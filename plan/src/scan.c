/** @file scan.c
 *  @brief Instrument scan characteristics and timing
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>

#include "bisect.h"
#include "solar.h"

#define CCD_NOMINAL_MODE_INTEGRATION_TIME 0.118  /* sec */

typedef struct
{
   double lon;
   double lat;
}
Surface_Point_Type;

typedef struct
{
   double mirror_x;  /* microradian, +toward East */
   double mirror_y;  /* microradian, +toward South */
}
Scan_Angle_Type;

typedef struct
{
   double integration_time;
   double position_dwell;
   double scan_reset;
   double scan_timing_margin;
   int scan_step_quantum;
   /* The number of steps in any scan is an integer multiple of SCAN_STEP_QUANTUM.
    * This quantization reduces the total number of unique command blocks that need
    * to be validated to support TEMPO operations. This approach limits the impact
    * to operations in case the hardware simulator stops working and no new command
    * blocks can be validated. By default scan_step_quantum=1.
    */

}
Step_Config_Type;

typedef struct
{
   Surface_Point_Type pt;
   Scan_Angle_Type ang;
   Step_Config_Type dt;
   double width;
   int num;
}
Surface_Region_Type;

#define SPLIT_SCAN_TYPE_PRIVATE_DATA \
   Surface_Point_Type scan_beg; \
   Surface_Point_Type scan_end; \
   Surface_Point_Type ctrl; \
   Scan_Angle_Type scan_beg_angle; \
   Scan_Angle_Type scan_end_angle; \
   Step_Config_Type dt; \
   double weight; \
   int base_scan_method_index; \
   int num_repeats_cbm;

#define TWILIGHT_SCAN_TYPE_PRIVATE_DATA \
   Surface_Region_Type east; \
   Surface_Region_Type west;

#define SCAN_TYPE_PRIVATE_DATA \
   double min_sun_angle; \
   double max_sza; \
   double short_scan_frac; \
   Surface_Point_Type scan_beg; \
   Surface_Point_Type scan_end; \
   Scan_Angle_Type scan_beg_angle; \
   Scan_Angle_Type scan_end_angle; \
   Surface_Point_Type day_beg; \
   Surface_Point_Type day_end; \
   Step_Config_Type dt; \
   double step_size; \
   int num_scan_steps; \
   uint16_t scan_type;
#include "scan.h"

int scan_quantize_num_steps_ceil (double fnum, int q)
{
   return q * ceil (fnum/q);
}

int scan_quantize_num_steps_floor (double fnum, int q)
{
   return q * floor (fnum/q);
}

static void free_scan_type (Scan_Type *st)
{
   if (NULL == st)
     return;
   FREE(st);
}

static Scan_Type *new_scan_type (uint16_t scan_type)
{
   Scan_Type *st = NULL;

   if (NULL == (st = (Scan_Type *)MALLOC (sizeof *st)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)st, 0, sizeof *st);

   st->scan_type = scan_type;

   return st;
}

static int read_limits (config_t *cfg, Scan_Type *st)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "limits_config")))
     goto return_error;

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "max_sza", &st->max_sza))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "min_sun_angle", &st->min_sun_angle)))
     goto return_error;

   /* hard-coded safety constraint */
   if ((st->min_sun_angle < 60.0) || (st->min_sun_angle > 180.0))
     {
        fprintf (stderr, "*** SAFETY CONSTRAINT: limits_config.min_sun_angle must be in the range [60.0, 180.0] degrees\n");
        goto return_error;
     }

   return 0;

return_error:
     tell_verror (TELL_INVALID_PARM_ERROR,
                  "%s: accessing limits_config in param file: %s",
                  __func__, config_error_file (cfg));
   return -1;
}

static int read_surface_point (config_setting_t *s, const char *name,
                               Surface_Point_Type *pt)
{
   config_setting_t *sub;

   if (NULL == (sub = config_setting_get_member (s, name)))
     return 0;

   if ((CONFIG_TRUE != config_setting_lookup_float (sub, "lon", &pt->lon))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "lat", &pt->lat)))
     return -1;

   return 1;
}

static int read_scan_angles (config_setting_t *s, const char *name,
                             Scan_Angle_Type *pt)
{
   config_setting_t *sub;
   double nan_value = nan("");

   pt->mirror_x = nan_value;
   pt->mirror_y = nan_value;

   if (NULL == (sub = config_setting_get_member (s, name)))
     return 0;

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "mirror_x", &pt->mirror_x))
     {
        pt->mirror_x = nan_value;
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "mirror_y", &pt->mirror_y))
     pt->mirror_y = 0.0;

   /* Using the law of reflection, tilting the mirror by an angle X
    * moves the (reflected) line of sight by an angle 2X.
    *
    * Obviously, the scan step size and start/end positions must be
    * in the same coordinates and have the same units.
    *
    * It's convenient to have the step size in field of regard
    * coordinates because that's directly comparable to the angular
    * size of features on the ground.
    *
    * Here, therefore, we use the law of reflection to convert
    * instrument mirror tilt angle coordinates to angular coordinates
    * in the field of regard.
    */

   pt->mirror_x *= 2;
   pt->mirror_y *= 2;

   return 1;
}

static int bsearch_f (float t, const float *x, int n)
{
   int n0, n1, n2;
   float xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

static int lookup_int_time_and_dwell_time (const char *ccdtiming_path, int num_coadds,
                                           Step_Config_Type *dt)
{
   float *int_time = NULL;
   size_t s_num, s_num_coadds;
   const char *grp_name;
   int ncid, dimid, start[2], count[2];
   int grp_long, grp_short, grp_nom, grp;
   int k, is_nominal_mode, num_int_lines;
   float integration_time, total_time;
   float max_short_int_time, min_long_int_time, nominal_int_time;
   int status = -1;

   integration_time = (float) dt->integration_time;

   if (0 != TIO_open (ccdtiming_path, NC_NOWRITE, &ncid))
     return -1;

   if (0 != TIO_inq_dim (ncid, "coadds", &dimid, &s_num_coadds))
     goto return_status;
   if (num_coadds > (int) s_num_coadds)
     {
        TIO_close (ncid);
        return 1;
     }

   if (0 != TIO_inq_dim (ncid, "int_lines", &dimid, &s_num))
     goto return_status;
   num_int_lines = s_num;

   count[0] = 1;

   if (0 != TIO_inq_grp (ncid, "nominal_mode", &grp_nom))
     goto return_status;
   start[0] = 0;
   if (0 != TIO_get_var_section (grp_nom, "integration_time", start, count, NC_FLOAT, &nominal_int_time))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, "short_mode", &grp_short))
     goto return_status;
   start[0] = num_int_lines-1;
   if (0 != TIO_get_var_section (grp_short, "integration_time", start, count, NC_FLOAT, &max_short_int_time))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, "long_mode", &grp_long))
     goto return_status;
   start[0] = 0;
   if (0 != TIO_get_var_section (grp_long, "integration_time", start, count, NC_FLOAT, &min_long_int_time))
     goto return_status;

   is_nominal_mode = ((max_short_int_time < integration_time)
                      && (integration_time < min_long_int_time));

   if (is_nominal_mode)
     {
        grp = grp_nom;
        grp_name = "nominal_mode";
     }
   else if (integration_time < nominal_int_time)
     {
        grp = grp_short;
        grp_name = "short_mode";
     }
   else
     {
        grp = grp_long;
        grp_name = "long_mode";
     }

   if (is_nominal_mode)
     {
        count[0] = 1;
        start[0] = num_coadds-1;
        if (0 != TIO_get_var_section (grp, "total_time", start, count, NC_FLOAT, &total_time))
          goto return_status;

        if (Plan_Verbose)
          {
             fprintf (stderr, "%d coadds x %f sec => %s => pos total_time = %f sec\n",
                      num_coadds, integration_time, grp_name, total_time);
          }

        dt->integration_time = (double) nominal_int_time;
        dt->position_dwell = (double) total_time;
     }
   else
     {
        if (NULL == (int_time = (float *)MALLOC (num_int_lines * sizeof(float))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto return_status;
          }
        start[0] = 0;
        count[0] = num_int_lines;
        if (0 != TIO_get_var_section (grp, "integration_time", start, count, NC_FLOAT, int_time))
          goto return_status;

        if (integration_time < int_time[0]) {
           k = 0;
        } else if (integration_time > int_time[num_int_lines-1]) {
           status = 1;
           goto quiet_return_status;
        } else {
           k = bsearch_f (integration_time, int_time, num_int_lines);
        }

        start[0] = k;               /* k = int_lines-1 */
        start[1] = num_coadds-1;
        count[0] = 1;
        count[1] = 1;
        if (0 != TIO_get_var_section (grp, "total_time", start, count, NC_FLOAT, &total_time))
          goto return_status;

        if (Plan_Verbose)
          {
             fprintf (stderr, "%d coadds x %f sec (%d int lines) => %s => pos total_time = %f sec\n",
                      num_coadds, int_time[k], k+1, grp_name, total_time);
          }

        dt->integration_time = (double) int_time[k];
        dt->position_dwell = (double) total_time;
     }

   status = 0;

return_status:
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: lookup failed", __func__);
     }
quiet_return_status:
   FREE(int_time);
   TIO_close (ncid);
   return status;
}

static int read_step_config (config_t *cfg, config_setting_t *s, Step_Config_Type *dt)
{
   config_setting_t *sub;
   const char *timing_file_var = "refdata_config.ccdtiming_path";
   const char *ccdtiming_path;
   char *path = NULL;
   double frame_transfer_time;
   double readout_time;
   int num_coadds, scan_step_quantum, status;

   if (NULL == (sub = config_setting_get_member (s, "step_config")))
     return -1;

   if ((CONFIG_TRUE != config_setting_lookup_float (sub, "integration_time", &dt->integration_time))
       || (CONFIG_TRUE != config_setting_lookup_int (sub, "num_coadds", &num_coadds))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "scan_reset", &dt->scan_reset))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "scan_timing_margin", &dt->scan_timing_margin)))
     return -1;

   /* optional user-specified value */
   if (CONFIG_TRUE == config_setting_lookup_int (sub, "scan_step_quantum", &scan_step_quantum))
     dt->scan_step_quantum = scan_step_quantum;
   else
     dt->scan_step_quantum = 1;

   /* Try to get the integration time and dwell time from the lookup table */
   if (CONFIG_TRUE == config_lookup_string (cfg, timing_file_var, &ccdtiming_path))
     {
        if (NULL == (path = expand_string (ccdtiming_path)))
          return -1;
        status = lookup_int_time_and_dwell_time (path, num_coadds, dt);
        FREE(path);
        if (status < 0 )
          {
             return status;
          }
        else if (status)
          {
             fprintf (stderr, "*** WARNING: %s: unable to use CCD timing table, scan timing is approximate\n", __func__);
          }
     }
   else
     {
        fprintf (stderr, "*** WARNING: %s: config file variable '%s' is not defined (CCD timing file)\n",
                 __func__, timing_file_var);
     }

   /* If we cannot use the lookup table, fall back to a crude approximation */
   if ((CONFIG_TRUE != config_setting_lookup_float (sub, "frame_transfer_time", &frame_transfer_time))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "readout_time", &readout_time)))
     return -1;

   dt->position_dwell = (num_coadds * (dt->integration_time + frame_transfer_time)
                         + frame_transfer_time + readout_time);

   /* For shorter than nominal integration times, the instrument
    * must generate a delay frame to avoid overloading the data formatter.
    * The delay frame increases the time spent at each mirror position.
    */
   if (dt->integration_time <= CCD_NOMINAL_MODE_INTEGRATION_TIME)
     {
        dt->position_dwell += (frame_transfer_time + readout_time);
     }

   return 0;
}

static int read_master_table_step_size (config_t *cfg, double *step_size)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "master_table_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing master_table_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "step_size", step_size))
     return -1;

   return 0;
}

static int read_scan_config (config_t *cfg, Scan_Type *st)
{
   config_setting_t *s;
   int have_scan_beg_angle, have_scan_end_angle;

   if (NULL == (s = config_lookup (cfg, "scan_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing scan_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_scan_beg_angle = read_scan_angles (s, "scan_beg_angle", &st->scan_beg_angle)) < 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading scan_config.scan_beg_angle': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_scan_end_angle = read_scan_angles (s, "scan_end_angle", &st->scan_end_angle)) < 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading scan_config.scan_end_angle': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_scan_beg_angle == 0)
       && (1 != read_surface_point (s, "scan_beg", &st->scan_beg)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'scan_beg': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   st->num_scan_steps = -1;
   if (CONFIG_TRUE == config_setting_lookup_int (s, "scan_num_steps", &st->num_scan_steps))
     {
        double nan_value = nan("");
        st->scan_end.lon = nan_value;
        st->scan_end.lat = nan_value;
     }

   if ((have_scan_end_angle == 0)
       && (1 != read_surface_point (s, "scan_end", &st->scan_end)))
     {
        if (st->num_scan_steps <= 0)
          {
             tell_verror (TELL_INVALID_PARM_ERROR, "%s: reading surface point 'scan_end': %s",
                          __func__, config_error_file (cfg));
             return -1;
          }
     }

   if (1 != read_surface_point (s, "day_beg", &st->day_beg))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'day_beg': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (1 != read_surface_point (s, "day_end", &st->day_end))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'day_end': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   /* Controls overlap of short scans; (0 <= f <= 1)
    * values >0.5 give non-zero overlap.
    * value = 0.5 gives gapless coverage with no overlap
    * values <0.5 probably aren't particularly useful
    */
   if (CONFIG_TRUE != config_setting_lookup_float (s, "short_scan_frac", &st->short_scan_frac))
     {
        st->short_scan_frac = 0.5;
     }

   if (0 != read_step_config (cfg, s, &st->dt))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading step_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_master_table_step_size (cfg, &st->step_size))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading step_size: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   /* Because scan plans are closely tied to the coordinates of points
    * on the Earth's surface, it is convenient for the planning calculations
    * to use the scan step as size in field of regard coordinates
    * (e.g. as an angular step across the field of regard).
    * By the law of reflection, the step size in the field of regard
    * is twice the mirror tilt angle. The value in the parameter file
    * is the mirror tilt angle, so we need to multiply that value by two:
    */
   st->step_size *= 2;

   return 0;
}

static int read_surface_region (config_t *cfg, config_setting_t *s, const char *name, Surface_Region_Type *reg)
{
   config_setting_t *sub;
   int have_angles;

   if (NULL == (sub = config_setting_get_member (s, name)))
     return -1;

   if ((have_angles = read_scan_angles (sub, "angles", &reg->ang)) < 0)
     return -1;

   if ((have_angles == 0)
       && (1 != read_surface_point (sub, "point", &reg->pt)))
     return -1;

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "width", &reg->width))
     return -1;
   if (CONFIG_TRUE != config_setting_lookup_int (sub, "num_subdivisions", &reg->num))
     return -1;

   if (reg->num <= 0) reg->num = 1;

   /* Step config may be region specific */
   if (0 != read_step_config (cfg, sub, &reg->dt))
     {
        if (0 != read_step_config (cfg, s, &reg->dt))
          return -1;
     }

   return 0;
}

static int read_twilight_scan_config (config_t *cfg, Twilight_Scan_Type *twilight_scan)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "twilight_scan_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing twilight_scan_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_region (cfg, s, "east_region", &twilight_scan->east))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'twilight_scan_config:east': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_region (cfg, s, "west_region", &twilight_scan->west))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'twilight_scan_config:west': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int read_split_scan_config (config_t *cfg, Split_Scan_Type *sst,
                                   const char *name)
{
   config_setting_t *s;
   int have_beg_angle, have_end_angle;

   if (NULL == (s = config_lookup (cfg, name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing split_scan_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_beg_angle = read_scan_angles (s, "scan_beg_angle", &sst->scan_beg_angle)) < 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading scan_config.scan_beg_angle': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_end_angle = read_scan_angles (s, "scan_end_angle", &sst->scan_end_angle)) < 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading scan_config.scan_end_angle': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_beg_angle == 0)
       && (1 != read_surface_point (s, "scan_beg", &sst->scan_beg)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'split_scan_config:east': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((have_end_angle == 0)
       && (1 != read_surface_point (s, "scan_end", &sst->scan_end)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'split_scan_config:west': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (1 != read_surface_point (s, "control", &sst->ctrl))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'split_scan_config:control': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "weight", &sst->weight))
     sst->weight = 1.0;

   if (0 != read_step_config (cfg, s, &sst->dt))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading step_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int read_params (config_t *cfg, Scan_Type *st)
{
   if (0 != read_limits (cfg, st))
     return -1;

   if (0 != read_scan_config (cfg, st))
     return -1;

   return 0;
}

typedef struct
{
   Solar_Geom_Type *sgt;
   Surface_Point_Type pt;
   double max_sza;
}
SZA_Bisect_Type;

static int dsza_vs_time (double jd_utc, double *dsza, void *v)
{
   SZA_Bisect_Type *b = (SZA_Bisect_Type *)v;
   Surface_Point_Type *pt = &b->pt;
   Solar_Geom_Type *sgt = b->sgt;
   double sza;

   if (0 != sgt->sgt_solar_zenith_angle (sgt, jd_utc, pt->lon, pt->lat, &sza))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: evaluating solar zenith angle",
                     __func__);
        return -1;
     }

   *dsza = sza - b->max_sza;

   return 0;
}

typedef struct
{
   Solar_Geom_Type *sgt;
   double target_sun_angle;
}
Sun_Angle_Bisect_Type;

static int sun_angle_vs_time (double jd_utc, double *dsa, void *v)
{
   Sun_Angle_Bisect_Type *b= (Sun_Angle_Bisect_Type *)v;
   Solar_Geom_Type *sgt = b->sgt;
   double sun_angle;

   if (0 != sgt->sgt_sat_sun_position (sgt, jd_utc, &sun_angle, NULL, NULL))
     return -1;
   *dsa = sun_angle - b->target_sun_angle;

   return 0;
}

static int find_safe_limit_time (Solar_Geom_Type *sgt,
                                 double min_sun_angle, int is_morning,
                                 double jd_utc, double *jd_utc_safe)
{
   Sun_Angle_Bisect_Type b;
   double jd_utc1, jd_utc2;

   b.sgt = sgt;
   b.target_sun_angle = min_sun_angle;

   if (is_morning)
     {
        jd_utc1 = jd_utc - 5.0/24;
        jd_utc2 = jd_utc + 1.0/24;
     }
   else
     {
        jd_utc1 = jd_utc - 1.0/24;;
        jd_utc2 = jd_utc + 5.0/24;
     }

   if (0 != bisection (sun_angle_vs_time, jd_utc1, jd_utc2, &b, &jd_utc))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: bisection failed (is_morning=%d, sun_angle=%0.2f not found between jd1=%0.2f jd2=%02f)",
                     __func__, is_morning, min_sun_angle, jd_utc1, jd_utc2);
        return -1;
     }

   *jd_utc_safe = jd_utc;

   return 0;
}

int scan_irradiance_time (Solar_Geom_Type *sgt, double irr_sun_angle, int after_midnight,
                          double jd_utc, double *jd_utc_irr)
{
   Sun_Angle_Bisect_Type b;
   double jd_utc1, jd_utc2;

   /* Assume: jd_utc = local midnight at the satellite's longitude */

   b.sgt = sgt;
   b.target_sun_angle = irr_sun_angle;

   if (0 == after_midnight)
     {
        /* Expect the nominal sun angle to occur within 3 hours BEFORE midnight */
        jd_utc1 = jd_utc - 3.0/24;
        jd_utc2 = jd_utc;
     }
   else
     {
        /* Expect the nominal sun angle to occur within 3 hours AFTER midnight */
        jd_utc1 = jd_utc;
        jd_utc2 = jd_utc + 3.0/24;
     }

   if (0 != bisection (sun_angle_vs_time, jd_utc1, jd_utc2, &b, &jd_utc))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: bisection failed (sun_angle=%0.2f not found between jd1=%0.2f jd2=%02f)",
                     __func__, irr_sun_angle, jd_utc1, jd_utc2);
        return -1;
     }

   *jd_utc_irr = jd_utc;

   return 0;
}

static double floor_sec (double t_days)
{
   return floor (t_days * SEC_PER_DAY) / SEC_PER_DAY;
}

static double ceil_sec (double t_days)
{
   return ceil (t_days * SEC_PER_DAY) / SEC_PER_DAY;
}

int scan_sza_time (Solar_Geom_Type *sgt, double max_sza, double jd_utc,
                   double lon, double lat, int is_start, double *jd_utc_sza)
{
   SZA_Bisect_Type b = {0};
   double jd_utc1, jd_utc2;

   /* Assume jd_utc = local midnight (beginning of 24 hour period to examine)
    *                 in the satellite time-zone */
   if (is_start)
     {
        jd_utc1 = jd_utc;
        jd_utc2 = jd_utc + 0.5;
     }
   else
     {
        jd_utc1 = jd_utc + 0.5;
        jd_utc2 = jd_utc + 1.0;
     }

   b.max_sza = max_sza;
   b.sgt = sgt;
   b.pt.lon = lon;
   b.pt.lat = lat;

   if (0 != bisection (dsza_vs_time, jd_utc1, jd_utc2, &b, jd_utc_sza))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: bisection failed", __func__);
        return -1;
     }

   return 0;
}

int scan_limit_times (const Scan_Type *st, double jd_utc,
                      Solar_Geom_Type *sgt,
                      Scan_Limit_Times_Type *slt)
{
   SZA_Bisect_Type b = {0};
   double jd_utc1, jd_utc2, jd_utc_midpoint;

   /* Assume: jd_utc  = local midnight
    *         jd_utc1 = local noon,
    *         jd_utc2 = local midnight,
    * where "local" means, e.g. at satellite longitude
    */
   jd_utc1 = jd_utc + 0.5;
   jd_utc2 = jd_utc1 + 0.5;

   b.max_sza = st->max_sza;
   b.sgt = sgt;

   /* Earliest morning time when SZA=max(SZA) at eastern point */
   b.pt = st->day_beg; /* struct copy */
   if (0 != bisection (dsza_vs_time, jd_utc, jd_utc1, &b, &slt->jd_utc_beg))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: bisection failed (jd_utc_beg)",
                     __func__);
        return -1;
     }
   /* Latest evening time when SZA=max(SZA) at western point */
   b.pt = st->day_end; /* struct copy */
   if (0 != bisection (dsza_vs_time, jd_utc1, jd_utc2, &b, &slt->jd_utc_end))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: bisection failed (jd_utc_end)",
                     __func__);
        return -1;
     }

   /* Find safe limit times (when it's ok to have the aperture open) */
   if (0 != find_safe_limit_time (sgt, st->min_sun_angle, 1, slt->jd_utc_beg, &slt->jd_utc_beg_safe))
     return -1;
   if (0 != find_safe_limit_time (sgt, st->min_sun_angle, 0, slt->jd_utc_end, &slt->jd_utc_end_safe))
     return -1;
   if (slt->jd_utc_end_safe <= slt->jd_utc_beg_safe)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid sun-angle safe interval",
                     __func__);
        return -1;
     }

   /* Make sure that radiance scan interval is inside the safe interval */
   if (slt->jd_utc_beg < slt->jd_utc_beg_safe)
     {
        slt->jd_utc_beg = slt->jd_utc_beg_safe;
     }
   if (slt->jd_utc_end > slt->jd_utc_end_safe)
     {
        slt->jd_utc_end = slt->jd_utc_end_safe;
     }
   if (slt->jd_utc_end <= slt->jd_utc_beg)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed sun angle safety constraint",
                     __func__);
        return -1;
     }

   if (Plan_Verbose > 1)
     {
        char buf[32];
        if (0 != mkjdtimestr (jd_utc, buf, sizeof(buf)))
          return -1;
        fprintf (stderr, "%s max scan %0.3f hr (max safe = %0.3f hr)\n", buf,
                 (slt->jd_utc_end      - slt->jd_utc_beg     ) * SEC_PER_DAY / 3600.0,
                 (slt->jd_utc_end_safe - slt->jd_utc_beg_safe) * SEC_PER_DAY / 3600.0);
     }

   jd_utc_midpoint = 0.5 * (slt->jd_utc_beg + slt->jd_utc_end);

   /* Full illumination is defined to begin when SZA=max(SZA)
    * first occurs at the western point.
    */
   b.pt.lat = st->day_end.lat;
   b.pt.lon = st->day_end.lon;
   if (0 != bisection (dsza_vs_time, slt->jd_utc_beg, jd_utc_midpoint, &b, &slt->jd_utc_beg_full))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: bisection failed (jd_utc_beg_full)",
                     __func__);
        return -1;
     }

   /* Full illumination is defined to end when SZA=max(SZA)
    * last occurs at the eastern point.
    */
   b.pt.lat = st->day_beg.lat;
   b.pt.lon = st->day_beg.lon;
   if (0 != bisection (dsza_vs_time, jd_utc_midpoint, slt->jd_utc_end, &b, &slt->jd_utc_end_full))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: bisection failed (jd_utc_end_full)",
                     __func__);
        return -1;
     }

   /* These limit times are somewhat arbitrary, and it's nice to start
    * scans on whole second boundaries, so adjust accordingly
    */
   slt->jd_utc_beg = ceil_sec (slt->jd_utc_beg);
   slt->jd_utc_end = floor_sec (slt->jd_utc_end);
   slt->jd_utc_beg_full = ceil_sec (slt->jd_utc_beg_full);
   slt->jd_utc_end_full = floor_sec (slt->jd_utc_end_full);
   slt->jd_utc_beg_safe = ceil_sec (slt->jd_utc_beg_safe);
   slt->jd_utc_end_safe = floor_sec (slt->jd_utc_end_safe);
   slt->user_imposed_start_time = 0;

   return 0;
}

static int scan_beg_angle (const Scan_Type *st, double *mirror_x, double *mirror_y)
{
   const Scan_Angle_Type *sa = &st->scan_beg_angle;
   if (isnan(sa->mirror_x) || isnan(sa->mirror_y))
     return -1;
   *mirror_x = sa->mirror_x;
   *mirror_y = sa->mirror_y;
   return 0;
}

static int scan_end_angle (const Scan_Type *st, double *mirror_x, double *mirror_y)
{
   const Scan_Angle_Type *sa = &st->scan_end_angle;
   if (isnan(sa->mirror_x) || isnan(sa->mirror_y))
     return -1;
   *mirror_x = sa->mirror_x;
   *mirror_y = sa->mirror_y;
   return 0;
}

static int scan_beg_point (const Scan_Type *st, double *lon, double *lat)
{
   *lon = st->scan_beg.lon;
   *lat = st->scan_beg.lat;
   return 0;
}

static int scan_end_point (const Scan_Type *st, double *lon, double *lat)
{
   *lon = st->scan_end.lon;
   *lat = st->scan_end.lat;
   return 0;
}

static int scan_day_beg_point (const Scan_Type *st, double *lon, double *lat)
{
   *lon = st->day_beg.lon;
   *lat = st->day_beg.lat;
   return 0;
}

static int scan_day_end_point (const Scan_Type *st, double *lon, double *lat)
{
   *lon = st->day_end.lon;
   *lat = st->day_end.lat;
   return 0;
}

static int scan_num_steps (const Scan_Type *st)
{
   return st->num_scan_steps;
}

static int twilight_scan_region_angles (const Twilight_Scan_Type *tst, int is_east, double *mirror_x, double *mirror_y,
                                        double *width, int *num)
{
   const Surface_Region_Type *reg = is_east ? &tst->east : &tst->west;
   const Scan_Angle_Type *ang = &reg->ang;
   if (isnan(ang->mirror_x) || isnan(ang->mirror_y))
     return -1;
   *mirror_x = reg->ang.mirror_x;
   *mirror_y = reg->ang.mirror_y;
   *width = reg->width;
   *num = reg->num;
   return 0;
}

static int twilight_scan_region (const Twilight_Scan_Type *tst, int is_east, double *lon, double *lat,
                                 double *width, int *num)
{
   const Surface_Region_Type *reg = is_east ? &tst->east : &tst->west;

   *lon = reg->pt.lon;
   *lat = reg->pt.lat;
   *width = reg->width;
   *num = reg->num;

   return 0;
}

static int split_scan_region_angles (const Split_Scan_Type *sst,
                                     double *beg_x, double *beg_y,
                                     double *end_x, double *end_y)
{
   const Scan_Angle_Type *sb = &sst->scan_beg_angle;
   const Scan_Angle_Type *se = &sst->scan_end_angle;
   if (isnan(sb->mirror_x) || isnan(sb->mirror_y) || isnan(se->mirror_x) || isnan(se->mirror_y))
     return -1;
   *beg_x = sb->mirror_x;
   *beg_y = sb->mirror_y;
   *end_x = se->mirror_x;
   *end_y = se->mirror_y;
   return 0;
}

static int split_scan_region (const Split_Scan_Type *sst,
                              double *beg_lon, double *beg_lat,
                              double *end_lon, double *end_lat)
{
   *beg_lon = sst->scan_beg.lon;
   *beg_lat = sst->scan_beg.lat;
   *end_lon = sst->scan_end.lon;
   *end_lat = sst->scan_end.lat;
   return 0;
}

static int split_scan_control (const Split_Scan_Type *sst,
                               double *lon, double *lat)
{
   *lon = sst->ctrl.lon;
   *lat = sst->ctrl.lat;
   return 0;
}

static double scan_min_sun_angle (const Scan_Type *st)
{
   return st->min_sun_angle;
}

static double scan_max_sza (const Scan_Type *st)
{
   return st->max_sza;
}

static double scan_integration_time (const Scan_Type *st)
{
   return st->dt.integration_time;
}

static double twilight_scan_integration_time (const Twilight_Scan_Type *tst, int is_east)
{
   const Surface_Region_Type *reg = is_east ? &tst->east : &tst->west;
   return reg->dt.integration_time;
}

static double split_scan_integration_time (const Split_Scan_Type *sst)
{
   return sst->dt.integration_time;
}

static double scan_step_size (const Scan_Type *st)
{
   return st->step_size;
}

static double short_scan_frac (const Scan_Type *st)
{
   return st->short_scan_frac;
}

/* scan duration [days] */
static double __scan_duration_days (const Step_Config_Type *dt, int num_steps)
{
   double duration = (dt->position_dwell * (num_steps + 1)
                      + 2*dt->scan_reset
                      + dt->scan_timing_margin) / SEC_PER_DAY;
   /* adjust scan duration to an integer number of seconds [in units of days] */
   return ceil_sec (duration);
}
static int __scan_num_steps_in_duration (const Step_Config_Type *dt, double duration_days)
{
   /* underestimate the number of steps to allow margin */
   double fnp1 = (floor_sec(duration_days) * SEC_PER_DAY
              - dt->scan_timing_margin - 2*dt->scan_reset) / dt->position_dwell;
   int n = scan_quantize_num_steps_floor (fnp1 - 1, dt->scan_step_quantum);
   return (n > 0) ? n : 0;
}

static double scan_duration (const Scan_Type *st, int num_steps)
{
   return __scan_duration_days (&st->dt, num_steps);
}
static int scan_num_steps_in_duration (const Scan_Type *st, double duration_days)
{
   return __scan_num_steps_in_duration (&st->dt, duration_days);
}

static double twilight_scan_duration (const Twilight_Scan_Type *tst, int is_east, int num_steps)
{
   const Surface_Region_Type *reg = is_east ? &tst->east : &tst->west;
   /* FIXME - should twilight scan duration be computed differently? */
   return __scan_duration_days (&reg->dt, num_steps);
}
static int twilight_scan_num_steps_in_duration (const Twilight_Scan_Type *tst, int is_east, double duration_days)
{
   const Surface_Region_Type *reg = is_east ? &tst->east : &tst->west;
   return __scan_num_steps_in_duration (&reg->dt, duration_days);
}

static int scan_print_params (const Scan_Type *st, const char *pprefix,
                              FILE *fp)
{
   const char *prefix = pprefix ? pprefix : "";
   (void) fprintf (fp, "%s %0.2f deg = max(SZA)\n",
                   prefix, st->max_sza);
   (void) fprintf (fp, "%s %0.2f deg = min(sun_angle)\n",
                   prefix, st->min_sun_angle);
   (void) fprintf (fp, "%s (%8.3f,%7.3f) = Eastern scan limit point\n",
                   prefix, st->scan_beg.lon, st->scan_beg.lat);
   if (0 == isnan (st->scan_end.lon))
     {
        (void) fprintf (fp, "%s (%8.3f,%7.3f) = Western scan limit point\n",
                        prefix, st->scan_end.lon, st->scan_end.lat);
     }
   else
     {
        (void) fprintf (fp, "%s %d = Number of scan steps\n",
                        prefix, st->num_scan_steps);
     }
   (void) fprintf (fp, "%s (%8.3f,%7.3f) = Eastern SZA control point\n",
                   prefix, st->day_beg.lon, st->day_beg.lat);
   (void) fprintf (fp, "%s (%8.3f,%7.3f) = Western SZA control point\n",
                   prefix, st->day_end.lon, st->day_end.lat);
   return 0;
}

static uint16_t query_scan_type (const Scan_Type *st)
{
   return st->scan_type;
}

static int scan_config_step_quantum (const Scan_Type *st)
{
   return st->dt.scan_step_quantum;
}

Scan_Type *scan_open (config_t *cfg, uint16_t scan_type)
{
   Scan_Type *st = NULL;

   if (NULL == (st = new_scan_type (scan_type)))
     return NULL;

   st->st_delete = free_scan_type;
   st->st_scan_duration = scan_duration;
   st->st_scan_num_steps_in_duration = scan_num_steps_in_duration;
   st->st_step_size = scan_step_size;
   st->st_short_scan_frac = short_scan_frac;
   st->st_integration_time = scan_integration_time;
   st->st_min_sun_angle = scan_min_sun_angle;
   st->st_max_sza = scan_max_sza;
   st->st_scan_beg_angle = scan_beg_angle;
   st->st_scan_end_angle = scan_end_angle;
   st->st_scan_beg = scan_beg_point;
   st->st_scan_end = scan_end_point;
   st->st_scan_day_beg = scan_day_beg_point;
   st->st_scan_day_end = scan_day_end_point;
   st->st_scan_num_steps = scan_num_steps;
   st->st_print_params = scan_print_params;
   st->st_scan_type = query_scan_type;
   st->st_scan_step_quantum = scan_config_step_quantum;

   if (0 != read_params (cfg, st))
     {
        free_scan_type (st);
        return NULL;
     }

   return st;
}

static void free_twilight_scan_type (Twilight_Scan_Type *tst)
{
   if (tst == NULL)
     return;
   FREE(tst);
}

Twilight_Scan_Type *twilight_scan_open (config_t *cfg)
{
   Twilight_Scan_Type *tst = NULL;

   if (NULL == (tst = (Twilight_Scan_Type *) MALLOC (sizeof *tst)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   tst->tst_delete = free_twilight_scan_type;
   tst->tst_twilight_scan_region = twilight_scan_region;
   tst->tst_twilight_scan_region_angles = twilight_scan_region_angles;
   tst->tst_twilight_scan_duration = twilight_scan_duration;
   tst->tst_twilight_scan_num_steps_in_duration = twilight_scan_num_steps_in_duration;
   tst->tst_twilight_integration_time = twilight_scan_integration_time;

   if (0 != read_twilight_scan_config (cfg, tst))
     {
        free_twilight_scan_type (tst);
        return NULL;
     }

   return tst;
}

static void free_split_scan_type (Split_Scan_Type *sst)
{
   if (sst == NULL)
     return;
   FREE(sst);
}

static int parse_scan_method_string (const char *pscan_method,
                                     char *base_scan_method, size_t size_bsm,
                                     char *config_group_name, size_t size_cgn,
                                     int *num_repeats_cbm)
{
   const char *delim = "-";
   char *scan_method = NULL;
   char *tok = NULL;
   int status = -1;

   *num_repeats_cbm = 0;

   if (NULL == (scan_method = strdup (pscan_method)))
     return -1;

   if (NULL == (tok = strtok (scan_method, delim)))
     goto free_and_return;
   if (0 != strcmp (tok, "split"))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unexpected scan method string: %s", __func__, scan_method);
        goto free_and_return;
     }

   if (NULL == (tok = strtok (NULL, delim)))
     goto free_and_return;
   strncpy (base_scan_method, tok, size_bsm-1);

   if (NULL == (tok = strtok (NULL, delim)))
     goto free_and_return;
   strncpy (config_group_name, tok, size_cgn-1);

   /* optional 3rd field to specify a custom scanning CBM */
   if (NULL != (tok = strtok (NULL, delim)))
     {
        if (1 != sscanf (tok, "%d", num_repeats_cbm))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: reading num_repeats_cbm field from scan method string: %s",
                          __func__, pscan_method);
             goto free_and_return;
          }
        if (*num_repeats_cbm < 0)
          {
             tell_verror (TELL_INVALID_PARM_ERROR, "%s: scan method string %s (num_repeats_cbm = %d)",
                          __func__, pscan_method, *num_repeats_cbm);
             goto free_and_return;
          }
     }

   status = 0;
free_and_return:
   FREE(scan_method);
   return status;
}

static int split_scan_base_scan_method (const Split_Scan_Type *sst)
{
   return sst->base_scan_method_index;
}

static int split_scan_num_repeats_cbm (const Split_Scan_Type *sst)
{
   return sst->num_repeats_cbm;
}

static double split_scan_weight (const Split_Scan_Type *sst)
{
   return sst->weight;
}

Split_Scan_Type *split_scan_open (config_t *cfg, const char *scan_method)
{
   Split_Scan_Type *sst = NULL;
#define MAX_SIZE_METHOD_NAME         (16)
#define MAX_SIZE_CONFIG_GROUP_NAME   (256)
   char base_scan_method[MAX_SIZE_METHOD_NAME];
   char config_group_name[MAX_SIZE_CONFIG_GROUP_NAME];
   int base_scan_method_index;
   int num_repeats_cbm;

   if (0 != parse_scan_method_string (scan_method,
                                      base_scan_method, MAX_SIZE_METHOD_NAME,
                                      config_group_name, MAX_SIZE_CONFIG_GROUP_NAME,
                                      &num_repeats_cbm))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: parsing scan_method string: %s",
                     __func__, scan_method);
        return NULL;
     }

   if (0 == strcmp (base_scan_method, "std"))
     base_scan_method_index = SCAN_SPLIT_STD;
   else if (0 == strcmp (base_scan_method, "opt1"))
     base_scan_method_index = SCAN_SPLIT_OPT1;
   else
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported base scan method: %s",
                     __func__, scan_method);
        return NULL;
     }

   if (NULL == (sst = (Split_Scan_Type *) MALLOC (sizeof *sst)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   sst->base_scan_method_index = base_scan_method_index;
   sst->num_repeats_cbm = num_repeats_cbm;
   sst->weight = 1.0;

   sst->sst_base_scan_method = split_scan_base_scan_method;
   sst->sst_num_repeats_cbm = split_scan_num_repeats_cbm;
   sst->sst_weight = split_scan_weight;
   sst->sst_delete = free_split_scan_type;
   sst->sst_scan_region_angles = split_scan_region_angles;
   sst->sst_scan_region = split_scan_region;
   sst->sst_scan_control = split_scan_control;
   sst->sst_scan_integration_time = split_scan_integration_time;

   if (0 != read_split_scan_config (cfg, sst, config_group_name))
     {
        free_split_scan_type (sst);
        return NULL;
     }

   return sst;
}
