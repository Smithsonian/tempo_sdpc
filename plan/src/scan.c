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

#include "bisect.h"
#include "solar.h"

typedef struct
{
   double lon;
   double lat;
}
Surface_Point_Type;

typedef struct
{
   double integration_time;
   double frame_transfer_time;
   double readout_time;
   double step_dwell;
   double scan_reset;
   double scan_timing_margin;
   int num_coadds;
}
Step_Config_Type;

typedef struct
{
   Surface_Point_Type pt;
   double width;
   int num;
}
Surface_Region_Type;

typedef struct
{
   Surface_Region_Type east;
   Surface_Region_Type west;
   Step_Config_Type dt;
}
Night_Scan_Type;

#define SCAN_TYPE_PRIVATE_DATA \
   double min_sun_angle; \
   double max_sza; \
   Surface_Point_Type scan_beg; \
   Surface_Point_Type scan_end; \
   Surface_Point_Type day_beg; \
   Surface_Point_Type day_end; \
   Step_Config_Type dt; \
   double step_size; \
   uint16_t scan_type; \
   Night_Scan_Type night_scan;
#include "scan.h"

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
     return -1;

   if ((CONFIG_TRUE != config_setting_lookup_float (sub, "lon", &pt->lon))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "lat", &pt->lat)))
     return -1;

   return 0;
}

static int read_step_config (config_setting_t *s, Step_Config_Type *dt)
{
   config_setting_t *sub;

   if (NULL == (sub = config_setting_get_member (s, "step_config")))
     return -1;

   if ((CONFIG_TRUE != config_setting_lookup_float (sub, "integration_time", &dt->integration_time))
       || (CONFIG_TRUE != config_setting_lookup_int (sub, "num_coadds", &dt->num_coadds))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "scan_reset", &dt->scan_reset))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "scan_timing_margin", &dt->scan_timing_margin))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "frame_transfer_time", &dt->frame_transfer_time))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "readout_time", &dt->readout_time)))
     return -1;

   dt->step_dwell =
     dt->num_coadds * (dt->integration_time + dt->frame_transfer_time)
       + dt->frame_transfer_time + dt->readout_time;

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

   if (NULL == (s = config_lookup (cfg, "scan_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing scan_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_point (s, "scan_beg", &st->scan_beg))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'scan_beg': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_point (s, "scan_end", &st->scan_end))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'scan_end': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_point (s, "day_beg", &st->day_beg))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'day_beg': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_point (s, "day_end", &st->day_end))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'day_end': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_step_config (s, &st->dt))
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

   return 0;
}

static int read_surface_region (config_setting_t *s, const char *name, Surface_Region_Type *reg)
{
   config_setting_t *sub;

   if (0 != read_surface_point (s, name, &reg->pt))
     return -1;

   if (NULL == (sub = config_setting_get_member (s, name)))
     return -1;

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "width", &reg->width))
     return -1;
   if (CONFIG_TRUE != config_setting_lookup_int (sub, "num_subdivisions", &reg->num))
     return -1;

   if (reg->num <= 0) reg->num = 1;

   return 0;
}

static int read_night_scan_config (config_t *cfg, Night_Scan_Type *night_scan)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "night_scan_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing night_scan_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_region (s, "east_region", &night_scan->east))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'night_scan_config:east': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_surface_region (s, "west_region", &night_scan->west))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading surface point 'night_scan_config:west': %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != read_step_config (s, &night_scan->dt))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading night_scan_config:step_config: %s",
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

   if (0 != read_night_scan_config (cfg, &st->night_scan))
     return -1;

   return 0;
}

typedef struct
{
   const Scan_Type *st;
   Solar_Geom_Type *sgt;
   Surface_Point_Type pt;
}
SZA_Bisect_Type;

static int dsza_vs_time (double jd_utc, double *dsza, void *v)
{
   SZA_Bisect_Type *b = (SZA_Bisect_Type *)v;
   Surface_Point_Type *pt = &b->pt;
   Solar_Geom_Type *sgt = b->sgt;
   const Scan_Type *st = b->st;
   double sza;

   if (0 != sgt->sgt_solar_zenith_angle (sgt, jd_utc, pt->lon, pt->lat, &sza))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: evaluating solar zenith angle",
                     __func__);
        return -1;
     }

   *dsza = sza - st->max_sza;

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

   if (0 != sgt->sgt_sat_sun_angles (sgt, jd_utc, &sun_angle, NULL))
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
        jd_utc1 = jd_utc - 4.0/24;
        jd_utc2 = jd_utc + 0.5/24;
     }
   else
     {
        jd_utc1 = jd_utc - 0.5/24;;
        jd_utc2 = jd_utc + 4.0/24;
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

int scan_irradiance_time (Solar_Geom_Type *sgt, double irr_sun_angle,
                          double jd_utc, double *jd_utc_irr)
{
   Sun_Angle_Bisect_Type b;
   double jd_utc1, jd_utc2;

   /* Assume: jd_utc = local midnight at the satellite's longitude */

   b.sgt = sgt;
   b.target_sun_angle = irr_sun_angle;

   /* Expect the nominal sun angle to occur within 3 hours before midnight */
   jd_utc1 = jd_utc - 3.0/24;
   jd_utc2 = jd_utc;

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

   b.st = st;
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

   return 0;
}

static int scan_beg_point (const Scan_Type *st, double *lon, double *lat)
{
   if (lon) *lon = st->scan_beg.lon;
   if (lat) *lat = st->scan_beg.lat;
   return 0;
}

static int scan_end_point (const Scan_Type *st, double *lon, double *lat)
{
   if (lon) *lon = st->scan_end.lon;
   if (lat) *lat = st->scan_end.lat;
   return 0;
}

static int night_scan_region (const Scan_Type *st, int is_east, double *lon, double *lat,
                              double *width, int *num)
{
   const Surface_Region_Type *reg = is_east ? &st->night_scan.east : &st->night_scan.west;

   *lon = reg->pt.lon;
   *lat = reg->pt.lat;
   *width = reg->width;
   *num = reg->num;

   return 0;
}

static double scan_min_sun_angle (const Scan_Type *st)
{
   return st->min_sun_angle;
}

static double scan_integration_time (const Scan_Type *st)
{
   return st->dt.integration_time;
}

static double night_scan_integration_time (const Scan_Type *st)
{
   return st->night_scan.dt.integration_time;
}

static double scan_step_size (const Scan_Type *st)
{
   return st->step_size;
}

/* scan duration [days] */
static double __scan_duration_days (const Step_Config_Type *dt, int num_steps)
{
   double duration = (dt->step_dwell * num_steps
                      + 2*dt->scan_reset
                      + dt->scan_timing_margin) / SEC_PER_DAY;
   /* adjust scan duration to an integer number of seconds [in units of days] */
   return ceil_sec (duration);
}

static double scan_duration (const Scan_Type *st, int num_steps)
{
   return __scan_duration_days (&st->dt, num_steps);
}

static double night_scan_duration (const Scan_Type *st, int num_steps)
{
   /* FIXME - should night scan duration be computed differently? */
   return __scan_duration_days (&st->night_scan.dt, num_steps);
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
   (void) fprintf (fp, "%s (%8.3f,%7.3f) = Western scan limit point\n",
                   prefix, st->scan_end.lon, st->scan_end.lat);
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

Scan_Type *scan_open (config_t *cfg, uint16_t scan_type)
{
   Scan_Type *st = NULL;

   if (NULL == (st = new_scan_type (scan_type)))
     return NULL;

   st->st_delete = free_scan_type;
   st->st_scan_duration = scan_duration;
   st->st_step_size = scan_step_size;
   st->st_integration_time = scan_integration_time;
   st->st_min_sun_angle = scan_min_sun_angle;
   st->st_scan_beg = scan_beg_point;
   st->st_scan_end = scan_end_point;
   st->st_print_params = scan_print_params;
   st->st_scan_type = query_scan_type;

   st->st_night_scan_region = night_scan_region;
   st->st_night_scan_duration = night_scan_duration;
   st->st_night_integration_time = night_scan_integration_time;

   if (0 != read_params (cfg, st))
     {
        free_scan_type (st);
        return NULL;
     }

   return st;
}
