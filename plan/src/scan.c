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

#define SCAN_TYPE_PRIVATE_DATA \
   double min_sun_angle; \
   double max_sza; \
   Surface_Point_Type scan_beg; \
   Surface_Point_Type scan_end; \
   Surface_Point_Type day_beg; \
   Surface_Point_Type day_end; \
   Step_Config_Type dt; \
   double step_size;
#include "scan.h"

static void free_scan_type (Scan_Type *st)
{
   if (NULL == st)
     return;
   FREE(st);
}

static Scan_Type *new_scan_type (void)
{
   Scan_Type *st = NULL;

   if (NULL == (st = (Scan_Type *)MALLOC (sizeof *st)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)st, 0, sizeof *st);

   return st;
}

static int read_limits (config_t *cfg, Scan_Type *st)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "limits_config")))
     return -1;

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "max_sza", &st->max_sza))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "min_sun_angle", &st->min_sun_angle)))
     return -1;

   return 0;
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

static int read_step_size (config_t *cfg, double *step_size)
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

static int read_params (config_t *cfg, Scan_Type *st)
{
   config_setting_t *s;

   if (0 != read_limits (cfg, st))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing scan_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

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

   if (0 != read_step_size (cfg, &st->step_size))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading step_size: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

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
   double min_sun_angle;
}
Sun_Angle_Bisect_Type;

static int sun_angle_vs_time (double jd_utc, double *dsa, void *v)
{
   Sun_Angle_Bisect_Type *b= (Sun_Angle_Bisect_Type *)v;
   Solar_Geom_Type *sgt = b->sgt;
   double sun_angle;

   if (0 != sgt->sgt_sat_sun_angle (sgt, jd_utc, &sun_angle))
     return -1;
   *dsa = sun_angle - b->min_sun_angle;

   return 0;
}

static int fixup_sun_angle (Solar_Geom_Type *sgt,
                            double min_sun_angle,
                            int is_morning, double *jd_utc_sza)
{
   Sun_Angle_Bisect_Type b;
   double jd_utc, jd_utc1, jd_utc2, sun_angle;

   jd_utc = *jd_utc_sza;

   /* If the sun-angle is ok, do nothing */
   if (0 != sgt->sgt_sat_sun_angle (sgt, jd_utc, &sun_angle))
     return -1;
   if (sun_angle > min_sun_angle)
     return 0;

   b.sgt = sgt;
   b.min_sun_angle = min_sun_angle;

   if (is_morning)
     {
        jd_utc1 = jd_utc;
        jd_utc2 = jd_utc + 2.0/24;
     }
   else
     {
        jd_utc1 = jd_utc - 2.0/24;
        jd_utc2 = jd_utc;
     }

   if (0 != bisection (sun_angle_vs_time, jd_utc1, jd_utc2, &b, &jd_utc))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: bisection failed (sun angle < %0.2f between jd1=%0.2f jd2=%02f)",
                     __func__, min_sun_angle, jd_utc1, jd_utc2);
        return -1;
     }

   *jd_utc_sza = jd_utc;

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

   /* Ensure sun angle constraints are met */
   if (0 != fixup_sun_angle (sgt, st->min_sun_angle, 1, &slt->jd_utc_beg))
     return -1;
   if (0 != fixup_sun_angle (sgt, st->min_sun_angle, 0, &slt->jd_utc_end))
     return -1;

   if (slt->jd_utc_beg >= slt->jd_utc_end)
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

static double scan_min_sun_angle (const Scan_Type *st)
{
   return st->min_sun_angle;
}

static double scan_integration_time (const Scan_Type *st)
{
   return st->dt.integration_time;
}

static double scan_step_size (const Scan_Type *st)
{
   return st->step_size;
}

/* scan duration [days] */
static double scan_duration (const Scan_Type *st, int num_steps)
{
   const Step_Config_Type *dt = &st->dt;
   double duration = (dt->step_dwell * num_steps
                      + 2*dt->scan_reset
                      + dt->scan_timing_margin) / SEC_PER_DAY;
   return duration;
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

Scan_Type *scan_open (config_t *cfg)
{
   Scan_Type *st = NULL;

   if (NULL == (st = new_scan_type ()))
     return NULL;

   st->st_delete = free_scan_type;
   st->st_scan_duration = scan_duration;
   st->st_step_size = scan_step_size;
   st->st_integration_time = scan_integration_time;
   st->st_min_sun_angle = scan_min_sun_angle;
   st->st_scan_beg = scan_beg_point;
   st->st_scan_end = scan_end_point;
   st->st_print_params = scan_print_params;

   if (0 != read_params (cfg, st))
     {
        free_scan_type (st);
        return NULL;
     }

   return st;
}
