/** @file scan_methods.c
 *  @brief Support different scan planning methods
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

#include <tempo_geo.h>

#include "scan.h"
#include "plan_list.h"
#include "scan_methods.h"
#include "vis.h"

#define DEGTORAD       (M_PI/180.0)
#define DEGTOMICRORAD  (1.e6*DEGTORAD)

typedef struct
{
   double azimuth;
   double elevation;
}
AziElev_Type;

static int compute_scan_angles (const EarthPoint *pt, double sat_lon, AziElev_Type *apt)
{
   TempoGeoErr error;

   /* need sat_lon in radians */
   sat_lon /= DEGTORAD;

   if ((error = computeScanAngles (pt, sat_lon, SCAN_AZ_FIRST,
                                   TEMPO_FIRST_CORRECTION,
                                   &apt->azimuth, &apt->elevation)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computeScanAngles failed (%s)",
                     __func__, tempoGeoErrorString (error));
        return -1;
     }

   /* want values in microradians */

   apt->azimuth *= DEGTOMICRORAD;
   apt->elevation *= DEGTOMICRORAD;

   return 0;
}

static int radiance_scan_endpoints (const Scan_Type *st,
                                    const Solar_Geom_Type *solar_geom,
                                    AziElev_Type *beg,
                                    AziElev_Type *end)
{
   EarthPoint beg_pt={0}, end_pt={0};
   double sat_lon;

   if (0 != solar_geom->sgt_geosat_longitude(solar_geom, &sat_lon))
     return -1;

   if (0 != st->st_scan_beg (st, &beg_pt.theLon, &beg_pt.theLat))
     return -1;
   if (0 != st->st_scan_end (st, &end_pt.theLon, &end_pt.theLat))
     return -1;

   if (0 != compute_scan_angles (&beg_pt, sat_lon, beg))
     return -1;
   if (0 != compute_scan_angles (&end_pt, sat_lon, end))
     return -1;

   return 0;
}

static int nightlights_scan_endpoints (const Scan_Type *st,
                                       const Solar_Geom_Type *solar_geom,
                                       const Scan_Limit_Times_Type *limit_times,
                                       int is_east,
                                       AziElev_Type *beg,
                                       AziElev_Type *end)
{
   EarthPoint pt={0};
   AziElev_Type p;
   double sat_lon, width;
   double sub_width, sub_offset;
   int i, num;

   if (0 != solar_geom->sgt_geosat_longitude(solar_geom, &sat_lon))
     return -1;

   if (0 != st->st_night_scan_region (st, is_east, &pt.theLon, &pt.theLat, &width, &num))
     return -1;

   if (0 != compute_scan_angles (&pt, sat_lon, &p))
     return -1;

   beg->elevation = p.elevation;
   end->elevation = p.elevation;

   /* Scan adjacent subdivisions on successive days */
   i = ((int) limit_times->jd_utc_beg) % num;
   sub_width = width / num;
   sub_offset = p.azimuth + i * sub_width;

   /* width>0 means the region is eastward of the point.
    * Scans always begin on the eastern side. */
   if (width > 0)
     {
        beg->azimuth = sub_offset + sub_width;
        end->azimuth = sub_offset;
     }
   else
     {
        beg->azimuth = sub_offset;
        end->azimuth = sub_offset + sub_width;
     }

   return 0;
}

static int std_scan_table (const Scan_Type *st,
                           const AziElev_Type *beg, const AziElev_Type *end,
                           double *xstart, double *ystart, int *num_steps)
{
   double step_size;
   int max_num_steps;

   /* In this context, we care only about the absolute value of the
    * mirror step size */
   step_size = fabs(st->st_step_size (st));
   if (step_size == 0.0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: invalid mirror step size = %g",
                     __func__, step_size);
        return -1;
     }

   /* Coordinate system should ensure beg->azimuth > end->azimuth */
   max_num_steps = (beg->azimuth - end->azimuth) / step_size;

   if (max_num_steps <= 0)
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: max_num_steps = %d (expected > 0)",
                     __func__, max_num_steps);
        return -1;
     }

   *xstart = beg->azimuth;
   *ystart = beg->elevation;
   *num_steps = max_num_steps;

   return 0;
}

static Plan_List_Type *
std_plan (const Scan_Type *st, Solar_Geom_Type *solar_geom,
          const Scan_Limit_Times_Type *limit_times)
{
   Plan_List_Type *entry = NULL;
   AziElev_Type beg={0}, end={0};
   double time_full_scan, xstart, ystart;
   int num_steps;

   if (0 != radiance_scan_endpoints (st, solar_geom, &beg, &end))
     return NULL;
   if (0 != std_scan_table (st, &beg, &end, &xstart, &ystart, &num_steps))
     return NULL;

   time_full_scan = st->st_scan_duration (st, num_steps);

   if (NULL == (entry = plan_list_entry_alloc ()))
     return NULL;

   entry->tstart = limit_times->jd_utc_beg;
   entry->xstart = xstart;
   entry->ystart = ystart;
   entry->num_steps = num_steps;
   entry->scan_duration = time_full_scan * SEC_PER_DAY;
   entry->integration_time = st->st_integration_time (st);
   entry->num_repeats = ceil ((limit_times->jd_utc_end - limit_times->jd_utc_beg)
                              / time_full_scan);
   return entry;
}

static Plan_List_Type *
nightlights_dawn_plan (const Scan_Type *st, Solar_Geom_Type *solar_geom,
                       const Scan_Limit_Times_Type *limit_times)
{
   Plan_List_Type *entry = NULL;
   AziElev_Type beg={0}, end={0};
   double time_full_scan, xstart, ystart;
   int num_steps, num_repeats;

   if (0 != nightlights_scan_endpoints (st, solar_geom, limit_times, 0, &beg, &end))
     return NULL;
   if (0 != std_scan_table (st, &beg, &end, &xstart, &ystart, &num_steps))
     return NULL;

   time_full_scan = st->st_night_scan_duration (st, num_steps);

   /* may be zero */
   num_repeats = floor ((limit_times->jd_utc_beg - limit_times->jd_utc_beg_safe)
                        / time_full_scan);

   if (NULL == (entry = plan_list_entry_alloc ()))
     return NULL;

   entry->tstart = limit_times->jd_utc_beg_safe;
   entry->xstart = xstart;
   entry->ystart = ystart;
   entry->num_steps = num_steps;
   entry->scan_duration = time_full_scan * SEC_PER_DAY;
   entry->integration_time = st->st_night_integration_time (st);
   entry->num_repeats = num_repeats;

   return entry;
}

static Plan_List_Type *
nightlights_dusk_plan (const Scan_Type *st, Solar_Geom_Type *solar_geom,
                       const Scan_Limit_Times_Type *limit_times)
{
   Plan_List_Type *entry = NULL;
   AziElev_Type beg={0}, end={0};
   double time_full_scan, xstart, ystart;
   int num_steps, num_repeats;

   if (0 != nightlights_scan_endpoints (st, solar_geom, limit_times, 1, &beg, &end))
     return NULL;
   if (0 != std_scan_table (st, &beg, &end, &xstart, &ystart, &num_steps))
     return NULL;

   time_full_scan = st->st_night_scan_duration (st, num_steps);

   /* may be zero */
   num_repeats = floor ((limit_times->jd_utc_end_safe - limit_times->jd_utc_end)
                        / time_full_scan);

   if (NULL == (entry = plan_list_entry_alloc ()))
     return NULL;

   entry->tstart = limit_times->jd_utc_end;
   entry->xstart = xstart;
   entry->ystart = ystart;
   entry->num_steps = num_steps;
   entry->scan_duration = time_full_scan * SEC_PER_DAY;
   entry->integration_time = st->st_night_integration_time (st);
   entry->num_repeats = num_repeats;

   return entry;
}

typedef struct
{
   double duration;
   double tstart;
   double xstart;
   double ystart;
   int num_steps;
   int num_repeats;
}
Table_Schedule;

static int append_entry (Plan_List_Type **lst, const Scan_Type *st,
                         const Table_Schedule *x)
{
   Plan_List_Type *entry = NULL;

   if (x->num_steps <= 0)
     return 0;

   if (NULL == (entry = plan_list_entry_alloc ()))
     return -1;

   entry->tstart = x->tstart;
   entry->xstart = x->xstart;
   entry->ystart = x->ystart;
   entry->num_steps = x->num_steps;
   entry->scan_duration = x->duration * SEC_PER_DAY;
   entry->integration_time = st->st_integration_time (st);
   entry->num_repeats = x->num_repeats;

   return plan_list_append (lst, entry);
}

static int opt1_scan_table (const Scan_Type *st, const AziElev_Type *beg, const AziElev_Type *end,
                            double *xstart, double *ystart, int *num_steps)
{
   double step_size;
   int max_num_steps;

   /* In this context, we care only about the absolute value of the
    * mirror step size */
   step_size = fabs(st->st_step_size (st));
   if (step_size == 0.0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: invalid mirror step size = %g",
                     __func__, step_size);
        return -1;
     }

   /* Coordinate system should ensure beg->azimuth > end->azimuth */
   max_num_steps = (beg->azimuth - end->azimuth) / step_size;

   if (max_num_steps <= 0)
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: max_num_steps = %d (expected > 0)",
                     __func__, max_num_steps);
        return -1;
     }

   /* optimize sunset/sunrise scans
    * 0=sunrise, 1=mid-day, 2=sunset
    */
   xstart[0] = beg->azimuth;
   xstart[1] = xstart[0];
   xstart[2] = (beg->azimuth + end->azimuth) * 0.5;

   *ystart = beg->elevation;

   num_steps[0] = max_num_steps/2;
   num_steps[1] = max_num_steps;
   num_steps[2] = max_num_steps - num_steps[0];

   return 0;
}

static Plan_List_Type *
opt1_plan (const Scan_Type *st, Solar_Geom_Type *solar_geom,
           const Scan_Limit_Times_Type *limit_times)
{
   Plan_List_Type *opt_scans = NULL;
   Table_Schedule rise={0}, full={0}, set={0};
   AziElev_Type beg={0}, end={0};
   double time_midpoint, sun_angle, min_sun_angle, jd_utc;
   double xstart[3], ystart;
   int num_steps[3];

   if (0 != radiance_scan_endpoints (st, solar_geom, &beg, &end))
     return NULL;
   if (0 != opt1_scan_table (st, &beg, &end, xstart, &ystart, num_steps))
     return NULL;

   rise.xstart = xstart[0];   rise.num_steps = num_steps[0];
   full.xstart = xstart[1];   full.num_steps = num_steps[1];
    set.xstart = xstart[2];    set.num_steps = num_steps[2];

   rise.ystart = ystart;
   full.ystart = ystart;
    set.ystart = ystart;

   rise.duration = st->st_scan_duration (st, rise.num_steps);
   full.duration = st->st_scan_duration (st, full.num_steps);
   set.duration  = st->st_scan_duration (st, set.num_steps);

   min_sun_angle = st->st_min_sun_angle (st);

   time_midpoint = 0.5 * (limit_times->jd_utc_beg + limit_times->jd_utc_end);

   /* number of full scans */
   full.num_repeats = ceil((limit_times->jd_utc_end_full
                            - (limit_times->jd_utc_beg_full - full.duration))
                           / full.duration);

   /* The first full scan actually starts at: */
   full.tstart = time_midpoint - 0.5 * full.num_repeats * full.duration;

   /* Fill out the morning with sunrise scans: */
   rise.tstart = full.tstart;   /* initialization */
   rise.num_repeats = ceil((full.tstart - limit_times->jd_utc_beg)
                          / rise.duration);
   /* impose sun_angle safety constraint */
   while (rise.num_repeats > 0)
     {
        rise.tstart = full.tstart - rise.num_repeats * rise.duration;
        if (0 != solar_geom->sgt_sat_sun_angle (solar_geom, rise.tstart, &sun_angle))
          return NULL;
        if (sun_angle > min_sun_angle)
          break;
        rise.num_repeats -= 1;
     }

   /* Fill out the afternoon with sunset scans: */
   set.tstart = full.tstart + full.num_repeats * full.duration;
   set.num_repeats = ceil((limit_times->jd_utc_end - set.tstart)
                         / set.duration);
   /* impose sun_angle safety constraint */
   while (set.num_repeats > 0)
     {
        double set_tend = set.tstart + set.num_repeats * set.duration;
        if (0 != solar_geom->sgt_sat_sun_angle (solar_geom, set_tend, &sun_angle))
          return NULL;
        if (sun_angle > min_sun_angle)
          break;
        set.num_repeats -= 1;
     }

   /* make sure start times are consistent */
   jd_utc = rise.tstart;
   jd_utc += rise.num_repeats * rise.duration;
   full.tstart = jd_utc;
   jd_utc += full.num_repeats * full.duration;
   set.tstart = jd_utc;

   /* append scan entries */
   if ((0 != append_entry (&opt_scans, st, &rise))
       || (0 != append_entry (&opt_scans, st, &full))
       || (0 != append_entry (&opt_scans, st, &set)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: appending scan plan entries",
                     __func__);
        plan_list_free (opt_scans);
        return NULL;
     }

   return opt_scans;
}

static int scan_vis (Vis_Type *v, const Plan_List_Type *lst, double step_size,
                     int num_days, int ncid)
{
   const Plan_List_Type *entry;
   char var_name[32];
   double *sza = NULL, jd_utc;
   int i, k, status = -1;

#define VAR_NAME_FMT "sza_%02d"

   if (0 != vis_write_grid (v, ncid))
     goto return_status;

   k = 0;
   for (entry = lst; entry != NULL; entry = entry->next)
     {
        for (i = 0; i < entry->num_repeats; i++)
          {
             jd_utc = entry->tstart + i * entry->scan_duration / SEC_PER_DAY;
             if (NULL == (sza = vis_sza (v, jd_utc, sza)))
               goto return_status;
             snprintf (var_name, sizeof(var_name), VAR_NAME_FMT, k);
             k++;
             if (0 != vis_write_value (v, ncid, jd_utc, var_name, sza, step_size, entry))
               goto return_status;
          }
        /* Process only num_days of the plan:
         * (the door-closed part of the orbit is always 1/3 of a day). */
        jd_utc = entry->tstart + entry->num_repeats * entry->scan_duration / SEC_PER_DAY;
        if ((entry->next != NULL)
            && (entry->next->tstart - jd_utc > 0.25))
          {
             num_days--;
             if (num_days == 0)
               break;
          }
     }

#if 0
   if (entry)
     {
        /* end of last scan */
        jd_utc = entry->tstart + entry->num_repeats * entry->scan_duration / SEC_PER_DAY;
        if (NULL == vis_sza (v, jd_utc, sza))
          goto return_status;
        snprintf (var_name, sizeof(var_name), VAR_NAME_FMT, k);
        if (0 != vis_write_value (v, ncid, jd_utc, var_name, sza, step_size, entry))
          goto return_status;
     }
#endif

   status = 0;
return_status:
   FREE(sza);
   return status;
}

typedef struct
{
   const char *name;
   Scan_Method_Type method;
}
Method_Entry;
#define METHOD_ENTRY(name,plan,vis) {name,{plan,vis}}
#define METHOD_TABLE_END  {NULL,{NULL,NULL}}

static Method_Entry Method_Table[] =
{
   METHOD_ENTRY("std", std_plan, scan_vis),
   METHOD_ENTRY("opt1", opt1_plan, scan_vis),
   METHOD_ENTRY("nl_dawn", nightlights_dawn_plan, scan_vis),
   METHOD_ENTRY("nl_dusk", nightlights_dusk_plan, scan_vis),
   METHOD_TABLE_END
};

const Scan_Method_Type *find_scan_method (const char *name)
{
   Method_Entry *e = NULL;

   for (e = Method_Table; e->name != NULL; e++)
     {
        if (0 == strcmp (e->name, name))
          {
             return &e->method;
          }
     }

   return NULL;
}
