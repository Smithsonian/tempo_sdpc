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

static int scan_table_params (const Scan_Type *st,
                              const Solar_Geom_Type *solar_geom,
                              int num_tables, double *xstart, int *num_steps)
{
   EarthPoint beg_pt={0}, end_pt={0};
   AziElev_Type beg, end;
   TempoGeoErr error;
   double sat_lon, step_size;
   int max_num_steps;

   if (0 != solar_geom->sgt_geosat_longitude(solar_geom, &sat_lon))
     return -1;
   if (0 != st->st_scan_beg (st, &beg_pt.theLon, &beg_pt.theLat))
     return -1;
   if (0 != st->st_scan_end (st, &end_pt.theLon, &end_pt.theLat))
     return -1;

   sat_lon /= DEGTORAD;

   if ((error = computeScanAngles (&beg_pt, sat_lon, SCAN_AZ_FIRST,
                                   TEMPO_FIRST_CORRECTION,
                                   &beg.azimuth, &beg.elevation)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computeScanAngles failed (%s)",
                     __func__, tempoGeoErrorString (error));
        return -1;
     }

   if ((error = computeScanAngles (&end_pt, sat_lon, SCAN_AZ_FIRST,
                                   TEMPO_FIRST_CORRECTION,
                                   &end.azimuth, &end.elevation)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computeScanAngles failed (%s)",
                     __func__, tempoGeoErrorString (error));
        return -1;
     }

   /* want xstart in microradians */

   beg.azimuth *= DEGTOMICRORAD;   beg.elevation *= DEGTOMICRORAD;
   end.azimuth *= DEGTOMICRORAD;   end.elevation *= DEGTOMICRORAD;

   if ((step_size = st->st_step_size (st)) <= 0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: invalid mirror step size = %g",
                     __func__, step_size);
        return -1;
     }

   /* Coordinate system should ensure beg.azimuth > end.azimuth */
   max_num_steps = (beg.azimuth - end.azimuth) / step_size;

   if (max_num_steps <= 0)
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: max_num_steps = %d (expected > 0)",
                     __func__, max_num_steps);
        return -1;
     }

   switch (num_tables)
     {
      case 1:
        /* baseline */
        *xstart = beg.azimuth;
        *num_steps = max_num_steps;
        break;

      case 3:
        /* optimize sunset/sunrise scans
         * 0=sunrise, 1=mid-day, 2=sunset
         */
        xstart[0] = beg.azimuth;
        xstart[1] = xstart[0];
        xstart[2] = (beg.azimuth + end.azimuth) * 0.5;

        num_steps[0] = max_num_steps/2;
        num_steps[1] = max_num_steps;
        num_steps[2] = max_num_steps - num_steps[0];
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: num_tables=%d not supported",
                     __func__, num_tables);
        return -1;
     }

   return 0;
}

static Plan_List_Type *
std_plan (const Scan_Type *st, Solar_Geom_Type *solar_geom,
          const Scan_Limit_Times_Type *limit_times)
{
   Plan_List_Type *entry = NULL;
   double time_full_scan, xstart;
   int num_steps;

   if (0 != scan_table_params (st, solar_geom, 1, &xstart, &num_steps))
     return NULL;

   if (NULL == (entry = plan_list_entry_alloc ()))
     return NULL;

   time_full_scan = st->st_scan_duration (st, num_steps);

   entry->tstart = limit_times->jd_utc_beg;
   entry->xstart = xstart;
   entry->num_steps = num_steps;
   entry->scan_duration = time_full_scan * SEC_PER_DAY;
   entry->integration_time = st->st_integration_time (st);
   entry->num_repeats = ceil ((limit_times->jd_utc_end - limit_times->jd_utc_beg)
                              / time_full_scan);
   return entry;
}

static int std_vis (Vis_Type *v, const Plan_List_Type *lst, int ncid)
{
   double *sza = NULL, jd_utc;
   int status = -1;

   if (0 != vis_write_grid (v, ncid))
     goto return_status;

   jd_utc = lst->tstart;
   if (NULL == (sza = vis_sza (v, jd_utc, NULL)))
     goto return_status;
   if (0 != vis_write_value (v, ncid, jd_utc, "sza_beg", sza))
     goto return_status;

   jd_utc = (lst->tstart + lst->num_repeats * lst->scan_duration / SEC_PER_DAY);
   if (NULL == vis_sza (v, jd_utc, sza))
     goto return_status;
   if (0 != vis_write_value (v, ncid, jd_utc, "sza_end", sza))
     goto return_status;

   status = 0;
return_status:
   FREE(sza);
   return status;
}

typedef struct
{
   double duration;
   double tstart;
   double xstart;
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
   entry->num_steps = x->num_steps;
   entry->scan_duration = x->duration * SEC_PER_DAY;
   entry->integration_time = st->st_integration_time (st);
   entry->num_repeats = x->num_repeats;

   return plan_list_append (lst, entry);
}

#define NUM_TABLES_OPT1 3

static Plan_List_Type *
opt1_plan (const Scan_Type *st, Solar_Geom_Type *solar_geom,
           const Scan_Limit_Times_Type *limit_times)
{
   Plan_List_Type *opt_scans = NULL;
   Table_Schedule rise={0}, full={0}, set={0};
   double time_midpoint, sun_angle, min_sun_angle, jd_utc;
   double xstart[NUM_TABLES_OPT1];
   int num_steps[NUM_TABLES_OPT1];

   if (0 != scan_table_params (st, solar_geom, NUM_TABLES_OPT1,
                               xstart, num_steps))
     return NULL;

   rise.xstart = xstart[0];   rise.num_steps = num_steps[0];
   full.xstart = xstart[1];   full.num_steps = num_steps[1];
    set.xstart = xstart[2];    set.num_steps = num_steps[2];

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

static int opt1_vis (Vis_Type *v, const Plan_List_Type *lst, int ncid)
{
   const Plan_List_Type *entry;
   const char *variable_names[] =
     {"sza_beg", "sza_full_beg", "sza_full_end", NULL};
   const char **name;
   double *sza = NULL, jd_utc;
   int status = -1;

   if (0 != vis_write_grid (v, ncid))
     goto return_status;

   entry = lst;
   for (name = variable_names; *name != NULL; name++)
     {
        jd_utc = entry->tstart;
        if (NULL == (sza = vis_sza (v, jd_utc, sza)))
          goto return_status;
        if (0 != vis_write_value (v, ncid, jd_utc, *name, sza))
          goto return_status;

        entry = entry->next;
        if (entry == NULL)
          {
             /* unexpectedly short list, but don't complain */
             status = 0;
             goto return_status;
          }
     }

   entry = lst->next->next;

   /* end of last scan */
   jd_utc = entry->tstart + entry->num_repeats * entry->scan_duration / SEC_PER_DAY;
   if (NULL == vis_sza (v, jd_utc, sza))
     goto return_status;
   if (0 != vis_write_value (v, ncid, jd_utc, "sza_end", sza))
     goto return_status;

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
   METHOD_ENTRY("std", std_plan, std_vis),
   METHOD_ENTRY("opt1", opt1_plan, opt1_vis),
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
