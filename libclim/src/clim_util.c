#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <wordexp.h>

#include <tell.h>

#define CLIM_PATH_PATTERN_CLIMATOLOGY \
   "$SDPC_REFDATA_DIR/climatology/monthly_by_hour/GEOS-CF_NRT.month%02d_%02dz_reorder.nc4"

#define CLIM_PATH_PATTERN_FORECAST \
   "$SDPC_ANCILLARY_ROOT/geos_cf/%4d/%d/GEOS-CF.v01.rpl.sat_inst_1hr_r720x361_v72.%4d%02d%02d_%02d00z_reorder.nc4"

#define CLIM_PATH_PATTERN_PRESSURE_ETA \
   "$SDPC_REFDATA_DIR/climatology/GEOS-Chem_72_layer_vertical_grid.nc"

#define CLIM_PATH_PATTERN_CLOUD \
   "$SDPC_REFDATA_DIR/climatology/omcldrr_pressure.nc"

extern int make_climatology_path (int month, int hour, char *buf, int bufsize);
extern long make_timet (int year, int month, int day, int hour);
extern int make_forecast_path (time_t tt, char *buf, int bufsize);
extern int have_forecast_files (time_t tt, int num_hours);
extern int make_pressure_eta_path (char *buf, int bufsize);
extern int make_cloud_climatology_path (char *buf, int bufsize);

static char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
     }

   return s_exp;
}

static int expand_buffer_in_place (char *buf, size_t bufsize)
{
   char *fmt = NULL;
   size_t n;

   if (NULL == (fmt = expand_string (buf)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: expand_string failed", __func__);
        return -1;
     }

   n = strlen(fmt) + 1;

   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %ld > %ld",
                     __func__, n, bufsize);
        free (fmt);
        return -1;
     }

   strncpy (buf, fmt, bufsize);
   free(fmt);
   return 0;
}

int make_climatology_path (int month, int hour, char *buf, int bufsize)
{
   const char fmt[] = CLIM_PATH_PATTERN_CLIMATOLOGY;
   int n;

   if ((n = snprintf (buf, bufsize, fmt, month, hour)) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        return -1;
     }
   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %d > %d",
                     __func__, n, bufsize);
        return -1;
     }

   return expand_buffer_in_place (buf, bufsize);
}

long make_timet (int year, int month, int day, int hour)
{
   struct tm tm = {0};

   tm.tm_year = year - 1900;
   tm.tm_mon = month - 1;
   tm.tm_mday = day;
   tm.tm_hour = hour;

   return (long) timegm (&tm);
}

int make_forecast_path (time_t tt, char *buf, int bufsize)
{
   const char fmt[] = CLIM_PATH_PATTERN_FORECAST;
   struct tm tm = {0};
   int n, year, dayofyear, month;

   if (NULL == gmtime_r (&tt, &tm))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gmtime_r failed (tt = %ld", __func__, tt);
        return -1;
     }

   year = tm.tm_year + 1900;
   dayofyear = tm.tm_yday + 1;
   month = tm.tm_mon + 1;

   if ((n = snprintf (buf, bufsize, fmt, year, dayofyear, year, month, tm.tm_mday, tm.tm_hour)) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        return -1;
     }
   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %d > %d",
                     __func__, n, bufsize);
        return -1;
     }

   return expand_buffer_in_place (buf, bufsize);
}

int have_forecast_files (time_t tt, int num_hours)
{
#define CF_BUFSIZE 1024
   char path[CF_BUFSIZE];
   struct tm tm = {0};
   int hour;

   if ((num_hours < 0) || (num_hours > 24))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid num_hours = %d", __func__, num_hours);
        return -1;
     }
   if (tt < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid start time = %ld (time_t)", __func__, tt);
        return -1;
     }

   if (NULL == gmtime_r (&tt, &tm))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gmtime_r failed (t=%ld)", __func__, tt);
        return -1;
     }

   for (hour = 0; hour <= num_hours; hour++)
     {
        if (0 != make_forecast_path (tt + hour * 3600, path, CF_BUFSIZE))
          return -1;

        if (0 != access (path, F_OK))
          return 0;
     }

   return 1;
}

int make_pressure_eta_path (char *buf, int bufsize)
{
   const char fmt[] = CLIM_PATH_PATTERN_PRESSURE_ETA;
   int n;

   if ((n = snprintf (buf, bufsize, fmt)) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        return -1;
     }
   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %d > %d",
                     __func__, n, bufsize);
        return -1;
     }

   return expand_buffer_in_place (buf, bufsize);
}

int make_cloud_climatology_path (char *buf, int bufsize)
{
   const char fmt[] = CLIM_PATH_PATTERN_CLOUD;
   int n;

   if ((n = snprintf (buf, bufsize, fmt)) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        return -1;
     }
   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %d > %d",
                     __func__, n, bufsize);
        return -1;
     }

   return expand_buffer_in_place (buf, bufsize);
}
