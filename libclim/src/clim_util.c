#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <wordexp.h>

#include <ioclib.h>
#include <tell.h>

typedef struct
{
   char *cf_climatology;
   char *cf_forecast;
   char *cf_pressure_grid;
   char *cloud_climatology;
}
Clim_File_Pattern_Type;
static Clim_File_Pattern_Type Clim_File_Patterns = {0};

#define FILE_PATTERN_PTR(s) (Clim_File_Patterns.s)

extern int read_config_file (const char *config_file);
extern int make_climatology_path (int month, int hour, char *buf, int bufsize);
extern int make_forecast_path (time_t tt, char *buf, int bufsize);
extern int make_pressure_eta_path (char *buf, int bufsize);
extern int make_cloud_climatology_path (char *buf, int bufsize);

extern int have_forecast_files (time_t tt, int num_hours);
extern long make_timet (int year, int month, int day, int hour);

static char *expand_string (const char *s, int quiet)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        if (!quiet) tell_verror (TELL_UNKNOWN_ERROR,
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

   if (NULL == (fmt = expand_string (buf, 0)))
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

static int replace_alloc_str (char **targ, char *src)
{
   ioclib_free (*targ);
   *targ = ioclib_strdup (src);
   return *targ ? 0 : -1;
}
#define REPLACE_ALLOC_STR(old,p,field) \
   do {if (replace_alloc_str (&old->field, p->field)) return -1; } while (0)

static int replace_file_patterns (Clim_File_Pattern_Type *old, Clim_File_Pattern_Type *p)
{
   if ((old == NULL) || (p == NULL))
     return -1;
   REPLACE_ALLOC_STR(old,p,cf_climatology);
   REPLACE_ALLOC_STR(old,p,cf_forecast);
   REPLACE_ALLOC_STR(old,p,cf_pressure_grid);
   REPLACE_ALLOC_STR(old,p,cloud_climatology);
   return 0;
}

int read_config_file (const char *config_file)
{
   Clim_File_Pattern_Type p = {0};
   IOCLib_String_Array_Obj_Type *sa = NULL;
   IOCLib_Config_String_Type tbl[] =
     {
        {"GEOSCF_Climatology_Files", &p.cf_climatology, IOCLIB_CONFIG_TYPE_STR},
        {"GEOSCF_Forecast_Files", &p.cf_forecast, IOCLIB_CONFIG_TYPE_STR},
        {"GEOSCF_Pressure_Grid", &p.cf_pressure_grid, IOCLIB_CONFIG_TYPE_STR},
        {"Cloud_Climatology", &p.cloud_climatology, IOCLIB_CONFIG_TYPE_STR},
        {NULL,NULL,0}
     };
   char *cfg_file = NULL;
   int status = -1;

   /* When a config file path is not provided, search
    * for the path in a few likely places, starting with
    * the closest, and expanding outward.
    */

   if ((config_file == NULL)
       && (NULL == (config_file = getenv ("SDPC_GEOSCF_CONFIG"))))   /* check local environment customization */
     {
        char *conf_file_paths[] =
          {
             "$SDPC_PIPE_DIR/etc/table_config.ini",  /* running pipeline context */
             "$SDPC_ROOT/etc/table_config.ini",            /* software installation context */
             NULL
          };
        char **cfp;
         for (cfp = conf_file_paths; *cfp != NULL; cfp++)
          {
             if (NULL != (cfg_file = expand_string (*cfp, 1)))
               {
                  if (0 == access (cfg_file, F_OK | R_OK))
                    break;
               }
          }
        if (cfg_file == NULL)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: config_file not specified, default file not found", __func__);
             return -1;
           }
        config_file = cfg_file;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: reading config file: %s", __func__, config_file);

   if (NULL == (sa = ioclib_config_get_strings (config_file, "GEOSCF", tbl)))
     goto return_status;

   replace_file_patterns (&Clim_File_Patterns, &p);

   status = 0;
return_status:
   ioclib_string_array_obj_free (sa);
   if (cfg_file) free(cfg_file);
   return status;
}

int make_climatology_path (int month, int hour, char *buf, int bufsize)
{
   const char *fmt = FILE_PATTERN_PTR(cf_climatology);
   int n;

   if (fmt == NULL)
     {
        if (read_config_file (NULL))
          return -1;
        fmt = FILE_PATTERN_PTR(cf_climatology);
     }

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
   const char *fmt = FILE_PATTERN_PTR(cf_forecast);
   struct tm tm = {0};
   struct tm rpl_tm = {0};
   time_t rpl_tt;
   int rpl_year, rpl_month, rpl_day;
   int i, n, year, month, day, dayofyear;

   if (fmt == NULL)
     {
        if (read_config_file (NULL))
          return -1;
        fmt = FILE_PATTERN_PTR(cf_forecast);
     }

   if (NULL == gmtime_r (&tt, &tm))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gmtime_r failed (tt = %ld", __func__, tt);
        return -1;
     }

   year = tm.tm_year + 1900;
   month = tm.tm_mon + 1;
   day = tm.tm_mday;
   dayofyear = tm.tm_yday + 1;

   /* rpl_tt is the date of the "replay" upon which the forecast is based.
    * We prefer forecasts based on the previous day's "replay", but
    * sometimes we must settle for a two day old forecast */
   rpl_tt = tt;

   for (i = 0; i < 2; i++)
     {
        rpl_tt -= 86400;
        if (NULL == gmtime_r (&rpl_tt, &rpl_tm))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: gmtime_r failed (rpl_tt = %ld", __func__, rpl_tt);
             return -1;
          }

        rpl_year = rpl_tm.tm_year + 1900;
        rpl_month = rpl_tm.tm_mon + 1;
        rpl_day = rpl_tm.tm_mday;

        if ((n = snprintf (buf, bufsize, fmt, year, dayofyear,
                           rpl_year, rpl_month, rpl_day,
                           year, month, day, tm.tm_hour)) < 0)
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

        if (0 != expand_buffer_in_place (buf, bufsize))
          return -1;

        /* If the preferred file doesn't exist, we'll look for an older forecast */
        if (0 == access (buf, F_OK | R_OK))
          return 0;
     }

   /* As long as the path was created successfully, we return "success",
    * even when the generated path doesn't exist. */
   return 0;
}

/* return 1 means we have the forecast file,
 * return 0 means we don't have the forecast file
 * return -1 means an error occurred
 */
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
   const char *fmt = FILE_PATTERN_PTR(cf_pressure_grid);
   int n;

   if (fmt == NULL)
     {
        if (read_config_file (NULL))
          return -1;
        fmt = FILE_PATTERN_PTR(cf_pressure_grid);
     }

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
   const char *fmt = FILE_PATTERN_PTR(cloud_climatology);
   int n;

   if (fmt == NULL)
     {
        if (read_config_file (NULL))
          return -1;
        fmt = FILE_PATTERN_PTR(cloud_climatology);
     }

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
