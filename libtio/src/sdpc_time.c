#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include "tio.h"
#include "_tio.h"
#include "tio_template.h"

enum
{
   TASK_UNKNOWN,
   TASK_PRINT_TIMES,
   TASK_PRINT_TIMESTAMP,
   TASK_PRINT_SAT_LOCAL_DAY,
   TASK_FIX_FILE
};

/* GPS epoch: 1980-01-06T00:00:00Z */
#define EPOCH_DEFAULT "1980-01-06T00:00:00Z"

#define BUFSIZE 64

static void usage (void)
{
   int sc_timezone_default;
   (void) _pTIO_get_sc_timezone (&sc_timezone_default);
   fprintf (stderr, "Usage: sdpc_time <TAI seconds since epoch>\n");
   fprintf (stderr, "   or: sdpc_time YYYY-MM-DDTHH:MM:SS.SSSZ\n");
   fprintf (stderr, "   or: sdpc_time YYYYMMDDTHHMMSS.SSSZ\n");
   fprintf (stderr, "   or: sdpc_time dDDDDDmMMMMMMMMuUUU\n");
   fprintf (stderr, "   or: sdpc_time -f FILE [-g PATH] [-v VARNAME]\n");
   fprintf (stderr, "Options:\n");
   fprintf (stderr, "  -e | --epoch TSTAMP Epoch defined as an ISO-8601 UTC timestamp string\n");
   fprintf (stderr, "                      [default: %s]\n", EPOCH_DEFAULT);
   fprintf (stderr, "  -z | --zone h       Spacecraft local time-zone offset from UTC. Must be in range [-12,12].\n");
   fprintf (stderr, "                      [default: %d]\n", sc_timezone_default);
   fprintf (stderr, "  -f | --fix FILE     Fix header timestamps in TEMPO netcdf4/HDF5 file\n");
   fprintf (stderr, "  -d | --day FILE     Print satellite-local day number for file\n");
   fprintf (stderr, "  -g | --grp PATH     File group containing time variable [default: /]\n");
   fprintf (stderr, "  -v | --var VARNAME  Name of time variable [default: /time]\n");
   fprintf (stderr, "  -w | --write        Write time_reference timestamp to netcdf4/HDF5 file header\n");
   fprintf (stderr, "  -T | --timestamp    Generate a UTC timestamp for a TEMPO filename\n");
   fprintf (stderr, "  ISO-8601 UTC timestamp format: YYYY-MM-DDTHH:MM:SSZ or YYYYMMDDTHHMMSSZ\n");
   exit (EXIT_SUCCESS);
}

static int fix_header_timestamp (const char *path, const char *grp_path,
                                 const char *var, int write_epoch)
{
   TIO_Var_Info_Type info = {0};
   double tstart, tend;
   int ncid, grp, start, count;

   if (0 != TIO_open (path, NC_WRITE, &ncid))
     return -1;

   if (write_epoch)
     {
        if (0 != tio_write_epoch_timestamp (ncid, NC_GLOBAL))
          return -1;
     }

   if (0 != TIO_inq_grp (ncid, grp_path, &grp))
     return -1;

   if (0 != TIO_inq_var (grp, var, &info))
     return -1;

   start = 0;
   count = 1;
   if (0 != TIO_get_var_section (grp, var, &start, &count, NC_DOUBLE, &tstart))
     return -1;

   start = info.dimlens[0]-1;
   count = 1;
   if (0 != TIO_get_var_section (grp, var, &start, &count, NC_DOUBLE, &tend))
     return -1;

   if ((0 != TIO_write_timestamp (ncid, NC_GLOBAL, "time_coverage_start", tstart))
       ||(0 != TIO_write_timestamp (ncid, NC_GLOBAL, "time_coverage_end", tend)))
     return -1;

   if (0 != TIO_close (ncid))
     return -1;

   return 0;
}

static int print_sat_local_day (const char *path)
{
   int ncid, status = -1;
   double tstart, flocal_day;

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     return -1;

   if (0 != tio_use_file_epoch (ncid))
     goto return_status;

   if (0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &tstart))
     goto return_status;

   if (0 != tio_time_sat_local_day_number (tstart, &flocal_day))
     goto return_status;

   fprintf (stdout, "%d\n", (int) flocal_day);

   status = 0;
return_status:
   TIO_close (ncid);
   return status;
}

static int print_ioc_string (double taix_sec)
{
   int day, msec, usec, sec_per_day = 86400;
   double f_msec;

   day = taix_sec / sec_per_day;
   f_msec = (taix_sec - day * sec_per_day) * 1000;

   msec = f_msec;
   usec = (f_msec - msec) * 1000;

   fprintf (stdout, "IOC: d%05dm%08du%03d\n", day, msec, usec);

   return 0;
}

static int make_iso8601_string (double taix, int want_delim, int want_frac,
                                char *buf, int bufsize)
{
   double hour, minf, sec;
   int year, month, day, hr, min;
   const char *dash;
   const char *colon;

   if (0 != tio_time_taix_to_utc_caldate (taix, &year, &month, &day, &hour))
     return -1;

   hr   = (int)hour;
   minf = (hour - hr)*60;
   min = (int)minf;
   sec = (minf - min)*60;

   /* Truncate numerical value to match output precision.
    * If we don't do this, rounding on output may show
    * a string timestamp with seconds=60.
    */
   sec = ((int)(sec * 1.e6))/1.e6;

   dash  = want_delim ? "-" : "";
   colon = want_delim ? ":" : "";

   if (want_frac)
     {
        return snprintf (buf, bufsize, "%4d%s%02d%s%02dT%02d%s%02d%s%09.6fZ",
                         year, dash, month, dash, day,
                         hr, colon, min, colon, sec);
     }
   else
     {
        return snprintf (buf, bufsize, "%4d%s%02d%s%02dT%02d%s%02d%s%02dZ",
                         year, dash, month, dash, day,
                         hr, colon, min, colon, (int) sec);
     }
}

static int print_timestamp (double taix)
{
   char buf[BUFSIZE];
   if (make_iso8601_string (taix, 0, 0, buf, sizeof(buf)) < 0)
     return -1;
   fprintf (stdout, "%s", buf);
   return 0;
}

static int print_taix_as_strings (double taix)
{
   double tai, utc;
   time_t tt;
   struct tm tm = {0};
   char buf[BUFSIZE];

   if (make_iso8601_string (taix, 1, 1, buf, sizeof(buf)) < 0)
     return -1;
   fprintf (stdout, "UTC: %s\n", buf);

   if (0 != tio_time_taix_to_utc (taix, &utc))
     return -1;

   tt = (time_t) utc;
   gmtime_r (&tt, &tm);
   fprintf (stdout, "     time_t: %ld    day of year: %d (1..366)\n",
            tt, 1+tm.tm_yday);

   /* taix is seconds since TEMPO epoch,
    * tai is seconds since Unix epoch
    */
   if (-1 == tio_time_taix_to_tai (taix, &tai))
     return -1;
   tt = (time_t)tai;
   gmtime_r (&tt, &tm);
   strftime (buf, BUFSIZE, "%Y-%m-%dT%H:%M:%S", &tm);
   fprintf (stdout, "TAI: %s.%06d+00\n", buf, (int)(round((tai-tt)*1e6)));

   return 0;
}

static int convert_utc_string_to_taix (const char *arg, double *ptaix)
{
   double taix, fsec = 0.0;
   char *dot;
   size_t len;
   char buf[BUFSIZE];

   if (arg == NULL)
     {
        fprintf (stderr, "%s: NULL string\n", __func__);
        return -1;
     }

   dot = strchr (arg, '.');

   if (dot)
     {
        len = dot - arg;
        if (1 != sscanf (dot, "%le", &fsec))
          goto error_return;
     }
   else len = strlen (arg);

   if (len >= BUFSIZE)
     goto error_return;
   memset (buf, 0, sizeof(buf));
   strncpy (buf, arg, len);

   if (buf[len-1] != 'Z')
     {
        if (len+1 >= BUFSIZE)
          goto error_return;
        buf[len] = 'Z';
        buf[len+1] = 0;
     }

   if (0 != tio_time_utcstr_to_taix (buf, &taix))
     return -1;
   taix += fsec;

   if (ptaix) *ptaix = taix;

   return 0;

error_return:
   fprintf (stderr, "*** ERROR: converting timestamp %s\n", arg);
   return -1;
}

static int convert_ioc_string_to_taix (const char *str, double *ptaix)
{
   int day, msec, usec;
   double taix;

   if (str == NULL)
     {
        fprintf (stderr, "%s: NULL string\n", __func__);
        return -1;
     }

   if (3 != sscanf (str, "d%5dm%8du%3d", &day, &msec, &usec))
     {
        fprintf (stderr, "*** Error: parsing timestamp: %s\n", str);
        return -1;
     }

   taix = day * 86400.0 + msec/1000.0 + usec/1.e6;
   if (ptaix) *ptaix = taix;

   return 0;
}

int main (int argc, char **argv)
{
   static struct option long_options[] =
     {
        {"day",   required_argument, 0, 'd'},
        {"write", no_argument,       0, 'w'},
        {"epoch", required_argument, 0, 'e'},
        {"zone",  required_argument, 0, 'z'},
        {"fix",   required_argument, 0, 'f'},
        {"grp",   required_argument, 0, 'g'},
        {"var",   required_argument, 0, 'v'},
        {"timestamp", no_argument,   0, 'T'},
        {0,0,0,0}
     };
   const char *timestamp_string = NULL;
   const char *epoch_string = EPOCH_DEFAULT;
   const char *path = NULL;
   const char *grp = "/";
   const char *var = "time";
   double taix = 0.0;
   double flocal_day;
   int exit_status = EXIT_FAILURE;
   int status = -1;
   int write_epoch = 0;
   int task = TASK_UNKNOWN;
   int sc_timezone = INT_MAX;
   int utc_day;
   int have_utc_string = 0;
   int have_ioc_string = 0;

   if (argc < 2)
     usage();

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "wd:e:f:g:Tv:z:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d??\n", c);
             break;
           case 'w':
             write_epoch++;
             break;
           case 'g':
             grp = optarg;
             break;
           case 'd':
             path = optarg;
	     task = TASK_PRINT_SAT_LOCAL_DAY;
             break;
           case 'f':
             path = optarg;
	     task = TASK_FIX_FILE;
             break;
           case 'v':
             var = optarg;
             break;
	   case 'e':
	     epoch_string = optarg;
	     break;
           case 'z':
	     if (1 != sscanf (optarg, "%d", &sc_timezone))
	       {
		  fprintf (stderr, "*** Error: setting spacecraft time zone\n");
		  exit(1);
	       }
             break;
           case 'T':
	     task = TASK_PRINT_TIMESTAMP;
             break;
          }
     }

   if ((optind == argc) && (task != TASK_FIX_FILE))
     usage();
   else
     {
        int len;

        timestamp_string = argv[optind++];
        len = strlen (timestamp_string);

        if (timestamp_string[len-1] == 'Z')
          have_utc_string = 1;
        else if (timestamp_string[0] == 'd')
          have_ioc_string = 1;
        else if (1 != sscanf (timestamp_string, "%le", &taix))
          {
             fprintf (stderr, "*** Error: invalid timestamp string: %s\n", timestamp_string);
             exit(1);
          }
        if (task == TASK_UNKNOWN) task = TASK_PRINT_TIMES;
     }

   if (sc_timezone != INT_MAX)
     {
        /* timezone set on command line */
        if (0 != _pTIO_set_sc_timezone (sc_timezone))
          goto error_return;
     }
   else if (0 != _pTIO_get_sc_timezone (&sc_timezone))
     {
        goto error_return;
     }

   if (task == TASK_PRINT_SAT_LOCAL_DAY)
     {
        status = print_sat_local_day (path);
        goto error_return;
     }

   if (0 != tio_time_set_taix_epoch (epoch_string))
     goto error_return;

   if (have_utc_string)
     {
        if (0 != convert_utc_string_to_taix (timestamp_string, &taix))
          goto error_return;
     }
   else if (have_ioc_string)
     {
        if (0 != convert_ioc_string_to_taix (timestamp_string, &taix))
          goto error_return;
     }

   switch (task)
     {
      case TASK_PRINT_TIMESTAMP:
        status = print_timestamp (taix);
        break;

      case TASK_PRINT_TIMES:
        fprintf (stdout, "SEC: %0.6f\n", taix);
	status = print_taix_as_strings (taix);
        print_ioc_string (taix);

        utc_day = taix / 86400.0;
        if (0 != tio_time_sat_local_day_number (taix, &flocal_day))
          goto error_return;

        fprintf (stdout, "DAY: %d UTC\n", utc_day);
        fprintf (stdout, "DAY: %f local at UTC%+03d\n", flocal_day, sc_timezone);
	break;

      case TASK_FIX_FILE:
        status = fix_header_timestamp (path, grp, var, write_epoch);
	break;

      default:
	fprintf (stderr, "*** Error: unsupported task\n");
	break;
     }

error_return:
   exit_status = status ? EXIT_FAILURE : EXIT_SUCCESS;

   return exit_status;
}
