#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>

#include "tio.h"
#include "_tio.h"
#include "tio_template.h"

enum
{
   TASK_UNKNOWN = 0,
   TASK_CONVERT_TAI = 1,
   TASK_CONVERT_UTC_STRING = 2,
   TASK_CONVERT_IOC_STRING = 3,
   TASK_FIX_FILE = 4
};

#define EPOCH_DEFAULT "2000-01-01T00:00:00Z"

static void usage (void)
{
   fprintf (stderr, "Usage: ttime -s SECONDS\n");
   fprintf (stderr, "   or: ttime -u YYYY-MM-DDTHH:MM:SSZ\n");
   fprintf (stderr, "   or: ttime -i dDDDDDmMMMMMMMMuUUU\n");
   fprintf (stderr, "   or: ttime -f FILE [-g PATH] [-v VARNAME]\n");
   fprintf (stderr, "Options:\n");
   fprintf (stderr, "  -e | --epoch TSTAMP Epoch defined as an ISO-8601 UTC timestamp string\n");
   fprintf (stderr, "                      [default: %s]\n", EPOCH_DEFAULT);
   fprintf (stderr, "  -f | --fix FILE     Fix header timestamps in TEMPO netcdf4/HDF5 file\n");
   fprintf (stderr, "  -g | --grp PATH     File group containing time variable [default: /]\n");
   fprintf (stderr, "  -v | --var VARNAME  Name of time variable [default: /time]\n");
   fprintf (stderr, "  -w | --write        Write time_reference timestamp to netcdf4/HDF5 file header\n");
   fprintf (stderr, "  -i | --ioc TSTAMP   Convert IOC timestamp string to TAI seconds since epoch\n");
   fprintf (stderr, "  -u | --utc TSTAMP   Convert UTC timestamp string to TAI seconds since epoch\n");
   fprintf (stderr, "  -s | --sec SECONDS  Convert TAI seconds since epoch to UTC timestamp string\n");
   fprintf (stderr, "  -d | --delim        Output UTC timestamp string omitting delimiters :-\n");
   fprintf (stderr, "  ISO-8601 UTC timestamp format: YYYY-MM-DDTHH:MM:SSZ\n");
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

static int print_ioc_string (double tai_sec)
{
   int day, msec, usec, sec_per_day = 86400;
   double f_msec;

   day = tai_sec / sec_per_day;
   f_msec = (tai_sec - day * sec_per_day) * 1000;

   msec = f_msec;
   usec = (f_msec - msec) * 1000;

   fprintf (stdout, "d%05dm%08du%03d\n", day, msec, usec);

   return 0;
}

static int convert_tai_to_utc_string (double elapsed_seconds, int omit_delimiters)
{
   double hour, minf, sec;
   int year, month, day, hr, min;

   if (0 != tio_time_taix_to_utc_caldate (elapsed_seconds, &year, &month, &day, &hour))
     return -1;

   hr   = (int)hour;
   minf = (hour - hr)*60;
   min = (int)minf;
   sec = (minf - min)*60;

   if (omit_delimiters)
     {
        fprintf (stdout, "%4d%02d%02dT%02d%02d%02.0f\n",
                 year, month, day, hr, min, sec);
     }
   else
     {
        fprintf (stdout, "%4d-%02d-%02dT%02d:%02d:%09.6fZ\n",
                 year, month, day, hr, min, sec);
     }

   return 0;
}

static int convert_utc_string_to_tai (const char *arg, double *ptai)
{
   double tempo, fsec = 0.0;
   char *dot;
   size_t len;
#define BUFSIZE 32
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

   if (0 != tio_time_utcstr_to_taix (buf, &tempo))
     return -1;
   tempo += fsec;

   if (ptai) *ptai = tempo;

   fprintf (stdout, "%0.6f\n", tempo);
   return 0;

error_return:
   fprintf (stderr, "*** ERROR: converting timestamp %s\n", arg);
   return -1;
}

static int convert_ioc_string_to_tai (const char *str, double *ptai)
{
   int day, msec, usec;
   double tai;

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

   tai = day * 86400.0 + msec/1000.0 + usec/1.e6;
   if (ptai) *ptai = tai;

   fprintf (stdout, "%0.6f\n", tai);

   return 0;
}

int main (int argc, char **argv)
{
   static struct option long_options[] =
     {
        {"utc",   required_argument, 0, 'u'},
        {"sec",   required_argument, 0, 's'},
        {"delim", no_argument,       0, 'd'},
        {"write", no_argument,       0, 'w'},
        {"epoch", required_argument, 0, 'e'},
        {"fix",   required_argument, 0, 'f'},
        {"grp",   required_argument, 0, 'g'},
        {"var",   required_argument, 0, 'v'},
        {0,0,0,0}
     };
   const char *utc_string = NULL;
   const char *ioc_string = NULL;
   const char *epoch_string = EPOCH_DEFAULT;
   const char *path = NULL;
   const char *grp = "/";
   const char *var = "time";
   double elapsed_seconds = 0.0;
   double tai;
   int exit_status = EXIT_FAILURE;
   int status = -1;
   int write_epoch = 0;
   int omit_delimiters = 0;
   int task = TASK_UNKNOWN;

   if (argc < 3)
     usage();

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "wde:f:g:i:s:u:v:", long_options, &option_index);
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
           case 'd':
             omit_delimiters++;
             break;
           case 'g':
             grp = optarg;
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
           case 's':
	     task = TASK_CONVERT_TAI;
	     if (1 != sscanf (optarg, "%le", &elapsed_seconds))
	       {
		  fprintf (stderr, "*** Error: converting %s to elapsed seconds since the epoch\n",
			   optarg ? optarg : "<null>");
		  exit(1);
	       }
             break;
           case 'u':
	     task = TASK_CONVERT_UTC_STRING;
	     utc_string = optarg;
             break;
           case 'i':
	     task = TASK_CONVERT_IOC_STRING;
	     ioc_string = optarg;
             break;
          }
     }

   if (0 != tio_time_set_taix_epoch (epoch_string))
     goto error_return;

   switch (task)
     {
      case TASK_CONVERT_TAI:
	status = convert_tai_to_utc_string (elapsed_seconds, omit_delimiters);
        print_ioc_string (elapsed_seconds);
	break;

      case TASK_CONVERT_UTC_STRING:
	status = convert_utc_string_to_tai (utc_string, &tai);
        print_ioc_string (tai);
	break;

      case TASK_CONVERT_IOC_STRING:
        if ((0 == convert_ioc_string_to_tai (ioc_string, &tai))
            && (0 == convert_tai_to_utc_string (tai, omit_delimiters)))
          {
             status = 0;
          }
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
