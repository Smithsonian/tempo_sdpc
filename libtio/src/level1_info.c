#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define __USE_XOPEN
#include <time.h>

#include <limits.h>
#include <stddef.h>
#include <unistd.h>
#include <getopt.h>

#include <netcdf.h>
#include <tell.h>

#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

#define MAX_PRODUCT_TYPE_LEN 16

enum
{
   TASK_UNKNOWN,
   TASK_PRINT_DIR,
   TASK_PRINT_MONTH,
   TASK_PRINT_SATDAY_DIR,
   TASK_PRINT_SCANID
};

static void usage (int argc, char **argv)
{
   (void) argc;
   fprintf (stderr, "Usage: %s <option> <level1-product-file>\n", argv[0]);
   fprintf (stderr, "Options:\n");
   fprintf (stderr, "  -d | --dir       Print the archive sub-directory path for this file\n");
   fprintf (stderr, "  -l | --localday  Print satellite-local day subdirectory name\n");
   fprintf (stderr, "  -m | --month     Print the calendar month when this file started\n");
   fprintf (stderr, "  -s | --scanid    Print the unique scan id number for this file\n");
}

static int print_sat_local_day_number (const _pTIO_Granule_Ident_Type *gid)
{
   double sat_day;

   if (0 != tio_time_sat_local_day_number (gid->tstart, &sat_day))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computing satellite-local day number", __func__);
        return -1;
     }

   if (fprintf (stdout, "D%05d", (int) sat_day) < 0)
     return -1;

   return 0;
}

static int print_archive_subdir (const _pTIO_Granule_Ident_Type *gid, const char *product_type)
{
   double sat_day;

   if (0 != tio_time_sat_local_day_number (gid->tstart, &sat_day))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computing satellite-local day number", __func__);
        return -1;
     }

   /* The archive is organized by major Level 1 product types:
    *    DRK, DRKL
    *    IRR, IRRL
    *    RAD, RADT
    * "DRK" and "IRR" products are organized by sat_day/tstart.
    * "RAD" products are organized by sat_day/scan/granule index.
    */

   if (gid->granule_num > 0)
     {
        /* All Level 1 radiance file variants should have RAD in the product_type.
         * If that's missing, then this file type is not supported.
         */
        if (NULL == strstr (product_type, "RAD"))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported product type: %s", __func__, product_type);
             return -1;
          }
        if (fprintf (stdout, "%s/D%05d/S%03d/G%02d", product_type,
                     (int) sat_day, gid->scan_num, gid->granule_num) < 0)
          return -1;
     }
   else  /* All other Level1 file types are organized by product_type/date/tstart */
     {
        char buf[MAX_ISOTIME_LEN];

        if (0 != TIO_mktimestamp_str (gid->tstart, 0, buf, sizeof(buf)))
          return -1;

        if (fprintf (stdout, "%s/D%05d/%s", product_type,
                     (int) sat_day, buf) < 0)
          return -1;
     }

   return 0;
}

static int print_scan_id (const _pTIO_Granule_Ident_Type *gid)
{
   double sat_day;
   long scan_id;

   /* All radiance products will have scan_num >= 0. */

   if (gid->scan_num < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported product: scan_num = %d", __func__, gid->scan_num);
        return -1;
     }

   if (0 != tio_time_sat_local_day_number (gid->tstart, &sat_day))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computing satellite-local day number", __func__);
        return -1;
     }

   scan_id = 1000L * ((long) sat_day) + gid->scan_num;

   if (fprintf (stdout, "%ld", scan_id) < 0)
     return -1;

   return 0;
}

static int print_product_month (const _pTIO_Granule_Ident_Type *gid)
{
   struct tm tstart;
   const char *month[] = {
      "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
      "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
   };

   if ((-1 == _pTIO_parse_timestr (gid->tstart_str, &tstart))
       || (tstart.tm_mon < 0 || tstart.tm_mon > 11))
     return -1;

   if (fprintf (stdout, "%s", month[tstart.tm_mon]) < 0)
     return -1;

   return 0;
}

int main (int argc, char **argv)
{
   _pTIO_Granule_Ident_Type gid = {0};
   char product_type[MAX_PRODUCT_TYPE_LEN];
   char *level1_file = NULL;
   int task = TASK_UNKNOWN;
   int status = 1;
   int ncid, xtype;
   size_t len;
   static struct option long_options[] =
     {
        {"dir",      no_argument, 0, 'd'},
        {"localday", no_argument, 0, 'l'},
        {"month",    no_argument, 0, 'm'},
        {"scanid",   no_argument, 0, 's'},
	{"help",     no_argument, 0, 'h'},
        {0,0,0,0}
     };

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "dlms", long_options, &option_index);
        if (c == -1) break;

        switch (c)
          {
           case 'd':
             task = TASK_PRINT_DIR;
             break;

           case 'l':
             task = TASK_PRINT_SATDAY_DIR;
             break;

           case 'm':
             task = TASK_PRINT_MONTH;
             break;

           case 's':
             task = TASK_PRINT_SCANID;
             break;

           case '?':
             fprintf (stderr, "Unknown option -%c'.\n", optopt);
             usage (argc, argv);
             return 1;

           default:
             usage (argc, argv);
             return 1;
          }
     }

   if (argc - optind < 1)
     {
        usage (argc, argv);
        return 1;
     }

   if (task == TASK_UNKNOWN)
     {
        usage (argc, argv);
        return 1;
     }

   level1_file = argv[optind];

   if (0 != TIO_open (level1_file, NC_NOWRITE, &ncid))
     return 1;

   if ((0 != tio_use_file_epoch (ncid))
       || (-1 == _pTIO_read_granule_ident (ncid, &gid)))
     {
        (void) TIO_close (ncid);
        return 1;
     }

   memset (product_type, 0, MAX_PRODUCT_TYPE_LEN);

   if ((NC_NOERR == nc_inq_att (ncid, NC_GLOBAL, "product_type", &xtype, &len))
       && (len < MAX_PRODUCT_TYPE_LEN))
     {
        if (0 != TIO_get_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, product_type))
          {
             (void) TIO_close (ncid);
             return 1;
          }
     }

   (void) TIO_close (ncid);

   switch (task)
     {
      case TASK_PRINT_DIR:
        status = print_archive_subdir (&gid, product_type);
        break;

      case TASK_PRINT_MONTH:
        status = print_product_month (&gid);
        break;

      case TASK_PRINT_SATDAY_DIR:
        status = print_sat_local_day_number (&gid);
        break;

      case TASK_PRINT_SCANID:
        status = print_scan_id (&gid);
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid task", __func__);
        break;
     }

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
