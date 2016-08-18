#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define __USE_XOPEN
#include <time.h>

#include <limits.h>
#include <stddef.h>
#include <unistd.h>

#include <netcdf.h>
#include <tell.h>

#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

#define BUFSIZE 1024

static void usage (int argc, char **argv)
{
   (void) argc;
   fprintf (stderr, "Usage: %s [options] <radiance-file>\n", argv[0]);
   fprintf (stderr, "Options:\n");
   fprintf (stderr, " -o <info-file> Output file name\n");
   fprintf (stderr, " -v <version>   Processing version number\n");
}

static int write_granule_ident_file (const char *file,
                                     _pTIO_Granule_Ident_Type *gid,
                                     int version)
{
   struct tm tstart, tend;
   FILE *fp;

   if ((-1 == _pTIO_parse_timestr (gid->tstart_str, &tstart))
       || (-1 == _pTIO_parse_timestr (gid->tend_str, &tend)))
     return -1;

   if (NULL == (fp = fopen (file, "w")))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: opening file %s for writing\b",
                     __func__, file);
        return -1;
     }

   fprintf (fp, "processing_version,%d\n", version);
   fprintf (fp, "scan_seq_num,%d\n", gid->scan_seq_num);
   fprintf (fp, "granule_seq_num,%d\n", gid->granule_seq_num);
   fprintf (fp, "granule_num,%d\n", gid->granule_num);

   fprintf (fp, "tstart_year,%d\n", tstart.tm_year + 1900);
   fprintf (fp, "tstart_month,%d\n", tstart.tm_mon + 1);
   fprintf (fp, "tstart_mday,%d\n", tstart.tm_mday);
   fprintf (fp, "tstart_hour,%d\n", tstart.tm_hour);
   fprintf (fp, "tstart_min,%d\n", tstart.tm_min);
   fprintf (fp, "tstart_sec,%d\n", tstart.tm_sec);
   fprintf (fp, "tstart_wday,%d\n", tstart.tm_wday);
   fprintf (fp, "tstart_yday,%d\n", tstart.tm_yday);

   fprintf (fp, "tend_year,%d\n", tend.tm_year + 1900);
   fprintf (fp, "tend_month,%d\n", tend.tm_mon + 1);
   fprintf (fp, "tend_mday,%d\n", tend.tm_mday);
   fprintf (fp, "tend_hour,%d\n", tend.tm_hour);
   fprintf (fp, "tend_min,%d\n", tend.tm_min);
   fprintf (fp, "tend_sec,%d\n", tend.tm_sec);
   fprintf (fp, "tend_wday,%d\n", tend.tm_wday);
   fprintf (fp, "tend_yday,%d\n", tend.tm_yday);

   errno = 0;
   if (0 != fclose (fp))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s; errno=%d (%s)\n",
                     __func__, file, errno, strerror(errno));
        return -1;
     }

   return 0;
}

int main (int argc, char **argv)
{
   _pTIO_Granule_Ident_Type gid;
   char *output_file = NULL;
   char *radiance_file = NULL;
   int c, ncid;
   int version = INT_MAX; /* version=1 default might cause trouble here */

   while ((c = getopt (argc, argv, "o:v:")) != -1)
     {
        switch (c)
          {
           case 'o':
             output_file = optarg;
             break;

           case 'v':
             if (1 != sscanf (optarg, "%d", &version))
               {
                  fprintf (stderr, "*** ERROR: invalid version: %s\n", optarg);
                  return 1;
               }
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

   if (output_file == NULL)
     {
        usage (argc, argv);
        return 1;
     }

   if (argc - optind < 1)
     {
        usage (argc, argv);
        return 1;
     }

   radiance_file = argv[optind];

   if (0 != TIO_open (radiance_file, NC_NOWRITE, &ncid))
     return 1;

   if (-1 == _pTIO_read_granule_ident (ncid, &gid))
     {
        (void) TIO_close (ncid);
        return -1;
     }

   (void) TIO_close (ncid);

   if (0 != write_granule_ident_file (output_file, &gid, version))
     return -1;

   return 0;
}
