#include <stdio.h>
#include <unistd.h>
#include <tell.h>

#include <netcdf.h>

#include "tio.h"
#include "tio_template.h"

#define BUFSIZE 1024

static void usage (int argc, char **argv)
{
   (void) argc;
   fprintf (stderr, "Usage: %s [options] -p <prod> <radiance-file>\n", argv[0]);
   fprintf (stderr, "Options:\n");
   fprintf (stderr, " -p <prod>      L2 product abbreviation, e.g. no2\n");
   fprintf (stderr, " -v <version>   Processing version number\n");
}

int main (int argc, char **argv)
{
   char buf[BUFSIZE];
   char *prod_abbrev = NULL;
   char *radiance_file = NULL;
   int c, n, ncid, version = 1;

   while ((c = getopt (argc, argv, "p:v:")) != -1)
     {
        switch (c)
          {
           case 'p':
             prod_abbrev = optarg;
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

   if (argc - optind < 1)
     {
        usage (argc, argv);
        return 1;
     }

   radiance_file = argv[optind];

   if (0 != TIO_open (radiance_file, NC_NOWRITE, &ncid))
     return 1;

   n = TIO_filename_from_granule (ncid, prod_abbrev, version,
                                  buf, BUFSIZE);
   (void) TIO_close (ncid);

   if (n >= BUFSIZE)
     {
        Tell_verror (TELL_RUNTIME_ERROR,
                     "filename truncated; length %d exceeds buffer size %d",
                     n, BUFSIZE);
        return 1;
     }

   fprintf (stdout, "%s\n", buf);
   return 0;
}
