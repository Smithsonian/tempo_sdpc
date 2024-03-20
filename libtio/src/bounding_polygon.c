#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>

#include <tell.h>

#include "tio.h"

static int process_file (const char *ncfile, int replace)
{
   TIO_Meta_Type *meta = NULL;
   float *lon=NULL, *lat=NULL;
   int ncid, grp, num, mode, i;
   int status = -1;

   mode = replace ? NC_WRITE : NC_NOWRITE;

   if ((0 != TIO_open (ncfile, mode, &ncid))
       || (0 != TIO_inq_grp (ncid, "band_290_490_nm", &grp)))
     {
        fprintf (stderr, "*** Error opening netcdf file metadata group: %s\n", ncfile);
        goto cleanup_and_exit;
     }

   if (replace)
     {
        if (NULL == (meta = tio_meta_open ()))
          goto cleanup_and_exit;
        if (0 != tio_meta_set_lev1_bounding_polygon (meta, grp))
          goto cleanup_and_exit;
        if (0 != tio_meta_write_ncattr (meta, ncid))
          goto cleanup_and_exit;
     }
   else
     {
        if (0 != __tio_make_lev1_bounding_polygon (grp, &num, &lon, &lat))
          goto cleanup_and_exit;
        fprintf (stderr, "bounding_polygon has %d vertices\n", num);

        for (i = 0; i < num; i++)
          {
             fprintf (stdout, "%5d %10.7e %10.7e\n", i, lon[i], lat[i]);
          }
     }

   status = 0;

cleanup_and_exit:
   tio_meta_close (meta);
   if (lon) free(lon);
   if (lat) free(lat);
   (void) TIO_close (ncid);

   return status;
}

static int usage (const char *pgm)
{
   fprintf (stdout, "Usage: %s [options] <level1b-radiance-file>\n", pgm);
   fprintf (stdout, "Options:\n");
   fprintf (stdout, " -h | --help              Print this listing\n");
   fprintf (stdout, " -r | --replace           Replace the bounding polygon in the specified file\n");
   fprintf (stdout, " -b | --band BAND_KM      Preserve boundary details larger than this length scale [km]\n");
   fprintf (stdout, " -l | --limb MAX_VZA_DEG  Maximum VZA to include [DEG]\n");
   exit(EXIT_SUCCESS);
}

int main (int argc, char **argv)
{
   const char *pgm = "bounding_polygon";
   static struct option long_options[] =
     {
        {"help",    no_argument,       0, 'h'},
        {"replace", no_argument,       0, 'r'},
        {"band",    required_argument, 0, 'b'},
        {"limb",    required_argument, 0, 'l'},
        {0,0,0,0}
     };
   int status, replace = 0;
   float band_km = -1.0;      /* band_km < 0 means use default */
   float vza_max_deg = -1.0;  /* vza_max_deg < 0 means use default */
   char *path = NULL;

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hrb:l:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d??\n", c);
             break;
           case 'h':
             usage (pgm);
             break;
           case 'r':
             replace++;
             break;
           case 'b':
             if (1 != sscanf (optarg, "%f", &band_km))
               {
                  fprintf (stderr, "*** Error: reading argument: %s\n", optarg);
                  return 1;
               }
             break;
           case 'l':
             if (1 != sscanf (optarg, "%f", &vza_max_deg))
               {
                  fprintf (stderr, "*** Error: reading argument: %s\n", optarg);
                  return 1;
               }
             break;
          }
     }

   if (argc - optind < 1)
     {
        usage (pgm);
        return 1;
     }

   path = argv[optind];

   tell_open (pgm, -1, 0);
   (void) __tio_set_bounding_polygon_controls (band_km, vza_max_deg);
   status = process_file (path, replace);
   tell_close();

   if (status)
     {
        fprintf (stderr, "*** %s failed\n", pgm);
     }

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
