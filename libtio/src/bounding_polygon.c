#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
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
       || (0 != TIO_def_grp (ncid, "band_290_490_nm", &grp)))
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
   fprintf (stdout, "Usage: %s [--replace] FILE\n", pgm);
   exit(EXIT_SUCCESS);
}

int main (int argc, char **argv)
{
   const char *pgm = "bounding_polygon";
   int status, replace = 0;
   char *path = NULL;

   if ((argc < 2)
       || (0 == strcmp ("-h", argv[1]))
       || (0 == strcmp ("--help", argv[1])))
     {
        usage (pgm);
     }

   switch (argc)
     {
      default:
        /* FALLTHRU */
      case 2:
        path = argv[1];
        break;

      case 3:
        if (0 == strcmp ("--replace", argv[1]))
          {
             replace++;
             path = argv[2];
          }
        else usage (pgm);
        break;
     }

   tell_open (pgm, -1, 0);
   status = process_file (path, replace);
   tell_close();

   if (status)
     {
        fprintf (stderr, "*** %s failed\n", pgm);
     }

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
