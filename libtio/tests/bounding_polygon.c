#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <tio.h>
#include <tell.h>

#include "simplify.h"

static int process_file (const char *ncfile)
{
   int ncid, grp, num, i;
   float *lon=NULL, *lat=NULL;
   int *indices = NULL;
   float band_km = 5.0;
   int num_kept, status = -1;

   if ((0 != TIO_open (ncfile, NC_NOWRITE, &ncid))
       || (0 != TIO_def_grp (ncid, "band_290_490_nm", &grp)))
     {
        fprintf (stderr, "*** Error opening netcdf file metadata group: %s\n", ncfile);
        goto cleanup_and_exit;
     }

   if (0 != __tio_make_lev1_bounding_polygon (grp, &num, &lon, &lat, NULL, NULL))
     goto cleanup_and_exit;

   if ((num_kept = simplify_dp (lon, lat, num, band_km, &indices)) < 0)
     goto cleanup_and_exit;

   for (i = 0; i < num_kept; i++)
     {
        int k = indices[i];
        fprintf (stdout, "%5d %10.7e %10.7e\n", i, lon[k], lat[k]);
     }

   fprintf (stderr, "bounding_polygon has %d vertices\n", num_kept);

   status = 0;

cleanup_and_exit:
   if (lon) free(lon);
   if (lat) free(lat);
   if (indices) free(indices);
   (void) TIO_close (ncid);

   return status;
}

int main (int argc, char **argv)
{
   int status;

   if (argc != 2)
     {
        fprintf (stdout, "Usage: %s FILE\n", argv[0]);
        return 0;
     }

   tell_open ("bounding_polygon", -1, 0);
   status = process_file(argv[1]);
   tell_close();

   if (status)
     {
        fprintf (stderr, "*** %s failed\n", argv[0]);
     }

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
