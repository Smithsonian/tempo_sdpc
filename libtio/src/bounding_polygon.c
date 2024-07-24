#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>

#include <tell.h>

#include "tio.h"

static int read_polygon (const char *srcpath, float **lon, float **lat, int *num)
{
   FILE *fp;
#define FILE_FMT "%d,%f,%f"
   char buf[128];
   float slon, slat;
   float *plon = NULL;
   float *plat = NULL;
   int n, i;

   if (NULL == (fp = fopen (srcpath, "r")))
     {
        fprintf (stderr, "*** Error opening file: %s\n", srcpath);
        return -1;
     }
   n = 0;
   for (;;)
     {
        memset ((char *)buf, 0, sizeof(buf));
        if (NULL == fgets (buf, sizeof(buf), fp))
          break;
        if (3 != sscanf (buf, FILE_FMT, &i,&slon,&slat))
          break;
        n++;
     }
   fclose (fp);

   if ((NULL == (plon = (float *)malloc (n * sizeof(float))))
       || (NULL == (plat = (float *)malloc (n * sizeof(float)))))
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        goto return_error;
     }

   if (NULL == (fp = fopen (srcpath, "r")))
     {
        fprintf (stderr, "*** Error opening file: %s\n", srcpath);
        goto return_error;
     }
   for (i = 0; i < n; i++)
     {
        int k;
        memset ((char *)buf, 0, sizeof(buf));
        if (NULL == fgets (buf, sizeof(buf), fp))
          break;
        if (3 != sscanf (buf, FILE_FMT, &k,&plon[i],&plat[i]))
          break;
     }
   fclose (fp);

   *lon = plon;
   *lat = plat;
   *num = n;
   return 0;

return_error:
   free (plon);
   free (plat);
   return -1;
}

static int process_source_polygon (TIO_Meta_Type *meta, const char *srcpath, float band_km)
{
   float *lon=NULL, *lat=NULL;
   float *smp_lon=NULL, *smp_lat=NULL;
   int *indices = NULL;
   int i, num, num_smp, num_kept, status = -1;

   if (0 != read_polygon (srcpath, &lon, &lat, &num))
     return -1;

   /* simplify the polygon */
   if ((num_kept = tio_meta_simplify_dp (lon, lat, num, band_km, &indices)) < 0)
     return -1;

   /* allocate room for an extra point if we need to close the polygon */
   num_smp = num_kept + 1;
   if (NULL == (smp_lon = (float *)malloc (2 * num_smp * sizeof(float))))
     {
        fprintf (stderr, "*** malloc failed\n");
        goto return_status;
     }
   smp_lat = smp_lon + num_smp;

   for (i = 0; i < num_kept; i++)
     {
        int k = indices[i];
        smp_lon[i] = lon[k];
        smp_lat[i] = lat[k];
     }

   /* If necessary, close the polygon */
   if ((smp_lon[num_kept-1] != smp_lon[0]) || (smp_lat[num_kept-1] != smp_lat[0]))
     {
        smp_lon[num_kept] = smp_lon[0];
        smp_lat[num_kept] = smp_lat[0];
        num_kept++;
     }

   /* set metadata entries */
   fprintf (stderr, "num_read=%d  num_kept=%d\n", num, num_kept);
   if (0 != tio_meta_set_odl_bounding_polygon (meta, smp_lon, smp_lat, num_kept))
     goto return_status;
   if (0 != tio_meta_set_acdd_geospatial_bounds (meta, smp_lon, smp_lat, num_kept))
     goto return_status;

   status = 0;
return_status:
   if (status)
     {
        free (lon);
        free (lat);
        free (indices);
        free (smp_lon);
        free (smp_lat);
     }

   return status;
}

static int process_file (const char *ncfile, int replace, const char *srcpath, float band_km)
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
        if (srcpath == NULL)
          {
             if (0 != tio_meta_set_lev1_bounding_polygon (meta, grp))
               goto cleanup_and_exit;
          }
        else
          {
             if (0 != process_source_polygon (meta, srcpath, band_km))
               goto cleanup_and_exit;
          }
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
             fprintf (stdout, "%5d,%10.7e,%10.7e\n", i, lon[i], lat[i]);
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
   fprintf (stdout, " -s | --source FILE       Read the bounding polygon from the specified file\n");
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
        {"source",  required_argument, 0, 's'},
        {"band",    required_argument, 0, 'b'},
        {"limb",    required_argument, 0, 'l'},
        {0,0,0,0}
     };
   int status, replace = 0;
   float band_km = -1.0;      /* band_km < 0 means use default */
   float vza_max_deg = -1.0;  /* vza_max_deg < 0 means use default */
   char *path = NULL;
   char *srcpath = NULL;

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hrs:b:l:", long_options, &option_index);
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
           case 's':
             srcpath = optarg;
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
   status = process_file (path, replace, srcpath, band_km);
   tell_close();

   if (status)
     {
        fprintf (stderr, "*** %s failed\n", pgm);
     }

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
