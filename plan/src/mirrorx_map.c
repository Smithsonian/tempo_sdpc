#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <wordexp.h>
#include <libconfig.h>

#include "tio.h"

#include "compute_scan_angles.h"

#define DEGTORAD       (M_PI/180.0)
#define SMA_MAX_CALIBRATED_MIRROR_X  49600.0

static int set_geometry_params (config_t *cfg, double *sat_lon, double *bs_lon, double *bs_lat)
{
   config_setting_t *s;
   config_setting_t *sub;
   double ewbias, nsbias, clockingbias, telescopeOffset;

   if (NULL == (s = config_lookup (cfg, "sat_config")))
     {
        fprintf (stderr, "*** Error: %s: accessing sat_config in param file: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "ewbias", &ewbias))
     {
        fprintf (stderr, "*** Error: %s: reading ewbias: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "nsbias", &nsbias))
     {
        fprintf (stderr, "*** Error: %s: reading nsbias: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "clockingbias", &clockingbias))
     {
        fprintf (stderr, "*** Error: %s: reading clockingbias: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "telescopeOffset", &telescopeOffset))
     {
        fprintf (stderr, "*** Error: %s: reading telescopeOffset: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "sat_lon", sat_lon))
     {
        fprintf (stderr, "*** Error: %s: reading sat_lon: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (sub = config_setting_get_member (s, "boresight")))
     return -1;

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "lon", bs_lon))
     {
        fprintf (stderr, "*** Error: %s: reading boresight longitude: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "lat", bs_lat))
     {
        fprintf (stderr, "*** Error: %s: reading boresight latitude: %s",
                 __func__, config_error_file (cfg));
        return -1;
     }

   return __set_geometry_params (ewbias, nsbias, clockingbias, telescopeOffset);
}

static char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        fprintf (stderr, "*** Error: %s: expanding path: %s\n", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        fprintf (stderr, "*** Error: %s: strdup failed\n", __func__);
     }

   return s_exp;
}

static int find_config_file (char **config_file)
{
   const char *defaults[] = {"plan.cfg", "$SDPC_ROOT/share/plan.cfg", NULL};
   const char **p;
   char *s;

   for (p = defaults; *p != NULL; p++)
     {
        if (NULL == (s = expand_string (*p)))
          return -1;
        if (0 == access (s, F_OK | R_OK))
          break;
        free (s);
     }

   *config_file = s;
   return s ? 1 : 0;
}

static void usage (void)
{
   fprintf (stderr, "Usage: mirrorx_map [options] <output-filename>\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -c | --config FILE       Configuration file path\n");
   fprintf (stderr, "   -o | --lon delta,min,num longitude grid [deg]\n");
   fprintf (stderr, "   -a | --lat delta,min,num latitude grid [deg]\n");
   fprintf (stderr, "   -S | --satlon LON        Host satellite longitude [deg]\n");
   exit (EXIT_SUCCESS);
}

static int write_mirrorx_table (const char *outfile, int num_lons, const float *lons,
                                int num_lats, const float *lats, const float *mirrorx)
{
   int dimids[2], start[2], count[2];
   int ncid, dimid_lon, dimid_lat;
   int varid_lon, varid_lat, varid_mirrorx;
   int status = -1;

   if (0 != TIO_create (outfile, NC_NETCDF4, &ncid))
     return -1;

   if ((0 != TIO_def_dim (ncid, "longitude", num_lons, &dimid_lon))
       || (0 != TIO_def_dim (ncid, "latitude", num_lats, &dimid_lat)))
     goto return_status;

   dimids[0] = dimid_lat;
   dimids[1] = dimid_lon;

   if ((0 != TIO_def_var (ncid, "longitude", NC_FLOAT, 1, &dimid_lon, &varid_lon))
       || (0 != TIO_def_var (ncid, "latitude", NC_FLOAT, 1, &dimid_lat, &varid_lat))
       || (0 != TIO_def_var (ncid, "mirrorx", NC_FLOAT, 2, dimids, &varid_mirrorx)))
     {
        goto return_status;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_lats;
   count[1] = num_lons;

   if ((0 != TIO_put_var_section (ncid, "latitude", &start[0], &count[0], NC_FLOAT, lats))
       || (0 != TIO_put_var_section (ncid, "longitude", &start[1], &count[1], NC_FLOAT, lons))
       || (0 != TIO_put_var_section (ncid, "mirrorx", start, count, NC_FLOAT, mirrorx)))
     goto return_status;

   status = 0;
return_status:
   (void) TIO_close (ncid);
   return status;
}

typedef struct
{
   double min, delta;
   int num;
}
Grid_Type;

int main (int argc, char **argv)
{
   config_t cfg = {0};
   static struct option long_options[] =
     {
        {"help",   no_argument,       0, 'h'},
        {"config", required_argument, 0, 'c'},
        {"satlon", required_argument, 0, 'S'},
        {"lon",    required_argument, 0, 'o'},
        {"lat",    required_argument, 0, 'a'},
        {0,0,0,0}
     };
   char *config_file = NULL;
   char *outfile = NULL;
   double nan_value = nan("");
   double sat_lon = nan_value;
   Grid_Type lon = {.min = -168.0, .delta = 0.1,/*0.02*/ .num = 1550 /*7750*/};
   Grid_Type lat = {.min =   14.0, .delta = 0.1,/*0.02*/ .num =  590 /*2950*/};
   float *lon_array = NULL;
   float *lat_array = NULL;
   float *mirrorx = NULL;
   float mirrorx_fill_value = -9999.0;
   int status = -1;
   int malloced_config_file = 0;
   EarthPoint earth_point = {0};
   AziElev_Type scan_angles = {0};
   AziElev_Type mirror_angles = {0};
   double cfg_sat_lon, bs_lon, bs_lat;
   int i, j, num_fail;

   config_init (&cfg);

   if ((malloced_config_file = find_config_file (&config_file)) < 0)
     goto exit_status;

   /* If we found a config file, read it now, otherwise, keep going
    * in case there's a config file on the command line */
   if (config_file)
     {
        if (0 == config_read_file (&cfg, config_file))
          {
             fprintf (stderr, "*** Error: reading %s:%d - %s",
                      config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
             goto exit_status;
          }
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hc:s:o:a:g:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "*** Error: getopt returned character %d??", c);
             usage();
             break;
           case 0:
             /* handle long-only options */
             break;
           case 'c':
             if (0 != access (optarg, F_OK | R_OK))
               {
                  fprintf (stderr, "*** Error: cannot access %s\n", optarg);
                  goto exit_status;
               }
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (malloced_config_file)
               {
                  free(config_file);
                  malloced_config_file = 0;
               }
             config_file = optarg;
             if (0 == config_read_file (&cfg, config_file))
               {
                  fprintf (stderr, "*** Error: reading %s:%d - %s",
                           config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
                  goto exit_status;
               }
             break;
           case 'h':
             usage();
             break;
           case 's':
             if (1 != sscanf (optarg, "%le", &sat_lon))
               {
                  usage();
               }
             break;
           case 'a':
             if (1 != sscanf (optarg, "%le,%le,%d", &lat.delta, &lat.min, &lat.num))
               {
                  usage();
               }
             break;
           case 'o':
             if (1 != sscanf (optarg, "%le,%le,%d", &lon.delta, &lon.min, &lon.num))
               {
                  usage();
               }
             break;
          }
     }

   if (config_file == NULL) usage();

   if (optind == argc)
     usage();

   outfile = argv[optind++];

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   if (0 != set_geometry_params (&cfg, &cfg_sat_lon, &bs_lon, &bs_lat))
     goto exit_status;

   if (isnan(sat_lon)) sat_lon = cfg_sat_lon;
   sat_lon *= DEGTORAD;

   if ((NULL == (lon_array = (float *)malloc (lon.num * sizeof(float))))
       || (NULL == (lat_array = (float *)malloc (lat.num * sizeof(float))))
       || (NULL == (mirrorx = (float *)malloc (lat.num * lon.num * sizeof(float)))))
     {
        fprintf (stderr, "*** malloc failed\n");
        goto exit_status;
     }

   for (j = 0; j < lat.num; j++)
     {
        lat_array[j] = lat.min + j * lat.delta;
     }
   for (i = 0; i < lon.num; i++)
     {
        lon_array[i] = lon.min + i * lon.delta;
     }

   earth_point.theAlt = 0.0;
   num_fail = 0;
   for (i = 0; i < lat.num; i++)
     {
        int offset = i * lon.num;
        earth_point.theLat = lat_array[i];
        fprintf (stderr, "%04d/%4d\r", i, lat.num);
        for (j = 0; j < lon.num; j++)
          {
             double mx;

             earth_point.theLon = lon_array[j];
             if (0 == __compute_scan_angles (&earth_point, sat_lon, &scan_angles))
               {
                  mirror_angles.azimuth = 0.5 * scan_angles.azimuth;
                  mirror_angles.elevation = 0.5 * scan_angles.elevation;
                  mx = mirror_angles.azimuth;
                  if (fabs(mx) > SMA_MAX_CALIBRATED_MIRROR_X) mx = mirrorx_fill_value;
               }
             else
               {
                  mx = mirrorx_fill_value;
                  num_fail++;
               }

             mirrorx[offset + j] = mx;
          }
     }

   if (0 != write_mirrorx_table (outfile, lon.num, lon_array, lat.num, lat_array, mirrorx))
     goto exit_status;

   fprintf (stderr, "done (num_fail=%d)\n", num_fail);

   status = 0;
exit_status:
   if (malloced_config_file) free(config_file);
   free(lon_array);
   free(lat_array);
   free(mirrorx);
   config_destroy (&cfg);
   return status;
}
