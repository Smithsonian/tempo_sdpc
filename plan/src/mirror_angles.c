#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <libconfig.h>

#include "compute_scan_angles.h"

#define DEGTORAD       (M_PI/180.0)

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

static void usage (void)
{
   fprintf (stderr, "Usage: mirror_angles [options]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -c | --config FILE       Configuration file path\n");
   fprintf (stderr, "   -o | --lon LON           Surface point longitude [deg]\n");
   fprintf (stderr, "   -a | --lat LAT           Surface point latitude [deg]\n");
   fprintf (stderr, "   -g | --gridlon n[,dlon]  Number of longitude steps; longitude step size [deg]\n");
   fprintf (stderr, "   -S | --satlon LON        Host satellite longitude [deg]\n");
   exit (EXIT_SUCCESS);
}

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
        {"gridlon", required_argument, 0, 'g'},
        {0,0,0,0}
     };
   char *config_file = "plan.cfg";
   double nan_value = nan("");
   double sat_lon = nan_value;
   double lon = nan_value;
   double lat = nan_value;
   double dlon = 0.25;
   int status = -1;
   int num_lons = 0;
   EarthPoint earth_point = {0};
   AziElev_Type scan_angles = {0};
   AziElev_Type mirror_angles = {0};
   double cfg_sat_lon, bs_lon, bs_lat;

   config_init (&cfg);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hs:o:a:g:", long_options, &option_index);
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
             config_file = optarg;
             break;
           case 'h':
             usage();
             break;
           case 'g':
             if (NULL == strchr (optarg, ','))
               {
                  if (1 != sscanf (optarg, "%d", &num_lons))
                    {
                       usage();
                    }
               }
             else if (2 != sscanf (optarg, "%d,%le", &num_lons, &dlon))
               {
                  usage();
               }
             break;
           case 's':
             if (1 != sscanf (optarg, "%le", &sat_lon))
               {
                  usage();
               }
             break;
           case 'a':
             if (1 != sscanf (optarg, "%le", &lat))
               {
                  usage();
               }
             break;
           case 'o':
             if (1 != sscanf (optarg, "%le", &lon))
               {
                  usage();
               }
             break;
          }
     }

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   if (0 == config_read_file (&cfg, config_file))
     {
        fprintf (stderr, "*** Error: reading %s:%d - %s",
                 config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        goto exit_status;
     }

   if (0 != set_geometry_params (&cfg, &cfg_sat_lon, &bs_lon, &bs_lat))
     goto exit_status;

   if (isnan(sat_lon)) sat_lon = cfg_sat_lon;
   if (isnan(lon)) lon = bs_lon;
   if (isnan(lat)) lat = bs_lat;

   sat_lon *= DEGTORAD;

   if ((num_lons < 1) || (dlon == 0.0))
     {
        earth_point.theLon = lon;
        earth_point.theLat = lat;
        earth_point.theAlt = 0.0;

        if (0 != __compute_scan_angles (&earth_point, sat_lon, &scan_angles))
          {
             fprintf (stderr, "*** Error: compute_scan_angles failed\n");
             goto exit_status;
          }

        mirror_angles.azimuth = 0.5 * scan_angles.azimuth;
        mirror_angles.elevation = 0.5 * scan_angles.elevation;

        fprintf (stdout, "     (lon,     lat) = (%9.4f, %9.4f) deg\n", lon, lat);
        fprintf (stdout, "(mirror_x,mirror_y) = (%9.2f, %9.2f) microradian\n", mirror_angles.azimuth, mirror_angles.elevation);
        fprintf (stdout, "  (scan_x,  scan_y) = (%9.2f, %9.2f) microradian\n", scan_angles.azimuth, scan_angles.elevation);
     }
   else
     {
        double scan_azimuth0 = 0.0;
        int i;

        earth_point.theLat = lat;
        earth_point.theAlt = 0.0;

        fprintf (stdout, "longitude,mirror_x,scan_x,delta_scan_x\n");
        for (i = 0; i < num_lons; i++)
          {
             earth_point.theLon = lon + i * dlon;
             if (0 != __compute_scan_angles (&earth_point, sat_lon, &scan_angles))
               {
                  fprintf (stderr, "*** Error: compute_scan_angles failed\n");
                  goto exit_status;
               }

             if (i == 0) scan_azimuth0 = scan_angles.azimuth;

             mirror_angles.azimuth = 0.5 * scan_angles.azimuth;
             mirror_angles.elevation = 0.5 * scan_angles.elevation;

             fprintf (stdout, "%9.4f,%9.2f,%9.2f,%9.2f\n",
                      earth_point.theLon, mirror_angles.azimuth,
                      scan_angles.azimuth, scan_angles.azimuth - scan_azimuth0);
          }
     }

   status = 0;
exit_status:
   config_destroy (&cfg);
   return status;
}
