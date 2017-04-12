/** @file l1_inr_post.c
 *  @brief Main program
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>

#include "config.h"
#include "granule.h"
#include "snow.h"
#include "elevation.h"
#include "land_cover.h"

static void usage (void)
{
   fprintf (stderr, "Usage: L1_inr_post [options] <input-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -s | --snow FILE       snow and ice mask file\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -c | --config FILE     configuration file\n");
   fprintf (stderr, "   -v | --verbose lev     logging verbosity\n");
   exit (EXIT_SUCCESS);
}

static int read_config_file (const char *config_file, config_t *cfg)
{
   if (0 == config_read_file (cfg, config_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: Reading %s:%d - %s",
                     __func__, config_error_file(cfg),
                     config_error_line(cfg), config_error_text(cfg));
        return -1;
     }

   return 0;
}

static int set_elevation (Granule_Type *gt, config_t *cfg)
{
   Elevation_Type *et;
   int status;

   if (NULL == (et = elevation_init (cfg)))
     return -1;

   status = gt->gt_set_elevation (gt, et);
   et->et_delete (et);

   return status;
}

static int read_angles (config_t *cfg, double *max_glint_angle_deg,
                       double *max_eclipse_angle_deg)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "ground_pixel_quality_flags")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ground_pixel_quality_flags in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "max_glint_angle_deg", max_glint_angle_deg))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading max_glint_angle_deg: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "max_eclipse_angle_deg", max_eclipse_angle_deg))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading max_eclipse_angle_deg: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int set_ground_pixel_quality_flags (Granule_Type *gt, config_t *cfg,
                                           const char *snow_file)
{
   Snow_Type *snow = NULL;
   Land_Cover_Type *land_cover = NULL;
   double max_glint_angle, max_eclipse_angle;
   int status = -1;

   if (0 != read_angles (cfg, &max_glint_angle, &max_eclipse_angle))
     return -1;

   if (NULL == (snow = snow_init (snow_file)))
     return -1;

   if (NULL == (land_cover = land_cover_init (cfg)))
     goto return_status;

   status = gt->gt_set_ground_pixel_flags (gt, max_glint_angle, max_eclipse_angle,
                                           snow, land_cover);

   status = 0;
return_status:
   if (snow) snow->sn_delete (snow);
   if (land_cover) land_cover->lc_delete (land_cover);

   return status;
}

static int process_inputs (Granule_Type *gt, config_t *cfg,
                           const char *snow_file)
{
   if (0 != set_elevation (gt, cfg))
     return -1;

   if (0 != gt->gt_set_object_angles (gt))
     return -1;

   if (0 != set_ground_pixel_quality_flags (gt, cfg, snow_file))
     return -1;

   return 0;
}

int main (int argc, char **argv)
{
   const char appname[] = "L1_inr_post";
   char *config_file = "l1_inr_post.cfg";
   config_t cfg;
   Granule_Type *gt = NULL;
   char *input_file = NULL;
   char *snow_file = NULL;
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"snow",    required_argument, 0, 's'},
        {"config",  optional_argument, 0, 'c'},
        {"verbose", optional_argument, 0, 'v'},
        {0,0,0,0}
     };

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);

   config_init (&cfg);

   /* Try reading the default config file, but if it doesn't exist,
    * keep going in case there's a config file on the command line */
   if (0 == access (config_file, F_OK | R_OK))
     {
        if (-1 == read_config_file (config_file, &cfg))
          goto return_status;
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "c:s:v:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: getopt returned character %d??",
                          __func__, c);
             goto return_status;
             break;
           case 's': snow_file = optarg;
             break;
           case 'c': config_file = optarg;
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (-1 == read_config_file (config_file, &cfg))
               goto return_status;
             break;
           case 'v':
               {
                  int log_level;
                  if (1 == sscanf (optarg, "%d", &log_level))
                    (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);
               }
             break;
          }
     }

   if (optind == argc)
     usage();

   input_file = argv[optind++];

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   tell_vlog (TELL_MSGTYPE_INFO, 0, "start %s", input_file);

   gt = granule_open (input_file);

   if (gt)
     {
        status = process_inputs (gt, &cfg, snow_file);
        gt->gt_close (gt);
     }

return_status:
   config_destroy (&cfg);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "status=%d, finished %s",
              status, input_file);

   return status;
}
