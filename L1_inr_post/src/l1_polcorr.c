/** @file l1_polcorr.c
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
#include <proj_api.h>
#include <gsl/gsl_errno.h>
#include <tell.h>
#include <tio_template.h>

#include "config.h"

#include "lps.h"

/* Update the fortran interface whenever this definition changes */
typedef struct
{
   const char *rad_file;
   const char *qu_file;
   Lps_Type *lps;
   int uv_beg, uv_end;
   int vis_beg, vis_end;
   int merge_bands;
   int use_mler;
   int debug_output;
   int step;
   int xtrack;
   double delta_pa;
}
Polcorr_Type;

static void free_polcorr_type (Polcorr_Type *pt)
{
   if (pt == NULL)
     return;
}

extern int polcorrect (const Polcorr_Type *pt);

static void usage (void)
{
   fprintf (stderr, "Usage: L1_polcorr [options] <input-file>\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -c | --config FILE     configuration file\n");
   fprintf (stderr, "   -s | --step s          step index to process (1 <= s <= num_steps_in_granule)\n");
   fprintf (stderr, "   -x | --xtrack x        xtrack index to process (1 <= x <= num_xtrack)\n");
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

static int
read_retrieval_limits (config_setting_t *s, const char *name,
                      int *beg, int *end)
{
   config_setting_t *sub;
   int num;

   if (NULL == (sub = config_setting_get_member (s, name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing config setting %s in %s",
                     __func__, name, config_setting_source_file (s));
        return -1;
     }

   if ((num = config_setting_length (sub)) != 2)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: invalid config setting %s in %s",
                     __func__, name, config_setting_source_file (s));
        return -1;
     }

   *beg = config_setting_get_int_elem (sub, 0);
   *end = config_setting_get_int_elem (sub, 1);

   return 0;
}

static int process_inputs (config_t *cfg, const char *rad_file,
                           int step, int xtrack)
{
   Polcorr_Type pt = {0};
   config_setting_t *s;
   int status = -1;

   pt.rad_file = rad_file;
   pt.step = step;
   pt.xtrack = xtrack;

   if (NULL == (s = config_lookup (cfg, "polcorr_tables")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing polcorr in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "qu_lut", &pt.qu_file))
       || (CONFIG_TRUE != config_setting_lookup_bool (s, "use_mler", &pt.use_mler))
       || (CONFIG_TRUE != config_setting_lookup_bool (s, "merge_bands", &pt.merge_bands))
       || (CONFIG_TRUE != config_setting_lookup_bool (s, "debug_output", &pt.debug_output))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "delta_pa", &pt.delta_pa)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing %s in param file: %s",
                     __func__, config_setting_name(s),
                     config_setting_source_file (s));
        return -1;
     }

   if ((0 != read_retrieval_limits (s, "waves_uv", &pt.uv_beg, &pt.uv_end))
       || (0 != read_retrieval_limits (s, "waves_vis", &pt.vis_beg, &pt.vis_end)))
     {
        goto return_status;
     }

   if (NULL == (pt.lps = lps_open (cfg)))
     goto return_status;

   if (0 != polcorrect (&pt))
     return -1;

   status = 0;
return_status:
   lps_close (pt.lps);
   free_polcorr_type (&pt);

   return status;
}

int main (int argc, char **argv)
{
   const char appname[] = "L1_polcorr";
   char *config_file = "l1_inr_post.cfg";
   config_t cfg;
   char *input_file = NULL;
   int status = EXIT_FAILURE;
   int step = 0;
   int xtrack = 0;
   static struct option long_options[] =
     {
        {"config",  required_argument, 0, 'c'},
        {"step",    required_argument, 0, 's'},
        {"xtrack",  required_argument, 0, 'x'},
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
        int c = getopt_long (argc, argv, "c:s:x:v:", long_options, &option_index);
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
           case 'c': config_file = optarg;
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (-1 == read_config_file (config_file, &cfg))
               goto return_status;
             break;
           case 's':
             if (1 != sscanf (optarg, "%d", &step))
               goto return_status;
             break;
           case 'x':
             if (1 != sscanf (optarg, "%d", &xtrack))
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
     {
        config_destroy(&cfg);
        tell_close();
        usage();
     }

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

   gsl_set_error_handler_off();
   tell_vlog (TELL_MSGTYPE_INFO, 0, "start %s", input_file);

   status = process_inputs (&cfg, input_file, step, xtrack);

return_status:
   config_destroy (&cfg);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "status=%d, finished %s",
              status, input_file);
   tell_close();

   return status;
}
