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
#include <wordexp.h>

#include <libconfig.h>
#include <proj_api.h>
#include <gsl/gsl_errno.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"

#include "lps.h"
#include "lps_apply.h"

static char *_pCommand_Line;

enum
{
   TASK_LPS_CORRECT = 0,
   TASK_LPS_APPLY = 1
};

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
   int diag_output;
   int step;
   int xtrack;
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
   fprintf (stderr, "   -h | --help            print this usage message\n");
   fprintf (stderr, "   -c | --config FILE     configuration file\n");
   fprintf (stderr, "   -d | --diag            generate diagnostic output\n");
   fprintf (stderr, "   -a | --apply           apply LPS error for Q,U in input file\n");
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

static int process_inputs (config_t *cfg, const char *rad_file, int task,
                           int step, int xtrack, int diag_output)
{
   Polcorr_Type pt = {0};
   wordexp_t we;
   const char *qu_file;
   config_setting_t *s;
   int status = -1;

   memset ((char *)&we, 0, sizeof(wordexp_t));

   pt.rad_file = rad_file;
   pt.step = step;
   pt.xtrack = xtrack;
   pt.diag_output = diag_output;

   if (NULL == (s = config_lookup (cfg, "polcorr_tables")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing polcorr in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "qu_lut", &qu_file))
       || (CONFIG_TRUE != config_setting_lookup_bool (s, "use_mler", &pt.use_mler))
       || (CONFIG_TRUE != config_setting_lookup_bool (s, "merge_bands", &pt.merge_bands)))
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

   if ((0 != wordexp (qu_file, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, qu_file);
        goto return_status;
     }
   pt.qu_file = we.we_wordv[0];

   if (NULL == (pt.lps = lps_open (cfg)))
     goto return_status;

   switch (task)
     {
      case TASK_LPS_CORRECT:
        status = polcorrect (&pt);
        break;

      case TASK_LPS_APPLY:
        status = lps_apply (pt.lps, rad_file, _pCommand_Line);
        break;

      default:
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: unsupported task id = %d", __func__, task);
        break;
     }

return_status:
   wordfree (&we);
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
   int task = TASK_LPS_CORRECT;
   int diag_output = 0;
   int step = 0;
   int xtrack = 0;
   static struct option long_options[] =
     {
        {"config",  required_argument, 0, 'c'},
        {"diag",    no_argument,       0, 'd'},
        {"apply",   no_argument,       0, 'a'},
        {"help",    no_argument,       0, 'h'},
        {"step",    required_argument, 0, 's'},
        {"xtrack",  required_argument, 0, 'x'},
        {"verbose", optional_argument, 0, 'v'},
        {0,0,0,0}
     };

   if (argc < 2)
     usage();

   /* NULL return is ok */
   _pCommand_Line = tio_concat_argv (argc, argv, NULL, 0);

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
        int c = getopt_long (argc, argv, "ahc:ds:x:v:", long_options, &option_index);
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
           case 'h':
             config_destroy(&cfg);
             tell_close();
             usage();
             break;
           case 'c': config_file = optarg;
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (-1 == read_config_file (config_file, &cfg))
               goto return_status;
             break;
           case 'd':
             diag_output = 1;
             break;
           case 'a':
             task = TASK_LPS_APPLY;
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

   status = process_inputs (&cfg, input_file, task, step, xtrack, diag_output);

return_status:
   config_destroy (&cfg);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "status=%d, finished %s",
              status, input_file);
   tell_close();

   return status;
}
