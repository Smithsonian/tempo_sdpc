#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>

#include "config.h"
#include "control.h"
#include "util.h"
#include "process.h"
#include "version.h"

static void usage (void)
{
   fprintf (stderr, "Usage: L0_ccd [options] <input-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -d | --dark FILE       Input corrected dark current file\n");
   fprintf (stderr, "   -o | --output FILE     Output file\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -s | --solar FILE      Solar irradiance file\n");
   fprintf (stderr, "   -b | --bpix FILE       Bad pixel file\n");
   fprintf (stderr, "   -c | --config FILE     Configuration file\n");
   fprintf (stderr, "   -i | --instr FILE      Instrument telemetry points file.\n");
   fprintf (stderr, "                          Provide a list using FILE=@DIR/paths.lis.\n");
   fprintf (stderr, "                          When FILE provides a directory path,\n");
   fprintf (stderr, "                          files matching hk_glob_pattern are examined\n");
   fprintf (stderr, "                          (hk_glob_pattern is defined in the config file)\n");
   fprintf (stderr, "   -n | --num N           Process <= N exposure records \n");
   fprintf (stderr, "   -t | --trend FILE      Trending parameter output file\n");
   fprintf (stderr, "   -v | --verbose         Verbosity (more instances means more verbose, e.g. -vvv)\n");
   fprintf (stderr, "   -V | --Version         Processing version number [default=%d]\n", process_get_version());
   fprintf (stderr, "   -h | --help            Print this usage message\n");
   exit (EXIT_SUCCESS);
}

static int read_config_file (const char *config_file,
                             config_t *cfg, Control_Type *ctrl)
{
   config_setting_t *setting;
   const char *template_dir;

   if (0 == config_read_file (cfg, config_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: Reading %s:%d - %s",
                     __func__, config_error_file(cfg),
                     config_error_line(cfg), config_error_text(cfg));
        return -1;
     }

   if (NULL == (setting = config_lookup (cfg, "metadata")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'template' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   ctrl->metadata_template_dir = NULL;

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "metadata_template_dir", &template_dir))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "metadata template path not found: skipping template expansion");
     }
   else if (NULL == (ctrl->metadata_template_dir = expand_string (template_dir)))
     {
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "hk_glob_pattern", &ctrl->instr_glob))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading hk_glob_pattern in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (setting = config_lookup (cfg, "calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'calibration' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string
       (setting, "badpix_file", &ctrl->bpix_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading badpix_file in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int check_env (Control_Type *ctrl)
{
   const char *env = "SDPC_DIAGNOSTIC_INDEX";
   char *val;

   ctrl->diagnostic_index = -1;

   if (NULL == (val = getenv (env)))
     return 0;

   /* setting it to "OFF" means "not set" */
   if (0 == strcasecmp (val, "OFF"))
     return 0;

   if (1 != sscanf (val, "%d", &ctrl->diagnostic_index))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading environment variable %s = %s",
                     __func__, env, val);
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s = %d", env, ctrl->diagnostic_index);

   return 0;
}

int main (int argc, char **argv)
{
   const char appname[] = "L0_ccd";
   char *config_file = "l0_ccd.cfg";
   config_t cfg;
   Control_Type ctrl = {0};
   int log_level = 0;
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"config",  required_argument, 0, 'c'},
        {"bpix",    required_argument, 0, 'b'},
        {"dark",    required_argument, 0, 'd'},
        {"solar",   required_argument, 0, 's'},
        {"instr",   required_argument, 0, 'i'},
        {"trend",   required_argument, 0, 't'},
        {"output",  required_argument, 0, 'o'},
        {"num",     required_argument, 0, 'n'},
        {"Version", required_argument, 0, 'V'},
        {"verbose", no_argument,       0, 'v'},
	{"help",    no_argument,       0, 'h'},
        {0,0,0,0}
     };

   ctrl.limit_num_granules = INT_MAX;
   ctrl.pge_version_string = L0CCD_VERSION_STRING;

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);

   config_init (&cfg);

   /* Try reading the default config file, but if it doesn't exist,
    * keep going in case there's a config file on the command line */
   if (0 == access (config_file, F_OK | R_OK))
     {
        if (-1 == read_config_file (config_file, &cfg, &ctrl))
          goto return_status;
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hvb:c:d:s:i:o:n:t:V:", long_options, &option_index);
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
             if (-1 == read_config_file (config_file, &cfg, &ctrl))
               goto return_status;
             break;
           case 'b': ctrl.bpix_file = optarg;
             break;
           case 'd': ctrl.dark_file = optarg;
             break;
           case 'h':
	     usage();
             break;
           case 's': ctrl.irr_file = optarg;
             break;
           case 'i': ctrl.instr_status_file = optarg;
             break;
           case 'o': ctrl.output_file = optarg;
             break;
           case 'n':
             if (1 != sscanf (optarg, "%u", &ctrl.limit_num_granules))
	       usage();
             break;
           case 't': ctrl.trend_file = optarg;
             break;
           case 'V':
               {
                  int version;
                  if (1 != sscanf (optarg, "%d", &version))
                    usage();
                  process_set_version (version);
               }
                  break;
           case 'v':
             log_level++;
          }
     }

   if ((optind == argc) || (ctrl.output_file == NULL))
     usage();

   ctrl.input_file = argv[optind++];

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   (void) tio_set_cmdline (argc, argv);

   (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);

   check_env (&ctrl);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "start %s", ctrl.input_file);

   if (-1 == process_inputs (&cfg, &ctrl))
     goto return_status;

   status = EXIT_SUCCESS;
return_status:
   config_destroy (&cfg);
   FREE(ctrl.metadata_template_dir);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "status=%d, finished %s",
              status, ctrl.input_file);
   tell_close();

   return status;
}
