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
#include "control.h"
#include "process.h"

static void usage (void)
{
   fprintf (stderr, "Usage: L0_ccd [options] <input-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -d | --dark FILE       input corrected dark current file\n");
   fprintf (stderr, "   -i | --instr FILE      instrument telemetry points file;\n");
   fprintf (stderr, "                          provide a list using FILE=@DIR/paths.lis\n");
   fprintf (stderr, "   -o | --output FILE     output file\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -b | --bpix FILE       bad pixel file\n");
   fprintf (stderr, "   -c | --config FILE     configuration file\n");
   fprintf (stderr, "   -n | --num N           process <= N exposure records \n");
   fprintf (stderr, "   -v | --verbose lev     logging level\n");
   exit (EXIT_SUCCESS);
}

static int read_config_file (const char *config_file,
                             config_t *cfg, Control_Type *ctrl)
{
   config_setting_t *setting;

   if (0 == config_read_file (cfg, config_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: Reading %s:%d - %s",
                     __func__, config_error_file(cfg),
                     config_error_line(cfg), config_error_text(cfg));
        return -1;
     }

   if (NULL == (setting = config_lookup (cfg, "ccd_calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_calibration in param file: %s",
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

   if (CONFIG_TRUE != config_setting_lookup_string
       (setting, "hk_glob_pattern", &ctrl->instr_glob))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading hk_glob_pattern in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

int main (int argc, char **argv)
{
   const char appname[] = "L0_ccd";
   char *config_file = "l0_ccd.cfg";
   config_t cfg;
   Control_Type ctrl = {0};
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"config",  optional_argument, 0, 'c'},
        {"bpix",    optional_argument, 0, 'b'},
        {"dark",    required_argument, 0, 'd'},
        {"instr",   required_argument, 0, 'i'},
        {"output",  required_argument, 0, 'o'},
        {"num",     optional_argument, 0, 'n'},
        {"verbose", optional_argument, 0, 'v'},
        {0,0,0,0}
     };

   ctrl.limit_num_granules = INT_MAX;

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
        int c = getopt_long (argc, argv, "b:c:d:i:o:v:n:", long_options, &option_index);
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
           case 'i': ctrl.instr_status_file = optarg;
             break;
           case 'o': ctrl.output_file = optarg;
             break;
           case 'n':
             if (1 != sscanf (optarg, "%u", &ctrl.limit_num_granules))
               usage();
           case 'v':
             {
                int log_level;
                if (1 == sscanf (optarg, "%d", &log_level))
                  (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);
             }
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

   tell_vlog (TELL_MSGTYPE_INFO, 0, "start %s", ctrl.input_file);

   if (-1 == process_inputs (&cfg, &ctrl))
     goto return_status;

   status = EXIT_SUCCESS;
return_status:
   config_destroy (&cfg);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "status=%d, finished %s",
              status, ctrl.input_file);
   tell_close();

   return status;
}
