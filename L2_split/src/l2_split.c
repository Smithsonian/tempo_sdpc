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
#include "process.h"

static void usage (void)
{
   fprintf (stderr, "Usage: L2_split [options] file1 file2 ...\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help            print this usage message\n");
   fprintf (stderr, "   -c | --config FILE     configuration file\n");
   fprintf (stderr, "   -v | --verbose lev     logging level\n");
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

int main (int argc, char **argv)
{
   const char appname[] = "L2_split";
   char *config_file = "l2_split.cfg";
   config_t cfg;
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"config",  required_argument, 0, 'c'},
        {"help",    no_argument,       0, 'h'},
        {"verbose", required_argument, 0, 'v'},
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
        int c = getopt_long (argc, argv, "hc:v:", long_options, &option_index);
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
           case 'h':
             usage();
             break;
           case 'v':
             {
                int log_level;
                if (1 == sscanf (optarg, "%d", &log_level))
                  (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);
             }
          }
     }

   if (optind == argc)
     usage();

   tio_set_cmdline (argc, argv);

   if (0 != process_files (&cfg, argc-optind, &argv[optind]))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s:  processing files", __func__);
        goto return_status;
     }

   status = EXIT_SUCCESS;
return_status:
   config_destroy (&cfg);
   tell_close();

   return status;
}
