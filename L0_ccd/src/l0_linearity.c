#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <tell.h>

#include "config.h"
#include "linearity.h"

static void usage (void)
{
   fprintf (stderr, "Usage: L0_linearity [options] <input-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -o | --output FILE     Output file\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -v | --verbose         Verbosity (more instances means more verbose, e.g. -vvv)\n");
   fprintf (stderr, "   -h | --help            Print this usage message\n");
   exit (EXIT_SUCCESS);
}

int main (int argc, char **argv)
{
   const char appname[] = "L0_linearity";
   const char *input_file = NULL;
   const char *output_file = NULL;
   int log_level = 0;
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"output",  required_argument, 0, 'o'},
        {"verbose", no_argument,       0, 'v'},
	{"help",    no_argument,       0, 'h'},
        {0,0,0,0}
     };

   if (argc < 2)
     usage();

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hvb:c:d:i:o:n:V:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d??", c);
             goto return_status;
             break;
           case 'h':
	     usage();
             break;
           case 'o': output_file = optarg;
             break;
           case 'v':
             log_level++;
             break;
          }
     }

   if ((optind == argc) || (output_file == NULL))
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

   tell_open (appname, -1, 0);

   (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "start %s", input_file);

   if (-1 == derive_linearity (input_file, output_file))
     goto return_status;

   status = EXIT_SUCCESS;
return_status:

   tell_vlog (TELL_MSGTYPE_INFO, 0, "status=%d, finished %s",
              status, input_file);
   tell_close();

   return status;
}
