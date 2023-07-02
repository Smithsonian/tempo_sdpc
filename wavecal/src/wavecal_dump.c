#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <ioclib.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal_dump FILE\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -b | --band NAME      Wavelength band to dump (group name)\n");
   fprintf (stderr, "   -s | --step STEP      Mirror step index\n");
   fprintf (stderr, "   -x | --xtrack XTRACK  Cross-track index\n");
   fprintf (stderr, "   -h | --help           Print this usage message\n");
   exit (EXIT_SUCCESS);
}

static int perform_dump (int grp, int step, int xtrack)
{
   TIO_Var_Info_Type info = {0};
   int start[3], count[3], i, num_waves;
   float *waves = NULL;
   float *nominal_waves = NULL;
   int status = -1;

   if (0 != TIO_inq_var (grp, TEMPO_VAR_WAVELEN_NOMINAL, &info))
     goto close_and_return;

   num_waves = info.dimlens[1];

   if (NULL == (waves = (float *)MALLOC (2 * num_waves * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }
   nominal_waves = waves + num_waves;

   start[0] = step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = num_waves;

   if (0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELENGTH, start, count, TIO_FLOAT, waves))
     goto close_and_return;

   start[0] = xtrack;
   start[1] = 0;
   start[2] = 0;

   count[0] = 1;
   count[1] = num_waves;
   count[2] = 0;

   if (0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELEN_NOMINAL, start, count, TIO_FLOAT, nominal_waves))
     goto close_and_return;

   fprintf (stdout, "index     wave  nominal wave-nominal\n");
   for (i = 0; i < num_waves; i++)
     {
        fprintf (stdout, "%5d %8.3f %8.3f %8.3f\n",
                 i, waves[i], nominal_waves[i], waves[i] - nominal_waves[i]);
     }

   status = 0;
close_and_return:
   FREE(waves);

   return status;
}

int main (int argc, char **argv)
{
   const char appname[] = "wavecal_dump";
   int status = EXIT_FAILURE;
   const char *file = NULL;
   const char *grpname = "band_290_490_nm";
   int ncid, grp;
   int step = 0;
   int xtrack = 0;
   static struct option long_options[] =
     {
        {"help",   no_argument, 0, 'h'},
        {"step",   required_argument, 0, 's'},
        {"xtrack", required_argument, 0, 'x'},
        {0,0,0,0}
     };

   if (argc == 0)
     usage();

   tell_open (appname, -1, 0);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hs:x:", long_options, &option_index);
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
           case 'b':
             grpname = optarg;
             break;
           case 's':
             if (1 != sscanf (optarg, "%d", &step))
               usage();
             break;
           case 'x':
             if (1 != sscanf (optarg, "%d", &xtrack))
               usage();
             break;
          }
     }

   if (optind == argc)
     usage();

   file = argv[optind];

   tio_set_cmdline (argc, argv);

   if (0 != TIO_open (file, NC_WRITE, &ncid))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, grpname, &grp))
     goto return_status;

   if (0 != perform_dump (grp, step, xtrack))
     goto return_status;

   status = 0;
return_status:
   (void) TIO_close (ncid);
   tell_close();
   return status;
}
