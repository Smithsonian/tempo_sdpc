#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_errno.h>
#include <libconfig.h>
#include <tell.h>
#include <tio.h>

#include "config.h"
#include "wavecal.h"

typedef struct
{
   double *rad;
   double *rad_err;
   double *y;
   size_t n;
}
Radiance_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal [options] <input-file>\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -o | --outpar FILE     output file for wavelength parameters\n");
   fprintf (stderr, "   -s | --silent          turn off file output (for timing)\n");
   fprintf (stderr, "   -x | --xtrack N        cross-track pixel index, 0 is northernmost\n");
   exit (EXIT_SUCCESS);
}

static void free_radiance (Radiance_Type *r)
{
   if (r == NULL)
     return;
   FREE(r->rad);
}

static int alloc_radiance (Radiance_Type *r, size_t n)
{
   if (NULL == (r->rad = (double *)MALLOC (3 * n * sizeof(double))))
     return -1;
   r->n = n;
   r->y = r->rad + n;
   r->rad_err = r->y + n;
   return 0;
}

static int read_radiance (int ncid, int xtrack, Radiance_Type *r)
{
   int start[2], count[2];
   double rmax, a, sn_max=2500.0;
   size_t i;
   int status = -1;

   start[0] = xtrack;
   start[1] = 0;

   count[0] = 1;
   count[1] = r->n;

   if (0 != TIO_get_var_section (ncid, "RadianceObserved", start, count, TIO_DOUBLE,
                                 r->rad))
     goto return_error;

   if (0 != TIO_get_var_section (ncid, "Wavelength", start, count, TIO_DOUBLE,
                                 r->y))
     goto return_error;

   /* To fake some plausible uncertainty values,
    * assume 1) radiance is proportional to counts
    *        2) uncertainties are Poisson, so sigma=sqrt(N)
    *        3) peak radiance has signal-to-noise ratio, sn_max
    * Therefore:
    *    sn_max = sqrt(nmax)   and rmax = a*nmax
    *   => a = rmax/sn_max^2
    * so that for any r, n(r) = r/a,
    * and r_err = a*sqrt(n(r)) = a*sqrt(r/a) = sqrt(a*r)
    */
   rmax = r->rad[0];
   for (i = 1; i < r->n; i++)
     {
        if (r->rad[i] > rmax) rmax = r->rad[i];
     }
   a = rmax / (sn_max * sn_max);
   for (i = 0; i < r->n; i++)
     {
        r->rad_err[i] = sqrt(a * r->rad[i]);
     }

   status = 0;
return_error:
   return status;
}

static int write_resid (size_t num_wave,
                        const double *wave0,
                        const double *model, const double *spec,
                        const double *resid)
{
   FILE *fp;
   size_t i;

   if (NULL == (fp = fopen ("resid.dat", "w")))
     return -1;

   for (i = 0; i < num_wave; i++)
     {
        fprintf (fp, "%5ld %15.6e %15.6f %15.6e %15.6e\n",
                 i, spec[i], wave0[i], model[i], resid[i]);
     }

   return fclose (fp);
}

static int write_term_info (const Wavecal_Term_Info_Type *info, const double *wave)
{
#define BUFSIZE 256
   char filename[BUFSIZE];
   FILE *fp;
   size_t i;

   if (info == NULL)
     return -1;

   snprintf (filename, BUFSIZE, "term_%s.dat", info->name);

   if (NULL == (fp = fopen (filename, "w")))
     {
        fprintf (fp, "*** Error opening %s for writing\n", filename);
        return -1;
     }

   for (i = 0; i < info->num_values; i++)
     {
        fprintf (fp, "%4d %15.6f %15.6e\n", i, wave[i], info->value[i]);
     }

   (void) fclose (fp);

   return 0;
}

static int write_fit_results (FILE *fp, int xtrack, int num_wave,
                              const Wavecal_Type *wct,
                              const Wavecal_Result_Type *wavecal_result)
{
   Wavecal_Term_Info_Type info = {0};
   size_t i;
   int nth;

   fprintf (fp, "%4d %12.4e %4d ", xtrack, wavecal_result->bestnorm,
            wavecal_result->nfev);
   for (i = 0; i < wavecal_result->num_wave_params; i++)
     {
        fprintf (fp, "%15.9e ", wavecal_result->wave_params[i]);
     }
   fprintf (fp, "\n");

   write_resid (num_wave, wavecal_result->wave,
                wavecal_result->model, wavecal_result->spec_scaled,
                wavecal_result->residuals);
   nth = 0;
   do
     {
        nth = wavecal_query_term (wct, nth, &info);
        write_term_info (&info, wavecal_result->wave);
     }
   while (nth > 0);

   return 0;
}

int main (int argc, char **argv)
{
   const char appname[] = "wavecal_driver";
   const char *config_file = "l0_ccd.cfg";
   const char *file = NULL;
   const char *params_outfile = NULL;
   FILE *fp = stderr;
   config_t cfg;
   Radiance_Type rad = {0};
   int status = EXIT_FAILURE;
   Wavecal_Type *wct = NULL;
   Wavecal_Config_Type wavecal_config = {0};
   Wavecal_Result_Type wavecal_result = {0};
   double *y0 = NULL;
   double nan_value = nan("");
   int ncid, xtrack = -1;
   int beg_xtrack, end_xtrack;
   int silent = 0;
   size_t i, len;
   static struct option long_options[] =
     {
        {"config",  optional_argument, 0, 'c'},
        {"outpar",  optional_argument, 0, 'o'},
        {"silent",  optional_argument, 0, 's'},
        {"xtrack",  optional_argument, 0, 'x'},
        {0,0,0,0}
     };

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);
   gsl_set_error_handler_off ();

   config_init (&cfg);

   /* Try reading the default config file, but if it doesn't exist,
    * keep going in case there's a config file on the command line */
   if (0 == access (config_file, F_OK | R_OK))
     {
        if (0 == config_read_file (&cfg, config_file))
          goto return_status;
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "c:x:o:s", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d??", c);
             goto return_status;
             break;
           case 'c': config_file = optarg;
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (0 == config_read_file (&cfg, config_file))
               goto return_status;
             break;
           case 'o': params_outfile = optarg;
             break;
           case 's': silent++;
             break;
           case 'x':
             if (1 != sscanf (optarg, "%d", &xtrack))
               usage();
             break;
          }
     }

   if (optind == argc)
     usage();

   file = argv[optind++];

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
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "Reading %s:%d - %s",
                     config_error_file(&cfg),
                     config_error_line(&cfg), config_error_text(&cfg));
        goto return_status;
     }

   wavecal_config.fill_value = nan_value;

   len = 1001;

   if (alloc_radiance (&rad, len))
     goto return_status;

   if (NULL == (y0 = (double *)malloc (len * sizeof(double))))
     goto return_status;

   for (i = 0; i < len; i++)
     {
        y0[i] = 290.0 + i*0.2;
     }

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if (NULL == (wct = wavecal_open (&cfg, rad.n, 1)))
     goto return_status;

   if (params_outfile)
     {
        if (NULL == (fp = fopen (params_outfile, "w")))
          {
             fprintf (stderr, "*** Error: opening file for writing: %s\n", params_outfile);
             goto return_status;
          }
     }

   /* default test xtrack = 196;  test file contains 301 */
   if (xtrack < 0)
     {
        beg_xtrack = 0;
        end_xtrack = 301;
     }
   else
     {
        beg_xtrack = xtrack;
        end_xtrack = xtrack+1;
     }

   for (xtrack = beg_xtrack; xtrack < end_xtrack; xtrack++)
     {
        int nth;
        if (read_radiance (ncid, xtrack, &rad))
          goto return_status;

        if (0 != wavecal_fit (wct, xtrack, y0, rad.rad, rad.rad_err,
                              &wavecal_config, &wavecal_result))
          goto return_status;

        if (!silent)
          {
             (void) write_fit_results (fp, xtrack, rad.n, wct, &wavecal_result);
          }
     }

   status = EXIT_SUCCESS;
return_status:
   FREE(y0);
   TIO_close (ncid);
   free_radiance (&rad);
   config_destroy (&cfg);
   tell_close();
   wavecal_close (wct);

   if ((fp != NULL) && (params_outfile != NULL))
     {
        (void) fclose (fp);
     }

   return status;
}
