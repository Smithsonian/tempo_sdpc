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
#include <tio_template.h>

#include "config.h"
#include "wavecal.h"

#define PIXEL_SIZE_NANOMETERS   0.2
#define MIN_WAVELENGTH_UV     288.0
#define MIN_WAVELENGTH_VIS    536.8

typedef struct
{
   double *spec;
   double *spec_err;
   size_t n;
}
Spectrum_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal [options] <input-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -g | --group NAME          name of netCDF4 file group containing spectra\n");
   fprintf (stderr, "   -S | --yStart WAVELENGTH   starting wavelength\n");
   fprintf (stderr, "   -D | --yDelta DELTA        wavelength step per pixel\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help            print this usage message\n");
   fprintf (stderr, "   -m | --mirror STEP     mirror step index\n");
   fprintf (stderr, "   -x | --xtrack N        cross-track pixel index, 0 is northernmost\n");
   fprintf (stderr, "   -o | --outpar FILE     output file for wavelength parameters\n");
   fprintf (stderr, "   -v | --verbose         turn on verbose output\n");
   exit (EXIT_SUCCESS);
}

static void free_spectrum (Spectrum_Type *r)
{
   if (r == NULL)
     return;
   FREE(r->spec);
}

static int alloc_spectrum (Spectrum_Type *r, size_t n)
{
   if (NULL == (r->spec = (double *)MALLOC (2 * n * sizeof(double))))
     return -1;
   r->n = n;
   r->spec_err = r->spec + n;
   return 0;
}

static void fake_spectrum_errors (Spectrum_Type *r, double snr_max)
{
   double rmax, a;
   size_t i;

   /* To fake some plausible uncertainty values,
    * assume 1) spectral radiance/irradiance is proportional to counts
    *        2) uncertainties are Poisson, so sigma=sqrt(N)
    *        3) peak has signal-to-noise ratio, snr_max
    * Therefore:
    *    snr_max = sqrt(nmax)   and rmax = a*nmax
    *   => a = rmax/snr_max^2
    * so that for any r, n(r) = r/a,
    * and r_err = a*sqrt(n(r)) = a*sqrt(r/a) = sqrt(a*r)
    */
   rmax = r->spec[0];
   for (i = 1; i < r->n; i++)
     {
        if (r->spec[i] > rmax) rmax = r->spec[i];
     }
   a = rmax / (snr_max * snr_max);
   for (i = 0; i < r->n; i++)
     {
        r->spec_err[i] = sqrt(a * r->spec[i]);
     }
}

static int read_spectrum (int ncid, int step, int xtrack, int is_irradiance,
                          Spectrum_Type *r)
{
   const char *var_spec;
   const char *var_err;
   int start[3], count[3];
   int status = -1;

   if (is_irradiance)
     {
        var_spec = TEMPO_VAR_IRRADIANCE;
        var_err = TEMPO_VAR_IRRADIANCE_ERROR;
     }
   else
     {
        var_spec = TEMPO_VAR_RADIANCE;
        var_err = TEMPO_VAR_RADIANCE_ERROR;
     }

   start[0] = step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = r->n;

   if (0 != TIO_get_var_section (ncid, var_spec, start, count, TIO_DOUBLE, r->spec))
     goto return_error;

   tell_push_queue();
   status = TIO_get_var_section (ncid, var_err, start, count, TIO_DOUBLE, r->spec_err);
   tell_pop_queue(1);
   if (status)
     {
        double snr_max = 2500.0;
        fprintf (stderr, "*** Faking spectrum errors assuming max(SNR) = %g\n", snr_max);
        fake_spectrum_errors (r, snr_max);
     }

   status = 0;
return_error:
   return status;
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
        fprintf (fp, "%4ld %15.6f %15.6e\n", i, wave[i], info->value[i]);
     }

   (void) fclose (fp);

   return 0;
}

static void write_fit (const Wavecal_Result_Type *r, size_t num_values)
{
   FILE *fp = NULL;
   size_t i;

   if (NULL == (fp = fopen ("result.dat", "w")))
     {
        fprintf (stderr, "*** %s: Error writing results\n", __func__);
        return;
     }
   for (i = 0; i < num_values; i++)
     {
        fprintf (fp, "%4ld %9.4f %15.4e %15.4e %15.4e\n", i,
                 r->wave[i],
                 r->spec_scaled[i],
                 r->model[i],
                 r->residuals[i]);
     }
   (void) fclose (fp);
}

static int write_fit_details (FILE *fp, int xtrack,
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

   nth = 0;
   do
     {
        nth = wavecal_query_term (wct, nth, &info);
        write_term_info (&info, wavecal_result->wave);
     }
   while (nth > 0);

   (void) wavecal_query_term (wct, 0, &info);
   write_fit (wavecal_result, info.num_values);

   return 0;
}

static int create_result_group (int parent_grp, const char *grp_name,
                                const TIO_Var_Info_Type *spectrum_info,
                                const Wavecal_Result_Type *wavecal_result,
                                int *pgrp)
{
   int grp, varid, param_dimids[3];
   size_t params_dimlen;

   if (0 != TIO_def_grp (parent_grp, grp_name, &grp))
     return -1;

   memcpy ((char *)param_dimids, (char *)spectrum_info->dimids,
           3 * sizeof (int));
   params_dimlen = wavecal_result->num_wave_params;
   if (0 != TIO_def_dim (grp, "params", params_dimlen, &param_dimids[2]))
     return -1;

   if ((0 != TIO_def_var (grp, "wavelength", TIO_FLOAT,
                          spectrum_info->ndims,
                          spectrum_info->dimids, &varid))
       || (0 != TIO_def_var (grp, "model", TIO_FLOAT,
                             spectrum_info->ndims,
                             spectrum_info->dimids, &varid))
       || (0 != TIO_def_var (grp, "residuals", TIO_FLOAT,
                             spectrum_info->ndims,
                             spectrum_info->dimids, &varid))
       || (0 != TIO_def_var (grp, "params", TIO_FLOAT,
                             spectrum_info->ndims,
                             param_dimids, &varid))
      )
     {
        return -1;
     }

   *pgrp = grp;

   return 0;
}

static int write_results (int parent_grp, const TIO_Var_Info_Type *spectrum_info,
                          int step, int xtrack, int num_wave,
                          const Wavecal_Type *wct,
                          const Wavecal_Result_Type *wavecal_result)
{
   const char grp_name[] = "wavecal_test";
   int grp, status, start[3], count[3];

   (void) wct;

   tell_push_queue ();
   status = TIO_inq_grp (parent_grp, grp_name, &grp);
   tell_pop_queue (1);
   if (status)
     {
        if (0 != create_result_group (parent_grp, grp_name, spectrum_info,
                                      wavecal_result, &grp))
          return -1;
     }

   start[0] = step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = num_wave;

   if ((0 != TIO_put_var_section (grp, "wavelength", start, count, TIO_DOUBLE,
                                  wavecal_result->wave))
       || (0 != TIO_put_var_section (grp, "model", start, count, TIO_DOUBLE,
                                  wavecal_result->model))
       || (0 != TIO_put_var_section (grp, "residuals", start, count, TIO_DOUBLE,
                                  wavecal_result->residuals)))
     {
        return -1;
     }

   count[2] = wavecal_result->num_wave_params;
   if (0 != TIO_put_var_section (grp, "params", start, count, TIO_DOUBLE,
                                 wavecal_result->wave_params))
     return -1;

   return 0;
}

int main (int argc, char **argv)
{
   const char appname[] = "wavecal_driver";
   const char *config_file = "l0_ccd.cfg";
   const char *file = NULL;
   const char *params_outfile = NULL;
   const char *group_name = NULL;
   FILE *fp = stderr;
   config_t cfg;
   Spectrum_Type spec = {0};
   int status = EXIT_FAILURE;
   Wavecal_Type *wct = NULL;
   TIO_Var_Info_Type spectrum_info = {0};
   Wavecal_Config_Type wavecal_config = {0};
   Wavecal_Result_Type wavecal_result = {0};
   double *y0 = NULL;
   double nan_value = nan("");
   double y_start = nan_value;
   double y_delta = nan_value;
   int ncid, grp, step = -1, xtrack = -1;
   int beg_xtrack, end_xtrack;
   int beg_step, end_step;
   int is_irradiance = 0;
   int verbose = 0;
   size_t i, len;
   static struct option long_options[] =
     {
        {"help",    optional_argument, 0, 'h'},
        {"config",  optional_argument, 0, 'c'},
        {"outpar",  optional_argument, 0, 'o'},
        {"xtrack",  optional_argument, 0, 'x'},
        {"mirror",  optional_argument, 0, 'm'},
        {"verbose", optional_argument, 0, 'v'},
        {"yDelta",  optional_argument, 0, 'D'},
        {"yStart",  optional_argument, 0, 'S'},
        {"group",   required_argument, 0, 'g'},
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
        int c = getopt_long (argc, argv, "hm:D:S:c:g:x:o:v", long_options, &option_index);
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
           case 'h':
             usage();
             break;
           case 'g': group_name = optarg;
             break;
           case 'o': params_outfile = optarg;
             break;
           case 'v': verbose++;
             break;

           case 'S':
             if (1 != sscanf (optarg, "%le", &y_start))
               usage();
             break;
           case 'D':
             if (1 != sscanf (optarg, "%le", &y_delta))
               usage();
             break;
           case 'x':
             if (1 != sscanf (optarg, "%d", &xtrack))
               usage();
             break;
           case 'm':
             if (1 != sscanf (optarg, "%d", &step))
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

   if (group_name == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "group name not specified");
        usage();
     }

   if (isnan(y_delta))
     {
        y_delta = PIXEL_SIZE_NANOMETERS;
     }

   if (isnan(y_start))
     {
        if (0 == strcasecmp (group_name, TEMPO_BAND_NAME_UV))
          y_start = MIN_WAVELENGTH_UV;
        else if (0 == strcasecmp (group_name, TEMPO_BAND_NAME_VIS))
          y_start = MIN_WAVELENGTH_VIS;
        else
          usage();
     }

   wavecal_config.fill_value = nan_value;

   if (0 != TIO_open (file, NC_WRITE, &ncid))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, group_name, &grp))
     goto return_status;

   is_irradiance = -1;
   tell_push_queue();
   if (0 == TIO_inq_var (grp, TEMPO_VAR_RADIANCE, &spectrum_info))
     {
        is_irradiance = 0;
     }
   else if (0 == TIO_inq_var (grp, TEMPO_VAR_IRRADIANCE, &spectrum_info))
     {
        is_irradiance = 1;
     }
   tell_pop_queue(1);
   if (is_irradiance < 0)
     {
        fprintf (stderr, "*** unsupported file type: %s\n", file);
        goto return_status;
     }

   /* expected dimensions are: [mirror_step, xtrack, spectral_channel] */
   len = spectrum_info.dimlens[2];

   if (alloc_spectrum (&spec, len))
     goto return_status;

   if (NULL == (y0 = (double *)malloc (len * sizeof(double))))
     goto return_status;

   for (i = 0; i < len; i++)
     {
        y0[i] = y_start + i*y_delta;
     }

   if (NULL == (wct = wavecal_open (&cfg, group_name, spec.n, is_irradiance)))
     goto return_status;

   if (params_outfile)
     {
        if (NULL == (fp = fopen (params_outfile, "w")))
          {
             fprintf (stderr, "*** Error: opening file for writing: %s\n", params_outfile);
             goto return_status;
          }
     }

   if (xtrack < 0)
     {
        beg_xtrack = 0;
        end_xtrack = spectrum_info.dimlens[1];
     }
   else
     {
        beg_xtrack = xtrack;
        end_xtrack = xtrack+1;
     }

   if (step < 0)
     {
        beg_step = 0;
        end_step = spectrum_info.dimlens[0];
     }
   else
     {
        beg_step = step;
        end_step = step+1;
     }

   for (step = beg_step; step < end_step; step++)
     {
        for (xtrack = beg_xtrack; xtrack < end_xtrack; xtrack++)
          {
             if (read_spectrum (grp, step, xtrack, is_irradiance, &spec))
               goto return_status;

             if (wavecal_fit (wct, xtrack, y0, spec.spec, spec.spec_err,
                              &wavecal_config, &wavecal_result))
               goto return_status;

             if (write_results (grp, &spectrum_info, step, xtrack, spec.n, wct,
                                &wavecal_result))
               goto return_status;

             if (verbose)
               {
                  fprintf (stderr, "%4d %4d %12.4e %4d\n", step, xtrack,
                           wavecal_result.bestnorm, wavecal_result.nfev);

                  if (params_outfile)
                    {
                       (void) write_fit_details (fp, xtrack, wct, &wavecal_result);
                    }
               }
          }
     }

   status = EXIT_SUCCESS;
return_status:
   FREE(y0);
   TIO_close (ncid);
   free_spectrum (&spec);
   config_destroy (&cfg);
   tell_close();
   wavecal_close (wct);

   if ((fp != NULL) && (params_outfile != NULL))
     {
        (void) fclose (fp);
     }

   return status;
}
