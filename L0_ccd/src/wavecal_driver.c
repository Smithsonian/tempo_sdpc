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

#define WAVECAL_PARAM_NAME      "wavecal_params"
#define WAVECAL_PARAM_DIM_NAME  "wavecal_par"

typedef struct
{
   double *spec;
   double *spec_err;
   unsigned int *pixel_quality_flag;
   size_t n;
}
Spectrum_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal_driver [options] <input-file> <output-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -g | --group NAME          name of netCDF4 file group containing spectra\n");
   fprintf (stderr, "   -S | --yStart WAVELENGTH   starting wavelength\n");
   fprintf (stderr, "   -D | --yDelta DELTA        wavelength step per pixel\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help                print this usage message\n");
   fprintf (stderr, "   -d | --debug               write diagnostic information to output file\n");
   fprintf (stderr, "   -b | --block i:num         mirror step blocking specification\n");
   fprintf (stderr, "   -m | --mirror STEP         mirror step index\n");
   fprintf (stderr, "   -x | --xtrack N            cross-track pixel index, 0 is northernmost\n");
   fprintf (stderr, "   -w | --wavepar FILE        output file for wavelength parameters\n");
   fprintf (stderr, "   -c | --config FILE         path to configuration file\n");
   fprintf (stderr, "   -v | --verbose             turn on verbose output\n");
   exit (EXIT_SUCCESS);
}

static void free_spectrum (Spectrum_Type *r)
{
   if (r == NULL)
     return;
   FREE(r->spec);
   FREE(r->pixel_quality_flag);
}

static int alloc_spectrum (Spectrum_Type *r, size_t n)
{
   if ((NULL == (r->spec = (double *)MALLOC (2 * n * sizeof(double))))
       || (NULL == (r->pixel_quality_flag = (unsigned int *)MALLOC (n * sizeof(unsigned int)))))
     return -1;
   memset ((char *) r->spec, 0, 2*n*sizeof(double));
   memset ((char *) r->pixel_quality_flag, 0, n * sizeof(unsigned int));
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

   if ((0 != TIO_get_var_section (ncid, var_spec, start, count, TIO_DOUBLE, r->spec))
       || (0 != TIO_get_var_section (ncid, TEMPO_VAR_PQF, start, count, TIO_UINT, r->pixel_quality_flag)))
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

static int *read_inr_quality_flag (int grp, size_t *dimlens)
{
   int *inr_quality_flag = NULL;
   int start[2], count[2];
   size_t len = dimlens[0] * dimlens[1];

   if (NULL == (inr_quality_flag = (int *)MALLOC (len * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)inr_quality_flag, 0, len * sizeof(int));

   start[0] = 0;
   start[1] = 0;
   count[0] = dimlens[0];
   count[1] = dimlens[1];

   if (0 != TIO_get_var_section (grp, TEMPO_VAR_INRQF, start, count, TIO_INT, inr_quality_flag))
     {
        FREE(inr_quality_flag);
        return NULL;
     }

   return inr_quality_flag;
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
                              const double *wave_params,
                              const Wavecal_Result_Type *wavecal_result)
{
   Wavecal_Term_Info_Type info = {0};
   size_t i;
   int nth;

   fprintf (fp, "%4d %12.4e %4d %4d ", xtrack, wavecal_result->bestnorm,
            wavecal_result->nfev, wavecal_result->opt_status);
   for (i = 0; i < wavecal_result->num_wave_params; i++)
     {
        fprintf (fp, "%15.9e ", wave_params[i]);
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

static int create_result_file (const char *path, const char *group_name,
                               size_t beg_step, size_t end_step, size_t step_dimlen,
                               size_t num_xtrack, int start_pix, int num_pix,
                               size_t params_dimlen)
{
   int ncid, varid, param_dimids[3], start, count;
   size_t i, num_steps = end_step - beg_step;
   int max_num_steps = step_dimlen;
   int *steps = NULL;

   if (0 != TIO_create (path, NC_NETCDF4, &ncid))
     return -1;

   if ((0 != TIO_put_att (ncid, NC_GLOBAL, "group_name", TIO_CHAR,
                          strlen(group_name), group_name))
       || (0 != TIO_put_att (ncid, NC_GLOBAL, "mirror_step_dimlen", TIO_INT, 1, &max_num_steps)))
     goto close_and_return;

   if ((0 != TIO_def_dim (ncid, TEMPO_DIM_STEP, num_steps, &param_dimids[0]))
       || (0 != TIO_def_dim (ncid, TEMPO_DIM_XTRACK, num_xtrack, &param_dimids[1]))
       || (0 != TIO_def_dim (ncid, WAVECAL_PARAM_DIM_NAME, params_dimlen, &param_dimids[2])))
     goto close_and_return;

   if (0 != TIO_def_var (ncid, WAVECAL_PARAM_NAME, TIO_FLOAT, 3, param_dimids, &varid))
     goto close_and_return;
   if ((0 != TIO_put_att (ncid, varid, "start_spectral_channel", TIO_INT, 1, &start_pix))
       ||(0 != TIO_put_att (ncid, varid, "num_spectral_channels", TIO_INT, 1, &num_pix)))
     goto close_and_return;

   if (0 != TIO_def_var (ncid, "bestnorm", TIO_FLOAT, 2, param_dimids, &varid))
     goto close_and_return;

   if (0 != TIO_def_var (ncid, TEMPO_DIM_STEP, TIO_INT, 1, &param_dimids[0], &varid))
     goto close_and_return;

   if (NULL == (steps = (int *)MALLOC (num_steps * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   for (i = beg_step; i < end_step; i++)
     {
        steps[i-beg_step] = i;
     }

   start = 0;
   count = num_steps;

   if (0 != TIO_put_var_section (ncid, TEMPO_DIM_STEP, &start, &count, TIO_INT, steps))
     goto close_and_return;

   FREE(steps);
   return ncid;

close_and_return:
   FREE(steps);
   (void) TIO_close (ncid);
   return -1;
}

static int write_result (int ncid, int beg_step, int step, int xtrack,
                         const double *wave_params,
                         const Wavecal_Result_Type *wavecal_result)
{
   int start[3], count[3];

   start[0] = step - beg_step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = wavecal_result->num_wave_params;

   if ((0 != TIO_put_var_section (ncid, WAVECAL_PARAM_NAME, start, count, TIO_DOUBLE,
                                  wave_params))
       ||(0 != TIO_put_var_section (ncid, "bestnorm", start, count, TIO_DOUBLE,
                                    &wavecal_result->bestnorm)))
     return -1;

   return 0;
}

static int create_diagnostic_group (int parent_grp, const char *grp_name,
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
   if (0 != TIO_def_dim (grp, WAVECAL_PARAM_DIM_NAME, params_dimlen, &param_dimids[2]))
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
       || (0 != TIO_def_var (grp, WAVECAL_PARAM_NAME, TIO_FLOAT,
                             spectrum_info->ndims,
                             param_dimids, &varid))
       || (0 != TIO_def_var (grp, "bestnorm", TIO_FLOAT, 2,
                             param_dimids, &varid))
      )
     {
        return -1;
     }

   *pgrp = grp;

   return 0;
}

static int write_diagnostics (int parent_grp, const TIO_Var_Info_Type *spectrum_info,
                              int step, int xtrack, const Wavecal_Type *wct,
                              const double *wave_params,
                              const Wavecal_Result_Type *wavecal_result)
{
   const char grp_name[] = "wavecal_diagnostics";
   int grp, status, start[3], count[3];

   (void) wct;

   tell_push_queue ();
   status = TIO_inq_grp (parent_grp, grp_name, &grp);
   tell_pop_queue (1);
   if (status)
     {
        if (0 != create_diagnostic_group (parent_grp, grp_name, spectrum_info,
                                          wavecal_result, &grp))
          return -1;
     }

   start[0] = step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = wavecal_result->num_fit;

   if ((0 != TIO_put_var_section (grp, "wavelength", start, count, TIO_DOUBLE,
                                  wavecal_result->wave))
       || (0 != TIO_put_var_section (grp, "model", start, count, TIO_DOUBLE,
                                  wavecal_result->model))
       || (0 != TIO_put_var_section (grp, "residuals", start, count, TIO_DOUBLE,
                                  wavecal_result->residuals))
       || (0 != TIO_put_var_section (grp, "bestnorm", start, count, TIO_DOUBLE,
                                     &wavecal_result->bestnorm)))
     {
        return -1;
     }

   count[2] = wavecal_result->num_wave_params;
   if (0 != TIO_put_var_section (grp, "params", start, count, TIO_DOUBLE,
                                 wave_params))
     return -1;

   return 0;
}

int main (int argc, char **argv)
{
   const char appname[] = "wavecal_driver";
   const char *config_file = "wavecal.cfg";
   const char *input_file = NULL;
   const char *result_file = NULL;
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
   int *inr_quality_flag = NULL;
   double *y0 = NULL;
   double nan_value = nan("");
   double y_start = nan_value;
   double y_delta = nan_value;
   double *wave_params = NULL;
   int grp, ncid = 0, step = -1, xtrack = -1;
   int beg_xtrack, end_xtrack;
   int beg_step, end_step;
   int use_blocking = 0, this_block, num_blocks;
   int is_irradiance = 0;
   int verbose = 0;
   int ncid_result = 0, num_wave_params, start_pix, num_pix;
   int debug = 0;
   size_t i, step_dimlen, xtrack_dimlen, channel_dimlen;
   int fit_status_code;
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"debug",    no_argument, 0, 'd'},
        {"config",  required_argument, 0, 'c'},
        {"wavepar", required_argument, 0, 'w'},
        {"xtrack",  required_argument, 0, 'x'},
        {"block",  required_argument, 0, 'b'},
        {"mirror",  required_argument, 0, 'm'},
        {"verbose", no_argument, 0, 'v'},
        {"yDelta",  required_argument, 0, 'D'},
        {"yStart",  required_argument, 0, 'S'},
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
        int c = getopt_long (argc, argv, "hdb:m:D:S:c:g:x:w:v", long_options, &option_index);
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
           case 'd':
             debug++;
             break;
           case 'h':
             usage();
             break;
           case 'g': group_name = optarg;
             break;
           case 'w': params_outfile = optarg;
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
           case 'b':
             if (2 != sscanf (optarg, "%d:%d", &this_block, &num_blocks))
               usage();
             if (this_block < 0 || num_blocks <= this_block)
               {
                  fprintf (stderr, "*** wavecal_driver: invalid blocking specification\n");
                  goto return_status;
               }
             use_blocking++;
             break;
          }
     }

   if (optind+2 < argc)
     usage();

   input_file = argv[optind++];
   result_file = argv[optind++];

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

   if (0 != TIO_open (input_file, NC_WRITE, &ncid))
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
        fprintf (stderr, "*** unsupported file type: %s\n", input_file);
        goto return_status;
     }

   /* expected dimensions are: [mirror_step, xtrack, spectral_channel] */
   step_dimlen = spectrum_info.dimlens[0];
   xtrack_dimlen = spectrum_info.dimlens[1];
   channel_dimlen = spectrum_info.dimlens[2];

   if (alloc_spectrum (&spec, channel_dimlen))
     goto return_status;

   if (NULL == (y0 = (double *)MALLOC (channel_dimlen * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   for (i = 0; i < channel_dimlen; i++)
     {
        y0[i] = y_start + i*y_delta;
     }

   if (0 == is_irradiance)
     {
        if (NULL == (inr_quality_flag = read_inr_quality_flag (grp, spectrum_info.dimlens)))
          goto return_status;
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
        end_xtrack = xtrack_dimlen;
     }
   else
     {
        beg_xtrack = xtrack;
        end_xtrack = xtrack+1;
     }

   if (use_blocking)
     {
        int block_size = step_dimlen / num_blocks;
        int residual = step_dimlen - num_blocks * block_size;
        beg_step = this_block * block_size;
        if (residual > 0)
          {
             if (this_block < residual)
               {
                  block_size += 1;
                  beg_step += this_block;
               }
             else beg_step += residual;
          }
        end_step = beg_step + block_size;
        if (end_step > (int) step_dimlen)
          end_step = step_dimlen;
     }
   else if (step < 0)
     {
        beg_step = 0;
        end_step = step_dimlen;
     }
   else
     {
        beg_step = step;
        end_step = step+1;
     }

   num_wave_params = wavecal_num_wave_params (wct);
   if (num_wave_params <= 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid number of wavelength scale parameters: num_wave_params=%d",
                     __func__, num_wave_params);
        goto return_status;
     }

   if (NULL == (wave_params = (double *)MALLOC (num_wave_params * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   (void) wavecal_feature_window (wct, &start_pix, &num_pix);

   ncid_result = create_result_file (result_file, group_name,
                                     beg_step, end_step, step_dimlen,
                                     xtrack_dimlen,
                                     start_pix, num_pix, num_wave_params);
   if (ncid_result <= 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: problem creating result file: %s",
                     __func__, result_file);
        goto return_status;
     }

   for (step = beg_step; step < end_step; step++)
     {
        int *inrqf_step = NULL;

        if (0 == is_irradiance)
          {
             inrqf_step = inr_quality_flag + step * xtrack_dimlen;
          }

        for (xtrack = beg_xtrack; xtrack < end_xtrack; xtrack++)
          {
             if ((inrqf_step != NULL) && (inrqf_step[xtrack] != 0))
               continue;

             if (read_spectrum (grp, step, xtrack, is_irradiance, &spec))
               goto return_status;

             fit_status_code = wavecal_fit (wct, xtrack, y0, spec.spec, spec.spec_err,
                                            spec.pixel_quality_flag, &wavecal_config,
                                            wave_params, &wavecal_result);
             if (fit_status_code == WAVECAL_FIT_ERROR)
               goto return_status;

             if (write_result (ncid_result, beg_step, step, xtrack, wave_params, &wavecal_result))
               goto return_status;

             if (debug)
               {
                  if (write_diagnostics (grp, &spectrum_info, step, xtrack, wct, wave_params,
                                         &wavecal_result))
                    goto return_status;
               }

             if (verbose)
               {
                  fprintf (stderr, "%4d %4d %12.4e %4d %4d\n", step, xtrack,
                           wavecal_result.bestnorm, wavecal_result.nfev,
                           wavecal_result.opt_status);

                  if (params_outfile)
                    {
                       (void) write_fit_details (fp, xtrack, wct, wave_params, &wavecal_result);
                    }
               }
          }
     }

   status = EXIT_SUCCESS;
return_status:
   FREE(y0);
   FREE(inr_quality_flag);
   FREE(wave_params);
   if (ncid) TIO_close (ncid);
   if (ncid_result) TIO_close (ncid_result);
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
