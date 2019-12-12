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
#include "wavecal_adj.h"

#define DIMNAME_WAVELEN "wavelen"

typedef struct
{
   double *spec;
   double *spec_err;
   unsigned int *pixel_quality_flag;
   size_t n;
}
Spectrum_Type;

typedef struct
{
   int *inr_quality_flag;
   double *solar_zenith_angle;
   size_t num_step;
   size_t num_xtrack;
}
Geoloc_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal_driver [options] <input-file> <output-file>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -g | --group NAME          Name of netCDF4 file group containing spectra\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help                Print this usage message\n");
   fprintf (stderr, "   -d | --debug               Write diagnostic information to output file\n");
   fprintf (stderr, "   -a | --adjust              Apply wavelength-shift adjustment\n");
   fprintf (stderr, "   -b | --block i:num         Mirror step blocking specification\n");
   fprintf (stderr, "   -m | --mirror STEP         Mirror step index\n");
   fprintf (stderr, "   -x | --xtrack N            Cross-track pixel index, 0 is northernmost\n");
   fprintf (stderr, "   -w | --wavepar FILE        Output file for wavelength parameters\n");
   fprintf (stderr, "   -c | --config FILE         Path to configuration file\n");
   fprintf (stderr, "   -v | --verbose             Turn on verbose output\n");
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

static double *read_nominal_wavelength (int grp, size_t num_waves)
{
   double *y = NULL;
   int start, count;

   if (NULL == (y = (double *)MALLOC (num_waves * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   start = 0;
   count = num_waves;
   if (0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELEN_NOMINAL, &start, &count, TIO_DOUBLE, y))
     {
        FREE(y);
        return NULL;
     }

   return y;
}

static void free_geoloc_type (Geoloc_Type *g)
{
   if (g == NULL)
     return;
   FREE(g->inr_quality_flag);
   FREE(g->solar_zenith_angle);
   FREE(g);
}

static Geoloc_Type *alloc_geoloc_type (size_t num_step, size_t num_xtrack)
{
   Geoloc_Type *g = NULL;
   size_t len;

   if (NULL == (g = (Geoloc_Type *)MALLOC (sizeof *g)))
     {
        tell_verror(TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   g->num_step = num_step;
   g->num_xtrack = num_xtrack;

   len = num_step * num_xtrack;

   if ((NULL == (g->inr_quality_flag = (int *)MALLOC (len * sizeof(int))))
       || (NULL == (g->solar_zenith_angle = (double *)MALLOC (len * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_geoloc_type (g);
        return NULL;
     }

   memset ((char *)g->inr_quality_flag, 0, len * sizeof(int));
   memset ((char *)g->solar_zenith_angle, 0, len * sizeof(double));

   return g;
}

static Geoloc_Type *read_geolocation_vars (int grp, const size_t *dimlens)
{
   Geoloc_Type *g = NULL;
   int start[2], count[2];

   if (NULL == (g = alloc_geoloc_type (dimlens[0], dimlens[1])))
     return NULL;

   start[0] = 0;
   start[1] = 0;
   count[0] = dimlens[0];
   count[1] = dimlens[1];

   if ((0 != TIO_get_var_section (grp, TEMPO_VAR_INRQF, start, count, TIO_INT, g->inr_quality_flag))
       ||(0 != TIO_get_var_section (grp, TEMPO_VAR_SZ_ANGLE, start, count, TIO_DOUBLE, g->solar_zenith_angle)))
     {
        free_geoloc_type (g);
        return NULL;
     }

   return g;
}

static int read_sza_max (config_t *cfg, double *sza_max)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "wavecal_control")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing wavecal_control in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "sza_max", sza_max))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading sza_max from param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int will_calibrate_radiance (const Geoloc_Type *g, int step, int xtrack, double sza_max)
{
   int *inrqf_step = NULL;
   double *sza_step = NULL;

   if (g == NULL)
     return 1;

   inrqf_step = g->inr_quality_flag + step * g->num_xtrack;
   sza_step = g->solar_zenith_angle + step * g->num_xtrack;

   return ((inrqf_step[xtrack] == 0) && (sza_step[xtrack] < sza_max));
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
        fprintf (stderr, "*** Error opening %s for writing\n", filename);
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
                              const double *wave_params, int num_wave_params,
                              const Wavecal_Result_Type *wavecal_result)
{
   Wavecal_Term_Info_Type info = {0};
   int i, nth;

   if (wavecal_result == NULL)
     {
        fprintf (fp, "# %4d [no result]\n", xtrack);
        return 0;
     }

   fprintf (fp, "%4d %12.4e %4d %4d ", xtrack, wavecal_result->bestnorm,
            wavecal_result->nfev, wavecal_result->opt_status);
   for (i = 0; i < num_wave_params; i++)
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
                               size_t num_xtrack,
                               int start_pix, int num_pix, int num_coefs)
{
   int ncid, varid, param_dimids[3], start, count;
   size_t params_dimlen = num_coefs;
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
       || (0 != TIO_def_dim (ncid, TEMPO_DIM_WAVECAL_PARAM, params_dimlen, &param_dimids[2])))
     goto close_and_return;

   if (0 != TIO_def_var (ncid, TEMPO_VAR_WAVECAL_PARAM, TIO_FLOAT, 3, param_dimids, &varid))
     goto close_and_return;
   if ((0 != TIO_put_att (ncid, varid, "num_coefficients", TIO_INT, 1, &num_coefs))
       ||(0 != TIO_put_att (ncid, varid, "start_spectral_channel", TIO_INT, 1, &start_pix))
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
                         const double *final_coeff, int num_final_coeff,
                         const Wavecal_Result_Type *wavecal_result)
{
   int start[3], count[3];

   start[0] = step - beg_step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = num_final_coeff;

   if (0 != TIO_put_var_section (ncid, TEMPO_VAR_WAVECAL_PARAM, start, count, TIO_DOUBLE,
                                 final_coeff))
     return -1;

   if (wavecal_result)
     {
        if (0 != TIO_put_var_section (ncid, "bestnorm", start, count, TIO_DOUBLE,
                                    &wavecal_result->bestnorm))
          return -1;
     }

   return 0;
}

static int def_diagnostic_vars (int grp, const Wavecal_Result_Type *wavecal_result)
{
   const char *dimname_wavelen = DIMNAME_WAVELEN;
   int varid, dimid_wavelen, dimid_xtrack, dimid_step;
   int dimids[3];
   size_t num_waves = wavecal_result->num_fit;

   if (0 == TIO_inq_dimid (grp, dimname_wavelen, &dimid_wavelen))
     return 0;

   if ((0 != TIO_inq_dimid (grp, TEMPO_DIM_XTRACK, &dimid_xtrack))
       || (0 != TIO_inq_dimid (grp, TEMPO_DIM_STEP, &dimid_step)))
     return -1;

   if (0 != TIO_def_dim (grp, dimname_wavelen, num_waves, &dimid_wavelen))
     return -1;

   dimids[0] = dimid_step;
   dimids[1] = dimid_xtrack;
   dimids[2] = dimid_wavelen;

   if ((0 != TIO_def_var (grp, dimname_wavelen, TIO_FLOAT, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "model", TIO_FLOAT, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "spec_scaled", TIO_FLOAT, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "weight", TIO_FLOAT, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "residuals", TIO_FLOAT, 3, dimids, &varid)))
     {
        return -1;
     }

   return 0;
}

static int write_diagnostics (int grp, int beg_step, int step, int xtrack,
                              const Wavecal_Result_Type *wavecal_result)
{
   int start[3], count[3];

   if (wavecal_result == NULL)
     return 0;

   /* quick return if variables already defined */
   (void) def_diagnostic_vars (grp, wavecal_result);

   start[0] = step - beg_step;
   start[1] = xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = wavecal_result->num_fit;

   if ((0 != TIO_put_var_section (grp, DIMNAME_WAVELEN, start, count, TIO_DOUBLE,
                                  wavecal_result->wave))
       || (0 != TIO_put_var_section (grp, "model", start, count, TIO_DOUBLE,
                                     wavecal_result->model))
       || (0 != TIO_put_var_section (grp, "spec_scaled", start, count, TIO_DOUBLE,
                                     wavecal_result->spec_scaled))
       || (0 != TIO_put_var_section (grp, "weight", start, count, TIO_DOUBLE,
                                     wavecal_result->weight))
       || (0 != TIO_put_var_section (grp, "residuals", start, count, TIO_DOUBLE,
                                     wavecal_result->residuals)))
     {
        return -1;
     }

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
   Wadj_Type *wadj = NULL;
   TIO_Var_Info_Type spectrum_info = {0};
   Wavecal_Config_Type wavecal_config = {0};
   Wavecal_Result_Type wavecal_result = {0};
   Geoloc_Type *geoloc = NULL;
   TIO_Meta_Type *meta = NULL;
   double *y0 = NULL;
   double nan_value = nan("");
   double *wave_params = NULL;
   double *final_coeff = NULL;
   double sza_max;
   int grp, ncid = 0, step = -1, xtrack = -1;
   int beg_xtrack, end_xtrack;
   int beg_step, end_step;
   int use_blocking = 0, this_block, num_blocks;
   int is_irradiance = 0;
   int verbose = 0;
   int ncid_result = 0;
   int debug = 0;
   int apply_shift_adjust = 0;
   size_t step_dimlen, xtrack_dimlen, channel_dimlen;
   int num_wave_params, start_pix, num_pix, grp_meta;
   int num_final_coeff, final_start_pix, final_num_pix;
   int fit_status_code;
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"debug",   no_argument, 0, 'd'},
        {"verbose", no_argument, 0, 'v'},
        {"adjust",  no_argument, 0, 'a'},
        {"config",  required_argument, 0, 'c'},
        {"wavepar", required_argument, 0, 'w'},
        {"xtrack",  required_argument, 0, 'x'},
        {"block",   required_argument, 0, 'b'},
        {"mirror",  required_argument, 0, 'm'},
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
        int c = getopt_long (argc, argv, "ahdb:m:c:g:x:w:v", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d??", c);
             goto return_status;
             break;
           case 'a':
             apply_shift_adjust++;
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

   if (NULL == (meta = tio_meta_open ()))
     goto return_status;

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

   /* The shift-adjustment table is used only
    * for radiance wavelength calibration */
   if (apply_shift_adjust && (0 == is_irradiance))
     {
        if (NULL == (wadj = wadj_open (&cfg, group_name, meta)))
          goto return_status;
     }

   /* expected dimensions are: [mirror_step, xtrack, spectral_channel] */
   step_dimlen = spectrum_info.dimlens[0];
   xtrack_dimlen = spectrum_info.dimlens[1];
   channel_dimlen = spectrum_info.dimlens[2];

   if (alloc_spectrum (&spec, channel_dimlen))
     goto return_status;

   if (NULL == (y0 = read_nominal_wavelength (grp, channel_dimlen)))
     goto return_status;

   if (0 == is_irradiance)
     {
        if (NULL == (geoloc = read_geolocation_vars (grp, spectrum_info.dimlens)))
          goto return_status;

        if (0 != read_sza_max (&cfg, &sza_max))
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

   /* Allocate space to store the fitted Chebyshev series
    * coefficients which represent the wavelength grid
    * within a single spectrum's fit window */
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

   (void) wavecal_query_feature_window (wct, &start_pix, &num_pix);

   /* Decide how many Chebyshev series coefficients (per-spectrum) will be
    * written to the output file, and allocate a working array to hold the
    * coefficients for a single spectrum
    * If we're using a narrow fit window and applying the wavelength
    * shift adjustment, make sure the wavelength shift adjustment
    * lookup table is consistent.
    */
   if (wadj)
     {
        Wadj_Cheb_Series_Type narrow = {0};
        Wadj_Cheb_Series_Type full = {0};
        int num_narrow, num_full;

        if (0 != wadj_num_final_coeff (wadj, &num_final_coeff))
          goto return_status;
        if (NULL == (final_coeff = (double *)MALLOC (num_final_coeff * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto return_status;
          }

        if ((0 != wadj_narrow_band_get_attr (wadj, &narrow))
            ||(0 != wadj_full_band_get_attr (wadj, &full)))
          goto return_status;

        num_narrow = narrow.pix_max - narrow.pix_min + 1;
        num_full = full.pix_max - full.pix_min + 1;

        if ((narrow.pix_min != start_pix)
            || (num_narrow != num_pix))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: shift_adjust_table: narrow window mismatch",
                          __func__);
             goto return_status;
          }

        final_start_pix = full.pix_min;
        final_num_pix = num_full;
     }
   else
     {
        num_final_coeff = num_wave_params;
        final_coeff = wave_params;
        final_start_pix = start_pix;
        final_num_pix = num_pix;
     }

   /* Create a netcdf output file to hold the wavelength grid coefficients */
   ncid_result = create_result_file (result_file, group_name,
                                     beg_step, end_step, step_dimlen, xtrack_dimlen,
                                     final_start_pix, final_num_pix, num_final_coeff);
   if (ncid_result <= 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: problem creating result file: %s",
                     __func__, result_file);
        goto return_status;
     }

   if (0 != TIO_def_grp (ncid_result, "metadata", &grp_meta))
     goto return_status;
   if (0 != tio_meta_write_ncattr (meta, grp_meta))
     goto return_status;

   for (step = beg_step; step < end_step; step++)
     {
        for (xtrack = beg_xtrack; xtrack < end_xtrack; xtrack++)
          {
             Wavecal_Result_Type *wrt;

             wrt = NULL;

             /* Calibrate each spectrum that meets the filter criteria */
             if ((0 != is_irradiance)
                 || (0 != will_calibrate_radiance (geoloc, step, xtrack, sza_max)))
               {
                  if (read_spectrum (grp, step, xtrack, is_irradiance, &spec))
                    goto return_status;

                  fit_status_code = wavecal_fit (wct, xtrack, y0, spec.spec, spec.spec_err,
                                                 spec.pixel_quality_flag, &wavecal_config,
                                                 wave_params, &wavecal_result);
                  if (fit_status_code == WAVECAL_FIT_ERROR)
                    goto return_status;

                  wrt = &wavecal_result;
               }
             else
               {
                  /* Uncalibrated radiance spectra get the default wavelength grid */
                  if (0 != wavecal_get_initial_params (wct, y0, wave_params))
                    goto return_status;
               }

             if (wadj)
               {
                  if (0 != wavecal_adjust (wct, wadj, xtrack, wave_params, final_coeff))
                    goto return_status;
               }

             if (write_result (ncid_result, beg_step, step, xtrack, final_coeff, num_final_coeff, wrt))
               goto return_status;

             if (debug)
               {
                  if (write_diagnostics (ncid_result, beg_step, step, xtrack, wrt))
                    goto return_status;
               }

             if (verbose && wrt)
               {
                  fprintf (stderr, "%4d %4d %12.4e %4d %4d\n", step, xtrack,
                           wrt->bestnorm, wrt->nfev, wrt->opt_status);

                  if (params_outfile)
                    {
                       (void) write_fit_details (fp, xtrack, wct, wave_params, num_wave_params, wrt);
                    }
               }
          }
     }

   status = EXIT_SUCCESS;
return_status:
   FREE(y0);
   FREE(wave_params);
   if (final_coeff != wave_params) FREE(final_coeff);
   free_geoloc_type (geoloc);
   if (ncid) TIO_close (ncid);
   if (ncid_result) TIO_close (ncid_result);
   free_spectrum (&spec);
   config_destroy (&cfg);
   tell_close();
   wavecal_close (wct);
   wadj_close (wadj);
   tio_meta_close (meta);

   if ((fp != NULL) && (params_outfile != NULL))
     {
        (void) fclose (fp);
     }

   return status;
}
