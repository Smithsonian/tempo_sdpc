#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

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
   fprintf (stderr, "   -S | --s_block i:num       Mirror step blocking specification\n");
   fprintf (stderr, "   -X | --x_block i:num       Xtrack blocking specification\n");
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

static int read_nominal_wavelength (int grp, int xtrack, size_t num_waves, double *y)
{
   TIO_Var_Info_Type info = {0};
   int start[2], count[2];

   if (0 != TIO_inq_var (grp, TEMPO_VAR_WAVELEN_NOMINAL, &info))
     return -1;

   if (1 == info.ndims)
     {
        start[0] = 0;
        count[0] = num_waves;
     }
   else if (2 == info.ndims)
     {
        start[0] = xtrack;
        start[1] = 0;
        count[0] = 1;
        count[1] = num_waves;
     }
   else
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported array shape: %s",
                     __func__, TEMPO_VAR_WAVELEN_NOMINAL);
        return -1;
     }

   if (0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELEN_NOMINAL, start, count, TIO_DOUBLE, y))
     return -1;

   return 0;
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
                               int beg_step, int end_step, size_t step_dimlen,
                               int beg_xtrack, int end_xtrack, size_t xtrack_dimlen,
                               size_t num_spectral_channels,
                               int start_pix, int num_pix, int num_coefs,
                               int fitting_sf_params)
{
   int ncid, varid, start, count;
   int dimids_wavecal_params[3], dimids_sf_params[3];
   size_t params_dimlen = num_coefs;
   size_t num_steps = end_step - beg_step;
   size_t num_xtrack = end_xtrack - beg_xtrack;
   int max_num_steps = step_dimlen;
   int max_num_xtrack = xtrack_dimlen;
   int *steps = NULL;
   int *xtracks = NULL;
   int i;

   if (0 != TIO_create (path, NC_NETCDF4, &ncid))
     return -1;

   if ((0 != TIO_put_att (ncid, NC_GLOBAL, "group_name", TIO_CHAR,
                          strlen(group_name), group_name))
       || (0 != TIO_put_att (ncid, NC_GLOBAL, "mirror_step_dimlen", TIO_INT, 1, &max_num_steps))
       || (0 != TIO_put_att (ncid, NC_GLOBAL, "xtrack_dimlen", TIO_INT, 1, &max_num_xtrack))
      )
     goto close_and_return;

   if ((0 != TIO_def_dim (ncid, TEMPO_DIM_STEP, num_steps, &dimids_wavecal_params[0]))
       || (0 != TIO_def_dim (ncid, TEMPO_DIM_XTRACK, num_xtrack, &dimids_wavecal_params[1]))
       || (0 != TIO_def_dim (ncid, TEMPO_DIM_WAVECAL_PARAM, params_dimlen, &dimids_wavecal_params[2])))
     goto close_and_return;

   /* slitfun parameters [mirror_step, xtrack, spectral_channel] */
   dimids_sf_params[0] = dimids_wavecal_params[0];
   dimids_sf_params[1] = dimids_wavecal_params[1];
   if (0 != TIO_def_dim (ncid, TEMPO_DIM_CHANNEL, num_spectral_channels, &dimids_sf_params[2]))
     goto close_and_return;

   /* wavecal parameters */
   if (0 != TIO_def_var (ncid, TEMPO_VAR_WAVECAL_PARAM, TIO_DOUBLE, 3, dimids_wavecal_params, &varid))
     goto close_and_return;
   if ((0 != TIO_put_att (ncid, varid, "num_coefficients", TIO_INT, 1, &num_coefs))
       ||(0 != TIO_put_att (ncid, varid, "start_spectral_channel", TIO_INT, 1, &start_pix))
       ||(0 != TIO_put_att (ncid, varid, "num_spectral_channels", TIO_INT, 1, &num_pix)))
     goto close_and_return;

   if (0 != TIO_def_var (ncid, "bestnorm", TIO_DOUBLE, 2, dimids_wavecal_params, &varid))
     goto close_and_return;

   if (0 != TIO_def_var (ncid, TEMPO_DIM_STEP, TIO_INT, 1, &dimids_wavecal_params[0], &varid))
     goto close_and_return;
   if (0 != TIO_def_var (ncid, TEMPO_DIM_XTRACK, TIO_INT, 1, &dimids_wavecal_params[1], &varid))
     goto close_and_return;

   /* slit-function parameters */
   if (fitting_sf_params)
     {
        if ((0 != TIO_def_var (ncid, "sf_asym", TIO_FLOAT, 3, dimids_sf_params, &varid))
            ||(0 != TIO_def_var (ncid, "sf_hw1e", TIO_FLOAT, 3, dimids_sf_params, &varid))
            ||(0 != TIO_def_var (ncid, "sf_shape", TIO_FLOAT, 3, dimids_sf_params, &varid)))
          goto close_and_return;
     }

   /* Write step coordinate variable */
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

   /* Write xtrack coordinate variable */
   if (NULL == (xtracks = (int *)MALLOC (num_xtrack * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   for (i = beg_xtrack; i < end_xtrack; i++)
     {
        xtracks[i-beg_xtrack] = i;
     }

   start = 0;
   count = num_xtrack;

   if (0 != TIO_put_var_section (ncid, TEMPO_DIM_XTRACK, &start, &count, TIO_INT, xtracks))
     goto close_and_return;

   FREE(xtracks);
   return ncid;

close_and_return:
   FREE(steps);
   FREE(xtracks);
   (void) TIO_close (ncid);
   return -1;
}

static void fill_array (float a, float *v, int n)
{
   int i;
   for (i = 0; i < n; i++) v[i] = a;
}

static int write_sf_params (int ncid, int beg_step, int step, int beg_xtrack, int xtrack,
                            const Wavecal_Result_Type *wavecal_result)
{
   int num_fit = wavecal_result->num_fit;
   int start_pix = wavecal_result->start_pix;
   float hw1e, shape, asym;
   float *tmp = NULL;
   int start[3], count[3];
   int status = -1;

   if (wavecal_result->sf_params == NULL)
     return 0;

   hw1e = wavecal_result->sf_params[0];
   shape = wavecal_result->sf_params[1];
   asym = wavecal_result->sf_params[2];

   if (NULL == (tmp = (float *)MALLOC (num_fit * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc error", __func__);
        return -1;
     }

   start[0] = step - beg_step;
   start[1] = xtrack - beg_xtrack;
   start[2] = start_pix;
   count[0] = 1;
   count[1] = 1;
   count[2] = num_fit;

   fill_array (hw1e, tmp, num_fit);
   if (0 != TIO_put_var_section (ncid, "sf_hw1e", start, count, TIO_FLOAT, tmp))
     goto return_status;

   fill_array (shape, tmp, num_fit);
   if (0 != TIO_put_var_section (ncid, "sf_shape", start, count, TIO_FLOAT, tmp))
     goto return_status;

   fill_array (asym, tmp, num_fit);
   if (0 != TIO_put_var_section (ncid, "sf_asym", start, count, TIO_FLOAT, tmp))
     goto return_status;

   status = 0;
return_status:
   FREE(tmp);

   return status;
}

static int write_result (int ncid, int beg_step, int step, int beg_xtrack, int xtrack,
                         const double *final_coeff, int num_final_coeff,
                         const Wavecal_Result_Type *wavecal_result)
{
   int start[3], count[3];

   start[0] = step - beg_step;
   start[1] = xtrack - beg_xtrack;
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

        if (0 != write_sf_params (ncid, beg_step, step, beg_xtrack, xtrack, wavecal_result))
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

   if ((0 != TIO_def_var (grp, dimname_wavelen, TIO_DOUBLE, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "model", TIO_DOUBLE, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "spec_scaled", TIO_DOUBLE, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "weight", TIO_DOUBLE, 3, dimids, &varid))
       || (0 != TIO_def_var (grp, "residuals", TIO_DOUBLE, 3, dimids, &varid)))
     {
        return -1;
     }

   return 0;
}

static int write_diagnostics (int grp, int beg_step, int step, int beg_xtrack, int xtrack,
                              const Wavecal_Result_Type *wavecal_result)
{
   int start[3], count[3];

   if (wavecal_result == NULL)
     return 0;

   /* quick return if variables already defined */
   (void) def_diagnostic_vars (grp, wavecal_result);

   start[0] = step - beg_step;
   start[1] = xtrack - beg_xtrack;
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
   int use_s_blocking = 0, this_s_block, num_s_blocks;
   int use_x_blocking = 0, this_x_block, num_x_blocks;
   int is_irradiance = 0;
   int verbose = 0;
   int ncid_result = 0;
   int debug = 0;
   int apply_shift_adjust = 0;
   size_t step_dimlen, xtrack_dimlen, channel_dimlen;
   int num_wave_params, fitting_sf_params, start_pix, num_pix;
   int num_final_coeff, final_start_pix, final_num_pix;
   int fit_status_code, residual, progress, num_expected_results;
   struct timeval tv0 = {0};
   struct timeval tv1 = {0};
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"debug",   no_argument, 0, 'd'},
        {"verbose", no_argument, 0, 'v'},
        {"adjust",  no_argument, 0, 'a'},
        {"config",  required_argument, 0, 'c'},
        {"wavepar", required_argument, 0, 'w'},
        {"xtrack",  required_argument, 0, 'x'},
        {"s_block", required_argument, 0, 'S'},
        {"x_block", required_argument, 0, 'X'},
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
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "Reading %s:%d - %s",
                          config_error_file(&cfg),
                          config_error_line(&cfg), config_error_text(&cfg));
             goto return_status;
          }
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "ahdS:X:m:c:g:x:w:v", long_options, &option_index);
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
               {
                  tell_verror (TELL_INVALID_PARM_ERROR,
                               "Reading %s:%d - %s",
                               config_error_file(&cfg),
                               config_error_line(&cfg), config_error_text(&cfg));
                  goto return_status;
               }
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

           case 'S':
             if (2 != sscanf (optarg, "%d:%d", &this_s_block, &num_s_blocks))
               usage();
             if (this_s_block < 0 || num_s_blocks <= this_s_block)
               {
                  fprintf (stderr, "*** wavecal_driver: invalid step blocking specification\n");
                  goto return_status;
               }
             use_s_blocking++;
             break;
           case 'X':
             if (2 != sscanf (optarg, "%d:%d", &this_x_block, &num_x_blocks))
               usage();
             if (this_s_block < 0 || num_x_blocks <= this_x_block)
               {
                  fprintf (stderr, "*** wavecal_driver: invalid xtrack blocking specification\n");
                  goto return_status;
               }
             use_x_blocking++;
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

   (void) tell_set_log_level (TELL_MSGTYPE_INFO, verbose);
   (void) tell_set_log_level (TELL_MSGTYPE_WARN, verbose);

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

   if (0 != TIO_open (input_file, NC_NOWRITE, &ncid))
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

   /* If too many blocks are requested, ignore the excess */
   if ((use_s_blocking != 0) &&
       (num_s_blocks > (int) step_dimlen))
     {
        if (this_s_block > (int) step_dimlen)
          {
             status = EXIT_SUCCESS;
             goto return_status;
          }
        num_s_blocks = step_dimlen;
     }
   if ((use_x_blocking != 0) &&
       (num_x_blocks > (int) xtrack_dimlen))
     {
        if (this_x_block > (int) xtrack_dimlen)
          {
             status = EXIT_SUCCESS;
             goto return_status;
          }
        num_x_blocks = xtrack_dimlen;
     }

   if (alloc_spectrum (&spec, channel_dimlen))
     goto return_status;

   if (NULL == (y0 = (double *)MALLOC (channel_dimlen * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   if (0 == is_irradiance)
     {
        if (NULL == (geoloc = read_geolocation_vars (grp, spectrum_info.dimlens)))
          goto return_status;

        if (0 != read_sza_max (&cfg, &sza_max))
          goto return_status;
     }

   if (NULL == (wct = wavecal_open (&cfg, group_name, meta, spec.n, is_irradiance)))
     goto return_status;

   if (params_outfile)
     {
        if (NULL == (fp = fopen (params_outfile, "w")))
          {
             fprintf (stderr, "*** Error: opening file for writing: %s\n", params_outfile);
             goto return_status;
          }
     }

   if (use_x_blocking)
     {
        int x_block_size = xtrack_dimlen / num_x_blocks;
        residual = xtrack_dimlen - num_x_blocks * x_block_size;
        beg_xtrack = this_x_block * x_block_size;
        if (residual > 0)
          {
             if (this_x_block < residual)
               {
                  x_block_size += 1;
                  beg_xtrack += this_x_block;
               }
             else beg_xtrack += residual;
          }
        end_xtrack = beg_xtrack + x_block_size;
        if (end_xtrack > (int) xtrack_dimlen)
          end_xtrack = xtrack_dimlen;
     }
   else if (xtrack < 0)
     {
        beg_xtrack = 0;
        end_xtrack = xtrack_dimlen;
     }
   else
     {
        beg_xtrack = xtrack;
        end_xtrack = xtrack+1;
     }

   if (use_s_blocking)
     {
        int s_block_size = step_dimlen / num_s_blocks;
        residual = step_dimlen - num_s_blocks * s_block_size;
        beg_step = this_s_block * s_block_size;
        if (residual > 0)
          {
             if (this_s_block < residual)
               {
                  s_block_size += 1;
                  beg_step += this_s_block;
               }
             else beg_step += residual;
          }
        end_step = beg_step + s_block_size;
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

   if ((fitting_sf_params = wavecal_fitting_sf_params (wct)) < 0)
     goto return_status;

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
                                     beg_step, end_step, step_dimlen,
                                     beg_xtrack, end_xtrack, xtrack_dimlen,
                                     channel_dimlen,
                                     final_start_pix, final_num_pix, num_final_coeff,
                                     fitting_sf_params);
   if (ncid_result <= 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: problem creating result file: %s",
                     __func__, result_file);
        goto return_status;
     }

   /* Try to report progress at regular intervals */
   num_expected_results = (end_step-beg_step) * (end_xtrack-beg_xtrack);
   progress = 0;
   (void) gettimeofday (&tv0, NULL);

   for (step = beg_step; step < end_step; step++)
     {
        for (xtrack = beg_xtrack; xtrack < end_xtrack; xtrack++)
          {
             Wavecal_Result_Type *wrt;

             wrt = NULL;

             if (0 != read_nominal_wavelength (grp, xtrack, channel_dimlen, y0))
               goto return_status;

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

             if (write_result (ncid_result, beg_step, step, beg_xtrack, xtrack, final_coeff, num_final_coeff, wrt))
               goto return_status;

             progress++;
             (void) gettimeofday (&tv1, NULL);
             if (tv1.tv_sec - tv0.tv_sec > 60.0)
               {
                  tell_vinfo (0, "finished %d/%d", progress, num_expected_results);
                  tv0 = tv1;  /* struct copy */
               }

             if (debug)
               {
                  if (write_diagnostics (ncid_result, beg_step, step, beg_xtrack, xtrack, wrt))
                    goto return_status;
               }

             if ((verbose > 1) && (wrt != NULL))
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

   if (0 != tio_meta_write_ncattr (meta, ncid_result))
     goto return_status;

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
