#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "config.h"

#define SENSORCAL_PRIVATE_DATA \
   double *waves; \
   double *rcoeffs; \
   double *plate_trans; \
   double btdf; \
   double diffuser_trend; \
   int num_waves;
#include "sensorcal.h"

typedef struct Cal_Data_Type Cal_Data_Type;

struct Cal_Data_Type
{
   double *waves;            /* wavelengths */
   double *rmetric_conv;     /* radiometric calibration conversion factor */
   double *diffuser_trans;   /* diffuser plate transmission  */
   int num_waves;
   int num_waves_per_ccd;
   double delta_wave;
   double min_wave_uv;
   double min_wave_vis;
};

static void free_cal_data (Cal_Data_Type *data)
{
   if (data == NULL)
     return;
   FREE(data->waves);
   FREE(data);
}

static Cal_Data_Type *new_cal_data (int num_waves)
{
   Cal_Data_Type *data = NULL;
   size_t sizeof_data = 3 * num_waves * sizeof(double);

   if (NULL == (data = (Cal_Data_Type *) MALLOC (sizeof *data)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   data->num_waves = num_waves;

   if (NULL == (data->waves = (double *) MALLOC (sizeof_data)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        FREE(data);
        return NULL;
     }
   memset ((char *)data->waves, 0, sizeof_data);

   data->rmetric_conv = data->waves + num_waves;
   data->diffuser_trans = data->waves + num_waves * 2;

   return data;
}

static Cal_Data_Type *read_cal_file (const char *file)
{
   Cal_Data_Type *data = NULL;
   int ncid, start, count, dimid;
   size_t num_waves;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if (0 != TIO_inq_dim (ncid, "wavelength", &dimid, &num_waves))
     {
        (void) TIO_close (ncid);
        return NULL;
     }

   if (NULL == (data = new_cal_data (num_waves)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        (void) TIO_close (ncid);
        return NULL;
     }

   start = 0;
   count = num_waves;

   if ((0 != TIO_get_var_section (ncid, "wavelength", &start, &count, TIO_DOUBLE,
                                  data->waves))
       ||(0 != TIO_get_var_section (ncid, "conversion_factor", &start, &count, TIO_DOUBLE,
                                  data->rmetric_conv))
       ||(0 != TIO_get_var_section (ncid, "transmission", &start, &count, TIO_DOUBLE,
                                  data->diffuser_trans)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading sensor calibration: %s",
                     __func__, file);
        (void) TIO_close (ncid);
        free_cal_data (data);
        return NULL;
     }

   (void) TIO_close (ncid);

   return data;
}

static void make_wavelength_grid (double *y, int ny, double dy,
                                  double y0_uv, double y0_vis)
{
   int i, num = ny/2;

   for (i = 0; i < num; i++)
     {
        double delta_i = dy * i;
        y[i] = y0_uv + delta_i;
        y[i+num] = y0_vis + delta_i;
     }
}

static int cal_apply_rcoeffs (const Calibration_Type *cal, Image_Type *img)
{
   int y, ny = img->num_rows;
   int x, nx = img->num_cols;

   if (ny != cal->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected wavelength grid size n=%d (expected n=%d)",
                     __func__, ny, cal->num_waves);
        return -1;
     }

   for (y = 0; y < ny; y++)
     {
        Image_Pixel_Type *pixels = img->pixels + y * nx;
        double r_y = cal->rcoeffs[y];
        for (x = 0; x < nx; x++)
          {
             pixels[x] *= r_y;
          }
     }

   return 0;
}

static int cal_apply_btdf (const Calibration_Type *cal,
                           double solar_phi, double solar_theta,
                           Image_Type *img)
{
   int y, ny = img->num_rows;
   int x, nx = img->num_cols;
   double factor = cal->btdf * cal->diffuser_trend;

   /* FIXME: BTDF angle interpolation not implemented yet */
   (void) solar_phi;  (void) solar_theta;

   if (ny != cal->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected wavelength grid size n=%d (expected n=%d)",
                     __func__, ny, cal->num_waves);
        return -1;
     }

   for (y = 0; y < ny; y++)
     {
        Image_Pixel_Type *pixels = img->pixels + y * nx;
        double btdf = factor * cal->plate_trans[y];
        for (x = 0; x < nx; x++)
          {
             pixels[x] /= btdf;
          }
     }

   return 0;
}

static int cal_apply_prnu (const Calibration_Type *cal, Image_Type *img)
{
   (void) cal;
   (void) img;

   /* FIXME: PRNU not implemented yet */

   return 0;
}

static int cal_wavecal (const Calibration_Type *cal, Image_Type *img,
                        Image_Type *img_waves)
{
   int y, ny = img->num_rows;
   int x, nx = img->num_cols;
   double *cal_waves = cal->waves;

   if (ny != cal->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected wavelength grid size n=%d (expected n=%d)",
                     __func__, ny, cal->num_waves);
        return -1;
     }

   /* FIXME: wavelength calibration isn't implemented yet */

   for (y = 0; y < ny; y++)
     {
        Image_Pixel_Type *waves_y = img_waves->pixels + y * nx;
        for (x = 0; x < nx; x++)
          {
             waves_y[x] = cal_waves[y];
          }
     }

   return 0;
}

static void cal_delete (Calibration_Type *cal)
{
   if (cal == NULL)
     return;
   FREE(cal->waves);
   FREE(cal);
}

static Calibration_Type *cal_alloc (int num_waves)
{
   Calibration_Type *cal = NULL;
   size_t sizeof_waves_array = 3 * num_waves * sizeof(double);

   if (NULL == (cal = (Calibration_Type *)MALLOC (sizeof *cal)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)cal, 0, sizeof *cal);

   if (NULL == (cal->waves = (double *) MALLOC (sizeof_waves_array)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        cal_delete (cal);
        return NULL;
     }
   memset ((char *)cal->waves, 0, sizeof_waves_array);

   /* storage for calibration curves vs wavelength */
   cal->num_waves = num_waves;
   cal->rcoeffs = cal->waves + num_waves;
   cal->plate_trans = cal->waves + num_waves * 2;

   /* FIXME: placeholder for BTDF and diffuser trend */
   cal->btdf = 0.1;
   cal->diffuser_trend = 1.0;

   cal->cal_delete = cal_delete;
   cal->cal_apply_rcoeffs = cal_apply_rcoeffs;
   cal->cal_apply_prnu = cal_apply_prnu;
   cal->cal_apply_btdf = cal_apply_btdf;
   cal->cal_wavecal = cal_wavecal;

   return cal;
}

static int spline_interp (const double *x0, const double *y0, size_t n0,
                          const double *x, size_t n, double *y)
{
   gsl_interp_accel *acc;
   gsl_spline *spline;
   int status, interp_status = -1;
   size_t i;

   gsl_set_error_handler_off();

   if ((NULL == (acc = gsl_interp_accel_alloc ()))
       || (NULL == (spline = gsl_spline_alloc (gsl_interp_cspline, n0))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_error;
     }

   if (0 != (status = gsl_spline_init (spline, x0, y0, n0)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: %s",
                     __func__, gsl_strerror(status));
        goto return_error;
     }

   for (i = 0; i < n; i++)
     {
        y[i] = gsl_spline_eval (spline, x[i], acc);
     }

   interp_status = 0;
return_error:
   gsl_spline_free (spline);
   gsl_interp_accel_free (acc);

   return interp_status;
}

static int init_rcoeffs (Calibration_Type *cal,
                         const Cal_Data_Type *data)
{
   double *y0 = NULL;
   double *rmetric_conv = data->rmetric_conv;
   int i, status, n0 = data->num_waves;

   if (NULL == (y0 = (double *) MALLOC (n0 * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < n0; i++)
     {
        /* FIXME:
         * Input rmetric_conv value has units e-/(ph/s/sr/cm^2/nm)
         * which is the conversion factor for 4 spatial pixels.
         * Divide by 4 to obtain the conversion factor per spatial pixel.
         * Prototype also divides by 0.118 sec, which is
         * the nominal value of exposure_time.
         * I think this scaling is needed to account for the specific
         * exposure time and binning used when deriving these coefficients.
         * Presumably we'll eventually get coefficients with documentation
         * that explains these details.
         */
        double y0_i = rmetric_conv[i] / 4 / 0.118;
        y0[i] = (y0_i > 0) ? 1.0/y0_i : 0.0;
     }

   status = spline_interp (data->waves, y0, data->num_waves,
                           cal->waves, cal->num_waves, cal->rcoeffs);
   FREE(y0);

   return status;
}

static int init_plate_trans (Calibration_Type *cal,
                             const Cal_Data_Type *data)
{
   double *y0 = NULL;
   double *diffuser_trans = data->diffuser_trans;
   int i, status, n0 = data->num_waves;

   if (NULL == (y0 = (double *) MALLOC (n0 * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   /* The diffuser has 2 plates, so there are 2 transmission
    * factors, one for each plate: */
   for (i = 0; i < n0; i++)
     {
        double plate_trans = diffuser_trans[i];
        y0[i] = plate_trans * plate_trans;
     }

   status = spline_interp (data->waves, y0, data->num_waves,
                           cal->waves, cal->num_waves, cal->plate_trans);
   FREE(y0);

   return status;
}

static Calibration_Type *cal_init (const Cal_Data_Type *data)
{
   Calibration_Type *cal = NULL;

   if (NULL == (cal = cal_alloc (2 * data->num_waves_per_ccd)))
     return NULL;

   make_wavelength_grid (cal->waves, cal->num_waves, data->delta_wave,
                         data->min_wave_uv, data->min_wave_vis);

   if ((0 != init_rcoeffs (cal, data))
       || (0 != init_plate_trans (cal, data)))
     {
        cal_delete (cal);
        return NULL;
     }

   return cal;
}

Calibration_Type *sensorcal_init (config_t *cfg)
{
   config_setting_t *s, *sub;
   const char *sensorcal_file;
   Calibration_Type *cal = NULL;
   Cal_Data_Type *data = NULL;

   if (NULL == (s = config_lookup (cfg, "ccd_calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_calibration in param file: %s",
                     __func__, config_error_file (cfg));
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "sensorcal_file", &sensorcal_file))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading sensorcal_file",
                     __func__);
        return NULL;
     }

   if (NULL == (data = read_cal_file (sensorcal_file)))
     return NULL;

   if ((NULL == (sub = config_setting_get_member (s, "default_wavelength_grid")))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "delta_wave", &data->delta_wave))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "min_wave_uv", &data->min_wave_uv))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "min_wave_vis", &data->min_wave_vis))
      )
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading default_wavelength_grid parameters",
                     __func__);
        free_cal_data (data);
        return NULL;
     }

   if (NULL == (s = config_lookup (cfg, "ccd_parameters")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_parameters in param file: %s",
                     __func__, config_error_file (cfg));
        free_cal_data (data);
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_parallel_active", &data->num_waves_per_ccd))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading num_parallel_active",
                     __func__);
        free_cal_data (data);
        return NULL;
     }

   cal = cal_init (data);
   free_cal_data (data);

   if (cal == NULL)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: sensorcal init failed", __func__);
     }

   return cal;
}
