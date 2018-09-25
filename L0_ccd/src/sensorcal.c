#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "config.h"
#include "util.h"

typedef struct
{
   const char *name;
   double *waves;
   double *rcoeffs;
   double *plate_trans;
   int num_waves;
   int ybeg;
   int yend;
}
CCD_Cal_Type;

#define CCD_UV(cal)  &(cal)->ccds[TEMPO_BAND_UV]
#define CCD_VIS(cal) &(cal)->ccds[TEMPO_BAND_VIS]

#define SENSORCAL_PRIVATE_DATA \
   CCD_Cal_Type ccds[2]; \
   double btdf; \
   double diffuser_trend; \
   int num_waves_per_ccd;
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

static void make_wavelength_grid (double *y, int ny, double dy, double y0)
{
   int i;

   for (i = 0; i < ny; i++)
     {
        y[i] = y0 + dy * i;
     }
}

static int cal_apply_rcoeffs (const Calibration_Type *cal, Image_Type *img)
{
   const CCD_Cal_Type *uv = CCD_UV(cal);
   const CCD_Cal_Type *vis = CCD_VIS(cal);
   int y, ny = img->num_rows;
   int x, nx = img->num_cols;
   int yoffset;

   if (ny != uv->num_waves + vis->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected wavelength grid size n=%d (expected n=%d)",
                     __func__, ny, uv->num_waves+vis->num_waves);
        return -1;
     }

   for (y = 0; y < uv->num_waves; y++)
     {
        Image_Pixel_Type *pixels = img->pixels + y * nx;
        double r_y = uv->rcoeffs[y];
        for (x = 0; x < nx; x++)
          {
             if (pixels[x] != IMAGE_PIXEL_FILL_VALUE)
               pixels[x] *= r_y;
          }
     }

   yoffset = uv->num_waves;

   for (y = 0; y < vis->num_waves; y++)
     {
        Image_Pixel_Type *pixels = img->pixels + (y + yoffset) * nx;
        double r_y = vis->rcoeffs[y];
        for (x = 0; x < nx; x++)
          {
             if (pixels[x] != IMAGE_PIXEL_FILL_VALUE)
               pixels[x] *= r_y;
          }
     }

   return 0;
}

static int cal_apply_btdf (const Calibration_Type *cal,
                           double solar_phi, double solar_theta,
                           Image_Type *img)
{
   const CCD_Cal_Type *uv = CCD_UV(cal);
   const CCD_Cal_Type *vis = CCD_VIS(cal);
   int y, ny = img->num_rows;
   int x, nx = img->num_cols;
   int yoffset;
   double factor = cal->btdf * cal->diffuser_trend;

   /* FIXME: BTDF angle interpolation not implemented yet */
   (void) solar_phi;  (void) solar_theta;

   if (ny != uv->num_waves + vis->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected wavelength grid size n=%d (expected n=%d)",
                     __func__, ny, uv->num_waves + vis->num_waves);
        return -1;
     }

   for (y = 0; y < uv->num_waves; y++)
     {
        Image_Pixel_Type *pixels = img->pixels + y * nx;
        double btdf = factor * uv->plate_trans[y];
        for (x = 0; x < nx; x++)
          {
             if (pixels[x] != IMAGE_PIXEL_FILL_VALUE)
               pixels[x] /= btdf;
          }
     }

   yoffset = uv->num_waves;

   for (y = 0; y < vis->num_waves; y++)
     {
        Image_Pixel_Type *pixels = img->pixels + (y + yoffset) * nx;
        double btdf = factor * vis->plate_trans[y];
        for (x = 0; x < nx; x++)
          {
             if (pixels[x] != IMAGE_PIXEL_FILL_VALUE)
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

static const CCD_Cal_Type *ccd_cal (const Calibration_Type *cal,
                                    int band_index)
{
   switch (band_index)
     {
      case TEMPO_BAND_UV:
        return CCD_UV(cal);

      case TEMPO_BAND_VIS:
        return CCD_VIS(cal);

      default:
        break;
     }

   tell_verror (TELL_INVALID_PARM_ERROR,
                "%s: invalid band index = %d", __func__, band_index);
   return NULL;
}

static int cal_wavecal (const Calibration_Type *cal, int band_index,
                        int nx, const double *pspec, const double *pspec_err,
                        double *pwaves)
{
   const CCD_Cal_Type *ccd = NULL;

   (void) nx; (void) pspec; (void) pspec_err;

   /* Wavelength calibration is performed elsewhere.
    * Here, we just copy the nominal wavelength grid.
    */

   if (NULL == (ccd = ccd_cal (cal, band_index)))
     return -1;

   if (cal->num_waves_per_ccd != ccd->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected wavelength grid size n=%d (expected n=%d)",
                     __func__, cal->num_waves_per_ccd, ccd->num_waves);
        return -1;
     }

   memcpy ((char *)pwaves, (char *)ccd->waves,
           cal->num_waves_per_ccd * sizeof(double));

   return 0;
}

static void free_ccd_cal (CCD_Cal_Type *ccd)
{
   if (ccd == NULL)
     return;
   FREE(ccd->waves);
   FREE(ccd->rcoeffs);
   FREE(ccd->plate_trans);
}

static int alloc_ccd_cal (CCD_Cal_Type *ccd, int num_waves)
{
   if (ccd == NULL)
     return -1;

   if ((NULL == (ccd->waves = alloc_doubles (num_waves)))
       || (NULL == (ccd->rcoeffs = alloc_doubles (num_waves)))
       || (NULL == (ccd->plate_trans = alloc_doubles (num_waves))))
     return -1;

   ccd->num_waves = num_waves;

   return 0;
}

static void cal_delete (Calibration_Type *cal)
{
   if (cal == NULL)
     return;
   free_ccd_cal (CCD_UV(cal));
   free_ccd_cal (CCD_VIS(cal));
   FREE(cal);
}

static Calibration_Type *cal_alloc (int num_waves_per_ccd)
{
   Calibration_Type *cal = NULL;
   CCD_Cal_Type *uv, *vis;

   if (NULL == (cal = (Calibration_Type *)MALLOC (sizeof *cal)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)cal, 0, sizeof *cal);

   cal->num_waves_per_ccd = num_waves_per_ccd;

   if ((0 != alloc_ccd_cal (CCD_UV(cal), num_waves_per_ccd))
       || (0 != alloc_ccd_cal (CCD_VIS(cal), num_waves_per_ccd)))
     {
        cal_delete (cal);
        return NULL;
     }

   uv = CCD_UV(cal);
   uv->name = TEMPO_BAND_NAME_UV;
   uv->ybeg = 0;
   uv->yend = uv->num_waves;

   vis = CCD_VIS(cal);
   vis->name = TEMPO_BAND_NAME_VIS;
   vis->ybeg = uv->num_waves;
   vis->yend = vis->ybeg + vis->num_waves;

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

static int init_rcoeffs (const Cal_Data_Type *data,
                         Calibration_Type *cal)
{
   CCD_Cal_Type *uv = CCD_UV(cal);
   CCD_Cal_Type *vis = CCD_VIS(cal);
   double *rmetric_conv = data->rmetric_conv;
   int i, n0 = data->num_waves;
   double *y0 = NULL;

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

   if ((0 != spline_interp (data->waves, y0, data->num_waves,
                           uv->waves, uv->num_waves, uv->rcoeffs))
       || (0 != spline_interp (data->waves, y0, data->num_waves,
                               vis->waves, vis->num_waves, vis->rcoeffs)))
     {
        FREE(y0);
        return -1;
     }

   FREE(y0);
   return 0;
}

static int init_plate_trans (const Cal_Data_Type *data,
                             Calibration_Type *cal)
{
   CCD_Cal_Type *uv = CCD_UV(cal);
   CCD_Cal_Type *vis = CCD_VIS(cal);
   double *diffuser_trans = data->diffuser_trans;
   int i, n0 = data->num_waves;
   double *y0 = NULL;

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

   if ((0 !=spline_interp (data->waves, y0, data->num_waves,
                           uv->waves, uv->num_waves, uv->plate_trans))
       ||(0 !=spline_interp (data->waves, y0, data->num_waves,
                             vis->waves, vis->num_waves, vis->plate_trans)))
     {
        FREE(y0);
        return -1;
     }

   FREE(y0);
   return 0;
}

static Calibration_Type *cal_init (const Cal_Data_Type *data)
{
   Calibration_Type *cal = NULL;
   CCD_Cal_Type *uv, *vis;

   if (NULL == (cal = cal_alloc (data->num_waves_per_ccd)))
     return NULL;

   uv = CCD_UV(cal);
   vis = CCD_VIS(cal);

   make_wavelength_grid (uv->waves, uv->num_waves,
                         data->delta_wave, data->min_wave_uv);
   make_wavelength_grid (vis->waves, vis->num_waves,
                         data->delta_wave, data->min_wave_vis);

   if ((0 != init_rcoeffs (data, cal))
       || (0 != init_plate_trans (data, cal)))
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
   char *path = NULL;

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

   if (NULL == (path = expand_path (sensorcal_file)))
     return NULL;
   data = read_cal_file (path);
   FREE(path);
   if (NULL == data)
     return NULL;

   if ((NULL == (sub = config_setting_get_member (s, "nominal_wavelength_grid")))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "delta_wave", &data->delta_wave))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "min_wave_uv", &data->min_wave_uv))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "min_wave_vis", &data->min_wave_vis))
      )
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading nominal_wavelength_grid parameters",
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

   return cal;
}

void sdt_free (Spectral_Data_Type *sdt)
{
   if (sdt == NULL)
     return;
   FREE(sdt->wave);
   FREE(sdt->img);
   FREE(sdt->pqf);
   FREE(sdt);
}

static Spectral_Data_Type *sdt_alloc (int num_xtrack, int num_channels)
{
   Spectral_Data_Type *sdt = NULL;
   size_t img_size = num_xtrack * num_channels;

   if (NULL == (sdt = (Spectral_Data_Type *)MALLOC (sizeof *sdt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sdt, 0, sizeof *sdt);

   sdt->num_xtrack = num_xtrack;
   sdt->num_channels = num_channels;

   if ((NULL == (sdt->img = (double *)MALLOC (2 * img_size * sizeof(double))))
       || (NULL == (sdt->pqf = (Image_Pqf_Bitmap_Type *)MALLOC (img_size * sizeof(Image_Pqf_Bitmap_Type))))
       || (NULL == (sdt->wave = (double *)MALLOC (num_channels * sizeof(double)))))
     {
        sdt_free (sdt);
        return NULL;
     }
   sdt->img_err = sdt->img + img_size;

   return sdt;
}

static void copy_image_pixels_to_wavelength_order
(const Image_Type *img, int ybeg, int yend, double *outbuf)
{
   Image_Pixel_Type *pixels_img = img->pixels;
   int x, nx = img->num_cols;
   int y, ny = yend - ybeg;

   /* Copy wavelength range [ybeg,yend) from img -> outbuf.
    * In img, x varies fastest. In outbuf, y varies fastest. */

   for (x = 0; x < nx; x++)
     {
        double *pixels_out = outbuf + x * ny;
        for (y = ybeg; y < yend; y++)
          {
             pixels_out[y-ybeg] = pixels_img[x + y * nx];
          }
     }
}

static void copy_image_pqf_to_wavelength_order
(const Image_Type *img, int ybeg, int yend, Image_Pqf_Bitmap_Type *outbuf)
{
   const Image_Pqf_Bitmap_Type *img_pqf = img->pixel_quality_flags;
   int x, nx = img->num_cols;
   int y, ny = yend - ybeg;

   /* Copy wavelength range [ybeg,yend) from img_pqf -> outbuf.
    * In img, x varies fastest. In outbuf, y varies fastest. */

   for (x = 0; x < nx; x++)
     {
        Image_Pqf_Bitmap_Type *pqf_out = outbuf + x * ny;
        for (y = ybeg; y < yend; y++)
          {
             pqf_out[y-ybeg] = img_pqf[x + y * nx];
          }
     }
}

Spectral_Data_Type *
sdt_extract_band (const Calibration_Type *cal, int band_id,
                  const Image_Type *img,
                  const Image_Type *img_err)
{
   Spectral_Data_Type *sdt = NULL;
   const CCD_Cal_Type *ccd;

   if (NULL == (ccd = ccd_cal (cal, band_id)))
     return NULL;

   if (NULL == (sdt = sdt_alloc (img->num_cols, ccd->num_waves)))
     return NULL;

   sdt->name = ccd->name;
   copy_image_pixels_to_wavelength_order (img, ccd->ybeg, ccd->yend, sdt->img);
   copy_image_pixels_to_wavelength_order (img_err, ccd->ybeg, ccd->yend, sdt->img_err);
   copy_image_pqf_to_wavelength_order (img, ccd->ybeg, ccd->yend, sdt->pqf);

   return sdt;
}
