#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_float.h>

#include "config.h"
#include "util.h"

#define SLIT_AOV_MAX_DEG (2.3)

#define SHADOW_TOP (1<<0)
#define SHADOW_BOT (1<<1)

typedef struct
{
   int is_reference_diffuser;

   int num_waves;
   float *waves;
   /* [nm] wavelength */

   int num_aov;
   float *aov;
   /* [deg] "angle of view" */

   int num_slope_aoi;
   float *slope_aoi;
   /* polynomial coefficients for interpolating the dependence on angle of incidence
    * (polar angle, theta) */

   float aoi_nom;
   /* [deg] nominal angle of incidence (polar angle, theta) */

   float *btdfe_lut;
   /* [dimensionless] (num_waves,num_aov) Effective BTDF, averaged over pixel FOV */

   float *trend;
   /* [dimensionless] (row,col) BTDF trend multiplicative correction factor
    * row = xtrack index, col = spectral index (fastest varying)
    */
}
BTDF_Type;

typedef struct
{
   int num_waves;
   int num_kernels;

   float *bb_kernels_transpose;
   /**< Transpose of Broad band kernels for stray light [num_kernels,num_waves] */

   float *bb_stray_light;
   /**< Broad band stray light for each kernel [num_waves,num_kernels] */

   float *bb_source_mean;
   /**< Broad band source mean for each kernel [num_kernels] */
}
BB_Kernel_Type;

typedef struct
{
   size_t num_waves;
   float *Ainv;
   /**< Ainv matrix, square, [num_waves,num_waves] */

   double scale_factor_vis_to_uv;
   int use_shadows;
}
PSF_Matrix_Type;

typedef struct
{
   float *spec;
   int *count;
   int num_rows;
   int first_col;
   int last_col;
}
Shadow_Type;

#define SENSORCAL_PRIVATE_DATA \
   BB_Kernel_Type *sl_bbk; \
   PSF_Matrix_Type *sl_psf; \
   BTDF_Type *diffuser_wrk; \
   BTDF_Type *diffuser_ref; \
   float *wavelength_grid; \
   float *radcal_coeffs; \
   int num_waves; \
   int num_xpos; \
   double btdf; \
   double diffuser_trend; \
   unsigned int straylight_shadow_method;
#include "sensorcal.h"

static int read_wavelength_grid (Calibration_Type *cal, const char *file)
{
   size_t num_waves, num_xpos, len;
   int start[2], count[2];
   int ncid, dimid, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if ((0 != TIO_inq_dim (ncid, "wave", &dimid, &num_waves))
       || (0 != TIO_inq_dim (ncid, "Xpos", &dimid, &num_xpos)))
     goto close_and_return;

   cal->num_waves = num_waves;
   cal->num_xpos = num_xpos;

   len = num_waves * num_xpos;

   if (NULL == (cal->wavelength_grid = (float *)MALLOC (len * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_xpos;

   if (0 != TIO_get_var_section (ncid, "wavelength_grid", start, count, TIO_FLOAT,
                                  cal->wavelength_grid))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading nominal wavelength grid: %s",
                     __func__, file);
        goto close_and_return;
     }

   status = 0;
close_and_return:
   (void) TIO_close (ncid);
   if (status)
     {
        FREE(cal->wavelength_grid);
        cal->wavelength_grid = NULL;
     }

   return 0;
}

static int apply_radcal_trend_correction (Calibration_Type *cal, const char *file)
{
   const char var_name[] = "radcal_coeffs_trend_correction";
   size_t num_waves, num_xpos, len, i;
   float *trend_corr = NULL;
   int start[2], count[2];
   int ncid, dimid, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s (%s)", file, var_name);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if ((0 != TIO_inq_dim (ncid, "wave", &dimid, &num_waves))
       || (0 != TIO_inq_dim (ncid, "Xpos", &dimid, &num_xpos)))
     goto close_and_return;

   if (((int) num_waves != cal->num_waves)
       || ((int) num_xpos != cal->num_xpos))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mismatched array dimensions: %s (expected wave=%d, Xpos=%d)",
                     __func__, file, cal->num_waves, cal->num_xpos);
        goto close_and_return;
     }

   len = num_waves * num_xpos;

   if (NULL == (trend_corr = (float *)MALLOC (len * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_xpos;

   if (0 != TIO_get_var_section (ncid, var_name, start, count, TIO_FLOAT, trend_corr))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %s: %s", __func__, var_name, file);
        goto close_and_return;
     }

   for (i = 0; i < len; i++)
     {
        cal->radcal_coeffs[i] *= trend_corr[i];
     }

   status = 0;
close_and_return:
   (void) TIO_close (ncid);
   FREE(trend_corr);

   return status;
}

static int read_radcal_coeffs (Calibration_Type *cal, const char *file,
                              const char *trend_file)
{
   size_t num_waves, num_xpos, len, i;
   int start[2], count[2];
   int ncid, dimid, status = -1;
   int count_invalid;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if ((0 != TIO_inq_dim (ncid, "wave", &dimid, &num_waves))
       || (0 != TIO_inq_dim (ncid, "Xpos", &dimid, &num_xpos)))
     goto close_and_return;

   cal->num_waves = num_waves;
   cal->num_xpos = num_xpos;

   len = num_waves * num_xpos;

   if (NULL == (cal->radcal_coeffs = (float *)MALLOC (len * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_xpos;

   if (0 != TIO_get_var_section (ncid, "radcal_coeffs", start, count, TIO_FLOAT,
                                  cal->radcal_coeffs))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading sensor calibration: %s",
                     __func__, file);
        goto close_and_return;
     }

   count_invalid = 0;
   for (i = 0; i < len; i++)
     {
        if (cal->radcal_coeffs[i] < 0.0)
          {
             count_invalid++;
             cal->radcal_coeffs[i] = 0.0;
          }
     }
   if (count_invalid)
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0, "read/zeroed %d invalid radcal_coeffs from %s",
                   count_invalid, file);
     }

   if (trend_file)
     {
        if (0 != apply_radcal_trend_correction (cal, trend_file))
          goto close_and_return;
     }

   status = 0;
close_and_return:
   (void) TIO_close (ncid);
   if (status)
     {
        FREE(cal->radcal_coeffs);
        cal->radcal_coeffs = NULL;
     }

   return status;
}

static int cal_apply_radcal_coeffs (const Calibration_Type *cal, Image_Type *img)
{
   Image_Pixel_Type *pixels = img->pixels;
   float *radcal_coeffs = cal->radcal_coeffs;
   size_t i, n;

   if ((cal->num_waves != img->num_rows)
       || (cal->num_xpos != img->num_cols))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: array size mismatch: img[%d,%d] (expected [%d,%d])", __func__,
                     img->num_rows, img->num_cols,
                     cal->num_waves, cal->num_xpos);
        return -1;
     }

   n = img->num_rows * img->num_cols;

   for (i = 0; i < n; i++)
     {
        if (pixels[i] != IMAGE_PIXEL_FILL_VALUE)
          {
             pixels[i] *= radcal_coeffs[i];
          }
     }

   return 0;
}

static void btdf_free (BTDF_Type *btdf)
{
   if (btdf == NULL)
     return;
   FREE(btdf->waves);
   FREE(btdf->aov);
   FREE(btdf->slope_aoi);
   FREE(btdf->btdfe_lut);
   FREE(btdf->trend);
   FREE(btdf);
}

static BTDF_Type *btdf_alloc (size_t num_waves, size_t num_aov, size_t num_slope_aoi)
{
   BTDF_Type *btdf = NULL;

   if (NULL == (btdf = (BTDF_Type *)MALLOC (sizeof *btdf)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)btdf, 0, sizeof (*btdf));

   if ((NULL == (btdf->waves = (float *)MALLOC (num_waves * sizeof(float))))
       || (NULL == (btdf->aov = (float *)MALLOC (num_aov * sizeof(float))))
       || (NULL == (btdf->slope_aoi = (float *)MALLOC (num_slope_aoi * sizeof(float))))
       || (NULL == (btdf->btdfe_lut = (float *)MALLOC (num_aov * num_waves * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        btdf_free (btdf);
        return NULL;
     }

   btdf->num_aov = num_aov;
   btdf->num_waves = num_waves;
   btdf->num_slope_aoi = num_slope_aoi;

   return btdf;
}

static int init_btdf_trend_correction (BTDF_Type *btdf, int num_xpos, int num_waves,
                                       const char *trend_file)
{
   const char *var_name = NULL;
   size_t num_waves_file, num_xpos_file, len, i;
   int start[2], count[2];
   int ncid, dimid, status = -1;

   len = num_waves * num_xpos;

   if (NULL == (btdf->trend = (float *)MALLOC (len * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   if (trend_file == NULL)
     {
        for (i = 0; i < len; i++)
          {
             btdf->trend[i] = 1.0;
          }
        return 0;
     }

   if (btdf->is_reference_diffuser)
     {
        var_name = "BTDF_ref_trend_correction";
     }
   else
     {
        var_name = "BTDF_work_trend_correction";
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s (%s)", trend_file, var_name);

   if (0 != TIO_open (trend_file, NC_NOWRITE, &ncid))
     return -1;

   if ((0 != TIO_inq_dim (ncid, "wave", &dimid, &num_waves_file))
       || (0 != TIO_inq_dim (ncid, "Xpos", &dimid, &num_xpos_file)))
     goto close_and_return;

   if ((num_waves != (int) num_waves_file)
       || (num_xpos != (int) num_xpos_file))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mismatched array dimensions: %s (expected wave=%d, Xpos=%d)",
                     __func__, trend_file, num_waves, num_xpos);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_xpos;

   if (0 != TIO_get_var_section (ncid, var_name, start, count, TIO_FLOAT, btdf->trend))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %s: %s", __func__, var_name, trend_file);
        goto close_and_return;
     }

   status = 0;
close_and_return:
   (void) TIO_close (ncid);

   return status;
}

static BTDF_Type *read_btdf_parameters (int is_reference_diffuser, const char *file)
{
   BTDF_Type *btdf = NULL;
   const char *slope_aoi_name;
   const char *lut_name;
   size_t num_waves, num_aov, num_slope_aoi;
   int start[2], count[2];
   int ncid, dimid, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if ((0 != TIO_inq_dim (ncid, "n_BTDF_w", &dimid, &num_waves))
       || (0 != TIO_inq_dim (ncid, "n_BTDF_aov", &dimid, &num_aov))
       || (0 != TIO_inq_dim (ncid, "n_BTDF_slope_aoi", &dimid, &num_slope_aoi)))
     goto close_and_return;

   if (NULL == (btdf = btdf_alloc (num_waves, num_aov, num_slope_aoi)))
     goto close_and_return;

   btdf->is_reference_diffuser = is_reference_diffuser;

   if (is_reference_diffuser)
     {
        slope_aoi_name = "BTDF_ref_slope_aoi";
        lut_name = "BTDFe_ref_lut";
     }
   else
     {
        slope_aoi_name = "BTDF_work_slope_aoi";
        lut_name = "BTDFe_work_lut";
     }

   start[0] = 0;
   count[0] = num_waves;

   if (0 != TIO_get_var_section (ncid, "BTDF_w", start, count, TIO_FLOAT, btdf->waves))
     goto close_and_return;

   start[0] = 0;
   count[0] = num_aov;

   if (0 != TIO_get_var_section (ncid, "BTDF_aov", start, count, TIO_FLOAT, btdf->aov))
     goto close_and_return;

   start[0] = 0;
   count[0] = 1;

   if (0 != TIO_get_var_section (ncid, "BTDF_aoi_nom", start, count, TIO_FLOAT, &btdf->aoi_nom))
     goto close_and_return;

   start[0] = 0;
   count[0] = num_slope_aoi;

   if (0 != TIO_get_var_section (ncid, slope_aoi_name, start, count, TIO_FLOAT, btdf->slope_aoi))
     goto close_and_return;

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_aov;

   if (0 != TIO_get_var_section (ncid, lut_name, start, count, TIO_FLOAT, btdf->btdfe_lut))
     goto close_and_return;

   status = 0;
close_and_return:
   if (ncid)
     {
        (void) TIO_close (ncid);
     }

   if (status)
     {
        btdf_free (btdf);
        btdf = NULL;
     }

   return btdf;
}

static int read_btdf (Calibration_Type *cal, const char *path, const char *trend_file)
{
   if ((NULL == (cal->diffuser_wrk = read_btdf_parameters (0, path)))
       || (0 != init_btdf_trend_correction (cal->diffuser_wrk, cal->num_xpos, cal->num_waves, trend_file)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading working BTDF parameters from %s",
                     __func__, path);
        return -1;
     }

   if ((NULL == (cal->diffuser_ref = read_btdf_parameters (1, path)))
       || (0 != init_btdf_trend_correction (cal->diffuser_ref, cal->num_xpos, cal->num_waves, trend_file)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading reference BTDF parameters from %s",
                     __func__, path);
        return -1;
     }

   return 0;
}

#define BTDF_SOLAR_THETA_NOMINAL     30.0  /* deg */
#define BTDF_SOLAR_THETA_WARN_DELTA   5.0  /* deg */

static int cal_apply_btdf (const Calibration_Type *cal,
                           int is_reference_diffuser,
                           double solar_phi_deg, double solar_theta_deg,
                           Image_Type *img)
{
   const BTDF_Type *bt;
   float aov_min, aov_step, hs;
   float wave_min, wave_step;
   int p, s;

   if (enable_state_query_bool (ENABLE_BTDF) < 1)
     return 0;

   tell_vlog (TELL_MSGTYPE_INFO, 2, "btdf solar angles: theta = %7.3f phi = %7.3f deg",
              solar_theta_deg, solar_phi_deg);

   if (fabs (solar_theta_deg - BTDF_SOLAR_THETA_NOMINAL) > BTDF_SOLAR_THETA_WARN_DELTA)
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "diffuser solar incidence angle = %0.3f deg (>%g deg away from nominal value = %0.3f deg)",
                   solar_theta_deg, BTDF_SOLAR_THETA_WARN_DELTA, BTDF_SOLAR_THETA_NOMINAL);
     }

   /* Total BTDF has no significant dependence on azimuth.
    * Polarization correction does have azimuthal dependence.
    */
   (void) solar_phi_deg;

   if (img->num_cols != cal->num_xpos)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: image size mismatch: img->num_cols=%d (expected %d)",
                     __func__, img->num_cols, cal->num_xpos);
        return -1;
     }
   if (img->num_rows != cal->num_waves)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: image size mismatch: img->num_rows=%d (expected %d)",
                     __func__, img->num_rows, cal->num_waves);
        return -1;
     }

   bt = is_reference_diffuser ? cal->diffuser_ref : cal->diffuser_wrk;

   /* This code assumes that bt->waves and bt->aov are
    * monotonic increasing, with fixed-width grid spacing.
    */

   wave_min = bt->waves[0];
   wave_step = bt->waves[1] - wave_min;

   aov_min = bt->aov[0];
   aov_step = bt->aov[1] - aov_min;

   hs = (img->num_cols - 1)/2.0;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *img_pixels = img->pixels + p * img->num_cols;
        float *waves = cal->wavelength_grid + p * cal->num_xpos;
        float *bt_trend = bt->trend + p * cal->num_xpos;
        for (s = 0; s < img->num_cols; s++)
          {
             float *lut_w0, *lut_w1;
             float aov_s, wave_s, fa0, fa1, fw0, fw1, b0, b1;
             float aoi_correction, btdfe_s;
             int iw, ia;

             if (img_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;

             aov_s = SLIT_AOV_MAX_DEG * (1.0 - s/hs);
             wave_s = waves[s];

             /* angle of incidendence correction */
             aoi_correction = ((bt->slope_aoi[0] + wave_s * bt->slope_aoi[1])
                               * (solar_theta_deg - bt->aoi_nom) / 100.0);

             /* Bilinear interpolation of effective BTDF, with no extrapolation */

             iw = (waves[s] - wave_min) / wave_step;
             ia = (aov_s - aov_min) / aov_step;

             if (iw < 0)
               {
                  iw = 0;
                  wave_s = bt->waves[iw];
               }
             else if (iw > bt->num_waves-2)
               {
                  iw = bt->num_waves-2;
                  wave_s = bt->waves[iw];
               }

             if (ia < 0)
               {
                  ia = 0;
                  aov_s = bt->aov[ia];
               }
             else if (ia > bt->num_aov-2)
               {
                  ia = bt->num_aov-2;
                  aov_s = bt->aov[ia];
               }

             /* LUT rows bracketing the wavelength */
             lut_w0 = bt->btdfe_lut + bt->num_aov * iw;
             lut_w1 = bt->btdfe_lut + bt->num_aov * (iw + 1);

             /* In each row, interpolate in aov: */
             fa0 = (bt->aov[ia+1] - aov_s) / aov_step;
             fa1 = (aov_s   - bt->aov[ia]) / aov_step;

             b0 = fa0 * lut_w0[ia] + fa1 * lut_w0[ia+1];
             b1 = fa0 * lut_w1[ia] + fa1 * lut_w1[ia+1];

             /* interpolate in wavelength: */
             fw0 = (bt->waves[iw+1] - wave_s) / wave_step;
             fw1 = (wave_s   - bt->waves[iw]) / wave_step;

             btdfe_s = (fw0 * b0 + fw1 * b1) * (1.0 + aoi_correction);

             img_pixels[s] /= btdfe_s * bt_trend[s];
          }
     }

   return 0;
}

static void free_bb_kernel_type (BB_Kernel_Type *sl)
{
   if (sl == NULL)
     return;
   FREE(sl->bb_kernels_transpose);
   FREE(sl->bb_stray_light);
   FREE(sl->bb_source_mean);
   FREE(sl);
}

static BB_Kernel_Type *alloc_bb_kernel_type (size_t num_waves, size_t num_kernels)
{
   BB_Kernel_Type *sl = NULL;

   if (NULL == (sl = (BB_Kernel_Type *)MALLOC (sizeof *sl)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sl, 0, sizeof *sl);

   if ((NULL == (sl->bb_kernels_transpose = (float *) MALLOC (num_waves * num_kernels * sizeof(float))))
       || (NULL == (sl->bb_stray_light = (float *) MALLOC (num_waves * num_kernels * sizeof(float))))
       || (NULL == (sl->bb_source_mean = (float *) MALLOC (num_kernels * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_bb_kernel_type (sl);
        return NULL;
     }

   sl->num_waves = num_waves;
   sl->num_kernels = num_kernels;

   return sl;
}

static int read_bb_kernels (Calibration_Type *cal, config_t *cfg)
{
   BB_Kernel_Type *sl = NULL;
   float *bb_kernels = NULL;
   const char *path_str;
   char *path = NULL;
   size_t i, k, num_waves, num_kernels;
   int start[2], count[2];
   int ncid, dimid, status = -1;

   if (CONFIG_TRUE != config_lookup_string (cfg, "calibration.straylight.bb_kernel", &path_str))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading ", __func__);
        return -1;
     }

   if (NULL == (path = expand_string (path_str)))
     return -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", path);

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     goto close_and_return;

   if ((0 != TIO_inq_dim (ncid, "wave", &dimid, &num_waves))
       || (0 != TIO_inq_dim (ncid, "n_BBSL_kernel", &dimid, &num_kernels)))
     goto close_and_return;

   if (NULL == (sl = alloc_bb_kernel_type (num_waves, num_kernels)))
     goto close_and_return;

   if (NULL == (bb_kernels = (float *)MALLOC (num_waves * num_kernels * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_kernels;

   if ((0 != TIO_get_var_section (ncid, "BB_kernels", start, count, TIO_FLOAT, bb_kernels))
       || (0 != TIO_get_var_section (ncid, "BB_SLs", start, count, TIO_FLOAT, sl->bb_stray_light)))
     goto close_and_return;

   count[0] = num_kernels;
   if (0 != TIO_get_var_section (ncid, "BB_SourceMean", start, count, TIO_FLOAT, sl->bb_source_mean))
     goto close_and_return;

   /* Apparently the kernels need to be normalized */
   for (i = 0; i < num_waves; i++)
     {
        float *kern = bb_kernels + i * num_kernels;
        for (k = 0; k < num_kernels; k++)
          {
             kern[k] /= num_waves * sl->bb_source_mean[k];
          }
     }

   /* Store the kernels transposed
    * [waves,b] -> [b,waves]
    */
   for (i = 0; i < num_waves; i++)
     {
        float *kern = bb_kernels + i * num_kernels;
        float *kern_Ti = sl->bb_kernels_transpose + i;
        for (k = 0; k < num_kernels; k++)
          {
             kern_Ti[k * num_waves] = kern[k];
          }
     }

   cal->sl_bbk = sl;
   status = 0;

close_and_return:
   (void) TIO_close (ncid);
   if (status)
     {
        free_bb_kernel_type (sl);
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading %s", __func__, path);
     }
   FREE(bb_kernels);
   FREE(path);

   return status;
}

static void free_psf_matrix_type (PSF_Matrix_Type *psf)
{
   if (psf == NULL)
     return;
   FREE(psf->Ainv);
   FREE(psf);
}

static PSF_Matrix_Type *alloc_psf_matrix_type (size_t num_waves)
{
   PSF_Matrix_Type *psf = NULL;

   if (NULL == (psf = (PSF_Matrix_Type *)MALLOC (sizeof *psf)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)psf, 0, sizeof (*psf));

   if (NULL == (psf->Ainv = (float *)MALLOC (num_waves * num_waves * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_psf_matrix_type (psf);
        return NULL;
     }

   psf->num_waves = num_waves;

   return psf;
}

static int read_sl_psf_matrix (Calibration_Type *cal, config_t *cfg)
{
   config_setting_t *s, *m;
   PSF_Matrix_Type *psf = NULL;
   const char *path_str;
   char *path = NULL;
   int ncid, dimid, start[2], count[2];
   int use_shadows;
   double scale_factor_vis_to_uv;
   size_t num_waves;
   int status = -1;

   if (CONFIG_TRUE != config_lookup_string (cfg, "calibration.straylight.psf", &path_str))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading ", __func__);
        return -1;
     }

   if (NULL == (path = expand_string (path_str)))
     return -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", path);

   if (NULL == (s = config_lookup (cfg, "straylight.psf_method")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'straylight.psf_method' in param file: %s",
                     __func__, config_error_file (cfg));
        goto return_status;
     }

   if (NULL == (m = config_setting_get_member (s, "use_shadows")))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file", __func__);
        goto return_status;
     }

   use_shadows = config_setting_get_bool (m);

   if (CONFIG_TRUE != config_setting_lookup_float (s, "scale_factor_vis_to_uv", &scale_factor_vis_to_uv))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file", __func__);
        goto return_status;
     }

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     goto return_status;

   if (0 != TIO_inq_dim (ncid, "phony_dim_0", &dimid, &num_waves))
     goto return_status;

   if (NULL == (psf = alloc_psf_matrix_type (num_waves)))
     goto return_status;

   psf->use_shadows = use_shadows;
   psf->scale_factor_vis_to_uv = scale_factor_vis_to_uv;

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_waves;

   if (0 != TIO_get_var_section (ncid, "AINV", start, count, TIO_FLOAT, psf->Ainv))
     goto return_status;

   free_psf_matrix_type (cal->sl_psf);
   cal->sl_psf = psf;

   status = 0;
return_status:
   TIO_close (ncid);
   if (status)
     {
        free_psf_matrix_type (psf);
        cal->sl_psf = NULL;
     }
   FREE(path);

   return status;
}

typedef struct
{
   double *good;
   double *good_idx;
   double *bad_idx;
   size_t num_good;
   size_t num_bad;
   gsl_interp_accel *acc;
}
Hole_Info_Type;

static void free_hole_info (Hole_Info_Type *h)
{
   if (h == NULL)
     return;
   FREE(h->good);
   gsl_interp_accel_free (h->acc);
}

static int alloc_hole_info (size_t n, Hole_Info_Type *h)
{
   if ((NULL == (h->good = (double *)MALLOC (3*n*sizeof(double))))
       || (NULL == (h->acc = gsl_interp_accel_alloc())))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_hole_info (h);
        return -1;
     }

   h->good_idx = h->good + n;
   h->bad_idx = h->good + n*2;
   h->num_good = 0;
   h->num_bad = 0;

   return 0;
}

static int interp_row (Hole_Info_Type *h, Image_Pixel_Type *pix, const Image_Pqf_Bitmap_Type *pqf,
                       Image_Pqf_Bitmap_Type mask, size_t n)
{
   gsl_interp *lin = NULL;
   size_t i, g, b;
   int status = -1;

   /* Use the good pixels to define a function, then linearly
    * interpolate values to fill in the bad pixels in this row
    */

   g = 0;
   b = 0;
   for (i = 0; i < n; i++)
     {
        if (pqf[i] & mask)
          {
             h->bad_idx[b] = i;
             b++;
          }

        else
          {
             h->good[g] = pix[i];
             h->good_idx[g] = i;
             g++;
          }
     }

   h->num_good = g;
   h->num_bad = b;

   if ((NULL == (lin = gsl_interp_alloc (gsl_interp_linear, h->num_good)))
       || (0 != gsl_interp_init (lin, h->good_idx, h->good, h->num_good)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gsl_interp_linear initialization failed", __func__);
        goto return_status;
     }

   for (b = 0; b < h->num_bad; b++)
     {
        double x = h->bad_idx[b];
        double val = gsl_interp_eval (lin, h->good_idx, h->good, x, h->acc);
        pix[(size_t)x] = (float) val;
     }

   status = 0;
return_status:

   gsl_interp_free (lin);

   return status;
}

static int fill_image_holes (Image_Type *img)
{
   Hole_Info_Type h = {0};
   Image_Pqf_Bitmap_Type mask = IMAGE_PQF_BAD_PIXEL | IMAGE_PQF_MISSING_DATA;
   size_t np = img->num_rows;
   size_t ns = img->num_cols;
   size_t s, p;

   /* Bad features are likely to extend along columns, so we'll fill in
    * bad pixels by interpolating along each row.
    */

   if (0 != alloc_hole_info (ns, &h))
     return -1;

   for (p = 0; p < np; p++)
     {
        Image_Pixel_Type *pixels = img->pixels + p * ns;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * ns;

        for (s = 0; s < ns; s++)
          {
             if (pqf[s] & mask)
               {
                  if (0 != interp_row (&h, pixels, pqf, mask, ns))
                    {
                       free_hole_info (&h);
                       return -1;
                    }
                  break;
               }
          }
     }

   free_hole_info (&h);

   return 0;
}

static int slcorr_using_bb_kernels (const Calibration_Type *cal, Image_Type *img)
{
   BB_Kernel_Type *sl_bbk = cal->sl_bbk;
   gsl_matrix_float_view Kt, S, I0, KtI0, D;
   Image_Type *img0 = NULL;
   Image_Type *kti0 = NULL;
   Image_Type *d = NULL;
   Image_Pixel_Type *pix;
   Image_Pixel_Type *pix0;
   Image_Pqf_Bitmap_Type *pqf;
   Image_Pqf_Bitmap_Type *pqf0;
   Image_Pqf_Bitmap_Type mask = IMAGE_PQF_MISSING_DATA | IMAGE_PQF_BAD_PIXEL;
   size_t num_rows = img->num_rows;
   size_t num_cols = img->num_cols;
   size_t i, num_pixels = num_rows * num_cols;
   int status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "straylight correction: broad-band kernels");

   if (NULL == (img0 = image_dup (img)))
     return -1;

   if (0 != fill_image_holes (img0))
     goto return_status;

   if ((NULL == (kti0 = image_new (sl_bbk->num_kernels, num_cols)))
       ||(NULL == (d = image_new (img->num_rows, num_cols))))
     goto return_status;

   Kt = gsl_matrix_float_view_array (sl_bbk->bb_kernels_transpose, sl_bbk->num_kernels, sl_bbk->num_waves);
   S  = gsl_matrix_float_view_array (sl_bbk->bb_stray_light, sl_bbk->num_waves, sl_bbk->num_kernels);
   I0 = gsl_matrix_float_view_array (img0->pixels, num_rows, num_cols);

   KtI0 = gsl_matrix_float_view_array (kti0->pixels, sl_bbk->num_kernels, num_cols);
   D    = gsl_matrix_float_view_array (d->pixels, num_rows, num_cols);

   /* We want to compute
    *   I = I0 - S * (transpose(K) * I0),
    * where '*' is matrix multiplication, and where the array dimensions are:
    *  [p,s] = [p,s] - [p,b] * ([b,p] * [p,s])    p=parallel, s=serial
    *        = [p,s] - [p,b] * [b,s]
    *        = [p,s] - [p,s]
    *        = [p,s]   (** showing that the dimensions work out **)
    * Note:
    *    *) K should be transposed on input so that it has dimensions [b,p].
    *    *) the original IDL expression looked different because IDL stores arrays
    *       differently and because the IDL code had additional transpose calls.
    *
    * SGEMM computes C = alpha op(A) op(B) + beta C
    */

   /* KtI0 = transpose(K) * I0 */
   gsl_blas_sgemm (CblasNoTrans, CblasNoTrans, 1.0, &Kt.matrix, &I0.matrix, 0.0, &KtI0.matrix);

   /* D = S * KtI0 */
   gsl_blas_sgemm (CblasNoTrans, CblasNoTrans, 1.0, &S.matrix, &KtI0.matrix, 0.0, &D.matrix);

   /* I0 -= D */
   gsl_matrix_float_sub (&I0.matrix, &D.matrix);

   if (0) (void) image_write_raw (kti0, "SLkti0");
   if (0) (void) image_write_raw (d, "SLdelta");

   /* Over-write the input image, leaving original bad pixel values in place. */
   pqf = img->pixel_quality_flags;
   pqf0 = img0->pixel_quality_flags;
   pix0 = img0->pixels;
   pix = img->pixels;
   for (i = 0; i < num_pixels; i++)
     {
        if (0 == (pqf0[i] & mask))
          {
             if ((pix[i] > 0) && (pix0[i] < 0))
               {
                  pqf[i] |= IMAGE_PQF_STRAYLIGHT_CORR_ERROR;
               }
             pix[i] = pix0[i];
          }
     }

   status = 0;
return_status:
   image_free (img0);
   image_free (kti0);
   image_free (d);
   return status;
}

static void shadow_free (Shadow_Type *sh)
{
   if (sh == NULL)
     return;
   FREE(sh->spec);
   FREE(sh->count);
}

static int shadow_alloc (Shadow_Type *sh, int num_rows, int first_col, int last_col)
{
   if (sh == NULL)
     return -1;

   sh->spec = NULL;
   sh->count = NULL;

   if ((NULL == (sh->spec = (float *)MALLOC (num_rows * sizeof(float))))
       || (NULL == (sh->count = (int *)MALLOC (num_rows * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        shadow_free(sh);
        return -1;
     }

   sh->first_col = first_col;
   sh->last_col = last_col + 1;   /* loop runs (i=first; i<last; i++) */
   sh->num_rows = num_rows;

   return 0;
}

static int shadow_mean_spectrum (Shadow_Type *sh, const Image_Type *img)
{
   int p, s, sb, se, np;
   float *spec = sh->spec;
   int *count = sh->count;

   np = img->num_rows;
   sb = sh->first_col;
   se = sh->last_col;

   memset ((char *)spec, 0, np * sizeof(float));

   for (p = 0; p < np; p++)
     {
        Image_Pixel_Type *img_pixel = img->pixels + p * img->num_cols;
        count[p] = 0;
        spec[p] = 0.0;
        for (s = sb; s < se; s++)
          {
             Image_Pixel_Type pix = img_pixel[s];
             if (pix != IMAGE_PIXEL_FILL_VALUE)
               {
                  spec[p] += pix;
                  count[p] += 1;
               }
          }
     }

   for (p = 0; p < np; p++)
     {
        if (count[p] > 0) spec[p] /= count[p];
     }

   return 0;
}

static int subtract_shadow1 (const Shadow_Type *sh, Image_Type *img)
{
   int p, s;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *img_pixel = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        float spec_p = sh->spec[p];
        for (s = 0; s < img->num_cols; s++)
          {
             Image_Pixel_Type pix = img_pixel[s];
             if (pix != IMAGE_PIXEL_FILL_VALUE)
               {
                  if ((pix >= 0.0) && (spec_p > pix))
                    {
                       pqf[s] |= IMAGE_PQF_STRAYLIGHT_CORR_ERROR;
                    }
                  img_pixel[s] -= spec_p;
               }
          }
     }

   return 0;
}

static int subtract_weighted_shadows (const Shadow_Type *top, const Shadow_Type *bot,
                                      Image_Type *img)
{
   float *col_weight = NULL;
   int p, s, num_not_shadowed;

   if (NULL == (col_weight = (float *)MALLOC (img->num_cols * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)col_weight, 0, img->num_cols * sizeof(float));

   num_not_shadowed = bot->first_col - top->last_col;

   for (s = top->first_col; s < top->last_col; s++)
     {
        col_weight[s] = 1.0;
     }
   for (s = top->last_col; s < bot->first_col; s++)
     {
        col_weight[s] = 1.0 - (s - top->last_col)/num_not_shadowed;
     }
   for (s = bot->first_col; s < bot->last_col; s++)
     {
        col_weight[s] = 0.0;
     }

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *img_pixel = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        float top_spec_p = top->spec[p];
        float bot_spec_p = bot->spec[p];

        for (s = 0; s < img->num_cols; s++)
          {
             Image_Pixel_Type pix = img_pixel[s];
             float spec_p, wt;

             if (pix == IMAGE_PIXEL_FILL_VALUE)
               continue;

             wt = col_weight[s];
             spec_p = top_spec_p * wt + bot_spec_p * (1.0 - wt);

             if ((pix >= 0.0) && (spec_p > pix))
               {
                  pqf[s] |= IMAGE_PQF_STRAYLIGHT_CORR_ERROR;
               }

             img_pixel[s] -= spec_p;
          }
     }

   FREE(col_weight);

   return 0;
}

static int subtract_shadows (unsigned int method, const Image_Type *img_src,
                             Image_Type *img)
{
   const Shadow_Type *sh = NULL;
   Shadow_Type top = {0};
   Shadow_Type bot = {0};
   int status = -1;

   if (method & SHADOW_TOP)
     {
        /* shadowed columns on the top/north end of the slit */
        if (0 != shadow_alloc (&top, img_src->num_rows, 0, 4))
          goto return_status;
        if (0 != shadow_mean_spectrum (&top, img_src))
          goto return_status;
        sh = &top;
     }

   if (method & SHADOW_BOT)
     {
        /* shadowed columns on the bottom/south end of the slit */
        if (0 != shadow_alloc (&bot, img_src->num_rows, 2043, 2047))
          goto return_status;
        if (0 != shadow_mean_spectrum (&bot, img_src))
          goto return_status;
        sh = &bot;
     }

   if (method == (SHADOW_BOT | SHADOW_TOP))
     {
        if (0 != subtract_weighted_shadows (&top, &bot, img))
          goto return_status;
     }
   else if (sh)
     {
        if (0 != subtract_shadow1 (sh, img))
          goto return_status;
     }
   else
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid straylight correction method", __func__);
        goto return_status;
     }

   status = 0;
return_status:
   shadow_free (&top);
   shadow_free (&bot);

   return status;
}

static int slcorr_using_shadows (const Calibration_Type *cal, Image_Type *img)
{
   unsigned int method = cal->straylight_shadow_method;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "straylight correction: shadow");
   return subtract_shadows (method, img, img);
}

static int slcorr_using_psf (const Calibration_Type *cal, Image_Type *img)
{
   PSF_Matrix_Type *psf = cal->sl_psf;
   gsl_matrix_float_view Ainv, I0, I;
   Image_Type *img0 = NULL;
   Image_Type *corr = NULL;
   Image_Pixel_Type *pix;
   Image_Pixel_Type *pix_corr;
   Image_Pqf_Bitmap_Type *pqf;
   Image_Pqf_Bitmap_Type mask = IMAGE_PQF_MISSING_DATA | IMAGE_PQF_BAD_PIXEL;
   size_t num_rows = img->num_rows;
   size_t num_cols = img->num_cols;
   size_t i, num_pixels = num_rows * num_cols;
   int status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "straylight correction: PSF");

   if (NULL == (img0 = image_dup (img)))
     goto return_status;

   if (0 != fill_image_holes (img0))
     goto return_status;

   if (NULL == (corr = image_new (img->num_rows, img->num_cols)))
     goto return_status;

   Ainv = gsl_matrix_float_view_array (psf->Ainv, psf->num_waves, psf->num_waves);
   I0 = gsl_matrix_float_view_array (img0->pixels, num_rows, num_cols);
   I  = gsl_matrix_float_view_array (corr->pixels, num_rows, num_cols);

   /* We want to compute:
    *    I = Ainv * I0,
    * where '*' is matrix multiplication, and where the array dimensions are:
    *  [p,s] = [p,p] * [p,s]      p=parallel, s=serial
    *        = [p,s]   (** showing that the dimensions work out **)
    *
    * SGEMM computes C = alpha op(A) op(B) + beta C
    */

   /* I = Ainv * I0 */
   gsl_blas_sgemm (CblasNoTrans, CblasNoTrans, 1.0, &Ainv.matrix, &I0.matrix, 0.0, &I.matrix);

   /* UV correction proportional to uncorrected VIS band signal: */
   if (psf->scale_factor_vis_to_uv > 0.0)
     {
        Image_Pixel_Type *pixels0 = img0->pixels;
        Image_Pixel_Type *pixels  = corr->pixels;
        double uv_corr, vis0_total;

        /* Compute uncorrected VIS signal */
        vis0_total = 0.0;
        for (i = 0; i < num_pixels/2; i++)
          {
             vis0_total += pixels0[i];
          }

        /* Scale UV correction */
        uv_corr = vis0_total * psf->scale_factor_vis_to_uv;

        /* Apply UV correction */
        for (i = num_pixels/2; i < num_pixels; i++)
          {
             pixels[i] -= uv_corr;
          }
     }

   /* Over-write the input image, leaving original bad pixel values in place. */
   pqf = img->pixel_quality_flags;
   pix = img->pixels;
   pix_corr = corr->pixels;
   for (i = 0; i < num_pixels; i++)
     {
        if (0 == (pqf[i] & mask))
          {
             if ((pix[i] > 0) && (pix_corr[i] < 0))
               {
                  pqf[i] |= IMAGE_PQF_STRAYLIGHT_CORR_ERROR;
               }
             pix[i] = pix_corr[i];
          }
     }

   if (psf->use_shadows)
     {
        unsigned int method = cal->straylight_shadow_method;
        if (0 != subtract_shadows (method, img0, img))
          goto return_status;
     }

   status = 0;
return_status:
   image_free (img0);
   image_free (corr);
   return status;
}

static int cal_nominal_wavelength_grid (const Calibration_Type *cal, int band_index,
                                        double *pwaves)
{
   int iw0, iw, ix, num_waves;

   /* Nominal wavelength grid is cal->wavelength_grid[waves,xpos].
    * xpos is the fastest varying index, and cal->wavelength_grid[0,*]
    * is the long wavelength end.
    */

   switch (band_index)
     {
      case TEMPO_BAND_UV:
        iw0 = cal->num_waves;   /* shortest UV */
        break;
      case TEMPO_BAND_VIS:
        iw0 = cal->num_waves/2; /* shortest VIS */
        break;
      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unrecognized band index=%d", __func__, band_index);
        return -1;
     }

   /* Define the "nominal" wavelength grid to be the wavelength grid
    * across the middle of the chip
    */
   ix = cal->num_xpos/2;
   num_waves = cal->num_waves/2;

   for (iw = 0; iw < num_waves; iw++)
     {
        float *cal_wavelen = cal->wavelength_grid + (iw0 - iw - 1) * cal->num_xpos;
        pwaves[iw] = cal_wavelen[ix];
     }

   return 0;
}

static void cal_delete (Calibration_Type *cal)
{
   if (cal == NULL)
     return;
   free_bb_kernel_type (cal->sl_bbk);
   free_psf_matrix_type (cal->sl_psf);
   btdf_free(cal->diffuser_wrk);
   btdf_free(cal->diffuser_ref);
   FREE(cal->radcal_coeffs);
   FREE(cal->wavelength_grid);
   FREE(cal);
}

#if 0
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
#endif

static int read_shadow_qualifiers (Calibration_Type *cal, config_t *cfg)
{
   config_setting_t *s;
   const char *which_shadows;

   if (NULL == (s = config_lookup (cfg, "straylight.shadow_method")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'straylight.shadow_method' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "which_shadows", &which_shadows))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file", __func__);
        return -1;
     }

   cal->straylight_shadow_method = 0;

   if (0 != strpbrk (which_shadows, "NT"))
     {
        cal->straylight_shadow_method |= SHADOW_TOP;
     }

   if (0 != strpbrk (which_shadows, "SB"))
     {
        cal->straylight_shadow_method |= SHADOW_BOT;
     }

   return 0;
}

static int config_straylight_method (Calibration_Type *cal, config_t *cfg)
{
   const char *sl_method = enable_state_query_enum (ENABLE_STRAYLIGHT);

   if (0 == strcmp (sl_method, "none"))
     {
        cal->cal_straylight_correction = NULL;
        return 0;
     }

   if (0 == strcmp (sl_method, "bb_kernel"))
     {
        if (0 != read_bb_kernels (cal, cfg))
          return -1;
        cal->cal_straylight_correction = slcorr_using_bb_kernels;
        return 0;
     }

   if (0 == strcmp (sl_method, "psf"))
     {
        if (0 != read_sl_psf_matrix (cal, cfg))
          return -1;
        cal->cal_straylight_correction = slcorr_using_psf;
        return read_shadow_qualifiers (cal, cfg);
     }

   if (0 == strncmp (sl_method, "shadow", 6))
     {
        cal->cal_straylight_correction = slcorr_using_shadows;
        return read_shadow_qualifiers (cal, cfg);
     }

   tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported straylight correction method: %s",
                __func__, sl_method);
   return -1;
}

static int optional_config_path (config_setting_t *s, const char *name, char **path)
{
   const char *str = NULL;

   *path = NULL;

   /* no such setting, or empty string is ok */
   if ((CONFIG_TRUE != config_setting_lookup_string (s, name, &str))
       || (*str == 0))
     return 0;

   if (NULL == (*path = expand_string (str)))
     return -1;

   return 0;
}

Calibration_Type *sensorcal_init (config_t *cfg, TIO_Meta_Type *meta)
{
   config_setting_t *s;
   const char *sensorcal_file = NULL;
   Calibration_Type *cal = NULL;
   char *path = NULL;
   char *radcal_trend_file = NULL;
   char *btdf_trend_file = NULL;

   int status = -1;

   if ((0 != enable_state_define (cfg, ENABLE_STRAYLIGHT))
       || (0 != enable_state_define (cfg, ENABLE_BTDF)))
     return NULL;

   if (NULL == (s = config_lookup (cfg, "calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'calibration' in param file: %s",
                     __func__, config_error_file (cfg));
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "sensorcal_file", &sensorcal_file))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading sensorcal_file", __func__);
        return NULL;
     }
   if (NULL == (path = expand_string (sensorcal_file)))
     return NULL;

   if ((0 != optional_config_path (s, "trend_rcoef", &radcal_trend_file))
       || (0 != optional_config_path (s, "trend_btdf", &btdf_trend_file)))
     goto free_and_return;

   if (NULL == (cal = (Calibration_Type *)MALLOC (sizeof *cal)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }
   memset ((char *)cal, 0, sizeof(*cal));

   cal->cal_delete = cal_delete;
   cal->cal_apply_radcal_coeffs = cal_apply_radcal_coeffs;
   cal->cal_apply_btdf = cal_apply_btdf;
   cal->cal_nominal_wavelength_grid = cal_nominal_wavelength_grid;

   if ((0 != read_radcal_coeffs (cal, path, radcal_trend_file))
       || (0 != read_wavelength_grid (cal, path))
       || (0 != meta_record_basename (meta, path))
       || (0 != meta_record_basename (meta, radcal_trend_file)))
     {
        goto free_and_return;
     }

   if (enable_state_query_bool (ENABLE_BTDF) > 0)
     {
        if ((0 != read_btdf (cal, path, btdf_trend_file))
            || (0 != meta_record_basename (meta, btdf_trend_file)))
          goto free_and_return;
     }

   if (0 != config_straylight_method (cal, cfg))
     goto free_and_return;

   status = 0;
free_and_return:
   FREE(path);
   FREE(radcal_trend_file);
   FREE(btdf_trend_file);
   if (status)
     {
        cal_delete(cal);
        cal = NULL;
     }

   return cal;
}

void sdt_free (Spectral_Data_Type *sdt)
{
   if (sdt == NULL)
     return;
   FREE(sdt->name);
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
             pixels_out[yend-y-1] = pixels_img[x + y * nx];
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
             pqf_out[yend-y-1] = img_pqf[x + y * nx];
          }
     }
}

Spectral_Data_Type *
sdt_extract_band (const Calibration_Type *cal, int band_id,
                  const Image_Type *img,
                  const Image_Type *img_err)
{
   Spectral_Data_Type *sdt = NULL;
   const char *band_name = NULL;
   int beg, end;

   switch (band_id)
     {
      case TEMPO_BAND_VIS:
        beg = 0;
        end = cal->num_waves/2;
        band_name = TEMPO_BAND_NAME_VIS;
        break;
      case TEMPO_BAND_UV:
        beg = cal->num_waves/2;
        end = cal->num_waves;
        band_name = TEMPO_BAND_NAME_UV;
        break;
      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported band index=%d", __func__, band_id);
        return NULL;
     }

   if (NULL == (sdt = sdt_alloc (img->num_cols, end-beg+1)))
     return NULL;

   if (NULL == (sdt->name = strdup (band_name)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
        sdt_free (sdt);
        return NULL;
     }

   copy_image_pixels_to_wavelength_order (img, beg, end, sdt->img);
   copy_image_pixels_to_wavelength_order (img_err, beg, end, sdt->img_err);
   copy_image_pqf_to_wavelength_order (img, beg, end, sdt->pqf);

   return sdt;
}
