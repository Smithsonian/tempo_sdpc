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

#include "granule.h"
#include "config.h"
#include "util.h"

#define DEGTORAD (M_PI/180.0)
#define RADTODEG (180.0/M_PI)

#define SLIT_AOV_STEP_RAD (41.49e-6)

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

   int do_scat_ang_corr;
   /* whether to perform scattering angle correction of BTDF */
   int do_extra_corr;
   /* whether to perform extra BTDF correction */

   int num_slope_aoi;
   float *slope_aoi;
   /* polynomial coefficients for interpolating the dependence on angle of incidence
    * (polar angle, theta) */
   float *slope_aoi_extra;
   /* Extra adjustment to BTDF sensitivity (minor correction for angle of incidence) */

   float aoi_nom;
   /* [deg] nominal angle of incidence (polar angle, theta) */

   float *btdfe_lut;
   /* [dimensionless] (num_waves,num_aov) Effective BTDF, averaged over pixel FOV */

   float lab_sf;
   /* [dimensionless] Lab results scaling factor */

   float *trend;
   /* [dimensionless] (row,col) BTDF trend multiplicative correction factor
    * row = xtrack index, col = spectral index (fastest varying)
    */
}
BTDF_Type;

typedef struct
{
   double *lpsens;   /**< linear polarization sensivity [num_xtrack, num_wave] */
   double *angmax;   /**< angle of maximum transmission [num_xtrack, num_wave] [rad] */
   int num_xtrack;
   int num_wave;
}
Lps_Table_Type;

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

static void free_lps_table (Lps_Table_Type *tbl);

#define SENSORCAL_PRIVATE_DATA \
   BB_Kernel_Type *sl_bbk; \
   PSF_Matrix_Type *sl_psf; \
   BTDF_Type *diffuser_wrk; \
   BTDF_Type *diffuser_ref; \
   Lps_Table_Type *lps_uv; \
   Lps_Table_Type *lps_vis; \
   double *diffuser_index; \
   float *wavelength_grid; \
   float *radcal_coeffs; \
   int num_waves; \
   int num_xpos; \
   unsigned int straylight_shadow_method; \
   int num_col_top; \
   int num_col_bot;
#include "sensorcal.h"

static int read_wavelength_grid_cal (Calibration_Type *cal, const char *file)
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

static int read_wavelength_grid_irr (Calibration_Type *cal, const char *irr_file)
{
   TIO_Var_Info_Type info_uv, info_vis;
   size_t num_waves_uv, num_waves_vis, num_xpos_uv, num_xpos_vis, len_uv, len_vis, len;
   float *wave_uv = NULL, *wave_vis = NULL;
   int grp_uv, grp_vis, start[3], count[3], ix, iw;
   int ncid, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", irr_file);

   if (0 != TIO_open (irr_file, NC_NOWRITE, &ncid))
     return -1;

   if ((0 != TIO_inq_grp (ncid, "band_290_490_nm", &grp_uv))
       || (0 != TIO_inq_grp (ncid, "band_540_740_nm", &grp_vis)))
     goto close_and_return;

   if ((0 != TIO_inq_var (grp_uv, "irradiance", &info_uv))
       || (0 != TIO_inq_var (grp_vis, "irradiance", &info_vis)))
     goto close_and_return;

   num_waves_uv  = info_uv.dimlens [2];
   num_xpos_uv   = info_uv.dimlens [1];
   num_waves_vis = info_vis.dimlens[2];
   num_xpos_vis  = info_vis.dimlens[1];

   if ((num_xpos_uv != num_xpos_vis)
       || (num_waves_uv != num_waves_vis))
     {
        tell_verror (TELL_IO_READ_ERROR,
                     "%s: the UV and VIS wavelength variables have different dimensions: %s",
                     __func__, irr_file);
        goto close_and_return;
     }
   else
     {
        cal->num_waves = num_waves_uv + num_waves_vis;
        cal->num_xpos  = num_xpos_uv;
     }

   len_uv  = num_waves_uv  * num_xpos_uv;
   len_vis = num_waves_vis * num_xpos_vis;

   if ((NULL == (wave_uv = (float *)MALLOC (len_uv * sizeof(float))))
       || (NULL == (wave_vis = (float *)MALLOC (len_vis * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = num_xpos_uv;
   count[2] = num_waves_uv;

   if (0 != TIO_get_var_section (grp_uv, "wavelength", start, count, TIO_FLOAT, wave_uv))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading UV wavelength grid: %s",
                     __func__, irr_file);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = num_xpos_vis;
   count[2] = num_waves_vis;

   if (0 != TIO_get_var_section (grp_vis, "wavelength", start, count, TIO_FLOAT, wave_vis))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading VIS wavelength grid: %s",
                     __func__, irr_file);
        goto close_and_return;
     }

   len = cal->num_waves * cal->num_xpos;

   if (NULL == (cal->wavelength_grid = (float *)MALLOC (len * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   for (ix = 0; ix < cal->num_xpos; ix++)
     {
        float *wave_uv_x  = wave_uv  + ix * num_waves_uv;
        float *wave_vis_x = wave_vis + ix * num_waves_vis;
        for (iw = 0; iw < cal->num_waves/2; iw++)
          {
             float *nom_uv  = cal->wavelength_grid + (cal->num_waves   - iw - 1) * cal->num_xpos;
             float *nom_vis = cal->wavelength_grid + (cal->num_waves/2 - iw - 1) * cal->num_xpos;
             nom_uv [ix] = wave_uv_x [iw];
             nom_vis[ix] = wave_vis_x[iw];
          }
     }

   status = 0;
close_and_return:
   (void) TIO_close (ncid);
   if (status)
     {
        FREE(cal->wavelength_grid);
        cal->wavelength_grid = NULL;
     }
   FREE(wave_uv);
   FREE(wave_vis);

   return 0;
}

static int read_wavelength_grid (Calibration_Type *cal, const char *file, const char *irr_file,
                                 int exposure_type)
{
   if (irr_file == NULL)
     {
        if (0 != read_wavelength_grid_cal (cal, file))
          return -1;
     }
   else
     {
        switch (exposure_type)
          {
           case EXPREC_TYPE_IRR_WRK:
           case EXPREC_TYPE_IRR_REF:
             if (0 != read_wavelength_grid_cal (cal, file))
               return -1;
             break;
           case EXPREC_TYPE_RAD:
           case EXPREC_TYPE_RAD_TWI:
             if (0 != read_wavelength_grid_irr (cal, irr_file))
               return -1;
             break;

           default:
             tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported exposure record type = %d", __func__, exposure_type);
             return -1;
          }
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
   FREE(btdf->slope_aoi_extra);
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
       || (NULL == (btdf->slope_aoi_extra = (float *)MALLOC (num_slope_aoi * sizeof(float))))
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

static BTDF_Type *read_btdf_parameters (int is_reference_diffuser, const char *file, config_t *cfg)
{
   config_setting_t *s, *m;
   BTDF_Type *btdf = NULL;
   const char *slope_aoi_name;
   const char *lut_name;
   size_t num_waves, num_aov, num_slope_aoi;
   int start[2], count[2];
   int ncid, dimid, status = -1;
   int do_extra_corr, do_scat_ang_corr;
   double lab_scaling_factor;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (NULL == (s = config_lookup (cfg, "btdf")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'btdf' in param file: %s",
                     __func__, config_error_file (cfg));
        goto close_and_return;
     }

   if (NULL == (m = config_setting_get_member (s, "do_scat_ang_corr")))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file",__func__);
        goto close_and_return;
     }

   do_scat_ang_corr = config_setting_get_bool (m);

   if (NULL == (m = config_setting_get_member (s, "do_extra_corr")))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file",__func__);
        goto close_and_return;
     }

   do_extra_corr = config_setting_get_bool (m);

   if (CONFIG_TRUE != config_setting_lookup_float (s, "lab_scaling_factor", &lab_scaling_factor))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file:", __func__);
        goto close_and_return;
     }

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if ((0 != TIO_inq_dim (ncid, "n_BTDF_w", &dimid, &num_waves))
       || (0 != TIO_inq_dim (ncid, "n_BTDF_aov", &dimid, &num_aov))
       || (0 != TIO_inq_dim (ncid, "n_BTDF_slope_aoi", &dimid, &num_slope_aoi)))
     goto close_and_return;

   if (NULL == (btdf = btdf_alloc (num_waves, num_aov, num_slope_aoi)))
     goto close_and_return;

   btdf->do_scat_ang_corr = do_scat_ang_corr;
   btdf->do_extra_corr = do_extra_corr;
   btdf->lab_sf = lab_scaling_factor;
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

   if (0 != TIO_get_var_section (ncid, "BTDF_slope_aoi_extra", start, count, TIO_FLOAT, btdf->slope_aoi_extra))
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

static int read_btdf (Calibration_Type *cal, config_t *cfg, const char *path, const char *trend_file)
{
   if ((NULL == (cal->diffuser_wrk = read_btdf_parameters (0, path, cfg)))
       || (0 != init_btdf_trend_correction (cal->diffuser_wrk, cal->num_xpos, cal->num_waves, trend_file)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading working BTDF parameters from %s",
                     __func__, path);
        return -1;
     }

   if ((NULL == (cal->diffuser_ref = read_btdf_parameters (1, path, cfg)))
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
                           Image_Type *img, Image_Type *img_diag)
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
             float aov_s_rad, wave_s, fa0, fa1, fw0, fw1, b0, b1;
             float aoi_correction, aoi_correction_extra, sa_correction, btdfe_s;
             float sf = bt->lab_sf;
             double solar_theta = solar_theta_deg * DEGTORAD;
             double sin_aoi = sin(solar_theta);
             double cos_aoi = cos(solar_theta);
             int iw, ia;

             if (img_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               {
                  if (img_diag)
                    {
                       Image_Pixel_Type *img_d = img_diag->pixels + p * img_diag->num_cols;
                       img_d[s] = IMAGE_PIXEL_FILL_VALUE;
                    }
                  continue;
               }

             aov_s_rad = SLIT_AOV_STEP_RAD * (hs - s);
             wave_s = waves[s];

             /* angle of incidendence correction */
             aoi_correction = ((bt->slope_aoi[0] + wave_s * bt->slope_aoi[1])
                               * (solar_theta_deg - bt->aoi_nom) / 100.0);

             aoi_correction_extra = ((bt->slope_aoi_extra[0] + wave_s * bt->slope_aoi_extra[1])
                                     * (bt->aoi_nom - solar_theta_deg) / 100.0);

             /* scattering angle correction */
             double aov_phi      = (aov_s_rad > 0) ? 0 : 180;
             double scat_ang     = acos(cos_aoi * cos(fabs(aov_s_rad))
                                    + sin_aoi * sin(fabs(aov_s_rad)) * cos((solar_phi_deg - aov_phi) * DEGTORAD))
                                    * RADTODEG;
             double scat_ang_nom = acos(cos_aoi * cos(fabs(aov_s_rad))
                                    + sin_aoi * sin(fabs(aov_s_rad)) * cos((        -90.0 - aov_phi) * DEGTORAD))
                                    * RADTODEG;

             sa_correction = -sf * ((bt->slope_aoi[0] + wave_s * bt->slope_aoi[1])
                                    * (scat_ang - scat_ang_nom) / 100.0);

             /* Bilinear interpolation of effective BTDF, with no extrapolation */
             double aov_s = aov_s_rad * RADTODEG;

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

             if (bt->do_scat_ang_corr)
               btdfe_s /= (1.0 + sa_correction);

             if (bt->do_extra_corr)
               btdfe_s /= (1.0 + aoi_correction_extra);

             img_pixels[s] /= btdfe_s * bt_trend[s];
             if (img_diag)
               {
                  Image_Pixel_Type *img_d = img_diag->pixels + p * img_diag->num_cols;
                  img_d[s] = 1.0 / (btdfe_s * bt_trend[s]);
               }
          }
     }

   return 0;
}

static int lps_lookup (const Calibration_Type *cal, int p, int s, double *lpsens, double *angmax)
{
   int iwave, end, offset;
   Lps_Table_Type *lps;

   /* LPS tables are split into UV, VIS bands with arrays dimensioned as
    * [num_xtrack, num_wave], so num_wave varies fastest, with wavelength
    * monotonic increasing.
    * Before splitting into bands, the CCD images are dimensioned as
    * [num_wave, num_xtrack], with wavelengths in reverse order. Ugh.
    */

   if (p < cal->num_waves/2)
     { /* VIS band */
        lps = cal->lps_vis;
        end = cal->num_waves/2;
     }
   else
     { /* UV band */
        lps = cal->lps_uv;
        end = cal->num_waves;
     }

   /* xtrack == s */
   iwave = end - p - 1;
   offset = s * lps->num_wave + iwave;

   *lpsens = lps->lpsens[offset];
   *angmax = lps->angmax[offset];

   return 0;
}

static int cal_apply_diffuser_polcorr (const Calibration_Type *cal,
                                       double solar_phi_deg, double solar_theta_deg,
                                       Image_Type *img, Image_Type *img_diag)
{
   double solar_phi = solar_phi_deg * DEGTORAD;
   double solar_theta = solar_theta_deg * DEGTORAD;
   double sin_aoi = sin(solar_theta);   /* angle of incidence (polar angle) */
   double cos_aoi = cos(solar_theta);
   int p, s;

   if (enable_state_query_bool (ENABLE_DIFF_POLCORR) < 1)
     return 0;

   /* To correct for polarization induced by the diffuser, consider the
    * Stokes vector and use Mueller matrices to model the effect of the
    * diffuser.
    * 1. Incident solar radiation is unpolarized, so the initial Stokes vector
    *    is (S0, 0, 0, 0).  (e.g. S1, S2, S3 are all zero).
    * 2. Mt' = Mueller matrix of the TEMPO instrument
    *    Md' = Mueller matrix of the diffuser
    *    Mueller matrix elements can be derived by considering the effect
    *    of the diffuser material on each component of the Stokes vector.
    *    Because the incident radiation is unpolarized, only a few matrix
    *    elements are needed.
    *    The matrix elements may be expressed in terms of the linear
    *    polarization sensitivity, LPS, of the diffuser material.
    *    The LPS may be written in terms of the Fresnel power transmission
    *    coefficients Ts, Tp, of each polarization mode (s, p).
    * 3. To account for the illumination geometry, the necessary Mueller
    *    matrices rotations are incorporated using:
    *       M' = R(x) M R(-x)
    *    where R(x) is the rotation matrix,
    *          Mt is rotated by 'chi' = angle of maximum transmission,
    *      and Md is rotated by 'az' = the azimuth angle of the incident solar
    *          irradiance.  'az' is the usual azimuthal angle in the coordinate
    *          system where:
    *          \z is the diffuser unit normal pointing toward earth,
    *          \x is a unit normal pointing northward along the slit, and
    *          \y = (\z x \x) completes the right-handed coordinate triad
    *          of unit vectors.
    * 4. The diffuser polarization correction factor is the ratio:
    *       factor = S0/S
    *    where S0 is the signal transmitted by a perfect non-polarizing diffuser,
    *      and S is the signal transmitted by the actual diffuser.
    *
    * Part of the derivation was checked with Mathematica using the following
    * expressions (Mathematica syntax):
    *
    * # Define unrotated Mueller matrices in terms of Lt,Ld which are the telescope
    * # and diffuser LPS values, respectively. Variables at,ad get multiplied by zero,
    * # so their structure is not important for the final result.
    *    Mt = {{1,Lt,0,0},{Lt,1,0,0},{0,0,at,0},{0,0,0,at}}
    *    Md = {{1,Ld,0,0},{Ld,1,0,0},{0,0,ad,0},{0,0,0,ad}}
    * # Define the rotation matrices R(+az),R(-az),R(+chi),R(-chi)
    * # using abbreviations z=azimuth, x=chi
    *    Rpz = {{1,0,0,0},{0,Cos[2z],-Sin[2z],0},{0, Sin[2z],Cos[2z],0},{0,0,0,1}}
    *    Rmz = {{1,0,0,0},{0,Cos[2z], Sin[2z],0},{0,-Sin[2z],Cos[2z],0},{0,0,0,1}}
    *    Rpx = {{1,0,0,0},{0,Cos[2x],-Sin[2x],0},{0, Sin[2x],Cos[2x],0},{0,0,0,1}}
    *    Rmx = {{1,0,0,0},{0,Cos[2x], Sin[2x],0},{0,-Sin[2x],Cos[2x],0},{0,0,0,1}}
    * # Compute rotated Mueller matrices
    *    tprime = Rpx.Mt.Rmx
    *    dprime = Rpz.Md.Rmz
    *    signal = tprime . dprime
    *    signal[[1,1]]
    *    Out[10]= 1 + Ld Lt Cos[2 x] Cos[2 z] + Ld Lt Sin[2 x] Sin[2 z]
    * Here, 'signal' is really the ratio -- e.g., the incident Stokes vector
    * is (1,0,0,0).
    *
    * Finally, we use cos(a-b) = cos(a)cos(b) + sin(a)*sin(b)
    */

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *img_pixels = img->pixels + p * img->num_cols;
        double *diffuser_index = cal->diffuser_index + p * img->num_cols;

        for (s = 0; s < img->num_cols; s++)
          {
             double sin_t, cos_t, r, x, lps_d, lps_t, ang_max, diffuser_polcorr;
             double n2 = diffuser_index[s];  /* index of refraction */

             if (img_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               {
                  if (img_diag)
                    {
                       Image_Pixel_Type *img_d = img_diag->pixels + p * img_diag->num_cols;
                       img_d[s] = IMAGE_PIXEL_FILL_VALUE;
                    }
                  continue;
               }

             /* Snell's law, angle of incidence -> angle of transmission */
             sin_t = sin_aoi / n2;
             cos_t = (fabs(sin_t) > 1.0) ? 0.0 : sqrt (1.0 - sin_t * sin_t);

             /* The two-plate diffuser linear polarization sensitivity, LPS, may
              * be expressed as lps_d = (Tp^2-Ts^2)/(Tp^2+Ts^2) where
              * Ts, Tp are power transmission coefficients for perpendicular (s),
              * and parallel polarization modes (p), respectively.
              * [Not to be confused with serial/parallel readout of the CCD].
              * For details on Ts,Tp, see e.g. the wikipedia article on Fresnel Equations.
              * We can write lps_d in terms of the ratio, x = Ts/Tp.
              * Ts and Tp have a common numerator and can be written in the form:
              *     Ts = A/a^2, Tp=A/b^2,
              * so that x = Ts/Tp = b^2/a^2. Defining r=b/a, we then have:
              */
             r = (cos_t + n2 * cos_aoi) / (cos_aoi + n2 * cos_t);
             x = r*r;
             lps_d = (1.0 - x*x) / (1.0 + x*x);

             /* lps_t = spectrometer linear polarization sensitivity */
             (void) lps_lookup (cal, p, s, &lps_t, &ang_max);

             /* Factor to correct for polarization caused by the diffuser. */
             diffuser_polcorr = 1.0/(1.0 + lps_d * lps_t * cos(2*(ang_max - solar_phi)));

             /* apply the correction */
             img_pixels[s] *= diffuser_polcorr;
             if (img_diag)
               {
                  Image_Pixel_Type *img_d = img_diag->pixels + p * img_diag->num_cols;
                  img_d[s] = diffuser_polcorr;
               }
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
   gsl_matrix_float_view Ainv;
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

   if (0 != TIO_inq_dim (ncid, "wave", &dimid, &num_waves))
     goto return_status;

   if (NULL == (psf = alloc_psf_matrix_type (num_waves)))
     goto return_status;

   psf->use_shadows = use_shadows;
   psf->scale_factor_vis_to_uv = scale_factor_vis_to_uv;

   start[0] = 0;
   start[1] = 0;
   count[0] = num_waves;
   count[1] = num_waves;

   /* After this read, psf->Ainv actually holds transpose(Ainv).
    * Immediately after the read, we'll transpose the array in-place. */
   if (0 != TIO_get_var_section (ncid, "PSF_1dSpectral_AInv", start, count, TIO_FLOAT, psf->Ainv))
     goto return_status;

   /* Transpose psf->Ainv in place */
   Ainv = gsl_matrix_float_view_array (psf->Ainv, psf->num_waves, psf->num_waves);
   if (0 != gsl_matrix_float_transpose (&Ainv.matrix))
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
    * interpolate values to fill in the bad pixels in this row.
    * In this process, exclude the first five and the last seven
    * columns since they are "shadow" pixels.
    */

   g = 0;
   b = 0;
   for (i = 5; i < (n-7); i++)
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
   int p, s;
   float num_not_shadowed;

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

static int subtract_shadows (unsigned int method, int num_col_top, int num_col_bot,
                             const Image_Type *img_src, Image_Type *img)
{
   const Shadow_Type *sh = NULL;
   Shadow_Type top = {0};
   Shadow_Type bot = {0};
   int status = -1;

   if (method & SHADOW_TOP)
     {
        /* shadowed columns on the top/north end of the slit */
        if (0 != shadow_alloc (&top, img_src->num_rows, 0, (num_col_top-1)))
          goto return_status;
        if (0 != shadow_mean_spectrum (&top, img_src))
          goto return_status;
        sh = &top;
     }

   if (method & SHADOW_BOT)
     {
        /* shadowed columns on the bottom/south end of the slit */
        if (0 != shadow_alloc (&bot, img_src->num_rows, (2048-num_col_bot), 2047))
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
   int num_col_top = cal->num_col_top;
   int num_col_bot = cal->num_col_bot;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "straylight correction: shadow");
   return subtract_shadows (method, num_col_top, num_col_bot, img, img);
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
        int num_col_top = cal->num_col_top;
        int num_col_bot = cal->num_col_bot;
        if (0 != subtract_shadows (method, num_col_top, num_col_bot, img0, img))
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

   num_waves = cal->num_waves/2;

   for (ix = 0; ix < cal->num_xpos; ix++)
     {
        double *pwaves_x = pwaves + ix * num_waves;
        for (iw = 0; iw < num_waves; iw++)
          {
             float *cal_wavelen = cal->wavelength_grid + (iw0 - iw - 1) * cal->num_xpos;
             pwaves_x[iw] = cal_wavelen[ix];
          }
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
   free_lps_table (cal->lps_uv);
   free_lps_table (cal->lps_vis);
   FREE(cal->radcal_coeffs);
   FREE(cal->wavelength_grid);
   FREE(cal->diffuser_index);
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

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_col_top", &cal->num_col_top))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file", __func__);
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_col_bot", &cal->num_col_bot))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading config file", __func__);
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

static void free_lps_table (Lps_Table_Type *tbl)
{
   if (tbl == NULL)
     return;
   FREE(tbl->lpsens);
   FREE(tbl->angmax);
   FREE(tbl);
}

static Lps_Table_Type *alloc_lps_table (int num_xtrack, int num_wave)
{
   Lps_Table_Type *tbl = NULL;
   size_t len = num_xtrack * num_wave * sizeof(double);

   if (NULL == (tbl = (Lps_Table_Type *)MALLOC (sizeof *tbl)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)tbl, 0, sizeof (*tbl));

   if ((NULL == (tbl->lpsens = (double *)MALLOC (len)))
       || (NULL == (tbl->angmax = (double *)MALLOC (len))))
     {
        free_lps_table (tbl);
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   tbl->num_xtrack = num_xtrack;
   tbl->num_wave = num_wave;

   return tbl;
}

static Lps_Table_Type *read_lps_table (int grp)
{
   Lps_Table_Type *tbl = NULL;
   TIO_Var_Info_Type info = {0};
   const char lps_var[] = "linear_polarization_sensitivity";
   const char angmax_var[] = "angle_of_maximum_transmission";
   int k, num_xtrack, num_wave;
   int start[3], count[3];

   if (0 != TIO_inq_var (grp, lps_var, &info))
     return NULL;

   /* dimlen[0] is slowest varying,
    * dimlen[ndims-1] is fastest varying */
   num_xtrack = info.dimlens[1];
   num_wave = info.dimlens[2];

   if (NULL == (tbl = alloc_lps_table (num_xtrack, num_wave)))
     return NULL;

   /* read the table for the middle mirror position */
   start[0] = 1;
   start[1] = 0;
   start[2] = 0;

   count[0] = 1;
   count[1] = num_xtrack;
   count[2] = num_wave;

   if ((0 != TIO_get_var_section (grp, lps_var, start, count, NC_DOUBLE, tbl->lpsens))
       ||(0 != TIO_get_var_section (grp, angmax_var, start, count, NC_DOUBLE, tbl->angmax)))
     goto return_error;

   /* convert deg -> radians */
   for (k = 0; k < num_xtrack * num_wave; k++)
     {
        tbl->angmax[k] *= DEGTORAD;
     }

   return tbl;
return_error:
   free_lps_table (tbl);
   return NULL;
}

static int read_lps (Calibration_Type *cal, config_t *cfg)
{
   const char *lps_file = NULL;
   char *lps_path = NULL;
   int status = -1;
   int ncid, grp;

   if (CONFIG_TRUE != config_lookup_string (cfg, "calibration.lps_file", &lps_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading calibration.lps_file: %s", __func__, config_error_file(cfg));
        return -1;
     }

   if (NULL == (lps_path = expand_string (lps_file)))
     return -1;

   if (0 != TIO_open (lps_path, NC_NOWRITE, &ncid))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, TEMPO_BAND_NAME_UV, &grp))
     goto return_status;

   if (NULL == (cal->lps_uv = read_lps_table (grp)))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, TEMPO_BAND_NAME_VIS, &grp))
     goto return_status;

   if (NULL == (cal->lps_vis = read_lps_table (grp)))
     goto return_status;

   if (   (cal->num_waves/2 != cal->lps_uv->num_wave)
       || (cal->num_waves/2 != cal->lps_vis->num_wave)
       || (cal->num_xpos != cal->lps_uv->num_xtrack)
       || (cal->num_xpos != cal->lps_vis->num_xtrack)
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: LPS table dimensions are inconsistent", __func__);
        goto return_status;
     }

   status = 0;
return_status:
   FREE(lps_path);
   return status;
}

static int diffuser_index_of_refraction (Calibration_Type *cal)
{
   /* Sellmeier dispersion equation parameters from:
    * "Interspecimen Comparison of the Refractive Index of Fused Silica",
    * Malitson, 1965, J. of Opt. Soc. Am., 55, 1205 */
   double a[] = {0.6961663, 0.4079426, 0.8974794};
   double y0[] = {0.0684043, 0.1162414, 9.896161};  /* microns */
   double micron_per_nm = 1.e-3;
   int n = sizeof(y0)/sizeof(*y0);
   int i, k, num = cal->num_waves * cal->num_xpos;

   if (NULL == (cal->diffuser_index = (double *)MALLOC (num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < num; i++)
     {
        double y = micron_per_nm * cal->wavelength_grid[i];
        double nsqr = 1.0;
        for (k = 0; k < n; k++)
          {
             double r = y0[k]/y;
             nsqr += a[k] / (1.0 - r * r);
          }
        cal->diffuser_index[i] = sqrt(nsqr);
     }

   return 0;
}

Calibration_Type *sensorcal_init (config_t *cfg, TIO_Meta_Type *meta, const char *irr_file,
                                  int exposure_type)
{
   config_setting_t *s;
   const char *sensorcal_file = NULL;
   Calibration_Type *cal = NULL;
   char *path = NULL;
   char *irr_path = NULL;
   char *radcal_trend_file = NULL;
   char *btdf_trend_file = NULL;

   int status = -1;

   if ((0 != enable_state_define (cfg, ENABLE_STRAYLIGHT))
       || (0 != enable_state_define (cfg, ENABLE_BTDF))
       || (0 != enable_state_define (cfg, ENABLE_DIFF_POLCORR)))
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
   cal->cal_apply_diffuser_polcorr = cal_apply_diffuser_polcorr;
   cal->cal_nominal_wavelength_grid = cal_nominal_wavelength_grid;

   if (irr_file != NULL)
     {
        if (NULL == (irr_path = expand_string (irr_file)))
          return NULL;
     }

   if ((0 != read_radcal_coeffs (cal, path, radcal_trend_file))
       || (0 != read_wavelength_grid (cal, path, irr_path, exposure_type)))
     {
        goto free_and_return;
     }

   if (0)
     {
        if ((0 != meta_record_basename (meta, path))
            || (0 != meta_record_basename (meta, radcal_trend_file)))
          goto free_and_return;
     }

   if (enable_state_query_bool (ENABLE_BTDF) > 0)
     {
        if (0 != read_btdf (cal, cfg, path, btdf_trend_file))
          goto free_and_return;

        if (enable_state_query_bool (ENABLE_DIFF_POLCORR) > 0)
          {
             if (0 != read_lps (cal, cfg))
               goto free_and_return;
             if (0 != diffuser_index_of_refraction (cal))
               goto free_and_return;
          }
     }

   if (0 != config_straylight_method (cal, cfg))
     goto free_and_return;

   status = 0;
free_and_return:
   FREE(path);
   FREE(irr_path);
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
       || (NULL == (sdt->wave = (double *)MALLOC (img_size * sizeof(double)))))
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

   if (NULL == (sdt = sdt_alloc (img->num_cols, end-beg)))
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
