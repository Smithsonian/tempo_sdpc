#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "image.h"
#include "util.h"

#define NUM_QUAD  4

#define DARK_PRIVATE_DATA \
   Image_Type *dc_image; \
   float ref_fpa_temp; \
   float dc_coeffs[NUM_QUAD];
#include "dark.h"

static void drk_close (Dark_Type *drk)
{
   if (drk == NULL)
     return;
   image_free(drk->dc_image);
   FREE(drk);
}

static int apply_quad_factor (Image_Type *img, int pb, int pe, int sb, int se, float fac)
{
   int s, p;

   for (p = pb; p < pe; p++)
     {
        Image_Pixel_Type *quad_pixel = img->pixels + p * img->num_cols;
        for (s = sb; s < se; s++)
          {
             if (quad_pixel[s] != IMAGE_PIXEL_FILL_VALUE)
               {
                  quad_pixel[s] *= fac;
               }
          }
     }

   return 0;
}

static int drk_image_Tfpa_corrected (const Dark_Type *drk, const Dark_Lookup_Type *dlt,
                                     Image_Type *img)
{
   const char *method = enable_state_query_enum (ENABLE_DARK);
   float delta_invt, fac[NUM_QUAD];
   int nr = img->num_rows;
   int nc = img->num_cols;
   int i;

   if (0 != image_copy (drk->dc_image, img))
     return -1;

   if (0 != strcmp (method, "mean_tfpa"))
     return 0;

   delta_invt = 1.0/dlt->fpa_temp - 1.0/drk->ref_fpa_temp;

   for (i = 0; i < NUM_QUAD; i++)
     {
        fac[i] = exp(drk->dc_coeffs[i] * delta_invt);
     }

   (void) apply_quad_factor (img,    0, nr/2,    0, nc/2, fac[0]);  /* A */
   (void) apply_quad_factor (img,    0, nr/2, nc/2,   nc, fac[1]);  /* B */
   (void) apply_quad_factor (img, nr/2,   nr, nc/2,   nc, fac[2]);  /* C */
   (void) apply_quad_factor (img, nr/2,   nr,    0, nc/2, fac[3]);  /* D */

   return 0;
}

static int read_image (int ncid, int k, Image_Type *img)
{
   int start[3], count[3];

   start[0] = k;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = img->num_rows;
   count[2] = img->num_cols;

   if ((0 != TIO_get_var_section (ncid, "image", start, count, TIO_FLOAT, img->pixels))
       || (0 != TIO_get_var_section (ncid, TEMPO_VAR_PQF, start, count, TIO_USHORT, img->pixel_quality_flags)))
     {
        return -1;
     }

   img->image_type = IMAGE_TYPE_ACTIVE;

   return 0;
}

static int image_add_weighted (Image_Type *img, Image_Type *weights, double weight, const Image_Type *tmp)
{
   const Image_Pqf_Bitmap_Type *pqf = tmp->pixel_quality_flags;
   const Image_Pixel_Type *pix = tmp->pixels;
   Image_Pixel_Type *wt = weights->pixels;
   Image_Pixel_Type *sum = img->pixels;
   size_t i, n = img->num_rows * img->num_cols;

   for (i = 0; i < n; i++)
     {
        if (pqf[i] == 0)
          {
             sum[i] += pix[i] * weight;
             wt[i] += weight;
          }
     }

   return 0;
}

static int image_divide (Image_Type *img, const Image_Type *denom)
{
   const Image_Pixel_Type *den = denom->pixels;
   Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags;
   Image_Pixel_Type *pix = img->pixels;
   size_t i, n = img->num_rows * img->num_cols;

   for (i = 0; i < n; i++)
     {
        if (pix[i] != IMAGE_PIXEL_FILL_VALUE)
          {
             if (isfinite(den[i]) && (den[i] != 0.0))
               {
                  pix[i] /= den[i];
               }
             else
               {
                  pix[i] = IMAGE_PIXEL_FILL_VALUE;
                  pqf[i] |= IMAGE_PQF_BAD_PIXEL;
               }
          }
     }

   return 0;
}

static int get_averaging_weights (int ncid, int num, double *wt)
{
   int k, start, count, num_zero_exposures;

   /* Assume all exposure times are >= 0.
    * If all exposure_times are zero, then weight equally when computing the mean.
    * If any exposure_times are non-zero, then weight by exposure time so that the
    *                               terms with exposure_times=0 contribute nothing.
    * If the data are packaged so that all images in a single file have the same exposure
    * time, then this should do the right thing whether the exposure time is zero or not.
    * If exposure times vary within a single file, then this also seems the right approach.
    */

   start = 0;
   count = num;
   if (0 != TIO_get_var_section (ncid, "exposure_time", &start, &count, TIO_DOUBLE, wt))
     return -1;

   num_zero_exposures = 0;
   for (k = 0; k < num; k++)
     {
        if (wt[k] == 0.0) num_zero_exposures++;
     }

   if (num == num_zero_exposures)
     {
        for (k = 0; k < num; k++)
          {
             wt[k] = 1.0;
          }
     }

   return 0;
}

static Image_Type *compute_image_mean (int ncid)
{
   TIO_Var_Info_Type info = {0};
   Image_Type *img = NULL;
   Image_Type *tmp = NULL;
   Image_Type *weights = NULL;
   double *wt = NULL;
   int k, num, status = -1;

   if (0 != TIO_inq_var (ncid, "image", &info))
     return NULL;

   if ((NULL == (img = image_new (info.dimlens[1], info.dimlens[2])))
       || (NULL == (weights = image_dup (img)))
       || (NULL == (tmp = image_dup (img))))
     goto return_status;

   num = info.dimlens[0];

   tell_vlog (TELL_MSGTYPE_INFO, 1, "averaging %d dark frames...", num);

   if (NULL == (wt = (double *)MALLOC (num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   if (0 != get_averaging_weights (ncid, num, wt))
     goto return_status;

   for (k = 0; k < num; k++)
     {
        if (0 != read_image (ncid, k, tmp))
          goto return_status;
        if (0 != image_add_weighted (img, weights, wt[k], tmp))
          goto return_status;
     }

   if (0 != image_divide (img, weights))
     goto return_status;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "done");

   status = 0;

return_status:
   image_free (tmp);
   image_free (weights);
   if (status)
     {
        image_free (img);
        img = NULL;
     }
   FREE(wt);

   return img;
}

static int drk_open (Dark_Type *drk, const char *path)
{
   char product_type[TIO_MAX_SHORT_NAME_LEN];
   int ncid, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", path);

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading",
                     __func__, path);
        return -1;
     }

   if (0 != TIO_get_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, product_type))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading product_type attribute from file %s\n",
                     __func__, path);
        goto close_and_return;
     }

   if (0 != strcmp (product_type, TEMPO_PROD_TYPE_DRK))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: no support for dark data with product_type=%s\n",
                     __func__, product_type);
        goto close_and_return;
     }

   image_free (drk->dc_image);
   if (NULL == (drk->dc_image = compute_image_mean (ncid)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computing mean dark current image from file %s\n",
                     __func__, path);
        goto close_and_return;
     }

   status = 0;
close_and_return:
   (void) TIO_close (ncid);

   return status;
}

static int read_Tfpa_coeffs (Dark_Type *drk, const char *path)
{
   int start[2], count[2], ncid, dimid, status = -1;
   size_t num_bands, num_coeffs;
   float coeffs[4];

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading",
                     __func__, path);
        return -1;
     }

   if ((0 != TIO_inq_dim (ncid, "band", &dimid, &num_bands))
       || (0 != TIO_inq_dim (ncid, "n_DC_Tfpa_coeff", &dimid, &num_coeffs)))
     goto return_status;

   if ((num_bands != 2) || (num_coeffs != 2))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unexpected array dimensions: band=%ld, n_DC_Tfpa_coeff=%ld",
                     __func__, num_bands, num_coeffs);
        goto return_status;
     }

   start[0] = 0;
   count[0] = 1;

   if (0 != TIO_get_var_section (ncid, "ref_Tfpa_4dc", start, count, TIO_FLOAT, &drk->ref_fpa_temp))
     goto return_status;

   start[0] = 0;
   start[1] = 0;
   count[0] = 2;
   count[1] = 2;

   if (0 != TIO_get_var_section (ncid, "DC_Tfpa_coeffs", start, count, TIO_FLOAT, coeffs))
     return -1;

   /* The current calibration file stores coefficients for a 2 parameter
    * correction:  a0 * exp(a1 * 1/T) for each of the two CCDs.
    * However, we will use the measured dark as the norm, so we use only the a1
    * coefficients from the file to apply the T-dependent part of the scaling.
    * Also, per-quadrant coefficients will be provided in the next iteration
    * of the calibration file, so we prepare for that by defining dc_coeffs[]
    * as an array of length 4.
    */
   drk->dc_coeffs[0] = coeffs[1];
   drk->dc_coeffs[1] = coeffs[1];
   drk->dc_coeffs[2] = coeffs[3];
   drk->dc_coeffs[3] = coeffs[3];

   status = 0;
return_status:
   (void) TIO_close (ncid);
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading %s",  __func__, path);
        /* drop */
     }

   return status;
}

static int read_dc_params (config_t *cfg, Dark_Type *drk)
{
   config_setting_t *s;
   const char *sensorcal_file = NULL;
   char *path = NULL;
   int status = -1;

   if (NULL == (s = config_lookup (cfg, "calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'calibration' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "sensorcal_file", &sensorcal_file))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading sensorcal_file", __func__);
        return -1;
     }

   if (NULL == (path = expand_string (sensorcal_file)))
     return -1;

   if (0 != read_Tfpa_coeffs (drk, path))
     goto return_status;

   status = 0;
return_status:
   FREE(path);

   return status;
}

Dark_Type *drk_init (config_t *cfg)
{
   Dark_Type *drk = NULL;

   if (enable_state_define (cfg, ENABLE_DARK) < 0)
     return NULL;

   if (NULL == (drk = (Dark_Type *)MALLOC (sizeof (*drk))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)drk, 0, sizeof (*drk));

   if (0 != read_dc_params (cfg, drk))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initialization failed", __func__);
        drk_close (drk);
        return NULL;
     }

   drk->drk_close = drk_close;
   drk->drk_open = drk_open;
   drk->drk_get_image = drk_image_Tfpa_corrected;

   return drk;
}
