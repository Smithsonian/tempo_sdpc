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
#include "granule.h"
#include "current.h"
#include "util.h"

#define NUM_QUAD  4

#define DARK_PRIVATE_DATA \
   Image_Type *dc_image; \
   float ref_fpa_temp; \
   float dc_coeffs[NUM_QUAD]; \
   float mean_sdc[NUM_QUAD];
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

static int apply_factors (Image_Type *img, const float *fac)
{
   int nr = img->num_rows;
   int nc = img->num_cols;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "dark scaling quadrants A,B,C,D:  %f %f %f %f",
              fac[0], fac[1], fac[2], fac[3]);

   (void) apply_quad_factor (img,    0, nr/2,    0, nc/2, fac[0]);  /* A */
   (void) apply_quad_factor (img,    0, nr/2, nc/2,   nc, fac[1]);  /* B */
   (void) apply_quad_factor (img, nr/2,   nr, nc/2,   nc, fac[2]);  /* C */
   (void) apply_quad_factor (img, nr/2,   nr,    0, nc/2, fac[3]);  /* D */

   return 0;
}

static int drk_image (const Dark_Type *drk, Image_Type *img)
{
   return image_copy (drk->dc_image, img);
}

static int drk_image_Tfpa_adj (const Dark_Type *drk, float fpa_temp, Image_Type *img)
{
   float delta_invt, fac[NUM_QUAD];
   int i;

   delta_invt = 1.0/fpa_temp - 1.0/drk->ref_fpa_temp;

   for (i = 0; i < NUM_QUAD; i++)
     {
        fac[i] = exp(drk->dc_coeffs[i] * delta_invt);
     }

   return apply_factors (img, fac);
}

static int drk_image_sdc_adj (const Dark_Type *drk, float *target_sdc, Image_Type *img)
{
   float fac[NUM_QUAD];
   int i;

   for (i = 0; i < NUM_QUAD; i++)
     {
        fac[i] = target_sdc[i] / drk->mean_sdc[i];
     }

   return apply_factors (img, fac);
}

static int drk_open (Dark_Type *drk, const char *path)
{
   TIO_Var_Info_Type info = {0};
   char product_type[TIO_MAX_SHORT_NAME_LEN];
   int ncid, start[2], count[2], status = -1;

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

   start[0] = 0;
   start[1] = 0;
   count[0] = 1;
   count[1] = 4;

   if (0 != TIO_get_var_section (ncid, "mean_sdc", start, count, TIO_FLOAT, drk->mean_sdc))
     goto close_and_return;

   if (0 != TIO_inq_var (ncid, "image", &info))
     goto close_and_return;

   if (NULL == (drk->dc_image = image_new (info.dimlens[1], info.dimlens[2])))
     goto close_and_return;

   if (0 != current_image_read (ncid, 0, drk->dc_image))
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
   drk->drk_image = drk_image;
   drk->drk_image_Tfpa_adj = drk_image_Tfpa_adj;
   drk->drk_image_sdc_adj = drk_image_sdc_adj;

   return drk;
}
