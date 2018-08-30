#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libconfig.h>
#include <tell.h>

#include "config.h"
#include "util.h"
#include "image.h"

typedef struct
{
   int num_serial;              /* total number of serial readout pixels */
   int num_parallel;            /* total number of parallel readout pixels */

   /* serial readout */
   int num_serial_active;       /* number of photo-active pixels per quadrant */
   int num_serial_leading;      /* number of leading buffer pixels per quadrant */
   int num_serial_trailing;     /* number of trailing buffer pixels per quadrant */

   /* parallel readout */
   int num_parallel_active;     /* number of photo-active pixels per quadrant */
   int num_parallel_oclock;     /* number of pixels overclocked for smear correction per quadrant */
   int num_parallel_sdc;        /* number of pixels containing storage region dark current per quadrant */

   /* pixel values */
   int num_readout_bits;        /* number of bits per pixel in CCD readout */
   int num_coadd_bits;          /* number of bits per pixel in coadded image */
   int num_full_well;           /* maximum number of e- a CCD pixel can hold */
}
CCD_Param_Type;

#define NUM_OCTANTS  8
#define MAX_NUM_NLCOEFFS 6

typedef struct
{
   double gain;                        /* [DN/e-] */
   double offset;                      /* [DN] */
   double cte;                         /* [dimensionless] */
   double readout_noise;               /* [e-] */
   double quant_noise;                 /* [e-] */
   double nlcoeffs[MAX_NUM_NLCOEFFS];
   int num_nlcoeffs;
}
Gain_Param_Type;

typedef int Smear_Corr_Method_Type
(const CCD_Param_Type *, const Image_Subset_Type *,
    int, int, const void *, const Image_Type *, Image_Pixel_Type *);

#define CCD_TYPE_PRIVATE_DATA \
   CCD_Param_Type params; \
   Image_Subset_Type *psubsets; \
   Image_Subset_Type *half; \
   Image_Subset_Type *quad; \
   Image_Subset_Type *oct; \
   Gain_Param_Type gain_params[NUM_OCTANTS]; \
   int num_octants; \
   Smear_Corr_Method_Type *ccd_smear_correction_method;
#include "ccd.h"

static void ccd_delete (CCD_Type *ccd)
{
   if (ccd == NULL)
     return;
   FREE(ccd->psubsets);
   FREE(ccd);
}

static int init_image_subsets (CCD_Type *ccd)
{
   CCD_Param_Type *p = &ccd->params;
   int nr, nc, hr, hc;

   /* Do this only once, because the image size is invariant. */
   if (ccd->psubsets != NULL)
     return 0;

   if (NULL == (ccd->psubsets = image_new_subsets (2+4+8)))
     return -1;

   ccd->half = ccd->psubsets;
   ccd->quad = ccd->half + 2;
   ccd->oct  = ccd->quad + 4;

   nr = p->num_parallel;
   nc = p->num_serial;
   hr = nr/2;
   hc = nc/2;

   /* Image subsets are specified with X_beg <= i < X_end, for X=row|col.
    * row =0, col=0 is upper left of quadrant A on UV CCD.
    * (e.g. row=0 is minimum wavelength, col=0 is north-most)
    * Sign of row_step indicates parallel readout direction:
    *     row_step < 0 means readout toward row zero
    *     row_step > 0 means readout away from row zero
    * Sign of col_step indicates serial readout direction only
    * for quadrants and octants; abs(col_step) = 2 for octants, 1 otherwise.
    *     col_step < 0 means readout toward col zero
    *     col_step > 0 means readout away from col zero
    */

   image_set_subset (&ccd->half[0],  0, hr, -1,    0, nc  , +1); /* UV */
   image_set_subset (&ccd->half[1], hr, nr, +1,    0, nc  , +1); /* VIS */

   image_set_subset (&ccd->quad[0],  0, hr, -1,    0, hc  , -1); /* UV-A */
   image_set_subset (&ccd->quad[1],  0, hr, -1,   hc, nc  , +1); /* UV-B */
   image_set_subset (&ccd->quad[2], hr, nr, +1,    0, hc  , -1); /* VIS-C */
   image_set_subset (&ccd->quad[3], hr, nr, +1,   hc, nc  , +1); /* VIS-D */

   image_set_subset (&ccd->oct[0],   0, hr, -1,    0, hc-1, -2); /* UV-Ae */
   image_set_subset (&ccd->oct[1],   0, hr, -1,    1, hc  , -2); /* UV-Ao */
   image_set_subset (&ccd->oct[2],   0, hr, -1, hc  , nc-1, +2); /* UV-Be */
   image_set_subset (&ccd->oct[3],   0, hr, -1, hc+1, nc  , +2); /* UV-Bo */

   image_set_subset (&ccd->oct[4],  hr, nr, +1,    0, hc-1, -2); /* VIS-Ce */
   image_set_subset (&ccd->oct[5],  hr, nr, +1,    1, hc  , -2); /* VIS-Co */
   image_set_subset (&ccd->oct[6],  hr, nr, +1, hc  , nc-1, +2); /* VIS-De */
   image_set_subset (&ccd->oct[7],  hr, nr, +1, hc+1, nc  , +2); /* VIS-Do */

   return 0;
}

static int ccd_correct_coadd (const CCD_Type *ccd, int num_coadds, Image_Type *img)
{
   const CCD_Param_Type *ccdp = &ccd->params;
   int saturation_level_readout = (1 << ccdp->num_readout_bits) - 1;
   int saturation_level_coadded = (1 << ccdp->num_coadd_bits) - 1;
   Image_Pixel_Type *pixels = img->pixels;
   Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags;
   int i, num_pixels = img->num_rows * img->num_cols;

   for (i = 0; i < num_pixels; i++)
     {
        if (pixels[i] == IMAGE_PIXEL_FILL_VALUE)
          continue;

        if (pixels[i] >= saturation_level_coadded)
          pixel_quality_flags[i] |= IMAGE_PQF_SATURATED;

        pixels[i] /= num_coadds;

        if (pixels[i] >= saturation_level_readout)
          pixel_quality_flags[i] |= IMAGE_PQF_SATURATED;
     }

   return 0;
}

static int correct_offset_oct (const CCD_Param_Type *ccdp,
                               const Image_Subset_Type *oct,
                               Image_Type *img)
{
   int unused_serial_trailing = 1;
   int oct_num_serial_trailing = (ccdp->num_serial_trailing/2
                                  - unused_serial_trailing);
   int s, sb, se, sb0, se0, p, pb0, pe0;

   if ((oct == NULL) || (img == NULL))
     return -1;

   pb0 = oct->row_beg;
   pe0 = oct->row_end;
   sb0 = oct->col_beg;
   se0 = oct->col_end;

   if (oct->col_step > 0)
     {
        sb = oct->col_beg;
        se = oct->col_beg + oct_num_serial_trailing;
     }
   else
     {
        sb = oct->col_end - oct_num_serial_trailing - 1;
        se = oct->col_end;
     }
   if (0) fprintf (stderr, "sb=%4d se=%4d\n", sb, se);

   for (p = pb0; p < pe0; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        double offset = 0.0;
        int count = 0;
        for (s = sb; s < se; s += 2)
          {
             if (oct_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             offset += oct_pixels[s];
             count += 1;
          }
        offset /= count;
        for (s = sb0; s < se0; s += 2)
          {
             if (oct_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             oct_pixels[s] -= offset;
          }
     }

   return 0;
}

static int ccd_correct_offset (const CCD_Type *ccd, Image_Type *img)
{
   int i;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        if (-1 == correct_offset_oct (&ccd->params, &ccd->oct[i], img))
          return -1;
     }

   return 0;
}

static int correct_nonlinearity_oct (const Gain_Param_Type *gpt,
                                     const Image_Subset_Type *oct,
                                     Image_Type *img)
{
   int s, sb0, se0, p, pb0, pe0;
   int num_coeffs = gpt->num_nlcoeffs;
   const double *coeffs = gpt->nlcoeffs;

   pb0 = oct->row_beg;
   pe0 = oct->row_end;
   sb0 = oct->col_beg;
   se0 = oct->col_end;

   for (p = pb0; p < pe0; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        for (s = sb0; s < se0; s += 2)
          {
             Image_Pixel_Type pixel_value = oct_pixels[s];
             Image_Pixel_Type tmp = 0.0;
             Image_Pixel_Type tmpx = 1.0;
             int k;
             if (pixel_value == IMAGE_PIXEL_FILL_VALUE)
               continue;
             for (k = 0; k < num_coeffs; k++)
               {
                  tmp += coeffs[k] * tmpx;
                  tmpx *= pixel_value;
               }
             oct_pixels[s] = tmp;
          }
     }

   return 0;
}

static int ccd_correct_nonlinearity (const CCD_Type *ccd, Image_Type *img)
{
   int i;

   for (i = 0; i < ccd->num_octants; i++)
     {
        if (-1 == correct_nonlinearity_oct (&ccd->gain_params[i], &ccd->oct[i], img))
          return -1;
     }

   return 0;
}

static int correct_gain_oct (const Gain_Param_Type *gpt,
                             const Image_Subset_Type *oct,
                             Image_Type *img)
{
   int s, sb0, se0, p, pb0, pe0;
   double offset = gpt->offset;
   double gain = gpt->gain;

   pb0 = oct->row_beg;
   pe0 = oct->row_end;
   sb0 = oct->col_beg;
   se0 = oct->col_end;

   for (p = pb0; p < pe0; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        for (s = sb0; s < se0; s += 2)
          {
             if (oct_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             oct_pixels[s] = (oct_pixels[s] - offset) / gain;
          }
     }

   return 0;
}

static int flag_gain_corrected_saturation (const CCD_Type *ccd,
                                           Image_Type *img)
{
   const CCD_Param_Type *ccdp = &ccd->params;
   Image_Pixel_Type *pixels = img->pixels;
   Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags;
   int num_pixels = img->num_rows * img->num_cols;
   int num_full_well = ccdp->num_full_well;
   int i;

   for (i = 0; i < num_pixels; i++)
     {
        if (pixels[i] >= num_full_well)
          pixel_quality_flags[i] |= IMAGE_PQF_SATURATED;
     }

   return 0;
}

static int ccd_correct_gain (const CCD_Type *ccd, Image_Type *img)
{
   int i;

   for (i = 0; i < ccd->num_octants; i++)
     {
        if (-1 == correct_gain_oct (&ccd->gain_params[i], &ccd->oct[i], img))
          return -1;
     }

   return flag_gain_corrected_saturation (ccd, img);
}

static int smear_corr_region (const CCD_Param_Type *ccdp,
                              const Image_Subset_Type *quad,
                              int *sb0p, int *se0p, int *pb0p, int *pe0p)
{
   int sb0, se0, pb0, pe0;

   /* Operate on serial readout columns [sb0,se0).
    * Apply smear correction to parallel readout rows [pb0, pe0).
    */

   sb0 = quad->col_beg;
   se0 = quad->col_end;
   if (quad->col_step < 0)
     {
        sb0 += ccdp->num_serial_leading;
        se0 -= ccdp->num_serial_trailing;
     }
   else
     {
        sb0 += ccdp->num_serial_trailing;
        se0 -= ccdp->num_serial_leading;
     }

   if (quad->row_step < 0)
     {
        pb0 = quad->row_beg + ccdp->num_parallel_sdc;
        pe0 = quad->row_end - ccdp->num_parallel_oclock;
     }
   else
     {
        pb0 = quad->row_beg + ccdp->num_parallel_oclock;
        pe0 = quad->row_end - ccdp->num_parallel_sdc;
     }

   if (sb0p) *sb0p = sb0;
   if (se0p) *se0p = se0;
   if (pb0p) *pb0p = pb0;
   if (pe0p) *pe0p = pe0;

   return 0;
}

static int smear_correction_using_oclocks (const CCD_Param_Type *ccdp,
                                           const Image_Subset_Type *quad,
                                           int sb0, int se0,
                                           const void *client_data,
                                           const Image_Type *img,
                                           Image_Pixel_Type *smear_corr)
{
   Image_Pixel_Type *quad_pixels = img->pixels;
   int img_num_cols = img->num_cols;
   int s, p, pb, pe;

   (void) client_data;

   /* Compute smear correction using parallel readout rows [pb, pe). */
   if (quad->row_step < 0)
     {
        pb  = quad->row_end - ccdp->num_parallel_oclock;
        pe  = quad->row_end;
     }
   else
     {
        pb  = quad->row_beg;
        pe  = pb + ccdp->num_parallel_oclock;
     }

   for (s = sb0; s < se0; s++)
     {
        Image_Pixel_Type pixsum = 0.0;
        int pixcount = 0;
        for (p = pb; p < pe; p++)
          {
             Image_Pixel_Type pixel_value = quad_pixels[s + p * img_num_cols];
             if (pixel_value == IMAGE_PIXEL_FILL_VALUE)
               continue;
             pixsum += pixel_value;
             pixcount += 1;
          }
        smear_corr[s] = pixsum / pixcount;
     }

   return 0;
}

static int smear_correction_using_timing (const CCD_Param_Type *ccdp,
                                          const Image_Subset_Type *quad,
                                          int sb0, int se0,
                                          const void *client_data,
                                          const Image_Type *img,
                                          Image_Pixel_Type *smear_corr)
{
   Image_Pixel_Type *quad_pixels = img->pixels;
   int img_num_cols = img->num_cols;
   int s, p, pb0, pe0;
   double smear_fraction;

   if (client_data == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: smear fraction not specified (NULL pointer)",
                     __func__);
        return -1;
     }

   smear_fraction = *(double *)client_data;

   if (-1 == smear_corr_region (ccdp, quad, NULL, NULL, &pb0, &pe0))
     return -1;

   for (s = sb0; s < se0; s++)
     {
        Image_Pixel_Type pixsum = 0.0;
        int pixcount = 0;
        for (p = pb0; p < pe0; p++)
          {
             Image_Pixel_Type pixel_value = quad_pixels[s + p * img_num_cols];
             if (pixel_value == IMAGE_PIXEL_FILL_VALUE)
               continue;
             pixsum += pixel_value;
             pixcount += 1;
          }
        smear_corr[s] = (pixsum / pixcount) * smear_fraction;
     }

   return 0;
}

static int correct_smear_quad (const CCD_Param_Type *ccdp,
                               const Image_Subset_Type *quad,
                               Smear_Corr_Method_Type *calc_correction,
                               const void *client_data, Image_Type *img)
{
   Image_Pixel_Type *quad_pixels = img->pixels;
   Image_Pixel_Type *smear_corr = NULL;
   int s, sb0, se0, p, pb0, pe0;

   if (-1 == smear_corr_region (ccdp, quad, &sb0, &se0, &pb0, &pe0))
     return -1;

   if (0) fprintf (stderr, "smear:  sb0=%4d se0=%4d pb0=%4d pe0=%4d\n",
                   sb0, se0, pb0, pe0);

   if (NULL == (smear_corr = (Image_Pixel_Type *) MALLOC (img->num_cols * sizeof(Image_Pixel_Type))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   if (-1 == (*calc_correction)(ccdp, quad, sb0, se0, client_data, img, smear_corr))
     {
        FREE(smear_corr);
        return -1;
     }

   for (p = pb0; p < pe0; p++)
     {
        quad_pixels = img->pixels + p * img->num_cols;
        for (s = sb0; s < se0; s++)
          {
             if (quad_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             quad_pixels[s] -= smear_corr[s];
          }
     }

   FREE(smear_corr);

   return 0;
}

static int ccd_correct_smear (const CCD_Type *ccd, const void *client_data,
                              Image_Type *img)
{
   int i;

   if (ccd->ccd_smear_correction_method == NULL)
     {
        tell_verror (TELL_INTERNAL_ERROR,
                     "%s: smear correction method not specified!",
                     __func__);
        return -1;
     }

   for (i = 0; i < 4; i++)
     {
        if (-1 == correct_smear_quad (&ccd->params, &ccd->quad[i],
                                      ccd->ccd_smear_correction_method,
                                      client_data, img))
          return -1;
     }

   return 0;
}

typedef struct
{
   const char *name;
   Smear_Corr_Method_Type *method;
}
Smear_Corr_Method_Entry_Type;
static Smear_Corr_Method_Entry_Type Smear_Corr_Methods[] =
{
   {"oclocks", smear_correction_using_oclocks},
   {"timing", smear_correction_using_timing},
   {NULL, NULL}
};

int __ccd_set_smear_corr_method (CCD_Type *ccd, const char *name)
{
   Smear_Corr_Method_Entry_Type *m;

   for (m = Smear_Corr_Methods; m->name != NULL; m++)
     {
        if (0 == strcasecmp (m->name, name))
          {
             tell_vlog (TELL_MSGTYPE_INFO, 1, "smear correction method: %s", name);
             ccd->ccd_smear_correction_method = m->method;
             return 0;
          }
     }

   tell_verror (TELL_INVALID_PARM_ERROR,
                "%s: unsupported smear correction method = %s",
                __func__, name);

   return -1;
}

static int mean_sdc_quad (const CCD_Param_Type *ccdp,
                          const Image_Subset_Type *quad,
                          const Image_Type *img,
                          float *mean_sdc_per_pixel)
{
   int s, sb0, se0, p, pb, pe, pixcount;
   float pixsum;

   if (-1 == smear_corr_region (ccdp, quad, &sb0, &se0, NULL, NULL))
     return -1;

   /* [**] Only one SDC of the rows is actually used */
   if (quad->row_step < 0)
     {
        pb = quad->row_beg;
        pe = quad->row_beg + 1;  /* [**] pe = quad->row_beg + ccdp->num_parallel_sdc; */
     }
   else
     {
        pb = quad->row_end - 1; /* [**] pb = quad->row_end - ccdp->num_parallel_sdc; */
        pe = quad->row_end;
     }

   pixsum = 0.0;
   pixcount = 0;
   for (p = pb; p < pe; p++)
     {
        const Image_Pixel_Type *quad_pixels = img->pixels + p * img->num_cols;
        for (s = sb0; s < se0; s++)
          {
             if (quad_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             pixsum += quad_pixels[s];
             pixcount += 1;
          }
     }

   *mean_sdc_per_pixel = (pixsum / pixcount) / ccdp->num_parallel_active;

   return 0;
}

static int ccd_mean_storage_region_dark (const CCD_Type *ccd,
                                         const Image_Type *img,
                                         float mean_sdc[4])
{
   int i;

   for (i = 0; i < 4; i++)
     {
        if (-1 == mean_sdc_quad (&ccd->params, &ccd->quad[i], img,
                                 &mean_sdc[i]))
          return -1;
     }

   return 0;
}

static int map_active_pixels (const CCD_Type *ccd, int elem_size,
                              char *active, int active_num_cols,
                              char *padded, int padded_num_cols,
                              int (*process_row)(void *, void *, size_t, int))
{
   const CCD_Param_Type *ccdp = &ccd->params;
   int padded_top_offset, padded_bottom_offset, padded_quad_row_offset;
   int padded_offset, active_offset;
   int num_active_rows, p;
   size_t num_serial_active;

   num_active_rows = 2 * ccdp->num_parallel_active;

   padded_top_offset = (ccdp->num_parallel_sdc * padded_num_cols
                        + ccdp->num_serial_leading);
   padded_bottom_offset = (padded_top_offset
                           + 2 * ccdp->num_parallel_oclock * padded_num_cols);
   padded_quad_row_offset = (ccdp->num_serial_active
                             + 2 * ccdp->num_serial_trailing);
   num_serial_active = ccdp->num_serial_active;

   /* top half */
   for (p = 0; p < num_active_rows/2; p++)
     {
        /* left quad */
        padded_offset = p * padded_num_cols + padded_top_offset;
        active_offset = p * active_num_cols;
        if (process_row (active + active_offset * elem_size,
                         padded + padded_offset * elem_size,
                         num_serial_active, elem_size)) return -1;
        /* right quad */
        padded_offset += padded_quad_row_offset;
        active_offset += ccdp->num_serial_active;
        if (process_row (active + active_offset * elem_size,
                         padded + padded_offset * elem_size,
                         num_serial_active, elem_size)) return -1;
     }

   /* bottom half */
   for (p = num_active_rows/2; p < num_active_rows; p++)
     {
        /* left quad */
        padded_offset = p * padded_num_cols + padded_bottom_offset;
        active_offset = p * active_num_cols;
        if (process_row (active + active_offset * elem_size,
                         padded + padded_offset * elem_size,
                         num_serial_active, elem_size)) return -1;
        /* right quad */
        padded_offset += padded_quad_row_offset;
        active_offset += ccdp->num_serial_active;
        if (process_row (active + active_offset * elem_size,
                         padded + padded_offset * elem_size,
                         num_serial_active, elem_size)) return -1;
     }

   return 0;
}

static int copy_active_from_padded (void *p_active, void *p_padded, size_t num_elem, int elem_size)
{
   size_t len = num_elem * elem_size;
   memcpy ((char *)p_active, (char *)p_padded, len);
   return 0;
};

static Image_Type *ccd_copy_active_pixels (const CCD_Type *ccd,
                                           const Image_Type *img)
{
   const CCD_Param_Type *ccdp = &ccd->params;
   int image_type, num_active_rows, num_active_cols;
   Image_Type *aimg = NULL;

   if (IMAGE_TYPE_PADDED != (image_type = image_get_type (img)))
     {
        tell_verror (TELL_INTERNAL_ERROR,
                     "%s: unexpected image type %d (expected type %d)",
                     __func__, image_type, IMAGE_TYPE_PADDED);
        return NULL;
     }

   num_active_rows = 2 * ccdp->num_parallel_active;
   num_active_cols = 2 * ccdp->num_serial_active;

   if (NULL == (aimg = image_new (num_active_rows, num_active_cols)))
     return NULL;
   image_set_type (aimg, IMAGE_TYPE_ACTIVE);

   (void) map_active_pixels (ccd, sizeof(*img->pixels),
                             (char *)aimg->pixels, aimg->num_cols,
                             (char *)img->pixels, img->num_cols,
                             copy_active_from_padded);
   (void) map_active_pixels (ccd, sizeof(*img->pixel_quality_flags),
                             (char *)aimg->pixel_quality_flags, aimg->num_cols,
                             (char *)img->pixel_quality_flags, img->num_cols,
                             copy_active_from_padded);

   return aimg;
}

static void ccd_active_image_dims (const CCD_Type *ccd,
                                   int *num_parallel_active_full,
                                   int *num_serial_active_full)
{
   const CCD_Param_Type *ccdp = &ccd->params;

   if (num_parallel_active_full)
     *num_parallel_active_full = 2*ccdp->num_parallel_active;

   if (num_serial_active_full)
     *num_serial_active_full = 2*ccdp->num_serial_active;
}

static int apply_active_flag_array_to_padded (void *p_active, void *p_padded,
                                              size_t num_elem, int size_elem)
{
   Image_Pqf_Bitmap_Type *active = (Image_Pqf_Bitmap_Type *)p_active;
   Image_Pqf_Bitmap_Type *padded = (Image_Pqf_Bitmap_Type *)p_padded;
   size_t i;

   (void) size_elem;

   for (i = 0; i < num_elem; i++)
     {
        padded[i] |= active[i];
     }

   return 0;
}

static int ccd_apply_pixel_quality_flags (const CCD_Type *ccd, Image_Type *img,
                                          const Image_Pqf_Bitmap_Type *flags,
                                          int num_rows, int num_cols)
{
   int num_parallel_active_full, num_serial_active_full;
   int image_type;

   if (IMAGE_TYPE_PADDED != (image_type = image_get_type (img)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected image, image_type=%d (expected a padded image, image_type=%d)",
                     __func__, image_type, IMAGE_TYPE_PADDED);
        return -1;
     }

   ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);

   if ((num_rows != num_parallel_active_full)
       || (num_cols != num_serial_active_full))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mismatched quality flag array dimensions: got (r:%d,c:%d), expected (r:%d,c:%d)",
                     __func__, num_rows, num_cols, num_parallel_active_full, num_serial_active_full);
        return -1;
     }

   return map_active_pixels (ccd, sizeof(Image_Pqf_Bitmap_Type),
                             (char *)flags, num_cols,
                             (char *)img->pixel_quality_flags, img->num_cols,
                             apply_active_flag_array_to_padded);
}

static void update_noisesq_quad (const CCD_Type *ccd, int quad,
                                 float sdc, Image_Type *noisesq,
                                 int pb, int pe, int pstep,
                                 int sb, int se, int sstep)
{
   const CCD_Param_Type *ccdp = &ccd->params;
   const Gain_Param_Type *gpt0 = &ccd->gain_params[2*quad];
   const Gain_Param_Type *gpt1 = &ccd->gain_params[2*quad+1];
   float cte, readout_noise, quant_noise, readnoise_sq;
   int p, s, num_xfer_smear;

   num_xfer_smear = ccdp->num_parallel_oclock;
   cte = 0.5 * (gpt0->cte + gpt1->cte);
   readout_noise = 0.5 * (gpt0->readout_noise + gpt1->readout_noise);
   quant_noise = 0.5 * (gpt0->quant_noise + gpt1->quant_noise);
   readnoise_sq = (readout_noise*readout_noise
                   + quant_noise*quant_noise);

   /* FIXME: This expression for CTE noise is an empirical fit to
    * OMPS data.  Presumably it will eventually be replaced by a
    * TEMPO-specific expression.
    */
   for (p = pb; p < pe; p++)
     {
        Image_Pixel_Type *pnsq = noisesq->pixels + p * noisesq->num_cols;
        Image_Pqf_Bitmap_Type *pqf = noisesq->pixel_quality_flags + p * noisesq->num_cols;
        float sdc_scaled, sdc_factor;
        int num_xfer_p = 1 + ((pstep < 0) ? p : (pe - p));
        sdc_scaled = (sdc * num_xfer_p) /(pe - pb);
        sdc_factor = pow (sdc_scaled, 0.715);
        for (s = sb; s < se; s++)
          {
             int num_xfer_s, num_xfer1, num_xfer2;
             float ctesq;
             if (pqf[s])
               {
                  /* Presumably we don't need uncertainties for bad pixels */
                  pnsq[s] = IMAGE_PIXEL_FILL_VALUE;
                  continue;
               }
             num_xfer_s = 1 + ((sstep < 0) ? s : (se - s));
             num_xfer1 = num_xfer_p + num_xfer_smear + num_xfer_s;
             num_xfer2 = num_xfer_p + num_xfer1;
             ctesq = (num_xfer1 * pow (pnsq[s], 0.715)
                      + num_xfer2 * sdc_factor);
             pnsq[s] += ctesq * 2.0 * (1.0 - cte) + readnoise_sq;
          }
     }
}

/* Assume img_noisesq has been initialized with a copy of the
 * active-region image, smear-corrected, with pixels in units of electrons */
static int ccd_update_noisesq (const CCD_Type *ccd, const float *sdc,
                               Image_Type *img_noisesq)
{
   int nr = img_noisesq->num_rows;
   int nc = img_noisesq->num_cols;
   int image_type;

   if (IMAGE_TYPE_ACTIVE != (image_type = image_get_type (img_noisesq)))
     {
        tell_verror (TELL_INTERNAL_ERROR,
                     "%s: unexpected image type %d (expected type %d)",
                     __func__, image_type, IMAGE_TYPE_ACTIVE);
        return -1;
     }

   update_noisesq_quad (ccd, 0, sdc[0], img_noisesq,    0, nr/2, -1,    0, nc/2, -1); /* UV-A */
   update_noisesq_quad (ccd, 1, sdc[1], img_noisesq,    0, nr/2, -1, nc/2,   nc, +1); /* UV-B */
   update_noisesq_quad (ccd, 2, sdc[2], img_noisesq, nr/2,   nr, +1,    0, nc/2, -1); /* VIS-C */
   update_noisesq_quad (ccd, 3, sdc[3], img_noisesq, nr/2,   nr, +1, nc/2,   nc, +1); /* VIS-D */

   return 0;
}

static CCD_Type *ccd_create (void)
{
   CCD_Type *ccd = NULL;

   if (NULL == (ccd = (CCD_Type *) MALLOC (sizeof *ccd)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: allocating CCD_Type", __func__);
        return NULL;
     }
   memset ((char *)ccd, 0, sizeof *ccd);

   ccd->ccd_delete = ccd_delete;
   ccd->ccd_correct_coadd = ccd_correct_coadd;
   ccd->ccd_correct_offset = ccd_correct_offset;
   ccd->ccd_correct_nonlinearity = ccd_correct_nonlinearity;
   ccd->ccd_correct_gain = ccd_correct_gain;
   ccd->ccd_correct_smear = ccd_correct_smear;
   ccd->ccd_mean_storage_region_dark = ccd_mean_storage_region_dark;
   ccd->ccd_copy_active_pixels = ccd_copy_active_pixels;
   ccd->ccd_apply_pixel_quality_flags = ccd_apply_pixel_quality_flags;
   ccd->ccd_active_image_dims = ccd_active_image_dims;
   ccd->ccd_update_noisesq = ccd_update_noisesq;

   /* default methods */
   ccd->ccd_smear_correction_method = smear_correction_using_oclocks;

   return ccd;
}

typedef struct
{
   const char *name;
   size_t name_offset;
}
Param_Table_Type;
#define PARAM_ENTRY(name) {#name, offsetof(CCD_Param_Type, name)}
#define PARAM_TABLE_END   {NULL, 0}
static Param_Table_Type Param_Table[] =
{
   PARAM_ENTRY(num_serial_active),
   PARAM_ENTRY(num_serial_leading),
   PARAM_ENTRY(num_serial_trailing),
   PARAM_ENTRY(num_parallel_active),
   PARAM_ENTRY(num_parallel_oclock),
   PARAM_ENTRY(num_parallel_sdc),
   PARAM_ENTRY(num_readout_bits),
   PARAM_ENTRY(num_coadd_bits),
   PARAM_ENTRY(num_full_well),
   PARAM_TABLE_END
};

static int read_ccd_params (config_setting_t *s, CCD_Param_Type *params)
{
   Param_Table_Type *p;

   for (p = Param_Table; p->name != NULL; p++)
     {
        int *addr = (int *)((char *)params + p->name_offset);
        if (CONFIG_TRUE != config_setting_lookup_int (s, p->name, addr))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading parameter %s",
                          __func__, p->name);
             return -1;
          }
     }

   return 0;
}

static int parse_param_file (config_t *cfg, CCD_Type *ccd)
{
   config_setting_t *setting;

   if (NULL == (setting = config_lookup (cfg, "ccd_parameters")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return read_ccd_params (setting, &ccd->params);
}

static int init_subsets (CCD_Type *ccd)
{
   CCD_Param_Type *p = &ccd->params;

   p->num_serial = 2 * (p->num_serial_active
                        + p->num_serial_leading + p->num_serial_trailing);

   p->num_parallel = 2 * (p->num_parallel_active
                          + p->num_parallel_oclock + p->num_parallel_sdc);

   if (-1 == init_image_subsets (ccd))
     return -1;

   return 0;
}

static int read_dn_to_charge_params (CCD_Type *ccd, const char *gain_file)
{
   FILE *fp;

   if (gain_file == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: gain_file=NULL", __func__);
        return -1;
     }

   if (NULL == (fp = fopen (gain_file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: reading %s",
                     __func__, gain_file);
        return -1;
     }
   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", gain_file);

   ccd->num_octants = 0;
   while (ccd->num_octants < NUM_OCTANTS)
     {
#define BUFSIZE  256
        char buf[BUFSIZE];
        Gain_Param_Type si;
        int id, num_fields;
        if (NULL == fgets (buf, sizeof(buf), fp))
          {
             if (feof(fp)) break;
             continue;
          }
        if (*buf == '#')
          continue;
        si.num_nlcoeffs = MAX_NUM_NLCOEFFS;
        num_fields = sscanf (buf, "%d,%le,%le,%le,%le,%le, %le,%le,%le,%le,%le,%le",
                             &id,
                             &si.gain,
                             &si.offset,
                             &si.cte,
                             &si.readout_noise,
                             &si.quant_noise,
                             &si.nlcoeffs[0],
                             &si.nlcoeffs[1],
                             &si.nlcoeffs[2],
                             &si.nlcoeffs[3],
                             &si.nlcoeffs[4],
                             &si.nlcoeffs[5]
                            );
        if (num_fields != 12)
          {
             tell_verror (TELL_INVALID_DATA_ERROR, "%s: parsing gain file: %s",
                          __func__, gain_file);
             break;
          }
        if (id < 0 || NUM_OCTANTS <= id)
          {
             tell_verror (TELL_INVALID_PARM_ERROR, "%s: invalid quadrant id=%d in gain file: %s",
                          __func__, id, gain_file);
             break;
          }
        ccd->gain_params[id] = si; /* struct copy */
        ccd->num_octants++;
     }

   (void) fclose (fp);

   if (ccd->num_octants < NUM_OCTANTS)
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading gain file: %s",
                     __func__, gain_file);
        return -1;
     }

   return 0;
}

static int init_gain (config_t *cfg, CCD_Type *ccd)
{
   config_setting_t *setting;
   const char *gain_file;
   char *path = NULL;
   int status;

   if (NULL == (setting = config_lookup (cfg, "ccd_calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_calibration in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "gain_file", &gain_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading gain_file in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (path = expand_path (gain_file)))
     return -1;
   status = read_dn_to_charge_params (ccd, path);
   FREE(path);

   return status;
}

static int init_methods (config_t *cfg, CCD_Type *ccd)
{
   config_setting_t *setting;
   const char *smear_method;

   if (NULL == (setting = config_lookup (cfg, "ccd_methods")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_methods in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "smear_corr", &smear_method))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading smear correction method in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (-1 == __ccd_set_smear_corr_method (ccd, smear_method))
     return -1;

   return 0;
}

CCD_Type *ccd_init (config_t *cfg)
{
   CCD_Type *ccd = NULL;

   if (NULL == (ccd = ccd_create ()))
     return NULL;

   if (-1 == parse_param_file (cfg, ccd))
     goto error_return;

   if (-1 == init_subsets (ccd))
     goto error_return;

   if (-1 == init_gain (cfg, ccd))
     goto error_return;

   if (-1 == init_methods (cfg, ccd))
     goto error_return;

   return ccd;

error_return:
   ccd_delete (ccd);
   return NULL;
}
