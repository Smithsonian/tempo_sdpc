/* -*- mode: C; mode: fold -*- */

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

#define NUM_QUAD     4
#define NUM_OCTANTS  8

#define SAT_FUDGE_FACTOR (0.999)

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
}
CCD_Param_Type;

typedef struct
{
   float *gain;               /* [e-/DN] gain */
   float *temp;               /* [C] temperature */
   int num_gain_tfpa;
}
Gain_LUT_Type;

typedef struct
{
   Gain_LUT_Type gain_Tfpa;   /* [e-/DN] gain variation w.r.t. FPA temperature (look-up table) */
   float *gain_Tfpe_coeffs;   /* coefficients for gain variation w.r.t FPE temperature */
   int num_gain_tfpe_coeff;
   float gain_at_Tref;        /* [e-/DN] Pre-launch gain at reference FPE and FPA temperatures */
   float Tref_fpe;            /* Reference FPE temperature for gain */
   float Tref_fpa;            /* Reference FPA temperature for gain */
}
Gain_Param_Type;

typedef struct
{
   Gain_Param_Type gain;      /* gain parameters */
   float *nonlinearity_lut;   /* nonlinearity lookup table for each DN */
   int num_dn;
   float ccd_gate_limit;      /* [e-] CCD effective full-well (serial readout gate saturation limit) */
}
Octant_Response_Type;

typedef struct
{
   float *value;             /* [dimensionless] Pixel Response Non-Uniformity */
   int num_cols;
   int num_rows;
}
PRNU_Type;

typedef struct
{
   float cte;                 /* [dimensionless] */
   float readnoise_sq;        /* [e-] */
   PRNU_Type prnu;
}
Response_Info_Type;

typedef struct
{
   float mean_eoffset0[NUM_OCTANTS];   /* Mean electronic offset at FPS testing */
   float mean_eoffset[NUM_OCTANTS];    /* Mean electronic offset in the data */
   int phase_change[NUM_QUAD];
}
Phase_Change_Type;

typedef int Smear_Corr_Method_Type
(const CCD_Param_Type *, const Image_Subset_Type *,
    int, int, const void *, const Image_Type *, Image_Pixel_Type *);

#define CCD_TYPE_PRIVATE_DATA \
   CCD_Param_Type params; \
   Image_Subset_Type *psubsets; \
   Image_Subset_Type *half; \
   Image_Subset_Type *quad; \
   Image_Subset_Type *oct; \
   Response_Info_Type resp_info; \
   Phase_Change_Type pct; \
   Octant_Response_Type oct_resp_data[NUM_OCTANTS]; \
   int oct_resp_index[NUM_OCTANTS]; \
   float crosstalk_matrix[NUM_QUAD * NUM_QUAD]; \
   Smear_Corr_Method_Type *ccd_smear_correction_method;
#include "ccd.h"

static void free_gain_param_type (Gain_Param_Type *gpt)
{
   if (gpt == NULL)
     return;
   FREE(gpt->gain_Tfpa.gain);
   FREE(gpt->gain_Tfpe_coeffs);
}

static int alloc_gain_param_type (Gain_Param_Type *gpt, size_t num_gain_tfpa,
                                  size_t num_gain_tfpe_coeff)
{
   Gain_LUT_Type *glt = &gpt->gain_Tfpa;

   if (NULL == (glt->gain = (float *) MALLOC (2 * num_gain_tfpa * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   glt->temp = glt->gain + num_gain_tfpa;
   glt->num_gain_tfpa = num_gain_tfpa;

   if (NULL == (gpt->gain_Tfpe_coeffs = (float *) MALLOC (num_gain_tfpe_coeff * sizeof(float))))
     {
        FREE(glt->gain);
        glt->gain = NULL;
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
     }
   gpt->num_gain_tfpe_coeff = num_gain_tfpe_coeff;

   return 0;
}

static void free_oct_resp (Octant_Response_Type *oct_resp)
{
   if (oct_resp == NULL)
     return;
   free_gain_param_type (&oct_resp->gain);
   FREE(oct_resp->nonlinearity_lut);
}

static void free_prnu (PRNU_Type *prnu)
{
   if (prnu == NULL)
     return;
   FREE(prnu->value);
}

static void ccd_delete (CCD_Type *ccd)
{
   int i;
   if (ccd == NULL)
     return;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        free_oct_resp (&ccd->oct_resp_data[i]);
     }

   free_prnu (&ccd->resp_info.prnu);

   FREE(ccd->psubsets);
   FREE(ccd);
}

static int init_ccd_params (CCD_Type *ccd)
{
   CCD_Param_Type *p = &ccd->params;

   /* Parallel readout (spectral) dimension
    * 1046 "rows" or "lines" per quadrant */
   p->num_parallel_active = 1028;   /* number of active lines (including 4 for alignment) */
   p->num_parallel_oclock = 16;     /* number of rows overclocked for smear correction */
   p->num_parallel_sdc = 2;         /* number of storage region dark current readout rows */

   /* serial readout (spatial) dimension
    * 1056 "columns" or pixels per line */
   p->num_serial_active = 1024;     /* number of active pixels */
   p->num_serial_leading = 10;      /* number of leading buffer pixels */
   p->num_serial_trailing = 22;     /* number of trailing buffer pixels */

   /* Total number of pixels within a Level 0 image frame:
    * (2092 rows=lines) x (2112 columns) */
   p->num_parallel = 2 * (p->num_parallel_active + p->num_parallel_oclock
                          + p->num_parallel_sdc);
   p->num_serial = 2 * (p->num_serial_active
                        + p->num_serial_leading + p->num_serial_trailing);

   p->num_readout_bits = 14;
   /* number of bits per pixel in CCD readout */

   p->num_coadd_bits = 20;
   /* number of bits per pixel in co-added image */

   return 0;
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
    * row =0, col=0 is the outer corner of quadrant A on UV CCD.
    * (e.g. row=0 is minimum wavelength, col=0 is north-most)
    * Sign of row_step indicates parallel readout direction:
    *     row_step < 0 means readout toward row zero
    *     row_step > 0 means readout away from row zero
    * Sign of col_step indicates serial readout direction only for
    * quadrants and octants; abs(col_step) = 2 for octants, 1 otherwise.
    *     col_step < 0 means readout toward col zero
    *     col_step > 0 means readout away from col zero
    *
    * Sensor labeling and array indexing:
    *      |----------|----------|
    *      |          |          |
    *      |          |          |     <-- VIS ccd
    *      |          |          |
    *      |D         |C         |
    *      |----------|----------|
    *      |          |          |
    *      |          |          |     <-- UV ccd
    *      |          |          |
    *      |A         |B         |
    *      |----------|----------|  <- row 0
    *      columns --->
    *
    * Serial readout at each outer corner.
    * UV chip (AB) parallel read-out is toward shorter wavelengths ("downward")
    * VIS chip (DC) parallel read-out is toward longer wavelengths ("upward")
    *
    * The netcdf4 calibration file organizes the relevant calibration data as
    * arrays dimensioned like TSOC_eOffset(ADC,quad), with quad the fastest
    * varying index.
    * The quadrants are stored in order: (quad=0,1,2,3  => A,B,C,D)
    * The ADCs are stored in order (ADC=0,1 => odd, even).
    * Therefore, the octants are organized in this order:
    *  (Ao,Bo,Co,Do, Ae,Be,Ce,De)
    * e.g. odds first, then evens.
    * We organize the image octant-subset arrays accordingly.
    */

   image_set_subset (&ccd->half[0],  0, hr, -1,    0, nc  , +1); /* UV */
   image_set_subset (&ccd->half[1], hr, nr, +1,    0, nc  , +1); /* VIS */

   image_set_subset (&ccd->quad[0],  0, hr, -1,    0, hc  , -1); /* UV-A */
   image_set_subset (&ccd->quad[1],  0, hr, -1,   hc, nc  , +1); /* UV-B */
   image_set_subset (&ccd->quad[2], hr, nr, +1,   hc, nc  , +1); /* VIS-C */
   image_set_subset (&ccd->quad[3], hr, nr, +1,    0, hc  , -1); /* VIS-D */

   image_set_subset (&ccd->oct[0],   0, hr, -1,    1, hc  , -2); /* UV-Ao */
   image_set_subset (&ccd->oct[1],   0, hr, -1, hc+1, nc  , +2); /* UV-Bo */
   image_set_subset (&ccd->oct[2],  hr, nr, +1, hc+1, nc  , +2); /* VIS-Co */
   image_set_subset (&ccd->oct[3],  hr, nr, +1,    1, hc  , -2); /* VIS-Do */

   image_set_subset (&ccd->oct[4],   0, hr, -1,    0, hc-1, -2); /* UV-Ae */
   image_set_subset (&ccd->oct[5],   0, hr, -1, hc  , nc-1, +2); /* UV-Be */
   image_set_subset (&ccd->oct[6],  hr, nr, +1, hc  , nc-1, +2); /* VIS-Ce */
   image_set_subset (&ccd->oct[7],  hr, nr, +1,    0, hc-1, -2); /* VIS-De */

   return 0;
}

static int mean_serial_trailing_oct (const CCD_Param_Type *ccdp,
                                     const Image_Subset_Type *oct,
                                     const Image_Type *img,
                                     int num_skip, int num_selected,
                                     double *mean)
{
   int s, sb, se, p, pb, pe, num;
   double sum;

   /* Consider the trailing serial pixels in quadrant A that reads
    * out toward the left.
    * This is a lovely diagram of those num_serial_trailing pixels:
    *    <--- readout ...aaaSSSSSSSSiiiiiiiiiiUUUU
    * 'a' represents the last few photo-active pixels.
    * 'S' represents the first num_skip pixels (not included in avg)
    * 'i' represents the next num_selected pixels that will be included
    *     in the average,
    * 'U' represents any remaining unused serial trailing pixels.
    */

   if ((oct == NULL) || (img == NULL))
     return -1;

   pb = oct->row_beg;
   pe = oct->row_end;

   if (oct->col_step > 0)
     {/* B, C */
        sb = (oct->col_beg
              + ccdp->num_serial_trailing
              - num_skip
              - num_selected);
     }
   else
     {/* A, D */
        sb = (oct->col_end
              - ccdp->num_serial_trailing - 1
              + num_skip);
     }
   se = sb + num_selected;

   num = 0;
   sum = 0.0;
   for (p = pb; p < pe; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        for (s = sb; s < se; s += 2)
          {
             if (oct_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             sum += oct_pixels[s];
             num += 1;
          }
     }

   *mean = sum / num;

   return 0;
}

static int compute_mean_eoffsets (CCD_Type *ccd, const Image_Type *img)
{
   Phase_Change_Type *pct = &ccd->pct;
   int num_skip = 8;
   int num_selected = 10;
   int i;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        double mean_eoffset;
        if (0 != mean_serial_trailing_oct (&ccd->params, &ccd->oct[i], img,
                                           num_skip, num_selected, &mean_eoffset))
          {
             return -1;
          }
        pct->mean_eoffset[i] = mean_eoffset;
     }

   return 0;
}

static int configure_using_octant_phase (CCD_Type *ccd)
{
   Phase_Change_Type *pct = &ccd->pct;
   int i;

   for (i = 0; i < NUM_QUAD; i++)
     {
        if (pct->phase_change[i])
          {
             /* phase change: even and odd phases are
              * opposite what they were at FPS testing,
              * so swap the corresponding calibration
              * parameters */
             ccd->oct_resp_index[i]   = i+4;
             ccd->oct_resp_index[i+4] = i;
          }
        else
          {
             /* No phase change: even and odd phases are
              * the same as they were at FPS testing */
             ccd->oct_resp_index[i]   = i;
             ccd->oct_resp_index[i+4] = i+4;
          }
     }

   return 0;
}

static int ccd_configure_using_octant_phase (CCD_Type *ccd, const Image_Type *img)
{
   Phase_Change_Type *pct = &ccd->pct;
   float *eoff0 = pct->mean_eoffset0;
   float *eoff = pct->mean_eoffset;
   int i;

   if (0 != compute_mean_eoffsets (ccd, img))
     return -1;

   /* (ccd->phase_change[i] != 0) means yes, the odd/even state
    * in quadrant 'i' is now opposite to the state observed
    * during FPS testing */

   for (i = 0; i < NUM_QUAD; i++)
     {
        float diff0 = eoff0[i+4] - eoff0[i];
        float diff = eoff[i+4] - eoff[i];
        pct->phase_change[i] = (diff0 * diff < 0) ? 1 : 0;
     }

   return configure_using_octant_phase (ccd);
}

static int ccd_correct_coadd (const CCD_Type *ccd, int num_coadds, Image_Type *img)
{
   const CCD_Param_Type *ccdp = &ccd->params;
   int saturation_level_coadded = (1 << ccdp->num_coadd_bits) - 1;
   float saturation_threshold_coadded = SAT_FUDGE_FACTOR * saturation_level_coadded;
   int saturation_level_readout = (1 << ccdp->num_readout_bits) - 1;
   float saturation_threshold_readout = SAT_FUDGE_FACTOR * saturation_level_readout;
   Image_Pixel_Type *pixels = img->pixels;
   Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags;
   int i, num_pixels = img->num_rows * img->num_cols;

   for (i = 0; i < num_pixels; i++)
     {
        if (pixels[i] == IMAGE_PIXEL_FILL_VALUE)
          continue;

        if (pixels[i] >= saturation_threshold_coadded)
          pixel_quality_flags[i] |= IMAGE_PQF_SATURATED;

        pixels[i] /= num_coadds;

        if (pixels[i] >= saturation_threshold_readout)
          pixel_quality_flags[i] |= IMAGE_PQF_SATURATED;
     }

   return 0;
}

static int correct_offset_oct (float mean_eoffset,
                               const Image_Subset_Type *oct,
                               Image_Type *img)
{
   int s, sb0, se0, p, pb0, pe0;

   if ((oct == NULL) || (img == NULL))
     return -1;

   pb0 = oct->row_beg;
   pe0 = oct->row_end;
   sb0 = oct->col_beg;
   se0 = oct->col_end;

   for (p = pb0; p < pe0; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb0; s < se0; s += 2)
          {
             if (oct_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             oct_pixels[s] -= mean_eoffset;
             if (oct_pixels[s] < 0)
               {
                  pixel_quality_flags[s] |= IMAGE_PQF_OFFSET_CORR_ERROR;
               }
          }
     }

   return 0;
}

static int ccd_correct_offset (const CCD_Type *ccd, Image_Type *img)
{
   const Phase_Change_Type *pct = &ccd->pct;
   int i;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        if (-1 == correct_offset_oct (pct->mean_eoffset[i], &ccd->oct[i], img))
          return -1;
     }

   return 0;
}

static int correct_nonlinearity_oct (const Octant_Response_Type *oct_resp,
                                     const Image_Subset_Type *oct,
                                     Image_Type *img)
{
   int s, sb0, se0, p, pb0, pe0;
   float *lut = oct_resp->nonlinearity_lut;
   int num_dn = oct_resp->num_dn;

   pb0 = oct->row_beg;
   pe0 = oct->row_end;
   sb0 = oct->col_beg;
   se0 = oct->col_end;

   for (p = pb0; p < pe0; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb0; s < se0; s += 2)
          {
             Image_Pixel_Type tmp, pixel_value = oct_pixels[s];
             int ipixel_value;
             float frac;

             if ((pixel_value == IMAGE_PIXEL_FILL_VALUE)
                 || (pixel_value < 0))
               continue;

             ipixel_value = floor(pixel_value);
             frac = 1.0 - (pixel_value - ipixel_value);

             if (ipixel_value < num_dn-2)
               {
                  tmp = lut[ipixel_value] * frac + lut[ipixel_value+1] * (1.0 - frac);
               }
             else
               {
                  tmp = pixel_value * lut[num_dn-1] / (num_dn-1.0);
               }

             oct_pixels[s] = tmp;
             if (oct_pixels[s] < 0)
               {
                  pixel_quality_flags[s] |= IMAGE_PQF_NONLINEAR_RANGE_ERROR;
                  /* FIXME: setting the right bit? */
               }
          }
     }

   return 0;
}

static int ccd_correct_nonlinearity (const CCD_Type *ccd, Image_Type *img)
{
   int i;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        int k = ccd->oct_resp_index[i];
        if (-1 == correct_nonlinearity_oct (&ccd->oct_resp_data[k], &ccd->oct[i], img))
          return -1;
     }

   return 0;
}

static int correct_crosstalk_half (const CCD_Type *ccd, int which_half, const Image_Type *img0,
                                   float *crosstalk_vector, Image_Type *img)
{
   const Image_Subset_Type *half = &ccd->half[which_half];
   const Image_Subset_Type *q0 = &ccd->quad[which_half*2];
   const Image_Subset_Type *q1 = &ccd->quad[which_half*2 + 1];
   int s0, sb0, se0, s1, sb1, se1, p, pb, pe;
   float frac_to_lhs, frac_to_rhs;

   pb = half->row_beg;
   pe = half->row_end;

   sb0 = q0->row_beg;
   se0 = q0->row_end;
   sb1 = q1->row_beg;
   se1 = q1->row_end;

   frac_to_rhs = crosstalk_vector[0];
   frac_to_lhs = crosstalk_vector[1];

   for (p = pb; p < pe; p++)
     {
        Image_Pixel_Type *img0_pixels = img0->pixels + p * img0->num_cols;
        Image_Pixel_Type *img_pixels = img->pixels + p * img->num_cols;

        for (s0 = sb0, s1 = se1-1; (s0 < se0) && (s1 >= sb1); s0++, s1--)
          {
             if ((img0_pixels[s0] != IMAGE_PIXEL_FILL_VALUE)
                 && (img0_pixels[s1] != IMAGE_PIXEL_FILL_VALUE))
               {
                  img_pixels[s0] = img0_pixels[s0] - frac_to_lhs * img0_pixels[s1];
                  img_pixels[s1] = img0_pixels[s1] - frac_to_rhs * img0_pixels[s0];
               }
          }
     }

   return 0;
}

static int ccd_correct_crosstalk (const CCD_Type *ccd, Image_Type *img)
{
   Image_Type *img0 = NULL;
   float crosstalk_AB[2];
   float crosstalk_DC[2];

   if (NULL == (img0 = image_dup (img)))
     return -1;

   /* crosstalk matrix elements are stored in order:
    * M(from, to) where the 'to' is the fastest varying index.
    * In offset order, the matrix elements are:
    *   [ M(0,0:3), M(1,0:3), M(2,0:3), M(3,0:3) ]
    * The only coupling that matters is A <-> B, and D <-> C,
    * so the only non-zero coefficients are:
    *    M(0,1) = (A -> B) : offset=0*4+1 =  1
    *    M(1,0) = (B -> A) : offset=1*4+0 =  4
    *    M(2,3) = (C -> D) : offset=2*4+3 = 11
    *    M(3,2) = (D -> C) : offset=3*4+2 = 14
    *
    * Remember: the physical orientation of the halves
    * are {DC} and {AB} so the LHS elements are {A,D}
    * and the RHS elements are {B,C}.
    */

   crosstalk_AB[0] = ccd->crosstalk_matrix[1];
   crosstalk_AB[1] = ccd->crosstalk_matrix[4];
   crosstalk_DC[0] = ccd->crosstalk_matrix[14];   /* <-- Note the index order!! */
   crosstalk_DC[1] = ccd->crosstalk_matrix[11];

   if ((0 != correct_crosstalk_half (ccd, 0, img0, crosstalk_AB, img))
       || (0 != correct_crosstalk_half (ccd, 1, img0, crosstalk_DC, img)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: crosstalk correction failed", __func__);
        image_free (img0);
        return -1;
     }

   image_free (img0);
   return 0;
}

static int interpolate_gain_vs_Tfpa (const Gain_Param_Type *gpt, float fpa_temp, float *pgain)
{
   const Gain_LUT_Type *lut = &gpt->gain_Tfpa;
   float *gain = lut->gain;
   float *temp = lut->temp;
   int k, n = lut->num_gain_tfpa;

   if (fpa_temp < temp[0])
     {
        *pgain = gain[0];
     }
   else if (fpa_temp >= temp[n-1])
     {
        *pgain = gain[n-1];
     }
   else
     {
        float f;
        if ((k = bsearch_f (fpa_temp, temp, n)) < 0)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: binary search failed, fpa_temp=%f", __func__, fpa_temp);
             return -1;
          }
        f = (fpa_temp - temp[k]) / (temp[k+1] - temp[k]);
        *pgain = (1.0 - f) * gain[k] + f * gain[k+1];
     }

   return 0;
}

static int gain_corr_Tfpa (const Gain_Param_Type *gpt, float fpa_temp, float *corr)
{
   float gain_T;

   if (fpa_temp == gpt->Tref_fpa)
     {
        *corr = 1.0;
        return 0;
     }

   if (0 != interpolate_gain_vs_Tfpa (gpt, fpa_temp, &gain_T))
     return -1;

   *corr = gain_T / gpt->gain_at_Tref;

   return 0;
}

static int gain_corr_Tfpe (const Gain_Param_Type *gpt, float fpe_temp, float *corr)
{
   float gain_T;
   const float *c = gpt->gain_Tfpe_coeffs;

   if (fpe_temp == gpt->Tref_fpe)
     {
        *corr = 1.0;
        return 0;
     }

   if (gpt->num_gain_tfpe_coeff != 4)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: got %d coefficients describing gain dependence on FPE temperature (expected 4)",
                     __func__, gpt->num_gain_tfpe_coeff);
        return -1;
     }

   gain_T = c[0] * exp (c[1] * fpe_temp) + c[2] * exp(c[3] * fpe_temp);

   *corr = gain_T / gpt->gain_at_Tref;

   return 0;
}

static int interpolate_gain (const Gain_Param_Type *gpt,
                             float fpa_temp, float fpe_temp, float *gain)
{
   float corr_fpa, corr_fpe;

   if ((0 != gain_corr_Tfpa (gpt, fpa_temp, &corr_fpa))
       || (0 != gain_corr_Tfpe (gpt, fpe_temp, &corr_fpe)))
     return -1;

   *gain = gpt->gain_at_Tref * corr_fpa * corr_fpe;

   return 0;
}

static int correct_gain_oct (const Octant_Response_Type *oct_resp,
                             const Image_Subset_Type *oct, Image_Type *img,
                             float fpa_temp, float fpe_temp)
{
   int s, sb, se, p, pb, pe;
   float gain, saturation_threshold_gate;

   if (0 != interpolate_gain (&oct_resp->gain, fpa_temp, fpe_temp, &gain))
     return -1;

   pb = oct->row_beg;
   pe = oct->row_end;
   sb = oct->col_beg;
   se = oct->col_end;

   /* serial readout gate saturation threshold [e-] */
   saturation_threshold_gate = SAT_FUDGE_FACTOR * oct_resp->ccd_gate_limit;

   for (p = pb; p < pe; p += 1)
     {
        Image_Pixel_Type *oct_pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb; s < se; s += 2)
          {
             if (oct_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             /* gain [e-/DN] */
             oct_pixels[s] *= gain;
             if (oct_pixels[s] < 0)
               {
                  pixel_quality_flags[s] |= IMAGE_PQF_PROCESSING_ERROR;
               }
             else if (oct_pixels[s] > saturation_threshold_gate)
               {
                  pixel_quality_flags[s] |= IMAGE_PQF_SATURATED;
               }
          }
     }

   return 0;
}

static int ccd_correct_gain (const CCD_Type *ccd, Image_Type *img,
                             float fpa_temp, float fpe_temp)
{
   int i;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        int k = ccd->oct_resp_index[i];
        if (-1 == correct_gain_oct (&ccd->oct_resp_data[k], &ccd->oct[i], img,
                                    fpa_temp, fpe_temp))
          {
             return -1;
          }
     }

   return 0;
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
   int num_skip = 8, num_include = 4;
   int s, p, pb, pe;

   (void) client_data;

   /* Compute smear correction using a subset of the overclocked
    * parallel readout rows [pb, pe). */
   if (quad->row_step < 0)
     {
        pb = quad->row_end-1 - ccdp->num_parallel_oclock + num_skip;
        pe = pb + num_include;
     }
   else
     {
        pe = quad->row_beg + ccdp->num_parallel_oclock - num_skip;
        pb = pe - num_include;
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
   Image_Pixel_Type *quad_pixels = NULL;
   Image_Pqf_Bitmap_Type *quad_pqf = NULL;
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
        quad_pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb0; s < se0; s++)
          {
             if (quad_pixels[s] == IMAGE_PIXEL_FILL_VALUE)
               continue;
             quad_pixels[s] -= smear_corr[s];
             if (quad_pixels[s] < 0)
               {
                  quad_pqf[s] |= IMAGE_PQF_SMEAR_CORR_ERROR;
               }
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

   for (i = 0; i < NUM_QUAD; i++)
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
                          int num_dg_rows, int num_tg_rows,
                          float *mean_sdc_per_pixel)
{
   int s, sb0, se0, p, pb, pe, pixcount;
   float pixsum;

   if (num_tg_rows <= 0)
     {
        *mean_sdc_per_pixel = 0.0;
        return 0;
     }

   if (-1 == smear_corr_region (ccdp, quad, &sb0, &se0, NULL, NULL))
     return -1;

   /* Apparently, only one imaging region row is used
    * to hold the sum over selected storage dark current rows.
    * Which row?  (FIXME?)
    */
   if (quad->row_step < 0)
     {
        pb = quad->row_beg;
        pe = quad->row_beg + 1;
     }
   else
     {
        pb = quad->row_end - 1;
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

   /* Storage region dark current (SDC) depends on the distance
    * of the selected rows from the serial read out.
    * (e.g., ~no SDC for 290/740 nm, while full exposure
    * to SDC at 490/540 nm).
    * To correct for this, we divide by the distance in pixels
    * from the serial readout, to the midpoint of the selected
    * range of rows.  The easiest way to compute that is to just
    * use the instrument command parameters, NUM_DG_ROWS, NUM_TG_ROWS,
    * described in section 3.3 of the TEMPO ConOps, DRD SE-13,
    * doc # 2418231 rev G.
    */

   if (pixcount > 0)
     {
        float sr_distance = num_dg_rows + num_tg_rows/2.0;
        *mean_sdc_per_pixel = (((pixsum / num_tg_rows) / pixcount) / sr_distance);
     }

   return 0;
}

static int ccd_mean_storage_region_dark (const CCD_Type *ccd,
                                         const Image_Type *img,
                                         int num_dg_rows, int num_tg_rows,
                                         float *mean_sdc)
{
   int i;

   for (i = 0; i < NUM_QUAD; i++)
     {
        if (-1 == mean_sdc_quad (&ccd->params, &ccd->quad[i], img, num_dg_rows, num_tg_rows,
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

static void update_noisesq_quad (const CCD_Type *ccd, float sdc, Image_Type *noisesq,
                                 int pb, int pe, int pstep,
                                 int sb, int se, int sstep)
{
   const Response_Info_Type *rit = &ccd->resp_info;
   int num_xfer_smear = ccd->params.num_parallel_oclock;
   float sdc_scaled, sdc_factor;
   int p, s, num_xfer_p;

   /* FIXME: This expression for CTE noise is an empirical fit to
    * OMPS data.  Presumably it will eventually be replaced by a
    * TEMPO-specific expression.
    */
   for (p = pb; p < pe; p++)
     {
        Image_Pixel_Type *pnsq = noisesq->pixels + p * noisesq->num_cols;
        Image_Pqf_Bitmap_Type *pqf = noisesq->pixel_quality_flags + p * noisesq->num_cols;
        num_xfer_p = (pstep < 0) ? (1+p) : (pe-p);
        sdc_scaled = (sdc * num_xfer_p) /(pe - pb);
        sdc_factor = pow (sdc_scaled, 0.715);
        for (s = sb; s < se; s++)
          {
             int num_xfer_s, num_xfer1, num_xfer2;
             float ctesq;
             if (pqf[s] != 0)
               {
                  /* Presumably we don't need uncertainties for bad pixels */
                  pnsq[s] = IMAGE_PIXEL_FILL_VALUE;
                  continue;
               }
             num_xfer_s = (sstep < 0) ? (1+s) : (se - s);
             num_xfer1 = num_xfer_p + num_xfer_smear + num_xfer_s;
             num_xfer2 = num_xfer_p + num_xfer1;
             ctesq = (num_xfer1 * pow (pnsq[s], 0.715)
                      + num_xfer2 * sdc_factor);
             pnsq[s] += ctesq * 2.0 * (1.0 - rit->cte) + rit->readnoise_sq;
          }
     }
}

/* Assume img_noisesq has been initialized with a copy of the
 * active-region image with pixels in units of electrons */
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

   update_noisesq_quad (ccd, sdc[0], img_noisesq,    0, nr/2, -1,    0, nc/2, -1); /* UV-A */
   update_noisesq_quad (ccd, sdc[1], img_noisesq,    0, nr/2, -1, nc/2,   nc, +1); /* UV-B */
   update_noisesq_quad (ccd, sdc[2], img_noisesq, nr/2,   nr, +1, nc/2,   nc, +1); /* VIS-C */
   update_noisesq_quad (ccd, sdc[3], img_noisesq, nr/2,   nr, +1,    0, nc/2, -1); /* VIS-D */

   return 0;
}

static int ccd_correct_prnu (const CCD_Type *ccd, Image_Type *img)
{
   const PRNU_Type *pt = &ccd->resp_info.prnu;
   int s, p, image_type;

   if (IMAGE_TYPE_ACTIVE != (image_type = image_get_type (img)))
     {
        tell_verror (TELL_INTERNAL_ERROR,
                     "%s: unexpected image type %d (expected type %d)",
                     __func__, image_type, IMAGE_TYPE_ACTIVE);
        return -1;
     }

   for (p = 0; p < pt->num_rows; p++)
     {
        Image_Pixel_Type *img_pixels = img->pixels + p * img->num_cols;
        float *prnu = pt->value + p * pt->num_cols;
        for (s = 0; s < pt->num_cols; s++)
          {
             if (img_pixels[s] != IMAGE_PIXEL_FILL_VALUE)
               {
                  img_pixels[s] /= prnu[s];
               }
          }
     }

   return 0;
}

static CCD_Type *ccd_create (void)
{
   CCD_Type *ccd = NULL;
   int i;

   if (NULL == (ccd = (CCD_Type *) MALLOC (sizeof *ccd)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: allocating CCD_Type", __func__);
        return NULL;
     }
   memset ((char *)ccd, 0, sizeof *ccd);

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        memset ((char *)&ccd->oct_resp_data[i], 0, sizeof (Octant_Response_Type));
     }

   ccd->ccd_delete = ccd_delete;
   ccd->ccd_correct_coadd = ccd_correct_coadd;
   ccd->ccd_configure_using_octant_phase = ccd_configure_using_octant_phase;
   ccd->ccd_correct_offset = ccd_correct_offset;
   ccd->ccd_correct_nonlinearity = ccd_correct_nonlinearity;
   ccd->ccd_correct_crosstalk = ccd_correct_crosstalk;
   ccd->ccd_correct_gain = ccd_correct_gain;
   ccd->ccd_correct_smear = ccd_correct_smear;
   ccd->ccd_mean_storage_region_dark = ccd_mean_storage_region_dark;
   ccd->ccd_copy_active_pixels = ccd_copy_active_pixels;
   ccd->ccd_apply_pixel_quality_flags = ccd_apply_pixel_quality_flags;
   ccd->ccd_active_image_dims = ccd_active_image_dims;
   ccd->ccd_update_noisesq = ccd_update_noisesq;
   ccd->ccd_correct_prnu = ccd_correct_prnu;

   /* default methods */
   ccd->ccd_smear_correction_method = smear_correction_using_oclocks;

   return ccd;
}

static int read_octant_gain (int ncid, int adc, int quad, Gain_Param_Type *gpt)
{
   Gain_LUT_Type *glt = &gpt->gain_Tfpa;
   size_t num_gain_tfpa, num_gain_tfpe_coeff;
   int start[3], count[3], dimid;

   if ((0 != TIO_inq_dim (ncid, "n_gain_Tfpa", &dimid, &num_gain_tfpa))
       || (0 != TIO_inq_dim (ncid, "n_gain_Tfpe_coeff", &dimid, &num_gain_tfpe_coeff)))
     return -1;

   if (0 != alloc_gain_param_type (gpt, num_gain_tfpa, num_gain_tfpe_coeff))
     return -1;

   start[0] = 0;
   count[0] = num_gain_tfpa;

   if (0 != TIO_get_var_section (ncid, "Tfpas_4gain_LUT", start, count, TIO_FLOAT, glt->temp))
     return -1;

   start[0] = adc;
   start[1] = quad;
   start[2] = 0;
   count[0] = 1;
   count[1] = 1;
   count[2] = num_gain_tfpa;

   if (0 != TIO_get_var_section (ncid, "gain_Tfpa_LUT", start, count, TIO_FLOAT, glt->gain))
     return -1;

   count[2] = num_gain_tfpe_coeff;

   if (0 != TIO_get_var_section (ncid, "gain_Tfpe_coeffs", start, count, TIO_FLOAT, gpt->gain_Tfpe_coeffs))
     return -1;

   if (0 != TIO_get_var_section (ncid, "gain_refTs", start, count, TIO_FLOAT, &gpt->gain_at_Tref))
     return -1;

   if (0 != TIO_get_var_section (ncid, "ref_Tfpa_4gain", start, count, TIO_FLOAT, &gpt->Tref_fpa))
     return -1;

   if (0 != TIO_get_var_section (ncid, "ref_Tfpe_4gain", start, count, TIO_FLOAT, &gpt->Tref_fpe))
     return -1;

   return 0;
}

static int read_octant_nonlinearity_lut (int ncid, int adc, int quad, Octant_Response_Type *oct_resp)
{
   int dimid, start[3], count[3];
   size_t num_dn;

   if (0 != TIO_inq_dim (ncid, "DN", &dimid, &num_dn))
     return -1;
   oct_resp->num_dn = num_dn;

   if (NULL == (oct_resp->nonlinearity_lut = (float *)MALLOC (num_dn * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = adc;
   start[1] = quad;
   start[2] = 0;
   count[0] = 1;
   count[1] = 1;
   count[2] = num_dn;

   if (0 != TIO_get_var_section (ncid, "nonlinearity_LUT", start, count, TIO_FLOAT,
                                 oct_resp->nonlinearity_lut))
     return -1;

   return 0;
}

static int read_octant_gate_limit (int ncid, int adc, int quad, float *ccd_gate_limit)
{
   int start[2], count[2];

   start[0] = adc;
   start[1] = quad;
   count[0] = 1;
   count[1] = 1;

   return TIO_get_var_section (ncid, "ccd_gate_limit", start, count, TIO_FLOAT, ccd_gate_limit);
}

static int read_octant_params1 (int ncid, int adc, int quad, Octant_Response_Type *oct_resp)
{
   if (0 != read_octant_gain (ncid, adc, quad, &oct_resp->gain))
     return -1;

   if (0 != read_octant_nonlinearity_lut (ncid, adc, quad, oct_resp))
     return -1;

   if (0 != read_octant_gate_limit (ncid, adc, quad, &oct_resp->ccd_gate_limit))
     return -1;

   return 0;
}

static int read_octant_params (CCD_Type *ccd, const char *cal_file)
{
   int i, ncid, adc, quad;
   int status = -1;

   if (0 != TIO_open (cal_file, NC_NOWRITE, &ncid))
     return -1;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        /* i is the octant index where octants are ordered like so:
         * (Ao,Bo,Co,Do, Ae,Be,Ce,De).
         * The cal file uses indexing (adc,quad)
         * where adc=(0,1) -> (odd,even)
         * and where quad=(0,1,2,3) -> (A,B,C,D).
         * So, map 'i' to (adc,quad):
         */

        adc  = i / 4;
        quad = i % 4;

        if (0 != read_octant_params1 (ncid, adc, quad, &ccd->oct_resp_data[i]))
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: reading calibration params for octant %d (adc=%d, quad=%d)",
                          __func__, i, adc, quad);
             goto return_status;
          }
     }

   status = 0;
return_status:
   (void) TIO_close (ncid);
   return status;
}

static int read_prnu (int ncid, int num_parallel_active_full, int num_serial_active_full,
                      PRNU_Type *pt)
{
   size_t num_rows, num_cols, len, i;
   int dimid, start[2], count[2];
   float *prnu = NULL;

   if (0 != TIO_inq_dim (ncid, "wave", &dimid, &num_rows))
     return -1;
   if (0 != TIO_inq_dim (ncid, "Xpos", &dimid, &num_cols))
     return -1;

   if ((num_rows != (size_t) num_parallel_active_full)
       || (num_cols != (size_t) num_serial_active_full))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: PRNU array size mismatch: got (wave,Xpos)=(%ld,%ld), expected (%d,%d)",
                     __func__, num_rows, num_cols, num_parallel_active_full, num_serial_active_full);
        return -1;
     }

   len = num_rows * num_cols;
   if (NULL == (prnu = (float *)MALLOC (len * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_rows;
   count[1] = num_cols;

   if (0 != TIO_get_var_section (ncid, "prnu", start, count, TIO_FLOAT, prnu))
     {
        FREE(prnu);
        return -1;
     }

   /* Validate */
   for (i = 0; i < len; i++)
     {
        if ((0 == isfinite(prnu[i]))
            || (prnu[i] <= 0.0))
          {
             tell_verror (TELL_INVALID_PARM_ERROR, "%s: read invalid PRNU(%ld,%ld)=%e",
                          __func__, i / num_cols, i % num_cols, prnu[i]);
             FREE(prnu);
             return -1;
          }
     }

   pt->num_rows = num_rows;
   pt->num_cols = num_cols;
   pt->value = prnu;

   return 0;
}

static int read_cal_params (CCD_Type *ccd, const char *path)
{
   Phase_Change_Type *pct = &ccd->pct;
   Response_Info_Type *rit = &ccd->resp_info;
   float readout_noise, quantization_noise;
   int start[2], count[2];
   int num_serial_active_full, num_parallel_active_full;
   int ncid, status = -1;

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     return -1;

   start[0] = 0;
   start[1] = 0;
   count[0] = 2;  /* NUM ADCs per quadrant */
   count[1] = NUM_QUAD;

   if (0 != TIO_get_var_section (ncid, "TSOC_eOffset", start, count, TIO_FLOAT, pct->mean_eoffset0))
     goto return_status;

   start[0] = 1;
   count[0] = 1;

   if (0 != TIO_get_var_section (ncid, "cte", start, count, TIO_FLOAT, &rit->cte))
     goto return_status;
   if (0 != TIO_get_var_section (ncid, "readout_noise", start, count, TIO_FLOAT, &readout_noise))
     goto return_status;
   if (0 != TIO_get_var_section (ncid, "quantization_noise", start, count, TIO_FLOAT, &quantization_noise))
     goto return_status;

   rit->readnoise_sq = hypotf (readout_noise, quantization_noise);

   ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);

   if (0 != read_prnu (ncid, num_parallel_active_full, num_serial_active_full,
                       &ccd->resp_info.prnu))
     goto return_status;

   status = 0;

return_status:
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading calibration file %s",
                     __func__, path);
     }
   (void) TIO_close (ncid);
   return status;
}

static int read_crosstalk_matrix (CCD_Type *ccd, const char *path)
{
   int ncid, start[2], count[2];
   int status = -1;

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     return -1;

   start[0] = 0;
   start[1] = 0;
   count[0] = NUM_QUAD;
   count[1] = NUM_QUAD;

   if (0 != TIO_get_var_section (ncid, "crosstalk_matrix", start, count, TIO_FLOAT,
                                 ccd->crosstalk_matrix))
     goto return_status;

   status = 0;
return_status:
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading calibration file %s",
                     __func__, path);
     }
   (void) TIO_close (ncid);
   return status;
}

static int init_ccd_cal_params (config_t *cfg, CCD_Type *ccd, TIO_Meta_Type *meta)
{
   config_setting_t *setting;
   const char *cal_param_file;
   char *path = NULL;
   int status = -1;

   if (NULL == (setting = config_lookup (cfg, "ccd_calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_calibration in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "cal_param_file", &cal_param_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading cal_param_file in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (path = expand_path (cal_param_file)))
     return -1;

   /* Letting each subroutine have the path makes it easier
    * to switch to reading different bits from different files
    * should we ever want to do that.  Sure, I could just
    * implement it that way from the beginning, but that might
    * needlessly clutter the parameter file.
    */

   if ((0 != read_cal_params (ccd, path))
       || (0 != read_octant_params (ccd, path))
       || (0 != read_crosstalk_matrix (ccd, path)))
     goto return_status;

   if (0 != meta_record_basename (meta, path))
     goto return_status;

   status = 0;
return_status:

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

CCD_Type *ccd_init (config_t *cfg, TIO_Meta_Type *meta)
{
   CCD_Type *ccd = NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: starting", __func__);

   if (NULL == (ccd = ccd_create ()))
     return NULL;

   if (0 != init_ccd_params (ccd))
     goto error_return;

   if (0 != init_image_subsets (ccd))
     goto error_return;

   if (0 != init_ccd_cal_params (cfg, ccd, meta))
     goto error_return;

   if (-1 == init_methods (cfg, ccd))
     goto error_return;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: succeeded", __func__);

   return ccd;

error_return:
   ccd_delete (ccd);
   return NULL;
}
