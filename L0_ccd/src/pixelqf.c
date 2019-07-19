#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <libconfig.h>
#include <tell.h>

#include "config.h"
#include "image.h"
#include "util.h"

#define PIXELQF_TYPE_PRIVATE_DATA \
   int num_sigmas_hot_threshold; \
   int num_sigmas_cold_threshold; \
   int transient_neighbor_hw_serial; \
   int transient_neighbor_hw_parallel; \
   double transient_threshold;
#include "pixelqf.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void pqf_delete (Pixelqf_Type *pt)
{
   if (pt == NULL)
     return;
   FREE(pt);
}

static void compute_goodpix_mean_and_stddev (const Image_Type *img,
                                             int pb, int pe, int sb, int se,
                                             double *meanp, double *stddevp)
{
   Image_Pixel_Type *pixels;
   Image_Pqf_Bitmap_Type *pqf;
   double sum, dev, sum_sqdev, mean, stddev;
   int p, s, count;

   sum = 0.0;
   count = 0;
   for (p = pb; p < pe; p++)
     {
        pixels = img->pixels + p * img->num_cols;
        pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb; s < se; s++)
          {
             if (pqf[s]) continue;
             sum += pixels[s];
             count += 1;
          }
     }
   mean = sum / count;

   sum_sqdev = 0.0;
   for (p = pb; p < pe; p++)
     {
        pixels = img->pixels + p * img->num_cols;
        pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb; s < se; s++)
          {
             if (pqf[s]) continue;
             dev = pixels[s] - mean;
             sum_sqdev += dev * dev;
             count += 1;
          }
     }
   stddev = sqrt (sum_sqdev / count);

   if (meanp) *meanp = mean;
   if (stddevp) *stddevp = stddev;
}

static int flag_hotcold (const Pixelqf_Type *pt, Image_Type *img,
                         int pb, int pe, int sb, int se)
{
   Image_Pixel_Type *pixels;
   Image_Pixel_Type hot_thresh, cold_thresh;
   double mean, stddev;
   Image_Pqf_Bitmap_Type *pqf;
   int p, s;

   compute_goodpix_mean_and_stddev (img, pb, pe, sb, se,
                                    &mean, &stddev);

   hot_thresh = mean + pt->num_sigmas_hot_threshold * stddev;

   for (p = pb; p < pe; p++)
     {
        pixels = img->pixels + p * img->num_cols;
        pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb; s < se; s++)
          {
             if ((pqf[s] == 0) && (pixels[s] > hot_thresh))
               pqf[s] |= IMAGE_PQF_HOT_PIXEL;
          }
     }

   /* Ignore hot pixels when searching for cold pixels */
   compute_goodpix_mean_and_stddev (img, pb, pe, sb, se,
                                    &mean, &stddev);

   cold_thresh = mean - pt->num_sigmas_cold_threshold * stddev;
   if (cold_thresh < 0) cold_thresh = 0.0;

   for (p = pb; p < pe; p++)
     {
        pixels = img->pixels + p * img->num_cols;
        pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb; s < se; s++)
          {
             if ((pqf[s] == 0) && (pixels[s] < cold_thresh))
               pqf[s] |= IMAGE_PQF_COLD_PIXEL;
          }
     }

   return 0;
}

static int pqf_flag_hotcold (const Pixelqf_Type *pt, Image_Type *img)
{
   int nr = img->num_rows;
   int nc = img->num_cols;

   if (enable_state_query_bool (ENABLE_HOTCOLD) < 1)
     return 0;

   flag_hotcold (pt, img,    0, nr/2,    0, nc/2);  /* A */
   flag_hotcold (pt, img, nr/2, nr  ,    0, nc/2);  /* B */
   flag_hotcold (pt, img, nr/2, nr  , nc/2, nc  );  /* C */
   flag_hotcold (pt, img,    0, nr/2, nc/2, nc  );  /* D */

   return 0;
}

typedef struct
{
   Image_Pqf_Bitmap_Type *pqf;
   int hw_serial;
   int hw_parallel;
   Image_Pqf_Bitmap_Type loc_mask;   /* mask that matches pixels to examine */
   Image_Pqf_Bitmap_Type set_mask;   /* mask indicating bits to set */
}
Neighbor_Type;

static int flag_neighbor1 (Neighbor_Type *nt, const Image_Type *img,
                           int pb, int pe, int sb, int se)
{
   Image_Pqf_Bitmap_Type loc_mask = nt->loc_mask;
   Image_Pqf_Bitmap_Type set_mask = nt->set_mask;
   int p, s;

   for (p = pb; p < pe; p++)
     {
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        Image_Pqf_Bitmap_Type *neighbor_pqf;

        for (s = sb; s < se; s++)
          {
             int pp, ss, nsb, nse, npb, npe;

             if (0 == (pqf[s] & loc_mask))
               continue;

             nsb = s - nt->hw_serial;
             nsb = MAX(sb,nsb);

             nse = s + nt->hw_serial;
             nse = MIN(se,nse);

             npb = p - nt->hw_parallel;
             npb = MAX(pb,npb);

             npe = p + nt->hw_parallel;
             npe = MIN(pe,npe);

             for (pp=npb; pp<npe; pp++)
               {
                  neighbor_pqf = nt->pqf + pp * img->num_cols;
                  for (ss=nsb; ss<nse; ss++)
                    {
                       neighbor_pqf[ss] |= set_mask;
                    }
               }
          }
     }

   return 0;
}

static int alloc_neighbor_type (const Image_Type *img, Neighbor_Type *nt)
{
   size_t pqf_array_size = (img->num_rows * img->num_cols
                            * sizeof(Image_Pqf_Bitmap_Type));

   if (NULL == (nt->pqf = (Image_Pqf_Bitmap_Type *)MALLOC (pqf_array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)nt->pqf, 0, pqf_array_size);

   return 0;
}

static void dealloc_neighbor_type (Neighbor_Type *nt)
{
   if (nt == NULL)
     return;
   FREE(nt->pqf);
}

static void flag_neighbor (const Image_Type *img, Neighbor_Type *nt)
{
   int nr = img->num_rows;
   int nc = img->num_cols;

   flag_neighbor1 (nt, img,    0, nr/2,    0, nc/2);  /* A */
   flag_neighbor1 (nt, img,    0, nr/2, nc/2, nc  );  /* B */
   flag_neighbor1 (nt, img, nr/2, nr  , nc/2, nc  );  /* C */
   flag_neighbor1 (nt, img, nr/2, nr  ,    0, nc/2);  /* D */
}

static int pqf_flag_neighbor (const Pixelqf_Type *pt, Image_Type *img,
                              int hw_serial, int hw_parallel,
                              Image_Pqf_Bitmap_Type loc_mask,
                              Image_Pqf_Bitmap_Type set_mask)
{
   Neighbor_Type nt;
   Image_Pqf_Bitmap_Type *pqf;
   int nr = img->num_rows;
   int nc = img->num_cols;
   int i;

   (void) pt;

   memset ((char *)&nt, 0, sizeof (nt));

   if (-1 == alloc_neighbor_type (img, &nt))
     return -1;

   nt.hw_serial = hw_serial;
   nt.hw_parallel = hw_parallel;
   nt.loc_mask = loc_mask;
   nt.set_mask = set_mask;

   flag_neighbor (img, &nt);

   pqf = img->pixel_quality_flags;
   for (i = 0; i < nr * nc; i++)
     {
        pqf[i] |= nt.pqf[i];
     }

   dealloc_neighbor_type (&nt);

   return 0;
}

static int flag_transients (const Pixelqf_Type *pt, Image_Type *img,
                            float *spikefs, int pb, int pe, int sb, int se)
{
   double threshold = pt->transient_threshold;
   int hw_serial = pt->transient_neighbor_hw_serial;
   int hw_parallel = pt->transient_neighbor_hw_parallel;
   Image_Pqf_Bitmap_Type *pqf;
   int s, p, num_pixels_flagged;

   num_pixels_flagged = 0;

   for (p = pb; p < pe; p++)
     {
        float *spikefs_row = spikefs + p * img->num_cols;
        pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = sb; s < se; s++)
          {
             int pp, ss, nsb, nse, npb, npe, count;
             double sum;

             if (spikefs_row[s] <= threshold)
               continue;

             nsb = s - hw_serial;
             nsb = MAX(sb,nsb);

             nse = s + hw_serial + 1;
             nse = MIN(se,nse);

             npb = p - hw_parallel;
             npb = MAX(pb,npb);

             npe = p + hw_parallel + 1;
             npe = MIN(pe,npe);

             sum = 0.0;
             count = 0;
             for (pp=npb; pp<npe; pp++)
               {
                  int o = pp * img->num_cols;
                  for (ss=nsb; ss<nse; ss++)
                    {
                       sum += spikefs[ss + o];
                       count += 1;
                    }
               }

             if (count > 1)
               {
                  double spikefs_adjacent_mean, contrast;
                  spikefs_adjacent_mean = ((sum - spikefs_row[s]) / (count - 1));
                  contrast = spikefs_row[s] / spikefs_adjacent_mean - 1.0;
                  if (contrast > threshold)
                    {
                       pqf[s] |= IMAGE_PQF_TRANSIENT_PIXEL;
                       num_pixels_flagged++;
                    }
               }
          }
     }

   if (num_pixels_flagged > 0)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "found %d transients", num_pixels_flagged);
     }

   return num_pixels_flagged;
}

static float *compute_spikefs (const Image_Pqf_Bitmap_Type *bpixmap,
                               const Image_Type *img_ref, Image_Type *img)
{
   size_t spikefs_array_size = (img->num_rows * img->num_cols
                                * sizeof(float));
   float *spikefs = NULL;
   int s, p;

   if (NULL == (spikefs = (float *) MALLOC (spikefs_array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)spikefs, 0, spikefs_array_size);

   for (p = 0; p < img->num_rows; p++)
     {
        const Image_Pqf_Bitmap_Type *bpix = bpixmap + p * img->num_cols;
        Image_Pixel_Type *pixels = img->pixels + p * img->num_cols;
        Image_Pixel_Type *ref_pixels = img_ref->pixels + p * img_ref->num_cols;
        float *spikefs_row = spikefs + p * img->num_cols;
        for (s = 0; s < img->num_cols; s++)
          {
             if ((pixels[s] > 0) && (ref_pixels[s] > 0)
                 && (bpix[s] == 0))
               {
                  spikefs_row[s] = pixels[s] / ref_pixels[s];
               }
          }
     }

   return spikefs;
}

/* Flag transient pixels by comparing to a reference image */
static int pqf_flag_transients (const Pixelqf_Type *pt,
                                const Image_Pqf_Bitmap_Type *bpixmap,
                                const Image_Type *img_ref,
                                Image_Type *img)
{
   float *spikefs;
   int nr = img->num_rows;
   int nc = img->num_cols;

   if (NULL == (spikefs = compute_spikefs (bpixmap, img_ref, img)))
     return -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "searching by quadrant for transient pixels:");

   flag_transients (pt, img, spikefs,    0, nr/2,    0, nc/2);
   flag_transients (pt, img, spikefs,    0, nr/2, nc/2, nc  );
   flag_transients (pt, img, spikefs, nr/2, nr  ,    0, nc/2);
   flag_transients (pt, img, spikefs, nr/2, nr  , nc/2, nc  );

   FREE(spikefs);

   return 0;
}

static Pixelqf_Type *pqf_create (void)
{
   Pixelqf_Type *pt;

   if (NULL == (pt = (Pixelqf_Type *) MALLOC (sizeof *pt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: allocating Pixelqf_Type", __func__);
        return NULL;
     }
   memset ((char *)pt, 0, sizeof *pt);

   pt->pqf_delete = pqf_delete;
   pt->pqf_flag_hotcold = pqf_flag_hotcold;
   pt->pqf_flag_neighbor = pqf_flag_neighbor;
   pt->pqf_flag_transients = pqf_flag_transients;

   return pt;
}

static int parse_param_file (config_t *cfg, Pixelqf_Type *pt)
{
   config_setting_t *s, *sub;

   if ((0 != enable_state_define (cfg, ENABLE_TRANSIENTS))
       || (0 != enable_state_define (cfg, ENABLE_HOTCOLD)))
     return -1;

   if (NULL == (s = config_lookup (cfg, "pixel_quality_flag_params")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing pqf_params in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "hot_thresh", &pt->num_sigmas_hot_threshold))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "cold_thresh", &pt->num_sigmas_cold_threshold)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading hot/cold pixel thresholds",
                     __func__);
        return -1;
     }

   if ((NULL == (sub = config_setting_get_member (s, "transients")))
       || (CONFIG_TRUE != config_setting_lookup_int (sub, "neighbor_hw_serial", &pt->transient_neighbor_hw_serial))
       || (CONFIG_TRUE != config_setting_lookup_int (sub, "neighbor_hw_parallel", &pt->transient_neighbor_hw_parallel))
       || (CONFIG_TRUE != config_setting_lookup_float (sub, "threshold", &pt->transient_threshold))
       )
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading transient flag parameters",
                     __func__);
        return -1;
     }

   return 0;
}

Pixelqf_Type *pixelqf_init (config_t *cfg)
{
   Pixelqf_Type *pt = NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: starting", __func__);

   if (NULL == (pt = pqf_create ()))
     return NULL;

   if (-1 == parse_param_file (cfg, pt))
     goto error_return;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: succeeded", __func__);

   return pt;
error_return:
   pqf_delete (pt);
   return NULL;
}
