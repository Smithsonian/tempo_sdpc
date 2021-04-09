#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libconfig.h>
#include <image.h>
#include <pixelqf.h>

#include "util.h"

#define MALLOC malloc
#define FREE free

#define NUM_ROWS 2000
#define NUM_COLS NUM_ROWS

static void image_set_diag (Image_Type *img, Image_Pixel_Type value)
{
   int i, n;

   n = (img->num_rows < img->num_cols) ? img->num_rows : img->num_cols;

   for (i = 0; i < n; i++)
     {
        Image_Pixel_Type *pixels = img->pixels + i * img->num_cols;
        pixels[i] = value;
     }
}

static int image_test_pqf_diag (Image_Type *img,
                                Image_Pqf_Bitmap_Type mask_diag)
{
   int r, c;

   for (r = 0; r < img->num_rows; r++)
     {
        Image_Pqf_Bitmap_Type *pqf =
          img->pixel_quality_flags + r * img->num_cols;
        for (c = 0; c < img->num_cols; c++)
          {
             if ((r != c) && (pqf[c] != 0))
               {
                  fprintf (stderr, "%s: r=%d c=%d pqf=%d (expected pqf=0)\n",
                           __func__, r, c, pqf[c]);
                  return 0;
               }
             else if ((r == c) && (pqf[c] != mask_diag))
               {
                  fprintf (stderr, "%s: r=%d c=%d pqf=%d (expected pqf=%d)\n",
                           __func__, r, c, pqf[c], mask_diag);
                  return 0;
               }
          }
     }

   return 1;
}

static int test_hotcold (Pixelqf_Type *pt)
{
   Image_Type *img;
   Image_Pixel_Type p0 = 1024.0;
   Image_Pixel_Type sigma = sqrt(p0);
   Image_Pixel_Type hot_value = p0 + 6*sigma;
   Image_Pixel_Type cold_value = p0 - 6*sigma;
   Dark_Trend_Type dtt = {0};
   int status = -1;

   if (NULL == (img = image_new (NUM_ROWS, NUM_COLS)))
     return -1;

   /* Note that there is a correlation between the
    * mean pixel value, the hot/cold pixel threshold,
    * and the size of the test image.  For small images,
    * the diagonal contains a larger fraction of the
    * total number of pixels, so a line of hot pixels
    * on the diagonal has a large effect on the mean value.
    */

   image_set (img, p0, 0);
   image_set_diag (img, hot_value);
   if (0 != pt->pqf_flag_hotcold (pt, img, &dtt))
     goto return_status;
   if (0 == image_test_pqf_diag (img, IMAGE_PQF_HOT_PIXEL))
     {
        fprintf (stderr, "*** hot pixel flagging failed\n");
        goto return_status;
     }

   image_set (img, p0, 0);
   image_set_diag (img, cold_value);
   if (0 != pt->pqf_flag_hotcold (pt, img, &dtt))
     goto return_status;
   if (0 == image_test_pqf_diag (img, IMAGE_PQF_COLD_PIXEL))
     {
        fprintf (stderr, "*** cold pixel flagging failed\n");
        goto return_status;
     }

   status = 0;
return_status:
   image_free (img);
   return status;
}

static int test_neighbor1 (Pixelqf_Type *pt, Image_Type *img,
                           int r0, int c0, int hw)
{
   Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags;
   Image_Pqf_Bitmap_Type pqf_rc;
   int num_rows = img->num_rows;
   int num_cols = img->num_cols;
   int loc_mask = IMAGE_PQF_MISSING_DATA;
   int set_mask = IMAGE_PQF_MISSING_DATA;
   int row_hw = hw;
   int col_hw = hw;
   int beg_quad_row, end_quad_row;
   int beg_quad_col, end_quad_col;
   int r, c, fail = 0;

   pqf[c0 + r0*num_cols] = loc_mask;

   if (0 != pt->pqf_flag_neighbor (pt, img, col_hw, row_hw, loc_mask, set_mask))
     return -1;

   /* quadrant boundaries are special */
   if (c0 < num_cols/2)
     {
        beg_quad_col = 0;
        end_quad_col = num_cols/2;
     }
   else
     {
        beg_quad_col = num_cols/2;
        end_quad_col = num_cols;
     }

   if (r0 < num_rows/2)
     {
        beg_quad_row = 0;
        end_quad_row = num_rows/2;
     }
   else
     {
        beg_quad_row = num_rows/2;
        end_quad_row = num_rows;
     }

   for (r = r0 - row_hw; (r < r0 + row_hw) && (r < end_quad_row); r++)
     {
        if (r < beg_quad_row)
          continue;

        for (c = c0 - col_hw; (c < c0 + col_hw) && (c < end_quad_col); c++)
          {
             if (c < beg_quad_col)
               continue;

             pqf_rc = pqf[c + r*num_cols];
             if (0 == (pqf_rc & set_mask))
               fail++;
#if 0
             fprintf (stderr, "*** %s: r0=%d c0=%d r=%d c=%d pqf=%d (expected %d) %s\n",
                      __func__, r0, c0, r, c, pqf_rc, set_mask,
                      (pqf_rc & set_mask) ? "OK" : "BAD");
#endif

             pqf[c + r*num_cols] = 0;
          }
     }

   pqf[c0 + r0*num_cols] = 0;

#if 0
   fprintf (stderr, "r0=%d c0=%d %s --------------------\n",
            r0, c0, fail ? "BAD" : "OK");
#endif

   return fail;
}

static int test_neighbor (Pixelqf_Type *pt)
{
   Image_Type *img;
   int hw, r, c, nr=10, nc=10;

   if (NULL == (img = image_new (nr, nc)))
     return -1;

   for (hw = 1; hw < 4; hw++)
     {
        image_set (img, 1.0, 0);

        for (r = 0; r < nr; r++)
          {
             for (c = 0; c < nc; c++)
               {
                  if (0 != test_neighbor1 (pt, img, r, c, hw))
                    {
                       image_free (img);
                       return -1;
                    }
               }
          }
     }

   image_free (img);
   return 0;
}

static int test_transient1 (Pixelqf_Type *pt,
                            const Image_Pqf_Bitmap_Type *bpixmap,
                            const Image_Type *img_ref,
                            Image_Type *img, int r0, int c0)
{
   Image_Pixel_Type pixel0;
   Image_Pqf_Bitmap_Type pqf0;
   int p0 = c0 + r0 * img->num_cols;
   int i, num_pixels = img->num_rows * img->num_cols;
   int fail = 0;

   pixel0 = img->pixels[p0];
   pqf0 = img->pixel_quality_flags[p0];

   img->pixels[p0] += 10.0;

   if (0 != pt->pqf_flag_transients (pt, bpixmap, img_ref, img))
     return -1;

   for (i = 0; i < num_pixels; i++)
     {
        Image_Pqf_Bitmap_Type pqf = img->pixel_quality_flags[i];
        if (i != p0)
          {
             if (0 != (pqf & IMAGE_PQF_TRANSIENT_PIXEL))
               {
                  fprintf (stderr, "*** %s: i=%d pixel=%g pqf=%d (expected 0)\n",
                           __func__, i, img->pixels[i], pqf);
                  fail++;
               }
          }
        else if (0 == (pqf & IMAGE_PQF_TRANSIENT_PIXEL))
          {
             fprintf (stderr, "*** %s: i=%d pixel=%g pqf=%d (expected %d)\n",
                      __func__, i, img->pixels[i], pqf, IMAGE_PQF_TRANSIENT_PIXEL);
             fail++;
          }
     }

   img->pixels[p0] = pixel0;
   img->pixel_quality_flags[p0] = pqf0;

   return fail;
}

static int test_transient (Pixelqf_Type *pt)
{
   Image_Type *img=NULL, *img_ref=NULL;
   Image_Pqf_Bitmap_Type *bpixmap = NULL;
   int r, c, nr=10, nc=10;
   size_t bpixmap_size = nr * nc * sizeof(*bpixmap);
   int status = -1;

   if ((NULL == (img = image_new (nr, nc)))
       || (NULL == (img_ref = image_new (nr, nc))))
     goto return_status;

   if (NULL == (bpixmap = (Image_Pqf_Bitmap_Type *)MALLOC (bpixmap_size)))
     goto return_status;

   memset ((char *)bpixmap, 0, bpixmap_size);
   image_set (img, 1.0, 0);
   image_set (img_ref, 1.0, 0);

   for (r = 0; r < nr; r++)
     {
        for (c = 0; c < nc; c++)
          {
             if (0 != test_transient1 (pt, bpixmap, img_ref, img, r, c))
               goto return_status;
          }
     }

   status = 0;
return_status:
   FREE(bpixmap);
   image_free (img);
   image_free (img_ref);
   return status;
}

static int perform_test (config_t *cfg)
{
   Pixelqf_Type *pt = NULL;
   typedef int test_fun_type (Pixelqf_Type *);
   test_fun_type *test_funs[] =
     {
        test_hotcold,
        test_neighbor,
        test_transient,
        NULL
     };
   test_fun_type **test_fun;

   if (NULL == (pt = pixelqf_init (cfg)))
     return -1;

   for (test_fun = test_funs; *test_fun != NULL; test_fun++)
     {
        if (0 != (*test_fun)(pt))
          {
             pt->pqf_delete (pt);
             return -1;
          }
     }

   pt->pqf_delete (pt);
   return 0;
}

int main (void)
{
   const char *config_file = "test_pixelqf.cfg";
   config_t cfg;
   int status = EXIT_FAILURE;

   config_init (&cfg);

   if (0 == config_read_file (&cfg, config_file))
     {
        fprintf (stderr, "%s: Reading %s:%d - %s\n",
                 __func__, config_error_file(&cfg),
                 config_error_line(&cfg), config_error_text(&cfg));
        goto return_status;
     }

   if (0 == perform_test (&cfg))
     status = EXIT_SUCCESS;

return_status:
   config_destroy (&cfg);
   if (status)
     fprintf (stderr, "*** ERROR: test_pixelqf failed\n");

   return status;
}

