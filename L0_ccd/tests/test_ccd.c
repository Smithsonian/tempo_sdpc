#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <libconfig.h>
#include <image.h>
#include <ccd.h>

#include "util.h"

#define TEST_DIM_SERIAL   32
#define TEST_DIM_PARALLEL 32
#define TEST_NUM_SERIAL_ACTIVE    8
#define TEST_NUM_PARALLEL_ACTIVE 10
#define TEST_NUM_READOUT_BITS 14
#define TEST_NUM_COADD_BITS   20
#define TEST_FULL_WELL  240000

static void init_rownum_image (Image_Type *img)
{
   int p, s;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        if (0) fprintf (stderr, "%2d ", p);
        for (s = 0; s < img->num_cols; s++)
          {
             pixels[s] = p + 1.0;
             pqf[s] = 0;
             if (0) fprintf (stderr, "%4.1f ", pixels[s]);
          }
        if (0) fprintf (stderr, "\n");
     }
}

static void init_colnum_image (Image_Type *img)
{
   int p, s;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        if (0) fprintf (stderr, "%2d ", p);
        for (s = 0; s < img->num_cols; s++)
          {
             pixels[s] = s + 1.0;
             pqf[s] = 0;
             if (0) fprintf (stderr, "%4.1f ", pixels[s]);
          }
        if (0) fprintf (stderr, "\n");
     }
}

static int test_ccd_correct_coadd (CCD_Type *ccd)
{
   Image_Type *img = NULL;
   int num_rows = TEST_DIM_PARALLEL;
   int num_cols = TEST_DIM_SERIAL;
   int num_coadds = 21;
   int saturation_level_readout = (1 << TEST_NUM_READOUT_BITS) - 1;
   int saturation_level_coadded = (1 << TEST_NUM_COADD_BITS) - 1;
   Image_Pixel_Type sat_test_value;
   int status = -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   image_set (img, num_coadds, 0);
   if ((0 != ccd->ccd_correct_coadd (ccd, num_coadds, img))
       || (0 == image_test (img, 1.0, 0)))
     {
        fprintf (stderr, "*** coaddition correction\n");
        goto return_status;
     }

   image_set (img, IMAGE_PIXEL_FILL_VALUE, 0);
   if ((0 != ccd->ccd_correct_coadd (ccd, num_coadds, img))
       || (0 == image_test (img, IMAGE_PIXEL_FILL_VALUE, 0)))
     {
        fprintf (stderr, "*** coaddition correction (fill value)\n");
        goto return_status;
     }

   sat_test_value = saturation_level_readout * num_coadds + 1;
   image_set (img, sat_test_value, 0);
   if ((0 != ccd->ccd_correct_coadd (ccd, num_coadds, img))
       || (0 == image_test (img, sat_test_value/num_coadds, IMAGE_PQF_SATURATED)))
     {
        fprintf (stderr, "*** coaddition correction (readout sat)\n");
        goto return_status;
     }

   sat_test_value = saturation_level_coadded * num_coadds + 1;
   image_set (img, sat_test_value, 0);
   if ((0 != ccd->ccd_correct_coadd (ccd, num_coadds, img))
       || (0 == image_test (img, sat_test_value/num_coadds, IMAGE_PQF_SATURATED)))
     {
        fprintf (stderr, "*** coaddition correction (coadd sat)\n");
        goto return_status;
     }

   sat_test_value = (saturation_level_readout - 1) * num_coadds;
   image_set (img, sat_test_value, 0);
   if ((0 != ccd->ccd_correct_coadd (ccd, num_coadds, img))
       || (0 == image_test (img, sat_test_value/num_coadds, 0)))
     {
        fprintf (stderr, "*** coaddition correction (sub readout sat)\n");
        goto return_status;
     }

   status = 0;
return_status:
   image_free (img);
   return status;
}

static int test_ccd_correct_offset (CCD_Type *ccd)
{
   Image_Type *img = NULL;
   int num_rows = TEST_DIM_PARALLEL;
   int num_cols = TEST_DIM_SERIAL;
   int itest = num_cols/2;
   int status = -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   /* To test offset correction, generate an image in which
    * each pixel contains its row number.
    * The offset correction subtracts the mean trailing serial
    * pixels from all pixels in that row.
    * After applying this correction, its expected that all
    * of the corrected pixels will be zero.
    */
   init_rownum_image (img);
   if ((0 != ccd->ccd_correct_offset (ccd, img))
       || (0 == image_test (img, 0.0, 0)))
     {
        fprintf (stderr, "*** offset correction\n");
        goto return_status;
     }

   /* Verify that fill values remain unmodified */
   init_rownum_image (img);
   img->pixels[itest] = IMAGE_PIXEL_FILL_VALUE;
   if ((0 != ccd->ccd_correct_offset (ccd, img))
       || (img->pixels[itest] != IMAGE_PIXEL_FILL_VALUE))
     {
        fprintf (stderr, "*** offset correction (fill value)\n");
        goto return_status;
     }

   status = 0;
return_status:
   image_free (img);
   return status;
}

static int test_ccd_correct_nonlinearity (CCD_Type *ccd)
{
   Image_Type *img = NULL;
   int num_rows = TEST_DIM_PARALLEL;
   int num_cols = TEST_DIM_SERIAL;
   Image_Pixel_Type test_value = 1.0;
   Image_Pixel_Type expected_value = 21.0;
   int itest = num_cols/2;
   int status = -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   /* To test nonlinearity correction, generate a constant image,
    * and apply a nonlinearity correction defined so that it changes
    * the image to a known constant value.  In particular,
    * for f(x) = \sum_k k * x^k, we have f(1)=21 for k=1,2,...6
    * After applying this correction, its expected that all
    * of the corrected pixels will have the expected value.
    */

   image_set (img, test_value, 0);
   if ((0 != ccd->ccd_correct_nonlinearity (ccd, img))
       || (0 == image_test (img, expected_value, 0)))
     {
        fprintf (stderr, "*** nonlinearity correction\n");
        goto return_status;
     }

   /* Verify that fill values remain unmodified */
   image_set (img, test_value, 0);
   img->pixels[itest] = IMAGE_PIXEL_FILL_VALUE;
   if ((0 != ccd->ccd_correct_nonlinearity (ccd, img))
       || (img->pixels[itest] != IMAGE_PIXEL_FILL_VALUE))
     {
        fprintf (stderr, "*** nonlinearity correction (fill value)\n");
        goto return_status;
     }

   status = 0;
return_status:
   image_free (img);
   return status;
}

static int test_ccd_correct_gain (CCD_Type *ccd)
{
   Image_Type *img = NULL;
   int num_rows = TEST_DIM_PARALLEL;
   int num_cols = TEST_DIM_SERIAL;
   Image_Pixel_Type test_value = 3.0;
   Image_Pixel_Type expected_value = 1.0;
   int itest = num_cols/2;
   int status = -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   /* To test gain, generate a constant image,
    * and apply a gain correction defined so that it changes
    * the image to a known constant value.
    * The gain correction is: (test_value - offset)/gain
    * so, for offset=1, gain=2, test_value=3 gives an expected
    * value of 1.
    * After applying this correction, its expected that all
    * of the corrected pixels will have the expected value.
    */

   image_set (img, test_value, 0);
   if ((0 != ccd->ccd_correct_gain (ccd, img))
       || (0 == image_test (img, expected_value, 0)))
     {
        fprintf (stderr, "*** gain correction\n");
        goto return_status;
     }

   /* Verify that fill values remain unmodified */
   image_set (img, test_value, 0);
   img->pixels[itest] = IMAGE_PIXEL_FILL_VALUE;
   if ((0 != ccd->ccd_correct_gain (ccd, img))
       || (img->pixels[itest] != IMAGE_PIXEL_FILL_VALUE))
     {
        fprintf (stderr, "*** gain correction (fill value)\n");
        goto return_status;
     }

   /* Verify that gain-corrected pixels that exceed full-well
    * are marked as saturated
    */
   image_set (img, test_value, 0);
   img->pixels[itest] = 2*(TEST_FULL_WELL+1) + 1;
   if ((0 != ccd->ccd_correct_gain (ccd, img))
       || (img->pixel_quality_flags [itest] != IMAGE_PQF_SATURATED))
     {
        fprintf (stderr, "*** gain correction (full well)\n");
        goto return_status;
     }

   status = 0;
return_status:
   image_free (img);
   return status;
}

static int test_ccd_correct_smear_method (CCD_Type *ccd,
                                          const char *method)
{
   Image_Type *img = NULL;
   Image_Type *img_active = NULL;
   int num_rows = TEST_DIM_PARALLEL;
   int num_cols = TEST_DIM_SERIAL;
   double smear_fraction = 1.0;
   Image_Pixel_Type expected_value = 0.0;
   int itest = num_cols*(num_rows/2);
   int status = -1;

   if (0 != __ccd_set_smear_corr_method (ccd, method))
     return -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   /* To test smear correction, generate an image in which
    * each pixel contains its column number.
    * The smear correction subtracts a column-dependent value
    * from all pixels in each column.
    * After applying this correction, its expected that all
    * of the corrected pixels in the _active_ image
    * will have the expected value. In this case:
    *   method=timing, smear_fraction=1.0 => expected_value=0.0
    *
    * Note that smear_fraction is used only by the 'timing' method
    */
   init_colnum_image (img);
   image_set_type (img, IMAGE_TYPE_PADDED);
   if ((0 != ccd->ccd_correct_smear (ccd, &smear_fraction, img))
       || (NULL == (img_active = ccd->ccd_copy_active_pixels (ccd, img)))
       || (0 == image_test (img_active, expected_value, 0)))
     {
        fprintf (stderr, "*** smear correction\n");
        goto return_status;
     }

   /* Verify that fill values remain unmodified */
   init_colnum_image (img);
   img->pixels[itest] = IMAGE_PIXEL_FILL_VALUE;
   if ((0 != ccd->ccd_correct_smear (ccd, &smear_fraction, img))
       || (img->pixels[itest] != IMAGE_PIXEL_FILL_VALUE))
     {
        fprintf (stderr, "*** smear correction (fill value)\n");
        goto return_status;
     }

   status = 0;
return_status:
   image_free (img);
   image_free (img_active);
   return status;
}

static int test_ccd_correct_smear (CCD_Type *ccd)
{
   const char *method_names[] = {"oclocks", "timing", NULL};
   const char **method;
   int status;

   for (method = method_names; *method != NULL; method++)
     {
        if (0 != (status = test_ccd_correct_smear_method (ccd, *method)))
          return status;
     }

   return 0;
}

static int test_ccd_mean_sdc (CCD_Type *ccd)
{
   Image_Type *img = NULL;
   int num_rows = TEST_DIM_PARALLEL;
   int num_cols = TEST_DIM_SERIAL;
   float mean_sdc[4];
   float sdc_row_value[4] = {1.0, 1.0, TEST_DIM_PARALLEL, TEST_DIM_PARALLEL};
   float expected_mean_sdc[4];
   int i, status = -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   init_rownum_image (img);
   for (i = 0; i < 4; i++)
     {
        expected_mean_sdc[i] = sdc_row_value[i] / TEST_NUM_PARALLEL_ACTIVE;
     }

   if (0 != ccd->ccd_mean_storage_region_dark (ccd, img, mean_sdc))
     {
        fprintf (stderr, "*** mean sdc\n");
        goto return_status;
     }

   for (i = 0; i < 4; i++)
     {
        if (fabs(mean_sdc[i] - expected_mean_sdc[i])
            > FLT_EPSILON * expected_mean_sdc[i])
          {
             fprintf (stderr, "*** mean_sdc[%d] = %7.4f (expected %7.4f)\n",
                      i, mean_sdc[i], expected_mean_sdc[i]);
             goto return_status;
          }
     }

   status = 0;
return_status:
   image_free (img);
   return status;
}

static int test_ccd_active_image_dims (CCD_Type *ccd)
{
   int num_parallel_active_full, num_serial_active_full;

   ccd->ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);
   if ((num_parallel_active_full != 2*TEST_NUM_PARALLEL_ACTIVE)
       || (num_serial_active_full != 2*TEST_NUM_SERIAL_ACTIVE))
     {
        fprintf (stderr, "*** active image dims\n");
        return -1;
     }

   return 0;
}

static int perform_test (int argc, char **argv, config_t *cfg)
{
   CCD_Type *ccd = NULL;
   TIO_Meta_Type *meta = NULL;
   typedef int test_fun_type (CCD_Type *);
   test_fun_type *test_funs[] =
     {
        test_ccd_correct_coadd,
        test_ccd_correct_offset,
        test_ccd_correct_nonlinearity,
        test_ccd_correct_gain,
        test_ccd_correct_smear,
        test_ccd_mean_sdc,
        test_ccd_active_image_dims,
        NULL
     };
   test_fun_type **fun;
   int status = -1;

   (void) argc; (void) argv;

   if (NULL == (meta = tio_meta_open ()))
     return -1;

   if (NULL == (ccd = ccd_init (cfg, meta)))
     goto return_status;

   for (fun = test_funs; *fun != NULL; fun++)
     {
        if (0 != (*fun)(ccd))
          goto return_status;
     }

   status = 0;
return_status:
   tio_meta_close (meta);
   if (ccd) ccd->ccd_delete (ccd);
   return status;
}

int main (int argc, char **argv)
{
   const char *config_file = "test_ccd.cfg";
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

   if (0 == perform_test (argc, argv, &cfg))
     status = EXIT_SUCCESS;

return_status:
   config_destroy (&cfg);
   if (status)
     fprintf (stderr, "*** ERROR: test_ccd failed\n");

   return status;
}
