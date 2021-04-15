#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <config.h>
#include <image.h>

#include "util.h"

static int test_write_raw (const Image_Type *img, const char *prefix)
{
   Image_Type *in = NULL;
   char buf[256];
   FILE *fp = NULL;
   int status = -1;
   size_t i, num_pixels, num_read;

   if (0 != image_write_raw (img, prefix))
     return -1;

   if (NULL == (in = image_new (img->num_rows, img->num_cols)))
     return -1;

   num_pixels = img->num_rows * img->num_cols;

   sprintf (buf, "%s_pixels.arr", prefix);
   if (NULL == (fp = fopen (buf, "r")))
     {
        fprintf (stderr, "*** %s: open failed: %s\n", __func__, buf);
        goto return_status;
     }

   num_read = fread (in->pixels, sizeof(*in->pixels), num_pixels, fp);
   (void) fclose (fp);

   if (num_read != num_pixels)
     {
        fprintf (stderr, "*** %s: read failed: %s\n", __func__, buf);
        goto return_status;
     }

   sprintf (buf, "%s_pqf.arr", prefix);
   if (NULL == (fp = fopen (buf, "r")))
     {
        fprintf (stderr, "*** %s: open failed: %s\n", __func__, buf);
        goto return_status;
     }

   num_read = fread (in->pixel_quality_flags,
                     sizeof(*in->pixel_quality_flags), num_pixels, fp);
   (void) fclose(fp);

   if (num_read != num_pixels)
     {
        fprintf (stderr, "*** %s: read failed: %s\n", __func__, buf);
        goto return_status;
     }

   for (i = 0; i < num_pixels; i++)
     {
        if ((in->pixels[i] != img->pixels[i])
            || (in->pixel_quality_flags[i] != img->pixel_quality_flags[i]))
          {
             fprintf (stderr, "*** %s: input image differs!\n", __func__);
             goto return_status;
          }
     }

   status = 0;
return_status:
   image_free (in);
   return status;
}

static int perform_test (int num_rows, int num_cols)
{
   Image_Type *img = NULL;
   Image_Type *dup = NULL;
   Image_Pixel_Type pixel_value = 2.0;
   Image_Pqf_Bitmap_Type pqf_value = 2;
   int status = -1;

   if (NULL == (img = image_new (num_rows, num_cols)))
     return -1;

   util_image_set (img, pixel_value, pqf_value);
   if (NULL == (dup = image_dup (img)))
     goto return_status;

   if (0 == util_image_test (dup, pixel_value, pqf_value))
     goto return_status;

   image_scale (dup, pixel_value);

   if (0 == util_image_test (dup, pixel_value*pixel_value, pqf_value))
     goto return_status;

   image_sqrt (dup);

   if (0 == util_image_test (dup, pixel_value, pqf_value))
     goto return_status;

   if (0 != test_write_raw (dup, "img_dup"))
     goto return_status;

   status = 0;
return_status:
   image_free (img);
   image_free (dup);
   return status;
}

int main (void)
{
   int num_rows = 10, num_cols = 5;
   int status = EXIT_FAILURE;

   if (0 == perform_test (num_rows, num_cols))
     status = EXIT_SUCCESS;

   if (status)
     fprintf (stderr, "*** ERROR: test_image failed\n");

   return status;
}

