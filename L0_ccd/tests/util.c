#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <image.h>
#include "util.h"

void util_image_set (Image_Type *img, Image_Pixel_Type c, Image_Pqf_Bitmap_Type b)
{
   int p, s;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = 0; s < img->num_cols; s++)
          {
             pixels[s] = c;
             pqf[s] = b;
          }
     }
}

int util_image_test (Image_Type *img, Image_Pixel_Type c, Image_Pqf_Bitmap_Type b)
{
   int p, s;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = 0; s < img->num_cols; s++)
          {
             if ((pqf[s] != b)
                 || ((pixels[s] != c)
                     && (fabs(pixels[s]-c) > FLT_EPSILON*(0.5*(pixels[s]+c)))))
               {
                  fprintf (stderr,
                           "*** Error: p=%d s=%d pixels=%g pqf=%d (expected %g,%d)\n",
                           p, s, pixels[s], pqf[s], c, b);
                  return 0;
               }
          }
     }

   return 1;
}

