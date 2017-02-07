#include <stdio.h>
#include <stdlib.h>

#include <image.h>
#include "util.h"

void image_set (Image_Type *img, Image_Pixel_Type c, Image_Pqf_Bitmap_Type b)
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

int image_test (Image_Type *img, Image_Pixel_Type c, Image_Pqf_Bitmap_Type b)
{
   int p, s;

   for (p = 0; p < img->num_rows; p++)
     {
        Image_Pixel_Type *pixels = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags + p * img->num_cols;
        for (s = 0; s < img->num_cols; s++)
          {
             if ((pixels[s] != c) || (pqf[s] != b))
               {
                  fprintf (stderr,
                           "*** Error: p=%d s=%d pixels=%f pqf=%d (expected %f,%d)\n",
                           p, s, pixels[s], pqf[s], c, b);
                  return 0;
               }
          }
     }

   return 1;
}

