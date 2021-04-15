#ifndef __TEST_UTIL_H__
#define __TEST_UTIL_H__ 1

#include <image.h>

extern void util_image_set (Image_Type *img, Image_Pixel_Type c, Image_Pqf_Bitmap_Type b);
extern int util_image_test (Image_Type *img, Image_Pixel_Type c, Image_Pqf_Bitmap_Type b);

#endif
