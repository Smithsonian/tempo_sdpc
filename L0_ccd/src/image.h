#ifndef __IMAGE_INCLUDE__
# define __IMAGE_INCLUDE__ 1

#include <tio.h>

typedef float Image_Pixel_Type;
typedef unsigned short Image_Pqf_Bitmap_Type;

typedef struct
{
   Image_Pixel_Type *pixels;
   Image_Pqf_Bitmap_Type *pixel_quality_flags;
   int num_cols;
   int num_rows;
   int image_type;
}
Image_Type;

enum
{
   IMAGE_TYPE_UNKNOWN = 0,
   IMAGE_TYPE_PADDED = 1,      /* includes photo-active, oclocks, sdc, leading trailing */
   IMAGE_TYPE_ACTIVE = 2       /* photo-active pixels only */
};

/* fill value must be negative! */
#define IMAGE_PIXEL_FILL_VALUE  (-TIO_FILL_FLOAT)

/* pixel quality flag bit assignments */
#define IMAGE_PQF_MISSING_DATA          (1<<0)
#define IMAGE_PQF_BAD_PIXEL             (1<<1)
#define IMAGE_PQF_PROCESSING_ERROR      (1<<2)
#define IMAGE_PQF_TRANSIENT_PIXEL       (1<<3)
#define IMAGE_PQF_RTS_PIXEL             (1<<4)
#define IMAGE_PQF_SATURATED             (1<<5)
#define IMAGE_PQF_NOISE_UNDERFLOW       (1<<6)
#define IMAGE_PQF_DARK_CORR_ERROR       (1<<7)
#define IMAGE_PQF_OFFSET_CORR_ERROR     (1<<8)
#define IMAGE_PQF_SMEAR_CORR_ERROR      (1<<9)
#define IMAGE_PQF_STRAYLIGHT_CORR_ERROR (1<<10)
#define IMAGE_PQF_NONLINEAR_RANGE_ERROR (1<<11)
/* bits 0-11 are documented in the ICD and are consistent with TROPOMI,
 * (or at least they were consistent at some point). */
#define IMAGE_PQF_HOT_PIXEL             (1<<12)
#define IMAGE_PQF_COLD_PIXEL            (1<<13)

typedef struct
{
   int row_beg, row_end, row_step;
   int col_beg, col_end, col_step;
}
Image_Subset_Type;

extern void image_free (Image_Type *img);
extern Image_Type *image_new (int num_rows, int num_cols);
extern Image_Type *image_dup (const Image_Type *img);
extern int image_copy (const Image_Type *src, Image_Type *dest);
extern void image_sqrt (Image_Type *img);
extern void image_scale (Image_Type *img, double s);

extern Image_Subset_Type *image_new_subsets (int num_subsets);
extern void image_set_subset (Image_Subset_Type *s,
                              int row_beg, int row_end, int row_step,
                              int col_beg, int col_end, int col_step);

extern int image_get_type (const Image_Type *img);
extern void image_set_type (Image_Type *img, int image_type);

extern int image_write_raw (const Image_Type *img, const char *file);

#endif
