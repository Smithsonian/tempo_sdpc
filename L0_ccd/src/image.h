#ifndef __IMAGE_INCLUDE__
#define __IMAGE_INCLUDE__ 1
/** @file image.h
 *  @brief Basic image data object
 *
 * This module implements generic operations on a simple
 * image object.
 */

#include <tio.h>

/** data type for image pixels */
typedef float Image_Pixel_Type;

/** data type for image pixel quality flags */
typedef unsigned short Image_Pqf_Bitmap_Type;

/** @brief Struct defining an image object */
typedef struct
{
   Image_Pixel_Type *pixels;
   Image_Pqf_Bitmap_Type *pixel_quality_flags;
   int num_cols;
   int num_rows;
   int image_type;
}
Image_Type;

/** @brief Enum identifying supported image types
 *  @anchor image_type_enum_list
 */
enum
{
   IMAGE_TYPE_UNKNOWN = 0,
   IMAGE_TYPE_PADDED = 1,      /**< includes photo-active, oclocks, sdc, leading trailing */
   IMAGE_TYPE_ACTIVE = 2       /**< photo-active pixels only */
};

/** image pixel fill value must be negative! */
#define IMAGE_PIXEL_FILL_VALUE  (-TIO_FILL_FLOAT)

/** @brief Pixel quality flag bit assignments
 *  @anchor pixel_quality_flag_bits
 */
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

/** @brief Struct defining an image subset.
 *
 * An image subset is the set of pixels with coordinates
 * @verbatim
    row = row_beg + n * row_step, n = 0,1,2, ...
    col = col_beg + m * col_step, m = 0,1,2, ...
  such that
    row_beg <= row < row_end
    col_beg <= col < col_end
  @endverbatim
 */
typedef struct
{
   int row_beg, row_end, row_step;
   int col_beg, col_end, col_step;
}
Image_Subset_Type;

/** Free an image object
 * @param   img    A pointer to an Image_Type object
 */
extern void image_free (Image_Type *img);

/** Create a new image object
 * @param   num_rows  Dimension of the slowest varying array index
 * @param   num_cols  Dimension of the quickest varying array index
 * @return a non-NULL pointer to a new Image_Type object on success,
 * NULL on error.
 *
 * Create a new \a Image_Type object with the specified dimensions.
 */
extern Image_Type *image_new (int num_rows, int num_cols);

/** Duplicate an image object
 * @param   img     A non-NULL pointer to an existing Image_Type object.
 * @return a non-NULL pointer to a duplicate Image_Type object on success,
 * NULL on error.
 *
 * The duplicate Image_Type object contains a copy of the pixel values
 * and pixel quality flags of the input Image_Type object.  The input
 * Image_Type object is not modified.
 */
extern Image_Type *image_dup (const Image_Type *img);

/** Copy data fields from one image object to another with the same dimensions
 * @param   src    A non-NULL pointer to an existing Image_Type object
 * @param   dest   A non-NULL pointer to an existing Image_Type object
 * @return 0 on success, non-zero on error
 *
 * The data fields of the \a src object are copied to the \a dest object.
 * Both \a src and \a dest objects must have the same dimensions.
 */
extern int image_copy (const Image_Type *src, Image_Type *dest);

/** Compute the square root of an image
 * @param   img    A non-NULL pointer to an existing Image_Type object
 *
 * Each pixel value in \a img is replaced by its square root.
 * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
 */
extern void image_sqrt (Image_Type *img);

/** Multiply an image by a scalar constant
 * @param   img    A non-NULL pointer to an existing Image_Type object
 * @param   s      A scalar constant
 *
 * Each pixel value in \a img is multipled by the scalar constant \a s.
 * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
 */
extern void image_scale (Image_Type *img, double s);

/** Allocate an array of image subset structs
 * @param num_subsets   The number of Image_Subset_Type struct to allocate
 * @return a non-NULL pointer to an array of Image_Subset_Type structs on success,
 * or NULL on error.
 */
extern Image_Subset_Type *image_new_subsets (int num_subsets);

/** Define an image subset
 * @param  s       A non-NULL pointer to an Image_Subset_Type struct
 * @param row_beg  The beginning row of the image subset.
 * @param row_end  The row that lies one row_step beyond the last row included
 *                 in the subset.
 * @param row_step The interval between rows included in the subset.
 * @param col_beg  The beginning column of the image subset.
 * @param col_end  The column that lies one col_step beyond the last
 *                 column included in the subset.
 * @param col_step The interval between columns included in the subset.
 *
 * An image subset is the set of all pixels with array indices:
 * @verbatim
      row = row_beg + n * row_step, n=0,1,2,...
      col = col_beg + m * col_step, m=0,1,2,...
    such that
      row_beg <= row < row_end
      col_beg <= col < col_end
   @endverbatim
 */
extern void image_set_subset (Image_Subset_Type *s,
                              int row_beg, int row_end, int row_step,
                              int col_beg, int col_end, int col_step);

/** Query the type of an image
 * @param   img    A non-NULL pointer to an existing Image_Type object
 * @return the \ref image_type_enum_list "image type" enum value
 */
extern int image_get_type (const Image_Type *img);

/** Set the type of an image
 * @param   img         A non-NULL pointer to an existing Image_Type object
 * @param   image_type  The \ref image_type_enum_list "image type" enum value
 */
extern void image_set_type (Image_Type *img, int image_type);

/** Write image data fields to raw binary-formatted files
 * @param   img         A non-NULL pointer to an existing Image_Type object
 * @param   file        A filename prefix.
 *
 * <em>This function is intended to be used only for debugging.</em>
 */
extern int image_write_raw (const Image_Type *img, const char *file);

#endif
