#ifndef __TEMPO_PIXELQF_INCLUDE__
#define __TEMPO_PIXELQF_INCLUDE__
/** @file pixelqf.h
 *  @brief Algorithms for defining selected pixel quality flags
 */

#include <libconfig.h>
#include "image.h"

typedef struct
{
   int num_hot_pixels[4];  /* per quadrant: A,B,C,D */
   int num_cold_pixels[4];
   float mean_dark_current[4];  /* electrons/sec */
}
Dark_Trend_Type;

typedef struct Pixelqf_Type Pixelqf_Type;

/** @brief Struct providing functions to set selected pixel quality flags */
struct Pixelqf_Type
{
   /** Free a Pixelqf_Type object
    * @param  pt  non-NULL pointer to a Pixelqf_Type object
    */
   void (*pqf_delete)(Pixelqf_Type *);

   /** Flag hot and cold pixels
    * @param pt  non-NULL pointer to a Pixelqf_Type object
    * @param img  non-NULL pointer to an Image_Type object
    * @param dark_trend  non-NULL pointer to a Dark_Trend_Type object
    * @return 0 on success, non-zero on error
    *
    * Hot pixels exceed the mean pixel value in a quadrant by more than
    * \a hot_thresh standard deviations.  Cold pixels fall
    * below the mean pixel value in a quadrant by more than \a cold_thresh
    * standard deviations.  The values of \a hot_thresh and
    * \a cold_thresh are obtained from a configuration file.
    * Hot pixels are flagged by setting the pixel quality flag
    * bit associated with \a IMAGE_PQF_HOT_PIXEL.
    * Cold pixels are flagged by setting the pixel quality flag
    * bit associated with \a IMAGE_PQF_COLD_PIXEL.
    */
   int (*pqf_flag_hotcold)(const Pixelqf_Type *, Image_Type *, Dark_Trend_Type *);

   /** Flag pixels adjacent to pixels that match a specified bitmask
    * @param pt  non-NULL pointer to a Pixelqf_Type object
    * @param img  non-NULL pointer to an Image_Type object
    * @param hw_serial  neighborhood half-width [pixels] in serial readout direction
    * @param hw_parallel neighborhood half-width [pixels] in parallel readout direction
    * @param loc_mask   mask used to select image pixels
    * @param set_mask   mask used to set bits in neighboring pixels
    * @return 0 on success, non-zero on error
    *
    * Each pixel in the provided image has a corresponding set of pixel
    * quality flags, \a pqf.  Pixels are selected by identifying those
    * for which pqf|loc_mask is non-zero.  The quality flags for neighbors
    * of each selected pixel are modified by setting the bits in \a set_mask.
    */
   int (*pqf_flag_neighbor)(const Pixelqf_Type *, Image_Type *, int, int,
                            Image_Pqf_Bitmap_Type, Image_Pqf_Bitmap_Type);

   /** Flag transient pixels by comparing with a reference image
    * @param pt       non-NULL pointer to a Pixelqf_Type object
    * @param bpixmap  non-NULL pointer to a bad pixel map
    * @param img_ref  non-NULL pointer to a reference image
    * @param img      non-NULL pointer to an Image_Type object
    *
    * This algorithm examines pixels that have positive values
    * in both the input image and the reference image and that
    * do not coincide with known bad pixels.  The reference image
    * is normally constructed by averaging exposure records acquired
    * immediately before and immediately after the image of interest.
    * When
    * @code
    *    pixel(img)/pixel(img_ref) > 1 + threshold
    * @endcode
    * the corresponding image pixel may be flagged as a transient
    * if the contrast with adjacent pixels is also sufficiently high.
    * Transient pixels are flagged by setting the bit associated
    * with IMAGE_PQF_TRANSIENT_PIXEL.
    */
   int (*pqf_flag_transients)(const Pixelqf_Type *,
                              const Image_Pqf_Bitmap_Type *,
                              const Image_Type *, Image_Type *);

#ifdef PIXELQF_TYPE_PRIVATE_DATA
   PIXELQF_TYPE_PRIVATE_DATA
#endif
};

/** Initialize a Pixelqf_Type object
 * @param cfg   Pointer to a config_t struct associated with an open configuration file
 * @return non-NULL pointer to a Pixelqf_Type struct on success, NULL on error
 */
extern Pixelqf_Type *pixelqf_init (config_t *cfg);

#endif
