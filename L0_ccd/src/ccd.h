#ifndef __TEMPO_CCD_INCLUDE__
#define __TEMPO_CCD_INCLUDE__
/** @file ccd.h
 *  @brief Basic CCD data conversion
 *
 * This module's primary function is to facilitate converting
 * raw pixel values to a number of electrons in each pixel.
 */

#include <libconfig.h>
#include "image.h"

typedef struct CCD_Type CCD_Type;

/** @brief Struct member functions facilitate basic CCD data conversion and correction.
 *
 * The functions that perform corrections are intended to be called
 * in the following order:
 *   -# ccd_correct_coadd
 *   -# ccd_decide_phase_change
 *   -# ccd_correct_offset
 *   -# ccd_correct_nonlinearity
 *   -# ccd_correct_gain
 *   -# ccd_correct_smear
 */
struct CCD_Type
{
   /** A destructor.
    * @param ccd         A non-null pointer to a CCD_Type object
    *
    * Call this function to release resources associated with the
    * CCD_Type object.
    *
    * For example:
    * @code
    *   CCD_Type *ccd = ccd_init (cfg);   // initialize CCD_Type object
    *    ...
    *   ccd->ccd_delete (ccd);            // free CCD_Type object
    * @endcode
    */
   void (*ccd_delete)(CCD_Type *);

   /** Co-addition correction
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param num_coadds  Number of co-added frames
    * @param img         A non-null pointer to an Image_Type object
    * @return 0 on success, non-zero on error
    *
    * The value of each pixel in the input image is divided by \a num_coadds.
    * The input image is modified in place.
    * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
    * Saturated pixels are defined as those that overflow 20 bits
    * before correction, or that overflow 14 bits after correction.
    * Saturated pixels are flagged by setting the \a IMAGE_PQF_SATURATED
    * bit in the pixel quality flag.
    */
   int (*ccd_correct_coadd)(const CCD_Type *, int, Image_Type *);

   int (*ccd_configure_using_octant_phase)(CCD_Type *, const Image_Type *);
   int (*ccd_correct_crosstalk)(const CCD_Type *, Image_Type *);

   /** Offset correction
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param img         A non-null pointer to an Image_Type object
    * @return 0 on success, non-zero on error
    *
    * Each pixel value includes an additive offset.  The value of the
    * offset for each octant is derived by computing the average of a
    * subset of the trailing serial readout pixels from that octant.
    * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are ignored and are
    * not modified.  For each octant, the offset computed for each row
    * is subtracted from the pixels in that row.  The input image is
    * modified in place.
    */
   int (*ccd_correct_offset)(const CCD_Type *, Image_Type *);

   /** Nonlinearity correction
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param img         A non-null pointer to an Image_Type object
    * @return 0 on success, non-zero on error
    *
    * Nonlinearity correction is a polynomial transformation
    * of the pixel values in the input image.  The image is modified
    * in place. Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not
    * modified. Each octant may have a separate nonlinearity correction.
    */
   int (*ccd_correct_nonlinearity)(const CCD_Type *, Image_Type *);

   /** Gain correction
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param img         A non-null pointer to an Image_Type object
    * @param fpa_temp    Focal Plane Assembly temperature
    * @param fpe_temp    Focal Plane Electronics temperature
    * @return 0 on success, non-zero on error
    *
    * Gain correction is a linear transformation of the pixels in the
    * input image.  Each octant may have a separate gain correction.
    * After this transformation, each pixel value represents
    * the number of electrons in the potential well of the pixel.
    * The image is modified in place. Pixels containing
    * \a IMAGE_PIXEL_FILL_VALUE are not modified.
    * Corrected pixels with values that exceed the expected full-well
    * number of electrons in a pixel are flagged by setting the
    * \a IMAGE_PQF_SATURATED bit in the pixel quality flag.
    */
   int (*ccd_correct_gain)(const CCD_Type *, Image_Type *, float, float);

   /** Smear correction
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param client_data A pointer to optional, method-dependent client data
    * @param img         A non-null pointer to an Image_Type object
    * @return 0 on success, non-zero on error
    *
    * A small amount of charge is subtracted from each image pixel to
    * account for photoelectrons generated during the readout interval
    * between exposures.  The correction is applied columnwise.
    * In each column, \f$k\f$, the size of the correction, \f$c_k\f$,
    * depends on which smear correction algorithm is selected.
    * The \a oclocks method derives
    * the correction for each column by averaging the parallel overclock
    * pixel values in that column.  The \a timing method computes the
    * correction, \f$c_k\f$, using the mean
    * value of each column \f$(\bar{p_k})\f$,
    * the frame transfer (readout) time \f$(t_r)\f$,
    * and the integration (exposure) time \f$(t_i)\f$:
    * \f[
    *    c_k = \frac{t_r}{t_r+t_i}\bar{p_k}.
    * \f]
    * The image is modified in place. Pixels containing
    * \a IMAGE_PIXEL_FILL_VALUE are ignored, and are not modified.
    */
   int (*ccd_correct_smear)(const CCD_Type *, const void *, Image_Type *);

   /** Compute the mean, per-pixel storage region dark count
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param img         A non-null pointer to an Image_Type object
    * @param num_dg_rows First storage region dark row summed.
    * @param num_tg_rows Number of storage region rows summed.
    * @param mean_sdc    Mean, per-pixel storage region dark count in each quadrant.
    * @return 0 on success, non-zero on error
    *
    * One row in each quadrant contains the total dark counts generated in
    * (a subset of) the storage-region pixels during the readout interval.
    * The mean per-pixel storage-region dark count accumulated during the
    * readout interval is obtained by dividing the pixel value by the
    * number of pixels it represents.
    * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are ignored.
    */
   int (*ccd_mean_storage_region_dark)(const CCD_Type *, const Image_Type *,
                                       int, int, float *mean_sdc);

   /** Retrieve the dimensions of the photo-active CCD image
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param num_parallel_active_full  Pointer to the number of active pixels
    *                                  along the parallel readout direction
    * @param num_serial_active_full    Pointer to the number of active pixels
    *                                  along the serial readout direction.
    *
    * The return values give the dimensions of the entire photo-active region
    * of the full CCD (including both UV and Vis bands).
    */
   void (*ccd_active_image_dims)(const CCD_Type *, int *, int *);

   /** Create a new image containing only the photo-active pixels
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param img         A non-null pointer to an Image_Type object
    * @return a non-NULL pointer to a new Image_Type object on success,
    * NULL on error.
    *
    * Create a new \a Image_Type object that contains only the photo-active
    * pixel values and the corresponding pixel quality flags.  The pixels
    * associated with leading/trailing serial readout, the parallel overclocks,
    * and the storage region dark current are omitted from the newly created
    * image.  The input Image_Type object is not modified.
    */
   Image_Type *(*ccd_copy_active_pixels)(const CCD_Type *, const Image_Type *);

   int (*ccd_apply_pixel_quality_flags)(const CCD_Type *, Image_Type *img,
                                        const Image_Pqf_Bitmap_Type *flags,
                                        int num_rows, int num_cols);

   /** Add noise contribution from CTE and quantization
    * @param ccd         A non-null pointer to the CCD_Type object
    * @param sdc         A non-null pointer to a float array of size 4
    *                    containing the mean per-pixel storage region
    *                    dark count in each quadrant.
    * @param img         A non-null pointer to an Image_Type object
    *                    initialized with a copy of the active-region
    *                    image, smear-corrected, with pixels in units of
    *                    electrons.
    *
    * The expression for CTE noise is an empirical fit to OMPS data
    * and will eventually be replaced by a TEMPO-specific expression.
    * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are ignored and are
    * not modified.
    */
   int (*ccd_update_noisesq)(const CCD_Type *, const float *, Image_Type *);

   int (*ccd_correct_prnu)(const CCD_Type *, Image_Type *);

#ifdef CCD_TYPE_PRIVATE_DATA
   CCD_TYPE_PRIVATE_DATA
#endif
};

/** Initialize a CCD_Type object
 * @param   cfg   A non-NULL pointer to a \ref config_t object associated
 *                with an open configuration file.
 * @param  meta   A non-NULL pointer to a \ref TIO_Meta_Type object
 * @return a non-NULL pointer to a CCD_Type object on success, NULL on error.
 *
 * User-configurable parameter values are provided in a configuration file.
 * The values of these parameters are used to configure the CCD_Type object.
 * that is returned.
 *
 * When the CCD_Type object is no longer needed, free the associated resources
 * by passing the CCD_Type pointer to the destructor function, \ref ccd_delete
 */
extern CCD_Type *ccd_init (config_t *cfg, TIO_Meta_Type *meta);

/** Select the smear correction method
 * @param ccd         A non-null pointer to the CCD_Type object
 * @param name        Pointer to the name of a valid smear correction method
 * @return 0 on success, non-zero on error
 *
 * <em> This function has global scope only to facilitate testing.
 * It is not intended for use in any other context. </em>
 *
 * It facilitates testing by enabling selection of a
 * smear correction method without reading a configuration file.
 */
extern int __ccd_set_smear_corr_method (CCD_Type *ccd, const char *name);

typedef struct CCD_Linearity_Type CCD_Linearity_Type;
typedef struct CCD_Select_Type CCD_Select_Type;

struct CCD_Linearity_Type
{
   void (*clt_delete)(CCD_Linearity_Type *);
   int (*clt_correct_coadd)(const CCD_Linearity_Type *, int, Image_Type *);
   int (*clt_correct_offset)(const CCD_Linearity_Type *, Image_Type *);
   int (*clt_trimmed_sample_mean)(const CCD_Linearity_Type *, const Image_Type *,
                                  CCD_Select_Type *, double *);
#ifdef CCD_LINEARITY_TYPE_PRIVATE_DATA
   CCD_LINEARITY_TYPE_PRIVATE_DATA
#endif
};

extern CCD_Select_Type *clt_select_alloc (size_t num_pixels);
extern void clt_select_free (CCD_Select_Type *sel);

extern CCD_Linearity_Type *ccd_linearity_init (void);

#endif
