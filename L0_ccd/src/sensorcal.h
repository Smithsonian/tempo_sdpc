#ifndef __SENSORCAL_INCLUDE__
#define __SENSORCAL_INCLUDE__ 1
/** @file sensorcal.h
 *  @brief Sensor calibration interface
 */

#include <libconfig.h>
#include "image.h"

enum
{
   TEMPO_BAND_UV,
   TEMPO_BAND_VIS
};

typedef struct Calibration_Type Calibration_Type;

/** @brief Struct member functions perform various sensor calibrations
 */
struct Calibration_Type
{
   /** Free a Calibration_Type object
    * @param cal  non-NULL pointer to a Calibration_Type object
    */
   void (*cal_delete)(Calibration_Type *);

   /** Apply radiometric calibration
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param img  non-NULL pointer to an uncalibrated image
    * @return 0 on success, non-zero on error
    *
    * Applying the radiometric calibration modifies the image in place.
    * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
    */
   int (*cal_apply_radcal_coeffs)(const Calibration_Type *, Image_Type *);

   /** Apply BTDF calibration
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param is_reference   0 for the working diffuser, non-zero for the reference diffuser
    * @param solar_phi    azimuthal angle coordinate of the solar position vector
    * @param solar_theta  polar angle coordinate of the solar position vector
    * @param img  non-NULL pointer to an uncalibrated image
    * @param img_diag    pointer to space for the BTDF values (diagnostic - NULL is ok)
    * @return 0 on success, non-zero on error
    *
    * The BTDF is the bidirectional transmission distribution function
    * describing the relevant diffuser.  The BTDF factor depends on the direction
    * of the incident solar irradiance.  Applying the BTDF factor modifies
    * the image in place. Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
    */
   int (*cal_apply_btdf)(const Calibration_Type *, int, double, double,
                         Image_Type *, Image_Type *);

   /** Apply correction for diffuser polarization
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param solar_phi    azimuthal angle coordinate of the solar position vector
    * @param solar_theta  polar angle coordinate of the solar position vector
    * @param img  non-NULL pointer to an uncorrected irradiance image
    * @param img_diag    pointer to space for the BTDF values (diagnostic - NULL is ok)
    * @return 0 on success, non-zero on error
    *
    * The input image is modified in place.
    * Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
    */
   int (*cal_apply_diffuser_polcorr)(const Calibration_Type *, double, double,
                                     Image_Type *, Image_Type *);

   int (*cal_straylight_correction)(const Calibration_Type *, Image_Type *);

   /** Define a nominal wavelength grid
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param band_id  integer band index (TEMPO_BAND_UV | TEMPO_BAND_VIS)
    * @param img_waves non-NULL pointer to the nominal wavelength grid
    * @return 0 on success, non-zero on error
    *
    * The initial Level 1 radiance and irradiance files contain a nominal
    * wavelength grid that's used to initialize wavelength calibration
    * (which happens later), and also, for radiances, to provide INR
    * with a nominal wavelength grid.
    */
   int (*cal_nominal_wavelength_grid)(const Calibration_Type *, int, double *);

#ifdef SENSORCAL_PRIVATE_DATA
   SENSORCAL_PRIVATE_DATA
#endif
};

/** Initialize a Calibration_Type object
 * @param cfg  non-NULL pointer to a config_t object associated with an open configuration file
 * @param meta non-NULL pointer to a \ref TIO_Meta_Type object
 * @return non-NULL pointer to a Calibration_Type object on success, NULL on error
 */
extern Calibration_Type *sensorcal_init (config_t *cfg, TIO_Meta_Type *meta, const char *irr_file, int exposure_type);

typedef struct
{
   double *img;
   double *img_err;
   double *wave;
   Image_Pqf_Bitmap_Type *pqf;
   char *name;
   int num_xtrack;
   int num_channels;
}
Spectral_Data_Type;
/**< Spectral_Data_Type arrays are ordered so that wavelength
 * is the fastest varying index to facilitate wavelength
 * calibration and spectral analysis.
 */

extern void sdt_free (Spectral_Data_Type *sdt);

extern Spectral_Data_Type *
sdt_extract_band (const Calibration_Type *cal, int band_id,
                  const Image_Type *img, const Image_Type *img_err);

#endif
