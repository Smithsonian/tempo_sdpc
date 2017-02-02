#ifndef __SENSORCAL_INCLUDE__
#define __SENSORCAL_INCLUDE__ 1
/** @file sensorcal.h
 *  @brief Sensor calibration interface
 */

#include <libconfig.h>
#include "image.h"

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
   int (*cal_apply_rcoeffs)(const Calibration_Type *, Image_Type *);

   /** Apply BTDF calibration
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param solar_phi    azimuthal angle coordinate of the solar position vector
    * @param solar_theta  polar angle coordinate of the solar position vector
    * @param img  non-NULL pointer to an uncalibrated image
    * @return 0 on success, non-zero on error
    *
    * The BTDF is the bidirectional transmission distribution function
    * describing the relevant diffuser.  The BTDF factor depends on the direction
    * of the incident solar irradiance.  Applying the BTDF factor modifies
    * the image in place. Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
    */
   int (*cal_apply_btdf)(const Calibration_Type *, double, double,
                         Image_Type *);

   /** Apply PRNU calibration
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param img  non-NULL pointer to an uncalibrated image
    * @return 0 on success, non-zero on error
    *
    * The PRNU is the pixel response non-uniformity.  Applying the PRNU factor modifies
    * the image in place. Pixels containing \a IMAGE_PIXEL_FILL_VALUE are not modified.
    *
    * Not implemented yet
    */
   int (*cal_apply_prnu)(const Calibration_Type *, Image_Type *);

   /** Perform wavelength calibration
    * @param cal  non-NULL pointer to a Calibration_Type object
    * @param img  non-NULL pointer to an uncalibrated image
    * @param img_waves non-NULL pointer to the output calibrated wavelength arrays
    * @return 0 on success, non-zero on error
    *
    * Not implemented yet
    */
   int (*cal_wavecal)(const Calibration_Type *, Image_Type *, Image_Type *);

#ifdef SENSORCAL_PRIVATE_DATA
   SENSORCAL_PRIVATE_DATA
#endif
};

/** Initialize a Calibration_Type object
 * @param cfg  non-NULL pointer to a config_t object associated with an open configuration file
 * @return non-NULL pointer to a Calibration_Type object on success, NULL on error
 */
extern Calibration_Type *sensorcal_init (config_t *cfg);

#endif
