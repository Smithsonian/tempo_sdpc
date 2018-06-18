#ifndef __POLCORR_LPS_INC__
#define __POLCORR_LPS_INC__ 1

/** @file lps.h
 *  @brief Interface for linear polarization sensitivity
 *         lookup table
 */

#include <libconfig.h>

/* FIXME - these indices should be provided by tio_template.h */
enum
{
   TEMPO_BAND_UV = 0,
   TEMPO_BAND_VIS = 1
};

typedef struct Lps_Type Lps_Type;

/** Evaluate linear polarization sensitivity vs wavelength
 * @param[in]  lps          Pointer to \a Lps_Type object allocated by \a lps_open
 * @param[in]  band_index   Integer index indicating the TEMPO band of interest
 *                          0 = UV band, 1 = VIS band
 * @param[in]  xtrack       Integer north-south index, 0 is northern-most
 * @param[in]  lon          Longitude coordinate [deg]
 * @param[in]  lat          Latitude coordinate [deg]
 * @param[in]  n            Number of wavelengths
 * @param[in]  wave         Pointer to array of \a n wavelengths [nm]
 * @param[out] lpsens       Linear polarization sensitivity [dimensionless] vs wavelength
 * @param[out] angmax       Angle [radians] of maximum transmission vs wavelength,
 *                          relative to the instrument reference plane (IRP)
 * @param[out] lmp_irp      Angle [radians] between the local meridian plane (LMP) for the
 *                          specified surface coordinate (lon,lat) and the IRP
 * @return 0 on success, -1 on error
 *
 * The linear polarization sensitivity (LPS) of the instrument is defined in terms
 * of polarization relative to an instrument reference plane, (IRP).  The LPS
 * at each wavelength is characterized by a dimensionless scalar, lpsens, along
 * with the angle of maximum transmission, angmax.
 *
 * For now, the instrument reference plane (IRP) has been arbtrarily defined
 * as the plane containing the instrument slit and the boresight vector.
 * If Ball chooses a different definition, this code may need to be modified
 * to reflect that choice.
 *
 * The look-up tables used to predict the polarization of the backscattered
 * radiance incident on the detector define the polarization relative to
 * the local meridian plane (LMP) which varies with position on the Earth's
 * surface.  The LMP is defined as the plane containing both the surface normal
 * at (lon,lat) and the vector from the surface point to the instrument.
 *
 * Because the polarization reference plane of the instrument differs from
 * the polarization reference plane of the (predicted) incident radiation,
 * the polarization correction requires the angle between these two planes.
 * The required angle, lmp_irp, is the angle between the westward normal of
 * the IRP and the westward normal of the LMP.
 */
extern int
lps_eval (Lps_Type *lps, int band_index, int xtrack,
          double lon, double lat, int n, const double *wave,
          double *linear_polarization_sensitivity,
          double *angle_of_max_transmission,
          double *lmp_irp_angle);

/** Initialize \a Lps_Type object
 * @param[in]  cfg   pointer to open \a config_t object
 * @return pointer to initialized \a Lps_Type on success, NULL on error
 *
 * Data file names and control parameters are retrieved from the
 * configuration file, and an \a Lps_Type is initialized.
 */
extern Lps_Type *lps_open (config_t *cfg);

/** Free resources associated with \a Lps_Type object
 * @param[in]  lps   pointer to \a Lps_Type object allocated by \a lps_open
 */
extern void lps_close (Lps_Type *lps);

#endif
