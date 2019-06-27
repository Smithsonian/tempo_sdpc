#ifndef __WAVECAL_ADJ_INCLUDE__
#define __WAVECAL_ADJ_INCLUDE__ 1

#include <libconfig.h>

/** @file wavecal_adj.h
 *  @brief Wavelength shift adjustment
 *
 * This module provides an interface for applying a wavelength shift adjustment
 * to a pre-determined set of wavelength grids defined by a set of Chebyshev
 * series coefficients.
 */

typedef struct Wadj_Type Wadj_Type;

typedef struct
{
   int pix_min;           /**< pixel index of spectral window lower boundary (inclusive) */
   int pix_max;           /**< pixel index of spectral window upper boundary (inclusive) */
   int num_series_coeff;  /**< number of Chebyshev series coefficients used to define wavelength grid in [pix_min,pix_max] */
}
Wadj_Cheb_Series_Type;

/** Read the wavelength adjustment lookup table for the specified netcdf file group (spectral band).
 * @param  cfg    Pointer to an open \a config_t structure initialized
 *                with an appropriate configuration file.
 * @param  group  Name of the netcdf file group containing the lookup table for the
 *                spectral band of interest.
 * @param  meta   Pointer to an open \a TIO_Meta_type structure
 * @return A valid \a Wadj_Type pointer on success, \a NULL on error
 */
extern Wadj_Type *wadj_open (config_t *cfg, const char *group, TIO_Meta_Type *meta);

/** Release resources associated with a wavelength adjustment lookup table
 * @param   wadj  Pointer to a \a Wadj_Type pointer initialized by \a wadj_open
 */
extern void wadj_close (Wadj_Type *wadj);

/** Query the narrow-band definition used in the wavelength adjustment lookup table
 * @param   wadj  Pointer to a \a Wadj_Type pointer initialized by \a wadj_open.
 * @param   cheb  Pointer to a target \a Wadj_Cheb_Series_Type structure to hold
 *                the result values.
 * @return 0 on success, -1 on error
 */
extern int wadj_narrow_band_get_attr (const Wadj_Type *wadj, Wadj_Cheb_Series_Type *cheb);

/** Query the full-band definition used in the wavelength adjustment lookup table
 * @param   wadj  Pointer to a \a Wadj_Type pointer initialized by \a wadj_open.
 * @param   cheb  Pointer to a target \a Wadj_Cheb_Series_Type structure to hold
 *                the result values.
 * @return 0 on success, -1 on error
 */
extern int wadj_full_band_get_attr (const Wadj_Type *wadj, Wadj_Cheb_Series_Type *cheb);

/** Obtain a read-only pointer to the narrow-band Chebyshev series coefficients
 * @param   wadj  Pointer to a \a Wadj_Type pointer initialized by \a wadj_open.
 * @param   xtrack  Cross-track pixel index
 * @return a valid, read-only pointer on success, NULL on error
 */
extern const double *wadj_narrow_band_coeff (const Wadj_Type *wadj, size_t xtrack);

/** Query the final number of Chebyshev series coefficients needed after applying the wavelength shift adjustment.
 * @param   wadj       Pointer to a \a Wadj_Type pointer initialized by \a wadj_open.
 * @param   num_coeff  Final number of Chebyshev series coefficients
 * @return 0 on success, -1 on error
 */
extern int wadj_num_final_coeff (const Wadj_Type *wadj, int *num_coeff);

/** Compute the final Chebyshev series coefficients, including the wavelength shift adjustment.
 * @param   wadj      Pointer to a \a Wadj_Type pointer initialized by \a wadj_open.
 * @param   xtrack    Cross-track pixel index
 * @param   shift     Measured wavelength shift, (fit_y0 - table_y0)
 * @param   coeff     Pointer to an array of size consistent with the return value from \a wadj_num_final_coeff
 * @param 0 on success, -1 on error
 */
extern int wadj_final_coeff (const Wadj_Type *wadj, size_t xtrack, double shift, double *coeff);

#endif
