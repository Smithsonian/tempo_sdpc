#ifndef __SLIT_FUNCTION_INCLUDE_H__
#define __SLIT_FUNCTION_INCLUDE_H__ 1

/** @file slit_function.h
 *  @brief Slit function convolution
 *
 * This module provides an interface for convolving a high resolution
 * spectrum with an instrumental slit function, a.k.a. instrumental spectral
 * response function (ISRF).
 */

#define SFT_NUM_PARAMS (3)

typedef struct Slit_Function_Type Slit_Function_Type;

/** Function interface for evaluating a slit function
 */
typedef int SFT_Eval_Type (const double *x, size_t nx, double *params,
                           double *value, double *param_step,
                           double *param_derivs[SFT_NUM_PARAMS]);

/** Function interface for retrieving wavelength[k] slit function parameters
 */
typedef int SFT_Param_Type (int k, int num_pars, double *pars, void *cl);

/** Free storage associated with \a Slit_Function_Type structure
 */
extern void sft_free (Slit_Function_Type *);

/** Allocate a new \a Slit_Function_Type structure
 * @param[in]  num_sf      Number of wavelength points in the computed slit function shape.
 */
extern Slit_Function_Type *sft_new (int num_sf);

/** Configure \a Slit_Function_Type structure
 * @param[in]  sft         Pointer to a \a Slit_Function_Type structure allocated by \a sft_new
 * @param[in]  sf_eval     Function to evaluate the slit function
 * @param[in]  dx          Wavelength grid spacing [nm].  For slit function convolution, the
 *                         wavelength grid is assumed to have a fixed spacing.
 * @param[in]  param_step  Parameter delta used in computing numerical derivatives of
 *                         the slit-function with respect to each parameter.
 */
extern int sft_config (Slit_Function_Type *sft, SFT_Eval_Type *sf_eval,
                       double dx, double *param_step);

/** Apply slit function to a padded spectrum array
 * @param[in] sft          Pointer to a \a Slit_Function_Type structure allocated by \a sft_new,
 *                         and configured by \a sft_config
 * @param[in] sf_params    Function to retrieve slit function parameters for each wavelength index
 * @param[in] cl           Client data to pass to function \a sf_params
 * @param[in] num_waves    Number of wavelengths in the target spectrum
 * @param[in] spec_padded  Target spectrum with \a num_waves wavelengths, padded with \a num_sf/2
 *                         zeros at each end, so that the total array size is \a (num_waves+num_sf)
 * @param[out] spec_convolved  Output convolved spectrum of length \a num_waves
 * @param[in]  compute_derivs  Integer parameter.  When non-zero, the output will include slit
 *                             function parameter derivatives convolved with the target spectrum.
 * @param[out] spec_derivs_convolved   Array of pointers to arrays of length \a num_waves;
 *                                     storage for output slit function parameter
 *                                     derivatives convolved with the target spectrum.
 */
extern int sft_apply (Slit_Function_Type *sft, SFT_Param_Type *sf_params, void *cl,
                      int num_waves, const double *spec_padded, double *spec_convolved,
                      int compute_derivs, double *spec_derivs_convolved[SFT_NUM_PARAMS]);

#endif
