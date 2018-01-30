#ifndef __INTERP_POLY_INCLUDE__
#define __INTERP_POLY_INCLUDE__ 1

/** @file interp_poly.h
 *  @brief Polynomial interpolation
 *
 * This module performs 1-D interpolation using a
 * polynomial expansion of the form:
 *
 *    f(x) = \sum c(k) (x-x0)^k  k = 0,1,..N-1
 *
 * where the constant coefficients, c(k), are derived by
 * a least-squares fit to a set of data points, X(i), Y(i), i=1..M
 */

typedef struct Poly_Type Poly_Type;

typedef struct
{
   double xexp;        /**< X-coordinate for series expansion */
   size_t num_coef;    /**< number of series coefficients */
}
Poly_Config_Type;

/** Free memory associated with a \a Poly_Type object
 * @param pt  A pointer to a \a Poly_Type object
 */
extern void poly_free (Poly_Type *pt);

/** Initialize a \a Poly_Type object
 * @param  n   The number of data points
 * @param  x   Read-only pointer to an array of X coordinates
 * @param  y   Read-only pointer to an array of Y coordinates
 * @param  config  Read-only pointer to an initialized struct of type \a Poly_Config_Type
 *
 * @return a non-NULL pointer to a \a Poly_Type object on success,
 * NULL on error.
 */
extern Poly_Type *
poly_interpol (size_t n, const double *x, const double *y,
               const Poly_Config_Type *config);

/** Interpolate using a least-squares polynomial.
 * @param  pt  Pointer to an initialized \a Poly_Type object.
 * @param  n   The number of data points
 * @param  x   Read-only pointer to an array of X coordinates
 * @param  yest   Pointer to a pre-allocated array to hold the interpolated
 *             Y values.
 *
 * @return 0 on success, -1 on error.
 */
extern int
poly_eval (Poly_Type *pt, size_t n, const double *x, double *yest);

/** Retrieve polynomial coefficients from a \a Poly_Type object.
 * @param pt           Pointer to an initialized \a Poly_Type object.
 * @param num_coeffs   The number of coefficients to be retrieved
 *                     (which must exactly match the value used to
 *                      initialize the \a Poly_Type object).
 * @param coeffs       Pointer to a pre-allocated array to hold the
 *                     polynomial coefficients.
 * @return 0 on success, -1 on error.
 *
 * The coefficients are stored so that \a coeffs[k] is the coefficient
 * of, \a (x-x0)^k, the polynomial term of order \a k.
 */
extern int poly_get_coeffs (Poly_Type *pt, size_t num_coeffs,
                            double *coeffs);

/** Set polynomial coefficients in a \a Poly_Type object.
 * @param pt           Pointer to an initialized \a Poly_Type object.
 * @param num_coeffs   The number of coefficients to be retrieved
 *                     (which must exactly match the value used to
 *                      initialize the \a Poly_Type object).
 * @param coeffs       Read-only pointer to an array holding the
 *                     polynomial coefficients.
 * @return 0 on success, -1 on error.
 *
 * The coefficients are stored so that \a coeffs[k] is the coefficient
 * of, \a (x-x0)^k, the polynomial term of order \a k.
 */
extern int poly_set_coeffs (Poly_Type *pt, size_t num_coeffs,
                            const double *coeffs);

#endif
