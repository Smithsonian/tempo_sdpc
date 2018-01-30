#ifndef __INTERP_CHEBYSHEV_INCLUDE__
#define __INTERP_CHEBYSHEV_INCLUDE__ 1

/** @file interp_cheb.h
 *  @brief Chebyshev interpolation
 *
 * This module performs 1-D interpolation using a Chebyshev
 * series expansion of the form:
 *
 *    f(t) = \sum c(k) T(t;k)  k = 0,1,..N-1
 *
 * where t is in the range [-1,1],
 *     T(t;k) is a Chebyshev polynomial of the first kind of order k
 * and where the constant coefficients, c(k), are derived by
 * a least-squares fit to a set of data points, F(i), t(i), i=1..M
 *
 * For functions defined on an interval x in [a,b], define
 *     t = (2*x - a - b)/(b - a)
 */

typedef struct Cheb_Type Cheb_Type;

typedef struct
{
   double xmin;     /**< lower bound of X-coordinate interval */
   double xmax;     /**< upper bound of X-coordinate interval */
   size_t num_coef; /**< number of Chebyshev coefficients */
}
Cheb_Config_Type;

/** Free memory associated with a \a Cheb_Type object
 * @param  ct  A pointer to a \a Cheb_Type object
 */
extern void cheb_free (Cheb_Type *ct);

/** Initialize a \a Cheb_Type object
 * @param  n   The number of data points
 * @param  x   Read-only pointer to an array of X coordinates
 * @param  y   Read-only pointer to an array of Y coordinates
 * @param  config  Read-only pointer to an initialized struct of type \a Cheb_Config_Type
 *
 * @return a non-NULL pointer to a \a Cheb_Type object on success,
 * NULL on error.
 */
extern Cheb_Type *
cheb_interpol (size_t n, const double *x, const double *y,
               const Cheb_Config_Type *config);

/** Interpolate using a least-squares Chebyshev series.
 * @param  ct  Pointer to an initialized \a Cheb_Type object.
 * @param  n   The number of data points
 * @param  x   Read-only pointer to an array of X coordinates
 * @param  y   Pointer to a pre-allocated array to hold the interpolated
 *             Y values.
 *
 * @return 0 on success, -1 on error.
 */
extern int
cheb_eval (Cheb_Type *ct, size_t n, const double *x, double *yest);

/** Retrieve Chebyshev coefficients from a \a Cheb_Type object.
 * @param ct           Pointer to an initialized \a Cheb_Type object.
 * @param num_coeffs   The number of coefficients to be retrieved
 *                     (which must exactly match the value used to
 *                      initialize the \a Cheb_Type object).
 * @param coeffs       Pointer to a pre-allocated array to hold the
 *                     Chebyshev coefficients.
 * @return 0 on success, -1 on error.
 *
 * The coefficients are stored so that \a coeffs[k] is the coefficient
 * of, \a T(t;k), the Chebyshev polynomial of order \a k.
 */
extern int cheb_get_coeffs (Cheb_Type *ct, size_t num_coeffs,
                            double *coeffs);

/** Set Chebyshev coefficients in a \a Cheb_Type object.
 * @param ct           Pointer to an initialized \a Cheb_Type object.
 * @param num_coeffs   The number of coefficient
 *                     (which must exactly match the value used to
 *                      initialize the \a Cheb_Type object).
 * @param coeffs       Read-only pointer to an array of
 *                     Chebyshev coefficients.
 * @return 0 on success, -1 on error.
 *
 * The coefficients are stored so that \a coeffs[k] is the coefficient
 * of, \a T(t;k), the Chebyshev polynomial of order \a k.
 */
extern int cheb_set_coeffs (Cheb_Type *ct, size_t num_coeffs,
                            const double *coeffs);

/** Evaluate a Chebyshev series at a point x using Clenshaw's recursion
 * @param  a   Interpolation domain X coordinate lower bound
 * @param  b   Interpolation domain X coordinate upper bound
 * @param  c   Read-only pointer to an array of Chebyshev series coefficients,
 *             with \a c[k] the coefficient of the Chebyshev polynomial of order k.
 * @param  m   The number of series coefficients.
 * @param  x   The X-coordinate where the series should be evaluated.
 * @param  value  The computed value of the series at the specified X-coordinate
 * @return 0 on success, -1 on error.
 *
 * As far as I know, Clenshaw's recursion is the most efficient
 * and numerically stable way to evaluate a Chebyshev series.
 */
extern int cheb_clenshaw_eval (double a, double b, const double *c,
                               size_t m, double x, double *value);

#endif
