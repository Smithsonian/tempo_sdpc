#ifndef __INTERP_CSPLINE_INCLUDE__
#define __INTERP_CSPLINE_INCLUDE__ 1

/** @file interp_cspline.h
 *  @brief Cubic spline interpolation
 *
 * This module performs 1-D cubic spline interpolation.
 */

typedef struct Cspline_Type Cspline_Type;
typedef struct Cspline_Config_Type Cspline_Config_Type;

struct Cspline_Config_Type
{
   int unused;
};

/** Free memory associated with a \a Cspline_Type object
 * @param   A pointer to a \a Cspline_Type object
 */
extern void cspline_free (Cspline_Type *ct);

/** Initialize a \a Cspline_Type object
 * @param  n   The number of data points
 * @param  x   Read-only pointer to an array of X coordinates
 * @param  y   Read-only pointer to an array of Y coordinates
 * @param  config  Read-only pointer to an initialized struct of type \a Cspline_Config_Type
 *
 * @return a non-NULL pointer to a \a Cspline_Type object on success,
 * NULL on error.
 */
extern Cspline_Type *
cspline_interpol (size_t n, const double *x, const double *y,
                  void *config);

/** Interpolate using cubic splines
 * @param  ct  Pointer to an initialized \a Cspline_Type object.
 * @param  n   The number of data points
 * @param  x   Read-only pointer to an array of X coordinates
 * @param  y   Pointer to a pre-allocated array to hold the interpolated
 *             Y values.
 *
 * @return 0 on success, -1 on error.
 */
extern int
cspline_eval (Cspline_Type *ct, size_t n, const double *x, double *yest);
#endif
