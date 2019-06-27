#ifndef __INTERP_INCLUDE__
#define __INTERP_INCLUDE__ 1

/** @file interp.h
 *  @brief Generic interpolation
 *
 * This module provides an interface for performing 1-D interpolation.
 */

enum
{
   INTERP_TYPE_CSPLINE = 0,
   INTERP_TYPE_CHEB    = 1,
   INTERP_TYPE_POLY    = 2
};

typedef struct Interp_Type Interp_Type;

struct Interp_Type
{
   /** Free memory associated with this \a Interp_Type object
    * @param   A pointer to this \a Interp_Type object
    */
   void (*it_delete)(Interp_Type *);

   /** Retrieve the method index for this interpolation method
    * @param   A pointer to this \a Interp_Type object
    * @return  A valid method index.
    */
   int (*it_method_id)(const Interp_Type *);

   /** Initialize the interpolation method using a set of data points
    * @param  it  A pointer to this \a Interp_Type object
    * @param  n   The number of (X,Y) data points
    * @param  x   A read-only pointer to an array of X coordinates.
    * @param  y   A read-only pointer to an array of Y coordinates.
    * @return  0 on success, -1 on error.
    */
   int (*it_interp_init)(Interp_Type *, size_t, const double *, const double *);

   /** Interpolate the function value at a set of X coordinates
    * @param  it  A pointer to this \a Interp_Type object
    * @param  n   The number of X coordinates
    * @param  x   A read-only pointer to an array of X coordinates.
    * @param  y   A pointer to a pre-allocated array to hold the
    *             interpolated Y values.
    * @return  0 on success, -1 on error.
    */
   int (*it_interp_eval)(Interp_Type *, size_t, const double *, double *);

   size_t num_coef;   /**< [method-dependent option] the number of interpolation coefficients */
   double xmin;       /**< [method-dependent option] lower bound on X coordinate domain */
   double xmax;       /**< [method-dependent option] upper bound on X coordinate domain */
   double xexp;       /**< [method-dependent option] X coordinate for polynomial expansion */

#ifdef INTERP_TYPE_PRIVATE_DATA
   INTERP_TYPE_PRIVATE_DATA
#endif
};

/** Create an \a Interp_Type object for a particular interpolation method
 * @param method_name   The name of the interpolation method
 * @return A valid \a Interp_Type pointer on success, NULL on error.
 */
extern Interp_Type *interp_create (const char *method_name);

#endif
