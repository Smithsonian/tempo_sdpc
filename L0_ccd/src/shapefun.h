#ifndef __SHAPEFUN_INCLUDE__
#define __SHAPEFUN_INCLUDE__ 1

/** @file shapefun.h
 *  @brief Wavelength-dependent shape functions for wavelength calibration
 *
 * This module facilitates wavelength calibration by providing
 * wavelength-dependent shape functions for representing various
 * model spectrum components.  For example, the effect of scattering
 * in the atmosphere may be accounted for by introducing a wavelength
 * dependent coefficient on the absorbing column.
 */

#include <libconfig.h>

enum
{
   SHAPEFUN_TYPE_CSPLINE,
   SHAPEFUN_TYPE_CHEB,
   SHAPEFUN_TYPE_POLY,
   SHAPEFUN_TYPE_POSITIVE_COEFF
};

typedef struct
{
   /** The initial state of a shape function may be specified
    * by a set of (x,y) pairs.
    */
   double *x;     /**< array of X coordinate values */
   double *y;     /**< array of Y coordinate values */
   size_t n;      /**< number of (X,Y) coordinate pairs */
   int malloced;  /**< optional flag to facilitate memory management */
}
Shapefun_Init_Type;

typedef struct Shapefun_Type Shapefun_Type;

struct Shapefun_Type
{
   /** A destructor.
    *  @param  sf   A non-null pointer to this \a Shapefun_Type object
    *
    * Call this function to release resources associated with this
    * \a Shapefun_Type object.
    */
   void (*st_delete)(Shapefun_Type *);

   /** Retrieve the number of parameters associated with this \a Shapefun_Type object
    * @param sf  A pointer to this \a Shapefun_Type object
    * @return The number of parameters
    */
   int (*st_num_params)(const Shapefun_Type *);

   /** Query the method id for this \a Shapefun_Type object
    * @param sf  A pointer to this \a Shapefun_Type object
    * @return The method id
    */
   int (*st_method)(const Shapefun_Type *);

   /** Initialize the parameters for this \a Shapefun_Type object
    * @param  sf  A read-only pointer to this \a Shapefun_Type object
    * @param  si  A read-only pointer to an initialized \a Shapefun_Init_Type object
    * @param  n   The number of parameters for this \a Shapefun_Type object
    *             (see the \a st_num_params method)
    * @param  p   An output array of parameter values.
    *
    * For example, if the \a Shapefun_Init_Type provides a set of (x,y)
    * pairs defining a shape, then those points will be used to initialize
    * an interpolation method which will smoothly interpolate to any
    * valid X coordinate in the domain.
    *
    * Note that, depending on the method, the parameters may be either
    * least-squares coefficients (e.g. for Chebyshev interpolation) or
    * they may be interpolated function values for a fixed set of nodes
    * (e.g. for the spline interpolation method).
    *
    * @return 0 on success, -1 on error
    */
   int (*st_init_params)(const Shapefun_Type *, const Shapefun_Init_Type *,
                         size_t, double *);

   /** Evaluate this \a Shapefun_Type object at a set of X coodinates.
    * @param sf  A read-only pointer to this \a Shapefun_Type object
    * @param p   A read-only array of parameter values
    * @param n   The number of X coordinates for which interpolated values
    *            are wanted.
    * @param x   A read-only array of X coordinate values.
    * @param yest  An array of interpolated values.
    *
    * @return 0 on success, -1 on error
    */
   int (*st_eval)(const Shapefun_Type *, const double *,
                  size_t, const double *, double *);

   int st_apply_external_scaling;  /**< if non-zero, result requires externally managed scale factor */

   double *node_x;        /**< [method-specific option] fixed, ordered set of x coordinates */
   size_t num_nodes;      /**< [method-specific option] number of fixed nodes */
   int malloced_node_x;   /**< [method-specific option] flag used to facilitate memory management */

   size_t num_coef;       /**< [method-specific option] number of least-squares coefficients */
   double xmin;           /**< [method-specific option] interpolation domain X coordinate lower bound */
   double xmax;           /**< [method-specific option] interpolation domain X coordinate upper bound */
   double xexp;           /**< [method-specific option] X-coordinate for polynomial expansion */

#ifdef SHAPEFUN_PRIVATE_DATA
   SHAPEFUN_PRIVATE_DATA
#endif
};

/** Create an \a Shapefun_Type object based on a particular interpolation method
 * @param method_name   The name of the interpolation method
 * @return A valid \a Shapefun_Type pointer on success, NULL on error.
 */
extern Shapefun_Type *shapefun_create (const char *method_name);

#endif
