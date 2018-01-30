#ifndef __WCI_LSQ_INCLUDE__
#define __WCI_LSQ_INCLUDE__ 1

/** @file lsq.h
 *  @brief Linear least-squares optimization
 *
 * This module provides an interface for solving reasonably
 * generic linear least-squares problems.
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct Lsq_Type Lsq_Type;
typedef struct Lsq_Method_Type Lsq_Method_Type;

struct Lsq_Method_Type
{
   /** Compute one row of elements for the pseudo-Vandermonde
    * matrix defining the least-squares linear system.
    * @param  m  A non-null pointer to this \a Lsq_Method_Type object
    * @param  X  A non-null pointer to a GSL matrix object to receive
    *            the row elements
    * @param  i  The index of the row to be populated.
    * @param xi  The value of the ith X-coordinate (corresponding to this row).
    *
    * @return 0 on success, -1 on error.
    */
   int (*lsqm_set_row)(Lsq_Method_Type *m, gsl_matrix *X, size_t i, double xi);

   /** Pointer to method-specific private data */
   void *client_data;
};

/** Free memory associated with an \a Lsq_Type object
 * @param  lsq  A pointer to an \a Lsq_Type object
 */
extern void lsq_free (Lsq_Type *lsq);

/** Allocate an \a Lsq_Type object for linear least-squares fitting
 * @param  num_data   The number of data points
 * @param  num_coef   The number of fit coefficients
 * @param  method     Pointer to a fully initialized \a Lsq_Method_Type
 *                    structure.
 * @return a non-NULL pointer to  an \a Lsq_Type object on success,
 * NULL on error.
 */
extern Lsq_Type *lsq_alloc (size_t num_data, size_t num_coef,
                            Lsq_Method_Type *method);

/** Provide (x,y) data points for a linear least-squares fit
 * @param  lsq    An \a Lsq_Type object created by \a lsq_alloc
 * @param  n      The number of data points
 * @param  x      A pointer to an array of \a n \a x[i] data values
 * @param  y      A pointer to an array of \a n \a y[i] data values
 *
 * @return 0 on success, -1 on error.
 */
extern int lsq_set_data (Lsq_Type *lsq, size_t n,
                         const double *x, const double *y);

/** Solve a linear least-squares system
 * @param  lsq    An \a Lsq_Type object created by \a lsq_alloc
 *                and populated with data using \a lsq_set_data.
 * @param  coef   A pointer to an output GSL vector containing the
 *                least-squares solution coefficients.  Storage
 *                for the array must have been previously allocated.
 * @param chisqr  The address in which the least-squares fit statistic
 *                is to be stored.
 *
 * @return 0 on success, -1 on error.
 */
extern int lsq_solve (Lsq_Type *lsq, gsl_vector *coef, double *chisqr);

#endif
