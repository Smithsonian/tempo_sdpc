#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>

#include "config.h"
#include "lsq.h"

/* Constructing a linear least squares fit to N data points (x_i, y_i),
 * with weights, w_i, means that we minimize the chi-square statistic:
 *       \chi^2 = \sum_i w_i * (y_i - X_ij c_j)^2
 * where
 *      y_i = data value
 *      w_i = data weight
 *      c_j = fit coefficient
 *     X_ij = matrix of "predictor variables"
 *
 * If we're fitting a polynomial, then the fitted y value is
 *    Y_i = X_ij c_j,
 * so that the corresponding "predictor variable" matrix is
 *     X_ij = (x_i)^j
 * where j is the jth order polynomial term
 *
 * Similarly, if we're fitting with B-splines, then
 *    X_ij = B_j(x_i)
 * where B_j(x) is the jth order B-spline function.
 * Recall that B-splines form a basis for the vector space of
 * all spline functions. The B-spline basis functions are
 * constructed using a user-specified set of knots within
 * each sub interval.
 *
 * If we're fitting with Chebyshev polynomials, then
 *    X_ij = T_j(x_i)
 * where T_j(x) is the jth order Chebyshev polynomial.
 * Recall that the Chebyshev polynomials are defined on
 * the range [-1,1] using symmetrically placed abcissae,
 * so an appropriate linear transformation is implied.
 *
 * To minimize chi-square, we use GSL routines that
 * take as input the vectors (x_i, y_i, w_i) and the matrix X_ij,
 * and compute the least-squares best-fit coefficients, c_j.
 *
 * NOTE: for the current application, we don't expect to use
 * the weights for linear least-squares, so I'm omitting support
 * for that.
 */

struct Lsq_Type
{
   gsl_vector *x;
   gsl_vector *y;
   gsl_matrix *X;
   gsl_matrix *cov;
   gsl_multifit_linear_workspace *mw;
   size_t num_data;
   size_t num_coef;
   Lsq_Method_Type *method;
};

void lsq_free (Lsq_Type *lsq)
{
   if (lsq == NULL)
     return;
   gsl_vector_free (lsq->x);
   gsl_vector_free (lsq->y);
   gsl_matrix_free (lsq->X);
   gsl_matrix_free (lsq->cov);
   gsl_multifit_linear_free (lsq->mw);
   FREE(lsq);
}

Lsq_Type *lsq_alloc (size_t num_data, size_t num_coef,
                     Lsq_Method_Type *method)
{
   Lsq_Type *lsq = NULL;

   if (NULL == (lsq = (Lsq_Type *)MALLOC (sizeof *lsq)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)lsq, 0, sizeof *lsq);

   lsq->num_data = num_data;
   lsq->num_coef = num_coef;
   lsq->method = method;

   if ((NULL == (lsq->x = gsl_vector_alloc(num_data)))
       || (NULL == (lsq->y = gsl_vector_alloc(num_data)))
       || (NULL == (lsq->X = gsl_matrix_alloc(num_data, num_coef)))
       || (NULL == (lsq->cov = gsl_matrix_alloc(num_coef, num_coef)))
       || (NULL == (lsq->mw = gsl_multifit_linear_alloc (num_data, num_coef)))
       )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        lsq_free (lsq);
        return NULL;
     }

   return lsq;
}

static int lsq_set_matrix (Lsq_Type *lsq)
{
   Lsq_Method_Type *m = lsq->method;
   size_t i;

   for (i = 0; i < lsq->num_data; i++)
     {
        double xi = gsl_vector_get (lsq->x, i);
        if (0 != m->lsqm_set_row (m, lsq->X, i, xi))
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: matrix element evaluation failed",
                          __func__);
             return -1;
          }
     }

   return 0;
}

int lsq_set_data (Lsq_Type *lsq, size_t n,
                  const double *x, const double *y)
{
   size_t i;

   if (n != lsq->num_data)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: number of data values = %ld (expected %ld)",
                     __func__, n, lsq->num_data);
        return -1;
     }

   for (i = 0; i < n; i++)
     gsl_vector_set (lsq->x, i, x[i]);

   for (i = 0; i < n; i++)
     gsl_vector_set (lsq->y, i, y[i]);

   return 0;
}

int lsq_solve (Lsq_Type *lsq, gsl_vector *coef, double *chisqr)
{
   int status;

   if (0 != lsq_set_matrix (lsq))
     return -1;

   if (0 != (status = gsl_multifit_linear (lsq->X, lsq->y, coef,
                                           lsq->cov, chisqr, lsq->mw)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: %s",
                     __func__, gsl_strerror(status));
        return -1;
     }

   return 0;
}
