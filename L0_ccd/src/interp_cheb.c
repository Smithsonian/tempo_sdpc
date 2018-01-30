#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <tell.h>

#include "config.h"

#include "interp_cheb.h"
#include "lsq.h"

struct Cheb_Type
{
   gsl_vector *T;
   gsl_vector *coef;
   size_t num_coef;
   double xmin;
   double xmax;
};

void cheb_free (Cheb_Type *ct)
{
   if (ct == NULL)
     return;
   gsl_vector_free (ct->T);
   gsl_vector_free (ct->coef);
   FREE(ct);
}

static Cheb_Type *cheb_alloc (int num_coef)
{
   Cheb_Type *ct = NULL;

   if (NULL == (ct = (Cheb_Type *)MALLOC (sizeof *ct)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)ct, 0, sizeof *ct);

   ct->num_coef = num_coef;

   if ((NULL == (ct->T = gsl_vector_alloc (num_coef)))
       || (NULL == (ct->coef = gsl_vector_alloc (num_coef))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        cheb_free (ct);
        return NULL;
     }

   return ct;
}

/* Assumes -1 <= x <= 1 */
static int cheb_series_terms (double x, size_t n, double *v)
{
   double x2, vjm1, vjm2;
   size_t j;

   /* Computes one row of the pseudo-Vandermonde matrix of degree
    * (num_coef-1) for Chebyshev polynomials:
    *    V[i,j] = T_j(x_i)
    * where i = row index = 0,...N-1
    *       j = column index = 0,1,...degree
    *
    * and where N = number of data points, (x_i, y_i)
    *
    * e.g.  Y_i = V_ij c_j is the least squares estimate
    * for the data (xi,yi) when the c_j are the linear
    * least squares best fit coefficients.
    *
    * For simplicity, we use variables x=x_i (scalar), v=v_i (vector)
    */

   vjm2 = x*0 + 1.0;
   v[0] = vjm2;
   if (n == 1)
     return 0;

   vjm1 = x;
   v[1] = vjm1;

   x2 = 2*x;

   for (j = 2; j < n; j++)
     {
        double vj = vjm1 * x2 - vjm2;
        v[j] = vj;
        vjm2 = vjm1;
        vjm1 = vj;
     }

   return 0;
}

static int cheb_series_terms_ct (Cheb_Type *ct, double x, gsl_vector *v)
{
   /* t is coordinate in [-1,1] interval */
   double t = (2.0 * x - ct->xmax - ct->xmin) / (ct->xmax - ct->xmin);
   return cheb_series_terms (t, ct->num_coef, v->data);
}

/* Use Clenshaw recursion to evaluate a truncated Chebyshev series,
 *       f(x) = \sum c(k) T(t;k), k = 0,1,2...order
 * at a particular coordinate x, in the interval [a,b], via the Chebyshev
 * coordinate t, on the interval [-1,1], is: t = (2*x-a-b)/(b-a).
 * In this expansion, T(t;k) is a Chebyshev polynomial of the first kind
 * of order k, and the c(k) are constant coefficients.
 */ 
inline int cheb_clenshaw_eval (double a, double b, const double *coef,
                               size_t num_coef, double x, double *value)
{
   size_t i, order = num_coef-1;
   double d1 = 0.0;
   double d2 = 0.0;
   double t, t2;

#if 0
   if (x < a || b < x)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: x=%g out of range (expected values in interval [%g,%g])",
                     __func__, x, a, b);
        return -1;
     }
#endif

   /* Clenshaw recursion: original reference is
    * Clenshaw, C.W., Math. Comp. 9 (1955), 118-120
    */

   /* t is coordinate in [-1,1] interval */
   t = (2.0 * x - a - b) / (b - a);
   t2 = 2.0 * t;

   for (i = order; i >= 1; i--)
     {
        double temp = d1;
        d1 = t2 * d1 - d2 + coef[i];
        d2 = temp;
     }

   /* Note that a common formulation of the Clenshaw recursion
    * has a factor 1/2 in the constant term of the final expression.
    * But in this application, we're generating our Chebyshev
    * series coefficients by solving a least-squares system, so
    * we don't have that factor of 1/2 in the constant term. */

   *value = t * d1 - d2 + coef[0];

   return 0;
}

int cheb_eval (Cheb_Type *ct,
               size_t n, const double *x, double *yest)
{
   double a = ct->xmin;
   double b = ct->xmax;
   double *c = ct->coef->data;
   size_t m = ct->coef->size;
   size_t i;

   for (i = 0; i < n; i++)
     {
        if (0 != cheb_clenshaw_eval (a, b, c, m, x[i], &yest[i]))
          return -1;
     }

   return 0;
}

int cheb_get_coeffs (Cheb_Type *ct, size_t num_coeffs, double *coeffs)
{
   size_t i;
   if (num_coeffs != ct->num_coef)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size=%ld, %ld is required",
                     __func__, num_coeffs, ct->num_coef);
        return -1;
     }

   for (i = 0; i < ct->num_coef; i++)
     {
        coeffs[i] = gsl_vector_get (ct->coef, i);
     }

   return 0;
}

int cheb_set_coeffs (Cheb_Type *ct, size_t num_coeffs, const double *coeffs)
{
   size_t i;
   if (num_coeffs != ct->num_coef)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size=%ld, %ld is required",
                     __func__, num_coeffs, ct->num_coef);
        return -1;
     }

   for (i = 0; i < ct->num_coef; i++)
     {
        gsl_vector_set (ct->coef, i, coeffs[i]);
     }

   return 0;
}

static int cheb_lsq_set_row (Lsq_Method_Type *m,
                             gsl_matrix *X, size_t i, double xi)
{
   Cheb_Type *ct = (Cheb_Type *)m->client_data;
   size_t j;

   if (0 != cheb_series_terms_ct (ct, xi, ct->T))
     return -1;

   for (j = 0; j < ct->num_coef; j++)
     {
        double Tij = gsl_vector_get (ct->T, j);
        gsl_matrix_set (X, i, j, Tij);
     }

   return 0;
}

static void cheb_lsq_config (Cheb_Type *ct, Lsq_Method_Type *m)
{
   m->client_data = ct;
   m->lsqm_set_row = cheb_lsq_set_row;
}

Cheb_Type *cheb_interpol (size_t num_data, const double *x, const double *y,
                          const Cheb_Config_Type *cct)
{
   size_t num_coef = cct->num_coef;
   Lsq_Method_Type m;
   Lsq_Type *lsq = NULL;
   Cheb_Type *ct = NULL;
   double chisqr;

   if (NULL == (ct = cheb_alloc (num_coef)))
     return NULL;

   if (cct->xmax > cct->xmin)
     {
        ct->xmin = cct->xmin;
        ct->xmax = cct->xmax;
     }
   else
     {
        ct->xmin = x[0];
        ct->xmax = x[num_data - 1];
     }

   cheb_lsq_config (ct, &m);

   if (NULL == (lsq = lsq_alloc (num_data, num_coef, &m)))
     goto error_return;

   if (0 != lsq_set_data (lsq, num_data, x, y))
     goto error_return;

   if (0 != lsq_solve (lsq, ct->coef, &chisqr))
     goto error_return;

   lsq_free (lsq);

   return ct;

error_return:
   cheb_free (ct);
   lsq_free (lsq);
   return NULL;
}
