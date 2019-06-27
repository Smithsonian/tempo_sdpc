#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>

#include <tell.h>

#include "config.h"

#include "interp_poly.h"
#include "lsq.h"

struct Poly_Type
{
   gsl_vector *coef;
   size_t num_coef;
   double xexp;
};

void poly_free (Poly_Type *pt)
{
   if (pt == NULL)
     return;
   gsl_vector_free (pt->coef);
   FREE(pt);
}

static Poly_Type *poly_alloc (size_t num_coef)
{
   Poly_Type *pt = NULL;

   if (NULL == (pt = (Poly_Type *)MALLOC (sizeof *pt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   pt->num_coef = num_coef;
   pt->coef = NULL;

   if (NULL == (pt->coef = gsl_vector_alloc (num_coef)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        poly_free (pt);
        return NULL;
     }

   return pt;
}

int poly_eval (Poly_Type *pt,
               size_t n, const double *x, double *yest)
{
   size_t i, num_coef = pt->num_coef;
   double x0 = pt->xexp;

   for (i = 0; i < n; i++)
     {
        yest[i] = gsl_poly_eval (pt->coef->data, num_coef, x[i]-x0);
     }

   return 0;
}

int poly_get_coeffs (Poly_Type *pt, size_t num_coeffs, double *coeffs)
{
   size_t i;
   if (num_coeffs != pt->num_coef)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size=%ld, %ld is required",
                     __func__, num_coeffs, pt->num_coef);
        return -1;
     }

   for (i = 0; i < pt->num_coef; i++)
     {
        coeffs[i] = gsl_vector_get (pt->coef, i);
     }

   return 0;
}

int poly_set_coeffs (Poly_Type *pt, size_t num_coeffs, const double *coeffs)
{
   size_t i;
   if (num_coeffs != pt->num_coef)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size=%ld, %ld is required",
                     __func__, num_coeffs, pt->num_coef);
        return -1;
     }

   for (i = 0; i < pt->num_coef; i++)
     {
        gsl_vector_set (pt->coef, i, coeffs[i]);
     }

   return 0;
}

static int poly_lsq_set_row (Lsq_Method_Type *m,
                             gsl_matrix *X, size_t i, double xi)
{
   Poly_Type *pt = (Poly_Type *)m->client_data;
   size_t j, num_coef = pt->num_coef;
   double X_ij;

   X_ij = 1.0;
   for (j = 0; j < num_coef; j++)
     {
        gsl_matrix_set (X, i, j, X_ij);
        X_ij *= xi;
     }

   return 0;
}

static void poly_lsq_config (Poly_Type *pt, Lsq_Method_Type *m)
{
   m->client_data = pt;
   m->lsqm_set_row = poly_lsq_set_row;
}

Poly_Type *poly_interpol (size_t num_data, const double *x, const double *y,
                          const Poly_Config_Type *pct)
{
   size_t num_coef = pct->num_coef;
   Lsq_Method_Type m;
   Lsq_Type *lsq = NULL;
   Poly_Type *pt = NULL;
   double chisqr;

   if (NULL == (pt = poly_alloc (num_coef)))
     return NULL;

   pt->xexp = pct->xexp;

   poly_lsq_config (pt, &m);

   if (NULL == (lsq = lsq_alloc (num_data, num_coef, &m)))
     goto error_return;

   if (pt->xexp == 0.0)
     {
        if (0 != lsq_set_data (lsq, num_data, x, y))
          goto error_return;
     }
   else
     {
        double x0 = pt->xexp;
        double *x_x0 = NULL;
        size_t i;
        int status;
        if (NULL == (x_x0 = MALLOC (num_data * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto error_return;
          }
        for (i = 0; i < num_data; i++)
          {
             x_x0[i] = x[i] - x0;
          }
        status = lsq_set_data (lsq, num_data, x_x0, y);
        FREE(x_x0);
        if (status < 0)
          goto error_return;
     }

   if (0 != lsq_solve (lsq, pt->coef, &chisqr))
     goto error_return;

   lsq_free (lsq);

   return pt;

error_return:
   poly_free (pt);
   lsq_free (lsq);
   return NULL;
}
