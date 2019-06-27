#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>

#include <tell.h>

#include "config.h"
#include "interp_bspline.h"
#include "lsq.h"
#include "util.h"

#define MAX(a,b) (((a)>(b)) ? (a) : (b))

struct Bspline_Type
{
   gsl_bspline_workspace *bw;
   gsl_vector *B;
   gsl_vector *coef;
   size_t num_coef;
   size_t nbreak;
   double xmin;
   double xmax;
};

/* k=(1+m) where m is the polynomial order of the B-spline basis functions.
 * For cubic splines, the order is m=3, so k=4 */
#define BSPLINE_K 4

void bspline_free (Bspline_Type *bt)
{
   if (bt == NULL)
     return;
   gsl_bspline_free (bt->bw);
   gsl_vector_free (bt->B);
   gsl_vector_free (bt->coef);
   FREE(bt);
}

static Bspline_Type *bspline_alloc (int nbreak)
{
   Bspline_Type *bt = NULL;

   if (NULL == (bt = (Bspline_Type *)MALLOC (sizeof *bt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)bt, 0, sizeof *bt);

   if (NULL == (bt->bw = gsl_bspline_alloc (BSPLINE_K, nbreak)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        bspline_free (bt);
        return NULL;
     }

   bt->nbreak = nbreak;
   bt->num_coef = gsl_bspline_ncoeffs (bt->bw);

   fprintf (stderr, "bspline:  ncoef=%ld nbreak=%ld\n",
            bt->num_coef, bt->nbreak);

   if ((NULL == (bt->B = gsl_vector_alloc (bt->num_coef)))
       || (NULL == (bt->coef = gsl_vector_alloc (bt->num_coef))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        bspline_free (bt);
        return NULL;
     }

   return bt;
}

static int bspline_define_knots (Bspline_Type *bt, const double *breakpts,
                                 size_t num_breakpts)
{
   gsl_vector *bpts = NULL;
   int status;
   size_t i;

   if (breakpts == NULL)
     {
        if (0 != (status = gsl_bspline_knots_uniform (bt->xmin, bt->xmax, bt->bw)))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: %s",
                          __func__, gsl_strerror (status));
             return -1;
          }
        return 0;
     }

   if (NULL == (bpts = gsl_vector_alloc (bt->nbreak)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < num_breakpts; i++)
     {
        gsl_vector_set (bpts, i, breakpts[i]);
     }

   status = gsl_bspline_knots (bpts, bt->bw);
   gsl_vector_free (bpts);

   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: %s",
                     __func__, gsl_strerror (status));
        return -1;
     }

   return 0;
}

int bspline_eval (Bspline_Type *bt,
                  size_t n, const double *x, double *yest)
{
   size_t i;

   for (i = 0; i < n; i++)
     {
        if (0 != gsl_bspline_eval (x[i], bt->B, bt->bw))
          return -1;
        if (0 != gsl_blas_ddot (bt->B, bt->coef, &yest[i]))
          return -1;
     }

   return 0;
}

static int bspline_lsq_set_row (Lsq_Method_Type *m,
                                gsl_matrix *X, size_t i, double xi)
{
   Bspline_Type *bt = (Bspline_Type *)m->client_data;
   size_t j;

   gsl_bspline_eval (xi, bt->B, bt->bw);

   for (j = 0; j < bt->num_coef; j++)
     {
        double Bj = gsl_vector_get (bt->B, j);
        gsl_matrix_set (X, i, j, Bj);
     }

   return 0;
}

static void bspline_lsq_config (Bspline_Type *bt, Lsq_Method_Type *m)
{
   m->client_data = bt;
   m->lsqm_set_row = bspline_lsq_set_row;
}

Bspline_Type *bspline_interpol (size_t num_data, const double *x, const double *y,
                                const Bspline_Config_Type *bct)
{
   Lsq_Method_Type m;
   Lsq_Type *lsq = NULL;
   Bspline_Type *bt = NULL;
   const double *breakpts = NULL;
   double chisqr;
   size_t num_breakpts = 0;

   if (bct)
     {
        breakpts = bct->breakpts;
        num_breakpts = bct->num_breakpts;
     }

   if (num_breakpts == 0)
     {
        num_breakpts = num_data + BSPLINE_K;
     }

   if (NULL == (bt = bspline_alloc (num_breakpts)))
        return NULL;

   if ((bct == NULL) || (bct->xmin == bct->xmax))
     {
        bt->xmin = x[0];
        bt->xmax = x[num_data-1];
     }
   else
     {
        bt->xmin = bct->xmin;
        bt->xmax = bct->xmax;
     }

   if (0 != bspline_define_knots (bt, breakpts, num_breakpts))
     goto error_return;

   bspline_lsq_config (bt, &m);

   if (NULL == (lsq = lsq_alloc (num_data, bt->num_coef, &m)))
     goto error_return;

   if (0 != lsq_set_data (lsq, num_data, x, y))
     goto error_return;

   if (0 != lsq_solve (lsq, bt->coef, &chisqr))
     goto error_return;

   lsq_free (lsq);

   return bt;

error_return:
   bspline_free (bt);
   lsq_free (lsq);
   return NULL;
}
