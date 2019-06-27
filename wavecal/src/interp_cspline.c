#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <tell.h>

#include "config.h"

#include "interp_cspline.h"
#include "lsq.h"

struct Cspline_Type
{
   gsl_spline *spline;
   gsl_interp_accel *acc;
   size_t n;
};

void cspline_free (Cspline_Type *ct)
{
   if (ct == NULL)
     return;
   gsl_interp_accel_free (ct->acc);
   gsl_spline_free (ct->spline);
   FREE(ct);
}

static Cspline_Type *cspline_alloc (size_t n)
{
   Cspline_Type *ct = NULL;

   if (NULL == (ct = (Cspline_Type *)MALLOC (sizeof *ct)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)ct, 0, sizeof *ct);

   ct->n = n;

   if ((NULL == (ct->acc = gsl_interp_accel_alloc()))
       || (NULL == (ct->spline = gsl_spline_alloc (gsl_interp_cspline, ct->n))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        cspline_free (ct);
        return NULL;
     }

   return ct;
}

int cspline_eval (Cspline_Type *ct,
                  size_t n, const double *x, double *yest)
{
   size_t i;
   int status;

   for (i = 0; i < n; i++)
     {
        status = gsl_spline_eval_e (ct->spline, x[i], ct->acc, &yest[i]);
        if (status)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: %s", __func__, gsl_strerror(status));
             return -1;
          }
     }

   return 0;
}

Cspline_Type *cspline_interpol (size_t num_data, const double *x, const double *y,
                                void *config)
{
   Cspline_Type *ct = NULL;
   int status;

   (void) config;

   if (NULL == (ct = cspline_alloc (num_data)))
     return NULL;

   if (0 != (status = gsl_spline_init (ct->spline, x, y, num_data)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: %s", __func__, gsl_strerror(status));
        cspline_free (ct);
        return NULL;
     }

   return ct;
}
