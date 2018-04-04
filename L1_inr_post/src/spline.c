/** @file spline.c
 *  @brief Cubic spline implementation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include <tell.h>

#include "spline.h"

int spline (const double *x, const double *y, int nx,
            const double *xs, double *ys, int nxs)
{
   double nan_value = nan("");
   gsl_interp *interp = NULL;
   gsl_interp_accel *acc = NULL;
   size_t n = nx;
   size_t i, ns = nxs;

   if ((NULL == (interp = gsl_interp_alloc (gsl_interp_cspline, n)))
       || (NULL == (acc = gsl_interp_accel_alloc ())))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gsl_interp alloc failed", __func__);
        gsl_interp_free (interp);
        return -1;
     }

   if (gsl_interp_init (interp, x, y, n))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gsl_interp init failed", __func__);
        gsl_interp_free (interp);
        gsl_interp_accel_free (acc);
        return -1;
     }

   for (i = 0; i < ns; i++)
     {
        if (gsl_interp_eval_e (interp, x, y, xs[i], acc, &ys[i]))
          ys[i] = nan_value;
     }

   gsl_interp_accel_free (acc);
   gsl_interp_free (interp);

   return 0;
}
