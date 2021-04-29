#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <tell.h>

#include "slit_function_asg.h"

/* Asymmetric super-Gaussian (asg)
 * Parameterization following Beirle et al 2017, Atmos. Meas. Tech., 10, 581 */
static inline double asg_exp (double x, void *params)
{
   double *par = (double *)params;
   double w   = par[0];   /* half-width at 1/e */
   double k   = par[1];   /* exponent of exp arg*/
   double a_w = par[2];   /* width asymmetry */
   double s = (x > 0) ? (w + a_w) : (w - a_w);

   if (s == 0.0)
     return 0.0;

   return exp (-pow(fabs(x/s), k));
}

typedef struct
{
   double norm;
   double params[3];
   gsl_integration_workspace *work;
   size_t limit;
}
Norm_Cache_Type;

static Norm_Cache_Type Norm_Cache =
{
   GSL_NAN,
   {GSL_POSINF, GSL_POSINF, GSL_POSINF},
   NULL,
   61*32   /* 32 copies of the 61-point Gauss-Kronod rule */
};

int asg_compute_norm (double *params, double *norm)
{
   Norm_Cache_Type *nc = &Norm_Cache;
   gsl_function fptr = {0};
   int status, compute_norm = 0;
   double a, b, epsabs, epsrel, abserr, w;

   if (nc->work == NULL)
     {
        compute_norm++;
        /* Allocate global work space on the first call,
         * and let it be freed when the process exits (ugh) */
        if (NULL == (nc->work = gsl_integration_workspace_alloc (nc->limit)))
          return -1;
     }
   else
     {
        int i;
        for (i = 0; i < 3; i++)
          {
             if (0 != gsl_fcmp (nc->params[i], params[i], DBL_EPSILON))
               {
                  compute_norm++;
                  break;
               }
          }
     }

   if (compute_norm == 0)
     {
        *norm = nc->norm;
        return 0;
     }

   memcpy ((char *)nc->params, (char *)params, 3 * sizeof(double));

   fptr.function = asg_exp;
   fptr.params = params;

   w = params[0];  /* half-width at 1/e */

   /* 2*w = full-width at 1/e; area is of order 2w * exp(0) = 2w.
    * => scale desired absolute error tolerance by 2w.
    */

   epsabs = 1.e-4 * (2*w);
   epsrel = 1.e-6;
   abserr = 0.0;
   a = -6 * w;
   b = +6 * w;

   /* Non-adaptive integration sometimes fails to meet the required accuracy */
   status = gsl_integration_qag (&fptr, a, b, epsabs, epsrel,
                                 nc->limit, GSL_INTEG_GAUSS51, nc->work,
                                 norm, &abserr);
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: GSL error: %s",
                     __func__, gsl_strerror (status));
        return -1;
     }

   if (abserr > 1.e-2 * (2*w))
     {
        tell_vwarn (0, "%s: WARNING: abserr = %g (norm=%g)",
                    __func__, abserr, *norm);
     }

   nc->norm = *norm;

   return 0;
}

static int asg_normed (const double *x, size_t nx, double *params, double asg_area, double *value)
{
   size_t i;

   /* ASG normalized to have unit area within finite integration limits.
    * If we're given a valid norm, use it.  Otherwise, compute it.
    */
   if (asg_area < 0.0)
     {
        if (0 != asg_compute_norm (params, &asg_area))
          return -1;
     }

   for (i = 0; i < nx; i++)
     {
        value[i] = asg_exp (x[i], params) / asg_area;
     }

   return 0;
}

int asg_normed_plus_derivs (const double *x, size_t nx, double *params, double norm,
                            double *value,
                            double *param_step,
                            double *param_derivs[3])
{
   int j;

   if (0 != asg_normed (x, nx, params, norm, value))
     return -1;

   /* optionally compute derivatives */
   if ((param_step == NULL) || (param_derivs == NULL))
     return 0;

   for (j = 0; j < 3; j++)
     {
        double *deriv = param_derivs[j];
        double par_copy[3];
        size_t i;

        if (deriv == NULL)
          continue;

        memcpy ((char *)par_copy, (char *)params, 3 * sizeof(double));
        par_copy[j] += param_step[j];

        if (0 != asg_normed (x, nx, par_copy, norm, deriv))
          return -1;

        for (i = 0; i < nx; i++)
          {
             deriv[i] = (deriv[i] - value[i]) / param_step[j];
          }
     }

   return 0;
}

#ifdef UNIT_TEST

#define NX 1025

int main (void)
{
   double params[3] = {1.0, 2, 0.0};
   double param_step[3] = {1.e-4, 1.e-3, 1.e-5};
   double x[NX], value[NX];
   double dvdp0[NX], dvdp1[NX], dvdp2[NX];
   double *derivs[3];
   double xb = -3.0;
   double xe = +3.0;
   double norm = -1.0;
   size_t i, nx = NX;

   for (i = 0; i < nx; i++)
     {
        x[i] = xb + i * (xe-xb) / (nx - 1);
     }

   derivs[0] = dvdp0;
   derivs[1] = dvdp1;
   derivs[2] = dvdp2;

   if (0 != asg_normed_plus_derivs (x, nx, params, norm, value, param_step, derivs))
     return 1;

   for (i = 0; i < nx; i++)
     {
        fprintf (stdout, "%22.14e %22.14e %22.14e %22.14e %22.14e\n",
                 x[i], value[i], derivs[0][i], derivs[1][i], derivs[2][i]);
     }

   return 0;
}
#endif
