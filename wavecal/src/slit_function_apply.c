#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <tell.h>

#include "slit_function.h"
#include "slit_function_asg.h"

struct Slit_Function_Type
{
   /* function to compute slit-function value and parameter derivatives */
   int (*sf_eval)(const double *x, size_t nx, double *params,
                  double *value, double *param_step, double *param_derivs[SFT_NUM_PARAMS]);

   double param_step[SFT_NUM_PARAMS];  /* delta step to use in compute derivs w.r.t. each parameter */
   double *params;                     /* parameter array [num_waves, k], param index k varies fastest */
   int num_waves;

   int num_sf;    /* number of equally spaced points at which slit-function will be evaluated */
   double dx;     /* spacing of slit-function eval points */
   double *x;     /* slit-function eval points */
   double *sf;    /* slit-function values */

   double *sf_derivs[SFT_NUM_PARAMS];  /* slit-function derivative w.r.t. each parameter */
};

void sft_free (Slit_Function_Type *sft)
{
   int i;

   if (sft == NULL)
     return;

   FREE(sft->x);
   FREE(sft->sf);
   for (i = 0; i < SFT_NUM_PARAMS; i++)
     {
        FREE(sft->sf_derivs[i]);
        sft->sf_derivs[i] = NULL;
     }
   FREE(sft);
}

Slit_Function_Type *sft_init (double dx, size_t num_sf, size_t num_waves,
                              double *params, double *param_step)
{
   Slit_Function_Type *sft = NULL;
   int j, nj, m;

   if (NULL == (sft = (Slit_Function_Type *)MALLOC (sizeof *sft)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sft, 0, sizeof(*sft));

   if ((NULL == (sft->x = (double *)MALLOC (num_sf * sizeof(double))))
       || (NULL == (sft->sf = (double *)MALLOC (num_sf * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        sft_free (sft);
        return NULL;
     }

   /* Assume spectrum has padding of length, m=num_sf/2 at
    * both ends, so the first real spectrum point is at spec[m],
    * and we can access spec[num_waves + m - 1].
    */

   m = num_sf/2;
   nj = num_sf;

   /* x[j] = wavelength grid relative to the slit-function center */
   for (j = 0; j < nj; j++)
     {
        sft->x[j] = (j-m)*dx;
     }

   memset ((char *)sft->sf, 0, num_sf * sizeof(double));
   memcpy ((char *)sft->param_step, (char *)param_step, SFT_NUM_PARAMS * sizeof(double));

   for (j = 0; j < SFT_NUM_PARAMS; j++)
     {
        if (NULL == (sft->sf_derivs[j] = (double *)MALLOC (num_sf * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             sft_free (sft);
             return NULL;
          }
        memset ((char *)sft->sf_derivs[j], 0, num_sf * sizeof(double));
     }

   sft->dx = dx;
   sft->num_sf = num_sf;
   sft->num_waves = num_waves;
   sft->params = params;

   sft->sf_eval = asg_normed_plus_derivs;

   return sft;
}

static int cached_params_differ (const double *p0, const double *p, int n, double tol)
{
   int i;

   for (i = 0; i < n; i++)
     {
        double diff = (p[i] - p0[i]);
        double mid  = (p[i] + p0[i]) * 0.5;
        if (fabs(diff) > tol * fabs(mid))
          return 1;
     }

   return 0;
}

/* When the slit-function shape is fixed, one can use FFTs to compute the
 * convolution more efficiently.  Because the slit-function shape may vary
 * we're probably stuck with direct integration.  To reduce the cost, we
 * cache the slit function, and re-evaluate it only when any parameter change
 * exceeds some tolerance.
 */
int sft_apply (Slit_Function_Type *sft, const double *spec, double *spec_convolved,
               int compute_derivs, double *spec_derivs_convolved[SFT_NUM_PARAMS])
{
   double prev_par[SFT_NUM_PARAMS];
   int k, m, nsf;

   nsf = sft->num_sf;
   m   = nsf / 2;

   memset ((char *)prev_par, 0, SFT_NUM_PARAMS * sizeof(double));

   for (k = m; k < sft->num_waves + m; k++)
     {
        double *par = sft->params + (k-m) * SFT_NUM_PARAMS;
        double s;
        int j;

        if (cached_params_differ (par, prev_par, SFT_NUM_PARAMS, DBL_EPSILON))
          {
             if (0 != sft->sf_eval (sft->x, sft->num_sf, par, sft->sf, sft->param_step, sft->sf_derivs))
               return -1;

             memcpy ((char *)prev_par, (char *)par, SFT_NUM_PARAMS * sizeof(double));
          }

        /* Evaluate convolution integral using trapezoid rule
         * for uniform grid spacing */
        s = 0.5 * (spec[k-m] * sft->sf[nsf-1] + spec[k+m-1] * sft->sf[0]);
        for (j = 1; j < nsf-1; j++)
          {
             s += spec[k-m+j] * sft->sf[nsf-j-1];
          }
        spec_convolved[k-m] = s * sft->dx;

        if (compute_derivs)
          {
             int n;
             for (n = 0; n < SFT_NUM_PARAMS; n++)
               {
                  double *deriv = sft->sf_derivs[n];
                  double *deriv_convolved = spec_derivs_convolved[n];

                  if (deriv_convolved == NULL)
                    continue;

                  /* Evaluate convolution integral using trapezoid rule
                   * for uniform grid spacing */
                  s = 0.5 * (spec[k-m]*deriv[nsf-1] + spec[k+m-1]*deriv[0]);
                  for (j = 1; j < nsf-1; j++)
                    {
                       s += spec[k-m+j] * deriv[nsf-j-1];
                    }
                  deriv_convolved[k-m] = s * sft->dx;
               }
          }
     }

   return 0;
}

#ifdef UNIT_TEST

#include <gsl/gsl_randist.h>

int main (void)
{
   Slit_Function_Type *sft = NULL;
   double params0[SFT_NUM_PARAMS] = {0.25, 2.0, 0.0};
   double param_step[SFT_NUM_PARAMS] = {1.e-4, 1.e-4, 1.e-4};
   double *tmp = NULL;
   double *params = NULL;
   double *spec_padded = NULL;
   double *spec_convolved = NULL;
   double *spec_derivs_convolved[3] = {NULL, NULL, NULL};
   double dx = 0.02;
   int num_sf = 12 * params0[0]/dx;
   int num_waves = num_sf * 2;
   int compute_derivs = 1;
   size_t offset;
   int i, i0, m, len_tmp;
   int status = -1;

   if (NULL == (params = (double *)MALLOC (num_waves * 3 * sizeof(double))))
     {
        fprintf (stderr, "%s: malloc failed", __func__);
        goto return_status;
     }

   for (i = 0; i < num_waves; i++)
     {
        double *par = params + i*3;
        memcpy ((char *)par, (char *)params0, 3 * sizeof(double));
     }

   len_tmp = 5*num_waves + num_sf;
   if (NULL == (tmp = (double *)MALLOC (len_tmp * sizeof(double))))
     {
        fprintf (stderr, "%s: malloc failed", __func__);
        goto return_status;
     }
   memset ((char *)tmp, 0, len_tmp * sizeof(double));

   offset = 0;
   spec_padded = tmp + offset;
   offset += num_waves + num_sf;

   spec_convolved = tmp + offset;
   offset += num_waves;

   for (i = 0; i < SFT_NUM_PARAMS; i++)
     {
        spec_derivs_convolved[i] = tmp + offset;
        offset += num_waves;
     }

   m = num_sf/2;
   i0 = num_waves/2;

   for (i = 0; i < num_sf; i++)
     {
        double x = (i-m) * dx;
        spec_padded[i0 + i] = gsl_ran_gaussian_pdf (x, params0[0]/sqrt(2));
     }

   if (NULL == (sft = sft_init (dx, num_sf, num_waves, params, param_step)))
     goto return_status;

   if (0 != sft_apply (sft, spec_padded, spec_convolved, compute_derivs, spec_derivs_convolved))
     goto return_status;

   /* In terms of our width parameter, w, the gaussian sigma = w/sqrt(2),
    * so when the target is a gaussian, sigma1,and we convolve with gaussian,
    * sigma2, the convolved result should have:
    *   sigma = hypot(sigma1,sigma2)
    *         = sqrt(0.5*w^2, 0.5*2^2)
    *         = w
    */

   fprintf (stdout, "# i s conv(s) conv(s_expect) conv(d0) conv(d1) conf(d2)\n");
   for (i = 0; i < num_waves; i++)
     {
        fprintf (stdout, "%4d %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n",
                 i, spec_padded[i+m], spec_convolved[i],
                 gsl_ran_gaussian_pdf ((i-i0-1)*dx, params0[0]),
                 spec_derivs_convolved[0][i],
                 spec_derivs_convolved[1][i],
                 spec_derivs_convolved[2][i]);
     }

   status = 0;
return_status:
   FREE(params);
   FREE(tmp);
   sft_free (sft);

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
#endif
