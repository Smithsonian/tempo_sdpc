#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <tell.h>

#include "slit_function.h"

struct Slit_Function_Type
{
   SFT_Eval_Type *sf_eval;
   /**< function to compute slit-function value and parameter derivatives */

   double param_step[SFT_NUM_PARAMS];
   /**< delta step to use in compute derivs w.r.t. each parameter */

   int num_sf;    /**< number of equally spaced points at which slit-function will be evaluated */
   double dx;     /**< spacing of slit-function eval points */
   double *x;     /**< slit-function eval points */
   double *sf;    /**< slit-function values */

   double *sf_derivs[SFT_NUM_PARAMS];
   /**< slit-function derivative w.r.t. each parameter */
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

Slit_Function_Type *sft_new (int num_sf)
{
   Slit_Function_Type *sft = NULL;
   int i;

   if (NULL == (sft = (Slit_Function_Type *)MALLOC (sizeof *sft)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sft, 0, sizeof(*sft));

   sft->num_sf = num_sf;

   if ((NULL == (sft->x = (double *)MALLOC (num_sf * sizeof(double))))
       || (NULL == (sft->sf = (double *)MALLOC (num_sf * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        sft_free (sft);
        return NULL;
     }

   memset ((char *)sft->x, 0, num_sf * sizeof(double));
   memset ((char *)sft->sf, 0, num_sf * sizeof(double));

   for (i = 0; i < SFT_NUM_PARAMS; i++)
     {
        if (NULL == (sft->sf_derivs[i] = (double *)MALLOC (num_sf * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             sft_free (sft);
             return NULL;
          }
        memset ((char *)sft->sf_derivs[i], 0, num_sf * sizeof(double));
     }

   return sft;
}

int sft_config (Slit_Function_Type *sft, SFT_Eval_Type *sf_eval,
                double dx, double *param_step)
{
   int j, m;

   sft->sf_eval = sf_eval;

   sft->dx = dx;

   m = sft->num_sf/2;

   /* x[j] = wavelength grid relative to the slit-function center */
   for (j = 0; j < sft->num_sf; j++)
     {
        sft->x[j] = (j-m)*dx;
     }

   memcpy ((char *)sft->param_step, (char *)param_step, SFT_NUM_PARAMS * sizeof(double));

   return 0;
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
int sft_apply (Slit_Function_Type *sft, SFT_Param_Type *sf_params, void *cl,
               int num_waves, const double *spec_padded, double *spec_convolved,
               double *spec_derivs_convolved[SFT_NUM_PARAMS])
{
   double prev_par[SFT_NUM_PARAMS];
   int k, m, nsf;

   /* Assume spec_padded has padding of length, m=num_sf/2 at
    * each end, so the first real spectrum point is at spec[m],
    * and we can access spec[num_waves + m - 1].
    */

   nsf = sft->num_sf;
   m   = nsf / 2;

   memset ((char *)prev_par, 0, SFT_NUM_PARAMS * sizeof(double));

   for (k = m; k < num_waves + m; k++)
     {
        double s, par[3];
        int j;

        if (0 != sf_params (k-m, SFT_NUM_PARAMS, par, cl))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: retrieving parameters for wavelength index = %d",
                          __func__, k-m);
             return -1;
          }

        if (cached_params_differ (par, prev_par, SFT_NUM_PARAMS, DBL_EPSILON))
          {
             if (0 != sft->sf_eval (sft->x, sft->num_sf, par, sft->sf, sft->param_step, sft->sf_derivs))
               return -1;

             memcpy ((char *)prev_par, (char *)par, SFT_NUM_PARAMS * sizeof(double));
          }

        /* Evaluate convolution integral using trapezoid rule
         * for uniform grid spacing */
        s = 0.5 * (spec_padded[k-m] * sft->sf[nsf-1] + spec_padded[k+m-1] * sft->sf[0]);
        for (j = 1; j < nsf-1; j++)
          {
             s += spec_padded[k-m+j] * sft->sf[nsf-j-1];
          }
        spec_convolved[k-m] = s * sft->dx;

        if (spec_derivs_convolved)
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
                  s = 0.5 * (spec_padded[k-m]*deriv[nsf-1] + spec_padded[k+m-1]*deriv[0]);
                  for (j = 1; j < nsf-1; j++)
                    {
                       s += spec_padded[k-m+j] * deriv[nsf-j-1];
                    }
                  deriv_convolved[k-m] = s * sft->dx;
               }
          }
     }

   return 0;
}

#ifdef UNIT_TEST

#include <gsl/gsl_randist.h>
#include "slit_function_asg.h"

static int get_params (int k, int num_pars, double *pars, void *cl)
{
   double params0[SFT_NUM_PARAMS] = {0.25, 2.0, 0.0};

   (void) k; (void) cl;

   if (num_pars != SFT_NUM_PARAMS)
     {
        fprintf (stderr, "%s: unexpected number of parameters requested (num_pars=%d, expected %d)\n",
                 __func__, num_pars, SFT_NUM_PARAMS);
        return -1;
     }
   memcpy ((char *)pars, (char *)params0, 3 * sizeof(double));
   return 0;
}

int main (void)
{
   Slit_Function_Type *sft = NULL;
   double param_step[SFT_NUM_PARAMS] = {1.e-4, 1.e-4, 1.e-4};
   double pars[3];
   double *tmp = NULL;
   double *spec_padded = NULL;
   double *spec_convolved = NULL;
   double *spec_derivs_convolved[3] = {NULL, NULL, NULL};
   double dx = 0.02;
   int num_sf, num_waves;
   size_t offset;
   int i, i0, m, len_tmp;
   int status = -1;

   (void) get_params (0, SFT_NUM_PARAMS, pars, NULL);
   num_sf = 12 * pars[0]/dx;
   num_waves = num_sf * 2;

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
        spec_padded[i0 + i] = gsl_ran_gaussian_pdf (x, pars[0]/sqrt(2));
     }

   if (NULL == (sft = sft_new (num_sf)))
     goto return_status;

   if (0 != sft_config (sft, asg_normed_plus_derivs, dx, param_step))
     goto return_status;

   if (0 != sft_apply (sft, get_params, NULL, num_waves, spec_padded,
                       spec_convolved, spec_derivs_convolved))
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
                 gsl_ran_gaussian_pdf ((i-i0-1)*dx, pars[0]),
                 spec_derivs_convolved[0][i],
                 spec_derivs_convolved[1][i],
                 spec_derivs_convolved[2][i]);
     }

   status = 0;
return_status:
   FREE(tmp);
   sft_free (sft);

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
#endif
