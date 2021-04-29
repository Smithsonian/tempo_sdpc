#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "slit_function_asg.h"
#include "util.h"

typedef struct
{
   double *var;         /* packed with wavelength index varying fastest */
   int num_waves;
   int num_xtrack;
}
Table_Var_Type;

#define SLIT_FUNCTION_TABLE_PRIVATE_DATA \
   Table_Var_Type waves; \
   Table_Var_Type width; \
   Table_Var_Type power; \
   Table_Var_Type asym; \
   double *norm; \
   int num_params;
#include "slit_function_table.h"

static void free_table_var (Table_Var_Type *tv)
{
   FREE(tv->var);
   tv->var = NULL;
}

static void sf_table_free (SF_Table_Type *stt)
{
   if (stt == NULL)
     return;
   free_table_var (&stt->waves);
   free_table_var (&stt->width);
   free_table_var (&stt->power);
   free_table_var (&stt->asym);
   FREE(stt->norm);
   FREE(stt);
}

static int stt_close (SF_Table_Type *stt)
{
   sf_table_free (stt);
   return 0;
}

static int stt_size (const SF_Table_Type *stt, int *num_xtrack, int *num_waves, int *num_params)
{
   if (num_xtrack) *num_xtrack = stt->width.num_xtrack;
   if (num_waves) *num_waves = stt->width.num_waves;
   if (num_params) *num_params = stt->num_params;
   return 0;
}

static int compute_sf_table_norms (SF_Table_Type *stt)
{
   Table_Var_Type *tv_width = &stt->width;
   double *width = stt->width.var;
   double *power = stt->power.var;
   double *asym = stt->asym.var;
   size_t i, n = tv_width->num_xtrack * tv_width->num_waves;

   tell_vlog (TELL_MSGTYPE_INFO, 0, "computing slit function table norms [start]");

   if (NULL == (stt->norm = (double *)MALLOC (n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < n; i++)
     {
        double norm, params[3];
        params[0] = width[i];
        params[1] = power[i];
        params[2] = asym[i];
        if (0 != asg_compute_norm (params, &norm))
          {
             /* Returning a value <= 0 would eventually trigger another attempt
              * to compute this, which would again fail.  To avoid that, we must
              * return a positive value.  Since we divide by the norm, returning
              * DBL_MAX would give zero contribution from this, while returning 1.0
              * is like hoping the profile is already normalized, even though we
              * failed to compute it. */
             norm = DBL_MAX;
          }
        stt->norm[i] = norm;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 0, "computing slit function table norms [done]");

   return 0;
}

static int stt_get_params (const SF_Table_Type *stt, int xtrack, double wave0, double *params, double *norm)
{
   double *waves_x, *width_x, *power_x, *asym_x, *norm_x;
   int x_offset, num_waves;

   if (xtrack < 0 || xtrack >= stt->waves.num_xtrack)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: xtrack = %d", __func__, xtrack);
        return -1;
     }

   num_waves = stt->waves.num_waves;
   x_offset  = xtrack * num_waves;

   waves_x = stt->waves.var + x_offset;
   width_x = stt->width.var + x_offset;
   power_x = stt->power.var + x_offset;
   asym_x  = stt->asym.var  + x_offset;
   norm_x  = stt->norm + x_offset;

   if (wave0 < waves_x[0])
     {
        params[0] = width_x[0];
        params[1] = power_x[0];
        params[2] =  asym_x[0];
        *norm = norm_x[0];
     }
   else if (wave0 >= waves_x[num_waves-1])
     {
        params[0] = width_x[num_waves-1];
        params[1] = power_x[num_waves-1];
        params[2] =  asym_x[num_waves-1];
        *norm = norm_x[num_waves-1];
     }
   else
     {
        int k = bsearch_d (wave0, waves_x, num_waves);
        double frac = (waves_x[k+1] - wave0) / (waves_x[k+1] - waves_x[k]);

        params[0] = width_x[k] * frac + width_x[k+1] * (1.0 - frac);
        params[1] = power_x[k] * frac + power_x[k+1] * (1.0 - frac);
        params[2] =  asym_x[k] * frac +  asym_x[k+1] * (1.0 - frac);
        *norm     =  norm_x[k] * frac +  norm_x[k+1] * (1.0 - frac);
     }

   return 0;
}

static SF_Table_Type *sf_table_new (int num_params)
{
   SF_Table_Type *stt = NULL;

   if (NULL == (stt = (SF_Table_Type *)MALLOC (sizeof *stt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)stt, 0, sizeof *stt);

   stt->stt_close = stt_close;
   stt->stt_get_params = stt_get_params;
   stt->stt_size = stt_size;

   stt->num_params = num_params;

   return stt;
}

static int read_table_var (int grp, const char *name, Table_Var_Type *tv)
{
   TIO_Var_Info_Type info = {0};
   int start[2], count[2];

   memset ((char *)tv, 0, sizeof (*tv));

   if (0 != TIO_inq_var (grp, name, &info))
     return -1;

   tv->num_xtrack = info.dimlens[0];
   tv->num_waves = info.dimlens[1]; /* fastest varying index in the file */

   if (NULL == (tv->var = (double *)MALLOC (tv->num_waves * tv->num_xtrack * sizeof(double))))
     {
        free_table_var (tv);
        tv = NULL;
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = tv->num_xtrack;
   count[1] = tv->num_waves;

   if (0 != TIO_get_var_section (grp, name, start, count, TIO_DOUBLE, tv->var))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading variable %s", __func__, name);
        free_table_var (tv);
        tv = NULL;
        return -1;
     }

   return 0;
}

static int read_sf_table (int grp, SF_Table_Type *stt)
{
   if ((0 != read_table_var (grp, "sf_hw1e", &stt->width))
       || (0 != read_table_var (grp, "sf_shape", &stt->power))
       || (0 != read_table_var (grp, "sf_asym", &stt->asym)))
     {
        return -1;
     }

   if (0 != read_table_var (grp, "sf_wavelength", &stt->waves))
     return -1;

   return 0;
}

static int read_irr_table_var (int grp, const char *name,
                               TIO_Var_Info_Type *info, Table_Var_Type *tv)
{
   int start[3], count[3];

   memset ((char *)tv, 0, sizeof (*tv));

   if (info->ndims == 0)
     {
        /* The wavelength grid won't be explicitly defined in the irradiance file
         * so for that variable, we take the dimension sizes from a previously
         * read SF parameter. libtio will automatically construct the grid upon request.
         */
        if (0 != TIO_inq_var (grp, name, info))
          return -1;
     }

   tv->num_xtrack = info->dimlens[1];
   tv->num_waves = info->dimlens[2]; /* fastest varying index in the file */

   if (NULL == (tv->var = (double *)MALLOC (tv->num_waves * tv->num_xtrack * sizeof(double))))
     {
        free_table_var (tv);
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = tv->num_xtrack;
   count[2] = tv->num_waves;

   if (0 != TIO_get_var_section (grp, name, start, count, TIO_DOUBLE, tv->var))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading variable %s", __func__, name);
        free_table_var (tv);
        return -1;
     }

   return 0;
}

static int read_irr_sf_table (int grp, SF_Table_Type *stt)
{
   TIO_Var_Info_Type info = {0};

   if ((0 != read_irr_table_var (grp, "sf_hw1e", &info, &stt->width))
       || (0 != read_irr_table_var (grp, "sf_shape", &info, &stt->power))
       || (0 != read_irr_table_var (grp, "sf_asym", &info, &stt->asym)))
     {
        return -1;
     }

   if (0 != read_irr_table_var (grp, TEMPO_VAR_WAVELENGTH, &info, &stt->waves))
     return -1;

   return 0;
}

SF_Table_Type *sf_table_open (const char *sf_file, const char *band_name)
{
   SF_Table_Type *stt = NULL;
   int num_params = 3;
   int grp, varid, status, is_sf_table, ncid = 0;

   if (0 != TIO_open (sf_file, NC_NOWRITE, &ncid))
     return NULL;

   if ((0 != TIO_inq_grp (ncid, band_name, &grp))
       || (NULL == (stt = sf_table_new (num_params))))
     {
        TIO_close (ncid);
        return NULL;
     }

   /* The slit function may come from in an irradiance file */
   if (0 == tio_inq_varid (grp, TEMPO_VAR_IRRADIANCE, &varid))
     is_sf_table = 0;
   else
     is_sf_table = 1;

   if (is_sf_table)
     status = read_sf_table (grp, stt);
   else
     status = read_irr_sf_table (grp, stt);

   if (status == 0)
     {
        status = compute_sf_table_norms (stt);
     }

   if (status)
     {
        sf_table_free (stt);
        stt = NULL;
     }

   TIO_close (ncid);

   return stt;
}

#ifdef UNIT_TEST

typedef struct
{
   int xtrack;
   double wave0;
   const double expected_params[3];
}
SF_Test_Type;

static int dbl_cmp (double a, double b, double tol)
{
   double diff =     (a-b);
   double mean = 0.5*(a+b);
   return fabs(diff) > (fabs(mean) * tol);
}

static int test_table (const char *sf_file,
                       const char *band_name, SF_Test_Type *test_cases)
{
   SF_Table_Type *stt = NULL;
   SF_Test_Type *ptest;
   double tol = 1.e-4;
   int status = -1;
   int num_errors, nw, nx, np;

   if (NULL == (stt = sf_table_open (sf_file, band_name)))
     return -1;

   if (0 != stt->stt_size (stt, &nx, &nw, &np))
     goto return_status;

   if ((nx != 2048) || (nw != 1028) || (np != 3))
     {
        fprintf (stderr, "*** Warning: unexpected table dimensions: xtrack=%d wave=%d params=%d\n",
                 nx, nw, np);
     }

   num_errors = 0;
   for (ptest = test_cases; (ptest != NULL) && (ptest->xtrack >= 0); ptest++)
     {
        double params[3];
        int i;
        if (0 != stt->stt_get_params (stt, ptest->xtrack, ptest->wave0, params))
          goto return_status;
        for (i = 0; i < 3; i++)
          {
             if (dbl_cmp(params[i], ptest->expected_params[i], tol))
               {
                  fprintf (stderr, "*** Warning: params value mismatch: expected p[%d] = %f got %f (tolerance=%e)\n",
                           i, ptest->expected_params[i], params[i], tol);
                  num_errors++;
                  continue;
               }
          }
     }

   if (num_errors == 0) status = 0;
return_status:
   if (stt) stt->stt_close (stt);
   return status;
}

int main (void)
{
   const char *sf_file = "/data/tempo/sdpc/refdata/instrument/TEMPO_Slit_Function_Nischal10292019_V1.nc";
   SF_Test_Type vis_tests[] =
     {
        {   0, 537.2, {0.30453, 3.4585,  0.0059415}},
        {2047, 537.2, {0.31403, 3.4443, -0.0025236}},
        {2047, 742.6, {0.28651, 2.4118, -0.0027844}},
        {-1,-1,{0,0,0}}
     };
   SF_Test_Type uv_tests[] =
     {
        {   0, 291.8, {0.30254, 3.8627, 7.46e-05}},
        {2047, 291.8, {0.30658, 3.6683, -0.00445}},
        {2047, 497.2, {0.30796, 3.5532, -0.00308}},
        {-1,-1,{0,0,0}}
     };

   if (0 != test_table (sf_file, TEMPO_BAND_NAME_VIS, vis_tests))
     return EXIT_FAILURE;

   if (0 != test_table (sf_file, TEMPO_BAND_NAME_UV, uv_tests))
     return EXIT_FAILURE;

   return EXIT_SUCCESS;
}

#endif
