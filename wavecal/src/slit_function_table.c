#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

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
   int num_params;
#include "slit_function_table.h"

static void free_table_var (Table_Var_Type *tv)
{
   FREE(tv->var);
}

static void sf_table_free (SF_Table_Type *stt)
{
   if (stt == NULL)
     return;
   free_table_var (&stt->waves);
   free_table_var (&stt->width);
   free_table_var (&stt->power);
   free_table_var (&stt->asym);
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

static int stt_get_params (const SF_Table_Type *stt, int xtrack, double wave0, double *params)
{
   double *waves_x, *width_x, *power_x, *asym_x;
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

   if (wave0 < waves_x[0])
     {
        params[0] = width_x[0];
        params[1] = power_x[0];
        params[2] =  asym_x[0];
     }
   else if (wave0 >= waves_x[num_waves-1])
     {
        params[0] = width_x[num_waves-1];
        params[1] = power_x[num_waves-1];
        params[2] =  asym_x[num_waves-1];
     }
   else
     {
        int k = bsearch_d (wave0, waves_x, num_waves);
        double frac = (waves_x[k+1] - wave0) / (waves_x[k+1] - waves_x[k]);

        params[0] = width_x[k] * frac + width_x[k+1] * (1.0 - frac);
        params[1] = power_x[k] * frac + power_x[k+1] * (1.0 - frac);
        params[2] =  asym_x[k] * frac +  asym_x[k+1] * (1.0 - frac);
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

static int read_table_var (int ncid, const char *band_name, const char *name,
                           Table_Var_Type *tv)
{
   TIO_Var_Info_Type info = {0};
   int start[2], count[2];
   double *vars_for_wave_i = NULL;
   int num_xtrack, num_waves;
   int beg_row, end_row, k, i;
   size_t len;

   memset ((char *)tv, 0, sizeof (*tv));

   if (0 != TIO_inq_var (ncid, name, &info))
     return -1;

   num_waves = info.dimlens[0];
   num_xtrack = info.dimlens[1];  /* fastest varying index in the file */

   /* Parameters are are stored like so:
    * -----------------------------  290 nm: wave_index=2056-1
    * |             |             |
    * |             |             |
    * |   D         |   C         |    <- UV
    * |             |             |
    * |             |             |  490 nm
    * -----------------------------
    * |             |             |  540 nm
    * |             |             |
    * |   A         |   B         |    <- VIS
    * |             |             |
    * |             |             |
    * -----------------------------  740 nm: wave_index=0
    * north                     south
    * xtrack=0                  xtrack=2048-1
    */
   if (0 == strcmp (band_name, TEMPO_BAND_NAME_VIS))
     {
        beg_row = 0;
        end_row = num_waves/2;
     }
   else if (0 == strcmp (band_name, TEMPO_BAND_NAME_UV))
     {
        beg_row = num_waves/2;
        end_row = num_waves;
     }
   else
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: unrecognized band name: %s", __func__, band_name);
        return -1;
     }

   tv->num_waves = num_waves/2;
   tv->num_xtrack = num_xtrack;

   len = tv->num_waves * tv->num_xtrack;

   if ((NULL == (tv->var = (double *)MALLOC (len * sizeof(double))))
       || (NULL == (vars_for_wave_i = (double *)MALLOC (num_xtrack * sizeof(double)))))
     {
        free_table_var (tv);
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   /* When using the parameters, we want to use a wavelength index that goes
    * from 0->(N-1) within a single band, with wavelengths in increasing order.
    * To achieve this, we re-order the array on input.
    */

   count[0] = 1;
   start[1] = 0;
   count[1] = num_xtrack;

   for (k=beg_row, i=tv->num_waves-1; (k<end_row) && (i >= 0); k++, i--)
     {
        int j;
        start[0] = k;
        if (0 != TIO_get_var_section (ncid, name, start, count, TIO_DOUBLE, vars_for_wave_i))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading variable %s", __func__, name);
             free_table_var (tv);
             FREE(vars_for_wave_i);
             return -1;
          }

        /* tv->var has dimensions [xtrack, wave] with wave varying fastest so */
        for (j = 0; j < num_xtrack; j++)
          {
             tv->var[i + j * tv->num_waves] = vars_for_wave_i[j];
          }
     }

   FREE(vars_for_wave_i);

   return 0;
}

SF_Table_Type *sf_table_open (const char *sf_file, const char *cal_file,
                              const char *band_name)
{
   SF_Table_Type *stt = NULL;
   int num_params = 3;
   int ncid = 0;

   if (NULL == (stt = sf_table_new (num_params)))
     return NULL;

   if (0 != TIO_open (sf_file, NC_NOWRITE, &ncid))
     {
        sf_table_free (stt);
        return NULL;
     }

   if ((0 != read_table_var (ncid, band_name, "W", &stt->width))
       || (0 != read_table_var (ncid, band_name, "K", &stt->power))
       || (0 != read_table_var (ncid, band_name, "A", &stt->asym)))
     {
        sf_table_free (stt);
        TIO_close (ncid);
        return NULL;
     }

   if (0 != TIO_open (cal_file, NC_NOWRITE, &ncid))
     {
        sf_table_free (stt);
        return NULL;
     }
   if (0 != read_table_var (ncid, band_name, "wavelength_grid", &stt->waves))
     {
        sf_table_free (stt);
        TIO_close (ncid);
        return NULL;
     }

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

static int test_table (const char *sf_file, const char *cal_file,
                       const char *band_name, SF_Test_Type *test_cases)
{
   SF_Table_Type *stt = NULL;
   SF_Test_Type *ptest;
   double tol = 1.e-4;
   int status = -1;
   int num_errors, nw, nx, np;

   if (NULL == (stt = sf_table_open (sf_file, cal_file, band_name)))
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
   const char *sf_file = "/data/tempo/sdpc/refdata/instrument/TEMPO_Slit_Function_10_29_2019_V1.h5";
   const char *cal_file = "/data/tempo/sdpc/refdata/instrument/TEMPO_senscal_larc_update_V2_04242019_v1.1.nc";
   SF_Test_Type vis_tests[] =
     {
        {   0, 537.2, {0.30453, 3.4585,  0.0059415}},  /* e.g.  h5dump -d W -s 1027,0 -c 1,2048 $FILE */
        {2047, 537.2, {0.31403, 3.4443, -0.0025236}},  /* e.g.  h5dump -d W -s 1027,0 -c 1,2048 $FILE */
        {2047, 742.6, {0.28651, 2.4118, -0.0027844}},  /* e.g.  h5dump -d W -s 0,0    -c 1,2048 $FILE */
        {-1,-1,{0,0,0}}
     };
   SF_Test_Type uv_tests[] =
     {
        {   0, 291.8, {0.30254, 3.8627, 7.46e-05}},  /* e.g.  h5dump -d W -s 2055,0 -c 1,2048 $FILE */
        {2047, 291.8, {0.30658, 3.6683, -0.00445}},  /* e.g.  h5dump -d W -s 2055,0 -c 1,2048 $FILE */
        {2047, 497.2, {0.30796, 3.5532, -0.00308}},  /* e.g.  h5dump -d W -s 1028,0 -c 1,2048 $FILE */
        {-1,-1,{0,0,0}}
     };

   if (0 != test_table (sf_file, cal_file, TEMPO_BAND_NAME_VIS, vis_tests))
     return EXIT_FAILURE;

   if (0 != test_table (sf_file, cal_file, TEMPO_BAND_NAME_UV, uv_tests))
     return EXIT_FAILURE;

   return EXIT_SUCCESS;
}

#endif
