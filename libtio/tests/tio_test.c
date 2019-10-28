/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "netcdf.h"
#include "tio.h"
#include "_tio.h"
#include "tio_template.h"

static float *generate_data (int n) /*{{{*/
{
   int i;
   float *y;

   if (NULL == (y = (float *) malloc (n * sizeof(*y))))
     {
        fprintf (stderr, "*** malloc failed\n");
        return NULL;
     }

   for (i = 0; i < n; i++)
     {
        y[i] = (float)i;
     }

   return y;
}

/*}}}*/

static float *generate_err (int n, float err) /*{{{*/
{
   int i;
   float *y;

   if (NULL == (y = (float *) malloc (n * sizeof(*y))))
     {
        fprintf (stderr, "*** malloc failed\n");
        return NULL;
     }

   for (i = 0; i < n; i++)
     {
        y[i] = err;
     }

   return y;
}

/*}}}*/

static int compare_data (int n, float *out, float *in) /*{{{*/
{
   int i;

   for (i = 0; i < n; i++)
     {
        if (isnan(out[i]) && (0 == isnan(in[i])))
          break;

        if (out[i] == 0.0 && abs(in[i]) > FLT_EPSILON)
          break;

        if (abs(in[i] - out[i]) > FLT_EPSILON * fabs(in[i]))
          break;
     }

   if (i < n)
     {
        fprintf (stderr, "*** error:  data[%d] in=%15.10e != out=%15.10e\n",
                 i, in[i], out[i]);
        return -1;
     }

   return 0;
}

/*}}}*/

static int test_def_grp (int ncid) /*{{{*/
{
   int grp, ignore_grp;

   if (-1 == TIO_def_grp (ncid, "xxx", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }

   if ((-1 == TIO_inq_grp (ncid, "xxx", &grp))
       || (grp != ignore_grp))
     {
        fprintf (stderr, "*** TIO_inq_grp failed\n");
        return -1;
     }
   if (-1 == TIO_def_grp (ncid, "xxx/a/b/c", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }
   /* duplicate and trailing slashes are ignored */
   if (-1 == TIO_def_grp (ncid, "/xxx/a//qqq/r/", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }
   /* no-op if path already exists */
   if (-1 == TIO_def_grp (ncid, "/xxx/a/qqq/r", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }
   if (-1 == TIO_def_grp (ncid, "/xxx/a/zzz//", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }
   if (-1 == TIO_def_grp (ncid, "/", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }
   if (-1 == TIO_def_grp (ncid, "///", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
        return -1;
     }

   return 0;
}

/*}}}*/

static int test_dims (int ncid) /*{{{*/
{
   const char dimname[] = "test_dim";
   size_t len, test_dimlen = 128;
   int test_dimid, id;
   char buf[TIO_MAX_NAME_LEN];

   if (-1 == TIO_def_dim (ncid, dimname, test_dimlen, &test_dimid))
     return -1;

   if (-1 == TIO_inq_dimid (ncid, dimname, &id))
     return -1;
   if (id != test_dimid)
     {
        fprintf (stderr, "*** TIO_inq_dimid failed\n");
        return -1;
     }

   if (-1 == TIO_inq_dimname (ncid, test_dimid, buf))
     return -1;
   if (0 != strcmp (buf, dimname))
     {
        fprintf (stderr, "*** TIO_inq_dimname failed\n");
        return -1;
     }

   if (-1 == TIO_inq_dim (ncid, dimname, &id, &len))
     return -1;
   if ((id != test_dimid) || (len != test_dimlen))
     {
        fprintf (stderr, "*** TIO_inq_dim failed\n");
        return -1;
     }

   return 0;
}

/*}}}*/

static int dontcopy_attr (const char *attr)
{
   return (0 == strcmp (attr, "_FillValue"));
}

static int test_def_var (int ncid, const char *name, int type) /*{{{*/
{
   TIO_Var_Info_Type vi, vi2;
   static TIO_Attr_Text_Type attrs[] =
     {
        {"test_attr1", "This is an attribute test"},
        {"test_attr2", "This is another attribute test"},
        {NULL,NULL}
     };
   const char test_name[] = "test_var";
   size_t chunksizes[TIO_MAX_VAR_DIMS];
   int i, test_id, dimids_ok;

   if (-1 == TIO_inq_var (ncid, name, &vi))
     return -1;

   if (-1 == TIO_def_var (ncid, test_name, type, vi.ndims, vi.dimids, &test_id))
     return -1;

   if (-1 == TIO_put_text_attrs (ncid, test_id, attrs))
     return -1;

   if (-1 == TIO_inq_var (ncid, name, &vi2))
     return -1;
   dimids_ok = 1;
   for (i = 0; i < vi.ndims; i++)
     {
        chunksizes[i] = vi.dimlens[i]/2;
        if (chunksizes[i] == 0) chunksizes[i] = 1;

        if (vi.dimids[i] != vi2.dimids[i])
          {
             dimids_ok = 0;
             break;
          }
     }
   if ((vi2.ndims != vi.ndims) || (dimids_ok == 0))
     {
        fprintf (stderr, "*** ERROR: TIO_def_var/TIO_inq_var are inconsistent!\n");
        return -1;
     }

   if (-1 == TIO_def_var_deflate (ncid, test_id, 1, 1, 1))
     return -1;
   if (-1 == TIO_def_var_chunking (ncid, test_id, NC_CHUNKED, chunksizes))
     return -1;

   if (-1 == TIO_def_var_fill (ncid, test_id, 1, NULL))
     return -1;

   if (-1 == TIO_copy_attrs (ncid, vi.varid, dontcopy_attr,
                             ncid, test_id))
     return -1;

   return 0;
}

/*}}}*/

static int test_numerical_attribute_io (int target_ncid) /*{{{*/
{
#define NUM_ATT_VALUE 121
   unsigned char ubi, ub = NUM_ATT_VALUE;
   unsigned short usi, us = NUM_ATT_VALUE;
   char bi, b = NUM_ATT_VALUE;
   short si, s = NUM_ATT_VALUE;
   float fi, f = NUM_ATT_VALUE;
   double di, d = NUM_ATT_VALUE;
   unsigned long long u64i, u64 = NUM_ATT_VALUE;

   if (   (0 != TIO_put_att (target_ncid, NC_GLOBAL, "byte_att", NC_BYTE, 1, &b))
       || (0 != TIO_put_att (target_ncid, NC_GLOBAL, "ubyte_att", NC_UBYTE, 1, &ub))
       || (0 != TIO_put_att (target_ncid, NC_GLOBAL, "short_att", NC_SHORT, 1, &s))
       || (0 != TIO_put_att (target_ncid, NC_GLOBAL, "ushort_att", NC_USHORT, 1, &us))
       || (0 != TIO_put_att (target_ncid, NC_GLOBAL, "uint64_att", NC_UINT64, 1, &u64))
       || (0 != TIO_put_att (target_ncid, NC_GLOBAL, "float_att", NC_FLOAT, 1, &f))
       || (0 != TIO_put_att (target_ncid, NC_GLOBAL, "double_att", NC_DOUBLE, 1, &d)))
     {
        fprintf (stderr, "*** error writing numerial attribute\n");
        return -1;
     }

   if ((0 != TIO_get_att (target_ncid, NC_GLOBAL, "byte_att", NC_BYTE, &bi))
       || (0 != TIO_get_att (target_ncid, NC_GLOBAL, "ubyte_att", NC_UBYTE, &ubi))
       || (0 != TIO_get_att (target_ncid, NC_GLOBAL, "short_att", NC_SHORT, &si))
       || (0 != TIO_get_att (target_ncid, NC_GLOBAL, "ushort_att", NC_USHORT, &usi))
       || (0 != TIO_get_att (target_ncid, NC_GLOBAL, "uint64_att", NC_UINT64, &u64i))
       || (0 != TIO_get_att (target_ncid, NC_GLOBAL, "float_att", NC_FLOAT, &fi))
       || (0 != TIO_get_att (target_ncid, NC_GLOBAL, "double_att", NC_DOUBLE, &di))
       || (bi != b) || (ubi != ub) || (si != s) || (usi != us) || (u64i != u64)
       || (fi != f) || (di != d))
     {
        fprintf (stderr, "*** error reading numerial attribute\n");
        return -1;
     }

   return 0;
}

/*}}}*/

static int create_float_waves (float w0, size_t num_waves, float **pwaves) /*{{{*/
{
   float *waves = NULL;
   size_t i;

   if (NULL == (waves = (float *) malloc (num_waves * sizeof(float))))
     {
        fprintf (stderr, "*** %s: malloc failed\n", __func__);
        return -1;
     }

   for (i = 0; i < num_waves; i++)
     {
        waves[i] = w0 + i;
     }

   *pwaves = waves;

   return 0;
}

/*}}}*/

/* Initial condition: group contains only 1D 'nominal_wavelength'
 * Read 'wavelength' 3D and verify that the result is 'nominal_wavelength'.
 */
static int test_wavelength_nominal (int ncid) /*{{{*/
{
   TIO_Var_Info_Type info = {0};
   const char grpname[] = "band_290_490_nm";
   float *nominal_waves = NULL;
   float *waves = NULL;
   float wave0 = 290.0;
   size_t i, step, xtrack, num_chan, num_xtrack, num_step;
   size_t waves_offset;
   int grp, varid, start[3], count[3];
   int status = -1;

   if (NC_NOERR != nc_inq_grp_ncid (ncid, grpname, &grp))
     {
        fprintf (stderr, "%s: Error accessing group %s\n", __func__, grpname);
        return -1;
     }

   if (NC_NOERR == nc_inq_varid (grp, TEMPO_VAR_WAVELENGTH, &varid))
     {
        fprintf (stderr, "%s: Error test assumes that variable %s does not exist\n", __func__, TEMPO_VAR_WAVELENGTH);
        return -1;
     }

   if (NC_NOERR != nc_inq_varid (grp, TEMPO_VAR_WAVELEN_NOMINAL, &varid))
     {
        fprintf (stderr, "%s: Error accessing variable %s\n", __func__, TEMPO_VAR_WAVELEN_NOMINAL);
        return -1;
     }

   if (0 != TIO_inq_var (grp, "radiance", &info))
     return -1;
   num_step = info.dimlens[0];
   num_xtrack = info.dimlens[1];
   num_chan = info.dimlens[2];

   if (0 != create_float_waves (wave0, num_chan, &nominal_waves))
     return -1;

   if (NULL == (waves = (float *)malloc (num_step * num_xtrack * num_chan * sizeof(float))))
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        return -1;
     }

   start[0] = 0;
   count[0] = num_chan;

   if (0 != TIO_put_var_section (grp, TEMPO_VAR_WAVELEN_NOMINAL, start, count, NC_FLOAT, nominal_waves))
     goto return_status;

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;

   count[0] = num_step;
   count[1] = num_xtrack;
   count[2] = num_chan;

   if (0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELENGTH, start, count, NC_FLOAT, waves))
     goto return_status;

   waves_offset = 0;
   for (step = 0; step < num_step; step++)
     {
        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             float *pwaves = waves + waves_offset;
             for (i = 0; i < num_chan; i++)
               {
                  float diff = pwaves[i] - nominal_waves[i];
                  float avg = (pwaves[i] + nominal_waves[i]) * 0.5;
                  if (fabs(diff) > FLT_EPSILON*fabs(avg))
                    {
                       fprintf (stderr, "*** %s: Error: wavelength mismatch: pwaves[%ld]=%g nominal_waves[%ld]=%g\n",
                                __func__, i, pwaves[i], i, nominal_waves[i]);
                       goto return_status;
                    }
               }
             waves_offset += num_chan;
          }
     }

   status = 0;
return_status:
   free(nominal_waves);
   free(waves);
   if (status)
     {
        fprintf (stderr, "%s: FAIL\n", __func__);
     }

   return status;
}

/*}}}*/

static int write_wavecal_params (int grp, float **pwaves, int *pnum_waves) /*{{{*/
{
   TIO_Var_Info_Type info = {0};
   float *waves = NULL;
   float *wavecal_params = NULL;
   float wave0 = 280.0;  /* intentionally different from nominal_waves[0] */
   float coefs[2];
   size_t n, params_dimlen, step, xtrack, num_step, num_xtrack, len;
   size_t wp_offset;
   int start_spectral_channel, num_spectral_channels, num_coefficients;
   int param_dimid, dimids[3], varid, start[3], count[3];
   int status = -1;

   if (0 != TIO_inq_var (grp, "radiance", &info))
     return -1;
   n = info.dimlens[2];

   if (0 != create_float_waves (wave0, n, &waves))
     return -1;

   *pwaves = waves;
   *pnum_waves = n;

   params_dimlen = 2;

   num_coefficients = params_dimlen;
   num_spectral_channels = n;
   start_spectral_channel = 0;

   if (0 != TIO_def_dim (grp, TEMPO_DIM_WAVECAL_PARAM, params_dimlen, &param_dimid))
     goto return_status;

   dimids[0] = info.dimids[0];
   dimids[1] = info.dimids[1];
   dimids[2] = param_dimid;

   if (0 != TIO_def_var (grp, TEMPO_VAR_WAVECAL_PARAM, TIO_FLOAT, 3, dimids, &varid))
     goto return_status;
   if ((0 != TIO_put_att (grp, varid, "num_coefficients", TIO_INT, 1, &num_coefficients))
       || (0 != TIO_put_att (grp, varid, "start_spectral_channel", TIO_INT, 1, &start_spectral_channel))
       || (0 != TIO_put_att (grp, varid, "num_spectral_channels", TIO_INT, 1, &num_spectral_channels)))
     goto return_status;

   /* In the linear case, the Chebyshev expansion coefficients are fairly obvious:
    *    \lambda(i) = c[0] * T0(x) + c[1] * T1(x)
    * where x(i) = (2*i - i0 - i1) / (i1 - i0)
    *   and  i \in [i0, i1]
    *   so that -1 <= x <= 1
    * Also, T0(x) = 1, T1(x) = x
    * Therefore:
    *    \lambda(i) = c[0] + c[1] * x
    * Here i0=0, i1=n-1, so the coefficients are determined by:
    *          lambda(0) = c[0] - c[1]
    *        lambda(n-1) = c[0] + c[1]
    *  =>  c[0] = (lambda(n-1) + lambda(0)) / 2
    *      c[1] = (lambda(n-1) - lambda(0)) / 2
    */
   coefs[0] = (waves[n-1] + waves[0]) * 0.5;
   coefs[1] = (waves[n-1] - waves[0]) * 0.5;

   num_step = info.dimlens[0];
   num_xtrack = info.dimlens[1];

   len = params_dimlen * num_step * num_xtrack;

   if (NULL == (wavecal_params = (float *)malloc (len * sizeof(float))))
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        goto return_status;
     }

   wp_offset = 0;
   for (step = 0; step < num_step; step++)
     {
        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             float *wp = wavecal_params + wp_offset;
             wp[0] = coefs[0];
             wp[1] = coefs[1];
             wp_offset += params_dimlen;
          }
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = num_step;
   count[1] = num_xtrack;
   count[2] = params_dimlen;

   if (0 != TIO_put_var_section (grp, TEMPO_VAR_WAVECAL_PARAM, start, count, NC_FLOAT, wavecal_params))
     goto return_status;

   status = 0;
return_status:
   free(wavecal_params);

   return status;
}

/*}}}*/

/* Initial condition: group contains only 'nominal_wavelength'.
 * Create wavecal_params and verify that reading 'wavelength' yields
 * wavelengths derived from that, and not from nominal_wavelength.
 */
static int test_wavelength_wavecal (int ncid) /*{{{*/
{
   TIO_Var_Info_Type info = {0};
   const char grpname[] = "band_290_490_nm";
   float *waves1 = NULL;
   float *waves = NULL;
   int grp, step, i, xtrack, num_step, num_xtrack, num_chan, start[3], count[3];
   size_t len, waves_offset;
   int status = -1;

   if (NC_NOERR != nc_inq_grp_ncid (ncid, grpname, &grp))
     {
        fprintf (stderr, "%s: Error accessing group %s\n", __func__, grpname);
        return -1;
     }

   if (0 != TIO_inq_var (grp, "radiance", &info))
     return -1;

   num_step = info.dimlens[0];
   num_xtrack = info.dimlens[1];

   if (0 != write_wavecal_params (grp, &waves1, &num_chan))
     goto return_status;

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;

   count[0] = num_step;
   count[1] = num_xtrack;
   count[2] = num_chan;

   len = num_step * num_xtrack * num_chan;
   if (NULL == (waves = (float *)malloc (len * sizeof(float))))
     {
        fprintf (stderr, "*** %s: malloc failed\n", __func__);
        goto return_status;
     }

   if (0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELENGTH, start, count, NC_FLOAT, waves))
     goto return_status;

   waves_offset = 0;
   for (step = 0; step < num_step; step++)
     {
        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             float *pwaves = waves + waves_offset;
             for (i = 0; i < num_chan; i++)
               {
                  float diff = pwaves[i] - waves1[i];
                  float avg = (pwaves[i] + waves1[i]) * 0.5;
                  if (fabs(diff) > FLT_EPSILON*fabs(avg))
                    {
                       fprintf (stderr, "*** %s: Error: wavelength mismatch: pwaves[%d]=%g waves1[%d]=%g\n",
                                __func__, i, pwaves[i], i, waves1[i]);
                       goto return_status;
                    }
               }
             waves_offset += num_chan;
          }
     }

   status = 0;
return_status:
   free(waves1);
   free(waves);
   return status;
}

/*}}}*/

/* Verify reading 'wavelength' succeeds when 'nominal_wavelength'
 * and 'wavecal_params' are not present.
 */
static int test_wavelength_literal (int ncid) /*{{{*/
{
   const char varname[] = "wavelength";
   size_t i, num_waves = 32;
   int dimid, varid, start, count;
   float w0 = 200.0;
   float *waves = NULL;
   float *in_waves = NULL;
   int status = -1;

   if (NULL == (in_waves = (float *)malloc (num_waves * sizeof(float))))
     {
        fprintf (stderr, "*** %s: malloc failed\n", __func__);
        return -1;
     }

   if (0 != create_float_waves (w0, num_waves, &waves))
     goto return_status;

   if ((NC_NOERR != nc_def_dim (ncid, "test_wavedim", num_waves, &dimid))
       || (NC_NOERR != nc_def_var (ncid, varname, NC_FLOAT, 1, &dimid, &varid)))
     {
        fprintf (stderr, "*** %s: Error defining test variable '%s'\n", __func__, varname);
        goto return_status;
     }

   start = 0;
   count = num_waves;
   if (0 != TIO_put_var_section (ncid, varname, &start, &count, NC_FLOAT, waves))
     goto return_status;

   if (0 != TIO_get_var_section (ncid, varname, &start, &count, NC_FLOAT, in_waves))
     goto return_status;

   for (i = 0; i < num_waves; i++)
     {
        if (in_waves[i] != waves[i])
          {
             fprintf (stderr, "*** %s: input/output mismatch: wrote waves[%ld] = %g  read in_waves[%ld] = %g\n",
                      __func__, i, waves[i], i, in_waves[i]);
             goto return_status;
          }
     }

   status = 0;
return_status:
   free(waves);
   free(in_waves);
   return status;
}

/*}}}*/

static int test_wavelength_input (const char *file) /*{{{*/
{
   int ncid;
   int status = -1;

   if (-1 == TIO_open (file, NC_WRITE, &ncid))
     return -1;

   if ((0 != test_wavelength_nominal (ncid))
       || (0 != test_wavelength_wavecal (ncid))
       || (0 != test_wavelength_literal (ncid)))
     goto close_and_return;

   status = 0;
close_and_return:
   if (0 != TIO_close (ncid))
     return -1;

   return status;
}

/*}}}*/

static int check_granule_ident (int ncid) /*{{{*/
{
   int scan_num, granule_num, granule_flag, itest, varid;
   int start=0, count=1;
   double tstart, tend, dtest;

   scan_num = 10;
   granule_num = 3;
   granule_flag = 1;

   if (0 != tio_write_granule_ident_indices (ncid, scan_num, granule_num))
     return -1;
   if (0 != tio_write_granule_flag_var (ncid, granule_flag))
     return -1;

   if ((0 != TIO_get_att (ncid, NC_GLOBAL, "scan_num", NC_INT, &itest)
        || (itest != scan_num)))
     {
        fprintf (stderr, "*** %s: Error expected scan_num=%d, got %d\n",
                 __func__, scan_num, itest);
        return -1;
     }

   if ((0 != TIO_get_att (ncid, NC_GLOBAL, "granule_num", NC_INT, &itest)
        || (itest != granule_num)))
     {
        fprintf (stderr, "*** %s: Error expected granule_num=%d, got %d\n",
                 __func__, granule_num, itest);
        return -1;
     }

   if (0 != tio_inq_varid (ncid, "granule_flag", &varid))
     return -1;

   if ((0 != TIO_get_var_section (ncid, "granule_flag", &start, &count, NC_INT, &itest)
        || (itest != granule_flag)))
     {
        fprintf (stderr, "*** %s: Error expected granule_flag=%d, got %d\n",
                 __func__, granule_flag, itest);
        return -1;
     }

   if ((0 != tio_time_utcstr_to_taix ("2018-09-26T18:00:00Z", &tstart))
       || (0 != tio_time_utcstr_to_taix ("2018-09-26T18:00:01Z", &tend)))
     return -1;

   if (0 != tio_write_granule_ident_times (ncid, tstart, tend))
     return -1;

   if ((0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &dtest)
        || (dtest != tstart)))
     {
        fprintf (stderr, "*** %s: Error expected tstart=%17.15e, got %17.15e\n",
                 __func__, tstart, dtest);
        return -1;
     }

   if ((0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, &dtest)
        || (dtest != tend)))
     {
        fprintf (stderr, "*** %s: Error expected tend=%17.15e, got %17.15e\n",
                 __func__, tend, dtest);
        return -1;
     }

   return 0;
}

/*}}}*/

static int write_lonlat_arrays (int ncid, int ntracks, int nxtrack) /*{{{*/
{
   float *lon=NULL, *lat=NULL, *vza=NULL;
   int *inrqf=NULL;
   float lon_0, lat_0, dlon, dlat;
   size_t len = nxtrack * ntracks;
   int i, j, grp, start[2], count[2], status = -1;

   if ((NULL == (lon = (float *)malloc (len * sizeof(float))))
       ||(NULL == (lat = (float *)malloc (len * sizeof(float))))
       ||(NULL == (vza = (float *)malloc (len * sizeof(float))))
       ||(NULL == (inrqf = (int *)malloc (len * sizeof(int))))
      )
     {
        fprintf (stderr, "*** %s: malloc failed\n", __func__);
        goto return_status;
     }

   memset ((char *)inrqf, 0, len * sizeof(int));

   /* At some point, these arrays may need to be flipped in some way */
   lon_0 = -90.0;   dlon = -0.08;
   lat_0 = +45.0;   dlat = -0.08;

   for (i = 0; i < ntracks; i++)
     {
        float *plon = lon + i * nxtrack;
        float *plat = lat + i * nxtrack;
        float *pvza = vza + i * nxtrack;
        float lon_i = lon_0 + i * dlon;
        for (j = 0; j < nxtrack; j++)
          {
             plon[j] = lon_i;
             plat[j] = lat_0 + j * dlat;
             pvza[j] = 0.0;
          }
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = ntracks;
   count[1] = nxtrack;

   if (0 != TIO_inq_grp (ncid, "band_290_490_nm", &grp))
     goto return_status;

   if ((0 != TIO_put_var_section (grp, TEMPO_VAR_LONGITUDE, start, count, TIO_FLOAT, lon))
       ||(0 != TIO_put_var_section (grp, TEMPO_VAR_LATITUDE, start, count, TIO_FLOAT, lat))
       ||(0 != TIO_put_var_section (grp, TEMPO_VAR_VZ_ANGLE, start, count, TIO_FLOAT, vza))
       ||(0 != TIO_put_var_section (grp, TEMPO_VAR_INRQF, start, count, TIO_INT, inrqf))
      )
     goto return_status;

   status = 0;
return_status:
   free(lon);
   free(lat);
   free(vza);
   free(inrqf);
   return status;
}

/*}}}*/

static int test_scan_type (void)
{
   uint16_t scan_label, scan_num, scan_type;

#define TEST_SCAN_NUM  ((uint16_t)10)
#define TEST_SCAN_TYPE  ((uint16_t)40)

   scan_num = TEST_SCAN_NUM;
   scan_type = TEST_SCAN_TYPE;

   if (0 != tio_make_scan_label (&scan_label, scan_num, scan_type))
     {
        fprintf (stderr, "*** tio_make_scan_label failed\n");
        return -1;
     }

   tio_parse_scan_label (scan_label, &scan_num, &scan_type);
   if ((scan_num != TEST_SCAN_NUM) || (scan_type != TEST_SCAN_TYPE))
     {
        fprintf (stderr, "*** tio_make_scan_type failed\n");
        return -1;
     }

   return 0;
}

static int test_l1_radiance (const char *file, int ntracks, int nxtrack, int ny) /*{{{*/
{
   int ncid, varid, status, grp, err=-1;
   char data_name[] = TEMPO_VAR_RADIANCE;
   char err_name[] = TEMPO_VAR_RADIANCE_ERROR;
   char attr_name[] = "foo";
#define BUFSIZE 1024
   char namebuf[BUFSIZE];
   int target_ncid, read_scan_num, scan_num;
   int field_type = TIO_FLOAT;
   int attr_type_in, attr_type = TIO_INT64, attr_type_conversion = TIO_UINT;
   int attr_len = 1, attr_len_in;
   long long attr = UINT_MAX;
   long long attr_in;
   unsigned int attr_in_conversion;
   char *grp_name;
   char *string_array[1] = {NULL};
   char *string_array_expected[1] = {"a string attribute"};
   float *data = NULL, *data_err = NULL;
   float *data_in = NULL, *data_err_in = NULL;
   double *dbl_err=NULL;
   int data_size = ntracks * nxtrack * ny;
   int start[3], count[3], sub_grp;
   int processing_level;
   int processing_level_type;
   TIO_Scan_Group_Type sgrps[] =
     {
        {"band_290_490_nm", 0, 0},
        {"band_540_740_nm", 0, 0},
     };
   int i, num_sgrps = sizeof (sgrps) / sizeof(sgrps[0]);
   TIO_Scan_Ident_Type *scan_ident = NULL;
   double test_timestamp_out, test_timestamp_in;

   if (0 != test_scan_type ())
     return -1;

   for (i = 0; i < num_sgrps; i++)
     {
        TIO_Scan_Group_Type *s = &sgrps[i];
        s->num_xtrack = nxtrack;
        s->num_channels = ny;
     }

   /* nc_set_log_level(3); */

   start[0] = 0;       start[1] = 0;       start[2] = 0;
   count[0] = ntracks; count[1] = nxtrack; count[2] = ny;

   if ((NULL == (data = generate_data (data_size)))
       || (NULL == (data_err = generate_err (data_size, 1.e-4))))
     {
        fprintf (stderr, "*** error generating data\n");
        return -1;
     }

   if (NULL == (data_in = (float *) malloc (2 * data_size * sizeof(float))))
     {
        fprintf (stderr, "*** malloc failed\n");
        goto cleanup;
     }
   memset ((char *)data_in, 0, 2 * data_size * sizeof(float));
   data_err_in = data_in + data_size;

   if (NULL == (dbl_err = (double *) malloc (data_size * sizeof(double))))
     {
        fprintf (stderr, "*** malloc failed\n");
        goto cleanup;
     }
   for (i = 0; i < data_size; i++)
     {
        dbl_err[i] = (double) data_err[i];
     }

   if (-1 == TIO_create (file, NC_NETCDF4, &ncid))
     goto cleanup;

   if (-1 == test_def_grp (ncid))
        goto cleanup;

   if (-1 == TIO_l1_radiance_template (ncid, ntracks, num_sgrps, sgrps))
     {
        fprintf (stderr, "*** failed creating L1 radiance template in %s\n", file);
        goto cleanup;
     }

   if (0 != check_granule_ident (ncid))
     goto cleanup;

   grp_name = sgrps[0].name;

   if (NC_NOERR != (status = nc_inq_grp_full_ncid (ncid, grp_name, &grp)))
     {
        fprintf (stderr, "*** error finding group %s in file %s (%s)\n",
                 grp_name, file, nc_strerror(status));
        goto cleanup;
     }

   if (-1 == test_dims (ncid))
     goto cleanup;

   if (-1 == TIO_put_var_section (grp, data_name, start, count, field_type, data))
     {
        fprintf (stderr, "*** failed writing %s in file %s\n",
                 data_name, file);
        goto cleanup;
     }
   /* write as a float */
   if (-1 == TIO_put_var_section (grp, err_name, start, count, field_type, data_err))
     {
        fprintf (stderr, "*** failed writing %s in file %s\n",
                 err_name, file);
        goto cleanup;
     }
   /* write again as a double, just to exercise the type conversion code */
   if (-1 == TIO_put_var_section (grp, err_name, start, count, TIO_DOUBLE, dbl_err))
     {
        fprintf (stderr, "*** failed writing %s in file %s\n",
                 err_name, file);
        goto cleanup;
     }

   if (-1 == test_def_var (grp, err_name, field_type))
     goto cleanup;

   /* test writing to attributes */
   if (NC_NOERR != (status = nc_inq_varid (grp, data_name, &varid)))
     {
        fprintf (stderr, "*** error finding variable %s in file %s (%s)\n",
                 data_name, file, nc_strerror(status));
        goto cleanup;
     }

   if (-1 == TIO_put_att (grp, varid, attr_name, attr_type, attr_len, &attr))
     {
        fprintf (stderr, "*** TIO_put_att failed\n");
        goto cleanup;
     }

   /* Having written an attribute, try writing a different type value.
    * The library should generate an error.
    */
     {
        int one=1;
        fprintf (stderr, "expect error here:\n");
        if (0 == TIO_put_att (grp, varid, attr_name, TIO_INT, attr_len, &one))
          {
             fprintf (stderr, "*** expected attribute type mismatch error \n");
             goto cleanup;
          }
     }

   /* test writing to enum attributes */
   processing_level = TIO_PROC_LEVEL_1A;
   if ((-1 == TIO_inq_att (ncid, NC_GLOBAL, "processing_level", &processing_level_type, NULL))
       || (-1 == TIO_put_att (ncid, NC_GLOBAL, "processing_level", processing_level_type, 1, &processing_level)))
     {
        fprintf (stderr, "*** error writing to enum attribute\n");
        goto cleanup;
     }

   /* test writing a timestamp */
   test_timestamp_out = 5.e8;
   if (0 != TIO_write_timestamp (ncid, NC_GLOBAL, "test_timestamp", test_timestamp_out))
     goto cleanup;
   if (0 != TIO_get_att (ncid, NC_GLOBAL, "test_timestamp_since_epoch", NC_DOUBLE, &test_timestamp_in))
     goto cleanup;
   if (test_timestamp_in != test_timestamp_out)
     {
        fprintf (stderr, "*** error reading timestamp\n");
        goto cleanup;
     }

   if (0 != write_lonlat_arrays (ncid, ntracks, nxtrack))
     goto cleanup;

   if ((0 != TIO_def_grp (grp, "subgroup", &sub_grp))
       || (0 != tio_def_l1_radiance_angle_vars (sub_grp))
       || (0 != tio_def_var_ground_pixel_quality_flag  (sub_grp))
       || (0 != tio_def_var_radiance_status_flag (sub_grp)))
     goto cleanup;

   if (0 != tio_sync (ncid))
     goto cleanup;

   if (-1 == TIO_close (ncid))
     goto cleanup;

   if (0 != test_wavelength_input (file))
     goto cleanup;

   if (-1 == TIO_open (file, NC_NOWRITE, &ncid))
     goto cleanup;

   if (0 != tio_use_file_epoch (ncid))
     goto cleanup;

   if (NC_NOERR != (status = nc_inq_grp_full_ncid (ncid, grp_name, &grp)))
     {
        fprintf (stderr, "*** error finding group %s in file %s (%s)\n",
                 grp_name, file, nc_strerror(status));
        goto cleanup;
     }

   /* check attribute properties */
   if (-1 == TIO_inq_att (grp, varid, attr_name, &attr_type_in, &attr_len_in))
     {
        fprintf (stderr, "*** TIO_inq_att failed\n");
        goto cleanup;
     }
   if ((attr_type != attr_type_in) || (attr_len != attr_len_in))
     {
        fprintf (stderr, "*** mismatched attr properties\n");
        fprintf (stderr, "attr_type = %d  attr_type_in=%d\n", attr_type, attr_type_in);
        fprintf (stderr, "attr_len  = %d  attr_len_in =%d\n", attr_len, attr_len_in);
        goto cleanup;
     }
   if (-1 == TIO_get_att (grp, varid, attr_name, attr_type, &attr_in))
     {
        fprintf (stderr, "*** TIO_get_att failed\n");
        goto cleanup;
     }
   if (attr != attr_in)
     {
        fprintf (stderr, "*** read wrong attribute value\n");
        fprintf (stderr, "attr = %lld  attr_in=%lld\n", attr, attr_in);
        goto cleanup;
     }

   /* test attribute type conversion */
   if (-1 == TIO_get_att (grp, varid, attr_name, attr_type_conversion, &attr_in_conversion))
     {
        fprintf (stderr, "*** TIO_get_att failed\n");
        goto cleanup;
     }
   if ((unsigned int) attr != attr_in_conversion)
     {
        fprintf (stderr, "*** conversion read wrong attribute value\n");
        fprintf (stderr, "attr = %lld  attr_in_conversion=%u\n", attr, attr_in_conversion);
        goto cleanup;
     }

   /* test reading enum attributes */
   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "processing_level", processing_level_type, &processing_level))
     {
        fprintf (stderr, "*** error reading enum attribute\n");
        goto cleanup;
     }
   if (processing_level != TIO_PROC_LEVEL_1A)
     {
        fprintf (stderr, "*** error:  processing_level=%u expected %u\n",
                 processing_level, TIO_PROC_LEVEL_1A);
        goto cleanup;
     }

   /* test variable input */
   if (-1 == TIO_get_var_section (grp, data_name, start, count, field_type, data_in))
     {
        fprintf (stderr, "*** error reading variable %s from file %s (%s)\n",
                 data_name, file, nc_strerror(status));
        goto cleanup;
     }
   if (compare_data (data_size, data, data_in))
     goto cleanup;

   if (-1 == TIO_get_var_section (grp, err_name, start, count, field_type, data_err_in))
     {
        fprintf (stderr, "*** error reading variable %s from file %s (%s)\n",
                 err_name, file, nc_strerror(status));
        goto cleanup;
     }
   /* Note that data_err is modified in-place on output
    * (to make compression more efficient), but data_err_in
    * is _not_ modified on input -- so the comparison is
    * expected to work.
    */
   if (compare_data (data_size, data_err, data_err_in))
     goto cleanup;
   /* read as a double just to exercise the conversion code */
   if (-1 == TIO_get_var_section (grp, err_name, start, count, TIO_DOUBLE, dbl_err))
     {
        fprintf (stderr, "*** error reading variable %s from file %s (%s)\n",
                 err_name, file, nc_strerror(status));
        goto cleanup;
     }
   /* exercise the I/O method enable/disable function */
   if ((0 != _TIO_set_io_method_enable (err_name, 0, 0))
       || (0 != _TIO_set_io_method_enable (err_name, 1, 1))
       || (-1 != _TIO_set_io_method_enable ("nonexistent", 1, 1)))
     {
        fprintf (stderr, "*** Error controlling per-variable I/O methods\n");
        goto cleanup;
     }

   /* test granule id functions */
   if (-1 == TIO_filename_from_granule (ncid, "test", 1, 1, namebuf, sizeof(namebuf)))
     {
        fprintf (stderr, "*** Error generating filename from granule id\n");
        goto cleanup;
     }

   if (NULL == (scan_ident = TIO_new_scan_ident ()))
     goto cleanup;

   if (-1 == TIO_attach_granule_ident (ncid, scan_ident))
     goto cleanup;

   if (-1 == TIO_create (namebuf, NC_NETCDF4, &target_ncid))
     {
        fprintf (stderr, "*** Error creating file %s\n", namebuf);
        goto cleanup;
     }

   if (-1 == TIO_copy_granule_ident (ncid, target_ncid))
     {
        fprintf (stderr, "*** Error copying granule id to %s\n", namebuf);
        goto cleanup;
     }
   if (-1 == tio_copy_granule_flag_var (ncid, target_ncid))
     goto cleanup;

   if (-1 == TIO_label_product (target_ncid, "just testing", 1))
     {
        fprintf (stderr, "*** Error labeling granule %s\n", namebuf);
        goto cleanup;
     }

   if ((0 != TIO_get_att (target_ncid, NC_GLOBAL, "scan_num", NC_INT, &read_scan_num))
       || (0 != TIO_get_att (ncid, NC_GLOBAL, "scan_num", NC_INT, &scan_num)))
     {
        fprintf (stderr, "*** Error reading scan_num\n");
        goto cleanup;
     }
   if (read_scan_num != scan_num)
     {
        fprintf (stderr,
                 "*** Error: value mismatch: read_scan_num=%d scan_num=%d\n",
                 read_scan_num, scan_num);
        goto cleanup;
     }

   if (-1 == TIO_write_scan_ident (target_ncid, scan_ident))
     goto cleanup;

   /* Test numerical attribute types */
   if (-1 == test_numerical_attribute_io (target_ncid))
     goto cleanup;

   /* Test string attribute I/O and free */
   if ((-1 == TIO_put_att (target_ncid, NC_GLOBAL, "a_string_att", NC_STRING, 1, string_array_expected))
       || (-1 == TIO_get_att (target_ncid, NC_GLOBAL, "a_string_att", NC_STRING, string_array)))
     {
        fprintf (stderr, "*** error reading time_reference attribute\n");
        goto cleanup;
     }
   if (0 != strcmp (string_array[0], string_array_expected[0]))
     {
        fprintf (stderr, "*** error:  string_array[0]=%s (expected %s)\n",
                 string_array[0], string_array_expected[0]);
        (void) TIO_free_string (1, string_array);
        goto cleanup;
     }
   (void) TIO_free_string (1, string_array);

   if (-1 == TIO_close (target_ncid))
     goto cleanup;

   (void) remove (namebuf);

   if (-1 == TIO_close (ncid))
     goto cleanup;

   err = 0;
cleanup:
   free(data);
   free(data_in);
   free(data_err);
   free(dbl_err);
   TIO_free_scan_ident (scan_ident);

   if (err) fprintf (stderr, "*** TEST FAILED (test_l1_radiance)\n");
   return err;
}

/*}}}*/

static int test_l1_irradiance (const char *file, int ntracks, int nxtrack, int ny) /*{{{*/
{
   int ncid, err=-1;
   TIO_Scan_Group_Type sgrps[] =
     {
        {"band_290_490_nm", 0, 0},
        {"band_540_740_nm", 0, 0},
     };
   int i, num_sgrps = sizeof (sgrps) / sizeof(sgrps[0]);

   for (i = 0; i < num_sgrps; i++)
     {
        TIO_Scan_Group_Type *s = &sgrps[i];
        s->num_xtrack = nxtrack;
        s->num_channels = ny;
     }

   /* nc_set_log_level(3); */

   if (-1 == TIO_create (file, NC_NETCDF4, &ncid))
     goto cleanup;

   if (-1 == TIO_l1_irradiance_template (ncid, ntracks, num_sgrps, sgrps))
     {
        fprintf (stderr, "*** failed creating L1 irradiance template in %s\n", file);
        goto cleanup;
     }

   if (-1 == TIO_close (ncid))
     goto cleanup;

   err = 0;
cleanup:

   if (err) fprintf (stderr, "*** TEST FAILED (test_l1_irradiance)\n");
   return err;
}

/*}}}*/

int main (void)
{
   int ntracks=8, nxtrack=6, ny=5;
   const char *epoch = "2000-01-01T12:00:00Z";
   const char *buf_expected = "TEMPO_" TEMPO_PROD_TYPE_RAD "_L1_V01_20000101T120000Z.nc";
   char buf[72];

   if (0 != tio_time_set_taix_epoch (epoch))
     return 1;

   if ( __tio_filename_string (buf, sizeof(buf), 0.0, TEMPO_PROD_TYPE_RAD, 1, 1) < 0)
     {
        fprintf (stderr, "*** Error: __tio_filename_string failed\n");
        return 1;
     }

   if (0 != strcmp (buf, buf_expected))
     {
        fprintf (stderr, "*** Error: filename mismatch: buf=%s, expected=%s\n",
                 buf, buf_expected);
        return 1;
     }

   if (test_l1_radiance ("delete_radiance.nc", ntracks, nxtrack, ny))
     return 1;

   if (test_l1_irradiance ("delete_irradiance.nc", ntracks, nxtrack, ny))
     return 1;

   return 0;
}
