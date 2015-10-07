#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "netcdf.h"
#include "tio.h"
#include "tio_template.h"

static float *generate_data (int n)
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

static float *generate_relerr (int n, float xp_min, float xp_max,
                               float u_missing)
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
        float xp = xp_min + i * (xp_max - xp_min) / (n-1.0);
        y[i] = pow (10.0, xp);
     }
   y[n/2] = u_missing;

   return y;
}

static int compare_data (int n, float *out, float *in)
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

static int test_def_grp (int ncid)
{
   int ignore_grp;

   if (-1 == TIO_def_grp (ncid, "xxx", &ignore_grp))
     {
        fprintf (stderr, "*** TIO_def_grp failed\n");
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

static int test_l1_radiance (const char *file, int ntracks, int nxtrack, int ny)
{
   int ncid, varid, status, grp, err=-1;
   char data_name[] = TEMPO_VAR_RADIANCE;
   char relerr_name[] = TEMPO_VAR_RADIANCE_ERROR;
   char attr_name[] = "foo";
   int field_type = TIO_FLOAT;
   int attr_type_in, attr_type = TIO_INT64, attr_type_conversion = TIO_UINT;
   size_t attr_len = 1, attr_len_in;
   long long attr = UINT_MAX;
   long long attr_in;
   unsigned int attr_in_conversion;
   char *grp_name;
   float *data = NULL, *relerr = NULL;
   float *data_in = NULL, *relerr_in = NULL;
   double *dbl_relerr=NULL;
   int data_size = ntracks * nxtrack * ny;
   int start[3], count[3];
   int processing_level;
   int processing_level_type;
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

   start[0] = 0;       start[1] = 0;       start[2] = 0;
   count[0] = ntracks; count[1] = nxtrack; count[2] = ny;

   if ((NULL == (data = generate_data (data_size)))
       || (NULL == (relerr = generate_relerr (data_size, -4.0, 2.0, -TIO_FILL_FLOAT))))
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
   relerr_in = data_in + data_size;

   if (NULL == (dbl_relerr = (double *) malloc (data_size * sizeof(double))))
     {
        fprintf (stderr, "*** malloc failed\n");
        goto cleanup;
     }
   for (i = 0; i < data_size; i++)
     {
        dbl_relerr[i] = (double) relerr[i];
     }

   if (NC_NOERR != (status = nc_create (file, NC_NETCDF4, &ncid)))
     {
        fprintf (stderr, "*** error opening %s (%s)\n",
                 file, nc_strerror(status));
        goto cleanup;
     }

   if (-1 == test_def_grp (ncid))
        goto cleanup;

   if (-1 == TIO_l1_radiance_template (ncid, ntracks, num_sgrps, sgrps))
     {
        fprintf (stderr, "*** failed creating L1 radiance template in %s\n", file);
        goto cleanup;
     }

   grp_name = sgrps[0].name;

   if (NC_NOERR != (status = nc_inq_grp_full_ncid (ncid, grp_name, &grp)))
     {
        fprintf (stderr, "*** error finding group %s in file %s (%s)\n",
                 grp_name, file, nc_strerror(status));
        goto cleanup;
     }

   if (-1 == TIO_put_var_section (grp, data_name, start, count, field_type, data))
     {
        fprintf (stderr, "*** failed writing %s in file %s\n",
                 data_name, file);
        goto cleanup;
     }
   /* write as a float */
   if (-1 == TIO_put_var_section (grp, relerr_name, start, count, field_type, relerr))
     {
        fprintf (stderr, "*** failed writing %s in file %s\n",
                 relerr_name, file);
        goto cleanup;
     }
   /* write again as a double, just to exercise the type conversion code */
   if (-1 == TIO_put_var_section (grp, relerr_name, start, count, TIO_DOUBLE, dbl_relerr))
     {
        fprintf (stderr, "*** failed writing %s in file %s\n",
                 relerr_name, file);
        goto cleanup;
     }

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

   if (NC_NOERR != (status = nc_close (ncid)))
     {
        fprintf (stderr, "*** error closing file %s\n", file);
        goto cleanup;
     }

   if (NC_NOERR != (status = nc_open (file, NC_NOWRITE, &ncid)))
     {
        fprintf (stderr, "*** error opening file %s\n", file);
        goto cleanup;
     }

   if (NC_NOERR != (status = nc_inq_grp_full_ncid (ncid, grp_name, &grp)))
     {
        fprintf (stderr, "*** error finding group %s in file %s (%s)\n",
                 grp_name, file, nc_strerror(status));
        goto cleanup;
     }

   /* check attribute propreties */
   if (-1 == TIO_inq_att (grp, varid, attr_name, &attr_type_in, &attr_len_in))
     {
        fprintf (stderr, "*** TIO_inq_att failed\n");
        goto cleanup;
     }
   if ((attr_type != attr_type_in) || (attr_len != attr_len_in))
     {
        fprintf (stderr, "*** mismatched attr properties\n");
        fprintf (stderr, "attr_type = %d  attr_type_in=%d\n", attr_type, attr_type_in);
        fprintf (stderr, "attr_len = %lu  attr_len_in=%lu\n", attr_len, attr_len_in);
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

   if (-1 == TIO_get_var_section (grp, relerr_name, start, count, field_type, relerr_in))
     {
        fprintf (stderr, "*** error reading variable %s from file %s (%s)\n",
                 relerr_name, file, nc_strerror(status));
        goto cleanup;
     }
   if (compare_data (data_size, relerr, relerr_in))
     goto cleanup;
   /* read as a double just to exercise the conversion code */
   if (-1 == TIO_get_var_section (grp, relerr_name, start, count, TIO_DOUBLE, dbl_relerr))
     {
        fprintf (stderr, "*** error reading variable %s from file %s (%s)\n",
                 relerr_name, file, nc_strerror(status));
        goto cleanup;
     }

   if (NC_NOERR != (status = nc_close (ncid)))
     {
        fprintf (stderr, "*** error closing file %s\n", file);
        goto cleanup;
     }

   err = 0;
cleanup:
   free(data);
   free(data_in);
   free(relerr);
   free(dbl_relerr);

   if (err) fprintf (stderr, "*** TEST FAILED (test_l1_radiance)\n");
   return err;
}

static int test_l1_irradiance (const char *file, int ntracks, int nxtrack, int ny)
{
   int ncid, status, err=-1;
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

   if (NC_NOERR != (status = nc_create (file, NC_NETCDF4, &ncid)))
     {
        fprintf (stderr, "*** error opening %s (%s)\n",
                 file, nc_strerror(status));
        goto cleanup;
     }

   if (-1 == TIO_l1_irradiance_template (ncid, ntracks, num_sgrps, sgrps))
     {
        fprintf (stderr, "*** failed creating L1 irradiance template in %s\n", file);
        goto cleanup;
     }

   if (NC_NOERR != (status = nc_close (ncid)))
     {
        fprintf (stderr, "*** error closing file %s\n", file);
        goto cleanup;
     }

   err = 0;
cleanup:

   if (err) fprintf (stderr, "*** TEST FAILED (test_l1_irradiance)\n");
   return err;
}

int main (void)
{
   int ntracks=8, nxtrack=6, ny=5;

   if (test_l1_radiance ("delete_radiance.nc", ntracks, nxtrack, ny))
     return 1;

   if (test_l1_irradiance ("delete_irradiance.nc", ntracks, nxtrack, ny))
     return 1;

   return 0;
}
