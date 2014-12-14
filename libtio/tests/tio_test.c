#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "netcdf.h"
#include "tio.h"

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

static int compare_data (int n, float *out, float *in)
{
   int i;

   for (i = 0; i < n; i++)
     {
        if (in[i] != out[i])
          {
             fprintf (stderr, "*** error:  data[%d] %15.10e != %15.10e\n",
                      i, in[i], out[i]);
             return -1;
          }
     }

   return 0;
}

static int test_l1_radiance (const char *file, int ntracks, int nxtrack, int ny)
{
   int ncid, varid, status, grp, err=-1;
   char field_name[] = TEMPO_VAR_RADIANCE;
   char attr_name[] = "foo";
   int field_type = TIO_FLOAT;
   int attr_type_in, attr_type = TIO_INT64, attr_type_conversion = TIO_UINT;
   size_t attr_len = 1, attr_len_in;
   long long attr = UINT_MAX;
   long long attr_in;
   unsigned int attr_in_conversion;
   char *grp_name;
   float *data = NULL;
   float *data_in = NULL;
   int data_size = ntracks * nxtrack * ny;
   int track, num_write;
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

   track = 0;
   num_write = ntracks;

   if (NULL == (data = generate_data (2 * data_size)))
     {
        fprintf (stderr, "*** error generating data\n");
        return -1;
     }
   data_in = data + data_size;

   if (NC_NOERR != (status = nc_create (file, NC_NETCDF4, &ncid)))
     {
        fprintf (stderr, "*** error opening %s (%s)\n",
                 file, nc_strerror(status));
        goto cleanup;
     }

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

   if (-1 == TIO_put_var_section (grp, field_name, track, num_write, field_type, data))
     {
        fprintf (stderr, "*** failed writing field %s in file %s\n",
                 field_name, file);
        goto cleanup;
     }

   /* test writing to attributes */
   if (NC_NOERR != (status = nc_inq_varid (grp, field_name, &varid)))
     {
        fprintf (stderr, "*** error finding variable %s in file %s (%s)\n",
                 field_name, file, nc_strerror(status));
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
   if (-1 == TIO_get_var_section (grp, field_name, track, num_write, field_type, data_in))
     {
        fprintf (stderr, "*** error reading variable %s from file %s (%s)\n",
                 field_name, file, nc_strerror(status));
        goto cleanup;
     }
   if (compare_data (data_size, data, data_in))
     goto cleanup;

   if (NC_NOERR != (status = nc_close (ncid)))
     {
        fprintf (stderr, "*** error closing file %s\n", file);
        goto cleanup;
     }

   err = 0;
cleanup:
   free(data);

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
   int ntracks=10, nxtrack=10, ny=10;

   if (test_l1_radiance ("delete_radiance.nc", ntracks, nxtrack, ny))
     return 1;

   if (test_l1_irradiance ("delete_irradiance.nc", ntracks, nxtrack, ny))
     return 1;

   return 0;
}
