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

int main (void)
{
   int ncid, status, grp;
   int ntracks=10, nxtrack=10, ny=10;
   char file[] = "delete_me.nc";
   char grp_name[] = TIO_GRP_NAME_BAND1;
   char field_name[] = TIO_VAR_NAME_RADIANCE;
   char attr_name[] = "foo";
   int field_type = TIO_FLOAT;
   int attr_type_in, attr_type = TIO_INT64, attr_type_conversion = TIO_UINT;
   size_t attr_len = 1, attr_len_in;
   long long attr = UINT_MAX;
   long long attr_in;
   unsigned int attr_in_conversion;
   float *data = NULL;
   float *data_in = NULL;
   int data_size = ntracks * nxtrack * ny;
   int err = 1;
   int track, num_write;

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

   if (-1 == TIO_create_l1b_template (ncid, ntracks, nxtrack, ny))
     {
        fprintf (stderr, "*** failed creating L1b template in %s\n", file);
        goto cleanup;
     }

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

   if (-1 == TIO_put_att (grp, field_name, attr_name, attr_type, attr_len, &attr))
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
        if (0 == TIO_put_att (grp, field_name, attr_name, TIO_INT, attr_len, &one))
          {
             fprintf (stderr, "*** expected attribute type mismatch error \n");
             goto cleanup;
          }
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
   if (-1 == TIO_inq_att (grp, field_name, attr_name, &attr_type_in, &attr_len_in))
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
   if (-1 == TIO_get_att (grp, field_name, attr_name, attr_type, &attr_in))
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
   if (-1 == TIO_get_att (grp, field_name, attr_name, attr_type_conversion, &attr_in_conversion))
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

   if (err) fprintf (stderr, "*** TEST FAILED\n");
   return err;
}
