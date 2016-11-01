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

static float *generate_err (int n, float err)
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

static int test_dims (int ncid)
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

static int dontcopy_attr (const char *attr)
{
   return (0 == strcmp (attr, "_FillValue"));
}

static int test_def_var (int ncid, const char *name, int type)
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

static int test_l1_radiance (const char *file, int ntracks, int nxtrack, int ny)
{
   int ncid, varid, status, grp, err=-1;
   char data_name[] = TEMPO_VAR_RADIANCE;
   char err_name[] = TEMPO_VAR_RADIANCE_ERROR;
   char attr_name[] = "foo";
#define BUFSIZE 1024
   char namebuf[BUFSIZE];
   int target_ncid, read_scan_seq_num, scan_seq_num;
   int field_type = TIO_FLOAT;
   int attr_type_in, attr_type = TIO_INT64, attr_type_conversion = TIO_UINT;
   int attr_len = 1, attr_len_in;
   long long attr = UINT_MAX;
   long long attr_in;
   unsigned int attr_in_conversion;
   char *grp_name;
   float *data = NULL, *data_err = NULL;
   float *data_in = NULL, *data_err_in = NULL;
   double *dbl_err=NULL;
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
   TIO_Scan_Ident_Type *scan_ident = NULL;

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

   if (-1 == TIO_close (ncid))
     goto cleanup;

   if (-1 == TIO_open (file, NC_NOWRITE, &ncid))
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
   if (-1 == TIO_filename_from_granule (ncid, "test", 1, namebuf, sizeof(namebuf)))
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

   if (-1 == TIO_label_product (target_ncid, "just testing", 1))
     {
        fprintf (stderr, "*** Error labeling granule %s\n", namebuf);
        goto cleanup;
     }

   if ((0 != TIO_get_att (target_ncid, NC_GLOBAL, "scan_seq_num", NC_INT, &read_scan_seq_num))
       || (0 != TIO_get_att (ncid, NC_GLOBAL, "scan_seq_num", NC_INT, &scan_seq_num)))
     {
        fprintf (stderr, "*** Error reading scan_seq_num\n");
        goto cleanup;
     }
   if (read_scan_seq_num != scan_seq_num)
     {
        fprintf (stderr,
                 "*** Error: value mismatch: read_scan_seq_num=%d scan_seq_num=%d\n",
                 read_scan_seq_num, scan_seq_num);
        goto cleanup;
     }

   if (-1 == TIO_write_scan_ident (target_ncid, scan_ident))
     goto cleanup;

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

static int test_l1_irradiance (const char *file, int ntracks, int nxtrack, int ny)
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

int main (void)
{
   int ntracks=8, nxtrack=6, ny=5;

   if (test_l1_radiance ("delete_radiance.nc", ntracks, nxtrack, ny))
     return 1;

   if (test_l1_irradiance ("delete_irradiance.nc", ntracks, nxtrack, ny))
     return 1;

   return 0;
}
