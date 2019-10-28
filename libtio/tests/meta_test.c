#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <tio.h>
#include <_tio.h>

static int try_date_time_range_scan (TIO_Meta_Type *meta)
{
   _pTIO_Granule_Ident_Type g2 =
     {
        .next = NULL,
        .scan_num = 1,
        .granule_num = 2,
        .tstart = 2.0,
        .tend = 3.0,
        .tstart_str = "2019-05-09T03:00:00Z",
        .tend_str = "2019-05-09T04:00:00Z"
     };
   _pTIO_Granule_Ident_Type g1 =
     {
        .next = &g2,
        .scan_num = 1,
        .granule_num = 1,
        .tstart = 1.0,
        .tend = 2.0,
        .tstart_str = "2019-05-09T01:00:00Z",
        .tend_str = "2019-05-09T02:00:00Z"
     };
   TIO_Scan_Ident_Type lst = {&g1};

   return tio_meta_set_datetime_range_scan (meta, &lst);
};

int main (void)
{
   TIO_Meta_Type *meta = NULL;
   int status = EXIT_FAILURE;
   int meta_int = 12345;
   unsigned meta_uint = UINT_MAX;
   double meta_dbl_array[] = {2.0e42, -4.0e-42, 6.0, 8.0};
   char *meta_str_array[] = {"The", "Sun", "Also", "Rises"};
   char *strings_to_append[] = {"append A", "append B", "append C"};
   const char *ncfile = "delete_radiance.nc";
   int ncid, grp, uv_grp, i;

   if (NULL == (meta = tio_meta_open ()))
     goto cleanup_and_exit;

   if (0 != tio_meta_set (meta, "INT_KEYWORD1", TIO_META_TYPE_INT, 1, &meta_int))
     goto cleanup_and_exit;

   if (0 != tio_meta_set (meta, "UINT_KEYWORD1", TIO_META_TYPE_UINT, 1, &meta_uint))
     goto cleanup_and_exit;

   if (0 != tio_meta_set (meta, "DBL_KEYWORD1", TIO_META_TYPE_DOUBLE,
                          sizeof(meta_dbl_array)/sizeof(double), meta_dbl_array))
     goto cleanup_and_exit;

   if (0 != tio_meta_set (meta, "STR_KEYWORD1", TIO_META_TYPE_STRING,
                          sizeof(meta_str_array)/sizeof(char *), meta_str_array))
     goto cleanup_and_exit;

   if (0 != tio_meta_set (meta, "STR_KEYWORD2", TIO_META_TYPE_STRING, 1, "Just one string"))
     goto cleanup_and_exit;

   if (0 != tio_meta_set (meta, "STR_APPEND", TIO_META_TYPE_STRING, 1, "The first string"))
     goto cleanup_and_exit;
   if ((0 != tio_meta_append_string (meta, "STR_APPEND", "append 1"))
       || (0 != tio_meta_append_string (meta, "STR_APPEND", "append 1"))  /* duplicates should be ignored */
       || (0 != tio_meta_append_string (meta, "STR_APPEND", "append 2")))
     {
        goto cleanup_and_exit;
     }
   for (i = 0; i < 3; i++)
     {
        if (0 != tio_meta_append_string (meta, "STR_APPEND", strings_to_append[i]))
          goto cleanup_and_exit;
     }

   if (0 != tio_meta_set (meta, "STR_KEYWORD3", TIO_META_TYPE_STRING, 1, "Another string"))
     goto cleanup_and_exit;

   if (0 != try_date_time_range_scan (meta))
     goto cleanup_and_exit;

   if (0 != tio_meta_set_standard (meta, "product_file_name.nc", "product_short_name", 1, "v1.2.3"))
     goto cleanup_and_exit;

   if (0 != tio_meta_expand_stream (meta, stdin, stdout))
     goto cleanup_and_exit;

   if (0 != tio_meta_expand_file (meta, "meta.in", "xxx"))
     goto cleanup_and_exit;

   if ((0 != TIO_open (ncfile, NC_WRITE, &ncid))
       || (0 != TIO_def_grp (ncid, "metadata", &grp)))
     {
        fprintf (stderr, "*** Error opening netcdf file metadata group: %s\n", ncfile);
        goto cleanup_and_exit;
     }

   if (0 != tio_meta_set_datetime_range (meta, ncid))
     goto cleanup_and_exit;

   if (0 != TIO_inq_grp (ncid, "band_290_490_nm", &uv_grp))
     goto cleanup_and_exit;

   if (0 != tio_meta_set_lev1_bounding_polygon_and_centroid (meta, uv_grp))
     goto cleanup_and_exit;

   if (0 != tio_meta_write_ncattr (meta, grp))
     goto cleanup_and_exit;

   if (0 != TIO_close (ncid))
     {
        fprintf (stderr, "*** Error closing netcdf file: %s\n", ncfile);
        goto cleanup_and_exit;
     }

   /* Test appending to a netcdf file string metadata keyword.
    * For this purpose, we don't care about the ascii .met file.
    */

   if ((0 != TIO_open (ncfile, NC_WRITE, &ncid))
       || (0 != TIO_inq_grp (ncid, "metadata", &grp)))
     {
        fprintf (stderr, "*** Error opening netcdf file metadata group: %s\n", ncfile);
        goto cleanup_and_exit;
     }

   if (0 != tio_meta_ncinit (meta, grp, "STR_KEYWORD2", TIO_META_TYPE_STRING))
     goto cleanup_and_exit;

   if (0 != tio_meta_append_string (meta, "STR_KEYWORD2", "and another one"))
     goto cleanup_and_exit;

   /* run this here just to improve code coverage */
   (void) tio_meta_set_noexpand (meta, "STR_KEYWORD2", 1);

   if (0 != tio_meta_write_ncattr (meta, grp))
     goto cleanup_and_exit;
   if (0 != TIO_close (ncid))
     {
        fprintf (stderr, "*** Error closing netcdf file: %s\n", ncfile);
        goto cleanup_and_exit;
     }

   status = EXIT_SUCCESS;
cleanup_and_exit:
   tio_meta_close (meta);
   return status;
}
