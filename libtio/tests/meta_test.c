#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <tio.h>

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
   int ncid, grp;

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
   if ((0 != tio_meta_append (meta, "STR_APPEND", TIO_META_TYPE_STRING, 1, "append 1"))
       || (0 != tio_meta_append (meta, "STR_APPEND", TIO_META_TYPE_STRING, 1, "append 2"))
       || (0 != tio_meta_append (meta, "STR_APPEND", TIO_META_TYPE_STRING, 3, strings_to_append))
      )
     goto cleanup_and_exit;
   
   if (0 != tio_meta_set (meta, "STR_KEYWORD3", TIO_META_TYPE_STRING, 1, "Another string"))
     goto cleanup_and_exit;

   if (0 != tio_meta_expand_stream (meta, stdin, stdout))
     goto cleanup_and_exit;

   if ((0 != TIO_open (ncfile, NC_WRITE, &ncid))
       || (0 != TIO_def_grp (ncid, "metadata", &grp)))
     {
        fprintf (stderr, "*** Error opening netcdf file metadata group: %s\n", ncfile);
        goto cleanup_and_exit;
     }

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
