#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"
#include "tio.h"

static int try_history (const char *filename, const char *str, const char *label)
{
   int ncid, status = -1;

   if (0 != TIO_create (filename, NC_NETCDF4, &ncid))
     goto return_status;

   /* test creating a global history attribute */
   if (0 != tio_append_history (ncid, str))
     goto return_status;
   /* test appending to an existing history attribute */
   if (0 != tio_append_history (ncid, str))
     goto return_status;

   if (0 != TIO_close (ncid))
     goto return_status;

   status = 0;
return_status:
   if (status)
     {
        fprintf (stderr, "*** ERROR: %s failed\n", label);
     }
   return status;
}

int main (int argc, char **argv)
{
   const char *file = "history_test.nc";
   char *str = NULL;
   size_t len;
   int status = EXIT_FAILURE;

   if (NULL == (str = tio_concat_argv (argc, argv, NULL, 0)))
     {
        fprintf (stderr, "*** ERROR: tio_concat_argv\n");
        goto exit_status;
     }

   len = strlen(str) + 1;

   if (0 != try_history (file, str, "tio_concat_argv alloc"))
     goto exit_status;

   free (str);
   str = NULL;

   if (NULL == (str = malloc (len * sizeof(char))))
     {
        fprintf (stderr, "*** ERROR: malloc failed (len = %ld)\n", len);
        goto exit_status;
     }
   if (NULL == tio_concat_argv (argc, argv, str, len))
     goto exit_status;

   if (0 != try_history (file, str, "preallocated"))
     goto exit_status;

   status = EXIT_SUCCESS;
exit_status:
   free(str);
   return status;
}
