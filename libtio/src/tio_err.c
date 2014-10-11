#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>

#include "netcdf.h"
#include "tio.h"
#include "_tio.h"

void _pTIO_err_verror (const char *fmt, ...)
{
   va_list ap;

   (void) fprintf (stderr, "**ERROR: ");
   va_start (ap, fmt);
   (void) vfprintf (stderr, fmt, ap);
   va_end (ap);
   (void) fprintf (stderr, "\n");
}

void _pTIO_err_verror_nc (int status, const char *fmt, ...)
{
   va_list ap;

   (void) fprintf (stderr, "**ERROR: ");
   va_start (ap, fmt);
   (void) vfprintf (stderr, fmt, ap);
   va_end (ap);
   (void) fprintf (stderr, "\n(%s)\n", nc_strerror(status));
}
