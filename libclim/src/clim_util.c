#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>

#include <tell.h>

extern int make_pressure_filename (int month, char *buf, int bufsize);
extern int make_climatology_filename (const char *species, int month, int hour,
                                      char *path, int path_bufsize);

static char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
     }

   return s_exp;
}

static int expand_buffer_in_place (char *buf, size_t bufsize)
{
   char *fmt = NULL;
   size_t n;

   if (NULL == (fmt = expand_string (buf)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: expand_string failed", __func__);
        return -1;
     }

   n = strlen(fmt) + 1;

   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %ld > %ld",
                     __func__, n, bufsize);
        free (fmt);
        return -1;
     }

   strncpy (buf, fmt, bufsize);
   free(fmt);
   return 0;
}

int make_climatology_filename (const char *species, int month, int hour,
                               char *buf, int bufsize)
{
   const char fmt[] = "$SDPC_REFDATA_DIR/climatology/%s/gcnr_TRC%s_2013%02d_%02d00.nc";
   int n;

   if ((n = snprintf (buf, bufsize, fmt, species, species, month, hour)) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        return -1;
     }
   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %d > %d",
                     __func__, n, bufsize);
        return -1;
     }

   return expand_buffer_in_place (buf, bufsize);
}

int make_pressure_filename (int month, char *buf, int bufsize)
{
   const char fmt[] = "$SDPC_REFDATA_DIR/climatology/PRES/gcnr_pressure_2013%02d.nc";
   int n;

   if ((n = snprintf (buf, bufsize, fmt, month)) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        return -1;
     }
   if (n >= bufsize)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: buffer size exceeded: string length = %d > %d",
                     __func__, n, bufsize);
        return -1;
     }

   return expand_buffer_in_place (buf, bufsize);
}
