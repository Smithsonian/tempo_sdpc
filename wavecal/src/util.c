#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>

#include <tell.h>
#include <tio.h>
#include "config.h"
#include "util.h"

int bsearch_d (double t, const double *x, int n)
{
   int n0, n1, n2;
   double xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

double *alloc_doubles (int n)
{
   double *a = NULL;

   if (NULL == (a = (double *)MALLOC (n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   return a;
}

int find_x (double x, const double *a, int na)
{
   if (x < a[0] || x >= a[na-1])
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: value=%g outside array bounds [%g,%g)",
                     __func__, x, a[0], a[na-1]);
        return -1;
     }
   return bsearch_d (x, a, na);
}

/* return 0 means absent, -1 means error, 1 means read successfully */
int read_config_float_array (config_setting_t *s, const char *name,
                             double **pa, size_t *pnum_a)
{
   config_setting_t *fs;
   double *a;
   unsigned int i, na;
   int len;

   *pa = NULL;
   if (pnum_a) *pnum_a = 0;

   /* absent/empty is ok */
   if ((NULL == (fs = config_setting_get_member (s, name)))
       || (0 == (len = config_setting_length (fs))))
     return 0;

   if (NULL == (a = (double *)MALLOC (len * sizeof (double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   na = len;

   for (i = 0; i < na; i++)
     {
        a[i] = config_setting_get_float_elem (fs, i);
     }

   *pa = a;
   if (pnum_a) *pnum_a = na;

   return 1;
}

char *expand_path (const char *path)
{
   wordexp_t we;
   char *path_exp = NULL;

   memset ((char *)&we, 0, sizeof(wordexp_t));

   if ((0 != wordexp (path, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, path);
        goto return_error;
     }

   if (NULL == (path_exp = strdup (we.we_wordv[0])))
     {
        tell_verror (TELL_MALLOC_ERROR,
                     "%s: strdup failed", __func__);
     }

return_error:
   wordfree (&we);
   return path_exp;
}

char *path_concat (const char *dir, const char *basename)
{
   int status, len;
   char *s;

   len = strlen (dir) + strlen(basename) + 2;
   if (NULL == (s = MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   status = snprintf (s, len, "%s/%s", dir, basename);
   if ((status < 0) || (status >= len))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        FREE(s);
        return NULL;
     }

   return s;
}

int meta_record_basename (TIO_Meta_Type *meta, const char *path)
{
   const char *path_basename;

   if (NULL != (path_basename = strrchr (path, '/')))
     {
        path_basename++;
     }
   else path_basename = path;

   return tio_meta_append_string (meta, "input_pointer", path_basename);
}
