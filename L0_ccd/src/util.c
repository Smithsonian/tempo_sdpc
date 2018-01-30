#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <tell.h>
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
