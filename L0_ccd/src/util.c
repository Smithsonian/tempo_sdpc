#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>

#include <tell.h>
#include <tio.h>
#include "config.h"
#include "util.h"

#define MAX_ENABLE_STATE_SETTINGS  100

typedef struct
{
   /* These const char * pointers are assumed to point to either
    * a string literal, or to memory managed by libconfig
    */
   const char *name;
   const char *enum_state;
   int boolean_state;
   int state_type;
}
Enable_State_Type;

static Enable_State_Type Enable_Settings[MAX_ENABLE_STATE_SETTINGS];
static int Num_Enable_Settings;

static int find_enable_setting (const char *name)
{
   int i;

   for (i = 0; i < Num_Enable_Settings; i++)
     {
        Enable_State_Type *s = &Enable_Settings[i];
        if (0 == strcmp (s->name, name))
          {
             return i;
             break;
          }
     }

   return -1;
}

static int enable_state_new (const char *name, int type, int boolean_state, const char *enum_state)
{
   Enable_State_Type *s;
   int n;

   if (Num_Enable_Settings == MAX_ENABLE_STATE_SETTINGS)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: enable state table is full", __func__);
        return -1;
     }

   if ((n = find_enable_setting (name)) < 0)
     {
        n = Num_Enable_Settings;
        Num_Enable_Settings++;
     }

   s = &Enable_Settings[n];
   s->state_type = type;
   s->boolean_state = boolean_state;
   s->enum_state = enum_state;
   s->name = name;

   return 0;
}

int enable_state_define (config_t *cfg, const char *state_name)
{
   config_setting_t *s, *m;
   const char *enum_state = NULL;
   int boolean_state = 0;
   int type;

   if ((NULL == (s = config_lookup (cfg, "enable")))
       || (NULL == (m = config_setting_get_member (s, state_name))))
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: reading %s from group 'enable' in %s",
                     __func__, state_name, config_error_file(cfg));
        return -1;
     }

   type = config_setting_type (m);

   switch (type)
     {
      case CONFIG_TYPE_BOOL:
        boolean_state = config_setting_get_bool (m);
        tell_vlog (TELL_MSGTYPE_INFO, 1, "SELECT %s: %s", state_name, boolean_state ? "ON" : "OFF");
        break;

      case CONFIG_TYPE_STRING:
        enum_state = config_setting_get_string (m);
        tell_vlog (TELL_MSGTYPE_INFO, 1, "SELECT %s: %s", state_name, enum_state ? enum_state : "NULL");
        break;

      default:
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: invalid setting: %s",  __func__, state_name);
        return -1;
     }

   return enable_state_new (state_name, type, boolean_state, enum_state);
}

const char *enable_state_query_enum (const char *name)
{
   Enable_State_Type *s;
   int n;

   if ((n = find_enable_setting (name)) < 0)
     return NULL;

   s = &Enable_Settings[n];

   if (s->state_type != CONFIG_TYPE_STRING)
     return NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: %s", name, s->enum_state);

   return s->enum_state;
}

int enable_state_query_bool (const char *name)
{
   Enable_State_Type *s;
   int n;

   if ((n = find_enable_setting (name)) < 0)
     return -1;

   s = &Enable_Settings[n];

   if (s->state_type != CONFIG_TYPE_BOOL)
     return -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: %s", name, s->boolean_state ? "ON" : "OFF");

   return s->boolean_state;
}

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

int bsearch_f (float t, const float *x, int n)
{
   int n0, n1, n2;
   float xt;

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

char *expand_string (const char *s)
{
   wordexp_t we;
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof(wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding string: %s", __func__, s);
        goto return_error;
     }

   if (NULL == (s_exp = strdup (we.we_wordv[0])))
     {
        tell_verror (TELL_MALLOC_ERROR,
                     "%s: strdup failed", __func__);
     }

return_error:
   wordfree (&we);
   return s_exp;
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

   if (path == NULL)
     return 0;

   if (NULL != (path_basename = strrchr (path, '/')))
     {
        path_basename++;
     }
   else path_basename = path;

   return tio_meta_append_string (meta, "input_files", path_basename);
}
