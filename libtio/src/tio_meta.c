/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

#define __USE_XOPEN   /* for strptime */
#include <time.h>

#include <tell.h>

#include "tio.h"
#include "_tio.h"
#include "tio_template.h"
#include "simplify.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2  (M_PI/2.0)
#endif

#define DEGTORAD  (M_PI/180.0)

#define GEOSPATIAL_BOUNDS_CRS "EPSG:4326"

struct TIO_Meta_Type
{
   TIO_Meta_Type *next;
   char *name;
   void *values;
   int num_values;
   int num_alloc;
   int value_type;
   int noexpand;
   int odl_only;
};

#define KEYNAME_VALID_CHARS "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_."

static void free_string_array (char **a, int na)
{
   int i;

   if (a == NULL)
     return;

   for (i = 0; i < na; i++)
     {
        TIO_FREE(a[i]);
     }

   TIO_FREE(a);
}

static char **new_string_array (int num_strings)
{
   char **a = NULL;

   if (NULL == (a = (char **)TIO_MALLOC (num_strings * sizeof (char *))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)a, 0, num_strings * sizeof (char *));
   return a;
}

static char **dup_string_array (char **strings, int num_strings)
{
   char **a = NULL;
   int i;

   if (NULL == (a = new_string_array (num_strings)))
     return NULL;

   for (i = 0; i < num_strings; i++)
     {
        if (NULL == (a[i] = strdup (strings[i])))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_string_array (a, num_strings);
             return NULL;
          }
     }

   return a;
}

static void meta_free1 (TIO_Meta_Type *meta)
{
   if (meta == NULL)
     return;

   if (meta->value_type == TIO_META_TYPE_STRING)
     {
        free_string_array (meta->values, meta->num_values);
     }
   else TIO_FREE(meta->values);

   TIO_FREE(meta->name);
   TIO_FREE(meta);
}

static void meta_free (TIO_Meta_Type *meta)
{
   if (meta == NULL)
     return;

   while (meta != NULL)
     {
        TIO_Meta_Type *next = meta->next;
        meta_free1(meta);
        meta = next;
     }
}

static TIO_Meta_Type *meta_alloc (void)
{
   TIO_Meta_Type *meta = NULL;

   if (NULL == (meta = (TIO_Meta_Type *)TIO_MALLOC (sizeof *meta)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)meta, 0, sizeof (*meta));

   return meta;
}

static int invalid_metadata_entry (const char *name, int value_type, int num_values, const void *values)
{
   if ((name == NULL) || (*name == 0))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid metadata name: %s",
                     __func__, name ? name : "NULL");
        return 1;
     }

   if (strspn (name, KEYNAME_VALID_CHARS) != strlen(name))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid metadata name: %s",
                     __func__, name);
        return 1;
     }

   if ((num_values < 1) || (values == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid metadata value: %s",
                     __func__, name);
        return 1;
     }

   switch (value_type)
     {
      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported metadata type: %s %d",
                     __func__, name, value_type);
        return 1;

      case TIO_META_TYPE_DOUBLE:
      case TIO_META_TYPE_FLOAT:
      case TIO_META_TYPE_INT:
      case TIO_META_TYPE_UINT:
      case TIO_META_TYPE_CHAR:
      case TIO_META_TYPE_STRING:
        break;
     }

   return 0;
}

static int meta_write_value_ascii (const TIO_Meta_Type *m, int want_num_values, FILE *fp)
{
   int i;
   char *cp = (char *)m->values;
   int *ip = (int *)m->values;
   unsigned int *uip = (unsigned int *)m->values;
   double *dp = (double *)m->values;
   float *sp = (float *)m->values;
   char **spp = (char **)m->values;
   int need_parens;

   if (want_num_values)
     {
        (void) fprintf (fp, "%d", m->num_values);
        return 0;
     }

   need_parens = ((m->num_values > 1)
                  && (m->value_type != TIO_META_TYPE_CHAR));

   if (need_parens)
     {
        fputs ("(", fp);
     }

   switch (m->value_type)
     {
      case TIO_META_TYPE_CHAR:
        /* This is assumed to be a text string of NC_CHAR type for fortran compatibility */
        fprintf (fp, "\"%s\"", cp);
        break;

      case TIO_META_TYPE_INT:
        fprintf (fp, "%d", ip[0]);
        for (i = 1; i < m->num_values; i++)
          {
             fprintf (fp, ", %d", ip[i]);
          }
        break;

      case TIO_META_TYPE_UINT:
        fprintf (fp, "%u", uip[0]);
        for (i = 1; i < m->num_values; i++)
          {
             fprintf (fp, ", %u", uip[i]);
          }
        break;

      case TIO_META_TYPE_DOUBLE:
        fprintf (fp, "%0.15g", dp[0]);
        for (i = 1; i < m->num_values; i++)
          {
             fprintf (fp, ", %0.15g", dp[i]);
          }
        break;

      case TIO_META_TYPE_FLOAT:
        fprintf (fp, "%0.7g", sp[0]);
        for (i = 1; i < m->num_values; i++)
          {
             fprintf (fp, ", %0.7g", sp[i]);
          }
        break;

      case TIO_META_TYPE_STRING:
        fprintf (fp, "\"%s\"", spp[0]);
        for (i = 1; i < m->num_values; i++)
          {
             fprintf (fp, ", \"%s\"", spp[i]);
          }
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid value type: %d",
                     __func__, m->value_type);
        return -1;
     }

   if (need_parens)
     {
        fputs (")", fp);
     }

   return 0;
}

static TIO_Meta_Type *meta_new (const char *name, int value_type, int num_values, const void *values)
{
   TIO_Meta_Type *meta;
   int value_num_bytes;

   if (NULL == (meta = meta_alloc()))
     return NULL;

   if (invalid_metadata_entry (name, value_type, num_values, values))
     return NULL;

   meta->value_type = value_type;
   meta->num_values = num_values;

   switch (value_type)
     {
      case TIO_META_TYPE_DOUBLE:
        value_num_bytes = 8;
        break;

      case TIO_META_TYPE_FLOAT:
      case TIO_META_TYPE_INT:
      case TIO_META_TYPE_UINT:
        value_num_bytes = 4;
        break;

      case TIO_META_TYPE_CHAR:
      case TIO_META_TYPE_STRING:
        value_num_bytes = 1;
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported value type %d",
                     __func__, value_type);
        goto return_error;
     }

   if (NULL == (meta->name = strdup (name)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
        goto return_error;
     }

    if (value_type == TIO_META_TYPE_STRING)
     {
        char *pv = (char *)values;
        char **ppv = (num_values == 1) ? &pv : (char **)values;

        if (NULL == (meta->values = (void *) dup_string_array (ppv, num_values)))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto return_error;
          }
     }
   else
     {
        size_t value_size = num_values * value_num_bytes;

        if (NULL == (meta->values = TIO_MALLOC (value_size)))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto return_error;
          }
        memcpy ((char *)meta->values, (char *)values, value_size);
     }

   meta->num_alloc = meta->num_values;

   return meta;

return_error:
   meta_free1(meta);
   return NULL;
}

void tio_meta_close (TIO_Meta_Type *meta)
{
   if (meta == NULL)
     return;
   meta_free (meta);
}

TIO_Meta_Type *tio_meta_open (void)
{
   return meta_alloc();
}

int tio_meta_set (TIO_Meta_Type *lst, const char *name,
                  int value_type, int num_values, const void *values)
{
   TIO_Meta_Type *meta = NULL;

   if (lst == NULL)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: metadata list is uninitialized", __func__);
        return -1;
     }

   if (NULL == (meta = meta_new (name, value_type, num_values, values)))
     return -1;

   /* If we already have an entry with this name, silently replace it.
    * Otherwise, append the new entry.
    */

   for ( ; lst->next != NULL; lst = lst->next)
     {
        if (0 == strcmp (lst->next->name, name))
          {
             meta->next = lst->next->next;
             meta_free1 (lst->next);
             break;
          }
     }

   lst->next = meta;

   return 0;
}

static TIO_Meta_Type *find_key_by_name (TIO_Meta_Type *lst, const char *name)
{
   TIO_Meta_Type *meta;

   for (meta = lst; meta != NULL; meta = meta->next)
     {
        if (meta->name)
          {
             if (0 == strcmp (meta->name, name))
               return meta;
          }
     }

   return NULL;
}

int tio_meta_append_string (TIO_Meta_Type *lst, const char *name, const char *str)
{
   TIO_Meta_Type *meta;
   char **ppc = NULL;
   char *s = NULL;
   int i;

   /* If the keyword isn't found, create it */
   if (NULL == (meta = find_key_by_name (lst, name)))
     {
        return tio_meta_set (lst, name, TIO_META_TYPE_STRING, 1, str);
     }

   if (meta->value_type == TIO_META_TYPE_CHAR)
     {
        /* 2 strings, plus a space delimiter, plus null char */
        int len = 2 + strlen(str) + strlen ((char *)meta->values);

        /* silently ignore duplicates */
        if (NULL != strstr ((char *)meta->values, str))
          return 0;

        if (NULL == (s = (char *) TIO_MALLOC (len * sizeof(char))))
          {
             tell_verror(TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return -1;
          }
        /* Strings appended to char array will be space delimited.
         * It's not totally general, but this is intended for a particular
         * application -- char arrays for metadata will normally contain text,
         * so inserting a space delimiter is a reasonable default. */
        snprintf (s, len, "%s %s", (char *)meta->values, str);
        TIO_FREE(meta->values);
        meta->values = s;
        meta->num_values = len;
        meta->num_alloc = len;
        return 0;
     }

   if (meta->value_type != TIO_META_TYPE_STRING)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: is not a string keyword (value_type=%d)",
                     __func__, meta->value_type);
        return -1;
     }

   /* silently ignore duplicates */
   ppc = (char **)meta->values;
   for (i = 0; i < meta->num_values; i++)
     {
        if (0 == strcmp (ppc[i], str))
          return 0;
     }

   if (meta->num_alloc < meta->num_values + 1)
     {
        int new_size = (meta->num_values + 1) + 10;
        void *v;
        if (NULL == (v = TIO_REALLOC (meta->values, new_size * sizeof(char *))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
             return -1;
          }
        meta->values = v;
        meta->num_alloc = new_size;
     }

   if (NULL == (s = strdup (str)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   ppc = (char **)meta->values;
   ppc[meta->num_values++] = s;

   return 0;
}

int tio_meta_set_noexpand (TIO_Meta_Type *lst, const char *name, int noexpand)
{
   TIO_Meta_Type *meta;

   if (NULL == (meta = find_key_by_name (lst, name)))
     return -1;

   meta->noexpand = noexpand;

   return 0;
}

static int append_att_chars (TIO_Meta_Type *meta, int grp, int varid,
                             const char *name, int att_len)
{
   char *att = NULL;
   int status = -1;

   if (NULL == (att = (char *)TIO_MALLOC (att_len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   if (NC_NOERR != nc_get_att_text (grp, varid, name, att))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading attribute %s", __func__, name);
        goto return_status;
     }

   if (0 != tio_meta_append_string (meta, name, att))
     goto return_status;

   status = 0;
return_status:
   TIO_FREE(att);
   return status;
}

static int append_att_strings (TIO_Meta_Type *meta, int grp, int varid,
                               const char *name, int att_len)
{
   char **atts = NULL;
   int i, status = -1;

   if (NULL == (atts = (char **)TIO_MALLOC (att_len * sizeof(char *))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   if (NC_NOERR != nc_get_att_string (grp, varid, name, atts))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading attribute %s", __func__, name);
        goto return_status;
     }

   for (i = 0; i < att_len; i++)
     {
        if (0 != tio_meta_append_string (meta, name, atts[i]))
          goto return_status;
     }

   status = 0;
return_status:
   /* free space allocated on string input */
   if (atts) (void) nc_free_string (att_len, atts);
   TIO_FREE(atts);
   return status;
}

int tio_meta_ncinit (TIO_Meta_Type *meta, int grp, const char *name, int value_type)
{
   int is_char, varid = NC_GLOBAL;
   nc_type att_type;
   size_t att_len;

   switch (value_type)
     {
      case TIO_META_TYPE_STRING:
        is_char = 0;
        break;

      case TIO_META_TYPE_CHAR:
        is_char = 1;
        break;

      default:
        tell_verror (TELL_UNSUPPORTED_ERROR, "%s: invalid metadata attribute type", __func__);
        return -1;
     }

   /* When the attribute doesn't exist, return silently */
   if (NC_NOERR != nc_inq_att (grp, varid, name, &att_type, &att_len))
     return 0;

   if ((is_char == 1) && (att_type == NC_CHAR))
     {
        return append_att_chars (meta, grp, varid, name, att_len);
     }
   else if ((is_char == 0) && (att_type == NC_STRING))
     {
        return append_att_strings (meta, grp, varid, name, att_len);
     }

   tell_verror (TELL_RUNTIME_ERROR,
                "%s: attribute type mismatch: %s (att_type=%d, value_type=%d)",
                __func__, name, att_type, value_type);
   return -1;
}

static int key_search (const char *s, size_t slen, char **pos, size_t *len)
{
   char *beg_key;
   const char *p;
   const char *end;
   const char *pk;
   size_t n;

   *pos = NULL;
   *len = 0;

   p = s;
   end = s + slen;

   while (p < end)
     {
        /* smallest keyword requires 4 characters:  ${X} */
        if ((NULL == (beg_key = strchr (p, '$')))
            || ((beg_key + 3) >= end))
          return 0;

        pk = beg_key + 1;

        if (*pk != '{')
          {
             p = pk;
             continue;
          }

        pk++;

        /* # is a keyword modifier, standing for num_values */
        if (*pk == '#') pk++;

        n = strspn (pk, KEYNAME_VALID_CHARS);

        if (pk + n >= end)
          return 0;

        pk += n;

        if (*pk != '}')
          {
             p = pk;
             continue;
          }

        /* Length includes keyword delimiters, ${} */
        *pos = beg_key;
        *len = (pk + 1) - beg_key;
        break;
     }

   return 0;
}

static int key_expand (const TIO_Meta_Type *meta, const char *keypos, size_t keylen, FILE *fp)
{
   char buf[128];
   const TIO_Meta_Type *k;
   const char *name;
   size_t namelen, len;
   int want_num_values;

   /* Keyword of the form ${KEYNAME} expands to keyword value. */
   name = &keypos[2];
   namelen = keylen - 3;

   /* Keyword of the form ${#KEYNAME} expands to the number of values. */
   if (name[0] == '#')
     {
        want_num_values = 1;
        name = &keypos[3];
        namelen--;
     }
   else want_num_values = 0;

   for (k = meta; k != NULL; k = k->next)
     {
        if (k->name == NULL)
          continue;

        len = strlen(k->name);

        if (namelen != len)
          continue;

        if (0 == strncmp (name, k->name, len))
          {
             if (k->noexpand == 0)
               return meta_write_value_ascii (k, want_num_values, fp);
             else break;
          }
     }

   /* Pass-through undefined or noexpand keywords */
   memset ((char *)buf, 0, sizeof(buf));
   len = keylen < sizeof(buf) ? keylen : sizeof(buf);
   strncpy (buf, keypos, len);
   tell_vwarn (0, "%s: undefined keyword: %s", __func__, buf);
   fprintf (fp, "%s", buf);

   return 0;
}

/* Expand keywords of the form ${[#]KEYWORDNAME} */
int tio_meta_expand_stream (const TIO_Meta_Type *meta, FILE *fp_template, FILE *fp_outfile)
{
   char *line = NULL;
   int return_status = -1;

   while (1)
     {
        char *beg;
        char *end;
        size_t linelen;
        int status;

        TIO_FREE(line);
        line = NULL;

        if ((status = tio_fgets (&line, &linelen, fp_template)) < 0)
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading metadata pattern file", __func__);
             return -1;
          }
        if (status == 0)
          break;

        beg = line;
        end = line + linelen;

        while (beg < end)
          {
             char *keypos;
             size_t keylen;

             if (-1 == key_search (beg, end-beg, &keypos, &keylen))
               goto return_error;

             if (keypos == NULL)
               {
                  /* No keywords, so just write the line */
                  (void) fwrite (beg, sizeof(char), end-beg, fp_outfile);
                  beg = end;
               }
             else
               {
                  /* Write out any text preceeding the keyword,
                   * expand the keyword, then continue scanning
                   */
                  if (keypos-beg > 0)
                    {
                       (void) fwrite (beg, sizeof(char), keypos-beg, fp_outfile);
                    }

                  if (0 != key_expand (meta, keypos, keylen, fp_outfile))
                    goto return_error;

                  beg = keypos + keylen;
               }
          }
     }

   return_status = 0;

return_error:
   TIO_FREE(line);

   return return_status;
}

static char *mkstrcat (const char *str, const char *suffix)
{
   char *s = NULL;
   int len;

   len = strlen(str) + strlen(suffix) + 1;
   if (NULL == (s = TIO_MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   strcpy (s, str);
   strcat (s, suffix);

   return s;
}

int tio_meta_expand_file (const TIO_Meta_Type *meta, const char *ptemplate_file,
                          const char *outfile_root)
{
   FILE *fp_template = NULL;
   FILE *fp_outfile = NULL;
   char *tmpl = NULL;
   char *outfile = NULL;
   int do_rename=0, status = -1;

   /* The output file root must always be specified */
   if (outfile_root == NULL)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: outfile_root == NULL", __func__);
        return -1;
     }

   if (NULL == (outfile = mkstrcat (outfile_root, ".met")))
     goto return_status;

   /* ptemplate_file==NULL means that we're updating an existing template
    * that's in the same directory as the input file, and that has the
    * same filename, apart from a .met extension.
    * In this case, we'll use the outfile name as the template,
    * write the expanded file to a temporary file, and then rename.
    */
   if (ptemplate_file)
     {
        tmpl = (char *)ptemplate_file;
     }
   else
     {
        tmpl = outfile;

        if (0 != access (tmpl, F_OK | R_OK))
          {
             tell_vwarn (0, "metadata template not found: %s", tmpl);
             TIO_FREE(outfile);
             return 0;
          }

        if (NULL == (outfile = mkstrcat (tmpl, ".tmp")))
          goto return_status;

        do_rename = 1;
     }

   if (NULL == (fp_template = fopen (tmpl, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading",
                     __func__, tmpl);
        goto return_status;
     }

   if (NULL == (fp_outfile = fopen (outfile, "w")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for writing", __func__, outfile);
        goto return_status;
     }

   if (0 != tio_meta_expand_stream (meta, fp_template, fp_outfile))
     goto return_status;

   if ((do_rename != 0)
       && (0 != rename (outfile, tmpl)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: rename failed: %s -> %s",
                     __func__, outfile, tmpl);
     }

   status = 0;
return_status:
   if (fp_template) (void) fclose (fp_template);
   if (fp_outfile)
     {
        if (0 != fclose (fp_outfile))
          {
             tell_verror (TELL_IO_ERROR, "%s: closing file: %s", __func__, outfile);
          }
     }

   TIO_FREE(outfile);
   if (ptemplate_file == NULL)
     {
        TIO_FREE(tmpl);
     }

   return status;
}

int tio_meta_write_ncattr (const TIO_Meta_Type *meta, int grp)
{
   if (meta == NULL)
     return -1;

   for (; meta != NULL; meta = meta->next)
     {
        int xtype;

        if (meta->odl_only) continue;

        switch (meta->value_type)
          {
           case TIO_META_TYPE_UNDEFINED:
             continue;
           case TIO_META_TYPE_INT: xtype = NC_INT;
             break;
           case TIO_META_TYPE_UINT: xtype = NC_UINT;
             break;
           case TIO_META_TYPE_FLOAT: xtype = NC_FLOAT;
             break;
           case TIO_META_TYPE_DOUBLE: xtype = NC_DOUBLE;
             break;
           case TIO_META_TYPE_STRING: xtype = NC_STRING;
             break;
           case TIO_META_TYPE_CHAR: xtype = NC_CHAR;
             break;
           default:
             tell_verror (TELL_RUNTIME_ERROR, "%s: invalid metadata type: %d",
                          __func__, meta->value_type);
             return -1;
          }

        if (0 != TIO_put_att (grp, NC_GLOBAL, meta->name, xtype, meta->num_values, meta->values))
          return -1;
     }

   return 0;
}

static inline float merge_coordinates (float corner1, float corner2, float center, float fill_value)
{
   if (corner1 == fill_value)
     {
        return (corner2 == fill_value) ? center : corner2;
     }

   return (corner2 == fill_value) ? corner1 : 0.5 * (corner1 + corner2);
}

static int flag_nonmonotonic_columns (const float *lon2d, const int *inrqf, int num_steps, int num_xtrack, int x,
                                      int increasing_eastward, float *delta_lon, int *exclude_column)
{
   float lon_a, lon_b;
   int s, sign;

   memset ((char *)delta_lon, 0, num_steps * sizeof(float));

   for (s = 0; s < num_steps; s++)
     {
        if ((inrqf + s * num_xtrack)[x] == 0)
          break;
     }
   if (s == num_steps)
     return 0;

   lon_a = (lon2d + s * num_xtrack)[x];
   s++;
   lon_b = lon_a;

   for (  ; s < num_steps; s++)
     {
        if ((inrqf + s * num_xtrack)[x] == 0)
          {
             lon_b = (lon2d + s * num_xtrack)[x];
             delta_lon[s] = lon_b - lon_a;
             lon_a = lon_b;
          }
     }

   sign = increasing_eastward ? -1 : +1;

   for (s = 0; s < num_steps; s++)
     {
        if (delta_lon[s] * sign < 0) exclude_column[s] = 1;
     }

   return 0;
}

static int intersect_lines (const double *seg1, const double *seg2,
                            double *x, double *y, double rtol)
{
   double x21, x31, x43, y21, y31, y43;
   double num_a, num_b, denom;
   double tol_a, tol_b, tol_d;
   int type;

   /* seg1 = [x1,y1, x2,y2]
    * seg2 = [x3,y3, x4,y4]
    */

   x21 = seg1[2] - seg1[0];
   x31 = seg2[0] - seg1[0];
   x43 = seg2[2] - seg2[0];

   y21 = seg1[3] - seg1[1];
   y31 = seg2[1] - seg1[1];
   y43 = seg2[3] - seg2[1];

   num_a = x43 * y31 - x31 * y43;
   num_b = x21 * y31 - x31 * y21;
   denom = x43 * y21 - x21 * y43;

   tol_a = fabs(x43*y31) + fabs(x31*y43);
   tol_b = fabs(x21*y31) + fabs(x31*y21);
   tol_d = fabs(x43*y21) + fabs(x21*y43);

   if (isnan(rtol)) rtol = DBL_EPSILON;
   tol_a = rtol * ((tol_a > 0) ? tol_a : 1.0);
   tol_b = rtol * ((tol_b > 0) ? tol_b : 1.0);
   tol_d = rtol * ((tol_d > 0) ? tol_d : 1.0);

   /*
    * type:
    * -1 = lines parallel, no intersection point
    *  0 = lines coincide
    *  1 = intersection point on segment 1, but not segment 2
    *  2 = intersection point on segment 2, but not segment 1
    *  3 = intersection point on both segment 1 and segment 2
    */

   if ((fabs(num_a) < tol_a)
       && (fabs(num_b) < tol_b)
       && (fabs(denom) < tol_d))
     {
        /* lines coincide */
        *x = seg1[0];
        *y = seg1[1];
        type = 0;
     }
   else if (fabs(denom) < tol_d)
     {
        /* lines parallel */
        *x = nan("");
        *y = nan("");
        type = -1;
     }
   else
     {
        /* lines intersect */
        double a = num_a / denom;
        double b = num_b / denom;
        *x = seg1[0] + a * x21;
        *y = seg1[1] + a * y21;
        type = 0;
        /* FIXME - floating point equality comparison... */
        if (0 <= a && a <= 1) type += 1;
        if (0 <= b && b <= 1) type += 2;
     }

   return type;
}

/* Brute force test for polygon self-intersection.
 * Shamos-Hoey algorithm is more general and more efficient, but also more complicated. */
static int polygon_self_intersects (const double *px, const double *py, int n)
{
   int a, b;

   for (a = 0; a < n-1; a++)
     {
        double aseg[4];

        aseg[0] = px[a  ];  aseg[1] = py[a  ];
        aseg[2] = px[a+1];  aseg[3] = py[a+1];

        for (b = 0; b < n-1; b++)
          {
             double bseg[4], x, y;

             /* compare sides that do not share a vertex */
             if ((a == b) || (a == (b+1)) || ((a+1) == b)
                 || ((a == 0) && (b == (n-2)))
                 || ((b == 0) && (a == (n-2))))
               continue;

             bseg[0] = px[b  ];  bseg[1] = py[b  ];
             bseg[2] = px[b+1];  bseg[3] = py[b+1];

             if (3 == intersect_lines (aseg, bseg, &x, &y, DBL_EPSILON))
               return 1;
          }
     }

   return 0;
}

/* output bounding polygon resolution */
#ifndef DOUGLAS_PEUCKER_BAND_WIDTH_KM
# define DOUGLAS_PEUCKER_BAND_WIDTH_KM 5.0
#endif
static float Douglas_Peucker_Band_Width_Km = DOUGLAS_PEUCKER_BAND_WIDTH_KM;

/* avoid pixels near the Earth's limb */
#ifndef BOUNDING_POLYGON_MAX_VZA_DEG
# define BOUNDING_POLYGON_MAX_VZA_DEG 80.0
#endif
static float Bounding_Polygon_Max_VZA_Deg = BOUNDING_POLYGON_MAX_VZA_DEG;

int __tio_set_bounding_polygon_controls (float band_km, float vza_max_deg)
{
   Douglas_Peucker_Band_Width_Km = (band_km > 0.0) ? band_km : DOUGLAS_PEUCKER_BAND_WIDTH_KM;
   Bounding_Polygon_Max_VZA_Deg = (vza_max_deg > 0.0) ? vza_max_deg : BOUNDING_POLYGON_MAX_VZA_DEG;
   return 0;
}

int __tio_make_lev1_bounding_polygon (int grp, int *num, float **plon, float **plat)
{
   TIO_Var_Info_Type info;
   double *tmp_lon = NULL, *tmp_lat = NULL;
   float *vza2d=NULL, *lon2d=NULL, *lat2d=NULL, *lon=NULL, *lat=NULL;
   float *lon2d_bnds=NULL, *lat2d_bnds=NULL, *delta_lon=NULL;
   int *inrqf=NULL, *indices=NULL, *exclude_column=NULL;
   int *bx1=NULL, *bx2=NULL, *bs1=NULL, *bs2=NULL, *bdry=NULL, *side=NULL;
   int start[3], count[3];
   int num_steps, num_xtrack, num_pixels, max_num_boundary, num_fixed;
   int s, x, i, n, x_first_ok, x_last_ok, s_first_ok, s_last_ok;
   int varid, no_fill, lon_bounds_status, num_kept, polygon_type;
   int status = -1, increasing_eastward=1;
   float fill_value = TIO_FILL_FLOAT;
   float band_km = Douglas_Peucker_Band_Width_Km;
   float vza_max_deg = Bounding_Polygon_Max_VZA_Deg;
   float lon_i, lat_i;
   int dx=1, ds=1;           /* set >1 to reduce final polygon point density */

   *num = 0;
   *plon = NULL;
   *plat = NULL;

   if (0 != TIO_inq_var (grp, TEMPO_VAR_LONGITUDE, &info))
     return -1;

   num_steps = info.dimlens[0];
   num_xtrack = info.dimlens[1];
   num_pixels = num_xtrack * num_steps;

   max_num_boundary = 2 * (num_steps + num_xtrack);

   if ((NULL == (inrqf = (int *)TIO_MALLOC (num_pixels * sizeof(int))))
       || (NULL == (lon2d = (float *)TIO_MALLOC (num_pixels * sizeof(float))))
       || (NULL == (lat2d = (float *)TIO_MALLOC (num_pixels * sizeof(float))))
       || (NULL == (vza2d = (float *)TIO_MALLOC (num_pixels * sizeof(float))))
       || (NULL == (side = (int *)TIO_MALLOC (max_num_boundary * sizeof(int))))
       || (NULL == (indices = (int *)TIO_MALLOC (2 * max_num_boundary * sizeof(int))))
       || (NULL == (exclude_column = (int *)TIO_MALLOC (num_steps * sizeof(int))))
       || (NULL == (delta_lon = (float *)TIO_MALLOC (num_steps * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_steps;
   count[1] = num_xtrack;

   if ((0 != TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE, start, count, TIO_FLOAT, lon2d))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_LATITUDE, start, count, TIO_FLOAT, lat2d))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_VZ_ANGLE, start, count, TIO_FLOAT, vza2d))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_INRQF, start, count, TIO_INT, inrqf)))
     {
        goto return_status;
     }

   /* Ensure that the INR QF is non-zero for all pixels that
    * have a fill value for a pixel center coordinate.
    */
   num_fixed = 0;
   for (i = 0; i < num_pixels; i++)
     {
        if (inrqf[i] != 0)
          continue;
        if ((lon2d[i] == fill_value) || (lat2d[i] == fill_value))
          {
             inrqf[i] = 1;
             num_fixed++;
          }
     }
   if (num_fixed > 0)
     {
        tell_vwarn (0, "boundary polygon construction found %d pixels with INR quality flag=0 (nominal), and a fill value in the geospatial coordinates", num_fixed);
     }

   /* If we have lon/lat bounds, we'll use them to define the polygon boundaries. */

   tell_push_queue ();
   lon_bounds_status = tio_inq_varid (grp, TEMPO_VAR_LONGITUDE_BOUNDS, &varid);
   tell_pop_queue (1);

   if (lon_bounds_status == 0)
     {
        if ((NULL == (lon2d_bnds = (float *)TIO_MALLOC (4 * num_pixels * sizeof(float))))
            || (NULL == (lat2d_bnds = (float *)TIO_MALLOC (4 * num_pixels * sizeof(float)))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto return_status;
          }

        if (0 != TIO_inq_var_fill (grp, varid, &no_fill, &fill_value))
          goto return_status;

        start[2] = 0;
        count[2] = 4;

        if ((0 != TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS, start, count, TIO_FLOAT, lon2d_bnds))
            || (0 != TIO_get_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS, start, count, TIO_FLOAT, lat2d_bnds)))
          goto return_status;
     }

   /* Make sure the longitude grid is monotonic in the step dimension.
    * Exclude any column with a non-monotonic longitude change
    */
   memset ((char *)exclude_column, 0, num_steps * sizeof(int));
   flag_nonmonotonic_columns (lon2d, inrqf, num_steps, num_xtrack, 0.2*num_xtrack, increasing_eastward, delta_lon, exclude_column);
   flag_nonmonotonic_columns (lon2d, inrqf, num_steps, num_xtrack, 0.8*num_xtrack, increasing_eastward, delta_lon, exclude_column);
   for (s = 0; s < num_steps; s++)
     {
        int *inrqf_step = inrqf + s * num_xtrack;
        for (x = 0; x < num_xtrack; x++)
          {
             if (exclude_column[s]) inrqf_step[x] = 1;
          }
     }

   /* Points with valid (lon,lat) coordinates have inrqf == 0.
    * We want the boundary of this region.
    * Find the boundary points by scanning each row and each column
    * to find the outermost endpoints of the region with inrqf==0.
    */

   bs1 = indices;
   bs2 = bs1 + num_xtrack;
   bx1 = bs2 + num_xtrack;
   bx2 = bx1 + num_steps;
   bdry = bx2 + num_steps;

   /* for each step, find the xtrack dimension boundaries: */
   for (s = 0; s < num_steps; s++)
     {
        float *vza_step = vza2d + s * num_xtrack;
        int *inrqf_step = inrqf + s * num_xtrack;
        x_first_ok = -1;
        x_last_ok = -1;
        for (x = 0; x < num_xtrack; x++)
          {
             if ((inrqf_step[x] == 0) && (vza_step[x] < vza_max_deg))
               {
                  x_last_ok = x;
                  if (x_first_ok < 0)
                    x_first_ok = x;
               }
          }
        if (x_first_ok >= 0)
          {
             bx1[s] = x_first_ok + s * num_xtrack;
             bx2[s] = x_last_ok + s * num_xtrack;
          }
        else
          {
             /* entire column was bad */
             bx1[s] = -1;
             bx2[s] = -1;
          }
     }

   /* for each xtrack, find the step dimension boundaries: */
   for (x = 0; x < num_xtrack; x++)
     {
        s_first_ok = -1;
        s_last_ok = -1;
        for (s = 0; s < num_steps; s++)
          {
             if ((inrqf[x + s * num_xtrack] == 0)
                 && (vza2d[x + s * num_xtrack] < vza_max_deg))
               {
                  s_last_ok = s;
                  if (s_first_ok < 0)
                    s_first_ok = s;
               }
          }
        if (s_first_ok >= 0)
          {
             bs1[x] = x + s_first_ok * num_xtrack;
             bs2[x] = x + s_last_ok * num_xtrack;
          }
        else
          {
             /* entire row was bad */
             bs1[x] = -1;
             bs2[x] = -1;
          }
     }

   /* The region boundary polygon we want is now this set of points
    * {bs2, bx2, bs1, bx1}, with some segments reversed to maintain
    * a consistent CCW ordering, and skipping -1 boundary indices.
    */

   n = 0;
   for (x = 0; x < num_xtrack; x += dx)
     {
        if (bs2[x] >= 0)
          {
             side[n] = 0;
             bdry[n] = bs2[x];
             n++;
          }
     }
   for (s = num_steps-1; s >= 0; s -= ds)
     {
        if (bx2[s] >= 0)
          {
             side[n] = 1;
             bdry[n] = bx2[s];
             n++;
          }
     }
   for (x = num_xtrack-1; x >= 0; x -= dx)
     {
        if (bs1[x] >= 0)
          {
             side[n] = 2;
             bdry[n] = bs1[x];
             n++;
          }
     }
   /* The northernmost boundary points are sometimes jagged and/or
    * confusing to deal with.  The simplest solution is to just omit the
    * the most troublesome of these points, letting any remaining boundary gap
    * close with a single line segment connecting the northernmost endpoints
    * of the eastern and western boundaries.
    */
   for (s = 0; s < num_steps; s += ds)
     {
        if ((bx1[s] % num_xtrack) == 0)
          {
             side[n] = 3;
             bdry[n] = bx1[s];
             n++;
          }
     }

   /* If necessary, we will explicitly close the polygon by appending
    * a copy of the first point */

   if ((NULL == (lon = (float *)TIO_MALLOC(n * sizeof(float))))
       ||(NULL == (lat = (float *)TIO_MALLOC(n * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   if (lon_bounds_status != 0)
     {
        /* pixel bounds not available - use the pixel centers */
        for (i = 0; i < n; i++)
          {
             int k = bdry[i];
             lon[i] = lon2d[k];
             lat[i] = lat2d[k];
          }
     }
   else
     {
        for (i = 0; i < n; i++)
          {
             int k = bdry[i];
             float lon_k = lon2d[k];
             float lat_k = lat2d[k];
             float *lonbnds = lon2d_bnds + k*4;
             float *latbnds = lat2d_bnds + k*4;
             switch (side[i])
               {
                case 0:
                  /* "west"  e.g. NW,SW corners, regardless of actual lon,lat coordinate values  */
                  lon_i = merge_coordinates (lonbnds[1], lonbnds[2], lon_k, fill_value);
                  lat_i = merge_coordinates (latbnds[1], latbnds[2], lat_k, fill_value);
                  break;
                case 1:
                  /* "south", e.g. SW,SE corners */
                  lon_i = merge_coordinates (lonbnds[2], lonbnds[3], lon_k, fill_value);
                  lat_i = merge_coordinates (latbnds[2], latbnds[3], lat_k, fill_value);
                  break;
                case 2:
                  /* "east", e.g. SE,NE corners */
                  lon_i = merge_coordinates (lonbnds[3], lonbnds[0], lon_k, fill_value);
                  lat_i = merge_coordinates (latbnds[3], latbnds[0], lat_k, fill_value);
                  break;
                case 3:
                  /* "north" e.g. NE,NW corners */
                  lon_i = merge_coordinates (lonbnds[0], lonbnds[1], lon_k, fill_value);
                  lat_i = merge_coordinates (latbnds[0], latbnds[1], lat_k, fill_value);
                  break;
                default:
                  tell_verror (TELL_RUNTIME_ERROR, "%s: this should never happen!!", __func__);
                  goto return_status;
               }
             lon[i] = lon_i;
             lat[i] = lat_i;
          }
     }

   /* simplify needs double arrays for temporary storage */
   if (NULL == (tmp_lon = (double *)TIO_MALLOC (2 * n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }
   tmp_lat = tmp_lon + n;

   /* Try to simplify the polygon, but if the result is self-intersecting,
    * keep all the boundary points..
    */
   for ( ; band_km > 0.0; band_km -= 1.0)
     {
        TIO_FREE(indices);
        indices = NULL;

        if ((num_kept = simplify_dp (lon, lat, n, band_km, &indices)) < 0)
          goto return_status;

        for (i = 0; i < num_kept; i++)
          {
             int k = indices[i];
             tmp_lon[i] = (double) lon[k];
             tmp_lat[i] = (double) lat[k];
          }

        /* If necessary, close the polygon */
        if ((tmp_lon[num_kept-1] != tmp_lon[0]) || (tmp_lat[num_kept-1] != tmp_lat[0]))
          {
             tmp_lon[num_kept] = (double) tmp_lon[0];
             tmp_lat[num_kept] = (double) tmp_lat[0];
             num_kept++;
          }

        if (0 == (polygon_type = polygon_self_intersects (tmp_lon, tmp_lat, num_kept)))
          break;
     }

   if (polygon_type == 0)
     {
        for (i = 0; i < num_kept; i++)
          {
             lon[i] = (float) tmp_lon[i];
             lat[i] = (float) tmp_lat[i];
          }
     }
   else num_kept = n;  /* e.g. boundary simplification failed */

   *num = num_kept;
   *plon = lon;
   *plat = lat;

   status = 0;
return_status:
   TIO_FREE(inrqf);
   TIO_FREE(lon2d);
   TIO_FREE(lat2d);
   TIO_FREE(lon2d_bnds);
   TIO_FREE(lat2d_bnds);
   TIO_FREE(vza2d);
   TIO_FREE(indices);
   TIO_FREE(side);
   TIO_FREE(tmp_lon);
   TIO_FREE(exclude_column);
   TIO_FREE(delta_lon);

   if (status)
     {
        TIO_FREE(lon);
        TIO_FREE(lat);
     }

   return status;
}

/*{{{ Fortran-callable interface for __tio_make_lev1_bounding_polygon */
struct Bounding_Polygon_Type
{
   float *lon;
   float *lat;
   int num;
};
void __free_lev1_bounding_polygon_struct (struct Bounding_Polygon_Type *bpt)
{
   if (bpt == NULL)
     return;
   TIO_FREE(bpt->lon);
   TIO_FREE(bpt->lat);
}
int __make_lev1_bounding_polygon_struct (const int *grp, struct Bounding_Polygon_Type *bpt)
{
   return __tio_make_lev1_bounding_polygon (*grp, &bpt->num, &bpt->lon, &bpt->lat);
}

/*}}}*/

int tio_meta_set_odl_bounding_polygon (TIO_Meta_Type *meta,
                                       const float *lon, const float *lat, int num)
{
   const char *name_lons = "polygon_longitudes";
   const char *name_lats = "polygon_latitudes";
   const char *name_seqs = "polygon_sequence";
   int *seq = NULL;
   int i, status = -1;

   if (NULL == (seq = (int *)TIO_MALLOC (num * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   for (i = 0; i < num; i++)
     {
        seq[i] = i+1;
     }

   /* These variables are for ODL only, and may be be obsolete */
   if ((0 != tio_meta_set (meta, name_lons, TIO_META_TYPE_FLOAT, num, lon))
       || (0 != tio_meta_set (meta, name_lats, TIO_META_TYPE_FLOAT, num, lat))
       || (0 != tio_meta_set (meta, name_seqs, TIO_META_TYPE_INT, num, seq)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: setting boundary polygon vertex arrays", __func__);
        goto return_status;
     }

   /* mark these as ODL only */
   find_key_by_name (meta, name_lons)->odl_only = 1;
   find_key_by_name (meta, name_lats)->odl_only = 1;
   find_key_by_name (meta, name_seqs)->odl_only = 1;

   status = 0;
return_status:

   TIO_FREE(seq);
   return status;
}

static int make_acdd_geospatial_bounds (const float *lon, const float *lat, int num,
                                        char **str, float lon_range[2], float lat_range[2])
{
   char *end_str = NULL, *p = NULL;
   int len, i, n, c;

   lon_range[0] = lat_range[0] = +FLT_MAX;
   lon_range[1] = lat_range[1] = -FLT_MAX;

   /* geospatial_bounds format:
    * geospatial_bounds = POLYGON((lat0 lon0, lat1 lon1, ... lat0 lon0))
    * Coordinates of each point _must_ be ordered as lat_i lon_i.
    * For each lat lon pair, we need 9+9+1 = 19 characters.  This allows
    * for 7 significant figures, plus sign, plus decimal for each value,
    * plus a space between them.  For N vertices, we'll also have N-1 commas.
    */

   len = num*20 + 12;
   if (NULL == (*str = (char *)TIO_MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   end_str = *str + len;

   strcpy (*str, "POLYGON(");
   p = *str + 8;
   c = '(';
   for (i = 0; (p < end_str) && (i < num); i++)
     {
        float lon_i = lon[i];
        float lat_i = lat[i];
        if (lon_i < lon_range[0]) lon_range[0] = lon_i;
        if (lon_i > lon_range[1]) lon_range[1] = lon_i;
        if (lat_i < lat_range[0]) lat_range[0] = lat_i;
        if (lat_i > lat_range[1]) lat_range[1] = lat_i;
        /* Coordinates of each point _must_ be ordered as lat_i lon_i */
        if ((n = snprintf (p, end_str-p, "%c%0.4f %0.4f", c, lat_i, lon_i)) < 0)
          {
             TIO_FREE(*str);
             *str = NULL;
             return -1;
          }
        c = ',';
        p += n;
     }
   strcpy (p, "))");

   return 0;
}

int _pTIO_write_acdd_geospatial_attrs (int grp, const float *lon, const float *lat, int num)
{
   char *str = NULL;
   char *crs = GEOSPATIAL_BOUNDS_CRS;
   float lon_range[2], lat_range[2];

   if (0 != make_acdd_geospatial_bounds (lon, lat, num, &str, lon_range, lat_range))
     return -1;

   if ((0 != TIO_put_att (grp, NC_GLOBAL, "geospatial_bounds", NC_CHAR, 1 + strlen(str), str))
       || (0 != TIO_put_att (grp, NC_GLOBAL, "geospatial_bounds_crs", NC_CHAR, 1 + strlen(crs), crs))
       || (0 != TIO_put_att (grp, NC_GLOBAL, "geospatial_lon_min", NC_FLOAT, 1, &lon_range[0]))
       || (0 != TIO_put_att (grp, NC_GLOBAL, "geospatial_lon_max", NC_FLOAT, 1, &lon_range[1]))
       || (0 != TIO_put_att (grp, NC_GLOBAL, "geospatial_lat_min", NC_FLOAT, 1, &lat_range[0]))
       || (0 != TIO_put_att (grp, NC_GLOBAL, "geospatial_lat_max", NC_FLOAT, 1, &lat_range[1])))
     {
        TIO_FREE(str);
        return -1;
     }

   TIO_FREE(str);

   return 0;
}

int tio_meta_set_acdd_geospatial_bounds (TIO_Meta_Type *meta,
                                         const float *lon, const float *lat, int num)
{
   char *str = NULL;
   float lon_range[2], lat_range[2];
   int len_crs = 1 + strlen(GEOSPATIAL_BOUNDS_CRS);

   if (0 != make_acdd_geospatial_bounds (lon, lat, num, &str, lon_range, lat_range))
     return -1;

   if ((0 != tio_meta_set (meta, "geospatial_bounds", TIO_META_TYPE_CHAR, 1+strlen(str), str))
       ||(0 != tio_meta_set (meta, "geospatial_bounds_crs", TIO_META_TYPE_CHAR, len_crs, GEOSPATIAL_BOUNDS_CRS))
       ||(0 != tio_meta_set (meta, "geospatial_lon_min", TIO_META_TYPE_FLOAT, 1, &lon_range[0]))
       ||(0 != tio_meta_set (meta, "geospatial_lon_max", TIO_META_TYPE_FLOAT, 1, &lon_range[1]))
       ||(0 != tio_meta_set (meta, "geospatial_lat_min", TIO_META_TYPE_FLOAT, 1, &lat_range[0]))
       ||(0 != tio_meta_set (meta, "geospatial_lat_max", TIO_META_TYPE_FLOAT, 1, &lat_range[1])))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: setting geospatial_bounds", __func__);
        TIO_FREE(str);
        return -1;
     }

   TIO_FREE(str);
   return 0;
}

int tio_meta_set_lev1_bounding_polygon (TIO_Meta_Type *meta, int grp)
{
   float *lon=NULL, *lat=NULL;
   int num, status = -1;

   if (0 != __tio_make_lev1_bounding_polygon (grp, &num, &lon, &lat))
     return -1;

   if (0 != tio_meta_set_odl_bounding_polygon (meta, lon, lat, num))
     goto return_status;

   if (0 != tio_meta_set_acdd_geospatial_bounds (meta, lon, lat, num))
     goto return_status;

   status = 0;
return_status:
   TIO_FREE(lon);
   TIO_FREE(lat);
   return status;
}

int tio_meta_simplify_dp (const float *lon_deg, const float *lat_deg, int num,
                          float band_km, int **pindex)
{
   return simplify_dp (lon_deg, lat_deg, num, band_km, pindex);
}

#define MAX_DATETIME_KEYLEN 72
static int meta_set_datetime (TIO_Meta_Type *meta, const char *str, const char *key_prefix)
{
   char date_str[MAX_ISOTIME_LEN], time_str[MAX_ISOTIME_LEN];
   char date_key[MAX_DATETIME_KEYLEN], time_key[MAX_DATETIME_KEYLEN];
   int len, status;

   if (2 != (status = sscanf (str, "%[^T]T%[^Z]Z", date_str, time_str)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: parsing timestamp string: %s (status=%d)",
                     __func__, str, status);
        return -1;
     }

   if (((len = snprintf (date_key, MAX_DATETIME_KEYLEN, "%sdate", key_prefix)) < 0)
       || (len == MAX_DATETIME_KEYLEN))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: constructing date keyword", __func__);
        return -1;
     }

   if (((len = snprintf (time_key, MAX_DATETIME_KEYLEN, "%stime", key_prefix)) < 0)
       || (len == MAX_DATETIME_KEYLEN))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: constructing time keyword", __func__);
        return -1;
     }

   if ((0 != tio_meta_set (meta, date_key, TIO_META_TYPE_CHAR, 1+strlen(date_str), date_str))
       ||(0 != tio_meta_set (meta, time_key, TIO_META_TYPE_CHAR, 1+strlen(time_str), time_str)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: setting date/time metadata keywords", __func__);
        return -1;
     }

   return 0;
}

int tio_meta_set_datetime_range (TIO_Meta_Type *meta, int ncid)
{
   char str[MAX_ISOTIME_LEN];

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, str))
     return -1;

   if (0 != meta_set_datetime (meta, str, "begin_"))
     return -1;

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, str))
     return -1;

   if (0 != meta_set_datetime (meta, str, "end_"))
     return -1;

   return 0;
}

int tio_meta_set_datetime_range_scan (TIO_Meta_Type *meta, const TIO_Scan_Ident_Type *lst)
{
   const _pTIO_Granule_Ident_Type *beg=NULL, *end=NULL, *gid;
   double t_beg, t_end;

   if (lst == NULL)
     return -1;

   beg = lst->granule_ident;
   end = beg;

   t_beg = beg->tstart;
   t_end = beg->tend;

   for (gid = lst->granule_ident; gid != NULL; gid = gid->next)
     {
        if (gid->tstart < t_beg)
          {
             t_beg = gid->tstart;
             beg = gid;
          }

        if (gid->tend > t_end)
          {
             t_end = gid->tend;
             end = gid;
          }
     }

   if (0 != meta_set_datetime (meta, beg->tstart_str, "begin_"))
     return -1;

   if (0 != meta_set_datetime (meta, end->tend_str, "end_"))
     return -1;

   return 0;
}

int tio_meta_set_datetime_production (TIO_Meta_Type *meta)
{
   char time_str[MAX_ISOTIME_LEN];
   struct tm tm = {0};
   time_t now;
   size_t len;

   now = time(NULL);
   gmtime_r (&now, &tm);

   len = strftime (time_str, sizeof(time_str), "%Y-%m-%dT%H:%M:%SZ", &tm);
   if ((len >= sizeof(time_str)) || (len == 0))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: strftime failed", __func__);
        return -1;
     }

   if (0 != tio_meta_set (meta, "production_date_time", TIO_META_TYPE_CHAR, 1+strlen(time_str), time_str))
     return -1;

   return 0;
}

int tio_meta_set_standard (TIO_Meta_Type *meta,
                           const char *product_file_name,
                           const char *product_short_name,
                           int product_versionid,
                           const char *pge_version_string)
{
   const char *basename;

   if (NULL != (basename = strrchr (product_file_name, '/')))
     {
        basename++;
     }
   else basename = product_file_name;

   if ((0 != tio_meta_set (meta, "local_granule_id", TIO_META_TYPE_CHAR, 1+strlen(basename), basename))
       || (0 != tio_meta_set (meta, "version_id", TIO_META_TYPE_INT, 1, &product_versionid))
       || (0 != tio_meta_set (meta, "pge_version", TIO_META_TYPE_CHAR, 1+strlen(pge_version_string), pge_version_string)))
     {
        return -1;
     }

   if (product_short_name != NULL)
     {
        if (0 != tio_meta_set (meta, "shortname", TIO_META_TYPE_CHAR, 1+strlen(product_short_name), product_short_name))
          return -1;
     }

   return 0;
}
