/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define __USE_XOPEN   /* for strptime */
#include <time.h>

#include <tell.h>

#include "tio.h"
#include "_tio.h"
#include "tio_template.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2  (M_PI/2.0)
#endif

#define DEGTORAD  (M_PI/180.0)

struct TIO_Meta_Type
{
   TIO_Meta_Type *next;
   char *name;
   void *values;
   int num_values;
   int num_alloc;
   int value_type;
   int noexpand;
};

#define KEYNAME_VALID_CHARS "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_."

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
      case TIO_META_TYPE_STRING:
        break;
     }

   return 0;
}

static int meta_write_value_ascii (const TIO_Meta_Type *m, int want_num_values, FILE *fp)
{
   int i;
   int *ip = (int *)m->values;
   unsigned int *uip = (unsigned int *)m->values;
   double *dp = (double *)m->values;
   float *sp = (float *)m->values;
   char **spp = (char **)m->values;

   if (want_num_values)
     {
        (void) fprintf (fp, "%d", m->num_values);
        return 0;
     }

   if (m->num_values > 1)
     {
        fputs ("(", fp);
     }

   switch (m->value_type)
     {
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

   if (m->num_values > 1)
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

int tio_meta_ncinit (TIO_Meta_Type *meta, int grp, const char *name, int value_type)
{
   int varid = NC_GLOBAL;
   nc_type att_type;
   size_t i, att_len;
   char **atts = NULL;
   int status = 0;

   if (value_type != TIO_META_TYPE_STRING)
     {
        tell_verror (TELL_UNSUPPORTED_ERROR, "%s: reading non-string metadata attribute from netcdf",
                     __func__);
        return -1;
     }

   /* It's ok if the attribute doesn't exist */
   if (NC_NOERR != nc_inq_att (grp, varid, name, &att_type, &att_len))
     return 0;

   if (att_type != NC_STRING)
     {
        tell_verror (TELL_UNSUPPORTED_ERROR, "%s: attribute %s is not of type NC_STRING (type=%d)",
                     __func__, name, att_type);
        return -1;
     }

   if (NULL == (atts = TIO_MALLOC (att_len * sizeof(char *))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   if (NC_NOERR != nc_get_att (grp, varid, name, atts))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading attribute %s", __func__, name);
        TIO_FREE(atts);
        return -1;
     }

   for (i = 0; i < att_len; i++)
     {
        if (0 != tio_meta_append_string (meta, name, atts[i]))
          {
             status = -1;
             break;
          }
     }

   /* free space allocated on string input */
   (void) nc_free_string (att_len, (char **)atts);
   TIO_FREE(atts);

   return status;
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

static int compute_centroid (const float *lon2d, const float *lat2d, const int *inrqf, int num_pixels,
                             float *lon_centroid, float *lat_centroid)
{
   float x, y, z, phi, theta;
   int i, n;

   /* NULL means these values aren't wanted */
   if ((lon_centroid == NULL) || (lat_centroid == NULL))
     return 0;

   x = y = z = 0.0;
   n = 0;

   for (i = 0; i < num_pixels; i++)
     {
        if (inrqf[i] == 0)
          {
             float sin_theta;
             phi = lon2d[i] * DEGTORAD;
             theta = M_PI_2 - lat2d[i] * DEGTORAD;
             sin_theta = sin(theta);

             x += sin_theta * cos(phi);
             y += sin_theta * sin(phi);
             z += cos(theta);
             n++;
          }
     }

   x /= n;
   y /= n;
   z /= n;

   *lon_centroid = atan2 (y, x) / DEGTORAD;
   *lat_centroid = atan2 (z, hypot (x, y)) / DEGTORAD;

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

int __tio_make_lev1_bounding_polygon (int grp, int *num, float **plon, float **plat,
                                      float *lon_centroid, float *lat_centroid)
{
   TIO_Var_Info_Type info;
   float *vza2d=NULL, *lon2d=NULL, *lat2d=NULL, *lon=NULL, *lat=NULL;
   float *lon2d_bnds=NULL, *lat2d_bnds=NULL;
   int *inrqf=NULL, *indices=NULL;
   int *bx1=NULL, *bx2=NULL, *bs1=NULL, *bs2=NULL, *bdry=NULL, *side=NULL;
   int start[3], count[3];
   int num_steps, num_xtrack, num_pixels, max_num_boundary;
   int s, x, i, n, x_first_ok, x_last_ok, s_first_ok, s_last_ok;
   int varid, no_fill, lon_bounds_status;
   int status = -1;
   float fill_value = TIO_FILL_FLOAT;
   float vza_max_deg = 85.0; /* avoid pixels near the Earth's limb */
   int dx=8, ds=2;            /* reduce final polygon point density */

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
       || (NULL == (indices = (int *)TIO_MALLOC (2 * max_num_boundary * sizeof(int)))) )
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

   if (0 != compute_centroid (lon2d, lat2d, inrqf, num_pixels, lon_centroid, lat_centroid))
     goto return_status;

   /* Points with valid (lon,lat) coordinates have inrqf == 0.
    * We want the boundary of this region.
    * Find the boundary points by scanning each row and each column
    * to find the outermost endpoints of the region with inrqf==0.
    *
    * NOTE: The centroid ignores VZA, so the bounding polygon is clipped
    * relative to the centroid
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
#if 0
   /* The northernmost boundary points are sometimes jagged and/or
    * confusing to deal with.  The simplest solution is to just omit
    * these boundary points, letting that side of the polygon close
    * with a single line segment connecting the northernmost endpoints
    * of the eastern and western boundaries.
    */
   for (s = num_steps-1; s >= 0; s -= ds)
     {
        if (bx1[s] >= 0)
          {
             side[n] = 3;
             bdry[n] = bx1[s];
             n++;
          }
     }
#endif

   /* The polygon is assumed to be a closed shape, so there's no need
    * to explicitly close it by appending a copy of the first point */

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
             float lon_i, lat_i;
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

   *num = n;
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
   float centroid_lon;
   float centroid_lat;
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
   return __tio_make_lev1_bounding_polygon (*grp, &bpt->num, &bpt->lon, &bpt->lat,
                                            &bpt->centroid_lon, &bpt->centroid_lat);
}

/*}}}*/

int tio_meta_set_lev1_bounding_polygon_and_centroid (TIO_Meta_Type *meta, int grp)
{
   float *lon=NULL, *lat=NULL;
   int *seq = NULL;
   float centroid_lon, centroid_lat;
   int i, num, status = -1;

   if (0 != __tio_make_lev1_bounding_polygon (grp, &num, &lon, &lat, &centroid_lon, &centroid_lat))
     return -1;

   if (NULL == (seq = (int *)TIO_MALLOC (num * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   for (i = 0; i < num; i++)
     {
        seq[i] = i+1;
     }

   if ((0 != tio_meta_set (meta, "GRINGPOINTLONGITUDE", TIO_META_TYPE_FLOAT, num, lon))
       || (0 != tio_meta_set (meta, "GRINGPOINTLATITUDE", TIO_META_TYPE_FLOAT, num, lat))
       || (0 != tio_meta_set (meta, "GRINGPOINTSEQUENCENO", TIO_META_TYPE_INT, num, seq))
       || (0 != tio_meta_set (meta, "CENTROID_MEAN_LONGITUDE", TIO_META_TYPE_FLOAT, 1, &centroid_lon))
       || (0 != tio_meta_set (meta, "CENTROID_MEAN_LATITUDE", TIO_META_TYPE_FLOAT, 1, &centroid_lat))
      )
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: setting boundary polygon metadata", __func__);
        goto return_status;
     }

   status = 0;
return_status:
   TIO_FREE(lon);
   TIO_FREE(lat);
   TIO_FREE(seq);
   return status;
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

   if (((len = snprintf (date_key, MAX_DATETIME_KEYLEN, "%sDATE", key_prefix)) < 0)
       || (len == MAX_DATETIME_KEYLEN))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: constructing date keyword", __func__);
        return -1;
     }

   if (((len = snprintf (time_key, MAX_DATETIME_KEYLEN, "%sTIME", key_prefix)) < 0)
       || (len == MAX_DATETIME_KEYLEN))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: constructing time keyword", __func__);
        return -1;
     }

   if ((0 != tio_meta_set (meta, date_key, TIO_META_TYPE_STRING, 1, date_str))
       ||(0 != tio_meta_set (meta, time_key, TIO_META_TYPE_STRING, 1, time_str)))
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

   if (0 != meta_set_datetime (meta, str, "RANGEBEGINNING"))
     return -1;

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, str))
     return -1;

   if (0 != meta_set_datetime (meta, str, "RANGEENDING"))
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

   if (0 != meta_set_datetime (meta, beg->tstart_str, "RANGEBEGINNING"))
     return -1;

   if (0 != meta_set_datetime (meta, end->tend_str, "RANGEENDING"))
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

   if (0 != tio_meta_set (meta, "PRODUCTIONDATETIME", TIO_META_TYPE_STRING, 1, time_str))
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

   if ((0 != tio_meta_set (meta, "LOCALGRANULEID", TIO_META_TYPE_STRING, 1, basename))
       || (0 != tio_meta_set (meta, "VERSIONID", TIO_META_TYPE_INT, 1, &product_versionid))
       || (0 != tio_meta_set (meta, "PGEVERSION", TIO_META_TYPE_STRING, 1, pge_version_string)))
     {
        return -1;
     }

   if (product_short_name != NULL)
     {
        if (0 != tio_meta_set (meta, "SHORTNAME", TIO_META_TYPE_STRING, 1, product_short_name))
          return -1;
     }

   return 0;
}
