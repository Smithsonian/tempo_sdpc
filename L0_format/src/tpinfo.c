/** @file tpinfo.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Search telemetry point definitions
 */

#include <errno.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <csv.h>

#include <tell.h>

#include "tpinfo.h"

#ifndef MALLOC
#define MALLOC malloc
#endif

#ifndef FREE
#define FREE free
#endif

#define CSV_BUFSIZE 1024
#define TPINFO_MALLOC_ERROR 1

typedef struct
{
   size_t fields;
   size_t rows;
}
Count_Type;

struct TPInfo_Type
{
   char **mnemonics;
   size_t *sort_index;
   TPFields_Type *fields;
   size_t num_rows;
   size_t row;
   size_t column;
   int error_status;
};

typedef struct
{
   const char *field_name;
   size_t column;
   size_t tpinfo_offset;
}
CSV_Map_Type;
#define CSV_FIELD(name)     {#name,-1,offsetof(TPFields_Type,name)}
#define CSV_FIELD_TABLE_END {NULL,-1,0}

static CSV_Map_Type CSV_Map[] =
{
   CSV_FIELD(mnemonic),
   CSV_FIELD(synopsis),
   CSV_FIELD(units),
   CSV_FIELD(enumlist),
   CSV_FIELD_TABLE_END
};

static void free_tpfields_type_fields (TPFields_Type *info)
{
   if (info == NULL)
     return;
   FREE(info->mnemonic);
   FREE(info->synopsis);
   FREE(info->units);
   FREE(info->enumlist);
}

static void free_tpinfo_type (TPInfo_Type *tp)
{
   size_t i;
   if (tp == NULL)
     return;
   for (i = 0; i < tp->num_rows; i++)
     {
        free_tpfields_type_fields (&tp->fields[i]);
     }
   FREE(tp->mnemonics);
   FREE(tp->fields);
   FREE(tp->sort_index);
   FREE(tp);
}

static TPInfo_Type *new_tpinfo_type (const Count_Type *ct)
{
   TPInfo_Type *tp = NULL;
   size_t nrows = ct->rows;

   if (NULL == (tp = (TPInfo_Type *) MALLOC (sizeof *tp)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)tp, 0, sizeof *tp);

   if ((NULL == (tp->mnemonics = (char **) MALLOC (nrows * sizeof(char *))))
       || (NULL == (tp->sort_index = (size_t *) MALLOC (nrows * sizeof(size_t)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_tpinfo_type (tp);
        return NULL;
     }
   memset ((char *)tp->mnemonics, 0, nrows * sizeof (char *));
   tp->mnemonics[0] = "";   /* No real mnemonic will appear at row 0 */

   if (NULL == (tp->fields = (TPFields_Type *) MALLOC (nrows * sizeof(TPFields_Type))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_tpinfo_type (tp);
        return NULL;
     }
   memset ((char *)tp->fields, 0, nrows * sizeof(TPFields_Type));

   tp->num_rows = nrows;
   tp->row = 0;
   tp->column = 0;
   tp->error_status = 0;

   return tp;
}

static char *strdup_plus_null (const char *s, size_t len)
{
   char *dup;
   if (s == NULL)
     return NULL;
   if (NULL == (dup = (char *) MALLOC (len+1)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memcpy (dup, s, len);
   *(dup + len) = 0;  /* ensure null-terminated */
   return dup;
}

static void csv_field_header (void *field, size_t len_field, void *client_data)
{
   TPInfo_Type *tp = (TPInfo_Type *)client_data;
   char *field_value = (char *)field;
   CSV_Map_Type *map;

   for (map = CSV_Map; map->field_name != NULL; map++)
     {
        if (strncmp (map->field_name, field_value, len_field) == 0)
          {
             map->column = tp->column;
             break;
          }
     }

   tp->column++;
}

static void csv_field_values (void *field, size_t len_field, void *client_data)
{
   TPInfo_Type *tp = (TPInfo_Type *)client_data;
   char *field_value = (char *)field;
   size_t mnemonic_column = CSV_Map[0].column;
   CSV_Map_Type *map;

   for (map = CSV_Map; map->field_name != NULL; map++)
     {
        if (map->column == tp->column)
          {
             TPFields_Type *info;
             char *value = strdup_plus_null (field_value, len_field);
             if ((value == NULL) && (field_value != NULL))
               {
                  tp->error_status = TPINFO_MALLOC_ERROR;
                  return;
               }
             info = &tp->fields[tp->row];
             *(char **)((char *)info + map->tpinfo_offset) = value;
             if (tp->column == mnemonic_column)
               {
                  tp->mnemonics[tp->row] = value;
               }
             break;
          }
     }

   tp->column++;
}

static void csv_field_callback (void *column, size_t len_column, void *client_data)
{
   TPInfo_Type *tp = (TPInfo_Type *)client_data;
   if (tp->row)
     csv_field_values (column, len_column, client_data);
   else
     csv_field_header (column, len_column, client_data);
}

static void csv_row_callback (int c, void *client_data)
{
   TPInfo_Type *tp = (TPInfo_Type *)client_data;
   (void) c;
   tp->column = 0;
   tp->row++;
}

static void csv_field_count (void *field, size_t len_field, void *cl)
{
   Count_Type *ct = (Count_Type *)cl;
   (void) field; (void) len_field;
   ct->fields++;
}

static void csv_row_count (int c, void *cl)
{
   Count_Type *ct = (Count_Type *)cl;
   (void) c;
   ct->rows++;
}

static int csv_parse1 (const char *file, void cb1(void *, size_t, void *),
                       void cb2(int, void *), void *client_data)
{
   FILE *fp;
   size_t bytes_read;
   struct csv_parser csvp;
   unsigned char options = CSV_EMPTY_IS_NULL;
   int status = -1;

   if (0 != csv_init (&csvp, options))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: csv_init failed", __func__);
        goto cleanup_and_return;
     }

   if (NULL == (fp = fopen (file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading (%s)",
                     __func__, file, strerror(errno));
        goto cleanup_and_return;
     }

   for (;;)
     {
        char buf[CSV_BUFSIZE];
        size_t bytes_parsed;

        bytes_read = fread (buf, 1, sizeof(buf), fp);
        if (bytes_read > 0)
          {
             bytes_parsed = csv_parse (&csvp, buf, bytes_read, cb1, cb2, client_data);
             if (bytes_parsed != bytes_read)
               {
                  tell_verror (TELL_APPLICATION_ERROR,
                               "%s: csv_parse failed: file=%s (%s)",
                               __func__, file, csv_strerror(csv_error (&csvp)));
                  (void) fclose (fp);
                  goto cleanup_and_return;
               }
          }

        if (feof(fp))
          break;
     }

   if (ferror(fp))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %s", __func__, file);
        goto cleanup_and_return;
     }

   (void) fclose (fp);

   status = 0;
cleanup_and_return:
   csv_fini (&csvp, cb1, cb2, client_data);
   csv_free (&csvp);

   return status;
}

static const TPInfo_Type *Temp_tp;

static int compar_mnemonics (const void *va, const void *vb)
{
   const int *ia = va;
   const int *ib = vb;
   return strcmp (Temp_tp->mnemonics[*ia], Temp_tp->mnemonics[*ib]);
}

static int define_sort_index (TPInfo_Type *tp)
{
   size_t i, nrows = tp->num_rows;

   for (i = 0; i < nrows; i++)
     {
        tp->sort_index[i] = i;
     }

   Temp_tp = tp;
   qsort (tp->sort_index, nrows, sizeof(size_t), compar_mnemonics);
   Temp_tp = NULL;

   return 0;
}

TPInfo_Type *tpinfo_init (const char *file)
{
   TPInfo_Type *tp = NULL;
   Count_Type ct;

   ct.rows = 0;
   ct.fields = 0;
   if (-1 == csv_parse1 (file, csv_field_count, csv_row_count, &ct))
     return NULL;

   if (NULL == (tp = new_tpinfo_type (&ct)))
     return NULL;

   if (-1 == csv_parse1 (file, csv_field_callback, csv_row_callback, tp))
     goto cleanup_and_return;

   if (tp->error_status)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: parsing %s", __func__, file);
        goto cleanup_and_return;
     }

   if (-1 == define_sort_index (tp))
     goto cleanup_and_return;

   return tp;

cleanup_and_return:
   free_tpinfo_type (tp);
   return NULL;
}

void tpinfo_free (TPInfo_Type *tp)
{
   free_tpinfo_type (tp);
}

void tpinfo_free_tpfields (TPFields_Type *info)
{
   FREE(info);
}

static TPFields_Type *new_tpfields_type (void)
{
   TPFields_Type *info = NULL;
   if (NULL == (info = (TPFields_Type *) MALLOC (sizeof *info)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)info, 0, sizeof *info);
   return info;
}

static int find_mnemonic (const void *vkey, const void *vmember)
{
   const char *key = vkey;
   const size_t *ip = vmember;
   const char *memb = Temp_tp->mnemonics[*ip];
   return strcmp (key, memb);
}

TPFields_Type *tpinfo_find (const TPInfo_Type *tp, const char *mnemonic,
                            TPFields_Type *fieldsp)
{
   TPFields_Type *fields = NULL;
   void *found;
   size_t index;

   if (tp == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: got NULL pointer", __func__);
        return NULL;
     }

   /* assume mnemonics are sorted in alpha-order */
   Temp_tp = tp;
   found = bsearch (mnemonic, tp->sort_index, tp->num_rows,
                    sizeof(size_t), find_mnemonic);
   Temp_tp = NULL;

   /* not found */
   if (found == NULL)
     {
        if (fieldsp)
          {
             memset ((char *)fieldsp, 0, sizeof *fieldsp);
          }
        return NULL;
     }

   index = *(size_t *)found;

   /* found, and result struct provided */
   if (fieldsp)
     {
        /* struct copy */
        *fieldsp = tp->fields[index];
        return fieldsp;
     }

   /* found, and result struct allocated */
   if (NULL == (fields = new_tpfields_type ()))
     return NULL;

   /* struct copy */
   *fields = tp->fields[index];

   return fields;
}
