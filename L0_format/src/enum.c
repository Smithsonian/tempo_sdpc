/** @file enum.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Manage enumerated data types
 */

#define _GNU_SOURCE
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <search.h>  /* for hsearch_r, etc */
#include <tio.h>
#include <tell.h>

#include "enum.h"

#ifndef MALLOC
#define MALLOC malloc
#endif

#ifndef FREE
#define FREE free
#endif

typedef struct Enum_Type Enum_Type;
struct Enum_Type
{
   Enum_Type *next;
   const char *enumlist;
   int type_id;
};

struct Enum_Lookup_Type
{
   Enum_Type *enum_defs;
   struct hsearch_data enum_hash_table;
};

static Enum_Type *new_enum_type (const char *enumlist, int type_id)
{
   Enum_Type *e = NULL;
   if (NULL == (e = (Enum_Type *) MALLOC (sizeof *e)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   e->next = NULL;
   e->enumlist = enumlist;
   e->type_id = type_id;
   return e;
}

static void free_enum_defs (Enum_Type *lst)
{
   while (lst != NULL)
     {
        Enum_Type *next = lst->next;
        FREE(lst);
        lst = next;
     }
}

static int elt_lookup (Enum_Lookup_Type *elt, const char *enumlist,
                       int *type_id)
{
   ENTRY enum_query, *enum_result = NULL;

   enum_query.key = (char *)enumlist;
   if (0 != hsearch_r (enum_query, FIND, &enum_result, &elt->enum_hash_table))
     {
        Enum_Type *er = (Enum_Type *) enum_result->data;
        *type_id = er->type_id;
        return 0;
     }

   return 1;
}

int elt_is_valid (Enum_Lookup_Type *elt, const char *enumlist)
{
   int num_values;
   const char *p;

   if ((elt == NULL) || (enumlist == NULL))
     return 0;

   num_values = 0;
   p = enumlist;

   for (;;)
     {
        if (NULL == (p = strchr (p, ':')))
          break;
        num_values++;
        p++;
     }

   return (num_values > 1);
}

static int elt_insert (Enum_Lookup_Type *elt, const char *enumlist,
                       int type_id)
{
   ENTRY enum_query, *enum_result = NULL;
   Enum_Type *er = NULL;

   if (NULL == (er = new_enum_type (enumlist, type_id)))
     return -1;
   er->next = elt->enum_defs;
   elt->enum_defs = er;

   enum_query.key = (char *)enumlist;
   enum_query.data = (void *)er;
   if (0 == hsearch_r (enum_query, ENTER, &enum_result, &elt->enum_hash_table))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: hsearch_r failed entering enum = %s",
                     __func__, enumlist ? enumlist : "(null)");
        return -1;
     }

   return 0;
}

static void free_enum_table (TIO_Enum_Type *enum_table, int num)
{
   if (enum_table == NULL)
     return;
   while (num-- > 0)
     {
        FREE(enum_table[num].name);
     }
   FREE(enum_table);
}

/* allocate N+1, and assume NULL-terminated (enum_table->name=NULL marks end) */
static TIO_Enum_Type *new_enum_table (int num)
{
   TIO_Enum_Type *enum_table = NULL;
   size_t size = (num + 1) * sizeof(*enum_table);

   if (NULL == (enum_table = (TIO_Enum_Type *) MALLOC (size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)enum_table, 0, size);
   return enum_table;
}

static int parse_enumlist (const char *enumlist, int *pnum_keys,
                           TIO_Enum_Type **entries)
{
   TIO_Enum_Type *enum_table = NULL;
   const char *p;
   int k, num_keys;
   int status = -1;

   *pnum_keys = 0;
   *entries = NULL;

   /* Parse enumlist of the form 'name1:int1, name2:int2, ..'
    */

   num_keys = 0;
   p = enumlist;
   for (;;)
     {
        if (NULL == (p = strchr (p, ':')))
          break;
        num_keys++;
        p++;
     }

   if (NULL == (enum_table = new_enum_table (num_keys)))
     return -1;

   p = enumlist;
   for (k = 0; k < num_keys; k++)
     {
        TIO_Enum_Type *entry = &enum_table[k];
        int len_spaces = strspn (p, " \t");
        const char *pname = p + len_spaces;
        const char *colon = strchr (pname, ':');
        int len_name;
        char *comma;

        if (colon == NULL)
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: invalid enum definition %s",
                          __func__, pname);
             goto return_status;
          }

        len_name = colon - pname;
        if ((len_name == 0) || (len_name >= NC_MAX_NAME))
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: invalid enum definition %s",
                          __func__, p);
             goto return_status;
          }

        if (NULL == (entry->name = MALLOC (len_name+1)))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto return_status;
          }
        strncpy (entry->name, pname, len_name);
        entry->name[len_name] = 0;

        if (1 != sscanf (colon, ":%d", &entry->value))
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: error parsing enum value %s",
                          __func__, colon);
             goto return_status;
          }

        if (NULL == (comma = strchr (p, ',')))
          break;
        p = comma + 1;
     }

   *pnum_keys = num_keys;
   *entries = enum_table;

   status = 0;
return_status:
   if (status)
     {
        free_enum_table (enum_table, num_keys);
     }
   return status;
}

static int define_enum (int ncid, const char *mnemonic, int type_raw,
                        const char *enumlist, int *ptypeid)
{
   char buf[NC_MAX_NAME];
   TIO_Enum_Type *enum_table = NULL;
   int num_keys;

   if (snprintf (buf, NC_MAX_NAME, "enum_%s", mnemonic) >= NC_MAX_NAME)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: snprintf output was truncated to %d characters: buf=%s",
                     __func__, NC_MAX_NAME, buf);
        return -1;
     }

   if (-1 == parse_enumlist (enumlist, &num_keys, &enum_table))
     return -1;

   if (-1 == TIO_define_enum_table (ncid, buf, type_raw, enum_table, ptypeid))
     return -1;

   free_enum_table (enum_table, num_keys);

   return 0;
}

int elt_define (Enum_Lookup_Type *elt, int ncid, const char *mnemonic,
                const char *enumlist, int base_type, int *type_id)
{
   if (0 == elt_lookup (elt, enumlist, type_id))
     return 0;

   if (-1 == define_enum (ncid, mnemonic, base_type, enumlist, type_id))
     return -1;
   if (-1 == elt_insert (elt, enumlist, *type_id))
     return -1;

   return 0;
}

void elt_close (Enum_Lookup_Type *elt)
{
   if (elt == NULL)
     return;
   free_enum_defs (elt->enum_defs);
   elt->enum_defs = NULL;
   hdestroy_r (&elt->enum_hash_table);
   FREE(elt);
}

Enum_Lookup_Type *elt_open (size_t size)
{
   Enum_Lookup_Type *elt = NULL;
   if (NULL == (elt = (Enum_Lookup_Type *) MALLOC (sizeof *elt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   elt->enum_defs = NULL;
   memset ((char *)&elt->enum_hash_table, 0, sizeof(struct hsearch_data));

   if (0 == hcreate_r (size, &elt->enum_hash_table))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: hcreate_r failed, size=%ld (%s)",
                     __func__, size, strerror(errno));
        elt_close (elt);
        return NULL;
     }

   return elt;
}
