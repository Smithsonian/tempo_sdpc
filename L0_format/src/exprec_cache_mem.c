#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <tell.h>
#include <ioclib.h>

#include "l0_format.h"

/* This module maintains a queue of exposure records, represented
 * by Exprec_Rec_Type. Incoming records are appended at the end of
 * the queue. Records written are out starting with the record at
 * the front of the queue. Records are removed from (front of) the
 * queue as they are written out.
 */

typedef struct Exprec_Rec_Type Exprec_Rec_Type;

struct Exprec_Rec_Type
{
   Exprec_Rec_Type *next;
   char *file;
   size_t file_index;
};

#define EXPREC_CACHE_METHOD_PRIVATE_DATA \
   size_t num_erecs_cached; \
   Exprec_Rec_Type *erec_list;
#include "exprec_cache.h"

static void free_rec1 (Exprec_Rec_Type *rec)
{
   if (rec == NULL)
     return;
   FREE(rec->file);
   FREE(rec);
}

static void free_rec_list (Exprec_Rec_Type *rec)
{
   while (rec)
     {
        Exprec_Rec_Type *next = rec->next;
        free_rec1 (rec);
        rec = next;
     }
}

static Exprec_Rec_Type *new_rec (const char *file, size_t file_index)
{
   Exprec_Rec_Type *rec = NULL;
   if (NULL == (rec = (Exprec_Rec_Type *)MALLOC (sizeof *rec)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)rec, 0, sizeof *rec);

   rec->file_index = file_index;

   if (NULL == (rec->file = strdup (file)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
        free_rec1 (rec);
        return NULL;
     }

   return rec;
}

static int append_rec (Exprec_Rec_Type *head, Exprec_Rec_Type *rec)
{
   Exprec_Rec_Type *r;

   for (r = head; r != NULL; r = r->next)
     {
        if (r->next == NULL)
          {
             r->next = rec;
             return 0;
          }
     }

   tell_verror (TELL_RUNTIME_ERROR, "%s: appending exposure record: file_index=%ld %s",
                __func__, rec->file_index, rec->file);

   return -1;
}

static int cache_erec (Exprec_Cache_Method_Type *cmt, const char *file, size_t file_index)
{
   Exprec_Rec_Type *rec = NULL;

   if (NULL == (rec = new_rec (file, file_index)))
     return -1;

   /* append new record to the end of the queue */
   if (cmt->erec_list == NULL)
     {
        cmt->erec_list = rec;
     }
   else if (0 != append_rec (cmt->erec_list, rec))
     {
        free_rec1(rec);
        return -1;
     }

   cmt->num_erecs_cached++;

   tell_vinfo (2, "cached erec: %ld file_index=%ld %s", cmt->num_erecs_cached, file_index, file);

   return 0;
}

static int cache_num_recs (Exprec_Cache_Method_Type *cmt, size_t *num_erecs_cached)
{
   *num_erecs_cached = cmt->num_erecs_cached;
   return 0;
}

static int cache_close (Exprec_Cache_Method_Type *cmt)
{
   (void) cmt;
   return 0;
}

static int cache_open (Exprec_Cache_Method_Type *cmt)
{
   (void) cmt;
   return 0;
}

static IOCSDPC_Exprec_Type *open_erec (Exprec_Rec_Type *rec)
{
   IOCSDPC_Common_Header_Type chdr = {0};
   IOCSDPC_Exprec_Type *erec = NULL;
   size_t k;
   int fd;

   if (-1 == (fd = iocsdpc_open_file_read (rec->file, 0, &chdr)))
     {
        tell_vlog (TELL_MSGTYPE_ERROR, 0, "%s: opening file: %s", __func__, rec->file);
        return NULL;
     }

   if (NULL == (erec = iocsdpc_exprec_fdopen_read (rec->file, fd, &chdr)))
     {
        tell_vlog (TELL_MSGTYPE_ERROR, 0, "%s: reading exprec from file: %s", __func__, rec->file);
        ioclib_fd_close (fd);
        return NULL;
     }

   if (rec->file_index == 0)
     return erec;

   for (k = 1; /* until EOF */ ; k++)
     {
        if (1 != iocsdpc_exprec_open_next (erec))
          return NULL;
        if (k == rec->file_index)
          break;
     }

   return erec;
}

/* get the first (oldest) record in the queue */
static IOCSDPC_Exprec_Type *cache_erec_get (Exprec_Cache_Method_Type *cmt)
{
   Exprec_Rec_Type *rec = cmt->erec_list;

   if (NULL == rec)
     return NULL;

   tell_vinfo (3, "get cached erec: file_index=%ld %s", rec->file_index, rec->file);

   return open_erec (rec);
}

/* retrieve the path to the file containing the first (oldest) record in the queue */
static int cache_erec_path (Exprec_Cache_Method_Type *cmt, char *buf, size_t buflen)
{
   Exprec_Rec_Type *rec = cmt->erec_list;
   size_t len;

   if (NULL == rec)
     return -1;

   strncpy (buf, rec->file, buflen);
   len = strlen (rec->file);
   if (len >= buflen) buf[buflen-1] = 0;

   return len;
}

/* The oldest record is bad.  All we can do is complain. */
static int cache_erec_bad (Exprec_Cache_Method_Type *cmt)
{
  Exprec_Rec_Type *rec = cmt->erec_list;

   if (NULL == rec)
     return -1;

   tell_vinfo (0, "bad exposure record: file=%s file_index=%ld",
               rec->file, rec->file_index);
   return 0;
}

/* We're done with the oldest record, so drop it from the queue */
static int cache_erec_done (Exprec_Cache_Method_Type *cmt)
{
   Exprec_Rec_Type *head = cmt->erec_list;
   Exprec_Rec_Type *r;

   if (head == NULL)
     return -1;

   r = head->next;
   tell_vinfo (2, "del cached erec: file_index=%ld %s",
               head->file_index, head->file);
   free_rec1 (head);
   cmt->erec_list = r;

   cmt->num_erecs_cached--;

   return 0;
}

static void cache_delete (Exprec_Cache_Method_Type *cmt)
{
   free_rec_list (cmt->erec_list);
   FREE(cmt);
}

Exprec_Cache_Method_Type *open_erec_cache_mem (config_t *cfg)
{
   Exprec_Cache_Method_Type *cmt = NULL;

   (void) cfg;

   if (NULL == (cmt = (Exprec_Cache_Method_Type *)MALLOC (sizeof *cmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)cmt, 0, sizeof *cmt);

   cmt->cache_erec = cache_erec;
   cmt->cache_num_recs = cache_num_recs;
   cmt->cache_open = cache_open;
   cmt->cache_close = cache_close;
   cmt->cache_erec_get = cache_erec_get;
   cmt->cache_erec_path = cache_erec_path;
   cmt->cache_erec_bad = cache_erec_bad;
   cmt->cache_erec_done = cache_erec_done;
   cmt->cache_delete = cache_delete;

   return cmt;
}

