#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tell.h>

#include "met.h"

typedef struct Met_File_List_Type Met_File_List_Type;
struct Met_File_List_Type
{
   Met_File_List_Type *next;
   Met_File_Type *mft;
};

struct Met_List_Type
{
   unsigned int flags;
   Met_File_List_Type *list;
};

static void free_file_list1 (Met_File_List_Type *list)
{
   Met_File_Type *mft;
   if (list == NULL)
     return;
   mft = list->mft;
   if (mft)
     {
        mft->mft_close (mft);
        mft = NULL;
     }
   FREE(list);
}

void met_list_free (Met_List_Type *mlt)
{
   Met_File_List_Type *list;

   if (mlt == NULL)
     return;

   list = mlt->list;

   while (list != NULL)
     {
        Met_File_List_Type *next = list->next;
        free_file_list1 (list);
        list = next;
     }

   FREE(mlt);
}

static int append_entry (Met_List_Type *mlt, Met_File_Type *mft)
{
   Met_File_List_Type *entry = NULL;
   Met_File_List_Type *p = NULL;

   if (NULL == (entry = (Met_File_List_Type *)MALLOC (sizeof *entry)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   entry->next = NULL;
   entry->mft = mft;

   if (mlt->list != NULL)
     {
        for (p = mlt->list; p->next != NULL; p = p->next)
          ;
        p->next = entry;
     }
   else mlt->list = entry;

   return 0;
}

int met_list_add_file (Met_List_Type *mlt, const char *path)
{
   Met_File_Type *mft = NULL;

   if (NULL == (mft = met_open_file_grib2 (path, mlt->flags)))
     return -1;

   if (0 != append_entry (mlt, mft))
     {
        mft->mft_close (mft);
        return -1;
     }

   return 0;
}

int met_list_interp (Met_List_Type *mlt, float lon, float lat, Met_Value_Type *mvt)
{
   Met_File_List_Type *entry = mlt->list;
   int status = MFT_INTERP_FAIL;

   /* When the interpolation succeeds, we're done -> success
    * On a domain error, we try the next dataset.
    * On any other error, we're done -> fail
    */

   while (entry)
     {
        Met_File_Type *mft = entry->mft;
        status = mft->mft_interp (mft, lon, lat, mvt);
        if (status != MFT_INTERP_DOMAIN_ERROR)
          break;
        entry = entry->next;
     }

   return status;
}

Met_List_Type *met_list_new (int flags)
{
   Met_List_Type *mlt = NULL;

   if (NULL == (mlt = (Met_List_Type *)MALLOC (sizeof *mlt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   mlt->flags = flags;
   mlt->list = NULL;

   return mlt;
}
