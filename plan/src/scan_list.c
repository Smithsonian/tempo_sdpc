/** @file scan_list.c
 *  @brief Manage a list of instrument scan parameters
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>

#include "scan_list.h"

void scan_list_entry_free (Scan_List_Entry *stt)
{
   if (stt == NULL)
     return;
   FREE(stt);
}

void scan_list_free (Scan_List_Entry *head)
{
   while (head != NULL)
     {
        Scan_List_Entry *stt = head->next;
        scan_list_entry_free (head);
        head = stt;
     }
}

Scan_List_Entry *scan_list_entry_alloc (void)
{
   Scan_List_Entry *stt = NULL;

   if (NULL == (stt = (Scan_List_Entry *) MALLOC (sizeof *stt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)stt, 0, sizeof *stt);

   return stt;
}

int scan_list_append (Scan_List_Entry **phead,
                      Scan_List_Entry *stt)
{
   Scan_List_Entry *head;

   if (phead == NULL)
     return -1;

   head = *phead;

   if (head == NULL)
     {
        *phead = stt;
        return 0;
     }

   for ( ; head != NULL; head = head->next)
     {
        if (head->next == NULL)
          {
             head->next = stt;
             return 0;
          }
     }

   return -1;
}

#define TIME_BUFSIZE 21

int scan_list_write (FILE *fp, const Scan_List_Entry *head)
{
   const Scan_List_Entry *entry;
   const char header_comment[] =
     "# tstart,duration,xstart,num_steps,step_exposure,num_repeats,timestamp\n";

   if (fprintf (fp, header_comment) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
        return -1;
     }

   for (entry = head; entry != NULL; entry = entry->next)
     {
        char buf[TIME_BUFSIZE];

        if (0 != mkjdtimestr (entry->tstart, buf, sizeof(buf)))
          return -1;

        if (fprintf (fp, "%0.6f,%0.3f,%d,%d,%0.3f,%d,\"%s\"\n",
                     entry->tstart,
                     entry->scan_duration,
                     (int) entry->xstart,
                     entry->num_steps,
                     entry->step_exposure,
                     entry->num_repeats,
                     buf) < 0)
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
             return -1;
          }
     }

   return 0;
}
