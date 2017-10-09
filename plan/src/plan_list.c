/** @file plan_list.c
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

#include "plan_list.h"

void plan_list_entry_free (Plan_List_Type *ple)
{
   if (ple == NULL)
     return;
   FREE(ple);
}

void plan_list_free (Plan_List_Type *head)
{
   while (head != NULL)
     {
        Plan_List_Type *ple = head->next;
        plan_list_entry_free (head);
        head = ple;
     }
}

Plan_List_Type *plan_list_entry_alloc (void)
{
   Plan_List_Type *ple = NULL;

   if (NULL == (ple = (Plan_List_Type *) MALLOC (sizeof *ple)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)ple, 0, sizeof *ple);

   return ple;
}

int plan_list_append (Plan_List_Type **phead,
                      Plan_List_Type *ple)
{
   Plan_List_Type *head;

   if (phead == NULL)
     return -1;

   head = *phead;

   if (head == NULL)
     {
        *phead = ple;
        return 0;
     }

   for ( ; head != NULL; head = head->next)
     {
        if (head->next == NULL)
          {
             head->next = ple;
             return 0;
          }
     }

   return -1;
}

#define TIME_BUFSIZE 21

int plan_list_write (FILE *fp, const Plan_List_Type *head)
{
   const Plan_List_Type *entry;
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
