/** @file plan_list.c
 *  @brief Manage a list of instrument scan parameters
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <libconfig.h>
#include <libnovas.h>
#include <tio.h>
#include <tio_template.h>
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

#define TIME_BUFSIZE 32

int plan_list_write (FILE *fp, double (*mirror_tilt)(double),
                     const Plan_List_Type *head)
{
   const Plan_List_Type *entry;
   const char header_comment[] =
     "time,duration,mirror_x,num_steps,integration_time,timestamp\n";
   double unix_epoch_jd;

   unix_epoch_jd = novas_julian_date (1970,1,1,0.0);

   if (fprintf (fp, header_comment) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
        return -1;
     }

   for (entry = head; entry != NULL; entry = entry->next)
     {
        double tstart_utc, tstart_tai, fsw_xstart;
        char buf[TIME_BUFSIZE];
        int i;

        tstart_utc = (entry->tstart - unix_epoch_jd) * SEC_PER_DAY;
        if (0 != tio_time_utc_to_tempo (tstart_utc, &tstart_tai))
          return -1;

        if (mirror_tilt)
          fsw_xstart = mirror_tilt (entry->xstart);
        else
          fsw_xstart = entry->xstart;

        for (i = 0; i < entry->num_repeats; i++)
          {
             double tstart_jd = (entry->tstart
                                 + i * entry->scan_duration/SEC_PER_DAY);

             if (0 != mkjdtimestr (tstart_jd, buf, sizeof(buf)))
               return -1;

             if (fprintf (fp, "%0.16e,%0.3f,%0.1f,%d,%0.3f,\"%s\"\n",
                          tstart_tai,
                          entry->scan_duration,
                          fsw_xstart,
                          entry->num_steps,
                          entry->integration_time,
                          buf) < 0)
               {
                  tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
                  return -1;
               }

             tstart_tai += entry->scan_duration;
          }
     }

   return 0;
}
