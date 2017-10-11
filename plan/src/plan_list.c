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
#include <tio_template.h>
#include <tell.h>

#include "plan_list.h"

static int jd_tempo_epoch (double *epoch_jd)
{
   struct tm tm;
   const char *epoch_str = TIO_TIME_REFERENCE_STRING;
   short int year, month, day;
   double hour = 0.0;

   if (0 != TIO_parse_timestr (epoch_str, &tm))
     return -1;

   year = 1900 + tm.tm_year;
   month = 1 + tm.tm_mon;
   day = tm.tm_mday;

   *epoch_jd = novas_julian_date (year, month, day, hour);
   return 0;
}

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

int plan_list_write (FILE *fp, double (*mirror_tilt)(double),
                     const Plan_List_Type *head)
{
   const Plan_List_Type *entry;
   const char header_comment[] =
     "time,duration,mirror_x,num_steps,integration_time,timestamp\n";
   double epoch;

   if (0 != jd_tempo_epoch (&epoch))
     return -1;

   if (fprintf (fp, header_comment) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
        return -1;
     }

   for (entry = head; entry != NULL; entry = entry->next)
     {
        char buf[TIME_BUFSIZE];
        int i;

        for (i = 0; i < entry->num_repeats; i++)
          {
             double fsw_xstart;
             double tstart = (entry->tstart
                              + i * entry->scan_duration/SEC_PER_DAY);

             if (0 != mkjdtimestr (tstart, buf, sizeof(buf)))
               return -1;

             if (mirror_tilt)
               fsw_xstart = mirror_tilt (entry->xstart);
             else
               fsw_xstart = entry->xstart;

             if (fprintf (fp, "%0.14e,%0.3f,%d,%d,%0.3f,\"%s\"\n",
                          (tstart - epoch)*SEC_PER_DAY,
                          entry->scan_duration,
                          (int) fsw_xstart,
                          entry->num_steps,
                          entry->integration_time,
                          buf) < 0)
               {
                  tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
                  return -1;
               }
          }
     }

   return 0;
}
