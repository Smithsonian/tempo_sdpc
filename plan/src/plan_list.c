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

#define SCAN_START_TIME_ERROR_TOLERANCE_SEC (1.e-2)
#define SEC_PER_HOUR  (3600.0)

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

Plan_List_Type *plan_list_entry_alloc (uint16_t scan_type)
{
   Plan_List_Type *ple = NULL;

   if (NULL == (ple = (Plan_List_Type *) MALLOC (sizeof *ple)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)ple, 0, sizeof *ple);

   ple->scan_type = scan_type;

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
     "label,time,duration,mirror_x,num_steps,integration_time,timestamp\n";
   double unix_epoch_jd;
   double previous_entry_tstop_tai, previous_entry_jd_utc_end;
   uint16_t last_scan_type, scan_num;

   unix_epoch_jd = novas_julian_date (1970,1,1,0.0);

   if (fprintf (fp, header_comment) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
        return -1;
     }

   previous_entry_tstop_tai = 0.0;
   previous_entry_jd_utc_end = 0.0;
   last_scan_type = head->scan_type;

   scan_num = 1;

   for (entry = head; entry != NULL; entry = entry->next)
     {
        double tstart_utc, tstart_tai, fsw_xstart;
        char buf[TIME_BUFSIZE];
        int new_scan_type = (last_scan_type != entry->scan_type);
        int new_day = (previous_entry_jd_utc_end < entry->tstart);
        int i;

        tstart_utc = (entry->tstart - unix_epoch_jd) * SEC_PER_DAY;
        if (0 != tio_time_utc_to_taix (tstart_utc, &tstart_tai))
          return -1;

        if (previous_entry_tstop_tai - tstart_tai > SCAN_START_TIME_ERROR_TOLERANCE_SEC)
          {
             tell_verror (TELL_INTERNAL_ERROR,
                          "%s: inconsistent scan start times: tstart(TAI) =%0.16e/(%f jd_utc) is %e sec BEFORE the end of the previous scan",
                          __func__, tstart_tai, entry->tstart, previous_entry_tstop_tai - tstart_tai);
             return -1;
          }

        if (mirror_tilt)
          fsw_xstart = mirror_tilt (entry->xstart);
        else
          fsw_xstart = entry->xstart;

        /* Restart scan numbering each day, and whenever scan_type changes.
         * (scan_num=0 is used as a fill value, so we number scans from 1) */
        if (new_day || new_scan_type)
          {
             scan_num = 1;
          }

        for (i = 0; i < entry->num_repeats; i++, scan_num++)
          {
             uint16_t scan_label;
             double tstart_jd = (entry->tstart
                                 + i * entry->scan_duration/SEC_PER_DAY);

             if (0 != mkjdtimestr (tstart_jd, buf, sizeof(buf)))
               return -1;

             if (0 != tio_make_scan_label (&scan_label, entry->scan_type, scan_num))
               return -1;

             if (fprintf (fp, "%d,%0.3f,%0.3f,%0.1f,%d,%0.3f,\"%s\"\n",
                          scan_label,
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

        previous_entry_tstop_tai = tstart_tai;
        previous_entry_jd_utc_end = entry->jd_utc_end_safe;
        last_scan_type = entry->scan_type;
     }

   return 0;
}
