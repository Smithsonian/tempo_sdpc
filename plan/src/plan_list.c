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

Plan_List_Type *plan_list_entry_alloc (uint16_t scan_type, int region_id)
{
   Plan_List_Type *ple = NULL;

   if (NULL == (ple = (Plan_List_Type *) MALLOC (sizeof *ple)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)ple, 0, sizeof *ple);

   ple->scan_type = scan_type;
   ple->region_id = region_id;

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

static double mirror_tilt (double azimuth)
{
   /* The FOR coordinate system refers to the azimuth and elevation angular
    * coordinates in the field of regard indicating the line of sight
    * from which we want to collect photons entering through the slit.
    * To command the instrument, the flight software wants the mirror tilt
    * angle needed to access that line of sight.  Using the law of reflection,
    * the mirror tilt angle is half the azimuthal angle coordinate in the
    * field of regard.
    *
    * The azimuthal angle coordinate in the field of regard increases toward
    * the east (+X axis in a right-handed coordinate system).  The elevation
    * coordinate increases toward the south (+Y axis).
    *
    * The C&THB documentation for the SMA_MOVE command says:
    *
    * "Neglecting alignment tolerances, a motion of the scan mirror in
    * the positive X-direction moves the line of sight in the Spacecraft
    * +X direction (East) . A motion of the scan mirror in the positive
    * Y-direction moves the line of sight in the Spacecraft +Y direction
    * (South)."
    *
    * Therefore, the mirror tilt angle +X coordinate has the same sign
    * as the +X azimuthal angle in the field of regard.
    *
    * Both angles are in microradians.
    */

   return 0.5 * azimuth;
}

#define TIME_BUFSIZE 32

int plan_list_write (FILE *fp, int mark_scan_seq_start, const Plan_List_Type *head)
{
   const Plan_List_Type *entry;
   const char header_comment[] =
     "label,time,duration,mirror_x,num_steps,integration_time,repeat,timestamp,comment\n";
   double unix_epoch_jd;
   double previous_entry_tstop_tai, previous_entry_jd_utc_end;
   uint16_t scan_num;
   int num_scan_csm, num_days, scan_num_to_inr;

   unix_epoch_jd = novas_julian_date (1970,1,1,0.0);

   if (fprintf (fp, header_comment) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
        return -1;
     }

   previous_entry_tstop_tai = 0.0;
   previous_entry_jd_utc_end = 0.0;

   scan_num = 1;
   num_scan_csm = 0;
   num_days = 0;
   scan_num_to_inr = 0;

   for (entry = head; entry != NULL; entry = entry->next)
     {
        double tstart_utc, tstart_tai, fsw_xstart;
        char buf[TIME_BUFSIZE];
        int is_twilight = (entry->scan_type & TEMPO_SCAN_TYPE_NIGHTLIGHTS);
        int new_day = (previous_entry_jd_utc_end < entry->tstart);
        int i, num_scans;

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

        /* The plan is generated using an azimuthal coordinate in the field of regard,
         * but for the IOC plan, we want to write out the mirror tilt angle */
        fsw_xstart = mirror_tilt (entry->xstart);

        /* Restart scan numbering each day, (scan_num=0 is used as a fill value,
         * so we number scans from 1).
         * IMPORTANT: Every scan must have a unique value of scan_num within each day.
         * This ensures that, by combining scan_num with the satellite-local day number
         * counter, we can construct an integer scan_id label that uniquely identifies
         * each scan throughout the entire mission.
         * DO NOT change this unless you fully understand the consequences.
         */
        if (new_day)
          {
             scan_num = 1;       /* <- WARNING: Don't change the scan numbering! */
             scan_num_to_inr = 0;
             num_days++;
          }

        num_scans = (entry->num_repeats_cbm > 0) ? entry->num_repeats_cbm : 1;

        for (i = 0; i < entry->num_repeats; i++, scan_num += num_scans)
          {
             uint16_t scan_label, scan_type;
             double tstart_jd = (entry->tstart
                                 + i * entry->scan_duration/SEC_PER_DAY);

             if (0 != mkjdtimestr (tstart_jd, buf, sizeof(buf)))
               return -1;

             /* mark only the first scan of each new sequence destined for INR */
             if (is_twilight == 0) scan_num_to_inr++;
             if (((scan_num_to_inr == 1) || (entry->post_maneuver != 0))
                 && (i == 0) && (mark_scan_seq_start != 0))
               {
                  scan_type = entry->scan_type | TEMPO_SCAN_TYPE_SCAN_SEQ_START;
               }
             else scan_type = entry->scan_type;

             if (0 != tio_make_scan_label (&scan_label, scan_type, scan_num))
               return -1;

             num_scan_csm += entry->num_repeats_cbm ? 2 : 3;

             if (fprintf (fp, "%d,%0.3f,%0.3f,%0.1f,%d,%0.6f,%d,\"%s\",\n",
                          scan_label,
                          tstart_tai,
                          entry->scan_duration / num_scans,
                          fsw_xstart,
                          entry->num_steps,
                          entry->integration_time,
                          entry->num_repeats_cbm,
                          buf) < 0)
               {
                  tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
                  return -1;
               }

             tstart_tai += entry->scan_duration;
          }

        previous_entry_tstop_tai = tstart_tai;
        previous_entry_jd_utc_end = entry->jd_utc_end_safe;
     }

   fprintf (stdout, "    %5d Plan days\n", num_days);
   fprintf (stdout, "    %5d Scan CSM entries\n", num_scan_csm);
   fprintf (stdout, " => %5d Total CSM entries (estimated)\n", num_scan_csm + 5 * num_days);

   return 0;
}

void plan_stats_list_free (Plan_Stats_Type *stats)
{
   while (stats != NULL)
     {
        Plan_Stats_Type *next = stats->next;
        FREE(stats);
        stats = next;
     }
}

Plan_Stats_Type *plan_stats_list_append (Plan_Stats_Type *head, double jd_utc_beg, double jd_utc_end,
                                         double jd_utc_beg_safe, double jd_utc_end_safe)
{
   Plan_Stats_Type *p, *stats;

   if (NULL == (stats = (Plan_Stats_Type *)MALLOC (sizeof *stats)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)stats, 0, sizeof (*stats));
   stats->jd_utc_beg = jd_utc_beg;
   stats->jd_utc_end = jd_utc_end;
   stats->jd_utc_beg_safe = jd_utc_beg_safe;
   stats->jd_utc_end_safe = jd_utc_end_safe;

   for (p = head; p != NULL; p = p->next)
     {
        if (p->next == NULL)
          break;
     }

   p->next = stats;

   return stats;
}

int plan_stats_set_scan_times (Plan_Stats_Type *stats, const Plan_List_Type *head)
{
   const Plan_List_Type *entry;
   double tbeg, tend;

   tbeg = 0.0;
   tend = 0.0;

   for (entry = head; entry != NULL; entry = entry->next)
     {
        int is_twilight = (entry->scan_type & TEMPO_SCAN_TYPE_NIGHTLIGHTS);
        if (is_twilight) continue;

        if (tbeg == 0.0)
          {
             tbeg = entry->tstart;
             tend = tbeg;
          }

        tend += entry->num_repeats * entry->scan_duration/SEC_PER_DAY;
     }

   stats->radiance_scan_first_start = tbeg;
   stats->radiance_scan_last_end = tend;

   return 0;
}

static double Unix_Epoch_JD;
static int jd_to_tai (double tstamp_jd, double *tstamp_tai)
{
   double tstamp_utc;

   if (Unix_Epoch_JD == 0.0)
     {
        Unix_Epoch_JD = novas_julian_date (1970,1,1,0.0);
     }

   tstamp_utc = (tstamp_jd - Unix_Epoch_JD) * SEC_PER_DAY;
   return tio_time_utc_to_taix (tstamp_utc, tstamp_tai);
}

int plan_stats_write (const Plan_Stats_Type *stats, double min_sun_angle, const char *filename)
{
   FILE *fp = NULL;
   const Plan_Stats_Type *p;
   const char hdr[] =
     "#  scan_beg/end = Radiance scanning start/end times (ignoring maneuvers)\n"
     "#   sza_beg/end = SZA constraint times\n"
     "#  safe_beg/end = Safety constraint times\n"
     "date,scan_beg,scan_end,sza_beg,sza_end,safe_beg,safe_end\n";
   struct time_bounds
     {
        double beg;
        double end;
     }
   safe, sza, scan, delta;
   time_t epoch;
   time_t now_tt = time(NULL);
   char epoch_str[32];
   int num, status = -1;

   epoch = tio_time_taix_epoch_timet();
   if (0 != TIO_mktimestamp_str (0.0, 1, epoch_str, sizeof(epoch_str)))
     return -1;

   if (NULL == (fp = fopen (filename, "w")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening file for writing: %s", __func__, filename);
        return -1;
     }

   (void) fprintf (fp, "# Created: %s", ctime(&now_tt));
   (void) fprintf (fp, "# epoch = %s = %ld (time_t)\n", epoch_str, epoch);
   (void) fprintf (fp, hdr);

   p = stats;
   /* head node may be only a handle */
   if (p->radiance_scan_first_start == 0.0)
     p = p->next;

   num = 0;
   delta.beg = 0.0;
   delta.end = 0.0;
   for ( ; p != NULL; p = p->next)
     {
        short year, month, day;
        double hour;

        novas_cal_date (p->radiance_scan_first_start, &year, &month, &day, &hour);

        if ((0 != jd_to_tai (p->jd_utc_beg_safe, &safe.beg))
            || (0 != jd_to_tai (p->jd_utc_end_safe, &safe.end))
            || (0 != jd_to_tai (p->jd_utc_beg, &sza.beg))
            || (0 != jd_to_tai (p->jd_utc_end, &sza.end))
            || (0 != jd_to_tai (p->radiance_scan_first_start, &scan.beg))
            || (0 != jd_to_tai (p->radiance_scan_last_end, &scan.end)))
          goto close_and_return;

        if (fprintf (fp, "%4d-%02d-%02d,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
                     year, month, day,
                     scan.beg, scan.end, sza.beg, sza.end, safe.beg, safe.end)<0)
          goto close_and_return;

        delta.beg += scan.beg - safe.beg;
        delta.end += safe.end - scan.end;
        num++;
     }

   delta.beg /= num;
   delta.end /= num;
   fprintf (stderr, "Mean safety margin [min]:  morning:%0.1f, evening:%0.1f  [min_sun_angle=%0.1f deg]\n",
            delta.beg/60.0, delta.end/60.0, min_sun_angle);

   status = 0;
close_and_return:
   fclose (fp);
   return status;
}
