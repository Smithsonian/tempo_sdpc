#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Adapted from iocdbtime.c by John E. Davis <jedavis@cfa.harvard.edu> */

#include <tell.h>
#include "tio.h"
#include "_tio.h"

/* To simplify anomaly resolution, choose the TEMPO epoch to be the
   same as the instrument on-board time epoch.  The CDS time code that
   appears in the secondary header of the CCSDS packets is referenced
   to the OBT epoch. Requirement GIF-297 in Ball document 2419038,
   "Remote Electronics Subsystem (RES) Specification for the
   Instrument Space Products Geostationary UV-VIS Dispersive Imaging
   Spectrometers", states that "The on-board time (OBT) epoch shall be
   12:00, January 1, 2000 UTC.".

   Note that the CDS time code epoch is offset from the J2000.0
   standard epoch by 64.184 seconds.

   Following USNO Circular No. 179, 2005 Oct 20, page 9, the
   J2000.0 standard epoch is defined as follows:
   1) J2000.0 is 2000 January 1, 12h TT (JD 2451545.0 TT) at the
      geocenter, where TT = Terrestrial Time.
   2) TT = TAI + 32.184 sec, where TAI = International Atomic Time.
   3) Between 1999-2005, TAI = UTC + 32 sec, where
       UTC = Coordinated Universal Time.
   Therefore, on Jan 1, 2000, the UTC of the epoch is:
    UTC = TAI - 32
        = (TT - 32.184) - 32
        = 12:00:00 - 64.184
        = 11:59:00 -  4.184
        = 11:58:55.816
*/

#ifndef LEAP_SECONDS_LIST
# define LEAP_SECONDS_LIST "/usr/share/zoneinfo/leap-seconds.list"
#endif
static const char *_pTIO_Leap_Sec_File = LEAP_SECONDS_LIST;

/* The code here assumes that the Unix time is synchronized to NTP servers.
 * The leapsecond table from
 *  <https://www.ietf.org/timezones/data/leap-seconds.list>
 * is used.  This table contains a number of comments, some metadata, and
 * two columns that give the UTC offset from TAI as a function of NTP.
 * Specifically, the second column gives the value of TAI-UTC at the NTP
 * time in the first column.  The NTP time differs from the synchronized
 * Unix or POSIX time by 2208988800 seconds.
 */

/* The NTP EPOCH is at 1 January 1900, 00:00:00. */
#define UNIX_EPOCH_NTP (2208988800LL)
/* Note: (70*365 + 17leapdays)*86400 = 2208988800 */

/* The TEMPO Epoch is at 12:00, January 1, 2000 UTC */
#define TEMPO_EPOCH_TIME_T ((time_t)946728000L)
#define TEMPO_EPOCH_TAI    ((time_t)946728032L)

/* The present assumption is that the TEMPO time represents the number of
 * actual elapsed seconds since the TEMPO epoch.  To convert this to UTC,
 * the number of elapsed leap seonds must be accounted for.
 */

typedef struct
{
   time_t unix_time;
   time_t tai_time;
   int tai_utc_offset;
}
Leap_Second_Entry_Type;

typedef struct
{
   Leap_Second_Entry_Type *entries;
   unsigned int num_entries;
   long long expiration_time;
}
Leap_Second_Table_Type;

static Leap_Second_Table_Type *Leap_Second_Table;

static void free_leap_second_table (Leap_Second_Table_Type *lstt)
{
   if (lstt == NULL) return;
   TIO_FREE (lstt->entries);
   TIO_FREE (lstt);
}

static int resize_leap_second_table (Leap_Second_Table_Type *lstt, unsigned int new_size)
{
   Leap_Second_Entry_Type *entries;

   entries = (Leap_Second_Entry_Type *)TIO_REALLOC(lstt->entries, new_size*sizeof(Leap_Second_Entry_Type));
   if (entries == NULL)
     return -1;
   lstt->entries = entries;

   return 0;
}

int tio_fgets (char **linep, size_t *lenp, FILE *fp)
{
   char *line = NULL;
   size_t len = 0;
   size_t dlen = 128;

   if (fp == NULL)
     {
        *linep = NULL;
        if (lenp != NULL) *lenp = 0;
        return -1;
     }

   if (feof (fp))
     {
        *linep = NULL;
        if (lenp != NULL) *lenp = 0;
        return 0;
     }

   while (1)
     {
        char *line1;
        size_t len1;

        line1 = (char *)TIO_REALLOC (line, len + dlen);
        if (line1 == NULL)
          goto return_error;
        line = line1;

        line1 += len;
        if (NULL == fgets (line1, dlen, fp))
          {
             if (feof (fp))
               {
                  if (lenp != NULL) *lenp = len;
                  if (len == 0)
                    {
                       TIO_FREE (line);
                       *linep = NULL;
                       return 0;
                    }
                  *linep = line;
                  return 1;
               }

             tell_verror (TELL_IO_READ_ERROR, "fgets failed: %s", strerror(errno));
             goto return_error;
          }

        len1 = strlen (line1);
        len += len1;

        if ((1 + len1 < dlen)
            || (line1[len1-1] == '\n'))
          {
             if (lenp != NULL) *lenp = len;
             *linep = line;
             return 1;
          }
     }

return_error:

   TIO_FREE (line);       /* NULL ok */
   if (lenp != NULL) *lenp = 0;
   *linep = NULL;
   return -1;
}

static int read_leapsecond_table (const char *file)
{
   Leap_Second_Table_Type *lstt;
   unsigned int num_allocated, lineno;
   FILE *fp;

   lstt = (Leap_Second_Table_Type *)TIO_MALLOC (sizeof(Leap_Second_Table_Type));
   if (lstt == NULL)
     return -1;
   memset (lstt, 0, sizeof(Leap_Second_Table_Type));

   fp = fopen (file, "r");
   if (fp == NULL)
     {
        free_leap_second_table (lstt);
        tell_verror (TELL_IO_OPEN_ERROR, "Unable to open leapsecond table %s", file);
        return -1;
     }

   num_allocated = 0;
   lineno = 0;
   while (1)
     {
        Leap_Second_Entry_Type *e;
        long long ntpsecs;
        char *line;
        int tai_utc_offset;
        int status;

        status = tio_fgets (&line, NULL, fp);
        if (status == -1)
          goto return_error;
        if (status == 0)
          break;

        lineno++;

        if ((*line == '#') || (*line == '\n') || (*line == 0))
          {
             if ((*line == '#') && (line[1] == '@'))
               {
                  if (1 != sscanf (line+2, "%lld", &lstt->expiration_time))
                    {
                       tell_verror (TELL_INVALID_DATA_ERROR, "Failed to parse the leapsecond file expiration time on line %u of %s",
                                    lineno, file);
                       TIO_FREE (line);
                       goto return_error;
                    }
                  lstt->expiration_time -= UNIX_EPOCH_NTP;
               }
             TIO_FREE (line);
             continue;
          }

        if (2 != sscanf (line, "%lld %d", &ntpsecs, &tai_utc_offset))
          {
             tell_verror (TELL_INVALID_DATA_ERROR, "The leapsecond file %s appears be be corrupt on line %u",
                          file, lineno);
             TIO_FREE (line);
             goto return_error;
          }

        if (num_allocated == lstt->num_entries)
          {
             num_allocated += 10;
             if (-1 == resize_leap_second_table (lstt, num_allocated))
               {
                  TIO_FREE (line);
                  goto return_error;
               }
          }

        e = lstt->entries + lstt->num_entries;
        e->unix_time = ntpsecs - UNIX_EPOCH_NTP;
        e->tai_time = e->unix_time + tai_utc_offset;
        e->tai_utc_offset = tai_utc_offset;

        lstt->num_entries++;
        TIO_FREE (line);
     }

   (void) fclose (fp);
   free_leap_second_table (Leap_Second_Table);
   Leap_Second_Table = lstt;
   return 0;

return_error:
   (void) fclose (fp);
   free_leap_second_table (lstt);
   return -1;
}

double tio_tempo_epoch_timet (void)
{
   return (double) TEMPO_EPOCH_TIME_T;
}

static int initialize (void)
{
   if (-1 == read_leapsecond_table (_pTIO_Leap_Sec_File))
     return -1;

   return 0;
}

int tio_time_tempo_to_utc_caldate
(double tempo_time, int *year, int *month, int *day, double *hour)
{
   struct tm tm;
   double utc;
   time_t tt;

   if (0 != tio_time_tempo_to_utc (tempo_time, &utc))
     return -1;

   tt = (time_t) utc;
   if (NULL == gmtime_r (&tt, &tm))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: gmtime_r failed: tt=%ld",
                     __func__, tt);
        return -1;
     }

   *year = 1900 + tm.tm_year;
   *month = 1 + tm.tm_mon;
   *day = tm.tm_mday;
   *hour = (tm.tm_hour
            + (1.0/60)*(tm.tm_min
                        + (1.0/60)*(tm.tm_sec + (utc - tt))));
   return 0;
}

int tio_time_tempo_to_utc (double tempo_time, double *utc_time)
{
   Leap_Second_Table_Type *lstt;
   Leap_Second_Entry_Type *entries;
   double tai;
   unsigned int ilo, ihi;

   lstt = Leap_Second_Table;
   if (lstt == NULL)
     {
        if (-1 == initialize ())
          return -1;
        lstt = Leap_Second_Table;
     }

   if (tempo_time < 0.0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "tio_time_tempo_to_utc: tempo time value is less than 0");
        return -1;
     }

   tai = (TEMPO_EPOCH_TAI + tempo_time);

   entries = lstt->entries;
   ilo = 0;
   ihi = lstt->num_entries;

   if ((tai < entries->tai_time)
       || (tai >= entries[ihi-1].tai_time))
     {
        if (tai < entries->tai_time)
          *utc_time = tai;
        else
          *utc_time = tai - entries[ihi-1].tai_utc_offset;

        if (*utc_time >= lstt->expiration_time)
          tell_vwarn (0, "time value out of range, consider updating the leapsecond table");

        return 0;
     }

   while (ilo + 1 < ihi)
     {
        unsigned int i;

        i = ilo + (ihi - ilo)/2;
        if (tai >= entries[i].tai_time)
          ilo = i;
        else
          ihi = i;
     }

   *utc_time = tai - entries[ilo].tai_utc_offset;
   return 0;
}

int tio_time_tempo_to_tai (double tempo_time, double *tai_time)
{
   *tai_time = tempo_time + TEMPO_EPOCH_TAI;
   return 0;
}

int tio_time_utc_to_tempo (double utc_time, double *tempo_time)
{
   Leap_Second_Table_Type *lstt;
   Leap_Second_Entry_Type *entries;
   double tai;
   unsigned int ilo, ihi;

   if (utc_time < TEMPO_EPOCH_TIME_T)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "tio_time_utc_to_tempo: UTC time value is less than the TEMPO epoch");
        return -1;
     }

   lstt = Leap_Second_Table;
   if (lstt == NULL)
     {
        if (-1 == initialize ())
          return -1;
        lstt = Leap_Second_Table;
     }

   entries = lstt->entries;
   ilo = 0;
   ihi = lstt->num_entries;

   if ((utc_time < entries->unix_time)
       || (utc_time >= entries[ihi-1].unix_time))
     {
        if (utc_time < entries->unix_time)
          tai = utc_time;
        else
          tai = utc_time + entries[ihi-1].tai_utc_offset;

        if (utc_time >= lstt->expiration_time)
          tell_vwarn (0, "time value out of range, consider updating the leapsecond table");

        *tempo_time = tai - TEMPO_EPOCH_TAI;
        return 0;
     }

   while (ilo + 1 < ihi)
     {
        unsigned int i;

        i = ilo + (ihi - ilo)/2;
        if (utc_time >= entries[i].unix_time)
          ilo = i;
        else
          ihi = i;
     }

   tai = utc_time + entries[ilo].tai_utc_offset;
   *tempo_time = tai - TEMPO_EPOCH_TAI;

   return 0;
}

