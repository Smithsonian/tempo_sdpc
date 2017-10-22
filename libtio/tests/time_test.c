#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tio.h>
#include <math.h>

#define BUFSIZE 32

static int check_time_conversions (double tempo)
{
   char buf[BUFSIZE];
   double other, tempo_from_utc, hour, minf, sec;
   time_t tt;
   struct tm tm;
   int year, month, day, hr, min;

  if (-1 == tio_time_tempo_to_utc (tempo, &other))
     return -1;
   if (-1 == tio_time_tempo_to_utc_caldate (tempo, &year, &month, &day, &hour))
     return -1;
   if (0 != tio_time_utc_to_tempo (other, &tempo_from_utc))
     return -1;

   tt = (time_t)other;
   gmtime_r (&tt, &tm);
   strftime (buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
   fprintf (stdout, "%s.%06d UTC : %06d usec\n", buf,
            (int)(round((other-tt)*1e6)),
            (int)(round((tempo-tempo_from_utc)*1e6)));

   hr   = (int)hour;
   minf = (hour - hr)*60;
   min = (int)minf;
   sec = (minf - min)*60;

   /* this string should exactly match the above UTC string */
   fprintf (stdout, "%4d-%02d-%02dT%02d:%02d:%09.6f\n",
            year, month, day, hr, min, sec);

   if (-1 == tio_time_tempo_to_tai (tempo, &other))
     return -1;
   tt = (time_t)other;
   gmtime_r (&tt, &tm);
   strftime (buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
   fprintf (stdout, "%s.%06d TAI\n", buf, (int)(round((other-tt)*1e6)));

   return 0;
}

int main (void)
{
   double tempo_time;
   int status;

   while (!feof(stdin))
     {
        status = fscanf (stdin, "%lf", &tempo_time);
        if (status != 1)
          {
             if (status == EOF) break;
             fprintf (stderr, "*** FAIL (corrupted input?)\n");
             exit(EXIT_FAILURE);
          }
        status = check_time_conversions (tempo_time);
        if (status)
          {
             fprintf (stderr, "*** FAIL (time scale conversion test)\n");
             exit(EXIT_FAILURE);
          }
     }

   return EXIT_SUCCESS;
}
