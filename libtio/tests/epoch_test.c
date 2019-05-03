#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tio.h>
#include <math.h>

static int check_epoch (const char *epoch_str, time_t tt)
{
   double tempo0, utc0;

   if (0 != tio_time_set_taix_epoch (epoch_str))
     {
        fprintf (stderr, "*** FAIL setting epoch = %s\n", epoch_str);
        return -1;
     }

   tempo0 = 0.0;
   if ((0 != tio_time_taix_to_utc (tempo0, &utc0))
       || (utc0 != (double) tt))
     {
        fprintf (stderr, "*** FAIL converting tempo time to UTC\n");
        return -1;
     }

   if ((0 != tio_time_utc_to_taix (utc0, &tempo0))
       || (tempo0 != 0.0))
     {
        fprintf (stderr, "*** FAIL converting UTC to tempo time\n");
        return -1;
     }

   return 0;
}

int main (void)
{
   time_t noon_epoch_timet = 946728000L;
   time_t midnight_epoch_timet = noon_epoch_timet - 12*3600L;

   if (check_epoch ("2000-01-01T12:00:00Z", noon_epoch_timet))
     exit(EXIT_FAILURE);

   if (check_epoch ("2000-01-01T00:00:00Z", midnight_epoch_timet))
     exit(EXIT_FAILURE);

   return EXIT_SUCCESS;
}
