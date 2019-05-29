/** @file filedb_ephemeris.c
 *  @brief ephemeris
 */

#include "config.h"
#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>

#include "filedb.h"

static int parse_timestamp (const char *basename, struct tm *tm)
{
   size_t len;
   const char *tstamp_start;
   const char *p;

   len = strlen(basename);
   if (len < 34)
     {
        fprintf (stderr, "*** Error: basename appears too short: %s (strlen = %ld)\n",
                 basename, len);
        return -1;
     }

   /* e.g. tempo_20190522T225358Z_lt_pred.eph */
   if (NULL == (tstamp_start = strchr (basename, '2')))
     goto parse_failed;

   if ((tstamp_start + 15 > basename + len)
       || (tstamp_start[8] != 'T')
       || (tstamp_start[15] != 'Z'))
     goto parse_failed;

   p = strptime (tstamp_start, "%Y%m%dT%H%M%S", tm);
   len = p - tstamp_start;
   if (len == 15)
     return 0;

parse_failed:
   fprintf (stderr, "*** Error: parsing %s\n", basename);
   return -1;
}

int config_ephemeris (Filedb_Type *fdb, config_t *cfg, const char *name)
{
   const char *cfg_name = name;
   const char *sc;

   /* name may have an embedded selector string */
   if (NULL != (sc = strchr (cfg_name, ':')))
     {
        cfg_name = sc + 1;
     }

   if (0 != read_config_common (fdb, cfg, cfg_name))
     return -1;

   fdb->parse_timestamp = parse_timestamp;

   return 0;
}
