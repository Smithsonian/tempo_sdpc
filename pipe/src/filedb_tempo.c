/** @file filedb_tempo.c
 *  @brief TEMPO data products
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
   char *p;

   len = strlen(basename);
   if (len != 36)
     {
        fprintf (stderr, "**** Error: basename appears too short: %s (strlen = %ld)\n",
                 basename, len);
        return -1;
     }

   /* e.g. TEMPO_irr_L1_V01_20181028T013009Z.nc */
   tstamp_start = basename + 17;
   p = strptime (tstamp_start, "%Y%m%dT%H%M%S", tm);
   len = p - tstamp_start;
   if (len == 15)
     return 0;

   fprintf (stderr, "*** Error: parsing %s\n", basename);
   return -1;
}

int config_tempo (Filedb_Type *fdb, config_t *cfg, const char *name)
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
