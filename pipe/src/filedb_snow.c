/** @file filedb_snow.c
 *  @brief Snow and ice cover database
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
   if (len < 22)
     {
        fprintf (stderr, "**** Error: basename appears too short: %s (strlen = %ld)\n",
                 basename, len);
        return -1;
     }

   /* e.g. ims2020182_1km_GIS_v1.3.tif */
   tstamp_start = basename + 3;
   p = strptime (tstamp_start, "%Y%j", tm);
   len = p - tstamp_start;
   if (len == 7)
     return 0;

   fprintf (stderr, "*** Error: parsing %s\n", basename);
   return -1;
}

int config_snow (Filedb_Type *fdb, config_t *cfg, const char *name)
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
