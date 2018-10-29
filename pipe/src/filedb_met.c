/** @file filedb_met.c
 *  @brief Weather forecast database
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
   char *p;

   /* e.g. YYYYMMDDHH.nam.tffz.conusnest.hiresfHH.tm00.grib2 */

   p = strptime (basename, "%Y%m%d%H", tm);
   len = p - basename;
   if (len == 10)
     return 0;

   fprintf (stderr, "*** Error: parsing %s\n", basename);
   return -1;
}

int config_met (Filedb_Type *fdb, config_t *cfg, const char *name)
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
