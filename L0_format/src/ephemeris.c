/** @file smc.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Dec 2021
 *  @brief Process ephemeris files
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tio.h>
#include <tio_template.h>
#include <tell.h>

#define PROCESS_METHOD_PRIVATE_DATA \
   char *out_dirname; \
   double outfile_timestamp_end;
#include "l0_format.h"

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        int unwind, const char *incoming_dir)
{
   (void) pmt; (void) tpinfo; (void) unwind; (void) incoming_dir;
   return 0;
}

static void delete_ephemeris (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   FREE(pmt->out_dirname);
   FREE(pmt);
}

static int process_ephemeris (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                             const char *file, void *client_data)
{
   (void) tpinfo; (void) client_data;

   if (0 != copy_file_to_dir (file, pmt->out_dirname, ioclib_basename (file)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed copying %s to %s", __func__,
                     file, pmt->out_dirname);
          return -1;
     }
   return 0;
}

static int parse_ephemeris_params (config_t *cfg, Process_Method_Type *pmt)
{
   const char *out_dirname;

   if (CONFIG_TRUE != config_lookup_string (cfg, "ephemeris.output_dir", &out_dirname))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'ephemeris' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (pmt->out_dirname = expand_string (out_dirname)))
     return -1;

   return 0;
}

Process_Method_Type *init_ephemeris_method (config_t *cfg)
{
   Process_Method_Type *pmt = NULL;

   if (NULL == (pmt = (Process_Method_Type *)MALLOC (sizeof *pmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)pmt, 0, sizeof *pmt);

   if (-1 == parse_ephemeris_params (cfg, pmt))
     {
        delete_ephemeris (pmt);
        return NULL;
     }

   pmt->pmt_process = process_ephemeris;
   pmt->pmt_delete = delete_ephemeris;
   pmt->pmt_flush_cache = flush_cache;

   /* Unused */
   pmt->pmt_query_latest_timestamp = NULL;

   return pmt;
}
