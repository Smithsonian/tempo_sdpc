/** @file smc.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Dec 2021
 *  @brief Process maneuver files
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

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo)
{
   (void) pmt; (void) tpinfo;
   return 0;
}

static void delete_maneuver (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   FREE(pmt->out_dirname);
   FREE(pmt);
}

static int symlink_latest (const char *dir)
{
   IOCLib_Listdir_Type *ld = NULL;
   const char *latest = "latest";
   const char *newest_file;
#define BUFSIZE 128
   char buf[BUFSIZE];
   ssize_t buflen;
   int status = -1;

   /* Get a list of files sorted in ascending order */
   if (NULL == (ld = ioclib_listdir (dir, IOCLIB_LISTDIR_SORT)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: ioclib_listdir failed: %s",
                     __func__, dir);
        return -1;
     }
   if (ld->num_files == 0)
     {
        ioclib_listdir_free (ld);
        return 0;
     }

   /* We assume the files have a consistent naming scheme so that,
    * when sorted in ascending order, the newest file appears last.
    */
   newest_file = ld->files[ld->num_files-1];

   /* Read the 'latest' symlink */
   buflen = readlinkat (ld->dirfd, latest, buf, BUFSIZE);

   /* The 'latest' symlink should exist and should point to the newest file.
    * If we don't have that symlink, then try to create it atomically.
    * First, create the symlink under a temporary name, then rename it.
    */
   if ((buflen < 0) || (buflen == BUFSIZE)
       || (0 != strcmp (buf, newest_file)))
     {
        const char *tmp = "new_latest";
        if (symlinkat (newest_file, ld->dirfd, tmp) < 0)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: error creating symlink %s in %s: (%s)",
                          __func__, tmp, dir, strerror(errno));
             goto return_error;
          }
        if (renameat (ld->dirfd, tmp, ld->dirfd, latest) < 0)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: error renaming symlink %s in %s: (%s)",
                          __func__, tmp, dir, strerror(errno));
             goto return_error;
          }
     }

   status = 0;
return_error:
   ioclib_listdir_free (ld);
   return status;
}

static int process_maneuver (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                             const char *file, void *client_data)
{
   (void) tpinfo; (void) client_data;

   if (0 != copy_file_to_dir (file, pmt->out_dirname, ioclib_basename (file)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed copying %s to %s", __func__,
                     file, pmt->out_dirname);
        return -1;
     }

   return symlink_latest (pmt->out_dirname);
}

static int parse_maneuver_params (config_t *cfg, Process_Method_Type *pmt)
{
   const char *out_dirname;

   if (CONFIG_TRUE != config_lookup_string (cfg, "maneuver.output_dir", &out_dirname))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'maneuver' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (pmt->out_dirname = expand_string (out_dirname)))
     return -1;

   return 0;
}

Process_Method_Type *init_maneuver_method (config_t *cfg)
{
   Process_Method_Type *pmt = NULL;

   if (NULL == (pmt = (Process_Method_Type *)MALLOC (sizeof *pmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)pmt, 0, sizeof *pmt);

   if (-1 == parse_maneuver_params (cfg, pmt))
     {
        delete_maneuver (pmt);
        return NULL;
     }

   pmt->pmt_process = process_maneuver;
   pmt->pmt_delete = delete_maneuver;
   pmt->pmt_flush_cache = flush_cache;

   /* Unused */
   pmt->pmt_query_latest_timestamp = NULL;

   return pmt;
}
