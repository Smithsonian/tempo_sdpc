#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <tell.h>
#include <ioclib.h>

#include "l0_format.h"

/* This module maintains a a disk cache of exposure record files,
 * each containing one exposure record.
 * The calling routine ensures that all records in the disk
 * cache are of the same type and could occupy the same granule.
 * For example, all cached radiances belong to a single scan.
 * Cached dark or irradiance exposure records all belong to a
 * single observing "session".
 * The cache is maintained as a FIFO queue.  Records are added
 * at the end and removed from the front.
 */

#define EXPREC_CACHE_METHOD_PRIVATE_DATA \
   char *cache_dirname; \
   size_t num_erecs_cached; \
   char **files; \
   size_t num_files; \
   size_t front;
#include "exprec_cache.h"

static int cache_erec (Exprec_Cache_Method_Type *cmt, const char *file, size_t file_index)
{
   const char *basename = ioclib_basename (file);
   char *cache_path;
   int status;

   /* For this method, we assume all files have 1 erec, so file_index==0, always */
   (void) file_index;

   if (NULL == (cache_path = ioclib_pathconcat (cmt->cache_dirname, basename)))
     return -1;

   status = ioclib_rename (file, cache_path);
   ioclib_free (cache_path);

   if (status == 0)
     {
        cmt->num_erecs_cached++;
     }

   return status;
}

static int cache_num_recs (Exprec_Cache_Method_Type *cmt, size_t *num_erecs_cached)
{
   /* num_files omits any known-bad files (in the .bad subdirectory) */
   if (cmt->num_files > 0)
     *num_erecs_cached = cmt->num_files;
   else
     *num_erecs_cached = cmt->num_erecs_cached;

   return 0;
}

static int cache_close (Exprec_Cache_Method_Type *cmt)
{
   ioclib_string_array_free (cmt->files, cmt->num_files);
   cmt->files = NULL;
   cmt->num_files = 0;
   cmt->front = 0;
   return 0;
}

static int cache_open (Exprec_Cache_Method_Type *cmt)
{
   cmt->files = ioclib_dir_list(cmt->cache_dirname, &cmt->num_files, IOCLIB_LISTDIR_SORT);
   cmt->front = 0;

   return cmt->files ? 0 : -1;
}

/* return the path at the front of the FIFO queue */
static char *get_oldest_file_path (Exprec_Cache_Method_Type *cmt)
{
   if (cmt->files == NULL)
     return NULL;

   if (cmt->front >= cmt->num_files)
     return NULL;

   return ioclib_pathconcat (cmt->cache_dirname, cmt->files[cmt->front]);
}

static IOCSDPC_Exprec_Type *cache_erec_get (Exprec_Cache_Method_Type *cmt)
{
   IOCSDPC_Common_Header_Type chdr = {0};
   IOCSDPC_Exprec_Type *erec = NULL;
   char *path = NULL;
   int fd;

   if (NULL == (path = get_oldest_file_path (cmt)))
     return NULL;

   if (-1 == (fd = iocsdpc_open_file_read (path, 0, &chdr)))
     {
        ioclib_free (path);
        return NULL;
     }

   if (NULL == (erec = iocsdpc_exprec_fdopen_read (path, fd, &chdr)))
     {
        ioclib_fd_close (fd);
        erec = NULL;
     }

   ioclib_free (path);

   return erec;
}

static int cache_erec_bad (Exprec_Cache_Method_Type *cmt)
{
   char *path = NULL;
   int status;

   if (NULL == (path = get_oldest_file_path (cmt)))
     return -1;

   cmt->front++;

   if (0 != (status = ioclib_rename_to_bad_file (path)))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: ioclib_rename_to_bad_file, failed: file=%s",
                     __func__, path);
     }
   ioclib_free (path);

   return status;
}

static int cache_erec_done (Exprec_Cache_Method_Type *cmt)
{
   char *path = NULL;
   int status;

   if (NULL == (path = get_oldest_file_path (cmt)))
     return -1;

   cmt->front++;

   if (0 != (status = ioclib_unlink (path)))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: unlink failed: file=%s", __func__, path);
     }
   ioclib_free (path);
   cmt->num_erecs_cached--;

   return status;
}

static int cache_erec_path (Exprec_Cache_Method_Type *cmt, char *buf, size_t buflen)
{
   char *path = NULL;
   size_t len;

   if (NULL == (path = get_oldest_file_path (cmt)))
     return -1;

   strncpy (buf, path, buflen);
   len = strlen (path);
   if (len >= buflen) buf[buflen-1] = 0;

   ioclib_free (path);

   return len;
}

static void cache_delete (Exprec_Cache_Method_Type *cmt)
{
   FREE(cmt->cache_dirname);
   ioclib_string_array_free (cmt->files, cmt->num_files);
   FREE(cmt);
}

static int valid_cache_directory_path (const char *path)
{
   struct stat st;

   /* nonexistent is ok - we'll create it */
   if (0 != stat (path, &st))
     return 1;

   /* an empty, accessible directory is ok */
   if (S_ISDIR(st.st_mode))
     {
        size_t num_files = 0;
        char **files = NULL;
        if (NULL == (files = ioclib_dir_list (path, &num_files, IOCLIB_LISTDIR_SORT)))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: ioclib_dir_list failed: %s", __func__, path);
             return 0;  /* inaccessible is not ok */
          }
        ioclib_string_array_free (files, num_files);
        if (num_files == 0) return 1;
        tell_verror (TELL_RUNTIME_ERROR, "%s: cache directory is not empty: %s", __func__, path);
     }

   /* 0 means invalid */
   return 0;
}

static int parse_exprec_cache_params (config_t *cfg, Exprec_Cache_Method_Type *cmt)
{
   config_setting_t *s;
   const char *cache_dirname;

   if (NULL == (s = config_lookup (cfg, "exprec")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'exprec' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "cache_dir", &cache_dirname))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'exprec' cache parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (cmt->cache_dirname = expand_string (cache_dirname)))
     return -1;

   if (0 == valid_cache_directory_path (cmt->cache_dirname))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid cache directory path: %s",
                     __func__, cmt->cache_dirname);
        return -1;
     }

   if (0 != ioclib_mkdir (cmt->cache_dirname, 0))
     return -1;

   return 0;
}

Exprec_Cache_Method_Type *open_erec_cache_disk (config_t *cfg)
{
   Exprec_Cache_Method_Type *cmt = NULL;

   if (NULL == (cmt = (Exprec_Cache_Method_Type *)MALLOC (sizeof *cmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)cmt, 0, sizeof *cmt);

   if (0 != parse_exprec_cache_params (cfg, cmt))
     {
        cache_delete (cmt);
        return NULL;
     }

   cmt->cache_erec = cache_erec;
   cmt->cache_num_recs = cache_num_recs;
   cmt->cache_open = cache_open;
   cmt->cache_close = cache_close;
   cmt->cache_erec_get = cache_erec_get;
   cmt->cache_erec_path = cache_erec_path;
   cmt->cache_erec_bad = cache_erec_bad;
   cmt->cache_erec_done = cache_erec_done;
   cmt->cache_delete = cache_delete;

   return cmt;
}
