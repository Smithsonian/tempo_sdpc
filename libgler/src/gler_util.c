#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <wordexp.h>

#include <ioclib.h>
#include <tell.h>

#include <tio.h>
#include "gler_util.h"

#define BUFSIZE 1024

typedef struct
{
   char **files;
   int num_files;
   int *yearday;
}
GLER_File_Map_Type;

struct GLER_Type
{
   GLER_File_Map_Type *land;
   GLER_File_Map_Type *snow;
   GLER_File_Map_Type *ocean;
};

static void free_gler_file_map_type (GLER_File_Map_Type *fmt)
{
   if (NULL == fmt)
     return;
   if (fmt->files)
     {
        int i;
        for (i = 0; i < fmt->num_files; i++)
          {
             free (fmt->files[i]);
             fmt->files[i] = NULL;
          }
        free (fmt->files);
        fmt->files = NULL;
     }
   free(fmt->yearday);
   free(fmt);
}

static GLER_File_Map_Type *alloc_gler_file_map_type (int num_files)
{
   GLER_File_Map_Type *fmt = NULL;

   if (NULL == (fmt = (GLER_File_Map_Type *)malloc (sizeof *fmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)fmt, 0, sizeof (*fmt));

   if ((NULL == (fmt->yearday = (int *)malloc (num_files * sizeof(int))))
       || (NULL == (fmt->files = (char **)malloc (num_files * sizeof(char *)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_gler_file_map_type (fmt);
        return NULL;
     }
   memset ((char *)fmt->files, 0, num_files * sizeof (char *));
   fmt->num_files = num_files;

   return fmt;
}

static char *expand_string (const char *s, int quiet)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        if (!quiet) tell_verror (TELL_UNKNOWN_ERROR,
                                 "%s: expanding path: %s", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
     }

   return s_exp;
}

static char *expand_file_map_path (const char *path, int iwave)
{
   char *dirpath = NULL;
   char *exp_dirpath = NULL;
   char *exp_path = NULL;
   char *filled_path = NULL;
   size_t len;
   int n, status = -1;

   len = 9 + strlen(path);
   if (NULL == (filled_path = ioclib_malloc (len)))
     return NULL;

   n = snprintf (filled_path, len, path, iwave);
   if ((n < 0) || (n >= (int) len))
     goto return_status;

   if (NULL == (dirpath = ioclib_dirname (filled_path)))
     goto return_status;

   if (NULL == (exp_dirpath = expand_string (dirpath, 0)))
     goto return_status;

   if (NULL == (exp_path = ioclib_pathconcat (exp_dirpath, ioclib_basename(filled_path))))
     goto return_status;

   status = 0;
return_status:
   ioclib_free (dirpath);
   ioclib_free (filled_path);
   free (exp_dirpath);
   if (status)
     {
        ioclib_free (exp_path);
        exp_path = NULL;
     }

   return exp_path;
}

static GLER_File_Map_Type *gler_glob_file_map (const char *pattern, int iwave)
{
   GLER_File_Map_Type *fmt = NULL;
   IOCLib_Glob_Type *g = NULL;
   char *exp_pattern = NULL;
   size_t i;
   int status = -1;

   if (pattern == NULL)
     return NULL;

   if (NULL == (exp_pattern = expand_file_map_path (pattern, iwave)))
     return NULL;

   /* Get sorted file list */
   if (NULL == (g = ioclib_glob (exp_pattern, 0)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: glob failed", __func__);
        goto return_error;
     }

   if (g->num_files == 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: no files match glob pattern: %s", __func__, pattern);
        goto return_error;
     }

   if (NULL == (fmt = alloc_gler_file_map_type (g->num_files)))
     goto return_error;

   for (i = 0; i < g->num_files; i++)
     {
        int ncid, start=0, count=1;

        if (0 != TIO_open (g->files[i], NC_NOWRITE, &ncid))
          goto return_error;

        if (0 != TIO_get_var_section (ncid, "doy", &start, &count, NC_INT, &fmt->yearday[i]))
          goto return_error;

        TIO_close (ncid);

        if (NULL == (fmt->files[i] = strdup (g->files[i])))
          goto return_error;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: got yday=%d file=%s", __func__, fmt->yearday[i], fmt->files[i]);
     }

   status = 0;
return_error:
   ioclib_glob_free (g);
   free(exp_pattern);
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed", __func__);
        free_gler_file_map_type (fmt);
        fmt = NULL;
     }

   return fmt;
}

void gler_close (GLER_Type *gt)
{
   if (NULL == gt)
     return;
   free_gler_file_map_type (gt->land);
   free_gler_file_map_type (gt->snow);
   free_gler_file_map_type (gt->ocean);
   free(gt);
}

GLER_Type *gler_open (int iwave, const char *config_file)
{
   GLER_Type *gt = NULL;
   IOCLib_String_Array_Obj_Type *sa = NULL;
   char *land_glob = NULL;
   char *snow_glob = NULL;
   char *ocean_glob = NULL;
   char *cfg_file = NULL;
   int status = -1;
   IOCLib_Config_String_Type tbl[] =
     {
        {"GLER_Land_File_Pattern", &land_glob, IOCLIB_CONFIG_TYPE_STR},
        {"GLER_Snow_File_Pattern", &snow_glob, IOCLIB_CONFIG_TYPE_STR},
        {"GLER_Ocean_File_Pattern", &ocean_glob, IOCLIB_CONFIG_TYPE_STR},
        {NULL,NULL,0}
     };

   /* When a config file path is not provided, search
    * for the path in a few likely places, starting with
    * the closest, and expanding outward.
    */

   if ((config_file == NULL)
       && (NULL == (config_file = getenv ("SDPC_GLER_CONFIG"))))  /* check local environment customization */
     {
        char *conf_file_paths[] =
          {
             "$SDPC_PIPE_DIR/etc/table_config.ini",  /* running pipeline context */
             "$SDPC_ROOT/etc/table_config.ini",            /* software installation context */
             NULL
          };
        char **cfp;
         for (cfp = conf_file_paths; *cfp != NULL; cfp++)
          {
             if (NULL != (cfg_file = expand_string (*cfp, 1)))
               {
                  if (0 == access (cfg_file, F_OK | R_OK))
                    break;
               }
          }
        if (cfg_file == NULL)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: config_file not specified, default file not found", __func__);
             goto return_error;
           }
        config_file = cfg_file;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: reading config file: %s", __func__, config_file);

   if (NULL == (sa = ioclib_config_get_strings (config_file, "GLER", tbl)))
     goto return_error;

   if (NULL == (gt = (GLER_Type *)malloc (sizeof *gt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc error", __func__);
        goto return_error;
     }
   memset ((char *)gt, 0, sizeof (*gt));

   if (NULL == (gt->land = gler_glob_file_map (land_glob, iwave)))
     goto return_error;

   if (NULL == (gt->snow = gler_glob_file_map (snow_glob, iwave)))
     goto return_error;

   if (NULL == (gt->ocean = gler_glob_file_map (ocean_glob, iwave)))
     goto return_error;

   status = 0;
return_error:
   ioclib_string_array_obj_free (sa);
   if (cfg_file) free(cfg_file);
   if (status)
     {
        gler_close (gt);
        gt = NULL;
     }
   return gt;
}

static int bsearch_i (int t, const int *x, int n)
{
   int n0, n1, n2;
   int xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

/* Given a time, t, return two indices, (a,b), and a weight, awt,
 * such that the interpolation between files should be:
 *   f(t) = awt * file[a] + (1 - awt) * file[b]
 */
static int file_lookup (const GLER_File_Map_Type *fmt, double taix, int *a, int *b, double *awt)
{
   int year, yday, n, ka, kb;
   double yb, ya;

   /* This returns yday in the range [0,365] */
   if (0 != tio_time_taix_to_yearday (taix, &year, &yday))
     return -1;

   /* number days from [1,365] and ignore leap days */
   if (yday < 365) yday++;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: yday=%d", __func__, yday);

   n = fmt->num_files;

   /* We may lack sample points near the beginning/end of the year,
    * but the table should be viewed as periodic, so that's ok.
    */
   if ((yday < fmt->yearday[0])
       || (yday >= fmt->yearday[n-1]))
     {
        ka = n-1;
        kb = 0;
        ya = (double) fmt->yearday[ka];
        yb = 365.0 + fmt->yearday[kb];
        if (yday < fmt->yearday[0]) yday += 365;
     }
   else
     {
        ka = bsearch_i (yday, fmt->yearday, fmt->num_files);
        kb = ka + 1;
        ya = (double) fmt->yearday[ka];
        yb = (double) fmt->yearday[kb];
     }

   *awt = (yb - yday) / (yb - ya);

   *a = ka;
   *b = kb;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: a=%d, b=%d, awt=%f", __func__, ka, kb, *awt);

   return 0;
}

static int copy_path (const GLER_File_Map_Type *fmt, int k, char *path, int pathlen)
{
   char **files = fmt->files;
   int num_files = fmt->num_files;
   int len;

   if (k > num_files)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid index: %d", __func__, k);
        return -1;
     }

   len = strlen(files[k]);
   if (len >= pathlen)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: string length %d exceeds target buffer size %d",
                     __func__, len, pathlen);
        return -1;
     }

   (void) strncpy (path, files[k], pathlen);

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: path -> %s", __func__, path);

   return 0;
}

int gler_land_lookup (const GLER_Type *gt, double taix, int *a, int *b, double *awt)
{
   if (gt == NULL)
     return -1;
   return file_lookup (gt->land, taix, a, b, awt);
}

int gler_snow_lookup (const GLER_Type *gt, double taix, int *a, int *b, double *awt)
{
   if (gt == NULL)
     return -1;
   return file_lookup (gt->snow, taix, a, b, awt);
}

int gler_ocean_lookup (const GLER_Type *gt, double taix, int *a, int *b, double *awt)
{
   if (gt == NULL)
     return -1;
   return file_lookup (gt->ocean, taix, a, b, awt);
}

int gler_land_file (const GLER_Type *gt, int k, char *path, int pathlen)
{
   if (gt == NULL)
     return -1;
   return copy_path (gt->land, k, path, pathlen);
}

int gler_snow_file (const GLER_Type *gt, int k, char *path, int pathlen)
{
   if (gt == NULL)
     return -1;
   return copy_path (gt->snow, k, path, pathlen);
}

int gler_ocean_file (const GLER_Type *gt, int k, char *path, int pathlen)
{
   if (gt == NULL)
     return -1;
   return copy_path (gt->ocean, k, path, pathlen);
}

