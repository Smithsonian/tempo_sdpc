/** @file filedb.c
 *  @brief Main program
 */

#include "config.h"
#define _XOPEN_SOURCE 500  /* for nftw */
#define _GNU_SOURCE        /* for timegm */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <getopt.h>
#include <wordexp.h>
#include <ftw.h>
#include <fnmatch.h>

#include <libconfig.h>

#include <ioclib.h>
#include <tio.h>
#include "filedb.h"

#define FILEDB_MODE (S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP)
#define FILEDB_LIST_SUFFIX   ".lst"

#define MAX_NUM_OPEN_DIRS 15

struct Filedb_Entry_Type
{
   Filedb_Entry_Type *next;
   char *fpath;
   double timestamp;
   long offset;
};

typedef struct
{
   double *timestamp;
   size_t *order;
   long *offset;
   size_t *size;
   size_t num_entries;
}
Filedb_Table_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: filedb DBNAME [options]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -c | --config FILE     configuration file\n");
   fprintf (stderr, "   -u | --update          update the lookup table\n");
   fprintf (stderr, "   -f | --find            search the lookup table\n");
   fprintf (stderr, "   -H | --header FILE     match header timestamp in TEMPO data product FILE\n");
   fprintf (stderr, "   -s | --sec SECONDS     UTC time elapsed since the Unix epoch [sec]\n");
   fprintf (stderr, "                          (e.g. a unix time_t value)\n");
   fprintf (stderr, "   -d | --delay SECONDS   delay closing database file (for testing only)\n\n");
   fprintf (stderr, "   -h | --help            print this usage message\n");
   fprintf (stderr, "WARNING: Because locking of network-mounted files is unreliable,\n");
   fprintf (stderr, "         lookup tables should reside on a local disk.\n");
   exit (EXIT_SUCCESS);
}

static int read_config_file (const char *config_file, config_t *cfg)
{
   if (0 == config_read_file (cfg, config_file))
     {
        fprintf (stderr, "%s: Reading %s:%d - %s\n",
                 __func__, config_error_file(cfg),
                 config_error_line(cfg), config_error_text(cfg));
        return -1;
     }

   return 0;
}

static char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        fprintf (stderr, "%s: expanding path: %s\n", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        fprintf (stderr, "%s: strdup failed\n", __func__);
     }

   return s_exp;
}

static void filedb_entry_free (Filedb_Entry_Type *et)
{
   if (et == NULL)
     return;
   FREE(et->fpath);
   FREE(et);
}

static Filedb_Entry_Type *filedb_entry_new (const char *fpath, double timestamp)
{
   Filedb_Entry_Type *et = NULL;

   if (NULL == (et = (Filedb_Entry_Type *)MALLOC (sizeof *et)))
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        return NULL;
     }
   memset ((char *)et, 0, sizeof(*et));

   if (NULL == (et->fpath = strdup (fpath)))
     {
        fprintf (stderr, "%s: strdup failed\n", __func__);
        filedb_entry_free (et);
        return NULL;
     }

   et->next = NULL;
   et->offset = 0;
   et->timestamp = timestamp;

   return et;
}

void filedb_list_free (Filedb_Entry_Type *lst)
{
   while (lst)
     {
        Filedb_Entry_Type *next = lst->next;
        filedb_entry_free (lst);
        lst = next;
     }
}

static void filedb_table_free (Filedb_Table_Type *ft)
{
   if (ft == NULL)
     return;
   FREE(ft->timestamp);
   FREE(ft->order);
   FREE(ft->size);
   FREE(ft->offset);
   FREE(ft);
}

static Filedb_Table_Type *filedb_table_alloc (size_t num_entries)
{
   Filedb_Table_Type *ft = NULL;

   if (NULL == (ft = (Filedb_Table_Type *)MALLOC (sizeof *ft)))
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        return NULL;
     }
   memset ((char *)ft, 0, sizeof *ft);

   if ((NULL == (ft->timestamp = (double *)MALLOC (num_entries * sizeof(*ft->timestamp))))
       || (NULL == (ft->order = (size_t *)MALLOC (num_entries * sizeof(*ft->order))))
       || (NULL == (ft->size = (size_t *)MALLOC (num_entries * sizeof(*ft->size))))
       || (NULL == (ft->offset = (long *)MALLOC (num_entries * sizeof(*ft->offset))))
      )
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        filedb_table_free (ft);
        return NULL;
     }

   ft->num_entries = num_entries;

   return ft;
}

static double *_pTimestamp;

static int compare_timestamps (const void *va, const void *vb)
{
   size_t ia = *(const size_t *)va;
   size_t ib = *(const size_t *)vb;
   double ta = _pTimestamp[ia];
   double tb = _pTimestamp[ib];
   if (ta < tb) return -1;
   else if (ta > tb) return +1;
   else return 0;
}

static int filedb_table_sort (Filedb_Table_Type *ft)
{
   _pTimestamp = ft->timestamp;
   qsort (ft->order, ft->num_entries, sizeof(size_t), compare_timestamps);
   return 0;
}

static int filedb_table_write (int fd, const Filedb_Table_Type *ft)
{
   size_t i, o, n = ft->num_entries;
   ssize_t num_bytes_written;

   num_bytes_written = write (fd, &ft->num_entries, sizeof(size_t));
   if (sizeof(size_t) != num_bytes_written)
     return -1;

   for (i = 0; i < n; i++)
     {
        o = ft->order[i];
        num_bytes_written = write (fd, &ft->timestamp[o], sizeof(double));
        if (sizeof(double) != num_bytes_written)
          return -1;
     }

   for (i = 0; i < n; i++)
     {
        o = ft->order[i];
        num_bytes_written = write (fd, &ft->offset[o], sizeof(long));
        if (sizeof(long) != num_bytes_written)
          return -1;
     }

   for (i = 0; i < n; i++)
     {
        o = ft->order[i];
        num_bytes_written = write (fd, &ft->size[o], sizeof(size_t));
        if (sizeof(size_t) != num_bytes_written)
          return -1;
     }

   return 0;
}

static Filedb_Table_Type *filedb_table_read (int fd)
{
   Filedb_Table_Type *ft = NULL;
   size_t num_entries, len;
   ssize_t num_bytes_read;

   len = sizeof(size_t);
   if (((num_bytes_read = read (fd, &num_entries, len)) < 0)
       || ((size_t) num_bytes_read != len))
     goto error_return;

   if (NULL == (ft = filedb_table_alloc (num_entries)))
     goto error_return;

   len = num_entries * sizeof(double);
   if (((num_bytes_read = read (fd, ft->timestamp, len)) < 0)
       || ((size_t) num_bytes_read != len))
     goto error_return;

   len = num_entries * sizeof(long);
   if (((num_bytes_read = read (fd, ft->offset, len)) < 0)
       || ((size_t) num_bytes_read != len))
     goto error_return;

   len = num_entries * sizeof(size_t);
   if (((num_bytes_read = read (fd, ft->size, len)) < 0)
       || ((size_t) num_bytes_read != len))
     goto error_return;

   /* ft->order is used only on output.  It's not used for table lookups */
   memset ((char *)ft->order, 0, num_entries * sizeof(size_t));

   return ft;
error_return:
   filedb_table_free (ft);
   return NULL;
}

static int open_with_lock (const char *filename, int readonly)
{
   struct stat st;
   mode_t mode = FILEDB_MODE;
   int lock_type;
   int fd;

   if (0 != stat (filename, &st))
     {
        if (readonly)
          {
             fprintf (stderr, "*** %s: nonexistent file %s (%s)\n",
                        __func__, filename, strerror(errno));
             return -1;
          }
        else if ((fd = creat (filename, mode)) < 0)
          {
             fprintf (stderr, "*** %s: error creating %s (%s)\n",
                        __func__, filename, strerror(errno));
             return -1;
          }
     }
   else
     {
        int flags = readonly ? O_RDONLY : O_RDWR;

        if ((fd = open (filename, flags, mode)) < 0)
          {
             fprintf (stderr, "*** %s: error opening %s (%s)\n",
                      __func__, filename, strerror(errno));
             return -1;
          }
     }

   lock_type = readonly ? LOCK_SH : LOCK_EX;

   /* may block if another process holds the lock */
   if (0 != flock (fd, lock_type))
     {
        fprintf (stderr, "*** %s: flock failed on %s (%s)\n",
                 __func__, filename, strerror(errno));
        (void) close (fd);
        return -1;
     }

   return fd;
}

static double Delay_Time;

static void delay_if_requested (void)
{
   struct timespec req, rem;
   time_t delay_sec;

   if (Delay_Time <= 0)
     return;

   delay_sec = (time_t) Delay_Time;

   req.tv_sec = delay_sec;
   req.tv_nsec = (long) (Delay_Time - delay_sec) * 1.e9;

   fprintf (stderr, "delay %g sec...\n", Delay_Time);

   if (0 == nanosleep (&req, &rem))
     return;

   fprintf (stderr, "delay interrupted: %g sec remaining\n", (rem.tv_sec + rem.tv_nsec/1.e9));
}

static int close_and_unlock (int fd)
{
   delay_if_requested ();
   (void) flock (fd, LOCK_UN);
   return close (fd);
}

static size_t list_size (const Filedb_Entry_Type *lst)
{
   size_t num_entries = 0;

   while (lst)
     {
        Filedb_Entry_Type *next = lst->next;
        num_entries++;
        lst = next;
     }

   return num_entries;
}

static Filedb_Table_Type *sorted_table_from_list (const Filedb_Entry_Type *lst)
{
   Filedb_Table_Type *ft = NULL;
   const Filedb_Entry_Type *et = NULL;
   size_t i, num_entries;

   if (0 == (num_entries = list_size (lst)))
     return NULL;

   if (NULL == (ft = filedb_table_alloc (num_entries)))
     return NULL;

   i = 0;

   for (et = lst; et != NULL; et = et->next)
     {
        ft->timestamp[i] = et->timestamp;
        ft->order[i] = i;
        ft->offset[i] = 0;
        ft->size[i] = strlen (et->fpath);
        i++;
     }

   if (0 != filedb_table_sort (ft))
     {
        filedb_table_free (ft);
        return NULL;
     }

   return ft;
}

static char *append_suffix (const char *s, const char *suffix)
{
   size_t len = strlen(s) + strlen(suffix) + 1;
   char *ns = NULL;

   if (NULL == (ns = (char *)MALLOC (len)))
     {
        fprintf (stderr, "*** %s: malloc failed\n", __func__);
        return NULL;
     }

   snprintf (ns, len, "%s%s", s, suffix);
   return ns;
}

static int write_list_and_sync_table (const char *filename, const Filedb_Entry_Type *lst,
                                     Filedb_Table_Type *ft)
{
   const Filedb_Entry_Type *et;
   char *lst_filename = NULL;
   FILE *fp = NULL;
   size_t i;
   int status = -1;

   if (NULL == (lst_filename = append_suffix (filename, FILEDB_LIST_SUFFIX)))
     return -1;

   if (NULL == (fp = fopen (lst_filename, "w")))
     {
        fprintf (stderr, "*** %s: error opening %s for writing\n", __func__, lst_filename);
        goto free_and_return;
     }

   for (et = lst, i = 0; et != NULL; et = et->next, i++)
     {
        ft->offset[i] = ftell (fp);
        if (fprintf (fp, "%s\n", et->fpath) < 0)
          {
             fprintf (stderr, "*** %s: error writing %s\n", __func__, lst_filename);
             goto free_and_return;
          }
     }

   status = 0;
free_and_return:
   if (0 != fclose (fp))
     {
        fprintf (stderr, "*** %s: error closing %s\n", __func__, lst_filename);
     }
   FREE(lst_filename);
   return status;
}

static int filedb_write_while_locked (int fd, const Filedb_Entry_Type *lst, const char *filename)
{
   Filedb_Table_Type *ft = NULL;
   int status = -1;

   if (NULL == (ft = sorted_table_from_list (lst)))
     goto free_and_return;

   if (0 != write_list_and_sync_table (filename, lst, ft))
     goto free_and_return;

   if (0 != filedb_table_write (fd, ft))
     goto free_and_return;

   status = 0;
free_and_return:
   filedb_table_free (ft);
   return status;
}

static int filedb_write (const char *filename, const Filedb_Entry_Type *lst)
{
   int fd, status;

   if (lst == NULL)
     {
        fprintf (stderr, "*** %s: empty filename list?\n", __func__);
        return -1;
     }

   if ((fd = open_with_lock (filename, 0)) < 0)
     return -1;

   status = filedb_write_while_locked (fd, lst, filename);

   if (0 != close_and_unlock (fd))
     {
        fprintf (stderr, "*** %s: failed writing %s\n", __func__, filename);
        status = -1;
        /* drop */
     }

   return status;
}

static size_t bsearch_d (double t, double *x, size_t n)
{
   size_t n0, n1, n2;
   double xt;

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

static char *read_filename_at_offset (const char *filename, long offset, size_t size)
{
   FILE *fp = NULL;
   char *lst_filename = NULL;
   char *db_file = NULL;

   if (NULL == (lst_filename = append_suffix (filename, FILEDB_LIST_SUFFIX)))
     return NULL;

   if (NULL == (fp = fopen (lst_filename, "r")))
     {
        fprintf (stderr, "*** %s: failed opening %s\n", __func__, lst_filename);
        goto return_status;
     }

   if (0 != fseek (fp, offset, SEEK_SET))
     {
        fprintf (stderr, "*** %s: illegal seek: %s\n", __func__, lst_filename);
        goto return_status;
     }

   /* leave space for the terminating null byte */
   if (NULL == (db_file = (char *)MALLOC (size + 1)))
     {
        fprintf (stderr, "*** %s: malloc failed\n", __func__);
        goto return_status;
     }

   if (NULL == fgets (db_file, size + 1, fp))
     {
        fprintf (stderr, "*** %s: failed reading from %s\n", __func__, lst_filename);
        goto return_status;
     }

return_status:
   if (fp) fclose (fp);
   FREE(lst_filename);

   return db_file;
}

static char *filedb_lookup (const char *filename, double query_time)
{
   Filedb_Table_Type *ft = NULL;
   char *db_file = NULL;
   size_t n, id;
   int fd, status = -1;

   if ((fd = open_with_lock (filename, 1)) < 0)
     return NULL;

   if (NULL == (ft = filedb_table_read (fd)))
     {
        fprintf (stderr, "*** %s: failed reading %s\n", __func__, filename);
        goto cleanup_and_return;
     }

   n = ft->num_entries;

   if (query_time < ft->timestamp[0])
     id = 0;
   else if (ft->timestamp[n-1] < query_time)
     id = n-1;
   else
     {
        id = bsearch_d (query_time, ft->timestamp, n);

        /* FIXME?
         * Maybe the lookup should be file-type specific?
         *
         * Ideally, each file would span a defined interval, and we would
         * pick the file that contained the query point.  But the files
         * aren't always convenient to read (grib2!), and an applicable
         * interval is not always defined (e.g. snow cover).
         *
         * For now, we pick the file that has the nearest timestamp.
         * For the ephemeris, we may have to pre-process the ephemeris data
         * files to remove any ambiguities.
         */

        if (id+1 < n)
          {
             double t0 = fabs(ft->timestamp[id] - query_time);
             double t1 = fabs(query_time - ft->timestamp[id+1]);
             if (t1 < t0) id++;
          }
     }

   if (NULL == (db_file = read_filename_at_offset (filename, ft->offset[id], ft->size[id])))
     goto cleanup_and_return;

   status = 0;
cleanup_and_return:
   filedb_table_free (ft);
   (void) close_and_unlock (fd);

   if (status)
     {
        FREE(db_file);
        db_file = NULL;
     }

   return db_file;
}

static void free_filedb_type (Filedb_Type *fdb)
{
   if (fdb == NULL)
     return;
   filedb_list_free (fdb->lst);
   FREE(fdb->root_dir);
   FREE(fdb->lookup_table);
   FREE(fdb);
}

static Filedb_Type *alloc_filedb_type (void)
{
   Filedb_Type *fdb = NULL;

   if (NULL == (fdb = (Filedb_Type *)MALLOC (sizeof *fdb)))
     {
        fprintf (stderr, "%s: malloc failed\n", __func__);
        return NULL;
     }

   memset ((char *)fdb, 0, sizeof *fdb);

   return fdb;
}

static int append_entry (Filedb_Type *fdb, const char *fpath, double timestamp)
{
   Filedb_Entry_Type *et = NULL;

   if (NULL == (et = filedb_entry_new (fpath, timestamp)))
     return -1;

   et->next = fdb->lst;
   fdb->lst = et;

   return 0;
}

static Filedb_Type *__pFDB;

static void set_global_filedb_ptr (Filedb_Type *fdb)
{
   __pFDB = fdb;
}

static Filedb_Type *get_global_filedb_ptr (void)
{
   return __pFDB;
}

static int direntry_handler (const char *fpath, const struct stat *sb, int typeflag,
                             struct FTW *ftwbuf)
{
   int fnm_flags = FNM_PATHNAME | FNM_PERIOD | FNM_NOESCAPE;
   Filedb_Type *fdb = get_global_filedb_ptr();
   const char *pbasename;
   struct tm tm = {0};
   double utc;

   (void) sb;

   /* ignore anything that isn't a regular file */
   if (typeflag != FTW_F)
     return 0;

   pbasename = fpath + ftwbuf->base;

   /* ignore filenames that don't match the expected pattern */
   if (0 != fnmatch (fdb->basename_pattern, pbasename, fnm_flags))
     return 0;

   if (0 != fdb->parse_timestamp (pbasename, &tm))
     return -1;

   /* timegm is a GNU extension, but this is easier than persuading
    * mktime to interpret the struct tm as UTC instead of local time.
    */
   utc = (double) timegm (&tm);

   return append_entry (fdb, fpath, utc);
}

static int create_lookup_table_dir (const char *path)
{
   char *dirname = NULL;
   int status = -1;

   if (NULL == (dirname = ioclib_dirname (path)))
     {
        fprintf (stderr, "*** %s: unable to extract directory path from %s\n",
                 __func__, path ? path : "(null)");
        return -1;
     }

   /* mode 07555 = drwxr-xr-x */
   if (0 != ioclib_mkdir (dirname, 0755))
     {
        fprintf (stderr, "*** %s: unable to create directory path %s\n",
                 __func__, dirname);
        goto cleanup_and_return;
     }

   status = 0;
cleanup_and_return:
   ioclib_free (dirname);
   return status;
}

static int fdb_initialize (Filedb_Type *fdb)
{
   int nopenfd = MAX_NUM_OPEN_DIRS;
   int nftw_flags = FTW_MOUNT | FTW_PHYS;

   if (0 != create_lookup_table_dir (fdb->lookup_table))
     return -1;

   set_global_filedb_ptr (fdb);

   if (0 != nftw (fdb->root_dir, direntry_handler, nopenfd, nftw_flags))
     {
        fprintf (stderr, "%s: nftw failed looking at root_dir=%s\n",
		 __func__, fdb->root_dir);
        return -1;
     }

   if (fdb->lst == NULL)
     {
        fprintf (stderr, "*** No matching files in directory %s\n",
                 fdb->root_dir);
        return -1;
     }

   if (0 != filedb_write (fdb->lookup_table, fdb->lst))
     return -1;

   return 0;
}

int read_config_common (Filedb_Type *fdb, config_t *cfg, const char *name)
{
   config_setting_t *s;
   const char *root_dir;
   const char *lookup_table;

   if (NULL == (s = config_lookup (cfg, name)))
     {
        fprintf (stderr, "%s: accessing config setting %s:%s\n",
                 __func__, name, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "root_dir", &root_dir))
       ||(CONFIG_TRUE != config_setting_lookup_string (s, "basename_pattern", &fdb->basename_pattern))
       ||(CONFIG_TRUE != config_setting_lookup_string (s, "lookup_table", &lookup_table)))
     {
        fprintf (stderr, "%s: reading config setting %s:%s\n",
                 __func__, name, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (fdb->root_dir = expand_string (root_dir)))
       || (NULL == (fdb->lookup_table = expand_string (lookup_table))))
     return -1;

   return 0;
}

static int query_file_timestamp (const char *file, double *timestamp_utc)
{
   int ncid, status = -1;
   double tempo_time;

   if (0 != access (file, F_OK | R_OK))
     {
        fprintf (stderr, "*** %s: cannot read file %s\n", __func__, file);
        return -1;
     }

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if (0 != tio_use_file_epoch (ncid))
     goto close_and_return;

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &tempo_time))
     goto close_and_return;

   if (0 != tio_time_taix_to_utc (tempo_time, timestamp_utc))
     goto close_and_return;

   status = 0;
close_and_return:
   (void) TIO_close (ncid);
   return status;
}

typedef struct
{
   const char *name;
   int (*config_method)(Filedb_Type *, config_t *, const char *);
}
Filedb_Method_Type;

#define FILEDB_METHOD(name) {#name,config_##name}
#define FILEDB_METHOD_TABLE_END {NULL,NULL}

static Filedb_Method_Type Filedb_Methods_List[] =
{
   FILEDB_METHOD(met),
   FILEDB_METHOD(snow),
   FILEDB_METHOD(tempo),
   FILEDB_METHOD(ephemeris),
   FILEDB_METHOD_TABLE_END
};

static Filedb_Type *select_method (config_t *cfg, const char *method_name)
{
   Filedb_Method_Type *mt = Filedb_Methods_List;
   Filedb_Type *fdb = NULL;

   if (mt == NULL)
     return 0;

   for (; mt->name != NULL; mt++)
     {
        size_t len = strlen (mt->name);
        if (0 == strncmp (mt->name, method_name, len))
          break;
     }

   if (mt->name == NULL)
     {
        fprintf (stderr, "*** %s: unsupported method: %s\n", __func__, method_name);
        return NULL;
     }

   if (NULL == (fdb = alloc_filedb_type ()))
     return NULL;

   if (0 != mt->config_method (fdb, cfg, method_name))
     {
        free_filedb_type (fdb);
        return NULL;
     }

   return fdb;
}

int main (int argc, char **argv)
{
   config_t cfg;
   char *config_file = "filedb.cfg";
   Filedb_Type *fdb = NULL;
   const char *dbname;
   char *result_filename = NULL;
   double timestamp_utc = DBL_MAX;
   int status = EXIT_FAILURE;
   int task;
   enum
     {
        TASK_UPDATE = 1,
        TASK_FIND = 2
     };
   static struct option long_options[] =
     {
        {"config", required_argument, 0, 'c'},
        {"header", required_argument, 0, 'H'},
        {"sec",    required_argument, 0, 's'},
        {"delay",  required_argument, 0, 'd'},
        {"help",   no_argument,       0, 'h'},
        {"update", no_argument,       0, 'u'},
        {"find",   no_argument,       0, 'f'},
        {0,0,0,0}
     };

   if (argc < 2)
     usage();

   config_init (&cfg);

   /* Try reading the default config file, but if it doesn't exist,
    * keep going in case there's a config file on the command line */
   if (0 == access (config_file, F_OK | R_OK))
     {
        if (-1 == read_config_file (config_file, &cfg))
          goto return_status;
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "c:H:s:d:huf", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "%s: getopt returned character %d??\n",
                      __func__, c);
             goto return_status;
             break;
           case 'c': config_file = optarg;
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (-1 == read_config_file (config_file, &cfg))
               goto return_status;
             break;

           case 'h':
             config_destroy (&cfg);
             usage();
             break;

           case 'd':
             if (1 != sscanf (optarg, "%le", &Delay_Time))
               goto return_status;
             break;

           case 'H':
             if (0 != query_file_timestamp (optarg, &timestamp_utc))
               goto return_status;
             break;

           case 's':
             if (1 != sscanf (optarg, "%le", &timestamp_utc))
               goto return_status;
             break;

           case 'u':
             task = TASK_UPDATE;
             break;
           case 'f':
             task = TASK_FIND;
             break;
          }
     }

   if (optind == argc)
     {
        config_destroy(&cfg);
        usage();
     }

   dbname = argv[optind++];

   if (optind < argc)
     {
        fprintf (stderr, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stderr, "%s ", argv[optind++]);
          }
        fprintf (stderr, "\n");
     }

   if (NULL == (fdb = select_method (&cfg, dbname)))
     goto return_status;

   switch (task)
     {
      case TASK_UPDATE:
        if (0 != fdb_initialize (fdb))
          goto return_status;
        break;

      case TASK_FIND:
        if (NULL == (result_filename = filedb_lookup (fdb->lookup_table, timestamp_utc)))
          goto return_status;
        fputs (result_filename, stdout);
        FREE(result_filename);
        break;

      default:
        fprintf (stderr, "*** Unsupported task: task=%d\n", task);
        goto return_status;
        break;
     }

   status = EXIT_SUCCESS;

return_status:

   free_filedb_type (fdb);
   config_destroy (&cfg);

   return status;
}
