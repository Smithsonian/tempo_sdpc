/** @file   l0_format.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date   Oct 2016
 *  @brief  Main program
 */

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <getopt.h>
#include <signal.h>
#include <wordexp.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "md5.h"
#include "l0_format.h"

#ifndef PROCESSED_FILE_LOG_BASENAME
# define PROCESSED_FILE_LOG_BASENAME "processed_file.log"
#endif

#ifndef PROCESSED_FILE_MAX_LOG_ENTRIES
# define PROCESSED_FILE_MAX_LOG_ENTRIES 100000
#endif

#ifndef CACHE_FLUSH_EXPREC_WAIT_SECS_MINIMUM
# define CACHE_FLUSH_EXPREC_WAIT_SECS_MINIMUM 60.0
#endif

#define SDPC_FILETYPE_MANEUVER   (-100)
#define SDPC_FILETYPE_EPHEMERIS  (-101)

#define MD5_NUM_BYTES 16

static int Have_Epoch;
static int Perform_Archive_Registration;
static int Disable_MD5_Checksum_Files;
static const char *Archive_Root_Dir = NULL;

static const char *Public_Mirror_Root_Dir = NULL;

static int Processing_Version = 1;

static Process_Method_Type *Exprec_Process_Method;
static Process_Method_Type *Iru_Process_Method;

static int caught_signal (void);
static void log_caught_signal (void);

typedef struct
{
   int filetype;
   Process_Method_Type *method;
   Process_Method_Type *(*init)(config_t *);
   Process_Method_Callback_Function *post_process_callback;
}
Process_Method_Table_Type;

typedef struct
{
   double max_iru_knowledge_gap_duration;
   /* Time between telemetry-only granules;
    * A negative value means "don't generate telemetry-only granules".
    */
   double latest_iru_timestamp_seen;
   double latest_radiance_timestamp_seen;
   double latest_iru_only_interval_end_time;
   int need_first_iru_timestamp;
   char *dir;
}
IRU_Interval_Type;

typedef struct
{
   const char *logdir;
   char *path;
   FILE *fp;
   size_t curr_num_file_entries;
}
Processed_File_Log_Type;

typedef struct
{
   const char *input_filename_glob_pattern;
   char *incoming_dir;
   char *tpinfo_file;
   double monitor_wait_secs;
   double start_time;
   double stop_time;
   double cache_flush_exprec_wait_secs;
   int exit_on_emptydir;
   IRU_Interval_Type iru_interval;
   Processed_File_Log_Type *log_incoming;
}
Control_Type;

static double First_Packet_Time;
static double Last_Packet_Time;

static void usage (void)
{
   fprintf (stderr, "Usage: L0_format [options] [config-file]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -e | --empty             Exit when the input directory is empty\n");
   fprintf (stderr, "   -a | --archive DIR       Archive files in directory DIR\n");
   fprintf (stderr, "   -m | --mirror DIR        Public mirror files in directory DIR\n");
   fprintf (stderr, "   -n | --nochecksum        Disable generation of product MD5 checksum files\n");
   fprintf (stderr, "   -L | --logdir DIR        Log processed files in directory DIR\n");
   fprintf (stderr, "   -r | --register          Perform database registration of archived files\n");
   fprintf (stderr, "   -c | --cache DIR         Process cached directories matching regex DIR\n");
   fprintf (stderr, "                            e.g. --cache 'd710[1-4]/h[0-2][0-9]'\n");
   fprintf (stderr, "                            If DIR begins with '@', it's the path\n");
   fprintf (stderr, "                            to a file containing a regex list\n");
   fprintf (stderr, "   -t | --tstart SEC        Process cache files newer than SEC since the TEMPO epoch\n");
   fprintf (stderr, "   -V | --Version N         Processing version number\n");
   fprintf (stderr, "   -v | --verbose           Increase verbosity (-vv is more verbose)\n");
   exit (EXIT_SUCCESS);
}

int get_processing_version (void)
{
   return Processing_Version;
}

char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
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

static void set_archive_root_dir (const char *dir)
{
   Archive_Root_Dir = dir;
}

static const char *get_archive_root_dir (void)
{
   return Archive_Root_Dir;
}

static void set_public_mirror_root_dir (const char *dir)
{
   Public_Mirror_Root_Dir = dir;
}

static int close_processed_file_log (Processed_File_Log_Type *flt)
{
   time_t now = {0};
   struct tm tm = {0};
   struct stat st = {0};
   char suffix[64];
   char *new_path = NULL;
   int status;
   size_t len;

   if (flt == NULL)
     return 0;

   if (flt->fp)
     {
        if (0 != fclose (flt->fp))
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: closing: %s ", __func__, flt->path);
             return -1;
          }
     }
   flt->fp = NULL;

   /* If no file exists we're done */
   if (0 != stat (flt->path, &st))
     return 0;

   /* Rename the existing file by appending a timestamp */
   time(&now);
   gmtime_r(&now, &tm);
   if (0 == strftime (suffix, sizeof(suffix), ".%Y%m%dT%H%M%S", &tm))
     return -1;

   len = strlen (flt->path) + strlen(suffix) + 1;
   if (NULL == (new_path = (char *)MALLOC (len)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   strcpy (new_path, flt->path);
   strcat (new_path, suffix);

   status = rename (flt->path, new_path);
   ioclib_free (new_path);
   if (status != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: rename failed: %s ", __func__, flt->path);
        return -1;
     }

   return 0;
}

static int finalize_processed_file_log (Processed_File_Log_Type *flt)
{
   int status = close_processed_file_log (flt);
   if (flt)
     {
        ioclib_free (flt->path);
     }
   return status;
}

static int new_processed_file_log (Processed_File_Log_Type *flt)
{
   if (flt->path == NULL)
     {
        if (NULL == (flt->path = ioclib_pathconcat (flt->logdir, PROCESSED_FILE_LOG_BASENAME)))
          return -1;
     }

   if (0 != close_processed_file_log (flt))
     return -1;

   if (NULL == (flt->fp = fopen (flt->path, "w")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening log file for writing: %s ", __func__, flt->path);
        return -1;
     }
   flt->curr_num_file_entries = 0;

   return 0;
}

static int open_processed_file_log (Processed_File_Log_Type *flt)
{
   if (flt == NULL)
     return 0;
   if (flt->logdir)
     {
        if (0 != ioclib_mkdir (flt->logdir, 0))
          return -1;
     }
   return new_processed_file_log (flt);
}

static int flush_processed_file_log (Processed_File_Log_Type *flt)
{
   if (flt == NULL)
     return 0;

   return fflush (flt->fp);
}

static int log_processed_file (Processed_File_Log_Type *flt, const char *path,
                               time_t mtime_tv_sec, int status)
{
   struct timeval tv = {0};

   if (flt == NULL)
     return 0;

   if (flt->curr_num_file_entries >= PROCESSED_FILE_MAX_LOG_ENTRIES)
     {
        if (0 != new_processed_file_log (flt))
          return -1;
     }

   if (flt->curr_num_file_entries == 0)
     {
        char *dirname;
        if (NULL == (dirname = ioclib_dirname (path)))
          return -1;
        (void) fprintf (flt->fp, "# %s\n", dirname);
        ioclib_free (dirname);
     }

   (void) gettimeofday (&tv, NULL);

   if (fprintf (flt->fp, "%ld,%06ld,%d,%s,%ld\n",
                tv.tv_sec, tv.tv_usec, status, ioclib_basename (path), mtime_tv_sec) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing to log file: %s", __func__, flt->path);
        return -1;
     }

   flt->curr_num_file_entries++;

   if (0 == (flt->curr_num_file_entries % 10))
     {
        (void) flush_processed_file_log (flt);
     }

   return 0;
}

static void free_control_type_fields (Control_Type *ctrl)
{
   FREE(ctrl->incoming_dir);
   FREE(ctrl->tpinfo_file);
   FREE(ctrl->iru_interval.dir);
}

static int read_main_params (config_t *cfg, Control_Type *ctrl)
{
   config_setting_t *s;
   const char *incoming_dir;
   const char *tpinfo_file;

   if (NULL == (s = config_lookup (cfg, "main")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'main' in param file: %s",
                     __func__, config_error_file(cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "incoming_dir", &incoming_dir))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "input_filename_glob_pattern", &ctrl->input_filename_glob_pattern))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "tpinfo_file", &tpinfo_file))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "monitor_wait_secs", &ctrl->monitor_wait_secs))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "start_time", &ctrl->start_time))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "cache_flush_exprec_wait_secs", &ctrl->cache_flush_exprec_wait_secs))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'main' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (ctrl->cache_flush_exprec_wait_secs < CACHE_FLUSH_EXPREC_WAIT_SECS_MINIMUM)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: %s: cache_flush_exprec_wait_secs = %0.3f sec (minimum valid value = %0.3f sec)",
                     __func__, config_error_file (cfg),
                     ctrl->cache_flush_exprec_wait_secs,
                     CACHE_FLUSH_EXPREC_WAIT_SECS_MINIMUM);
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "stop_time", &ctrl->stop_time))
     ctrl->stop_time = 0;

   if ((NULL == (ctrl->incoming_dir = expand_string (incoming_dir)))
       || (NULL == (ctrl->tpinfo_file = expand_string (tpinfo_file))))
     {
        return -1;
     }

   return 0;
}

static int read_exprec_params (config_t *cfg, Control_Type *ctrl)
{
   config_setting_t *s;
   const char *exprec_out_dirname;

   if (NULL == (s = config_lookup (cfg, "exprec")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'exprec' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &exprec_out_dirname))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading parameter 'exprec:output_dir' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (ctrl->iru_interval.dir = expand_string (exprec_out_dirname)))
     return -1;

   return 0;
}

static int read_iru_params (config_t *cfg, Control_Type *ctrl)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "iru")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'iru' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "inr_update_interval", &ctrl->iru_interval.max_iru_knowledge_gap_duration))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading parameter 'iru:inr_update_interval' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (ctrl->iru_interval.max_iru_knowledge_gap_duration >= 0)
     ctrl->iru_interval.need_first_iru_timestamp = 1;
   else
     ctrl->iru_interval.need_first_iru_timestamp = 0;

   return 0;
}

static int parse_param_file (config_t *cfg, const char *cfg_file,
                             Control_Type *ctrl)
{
   if (0 == config_read_file (cfg, cfg_file))
    {
       tell_verror (TELL_INVALID_PARM_ERROR, "%s: Reading %s:%d - %s",
                    __func__, config_error_file(cfg),
                    config_error_line(cfg), config_error_text(cfg));
       return -1;
    }

   if (0 != read_main_params (cfg, ctrl))
     return -1;

   if (0 != read_exprec_params (cfg, ctrl))
     return -1;

   if (0 != read_iru_params (cfg, ctrl))
     return -1;

   return 0;
}

static Process_Method_Type *
find_process_method (const Process_Method_Table_Type *tbl, int filetype)
{
   for ( ; tbl->init != NULL; tbl++)
     {
        if (tbl->filetype == filetype)
          return tbl->method;
     }

   tell_verror (TELL_APPLICATION_ERROR,
                "%s: unsupported file type: filetype=%d", __func__, filetype);

   return NULL;
}

static int make_iru_only_path (const char *dir, double tstart, char *path, size_t pathsize)
{
   size_t n, strsize;
   char *str = path;

   if (dir == NULL)
     {
        str = path;
        strsize = pathsize;
     }
   else
     {
        n = snprintf (path, pathsize, "%s/", dir);
        if (n >= pathsize)
          {
             tell_verror (TELL_APPLICATION_ERROR,
                          "%s: path length %ld truncated to buffer size %ld)",
                          __func__, n, pathsize);
             return -1;
          }
        str = path + n;
        strsize = pathsize - n;
     }

   n = __tio_filename_string (str, strsize, tstart, TEMPO_PROD_TYPE_INR, 0, 0);
   if (n >= strsize)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: basename length %ld truncated to buffer size %ld)",
                     __func__, n, strsize);
        return -1;
     }

   return 0;
}

static int write_iru_only_interval (const char *dir, double tbeg, double tend)
{
   char path[MAX_PATHLEN];
   time_t epoch;
   double tbeg_utc;
   FILE *fp;

   if (0 != tio_time_taix_to_utc (tbeg, &tbeg_utc))
     return -1;

   if (0 != make_iru_only_path (dir, tbeg, path, sizeof(path)))
     return -1;

   tell_vinfo (0, "creating file %s", path);

   if (NULL == (fp = fopen (path, "w")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: error opening file %s", __func__, path);
        return -1;
     }

   tell_vinfo (0, "tbeg=%0.6f tend=%0.6f", tbeg, tend);

   epoch = (time_t) tio_time_taix_epoch_timet();

   if (fprintf (fp, "%0.6f,%0.6f,%ld,%ld\n", tbeg, tend, epoch, (time_t) tbeg_utc) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: error writing to file %s", __func__, path);
        (void) fclose (fp);
        return -1;
     }

   if (0 != fclose (fp))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: error closing file %s", __func__, path);
        return -1;
     }

   return 0;
}

static int ensure_iru_coverage_for_inr (IRU_Interval_Type *iru_interval, double t_beg, double t_end)
{
   double dt_max = iru_interval->max_iru_knowledge_gap_duration;
   double dt_gap, dt_file;
   int i, n;

   /* Decide how many files */
   dt_gap = t_end - t_beg;
   n = (dt_gap < dt_max) ? 1 : (dt_gap / dt_max);
   dt_file = dt_gap / n;

   for (i = 0; i < n; i++)
     {
        double tbeg = t_beg + i*dt_file;
        double tend = tbeg + dt_file;
        if (0 != write_iru_only_interval (iru_interval->dir, tbeg, tend))
          return -1;
        iru_interval->latest_iru_only_interval_end_time = tend;
     }

   return 0;
}

static int exprec_post_process_callback (Process_Method_Type *pmt, void *client_data)
{
   IRU_Interval_Type *iru_interval = (IRU_Interval_Type *)client_data;
   double dt_max = iru_interval->max_iru_knowledge_gap_duration;
   double t_only = iru_interval->latest_iru_only_interval_end_time;
   double t_rad_prev = iru_interval->latest_radiance_timestamp_seen;
   double t_rad, t_last;

   /* dt_max < 0 means "don't generate telemetry-only granules" */
   if (dt_max < 0.0)
     return 0;

   t_last = (t_rad_prev > t_only) ? t_rad_prev : t_only;

   if (0 != pmt->pmt_query_latest_timestamp (pmt, IOCSDPC_EXPREC_TYPE_RAD, &t_rad))
     return -1;

   if (t_last > 0)
     {
        double dt_only_min = 3;           /* [sec] min duration desired for a single telemetry-only file */
        double dt_rad_pad = 10;           /* [sec] padding that preceeds every radiance granule */
        double dt_only = t_rad - t_only;
        double dt_rad = t_rad - t_rad_prev;

        if ((dt_rad > dt_rad_pad)
            && (dt_only_min < dt_only) && (dt_only < dt_max))
          {
             /* Fill small gap before first radiance scan */
             if (0 != ensure_iru_coverage_for_inr (iru_interval, t_only, t_rad))
               return -1;
          }
        else if ((dt_rad > dt_rad_pad) && (dt_only > dt_max))
          {
             /* Fill largish gap between radiance scans */
             if (0 != ensure_iru_coverage_for_inr (iru_interval, t_rad_prev, t_rad))
               return -1;
          }
        else if (t_rad - t_last > dt_max)
          {
             if (0 != ensure_iru_coverage_for_inr (iru_interval, t_last, t_rad))
               return -1;
          }
     }

   iru_interval->latest_radiance_timestamp_seen = t_rad;

   return 0;
}

static int iru_post_process_callback (Process_Method_Type *pmt, void *client_data)
{
   IRU_Interval_Type *iru_interval = (IRU_Interval_Type *)client_data;
   double dt_max = iru_interval->max_iru_knowledge_gap_duration;
   double t_iru, t_last, t_only, t_rad;

   /* dt_max < 0 means "don't generate telemetry-only granules" */
   if (dt_max < 0.0)
     return 0;

   if (0 != pmt->pmt_query_latest_timestamp (pmt, 0, &t_iru))
     return -1;
   iru_interval->latest_iru_timestamp_seen = t_iru;

   /* iru_interval->latest_iru_only_interval_end_time is initialized
    * when the first IRU file is opened. Therefore, when we get to this
    * point, we're guaranteed that its value is >= 0, so we have t_only >= 0.
    * If no radiance data has been seen yet, then we might have t_rad <=0.
    * Avoid gaps by padding the IRU coverage in every granule.
    */
   t_only = iru_interval->latest_iru_only_interval_end_time;
   t_rad = iru_interval->latest_radiance_timestamp_seen;

   /* What's the latest IRU timestamp headed for the INR subsystem? */
   t_last = (t_rad > t_only) ? t_rad : t_only;

   /* Does the INR subsystem need an IRU update? */
   if (t_iru - t_last > dt_max)
     {
        double t_end = t_last + dt_max;
        if (0 != ensure_iru_coverage_for_inr (iru_interval, t_last, t_end))
          return -1;
     }

   return 0;
}

static int get_first_iru_sample_time (const char *file, int fd, IOCSDPC_Common_Header_Type *chdr,
                                      double *sample_time)
{
   IOCSDPC_IRU_Type *iru = NULL;
   IOCSDPC_IRU_Record_Type rec = {0};

   *sample_time = -1.0;

   if (NULL == (iru = iocsdpc_iru_fdopen_read (file, fd, chdr)))
     return -1;

   for (;;)
     {
        unsigned int num_read;
        if ((0 != iocsdpc_iru_read (iru, &rec, 1, &num_read))
            || (num_read != 1))
          return -1;
        if (rec.sample_time > 0.0)
          break;
     }

   *sample_time = rec.sample_time;

   return 0;
}

static int classify_file (const char *file, Control_Type *ctrl, int *filetype, int *skip)
{
   IOCSDPC_Common_Header_Type chdr = {0};
   char *ext;
   int fd;

   /* The file stream may contain:
    *  - maneuver files: <prefix>_maneuver.csv
    *  - ephemeris files: <prefix>_ephemeris.csv
    *  - L0 sciextract-produced data products <prefix>_*.*
    * where prefix looks like tempo_dDDDDDmMMMMMMMMuUUU_rR
    */

   *skip = 0;

   if ((NULL != (ext = ioclib_extname (file)))
       && (0 == strcmp (ext, ".csv")))
     {
        char *basename = ioclib_basename (file);
        if (NULL != strstr (basename, "_maneuver.csv"))
          {
             *filetype = SDPC_FILETYPE_MANEUVER;
             return 0;
          }
        else if (NULL != strstr (basename, "_ephemeris.csv"))
          {
             *filetype = SDPC_FILETYPE_EPHEMERIS;
             return 0;
          }
        /* reject unrecognized file */
        tell_vinfo (0, "%s: unrecognized CSV file type: %s", __func__, file);
        return -1;
     }

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;
   First_Packet_Time = chdr.first_packet_time;
   Last_Packet_Time = chdr.last_packet_time;

   /* To handle a corner case when generating telemetry-only radiance
    * files, initialize iru_interval.latest_iru_only_interval_end_time
    * to the first IRU sample time >= 0.
    */
   if ((chdr.filetype == IOCSDPC_FILETYPE_IRU)
       && (ctrl->iru_interval.need_first_iru_timestamp != 0))
     {
        double sample_time;
        if ((0 != get_first_iru_sample_time (file, fd, &chdr, &sample_time))
            || (sample_time < 0.0))
          {
             (void) ioclib_fd_close (fd);
             return -1;
          }
        /* At this point, we know that we've never sent any IRU data to the
         * INR subsystem, so the current IRU knowledge gap starts here. */
        ctrl->iru_interval.latest_iru_only_interval_end_time = sample_time;
        ctrl->iru_interval.need_first_iru_timestamp = 0;
     }

   (void) ioclib_fd_close (fd);

   *filetype = chdr.filetype;

   /* We may want to skip files prior to some user-specified time */
   if (ctrl->start_time > 0)
     {
        if (0 != verify_epoch (chdr.epoch))
          return -1;
        *skip = (chdr.last_packet_time < ctrl->start_time);
     }

   return 0;
}

static int process_file (const Process_Method_Table_Type *tbl, const TPInfo_Type *tpinfo,
                         Control_Type *ctrl, const char *file)
{
   Process_Method_Type *pmt;
   struct stat st = {0};
   int filetype, skip, status;

   if (0 != lstat (file, &st))
     return -1;

   if (0 != classify_file (file, ctrl, &filetype, &skip))
     return -1;

   if (skip)
     {
        tell_vinfo (1, "skipped: %s", file);
        return 0;
     }

   if (NULL == (pmt = find_process_method (tbl, filetype)))
     return -1;

   tell_vinfo (1, "processing: %s", file);
   status = pmt->pmt_process (pmt, tpinfo, file, &ctrl->iru_interval);

   /* Complain when logging fails, but don't stop processing */
   (void) log_processed_file (ctrl->log_incoming, file, st.st_mtim.tv_sec, status);

   return status;
}

static int maybe_flush_exprec_cache (const TPInfo_Type *tpinfo, Control_Type *ctrl);

static int process_live_stream_dir_files (const Process_Method_Table_Type *tbl,
                                          const TPInfo_Type *tpinfo, Control_Type *ctrl,
                                          char **file_list, size_t num_files)
{
   size_t i;

   for (i = 0; i < num_files; i++)
     {
        char *file = file_list[i];

        if (caught_signal())
          break;

        if (0 != process_file (tbl, tpinfo, ctrl, file))
          {
             tell_vinfo (0, "%s: bad file: %s", __func__, file);
             if (0 != ioclib_rename_to_bad_file (file))
               {
                  tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_rename_to_bad_file failed, file=%s",
                               __func__, file ? file : "(null)");
                  return -1;
               }
             continue;
          }

        /* If processing involved a rename, then deletion
         * will be handled elsewhere. */
        if (ioclib_isfile (file, NULL))
          {
             if (0 != ioclib_unlink (file))
               return -1;
          }

        if (0 != maybe_flush_exprec_cache (tpinfo, ctrl))
          return -1;
     }

   return 0;
}

static int flush_caches (const Process_Method_Table_Type *tbl, const TPInfo_Type *tpinfo,
                         int unwind, const char *incoming_dir)
{
   Process_Method_Type *pmt;
   int num_failed = 0;

   for (; tbl->init != NULL; tbl++)
     {
        pmt = tbl->method;
        if (pmt->pmt_flush_cache)
          {
             if (0 != pmt->pmt_flush_cache (pmt, tpinfo, unwind, incoming_dir))
               num_failed++;
          }
     }

   return num_failed;
}

static int flush_iru_coverage_for_inr (IRU_Interval_Type *iru_interval)
{
   Process_Method_Type *pmt = Iru_Process_Method;
   double dt_max = iru_interval->max_iru_knowledge_gap_duration;
   double t_iru, t_last, t_only, t_rad;

   /* dt_max < 0 means "don't generate telemetry-only granules" */
   if (dt_max < 0.0)
     return 0;

   if (0 != pmt->pmt_query_latest_timestamp (pmt, 0, &t_iru))
     return -1;
   iru_interval->latest_iru_timestamp_seen = t_iru;

   /* If we haven't seen any IRU data, there's nothing to flush */
   if (t_iru < 0.0)
     return 0;

   /* We've seen IRU data, so we have t_only >= 0.0 */
   t_only = iru_interval->latest_iru_only_interval_end_time;
   t_rad = iru_interval->latest_radiance_timestamp_seen;

   /* What's the latest IRU timestamp headed for the INR subsystem? */
   t_last = (t_rad > t_only) ? t_rad : t_only;

   /* Flush whatever we have */
   if (t_iru > t_last)
     {
        if (0 != ensure_iru_coverage_for_inr (iru_interval, t_last, t_iru))
          return -1;
     }

   return 0;
}

static int maybe_flush_exprec_cache (const TPInfo_Type *tpinfo, Control_Type *ctrl)
{
   Process_Method_Type *pmt = Exprec_Process_Method;
   double last_erec_cached_timestamp, age_secs;

   /* If exposure records are cached, but there's a gap in telemetry since we've seen one,
    * that likely indicates that it's time to flush the cache and close that file.
    * This is an important mechanism for triggering DRK processing soon after DRK records
    * stop arriving, but a similar situation can arrive for other exposure record types.
    * It's important to check after every successfully processed file
    * to ensure that we catch a telemetry gap between any two files.
    */

   /* If we're not processing exposure records, there's nothing more to do here */
   if (pmt == NULL)
     return 0;

   if (0 != pmt->pmt_query_last_erec_cached_timestamp (pmt, &last_erec_cached_timestamp))
     return -1;

   /* If the cache is empty, we're done */
   if (last_erec_cached_timestamp <= 0.0)
     return 0;

   /* Incoming data stream is time-ordered, so once we receive a file beginning
    * with timestamp, T, all subsequent timestamps, t, should be >=T */
   age_secs = First_Packet_Time - last_erec_cached_timestamp;

   /* If exposure records arrived recently, then it's too soon to flush the exprec cache. */
   if (age_secs < ctrl->cache_flush_exprec_wait_secs)
     return 0;

   tell_vinfo (0, "flush exprec cache (first_packet_time - last_erec_time = %0.3f sec > cache_flush_exprec_wait_secs = %0.3f sec)",
               age_secs, ctrl->cache_flush_exprec_wait_secs);

   return pmt->pmt_flush_cache (pmt, tpinfo, 0, NULL);
}

#define PROCESS_METHOD(filetype,init,callback) {filetype,NULL,init,callback}
#define PROCESS_METHODS_TABLE_END {-1,NULL,NULL,NULL}

static Process_Method_Table_Type Method_Table[] =
{
   PROCESS_METHOD(IOCSDPC_FILETYPE_EXPREC, init_exprec_method, exprec_post_process_callback),
   PROCESS_METHOD(IOCSDPC_FILETYPE_TPSEC, init_tpsec_method, NULL),
   PROCESS_METHOD(IOCSDPC_FILETYPE_IRU, init_iru_method, iru_post_process_callback),
   PROCESS_METHOD(IOCSDPC_FILETYPE_SMC, init_smc_method, NULL),
   PROCESS_METHOD(SDPC_FILETYPE_MANEUVER, init_maneuver_method, NULL),
   PROCESS_METHOD(SDPC_FILETYPE_EPHEMERIS, init_ephemeris_method, NULL),
   PROCESS_METHODS_TABLE_END
};

static int init_methods_table (Process_Method_Table_Type *tbl,
                               config_t *cfg)
{
   if (tbl == NULL)
     return -1;

   for ( ; tbl->init != NULL; tbl++)
     {
        tbl->method = (*tbl->init)(cfg);
        if (NULL == tbl->method)
          return -1;

        tbl->method->pmt_post_process_callback = tbl->post_process_callback;

        if (tbl->filetype == IOCSDPC_FILETYPE_EXPREC)
          Exprec_Process_Method = tbl->method;
        else if (tbl->filetype == IOCSDPC_FILETYPE_IRU)
          Iru_Process_Method = tbl->method;
     }

   return 0;
}

static void delete_methods_table (Process_Method_Table_Type *tbl)
{
   if (tbl == NULL)
     return;

   for ( ; tbl->init != NULL; tbl++)
     {
        Process_Method_Type *m = tbl->method;
        if ((m != NULL) && (m->pmt_delete != NULL))
          {
             m->pmt_delete (m);
          }
     }
}

static volatile int Sigterm_Received;
static void sigterm_handler (int sig)
{
   (void) sig;
   Sigterm_Received++;
}
static void catch_sigterm (void)
{
   struct sigaction new_action;
   new_action.sa_handler = sigterm_handler;
   sigemptyset (&new_action.sa_mask);
   new_action.sa_flags = 0;
   sigaction (SIGTERM, &new_action, NULL);
}

static volatile int Sigint_Received;
static void sigint_handler (int sig)
{
   (void) sig;
   Sigint_Received++;
}
static void catch_sigint (void)
{
   struct sigaction new_action;
   new_action.sa_handler = sigint_handler;
   sigemptyset (&new_action.sa_mask);
   new_action.sa_flags = 0;
   sigaction (SIGINT, &new_action, NULL);
}

static void log_caught_signal (void)
{
   const char *signame;

   if (Sigterm_Received)
     signame = "SIGTERM";
   else if (Sigint_Received)
     signame = "SIGINT";
   else
     signame = "unknown";

   tell_vinfo (0, "caught signal: %s (pid=%d)", signame, getpid());
}

static int caught_signal (void)
{
   return (Sigterm_Received || Sigint_Received);
}

static int process_live_stream (Process_Method_Table_Type *tbl,
                                const TPInfo_Type *tpinfo, Control_Type *ctrl)
{
   IOCLib_Glob_Type *gt = NULL;
   char *pattern = NULL;
   int received_signal;
   int status = -1;

   pattern = ioclib_pathconcat (ctrl->incoming_dir,
                                ctrl->input_filename_glob_pattern);
   if (NULL == pattern)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed",
                    __func__);
        return -1;
     }

   tell_vinfo (0, "processing %s", pattern);
   if (ctrl->start_time > 0)
     {
        tell_vinfo (0, "start time = %f sec since the epoch", ctrl->start_time);
     }

   while (0 == (received_signal = caught_signal()))
     {
        ioclib_glob_free (gt);
        gt = NULL;
        if (NULL == (gt = ioclib_glob (pattern, 0)))
          {
             tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_glob failed", __func__);
             goto return_status;
          }

        if (ctrl->exit_on_emptydir && (gt->num_files == 0))
          break;

        if (-1 == process_live_stream_dir_files (tbl, tpinfo, ctrl, gt->files, gt->num_files))
          goto return_status;

        (void) ioclib_sleep (ctrl->monitor_wait_secs);
     }

   if (received_signal != 0)
     log_caught_signal();

   tell_vinfo (0, "flush caches on exit");
   if (0 != flush_caches (tbl, tpinfo, received_signal, ctrl->incoming_dir))
     goto return_status;
   if (0 != flush_iru_coverage_for_inr (&ctrl->iru_interval))
     goto return_status;
   (void) flush_processed_file_log (ctrl->log_incoming);

   status = 0;
return_status:
   ioclib_free (pattern);
   ioclib_glob_free (gt);

   return status;
}

static int process_cache_dir_pattern (Process_Method_Table_Type *tbl,
                                      const TPInfo_Type *tpinfo, Control_Type *ctrl,
                                      const char *cache_dir_pattern)
{
   IOCLib_Glob_Type *gt = NULL;
   char **files = NULL;
   char *path = NULL;
   int status = -1;
   size_t i, k, num_files;

   if (NULL == (gt = ioclib_glob (cache_dir_pattern, 0)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_glob failed", __func__);
        return -1;
     }

   for (i = 0; i < gt->num_files; i++)
     {
        const char *dir = gt->files[i];

        if (NULL == (files = ioclib_dir_list(dir, &num_files, IOCLIB_LISTDIR_SORT)))
          goto return_status;

        tell_vinfo (0, "begin processing %ld files from directory %s", num_files, dir);
        if (ctrl->start_time > 0)
          {
             tell_vinfo (0, "start time = %f sec since the epoch", ctrl->start_time);
          }

        for (k = 0; k < num_files; k++)
          {
             if (NULL == (path = ioclib_pathconcat (dir, files[k])))
               goto return_status;

             if (caught_signal())
               {
                  log_caught_signal();
                  status = 0;
                  goto return_status;
               }

             (void) process_file (tbl, tpinfo, ctrl, path);

             if (0 != maybe_flush_exprec_cache (tpinfo, ctrl))
               goto return_status;

             if ((ctrl->stop_time > 0) && (Last_Packet_Time > ctrl->stop_time))
               {
                  tell_vinfo (0, "stopping:  last packet time=%f  exceeds specified stop_time=%f",
                              Last_Packet_Time, ctrl->stop_time);
                  goto last_packet_time_exceeds_stop_time;
               }

             ioclib_free (path);
             path = NULL;
          }

        tell_vinfo (0, "end processing %ld files from directory %s", num_files, dir);

        ioclib_string_array_free (files, num_files);
        files = NULL;
        num_files = 0;
     }

last_packet_time_exceeds_stop_time:
   status = 0;
return_status:
   ioclib_glob_free (gt);
   ioclib_string_array_free (files, num_files);
   ioclib_free (path);

   return status;
}

static int process_cache_dirs (Process_Method_Table_Type *tbl,
                               const TPInfo_Type *tpinfo, Control_Type *ctrl,
                               const char *cache_dir_arg)
{
   int received_signal = 0;

   if (*cache_dir_arg != '@')
     {
        if (0 != process_cache_dir_pattern (tbl, tpinfo, ctrl, cache_dir_arg))
          return -1;
     }
   else
     {
        FILE *fp = NULL;
        const char *path = cache_dir_arg + 1;

        if (NULL == (fp = fopen (path, "r")))
          {
             tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s", __func__, path);
             return -1;
          }

        for (;;)
          {
             char *newline;
             char buf[1024];

             if (0 != (received_signal = caught_signal()))
               break;

             if (NULL == fgets (buf, sizeof(buf), fp))
               break;
             /* Assume no leading whitespace, no line-breaks within regex,
              * each line ends with a single newline
              */
             if (NULL != (newline = strchr (buf, '\n')))
               *newline = 0;

             /* leading # is a comment character */
             if (buf[0] == '#')
               continue;

             if (0 != process_cache_dir_pattern (tbl, tpinfo, ctrl, buf))
               {
                  fclose (fp);
                  return -1;
               }
          }

        fclose (fp);
     }

   tell_vinfo (0, "flush caches on exit");
   /* If processing was interrupted, we may want to delete an incomplete image granule */
   if (0 != flush_caches (tbl, tpinfo, received_signal, NULL))
     return -1;
   if (0 != flush_iru_coverage_for_inr (&ctrl->iru_interval))
     return -1;
   (void) flush_processed_file_log (ctrl->log_incoming);
   return 0;
}

int verify_epoch (time_t epoch)
{
   if (Have_Epoch)
     {
        double current_epoch = tio_time_taix_epoch_timet ();
        if (epoch != current_epoch)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: epoch mismatch: old=%f  new=%ld",
                          __func__, current_epoch, epoch);
             return -1;
          }
     }
   else
     {
        if (0 != tio_time_set_taix_epoch_timet (epoch))
          return -1;
        Have_Epoch = 1;
     }

   return 0;
}

int make_level0_archdir_path (char **archdir_path, double sec_since_epoch, int scan_num,
                              const char *suffix)
{
   char buf[MAX_PATHLEN];
   size_t bufsize = sizeof(buf);
   const char *root_path;
   char *path = NULL;
   double sat_day;
   size_t n;

   /* NULL means don't perform archiving */
   if (NULL == (root_path = get_archive_root_dir ()))
     {
        *archdir_path = NULL;
        return 0;
     }

   /* Number of days since the TEMPO epoch, spacecraft local time.
    * Spacecraft local time is used because it makes the archive organization
    * more intuitive.  To force UTC time in the archive, set SC_Timezone=0.
    */
   if (0 != tio_time_sat_local_day_number (sec_since_epoch, &sat_day))
     return -1;

   /* e.g. ${SDPC_ARCHIVE_DIR}/L0/Ddddd/${file_type}
    *   or ${SDPC_ARCHIVE_DIR}/L0/Ddddd/${file_type}/Sddd
    */
   if (scan_num < 0)
     {
        n = snprintf (buf, bufsize, "%s/L0/D%05d/%s",
                      root_path, (int) sat_day, suffix);
     }
   else
     {
        n = snprintf (buf, bufsize, "%s/L0/D%05d/%s/S%03d",
                      root_path, (int) sat_day, suffix, scan_num);
     }

   if (n >= bufsize)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: basename length %ld truncated to buffer size %ld)",
                     __func__, n, bufsize);
        return -1;
     }

   if (NULL == (path = strdup (buf)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
        return -1;
     }

   *archdir_path = path;

   return 0;
}

int make_level0_basename (char *buf, int bufsize,
                          double sec_since_epoch, int processing_version,
                          const char *suffix, const Radiance_Ident_Type *identp)
{
   int n, level = 0;

   if (identp)
     {
        n = __tio_filename_string_indexed (buf, bufsize,
                                           sec_since_epoch, suffix, level, processing_version,
                                           identp->scan_num, identp->granule_num);
     }
   else
     {
        n = __tio_filename_string (buf, bufsize,
                                   sec_since_epoch, suffix, level, processing_version);
     }

   if (n >= bufsize)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: basename length %d truncated to buffer size %d)",
                     __func__, n, bufsize);
        return -1;
     }

   return 0;
}

static int make_hidden_basename (const char *basename, char *buf, int bufsize)
{
   if (snprintf (buf, bufsize, ".%s", basename) >= bufsize)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: basename exceeds %d char buffer size",
                     __func__, bufsize);
        return -1;
     }
   return 0;
}

int write_attr_global_timestamp (int ncid, const char *tstamp_name,
                                 double tstamp_value)
{
   return TIO_write_timestamp (ncid, NC_GLOBAL, tstamp_name, tstamp_value);
}

int write_std_global_metadata (int ncid, const IOCSDPC_Common_Header_Type *chdr)
{
   int chdr_content_version = (chdr != NULL) ? chdr->content_version : 0;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "content_version", NC_INT, 1, &chdr_content_version))
     return -1;
   return tio_write_epoch_timestamp (ncid, NC_GLOBAL);
}

int create_hidden (const char *dirname, const char *basename, int *ncid)
{
   char hidden_basename[MAX_BASENAME_SIZE];
   char *path = NULL;
   int status = 0;

   *ncid = INT_MAX;

   if (-1 == make_hidden_basename (basename, hidden_basename, MAX_BASENAME_SIZE))
     return -1;

   tell_vinfo (0, "creating file %s/%s", dirname, hidden_basename);

   if ((NULL == (path = ioclib_pathconcat (dirname, hidden_basename)))
       || (-1 == TIO_create (path, NC_NETCDF4, ncid))
      )
     {
        status = -1;
     }

   FREE(path);
   return status;
}

/* From the ioblksize.h header in GNU coreutils, a 128 kiB block size
 * minimizes system overhead when copying files on a wide variety
 * of computer systems.  GNU cp uses this.
 */
enum {IO_BUFSIZE = 128*1024};
static inline size_t io_blksize (struct stat *sb)
{
   size_t blksize = sb->st_blksize;
   return (IO_BUFSIZE > blksize) ? IO_BUFSIZE : blksize;
}

static int copy_file (const char *from, const char *to)
{
   mode_t mode_create = 00644;  /* rw-r--r-- */
   struct stat st = {0};
   char *buf = NULL;
   size_t bufsize;
   ssize_t nread;
   int fd_from, fd_to = -1, status = -1;

   if (0 != stat (from, &st))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: cannot stat %s", __func__, from);
        return -1;
     }

   bufsize = io_blksize(&st);

   if ((fd_from = open (from, O_RDONLY)) < 0)
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: opening %s", __func__, from);
        return -1;
     }
   if (0 != posix_fadvise(fd_from, 0, 0, POSIX_FADV_SEQUENTIAL))
     {
        tell_vwarn (0, "%s: posix_fadvise failed: %s", __func__, strerror(errno));
     }

   if ((fd_to = open (to, O_WRONLY | O_CREAT | O_EXCL, mode_create)) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: opening %s (%s)", __func__, to, strerror (errno));
        goto return_status;
     }
   if (0 != posix_fallocate (fd_to, 0, st.st_size))
     {
        tell_vwarn (0, "%s: posix_fallocate failed: %s", __func__, strerror(errno));
     }

   if (NULL == (buf = (char *) MALLOC (bufsize)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   while ((nread = read (fd_from, buf, bufsize)) > 0)
     {
        char *pbuf = buf;
        ssize_t nwritten;

        do
          {
             if ((nwritten = write (fd_to, pbuf, nread)) >= 0)
               {
                  nread -= nwritten;
                  pbuf += nwritten;
               }
             else if (errno != EINTR)
               {
                  tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s (%s)", __func__, to, strerror(errno));
                  goto return_status;
               }
          } while (nread > 0);
     }

   if (nread == 0)
     {
        if (close (fd_to) < 0)
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: closing %s (%s)", __func__, to, strerror(errno));
             fd_to = -1;
             goto return_status;
          }
        /* Success! */
        status = 0;
     }
   else
     {
        /* Error: nread < 0 */
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %s (%s)", __func__, from, strerror(errno));
     }

return_status:
   FREE(buf);
   (void) close (fd_from);
   if (fd_to > 0)
     (void) close (fd_to);

   return status;
}

static int enforce_path_uniqueness (char **path)
{
   char timestr[] = "yyyymmddThhmmssZ";
   size_t len = strlen(timestr);
   struct stat st = {0};
   struct tm tm = {0};
   char *newpath = NULL;
   time_t now;
   size_t n;
   int status;

   if (path == NULL)
     return -1;

   if (0 != stat (*path, &st))
     return 0;

   time(&now);
   gmtime_r (&now, &tm);
   if (len != (n = strftime (timestr, sizeof(timestr), "%Y%m%dT%H%M%SZ", &tm)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected strftime return value = %ld (expected %ld)",
                     __func__, n, len);
        return -1;
     }

   n += strlen(*path) + 2;  /* two strings, plus '.' plus terminating null char */

   if (NULL == (newpath = (char *)MALLOC (n)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   status = snprintf (newpath, n, "%s.%s", *path, timestr);
   if ((status < 0) || ((size_t) status >= n))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: adding suffix %s to path %s", __func__, timestr, *path);
        FREE(newpath);
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "avoiding filename collision: newpath=%s", newpath);

   FREE(*path);
   *path = newpath;

   return 0;
}

int copy_file_to_dir (const char *path, const char *copydir, const char *basename)
{
   char *copypath = NULL;
   int status = -1;

   if (0 != ioclib_mkdir (copydir, 0))
     return -1;

   if (NULL == (copypath = ioclib_pathconcat (copydir, basename)))
     return -1;

   if (0 != enforce_path_uniqueness (&copypath))
     goto return_status;

   tell_vinfo (0, "copying %s %s", path, copypath);

   if (0 != copy_file (path, copypath))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: copying from %s to %s", __func__, path, copypath);
        goto return_status;
     }

   status = 0;
return_status:

   FREE(copypath);
   return status;
}

static int compute_file_md5sum (const char *path, char *md5sum, int md5sum_size)
{
   FILE *fp = NULL;
   uint8_t result[MD5_NUM_BYTES];
   int i;

   if (md5sum_size < 2*MD5_NUM_BYTES+1)
     return -1;

   if (NULL == (fp = fopen (path, "r")))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %s", __func__, path);
        return -1;
     }
   md5File (fp, result);
   (void) fclose (fp);

   /* md5sum is the hexadecimal string representation of the MD5 checksum */
   for (i = 0; i < MD5_NUM_BYTES; i++)
     {
        sprintf (&md5sum[2*i], "%02x", result[i]);
     }
   md5sum[2*MD5_NUM_BYTES] = 0;

   return 0;
}

static int write_md5sum_file (const char *path)
{
   FILE *fp = NULL;
   char *path_md5 = NULL;
   char md5[2*MD5_NUM_BYTES+1];
   int status = -1;

   if (Disable_MD5_Checksum_Files)
     return 0;

   if (NULL == (path_md5 = ioclib_strcat (path, ".md5")))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: ioclib_strcat failed", __func__);
        return status;
     }

   if (NULL == (fp = fopen (path_md5, "w")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: cannot open %s", __func__, path_md5);
        goto return_status;
     }

   if (0 != compute_file_md5sum (path, md5, sizeof(md5)))
     goto return_status;

   if (fprintf (fp, "%s\n", md5) < 0)
     {
        (void) fclose (fp);
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s", __func__, path_md5);
        goto return_status;
     }

   if (0 != fclose (fp))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: closing %s", __func__, path_md5);
        goto return_status;
     }

   status = 0;
return_status:
   ioclib_free (path_md5);
   return status;
}

static int register_with_symlink (const char *dir, const char *basename)
{
   char *archived_path = NULL;
   char *registry_dir = NULL;
   char *symlink_path = NULL;
   const char *root_path;
   int status = -1;

   /* NULL means don't perform archiving */
   if (NULL == (root_path = get_archive_root_dir ()))
     return 0;

   if (Perform_Archive_Registration == 0)
     return 0;

   if ((NULL == (archived_path = ioclib_pathconcat (dir, basename)))
       || (NULL == (registry_dir = ioclib_pathconcat (root_path, "registry/incoming")))
       || (NULL == (symlink_path = ioclib_pathconcat (registry_dir, basename)))
      )
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed", __func__);
        goto return_status;
     }

   if (0 != write_md5sum_file (archived_path))
     goto return_status;

   if (0 != ioclib_mkdir (registry_dir, 0))
     goto return_status;

   if (0 != symlink (archived_path, symlink_path))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: symlink failed: %s -> %s",
                     __func__, archived_path, symlink_path);
        goto return_status;
     }

   status = 0;
return_status:
   ioclib_free (archived_path);
   ioclib_free (registry_dir);
   ioclib_free (symlink_path);
   return status;
}

static int read_time_coverage_start (const char *path, double *timestamp_start)
{
   int ncid, status;
   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     return -1;
   status = TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, timestamp_start);
   (void) TIO_close (ncid);
   return status ? -1 : 0;
}

static int release_with_symlink (const char *dir, const char *basename)
{
   char *archived_path = NULL;
   char *symlink_path = NULL;
   char *symlink_dir = NULL;
   struct stat st = {0};
   int len, n, status = -1;
   double sec_since_epoch, sat_day;

   if (Public_Mirror_Root_Dir == NULL)
     return 0;

   /* Silent return when directory Public_Mirror_Root_Dir does not exist */
   if ((0 != stat (Public_Mirror_Root_Dir, &st))
       || (0 == S_ISDIR(st.st_mode)))
     return 0;

   if (NULL == (archived_path = ioclib_pathconcat (dir, basename)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed", __func__);
        goto return_status;
     }

   if (0 != read_time_coverage_start (archived_path, &sec_since_epoch))
     goto return_status;

   if (0 != tio_time_sat_local_day_number (sec_since_epoch, &sat_day))
     goto return_status;

   len = strlen(Public_Mirror_Root_Dir) + 11;

   if (NULL == (symlink_dir = (char *)MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   n = snprintf (symlink_dir, len, "%s/D%05d/L0", Public_Mirror_Root_Dir, (int) sat_day);
   if ((n < 0) || (n >= len))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: error generating path to public mirror copy", __func__);
        return -1;
     }

   if (0 != ioclib_mkdir (symlink_dir, 0))
     goto return_status;

   if (NULL == (symlink_path = ioclib_pathconcat (symlink_dir, basename)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed", __func__);
        goto return_status;
     }

   if (0 != symlink (archived_path, symlink_path))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: symlink failed: %s -> %s",
                     __func__, archived_path, symlink_path);
        goto return_status;
     }

   status = 0;
return_status:
   ioclib_free (archived_path);
   ioclib_free (symlink_dir);
   ioclib_free (symlink_path);
   return status;
}

int copy_hidden (const char *dirname, const char *basename, const char *copydir)
{
   char hidden_basename[MAX_BASENAME_SIZE];
   char *oldpath = NULL;
   char *newpath = NULL;
   int status = -1;

   /* copydir == NULL means "don't copy" */
   if (copydir == NULL)
     return 0;

   if (-1 == make_hidden_basename (basename, hidden_basename, MAX_BASENAME_SIZE))
     return -1;

   if ((NULL == (oldpath = ioclib_pathconcat (dirname, hidden_basename)))
       || (NULL == (newpath = ioclib_pathconcat (dirname, basename))))
     goto return_status;

   if (0 != copy_file_to_dir (oldpath, copydir, basename))
     goto return_status;

   if (0 != register_with_symlink (copydir, basename))
     goto return_status;

   if (0 != release_with_symlink (copydir, basename))
     goto return_status;

   status = 0;
return_status:

   FREE(oldpath);
   FREE(newpath);
   return status;
}

int rename_hidden (const char *dirname, const char *basename)
{
   char hidden_basename[MAX_BASENAME_SIZE];
   char *oldpath = NULL;
   char *newpath = NULL;
   int status = -1;

   if (-1 == make_hidden_basename (basename, hidden_basename, MAX_BASENAME_SIZE))
     return -1;

   if ((NULL == (oldpath = ioclib_pathconcat (dirname, hidden_basename)))
       || (NULL == (newpath = ioclib_pathconcat (dirname, basename))))
     goto return_status;

   if (-1 == ioclib_rename (oldpath, newpath))
     goto return_status;

   status = 0;
return_status:

   FREE(oldpath);
   FREE(newpath);
   return status;
}

int remove_hidden (const char *dirname, const char *basename)
{
   char hidden_basename[MAX_BASENAME_SIZE];
   char *path = NULL;
   int status;

   if (-1 == make_hidden_basename (basename, hidden_basename, MAX_BASENAME_SIZE))
     return -1;

   if (NULL == (path = ioclib_pathconcat (dirname, hidden_basename)))
     return -1;

   status = ioclib_unlink (path);
   FREE(path);
   return status;
}

int annotate_var (int grp, int varid, const char *descr, const char *units)
{
   int len;

   if (descr)
     {
        len = strlen(descr);
        if (-1 == TIO_put_att (grp, varid, "comment", NC_CHAR, len, descr))
          return -1;
     }

   if (units)
     {
        len = strlen(units);
        if (-1 == TIO_put_att (grp, varid, "units", NC_CHAR, len, units))
          return -1;
     }

   return 0;
}

int record_source_file_contribution (int grp, const char *source_file, int num_records)
{
   TIO_Meta_Type *meta = NULL;
   const char *varname = "level0_source";
   char *basename = NULL;
   char buf[128];
   int bufsize = sizeof(buf);
   int n, status = -1;

   if (NULL == (basename = ioclib_basename (source_file)))
     return -1;

   n = snprintf (buf, bufsize, "%s|%d", basename, num_records);
   if ((n < 0) || (n >= bufsize))
     {
        tell_vwarn (0, "%s: updating metadata variable %s", __func__, varname);
        return -1;
     }

   if (NULL == (meta = tio_meta_open ()))
     return -1;

   if (0 != tio_meta_ncinit (meta, grp, varname, TIO_META_TYPE_STRING))
     goto return_status;

   if (0 != tio_meta_append_string (meta, varname, buf))
     goto return_status;

   if (0 != tio_meta_write_ncattr (meta, grp))
     goto return_status;

   status = 0;
return_status:
   tio_meta_close (meta);
   return status;
}

int main (int argc, char **argv)
{
   const char *appname = "L0_format";
   const char *param_file = "l0_format.cfg";
   Process_Method_Table_Type *tbl = Method_Table;
   Processed_File_Log_Type incoming_log_info = {0};
   Control_Type ctrl = {0};
   config_t cfg = {0};
   TPInfo_Type *tp = NULL;
   int verbose = 0;
   int cache_method = EXPREC_CACHE_DISK;
   const char *cache_dir_pattern = NULL;
   int status = EXIT_FAILURE;
   double start_time = 0.0;
   static struct option long_options[] =
     {
        {"help",     no_argument,       0, 'h'},
        {"archive",  required_argument, 0, 'a'},
        {"mirror",   required_argument, 0, 'm'},
        {"cache",    required_argument, 0, 'c'},
        {"logdir",   required_argument, 0, 'L'},
        {"tstart",   required_argument, 0, 't'},
        {"empty",    no_argument,       0, 'e'},
        {"register", no_argument,       0, 'r'},
        {"nochecksum", no_argument,     0, 'n'},
        {"verbose",  no_argument,       0, 'v'},
        {"Version",  required_argument, 0, 'V'},
        {0,0,0,0}
     };

   memset ((char *)&ctrl, 0, sizeof ctrl);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "ha:m:nc:eL:rvV:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d?\n", c);
             goto return_status;
             break;
           case 'h':
             usage();
             break;
           case 'a':
             set_archive_root_dir (optarg);
             break;
           case 'm':
             set_public_mirror_root_dir (optarg);
             break;
           case 'n':
             Disable_MD5_Checksum_Files++;
             break;
           case 'L':
             incoming_log_info.logdir = optarg;
             ctrl.log_incoming = &incoming_log_info;
             break;
           case 'c':
             cache_dir_pattern = optarg;
             cache_method = EXPREC_CACHE_MEM;
             break;
           case 'e':
             ctrl.exit_on_emptydir = 1;
             break;
           case 'r':
             Perform_Archive_Registration++;
             break;
           case 't':
             if (1 != sscanf (optarg, "%le", &start_time))
               {
                  fprintf (stderr, "*** Error parsing start time: %s\n", optarg);
                  goto return_status;
               }
             break;
           case 'v':
             verbose++;
             break;
           case 'V':
             if (1 != sscanf (optarg, "%d", &Processing_Version))
               {
                  fprintf (stderr, "*** Error parsing processing version: %s\n", optarg);
                  goto return_status;
               }
             break;
          }
     }

   if (optind < argc)
     {
        param_file = argv[optind++];
     }

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   if (0 != access (param_file, F_OK | R_OK))
     {
        fprintf (stderr, "*** Cannot access param file: %s\n", param_file);
        usage();
     }

   tell_open (appname, -1, 0);
   tell_set_log_level (TELL_MSGTYPE_INFO, verbose);
   config_init (&cfg);

   Have_Epoch = 0;

   if (-1 == parse_param_file (&cfg, param_file, &ctrl))
     goto return_status;

   /* command-line start_time overrides the control file */
   if (start_time != 0.0)
     {
        ctrl.start_time = start_time;
     }

   if (NULL == (tp = tpinfo_init (ctrl.tpinfo_file)))
     goto return_status;

   /* must precede init_methods_table call */
   set_exprec_cache_method (cache_method);

   if (-1 == init_methods_table (tbl, &cfg))
     goto return_status;

   if (0 != open_processed_file_log (ctrl.log_incoming))
     goto return_status;

   tell_vinfo (0, "started (pid=%d)", getpid());

   catch_sigterm ();
   catch_sigint ();

   if (cache_dir_pattern)
     status = process_cache_dirs (tbl, tp, &ctrl, cache_dir_pattern);
   else
     status = process_live_stream (tbl, tp, &ctrl);

   tell_vinfo (0, "last packet time = %f", Last_Packet_Time);

   delete_methods_table (tbl);

   status = (status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

return_status:
   tell_vinfo (0, "exiting: status = %d (pid=%d)", status, getpid());
   finalize_processed_file_log (ctrl.log_incoming);
   free_control_type_fields (&ctrl);
   tpinfo_free (tp);
   config_destroy (&cfg);
   tell_close();
   return status;
}
