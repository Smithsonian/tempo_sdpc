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
#include <sys/types.h>
#include <sys/stat.h>
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

#include "l0_format.h"

static int Have_Epoch;
static int Perform_Archive_Registration;
static const char *Archive_Root_Dir = NULL;
static double Cache_Start_Time = 0.0;

static Process_Method_Type *Exprec_Process_Method;

typedef struct
{
   int filetype;
   Process_Method_Type *method;
   Process_Method_Type *(*init)(config_t *);
   Process_Method_Callback_Function *post_process_callback;
}
Process_Method_Table_Type;

typedef struct IRU_Gap_Type IRU_Gap_Type;
struct IRU_Gap_Type
{
   IRU_Gap_Type *next;
   double tbeg;
   double tend;
};

typedef struct
{
   double max_iru_knowledge_gap_duration;
   double latest_iru_timestamp_seen;
   double latest_radiance_timestamp_seen;
   double latest_iru_only_interval_end_time;
   char *dir;
   IRU_Gap_Type *gap_list;
}
IRU_Interval_Type;

typedef struct
{
   const char *input_filename_glob_pattern;
   char *incoming_dir;
   char *tpinfo_file;
   double monitor_wait_secs;
   double cache_flush_idle_wait_secs;
   double cache_flush_exprec_wait_secs;
   int exit_on_emptydir;
   IRU_Interval_Type iru_interval;
}
Control_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: L0_format [options] [config-file]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -e | --empty             Exit when the input directory is empty\n");
   fprintf (stderr, "   -a | --archive DIR       Archive files in directory DIR\n");
   fprintf (stderr, "   -c | --cache DIR         Process cached directories matching regex DIR\n");
   fprintf (stderr, "   -t | --tstart SEC        Process cache files newer than SEC since the TEMPO epoch\n");
   fprintf (stderr, "   -v | --verbose           Increase verbosity (-vv is more verbose)\n");
   exit (EXIT_SUCCESS);
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

static void iru_gap_free_list (IRU_Gap_Type *glst)
{
   while (glst)
     {
        IRU_Gap_Type *next = glst->next;
        FREE(glst);
        glst = next;
     }
}

/* add a new interval to a sorted list */
static int iru_gap_append (IRU_Gap_Type **glst, double tbeg, double tend)
{
   IRU_Gap_Type *g = NULL;
   IRU_Gap_Type *p;

   if (tbeg >= tend)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: invalid IRU gap: tbeg=%f tend=%f",
                     __func__, tbeg, tend);
        return -1;
     }

   if (NULL == (g = (IRU_Gap_Type *)MALLOC (sizeof *g)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   g->next = NULL;
   g->tbeg = tbeg;
   g->tend = tend;

   /* If this is the first entry, we're done */
   if (*glst == NULL)
     {
        *glst = g;
        return 0;
     }

   /* Append new entry */
   for (p = *glst; p != NULL; p = p->next)
     {
        if (p->next != NULL)
          continue;

        if (p->tend > g->tbeg)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: IRU gap out of order: list tend=%f gap tbeg=%f",
                          __func__, p->tend, g->tbeg);
             FREE(g);
             return -1;
          }

        p->next = g;
        break;
     }

   return 0;
}

static void free_control_type_fields (Control_Type *ctrl)
{
   FREE(ctrl->incoming_dir);
   FREE(ctrl->tpinfo_file);
   FREE(ctrl->iru_interval.dir);
   iru_gap_free_list (ctrl->iru_interval.gap_list);
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
       || (CONFIG_TRUE != config_setting_lookup_float (s, "cache_flush_idle_wait_secs", &ctrl->cache_flush_idle_wait_secs))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "cache_flush_exprec_wait_secs", &ctrl->cache_flush_exprec_wait_secs))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'main' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

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

   ctrl->iru_interval.gap_list = NULL;

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

   n = __tio_filename_string (str, strsize, tstart, "inr", 0, 0);
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

static int exprec_post_process_callback (Process_Method_Type *pmt, void *client_data)
{
   IRU_Interval_Type *iru_interval = (IRU_Interval_Type *)client_data;
   double t_max = iru_interval->max_iru_knowledge_gap_duration;
   double t_rad_prev = iru_interval->latest_radiance_timestamp_seen;
   double t_rad;

   if (0 != pmt->pmt_query_latest_timestamp (pmt, IOCSDPC_EXPREC_TYPE_RADIANCE, &t_rad))
     return -1;

   /* record any large gaps between radiance exposure records */
   if ((t_rad_prev > 0) && (t_rad - t_rad_prev > t_max))
     {
        if (0 != iru_gap_append (&iru_interval->gap_list, t_rad_prev, t_rad))
          return -1;
     }

   iru_interval->latest_radiance_timestamp_seen = t_rad;

   return 0;
}

static int resolve_iru_gaps (IRU_Interval_Type *iru_interval)
{
   IRU_Gap_Type *head = iru_interval->gap_list;
   double t_max = iru_interval->max_iru_knowledge_gap_duration;

   /* If there are no gaps, do nothing */
   if (head == NULL)
     return 0;

   /* Try to resolve previously accumulated gaps */
   while (head != NULL)
     {
        IRU_Gap_Type *next;
        double dt_gap, dt;
        int i, n;

        /* Is the oldest gap in radiance data covered by the available IRU data? */
        if (iru_interval->latest_iru_timestamp_seen < head->tend)
          return 0;

        /* Write out time intervals for IRU-only files to span this gap */
        dt_gap = head->tend - head->tbeg;
        n = (dt_gap < t_max) ? 1 : (dt_gap / t_max);
        dt = dt_gap / n;

        for (i = 0; i < n; i++)
          {
             double tbeg = head->tbeg + i*dt;
             double tend = tbeg + dt;
             if (0 != write_iru_only_interval (iru_interval->dir, tbeg, tend))
               return -1;
             iru_interval->latest_iru_only_interval_end_time = tend;
          }

        /* Gap issue is resolved -- no need to track it further */
        next = head->next;
        FREE(head);
        head = next;
        iru_interval->gap_list = head;
     }

   return 0;
}

static int iru_post_process_callback (Process_Method_Type *pmt, void *client_data)
{
   IRU_Interval_Type *iru_interval = (IRU_Interval_Type *)client_data;
   double t_max = iru_interval->max_iru_knowledge_gap_duration;
   double t_iru, t_last, t_end, t_only, t_rad;

   if (0 != pmt->pmt_query_latest_timestamp (pmt, 0, &t_iru))
     return -1;
   iru_interval->latest_iru_timestamp_seen = t_iru;

   if (0 != resolve_iru_gaps (iru_interval))
     return -1;

   t_only = iru_interval->latest_iru_only_interval_end_time;
   t_rad = iru_interval->latest_radiance_timestamp_seen;

   /* Initialization: If we've never sent any IRU data to the INR subsystem,
    * we'll assume the current IRU knowledge gap starts with the latest
    * IRU timestamp.  Avoid gaps by padding the IRU coverage in every granule.
    */
   if (t_only <= 0.0 && t_rad <= 0.0)
     {
        iru_interval->latest_iru_only_interval_end_time = t_iru;
        return 0;
     }

   /* What's the latest IRU timestamp headed for the INR subsystem? */
   t_last = (t_rad > t_only) ? t_rad : t_only;

   /* Does the INR subsystem need an IRU update? */
   if (t_iru - t_last < t_max)
     return 0;

   t_end = t_last + t_max;

   return write_iru_only_interval (iru_interval->dir, t_last, t_end);
}

static int process_file (const Process_Method_Table_Type *tbl, const TPInfo_Type *tpinfo,
                         Control_Type *ctrl, const char *file)
{
   Process_Method_Type *pmt;
   IOCSDPC_Common_Header_Type chdr = {0};
   int fd, filetype;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;
   (void) ioclib_fd_close (fd);

   filetype = chdr.filetype;

   /* If we're processing the science data cache, we may want to skip files
    * prior to some user-specified time */
   if (Cache_Start_Time > 0)
     {
        if (0 != verify_epoch (chdr.epoch))
          return -1;
        if (chdr.last_packet_time < Cache_Start_Time)
          return 0;
     }

   if (NULL == (pmt = find_process_method (tbl, filetype)))
     return -1;

   return pmt->pmt_process (pmt, tpinfo, file, &ctrl->iru_interval);
}

static int process_dir_files (const Process_Method_Table_Type *tbl,
                              const TPInfo_Type *tpinfo, Control_Type *ctrl,
                              char **file_list, size_t num_files)
{
   size_t i;

   for (i = 0; i < num_files; i++)
     {
        char *file = file_list[i];
        if (0 == process_file (tbl, tpinfo, ctrl, file))
          {
             tell_vinfo (1, "processed: %s", file);
             /* If processing involved a rename, then deletion
              * will be handled elsewhere. */
             if (ioclib_isfile (file, NULL))
               {
                  if (0 != ioclib_unlink (file))
                    return -1;
               }
          }
        else
          {
             tell_vinfo (0, "%s: bad file: %s", __func__, file);
             if (0 != ioclib_rename_to_bad_file (file))
               {
                  tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_rename_to_bad_file failed, file=%s",
                               __func__, file ? file : "(null)");
                  return -1;
               }
          }
     }

   return 0;
}

static int flush_caches (const Process_Method_Table_Type *tbl,
                         const TPInfo_Type *tpinfo)
{
   Process_Method_Type *pmt;
   int num_failed = 0;

   for (; tbl->init != NULL; tbl++)
     {
        pmt = tbl->method;
        if (pmt->pmt_flush_cache)
          {
             if (0 != pmt->pmt_flush_cache (pmt, tpinfo))
               num_failed++;
          }
     }

   return num_failed;
}

static int maybe_flush_exprec_cache (const TPInfo_Type *tpinfo, Control_Type *ctrl)
{
   Process_Method_Type *pmt = Exprec_Process_Method;
   time_t when_last_exprec_cached, age_secs;

   /* If we're not processing exposure records, there's nothing more to do here */
   if (pmt == NULL)
     return 0;

   if (0 != pmt->pmt_query_when_last_rec_cached (pmt, &when_last_exprec_cached))
     return -1;

   /* If the cache is empty, we're done */
   if (when_last_exprec_cached <= 0)
     return 0;

   age_secs = time(NULL) - when_last_exprec_cached;

   if (age_secs < ctrl->cache_flush_exprec_wait_secs)
     return 0;

   tell_vinfo (0, "flush exprec cache (newest cached data is %ld sec old)", age_secs);

   return pmt->pmt_flush_cache (pmt, tpinfo);
}

#define PROCESS_METHOD(filetype,init,callback) {filetype,NULL,init,callback}
#define PROCESS_METHODS_TABLE_END {-1,NULL,NULL,NULL}

static Process_Method_Table_Type Method_Table[] =
{
   PROCESS_METHOD(IOCSDPC_FILETYPE_EXPREC, init_exprec_method, exprec_post_process_callback),
   PROCESS_METHOD(IOCSDPC_FILETYPE_TPSEC, init_tpsec_method, NULL),
   PROCESS_METHOD(IOCSDPC_FILETYPE_IRU, init_iru_method, iru_post_process_callback),
   PROCESS_METHOD(IOCSDPC_FILETYPE_SMC, init_smc_method, NULL),
   /* PROCESS_METHOD(IOCSDPC_FILETYPE_TLMRAW, init_tlmraw_method), */
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
   int may_have_cached_files, status = -1;
   double time_since_last_file;

   pattern = ioclib_pathconcat (ctrl->incoming_dir,
                                ctrl->input_filename_glob_pattern);

   tell_vinfo (0, "processing %s", pattern);

   if (NULL == pattern)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed",
                    __func__);
        return -1;
     }

   time_since_last_file = 0.0;
   may_have_cached_files = 0;

   while (0 == caught_signal())
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

        if (-1 == process_dir_files (tbl, tpinfo, ctrl, gt->files, gt->num_files))
          goto return_status;

        (void) ioclib_sleep (ctrl->monitor_wait_secs);

        /* If we haven't seen an exprec in a while, maybe it's time to flush the exprec cache */
        if (0 != maybe_flush_exprec_cache (tpinfo, ctrl))
          goto return_status;

        if (gt->num_files)
          {
             time_since_last_file = 0.0;
             may_have_cached_files = 1;
          }
        else
          {
             time_since_last_file += ctrl->monitor_wait_secs;
          }

        /* If we haven't seen any files in a while, all the caches may need flushing. */
        if ((may_have_cached_files != 0) &&
            (time_since_last_file > ctrl->cache_flush_idle_wait_secs))
          {
             tell_vinfo (0, "flush caches (%g sec since last file)", time_since_last_file);
             if (0 != flush_caches (tbl, tpinfo))
               goto return_status;
             may_have_cached_files = 0;
          }
     }

   if (caught_signal())
     log_caught_signal();

   tell_vinfo (0, "flush caches on exit");
   if (0 != flush_caches (tbl, tpinfo))
     goto return_status;

   status = 0;
return_status:
   ioclib_free (pattern);
   ioclib_glob_free (gt);

   return status;
}

static int process_cache_dirs (Process_Method_Table_Type *tbl,
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
        if (Cache_Start_Time > 0)
          {
             tell_vinfo (0, "start time = %f sec since the epoch", Cache_Start_Time);
          }

        for (k = 0; k < num_files; k++)
          {
             if (NULL == (path = ioclib_pathconcat (dir, files[k])))
               goto return_status;

             if (caught_signal())
               {
                  log_caught_signal();
                  goto return_status;
               }

             if (-1 == process_file (tbl, tpinfo, ctrl, path))
               goto return_status;

             ioclib_free (path);
             path = NULL;
          }

        tell_vinfo (0, "end processing %ld files from directory %s", num_files, dir);

        ioclib_string_array_free (files, num_files);
        files = NULL;
        num_files = 0;
     }

   tell_vinfo (0, "flush caches on exit");
   if (0 != flush_caches (tbl, tpinfo))
     goto return_status;

   status = 0;
return_status:
   ioclib_glob_free (gt);
   ioclib_string_array_free (files, num_files);
   ioclib_free (path);

   return status;
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
                              int processing_version, const char *suffix)
{
   char buf[MAX_PATHLEN];
   size_t bufsize = sizeof(buf);
   const char *root_path;
   char *path = NULL;
   int sat_day;
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

   /* e.g. ${SDPC_ARCHIVE_DIR}/L0/${version}/ddddd/${file_type}
    *   or ${SDPC_ARCHIVE_DIR}/L0/${version}/ddddd/${file_type}/scan
    */
   if (scan_num < 0)
     {
        n = snprintf (buf, bufsize, "%s/L0/%d/%d/%s",
                      root_path, processing_version, sat_day, suffix);
     }
   else
     {
        n = snprintf (buf, bufsize, "%s/L0/%d/%d/%s/%d",
                      root_path, processing_version, sat_day, suffix, scan_num);
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

int write_attr_global_product_type (int ncid, const char *product_type)
{
   int len;
   if ((product_type == NULL) || (*product_type == 0))
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: null product_type", __func__);
        return -1;
     }
   len = strlen(product_type);
   return TIO_put_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, len, product_type);
}

int write_attr_global_timestamp (int ncid, const char *tstamp_name,
                                 double tstamp_value)
{
   return TIO_write_timestamp (ncid, NC_GLOBAL, tstamp_name, tstamp_value);
}

int write_std_global_metadata (int ncid, const IOCSDPC_Common_Header_Type *chdr)
{
   (void) chdr;
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

   if ((NULL == (path = ioclib_pathconcat (dirname, hidden_basename)))
       || (-1 == TIO_create (path, NC_NETCDF4, ncid))
      )
     {
        status = -1;
     }

   FREE(path);
   return status;
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

   bufsize = st.st_blksize;

   if ((fd_from = open (from, O_RDONLY)) < 0)
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: opening %s", __func__, from);
        return -1;
     }

   if ((fd_to = open (to, O_WRONLY | O_CREAT | O_EXCL, mode_create)) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: opening %s (%s)", __func__, to, strerror (errno));
        goto return_status;
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

static int perform_copy (const char *path, const char *copydir, const char *basename)
{
   char *copypath = NULL;
   int status = -1;

   /* NULL means "don't copy" */
   if (copydir == NULL)
     return 0;

   if (0 != ioclib_mkdir (copydir, 0))
     return -1;

   if (NULL == (copypath = ioclib_pathconcat (copydir, basename)))
     return -1;

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

static int register_with_symlink (const char *copydir, const char *basename)
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

   if ((NULL == (archived_path = ioclib_pathconcat (copydir, basename)))
       || (NULL == (registry_dir = ioclib_pathconcat (root_path, "registry/incoming")))
       || (NULL == (symlink_path = ioclib_pathconcat (registry_dir, basename)))
      )
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed", __func__);
        goto return_status;
     }

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

/* Close hidden file $dirname/.${basename}
 * and rename to $dirname/$basename.
 * Optionally, if copydir != NULL, put a copy in $copydir/$basename
 * before performing the rename.
 */
int close_hidden (int ncid, const char *dirname, const char *basename,
                  const char *copydir)
{
   char hidden_basename[MAX_BASENAME_SIZE];
   char *oldpath = NULL;
   char *newpath = NULL;
   int status = -1;

   if (-1 == TIO_close (ncid))
     return -1;

   if (-1 == make_hidden_basename (basename, hidden_basename, MAX_BASENAME_SIZE))
     return -1;

   if ((NULL == (oldpath = ioclib_pathconcat (dirname, hidden_basename)))
       || (NULL == (newpath = ioclib_pathconcat (dirname, basename))))
     goto return_status;

   if (0 != perform_copy (oldpath, copydir, basename))
     goto return_status;

   if (0 != register_with_symlink (copydir, basename))
     goto return_status;

   if (-1 == ioclib_rename (oldpath, newpath))
     goto return_status;

   status = 0;
return_status:

   FREE(oldpath);
   FREE(newpath);
   return status;
}

int remove_file (const char *dirname, const char *basename)
{
   char *path = NULL;
   int status;

   if (NULL == (path = ioclib_pathconcat (dirname, basename)))
     return -1;
   status = ioclib_unlink (path);
   ioclib_free (path);
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

int main (int argc, char **argv)
{
   const char *appname = "L0_format";
   const char *param_file = "l0_format.cfg";
   Process_Method_Table_Type *tbl = Method_Table;
   Control_Type ctrl = {0};
   config_t cfg = {0};
   TPInfo_Type *tp = NULL;
   int verbose = 0;
   int cache_method = EXPREC_CACHE_DISK;
   const char *cache_dir_pattern = NULL;
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"help",     no_argument,       0, 'h'},
        {"archive",  required_argument, 0, 'a'},
        {"cache",    required_argument, 0, 'c'},
        {"tstart",   required_argument, 0, 't'},
        {"empty",    no_argument,       0, 'e'},
        {"register", no_argument,       0, 'r'},
        {"verbose",  no_argument,       0, 'v'},
        {0,0,0,0}
     };

   memset ((char *)&ctrl, 0, sizeof ctrl);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "ha:c:erv", long_options, &option_index);
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
             if (1 != sscanf (optarg, "%le", &Cache_Start_Time))
               {
                  fprintf (stderr, "*** Error parsing start time: %s\n", optarg);
                  goto return_status;
               }
             break;
           case 'v':
             verbose++;
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

   tell_open (appname, -1, 0);
   tell_set_log_level (TELL_MSGTYPE_INFO, verbose);
   config_init (&cfg);

   Have_Epoch = 0;

   if (-1 == parse_param_file (&cfg, param_file, &ctrl))
     goto return_status;

   if (NULL == (tp = tpinfo_init (ctrl.tpinfo_file)))
     goto return_status;

   /* must precede init_methods_table call */
   set_exprec_cache_method (cache_method);

   if (-1 == init_methods_table (tbl, &cfg))
     goto return_status;

   tell_vinfo (0, "started (pid=%d)", getpid());

   catch_sigterm ();
   catch_sigint ();

   if (cache_dir_pattern)
     status = process_cache_dirs (tbl, tp, &ctrl, cache_dir_pattern);
   else
     status = process_live_stream (tbl, tp, &ctrl);

   delete_methods_table (tbl);

   status = (status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

return_status:
   tell_vinfo (0, "exiting: status = %d (pid=%d)", status, getpid());
   free_control_type_fields (&ctrl);
   tpinfo_free (tp);
   config_destroy (&cfg);
   tell_close();
   return status;
}
