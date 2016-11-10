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
#include <unistd.h>
#include <signal.h>

#define __USE_XOPEN  /* for strptime */
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "l0_format.h"
#include "daemon.h"

enum
{
   NODELIM_TIMESTAMP = 0,
     DELIM_TIMESTAMP = 1
};

typedef struct
{
   int filetype;
   Process_Method_Type *method;
   Process_Method_Type *(*init)(config_t *);
}
Process_Method_Table_Type;

typedef struct
{
   const char *incoming_dir;
   const char *tpinfo_file;
   const char *input_filename_glob_pattern;
   const char *daemon_logfile_path;
   double monitor_wait_secs;
   int exit_on_emptydir;
   int daemon;
}
Control_Type;

static int parse_param_file (config_t *cfg, const char *cfg_file,
                             Control_Type *ctrl)
{
   config_setting_t *s;

   if (0 == config_read_file (cfg, cfg_file))
    {
       tell_verror (TELL_INVALID_PARM_ERROR,
                    "%s: Reading %s: %s:%d - %s",
                    __func__, cfg_file, config_error_file(cfg),
                    config_error_line(cfg), config_error_text(cfg));
       return -1;
    }

   if (NULL == (s = config_lookup (cfg, "main")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'main' in param file: %s",
                     __func__, cfg_file);
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "incoming_dir", &ctrl->incoming_dir))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "input_filename_glob_pattern", &ctrl->input_filename_glob_pattern))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "tpinfo_file", &ctrl->tpinfo_file))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "daemon_logfile_path", &ctrl->daemon_logfile_path))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "monitor_wait_secs", &ctrl->monitor_wait_secs))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'main' parameters in param file: %s",
                     __func__, cfg_file);
        return -1;
     }

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

static int process_file (const Process_Method_Table_Type *tbl,
                         const TPInfo_Type *tpinfo,
                         const char *file)
{
   Process_Method_Type *pmt;
   IOCSDPC_Common_Header_Type chdr;
   int fd, status;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;

   if (NULL == (pmt = find_process_method (tbl, chdr.filetype)))
     return -1;

   status = pmt->process (pmt, tpinfo, file, fd, &chdr);

   (void) ioclib_fd_close (fd);

   return status;
}

static int process_dir_files (const Process_Method_Table_Type *tbl,
                              const TPInfo_Type *tpinfo,
                              char **file_list, size_t num_files)
{
   size_t i;

   for (i = 0; i < num_files; i++)
     {
        char *file = file_list[i];
        if (0 == process_file (tbl, tpinfo, file))
          {
             tell_vlog (TELL_MSGTYPE_INFO, 0, "processed: %s", file);
             if (0 != ioclib_unlink (file))
               return -1;
          }
        else
          {
             tell_vlog (TELL_MSGTYPE_INFO, 0, "bad file: %s", file);
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

#define PROCESS_METHOD(filetype,init) {filetype,NULL,init}
#define PROCESS_METHODS_TABLE_END {-1,NULL,NULL}

static Process_Method_Table_Type Method_Table[] =
{
   PROCESS_METHOD(IOCSDPC_FILETYPE_EXPREC, init_exprec_method),
   PROCESS_METHOD(IOCSDPC_FILETYPE_TPSEC, init_tpsec_method),
   /* PROCESS_METHOD(IOCSDPC_FILETYPE_IRU, init_iru_method), */
   /* PROCESS_METHOD(IOCSDPC_FILETYPE_SMC, init_smc_method), */
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
        if ((m != NULL) && (m->delete != NULL))
          {
             m->delete (m);
          }
     }
}

static volatile int Sighup_Received;
static void sighup_handler (int sig)
{
   (void) sig;
   Sighup_Received++;
}
static void catch_sighup (void)
{
   struct sigaction new_action;
   new_action.sa_handler = sighup_handler;
   sigemptyset (&new_action.sa_mask);
   new_action.sa_flags = 0;
   sigaction (SIGHUP, &new_action, NULL);
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

   if (Sighup_Received)
     signame = "SIGHUP";
   else if (Sigint_Received)
     signame = "SIGINT";
   else
     signame = "unknown";

   tell_vlog (TELL_MSGTYPE_INFO, 0, "caught signal: %s (pid=%d)",
              signame, getpid());
}

#define NO_CAUGHT_SIGNALS ((Sighup_Received == 0) && (Sigint_Received == 0))

static int monitor_dir (Process_Method_Table_Type *tbl,
                        const TPInfo_Type *tpinfo, Control_Type *ctrl)
{
   IOCLib_Glob_Type *gt = NULL;
   char *pattern = NULL;
   int status = -1;

   pattern = ioclib_pathconcat (ctrl->incoming_dir,
                                ctrl->input_filename_glob_pattern);

   tell_vlog (TELL_MSGTYPE_INFO, 0, "processing %s", pattern);

   if (NULL == pattern)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_pathconcat failed",
                    __func__);
        return -1;
     }

   while (NO_CAUGHT_SIGNALS)
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

        if (-1 == process_dir_files (tbl, tpinfo, gt->files, gt->num_files))
          goto return_status;

        (void) ioclib_sleep (ctrl->monitor_wait_secs);
     }

   status = 0;
return_status:
   log_caught_signal();
   ioclib_free (pattern);
   ioclib_glob_free (gt);

   return status;
}

static time_t Epoch_Time_t;

static int init_epoch_time_t (void)
{
   struct tm tm;

   memset ((char *)&tm, 0, sizeof(struct tm));
   if (NULL == strptime (TIO_TIME_REFERENCE_STRING,
                         TIO_DELIM_TIMESTAMP_FORMAT, &tm))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: strptime failed: s=%s",
                     __func__, TIO_TIME_REFERENCE_STRING);
        return -1;
     }

   Epoch_Time_t = timegm (&tm);

   return 0;
}

static int mktimestamp_str (double sec_since_epoch, int delim, char *buf, int bufsize)
{
   struct tm tm;
   time_t tt;
   int status;

   tt = Epoch_Time_t + sec_since_epoch;

   memset ((char *)&tm, 0, sizeof(struct tm));

   if (NULL == gmtime_r (&tt, &tm))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: gmtime_r failed: tt=%ld",
                     __func__, tt);
        return -1;
     }

   if (delim == 0)
     status = strftime (buf, bufsize, TIO_NODELIM_TIMESTAMP_FORMAT, &tm);
   else
     status = strftime (buf, bufsize, TIO_DELIM_TIMESTAMP_FORMAT, &tm);

   if (0 == status)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: strftime failed, tt=%ld",
                     __func__, tt);
        return -1;
     }

   return 0;
}

int make_level0_basename (double sec_since_epoch, int processing_version,
                          const char *suffix, char *buf, int bufsize)
{
   char tstr[MAX_ISOTIME_LEN];
   int n;

   if (-1 == mktimestamp_str (sec_since_epoch, NODELIM_TIMESTAMP, tstr, MAX_ISOTIME_LEN))
     return -1;

   n = snprintf (buf, bufsize, "tempo_%s_v%d_%s.nc",
                 tstr, processing_version, suffix);
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
   char buf[MAX_ISOTIME_LEN];
   int len;

   if (-1 == mktimestamp_str (tstamp_value, DELIM_TIMESTAMP, buf, sizeof(buf)))
     return -1;

   len = strlen (buf);
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, tstamp_name, NC_CHAR, len, buf))
     return -1;

   return 0;
}

static int write_std_global_metadata (int ncid)
{
   const char *time_ref = TIO_TIME_REFERENCE_STRING;
   int len = strlen (time_ref);

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_reference", NC_CHAR, len, time_ref))
     return -1;

   return 0;
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
       || (-1 == write_std_global_metadata (*ncid))
      )
     {
        status = -1;
     }

   FREE(path);
   return status;
}

/* Close hidden file $dirname/.${basename} and
 * and rename to $dirname/$basename */
int close_hidden (int ncid, const char *dirname, const char *basename)
{
   char hidden_basename[MAX_BASENAME_SIZE];
   char *oldpath = NULL;
   char *newpath = NULL;
   int status = 0;

   if (-1 == TIO_close (ncid))
     return -1;

   if (-1 == make_hidden_basename (basename, hidden_basename, MAX_BASENAME_SIZE))
     return -1;

   if ((NULL == (oldpath = ioclib_pathconcat (dirname, hidden_basename)))
       || (NULL == (newpath = ioclib_pathconcat (dirname, basename)))
       || (-1 == ioclib_rename (oldpath, newpath)))
     {
        status = -1;
     }

   FREE(oldpath);
   FREE(newpath);
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
   Control_Type ctrl;
   config_t cfg;
   TPInfo_Type *tp = NULL;
   int status = EXIT_FAILURE;

   memset ((char *)&ctrl, 0, sizeof ctrl);
   if (0 == strcmp (argv[1], "--empty"))
     {
        ctrl.exit_on_emptydir = 1;
        argc--;
        argv++;
     }

   if (0 == strcmp (argv[1], "--daemon"))
     {
        ctrl.daemon = 1;
        argc--;
        argv++;
     }

   if (argc > 1)
     param_file = argv[1];

   tell_open (appname, -1, 0);

   if (-1 == init_epoch_time_t ())
     goto return_status;

   config_init (&cfg);

   if (-1 == parse_param_file (&cfg, param_file, &ctrl))
     goto return_status;

   if (NULL == (tp = tpinfo_init (ctrl.tpinfo_file)))
     goto return_status;

   if (-1 == init_methods_table (tbl, &cfg))
     goto return_status;

   if (ctrl.daemon)
     {
        /* daemonize calls tell_open for logfile */
        tell_close();
        if (0 != daemonize (appname, ctrl.daemon_logfile_path))
          goto return_status;
        tell_vlog (TELL_MSGTYPE_INFO, 0, "daemon started (pid=%d)", getpid());
     }
   catch_sighup ();
   catch_sigint ();

   status = monitor_dir (tbl, tp, &ctrl);
   delete_methods_table (tbl);

   status = (status == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

return_status:
   if (ctrl.daemon)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 0, "daemon exiting: status = %d (pid=%d)",
                   status, getpid());
     }
   tpinfo_free (tp);
   config_destroy (&cfg);
   tell_close();
   return status;
}
