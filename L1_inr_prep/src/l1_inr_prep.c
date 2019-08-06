/** @file l1_inr_prep.c
 *  @brief Main program
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <wordexp.h>

#include <libconfig.h>

#include <ioclib.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "row_select.h"
#include "radiance.h"
#include "ephem.h"

#define BASENAME_SIZE 256

typedef struct
{
   const char *file_glob_pattern;
   int num_pad;
   Row_Select_Type *rst;
}
Selection_Type;

typedef struct
{
   char *tmp_path;
   char *final_path;
   char *target_dir;
   int processing_version;
}
Rename_Path_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: L1_inr_prep [options] [FILE]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -b | --begin <start-time>  start time (sec since epoch)\n");
   fprintf (stderr, "   -e | --end <end-time>      stop time (sec since epoch)\n");
   fprintf (stderr, "   -E | --epoch SEC           epoch (UTC sec since Unix epoch, e.g. a time_t value)\n");
   fprintf (stderr, "   -p | --ephemeris FILE      ephemeris file\n");
   fprintf (stderr, "   -d | --delay SEC           delay start (to wait for all telemetry to arrive)\n");
   fprintf (stderr, "   -c | --config FILE         configuration file\n");
   fprintf (stderr, "   -v | --verbose lev         logging verbosity\n");
   exit (EXIT_SUCCESS);
}

static int read_config_file (const char *config_file, config_t *cfg)
{
   if (0 == config_read_file (cfg, config_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: Reading %s:%d - %s",
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

static int read_common_params (config_t *cfg, const char *setting_name,
                               Selection_Type *st)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, setting_name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing '%s' in param file: %s",
                     __func__, setting_name, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "file_glob", &st->file_glob_pattern))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num_pad", &st->num_pad)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading '%s' parameters in param file: %s",
                     __func__, setting_name, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int copy_iru (Radiance_Type *r, config_t *cfg,
                     double time_beg, double time_end, int pad_enable)
{
   Selection_Type iru = {0};
   int status = -1;

   if (0 != read_common_params (cfg, "iru_config", &iru))
     return -1;

   if (0 != row_select_scan (time_beg, time_end,
                             pad_enable ? iru.num_pad : 0,
                             iru.file_glob_pattern, &iru.rst))
     goto return_status;

   if (iru.rst)
     {
        if (0 != radiance_copy_iru (r, iru.rst))
          goto return_status;
     }

   status = 0;
return_status:
   row_select_free (iru.rst);
   return status;
}

static int copy_smc (Radiance_Type *r, config_t *cfg,
                     double time_beg, double time_end, int pad_enable)
{
   Selection_Type smc = {0};
   int status = -1;

   if (0 != read_common_params (cfg, "smc_config", &smc))
     return -1;

   if (0 != row_select_scan (time_beg, time_end,
                             pad_enable ? smc.num_pad : 0,
                             smc.file_glob_pattern, &smc.rst))
     goto return_status;

   if (smc.rst)
     {
        if (0 != radiance_copy_smc (r, smc.rst))
          goto return_status;
     }

   status = 0;
return_status:
   row_select_free (smc.rst);
   return status;
}

static int rename_radiance_file (Rename_Path_Type *rpt)
{
   char basename[BASENAME_SIZE];
   int ncid, status = -1;

   if (0 != TIO_open (rpt->tmp_path, NC_NOWRITE, &ncid))
     return -1;

   if (TIO_filename_from_granule (ncid, TEMPO_PROD_TYPE_RAD, 1, rpt->processing_version,
                                  basename, BASENAME_SIZE) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: generating filename", __func__);
        goto return_status;
     }

   if (NULL == (rpt->final_path = ioclib_pathconcat (rpt->target_dir, basename)))
     goto return_status;

   if (0 != rename (rpt->tmp_path, rpt->final_path))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: renaming %s to %s",
                     __func__, rpt->tmp_path, rpt->final_path);
        goto return_status;
     }

   status = 0;
return_status:
   (void) TIO_close (ncid);
   return status;
}

static char *temp_radiance_path (const char *target_dir)
{
   const char tmp_fmt[] = ".l1_inr_prep_%d_tmprad1.nc";
   char buf[BASENAME_SIZE];

   if (snprintf (buf, BASENAME_SIZE, tmp_fmt, getpid()) >= BASENAME_SIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: filename truncated (bufsize=%d is too small)",
                     __func__, BASENAME_SIZE);
        return NULL;
     }

   return ioclib_pathconcat (target_dir, buf);
}

static int read_rename_config (config_t *cfg, Rename_Path_Type *rpt)
{
   config_setting_t *s;
   const char *target_dir;

   if (NULL == (s = config_lookup (cfg, "telemetry_only_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'telemetry_only_config' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "target_dir", &target_dir))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'target_dir' from param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &rpt->processing_version))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'processing_version' from param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (rpt->target_dir = expand_string (target_dir)))
       || (NULL == (rpt->tmp_path = temp_radiance_path (rpt->target_dir))))
     return -1;

   return 0;
}

static void free_rename_path_type (Rename_Path_Type *rpt)
{
   if (rpt == NULL)
     return;
   ioclib_free (rpt->tmp_path);
   ioclib_free (rpt->final_path);
   FREE(rpt->target_dir);
}

static int copy_ephem (Radiance_Type *r, config_t *cfg,
                       double time_beg, double time_end, int pad_enable,
                       const char *ephemeris_file)
{
   Eph_Type eph = {0};
   config_setting_t *s;
   int num_pad, status = -1;

   if (NULL == (s = config_lookup (cfg, "ephemeris_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'ephemeris_config' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_pad", &num_pad))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'ephemeris_config' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 != eph_read_subset (&eph, ephemeris_file, time_beg, time_end,
                             pad_enable ? num_pad : 0))
     goto return_status;

   if (eph.n > 0)
     {
        if (0 != radiance_write_eph (r, &eph))
          goto return_status;
     }

   status = 0;
return_status:
   eph_free (&eph);

   return status;
}

static int process_inputs (config_t *cfg,
                           const char *radiance_file,
                           double time_beg, double time_end,
                           const char *ephemeris_file)
{
   Radiance_Type *r = NULL;
   Rename_Path_Type rpt = {0};
   const char *logmsg_filename = NULL;
   int radiance_is_telemetry_only = 0;
   int pad_enable = 1;
   int status = -1;

   if (radiance_file)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 0, "processing file: %s", radiance_file);
        logmsg_filename = radiance_file;
        if (NULL == (r = radiance_open (radiance_file)))
          return -1;
        if (0 != radiance_interval (r, &time_beg, &time_end))
          goto return_status;
     }
   else
     {
        /* If necessary, create a telemetry-only radiance file */
        tell_vlog (TELL_MSGTYPE_INFO, 0, "create telemetry-only radiance file: [%f, %f]", time_beg, time_end);
        radiance_is_telemetry_only = 1;
        if (0 != read_rename_config (cfg, &rpt))
          goto return_status;
        if (NULL == (r = radiance_create (rpt.tmp_path, rpt.processing_version)))
          goto return_status;
        logmsg_filename = rpt.tmp_path;
     }

   if (isnan(time_beg)
       || isnan(time_end)
       || (time_beg >= time_end))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: invalid time interval [%g, %g)",
                     __func__, time_beg, time_end);
        goto return_status;
     }

   /* Copy IRU, SMC time series spanning [time_beg,time_end),
    * handling the case where padding forces consideration
    * of additional files.
    */

   if (0 != copy_smc (r, cfg, time_beg, time_end, pad_enable))
     goto return_status;

   if (0 != copy_iru (r, cfg, time_beg, time_end, pad_enable))
     goto return_status;

   /* Copy subset of ephemeris */
   if (0 != copy_ephem (r, cfg, time_beg, time_end, pad_enable, ephemeris_file))
     goto return_status;

   if (radiance_is_telemetry_only)
     {
        /* Finalize the temporary radiance file */
        if (0 != radiance_update_coverage_times (r, time_beg, time_end))
          goto return_status;

        /* close before renaming */
        radiance_close (r);
        r = NULL;
        if (0 != rename_radiance_file (&rpt))
          goto return_status;
        logmsg_filename = rpt.final_path;
     }

   status = 0;
return_status:
   tell_vlog (TELL_MSGTYPE_INFO, 0, "exit status=%d, file=%s",
              status, logmsg_filename ? logmsg_filename : "(null)");
   radiance_close (r);
   free_rename_path_type (&rpt);

   return status;
}

static void delay_start (time_t n)
{
   struct timespec req;

   req.tv_sec = n;
   req.tv_nsec = 0;

   for (;;)
     {
        struct timespec rem;
        if (0 == nanosleep (&req, &rem))
          return;
        if (errno != EINTR) break;
        req = rem;
     }
}

int main (int argc, char **argv)
{
   const char appname[] = "L1_inr_prep";
   char *config_file = "l1_inr_prep.cfg";
   config_t cfg;
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"help",    no_argument,       0, 'h'},
        {"begin",   required_argument, 0, 'b'},
        {"end",     required_argument, 0, 'e'},
        {"epoch",   required_argument, 0, 'E'},
        {"delay",   required_argument, 0, 'd'},
        {"config",  required_argument, 0, 'c'},
        {"verbose", required_argument, 0, 'v'},
        {"ephemeris", required_argument, 0, 'p'},
        {0,0,0,0}
     };

   double nan_value = nan("");
   double time_beg = nan_value;
   double time_end = nan_value;
   time_t delay_sec = 0;
   char *radiance_file = NULL;
   char *ephemeris_file = NULL;
   int print_usage = 0;
   int have_epoch = 0;

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);

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
        int c = getopt_long (argc, argv, "hb:c:d:e:E:p:v:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "%s: getopt returned character %d??\n", __func__, c);
             goto return_status;
             break;
           case 'h':
             print_usage = 1;
             goto return_status;
             break;
           case 'b':
             if (1 != sscanf (optarg, "%le", &time_beg))
               goto return_status;
             break;
           case 'e':
             if (1 != sscanf (optarg, "%le", &time_end))
               goto return_status;
             break;
           case 'E':
               {
                  time_t epoch;
                  if (1 != sscanf (optarg, "%ld", &epoch))
                    goto return_status;
                  if (0 != tio_time_set_taix_epoch_timet (epoch))
                    goto return_status;
                  have_epoch++;
               }
             break;
           case 'd':
             if (1 != sscanf (optarg, "%ld", &delay_sec))
               goto return_status;
             break;
           case 'p':
             ephemeris_file = optarg;
             break;
           case 'c':
             config_file = optarg;
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             if (-1 == read_config_file (config_file, &cfg))
               goto return_status;
             break;
           case 'v':
               {
                  int log_level;
                  if (1 == sscanf (optarg, "%d", &log_level))
                    (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);
               }
             break;
          }
     }

   if (optind == 0)
     {
        print_usage = 1;
        goto return_status;
     }

   if (optind < argc)
     {
        radiance_file = argv[optind++];
     }

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored: \n");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   if (ephemeris_file == NULL)
     {
        fprintf (stderr, "*** Ephemeris file not specified\n");
        print_usage = 1;
        goto return_status;
     }

   if (delay_sec > 0)
     {
        delay_start (delay_sec);
     }

   /* When either time boundary is set, the epoch must also be set */
   if ((have_epoch == 0)
       && ((0 == isnan(time_beg)) || (0 == isnan(time_end))))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: epoch is not set", __func__);
        goto return_status;
     }

   if (0 != process_inputs (&cfg, radiance_file, time_beg, time_end, ephemeris_file))
     goto return_status;

   status = 0;
return_status:
   config_destroy (&cfg);
   tell_close ();

   if (print_usage) usage();

   return status;
}
