#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>

#include <ioclib.h>
#include <skalibs/tai.h>

/* Somewhat arbitrary upper limit on the maximum valid time_t value */
#define MAX_VALID_TIMET (100*365L*86400L)

typedef struct
{
   regex_t reg;
   tain tbeg;
   tain tend;
   int have_regex;
   int have_time_window;
   int case_insensitive;
   int exclude_match;
}
Filter_Type;

static void filter_free (Filter_Type *filter)
{
   if (filter == NULL)
     return;
   if (filter->have_regex)
     {
        regfree (&filter->reg);
     }
}

static int filter_init_regex (Filter_Type *filter, const char *regex)
{
   int cflags = REG_EXTENDED | REG_NOSUB | REG_NEWLINE;
   int errcode;

   if (regex == NULL)
     {
        filter->have_regex = 0;
        return 0;
     }

   if (filter->case_insensitive)
     {
        cflags |= REG_ICASE;
     }
   if (0 != (errcode = regcomp (&filter->reg, regex, cflags)))
     {
        fprintf (stderr, "*** ERROR: compiling regular expression (err=%d): %s\n", errcode, regex);
        return -1;
     }
   filter->have_regex = 1;

   return 0;
}

static int filter_init_time_window (Filter_Type *filter,
                                    const struct timeval *tbeg, const struct timeval *tend)
{
   if (0 == tain_from_timeval_sysclock (&filter->tbeg, tbeg))
     {
        fprintf (stderr, "*** Error: tain_from_timeval_sysclock conversion failed: tbeg=%ld\n", tbeg->tv_sec);
        return -1;
     }

   if (0 == tain_from_timeval_sysclock (&filter->tend, tend))
     {
        fprintf (stderr, "*** Error: tain_from_timeval_sysclock conversion failed tend=%ld\n", tend->tv_sec);
        return -1;
     }
   filter->have_time_window = (0 != tain_less (&filter->tbeg, &filter->tend));

   return 0;
}

static int filter_match_line (const Filter_Type *filter, const char *linep)
{
   /* time filter */
   if (filter->have_time_window)
     {
        tain tn;
        /* When linep starts with a valid TAI64N timestamp,
         * timestamp_scan parses that value, and returns 25.
         * Otherwise, it returns 0.
         */
        if (0 == timestamp_scan (linep, &tn))
          return 0;
        if ((0 != tain_less (&tn, &filter->tbeg))
            || (0 != tain_less (&filter->tend, &tn)))
          {
             return 0;
          }
     }

   /* regular expression filter */
   if (filter->have_regex)
     {
        if (0 != regexec (&filter->reg, linep, 0, NULL, 0))
          return 0;
     }

   /* Get to here only when all filters match */
   return 1;
}

static int process_logfile (const Filter_Type *filter, const char *dirpath, const char *file)
{
   FILE *fp = NULL;
   char *path = NULL;
   int status = -1;

   if (*file == '-')
     {
        fp = stdin;
     }
   else
     {
        /* Log files with names that correspond to a valid TAIN64 timestamp
         * are named using the last timestamp entry in the file.
         * When that last entry is earlier than our filter's time window,
         * there's no point in opening the file, so we return immediately.
         * When the log file name begins with '@' but is not a valid
         * TAIN64 timestamp, we open the file and process it normally.
         */
        if ((*file == '@') && (filter->have_time_window != 0))
          {
             tain tn;
             if ((25 == timestamp_scan (file, &tn))
                 && (0 != tain_less (&tn, &filter->tbeg)))
               {
                  return 0;
               }
          }

        if (NULL == (path = ioclib_pathconcat (dirpath, file)))
          {
             fprintf (stderr, "*** Error: %s: ioclib_pathconcat failed\n", __func__);
             return -1;
          }

        if (NULL == (fp = fopen (path, "r")))
          {
             fprintf (stderr, "*** Error: failed opening file for reading: %s\n", path);
             return -1;
          }
     }

   for (;;)
     {
        char *linep;
        size_t lenp;
        int ret;
        if ((ret = ioclib_fgets (&linep, &lenp, fp)) < 0)
          {
             fprintf (stderr, "*** Error: %s: read failed\n", __func__);
             goto return_status;
          }
        else if (ret == 0)
          break; /* EOF */

        if (filter_match_line (filter, linep) != filter->exclude_match)
          {
             ret = fputs (linep, stdout);
          }
        else ret = 0;

        ioclib_free (linep);

        if (ret == EOF) /* EOF on stdout */
          break;
     }

   status = 0;
return_status:
   ioclib_free (path);
   if (fp != stdin)
     {
        fclose (fp);
     }
   return status;
}

/* file_cbfun callback returns:
 *   1 to continue processing;
 *   0 to stop processing so that ioclib_process_dir returns 0;
 *  -1 to stop processing so that ioclib_process_dir returns -1.
 */
static int file_cbfun (int dirfd, const char *dirpath, const char *file, unsigned int idx, void *cd)
{
   int status;
   (void) dirfd; (void) idx;
   status = process_logfile ((const Filter_Type *)cd, dirpath, file);
   return (status == 0) ? 1 : status;
}

static int process_logdir (Filter_Type *filter, const char *dirpath)
{
   const char *globs[] = {"@*.s", "current"};
   unsigned int nglobs = sizeof(globs)/sizeof(*globs);
   if (ioclib_process_dir2 (dirpath, globs, nglobs, IOCLIB_LISTDIR_SORT, file_cbfun, filter) < 0)
     return -1;
   return 0;
}

static int classify_path (const char *path, char **dirpath, char **file)
{
   if (0 != ioclib_isdir (path, "rx"))
     {
        if (NULL == (*dirpath = strdup (path)))
          {
             fprintf (stderr, "*** Error: strup failed\n");
             return -1;
          }
        *file = NULL;
        return 0;
     }
   else if (0 != ioclib_isfile (path, "r"))
     {
        if (NULL == (*dirpath = ioclib_dirname (path)))
          {
             fprintf (stderr, "*** Error: ioclib_dirname failed\n");
             return -1;
          }
        *file = ioclib_basename (path);
        return 0;
     }

   fprintf (stderr, "*** Error: cannot read from %s\n", path);
   return -1;
}

static int usage (void)
{
   (void) fprintf (stderr, "Usage: filter_s6_log [options] PATH\n");
   (void) fprintf (stderr, "Options:\n");
   (void) fprintf (stderr, " -h|--help          Print this usage message\n");
   (void) fprintf (stderr, " -b|--begin TIME    Begin time (time_t)\n");
   (void) fprintf (stderr, " -e|--end TIME      End time (time_t)\n");
   (void) fprintf (stderr, " -r|--regex FILTER  Regular expression filter\n");
   (void) fprintf (stderr, " -i|--insensitive   Case insensitive matching\n");
   (void) fprintf (stderr, " -x|--exclude       Exclude matching lines\n");
   (void) fprintf (stderr, "\n");
   (void) fprintf (stderr, "PATH may be either a log file, or a log directory.\n");
   (void) fprintf (stderr, "To filter a log file piped from stdin, use PATH='-'\n");
   (void) fprintf (stderr, "\n");
   exit(0);
}

int main (int argc, char **argv)
{
   int exit_status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"begin",       required_argument, 0, 'b'},
        {"end",         required_argument, 0, 'e'},
        {"help",        no_argument,       0, 'h'},
        {"insensitive", no_argument,       0, 'i'},
        {"regex",       required_argument, 0, 'r'},
        {"exclude",     no_argument,       0, 'x'},
        {0,0,0,0}
     };
   char *regex = NULL;
   char *path = NULL;
   char *dirpath = NULL;
   char *file = "-";
   struct timeval tbeg = {0};
   struct timeval tend = {MAX_VALID_TIMET, 0};
   Filter_Type filter = {0};

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "b:e:hir:x", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "*** unrecognized option: getopt_long returned character %d??\n", c);
             goto return_status;
             break;
           case 'h':
             usage ();
             break;

           case 'b':  /* filter begin time */
             if (1 != sscanf (optarg, "%ld", &tbeg.tv_sec))
               {
                  fprintf (stderr, "*** ERROR: could not parse begin time = %s\n", optarg ? optarg : "<null>");
                  goto return_status;
               }
             break;
           case 'e':  /* filter end time */
             if (1 != sscanf (optarg, "%ld", &tend.tv_sec))
               {
                  fprintf (stderr, "*** ERROR: could not parse end time = %s\n", optarg ? optarg : "<null>");
                  goto return_status;
               }
             break;

           case 'i':
             filter.case_insensitive = 1;
             break;
           case 'r':  /* regular expression */
             regex = optarg;
             break;
           case 'x':  /* exclude matching lines */
             filter.exclude_match = 1;
             break;
          }
     }

   if (optind < argc)
     {
        path = argv[optind++];
     }

   if (optind < argc)
     {
        fprintf (stderr, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stderr, "%s ", argv[optind++]);
          }
        fprintf (stderr, "\n");
     }

   if (0 != filter_init_regex (&filter, regex))
     goto return_status;

   if (0 != filter_init_time_window (&filter, &tbeg, &tend))
     goto return_status;

   if (path)
     {
        if (0 != classify_path (path, &dirpath, &file))
          goto return_status;
     }

   if (file == NULL)
     {
        if (0 != process_logdir (&filter, dirpath))
          goto return_status;
     }
   else
     {
        if (0 != process_logfile (&filter, dirpath, file))
          goto return_status;
     }

   exit_status = EXIT_SUCCESS;
return_status:
   filter_free (&filter);
   ioclib_free (dirpath);
   return exit_status;
}
