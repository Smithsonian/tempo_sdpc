/** @file daemon.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Daemon process initialization
 */

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <ioclib.h>
#include <tell.h>

/* stdio isn't optimal for signal handling, but this
 * function is only called if something goes wrong
 * during initialization.
 */
static void print_errmsg (const char *fmt, ...)
{
   va_list ap;
   va_start (ap, fmt);
   (void) vfprintf (stderr, fmt, ap);
   va_end (ap);
   (void) fputc ('\n', stderr);
}

static int open_logfile (const char *appname, const char *path)
{
   int fd, flags = O_WRONLY | O_CREAT;
   mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

   if ((fd = open (path, flags, mode)) < 0)
     {
        print_errmsg ("%s: open failed: %s", appname, path);
        return -1;
     }

   tell_open (appname, fd, 0);
   return 0;
}

static const char *Temp_Pid_File = NULL;
static void delete_pidfile (void)
{
   if (Temp_Pid_File)
     {
        (void) ioclib_unlink (Temp_Pid_File);
     }
}

static int make_pidfile (const char *appname)
{
#define BUFSIZE 256
   char dirname[BUFSIZE];
   char basename[BUFSIZE];
   const char *user = NULL;
   char *path = NULL;
   FILE *fp = NULL;

   if (appname == NULL)
     return -1;

   if (NULL == (user = getenv ("USER")))
     user = "";

   snprintf (dirname, sizeof(dirname), "/var/tmp/%s/%s", user, appname);
   if (0 != ioclib_mkdir (dirname, 0))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: creating %s",
                     __func__, dirname);
        return -1;
     }

   snprintf (basename, sizeof(basename), "%d", getpid());
   if (NULL == (path = ioclib_pathconcat (dirname, basename)))
     return -1;

   if ((NULL == (fp = fopen (path, "w")))
       || (0 != fclose (fp)))
     {
        ioclib_free(path);
        tell_verror (TELL_IO_WRITE_ERROR, "%s: creating %s (%s)",
                     __func__, path, strerror(errno));
        return -1;
     }

   Temp_Pid_File = path;
   atexit (&delete_pidfile);

   return 0;
}

/* Initial version from Stevens APUE */
int daemonize (const char *appname, const char *logfile_path)
{
   unsigned int i;
   int fd0, fd1, fd2;
   pid_t pid;
   struct rlimit rl;
   struct sigaction sa;

   /* Clear file creation mask. */
   umask(0);

   /* Get maximum number of file descriptors. */
   if (getrlimit(RLIMIT_NOFILE, &rl) < 0)
     {
        print_errmsg ("%s: can't get file limit", appname);
        return -1;
     }

   /* Become a session leader to lose controlling TTY. */
   if ((pid = fork()) < 0)
     {
        print_errmsg ("%s: can't fork", appname);
        return -1;
     }
   else if (pid != 0) /* parent */
     {
        _exit(0);
     }
   setsid();

   /* Ensure future opens won't allocate controlling TTYs. */
   sa.sa_handler = SIG_IGN;
   sigemptyset(&sa.sa_mask);
   sa.sa_flags = 0;
   if (sigaction (SIGHUP, &sa, NULL) < 0)
     {
        print_errmsg ("%s: can't ignore SIGHUP", appname);
        return -1;
     }
   if ((pid = fork()) < 0)
     {
        print_errmsg ("%s: can't fork", appname);
        return -1;
     }
   else if (pid != 0) /* parent */
     {
        fprintf (stdout, "%s started pid= %d log= %s\n",
                 appname, pid, logfile_path);
        _exit(0);
     }

   /* Change the current working directory to the root so
    * we won't prevent file systems from being unmounted. */
   if (chdir("/") < 0)
     {
        print_errmsg ("%s: can't change directory to /", appname);
        return -1;
     }

   /* Close all open file descriptors. */
   if (rl.rlim_max == RLIM_INFINITY)
     rl.rlim_max = 1024;
   for (i = 0; i < rl.rlim_max; i++)
     {
        close(i);
     }

   /* Attach file descriptors 0, 1, and 2 to /dev/null. */
   fd0 = open("/dev/null", O_RDWR);
   fd1 = dup(0);
   fd2 = dup(0);

   if (-1 == open_logfile (appname, logfile_path))
     return -1;

   if (fd0 != 0 || fd1 != 1 || fd2 != 2)
     {
        tell_vlog (TELL_MSGTYPE_ERROR, -1,
                   "unexpected file descriptors %d %d %d", fd0, fd1, fd2);
        return -1;
     }

   return make_pidfile (appname);
}
