#include "config.h"

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <slang.h>
#include <tell.h>
#include "version.h"

SLANG_MODULE(daemon);

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
        SLang_set_error (SL_INTRINSIC_ERROR);
        return -1;
     }

   Tell_open (appname, fd, 0);
   return 0;
}

/* Initial version from Stevens APUE */
static void daemonize_intrin (const char *appname, const char *logfile_path)
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
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   /* Become a session leader to lose controlling TTY. */
   if ((pid = fork()) < 0)
     {
        print_errmsg ("%s: can't fork", appname);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
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
   if (sigaction(SIGHUP, &sa, NULL) < 0)
     {
        print_errmsg ("%s: can't ignore SIGHUP", appname);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }
   if ((pid = fork()) < 0)
     {
        print_errmsg ("%s: can't fork", appname);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }
   else if (pid != 0) /* parent */
     {
        _exit(0);
     }

   /* Change the current working directory to the root so
    * we won't prevent file systems from being unmounted. */
   if (chdir("/") < 0)
     {
        print_errmsg ("%s: can't change directory to /", appname);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
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
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   if (fd0 != 0 || fd1 != 1 || fd2 != 2)
     {
        tell_vlog (TELL_MSGTYPE_ERROR, -1,
                   "unexpected file descriptors %d %d %d", fd0, fd1, fd2);
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }
}

static void daemon_log_intrin (int *msg_type, const char *msg)
{
   tell_vlog (*msg_type, -1, "%s", msg);
}

static SLang_Intrin_Fun_Type Module_Intrinsics [] =
{
   MAKE_INTRINSIC_2 ("daemonize", daemonize_intrin, SLANG_VOID_TYPE, SLANG_STRING_TYPE, SLANG_STRING_TYPE),
   MAKE_INTRINSIC_2 ("daemon_log", daemon_log_intrin, SLANG_VOID_TYPE, SLANG_INT_TYPE, SLANG_STRING_TYPE),
   SLANG_END_INTRIN_FUN_TABLE
};

static SLang_IConstant_Type Module_IConstants [] =
{
   MAKE_ICONSTANT("_daemon_module_version", MODULE_VERSION_NUMBER),
   MAKE_ICONSTANT("LOG_ERR", TELL_MSGTYPE_ERROR),
   MAKE_ICONSTANT("LOG_WARN", TELL_MSGTYPE_WARN),
   MAKE_ICONSTANT("LOG_INFO", TELL_MSGTYPE_INFO),
   SLANG_END_ICONST_TABLE
};

static SLang_Intrin_Var_Type Module_Variables [] =
{
   MAKE_VARIABLE("_daemon_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

int init_daemon_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns = SLns_create_namespace (ns_name);
   if (ns == NULL)
     return -1;

   if (-1 == SLns_add_intrin_fun_table (ns, Module_Intrinsics, NULL)
       || -1 == SLns_add_iconstant_table (ns, Module_IConstants, NULL)
       || -1 == SLns_add_intrin_var_table (ns, Module_Variables, NULL))
     {
        return -1;
     }

   return 0;
}

/* This function is optional */
void deinit_daemon_module (void)
{
   Tell_close();
}

