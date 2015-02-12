#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "tell.h"

#define MSG_ENV "PGSMSG"
#ifndef MSG_DIR
# define MSG_DIR ""
#endif

/* this stuff really ought to go into an include file */
extern void unbufferstdout_(void);
extern void c_exit_ (int *);
extern void maybe_setenv_msgenv_ (int *);

void unbufferstdout_(void)
{
   if (getenv ("OMIABORT") != NULL)
     atexit (abort);
   setbuf(stdout,NULL);
}

void c_exit_ (int *status)
{
   if (status)
     exit (*status);
   else
     exit (0);
}

static int __maybe_setenv_msgenv (void)
{
   struct stat st;
   const char msg_env[] = MSG_ENV;
   const char msg_dir[] = MSG_DIR;

   /* if the environment variable is already set, don't override */
   if (NULL != getenv (msg_env))
     return 0;

   if (-1 == stat (msg_dir, &st))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "cannot stat message directory: %s",
                     msg_dir);
        return -1;
     }

   if (0 == (S_ISDIR(st.st_mode))
       || (0 == (st.st_mode & S_IRUSR)))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "path is not a readable directory: %s",
                     msg_dir);
        return -1;
     }

   if (-1 == putenv (MSG_ENV "=" MSG_DIR))
     {
        Tell_verror (TELL_IO_ERROR, "failed setting environment variable %s='%s'",
                     msg_env, msg_dir);
        return -1;
     }

   return 0;
}

void maybe_setenv_msgenv_ (int *errstat)
{
   if (errstat == NULL) return;
   *errstat = __maybe_setenv_msgenv ();
}

