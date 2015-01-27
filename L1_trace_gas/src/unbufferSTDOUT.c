#include <stdio.h>
#include <stdlib.h>
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

