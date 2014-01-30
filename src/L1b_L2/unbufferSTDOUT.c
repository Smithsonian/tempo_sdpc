#include <stdio.h>
#include <stdlib.h>
void unbufferstdout_(void)
{
   if (getenv ("OMIABORT") != NULL)
     atexit (abort);
   setbuf(stdout,NULL);
}
