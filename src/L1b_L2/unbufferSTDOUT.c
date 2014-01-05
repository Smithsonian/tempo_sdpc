#include <stdio.h>
#include <stdlib.h>
void unbufferstdout_()
{
   /* atexit (abort); */
   setbuf(stdout,NULL);
}
