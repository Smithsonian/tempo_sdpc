/*------------------------------------------------------------------------------*\
** unbufferSTDOUT.c
**
** Causes STDOUT to behave similarly to STDERR in that it will not be buffered.
\*------------------------------------------------------------------------------*/

#include <stdio.h>

void unbufferstdout_() { setbuf( stdout, NULL ); }

/*------------------------------------------------------------------------------*/
