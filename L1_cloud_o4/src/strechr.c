/*------------------------------------------------------------------------------*\
** strechr.c
**
** Locates the first occurrence from the end that is not c in string pointed
** to by s.  Returns pointer to the last c.  For example,
** strechr("abcdefg    \0", ' ') would find "g" and point to the space just
** following the "g".  Returns NULL if string has no c characters.
** This is useful in removing trailing blanks in a string.
\*------------------------------------------------------------------------------*/

#include <string.h>
#include "O3_PEATE_Common.h"

char *strechr(char *s, char c)
{
     char *x = strchr(s, '\0') - 1;
     while ( x >= s )
     {
        if ( *x != c ) return x+1;
        x--;
     }
     return NULL;
}

/*------------------------------------------------------------------------------*/
