/*------------------------------------------------------------------------------*\
** str2lwr.c
**
** Converts a string to all lower case.  For example, str2lwr("A1b=C2d\0", out)
** would result in "a1b=c2d\0".  Strings s and r may be the same string or
** different strings.  Returns the number of characters in the string (including
** the terminating null).
\*------------------------------------------------------------------------------*/

#include "O3_PEATE_Common.h"

int str2lwr(char *s, char *r)
{
    int i = 1; while( *s )
    {
       if ( isalpha( (int)(*s) ) ) { *r = (char) tolower( (int)(*s) ); }
       else                        { *r =                       *s;    }
       s++; r++; i++;
    }
    *r = '\0';
    return i;
}

/*------------------------------------------------------------------------------*/
