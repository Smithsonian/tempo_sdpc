/*------------------------------------------------------------------------------*\
** strxchr.c
**
** Locates the first occurrence that is not c in string pointed to by s.
** For example, strxchr("   abcd   \0", ' ') would return a pointer to the "a".
** Returns NULL if string is all c characters terminated by a null.
** This is useful in removing leading blanks or zeroes in a string or finding
** the beginning of useful information in a character stream.
\*------------------------------------------------------------------------------*/

#include "O3_PEATE_Common.h"

char *strxchr(char *s, char c)
{
     char *x = s; while ( *x ) { if ( *x != c ) return x; x++; } return NULL;
}

/*------------------------------------------------------------------------------*/
