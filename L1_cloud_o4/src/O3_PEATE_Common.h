/*------------------------------------------------------------------------------*\
** O3_PEATE_Common.h
**
** Includes definitions for a variety of useful things that can be shared.
\*------------------------------------------------------------------------------*/

#include <stddef.h>		/* defines NULL, a generic "no value" pointer	*/
#include <ctype.h>		/* defines isalpha(), etc.			*/

#define bool  int		/* boolean is same as int			*/
#define TRUE  (1==1)
#define FALSE (1==0)
#define UNDEF ( -2 )

#define ReturnIfFatal(rco) if ( rco < 0 ) { return rco; }

char *strechr(char *s, char  c);
char *strxchr(char *s, char  c);
int   str2upr(char *s, char *r);
int   str2lwr(char *s, char *r);

bool  SetTrueIfLittleEndian();

/*------------------------------------------------------------------------------*/
