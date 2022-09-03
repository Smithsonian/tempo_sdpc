/*------------------------------------------------------------------------------*\
** ymlNrmKey.c
**
** - removes white space from input string "a"
** - changes all characters to upper-case
**
**   LDx = level of detail
**   a[] = input  string, null-terminated
**   b[] = output string, null-terminated
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlNrmKey(int LDx, char a[], char b[])
{
    int len = strlen(a);			/* length of input string	*/
    int max = YML_KEY_LEN - 1;			/* one less due to terminator	*/
    int rtc = 0;				/* return code			*/

   (void) LDx;

    if ( len > max )
    {
       Error("key length %d exceeds %d chars; truncating to %d.", len, max, max);
       Info(0, "<%s>", a);
       rtc = 4;
       /* not fatal; string truncated, user informed				*/
    }

    int j = 0; int k = 0; for( k = 0; k < YML_KEY_LEN; k++ )
    {
        if ( a[k] ==     ' ' ) continue;	/* remove white space		*/
        if ( a[k] ==    '\0' ) break;		/* end of source string		*/

        b[j++] = toupper(a[k]);			/* convert to upper case	*/
    }

    b[j] = '\0';				/* terminate destination	*/

    return(rtc);
}

/*------------------------------------------------------------------------------*/
