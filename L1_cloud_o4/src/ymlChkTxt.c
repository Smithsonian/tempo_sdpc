/*------------------------------------------------------------------------------*\
** ymlChkTxt.c
**
** Points out everything in the text that breaks the YAML parser.
**
**   LDx = level of detail
**   str = string of ASCII characters, null-terminated
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlChkTxt(int LDx, char str[])
{
    int rtc = 0;				/* return code			*/
    int ix1 = 1;				/* counter of lines in file	*/
    int ix2 = 0;				/* counter of columns per line	*/

    size_t cnt = strlen(str);			/* count of chars in string	*/

    int ix3; for (ix3 = 0; ix3 < cnt; ix3++)	/* examine every char in string	*/
    {
        if ( str[ix3] == '\n' ) { ix1++; ix2 = 0; }	/* new line begins	*/

        if ( str[ix3] == '\t' )				/* YAML outlawed TABs	*/
        {
           Error("YAML file has TAB character: line %d column %d", ix1, ix2);
           rtc = -7;
        }
        ix2++;					/* counter of columns per line	*/
    }

    return(rtc);
}

/*------------------------------------------------------------------------------*/
