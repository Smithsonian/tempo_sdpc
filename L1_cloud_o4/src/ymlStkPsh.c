/*------------------------------------------------------------------------------*\
** ymlStkPsh.c
**
** Push an entry onto an internal stack.
**
**   LDx = level of detail
**   str = string to add
**   lvl = index in stk[]
**   stk = stack to add onto
**   max = size of stk[]
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlStkPsh(int LDx, char* str, int lvl, char* stk[], int max)
{
    Info(LDx, "OK --------------------- ymlStkPsh()  lvl=%d  str=<%s>", lvl, str);

    int rtc = 0;				/* return code			*/

    if ( lvl < max ) { stk[lvl] = str; }
    else
    {
       Error("stk[] overflow; %d > %d", lvl+1, max);
       rtc = 5;
       /* not fatal; want to find the max stack size needed			*/
    }

    return(rtc);
}

/*------------------------------------------------------------------------------*/
