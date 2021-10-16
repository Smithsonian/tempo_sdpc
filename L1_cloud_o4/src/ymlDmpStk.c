/*------------------------------------------------------------------------------*\
** ymlDmpStk.c
**
** Dump internal stack (stack of pointers to parts of key).
**
**   LDx = level of detail
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlDmpStk(int LDx, yTB tbl[])
{
    for( stx = 0; stx < YML_MAX_STK; stx++ )			  /* loop stk[]	*/
    {
        Info(LDx, "OK stk[%02d] = %s", stx, tbl[cfx].stk[stx]);   /* one line	*/

        if ( tbl[cfx].stk[stx] == 0 ) { break; }		  /* end	*/
    }

    return(0);
}

/*------------------------------------------------------------------------------*/
