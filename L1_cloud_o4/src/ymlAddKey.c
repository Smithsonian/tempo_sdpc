/*------------------------------------------------------------------------------*\
** ymlAddKey.c
**
** Concatenates parts of key into one normalized string.
**
**   LDx = level of detail
**   nry = number of elements in tbl[]
**   tbl = table for YAML data
**   lvl = number of elements in stk[]
**   stk = stack of key parts
**   max = size of stk[]
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlAddKey(int LDx, int nry, yTB tbl[], int lvl, char* stk[], int max)
{
    if ( lvl >= max ) { return(0); }		/* error reported previously	*/
    if ( cfx >= nry ) { return(0); }		/* error reported previously	*/

    if ( max > YML_MAX_STK )
    {
       Error("max > YML_MAX_STK and that should never happen");
       return(-14);
    }

    int rtc = 0;				/* return code			*/

    int cnt; for( cnt = 0; cnt <= lvl; cnt++ )	/* copy stk[] to tbl[].stk	*/
    {
        if ( cnt < max )
        {
           tbl[cfx].stk[cnt] = stk[cnt];	/* pointers to key parts	*/
        }
    }

    ymlDmpStk(LDx, tbl);			/* show stack (all to add up)	*/

    char buf[YML_KEY_LEN];			/* temp buffer related to stk[]	*/

    size_t ptr = 0;				/* index in buf[]		*/

    int idx; for (idx = 0; idx < max; idx++)	/* walk items on stack		*/
    {
        char* itm = tbl[cfx].stk[idx];		/* convenience variable		*/

        Info(LDx, "OK considering <%s>", itm);

        if ( itm == '\0' ) { break; }		/* end of parts list		*/

        size_t len = strlen(itm);		/* size of item on stk		*/
        size_t end = ptr + len;			/* projected end after copy	*/

        size_t stp = YML_KEY_LEN - 1;		/* one less due to terminator	*/

        if ( end > stp )			/* projected overflow?		*/
        {
           int dif = end - stp;			/* amount of overflow		*/

           Error("buf[YML_KEY_LEN] overflow: %d > %d by %d", end, stp, dif);
           Info(L0, "<%s>", itm);

           rtc = 2;

           /* not fatal since buf[] will be truncated if it is too long		*/
        }

        strcpy(&buf[ptr], itm);			/* add part to accumulation	*/
        ptr = end;				/* move pointer to new end	*/
    }

    /*--------------------------------------------------------------------------*/

    if ( tbl[cfx].stk[0] == '\0' )		/* no key; value part of array	*/
    {
        strcpy(&buf[ptr], tbl[cfx-1].key);	/* use previous value's key	*/
    }

    /*--------------------------------------------------------------------------*/

    ymlNrmKey(LDx, buf, tbl[cfx].key);		/* normalize; cp to destination	*/

    Info(LDx, "OK added <%s>", tbl[cfx].key);

    return(rtc);
}

/*------------------------------------------------------------------------------*/
