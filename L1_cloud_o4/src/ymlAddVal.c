/*------------------------------------------------------------------------------*\
** ymlAddVal.c
**
** Add one entry into the internal table of parameter-value pairs.
**
**   LDx = level of detail
**   nry = number of elements in tbl[]
**   tbl = table for YAML data
**   lvl = YAML level number
**   blk = YAML block number
**   val = value (string) to add
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlAddVal(int LDx, int nry, yTB tbl[], int lvl, int blk, char* val)
{
    Info(LDx, "OK #################### ymlAddVal() lvl=%d blk=%d cfx=%d val=<%s>",
                                                   lvl,   blk,   cfx,   val);

    int rtc = 0;				/* return code			*/

    /*--------------------------------------------------------------------------*/
    /* add value								*/

    size_t len = strlen(val);			/* length of val (no end null)	*/

    size_t stp = YML_VAL_LEN - 1;		/* one less due to terminator	*/

    if ( len > stp )
    {
       Error("yTB[].val[YML_VAL_LEN] overflow: %ld > %ld", len, stp);
       rtc = 8;
       /* not fatal; value will be truncated if too long			*/
    }

    if ( cfx < nry )
    {
       strncpy(tbl[cfx].val, val, stp);		/* cp val to final destination	*/

       tbl[cfx].val[stp] = '\0';		/* terminate if too long	*/

       tbl[cfx].lvl = lvl;			/* YAML level number		*/
       tbl[cfx].blk = blk;			/* YAML block number		*/
       tbl[cfx].usd = 0;			/* item use count (init)	*/

       Info(LDx, "OK ymlAddVal() yTB[%d] lvl=%d blk=%d key=<%s> val=<%s>",
                     cfx, tbl[cfx].lvl, tbl[cfx].blk, tbl[cfx].key, tbl[cfx].val);

       /*-----------------------------------------------------------------------*/
       /* if item is "Maximum Level Of Detail", set the global value		*/

       char* key = tbl[cfx].key;			/* convenience variable	*/

       rtc = ymlMaxDtl(LDx, key, val);	if ( rtc < 0 ) { return(rtc); }
    }

    /*--------------------------------------------------------------------------*/
    /* test if there is space to add the next entry				*/

    cfx++;				/* add one for entry to be added next	*/
    stx = 0;				/* re-init index for writing stk[]	*/

    if ( cfx >= nry )
    {
       Error("yTB[%d] overflow; increase dimension to %d", nry, cfx+1);
       /* not fatal; want to find needed max dimension				*/
       rtc = 3;
    }

    return(rtc);
}

/*------------------------------------------------------------------------------*/
