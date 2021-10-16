/*------------------------------------------------------------------------------*\
** ymlIniTbl.c
**
** Initializes the main memory structure tbl[].
**
**   LDx = level of detail
**   nrp = pointer to number of elements in tbl[]
**   ptr = pointer to tbl[]
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int cfx;					/* index for writing tbl[]	*/
int stx;					/* index for writing stk[]	*/

int ymlIniTbl(int LDx, int* nrp, yTB** ptr)
{
    int   msz = EWI_LEN;			/* error msg[] buffer size	*/
    char  msg  [EWI_LEN];			/* error or progress message	*/
    int   ret = 0;				/* return code			*/

    int   nry = *nrp;				/* number of elements in tbl[]  */

    /*--------------------------------------------------------------------------*/

    size_t esz = (size_t)( sizeof(yTB) );	/* number of bytes in tbl[0]	*/
    size_t fsz = (size_t)nry * esz;		/* number of bytes in tbl[nry]	*/

    Info(LDx, "OK each element in tbl[] is %d bytes", esz);

    snprintf(msg, msz, "allocating %d bytes of memory for tbl[%d]", fsz, nry);

    *ptr = (yTB*)malloc(fsz);

    yTB* tbl = *ptr;				/* tbl[] is what ptr points to	*/

    if ( tbl == NULL ) { Error(msg); ret = -1; return(ret); }

    Info(LDx, "OK %s", msg);

    /*--------------------------------------------------------------------------*/

    for( cfx = 0; cfx < nry; cfx++ )
    {
        tbl[cfx].key[0] = '\0';			/* parameter name		*/
        tbl[cfx].val[0] = '\0';			/* parameter value		*/
        tbl[cfx].lvl    =   0 ;			/* YAML level number		*/
        tbl[cfx].blk    =   0 ;			/* YAML block number		*/
        tbl[cfx].usd    =   0 ;			/* item use count		*/

        for( stx = 0; stx < YML_MAX_STK; stx++ )
        {
           tbl[cfx].stk[stx] = NULL;		/* stack of key part pointers	*/
        }
    }

    cfx = 0;					/* index for writing tbl[]	*/
    stx = 0;					/* index for writing stk[]	*/

    /*--------------------------------------------------------------------------*/

    return(ret);
}

/*------------------------------------------------------------------------------*/
