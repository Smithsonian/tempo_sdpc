/*------------------------------------------------------------------------------*\
** ymlClrUsdFlg.c
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlClrUsdFlg(int nry, yTB tbl[], char* key)
{
    int LDa = L7;				/* level of detail		*/

    if ( tbl[0].key[0] == '\0' )		/* tbl[] not initialized	*/
    {
        Error("tbl[] has not been initialized; call ymlOpen() first");
        return(-10);
    }

    char nky[YML_KEY_LEN]; size_t len = YML_KEY_LEN; ymlNrmKey(LDa, key, nky);

    for( cfx = 0; cfx < nry; cfx++ )		/* loop all entries in tbl[]	*/
    {
        if ( 0 == strncmp(nky, tbl[cfx].key, len) ) { tbl[cfx].usd = 0; }
    }
    return(0);
}

/*------------------------------------------------------------------------------*/
