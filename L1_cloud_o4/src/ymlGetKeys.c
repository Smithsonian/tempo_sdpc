/*------------------------------------------------------------------------------*\
** ymlGetKeys.c
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlGetKeys(int nry, yTB tbl[], int lvl, int max, char arr[max][YML_VAL_LEN])
{
    int LDa = L7;					/* level of detail	*/

    if ( tbl[0].key[0] == '\0' )
    {
        Error("tbl[] has not been initialized; call ymlOpen() first");
        return(-10);
    }

    int rtc = 0;					/* return code		*/

    int num = 0;			/* counter of matching keys found	*/

    int ix1; for( ix1 = 0; ix1 < max; ix1++ ) { arr[ix1][0] = '\0'; }	/* init	*/

    int ixo = 0;					/* nr of items found	*/

    /* empty entry is an early end of list, but arr[] may be filled completely	*/

    int ix2; for( ix2 = 0; ix2 < nry; ix2++ )		/* searching tbl[]	*/
    {
        char* key = tbl[ix2].key;			/* convenience variable	*/

        if ( tbl[ix2].key[0] == '\0' ) {    break; }	/* done; end of tbl[]	*/
        if ( tbl[ix2].lvl    !=  lvl ) { continue; }	/* skip: do not want it	*/

        if ( strcmp(key, tbl[ix2-1].key) != 0 )		/* not equal previous	*/
        {
           if ( ixo < max )				/* prevents overflow	*/
           {
              strncpy(arr[ixo], key, YML_VAL_LEN);	/* grab it, cp to arr[]	*/

              arr[ixo][YML_VAL_LEN-1] = '\0';		/* terminate to be safe	*/
           }

           ixo++;					/* next slot in arr[]	*/
        }
    }

    if ( ixo > max )
    {
       Error("arr[] overflow; ixo=%d > max=%d; too many keys", ixo, max);
       rtc = 7;
       /* overflow is not fatal since arr[] was truncated; no damage was done	*/
    }

    return(rtc);
}

/*------------------------------------------------------------------------------*/
