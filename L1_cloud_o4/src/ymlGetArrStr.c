/*------------------------------------------------------------------------------*\
** ymlGetArrStr.c
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlGetArrStr(int nry, yTB tbl[], char* key, int cnt, int inb,
				int* blk, int arm, char arr[arm][YML_VAL_LEN])
{
    int LDa = L7;				/* level of detail		*/

    if ( cnt < 0 ) { Error("ymlGetArrStr() count %d may not be negative", cnt);
       return(-13);
    }

    int rtc = 0;				/* return code			*/

    if ( cnt > nry )
    {
       /* although this is an error, it is not fatal; always display (i.e., L0)	*/

       Warn(L0, "ymlGetArrStr() count %d > %d", cnt, nry);
       Info(L0, "Searching past the table's end should not return anything.");
       rtc = 6;
    }

    char tmp[YML_KEY_LEN]; size_t len = YML_KEY_LEN; ymlNrmKey(LDa, key, tmp);

    if ( tbl[0].key[0] == '\0' )		/* tbl[] not initialized	*/
    {
        Error("tbl[] has not been initialized; call ymlOpen() first");
        return(-10);
    }

    int num = 0;			/* counter of matching keys found	*/
    int hit = 0;			/* flag set when key exists		*/
    int ixo = 0;			/* number of values found		*/

    int ht1 = 0;		/* flag set when array found, count not reached	*/
    int ht2 = 0;		/* flag set when array found and count reached	*/
    int cur = INT_MIN;		/* block number we're reading			*/

    int ix1; for( ix1 = 0; ix1 < arm; ix1++ ) { arr[ix1][0] = '\0'; }   /* init */

    int ix2; for( ix2 = 0; ix2 < nry; ix2 ++ )		/* searching tbl[]	*/
    {
       if ( tbl[ix2].key[0] == '\0' ) { break; }	/* done; end of tbl[]	*/

       if ( strncmp(tmp, tbl[ix2].key, len) == 0 )	/* found one		*/
       {
          hit = 1;					/* key exists		*/

          Info(LDa, "key <%s> found; used = %d", tbl[ix2].key, tbl[ix2].usd);
       }

       if ( tbl[ix2].blk <     inb  ) { continue; }	/* skip: block nr low	*/
       if ( tbl[ix2].blk == cur + 1 ) { break;    }	/* next block started	*/

       if ( strncmp(tmp, tbl[ix2].key, len) == 0 )	/* found one		*/
       {
          if ( ht1 == 0 ) { cur = tbl[ix2].blk; num++; }   /* catch block nr	*/

          ht1 = 1;					/* may not be the one	*/

          if ( num == cnt || ( ! cnt && ! tbl[ix2].usd ) ) { ht2 = 1; }	/* yes!	*/

          if ( ht2 == 1 )				/* we want it		*/
          {
             if ( ixo < arm )				/* prevents overflow	*/
             {
                strncpy(arr[ixo], tbl[ix2].val, YML_VAL_LEN);	/* grab it	*/
                arr[ixo][YML_VAL_LEN-1] = '\0';			/* if too long	*/
             }
             ixo ++;					/* count past overflow	*/
          }
       }
       else						/* not of interest	*/
       {
          if ( ht2 == 1 ) { break; }			/* done			*/
               ht1  = 0;				/* reset flag		*/
       }
    }

    if ( ixo > arm )
    {
       Error("arr[%d] overflow getting <%s>; increase to %d", arm, key, ixo);
       rtc = 7;
       /* overflow is not fatal since arr[] was truncated; no damage was done	*/
    }

    *blk = cur;					/* block number			*/

    if ( hit == 0 ) { return(-12); }	/* key not found			*/
    if ( ixo == 0 ) { return(  1); }	/* key found, but there is no value	*/
                      return(rtc);	/* at lease one value found		*/
}

/*------------------------------------------------------------------------------*/
