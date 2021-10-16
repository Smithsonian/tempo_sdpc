/*------------------------------------------------------------------------------*\
** ymlGetOneStr.c
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlGetOneStr(int nry, yTB tbl[], char* key, int cnt, char* val)
{
    int LDa = L7;				/* level of detail		*/

    if ( cnt < 0 ) { Error("ymlGetOneStr() count %d may not be negative", cnt);
       return(-13);
    }

    int rtc = 0;				/* return code			*/

    if ( cnt > nry )
    {
       /* although this is an error, it is not fatal; always display (i.e., L0)	*/

       Warn(L0, "ymlGetOneStr() count %d > %d", cnt, nry);
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

    *val = '\0';			/* initialize returned value		*/

    int ix2; for( ix2 = 0; ix2 < nry; ix2 ++ )		/* loop tbl[]		*/
    {
       if ( tbl[ix2].key[0] == '\0' ) { break; }	/* done; end of tbl[]	*/

       if ( strncmp(tmp, tbl[ix2].key, len) == 0 )	/* found one		*/
       {
          hit = 1;					/* key exists		*/
          num ++;					/* count of finds	*/

          Info(LDa, "key <%s> found; used = %d", tbl[ix2].key, tbl[ix2].usd);

          if ( cnt == num || ( ! cnt && ! tbl[ix2].usd ) )	/* it is it!	*/
          {
             strncpy(val, tbl[ix2].val, YML_VAL_LEN);		/* grab it	*/
             val[YML_VAL_LEN-1] = '\0';				/* if too long	*/
             tbl[ix2].usd = 1;					/* row now used	*/
             ixo = 1;						/* value found	*/
             break;						/* done		*/
          }
       }
    }

    if ( hit == 0 ) { return(-12); }	/* key not found			*/
    if ( ixo == 0 ) { return(  1); }	/* key found, but there is no value	*/
                      return(rtc);	/* value found				*/
}

/*------------------------------------------------------------------------------*/
