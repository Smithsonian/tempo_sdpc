/*------------------------------------------------------------------------------*\
** GetConfigInteger.c -- see GetConfig.h for usage
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN3 ( INT, GetConfigInteger, GETCONFIGINTEGER, getconfiginteger,
              STRING, STRING, PINT )

/*------------------------------------------------------------------------------*/

int GetConfigInteger(char* flg, char* key, int* val)
{
    char tmp[CFG_VAL_LEN];

    int ier = GetConfigString(flg, key, tmp);

    if ( ier == 0 ) { *val = atoi(tmp); } else { *val = 0; }

    return(ier);
}

/*------------------------------------------------------------------------------*/
