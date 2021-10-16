/*------------------------------------------------------------------------------*\
** GetConfigDouble.c -- see GetConfig.h for usage
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN3 ( INT, GetConfigDouble, GETCONFIGDOUBLE, getconfigdouble,
              STRING, STRING, PDOUBLE )

/*------------------------------------------------------------------------------*/

int GetConfigDouble(char* flg, char* key, double* val)
{
    char tmp[CFG_VAL_LEN];

    int ier = GetConfigString(flg, key, tmp);

    if ( ier == 0 ) { *val = (double) atof(tmp); } else { *val = 0.0; }

    return(ier);
}

/*------------------------------------------------------------------------------*/
