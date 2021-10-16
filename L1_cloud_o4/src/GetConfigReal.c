/*------------------------------------------------------------------------------*\
** GetConfigReal.c -- see GetConfig.h for usage
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN3 ( INT, GetConfigReal, GETCONFIGREAL, getconfigreal,
              STRING, STRING, PFLOAT )

/*------------------------------------------------------------------------------*/

int GetConfigReal(char* flg, char* key, float* val)
{
    char tmp[CFG_VAL_LEN];

    int ier = GetConfigString(flg, key, tmp);

    if ( ier == 0 ) { *val = (float) atof(tmp); } else { *val = 0.0; }

    return(ier);
}

/*------------------------------------------------------------------------------*/
