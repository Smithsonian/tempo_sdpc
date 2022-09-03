/*------------------------------------------------------------------------------*\
** GetConfigRealArr.c
**
** Parses a list of reals from an ASCII control file.
** It is up to the user of this function to ensure proper allocation of vls[].
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN4 ( INT, GetConfigRealArr, GETCONFIGREALARR, getconfigrealarr,
              STRING, STRING, PFLOAT, INT )

/*------------------------------------------------------------------------------*/

int GetConfigRealArr (char* flg, char* key, float* vls, const int cnt)
{

    char val_str[CFG_VAL_LEN];
    int  err=0, i;

    for ( i = 0 ; i < cnt; ++i )
    {
        err = GetConfigString(flg, key, val_str);

        if ( ! err ) { vls[i] = atof(val_str); }
    }

    return(err);
}

/*------------------------------------------------------------------------------*/
