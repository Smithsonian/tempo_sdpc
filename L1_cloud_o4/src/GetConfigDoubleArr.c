/*------------------------------------------------------------------------------*\
** GetConfigDoubleArr.c -- see GetConfig.h for usage
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN4 ( INT, GetConfigDoubleArr, GETCONFIGDOUBLEARR, getconfigdoublearr,
              STRING, STRING, PDOUBLE, INT )

/*------------------------------------------------------------------------------*/

int GetConfigDoubleArr(char* flg, char* key, double* vls, const int cnt)
{
    char tmp[CFG_VAL_LEN];

    int ier = 0;

    int i; for ( i = 0; i < cnt; ++i )
    {
        int err = GetConfigString(flg, key, tmp);
        if ( err == 0 ) { vls[i] = atof(tmp);      }
        else            { vls[i] = 0.0; ier = err; }
    }

    return(ier);
}

/*------------------------------------------------------------------------------*/
