/*------------------------------------------------------------------------------*\
** GetConfigString.c -- see GetConfig.h for usage
**
** All the other GetConfig_X_() routines use this one.
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN3 ( INT, GetConfigString, GETCONFIGSTRING, getconfigstring,
              STRING, STRING, PSTRING )

/*------------------------------------------------------------------------------*/

extern int LoadControlFile(void);

int GetConfigString(char* flg, char* key, char* val)
{
    int LDa = L8;				/* level of detail		*/
    int LDb = L9;				/* level of detail		*/

    /*--------------------------------------------------------------------------*/

    if ( nry == 0 )				/* if internal table empty...	*/
    {
        Info(LDb, "OK need to read and parse control file");

        int er1 = LoadControlFile();	ReturnIfFatal(er1);
    }

    /*--------------------------------------------------------------------------*/

    int er2 = ymlGetOneStr(nry, cfg, key, 0, val);

    if ( ! *val )
    {
       if ( *flg == 'E' || *flg == 'e' )
       {
          Error("no value for config <%s>", key);
          er2 = -14;
       }
       else
       {
          Warn(LDa, "no value for config <%s>", key);
          er2 = 1;
       }
    }

    Info(LDb, "value for <%s> is <%s>", key, val);

    return(er2);
}

/*------------------------------------------------------------------------------*/
