/*------------------------------------------------------------------------------*\
** ymlMaxDtl.c
**
** Sets Maximum Level Of Detail
**
**   LDx = level of detail
**   key = parameter name
**   val = parameter value
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlMaxDtl(int LDx, char* key, char* val)
{
    int rtc = 0;					/* return code		*/

    if ( strcmp(key, "RUNTIMEPARAMETERSMAXIMUMLEVELOFDETAIL") == 0 )
    {
       int num = atoi(val);

       if ( num == 0 && *val != '0' )
       {
          Error("Maximum Level Of Detail <%s> not an integer", val);
          rtc = 9;
       }
       else
       {
          Info(LDx, "OK ++++++++++++++++++++ setting level of detail to %d", num);
          SetMaximumLevelOfDetail(num);
       }
    }

    return(rtc);
}

/*------------------------------------------------------------------------------*/
