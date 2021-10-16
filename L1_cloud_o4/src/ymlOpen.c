/*------------------------------------------------------------------------------*\
** ymlOpen.c
**
**   fnm = file name
**   nrp = pointer to number of elements in tbl[]
**   ptr = pointer to table for YAML data
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlOpen(char* fnm, int* nrp, yTB** ptr)
{
    int Lx = L7;				/* level of detail		*/

    char* raw;					/* raw copy of YAML file in mem	*/

    int er1 = ymlRdFile(Lx, fnm, &raw);			ReturnIfFatal(er1);
    int er2 = ymlEstNrE(Lx, raw, nrp);			ReturnIfFatal(er2);
    int er4 = ymlParTxt(Lx, 'E', raw, nrp, ptr);	ReturnIfFatal(er4);

    return(0);
}

/*------------------------------------------------------------------------------*/
