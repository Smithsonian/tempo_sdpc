/*------------------------------------------------------------------------------*\
** LoadControlFile.c
**
** Allocates memory for the control file, reads it, and YAML-parses it.
**
** This routine is not intended to be called directly by programmers; call
** the GetConfig_X_() routines instead.  See GetConfig.h.
**
** Use MetaConfig() to write the control file content into H5 files.
**
** The control file is also made available to users as a single null-terminated
** string via the global ControlFile variable.  Here is how another C function
** could use it:
**
**   #include "GetConfig.h"
**   printf("%s", ControlFile);
**
** Please note that the ControlFile variable is not initialized until the user
** makes a call to any of the GetConfig_X_() functions.  Once initialized, the
** memory is allocated using malloc(), and never freed; that's on purpose so
** that the user may access the allocated information at any time, right
** through the main program's termination.
\*------------------------------------------------------------------------------*/

#include "Config.h"

char* ControlFile;			/* control file (null-terminated str)	*/

yTB* cfg;				/* CONFIG_READER main data structure	*/
int  nry = 0;				/* number of elements in cfg[]		*/

/*------------------------------------------------------------------------------*/

int LoadControlFile(void)
{
    int LDa = L7;				/* level of detail		*/

    /*--------------------------------------------------------------------------*/
    /* read control file name from command line					*/

    int  num = 1;				/* arg number on command line	*/
    int  siz = VAL_LEN;				/* control file name max size	*/
    char fnm[C_VAL_LEN];			/* control file name buffer	*/

    int i; for ( i = 0; i <= siz; i++ ) { fnm[i] = '\0'; }	/* name init	*/

    /*--------------------------------------------------------------------------*/

#ifdef pgf
    getarg_(&num, &fnm, siz);			/* control file name retrieved	*/
#endif

#ifdef ifort
    getarg_(&num, &fnm, siz);			/* control file name retrieved	*/
#endif

#ifdef g77
    getarg_(&num, &fnm, siz);			/* control file name retrieved	*/
#endif

#ifdef gfortran
    _gfortran_getarg_i4(&num, &fnm, siz);	/* control file name retrieved	*/
#endif

    /*--------------------------------------------------------------------------*/

    char *ptr = strechr(fnm, ' ');			/* find string end	*/

    if ( ptr ) { *ptr = '\0'; } else { *fnm = '\0'; }	/* zap trailing blanks	*/

    Info(LDa, "OK control file name: <%s>", fnm);

    /*--------------------------------------------------------------------------*/

    char* raw;						/* copy of text in mem	*/

    int er1 = ymlRdFile(LDa, fnm, &raw);		ReturnIfFatal(er1);

    ControlFile = raw;					/* pointer for users	*/

    Info(LDa, "OK read control file <%s>", fnm);

    /*--------------------------------------------------------------------------*/

    int er2 = ymlEstNrE(LDa, raw, &nry);		ReturnIfFatal(er2);

    Info(LDa, "OK estimated %d elements in cfg[]", nry);

    int er3 = ymlIniTbl(LDa, &nry, &cfg);		ReturnIfFatal(er3);

    Info(LDa, "OK initialized YAML parser");

    int er4 = ymlParTxt(LDa, 'E', raw, &nry, &cfg);	ReturnIfFatal(er4);

    Info(LDa, "OK YAML-parsed <%s>", fnm);

    return er4;
}

/*------------------------------------------------------------------------------*/
