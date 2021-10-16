/*------------------------------------------------------------------------------*\
** ymlRdFile.c
**
** This routine does the following:
** - determines the file's size
** - allocates that much memory
** - opens the file
** - reads the file into the allocated memory
** - closes the file
** - returns a pointer to the allocated and populated memory
**
**   LDx = level of detail
**   fnm = file name
**   mem = pointer to memory (returned)
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

char* svp = NULL;				/* saved pointer for ymlClose()	*/

void ymlClose(void) { free(svp); svp = NULL; }	/* release allocated memory	*/

/*------------------------------------------------------------------------------*/

int ymlRdFile(int LDx, char* fnm, char** mem)
{
    int    mss = EWI_LEN;			/* error msg[] buffer size	*/
    char   msg  [EWI_LEN];			/* error or progress message	*/

    /*--------------------------------------------------------------------------*/

    size_t stp = YML_KEY_LEN - 1;		/* one less due to terminator	*/

    char   fn2[YML_KEY_LEN]; strncpy(fn2, fnm, stp); fn2[stp] = '\0';

    char*  fn3 = strxchr(fn2, ' ');		/* find string's start		*/
    char*  end = strechr(fn2, ' ');		/* find string's end		*/

    if ( ! fn3 ) { fn3 = fn2; }			/* start not found; set to init	*/
    if ( ! end ) { end = fn2; }			/* end   not found; set to init	*/

          *end =  '\0';				/* remove trailing blanks	*/

    if (  *fn2 == '\0' ) { Error("no file; file name is blank"); return(-1); }

    Info(LDx, "OK file name <%s>", fn3);

    /*--------------------------------------------------------------------------*/

    snprintf(msg, mss, "getting size of <%s>", fn3);

    struct stat stt; int rc1 = stat(fn3, &stt); if ( rc1 < 0 ) { Err(-19);} Err(0);

    size_t fsz = stt.st_size + 1;		/* plus one for terminator	*/

    /*--------------------------------------------------------------------------*/

    snprintf(msg, mss, "allocating %d bytes of memory for file", fsz);

    *mem = malloc(fsz); if ( *mem == NULL ) { Err(-2);} Err(0);

    /*--------------------------------------------------------------------------*/

    snprintf(msg, mss, "scheduling a call to free memory at exit");

    svp = *mem;					/* ptr needed by ymlClose()	*/

    int rc2 = atexit(ymlClose); if ( rc2 != 0 ) { Err(-11);} Err(0);

    /*--------------------------------------------------------------------------*/

    snprintf(msg, mss, "opening <%s>", fn3);

    FILE* fid = fopen(fn3, "r");

    if ( fid == NULL )
    {
       Error("errno = %d: %s", errno, strerror(errno) );

       Err(-3);			/* errno is often positive, but this is fatal	*/
    }

    Err(0);

    /*--------------------------------------------------------------------------*/

    size_t len = fread(*mem, 1, stt.st_size, fid);

    if ( len != stt.st_size )
    {
       Error("%s is %d bytes != %d bytes read", fn3, stt.st_size, len);
       return(-4);
    }

    if ( len  < 0 ) { Error("reading <%s>", fn3); return(-5); }

    if ( len == 0 ) { Error("<%s> is empty", fn3); return(-6); }

    Info(LDx, "OK read %d bytes from <%s>", len, fn3);

    (*mem)[stt.st_size] = '\0';			/* terminate			*/

    /*--------------------------------------------------------------------------*/

    snprintf(msg, mss, "closing <%s>", fn3); int rc3 = fclose(fid); Err(rc3);

    /*--------------------------------------------------------------------------*/

    return(0);
}

/*------------------------------------------------------------------------------*/
