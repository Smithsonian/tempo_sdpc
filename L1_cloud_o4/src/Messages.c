/*------------------------------------------------------------------------------*\
* Messages.c -- Error Message Display System
*
* See Messages.h or Messages.inc for usage information.
\*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include "Messages.h"

#define OUT stderr

int  max = 2;				/* default "Maximum Level Of Detail"	*/

int  thr = 7;				/* threshold for line number display	*/

char ewibuf[];				/* see include file			*/

/*------------------------------------------------------------------------------*/

int SetMaximumLevelOfDetail(int val) { max = val; }
int GetMaximumLevelOfDetail(int *Lx) { *Lx = max; }

/*------------------------------------------------------------------------------*/

int MsgErr(       char *f, int n, char *m) { MsgCmn("ERROR", 0, f, n, m); }
int MsgWrn(int v, char *f, int n, char *m) { MsgCmn( "WARN", v, f, n, m); }
int MsgInf(int v, char *f, int n, char *m) { MsgCmn(     "", v, f, n, m); }

/*------------------------------------------------------------------------------*\
* t = message type			(string "ERROR", "WARN", "")
* v = message's level of detail 	(integer 0, 1, 2, ...)
* f = filename				(short string)
* n = line number			(integer)
* m = message text			(string)
\*------------------------------------------------------------------------------*/

int MsgCmn(char *t, int v, char *f, int n, char *m)
{
    if ( v <= max )
    {
       ConvertNLto0x0a(m);				/* convert \n to 0x0a	*/
       RemoveNLfromEnd(m);				/* remove  \n from end	*/

       if ( *(m) == ' ' ) { m = m + 1; }		/* remove leading blank	*/

       if ( max > thr )
       {
          fprintf(OUT, "%7s %-72s (%d:%s)\n", t, m, n, f);	/* with name	*/
       } else {
          fprintf(OUT, "%7s %s\n", t, m);			/*  no  name	*/
       }
    }
}

/*------------------------------------------------------------------------------*/

int ConvertNLto0x0a(char *str)				/* convert \n to 0x0a	*/
{
    int p = 0;
    while ( p < EWI_LEN )
    {
        if ( str[p] == '\\' && str[p+1] == 'n' )
        {
           str[p++] = ' '; str[p] = '\n';
        }
        p++;
    }
}

/*------------------------------------------------------------------------------*/

int RemoveNLfromEnd(char *str)				/* remove \n from end	*/
{
    char *BufPtr = strrchr(str, '\n');			/* get ptr to last \n	*/

    if ( BufPtr != 0 && *(BufPtr + 1) == '\0' )		/* if \n at end of str	*/
    {
       *BufPtr = '\0';					/* replace \n with \0	*/
    }
}

/*------------------------------------------------------------------------------*/
/* FORTRAN bindings								*/

#include <cfortran.h>

FCALLSCFUN3( INT, MsgErr, MSGERR, msgerr,      STRING, INT, STRING )
FCALLSCFUN4( INT, MsgWrn, MSGWRN, msgwrn, INT, STRING, INT, STRING )
FCALLSCFUN4( INT, MsgInf, MSGINF, msginf, INT, STRING, INT, STRING )

FCALLSCFUN1( INT, SetMaxLoD, SETMAXLOD, setmaxlod,  INT )
FCALLSCFUN1( INT, SetMaximumLevelOfDetail,
                  SETMAXIMUMLEVELOFDETAIL, setmaximumlevelofdetail,  INT )

FCALLSCFUN1( INT, GetMaxLoD, GETMAXLOD, getmaxlod, PINT )
FCALLSCFUN1( INT, GetMaximumLevelOfDetail,
                  GETMAXIMUMLEVELOFDETAIL, getmaximumlevelofdetail, PINT )

/*------------------------------------------------------------------------------*/
