/*------------------------------------------------------------------------------*\
** SetTrueIfLittleEndian.c
**
** Returns TRUE if our machine is of LittleEndian architecture.
** Returns FALSE if BigEndian.
** Returns UNDEF if error.
** From 'C: A Reference Manual', by Samuel P. Harbison and Guy L. Steele Jr.
\*------------------------------------------------------------------------------*/

#include "O3_PEATE_Common.h"

int SetTrueIfLittleEndian()
{
    union { long l; char c[sizeof(long)]; } u;
    u.l = 1;
    if ( u.c[0] == 1)                    return TRUE;
    else if (u.c[sizeof(long) - 1] == 1) return FALSE;
    else                                 return UNDEF;
}

/*------------------------------------------------------------------------------*/
/* FORTRAN bindings								*/

#include <cfortran.h>

FCALLSCFUN0(INT, SetTrueIfLittleEndian, SETTRUEIFLITTLEENDIAN, settrueiflittleendian)

/*------------------------------------------------------------------------------*/

