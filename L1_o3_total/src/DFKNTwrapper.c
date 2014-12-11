#include <mfhdf.h>
#include <cfortHdf.h>

int
DFKNTwrapper( int32  datatype )
/*****************************************************************************
!C

!Description:

!Input Parameters:
  datatype 
!Output Parameters:
  None
  
!Return
  size in bytes for the input data type on successful return
  -1      when unknow datatype encountered 
 
!Revision History:
  Revision 0.1  11/26/2001  Kai Yang/UMBC
  Adopted smf.c from MODIS code.

!Team-unique Header:
  This software was developed by the OMI Science Team Support
  Group for the National Aeronautics and Space Administration, Goddard
  Space Flight Center, under NASA Task 916-003-1

!References and Credits
  Written by 
  Kai Yang 
  University of Maryland Baltimore County
  email: Kai.Yang-1@lnasa.gov
  
!Design Notes

!END
*****************************************************************************/
{
  int size;

  size = DFKNTsize( datatype );
  return size;
}

/* FORTRAN bindings */

FCALLSCFUN1( INT, DFKNTwrapper, DFKNTSIZE, dfkntsize, INT )
