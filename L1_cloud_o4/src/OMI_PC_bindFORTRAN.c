/*******************************************************************************
BEGIN_FILE_PROLOG:

FILENAME:
  OMI_PC_bindFORTRAN.c

DESCRIPTION:
  This file contains FORTRAN bindings for the OMI Get Reference tool.

AUTHOR:
  Ellyne Kinney / SSAI

HISTORY:
  17-JUL-2003  EKK  Initial version

END_FILE_PROLOG:
*******************************************************************************/
#include <PGS_PC.h>
#include <cfortran.h>
#include <cfortHdf.h>

/*********************
 * cfortran.h MACROS *
 *********************/

/* The following macros are expanded by cfortran.h to C functions which allow
 *    FORTRAN users to access the SDP Toolkit relatively painlessly. */

/* OMI_PC_GetReference() */ 

FCALLSCFUN4(INT, OMI_PC_GetReference, OMI_PC_GETREFERENCE, omi_pc_getreference,\
		            INT, INT, PSTRING, PSTRING)

/* OMI_PC_RefIDsplit() */ 
 
FCALLSCFUN3(INT, OMI_PC_RefIDsplit, OMI_PC_REFIDSPLIT, omi_pc_refidsplit, \
		            PSTRING, PSTRING, PSTRING)
