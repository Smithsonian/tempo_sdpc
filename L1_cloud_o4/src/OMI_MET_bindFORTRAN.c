/*******************************************************************************
BEGIN_FILE_PROLOG:

FILENAME:
  OMI_MET_bindFORTRAN.c

DESCRIPTION:
  This file contains FORTRAN bindings for the OMI Toolkit Metadata (MET) tools.

AUTHOR:
  Ellyne Kinney / SSAI

HISTORY:
  24-APR-2003  EKK  Initial version

END_FILE_PROLOG:
*******************************************************************************/
#include <PGS_MET.h>
#include <cfortran.h>
#include <cfortHdf.h>

/*********************
 * cfortran.h MACROS *
 *********************/

/* The following macros are expanded by cfortran.h to C functions which allow
 *    FORTRAN users to access the SDP Toolkit relatively painlessly. */

/* OMI_LocalGranuleID() */ 
 
#define omi_localgranuleid_STRV_A2 NUM_ELEMS(PGSd_MET_NUM_OF_GROUPS)

FCALLSCFUN2(INT, OMI_localGranuleID, OMIR_LOCALGRANULEID, omi_localgranuleid, \
	    INT, PSTRINGV)

