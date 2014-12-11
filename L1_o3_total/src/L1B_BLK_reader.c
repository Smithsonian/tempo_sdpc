#include "PGS_OZT_52050.h"
#include "l1b_reader.h"
#include "hdf.h"
#include "cfortHdf.h"
#include "HdfEosDef.h"
#include "omi_smf.h"

#define FUNCTION_NAME "Get_L1B_Data_Block()"
#define MAX_LEN_DIM_LIST 256
#define MAX_RANK         4

PGSt_SMF_status
Get_L1B_Data_Block( char *filename,       /* input  */
                    char *swathname,      /* input  */
                    char *fieldname,      /* input  */
                    int   numberType,     /* input  */
                    int   Line_start,     /* input  */
                    int   nLines,         /* input  */
                    int  *rank,           /* output */
                    int  *dims,           /* output */
                    void *data_buffer )   /* output */
/*****************************************************************************
!C

!Description:
  Get_L1B_Data_Block reads from a swath in an HE4 file a number of
  lines of data in the datafield. Here line refers to first dimension for
  the datafield in the HE4 files, and index in HDF file is 0 based. 

!Input Parameters:
  char *filename    HE4 filename
  char *swathname   swath name in the HE4 file
  char *fieldname   data set name in the swath
  int   numberType  data type of the data set
  int   Line_start  starting index (0 based) of the first dimension of
                    this data set wheren the data will be retrieved
                    from the swath in the HE4 file 
  int   nLines      the number of lines of data will be retrieved

!Output Parameters:
  int  *rank,       the rank, i.e., the number of dimensions of this data set
  int  *dims,       the extent of each dimension
  void *data_buffer pointer to the memory location where the nLines of data 
                    is stored.
!Return
   OZT_S_SUCCESS  successful return
   OZT_E_FAILURE  error encountered

!Revision History:
  Revision 0.1  11/26/2001  Kai Yang/UMBC

!Team-unique Header:
  This software was developed by the OMI Science Team Support
  Group for the National Aeronautics and Space Administration, Goddard
  Space Flight Center, under NASA Task 916-003-1

!References and Credits
  Written by
  Kai Yang
  University of Maryland Baltimore County 
  email: Kai.Yang-1@nasa.gov

!Design Notes

  A fortran binding of this function is provided here. Fortran code
  can link and use this function directly.

  This functions can read a data set with a maximum rank of 4. The first
  dimension is refered to as line which starts from 0 in this function.

!END
*****************************************************************************/
{  int32  swfid_he4 = FAIL;
   int32  SWid_he4  = FAIL;
   intn   status_he4= FAIL; 
   int32  numberType_local;
   int32  rank_local;
   int32  dims_local[MAX_RANK];
   int32  start[MAX_RANK] = {0,0,0,0};
   int32  edge[MAX_RANK];
   int    di;
   static PGSt_integer s_code = 0;
   char   dimlist[MAX_LEN_DIM_LIST] = {'\0'};
   char   msg[PGS_SMF_MAX_MSG_SIZE];  /* PGS_SMF_MAX_MSG_SIZE is
                                         defined in PGS_SMF.h */
   /* check input names to make sure they are not empty */
   if( filename == NULL )
   {  OMI_SMF_setmsg( OZT_E_INPUT, "empty filename", FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }
   else if( swathname == NULL )
   {  OMI_SMF_setmsg( OZT_E_INPUT, "empty swathname", FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }
   else if( fieldname == NULL )
   {  OMI_SMF_setmsg( OZT_E_INPUT, "empty fieldname", FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   if( Line_start < 0 || nLines < 0 )
   {  sprintf( msg, "Line_start = %d, nLines = %d.", Line_start, nLines );
      OMI_SMF_setmsg( OZT_E_VALID_RANGE, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }
   
   /* open the swath file for read */
   swfid_he4 = SWopen( filename, DFACC_READ );
   if( swfid_he4 == FAIL )
   {  sprintf( msg, "open file %s for read failed.", filename );
      OMI_SMF_setmsg( OZT_E_FILE_OPEN, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   /* attach the swath */
   SWid_he4 = SWattach( swfid_he4, swathname );
   if( SWid_he4 == FAIL )
   {  sprintf( msg, "attach swath %s in file %s failed.", swathname, filename );
      OMI_SMF_setmsg( OZT_E_SWATH_ATTACH, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   status_he4 = SWfieldinfo( SWid_he4, fieldname, 
                            &rank_local, dims_local,
                            &numberType_local, dimlist );
   if( status_he4 == FAIL )
   {  sprintf( msg, 
              "%s: retrieve field info for %s failed in swath %s in file %s",
              "SWfieldinfo()", fieldname, swathname, filename );
      OMI_SMF_setmsg( OZT_E_HDFEOS, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   if( numberType_local != (int32) numberType )
   {  sprintf( msg, "pre-set number type different from that in L1B file." ); 
      OMI_SMF_setmsg( OZT_E_INPUT, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   *rank = rank_local;

   for( di = 0; di < rank_local; di++ )
   {  dims[di] = (int) dims_local[di];
      edge[di] = (int) dims_local[di];
   }

   if( Line_start >= dims[0] || (Line_start + nLines) > dims[0] )
   {  sprintf( msg, "Line_start = %d, Line_start+ nLines = %d, dims[0] = %d.", 
               Line_start, Line_start + nLines, dims[0] );
      OMI_SMF_setmsg( OZT_E_VALID_RANGE, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   if( nLines > 0 )
   {  start[0] = Line_start;
      edge[0]  = nLines;

      status_he4 = SWreadfield( SWid_he4, fieldname, 
                                start, NULL, edge,
                                (VOIDP) data_buffer );
      if( status_he4 == FAIL )
      {  sprintf( msg, 
              "%s: read field data for %s failed in swath %s in file %s",
              "SWreadfield()", fieldname, swathname, filename );
         OMI_SMF_setmsg( OZT_E_HDFEOS, msg, FUNCTION_NAME, s_code );
         return OZT_E_FAILURE;
   }  }

   /* detach the attached swath */
   if( SWdetach( SWid_he4 ) == FAIL )
   {  sprintf( msg, "detach swath %s in file %s failed.", swathname, filename );
      OMI_SMF_setmsg( OZT_E_SWATH_ATTACH, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }

   /* close the opened swath file */
   if( SWclose( swfid_he4 ) == FAIL )
   {  sprintf( msg, "close file %s failed.", filename );
      OMI_SMF_setmsg( OZT_E_FILE_OPEN, msg, FUNCTION_NAME, s_code );
      return OZT_E_FAILURE;
   }
   return OZT_S_SUCCESS;
}

/* HDF types used in FORTRAN bindings */

#if defined(DEC_ALPHA) || defined(IRIX) || defined(UNICOS)

#define INT32  INT
#define INT32V INTV
#define PINT32 PINT

#else

#define INT32  LONG
#define INT32V LONGV
#define PINT32 PLONG

#endif

/* FORTRAN bindings */

FCALLSCFUN9(INT, Get_L1B_Data_Block, GETL1BDBLK, getl1bdblk, STRING, STRING,
            STRING, INT, INT, INT, PINT, INTV, PVOID )
