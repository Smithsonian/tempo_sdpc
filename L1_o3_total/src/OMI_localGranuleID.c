#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "PGS_TD.h"
#include <cfortran.h>
#include "PGS_OZT_52050.h"
#include "omi_smf.h"

#define FUNCTION_NAME "OMI_localGranuleID()"
#define zero 0

PGSt_SMF_status 
OMI_localGranuleID( char          *RangeBeginningDateTime, 
                    char          *InstPltfm,
                    char          *product_level,
                    char          *ShortName,
                    PGSt_integer   orbitNumber,
                    PGSt_integer   VersionID,
                    char          *filetype,
                    char          *localGranuleID )
{
  PGSt_SMF_status returnStatus;
  /* PGSt_SMF_status code; */         /* status code returned by
                                         PGS_SMF_GetMsg() */
  char msg[PGS_SMF_MAX_MSGBUF_SIZE]={0};
                          /* holds the message string associated with
                             the error code returned by GetMsg */
  char     year[5] = {'\0'},
           month[3]= {'\0'},
           day[4]  = {'\0'},
           hh[3] = {'\0'},
           mm[3] = {'\0'};

  struct tm *time_ptr = NULL;
  time_t     time_sec;

  char     proc_year[5], proc_month[3], proc_day[3],
           proc_hour[3], proc_minute[3], proc_second[3];

  returnStatus = OZT_S_SUCCESS;

  returnStatus = PGS_TD_timeCheck( RangeBeginningDateTime );
  if( returnStatus != PGS_S_SUCCESS )
  {  if( returnStatus == PGSTD_M_ASCII_TIME_FMT_B )
     {  PGS_TD_ASCIItime_BtoA( RangeBeginningDateTime, msg ); 
        sprintf(  RangeBeginningDateTime, "%s", msg );
     }
     else
     {  sprintf(msg, "%s%s", "Error in RangeBeginningDateTime:", 
                      RangeBeginningDateTime );
        OMI_SMF_setmsg( OZT_E_INPUT, msg, FUNCTION_NAME, zero );
     }
  }

  sscanf( RangeBeginningDateTime, "%4s-%2s-%2sT%2s:%2s",
          year, month, day, hh, mm );

  time_sec = time(NULL);
  time_ptr = (struct tm *) gmtime(&time_sec);

  if( time_ptr == NULL )
  {  sprintf(msg, "%s", "gmtime convert to local UTC failed");
     OMI_SMF_setmsg( OZT_E_INPUT, msg, FUNCTION_NAME, zero );
     sprintf( localGranuleID, "NOT_SET" );
     return  OZT_E_INPUT;
  }

  sprintf( proc_year,   "%.4d", 1900+(time_ptr->tm_year) );
  sprintf( proc_month,  "%.2d", time_ptr->tm_mon  + 1    );
  sprintf( proc_day,    "%.2d", time_ptr->tm_mday        );
  sprintf( proc_hour,   "%.2d", time_ptr->tm_hour        );
  sprintf( proc_minute, "%.2d", time_ptr->tm_min         );
  sprintf( proc_second, "%.2d", time_ptr->tm_sec         );

  sprintf( localGranuleID,
           "%s_%s-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%sm%s%st%s%s%s.%s",
           InstPltfm, product_level,
           ShortName, year, month, day, hh, mm,
           orbitNumber, VersionID,
           proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second,
           filetype );

  return OZT_S_SUCCESS; 
}

/* FORTRAN bindings */

FCALLSCFUN8( INT, OMI_localGranuleID, OMI_LOCALGRANULEID, omi_localgranuleid,
             STRING, STRING, STRING, STRING, INT, INT, STRING, PSTRING )
