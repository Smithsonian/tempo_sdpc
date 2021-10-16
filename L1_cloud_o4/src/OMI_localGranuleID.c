/*******************************************************************************
BEGIN_FILE_PROLOG:

FILENAME:
        OMI_localGranuleID

DESCRIPTION:
	Reads the  following metadata parameters from the inventory metadata:
	SHORTNAME, RANGEBEGINDATE. RANGEBEGINTIME, VERSIONID, & ORBITNUMBER.
	Constructs the LocalGranuleID from these metadata and from the 
	processing level.  Writes the LocalGranuleID metadata parameter 
	to the inventory metadata.

AUTHOR:
        Ellyne Kinney / SSAI

HISTORY:
        24-APR-03      EKK     Initial version
        15-JUN-06      EKK     Added L2G to processing level

END_FILE_PROLOG
*******************************************************************************/

/* include files */

#include <string.h>
#include <time.h>
#include <PGS_TD.h>
#include "omi_smf.h" /* include PGS_PC.h PGS_SMF.h */
#include "OMI_MET_Tools.h"
#include "PGS_OMI_52090.h"

/***************************************************************************
BEGIN_PROLOG:

TITLE:
        Populates the LOCALGRANULEID metadata parameter 

NAME:
        OMI_localGranuleID()

SYNOPSIS:
C:
        #include "PGS_MET.h"

        PGSt_SMF_status
        OMI_localGranuleID( enum processing_level ProcessingLevel,
                            PGSt_MET_all_handles mdHandles)

FORTRAN:
        include "PGS_MET_13.f"
        include "PGS_MET.f"
        include "PGS_SMF.h"

        integer function omi_localgranuleid(ProcessingLevel, mdHandles)

        integer          ProcessingLevel 
        character*(49)   mdHandles(20)

DESCRIPTION:
	Reads the  following metadata parameters from the inventory metadata:
	SHORTNAME, RANGEBEGINDATE. RANGEBEGINTIME, VERSIONID, & ORBITNUMBER.
	Constructs the LocalGranuleID from these metadata and from the 
	processing level.  Writes the LocalGranuleID metadata parameter 
	to the inventory metadata.

        Processing levels: L1  = 1
                           L2  = 2
			   L3  = 3
			   L2A = 4 
			   L2B = 5
			   L2B = 6

INPUTS:
        Name            Description            Units           Min          Max
        ----            -----------            -----           ---          ---
        ProcessingLevel Processing level       none             1            6

        mdHandles       metadata groups        none            N/A          N/A
                        in MCF
			
OUTPUTS:
        None

RETURNS:
        OMI_S_SUCCESS
	OMI_E_GENERAL            General error message
	OMI_E_MEM_ALLOC          Memory allocation error
	OMI_E_MET_READ           Metadata read failure
	OMI_E_MET_WRITE          Metadata write failure
	OMI_E_GRAN_ID            Local Granule ID could not be set

EXAMPLES:
C:
        #include <PGS_SMF.h>
        #include <PGS_MET.h>

        #define MCF_FILE       209040

	PGSt_MET_all_handles mdhandles;
	PGSt_integer         HE4fid;

	/
         * Read Metadata Configuration File (MCF) file into memory
	 * The syntax of the MCF is checked
	 * The SHORTNAME and VERSIONID metadata parameters are defined in
	 * the MCF.  
         *
         * The ORBITNUMBER metadata parameter is defined the the PCF.
         * PCF file must define the ORBITNUMBER and SMF_VERBOSITY parameters
         * in the USER DEFINED RUNTIME PARAMETERS section as follows:
         *         200105|ORBITNUMBER|10001
         *         200100|SMF_VERBOSITY|1
	 *
	 * MCF_FILE is the LUN of the MCF in the PCF.
	 * mdhandles is the structure that contains the Metadata handles 
	  /

        PGS_MET_Init( MCF_FILE, mdhandles );

	/ 
         * For objects in the MCF with data location of PGE use the following
	 * functions to set the value for Metadata determined by the PGE.
         *
         * mdhandles[1] designates the metadata should be written to the 
	 * inventory metadata handle.
          /

	PGS_MET_SetAttr( mdhandles[1], "RangeBeginningDate", "2002-01-01");
        PGS_MET_SetAttr( mdhandles[1], "RangeBeginningTime", "01:01:59.9" );

	OMI_localGranuleID( L1, mdhandles);

	/ 
        *	Get the metadata ID  
	 /

	PGS_MET_SDstart( "SampleHDF.he4", HDF4_ACC_RDWR, &HE4fid );

	/
         * Write the the metadata
	 * mdhandles[1] is the metadata handle
         * 'CoreMetadata.0' is the name of the metadata
          /

	returnStatus = PGS_MET_Write( mdhandles[1], "CoreMetadata.0", HE4fid );

FORTRAN:
        include "PGS_SMF.f"
        include "PGS_MET_13.f"
        include "PGS_MET.f"

	include "PGS_OMI_52090.f"	!OMI error messages

	!  File Access Flags for FORTRAN
	INTEGER (KIND = 4), PARAMETER :: HE5_HDFE_RDWR   = 0
	INTEGER (KIND = 4), PARAMETER :: HE5_HDFE_RDONLY = 1
	INTEGER (KIND = 4), PARAMETER :: HE5_HDFE_TRUNC  = 2

        !  LUN of the MCF file and output file
        INTEGER, PARAMETER :: MCF_FILE = 10241

        !  The string that holds the Metadata handles
        CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: GROUPS(PGSd_MET_NUM_OF_GROUPS)

	!  Types of Metadata
        INTEGER, PARAMETER :: INVENTORY = 2
        INTEGER, PARAMETER :: ARCHIVE   = 3

        !  Metadata ID
        INTEGER (KIND = 4) :: metID

	!  Functions called by this program
        INTEGER ::pgs_met_init, pgs_met_setattr_s, omi_localgranuleid
        INTEGER ::pgs_met_sfstart, pgs_met_write, pgs_met_remove

        !  Read Metadata Configuration File (MCF) file into memory
        !  The syntax of the MCF is checked
        !  The SHORTNAME and VERSIONID metadata parameters are defined in 
	!  the MCF.  

        !  The ORBITNUMBER metadata parameter is defined the the PCF.
        !  PCF file must define the ORBITNUMBER and SMF_VERBOSITY parameters
        !  in the USER DEFINED RUNTIME PARAMETERS section as follows:
        !          200105|ORBITNUMBER|10001
        !          200100|SMF_VERBOSITY|1
        !
	!  MCF_FILE is the LUN of the MCF in the PCF.
        !  GROUPS is the structure that contains the Metadata handles
        
        status = pgs_met_init(MCF_FILE,GROUPS)

	!  For Objects in the MCF with data location of PGE use the following
	!  functions to set the value for Metadata determined by the PGE.
	!  The pgs_met_setattr_s function sets values with the data type of
	!  character string.
	!
	!  GROUPS(INVENTORY) designates to which metadata handle the metadata
	!  should be written.
	!

	status = pgs_met_setattr_s(GROUPS(INVENTORY),"RANGEBEGINNINGDATE", \
				   "2002-01-01")
	status = pgs_met_setattr_s(GROUPS(INVENTORY),"RANGEBEGINNINGTIME", \
		       	 	   "01:01:59.9")

	!  Construct and write the localgranuleid to the inventory metadata

	status = omi_localgranuleid(2,GROUPS);

	!  Get the metadata ID

	status =  pgs_met_sfstart("SampleHDF.he5",HDF5_ACC_RDWR,metID)

	!  Write the the metadata
	!  groups(INVENTORY) is the metadata handle
	!  'coremetadata' is the name of the metadata

	status = pgs_met_write(GROUPS(INVENTORY),'coremetadata',metID)

NOTES:
        MCF file must be in the format described in the MET userguide
        
DETAILS:
        This routine reads the SHORTNAME, RANGEBEGINDATE, RANGEBEGINTIME,
	VERSIONID and ORBITNUMBER metadata parameters from the inventory
        metadata.  The localGranuleID string is then constructed with these 
	values and the ProcessingLevel argument.  The localGranuleID string
	in then written to the localGranuleID metadata parameter.
	The SHORTNAME and VERSIONID are defined the the MCF and the ORBITNUMBER
	is defined the the PCF.  The PGE must call the PGS_MET_Init() function 
	and set the RANGEBEGINDATE and RANGEBEGINTIME metadata parameters 
	before calling the OMI_localGranuleID() function.

GLOBALS:
        None  

FILES:
	MCF, PCF

FUNCTIONS_CALLED:
        OMI_SMF_setmsg
	PGS_SMF_GetMsg
        PGS_MET_GetSetAttr
	PGS_MET_SetAttr

END_PROLOG:
***************************************************************************/

#define FUNCTION_NAME "OMI_localGranuleID()"

PGSt_SMF_status 
OMI_localGranuleID( enum processing_level ProcessingLevel, 
                    PGSt_MET_all_handles  mdhandles)
{
  char  *ShortName;
  char  *RangeBeginDate;
  char  *RangeBeginTime;
  int    OrbitNumber[1] = {0};
  char   VID[PGSd_PC_VALUE_LENGTH_MAX] = "-1\0";
  int    VersionID = -1;
  int    VersionID_LUN = 200205;

  char     year[5] = "\0", month[3] = "\0", day[3] = "\0";
  char     hour[3] = "\0", minute[3] = "\0";
		
  time_t time_sec;
  struct tm
        *time_ptr = NULL;

  char   proc_time[PGSd_PC_VALUE_LENGTH_MAX] = "\0",
         proc_year[PGSd_PC_VALUE_LENGTH_MAX] = "\0",
         proc_month[PGSd_PC_VALUE_LENGTH_MAX] = "\0",
         proc_day[PGSd_PC_VALUE_LENGTH_MAX] = "\0",
         proc_hour[PGSd_PC_VALUE_LENGTH_MAX] = "\0",
         proc_minute[PGSd_PC_VALUE_LENGTH_MAX] = "\0",
         proc_second[PGSd_PC_VALUE_LENGTH_MAX] = "\0";

  char   localGranuleID[PGSd_PC_FILE_PATH_MAX]; /* granule Id */
  void  *local; 

  int version;
  PGSt_SMF_status      returnStatus = PGS_S_SUCCESS;
  char                 mnemonic[PGS_SMF_MAX_MNEMONIC_SIZE];
                                                   /* status mnemonic
                                                      returned by
                                                      PGS_SMF_GetMsg() */
  char                 msg[PGS_SMF_MAX_MSG_SIZE];  /* status messsage
                                                      returned by call to
                                                      PGS_SMF_GetMsg() */
  PGSt_SMF_status      code;                       /* status code returned
                                                      by PGS_SMF_GetMsg() */

  /* Get SHORTNAME from Metadata */
  ShortName = (void *)NULL;
  ShortName = malloc(sizeof(char) * PGSd_PC_VALUE_LENGTH_MAX);

  if (ShortName == NULL)
  {  sprintf(msg, "unable to allocate memory: ShortName\n\n");
     OMI_SMF_setmsg(OMI_E_MEM_ALLOC,msg,FUNCTION_NAME, zero);
     return OMI_E_MEM_ALLOC;
  }   

  returnStatus = PGS_MET_GetSetAttr(mdhandles[1], "SHORTNAME", &ShortName);

  if( returnStatus != PGS_S_SUCCESS)
  {  PGS_SMF_GetMsg(&code,mnemonic,msg);
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  } else {
     sprintf(msg,"ShortName is: %s\n", ShortName);
     OMI_SMF_setmsg(OMI_S_MET_READ,msg,FUNCTION_NAME, two);
  }

  /* Get RANGEBEGINNINGDATE from Metadata */
  RangeBeginDate = (void *)NULL;
  RangeBeginDate = malloc(sizeof(char) * PGSd_PC_VALUE_LENGTH_MAX);

  if (RangeBeginDate == NULL)
  {  sprintf(msg, "unable to allocate memory: RangeBeginDate\n\n");
     OMI_SMF_setmsg(OMI_E_MEM_ALLOC,msg,FUNCTION_NAME, zero);
     return OMI_E_MEM_ALLOC;
  }   

  returnStatus = PGS_MET_GetSetAttr(mdhandles[1], "RANGEBEGINNINGDATE", 
		  &RangeBeginDate);

  if( returnStatus != PGS_S_SUCCESS)
  {  PGS_SMF_GetMsg(&code,mnemonic,msg);
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  }

  returnStatus = sscanf( RangeBeginDate, "%4s-%2s-%2s", year, month, day);

  if (returnStatus != 3)
  {  sprintf(msg, "unable to read: RangeBeginDate\n\n");
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  } else {
     sprintf(msg, "RangeBeginDate is: %4s-%2s-%2s\n", year, month, day);
     OMI_SMF_setmsg(OMI_S_MET_READ,msg,FUNCTION_NAME, two);
  }

  /* Get RANGEBEGINNINGTIME from Metadata */
  RangeBeginTime = (void *)NULL;
  RangeBeginTime = malloc(sizeof(char) * PGSd_PC_VALUE_LENGTH_MAX);

  if (RangeBeginTime == NULL)
  {  sprintf(msg, "unable to allocate memory: RangeBeginDate\n\n");
     OMI_SMF_setmsg(OMI_E_MEM_ALLOC,msg,FUNCTION_NAME, zero);
     return OMI_E_MEM_ALLOC;
  }   

  returnStatus = PGS_MET_GetSetAttr(mdhandles[1], "RANGEBEGINNINGTIME", 
		  &RangeBeginTime);

  if( returnStatus != PGS_S_SUCCESS)
  {  PGS_SMF_GetMsg(&code,mnemonic,msg);
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  }

  returnStatus = sscanf( RangeBeginTime, "%2s:%2s", hour, minute);

  if (returnStatus != 2)
  {  sprintf(msg, "unable to read: RangeBeginTime\n\n");
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  } else {
     sprintf(msg, "RangeBeginTime is: %2s:%2s\n", hour, minute);
     OMI_SMF_setmsg(OMI_S_MET_READ,msg,FUNCTION_NAME, two);
  }

  /* Get VERSIONID from PCF */
  version = 1;
  returnStatus = PGS_PC_GetConfigData(VersionID_LUN, VID);

  VersionID = atoi(VID);

  if( returnStatus != PGS_S_SUCCESS)
  {  PGS_SMF_GetMsg(&code,mnemonic,msg);
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  } else {
     sprintf(msg, "VersionID is: %d\n", VersionID);
     OMI_SMF_setmsg(OMI_S_MET_READ,msg,FUNCTION_NAME, two);
  }

  /* Get ORBITNUMBER from Metadata */
  returnStatus = PGS_MET_GetSetAttr(mdhandles[1], "ORBITNUMBER.1", &OrbitNumber);

  if( returnStatus != PGS_S_SUCCESS)
  {  PGS_SMF_GetMsg(&code,mnemonic,msg);
     OMI_SMF_setmsg(OMI_E_MET_READ,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_READ;
  } else {
     sprintf(msg, "OrbitNumber is: %d\n", OrbitNumber[0]);
     OMI_SMF_setmsg(OMI_S_MET_READ,msg,FUNCTION_NAME, two);
  }

  /* Get current processing time and convert to GMT */
  time(&time_sec);
  time_ptr = gmtime(&time_sec);
  if( time_ptr == NULL )
  {  sprintf(msg, "%s", "gmtime convert to local UTC failed");
     OMI_SMF_setmsg(OMI_E_GENERAL, msg, FUNCTION_NAME, zero );
     sprintf( localGranuleID, "NOT_SET" );
     return OMI_E_GENERAL;
  } else {
     sprintf( proc_year, "%.4d", 1900+(time_ptr->tm_year) );
     sprintf( proc_month, "%.2d", time_ptr->tm_mon+1 );
     sprintf( proc_day, "%.2d", time_ptr->tm_mday );
     sprintf( proc_hour, "%.2d", time_ptr->tm_hour );
     sprintf( proc_minute, "%.2d", time_ptr->tm_min );
     sprintf( proc_second, "%.2d", time_ptr->tm_sec );
  }

  /* Construct localGranuleID string based on the ProcessingLevel */
  switch (ProcessingLevel)
  {
    case L1:
	    
      sprintf( localGranuleID, 
        "OMI-Aura_L1-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he4",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);

      break;

    case L2:
	    
      sprintf( localGranuleID, 
        "OMI-Aura_L2-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he5",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);

      if (strcmp(ShortName,"N7AERUV") == 0)
      { sprintf( localGranuleID, 
        "TOMS-N7_L2-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he5",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);
      }

      break;

    case L3:
	    
      sprintf( localGranuleID, 
        "OMI-Aura_L3-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he5",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);

      break;

    case L2A:
	    
      sprintf( localGranuleID, 
        "OMI-Aura_L2A-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he5",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);

      break;

    case L2B:
	    
      sprintf( localGranuleID, 
        "OMI-Aura_L2B-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he5",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);

      break;

    case L2G:
	    
      sprintf( localGranuleID, 
        "OMI-Aura_L2G-%s_%4sm%2s%2st%2s%2s-o%05d_v%03d-%4sm%2s%2st%2s%2s%2s.he5",
        ShortName, year, month, day, hour, minute, OrbitNumber[0], VersionID, 
	proc_year, proc_month, proc_day, proc_hour, proc_minute, proc_second);

      break;

    default:

      sprintf(msg, "%s", "Processing Level for Granule ID Invalid");
      OMI_SMF_setmsg(OMI_E_GRAN_ID, msg, FUNCTION_NAME, zero );
      return OMI_E_GRAN_ID;
   
  }      

  /* Set the LOCALGRANULEID metadata parameter */
  local = (void *) localGranuleID;
  returnStatus = PGS_MET_SetAttr( mdhandles[1], "LOCALGRANULEID", &local );

  if( returnStatus != PGS_S_SUCCESS)
  {  PGS_SMF_GetMsg(&code,mnemonic,msg);
     OMI_SMF_setmsg(OMI_E_MET_WRITE,msg,FUNCTION_NAME, zero);
     return OMI_E_MET_WRITE;
  } else {
     sprintf(msg, "Local Granule ID is: %s\n", localGranuleID);
     OMI_SMF_setmsg(OMI_S_GRAN_ID,msg,FUNCTION_NAME, two);
  }

  return OMI_S_SUCCESS; 
}
