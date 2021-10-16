/*******************************************************************************
BEGIN_FILE_PROLOG:

FILENAME:
        OMI_PC_GetReference

DESCRIPTION:
        Calls the SDP Toolkit function PGS_PC_GetReference and splits the 
	file reference into a path and basename.

AUTHOR:
        Ellyne Kinney / SSAI

HISTORY:
        17-JUL-03      EKK     Initial version

END_FILE_PROLOG
*******************************************************************************/

/* include files */

#include <stdio.h>
#include <PGS_PC.h>
#include "omi_smf.h"          /* include PGS_PC.h PGS_SMF.h */
#include "PGS_OMI_52090.h"

/***************************************************************************
BEGIN_PROLOG:

TITLE:
       Gets path and basename from the Process Control file

NAME:
       OMI_PC_GetReference()

SYNOPSIS:
C:
       #include <PGS_PC.h>
       #include "omi_smf.h"

       PGSt_SMF_status
       OMI_PC_GetReference ( PGSt_PC_Logical  prodID,
                             PGSt_integer    *version, 
			     char            *path,
                             char            *basename)

FORTRAN:
       include "PGS_PC_9.f"
       include "PGS_PC.f"
       include "PGS_SMF.h"

       integer function omi_pc_getreference(prodID, version, path, basename)

       integer          prodID
       integer          version
       character*1025   path
       character*1025   basename

DESCRIPTION:
        This tool may be used to obtain a physical reference (file name or
        universal identifier) from a logical identifier. The file basename
	and path name are returned as separate variables.

INPUTS:
        Name            Description                     Units   Min     Max

        prodID          User defined constant identifier
                        that internally represents the
                        current product.

        version         Version of reference to get.
                        Remember,  for Product Input
                        files and Product Output files
                        there can be a many-to-one
                        relationship.  See OUTPUTS.

OUTPUTS:
        Name            Description                     Units   Min     Max

        path            The actual file path     
                        returned as a string.

        basename        The actual filename 
                        returned as a string.

        version         The number of versions remaining
                        for the requested Product ID.

RETURNS:
        OMI_S_SUCCESS               successful execution
        PGSPC_W_NO_REFERENCE_FOUND  link number does not have the data that
                                    mode is requesting
        PGSPC_E_DATA_ACCESS_ERROR   problem while accessing PCS data
	OMI_E_GENERAL               could not split reference 
	                            into path and basename

EXAMPLES:

C:
        #define         SML1BGEO 2530

        PGSt_integer    version;
        char            path[PGSd_PC_FILE_PATH_MAX];
        char            basename[PGSd_PC_FILE_PATH_MAX];
        PGSt_SMF_status returnStatus;

        /# Get first version of the file #/
        version = 1;

        returnStatus = OMI_PC_GetReference(SML1BGEO,&version,path,basename);

        /# version now contains the number of versions remaining #/

        if (returnStatus != OMI_S_SUCCESS)
                goto EXCEPTION;
        else
        { /# perform necessary operations on file #/ }
                        .
                        .
                        .
        EXCEPTION:
                return returnStatus;

FORTRAN:
                IMPLICIT NONE

                integer         version
                character*1025  path
                character*1025  basename
                integer         returnstatus
                integer         omi_pc_getreference
                integer         sml1bgeo
                parameter       (modis1a = 2530)

C               Get the first version of the file
                version = 1

                returnstatus = getreference(sml1bgeo,version,path,basename)

                if (returnstatus .ne. omi_s_success)
                        goto 9999
                else
C                       perform necessary operations on file
                        .
                        .
                        .
        9999    return

NOTES:
        All reference identifier strings are guaranteed to be no greater
        than PGSd_PC_FILE_PATH_MAX characters in length (see PGS_PC.h).

        The version returns the number of files remaining for the product
        group.  For example, if there are eight (8) versions of a file
        when the user requests version one (1) the value seven (7) is
        returned in version.  When the user requests version two (2) the
        value six (6) is returned in version, etc.  Therefore, it is not
        recommended to use version as a loop counter that is also passed
        into PGS_PC_GetReference().

        The one-to-many relationship is only supported with Product Input
        and Product Output files.  Files listed in other sections of the
        Process Control File are only supported on a one-to-one relationship.
        Therefore, the variable version in not used when searching for files
        in sections of the PCF other than the Product Input and Product Output
        files.  If the user is knowingly searching for a file that is not in
        the Product Input or Product Output files of the PCF the user should
        create a version variable, set it equal to one, and pass in the address
        of that variable.


REQUIREMENTS:
       PGSTK-1290

DETAILS:
        NONE

GLOBALS:
        NONE

FILES:
        NONE

FUNCTIONS_CALLED:
        OMI_PC_GetReference              Get process control information.
        PGS_SMF_GetMsg                   Get error/status message.
	OMI_SMF_setmsg                   Set error/status message.
	OMI_PC_RefIDsplit                Split Reference into path and 
	                                 basename

END_PROLOG:
***************************************************************************/

#define FUNCTION_NAME "OMI_PC_GetReference()"

PGSt_SMF_status
OMI_PC_GetReference(                    /* get physical reference */
    PGSt_PC_Logical   prodID,           /* logical value */
    PGSt_integer      *version,         /* version of reference requested */
    char              *path,            /* actual path string */
    char              *basename)        /* actual basename string */

{
    enum err_level { zero, one, two, three};

    char referenceID[PGSd_PC_FILE_PATH_MAX] = "\0"; 

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

    /* Get the Reference ID using the standard SDP Toolkit routine */
    returnStatus = PGS_PC_GetReference(prodID, version, referenceID);

    /* Test for error condition and return error */
    if( returnStatus != PGS_S_SUCCESS )
    {  PGS_SMF_GetMsg(&code,mnemonic,msg);
       OMI_SMF_setmsg(OMI_E_PCF_REFID,msg,FUNCTION_NAME, zero);
       return OMI_E_PCF_REFID;
    }

    /* Split the reference ID into path and basename strings */
    returnStatus = OMI_PC_RefIDsplit(referenceID, path, basename);

    /* Test for error condition and return error */
    if( returnStatus != OMI_S_SUCCESS )
    { sprintf(msg, "Could not split reference into path and basename.");
      OMI_SMF_setmsg(OMI_E_PCF_SPLIT,msg,FUNCTION_NAME, zero);
      return OMI_E_PCF_SPLIT;
    }

    /* return success */
    return OMI_S_SUCCESS;
}

/***************************************************************************
BEGIN_PROLOG:

TITLE:
       Gets path and basename from the reference ID

NAME:
       OMI_PC_RefIDsplit()

SYNOPSIS:
C:
       #include "omi_smf.h"

       PGSt_SMF_status
       OMI_PC_RefIDsplit( char    *referenceID, 
			  char    *path,
                          char    *basename)

FORTRAN:
       include "PGS_SMF.h"

       integer function omi_pc_refidsplit(referenceID, path, basename)

       character*1025   referenceID
       character*1025   path
       character*1025   basename

DESCRIPTION:
        This tool may be used to split the reference ID into separate 
	strings for path name and basename.

INPUTS:
        Name            Description                     Units   Min     Max

        referenceID     

OUTPUTS:
        Name            Description                     Units   Min     Max

        path            The actual file path     
                        returned as a string.

        basename        The actual filename 
                        returned as a string.


RETURNS:
        OMI_S_SUCCESS               successful execution
	OMI_E_GENERAL               ReferenceID is an empty string

EXAMPLES:

C:

        char            referenceid[PGSd_PC_FILE_PATH_MAX];
        char            path[PGSd_PC_FILE_PATH_MAX];
        char            basename[PGSd_PC_FILE_PATH_MAX];
        PGSt_SMF_status returnStatus;

        /# Split the referenceID into path and basename strings #/ 
        returnStatus = OMI_PC_RefIDsplit(referenceid,path,basename);

        if (returnStatus != OMI_S_SUCCESS)
                goto EXCEPTION;
        else
        { /# perform necessary operations on file #/ }
                        .
                        .
                        .
        EXCEPTION:
                return returnStatus;

FORTRAN:
                IMPLICIT NONE

                character*1025  referenceid
                character*1025  path
                character*1025  basename
                integer         returnstatus
                integer         omi_pc_refidsplit

C               Split the referenceID into path and basename strings

                returnstatus = refidsplit(referenceid,path,basename)

                if (returnstatus .ne. omi_s_success)
                        goto 9999
                else
C                       perform necessary operations on file
                        .
                        .
                        .
        9999    return

NOTES:
        NONE

DETAILS:
        NONE

GLOBALS:
        NONE

FILES:
        NONE

FUNCTIONS_CALLED:
	OMI_SMF_setmsg                   Set error/status message.


END_PROLOG:
***************************************************************************/

#define FUNCTION_NAME2 "OMI_PC_RefIDsplit()"

PGSt_SMF_status
OMI_PC_RefIDsplit(char referenceID[],   /* actual reference ID string */
		  char path[],          /* actual path string */
		  char basename[])      /* actual basename string */
{
    enum err_level { zero, one, two, three};
    int   i,j;
    int   size = 0;

    char                 mnemonic[PGS_SMF_MAX_MNEMONIC_SIZE];
                                                    /* status mnemonic
                                                       returned by
                                                       PGS_SMF_GetMsg() */
    char                 msg[PGS_SMF_MAX_MSG_SIZE]; /* status messsage
                                                       returned by call to
                                                       PGS_SMF_GetMsg() */
    PGSt_SMF_status      code;                      /* status code returned
                                                       by PGS_SMF_GetMsg() */


    /* Get reference ID string length */
    size = strlen(referenceID);

    /* Test for zero length string, return error */
    if(size<=0)
    {  sprintf( msg,"ReferenceID is an empty string!");
       OMI_SMF_setmsg(OMI_E_GENERAL,msg,FUNCTION_NAME2, zero);
       return OMI_E_GENERAL;
    }

    /* Find the last path deliminator */
    for(i=size; i>=0; --i)
	if(referenceID[i] == '/') 
		break;

    /* Copy the contents left of the delimintor into the path string */
    for(j=0; j<i; j++)
	path[j] = referenceID[j];
    path[i] = '\0';

    /* Copy the contents right of the delimintor into the basename string */
    for(j=i; j<size; j++)
	basename[j-i] = referenceID[j+1];
    basename[size-i] = '\0';

    /* return success */
    return OMI_S_SUCCESS;
}
