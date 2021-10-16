/****************************************************************************
BEGIN_FILE_PROLOG:

FILENAME:

        OMI_MET_Tools.h

DESCRIPTION:
        This file contains struct and function prototype information for 
	the OMI_MET_Tools module. 
AUTHOR:
        Ellyne Kinney / SSAI

HISTORY:
        24-APR-03 EKK Initial version

END_FILE_PROLOG:
 *****************************************************************************/
#ifndef OMI_MET_TOOLS_H        /* avoid re-inclusion */
#define OMI_MET_TOOLS_H

#include <PGS_MET.h>
#include <PGS_SMF.h>

enum err_level { zero, one, two, three};

enum processing_level { L1=1, L2, L3, L2A, L2B, L2G };

/* function prototypes */

PGSt_SMF_status           
OMI_localGranuleID( enum processing_level ProcessingLevel,
		    PGSt_MET_all_handles  mdhandles);

#endif
