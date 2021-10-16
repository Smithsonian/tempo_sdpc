/*------------------------------------------------------------------------------*\
** Config.h -- not for users
**
** If you wish to call any of the GetConfig_X_() functions, include GetConfig.h
** in your source (and not this file).
\*------------------------------------------------------------------------------*/

#include <cfortran.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "O3_PEATE_Common.h"
#include "Messages.h"
#include "GetConfig.h"

extern char* ControlFile;		/* control file (null-terminated str)	*/
extern yTB*  cfg;			/* CONFIG_READER main data structure	*/
extern int   nry;			/* number of elements in cfg[]		*/

/*------------------------------------------------------------------------------*/
