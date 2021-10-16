/*------------------------------------------------------------------------------*\
** MetaConfig.c
**
** Writes a copy of the control file into the specified H5 file.
**
** See ReadControlFile.c for more information.
\*------------------------------------------------------------------------------*/

#include "Config.h"

/*------------------------------------------------------------------------------*/
/* FORTRAN binding								*/

FCALLSCFUN1 ( INT, MetaConfig, METACONFIG, metaconfig, PINT )

/*------------------------------------------------------------------------------*/
/*  fid = file ID identifies the H5 file to write into				*/

int MetaConfig(hid_t *fid)
{
   int    ierr = 0;
   herr_t herr;
   char   name[CFG_VAL_LEN];

   sprintf(name, "ControlFileContents");

   herr = H5LTmake_dataset_string(*fid, name, ControlFile);

   if ( herr < 0 )
   {
      Error("MetaConfig: Error writing <%s> input pointer data", name);
      ierr = (int)herr;
   }

   return(ierr);
}

/*------------------------------------------------------------------------------*/
