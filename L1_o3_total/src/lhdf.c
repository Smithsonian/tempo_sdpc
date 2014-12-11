#include <stdio.h>
#include "lhdf.h"

int
Lhdf_selectSDS( int32  sd_fid,
                sds_t *sds_p )
/******************************************************************************
!C
!Description: Lhdf_selectSDS-- select an SDS (Science Data Set no
 more than rank 3) i.e., obtaing the sds_p->id, from the HDF or 
 EOS_HDF file defined by sd_fid.  this function compares input 
 rank, data type with those of the selected SDS to make sure 
 they are exactly the same.

!Input Parameters:
 sd_fid         Science data file id
 sds_p          Data structure defining the SDS

!Output Parameters:
 sds_p          Data structure defining the SDS with
                member 'id' updated.
                dimensions filled.

 (returns)      status:
                  LH2G_SUCCESS - normal return
                  LH2G_ERROR - error return

!Revision History:
   Revision 00.0  02/08/1998
   Kai Yang,  UMBC
   Original version.

   Revision 00.0  03/14/2000
   Kai Yang,  UMBC
   modified Lhdf_Get so that it can read 1-d sds.  Before it was a 2-d sds
   reader only.

!Team-unique Header:
!References and Credits:

!Design Notes:

!END
***************************************************************************/
{
  int32 sds_index;
  int   irank;
  int32 dims[MAXRANK];
  char  name[H4_MAX_NC_NAME];
  int32 rank,
        num_type,
        attributes;

  if( (sds_index = SDnametoindex( sd_fid, sds_p->name )) == FAIL )
  {  fprintf( stderr, "SDnametoindex failed %s", 
             "Lhdf_selectSDS, lhdf.c" );  
     return FAIL;
  }

  if( (sds_p->id = SDselect( sd_fid, sds_index )) == FAIL )
  {  fprintf( stderr, "index-to-sds_id (SDselect) failed, %s", 
             "Lhdf_selectSDS, lhdf.c" );  
     return FAIL;
  }

   /* get the rank and dim_sizes */
  if(  SDgetinfo( sds_p->id,
                  name,
                  &rank,
                  dims,
                  &num_type,
                  &attributes ) == FAIL )
  {  fprintf( stderr, "SDgetinfo failed %s", 
             "Lhdf_selectSDS, lhdf.c" );  
     return FAIL;
  }

  if( num_type != sds_p->type )
  {  fprintf( stderr, "data type not matched type_out=%d, type_in=%d, %s\n", 
             (int) num_type, (int) sds_p->type, "Lhdf_selectSDS, lhdf.c" );  
     return FAIL;
  }

  if( rank != sds_p->rank )
  {  fprintf( stderr, "rank not matched %s\n", 
             "Lhdf_selectSDS, lhdf.c" );  
     return FAIL;
  }

  for( irank = 0; irank < sds_p->rank; irank++ )
     sds_p->dim[irank].nval = dims[irank];

  return SUCCEED;
}

int
Lhdf_createSDS( int32  sd_fid,
                sds_t *sds )
{
  int irank;
  int32 dims[MAXRANK];

  for (irank = 0; irank < sds->rank; irank++)
    dims[irank] = sds->dim[irank].nval;

  if( (sds->id = SDcreate(sd_fid, sds->name, sds->type,
                          sds->rank, dims)) == FAIL )
  {  fprintf( stderr, "Lhdf_createSDS %s", sds->name );
     return FAIL;
  }

  for( irank = 0; irank < sds->rank; irank++ )
  {
    if( (sds->dim[irank].id = SDgetdimid(sds->id, irank)) == FAIL )
    {  fprintf( stderr, "Lh2g_createSDS %s", "dim_id" );
       return FAIL;
    }

/*  if( SDsetdimval_comp(sds->dim[irank].id, SD_DIMVAL_BW_INCOMP) == FAIL )
 *  {  fprintf( stderr, "Lh2g_createSDS %s", "bad compression" );
 *     return FAIL;
 *  }
 */
    if( SDsetdimname(sds->dim[irank].id, sds->dim[irank].name) == FAIL )
    {  fprintf( stderr, "Lh2g_createSDS %s", "set dimension failure" );
       return FAIL;
    }
  }
  return SUCCEED;
}


int 
Lhdf_Get( sds_t       *sds_p,
          void        *ptr )
/*
!C******************************************************************************
!Description: Lhdf_Get (get parameter) read the data for one set of 
 observations in one or more lines of an parameter stored in an HDF file.

!Input Parameters:
 file           pointer to Lh2g file definition structure
*sds_p          pointer of the sds structures to be read from the 2g file
 iline          start line to read (0 based)
 nline_read     number of lines to read

!Output Parameters:
 ptr            pointer to memory space ( allocated outside this function, 
                usually in the calling function) where retrived data are
                stored. 
 (returns)      status:
                  LH2G_SUCCESS - normal return
                  LH2G_ERROR - error return

!Revision History:
!Team-unique Header:
!References and Credits:

!Design Notes:

!END
******************************************************************************/
{
  int     irank;
  int32   start[MAXRANK], nval[MAXRANK];

  if( sds_p->rank > MAXRANK )
  {  fprintf( stderr, "input rank is greater than MAXRANK %s", 
             "Lhdf_Get, lhdf.c" );  
     return FAIL;
  }

  for( irank = 0; irank < sds_p->rank; irank++ )
  {  start[irank] = 0;
     nval[irank]  = sds_p->dim[irank].nval;
  }

  if( SDreaddata(sds_p->id, start, NULL, nval, ptr ) == FAIL )
  {  fprintf( stderr, "SDreaddata failed %s", 
             "Lhdf_selectSDS, lhdf.c" );  
     return FAIL;
  }

  return SUCCEED;
}

int 
Lhdf_Put( sds_t       *sds_p,
          void        *ptr )
{
  int     irank;
  int32   start[MAXRANK], nval[MAXRANK];


  for( irank = 0; irank < sds_p->rank; irank++ )
  {  start[irank] = 0;
     nval[irank]  = sds_p->dim[irank].nval;
  }

  if( SDwritedata(sds_p->id, start, NULL, nval, ptr ) == FAIL )
  {  fprintf( stderr, "SDreaddata failed %s",
             "Lhdf_selectSDS, lhdf.c" );
     return FAIL;
  }

  return SUCCEED;
}

