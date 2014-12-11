/******************************************************************************
!C-INC
!Description:
   include file lhdf library.

!Input Parameters:

!Output Parameters:

!Revision History:
 Revision 00.1  02/08/1998
 Kai Yang,  UMBC

!Team-unique Header:
  This software was developed by:
    MODIS Land Science Team for the National Aeronautics and
    Space Administration, Goddard Space Flight Center, under
    NASA <contract or task number>

!References and Credits:
    Developers:
      Kai Yang
      University of Maryland Baltimore County
!Design Notes:

!END
****************************************************************************/

#ifndef LHDF_H
#define LHDF_H
#include "mfhdf.h"

#define MAXRANK 6
typedef struct
{ int32 index,
        id,
        rank;
  struct
  { int32 nval,
          id;
    char  name[80];
  } dim[MAXRANK];

  int32     type,
            nattr;
  char      name[80];
  char      long_name[100];  /* long name */
  char      units[80];      /* units */
  long int  vrange[2];      /* valid range */
  long int  fillv;          /* fill value */
} sds_t;

int
Lhdf_selectSDS( int32  sd_fid,
                sds_t *sds_p );

int
Lhdf_createSDS( int32  sd_fid,
                sds_t *sds );
int
Lhdf_Get( sds_t       *sds_p,
          void        *ptr );

int
Lhdf_Put( sds_t       *sds_p,
          void        *ptr );
#endif
