/*
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
*/ 

#include <PGS_PC.h>
#include "PGS_OZT_52050.h"
#include "pcf.h"
#include "lhdf.h"
#include "omi_smf.h"
#include <cfortran.h>

sds_t
sds_terrainPressures =
{ -1, -1, 2,
  { {-1,  -1, ""}, {-1,  -1, ""}, {-1,  -1, ""},
    {-1,  -1, ""}, {-1,  -1, ""}, {-1,  -1, ""} },
  DFNT_FLOAT32, -1, "TerrainPressure", "TerrainPressure", "", 
  {0, 1014}, 9999
};

#define XDim 1080
#define YDim 540

static float32 terP[YDim][XDim];


#define FUNCTION_NAME "OMI_pixGetTerPres"
#define zero 0

PGSt_SMF_status 
OMI_pixGetTerPres( float  latitude,
                   float  longitude,
                   float *pterrain ) 
{  
   static int firsttime = -1;

   static float dlat, dlon;
   float  fx, fy;
   float  foo;
   int    lg_y, lg_x;
   /* float  frac; */
  
   if( firsttime == -1 )
   {   
      PGSt_integer    file_v;
      PGSt_SMF_status returnStatus;
      char            terrainpresfn[PGSd_PC_VALUE_LENGTH_MAX];
      char            msg[PGS_SMF_MAX_MSG_SIZE];
      int32           sd_fid;
      int32           /* im, */ status = -1;

      file_v = 1;
      /* get the Land Cover names */
      returnStatus = PGS_PC_GetReference( TERRAINPRES_LUN, &file_v, 
                                          terrainpresfn);
      if( returnStatus != PGS_S_SUCCESS )
      {  OMI_SMF_setmsg( OZT_E_INPUT,"get terrain pressure fn failed",
                         FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }

      sd_fid = SDstart( terrainpresfn, DFACC_RDONLY );
      if( sd_fid == FAIL )
      {  sprintf( msg, "Terrain pressure file not exist: %s", terrainpresfn );
         OMI_SMF_setmsg(OZT_E_INPUT,msg,FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }

      if( Lhdf_selectSDS( sd_fid, &sds_terrainPressures ) != SUCCEED )
      {  sprintf( msg, "selecting %s from %s failed", 
                      sds_terrainPressures.name, terrainpresfn );
         OMI_SMF_setmsg(OZT_E_HDF,msg,FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }

      if( YDim != sds_terrainPressures.dim[0].nval ||
          XDim != sds_terrainPressures.dim[1].nval )
      {  OMI_SMF_setmsg( OZT_E_INPUT,"Dimensions not matched",
                         FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }

      status =  Lhdf_Get( &sds_terrainPressures, terP );

      if( status != SUCCEED )
      {  sprintf( msg, "reading terrain pressure from file %s failed", 
                  terrainpresfn );
         OMI_SMF_setmsg(OZT_E_HDF,msg,FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }
      SDendaccess( sds_terrainPressures.id );

      SDend( sd_fid );

      dlat = 180.0/YDim;
      dlon = 360.0/XDim;
      firsttime = 0;
   }  

   fy   = ( 90.0 - latitude - dlat/2 )/dlat ;
   fx   = ( longitude + 180.0 )/dlon;

   lg_y = (int) fy;
   lg_x = (int) fx;
   fy   = fy - lg_y;
   fx   = fx - lg_x;


   if( lg_y < 0 )
      lg_y = 0;
   else if( lg_y >= YDim )
      lg_y = YDim-1;

   if( lg_x < 0 )
      lg_x = 0;
   else if( lg_x >= XDim )
      lg_x = XDim-1;

   if( lg_y == 0 || lg_y == YDim-1 || lg_x == 0 || lg_x == XDim-1)
   {  foo = terP[lg_y ][ lg_x ];
   }
   else 
   {  foo  = (1.0-fy)*(1.0-fx)*terP[lg_y  ][lg_x  ] \
           +      fy *(1.0-fx)*terP[lg_y+1][lg_x  ] \
           +      fy *     fx *terP[lg_y+1][lg_x+1] \
           + (1.0-fy)*     fx *terP[lg_y  ][lg_x+1];

   }
  
   *pterrain = foo;
   return OZT_S_SUCCESS;
}

/* FORTRAN bindings */

FCALLSCFUN3( INT, OMI_pixGetTerPres, OMI_PIXgETTERPRES, omi_pixgetterpres,\
             FLOAT, FLOAT, PFLOAT ) 
