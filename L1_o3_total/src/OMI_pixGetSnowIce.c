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
sds_snowIcePercent =
{ -1, -1, 2,
  { {-1,  -1, ""}, {-1,  -1, ""}, {-1,  -1, ""},
    {-1,  -1, ""}, {-1,  -1, ""}, {-1,  -1, ""} },
  DFNT_FLOAT32, -1, "snow_ice00", "snow_ice00", "", 
  {0, 100}, 9999
};

#define XDim 360
#define YDim 180

static float32 snowIceP1[YDim][XDim];
static float32 snowIceP2[YDim][XDim];

float
OMI_mmddInterp( int  year,
                int  month,
                int  day,
                int *dday,
                int *mm_cur,
                int *mm_pre );

#define FUNCTION_NAME "OMI_pixGetSnowIce"
#define zero 0

PGSt_SMF_status 
OMI_pixGetSnowIce( float  latitude,
                   float  longitude,
                   int    year, 
                   int    month,
                   int    day,
                   float *snowicePercent ) 
{  
   static int mm_current = -1;
   static int yy = -1;

   static float dlat, dlon;
   float  fx, fy;
   float  foo1, foo2;
   int    lg_y, lg_x;
   int    dday, mm_cur, mm_pre; 
   float  frac;
  
   frac = OMI_mmddInterp( year, month, day, &dday, &mm_cur, &mm_pre ); 

   if( mm_current != mm_cur || yy != year )
   {   
      PGSt_integer    file_v;
      PGSt_SMF_status returnStatus;
      char            snowicefn[PGSd_PC_VALUE_LENGTH_MAX];
      char            msg[PGS_SMF_MAX_MSG_SIZE];
      int32           sd_fid;
      int32           im, status = -1;

      mm_current = mm_cur;
      yy         = year;

      file_v = 1;
      /* get the snow/ice Cover names */
      returnStatus = PGS_PC_GetReference( CLIMSNOWICE_LUN, &file_v, snowicefn);
      if( returnStatus != PGS_S_SUCCESS )
      {  OMI_SMF_setmsg( OZT_E_INPUT,"get clim snow/ice fn failed",
                         FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }

      sd_fid = SDstart( snowicefn, DFACC_RDONLY );
      if( sd_fid == FAIL )
      {  sprintf( msg, "Snow/ice file not exist: %s", snowicefn );
         OMI_SMF_setmsg(OZT_E_INPUT,msg,FUNCTION_NAME, zero);
         return OZT_E_FAILURE;
      }

      for( im = 0; im < 2; im++ )
      {   if( im == 0 ) sprintf( sds_snowIcePercent.name, "%s%02d", 
                                 "snow_ice", mm_cur );
          if( im == 1 ) sprintf( sds_snowIcePercent.name, "%s%02d", 
                                 "snow_ice", mm_pre );

          if( Lhdf_selectSDS( sd_fid, &sds_snowIcePercent ) != SUCCEED )
          {  sprintf( msg, "selecting %s from %s failed", 
                      sds_snowIcePercent.name, snowicefn );
             OMI_SMF_setmsg(OZT_E_HDF,msg,FUNCTION_NAME, zero);
             return OZT_E_FAILURE;
          }

          if( YDim != sds_snowIcePercent.dim[0].nval ||
              XDim != sds_snowIcePercent.dim[1].nval )
          {  OMI_SMF_setmsg( OZT_E_INPUT,"Dimensions not matched",
                             FUNCTION_NAME, zero);
             return OZT_E_FAILURE;
          }

          if( im == 0 ) status =  Lhdf_Get( &sds_snowIcePercent, snowIceP1 );
          if( im == 1 ) status =  Lhdf_Get( &sds_snowIcePercent, snowIceP2 );

          if( status != SUCCEED )
          {  sprintf( msg, "reading snow/ice percent from file %s failed", 
                      snowicefn );
             OMI_SMF_setmsg(OZT_E_HDF,msg,FUNCTION_NAME, zero);
             return OZT_E_FAILURE;
          }
          SDendaccess( sds_snowIcePercent.id );
      } /* end for( im = 0; im < 2; im++ ) */

      SDend( sd_fid );

      dlat = 180.0/YDim;
      dlon = 360.0/XDim;
   }  

   fy   = ( 90.0 - latitude  )/dlat ;
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

   foo1 = snowIceP1[lg_y ][ lg_x ];
   foo2 = snowIceP2[lg_y ][ lg_x ];

/* commented out interpolation in (lat,lon), interploation is
   only done in day of month to match TOMS V8.

   if( lg_y == 0 || lg_y == YDim-1 || lg_x == 0 || lg_x == XDim-1)
   {  foo1 = snowIceP1[lg_y ][ lg_x ];
      foo2 = snowIceP2[lg_y ][ lg_x ];
   }
   else 
   {  foo1 = (1.0-fy)*(1.0-fx)*snowIceP1[lg_y  ][lg_x  ] \
           +      fy *(1.0-fx)*snowIceP1[lg_y+1][lg_x  ] \
           +      fy *     fx *snowIceP1[lg_y+1][lg_x+1] \
           + (1.0-fy)*     fx *snowIceP1[lg_y  ][lg_x+1];

      foo2 = (1.0-fy)*(1.0-fx)*snowIceP2[lg_y  ][lg_x  ] \
           +      fy *(1.0-fx)*snowIceP2[lg_y+1][lg_x  ] \
           +      fy *     fx *snowIceP2[lg_y+1][lg_x+1] \
           + (1.0-fy)*     fx *snowIceP2[lg_y  ][lg_x+1];
   }
 */
   *snowicePercent = foo1+(foo2-foo1)*frac;
   return OZT_S_SUCCESS;
}

/* FORTRAN bindings */

FCALLSCFUN6( INT, OMI_pixGetSnowIce, OMI_PIXGETSNOWICE, omi_pixgetsnowice,\
             FLOAT, FLOAT, INT, INT, INT, PFLOAT ) 
