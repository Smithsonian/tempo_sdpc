#ifndef __PRIVATE_PROCESS_INCLUDE_H__
#define __PRIVATE_PROCESS_INCLUDE_H__ 1

#include "granule.h"
#include "pixelqf.h"

typedef struct
{
   Granule_Exprec_Type *exprec;
   Image_Type *img_err;
   Dark_Trend_Type dark_trend;
   float storage_region_dark[4];
   float fpa_temp;
   float fpe_temp;
   double earth_sun_distance;
   int index;
}
Exprec_Meta_Type;

#endif
