
#ifndef __COMPUTE_SCAN_ANGLES_H__
#define __COMPUTE_SCAN_ANGLES_H__

#include <tempo_geo.h>

typedef struct
{
   double azimuth;
   double elevation;
}
AziElev_Type;

typedef struct
{
   double ewbias_rad;      /**< East-West bias [rad] */
   double nsbias_rad;      /**< North-South bias [rad] */
   double clockbias_rad;   /**< clock angle [rad] */
   double theta0_rad;      /**< mirror incidence angle [rad] */
}
Geometry_Param_Type;

extern int __compute_scan_angles (const EarthPoint *c_pt, double sat_lon, AziElev_Type *apt);
extern int __set_geometry_params (double ewbias, double nsbias, double clockingbias, double telescopeOffset);

#endif
