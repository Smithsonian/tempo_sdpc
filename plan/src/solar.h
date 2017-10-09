#ifndef __PLAN_SGT_H__
#define __PLAN_SGT_H__ 1

/** @file solar.h
 *  @brief Interface for solar illumination geometry
 */

#include <libconfig.h>

typedef struct Solar_Geom_Type Solar_Geom_Type;

struct Solar_Geom_Type
{
   void (*sgt_delete)(Solar_Geom_Type *sgt);
   /**< delete this \c Solar_Geom_Type object */

   int (*sgt_geosat_longitude)(const Solar_Geom_Type *sgt, double *lon);
   /**< retrieve the geostationary satellite longitude */

   int (*sgt_solar_zenith_angle)(const Solar_Geom_Type *sgt,
                                 double jd_utc, double lon, double lat,
                                 double *psza);
   /**< compute the solar zenith angle for the specified place and time */

   int (*sgt_sat_sun_angle)(Solar_Geom_Type *sgt, double jd_utc,
                            double *psun_angle);
   /**< compute the angle between the sun and the instrument boresight */

   int (*sgt_print_params)(Solar_Geom_Type *sgt, const char *, FILE *);
   /**< print selected parameters to the specified open FILE pointer */

#ifdef SGT_PRIVATE_DATA
   SGT_PRIVATE_DATA
#endif
};

/** Initialize an object of \ref Solar_Geom_Type
 * @param[in] cfg  Pointer to the open configuration file
 * @return on success, a pointer to an initialized object of \c Solar_Geom_Type;
 *         on error, a NULL pointer
 */
extern Solar_Geom_Type *solar_geom_init (config_t *cfg);

#endif
