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
   int (*sgt_geosat_longitude)(const Solar_Geom_Type *sgt, double *lon);
   int (*sgt_solar_zenith_angle)(const Solar_Geom_Type *sgt,
                                 double jd_utc, double lon, double lat,
                                 double *psza);
   int (*sgt_sat_sun_angle)(Solar_Geom_Type *sgt, double jd_utc,
                            double *psun_angle);
   int (*sgt_print_params)(Solar_Geom_Type *sgt, const char *, FILE *);

#ifdef SGT_PRIVATE_DATA
   SGT_PRIVATE_DATA
#endif
};

extern Solar_Geom_Type *solar_geom_init (config_t *cfg);

#endif
