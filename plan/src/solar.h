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
   int (*sgt_boresight_angles)(const Solar_Geom_Type *sgt, double *elev_deg, double *azi_deg);
   /**< retrieve the boresight elevation and azimuth in degrees */

   int (*sgt_solar_zenith_angle)(const Solar_Geom_Type *sgt,
                                 double jd_utc, double lon, double lat,
                                 double *psza);
   /**< compute the solar zenith angle for the specified place and time */

   int (*sgt_solar_xyz) (Solar_Geom_Type *sgt, double jd_utc, double sun_itrs[3]);
   /**< compute the ITRS coordinates of the sun (e.g. WGS84 XYZ) */

   int (*sgt_sat_sun_position)(Solar_Geom_Type *sgt, double jd_utc,
                               double *ptheta, double *pphi, double *earth_sun_distance);
   /**< Compute the position angles of the sun relative to the instrument boresight.
    * @param[in] sgt     Pointer to a struct of type \a Solar_Geom_Type
    * @param[in] jd_utc  Julian date for the UTC time of interest.
    * @param[out] ptheta  Polar angle of the sun [deg] (must be non-NULL)
    * @param[out] pphi    Azimuth angle of the sun [deg], CCW from the northern end of
    *                     the instrument slit.  (NULL is ok)
    * @param[out] earth_sun_distance   Distance of the sun from geocenter [km] (NULL is ok)
    */

   int (*sgt_print_params)(const Solar_Geom_Type *sgt, const char *, FILE *);
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
