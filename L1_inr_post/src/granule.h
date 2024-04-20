#ifndef __TEMPO_GRANULE_INCLUDE__
#define __TEMPO_GRANULE_INCLUDE__ 1

/** @file granule.h
 *  @brief Load granule geolocation information,
 *         set geolocation-dependent variables
 */

#include <libconfig.h>

#include "elevation.h"
#include "land_cover.h"
#include "snow.h"

typedef struct Granule_Type Granule_Type;

struct Granule_Type
{
   void (*gt_close)(Granule_Type *);
   /**<  Free resources associated with Granule_Type object
    * @param[in]  et  pointer to Granule_Type object from granule_open
    */

   int (*gt_set_snow_ice_fraction)(Granule_Type *, Snow_Type *);
   /**< Set snow_ice_fraction variable
    * @param[in]  gt   pointer to Granule_Type object from granule_open
    * @param[in]  sn   pointer to Snow_Type object from snow_init
    * @return 0 on success, -1 on error
    */

   int (*gt_set_elevation)(Granule_Type *, const Elevation_Type *);
   /**<  Set variables associated with elevation
    * @param[in]  gt  pointer to Granule_Type object from granule_open
    * @param[in]  et  pointer to Elevation_Type object from elevation_init
    * @return 0 on success, -1 on error
    */

   int (*gt_set_object_angles)(Granule_Type *, int);
   /**<  Set solar and satellite viewing angles
    * @param[in]  gt  pointer to Granule_Type object from granule_open
    * @param[in]  is_radt  integer, non-zero if this is a RADT product
    * @return 0 on success, -1 on error
    */

   int (*gt_set_ground_pixel_flags)(Granule_Type *, double, double,
                                    const Land_Cover_Type *);
   /**<  Set ground pixel flags
    * @param[in]  gt  pointer to Granule_Type object from granule_open
    * @param[in]  max_glint_angle  Maximum angle between satellite viewing
    *                  vector and nominal solar reflection vector for which
    *                  the glint possibility flag should be set.
    * @param[in]  max_eclipse_angle  Maximum angle between sun and moon
    *                  for which the eclipse flag should be set.
    * @param[in]  lc  pointer to Land_Cover_Type object from land_cover_init
    * @return 0 on success, -1 on error
    */

   int (*gt_set_earth_sun_distance) (Granule_Type *);
   /**< Set Earth-Sun distance
    * @param[in] gt  pointer to Granule_Type object from granule_open
    * @return 0 on success, -1 on error
    *
    * This function uses the (X,Y,Z) coordinates of the Sun provided
    * by INR processing to define an Earth-Sun distance global variable
    * in the output NetCDF file.
    */

   int (*gt_is_twilight_granule) (Granule_Type *, int *);
   int (*gt_set_exposure_time_valid_max) (Granule_Type *gt, float);

#ifdef GRANULE_PRIVATE_DATA
   GRANULE_PRIVATE_DATA
#endif
};

/**  Initialize Granule_Type object
 * @param[in] file  path to geolocated radiance file
 * @param[in] correct_parallax  if non-zero, apply the parallax correction
 * @param[in] meta              pointer to an initialized instance of
 *                              \a TIO_Meta_Type
 * @param[in] cfg               pointer to an initialized instance of
 *                              \a config_t
 * @return pointer to initialized Granule_Type on success, NULL on error
 */
extern Granule_Type *granule_open (const char *file, int correct_parallax,
                                   TIO_Meta_Type *meta, config_t *cfg);

#endif
