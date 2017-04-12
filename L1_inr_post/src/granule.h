#ifndef __TEMPO_GRANULE_INCLUDE__
#define __TEMPO_GRANULE_INCLUDE__ 1

/** @file granule.h
 *  @brief Load granule geolocation information,
 *         set geolocation-dependent variables
 */

#include "elevation.h"
#include "snow.h"
#include "land_cover.h"

typedef struct Granule_Type Granule_Type;

struct Granule_Type
{
   void (*gt_close)(Granule_Type *);
   /**<  Free resources associated with Granule_Type object
    * @param[in]  et  pointer to Granule_Type object from granule_open
    */

   int (*gt_set_elevation)(Granule_Type *, const Elevation_Type *);
   /**<  Set variables associated with elevation
    * @param[in]  gt  pointer to Granule_Type object from granule_open
    * @param[in]  et  pointer to Elevation_Type object from elevation_init
    * @return 0 on success, -1 on error
    */

   int (*gt_set_object_angles)(Granule_Type *);
   /**<  Set solar and satellite viewing angles
    * @param[in]  gt  pointer to Granule_Type object from granule_open
    * @return 0 on success, -1 on error
    */

   int (*gt_set_ground_pixel_flags)(Granule_Type *, double, double,
                                    const Snow_Type *,
                                    const Land_Cover_Type *);
   /**<  Set ground pixel flags
    * @param[in]  gt  pointer to Granule_Type object from granule_open
    * @param[in]  max_glint_angle  Maximum angle between satellite viewing
    *                  vector and nominal solar reflection vector for which
    *                  the glint possibility flag should be set.
    * @param[in]  max_eclipse_angle  Maximum angle between sun and moon
    *                  for which the eclipse flag should be set.
    * @param[in]  sn  pointer to Snow_Type object from snow_init
    * @param[in]  lc  pointer to Land_Cover_Type object from land_cover_init
    * @return 0 on success, -1 on error
    */

#ifdef GRANULE_PRIVATE_DATA
   GRANULE_PRIVATE_DATA
#endif
};

/**  Initialize Granule_Type object
 * @param[in] file  path to geolocated radiance file
 * @return pointer to initialized Granule_Type on success, NULL on error
 */
extern Granule_Type *granule_open (const char *file);

#endif
