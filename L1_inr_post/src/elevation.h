#ifndef __TEMPO_ELEVATION_H__
#define __TEMPO_ELEVATION_H__ 1

/** @file elevation.h
 *  @brief Load elevation data, perform coordinate-based nearest-neighbor
 *         elevation lookup, compute overlap area-weighted mean elevation
 *         within pixel boundaries.
 */

#include <libconfig.h>

typedef struct Elevation_Type Elevation_Type;

struct Elevation_Type
{
   void (*et_delete)(Elevation_Type *);
   /**<  Free resources associated with Elevation_Type object
    * @param[in]  et  pointer to Elevation_Type object from \c elevation_init
    */

   int (*et_lookup)(const Elevation_Type *, unsigned int,
                    const double *, const double *, short *, short *);
   /**<  Nearest-neighbor elevation lookup
    * @param[in]  et   pointer to Elevation_Type object from \c elevation_init
    * @param[in]  num  number of lon,lat points to process
    * @param[in]  lon  pointer to longitude coordinates
    * @param[in]  lat  pointer to latitude coordinates
    * @param[out] elev  pointer to nearest-neighbor elevation mean values [meters]
    * @param[out] elev_stddev  pointer to nearest-neighbor elevation standard deviation values [meters]
    * @return 0 on success, -1 on error
    */

   int (*et_regrid)(const Elevation_Type *, unsigned int,
                    const double *, const double *, short *, short *);
   /**<  Compute overlap area-weighted mean elevation within pixel boundaries
    * @param[in]  et   pointer to Elevation_Type object from \c elevation_init
    * @param[in]  num_pixels  number of lon,lat pixels to process
    * @param[in]  lon_cnrp  pointer to pixel corner longitude coordinates, packed
    *                  in groups of 4, in CCW order around the pixel boundary.
    * @param[in]  lat_cnrp  pointer to pixel corner latitude coordinates, packed
    *                  in groups of 4, in CCW order around the pixel boundary.
    * @param[out] elev  pointer to overlap area-weighted elevation mean values [meters]
    * @param[out] elev_stddev  pointer to overlap area-weighted elevation standard deviation values [meters]
    * @return 0 on success, -1 on error
    */

#ifdef ELEVATION_TYPE_PRIVATE_DATA
   ELEVATION_TYPE_PRIVATE_DATA
#endif
};

/**  Initialize Elevation_Type object
 * @param[in] cfg   pointer to \c config_t object associated with
 *                  \c L1_inr_post parameter file
 * @return pointer to initialized Elevation_Type on success, NULL on error
 */
extern Elevation_Type *elevation_init (config_t *cfg);

#endif
