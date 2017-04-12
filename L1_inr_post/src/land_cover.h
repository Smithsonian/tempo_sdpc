#ifndef __TEMPO_LAND_COVER_H__
#define __TEMPO_LAND_COVER_H__ 1

/** @file land_cover.h
 *  @brief Load MODIS land cover data, perform coordinate-based
 *         nearest-neighbor lookups
 */

#include <libconfig.h>

typedef struct Land_Cover_Type Land_Cover_Type;

struct Land_Cover_Type
{
   void (*lc_delete)(Land_Cover_Type *);
   /**<  Free resources associated with Land_Cover_Type object
    * @param[in]  lc  pointer to Land_Cover_Type object from land_cover_init
    */

   int (*lc_lookup_type1)(const Land_Cover_Type *, unsigned int,
                          const double *, const double *,
                          unsigned char *);
   /**<  Nearest-neighbor lookup of MODIS Type 1 land cover mask
    * @param[in]  lc   pointer to Land_Cover_Type object from land_cover_init
    * @param[in]  num  number of lon,lat points to process
    * @param[in]  lon  pointer to longitude coordinates
    * @param[in]  lat  pointer to latitude coordinates
    * @param[out] mask  pointer to nearest-neighbor land cover mask value
    * @return 0 on success, -1 on error
    */

   int (*lc_lookup_typeqc)(const Land_Cover_Type *, unsigned int,
                           const double *, const double *,
                           unsigned char *);
   /**<  Nearest-neighbor lookup of MODIS Type QC land cover mask
    * @param[in]  lc   pointer to Land_Cover_Type object from land_cover_init
    * @param[in]  num  number of lon,lat points to process
    * @param[in]  lon  pointer to longitude coordinates
    * @param[in]  lat  pointer to latitude coordinates
    * @param[out] mask  pointer to nearest-neighbor land cover mask value
    * @return 0 on success, -1 on error
    */

#ifdef LAND_COVER_PRIVATE_DATA
   LAND_COVER_PRIVATE_DATA
#endif
};

/**  Initialize Land_Cover_Type object
 * @param[in] cfg   pointer to config_t object associated with
 *                  \c L1_inr_post parameter file
 * @return pointer to initialized Land_Cover_Type on success, NULL on error
 */
extern Land_Cover_Type *land_cover_init (config_t *cfg);

#endif
