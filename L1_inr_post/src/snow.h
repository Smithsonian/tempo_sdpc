#ifndef __TEMPO_SNOW_ICE_H__
#define __TEMPO_SNOW_ICE_H__ 1

/** @file snow.h
 *  @brief Load snow and ice cover data,
 *         perform coordinate-based nearest-neighbor lookup
 */

typedef struct Snow_Type Snow_Type;

struct Snow_Type
{
   void (*sn_delete)(Snow_Type *);
   /**<  Free resources associated with Snow_Type object
    * @param[in]  sn  pointer to Snow_Type object from \c snow_init
    */

   int (*sn_regrid)(Snow_Type *, unsigned int,
                    const double *, const double *, float *);
   /**<  Regrid snow and ice mask to derive the fractional area covered by snow and/or ice
    * @param[in]  sn   pointer to Snow_Type object from \c snow_init
    * @param[in]  num  number of lon,lat points to process
    * @param[in]  lon  pointer to longitude coordinates
    * @param[in]  lat  pointer to latitude coordinates
    * @param[out] frac  fraction of pixel area covered by snow and/or ice
    * @return 0 on success, -1 on error
    */

#ifdef SNOW_TYPE_PRIVATE_DATA
   SNOW_TYPE_PRIVATE_DATA
#endif
};

/**  Initialize Snow_Type object
 * @param[in] file  path to IMS product file
 * @return pointer to initialized Snow_Type on success, NULL on error
 */
extern Snow_Type *snow_init (const char *file);

#endif
