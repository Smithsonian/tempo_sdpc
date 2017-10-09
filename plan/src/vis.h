#ifndef __PLAN_VIS_H__
#define __PLAN_VIS_H__ 1

/** @file vis.h
 *  @brief Map solar zenith angle vs position to help visualize plan
 */

#include <libconfig.h>

typedef struct Vis_Type Vis_Type;

/** Free a \c Vis_Type object
 * @param[in]  v   pointer to a \c Vis_Type object allocated by \ref vis_init
 */
extern void vis_free(Vis_Type *v);

/** Initialize a \c Vis_Type object
 * @param[in]  cfg   pointer to an open configuration file
 * @param[in]  sgt   pointer to an initialized object of type \ref Solar_Geom_Type
 * @return on success, an initialized \c Vis_Type object
 */
extern Vis_Type *vis_init (config_t *cfg, Solar_Geom_Type *sgt);

/** Compute a solar zenith angle map for the specified Julian date.
 * @param[in]  v        Pointer to a \c Vis_Type object initialized by \c vis_init
 * @param[in]  jd_utc   The Julian date
 * @param[in]  psza      Optional pointer to an array of doubles large enough
 *                      to hold a square image of size \c img_size.
 *                      If NULL, then an array will be allocated
 * @return on success, either a pointer to an allocated array (if psza=NULL)
 *         or the input pointer, \c psza;
 *         on error, a NULL pointer.
 */
extern double *vis_sza (const Vis_Type *v, double jd_utc, double *psza);

/** Write lon,lat grid coordinates to a netCDF file
 * @param[in]  v       Pointer to a \c Vis_Type object initialized by \c vis_init
 * @param[in]  ncid    netCDF file identifier
 * @return 0 on success, -1 on error
 */
extern int vis_write_grid (Vis_Type *v, int ncid);

/** Write a solar zenith angle map to a named variable in a netCDF file
 * @param[in]  v       Pointer to a \c Vis_Type object initialized by \c vis_init
 * @param[in]  ncid    netCDF file identifier
 * @param[in]  jd_utc  Julian date of the SZA map
 * @param[in]  name    netCDF file variable name for the SZA map
 * @param[in]  value   Pointer to the SZA map array
 * @return 0 on success, -1 on error
 */
extern int vis_write_value (const Vis_Type *v, int ncid, double jd_utc,
                            const char *name, const double *value);

#endif
