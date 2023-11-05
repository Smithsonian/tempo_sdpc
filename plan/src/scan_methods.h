#ifndef __PLAN_SCAN_METHODS_H__
#define __PLAN_SCAN_METHODS_H__ 1

/** @file scan_methods.h
 *  @brief Support different scan planning methods
 */

#include "scan.h"
#include "plan_list.h"
#include "vis.h"

#define SMA_MAX_CALIBRATED_MIRROR_X  49600.0
#define SMA_MAX_CALIBRATED_MIRROR_Y   4400.0
/* The maximum range over which the scan mirror should be commanded is:
 *    |X| <= SMA_MAX_CALIBRATED_MIRROR_X,
 *    |Y| <= SMA_MAX_CALIBRATED_MIRROR_Y.
 * These coordinates refer to the scan mirror tilt angle in microradians.
 * (From TEMPO ConOps, Ball doc 2418231, Rev E, section 12.3.2, page 75,
 *  10/12/2017)
 */

typedef struct
{
   Plan_List_Type *(*sm_plan)(const Scan_Type *, Solar_Geom_Type *,
                              const Scan_Limit_Times_Type *, void *);
   /**< generate a scan plan */

   int (*sm_vis)(Vis_Type *v, const Plan_List_Type *, double, int, double *, int);
   /**< compute solar zenith angle maps and write them to the specified netCDF file */
}
Scan_Method_Type;

/** Find the named scan method
 * @param[in] name   Name of the scan method
 * @return on success, a pointer to a Scan_Method_Type object;
 *         on error, a NULL pointer
 */
extern const Scan_Method_Type *find_scan_method (const char *name);

/** Use a lookup table to invert the mapping from FOR scan angles (x,y)
 *  to surface coordinates (lon,lat)
 */
extern int scan_xy_to_lonlat (const double *x_urad, const double *y_urad, int n,
                              double *lon_deg, double *lat_deg, double sat_lon);
extern void scan_set_lonlat_bounding_box (double lon_min, double lon_max,
                                          double lat_min, double lat_max, double yoffset);

extern int Set_Geometry_Params (double ewbias, double nsbias, double clockingbias, double telescopeOffset);
#endif
