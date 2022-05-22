#ifndef __PLAN_SCAN_METHODS_H__
#define __PLAN_SCAN_METHODS_H__ 1

/** @file scan_methods.h
 *  @brief Support different scan planning methods
 */

#include "scan.h"
#include "plan_list.h"
#include "vis.h"

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

#endif
