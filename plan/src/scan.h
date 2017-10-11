#ifndef __PLAN_SCAN_H__
#define __PLAN_SCAN_H__ 1

#include <libconfig.h>
#include "solar.h"

/** @file scan.h
 *  @brief Instrument scan characteristics and timing
 */

typedef struct Scan_Type Scan_Type;

struct Scan_Type
{
   void (*st_delete)(Scan_Type *);
   /**< delete an object of type \c Scan_Type */

   double (*st_scan_duration)(const Scan_Type *, int);
   /**< compute the time [sec] required to complete a scan with N steps */

   double (*st_integration_time)(const Scan_Type *);
   /**< retrieve the integration time [sec] for a single exposure in a co-add */

   int (*st_step_size)(const Scan_Type *);
   /**< retrieve the step size [microradian] for each scan step */

   double (*st_min_sun_angle)(const Scan_Type *);
   /**< retrieve the minimum allowed sun angle [deg] */

   int (*st_scan_beg)(const Scan_Type *, double *, double *);
   /**< retrieve the (lon,lat) coordinates of the scan's eastern limit [deg] */

   int (*st_scan_end)(const Scan_Type *, double *, double *);
   /**< retrieve the (lon,lat) coordinates of the scan's western limit [deg] */

   int (*st_print_params)(const Scan_Type *, const char *, FILE *);
   /**< print the scan parameters to an open FILE pointer */

#ifdef SCAN_TYPE_PRIVATE_DATA
   SCAN_TYPE_PRIVATE_DATA
#endif
};

/** Initialize a \c Scan_Type object
 * @param[in]  cfg   Pointer to an open configuration file
 * @return on success, a pointer to an initialized \c Scan_Type object;
 *         on error, a NULL pointer.
*/
extern Scan_Type *scan_open (config_t *cfg);

typedef struct
{
   double jd_utc_beg;
   /**< Earliest time when scanning could begin [days] */

   double jd_utc_end;
   /**< Time when scanning should end [days] */

   double jd_utc_beg_full;
   /**< Time when max ROI illumination begins [days] */

   double jd_utc_end_full;
   /**< Time when max ROI illumination ends [days] */
}
Scan_Limit_Times_Type;

/** Compute the fiducial times of day that may limit scanning activities
 * @param[in]  st      Pointer to an initialized \c Scan_Type object
 * @param[in]  jd_utc  Julian date of local midnight at the beginning of
 *                     the day of interest.
 * @param[in]  sgt     Pointer to an initialized \c Solar_Geom_Type object
 * @param[out] slt     Pointer to the destination Scan_Limit_Times_Type object
 * @return 0 on success, -1 on error
*/
extern int scan_limit_times (const Scan_Type *st, double jd_utc,
                             Solar_Geom_Type *sgt,
                             Scan_Limit_Times_Type *slt);
#endif
