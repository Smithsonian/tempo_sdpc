#ifndef __PLAN_SCAN_H__
#define __PLAN_SCAN_H__ 1

#include <stdint.h>
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

   uint16_t (*st_scan_type)(const Scan_Type *);
   /**< retrieve the scan type value */

   double (*st_scan_duration)(const Scan_Type *, int);
   /**< compute the time [sec] required to complete a scan with N steps */

   double (*st_integration_time)(const Scan_Type *);
   /**< retrieve the integration time [sec] for a single exposure in a co-add */

   double (*st_step_size)(const Scan_Type *);
   /**< retrieve the step size [microradian] for each scan step */

   int (*st_scan_num_steps)(const Scan_Type *);
   /**< retrieve the requested number of scan steps */

   double (*st_min_sun_angle)(const Scan_Type *);
   /**< retrieve the minimum allowed sun angle [deg] */

   int (*st_scan_beg_angle)(const Scan_Type *, double *, double *);
   /**< retrieve the (mirror_x,mirror_y) coordinates of the scan's eastern limit [urad] */

   int (*st_scan_end_angle)(const Scan_Type *, double *, double *);
   /**< retrieve the (mirror_x,mirror_y) coordinates of the scan's western limit [urad] */

   int (*st_scan_beg)(const Scan_Type *, double *, double *);
   /**< retrieve the (lon,lat) coordinates of the scan's eastern limit [deg] */

   int (*st_scan_end)(const Scan_Type *, double *, double *);
   /**< retrieve the (lon,lat) coordinates of the scan's western limit [deg] */

   int (*st_scan_day_beg)(const Scan_Type *, double *, double *);
   /**< retrieve the (lon,lat) coordinates of the day-begin control point [deg] */

   int (*st_scan_day_end)(const Scan_Type *, double *, double *);
   /**< retrieve the (lon,lat) coordinates of the day-end control point [deg] */

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
extern Scan_Type *scan_open (config_t *cfg, uint16_t scan_type);

typedef struct Twilight_Scan_Type Twilight_Scan_Type;

struct Twilight_Scan_Type
{
   void (*tst_delete)(Twilight_Scan_Type *);
   /**< delete an object of type \c Twilight_Scan_Type */

   int (*tst_twilight_scan_region_angles)(const Twilight_Scan_Type *, int, double *, double *, double *, int *);
   /**< retrieve twilight scan region as (mirror_x,mirror_y) of one boundary [urad], and an eastward or westward extent [urad]  */

   int (*tst_twilight_scan_region)(const Twilight_Scan_Type *, int, double *, double *, double *, int *);
   /**< retrieve twilight scan region as (lon,lat) of one boundary [deg], and an eastward or westward extent [urad]  */

   double (*tst_twilight_scan_duration)(const Twilight_Scan_Type *, int);
   /**< compute the time [sec] required to complete a scan with N steps */

   double (*tst_twilight_integration_time)(const Twilight_Scan_Type *);
   /**< retrieve the integration time [sec] for a single twilight exposure in a co-add */

#ifdef TWILIGHT_SCAN_TYPE_PRIVATE_DATA
   TWILIGHT_SCAN_TYPE_PRIVATE_DATA
#endif
};

extern Twilight_Scan_Type *twilight_scan_open (config_t *cfg);

typedef struct Split_Scan_Type Split_Scan_Type;

enum
{
   SCAN_SPLIT_STD = 0,
   SCAN_SPLIT_OPT1 = 1
};

struct Split_Scan_Type
{
   void (*sst_delete)(Split_Scan_Type *);

   int (*sst_scan_region_angles)(const Split_Scan_Type *, double *, double *, double *, double *);
   /* beg_x, beg_y, end_x, end_y */

   int (*sst_scan_region)(const Split_Scan_Type *, double *, double *, double *, double *);
   /* beg_lon, beg_lat, end_lon, end_lat */

   double (*sst_scan_integration_time) (const Split_Scan_Type *);

   int (*sst_base_scan_method)(const Split_Scan_Type *);
   int (*sst_num_repeats_cbm)(const Split_Scan_Type *);
   double (*sst_weight)(const Split_Scan_Type *);

#ifdef SPLIT_SCAN_TYPE_PRIVATE_DATA
   SPLIT_SCAN_TYPE_PRIVATE_DATA
#endif
};

extern Split_Scan_Type *split_scan_open (config_t *cfg, const char *scan_method);

typedef struct
{
   double jd_utc_beg;
   /**< Earliest time when radiance scanning could begin [days] */

   double jd_utc_end;
   /**< Time when radiance scanning should end [days] */

   double jd_utc_beg_full;
   /**< Time when ROI max illumination begins [days] */

   double jd_utc_end_full;
   /**< Time when ROI max illumination ends [days] */

   double jd_utc_beg_safe;
   /**< Earliest time when the aperture may safely open [days] */

   double jd_utc_end_safe;
   /**< Time when the aperture must close [days] */
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

/** Compute the time at which an irradiance measurement may be taken with a specified sun angle
 * @param[in]  sgt     Pointer to an initialized \c Solar_Geom_Type object
 * @param[in]  irr_sun_angle   Desired angle of incidence of sunlight [deg]
 * @param[in]  after_midnight  Non-zero selects times after midnight.
 * @param[in]  jd_utc  Julian date of local midnight on the day of interest.
 * @param[out] jd_utc_irr  Julian date for the irradiance observation
 * @return 0 on success, -1 on error
 */
extern int scan_irradiance_time (Solar_Geom_Type *sgt, double irr_sun_angle, int after_midnight,
                                 double jd_utc, double *jd_utc_irr);
#endif
