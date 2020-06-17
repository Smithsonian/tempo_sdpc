#ifndef __L1_INR_PREP_RADIANCE__
#define __L1_INR_PREP_RADIANCE__ 1
/** @file radiance.h
 *  @brief Manage radiance file I/O
 */

#include <libconfig.h>
#include "row_select.h"
#include "ephem.h"

typedef struct Radiance_Type Radiance_Type;

struct Radiance_Type
{
   int ncid;               /**< open netcdf file id */
   char *file;             /**< file path */

#ifdef RADIANCE_PRIVATE_DATA
   RADIANCE_PRIVATE_DATA
#endif
};

/** Free memory allocated by @c radiance_open or @c radiance_create
 *
 * @param[in] r   Pointer to @c Radiance_Type struct
 */
extern void radiance_close (Radiance_Type *r);

/** Open an existing Level 1 radiance file
 *
 * @param[in] file   Path to existing Level 1 radiance file
 * @return On sucess, a pointer to a @c Radiance_Type struct
 *         On error, a NULL pointer
 */
extern Radiance_Type *radiance_open (const char *file);

/** Create a telemetry-only Level 1 radiance file for INR
 *
 * @param[in] file   Temporary path to Level 1 radiance file
 *                   to be created
 * @param[in] processing_version  Processing version number for
 *                                the new Level 1 radiance file
 * @return On sucess, a pointer to a @c Radiance_Type struct
 *         On error, a NULL pointer
 *
 * The purpose of a telemetry-only Level 1 radiance file is
 * to provide gyroscope (IRU) and scan mechanism (SMC) telemetry,
 * and ephemeris data for INR processing during time intervals
 * during time intervals when radiance spectra are not being acquired.
 */
extern Radiance_Type *radiance_create (const char *file,
                                       int processing_version);

/** Update header coverage time keywords in telemetry-only radiance file
 *
 *  @param[in] r      Pointer to a @c Radiance_Type structure
 *                    associated with a telemetry-only radiance file.
 *  @param[in] tstart Effective coverage start time in TAI seconds since the TEMPO epoch
 *  @param[in] tstop  Effective coverage stop  time in TAI seconds since the TEMPO epoch
 *  @return 0 on success, -1 on error
 *
 * In a telemetry-only radiance file the header keywords,
 * @c time_coverage_start, and @c time_coverage_end, are determined
 * by the processing context.  The endpoints of the IRU, SMC and ephemeris
 * time series that are stored in the file may extend a little beyond the time
 * interval determined from the processing context.  This padding is intended
 * to ensure that the INR software receives complete time coverage of each
 * relevant data stream.  Note that the precise values
 * of the actual time series endpoints may differ slightly from the endpoints
 * of the time interval specified on the command line when @c L1_inr_prep
 * is run.
 */
extern int radiance_update_coverage_times (Radiance_Type *r, double tstart, double tstop);

/** Query the coverage time interval associated with a @c Radiance_Type struct
 *
 * @param[in] r         Pointer to a @c Radiance_Type struct
 * @param[out] tstart   Beginning of the coverage time interval, expressed
 *                      as the number of seconds elapsed since the TEMPO epoch.
 * @param[out] tstop    End of the coverage time interval, expressed
 *                      as the number of seconds elapsed since the TEMPO epoch.
 * @return 0 on success, -1 on error
 */
extern int radiance_interval (Radiance_Type *r,
                              double *tstart, double *tstop);

/** Copy gyroscope (IRU) time series data into a Level 1 radiance file
 *
 * @param[in] r    Pointer to a @c Radiance_Type struct associated with an
 *                 an open Level 1 radiance file
 * @param[in] meta Pointer to @c TIO_Meta_Type struct associated with an open
 *                 Level 1 radiance file
 * @param[in] rst  Pointer to a @c Row_Select_Type object referencing
 *                 the relevant gyroscope time series data
 * @return 0 on success, -1 on error
 */
extern int radiance_copy_iru (Radiance_Type *r, TIO_Meta_Type *meta,
                              const Row_Select_Type *rst);

/** Copy scan mechanism controller (SMC) time series data into a Level 1 radiance file
 *
 * @param[in] r    Pointer to a @c Radiance_Type struct associated with an
 *                 an open Level 1 radiance file
 * @param[in] meta Pointer to @c TIO_Meta_Type struct associated with an open
 *                 Level 1 radiance file
 * @param[in] rst  Pointer to a @c Row_Select_Type object referencing
 *                 the relevant SMC time series data
 * @return 0 on success, -1 on error
 */
extern int radiance_copy_smc (Radiance_Type *r, TIO_Meta_Type *meta,
                              const Row_Select_Type *rst);

/** Copy ephemeris time series data into a Level 1 radiance file
 *
 * @param[in] r    Pointer to a @c Radiance_Type struct associated with an
 *                 an open Level 1 radiance file
 * @param[in] rst  Pointer to an @c Eph_Type object containing
 *                 the relevant ephemeris time series data
 * @return 0 on success, -1 on error
 */
extern int radiance_write_eph (Radiance_Type *r,
                               const Eph_Type *eph);

#endif
