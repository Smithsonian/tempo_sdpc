#ifndef __L1_INR_PREP_EPHEM__
#define __L1_INR_PREP_EPHEM__ 1
/** @file ephem.h
 *  @brief Manage ephemeris subsetting
 */

typedef struct
{
   double *x;
   double *y;
   double *z;
}
Eph_Vector_Type;

typedef struct
{
   Eph_Vector_Type r;  /**< position in WGS84 cartesian coordinates [km] */
   Eph_Vector_Type v;  /**< velocity [km/sec] */
   double *t;          /**< time [sec since TEMPO epoch] */
   size_t n;           /**< number of ephemeris time samples */
   size_t num_alloc;   /**< number of samples allocated */
}
Eph_Type;

/** Free memory allocated by @c eph_read_subset
 *
 *  @param[in] eph  Pointer to @c Eph_Type struct with fields allocated by @c eph_read_subset
 */
extern void eph_free (Eph_Type *eph);

/** Read ephemeris data for a specified time interval
 *
 * @param[in] eph   Pointer to @c Eph_Type struct.  On successful
 *                  return, the fields of this structure will be
 *                  allocated and populated with ephemeris data.
 * @param[in] file  Path to the ephemeris data file, formatted
 *                  in accordance with the SOC/IOC IDD.
 * @param[in] time_beg   Beginning of the time interval for
 *                       which ephemeris data will be read in,
 *                       expressed as the number of seconds elapsed
 *                       since the TEMPO epoch.
 * @param[in] time_end   End of the time interval for
 *                       which ephemeris data will be read in,
 *                       expressed as the number of seconds elapsed
 *                       since the TEMPO epoch.
 * @param[in] num_pad    Number of additional ephemeris points to
 *                       include as padding in addition to the
 *                       primary time interval from @c time_beg to
 *                       @c time_end.
 * @return 0 on success, -1 on error.
 */
extern int eph_read_subset (Eph_Type *eph, const char *file,
                            double time_beg, double time_end,
                            int num_pad);

#endif
