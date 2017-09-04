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

extern void eph_free (Eph_Type *eph);

extern int eph_read_subset (Eph_Type *eph, const char *file,
                            double time_beg, double time_end,
                            int num_pad);

#endif
