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
   double (*st_scan_duration)(const Scan_Type *, int);
   double (*st_step_exposure)(const Scan_Type *);
   int (*st_step_size)(const Scan_Type *);
   double (*st_min_sun_angle)(const Scan_Type *);
   int (*st_scan_beg)(const Scan_Type *, double *, double *);
   int (*st_scan_end)(const Scan_Type *, double *, double *);
   int (*st_print_params)(const Scan_Type *, const char *, FILE *);

#ifdef SCAN_TYPE_PRIVATE_DATA
   SCAN_TYPE_PRIVATE_DATA
#endif
};

extern Scan_Type *scan_open (config_t *cfg);

typedef struct
{
   double jd_utc_beg;
   double jd_utc_end;
   double jd_utc_beg_full;
   double jd_utc_end_full;
}
Scan_Limit_Times_Type;

extern int scan_limit_times (const Scan_Type *st, double jd_utc,
                             Solar_Geom_Type *sgt,
                             Scan_Limit_Times_Type *slt);
#endif
