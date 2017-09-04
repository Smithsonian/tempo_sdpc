#ifndef __L1_INR_PREP_RADIANCE__
#define __L1_INR_PREP_RADIANCE__ 1

#include <libconfig.h>
#include "row_select.h"
#include "ephem.h"

typedef struct Radiance_Type Radiance_Type;

struct Radiance_Type
{
   int ncid;
   char *file;

#ifdef RADIANCE_PRIVATE_DATA
   RADIANCE_PRIVATE_DATA
#endif
};

extern void radiance_close (Radiance_Type *r);
extern Radiance_Type *radiance_open (const char *file);
extern Radiance_Type *radiance_create (const char *file,
                                       int processing_version);
extern int radiance_update_coverage_times (Radiance_Type *r);

extern int radiance_interval (Radiance_Type *r,
                              double *tstart, double *tstop);

extern int radiance_copy_iru (Radiance_Type *r,
                              const Row_Select_Type *rst);

extern int radiance_copy_smc (Radiance_Type *r,
                              const Row_Select_Type *rst);

extern int radiance_write_eph (Radiance_Type *r,
                               const Eph_Type *eph);

#endif
