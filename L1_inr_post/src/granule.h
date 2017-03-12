#ifndef __TEMPO_GRANULE_INCLUDE__
#define __TEMPO_GRANULE_INCLUDE__ 1

#include "elevation.h"
#include "snow.h"
#include "land_cover.h"

typedef struct Granule_Type Granule_Type;

struct Granule_Type
{
   void (*gt_close)(Granule_Type *);
   int (*gt_set_elevation)(Granule_Type *, const Elevation_Type *);
   int (*gt_set_object_angles)(Granule_Type *);
   int (*gt_set_ground_pixel_flags)(Granule_Type *, double, double,
                                    const Snow_Type *,
                                    const Land_Cover_Type *);

#ifdef GRANULE_PRIVATE_DATA
   GRANULE_PRIVATE_DATA
#endif
};

extern Granule_Type *granule_open (const char *file);

#endif
