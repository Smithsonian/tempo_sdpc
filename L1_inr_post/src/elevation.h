#ifndef __TEMPO_ELEVATION_H__
#define __TEMPO_ELEVATION_H__ 1

#include <libconfig.h>

typedef struct Elevation_Type Elevation_Type;

struct Elevation_Type
{
   void (*et_delete)(Elevation_Type *);
   int (*et_lookup)(const Elevation_Type *, unsigned int,
                    const double *, const double *, short *, short *);
   int (*et_regrid)(const Elevation_Type *, unsigned int,
                    const double *, const double *, short *, short *);

#ifdef ELEVATION_TYPE_PRIVATE_DATA
   ELEVATION_TYPE_PRIVATE_DATA
#endif
};

extern Elevation_Type *elevation_init (config_t *cfg);

#endif
