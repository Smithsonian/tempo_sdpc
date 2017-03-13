#ifndef __TEMPO_LAND_COVER_H__
#define __TEMPO_LAND_COVER_H__ 1

#include <libconfig.h>

typedef struct Land_Cover_Type Land_Cover_Type;

struct Land_Cover_Type
{
   void (*lc_delete)(Land_Cover_Type *);
   int (*lc_lookup_type1)(const Land_Cover_Type *, unsigned int,
                          const double *, const double *,
                          unsigned char *);
   int (*lc_lookup_typeqc)(const Land_Cover_Type *, unsigned int,
                           const double *, const double *,
                           unsigned char *);

#ifdef LAND_COVER_PRIVATE_DATA
   LAND_COVER_PRIVATE_DATA
#endif
};

extern Land_Cover_Type *land_cover_init (config_t *cfg);

#endif
