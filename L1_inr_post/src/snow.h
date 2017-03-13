#ifndef __TEMPO_SNOW_ICE_H__
#define __TEMPO_SNOW_ICE_H__ 1

typedef struct Snow_Type Snow_Type;

struct Snow_Type
{
   void (*sn_delete)(Snow_Type *);
   int (*sn_lookup)(const Snow_Type *, unsigned int,
                    const double *, const double *, unsigned char *);

#ifdef SNOW_TYPE_PRIVATE_DATA
   SNOW_TYPE_PRIVATE_DATA
#endif
};

extern Snow_Type *snow_init (const char *file);

#endif
