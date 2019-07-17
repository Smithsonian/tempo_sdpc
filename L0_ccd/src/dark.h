#ifndef __TEMPO_DARK_INCLUDE_H__
#define __TEMPO_DARK_INCLUDE_H__ 1

#include <libconfig.h>
#include "image.h"

typedef struct Dark_Type Dark_Type;

typedef struct
{
   float exposure_time;
   float fpa_temp;
   float *storage_region_dark;
}
Dark_Lookup_Type;

struct Dark_Type
{
   void (*drk_close)(Dark_Type *);

   int (*drk_open)(Dark_Type *, const char *path);
   /* path may be either a dark granule, or a lookup table */

   int (*drk_get_image)(const Dark_Type *, const Dark_Lookup_Type *, Image_Type *img);

#ifdef DARK_PRIVATE_DATA
   DARK_PRIVATE_DATA
#endif
};

extern Dark_Type *drk_init (config_t *cfg);

#endif
