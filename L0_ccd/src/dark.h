#ifndef __TEMPO_DARK_INCLUDE_H__
#define __TEMPO_DARK_INCLUDE_H__ 1

#include <libconfig.h>
#include "image.h"

typedef struct Dark_Type Dark_Type;

struct Dark_Type
{
   void (*drk_close)(Dark_Type *);

   int (*drk_open)(Dark_Type *, const char *path);
   /* path may be either a dark granule, or a lookup table */

   int (*drk_image)(const Dark_Type *, Image_Type *img);

   int (*drk_image_Tfpa_adj)(const Dark_Type *, float fpa_temp, Image_Type *img);

   int (*drk_image_sdc_adj)(const Dark_Type *, float *sdc, Image_Type *img);

#ifdef DARK_PRIVATE_DATA
   DARK_PRIVATE_DATA
#endif
};

extern Dark_Type *drk_init (config_t *cfg);

#endif
