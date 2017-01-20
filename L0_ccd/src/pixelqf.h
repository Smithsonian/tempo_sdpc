#ifndef __TEMPO_PIXELQF_INCLUDE__
#define __TEMPO_PIXELQF_INCLUDE__

#include <libconfig.h>
#include "image.h"

typedef struct Pixelqf_Type Pixelqf_Type;

struct Pixelqf_Type
{
   void (*pqf_delete)(Pixelqf_Type *);
   int (*pqf_flag_hotcold)(const Pixelqf_Type *, Image_Type *);
   int (*pqf_flag_neighbor)(const Pixelqf_Type *, Image_Type *,
                            Image_Pqf_Bitmap_Type, Image_Pqf_Bitmap_Type);
   int (*pqf_flag_transients)(const Pixelqf_Type *,
                              const Image_Pqf_Bitmap_Type *,
                              const Image_Type *, Image_Type *);

#ifdef PIXELQF_TYPE_PRIVATE_DATA
   PIXELQF_TYPE_PRIVATE_DATA
#endif
};

extern Pixelqf_Type *pixelqf_init (config_t *cfg);

#endif
