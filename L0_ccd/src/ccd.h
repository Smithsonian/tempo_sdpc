#ifndef __TEMPO_CCD_INCLUDE__
#define __TEMPO_CCD_INCLUDE__

#include <libconfig.h>
#include "image.h"

typedef struct CCD_Type CCD_Type;

struct CCD_Type
{
   void (*ccd_delete)(CCD_Type *);
   int (*ccd_correct_coadd)(const CCD_Type *, int, Image_Type *);
   int (*ccd_correct_offset)(const CCD_Type *, Image_Type *);
   int (*ccd_correct_nonlinearity)(const CCD_Type *, Image_Type *);
   int (*ccd_correct_gain)(const CCD_Type *, Image_Type *);
   int (*ccd_correct_smear)(const CCD_Type *, const void *, Image_Type *);
   int (*ccd_mean_storage_region_dark)(const CCD_Type *, const Image_Type *, float mean_sdc[4]);

   Image_Type *(*ccd_select_active_pixels)(const CCD_Type *, const Image_Type *);
   int (*ccd_update_noisesq)(const CCD_Type *, const float *, Image_Type *);

#ifdef CCD_TYPE_PRIVATE_DATA
   CCD_TYPE_PRIVATE_DATA
#endif
};

extern CCD_Type *ccd_init (config_t *cfg);

#endif
