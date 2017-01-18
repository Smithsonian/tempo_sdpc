#ifndef __SENSORCAL_INCLUDE__
#define __SENSORCAL_INCLUDE__ 1

#include <libconfig.h>
#include "image.h"

typedef struct Calibration_Type Calibration_Type;

struct Calibration_Type
{
   void (*cal_delete)(Calibration_Type *);
   int (*cal_apply_rcoeffs)(const Calibration_Type *, Image_Type *);
   int (*cal_apply_btdf)(const Calibration_Type *, double, double,
                         Image_Type *);
   int (*cal_apply_prnu)(const Calibration_Type *, Image_Type *);
   int (*cal_wavecal)(const Calibration_Type *, Image_Type *, Image_Type *);

#ifdef SENSORCAL_PRIVATE_DATA
   SENSORCAL_PRIVATE_DATA
#endif
};

extern Calibration_Type *sensorcal_init (config_t *cfg);

#endif
