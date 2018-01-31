#ifndef __WAVECAL_INCLUDE__
#define __WAVECAL_INCLUDE__ 1

#include <libconfig.h>

typedef struct Wavecal_Type Wavecal_Type;

typedef struct
{
   double fill_value;
   int *index_lim;
}
Wavecal_Config_Type;

typedef struct
{
   const double *wave;
}
Wavecal_Result_Type;

extern void wavecal_close (Wavecal_Type *wct);
extern Wavecal_Type *wavecal_open (config_t *cfg, int mode);
/* (mode = 0) => irradiance
 * (mode = 1) => radiance */

extern int wavecal_fit (Wavecal_Type *wct, int xtrack,
                        int num_wave, const double *wave,
                        const double *spec, const double *specerr,
                        const Wavecal_Config_Type *config,
                        Wavecal_Result_Type *result);

#endif
