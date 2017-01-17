#ifndef __GRANULE_INCLUDE__
#define __GRANULE_INCLUDE__ 1

#include "image.h"

/* symbols defined for consistency with IOCSDPC definitions */
enum
{
   EXPREC_TYPE_RADIANCE   = 0,
   EXPREC_TYPE_DARK       = 1,
   EXPREC_TYPE_IRRADIANCE = 2,
   EXPREC_TYPE_LIN_IRR    = 3,
   EXPREC_TYPE_LIN_DARK   = 4,
   EXPREC_TYPE_UNKNOWN    = -1
};

#define EXPREC_TYPE_IS_LINEARITY(type) \
   (((type) == EXPREC_TYPE_LIN_DARK) || ((type) == EXPREC_TYPE_LIN_IRR))

#define EXPREC_TYPE_IS_DARK(type) \
   (((type) == EXPREC_TYPE_DARK) || ((type) == EXPREC_TYPE_LIN_DARK))

#define EXPREC_TYPE_IS_IRRADIANCE(type) \
   (((type) == EXPREC_TYPE_IRRADIANCE) || ((type) == EXPREC_TYPE_LIN_IRR))

typedef struct
{
   double start_time;
   double exposure_time;
   double frame_transfer_time;  /* [sec] frame transfer time */
   double readout_time;         /* [sec] storage region readout time */
   int exposure_type;
   int num_coadds;
   Image_Type *img;
}
Granule_Exprec_Type;

typedef struct Granule_Type Granule_Type;

struct Granule_Type
{
   void (*granule_close) (Granule_Type *);
   int (*granule_num_exprecs)(const Granule_Type *);
   Granule_Exprec_Type *(*granule_read_exprec_by_index) (const Granule_Type *, int, Granule_Exprec_Type **);
   void (*granule_free_exprec) (Granule_Exprec_Type *);

#ifdef GRANULE_TYPE_PRIVATE_DATA
   GRANULE_TYPE_PRIVATE_DATA
#endif
};

extern Granule_Type *granule_open (const char *file);

#endif
