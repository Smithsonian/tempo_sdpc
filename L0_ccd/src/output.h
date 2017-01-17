#ifndef __TEMPO_OUTPUT_H__
#define __TEMPO_OUTPUT_H__ 1

#include <libconfig.h>
#include "granule.h"

typedef struct
{
   const Granule_Exprec_Type *exprec;
   const Image_Type *img_err;
}
Output_Exprec_Type;

typedef struct Output_Type Output_Type;

struct Output_Type
{
   int (*out_set_file)(Output_Type *, const char *);
   int (*out_set_dims)(Output_Type *, int, int, int);
   int (*out_create)(Output_Type *);
   int (*out_write_rec)(Output_Type *, int, const Output_Exprec_Type *);
   int (*out_copy_metadata)(Output_Type *, int);
   int (*out_close)(Output_Type *);
   void (*out_free)(Output_Type *);

#ifdef OUTPUT_PRIVATE_DATA
   OUTPUT_PRIVATE_DATA
#endif
};

extern Output_Type *output_alloc (config_t *cfg, int exposure_type);

#endif
