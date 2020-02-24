#ifndef __SLIT_FUNCTION_INCLUDE_H__
#define __SLIT_FUNCTION_INCLUDE_H__ 1

#define SFT_NUM_PARAMS (3)

typedef struct Slit_Function_Type Slit_Function_Type;

extern void sft_free (Slit_Function_Type *);

extern Slit_Function_Type *sft_init (double dx, size_t num_sf, size_t num_waves,
                                     double *params, double *param_step);

extern int sft_apply (Slit_Function_Type *sft, const double *spec, double *spec_convolved,
                      int compute_derivs, double *spec_derivs_convolved[SFT_NUM_PARAMS]);

#endif
