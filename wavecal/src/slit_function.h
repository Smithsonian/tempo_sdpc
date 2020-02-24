#ifndef __SLIT_FUNCTION_INCLUDE_H__
#define __SLIT_FUNCTION_INCLUDE_H__ 1

#define SFT_NUM_PARAMS (3)

typedef struct Slit_Function_Type Slit_Function_Type;
typedef int SFT_Eval_Type (const double *x, size_t nx, double *params,
                           double *value, double *param_step, double *param_derivs[SFT_NUM_PARAMS]);

extern void sft_free (Slit_Function_Type *);

extern Slit_Function_Type *sft_new (size_t num_sf, size_t num_waves);

extern int sft_config (Slit_Function_Type *sft, SFT_Eval_Type *sf_eval,
                       double dx, double *param_step, double *params);

extern int sft_apply (Slit_Function_Type *sft, const double *spec, double *spec_convolved,
                      int compute_derivs, double *spec_derivs_convolved[SFT_NUM_PARAMS]);

#endif
