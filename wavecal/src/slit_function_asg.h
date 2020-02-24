#ifndef __SLIT_FUNCTION_ASG_INCLUDE_H__
#define __SLIT_FUNCTION_ASG_INCLUDE_H__ 1

extern int asg_normed_plus_derivs (const double *x, size_t nx, double *params,
                                   double *value,
                                   double *param_step,
                                   double *param_derivs[3]);

#endif
