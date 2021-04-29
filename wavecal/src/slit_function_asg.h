#ifndef __SLIT_FUNCTION_ASG_INCLUDE_H__
#define __SLIT_FUNCTION_ASG_INCLUDE_H__ 1

/* When norm is unknown, set norm <= 0.0 */
extern int asg_normed_plus_derivs (const double *x, size_t nx, double *params, double norm,
                                   double *value,
                                   double *param_step,
                                   double *param_derivs[3]);

extern int asg_compute_norm (double *params, double *norm);

#endif
