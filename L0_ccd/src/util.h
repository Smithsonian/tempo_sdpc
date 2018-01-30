#ifndef __UTIL_INCLUDE__
#define __UTIL_INCLUDE__ 1

extern int bsearch_d (double t, const double *x, int n);
extern double *alloc_doubles (int n);
extern int find_x (double x, const double *a, int na);

#include <libconfig.h>
extern int read_config_float_array (config_setting_t *s, const char *name,
                                    double **pa, size_t *pnum_a);

#endif
