#ifndef __UTIL_INCLUDE__
#define __UTIL_INCLUDE__ 1

extern int bsearch_f (float t, const float *x, int n);
extern int bsearch_d (double t, const double *x, int n);
extern double *alloc_doubles (int n);
extern int find_x (double x, const double *a, int na);
extern char *expand_path (const char *path);
extern char *path_concat (const char *dir, const char *basename);

#include <libconfig.h>
extern int read_config_float_array (config_setting_t *s, const char *name,
                                    double **pa, size_t *pnum_a);

#include <tio.h>
extern int meta_record_basename (TIO_Meta_Type *meta, const char *path);

#endif
