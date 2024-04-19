#ifndef __UTIL_INCLUDE__
#define __UTIL_INCLUDE__ 1

extern int bsearch_f (float t, const float *x, int n);
extern int bsearch_d (double t, const double *x, int n);
extern double *alloc_doubles (int n);
extern int find_x (double x, const double *a, int na);
extern char *expand_string (const char *s);
extern char *path_concat (const char *dir, const char *basename);

#include <libconfig.h>
extern int read_config_float_array (config_setting_t *s, const char *name,
                                    double **pa, size_t *pnum_a);

extern int enable_state_define (config_t *cfg, const char *state_name);
extern const char *enable_state_query_enum (const char *name);
extern int enable_state_query_bool (const char *name);

#include "enable_states.h"

#include <tio.h>
extern int meta_record_basename (TIO_Meta_Type *meta, const char *path);

#include "image.h"
typedef struct Trend_File_Type Trend_File_Type;
typedef struct Trend_Record_Type Trend_Record_Type;

extern Trend_File_Type *trend_collect_open (const char *trend_file, int exposure_type);
extern int trend_collect_close (Trend_File_Type *tft);

extern Trend_Record_Type *trend_collect_new_record (Trend_File_Type *tft);
extern void trend_collect_free_record (Trend_Record_Type *tr);
extern Trend_Record_Type *trend_collect_set_active_record (Trend_Record_Type *tr);

extern int trend_collect_write_record (const Trend_Record_Type *tr);

extern int trend_collect_time (double start_time, int index);
extern int trend_collect_eoffsets (const float *eoffsets, const int *phase_change);
extern int trend_collect_temp (float spec_temp, float tele_temp, float bench_temp);
extern int trend_collect_gain (float fpa_temp, float fpe_temp, const float *gain);
extern int trend_collect_readnoise (const float *readnoise);
extern int trend_collect_sdc (int num_dg_rows, int num_tg_rows, const float *sdc);
extern int trend_collect_dc_mean (const float *mean_dc, const float *stddev_dc);
extern int trend_collect_solar_angles (double solar_theta, double solar_phi, int use_reference_diffuser);
extern int trend_collect_pqf (const Image_Type *img);
extern int trend_collect_pqf_uv (Image_Pqf_Bitmap_Type *img, int num_rows, int num_cols);
extern int trend_collect_pqf_vis (Image_Pqf_Bitmap_Type *img, int num_rows, int num_cols);

#endif
