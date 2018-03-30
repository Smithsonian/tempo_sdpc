#ifndef __POLCORR_LPS_INC__
#define __POLCORR_LPS_INC__ 1

typedef struct Lps_Type Lps_Type;

extern int
lps_eval (Lps_Type *lps, int band_index, int xtrack,
          double lon, double lat, int n, const double *wave,
          double *linear_polarization_sensitivity,
          double *angle_of_max_transmission,
          double *lmp_irp_angle);

extern Lps_Type *lps_open (config_t *cfg);
extern void lps_close (Lps_Type *);

#endif
