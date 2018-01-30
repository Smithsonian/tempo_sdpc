#ifndef __INTERP_CHEBYSHEV_INCLUDE__
#define __INTERP_CHEBYSHEV_INCLUDE__ 1

typedef struct Cheb_Type Cheb_Type;

typedef struct
{
   double xmin;
   double xmax;
   size_t num_coef;
}
Cheb_Config_Type;

extern void cheb_free (Cheb_Type *ct);

extern Cheb_Type *
cheb_interpol (size_t n, const double *x, const double *y,
               const Cheb_Config_Type *config);

extern int
cheb_eval (Cheb_Type *ct, size_t n, const double *x, double *yest);

extern int cheb_get_coeffs (Cheb_Type *ct, size_t num_coeffs, double *coeffs);
extern int cheb_set_coeffs (Cheb_Type *ct, size_t num_coeffs, const double *coeffs);

extern int cheb_clenshaw_eval (double a, double b, const double *c,
                               size_t m, double x, double *value);

#endif
