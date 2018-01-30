#ifndef __INTERP_POLY_INCLUDE__
#define __INTERP_POLY_INCLUDE__ 1

typedef struct Poly_Type Poly_Type;

typedef struct
{
   double xexp;
   size_t num_coef;
}
Poly_Config_Type;

extern void poly_free (Poly_Type *pt);

extern Poly_Type *
poly_interpol (size_t n, const double *x, const double *y,
               const Poly_Config_Type *config);

extern int
poly_eval (Poly_Type *pt, size_t n, const double *x, double *yest);

extern int poly_get_coeffs (Poly_Type *pt, size_t num_coeffs, double *coeffs);
extern int poly_set_coeffs (Poly_Type *pt, size_t num_coeffs, const double *coeffs);

#endif
