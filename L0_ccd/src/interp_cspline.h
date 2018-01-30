#ifndef __INTERP_CSPLINE_INCLUDE__
#define __INTERP_CSPLINE_INCLUDE__ 1

typedef struct Cspline_Type Cspline_Type;
typedef struct Cspline_Config_Type Cspline_Config_Type;

struct Cspline_Config_Type
{
   int unused;
};

extern void cspline_free (Cspline_Type *ct);

extern Cspline_Type *
cspline_interpol (size_t n, const double *x, const double *y,
                  void *config);

extern int
cspline_eval (Cspline_Type *ct, size_t n, const double *x, double *yest);
#endif
