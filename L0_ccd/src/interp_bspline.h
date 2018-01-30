#ifndef __INTERP_BSPLINE_INCLUDE__
#define __INTERP_BSPLINE_INCLUDE__ 1

typedef struct Bspline_Type Bspline_Type;

/* xmin==xmax  means use full data range
 * breakpts=NULL means spline through all data points
 *   (instead of least-squares spline in sub-intervals)
 */
typedef struct
{
   double xmin;
   double xmax;
   const double *breakpts;
   size_t num_breakpts;
}
Bspline_Config_Type;

extern void bspline_free (Bspline_Type *bt);

extern Bspline_Type *
bspline_interpol (size_t n, const double *x, const double *y,
                  const Bspline_Config_Type *config);

extern int
bspline_eval (Bspline_Type *bt, size_t n, const double *x, double *yest);

#endif
