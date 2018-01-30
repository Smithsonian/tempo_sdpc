#ifndef __WCI_LSQ_INCLUDE__
#define __WCI_LSQ_INCLUDE__ 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct Lsq_Type Lsq_Type;
typedef struct Lsq_Method_Type Lsq_Method_Type;

struct Lsq_Method_Type
{
   int (*lsqm_set_row)(Lsq_Method_Type *m, gsl_matrix *X, size_t i, double xi);
   void *client_data;
};

extern void lsq_free (Lsq_Type *lsq);

extern Lsq_Type *lsq_alloc (size_t num_data, size_t num_coef,
                            Lsq_Method_Type *method);

extern int lsq_set_data (Lsq_Type *lsq, size_t n,
                         const double *x, const double *y);

extern int lsq_solve (Lsq_Type *lsq, gsl_vector *coef, double *chisqr);

#endif
