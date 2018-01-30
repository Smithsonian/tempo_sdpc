#ifndef __INTERP_INCLUDE__
#define __INTERP_INCLUDE__ 1

enum
{
   INTERP_TYPE_CSPLINE = 0,
   INTERP_TYPE_CHEB = 1,
   INTERP_TYPE_POLY = 2
};

typedef struct Interp_Type Interp_Type;

struct Interp_Type
{
   void (*it_delete)(Interp_Type *);
   int (*it_method_id)(const Interp_Type *);
   int (*it_interp_init)(Interp_Type *, size_t, const double *, const double *);
   int (*it_interp_eval)(Interp_Type *, size_t, const double *, double *);

   size_t num_coef;
   double xmin;
   double xmax;
   double xexp;
   const double *breakpts;
   size_t num_breakpts;

#ifdef INTERP_TYPE_PRIVATE_DATA
   INTERP_TYPE_PRIVATE_DATA
#endif
};

extern Interp_Type *interp_create (const char *method_name);

#endif
