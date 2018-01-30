#ifndef __SHAPEFUN_INCLUDE__
#define __SHAPEFUN_INCLUDE__ 1

#include <libconfig.h>

typedef struct
{
   double *x;
   double *y;
   size_t n;
   int malloced;
}
Shapefun_Init_Type;

typedef struct Shapefun_Type Shapefun_Type;

struct Shapefun_Type
{
   void (*st_delete)(Shapefun_Type *);
   int (*st_method_id)(const Shapefun_Type *);
   int (*st_num_params)(const Shapefun_Type *);
   int (*st_init_params)(const Shapefun_Type *, const Shapefun_Init_Type *,
                         size_t, double *);
   int (*st_eval)(const Shapefun_Type *, const double *,
                  size_t, const double *, double *);

   double *node_x;      /**< fixed, ordered set of x coordinates */
   size_t num_nodes;
   int malloced_node_x;

   size_t num_coef;
   double xmin;
   double xmax;
   double xexp;

#ifdef SHAPEFUN_PRIVATE_DATA
   SHAPEFUN_PRIVATE_DATA
#endif
};

extern Shapefun_Type *shapefun_create (const char *method_name);

#endif
