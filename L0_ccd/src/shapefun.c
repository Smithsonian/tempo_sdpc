/* -*- mode: C; mode: fold -*- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <gsl/gsl_poly.h>

#include "config.h"
#include "interp_cspline.h"
#include "interp_cheb.h"
#include "interp_poly.h"
#include "interp.h"
#include "util.h"

/**< An adjustable parameter is associated with each shapefun node.
 * Storage for the parameter vector is provided by the calling
 * routine. The initial parameter values are derived from
 * the input function (x0,y0) in a way that depends on
 * the interpolation method and on how many (x0,y0) initialization
 * points are provided.
 * The scaling function will be evaluated at a fixed set of
 * grid coordinates, node_x.
 */

#define SHAPEFUN_PRIVATE_DATA \
   int method_id;
#include "shapefun.h"

enum
{
   SHAPEFUN_TYPE_CSPLINE,
   SHAPEFUN_TYPE_CHEB,
   SHAPEFUN_TYPE_POLY,
   SHAPEFUN_TYPE_POSITIVE_COEFF
};

static void free_shapefun_type (Shapefun_Type *st)
{
   if (st == NULL)
     return;
   FREE(st);
}

static int st_num_params (const Shapefun_Type *st)
{
   int num_params;

   switch (st->method_id)
     {
      case SHAPEFUN_TYPE_CSPLINE:
        num_params = st->num_nodes;
        break;

      case SHAPEFUN_TYPE_CHEB:
      case SHAPEFUN_TYPE_POLY:
        num_params = st->num_coef;
        break;

      case SHAPEFUN_TYPE_POSITIVE_COEFF:
        num_params = 1;
        break;

      default:
        num_params = 0;
        break;
     }

   return num_params;
}

static int cspline_config (const Shapefun_Type *tf, Cspline_Config_Type *c)
{
   (void) tf; (void) c;
   return 0;
}

static int cheb_config (const Shapefun_Type *st, Cheb_Config_Type *c)
{
   if (st == NULL || c == NULL)
     return -1;
   c->xmin = st->xmin;
   c->xmax = st->xmax;
   c->num_coef = st->num_coef;
   return 0;
}

static int poly_config (const Shapefun_Type *st, Poly_Config_Type *c)
{
   if (st == NULL || c == NULL)
     return -1;
   c->num_coef = st->num_coef;
   c->xexp = st->xexp;
   return 0;
}

static int cspline_onto_control_nodes (const Shapefun_Type *st,
                                       const Shapefun_Init_Type *init,
                                       double *node_coeffs)
{
   Cspline_Config_Type cct = {0};
   Cspline_Type *ct = NULL;
   int status;
   if (cspline_config (st, &cct)<0)
     return -1;
   if (NULL == (ct = cspline_interpol (init->n, init->x, init->y, &cct)))
     return -1;
   status = cspline_eval (ct, st->num_nodes, st->node_x, node_coeffs);
   cspline_free (ct);
   return status;
}

static int
init_cspline_shapefun_params (const Shapefun_Type *st,
                              const Shapefun_Init_Type *init,
                              size_t num_nodes, double *node_coeffs)
{
   size_t n = st->num_nodes;

   if ((init == NULL) || (init->x == NULL) || (init->y == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: got NULL pointer", __func__);
        return -1;
     }

   if (n == 0)
     return 0;

   /* For spline scale functions, the node_coeffs[i] values are interpreted
    * as scale function values at each node_x[i] coordinate */

   if (num_nodes < n)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: array size too small (size=%ld, %ld is required)",
                     __func__, num_nodes, n);
        return -1;
     }

   if (init->n < 2)
     {
        double s0 = (init->n == 0) ? 0.0 : init->y[0];
        size_t i;

        for (i = 0; i < n; i++)
          {
             node_coeffs[i] = s0;
          }
     }
   else
     {
        return cspline_onto_control_nodes(st, init, node_coeffs);
     }

   return 0;
}

static int eval_cspline_shapefun (const Shapefun_Type *st, const double *node_coeffs,
                                  size_t n, const double *x, double *y)
{
   Cspline_Config_Type cct = {0};
   Cspline_Type *ct = NULL;
   int status;

   if (cspline_config (st, &cct) < 0)
     return -1;

   /* For spline scale functions, the node_coeffs[i] values are interpreted
    * as scale function values at each node_x[i] coordinate */

   if (NULL == (ct = cspline_interpol (st->num_nodes, st->node_x, node_coeffs, &cct)))
     return -1;

   status = cspline_eval (ct, n, x, y);
   cspline_free (ct);

   if (status)
     return -1;

   return 0;
}

static int init_cheb_shapefun_params (const Shapefun_Type *st,
                                      const Shapefun_Init_Type *init,
                                      size_t num_nodes, double *node_coeffs)
{
   size_t n = st->num_coef;

   /* For Chebyshev scale functions, the node_coeffs[i] values are interpreted
    * as the Chebyshev series coefficients (contrast with usage for spline
    * scale functions) */

   if (num_nodes != n)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size (size=%ld, %ld is required)",
                     __func__, num_nodes, n);
        return -1;
     }

   if (n == 0)
     return 0;

   if ((init == NULL) || (init->x == NULL) || (init->y == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: got NULL pointer", __func__);
        return -1;
     }

   if (init->n < 2)
     {
        double s0 = (init->n == 0) ? 0.0 : init->y[0];
        size_t i;

        node_coeffs[0] = s0;
        for (i = 1; i < n; i++) node_coeffs[i] = 0.0;
     }
   else
     {
        Cheb_Config_Type cct = {0};
        Cheb_Type *ct = NULL;
        int status;

        if (cheb_config (st, &cct) < 0)
          return -1;

        if (NULL == (ct = cheb_interpol (init->n, init->x, init->y, &cct)))
          return -1;
        status = cheb_get_coeffs (ct, n, node_coeffs);
        cheb_free (ct);
        if (status != 0)
          return -1;
     }

   return 0;
}

static int eval_cheb_shapefun (const Shapefun_Type *st, const double *node_coeffs,
                               size_t n, const double *x, double *y)
{
   double a = st->xmin;
   double b = st->xmax;
   size_t i, num_coef = st->num_coef;

   for (i = 0; i < n; i++)
     {
        if (0 != cheb_clenshaw_eval (a, b, node_coeffs, num_coef, x[i], &y[i]))
          return -1;
     }

   return 0;
}

static int init_poly_shapefun_params (const Shapefun_Type *st,
                                      const Shapefun_Init_Type *init,
                                      size_t num_nodes,
                                      double *node_coeffs)
{
   size_t n = st->num_coef;

   /* For polynomial scale functions, the node_coeffs[i] values are interpreted
    * as the polynomial coefficients (contrast with usage for spline
    * scale functions) */

   if (num_nodes != n)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size (size=%ld, %ld is required)",
                     __func__, num_nodes, n);
        return -1;
     }

   if (n == 0)
     return 0;

   if ((init == NULL) || (init->x == NULL) || (init->y == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: got NULL pointer", __func__);
        return -1;
     }

   if (init->n < 2)
     {
        double s0 = (init->n == 0) ? 0.0 : init->y[0];
        size_t i;

        node_coeffs[0] = s0;
        for (i = 1; i < n; i++) node_coeffs[i] = 0.0;
     }
   else
     {
        Poly_Config_Type pct = {0};
        Poly_Type *pt = NULL;
        int status;

        if (poly_config (st, &pct) < 0)
          return -1;

        if (NULL == (pt = poly_interpol (init->n, init->x, init->y, &pct)))
          return -1;

        status = poly_get_coeffs (pt, n, node_coeffs);
        poly_free (pt);

        if (0 != status)
          return -1;
     }

   return 0;
}

static int eval_poly_shapefun (const Shapefun_Type *st, const double *node_coeffs,
                               size_t n, const double *x, double *y)
{
   double x0 = st->xexp;
   int len = st->num_coef;
   size_t i;

   for (i = 0; i < n; i++)
     {
        y[i] = gsl_poly_eval (node_coeffs, len, x[i] - x0);
     }

   return 0;
}

static int
init_positive_coeff_shapefun_params (const Shapefun_Type *st,
                                     const Shapefun_Init_Type *init,
                                     size_t num_nodes,
                                     double *node_coeffs)
{
   (void) st;

   /* For 'positive_coeff' scale functions, node_coeffs[0]
    * is interpreted as the parameter to be squared */

   if (num_nodes != 1)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: incorrect destination array size (size=%ld, 1 is required)",
                     __func__, num_nodes);
        return -1;
     }

   if ((init == NULL) || (init->y == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: got NULL pointer", __func__);
        return -1;
     }

   node_coeffs[0] = init->y[0];

   return 0;
}

static int
eval_positive_coeff_shapefun (const Shapefun_Type *st, const double *node_coeffs,
                              size_t n, const double *x, double *y)
{
   double s2, s = node_coeffs[0];
   size_t i;

   (void) st; (void) x;

   s2 = s*s;

   for (i = 0; i < n; i++)
     {
        y[i] = s2;
     }

   return 0;
}

typedef struct
{
   const char *name;
   int method_id;

   int (*st_init_params)(const Shapefun_Type *,
                         const Shapefun_Init_Type *,
                         size_t, double *);
   int (*st_eval)(const Shapefun_Type *, const double *,
                  size_t, const double *, double *);
}
Shapefun_Method;
#define TF_METHOD(n,typ) \
   {#n,typ, \
        init_##n##_shapefun_params, \
        eval_##n##_shapefun}
#define TF_METHOD_TABLE_END {NULL,0,NULL,NULL}

static Shapefun_Method Tf_Method_Table[] =
{
   TF_METHOD(cspline,SHAPEFUN_TYPE_CSPLINE),
   TF_METHOD(cheb,SHAPEFUN_TYPE_CHEB),
   TF_METHOD(poly,SHAPEFUN_TYPE_POLY),
   TF_METHOD(positive_coeff,SHAPEFUN_TYPE_POSITIVE_COEFF),
   TF_METHOD_TABLE_END
};

static Shapefun_Method *find_shapefun_method (const char *name)
{
   Shapefun_Method *m = NULL;

   for (m = Tf_Method_Table; m->name != NULL; m++)
     {
        if (0 == strcmp (name, m->name))
          return m;
     }

   return NULL;
}

static int shapefun_config (Shapefun_Type *st, const char *method_name)
{
   Shapefun_Method *m = NULL;

   if (NULL == (m = find_shapefun_method (method_name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: unsupported method = %s",
                     __func__, method_name);
        return -1;
     }

   st->method_id = m->method_id;
   st->st_delete = free_shapefun_type;
   st->st_init_params = m->st_init_params;
   st->st_eval = m->st_eval;
   st->st_num_params = st_num_params;

   return 0;
}

static Shapefun_Type *alloc_shapefun_type (void)
{
   Shapefun_Type *st = NULL;

   if (NULL == (st = (Shapefun_Type *)MALLOC (sizeof *st)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)st, 0, sizeof *st);

   return st;
}

Shapefun_Type *shapefun_create (const char *method_name)
{
   Shapefun_Type *st = NULL;

   if (NULL == (st = alloc_shapefun_type ()))
     return NULL;

   if (0 != shapefun_config (st, method_name))
     {
        free_shapefun_type (st);
        return NULL;
     }

   return st;
}
