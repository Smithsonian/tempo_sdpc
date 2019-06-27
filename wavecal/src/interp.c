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
#include "util.h"

#define INTERP_TYPE_PRIVATE_DATA \
   int method_id; \
   void *method;
#include "interp.h"

static void free_interp_type (Interp_Type *it)
{
   if (it == NULL)
     return;
   FREE(it);
}

static Interp_Type *alloc_interp_type (void)
{
   Interp_Type *it = NULL;

   if (NULL == (it = (Interp_Type *)MALLOC (sizeof *it)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)it, 0, sizeof *it);

   return it;
}

static int it_method_id (const Interp_Type *it)
{
   return it->method_id;
}

static int cspline_config (Interp_Type *it, Cspline_Config_Type *c)
{
   (void) it; (void) c;
   return 0;
}

static int cheb_config (Interp_Type *it, Cheb_Config_Type *c)
{
   if (c == NULL || it == NULL)
     return -1;
   c->xmin = it->xmin;
   c->xmax = it->xmax;
   c->num_coef = it->num_coef;
   return 0;
}

static int poly_config (Interp_Type *it, Poly_Config_Type *c)
{
   if (c == NULL || it == NULL)
     return -1;
   c->num_coef = it->num_coef;
   c->xexp = it->xexp;
   return 0;
}

#define INTERP_DELETE_METHOD(typ_name,sym_name) \
static void delete_##sym_name (Interp_Type *it) \
{ \
   if (it == NULL) \
     return; \
   sym_name##_free ((typ_name##_Type *)it->method); \
   free_interp_type (it); \
}

#define INTERP_EVAL_METHOD(typ_name,sym_name) \
static int eval_##sym_name (Interp_Type *it, size_t n, const double *x, double *y) \
{ \
   typ_name##_Type *mt = (typ_name##_Type *)it->method; \
   if (mt == NULL) \
     { \
        tell_verror (TELL_RUNTIME_ERROR, "%s: method %s not initialized", __func__, #sym_name); \
     } \
   return sym_name##_eval (mt, n, x, y); \
}

#define INTERP_INIT_METHOD(typ_name,sym_name) \
static int init_##sym_name (Interp_Type *it, \
                            size_t n, const double *x, const double *y) \
{ \
   typ_name##_Config_Type m_cfg = {0}; \
 \
   sym_name##_free ((typ_name##_Type *)it->method); \
   if (sym_name##_config (it, &m_cfg) < 0) \
     return -1; \
 \
   if (NULL == (it->method = sym_name##_interpol (n, x, y, &m_cfg))) \
     return -1; \
 \
   return 0; \
}

#define INTERP_METHOD_SRC(typ_name,sym_name) \
   INTERP_DELETE_METHOD(typ_name,sym_name) \
   INTERP_INIT_METHOD(typ_name,sym_name) \
   INTERP_EVAL_METHOD(typ_name,sym_name)

INTERP_METHOD_SRC(Cspline,cspline)
INTERP_METHOD_SRC(Cheb,cheb)
INTERP_METHOD_SRC(Poly,poly)

typedef struct
{
   const char *name;
   int method_id;

   void (*it_delete)(Interp_Type *);
   int (*it_interp_init)(Interp_Type *, size_t, const double *, const double *);
   int (*it_interp_eval)(Interp_Type *, size_t, const double *, double *);
}
Interp_Method;
#define INTERP_METHOD(n,typ) \
   {#n,typ, \
        delete_##n, \
        init_##n, \
        eval_##n}
#define INTERP_METHOD_TABLE_END {NULL,0,NULL,NULL,NULL}

static Interp_Method Interp_Method_Table[] =
{
   INTERP_METHOD(cspline,INTERP_TYPE_CSPLINE),
   INTERP_METHOD(cheb,INTERP_TYPE_CHEB),
   INTERP_METHOD(poly,INTERP_TYPE_POLY),
   INTERP_METHOD_TABLE_END
};

static Interp_Method *find_interp_method (const char *name)
{
   Interp_Method *m = NULL;

   for (m = Interp_Method_Table; m->name != NULL; m++)
     {
        if (0 == strcmp (name, m->name))
          return m;
     }

   return NULL;
}

static int interp_config (Interp_Type *it, const char *method_name)
{
   Interp_Method *m = NULL;

   if (NULL == (m = find_interp_method (method_name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: unsupported method = %s",
                     __func__, method_name);
        return -1;
     }

   it->method_id = m->method_id;
   it->it_method_id = it_method_id;

   it->it_delete = m->it_delete;
   it->it_interp_init = m->it_interp_init;
   it->it_interp_eval = m->it_interp_eval;

   return 0;
}

Interp_Type *
interp_create (const char *method_name)
{
   Interp_Type *it = NULL;

   if (NULL == (it = alloc_interp_type ()))
     return NULL;

   if (0 != interp_config (it, method_name))
     {
        free_interp_type (it);
        return NULL;
     }

   return it;
}
