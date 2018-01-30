/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <tell.h>
#include <tio.h>
#include <gsl/gsl_errno.h>

#include "config.h"
#include "util.h"
#include "interp_cspline.h"
#include "interp.h"
#include "shapefun.h"
#include "wavecal.h"
#include "mpfit.h"

#define NCID_UNINITIALIZED -1

enum
{
   TERM_TYPE_AD1,
   TERM_TYPE_AD2,
   TERM_TYPE_LBE,
   TERM_TYPE_SC,
   TERM_TYPE_BL
};

#define NUM_TERMS 5

typedef struct
{
   int xtrack;         /**< cross-track index to be fitted */
   int num_wave;       /**< total number of wavelength points in measured spectrum */
   double *wave0;      /**< initial guess at wavelength grid for measured spectrum */
   double wave_pad;
   int index_lim[2];   /**< index limits on the wave0[] grid to be fitted */

   Shapefun_Type *shapefun;
   double *pindex;
   double *wave_params;
   size_t num_wave_params;
}
Window_Type;

typedef struct
{
   const char *path;
   const char *name_x;
   const char *name_y;
   double scale_factor;
}
File_Type;

typedef struct
{
   File_Type file;
   Interp_Type *interp;
}
Refspec_Type;

typedef struct Term_Type Term_Type;
struct Term_Type
{
   Term_Type *next;
   char *term_name;
   int term_type;
   Refspec_Type refspec;
   Shapefun_Type *shapefun;
   Shapefun_Init_Type shapefun_init;
   double *params;
   size_t num_params;
   double *value_workspace;
   size_t num_values;
};

typedef struct
{
   File_Type file;
   int ncid;
   Cspline_Type *cspline;
   double *irr0_workspace;
   double *irradiance;
   double *wavelen;
   size_t num_wavelen;
}
Reference_Irr_Type;

typedef struct
{
   double xtol;
   double ftol;
   int maxiter;
   int maxfev;
}
Fit_Control_Type;

struct Wavecal_Type
{
   Reference_Irr_Type irr;
   Window_Type window;
   Term_Type *terms;
   Fit_Control_Type fit_ctrl;
   int index_lim_ref[2];
   int mode;
};

static void free_interp_type (Interp_Type *it)
{
   if (it == NULL)
     return;
   if (it->it_delete)
     {
        it->it_delete (it);
     }
}

static void free_shapefun_type (Shapefun_Type *st)
{
   if (st == NULL)
     return;
   if (st->malloced_node_x)
     {
        FREE(st->node_x);
     }
   if (st->st_delete)
     {
        st->st_delete (st);
     }
}

static void free_shapefun_init_type (Shapefun_Init_Type *it)
{
   if (it == NULL)
     return;
   if (it->malloced)
     {
        FREE(it->x);
        FREE(it->y);
     }
}

static void free_term1 (Term_Type *tt)
{
   if (tt == NULL)
     return;
   FREE(tt->term_name);
   FREE(tt->params);
   FREE(tt->value_workspace);
   free_interp_type (tt->refspec.interp);
   free_shapefun_type (tt->shapefun);
   free_shapefun_init_type (&tt->shapefun_init);
   FREE(tt);
}

static void free_terms (Term_Type *tt)
{
   while (tt)
     {
        Term_Type *next = tt->next;
        free_term1 (tt);
        tt = next;
     }
}

static Term_Type *term_alloc (void)
{
   Term_Type *tt = NULL;

   if (NULL == (tt = (Term_Type *)MALLOC (sizeof *tt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)tt, 0, sizeof *tt);

   return tt;
}

static int read_filepar (config_setting_t *s, File_Type *ct)
{
   if ((CONFIG_TRUE != config_setting_lookup_string (s, "path", &ct->path))
       ||(CONFIG_TRUE != config_setting_lookup_string (s, "wavegrid", &ct->name_x))
       ||(CONFIG_TRUE != config_setting_lookup_string (s, "value", &ct->name_y)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading config setting 'file' (%s:%d)",
                     __func__, config_setting_source_file (s),
                     config_setting_source_line (s));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "scale_factor", &ct->scale_factor))
     ct->scale_factor = 1.0;

   return 0;
}

static void tell_config_error (config_setting_t *s, const char *func)
{
   tell_verror (TELL_INVALID_PARM_ERROR, "%s: reading from %s (%s:%d)",
                func, config_setting_name(s),
                config_setting_source_file(s),
                config_setting_source_line(s));
}

static int config_interp_method (config_setting_t *s, Interp_Type *it)
{
   int num_coef;

   /* FIXME - should implement method-specific validation */

   (void) config_setting_lookup_float (s, "xmin", &it->xmin);
   (void) config_setting_lookup_float (s, "xmax", &it->xmax);
   (void) config_setting_lookup_float (s, "xexp", &it->xexp);
   (void) config_setting_lookup_int (s, "num_coef", &num_coef);
   it->num_coef = num_coef;  /* note type conversion */

   return 0;
}

static int config_shapefun_method (config_setting_t *s,
                                   Shapefun_Type *st,
                                   Shapefun_Init_Type *st_init)
{
   config_setting_t *ss;
   int num_coef;

   /* FIXME - should implement method-specific validation */

   (void) config_setting_lookup_float (s, "xmin", &st->xmin);
   (void) config_setting_lookup_float (s, "xmax", &st->xmax);
   (void) config_setting_lookup_float (s, "xexp", &st->xexp);
   (void) config_setting_lookup_int (s, "num_coef", &num_coef);
   st->num_coef = num_coef;  /* note type conversion */

   if (NULL == (ss = config_setting_get_member (s, "function")))
     return 0;

   if (read_config_float_array (ss, "nodes", &st->node_x, &st->num_nodes) < 0)
     return -1;
   st->malloced_node_x = 1;

   if (st_init)
     {
        size_t ny;
        if (read_config_float_array (ss, "x0", &st_init->x, &st_init->n) < 0)
          return -1;
        st_init->malloced = 1;
        if (read_config_float_array (ss, "y0", &st_init->y, &ny) < 0)
          return -1;

        if (ny != st_init->n)
          {
             tell_config_error (ss, __func__);
          }
     }

   return 0;
}

static int open_refspec (config_setting_t *s,
                         File_Type *file, Interp_Type **pinterp)
{
   const char *method_name;
   config_setting_t *ss;
   Interp_Type *interp = NULL;

   *pinterp = NULL;

   /* optional */
   if (NULL == (ss = config_setting_get_member (s, "data")))
     return 0;

   if (0 != read_filepar (ss, file))
     return -1;

   if (CONFIG_FALSE == config_setting_lookup_string (s, "method", &method_name))
     {
        tell_config_error (s, __func__);
        return -1;
     }

   if (NULL == (interp = interp_create (method_name)))
     return -1;

   if ((interp != NULL)
       && (0 != config_interp_method (s, interp)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: failed configuring interp method: %s",
                     __func__, method_name);
        interp->it_delete (interp);
        return -1;
     }

   *pinterp = interp;

   return 0;
}

static int open_shapefun (config_setting_t *s,
                          Shapefun_Init_Type *shapefun_init,
                          Shapefun_Type **pshapefun)
{
   Shapefun_Type *shapefun = NULL;
   const char *method_name;

   *pshapefun = NULL;

   if (CONFIG_FALSE == config_setting_lookup_string (s, "method", &method_name))
     {
        tell_config_error (s, __func__);
        return -1;
     }

   if (NULL == (shapefun = shapefun_create (method_name)))
     return -1;

   if ((shapefun != NULL)
       && (0 != config_shapefun_method (s, shapefun, shapefun_init)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: failed configuring shapefun method: %s",
                     __func__, method_name);
        shapefun->st_delete (shapefun);
        return -1;
     }

   *pshapefun = shapefun;

   return 0;
}

static int get_term_type_id (const char *term_type)
{
   typedef struct
     {
        const char *name;
        int id;
     }
   symbol_type;
   static symbol_type symbol_table[] =
     {
        {"ad1", TERM_TYPE_AD1},
        {"ad2", TERM_TYPE_AD2},
        {"lbe", TERM_TYPE_LBE},
        {"sc",  TERM_TYPE_SC},
        {"bl",  TERM_TYPE_BL},
        {NULL, -1}
     };
   symbol_type *s;

   for (s = symbol_table; s->name != NULL; s++)
     {
        if (0 == strcasecmp (s->name, term_type))
          return s->id;
     }

   tell_verror (TELL_INVALID_PARM_ERROR, "%s: unknown term type: %s",
                __func__, term_type);

   return -1;
}

static Term_Type *term_open (config_setting_t *s)
{
   Term_Type *tt = NULL;
   const char *term_type;

   if (NULL == (tt = term_alloc ()))
     return NULL;

   if (NULL == (tt->term_name = strdup (config_setting_name(s))))
     goto return_error;

   if (CONFIG_FALSE == config_setting_lookup_string (s, "term", &term_type))
     goto return_error;

   if ((tt->term_type = get_term_type_id (term_type)) < 0)
     goto return_error;

   if (open_refspec (s, &tt->refspec.file, &tt->refspec.interp) < 0)
     goto return_error;

   if (open_shapefun (s, &tt->shapefun_init, &tt->shapefun) < 0)
     goto return_error;

   if ((tt->refspec.interp == NULL) && (tt->shapefun == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: NULL configuration encountered in %s:%s:%d",
                     __func__, config_setting_name(s),
                     config_setting_source_file(s),
                     config_setting_source_line(s));
        goto return_error;
     }

   return tt;
return_error:
   free_term1 (tt);
   return NULL;
}

static void free_window (Window_Type *win)
{
   if (win == NULL)
     return;
   free_shapefun_type (win->shapefun);
   FREE(win->wave_params);
   FREE(win->wave0);
   FREE(win->pindex);
}

static void free_reference_irr_type (Reference_Irr_Type *irr)
{
   if (irr == NULL)
     return;
   TIO_close (irr->ncid);
   cspline_free (irr->cspline);
   FREE(irr->wavelen);
   FREE(irr->irradiance);
   FREE(irr->irr0_workspace);
}

void wavecal_close (Wavecal_Type *wct)
{
   if (wct == NULL)
     return;
   free_reference_irr_type (&wct->irr);
   free_terms (wct->terms);
   free_window (&wct->window);
   FREE(wct);
}

static Wavecal_Type *alloc_wavecal (void)
{
   Wavecal_Type *wct = NULL;

   if (NULL == (wct = (Wavecal_Type *)MALLOC (sizeof *wct)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)wct, 0, sizeof *wct);

   return wct;
}

static int config_model_components (Wavecal_Type *wct, config_setting_t *s)
{
   Term_Type *term;
   config_setting_t *ss;
   unsigned int i, num;

   num = config_setting_length (s);

   for (i = 0; i < num; i++)
     {
        const char *name;
        if (NULL == (ss = config_setting_get_elem (s, i)))
          return -1;

        name = config_setting_name (ss);
        if (name[0] == '*') continue;
        if (1) fprintf (stderr, "add component: %s\n", name);

        if (NULL == (term = term_open (ss)))
          return -1;
        term->next = wct->terms;
        wct->terms = term;
     }

   return 0;
}

static int config_irr_reference (config_setting_t *s, Reference_Irr_Type *irr)
{
   File_Type *file = &irr->file;
   config_setting_t *ss;

   if (NULL == (ss = config_setting_get_member (s, "data")))
     return -1;
   if (0 != read_filepar (ss, file))
     return -1;

   if (0 != TIO_open (file->path, NC_NOWRITE, &irr->ncid))
     return -1;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file->path);

   return 0;
}

static int config_fit_window (config_setting_t *s, Window_Type *win)
{
   Shapefun_Type *st;
   const char *method_name;

   if (CONFIG_FALSE == config_setting_lookup_string (s, "method", &method_name))
     {
        tell_config_error (s, __func__);
        return -1;
     }

   if (NULL == (win->shapefun = shapefun_create (method_name)))
     return -1;

   if (0 != config_shapefun_method (s, win->shapefun, NULL))
     return -1;

   st = win->shapefun;
   win->num_wave_params = st->st_num_params (st);
   if (win->num_wave_params <= 0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: irradiance wavelength grid has no variable parameters",
                     __func__);
        return -1;
     }
   if (NULL == (win->wave_params = alloc_doubles (win->num_wave_params)))
     return -1;

   return 0;
}

static int config_control (config_setting_t *s, Wavecal_Type *wct)
{
   Fit_Control_Type *c = &wct->fit_ctrl;

   /* set defaults */
   c->xtol = 1.e-8;
   c->ftol = 1.e-8;
   c->maxiter = 200;
   c->maxfev = 0;

   /* defaults are over-ridden only if the setting is present
    * and provides the correct data type */
   (void) config_setting_lookup_float (s, "xtol", &c->xtol);
   (void) config_setting_lookup_float (s, "ftol", &c->ftol);
   (void) config_setting_lookup_int (s, "maxiter", &c->maxiter);
   (void) config_setting_lookup_int (s, "maxfev", &c->maxfev);

   return 0;
}

Wavecal_Type *wavecal_open (config_t *cfg, int mode)
{
   Wavecal_Type *wct = NULL;
   config_setting_t *s;

   if (NULL == (wct = alloc_wavecal ()))
     goto error_return;

   wct->mode = mode;

   if (NULL == (s = config_lookup (cfg, "wavecal_control")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing wavecal_control in param file: %s",
                     __func__, config_error_file (cfg));
        goto error_return;
     }

   if (0 != config_control (s, wct))
     goto error_return;

   if (NULL == (s = config_lookup (cfg, "wavecal_irr_reference")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing wavecal_irr_reference in param file: %s",
                     __func__, config_error_file (cfg));
        goto error_return;
     }

   if (0 != config_fit_window (s, &wct->window))
     goto error_return;

   if (0 != config_irr_reference (s, &wct->irr))
     goto error_return;

   if (mode == 1)
     {
        if (NULL == (s = config_lookup (cfg, "wavecal_rad_model")))
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: accessing wavecal_rad_model in param file: %s",
                          __func__, config_error_file (cfg));
             goto error_return;
          }

        if (0 != config_model_components (wct, s))
          goto error_return;
     }

   return wct;

error_return:
   wavecal_close(wct);
   return NULL;
}

static int read_var (int ncid, const char *name, int xtrack,
                     double **px, size_t *pnum)
{
   TIO_Var_Info_Type info = {0};
   int start[TIO_MAX_VAR_DIMS] = {0};
   int count[TIO_MAX_VAR_DIMS] = {0};
   double *x = NULL;
   size_t n;
   int status = -1;

   *px = NULL;
   if (pnum) *pnum = 0;

   if (0 != TIO_inq_var (ncid, name, &info))
     return -1;

   switch (info.ndims)
     {
      case 1:
        n = info.dimlens[0];
        start[0] = 0;
        count[0] = n;
        break;
      case 2:
        n = info.dimlens[1];
        start[0] = xtrack;
        start[1] = 0;
        count[0] = 1;
        count[1] = n;
        break;
      default:
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unexpected variable dimension (num_dims=%d) (%s)",
                     __func__, info.ndims, name);
        goto error_return;
     }

   if (NULL == (x = alloc_doubles (n)))
     goto error_return;

   if (0 != TIO_get_var_section (ncid, name, start, count, TIO_DOUBLE, x))
     goto error_return;

   *px = x;
   if (pnum) *pnum = n;

   status = 0;
error_return:
   if (status) FREE(x);
   return status;
}

static int read_irr_reference (Reference_Irr_Type *irr, int xtrack)
{
   File_Type *file = &irr->file;

   xtrack = 0;
   fprintf (stderr, "FIXME: forcing xtrack=%d on irradiance input\n", xtrack);

   /* FIXME better to allocate a buffer once, and re-use it */
   FREE(irr->wavelen);
   FREE(irr->irradiance);

   /* During operations, the reference spectrum file is expected
    * to be fairly large (e.g. several gigabytes) so rather than
    * read the whole file at once, we'll read reference spectra
    * as needed.  If I/O is a problem, we can consider a caching
    * layer. (Does netcdf already provide caching on input?)
    */
   if (0 != read_var (irr->ncid, file->name_x, xtrack,
                      &irr->wavelen, &irr->num_wavelen))
     return -1;

   if (0 != read_var (irr->ncid, file->name_y, xtrack,
                      &irr->irradiance, NULL))
     return -1;

   if (file->scale_factor != 1.0)
     {
        double *irr_i = irr->irradiance;
        double scale = file->scale_factor;
        size_t i;
        for (i = 0; i < irr->num_wavelen; i++)
          {
             irr_i[i] *= scale;
          }
     }

   FREE(irr->irr0_workspace); /* FIXME - preallocate this */
   if (NULL == (irr->irr0_workspace = alloc_doubles (irr->num_wavelen)))
     return -1;

   return 0;
}

static int init_shapefun (Term_Type *term)
{
   Shapefun_Type *st = term->shapefun;
   Shapefun_Init_Type *st_init = &term->shapefun_init;

   if (st == NULL)
     return 0;

   term->num_params = st->st_num_params (st);
   if (term->num_params > 0)
     {
        FREE(term->params);  /* FIXME - preallocate this */
        if (NULL == (term->params = alloc_doubles (term->num_params)))
          return -1;
        if (0 != st->st_init_params (st, st_init, term->num_params, term->params))
          return -1;
     }
   else term->params = NULL;

   return 0;
}

static int interp_filepar_init (Interp_Type *interp,
                                const File_Type *file, int xtrack)
{
   double *x = NULL;
   double *y = NULL;
   size_t nx, ny, n;
   int ncid, status = -1;

   if (0 != TIO_open (file->path, NC_NOWRITE, &ncid))
     return -1;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file->path);

   if ((0 != read_var (ncid, file->name_x, xtrack, &x, &nx))
       ||(0 != read_var (ncid, file->name_y, xtrack, &y, &ny)))
     goto error_return;

   if (nx == ny)
     {
        n = nx;
     }
   else
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: grid mismatch (nx=%ld, ny=%ld)",
                     __func__, nx, ny);
        goto error_return;
     }

   if (file->scale_factor != 1.0)
     {
        double scale_factor = file->scale_factor;
        size_t i;
        for (i = 0; i < n; i++)
          {
             y[i] *= scale_factor;
          }
     }

   if (0 != interp->it_interp_init (interp, n, x, y))
     goto error_return;

   status = 0;
error_return:
   TIO_close(ncid);
   FREE(x);
   FREE(y);
   return status;
}

static int read_rad_reference (Term_Type *terms, int xtrack)
{
   Term_Type *term;

   for (term = terms; term != NULL; term = term->next)
     {
        if (init_shapefun (term) < 0)
          return -1;

        if (term->refspec.interp == NULL)
          continue;

        if (0 != interp_filepar_init (term->refspec.interp,
                                      &term->refspec.file, xtrack))
          return -1;
     }

   return 0;
}

static int init_window (Wavecal_Type *wct, int xtrack,
                        int num_wave, const double *wave,
                        const Wavecal_Config_Type *config)
{
   Reference_Irr_Type *irr = &wct->irr;
   Window_Type *win = &wct->window;
   Shapefun_Type *shapefun = win->shapefun;
   Shapefun_Init_Type shapefun_init = {0};
   double beg_wave, end_wave;
   double *wave_start;
   double *irr_start;
   int i, num;

   if (0 != read_irr_reference (&wct->irr, xtrack))
     return -1;

   if (0 != read_rad_reference (wct->terms, xtrack))
     return -1;

   win->xtrack = xtrack;
   win->wave_pad = config->wave_pad;
   win->num_wave = num_wave;

   if (config->index_lim)
     {
        win->index_lim[0] = config->index_lim[0];
        win->index_lim[1] = config->index_lim[1];
     }
   else
     {
        win->index_lim[0] = 0;
        win->index_lim[1] = win->num_wave-1;
     }

   FREE(win->wave0); /* FIXME - preallocate this */
   if (NULL == (win->wave0 = alloc_doubles (win->num_wave)))
     return -1;
   memcpy ((char *)win->wave0, (char *)wave, win->num_wave * sizeof(double));

   FREE(win->pindex); /* FIXME - preallocate this */
   if (NULL == (win->pindex = alloc_doubles (win->num_wave)))
     return -1;
   for (i = 0; i < win->num_wave; i++)
     {
        win->pindex[i] = (double)i;
     }

   shapefun_init.x = win->pindex;
   shapefun_init.y = win->wave0;
   shapefun_init.n = win->num_wave;

   shapefun->xmin = win->index_lim[0];
   shapefun->xmax = win->index_lim[1];

   if (0 != shapefun->st_init_params (shapefun, &shapefun_init,
                                      win->num_wave_params, win->wave_params))
     return -1;

   /* derive padded wavelength range for reference spectra */
   beg_wave = win->wave0[ win->index_lim[0] ] - win->wave_pad;
   end_wave = win->wave0[ win->index_lim[1] ] + win->wave_pad;

   if ((wct->index_lim_ref[0] = find_x (beg_wave, irr->wavelen, irr->num_wavelen)) < 0)
     return -1;
   if ((wct->index_lim_ref[1] = find_x (end_wave, irr->wavelen, irr->num_wavelen)) < 0)
     return -1;

   num = wct->index_lim_ref[1] - wct->index_lim_ref[0] + 1;
   wave_start = irr->wavelen + wct->index_lim_ref[0];
   irr_start = irr->irradiance + wct->index_lim_ref[0];

   cspline_free (irr->cspline); /* FIXME - preallocate this */
   if (NULL == (irr->cspline = cspline_interpol (num, wave_start, irr_start, NULL)))
     return -1;

   return 0;
}

static int collect_params (Wavecal_Type *wct, size_t *pnum, double **pparams)
{
   Window_Type *win = &wct->window;
   Term_Type *term;
   double *params = NULL;
   double *pnext;
   int num;

   /* count parameters */
   num = win->num_wave_params;
   for (term = wct->terms; term != NULL; term = term->next)
     {
        if (term->shapefun)
          {
             num += term->num_params;
          }
     }

   /* allocate storage and copy initial values */
   if (NULL == (params = alloc_doubles (num)))
     return -1;

   memcpy ((char *)params, (char *)win->wave_params,
           win->num_wave_params * sizeof(double));
   pnext = params + win->num_wave_params;

   for (term = wct->terms; term != NULL; term = term->next)
     {
        if (term->shapefun)
          {
             memcpy ((char *)pnext, (char *)term->params,
                     term->num_params * sizeof(double));
             pnext += term->num_params;
          }
     }

   *pnum = num;
   *pparams = params;

   return 0;
}

static int evaluate_term (Term_Type *term, Window_Type *win,
                          const double *params)
{
   Refspec_Type *ref = &term->refspec;
   size_t i, offset, n = win->num_wave;
   double *v = NULL;

   if (term->value_workspace == NULL)
     {
        /* we may have 2 factors, so allocate space for both */
        if (NULL == (v = alloc_doubles (2 * n)))
          return -1;
        term->value_workspace = v;
        term->num_values = n;
     }

   v = term->value_workspace;
   memset ((char *)v, 0, 2*n * sizeof(double));

   offset = 0;

   if (term->shapefun)
     {
        Shapefun_Type *st = term->shapefun;
        offset = n;
        if (st->st_eval (st, params, n, win->wave0, v))
          return -1;
        if (0) fprintf (stderr, "evaluated term %s (shapefun)\n", term->term_name);
     }

   if (ref->interp)
     {
        Interp_Type *it = ref->interp;
        if (it->it_interp_eval (it, n, win->wave0, v + offset) < 0)
          return -1;
        if (offset)
          {
             for (i = 0; i < n; i++) v[i] *= v[i+offset];
          }
     }

   return 0;
}

typedef struct
{
   double *v;
   size_t n;
}
Array_Type;

static void array_free (Array_Type *a)
{
   FREE(a->v);
   a->v = NULL;
   a->n = 0;
}

static int array_alloc (Array_Type *a, size_t n)
{
   if (NULL == (a->v = (double *)MALLOC (n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)a->v, 0, n * sizeof(double));
   a->n = n;
   return 0;
}

static void free_term_arrays (Array_Type *a, size_t na)
{
   size_t i;
   for (i = 0; i < na; i++)
     {
        array_free (&a[i]);
     }
}
static int alloc_term_arrays (Array_Type *a, size_t na, size_t len)
{
   size_t i;
   for (i = 0; i < na; i++)
     {
        if (0 != array_alloc (&a[i], len))
          {
             free_term_arrays (a, na);
             return -1;
          }
     }
   return 0;
}

static int Write_Output = 0;
static int write_term (const Term_Type *t, const Window_Type *win)
{
#define BUFSIZE 256
   char filename[BUFSIZE];
   FILE *fp;
   size_t i, n;

   if (Write_Output == 0)
     return 0;

   if (snprintf (filename, BUFSIZE, "term_%s.dat", t->term_name) > BUFSIZE)
     return -1;

   if (NULL == (fp = fopen (filename, "w")))
     return -1;

   n = win->num_wave;
   for (i = 0; i < n; i++)
     {
        fprintf (fp, "%5ld %15.6f %15.6e\n",
                 i, win->wave0[i], t->value_workspace[i]);
     }

   return fclose (fp);
}

static int combine_terms (Wavecal_Type *wct, double *model)
{
   Window_Type *win = &wct->window;
   Reference_Irr_Type *irr = &wct->irr;
   Array_Type term_sum[NUM_TERMS];
   int count[NUM_TERMS];
   size_t num_terms = NUM_TERMS;
   size_t i, n = win->num_wave;
   double *ad1, *ad2, *i0, *lbe, *bl;
   Term_Type *t;

   if (0 != alloc_term_arrays (term_sum, num_terms, n))
     return -1;

   /* sum over terms of each type */

   memset ((char *)count, 0, num_terms * sizeof(int));

   for (t = wct->terms; t != NULL; t = t->next)
     {
        double *term_value = t->value_workspace;
        double *sum = term_sum[t->term_type].v;
        for (i = 0; i < n; i++)
          {
             sum[i] += term_value[i];
          }
        count[t->term_type] += 1;
        if (0 != write_term (t, win)) return -1; /* FIXME */
     }

   /* model = sc*(i0 + ad1)*exp(-lbe) + bl + ad2 */

   i0  = irr->irr0_workspace;
   ad1 = term_sum[TERM_TYPE_AD1].v;
   ad2 = term_sum[TERM_TYPE_AD2].v;
   lbe = term_sum[TERM_TYPE_LBE].v;
   bl  = term_sum[TERM_TYPE_BL].v;

   if (count[TERM_TYPE_SC] > 0)
     {
        double *sc  = term_sum[TERM_TYPE_SC].v;
        for (i = 0; i < n; i++)
          {
             model[i] = (sc[i] * (i0[i] + ad1[i]) * exp(-lbe[i])
                         + bl[i] + ad2[i]);
          }
     }
   else
     {
        for (i = 0; i < n; i++)
          {
             model[i] = ((i0[i] + ad1[i]) * exp(-lbe[i])
                         + bl[i] + ad2[i]);
          }
     }

   free_term_arrays (term_sum, num_terms);

   return 0;
}

static int forward_model (Wavecal_Type *wct, const double *params, double *model)
{
   Window_Type *win = &wct->window;
   Shapefun_Type *wl = win->shapefun;
   Reference_Irr_Type *irr = &wct->irr;
   Cspline_Type *cspline = irr->cspline;
   Term_Type *term;
   const double *par;

   par = params;

   /* compute wavelength as a function of pixel index */
   if (wl->st_eval (wl, par, win->num_wave, win->pindex, win->wave0) < 0)
     return -1;
   par += win->num_wave_params;

   /* evaluate the reference irradiance on the new wavelength grid */
   if (cspline_eval (cspline, win->num_wave, win->wave0, irr->irr0_workspace))
     return -1;

   /* evaluate all model terms on the new wavelength grid */
   for (term = wct->terms; term != NULL; term = term->next)
     {
        if (evaluate_term (term, win, par) < 0)
          return -1;
        par += term->num_params;
     }

   /* combine terms to construct the updated model spectrum */
   if (0 != combine_terms (wct, model))
     return -1;

   return 0;
}

static int write_resid (const Window_Type *win,
                        const double *model, const double *spec,
                        const double *resid)
{
   FILE *fp;
   size_t i, n;

   if (Write_Output == 0)
     return 0;

   if (NULL == (fp = fopen ("resid.dat", "w")))
     return -1;

   n = win->num_wave;
   for (i = 0; i < n; i++)
     {
        fprintf (fp, "%5ld %15.6e %15.6f %15.6e %15.6e\n",
                 i, spec[i], win->wave0[i], model[i], resid[i]);
     }

   return fclose (fp);
}

typedef struct
{
   Wavecal_Type *wct;
   const double *spec;
   const double *weight;
   double *model;
   size_t counter;
}
Mpfit_Interface_Type;

/* m = number of functions
 * n = number of parameters
 * x = vector of parameters
 * fvec = vector of function values (residuals)
 * dvec = function derivatives (optional)
 */
static int mpfit_objective_function
(int m, int n, double *x, double *fvec, double **dvec, void *private_data)
{
   Mpfit_Interface_Type *p = (Mpfit_Interface_Type *)private_data;
   const double *spec = p->spec;
   const double *weight = p->weight;
   double *model = p->model;
   int i;

   (void) n; (void) dvec;

   if (forward_model (p->wct, x, model) < 0)
     return -1;

   for (i = 0; i < m; i++)
     {
        fvec[i] = (model[i] - spec[i]) * weight[i];
     }

   if (Write_Output)
     {
        double sumsq = 0.0;
        for (i = 0; i < m; i++)
          {
             double diff = fvec[i];
             sumsq += diff * diff;
          }
        p->counter += 1;
        fprintf(stderr, "sumsq= %12.5e\n", sumsq);
        for (i = 0; i < n; i++)
          {
             fprintf (stderr, "%2d:%12.5e ", i, x[i]);
             if ((i+1)%5 == 0) fputs("\n", stderr);
          }
        fprintf(stderr, "\n");
     }

   return 0;
}

int wavecal_fit (Wavecal_Type *wct, int xtrack,
                 int num_wave, const double *wave,
                 const double *spec, const double *specerr,
                 const Wavecal_Config_Type *config,
                 Wavecal_Result_Type *result)
{
   Mpfit_Interface_Type mp = {0};
   struct mp_config_struct fit_config = {0};
   struct mp_result_struct fit_result = {0};
   Fit_Control_Type *fit_ctrl = &wct->fit_ctrl;
   Window_Type *win = NULL;
   double *params = NULL;
   double *model = NULL;
   double *spec_scaled = NULL;
   double *weight = NULL;
   double *residuals = NULL;
   double scale_factor;
   size_t num;
   int status = -1;
   int i, num_residuals, num_params, mp_status;

   if (0 != init_window (wct, xtrack, num_wave, wave, config))
     return -1;

   win = &wct->window;

   if (collect_params (wct, &num, &params) < 0)
     return -1;
   if (1) fprintf (stderr, "fit: %ld parameters\n", num);

   if ((NULL == (model = alloc_doubles (win->num_wave)))
       || (NULL == (spec_scaled = alloc_doubles (win->num_wave)))
       || (NULL == (weight = alloc_doubles (win->num_wave)))
       || (NULL == (residuals = alloc_doubles (win->num_wave)))) /* FIXME? */
     goto return_error;

   /* Scale the radiance and irradiance by the same factor */
   scale_factor = wct->irr.file.scale_factor;
   for (i = 0; i < win->num_wave; i++)
     {
        double err = specerr[i];
        spec_scaled[i] = spec[i] * scale_factor;
        if ((err != 0.0) && (isfinite(err) != 0))
          {
             weight[i] = 1.0/(fabs(err) * scale_factor);
          }
        else weight[i] = 0.0;
     }

   mp.wct = wct;
   mp.spec = spec_scaled;
   mp.weight = weight;
   mp.model = model;
   mp.counter = 0;

   num_params = num;
   num_residuals = win->num_wave;

   fit_config.xtol = fit_ctrl->xtol;
   fit_config.ftol = fit_ctrl->ftol;
   fit_config.maxiter = fit_ctrl->maxiter;
   fit_config.maxfev = (fit_ctrl->maxfev ?
                        fit_ctrl->maxfev : fit_ctrl->maxiter * num_params);

   fit_result.resid = residuals; /* FIXME? */

   mp_status = mpfit (&mpfit_objective_function, num_residuals,
                      num_params, params, NULL, &fit_config, &mp, &fit_result);
   fprintf (stderr, "mpfit returned mp_status = %d\n", mp_status);

   fprintf (stderr, "result:\n bestnorm=%g\n orignorm=%g\n niter=%d nfev=%d\n",
            fit_result.bestnorm,
            fit_result.orignorm,
            fit_result.niter,
            fit_result.nfev
           );

   if (1)
   {  /* FIXME - making sure the model is evaluated at best-fit */
      Write_Output = 1;
      if (0 != mpfit_objective_function (num_residuals, num_params,
                                         params, residuals, NULL, &mp))
        goto return_error;
      /* FIXME */
      if (write_resid (win, model, spec_scaled, residuals))
        return -1;
   }

   /* FIXME - what to return?
    * updated wavelengths: win->wave0?
    * best-fit wavelength parameters?
    */
   if (result)
     {
        result->wave = win->wave0;
     }

   status = 0;
return_error:
   FREE(params);
   FREE(model);
   FREE(spec_scaled);
   FREE(weight);
   FREE(residuals);

   return status;
}
