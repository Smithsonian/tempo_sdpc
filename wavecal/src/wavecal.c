/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <tell.h>
#include <tio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <mpfit.h>

#include "config.h"
#include "util.h"
#include "interp_cspline.h"
#include "interp.h"
#include "shapefun.h"
#include "wavecal.h"
#include "slit_function.h"
#include "slit_function_asg.h"
#include "slit_function_table.h"

#define NCID_UNINITIALIZED -1

#define MIN_ACCEPTABLE_GOOD_PIXEL_FRACTION  0.8

enum
{
   TERM_TYPE_AD1,
   TERM_TYPE_AD2,
   TERM_TYPE_LBE,
   TERM_TYPE_SC,
   TERM_TYPE_BL
};

#define NUM_TERM_TYPES 5

enum
{
   SF_MODE_NONE,
   SF_MODE_APPLY,
   SF_MODE_FIT
};

typedef struct
{
   char *sf_path;           /**< slit function lookup table */
   double *initial_params;  /**< initial slit function parameter values for fit. (optional - null is ok) */
   double *param_step;      /**< slit function parameter step values for numerical derivatives when fitting. (optional - null is ok) */
   int num_params;
   int num_pad_half_widths;
   int mode;
}
SF_Control_Type;

typedef struct
{
   double *model_padded;
   double *model_convolved;
   double *derivs_convolved[SFT_MAX_NUM_PARAMS];
   double *params;   /* storage for fitted slit function parameter result */
   double *params_error;  /* storage for parameter uncertainty estimate */
   int num_params;
   int num_alloc;    /* allocated size of model_padded/model_convolved arrays */
   int num_waves;    /* number of wavelengths in the unpadded model */
   int num_pad;      /* number of wavelengths zero-padding on one end */
}
SF_Convolution_Type;

typedef struct
{
   int *sample_indices;
   int num_samples;
}
Sample_Type;

typedef struct
{
   Sample_Type *sample_model;
   Sample_Type *sample_data;
   double *wave0_samp;
   double *spec_scaled_samp;
   double *weight_samp;
   double *waves_eval_samp;
   double *combined_model_samp;
   double *irr_value_samp;
   double *convolved_model;
   double isrf_pad_frac;
}
Window_Sample_Type;

typedef struct
{
   int xtrack;      /**< cross-track index to be fitted */
   size_t num_wave; /**< total number of wavelength points in measured spectrum */
   double *wave0;   /**< initial guess at wavelength grid for measured spectrum */

   double *xerror;  /**< storage to hold parameter uncertainty diagnostic */

   Cspline_Type *waves_nominal_cspline;   /**< nominal wavelength grid */
   double *waves_nominal_scratch;         /**< scratch space for computing grids */
   int adjust_nominal_wavelength;         /**< zero means NO, non-zero means YES */

   Shapefun_Type *wavegrid_shapefun;    /**< shape function for computing the wavelength grid */
   double *pindex;          /**< pixel index array for measured spectrum */
   double *wave_params;     /**< wavelength grid parameters */
   size_t num_wave_params;  /**< number of wavelength grid parameters */

   int start_pix;               /**< offset to sub-window to be fitted */
   double delta_wavelength;     /**< bin width */
   double feature_wavelength;   /**< wavelength of fiducial feature */
   double rad_mean_ratio;       /**< (mean_rad)/(mean_irr) within target wavelength band */

   /* storage for use during the fit iteration */
   double *model;           /**< computed model */
   double *spec_scaled;     /**< spectrum to calibrate, scaled by a constant */
   double *weight;          /**< weight for each spectrum pixel */
   double *residuals;       /**< weighted fit residual, (model - spec_scaled)*weight */

   SF_Convolution_Type *sfct;   /**< extra storage to perform the slit-function convolution */

   Window_Sample_Type samp;     /**< structure to support sub-sampling the data */
}
Window_Type;

typedef struct
{
   char *path;             /**< path to the file */
   const char *name_x;     /**< name of variable containing "X" */
   const char *name_y;     /**< name of variable containing "Y" */
   double scale_factor;    /**< scale factor to be applied to "Y" values */
}
File_Type;

typedef struct
{
   File_Type file;         /**< file specification */
   Interp_Type *interp;    /**< generic interpolation object */
}
Refspec_Type;

typedef struct Term_Type Term_Type;
struct Term_Type
{
   Term_Type *next;
   char *term_name;                   /**< user-specified name string */
   int term_type;                     /**< term type (enum) */
   Refspec_Type refspec;              /**< reference spectrum, if any */
   Shapefun_Type *shapefun;           /**< shape function, or NULL */
   Shapefun_Init_Type shapefun_init;  /**< shape function initialization data */
   double *params;                    /**< initial fit parameters for this term, if any */
   double *eval_params;               /**< params from last evaluation */
   size_t num_params;                 /**< number of fit parameters for this term */
   double *value;                     /**< term value vs wavelength */
   size_t num_values;                 /**< number of wavelengths in value array */
};

typedef struct
{
   File_Type file;                   /**< location of reference irradiance (file, variables) */
   int ncid;                         /**< netcdf file descriptor */
   Cspline_Type *cspline;            /**< spline interpolation object */
   double *irradiance;               /**< reference irradiance values */
   double *wavelen;                  /**< reference irradiance wavelength grid */
   size_t num_wavelen;               /**< number of wavelength points in grid */
   double delta_wavelength;          /**< reference irradiance wavelength grid spacing */
}
Reference_Irr_Type;

typedef struct
{
   double xtol;        /**< fit parameter value convergence criterion */
   double ftol;        /**< objective function convergence criterion */
   int maxiter;        /**< max number of iterations of the optimization algorithm */
   int maxfev;         /**< max number of evaluations of the objective function */
}
Fit_Control_Type;

struct Wavecal_Type
{
   Reference_Irr_Type irr;      /**< reference irradiance */
   Window_Type window;          /**< fit window */
   Fit_Control_Type fit_ctrl;   /**< optimization control parameters */
   TIO_Meta_Type *meta;         /**< metadata recording object */

   Slit_Function_Type *sft;     /**< slit function object */
   SF_Table_Type *sf_table;     /**< slit function lookup table (optional) */
   SF_Control_Type sf_ctrl;

   Term_Type *terms;            /**< terms in the model being fitted */
   double *term_sums[NUM_TERM_TYPES];   /**< sum over terms within each term type */
   double *irr0;                /**< reference irradiance interpolated onto target spectrum wavelength grid */
   int is_irradiance;
   int xtrack;                  /**< slit function lookup table requires xtrack index */
};

static int select_ref_irr_subset (const Reference_Irr_Type *irr, double wave_beg, double wave_end, int num_pad,
                                  int *pindex_irr_beg, int *pindex_irr_end);

static void sample_free (Sample_Type *s)
{
   if (s == NULL) return;
   FREE(s->sample_indices);
   FREE(s);
}

static Sample_Type *sample_alloc (int n)
{
   Sample_Type *s = NULL;

   if (NULL == (s = (Sample_Type *)MALLOC (sizeof *s)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   s->num_samples = n;
   s->sample_indices = NULL;

   if (NULL == (s->sample_indices = (int *)MALLOC (n * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        sample_free (s);
        return NULL;
     }

   return s;
}

static Sample_Type *sample_dup (const Sample_Type *s)
{
   Sample_Type *dup = NULL;
   if (NULL == (dup = sample_alloc (s->num_samples)))
     return NULL;
   memcpy ((char *)dup->sample_indices, (char *)s->sample_indices, s->num_samples * sizeof(int));
   return dup;
}

#if 0
static Sample_Type *sample_interval (int num, int first, int interval)
{
   Sample_Type *s = NULL;
   int i, nsamples = num / interval;
   int offset = first + (num/2) % interval;

   if (NULL == (s = sample_alloc (nsamples)))
     return NULL;

   for (i = 0; i < nsamples; i++)
     {
        s->sample_indices[i] = offset + i * interval;
     }

   return s;
}
#endif

static Sample_Type *sample_threshold (int num, const double *x, double thresh)
{
   Sample_Type *s = NULL;
   int i, k, nsamples;

   nsamples = 0;
   for (i = 0; i < num; i++)
     {
        if (x[i] > thresh) nsamples++;
     }

   if (NULL == (s = sample_alloc (nsamples)))
     return NULL;

   k = 0;
   for (i = 0; i < num; i++)
     {
        if (x[i] > thresh)
          {
             s->sample_indices[k++] = i;
          }
     }

   return s;
}

static int sample_gather (const Sample_Type *s, const double *x, double *xsamp, int *num_xsamp)
{
   int i;
   if (num_xsamp) *num_xsamp = s->num_samples;
   for (i = 0; i < s->num_samples; i++)
     {
        int k = s->sample_indices[i];
        xsamp[i] = x[k];
     }
   return 0;
}

static int sample_scatter (const Sample_Type *s, const double *xsamp, double *x)
{
   int i;
   for (i = 0; i < s->num_samples; i++)
     {
        int k = s->sample_indices[i];
        x[k] = xsamp[i];
     }
   return 0;
}

static int init_data_sampling (Window_Type *win, int interval, int half_width)
{
   Window_Sample_Type *samp = &win->samp;
   double *mask = NULL;
   int i0, i, num_samples;
   int num = win->num_wave;

   if (interval <= 1)
     return 0;

   if (half_width < 0) half_width = 0;

   if (NULL == (mask = (double *)MALLOC (num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)mask, 0, num * sizeof(double));

   /* Construct a mask that defines sample windows each spanning
    * (2*half_width+1) pixels and with window centers spaced a
    *  distance 'interval' apart.
    * The sample windows are offset from the end of the mask to
    * center the entire mask pattern in the user-specified data range.
    */

   i0 = win->start_pix + (num/2) % interval;
   for (i = 0; (i + i0) < num; i++)
     {
        int k = abs((i % interval) - half_width);
        if (k <= half_width) mask[i+i0] = 1.0;
     }

   samp->sample_data = sample_threshold (num, mask, 0.0);
   FREE(mask);
   if (NULL == samp->sample_data)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initializing data subsampling", __func__);
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1,
              "data sampling: interval=%d half_width=%d num_samples=%d",
              interval, half_width, samp->sample_data->num_samples);

   /* allocate additional temporary workspace based on the number of samples */
   num_samples = samp->sample_data->num_samples;

   if ((NULL == (samp->wave0_samp = alloc_doubles (num_samples)))
       || (NULL == (samp->spec_scaled_samp = alloc_doubles (num_samples)))
       || (NULL == (samp->weight_samp = alloc_doubles (num_samples))))
     return -1;

   return 0;
}

static int init_model_sampling (Window_Type *win, int sf_mode, const double *waves,
                                const Reference_Irr_Type *irr)
{
   SF_Convolution_Type *sfct = win->sfct;
   Window_Sample_Type *samp = &win->samp;
   double *irr_wavelen;
   double *mask = NULL;
   double wave_beg, wave_end;
   int i, index_irr_beg, index_irr_end, num_irr_waves, num_samples, num_pad;
   int *sample_indices = NULL;
   int status = -1;

   /* Quick return when slit function convolution is turned off --
    * Model and data are on the same grid */
   if (sf_mode == SF_MODE_NONE)
     {
        if ((NULL == (samp->sample_model = sample_dup (samp->sample_data)))
            || (NULL == (samp->waves_eval_samp = alloc_doubles (samp->sample_model->num_samples))))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: initializing model subsampling", __func__);
             return -1;
          }
        tell_vlog (TELL_MSGTYPE_INFO, 1,
                   "model sampling: num_samples=%d", samp->sample_model->num_samples);
        return 0;
     }

   /* Select the relevant subset of the reference irradiance wavelength grid */
   if (0 != select_ref_irr_subset (irr, win->wave0[0], win->wave0[win->num_wave-1], sfct->num_pad, &index_irr_beg, &index_irr_end))
     return -1;

   /* define irradiance subset pointers */
   num_irr_waves = index_irr_end - index_irr_beg + 1;
   irr_wavelen = irr->wavelen + index_irr_beg;

   wave_beg = irr_wavelen[0];
   wave_end = irr_wavelen[num_irr_waves-1];

   /* Ideally we should use the slit function to determine which model wavelengths
    * contribute to each sampled data pixel for the current slit function, but the
    * "exact" mapping could change with every fit iteration.  However, that's
    * needlessly complicated, because an approximate mapping is good enough.
    *
    * Instead, we sample the model in a (user adjustable) fixed-size window around
    * each sampled data wavelength. Increasing the window size can reduce the sensitivity
    * of the final result to this assumption.
    */
   if (NULL == (mask = alloc_doubles (num_irr_waves)))
     goto return_status;
   memset ((char *)mask, 0, num_irr_waves * sizeof(double));

   num_samples = samp->sample_data->num_samples;
   sample_indices = samp->sample_data->sample_indices;

   num_pad = sfct->num_pad;
   if (samp->isrf_pad_frac > 0.0) num_pad *= samp->isrf_pad_frac;

   /* Mark each model wavelength that contributes to a sampled data pixel */
   for (i = 0; i < num_samples; i++)
     {
        int k, i_irr, jbeg, jend, j;
        double wave_k;

        k = sample_indices[i];
        wave_k = waves[k];

        if ((wave_k < wave_beg) || (wave_end < wave_k))
          continue;
        i_irr = bsearch_d (wave_k, irr_wavelen, num_irr_waves);

        jbeg = i_irr - num_pad;
        jend = i_irr + num_pad;

        if (jbeg < 0) jbeg = 0;
        if (jend > num_irr_waves) jend = num_irr_waves;

        for (j = jbeg; j < jend; j++)
          {
             mask[j] += 1;
          }
     }

   /* The mask now marks all model wavelengths that map to the sampled data pixels */
   if (NULL == (samp->sample_model = sample_threshold (num_irr_waves, mask, 0.0)))
     goto return_status;

   tell_vlog (TELL_MSGTYPE_INFO, 1,
              "model sampling: isrf_pad_frac=%0.2g num_pad=%d num_samples=%d",
              samp->isrf_pad_frac, num_pad, samp->sample_model->num_samples);

   /* allocate additional temporary workspace based on the number of samples */
   if ((NULL == (samp->waves_eval_samp = alloc_doubles (samp->sample_model->num_samples)))
       || (NULL == (samp->combined_model_samp = alloc_doubles (samp->sample_model->num_samples)))
       || (NULL == (samp->irr_value_samp = alloc_doubles (sfct->num_alloc)))
       || (NULL == (samp->convolved_model = alloc_doubles (win->num_wave))))
     goto return_status;

   status = 0;
return_status:
   FREE(mask);

   return status;
}

static void free_window_sample (Window_Sample_Type *samp)
{
   if (samp == NULL)
     return;
   sample_free (samp->sample_model);
   sample_free (samp->sample_data);
   FREE(samp->wave0_samp);
   FREE(samp->spec_scaled_samp);
   FREE(samp->weight_samp);
   FREE(samp->waves_eval_samp);
   FREE(samp->combined_model_samp);
   FREE(samp->irr_value_samp);
   FREE(samp->convolved_model);
}

static int init_sampling (Window_Type *win, config_setting_t *s_band)
{
   config_setting_t *s_samp;
   double isrf_pad_frac;
   int interval, half_width;

   if (NULL == (s_samp = config_setting_get_member (s_band, "sampling")))
     return 0;

   if (CONFIG_TRUE != config_setting_lookup_int (s_samp, "interval", &interval))
     interval = 1;

   if (CONFIG_TRUE != config_setting_lookup_int (s_samp, "half_width", &half_width))
     half_width = 0;

   if (CONFIG_TRUE != config_setting_lookup_float (s_samp, "isrf_pad_frac", &isrf_pad_frac))
     isrf_pad_frac = 1.0;

   win->samp.isrf_pad_frac = isrf_pad_frac;

   /* Data sampling is independent of the slit function, so we can
    * initialize the necessary structures here.
    * Model sampling must be deferred until after the slit function
    * is initialized. */
   return init_data_sampling (win, interval, half_width);
}

static void free_file_type (File_Type *file)
{
   if (file == NULL)
     return;
   FREE(file->path);
}

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
   FREE(tt->eval_params);
   FREE(tt->value);
   free_file_type(&tt->refspec.file);
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
   const char *path = NULL;

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "path", &path))
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

   if (NULL == (ct->path = expand_path (path)))
     return -1;

   return 0;
}

static void tell_config_error (config_setting_t *s, const char *func)
{
   tell_verror (TELL_INVALID_PARM_ERROR, "%s: reading from %s (%s:%d)",
                func, config_setting_name(s),
                config_setting_source_file(s),
                config_setting_source_line(s));
}

typedef struct
{
   const char *mode_name;
   int mode_index;
}
SF_Mode_Table_Entry;
static SF_Mode_Table_Entry SF_Mode_Table[] =
{
   {"none",  SF_MODE_NONE},
   {"apply", SF_MODE_APPLY},
   {"fit",   SF_MODE_FIT},
   {NULL, -1}
};

static int config_slit_function (config_setting_t *s, SF_Control_Type *sf_ctrl)
{
   SF_Mode_Table_Entry *m = NULL;
   const char *mode;
   const char *path;
   int num_params;
   size_t ns = 0;

   if (CONFIG_TRUE != config_setting_lookup_string (s, "mode", &mode))
     return -1;

   for (m = SF_Mode_Table; m->mode_name != NULL; m++)
     {
        if (0 == strcmp (m->mode_name, mode))
          {
             sf_ctrl->mode = m->mode_index;
             break;
          }
     }
   if (m->mode_name == NULL)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported slit function model: %s",
                     __func__, mode);
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "sf_path", &path))
     {
        sf_ctrl->sf_path = NULL;
     }
   else
     {
        if (NULL == (sf_ctrl->sf_path = expand_path (path)))
          return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_isrf_half_widths", &sf_ctrl->num_pad_half_widths))
     return -1;

   /* When the SF is a lookup table, the params and param_step arrays are irrelevant,
    * but we still need to know how many parameters are in the table before we actually
    * read the table (because it's convenient to pre-allocate some structures and working
    * space before we actually read the SF parameter table) */
   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_params", &num_params))
     return -1;
   sf_ctrl->num_params = num_params;

   /* initial parameter values (FIXME - make this optional?) */
   if (read_config_float_array (s, "params", &sf_ctrl->initial_params, &ns) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error reading slit function params array from %s",
                     __func__, config_setting_source_file(s));
        return -1;
     }
   if ((ns > 0) && (ns != (size_t) sf_ctrl->num_params))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: inconsistent array size: params (in %s)",
                     __func__, config_setting_source_file (s));
        return -1;
     }

   /* param step sizes for numerical derivative calculation */
   if (read_config_float_array (s, "param_step", &sf_ctrl->param_step, &ns) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error reading slit function param_step array from %s",
                     __func__, config_setting_source_file(s));
        return -1;
     }

   if ((ns > 0) && (ns != (size_t) sf_ctrl->num_params))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: inconsistent array size: param_step (in %s)",
                     __func__, config_setting_source_file (s));
        return -1;
     }

   return 0;
}

static int config_interp_method (config_setting_t *s, Interp_Type *it)
{
   int num_coef;

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

   (void) config_setting_lookup_bool (s, "scale_by_mean_radiance_over_mean_irradiance",
                                      &st->st_apply_external_scaling);

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
   const char *method_name = "cspline";
   config_setting_t *ss;
   Interp_Type *interp = NULL;

   *pinterp = NULL;

   /* optional */
   if (NULL == (ss = config_setting_get_member (s, "data")))
     return 0;

   if (0 != read_filepar (ss, file))
     return -1;

#if 0
   if (CONFIG_FALSE == config_setting_lookup_string (s, "method", &method_name))
     {
        tell_config_error (s, __func__);
        return -1;
     }
#endif

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

static void free_sf_convolution_type (SF_Convolution_Type *sfct)
{
   int i;
   if (NULL == sfct)
     return;
   FREE(sfct->model_padded);
   FREE(sfct->model_convolved);
   FREE(sfct->params);
   FREE(sfct->params_error);
   if (sfct->derivs_convolved)
     {
        for (i = 0; i < sfct->num_params; i++)
          {
             FREE(sfct->derivs_convolved[i]);
          }
     }

   FREE(sfct);
}

static SF_Convolution_Type *alloc_sf_convolution_type (int num_waves, int num_pad, int num_params)
{
   SF_Convolution_Type *sfct = NULL;
   size_t len = num_waves + 2*num_pad;
   int i;

   if (NULL == (sfct = (SF_Convolution_Type *)MALLOC (sizeof *sfct)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   sfct->num_waves = 0;
   sfct->num_pad = num_pad;
   sfct->num_alloc = len;
   sfct->num_params = num_params;

   if ((NULL == (sfct->model_padded = (double *)MALLOC (len * sizeof(double))))
       || (NULL == (sfct->model_convolved = (double *)MALLOC (len * sizeof(double))))
       || (NULL == (sfct->params = (double *)MALLOC (num_params * sizeof(double))))
       || (NULL == (sfct->params_error = (double *)MALLOC (num_params * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_sf_convolution_type (sfct);
        return NULL;
     }

   for (i = 0; i < sfct->num_params; i++)
     {
        if ((NULL == (sfct->derivs_convolved[i] = (double *)MALLOC (len * sizeof(double)))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_sf_convolution_type (sfct);
             return NULL;
          }
        memset ((char *)sfct->derivs_convolved[i], 0, len * sizeof(double));
     }

   memset ((char *)sfct->model_padded, 0, len * sizeof(double));
   memset ((char *)sfct->model_convolved, 0, len * sizeof(double));

   return sfct;
}

static void free_window (Window_Type *win)
{
   if (win == NULL)
     return;
   FREE(win->wave0);
   FREE(win->xerror);
   cspline_free (win->waves_nominal_cspline);
   FREE(win->waves_nominal_scratch);
   FREE(win->pindex);
   FREE(win->model);
   FREE(win->spec_scaled);
   FREE(win->weight);
   FREE(win->residuals);
   FREE(win->wave_params);
   free_window_sample (&win->samp);
   free_shapefun_type (win->wavegrid_shapefun);
   free_sf_convolution_type (win->sfct);
}

static int alloc_window (Window_Type *win, int num_data_waves, int num_model_waves,
                         int num_pad, int num_sf_params)
{
   win->num_wave = (size_t) num_data_waves;

   if ((NULL == (win->wave0 = alloc_doubles (num_data_waves)))
       || (NULL == (win->waves_nominal_scratch = alloc_doubles (num_data_waves)))
       || (NULL == (win->pindex = alloc_doubles (num_data_waves)))
       || (NULL == (win->model = alloc_doubles (num_data_waves)))
       || (NULL == (win->spec_scaled = alloc_doubles (num_data_waves)))
       || (NULL == (win->weight = alloc_doubles (num_data_waves)))
       || (NULL == (win->residuals = alloc_doubles (num_data_waves))))
     return -1;

   win->xerror = NULL;        /* set this later */

   memset ((char *)&win->samp, 0, sizeof (win->samp));

   if (num_pad)
     {
        if (NULL == (win->sfct = alloc_sf_convolution_type (num_model_waves, num_pad, num_sf_params)))
          return -1;
     }

   return 0;
}

static void free_reference_irr_type (Reference_Irr_Type *irr)
{
   if (irr == NULL)
     return;
   if (irr->ncid) TIO_close (irr->ncid);
   cspline_free (irr->cspline);
   free_file_type (&irr->file);
   FREE(irr->wavelen);
   FREE(irr->irradiance);
}

static void free_term_sums (double **ts, size_t num_term_types)
{
   size_t i;

   if (ts == NULL)
     return;

   for (i = 0; i < num_term_types; i++)
     {
        FREE(ts[i]);
     }
}

static int alloc_term_sums (double **a, size_t na, size_t num)
{
   size_t i;

   for (i = 0; i < na; i++)
     {
        if (NULL == (a[i] = alloc_doubles (num)))
          return -1;
     }

   return 0;
}

static void zero_term_sums (double **a, size_t na, size_t num)
{
   size_t i, len = num * sizeof(double);
   for (i = 0; i < na; i++)
     {
        memset ((char *)a[i], 0, len);
     }
}

static int alloc_term_storage (Term_Type *term, int num_wave)
{
   Term_Type *t;

   for (t = term; t != NULL; t = t->next)
     {
        double *v;
        /* we may have 2 factors, so allocate space for both */
        if (NULL == (v = alloc_doubles (2 * num_wave)))
          return -1;
        t->value = v;
        t->num_values = num_wave;
     }

   return 0;
}

static void free_sf_ctrl (SF_Control_Type *sf_ctrl)
{
   if (sf_ctrl == NULL)
     return;
   FREE(sf_ctrl->initial_params);
   FREE(sf_ctrl->param_step);
   FREE(sf_ctrl->sf_path);
}

static void free_wavecal (Wavecal_Type *wct)
{
   size_t num_term_types = NUM_TERM_TYPES;
   if (wct == NULL)
     return;
   free_reference_irr_type (&wct->irr);
   free_terms (wct->terms);
   free_term_sums (wct->term_sums, num_term_types);
   free_window (&wct->window);
   sft_free (wct->sft);
   free_sf_ctrl (&wct->sf_ctrl);
   if (wct->sf_table)
     {
        SF_Table_Type *stt = wct->sf_table;
        stt->stt_close (stt);
   }
   FREE(wct->irr0);
   FREE(wct);
}

void wavecal_close (Wavecal_Type *wct)
{
   free_wavecal (wct);
}

static Wavecal_Type *alloc_wavecal (int num_data_waves, int num_model_waves,
                                    int num_pad, int num_sf_params)
{
   Wavecal_Type *wct = NULL;
   size_t num_term_types = NUM_TERM_TYPES;

   if (NULL == (wct = (Wavecal_Type *)MALLOC (sizeof *wct)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)wct, 0, sizeof *wct);

   if (NULL == (wct->irr0 = alloc_doubles (num_data_waves)))
     {
        free_wavecal (wct);
        return NULL;
     }

   if (0 != alloc_term_sums (wct->term_sums, num_term_types, num_model_waves))
     {
        free_wavecal (wct);
        return NULL;
     }

   if (0 != alloc_window (&wct->window, num_data_waves, num_model_waves, num_pad, num_sf_params))
     {
        free_wavecal (wct);
        return NULL;
     }

   if (num_pad)
     {
        if (NULL == (wct->sft = sft_new (num_pad * 2, num_sf_params)))
          {
             free_wavecal (wct);
             return NULL;
          }
     }

   /* default this parameter to irradiance calibration case */
   wct->window.rad_mean_ratio = 1.0;

   return wct;
}

static int config_model_components (Wavecal_Type *wct, config_setting_t *s)
{
   Term_Type *term;
   config_setting_t *s_terms;
   config_setting_t *ss;
   unsigned int i, num;

   if (NULL == (s_terms = config_setting_lookup (s, "model_components")))
     return 0;

   num = config_setting_length (s_terms);

   for (i = 0; i < num; i++)
     {
        const char *name;
        if (NULL == (ss = config_setting_get_elem (s_terms, i)))
          return -1;

        name = config_setting_name (ss);
        if (name[0] == '*') continue;

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
   TIO_Var_Info_Type info = {0};
   config_setting_t *ss;
   int start=0, count=1;

   if (NULL == (ss = config_setting_get_member (s, "data")))
     return -1;
   if (0 != read_filepar (ss, file))
     return -1;

   if (0 != TIO_open (file->path, NC_NOWRITE, &irr->ncid))
     return -1;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file->path);

   if (0 != TIO_inq_var (irr->ncid, file->name_y, &info))
     return -1;
   irr->num_wavelen = (info.ndims == 1) ? info.dimlens[0] : info.dimlens[1];

   if (0 != TIO_get_var_section (irr->ncid, "dWavelength", &start, &count, TIO_DOUBLE, &irr->delta_wavelength))
     return -1;

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

   (void) config_setting_lookup_bool (s, "adjust_nominal_wavelength", &win->adjust_nominal_wavelength);

   if (NULL == (win->wavegrid_shapefun = shapefun_create (method_name)))
     return -1;

   if (0 != config_shapefun_method (s, win->wavegrid_shapefun, NULL))
     return -1;

   st = win->wavegrid_shapefun;
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

static int read_irr_reference (Reference_Irr_Type *irr, TIO_Meta_Type *meta, int xtrack)
{
   File_Type *file = &irr->file;
   double dx0;
   size_t i;

   if (0)
     {
        if (0 != meta_record_basename (meta, file->path))
          return -1;
     }

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
        for (i = 0; i < irr->num_wavelen; i++)
          {
             irr_i[i] *= scale;
          }
     }

   /* Check for uniform grid spacing */
   dx0 = irr->wavelen[1] - irr->wavelen[0];
   for (i = 2; i < irr->num_wavelen; i++)
     {
        double dx = irr->wavelen[i] - irr->wavelen[i-1];
        /* FIXME - fix reference data to enforce tighter tolerances!? */
        if (gsl_fcmp (dx, dx0, 1.e-2) != 0)
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: reference irradiance has variable grid spacing: dx0=%e dx[%ld]=%e",
                          __func__, dx0, i, dx);
             return -1;
          }
     }
   irr->delta_wavelength = dx0;

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
        if (NULL == (term->eval_params = alloc_doubles (term->num_params)))
          return -1;
        /* we haven't evaluated the function yet, so init this array to zero */
        memset ((char *)term->eval_params, 0, term->num_params * sizeof(double));
     }
   else
     {
        term->params = NULL;
        term->eval_params = NULL;
     }

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

static int read_rad_reference (Term_Type *terms, TIO_Meta_Type *meta, int xtrack)
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

        if (0)
          {
             if (0 != meta_record_basename (meta, term->refspec.file.path))
               return -1;
          }
     }

   return 0;
}

static int init_window_reference_spectra (Wavecal_Type *wct, int xtrack)
{
   Reference_Irr_Type *irr = &wct->irr;
   Window_Type *win = &wct->window;

   win->xtrack = xtrack;

   if (0 != read_rad_reference (wct->terms, wct->meta, xtrack))
     return -1;

   if (0 != read_irr_reference (&wct->irr, wct->meta, xtrack))
     return -1;

   cspline_free (irr->cspline); /* FIXME - preallocate this */
   if (NULL == (irr->cspline = cspline_interpol (irr->num_wavelen, irr->wavelen,
                                                 irr->irradiance, NULL)))
     return -1;

   return 0;
}

static int init_window_shapefun (Wavecal_Type *wct, const double *wave)
{
   Window_Type *win = &wct->window;
   Shapefun_Type *shapefun = win->wavegrid_shapefun;
   Shapefun_Init_Type shapefun_init = {0};

   memcpy ((char *)win->wave0, (char *)wave, win->num_wave * sizeof(double));

   shapefun_init.x = win->pindex;
   shapefun_init.y = win->wave0;
   shapefun_init.n = win->num_wave;

   shapefun->xmin = shapefun_init.x[0];
   shapefun->xmax = shapefun_init.x[win->num_wave-1];

   if (win->adjust_nominal_wavelength)
     {
        /* FIXME - may want to support user-specified initial parameters.
         * But the model is usually a polynomial, so starting with all coefficients zero
         * should usually be the same as starting with the nominal wavelength grid.
         */
        memset ((char *)win->wave_params, 0, win->num_wave_params * sizeof(double));
        return 0;
     }
   else
     {
        return shapefun->st_init_params (shapefun, &shapefun_init,
                                         win->num_wave_params, win->wave_params);
     }
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

   /* put slit function parameters at the very end */
   if (wct->sf_ctrl.mode == SF_MODE_FIT)
     {
        num += wct->sf_ctrl.num_params;
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

   if ((wct->sf_ctrl.mode == SF_MODE_FIT)
       && (wct->sf_ctrl.initial_params != NULL))
     {
        /* set initial slit-function parameters  */
        memcpy ((char *)pnext, (char *)wct->sf_ctrl.initial_params,
                wct->sf_ctrl.num_params * sizeof(double));
     }

   *pnum = num;
   *pparams = params;

   return 0;
}

typedef struct
{
   double feature_wavelength;
   double delta_wavelength;
   int start_pix;
   int num_pix;
}
Feature_Window_Type;

static int read_feature_window (config_setting_t *s_band, Feature_Window_Type *fwin)
{
   config_setting_t *s;

   if (NULL == (s = config_setting_get_member (s_band, "feature_window")))
     return 1;

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "num_pix_fit", &fwin->num_pix))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "start_pix", &fwin->start_pix)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading feature_window in param file, %s:%d",
                     __func__, config_setting_source_file (s), config_setting_source_line (s));
          return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "delta_wavelength", &fwin->delta_wavelength))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "fid_wavelength", &fwin->feature_wavelength)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading feature_window in param file, %s:%d",
                     __func__, config_setting_source_file (s), config_setting_source_line (s));
          return -1;
     }

   return 0;
}

int wavecal_num_wave_params (const Wavecal_Type *wct)
{
   if (wct == NULL)
     return -1;
   return wct->window.num_wave_params;
}

int wavecal_fitting_sf_params (const Wavecal_Type *wct)
{
   if (wct == NULL)
     return -1;
   return (wct->sf_ctrl.mode == SF_MODE_FIT);
}

int wavecal_query_feature_window (const Wavecal_Type *wct, int *start_pix, int *num_pix,
                                  int *adjust_nominal_wavelength)
{
   const Window_Type *win = NULL;
   if (wct == NULL)
     return -1;
   win = &wct->window;
   if (start_pix) *start_pix = win->start_pix;
   if (num_pix) *num_pix = win->num_wave;
   if (adjust_nominal_wavelength) *adjust_nominal_wavelength = win->adjust_nominal_wavelength;
   return 0;
}

Wavecal_Type *wavecal_open (config_t *cfg, const char *cfg_name, TIO_Meta_Type *meta,
                            int max_num_data_waves, int is_irradiance)
{
   char *slit_function_setting = "slit_function";
   Wavecal_Type *wct = NULL;
   Window_Type *win = NULL;
   Feature_Window_Type fwin = {0};
   SF_Control_Type sf_ctrl = {0};
   Reference_Irr_Type ref_irr = {0};
   config_setting_t *s, *s_band, *s_irr, *s_slit;
   config_setting_t *s_rad = NULL;
   int i, num_data_waves, num_model_waves, num_pad;

   /* By default, operate on the entire spectrum. */
   fwin.num_pix = max_num_data_waves;
   fwin.start_pix = 0;
   fwin.delta_wavelength = 0.0;
   fwin.feature_wavelength = 0.0;

   /* Locate the configuration parameters for this spectral band, and target spectrum type */

   if (NULL == (s_band = config_lookup (cfg, cfg_name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing %s in param file: %s",
                     __func__, cfg_name, config_error_file (cfg));
        return NULL;
     }

   if (NULL == (s_irr = config_setting_get_member (s_band, "wavecal_irradiance")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing wavecal_irradiance in param file: %s",
                     __func__, config_error_file (cfg));
        goto error_return;
     }

   if (is_irradiance)
     {
        s_slit = config_setting_get_member (s_irr, slit_function_setting);  /* NULL is ok */
     }
   else
     {
        if (NULL == (s_rad = config_setting_get_member (s_band, "wavecal_radiance")))
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: accessing wavecal_radiance in param file: %s",
                          __func__, config_error_file (cfg));
             goto error_return;
          }
        s_slit = config_setting_get_member (s_rad, slit_function_setting);  /* NULL is ok */
     }

   if (read_feature_window (s_band, &fwin) < 0)
     goto error_return;

   if (((fwin.num_pix <= 0) || (fwin.num_pix > max_num_data_waves))
       || ((fwin.start_pix < 0) || (fwin.start_pix >= max_num_data_waves)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: invalid wavelength calibration window size: max_num_data_waves=%d num_pix=%d start_pix=%d",
                     __func__, max_num_data_waves, fwin.num_pix, fwin.start_pix);
        return NULL;
     }

   /* We always use a reference irradiance spectrum */
   if (0 != config_irr_reference (s_irr, &ref_irr))
     goto error_return;

   /* The slit function is optional */
   if (s_slit)
     {
        if (0 != config_slit_function (s_slit, &sf_ctrl))
          {
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: reading slit_function setting in %s",
                          __func__, config_setting_source_file (s_slit));
             goto error_return;
          }
     }
   else sf_ctrl.mode = SF_MODE_NONE;

   /* The slit function mode and the fit window size together determine
    * the size of the arrays used to compute model components.
    */

   num_data_waves = fwin.num_pix;

   if (sf_ctrl.mode == SF_MODE_NONE)
     {
        num_pad = 0;
        num_model_waves = fwin.num_pix;
     }
   else
     {
        /* Scale the number of padding points with the wavelength spacing in the
         * reference irradiance spectrum, and the anticipated width of the slit-function
         * (half-width at 1/e in nm).
         * For this purpose, we allow for a relatively large hw1e which will yield
         * a generous amount of padding.
         */
        double hw1e = 0.5;
        num_pad = (sf_ctrl.num_pad_half_widths
                   * (hw1e / ref_irr.delta_wavelength));

        /* We'll allocate the maximum size, but probably use
         * only the relevant subset. */
        num_model_waves = ref_irr.num_wavelen;
     }

   /* Allocate the wavecal structure for this configuration */

   if (NULL == (wct = alloc_wavecal (num_data_waves, num_model_waves,
                                     num_pad, sf_ctrl.num_params)))
     goto error_return;

   wct->is_irradiance = is_irradiance;
   wct->meta = meta;
   wct->irr = ref_irr;       /* struct copy */
   wct->sf_ctrl = sf_ctrl;   /* struct copy */

   if (sf_ctrl.sf_path)
     {
        if (NULL == (wct->sf_table = sf_table_open (sf_ctrl.sf_path, cfg_name)))
          goto error_return;
        if (0)
          {
             if (0 != meta_record_basename (meta, sf_ctrl.sf_path))
               goto error_return;
          }
     }

   if (0 != config_model_components (wct, is_irradiance ? s_irr : s_rad))
     goto error_return;

   if (0 != alloc_term_storage (wct->terms, num_model_waves))
     goto error_return;

   if (0 != config_fit_window (s_irr, &wct->window))
     goto error_return;

   win = &wct->window;
   win->start_pix = fwin.start_pix;
   win->delta_wavelength = fwin.delta_wavelength;
   win->feature_wavelength = fwin.feature_wavelength;

   for (i = 0; i < fwin.num_pix; i++)
     {
        win->pindex[i] = (double)(i + fwin.start_pix);
     }

   if (NULL == (s = config_lookup (cfg, "wavecal_control")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing wavecal_control in param file: %s",
                     __func__, config_error_file (cfg));
        goto error_return;
     }

   if (0 != config_control (s, wct))
     goto error_return;

   if (0 != init_sampling (win, s_band))
     goto error_return;

   return wct;

error_return:
   wavecal_close(wct);
   return NULL;
}

static int evaluate_term (Term_Type *term, size_t num_wave, double *waves,
                          double scale_factor, const double *params)
{
   Refspec_Type *ref = &term->refspec;
   size_t i, offset;
   double *v = term->value;

   memset ((char *)v, 0, 2*num_wave * sizeof(double));

   offset = 0;

   if (term->shapefun)
     {
        Shapefun_Type *st = term->shapefun;
        offset = num_wave;
        if (st->st_eval (st, params, num_wave, waves, v))
          return -1;
        if (st->st_apply_external_scaling)
          {
             for (i = 0; i < num_wave; i++) v[i] *= scale_factor;
          }
        if (term->num_params > 0)
          {
             /* save parameters where term was last evaluated */
             memcpy ((char *)term->eval_params, (char *)params, term->num_params * sizeof(double));
          }
     }

   if (ref->interp)
     {
        Interp_Type *it = ref->interp;
        if (it->it_interp_eval (it, num_wave, waves, v + offset) < 0)
          return -1;
        if (offset)
          {
             for (i = 0; i < num_wave; i++) v[i] *= v[i+offset];
          }
     }

   return 0;
}

int wavecal_query_term (const Wavecal_Type *wct, int nth,
                        Wavecal_Term_Info_Type *info)
{
   const Term_Type *t;
   const Window_Type *win;
   int i;

   if ((wct == NULL) || (info == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: got NULL pointer", __func__);
        return -1;
     }

   win = &wct->window;

   i = 0;
   for (t = wct->terms; t != NULL; t = t->next)
     {
        if (i == nth)
          {
             info->name = t->term_name;
             info->value = t->value;
             info->num_values = win->num_wave;
             return i+1;
          }
        i++;
     }

   return 0;
}

static int combine_terms (Wavecal_Type *wct, double *i0, size_t num_waves, double *model)
{
   int count[NUM_TERM_TYPES];
   size_t num_term_types = NUM_TERM_TYPES;
   size_t i;
   double *ad1, *ad2, *lbe, *bl;
   Term_Type *t;

   /* sum over terms of each type */
   zero_term_sums (wct->term_sums, num_term_types, num_waves);
   memset ((char *)count, 0, num_term_types * sizeof(int));

   for (t = wct->terms; t != NULL; t = t->next)
     {
        double *term_value = t->value;
        double *sum = wct->term_sums[t->term_type];
        for (i = 0; i < num_waves; i++)
          {
             sum[i] += term_value[i];
          }
        count[t->term_type] += 1;
     }

   /* model = sc*(i0 + ad1)*exp(-lbe) + bl + ad2 */

   ad1 = wct->term_sums[TERM_TYPE_AD1];
   ad2 = wct->term_sums[TERM_TYPE_AD2];
   lbe = wct->term_sums[TERM_TYPE_LBE];
   bl  = wct->term_sums[TERM_TYPE_BL];

   if (count[TERM_TYPE_SC] > 0)
     {
        double *sc  = wct->term_sums[TERM_TYPE_SC];
        for (i = 0; i < num_waves; i++)
          {
             model[i] = (sc[i] * (i0[i] + ad1[i]) * exp(-lbe[i])
                         + bl[i] + ad2[i]);
          }
     }
   else
     {
        for (i = 0; i < num_waves; i++)
          {
             model[i] = ((i0[i] + ad1[i]) * exp(-lbe[i])
                         + bl[i] + ad2[i]);
          }
     }

   return 0;
}

static int diagnostic_file_group_id_for_term (int grp, int term_type, const char *term_name)
{
   int len, grp_term;
   const char *type_name;
#define TERM_PATH_BUFSIZE 1024
   char path[TERM_PATH_BUFSIZE];

   switch (term_type)
     {
      case TERM_TYPE_AD1: type_name = "ad1"; break;
      case TERM_TYPE_AD2: type_name = "ad2"; break;
      case TERM_TYPE_LBE: type_name = "lbe"; break;
      case TERM_TYPE_SC:  type_name = "sc"; break;
      case TERM_TYPE_BL:  type_name = "bl"; break;
      default:
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: unknown term type id: %d",
                     __func__, term_type);
        return -1;
        break;
     }

   len = snprintf (path, TERM_PATH_BUFSIZE, "terms/%s/%s", type_name, term_name);
   if (len < 0 || len >= TERM_PATH_BUFSIZE)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: group path name exceeds buffer size (%d): term_name=%s",
                     __func__, TERM_PATH_BUFSIZE, term_name);
        return -1;
     }

   if (0 != TIO_def_grp (grp, path, &grp_term))
     return -1;

   return grp_term;
}

int wavecal_def_term_vars (const Wavecal_Type *wct, int grp,
                           int dimid_step, int dimid_xtrack)
{
   Term_Type *t;
   int dimids[3], dimid_params, dimid_waves;

   dimids[0] = dimid_step;
   dimids[1] = dimid_xtrack;

   for (t = wct->terms; t != NULL; t = t->next)
     {
        int grp_term, num_params, num_waves, varid;

        if ((grp_term = diagnostic_file_group_id_for_term (grp, t->term_type, t->term_name)) < 0)
          return -1;

        num_params = t->num_params;
        if (0 != TIO_def_dim (grp_term, "params", num_params, &dimid_params))
          return -1;

        num_waves = t->num_values;
        if (0 != TIO_def_dim (grp_term, "waves", num_waves, &dimid_waves))
          return -1;

        dimids[2] = dimid_params;
        if (0 != TIO_def_var (grp_term, "params", TIO_DOUBLE, 3, dimids, &varid))
          return -1;
        if (0 != TIO_def_var (grp_term, "params_error", TIO_DOUBLE, 3, dimids, &varid))
          return -1;

        dimids[2] = dimid_waves;
        if (0 != TIO_def_var (grp_term, "value", TIO_DOUBLE, 3, dimids, &varid))
          return -1;
     }

   return 0;
}

int wavecal_write_term_vars (const Wavecal_Type *wct, int grp,
                             int beg_step, int step, int beg_xtrack, int xtrack)
{
   const Window_Type *win = &wct->window;
   Term_Type *t;
   double *xerror = NULL;
   int start[3], count[3];

   start[0] = step - beg_step;
   start[1] = xtrack - beg_xtrack;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;

   /* offset past wavelength grid parameters */
   if (win->xerror)
     {
        xerror = win->xerror + win->num_wave_params;
     }

   for (t = wct->terms; t != NULL; t = t->next)
     {
        int grp_term;

        if ((grp_term = diagnostic_file_group_id_for_term (grp, t->term_type, t->term_name)) < 0)
          return -1;

        count[2] = t->num_params;
        if (0 != TIO_put_var_section (grp_term, "params", start, count, TIO_DOUBLE, t->eval_params))
          return -1;

        if (xerror)
          {
             if (0 != TIO_put_var_section (grp_term, "params_error", start, count, TIO_DOUBLE, xerror))
               return -1;
             xerror += t->num_params;
          }

        count[2] = t->num_values;
        if (0 != TIO_put_var_section (grp_term, "value", start, count, TIO_DOUBLE, t->value))
          return -1;
     }

   return 0;
}

static int monotonic_increasing (const double *x, size_t n)
{
   size_t i;

   for (i = 1; i < n; i++)
     {
        if (x[i-1] >= x[i])
          return 0;
     }

   return 1;
}

static int compute_wavelengths (const Window_Type *win, const double *par, size_t num_wave, double *pindex, double *waves)
{
   Shapefun_Type *wl = win->wavegrid_shapefun;

   if (wl->st_eval (wl, par, num_wave, pindex, waves) < 0)
     return -1;

   if (win->adjust_nominal_wavelength)
     {
        size_t i;
        if (cspline_eval (win->waves_nominal_cspline, num_wave, pindex, win->waves_nominal_scratch))
          return -1;
        for (i = 0; i < num_wave; i++)
          {
             waves[i] += win->waves_nominal_scratch[i];
          }
     }

   return 0;
}

static int pre_convolved_forward_model (Wavecal_Type *wct, const double *params, double *model, double **derivs)
{
   Window_Type *win = &wct->window;
   Window_Sample_Type *samp = &win->samp;
   Reference_Irr_Type *irr = &wct->irr;
   Term_Type *term;
   const double *par;
   double *waves_eval = NULL;
   int num_waves_eval;

   (void) derivs;

   par = params;

   /* compute wavelength as a function of pixel index */
   if (compute_wavelengths (win, par, win->num_wave, win->pindex, win->wave0) < 0)
     return -1;
   par += win->num_wave_params;

   /* validate grid */
   if (0 == monotonic_increasing (win->wave0, win->num_wave))
     {
        tell_vinfo (2, "%s: non-monotonic wavelength grid", __func__);
        return -1;
     }

   if (samp->sample_model == NULL)
     {
        num_waves_eval = win->num_wave;
        waves_eval = win->wave0;
     }
   else
     {
        if (0 != sample_gather (samp->sample_model, win->wave0, samp->waves_eval_samp, &num_waves_eval))
          return -1;
        waves_eval = samp->waves_eval_samp;
     }

   /* evaluate the reference irradiance on the new wavelength grid */
   if (cspline_eval (irr->cspline, num_waves_eval, waves_eval, wct->irr0))
     return -1;

   /* evaluate all model terms on the new wavelength grid */
   for (term = wct->terms; term != NULL; term = term->next)
     {
        double scale_factor = win->rad_mean_ratio;
        if (evaluate_term (term, num_waves_eval, waves_eval, scale_factor, par) < 0)
          return -1;
        par += term->num_params;
     }

   /* combine terms to construct the updated model spectrum */
   if (0 != combine_terms (wct, wct->irr0, num_waves_eval, model))
     return -1;

   return 0;
}

static int interp_array (const double *xa, const double *ya, size_t na, const double *x, double *y, size_t n)
{
   gsl_interp *interp = NULL;
   gsl_interp_accel *accel = NULL;
   size_t i;
   int nan, status = -1;

   if ((NULL == (accel = gsl_interp_accel_alloc ()))
       || (NULL == (interp = gsl_interp_alloc (gsl_interp_linear, na))))
     goto free_and_return;
   if (0 != gsl_interp_init (interp, xa, ya, na))
     goto free_and_return;

   nan = 0;
   for (i = 0; i < n; i++)
     {
        y[i] = gsl_interp_eval (interp, xa, ya, x[i], accel);
        if (gsl_isnan(y[i])) nan++;
     }
   if (nan)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: interpolation yielded %d NaNs", __func__, nan);
        goto free_and_return;
     }

   status = 0;
free_and_return:
   gsl_interp_free (interp);
   gsl_interp_accel_free (accel);

   return status;
}

static int interpolate_sf_convolved (double *wavelen, int num_wavelen,
                                     Window_Type *win, SF_Convolution_Type *sfct,
                                     double *model, double **derivs)
{
   int i;

   if (0 != interp_array (wavelen, sfct->model_convolved, num_wavelen, win->wave0, model, win->num_wave))
     return -1;

   if (derivs == NULL)
     return 0;

   for (i = 0; i < sfct->num_params; i++)
     {
        if (derivs[i] == NULL)
          continue;
        if (0 != interp_array (wavelen, sfct->derivs_convolved[i], num_wavelen, win->wave0, derivs[i], win->num_wave))
          return -1;
     }

   return 0;
}

typedef struct
{
   double params[SFT_MAX_NUM_PARAMS];
   SF_Table_Type *sf_table;
   const double *waves;
   int num_waves;
   int xtrack;
   int mode;
}
SF_Param_Lookup_Type;

static int get_sf_params (int wave_index, int num_pars, double *pars, double *norm, void *cl)
{
   SF_Param_Lookup_Type *lt = (SF_Param_Lookup_Type *)cl;
   SF_Table_Type *stt = lt->sf_table;
   int status;

   if ((stt != NULL) && (lt->mode == SF_MODE_APPLY))
     {
        /* When we have a lookup table, we can (optionally) pre-compute the norms and
         * then interpolate in wavelength.  If the norms were not pre-computed, this
         * can return norm<0 and that will trigger numerical integration to compute
         * the norm.  The latter approach is more accurate, but is extremely slow.
         */
        status = stt->stt_get_params (stt, lt->xtrack, lt->waves[wave_index], pars, norm);
     }
   else
     {
        memcpy ((char *)pars, (char *)lt->params, num_pars * sizeof(double));
        /* norm < 0 will trigger numerical integration to compute the norm. */
        *norm = -1.0;
        status = 0;
     }

   return status;
}

static int select_ref_irr_subset (const Reference_Irr_Type *irr, double wave_beg, double wave_end, int num_pad,
                                  int *pindex_irr_beg, int *pindex_irr_end)
{
   int index_irr_beg, index_irr_end;

   /* Select the relevant subset of the reference irradiance wavelength grid */
   index_irr_beg = bsearch_d (wave_beg, irr->wavelen, irr->num_wavelen);
   index_irr_end = bsearch_d (wave_end, irr->wavelen, irr->num_wavelen);
   index_irr_end++;

   /* evaluate the irradiance spectrum a bit beyond the required wavelength range */
   index_irr_beg -= num_pad;
   index_irr_end += num_pad;

   /* Make sure the indices are still valid */
   if (index_irr_beg < 0)
     index_irr_beg = 0;
   if ((int) irr->num_wavelen <= index_irr_end)
     index_irr_end = irr->num_wavelen - 1;

   /* Check for the impossible */
   if (index_irr_end < index_irr_beg)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid interval!!: %f @index_irr_beg=%d, %f @index_irr_end=%d",
                     __func__, wave_beg, index_irr_beg, wave_end, index_irr_end);
        return -1;
     }

   *pindex_irr_beg = index_irr_beg;
   *pindex_irr_end = index_irr_end;

   return 0;
}

static int convolve_forward_model (Wavecal_Type *wct, const double *params, double *model, double **derivs)
{
   Window_Type *win = &wct->window;
   Window_Sample_Type *samp = &win->samp;
   Reference_Irr_Type *irr = &wct->irr;
   SF_Convolution_Type *sfct = win->sfct;
   Term_Type *term;
   SF_Param_Lookup_Type sf_lookup = {0};
   const double *par;
   double *irr_wavelen = NULL;
   double *irr_value = NULL;
   double *waves_eval;
   double *model_tmp;
   size_t num_irr_waves;
   int index_irr_beg, index_irr_end;
   int status, index_slit_param0;
   int num_waves_eval;
   int use_derivs=0;

   par = params;

   /* compute wavelength as a function of pixel index
    * (wavelength grid parameters are at the front of the full param array) */
   if (compute_wavelengths (win, params, win->num_wave, win->pindex, win->wave0) < 0)
     return -1;

   /* validate grid */
   if (0 == monotonic_increasing (win->wave0, win->num_wave))
     {
        tell_vinfo (2, "%s: non-monotonic wavelength grid", __func__);
        return -1;
     }

   /* Select the relevant subset of the reference irradiance wavelength grid */
   if (0 != select_ref_irr_subset (irr, win->wave0[0], win->wave0[win->num_wave-1], sfct->num_pad, &index_irr_beg, &index_irr_end))
     return -1;

   /* define irradiance subset pointers */
   num_irr_waves = index_irr_end - index_irr_beg + 1;
   irr_wavelen = irr->wavelen + index_irr_beg;
   irr_value  = irr->irradiance + index_irr_beg;

   /* set the number of wavelengths in the unpadded model */
   sfct->num_waves = num_irr_waves;

   /* skip wavelength grid parameters */
   par += win->num_wave_params;

   if (samp->sample_model == NULL)
     {
        num_waves_eval = num_irr_waves;
        waves_eval = irr_wavelen;
     }
   else
     {
        if ((0 != sample_gather (samp->sample_model, irr_wavelen, samp->waves_eval_samp, &num_waves_eval))
            || (0 != sample_gather (samp->sample_model, irr_value, samp->irr_value_samp, NULL)))
          return -1;
        waves_eval = samp->waves_eval_samp;
     }

   /* Construct the model on the high-resolution reference irradiance wavelength grid */
   for (term = wct->terms; term != NULL; term = term->next)
     {
        double scale_factor = win->rad_mean_ratio;
        if (evaluate_term (term, num_waves_eval, waves_eval, scale_factor, par) < 0)
          return -1;
        par += term->num_params;
     }

   /* Combine terms to construct the updated model spectrum.
    * After the call, the zero-padded array sfct->model_padded looks like this:
    * [<num_pad zeros>|<num_irr_waves model values>|<num_pad zeros><more zeros>]
    */
   if (samp->sample_model == NULL)
     {
        if (0 != combine_terms (wct, irr_value, num_irr_waves, sfct->model_padded + sfct->num_pad))
          return -1;
        model_tmp = model;
     }
   else
     {
        if (0 != combine_terms (wct, samp->irr_value_samp, num_waves_eval, samp->combined_model_samp))
          return -1;
        memset ((char *)sfct->model_padded, 0, sfct->num_alloc * sizeof(double));
        if (0 != sample_scatter (samp->sample_model, samp->combined_model_samp, sfct->model_padded + sfct->num_pad))
          return -1;
        model_tmp = samp->convolved_model;
     }

   /* slit-function parameter lookup needs some context and reference data */
   sf_lookup.mode = wct->sf_ctrl.mode;
   sf_lookup.xtrack = wct->xtrack;
   sf_lookup.sf_table = wct->sf_table;
   sf_lookup.waves = irr_wavelen;
   sf_lookup.num_waves = sfct->num_waves;

   if (wct->sf_ctrl.mode == SF_MODE_FIT)
     {
        /* When fitting the slit function parameters are at the end of the array. */
        if (derivs) use_derivs++;
        index_slit_param0 = par - params;
        memcpy ((char *)sf_lookup.params, (char *)par, wct->sf_ctrl.num_params * sizeof(double));
     }
   else if (wct->sf_ctrl.initial_params)
     {
        memcpy ((char *)sf_lookup.params, (char *)wct->sf_ctrl.initial_params,
                wct->sf_ctrl.num_params * sizeof(double));
     }

   /* Convolve the model with the slit-function. */
   status = sft_apply (wct->sft, get_sf_params, &sf_lookup, sfct->num_waves, sfct->model_padded,
                       sfct->model_convolved,
                       use_derivs ? sfct->derivs_convolved : NULL);
   if (status) return -1;

   /* Interpolate the convolved model and derivatives onto the parametrized wavelength grid
    * (sfct->hr_model_convolved -> model) */
   status = interpolate_sf_convolved (irr_wavelen, num_irr_waves, win, sfct, model_tmp,
                                      use_derivs ? &derivs[index_slit_param0] : NULL);
   if (status) return -1;

   /* note that model_tmp = model when we're not sampling (see above) */
   if (samp->sample_data)
     {
        if (0 != sample_gather (samp->sample_data, model_tmp, model, NULL))
          return -1;
     }

   return 0;
}

static void write_params (FILE *fp, const double *x, int n)
{
   int i;
   fprintf (fp, "params:\n");
   for (i = 0; i < n; i++)
     {
        fprintf (fp, "%15.6e ", x[i]);
        if (0 == ((i+1) % 6)) fprintf (fp, "\n");
     }
   fprintf (fp, "\n");
}

static void write_statistic (FILE *fp, const double *f, int m)
{
   double sumsq;
   int i;

   sumsq = 0.0;
   for (i = 0; i < m; i++)
     {
        sumsq += f[i] * f[i];
     }
   fprintf (fp, "sumsq = %15.6e\n", sumsq);
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
   Window_Type *win = &p->wct->window;
   Window_Sample_Type *samp = &win->samp;
   const double *spec = p->spec;
   const double *weight = p->weight;
   double *model = p->model;
   int i, num_resid, eval_status;

   (void) n;

   if (0) write_params (stderr, x, n);

   num_resid = samp->sample_data ? samp->sample_data->num_samples : m;

   if (p->wct->sf_ctrl.mode == SF_MODE_NONE)
     {
        eval_status = pre_convolved_forward_model (p->wct, x, model, dvec);
     }
   else
     {
        eval_status = convolve_forward_model (p->wct, x, model, dvec);
     }

   if (eval_status)
     {
        /* Model evaluation failed, but returning a negative value
         * would halt the optimizer. Zero the model to penalize
         * this parameter set and keep going. */
        memset ((char *)model, 0, num_resid * sizeof(double));
     }

   for (i = 0; i < num_resid; i++)
     {
        fvec[i] = (model[i] - spec[i]) * weight[i];
     }

   if (dvec)
     {
        int k;
        for (k = 0; k < n; k++)
          {
             double *dvec_k = dvec[k];
             if (dvec_k == NULL)
               continue;
             if (eval_status == 0)
               {
                  for (i = 0; i < num_resid; i++)
                    {
                       dvec_k[i] *= weight[i];
                    }
               }
             else
               {
                  /* Handle failure of the model evaluation */
                  memset ((char *)dvec_k, 0, num_resid * sizeof(double));
               }
          }
     }

   if (0) write_statistic (stderr, fvec, num_resid);

   return 0;
}

static int compute_rad_mean_ratio (Wavecal_Type *wct)
{
   Window_Type *win = &wct->window;
   Reference_Irr_Type *irr = &wct->irr;
   double sum_irr, sum_rad;
   size_t i;

   /* compute wavelength as a function of pixel index */
   if (compute_wavelengths (win, win->wave_params, win->num_wave, win->pindex, win->wave0) < 0)
     return -1;

   /* evaluate the reference irradiance on the target wavelength grid */
   if (cspline_eval (irr->cspline, win->num_wave, win->wave0, wct->irr0))
     return -1;

   sum_irr = 0.0;
   sum_rad = 0.0;
   for (i = 0; i < win->num_wave; i++)
     {
        sum_irr += wct->irr0[i];
        sum_rad += win->spec_scaled[i];
     }

   if (sum_irr == 0.0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: irradiance sum is zero!!", __func__);
        return -1;
     }

   /* rad_mean_ratio = (mean_radiance)/(mean_irradiance)
    *                = (sum_rad[*]/num_rad) / (sum_irr[*]/num_irr)
    *                = sum_rad[*]/sum_irr[*]
    * because the wavelength grids are identical.
    */
   win->rad_mean_ratio = sum_rad / sum_irr;

   return 0;
}

/* When the wavelength grid is defined by Chebyshev polynomial
 * coefficients, the coefficients c[0] and c[1] correspond to:
 *
 *   c[0] = the midpoint wavelength
 *        = \lambda( (k0+k1)/2 )
 *   c[1] = the half-width of the calibration band
 *        = 0.5 * (\lambda(k1) - \lambda(k0))
 *
 * where the calibration band is:
 *    \lambda(k0) <= \lambda <= \lambda(k1)
 * and k0, k1 are pixel indices.
 */
static void estimate_midpoint_wavelength (Window_Type *win, const double *spec,
                                          double *wavelength_of_window_midpoint)
{
   Shapefun_Type *st = win->wavegrid_shapefun;
   int num_spec = win->num_wave;
   double bin_width = win->delta_wavelength;
   double wavelength_of_window_minimum = win->feature_wavelength;
   int i, index_of_window_minimum = 0;
   double window_minimum = spec[0];

   /* This estimate assumes the wavelength scale is
    * parameterized by a Chebyshev polynomial.
    * Don't use it for anything else.
    */
   if (st->st_method (st) != SHAPEFUN_TYPE_CHEB)
     return;

   /* If required parameters haven't been provided, do nothing */
   if ((bin_width <= 0.0)
       || (wavelength_of_window_minimum <= 0.0))
     return;

   for (i = 1; i < num_spec; i++)
     {
        if (0 < spec[i] && spec[i] < window_minimum)
          {
             window_minimum = spec[i];
             index_of_window_minimum = i;
          }
     }

   *wavelength_of_window_midpoint =
     (wavelength_of_window_minimum
       + bin_width * (0.5*(num_spec-1) - index_of_window_minimum));
}

int wavecal_get_initial_params (Wavecal_Type *wct, const double *p_wave,
                                double *wave_params)
{
   Window_Type *win = &wct->window;
   const double *wave = p_wave + win->start_pix;

   if (0 != init_window_shapefun (wct, wave))
     return -1;

   memcpy ((char *)wave_params, (char *)win->wave_params,
           win->num_wave_params * sizeof(double));

   return 0;
}

int wavecal_adjust (const Wavecal_Type *wct, const Wadj_Type *wadj, int xtrack,
                    const double *narrow_band_wave_params, double *final_wave_params)
{
   Wadj_Cheb_Series_Type attr_narrow;
   double narrow_band_mid_wl_shift = 0.0;

   /* The calling routine should make sure that the various
    * parameter array lengths are correct */
   if (0 != wadj_narrow_band_get_attr (wadj, &attr_narrow))
     return -1;

   if (attr_narrow.num_series_coeff > 0)
     {
        const Window_Type *win = &wct->window;
        const double *tbl_narrow_band_nwave_params;
        double mid_pix, mid_wl_fit, mid_wl_tbl;

        mid_pix = 0.5*(attr_narrow.pix_min + attr_narrow.pix_max);

        if (NULL == (tbl_narrow_band_nwave_params = wadj_narrow_band_coeff (wadj, xtrack)))
          return -1;

        if ((compute_wavelengths (win, narrow_band_wave_params, 1, &mid_pix, &mid_wl_fit) < 0)
            ||(compute_wavelengths (win, tbl_narrow_band_nwave_params, 1, &mid_pix, &mid_wl_tbl) < 0))
          return -1;

        narrow_band_mid_wl_shift = mid_wl_fit - mid_wl_tbl;
     }

   return wadj_final_coeff (wadj, xtrack, narrow_band_mid_wl_shift, final_wave_params);
}

static int init_slit_function (Wavecal_Type *wct, int xtrack, size_t num,
                               struct mp_par_struct **p_param_ctrl)
{
   struct mp_par_struct *param_ctrl = NULL;
   int k, k0;

   if (wct->sft == NULL)
     return 0;

   if (0 != sft_config (wct->sft, asg_normed_plus_derivs, wct->irr.delta_wavelength,
                        wct->sf_ctrl.param_step))
     {
        return -1;
     }

   /* slit-function lookup tables may require this */
   wct->xtrack = xtrack;

   if (wct->sf_ctrl.mode != SF_MODE_FIT)
     return 0;

   /* When fitting the slit function, the objective function
    * may compute parameter derivatives for the slit function
    * parameters. */
   if (NULL == (param_ctrl = (struct mp_par_struct *)MALLOC (num * sizeof (*param_ctrl))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)param_ctrl, 0, num * sizeof(*param_ctrl));

   /* Assume the slit function parameters are at the end of the full parameter array. */
   k0 = num - wct->sf_ctrl.num_params;
#if 0
   /* freeze the asymmetry parameter at zero (FIXME - control this from config file?) */
   param_ctrl[k0 + 2].fixed = 1;
#endif
   for (k = 0; k < wct->sf_ctrl.num_params; k++)
     {
        param_ctrl[k0 + k].side = 3;
        /* param_ctrl[k0 + k].deriv_debug = 1; */
     }

   *p_param_ctrl = param_ctrl;

   return 0;
}

int wavecal_fit (Wavecal_Type *wct, int xtrack,
                 const double *p_wave, const double *p_spec, const double *p_specerr,
                 const unsigned int *p_pixel_quality_flag, const Wavecal_Config_Type *config,
                 double *wave_params, Wavecal_Result_Type *result)
{
   Mpfit_Interface_Type mp = {0};
   struct mp_config_struct fit_config = {0};
   struct mp_result_struct fit_result = {0};
   struct mp_par_struct *param_ctrl = NULL;
   Fit_Control_Type *fit_ctrl = &wct->fit_ctrl;
   Window_Type *win = &wct->window;
   Window_Sample_Type *samp = &win->samp;
   double fill_value = config->fill_value;
   const double *wave = p_wave + win->start_pix;
   const double *spec = p_spec + win->start_pix;
   const double *specerr = p_specerr + win->start_pix;
   const unsigned int *pqf = p_pixel_quality_flag + win->start_pix;
   double *spec_scaled = win->spec_scaled;
   double *weight = win->weight;
   double *win_spec_scaled = NULL;
   double *win_weight = NULL;
   double *params = NULL;
   double scale_factor;
   size_t i, num;
   int status = WAVECAL_FIT_ERROR;
   int num_residuals, num_params, mp_status, num_good;

   if (0 != init_window_shapefun (wct, wave))
     return WAVECAL_FIT_ERROR;

   if (win->adjust_nominal_wavelength)
     {
        cspline_free (win->waves_nominal_cspline); /* FIXME - preallocate this */
        if (NULL == (win->waves_nominal_cspline = cspline_interpol (win->num_wave, win->pindex, wave, NULL)))
          return WAVECAL_FIT_ERROR;
     }

   /* If the fit fails, the initial wavelength parameter guess
    * will serve as the result */
   memcpy ((char *)wave_params, (char *)win->wave_params,
           win->num_wave_params * sizeof(double));

   if (0 != init_window_reference_spectra (wct, xtrack))
     return WAVECAL_FIT_ERROR;

   num_good = 0;
   /* Scale the radiance and irradiance by the same factor */
   scale_factor = wct->irr.file.scale_factor;
   for (i = 0; i < win->num_wave; i++)
     {
        double spec_i = spec[i];
        double err_i = specerr[i];

        if ((pqf[i] == 0) && (0 != isfinite(spec_i)) && (spec_i != fill_value))
          {
             spec_scaled[i] = spec_i * scale_factor;
             num_good++;
          }
        else spec_scaled[i] = 0.0;

        if ((pqf[i] == 0) && (0 != isfinite(err_i)) && (err_i != fill_value) && (err_i != 0.0)
            && (spec_i > 0))
          {
             weight[i] = 1.0/(fabs(err_i) * scale_factor);
          }
        else weight[i] = 0.0;

        // adding by WHou
        if (win->num_wave == 1028)
          {
             if ( (i < 10) || (i > 1017) )
               {
                  weight[i] = weight[i] / 100.0;
               }
          }
        // end adding
     }

   /* Too many bad pixels? */
   if (num_good < MIN_ACCEPTABLE_GOOD_PIXEL_FRACTION * win->num_wave)
     return WAVECAL_FIT_BAD;

   if (collect_params (wct, &num, &params) < 0)
     goto return_error;

   estimate_midpoint_wavelength (win, spec_scaled, &params[0]);

   if (wct->is_irradiance == 0)
     {
        if (0 != compute_rad_mean_ratio (wct))
          goto return_error;
     }

   /* If we're using the slit-function, initialize it here */
   if (0 != init_slit_function (wct, xtrack, num, &param_ctrl))
     goto return_error;

   if (samp->sample_data == NULL)
     {
        /* fit all data pixels in the specified range */
        win_spec_scaled = win->spec_scaled;
        win_weight = win->weight;
     }
   else
     {
        /* select a sample of data pixels to fit */
        if ((0 != sample_gather (samp->sample_data, win->spec_scaled, samp->spec_scaled_samp, NULL))
            || (0 != sample_gather (samp->sample_data, win->weight, samp->weight_samp, NULL))
            || (0 != sample_gather (samp->sample_data, win->wave0, samp->wave0_samp, NULL))
           )
          goto return_error;
        win_spec_scaled = samp->spec_scaled_samp;
        win_weight = samp->weight_samp;

        /* the slit function determines which model wavelengths
         * contribute to the sampled data pixels */
        if (0 != init_model_sampling (win, wct->sf_ctrl.mode, wave, &wct->irr))
          goto return_error;
     }

   mp.wct = wct;
   mp.spec = win_spec_scaled;
   mp.weight = win_weight;
   mp.model = win->model;
   mp.counter = 0;

   num_params = num;
   num_residuals = win->num_wave;

   if (NULL == (win->xerror = alloc_doubles (num_params)))
     goto return_error;
   memset ((char *)win->xerror, 0, num_params * sizeof(double));

   fit_config.xtol = fit_ctrl->xtol;
   fit_config.ftol = fit_ctrl->ftol;
   fit_config.maxiter = fit_ctrl->maxiter;
   fit_config.maxfev = (fit_ctrl->maxfev ?
                        fit_ctrl->maxfev : fit_ctrl->maxiter * num_params);

   fit_result.resid = win->residuals; /* FIXME: make this a debug option? */
   fit_result.xerror = win->xerror;

   if (0) write_params (stderr, params, num_params);

   mp_status = mpfit (&mpfit_objective_function, num_residuals,
                      num_params, params, param_ctrl, &fit_config, &mp, &fit_result);

   /* making sure the model is evaluated at best-fit */
   if (0 != mpfit_objective_function (num_residuals, num_params,
                                      params, win->residuals, NULL, &mp))
     goto return_error;

   if (result)
     {
        /* Warning: pointers are valid only until wavecal_close gets called */
        memcpy ((char *)win->wave_params, (char *)params,
                win->num_wave_params * sizeof(double));
        result->wave_params = win->wave_params;
        result->wave_params_error = win->xerror;
        result->num_wave_params = win->num_wave_params;
        if (wct->sf_ctrl.mode == SF_MODE_FIT)
          {
             SF_Convolution_Type *sfct = win->sfct;
             double *sf_params = params + (num_params - sfct->num_params);
             double *sf_params_error = win->xerror + (num_params - sfct->num_params);
             memcpy ((char *)sfct->params, (char *)sf_params, sfct->num_params * sizeof(double));
             memcpy ((char *)sfct->params_error, (char *)sf_params_error, sfct->num_params * sizeof(double));
             result->sf_params = sfct->params;
             result->sf_params_error = sfct->params_error;
             result->num_sf_params = sfct->num_params;
          }
        else
          {
             result->sf_params = NULL;
             result->num_sf_params = 0;
          }
        result->wave = samp->sample_data ? samp->wave0_samp : win->wave0;
        result->model = win->model;
        result->spec_scaled = win_spec_scaled;
        result->weight = win_weight;
        result->residuals = win->residuals;
        result->bestnorm = fit_result.bestnorm;
        result->num_fit = samp->sample_data ? (size_t) samp->sample_data->num_samples : win->num_wave;
        result->start_pix = win->start_pix;
        result->nfev = fit_result.nfev;
        result->niter = fit_result.niter;
        result->opt_status = mp_status;
        if (0) write_params (stderr, params, num_params);
        if (0) write_statistic (stderr, win->residuals, num_residuals);
     }

   status = ((0 < mp_status) && (mp_status < 4)) ? WAVECAL_FIT_GOOD : WAVECAL_FIT_BAD;

   if (status == WAVECAL_FIT_GOOD)
     {
        memcpy ((char *)wave_params, (char *)params,
                win->num_wave_params * sizeof(double));
     }

return_error:
   FREE(params);
   FREE(param_ctrl);

   return status;
}
