#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>
#include "config.h"
#include "granule.h"

#include "util.h"

#define MAX_ENABLE_STATE_SETTINGS  100

typedef struct
{
   /* These const char * pointers are assumed to point to either
    * a string literal, or to memory managed by libconfig
    */
   const char *name;
   const char *enum_state;
   int boolean_state;
   int state_type;
}
Enable_State_Type;

static Enable_State_Type Enable_Settings[MAX_ENABLE_STATE_SETTINGS];
static int Num_Enable_Settings;

static int find_enable_setting (const char *name)
{
   int i;

   for (i = 0; i < Num_Enable_Settings; i++)
     {
        Enable_State_Type *s = &Enable_Settings[i];
        if (0 == strcmp (s->name, name))
          {
             return i;
             break;
          }
     }

   return -1;
}

static int enable_state_new (const char *name, int type, int boolean_state, const char *enum_state)
{
   Enable_State_Type *s;
   int n;

   if (Num_Enable_Settings == MAX_ENABLE_STATE_SETTINGS)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: enable state table is full", __func__);
        return -1;
     }

   if ((n = find_enable_setting (name)) < 0)
     {
        n = Num_Enable_Settings;
        Num_Enable_Settings++;
     }

   s = &Enable_Settings[n];
   s->state_type = type;
   s->boolean_state = boolean_state;
   s->enum_state = enum_state;
   s->name = name;

   return 0;
}

int enable_state_define (config_t *cfg, const char *state_name)
{
   config_setting_t *s, *m;
   const char *enum_state = NULL;
   int boolean_state = 0;
   int type;

   if ((NULL == (s = config_lookup (cfg, "enable")))
       || (NULL == (m = config_setting_get_member (s, state_name))))
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: reading %s from group 'enable' in %s",
                     __func__, state_name, config_error_file(cfg));
        return -1;
     }

   type = config_setting_type (m);

   switch (type)
     {
      case CONFIG_TYPE_BOOL:
        boolean_state = config_setting_get_bool (m);
        tell_vlog (TELL_MSGTYPE_INFO, 1, "SELECT %s: %s", state_name, boolean_state ? "ON" : "OFF");
        break;

      case CONFIG_TYPE_STRING:
        enum_state = config_setting_get_string (m);
        tell_vlog (TELL_MSGTYPE_INFO, 1, "SELECT %s: %s", state_name, enum_state ? enum_state : "NULL");
        break;

      default:
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: invalid setting: %s",  __func__, state_name);
        return -1;
     }

   return enable_state_new (state_name, type, boolean_state, enum_state);
}

const char *enable_state_query_enum (const char *name)
{
   Enable_State_Type *s;
   int n;

   if ((n = find_enable_setting (name)) < 0)
     return NULL;

   s = &Enable_Settings[n];

   if (s->state_type != CONFIG_TYPE_STRING)
     return NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: %s", name, s->enum_state);

   return s->enum_state;
}

int enable_state_query_bool (const char *name)
{
   Enable_State_Type *s;
   int n;

   if ((n = find_enable_setting (name)) < 0)
     return -1;

   s = &Enable_Settings[n];

   if (s->state_type != CONFIG_TYPE_BOOL)
     return -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: %s", name, s->boolean_state ? "ON" : "OFF");

   return s->boolean_state;
}

int bsearch_d (double t, const double *x, int n)
{
   int n0, n1, n2;
   double xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

int bsearch_f (float t, const float *x, int n)
{
   int n0, n1, n2;
   float xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

double *alloc_doubles (int n)
{
   double *a = NULL;

   if (NULL == (a = (double *)MALLOC (n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   return a;
}

int find_x (double x, const double *a, int na)
{
   if (x < a[0] || x >= a[na-1])
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: value=%g outside array bounds [%g,%g)",
                     __func__, x, a[0], a[na-1]);
        return -1;
     }
   return bsearch_d (x, a, na);
}

/* return 0 means absent, -1 means error, 1 means read successfully */
int read_config_float_array (config_setting_t *s, const char *name,
                             double **pa, size_t *pnum_a)
{
   config_setting_t *fs;
   double *a;
   unsigned int i, na;
   int len;

   *pa = NULL;
   if (pnum_a) *pnum_a = 0;

   /* absent/empty is ok */
   if ((NULL == (fs = config_setting_get_member (s, name)))
       || (0 == (len = config_setting_length (fs))))
     return 0;

   if (NULL == (a = (double *)MALLOC (len * sizeof (double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   na = len;

   for (i = 0; i < na; i++)
     {
        a[i] = config_setting_get_float_elem (fs, i);
     }

   *pa = a;
   if (pnum_a) *pnum_a = na;

   return 1;
}

char *expand_string (const char *s)
{
   wordexp_t we;
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof(wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding string: %s", __func__, s);
        goto return_error;
     }

   if (NULL == (s_exp = strdup (we.we_wordv[0])))
     {
        tell_verror (TELL_MALLOC_ERROR,
                     "%s: strdup failed", __func__);
     }

return_error:
   wordfree (&we);
   return s_exp;
}

char *path_concat (const char *dir, const char *basename)
{
   int status, len;
   char *s;

   len = strlen (dir) + strlen(basename) + 2;
   if (NULL == (s = MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   status = snprintf (s, len, "%s/%s", dir, basename);
   if ((status < 0) || (status >= len))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed", __func__);
        FREE(s);
        return NULL;
     }

   return s;
}

int meta_record_basename (TIO_Meta_Type *meta, const char *path)
{
   const char *path_basename;

   if (path == NULL)
     return 0;

   if (NULL != (path_basename = strrchr (path, '/')))
     {
        path_basename++;
     }
   else path_basename = path;

   return tio_meta_append_string (meta, "input_files", path_basename);
}

/* Trend parameter collection */

#define NUM_QUAD     4
#define NUM_OCTANTS  8
#define NUM_PQF_BITS  ((int) (8*sizeof(Image_Pqf_Bitmap_Type)))

struct Trend_Record_Type
{
   const Trend_File_Type *tft;
   int index;
   int exposure_type;
   double start_time;
   float eoffsets[NUM_OCTANTS];
   int phase_change[NUM_QUAD];
   float gain[NUM_OCTANTS];
   float readnoise[NUM_OCTANTS];
   float fpe_temp;
   float fpa_temp;
   int num_dg_rows;
   int num_tg_rows;
   float storage_region_dark[NUM_QUAD];
   float mean_dc[NUM_QUAD];
   float stddev_dc[NUM_QUAD];
   int pqf_bits[NUM_PQF_BITS];
   int pqf_bits_uv[NUM_PQF_BITS];
   int pqf_bits_vis[NUM_PQF_BITS];
   double solar_theta;
   double solar_phi;
   int use_reference_diffuser;
};

struct Trend_File_Type
{
   int ncid;
   int exposure_type;
   Trend_Record_Type *tr;
};

static Trend_Record_Type *Active_Record = NULL;
static int Have_Trend_File = 0;

Trend_Record_Type *trend_collect_set_active_record (Trend_Record_Type *tr)
{
   Trend_Record_Type *old_tr = Active_Record;
   Active_Record = tr;
   return old_tr;
}

typedef struct
{
   const char *name;
   const char *text;
}
Text_Attr_Type;

static int define_text_attrs (int grp, int varid, const Text_Attr_Type *attrs)
{
   const Text_Attr_Type *a;

   for (a = attrs; a->name != NULL; a++)
     {
        size_t len = strlen(a->text) + 1;
        if (0 != TIO_put_att (grp, varid, a->name, TIO_CHAR, len, a->text))
          return -1;
     }

   return 0;
}

static int enable_var_deflation (int grp, int varid)
{
   int shuffle = 1, deflate = 1, deflate_level = 1;
   return TIO_def_var_deflate (grp, varid, shuffle, deflate, deflate_level);
}

static void trend_type_free (Trend_File_Type *tft)
{
   if (tft == NULL)
     return;
   FREE(tft);
}

static Trend_File_Type *trend_type_alloc (int exposure_type)
{
   Trend_File_Type *tft = NULL;

   if (NULL == (tft = (Trend_File_Type *)MALLOC (sizeof *tft)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)tft, 0, sizeof(*tft));
   tft->exposure_type = exposure_type;
   tft->tr = NULL;

   return tft;
}

int trend_collect_close (Trend_File_Type *tft)
{
   int status = 0;

   if (tft)
     {
        if (tft->ncid != 0)
          {
             status = TIO_close (tft->ncid);
          }
        trend_type_free (tft);
     }

   return status;
}

Trend_File_Type *trend_collect_open (const char *trend_file, int exposure_type)
{
   Trend_File_Type *tft = NULL;
   const char *product_type = NULL;
   int dimid_time, dimid_quad, dimid_oct, dimid_bit, varid;
   int dimid_time_quad[2];
   int dimid_time_oct[2];
   int dimid_time_bit[2];
   const Text_Attr_Type global_attrs[] =
     {
        {"comment",
"Dimension quad = CCD quadrants; quad={0,1,2,3} => {A,B,C,D}\n" \
"Dimension oct = CCD octants; ADC {odd,even} in each quadrant => oct={Ao,Bo,Co,Do, Ae,Be,Ce,De};\n" \
"Same storage order as (adc,quad) arrays in calibration key data file."},
        {NULL, NULL}
     };
   const Text_Attr_Type time_attrs[] =
     {
        {"long_name", "exposure start time"},
        {NULL, NULL}
     };
   const Text_Attr_Type eoffset_attrs[] =
     {
        {"long_name", "electronic offset"},
        {NULL, NULL}
     };
   const Text_Attr_Type phase_change_attrs[] =
     {
        {"long_name", "quadrant parity phase change"},
        {"comment", "0:quadrant even/odd column ADC configured as in FPS testing, 1: the ADC configuration is swapped"},
        {NULL, NULL}
     };
   const Text_Attr_Type gain_attrs[] =
     {
        {"long_name", "gain"},
        {"units", "electrons/DN"},
        {NULL, NULL}
     };
   const Text_Attr_Type fpa_temp_attrs[] =
     {
        {"long_name", "focal plane array temperature"},
        {"units", "C"},
        {NULL, NULL}
     };
   const Text_Attr_Type fpe_temp_attrs[] =
     {
        {"long_name", "focal plane electronics temperature"},
        {"units", "C"},
        {NULL, NULL}
     };
   const Text_Attr_Type readnoise_attrs[] =
     {
        {"long_name", "read-out noise"},
        {"units", "electrons"},
        {NULL, NULL}
     };
   const Text_Attr_Type num_dg_rows_attrs[] =
     {
        {"long_name", "number of rows offset from readout for storage region dark sum"},
        {NULL, NULL}
     };
   const Text_Attr_Type num_tg_rows_attrs[] =
     {
        {"long_name", "number of rows in storage region dark sum"},
        {NULL, NULL}
     };
   const Text_Attr_Type sdc_attrs[] =
     {
        {"long_name", "quadrant storage region dark current mean"},
        {"units", "electrons/s"},
        {NULL, NULL}
     };
   const Text_Attr_Type mean_dc_attrs[] =
     {
        {"long_name", "quadrant dark current mean"},
        {"units", "electrons/s"},
        {NULL, NULL}
     };
   const Text_Attr_Type stddev_dc_attrs[] =
     {
        {"long_name", "quadrant dark current standard deviation"},
        {"units", "electrons/s"},
        {NULL, NULL}
     };
   const Text_Attr_Type solar_theta_attrs[] =
     {
        {"long_name", "solar boresight angle"},
        {"units", "degrees"},
        {NULL, NULL}
     };
   const Text_Attr_Type solar_phi_attrs[] =
     {
        {"long_name", "solar azimuthal angle"},
        {"units", "degrees"},
        {NULL, NULL}
     };
   const Text_Attr_Type ref_diffuser_attrs[] =
     {
        {"long_name", "reference_diffuser"},
        {"comment", "0:working diffuser, 1:reference diffuser"},
        {NULL, NULL}
     };
   const Text_Attr_Type pqf_bits_attrs[] =
     {
        {"long_name", "count of pixels with each pixel_quality_flag bit set"},
        {NULL, NULL}
     };
   const Text_Attr_Type pqf_bits_uv_attrs[] =
     {
        {"long_name", "count of UV pixels with each pixel_quality_flag bit set"},
        {NULL, NULL}
     };
   const Text_Attr_Type pqf_bits_vis_attrs[] =
     {
        {"long_name", "count of VIS pixels with each pixel_quality_flag bit set"},
        {NULL, NULL}
     };
   int status = -1;

   if (trend_file == NULL)
     return NULL;

   if (NULL == (tft = trend_type_alloc (exposure_type)))
     return NULL;

   if (0 != TIO_create (trend_file, NC_NETCDF4, &tft->ncid))
     goto return_status;
   if (0 != tio_write_epoch_timestamp (tft->ncid, NC_GLOBAL))
     goto return_status;

   if (EXPREC_TYPE_IS_DARK(exposure_type))
     {
        product_type = "TREND_DRK";
     }
   else if (EXPREC_TYPE_IS_IRRADIANCE(exposure_type))
     {
        product_type = "TREND_IRR";
     }
   else
     {
        product_type = "TREND_RAD";
     }
   if (0 != TIO_label_product (tft->ncid, product_type, 0, 0))
     goto return_status;

   if (0 != define_text_attrs (tft->ncid, NC_GLOBAL, global_attrs))
     goto return_status;

   if ((0 != TIO_def_dim (tft->ncid, TEMPO_VAR_TIME, NC_UNLIMITED, &dimid_time))
       || (0 != TIO_def_dim (tft->ncid, "quad", 4, &dimid_quad))
       || (0 != TIO_def_dim (tft->ncid, "oct", 8, &dimid_oct))
       || (0 != TIO_def_dim (tft->ncid, "bit", NUM_PQF_BITS, &dimid_bit))
      )
     goto return_status;

   dimid_time_quad[0] = dimid_time;
   dimid_time_quad[1] = dimid_quad;

   dimid_time_oct[0] = dimid_time;
   dimid_time_oct[1] = dimid_oct;

   dimid_time_bit[0] = dimid_time;
   dimid_time_bit[1] = dimid_bit;

   if ((0 != TIO_def_var (tft->ncid, TEMPO_VAR_TIME, NC_DOUBLE, 1, &dimid_time, &varid))
       || (0 != define_text_attrs (tft->ncid, varid, time_attrs))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != tio_write_timestamp_unit_string (tft->ncid, TEMPO_VAR_TIME)))
     goto return_status;

   if ((0 != TIO_def_var (tft->ncid, "eoffsets", NC_FLOAT, 2, dimid_time_oct, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, eoffset_attrs)))
     goto return_status;
   if ((0 != TIO_def_var (tft->ncid, "phase_change", NC_INT, 2, dimid_time_quad, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, phase_change_attrs)))
     goto return_status;

   if ((0 != TIO_def_var (tft->ncid, "gain", NC_FLOAT, 2, dimid_time_oct, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, gain_attrs)))
     goto return_status;
   if ((0 != TIO_def_var (tft->ncid, "fpa_temp", NC_FLOAT, 1, &dimid_time, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, fpa_temp_attrs)))
     goto return_status;
   if ((0 != TIO_def_var (tft->ncid, "fpe_temp", NC_FLOAT, 1, &dimid_time, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, fpe_temp_attrs)))
     goto return_status;

   if ((0 != TIO_def_var (tft->ncid, "readout_noise", NC_FLOAT, 2, dimid_time_oct, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, readnoise_attrs)))
     goto return_status;

   if ((0 != TIO_def_var (tft->ncid, "num_dg_rows", NC_INT, 1, &dimid_time, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, num_dg_rows_attrs)))
     goto return_status;
   if ((0 != TIO_def_var (tft->ncid, "num_tg_rows", NC_INT, 1, &dimid_time, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, num_tg_rows_attrs)))
     goto return_status;
   if ((0 != TIO_def_var (tft->ncid, "dc_storage_region", NC_FLOAT, 2, dimid_time_quad, &varid))
       || (0 != enable_var_deflation (tft->ncid, varid))
       || (0 != define_text_attrs (tft->ncid, varid, sdc_attrs)))
     goto return_status;

   if (EXPREC_TYPE_IS_DARK(exposure_type))
     {
        if ((0 != TIO_def_var (tft->ncid, "dc_mean", NC_FLOAT, 2, dimid_time_quad, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, mean_dc_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "dc_stddev", NC_FLOAT, 2, dimid_time_quad, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, stddev_dc_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "pqf_bits", NC_INT, 2, dimid_time_bit, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, pqf_bits_attrs)))
          goto return_status;
     }
   else if (EXPREC_TYPE_IS_IRRADIANCE(exposure_type))
     {
        if ((0 != TIO_def_var (tft->ncid, "solar_theta", NC_DOUBLE, 1, &dimid_time, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, solar_theta_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "solar_phi", NC_DOUBLE, 1, &dimid_time, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, solar_phi_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "use_reference_diffuser", NC_INT, 1, &dimid_time, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, ref_diffuser_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "pqf_bits_uv", NC_INT, 2, dimid_time_bit, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, pqf_bits_uv_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "pqf_bits_vis", NC_INT, 2, dimid_time_bit, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, pqf_bits_vis_attrs)))
          goto return_status;
     }
   else /* radiance */
     {
        if ((0 != TIO_def_var (tft->ncid, "pqf_bits_uv", NC_INT, 2, dimid_time_bit, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, pqf_bits_uv_attrs)))
          goto return_status;
        if ((0 != TIO_def_var (tft->ncid, "pqf_bits_vis", NC_INT, 2, dimid_time_bit, &varid))
            || (0 != enable_var_deflation (tft->ncid, varid))
            || (0 != define_text_attrs (tft->ncid, varid, pqf_bits_vis_attrs)))
          goto return_status;
     }

   Have_Trend_File = 1;
   status = 0;
return_status:
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed\n", __func__);
        trend_collect_close (tft);
        tft = NULL;
     }

   return tft;
}

void trend_collect_free_record (Trend_Record_Type *tr)
{
   if (tr == NULL)
     return;
   FREE(tr);
}

static void trend_collect_clear_record (Trend_Record_Type *tr)
{
   int i;

   if (tr == NULL)
     return;

   tr->index = -1;

   tr->exposure_type = TIO_FILL_INT;
   tr->start_time = TIO_FILL_FLOAT;
   tr->fpa_temp = TIO_FILL_FLOAT;
   tr->fpe_temp = TIO_FILL_FLOAT;
   tr->num_dg_rows = TIO_FILL_INT;
   tr->num_tg_rows = TIO_FILL_INT;
   tr->solar_theta = TIO_FILL_DOUBLE;
   tr->solar_phi = TIO_FILL_DOUBLE;
   tr->use_reference_diffuser = TIO_FILL_INT;

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        tr->eoffsets[i] = TIO_FILL_FLOAT;
        tr->gain[i] = TIO_FILL_FLOAT;
        tr->readnoise[i] = TIO_FILL_FLOAT;
     }
   for (i = 0; i < NUM_QUAD; i++)
     {
        tr->phase_change[i] = TIO_FILL_INT;
        tr->storage_region_dark[i] = TIO_FILL_FLOAT;
        tr->mean_dc[i] = TIO_FILL_FLOAT;
        tr->stddev_dc[i] = TIO_FILL_FLOAT;
     }
   for (i = 0; i < NUM_PQF_BITS; i++)
     {
        tr->pqf_bits[i] = TIO_FILL_INT;
        tr->pqf_bits_uv[i] = TIO_FILL_INT;
        tr->pqf_bits_vis[i] = TIO_FILL_INT;
     }
}

Trend_Record_Type *trend_collect_new_record (Trend_File_Type *tft)
{
   Trend_Record_Type *tr = NULL;

   if (NULL == (tr = (Trend_Record_Type *)MALLOC (sizeof *tr)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   trend_collect_clear_record (tr);
   tr->tft = tft;

   return tr;
}

int trend_collect_write_record (const Trend_Record_Type *tr)
{
   int start[2], count[2];
   int ncid, exposure_type;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   ncid = tr->tft->ncid;
   exposure_type = tr->tft->exposure_type;

   start[0] = tr->index;
   start[1] = 0;
   count[0] = 1;
   count[1] = 0;

   if ((0 != TIO_put_var_section (ncid, TEMPO_VAR_TIME, start, count, TIO_DOUBLE, &tr->start_time))
       || (0 != TIO_put_var_section (ncid, "fpa_temp", start, count, TIO_FLOAT, &tr->fpa_temp))
       || (0 != TIO_put_var_section (ncid, "fpe_temp", start, count, TIO_FLOAT, &tr->fpe_temp))
       || (0 != TIO_put_var_section (ncid, "num_dg_rows", start, count, TIO_INT, &tr->num_dg_rows))
       || (0 != TIO_put_var_section (ncid, "num_tg_rows", start, count, TIO_INT, &tr->num_tg_rows)))
     return -1;

   count[1] = NUM_OCTANTS;
   if ((0 != TIO_put_var_section (ncid, "eoffsets", start, count, TIO_FLOAT, tr->eoffsets))
       || (0 != TIO_put_var_section (ncid, "gain", start, count, TIO_FLOAT, tr->gain))
       || (0 != TIO_put_var_section (ncid, "readout_noise", start, count, TIO_FLOAT, tr->readnoise)))
     return -1;

   count[1] = NUM_QUAD;
   if ((0 != TIO_put_var_section (ncid, "phase_change", start, count, TIO_FLOAT, tr->phase_change))
       || (0 != TIO_put_var_section (ncid, "dc_storage_region", start, count, TIO_FLOAT, tr->storage_region_dark)))
     return -1;

   if (EXPREC_TYPE_IS_DARK(exposure_type))
     {
        count[1] = NUM_QUAD;
        if ((0 != TIO_put_var_section (ncid, "dc_mean", start, count, TIO_FLOAT, tr->mean_dc))
            || (0 != TIO_put_var_section (ncid, "dc_stddev", start, count, TIO_FLOAT, tr->stddev_dc)))
          return -1;

        count[1] = NUM_PQF_BITS;
        if (0 != TIO_put_var_section (ncid, "pqf_bits", start, count, TIO_INT, tr->pqf_bits))
          return -1;
     }
   else if (EXPREC_TYPE_IS_IRRADIANCE(exposure_type))
     {
        count[1] = 0;
        if ((0 != TIO_put_var_section (ncid, "solar_theta", start, count, TIO_DOUBLE, &tr->solar_theta))
            || (0 != TIO_put_var_section (ncid, "solar_phi", start, count, TIO_DOUBLE, &tr->solar_phi))
            || (0 != TIO_put_var_section (ncid, "use_reference_diffuser", start, count, TIO_INT, &tr->use_reference_diffuser)))
          return -1;

        count[1] = NUM_PQF_BITS;
        if ((0 != (TIO_put_var_section (ncid, "pqf_bits_uv", start, count, TIO_INT, tr->pqf_bits_uv)))
            || (0 != (TIO_put_var_section (ncid, "pqf_bits_vis", start, count, TIO_INT, tr->pqf_bits_vis))))
          return -1;
     }
   else /* radiance */
     {
        count[1] = NUM_PQF_BITS;
        if ((0 != (TIO_put_var_section (ncid, "pqf_bits_uv", start, count, TIO_INT, tr->pqf_bits_uv)))
            || (0 != (TIO_put_var_section (ncid, "pqf_bits_vis", start, count, TIO_INT, tr->pqf_bits_vis))))
          return -1;
     }

   return 0;
}

int trend_collect_time (double start_time, int index)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   tr->index = index;
   tr->start_time = start_time;

   return 0;
}

/*  electronic offset for each octant => eoffsets[8]
 * phase determined for each quadrant => phase_change[4]
 */
int trend_collect_eoffsets (const float *eoffsets, const int *phase_change)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   memcpy ((char *)tr->eoffsets, (char *)eoffsets, NUM_OCTANTS * sizeof(float));
   memcpy ((char *)tr->phase_change, (char *)phase_change, NUM_QUAD * sizeof(int));

   return 0;
}

/* Gain, by octant => gain[8] */
int trend_collect_gain (float fpa_temp, float fpe_temp, const float *gain)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   tr->fpa_temp = fpa_temp;
   tr->fpe_temp = fpe_temp;
   memcpy ((char *)tr->gain, (char *)gain, NUM_OCTANTS * sizeof(float));

   return 0;
}

/* Read-out noise, by octant => readnoise[8] */
int trend_collect_readnoise (const float *readnoise)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   memcpy ((char*)tr->readnoise, (char *)readnoise, NUM_OCTANTS * sizeof(float));

   return 0;
}

/* Storage region dark current, by quadrant => sdc[4] */
int trend_collect_sdc (int num_dg_rows, int num_tg_rows, const float *sdc)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   tr->num_dg_rows = num_dg_rows;
   tr->num_tg_rows = num_tg_rows;
   memcpy ((char *)tr->storage_region_dark, (char *)sdc, NUM_QUAD * sizeof(float));

   return 0;
}

/* Mean dark current, by quadrant => mean_dc[4], stddev_dc[4] */
int trend_collect_dc_mean (const float *mean_dc, const float *stddev_dc)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   memcpy ((char *)tr->mean_dc, (char *)mean_dc, NUM_QUAD * sizeof(float));
   memcpy ((char *)tr->stddev_dc, (char *)stddev_dc, NUM_QUAD * sizeof(float));

   return 0;
}

int trend_collect_solar_angles (double solar_theta, double solar_phi, int use_reference_diffuser)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   tr->solar_theta = solar_theta;
   tr->solar_phi = solar_phi;
   tr->use_reference_diffuser = use_reference_diffuser;

   return 0;
}

static void count_bits (const Image_Type *img, int nbits, int *bit_count)
{
   int n;

   for (n = 0; n < nbits; n++)
     {
        bit_count[n] = image_count_mask_pixels (img, 1 << n);
     }
}

int trend_collect_pqf (const Image_Type *img)
{
   Trend_Record_Type *tr = Active_Record;

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   count_bits (img, NUM_PQF_BITS, tr->pqf_bits);

   return 0;
}

int trend_collect_pqf_uv (Image_Pqf_Bitmap_Type *pqf, int num_rows, int num_cols)
{
   Trend_Record_Type *tr = Active_Record;
   Image_Type img = {0};

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   img.pixel_quality_flags = pqf;
   img.num_rows = num_rows;
   img.num_cols = num_cols;

   count_bits (&img, NUM_PQF_BITS, tr->pqf_bits_uv);

   return 0;
}

int trend_collect_pqf_vis (Image_Pqf_Bitmap_Type *pqf, int num_rows, int num_cols)
{
   Trend_Record_Type *tr = Active_Record;
   Image_Type img = {0};

   if (!Have_Trend_File)
     return 0;
   if (tr == NULL)
     return -1;

   img.pixel_quality_flags = pqf;
   img.num_rows = num_rows;
   img.num_cols = num_cols;

   count_bits (&img, NUM_PQF_BITS, tr->pqf_bits_vis);

   return 0;
}
