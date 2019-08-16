#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"

#define NUM_BANDS 2

#define OUTPUT_PRIVATE_DATA \
   char *file; \
   char *metadata_template_dir; \
   int exposure_type; \
   int have_dims; \
   int num_recs; \
   int num_xtrack; \
   int num_waves; \
   int formatted_file_exists; \
   int ncid; \
   double tstart; \
   double tend;
#include "output.h"
#include "util.h"

static int out_set_file (Output_Type *out, const char *file)
{
   if (NULL == (out->file = strdup (file)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   return 0;
}

static int out_set_dims (Output_Type *out, int num_recs,
                         int num_xtrack, int num_waves)
{
   if ((num_recs < 1) || (num_xtrack < 1) || (num_waves < 1))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: invalid dimensions num_recs=%d, num_xtrack=%d, num_waves=%d",
                     __func__, num_recs, num_xtrack, num_waves);
        return -1;
     }

   out->num_recs = num_recs;
   out->num_xtrack = num_xtrack;
   out->num_waves = num_waves;
   out->have_dims = 1;

   return 0;
}

static int out_close (Output_Type *out)
{
   int status;

   if ((0 != TIO_write_timestamp (out->ncid, NC_GLOBAL, "time_coverage_start", out->tstart))
       ||(0 != TIO_write_timestamp (out->ncid, NC_GLOBAL, "time_coverage_end", out->tend)))
     {
        tell_vwarn (0, "%s: writing coverage time stamps", __func__);
     }

   status = TIO_close (out->ncid);
   out->ncid = 0;
   return status;
}

static int create_file (Output_Type *out)
{
   if (out->file == NULL)
     {
        tell_verror (TELL_USAGE_ERROR,
                     "%s: output file not specified", __func__);
        return -1;
     }

   if (out->have_dims == 0)
     {
        tell_verror (TELL_USAGE_ERROR,
                     "%s: output granule dimensions not specified", __func__);
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "creating %s", out->file);

   if (0 != TIO_create (out->file, NC_NETCDF4, &out->ncid))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: creating %s",
                     __func__, out->file);
        return -1;
     }

   return 0;
}

static void
define_tio_bands (Output_Type *out,
                  TIO_Scan_Group_Type tio_bands[NUM_BANDS])
{
   tio_bands[0].name = TEMPO_BAND_NAME_UV;
   tio_bands[0].num_xtrack = out->num_xtrack;
   tio_bands[0].num_channels = out->num_waves;

   tio_bands[1].name = TEMPO_BAND_NAME_VIS;
   tio_bands[1].num_xtrack = out->num_xtrack;
   tio_bands[1].num_channels = out->num_waves;
}

typedef int TIO_Template_Method (int, size_t, int, TIO_Scan_Group_Type *);

static int
create_file_of_type (Output_Type *out,
                     TIO_Template_Method *tmpl_method)
{
   TIO_Scan_Group_Type tio_bands[NUM_BANDS];

   if (0 != create_file (out))
     return -1;

   define_tio_bands (out, tio_bands);

   if (0 != (*tmpl_method) (out->ncid, out->num_recs, NUM_BANDS, tio_bands))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: creating %s",
                     __func__, out->file);
        out->out_close (out);
        return -1;
     }

   out->formatted_file_exists = 1;

   return 0;
}

static int create_irr_file (Output_Type *out)
{
   return create_file_of_type (out, TIO_l1_irradiance_template);
}

static int create_rad_file (Output_Type *out)
{
   return create_file_of_type (out, TIO_l1_radiance_template);
}

static void out_free (Output_Type *out)
{
   if (out == NULL)
     return;
   FREE(out->metadata_template_dir);
   FREE(out->file);
   FREE(out);
}

static int out_file_exists (const Output_Type *out)
{
   return out->formatted_file_exists ? 1 : 0;
}

static int
write_rec_band1 (Output_Type *out, int index, int w_nwg,
                 const Spectral_Data_Type *sdt,
                 const char *name_var, const char *name_var_err)
{
   int grp, start[3], count[3];

   if (sdt == NULL)
     return 0;

   if (0 != TIO_inq_grp (out->ncid, sdt->name, &grp))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: accessing group %s in %s",
                     __func__, sdt->name, out->file);
        return -1;
     }

   start[0] = index;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = out->num_xtrack;
   count[2] = out->num_waves;

   if (0 != TIO_put_var_section (grp, name_var, start, count, TIO_DOUBLE,
                                 sdt->img))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, name_var, out->file);
        return -1;
     }

   if (0 != TIO_put_var_section (grp, TEMPO_VAR_PQF, start, count, TIO_USHORT,
                                 sdt->pqf))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_PQF, out->file);
        return -1;
     }

   if (0 != TIO_put_var_section (grp, name_var_err, start, count, TIO_DOUBLE,
                                 sdt->img_err))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, name_var_err, out->file);
        return -1;
     }

   if (w_nwg)
     {
        if (0 != TIO_put_var_section (grp, TEMPO_VAR_WAVELEN_NOMINAL, &start[2], &count[2], TIO_DOUBLE,
                                      sdt->wave))
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                          __func__, TEMPO_VAR_WAVELEN_NOMINAL, out->file);
             return -1;
          }
     }

   return 0;
}

static void update_coverage_time_range (Output_Type *out,
                                        const Output_Metadata_Type *meta)
{
   double tend;

   if (isnan(meta->start_time) || (meta->start_time < 0))
     return;

   if (isnan(out->tstart)
       || (meta->start_time < out->tstart))
     {
        out->tstart = meta->start_time;
     }

   if (isnan(meta->exposure_time) || (meta->exposure_time < 0))
     return;

   tend = meta->start_time + meta->exposure_time;

   if (isnan(out->tend)
       || (tend > out->tend))
     {
        out->tend = tend;
     }
}

static int write_rec_meta (Output_Type *out, int index,
                           const Output_Metadata_Type *meta)
{
   int grp, start, count;

   update_coverage_time_range (out, meta);

   /* Assuming out->tstart has been set in update_coverage_time_range,
    * set the header earth-sun distance to the value for the first record.
    */
   if (out->tstart == meta->start_time)
     {
        /* header value is in meters, computed value is in km */
        if (0 != tio_set_earth_sun_distance (out->ncid, meta->earth_sun_distance * 1.e3))
          return -1;
     }

   if (0 != TIO_inq_grp (out->ncid, "/", &grp))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: accessing group / in %s",
                     __func__, out->file);
        return -1;
     }

   start = index;
   count = 1;

   if (0 != TIO_put_var_section (grp, TEMPO_VAR_TIME, &start, &count, TIO_DOUBLE,
                                 &meta->start_time))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_TIME, out->file);
        return -1;
     }

   if (0 != TIO_put_var_section (grp, TEMPO_VAR_EXPOSURE_TIME, &start, &count, TIO_DOUBLE,
                                 &meta->exposure_time))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_EXPOSURE_TIME, out->file);
        return -1;
     }

   if (0 != TIO_put_var_section (out->ncid, TEMPO_DIM_STEP, &start, &count, TIO_INT,
                                 &meta->mirror_step))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_DIM_STEP, out->file);
        return -1;
     }

   return 0;
}

static int
write_rec_bands (Output_Type *out, int index,
                 const Output_Exprec_Type *rec,
                 const char *name_var, const char *name_var_err)
{
   int w_nwg = rec->write_nominal_wavelength_grid;

   if (0 == out_file_exists (out))
     {
        tell_verror (TELL_USAGE_ERROR, "%s: no open output file",
                     __func__);
        return -1;
     }

   if (0 != write_rec_meta (out, index, &rec->meta))
     return -1;

   if (0 != write_rec_band1 (out, index, w_nwg, rec->uv, name_var, name_var_err))
     return -1;

   if (0 != write_rec_band1 (out, index, w_nwg, rec->vis, name_var, name_var_err))
     return -1;

   return 0;
}

static int write_irr_rec (Output_Type *out,
                          int index, const Output_Exprec_Type *rec)
{
   return write_rec_bands (out, index, rec,
                           TEMPO_VAR_IRRADIANCE,
                           TEMPO_VAR_IRRADIANCE_ERROR);
}

static int write_rad_rec (Output_Type *out,
                          int index, const Output_Exprec_Type *rec)
{
   return write_rec_bands (out, index, rec,
                           TEMPO_VAR_RADIANCE,
                           TEMPO_VAR_RADIANCE_ERROR);
}

static int out_get_ncid (const Output_Type *out)
{
   return out->ncid;
}

static int out_std_metadata (Output_Type *out, TIO_Meta_Type *meta, int ncid_from)
{
#define SHORTNAME_BUFSIZE 32
   char shortname[SHORTNAME_BUFSIZE];
   const char *prod_name = NULL;
   const char *template_basename = NULL;
   char *template_path = NULL;
   int grp_meta, n, status = -1;

   switch (out->exposure_type)
     {
      case EXPREC_TYPE_IRR_WRK:
        prod_name = TEMPO_PROD_TYPE_IRR;
        template_basename = "irradiance.met.template";
        break;

      case EXPREC_TYPE_IRR_REF:
        prod_name = TEMPO_PROD_TYPE_IRR_REF;
        template_basename = "irradiance.met.template";
        break;

      case EXPREC_TYPE_RAD:
        prod_name = TEMPO_PROD_TYPE_RAD;
        template_basename = "radiance.met.template";
        break;

      default:
        tell_vwarn (0, "%s: no metadata template expansion support for exposure records of type %d",
                    __func__, out->exposure_type);
        break;
     }

   n = snprintf (shortname, SHORTNAME_BUFSIZE, "TEMPO_%s", prod_name);
   if (n < 0 || n >= SHORTNAME_BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error generating shortname for %s", __func__, prod_name);
        goto return_status;
     }

   /* FIXME: set version numbers */
   if (0 != tio_meta_set_standard (meta, out->file, shortname, 1, "0.1.0"))
     goto return_status;

   if (0 != tio_meta_set_datetime_range (meta, ncid_from))
     goto return_status;

   if ((0 != TIO_def_grp (out->ncid, "metadata", &grp_meta))
       || (0 != tio_meta_write_ncattr (meta, grp_meta)))
     goto return_status;

   if ((out->exposure_type == EXPREC_TYPE_RAD)
       || (out->exposure_type == EXPREC_TYPE_IRR_WRK)
       || (out->exposure_type == EXPREC_TYPE_IRR_REF))
     {
        /* For radiance files, INPUTPOINTER gets expanded only in the
         * last processing step of Level 0-1, e.g. post-INR
         * For irradiance files, INPUTPOINTER gets expanded after
         * wavelength calibration.
         */
        tio_meta_set_noexpand (meta, "INPUTPOINTER", 1);
     }

   if ((out->metadata_template_dir != NULL)
       && (template_basename != NULL))
     {
        if (NULL == (template_path = path_concat (out->metadata_template_dir, template_basename)))
          goto return_status;
        tell_vlog (TELL_MSGTYPE_INFO, 1, "Expanding metadata template: %s", template_path);
        if (0 != tio_meta_expand_file (meta, template_path, out->file))
          goto return_status;
     }

   status = 0;
return_status:
   FREE(template_path);
   return status;
}

static int read_params (Output_Type *out, config_t *cfg)
{
   config_setting_t *setting;
   const char *template_dir;

   if (NULL == (setting = config_lookup (cfg, "metadata")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing group 'template' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   out->metadata_template_dir = NULL;

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "metadata_template_dir", &template_dir))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "metadata template path not found: skipping template expansion");
        return 0;
     }

   if (NULL == (out->metadata_template_dir = expand_path (template_dir)))
     return -1;

   return 0;
}

Output_Type *output_alloc (config_t *cfg, int exposure_type)
{
   Output_Type *out = NULL;
   double nan_value = nan("");

   (void) cfg;

   if (NULL == (out = (Output_Type *) MALLOC (sizeof *out)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)out, 0, sizeof *out);

   out->exposure_type = exposure_type;

   out->out_free = out_free;
   out->out_set_file = out_set_file;
   out->out_set_dims = out_set_dims;
   out->out_close = out_close;
   out->out_file_exists = out_file_exists;
   out->tstart = nan_value;
   out->tend = nan_value;
   out->out_ncid = out_get_ncid;
   out->out_std_metadata = out_std_metadata;

   if (0 != read_params (out, cfg))
     {
        out_free (out);
        return NULL;
     }

   switch (exposure_type)
     {
      default:
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR,
                     "%s: no support for exposure type %d",
                     __func__, exposure_type);
        FREE(out);
        return NULL;

      case EXPREC_TYPE_RAD:
        out->out_create = create_rad_file;
        out->out_write_rec = write_rad_rec;
        break;

      case EXPREC_TYPE_IRR_WRK: /* drop */
      case EXPREC_TYPE_IRR_REF:
        out->out_create = create_irr_file;
        out->out_write_rec = write_irr_rec;
        break;
     }

   return out;
}
