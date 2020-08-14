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
#include "process.h"
#include "image.h"

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
   int ncid_irr_frames; \
   double tstart; \
   double tend;
#include "output.h"
#include "util.h"

static int out_frames_ncid (const Output_Type *out)
{
   if (out->exposure_type == EXPREC_TYPE_IRR_WRK)
     return out->ncid_irr_frames;
   else return out->ncid;
}

static int out_root_ncid (const Output_Type *out)
{
   return out->ncid;
}

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

   return tio_history_append_cmdline (out->ncid);
}

static void
define_tio_bands (const Output_Type *out,
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
create_file_of_type (Output_Type *out, int ncid, int num_recs,
                     TIO_Template_Method *tmpl_method)
{
   TIO_Scan_Group_Type tio_bands[NUM_BANDS];

   define_tio_bands (out, tio_bands);

   if (0 != (*tmpl_method) (ncid, num_recs, NUM_BANDS, tio_bands))
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
   if (0 != create_file (out))
     return -1;

   switch (out->exposure_type)
     {
      case EXPREC_TYPE_IRR_REF:
        return create_file_of_type (out, out->ncid, out->num_recs, TIO_l1_ref_irradiance_template);

      case EXPREC_TYPE_IRR_WRK:
        /* The irradiance file used for processing has a single frame at the
         * top level to store the average irradiance.  The individual frames
         * are stored in the /frames group, which has the same uv/vis band structure
         */
        if (0 != create_file_of_type (out, out->ncid, 1, TIO_l1_wrk_irradiance_template))
          return -1;
        if (0 != TIO_def_grp (out->ncid, "frames", &out->ncid_irr_frames))
          return -1;
        return create_file_of_type (out, out->ncid_irr_frames, out->num_recs, TIO_l1_wrk_irradiance_template);

      default:
        /* No other Level 1 irradiance types are expected.
         * Irradiance data acquired for checking/correcting linearity
         * will not be processed beyond Level 0.
         */
        break;
     }

   tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported exposure type: %d",
                __func__, out->exposure_type);
   return -1;
}

static int create_rad_file (Output_Type *out)
{
   if (0 != create_file (out))
     return -1;
   return create_file_of_type (out, out->ncid, out->num_recs, TIO_l1_radiance_template);
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
   int grp, ncid, start[3], count[3];

   if (sdt == NULL)
     return 0;

   ncid = out_frames_ncid (out);

   if (0 != TIO_inq_grp (ncid, sdt->name, &grp))
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

        /* Yes, this is a particularly ugly hack */
        if (out->exposure_type == EXPREC_TYPE_IRR_WRK)
          {
             int avg_grp;
             if ((0 != TIO_inq_grp (out->ncid, sdt->name, &avg_grp))
                 || (0 != TIO_put_var_section (avg_grp, TEMPO_VAR_WAVELEN_NOMINAL, &start[2], &count[2], TIO_DOUBLE,
                                               sdt->wave)))
               {
                  tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                               __func__, TEMPO_VAR_WAVELEN_NOMINAL, out->file);
                  return -1;
               }
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
   int ncid, start, count;

   ncid = out_frames_ncid (out);

   update_coverage_time_range (out, meta);

   /* Assuming out->tstart has been set in update_coverage_time_range,
    * set the header earth-sun distance to the value for the first record.
    */
   if (out->tstart == meta->start_time)
     {
        /* header value is in meters, computed value is in km */
        if (0 != tio_set_earth_sun_distance (ncid, meta->earth_sun_distance * 1.e3))
          return -1;
     }

   start = index;
   count = 1;

   if (0 != TIO_put_var_section (ncid, TEMPO_VAR_TIME, &start, &count, TIO_DOUBLE,
                                 &meta->start_time))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_TIME, out->file);
        return -1;
     }

   if (0 != TIO_put_var_section (ncid, TEMPO_VAR_EXPOSURE_TIME, &start, &count, TIO_DOUBLE,
                                 &meta->exposure_time))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_EXPOSURE_TIME, out->file);
        return -1;
     }

   if (0 != TIO_put_var_section (ncid, TEMPO_DIM_STEP, &start, &count, TIO_INT,
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

static int out_std_metadata (Output_Type *out, TIO_Meta_Type *meta, int ncid_from)
{
#define SHORTNAME_BUFSIZE 32
   char shortname[SHORTNAME_BUFSIZE];
   const char *prod_name = NULL;
   const char *template_basename = NULL;
   char *template_path = NULL;
   int n, status = -1;

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

   n = snprintf (shortname, SHORTNAME_BUFSIZE, "TEMPO_%s_L1", prod_name);
   if (n < 0 || n >= SHORTNAME_BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error generating shortname for %s", __func__, prod_name);
        goto return_status;
     }

   /* FIXME: set version numbers */
   if (0 != tio_meta_set_standard (meta, out->file, shortname, process_get_version(), "0.1.0"))
     goto return_status;

   if (0 != tio_meta_set_datetime_range (meta, ncid_from))
     goto return_status;

   if (0 != tio_meta_write_ncattr (meta, out->ncid))
     goto return_status;

   if ((out->exposure_type == EXPREC_TYPE_RAD)
       || (out->exposure_type == EXPREC_TYPE_IRR_WRK)
       || (out->exposure_type == EXPREC_TYPE_IRR_REF))
     {
        /* For radiance files, input_pointer gets expanded only in the
         * last processing step of Level 0-1, e.g. post-INR
         * For irradiance files, input_pointer gets expanded after
         * wavelength calibration.
         */
        tio_meta_set_noexpand (meta, "input_files", 1);
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

static int insert_comment (int grp, const char *varname, const char *comment)
{
   int varid;
   if (0 != tio_inq_varid (grp, varname, &varid))
     return -1;
   if (0 != TIO_put_att (grp, varid, "comment", NC_CHAR, 1+strlen(comment), comment))
     return -1;
   return 0;
}

static int band_average_irradiance (Output_Type *out, const char *band_name)
{
   int src_grp, dest_grp;
   float *irr = NULL;
   float *irr_err = NULL;
   float *tmp_irr = NULL;
   float *tmp_irr_err = NULL;
   float *exptime = NULL;
   double *obstime = NULL;
   double obstime_avg, tot_obstime, tot_exptime;
   float exptime_avg, earth_sun_distance;
   unsigned short *pqf = NULL;
   unsigned short *tmp_pqf = NULL;
   char comment_string[1024];
   int i, k, len, start[3], count[3], n_obstime, n_exptime;
   int processing_version;
   int *num = NULL;
   int status = -1;

   /* average over frames stored in /frames/$band_name */
   if (0 != TIO_inq_grp (out->ncid_irr_frames, band_name, &src_grp))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: accessing group %s in %s",
                     __func__, band_name, out->file);
        return -1;
     }

   /* store averages in /$band_name */
   if (0 != TIO_inq_grp (out->ncid, band_name, &dest_grp))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: accessing group %s in %s",
                     __func__, band_name, out->file);
        return -1;
     }

   len = out->num_xtrack * out->num_waves;
   if ((NULL == (irr = (float *)MALLOC (2*len * sizeof(float))))
       || (NULL == (tmp_irr = (float *)MALLOC (2*len * sizeof(float))))
       || (NULL == (pqf = (unsigned short *)MALLOC (2*len * sizeof(unsigned short))))
       || (NULL == (num = (int *)MALLOC (len * sizeof(int))))
       || (NULL == (obstime = (double *)MALLOC (out->num_recs * sizeof(double))))
       || (NULL == (exptime = (float *)MALLOC (out->num_recs * sizeof(float))))
      )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc error", __func__);
        goto return_status;
     }
   memset ((char *)irr, 0, 2*len * sizeof(float));
   memset ((char *)tmp_irr, 0, 2*len * sizeof(float));
   memset ((char *)pqf, 0, 2*len * sizeof(unsigned short));
   memset ((char *)num, 0, len * sizeof(int));

   irr_err = irr + len;
   tmp_irr_err = tmp_irr + len;
   tmp_pqf = pqf + len;

   for (i = 0; i < out->num_recs; i++)
     {
        start[0] = i;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = out->num_xtrack;
        count[2] = out->num_waves;

        if ((0 != TIO_get_var_section (src_grp, TEMPO_VAR_PQF, start, count, TIO_USHORT, tmp_pqf))
            || (0 != TIO_get_var_section (src_grp, TEMPO_VAR_IRRADIANCE, start, count, TIO_FLOAT, tmp_irr))
            || (0 != TIO_get_var_section (src_grp, TEMPO_VAR_IRRADIANCE_ERROR, start, count, TIO_FLOAT, tmp_irr_err)))
          goto return_status;

        for (k = 0; k < len; k++)
          {
             if (tmp_pqf[k] == 0)
               {
                  num[k]     += 1;
                  irr[k]     += tmp_irr[k];
                  irr_err[k] += tmp_irr_err[k] * tmp_irr_err[k];
               }
             else
               {
                  pqf[k] |= tmp_pqf[k];
               }
          }
     }

   for (k = 0; k < len; k++)
     {
        if (num[k] > 0)
          {
             irr[k] /= num[k];
             irr_err[k] = sqrt(irr_err[k] /num[k]);
          }
        else
          {
             irr[k] = IMAGE_PIXEL_FILL_VALUE;
             irr_err[k] = IMAGE_PIXEL_FILL_VALUE;
          }
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = out->num_xtrack;
   count[2] = out->num_waves;

   if (0 != TIO_put_var_section (dest_grp, TEMPO_VAR_PQF, start, count, TIO_USHORT, pqf))
     goto return_status;
   if (0 != TIO_put_var_section (dest_grp, TEMPO_VAR_IRRADIANCE, start, count, TIO_FLOAT, irr))
     goto return_status;
   if (0 != TIO_put_var_section (dest_grp, TEMPO_VAR_IRRADIANCE_ERROR, start, count, TIO_FLOAT, irr_err))
     goto return_status;

   snprintf (comment_string, sizeof(comment_string), "Logical OR of /frames/%s/%s", band_name, TEMPO_VAR_PQF);
   if (0 != insert_comment (dest_grp, TEMPO_VAR_PQF, comment_string))
     goto return_status;
   snprintf (comment_string, sizeof(comment_string), "Average of /frames/%s/%s", band_name, TEMPO_VAR_IRRADIANCE);
   if (0 != insert_comment (dest_grp, TEMPO_VAR_IRRADIANCE, comment_string))
     goto return_status;
   snprintf (comment_string, sizeof(comment_string), "Root mean square of /frames/%s/%s", band_name, TEMPO_VAR_IRRADIANCE_ERROR);
   if (0 != insert_comment (dest_grp, TEMPO_VAR_IRRADIANCE_ERROR, comment_string))
     goto return_status;

   start[0] = 0;
   count[0] = out->num_recs;
   if ((0 != TIO_get_var_section (out->ncid_irr_frames, TEMPO_VAR_TIME, start, count, TIO_DOUBLE, obstime))
       ||(0 != TIO_get_var_section (out->ncid_irr_frames, TEMPO_VAR_EXPOSURE_TIME, start, count, TIO_FLOAT, exptime)))
     goto return_status;

   tot_obstime = 0.0;
   tot_exptime = 0.0;
   n_obstime = 0;
   n_exptime = 0;
   for (i = 0; i < out->num_recs; i++)
     {
        if (obstime[i] != TIO_FILL_DOUBLE)
          {
             tot_obstime += obstime[i];
             n_obstime++;
          }
        if (exptime[i] != TIO_FILL_FLOAT)
          {
             tot_exptime += exptime[i];
             n_exptime++;
          }
     }
   obstime_avg = n_obstime ? (tot_obstime / n_obstime) : TIO_FILL_DOUBLE;
   exptime_avg = n_exptime ? (tot_exptime / n_exptime) : TIO_FILL_FLOAT;

   start[0] = 0;
   count[0] = 1;
   if ((0 != TIO_put_var_section (out->ncid, TEMPO_VAR_TIME, start, count, TIO_DOUBLE, &obstime_avg))
       ||(0 != TIO_put_var_section (out->ncid, TEMPO_VAR_EXPOSURE_TIME, start, count, TIO_FLOAT, &exptime_avg)))
     goto return_status;

   snprintf (comment_string, sizeof(comment_string), "Average of /frames/%s", TEMPO_VAR_TIME);
   if (0 != insert_comment (out->ncid, TEMPO_VAR_TIME, comment_string))
     goto return_status;
   snprintf (comment_string, sizeof(comment_string), "Average of /frames/%s", TEMPO_VAR_EXPOSURE_TIME);
   if (0 != insert_comment (out->ncid, TEMPO_VAR_EXPOSURE_TIME, comment_string))
     goto return_status;

   i = 1;
   if (0 != TIO_put_var_section (out->ncid, TEMPO_DIM_STEP, start, count, TIO_INT, &i))
     goto return_status;

   if (0 != TIO_get_var_section (out->ncid_irr_frames, TEMPO_VAR_EARTH_SUN_DISTANCE, start, count, TIO_FLOAT, &earth_sun_distance))
     goto return_status;
   if (0 != TIO_put_var_section (out->ncid, TEMPO_VAR_EARTH_SUN_DISTANCE, start, count, TIO_FLOAT, &earth_sun_distance))
     goto return_status;

   if ((0 != TIO_write_timestamp (out->ncid_irr_frames, NC_GLOBAL, "time_coverage_start", out->tstart))
       ||(0 != TIO_write_timestamp (out->ncid_irr_frames, NC_GLOBAL, "time_coverage_end", out->tend)))
     {
        tell_vwarn (0, "%s: writing coverage time stamps to frames group", __func__);
     }

   /* We don't need this, but since it's in the file, let's try to get it right */
   processing_version = process_get_version();
   if (0 != TIO_put_att (out->ncid_irr_frames, NC_GLOBAL, "processing_version", NC_INT, 1, &processing_version))
     {
        tell_vwarn (0, "%s: writing processing_version to frames group", __func__);
     }

   status = 0;
return_status:
   FREE(irr);
   FREE(tmp_irr);
   FREE(pqf);
   FREE(num);
   FREE(obstime);
   FREE(exptime);

   return status;
}

static int generate_average_irradiance (Output_Type *out)
{
   if ((0 != band_average_irradiance (out, TEMPO_BAND_NAME_UV))
       || (0 != band_average_irradiance (out, TEMPO_BAND_NAME_VIS)))
     return -1;

   return 0;
}

static int out_finalize (Output_Type *out)
{
   switch (out->exposure_type)
     {
      case EXPREC_TYPE_IRR_WRK:
        return generate_average_irradiance (out);

      default:
        break;
     }

   return 0;
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

   if (NULL == (out->metadata_template_dir = expand_string (template_dir)))
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
   out->out_root_ncid = out_root_ncid;
   out->out_std_metadata = out_std_metadata;
   out->out_finalize = out_finalize;

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
