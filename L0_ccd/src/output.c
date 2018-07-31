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
   FREE(out->file);
   FREE(out);
}

static int out_file_exists (const Output_Type *out)
{
   return out->formatted_file_exists ? 1 : 0;
}

static int
write_rec_band1 (Output_Type *out, int index,
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

   if (0 != TIO_put_var_section (grp, TEMPO_VAR_WAVELENGTH, start, count, TIO_DOUBLE,
                                 sdt->wave))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_WAVELENGTH, out->file);
        return -1;
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
   if (0 == out_file_exists (out))
     {
        tell_verror (TELL_USAGE_ERROR, "%s: no open output file",
                     __func__);
        return -1;
     }

   if (0 != write_rec_meta (out, index, &rec->meta))
     return -1;

   if (0 != write_rec_band1 (out, index, rec->uv, name_var, name_var_err))
     return -1;

   if (0 != write_rec_band1 (out, index, rec->vis, name_var, name_var_err))
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

   out->out_free = out_free;
   out->out_set_file = out_set_file;
   out->out_set_dims = out_set_dims;
   out->out_close = out_close;
   out->out_file_exists = out_file_exists;
   out->tstart = nan_value;
   out->tend = nan_value;
   out->out_ncid = out_get_ncid;

   switch (exposure_type)
     {
      default:
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR,
                     "%s: no support for exposure type %d",
                     __func__, exposure_type);
        FREE(out);
        return NULL;

      case EXPREC_TYPE_RADIANCE:
        out->out_create = create_rad_file;
        out->out_write_rec = write_rad_rec;
        break;

      case EXPREC_TYPE_LIN_IRR:
        /* drop */
      case EXPREC_TYPE_IRRADIANCE:
        out->out_create = create_irr_file;
        out->out_write_rec = write_irr_rec;
        break;
     }

   return out;
}
