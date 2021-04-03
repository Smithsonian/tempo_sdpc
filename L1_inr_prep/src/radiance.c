/** @file radiance.c
 *  @brief Manage radiance file I/O
 */
#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libconfig.h>
#include <ioclib.h>
#include <tell.h>

#include <tio.h>
#include <tio_template.h>

#define RADIANCE_PRIVATE_DATA \
   double tstart; \
   double tstop; \
   int processing_version; \
   int created_for_inr_status_update;
#include "radiance.h"

static int meta_record_basename (TIO_Meta_Type *meta, const char *path)
{
   const char *path_basename;

   if (meta == NULL)
     return 0;

   if (path == NULL)
     return 0;

   if (NULL != (path_basename = strrchr (path, '/')))
     {
        path_basename++;
     }
   else path_basename = path;

   return tio_meta_append_string (meta, "input_files", path_basename);
}

static void free_radiance_type (Radiance_Type *r)
{
   if (r == NULL)
     return;
   ioclib_free (r->file);
   FREE(r);
}

static Radiance_Type *alloc_radiance_type (const char *file)
{
   Radiance_Type *r = NULL;
   double nan_value = nan("");

   if (NULL == (r = (Radiance_Type *)MALLOC (sizeof *r)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)r, 0, sizeof (*r));

   r->tstart = nan_value;
   r->tstop = nan_value;

   if (NULL == (r->file = ioclib_strdup (file)))
     {
        FREE(r);
        return NULL;
     }

   return r;
}

int radiance_update_coverage_times (Radiance_Type *r, double tstart, double tstop)
{
   if (0 == r->created_for_inr_status_update)
     return 0;

   if ((0 != TIO_write_timestamp (r->ncid, NC_GLOBAL, "time_coverage_start", tstart))
       || (0 != TIO_write_timestamp (r->ncid, NC_GLOBAL, "time_coverage_end", tstop)))
     return -1;

   return 0;
}

void radiance_close (Radiance_Type *r)
{
   if (r == NULL)
     return;

   (void) TIO_close (r->ncid);
   free_radiance_type (r);
}

Radiance_Type *radiance_open (const char *file)
{
   Radiance_Type *r = NULL;

   if (NULL == (r = alloc_radiance_type (file)))
     return NULL;

   if (0 != TIO_open (file, NC_WRITE, &r->ncid))
     goto return_error;

   if (0 != tio_history_append_cmdline (r->ncid))
     goto return_error;

   if (0 != tio_use_file_epoch (r->ncid))
     goto return_error;

   if (0 != TIO_get_att (r->ncid, NC_GLOBAL, "time_coverage_start_since_epoch",
                         NC_DOUBLE, &r->tstart))
     goto return_error;

   if (0 != TIO_get_att (r->ncid, NC_GLOBAL, "time_coverage_end_since_epoch",
                         NC_DOUBLE, &r->tstop))
     goto return_error;

   return r;
return_error:
   free_radiance_type (r);
   return NULL;
}

static TIO_Scan_Group_Type Scan_Groups[] =
{
   {TEMPO_BAND_NAME_UV,1,1},
   {TEMPO_BAND_NAME_VIS,1,1}
};

#define IS_TELEMETRY_ONLY  4   /* granule_flag */

Radiance_Type *radiance_create (const char *file, int processing_version)
{
   Radiance_Type *r = NULL;
   int num_scan_groups, num_steps;
   int start, count, granule_flag = IS_TELEMETRY_ONLY;

   if (NULL == (r = alloc_radiance_type (file)))
     return NULL;

   r->processing_version = processing_version;
   r->created_for_inr_status_update = 1;

   if (0 != TIO_create (file, NC_NETCDF4, &r->ncid))
     goto return_error;

   num_scan_groups = sizeof(Scan_Groups)/sizeof(TIO_Scan_Group_Type);
   num_steps = 0;

   if (0 != TIO_l1_radiance_template (r->ncid, num_steps,
                                      num_scan_groups, Scan_Groups))
     goto return_error;

   start = 0;
   count = 1;
   if (0 != TIO_put_var_section (r->ncid, TEMPO_VAR_GRANULE_FLAG, &start, &count,
                                 NC_INT, &granule_flag))
     goto return_error;

   if (0 != TIO_put_att (r->ncid, NC_GLOBAL, "processing_version",
                         NC_INT, 1, &r->processing_version))
     goto return_error;

   return r;
return_error:
   free_radiance_type (r);
   return NULL;
}

int radiance_interval (Radiance_Type *r, double *tstart, double *tstop)
{
   if (r == NULL)
     return -1;

   *tstart = r->tstart;
   *tstop = r->tstop;
   return 0;
}

static int radiance_write_times (Radiance_Type *r,
                                 const Row_Select_Type *rst_head, int grp,
                                 const char *time_var_name)
{
   const Row_Select_Type *rst;
   double tstart_global, tstop_global;
   double *rst_times;
   int count, start;

   tstart_global = ((rst_head->start < rst_head->num_times)
                    ? rst_head->times[rst_head->start]
                    : nan(""));
   tstop_global = tstart_global;

   start = 0;
   for (rst = rst_head; rst != NULL; rst = rst->next)
     {
        if (rst->start >= rst->num_times)
          continue;

        rst_times = rst->times + rst->start;
        count = rst->count;
        if (isnan(tstart_global)) tstart_global = rst_times[0];
        tstop_global = rst_times[count-1];

        if (0 != TIO_put_var_section (grp, time_var_name, &start, &count,
                                      NC_DOUBLE, rst_times))
          return -1;
        start += count;
     }

   if (isnan(tstart_global))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: empty time series! (this should never happen)",
                     __func__);
        return -1;
     }

   if (r->created_for_inr_status_update)
     {
        if (isnan(r->tstart) || (tstart_global < r->tstart))
          r->tstart = tstart_global;
        if (isnan(r->tstop) || (tstop_global > r->tstop))
          r->tstop = tstop_global;
     }

   return 0;
}

typedef struct
{
   void *bytes;
   size_t num_bytes;
}
Buffer_Type;

static int copy_file_var (const Row_Select_Type *rst, Buffer_Type *buf,
                          const char *from_var, int from_ncid, int copy_all,
                          const char *to_var, int to_grp, int to_start)
{
   TIO_Var_Info_Type info;
   size_t num_elem;
   int i, size_elem;
   int start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];

   if (0 != TIO_inq_var (from_ncid, from_var, &info))
     return -1;

   switch (info.type)
     {
      case NC_DOUBLE:
        size_elem = 8;
        break;

      case NC_FLOAT:
      case NC_INT:
      case NC_UINT:
        size_elem = 4;
        break;

      case NC_SHORT:
      case NC_USHORT:
        size_elem = 2;
        break;

      case NC_CHAR:
      case NC_BYTE:
      case NC_UBYTE:
        size_elem = 1;
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported type = %d",
                     __func__, info.type);
        return -1;
     }

   if (copy_all)
     {
        start[0] = 0;
        count[0] = info.dimlens[0];
     }
   else
     {
        start[0] = rst->start;
        count[0] = rst->count;
     }

   if (start[0] >= rst->num_times)
     return 0;

   num_elem = info.dimlens[0];
   for (i = 1; i < info.ndims; i++)
     {
        start[i] = 0;
        count[i] = info.dimlens[i];
        num_elem *= info.dimlens[i];
     }

   if (buf->num_bytes < num_elem * size_elem)
     {
        size_t num_bytes = num_elem * size_elem;
        void *v;
        if (NULL == (v = REALLOC(buf->bytes, num_bytes)))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return -1;
          }
        buf->bytes = v;
        buf->num_bytes = num_bytes;
     }

   if (0 != TIO_get_var_section (from_ncid, from_var, start, count,
                                 info.type, buf->bytes))
     return -1;

   start[0] = to_start;

   if (0 != TIO_put_var_section (to_grp, to_var, start, count,
                                 info.type, buf->bytes))
     return -1;

   return count[0];
}

typedef struct
{
   const char *from;
   const char *to;
}
Var_Name_Type;
#define VAR_TABLE_END {NULL,NULL}

static Var_Name_Type SMC_Vars[] =
{
   {TEMPO_VAR_SMADIT_SCANX, TEMPO_VAR_SMADIT_SCANX},
   {TEMPO_VAR_SMADIT_SCANY, TEMPO_VAR_SMADIT_SCANY},
   {TEMPO_VAR_SMADIT_RAWX, TEMPO_VAR_SMADIT_RAWX},
   {TEMPO_VAR_SMADIT_RAWY, TEMPO_VAR_SMADIT_RAWY},
   VAR_TABLE_END
};

static Var_Name_Type IRU_Vars[] =
{
   {TEMPO_VAR_GYRO_OUTPUT, TEMPO_VAR_GYRO_OUTPUT},
   VAR_TABLE_END
};

static Var_Name_Type IRU_Bias_Vars[] =
{
   {"time_of_scale_or_bias_update", TEMPO_VAR_TIME_GYRO_BIAS},
   {"bias", TEMPO_VAR_GYRO_BIAS},
   VAR_TABLE_END
};

static Var_Name_Type EPH_Vars[] =
{
   {"anc_satx", TEMPO_VAR_SAT_X},
   {"anc_saty", TEMPO_VAR_SAT_Y},
   {"anc_satz", TEMPO_VAR_SAT_Z},
   {"anc_satvx", TEMPO_VAR_SAT_VX},
   {"anc_satvy", TEMPO_VAR_SAT_VY},
   {"anc_satvz", TEMPO_VAR_SAT_VZ},
   VAR_TABLE_END
};

static int radiance_copy_vars (const Row_Select_Type *rst_head, TIO_Meta_Type *meta,
                               Var_Name_Type *v, const char *from_group_path,
                               int to_grp, int copy_all)
{
   const Row_Select_Type *rst;
   Buffer_Type buf = {0};
   int from_ncid = 0;
   int to_start = 0;
   int status = -1;

   for (rst = rst_head; rst != NULL; rst = rst->next)
     {
        int from_grp, rows_copied;

        if (0 != TIO_open (rst->file, NC_NOWRITE, &from_ncid))
          goto free_and_return;
        if (from_group_path)
          {
             if (0 != TIO_inq_grp (from_ncid, from_group_path, &from_grp))
               goto free_and_return;
          }
        else from_grp = from_ncid;

        if (0 != meta_record_basename (meta, rst->file))
          goto free_and_return;
        rows_copied = copy_file_var (rst, &buf, v->from, from_grp, copy_all,
                                     v->to, to_grp, to_start);
        (void) TIO_close (from_ncid);
        from_ncid = 0;
        if (rows_copied < 0)
          goto free_and_return;
        to_start += rows_copied;
     }

   status = 0;
free_and_return:
   FREE(buf.bytes);
   if (from_ncid != 0)
     {
        (void) TIO_close (from_ncid);
     }
   return status;
}

int radiance_copy_smc (Radiance_Type *r, TIO_Meta_Type *meta, const Row_Select_Type *rst_head)
{
   Var_Name_Type *v;
   int grp;

   if (0 != TIO_inq_grp (r->ncid, "/inr_input/telemetry/mirror", &grp))
     return -1;

   if (0 != radiance_write_times (r, rst_head, grp, TEMPO_VAR_TIME_SMA))
     return -1;

   for (v = SMC_Vars; v->from != NULL; v++)
     {
        if (0 != radiance_copy_vars (rst_head, meta, v, NULL, grp, 0))
          return -1;
     }

   return 0;
}

int radiance_copy_iru (Radiance_Type *r, TIO_Meta_Type *meta, const Row_Select_Type *rst_head)
{
   Var_Name_Type *v;
   int grp;

   if (0 != TIO_inq_grp (r->ncid, "/inr_input/telemetry/gyroscope", &grp))
     return -1;

   if (0 != radiance_write_times (r, rst_head, grp, TEMPO_VAR_TIME_GYRO))
     return -1;

   for (v = IRU_Vars; v->from != NULL; v++)
     {
        if (0 != radiance_copy_vars (rst_head, meta, v, NULL, grp, 0))
          return -1;
     }

   for (v = IRU_Bias_Vars; v->from != NULL; v++)
     {
        if (0 != radiance_copy_vars (rst_head, meta, v, NULL, grp, 1))
          return -1;
     }

   return 0;
}

int radiance_copy_eph (Radiance_Type *r, TIO_Meta_Type *meta, const char *from_group_path,
                       const Row_Select_Type *eph_head)
{
   Var_Name_Type *v;
   int grp;

   if (0 != TIO_inq_grp (r->ncid, "/inr_input/ephemeris", &grp))
     return -1;

   if (0 != radiance_write_times (r, eph_head, grp, TEMPO_VAR_TIME_EPHEM))
     return -1;

   for (v = EPH_Vars; v->from != NULL; v++)
     {
        if (0 != radiance_copy_vars (eph_head, meta, v, from_group_path, grp, 0))
          return -1;
     }

   return 0;
}

