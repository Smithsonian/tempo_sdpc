/** @file radiance.c
 *  @brief Manage radiance file I/O
 */
#include "config.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libconfig.h>
#include <ioclib.h>
#include <iocsdpc.h>
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

static int fixup_granule_flag (int ncid)
{
   int start, count, granule_num, granule_flag;

   /* Early in the TEMPO mission, the scan_seq_start bit was set in all
    * granules of a scan for which an INR warm restart was needed
    * (e.g. first scan of each day, or first scan following a maneuver).
    * On the other hand, the INR server triggered a warm restart only
    * when it encountered a subsequent granule with the scan_seq_start
    * bit not set. There was a misunderstanding here. The INR server implementation
    * assumed the scan_seq_start bit would only be set in the first granule
    * of the scan (handling of a single-granule scan might have been a problematic
    * corner case, but let's not go there right now). Later, the INR SW was modified
    * to perform an automatic restart after any sufficiently long idle period,
    * so setting the scan_seq_start bit was mostly not needed, and this mechanism
    * became much less important. To facilitate reprocessing of radiance scans collected
    * early in the mission, we fixup the misunderstanding by adjusting the bit setting.
    * We ensure that the scan_seq_start bit is set only in the first granule of any scan.
    * Everywhere else, we un-set it.
    */
   if (0 != TIO_get_att (ncid, NC_GLOBAL, "granule_num", NC_INT, &granule_num))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unable to read granule_num", __func__);
        return -1;
     }

   /* We make changes only in granules that aren't the first */
   if (granule_num == 1)
     return 0;

   start = 0;
   count = 1;
   if (0 != TIO_get_var_section (ncid, TEMPO_VAR_GRANULE_FLAG, &start, &count, NC_INT, &granule_flag))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unable to read granule_flag", __func__);
        return -1;
     }

   /* We make no change when the scan_seq_start bit is not set */
   if ((granule_flag & TEMPO_GRANULE_FLAG_SCAN_SEQ_START) == 0)
     return 0;

   /* When the scan_seq_start bit is set, we clear it, and update
    * the file before delivery to INR */
   granule_flag &= ~TEMPO_GRANULE_FLAG_SCAN_SEQ_START;
   if (0 != TIO_put_var_section (ncid, TEMPO_VAR_GRANULE_FLAG, &start, &count, NC_INT, &granule_flag))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unable to write granule_flag", __func__);
        return -1;
     }

   return 0;
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

   if (0 != fixup_granule_flag (r->ncid))
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

static int apply_scale_factor (double factor, int type, void *bytes, size_t num_elem)
{
   float *fval = NULL;
   double *dval = NULL;
   size_t i;

   switch (type)
     {
      case NC_FLOAT:
        fval = (float *)bytes;
        for (i = 0; i < num_elem; i++)
          {
             if ((0 == isnan(fval[i])) && (fval[i] != NC_FILL_FLOAT))
               {
                  fval[i] *= factor;
               }
          }
        break;

      case NC_DOUBLE:
        dval = (double *)bytes;
        for (i = 0; i < num_elem; i++)
          {
             if ((0 == isnan(dval[i])) && (dval[i] != NC_FILL_DOUBLE))
               {
                  dval[i] *= factor;
               }
          }
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported type = %d", __func__, type);
        return -1;
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
                          const double *factor,
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

   if (factor != NULL)
     {
        if (0 != apply_scale_factor (*factor, info.type, buf->bytes, num_elem))
          return -1;
     }

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
   int (*annotate)(int);
   const double *factor;
}
Var_Name_Type;
#define VAR_TABLE_END {NULL,NULL,NULL,NULL}

static Var_Name_Type SMC_Vars[] =
{
   {TEMPO_VAR_SMADIT_SCANX, TEMPO_VAR_SMADIT_SCANX, NULL, NULL},
   {TEMPO_VAR_SMADIT_SCANY, TEMPO_VAR_SMADIT_SCANY, NULL, NULL},
   {TEMPO_VAR_SMADIT_RAWX, TEMPO_VAR_SMADIT_RAWX, NULL, NULL},
   {TEMPO_VAR_SMADIT_RAWY, TEMPO_VAR_SMADIT_RAWY, NULL, NULL},
   VAR_TABLE_END
};

static int annotate_gyro_dqf (int grp);

static Var_Name_Type IRU_Vars[] =
{
   {TEMPO_VAR_GYRO_OUTPUT, TEMPO_VAR_GYRO_OUTPUT, NULL, NULL},
   {TEMPO_VAR_GYRO_DQF, TEMPO_VAR_GYRO_DQF, annotate_gyro_dqf, NULL},
   VAR_TABLE_END
};

static double _pGyro_Bias_Convert_Microradian_to_Radian = 1.e-6;
/* _pGyro_Bias_Convert_Microradian_to_Radian is the conversion factor:
 *    (microrad/sec)*(1.e-6 rad/microrad) -> (radian/sec)
 * Incoming IOC data stream provides gyro_bias values in microrad/sec.
 * Downstream, the INR software expects radian/sec.
 */

static Var_Name_Type IRU_Bias_Vars[] =
{
   {TEMPO_VAR_TIME_GYRO_BIAS, TEMPO_VAR_TIME_GYRO_BIAS, NULL, NULL},
   {"bias", TEMPO_VAR_GYRO_BIAS, NULL, &_pGyro_Bias_Convert_Microradian_to_Radian},
   VAR_TABLE_END
};

static Var_Name_Type EPH_Vars[] =
{
   {"anc_satx", TEMPO_VAR_SAT_X, NULL, NULL},
   {"anc_saty", TEMPO_VAR_SAT_Y, NULL, NULL},
   {"anc_satz", TEMPO_VAR_SAT_Z, NULL, NULL},
   {"anc_satvx", TEMPO_VAR_SAT_VX, NULL, NULL},
   {"anc_satvy", TEMPO_VAR_SAT_VY, NULL, NULL},
   {"anc_satvz", TEMPO_VAR_SAT_VZ, NULL, NULL},
   VAR_TABLE_END
};

static int annotate_gyro_dqf (int grp)
{
   uint16_t flag_masks[2] = {IOCSDPC_IRU_DATA_INVALID, IOCSDPC_IRU_DATA_SATURATED};
   char flag_meanings[] = "invalid saturated";
   int varid, status, len, num_values = 2;

   status = nc_inq_varid (grp, TEMPO_VAR_GYRO_DQF, &varid);
   if (NC_NOERR != status)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: cannot find variable %s", __func__, TEMPO_VAR_GYRO_DQF);
        return -1;
     }

   len = strlen(flag_meanings) + 1;
   status = nc_put_att_text (grp, varid, "flag_meanings", len, flag_meanings);
   if (NC_NOERR != status)
     {
        tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining ushort attribute %s (%s)",
                     __func__, "flag_meanings", nc_strerror(status));
        return -1;
     }

   status = nc_put_att_ushort (grp, varid, "flag_masks", NC_USHORT, num_values, flag_masks);
   if (NC_NOERR != status)
     {
        tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining ushort attribute %s (%s)",
                     __func__, "flag_masks", nc_strerror(status));
        return -1;
     }

   return 0;
}

static int radiance_copy_vars (const Row_Select_Type *rst_head, TIO_Meta_Type *meta,
                               Var_Name_Type *v, const char *from_group_path,
                               int to_grp, int copy_all)
{
   const Row_Select_Type *rst;
   Buffer_Type buf = {0};
   int from_ncid = 0;
   int to_start = 0;
   int status = -1;

   if (v->annotate)
     {
        if (0 != v->annotate (to_grp))
          return -1;
     }

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
        rows_copied = copy_file_var (rst, &buf, v->from, from_grp, copy_all, v->factor,
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

static int convert_ioc_string_to_taix (const char *str, double *ptaix)
{
   int day, msec, usec;
   double taix;

   if (str == NULL)
     {
        fprintf (stderr, "%s: NULL string\n", __func__);
        return -1;
     }

   if (3 != sscanf (str, "d%5dm%8du%3d", &day, &msec, &usec))
     {
        fprintf (stderr, "*** Error: parsing timestamp: %s\n", str);
        return -1;
     }

   taix = day * 86400.0 + msec/1000.0 + usec/1.e6;
   if (ptaix) *ptaix = taix;

   return 0;
}

typedef struct
{
   double *eph_time;
   double *satx;
   double *saty;
   double *satz;
   double *satvx;
   double *satvy;
   double *satvz;
   int num;
   int num_alloc;
}
Eph_Predicted_Type;

typedef struct
{
   double eph_time;
   double satx, saty, satz;
   double satvx, satvy, satvz;
}
Eph_Point_Type;

static void free_ephem (Eph_Predicted_Type *pred)
{
   FREE(pred->eph_time);
   FREE(pred->satx);
   FREE(pred->saty);
   FREE(pred->satz);
   FREE(pred->satvx);
   FREE(pred->satvy);
   FREE(pred->satvz);
}

static int realloc_dbl (double **x, int new_num)
{
   double *new_x;
   if (NULL == (new_x = (double *)REALLOC (*x, new_num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }
   *x = new_x;
   return 0;
}

static int append_ephem_point (Eph_Predicted_Type *pred, const Eph_Point_Type *pt)
{
   int i;

   /* assumes pred->num, pred->num_alloc are both initialized to zero */
   if (pred->num == pred->num_alloc)
     {
        int new_num = (pred->num_alloc > 0) ? (pred->num_alloc*2) : 10;
        if ((0 != realloc_dbl (&pred->eph_time, new_num))
            || (0 != realloc_dbl (&pred->satx, new_num))
            || (0 != realloc_dbl (&pred->saty, new_num))
            || (0 != realloc_dbl (&pred->satz, new_num))
            || (0 != realloc_dbl (&pred->satvx, new_num))
            || (0 != realloc_dbl (&pred->satvy, new_num))
            || (0 != realloc_dbl (&pred->satvz, new_num)))
          return -1;
        pred->num_alloc = new_num;
     }

   i = pred->num;
   pred->eph_time[i] = pt->eph_time;
   pred->satx[i] = pt->satx;
   pred->saty[i] = pt->saty;
   pred->satz[i] = pt->satz;
   pred->satvx[i] = pt->satvx;
   pred->satvy[i] = pt->satvy;
   pred->satvz[i] = pt->satvz;
   pred->num++;

   return 0;
}

static int parse_ephem_point (char ***data, unsigned int row, Eph_Point_Type *pt)
{
   if ((1 != sscanf (data[0][row], "%lf", &pt->eph_time))
       ||(1 != sscanf (data[1][row], "%lf", &pt->satx))
       ||(1 != sscanf (data[2][row], "%lf", &pt->saty))
       ||(1 != sscanf (data[3][row], "%lf", &pt->satz))
       ||(1 != sscanf (data[4][row], "%lf", &pt->satvx))
       ||(1 != sscanf (data[5][row], "%lf", &pt->satvy))
       ||(1 != sscanf (data[6][row], "%lf", &pt->satvz)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: parsing ephemeris data", __func__);
        return -1;
     }
   return 0;
}

static int transform_velocity_vector (Eph_Point_Type *pt)
{
   double v_geo = 3.074666284127684;  /* km/sec */
   double phi = atan2 (pt->saty, pt->satx);

   /* While the HGS/IOC ICD says that the predicted ephemeris velocities
    * are to be in km/sec, we seem to be getting meters/sec, relative to
    * the idealized geostationary orbital station.
    * INR wants something slightly different, so we make the (apparently)
    * necessary changes here.
    */

   /* convert m/sec to km/sec */
   pt->satvx *= 1.e-3;
   pt->satvy *= 1.e-3;
   pt->satvz *= 1.e-3;

   /* Add geostationary orbital velocity:
    * x = R cos(phi)  -> dx/dt = - R sin(phi) dphi/dt = - v_geo * sin(phi)
    * y = R sin(phi)  -> dy/dt =   R cos(phi) dphi/dt =   v_geo * cos(phi)
    */
   pt->satvx = - v_geo * sin(phi);
   pt->satvy =   v_geo * cos(phi);

   return 0;
}

static int append_predicted_ephemeris (Eph_Predicted_Type *pred, double time_beg, double time_end,
                                       int enable_adjust_velocity, const char *path)
{
   IOCLib_String_Table_Type *tbl = NULL;
   const char *eph_columns[] = {"time","sat_x","sat_y","sat_z","sat_vx","sat_vy","sat_vz"};
   int num_columns = sizeof(eph_columns)/sizeof(*eph_columns);
   unsigned int row;
   int loaded_points, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: reading %s", __func__, path);

   if (NULL == (tbl = ioclib_csv_read_string_table (path, eph_columns, num_columns)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading CSV file: %s", __func__, path);
        goto return_status;
     }

   loaded_points = 0;
   for (row = 0; row < tbl->num_rows; row++)
     {
        Eph_Point_Type pt;
        char *s0 = tbl->data[0][row];
        if (*s0 == ':') continue;
        if (0 != parse_ephem_point (tbl->data, row, &pt))
          goto return_status;
        if (pt.eph_time < time_beg)
          continue;
        if (pt.eph_time > time_end)
          break;
        if (enable_adjust_velocity)
          {
             (void) transform_velocity_vector (&pt);
          }
        if (0 != append_ephem_point (pred, &pt))
          goto return_status;
        loaded_points++;
     }

   status = loaded_points;
return_status:
   ioclib_free_string_table (tbl);
   return status;
}

static int write_predicted_ephemeris (Radiance_Type *r, const Eph_Predicted_Type *pred)
{
   int start = 0;
   int count = pred->num;
   int grp;

   if ((0 != TIO_inq_grp (r->ncid, "/inr_input/ephemeris", &grp))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_TIME_EPHEM, &start, &count, NC_DOUBLE, pred->eph_time))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_SAT_X, &start, &count, NC_DOUBLE, pred->satx))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_SAT_Y, &start, &count, NC_DOUBLE, pred->saty))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_SAT_Z, &start, &count, NC_DOUBLE, pred->satz))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_SAT_VX, &start, &count, NC_DOUBLE, pred->satvx))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_SAT_VY, &start, &count, NC_DOUBLE, pred->satvy))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_SAT_VZ, &start, &count, NC_DOUBLE, pred->satvz)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing ephemeris to file: %s", __func__, r->file);
        return -1;
     }

   return 0;
}

static int parse_filename (const char *filename, double *taix)
{
   const char *basename = ioclib_basename (filename);
   const char *timestamp;

   if ((NULL == (timestamp = strchr (basename, 'd')))
       || (0 != convert_ioc_string_to_taix (timestamp, taix)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error parsing timestamp: %s", __func__, timestamp);
        return -1;
     }
   return 0;
}

int radiance_copy_eph_predicted (Radiance_Type *r, double time_beg, double time_end,
                                 int enable_adjust_velocity, TIO_Meta_Type *meta, const char *eph_dir)
{
   IOCLib_Glob_Type *g = NULL;
   Eph_Predicted_Type pred = {0};
   char *eph_glob = NULL;
   double taix, delta, delta_min;
   unsigned int i;
   int k, num_points, status = -1;

   (void) meta;

   tell_vwarn (0, "%s: looking for predicted ephemeris files in dir: %s", __func__, eph_dir);

   if (NULL == (eph_glob = ioclib_pathconcat (eph_dir, "tempo_d*ephemeris.csv")))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: ioclib_pathconcat failed", __func__);
        return -1;
     }

   if ((NULL == (g = ioclib_glob (eph_glob, 0)))
       || (g->num_files == 0))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: no files found using glob pattern: %s", __func__, eph_glob);
        goto return_status;
     }

   /* Pick the ephemeris file received just prior to the start of the observation.
    * The ephemeris files span 21 days, are normally updated weekly, and have only
    * relatively coarse time resolution, so there's little point in reading more
    * than one file, especially since this is the 2nd or 3rd fallback option,
    * and should be used only rarely.
    */

   time_beg -= 600.0; /* pad the time interval [sec] */

   k = -1;
   delta_min = DBL_MAX;

   for (i = 0; i < g->num_files; i++)
     {
        /* The first ephemeris point in the file is the same as the filename timestamp */
        if (0 != parse_filename (g->files[i], &taix))
          goto return_status;
        if (taix > time_end)
          break;
        delta = time_beg - taix;
        if ((0 <= delta) && (delta < delta_min))
          {
             delta_min = delta;
             k = i;
          }
     }

   if (k >= 0)
     {
        if ((num_points = append_predicted_ephemeris (&pred, time_beg, time_end, enable_adjust_velocity, g->files[k])) < 0)
          goto return_status;
     }
   else num_points = 0;

   if (num_points > 0)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: read %d ephemeris points from %s", __func__, num_points, g->files[k]);
        if (0 != write_predicted_ephemeris (r, &pred))
          goto return_status;
     }
   else
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: found no predicted ephemeris data", __func__);
        goto return_status;
     }

   status = 0;
return_status:
   free_ephem (&pred);
   ioclib_free (eph_glob);
   ioclib_glob_free (g);
   return status;
}

