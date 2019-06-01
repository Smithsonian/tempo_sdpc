/** @file smc.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Aug 2017
 *  @brief Process inertial reference unit (iru) files
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#define PRODUCT_TYPE_IRU "iru"

#define PROCESS_METHOD_PRIVATE_DATA \
   char *out_dirname; \
   char *out_basename; \
   char *archdir_path; \
   int ncid; \
   int processing_version; \
   double outfile_timestamp_start; \
   double outfile_timestamp_end; \
   double outfile_deltat_sec; \
   double gyro_bias_time; \
   unsigned int gyro_dimension; \
   unsigned int bias_dimension; \
   void *outbuf; \
   size_t outbuf_num_bytes; \
   int num_written; \
   int num_files;
#include "l0_format.h"

#define IRU_FILE_VAR_TIME         "time"
#define IRU_FILE_VAR_SCALE_FACTOR "scale_factor"
#define IRU_FILE_VAR_BIAS         "bias"
#define IRU_FILE_VAR_BIAS_TIME    "time_of_scale_or_bias_update"

#define IRU_BIAS_DIMLEN   3

typedef struct
{
   const char *varname;    /**< variable name in netcdf file */
   int type;               /**< variable data type in netcdf file */
   int ndims;
   const char *units;      /**< unit string */
   const char *description;  /**< description of variable */
}
IRU_Var_Type;
#define IRU_VARS_END {NULL,-1,-1,NULL,NULL}
#define IRU_VAR(vn,ndims,type,units,descr) {vn,type,ndims,units,descr}

static IRU_Var_Type IRU_Var_Table[] =
{
   IRU_VAR(IRU_FILE_VAR_TIME,1,NC_DOUBLE,"sec","Time of sample (secs since TEMPO epoch)"),
   IRU_VAR(TEMPO_VAR_GYRO_OUTPUT,2,NC_INT,"DN","Gyro output data values"),
   IRU_VARS_END
};

static int close_iru_outfile (Process_Method_Type *pmt)
{
   if (pmt->ncid != INT_MAX)
     {
        if (-1 == write_attr_global_timestamp (pmt->ncid, "time_coverage_end",
                                               pmt->outfile_timestamp_end))
          return -1;
        if (-1 == close_hidden (pmt->ncid, pmt->out_dirname, pmt->out_basename, pmt->archdir_path))
          return -1;
     }
   pmt->ncid = INT_MAX;
   pmt->num_written = 0;
   pmt->num_files = 0;
   return 0;
}

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo)
{
   (void) tpinfo;
   return close_iru_outfile (pmt);
}

static void delete_iru (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   (void) close_iru_outfile (pmt);
   ioclib_free (pmt->out_basename);
   FREE(pmt->out_dirname);
   FREE(pmt->outbuf);
   FREE(pmt->archdir_path);
   FREE(pmt);
}

static int query_latest_timestamp (Process_Method_Type *pmt, int notused, double *timestamp)
{
   (void) notused;
   *timestamp = pmt->outfile_timestamp_end;
   return 0;
}

static int define_iru_vars (Process_Method_Type *pmt,
                            const IOCSDPC_IRU_Type *iru)
{
   IRU_Var_Type *vt;
   int storage = NC_CHUNKED;
   size_t chunk_sizes[2];
   int varid, dimid_time, dimid_bias_time, dimid_gyro_axis, dimid_gyro_bias;
   int dimids_gyro_bias[2], dimids_gyro_output[2], dimids_gyro_scale[2];
   const char *descr_bias_time =
     "Time when gyroscope scale or bias was updated (secs since TEMPO epoch)";

   pmt->gyro_bias_time = iru->gyro_bias_time;
   pmt->gyro_dimension = iru->gyro_dimension ? iru->gyro_dimension : 4;

   if ((0 != write_attr_global_product_type (pmt->ncid, PRODUCT_TYPE_IRU))
       || (0 != write_attr_global_timestamp (pmt->ncid, "time_coverage_start",
                                             pmt->outfile_timestamp_start)))
     return -1;

   /* FIXME?:
    * The dimension dimid_bias_time is being used somewhat inconsistently.
    * Each iru file incoming from the IOC has a single value of
    *     time_of_scale_or_bias_update
    *     scale_factor[gyro_axis]
    *     bias[gyro_bias]
    * Instead, the instrument will provide a bias time series, which
    * the IOC may ultimately provide as a separate file.  For now,
    * the bias_time dimension is effectively a counter indicating the
    * number of IRU products from the IOC merged to make the iru file.
    * This configuration is likely to change when we complete NDAs
    * with the gyroscope vendors and finally learn what's in the IRU
    * packets.
    */

   /* dimensions */
   if (-1 == TIO_def_dim (pmt->ncid, IRU_FILE_VAR_TIME, NC_UNLIMITED, &dimid_time))
     return -1;
   if (-1 == TIO_def_dim (pmt->ncid, TEMPO_VAR_TIME_GYRO_BIAS, NC_UNLIMITED, &dimid_bias_time))
     return -1;
   if (-1 == TIO_def_dim (pmt->ncid, TEMPO_DIM_GYRO_AXIS, pmt->gyro_dimension, &dimid_gyro_axis))
     return -1;
   if (-1 == TIO_def_dim (pmt->ncid, TEMPO_VAR_GYRO_BIAS, pmt->bias_dimension, &dimid_gyro_bias))
     return -1;

   /* bias time */
   if (0 != TIO_def_var (pmt->ncid, IRU_FILE_VAR_BIAS_TIME, NC_DOUBLE, 1, &dimid_bias_time, &varid))
     return -1;
   if (0 != annotate_var (pmt->ncid, varid, descr_bias_time, "sec"))
     return -1;

   /* scale factor */
   dimids_gyro_scale[0] = dimid_bias_time;
   dimids_gyro_scale[1] = dimid_gyro_axis;
   if (0 != TIO_def_var (pmt->ncid, IRU_FILE_VAR_SCALE_FACTOR, NC_FLOAT, 2, dimids_gyro_scale, &varid))
     return -1;
   if (0 != annotate_var (pmt->ncid, varid, "gyro scale factor", NULL))
     return -1;

   /* gyro bias */
   dimids_gyro_bias[0] = dimid_bias_time;
   dimids_gyro_bias[1] = dimid_gyro_bias;
   if (0 != TIO_def_var (pmt->ncid, IRU_FILE_VAR_BIAS, NC_FLOAT, 2, dimids_gyro_bias, &varid))
     return -1;
   if (0 != annotate_var (pmt->ncid, varid, "gyro bias", "microrad/sec"))
     return -1;

   /* gyro output time series */
   dimids_gyro_output[0] = dimid_time;
   dimids_gyro_output[1] = dimid_gyro_axis;

   chunk_sizes[0] = 1024;
   chunk_sizes[1] = pmt->gyro_dimension;

   for (vt = IRU_Var_Table; vt->varname != NULL; vt++)
     {
        if (0 != TIO_def_var (pmt->ncid, vt->varname, vt->type, vt->ndims, dimids_gyro_output, &varid))
          return -1;
        if (0 != TIO_def_var_chunking (pmt->ncid, varid, storage, chunk_sizes))
          return -1;
        if (0 != annotate_var (pmt->ncid, varid, vt->description, vt->units))
          return -1;
     }

   return 0;
}

static int realloc_outbuf (Process_Method_Type *pmt, size_t num_bytes)
{
   if (pmt->outbuf == NULL)
     {
        if (NULL == (pmt->outbuf = MALLOC (num_bytes)))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return -1;
          }
        pmt->outbuf_num_bytes = num_bytes;
     }
   else if (pmt->outbuf_num_bytes < num_bytes)
     {
        void *outbuf = NULL;
        if (NULL == (outbuf = REALLOC (pmt->outbuf, num_bytes)))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return -1;
          }
        pmt->outbuf = outbuf;
        pmt->outbuf_num_bytes = num_bytes;
     }

   return 0;
}

static int write_avg_sample_freq_attr (Process_Method_Type *pmt,
                                       double *sample_time,
                                       unsigned int num_records)
{
   TIO_Var_Info_Type info;
   double t0, t1, average_sample_frequency_hz;
   double delta_sum = 0.0;
   unsigned int i, num_deltas = 0;

   if (0 != TIO_inq_var (pmt->ncid, IRU_FILE_VAR_TIME, &info))
     return -1;

   t0 = sample_time[0];
   for (i = 1; i < num_records; i++)
     {
        t1 = sample_time[i];
        if ((0 == isnan(t0)) && (0 == isnan(t1)))
          {
             delta_sum += (t1 - t0);
             num_deltas += 1;
          }
        t0 = t1;
     }
   if (delta_sum > 0.0)
     {
        average_sample_frequency_hz = num_deltas / delta_sum;
     }
   else average_sample_frequency_hz = 0.0;

   if (0 != TIO_put_att (pmt->ncid, info.varid, "average_sample_frequency_hz",
                         NC_DOUBLE, 1, &average_sample_frequency_hz))
     return -1;

   return 0;
}

static int write_iru_records (Process_Method_Type *pmt,
                              IOCSDPC_IRU_Type *iru,
                              const IOCSDPC_IRU_Record_Type *rec_array)
{
   unsigned int num_records = iru->num_records;
   double time_last_sample = rec_array[num_records-1].sample_time;
   double *sample_time;
   int start[2], count[2], start_per_file[2], count_per_file[2];
   unsigned int i, ng;
   int32_t *gyro_data;
   size_t num_bytes;

   if (time_last_sample > pmt->outfile_timestamp_end)
     pmt->outfile_timestamp_end = time_last_sample;

   num_bytes = num_records * IOCSDPC_IRU_RECORD_SIZE;
   if (0 != realloc_outbuf (pmt, num_bytes))
     return -1;

   start[0] = pmt->num_written;
   start[1] = 0;
   count[0] = num_records;
   count[1] = pmt->gyro_dimension;

   sample_time = (double *)pmt->outbuf;
   for (i = 0; i < num_records; i++)
     {
        sample_time[i] = rec_array[i].sample_time;
     }

   if (0 != write_avg_sample_freq_attr (pmt, sample_time, num_records))
     return -1;

   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_TIME, start, count, NC_DOUBLE, sample_time))
     return -1;

   ng = pmt->gyro_dimension;
   gyro_data = (int32_t *)pmt->outbuf;
   for (i = 0; i < num_records; i++)
     {
        memcpy ((char *)&gyro_data[ng*i], (char *)rec_array[i].gyro_data, ng*sizeof(int32_t));
     }

   if (0 != TIO_put_var_section (pmt->ncid, TEMPO_VAR_GYRO_OUTPUT, start, count, NC_INT, gyro_data))
     return -1;

   pmt->num_written += num_records;

   start_per_file[0] = pmt->num_files;
   start_per_file[1] = 0;
   count_per_file[0] = 1;

   count_per_file[1] = pmt->bias_dimension;
   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_BIAS_TIME, start_per_file, count_per_file,
                                 NC_DOUBLE, &iru->gyro_bias_time))
     return -1;

   count_per_file[1] = pmt->bias_dimension;
   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_BIAS, start_per_file, count_per_file,
                                 NC_FLOAT, &iru->gyro_bias))
     return -1;

   count_per_file[1] = pmt->gyro_dimension;
   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_SCALE_FACTOR, start_per_file, count_per_file,
                                 NC_FLOAT, &iru->gyro_scale_factor))
     return -1;

   pmt->num_files += 1;

   return 0;
}

static int new_iru_outfile (Process_Method_Type *pmt, double timestamp,
                            const IOCSDPC_Common_Header_Type *chdr)
{
   char basename[MAX_BASENAME_SIZE];

   pmt->outfile_timestamp_start = timestamp;
   pmt->outfile_timestamp_end = timestamp;

   if (0 != verify_epoch (chdr->epoch))
     return -1;

   FREE(pmt->archdir_path);
   pmt->archdir_path = NULL;
   if (0 != make_level0_archdir_path (&pmt->archdir_path, timestamp,
                                      pmt->processing_version, PRODUCT_TYPE_IRU))
     return -1;

   if (0 != make_level0_basename (basename, sizeof(basename), timestamp,
                                  pmt->processing_version, PRODUCT_TYPE_IRU, NULL))
     return -1;

   if (pmt->ncid != INT_MAX)
     {
        if (0 != close_iru_outfile (pmt))
          return -1;
     }

   tell_vinfo (0, "creating file %s/%s", pmt->out_dirname, basename);

   if ((-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
       || (-1 == write_std_global_metadata (pmt->ncid, chdr)))
     return -1;
   ioclib_free (pmt->out_basename);
   if (NULL == (pmt->out_basename = ioclib_strdup (basename)))
     return -1;

   return 0;
}

static int select_iru_outfile (Process_Method_Type *pmt,
                               const IOCSDPC_IRU_Type *iru, double timestamp)
{
   /* FIXME - support splitting iru file across multiple
    * netcdf files? */
   if ((pmt->ncid != INT_MAX)
       && (0 < pmt->outfile_timestamp_start)
       && ((timestamp - pmt->outfile_timestamp_start)
           < pmt->outfile_deltat_sec)
       && (iru->gyro_bias_time == pmt->gyro_bias_time))
     return 0;

   if (0 != new_iru_outfile (pmt, timestamp, &iru->common_header))
     return -1;

   return define_iru_vars (pmt, iru);
}

static int process_iru (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        const char *file)
{
   IOCSDPC_Common_Header_Type chdr;
   IOCSDPC_IRU_Type *iru = NULL;
   IOCSDPC_IRU_Record_Type *rec_array = NULL;
   unsigned int num_read;
   size_t rec_array_size;
   int fd;

   (void) tpinfo;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;

   if (NULL == (iru = iocsdpc_iru_fdopen_read (file, fd, &chdr)))
     goto return_status;

   rec_array_size = iru->num_records * IOCSDPC_IRU_RECORD_SIZE;
   if (NULL == (rec_array = (IOCSDPC_IRU_Record_Type *)MALLOC (rec_array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   num_read = 0;
   if ((0 != iocsdpc_iru_read (iru, rec_array, iru->num_records, &num_read))
       || (num_read != iru->num_records))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %d records from %s",
                     __func__, iru->num_records, file);
        goto return_status;
     }

   if (0 != select_iru_outfile (pmt, iru, rec_array[0].sample_time))
     goto return_status;

   if (0 != write_iru_records (pmt, iru, rec_array))
     goto return_status;

   iocsdpc_iru_close (iru);
   ioclib_fd_close (fd);
   FREE(rec_array);
   return 0;

return_status:
   iocsdpc_iru_close (iru);
   ioclib_fd_close (fd);
   FREE(rec_array);
   return -1;
}

static int parse_iru_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;
   const char *out_dirname;

   if (NULL == (s = config_lookup (cfg, "iru")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'iru' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &pmt->processing_version))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "outfile_deltat_sec", &pmt->outfile_deltat_sec))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'iru' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (pmt->out_dirname = expand_string (out_dirname)))
     return -1;

   return 0;
}

Process_Method_Type *init_iru_method (config_t *cfg)
{
   Process_Method_Type *pmt = NULL;

   if (NULL == (pmt = (Process_Method_Type *)MALLOC (sizeof *pmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)pmt, 0, sizeof *pmt);

   if (-1 == parse_iru_params (cfg, pmt))
     {
        delete_iru (pmt);
        return NULL;
     }

   pmt->pmt_process = process_iru;
   pmt->pmt_delete = delete_iru;
   pmt->pmt_flush_cache = flush_cache;
   pmt->pmt_query_latest_timestamp = query_latest_timestamp;

   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;
   pmt->bias_dimension = IRU_BIAS_DIMLEN;

   return pmt;
}
