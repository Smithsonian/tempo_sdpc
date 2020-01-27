/** @file smc.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Aug 2017
 *  @brief Process scan mechanism controller (smc) files
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tio.h>
#include <tio_template.h>
#include <tell.h>

#define PROCESS_METHOD_PRIVATE_DATA \
   char *out_dirname; \
   char *out_basename; \
   char *archdir_path; \
   int ncid; \
   int processing_version; \
   double outfile_timestamp_start; \
   double outfile_timestamp_end; \
   double outfile_deltat_sec; \
   void *outbuf; \
   size_t outbuf_num_bytes; \
   int num_written;
#include "l0_format.h"

typedef struct
{
   const char *varname;    /**< variable name in netcdf file */
   unsigned int offset;    /**< offset in bytes of relevant field in IOCSDPC_SMC_Record_Type */
   int type;               /**< variable data type in netcdf file */
   int varid;              /**< variable id in netcdf file */
   const char *units;      /**< unit string */
   const char *description;  /**< description of variable */
}
SMC_Var_Type;
#define SMC_VARS_END {NULL,0,-1,-1,NULL,NULL}
#define SMC_VAR1(field,vn,type,units,descr) \
                 {vn,offsetof(IOCSDPC_SMC_Record_Type,field),type,-1,units,descr}
#define SMC_VAR(field,type,units,descr) SMC_VAR1(field,#field,type,units,descr)

static SMC_Var_Type SMC_Var_Table[] =
{
   SMC_VAR1(sample_time,"time",NC_DOUBLE,"sec",
            "Time of sample (secs since TEMPO epoch)"),
   SMC_VAR(dit_raw_x,NC_FLOAT,NULL,
          "Filtered raw DIT readout for sensor axis-x"),
   SMC_VAR(dit_raw_y,NC_FLOAT,NULL,
          "Filtered raw DIT readout for sensor axis-y"),
   SMC_VAR(proc_meas_x,NC_FLOAT,"urad",
          "Filtered estimated control axis-x scan pos"),
   SMC_VAR(proc_meas_y,NC_FLOAT,"urad",
          "Filtered estimated control axis-y scan pos"),
   SMC_VAR(proc_meas_x_err,NC_FLOAT,"urad",
          "Error in the proc_meas_x value"),
   SMC_VAR(proc_meas_y_err,NC_FLOAT,"urad",
          "Error in the proc_meas_y value"),
   SMC_VAR(x_motor_current,NC_FLOAT,"TBD",
          "motor current (X dir)"),
   SMC_VAR(y_motor_current,NC_FLOAT,"TBD",
          "motor current (Y dir)"),
   SMC_VARS_END
};

static int close_smc_outfile (Process_Method_Type *pmt)
{
   if (pmt->ncid != INT_MAX)
     {
        if (-1 == write_attr_global_timestamp (pmt->ncid, "time_coverage_end",
                                               pmt->outfile_timestamp_end))
          return -1;

        /* close the file */
        if (0 != TIO_close (pmt->ncid))
          return -1;
        pmt->ncid = INT_MAX;

        if (pmt->archdir_path)
          {
             /* Put a copy in the archive and remove the original */
             if (0 != copy_hidden (pmt->out_dirname, pmt->out_basename, pmt->archdir_path))
               return -1;
             if (0 != remove_hidden (pmt->out_dirname, pmt->out_basename))
               return -1;
          }
        else
          {
             /* If we're not archiving, un-hide the original file */
             if (0 != rename_hidden (pmt->out_dirname, pmt->out_basename))
               return -1;
          }
     }
   pmt->ncid = INT_MAX;
   pmt->num_written = 0;
   return 0;
}

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo)
{
   (void) tpinfo;
   return close_smc_outfile (pmt);
}

static void delete_smc (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   (void) close_smc_outfile (pmt);
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

static int define_smc_vars (Process_Method_Type *pmt)
{
   SMC_Var_Type *vt;
   int storage = NC_CHUNKED;
   size_t chunk_size = 1024;
   int dimid_time;

   if ((0 != TIO_label_product (pmt->ncid, TEMPO_PROD_TYPE_SMC, pmt->processing_version))
       || (0 != write_attr_global_timestamp (pmt->ncid, "time_coverage_start",
                                             pmt->outfile_timestamp_start)))
     return -1;

   if (-1 == TIO_def_dim (pmt->ncid, "time", NC_UNLIMITED, &dimid_time))
     return -1;

   for (vt = SMC_Var_Table; vt->varname != NULL; vt++)
     {
        if (0 != TIO_def_var (pmt->ncid, vt->varname, vt->type, 1, &dimid_time, &vt->varid))
          return -1;
        if (0 != TIO_def_var_chunking (pmt->ncid, vt->varid, storage, &chunk_size))
          return -1;
        if (0 != annotate_var (pmt->ncid, vt->varid, vt->description, vt->units))
          return -1;
     }

   return 0;
}

static int extract_float_field (const IOCSDPC_SMC_Record_Type *rec_array,
                                unsigned int num_values,
                                unsigned int offset, float *values)
{
   unsigned int i;

   for (i = 0; i < num_values; i++)
     {
        const IOCSDPC_SMC_Record_Type *rec_i = &rec_array[i];
        values[i] = *(const float *)((char *)rec_i + offset);
     }

   return 0;
}

static int extract_double_field (const IOCSDPC_SMC_Record_Type *rec_array,
                                 unsigned int num_values,
                                 unsigned int offset, double *values)
{
   unsigned int i;

   for (i = 0; i < num_values; i++)
     {
        const IOCSDPC_SMC_Record_Type *rec_i = &rec_array[i];
        values[i] = *(const double *)((char *)rec_i + offset);
     }

   return 0;
}

static int write_smc_field (const SMC_Var_Type *vt,
                            Process_Method_Type *pmt,
                            unsigned int num_records,
                            const IOCSDPC_SMC_Record_Type *rec_array)
{
   int start, count;

   switch (vt->type)
     {
      case NC_FLOAT:
        if (0 != extract_float_field (rec_array, num_records, vt->offset, pmt->outbuf))
          return -1;
        break;
      case NC_DOUBLE:
        if (0 != extract_double_field (rec_array, num_records, vt->offset, pmt->outbuf))
          return -1;
        break;
      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported field type = %d",
                     __func__, vt->type);
        return -1;
        break;
     }

   start = pmt->num_written;
   count = num_records;

   return TIO_put_var_section (pmt->ncid, vt->varname, &start, &count, vt->type, pmt->outbuf);
}

static int realloc_outbuf (Process_Method_Type *pmt, unsigned int num_records)
{
   size_t num_bytes = num_records * sizeof(double);

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
                                       unsigned int num_records,
                                       const IOCSDPC_SMC_Record_Type *rec_array)
{
   TIO_Var_Info_Type info;
   double t0, t1, average_sample_frequency_hz;
   double delta_sum = 0.0;
   unsigned int i, num_deltas = 0;

   if (0 != TIO_inq_var (pmt->ncid, "time", &info))
     return -1;

   t0 = rec_array[0].sample_time;
   for (i = 1; i < num_records; i++)
     {
        t1 = rec_array[i].sample_time;
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

static int write_smc_records (Process_Method_Type *pmt,
                              unsigned int num_records,
                              const IOCSDPC_SMC_Record_Type *rec_array)
{
   double time_last_sample = rec_array[num_records-1].sample_time;
   SMC_Var_Type *vt;

   if (time_last_sample > pmt->outfile_timestamp_end)
     pmt->outfile_timestamp_end = time_last_sample;

   if (0 != realloc_outbuf (pmt, num_records))
     return -1;

   for (vt = SMC_Var_Table; vt->varname != NULL; vt++)
     {
        if (0 != write_smc_field (vt, pmt, num_records, rec_array))
          return -1;
     }

   pmt->num_written += num_records;

   return 0;
}

static int new_smc_outfile (Process_Method_Type *pmt,
                            const IOCSDPC_Common_Header_Type *chdr,
                            double timestamp)
{
   char basename[MAX_BASENAME_SIZE];

   pmt->outfile_timestamp_start = timestamp;
   pmt->outfile_timestamp_end = timestamp;

   if (0 != verify_epoch (chdr->epoch))
     return -1;

   FREE(pmt->archdir_path);
   pmt->archdir_path = NULL;
   if (0 != make_level0_archdir_path (&pmt->archdir_path, timestamp, -1, TEMPO_PROD_TYPE_SMC))
     return -1;

   if (0 != make_level0_basename (basename, sizeof(basename), timestamp,
                                  pmt->processing_version, TEMPO_PROD_TYPE_SMC, NULL))
     return -1;

   if (pmt->ncid != INT_MAX)
     {
        if (0 != close_smc_outfile (pmt))
          return -1;
     }

   if ((-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
       || (-1 == write_std_global_metadata (pmt->ncid, chdr)))
     return -1;
   ioclib_free (pmt->out_basename);
   if (NULL == (pmt->out_basename = ioclib_strdup (basename)))
     return -1;

   return 0;
}

static int select_smc_outfile (Process_Method_Type *pmt,
                               const IOCSDPC_Common_Header_Type *chdr,
                               double timestamp)
{
   /* FIXME - support splitting smc file across multiple
    * netcdf files? */
   if ((pmt->ncid != INT_MAX)
       && (0 < pmt->outfile_timestamp_start)
       && ((timestamp - pmt->outfile_timestamp_start)
           < pmt->outfile_deltat_sec))
     return 0;

   if (0 != new_smc_outfile (pmt, chdr, timestamp))
     return -1;

   return define_smc_vars (pmt);
}

static int process_smc (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        const char *file, void *client_data)
{
   IOCSDPC_Common_Header_Type chdr;
   IOCSDPC_SMC_Type *smc = NULL;
   IOCSDPC_SMC_Record_Type *rec_array = NULL;
   unsigned int num_read;
   size_t rec_array_size;
   int fd;

   (void) tpinfo; (void) client_data;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     {
        tell_vlog (TELL_MSGTYPE_ERROR, 0, "%s: opening SMC file: %s", __func__, file);
        return -1;
     }

   if (NULL == (smc = iocsdpc_smc_fdopen_read (file, fd, &chdr)))
     goto return_status;

   rec_array_size = smc->num_records * IOCSDPC_SMC_RECORD_SIZE;
   if (NULL == (rec_array = (IOCSDPC_SMC_Record_Type *)MALLOC (rec_array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   num_read = 0;
   if ((0 != iocsdpc_smc_read (smc, rec_array, smc->num_records, &num_read))
       || (num_read != smc->num_records))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %d records from %s",
                     __func__, smc->num_records, file);
        goto return_status;
     }

   if (0 != select_smc_outfile (pmt, &smc->common_header, rec_array[0].sample_time))
     goto return_status;

   if (0 != write_smc_records (pmt, smc->num_records, rec_array))
     goto return_status;

   if (0 != write_avg_sample_freq_attr (pmt, smc->num_records, rec_array))
     goto return_status;

   iocsdpc_smc_close (smc);
   ioclib_fd_close (fd);
   FREE(rec_array);

   return tio_sync (pmt->ncid);

return_status:
   tell_vlog (TELL_MSGTYPE_ERROR, 0, "%s: processing SMC file: %s", __func__, file);
   iocsdpc_smc_close (smc);
   ioclib_fd_close (fd);
   FREE(rec_array);
   return -1;
}

static int parse_smc_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;
   const char *out_dirname;

   if (NULL == (s = config_lookup (cfg, "smc")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'smc' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "outfile_deltat_sec", &pmt->outfile_deltat_sec))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'smc' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((pmt->processing_version = get_processing_version()) < 0)
     return -1;

   if (NULL == (pmt->out_dirname = expand_string (out_dirname)))
     return -1;

   return 0;
}

Process_Method_Type *init_smc_method (config_t *cfg)
{
   Process_Method_Type *pmt = NULL;

   if (NULL == (pmt = (Process_Method_Type *)MALLOC (sizeof *pmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)pmt, 0, sizeof *pmt);

   if (-1 == parse_smc_params (cfg, pmt))
     {
        delete_smc (pmt);
        return NULL;
     }

   pmt->pmt_process = process_smc;
   pmt->pmt_delete = delete_smc;
   pmt->pmt_flush_cache = flush_cache;
   pmt->pmt_query_latest_timestamp = query_latest_timestamp;

   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;

   return pmt;
}
