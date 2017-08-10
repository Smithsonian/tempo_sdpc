/** @file smc.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Aug 2017
 *  @brief Process inertial reference unit (iru) files
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tio.h>
#include <tell.h>

#define PROCESS_METHOD_PRIVATE_DATA \
   const char *out_dirname; \
   char *out_basename; \
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
   int num_written;
#include "l0_format.h"

#define IRU_FILE_VAR_TIME         "time"
#define IRU_FILE_VAR_GYRO_OUTPUT  "gyro_output"
#define IRU_FILE_VAR_SCALE_FACTOR "scale_factor"
#define IRU_FILE_VAR_BIAS         "bias"
#define IRU_FILE_VAR_BIAS_TIME    "time_of_scale_or_bias_update"

#define IRU_FILE_DIM_GYRO_AXIS    "gyro_axis"
#define IRU_FILE_DIM_GYRO_BIAS    "gyro_bias"

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
   IRU_VAR(IRU_FILE_VAR_GYRO_OUTPUT,2,NC_INT,"DN","Gyro output data values"),
   IRU_VARS_END
};

static int close_iru_outfile (Process_Method_Type *pmt)
{
   if (pmt->ncid != INT_MAX)
     {
        if (-1 == write_attr_global_timestamp (pmt->ncid, "time_coverage_end",
                                               pmt->outfile_timestamp_end))
          return -1;
        if (-1 == close_hidden (pmt->ncid, pmt->out_dirname, pmt->out_basename))
          return -1;
     }
   pmt->ncid = INT_MAX;
   return 0;
}

static void delete_iru (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   (void) close_iru_outfile (pmt);
   ioclib_free (pmt->out_basename);
   FREE(pmt->outbuf);
   FREE(pmt);
}

static int define_iru_vars (Process_Method_Type *pmt,
                            const IOCSDPC_IRU_Type *iru)
{
   IRU_Var_Type *vt;
   int storage = NC_CHUNKED;
   size_t chunk_sizes[2];
   int start, count, varid;
   int dimid_time, dimid_gyro_axis, dimid_gyro_bias, dimids[2];
   const char *descr_bias_time =
     "Time when gyroscope scale or bias was updated (secs since TEMPO epoch)";

   pmt->gyro_bias_time = iru->gyro_bias_time;
   pmt->gyro_dimension = iru->gyro_dimension ? iru->gyro_dimension : 4;

   if (0 != write_attr_global_timestamp (pmt->ncid, "time_coverage_start",
                                         pmt->outfile_timestamp_start))
     return -1;

   /* dimensions */
   if (-1 == TIO_def_dim (pmt->ncid, IRU_FILE_VAR_TIME, NC_UNLIMITED, &dimid_time))
     return -1;
   if (-1 == TIO_def_dim (pmt->ncid, IRU_FILE_DIM_GYRO_AXIS, pmt->gyro_dimension, &dimid_gyro_axis))
     return -1;
   if (-1 == TIO_def_dim (pmt->ncid, IRU_FILE_DIM_GYRO_BIAS, pmt->bias_dimension, &dimid_gyro_bias))
     return -1;

   start = 0;

   /* bias time */
   if (0 != TIO_def_var (pmt->ncid, IRU_FILE_VAR_BIAS_TIME, NC_DOUBLE, 0, NULL, &varid))
     return -1;
   if (0 != annotate_var (pmt->ncid, varid, descr_bias_time, "sec"))
     return -1;

   count = 1;
   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_BIAS_TIME, &start, &count,
                                 NC_DOUBLE, &iru->gyro_bias_time))
     return -1;

   /* scale factor */
   if (0 != TIO_def_var (pmt->ncid, IRU_FILE_VAR_SCALE_FACTOR, NC_FLOAT, 1, &dimid_gyro_axis, &varid))
     return -1;
   if (0 != annotate_var (pmt->ncid, varid, "gyro scale factor", NULL))
     return -1;

   count = pmt->gyro_dimension;
   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_SCALE_FACTOR, &start, &count,
                                 NC_FLOAT, &iru->gyro_scale_factor))
     return -1;

   /* bias */
   if (0 != TIO_def_var (pmt->ncid, IRU_FILE_VAR_BIAS, NC_FLOAT, 1, &dimid_gyro_bias, &varid))
     return -1;
   if (0 != annotate_var (pmt->ncid, varid, "gyro bias", "microrad/sec"))
     return -1;

   count = pmt->bias_dimension;
   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_BIAS, &start, &count,
                                 NC_FLOAT, &iru->gyro_bias))
     return -1;

   /* gyro output time series */
   dimids[0] = dimid_time;
   dimids[1] = dimid_gyro_axis;

   chunk_sizes[0] = 1024;
   chunk_sizes[1] = pmt->gyro_dimension;

   for (vt = IRU_Var_Table; vt->varname != NULL; vt++)
     {
        if (0 != TIO_def_var (pmt->ncid, vt->varname, vt->type, vt->ndims, dimids, &varid))
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

static int write_iru_records (Process_Method_Type *pmt,
                              unsigned int num_records,
                              const IOCSDPC_IRU_Record_Type *rec_array)
{
   double time_last_sample = rec_array[num_records-1].sample_time;
   double *sample_time;
   int start[2], count[2];
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

   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_TIME, start, count, NC_DOUBLE, sample_time))
     return -1;

   ng = pmt->gyro_dimension;
   gyro_data = (int32_t *)pmt->outbuf;
   for (i = 0; i < num_records; i++)
     {
        memcpy ((char *)&gyro_data[ng*i], (char *)rec_array[i].gyro_data, ng*sizeof(int32_t));
     }

   if (0 != TIO_put_var_section (pmt->ncid, IRU_FILE_VAR_GYRO_OUTPUT, start, count, NC_INT, gyro_data))
     return -1;

   pmt->num_written += num_records;

   return 0;
}

static int new_iru_outfile (Process_Method_Type *pmt, double timestamp)
{
   char basename[MAX_BASENAME_SIZE];

   pmt->outfile_timestamp_start = timestamp;
   pmt->outfile_timestamp_end = timestamp;
   if (0 != make_level0_basename (timestamp, pmt->processing_version, "iru0",
                                  basename, sizeof(basename)))
     return -1;

   if (pmt->ncid != INT_MAX)
     {
        if (0 != close_iru_outfile (pmt))
          return -1;
     }

   if (-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
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
   if ((0 < pmt->outfile_timestamp_start)
       && ((timestamp - pmt->outfile_timestamp_start)
           < pmt->outfile_deltat_sec)
       && (iru->gyro_bias_time == pmt->gyro_bias_time))
     return 0;

   if (0 != new_iru_outfile (pmt, timestamp))
     return -1;

   return define_iru_vars (pmt, iru);
}

static int process_iru (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        const char *file, int fd, IOCSDPC_Common_Header_Type *chdrp)
{
   IOCSDPC_IRU_Type *iru = NULL;
   IOCSDPC_IRU_Record_Type *rec_array = NULL;
   unsigned int num_read;
   size_t rec_array_size;

   (void) tpinfo;

   if (NULL == (iru = iocsdpc_iru_fdopen_read (file, fd, chdrp)))
     return -1;

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

   if (0 != write_iru_records (pmt, iru->num_records, rec_array))
     goto return_status;

   iocsdpc_iru_close (iru);
   FREE(rec_array);
   return 0;

return_status:
   iocsdpc_iru_close (iru);
   FREE(rec_array);
   return -1;
}

static int parse_iru_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "iru")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'iru' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &pmt->processing_version))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &pmt->out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "outfile_deltat_sec", &pmt->outfile_deltat_sec))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'iru' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

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

   pmt->process = process_iru;
   pmt->delete = delete_iru;
   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;
   pmt->bias_dimension = IRU_BIAS_DIMLEN;

   return pmt;
}
