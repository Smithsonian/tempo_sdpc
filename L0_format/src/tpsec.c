/** @file tpsec.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Process telemetry point section files
 */

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "enum.h"

#define PRODUCT_TYPE_TPSEC "hk"

#define PROCESS_METHOD_PRIVATE_DATA \
   Enum_Lookup_Type *enum_lookup; \
   char *out_dirname; \
   char *out_basename; \
   char *archdir_path; \
   int ncid; \
   int processing_version; \
   double outfile_timestamp_start; \
   double outfile_timestamp_end; \
   double outfile_deltat_sec;
#include "l0_format.h"

#define L0_ENUM_TABLE_SIZE 256
#define L0_TPSEC_CHUNKSIZE 256

typedef struct
{
   double *raw;
   double *eng;
   unsigned int *limit_flags;
}
TPSec_Col_Type;

static void free_tpsec_col_type (TPSec_Col_Type *ct)
{
   if (ct == NULL)
     return;
   FREE(ct->limit_flags);
   FREE(ct->raw);
   FREE(ct);
}

static TPSec_Col_Type *alloc_tpsec_col_type (int nrows)
{
   TPSec_Col_Type *ct = NULL;
   if (NULL == (ct = (TPSec_Col_Type *)MALLOC (sizeof *ct)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   if ((NULL == (ct->raw = (double *) MALLOC (2 * nrows * sizeof(double))))
       ||(NULL == (ct->limit_flags = (unsigned int *) MALLOC (nrows * sizeof(unsigned int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_tpsec_col_type (ct);
        return NULL;
     }
   ct->eng = ct->raw + nrows;
   return ct;
}

static int close_outfile (Process_Method_Type *pmt)
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
   elt_close (pmt->enum_lookup);
   pmt->enum_lookup = NULL;
   return 0;
}

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo)
{
   (void) tpinfo;
   return close_outfile (pmt);
}

/* Create a new netCDF output file, optionally closing the current one. */
static int new_outfile (Process_Method_Type *pmt,
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
   if (0 != make_level0_archdir_path (&pmt->archdir_path, timestamp, -1,
                                      pmt->processing_version, PRODUCT_TYPE_TPSEC))
     return -1;

   if (-1 == make_level0_basename (basename, sizeof(basename), timestamp,
                                   pmt->processing_version, PRODUCT_TYPE_TPSEC, NULL))
     return -1;

   if (pmt->ncid != INT_MAX)
     {
        if (-1 == close_outfile (pmt))
          return -1;
     }

   tell_vinfo (0, "creating file %s/%s", pmt->out_dirname, basename);

   if ((-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
       || (-1 == write_std_global_metadata (pmt->ncid, chdr)))
     return -1;
   ioclib_free (pmt->out_basename);
   if (NULL == (pmt->out_basename = ioclib_strdup (basename)))
     return -1;

   if ((0 != write_attr_global_product_type (pmt->ncid, PRODUCT_TYPE_TPSEC))
       || (-1 == write_attr_global_timestamp (pmt->ncid, "time_coverage_start",
                                              pmt->outfile_timestamp_start)))
     return -1;

   if (NULL == (pmt->enum_lookup = elt_open (L0_ENUM_TABLE_SIZE)))
     return -1;

   return 0;
}

static int l0_nctype_from_ioctype (int ioc_data_type, int *nc_data_type)
{
   int type;
   switch (ioc_data_type)
     {
      case IOCDB_DTYPE_I8:  type = NC_BYTE; break;
      case IOCDB_DTYPE_U8:  type = NC_UBYTE; break;
      case IOCDB_DTYPE_I16: type = NC_SHORT; break;
      case IOCDB_DTYPE_U16: type = NC_USHORT; break;
      case IOCDB_DTYPE_I32: type = NC_INT; break;
      case IOCDB_DTYPE_U32: type = NC_UINT; break;
      case IOCDB_DTYPE_I64: type = NC_INT64; break;
      case IOCDB_DTYPE_U64: type = NC_UINT64; break;
      case IOCDB_DTYPE_F32: type = NC_FLOAT; break;
      case IOCDB_DTYPE_F64: type = NC_DOUBLE; break;
      case IOCDB_DTYPE_ENUM:
        type = NC_INT; /* AKA "STATE", uint32 in FS */
        break;
      default:
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: unsupported IOC data type %d",
                     __func__, ioc_data_type);
        return -1;
     }

   *nc_data_type = type;

   return 0;
}

static int append_suffix (const char *str, const char *suffix,
                          char *buf, int bufsize)
{
   return snprintf (buf, bufsize, "%s%s", str, suffix);
}

static int create_group
(const Process_Method_Type *pmt, const IOCSDPC_TPSec_Type *s,
    const TPInfo_Type *tpinfo, const char *name, int *pgrp)
{
   int grp, dimid_time, varid_time;
   IOCSDPC_TPSec_Col_Info_Type *col_info;
   int shuffle=1, deflate=1, deflate_level=1;
   size_t chunksizes = L0_TPSEC_CHUNKSIZE;
   char tmpname[72];
   int i, ncols;

   if (-1 == TIO_def_grp (pmt->ncid, name, &grp))
     return -1;
   if (-1 == TIO_def_dim (grp, "time", NC_UNLIMITED, &dimid_time))
     return -1;

   if (-1 == TIO_def_var (grp, "time", NC_DOUBLE, 1, &dimid_time, &varid_time))
     return -1;
   if (-1 == TIO_def_var_deflate (grp, varid_time, shuffle, deflate, deflate_level))
     return -1;
   if (-1 == TIO_def_var_chunking (grp, varid_time, NC_CHUNKED, &chunksizes))
     return -1;
   col_info = s->col_info;
   ncols = s->num_cols;

   for (i = 0; i < ncols; i++)
     {
        TPFields_Type tpfields;
        int is_enum_type, varid_eng, type_raw;
        int var_typeid;

        if (-1 == l0_nctype_from_ioctype (col_info->data_type, &type_raw))
          return -1;

        /* Is this an enum? */
        memset ((char *)&tpfields, 0, sizeof(TPFields_Type));
        if (NULL != tpinfo_find (tpinfo, col_info->mnemonic, &tpfields))
          {
             is_enum_type = elt_is_valid (pmt->enum_lookup, tpfields.enumlist);
             if (tpfields.enumlist && (is_enum_type == 0))
               {
                  tell_vwarn (0, "%s: invalid enum %s, enumlist=%s",
                              __func__, col_info->mnemonic, tpfields.enumlist);
               }
          }
        else is_enum_type = 0;

        /* Determine the typeid */
        if (is_enum_type)
          {
             if (-1 == elt_define (pmt->enum_lookup, pmt->ncid, col_info->mnemonic,
                                   tpfields.enumlist, type_raw, &var_typeid))
               return -1;
          }
        else var_typeid = NC_DOUBLE;

        /* engineering units */
        if (-1 == TIO_def_var (grp, col_info->mnemonic, var_typeid, 1, &dimid_time, &varid_eng))
          return -1;
        if (-1 == annotate_var (grp, varid_eng, tpfields.synopsis, tpfields.units))
          return -1;
        if (-1 == TIO_def_var_deflate (grp, varid_eng, shuffle, deflate, deflate_level))
          return -1;
        if (-1 == TIO_def_var_chunking (grp, varid_eng, NC_CHUNKED, &chunksizes))
          return -1;

        if (is_enum_type == 0)
          {
             int varid_raw, varid_lim;

             /* raw */
             append_suffix (col_info->mnemonic, "_raw", tmpname, sizeof(tmpname));
             if (-1 == TIO_def_var (grp, tmpname, type_raw, 1, &dimid_time, &varid_raw))
               return -1;
             if (-1 == TIO_def_var_deflate (grp, varid_raw, shuffle, deflate, deflate_level))
               return -1;
             if (-1 == nc_def_var_chunking (grp, varid_raw, NC_CHUNKED, &chunksizes))
               return -1;

             /* limit_flags */
             append_suffix (col_info->mnemonic, "_lim", tmpname, sizeof(tmpname));
             if (-1 == TIO_def_var (grp, tmpname, NC_UINT, 1, &dimid_time, &varid_lim))
               return -1;
             if (-1 == TIO_def_var_deflate (grp, varid_lim, shuffle, deflate, deflate_level))
               return -1;
             if (-1 == TIO_def_var_chunking (grp, varid_lim, NC_CHUNKED, &chunksizes))
               return -1;
          }

        col_info++;
     }

   *pgrp = grp;

   return 0;
}

/* Ensure that the right netCDF output file is open, and
 * return the index of the correct group in that file.
 */
static int select_dest_group
(Process_Method_Type *pmt, const IOCSDPC_TPSec_Type *s,
    const TPInfo_Type *tpinfo, double timestamp, int *pgrp)
{
#define BUFSIZE 72
   char buf[BUFSIZE];
   int group_exists;

   /* FIXME - support splitting a sciextract tpsec file
    * across multiple netcdf files? */
   if ((pmt->ncid == INT_MAX)
       || (pmt->outfile_timestamp_start < 0)
       || ((timestamp - pmt->outfile_timestamp_start)
           > pmt->outfile_deltat_sec))
     {
        if (-1 == new_outfile (pmt, &s->common_header, timestamp))
          return -1;
     }

   snprintf (buf, sizeof(buf), "/%s", s->section_name);

   tell_push_queue();
   group_exists = (0 == TIO_inq_grp (pmt->ncid, buf, pgrp));
   tell_pop_queue (1);
   if (group_exists)
     {
        return 0;
     }

   return create_group (pmt, s, tpinfo, buf, pgrp);
}

static double *extract_timestamps (const IOCSDPC_TPSec_Type *s,
                                   IOCSDPC_TPSec_Row_Type **row_list)
{
   int i, nrows = s->num_rows;
   double *t = NULL;

   if (NULL == (t = (double *) MALLOC (nrows * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   for (i = 0; i < nrows; i++)
     {
        IOCSDPC_TPSec_Row_Type *rt = row_list[i];
        t[i] = rt->timestamp;
     }

   return t;
}

static int extract_column (const IOCSDPC_TPSec_Type *s,
                           IOCSDPC_TPSec_Row_Type **row_list,
                           unsigned int col_index, TPSec_Col_Type *ct)
{
   IOCSDPC_TPSec_Item_Type *item;
   int i, nrows = s->num_rows;

   if (col_index >= s->num_cols)
     return -1;

   for (i = 0; i < nrows; i++)
     {
        IOCSDPC_TPSec_Row_Type *rt = row_list[i];
        item = rt->data_items + col_index;
        ct->raw[i] = item->raw;
        ct->eng[i] = item->eng;
        ct->limit_flags[i] = item->limit_flags;
     }

   return 0;
}

static int write_tpsec_row_list
(Process_Method_Type *pmt, const IOCSDPC_TPSec_Type *s,
    const TPInfo_Type *tpinfo, int grp, IOCSDPC_TPSec_Row_Type **row_list)
{
   IOCSDPC_TPSec_Col_Info_Type *col_info;
   TPSec_Col_Type *ct = NULL;
   double *timestamp = NULL;
   int *enum_tmp_space = NULL;
   int i, k, nrows, ncols, start, count, time_dimid;
   size_t time_dimlen;
   char tmpname[72];

   (void) pmt;

   col_info = s->col_info;
   ncols = s->num_cols;
   nrows = s->num_rows;

   if (NULL == (ct = alloc_tpsec_col_type (nrows)))
     return -1;

   if (NULL == (timestamp = extract_timestamps (s, row_list)))
     goto return_status;

   if (timestamp[nrows-1] > pmt->outfile_timestamp_end)
     pmt->outfile_timestamp_end = timestamp[nrows-1];

   if (-1 == TIO_inq_dim (grp, "time", &time_dimid, &time_dimlen))
     goto return_status;

   start = time_dimlen;
   count = nrows;

   if (-1 == TIO_put_var_section (grp, "time", &start, &count, NC_DOUBLE, timestamp))
     goto return_status;

   for (i = 0; i < ncols; i++)
     {
        TPFields_Type tpfields;
        int is_enum_type, eng_var_typeid;
        void *eng_values = NULL;

        if (-1 == extract_column (s, row_list, i, ct))
          goto return_status;

        /* Is this an enum? */
        memset ((char *)&tpfields, 0, sizeof(TPFields_Type));
        if (NULL != tpinfo_find (tpinfo, col_info->mnemonic, &tpfields))
          {
             is_enum_type = elt_is_valid (pmt->enum_lookup, tpfields.enumlist);
          }
        else is_enum_type = 0;

        /* engineering units: */
        if (is_enum_type == 0)
          {
             eng_var_typeid = NC_DOUBLE;
             eng_values = ct->eng;
          }
        else
          {
             TIO_Var_Info_Type var_info;
             if (-1 == TIO_inq_var (grp, col_info->mnemonic, &var_info))
               goto return_status;
             eng_var_typeid = var_info.type;
             if (enum_tmp_space == NULL)
               {
                  if (NULL == (enum_tmp_space = (int *) MALLOC (count * sizeof(int))))
                    {
                       tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
                       goto return_status;
                    }
               }
             for (k = 0; k < count; k++)
               {
                  enum_tmp_space[k] = ct->eng[k];
               }
             eng_values = enum_tmp_space;
          }
        if (-1 == TIO_put_var_section (grp, col_info->mnemonic, &start, &count, eng_var_typeid, eng_values))
          goto return_status;

        if (is_enum_type == 0)
          {
             /* raw */
             append_suffix (col_info->mnemonic, "_raw", tmpname, sizeof(tmpname));
             if (-1 == TIO_put_var_section (grp, tmpname, &start, &count, NC_DOUBLE, ct->raw))
               goto return_status;

             /* limit_flags */
             append_suffix (col_info->mnemonic, "_lim", tmpname, sizeof(tmpname));
             if (-1 == TIO_put_var_section (grp, tmpname, &start, &count, NC_UINT, ct->limit_flags))
               goto return_status;
          }

        col_info++;
     }

   free_tpsec_col_type (ct);
   FREE(timestamp);
   FREE(enum_tmp_space);
   return 0;

return_status:
   free_tpsec_col_type (ct);
   FREE(timestamp);
   FREE(enum_tmp_space);
   return -1;
}

static int process_tpsec_row_list
(Process_Method_Type *pmt, const IOCSDPC_TPSec_Type *s,
    const TPInfo_Type *tpinfo, IOCSDPC_TPSec_Row_Type **row_list)
{
   int grp;

   if (-1 == select_dest_group (pmt, s, tpinfo, row_list[0]->timestamp, &grp))
     return -1;

   if (-1 == write_tpsec_row_list (pmt, s, tpinfo, grp, row_list))
     return -1;

   return 0;
}

static IOCSDPC_TPSec_Row_Type **alloc_tpsec_row_list (int nrows)
{
   IOCSDPC_TPSec_Row_Type **row_list = NULL;
   if (NULL == (row_list = (IOCSDPC_TPSec_Row_Type **) MALLOC (nrows * sizeof(*row_list))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)row_list, 0, nrows * sizeof(*row_list));
   return row_list;
}

static void free_tpsec_row_list (IOCSDPC_TPSec_Row_Type **row_list, int nrows)
{
   if (row_list == NULL)
     return;
   while (nrows-- > 0)
     {
        iocsdpc_tpsec_free_row_type (row_list[nrows]);
     }
   FREE(row_list);
}

static int process_tpsec_file
(Process_Method_Type *pmt, const TPInfo_Type *tpinfo, const char *file)
{
   IOCSDPC_Common_Header_Type chdr;
   IOCSDPC_TPSec_Type *s;
   IOCSDPC_TPSec_Row_Type **row_list = NULL;
   unsigned int i, nrows;
   int fd;

   (void) pmt;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;

   if (NULL == (s = iocsdpc_tpsec_fdopen_read (file, fd, &chdr)))
     goto return_error;

   nrows = s->num_rows;

   if (NULL == (row_list = alloc_tpsec_row_list (nrows)))
     goto return_error;

   for (i = 0; i < nrows; i++)
     {
        if (-1 == iocsdpc_tpsec_read_row (s, &row_list[i]))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading row %u of %s",
                          __func__, i+1, file);
             goto return_error;
          }
     }

   if (-1 == process_tpsec_row_list (pmt, s, tpinfo, row_list))
     goto return_error;

   free_tpsec_row_list (row_list, nrows);
   iocsdpc_tpsec_close (s);
   ioclib_fd_close (fd);

   return tio_sync (pmt->ncid);

return_error:
   free_tpsec_row_list (row_list, nrows);
   iocsdpc_tpsec_close (s);
   ioclib_fd_close (fd);
   return -1;
}

static int parse_tpsec_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;
   const char *out_dirname;

   if (NULL == (s = config_lookup (cfg, "tpsec")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'tpsec' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &pmt->processing_version))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "outfile_deltat_sec", &pmt->outfile_deltat_sec)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'tpsec' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (pmt->out_dirname = expand_string (out_dirname)))
     return -1;

   return 0;
}

static void delete_tpsec (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   (void) close_outfile (pmt);
   ioclib_free (pmt->out_basename);
   FREE(pmt->out_dirname);
   FREE(pmt->archdir_path);
   FREE(pmt);
}

static int query_latest_timestamp (Process_Method_Type *pmt, int notused, double *timestamp)
{
   (void) notused;
   *timestamp = pmt->outfile_timestamp_end;
   return 0;
}

Process_Method_Type *init_tpsec_method (config_t *cfg)
{
   Process_Method_Type *pmt = NULL;

   if (NULL == (pmt = (Process_Method_Type *) MALLOC(sizeof *pmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)pmt, 0, sizeof *pmt);

   if (-1 == parse_tpsec_params (cfg, pmt))
     {
        delete_tpsec (pmt);
        return NULL;
     }

   pmt->pmt_process = process_tpsec_file;
   pmt->pmt_delete = delete_tpsec;
   pmt->pmt_flush_cache = flush_cache;
   pmt->pmt_query_latest_timestamp = query_latest_timestamp;

   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;
   pmt->enum_lookup = NULL;

   return pmt;
}
