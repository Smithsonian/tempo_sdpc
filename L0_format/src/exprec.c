/** @file exprec.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Process exposure record files
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tio.h>
#include <tell.h>

#define PROCESS_METHOD_PRIVATE_DATA \
   Enum_Lookup_Type *enum_lookup; \
   const char *out_dirname; \
   char *out_basename; \
   int ncid; \
   int processing_version; \
   double outfile_timestamp_start; \
   double outfile_timestamp_end; \
   char *exprec_type_string; \
   int exprec_type; \
   int outfile_erec_capacity; \
   int outfile_erec_count;
#include "l0_format.h"

#define EXPREC_ENUM_TABLE_SIZE 16

static int define_exprec_vars (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                               IOCSDPC_Exprec_Type *erec)
{
   const IOCSDPC_Image_Info_Item_Type *item;
   int ncid = pmt->ncid;
   int num_erec = pmt->outfile_erec_capacity;
   int shuffle=1, deflate=1, deflate_level=1;
   int dimid_time, dimid_row, dimid_col;
   int varid_image_start_time, varid_exposure_time, varid_readout_time, varid_frame_transfer_time;
   int varid_exprec, len;
   unsigned int nth, value;
   int dimids_exprec[3];
   size_t chunksizes[3];

   len = strlen(pmt->exprec_type_string)+1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "exprec_type",
                          NC_CHAR, len, pmt->exprec_type_string))
     return -1;
   if (-1 == write_attr_global_timestamp (pmt->ncid, "time_coverage_start",
                                          pmt->outfile_timestamp_start))
     return -1;

   if (-1 == TIO_def_dim (ncid, "time", num_erec, &dimid_time))
     return -1;
   if (-1 == TIO_def_dim (ncid, "row", erec->num_rows, &dimid_row))
     return -1;
   if (-1 == TIO_def_dim (ncid, "col", erec->num_cols, &dimid_col))
     return -1;

   if (-1 == TIO_def_var (ncid, "image_start_time", NC_DOUBLE, 1, &dimid_time, &varid_image_start_time))
     return -1;
   if (-1 == TIO_def_var (ncid, "exposure_time", NC_DOUBLE, 1, &dimid_time, &varid_exposure_time))
     return -1;
   if (-1 == TIO_def_var (ncid, "readout_time", NC_DOUBLE, 1, &dimid_time, &varid_readout_time))
     return -1;
   if (-1 == TIO_def_var (ncid, "frame_transfer_time", NC_DOUBLE, 1, &dimid_time, &varid_frame_transfer_time))
     return -1;

   nth = 0;
   while (NULL != (item = iocsdpc_image_info_get_value_by_index (erec, nth, &value)))
     {
        TPFields_Type tpfields;
        int varid_item, is_enum_type;

        /* Is this an enum? */
        memset ((char *)&tpfields, 0, sizeof(TPFields_Type));
        if (NULL != tpinfo_find (tpinfo, item->mnemonic, &tpfields))
          {
             is_enum_type = elt_is_valid (pmt->enum_lookup, tpfields.enumlist);
          }
        else is_enum_type = 0;

        /* Determine the typeid */
        if (is_enum_type)
          {
             if (-1 == elt_define (pmt->enum_lookup, ncid, item->mnemonic,
                                   tpfields.enumlist, NC_INT, &varid_item))
               return -1;
          }
        else varid_item = NC_UINT;

        if (-1 == TIO_def_var (ncid, item->mnemonic, varid_item, 1, &dimid_time, &varid_item))
          return -1;
        if (-1 == annotate_var (ncid, varid_item, tpfields.synopsis, tpfields.units))
          return -1;
        nth++;
     }

   dimids_exprec[0] = dimid_time;
   dimids_exprec[1] = dimid_row;
   dimids_exprec[2] = dimid_col;
   if (-1 == TIO_def_var (ncid, "image", NC_INT, 3, dimids_exprec, &varid_exprec))
     return -1;
   if (-1 == TIO_def_var_deflate (ncid, varid_exprec, shuffle, deflate, deflate_level))
     return -1;
   chunksizes[0] = 1;
   chunksizes[1] = 10;
   chunksizes[2] = erec->num_cols;
   if (-1 == TIO_def_var_chunking (ncid, varid_exprec, NC_CHUNKED, chunksizes))
     return -1;

   return 0;
}

static int close_outfile (Process_Method_Type *pmt)
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
   elt_close (pmt->enum_lookup);
   pmt->enum_lookup = NULL;
   return 0;
}

/* Create a new netCDF output file, optionally closing the current one. */
static int new_outfile (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        IOCSDPC_Exprec_Type *erec)
{
   char basename[MAX_BASENAME_SIZE];
   char *exprec_type_suffix;
   char *path = NULL;

   pmt->outfile_timestamp_start = erec->image_start_time;
   pmt->outfile_timestamp_end = erec->image_start_time;
   pmt->exprec_type = erec->exprec_type;
   pmt->outfile_erec_count = 0;

   ioclib_free (pmt->exprec_type_string);
   pmt->exprec_type_string = NULL;

   switch (pmt->exprec_type)
     {
      case IOCSDPC_EXPREC_TYPE_RADIANCE:
        exprec_type_suffix = "rad0";
        pmt->exprec_type_string = ioclib_strdup("radiance");
        break;
      case IOCSDPC_EXPREC_TYPE_DARK:
        exprec_type_suffix = "drk0";
        pmt->exprec_type_string = ioclib_strdup("dark");
        break;
      case IOCSDPC_EXPREC_TYPE_IRRADIANCE:
        exprec_type_suffix = "irr0";
        pmt->exprec_type_string = ioclib_strdup("irradiance");
        break;
      case IOCSDPC_EXPREC_TYPE_LIN_IRR:
        exprec_type_suffix = "irrlin0";
        pmt->exprec_type_string = ioclib_strdup("irradiance,linearity");
        break;
      case IOCSDPC_EXPREC_TYPE_LIN_DARK:
        exprec_type_suffix = "drklin0";
        pmt->exprec_type_string = ioclib_strdup("dark,linearity");
        break;
      case IOCSDPC_EXPREC_TYPE_UNKNOWN:
        /* drop */
      default:
        exprec_type_suffix = "unk0";
        pmt->exprec_type_string = ioclib_strdup("unknown");
        break;
     }

   /* if ioclib_strdup failed, it already produced a log message */
   if (pmt->exprec_type_string == NULL)
     return -1;

   if (-1 == make_level0_basename (erec->image_start_time, pmt->processing_version,
                                   exprec_type_suffix, basename, sizeof(basename)))
     return -1;

   if (pmt->ncid != INT_MAX)
     {
        if (-1 == close_outfile (pmt))
          return -1;
     }

   if (-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
     return -1;
   FREE(pmt->out_basename);
   if (NULL == (pmt->out_basename = ioclib_strdup (basename)))
     return -1;

   if (NULL == (pmt->enum_lookup = elt_open (EXPREC_ENUM_TABLE_SIZE)))
     goto return_status;

   return define_exprec_vars (pmt, tpinfo, erec);

return_status:
   FREE(path);
   return -1;
}

static int select_outfile (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                           IOCSDPC_Exprec_Type *erec)
{
   /* A new file will be opened when either:
    *  a) a new exposure record type is encountered, or
    *  b) when the output file max capacity reached.
    * FIXME? - should the output file also closed after a timeout?
    */
   if ((pmt->ncid != INT_MAX)
       && (pmt->exprec_type == erec->exprec_type)
       && (pmt->outfile_erec_count < pmt->outfile_erec_capacity))
     return 0;

   return new_outfile (pmt, tpinfo, erec);
}

static int write_exprec (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                         IOCSDPC_Exprec_Type *erec)
{
   const IOCSDPC_Image_Info_Item_Type *item;
   int32_t *data = NULL;
   unsigned int nth, value;
   int ncid, start[3], count[3];

   ncid = pmt->ncid;

   start[0] = pmt->outfile_erec_count;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = 0;
   count[2] = 0;

   if (erec->image_start_time > pmt->outfile_timestamp_end)
     pmt->outfile_timestamp_end = erec->image_start_time;

   if (-1 == TIO_put_var_section (ncid, "image_start_time", start, count,
                                  NC_DOUBLE, &erec->image_start_time))
     goto return_status;
   if (-1 == TIO_put_var_section (ncid, "exposure_time", start, count,
                                  NC_DOUBLE, &erec->exposure_time))
     goto return_status;
   if (-1 == TIO_put_var_section (ncid, "readout_time", start, count,
                                  NC_DOUBLE, &erec->readout_time))
     goto return_status;
   if (-1 == TIO_put_var_section (ncid, "frame_transfer_time", start, count,
                                  NC_DOUBLE, &erec->frame_transfer_time))
     goto return_status;

   nth = 0;
   while (NULL != (item = iocsdpc_image_info_get_value_by_index (erec, nth, &value)))
     {
        TPFields_Type tpfields;
        int is_enum_type, type_id;
        /* Is this an enum? */
        memset ((char *)&tpfields, 0, sizeof(TPFields_Type));
        if (NULL != tpinfo_find (tpinfo, item->mnemonic, &tpfields))
          {
             is_enum_type = elt_is_valid (pmt->enum_lookup, tpfields.enumlist);
          }
        else is_enum_type = 0;

        if (is_enum_type)
          {
             TIO_Var_Info_Type var_info;
             if (-1 == TIO_inq_var (ncid, item->mnemonic, &var_info))
               goto return_status;
             type_id = var_info.type;
          }
        else type_id = NC_INT;

        if (-1 == TIO_put_var_section (ncid, item->mnemonic, start, count, type_id, &value))
          goto return_status;
        nth++;
     }

   if (NULL == (data = iocsdpc_exprec_read_image (erec, NULL)))
     goto return_status;

   count[0] = 1;
   count[1] = erec->num_rows;
   count[2] = erec->num_cols;
   if (-1 == TIO_put_var_section (ncid, "image", start, count, NC_INT, data))
     goto return_status;

   ioclib_free (data);

   pmt->outfile_erec_count++;
   return 0;

return_status:
   ioclib_free (data);
   return -1;
}

static int process_exprec (Process_Method_Type *pmt,
                           const TPInfo_Type *tpinfo, const char *file,
                           int fd, IOCSDPC_Common_Header_Type *chdrp)
{
   IOCSDPC_Exprec_Type *erec;

   if (NULL == (erec = iocsdpc_exprec_fdopen_read (file, fd, chdrp)))
     return -1;

   if (-1 == select_outfile (pmt, tpinfo, erec))
     goto return_status;

   if (-1 == write_exprec (pmt, tpinfo, erec))
     goto return_status;

   iocsdpc_exprec_close (erec);
   return 0;

return_status:
   iocsdpc_exprec_close (erec);
   return -1;
}

static int parse_exprec_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "exprec")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'exprec' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &pmt->processing_version))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &pmt->out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "outfile_capacity", &pmt->outfile_erec_capacity))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'exprec' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static void delete_exprec (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   (void) close_outfile (pmt);
   ioclib_free (pmt->exprec_type_string);
   FREE(pmt->out_basename);
   FREE(pmt);
}

Process_Method_Type *init_exprec_method (config_t *cfg)
{
   Process_Method_Type *pmt = NULL;

   if (NULL == (pmt = (Process_Method_Type *)MALLOC (sizeof *pmt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if (-1 == parse_exprec_params (cfg, pmt))
     {
        delete_exprec (pmt);
        return NULL;
     }

   pmt->process = process_exprec;
   pmt->delete = delete_exprec;
   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;
   pmt->enum_lookup = NULL;

   pmt->exprec_type = INT_MAX;
   pmt->exprec_type_string = NULL;
   pmt->outfile_erec_count = 0;

   return pmt;
}
