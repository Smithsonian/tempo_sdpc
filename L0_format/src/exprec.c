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

/* crude but generous out-of-range values for TEMPO */
#define MAX_NUM_MIRROR_STEPS 1500
#define MAX_NUM_GRANULES 16

typedef struct
{
   unsigned int granule_sizes[MAX_NUM_GRANULES];
   unsigned int granule_ubound[MAX_NUM_GRANULES];
   unsigned int num_granules;
   unsigned int num_recs;
   unsigned int size_default;
   unsigned int size_max;
   unsigned int curr_granule;
}
Granule_Schedule_Type;

typedef struct
{
   int exprec_type;
   unsigned int curr_mirror_step;
   unsigned int num_mirror_steps;
}
File_Info_Type;

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
   int granule_size; \
   unsigned int curr_mirror_step; \
   Granule_Schedule_Type radiance_sched; \
   const char *cache_dirname;
#include "l0_format.h"

#define EXPREC_ENUM_TABLE_SIZE 16

static int compute_granule_sizes (unsigned int num_recs,
                                  unsigned int num_granules,
                                  Granule_Schedule_Type *sched)
{
   unsigned int i, n, remainder, ubound;

   if (num_granules > MAX_NUM_GRANULES)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: internal buffer size exceeded", __func__);
        return -1;
     }

   n = num_recs / num_granules;
   if (n < 1) n = 1;
   remainder = num_recs - n * num_granules;

   memset ((char *)sched->granule_sizes, 0, MAX_NUM_GRANULES * sizeof(unsigned int));
   memset ((char *)sched->granule_ubound, 0, MAX_NUM_GRANULES * sizeof(unsigned int));

   for (i = 0; i < num_granules; i++)
     {
        sched->granule_sizes[i] = n;
     }
   for (i = 0; i < remainder; i++)
     {
        sched->granule_sizes[i] += 1;
     }

   ubound = 0;
   for (i = 0; i < num_granules; i++)
     {
        ubound += sched->granule_sizes[i];
        sched->granule_ubound[i] = ubound;
     }

   sched->num_granules = num_granules;
   sched->num_recs = num_recs;
   sched->curr_granule = 0;

   return 0;
}

static int schedule_granules (unsigned int num_recs,
                              Granule_Schedule_Type *sched)
{
   unsigned int num_granules;

   num_granules = num_recs / sched->size_default;
   if (num_granules < 1) num_granules = 1;

   for (;;)
     {
        if (0 != compute_granule_sizes (num_recs, num_granules, sched))
          return -1;
        if (sched->granule_sizes[0] < sched->size_max)
          break;
        num_granules++;
     }

   return 0;
}

static double image_end_time (IOCSDPC_Exprec_Type *erec)
{
   return (erec->image_start_time
           + erec->exposure_time
           + erec->frame_transfer_time);
}

static int define_file_vars (Process_Method_Type *pmt,
                             const TPInfo_Type *tpinfo,
                             IOCSDPC_Exprec_Type *erec)
{
   const IOCSDPC_Image_Info_Item_Type *item;
   int ncid = pmt->ncid;
   int shuffle=1, deflate=1, deflate_level=1;
   int dimid_time, dimid_row, dimid_col;
   int varid_image_start_time, varid_exposure_time;
   int varid_readout_time, varid_frame_transfer_time;
   int varid_exprec, len;
   unsigned int nth, value;
   int dimids_exprec[3];
   size_t chunksizes[3];

   len = strlen(pmt->exprec_type_string)+1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "exprec_type",
                          NC_CHAR, len, pmt->exprec_type_string))
     return -1;
   if (-1 == write_attr_global_timestamp (ncid, "time_coverage_start",
                                          pmt->outfile_timestamp_start))
     return -1;

   if (-1 == TIO_def_dim (ncid, "time", pmt->granule_size, &dimid_time))
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
   if (pmt->enum_lookup)
     {
        elt_close (pmt->enum_lookup);
        pmt->enum_lookup = NULL;
     }
   return 0;
}

/* Create a new netCDF output file  */
static int new_outfile (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        int granule_size, IOCSDPC_Exprec_Type *erec)
{
   char basename[MAX_BASENAME_SIZE];
   char *exprec_type_suffix;

   pmt->outfile_timestamp_start = erec->image_start_time;
   pmt->outfile_timestamp_end = image_end_time (erec);
   pmt->exprec_type = erec->exprec_type;

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

   tell_vinfo (0, "creating file %s/%s", pmt->out_dirname, basename);

   if (-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
     return -1;
   FREE(pmt->out_basename);
   if (NULL == (pmt->out_basename = ioclib_strdup (basename)))
     return -1;

   if (NULL == (pmt->enum_lookup = elt_open (EXPREC_ENUM_TABLE_SIZE)))
     return -1;

   pmt->granule_size = granule_size;

   return define_file_vars (pmt, tpinfo, erec);
}

static int write_exprec (Process_Method_Type *pmt,
                         const TPInfo_Type *tpinfo,
                         int erec_start, IOCSDPC_Exprec_Type *erec)
{
   const IOCSDPC_Image_Info_Item_Type *item;
   int32_t *data = NULL;
   unsigned int nth, value;
   int ncid, start[3], count[3];

   if (erec_start >= pmt->granule_size)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: start index = %d out of range (dimension size=%d)",
                     __func__, erec_start, pmt->granule_size);
        return -1;
     }

   ncid = pmt->ncid;

   start[0] = erec_start;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = 0;
   count[2] = 0;

   if (erec->image_start_time > pmt->outfile_timestamp_end)
     pmt->outfile_timestamp_end = image_end_time (erec);

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

   return 0;

return_status:
   ioclib_free (data);
   return -1;
}

static int process_cache (Process_Method_Type *pmt,
                          const TPInfo_Type *tpinfo,
                          int process_all_files)
{
   const char *cache_dirname = pmt->cache_dirname;
   Granule_Schedule_Type other_sched = {0};
   Granule_Schedule_Type *sched = NULL;
   IOCSDPC_Exprec_Type *erec = NULL;
   char **files = NULL;
   char *path = NULL;
   unsigned int file_erec_count, total_erec_count;
   size_t i, num_files, num_removed=0;
   int fd, exprec_type, status = -1;

   /* If the cache directory is empty, do nothing.
    * Otherwise, the cache directory contains cached exposure records
    * of some type, and no output file is open.
    * Depending on the value of process_all_files, the goal is to pack
    * either some or all of the cached exposure records into netcdf granules.
    */

   if (NULL == (files = ioclib_dir_list(cache_dirname, &num_files, IOCLIB_LISTDIR_SORT)))
     return -1;

   /* If the cache is empty, do nothing */
   if (num_files == 0)
     {
        ioclib_free (files);
        return 0;
     }

   tell_vinfo (2, "processing exprec cache: %ld files", num_files);

   exprec_type = pmt->exprec_type;

   /* The schedule tells us how the exposure records should be packed
    * into granules. For radiances, the header tells us how many scan
    * steps were commanded, so scheduling is easier.
    * For other types of exposure records, we derive a packing schedule
    * based on the number of files we actually received. In most cases,
    * that will be a single granule, but we'll also handle the case when
    * multiple files might be needed.
    */
   if (exprec_type == IOCSDPC_EXPREC_TYPE_RADIANCE)
     {
        sched = &pmt->radiance_sched;
     }
   else
     {
        /* struct copy: field defaults from radiance_sched */
        other_sched = pmt->radiance_sched;
        if (0 != schedule_granules (num_files, &other_sched))
          return -1;
        sched = &other_sched;
     }

   file_erec_count = 0;
   total_erec_count = 0;

   for (i = 0; i < num_files; i++)
     {
        IOCSDPC_Common_Header_Type chdr;
        unsigned int next_erec;
        int need_new_outfile;

        if (NULL == (path = ioclib_pathconcat (cache_dirname, files[i])))
          goto return_status;

        if (-1 == (fd = iocsdpc_open_file_read (path, 0, &chdr)))
          {
             ioclib_free (path);
             path = NULL;
             goto return_status;
          }

        if (NULL == (erec = iocsdpc_exprec_fdopen_read (path, fd, &chdr)))
          {
             ioclib_fd_close (fd);
             ioclib_free (path);
             path = NULL;
             goto return_status;
          }

        /* Assert: exprec_type == erec->exprec_type */

        if (exprec_type == IOCSDPC_EXPREC_TYPE_RADIANCE)
          {
             unsigned int curr_mirror_step;
             if (NULL == iocsdpc_image_info_get_value (erec, "curr_mirror_step",
                                                       &curr_mirror_step))
               return -1;
             next_erec = curr_mirror_step;
          }
        else next_erec = total_erec_count;

        need_new_outfile = ((pmt->ncid == INT_MAX)
                            || (file_erec_count == (unsigned int) pmt->granule_size)
                            || (next_erec >= sched->granule_ubound[sched->curr_granule]));

        if (need_new_outfile)
          {
             unsigned int k, granule_size;

             if (-1 == close_outfile (pmt))
               return -1;

             /* From granule schedule, find out which file
              * should hold the next exposure record.
              */
             for (k = sched->curr_granule; k < sched->num_granules; k++)
               {
                  if (next_erec < sched->granule_ubound[k])
                    break;
               }

             if (k >= sched->num_granules)
               {
                  tell_vinfo (0, "bad file: %s", path);
                  if (0 != ioclib_rename_to_bad_file (path))
                    {
                       tell_verror (TELL_APPLICATION_ERROR,
                                    "%s: ioclib_rename_to_bad_file, failed: file=%s",
                                    __func__, path);
                       ioclib_free (path);
                       continue;
                    }
               }

             sched->curr_granule = k;
             granule_size = sched->granule_sizes[k];

             /* We prefer to open a file only when we can fill it. */
             if ((num_files - i < granule_size)
                 && (process_all_files == 0))
               {
                  /* No need to open yet: wait for more data */
                  break;
               }

             if ((process_all_files != 0)
                 && (granule_size > num_files-i))
               {
                  /* Must open now: size the granule to hold what's here */
                  granule_size = num_files - i;
               }

             if (0 != new_outfile (pmt, tpinfo, granule_size, erec))
               goto return_status;
             file_erec_count = 0;
          }

        if (0 != write_exprec (pmt, tpinfo, file_erec_count, erec))
          goto return_status;
        file_erec_count++;
        total_erec_count++;

        iocsdpc_exprec_close (erec);
        erec = NULL;
        if (0 != ioclib_unlink (path))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: unlink failed: path=%s",
                          __func__, path);
          }
        num_removed++;
        ioclib_free (path);
        path = NULL;
     }

   status = 0;
return_status:
   tell_vinfo (2, "done processing exprec cache: i=%ld num_files=%ld num_removed=%ld (status=%d)",
               i, num_files, num_removed, status);

   ioclib_free (path);
   ioclib_string_array_free (files, num_files);
   if (erec) iocsdpc_exprec_close (erec);
   (void) close_outfile (pmt);
   return status;
}

static int cache_file (const char *cache_dirname, const char *file)
{
   const char *basename = ioclib_basename (file);
   char *cache_path;
   int status;

   if (NULL == (cache_path = ioclib_pathconcat (cache_dirname, basename)))
     return -1;

   status = ioclib_rename (file, cache_path);
   ioclib_free (cache_path);
   return status;
}

static int radiance_belongs_to_later_granule (const Process_Method_Type *pmt,
                                              const File_Info_Type *file_info)
{
   const Granule_Schedule_Type *sched = &pmt->radiance_sched;
   unsigned int step_ubound = sched->granule_ubound[sched->curr_granule];
   if (file_info->exprec_type != IOCSDPC_EXPREC_TYPE_RADIANCE)
     return 0;
   return (file_info->curr_mirror_step >= step_ubound);
}

static int classify_file (const char *file, File_Info_Type *info)
{
   IOCSDPC_Common_Header_Type chdr;
   IOCSDPC_Exprec_Type *erec;
   int fd;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;
   if (NULL == (erec = iocsdpc_exprec_fdopen_read (file, fd, &chdr)))
     {
        ioclib_fd_close (fd);
        return -1;
     }

   info->exprec_type = erec->exprec_type;

   if ((NULL == iocsdpc_image_info_get_value (erec, "curr_mirror_step",
                                             &info->curr_mirror_step))
       || (NULL == iocsdpc_image_info_get_value (erec, "num_mirror_steps",
                                                 &info->num_mirror_steps)))
     {
        iocsdpc_exprec_close (erec);
        return -1;
     }

   iocsdpc_exprec_close (erec);

   if (info->exprec_type != IOCSDPC_EXPREC_TYPE_RADIANCE)
     return 0;

   if ((info->num_mirror_steps == 0)
       || (info->num_mirror_steps > MAX_NUM_MIRROR_STEPS)
       || (info->curr_mirror_step >= info->num_mirror_steps))
     {
        tell_vlog (TELL_MSGTYPE_INFO, 0,
                   "%s: invalid radiance exposure record header: curr_mirror_step=%d num_mirror_steps=%d file=%s",
                   __func__, info->curr_mirror_step,
                   info->num_mirror_steps, file);
        return -1;
     }

   return 0;
}

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo)
{
   return process_cache (pmt, tpinfo, 1);
}

static int process_exprec (Process_Method_Type *pmt,
                           const TPInfo_Type *tpinfo, const char *file)
{
   File_Info_Type file_info = {0};
   int is_radiance, is_new_type, is_new_scan, is_new_scan_length;

   if (0 != classify_file (file, &file_info))
     return -1;

   /* Invariant:  the cache always contains only records
    * that belong together and could share the same granule.
    */

   is_new_type = (file_info.exprec_type != pmt->exprec_type);
   is_radiance = (file_info.exprec_type == IOCSDPC_EXPREC_TYPE_RADIANCE);
   is_new_scan = (is_radiance
                  && (file_info.curr_mirror_step < pmt->curr_mirror_step));

   /* The cache contains records of type pmt->exprec_type.
    * Process the cache before changing the value of pmt->exprec_type. */
   if (is_new_type || is_new_scan)
     {
        if (0 != process_cache (pmt, tpinfo, 1))
          return -1;
     }

   pmt->exprec_type = file_info.exprec_type;
   if (0 != cache_file (pmt->cache_dirname, file))
     return -1;

   /* Non-radiance exposure records are cached until something
    * triggers cache processing; either arrival of a different
    * exposure type, or expiration of a timer. */

   if (is_radiance == 0)
     return 0;

   /* For radiance exposure records, we know how many to expect,
    * so we generate granules according to a pre-planned schedule */

   pmt->curr_mirror_step = file_info.curr_mirror_step;

   is_new_scan_length = (file_info.num_mirror_steps
                         != pmt->radiance_sched.num_recs);

   if (is_new_scan || is_new_scan_length)
     {
        /* new schedule upon new scan resets sched->curr_granule to 0 */
        if (0 != schedule_granules (file_info.num_mirror_steps,
                                    &pmt->radiance_sched))
          return -1;
     }

   if (radiance_belongs_to_later_granule (pmt, &file_info))
     {
        return process_cache (pmt, tpinfo, 0);
     }

   return 0;
}

static int parse_exprec_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;
   int size_default, size_max;

   if (NULL == (s = config_lookup (cfg, "exprec")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'exprec' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &pmt->processing_version))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &pmt->out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "cache_dir", &pmt->cache_dirname))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "granule_size_default", &size_default))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "granule_size_max", &size_max)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'exprec' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   pmt->radiance_sched.size_default = size_default;
   pmt->radiance_sched.size_max = size_max;

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
   memset ((char *)pmt, 0, sizeof *pmt);

   if (-1 == parse_exprec_params (cfg, pmt))
     {
        delete_exprec (pmt);
        return NULL;
     }

   if (0 != ioclib_mkdir (pmt->cache_dirname, 0))
     {
        delete_exprec (pmt);
        return NULL;
     }

   pmt->process = process_exprec;
   pmt->delete = delete_exprec;
   pmt->flush_cache = flush_cache;

   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;
   pmt->enum_lookup = NULL;

   pmt->exprec_type = INT_MAX;
   pmt->exprec_type_string = NULL;
   pmt->curr_mirror_step = UINT_MAX;

   return pmt;
}
