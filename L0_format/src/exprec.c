/** @file exprec.c
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Process exposure record files
 */

#include "config.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include <ioclib.h>
#include <iocsdpc.h>
#include <tio.h>
#include <tio_template.h>
#include <tell.h>

#include "exprec_cache.h"

#ifdef ENABLE_IMAGE_COMPRESSION
# define LEVEL0_IMAGE_SHUFFLE_SELECT 1
# define LEVEL0_IMAGE_DEFLATE_SELECT 1
# define LEVEL0_IMAGE_DEFLATE_LEVEL  1
#else
# define LEVEL0_IMAGE_SHUFFLE_SELECT 0
# define LEVEL0_IMAGE_DEFLATE_SELECT 0
# define LEVEL0_IMAGE_DEFLATE_LEVEL  0
#endif

typedef struct
{
   unsigned int *granule_sizes;
   unsigned int *granule_ubound;
   unsigned int num_granules_alloc;
   unsigned int num_granules;
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
   double image_end_time;
}
Exprec_Info_Type;

#define PROCESS_METHOD_PRIVATE_DATA \
   Enum_Lookup_Type *enum_lookup; \
   Exprec_Cache_Method_Type *cache_method; \
   char *out_dirname; \
   char *out_basename; \
   char *archdir_path; \
   int ncid; \
   int processing_version; \
   double latest_radiance_timestamp_seen; \
   double outfile_timestamp_start; \
   double outfile_timestamp_end; \
   const char *product_type; \
   const char *exprec_type_string; \
   int exprec_type; \
   int granule_size; \
   unsigned int curr_mirror_step; \
   time_t when_last_erec_cached; \
   Granule_Schedule_Type sched;
#include "l0_format.h"

#define EXPREC_ENUM_TABLE_SIZE 16

static int Exprec_Cache_Method;

void set_exprec_cache_method (int method)
{
   Exprec_Cache_Method = method;
}

static void sched_free (Granule_Schedule_Type *sched)
{
   FREE(sched->granule_sizes);
   memset ((char *)sched, 0, sizeof(Granule_Schedule_Type));
}

static int sched_realloc (Granule_Schedule_Type *sched, unsigned int n)
{
   unsigned int *m = NULL;
   size_t len;

   if (sched == NULL)
     return -1;
   if (n < sched->num_granules_alloc)
     return 0;

   /* allocate more than we need */
   sched->num_granules_alloc = n + n/2;
   len = sched->num_granules_alloc * sizeof(unsigned int);

   if (NULL == (m = (unsigned int *)REALLOC (sched->granule_sizes, 2 * len)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }

   sched->granule_sizes = m;
   sched->granule_ubound = m + sched->num_granules_alloc;
   sched->num_granules = n;

   return 0;
}

static int compute_granule_sizes (unsigned int num_recs,
                                  unsigned int num_granules,
                                  Granule_Schedule_Type *sched)
{
   unsigned int i, n, remainder, ubound;
   size_t len = sched->num_granules_alloc * sizeof(unsigned int);

   if (0 != sched_realloc (sched, num_granules))
     return -1;

   memset ((char *)sched->granule_sizes, 0, len);
   memset ((char *)sched->granule_ubound, 0, len);

   n = num_recs / num_granules;
   if (n < 1) n = 1;
   remainder = num_recs - n * num_granules;

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
          {
             tell_vlog (TELL_MSGTYPE_ERROR, 0, "%s: computing granule sizes", __func__);
             return -1;
          }
        if (sched->granule_sizes[0] < sched->size_max)
          break;
        num_granules++;
     }

   return 0;
}

static double image_end_time (const IOCSDPC_Exprec_Type *erec)
{
   return (erec->image_start_time
           + erec->exposure_time
           + erec->frame_transfer_time);
}

static int define_outfile_vars (Process_Method_Type *pmt,
                                const TPInfo_Type *tpinfo,
                                IOCSDPC_Exprec_Type *erec,
                                const Radiance_Ident_Type *identp)
{
   const IOCSDPC_Image_Info_Item_Type *item;
   int ncid = pmt->ncid;
   int shuffle=LEVEL0_IMAGE_SHUFFLE_SELECT;
   int deflate=LEVEL0_IMAGE_DEFLATE_SELECT;
   int deflate_level=LEVEL0_IMAGE_DEFLATE_LEVEL;
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

   if (0 != write_attr_global_product_type (ncid, pmt->product_type))
     return -1;

   if (identp)
     {
        if (0 != tio_write_granule_ident_indices (ncid, identp->scan_num,
                                                  identp->granule_num))
          return -1;
        if (0 != tio_define_granule_flag_var (ncid))
          return -1;
        if (0 != tio_write_granule_flag_var (ncid, identp->granule_flag))
          return -1;
	/* (scan_type==0) means "standard radiance scan".
	 * (scan_type!=0) means the data will probably require custom processing.
	 * It is assumed that all exposures within a given scan_num will have the
	 * same value of scan_type, but this is not checked.
	 */
	if (0 != TIO_put_att (ncid, NC_GLOBAL, "scan_type", NC_INT, 1, &identp->scan_type))
	  return -1;
     }

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
        if (-1 == close_hidden (pmt->ncid, pmt->out_dirname, pmt->out_basename, pmt->archdir_path))
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

static void set_radiance_granule_flag (int curr_granule, int num_granules,
                                       int *pgranule_flag)
{
   int granule_flag = 0;
   if (curr_granule == 0)
     {
        granule_flag |= TEMPO_GRANULE_FLAG_IS_FIRST;
     }
   if (curr_granule == num_granules-1)
     {
        granule_flag |= TEMPO_GRANULE_FLAG_IS_LAST;
     }
   *pgranule_flag = granule_flag;
}

/* Create a new netCDF output file  */
static int new_outfile (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                        int granule_size, IOCSDPC_Exprec_Type *erec,
                        int curr_granule, int num_granules)
{
   char basename[MAX_BASENAME_SIZE];
   Radiance_Ident_Type radiance_ident = {0};
   Radiance_Ident_Type *identp = NULL;
   uint16_t sdpc_scan_label, scan_num, scan_type;
   int scan_num_int;

   pmt->outfile_timestamp_start = erec->image_start_time;
   pmt->outfile_timestamp_end = image_end_time (erec);
   pmt->exprec_type = erec->exprec_type;

   pmt->exprec_type_string = NULL;
   scan_num_int = -1;

   switch (pmt->exprec_type)
     {
      case IOCSDPC_EXPREC_TYPE_RAD:
	/* lowest 16 bits of erec->scan_label contain the 16 bit label
	 * for the scan provided in the SDPC-generated radiance scan plan. */
	sdpc_scan_label = erec->scan_label & 0x0000ffff;
	tio_parse_scan_label (sdpc_scan_label, &scan_type, &scan_num);
        if (scan_type == 0)
          {
             pmt->product_type = TEMPO_PROD_TYPE_RAD;
             pmt->exprec_type_string = TEMPO_PROD_TYPESTR_RAD;
          }
        else
          {
             pmt->product_type = TEMPO_PROD_TYPE_RAD_TWI;
             pmt->exprec_type_string = TEMPO_PROD_TYPESTR_RAD_TWI;
          }
        scan_num_int = scan_num;
	radiance_ident.scan_num = scan_num;
	radiance_ident.scan_type = scan_type;
        radiance_ident.granule_num = curr_granule + 1;
        set_radiance_granule_flag (curr_granule, num_granules,
                                   &radiance_ident.granule_flag);
        identp = &radiance_ident;
        break;
      case IOCSDPC_EXPREC_TYPE_DARK:
        pmt->product_type = TEMPO_PROD_TYPE_DRK;
        pmt->exprec_type_string = TEMPO_PROD_TYPESTR_DRK;
        break;
      case IOCSDPC_EXPREC_TYPE_IRR_WORK:
        pmt->product_type = TEMPO_PROD_TYPE_IRR;
        pmt->exprec_type_string = TEMPO_PROD_TYPESTR_IRR;
        break;
      case IOCSDPC_EXPREC_TYPE_IRR_REF:
        pmt->product_type = TEMPO_PROD_TYPE_IRR_REF;
        pmt->exprec_type_string = TEMPO_PROD_TYPESTR_IRR_REF;
        break;
      case IOCSDPC_EXPREC_TYPE_LIN_IRR:
        pmt->product_type = TEMPO_PROD_TYPE_IRR_LIN;
        pmt->exprec_type_string = TEMPO_PROD_TYPESTR_IRR_LIN;
        break;
      case IOCSDPC_EXPREC_TYPE_LIN_DARK:
        pmt->product_type = TEMPO_PROD_TYPE_DRK_LIN;
        pmt->exprec_type_string = TEMPO_PROD_TYPESTR_DRK_LIN;
        break;
      case IOCSDPC_EXPREC_TYPE_UNKNOWN:
        /* drop */
      default:
        pmt->product_type = "UNK";
        pmt->exprec_type_string = "UNK";
        break;
     }

   /* if ioclib_strdup failed, it already produced a log message */
   if (pmt->exprec_type_string == NULL)
     return -1;

   if (0 != verify_epoch (erec->common_header.epoch))
     return -1;

   FREE(pmt->archdir_path);
   pmt->archdir_path = NULL;
   if (0 != make_level0_archdir_path (&pmt->archdir_path, erec->image_start_time, scan_num_int,
                                      pmt->processing_version, pmt->product_type))
     return -1;

   if (-1 == make_level0_basename (basename, sizeof(basename), erec->image_start_time,
                                   pmt->processing_version, pmt->product_type, identp))
     return -1;

   tell_vinfo (0, "creating file %s/%s", pmt->out_dirname, basename);

   if ((-1 == create_hidden (pmt->out_dirname, basename, &pmt->ncid))
       || (-1 == write_std_global_metadata (pmt->ncid, &erec->common_header)))
     return -1;
   FREE(pmt->out_basename);
   if (NULL == (pmt->out_basename = ioclib_strdup (basename)))
     return -1;

   if (NULL == (pmt->enum_lookup = elt_open (EXPREC_ENUM_TABLE_SIZE)))
     return -1;

   pmt->granule_size = granule_size;

   return define_outfile_vars (pmt, tpinfo, erec, identp);
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

   return tio_sync (ncid);

return_status:
   ioclib_free (data);
   return -1;
}

static int process_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                          int process_all_erecs)
{
   Exprec_Cache_Method_Type *cmt = pmt->cache_method;
   Granule_Schedule_Type *sched = NULL;
   IOCSDPC_Exprec_Type *erec = NULL;
   char path[MAX_PATHLEN];
   unsigned int outfile_erec_count, outfile_cumulative_erec_count;
   size_t cache_index, num_erecs;
   int exprec_type, is_radiance, status = -1;

   /* If the cache is empty, do nothing.
    * Otherwise, the cache contains exposure records
    * of some type, and no output file is open.
    * Depending on the value of process_all_erecs, the goal is to pack
    * either some or all of the cached exposure records into netcdf granules.
    */

   if (0 != cmt->cache_open (cmt))
     return -1;

   (void) cmt->cache_num_recs (cmt, &num_erecs);

   /* If the cache is empty, do nothing */
   if (num_erecs == 0)
     {
        pmt->when_last_erec_cached = 0;
        return cmt->cache_close (cmt);
     }

   exprec_type = pmt->exprec_type;
   is_radiance = (exprec_type == IOCSDPC_EXPREC_TYPE_RAD);

   /* The schedule tells us how the exposure records should be packed
    * into granules. For radiances, the header tells us how many scan
    * steps were commanded, so a schedule based on that was defined
    * when the first frame of the scan arrived.
    * For other types of exposure records, we derive a packing schedule
    * based on the number of records we actually received. In most cases,
    * that will be a single granule, but we'll also handle the case when
    * multiple granules might be needed.
    */
   if (is_radiance == 0)
     {
        if (0 != schedule_granules (num_erecs, &pmt->sched))
          goto return_status;
     }

   sched = &pmt->sched;

   outfile_erec_count = 0;
   outfile_cumulative_erec_count = 0;

   for (cache_index = 0; cache_index < num_erecs; cache_index++)
     {
        unsigned int next_erec, curr_mirror_step, index;
        int need_new_outfile;

        if  (NULL == (erec = cmt->cache_erec_get (cmt)))
          goto return_status;

        path[0] = 0;
        (void) cmt->cache_erec_path (cmt, path, sizeof(path));

        /* Assert: exprec_type == erec->exprec_type */
        if (exprec_type != erec->exprec_type)
          {
             tell_verror (TELL_INTERNAL_ERROR,
                          "%s: unexpected exposure record type %d (expected %d) file=%s",
                          __func__, erec->exprec_type, exprec_type, path);
             goto return_status;
          }

        if (is_radiance)
          {
             if (NULL == iocsdpc_image_info_get_value (erec, "curr_mirror_step",
                                                       &curr_mirror_step))
               goto return_status;
             next_erec = curr_mirror_step;
          }
        else next_erec = outfile_cumulative_erec_count;

        need_new_outfile = ((pmt->ncid == INT_MAX)
                            || (outfile_erec_count == (unsigned int) pmt->granule_size)
                            || (next_erec >= sched->granule_ubound[sched->curr_granule]));

        if (need_new_outfile)
          {
             unsigned int k, granule_size;

             if (-1 == close_outfile (pmt))
               goto return_status;

             /* From granule schedule, find out which granule file
              * should hold the next exposure record.
              */
             for (k = sched->curr_granule; k < sched->num_granules; k++)
               {
                  if (next_erec < sched->granule_ubound[k])
                    break;
               }

             if (k >= sched->num_granules)
               {
                  tell_vinfo (0, "%s: bad erec: %s (k >= num_granules: k=%d, num_granules=%d)",
                              __func__, path, k, sched->num_granules);
                  if (0 != cmt->cache_erec_bad (cmt))
                    goto return_status;
                  continue;
               }

             sched->curr_granule = k;
             granule_size = sched->granule_sizes[k];

             /* We prefer to open an output file only when we can fill it. */
             if ((num_erecs - cache_index < granule_size)
                 && (process_all_erecs == 0))
               {
                  /* No need to open yet: wait for more data */
                  break;
               }

             if ((process_all_erecs != 0)
                 && (granule_size > num_erecs-cache_index)
                 && (is_radiance == 0))
               {
                  /* Must open now: size the granule to hold what's here */
                  granule_size = num_erecs - cache_index;
               }

             if (0 != new_outfile (pmt, tpinfo, granule_size, erec,
                                   sched->curr_granule,
                                   sched->num_granules))
               goto return_status;
             outfile_erec_count = 0;
          }

        if (is_radiance)
          {
             index = curr_mirror_step;
             if (sched->curr_granule > 0)
               {
                  index -= sched->granule_ubound[sched->curr_granule-1];
               }
          }
        else index = outfile_erec_count;

        tell_vinfo (2, "writing exprec: cache_index=%ld", cache_index);

        if (0 != write_exprec (pmt, tpinfo, index, erec))
          goto return_status;
        outfile_erec_count++;
        outfile_cumulative_erec_count++;

        iocsdpc_exprec_close (erec);
        erec = NULL;

        (void) cmt->cache_erec_done (cmt);
     }

   status = 0;
return_status:
   if (status)
     {
        (void) cmt->cache_num_recs (cmt, &num_erecs);
        tell_vlog (TELL_MSGTYPE_ERROR, 0, "processing exprec cache: num_erecs=%ld",
                   num_erecs);
     }

   (void) cmt->cache_close (cmt);

   if (erec) iocsdpc_exprec_close (erec);
   (void) close_outfile (pmt);
   return status;
}

static int radiance_belongs_to_curr_granule (const Process_Method_Type *pmt,
                                             const Exprec_Info_Type *exprec_info)
{
   const Granule_Schedule_Type *sched = &pmt->sched;
   unsigned int step_ubound = sched->granule_ubound[sched->curr_granule];
   if (exprec_info->exprec_type != IOCSDPC_EXPREC_TYPE_RAD)
     return 0;
   return (exprec_info->curr_mirror_step < step_ubound);
}

static int classify_erec (IOCSDPC_Exprec_Type *erec, Exprec_Info_Type *info)
{
   info->exprec_type = erec->exprec_type;
   info->image_end_time = image_end_time (erec);

   if ((NULL == iocsdpc_image_info_get_value (erec, "curr_mirror_step", &info->curr_mirror_step))
       || (NULL == iocsdpc_image_info_get_value (erec, "num_mirror_steps", &info->num_mirror_steps)))
     return -1;

   if (info->exprec_type == IOCSDPC_EXPREC_TYPE_RAD)
     {
        if ((info->num_mirror_steps == 0)
            || (info->curr_mirror_step > info->num_mirror_steps))
          {
             tell_vlog (TELL_MSGTYPE_INFO, 0,
                        "%s: invalid radiance exposure record header: %s=%d %s=%d", __func__,
                        "curr_mirror_step", info->curr_mirror_step,
                        "num_mirror_steps", info->num_mirror_steps);
             return -1;
          }
     }

   return 0;
}

static int flush_cache (Process_Method_Type *pmt, const TPInfo_Type *tpinfo)
{
   return process_cache (pmt, tpinfo, 1);
}

static int process_exprec1 (Process_Method_Type *pmt,
                            const TPInfo_Type *tpinfo,
                            const Exprec_Info_Type *exprec_info,
                            const char *file, size_t exprec_index)
{
   Exprec_Cache_Method_Type *cmt = pmt->cache_method;
   int is_radiance, is_new_type, is_radiance_new_scan;

   /* Invariant:  the cache always contains only records
    * that belong together and could share the same granule.
    */

   is_new_type = (exprec_info->exprec_type != pmt->exprec_type);
   is_radiance = (exprec_info->exprec_type == IOCSDPC_EXPREC_TYPE_RAD);
   is_radiance_new_scan = (is_radiance && (exprec_info->curr_mirror_step
                                           < pmt->curr_mirror_step));

   /* Reject duplicate radiance */
   if (is_radiance && (exprec_info->curr_mirror_step == pmt->curr_mirror_step))
     return -1;

   if (is_radiance)
     {
        pmt->latest_radiance_timestamp_seen = exprec_info->image_end_time;
     }

   /* The cache contains records of type pmt->exprec_type.
    * Process the cache before changing the value of pmt->exprec_type. */
   if (is_new_type || is_radiance_new_scan)
     {
        if (0 != process_cache (pmt, tpinfo, 1))
          return -1;
     }

   pmt->exprec_type = exprec_info->exprec_type;
   if (0 != cmt->cache_erec (cmt, file, exprec_index))
     return -1;
   pmt->when_last_erec_cached = time(NULL);

   /* Non-radiance exposure records are cached until something
    * triggers cache processing: 1) arrival of a different
    * exposure type, 2) exceeding a cache size threshold,
    * or 3) expiration of a timer. */

   if (is_radiance == 0)
     {
        size_t num_erecs;
        (void) cmt->cache_num_recs (cmt, &num_erecs);
        if (num_erecs < 1.5*pmt->sched.size_default)
          {
             return 0;
          }
        return process_cache (pmt, tpinfo, 0);
     }

   /* For radiance exposure records, we know how many to expect,
    * so we generate granules according to a pre-planned schedule.
    */

   pmt->curr_mirror_step = exprec_info->curr_mirror_step;

   if (is_radiance_new_scan)
     {
        /* New schedule upon new scan resets sched->curr_granule to 0 */
        if (0 != schedule_granules (exprec_info->num_mirror_steps + 1,
                                    &pmt->sched))
          return -1;
     }

   if (radiance_belongs_to_curr_granule (pmt, exprec_info))
     return 0;

   /* When we see a frame that doesn't belong to the current
    * radiance granule, that granule is complete. */
   return process_cache (pmt, tpinfo, 0);
}

static int process_exprec (Process_Method_Type *pmt, const TPInfo_Type *tpinfo,
                           const char *file, void *client_data)
{
   IOCSDPC_Common_Header_Type chdr = {0};
   IOCSDPC_Exprec_Type *erec = NULL;
   size_t exprec_index;
   int fd, status;

   if (-1 == (fd = iocsdpc_open_file_read (file, 0, &chdr)))
     return -1;
   if (NULL == (erec = iocsdpc_exprec_fdopen_read (file, fd, &chdr)))
     {
        ioclib_fd_close (fd);
        return -1;
     }

   for (exprec_index=0; /* until EOF */ ;exprec_index++)
     {
        Exprec_Info_Type exprec_info = {0};

        if (0 != classify_erec (erec, &exprec_info))
          {
             tell_vlog (TELL_MSGTYPE_ERROR, 0, "%s: classifying exposure record: %s",
                        __func__, file);
             goto return_status;
          }

        if (0 != process_exprec1 (pmt, tpinfo, &exprec_info, file, exprec_index))
          goto return_status;

        if (pmt->pmt_post_process_callback != NULL)
          {
             if (0 != pmt->pmt_post_process_callback (pmt, client_data))
               goto return_status;
          }

        /* status=1  => have next record
         * status=0  => EOF
         * status=-1 => error occurred
         */
        if (1 != (status = iocsdpc_exprec_open_next (erec)))
          break;
     }

return_status:
   iocsdpc_exprec_close (erec);

   return status;
}

static int parse_exprec_params (config_t *cfg, Process_Method_Type *pmt)
{
   config_setting_t *s;
   const char *out_dirname;
   int size_default, size_max;

   if (NULL == (s = config_lookup (cfg, "exprec")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing 'exprec' in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "processing_version", &pmt->processing_version))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "output_dir", &out_dirname))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "granule_size_default", &size_default))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "granule_size_max", &size_max)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading 'exprec' parameters in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (pmt->out_dirname = expand_string (out_dirname)))
     return -1;

   pmt->sched.size_default = size_default;
   pmt->sched.size_max = size_max;

   return 0;
}

static void delete_exprec (Process_Method_Type *pmt)
{
   if (pmt == NULL)
     return;
   (void) close_outfile (pmt);
   sched_free (&pmt->sched);
   FREE(pmt->out_basename);
   FREE(pmt->out_dirname);
   FREE(pmt->archdir_path);

   if (pmt->cache_method)
     {
        Exprec_Cache_Method_Type *cmt = pmt->cache_method;
        cmt->cache_delete (cmt);
     }

   FREE(pmt);
}

static int query_latest_timestamp (Process_Method_Type *pmt, int exprec_type, double *timestamp)
{
   switch (exprec_type)
     {
      case IOCSDPC_EXPREC_TYPE_RAD:
        *timestamp = pmt->latest_radiance_timestamp_seen;
        break;

      default:
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR, "%s: exprec_type = %d", __func__, exprec_type);
        return -1;
        break;
     }

   return 0;
}

static int query_when_last_erec_cached (const Process_Method_Type *pmt, time_t *when)
{
   *when = pmt->when_last_erec_cached;
   return 0;
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

   pmt->pmt_process = process_exprec;
   pmt->pmt_delete = delete_exprec;
   pmt->pmt_flush_cache = flush_cache;
   pmt->pmt_query_latest_timestamp = query_latest_timestamp;
   pmt->pmt_query_when_last_rec_cached = query_when_last_erec_cached;

   pmt->latest_radiance_timestamp_seen = -1.0;

   pmt->out_basename = NULL;
   pmt->ncid = INT_MAX;
   pmt->outfile_timestamp_start = -1.0;
   pmt->outfile_timestamp_end = -1.0;
   pmt->enum_lookup = NULL;

   pmt->exprec_type = INT_MAX;
   pmt->exprec_type_string = NULL;
   pmt->curr_mirror_step = UINT_MAX;

   switch (Exprec_Cache_Method)
     {
      case EXPREC_CACHE_DISK:
        pmt->cache_method = open_erec_cache_disk (cfg);
        break;

      case EXPREC_CACHE_MEM:
        pmt->cache_method = open_erec_cache_mem (cfg);
        break;

      default:
        pmt->cache_method = NULL;
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: unsupported exprec cache method: %d",
                     __func__, Exprec_Cache_Method);
        break;
     }

   if (NULL == pmt->cache_method)
     {
        delete_exprec (pmt);
        pmt = NULL;
     }

   return pmt;
}
