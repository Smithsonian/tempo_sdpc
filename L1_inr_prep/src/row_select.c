/** @file row_select.c
 *  @brief Manage SMC, IRU time series subsetting
 */
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ioclib.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "row_select.h"
#include "bsearch.h"

enum
{
   FILE_ERROR_OCCURRED    = -1,
   FILE_FOLLOWS_INTERVAL  =  1,
   FILE_OVERLAPS_INTERVAL =  2,
   FILE_PRECEDES_INTERVAL =  3
};

static void row_select_free1 (Row_Select_Type *rst)
{
   if (rst == NULL)
     return;
   ioclib_free(rst->file);
   FREE(rst->times);
   FREE(rst);
}

void row_select_free (Row_Select_Type *rst)
{
   if (rst == NULL)
     return;

   while (rst)
     {
        Row_Select_Type *next = rst->next;
        row_select_free1(rst);
        rst = next;
     }
}

static Row_Select_Type *alloc_row_select (const char *file)
{
   Row_Select_Type *rst = NULL;

   if (NULL == (rst = (Row_Select_Type *)MALLOC (sizeof *rst)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)rst, 0, sizeof *rst);

   if (NULL == (rst->file = ioclib_strdup (file)))
     return NULL;

   return rst;
}

static int read_times (Row_Select_Type *rst, const char *time_var, int ncid)
{
   TIO_Var_Info_Type info;
   size_t size_times;

   if (time_var == NULL)
     {
        time_var = "time";
     }

   if (0 != TIO_inq_var (ncid, time_var, &info))
     return -1;

   rst->num_times = info.dimlens[0];

   size_times = rst->num_times * sizeof(double);
   if (NULL == (rst->times = (double *) MALLOC (size_times)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   if (0 != TIO_get_var_section (ncid, time_var, &rst->start, &rst->num_times,
                                 NC_DOUBLE, rst->times))
     {
        return -1;
     }

   return 0;
}

static int apply_selection (Row_Select_Type *rst,
                            double time_beg, double time_end)
{
   double frac;

   if (time_beg < rst->times[0])
     rst->start = 0;
   else if (time_beg > rst->times[rst->num_times-1])
     rst->start = rst->num_times;
   else
     {
        int b = bsearch_d (time_beg, rst->times, rst->num_times);
        if (b < 0) return -1;
        rst->start = b;
     }

   if (time_end < rst->times[0])
     rst->count = 0;
   else if (time_end > rst->times[rst->num_times-1])
     rst->count = rst->num_times - rst->start;
   else
     {
        int e = bsearch_d (time_end, rst->times, rst->num_times);
        if (e < 0) return -1;
        if ((e + 1) < rst->num_times) e += 1;
        rst->count = e - rst->start + 1;
     }

   if ((rst->start < rst->num_times)
       && (rst->count > 0))
     {
        double t0 = rst->times[rst->start];
        double t1 = rst->times[rst->start + rst->count - 1];
        if (time_end != time_beg)
          frac = (t1 - t0) / (time_end - time_beg);
        else
          frac = 1.0;
     }
   else frac = 0.0;

   tell_vlog (TELL_MSGTYPE_INFO, 1,
              "file has %d samples spanning %g of time interval",
              rst->count, frac);

   return 0;
}

static int estimate_padding (int ncid, int num_pad,
                             double *time_beg, double *time_end)
{
   TIO_Var_Info_Type time_var_info;
   const char *sample_hz_attname = "average_sample_frequency_hz";
   const char *time_varname = "time";
   double average_sample_frequency_hz, dt;
   size_t len_size_t;
   int xtype;

   if (0 != TIO_inq_var (ncid, time_varname, &time_var_info))
     return -1;

   if (NC_NOERR == nc_inq_att (ncid, time_var_info.varid, sample_hz_attname, &xtype, &len_size_t))
     {
        if (0 != TIO_get_att (ncid, time_var_info.varid, sample_hz_attname,
                              NC_DOUBLE, &average_sample_frequency_hz))
          {
             return -1;
          }
        dt = num_pad / average_sample_frequency_hz;
     }
   else
     {
#define MAX_NUM_TIME_SAMPLES 32

        double t_array[MAX_NUM_TIME_SAMPLES], delta_t;
        int size = time_var_info.dimlens[0];
        int i, num, start, count;

        start = 0;
        count = (size < MAX_NUM_TIME_SAMPLES) ? size : MAX_NUM_TIME_SAMPLES;

        if (0 != TIO_get_var_section (ncid, time_varname, &start, &count, NC_DOUBLE, t_array))
          return -1;

        delta_t = 0.0;
        num = 0;
        for (i = 1; i < count; i++)
          {
             double t0 = t_array[i-1];
             double t1 = t_array[i];
             if ((isfinite(t0) && (t0 != TIO_FILL_DOUBLE))
                 && (isfinite(t1) && (t1 != TIO_FILL_DOUBLE)))
               {
                  delta_t += t1-t0;
                  num++;
               }
          }
        if (num > 0) delta_t /= num;
        dt = num_pad * delta_t;
     }

   *time_beg -= dt;

   if (0) /* No need to pad beyond the endpoint */
     {
        *time_end += dt;
     }

   return 0;
}

static int examine_file (Row_Select_Type **rstp, int ncid, const char *file,
                         double time_beg, double time_end)
{
   double time_coverage_start, time_coverage_end;
   int status = FILE_ERROR_OCCURRED;

   *rstp = NULL;

   if ((0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch",
                          NC_DOUBLE, &time_coverage_start))
       ||(0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch",
                            NC_DOUBLE, &time_coverage_end)))
     {
        return status;
     }

   if (time_end < time_coverage_start)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "file follows interval: time_end:%f < time_coverage_start:%f",
                   time_end, time_coverage_start);
        status = FILE_FOLLOWS_INTERVAL;
     }
   else if (time_coverage_end < time_beg)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "file precedes interval: time_coverage_end:%f < time_beg:%f",
                   time_coverage_end, time_beg);
        status = FILE_PRECEDES_INTERVAL;
     }
   else
     {
        Row_Select_Type *rst = NULL;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "file overlaps interval");

        if (NULL == (rst = alloc_row_select (file)))
          return status;

        *rstp = rst;
        status = FILE_OVERLAPS_INTERVAL;
     }

   return status;
}

int row_select_scan (double time_beg, double time_end, int num_pad,
                     int num_files, char **file_list,
                     const char *group_path, const char *time_var,
                     Row_Select_Type **rstp)
{
   Row_Select_Type *rst_head = NULL;
   int return_status = -1;
   int i, ncid = 0;

   for (i = 0; i < num_files; i++)
     {
        Row_Select_Type *rst;
        const char *file = file_list[i];
        double time_beg_pad = time_beg;
        double time_end_pad = time_end;
        int grp, status;

        if (0 == strcasecmp (file, "NONE"))
          continue;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "Examining %s", file);

        ncid = 0;
        if (0 != TIO_open (file, NC_NOWRITE, &ncid))
          goto cleanup_and_return;
        if (group_path)
          {
             if (NC_NOERR != nc_inq_grp_full_ncid (ncid, group_path, &grp))
               {
                  (void) TIO_close(ncid);
                  ncid = 0;
                  continue;
               }
          }
        else grp = ncid;

        if (num_pad > 0)
          {
             /* When padding, expand the time interval to make sure
              * we include all relevant files */
             if (0 != estimate_padding (grp, num_pad, &time_beg_pad, &time_end_pad))
               goto cleanup_and_return;
          }

        status = examine_file (&rst, ncid, file, time_beg_pad, time_end_pad);
        if (status == FILE_ERROR_OCCURRED)
          {
             goto cleanup_and_return;
          }
        else if (status == FILE_PRECEDES_INTERVAL)
          {
             tell_vlog (TELL_MSGTYPE_INFO, 2, "file precedes interval: %s", file);
             (void) TIO_close (ncid);
             ncid = 0;
             continue;
          }
        else if (status == FILE_FOLLOWS_INTERVAL)
          {
             tell_vlog (TELL_MSGTYPE_INFO, 2, "file follows interval: %s", file);
             (void) TIO_close (ncid);
             ncid = 0;
             break;
          }
        else if (status == FILE_OVERLAPS_INTERVAL)
          {
             tell_vlog (TELL_MSGTYPE_INFO, 1, "file overlaps interval: %s", file);
             if (0 != read_times (rst, time_var, grp))
               goto cleanup_and_return;
             (void) TIO_close (ncid);
             ncid = 0;

             tell_vlog (TELL_MSGTYPE_INFO, 1, "read times: %s", file);
             if (0 != apply_selection (rst, time_beg_pad, time_end_pad))
               goto cleanup_and_return;

             if (rst_head == NULL)
               rst_head = rst;
             else
               {
                  Row_Select_Type *p;
                  for (p = rst_head; p != NULL; p = p->next)
                    {
                       if (p->next == NULL)
                         {
                            p->next = rst;
                            break;
                         }
                    }
               }
             tell_vlog (TELL_MSGTYPE_INFO, 1, "appended time samples from: %s", file);
          }
     }

   if (rst_head == NULL)
     {
        return_status = 0;
        tell_vwarn (0, "%s: no samples in time interval [%0.4f, %0.4f) with num_pad=%d",
                    __func__, time_beg, time_end, num_pad);
        goto cleanup_and_return;
     }

   return_status = 0;
cleanup_and_return:
   if (ncid != 0)
     {
        (void) TIO_close (ncid);
     }

   if (return_status)
     {
        row_select_free (rst_head);
        rst_head = NULL;
     }

   *rstp = rst_head;

   return return_status;
}
