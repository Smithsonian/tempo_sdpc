/** @file
 *  @brief TEMPO-specific utility functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define __USE_XOPEN   /* for strptime */
#include <time.h>

#include <stdarg.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>

#include <cfortran.h>
#include <netcdf.h>
#include <tell.h>

#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

#define DELIM_TIMESTAMP_FORMAT   "%Y-%m-%dT%H:%M:%SZ"
#define NODELIM_TIMESTAMP_FORMAT "%Y%m%dT%H%M%SZ"

struct TIO_Scan_Ident_Type
{
   _pTIO_Granule_Ident_Type *granule_ident;
};

int _pTIO_define_processing_level (int grp, int level)
{
   int status, enum_typeid;
   static TIO_Enum_Type enum_table[] =
     {
        {"level-0",  TIO_PROC_LEVEL_0},
        {"level-1a", TIO_PROC_LEVEL_1A},
        {"level-1b", TIO_PROC_LEVEL_1B},
        {"level-2",  TIO_PROC_LEVEL_2},
        {"level-3",  TIO_PROC_LEVEL_3},
        TIO_ENUM_TABLE_END
     };

   if (-1 == TIO_define_enum_table (grp, "processing_level_enum", NC_INT, enum_table, &enum_typeid))
     return -1;
   status = nc_put_att (grp, NC_GLOBAL,
                        "processing_level", enum_typeid, 1, &level);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining processing_level attribute (%s)",
                     __func__, nc_strerror(status));
        return -1;
     }

   return 0;
}

static int read_granule_ident_indices (int ncid, _pTIO_Granule_Ident_Type *gid)
{
   int start, count, attid;

   /* When granule ident indices aren't present, we assume this granule doesn't need them.
    * In this case, we set the indices to zero to indicate that they aren't used.
    */
   if (NC_ENOTATT == nc_inq_attid (ncid, NC_GLOBAL, "granule_num", &attid))
     {
        gid->scan_num = 0;
        gid->granule_num = 0;
        gid->granule_flag = 0;
        return 0;
     }

   if ((-1 == TIO_get_att (ncid, NC_GLOBAL, "scan_num", NC_INT, &gid->scan_num))
       ||(-1 == TIO_get_att (ncid, NC_GLOBAL, "granule_num", NC_INT, &gid->granule_num)))
     {
        return -1;
     }

   start = 0;
   count = 1;
   if (0 != TIO_get_var_section (ncid, TEMPO_VAR_GRANULE_FLAG, &start, &count, NC_INT, &gid->granule_flag))
     return -1;

   return 0;
}

int _pTIO_read_granule_ident (int ncid, _pTIO_Granule_Ident_Type *gid)
{
   if (0 != read_granule_ident_indices (ncid, gid))
     return -1;

   memset (gid->tstart_str, 0, MAX_ISOTIME_LEN);
   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, gid->tstart_str))
     return -1;

   memset (gid->tend_str, 0, MAX_ISOTIME_LEN);
   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, gid->tend_str))
     return -1;

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &gid->tstart))
     return -1;
   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, &gid->tend))
     return -1;

   return 0;
}

int tio_define_granule_flag_var (int ncid)
{
   static _pText_Attr_Type granule_flag_attrs[] =
     {
        {"flag_meanings",
             "is_first_granule_of_scan, is_last_granule_of_scan, is_telemetry_only"},
        _pTEXT_ATTRS_END
     };
   int flag_masks[] = {
      TEMPO_GRANULE_FLAG_IS_FIRST,
      TEMPO_GRANULE_FLAG_IS_LAST,
      TEMPO_GRANULE_FLAG_IS_TELEMETRY_ONLY
   };
   int num_masks = sizeof(flag_masks)/sizeof(*flag_masks);
   int varid, status;

   if (-1 == _pTIO_define_var_with_text_attrs (ncid, TEMPO_VAR_GRANULE_FLAG, NC_INT, 0, NULL, granule_flag_attrs, &varid))
     return -1;
   if (NC_NOERR != (status = nc_put_att_int (ncid, varid, "flag_masks", NC_INT,
                                             num_masks, flag_masks)))
     {
        tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining int attribute %s (%s)",
                     __func__, "flag_masks", nc_strerror(status));
        return -1;
     }

   return 0;
}

static int write_granule_ident_indices (int ncid, int scan_num,
                                        int granule_num, int granule_flag)
{
   int status, start, count, varid;

   /* In use, granule_num and scan_num are both positive values.
    * A value of zero indicates that the index will not be used.
    */
   if ((granule_num == 0) || (scan_num == 0))
     return 0;

   status = ((-1 == TIO_put_att (ncid, NC_GLOBAL, "scan_num", NC_INT, 1, &scan_num))
             ||(-1 == TIO_put_att (ncid, NC_GLOBAL, "granule_num", NC_INT, 1, &granule_num)));
   if (status)
     return -1;

   status = nc_inq_varid (ncid, TEMPO_VAR_GRANULE_FLAG, &varid);
   if (NC_NOERR != status)
     {
        if (status != NC_ENOTVAR)
          return -1;
        if (0 != tio_define_granule_flag_var (ncid))
          return -1;
        /* drop */
     }

   start = 0;
   count = 1;
   status = TIO_put_var_section (ncid, TEMPO_VAR_GRANULE_FLAG, &start, &count,
                                 NC_INT, &granule_flag);

   return status ? -1 : 0;
}

static int write_granule_ident_times (int ncid, const _pTIO_Granule_Ident_Type *gid)
{
   size_t len;

   len = strlen (gid->tstart_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, len, gid->tstart_str))
     return -1;

   len = strlen (gid->tend_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, len, gid->tend_str))
     return -1;

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, 1, &gid->tstart))
     return -1;

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, 1, &gid->tend))
     return -1;

   return 0;
}

int _pTIO_write_granule_ident (int ncid, const _pTIO_Granule_Ident_Type *gid)
{
   if (0 != write_granule_ident_indices (ncid, gid->scan_num, gid->granule_num,
                                         gid->granule_flag))
     return -1;

   return write_granule_ident_times (ncid, gid);
}

int tio_write_granule_ident_indices (int ncid, int scan_num, int granule_num,
                                    int granule_flag)
{
   return write_granule_ident_indices (ncid, scan_num, granule_num,
                                       granule_flag);
}

int tio_write_granule_ident_times (int ncid, double tstart, double tend)
{
   _pTIO_Granule_Ident_Type gid = {0};

   if (0 != tio_write_epoch_timestamp (ncid, NC_GLOBAL))
     return -1;

   gid.tstart = tstart;
   gid.tend = tend;

   if (0 != TIO_mktimestamp_str (tstart, 0, gid.tstart_str, sizeof(gid.tstart_str)))
     return -1;
   if (0 != TIO_mktimestamp_str (tend, 0, gid.tend_str, sizeof(gid.tend_str)))
     return -1;

   return write_granule_ident_times (ncid, &gid);
}

int _pTIO_parse_timestr (const char *timestr, struct tm *ptm)
{
   memset ((char *)ptm, 0, sizeof (struct tm));
   if (NULL == strptime (timestr, DELIM_TIMESTAMP_FORMAT, ptm))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: strptime failed: %s",
                     __func__, timestr);
        return -1;
     }

   return 0;
}

int __tio_filename_string_indexed (char *buf, int bufsize,
                                   double tstart, const char *label, int level, int version,
                                   int scan_num, int granule_num)
{
   char timestr[MAX_ISOTIME_LEN];

   if (0 != TIO_mktimestamp_str (tstart, 0, timestr, sizeof(timestr)))
     return -1;

   /* TEMPO_<label>_Ld_Vdd_<time>_SdddGdd.nc */
   return snprintf (buf, bufsize,
                    "TEMPO_%s_L%d_V%02d_%s_S%03dG%02d.nc",
                    label, level, version, timestr,
                    scan_num, granule_num);
}

int __tio_filename_string (char *buf, int bufsize,
                           double tstart, const char *label, int level, int version)
{
   char timestr[MAX_ISOTIME_LEN];

   if (0 != TIO_mktimestamp_str (tstart, 0, timestr, sizeof(timestr)))
     return -1;

   /* TEMPO_<label>_Ld_Vdd_<time>.nc */
   return snprintf (buf, bufsize,
                    "TEMPO_%s_L%d_V%02d_%s.nc",
                    label, level, version, timestr);
}

static int
_pTIO_filename_from_granule_ident (const _pTIO_Granule_Ident_Type *gid,
                                   const char *label, int level, int version,
                                   char *buf, int bufsize)
{
   int meaningful_granule_indices;
   int n;

   meaningful_granule_indices = ((gid->granule_num > 0)
                                 && (gid->scan_num > 0));

   if (meaningful_granule_indices)
     {
        n = __tio_filename_string_indexed (buf, bufsize, gid->tstart, label, level, version,
                                           gid->scan_num, gid->granule_num);
     }
   else
     {
        n = __tio_filename_string (buf, bufsize, gid->tstart, label, level, version);
     }

   if (n >= bufsize)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: filename length %d truncated to buffer size %d)",
                     __func__, n, bufsize);
        return -1;
     }

   return 0;
}

int TIO_copy_granule_ident (int ncid_from, int ncid_to)
{
   _pTIO_Granule_Ident_Type gid = {0};

   if ((-1 == _pTIO_read_granule_ident (ncid_from, &gid))
       || (-1 == _pTIO_write_granule_ident (ncid_to, &gid)))
     return -1;

   return 0;
}

static int same_granule_ident (const _pTIO_Granule_Ident_Type *gid1,
                               const _pTIO_Granule_Ident_Type *gid2)
{
   return ((gid1->scan_num == gid2->scan_num)
           && (gid1->granule_num == gid2->granule_num)
           && (gid1->granule_flag == gid2->granule_flag)
           && (gid1->tstart == gid2->tstart)
           && (gid1->tend == gid2->tend)
           && (0 == strncmp (gid1->tstart_str, gid2->tstart_str, MAX_ISOTIME_LEN))
           && (0 == strncmp (gid1->tend_str, gid2->tend_str, MAX_ISOTIME_LEN)));
}

int TIO_same_granule_ident (int ncid1, int ncid2)
{
   _pTIO_Granule_Ident_Type gid1 = {0};
   _pTIO_Granule_Ident_Type gid2 = {0};

   if ((-1 == _pTIO_read_granule_ident (ncid1, &gid1))
       || (-1 == _pTIO_read_granule_ident (ncid2, &gid2)))
     return -1;

   return same_granule_ident (&gid1, &gid2);
}

int TIO_filename_from_granule (int ncid, const char *label, int level, int version,
                               char *buf, int bufsize)
{
   _pTIO_Granule_Ident_Type gid = {0};

   if (-1 == _pTIO_read_granule_ident (ncid, &gid))
     return -1;

   return _pTIO_filename_from_granule_ident (&gid, label, level, version, buf, bufsize);
}

int TIO_label_product (int ncid, const char *product_type, int version)
{
   if (product_type != NULL)
     {
        size_t len = strlen (product_type) + 1;
        if (-1 == TIO_put_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, len, product_type))
          return -1;
     }

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "processing_version", NC_INT, 1, &version))
     return -1;

   return 0;
}

static void free_granule_ident (_pTIO_Granule_Ident_Type *gid)
{
   if (gid == NULL)
     return;
   TIO_FREE(gid);
}

static _pTIO_Granule_Ident_Type *new_granule_ident (void)
{
   _pTIO_Granule_Ident_Type *gid = NULL;
   if (NULL == (gid = (_pTIO_Granule_Ident_Type *) TIO_MALLOC (sizeof *gid)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)gid, 0, sizeof *gid);
   return gid;
}

static void free_scan_ident (TIO_Scan_Ident_Type *lst)
{
   _pTIO_Granule_Ident_Type *gid;
   if (NULL == lst)
     return;
   gid = lst->granule_ident;
   while (gid != NULL)
     {
        _pTIO_Granule_Ident_Type *next = gid->next;
        free_granule_ident (gid);
        gid = next;
     }
   TIO_FREE(lst);
}

TIO_Scan_Ident_Type *TIO_new_scan_ident (void)
{
   TIO_Scan_Ident_Type *lst = NULL;

   if (NULL == (lst = (TIO_Scan_Ident_Type *)TIO_MALLOC (sizeof *lst)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   lst->granule_ident = NULL;

   return lst;
}

void TIO_free_scan_ident (TIO_Scan_Ident_Type *lst)
{
   free_scan_ident (lst);
}

int TIO_attach_granule_ident (int ncid, TIO_Scan_Ident_Type *lst)
{
   _pTIO_Granule_Ident_Type *item = NULL;
   int status;

   if (lst == NULL)
     return -1;

   if (NULL == (item = new_granule_ident ()))
     return -1;
   if (-1 == _pTIO_read_granule_ident (ncid, item))
     {
        free_granule_ident (item);
        return -1;
     }

   if (lst->granule_ident == NULL)
     {
        lst->granule_ident = item;
        return 0;
     }

   if (item->scan_num == lst->granule_ident->scan_num)
     {
        status = 0;
     }
   else
     {
        Tell_verror (TELL_APPLICATION_ERROR,
                     "%s: scan_num mismatch: new item has scan_num=%d  granule list has scan_num=%d",
                     __func__,
                     item->scan_num,
                     lst->granule_ident->scan_num);
        status = -1;
     }

   free_granule_ident (item);

   return status;
}

enum
{
   NODELIM_TIMESTAMP = 0,
     DELIM_TIMESTAMP = 1
};

int _pTIO_tempo_time_from_utc_timestr (const char *str, double *tai_sec)
{
   struct tm tm;
   double utc;

   /* Note that parsing a time stamp string is likely to yield
    * only a low precision time stamp, because such a string
    * most likely has only integer seconds.  For a high precision
    * time stamp, it's better to read a floating point value,
    * if one is available */

   if (0 != _pTIO_parse_timestr (str, &tm))
     return -1;

   if ((utc = timegm (&tm)) < 0)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: timegm failed: s=%s",
                     __func__, str);
        return -1;
     }

   return tio_time_utc_to_tempo (utc, tai_sec);
}

int TIO_mktimestamp_str (double tempo_tai_offset,
                         int delim, char *buf, int bufsize)
{
   struct tm tm;
   double utc;
   time_t tt;
   int status;

   if (0 != tio_time_tempo_to_utc (tempo_tai_offset, &utc))
     return -1;
   tt = (time_t)utc;

   memset ((char *)&tm, 0, sizeof(struct tm));
   if (NULL == gmtime_r (&tt, &tm))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: gmtime_r failed: tt=%ld",
                     __func__, tt);
        return -1;
     }

   if (delim == 0)
     status = strftime (buf, bufsize, NODELIM_TIMESTAMP_FORMAT, &tm);
   else
     status = strftime (buf, bufsize, DELIM_TIMESTAMP_FORMAT, &tm);

   if (0 == status)
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: strftime failed, tt=%ld",
                     __func__, tt);
        return -1;
     }

   return 0;
}

int tio_write_epoch_timestamp (int ncid, int varid)
{
   char buf[MAX_ISOTIME_LEN];
   int len;

   if (0 != TIO_mktimestamp_str (0.0, 1, buf, sizeof(buf)))
     return -1;

   len = strlen(buf);

   if (-1 == TIO_put_att (ncid, varid, "time_reference", NC_CHAR, len, buf))
     return -1;

   return 0;
}

int TIO_write_timestamp (int ncid, int varid, const char *attr_name,
                         double secs_since_tempo_epoch)
{
   char epoch_str[MAX_ISOTIME_LEN];
   char timestamp_str[MAX_ISOTIME_LEN];
   char attr_name_since_epoch[TIO_MAX_NAME_LEN];
   double file_epoch_tai, secs_since_file_epoch, tempo_epoch_utc;
   time_t epoch_tt, tempo_epoch_tt;
   struct tm tm;
   size_t len;
   int n;

   /* Generate the UTC string for the timestamp */
   if (0 != TIO_mktimestamp_str (secs_since_tempo_epoch, DELIM_TIMESTAMP, timestamp_str, MAX_ISOTIME_LEN))
     return -1;

   /* The floating point TAI timestamp should be consistent
    * with the file's defined reference epoch, which may not
    * be the TEMPO epoch. */
   memset ((char *)epoch_str, 0, MAX_ISOTIME_LEN);
   if (0 != TIO_get_att (ncid, NC_GLOBAL, "time_reference", NC_CHAR, epoch_str))
     return -1;
   if (-1 == _pTIO_parse_timestr (epoch_str, &tm))
     return -1;
   epoch_tt = timegm (&tm);

   tempo_epoch_utc = tio_tempo_epoch_timet ();
   tempo_epoch_tt = (time_t)tempo_epoch_utc;

   if (epoch_tt == tempo_epoch_tt)
     {
        secs_since_file_epoch = secs_since_tempo_epoch;
     }
   else
     {
        double tempo, file_epoch_utc, tempo_epoch_tai, tai_timestamp;

        if (0 != tio_time_tempo_to_tai (0.0, &tempo_epoch_tai))
          return -1;
        tai_timestamp = tempo_epoch_tai + secs_since_tempo_epoch;

        file_epoch_utc = (double) epoch_tt;
        if ((0 != tio_time_utc_to_tempo (file_epoch_utc, &tempo))
            || (0 != tio_time_tempo_to_tai (tempo, &file_epoch_tai)))
          return -1;

        secs_since_file_epoch = tai_timestamp - file_epoch_tai;
     }

   n = snprintf (attr_name_since_epoch, TIO_MAX_NAME_LEN, "%s_since_epoch", attr_name);
   if (n >= TIO_MAX_NAME_LEN)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: snprintf failed (output variable name exceeds buffer size=%d characters)",
                     __func__, TIO_MAX_NAME_LEN);
        return -1;
     }

   len = strlen(timestamp_str) + 1;
   if ((0 != TIO_put_att (ncid, varid, attr_name, NC_CHAR, len, timestamp_str))
       ||(0 != TIO_put_att (ncid, varid, attr_name_since_epoch, NC_DOUBLE, 1, &secs_since_file_epoch)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing timestamp %s",
                     __func__, attr_name);
        return -1;
     }

   return 0;
}

int TIO_write_scan_ident (int ncid, TIO_Scan_Ident_Type *lst)
{
   _pTIO_Granule_Ident_Type *beg=NULL, *end=NULL, *gid;
   double t_beg, t_end;
   size_t len;

   if (lst == NULL)
     return -1;

   beg = lst->granule_ident;
   end = beg;

   t_beg = beg->tstart;
   t_end = beg->tend;

   for (gid = lst->granule_ident; gid != NULL; gid = gid->next)
     {
        if (gid->tstart < t_beg)
          {
             t_beg = gid->tstart;
             beg = gid;
          }

        if (gid->tend > t_end)
          {
             t_end = gid->tend;
             end = gid;
          }
     }

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "scan_num", NC_INT, 1, &beg->scan_num))
     return -1;

   len = strlen (beg->tstart_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, len, beg->tstart_str))
     return -1;

   len = strlen (end->tend_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, len, end->tend_str))
     return -1;

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, 1, &t_beg))
     return -1;

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, 1, &t_end))
     return -1;

   return 0;
}

FCALLSCFUN2(INT, TIO_copy_granule_ident, TIO_F_COPY_GRANULE_IDENT, tio_f_copy_granule_ident,
            INT, INT)
FCALLSCFUN2(INT, TIO_same_granule_ident, TIO_F_SAME_GRANULE_IDENT, tio_f_same_granule_ident,
            INT, INT)
FCALLSCFUN6(INT, TIO_filename_from_granule, TIO_F_FILENAME_FROM_GRANULE, tio_f_filename_from_granule,
            INT, STRING, INT, INT, PSTRING, INT)
FCALLSCFUN3(INT, TIO_label_product, TIO_F_LABEL_PRODUCT, tio_f_label_product,
            INT, STRING, INT)
FCALLSCFUN5(INT, tio_time_tempo_to_utc_caldate, TIO_F_TEMPO_TIME_TO_UTC_CALDATE, tio_f_tempo_time_to_utc_caldate,
            DOUBLE, PINT,PINT,PINT,PDOUBLE)
