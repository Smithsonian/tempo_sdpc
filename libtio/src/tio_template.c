/** @file
 *  @brief TEMPO-specific utility functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define __USE_XOPEN
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

int _pTIO_read_granule_ident (int ncid, _pTIO_Granule_Ident_Type *gid)
{
   if ((-1 == TIO_get_att (ncid, NC_GLOBAL, "scan_seq_num", NC_INT, &gid->scan_seq_num))
       ||(-1 == TIO_get_att (ncid, NC_GLOBAL, "granule_seq_num", NC_INT, &gid->granule_seq_num))
       ||(-1 == TIO_get_att (ncid, NC_GLOBAL, "granule_num", NC_INT, &gid->granule_num)))
     return -1;

   memset (gid->tstart_str, 0, MAX_ISOTIME_LEN);
   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, gid->tstart_str))
     return -1;

   memset (gid->tend_str, 0, MAX_ISOTIME_LEN);
   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, gid->tend_str))
     return -1;

   return 0;
}

static int
_pTIO_write_granule_ident (int ncid, const _pTIO_Granule_Ident_Type *gid)
{
   size_t len;

   if ((-1 == TIO_put_att (ncid, NC_GLOBAL, "scan_seq_num", NC_INT, 1, &gid->scan_seq_num))
       ||(-1 == TIO_put_att (ncid, NC_GLOBAL, "granule_seq_num", NC_INT, 1, &gid->granule_seq_num))
       ||(-1 == TIO_put_att (ncid, NC_GLOBAL, "granule_num", NC_INT, 1, &gid->granule_num)))
     return -1;

   len = strlen (gid->tstart_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, len, gid->tstart_str))
     return -1;

   len = strlen (gid->tend_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, len, gid->tend_str))
     return -1;

   return 0;
}

int _pTIO_parse_timestr (const char *timestr, struct tm *ptm)
{
   memset ((char *)ptm, 0, sizeof (struct tm));
   if (NULL == strptime (timestr, TIO_DELIM_TIMESTAMP_FORMAT, ptm))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: strptime failed: %s",
                     __func__, timestr);
        return -1;
     }

   return 0;
}

static int
_pTIO_filename_from_granule_ident (const _pTIO_Granule_Ident_Type *gid,
                                   const char *label, int version,
                                   char *buf, int bufsize)
{
   char timestr[MAX_ISOTIME_LEN];
   struct tm tm;
   int status;

   if (-1 == _pTIO_parse_timestr (gid->tstart_str, &tm))
     return -1;

   if (0 == strftime (timestr, MAX_ISOTIME_LEN,
                      TIO_NODELIM_TIMESTAMP_FORMAT, &tm))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: strftime failed",
                     __func__);
        return -1;
     }

   /* 6 digit sequence number should be enough for a 10 year mission:
    * (10 years)*(365 days/year)*(100 scans/day) = 3.65e5 scans
    * which fits in a 6 digit int.
    */
   status = snprintf (buf, bufsize,
                      "tempo_%s_%06d_%02d_v%d_%s.nc",
                      timestr,
                      gid->scan_seq_num,
                      gid->granule_num,
                      version,
                      label);
   return status;
}

int TIO_copy_granule_ident (int ncid_from, int ncid_to)
{
   _pTIO_Granule_Ident_Type gid;

   if ((-1 == _pTIO_read_granule_ident (ncid_from, &gid))
       || (-1 == _pTIO_write_granule_ident (ncid_to, &gid)))
     return -1;

   return 0;
}

static int same_granule_ident (const _pTIO_Granule_Ident_Type *gid1,
                               const _pTIO_Granule_Ident_Type *gid2)
{
   return ((gid1->scan_seq_num == gid2->scan_seq_num)
           && (gid1->granule_seq_num == gid2->granule_seq_num)
           && (gid1->granule_num == gid2->granule_num)
           && (0 == strncmp (gid1->tstart_str, gid2->tstart_str, MAX_ISOTIME_LEN))
           && (0 == strncmp (gid1->tend_str, gid2->tend_str, MAX_ISOTIME_LEN)));
}

int TIO_same_granule_ident (int ncid1, int ncid2)
{
   _pTIO_Granule_Ident_Type gid1, gid2;

   if ((-1 == _pTIO_read_granule_ident (ncid1, &gid1))
       || (-1 == _pTIO_read_granule_ident (ncid2, &gid2)))
     return -1;

   return same_granule_ident (&gid1, &gid2);
}

int TIO_filename_from_granule (int ncid, const char *label, int version,
                               char *buf, int bufsize)
{
   _pTIO_Granule_Ident_Type gid;

   if (-1 == _pTIO_read_granule_ident (ncid, &gid))
     return -1;

   return _pTIO_filename_from_granule_ident (&gid, label, version, buf, bufsize);
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

   if (item->scan_seq_num == lst->granule_ident->scan_seq_num)
     {
        status = 0;
     }
   else
     {
        Tell_verror (TELL_APPLICATION_ERROR,
                     "%s: scan_seq_num mismatch: new item has scan_seq_num=%d  granule list has scan_seq_num=%d",
                     __func__,
                     item->scan_seq_num,
                     lst->granule_ident->scan_seq_num);
        status = -1;
     }

   free_granule_ident (item);

   return status;
}

static int timet_from_timestr (const char *timestr, time_t *ptimet)
{
   const char tz_env[] = "TZ";
   const char *tz = NULL;
   struct tm tm;

   if (-1 == _pTIO_parse_timestr (timestr, &tm))
     return -1;

   tz = getenv (tz_env);
   setenv (tz_env, "", 1);     /* force UTC */
   *ptimet = mktime(&tm);
   if (tz)
     setenv (tz_env, tz, 1);
   else
     unsetenv (tz_env);

   return 0;
}

int TIO_write_scan_ident (int ncid, TIO_Scan_Ident_Type *lst)
{
   _pTIO_Granule_Ident_Type *beg=NULL, *end=NULL, *gid;
   time_t tt_beg=LONG_MAX, tt_end=0;
   size_t len;

   if (lst == NULL)
     return -1;

   beg = lst->granule_ident;
   end = beg;

   for (gid = lst->granule_ident; gid != NULL; gid = gid->next)
     {
        time_t tt;

        if (-1 == timet_from_timestr (gid->tstart_str, &tt))
          return -1;
        if (tt < tt_beg)
          {
             tt_beg = tt;
             beg = gid;
          }

        if (-1 == timet_from_timestr (gid->tend_str, &tt))
          return -1;
        if (tt > tt_end)
          {
             tt_end = tt;
             end = gid;
          }
     }

   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "scan_seq_num", NC_INT, 1, &beg->scan_seq_num))
     return -1;

   len = strlen (beg->tstart_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_start", NC_CHAR, len, beg->tstart_str))
     return -1;

   len = strlen (end->tend_str) + 1;
   if (-1 == TIO_put_att (ncid, NC_GLOBAL, "time_coverage_end", NC_CHAR, len, end->tend_str))
     return -1;

   return 0;
}

FCALLSCFUN2(INT, TIO_copy_granule_ident, TIO_F_COPY_GRANULE_IDENT, tio_f_copy_granule_ident,
            INT, INT)
FCALLSCFUN2(INT, TIO_same_granule_ident, TIO_F_SAME_GRANULE_IDENT, tio_f_same_granule_ident,
            INT, INT)
FCALLSCFUN5(INT, TIO_filename_from_granule, TIO_F_FILENAME_FROM_GRANULE, tio_f_filename_from_granule,
            INT, STRING, INT, PSTRING, INT)
FCALLSCFUN3(INT, TIO_label_product, TIO_F_LABEL_PRODUCT, tio_f_label_product,
            INT, STRING, INT)
