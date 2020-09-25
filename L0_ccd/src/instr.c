#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <tell.h>
#include <tio.h>
#include <ioclib.h>

#include "config.h"
#include "util.h"

typedef struct TP_Type TP_Type;

struct TP_Type
{
   TP_Type *next;
   double *timestamp;
   float *value;
   size_t num_times;
};

#define INSTR_PRIVATE_DATA \
   TP_Type *adc_temp0_derived; \
   TP_Type *fpe_temp1;
#include "instr.h"

typedef struct
{
   const char *glob_basename;
   double tstart;
   double tend;
}
Instr_Filter_Type;

static int find_entry1 (const TP_Type *tp, double t, int *entry)
{
   double *ti = tp->timestamp;
   int n = tp->num_times;

   if (ti[0] <= t && t <= ti[n-1])
     {
        *entry = bsearch_d (t, ti, n);
        return 0;
     }
   else if (t < ti[0])
     {
        *entry = 0;
        return -1;
     }
   else
     {
        *entry = n-1;
        return +1;
     }
}

static const TP_Type *find_entry (const TP_Type *tp, double t,
                                  int *entry, int *entry_status)
{
   const TP_Type *closest_tp = NULL;
   int closest_status = -1;
   int closest_index = 0;
   double min_delta = DBL_MAX;

   /* We could merge the lists, but this is probably fast enough */

   for ( ; tp != NULL; tp = tp->next)
     {
        int index, status;
        status = find_entry1 (tp, t, &index);
        if (status == 0)
          {
             *entry = index;
             *entry_status = 0;
             return tp;
          }
        else if (fabs(t-tp->timestamp[index]) < min_delta)
          {
             closest_index = index;
             closest_tp = tp;
             closest_status = status;
          }
     }

   *entry_status = closest_status;
   *entry = closest_index;

   return closest_tp;
}

/* FIXME: should get FPE and FPA temps separately, so the status can be returned separately for each */
static int instr_temps (const Instr_Type *instr, double timestamp,
                        float *fpa_temp, float *fpe_temp)
{
   const TP_Type *tp = NULL;
   int index, index_status;

   if (NULL == (tp = find_entry (instr->adc_temp0_derived, timestamp, &index, &index_status)))
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: FPA temperature lookup failed", __func__);
        return -2;
     }
   *fpa_temp = tp->value[index];
   tell_vlog (TELL_MSGTYPE_INFO, 2, "FPA temp lookup: time: %0.3f => FPA temp: %0.2f C", timestamp, *fpa_temp);

   if (NULL == (tp = find_entry (instr->fpe_temp1, timestamp, &index, &index_status)))
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: FPE temperature lookup failed", __func__);
        return -2;
     }
   *fpe_temp = tp->value[index];
   tell_vlog (TELL_MSGTYPE_INFO, 2, "FPE temp lookup: time: %0.3f => FPE temp: %.2f C", timestamp, *fpe_temp);

   return 0;
}

static void free_tp1 (TP_Type *tp)
{
   if (tp == NULL)
     return;
   FREE(tp->timestamp);
   FREE(tp->value);
   FREE(tp);
}

static void free_tp (TP_Type *tp)
{
   while (tp)
     {
        TP_Type *next = tp->next;
        free_tp1(tp);
        tp = next;
     }
}

static TP_Type *alloc_tp (size_t num_times)
{
   TP_Type *tp = NULL;

   if (NULL == (tp = (TP_Type *)MALLOC (sizeof *tp)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if ((NULL == (tp->timestamp = (double *) MALLOC (num_times * sizeof(double))))
       || (NULL == (tp->value = (float *) MALLOC (num_times * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_tp1 (tp);
        return NULL;
     }

   tp->num_times = num_times;
   tp->next = NULL;

   return tp;
}

static int append_tp (TP_Type **head, TP_Type *tp)
{
   TP_Type *p;

   if (tp == NULL)
     return -1;

   if (head == NULL)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: head pointer is NULL (should never happen!)", __func__);
        return -1;
     }

   if (*head == NULL)
     {
        *head = tp;
        return 0;
     }

   for (p = *head; p->next != NULL; p = p->next)
     {}

   p->next = tp;

   return 0;
}

static void free_instr (Instr_Type *instr)
{
   if (instr == NULL)
     return;
   free_tp (instr->adc_temp0_derived);
   free_tp (instr->fpe_temp1);
   FREE(instr);
}

static Instr_Type *new_instr_type (void)
{
   Instr_Type *instr = NULL;

   if (NULL == (instr = (Instr_Type *) MALLOC (sizeof *instr)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)instr, 0, sizeof *instr);

   instr->instr_delete = free_instr;
   instr->instr_temps = instr_temps;

   return instr;
}

static double *Qsort_Keys;
static int *Qsort_Indices;
static int Qsort_Index_Compare (const void *a, const void *b)
{
   double va = Qsort_Keys[Qsort_Indices[*(int *)a]];
   double vb = Qsort_Keys[Qsort_Indices[*(int *)b]];
   if (va < vb) return -1;
   else if (va > vb) return +1;
   else return 0;
}

static int time_sort1 (TP_Type *tp)
{
   int *sort_index = NULL;
   size_t i, num_keys = tp->num_times;

   if (NULL == (sort_index = (int *) MALLOC (tp->num_times * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < num_keys; i++)
     sort_index[i] = i;

   Qsort_Keys = tp->timestamp;
   Qsort_Indices = sort_index;
   num_keys = tp->num_times;
   qsort (sort_index, num_keys, sizeof(int), Qsort_Index_Compare);

   FREE(sort_index);
   sort_index = NULL;

   return 0;
}

static TIO_Meta_Type *_pMeta_Ptr;

static TP_Type *read_telemetry_point (int ncid, const char *grp_name,
                                      const char *var_name)
{
   TP_Type *tp = NULL;
   int group_exists, grp, dimid, start, count, varid;
   size_t num_times;

   tell_push_queue();
   group_exists = (0 == TIO_inq_grp (ncid, grp_name, &grp));
   tell_pop_queue (1);
   if (0 == group_exists)
     return NULL;

   if (NC_NOERR != tio_inq_varid (grp, var_name, &varid))
     {
        tell_vwarn (0, "%s: section %s: variable %s not present\n", __func__, grp_name, var_name);
        return NULL;
     }

   if (0 != TIO_inq_dim (grp, "time", &dimid, &num_times))
     return NULL;

   if (NULL == (tp = alloc_tp (num_times)))
     return NULL;

   start = 0;
   count = num_times;

   if ((0 != TIO_get_var_section (grp, "time", &start, &count, TIO_DOUBLE, tp->timestamp))
       || (0 != TIO_get_var_section (grp, var_name, &start, &count, TIO_FLOAT, tp->value)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading instrument telemetry point %s",  __func__, var_name);
        goto free_and_return;
     }

   /* FIXME: probably don't need to do this */
   if (0 != time_sort1 (tp))
     goto free_and_return;

   return tp;
free_and_return:
   free_tp1 (tp);
   return NULL;
}

static int read_instr1 (Instr_Type *instr, const char *file)
{
   /* FIXME - after fixing test data, stop looking for adc_temp0_derived in fpe_analog_tlm2
    * Some of the level 0 test data erroneously put adc_temp0_derived in fpe_analog_tlm2
    * instead of in adc_tlm (mea culpa).  In the near term, it's simplest to have the code
    * look in adc_tlm first, and then check for fpe_analog_tlm2 as a fallback.
    * With real data this should not be necessary.  The problem should be corrected in
    * the next set of test data.  Once this workaround is no longer needed, remove it!
    */
   char *adc_temp0_sections[] = {"adc_tlm", "fpe_analog_tlm2", NULL};
   char **section;
   int ncid, found_values = 0;
   TP_Type *tp;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   /* For the FPA temperature, the Command and Telemetry Handbook contains
    * at least 5 telemetry points that are related. Dave Flittner (NASA/LARC)
    * tells me that the ADC_TEMP0 is the one to use.
    * ADC_TEMP0 is used as the control point in the thermal control of the CCD.
    * It is located on the conduction bar between the S-Link, which is connected
    * to the FPA, and the FPA thermal interface with the Host S/C. In ground
    * testing, the on-die temperature sense resistors, CCD_TEMP1 and CCD_TEMP2
    * were noisy and sometimes erratic, while ADC_TEMP0 was most reliable.
    *
    * For the FPE temperature, FPE_TEMP1 is the only option.
    *
    * Regarding units, ADC_TEMP0 is actually in Ohms, so for a Celsius temperature,
    * we use ADC_TEMP0_DERIVED, which is provided by the IOC.  FPE_TEMP1 is already
    * in Celsius so we use that directly.
    */

   for (section = adc_temp0_sections; section != NULL; section++)
     {
        if (NULL != (tp = read_telemetry_point (ncid, *section, "adc_temp0_derived")))
          {
             if (0 != append_tp (&instr->adc_temp0_derived, tp))
               goto return_error;
             found_values++;
             break;
          }
     }

   if (NULL != (tp = read_telemetry_point (ncid, "fpe_analog_tlm2", "fpe_temp1")))
     {
        if (0 != append_tp (&instr->fpe_temp1, tp))
          goto return_error;
        found_values++;
     }

   if (found_values)
     {
        if (0 != meta_record_basename (_pMeta_Ptr, file))
          goto return_error;
     }

   TIO_close (ncid);

   return 0;

return_error:
   TIO_close (ncid);
   return -1;
}

static int read_instr_list (Instr_Type *instr, const char *path)
{
   FILE *fp = NULL;

   if (NULL == (fp = fopen (path, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s", __func__, path);
        return -1;
     }

   for (;;)
     {
        char *newline;
        char buf[1024];

        if (NULL == fgets (buf, sizeof(buf), fp))
          break;
        /* Assume no leading whitespace, no line-breaks within file path,
         * each line ends with a single newline.
         */
        if (NULL != (newline = strchr (buf, '\n')))
          *newline = 0;

        if (0 != read_instr1 (instr, buf))
          {
             fclose (fp);
             return -1;
          }
     }

   fclose (fp);

   return 0;
}

enum
{
   FILTER_ERROR = -1,
   FILTER_INCLUDES_FILE = 0,
   FILTER_EXCLUDES_FILE = 1,
   FILTER_HALTS_SEARCH = 2
};

static int filter_excludes_file (const Instr_Filter_Type *flt, const char *file,
                                 int *flag, double *delta_t)
{
   double file_tstart, file_tend;
   int ncid;

   tell_vlog (TELL_MSGTYPE_INFO, 2, "%s: looking at %s", __func__, file);

   *flag = FILTER_ERROR;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &file_tstart))
     {
        TIO_close(ncid);
        return -1;
     }

   if (flt->tend < file_tstart)
     {
        *flag = FILTER_HALTS_SEARCH;
        *delta_t = file_tstart - flt->tend;
        TIO_close(ncid);
        return 0;
     }

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, &file_tend))
     {
        TIO_close(ncid);
        return -1;
     }

   if (file_tend < flt->tstart)
     {
        *flag = FILTER_EXCLUDES_FILE;
        *delta_t = flt->tstart - file_tend;
     }
   else
     {
        *flag = FILTER_INCLUDES_FILE;
        *delta_t = 0.0;
     }

   (void) TIO_close (ncid);
   return 0;
}

static int read_instr_glob (Instr_Type *instr, const char *path, const Instr_Filter_Type *flt)
{
   IOCLib_Glob_Type *g = NULL;
   char *glob_path = NULL;
   unsigned int i;
   double min_delta_t = DBL_MAX;
   int nearest = -1;
   int status = -1;
   int num_files_selected = 0;

   if (flt->glob_basename == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: glob_basename == NULL", __func__);
        return -1;
     }

   if (NULL == (glob_path = ioclib_pathconcat (path, flt->glob_basename)))
     goto return_status;

   /* The globbing pattern is assumed to yield a time-ordered list of files */
   if ((NULL == (g = ioclib_glob (glob_path, 0)))
       || (g->num_files == 0))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: no files match pattern: %s", __func__, glob_path);
	goto return_status;
     }

   for (i = 0; i < g->num_files; i++)
     {
        double delta_t;
        int flag;

        if (filter_excludes_file (flt, g->files[i], &flag, &delta_t) < 0)
          goto return_status;

        if (delta_t < min_delta_t)
          {
             min_delta_t = delta_t;
             nearest = i;
          }

        if (flag == FILTER_ERROR)
          goto return_status;
        else if (flag == FILTER_EXCLUDES_FILE)
          continue;
        else if (flag == FILTER_HALTS_SEARCH)
          break;

        num_files_selected++;

        if (0 != read_instr1 (instr, g->files[i]))
          goto return_status;
     }

   if (num_files_selected == 0)
     {
        tell_vwarn (0, "%s: no files match selection criteria: %s",
                    __func__, glob_path);
        if (nearest >= 0)
          {
             if (0 != read_instr1 (instr, g->files[nearest]))
               goto return_status;
             tell_vwarn (0, "%s: using nearest match (delta_t=%f sec): %s",
                         __func__, min_delta_t, g->files[nearest]);
          }
     }

   status = 0;
return_status:
   if (status)
     {
	tell_verror (TELL_RUNTIME_ERROR, "%s: failed", __func__);
     }

   ioclib_free (glob_path);
   ioclib_glob_free (g);

   return status;
}

static char *make_hk_dir_path (int sat_day)
{
   char buf[1024];
   size_t bufsize = sizeof(buf);
   int len;

   /* If the environment variable is not set, look in the current directory */
   if (NULL == getenv ("SDPC_ARCHIVE_DIR"))
     {
        if (NULL == getcwd (buf, bufsize))
          return strdup (buf);
        else return NULL;
     }

   /* Construct path to archive directory containing the telemetry point
    * stream for this day.
    */
   len = snprintf (buf, bufsize, "$SDPC_ARCHIVE_DIR/L0/D%05d/HK", sat_day);
   if ((len < 0) || ((size_t) len >= bufsize))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: path exceeded buffer size", __func__);
        return NULL;
     }

   return expand_string (buf);
}

static int read_instr_default_paths (Instr_Type *instr, const Instr_Filter_Type *flt)
{
   char *dir = NULL;
   double day_beg, day_end;
   int iday_beg, iday_end;
   int status_beg, status_end;

   if ((0 != tio_time_sat_local_day_number (flt->tstart, &day_beg))
       ||(0 != tio_time_sat_local_day_number (flt->tend, &day_end)))
     return -1;

   iday_beg = (int) day_beg;
   iday_end = (int) day_end;

   /* Observations occuring near local midnight may have instrument
    * telemetry point files spread across two archive directories.
    * Radiance observations near midnight should never happen (violates mission safety rules).
    * Irradiance observations near midnight are possible, but unlikely.
    * Dark observations near midnight are allowed and must be supported.
    */

   if ((iday_end < iday_beg)
       || ((iday_end - iday_beg) > 1))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: invalid/unsupported observation interval: tstart=%f, tend=%f",
                     __func__, flt->tstart, flt->tend);
        return -1;
     }

   if (NULL == (dir = make_hk_dir_path (iday_beg)))
     return -1;

   status_beg = read_instr_glob (instr, dir, flt);
   FREE(dir);

   if (iday_end > iday_beg)
     {
        if (NULL == (dir = make_hk_dir_path (iday_end)))
          {
             free_instr (instr);
             return -1;
          }
        status_end = read_instr_glob (instr, dir, flt);
        FREE(dir);
     }

   if ((status_beg == 0) || (status_end == 0))
     return 0;

   return -1;
}

static int read_instr (Instr_Type *instr, const char *path, const Instr_Filter_Type *flt)
{
   struct stat st = {0};

   /* The input path may represent one of the following
    * alternatives:
    * 0) path = NULL
    *           => use the archive default, or, if SDPC_ARCHIVE_DIR
    *              is not set, look in the current directory for
    *              files matching the globbing expression
    * 1) path = '@LISTFILE'
    *           where LISTFILE is the path to a file containing
    *           a time-ordered list of filenames
    * 2) path = 'DIRPATH'
    *           where DIRPATH is the path to a directory containing
    *           files matching a globbing expression
    * 3) path = 'FILENAME'
    *           where FILENAME is the path to a single file
    */

   if (path == NULL)
     {
        return read_instr_default_paths (instr, flt);
     }

   if (*path == '@')
     {
        path++;
        return read_instr_list (instr, path);
     }

   if (stat(path, &st) == -1)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: cannot stat %s",
                     __func__, path ? path : "<null>");
        return -1;
     }

   if (S_ISDIR(st.st_mode))
     {
        return read_instr_glob (instr, path, flt);
     }

   return read_instr1 (instr, path);
}

Instr_Type *instr_open (const char *file, const char *glob_basename,
                        double tstart, double tend, TIO_Meta_Type *meta)
{
   Instr_Filter_Type flt = {0};
   Instr_Type *instr = NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: starting", __func__);

   flt.glob_basename = glob_basename;
   flt.tstart = tstart;
   flt.tend = tend;

   _pMeta_Ptr = meta;

   if (NULL == (instr = new_instr_type ()))
     return NULL;

   if (0 != read_instr (instr, file, &flt))
     {
        free_instr (instr);
        return NULL;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: succeeded", __func__);

   return instr;
}
