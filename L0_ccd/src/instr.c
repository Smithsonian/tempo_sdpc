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
   double *timestamp;
   double *value;
   size_t num_times;
   size_t num_allocated;
};

#define INSTR_PRIVATE_DATA \
   TP_Type adc_temp0_derived; \
   TP_Type adc_temp5_derived; \
   TP_Type adc_temp8_derived; \
   TP_Type adc_temp16_derived; \
   TP_Type fpe_temp1;
#include "instr.h"

typedef struct
{
   const char *glob_basename;
   double tstart;
   double tend;
}
Instr_Filter_Type;

static int lookup_value (const TP_Type *tp, double t, double *value)
{
   double *ti = tp->timestamp;
   int n = tp->num_times;

   if (t <= ti[0])
     {
        *value = tp->value[0];
        return -1;
     }
   else if (ti[n-1] <= t)
     {
        *value = tp->value[n-1];
        return +1;
     }
   else
     {
        int k = bsearch_d (t, ti, n);
        *value = tp->value[k];
     }

   return 0;
}

static int instr_fpa_temp (const Instr_Type *instr, double timestamp, float *fpa_temp)
{
   double value;

   if (0 != lookup_value (&instr->adc_temp0_derived, timestamp, &value))
     {
        tell_vwarn (0, "%s: FPA temperature lookup failed", __func__);
     }
   *fpa_temp = (float) value;

   return 0;
}

static int instr_fpe_temp (const Instr_Type *instr, double timestamp, float *fpe_temp)
{
   double value;

   if (0 != lookup_value (&instr->fpe_temp1, timestamp, &value))
     {
        tell_vwarn (0, "%s: FPE temperature lookup failed", __func__);
     }
   *fpe_temp = (float) value;

   return 0;
}

static int instr_spec_temp (const Instr_Type *instr, double timestamp, float *spec_temp)
{
   double value;

   if (0 != lookup_value (&instr->adc_temp5_derived, timestamp, &value))
     {
        tell_vwarn (0, "%s: Spectrometer temperature lookup failed", __func__);
     }
   *spec_temp = (float) value;

   return 0;
}

static int instr_tele_temp (const Instr_Type *instr, double timestamp, float *tele_temp)
{
   double value;

   if (0 != lookup_value (&instr->adc_temp8_derived, timestamp, &value))
     {
        tell_vwarn (0, "%s: Telescope temperature lookup failed", __func__);
     }
   *tele_temp = (float) value;

   return 0;
}

static int instr_bench_temp (const Instr_Type *instr, double timestamp, float *bench_temp)
{
   double value;

   if (0 != lookup_value (&instr->adc_temp16_derived, timestamp, &value))
     {
        tell_vwarn (0, "%s: Optical bench temperature lookup failed", __func__);
     }
   *bench_temp = (float) value;

   return 0;
}

static void free_tp (TP_Type *tp)
{
   if (tp == NULL)
     return;
   FREE(tp->timestamp);
   FREE(tp->value);
}

static int realloc_tp (TP_Type *tp, size_t num_times)
{
   double *tstamp = NULL;
   double *value = NULL;
   size_t num_new = tp->num_allocated + num_times;

   if (NULL == (tstamp = (double *)REALLOC (tp->timestamp, num_new * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   tp->timestamp = tstamp;

   if (NULL == (value = (double *)REALLOC (tp->value, num_new * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   tp->value = value;
   tp->num_allocated = num_new;

   return 0;
}

static void free_instr (Instr_Type *instr)
{
   if (instr == NULL)
     return;
   free_tp (&instr->adc_temp0_derived);
   free_tp (&instr->fpe_temp1);
   free_tp (&instr->adc_temp5_derived);
   free_tp (&instr->adc_temp8_derived);
   free_tp (&instr->adc_temp16_derived);
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
   instr->instr_fpa_temp = instr_fpa_temp;
   instr->instr_fpe_temp = instr_fpe_temp;
   instr->instr_spec_temp = instr_spec_temp;
   instr->instr_tele_temp = instr_tele_temp;
   instr->instr_bench_temp = instr_bench_temp;

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

static int time_sort (TP_Type *tp)
{
   int *sort_index = NULL;
   double *tmp = NULL;
   size_t i, k, num_keys = tp->num_times;

   if ((NULL == (sort_index = (int *) MALLOC (tp->num_times * sizeof(int))))
       || (NULL == (tmp = (double *) MALLOC (tp->num_times * sizeof(double)))))
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

   memcpy ((char *)tmp, (char *)tp->timestamp, num_keys * sizeof(double));
   for (i = 0; i < num_keys; i++)
     {
        k = sort_index[i];
        tp->timestamp[i] = tmp[k];
     }

   memcpy ((char *)tmp, (char *)tp->value, num_keys * sizeof(double));
   for (i = 0; i < num_keys; i++)
     {
        k = sort_index[i];
        tp->value[i] = tmp[k];
     }

   FREE(sort_index);
   FREE(tmp);

   return 0;
}

static TIO_Meta_Type *_pMeta_Ptr;

static int read_telemetry_point (TP_Type *tp, int ncid, const char *grp_name,
                                 const char *var_name)
{
   int group_exists, grp, dimid, start, count, varid;
   double *timestamp;
   double *value;
   size_t num_times;

   tell_push_queue();
   group_exists = (0 == TIO_inq_grp (ncid, grp_name, &grp));
   tell_pop_queue (1);
   if (0 == group_exists)
     return 1;

   if (NC_NOERR != tio_inq_varid (grp, var_name, &varid))
     {
        tell_vwarn (0, "%s: section %s: variable %s not present\n", __func__, grp_name, var_name);
        return -1;
     }

   if (0 != TIO_inq_dim (grp, "time", &dimid, &num_times))
     return -1;

   if (0 != realloc_tp (tp, num_times))
     return -1;

   timestamp = tp->timestamp + tp->num_times;
   value = tp->value + tp->num_times;

   start = 0;
   count = num_times;

   if ((0 != TIO_get_var_section (grp, "time", &start, &count, TIO_DOUBLE, timestamp))
       || (0 != TIO_get_var_section (grp, var_name, &start, &count, TIO_DOUBLE, value)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading instrument telemetry point %s",  __func__, var_name);
        return -1;
     }

   tp->num_times += num_times;

   return 0;
}

static int read_instr1 (Instr_Type *instr, const char *file)
{
   int ncid, status, found_values = 0;

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
    *
    * [April 18, 2024]
    * Three variables were added for trending purposes:
    * adc_temp5_derived, adc_temp8_derived, and adc_temp16_derived.
    * They represent temperatures of Spectrometer -X, Telescope T4, and
    * TEL -X-Y (bench heater set points, according to Xiong), respectively.
    */

   if ((status = read_telemetry_point (&instr->adc_temp0_derived, ncid, "adc_tlm", "adc_temp0_derived")) < 0)
     goto return_error;
   if (status == 0)
     {
        found_values++;
     }

   if ((status = read_telemetry_point (&instr->fpe_temp1, ncid, "fpe_analog_tlm2", "fpe_temp1")) < 0)
     goto return_error;
   if (status == 0)
     {
        found_values++;
     }

   if ((status = read_telemetry_point (&instr->adc_temp5_derived, ncid, "adc_tlm", "adc_temp5_derived")) < 0)
     goto return_error;
   if (status == 0)
     {
        found_values++;
     }

   if ((status = read_telemetry_point (&instr->adc_temp8_derived, ncid, "adc_tlm", "adc_temp8_derived")) < 0)
     goto return_error;
   if (status == 0)
     {
        found_values++;
     }

   if ((status = read_telemetry_point (&instr->adc_temp16_derived, ncid, "adc_tlm", "adc_temp16_derived")) < 0)
     goto return_error;
   if (status == 0)
     {
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
   else status_end = -1;

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

#define FILTER_WINDOW_HALF  (10)
#define FILTER_WINDOW_SIZE  (2*FILTER_WINDOW_HALF + 1)

static int fill_window (double *w, size_t nw, size_t i, const double *x, size_t n)
{
   size_t k, h = (nw-1)/2;

   for (k = 0; k < nw; k++)
     {
        int j = k + i-h;
        if (j < 0)
          {
             j = 0;
          }
        else if ((int)n <= j)
          {
             j = n-1;
          }
        w[k] = x[j];
     }

   return 0;
}

static int dbl_compare (const void *a, const void *b)
{
   double va = *(const double *)a;
   double vb = *(const double *)b;
   if (va < vb) return -1;
   else if (va > vb) return +1;
   else return 0;
}

static double median (double *w, size_t nw)
{
   qsort ((char *)w, nw, sizeof(double), dbl_compare);
   return w[(nw-1)/2];
}

static int median_filter (const double *x, size_t n, double *x_med)
{
   double w[FILTER_WINDOW_SIZE];
   size_t i, nw = FILTER_WINDOW_SIZE;

   for (i = 0; i < n; i++)
     {
        fill_window (w, nw, i, x, n);
        x_med[i] = median (w, nw);
     }

   return 0;
}

static int apply_median_filter (TP_Type *tp)
{
   double *tmp = NULL;
   size_t n = tp->num_times;

   if (NULL == (tmp = (double *) MALLOC (n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   (void) median_filter (tp->value, n, tmp);
   memcpy ((char *)tp->value, (char *)tmp, n * sizeof(double));

   FREE(tmp);

   return 0;
}

static int apply_filters (TP_Type *tp)
{
   if ((tp->timestamp == NULL)
       || (tp->value == NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: NULL telemetry point array", __func__);
        return -1;
     }

   /* In principle, sorting should not be necessary */
   if (0 != time_sort (tp))
     return -1;

   if (0 != apply_median_filter (tp))
     return -1;

   return 0;
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

   if ((0 != apply_filters (&instr->adc_temp0_derived))
       || (0 != apply_filters (&instr->fpe_temp1))
       || (0 != apply_filters (&instr->adc_temp5_derived))
       || (0 != apply_filters (&instr->adc_temp8_derived))
       || (0 != apply_filters (&instr->adc_temp16_derived)))
     {
        free_instr (instr);
        return NULL;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: succeeded", __func__);

   return instr;
}
