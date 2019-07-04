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

#define INSTR_PRIVATE_DATA \
   Instr_Type *next; \
   double *timestamp; \
   float *adc_temp0_derived; \
   float *fpe_temp1; \
   size_t num_times;
#include "instr.h"

typedef struct
{
   const char *glob_basename;
   double tstart;
   double tend;
}
Instr_Filter_Type;

static int find_entry1 (const Instr_Type *instr, double t, int *entry)
{
   double *ti = instr->timestamp;
   int n = instr->num_times;

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

static const Instr_Type *find_entry (const Instr_Type *instr, double t,
                                     int *entry, int *entry_status)
{
   const Instr_Type *closest_instr = NULL;
   int closest_status = -1;
   int closest_index = 0;
   double min_delta = DBL_MAX;

   /* We could merge the lists, but this is probably fast enough */

   for ( ; instr != NULL; instr = instr->next)
     {
        int index, status;
        status = find_entry1 (instr, t, &index);
        if (status == 0)
          {
             *entry = index;
             *entry_status = 0;
             return instr;
          }
        else if (fabs(t-instr->timestamp[index]) < min_delta)
          {
             closest_index = index;
             closest_instr = instr;
             closest_status = status;
          }
     }

   *entry_status = closest_status;
   *entry = closest_index;

   return closest_instr;
}

static int instr_temps (const Instr_Type *instr, double timestamp,
                        float *fpa_temp, float *fpe_temp)
{
   const Instr_Type *this_instr = NULL;
   int index, index_status;
   this_instr = find_entry (instr, timestamp, &index, &index_status);
   if (this_instr)
     {
        *fpa_temp = this_instr->adc_temp0_derived[index];
        *fpe_temp = this_instr->fpe_temp1[index];
        return index_status;
     }
   tell_verror (TELL_UNKNOWN_ERROR, "%s: temperature lookup failed", __func__);
   return -2;
}

static void free_instr1 (Instr_Type *instr)
{
   if (instr == NULL)
     return;
   FREE(instr->timestamp);
   FREE(instr->adc_temp0_derived);
   FREE(instr->fpe_temp1);
   FREE(instr);
}

static void free_instr (Instr_Type *instr)
{
   while (instr)
     {
        Instr_Type *next = instr->next;
        free_instr1 (instr);
        instr = next;
     }
}

static Instr_Type *new_instr_type (size_t num_times)
{
   Instr_Type *instr = NULL;

   if (NULL == (instr = (Instr_Type *) MALLOC (sizeof *instr)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)instr, 0, sizeof *instr);

   if ((NULL == (instr->timestamp = (double *) MALLOC (num_times * sizeof(double))))
       || (NULL == (instr->adc_temp0_derived = (float *) MALLOC (num_times * sizeof(float))))
       || (NULL == (instr->fpe_temp1 = (float *) MALLOC (num_times * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_instr1 (instr);
        return NULL;
     }

   instr->num_times = num_times;
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

static int time_sort1 (Instr_Type *instr)
{
   static int *sort_index = NULL;
   size_t i, num_keys = instr->num_times;

   if (NULL == (sort_index = (int *) MALLOC (instr->num_times * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < num_keys; i++)
     sort_index[i] = i;

   Qsort_Keys = instr->timestamp;
   Qsort_Indices = sort_index;
   num_keys = instr->num_times;
   qsort (sort_index, num_keys, sizeof(int), Qsort_Index_Compare);

   return 0;
}

static TIO_Meta_Type *_pMeta_Ptr;

static Instr_Type *read_instr1 (const char *file)
{
   Instr_Type *instr = NULL;
   int ncid, grp, dimid, start, count;
   size_t num_times;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if (0 != TIO_inq_grp (ncid, "fpe_analog_tlm2", &grp))
     goto return_error;

   if (0 != TIO_inq_dim (grp, "time", &dimid, &num_times))
     goto return_error;

   if (NULL == (instr = new_instr_type (num_times)))
     goto return_error;

   start = 0;
   count = num_times;

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

   if ((0 != TIO_get_var_section (grp, "time", &start, &count, TIO_DOUBLE,
                                  instr->timestamp))
       || (0 != TIO_get_var_section (grp, "adc_temp0_derived", &start, &count, TIO_FLOAT,
                                     instr->adc_temp0_derived))
       || (0 != TIO_get_var_section (grp, "fpe_temp1", &start, &count, TIO_FLOAT,
                                     instr->fpe_temp1))
      )
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading instrument status: %s",
                     __func__, file);
        goto return_error;
     }

   if (-1 == time_sort1 (instr))
     goto return_error;

   if (0 != meta_record_basename (_pMeta_Ptr, file))
     goto return_error;

   TIO_close (ncid);

   return instr;

return_error:
   TIO_close (ncid);
   free_instr1 (instr);
   return NULL;
}

static Instr_Type *read_instr_list (const char *path)
{
   FILE *fp = NULL;
   Instr_Type *head = NULL;

   if (NULL == (fp = fopen (path, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s", __func__, path);
        return NULL;
     }

   for (;;)
     {
        Instr_Type *instr, **tail;
        char *newline;
        char buf[1024];

        if (NULL == fgets (buf, sizeof(buf), fp))
          break;
        /* Assume no leading whitespace, no line-breaks within file path,
         * each line ends with a single newline.
         */
        if (NULL != (newline = strchr (buf, '\n')))
          *newline = 0;

        if (NULL == (instr = read_instr1 (buf)))
          {
             free_instr (instr);
             break;
          }

        /* Preserve the order, in case the files are sorted */
        if (head == NULL)
          head = instr;
        else
          {
             *tail = instr;
          }
        tail = &instr->next;
     }

   fclose (fp);

   return head;
}

enum
{
   FILTER_ERROR = -1,
   FILTER_INCLUDES_FILE = 0,
   FILTER_EXCLUDES_FILE = 1
};

static int filter_excludes_file (const Instr_Filter_Type *flt, const char *file)
{
   double file_tstart, file_tend;
   int ncid, status;

   status = FILTER_ERROR;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return FILTER_ERROR;

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &file_tstart))
     goto return_status;

   if (flt->tend < file_tstart)
     {
        status = FILTER_EXCLUDES_FILE;
        goto return_status;
     }

   if (-1 == TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, &file_tend))
     goto return_status;

   if (file_tend < flt->tstart)
     {
        status = FILTER_EXCLUDES_FILE;
        goto return_status;
     }

   status = FILTER_INCLUDES_FILE;

return_status:
   (void) TIO_close (ncid);
   return status;
}

static Instr_Type *read_instr_glob (const char *path, const Instr_Filter_Type *flt)
{
   IOCLib_Glob_Type *g = NULL;
   Instr_Type *head = NULL;
   char *glob_path = NULL;
   unsigned int i;
   int status_flag = -1;

   if (flt->glob_basename == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: glob_basename == NULL", __func__);
        return NULL;
     }

   if (NULL == (glob_path = ioclib_pathconcat (path, flt->glob_basename)))
     goto return_status;

   /* The globbing pattern is assumed to yield a time-ordered list of files */
   if ((NULL == (g = ioclib_glob (glob_path, 0)))
       || (g->num_files == 0))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: no matching files: %s", __func__, glob_path);
	goto return_status;
     }

   for (i = 0; i < g->num_files; i++)
     {
        Instr_Type *instr, **tail;
        int status = filter_excludes_file (flt, g->files[i]);

        if (status < 0)
          goto return_status;
        else if (status > 0)
          continue;

        if (NULL == (instr = read_instr1 (g->files[i])))
          goto return_status;

        /* Preserve the order */
        if (head == NULL)
          head = instr;
        else
          {
             *tail = instr;
          }
        tail = &instr->next;
     }

   status_flag = 0;
return_status:
   if (status_flag)
     {
	tell_verror (TELL_RUNTIME_ERROR, "%s: failed", __func__);
     }

   ioclib_free (glob_path);
   ioclib_glob_free (g);

   if (status_flag)
     {
        free_instr (head);
        head = NULL;
     }

   return head;
}

static Instr_Type *read_instr (const char *path, const Instr_Filter_Type *flt)
{
   struct stat st = {0};

   /* The input path may represent one of the following
    * alternatives:
    * 1) path = '@LISTFILE'
    *           where LISTFILE is the path to a file containing
    *           a time-ordered list of filenames
    * 2) path = 'DIRPATH'
    *           where DIRPATH is the path to a directory containing
    *           files matching a globbing expression
    * 3) path = 'FILENAME'
    *           where FILENAME is the path to a single file
    */

   if (*path == '@')
     {
        path++;
        return read_instr_list (path);
     }

   if (stat(path, &st) == -1)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: cannot stat %s",
                     __func__, path ? path : "<null>");
        return NULL;
     }

   if (S_ISDIR(st.st_mode))
     {
        return read_instr_glob (path, flt);
     }

   return read_instr1 (path);
}

Instr_Type *instr_open (const char *file, const char *glob_basename,
                        double tstart, double tend, TIO_Meta_Type *meta)
{
   Instr_Filter_Type flt = {0};
   Instr_Type *instr = NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: starting", __func__);

   if (file == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: received file == NULL",
                     __func__);
        return NULL;
     }

   flt.glob_basename = glob_basename;
   flt.tstart = tstart;
   flt.tend = tend;

   _pMeta_Ptr = meta;

   if (NULL == (instr = read_instr (file, &flt)))
     return NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: succeeded", __func__);

   return instr;
}
