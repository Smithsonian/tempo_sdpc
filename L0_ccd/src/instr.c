#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <tell.h>
#include <tio.h>

#include "config.h"
#include "util.h"

#define INSTR_PRIVATE_DATA \
   Instr_Type *next; \
   double *timestamp; \
   float *ccd_temp1; \
   float *ccd_temp2; \
   size_t num_times;
#include "instr.h"

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

static int instr_ccd_temp1 (const Instr_Type *instr, double timestamp,
                            float *ccd_temp1)
{
   const Instr_Type *this_instr = NULL;
   int index, index_status;
   this_instr = find_entry (instr, timestamp, &index, &index_status);
   if (this_instr)
     {
        *ccd_temp1 = this_instr->ccd_temp1[index];
        return index_status;
     }
   tell_verror (TELL_UNKNOWN_ERROR, "%s: CCD_TEMP1 lookup failed", __func__);
   return -2;
}

static int instr_ccd_temp2 (const Instr_Type *instr, double timestamp,
                            float *ccd_temp2)
{
   const Instr_Type *this_instr = NULL;
   int index, index_status;
   this_instr = find_entry (instr, timestamp, &index, &index_status);
   if (this_instr)
     {
        *ccd_temp2 = this_instr->ccd_temp2[index];
        return index_status;
     }
   tell_verror (TELL_UNKNOWN_ERROR, "%s: CCD_TEMP2 lookup failed", __func__);
   return -2;
}

static void free_instr1 (Instr_Type *instr)
{
   if (instr == NULL)
     return;
   FREE(instr->timestamp);
   FREE(instr->ccd_temp1);
   FREE(instr->ccd_temp2);
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
       || (NULL == (instr->ccd_temp1 = (float *) MALLOC (num_times * sizeof(float))))
       || (NULL == (instr->ccd_temp2 = (float *) MALLOC (num_times * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_instr1 (instr);
        return NULL;
     }

   instr->num_times = num_times;
   instr->instr_delete = free_instr;
   instr->instr_ccd_temp1 = instr_ccd_temp1;
   instr->instr_ccd_temp2 = instr_ccd_temp2;

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

   if ((0 != TIO_get_var_section (grp, "time", &start, &count, TIO_DOUBLE,
                                  instr->timestamp))
       || (0 != TIO_get_var_section (grp, "ccd_temp1", &start, &count, TIO_FLOAT,
                                     instr->ccd_temp1))
       || (0 != TIO_get_var_section (grp, "ccd_temp2", &start, &count, TIO_FLOAT,
                                     instr->ccd_temp2)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading instrument status: %s",
                     __func__, file);
        goto return_error;
     }

   if (-1 == time_sort1 (instr))
     goto return_error;

   TIO_close (ncid);
   return instr;

return_error:
   TIO_close (ncid);
   free_instr1 (instr);
   return NULL;
}

static Instr_Type *read_instr (const char *path)
{
   FILE *fp = NULL;
   Instr_Type *head = NULL;

   if (*path != '@')
     {
        return read_instr1 (path);
     }

   path++;

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

Instr_Type *instr_open (const char *file)
{
   if (file == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: received file == NULL",
                     __func__);
        return NULL;
     }

   return read_instr (file);
}
