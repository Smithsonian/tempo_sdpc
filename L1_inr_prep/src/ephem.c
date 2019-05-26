/** @file ephem.c
 *  @brief Manage ephemeris subsetting
 */
#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>

#include <tell.h>

#include "ephem.h"
#include "bsearch.h"

#define BUFSIZE      256
#define COMMENT_CHAR '#'

static void free_vector (Eph_Vector_Type *v)
{
   FREE(v->x); v->x = NULL;
   FREE(v->y); v->y = NULL;
   FREE(v->z); v->z = NULL;
}

static double *alloc_doubles (size_t num_alloc)
{
   double *x = NULL;
   if (NULL == (x = (double *)MALLOC (num_alloc * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   return x;
}

static int alloc_vector (Eph_Vector_Type *v, size_t num_alloc)
{
   v->x = NULL; v->y = NULL; v->z = NULL;

   if ((NULL == (v->x = alloc_doubles (num_alloc)))
       ||(NULL == (v->y = alloc_doubles (num_alloc)))
       ||(NULL == (v->z = alloc_doubles (num_alloc))))
     {
        free_vector (v);
        return -1;
     }

   return 0;
}

void eph_free (Eph_Type *eph)
{
   free_vector (&eph->r);
   free_vector (&eph->v);
   FREE(eph->t);
   eph->t = NULL;
   eph->n = 0;
   eph->num_alloc = 0;
}

static int alloc_eph (Eph_Type *eph, size_t num_alloc)
{
   eph->t = NULL;

   if ((0 != alloc_vector (&eph->r, num_alloc))
       || (0 != alloc_vector (&eph->v, num_alloc))
       || (NULL == (eph->t = alloc_doubles (num_alloc))))
     {
        eph_free (eph);
        return -1;
     }

   eph->n = 0;
   eph->num_alloc = num_alloc;

   return 0;
}

static int count_times (const char *file)
{
   FILE *fp = NULL;
   int num_times;
   double t0 = 0.0;

   if (NULL == (fp = fopen (file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading",
                     __func__, file);
        return -1;
     }

   num_times = 0;
   while (0 == feof(fp))
     {
        char line_buf[BUFSIZE];
        double t;
        if (NULL == fgets (line_buf, BUFSIZE, fp))
          break;
        if ((line_buf[0] == COMMENT_CHAR) || (line_buf[0] == '\n'))
          continue;
        if (NULL == strchr (line_buf, '\n'))
          {
             tell_verror (TELL_IO_READ_ERROR,
                          "%s: input lines exceed internal fixed buffer size!",
                          __func__);
             return -1;
          }
        if (1 != sscanf (line_buf, "%le", &t))
          {
             tell_verror (TELL_IO_READ_ERROR,
                          "%s: unexpected value on non-comment line: %s",
                          __func__, line_buf);
             return -1;
          }
        if (t0 >= t)
          {
             tell_verror (TELL_INVALID_DATA_ERROR,
                          "%s: unexpected ephemeris time ordering (expected monotonic increasing time order)",
                          __func__);
             return -1;
          }
        else t0 = t;
        num_times++;
     }

   (void) fclose (fp);

   return num_times;
}

static int read_eph (const char *file, size_t num_times, Eph_Type *eph)
{
   const char fmt_csv_doubles[] = " %le , %le , %le , %le , %le , %le , %le";
   Eph_Vector_Type *r, *v;
   FILE *fp = NULL;
   double *t;

   if (NULL == (fp = fopen (file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading",
                     __func__, file);
        return -1;
     }

   t = eph->t;
   r = &eph->r;
   v = &eph->v;

   eph->n = 0;
   while ((eph->n < eph->num_alloc) && (0 == feof(fp)))
     {
        char line_buf[BUFSIZE];
        double tn, rx,ry,rz, vx,vy,vz;
        size_t n = eph->n;

        if (NULL == fgets (line_buf, BUFSIZE, fp))
          break;

        if ((line_buf[0] == COMMENT_CHAR) || (line_buf[0] == '\n'))
          continue;

        if (7 != sscanf (line_buf, fmt_csv_doubles,
                         &tn, &rx, &ry, &rz, &vx, &vy, &vz))
          {
             tell_verror (TELL_IO_READ_ERROR,
                          "%s: unexpected value on non-comment line: %s",
                          __func__, line_buf);
             return -1;
          }

        t[n] = tn;
        r->x[n] = rx;
        r->y[n] = ry;
        r->z[n] = rz;
        v->x[n] = vx;
        v->y[n] = vy;
        v->z[n] = vz;
        eph->n++;
     }

   (void) fclose (fp);

   if (eph->n != num_times)
     {
        tell_verror (TELL_IO_READ_ERROR,
                     "%s: read %ld values, expecting %ld",
                     __func__, eph->n, num_times);
        return -1;
     }

   return 0;
}

static int select_interval (Eph_Type *eph,
                            double time_beg, double time_end,
                            int num_pad)
{
   double eph_beg, eph_end;
   int beg, end;

   if (eph->n == 0)
     {
        tell_vwarn (0, "%s: No ephemeris coverage", __func__);
        return 0;
     }

   eph_beg = eph->t[0];
   eph_end = eph->t[eph->n-1];

   if ((time_end < eph_beg) || (eph_end < time_beg))
     {
        tell_vwarn (0, "%s: No ephemeris coverage", __func__);
        eph->n = 0;
        return 0;
     }

   beg = bsearch_d (time_beg, eph->t, eph->n);
   end = bsearch_d (time_end, eph->t, eph->n);

   if ((beg < 0) || (end < 0))
     {
        tell_vwarn (0, "%s: Incomplete ephemeris coverage for interval [%f,%f]",
                    __func__, time_beg, time_end);
        if (beg < 0) beg = 0;
        if (end < 0) end = eph->n-1;
     }

   if (num_pad > 0)
     {
        if (num_pad < beg)
          beg -= num_pad;
        else beg = 0;

        if (end + num_pad < (int) eph->n)
          end += num_pad;
        else end = eph->n-1;
     }

   eph->n = (end - beg + 1);

   if (beg > 0)
     {
        size_t len = eph->n * sizeof(double);
        memmove (eph->t,   eph->t   + beg, len);
        memmove (eph->r.x, eph->r.x + beg, len);
        memmove (eph->r.y, eph->r.y + beg, len);
        memmove (eph->r.z, eph->r.z + beg, len);
        memmove (eph->v.x, eph->v.x + beg, len);
        memmove (eph->v.y, eph->v.y + beg, len);
        memmove (eph->v.z, eph->v.z + beg, len);
     }

   return 0;
}

int eph_read_subset (Eph_Type *eph, const char *file,
                     double time_beg, double time_end,
                     int num_pad)
{
   wordexp_t we = {0};
   int num_times, status = -1;

   if (0 != wordexp (file, &we, WRDE_NOCMD | WRDE_UNDEF))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, file);
        return -1;
     }

   if (we.we_wordc != 1)
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: ephemeris file pattern match is non-unique: %s",
                     __func__, file);
        goto cleanup_and_return;
     }

   file = we.we_wordv[0];

   if ((num_times = count_times (file)) < 0)
     goto cleanup_and_return;

   if (0 != alloc_eph (eph, num_times))
     goto cleanup_and_return;

   if (0 != read_eph (file, num_times, eph))
     {
        eph_free (eph);
        goto cleanup_and_return;
     }

   if (0 != select_interval (eph, time_beg, time_end, num_pad))
     {
        eph_free (eph);
        goto cleanup_and_return;
     }

   status = 0;
cleanup_and_return:
   wordfree (&we);

   return status;
}
