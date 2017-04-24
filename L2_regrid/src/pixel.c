#include "defs.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tell.h>
#include "poly.h"
#include "pixel.h"

#define NUM_OVERLAPS_HINT  4

typedef struct
{
   double *area;
   int *src_index;
   int num_overlaps;
   int num_alloc;
}
Pixel_Overlap_Type;

struct Pixel_List_Type
{
   Polygon_Type **poly;
   int *src_index;      /* NULL if no indirection needed */
   int num_polys;
};

struct Pixel_Regrid_Type
{
   Pixel_Overlap_Type **overlap;
   const Pixel_Grid_Param_Type *dest;
   const Pixel_List_Type *dest_area;
   int num_dest_pixels;
   int num_src_step;
   int num_src_xtrack;
};

static void free_overlap (Pixel_Overlap_Type *o)
{
   if (NULL == o)
     return;
   FREE(o->area);
   FREE(o->src_index);
   FREE(o);
}

static Pixel_Overlap_Type *new_overlap (void)
{
   Pixel_Overlap_Type *o = NULL;
   int n = NUM_OVERLAPS_HINT;

   if (NULL == (o = (Pixel_Overlap_Type *)MALLOC (sizeof *o)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if ((NULL == (o->area = (double *)MALLOC (n * sizeof(double))))
       || (NULL == (o->src_index = (int *)MALLOC (n * sizeof(int)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_overlap (o);
        return NULL;
     }

   o->num_overlaps = 0;
   o->num_alloc = n;

   return o;
}

static int realloc_overlap (Pixel_Overlap_Type *o, int new_n)
{
   double *area = NULL;
   int *src_index = NULL;

   if ((NULL == (area = (double *) REALLOC (o->area, new_n * sizeof(double))))
       || (NULL == (src_index = (int *) REALLOC (o->src_index, new_n * sizeof(int)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        FREE(area);
        return -1;
     }

   o->area = area;
   o->src_index = src_index;
   o->num_alloc = new_n;

   return 0;
}

static int push_overlap (Pixel_Overlap_Type *o, double area, int src_index)
{
   int n = o->num_overlaps;

   if (n == o->num_alloc)
     {
        if (-1 == realloc_overlap (o, 2*o->num_alloc))
          return -1;
     }

   o->area[n] = area;
   o->src_index[n] = src_index;
   o->num_overlaps++;

   return 0;
}

void Pixel_list_free (Pixel_List_Type *lst)
{
   if (lst == NULL)
     return;

   if (lst->poly != NULL)
     {
        int i;
        for (i = 0; i < lst->num_polys; i++)
          {
             Polygon_free (lst->poly[i]);
          }
        FREE(lst->poly);
     }
   FREE(lst->src_index);
   FREE(lst);
}

Pixel_List_Type *Pixel_list_new (int num_polys, int num_sides)
{
   Pixel_List_Type *lst = NULL;
   int i;

   if (NULL == (lst = (Pixel_List_Type *) MALLOC (sizeof (*lst))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   lst->poly = (Polygon_Type **) MALLOC (num_polys * sizeof (Polygon_Type *));
   if (lst->poly == NULL)
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   lst->num_polys = num_polys;
   for (i = 0; i < num_polys; i++)
     {
        if (NULL == (lst->poly[i] = Polygon_new (num_sides)))
          goto free_and_return;
     }

   lst->src_index = NULL;

   return lst;
free_and_return:
   Pixel_list_free (lst);
   return NULL;
}

int Pixel_list_use_src_index (Pixel_List_Type *lst)
{
   lst->src_index = (int *) MALLOC (lst->num_polys * sizeof(int));
   return lst->src_index ? 0 : -1;
}

int Pixel_list_set_src_index (Pixel_List_Type *lst, int i, int src_index)
{
   if (lst->src_index == NULL)
     return -1;
   lst->src_index[i] = src_index;
   return 0;
}

int Pixel_list_set_vertices (Pixel_List_Type *lst, int pix, int n,
                             const double *x, const double *y)
{
   return Polygon_set (lst->poly[pix], n, x, y);
}

int Pixel_grid_arrays (const Pixel_Grid_Param_Type *g,
                       double **x_corners, double **y_corners)
{
   double xsize = g->xmax - g->xmin;
   double ysize = g->ymax - g->ymin;
   double dx = xsize / g->nx;
   double dy = ysize / g->ny;
   double xmin = g->xmin;
   double ymin = g->ymin;
   int num_pixels = g->nx * g->ny;
   int nx = g->nx;
   double *xs=NULL, *ys=NULL;
   double *x, *y;
   int k, status = -1;

   if ((NULL == (xs = (double *) MALLOC (4*num_pixels * sizeof(double))))
       || (NULL == (ys = (double *) MALLOC (4*num_pixels * sizeof(double)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   for (k = 0; k < num_pixels; k++)
     {
        int ix = k % nx;
        int iy = k / nx;

        x = xs + 4*k;
        y = ys + 4*k;
        x[0] = xmin + ix * dx;
        y[0] = ymin + iy * dy;
        x[1] = x[0] + dx;
        y[1] = y[0];
        x[2] = x[1];
        y[2] = y[1] + dy;
        x[3] = x[0];
        y[3] = y[2];
     }

   *x_corners = xs;
   *y_corners = ys;

   status = 0;
free_and_return:
   if (status)
     {
        FREE(xs);
        FREE(ys);
     }

   return status;
}

static int poly_fill (const Pixel_Grid_Param_Type *dest,
                      int ix, int iy, double dx, double dy,
                      Polygon_Type *p)
{
   double x[4], y[4];

   x[0] = dest->xmin + ix * dx;
   x[1] = x[0] + dx;
   x[2] = x[1];
   x[3] = x[0];

   y[0] = dest->ymin + iy * dy;
   y[1] = y[0];
   y[2] = y[0] + dy;
   y[3] = y[2];

   return Polygon_set (p, 4, x, y);
}

static double my_round (double x)
{
   double xf, xi;

   xf = modf (x, &xi);                 /* x = xi + xf */
   if (xi > 0)
     {
        if (xf >= 0.5)
          return xi + 1.0;
     }
   else if (xi < 0)
     {
        if (xf <= -0.5)
          return xi - 1.0;
     }
   else if (xf < 0)                    /* xi=0 */
     {
        if (xf <= -0.5)
          return -1.0;
     }
   else if (xf >= 0.5)                 /* xi=0 */
     return 1.0;

   return xi;
}
/* C99 provides round, but for now, avoid requiring C99 */
#define ROUND(x) my_round(x)

/* assume array overlap[dest->num_polys] already allocated */
int Pixel_find_overlaps (Pixel_Regrid_Type *r,
                         const Pixel_List_Type *src_area,
                         const Pixel_List_Type *src_lookup)
{
   const Pixel_Grid_Param_Type *dest = r->dest;
   const Pixel_List_Type *dest_area = r->dest_area;
   Pixel_Overlap_Type **overlap = r->overlap;
   Polygon_Type *dest_poly = NULL;
   Polygon_Clip_Type *cl = NULL;
   double xsize, ysize, dx, dy;
   int have_dest_polygons;
   int k, num_src_dest_overlap = 0;

   xsize = dest->xmax - dest->xmin;
   ysize = dest->ymax - dest->ymin;

   dx = xsize / dest->nx;
   dy = ysize / dest->ny;

   have_dest_polygons = (dest_area != NULL);

   if (have_dest_polygons == 0)
     {
        if (NULL == (dest_poly = Polygon_new (4)))
          return -1;
     }

   if (NULL == (cl = Polygon_open_clip ()))
     return -1;

   if (src_lookup == NULL)
     src_lookup = src_area;

   /* partition each source polygon
    * among overlapping destination polygons
    */
   for (k = 0; k < src_area->num_polys; k++)
     {
        Polygon_Type *src_poly_area = src_area->poly[k];
        Polygon_Type *src_poly_lookup = src_lookup->poly[k];
        double xmn, xmx, ymn, ymx;
        int i, j, imn, imx, jmn, jmx;

        Polygon_bbox (src_poly_lookup, &xmn, &xmx, &ymn, &ymx);

        /* Does source polygon bbox overlap destination grid anywhere? */
        if ((xmn > dest->xmax) || (xmx < dest->xmin)
            || (ymn > dest->ymax) || (ymx < dest->ymin))
          {
             continue;
          }

        num_src_dest_overlap++;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

        /* clip bounding box */
        xmn = MAX(xmn, dest->xmin);
        xmx = MIN(xmx, dest->xmax);
        ymn = MAX(ymn, dest->ymin);
        ymx = MIN(ymx, dest->ymax);

        /* find destination cell index range */
        imn = (int) floor((xmn - dest->xmin) / dx);
        imx = (int) ROUND((xmx - dest->xmin) / dx);
        jmn = (int) floor((ymn - dest->ymin) / dy);
        jmx = (int) ROUND((ymx - dest->ymin) / dy);

        for (j = jmn; j < jmx; j++)
          {
             for (i = imn; i < imx; i++)
               {
                  Polygon_Type *p = NULL;
                  int dest_poly_index = i + j * dest->nx;
                  double overlap_area;

                  if (have_dest_polygons)
                    {
                       dest_poly = dest_area->poly[dest_poly_index];
                    }
                  else
                    {
                       if (-1 == poly_fill (dest, i, j, dx, dy, dest_poly))
                         return -1;
                    }

                  if (NULL == (p = Polygon_clip (cl, src_poly_area, dest_poly)))
                    return -1;
                  overlap_area = Polygon_area (p);
                  Polygon_free (p);

                  if (overlap_area > 0.0)
                    {
                       int src_index;
                       Pixel_Overlap_Type *o = overlap[dest_poly_index];
                       if (NULL == o)
                         {
                            if (NULL == (o = new_overlap()))
                              return -1;
                            overlap[dest_poly_index] = o;
                         }
                       src_index = ((src_area->src_index == NULL) ? k :
                                    src_area->src_index[k]);
                       if (-1 == push_overlap (o, overlap_area, src_index))
                         return -1;
                    }
               }
          }
     }

   Polygon_close_clip (cl);
   if (have_dest_polygons == 0)
     {
        Polygon_free (dest_poly);
     }

   return num_src_dest_overlap;
}

void Pixel_close_regrid (Pixel_Regrid_Type *r)
{
   if (r == NULL)
     return;

   if (r->overlap)
     {
        int n = r->num_dest_pixels;
        while (n-- > 0)
          {
             free_overlap (r->overlap[n]);
          }
        FREE(r->overlap);
     }

   FREE(r);
}

void Pixel_regrid_grow_srcdims (Pixel_Regrid_Type *r,
                                int src_max_step, int src_max_xtrack)
{
   int src_num_step = src_max_step + 1;
   int src_num_xtrack = src_max_xtrack + 1;
   if (src_num_step > r->num_src_step)
     r->num_src_step = src_num_step;
   if (src_num_xtrack > r->num_src_xtrack)
     r->num_src_xtrack = src_num_xtrack;
}

void Pixel_regrid_get_srcdims (const Pixel_Regrid_Type *r,
                               int *num_src_step, int *num_src_xtrack)
{
   *num_src_step = r->num_src_step;
   *num_src_xtrack = r->num_src_xtrack;
}

Pixel_Regrid_Type *
Pixel_open_regrid (const Pixel_Grid_Param_Type *dest,
                   const Pixel_List_Type *dest_area)
{
   Pixel_Regrid_Type *r = NULL;
   int len_overlap, num_dest;

   /* dest_area == NULL is ok */

   if (NULL == (r = (Pixel_Regrid_Type *) MALLOC (sizeof *r)))
     return NULL;

   num_dest = dest->nx * dest->ny;
   len_overlap = num_dest * sizeof (Pixel_Overlap_Type *);
   if (NULL == (r->overlap = (Pixel_Overlap_Type **) MALLOC (len_overlap)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        Pixel_close_regrid (r);
        return NULL;
     }

   memset ((char *)r->overlap, 0, len_overlap);

   r->dest = dest;
   r->dest_area = dest_area;
   r->num_dest_pixels = num_dest;

   r->num_src_step = 0;
   r->num_src_xtrack = 0;

   return r;
}

void Pixel_free_regrid_stats (Pixel_Regrid_Stats_Type *rs)
{
   if (rs == NULL)
     return;
   FREE(rs->min);
   FREE(rs->max);
   FREE(rs->num);
   FREE(rs);
}

Pixel_Regrid_Stats_Type *
Pixel_alloc_regrid_stats (int num_pixels, int num_values_per_pixel)
{
   Pixel_Regrid_Stats_Type *rs = NULL;
   int i, len;

   if (NULL == (rs = (Pixel_Regrid_Stats_Type *) MALLOC (sizeof *rs)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)rs, 0, sizeof *rs);

   len = num_pixels * num_values_per_pixel;
   if ((NULL == (rs->min = (double *)MALLOC (len * sizeof(double))))
       || (NULL == (rs->max = (double *)MALLOC (len * sizeof(double))))
       || (NULL == (rs->num = (int *)MALLOC (len * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        Pixel_free_regrid_stats (rs);
        return NULL;
     }

   for (i = 0; i < len; i++)
     {
        rs->min[i] = PIXEL_INIT_MIN_SAMPLE;
        rs->max[i] = PIXEL_INIT_MAX_SAMPLE;
        rs->num[i] = 0;
     }

   return rs;
}

int Pixel_regrid (const Pixel_Regrid_Type *r, const int *src_mask,
                  double fill_value, const double *src, double *dest,
                  Pixel_Regrid_Stats_Type *rs)
{
   int i, j;

   /* Quick return if source and destination grids don't overlap. */
   if (r->overlap == NULL)
     return 0;

   for (i = 0; i < r->num_dest_pixels; i++)
     {
        Pixel_Overlap_Type *o = r->overlap[i];
        double a_sum, awt_sum;
        double mn = PIXEL_INIT_MIN_SAMPLE;
        double mx = PIXEL_INIT_MAX_SAMPLE;
        int num = 0;

        dest[i] = fill_value;

        if (rs)
          {
             rs->min[i] = mn;
             rs->max[i] = mx;
             rs->num[i] = num;
          }

        if (o == NULL)
          continue;

        a_sum = awt_sum = 0.0;

        for (j = 0; j < o->num_overlaps; j++)
          {
             int k = o->src_index[j];
             if (src_mask[k] == 0)
               {
                  double a = o->area[j];
                  double src_k = src[k];
                  awt_sum += a * src_k;
                  a_sum   += a;

                  /* keep it simple -- always gather stats */
                  if (src_k > mx) mx = src_k;
                  if (src_k < mn) mn = src_k;
                  num += 1;
               }
          }

        if (a_sum > 0.0)
          {
             dest[i] = awt_sum / a_sum;
          }

        if (rs)
          {
             rs->min[i] = mn;
             rs->max[i] = mx;
             rs->num[i] = num;
          }
     }

   return 0;
}

#define REGRID_BYTES(typestr, type) \
static int regrid_bytes_##typestr (const Pixel_Regrid_Type *r, const int *src_mask, \
                                   const type *fill_value, const type *src, type *dest) \
{ \
   type or_all; \
   int i; \
 \
   /* Quick return if source and destination grids don't overlap. */ \
   if (r->overlap == NULL) \
     return 0; \
 \
   for (i = 0; i < r->num_dest_pixels; i++) \
     { \
        Pixel_Overlap_Type *o = r->overlap[i]; \
        int j; \
 \
        dest[i] = *fill_value; \
 \
        if (o == NULL) \
          continue; \
 \
        or_all = (type) 0; \
 \
        for (j = 0; j < o->num_overlaps; j++) \
          { \
             int k = o->src_index[j]; \
             if (src_mask[k] == 0) \
               { \
                  or_all |= src[k]; \
               } \
          } \
        dest[i] = or_all; \
     } \
 \
   return 0; \
}

REGRID_BYTES(uint64, unsigned long long)
REGRID_BYTES(uint,   unsigned int)
REGRID_BYTES(ushort, unsigned short)
REGRID_BYTES(ubyte,  unsigned char)
REGRID_BYTES(int64, long long)
REGRID_BYTES(int,   int)
REGRID_BYTES(short, short)
REGRID_BYTES(byte,  char)

int Pixel_regrid_bytes (const Pixel_Regrid_Type *r, const int *src_mask,
                        int value_type, const void *fill_value,
                        const void *src, void *dest)
{
   switch (value_type)
     {
      case VALUE_IS_UINT64:
        return regrid_bytes_uint64 (r, src_mask, fill_value, src, dest);
      case VALUE_IS_UINT:
        return regrid_bytes_uint   (r, src_mask, fill_value, src, dest);
      case VALUE_IS_USHORT:
        return regrid_bytes_ushort (r, src_mask, fill_value, src, dest);
      case VALUE_IS_UBYTE:
        return regrid_bytes_ubyte  (r, src_mask, fill_value, src, dest);
      case VALUE_IS_INT64:
        return regrid_bytes_int64 (r, src_mask, fill_value, src, dest);
      case VALUE_IS_INT:
        return regrid_bytes_int   (r, src_mask, fill_value, src, dest);
      case VALUE_IS_SHORT:
        return regrid_bytes_short (r, src_mask, fill_value, src, dest);
      case VALUE_IS_BYTE:
        return regrid_bytes_byte  (r, src_mask, fill_value, src, dest);
     }

   return -1;
}
