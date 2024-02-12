/* -*- mode: C; mode: fold -*- */
#include "defs.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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

/* Support diagnostic polygon output */ /*{{{*/
static int Write_Diagnostic_Polygons;

typedef struct
{
   FILE *fp;
   int wrote_first_line;
}
Diagnostic_Output_Type;

static int diagnostic_open (const char *filename, Diagnostic_Output_Type *v)
{
   char buf[1024];
   if (v == NULL)
     return 0;
   if (0 == access (filename, F_OK))
     {
        int n = 0;
        do
          {
             sprintf (buf, "%s.%d", filename, ++n);
          }
        while (0 == access (buf, F_OK));
        filename = buf;
     }
   if (NULL == (v->fp = fopen (filename, "w")))
     return -1;
   fputs ("{ \"type\" : \"GeometryCollection\", \"geometries\" : [\n", v->fp);
   v->wrote_first_line = 0;
   return 0;
}

static void diagnostic_close (Diagnostic_Output_Type *v)
{
   if ((v == NULL) || (v->fp == NULL))
     return;
   fputs ("\n]}\n", v->fp);
   (void) fclose (v->fp);
   v->fp = NULL;
}

static struct diagnostic_window
{
   double x0, y0;
   double dx, dy;
}
Diagnostic_Window =
{
   .x0 = 0.0,
   .y0 = 0.0,
   .dx = DBL_MAX,
   .dy = DBL_MAX,
};

void __Pixel_diagnostic_output (int i)
{
   Write_Diagnostic_Polygons = i;
}

void __Pixel_diagnostic_window (double x, double y, double dx, double dy)
{
   struct diagnostic_window *w = &Diagnostic_Window;
   w->x0 = x;
   w->y0 = y;
   w->dx = dx;
   w->dy = dy;
}

static void diagnostic_poly_write (Diagnostic_Output_Type *v, Polygon_Type *p, int index)
{
   FILE *fp = v->fp;
   struct diagnostic_window  *w = &Diagnostic_Window;
   const char *prefix = "{\"type\": \"Polygon\", \"coordinates\": [";
   double x, y;
   int i, n;

   if ((fp == NULL) || (p == NULL))
     return;

   n = Polygon_length (p);

   /* write out polygons with any vertex inside the window */
   for (i = 0; i < n; i++)
     {
        if (0 != Polygon_vertex (p, i, &x, &y))
          return;
        if ((fabs (x - w->x0) < w->dx)
            && (fabs (y - w->y0) < w->dy))
          break;
     }
   /* exclude polygons entirely outside the window */
   if (i == n)
     return;

   if (0 != Polygon_vertex (p, 0, &x, &y))
     return;

   if (v->wrote_first_line) fputs (",\n", fp);
   fprintf (fp, "%s[", prefix);
   fprintf (fp, "[%12.8e,%12.8e]", x, y);
   for (i = 1; i < n; i++)
     {
        if (0 != Polygon_vertex (p, i, &x, &y))
          break;
        fprintf (fp, ",[%12.8e,%12.8e]", x, y);
     }
   fprintf (fp, "]], \"id\": %d}", index);
   v->wrote_first_line = 1;
}

/*}}}*/

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

#if 0
int __Pixel_print_overlap (const Pixel_Regrid_Type *r, int dest_pixel_index)
{
   Pixel_Overlap_Type *o;
   int j;

   if (r->overlap == NULL)
     return -1;

   o = r->overlap[dest_pixel_index];

   fprintf (stderr, "src_index  overlap_area\n");
   for (j = 0; j < o->num_overlaps; j++)
     {
        fprintf (stderr, "%6d  %15.8e\n", o->src_index[j], o->area[j]);
     }
   return 0;
}
#endif

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

#define IS_FINITE(x) (((x)==(x)) && (fabs(x)<DBL_MAX))
#define VALID_POINT(x,y) (IS_FINITE(x) && IS_FINITE(y))

int Pixel_list_pack (Pixel_List_Type *pixel_list,
                     double *xs, double *ys, int num_pixels,
                     int num_pixel_vertices,
                     int *step, int *xtrack, int num_xtrack,
                     int xtrack_dimlen)
{
   int i;

   for (i = 0; i < num_pixels; i++)
     {
        int pix_xtrack = xtrack[i % num_xtrack];
        int pix_step   =   step[i / num_xtrack];
        /* pixel index in target full-scan array */
        int pix = pix_xtrack + pix_step * xtrack_dimlen;
        double *x = xs + i * num_pixel_vertices;
        double *y = ys + i * num_pixel_vertices;
        int j;

        for (j = 0; j < num_pixel_vertices; j++)
          {
             if (0 == VALID_POINT(x[j],y[j]))
               break;
          }
        if (j == num_pixel_vertices)
          {
             if ((-1 == Pixel_list_set_vertices (pixel_list, i, num_pixel_vertices, x, y))
                 || (-1 == Pixel_list_set_src_index (pixel_list, i, pix)))
               return -1;
          }
     }

   return 0;
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
   int num_pixel_vertices =
     (4 + 2 * g->num_extra_xpoints + 2 * g->num_extra_ypoints);
   int num_vertices = num_pixels * num_pixel_vertices;
   int nx = g->nx;
   double *xs=NULL, *ys=NULL;
   double *x, *y;
   double dx_j, dy_j;
   int k, status = -1;

   if ((NULL == (xs = (double *) MALLOC (num_vertices * sizeof(double))))
       || (NULL == (ys = (double *) MALLOC (num_vertices * sizeof(double)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   dx_j = dx / (1 + g->num_extra_xpoints);
   dy_j = dy / (1 + g->num_extra_ypoints);

   for (k = 0; k < num_pixels; k++)
     {
        double x0, y0, x1, y1;
        int i, j, ix, iy;

        ix = k % nx;
        iy = k / nx;

        x0 = xmin + ix * dx;
        y0 = ymin + iy * dy;
        x1 = x0 + dx;
        y1 = y0 + dy;

        x = xs + k * num_pixel_vertices;
        y = ys + k * num_pixel_vertices;

        j = 0;

        for (i = 0; i <= g->num_extra_xpoints; i++, j++)
          {
             x[j] = x0 + i * dx_j;
             y[j] = y0;
          }
        for (i = 0; i <= g->num_extra_ypoints; i++, j++)
          {
             x[j] = x1;
             y[j] = y0 + i * dy_j;
          }
        for (i = 0; i <= g->num_extra_xpoints; i++, j++)
          {
             x[j] = x1 - i * dx_j;
             y[j] = y1;
          }
        for (i = 0; i <= g->num_extra_ypoints; i++, j++)
          {
             x[j] = x0;
             y[j] = y1 - i * dy_j;
          }
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
   int num_negative_overlap_areas = 0;
   /* Support for diagnostic polygon output */
   Diagnostic_Output_Type v_src = {0};
   Diagnostic_Output_Type v_overlap = {0};

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

   if (Write_Diagnostic_Polygons != 0)
     {
        diagnostic_open ("polygons_overlap.json", &v_overlap);
        diagnostic_open ("polygons_src.json", &v_src);
        if (have_dest_polygons)
          {
             Diagnostic_Output_Type v_dest = {0};
             diagnostic_open ("polygons_dest.json", &v_dest);
             for (k = 0; k < dest_area->num_polys; k++)
               {
                  diagnostic_poly_write (&v_dest, dest_area->poly[k], k);
               }
             diagnostic_close (&v_dest);
          }
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
        int i, j, imn, imx, jmn, jmx, num_invalid;

        num_invalid = Polygon_bbox (src_poly_lookup, &xmn, &xmx, &ymn, &ymx);
        /* Does source polygon bbox have any invalid vertex coordinates ? */
        if (num_invalid) continue;

        /* Does source polygon bbox overlap destination grid anywhere? */
        if ((xmn > dest->xmax) || (xmx < dest->xmin)
            || (ymn > dest->ymax) || (ymx < dest->ymin))
          {
             continue;
          }

        if (Write_Diagnostic_Polygons != 0)
          {
             diagnostic_poly_write (&v_src, src_poly_area, k);
          }

        num_src_dest_overlap++;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

        /* clip bounding box */
        xmn = MAX(xmn, dest->xmin);
        xmx = MIN(xmx, dest->xmax);
        ymn = MAX(ymn, dest->ymin);
        ymx = MIN(ymx, dest->ymax);

        /* find destination cell index range. */
        imn = (int) floor((xmn - dest->xmin) / dx);
        imx = (int)  ceil((xmx - dest->xmin) / dx);
        jmn = (int) floor((ymn - dest->ymin) / dy);
        jmx = (int)  ceil((ymx - dest->ymin) / dy);

        for (j = jmn; j < jmx; j++)
          {
             for (i = imn; i < imx; i++)
               {
                  Polygon_Type *p = NULL;
                  int dest_poly_index = i + j * dest->nx;
                  double overlap_area, xp0, yp0;

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
                  if (overlap_area < 0.0)
                    {
                       num_negative_overlap_areas++;
                       if (0 == Polygon_vertex (p, 0, &xp0, &yp0))
                         {
                            tell_vwarn (1, "(x,y)=(%0.1f, %0.1f) overlap_area=%10.3e km^2 ",
                                        xp0, yp0, overlap_area/1.e6);
                         }
                    }
                  if ((Write_Diagnostic_Polygons != 0) && (overlap_area > 0.0))
                    {
                       diagnostic_poly_write (&v_overlap, p, dest_poly_index);
                    }
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

   if (num_negative_overlap_areas > 0)
     {
        tell_vwarn (0, "found %d negative pixel overlap areas", num_negative_overlap_areas);
     }

   Polygon_close_clip (cl);
   if (have_dest_polygons == 0)
     {
        Polygon_free (dest_poly);
     }

   if (Write_Diagnostic_Polygons)
     {
        diagnostic_close (&v_overlap);
        diagnostic_close (&v_src);
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
                                int src_num_step, int src_num_xtrack)
{
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
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

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
        rs->num[i] = PIXEL_INIT_NUM_SAMPLES;
     }

   return rs;
}

/* How many source pixels map to each destination mesh pixel?
 * This is primarily for diagnostic purposes */
int *Pixel_regrid_overlap_map (const Pixel_Regrid_Type *r)
{
   int *map = NULL;
   int i, j, num_src_pixels;

   if (r->overlap == NULL)
     return NULL;

   num_src_pixels = r->num_src_step * r->num_src_xtrack;
   if (NULL == (map = (int *)MALLOC (num_src_pixels * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)map, 0, num_src_pixels * sizeof(int));

   for (i = 0; i < r->num_dest_pixels; i++)
     {
        Pixel_Overlap_Type *o = r->overlap[i];
        if (o == NULL)
          continue;
        for (j = 0; j < o->num_overlaps; j++)
          {
             map[o->src_index[j]] += 1;
          }
     }

   return map;
}

double *Pixel_regrid_area_weight_sum (const Pixel_Regrid_Type *r,
                                      double scale, double fill_value)
{
   double *a_sum = NULL;
   double a_sum_i;
   int i, j;

   if (r->overlap == NULL)
     return NULL;

   if (scale <= 0.0) scale = 1.0;

   if (NULL == (a_sum = (double *)MALLOC (r->num_dest_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   for (i = 0; i < r->num_dest_pixels; i++)
     {
        if (r->overlap[i])
          {
             Pixel_Overlap_Type *o = r->overlap[i];
             a_sum_i = 0.0;
             for (j = 0; j < o->num_overlaps; j++)
               {
                  a_sum_i += o->area[j];
               }
             a_sum[i] = a_sum_i * scale;
          }
        else a_sum[i] = fill_value;
     }

   return a_sum;
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
             rs->min[i] = PIXEL_INIT_MIN_SAMPLE;
             rs->max[i] = PIXEL_INIT_MAX_SAMPLE;
             rs->num[i] = PIXEL_INIT_NUM_SAMPLES;
          }

        if (o == NULL)
          continue;

        a_sum = awt_sum = 0.0;

        for (j = 0; j < o->num_overlaps; j++)
          {
             int k = o->src_index[j];
             if ((src_mask == NULL)
                 || ((src_mask != NULL) && (src_mask[k] == 0)))
               {
                  double a = o->area[j];
                  double src_k = src[k];
                  if ((src_k == fill_value) || (0 == isfinite(src_k))) continue;
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

int Pixel_regrid_dest_boundary (const Pixel_Regrid_Type *r, const Pixel_Grid_Param_Type *dest,
                                int **pindices, int *num)
{
   int *mask=NULL, *indices=NULL;
   int *byb=NULL, *byt=NULL, *bdry=NULL;
   int max_num_boundary, y_first_ok, y_last_ok;
   int min_dy = (dest->ny > 100) ? (dest->ny/100) : 5;
   int i, j, n, status = -1;

   *pindices = NULL;
   *num = 0;

   /* Quick return if source and destination grids don't overlap. */
   if (r->overlap == NULL)
     return 0;

   max_num_boundary = 2 * (dest->nx + dest->ny);

   if ((NULL == (mask = (int *)MALLOC (r->num_dest_pixels * sizeof(int))))
       || (NULL == (bdry = (int *)MALLOC (max_num_boundary * sizeof(int))))
       || (NULL == (indices = (int *)MALLOC (max_num_boundary * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   /* create a mask showing which destination pixels
    * receive source contributions */
   for (i = 0; i < r->num_dest_pixels; i++)
     {
        mask[i] = (r->overlap[i] != NULL);
     }

   byt = indices;
   byb = byt + max_num_boundary/2;

   /* for each x, find the y boundaries */
   for (i = 0; i < dest->nx; i++)
     {
        y_first_ok = -1;
        y_last_ok = -1;
        for (j = 0; j < dest->ny; j++)
          {
             if (mask[i + j * dest->nx] != 0)
               {
                  y_last_ok = j;
                  if (y_first_ok < 0)
                    y_first_ok = j;
               }
          }
        /* require segment size to exceed some threshold */
        if ((y_first_ok >= 0) && ((y_last_ok-y_first_ok) > min_dy))
          {
             byb[i] = i + y_first_ok * dest->nx;
             byt[i] = i + y_last_ok * dest->nx;
          }
        else
          {
             byb[i] = -1;
             byt[i] = -1;
          }
     }

   /* The region boundary polygon we want is now this set of points
    * {byb, rev(byt)}, skipping -1 boundary indices.
    */
   n = 0;
   for (i = 0; i < dest->nx; i++)
     {
        if (byb[i] >= 0)
          {
             bdry[n] = byb[i];
             n++;
          }
     }
   for (i = dest->nx-1; i >= 0; i--)
     {
        if (byt[i] >= 0)
          {
             bdry[n] = byt[i];
             n++;
          }
     }

   *pindices = bdry;
   *num = n;

   status = 0;
return_status:
   FREE(mask);
   FREE(indices);
   return status;
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
             if ((src_mask == NULL) \
                 || ((src_mask != NULL) && (src_mask[k] == 0))) \
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

/* Regrid in the "reverse" direction:
 * FROM rectangular, non-overlapping mesh TO pixel list.
 * The overlap data structure is defined in the same way as before,
 * and we loop over the mesh pixels in the same order, but the
 * value assignment is in the opposite direction.
 *
 *  1) Pixel_Regrid_Type contains overlap info on all granules
 *  2) mesh_values covers the complete field of interest
 *     (e.g. enclosing all granules) for the mesh grid, whatever
 *     it may be, e.g. longitude-latitude.
 *     mesh_mask has the same dimensions as mesh_values.
 *  3) values covers all granules for the (mirror_step, xtrack)
 *     coordinate system.
 */
int Pixel_regrid_from_mesh (const Pixel_Regrid_Type *r, const int *mesh_mask,
                            double fill_value, const double *mesh_values,
                            double *values)
{
   int mesh_pixel, num_mesh_pixels;
   int list_pixel, num_list_pixels;
   double *wrk = NULL;
   double *a_sum = NULL;
   double *awt_sum = NULL;

   /* Quick return if source and destination grids don't overlap. */
   if (r->overlap == NULL)
     return 0;

   /* full field-of-regard dimensions */
   num_list_pixels = r->num_src_step * r->num_src_xtrack;
   num_mesh_pixels = r->num_dest_pixels;

   if (NULL == (wrk = (double *) MALLOC (2*num_list_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   awt_sum = wrk;
   a_sum = wrk + num_list_pixels;

   for (list_pixel = 0; list_pixel < num_list_pixels; list_pixel++)
     {
        values[list_pixel] = fill_value;
        awt_sum[list_pixel] = 0.0;
        a_sum[list_pixel] = 0.0;
     }

   for (mesh_pixel = 0; mesh_pixel < num_mesh_pixels; mesh_pixel++)
     {
        Pixel_Overlap_Type *o = r->overlap[mesh_pixel];
        double val_m;
        int j;

        if ((o == NULL)
            || ((mesh_mask != NULL) && (mesh_mask[mesh_pixel] != 0)))
          continue;

        val_m = mesh_values[mesh_pixel];
        if ((val_m == fill_value) || (val_m != val_m))
          continue;

        for (j = 0; j < o->num_overlaps; j++)
          {
             double a = o->area[j];
             list_pixel = o->src_index[j];
             awt_sum[list_pixel] += a * val_m;
             a_sum[list_pixel]   += a;
          }
     }

   for (list_pixel = 0; list_pixel < num_list_pixels; list_pixel++)
     {
        if (a_sum[list_pixel] != 0.0)
          {
             values[list_pixel] = awt_sum[list_pixel] / a_sum[list_pixel];
          }
     }

   FREE(wrk);

   return 0;
}
