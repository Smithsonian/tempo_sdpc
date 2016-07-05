#include "defs.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <tell.h>
#include "poly.h"

struct Polygon_Type
{
   double *v;
   int n;
   int n_alloc;
};

struct Polygon_Clip_Type
{
   Polygon_Type *p1;
};

void Polygon_free (Polygon_Type *p)
{
   if (p == NULL)
     return;
   FREE(p->v);
   FREE(p);
}

Polygon_Type *Polygon_new (int n)
{
   Polygon_Type *p = NULL;

   if (NULL == (p = (Polygon_Type *) MALLOC (sizeof *p)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   n = (n > 0) ? n : 1;

   p->n_alloc = n;
   p->n = 0;

   if (NULL == (p->v = (double *) MALLOC (2*n * sizeof(double))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        Polygon_free (p);
        return NULL;
     }

   return p;
}

static int poly_realloc (Polygon_Type *p, int n_new)
{
   double *v;

   if (n_new < 0)
     return -1;

   if (NULL == (v = (double *) REALLOC (p->v, 2*n_new * sizeof(double))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }

   p->v = v;
   p->n_alloc = n_new;

   return 0;
}

int Polygon_add (Polygon_Type *p, double x, double y)
{
   int n = p->n;

   if (n == p->n_alloc)
     {
        if (-1 == poly_realloc (p, 2*n))
          return -1;
     }

   /* don't add duplicate points */
   if ((n > 0)
       && ((fabs(p->v[2*n-2]-x) < DBL_EPSILON)
           && (fabs(p->v[2*n-1]-y) < DBL_EPSILON)))
     {
        return 0;
     }

   p->v[2*n  ] = x;
   p->v[2*n+1] = y;
   p->n++;

   return 0;
}

int Polygon_set (Polygon_Type *p, int n, const double *x, const double *y)
{
   double *p_v;
   int i;

   if (n > p->n_alloc)
     {
        if (-1 == poly_realloc (p, n))
          return -1;
     }

   p_v  = p->v;
   p->n = n;
   for (i = 0; i < n; i++)
     {
        p_v[2*i    ] = x[i];
        p_v[2*i + 1] = y[i];
     }

   return 0;
}

int Polygon_length (const Polygon_Type *p)
{
   return p->n;
}

int Polygon_vertex (const Polygon_Type *p, int i, double *x, double *y)
{
   if (i < 0 || i > p->n-1)
     return -1;
   *x = p->v[2*i  ];
   *y = p->v[2*i+1];
   return 0;
}

void Polygon_bbox (const Polygon_Type *p,
                   double *xmin, double *xmax, double *ymin, double *ymax)
{
   double *p_v, *p_vend = p->v + 2*p->n;
   double xmn, xmx, ymn, ymx;

   xmn = ymn = DBL_MAX;
   xmx = ymx = -DBL_MAX;

   for (p_v = p->v; p_v + 1 < p_vend; p_v += 2)
     {
        double x = p_v[0];
        double y = p_v[1];
        if (x < xmn) xmn = x;
        else if (x > xmx) xmx = x;
        if (y < ymn) ymn = y;
        else if (y > ymx) ymx = y;
     }

   *xmin = xmn;
   *xmax = xmx;
   *ymin = ymn;
   *ymax = ymx;
}

/* Area of planar polygon, ccw area is positive */
double Polygon_area (const Polygon_Type *p)
{
   int i, j, n = p->n;
   double *v = p->v;
   double area = 0.0;

   j = n-1;
   for (i = 0; i < n; i++)
     {
        area += (v[2*i] + v[2*j]) * (v[2*i+1] - v[2*j+1]);
        j = i;
    }

  return 0.5 * area;
}

/* returns >0 if pa,pb,pc are in CCW order */
static double orient (const double *pa, const double *pb, const double *pc)
{
   double acx, bcx, acy, bcy;
   acx = pa[0] - pc[0];
   bcx = pb[0] - pc[0];
   acy = pa[1] - pc[1];
   bcy = pb[1] - pc[1];
   return acx * bcy - acy * bcx;
}

/* returns 1 if line segments intersect,
          -1 if line segments co-linear
           0 if no intersection
 */
static int intersect (const double *a0, const double *a1,
                      const double *b0, const double *b1, double *xp)
{
   double da[2], db[2], d0[2], dbXda, g;
   da[0] = a1[0] - a0[0];
   da[1] = a1[1] - a0[1];
   db[0] = b1[0] - b0[0];
   db[1] = b1[1] - b0[1];
   d0[0] = a0[0] - b0[0];
   d0[1] = a0[1] - b0[1];

   /* a0 + f da = b0 + g db ->
    * a0 X da = b0 X da + g db X da ->
    * g = (a0 - b0) X da / (db X da)
    */
   dbXda = db[0] * da[1] - db[1] * da[0];
   if (fabs(dbXda) < DBL_EPSILON)
     return -1;  /* line segments co-linear */

   g = (d0[0] * da[1] - d0[1] * da[0]) / dbXda;
   if (g < 0.0 || 1.0 < g)
     return 0;  /* line segments don't intersect */

   /* line segments intersect at xp */
   xp[0] = b0[0] + g * db[0];
   xp[1] = b0[1] + g * db[1];

   return 1;
}

/* Sutherland-Hodgman polygon edge clipping.
 * The edge half-thickness parameter is necessary to handle the
 * case when one or more polygon edges coincide along all or
 * part of their lengths.  Its value should be "small" and positive.
 */
static int poly_edge_clip (const Polygon_Type *in,
                           double half_thickness,
                           const double *c0, const double *c1,
                           Polygon_Type *out)
{
   double *s = in->v + 2 * (in->n - 1);
   double s_loc, p_loc;
   int i;

   out->n = 0;

   s_loc = orient (c0, c1, s);
   if (s_loc > -half_thickness)
     {
        if (-1 == Polygon_add (out, s[0], s[1]))
          return -1;
     }

   for (i = 0; i < in->n; i++)
     {
        double xp[2];
        double *p = in->v + 2*i;

        p_loc = orient (c0, c1, p);
        if (((s_loc > half_thickness) && (p_loc < half_thickness))
            ||((s_loc < -half_thickness) && (p_loc > -half_thickness)))
          {
             if (intersect (c0, c1, s, p, xp) > 0)
               {
                  if (-1 == Polygon_add (out, xp[0], xp[1]))
                    return -1;
               }
          }

        if (i == in->n - 1)
          break;

        if ((p_loc > -half_thickness)
            || (fabs(s_loc) < half_thickness))
          {
             if (-1 == Polygon_add (out, p[0], p[1]))
               return -1;
          }

        s_loc = p_loc;
        s = p;
     }

   return 0;
}

/* Sutherland-Hodgman (convex) polygon clipping algorithm */
Polygon_Type *Polygon_clip (Polygon_Clip_Type *cl,
                            const Polygon_Type *in,
                            const Polygon_Type *clip)
{
   Polygon_Type *p1 = NULL;
   Polygon_Type *p2 = NULL;
   double h = DBL_EPSILON; /* half-thickness of polygon edges */
   double *cv, *cv_end;
   int status = -1;

   /* Use temp space allocated by Polygon_open_clip */
   p1 = cl->p1;

   if (NULL == (p2 = Polygon_new (in->n)))
     goto free_and_return;

   cv     = clip->v;
   cv_end = cv + 2*(clip->n-1);

   if (-1 == poly_edge_clip (in, h, cv_end, cv, p2))
     goto free_and_return;

   while (cv < cv_end)
     {
        Polygon_Type *tmp = p2;
        p2 = p1;
        p1 = tmp;

        if (p1->n < 3)
          {
             p2->n = 0;
             break;
          }

        if (-1 == poly_edge_clip (p1, h, cv, cv+2, p2))
          goto free_and_return;

        cv += 2;
     }

   status = 0;
free_and_return:

   /* preserve temp space for future calls */
   cl->p1 = p1;

   if (status)
     {
        Polygon_free (p2);
        p2 = NULL;
     }

   return p2;
}

void Polygon_close_clip (Polygon_Clip_Type *cl)
{
   if (cl == NULL)
     return;
   Polygon_free (cl->p1);
   FREE(cl);
}

Polygon_Clip_Type *Polygon_open_clip (void)
{
   Polygon_Clip_Type *cl;

   if (NULL == (cl = (Polygon_Clip_Type *)MALLOC (sizeof *cl)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if (NULL == (cl->p1 = Polygon_new (4)))
     {
        Polygon_close_clip (cl);
        return NULL;
     }

   return cl;
}
