#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <poly.h>

#define DEBUG 0

static int debug_print (FILE *fp, const char *fmt, ...)
{
   int status = 0;
#if DEBUG==0
   (void) fp; (void) fmt;
#else
   va_list ap;
   va_start (ap, fmt);
   status = vfprintf (fp, fmt, ap);
   va_end (ap);
#endif
   return status;
}

static int test_poly_area (void)
{
   Polygon_Type *p;
   double x[] = {0.0, 1.0, 1.0};
   double y[] = {0.0, 0.0, 1.0};
   double xx[3], yy[3];
   double area, xt, yt;
   int i;
   /* exercise realloc */
   if (NULL == (p = Polygon_new (0)))
     return -1;
   (void) Polygon_set (p, 3, x, y);
   if (0 != Polygon_get (p, 3, xx, yy))
     return -1;
   for (i = 0; i < 3; i++)
     {
        if (xx[i] != x[i])
          {
             fprintf (stderr, "***ERROR: Polygon_get:  xx[%d] = %g (expected %g)\n",
                     i, xx[i], x[i]);
             return -1;
          }
     }
   area = Polygon_area (p);
   if (fabs (area - 0.5) > DBL_EPSILON)
     {
        fprintf (stderr, "*** ERROR: Polygon_area: wrong triangle area\n");
        return -1;
     }

   /* make trangle into a unit square
    * (so test hits realloc in Polygon_area) */
   (void) Polygon_add (p, 0.0, 1.0);

   /* Test Polygon_vertex and hit error handling code */
   if ((0 != Polygon_vertex (p, 3, &xt, &yt))
       || (fabs (xt) > DBL_EPSILON) || (fabs(yt-1.0) > DBL_EPSILON))
     {
        fprintf (stderr, "*** ERROR: Polygon_vertex failed!\n");
        return -1;
     }

   if (-1 != Polygon_vertex (p, 5, &xt, &yt))
     {
        fprintf (stderr, "*** ERROR: Polygon_vertex should have failed here\n");
        return -1;
     }

   Polygon_free (p);
   return 0;
}

static void print_polygon (const Polygon_Type *p, const char *msg)
{
   int i, n = Polygon_length (p);
   if (msg) debug_print (stdout, "%s\n", msg);
   for (i = 0; i < n; i++)
     {
        double x, y;
        if (-1 == Polygon_vertex (p, i, &x, &y))
          return;
        debug_print (stdout, "%2d: (%15.7e,%15.7e)\n", i, x, y);
     }
}

static int check_area (Polygon_Type *p, double area_expected)
{
   double area;

   area = Polygon_area (p);
   if (fabs (area - area_expected) > DBL_EPSILON)
     {
        fprintf (stderr, "*** ERROR: clipped area = %g (expected %g)\n",
                 area, area_expected);
        return -1;
     }
   debug_print (stdout, "clipped area OK\n");

   return 0;
}

static int test_poly_clip (void)
{
   Polygon_Type *p=NULL, *p1=NULL, *p2=NULL, *pt=NULL;
   Polygon_Clip_Type *cl = NULL;
   /* triangle */
   double tx1[] = {0.0,  1.0, 1.0};
   double ty1[] = {0.0, -1.0, 1.0};
   /* unit square centered at origin */
   double x1[] = {-0.5,  0.5, 0.5, -0.5};
   double y1[] = {-0.5, -0.5, 0.5,  0.5};
   /* vector offsets for the clip polygon */
   double dx[] =
     {0.5, -0.5, -0.5,  0.5,   /* polygon overlap with no coincident edges */
      -0.5, 0.5,  0.0,  0.0,   /* polygon overlap with 2 coincident edges */
      0.0,                     /* all polygon edges coincident */
      2.0                      /* non-overlapping polygons */
     };
   double dy[] =
     {0.5,  0.5, -0.5, -0.5,
       0.0, 0.0, -0.5,  0.5,
     0.0,
     0.0
     };
   double area_expected[] =
     {0.25, 0.25, 0.25, 0.25,
      0.5,  0.5,  0.5,  0.5,
     1.0,
     0.0
     };
   double x2[4], y2[4];
   int j, num_cases, status = -1;

   /* exercise realloc */
   if ((NULL == (p1 = Polygon_new (0)))
       || (NULL == (p2 = Polygon_new (0))))
     goto return_status;

   (void) Polygon_set (p1, 4, x1, y1);

   if (NULL == (cl = Polygon_open_clip ()))
     goto return_status;

   num_cases = sizeof(dx)/sizeof(*dx);

   for (j = 0; j < num_cases; j++)
     {
        int k;

        debug_print (stdout, "shift clip polygon by dx=%g dy=%g\n", dx[j], dy[j]);
        for (k = 0; k < 4; k++)
          {
             x2[k] = x1[k] + dx[j];
             y2[k] = y1[k] + dy[j];
          }

        (void) Polygon_set (p2, 4, x2, y2);

        print_polygon (p1, "target polygon");
        print_polygon (p2, "clipper polygon");

        if (NULL == (p = Polygon_clip (cl, p1, p2)))
          goto return_status;

        print_polygon (p, "CLIPPED polygon");

        if (-1 == check_area (p, area_expected[j]))
          goto return_status;

        Polygon_free (p);
        p = NULL;
     }

   /* Test clipping polygon with slanted lines */
   if (NULL == (pt = Polygon_new(3)))
     goto return_status;
   (void) Polygon_set (pt, 3, tx1, ty1);
   print_polygon (p1, "target polygon");
   print_polygon (pt, "clipper triangle");
   if (NULL == (p = Polygon_clip (cl, p1, pt)))
     goto return_status;
   print_polygon (p, "CLIPPED polygon");
   if (-1 == check_area (p, 0.25))
     goto return_status;

   status = 0;
return_status:
   Polygon_close_clip (cl);
   Polygon_free (p1);
   Polygon_free (p2);
   Polygon_free (p);
   Polygon_free (pt);
   return status;
}

int main (void)
{
   if (test_poly_area ())
     return 1;
   if (test_poly_clip ())
     return 1;
   return 0;
}
