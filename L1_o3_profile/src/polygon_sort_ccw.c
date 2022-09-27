#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

extern int polygon_sort_ccw (double *x, double *y, int n);

/* This code will only ever process polygons with 4 vertices,
 * so use static temporary arrays for the sort. */
static struct
{
   double x[4];
   double y[4];
   double angles[4];
   int index[4];
}
Angle_Sort;

static int compare_angle_indices (const void *va, const void *vb)
{
   int ia = *(int *)va;
   int ib = *(int *)vb;
   double angle_a = Angle_Sort.angles[ia];
   double angle_b = Angle_Sort.angles[ib];
   if (angle_a < angle_b) return -1;
   if (angle_a > angle_b) return +1;
   return 0;
}

/* This is a very simple-minded sort to get the pixel corner
 * points into CCW order.
 * It makes no attempt to deal with coordinate system anomalies
 * such as poles or periodicity (e.g. no dateline crossing).
 * Since we're only sorting the points, we assume there's no
 * need to worry about physical units.
 */
int polygon_sort_ccw (double *x, double *y, int n)
{
   double xc, yc;
   int i;

   if (n != 4)
     {
        fprintf (stderr, "*** USAGE ERROR: %s: expected polygon with 4 vertices!! num_vertices=%d", __func__, n);
        return -1;
     }

   xc = 0.0;
   yc = 0.0;

   for (i = 0; i < n; i++)
     {
        xc += x[i];
        yc += y[i];
     }

   xc /= n;
   yc /= n;

   /* define angles so that the sort yields points ordered NE,NW,SW,SE */
   for (i = 0; i < n; i++)
     {
        double theta = atan2 (-(y[i]-yc), -(x[i]-xc)) + M_PI;
        Angle_Sort.angles[i] = theta;
        Angle_Sort.index[i] = i;
        Angle_Sort.x[i] = x[i];
        Angle_Sort.y[i] = y[i];
     }

   qsort (Angle_Sort.index, (size_t)n, sizeof(int), compare_angle_indices);

   for (i = 0; i < n; i++)
     {
        int k = Angle_Sort.index[i];
        x[i] = Angle_Sort.x[k];
        y[i] = Angle_Sort.y[k];
        /* fprintf (stderr, "%f %f %f\n", x[i], y[i], Angle_Sort.angles[k] * 180.0/M_PI); */
     }

   return 0;
}
