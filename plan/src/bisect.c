/** @file bisect.c
 *  @brief Bisection method root solver.
 */

#include "config.h"
#include <math.h>
#include <stdio.h>

#include "bisect.h"

int bisection (int (*func)(double, double *, void *),
               double a, double b, void *cd, double *xp)
{
   unsigned int count = 1;
   unsigned int bisect_count = 5;
   double fa, fb;
   double x = a;

   if (a > b)
     {
        double tmp = a; a = b; b = tmp;
     }

   a = nextafter (a, b);
   b = nextafter (b, a);

   if (0 != (*func)(a, &fa, cd))
     return -1;
   if (0 != (*func)(b, &fb, cd))
     return -1;

   if (fa * fb > 0)
     {
        /* "bisection: the interval may not bracket a root: f(a)*f(b)>0" */
        *xp = a;
        return -1;
     }

   while (b > a)
     {
        double fx;

        if (fb == 0.0)
          {
             x = b;
             break;
          }
        if (fa == 0.0)
          {
             x = a;
             break;
          }

        if (count % bisect_count)
          {
             x = (a*fb - b*fa) / (fb - fa);
             if (x < a || b < x)
               x = 0.5 * (a + b);
          }
        else x = 0.5 * (a + b);

        x = nextafter (x, a);

        if (x <= a || b <= x)
          break;

        if (0 != (*func)(x, &fx, cd))
          return -1;
        count++;

        if (fx*fa < 0.0)
          {
             fb = fx;
             b = x;
          }
        else
          {
             fa = fx;
             a = x;
          }
     }

   *xp = x;
   return 0;
}

#ifdef TESTING

#include <float.h>
#include <stdio.h>

static int f (double x, double *pf, void *cd)
{
   (void) cd;

   *pf = x - cos(x);

   return 0;
}

int main (void)
{
   double x, fx;

   (void) bisection (&f, 0.0, 1.0, NULL, &x);

   (void) f(x,&fx,NULL);

   if (fabs(fx) > DBL_EPSILON)
     fprintf (stdout, "Error\n");

   fprintf (stdout, "x = %15.13f\n", x);

   return 0;
}

#endif
