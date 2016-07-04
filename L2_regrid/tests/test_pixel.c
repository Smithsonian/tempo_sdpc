#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tell.h>
#include <poly.h>
#include <pixel.h>

#ifndef MALLOC
# define MALLOC malloc
#endif

#ifndef FREE
# define FREE free
#endif

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

static int test_regrid (int nx_src, int ny_src, float bin_factor,
                        double pixel_xoverlap,
                        double xshift, double yshift,
                        int expect_no_overlaps)
{
   Pixel_Grid_Param_Type dest_grid_params;
   Pixel_List_Type *src_pixel_list = NULL;
   Pixel_Regrid_Type *r = NULL;
   Pixel_Overlap_Info_Type *info = NULL;
   double *src_values=NULL, *dest_values=NULL;
   double x0 = 0.0, dx = 1.0/nx_src;
   double y0 = 0.0, dy = 1.0/ny_src;
   double src_sum, dest_sum, expected_dest_sum;
   double *xcnr = NULL, *ycnr = NULL;
   int num_src, nx_dest, ny_dest, num_dest, num_overlaps;
   int i, *src_mask= NULL;
   int status = -1;

   num_src = nx_src * ny_src;
   nx_dest = nx_src / bin_factor;
   ny_dest = ny_src / bin_factor;
   num_dest = nx_dest * ny_dest;

   if ((NULL == (src_mask = (int *)MALLOC (num_src*sizeof(int))))
       || (NULL == (src_values = (double *)MALLOC (num_src*sizeof(double))))
       || (NULL == (dest_values = (double *)MALLOC (num_dest*sizeof(double)))))
     {
        goto return_status;
     }

   if (NULL == (src_pixel_list = Pixel_list_new (num_src, 4)))
     goto return_status;

   if (-1 == Pixel_list_use_src_index (src_pixel_list))
     goto return_status;

   for (i = 0; i < num_src; i++)
     {
        double x[4], y[4];
        int ix = i % nx_src;
        int iy = i / nx_src;
        x[0] = x0 + ix * dx - pixel_xoverlap;
        x[1] = x[0] + (dx + 2*pixel_xoverlap);
        x[2] = x[1];
        x[3] = x[0];
        y[0] = y0 + iy * dy;
        y[1] = y[0];
        y[2] = y[0] + dy;
        y[3] = y[2];

        src_values[i] = 1.0;
        src_mask[i] = 0;

        if ((-1 == Pixel_list_set_vertices (src_pixel_list, i, 4, x, y))
            || (-1 == Pixel_list_set_src_index (src_pixel_list, i, i)))
          goto return_status;
     }

   dest_grid_params.xmin = xshift;
   dest_grid_params.xmax = xshift + 1.0;
   dest_grid_params.ymin = yshift;
   dest_grid_params.ymax = yshift + 1.0;
   dest_grid_params.nx = nx_dest;
   dest_grid_params.ny = ny_dest;

   if (-1 == Pixel_grid_arrays (&dest_grid_params, &xcnr, &ycnr))
     goto return_status;

   if (NULL == (r = Pixel_open_regrid (&dest_grid_params, NULL)))
     goto return_status;

   if (-1 == (num_overlaps = Pixel_find_overlaps (r, src_pixel_list, NULL)))
     goto return_status;

   if (expect_no_overlaps && (num_overlaps == 0))
     {
        status = 0;
        goto return_status;
     }

   if (NULL == (info = (Pixel_Overlap_Info_Type *) MALLOC (num_dest * sizeof(*info))))
     goto return_status;

   /* Pixel_regrid assumes that dest_values initialized to fill value */
   for (i = 0; i < num_dest; i++)
     {
        dest_values[i] = DBL_MAX;
     }

   if (-1 == Pixel_regrid (r, src_values, src_mask, dest_values, info))
     goto return_status;

   src_sum = 0.0;
   for (i = 0; i < num_src; i++)
     {
        debug_print (stderr, "src[%2d] = %15.6e\n", i, src_values[i]);
        src_sum += src_values[i];
     }

   dest_sum = 0.0;
   if (info != NULL)
     {
        for (i = 0; i < num_dest; i++)
          {
             Pixel_Overlap_Info_Type *oi = info + i;
             debug_print (stderr, "dest[%2d] = %15.6e  %d overlaps\n",
                          i, dest_values[i],
                          oi->num_overlaps);
             if (dest_values[i] != DBL_MAX)
               {
                  dest_sum += dest_values[i];
               }
          }
     }

   expected_dest_sum = ((num_dest * (src_sum/num_src)
                         * (1.0 - fabs(xshift)) * (1.0 - fabs(yshift))));
   if ((expected_dest_sum < 0.0) || (info == NULL))
     expected_dest_sum = 0.0;
   debug_print (stdout, "dest_sum = %g  expected_dest_sum = %g\n",
                dest_sum, expected_dest_sum);
   if (fabs(dest_sum - expected_dest_sum) > DBL_EPSILON * expected_dest_sum)
     {
        fprintf (stderr, "*** FAIL: dest_sum=%g (expected %g)\n",
                 dest_sum, expected_dest_sum);
     }

   status = 0;
return_status:

   FREE(src_values);
   FREE(src_mask);
   FREE(dest_values);
   FREE(info);
   FREE(xcnr);
   FREE(ycnr);
   Pixel_list_free (src_pixel_list);
   Pixel_close_regrid (r);
   return status;
}

#define N1 32
#define N2 64

int main (void)
{
   double pixel_xoverlap;

   /* src/dest grids are identical */
   if (test_regrid (N1, N1, 1.0, 0.0,  0.0, 0.0, 0))
     return 1;

   /* dest grid is finer */
   if (test_regrid (N1, N1, 0.25, 0.0,  0.0, 0.0, 0))
     return 1;

   /* source grid pixels overlap in X, dest grid is finer/coarser */
   pixel_xoverlap = (4.0/116.0)*(1.0/N1);
   if ((0 != test_regrid (N1, N1, 0.25, pixel_xoverlap,  0.0, 0.0, 0))
       ||(0 != test_regrid (N1, N1, 4.0, pixel_xoverlap,  0.0, 0.0, 0)))
     return 1;

   /* dest grid is coarser (test realloc_overlap) */
   if ((0 != test_regrid (N1, N2, 4.0, 0.0,  0.0, 0.0, 0))
       || (0 != test_regrid (N2, N1, 4.0, 0.0,  0.0, 0.0, 0)))
     return 1;

   /* map to offset grids to improve test coverage */
   if ((0 != test_regrid (N1, N1, 1.0, 0.0,  0.25, 0.25, 0))  /* overlap */
       || (0 != test_regrid (N1, N1, 1.0, 0.0,  2.0, 2.0, 1))) /* no grid overlap */
     return 1;

   return 0;
}
