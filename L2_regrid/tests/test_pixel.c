#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netcdf.h>

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

/* Note that we don't check the output of Pixel_regrid_bytes...
 * These bitfield tests are mainly intended to improve code coverage.
 * The floating point regridding logic is verified by other tests,
 * and the bitfields share a lot of that.
 */
#define BITFIELD_TEST_TYPE(type,fillsym,typestr) \
static int type##_bitfield_test (const Pixel_Regrid_Type *r, int num_src, int num_dest, \
                                 const int *src_mask) \
{ \
   typestr *a = NULL, *src, *dest; \
   typestr fill_value = NC_FILL_##fillsym; \
   int i, num = num_src + num_dest; \
 \
   if (NULL == (a = (typestr *)MALLOC (num * sizeof(*a)))) \
     return -1; \
   src = a; \
   dest = a + num_src; \
 \
   for (i = 0; i < num_src; i++) \
     { \
        src[i] = 1; \
     } \
   for (i = 0; i < num_dest; i++) \
     { \
        dest[i] = fill_value; \
     } \
 \
   if (-1 == Pixel_regrid_bytes (r, src_mask, VALUE_IS_##fillsym, &fill_value, \
                                 src, dest)) \
     { \
        FREE(a); \
        return -1; \
     } \
 \
   FREE(a); \
   return 0; \
}

BITFIELD_TEST_TYPE(ul,UINT64,unsigned long long)
BITFIELD_TEST_TYPE(ui,UINT,unsigned int)
BITFIELD_TEST_TYPE(us,USHORT,unsigned short)
BITFIELD_TEST_TYPE(uc,UBYTE,unsigned char)
BITFIELD_TEST_TYPE(l,INT64,long long)
BITFIELD_TEST_TYPE(i,INT,int)
BITFIELD_TEST_TYPE(s,SHORT,short)
BITFIELD_TEST_TYPE(c,BYTE,char)

typedef int Bitfield_Test_Type (const Pixel_Regrid_Type *, int, int, const int *);

static Bitfield_Test_Type *Bitfield_Tests[] = {
   &ul_bitfield_test,
   &ui_bitfield_test,
   &us_bitfield_test,
   &uc_bitfield_test,
   &l_bitfield_test,
   &i_bitfield_test,
   &s_bitfield_test,
   &c_bitfield_test,
   NULL
};

static int test_bitfields (const Pixel_Regrid_Type *r, int num_src, int num_dest,
                           const int *src_mask)
{
   Bitfield_Test_Type **btest;
   for (btest = Bitfield_Tests; *btest != NULL; btest++)
     {
        if (-1 == (*btest) (r, num_src, num_dest, src_mask))
          return -1;
     }
   return 0;
}

static int test_regrid (int nx_src, int ny_src, float bin_factor,
                        double pixel_xoverlap,
                        double xshift, double yshift,
                        int expect_no_overlaps)
{
   Pixel_Grid_Param_Type dest_grid_params;
   Pixel_List_Type *src_pixel_list = NULL;
   Pixel_Regrid_Type *r = NULL;
   double *src_values=NULL, *dest_values=NULL;
   double x0 = 0.0, dx = 1.0/nx_src;
   double y0 = 0.0, dy = 1.0/ny_src;
   double src_sum, dest_sum, expected_dest_sum;
   double *xcnr = NULL, *ycnr = NULL;
   int num_src, nx_dest, ny_dest, num_dest, num_overlaps;
   int new_num_step, new_num_xtrack;
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

   if (-1 == test_bitfields (r, num_src, num_dest, src_mask))
     {
        fprintf (stderr, "*** Error: bitfield test failed\n");
        goto return_status;
     }

   if (-1 == Pixel_regrid (r, src_mask, DBL_MAX, src_values, dest_values, NULL))
     goto return_status;

   src_sum = 0.0;
   for (i = 0; i < num_src; i++)
     {
        debug_print (stderr, "src[%2d] = %15.6e\n", i, src_values[i]);
        src_sum += src_values[i];
     }

   dest_sum = 0.0;
   for (i = 0; i < num_dest; i++)
     {
        debug_print (stderr, "dest[%2d] = %15.6e\n",
                     i, dest_values[i]);
        if (dest_values[i] != DBL_MAX)
          {
             dest_sum += dest_values[i];
          }
     }

   expected_dest_sum = ((num_dest * (src_sum/num_src)
                         * (1.0 - fabs(xshift)) * (1.0 - fabs(yshift))));
   if (expected_dest_sum < 0.0)
     expected_dest_sum = 0.0;
   debug_print (stdout, "dest_sum = %g  expected_dest_sum = %g\n",
                dest_sum, expected_dest_sum);
   if (fabs(dest_sum - expected_dest_sum) > DBL_EPSILON * expected_dest_sum)
     {
        fprintf (stderr, "*** FAIL: dest_sum=%g (expected %g)\n",
                 dest_sum, expected_dest_sum);
        goto return_status;
     }

   Pixel_regrid_grow_srcdims (r, INT_MAX-1, INT_MAX-1);
   Pixel_regrid_get_srcdims (r, &new_num_step, &new_num_xtrack);
   if ((new_num_step != INT_MAX) || (new_num_xtrack != INT_MAX))
     {
        fprintf (stderr, "*** FAIL: new_num_step=%d new_num_xtrack=%d, expected INT_MAX\n",
                 new_num_step, new_num_xtrack);
        goto return_status;
     }

   status = 0;
return_status:

   FREE(src_values);
   FREE(src_mask);
   FREE(dest_values);
   FREE(xcnr);
   FREE(ycnr);
   Pixel_list_free (src_pixel_list);
   Pixel_close_regrid (r);
   return status;
}

static int test_regrid_stat (int nx_src, int ny_src)
{
   Pixel_Grid_Param_Type dest_grid_params;
   Pixel_List_Type *src_pixel_list = NULL;
   Pixel_Regrid_Type *r = NULL;
   Pixel_Regrid_Stats_Type *regrid_stats = NULL;
   double *src_values=NULL, *dest_values=NULL;
   double x0 = 0.0, dx = 1.0/nx_src;
   double y0 = 0.0, dy = 1.0/ny_src;
   double *xcnr = NULL, *ycnr = NULL;
   /* These values define the test case -- if you change these
    * values the expected test values must be updated!! */
   float bin_factor = 2.0;          /* dest grid is coarser */
   double xshift = 0.5 * dx;
   double yshift = 0.0;
   double pixel_xoverlap = 0.0;
   int expect_no_overlaps = 0;
   /* end of test case parameters */
   int num_src, nx_dest, ny_dest, num_dest, num_overlaps;
   int i, *src_mask= NULL;
   int status = -1;

   if (DEBUG) fprintf (stderr, "nx_src = %d\n", nx_src);

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

        src_values[i] = 1.0 + (ix % 2);
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

   if (NULL == (regrid_stats = Pixel_alloc_regrid_stats (num_dest, 1)))
     goto return_status;

   if (-1 == Pixel_regrid (r, src_mask, DBL_MAX, src_values, dest_values, regrid_stats))
     goto return_status;

   for (i = 0; i < num_dest; i++)
     {
        struct {double min, max; int num;} expect;
        int ix = i % nx_dest;
        int iy = i / nx_dest;
        int expected;

        if (DEBUG) fprintf (stderr, "(%2d,%2d): %5.3e %5.3e %d\n",
                            ix, iy,
                            regrid_stats->min[i],
                            regrid_stats->max[i],
                            regrid_stats->num[i]);

        if (ix < nx_dest-1)
          {
             /* Apart from edge effects, we expect each destination
              * pixel to have contributions from 4 source pixels,
              * with the minimum source pixel value = 1 */
             expect.min = 1;
             expect.max = 2;
             expect.num = 4;
          }
        else
          {
             /* Handle edge cases:
              * When nx_src is even, the last column will have value=2,
              * so when the dest grid has a 1/2 pixel shift in the +X direction
              * the last column will see only 2 contributions with value=2,
              * and no contributions with value=1.
              * When nx_src is odd, vice-versa.
              */
             expect.num = 2;
             if ((nx_src / 2)*2 == nx_src)
               {
                  /* nx_src even */
                  expect.min = 2;
                  expect.max = 2;
               }
             else
               {
                  /* nx_src odd */
                  expect.min = 1;
                  expect.max = 1;
               }
          }

        expected = ((regrid_stats->min[i] == expect.min)
                    && (regrid_stats->max[i] == expect.max)
                    && (regrid_stats->num[i] == expect.num));
        if (!expected)
          {
             fprintf (stderr,
                      "*** unexpected regrid_stats: i=%d min=%g\tmax=%g\tnum=%d\n"
                      "                   expected:      min=%g\tmax=%g\tnum=%d\n",
                      i,
                      regrid_stats->min[i],
                      regrid_stats->max[i],
                      regrid_stats->num[i],
                      expect.min, expect.max, expect.num);
             goto return_status;
          }
     }

   status = 0;
return_status:

   FREE(src_values);
   FREE(src_mask);
   FREE(dest_values);
   FREE(xcnr);
   FREE(ycnr);
   Pixel_list_free (src_pixel_list);
   Pixel_close_regrid (r);
   Pixel_free_regrid_stats (regrid_stats);
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

   if (test_regrid_stat (N1, N1))
     return 1;
   if (test_regrid_stat (N1+1, N1))
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
