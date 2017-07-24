#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <proj_api.h>
#include <geodesic.h>

#include <libconfig.h>
#include <tell.h>

typedef struct
{
   int num_dim0;
   int num_dim1;
}
Window_Type;

#define FILTER_PRIVATE_DATA \
   Window_Type block; \
   Window_Type smooth1; \
   Window_Type smooth2; \
   Window_Type clip; \
   double clip_stddev;
#include "filter.h"

#include "proj.h"

#define MAX(a,b)  (((a)>(b))?(a):(b))
#define MIN(a,b)  (((a)<(b))?(a):(b))

/* WGS84 ellipsoid definition */
#define WGS84_SEMIMAJOR_AXIS     6378137.0        /* meters */
#define WGS84_FLATTENING_FACTOR  (1/298.257223563)

typedef struct
{
   gsl_interp *obj;
   gsl_interp_accel *acc;
   const gsl_interp_type *type;
   double *x;
   double *f;
   int num_finite;
   int num_alloc;
}
Interp_Type;

static double *pixel_area_vs_latitude (const Pixel_Grid_Param_Type *mesh)
{
   Pixel_Grid_Param_Type tmp_mesh;
   struct geod_geodesic g;
   double *areas = NULL;
   double *lons = NULL;
   double *lats = NULL;
   double dx;
   int i, num_pixel_vertices, status = -1;

   if (NULL == (areas = (double *)MALLOC (mesh->ny * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   /* struct copy */
   tmp_mesh = *mesh;

   /* Define new grid with dimensions (ny x 1).
    * Densify pixel edges with constant latitude (num_extra_xpoints > 0).
    * Pixel edges with constant longitude are geodesics and
    * need not be densified (num_extra_ypoints=0) */
   dx = (mesh->xmax - mesh->xmin) / mesh->nx;
   tmp_mesh.xmax = tmp_mesh.xmin + dx;
   tmp_mesh.nx = 1;
   tmp_mesh.num_extra_xpoints = 32;  /* FIXME: justify this */
   tmp_mesh.num_extra_ypoints =  0;

   num_pixel_vertices =
     (4 + 2 * tmp_mesh.num_extra_xpoints + 2 * tmp_mesh.num_extra_ypoints);

   /* compute pixel corner coordinates for (ny x 1) grid */
   if (0 != Pixel_grid_arrays (&tmp_mesh, &lons, &lats))
     goto free_and_return;

   /* Compute pixel areas.
    * What's the best way to do this?
    * The Earth is an ellipsoid, not a sphere.
    * The northern and southern pixel boundaries are not geodesics.
    * To compute the pixel area weights, we have at least 4 choices:
    *   - Approximate the earth as a sphere and use the analytic pixel areas.
    *   - Approximate the north and south pixel boundaries as piecewise
    *     geodesics, and use the proj library function.
    *   - Evaluate a 2D numerical integral for each pixel to compute
    *     the pixel area.
    *   - Perform the area-weighted regridding in Albers equal-area coordinates
    *     (e.g. use Pixel_regrid)
    * Considering the trade-off of accuracy, speed, efficiency,
    * and ease of implementation -- which method is best?
    *
    * 2D numerical integrals are a bit pricey, but we don't need many.
    * Doing the calculation in Albers coordinates is relatively slow
    * and memory intensive.
    * For small pixels, getting the pixel boundary a little bit wrong is
    * typically a 2nd order effect, and we can dramatically reduce that
    * error by densifying the polygon boundaries.
    *
    * For now, let's use the OTS library and approximate the
    * (appropriately densified) pixel boundaries as geodesics.
    * This should be fast, accurate and simple when integer
    * binning factors are used.
    */

   geod_init (&g, WGS84_SEMIMAJOR_AXIS, WGS84_FLATTENING_FACTOR);

   for (i = 0; i < mesh->ny; i++)
     {
        double geodesic_area, geodesic_perimeter;
        double *lats_i = lats + i * num_pixel_vertices;
        double *lons_i = lons + i * num_pixel_vertices;
        geod_polygonarea (&g, lats_i, lons_i, num_pixel_vertices,
                          &geodesic_area, &geodesic_perimeter);
        areas[i] = geodesic_area / 1.e6;  /* km^2 */
     }

   status = 0;
free_and_return:
   FREE(lons);
   FREE(lats);
   if (status)
     {
        FREE(areas);
        areas = NULL;
     }

   return areas;
}

static double *block_regrid_2d (const double *a,
                                const Pixel_Grid_Param_Type *mesh,
                                int dim0_blk, int dim1_blk)
{
   double nan_value = nan("");
   double *a_avg = NULL;
   double *area_sum = NULL;
   double *areas = NULL;
   int dim0, dim1, num_blk_pixels, num_bin0, num_bin1;
   int i, j, ii, jj;
   int status = -1;

   dim0 = mesh->ny;
   dim1 = mesh->nx;

   num_bin0 = dim0 / dim0_blk;
   num_bin1 = dim1 / dim1_blk;

   if ((num_bin0 * dim0_blk != dim0)
       || (num_bin1 * dim1_blk != dim1))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: array dimensions of source (%d,%d) are not integer multiples of target dimensions (%d, %d)",
                     __func__, dim0, dim1, dim0_blk, dim1_blk);
        return NULL;
     }

   if (NULL == (areas = pixel_area_vs_latitude (mesh)))
     return NULL;

   /* compute block sums */

   num_blk_pixels = dim0_blk * dim1_blk;
   if ((NULL == (a_avg = (double *)MALLOC (num_blk_pixels * sizeof(double))))
       ||(NULL == (area_sum = (double *)MALLOC (num_blk_pixels * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }
   memset ((char *)a_avg, 0, num_blk_pixels * sizeof(double));
   memset ((char *)area_sum, 0, num_blk_pixels * sizeof(double));

   ii = 0;
   jj = 0;

   for (i = 0; i < dim0; i++)
     {
        double pixel_area = areas[i];
        const double *a_row = a + i * dim1;
        double *area_sum_row;
        double *a_avg_row;

        ii = i / num_bin0;
        a_avg_row = a_avg + ii * dim1_blk;
        area_sum_row = area_sum + ii * dim1_blk;
        for (j = 0; j < dim1; j++)
          {
             double a_ij = a_row[j];
             if (isfinite(a_ij) && (a_ij != DBL_MAX))
               {
                  jj = j / num_bin1;
                  a_avg_row[jj] += a_ij * pixel_area;
                  area_sum_row[jj] += pixel_area;
               }
          }
     }

   /* compute area-weighted block averages */

   for (ii = 0; ii < num_blk_pixels; ii++)
     {
        if ((area_sum[ii] != 0) && (a_avg[ii] != 0))
          {
             a_avg[ii] /= area_sum[ii];
          }
        else a_avg[ii] = nan_value;
     }

   status = 0;
free_and_return:
   FREE(area_sum);
   FREE(areas);
   if (status)
     {
        FREE(a_avg);
        a_avg = NULL;
     }

   return a_avg;
}

static double *boxcar_smooth (const double *a, int dim0, int dim1,
                              const Window_Type *smo, double *psmoothed)
{
   double nan_value = nan("");
   double *smoothed = NULL;
   double sij;
   int i, nij;
   int hw0 = smo->num_dim0/2;
   int hw1 = smo->num_dim1/2;

   if (psmoothed == NULL)
     {
        if (NULL == (smoothed = (double *) MALLOC (dim0 * dim1 * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return NULL;
          }
     }
   else if (a == psmoothed)
     {
        tell_verror (TELL_USAGE_ERROR,
                     "%s: input/output array addresses must be distinct",
                     __func__);
        return NULL;
     }
   else smoothed = psmoothed;

   for (i = 0; i < dim0; i++)
     {
        int j, ib, ie;

        ib = MAX(i - hw0, 0);
        ie = MIN(i + hw0, dim0);

        for (j = 0; j < dim1; j++)
          {
             int jb, je, ii, jj;
             double a_ij;

             jb = MAX(j - hw1, 0);
             je = MIN(j + hw1, dim1);

             nij = 0;
             sij = 0.0;

             for (ii = ib; ii < ie; ii++)
               {
                  const double *a_row = a + ii * dim1;
                  for (jj = jb; jj < je; jj++)
                    {
                       a_ij = a_row[jj];
                       if (isfinite(a_ij) && (a_ij != DBL_MAX))
                         {
                            sij += a_ij;
                            nij++;
                         }
                    }
               }

             if (nij > 0)
               smoothed[j + i * dim1] = sij / nij;
             else
               smoothed[j + i * dim1] = nan_value;
          }
     }

   return smoothed;
}

static void replace_masked_pixels (const double *src, int dim0, int dim1,
                                   double *dest)
{
   int i, n = dim0 * dim1;

   for (i = 0; i < n; i++)
     {
        if (isnan(dest[i]))
          {
             dest[i] = src[i];
          }
     }
}

static double *clip_hotspots (const double *a, int dim0, int dim1,
                              const Window_Type *clp, double a_clip_stddev,
                              double *pclipped)
{
   double *clipped = NULL;
   double nan_value = nan("");
   int hw0 = clp->num_dim0/2;
   int hw1 = clp->num_dim1/2;
   int i, nij;

   if (pclipped == NULL)
     {
        if (NULL == (clipped = (double *) MALLOC (dim0 * dim1 * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return NULL;
          }
     }
   else if (a != pclipped)
     {
        clipped = pclipped;
        memcpy ((char *)clipped, (char *)a, dim0 * dim1 * sizeof(double));
     }
   else
     {
        tell_verror (TELL_USAGE_ERROR,
                     "%s: input/output array addresses must be distinct",
                     __func__);
        return NULL;
     }

   for (i = 0; i < dim0; i++)
     {
        const double *a_row;
        double sij, a_avg, a_sum_sqdev, a_stddev;
        int j, ib, ie;

        ib = MAX(i - hw0, 0);
        ie = MIN(i + hw0, dim0);

        for (j = 0; j < dim1; j++)
          {
             int jb, je, ii, jj;
             double a_ij;

             jb = MAX(j - hw1, 0);
             je = MIN(j + hw1, dim1);

             /* compute mean inside the window */

             sij = 0.0;
             nij = 0;
             for (ii = ib; ii < ie; ii++)
               {
                  a_row = a + ii * dim1;
                  for (jj = jb; jj < je; jj++)
                    {
                       a_ij = a_row[jj];
                       if (isfinite(a_ij) && (a_ij != DBL_MAX))
                         {
                            sij += a_ij;
                            nij++;
                         }
                    }
               }

             /* When the window contains too few valid values,
              * we have nothing to clip */
             if (nij < 2)
               continue;

             a_avg = sij / nij;

             /* compute standard deviation inside the window */

             a_sum_sqdev = 0.0;
             for (ii = ib; ii < ie; ii++)
               {
                  a_row = a + ii * dim1;
                  for (jj = jb; jj < je; jj++)
                    {
                       a_ij = a_row[jj];
                       if (isfinite(a_ij) && (a_ij != DBL_MAX))
                         {
                            double dev = a_ij - a_avg;
                            a_sum_sqdev += dev * dev;
                         }
                    }
               }

             a_stddev = sqrt (a_sum_sqdev / (nij - 1));

             /* If center pixel is high, replace with NaN */
             if (a[j + i * dim1] > (a_avg + a_clip_stddev * a_stddev))
               {
                  clipped[j + i * dim1] = nan_value;
               }
          }
     }

   return clipped;
}

static void interp_free (Interp_Type *it)
{
   if (it == NULL)
     return;
   gsl_interp_free (it->obj);
   gsl_interp_accel_free (it->acc);
   FREE(it->x);
   FREE(it->f);
   FREE(it);
}

static Interp_Type *interp_alloc (int n)
{
   Interp_Type *it = NULL;

   if (NULL == (it = (Interp_Type *)MALLOC (sizeof *it)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)it, 0, sizeof *it);

   it->type = gsl_interp_linear;
   it->num_alloc = n;

   if ((NULL == (it->x = (double *) MALLOC (n * sizeof(double))))
       || (NULL == (it->f = (double *) MALLOC (n * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        interp_free (it);
        return NULL;
     }

   if ((NULL == (it->obj = gsl_interp_alloc (it->type, n)))
       || (NULL == (it->acc = gsl_interp_accel_alloc ())))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        interp_free (it);
        return NULL;
     }

   return it;
}

static int interp_init (Interp_Type **p_it, const double *x, const double *f, int n)
{
   Interp_Type *it = *p_it;
   int i, k, num_finite;
   int gsl_errno;

   /* Filter nans and infs, then initialize the interpolation using
    * only the finite f[i] values.
    */

   num_finite = 0;
   for (i = 0; i < n; i++)
     {
        if (isfinite (f[i]))
          num_finite++;
     }

   it->num_finite = num_finite;

   if (num_finite < 2)
     return 0;

   /* Ideally, we would reallocate the Interp_Type only when we needed
    * more space (e.g. num_finite > it->num_alloc).  Unfortunately,
    * gsl_interp_init returns an error when called with fewer points
    * than were allocated by gsl_interp_alloc.  Why?  This means that
    * we can't separate allocation from initialization in an efficient way,
    * and we're forced to reallocate whenever (num_finite != it->num_alloc).
    * Sigh.
    */
   if (num_finite != it->num_alloc)
     {
        Interp_Type *q;
        if (NULL == (q = interp_alloc (num_finite)))
          return -1;
        q->num_finite = num_finite;
        interp_free (it);
        it = q;
        *p_it = it;
     }

   k = 0;
   for (i = 0; i < n; i++)
     {
        if (isfinite(f[i]))
          {
             it->x[k] = x[i];
             it->f[k] = f[i];
             k += 1;
          }
     }

   if (0 != (gsl_errno = gsl_interp_init (it->obj, it->x, it->f, it->num_finite)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: preparing to interpolate (%s)",
                     __func__, gsl_strerror (gsl_errno));
        return -1;
     }

   return 0;
}

static int interp_values (Interp_Type *it, double *xv, double *fv, int nv)
{
   int i;

   if (it->num_finite < 2)
     {
        double nan_value = nan("");
        for (i = 0; i < nv; i++)
          {
             fv[i] = nan_value;
          }
        return 0;
     }

   for (i = 0; i < nv; i++)
     {
        int gsl_errno =  gsl_interp_eval_e (it->obj, it->x, it->f,
                                            xv[i], it->acc, &fv[i]);
        if ((0 != gsl_errno) && (gsl_errno != GSL_EDOM))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: interpolation failed (%s)",
                          __func__, gsl_strerror (gsl_errno));
             return -1;
          }
     }

   return 0;
}

static int unblock_avg_2d (const double *a_blk, int dim0_blk, int dim1_blk,
                           const Pixel_Grid_Param_Type *mesh, double *a)
{
   Interp_Type *it = NULL;
   double *xc_blk = NULL;
   double *yc_blk = NULL;
   double *xx = NULL;
   double *yy = NULL;
   double *v0 = NULL;
   double *a_tmp = NULL;
   double dx_blk, dy_blk, xc_min_blk, yc_min_blk;
   double dx, dy, xc_min, yc_min;
   int i, k, dim0, dim1, n0, n1;
   int status = -1;

   /* y is index 0,
    * x is index 1 */
   dim0 = mesh->ny;
   dim1 = mesh->nx;

   n0 = dim0 / dim0_blk;
   n1 = dim1 / dim1_blk;

   if ((dim0 != n0 * dim0_blk) || (dim1 != n1 * dim1_blk))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: result dimensions (%d, %d) are not integer multiples of source dimensions (%d, %d)",
                     __func__, dim0, dim1, dim0_blk, dim1_blk);
        return -1;
     }

   if ((NULL == (xc_blk = (double *)MALLOC (dim1_blk * sizeof(double))))
       || (NULL == (yc_blk = (double *)MALLOC (dim0_blk * sizeof(double))))
       || (NULL == (xx = (double *)MALLOC (dim1 * sizeof(double))))
       || (NULL == (yy = (double *)MALLOC (dim0 * sizeof(double))))
       || (NULL == (v0 = (double *)MALLOC (dim0 * sizeof(double))))
       || (NULL == (a_tmp = (double *)MALLOC (dim0 * dim1_blk * sizeof(double))))
      )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   /* Generate pixel center coordinate arrays for coarse grid */
   dx_blk = (mesh->xmax - mesh->xmin) / dim1_blk;
   dy_blk = (mesh->ymax - mesh->ymin) / dim0_blk;
   xc_min_blk = mesh->xmin + 0.5 * dx_blk;
   yc_min_blk = mesh->ymin + 0.5 * dy_blk;
   for (i = 0; i < dim0_blk; i++)
     {
        yc_blk[i] = yc_min_blk + i * dy_blk;
     }
   for (i = 0; i < dim1_blk; i++)
     {
        xc_blk[i] = xc_min_blk + i * dx_blk;
     }

   /* Generate pixel center coordinate arrays for fine grid */
   dx = (mesh->xmax - mesh->xmin) / mesh->nx;
   dy = (mesh->ymax - mesh->ymin) / mesh->ny;
   xc_min = mesh->xmin + 0.5 * dx;
   yc_min = mesh->ymin + 0.5 * dy;
   for (i = 0; i < dim0; i++)
     {
        yy[i] = yc_min + i * dy;
     }
   for (i = 0; i < dim1; i++)
     {
        xx[i] = xc_min + i * dx;
     }

   /* Following the IDL prototype, we implement the 2D interpolation
    * as a series of 1D interpolations.
    * Where data is missing, we interpolate across the gaps.
    */

   if (NULL == (it = interp_alloc (MAX(dim0_blk, dim1_blk))))
     goto free_and_return;

   gsl_set_error_handler_off();

   /* 1D interpolation along latitude:
    *       a_blk[dim0_blk,dim1_blk] -> a_tmp[dim0,dim1_blk]
    * where the dimension order is [column, row], and
    * the rightmost index (dimension dim1_blk) changes fastest.
    */

   for (k = 0; k < dim1_blk; k++)
     {
        for (i = 0; i < dim0_blk; i++)
          {
             v0[i] = a_blk[k + i * dim1_blk];
          }

        if (0 != interp_init (&it, yc_blk, v0, dim0_blk))
          goto free_and_return;

        /* overwrite v0 with interpolated values */
        if (0 != interp_values (it, yy, v0, dim0))
          goto free_and_return;

        for (i = 0; i < dim0; i++)
          {
             a_tmp[k + i * dim1_blk] = v0[i];
          }
     }

   /* 1D interpolation along longitude:
    *       a_tmp[dim0,dim1_blk] -> a[dim0, dim1]
    * where the dimension order is [column, row], and
    * the rightmost index (dimension dim1_blk) changes fastest.
    */

   for (i = 0; i < dim0; i++)
     {
        const double *a_tmp_i = a_tmp + i * dim1_blk;
        double *a_i = a + i * dim1;

        if (0 != interp_init (&it, xc_blk, a_tmp_i, dim1_blk))
          goto free_and_return;

        if (0 != interp_values (it, xx, a_i, dim1))
          goto free_and_return;
     }

   status = 0;
free_and_return:

   interp_free (it);
   FREE(xc_blk);
   FREE(yc_blk);
   FREE(xx);
   FREE(yy);
   FREE(v0);
   FREE(a_tmp);

   return status;
}

static int filter_vert_strat (const Filter_Type *flt,
                              const Pixel_Grid_Param_Type *mesh,
                              double *vert_strat)
{
   const Window_Type *blk = &flt->block;
   double *tmp1 = NULL;
   double *tmp2 = NULL;
   int nb0, nb1;

   /* Regrid from the fine mesh to a coarse mesh */

   nb0 = mesh->ny / blk->num_dim0;
   nb1 = mesh->nx / blk->num_dim1;

   if (NULL == (tmp1 = block_regrid_2d (vert_strat, mesh, nb0, nb1)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: block_regrid_2d failed", __func__);
        goto return_error;
     }

   if (NULL == (tmp2 = boxcar_smooth (tmp1, nb0, nb1, &flt->smooth1, NULL)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: boxcar_smooth failed (#1)", __func__);
        goto return_error;
     }

   /* result over-writes tmp1 */
   replace_masked_pixels (tmp2, nb0, nb1, tmp1);

   /* result over-writes tmp2 */
   if (NULL == clip_hotspots (tmp1, nb0, nb1, &flt->clip, flt->clip_stddev, tmp2))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: clipping hot spots (#1)", __func__);
        goto return_error;
     }

   /* result over-writes tmp1 */
   if (NULL == clip_hotspots (tmp2, nb0, nb1, &flt->clip, flt->clip_stddev, tmp1))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: clipping hot spots (#2)", __func__);
        goto return_error;
     }

   /* result over-writes tmp2 */
   if (NULL == boxcar_smooth (tmp1, nb0, nb1, &flt->smooth2, tmp2))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: boxcar_smooth failed (#2)", __func__);
        goto return_error;
     }

   /* Interpolate result back to the original fine mesh */

   if (0 != unblock_avg_2d (tmp2, nb0, nb1, mesh, vert_strat))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unblocking array", __func__);
        goto return_error;
     }

   FREE(tmp1);
   FREE(tmp2);
   return 0;

return_error:
   FREE(tmp1);
   FREE(tmp2);
   return -1;
}

static void filter_delete (Filter_Type *flt)
{
   if (flt == NULL)
     return;
   FREE(flt);
}

static Filter_Type *new_filter (void)
{
   Filter_Type *flt = NULL;

   if (NULL == (flt = (Filter_Type *)MALLOC (sizeof *flt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)flt, 0, sizeof *flt);

   flt->filter_delete = filter_delete;
   flt->filter_apply = filter_vert_strat;

   return flt;
}

static int lookup_filter_window (const config_setting_t *s,
                                 int *num_dim0, int *num_dim1)
{
   if ((CONFIG_TRUE != config_setting_lookup_int (s, "num_dim0", num_dim0))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num_dim1", num_dim1)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining filter window", __func__);
        return -1;
     }
   return 0;
}

static int read_config_file (config_t *cfg, Filter_Type *flt)
{
   config_setting_t *setting, *s;

   if (NULL == (setting = config_lookup (cfg, "filter")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing filter definition in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (s = config_setting_get_member (setting, "block")))
       || (-1 == lookup_filter_window (s, &flt->block.num_dim0, &flt->block.num_dim1)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining blocking window", __func__);
        return -1;
     }
   if ((NULL == (s = config_setting_get_member (setting, "smooth1")))
       || (-1 == lookup_filter_window (s, &flt->smooth1.num_dim0, &flt->smooth1.num_dim1)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining smoothing window", __func__);
        return -1;
     }
   if ((NULL == (s = config_setting_get_member (setting, "smooth2")))
       || (-1 == lookup_filter_window (s, &flt->smooth2.num_dim0, &flt->smooth2.num_dim1)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining smoothing window", __func__);
        return -1;
     }
   if ((NULL == (s = config_setting_get_member (setting, "clip")))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "stddev", &flt->clip_stddev))
       || (-1 == lookup_filter_window (s, &flt->clip.num_dim0, &flt->clip.num_dim1)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining clipping window parameters", __func__);
        return -1;
     }

   return 0;
}

Filter_Type *filter_open (config_t *cfg)
{
   Filter_Type *flt = NULL;

   if (NULL == (flt = new_filter ()))
     return NULL;

   if (0 != read_config_file (cfg, flt))
     {
        filter_delete (flt);
        return NULL;
     }

   return flt;
}
