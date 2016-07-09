#include "defs.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tell.h>
#include <netcdf.h>
#include <tio.h>
#include <tio_template.h>
#include <proj_api.h>

#include "poly.h"
#include "pixel.h"
#include "regrid.h"

typedef struct
{
   /* polygon vertices in coordinates that facilitate coordinate lookup */
   Pixel_List_Type *pixel_lookup;
   /* polygon vertices in coordinates that facilitate area calculation */
   Pixel_List_Type *pixel_area;
   int max_step;
   int max_xtrack;
}
Source_Pixel_List_Type;

/* The Albers equal-area conic projection preserves areas
 * but in general, polygon edges do not project into
 * straight lines.  For "sufficiently small" polygons,
 * the polygon area computed _assuming_ straight edges
 * is "sufficiently accurate".
 */
static int longlat_to_albers (double *lon, double *lat, int n)
{
   projPJ albers = NULL, longlat = NULL;
   /* spatialreference.org ESRI Projection 102003 - USA Contiguous Albers Equal Area Conic */
   const char ctl_albers[] =
     "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs";
   /* spatialreference.org EPSG Projection 4326 - WGS 84  */
   const char ctl_longlat[] =
     "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
   int i, status = -1;

   if ((NULL == (albers = pj_init_plus (ctl_albers)))
       || (NULL == (longlat = pj_init_plus (ctl_longlat))))
     {
        Tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus failed", __func__);
        goto free_and_return;
     }

   for (i = 0; i < n; i++)
     {
        lon[i] *= DEG_TO_RAD;
        lat[i] *= DEG_TO_RAD;
     }

   status = pj_transform (longlat, albers, n, 1, lon, lat, NULL);
   if (status)
     {
        Tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        goto free_and_return;
     }

   status = 0;
free_and_return:
   pj_free (longlat);
   pj_free (albers);
   return status ? -1 : 0;
}

static int pack_pixel_list (Pixel_List_Type *pixel_list,
                            double *xs, double *ys, int num_pixels,
                            int *step, int num_xtrack)
{
   int i;

   for (i = 0; i < num_pixels; i++)
     {
        int pix_xtrack = i % num_xtrack;
        int pix_step_index = i / num_xtrack;
        int pix = pix_xtrack + step[pix_step_index] * num_xtrack;
        double *x = xs + 4*i;
        double *y = ys + 4*i;
        int j;

        for (j = 0; j < 4; j++)
          {
             if ((0 == isfinite(x[j])) || (0 == isfinite(y[j])))
               break;
          }
        if (j == 4)
          {
             if ((-1 == Pixel_list_set_vertices (pixel_list, i, 4, x, y))
                 || (-1 == Pixel_list_set_src_index (pixel_list, i, pix)))
               return -1;
          }
     }

   return 0;
}

static void free_source_pixel_list (Source_Pixel_List_Type *src)
{
   if (src == NULL)
     return;
   Pixel_list_free(src->pixel_area);
   Pixel_list_free(src->pixel_lookup);
   FREE(src);
}

static Source_Pixel_List_Type *new_source_pixel_list (void)
{
   Source_Pixel_List_Type *src = NULL;

   if (NULL == (src = (Source_Pixel_List_Type *) MALLOC (sizeof *src)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   src->pixel_area = NULL;
   src->pixel_lookup = NULL;
   src->max_step = 0;
   src->max_xtrack = 0;

   return src;
}

static int record_max_xtrack (int ncid, int num_xtrack, int *max_xtrack)
{
   int i, mx, start[2], count[2];
   int *xtrack = NULL;

   if (NULL == (xtrack = (int *) MALLOC (num_xtrack * sizeof (int))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = 0;
   count[0] = num_xtrack;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_XTRACK,
                                  start, count, TIO_INT, xtrack))
     {
        FREE(xtrack);
        return -1;
     }

   mx = 0;
   for (i = 0; i < num_xtrack; i++)
     {
        if (xtrack[i] > mx)
          mx = xtrack[i];
     }
   FREE(xtrack);

   *max_xtrack = mx;

   return 0;
}

static Source_Pixel_List_Type *
read_pixel_vertices (const char *file, const char *lonlat_grp)
{
   Source_Pixel_List_Type *src = NULL;
   TIO_Var_Info_Type vi;
   int ncid, grp, start[3], count[3];
   int num_steps, num_xtrack, num_pixels, num_sides, len_bounds;
   double *lon_bounds = NULL, *lat_bounds = NULL;
   double *albers_x_bounds, *albers_y_bounds;
   int *step = NULL;
   int i, status;

   if (-1 == TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if (-1 == TIO_inq_grp (ncid, lonlat_grp, &grp))
     return NULL;

   status = -1;

   /* assume xtrack, mirror_step dimensions and coordinate variables
    * are file-global */

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_XTRACK, &vi))
     goto free_and_return;
   num_xtrack = vi.dimlens[0];

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_STEP, &vi))
     goto free_and_return;
   num_steps = vi.dimlens[0];

   if (NULL == (src = new_source_pixel_list ()))
     goto free_and_return;

   if (NULL == (step = (int *) MALLOC (num_steps * sizeof (int))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   start[0] = 0;
   count[0] = num_steps;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_STEP,
                                  start, count, TIO_INT, step))
     goto free_and_return;

   for (i = 0; i < num_steps; i++)
     {
        if (step[i] > src->max_step)
          src->max_step = step[i];
     }

   /* this is a bit paranoid */
   if (-1 == record_max_xtrack (ncid, num_xtrack, &src->max_xtrack))
     goto free_and_return;

   /* read lon/lat bounds arrays */

   num_pixels = num_steps * num_xtrack;
   len_bounds = 4 * num_pixels * sizeof(double);

   if ((NULL == (lon_bounds = (double *) MALLOC (len_bounds)))
       || (NULL == (lat_bounds = (double *) MALLOC (len_bounds))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = num_steps;
   count[1] = num_xtrack;
   count[2] = 4;

   if ((-1 == TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS,
                                   start, count, TIO_DOUBLE, lon_bounds))
       || (-1 == TIO_get_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS,
                                      start, count, TIO_DOUBLE, lat_bounds)))
     {
        goto free_and_return;
     }

   /* Pack pixel vertices into pixel list structures */

   /* Zero-length polygons will indicate lines of sight that
    * have invalid lon-lat coordinates -- usually because
    * they don't intersect the earth. */
   num_sides = 0;

   if ((NULL == (src->pixel_area = Pixel_list_new (num_pixels, num_sides)))
       || (NULL == (src->pixel_lookup = Pixel_list_new (num_pixels, num_sides))))
     goto free_and_return;

   if ((-1 == Pixel_list_use_src_index (src->pixel_area))
       || (-1 == Pixel_list_use_src_index (src->pixel_lookup)))
     goto free_and_return;

   /* pixel_lookup coordinates will be used to determine which
    * destination pixels overlap each source pixel */
   if (-1 == pack_pixel_list (src->pixel_lookup,
                              lon_bounds, lat_bounds, num_pixels,
                              step, num_xtrack))
     {
        goto free_and_return;
     }

   if (-1 == longlat_to_albers (lon_bounds, lat_bounds, 4*num_pixels))
     goto free_and_return;
   /* NOTE: coordinate projection is done in place, so after the call,
    * (lon_bounds, lat_bounds) [deg] is really Albers (x,y) [meters] */
   albers_x_bounds = lon_bounds;
   albers_y_bounds = lat_bounds;

   /* pixel_area coordinates will be used to compute pixel overlap areas */
   if (-1 == pack_pixel_list (src->pixel_area,
                              albers_x_bounds, albers_y_bounds, num_pixels,
                              step, num_xtrack))
     {
        goto free_and_return;
     }

   status = 0;
free_and_return:
   (void) TIO_close (ncid);
   FREE(lon_bounds);
   FREE(lat_bounds);
   FREE(step);
   if (status)
     {
        free_source_pixel_list (src);
        src = NULL;
     }

   return src;
}

static int
find_all_pixel_overlaps (Pixel_Regrid_Type *r, char **files, int num_files,
                         const char *lonlat_grp)
{
   Source_Pixel_List_Type *src = NULL;
   int i;

   /* Loop over granule files, and accumulate contributions
    * to the pixel overlap array in each destination pixel.
    */

   for (i = 0; i < num_files; i++)
     {
        int num_overlaps;

        free_source_pixel_list (src);
        if (NULL == (src = read_pixel_vertices (files[i], lonlat_grp)))
          break;

        Pixel_regrid_grow_srcdims (r, src->max_step, src->max_xtrack);

        num_overlaps = Pixel_find_overlaps (r, src->pixel_area,
                                            src->pixel_lookup);
        if (num_overlaps < 0)
          break;
        else if (num_overlaps == 0)
          {
             /* FIXME:  should this be a warning message? */
             Tell_verror (TELL_UNKNOWN_ERROR,
                          "%s: no contribution to target grid from %s",
                          __func__, files[i]);
          }
     }

   free_source_pixel_list (src);

   if (i != num_files)
     {
        Tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: unexpected error on input %s",
                     __func__, files[i]);
        return -1;
     }

   return 0;
}

static Pixel_List_Type *
dest_pixel_area_coords (const Pixel_Grid_Param_Type *dest)
{
   Pixel_List_Type *pixel_list = NULL;
   int num_pixels = dest->nx * dest->ny;
   double *xs = NULL, *ys = NULL;
   double *x, *y;
   int i, status = -1;

   if (-1 == Pixel_grid_arrays (dest, &xs, &ys))
     return NULL;

   if (-1 == longlat_to_albers (xs, ys, 4*num_pixels))
     goto free_and_return;

   if (NULL == (pixel_list = Pixel_list_new (num_pixels, 4)))
     goto free_and_return;

   x = xs;
   y = ys;

   for (i = 0; i < num_pixels; i++)
     {
        if (-1 == Pixel_list_set_vertices (pixel_list, i, 4, x, y))
          goto free_and_return;
        x += 4;
        y += 4;
     }

   status = 0;
free_and_return:
   FREE(xs);
   FREE(ys);
   if (status)
     {
        Pixel_list_free (pixel_list);
        return NULL;
     }

   return pixel_list;
}

Pixel_Regrid_Type *
Regrid_open (const Pixel_Grid_Param_Type *dest,
             char **files, int num_files, const char *lonlat_grp)
{
   Pixel_Regrid_Type *r = NULL;
   Pixel_List_Type *dest_pixel_area = NULL;
   int status = -1;

   if (NULL == (dest_pixel_area = dest_pixel_area_coords (dest)))
     goto free_and_return;

   if (NULL == (r = Pixel_open_regrid (dest, dest_pixel_area)))
     goto free_and_return;

   if (-1 == find_all_pixel_overlaps (r, files, num_files, lonlat_grp))
     goto free_and_return;

   status = 0;
free_and_return:
   Pixel_list_free (dest_pixel_area);
   if (status)
     {
        Pixel_close_regrid (r);
        r = NULL;
     }

   return r;
}

void Regrid_close (Pixel_Regrid_Type *r)
{
   Pixel_close_regrid (r);
}
