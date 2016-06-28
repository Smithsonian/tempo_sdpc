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

#ifndef REALLOC
# define REALLOC realloc
#endif

#ifndef MALLOC
# define MALLOC malloc
#endif

#ifndef FREE
# define FREE free
#endif

typedef struct
{
   Pixel_List_Type *pixel_lookup;
   Pixel_List_Type *pixel_area;
   int num_step;
   int num_xtrack;
}
Source_Pixel_List_Type;

int map_strings (const char **str, int n,
                 int (*do_task)(const char *, void *),
                 void *client_data)
{
   int i;
   for (i = 0; i < n; i++)
     {
        if (-1 == do_task (str[i], client_data))
          return -1;
     }

   return 0;
}

static int count_pixels (const char *file, void *cl)
{
   Source_Pixel_List_Type *src = (Source_Pixel_List_Type *)cl;
   TIO_Var_Info_Type info_step, info_xtrack;
   int status, ncid;

   if (NC_NOERR != (status = nc_open (file, NC_NOWRITE, &ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading (%s)",
                     __func__, file, nc_strerror(status));
        return -1;
     }

   if ((-1 == TIO_inq_var (ncid, TEMPO_DIM_STEP, &info_step))
       || (-1 == TIO_inq_var (ncid, TEMPO_DIM_XTRACK, &info_xtrack)))
     return -1;

   src->num_step += info_step.dimlens[0];
   src->num_xtrack = info_xtrack.dimlens[0];

   (void) nc_close (ncid);

   return 0;
}

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
        Polygon_Type *p = pixel_list->poly[pix];
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
             if (-1 == Polygon_set (p, 4, x, y))
               return -1;
          }
     }

   return 0;
}

static int read_pixel_bounds (int ncid, Source_Pixel_List_Type *src)
{
   TIO_Var_Info_Type vi;
   int start[3], count[3];
   int num_steps, num_pixels, len_bounds;
   double *lon_bounds = NULL, *lat_bounds = NULL;
   double *albers_x_bounds, *albers_y_bounds;
   int *step = NULL;
   int status = -1;

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_STEP, &vi))
     return -1;

   num_steps = vi.dimlens[0];

   if (NULL == (step = (int *) MALLOC (num_steps * sizeof (int))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = 0;
   count[0] = num_steps;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_STEP,
                                  start, count, TIO_INT, step))
     goto free_and_return;

   num_pixels = num_steps * src->num_xtrack;
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
   count[1] = src->num_xtrack;
   count[2] = 4;

   if ((-1 == TIO_get_var_section (ncid, TEMPO_VAR_LONGITUDE_BOUNDS,
                                   start, count, TIO_DOUBLE, lon_bounds))
       || (-1 == TIO_get_var_section (ncid, TEMPO_VAR_LATITUDE_BOUNDS,
                                      start, count, TIO_DOUBLE, lat_bounds)))
     {
        goto free_and_return;
     }

   if (-1 == pack_pixel_list (src->pixel_lookup,
                              lon_bounds, lat_bounds, num_pixels,
                              step, src->num_xtrack))
     {
        goto free_and_return;
     }

   if (-1 == longlat_to_albers (lon_bounds, lat_bounds, 4*num_pixels))
     goto free_and_return;
   /* NOTE: coordinate projection is done in place, so after the call,
    * (lon_bounds, lat_bounds) [deg] is really Albers (x,y) [meters] */
   albers_x_bounds = lon_bounds;
   albers_y_bounds = lat_bounds;

   if (-1 == pack_pixel_list (src->pixel_area,
                              albers_x_bounds, albers_y_bounds, num_pixels,
                              step, src->num_xtrack))
     {
        goto free_and_return;
     }

   status = 0;
free_and_return:
   FREE(lon_bounds);
   FREE(lat_bounds);
   FREE(step);

   return status;
}

static int read_file_pixels (const char *file, void *cl)
{
   Source_Pixel_List_Type *src = (Source_Pixel_List_Type *)cl;
   int ncid, status;

   if (NC_NOERR != (status = nc_open (file, NC_NOWRITE, &ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR,
                     "%s: opening %s for reading (%s)",
                     __func__, file, nc_strerror(status));
        return -1;
     }

   status = read_pixel_bounds (ncid, src);
   (void) nc_close (ncid);

   return status;
}

static void free_source_pixel_list (Source_Pixel_List_Type *src)
{
   if (src == NULL)
     return;
   Pixel_list_free(src->pixel_area);
   Pixel_list_free(src->pixel_lookup);
   FREE(src);
}

static Source_Pixel_List_Type *
read_source_pixel_list (const char **files, int num_files)
{
   Source_Pixel_List_Type *src = NULL;
   int num_src_pixels, num_sides=0, status = -1;

   if (NULL == (src = (Source_Pixel_List_Type *) MALLOC (sizeof *src)))
     return NULL;
   src->num_step = 0;
   src->num_xtrack = 0;
   src->pixel_lookup = NULL;
   src->pixel_area = NULL;

   if (-1 == map_strings (files, num_files, count_pixels, src))
     goto free_and_return;

   /* Zero-length polygons will indicate lines of sight that
    * have invalid lon-lat coordinates -- usually because
    * they don't intersect the earth. */
   num_sides = 0;
   num_src_pixels = src->num_step * src->num_xtrack;

   if ((NULL == (src->pixel_area = Pixel_list_new (num_src_pixels, num_sides)))
       || (NULL == (src->pixel_lookup = Pixel_list_new (num_src_pixels, num_sides))))
     goto free_and_return;

   if (-1 == map_strings (files, num_files, read_file_pixels, src))
     goto free_and_return;

   status = 0;
free_and_return:
   if (status)
     {
        free_source_pixel_list (src);
        src = NULL;
     }

   return src;
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
        Polygon_Type *p = pixel_list->poly[i];
        if (-1 == Polygon_set (p, 4, x, y))
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
Regrid_open (const char **files, int num_files, int *src_dims,
             const Pixel_Grid_Param_Type *dest)
{
   Pixel_Regrid_Type *r = NULL;
   Pixel_List_Type *dest_pixel_area = NULL;
   Source_Pixel_List_Type *src = NULL;
   int status = -1;

   if (NULL == (dest_pixel_area = dest_pixel_area_coords (dest)))
     goto free_and_return;

   if (NULL == (src = read_source_pixel_list (files, num_files)))
     goto free_and_return;

   src_dims[0] = src->num_step;
   src_dims[1] = src->num_xtrack;

   if (NULL == (r = Pixel_open_regrid (src->pixel_area, dest,
                                       src->pixel_lookup, dest_pixel_area)))
     goto free_and_return;

   status = 0;
free_and_return:
   Pixel_list_free (dest_pixel_area);
   free_source_pixel_list (src);
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
