#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <proj_api.h>
#include <tell.h>

#include "proj.h"

/* The Albers equal-area conic projection preserves areas
 * but in general, polygon edges do not project into
 * straight lines.  For "sufficiently small" polygons,
 * the polygon area computed _assuming_ straight edges
 * is "sufficiently accurate".
 */
int proj_longlat_to_albers (double *lon, double *lat, int n)
{
   projPJ albers = NULL, longlat = NULL;
   /* spatialreference.org USA Contiguous Albers Equal Area Conic */
   const char ctl_albers[] =
     "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
   /* spatialreference.org EPSG Projection 4326 - WGS 84  */
   const char ctl_longlat[] =
     "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
   int i, status = -1;

   if ((NULL == (albers = pj_init_plus (ctl_albers)))
       || (NULL == (longlat = pj_init_plus (ctl_longlat))))
     {
        Tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus failed",
                     __func__);
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

Pixel_List_Type *
proj_pixel_list (const Pixel_Grid_Param_Type *mesh, Proj_cvt_type *cvt)
{
   Pixel_List_Type *pixel_list = NULL;
   int num_pixels = mesh->nx * mesh->ny;
   int num_pixel_vertices =
     (4 + 2*mesh->num_extra_xpoints + 2*mesh->num_extra_ypoints);
   double *xs = NULL, *ys = NULL;
   double *x, *y;
   int i, status = -1;

   if (-1 == Pixel_grid_arrays (mesh, &xs, &ys))
     return NULL;

   if (cvt != NULL)
     {
        if (-1 == (*cvt)(xs, ys, num_pixels * num_pixel_vertices))
          goto free_and_return;
     }

   if (NULL == (pixel_list = Pixel_list_new (num_pixels, num_pixel_vertices)))
     goto free_and_return;

   x = xs;
   y = ys;

   for (i = 0; i < num_pixels; i++)
     {
        if (-1 == Pixel_list_set_vertices (pixel_list, i, num_pixel_vertices, x, y))
          goto free_and_return;
        x += num_pixel_vertices;
        y += num_pixel_vertices;
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
