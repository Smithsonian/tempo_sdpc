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

#define BUFSIZE 1024
static int _pDiagnostic_Output;

typedef struct
{
   double *lon_bounds;
   double *lat_bounds;
   int *step;
   int *xtrack;
   int num_steps;
   int num_xtrack;
   int num_pixels;
}
Source_Pixel_Vertices_Type;

/* Valid longitude, latitude values */
#define INVALID_LONGITUDE(b) (((b) < -180.0) || (360.0 < (b)) || (0 == isfinite(b)))
#define INVALID_LATITUDE(b)  (((b) <  -90.0) || ( 90.0 < (b)) || (0 == isfinite(b)))

/* The Albers equal-area conic projection preserves areas
 * but in general, polygon edges do not project into
 * straight lines.  For "sufficiently small" polygons,
 * the polygon area computed _assuming_ straight edges
 * is "sufficiently accurate".
 */
static int longlat_to_albers (double *lon, double *lat, int n)
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

static void free_spv_type (Source_Pixel_Vertices_Type *spv)
{
   if (NULL == spv)
     return;
   FREE(spv->lon_bounds);
   FREE(spv->lat_bounds);
   FREE(spv->step);
   FREE(spv->xtrack);
   FREE(spv);
}

static Source_Pixel_Vertices_Type *new_spv_type (int num_steps, int num_xtrack)
{
   Source_Pixel_Vertices_Type *spv = NULL;
   int len_bounds;

   if (NULL == (spv = (Source_Pixel_Vertices_Type *)MALLOC (sizeof *spv)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)spv, 0, sizeof *spv);
   spv->num_steps = num_steps;
   spv->num_xtrack = num_xtrack;

   if ((NULL == (spv->step = (int *) MALLOC (num_steps * sizeof (int))))
       || (NULL == (spv->xtrack = (int *) MALLOC (num_xtrack * sizeof (int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_spv_type (spv);
        return NULL;
     }

   spv->num_pixels = num_steps * num_xtrack;
   len_bounds = 4 * spv->num_pixels * sizeof(double);

   if ((NULL == (spv->lon_bounds = (double *) MALLOC (len_bounds)))
       || (NULL == (spv->lat_bounds = (double *) MALLOC (len_bounds))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_spv_type (spv);
        return NULL;
     }

   return spv;
}

static Source_Pixel_Vertices_Type *
read_pixel_vertices (const char *file, const char *lonlat_grp)
{
   Source_Pixel_Vertices_Type *spv = NULL;
   TIO_Var_Info_Type vi;
   double *lon_bounds, *lat_bounds;
   double nan_value = nan("");
   int ncid, grp, start[3], count[3];
   int num_steps, num_xtrack, num_vertices;
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

   if (NULL == (spv = new_spv_type (num_steps, num_xtrack)))
     goto free_and_return;

   start[0] = 0;
   count[0] = num_steps;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_STEP,
                                  start, count, TIO_INT, spv->step))
     goto free_and_return;

   start[0] = 0;
   count[0] = num_xtrack;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_XTRACK,
                                  start, count, TIO_INT, spv->xtrack))
     goto free_and_return;

   /* read lon/lat bounds arrays */

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = num_steps;
   count[1] = num_xtrack;
   count[2] = 4;

   if ((-1 == TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS,
                                   start, count, TIO_DOUBLE, spv->lon_bounds))
       || (-1 == TIO_get_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS,
                                      start, count, TIO_DOUBLE, spv->lat_bounds)))
     {
        goto free_and_return;
     }

   num_vertices = 4 * num_steps * num_xtrack;
   lon_bounds = spv->lon_bounds;
   lat_bounds = spv->lat_bounds;

   for (i = 0; i < num_vertices; i++)
     {
        double lonb_i = lon_bounds[i];
        double latb_i = lat_bounds[i];
        if (INVALID_LONGITUDE(lonb_i)) lon_bounds[i] = nan_value;
        if (INVALID_LATITUDE(latb_i)) lat_bounds[i] = nan_value;
     }

   status = 0;
free_and_return:
   (void) TIO_close (ncid);
   if (status)
     {
        free_spv_type (spv);
        spv = NULL;
     }

   return spv;
}

/* Pack pixel vertices into pixel list structures */
static Pixel_List_Type *make_pixel_list (Source_Pixel_Vertices_Type *spv, int xtrack_dimlen)
{
   Pixel_List_Type *plt = NULL;
   int num_sides;

   /* Zero-length polygons will indicate lines of sight that
    * have invalid lon-lat coordinates -- usually because
    * they don't intersect the earth. */
   num_sides = 0;

   if ((NULL == (plt = Pixel_list_new (spv->num_pixels, num_sides)))
       || (-1 == Pixel_list_use_src_index (plt)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initializing pixel list", __func__);
        goto free_and_return;
     }

   if (-1 == Pixel_list_pack (plt, spv->lon_bounds, spv->lat_bounds,
                              spv->num_pixels, 4, spv->step, spv->xtrack, spv->num_xtrack,
                              xtrack_dimlen))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: packing pixel list", __func__);
        goto free_and_return;
     }

   return plt;

free_and_return:
   Pixel_list_free (plt);
   return NULL;
}

static void start_diagnostic_output (const char *file)
{
   FILE *fp;
   char buf[BUFSIZE];
   const char *base;
   double lon, lat;
   int n;

   /* Do nothing when not in diagnostic mode */
   if (_pDiagnostic_Output == 0) return;

   /* Look in the current directory for a .p file matching this granule's basename */
   if (NULL == (base = strrchr (file, '/')))
     {
        base = file;
     }
   else base++;

   n = snprintf (buf, BUFSIZE, "%s.p", base);
   if ((n < 0) || (n >= BUFSIZE))
     {
        fprintf (stderr, "%s: buffer size %d exceeded, file basename=%s\n", __func__, BUFSIZE, base);
        return;
     }
   /* When the file does not exist, return silently.
    * We may not want diagnostic output for this granule.*/
   if (NULL == (fp = fopen (buf, "r")))
     return;
   /* Read (lon,lat) in degrees */
   n = fscanf (fp, "%lf %lf", &lon, &lat);
   fclose (fp);
   if (n != 2)
     {
        fprintf (stderr, "%s: error reading diagnostic point from file: %s\n", __func__, buf);
        return;
     }

   /* After the call (lon,lat) is really Albers (x,y) in [meters] */
   if (0 != longlat_to_albers (&lon, &lat, 1))
     return;

   __Pixel_diagnostic_output (1);
   __Pixel_diagnostic_window (lon, lat, 50.0e3, 50.0e3);
}

static void stop_diagnostic_output (void)
{
   if (_pDiagnostic_Output == 0) return;
   __Pixel_diagnostic_output (0);
}

void Regrid_diagnostic_output (int b)
{
   _pDiagnostic_Output = b;
}

static int find_max_index (int ncid, const char *varname, size_t count_s, int *pmax)
{
   int *values = NULL;
   int i, max, start = 0, count = count_s;

   if (NULL == (values = (int *)MALLOC (count * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   if (0 != TIO_get_var_section (ncid, varname, &start, &count, TIO_INT, values))
     {
        FREE(values);
        return -1;
     }

   max = 0;
   for (i = 0; i < count; i++)
     {
        if (values[i] > max)
          max = values[i];
     }
   FREE(values);

   *pmax = max;
   return 0;
}

static int find_target_grid_size (char **files, int num_files,
                                  int *pmax_num_step, int *pmax_num_xtrack)
{
   int i, mx, ncid=0, max_step=0, max_xtrack=0;
   int status = -1;

   for (i = 0; i < num_files; i++)
     {
        TIO_Var_Info_Type vi = {0};

        if (0 != TIO_open (files[i], NC_NOWRITE, &ncid))
          return -1;

        if (-1 == TIO_inq_var (ncid, TEMPO_DIM_XTRACK, &vi))
          goto return_status;
        if (0 != find_max_index (ncid, TEMPO_DIM_XTRACK, vi.dimlens[0], &mx))
          goto return_status;
        if (mx > max_xtrack)
          max_xtrack = mx;

        if (-1 == TIO_inq_var (ncid, TEMPO_DIM_STEP, &vi))
          goto return_status;
        if (0 != find_max_index (ncid, TEMPO_DIM_STEP, vi.dimlens[0], &mx))
          goto return_status;
        if (mx > max_step)
          max_step = mx;

        (void) TIO_close (ncid);
        ncid = 0;
     }

   *pmax_num_step = max_step + 1;
   *pmax_num_xtrack = max_xtrack + 1;

   status = 0;
return_status:
   if ((status != 0) && (ncid != 0))
     {
        TIO_close (ncid);
     }
   return status;
}

static int
find_all_pixel_overlaps (Pixel_Regrid_Type *r, char **files, int num_files,
                         const char *lonlat_grp)
{
   Source_Pixel_Vertices_Type *spv = NULL;
   Pixel_List_Type *src_area = NULL;
   Pixel_List_Type *src_lookup = NULL;
   int i, max_num_step, max_num_xtrack;

   /* Loop over granule files, and accumulate contributions
    * to the pixel overlap array in each destination pixel.
    */

   if (0 != find_target_grid_size (files, num_files, &max_num_step, &max_num_xtrack))
     return -1;
   Pixel_regrid_grow_srcdims (r, max_num_step, max_num_xtrack);

   for (i = 0; i < num_files; i++)
     {
        int num_overlaps;

        free_spv_type (spv);
        if (NULL == (spv = read_pixel_vertices (files[i], lonlat_grp)))
          break;

        /* lookup pixel list uses coordinates that simplify
         * determining which destination pixels overlap each source pixel */
        Pixel_list_free (src_lookup);
        if (NULL == (src_lookup = make_pixel_list (spv, max_num_xtrack)))
          break;

        /* WARNING! coordinate projection is done in place, so after the call,
         * (lon_bounds, lat_bounds) [deg] is really Albers (x,y) [meters] */
        if (-1 == longlat_to_albers (spv->lon_bounds, spv->lat_bounds, 4*spv->num_pixels))
          break;

        /* area pixel list uses coordinates that are appropriate
         * for computing pixel overlap areas
         */
        Pixel_list_free (src_area);
        if (NULL == (src_area = make_pixel_list (spv, max_num_xtrack)))
          break;

        start_diagnostic_output (files[i]);
        num_overlaps = Pixel_find_overlaps (r, src_area, src_lookup);
        stop_diagnostic_output ();
        if (num_overlaps < 0)
          break;
        else if (num_overlaps == 0)
          {
             Tell_verror (TELL_UNKNOWN_ERROR,
                          "%s: no contribution to target grid from %s",
                          __func__, files[i]);
          }
     }

   free_spv_type (spv);
   Pixel_list_free (src_area);
   Pixel_list_free (src_lookup);

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
   int num_pixel_vertices =
     (4 + 2*dest->num_extra_xpoints + 2*dest->num_extra_ypoints);
   double *xs = NULL, *ys = NULL;
   double *x, *y;
   int i, status = -1;

   if (-1 == Pixel_grid_arrays (dest, &xs, &ys))
     return NULL;

   if (-1 == longlat_to_albers (xs, ys, num_pixels * num_pixel_vertices))
     goto free_and_return;

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
