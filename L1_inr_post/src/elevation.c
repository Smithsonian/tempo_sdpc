#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include <libconfig.h>
#include <proj_api.h>
#include <tell.h>
#include <tio.h>
#include <poly.h>

#include "config.h"
#include "map.h"

#define ELEVATION_TYPE_PRIVATE_DATA \
      Map_Type *elevation; \
      Map_Type *elevation_stddev;
#include "elevation.h"

/* FIXME?  or just convert file-specified fill value to DBL_MAX on input? */
#define INVALID_LONGLAT_COORD(c)  (((c) != (c)) || ((c) == DBL_MAX))

static int pixel_bbox (const double lon_cnr[4], const double lat_cnr[4],
                       double lon_bbox[2], double lat_bbox[2])
{
   int k;

   if (INVALID_LONGLAT_COORD(lon_cnr[0])
       || INVALID_LONGLAT_COORD(lat_cnr[0]))
     return -1;

   lon_bbox[0] = lon_cnr[0];
   lon_bbox[1] = lon_cnr[0];
   lat_bbox[0] = lat_cnr[0];
   lat_bbox[1] = lat_cnr[0];

   for (k = 1; k < 4; k++)
     {
        double lon_k = lon_cnr[k];
        double lat_k = lat_cnr[k];

        /* no bounding box when vertices out of range */
        if (INVALID_LONGLAT_COORD(lon_k)
            || INVALID_LONGLAT_COORD(lat_k))
          return -1;

        if (lon_k < lon_bbox[0]) lon_bbox[0] = lon_k;
        else if (lon_k > lon_bbox[1]) lon_bbox[1] = lon_k;

        if (lat_k > lat_bbox[0]) lat_bbox[0] = lat_k;
        else if (lat_k < lat_bbox[1]) lat_bbox[1] = lat_k;
     }

   return 0;
}

typedef struct
{
   projPJ albers;
   projPJ longlat;
}
Transform_Type;

static void tform_close (Transform_Type *tform)
{
   if (tform == NULL)
     return;
   pj_free (tform->albers);
   pj_free (tform->longlat);
}

/* The Albers equal-area conic projection preserves areas
 * but in general, polygon edges do not project into
 * straight lines.  For "sufficiently small" polygons,
 * the polygon area computed _assuming_ straight edges
 * is "sufficiently accurate".
 */
static int tform_open (Transform_Type *tform)
{
   /* USA Contiguous Albers Equal Area Conic */
   const char ctl_albers[] =
     "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
   /* spatialreference.org EPSG Projection 4326 - WGS 84  */
   const char ctl_longlat[] =
     "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

   tform->albers = NULL;
   tform->longlat = NULL;

   if ((NULL == (tform->albers = pj_init_plus (ctl_albers)))
       || (NULL == (tform->longlat = pj_init_plus (ctl_longlat))))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_init_plus failed", __func__);
        tform_close (tform);
        return -1;
     }

   return 0;
}

static int tform_apply (Transform_Type *tform,
                        int n, double *lon, double *lat)
{
   int i, status;

   for (i = 0; i < n; i++)
     {
        lon[i] *= DEG_TO_RAD;
        lat[i] *= DEG_TO_RAD;
     }

   status = pj_transform (tform->longlat, tform->albers, n, 1, lon, lat, NULL);
   if (status)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return -1;
     }

   return 0;
}

typedef struct
{
   double *values;
   int num_alloc;
}
Array_Type;

static void array_free (Array_Type *a)
{
   if (a == NULL)
     return;
   FREE(a->values);
}

static int array_alloc (Array_Type *a, int num)
{
   if (NULL == (a->values = (double *) MALLOC (num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   a->num_alloc = num;

   return 0;
}

static int array_realloc (Array_Type *a, int num)
{
   double *v;
   if (NULL == (v = (double *) REALLOC (a->values, num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }

   a->values = v;
   a->num_alloc = num;

   return 0;
}

/* Map rectangular tile grid onto an arbitrary pixel grid.
 * Assume the source tiles occupy a rectangular (lon,lat) grid.
 * Loop over destination grid pixels:
 *    - find source tiles/pixels that overlap destination pixel bounding box
 *    - for each overlapping source pixel:
 *       . compute pixel overlap area
 *       . increment destination pixel area-weighted value sums
 *    - compute area-weighted value(s)
 */
static int regrid_rect_src (const Elevation_Type *et, unsigned int num_pixels,
                            const double *lon_cnrp, const double *lat_cnrp,
                            short *elev, short *elev_stddev)
{
   Map_Type *map_elev = et->elevation;
   Map_Type *map_stddev = et->elevation_stddev;
   Map_Subset_Type *s = NULL;
   Polygon_Clip_Type *cl = NULL;
   Polygon_Type *clip = NULL;
   Polygon_Type *overlap = NULL;
   Array_Type stddev = {NULL,0};
   Transform_Type tform;
   unsigned int i, cnr;
   int status = -1;

   if (0 != tform_open (&tform))
     return -1;

   /* lon_cnr, lat_cnr are taken to have length 4*num_pixels,
    * with each pixel's corners packed in groups of 4,
    * CCW around the boundary.
    */

   if (NULL == (cl = Polygon_open_clip ()))
     return -1;

   if (NULL == (clip = Polygon_new (4)))
     goto return_status;

   /* Allocation will expand as needed */
   if (0 != array_alloc (&stddev, 128))
     goto return_status;

   for (i=0, cnr=0; i < num_pixels; i++, cnr+=4)
     {
        double lon_bbox[2], lat_bbox[2], x0[4], y0[4];
        double sum_area_elev, sum_area_stddev, sum_area;
        const double *lon_cnr, *lat_cnr;
        int k;

        /* default elevation to the fill value */
        elev[i] = TIO_FILL_SHORT;

        /* Next 4 (lon,lat) pixel corners: */
        lon_cnr = lon_cnrp + cnr;
        lat_cnr = lat_cnrp + cnr;

        /* Construct a rectangular bounding box for the pixel */
        if (0 != pixel_bbox (lon_cnr, lat_cnr, lon_bbox, lat_bbox))
          continue;

        /* Find elevation database pixels that overlap the bounding box. */
        if (NULL == map_subset (map_elev, lon_bbox, lat_bbox, &s))
          goto return_status;

        /* Nothing more to do when there's no elevation data */
        if (s->num_pixels == 0)
          continue;

        /* Collect elevation standard deviation values */
        if (s->num_pixels > stddev.num_alloc)
          {
             if (0 != array_realloc (&stddev, s->num_pixels))
               goto return_status;
          }
        if (0 != map_apply_subset (map_stddev, s, stddev.values))
          goto return_status;

        /* Make a modifiable copy of the pixel's vertex coordinates */
        memcpy ((char *)x0, (char *)lon_cnr, 4 * sizeof(double));
        memcpy ((char *)y0, (char *)lat_cnr, 4 * sizeof(double));

        /* Use Albers equal-area coordinates to compute overlap areas */
        if ((0 != tform_apply (&tform, 4, x0, y0))
            || (0 != Polygon_set (clip, 4, x0, y0)))
          goto return_status;

        /* Accumulate sums weighted by overlap area */
        sum_area = 0.0;
        sum_area_elev = 0.0;
        sum_area_stddev = 0.0;

        for (k = 0; k < s->num_pixels; k++)
          {
             Polygon_Type *subject = s->poly[k];
             double x[4], y[4], overlap_area;

             /* No contribution when elevation data is missing */
             if (SHRT_MIN == (short) s->values[k])
               continue;

             /* Use Albers equal-area coordinates for overlap calculation */
             if ((0 != Polygon_get (subject, 4, x, y))
                 || (0 != tform_apply (&tform, 4, x, y))
                 || (0 != Polygon_set (subject, 4, x, y)))
               goto return_status;

             /* Find overlap polygon and compute its area */
             Polygon_free (overlap);
             if (NULL == (overlap = Polygon_clip (cl, clip, subject)))
               goto return_status;
             /* zero length indicates no overlap */
             if (0 == Polygon_length (overlap))
               continue;

             overlap_area = Polygon_area (overlap);

             sum_area += overlap_area;
             sum_area_elev += overlap_area * s->values[k];
             sum_area_stddev += overlap_area * stddev.values[k];
          }

        if (sum_area > 0)
          {
             elev[i] = sum_area_elev / sum_area;
             elev_stddev[i] = sum_area_stddev / sum_area;
          }
     }

   status = 0;
return_status:
   map_free_subset (s);
   Polygon_free (clip);
   Polygon_free (overlap);
   Polygon_close_clip (cl);
   tform_close (&tform);
   array_free (&stddev);

   return status;
}

static int lookup_nearest_neighbor (const Elevation_Type *et, unsigned int num,
                                    const double *lon, const double *lat,
                                    short *elev, short *elev_stddev)
{
   if (0 != map_lookup_short (et->elevation, num, lon, lat, elev))
     return -1;
   if (0 != map_lookup_short (et->elevation_stddev, num, lon, lat, elev_stddev))
     return -1;
   return 0;
}

static Map_Type *load_map_tiles (config_setting_t *s,
                                 const char *tile_list, int map_layout)
{
   config_setting_t *sub;
   Map_Type *map = NULL;
   int i, num_tiles;

   if (NULL == (sub = config_setting_get_member (s, tile_list)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing map tile list '%s' in param file",
                     __func__, tile_list);
        return NULL;
     }

   if ((num_tiles = config_setting_length (sub)) <= 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: no map tiles?", __func__);
        return NULL;
     }

   if (NULL == (map = map_new (num_tiles, MAP_PIXEL_SHORT, map_layout)))
     return NULL;

   for (i = 0; i < num_tiles; i++)
     {
        const char *tile_file = config_setting_get_string_elem (sub, i);
        if (0 != map_add_tile (map, i, tile_file))
          {
             map_free (map);
             return NULL;
          }
     }

   return map;
}

static int init_tiles (Elevation_Type *et, config_t *cfg)
{
   config_setting_t *s;
   const char *map_layout_str;
   int map_layout;

   if (NULL == (s = config_lookup (cfg, "elevation")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing elevation data in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "tile_layout", &map_layout_str))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading elevation data tile_layout: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (0 == strcmp (map_layout_str, "latitude_bands"))
     map_layout = MAP_LAYOUT_LATITUDE_BANDS;
   else if (0 == strcmp (map_layout_str, "grid"))
     map_layout = MAP_LAYOUT_GRID;
   else
     {
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR, "%s: unsupported tile layout '%s'",
                     __func__, map_layout_str);
        return -1;
     }

   if (NULL == (et->elevation = load_map_tiles (s, "tiles_mean", map_layout)))
     return -1;

   if (NULL == (et->elevation_stddev = load_map_tiles (s, "tiles_stddev", map_layout)))
     return -1;

   return 0;
}

static void free_elevation_type (Elevation_Type *et)
{
   if (et == NULL)
     return;
   map_free (et->elevation);
   map_free (et->elevation_stddev);
   FREE(et);
}

static Elevation_Type *new_elevation_type (void)
{
   Elevation_Type *et = NULL;

   if (NULL == (et = (Elevation_Type *)MALLOC (sizeof *et)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)et, 0, sizeof *et);

   et->et_delete = free_elevation_type;
   et->et_lookup = lookup_nearest_neighbor;
   et->et_regrid = regrid_rect_src;

   return et;
}

Elevation_Type *elevation_init (config_t *cfg)
{
   Elevation_Type *et = NULL;

   if (NULL == (et = new_elevation_type ()))
     return NULL;

   if (0 != init_tiles (et, cfg))
     {
        free_elevation_type (et);
        return NULL;
     }

   return et;
}
