/** @file snow.c
 *  @brief Load snow and ice cover data
 *         Perform area-weighted regridding
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <tiffio.h>
#include <geotiff.h>
#include <xtiffio.h>

#include <proj_api.h>
#include <tell.h>
#include <tio.h>
#include <poly.h>

#include "config.h"

#define SNOW_TYPE_PRIVATE_DATA \
   unsigned int num_rows; \
   unsigned int num_cols; \
   unsigned char *pixels;
#include "snow.h"

typedef struct
{
   projPJ albers;
   projPJ longlat;
   projPJ polar_stereo;
}
Transform_Type;

/* IMS images have a fixed size, using polar stereographic projection */
static struct IMS_Image_Definition_Type
{
   double origin_x;     /* distances in meters */
   double origin_y;
   double pixel_dx;
   double pixel_dy;
   unsigned int num_cols;  /* along x */
   unsigned int num_rows;  /* along y */
}
IMS_Image_Def =
{
   .origin_x = -12288000.00,
   .origin_y = +12288000.00,
   .pixel_dx = +1000.0,
   .pixel_dy = -1000.0,
   .num_cols = 24576,
   .num_rows = 24576
};

static int pstereo_to_image (int n, const double *x, const double *y,
                             double *col, double *row)
{
   int i;

   for (i = 0; i < n; i++)
     {
        col[i] = (x[i] - IMS_Image_Def.origin_x) / IMS_Image_Def.pixel_dx;
        row[i] = (y[i] - IMS_Image_Def.origin_y) / IMS_Image_Def.pixel_dy;
     }

   return 0;
}

static int image_to_pstereo (int n, double *x, double *y,
                             const double *col, const double *row)
{
   int i;

   for (i = 0; i < n; i++)
     {
        x[i] = IMS_Image_Def.origin_x + col[i] * IMS_Image_Def.pixel_dx;
        y[i] = IMS_Image_Def.origin_y + row[i] * IMS_Image_Def.pixel_dy;
     }

   return 0;
}

static void tform_close (Transform_Type *tform)
{
   if (tform == NULL)
     return;
   pj_free (tform->albers);
   pj_free (tform->longlat);
   pj_free (tform->polar_stereo);
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
   /* IMS snow/ice data product uses this specific Polar Stereographic projection */
   const char ctl_polar_stereo[] =
     "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs";

   tform->albers = NULL;
   tform->longlat = NULL;
   tform->polar_stereo = NULL;

   if ((NULL == (tform->albers = pj_init_plus (ctl_albers)))
       || (NULL == (tform->longlat = pj_init_plus (ctl_longlat)))
       || (NULL == (tform->polar_stereo = pj_init_plus (ctl_polar_stereo))))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_init_plus failed", __func__);
        tform_close (tform);
        return -1;
     }

   return 0;
}

static int tform_lonlat_to_albers (Transform_Type *tform,
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

static int tform_lonlat_to_pstereo (Transform_Type *tform,
                                    int n, double *lon, double *lat)
{
   int i, status;

   for (i = 0; i < n; i++)
     {
        lon[i] *= DEG_TO_RAD;
        lat[i] *= DEG_TO_RAD;
     }

   status = pj_transform (tform->longlat, tform->polar_stereo, n, 1, lon, lat, NULL);
   if (status)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return -1;
     }

   return 0;
}

static int tform_pstereo_to_lonlat (Transform_Type *tform,
                                    int n, double *x, double *y)
{
   int i, status;

   status = pj_transform (tform->polar_stereo, tform->longlat, n, 1, x, y, NULL);
   if (status)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return -1;
     }

   for (i = 0; i < n; i++)
     {
        x[i] /= DEG_TO_RAD;  /* lon -> deg */
        y[i] /= DEG_TO_RAD;  /* lat -> deg */
     }

   return 0;
}

static void free_snow_type (Snow_Type *sn)
{
   if (sn == NULL)
     return;
   FREE(sn->pixels);
   FREE(sn);
}

static void find_start_count (const double *indices, unsigned int num,
                              unsigned int *start, unsigned int *count)
{
   unsigned int i, min, max;

   min = indices[0];
   max = min;

   for (i = 0; i < num; i++)
     {
        if (indices[i] < min) min = indices[i];
        if (indices[i] > max) max = indices[i];
     }

   *start = min;
   *count = max - min + 1;
}

#define INVALID_COORDINATES(lon,lat) ((0 == isfinite(lon)) || (0 == isfinite(lat)) \
                                      || (fabs(lon) > 360.0) || (fabs(lat) > 90.0))

static int define_bounding_box (Transform_Type *tform,
                                const double *lon_cnr, const double *lat_cnr,
                                unsigned int *col_start, unsigned int *row_start,
                                unsigned int *num_cols, unsigned int *num_rows)
{
   double row[4], col[4];
   double lon[4], lat[4];
   int i;

   for (i = 0; i < 4; i++)
     {
        if (INVALID_COORDINATES(lon_cnr[i], lat_cnr[i]))
          return 1;
        lon[i] = lon_cnr[i];
        lat[i] = lat_cnr[i];
     }

   if (0 != tform_lonlat_to_pstereo (tform, 4, lon, lat))
     return -1;

   if (0 != pstereo_to_image (4, lon, lat, col, row))
     return -1;

   find_start_count (col, 4, col_start, num_cols);
   find_start_count (row, 4, row_start, num_rows);

   return 0;
}

static int get_pixel_lonlat_corners (Transform_Type *tform, unsigned int col, unsigned int row,
                                     double *lon, double *lat)
{
   double dx[] = {-0.5, +0.5, +0.5, -0.5};
   double dy[] = {+0.5, +0.5, -0.5, -0.5}; /* CCW since pixel_dy < 0 */
   double cnr_col[4], cnr_row[4];
   int i;

   for (i = 0; i < 4; i++)
     {
        cnr_col[i] = col + dx[i];
        cnr_row[i] = row + dy[i];
     }

   /* after the call, lon <- x and lat <- y */
   if (0 != image_to_pstereo (4, lon, lat, cnr_col, cnr_row))
     return -1;

   return tform_pstereo_to_lonlat (tform, 4, lon, lat);
}

/* Map polar stereographic grid onto an arbitrary pixel grid.
 * Assume the source data occupies a polar stereographic (lon,lat) grid.
 * Loop over destination grid pixels:
 *    - find source pixels that overlap destination pixel bounding box
 *    - for each overlapping source pixel containing snow or sea ice:
 *       . compute pixel overlap area
 *       . increment overlap area sum
 *    - compute fraction of total pixel area covered by snow or sea ice
 */
static int regrid_polar_stereographic_src (Snow_Type *sn, unsigned int num_pixels,
                                           const double *lon_cnrp, const double *lat_cnrp,
                                           float *snow_ice_fraction)
{
   Transform_Type tform = {0};
   Polygon_Clip_Type *cl = NULL;
   Polygon_Type *target = NULL;
   Polygon_Type *overlap = NULL;
   Polygon_Type *subject = NULL;
   struct bbox_type
     {
        unsigned int row_start, num_rows;
        unsigned int col_start, num_cols;
     }
   bbox;
   unsigned int i, cnr;
   int status = -1;

   if (0 != tform_open (&tform))
     return -1;

   /* lon_cnr, lat_cnr are taken to have length 4*num_pixels,
    * with each pixel's corners packed in groups of 4,
    * CCW around the boundary.
    */

   if (NULL == (cl = Polygon_open_clip ()))
     goto return_status;

   if ((NULL == (target = Polygon_new (4)))
       || (NULL == (subject = Polygon_new (4))))
     goto return_status;

   for (i=0, cnr=0; i < num_pixels; i++, cnr+=4)
     {
        double x0[4], y0[4];
        const double *lon_cnr, *lat_cnr;
        unsigned int pix, num_bb_pixels;
        double sum_area;
        int pixel_status;

        snow_ice_fraction[i] = 0.0;

        /* Next 4 (lon,lat) pixel corners: */
        lon_cnr = lon_cnrp + cnr;
        lat_cnr = lat_cnrp + cnr;

        /* Find bounding box for this pixel in the projected image */
        if ((pixel_status = define_bounding_box (&tform, lon_cnr, lat_cnr,
                                                 &bbox.col_start, &bbox.row_start,
                                                 &bbox.num_cols, &bbox.num_rows)) < 0)
          goto return_status;

        /* skip invalid coordinates */
        if (pixel_status != 0)
          continue;

        /* Make a modifiable copy of the pixel's vertex coordinates */
        memcpy ((char *)x0, (char *)lon_cnr, 4 * sizeof(double));
        memcpy ((char *)y0, (char *)lat_cnr, 4 * sizeof(double));

        /* Use Albers equal-area coordinates to compute snow/ice overlap areas */
        if ((0 != tform_lonlat_to_albers (&tform, 4, x0, y0))
            || (0 != Polygon_set (target, 4, x0, y0)))
          goto return_status;

        num_bb_pixels = bbox.num_rows * bbox.num_cols;
        sum_area = 0.0;

        for (pix = 0; pix < num_bb_pixels; pix++)
          {
             unsigned int row = bbox.row_start + (pix / bbox.num_cols);
             unsigned int col = bbox.col_start + (pix % bbox.num_cols);
             int value = sn->pixels[col + row * sn->num_cols];
             double x[4], y[4];

             /* No contribution when there's no snow or ice.
              * IMS product pixel values are:
              *    0 = outside coverage area
              *    1 = Sea
              *    2 = Land (without snow)
              *    3 = Sea ice
              *    4 = Snow covered land
              */
             if (value < 3)
               continue;

             if (0 != get_pixel_lonlat_corners (&tform, col, row, x, y))
               goto return_status;

             /* Get Albers equal-area coordinates */
             if ((0 != tform_lonlat_to_albers (&tform, 4, x, y))
                 || (0 != Polygon_set (subject, 4, x, y)))
               goto return_status;

             /* Find overlap polygon and compute its area */
             Polygon_free (overlap);
             if (NULL == (overlap = Polygon_clip (cl, target, subject)))
               goto return_status;
             /* zero length indicates no overlap */
             if (0 == Polygon_length (overlap))
               continue;

             sum_area += Polygon_area (overlap);
          }

        if (sum_area > 0)
          {
             double total_area = Polygon_area (target);
             snow_ice_fraction[i] = sum_area / total_area;
          }
     }

   status = 0;
return_status:
   Polygon_free (subject);
   Polygon_free (target);
   Polygon_free (overlap);
   Polygon_close_clip (cl);
   tform_close (&tform);

   return status;
}

static Snow_Type *new_snow_type (void)
{
   Snow_Type *sn = NULL;

   if (NULL == (sn = (Snow_Type *)MALLOC (sizeof *sn)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sn, 0, sizeof *sn);

   sn->sn_delete = free_snow_type;
   sn->sn_regrid = regrid_polar_stereographic_src;

   return sn;
}

static int read_snow_file (Snow_Type *sn, const char *file)
{
   TIFF *tiff = NULL;
   unsigned int row;
   size_t image_size, scan_line_size;
   int pixel_size_bytes, status = -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "opening %s", file);

   if (NULL == (tiff = XTIFFOpen (file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s", __func__, file);
        return -1;
     }

   if (!TIFFGetField (tiff, TIFFTAG_IMAGELENGTH, &sn->num_rows)
       || !TIFFGetField (tiff, TIFFTAG_IMAGEWIDTH, &sn->num_cols))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading TIFF image dimensions",
                    __func__);
        goto close_and_return;
     }

   scan_line_size = TIFFScanlineSize (tiff);
   pixel_size_bytes = scan_line_size / sn->num_cols;

   if (sn->num_cols != scan_line_size)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: type mismatch: expected 8-bit pixels", __func__);
        goto close_and_return;
     }

   image_size = sn->num_rows * sn->num_cols * pixel_size_bytes;

   if (NULL == (sn->pixels = (unsigned char *) MALLOC (image_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   for (row = 0; row < sn->num_rows; row++)
     {
        void *pixel_row = sn->pixels + row * sn->num_cols;
        if (!TIFFReadScanline (tiff, pixel_row, row, 0))
          {
             tell_verror (TELL_IO_READ_ERROR,
                          "%s: reading TIFF image, scan line %d",
                          __func__, row);
             goto close_and_return;
          }
     }

   status = 0;
close_and_return:
   TIFFClose (tiff);

   return status;
}

Snow_Type *snow_init (const char *file)
{
   Snow_Type *sn = NULL;

   if (NULL == (sn = new_snow_type ()))
     return NULL;

   if (0 != read_snow_file (sn, file))
     {
        free_snow_type (sn);
        return NULL;
     }

   return sn;
}
