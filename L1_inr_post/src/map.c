#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>

#include <tiffio.h>
#include <geotiff/geotiff.h>
#include <geotiff/xtiffio.h>
#include <geotiff/geo_normalize.h>
#include <geotiff/geo_simpletags.h>
#include <geotiff/geovalues.h>

#include <libconfig.h>
#include <tell.h>

#include "config.h"
#include "map.h"

typedef struct
{
   TIFF *tiff;
   GTIF *gtiff;
   GTIFDefn defn;
   Map_Coord_Type lon_min;
   Map_Coord_Type lon_max;
   Map_Coord_Type lat_min;
   Map_Coord_Type lat_max;
   unsigned int num_rows;
   unsigned int num_cols;
   int pixel_type;
   union
     {
        void *v;
        short *s;
        unsigned char *uc;
     }
   pixels;
}
Tile_Type;

struct Map_Type
{
   Tile_Type *tiles;
   int (*find_tile)(Map_Type *, Map_Coord_Type, Map_Coord_Type);
   int num_tiles;
   int pixel_type;
   int layout_type;
   int last_tile_used;
};

static void tiff_warning_handler (const char *module, const char *fmt, va_list ap)
{
#define BUFSIZE 1024
   char buf[BUFSIZE];
   char **pw;
   char *pointless_warnings[] =
     {
        "Unknown field with tag 42113", /* GDAL_NODATA */
        NULL
     };

   if (vsnprintf (buf, sizeof(buf), fmt, ap) < 0)
     return;

   for (pw = pointless_warnings; *pw != NULL; pw++)
     {
        if (NULL != strstr (buf, *pw))
          return;
     }

   tell_vlog (TELL_MSGTYPE_WARN, 0, "%s%sWarning, %s\n",
              module ? module : "",
              module ? ": " : "", buf);
}

static TIFF *open_tiff (const char *file, const char *mode)
{
   TIFFErrorHandler warn;
   TIFF *tiff;

   /* Use a custom warning handler to filter pointless warnings
    * when the file is opened.
    */

   warn = TIFFSetWarningHandler (tiff_warning_handler);
   tiff = XTIFFOpen (file, mode);
   TIFFSetWarningHandler (warn);

   return tiff;
}

static int find_tile_grid (Map_Type *map, Map_Coord_Type lon, Map_Coord_Type lat)
{
   Tile_Type *tile;
   int i;

   /* If this is noticably inefficient, we may need to
    * implement a more clever search.
    */

   if ((0 <= map->last_tile_used)
       && (map->last_tile_used < map->num_tiles))
     {
        tile = &map->tiles[map->last_tile_used];
        if (((tile->lon_min <= lon) && (lon < tile->lon_max))
            && ((tile->lat_min <= lat) && (lat < tile->lat_max)))
          {
             return map->last_tile_used;
          }
     }

   for (i = 0; i < map->num_tiles; i++)
     {
        tile = &map->tiles[i];
        if (((tile->lon_min <= lon) && (lon < tile->lon_max))
            && ((tile->lat_min <= lat) && (lat < tile->lat_max)))
          {
             map->last_tile_used = i;
             return i;
          }
     }

   map->last_tile_used = -1;
   return -1;
}

static int find_tile_latitude_band (Map_Type *map, Map_Coord_Type lon, Map_Coord_Type lat)
{
   Tile_Type *tile;
   int i;

   /* Some latitude band tiles may contain pixels that don't map
    * onto the earth, so the corresponding longitude range may
    * be set to "invalid longitude".  For this reason, we don't
    * use the longitude value when selecting a tile.  If the
    * lon value isn't on the selected tile, then that will be
    * determined later.
    */
   (void) lon;

   /* If this is noticably inefficient, we may need to
    * implement a more clever search.
    */

   if ((0 <= map->last_tile_used) &&
       (map->last_tile_used < map->num_tiles))
     {
        tile = &map->tiles[map->last_tile_used];
        if ((tile->lat_min <= lat) && (lat < tile->lat_max))
          {
             return map->last_tile_used;
          }
     }

   for (i = 0; i < map->num_tiles; i++)
     {
        tile = &map->tiles[i];
        if ((tile->lat_min <= lat) && (lat < tile->lat_max))
          {
             map->last_tile_used = i;
             return i;
          }
     }

   map->last_tile_used = -1;
   return -1;
}

static void free_tile_members (Tile_Type *tile)
{
   if (tile == NULL)
     return;
   FREE(tile->pixels.v);
   GTIFFree (tile->gtiff);
   XTIFFClose (tile->tiff); /* XTIFFClose must follow GTIFFree */
}

static void free_tile_array (Tile_Type *tiles, int num_tiles)
{
   if (tiles == NULL)
     return;

   while (num_tiles-- > 0)
     {
        free_tile_members (&tiles[num_tiles]);
     }

   FREE(tiles);
}

static Tile_Type *new_tile_array (int num_tiles)
{
   Tile_Type *tiles = NULL;
   size_t array_size = num_tiles * sizeof (*tiles);

   if (NULL == (tiles = (Tile_Type *) MALLOC (array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)tiles, 0, array_size);

   return tiles;
}

void map_free (Map_Type *map)
{
   if (map == NULL)
     return;
   free_tile_array (map->tiles, map->num_tiles);
   FREE(map);
}

Map_Type *map_new (int num_tiles, int pixel_type, int layout_type)
{
   Map_Type *map = NULL;

   if (NULL == (map = (Map_Type *)MALLOC (sizeof *map)))
     return NULL;
   memset ((char *)map, 0, sizeof (*map));

   map->num_tiles = num_tiles;
   map->pixel_type = pixel_type;
   map->layout_type = layout_type;
   map->last_tile_used = -1;

   switch (layout_type)
     {
      case MAP_LAYOUT_GRID:
        map->find_tile = find_tile_grid;
        break;
      case MAP_LAYOUT_LATITUDE_BANDS:
        map->find_tile = find_tile_latitude_band;
        break;
      default:
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: invalid map layout_type=%d",
                     __func__, layout_type);
        map_free (map);
        return NULL;
     }

   if (NULL == (map->tiles = new_tile_array (num_tiles)))
     {
        map_free (map);
        return NULL;
     }

   return map;
}

static int read_tile (Tile_Type *tile, int pixel_type)
{
   unsigned int row;
   size_t image_size, scan_line_size;
   int pixel_size_bytes;

   if (!TIFFGetField (tile->tiff, TIFFTAG_IMAGELENGTH, &tile->num_rows)
       || !TIFFGetField (tile->tiff, TIFFTAG_IMAGEWIDTH, &tile->num_cols))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading TIFF image dimensions",
                    __func__);
        return -1;
     }

   scan_line_size = TIFFScanlineSize (tile->tiff);

   pixel_size_bytes = scan_line_size / tile->num_cols;

   switch (pixel_size_bytes)
     {
      case 1:
        tile->pixel_type = MAP_PIXEL_UBYTE;
        break;
      case 2:
        tile->pixel_type = MAP_PIXEL_SHORT;
        break;
      default:
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR,
                     "%s: unsupported pixel type: %d bytes",
                     __func__, pixel_size_bytes);
        return -1;
        break;
     }

   if (tile->pixel_type != pixel_type)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: type mismatch: tile pixel_type=%d (expected %d)",
                     __func__, tile->pixel_type, pixel_type);
        return -1;
     }

   if (NULL == (tile->gtiff = GTIFNew (tile->tiff)))
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: GTIFNew failed", __func__);
        return -1;
     }

   if (0 == GTIFGetDefn (tile->gtiff, &tile->defn))
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: GTIFGetDefn failed", __func__);
        return -1;
     }

   image_size = tile->num_rows * tile->num_cols * pixel_size_bytes;

   if (NULL == (tile->pixels.v = MALLOC (image_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (row = 0; row < tile->num_rows; row++)
     {
        void *pixel_row;
        switch (tile->pixel_type)
          {
           case MAP_PIXEL_SHORT:
             pixel_row = tile->pixels.s + row * tile->num_cols;
             break;
           case MAP_PIXEL_UBYTE:
             pixel_row = tile->pixels.uc + row * tile->num_cols;
             break;
          }
        if (!TIFFReadScanline (tile->tiff, pixel_row, row, 0))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading TIFF image, scan line %d",
                          __func__, row);
             return -1;
          }
     }

   return 0;
}

static int image_to_lonlat (Tile_Type *tile,
                            unsigned int col, unsigned int row,
                            Map_Coord_Type *lon, Map_Coord_Type *lat)
{
   GTIFDefn *defn = &tile->defn;
   double x, y;

   if (row > tile->num_rows || col > tile->num_cols)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: col=%d row=%d is outside the tile boundary num_cols=%d num_rows=%d",
                     __func__, col, row, tile->num_cols, tile->num_rows);
        return -1;
     }

   x = col;
   y = row;
   if (!GTIFImageToPCS(tile->gtiff, &x, &y))
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: computing projected coordinates",
                     __func__);
        return -1;
     }

   if (defn->Model != ModelTypeGeographic)
     {
        if (!GTIFProj4ToLatLong (defn, 1, &x, &y))
          {
             tell_verror (TELL_UNKNOWN_ERROR, "%s: computing (lon,lat) coordinates from projected coordinates",
                          __func__);
             return -1;
          }
     }

   *lon = x;
   *lat = y;
   return 0;
}

static int lonlat_to_image (Tile_Type *tile,
                            unsigned int *col, unsigned int *row,
                            Map_Coord_Type lon, Map_Coord_Type lat)
{
   GTIFDefn *defn = &tile->defn;
   double x, y;

#if 0
   if ((lat < tile->lat_min) || (tile->lat_max < lat)
       || (lon < tile->lon_min)|| (tile->lon_max < lon))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: lon=%10.6f lat=%10.6f is outside the tile boundary Lon=%10.6f:%10.6f Lat=%10.6f:%10.6f",
                     __func__, lon, lat,
                     tile->lon_min, tile->lon_max,
                     tile->lat_min, tile->lat_max);
        return -1;
     }
#endif

   x = lon;
   y = lat;

   if (defn->Model != ModelTypeGeographic)
     {
        if (!GTIFProj4FromLatLong (defn, 1, &x, &y))
          {
             tell_verror (TELL_UNKNOWN_ERROR, "%s: computing projected coordinates from (lon,lat)",
                          __func__);
             return -1;
          }
     }

   if (!GTIFPCSToImage(tile->gtiff, &x, &y))
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: computing image coordinates",
                     __func__);
        return -1;
     }

   *col = x;
   *row = y;

   if (*col > tile->num_cols)
     return -1;
   else if (*col == tile->num_cols)
     *col -= 1;

   if (*row > tile->num_rows)
     return -1;
   else if (*row == tile->num_rows)
     *row -= 1;

   return 0;
}

static int init_tile1 (const char *tile_file, int pixel_type,
                       Tile_Type *tile)
{
   if (NULL == (tile->tiff = open_tiff (tile_file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s", __func__, tile_file);
        return -1;
     }

   if (0 != read_tile (tile, pixel_type))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading %s",
                     __func__, tile_file);
        return -1;
     }

   /* NW corner is (0,0)
    * SE corner is (num_rows, num_cols)
    */
   if ((0 != image_to_lonlat (tile, 0, 0, &tile->lon_min, &tile->lat_max))
       || (0 != image_to_lonlat (tile, tile->num_cols, tile->num_rows,
                                 &tile->lon_max, &tile->lat_min)))
     {
        return -1;
     }

   return 0;
}

#define MAP_LOOKUP_TYPE(type,name,NAME,union_field,missing_value) \
int map_lookup_##name (Map_Type *map, unsigned int num, \
                       const Map_Coord_Type *lon, const Map_Coord_Type *lat, \
                       type *values) \
{ \
   unsigned int i; \
   unsigned int num_off_tiles, num_not_found; \
 \
   if (map->pixel_type != MAP_PIXEL_##NAME) \
     {\
        tell_verror (TELL_INVALID_PARM_ERROR, \
                     "%s: type mismatch, map pixel type = %d, expected %d", \
                     __func__, map->pixel_type, MAP_PIXEL_##NAME); \
        return -1; \
     } \
 \
   num_off_tiles = 0; \
   num_not_found = 0; \
 \
   for (i = 0; i < num; i++) \
     { \
        Tile_Type *tile; \
        unsigned int col, row; \
        int k; \
 \
        if ((k = map->find_tile (map, lon[i], lat[i])) < 0) \
          { \
             values[i] = missing_value; \
             num_off_tiles++; \
             continue; \
          } \
 \
        tile = &map->tiles[k]; \
 \
        if (0 != lonlat_to_image (tile, &col, &row, lon[i], lat[i])) \
          { \
             values[i] = missing_value; \
             num_not_found++; \
             continue; \
          } \
 \
        values[i] = tile->pixels.union_field[col + row * tile->num_cols]; \
     } \
 \
   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: num_not_found=%d num_off_tiles=%d\n", \
            __func__, num_not_found, num_off_tiles); \
 \
   return 0; \
}

MAP_LOOKUP_TYPE(short,short,SHORT,s,SHRT_MIN)
MAP_LOOKUP_TYPE(unsigned char,ubyte,UBYTE,uc,UCHAR_MAX)

int map_add_tile (Map_Type *map, int i, const char *file)
{
   Tile_Type *tile = &map->tiles[i];
   return init_tile1 (file, map->pixel_type, tile);
}

