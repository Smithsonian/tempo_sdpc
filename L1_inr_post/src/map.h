#ifndef __TEMPO_MAP_INCLUDE__
#define __TEMPO_MAP_INCLUDE__ 1

/** @file map.h
 *  @brief Provide an interface to a collection of non-overlapping map tiles,
 *   supporting coordinate-based lookups and rectangular bounding-box subsets.
 */

#include <poly.h>

typedef struct Map_Type Map_Type;
typedef double Map_Coord_Type;

enum
{
   MAP_PIXEL_SHORT = 0,
   MAP_PIXEL_UBYTE = 1
};

enum
{
   MAP_LAYOUT_GRID = 0,
   MAP_LAYOUT_LATITUDE_BANDS = 1
};

/** Allocate a new Map_Type object
 * @param[in]  num_tiles   Number of non-overlapping map tiles
 * @param[in]  pixel_type  Pixel data type (e.g. MAP_PIXEL_SHORT | MAP_PIXEL_UBYTE)
 * @param[in]  layout_type  Map tile arrangement (e.g. MAP_LAYOUT_GRID | MAP_LAYOUT_LATITUDE_BANDS)
 * @return pointer to initialized Map_Type on success, NULL on error
 */
extern Map_Type *map_new (int num_tiles, int pixel_type, int layout_type);

/** Add a map tile to a Map_Type object
 * @param[in] map  Map_Type object created by \c map_new
 * @param[in] i    index of the map tile to be added
 * @param[in] file path to the \c geotiff-formatted file containing the map tile
 * @return 0 on success, -1 on error
 */
extern int map_add_tile (Map_Type *map, int i, const char *file);

/** Free resources associated with a Map_Type object
 * @param[in]  map   Map_Type object created by \c map_new
 */
extern void map_free (Map_Type *map);

/** Perform coordinate-based nearest-neighbor lookup of short data type
 * @param[in]  map   Map_Type object created by \c map_new
 * @param[in]  num   Number of points to look up
 * @param[in]  lon   Pointer to longitude coordinates
 * @param[in]  lat   Pointer to longitude coordinates
 * @param[out] values  nearest-neighbor map pixel values
 * @return 0 on success, -1 on error
 */
extern int map_lookup_short (Map_Type *map, unsigned int num,
                             const Map_Coord_Type *lon, const Map_Coord_Type *lat,
                             short *values);

/** Perform coordinate-based nearest-neighbor lookup of unsigned char data type
 * @param[in]  map   Map_Type object created by \c map_new
 * @param[in]  num   Number of points to look up
 * @param[in]  lon   Pointer to longitude coordinates
 * @param[in]  lat   Pointer to longitude coordinates
 * @param[out] values  nearest-neighbor map pixel values
 * @return 0 on success, -1 on error
 */
extern int map_lookup_ubyte (Map_Type *map, unsigned int num,
                             const Map_Coord_Type *lon, const Map_Coord_Type *lat,
                             unsigned char *values);

typedef struct
{
   Polygon_Type **poly;   /*!< corners of each tile pixel that overlaps the bounding box */
   double *values;        /*!< value of each tile pixel that overlaps the bounding box */
   int num_pixels;        /*!< number of tile pixels that overlap the bounding box */

#ifdef MAP_SUBSET_PRIVATE_DATA
   MAP_SUBSET_PRIVATE_DATA
#endif
}
Map_Subset_Type;

/** Free resources associated with a Map_Subset_Type
 * @param[in] s   pointer to Map_Subset_Type created by \c map_subset
 */
extern void map_free_subset (Map_Subset_Type *s);

/** Create a Map_Subset_Type for a specified (lon,lat) bounding box
 * @param[in]  map   pointer to a fully initialized Map_Type object
 * @param[in]  lon_bbox  Bounding box corner longitude coordinates [NW, SE]
 * @param[in]  lat_bbox  Bounding box corner latitude coordinates [NW, SE]
 * @param[in]  psubset  Optional pointer to a Map_Subset_Type pointer.
 *                      When (psubset == NULL), a new Map_Subset_Type object
 *                      is allocated.  When (psubset != NULL), the storage
 *                      associated with psubset is re-used.
 * @return Map_Subset_Type pointer on success, NULL on error.
 */
extern Map_Subset_Type *map_subset (Map_Type *map,
                                    const Map_Coord_Type lon_bbox[2],
                                    const Map_Coord_Type lat_bbox[2],
                                    Map_Subset_Type **psubset);

/** Use a Map_Subset_Type to retrieve pixel values from a different Map_Type
 with identically gridded tiles
 * @param[in]  map   pointer to a fully initialized Map_Type object
 * @param[in]  s     pointer to a Map_Subset_Type object (constructed for a different
 *                   Map_Type object with identically gridded tiles)
 * @param[out] values  value of each tile pixel of map that overlaps the bounding box
 *                     for which s was defined.
 * @return 0 on success, -1 on error.
 */
extern int map_apply_subset (Map_Type *map, Map_Subset_Type *s,
                             double *values);

#endif
