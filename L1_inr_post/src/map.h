#ifndef __TEMPO_MAP_INCLUDE__
#define __TEMPO_MAP_INCLUDE__ 1

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

extern Map_Type *map_new (int num_tiles, int pixel_type, int layout_type);
extern int map_add_tile (Map_Type *map, int i, const char *file);
extern void map_free (Map_Type *map);

extern int map_lookup_short (Map_Type *map, unsigned int num,
                             const Map_Coord_Type *lon, const Map_Coord_Type *lat,
                             short *values);
extern int map_lookup_ubyte (Map_Type *map, unsigned int num,
                             const Map_Coord_Type *lon, const Map_Coord_Type *lat,
                             unsigned char *values);

typedef struct
{
   Polygon_Type **poly;
   double *values;
   int num_pixels;

#ifdef MAP_SUBSET_PRIVATE_DATA
   MAP_SUBSET_PRIVATE_DATA
#endif
}
Map_Subset_Type;

extern void map_free_subset (Map_Subset_Type *s);

extern Map_Subset_Type *map_subset (Map_Type *map,
                                    const Map_Coord_Type lon_bbox[2],
                                    const Map_Coord_Type lat_bbox[2],
                                    Map_Subset_Type **psubset);

extern int map_apply_subset (Map_Type *map, Map_Subset_Type *s,
                             double *values);

#endif
