#ifndef __TEMPO_MAP_INCLUDE__
#define __TEMPO_MAP_INCLUDE__ 1

typedef struct Map_Type Map_Type;

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
                             const float *lon, const float *lat,
                             short *values);
extern int map_lookup_ubyte (Map_Type *map, unsigned int num,
                             const float *lon, const float *lat,
                             unsigned char *values);

#endif
