#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>

#include "config.h"
#include "map.h"

#define ELEVATION_TYPE_PRIVATE_DATA \
      Map_Type *elevation;
#include "elevation.h"

static int elevation_lookup (const Elevation_Type *et, unsigned int num,
                             const double *lon, const double *lat,
                             short *elevation)
{
   Map_Type *map = et->elevation;
   return map_lookup_short (map, num, lon, lat, elevation);
}

static int init_tiles (Elevation_Type *et, config_t *cfg)
{
   config_setting_t *s, *sub;
   const char *map_layout_str;
   int i, num_tiles, map_layout;

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

   if (NULL == (sub = config_setting_get_member (s, "tiles")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing terrain height tiles in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   num_tiles = config_setting_length (sub);
   if (num_tiles <= 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: no terrain height files?", __func__);
        return -1;
     }

   et->elevation = map_new (num_tiles, MAP_PIXEL_SHORT, map_layout);
   if (NULL == et->elevation)
     return -1;

   for (i = 0; i < num_tiles; i++)
     {
        const char *tile_file = config_setting_get_string_elem (sub, i);
        if (0 != map_add_tile (et->elevation, i, tile_file))
          return -1;
     }

   return 0;
}

static void free_elevation_type (Elevation_Type *et)
{
   if (et == NULL)
     return;
   map_free (et->elevation);
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

   et->et_delete = free_elevation_type;
   et->et_lookup = elevation_lookup;

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
