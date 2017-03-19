#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>

#include "config.h"
#include "map.h"

#define LAND_COVER_PRIVATE_DATA \
      Map_Type *land_cover_type1; \
      Map_Type *land_cover_typeqc;
#include "land_cover.h"

static int land_cover_lookup_type1 (const Land_Cover_Type *lc,
                                    unsigned int num, const double *lon, const double *lat,
                                    unsigned char *mask)
{
   return map_lookup_ubyte (lc->land_cover_type1, num, lon, lat, mask);
}

static int land_cover_lookup_typeqc (const Land_Cover_Type *lc,
                                     unsigned int num, const double *lon, const double *lat,
                                     unsigned char *mask)
{
   return map_lookup_ubyte (lc->land_cover_typeqc, num, lon, lat, mask);
}

static Map_Type *load_map (config_setting_t *s, const char *name,
                           int pixel_type, int map_layout)
{
   config_setting_t *sub;
   Map_Type *map = NULL;
   int i, num_tiles;

   if (NULL == (sub = config_setting_get_member (s, name)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing config file setting %s", __func__, name);
        return NULL;
     }

   num_tiles = config_setting_length (sub);
   if (num_tiles <= 0)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: empty tile list?", __func__);
        return NULL;
     }

   if (NULL == (map = map_new (num_tiles, pixel_type, map_layout)))
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

static int init_tiles (Land_Cover_Type *lc, config_t *cfg)
{
   config_setting_t *s;
   const char *map_layout_str;
   int map_layout;

   if (NULL == (s = config_lookup (cfg, "land_cover")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing land_cover data in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "tile_layout", &map_layout_str))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading land_cover data tile_layout: %s",
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

   if (NULL == (lc->land_cover_type1 = load_map (s, "tiles_type1", MAP_PIXEL_UBYTE, map_layout)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: loading land cover type 1 tiles in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (lc->land_cover_typeqc = load_map (s, "tiles_typeQC", MAP_PIXEL_UBYTE, map_layout)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: loading land cover type QC tiles in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static void free_land_cover_type (Land_Cover_Type *lc)
{
   if (lc == NULL)
     return;
   map_free (lc->land_cover_type1);
   map_free (lc->land_cover_typeqc);
   FREE(lc);
}

static Land_Cover_Type *new_land_cover_type (void)
{
   Land_Cover_Type *lc = NULL;

   if (NULL == (lc = (Land_Cover_Type *)MALLOC (sizeof *lc)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)lc, 0, sizeof *lc);

   lc->lc_delete = free_land_cover_type;
   lc->lc_lookup_type1 = land_cover_lookup_type1;
   lc->lc_lookup_typeqc = land_cover_lookup_typeqc;

   return lc;
}

Land_Cover_Type *land_cover_init (config_t *cfg)
{
   Land_Cover_Type *lc = NULL;

   if (NULL == (lc = new_land_cover_type ()))
     return NULL;

   if (0 != init_tiles (lc, cfg))
     {
        free_land_cover_type (lc);
        return NULL;
     }

   return lc;
}
