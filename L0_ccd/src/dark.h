#ifndef __TEMPO_DARK_INCLUDE_H__
#define __TEMPO_DARK_INCLUDE_H__ 1

#include <libconfig.h>
#include "image.h"

typedef struct Dark_Array_Type Dark_Array_Type;

extern void dark_array_free (Dark_Array_Type *da);
extern Dark_Array_Type *dark_array_alloc (int num_darks);
extern int dark_array_elem_init (Dark_Array_Type *da, int i,
                                 Image_Type *img, double sdc, double fp_temp,
                                 double exposure_time);
extern int dark_array_write (const Dark_Array_Type *dark_array,
                             const char *file);

typedef struct Dark_Table_Type Dark_Table_Type;
typedef struct Dark_Config_Type Dark_Config_Type;

enum
{
   DARK_TABLE_ORDERED_BY_TEMP = 1,
   DARK_TABLE_ORDERED_BY_SDC = 2,
   DARK_TABLE_ORDERED_BY_EXPTIME = 3
};

struct Dark_Table_Type
{
   void (*dtt_delete)(Dark_Table_Type *);
   int (*dtt_ordering)(const Dark_Table_Type *);
   int (*dtt_write)(const Dark_Table_Type *, const char *);
   int (*dtt_interp)(const Dark_Table_Type *, double, Image_Type *);
   void (*dtt_domain)(const Dark_Table_Type *, double *, double *);

#ifdef DARK_TABLE_PRIVATE_DATA
   DARK_TABLE_PRIVATE_DATA
#endif
};

extern Dark_Table_Type *dark_table_read (const char *file);

extern void dark_table_config_free (Dark_Config_Type *dcfg);
extern Dark_Config_Type *dark_table_config (config_t *cfg, int is_linearity);
extern Dark_Table_Type *dark_table_create (const Dark_Config_Type *dcfg,
                                           const Dark_Array_Type *);
#endif
