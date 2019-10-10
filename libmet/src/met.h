#ifndef __MET_LIB_INTERFACE__
#define __MET_LIB_INTERFACE__ 1

typedef struct
{
   float pressure_surface;         /* pressure interpolated at (lon,lat) surface [hPa] */
   float pressure_tropopause;      /* pressure interpolated at (lon,lat) tropopause [hPa] */
   float *temperature_on_isobar;   /* temperature_on_isobar[i] interpolated at (lon,lat,isobars[i]) [K] */
   float *isobars;                 /* optional input vector of pressures, monotonic decreasing, [hPa] */
   int num_isobars;
}
Met_Value_Type;

/* Low level interface: */

enum
{
   MFT_INTERP_FAIL = -1,
   MFT_INTERP_SUCCESS = 0,
   MFT_INTERP_DOMAIN_ERROR = 1
};

typedef struct Met_File_Type Met_File_Type;

struct Met_File_Type
{
   void (*mft_close)(Met_File_Type *mft);
   int (*mft_interp)(Met_File_Type *mft, float lon, float lat, Met_Value_Type *mvt);

#ifdef MFT_PRIVATE_DATA
   MFT_PRIVATE_DATA
#endif
};

#define MET_READ_PRESSURE_SURFACE       (1<<0)
#define MET_READ_PRESSURE_TROPOPAUSE    (1<<1)
#define MET_READ_TEMPERATURE_ON_ISOBARS (1<<2)

extern Met_File_Type *met_open_file_grib2 (const char *path, int flags);

/* High level interface */

typedef struct Met_List_Type Met_List_Type;

extern Met_List_Type *met_list_new (int flags);
extern void met_list_free (Met_List_Type *met_list);
extern int met_list_add_file (Met_List_Type *met_list, const char *path);
extern int met_list_interp (Met_List_Type *met_list, float lon, float lat, Met_Value_Type *mvt);

#endif
