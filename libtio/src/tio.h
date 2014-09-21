#ifndef __TIO_INCLUDE_H__
#define __TIO_INCLUDE_H__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <stddef.h>

#define TIO_MAX_VAR_DIMS   7

#define TIO_BYTE   NC_BYTE
#define TIO_CHAR   NC_CHAR
#define TIO_UBYTE  NC_UBYTE
#define TIO_SHORT  NC_SHORT
#define TIO_USHORT NC_USHORT
#define TIO_INT    NC_INT
#define TIO_UINT   NC_UINT
#define TIO_INT64  NC_INT64
#define TIO_UINT64 NC_UINT64
#define TIO_FLOAT  NC_FLOAT
#define TIO_DOUBLE NC_DOUBLE

#define TIO_CF_CONVENTION_VERSION   "CF-1.6"
#define TIO_TIME_REFERENCE_STRING   "1970-01-01T00:00:00.0"

#define TIO_DIM_NAME_CORNER         "corner"
#define TIO_DIM_NAME_COV            "cov"
#define TIO_DIM_NAME_STEP           "mirror_step"
#define TIO_DIM_NAME_CHANNEL        "spectral_channel"
#define TIO_DIM_NAME_XTRACK         "xtrack"
#define TIO_DIM_NAME_XYDET          "xy_det"
#define TIO_DIM_NAME_XYSMA          "xy_sma"
#define TIO_DIM_NAME_XYZ            "xyz"
#define TIO_DIM_NAME_XYZSAT         "xyz_sat"

#define TIO_GRP_NAME_BAND1          "band_290_490_nm"
#define TIO_GRP_NAME_BAND2          "band_540_740_nm"
#define TIO_GRP_NAME_GEOMETRY       "geometry"
#define TIO_GRP_NAME_INRINPUT       "inr_input"
#define TIO_GRP_NAME_EPHEMERIS      "ephemeris"
#define TIO_GRP_NAME_MANEUVERS      "maneuvers"
#define TIO_GRP_NAME_TELEMETRY      "telemetry"
#define TIO_GRP_NAME_GYROSCOPE      "gyroscope"
#define TIO_GRP_NAME_MIRROR         "mirror"

#define TIO_VAR_NAME_CLOUDTOPHEIGHT "cloud_top_height"
#define TIO_VAR_NAME_COVARIANCE     "covariance"
#define TIO_VAR_NAME_DELTAV         "delta_v"
#define TIO_VAR_NAME_DQF            "data_quality_flag"
#define TIO_VAR_NAME_EXPTIME        "exptime"
#define TIO_VAR_NAME_GYROBIAS       "gyro_bias"
#define TIO_VAR_NAME_GYRORAW        "gyro_raw"
#define TIO_VAR_NAME_GYROSCALE      "gyro_scale"
#define TIO_VAR_NAME_MOONPOS        "moon_position"
#define TIO_VAR_NAME_PIXELSIZE      "pixel_size"
#define TIO_VAR_NAME_PIXELSCALE     "pixel_scale"
#define TIO_VAR_NAME_RADIANCE       "radiance"
#define TIO_VAR_NAME_SATPOS         "satellite_position"
#define TIO_VAR_NAME_SATVEL         "satellite_velocity"
#define TIO_VAR_NAME_SMADIT         "dit"
#define TIO_VAR_NAME_SRP            "solar_radiation_pressure"
#define TIO_VAR_NAME_SUNPOS         "sun_position"
#define TIO_VAR_NAME_TIME           "time"
#define TIO_VAR_NAME_TIME_EPHEM     "time"
#define TIO_VAR_NAME_TIME_GYRO      "time"
#define TIO_VAR_NAME_TIME_MANEUVER  "time"
#define TIO_VAR_NAME_TIME_SMA       "time"
#define TIO_VAR_NAME_WAVELENGTH     "wavelength"

typedef struct
{
   int varid;
   int ndims;
   int dimids[TIO_MAX_VAR_DIMS];
   size_t dimlens[TIO_MAX_VAR_DIMS];
}
TIO_Var_Info_Type;

extern int TIO_inq_var (int grp, const char *name, TIO_Var_Info_Type *info);
extern int TIO_put_var_section (int grp, const char *name,
                                size_t start0, size_t count0, int xtype,
                                const void *data);
extern int TIO_get_var_section (int grp, const char *name,
                                size_t start0, size_t count0, int xtype,
                                void *data);

extern int TIO_inq_att (int grp, const char *varname, const char *attname,
                        int *xtype, size_t *len);
extern int TIO_put_att (int grp, const char *varname, const char *attname,
                        int xtype, size_t len, const void *att);
extern int TIO_get_att (int grp, const char *varname, const char *attname,
                        int xtype, void *att);

extern int TIO_create_l1b_template (int ncid, size_t num_steps,
                                    size_t num_xtrack, size_t num_channels);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
