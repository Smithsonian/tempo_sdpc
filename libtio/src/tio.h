#ifndef __TIO_INCLUDE_H__
#define __TIO_INCLUDE_H__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <stddef.h>

/* L1 File format version number */
#define TIO_L1_FORMAT_VERSION  "0.1"

#define TIO_CF_CONVENTION_VERSION   "CF-1.6"
#define TIO_TIME_REFERENCE_STRING   "1970-01-01T00:00:00 UTC"

/* Maximum number of array dimensions */
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

#define TIO_DIM_NAME_CORNER         "corner"
#define TIO_DIM_NAME_COV            "cov"
#define TIO_DIM_NAME_STEP           "mirror_step"
#define TIO_DIM_NAME_CHANNEL        "spectral_channel"
#define TIO_DIM_NAME_XTRACK         "xtrack"
#define TIO_DIM_NAME_XYDET          "xy_det"
#define TIO_DIM_NAME_XYSMA          "xy_sma"
#define TIO_DIM_NAME_XYZ            "xyz"
#define TIO_DIM_NAME_XYZSAT         "xyz_sat"

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
#define TIO_VAR_NAME_INRQF          "inr_quality_flag"
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

enum TIO_INR_Status
{
   TIO_INR_NONE    = 0,
   TIO_INR_INITIAL = 1,
   TIO_INR_FINAL   = 2
};

enum TIO_Processing_Level
{
   TIO_PROC_LEVEL_0,
   TIO_PROC_LEVEL_1A,
   TIO_PROC_LEVEL_1B,
   TIO_PROC_LEVEL_1C,
   TIO_PROC_LEVEL_2,
   TIO_PROC_LEVEL_3
};

typedef struct
{
   char *name;
   size_t num_xtrack;
   size_t num_channels;
}
TIO_Radiance_Group_Type;

typedef struct
{
   int varid;   /**< variable ids assigned by nc_def_var */
   int ndims;   /**< number of dimensions */
   int dimids[TIO_MAX_VAR_DIMS];      /**< dimension ids assigned by nc_def_dim */
   size_t dimlens[TIO_MAX_VAR_DIMS];  /**< dimension sizes */
}
TIO_Var_Info_Type;

/** Get information about a variable using its name.
 * @param  grp    Index of group containing the variable
 * @param  name   Variable name string
 * @param  info   Pointer to the structure that will receive the information.
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_var (int grp, const char *name, TIO_Var_Info_Type *info);

/** Write a section of an N-dimensional variable
 * @param  grp     Index of group containing the variable
 * @param  name    Variable name string
 * @param  start0  Write data to the file starting at this value of the
 *                 slowest varying array index
 * @param  count0  Count of slowest varying array index to write
 * @param  xtype   Type of values to write
 * @param  data    Pointer to the array to be written
 * @return 0 on success, -1 on error
 */
extern int TIO_put_var_section (int grp, const char *name,
                                size_t start0, size_t count0, int xtype,
                                const void *data);

/** Read a section of an N-dimensional variable
 * @param  grp     Index of group containing the variable
 * @param  name    Variable name string
 * @param  start0  Read data from the file starting at this value of the
 *                 slowest varying array index
 * @param  count0  Count of slowest varying array index to read
 * @param  xtype   Type of values to be read
 * @param  data    Pointer to the array that will receive the input
 * @return 0 on success, -1 on error
 */
extern int TIO_get_var_section (int grp, const char *name,
                                size_t start0, size_t count0, int xtype,
                                void *data);

/** Query the type and size of an attribute
 * @param  grp      Index of group containing the attribute
 * @param  varname  Name of the corresponding variable.
 *                  Use varname=NULL to query global attributes.
 * @param  attname  Name of the attribute
 * @param  xtype    Pointer to attribute type
 * @param  len      Pointer to attribute length
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_att (int grp, const char *varname, const char *attname,
                        int *xtype, size_t *len);

/** Write an attribute value
 * @param  grp      Index of group containing the attribute
 * @param  varname  Name of the corresponding variable.
 *                  Use varname=NULL to query global attributes.
 * @param  attname  Name of the attribute
 * @param  xtype    Pointer to attribute type
 * @param  len      Pointer to attribute length
 * @param  att      Pointer to attribute value
 * @return 0 on success, -1 on error
 */
extern int TIO_put_att (int grp, const char *varname, const char *attname,
                        int xtype, size_t len, const void *att);

/** Read an attribute value
 * @param  grp      Index of group containing the attribute
 * @param  varname  Name of the corresponding variable.
 *                  Use varname=NULL to query global attributes.
 * @param  attname  Name of the attribute
 * @param  xtype    Pointer to attribute type
 * @param  att      Pointer to attribute value
 * @return 0 on success, -1 on error
 */
extern int TIO_get_att (int grp, const char *varname, const char *attname,
                        int xtype, void *att);

/** Create a template Level 1 data file
 * @param  ncid          Index returned by nc_create
 * @param  num_steps     Number of mirror steps
 * @param  num_rgrps     Number of radiance groups
 * @param  rgrps         Array of TIO_Radiance_Group_Type structs
 * @return 0 on success, -1 on error
 */
extern int TIO_create_l1_template (int ncid, size_t num_steps, int num_rgrps,
                                   TIO_Radiance_Group_Type *rgrps);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
