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

/* File dimension names */
#define TIO_DIM_NAME_CORNER         "corner"
#define TIO_DIM_NAME_COV            "cov"
#define TIO_DIM_NAME_STEP           "mirror_step"
#define TIO_DIM_NAME_CHANNEL        "spectral_channel"
#define TIO_DIM_NAME_XTRACK         "xtrack"

/* File group names */
#define TIO_GRP_NAME_GEOMETRY       "geometry"
#define TIO_GRP_NAME_INRINPUT       "inr_input"
#define TIO_GRP_NAME_EPHEMERIS      "ephemeris"
#define TIO_GRP_NAME_MANEUVERS      "maneuvers"
#define TIO_GRP_NAME_TELEMETRY      "telemetry"
#define TIO_GRP_NAME_GYROSCOPE      "gyroscope"
#define TIO_GRP_NAME_MIRROR         "mirror"

/* File variable names */
#define TIO_VAR_NAME_AMR                 "effective_area_to_mass_ratio"
#define TIO_VAR_NAME_CLOUDTOPHEIGHT      "cloud_top_height"
#define TIO_VAR_NAME_DELTAV_X            "delta_v_x"
#define TIO_VAR_NAME_DELTAV_Y            "delta_v_y"
#define TIO_VAR_NAME_DELTAV_Z            "delta_v_z"
#define TIO_VAR_NAME_DQF                 "data_quality_flag"
#define TIO_VAR_NAME_EXPOSURE_TIME       "exposure_time"
#define TIO_VAR_NAME_GRANULE_FLAG        "granule_flag"
#define TIO_VAR_NAME_GYRO_ROLL           "gyro_roll"
#define TIO_VAR_NAME_GYRO_PITCH          "gyro_pitch"
#define TIO_VAR_NAME_GYRO_YAW            "gyro_yaw"
#define TIO_VAR_NAME_INR_COVARIANCE      "inr_covariance"
#define TIO_VAR_NAME_INRQF               "inr_quality_flag"
#define TIO_VAR_NAME_IRRADIANCE          "irradiance"
#define TIO_VAR_NAME_MIRROR_STEP_SIZE    "mirror_step_size"
#define TIO_VAR_NAME_MOON_X              "moon_X"
#define TIO_VAR_NAME_MOON_Y              "moon_Y"
#define TIO_VAR_NAME_MOON_Z              "moon_Z"
#define TIO_VAR_NAME_PIXEL_SCALE_COLUMN  "pixel_scale_column"
#define TIO_VAR_NAME_PIXEL_SCALE_ROW     "pixel_scale_row"
#define TIO_VAR_NAME_PQF                 "pixel_quality_flag"
#define TIO_VAR_NAME_RADIANCE            "radiance"
#define TIO_VAR_NAME_SAT_X               "satellite_X"
#define TIO_VAR_NAME_SAT_Y               "satellite_Y"
#define TIO_VAR_NAME_SAT_Z               "satellite_Z"
#define TIO_VAR_NAME_SAT_VX              "satellite_velocity_X"
#define TIO_VAR_NAME_SAT_VY              "satellite_velocity_Y"
#define TIO_VAR_NAME_SAT_VZ              "satellite_velocity_Z"
#define TIO_VAR_NAME_SMADIT_EAST         "dit_east"
#define TIO_VAR_NAME_SMADIT_NORTH        "dit_north"
#define TIO_VAR_NAME_SUN_X               "sun_X"
#define TIO_VAR_NAME_SUN_Y               "sun_Y"
#define TIO_VAR_NAME_SUN_Z               "sun_Z"
#define TIO_VAR_NAME_TIME                "time"
#define TIO_VAR_NAME_TIME_EPHEM          "ephemeris_time"
#define TIO_VAR_NAME_TIME_GYRO           "gyro_time"
#define TIO_VAR_NAME_TIME_MANEUVER       "maneuver_time"
#define TIO_VAR_NAME_TIME_SMA            "sma_time"
#define TIO_VAR_NAME_WAVELENGTH          "wavelength"

/* fill values */
#define TIO_FILL_BYTE    ((signed char)-127)
#define TIO_FILL_CHAR    ((char)0)
#define TIO_FILL_SHORT   ((short)-32767)
#define TIO_FILL_INT     (-2147483647L)
#define TIO_FILL_FLOAT   (9.9692099683868690e+36f) /* near 15 * 2^119 */
#define TIO_FILL_DOUBLE  (9.9692099683868690e+36)
#define TIO_FILL_UBYTE   (255)
#define TIO_FILL_USHORT  (65535)
#define TIO_FILL_UINT    (4294967295U)
#define TIO_FILL_INT64   ((long long)-9223372036854775806LL)
#define TIO_FILL_UINT64  ((unsigned long long)18446744073709551614ULL)
#define TIO_FILL_STRING  ""

#define TIO_FILL_RADIANCE  (-TIO_FILL_FLOAT)
#define TIO_FILL_IRRADIANCE  (-TIO_FILL_FLOAT)

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
   TIO_PROC_LEVEL_2,
   TIO_PROC_LEVEL_3
};

typedef struct
{
   char *name;
   size_t num_xtrack;
   size_t num_channels;
}
TIO_Scan_Group_Type;

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

/** Create a template Level 1 radiance data file
 * @param  ncid          Index returned by nc_create
 * @param  num_steps     Number of mirror steps
 * @param  num_sgrps     Number of scan groups
 * @param  sgrps         Array of TIO_Scan_Group_Type structs
 * @return 0 on success, -1 on error
 */
extern int TIO_l1_radiance_template (int ncid, size_t num_steps, int num_sgrps,
                                     TIO_Scan_Group_Type *sgrps);

/** Create a template Level 1 irradiance data file
 * @param  ncid          Index returned by nc_create
 * @param  num_steps     Number of mirror steps
 * @param  num_sgrps     Number of scan groups
 * @param  sgrps         Array of TIO_Scan_Group_Type structs
 * @return 0 on success, -1 on error
 */
extern int TIO_l1_irradiance_template (int ncid, size_t num_steps, int num_sgrps,
                                       TIO_Scan_Group_Type *sgrps);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
