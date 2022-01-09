/*! @file
 *  @brief C public interface
 */
#ifndef __TIO_INCLUDE_H__
#define __TIO_INCLUDE_H__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <stddef.h>
#include <netcdf.h>

/* Number of dimensions in the nominal_wavelength array.
 * In SDPCv2, num_dims=1 -> nominal_wavelength(spectral_channel)
 * In SDPCv3, num_dims=2 -> nominal_wavelength(xtrack,spectral_channel)
 */
#ifndef TIO_NOMINAL_WAVELEN_NUM_DIMS
# define TIO_NOMINAL_WAVELEN_NUM_DIMS 2
#endif

/* Maximum number of array dimensions */
#define TIO_MAX_VAR_DIMS   7

/* Maximum name length */
#define TIO_MAX_NAME_LEN  128

/* Maximum product short name length */
#define TIO_MAX_SHORT_NAME_LEN  16

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
#define TIO_STRING NC_STRING

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

typedef struct
{
   int varid;   /**< variable ids assigned by nc_def_var */
   int type;    /**< storage data type */
   int ndims;   /**< number of dimensions */
   int dimids[TIO_MAX_VAR_DIMS];      /**< dimension ids assigned by nc_def_dim */
   size_t dimlens[TIO_MAX_VAR_DIMS];  /**< dimension sizes */
}
TIO_Var_Info_Type;

typedef struct
{
   char *name;
   char *text;
}
TIO_Attr_Text_Type;

typedef struct TIO_Scan_Ident_Type TIO_Scan_Ident_Type;

/** Create a file
 * @param[in]  path   Path to file
 * @param[in]  cmode  Create mode
 * @param[out] ncid   Pointer to id number
 * @return 0 on success, -1 on error
 */
extern int TIO_create (const char *path, int cmode, int *ncid);

/** Open an existing file
 * @param[in]  path   Path to existing file
 * @param[in]  omode  Open mode
 * @param[out] ncid   Pointer to id number of opened file
 * @return 0 on success, -1 on error
 */
extern int TIO_open (const char *path, int omode, int *ncid);

/** Close an open file
 * @param[in] ncid   Id number of currently open file
 * @return 0 on success, -1 on error
 */
extern int TIO_close (int ncid);

/** Sync a file to disk
 * @param[in] ncid   Id number of currently open file
 * @return 0 on success, -1 on error
 */
extern int tio_sync (int ncid);

/** Get id of a group given the full path to it
 * @param[in]  parent_ncid   File index
 * @param[in]  path          Full name path to the group
 * @param[out] grp           group index
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_grp (int parent_ncid, const char *path, int *grp);

/** Get the variable index for a given variable name.
 * @param  grp    Index of group containing the variable
 * @param  name   Variable name string
 * @param  varid  Pointer to the integer that will receive the variable index
 * @return 0 on success, -1 on error
 */
extern int tio_inq_varid (int grp, const char *name, int *varid);

/** Get information about a variable using its name.
 * @param  grp    Index of group containing the variable
 * @param  name   Variable name string
 * @param  info   Pointer to the structure that will receive the information.
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_var (int grp, const char *name, TIO_Var_Info_Type *info);

/*! Write a block of values to an N-dimensional array variable.
 * @param[in] grp  Group identifier
 * @param[in] name  Variable name.
 * @param[in] start  Starting offset in each dimension of the block to be
 *                    written to, with dimensions specified in C-order, so that
 *                    \a start[0] varies slowest and \a start[num_dims-1]
 *                    varies fastest.
 * @param[in] count Count of values in each dimension of the block to be written to.
 * @param[in] type   Internal data type in which the variable values are stored.
 *                   If this differs from the external type stored in the
 *                   file, type conversion will be attempted.
 * @param[in] data  Address where block of variable values will be stored.
 * @return 0 indicates success, -1 indicates failure.
 */
extern int TIO_put_var_section (int grp, const char *name,
                                const int *start, const int *count, int type,
                                const void *data);

/*! Read a block of values from an N-dimensional array variable.
 * @param[in] grp  Group identifier
 * @param[in] name  Variable name.
 * @param[in] start  Starting offset in each dimension of the block to be read,
 *                    with dimensions specified in C-order, so that
 *                    \a start[0] varies slowest and \a start[N-1]
 *                    varies fastest.
 * @param[in] count Count of values in each dimension of the block to be read.
 * @param[in] type   Internal data type in which to store variable values.
 *                   If this differs from the external type stored in the
 *                   file, type conversion will be attempted.
 * @param[out] data  Address where block of variable values will be stored.
 * @return 0 indicates success, -1 indicates failure.
 */
extern int TIO_get_var_section (int grp, const char *name,
                                const int *start, const int *count, int type,
                                void *data);

/** Enable or disable variable-specific I/O methods
 * @param  name        Name of the variable
 * @param  get_enable  Control the variable-specific input method:
 *                     non-zero to enable, zero to disable.
 * @param  put_enable  Control the variable-specific output method:
 *                     non-zero to enable, zero to disable.
 * @return 0 indicates success, -1 indicates the variable name is unrecognized
 */
extern int _TIO_set_io_method_enable (const char *name,
                                      int get_enable, int put_enable);

/** Query the type and size of an attribute
 * @param  grp      Index of group containing the attribute
 * @param  varid    Index of the corresponding variable.
 * @param  attname  Name of the attribute
 * @param  xtype    Pointer to attribute type
 * @param  len      Pointer to attribute length
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_att (int grp, int varid, const char *attname,
                        int *xtype, int *len);

/** Write an attribute value
 * @param  grp      Index of group containing the attribute
 * @param  varid    Index of the corresponding variable.
 * @param  attname  Name of the attribute
 * @param  xtype    Attribute type
 * @param  len      Attribute length
 * @param  att      Pointer to attribute value
 * @return 0 on success, -1 on error
 */
extern int TIO_put_att (int grp, int varid, const char *attname,
                        int xtype, int len, const void *att);

/** Write multiple text attributes
 * @param[in]  grp      Index of group containing the attribute
 * @param[in]  varid    Index of the corresponding variable.
 * @param[in]  attrs    Pointer to array of TIO_Attr_Text_Type
 *                      (last element has attr->name=NULL)
 * @return 0 on success, -1 on error
 */
extern int TIO_put_text_attrs (int grp, int varid,
                               const TIO_Attr_Text_Type *attrs);

/** Read an attribute value
 * @param  grp      Index of group containing the attribute
 * @param  varid    Index of the corresponding variable.
 * @param  attname  Name of the attribute
 * @param  xtype    Attribute type
 * @param  att      Pointer to attribute value
 * @return 0 on success, -1 on error
 */
extern int TIO_get_att (int grp, int varid, const char *attname,
                        int xtype, void *att);

/** Free string array returned by TIO_get_att
 * @param len    The number of character arrays
 * @param data   Pointer to the array of character arrays
 * @return 0 on success, -1 on error
 */
extern int TIO_free_string (size_t len, char **data);

/** Copy variable attributes between ncids
 * @param  ncid_infile     Input ncid
 * @param  id_var_infile   Input variable id (in input ncid)
 * @param  dontcopy_attr   (Optional) Pointer to function that
 *                         indicates whether or not an attribute
 *                         should be copied.
 * @param  ncid            Output ncid.
 * @param  id_var          Output variable id (in output ncid).
 * @return 0 on success, -1 on error
 */
extern int TIO_copy_attrs (int ncid_infile, int id_var_infile,
                           int (*dontcopy_attr)(const char *),
                           int ncid, int id_var);

/** Read the name of the dimension with a given dimid.
 * @param[in]   ncid     File ncid
 * @param[in]   dimid    Id number of the dimension
 * @param[out]  dimname  dimension name
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_dimname (int grp, int dimid, char *dimname);

/** Read the id number of a named dimension.
 * @param[in]   ncid     File ncid
 * @param[in]   dimname  dimension name
 * @param[out]  dimid    Id number of the dimension
 * @return same as nc_inq_dimid
 */
extern int TIO_inq_dimid (int grp, const char *dimname, int *dimid);

/** Read the id number and size of a named dimension.
 * @param[in]   ncid     File ncid
 * @param[in]   dimname  dimension name
 * @param[out]  dimid    Id number of the dimension
 * @param[out]  dimlen    size of the dimension
 * @return same as nc_inq_dimid and nc_inq_dimlen
 */
extern int TIO_inq_dim (int grp, const char *dimname,
                        int *dimid, size_t *dimlen);

/** Write a global attribute containing the git commit hash
 * @param  ncid     File ncid
 * @param  attname  Name of the attribute (uses default if NULL)
 * @return 0 on success, -1 on error
 */
extern int TIO_put_git_commit_hash (int ncid, const char *attname);

/** Define a group
 * @param  parent_ncid  The ncid that will contain the new group
 * @param  path         The path to the new group
 * @param  new_ncid     Pointer to ncid of the new group
 * @return 0 on success, -1 on error
 */
extern int TIO_def_grp (int parent_ncid, const char *path, int *new_ncid);

/** Define a dimension
 * @param[in]   ncid         The ncid that will contain the dimension
 * @param[in]   dimname      The name of the dimension
 * @param[in]   dimlen       The size_t size of the dimension, or NC_UNLIMITED
 * @param[out]  dimid        The id number of the new dimension
 * @return 0 on success, -1 on error
 */
extern int TIO_def_dim (int ncid, const char *dimname, size_t dimlen, int *dimid);

/** Define a variable
 * @param[in]   ncid      The ncid that will contain the variable
 * @param[in]   name      The name of the variable
 * @param[in]   type      The type of the variable
 * @param[in]   num_dims  The number of dimensions
 * @param[in]   dimids    The id number of each of the dimensions
 * @param[out]  varid     The id number of the new variable
 * @return 0 on success, -1 on error
 */
extern int TIO_def_var (int ncid, const char *name, int type,
                        int num_dims, const int *dimids, int *varid);

/** Define a variable's fill value
 * @param  grp         The group containing the variable
 * @param  varid       Variable id number
 * @param  no_fill     When no_fill mode is non-zero, fill values will not be written
 * @param  fill_value  Pointer to the fill value
 * @return 0 on success, -1 on error
 */
extern int TIO_def_var_fill (int grp, int varid, int no_fill, const void *fill_value);

/** Query a variable's fill value
 * @param  grp         The group containing the variable
 * @param  varid       Variable id number
 * @param  no_fill     When no_fill mode is non-zero, fill values will not be written
 * @param  fill_value  Pointer to the fill value
 * @return 0 on success, -1 on error
 */
extern int TIO_inq_var_fill (int grp, int varid, int *no_fill, void *fill_value);

/** Define a variable's deflate values
 * @param  grp            The group containing the variable
 * @param  varid          Variable id number
 * @param  shuffle        Non-zero to use shuffle
 * @param  deflate        Non-zero to use compresson
 * @param  deflate_level  Compression level to use
 * @return 0 on success, -1 on error
 */
extern int TIO_def_var_deflate (int grp, int varid,
                                int shuffle, int deflate, int deflate_level);

/** Define a variable's chunking configuration
 * @param grp          The group containing the variable
 * @param varid        Variable id number
 * @param storage      Storage class specifier NC_CHUNKED | NC_CONTIGUOUS
 * @param chunksizep   Array of chunk sizes, one per dimension
 * @return 0 on success, -1 on error
 */
extern int TIO_def_var_chunking (int grp, int varid,
                                 int storage, size_t *chunksizep);

/** Read the fill value of a named variable and convert it to a specified type
 * @param  grp    The group containing the variable
 * @param  name   The variable name
 * @param  type   The desired netcdf variable type, e.g. NC_INT, NC_FLOAT
 * @param  value  Pointer to the returned fill value
 * @return 0 on success, -1 on error
 */
extern int TIO_get_fill_value (int grp, const char *name, int type, void *value);

typedef struct
{
   char *name;
   int value;
}
TIO_Enum_Type;
#define TIO_ENUM_TABLE_END {NULL,0}

/** Define an enum data type using a provided array of name-integer pairs.
 * @param  grp    The group that is to contain the enum data type
 * @param  name   The name of the enum data type
 * @param  base_type  The desired netcdf integer base-type, e.g. NC_INT
 * @param  enum_table  Pointer to a NULL-terminated array of TIO_Enum_Type structs
 * @param[out] enum_typeid   The data type index of the newly defined enum data type
 * @return 0 on success, -1 on error
 */
extern int TIO_define_enum_table (int grp, const char *name, int base_type,
                                  const TIO_Enum_Type *enum_table, int *enum_typeid);

/** Convert TAI sec since X epoch to UTC sec since Unix epoch
 * @param[in]  taix_time     TAI seconds since the X epoch
 * @param[out] utc_time      UTC seconds since the Unix epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, while the UTC
 * time scale includes leap second corrections at irregular
 * intervals.
 */
extern int tio_time_taix_to_utc (double taix_time, double *utc_time);

/** Convert UTC sec since the Unix epoch to TAI sec since the X epoch
 * @param[in]  utc_time     UTC seconds since the Unix epoch
 * @param[out] taix_time    TAI seconds since the X epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, while the UTC
 * time scale includes leap second corrections at irregular
 * intervals.
 */
extern int tio_time_utc_to_taix (double utc_time, double *taix_time);

/** Convert UTC time string to TAI sec since the X epoch
 * @param[in]  utc_str    UTC time stamp of the form YYYY-MM-DDThh:mm:ssZ
 * @param[out] taix_time  TAI seconds since the X epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, while the UTC
 * time scale includes leap second corrections at irregular
 * intervals.
 */
extern int tio_time_utcstr_to_taix (const char *utc_str, double *taix_sec);

/** Convert UTC sec since the Unix epoch to TAI sec since the Unix epoch
 * @param[in]  utc_time    UTC seconds since the Unix epoch
 * @param[out] tai_time    TAI seconds since the Unix epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, while the UTC
 * time scale includes leap second corrections at irregular
 * intervals.
 */
extern int tio_time_utc_to_tai (double utc_time, double *tai_time);

/** Convert TAI sec since the X epoch to TAI sec since the Unix epoch
 * @param[in]  taix_time     TAI seconds since the X epoch
 * @param[out] tai_time      TAI seconds since the Unix epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, with no leap second
 * corrections.
 */
extern int tio_time_taix_to_tai (double taix_time, double *tai_time);

/** Unix time_t value for the X epoch
 * @return Unix time_t value for the X epoch
 *
 * The time_t value may be used for direct comparison with the Unix
 * NTP-synchronized clock value, which is a count of seconds
 * elapsed since the Unix epoch, January 1, 1970, 00:00:00 (UTC).
 * The NTP-synchronized Unix clock includes leap second corrections
 * and is therefore non-monotonic.
 */
extern double tio_time_taix_epoch_timet (void);

/** Set the X epoch using a UTC time string
 * @param[in]  utc   Epoch specification of the form YYYY-MM-DDThh:mm:ssZ
 * @return 0 on success, -1 on error
 */
extern int tio_time_set_taix_epoch (const char *utc_string);

/** Set the X epoch defined as a count of seconds (UTC) since the Unix epoch
 * @param[in]  tt    Epoch specification expressed as a count of seconds since
 *                   the Unix epoch (UTC)
 * @return 0 on success, -1 on error
 */
extern int tio_time_set_taix_epoch_timet (time_t tt);

/** Set the X epoch using the time_reference attribute stored in a file
 * @param[in]  ncid   Device id of a netCDF file open for reading
 * @return 0 on success, -1 on error
 */
extern int tio_use_file_epoch (int ncid);

/** Convert TAI sec since X epoch to calendar date and time, UTC
 * @param[in] taix_time   Elapsed TAI seconds since the X epoch
 * @param[out]  year       4-digit year
 * @param[out]  month      month, 1-12
 * @param[out]  day        day of month, 1-31
 * @param[out]  hour       hour of the day, 0-23.999...
 * @return 0 on success, -1 on error
 */
extern int tio_time_taix_to_utc_caldate
(double taix_time, int *year, int *month, int *day, double *hour);

/** Convert TAI sec since X epoch to calendar date and day of year, UTC
 * @param[in] taix_time   Elapsed TAI seconds since the X epoch
 * @param[out]  year      4-digit year
 * @param[out]  yday      number of days since 1 Jan, [0..365]
 * @return 0 on success, -1 on error
 */
extern int tio_time_taix_to_yearday (double taix_time, int *year, int *yday);

/** Compute satellite local-time day number given a count of TAI sec since X epoch
 * @param[in]   taix_time   Elapsed TAI seconds since the X epoch
 * @param[out]  sat_day     The day number in the satellite-local time zone.
 * @return 0 on success, -1 on error
 */
extern int tio_time_sat_local_day_number (double taix, double *sat_day);

/** Read one line from an ASCII text file
 * @param[out] linep   Pointer to an allocated string
 * @param[out] lenp    The length of the allocated string (NULL is ok)
 * @param[in]  fp      An open file pointer.
 * @return 0 on success, -1 on error
 *
 * The string returned will include the terminating newline if
 * one is present in the file.
 */
extern int tio_fgets (char **linep, size_t *lenp, FILE *fp);

extern int tio_def_var_radiance_status_flag (int grp);

/** Append a string to the global history attribute
 * @param[in]  ncid   netcdf file id
 * @param[in]  str    pointer to the history string to be appended
 * @return 0 on success, -1 on error
 */
extern int tio_append_history (int ncid, const char *str);
extern int tio_set_cmdline (int argc, char **argv);
extern int tio_history_append_cmdline (int ncid);

/** Concatenate string tokens
 * @param[in]  argc   number of string tokens
 * @param[in]  argv   array of pointers to string tokens
 * @param[in]  pstr   optional pointer to space for result string (NULL is ok)
 * @param[in]  len_pstr  length of space pointed to by pstr (ignored when pstr=NULL)
 * @return NULL on error
 *         When pstr=NULL, the return value is an allocated string.
 *         When pstrt!=NULL, the return value is pstr.
 */
extern char *tio_concat_argv (int argc, char **argv, char *pstr, size_t len_pstr);

/* Metadata interface */

enum
{
   TIO_META_TYPE_UNDEFINED = 0,
   TIO_META_TYPE_DOUBLE,
   TIO_META_TYPE_FLOAT,
   TIO_META_TYPE_INT,
   TIO_META_TYPE_UINT,
   TIO_META_TYPE_CHAR,
   TIO_META_TYPE_STRING
};

typedef struct TIO_Meta_Type TIO_Meta_Type;

/** Initialize a data structure for storing metadata keyword values
 * @return NULL on error, otherwise the return value is an opaque
 *         pointer of type \a TIO_Meta_Type
 */
extern TIO_Meta_Type *tio_meta_open (void);

/** Free memory associated with an instance of \a TIO_Meta_Type
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 */
extern void tio_meta_close (TIO_Meta_Type *meta);

/** Assign a scalar or array value to a metadata keyword
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  name   Metadata keyword name
 * @param[in]  value_type  Data type of the metadata keyword value.  Must be one of
 *                         \a TIO_META_TYPE_${t} where \a ${t} is \a DOUBLE|FLOAT|INT|UINT|STRING
 * @param[in]  num_values  Number of values to be assigned
 * @param[in]  values      \a{void *} pointer to the keyword values
 * @return 0 on success, -1 on error
 *
 * If a keyword with the same name has already been defined, the previous keyword
 * will be silently replaced.
 *
 * @code
 *  status = tio_meta_set (meta, "FOO", TIO_META_TYPE_STRING, 1, "Example scalar");
 *  status = tio_meta_set (meta, "BAR", TIO_META_TYPE_STRING, 5, array_of_strings);
 *  status = tio_meta_set (meta, "BAZ", TIO_META_TYPE_DOUBLE, 1, &a_double);
 * @endcode
 *
 * @see tio_meta_append_string
 */
extern int tio_meta_set (TIO_Meta_Type *meta, const char *name,
                         int value_type, int num_values, const void *values);

extern int tio_meta_set_acdd_geospatial_bounds (TIO_Meta_Type *meta,
                                                const float *lon, const float *lat, int num);
extern int tio_meta_set_odl_bounding_polygon (TIO_Meta_Type *meta,
                                              const float *lon, const float *lat, int num);

/** Read a keyword value from an attribute in the \a metadata group.
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  grp    netCDF group index, open for reading
 * @param[in]  name   Metadata keyword name
 * @param[in]  value_type  Data type of the metadata keyword value.  Must be one of
 *                         \a TIO_META_TYPE_${t} where \a ${t} is \a DOUBLE|FLOAT|INT|UINT|STRING
 *
 * @warning The current implementation requires @code value_type = TIO_META_TYPE_STRING @endcode
 *
 * @see tio_meta_expand_file
 * @see tio_meta_write_ncattr
 */
extern int tio_meta_ncinit (TIO_Meta_Type *meta, int grp, const char *name,
                            int value_type);

/** Append a value to metadata string keyword
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  name   Metadata keyword name
 * @param[in]  str    Pointer to a string
 * @return 0 on success, -1 on error
 *
 * If the named keyword does not exist, it will be initialized using the value
 * provided. String values are appended as new items in an array of strings
 * (strings are not concatenated).
 *
 * @see tio_meta_set
 */
extern int tio_meta_append_string (TIO_Meta_Type *lst, const char *name, const char *str);

/** Enable/disable keyword-value substitution for a named keyword
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  name   Metadata keyword name
 * @param[in]  noexpand  If non-zero, disable keyword expansion for the named keyword.
 *                       Otherwise, enable keyword expansion for the named keyword.
 *
 * This has no effect on how metadata keywords are written to netCDF files.
 *
 * @see tio_meta_expand_file
 * @see tio_meta_write_ncattr
 */
extern int tio_meta_set_noexpand (TIO_Meta_Type *meta, const char *name, int noexpand);

/** Expand keywords appearing in an ASCII input stream
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  fp_template  \a FILE pointer to an input stream containing keyword variables
 * @param[in]  fp_outfile  \a FILE pointer to an output stream with keyword variables
 *                         replaced with assigned keyword values.
 * @return 0 on success, -1 on error
 *
 * Keywords in the input stream must have one of the following forms:
 * @verbatim
 *       ${KEYWORDNAME}      expands to the keyword value
 *       ${#KEYWORDNAME}     expands to the number of keyword values
 * @endverbatim
 *
 * @see tio_meta_expand_file
 * @see tio_meta_write_ncattr
 */
extern int tio_meta_expand_stream (const TIO_Meta_Type *meta, FILE *fp_template,
                                   FILE *fp_outfile);

/** Expand keywords appearing in an ASCII template file
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  infile  If non-NULL, this provides the path to a template file containing
 *                     keyword variables.
 * @param[in]  outfile_root  The output filename will be constructed by
 *                           appending a @<".met"@> extension to this string.
 * @return 0 on success, -1 on error
 *
 * When @code infile==NULL @endcode, it is assumed that we are expanding a template
 * named @code ${outfile_root}.met @endcode.  In this case, a temporary file is
 * created to hold the expanded output, and then that temporary file is renamed
 * to the original template name.  Effectively, the original file is over-written
 * by a keyword-expanded version.
 *
 * Keywords in the input stream must have one of the following forms:
 * @verbatim
 *       ${KEYWORDNAME}      expands to the keyword value
 *       ${#KEYWORDNAME}     expands to the number of keyword values
 * @endverbatim
 *
 * @see tio_meta_expand_stream
 * @see tio_meta_write_ncattr
 */
extern int tio_meta_expand_file (const TIO_Meta_Type *meta, const char *infile,
                                 const char *outfile_sans_extname);

/** Set values for several standard metadata keywords
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  product_file_name   Basename of the data product file
 * @param[in]  product_short_name  Data product shortname (NULL is ok)
 * @param[in]  product_versionid   Data product versionid
 * @param[in]  pge_version_string  Version number string of the program used to generate
 *                                 the data product
 * @return 0 on success, -1 on error
 *
 * The provided keyword values are substituted for the following keyword names:
 * @verbatim
 * local_granule_id
 * shortname
 * version_id
 * pge_version
 * @endverbatim
 */
extern int tio_meta_set_standard (TIO_Meta_Type *meta,
                                  const char *product_file_name,
                                  const char *product_short_name,
                                  int product_versionid,
                                  const char *pge_version_string);

/** Set values for standard observation time interval keywords
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  ncid   Integer device code for a netCDF4/HDF5 data product file open for reading
 * @return 0 on success, -1 on error
 *
 * The following timestamps are parsed:
 * @verbatim
 *    time_coverage_start
 *    time_coverage_end
 * @endverbatim
 * to assign values to the following standard keywords
 * @verbatim
 *    RANGEBEGINNINGDATE
 *    RANGEBEGINNINGTIME
 *    RANGEENDINGDATE
 *    RANGEENDINGTIME
 * @endverbatim
 */
extern int tio_meta_set_datetime_range (TIO_Meta_Type *meta, int ncid);

/** Set values for standard observation time interval keywords for a full scan
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  lst    Pointer of type \a TIO_Scan_Ident_Type describing a full scan
 * @return 0 on success, -1 on error
 *
 * The following timestamps are parsed:
 * @verbatim
 *    time_coverage_start
 *    time_coverage_end
 * @endverbatim
 * to assign values to the following standard keywords
 * @verbatim
 *    RANGEBEGINNINGDATE
 *    RANGEBEGINNINGTIME
 *    RANGEENDINGDATE
 *    RANGEENDINGTIME
 * @endverbatim
 * where the start time is the beginning of the scan and the end time is the
 * end of the scan.
 */
extern int tio_meta_set_datetime_range_scan (TIO_Meta_Type *meta, const TIO_Scan_Ident_Type *lst);

/** Set value for standard production date keywords
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @return 0 on success, -1 on error
 *
 * The current time (UTC) is used to assign a value to the standard keyword
 * @c production_date_time
 */
extern int tio_meta_set_datetime_production (TIO_Meta_Type *meta);

/** Set value for geospatial bounding polygon keywords
 * @param[in]  meta  Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  grp   Integer device code for reading from the group in a netCDF4/HDF5
 *                   data product file containing \a longitude and \a latitude
 *                   pixel coordinates
 * @return 0 on success, -1 on error
 *
 * The following keyword values are set:
 * @verbatim
 *   polygon_longitudes      boundary polygon longitudes
 *   polygon_latitudes       boundary polygon latitudes
 *   polygon_sequence        integer indices giving the sequence in which the
 *                           (lon,lat) points trace the boundary in CCW order
 * @endverbatim
 *
 * @see __tio_make_lev1_bounding_polygon
 */
extern int tio_meta_set_lev1_bounding_polygon (TIO_Meta_Type *meta, int grp);

/** Construct the geospatial bounding polygon
 * @param[in]  grp   Integer device code for reading from the group in a netCDF4/HDF5
 *                   data product file containing \a longitude and \a latitude
 *                   pixel coordinates
 * @param[out] num   Number of polygon vertices
 * @param[out] plon  Pointer to an allocated array of longitude coordinates, in CCW
 *                   order around the polygon boundary.
 * @param[out] plat  Pointer to an allocated array of latitude coordinates, in CCW
 *                   order around the polygon boundary.
 * @return 0 on success, -1 on error
 *
 * The \a longitude and \a latitude coordinates in the file are processed to derive
 * the coordinates of a polygon that bounds the region with valid coordinates.
 *
 * @see tio_meta_set_lev1_bounding_polygon
 */
extern int __tio_make_lev1_bounding_polygon (int grp, int *num, float **plon, float **plat);

/** Write metadata keywords as attributes in a specified netCDF4/HDF5 group
 * @param[in]  meta   Pointer of type \a TIO_Meta_Type allocated by \a tio_meta_open
 * @param[in]  grp   Integer device code for writing to a group in a netCDF4/HDF5
 *                   data product file.
 * @return 0 on success, -1 on error
 *
 * @see tio_meta_expand_stream
 */
extern int tio_meta_write_ncattr (const TIO_Meta_Type *meta, int grp);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
