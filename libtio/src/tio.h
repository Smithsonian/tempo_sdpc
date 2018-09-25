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

/** Convert TAI sec since TEMPO epoch to UTC sec since Unix epoch
 * @param[in]  tempo_time    TAI seconds since the TEMPO epoch
 * @param[out] utc_time      UTC seconds since the Unix epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, while the UTC
 * time scale includes leap second corrections at irregular
 * intervals.
 */
extern int tio_time_tempo_to_utc (double tempo_time, double *utc_time);

/** Convert UTC sec since the Unix epoch to TAI sec since the TEMPO epoch
 * @param[in]  utc_time      UTC seconds since the Unix epoch
 * @param[out] tempo_time    TAI seconds since the TEMPO epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, while the UTC
 * time scale includes leap second corrections at irregular
 * intervals.
 */
extern int tio_time_utc_to_tempo (double utc_time, double *tempo_time);

/** Convert TAI sec since the TEMPO epoch to TAI sec since the Unix epoch
 * @param[in]  tempo_time    TAI seconds since the TEMPO epoch
 * @param[out] tai_time      TAI seconds since the Unix epoch
 * @return 0 on success, -1 on error
 *
 * The TAI time scale advances monotonically, with no leap second
 * corrections.
 */
extern int tio_time_tempo_to_tai (double tempo_time, double *tai_time);

/** Unix time_t value for the TEMPO epoch
 * @return Unix time_t value for the TEMPO epoch
 *
 * The time_t value may be used for direct comparison with the Unix
 * NTP-synchronized clock value, which is a count of seconds
 * elapsed since the Unix epoch, January 1, 1970, 00:00:00 (UTC).
 * The NTP-synchronized Unix clock includes leap second corrections
 * and is therefore non-monotonic.
 */
extern double tio_tempo_epoch_timet (void);

/** Convert TAI sec since TEMPO epoch to calendar date and time, UTC
 * @param[in] tempo_time   Elapsed TAI seconds since the TEMPO epoch
 * @param[out]  year       4-digit year
 * @param[out]  month      month, 1-12
 * @param[out]  day        day of month, 1-31
 * @param[out]  hour       UTC hour of the day, 0-23.999...
 * @return 0 on success, -1 on error
 */
extern int tio_time_tempo_to_utc_caldate
(double tempo_time, int *year, int *month, int *day, double *hour);

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

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
