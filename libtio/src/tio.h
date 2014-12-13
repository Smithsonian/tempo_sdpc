#ifndef __TIO_INCLUDE_H__
#define __TIO_INCLUDE_H__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <stddef.h>

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

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
