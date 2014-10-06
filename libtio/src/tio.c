#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "cfortran.h"

#include "netcdf.h"
#include "tio.h"
#include "_tio.h"

#define EMPTY()

int _pTIO_define_dims_using_offsets (int grp, const _pDim_Offsets_Type *offsets,
                                     _pDim_Table_Type *dim_table)
{
   const _pDim_Offsets_Type *o;
   char *p = (char *)dim_table;
   int status;

   for (o = offsets; o->name != NULL; o++)
     {
        status = nc_def_dim (grp, o->name, *(size_t *)(p + o->len_offset),
                             (int *)(p + o->id_offset));
        if (NC_NOERR != status)
          {
             _pTIO_err_verror_nc (status, "%s: defining dimension %s",
                                  __func__, o->name);
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_enum (int grp, const char *name,
                       const _pEnum_Type *enum_table, int *enum_typeid)
{
   const _pEnum_Type *e;
   int status;

   if (NC_NOERR != (status = nc_def_enum (grp, NC_INT, name, enum_typeid)))
     {
        _pTIO_err_verror_nc (status, "%s: defining enum %s", __func__, name);
        return -1;
     }

   for (e = enum_table; e->name != NULL; e++)
     {
        status = nc_insert_enum (grp, *enum_typeid, e->name, &e->value);
        if (status != NC_NOERR)
          {
             _pTIO_err_verror_nc (status, "%s: inserting value %s=%d for enum %s",
                                  __func__, e->name ? e->name : "(null)",
                                  e->value, name);
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_int_attrs (int grp, int varid, const _pInt_Attr_Type *attrs)
{
   const _pInt_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_int (grp, varid, a->name, NC_INT, 1, &a->value);
        if (NC_NOERR != status)
          {
             _pTIO_err_verror_nc (status, "%s: defining int attribute %s",
                                  __func__, a->name);
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_float_attrs (int grp, int varid, const _pFloat_Attr_Type *attrs)
{
   const _pFloat_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_float (grp, varid, a->name, NC_FLOAT, 1, &a->value);
        if (NC_NOERR != status)
          {
             _pTIO_err_verror_nc (status, "%s: defining float attribute %s",
                                  __func__, a->name);
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_text_attrs (int grp, int varid, const _pText_Attr_Type *attrs)
{
   const _pText_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        size_t len = strlen(a->text) + 1;
        status = nc_put_att_text (grp, varid, a->name, len, a->text);
        if (NC_NOERR != status)
          {
             _pTIO_err_verror_nc (status, "%s: defining text attribute %s",
                                  __func__, a->name);
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_var_with_text_attrs (int grp, const char *var_name, nc_type xtype,
                                      int num_dims, const int *dimids,
                                      const _pText_Attr_Type *text_attrs,
                                      int *pvarid)
{
   int status, varid;

   status = nc_def_var (grp, var_name, xtype, num_dims, dimids, &varid);
   if (NC_NOERR != status)
     {
        _pTIO_err_verror_nc (status, "%s: defining variable %s",
                             __func__, var_name);
        return -1;
     }

   if (text_attrs != NULL)
     {
        if (-1 == _pTIO_define_text_attrs (grp, varid, text_attrs))
          return -1;
     }

   if (pvarid != NULL)
     {
        *pvarid = varid;
     }

   return 0;
}

int _pTIO_put_fillvalue_attr (int grp, int varid, nc_type xtype)
{
   int status;
   void *pfill_value;
   char fill_char = TIO_FILL_CHAR;
   short fill_short = TIO_FILL_SHORT;
   int fill_int = TIO_FILL_INT;
   unsigned int fill_uint = (unsigned int) -1;
   float fill_float = TIO_FILL_FLOAT;
   double fill_double = TIO_FILL_DOUBLE;

   switch (xtype)
     {
      case NC_CHAR: pfill_value = &fill_char;
        break;
      case NC_SHORT: pfill_value = &fill_short;
        break;
      case NC_INT: pfill_value = &fill_int;
        break;
      case NC_UINT: pfill_value = &fill_uint;
        break;
      case NC_FLOAT: pfill_value = &fill_float;
        break;
      case NC_DOUBLE: pfill_value = &fill_double;
        break;
      default:
        _pTIO_err_verror ("%s: invalid fill value type xtype=%d", xtype);
        return -1;
     }

   status = nc_put_att (grp, varid, _FillValue, xtype, 1, pfill_value);
   if (NC_NOERR != status)
     {
        _pTIO_err_verror_nc (status, "writing fill value to grp=%d varid=%d",
                             grp, varid);
        return -1;
     }

   return 0;
}

int _pTIO_define_processing_level (int grp, int level)
{
   int status, enum_typeid;
   static _pEnum_Type enum_table[] =
     {
        {"level-0",  TIO_PROC_LEVEL_0},
        {"level-1a", TIO_PROC_LEVEL_1A},
        {"level-1b", TIO_PROC_LEVEL_1B},
        {"level-2",  TIO_PROC_LEVEL_2},
        {"level-3",  TIO_PROC_LEVEL_3},
        _pENUM_TABLE_END
     };

   if (-1 == _pTIO_define_enum (grp, "processing_level_enum", enum_table, &enum_typeid))
     return -1;
   status = nc_put_att (grp, NC_GLOBAL,
                        "processing_level", enum_typeid, 1, &level);
   if (NC_NOERR != status)
     {
        _pTIO_err_verror_nc (status, "%s: defining processing_level attribute", __func__);
        return -1;
     }

   return 0;
}

int TIO_inq_var (int grp, const char *name, TIO_Var_Info_Type *info)
{
   int status, i;

   if ((name == NULL) || (info == NULL))
     {
        _pTIO_err_verror ("%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_varid (grp, name, &info->varid)))
     {
        _pTIO_err_verror_nc (status, "%s: accessing variable %s in group %d",
                             __func__, name, grp);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_varndims (grp, info->varid, &info->ndims)))
     {
        _pTIO_err_verror_nc (status, "%s: accessing variable %s in group %d",
                             __func__, name, grp);
        return -1;
     }

   /* convenient default for the ndims=0 case */
   info->dimlens[0] = 1;

   if (info->ndims == 0)
     return 0;

   if (NC_NOERR != (status = nc_inq_vardimid (grp, info->varid, info->dimids)))
     {
        _pTIO_err_verror_nc (status,
                             "%s: accessing variable %s dimids in group %d",
                             __func__, name, grp);
        return -1;
     }
   for (i = 0; i < info->ndims; i++)
     {
        status = nc_inq_dimlen  (grp, info->dimids[i], &info->dimlens[i]);
        if (NC_NOERR != status)
          {
             _pTIO_err_verror_nc (status,
                                  "%s: accessing dimension %d length in group %d",
                                  __func__, info->dimids[i], grp);
             return -1;
          }
     }

   return 0;
}

/* #define TEST_WAVELENGTH_METHODS 1 */
#ifdef TEST_WAVELENGTH_METHODS
static int get_wavelengths (int grp, size_t start0, size_t count0, int xtype,
                            void *data)
{
   fprintf (stderr, "=====> called get_wavelengths\n");
   return 0;
}
static int put_wavelengths (int grp, size_t start0, size_t count0, int xtype,
                            const void *data)
{
   fprintf (stderr, "=====> called put_wavelengths\n");
   return 0;
}
#endif

typedef struct
{
   char *name;
   int (*get)(int, size_t, size_t, int,       void *);
   int (*put)(int, size_t, size_t, int, const void *);
}
IO_Methods_Type;
#define IO_METHODS_END {NULL,NULL,NULL}

static const
IO_Methods_Type *find_io_methods (const char *name,
                                  const IO_Methods_Type *io_methods)
{
   const IO_Methods_Type *m;

   if ((name == NULL) || (io_methods == NULL))
     return NULL;

   for (m = io_methods; m->name != NULL; m++)
     {
        if (0 == strcmp (name, m->name))
          return m;
     }

   return NULL;
}

static IO_Methods_Type IO_Methods[] =
{
#ifdef TEST_WAVELENGTH_METHODS
   {"wavelength", get_wavelengths, put_wavelengths},
#endif
   IO_METHODS_END
};

#define TIO_IO_VAR_SECTION(action,const_qual) \
int TIO_##action##_var_section (int grp, const char *name, \
                                size_t start0, size_t count0, int xtype, \
                                const_qual void *data) \
{ \
   TIO_Var_Info_Type info; \
   int status, varid, ndims; \
   size_t start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS]; \
   const IO_Methods_Type *io_method; \
 \
   io_method = find_io_methods (name, (const IO_Methods_Type *)&IO_Methods); \
   if (NULL != io_method) \
     { \
        return io_method->action (grp, start0, count0, xtype, data); \
     } \
 \
   start[0] = start0; \
   count[0] = count0; \
 \
   if (-1 == TIO_inq_var (grp, name, &info)) \
     return -1; \
 \
   ndims = info.ndims; \
   varid = info.varid; \
 \
   if (ndims <= 0) \
     { \
        _pTIO_err_verror ("%s: variable %s has ndims=%d", \
                          __func__, name, ndims); \
        return -1; \
     } \
   else if (ndims > 1) \
     { \
        int i; \
        for (i = 1; i < ndims; i++) \
          { \
             start[i] = 0; \
             count[i] = info.dimlens[i]; \
          } \
     } \
 \
   switch (xtype) \
     { \
      case NC_BYTE: \
        /* drop */ \
      case NC_CHAR: \
        status = nc_##action##_vara_text (grp, varid, start, count, (const_qual char *)data); \
        break; \
      case NC_UBYTE: \
        status = nc_##action##_vara_uchar (grp, varid, start, count, (const_qual unsigned char *)data); \
        break; \
      case NC_SHORT: \
        status = nc_##action##_vara_short (grp, varid, start, count, (const_qual short *)data); \
        break; \
      case NC_USHORT: \
        status = nc_##action##_vara_ushort (grp, varid, start, count, (const_qual unsigned short *)data); \
        break; \
      case NC_INT: \
        status = nc_##action##_vara_int (grp, varid, start, count, (const_qual int *)data); \
        break; \
      case NC_UINT: \
        status = nc_##action##_vara_uint (grp, varid, start, count, (const_qual unsigned int *)data); \
        break; \
      case NC_INT64: \
        status = nc_##action##_vara_longlong (grp, varid, start, count, (const_qual long long *)data); \
        break; \
      case NC_UINT64: \
        status = nc_##action##_vara_ulonglong (grp, varid, start, count, (const_qual unsigned long long *)data); \
        break; \
      case NC_FLOAT: \
        status = nc_##action##_vara_float (grp, varid, start, count, (const_qual float *)data); \
        break; \
      case NC_DOUBLE: \
        status = nc_##action##_vara_double (grp, varid, start, count, (const_qual double *)data); \
        break; \
      default: \
        _pTIO_err_verror ("%s: accessing variable %s using invalid type (xtype=%d)", \
                          __func__, name, xtype); \
        return -1; \
     } \
 \
   if (status != NC_NOERR) \
     { \
        _pTIO_err_verror_nc (status, "%s: accessing variable %s in group %d", \
                             __func__, name, grp); \
        return -1; \
     } \
 \
   return 0; \
}

TIO_IO_VAR_SECTION(get,EMPTY())
TIO_IO_VAR_SECTION(put,const)
#if 0
}
#endif

int TIO_inq_att (int grp, const char *varname, const char *attname,
                 int *xtype, size_t *len)
{
   int status, varid;

   if (NULL == attname)
     {
        _pTIO_err_verror ("%s: got a NULL pointer", __func__);
        return -1;
     }

   if (varname == NULL)
     varid = NC_GLOBAL;
   else if (NC_NOERR != (status = nc_inq_varid (grp, varname, &varid)))
     {
        _pTIO_err_verror_nc (status, "%s: accessing variable %s",
                             __func__, varname);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_att (grp, varid, attname, xtype, len)))
     {
        _pTIO_err_verror_nc (status, "%s: accessing attribute %s (varid=%d)",
                             __func__, attname, varid);
        return -1;
     }

   return 0;
}

int TIO_put_att (int grp, const char *varname, const char *attname,
                 int xtype, size_t len, const void *att)
{
   int status, varid;
   int file_atttype;

   if ((NULL == attname) || (att == NULL))
     {
        _pTIO_err_verror ("%s: got a NULL pointer", __func__);
        return -1;
     }

   if (varname == NULL)
     varid = NC_GLOBAL;
   else if (NC_NOERR != (status = nc_inq_varid (grp, varname, &varid)))
     {
        _pTIO_err_verror_nc (status, "%s: accessing variable %s",
                             __func__, varname);
        return -1;
     }

   /* Enforce attribute type-checking.
    * By default, netCDF will happily overwrite an existing
    * attribute with a different type value, but I don't think
    * we want that behavior.
    */
   status = nc_inq_atttype (grp, varid, attname, &file_atttype);
   if (NC_NOERR == status)
     {
        /* if the attribute exists, make sure we're writing
         * a value of the same type.
         */
        if (xtype != file_atttype)
          {
             _pTIO_err_verror ("%s: writing attribute %s: value type=%d, file value type=%d",
                               __func__, attname, xtype, file_atttype);
             return -1;
          }
     }
   else if (NC_ENOTATT != status)
     {
        _pTIO_err_verror_nc (status, "%s: attribute %s (varid=%d) query failed",
                             __func__, attname, varid);
        return -1;
     }

   status = nc_put_att (grp, varid, attname, xtype, len, att);
   if (NC_NOERR != status)
     {
        _pTIO_err_verror_nc (status, "%s: writing attribute %s (varid=%d)",
                             __func__, attname, varid);
        return -1;
     }

   return 0;
}

int TIO_get_att (int grp, const char *varname, const char *attname,
                 int xtype, void *att)
{
   int status, varid;

   if ((NULL == attname) || (att == NULL))
     {
        _pTIO_err_verror ("%s: got a NULL pointer", __func__);
        return -1;
     }

   if (varname == NULL)
     varid = NC_GLOBAL;
   else if (NC_NOERR != (status = nc_inq_varid (grp, varname, &varid)))
     {
        _pTIO_err_verror_nc (status, "%s: accessing variable %s",
                             __func__, varname);
        return -1;
     }

   /* Use the netCDF type-safe interface so that netCDF
    * performs any needed type conversion for built-in types.
    */
   switch (xtype)
     {
      case NC_BYTE:
        /* drop */
      case NC_CHAR:
        status = nc_get_att_text (grp, varid, attname, (char *)att);
        break;
      case NC_UBYTE:
        status = nc_get_att_uchar (grp, varid, attname, (unsigned char *)att);
        break;
      case NC_SHORT:
        status = nc_get_att_short (grp, varid, attname, (short *)att);
        break;
      case NC_USHORT:
        status = nc_get_att_ushort (grp, varid, attname, (unsigned short *)att);
        break;
      case NC_INT:
        status = nc_get_att_int (grp, varid, attname, (int *)att);
        break;
      case NC_UINT:
        status = nc_get_att_uint (grp, varid, attname, (unsigned int *)att);
        break;
      case NC_INT64:
        status = nc_get_att_longlong (grp, varid, attname, (long long *)att);
        break;
      case NC_UINT64:
        status = nc_get_att_ulonglong (grp, varid, attname, (unsigned long long *)att);
        break;
      case NC_FLOAT:
        status = nc_get_att_float (grp, varid, attname, (float *)att);
        break;
      case NC_DOUBLE:
        status = nc_get_att_double (grp, varid, attname, (double *)att);
        break;
      default:
        if (NC_NOERR != (status = nc_get_att (grp, varid, attname, att)))
          {
             _pTIO_err_verror_nc (status, "%s: reading attribute %s in group %d using unrecognized typeid=%d",
                                  __func__, attname, grp, xtype);
             return -1;
          }
        /* drop */
     }

   if (NC_NOERR != status)
     {
        _pTIO_err_verror_nc (status, "%s: reading attribute %s in group %d",
                             __func__, attname, grp);
        return -1;
     }

   return 0;
}

/* Fortran bindings */

FCALLSCFUN6(INT, TIO_get_var_section, TIOF_GET_L1BVAR, tiof_get_l1bvar,
            INT, STRING, INT, INT, INT, PVOID)
