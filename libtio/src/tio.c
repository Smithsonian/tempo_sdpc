/** @file
 *  @brief Core C functions
 */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include <cfortran.h>
#include <netcdf.h>
#include <tell.h>

#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

#define EMPTY()
#define _pTIO_STR(s) #s

int _pTIOMake_Name_Int_Arrays (_pName_Int_Pair_Type *array,
                               int *pnum_values, char **pnames, int **pvalues)
{
   _pName_Int_Pair_Type *p;
   int total_len, num, i;
   int *values = NULL;
   char *names = NULL;
   char *s, *end;

   num = total_len = 0;
   p = array;

   while (1)
     {
        if (p->name[0] == 0)
          break;
        total_len += strlen(p->name);
        num += 1;
        p++;
     }

   /* add room for num-1 spaces, plus 1 trailing null character */
   total_len += num;

   if ((NULL == (values = (int *) TIO_MALLOC (num * sizeof(int))))
       || (NULL == (names = (char *) TIO_MALLOC (total_len * sizeof(char)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        TIO_FREE(values);
        TIO_FREE(names);
        return -1;
     }

   s = names;
   end = names + total_len;

   for (i = 0; i < num; i++)
     {
        int len;
        p = &array[i];

        values[i] = p->value;

        len = strlen(p->name);
        strcpy (s, p->name);
        s += len;

        if (s+1 < end)
          {
             *s = ' ';
             s++;
          }
        else *s = 0;
     }

   *pnum_values = num;
   *pnames = names;
   *pvalues = values;

   return 0;
}

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
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dimension %s (%s)",
                          __func__, o->name, nc_strerror(status));
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
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining enum %s (%s)",
                     __func__, name, nc_strerror(status));
        return -1;
     }

   for (e = enum_table; e->name != NULL; e++)
     {
        status = nc_insert_enum (grp, *enum_typeid, e->name, &e->value);
        if (status != NC_NOERR)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: inserting value %s=%d for enum %s (%s)",
                          __func__, e->name ? e->name : "(null)",
                          e->value, name, nc_strerror(status));
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
        status = nc_put_att_int (grp, varid, a->name, NC_INT, a->num_values, a->value);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining int attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
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
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining float attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
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
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining text attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
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
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining variable %s (%s)",
                     __func__, var_name, nc_strerror(status));
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
   unsigned short fill_ushort = TIO_FILL_USHORT;
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
      case NC_USHORT: pfill_value = &fill_ushort;
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
        Tell_verror (TELL_INVALID_PARM,
                     "%s: invalid fill value type xtype=%d", __func__, xtype);
        return -1;
     }

   status = nc_put_att (grp, varid, _FillValue, xtype, 1, pfill_value);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "writing fill value to grp=%d varid=%d (%s)",
                     grp, varid, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_create (const char *path, int cmode, int *ncid)
{
   int status;

   if ((path == NULL) || (ncid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_create (path, cmode, ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: creating %s (%s)",
                     __func__, path, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_open (const char *path, int omode, int *ncid)
{
   int status;

   if ((path == NULL) || (ncid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_open (path, omode, ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s (%s)",
                     __func__, path, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_close (int ncid)
{
   int status;

   if (NC_NOERR != (status = nc_close (ncid)))
     {
        Tell_verror (TELL_IO_ERROR, "%s: closing file %d (%s)",
                     __func__, ncid, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_inq_var (int grp, const char *name, TIO_Var_Info_Type *info)
{
   int status, i;

   if ((name == NULL) || (info == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_varid (grp, name, &info->varid)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_vartype (grp, info->varid, &info->type)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_varndims (grp, info->varid, &info->ndims)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }

   /* convenient default for the ndims=0 case */
   info->dimlens[0] = 1;

   if (info->ndims == 0)
     return 0;

   if (NC_NOERR != (status = nc_inq_vardimid (grp, info->varid, info->dimids)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s dimids in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }
   for (i = 0; i < info->ndims; i++)
     {
        status = nc_inq_dimlen  (grp, info->dimids[i], &info->dimlens[i]);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: accessing dimension %d length in group %d (%s)",
                          __func__, info->dimids[i], grp, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

/* #define TEST_WAVELENGTH_METHODS 1 */
#ifdef TEST_WAVELENGTH_METHODS
static int get_wavelengths (int grp, int *start, int *count, int type,
                            void *data)
{
   fprintf (stderr, "=====> called get_wavelengths\n");
   return 0;
}
static int put_wavelengths (int grp, int *start, int *count, int type,
                            const void *data)
{
   fprintf (stderr, "=====> called put_wavelengths\n");
   return 0;
}
#endif

static int cvt_float_to_type (int num, const float *f, int type, void *v)
{
   double *dbl = (double *)v;
   int i;

   if (type != NC_DOUBLE)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "%s: Conversion from float to type=%d is not supported",
                     __func__, type);
        return -1;
     }

   for (i = 0; i < num; i++)
     {
        dbl[i] = (double) f[i];
     }

   return 0;
}

static int cvt_type_to_float (int num, float *f, int type, const void *v)
{
   double *dbl = (double *)v;
   int i;

   if (type != NC_DOUBLE)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "%s: Conversion from type=%d to float is not supported",
                     __func__, type);
        return -1;
     }

   for (i = 0; i < num; i++)
     {
        f[i] = (float) dbl[i];
     }

   return 0;
}

typedef struct
{
   char *name;
   int (*get)(int, const int *, const int *, int,       void *);
   int (*put)(int, const int *, const int *, int, const void *);
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

#define TIO_IO_VAR_SECTION(action,error_num,const_qual) \
int TIO_##action##_var_section (int grp, const char *name, \
                                int *istart, int *icount, int type, \
                                const_qual void *data) \
{ \
   TIO_Var_Info_Type info; \
   int status, varid, ndims, i; \
   size_t start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS]; \
   const IO_Methods_Type *io_method; \
 \
   io_method = find_io_methods (name, (const IO_Methods_Type *)&IO_Methods); \
   if (NULL != io_method) \
     { \
        return io_method->action (grp, istart, icount, type, data); \
     } \
 \
   if (-1 == TIO_inq_var (grp, name, &info)) \
     return -1; \
 \
   ndims = info.ndims; \
   varid = info.varid; \
 \
   if (ndims < 0) \
     { \
        Tell_verror (TELL_INVALID_PARM, "%s: variable %s has ndims=%d", \
                     __func__, name, ndims); \
        return -1; \
     } \
 \
   for (i = 0; i < TIO_MAX_VAR_DIMS; i++) \
     { \
        start[i] = 0; \
        count[i] = 0; \
     } \
 \
   for (i = 0; i < ndims; i++) \
     { \
        start[i] = (size_t) istart[i]; \
        count[i] = (icount[i] >= 0) ? (size_t) icount[i] : info.dimlens[i]; \
     } \
 \
   switch (type) \
     { \
      case NC_CHAR: \
        status = nc_##action##_vara_text (grp, varid, start, count, (const_qual char *)data); \
        break; \
      case NC_BYTE: \
        status = nc_##action##_vara_schar (grp, varid, start, count, (const_qual signed char *)data); \
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
      case NC_STRING: \
        status = nc_##action##_vara_string (grp, varid, start, count, (const_qual char **)data); \
        break; \
      default: \
        Tell_verror (TELL_INVALID_PARM, \
                     "%s: accessing variable %s using invalid type (type=%d)", \
                     __func__, name, type); \
        return -1; \
     } \
 \
   if (status != NC_NOERR) \
     { \
        Tell_verror (error_num, \
                     "%s: accessing variable %s in group %d (%s)", \
                     __func__, name, grp, nc_strerror(status)); \
        return -1; \
     } \
 \
   return 0; \
}

TIO_IO_VAR_SECTION(get,TELL_IO_READ_ERROR,EMPTY())
TIO_IO_VAR_SECTION(put,TELL_IO_WRITE_ERROR,const)

#if 0
}
#endif

int TIO_inq_att (int grp, int varid, const char *attname,
                 int *xtype, int *len)
{
   size_t len_size_t;
   int status;

   if (NULL == attname)
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_att (grp, varid, attname, xtype, &len_size_t)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing attribute %s (varid=%d) (%s)",
                     __func__, attname, varid, nc_strerror(status));
        return -1;
     }

   if (len) *len = len_size_t;

   return 0;
}

int TIO_put_att (int grp, int varid, const char *attname,
                 int xtype, int len, const void *att)
{
   int status;
   int file_atttype;
   size_t len_size_t;

   if ((NULL == attname) || (att == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
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
             Tell_verror (TELL_INVALID_PARM,
                          "%s: type mismatch: writing attribute %s: value type=%d, file value type=%d",
                          __func__, attname, xtype, file_atttype);
             return -1;
          }
     }
   else if (NC_ENOTATT != status)
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: attribute %s (varid=%d) query failed (%s)",
                     __func__, attname, varid, nc_strerror(status));
        return -1;
     }

   len_size_t = len;
   status = nc_put_att (grp, varid, attname, xtype, len_size_t, att);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: writing attribute %s (varid=%d) (%s)",
                     __func__, attname, varid, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_get_att (int grp, int varid, const char *attname,
                 int xtype, void *att)
{
   int status;

   if ((NULL == attname) || (att == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* Use the netCDF type-safe interface so that netCDF
    * performs any needed type conversion for built-in types.
    */
   switch (xtype)
     {
      case NC_CHAR:
        status = nc_get_att_text (grp, varid, attname, (char *)att);
        break;
      case NC_BYTE:
        status = nc_get_att_schar (grp, varid, attname, (signed char *)att);
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
      case NC_STRING:
        status = nc_get_att_string (grp, varid, attname, (char **)att);
        break;
      default:
        if (NC_NOERR != (status = nc_get_att (grp, varid, attname, att)))
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: reading attribute %s in group %d using unrecognized typeid=%d (%s)",
                          __func__, attname, grp, xtype, nc_strerror(status));
             return -1;
          }
        /* drop */
     }

   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_READ_ERROR, "%s: reading attribute %s in group %d (%s)",
                     __func__, attname, grp, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_inq_grp (int parent_ncid, const char *path, int *grp)
{
   int status;

   if ((path == NULL) || (grp == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_grp_full_ncid (parent_ncid, path, grp)))
     {
        Tell_verror (TELL_IO_READ_ERROR, "%s: reading group path %s (%s)",
                     __func__, path, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_def_grp (int parent_ncid, const char *path, int *new_ncid)
{
   char delim = '/';
   const char *p;
   unsigned int len;

   if ((path == NULL) || (new_ncid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if ((*path == delim) && (*(path+1) == 0))
     {
        int status = nc_inq_grp_full_ncid (parent_ncid, path, new_ncid);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: accessing id of root group (%s)",
                          __func__, nc_strerror (status));
             return -1;
          }
        return 0;
     }

   for (p = path; *p != 0; p += len)
     {
        char buf[TIO_MAX_NAME_LEN];
        int status, ncid;
        const char *end;

        while (*p == delim) p++;

        end = strchr (p, delim);
        if (end == NULL)
          {
             end = p + strlen(p);
          }
        len = end - p;

        if (len == 0) break;

        if (len + 1 > sizeof(buf))
          {
             Tell_verror (TELL_INVALID_PARM, "%s: group name is too long: %s", __func__, p);
             return -1;
          }
        strncpy (buf, p, len);
        buf[len] = 0;

        /* Does this group exist already? */
        status = nc_inq_grp_ncid (parent_ncid, buf, &ncid);
        if (NC_ENOGRP == status)
          {
             /* If the group doesn't exist, create it */
             status = nc_def_grp (parent_ncid, buf, &ncid);
             if (NC_NOERR != status)
               {
                  Tell_verror (TELL_IO_WRITE_ERROR, "%s: creating group %s (%s)",
                               __func__, buf, nc_strerror (status));
                  return -1;
               }
          }
        else if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR, "%s: accessing group %s (%s)",
                          __func__, buf, nc_strerror (status));
             return -1;
          }

        parent_ncid = ncid;
        *new_ncid = ncid;
     }

   return 0;
}

int TIO_def_dim (int ncid, const char *name, size_t len, int *id)
{
   int status;

   if ((name == NULL) || (id == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_dim (ncid, name, len, id)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dimension %s (%s)",
                     __func__, name, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var (int ncid, const char *name, int type,
                 int num_dims, const int *dimids, int *varid)
{
   int status;

   if ((name == NULL) || (dimids == NULL) || (varid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_var (ncid, name, type, num_dims, dimids, varid)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining variable %s (%s)",
                     __func__, name, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var_fill (int grp, int varid, int no_fill, const void *fill_value)
{
   int status;

   status = nc_def_var_fill (grp, varid, no_fill, fill_value);
   if (status != NC_NOERR)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: setting fill value for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_inq_var_fill (int grp, int varid, int *no_fill, void *fill_value)
{
   int status;

   status = nc_inq_var_fill (grp, varid, no_fill, fill_value);
   if (status != NC_NOERR)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: getting fill value for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var_deflate (int grp, int varid,
                         int shuffle, int deflate, int deflate_level)
{
   int status;

   status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level);
   if (status != NC_NOERR)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: setting deflate values for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_put_text_attrs (int grp, int varid, const TIO_Attr_Text_Type *attrs)
{
   const TIO_Attr_Text_Type *a;
   int status;

   if (attrs == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_text (grp, varid, a->name, strlen(a->text), a->text);
        if (status != NC_NOERR)
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing text attribute %s (%s)",
                          __func__, a->name, nc_strerror (status));
             return -1;
          }
     }

   return 0;
}

int TIO_copy_attrs (int ncid_infile, int id_var_infile,
                    int (*dontcopy_attr)(const char *),
                    int ncid, int id_var)
{
   char attname[TIO_MAX_NAME_LEN];
   int status, attnum, num_atts;

   status = nc_inq_varnatts (ncid_infile, id_var_infile, &num_atts);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: reading number of attributes for variable id=%d (%s)",
                     __func__, id_var_infile, nc_strerror (status));
        return -1;
     }

   for (attnum = 0; attnum < num_atts; attnum++)
     {
        status = nc_inq_attname (ncid_infile, id_var_infile, attnum, attname);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: reading attribute %d=%s (%s)",
                          __func__, attnum, attname,
                          nc_strerror (status));
             return -1;
          }
        if ((NULL != dontcopy_attr) && (1 == dontcopy_attr (attname)))
          continue;
        status = nc_copy_att (ncid_infile, id_var_infile, attname,
                              ncid, id_var);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: copying attribute %s (%s)",
                          __func__, attname, nc_strerror (status));
             return -1;
          }
     }

   return 0;
}

int TIO_inq_dimname (int grp, int dimid, char *dimname)
{
   int status;

   if (dimname == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_dimname (grp, dimid, dimname)))
     {
        Tell_verror (TELL_IO_READ_ERROR, "%s: reading dimension %d name (%s)",
                     __func__, dimid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_inq_dimid (int grp, const char *dimname, int *dimid)
{
   if ((dimname == NULL) || (dimid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* For such a low-level operation, it seems better to
    * _not_ generate an error message by default */
   return nc_inq_dimid (grp, dimname, dimid);
}

int TIO_inq_dim (int grp, const char *dimname, int *dimid, size_t *dimlen)
{
   int status;
   if ((dimname == NULL) || (dimid == NULL) || (dimlen == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* For such a low-level operation, it seems better to
    * _not_ generate an error message by default */

   status = nc_inq_dimid (grp, dimname, dimid);
   if (NC_NOERR != status)
     return status;

   return nc_inq_dimlen (grp, *dimid, dimlen);
}

int TIO_put_git_commit_hash (int grp, const char *attname)
{
   const char hash[] = GIT_COMMIT_HASH ;
   const char *name = "tio_commit";
   int status;

   /* Use a default attribute name when the supplied attribute
    * name is either NULL, empty, or begins with a space.
    */
   if ((attname != NULL) && (*attname != 0) && (*attname != ' '))
     name = attname;

   status = nc_put_att_text (grp, NC_GLOBAL, name, strlen(hash), hash);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: writing git commit hash (%s)",
                     __func__, nc_strerror(status));
        return -1;
     }

   return 0;
}

#define CONVERT_TO(to_type_name,to_type) \
static int convert_to_##to_type (int from_type, void *from_value, to_type_name *to_value) \
{ \
   to_type_name value; \
 \
   switch (from_type) \
     { \
      case NC_BYTE: value = *(signed char *)from_value; break; \
      case NC_SHORT: value = *(short *)from_value; break; \
      case NC_INT: value = *(int *)from_value; break; \
      case NC_INT64: value = *(long long *)from_value; break; \
      case NC_UBYTE: value = *(unsigned char *)from_value; break; \
      case NC_USHORT: value = *(unsigned short *)from_value; break; \
      case NC_UINT: value = *(unsigned int *)from_value; break; \
      case NC_UINT64: value = *(unsigned long long *)from_value; break; \
      case NC_FLOAT: value = *(float *)from_value; break; \
      case NC_DOUBLE: value = *(double *)from_value; break; \
      default: \
        Tell_verror (TELL_INVALID_PARM, "%s: unsupported type from_type=%d\n", \
                     __func__, from_type); \
        return -1; \
     } \
 \
   *to_value = value; \
 \
   return 0; \
}
CONVERT_TO(signed char,byte)
CONVERT_TO(unsigned char,ubyte)
CONVERT_TO(short,short)
CONVERT_TO(unsigned short,ushort)
CONVERT_TO(int,int)
CONVERT_TO(unsigned int,uint)
CONVERT_TO(long long,int64)
CONVERT_TO(unsigned long long,uint64)
CONVERT_TO(float,float)
CONVERT_TO(double,double)

int TIO_get_fill_value (int grp, const char *name, int type, void *value)
{
   int varid, file_type, no_fill, status, return_status = -1;
   void *fill = NULL;

   if ((name == NULL) || (value == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NULL == (fill = TIO_MALLOC (8)))
     {
        Tell_set_error (TELL_MALLOC_ERROR);
        return -1;
     }

   if ((NC_NOERR != (status = nc_inq_varid (grp, name, &varid)))
       || (NC_NOERR != (status = nc_inq_vartype (grp, varid, &file_type)))
       || (NC_NOERR != (status = nc_inq_var_fill (grp, varid, &no_fill, fill))))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: getting fill value for %s (%s)\n",
                     __func__, name, nc_strerror(status));
        goto cleanup;
     }

   switch (type)
     {
      case NC_BYTE:
        status = convert_to_byte (file_type, fill, (signed char *)value);
        break;
      case NC_UBYTE:
        status = convert_to_ubyte (file_type, fill, (unsigned char *)value);
        break;
      case NC_SHORT:
        status = convert_to_short (file_type, fill, (short *)value);
        break;
      case NC_USHORT:
        status = convert_to_ushort (file_type, fill, (unsigned short *)value);
        break;
      case NC_INT:
        status = convert_to_int (file_type, fill, (int *)value);
        break;
      case NC_UINT:
        status = convert_to_uint (file_type, fill, (unsigned int *)value);
        break;
      case NC_INT64:
        status = convert_to_int64 (file_type, fill, (long long *)value);
        break;
      case NC_UINT64:
        status = convert_to_uint64 (file_type, fill, (unsigned long long *)value);
        break;
      case NC_FLOAT:
        status = convert_to_float (file_type, fill, (float *)value);
        break;
      case NC_DOUBLE:
        status = convert_to_double (file_type, fill, (double *)value);
        break;
      default:
        Tell_verror (TELL_INVALID_PARM, "%s: unsupported destination type = %d\n",
                     __func__, type);
        goto cleanup;
     }

   if (status)
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: converting %s fill value to type = %d\n",
                     __func__, name, type);
        goto cleanup;
     }

   return_status = 0;
cleanup:
   free (fill);
   return return_status;
}

/* Fortran bindings */

FCALLSCFUN6(INT, TIO_get_var_section, TIOF_GET_VAR_SECTION, tiof_get_var_section,
            INT, STRING, PINT, PINT, INT, PVOID)
FCALLSCFUN6(INT, TIO_put_var_section, TIOF_PUT_VAR_SECTION, tiof_put_var_section,
            INT, STRING, PINT, PINT, INT, PVOID)
FCALLSCFUN2(INT, TIO_put_git_commit_hash, TIO_F_PUT_GIT_HASH, tio_f_put_git_hash,
            INT, STRING)
FCALLSCFUN3(INT, TIO_def_grp, TIO_F_DEF_GRP, tio_f_def_grp,
            INT, STRING, PINT)
FCALLSCFUN4(INT, TIO_get_fill_value, TIO_F_GET_FILL_VALUE, tio_f_get_fill_value,
            INT, STRING, INT, PVOID)
