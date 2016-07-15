#include "defs.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tell.h>
#include <netcdf.h>
#include <tio.h>
#include <tio_template.h>

#include "poly.h"
#include "pixel.h"
#include "regrid.h"
#include "var.h"

typedef union
{
   double *d;
   unsigned long long *ul;
   unsigned int   *ui;
   unsigned short *us;
   unsigned char  *uc;
   long long *l;
   int   *i;
   short *s;
   char  *c;
}
Value_Ptr_Type;

typedef union
{
   double d;
   unsigned long long ul;
   unsigned int   ui;
   unsigned short us;
   unsigned char  uc;
   long long l;
   int   i;
   short s;
   char  c;
}
Fill_Value_Type;

struct Var_Value_Buffer_Type
{
   Value_Ptr_Type src_values;
   Value_Ptr_Type dest_values;
   Fill_Value_Type fill_value;
   int value_type;
   int *src_mask;
   int num_dims, dimlens[TIO_MAX_VAR_DIMS];
   int num_step;                 /* Number of mirror steps in this granule */
   int num_xtrack;               /* Number of pixels along slit */
   int num_src_pixels;           /* number of spatial pixels in this granule */
   int num_dest_pixels;          /* number of lon/lat pixels in destination grid */
   int num_values_per_pixel;
   int bytes_per_value;
};

static int value_num_bytes (int value_type)
{
   switch (value_type)
     {
      case VALUE_IS_DOUBLE:
      case VALUE_IS_UINT64:
      case VALUE_IS_INT64:
        return 8;
      case VALUE_IS_UINT:
      case VALUE_IS_INT:
        return 4;
      case VALUE_IS_USHORT:
      case VALUE_IS_SHORT:
        return 2;
      case VALUE_IS_UBYTE:
      case VALUE_IS_BYTE:
        return 1;
      default:
        Tell_verror (TELL_APPLICATION_ERROR, "%s: invalid value_type=%d",
                     __func__, value_type);
        break;
     }

   return -1;
}

static int value_io_type (int value_type)
{
   switch (value_type)
     {
      case VALUE_IS_DOUBLE: return TIO_DOUBLE;
      case VALUE_IS_UINT64: return TIO_UINT64;
      case VALUE_IS_UINT:   return TIO_UINT;
      case VALUE_IS_USHORT: return TIO_USHORT;
      case VALUE_IS_UBYTE:  return TIO_UBYTE;
      case VALUE_IS_INT64:  return TIO_INT64;
      case VALUE_IS_INT:    return TIO_INT;
      case VALUE_IS_SHORT:  return TIO_SHORT;
      case VALUE_IS_BYTE:   return TIO_BYTE;
      default:
        Tell_verror (TELL_APPLICATION_ERROR, "%s: invalid value_type=%d",
                     __func__, value_type);
        break;
     }
   return -1;
}

static void value_default_fill (Var_Value_Buffer_Type *vb, int value_type)
{
   /* zero all the union bytes first */
   vb->fill_value.ul = 0;

   switch (value_type)
     {
      case NC_DOUBLE: vb->fill_value.d  = NC_FILL_DOUBLE; break;
      case NC_FLOAT:  vb->fill_value.d  = NC_FILL_FLOAT; break;
      case NC_UINT64: vb->fill_value.ul = NC_FILL_UINT64; break;
      case NC_UINT:   vb->fill_value.ui = NC_FILL_UINT; break;
      case NC_USHORT: vb->fill_value.us = NC_FILL_USHORT; break;
      case NC_UBYTE:  vb->fill_value.uc = NC_FILL_UBYTE; break;
      case NC_INT64:  vb->fill_value.l  = NC_FILL_INT64; break;
      case NC_INT:    vb->fill_value.i  = NC_FILL_INT; break;
      case NC_SHORT:  vb->fill_value.s  = NC_FILL_SHORT; break;
      case NC_BYTE:   vb->fill_value.c  = NC_FILL_BYTE; break;
     }
}

void Var_free_value_buffer (Var_Value_Buffer_Type *vb)
{
   if (vb == NULL)
     return;
   FREE(vb->src_values.d);
   FREE(vb->dest_values.d);
   FREE(vb->src_mask);
   FREE(vb);
}

Var_Value_Buffer_Type *
Var_new_value_buffer (int dest_nx, int dest_ny,
                      int src_num_steps, int src_num_xtrack)
{
   Var_Value_Buffer_Type *vb = NULL;
   int len_src, len_dest, len_mask;

   if (NULL == (vb = (Var_Value_Buffer_Type *)MALLOC (sizeof *vb)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   vb->src_mask = NULL;
   vb->src_values.d = NULL;
   vb->dest_values.d = NULL;

   vb->num_step = src_num_steps;
   vb->num_xtrack = src_num_xtrack;

   vb->num_src_pixels = vb->num_step * vb->num_xtrack;
   vb->num_dest_pixels = dest_nx * dest_ny;

   /* The most common level 2 variable to regrid is
    * a 2D array with one value per spatial pixel */
   vb->num_values_per_pixel = 1;
   vb->num_dims = 2;
   vb->dimlens[0] = src_num_steps;
   vb->dimlens[1] = src_num_xtrack;

   len_mask = vb->num_src_pixels * sizeof(int);
   vb->bytes_per_value = sizeof(double);
   len_src = vb->num_src_pixels * vb->bytes_per_value;
   len_dest = vb->num_dest_pixels * vb->bytes_per_value;

   /* By default, treat input values as numbers to be averaged,
    * and handle them internally as doubles.  Bit-fields will
    * be handled differently. */
   vb->value_type = VALUE_IS_DOUBLE;

   if ((NULL == (vb->src_values.d = (double *)MALLOC (len_src)))
       || (NULL == (vb->dest_values.d = (double *)MALLOC (len_dest)))
       || (NULL == (vb->src_mask = (int *) MALLOC (len_mask))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        Var_free_value_buffer (vb);
        return NULL;
     }

   memset ((char *)vb->src_mask, 0, len_mask);

   return vb;
}

static int maybe_realloc_value_buf (int ncid, Var_Value_Buffer_Type *vb,
                                    const char *var_name)
{
   TIO_Var_Info_Type vi;
   void *tmp;
   int bytes_per_value, have_num_src_bytes, have_num_dest_bytes;
   int i, num_src_values, num_dest_values;

   if (-1 == TIO_inq_var (ncid, var_name, &vi))
     return -1;

   have_num_src_bytes = ((vb->num_src_pixels * vb->num_values_per_pixel)
                         * vb->bytes_per_value);
   have_num_dest_bytes = ((vb->num_dest_pixels * vb->num_values_per_pixel)
                          * vb->bytes_per_value);

   vb->num_dims = vi.ndims;
   for (i = 0; i < vi.ndims; i++)
     {
        vb->dimlens[i] = vi.dimlens[i];
     }
   vb->num_values_per_pixel = 1;
   for (i = 2; i < vi.ndims; i++)
     {
        vb->num_values_per_pixel *= vi.dimlens[i];
     }

   num_src_values = vb->num_src_pixels * vb->num_values_per_pixel;
   num_dest_values = vb->num_dest_pixels * vb->num_values_per_pixel;

   if (-1 == (bytes_per_value = value_num_bytes (vb->value_type)))
     return -1;

   if ((num_src_values * bytes_per_value < have_num_src_bytes)
       && (num_dest_values * bytes_per_value < have_num_dest_bytes))
     return 0;

   if (NULL == (tmp = REALLOC (vb->src_values.d, num_src_values * bytes_per_value)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }
   vb->src_values.d = (double *)tmp;

   if (NULL == (tmp = REALLOC (vb->dest_values.d, num_dest_values * bytes_per_value)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }
   vb->dest_values.d = (double *)tmp;

   vb->bytes_per_value = bytes_per_value;

   return 0;
}

static int read_var_values1 (int ncid, int var_grp, Var_Value_Buffer_Type *vb,
                            const char *var_name)
{
   TIO_Var_Info_Type vi;
   int start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   int i, num_steps, num_pixels, num_values, bytes_per_value, in_type;
   int *step = NULL;
   unsigned char *var = NULL;
   int status = -1;

   num_steps = vb->dimlens[0];

   if (NULL == (step = (int *) MALLOC (num_steps * sizeof (int))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_return;
     }

   start[0] = 0;
   count[0] = num_steps;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_STEP,
                                  start, count, TIO_INT, step))
     goto cleanup_and_return;

   /* Set a default fill-value then let any fill-value in the file
    * override it.  If there's no fill-value in the file, the
    * default value won't be over-written. */
   if (-1 == TIO_inq_var (var_grp, var_name, &vi))
     goto cleanup_and_return;
   value_default_fill (vb, vi.type);
   if (vb->value_type == VALUE_IS_DOUBLE)
     {
        if (-1 == TIO_get_fill_value (var_grp, var_name, TIO_DOUBLE, &vb->fill_value.d))
          goto cleanup_and_return;
     }
   else
     {
        if (-1 == TIO_inq_var_fill (var_grp, vi.varid, NULL, &vb->fill_value.d))
          goto cleanup_and_return;
     }

   num_pixels = num_steps * vb->num_xtrack;

   for (i = 0; i < vb->num_dims; i++)
     {
        start[i] = 0;
     }

   count[0] = num_steps;
   count[1] = vb->num_xtrack;
   for (i = 2; i < vb->num_dims; i++)
     {
        count[i] = vb->dimlens[i];
     }

   if ((-1 == (bytes_per_value = value_num_bytes (vb->value_type)))
       || (-1 == (in_type = value_io_type (vb->value_type))))
     goto cleanup_and_return;

   num_values = num_pixels * vb->num_values_per_pixel;
   if (NULL == (var = (unsigned char *) MALLOC (num_values * bytes_per_value)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_return;
     }

   if (-1 == TIO_get_var_section (var_grp, var_name,
                                  start, count, in_type, var))
     goto cleanup_and_return;

   for (i = 0; i < num_pixels; i++)
     {
        int pix_xtrack = i % vb->num_xtrack;
        int pix_step_index = i / vb->num_xtrack;
        int pix = pix_xtrack + step[pix_step_index] * vb->num_xtrack;
        int bpp = bytes_per_value * vb->num_values_per_pixel;
        memcpy (vb->src_values.uc + pix * bpp, var + i * bpp, bpp);
     }

   status = 0;
cleanup_and_return:
   FREE(step);
   FREE(var);

   return status;
}

/* parse_var_path must be able to write to var_path */
static int parse_var_path (char *var_path,
                           char **grp_path, char **var_name)
{
   char *p;

   if (NULL == (p = strrchr (var_path, '/')))
     {
        *grp_path = NULL;
        *var_name = var_path;
     }
   else
     {
        *p = 0;
        *var_name = p + 1;
        *grp_path = var_path;
     }

   return 0;
}

static int read_var_values (Var_Value_Buffer_Type *vb, const char *var_path,
                            char **files, int num_files)
{
   char *var_path_copy = NULL;
   char *grp_path = NULL;
   char *var_name = NULL;
   int i, len, ncid, status = -1;

   len = strlen (var_path) + 1;
   if (NULL == (var_path_copy = (char *) MALLOC (len * sizeof(char))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memcpy (var_path_copy, var_path, len);

   if (-1 == parse_var_path (var_path_copy, &grp_path, &var_name))
     goto free_and_return;

   for (i = 0; i < num_files; i++)
     {
        int grp;

        if (-1 == TIO_open (files[i], NC_NOWRITE, &ncid))
          goto free_and_return;

        if (grp_path)
          {
             if (-1 == TIO_inq_grp (ncid, grp_path, &grp))
               goto free_and_return;
          }
        else grp = ncid;

        if (i == 0)
          {
             if (-1 == maybe_realloc_value_buf (grp, vb, var_name))
               goto free_and_return;
          }

        if (-1 == read_var_values1 (ncid, grp, vb, var_name))
          goto free_and_return;
     }

   status = 0;
free_and_return:
   (void) TIO_close (ncid);
   FREE(var_path_copy);

   return status;
}

#define COPY_FROM_STRIDED(typestr, type) \
static void copy_from_strided_##typestr (int num, int stride, \
                                         type fill_value, int *src_mask, \
                                         const type *strided_values, \
                                         type *packed_values) \
{ \
   int k; \
   if (strided_values == packed_values) \
     return; \
   for (k = 0; k < num; k++) \
     { \
        packed_values[k] = strided_values[k*stride]; \
        src_mask[k] = (packed_values[k] == fill_value); \
     } \
}

#define COPY_TO_STRIDED(typestr, type) \
static void copy_to_strided_##typestr (int num, int stride, \
                                       const type *packed_values, \
                                       type *strided_values) \
{ \
   int k; \
   if (strided_values == packed_values) \
     return; \
   for (k = 0; k < num; k++) \
     { \
        strided_values[k*stride] = packed_values[k]; \
     } \
}

COPY_FROM_STRIDED(double,double)
COPY_FROM_STRIDED(uint64,unsigned long long)
COPY_FROM_STRIDED(uint,  unsigned int)
COPY_FROM_STRIDED(ushort,unsigned short)
COPY_FROM_STRIDED(ubyte, unsigned char)
COPY_FROM_STRIDED(int64,long long)
COPY_FROM_STRIDED(int,  int)
COPY_FROM_STRIDED(short,short)
COPY_FROM_STRIDED(byte, char)

COPY_TO_STRIDED(double,double)
COPY_TO_STRIDED(uint64,unsigned long long)
COPY_TO_STRIDED(uint,  unsigned int)
COPY_TO_STRIDED(ushort,unsigned short)
COPY_TO_STRIDED(ubyte, unsigned char)
COPY_TO_STRIDED(int64,long long)
COPY_TO_STRIDED(int,  int)
COPY_TO_STRIDED(short,short)
COPY_TO_STRIDED(byte, char)

static void copy_from_strided_src (Var_Value_Buffer_Type *vb, int i,
                                   Value_Ptr_Type src_values)
{
   switch (vb->value_type)
     {
      case VALUE_IS_DOUBLE:
        copy_from_strided_double (vb->num_src_pixels,
                                  vb->num_values_per_pixel,
                                  vb->fill_value.d, vb->src_mask,
                                  vb->src_values.d + i, src_values.d);
        break;
      case VALUE_IS_UINT64:
        copy_from_strided_uint64 (vb->num_src_pixels,
                                  vb->num_values_per_pixel,
                                  vb->fill_value.ul, vb->src_mask,
                                  vb->src_values.ul + i, src_values.ul);
        break;
      case VALUE_IS_UINT:
        copy_from_strided_uint   (vb->num_src_pixels,
                                  vb->num_values_per_pixel,
                                  vb->fill_value.ui, vb->src_mask,
                                  vb->src_values.ui + i, src_values.ui);
        break;
      case VALUE_IS_USHORT:
        copy_from_strided_ushort (vb->num_src_pixels,
                                  vb->num_values_per_pixel,
                                  vb->fill_value.us, vb->src_mask,
                                  vb->src_values.us + i, src_values.us);
        break;
      case VALUE_IS_UBYTE:
        copy_from_strided_ubyte  (vb->num_src_pixels,
                                  vb->num_values_per_pixel,
                                  vb->fill_value.uc, vb->src_mask,
                                  vb->src_values.uc + i, src_values.uc);
        break;
      case VALUE_IS_INT64:
        copy_from_strided_int64 (vb->num_src_pixels,
                                 vb->num_values_per_pixel,
                                 vb->fill_value.l, vb->src_mask,
                                 vb->src_values.l + i, src_values.l);
        break;
      case VALUE_IS_INT:
        copy_from_strided_int   (vb->num_src_pixels,
                                 vb->num_values_per_pixel,
                                 vb->fill_value.i, vb->src_mask,
                                 vb->src_values.i + i, src_values.i);
        break;
      case VALUE_IS_SHORT:
        copy_from_strided_short (vb->num_src_pixels,
                                 vb->num_values_per_pixel,
                                 vb->fill_value.s, vb->src_mask,
                                 vb->src_values.s + i, src_values.s);
        break;
      case VALUE_IS_BYTE:
        copy_from_strided_byte  (vb->num_src_pixels,
                                 vb->num_values_per_pixel,
                                 vb->fill_value.c, vb->src_mask,
                                 vb->src_values.c + i, src_values.c);
        break;
     }
}

static void copy_to_strided_dest (Var_Value_Buffer_Type *vb, int i, Value_Ptr_Type dest_values)
{
   switch (vb->value_type)
     {
      case VALUE_IS_DOUBLE:
        copy_to_strided_double (vb->num_dest_pixels,
                                vb->num_values_per_pixel,
                                dest_values.d, vb->dest_values.d + i);
        break;
      case VALUE_IS_UINT64:
        copy_to_strided_uint64 (vb->num_dest_pixels,
                                vb->num_values_per_pixel,
                                dest_values.ul, vb->dest_values.ul + i);
        break;
      case VALUE_IS_UINT:
        copy_to_strided_uint   (vb->num_dest_pixels,
                                vb->num_values_per_pixel,
                                dest_values.ui, vb->dest_values.ui + i);
        break;
      case VALUE_IS_USHORT:
        copy_to_strided_ushort (vb->num_dest_pixels,
                                vb->num_values_per_pixel,
                                dest_values.us, vb->dest_values.us + i);
        break;
      case VALUE_IS_UBYTE:
        copy_to_strided_ubyte  (vb->num_dest_pixels,
                                vb->num_values_per_pixel,
                                dest_values.uc, vb->dest_values.uc + i);
        break;
      case VALUE_IS_INT64:
        copy_to_strided_int64 (vb->num_dest_pixels,
                               vb->num_values_per_pixel,
                               dest_values.l, vb->dest_values.l + i);
        break;
      case VALUE_IS_INT:
        copy_to_strided_int   (vb->num_dest_pixels,
                               vb->num_values_per_pixel,
                               dest_values.i, vb->dest_values.i + i);
        break;
      case VALUE_IS_SHORT:
        copy_to_strided_short (vb->num_dest_pixels,
                               vb->num_values_per_pixel,
                               dest_values.s, vb->dest_values.s + i);
        break;
      case VALUE_IS_BYTE:
        copy_to_strided_byte  (vb->num_dest_pixels,
                               vb->num_values_per_pixel,
                               dest_values.c, vb->dest_values.c + i);
        break;
     }
}

int Var_apply_regrid (const Pixel_Regrid_Type *r,
                      Var_Value_Buffer_Type *vb,
                      int value_type, const char *var_path,
                      char **files, int num_files)
{
   Value_Ptr_Type src_values, dest_values;
   int i, regrid_by_averaging, status = -1;

   vb->value_type = value_type;

   src_values.d = NULL;
   dest_values.d = NULL;

   if (-1 == read_var_values (vb, var_path, files, num_files))
     return -1;

   /* Regridded result will be returned in vb->dest_values.
    * For 2D data, this is straightforward but, for
    * higher-dimensional variables, we regrid one slice
    * (dest->nx * dest->ny values) at a time.
    */
   if (vb->num_values_per_pixel == 1)
     {
        src_values.d = vb->src_values.d;
        dest_values.d = vb->dest_values.d;
     }
   else
     {
        int len_src = vb->num_src_pixels * vb->bytes_per_value;
        int len_dest = vb->num_dest_pixels * vb->bytes_per_value;
        if ((NULL == (src_values.uc = (unsigned char *)MALLOC (len_src)))
            || (NULL == (dest_values.uc = (unsigned char *)MALLOC (len_dest))))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto free_and_return;
          }
     }

   regrid_by_averaging = (vb->value_type == VALUE_IS_DOUBLE);

   for (i = 0; i < vb->num_values_per_pixel; i++)
     {
        int rstatus;
        copy_from_strided_src (vb, i, src_values);
        if (regrid_by_averaging)
          {
             rstatus = Pixel_regrid (r, vb->src_mask, vb->fill_value.d,
                                     src_values.d, dest_values.d);
          }
        else
          {
             rstatus = Pixel_regrid_bytes (r, vb->src_mask, vb->value_type,
                                           &vb->fill_value.uc,
                                           src_values.uc, dest_values.uc);
          }
        if (rstatus)
          break;
        copy_to_strided_dest (vb, i, dest_values);
     }

   if (i == vb->num_values_per_pixel)
     {
        status = 0;
     }

free_and_return:

   if (vb->num_values_per_pixel > 1)
     {
        FREE(src_values.uc);
        FREE(dest_values.uc);
     }

   return status;
}

static int write_grid_1d (int ncid, int grp, const char *name,
                          double xmin, double xmax, int num,
                          const TIO_Attr_Text_Type *attrs)
{
   int i, dim_id, var_id, start, count;
   double *x = NULL;
   double dx;

   /* assume dimensions are global even when lon-lat variables
    * are in a group */
   if (-1 == TIO_def_dim (ncid, name, num, &dim_id))
     return -1;

   if ((-1 == TIO_def_var (grp, name, NC_FLOAT, 1, &dim_id, &var_id))
       || (-1 == TIO_put_text_attrs (grp, dim_id, attrs)))
     return -1;

   if (NULL == (x = (double *)MALLOC (num * sizeof(double))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   dx = (xmax - xmin) / num;
   for (i = 0; i < num; i++)
     {
        x[i] = xmin + (i + 0.5) * dx;
     }

   start = 0;
   count = num;
   if (-1 == TIO_put_var_section (grp, name, &start, &count, TIO_DOUBLE, x))
     {
        FREE(x);
        return -1;
     }

   FREE(x);
   return 0;
}

int Var_write_lonlat_grid (int ncid, const char *lonlat_grp,
                           const Pixel_Grid_Param_Type *dest)
{
   static TIO_Attr_Text_Type lon_attrs[] =
     {
        {"units", "degrees_east"},
        {NULL,NULL}
     };
   static TIO_Attr_Text_Type lat_attrs[] =
     {
        {"units", "degrees_north"},
        {NULL,NULL}
     };
   int grp = ncid;

   if (lonlat_grp)
     {
        if (-1 == TIO_def_grp (ncid, lonlat_grp, &grp))
          return -1;
     }

   if (-1 == write_grid_1d (ncid, grp, TEMPO_VAR_LONGITUDE,
                             dest->xmin, dest->xmax, dest->nx, lon_attrs))
     return -1;

   if (-1 == write_grid_1d (ncid, grp, TEMPO_VAR_LATITUDE,
                             dest->ymin, dest->ymax, dest->ny, lat_attrs))
     return -1;

   return 0;
}

static int dontcopy_attr (const char *attname)
{
   const char *lst[] = {
      "bounds"
   };
   int i, n = sizeof(lst)/sizeof(*lst);

   for (i = 0; i < n; i++)
     {
        if (0 == strcmp (lst[i], attname))
          return 1;
     }

   return 0;
}

static int copy_extra_dims (int ncid_infile, const TIO_Var_Info_Type *vi,
                            int ncid, int *dims)
{
   char dimname[TIO_MAX_NAME_LEN];
   int i;
   for (i = 2; i < vi->ndims; i++)
     {
        if (-1 == TIO_inq_dimname (ncid_infile, vi->dimids[i], dimname))
          return -1;
        /* If this dimension already exists in the output file,
         * record the dimid and continue.  Otherwise, create it. */
        if (0 == TIO_inq_dimid (ncid, dimname, &dims[i]))
          continue;
        if (-1 == TIO_def_dim (ncid, dimname, vi->dimlens[i], &dims[i]))
          return -1;
     }

   return 0;
}

int Var_write_values (int ncid, const Var_Value_Buffer_Type *vb,
                      const char *out_var_path,
                      int ncid_infile, const char *in_var_path)
{
   TIO_Var_Info_Type vi;
   char *in_var_path_copy=NULL, *out_var_path_copy=NULL;
   char *in_var_name=NULL, *out_var_name=NULL;
   char *in_grp_path = NULL, *out_grp_path = NULL;
   int start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   size_t lon_dimlen, lat_dimlen;
   int lon_dimid, lat_dimid, len_in_var_path, len_out_var_path;
   int i, dims[TIO_MAX_VAR_DIMS];
   int in_grp, in_varid, out_grp, out_varid, out_type;
   int in_no_fill, shuffle=1, deflate=1, deflate_level=1;
   int status = -1;

   /* To facilitate parsing group/var names, we first make
    * a copy of each path string that we can safely write to */
   len_in_var_path = strlen(in_var_path) + 1;
   if (NULL == (in_var_path_copy = (char *) MALLOC (len_in_var_path)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memcpy (in_var_path_copy, in_var_path, len_in_var_path);
   if (-1 == parse_var_path (in_var_path_copy, &in_grp_path, &in_var_name))
     goto free_and_return;

   len_out_var_path = strlen(out_var_path) + 1;
   if (NULL == (out_var_path_copy = (char *) MALLOC (len_out_var_path)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }
   memcpy (out_var_path_copy, out_var_path, len_out_var_path);
   if (-1 == parse_var_path (out_var_path_copy, &out_grp_path, &out_var_name))
     goto free_and_return;

   if (in_grp_path)
     {
        if (-1 == TIO_inq_grp (ncid_infile, in_grp_path, &in_grp))
          goto free_and_return;
     }
   else in_grp = ncid_infile;

   if (out_grp_path)
     {
        if (-1 == TIO_def_grp (ncid, out_grp_path, &out_grp))
          goto free_and_return;
     }
   else out_grp = ncid;

   if (-1 == TIO_inq_var (in_grp, in_var_name, &vi))
     goto free_and_return;

   in_varid = vi.varid;

   if (vi.ndims > 2)
     {
        if (-1 == copy_extra_dims (ncid_infile, &vi, ncid, dims))
          goto free_and_return;
     }

   if ((0 != TIO_inq_dim (ncid, TEMPO_VAR_LONGITUDE, &lon_dimid, &lon_dimlen))
       || (0 != TIO_inq_dim (ncid, TEMPO_VAR_LATITUDE, &lat_dimid, &lat_dimlen)))
     return -1;

   dims[0] = lat_dimid;
   dims[1] = lon_dimid;

   if ((-1 == TIO_def_var (out_grp, out_var_name, vi.type, vb->num_dims, dims, &out_varid))
       || (-1 == TIO_inq_var_fill (in_grp, in_varid, &in_no_fill, NULL))
       || (-1 == TIO_def_var_fill (out_grp, out_varid, in_no_fill, &vb->fill_value.d))
       || (-1 == TIO_def_var_deflate (out_grp, out_varid, shuffle, deflate, deflate_level))
       || (-1 == TIO_copy_attrs (in_grp, in_varid, dontcopy_attr, out_grp, out_varid)))
     {
        goto free_and_return;
     }

   for (i = 0; i < vb->num_dims; i++)
     {
        start[i] = 0;
     }
   count[0] = lat_dimlen;
   count[1] = lon_dimlen;
   for (i = 2; i < vb->num_dims; i++)
     {
        count[i] = vb->dimlens[i];
     }
   if ((-1 == (out_type = value_io_type (vb->value_type)))
       || (-1 == TIO_put_var_section (out_grp, out_var_name, start, count,
                                      out_type, vb->dest_values.d)))
     goto free_and_return;

   status = 0;
free_and_return:
   FREE(in_var_path_copy);
   FREE(out_var_path_copy);

   return status;
}
