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

#ifndef REALLOC
# define REALLOC realloc
#endif

#ifndef MALLOC
# define MALLOC malloc
#endif

#ifndef FREE
# define FREE free
#endif

typedef struct
{
   double *src_values;
   double *dest_values;
   int *src_mask;
   int num_dims, dimlens[TIO_MAX_VAR_DIMS];
   int num_step;                 /* Number of mirror steps in this granule */
   int num_xtrack;               /* Number of pixels along slit */
   int num_src_pixels;           /* number of spatial pixels in this granule */
   int num_dest_pixels;          /* number of lon/lat pixels in destination grid */
   int num_values_per_pixel;
   int num_alloc_values;
   Pixel_Overlap_Info_Type *overlap_info;
}
Var_Value_Buffer_Type;

void Var_free_value_buffer (Var_Value_Buffer_Type *vb)
{
   if (vb == NULL)
     return;
   FREE(vb->src_values);
   FREE(vb->src_mask);
   FREE(vb->overlap_info);
   FREE(vb);
}

Var_Value_Buffer_Type *
Var_new_value_buffer (int dest_nx, int dest_ny,
                      int src_num_step, int src_num_xtrack)
{
   Var_Value_Buffer_Type *vb = NULL;
   int len, len_mask;

   if (NULL == (vb = (Var_Value_Buffer_Type *)MALLOC (sizeof *vb)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   vb->src_values = NULL;
   vb->src_mask = NULL;

   vb->num_step = src_num_step;
   vb->num_xtrack = src_num_xtrack;

   vb->num_src_pixels = vb->num_step * vb->num_xtrack;
   vb->num_dest_pixels = dest_nx * dest_ny;

   /* The most common level 2 variable to regrid is
    * a 2D array with one value per spatial pixel */
   vb->num_values_per_pixel = 1;
   vb->num_dims = 2;
   vb->dimlens[0] = src_num_step;
   vb->dimlens[1] = src_num_xtrack;

   /* Note that vb->src_values and vb->dest_values
    * share a single malloced space */
   vb->num_alloc_values = vb->num_src_pixels + vb->num_dest_pixels;
   len = vb->num_alloc_values * sizeof(double);
   len_mask = vb->num_src_pixels * sizeof(int);

   if ((NULL == (vb->src_values = (double *)MALLOC (len)))
       || (NULL == (vb->src_mask = (int *) MALLOC (len_mask))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        Var_free_value_buffer (vb);
        return NULL;
     }
   vb->dest_values = vb->src_values + vb->num_src_pixels;

   memset ((char *)vb->src_mask, 0, len_mask);

#if 0
   len = vb->num_dest_pixels * sizeof(Pixel_Overlap_Info_Type);
   if (NULL == (vb->overlap_info = (Pixel_Overlap_Info_Type *)MALLOC (len)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        Var_free_value_buffer (vb);
        return NULL;
     }
#else
   vb->overlap_info = NULL;
#endif

   return vb;
}

static int maybe_realloc_value_buf (int ncid, Var_Value_Buffer_Type *vb,
                                    const char *var_name)
{
   TIO_Var_Info_Type vi;
   double *tmp;
   int i, need_num, len, num_src_values;

   if (-1 == TIO_inq_var (ncid, var_name, &vi))
     return -1;

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

   need_num = ((vb->num_src_pixels + vb->num_dest_pixels)
               * vb->num_values_per_pixel);
   if (need_num < vb->num_alloc_values)
     return 0;

   len = need_num * sizeof(double);
   if (NULL == (tmp = (double *)REALLOC (vb->src_values, len)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }

   num_src_values = vb->num_src_pixels * vb->num_values_per_pixel;
   vb->src_values = tmp;
   vb->dest_values = tmp + num_src_values;

   return 0;
}

static int read_var_values (int ncid, int var_grp, Var_Value_Buffer_Type *vb,
                            const char *var_name)
{
   int start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   int i, num_steps, num_pixels, num_values;
   int *step = NULL;
   double *var = NULL;
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

   num_values = num_pixels * vb->num_values_per_pixel;
   if (NULL == (var = (double *) MALLOC (num_values * sizeof(double))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_return;
     }

   if (-1 == TIO_get_var_section (var_grp, var_name,
                                  start, count, TIO_DOUBLE, var))
     goto cleanup_and_return;

   for (i = 0; i < num_pixels; i++)
     {
        int pix_xtrack = i % vb->num_xtrack;
        int pix_step_index = i / vb->num_xtrack;
        int pix = pix_xtrack + step[pix_step_index] * vb->num_xtrack;
        int k, nvpp = vb->num_values_per_pixel;
        for (k = 0; k < nvpp; k++)
          {
             vb->src_values[pix*nvpp + k] = var[i*nvpp + k];
          }
     }

   status = 0;
cleanup_and_return:
   FREE(step);
   FREE(var);

   return status;
}

static void copy_from_strided (int num, int stride,
                               const double *strided_values,
                               double *packed_values)
{
   int k;
   if (strided_values == packed_values)
     return;
   for (k = 0; k < num; k++)
     {
        packed_values[k] = strided_values[k*stride];
     }
}

static void copy_to_strided (int num, int stride,
                             const double *packed_values,
                             double *strided_values)
{
   int k;
   if (strided_values == packed_values)
     return;
   for (k = 0; k < num; k++)
     {
        strided_values[k*stride] = packed_values[k];
     }
}

static int parse_var_path (const char *var_path,
                           char **pgrp_path, char **pvar_name)
{
   char *grp_path = NULL, *var_name = NULL;
   const char *p;

   if (NULL == (p = strrchr (var_path, '/')))
     {
        grp_path = NULL;
        var_name = strdup (var_path);
     }
   else
     {
        if (NULL == (grp_path = strndup (var_path, p-var_path)))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return -1;
          }
        var_name = strdup (p + 1);
     }

   if (var_name == NULL)
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        FREE(grp_path);
        grp_path = NULL;
        return -1;
     }

   *pgrp_path = grp_path;
   *pvar_name = var_name;

   return 0;
}

int Var_apply_regrid (const Pixel_Regrid_Type *r, Var_Value_Buffer_Type *vb,
                      const char *var_path, char **files, int num_files)
{
   double *src_values = NULL;
   double *dest_values;
   char *var_name = NULL;
   char *grp_path = NULL;
   int i, ncid, status = -1;

   if (-1 == parse_var_path (var_path, &grp_path, &var_name))
     return -1;

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
               {
                  (void) TIO_close (ncid);
                  goto free_and_return;
               }
          }

        if (-1 == read_var_values (ncid, grp, vb, var_name))
          {
             (void) TIO_close (ncid);
             goto free_and_return;
          }

        (void) TIO_close (ncid);
     }

   /* Regridded result will be returned in vb->dest_values.
    * For 2D data, this is straightforward but, for
    * higher-dimensional variables, we regrid one slice
    * (dest->nx * dest->ny values) at a time.
    */
   if (vb->num_values_per_pixel == 1)
     {
        src_values = vb->src_values;
        dest_values = vb->dest_values;
     }
   else
     {
        int len = vb->num_src_pixels + vb->num_dest_pixels;
        if (NULL == (src_values = (double *)MALLOC (len * sizeof(double))))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto free_and_return;
          }
        dest_values = src_values + vb->num_src_pixels;
     }

   for (i = 0; i < vb->num_values_per_pixel; i++)
     {
        int j;
        copy_from_strided (vb->num_src_pixels,
                           vb->num_values_per_pixel,
                           vb->src_values + i, src_values);
        /* Pixel_regrid assumes dest_values initialized to
         * caller-defined INDEF */
        for (j = 0; j < vb->num_dest_pixels; j++)
          {
             dest_values[j] = HUGE_VAL;
          }
        if (-1 == Pixel_regrid (r, src_values, vb->src_mask,
                                dest_values, vb->overlap_info))
          {
             break;
          }
        copy_to_strided (vb->num_dest_pixels,
                         vb->num_values_per_pixel,
                         dest_values, vb->dest_values + i);
     }

   if (i == vb->num_values_per_pixel)
     {
        status = 0;
     }

free_and_return:
   FREE(grp_path);
   FREE(var_name);

   if (vb->num_values_per_pixel > 1)
     {
        FREE(src_values);
     }

   return status;
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
   double *lon=NULL, *lat=NULL;
   double dlon, dlat;
   int dim_lon, id_lon;
   int dim_lat, id_lat;
   int grp, start, count;
   int i, status = -1;

   if ((NULL == (lon = (double *)MALLOC (dest->nx * sizeof(double))))
       || (NULL == (lat = (double *)MALLOC (dest->ny * sizeof(double)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_exit;
     }

   dlon = (dest->xmax - dest->xmin) / dest->nx;
   dlat = (dest->ymax - dest->ymin) / dest->ny;
   for (i = 0; i < dest->nx; i++)
     {
        lon[i] = dest->xmin + (i + 0.5) * dlon;
     }
   for (i = 0; i < dest->ny; i++)
     {
        lat[i] = dest->ymin + (i + 0.5) * dlat;
     }

   /* assume dimensions are global even when lon-lat variables
    * are in a group */

   if ((-1 == TIO_def_dim (ncid, TEMPO_VAR_LONGITUDE, dest->nx, &dim_lon))
       || (-1 == TIO_def_dim (ncid, TEMPO_VAR_LATITUDE, dest->ny, &dim_lat)))
     goto cleanup_and_exit;

   if (lonlat_grp)
     {
        if (-1 == TIO_def_grp (ncid, lonlat_grp, &grp))
          goto cleanup_and_exit;
     }
   else grp = ncid;

   if ((-1 == TIO_def_var (grp, TEMPO_VAR_LONGITUDE, NC_FLOAT, 1, &dim_lon, &id_lon))
       || (-1 == TIO_put_text_attrs (grp, id_lon, lon_attrs)))
     goto cleanup_and_exit;

   if ((-1 == TIO_def_var (grp, TEMPO_VAR_LATITUDE, NC_FLOAT, 1, &dim_lat, &id_lat))
       || (-1 == TIO_put_text_attrs (grp, id_lat, lat_attrs)))
     goto cleanup_and_exit;

   start = 0;
   count = dest->nx;
   if (-1 == TIO_put_var_section (grp, TEMPO_VAR_LONGITUDE,
                                  &start, &count, TIO_DOUBLE, lon))
     goto cleanup_and_exit;

   start = 0;
   count = dest->ny;
   if (-1 == TIO_put_var_section (grp, TEMPO_VAR_LATITUDE,
                                  &start, &count, TIO_DOUBLE, lat))
     goto cleanup_and_exit;

   status = 0;
cleanup_and_exit:

   FREE(lon);
   FREE(lat);

   return status;
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
   char *in_var_name=NULL, *out_var_name=NULL;
   char *in_grp_path = NULL, *out_grp_path = NULL;
   int start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   size_t lon_dimlen, lat_dimlen;
   int lon_dimid, lat_dimid;
   int i, dims[TIO_MAX_VAR_DIMS];
   int in_grp, in_varid, out_grp, out_varid, num_dest_values;
   float fill_float = -NC_FILL_FLOAT;
   int shuffle=1, deflate=1, deflate_level=1;
   int status = -1;

   if ((-1 == parse_var_path (in_var_path, &in_grp_path, &in_var_name))
       || (-1 == parse_var_path (out_var_path, &out_grp_path, &out_var_name)))
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

   if ((-1 == TIO_def_var (out_grp, out_var_name, NC_FLOAT, vb->num_dims, dims, &out_varid))
       || (-1 == TIO_def_var_fill (out_grp, out_varid, 0, &fill_float))
       || (-1 == TIO_def_var_deflate (out_grp, out_varid, shuffle, deflate, deflate_level))
       || (-1 == TIO_copy_attrs (in_grp, in_varid, dontcopy_attr, out_grp, out_varid)))
     {
        goto free_and_return;
     }

   num_dest_values = vb->num_dest_pixels * vb->num_values_per_pixel;

   for (i = 0; i < num_dest_values; i++)
     {
        if (0 == isfinite(vb->dest_values[i]))
          vb->dest_values[i] = fill_float;
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
   if (-1 == TIO_put_var_section (out_grp, out_var_name, start, count,
                                  TIO_DOUBLE, vb->dest_values))
     goto free_and_return;

   status = 0;
free_and_return:
   FREE(in_var_name);
   FREE(in_grp_path);
   FREE(out_var_name);
   FREE(out_grp_path);
   return status;
}
