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
   const char *var_name;
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
}
Var_Value_Buffer_Type;

void Var_free_value_buffer (Var_Value_Buffer_Type *vb)
{
   if (vb == NULL)
     return;
   FREE(vb->src_values);
   FREE(vb->src_mask);
   FREE(vb);
}

Var_Value_Buffer_Type *
Var_new_value_buffer (int dest_nx, int dest_ny, int src_num_step, int src_num_xtrack)
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

   return vb;
}

static int maybe_realloc_value_buf (int ncid, Var_Value_Buffer_Type *vb)
{
   TIO_Var_Info_Type vi;
   double *tmp;
   int i, need_num, len, num_src_values;

   if (-1 == TIO_inq_var (ncid, vb->var_name, &vi))
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

static int read_var_values (int ncid, Var_Value_Buffer_Type *vb)
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

   if (-1 == TIO_get_var_section (ncid, vb->var_name,
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

int Var_apply_regrid (const Pixel_Regrid_Type *r, Var_Value_Buffer_Type *vb,
                      const char *var_name, const char **files, int num_files)
{
   double *src_values, *dest_values;
   int i, ncid, status = -1;

   vb->var_name = var_name;

   for (i = 0; i < num_files; i++)
     {
        int open_status;
        if (NC_NOERR != (open_status = nc_open (files[i], NC_NOWRITE, &ncid)))
          {
             Tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s (%s)",
                          __func__, files[i], nc_strerror(open_status));
             return -1;
          }

        if (i == 0)
          {
             if (-1 == maybe_realloc_value_buf (ncid, vb))
               return -1;
          }

        if (-1 == read_var_values (ncid, vb))
          return -1;

        (void) nc_close (ncid);
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
             return -1;
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
                                dest_values, NULL))
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

   if (vb->num_values_per_pixel > 1)
     {
        FREE(src_values);
     }

   return status;
}

#define NC_CHECK_STATUS(s) \
   do {if (NC_NOERR != (s)) goto cleanup_and_exit; } while (0);

#define RETURN_STATUS_ONERR(s) \
   do {if (NC_NOERR != (s)) { \
         fprintf (stderr, "*** ERROR: %s\n", nc_strerror(s)); \
         return s; \
       }} while (0);

int Var_write_lonlat_grid (int ncid, const Pixel_Grid_Param_Type *dest)
{
   const char units_lon[] = "degrees_east";
   const char units_lat[] = "degrees_north";
   double *lon=NULL, *lat=NULL;
   double dlon, dlat;
   int dim_lon, id_lon;
   int dim_lat, id_lat;
   size_t start, count;
   int i, status;

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

   status = nc_def_dim (ncid, TEMPO_VAR_LONGITUDE, dest->nx, &dim_lon);
   NC_CHECK_STATUS(status);
   status = nc_def_dim (ncid, TEMPO_VAR_LATITUDE, dest->ny, &dim_lat);
   NC_CHECK_STATUS(status);

   status = nc_def_var (ncid, TEMPO_VAR_LONGITUDE, NC_FLOAT, 1, &dim_lon, &id_lon);
   NC_CHECK_STATUS(status);
   status = nc_put_att_text (ncid, id_lon, "units", strlen(units_lon), units_lon);
   NC_CHECK_STATUS(status);

   status = nc_def_var (ncid, TEMPO_VAR_LATITUDE, NC_FLOAT, 1, &dim_lat, &id_lat);
   NC_CHECK_STATUS(status);
   status = nc_put_att_text (ncid, id_lat, "units", strlen(units_lat), units_lat);
   NC_CHECK_STATUS(status);

   start = 0;
   count = dest->nx;
   status = nc_put_vara_double (ncid, id_lon, &start, &count, lon);
   NC_CHECK_STATUS(status);

   start = 0;
   count = dest->ny;
   status = nc_put_vara_double (ncid, id_lat, &start, &count, lat);
   NC_CHECK_STATUS(status);

cleanup_and_exit:

   FREE(lon);
   FREE(lat);

   return 0;
}

static int dontcopy_attribute (const char *attname)
{
   const char *do_not_copy_list[] = {"bounds", NULL};
   const char **a;
   for (a = do_not_copy_list; *a != NULL; a++)
     {
        if (0 == strcmp (*a, attname))
          return 1;
     }

   return 0;
}

static int copy_var_atts (int ncid_infile, int id_var_infile,
                          int ncid, int id_var)
{
   char attname[TIO_MAX_NAME_LEN];
   int status, attnum, num_atts;

   status = nc_inq_varnatts (ncid_infile, id_var_infile, &num_atts);
   RETURN_STATUS_ONERR(status);
   for (attnum = 0; attnum < num_atts; attnum++)
     {
        status = nc_inq_attname (ncid_infile, id_var_infile, attnum, attname);
        RETURN_STATUS_ONERR(status);
        if (dontcopy_attribute (attname))
          continue;
        status = nc_copy_att (ncid_infile, id_var_infile, attname,
                              ncid, id_var);
        RETURN_STATUS_ONERR(status);
     }

   return 0;
}

static int check_extra_dims (int ncid_infile, const TIO_Var_Info_Type *vi,
                             int ncid, int *dims)
{
   char dimname[TIO_MAX_NAME_LEN];
   int i, status;
   for (i = 2; i < vi->ndims; i++)
     {
        status = nc_inq_dimname (ncid_infile, vi->dimids[i], dimname);
        RETURN_STATUS_ONERR(status);
        /* does this dimension already exist in the output file? */
        status = nc_inq_dimid (ncid, dimname, &dims[i]);
        if (status == NC_NOERR)
          continue;
        status = nc_def_dim (ncid, dimname, vi->dimlens[i], &dims[i]);
        RETURN_STATUS_ONERR(status);
     }

   return 0;
}

int Var_write_values (int ncid, const Var_Value_Buffer_Type *vb,
               int ncid_infile, const char *var_name)
{
   TIO_Var_Info_Type vi;
   size_t start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   int lon_dimlen, lon_dimid;
   int lat_dimlen, lat_dimid;
   int status, i, dims[TIO_MAX_VAR_DIMS];
   int id_var, num_dest_values;
   float fill_float = -NC_FILL_FLOAT;
   int shuffle=1, deflate=1, deflate_level=1;

   if (-1 == TIO_inq_var(ncid, TEMPO_VAR_LONGITUDE, &vi))
     return -1;
   lon_dimlen = vi.dimlens[0];
   lon_dimid = vi.dimids[0];
   if (-1 == TIO_inq_var(ncid, TEMPO_VAR_LATITUDE, &vi))
     return -1;
   lat_dimlen = vi.dimlens[0];
   lat_dimid = vi.dimids[0];

   dims[0] = lat_dimid;
   dims[1] = lon_dimid;

   if (-1 == TIO_inq_var (ncid_infile, var_name, &vi))
     return -1;

   if (vi.ndims > 2)
     {
        status = check_extra_dims (ncid_infile, &vi, ncid, dims);
        NC_CHECK_STATUS(status);
     }

   status = nc_def_var (ncid, var_name, NC_FLOAT, vi.ndims, dims, &id_var);
   NC_CHECK_STATUS(status);
   status = nc_def_var_fill(ncid, id_var, 0, &fill_float);
   NC_CHECK_STATUS(status);
   status = nc_def_var_deflate (ncid, id_var, shuffle, deflate, deflate_level);
   NC_CHECK_STATUS(status);
   status = copy_var_atts (ncid_infile, vi.varid, ncid, id_var);
   NC_CHECK_STATUS(status);

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
   status = nc_put_vara_double (ncid, id_var, start, count, vb->dest_values);
   NC_CHECK_STATUS(status);

cleanup_and_exit:
   if (status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: writing variable '%s' (%s)\n",
                     __func__, var_name, nc_strerror(status));
     }

   return status ? -1 : 0;
}
