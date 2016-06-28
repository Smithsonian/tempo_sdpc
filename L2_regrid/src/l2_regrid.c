#include <float.h>
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
   int num_step;
   int num_xtrack;
   int num_src;
   int num_dest;
}
Regrid_Buffer_Type;

static void free_regrid_buffer (Regrid_Buffer_Type *rb)
{
   if (rb == NULL)
     return;
   FREE(rb->src_values);
   FREE(rb->src_mask);
   FREE(rb);
}

static Regrid_Buffer_Type *
new_regrid_buffer (const int *src_dims, const Pixel_Grid_Param_Type *dest)
{
   Regrid_Buffer_Type *rb = NULL;
   int len, len_mask;

   if (NULL == (rb = (Regrid_Buffer_Type *)MALLOC (sizeof *rb)))
     return NULL;

   rb->src_values = NULL;
   rb->src_mask = NULL;

   rb->num_step = src_dims[0];
   rb->num_xtrack = src_dims[1];

   rb->num_src = rb->num_step * rb->num_xtrack;
   rb->num_dest = dest->nx * dest->ny;

   len = rb->num_src * rb->num_dest * sizeof(double);
   len_mask = rb->num_src * sizeof(int);

   if ((NULL == (rb->src_values = (double *)MALLOC (len)))
       || (NULL == (rb->src_mask = (int *) MALLOC (len_mask))))
     {
        free_regrid_buffer (rb);
        return NULL;
     }
   rb->dest_values = rb->src_values + rb->num_src;

   memset ((char *)rb->src_mask, 0, len_mask);

   return rb;
}

static int read_var_values (const char *file, void *cl)
{
   Regrid_Buffer_Type *rb = (Regrid_Buffer_Type *)cl;
   TIO_Var_Info_Type vi;
   int i, start[3], count[3], ncid, num_steps, num_pixels;
   int *step = NULL;
   double *var = NULL;
   int status = -1;

   if (NC_NOERR != (status = nc_open (file, NC_NOWRITE, &ncid)))
     return -1;

   if (-1 == TIO_inq_var (ncid, rb->var_name, &vi))
     goto cleanup_and_return;

   num_steps = vi.dimlens[0];

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

   num_pixels = num_steps * rb->num_xtrack;

   if (NULL == (var = (double *) MALLOC (num_pixels * sizeof(double))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_steps;
   count[1] = rb->num_xtrack;

   if (-1 == TIO_get_var_section (ncid, rb->var_name,
                                  start, count, TIO_DOUBLE, var))
     goto cleanup_and_return;

   for (i = 0; i < num_pixels; i++)
     {
        int pix_xtrack = i % rb->num_xtrack;
        int pix_step_index = i / rb->num_xtrack;
        int pix = pix_xtrack + step[pix_step_index] * rb->num_xtrack;
        rb->src_values[pix] = var[i];
     }

   status = 0;
cleanup_and_return:
   FREE(step);
   FREE(var);
   nc_close (ncid);

   return status;
}

static int
regrid_variable (Regrid_Buffer_Type *rb, const char *var_name,
                 const Pixel_Regrid_Type *r, const char **files, int num_files)
{
   double *rb_dest_values = rb->dest_values;
   int i, rb_num_dest = rb->num_dest;

   rb->var_name = var_name;

   if (-1 == map_strings (files, num_files, read_var_values, rb))
     return -1;

   /* Pixel_regrid assumes dest_values is initialized to <invalid> */
   for (i = 0; i < rb_num_dest; i++)
     {
        rb_dest_values[i] = HUGE_VAL;
     }

   /* result returned in rb->dest_values */
   return Pixel_regrid (r, rb->src_values,
                        rb->src_mask, rb->dest_values, NULL);
}

static int read_dest_grid_params (FILE *fp, Pixel_Grid_Param_Type *dest)
{
   double pixel_size_deg = 0.05;

   /* FIXME:  should read grid from input parameter file */
   (void) fp;

   /* longitude [deg] */
   dest->xmin = -150.0;
   dest->nx = 2000;
   dest->xmax = dest->xmin + dest->nx * pixel_size_deg;
   /* latitude [deg] */
   dest->ymin = 17.0;
   dest->ny = 860;
   dest->ymax = dest->ymin + dest->ny * pixel_size_deg;

   return 0;
}

typedef struct Product_Type Product_Type;
struct Product_Type
{
   int type;
   const char *outfile;
   int num_var_names;
   const char **var_names;
   int num_input_files;
   const char **input_files;
};

enum {
   PRODUCT_TYPE_TEST=1
};

/* FIXME:  These product param definitions are temporary,
 * and are only intended to facilitate early development and
 * testing. (I hope this was obvious.)
 */

static const char *__Test_Product_Var_Names[] = {
   "column"
};
static const char *__Test_Product_File_Names[] = {
   "/tmp/test_l2l3_g00_grid.nc",
   "/tmp/test_l2l3_g01_grid.nc",
   "/tmp/test_l2l3_g02_grid.nc",
   "/tmp/test_l2l3_g03_grid.nc",
   "/tmp/test_l2l3_g04_grid.nc",
   "/tmp/test_l2l3_g05_grid.nc",
   "/tmp/test_l2l3_g06_grid.nc",
   "/tmp/test_l2l3_g07_grid.nc",
   "/tmp/test_l2l3_g08_grid.nc",
   "/tmp/test_l2l3_g09_grid.nc"
};

static struct Product_Type Product_List[] =
{
   {PRODUCT_TYPE_TEST,
        "test_l3_out.nc",
        1, __Test_Product_Var_Names,
       10, __Test_Product_File_Names},
   {0,NULL,0,NULL,0,NULL}
};

#define NC_CHECK_STATUS(s) \
   do {if (NC_NOERR != (s)) goto cleanup_and_exit; } while (0);

static int write_lonlat_grid (int ncid, const Pixel_Grid_Param_Type *dest)
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

   /* FIXME!! should be using libtio!!  */
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

static int write_variable (int ncid, const char *var_name,
                           const Regrid_Buffer_Type *regrid_buf)
{
   TIO_Var_Info_Type vi;
   const char coord_lonlat[] = "longitude latitude";
   size_t start[2], count[2];
   int lon_dimlen, lon_dimid;
   int lat_dimlen, lat_dimid;
   int status, id_var, i, dims[2];
   double *values = regrid_buf->dest_values;
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
   status = nc_def_var (ncid, var_name, NC_FLOAT, 2, dims, &id_var);
   NC_CHECK_STATUS(status);
   status = nc_put_att_text (ncid, id_var, "coordinates", strlen(coord_lonlat), coord_lonlat);
   NC_CHECK_STATUS(status);
   status = nc_def_var_fill(ncid, id_var, 0, &fill_float);
   NC_CHECK_STATUS(status);
   status = nc_def_var_deflate (ncid, id_var, shuffle, deflate, deflate_level);
   NC_CHECK_STATUS(status);

   for (i = 0; i < regrid_buf->num_dest; i++)
     {
        if (0 == isfinite(values[i]))
          values[i] = fill_float;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = lat_dimlen;
   count[1] = lon_dimlen;
   status = nc_put_vara_double (ncid, id_var, start, count, regrid_buf->dest_values);
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

static int make_l3_product (const Product_Type *prod, const Pixel_Grid_Param_Type *dest,
                            const Pixel_Regrid_Type *r, Regrid_Buffer_Type *regrid_buf)
{
   int ncid, i, close_status, status = -1;

   if (NC_NOERR != (status = nc_create (prod->outfile, NC_NETCDF4, &ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: creating %s (%s)",
                     __func__, prod->outfile, nc_strerror(status));
        return -1;
     }

   if (-1 == write_lonlat_grid (ncid, dest))
     goto return_status;

   for (i = 0; i < prod->num_var_names; i++)
     {
        if (-1 == regrid_variable (regrid_buf, prod->var_names[i], r,
                                   prod->input_files, prod->num_input_files))
          {
             goto return_status;
          }
        if (-1 == write_variable (ncid, prod->var_names[i], regrid_buf))
          goto return_status;
     }

   status = 0;
return_status:
   if (NC_NOERR != (close_status = nc_close(ncid)))
     {
        Tell_verror (TELL_IO_ERROR, "%s: closing %s (%s)",
                     __func__, prod->outfile, nc_strerror(close_status));
        return -1;
     }

   return status;
}

int main (int argc, const char **argv)
{
   Pixel_Grid_Param_Type dest;
   Pixel_Regrid_Type *r = NULL;
   Regrid_Buffer_Type *regrid_buf = NULL;
   Product_Type *prod;
   int src_dims[2], status = 1;
   const char **grid_files;
   int num_grid_files;

   num_grid_files = argc-1;
   grid_files     = argv+1;

   if (-1 == read_dest_grid_params (NULL, &dest))
     goto return_status;

   if (NULL == (r = Regrid_open (grid_files, num_grid_files, src_dims, &dest)))
     goto return_status;

   /* regrid_buf gets re-used in looping over variables/files */
   if (NULL == (regrid_buf = new_regrid_buffer (src_dims, &dest)))
     goto return_status;

   /* FIXME:  The same area weights should work for all L2 products produced
    * in the same scan so, in principle, we could loop over all the L2 products
    * and not just the set we used to derive the weights.
    * Each L2 product implies a set of variables to be regridded.
    */
   for (prod = Product_List; prod->var_names != NULL; prod++)
     {
        if (-1 == make_l3_product (prod, &dest, r, regrid_buf))
          goto return_status;
     }

   status = 0;
return_status:
   free_regrid_buffer (regrid_buf);
   Regrid_close (r);

   return status;
}
