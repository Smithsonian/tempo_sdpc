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

#ifndef REALLOC
# define REALLOC realloc
#endif

#ifndef MALLOC
# define MALLOC malloc
#endif

#ifndef FREE
# define FREE free
#endif

static int read_dest_grid_params (FILE *fp, Pixel_Grid_Param_Type *dest)
{
   double pixel_size_deg = 0.05;

   /* FIXME:  should read grid from input parameter file */
   (void) fp;

   /* longitude [deg] */
   dest->xmin = -150.0;
   dest->nx = 2200;
   dest->xmax = dest->xmin + dest->nx * pixel_size_deg;
   /* latitude [deg] */
   dest->ymin = 17.0;
   dest->ny = 900;
   dest->ymax = dest->ymin + dest->ny * pixel_size_deg;

   return 0;
}

typedef struct Product_Type Product_Type;
struct Product_Type
{
   int type;
   const char *outfile;

   const char *in_lonlat_grp;
   const char *out_lonlat_grp;

   int num_var_names;
   const char **in_var_names;
   const char **out_var_names;

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

static const char *__Test_Product_In_Var_Names[] = {
   "var0", "var1", "/test_group/var2", "/test_group/var3"
};
static const char *__Test_Product_Out_Var_Names[] = {
   "var0", "var1", "/foo/bar/var2", "/foo/baz/var3"
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
        "/",    /* input file group containing lon-lat variables */
        "/geometry",   /* output file group containing lon-lat variables */
        4, __Test_Product_In_Var_Names, __Test_Product_Out_Var_Names,
       10, __Test_Product_File_Names},
   {0,NULL,NULL,NULL,0,NULL,NULL,0,NULL}
};

static int make_l3_product (const Product_Type *prod,
                            const Pixel_Grid_Param_Type *dest,
                            const Pixel_Regrid_Type *r, Var_Value_Buffer_Type *vb)
{
   int ncid=INT_MAX, ncid_infile=INT_MAX, i;
   int status = -1;

   if (-1 == TIO_create (prod->outfile, NC_NETCDF4, &ncid))
     return -1;

   if (-1 == Var_write_lonlat_grid (ncid, prod->out_lonlat_grp, dest))
     goto return_status;

   /* The first input file establishes each variable's dimensionality */
   if (-1 == TIO_open (prod->input_files[0], NC_NOWRITE, &ncid_infile))
     goto return_status;

   for (i = 0; i < prod->num_var_names; i++)
     {
        if (-1 == Var_apply_regrid (r, vb, prod->in_var_names[i],
                                    prod->input_files, prod->num_input_files))
          goto return_status;
        if (-1 == Var_write_values (ncid, vb, prod->out_var_names[i],
                                    ncid_infile, prod->in_var_names[i]))
          goto return_status;
     }

   status = 0;
return_status:
   if (ncid_infile != INT_MAX)
     {
        (void) TIO_close (ncid_infile);
     }
   if (-1 == TIO_close(ncid))
     return -1;

   return status;
}

int main (int argc, const char **argv)
{
   Pixel_Grid_Param_Type dest;
   Pixel_Regrid_Type *r = NULL;
   Var_Value_Buffer_Type *vb = NULL;
   Product_Type *prod;
   const char **grid_files;
   int num_grid_files, src_num_step, src_num_xtrack;
   int status = 1;

   num_grid_files = argc-1;
   grid_files     = argv+1;

   if (-1 == read_dest_grid_params (NULL, &dest))
     goto return_status;

   if (NULL == (r = Regrid_open (&dest, grid_files, num_grid_files,
                                 Product_List->in_lonlat_grp,
                                 &src_num_step, &src_num_xtrack)))
     goto return_status;

   if (NULL == (vb = Var_new_value_buffer (dest.nx, dest.ny,
                                           src_num_step, src_num_xtrack)))
     goto return_status;

   for (prod = Product_List; prod->in_var_names != NULL; prod++)
     {
        if (-1 == make_l3_product (prod, &dest, r, vb))
          goto return_status;
     }

   status = 0;
return_status:
   Var_free_value_buffer (vb);
   Regrid_close (r);

   return status;
}
