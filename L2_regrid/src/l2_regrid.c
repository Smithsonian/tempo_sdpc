/** @file l2_regrid.c
 *  @brief Main program; parameter file parsing
 */

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

#include <libconfig.h>

#include "poly.h"
#include "pixel.h"
#include "regrid.h"
#include "var.h"

#define DEFAULT_PARAM_FILE  "l2_regrid.cfg"

typedef struct Product_Type Product_Type;
struct Product_Type
{
   Product_Type *next;

   char *name;
   char *outfile;

   char *in_lonlat_grp;
   char *out_lonlat_grp;

   int num_var_names;
   char **in_var_names;
   char **out_var_names;

   int num_input_files;
   char **input_files;
};

static void free_product_type (Product_Type *p)
{
   int i;

   if (p->in_var_names)
     {
        int num_shared = 2 * p->num_var_names + p->num_input_files;
        for (i = 0; i < num_shared; i++)
          {
             FREE(p->in_var_names[i]);
          }
        FREE(p->in_var_names);
     }

   FREE(p->name);
   FREE(p->outfile);
   FREE(p->in_lonlat_grp);
   FREE(p->out_lonlat_grp);
   FREE(p);
}

static Product_Type *new_product_type (int num_var_names,
                                       int num_input_files)
{
   Product_Type *p = NULL;
   int num_strings;

   if (NULL == (p = (Product_Type *) MALLOC (sizeof *p)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)p, 0, sizeof (*p));

   p->next = NULL;
   p->num_var_names = num_var_names;
   p->num_input_files = num_input_files;

   num_strings = 2 * num_var_names + num_input_files;

   p->in_var_names = (char **) MALLOC (num_strings * sizeof (char *));
   if (NULL == p->in_var_names)
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_product_type (p);
        return NULL;
     }
   memset ((char *)p->in_var_names, 0, num_strings * sizeof(char *));

   p->out_var_names = p->in_var_names + num_var_names;
   p->input_files = p->out_var_names + num_var_names;

   return p;
}

static void free_product_list (Product_Type *plist)
{
   while (plist != NULL)
     {
        Product_Type *p = plist->next;
        free_product_type (plist);
        plist = p;
     }
}

static int lookup_grid_spec (const config_setting_t *s,
                             int *num, double *min, double *max)
{
   double delta;
   if ((CONFIG_TRUE != config_setting_lookup_float (s, "min", min))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num", num))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "delta", &delta)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining target grid", __func__);
        return -1;
     }
   *max = *min + (*num) * delta;
   return 0;
}

static int init_dest_grid (const config_setting_t *setting,
                           Pixel_Grid_Param_Type *dest)
{
   config_setting_t *s;

   if ((NULL == (s = config_setting_get_member (setting, "longitude")))
       || (-1 == lookup_grid_spec (s, &dest->nx, &dest->xmin, &dest->xmax)))
     return -1;

   if ((NULL == (s = config_setting_get_member (setting, "latitude")))
       || (-1 == lookup_grid_spec (s, &dest->ny, &dest->ymin, &dest->ymax)))
     return -1;

   return 0;
}

static Product_Type *init_product_type (const config_setting_t *setting)
{
   Product_Type *prod = NULL;
   config_setting_t *s, *vars, *longlat_group, *input_files;
   const char *name, *outfile, *in_grp, *out_grp;
   int i, num_vars, num_input_files;

   if (setting == NULL)
     return NULL;

   input_files = config_setting_get_member (setting, "input_files");
   if (NULL == input_files)
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing input_files list", __func__);
        return NULL;
     }
   num_input_files = config_setting_length (input_files);

   vars = config_setting_get_member (setting, "vars");
   if (NULL == vars)
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing vars list", __func__);
        return NULL;
     }
   num_vars = config_setting_length (vars);

   if (NULL == (prod = new_product_type (num_vars, num_input_files)))
     return NULL;

   for (i = 0; i < num_vars; i++)
     {
        const char *in_name, *out_name;
        if ((NULL == (s = config_setting_get_elem (vars, i)))
            || (CONFIG_TRUE != config_setting_lookup_string (s, "in", &in_name))
            || (CONFIG_TRUE != config_setting_lookup_string (s, "out", &out_name)))
          {
             Tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: accessing vars element i=%d", __func__, i);
             free_product_type (prod);
             return NULL;
          }
        if ((NULL == (prod->in_var_names[i] = strdup (in_name)))
            || (NULL == (prod->out_var_names[i] = strdup (out_name))))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_product_type (prod);
             return NULL;
          }
     }

   for (i = 0; i < num_input_files; i++)
     {
        const char *file;
        if (NULL == (file = config_setting_get_string_elem (input_files, i)))
          {
             Tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: accessing input_files element i=%d", __func__, i);
             free_product_type (prod);
             return NULL;
          }
        if (NULL == (prod->input_files[i] = strdup (file)))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_product_type (prod);
             return NULL;
          }
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (setting, "name", &name))
       || (CONFIG_TRUE != config_setting_lookup_string (setting, "output_file", &outfile)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR, "%s", __func__);
        free_product_type (prod);
        return NULL;
     }

   if ((NULL == (longlat_group = config_setting_get_member (setting, "longlat_group")))
       || (CONFIG_TRUE != config_setting_lookup_string (longlat_group, "in", &in_grp))
       || (CONFIG_TRUE != config_setting_lookup_string (longlat_group, "out", &out_grp)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing longlat_group", __func__);
        free_product_type (prod);
        return NULL;
     }

   if ((NULL == (prod->name = strdup (name)))
       || (NULL == (prod->outfile = strdup (outfile)))
       || (NULL == (prod->in_lonlat_grp = strdup (in_grp)))
       || (NULL == (prod->out_lonlat_grp = strdup (out_grp))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_product_type (prod);
        return NULL;
     }

   return prod;
}

static int parse_param_file (const char *cfg_file,
                             Pixel_Grid_Param_Type *dest,
                             Product_Type **product_list)
{
   Product_Type *plist = NULL;
   config_t cfg;
   config_setting_t *s;
   int i, num_products, status = -1;

   *product_list = NULL;
   config_init (&cfg);

   if(0 == config_read_file (&cfg, cfg_file))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: Reading %s: %s:%d - %s",
                     __func__, cfg_file, config_error_file(&cfg),
                     config_error_line(&cfg), config_error_text(&cfg));
        goto cleanup_and_return;
     }

   if (NULL == (s = config_lookup (&cfg, "target_grid")))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing target_grid in param file: %s",
                     __func__, cfg_file);
        goto cleanup_and_return;
     }

   if (-1 == init_dest_grid (s, dest))
        goto cleanup_and_return;

   if (NULL == (s = config_lookup (&cfg, "supported_data_products")))
     goto cleanup_and_return;

   num_products = config_setting_length (s);

   for (i = 0; i < num_products; i++)
     {
        Product_Type *prod;
        config_setting_t *prod_cfg = config_setting_get_elem (s, i);

        if (NULL == (prod = init_product_type (prod_cfg)))
          goto cleanup_and_return;

        prod->next = plist;
        plist = prod;
     }

   status = 0;
   *product_list = plist;

cleanup_and_return:
   if (status)
     {
        free_product_list (plist);
        plist = NULL;
     }

   config_destroy(&cfg);
   return status;
}

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

int main (int argc, char **argv)
{
   const char *param_file = DEFAULT_PARAM_FILE;
   Var_Value_Buffer_Type *vb = NULL;
   Pixel_Regrid_Type *r = NULL;
   Product_Type *product_list = NULL;
   Product_Type *prod = NULL;
   Pixel_Grid_Param_Type dest;
   int src_num_step, src_num_xtrack;
   int status = 1;

   Tell_open ("L2_regrid", -1, -1);

   if (argc > 1)
     param_file = argv[1];

   if (-1 == parse_param_file (param_file, &dest, &product_list))
     return 1;

   /* Compute pixel overlaps using the first set of products */
   prod = product_list;

   r = Regrid_open (&dest, prod->input_files, prod->num_input_files,
                    prod->in_lonlat_grp, &src_num_step, &src_num_xtrack);
   if (NULL == r)
     goto return_status;

   vb = Var_new_value_buffer (dest.nx, dest.ny,
                              src_num_step, src_num_xtrack);
   if (NULL == vb)
     goto return_status;

   for (prod = product_list; prod != NULL; prod = prod->next)
     {
        if (-1 == make_l3_product (prod, &dest, r, vb))
          goto return_status;
     }

   status = 0;
return_status:
   free_product_list (product_list);
   Var_free_value_buffer (vb);
   Regrid_close (r);

   return status;
}
