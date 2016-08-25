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
#include <strings.h>

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
   int processing_version;

   char *in_lonlat_grp;
   char *out_lonlat_grp;

   int num_var_names;
   char **in_var_names;
   char **out_var_names;
   char **var_qa_labels;
   int *value_types;

   int num_input_files;
   char **input_files;
};

static void free_product_type (Product_Type *p)
{
   if (p == NULL)
     return;

   if (p->in_var_names)
     {
        int i, num_shared = 3 * p->num_var_names + p->num_input_files;
        for (i = 0; i < num_shared; i++)
          {
             FREE(p->in_var_names[i]);
          }
        FREE(p->in_var_names);
     }

   FREE(p->value_types);
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
   int i, num_strings;

   if (NULL == (p = (Product_Type *) MALLOC (sizeof *p)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)p, 0, sizeof (*p));

   p->next = NULL;
   p->num_var_names = num_var_names;
   p->num_input_files = num_input_files;

   num_strings = 3 * num_var_names + num_input_files;

   p->in_var_names = (char **) MALLOC (num_strings * sizeof (char *));
   if (NULL == p->in_var_names)
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_product_type (p);
        return NULL;
     }
   memset ((char *)p->in_var_names, 0, num_strings * sizeof(char *));

   if (NULL == (p->value_types = (int *) MALLOC (num_var_names * sizeof(int))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_product_type (p);
        return NULL;
     }
   for (i = 0; i < num_var_names; i++)
     {
        p->value_types[i] = VALUE_IS_DOUBLE;
     }

   p->out_var_names = p->in_var_names + num_var_names;
   p->var_qa_labels = p->out_var_names + num_var_names;
   p->input_files = p->var_qa_labels + num_var_names;

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

static int set_value_type (int bitfield_type, int *value_type)
{
   switch (bitfield_type)
     {
      case  8: *value_type = VALUE_IS_UINT64; break;
      case -8: *value_type = VALUE_IS_INT64; break;
      case  4: *value_type = VALUE_IS_UINT; break;
      case -4: *value_type = VALUE_IS_INT; break;
      case  2: *value_type = VALUE_IS_USHORT; break;
      case -2: *value_type = VALUE_IS_SHORT; break;
      case  1: *value_type = VALUE_IS_UBYTE; break;
      case -1: *value_type = VALUE_IS_BYTE; break;
      case  0: *value_type = VALUE_IS_DOUBLE; break;
      default:
        Tell_verror (TELL_APPLICATION_ERROR,
                     "%s: unsupported value bitfield_type=%d", __func__, bitfield_type);
        return -1;
     }

   return 0;
}

static char *malloc_strcpy (const char *s)
{
   char *cpy = NULL;
   int len = strlen(s) + 1;
   if (NULL == (cpy = (char *) MALLOC (len)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memcpy (cpy, s, len);
   return cpy;
}

static Product_Type *init_product_type (const config_setting_t *setting)
{
   Product_Type *prod = NULL;
   config_setting_t *s, *vars, *longlat_group, *input_files;
   const char *name, *outfile, *in_grp, *out_grp;
   int i, num_vars, num_input_files, processing_version;

   if (setting == NULL)
     return NULL;

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "name", &name))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR, "%s: accessing name", __func__);
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "output_file", &outfile))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR, "%s: accessing output_file", __func__);
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_int (setting, "processing_version", &processing_version))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR, "%s: accessing processing_version", __func__);
        return NULL;
     }

   input_files = config_setting_get_member (setting, "input_files");
   if (NULL == input_files)
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing input_files list", __func__);
        return NULL;
     }
   num_input_files = config_setting_length (input_files);

   if (num_input_files == 0)
     {
        Tell_verror (TELL_USAGE_ERROR, "product %s: No input files", name);
        return NULL;
     }

   vars = config_setting_get_member (setting, "vars");
   if (NULL == vars)
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing vars list", __func__);
        return NULL;
     }
   num_vars = config_setting_length (vars);

   if (num_vars == 0)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "product %s: No variables to regrid", name);
        return NULL;
     }

   if (NULL == (prod = new_product_type (num_vars, num_input_files)))
     return NULL;

   prod->processing_version = processing_version;

   for (i = 0; i < num_vars; i++)
     {
        const char *in_name, *out_name, *var_qa_label;
        int bitfield_status, bitfield_type;
        if ((NULL == (s = config_setting_get_elem (vars, i)))
            || (CONFIG_TRUE != config_setting_lookup_string (s, "in", &in_name)))
          {
             Tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: accessing vars element i=%d", __func__, i);
             free_product_type (prod);
             return NULL;
          }

        if (CONFIG_TRUE != config_setting_lookup_string (s, "out", &out_name))
          out_name = in_name;

        if (CONFIG_TRUE != config_setting_lookup_string (s, "qa", &var_qa_label))
          var_qa_label = NULL;

        bitfield_status = config_setting_lookup_int (s, "bitfield_type", &bitfield_type);
        if (bitfield_status != CONFIG_TRUE)
          prod->value_types[i] = VALUE_IS_DOUBLE;
        else if (-1 == set_value_type (bitfield_type, &prod->value_types[i]))
          {
             free_product_type(prod);
             return NULL;
          }

        if ((NULL == (prod->in_var_names[i] = malloc_strcpy (in_name)))
            || (NULL == (prod->out_var_names[i] = malloc_strcpy (out_name))))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_product_type (prod);
             return NULL;
          }

        if (var_qa_label == NULL)
          prod->var_qa_labels[i] = NULL;
        else if (NULL == (prod->var_qa_labels[i] = malloc_strcpy (var_qa_label)))
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
        if (NULL == (prod->input_files[i] = malloc_strcpy (file)))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_product_type (prod);
             return NULL;
          }
     }

   if ((NULL == (longlat_group = config_setting_get_member (setting, "longlat_group")))
       || (CONFIG_TRUE != config_setting_lookup_string (longlat_group, "in", &in_grp)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing longlat_group", __func__);
        free_product_type (prod);
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (longlat_group, "out", &out_grp))
     out_grp = in_grp;

   if ((NULL == (prod->name = malloc_strcpy (name)))
       || (NULL == (prod->outfile = malloc_strcpy (outfile)))
       || (NULL == (prod->in_lonlat_grp = malloc_strcpy (in_grp)))
       || (NULL == (prod->out_lonlat_grp = malloc_strcpy (out_grp))))
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

   if (NULL == (s = config_lookup (&cfg, "data_products")))
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
                            const Pixel_Regrid_Type *r,
                            Var_Value_Buffer_Type *vb,
                            TIO_Scan_Ident_Type *lst)
{
   int ncid=INT_MAX, ncid_infile=INT_MAX, i;
   int status = -1;

   if (-1 == TIO_create (prod->outfile, NC_NETCDF4, &ncid))
     return -1;

   if (-1 == Var_write_lonlat_grid (ncid, prod->out_lonlat_grp, dest))
     goto return_status;

   if ((lst != NULL)
       && (-1 == TIO_write_scan_ident (ncid, lst)))
     goto return_status;

   if (-1 == TIO_label_product (ncid, prod->name, prod->processing_version))
     goto return_status;

   /* The first input file establishes each variable's dimensionality */
   if (-1 == TIO_open (prod->input_files[0], NC_NOWRITE, &ncid_infile))
     goto return_status;

   for (i = 0; i < prod->num_var_names; i++)
     {
        int want_qa = (prod->var_qa_labels[i] != NULL);
        if (-1 == Var_apply_regrid (r, vb, prod->value_types[i],
                                    prod->in_var_names[i], want_qa,
                                    prod->input_files, prod->num_input_files))
          goto return_status;
        if (-1 == Var_write_values (ncid, vb, prod->out_var_names[i],
                                    prod->var_qa_labels[i],
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

static TIO_Scan_Ident_Type *
read_scan_ident (char **input_files, int num_input_files, const char *name)
{
   TIO_Scan_Ident_Type *lst = NULL;
   int i;

   if (NULL == (lst = TIO_new_scan_ident ()))
     return NULL;

   for (i = 0; i < num_input_files; i++)
     {
        int ncid, status;
        if (-1 == TIO_open (input_files[i], NC_NOWRITE, &ncid))
          goto free_and_return;
        if (name != NULL)
          {
             /* If we have a product_type name,
              * then every granule should match it */
             char buf[TIO_MAX_NAME_LEN];
             if (-1 == TIO_get_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, buf))
               goto free_and_return;
             if (0 != strcasecmp (buf, name))
               {
                  Tell_verror (TELL_APPLICATION_ERROR,
                               "%s: product_type mismatch: expected %s got %s",
                               __func__, name, buf);
                  goto free_and_return;
               }
          }
        status = TIO_attach_granule_ident (ncid, lst);
        (void) TIO_close (ncid);
        if (status < 0)
          goto free_and_return;
     }

   return lst;

free_and_return:
   TIO_free_scan_ident (lst);
   return NULL;
}

int main (int argc, char **argv)
{
   const char *param_file = DEFAULT_PARAM_FILE;
   Var_Value_Buffer_Type *vb = NULL;
   TIO_Scan_Ident_Type *lst = NULL;
   Pixel_Regrid_Type *r = NULL;
   Product_Type *product_list = NULL;
   Product_Type *prod = NULL;
   Pixel_Grid_Param_Type dest;
   int src_num_steps, src_num_xtrack;
   int expect_scan_ident = 1;
   int status = 1;

   Tell_open ("L2_regrid", -1, -1);

   if (0 == strcmp (argv[1], "--noident"))
     {
        expect_scan_ident = 0;
        argv++;
        argc--;
     }

   if (argc > 1)
     param_file = argv[1];

   if (-1 == parse_param_file (param_file, &dest, &product_list))
     return 1;

   /* Compute pixel overlaps using the first set of products */
   prod = product_list;

   r = Regrid_open (&dest, prod->input_files, prod->num_input_files,
                    prod->in_lonlat_grp);
   if (NULL == r)
     goto return_status;

   Pixel_regrid_get_srcdims (r, &src_num_steps, &src_num_xtrack);

   vb = Var_new_value_buffer (dest.nx, dest.ny,
                              src_num_steps, src_num_xtrack);
   if (NULL == vb)
     goto return_status;

   for (prod = product_list; prod != NULL; prod = prod->next)
     {
        if ((expect_scan_ident != 0)
            && (NULL == (lst = read_scan_ident (prod->input_files, prod->num_input_files, prod->name))))
          goto return_status;
        if (-1 == make_l3_product (prod, &dest, r, vb, lst))
          goto return_status;
        TIO_free_scan_ident (lst);
        lst = NULL;
     }

   status = 0;
return_status:
   free_product_list (product_list);
   Var_free_value_buffer (vb);
   Regrid_close (r);
   TIO_free_scan_ident (lst);

   return status;
}
