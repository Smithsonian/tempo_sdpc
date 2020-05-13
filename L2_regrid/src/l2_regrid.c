/** @file l2_regrid.c
 *  @brief Main program; parameter file parsing
 */

#include "defs.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <wordexp.h>

#include <tell.h>
#include <netcdf.h>
#include <tio.h>
#include <tio_template.h>

#include <libconfig.h>

#include "poly.h"
#include "pixel.h"
#include "regrid.h"
#include "var.h"

typedef struct Product_Type Product_Type;
struct Product_Type
{
   Product_Type *next;

   char *name;
   char *outfile;
   char *shortname;
   char *metadata_template;
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

   int __num_files;
   char **__filenames;
};

static void usage (void)
{
   fprintf (stderr, "Usage: L2_regrid [options] [config-file]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help            Print this usage message\n");
   fprintf (stderr, "   -c | --config FILE     Use this configuration file\n");
   fprintf (stderr, "   -i | --ignore          Ignore scan metadata fields\n");
   fprintf (stderr, "   -d | --diagnostic      Generate diagnostic output\n");
   exit (EXIT_SUCCESS);
}

static void free_product_type (Product_Type *p)
{
   int i;

   if (p == NULL)
     return;

   if (p->in_var_names)
     {
        int num_shared = 3 * p->num_var_names;
        for (i = 0; i < num_shared; i++)
          {
             FREE(p->in_var_names[i]);
          }
        FREE(p->in_var_names);
     }

   if (p->__filenames)
     {
        for (i = 0; i < p->__num_files; i++)
          {
             FREE(p->__filenames[i]);
          }
        FREE(p->__filenames);
     }

   FREE(p->value_types);
   FREE(p->name);
   FREE(p->shortname);
   FREE(p->metadata_template);
   FREE(p->in_lonlat_grp);
   FREE(p->out_lonlat_grp);
   FREE(p);
}

static Product_Type *new_product_type (int num_var_names)
{
   Product_Type *p = NULL;
   int i, num_strings;

   if (NULL == (p = (Product_Type *) MALLOC (sizeof *p)))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)p, 0, sizeof (*p));

   p->num_var_names = num_var_names;

   num_strings = 3 * num_var_names;

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
                             int *num, double *min, double *max,
                             int *num_pixel_sub)
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

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_pixel_sub", num_pixel_sub))
     *num_pixel_sub = 0;

   *max = *min + (*num) * delta;
   return 0;
}

static int init_dest_grid (const config_setting_t *setting,
                           Pixel_Grid_Param_Type *dest)
{
   config_setting_t *s;

   if ((NULL == (s = config_setting_get_member (setting, "longitude")))
       || (-1 == lookup_grid_spec (s, &dest->nx, &dest->xmin, &dest->xmax, &dest->num_extra_xpoints)))
     return -1;

   if ((NULL == (s = config_setting_get_member (setting, "latitude")))
       || (-1 == lookup_grid_spec (s, &dest->ny, &dest->ymin, &dest->ymax, &dest->num_extra_ypoints)))
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

static int read_filename_list (const char *list_file, int *num_filesp,
                               char ***filenamesp)
{
   FILE *fp;
   char **filenames = NULL;
   int num_files = 0;
   int num_allocated = 0;
   int return_status = -1;

   *num_filesp = 0;
   *filenamesp = NULL;

   if (list_file == NULL)
     return -1;

/* Using a low initial value to exercise the realloc. */
#define DEFAULT_NUM_FILES 2

   if (NULL == (filenames = (char **)MALLOC (DEFAULT_NUM_FILES * sizeof (char *))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   num_allocated = DEFAULT_NUM_FILES;

   if (NULL == (fp = fopen (list_file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s", __func__, list_file);
        goto return_error;
     }

   while (1)
     {
        char *newline;
        char *line;
        int status;

        status = tio_fgets (&line, NULL, fp);
        if (status == -1)
          goto return_error;
        if (status == 0)
          break;

        if ((*line == '#') || (*line == '\n') || (*line == 0))
          {
             FREE(line);
             line = NULL;
             continue;
          }

        newline = strchr(line, '\n');
        if (newline) *newline = 0;

        filenames[num_files++] = line;

        if (num_files == num_allocated)
          {
             char **more_files = NULL;
             int new_num = 2*num_allocated;
             if (NULL == (more_files = REALLOC(filenames, new_num * sizeof(char *))))
               {
                  tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
                  goto return_error;
               }
             filenames = more_files;
             num_allocated = new_num;
          }
     }

   return_status = 0;
return_error:
   (void) fclose (fp);
   if (return_status)
     {
        FREE(filenames);
        filenames = NULL;
        num_files = 0;
     }

   *filenamesp = filenames;
   *num_filesp = num_files;

   return return_status;
}

static char *expand_path (const char *path)
{
   wordexp_t we;
   char *path_exp = NULL;

   memset ((char *)&we, 0, sizeof(wordexp_t));

   if ((0 != wordexp (path, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, path);
        goto return_error;
     }

   if (NULL == (path_exp = strdup (we.we_wordv[0])))
     {
        tell_verror (TELL_MALLOC_ERROR,
                     "%s: strdup failed", __func__);
     }

return_error:
   wordfree (&we);
   return path_exp;
}

static int init_product_type_metadata (const config_setting_t *s,
                                       Product_Type *prod)
{
   const char *shortname;
   const char *metadata_template;

   if (CONFIG_TRUE != config_setting_lookup_string (s, "shortname", &shortname))
     return -1;
   if (CONFIG_TRUE != config_setting_lookup_string (s, "template_file", &metadata_template))
     return -1;

   if (NULL == (prod->shortname = malloc_strcpy (shortname)))
     return -1;
   if (NULL == (prod->metadata_template = expand_path (metadata_template)))
     return -1;

   return 0;
}

static int init_product_type (const config_setting_t *setting,
                              Product_Type **prodp)
{
   Product_Type *prod = NULL;
   config_setting_t *s, *vars, *longlat_group, *s_meta;
   const char *name, *in_grp, *out_grp, *list_file;
   int i, num_vars;

   *prodp = NULL;

   if (setting == NULL)
     return -1;

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "filename_list", &list_file))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR, "%s: accessing filename_list", __func__);
        return -1;
     }

   /* silently ignore a missing list of filenames */
   if (0 != access (list_file, F_OK))
     return 0;

   if (CONFIG_TRUE != config_setting_lookup_string (setting, "name", &name))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR, "%s: accessing name", __func__);
        return -1;
     }

   vars = config_setting_get_member (setting, "vars");
   if (NULL == vars)
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing vars list", __func__);
        return -1;
     }
   num_vars = config_setting_length (vars);

   if (num_vars == 0)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "product %s: No variables to regrid", name);
        return -1;
     }

   if (NULL == (prod = new_product_type (num_vars)))
     return -1;

   if (0 != read_filename_list (list_file, &prod->__num_files, &prod->__filenames))
     {
        free_product_type (prod);
        return -1;
     }

   if (prod->__num_files < 2)
     {
        free_product_type (prod);
        /* silently ignore a missing list of input files */
        return 0;
     }

   prod->outfile = prod->__filenames[0];
   prod->input_files = prod->__filenames + 1;
   prod->num_input_files = prod->__num_files - 1;

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
             return -1;
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
             return -1;
          }

        if ((NULL == (prod->in_var_names[i] = malloc_strcpy (in_name)))
            || (NULL == (prod->out_var_names[i] = malloc_strcpy (out_name))))
          {
             free_product_type (prod);
             return -1;
          }

        if (var_qa_label == NULL)
          prod->var_qa_labels[i] = NULL;
        else if (NULL == (prod->var_qa_labels[i] = malloc_strcpy (var_qa_label)))
          {
             free_product_type (prod);
             return -1;
          }
     }

   if ((NULL == (longlat_group = config_setting_get_member (setting, "longlat_group")))
       || (CONFIG_TRUE != config_setting_lookup_string (longlat_group, "in", &in_grp)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing longlat_group", __func__);
        free_product_type (prod);
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (longlat_group, "out", &out_grp))
     out_grp = in_grp;

   /* may be absent */
   if (NULL != (s_meta = config_setting_get_member (setting, "metadata")))
     {
        if (0 != init_product_type_metadata (s_meta, prod))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: reading metadata settings for %s product",
                          __func__, name);
             free_product_type (prod);
             return -1;
          }
     }

   if ((NULL == (prod->name = malloc_strcpy (name)))
       || (NULL == (prod->in_lonlat_grp = malloc_strcpy (in_grp)))
       || (NULL == (prod->out_lonlat_grp = malloc_strcpy (out_grp))))
     {
        free_product_type (prod);
        return -1;
     }

   *prodp = prod;
   return 0;
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
        config_setting_t *prod_cfg = config_setting_get_elem (s, i);
        Product_Type *prod;

        if (0 != init_product_type (prod_cfg, &prod))
          goto cleanup_and_return;

        if (prod)
          {
             prod->next = plist;
             plist = prod;
          }
     }

   if (plist)
     {
        status = 0;
        *product_list = plist;
     }
   else
     {
        tell_verror (TELL_UNKNOWN_ERROR, "%s: no product files to regrid?",
                     __func__);
     }

cleanup_and_return:
   if (status)
     {
        free_product_list (plist);
        plist = NULL;
     }

   config_destroy(&cfg);
   return status;
}

static int meta_set_bounding_polygon (TIO_Meta_Type *meta,
                                      const Pixel_Grid_Param_Type *dest)
{
   float lon[4], lat[4];
   float centroid_lon, centroid_lat;

   /* The following standard keyword values are set:
    * @verbatim
    *   polygon_longitudes      boundary polygon longitudes
    *   polygon_latitudes       boundary polygon latitudes
    *   polygon_sequence        integer indices giving the sequence in which the
    *                           (lon,lat) points trace the boundary in CCW order
    *   centroid_mean_longitude longitude of the polygon centroid
    *   centroid_mean_latitude  latitude of the polygon centroid
    * @endverbatim
    */

   lon[0] = dest->xmin;  lat[0] = dest->ymin;
   lon[1] = dest->xmax;  lat[1] = dest->ymin;
   lon[2] = dest->xmax;  lat[2] = dest->ymax;
   lon[3] = dest->xmin;  lat[3] = dest->ymax;

   if ((0 != tio_meta_set_acdd_geospatial_bounds (meta, lon, lat, 4))
       || (0 != tio_meta_set_odl_bounding_polygon (meta, lon, lat, 4)))
     return -1;

   centroid_lon = 0.5 * (dest->xmin + dest->xmax);
   centroid_lat = 0.5 * (dest->ymin + dest->ymax);

   if ((0 != tio_meta_set (meta, "centroid_mean_longitude",  TIO_META_TYPE_FLOAT, 1, &centroid_lon))
       || (0 != tio_meta_set (meta, "centroid_mean_latitude",  TIO_META_TYPE_FLOAT, 1, &centroid_lat)))
     {
        return -1;
     }

   return 0;
}

static int write_metadata (TIO_Meta_Type *meta, int ncid,
                           const Product_Type *prod,
                           const Pixel_Grid_Param_Type *dest,
                           const TIO_Scan_Ident_Type *lst)
{
   const char *version_string = "0.1.0"; /* FIXME */
   int i;

   if (0 != tio_meta_set_datetime_production (meta))
     return -1;

   if (lst)
     {
        if (0 != tio_meta_set_datetime_range_scan (meta, lst))
          return -1;
     }

   if (0 != tio_meta_set_standard (meta, prod->outfile, prod->shortname,
                                   prod->processing_version, version_string))
     return -1;

   for (i = 0; i < prod->num_input_files; i++)
     {
        const char *basename = prod->input_files[i];
        const char *p;
        if (NULL != (p = strrchr (basename, '/')))
          {
             basename = p + 1;
          }
        if (0 != tio_meta_append_string (meta, "input_pointer", basename))
          return -1;
     }

   if (0 != meta_set_bounding_polygon (meta, dest))
     return -1;

   if (0 != tio_meta_write_ncattr (meta, ncid))
     return -1;

   if (prod->metadata_template)
     {
        if (0 != tio_meta_expand_file (meta, prod->metadata_template, prod->outfile))
          return -1;
     }

   return 0;
}

static int make_l3_product (const Product_Type *prod,
                            const Pixel_Grid_Param_Type *dest,
                            const Pixel_Regrid_Type *r,
                            Var_Value_Buffer_Type *vb,
                            TIO_Scan_Ident_Type *lst)
{
   TIO_Meta_Type *meta = NULL;
   int ncid=INT_MAX, ncid_infile=INT_MAX, i;
   int status = -1;

   if (-1 == TIO_create (prod->outfile, NC_NETCDF4, &ncid))
     return -1;

   if (NULL == (meta = tio_meta_open ()))
     goto return_status;

   if (-1 == Var_write_lonlat_grid (ncid, prod->out_lonlat_grp, dest))
     goto return_status;

   if ((lst != NULL)
       && (-1 == TIO_write_scan_ident (ncid, lst)))
     goto return_status;

   if (-1 == TIO_label_product (ncid, prod->name, 3, prod->processing_version))
     goto return_status;

   if (0 != write_metadata (meta, ncid, prod, dest, lst))
     goto return_status;

   /* The first input file establishes each variable's dimensionality */
   if (-1 == TIO_open (prod->input_files[0], NC_NOWRITE, &ncid_infile))
     goto return_status;

   if (0 != tio_use_file_epoch (ncid_infile))
     goto return_status;
   if (0 != tio_write_epoch_timestamp (ncid, NC_GLOBAL))
     goto return_status;

   for (i = 0; i < prod->num_var_names; i++)
     {
        int want_qa = (prod->var_qa_labels[i] != NULL);
        int apply_regrid_status =
          Var_apply_regrid (r, vb, prod->value_types[i],
                            prod->in_var_names[i], want_qa,
                            prod->input_files, prod->num_input_files);
        if (apply_regrid_status != 0)
          {
             if (apply_regrid_status > 0)
               continue;
             goto return_status;
          }
        if (-1 == Var_write_values (ncid, vb, prod->out_var_names[i],
                                    prod->var_qa_labels[i],
                                    ncid_infile, prod->in_var_names[i]))
          {
             goto return_status;
          }
     }

   status = 0;
return_status:
   if (ncid_infile != INT_MAX)
     {
        (void) TIO_close (ncid_infile);
     }
   if (-1 == TIO_close(ncid))
     return -1;
   tio_meta_close (meta);

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

static int read_processing_version (const char *file, int *processing_version)
{
   int ncid, status;
   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;
   status = TIO_get_att (ncid, NC_GLOBAL, "processing_version", NC_INT, processing_version);
   TIO_close (ncid);
   return status;
}

int main (int argc, char **argv)
{
   const char appname[] = "L2_regrid";
   const char *param_file = "l2_regrid.cfg";
   Var_Value_Buffer_Type *vb = NULL;
   TIO_Scan_Ident_Type *lst = NULL;
   Pixel_Regrid_Type *r = NULL;
   Product_Type *product_list = NULL;
   Product_Type *prod = NULL;
   Pixel_Grid_Param_Type dest;
   int src_num_steps, src_num_xtrack;
   int expect_scan_ident = 1;
   int status = 1;
   int want_diagnostic_output = 0;
   static struct option long_options[] =
     {
        {"help",       no_argument,       0, 'h'},
        {"config",     required_argument, 0, 'c'},
        {"ignore",     no_argument,       0, 'i'},
        {"diagnostic", no_argument,       0, 'd'},
        {0,0,0,0}
     };

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "dhic:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "*** Error: getopt returned character %d??", c);
             goto return_status;
             break;
           case 'c':
             param_file = optarg;
             break;
           case 'h':
             usage();
             break;
           case 'i':
             expect_scan_ident = 0;
             break;
           case 'd':
             want_diagnostic_output++;
             break;
          }
     }

   tell_open (appname, -1, -1);

   if (optind < argc)
     {
        param_file = argv[optind++];
     }

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   if (-1 == parse_param_file (param_file, &dest, &product_list))
     goto return_status;

   if (NULL != getenv ("SDPC_REGRID_DIAGNOSTICS"))
     want_diagnostic_output++;

   /* Compute pixel overlaps using the first set of products */
   prod = product_list;

   Regrid_diagnostic_output (want_diagnostic_output);

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
        if (0 != read_processing_version (prod->input_files[0], &prod->processing_version))
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
   tell_close ();

   return status;
}
