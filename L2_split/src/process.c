#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libconfig.h>
#include <tell.h>
#include <netcdf.h>
#include <tio.h>
#include <tio_template.h>

#include "scan.h"
#include "filter.h"
#include "process.h"

typedef struct
{
   double *slant_column;
   double *amf_trop;
   double *amf_strat;
   int dim0;
   int dim1;
}
Product_Var_Type;

static int lookup_grid_spec (const config_setting_t *s,
                             int *num, double *min, double *max,
                             int *num_extra_points)
{
   double delta;

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "min", min))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "delta", &delta))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num", num)))
     {
        Tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: defining target grid", __func__);
        return -1;
     }
   *max = *min + *num * delta;

   if (CONFIG_TRUE != config_setting_lookup_int (s, "num_extra_points", num_extra_points))
     *num_extra_points = 0;

   return 0;
}

static int init_mesh (config_t *cfg, Pixel_Grid_Param_Type *mesh)
{
   config_setting_t *setting, *s;

   if (NULL == (setting = config_lookup (cfg, "mesh")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing mesh definition in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (s = config_setting_get_member (setting, "longitude")))
       || (-1 == lookup_grid_spec (s, &mesh->nx, &mesh->xmin, &mesh->xmax, &mesh->num_extra_xpoints)))
     return -1;

   if ((NULL == (s = config_setting_get_member (setting, "latitude")))
       || (-1 == lookup_grid_spec (s, &mesh->ny, &mesh->ymin, &mesh->ymax, &mesh->num_extra_ypoints)))
     return -1;

   return 0;
}

static double *
read_apriori_vertical_trop_column (const char *file,
                                   const Pixel_Grid_Param_Type *cm)
{
   TIO_Var_Info_Type vi;
   int ncid, num_pixels, start[2], count[2];
   const char var_name[] = "no2";
   double *vtrop = NULL;
   int i;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if (0 != TIO_inq_var (ncid, var_name, &vi))
     goto return_error;

   /* FIXME:  A more general implementation might read grid information,
    *         and then regrid as needed.
    */
   if ((cm->ny != (int) vi.dimlens[0])
       || (cm->nx != (int) vi.dimlens[1]))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: array size mismatch: input apriori vertical trop column (%s)",
                     __func__, file);
        goto return_error;
     }

   num_pixels = vi.dimlens[0] * vi.dimlens[1];

   if (NULL == (vtrop = (double *) MALLOC (num_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_error;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = vi.dimlens[0];
   count[1] = vi.dimlens[1];

   if (0 != TIO_get_var_section (ncid, var_name, start, count, TIO_DOUBLE, vtrop))
     goto return_error;

   (void) TIO_close(ncid);

   for (i = 0; i < num_pixels; i++)
     {
        vtrop[i] *= 1.e15;
     }

   return vtrop;

return_error:
   (void) TIO_close(ncid);
   FREE(vtrop);
   return NULL;
}

static double *
vertical_strat_column (const Pixel_Grid_Param_Type *mesh,
                       const Product_Var_Type *mesh_vars,
                       const double *apriori_vert_trop, double trop_thresh)
{
   double nan_value = nan("");
   const double *slant_column = mesh_vars->slant_column;
   const double *amf_trop = mesh_vars->amf_trop;
   const double *amf_strat = mesh_vars->amf_strat;
   double *vstrat = NULL;
   int i, num_pixels;

   num_pixels = mesh->nx * mesh->ny;

   if (NULL == (vstrat = (double *)MALLOC (num_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   for (i = 0; i < num_pixels; i++)
     {
        double slant_column_i = slant_column[i];
        double s_trop, v_strat;

        vstrat[i] = nan_value;

        if ((slant_column_i != slant_column_i)
            || (slant_column_i == DBL_MAX))
          continue;

        s_trop = apriori_vert_trop[i] * amf_trop[i];
        v_strat = (slant_column_i - s_trop) / amf_strat[i];
        if (s_trop < trop_thresh * amf_strat[i])
          vstrat[i] = v_strat;
     }

   return vstrat;
}

static double *
vertical_trop_column (const Pixel_Grid_Param_Type *mesh,
                      const Product_Var_Type *mesh_vars,
                      const double *vert_strat)
{
   double nan_value = nan("");
   const double *slant_column = mesh_vars->slant_column;
   const double *amf_trop = mesh_vars->amf_trop;
   const double *amf_strat = mesh_vars->amf_strat;
   double *vtrop = NULL;
   int i, num_pixels;

   num_pixels = mesh->nx * mesh->ny;

   if (NULL == (vtrop = (double *)MALLOC (num_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   for (i = 0; i < num_pixels; i++)
     {
        if (amf_trop[i] != 0)
          {
             vtrop[i] = ((slant_column[i]
                          - vert_strat[i] * amf_strat[i]) /amf_trop[i]);
          }
        else vtrop[i] = nan_value;
     }

   return vtrop;
}

static void free_params (Config_Type *params)
{
   if (params == NULL)
     return;
   FREE(params->testdata_file);
   FREE(params->apriori_trop_file);
}

static int init_params (config_t *cfg, Config_Type *params)
{
   config_setting_t *s;
   const char *testdata_file;
   const char *apriori_trop_file;

   if (CONFIG_TRUE != config_lookup_float (cfg, "trop_thresh", &params->trop_thresh))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading trop_thresh in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (s = config_lookup (cfg, "refdata")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing refdata in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "testdata_file", &testdata_file))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "apriori_trop_file", &apriori_trop_file)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading refdata filenames in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   FREE(params->testdata_file);
   FREE(params->apriori_trop_file);
   if ((NULL == (params->testdata_file = strdup (testdata_file)))
       || (NULL == (params->apriori_trop_file = strdup (apriori_trop_file))))
     {
        free_params (params);
        return -1;
     }

   return 0;
}

static void free_var_type (Product_Var_Type *pvt)
{
   if (pvt == NULL)
     return;
   FREE(pvt->slant_column);
   FREE(pvt->amf_trop);
   FREE(pvt->amf_strat);
   FREE(pvt);
}

static Product_Var_Type *alloc_var_type (int dim0, int dim1)
{
   Product_Var_Type *pvt = NULL;
   int num_pixels = dim0 * dim1;

   if (NULL == (pvt = (Product_Var_Type *)MALLOC (sizeof *pvt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)pvt, 0, sizeof *pvt);

   if ((NULL == (pvt->slant_column = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (pvt->amf_strat = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (pvt->amf_trop = (double *)MALLOC (num_pixels * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_var_type (pvt);
        return NULL;
     }

   pvt->dim0 = dim0;
   pvt->dim1 = dim1;

   return pvt;
}

static Product_Var_Type *
regrid_product_vars (Scan_Vars_Type *sv, const Pixel_Grid_Param_Type *mesh,
                     const Pixel_Regrid_Type *r_mesh)
{
   Product_Var_Type *mesh_grid = NULL;
   double fill_value = nan("");
   int *mesh_mask = NULL;   /* FIXME? */
   int status = -1;

   /* Map variables from Level 2 product scan grid to the uniform mesh grid
    * that will be used for computations.
    */
   if (NULL == (mesh_grid = alloc_var_type (mesh->nx, mesh->ny)))
     goto free_and_return;

   if ((0 != Pixel_regrid (r_mesh, mesh_mask, fill_value, sv->slant_column, mesh_grid->slant_column, NULL))
       || (0 != Pixel_regrid (r_mesh, mesh_mask, fill_value, sv->amf_trop, mesh_grid->amf_trop, NULL))
       || (0 != Pixel_regrid (r_mesh, mesh_mask, fill_value, sv->amf_strat, mesh_grid->amf_strat, NULL)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mapping from scan grid to mesh grid", __func__);
        goto free_and_return;
     }

   status = 0;
free_and_return:
   if (status)
     {
        free_var_type (mesh_grid);
        mesh_grid = NULL;
     }

   return mesh_grid;
}

static void replace_fill (double *a, int n, double old_fill_value,
                          double new_fill_value)
{
   int i;

   for (i = 0; i < n; i++)
     {
        double v = a[i];
        if ((v == old_fill_value) || (0 == isfinite(v)))
          {
             a[i] = new_fill_value;
          }
     }
}

static int write_split (const Scan_Type *st,
                        const Scan_Vars_Type *sv,
                        const Pixel_Regrid_Type *r_mesh,
                        const double *mesh_vert_trop,
                        const double *mesh_vert_strat)
{
   int *dqf = sv->data_quality_flag;
   int *mesh_mask = NULL; /* FIXME? */
   double *scan_values = NULL;
   double *vert_strat = NULL;
   double *vert_trop = NULL;
   double nan_value = nan("");
   double fill_value = nan_value;
   double fill_value_output = -1.0e30;  /* (should not be a NaN!) */
   int i, num_scan_pixels;
   int status = -1;

   /* Regrid strat/trop column values back to Level 2 product scan grid */
   num_scan_pixels = sv->num_steps * sv->num_xtrack;
   if (NULL == (scan_values = (double *) MALLOC (2 * num_scan_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   vert_trop = scan_values;
   vert_strat = scan_values + num_scan_pixels;

   if ((0 != Pixel_regrid_from_mesh (r_mesh, mesh_mask, fill_value, mesh_vert_strat, vert_strat))
       || (0 != Pixel_regrid_from_mesh (r_mesh, mesh_mask, fill_value, mesh_vert_trop, vert_trop)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mapping from mesh grid to scan grid", __func__);
        goto free_and_return;
     }

   /* Mask strat/trop variables using original data quality flag */
   for (i = 0; i < num_scan_pixels; i++)
     {
        if (dqf[i] != 0)
          {
             vert_strat[i] = nan_value;
             vert_trop[i] = nan_value;
          }
     }

   /* For output, replace fill_value, NaN, Inf with standard fill value */
   if (isnan(fill_value) || (fill_value_output != fill_value))
     {
        replace_fill (vert_trop, num_scan_pixels, fill_value, fill_value_output);
        replace_fill (vert_strat, num_scan_pixels, fill_value, fill_value_output);
     }

   /* Write strat/trop vertical columns on Level 2 product scan grid */
   if (0 != scan_write_split (st, fill_value_output, vert_trop, vert_strat))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: writing strat/trop vertical column values",
                     __func__);
        goto free_and_return;
     }

   status = 0;
free_and_return:
   FREE(scan_values);
   return status;
}

static int filter_vert_strat (const Pixel_Grid_Param_Type *mesh,
                              double *mesh_vert_strat, config_t *cfg)
{
   Filter_Type *flt = NULL;
   int status;

   /* read filter parameters */
   if (NULL == (flt = filter_open (cfg)))
     return -1;

   status = flt->filter_apply (flt, mesh, mesh_vert_strat);

   flt->filter_delete (flt);

   return status;
}

int process_files (config_t *cfg, int num_files, char **files)
{
   Config_Type params = {0};
   Pixel_Grid_Param_Type mesh = {0};
   Scan_Type *st = NULL;
   Scan_Vars_Type *sv = NULL;
   Pixel_Regrid_Type *r_mesh = NULL;
   Product_Var_Type *mesh_vars = NULL;
   double *mesh_apriori_vert_trop = NULL;
   double *mesh_vert_strat = NULL;
   double *mesh_vert_trop = NULL;
   int num_steps, num_xtrack;
   int status = -1;

   if (0 != init_params (cfg, &params))
     return -1;

   /* Read mesh for calculations */
   if (0 != init_mesh (cfg, &mesh))
     goto free_and_return;

   /* Read scan grid */
   if (NULL == (st = scan_read_grids (num_files, files)))
     goto free_and_return;

   /* Prepare for regridding scan variables onto mesh grid */
   if (NULL == (r_mesh = scan_init_regrid (st, &mesh)))
     goto free_and_return;
   Pixel_regrid_get_srcdims (r_mesh, &num_steps, &num_xtrack);

   /* Read product variables on scan grid */
   if (0 != __scan_enable_testdata (st, params.testdata_file, r_mesh))
     goto free_and_return;
   if (NULL == (sv = scan_vars_alloc (num_steps, num_xtrack)))
     goto free_and_return;
   if (0 != scan_vars_read (st, sv))
     goto free_and_return;

   /* Regrid product variables onto mesh grid */
   if (NULL == (mesh_vars = regrid_product_vars (sv, &mesh, r_mesh)))
     goto free_and_return;

   /* Read apriori tropospheric vertical column and map to mesh grid */
   mesh_apriori_vert_trop =
     read_apriori_vertical_trop_column (params.apriori_trop_file, &mesh);
   if (NULL == mesh_apriori_vert_trop)
     goto free_and_return;

   /* Estimate stratospheric vertical column */
   mesh_vert_strat = vertical_strat_column (&mesh, mesh_vars, mesh_apriori_vert_trop,
                                            params.trop_thresh);
   if (mesh_vert_strat == NULL)
     goto free_and_return;

   /* Filter stratospheric vertical column estimate */
   if (0 != filter_vert_strat (&mesh, mesh_vert_strat, cfg))
     goto free_and_return;

   /* Compute tropospheric vertical column */
   mesh_vert_trop = vertical_trop_column (&mesh, mesh_vars, mesh_vert_strat);
   if (mesh_vert_trop == NULL)
     goto free_and_return;

   /* Write strat/trop variables to corresponding granule files */
   if (0 != write_split (st, sv, r_mesh, mesh_vert_trop, mesh_vert_strat))
     goto free_and_return;

   status = 0;
free_and_return:
   free_params (&params);
   scan_free (st);
   scan_vars_free (sv);
   Pixel_close_regrid (r_mesh);
   free_var_type (mesh_vars);
   FREE(mesh_apriori_vert_trop);
   FREE(mesh_vert_strat);
   FREE(mesh_vert_trop);

   return status;
}
