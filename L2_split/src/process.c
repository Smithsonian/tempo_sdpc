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
   double *vert_trop;
   double *vert_strat;
   int num_steps;
   int num_xtrack;
}
Split_Type;

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

   tell_vlog (TELL_MSGTYPE_INFO, 1, "working mesh longitude: [%7.3f, %7.3f, %4d]", mesh->xmin, mesh->xmax, mesh->nx);
   tell_vlog (TELL_MSGTYPE_INFO, 1, "working mesh latitude:  [%7.3f, %7.3f, %4d]", mesh->ymin, mesh->ymax, mesh->ny);

   return 0;
}

static double *map_vert_strat_to_mesh (Scan_Vars_Type *sv, const Pixel_Grid_Param_Type *mesh,
                                       const Pixel_Regrid_Type *r_mesh)
{
   double fill_value = nan("");
   double *vert_strat = NULL;
   int *mesh_mask = NULL;   /* FIXME? */
   int num_pixels, status = -1;

   /* Map from the Level 2 product scan grid to a uniform mesh grid
    * that will be used for computations.
    */
   num_pixels = mesh->nx * mesh->ny;
   if (NULL == (vert_strat = (double *)MALLOC (num_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if (0 != Pixel_regrid (r_mesh, mesh_mask, fill_value, sv->vert_strat, vert_strat, NULL))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mapping from scan grid to mesh grid", __func__);
        goto free_and_return;
     }

   status = 0;
free_and_return:
   if (status)
     {
        FREE(vert_strat);
        vert_strat = NULL;
     }

   return vert_strat;
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

static void free_split_type (Split_Type *split)
{
   if (split == NULL)
     return;
   FREE(split->vert_strat);
   FREE(split->vert_trop);
   FREE(split);
}

static Split_Type *alloc_split_type (int num_steps, int num_xtrack)
{
   Split_Type *split = NULL;
   int num_scan_pixels;

   if (NULL == (split = (Split_Type *)MALLOC (sizeof *split)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   split->num_steps = num_steps;
   split->num_xtrack = num_xtrack;

   num_scan_pixels = num_steps * num_xtrack;
   if ((NULL == (split->vert_trop = (double *) MALLOC (num_scan_pixels * sizeof(double))))
       || (NULL == (split->vert_strat = (double *) MALLOC (num_scan_pixels * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_split_type (split);
        return NULL;
     }

   return split;
}

static Split_Type *perform_split (const Scan_Vars_Type *sv,
                                  const Pixel_Regrid_Type *r_mesh,
                                  const double *mesh_vert_strat)
{
   int *dqf = sv->data_quality_flag;
   Split_Type *split = NULL;
   double *vert_strat = NULL;
   double *vert_trop = NULL;
   int *mesh_mask = NULL;   /* FIXME? */
   double nan_value = nan("");
   double fill_value = nan_value;
   int i, num_scan_pixels;

   if (NULL == (split = alloc_split_type (sv->num_steps, sv->num_xtrack)))
     return NULL;

   vert_strat = split->vert_strat;
   vert_trop = split->vert_trop;

   /* Map vert_strat column values from mesh back to Level 2 product scan grid */
   if (0 != Pixel_regrid_from_mesh (r_mesh, mesh_mask, fill_value, mesh_vert_strat, vert_strat))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: mapping from mesh grid to scan grid", __func__);
        free_split_type (split);
        return NULL;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "mapped filtered vertical strat column onto granule pixels");

   num_scan_pixels = sv->num_steps * sv->num_xtrack;

   /* Compute the tropospheric vertical column on the native scan pixel grid */
   for (i = 0; i < num_scan_pixels; i++)
     {
        if (sv->amf_trop[i] != 0)
          {
             vert_trop[i] = (sv->slant_column[i] -
                          vert_strat[i] * sv->amf_strat[i]) / sv->amf_trop[i];
          }
        else vert_trop[i] = nan_value;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "computed vertical trop column in granule pixels");

   /* Mask strat/trop variables using original data quality flag */
   for (i = 0; i < num_scan_pixels; i++)
     {
        if (dqf[i] != 0)
          {
             vert_strat[i] = nan_value;
             vert_trop[i] = nan_value;
          }
     }

   return split;
}

static int write_split (const Scan_Type *st, const Split_Type *split)
{
   double nan_value = nan("");
   double fill_value = nan_value;
   double fill_value_output = -1.0e30;  /* (should not be a NaN!) */
   int num_scan_pixels;

   num_scan_pixels = split->num_steps * split->num_xtrack;

   /* For output, replace fill_value, NaN, Inf with standard fill value */
   if (isnan(fill_value) || (fill_value_output != fill_value))
     {
        replace_fill (split->vert_strat, num_scan_pixels, fill_value, fill_value_output);
        replace_fill (split->vert_trop, num_scan_pixels, fill_value, fill_value_output);
     }

   /* Write strat/trop vertical columns on Level 2 product scan grid */
   if (0 != scan_write_split (st, fill_value_output, split->vert_trop, split->vert_strat))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: writing strat/trop vertical column values",
                     __func__);
        return -1;
     }

   return 0;
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
   Pixel_Grid_Param_Type mesh = {0};
   Scan_Type *st = NULL;
   Scan_Vars_Type *sv = NULL;
   Pixel_Regrid_Type *r_mesh = NULL;
   Split_Type *split = NULL;
   double *mesh_vert_strat = NULL;
   int num_steps, num_xtrack;
   int status = -1;

   /* Read mesh for calculations */
   if (0 != init_mesh (cfg, &mesh))
     goto free_and_return;

   /* Read scan granules, and compute initial estimate of
    * the stratospheric vertical column on the native scan grid.
    */
   if (NULL == (st = scan_read_granules (num_files, files, cfg)))
     goto free_and_return;

   /* Prepare for regridding scan variables onto mesh grid */
   if (NULL == (r_mesh = scan_init_regrid (st, &mesh)))
     goto free_and_return;
   Pixel_regrid_get_srcdims (r_mesh, &num_steps, &num_xtrack);

   if (NULL == (sv = scan_vars_alloc (num_steps, num_xtrack)))
     goto free_and_return;

   /* Pack product variables on scan grid-sized arrays */
   if (0 != scan_vars_pack (st, sv))
     goto free_and_return;

   /* Regrid stratospheric vertical column estimate onto mesh grid */
   if (NULL == (mesh_vert_strat = map_vert_strat_to_mesh (sv, &mesh, r_mesh)))
     goto free_and_return;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "regridded vertical strat column onto working mesh");

   /* Filter stratospheric vertical column estimate */
   if (0 != filter_vert_strat (&mesh, mesh_vert_strat, cfg))
     goto free_and_return;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "filtered vertical strat column (on working mesh)");

   /* Perform strat/trop separation on native scan grid */
   if (NULL == (split = perform_split (sv, r_mesh, mesh_vert_strat)))
     goto free_and_return;

   /* Write strat/trop variables to corresponding granule files */
   if (0 != write_split (st, split))
     goto free_and_return;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "done");

   status = 0;
free_and_return:
   scan_free (st);
   scan_vars_free (sv);
   Pixel_close_regrid (r_mesh);
   free_split_type (split);
   FREE(mesh_vert_strat);

   return status;
}
