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

#include "proj.h"
#include "scan.h"

static int _pEpoch_Set = 0;

/* Valid longitude, latitude values */
#define INVALID_LONGITUDE(b) ((0 == isfinite(b)) || ((b) < -180.0) || (360.0 < (b)))
#define INVALID_LATITUDE(b)  ((0 == isfinite(b)) || ((b) <  -90.0) || ( 90.0 < (b)))

typedef struct Granule_Type Granule_Type;

struct Granule_Type
{
   char *file;
   double *lon_bounds;
   double *lat_bounds;
   double *slant_column;
   double *amf_trop;
   double *amf_strat;
   double *vert_strat;
   double tstart;
   double tend;
   int *data_quality_flag;
   int *steps;
   int *xtrack;
   int num_xtrack;
   int num_steps;
};

struct Scan_Type
{
   Granule_Type **granules;
   int num_granules;
   int min_step;
   int max_step;
   int max_xtrack;
};

typedef struct
{
   double trop_thresh;
}
Params_Type;

static void granule_free (Granule_Type *gr)
{
   if (NULL == gr)
     return;

   FREE(gr->file);
   FREE(gr->steps);
   FREE(gr->xtrack);
   FREE(gr->lon_bounds);
   FREE(gr->lat_bounds);
   FREE(gr->slant_column);
   FREE(gr->amf_trop);
   FREE(gr->amf_strat);
   FREE(gr->vert_strat);
   FREE(gr->data_quality_flag);
   FREE(gr);
}

static Granule_Type *granule_new (const char *file)
{
   Granule_Type *gr = NULL;

   if (file == NULL)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: got NULL pointer", __func__);
        return NULL;
     }

   if (NULL == (gr = (Granule_Type *)MALLOC (sizeof *gr)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)gr, 0, sizeof *gr);

   if (NULL == (gr->file = strdup (file)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        granule_free (gr);
        return NULL;
     }

   return gr;
}

static int granule_alloc_data_arrays (Granule_Type *gr)
{
   int num_pixels;
   size_t len_bounds, len_doubles;

   if ((NULL == (gr->steps = (int *)MALLOC (gr->num_steps * sizeof(int))))
       || (NULL == (gr->xtrack = (int *)MALLOC (gr->num_xtrack * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   num_pixels = gr->num_steps * gr->num_xtrack;
   len_doubles = num_pixels * sizeof(double);
   len_bounds = 4 * len_doubles;

   if ((NULL == (gr->lon_bounds = (double *) MALLOC (len_bounds)))
       || (NULL == (gr->lat_bounds = (double *) MALLOC (len_bounds)))
       || (NULL == (gr->slant_column = (double *) MALLOC (len_doubles))
       || (NULL == (gr->amf_trop = (double *) MALLOC (len_doubles)))
       || (NULL == (gr->amf_strat = (double *) MALLOC (len_doubles)))
       || (NULL == (gr->vert_strat = (double *) MALLOC (len_doubles)))
       || (NULL == (gr->data_quality_flag = (int *) MALLOC (num_pixels * sizeof(int)))))
      )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   return 0;
}

static int read_pixel_vertices (Granule_Type *gr, int ncid)
{
   TIO_Var_Info_Type vi;
   int i, grp, num_pixels, start[3], count[3];
   double nan_value = nan("");
   double *lonb, *latb;

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_XTRACK, &vi))
     return -1;
   gr->num_xtrack = vi.dimlens[0];

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_STEP, &vi))
     return -1;
   gr->num_steps = vi.dimlens[0];

   if (0 != granule_alloc_data_arrays (gr))
     return -1;

   start[0] = 0;
   count[0] = gr->num_steps;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_STEP,
                                  start, count, TIO_INT, gr->steps))
     return -1;

   start[0] = 0;
   count[0] = gr->num_xtrack;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_XTRACK,
                                  start, count, TIO_INT, gr->xtrack))
     return -1;

   /* read lon/lat bounds arrays */

   if (-1 == TIO_inq_grp (ncid, "geolocation", &grp))
     return -1;

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = gr->num_steps;
   count[1] = gr->num_xtrack;
   count[2] = 4;

   if ((-1 == TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS,
                                   start, count, TIO_DOUBLE, gr->lon_bounds))
       || (-1 == TIO_get_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS,
                                      start, count, TIO_DOUBLE, gr->lat_bounds)))
     {
        return -1;
     }

   num_pixels = gr->num_steps * gr->num_xtrack;

   /* filter invalid values */
   lonb = gr->lon_bounds;
   latb = gr->lat_bounds;
   for (i = 0; i < 4 * num_pixels; i++)
     {
        double lonb_i = lonb[i];
        double latb_i = latb[i];
        if (INVALID_LONGITUDE(lonb_i)) lonb[i] = nan_value;
        if (INVALID_LATITUDE(latb_i)) latb[i] = nan_value;
     }

   return 0;
}

static void dbl_replace_fill_with_nan (double *a, int n, double fill_value)
{
   double nan_value = nan("");
   int i;
   for (i = 0; i < n; i++)
     {
        if (a[i] == fill_value) a[i] = nan_value;
     }
}

static int read_dbl_and_replace_fill (int grp, const char *name,
                                      int *count, double *value)
{
   TIO_Var_Info_Type info;
   int no_fill, start[] = {0, 0};
   double fill_value;

   if (0 != TIO_get_var_section (grp, name,
                                 start, count, TIO_DOUBLE, value))
     return -1;

   if (0 != TIO_inq_var (grp, name, &info))
     return -1;

   /* initialize to NaN, then read and see if the initial value changed */
   fill_value = nan("");

   if (info.type == TIO_FLOAT)
     {
        float float_fill = nan("");
        if (0 != TIO_inq_var_fill (grp, info.varid, &no_fill, &float_fill))
          return -1;
        fill_value = (double) float_fill;
     }
   else if (info.type == TIO_DOUBLE)
     {
        if (0 != TIO_inq_var_fill (grp, info.varid, &no_fill, &fill_value))
          return -1;
     }
   else
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported data type = %d",
                     __func__, info.type);
        return -1;
     }

   /* If the initial value didn't change, it appears that
    * there's no fill value to filter */
   if (0 == isnan(fill_value))
     {
        int num_pixels = count[0] * count[1];
        dbl_replace_fill_with_nan (value, num_pixels, fill_value);
     }

   return 0;
}

typedef struct
{
   double *psurf;
   double *ap;
   double *bp;
   int num_pressures;
   int num_pixels;
}
Pressure_Param_Type;

static void free_pressure_params (Pressure_Param_Type *p)
{
   if (p == NULL)
     return;
   FREE(p->psurf);
   FREE(p->ap);
   FREE(p->bp);
   FREE(p);
}

static Pressure_Param_Type *alloc_pressure_params (int num_pixels, int num_pressures)
{
   Pressure_Param_Type *p = NULL;

   if (NULL == (p = (Pressure_Param_Type *)MALLOC (sizeof *p)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)p, 0, sizeof *p);

   if ((NULL == (p->psurf = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (p->ap = (double *)MALLOC (num_pressures * sizeof(double))))
       || (NULL == (p->bp = (double *)MALLOC (num_pressures * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_pressure_params (p);
        return NULL;
     }

   p->num_pixels = num_pixels;
   p->num_pressures = num_pressures;

   return p;
}

static Pressure_Param_Type *read_pressure_params (const Granule_Type *gr, int ncid)
{
   TIO_Var_Info_Type info = {0};
   Pressure_Param_Type *p = NULL;
   const char *name = "surface_pressure";
   int grp, start[2], count[2], xtype, num_pressures, num_pixels;

   if (-1 == TIO_inq_grp (ncid, "support_data", &grp))
     return NULL;

   if ((0 != TIO_inq_var (grp, name, &info))
       || (0 != TIO_inq_att (grp, info.varid, "Eta_A", &xtype, &num_pressures)))
     return NULL;

   num_pixels = gr->num_steps * gr->num_xtrack;

   if (NULL == (p = alloc_pressure_params (num_pixels, num_pressures)))
     return NULL;

   start[0] = 0;
   start[1] = 0;
   count[0] = gr->num_steps;
   count[1] = gr->num_xtrack;

   if ((0 != TIO_get_var_section (grp, name, start, count, TIO_DOUBLE, p->psurf))
       || (0 != TIO_get_att (grp, info.varid, "Eta_A", TIO_DOUBLE, p->ap))
       || (0 != TIO_get_att (grp, info.varid, "Eta_B", TIO_DOUBLE, p->bp)))
     {
        free_pressure_params (p);
        return NULL;
     }

   return p;
}

static int find_tropopause (int pix, double ptrop, const Pressure_Param_Type *pt)
{
   double *ap = pt->ap;
   double *bp = pt->bp;
   double p0 = pt->psurf[pix];
   int n = pt->num_pressures;

   /* Is it faster to compute only the upper part of the
    * profile and search it linearly, or is it better to
    * compute the whole profile, and then do a binary search?
    * Let's try the linear approach first */
   while (n-- > 0)
     {
        double p = ap[n] + bp[n] * p0;
        if (p > ptrop)
          return n;
     }

   tell_verror (TELL_RUNTIME_ERROR,
                "%s: cannot find tropopause for pixel %d: P(trop) = %g",
                __func__, pix, ptrop);

   return -1;
}

static int compute_vstrat_from_file_data (Granule_Type *gr, int ncid,
                                          const Params_Type *params)
{
   Pressure_Param_Type *pt = NULL;
   int i, grp, start[3], count[3], dimid_levels, num_pixels;
   size_t dimlen_levels;
   double nan_value = nan("");
   double trop_thresh = params->trop_thresh;
   double *tropopause_pressure = NULL;
   double *gas_profile = NULL;
   int status = -1;

   if (0 != TIO_inq_dim (ncid, "swt_level", &dimid_levels, &dimlen_levels))
     return -1;

   if (-1 == TIO_inq_grp (ncid, "support_data", &grp))
     return -1;

   num_pixels = gr->num_steps * gr->num_xtrack;

   if  ((NULL == (tropopause_pressure = (double *)MALLOC (num_pixels * sizeof(double))))
        || (NULL == (gas_profile = (double *)MALLOC (num_pixels * dimlen_levels * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   count[0] = gr->num_steps;
   count[1] = gr->num_xtrack;
   count[2] = 0;

   if (0 != read_dbl_and_replace_fill (grp, "tropopause_pressure", count, tropopause_pressure))
     goto free_and_return;

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = gr->num_steps;
   count[1] = gr->num_xtrack;
   count[2] = dimlen_levels;

   if (0 != TIO_get_var_section (grp, "gas_profile", start, count, TIO_DOUBLE, gas_profile))
     goto free_and_return;

   if (NULL == (pt = read_pressure_params (gr, ncid)))
     goto free_and_return;

   for (i = 0; i < num_pixels; i++)
     {
        double *gas_profile_i = gas_profile + i * dimlen_levels;
        double vtrop_apriori, trop_slant;
        int k, ktrop;

        gr->vert_strat[i] = nan_value;

        if ((gr->data_quality_flag[i] != 0)
            || isnan(gr->amf_strat[i])
            || isnan(gr->slant_column[i])
            || isnan(tropopause_pressure[i]))
          {
             continue;
          }

        if ((ktrop = find_tropopause (i, tropopause_pressure[i], pt)) < 0)
          goto free_and_return;

        vtrop_apriori = 0.0;
        for (k = 0; k < ktrop; k++)
          {
             vtrop_apriori += gas_profile_i[k];
          }

        trop_slant = vtrop_apriori * gr->amf_trop[i];

        if (trop_slant < trop_thresh * gr->amf_strat[i])
          {
             gr->vert_strat[i] = (gr->slant_column[i] - trop_slant) / gr->amf_strat[i];
          }
     }

   status = 0;
free_and_return:
   free_pressure_params (pt);
   FREE(tropopause_pressure);
   FREE(gas_profile);

   return status;
}

static int read_data_arrays (Granule_Type *gr, int ncid)
{
   int grp, start[3], count[3];

   start[0] = 0;
   start[1] = 0;
   count[0] = gr->num_steps;
   count[1] = gr->num_xtrack;

   if (-1 == TIO_inq_grp (ncid, "product", &grp))
     return -1;

   if (0 != TIO_get_var_section (grp, "main_data_quality_flag",
                                 start, count, TIO_INT, gr->data_quality_flag))
     return -1;

   if (-1 == TIO_inq_grp (ncid, "support_data", &grp))
     return -1;

   if (0 != read_dbl_and_replace_fill (grp, "fitted_slant_column", count, gr->slant_column))
     return -1;

   if (0 != read_dbl_and_replace_fill (grp, "amf_troposphere", count, gr->amf_trop))
     return -1;

   if (0 != read_dbl_and_replace_fill (grp, "amf_stratosphere", count, gr->amf_strat))
     return -1;

   if (1)  /* FIXME - why are input slant columns < 0? */
   {
      double *slant_column = gr->slant_column;
      double nan_value = nan("");
      int i, num_pixels = gr->num_steps * gr->num_xtrack;
      for (i = 0; i < num_pixels; i++)
        {
           double slant_column_i = slant_column[i];
           if (slant_column_i < 0)
             {
                slant_column[i] = nan_value;
             }
        }
   }

   return 0;
}

static Granule_Type *granule_init (const char *file, const Params_Type *params)
{
   Granule_Type *gr = NULL;
   int ncid;

   if (NULL == (gr = granule_new (file)))
     return NULL;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     {
        granule_free (gr);
        return NULL;
     }

   if (0 == _pEpoch_Set)
     {
        if (0 != tio_use_file_epoch (ncid))
          goto free_and_return;
        _pEpoch_Set = 1;
     }

   if ((0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &gr->tstart))
       || (0 != TIO_get_att (ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, &gr->tend)))
     goto free_and_return;

   if (0 != read_pixel_vertices (gr, ncid))
     goto free_and_return;

   if (0 != read_data_arrays (gr, ncid))
     goto free_and_return;

   if (0 != compute_vstrat_from_file_data (gr, ncid, params))
     goto free_and_return;

   (void) TIO_close (ncid);
   return gr;

free_and_return:
   (void) TIO_close (ncid);
   granule_free (gr);
   return NULL;
}

void scan_free (Scan_Type *st)
{
   if (st == NULL)
     return;

   if (st->granules)
     {
        int i;
        for (i = 0; i < st->num_granules; i++)
          {
             granule_free (st->granules[i]);
          }
        FREE(st->granules);
     }
   FREE(st);
}

static Scan_Type *new_scan (int num_files)
{
   Scan_Type *st = NULL;
   size_t len_granule_array;

   if (NULL == (st = (Scan_Type *)MALLOC (sizeof *st)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)st, 0, sizeof *st);

   len_granule_array = num_files * sizeof (Granule_Type *);

   if (NULL == (st->granules = (Granule_Type **) MALLOC (len_granule_array)))
     {
        FREE(st);
        return NULL;
     }
   memset ((char *)st->granules, 0, len_granule_array);

   st->num_granules = num_files;

   return st;
}

static int init_params (config_t *cfg, Params_Type *params)
{
   if (CONFIG_TRUE != config_lookup_float (cfg, "trop_thresh", &params->trop_thresh))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading trop_thresh in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

Scan_Type *scan_read_granules (int num_files, char **files, config_t *cfg)
{
   Scan_Type *st = NULL;
   Params_Type params = {0};
   int i;

   if (0 != init_params (cfg, &params))
     return NULL;

   if (NULL == (st = new_scan (num_files)))
     return NULL;

   st->min_step = 0;
   st->max_step = 0;
   st->max_xtrack = 0;

   for (i = 0; i < st->num_granules; i++)
     {
        Granule_Type *gr = granule_init (files[i], &params);
        int j;

        if (gr == NULL)
          {
             scan_free (st);
             return NULL;
          }
        st->granules[i]  = gr;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "read %s", files[i]);

        /* record scan dimensions */
        for (j = 0; j < gr->num_xtrack; j++)
          {
             if (gr->xtrack[j] > st->max_xtrack)
               st->max_xtrack = gr->xtrack[j];
          }
        for (j = 0; j < gr->num_steps; j++)
          {
             if (gr->steps[j] > st->max_step)
               st->max_step = gr->steps[j];
             if (gr->steps[j] < st->min_step)
               st->min_step = gr->steps[j];
          }
     }

   return st;
}

int scan_time_interval (const Scan_Type *st, double *ptstart, double *ptend)
{
   double tbeg, tend;
   int i;

   tbeg = st->granules[0]->tstart;
   tend = st->granules[0]->tend;

   for (i = 1; i < st->num_granules; i++)
     {
        Granule_Type *gr = st->granules[i];
        if (gr->tstart < tbeg) tbeg = gr->tstart;
        if (gr->tend > tend) tend = gr->tstart;
     }

   *ptstart = tbeg;
   *ptend = tend;

   return 0;
}

static Pixel_List_Type *make_lonlat_pixel_list (const Granule_Type *gr)
{
   Pixel_List_Type *lonlat = NULL;
   int num_sides, num_pixels;

   /* Zero-length polygons will indicate lines of sight that
    * have invalid lon-lat coordinates -- usually because
    * they don't intersect the earth. */
   num_sides = 0;
   num_pixels = gr->num_steps * gr->num_xtrack;

   if (NULL == (lonlat = Pixel_list_new (num_pixels, num_sides))
       || (-1 == Pixel_list_use_src_index (lonlat)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initializing pixel list",
                     __func__);
        Pixel_list_free (lonlat);
        return NULL;
     }

   /* Assume all granules have the same xtrack size so that xtrack_dimlen = gr->num_xtrack */
   if (-1 == Pixel_list_pack (lonlat, gr->lon_bounds, gr->lat_bounds,
                              num_pixels, 4, gr->steps, gr->xtrack, gr->num_xtrack,
                              gr->num_xtrack))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: packing pixel list", __func__);
        Pixel_list_free (lonlat);
     }

   return lonlat;
}

static Pixel_List_Type *make_eqarea_pixel_list (const Granule_Type *gr)
{
   Pixel_List_Type *eqarea = NULL;
   double *bounds = NULL;
   double *lon_bounds, *lat_bounds;
   double *albers_x_bounds, *albers_y_bounds;
   int num_sides, num_pixels = gr->num_steps * gr->num_xtrack;
   size_t len_bounds = 4 * num_pixels * sizeof(double);
   int status = -1;

   if (NULL == (bounds = (double *)MALLOC (2 * len_bounds)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   lon_bounds = bounds;
   lat_bounds = bounds + 4*num_pixels;

   memcpy ((char *)lon_bounds, (char *)gr->lon_bounds, len_bounds);
   memcpy ((char *)lat_bounds, (char *)gr->lat_bounds, len_bounds);

   /* NOTE: coordinate projection is done in place, so after the call,
    * (lon_bounds, lat_bounds) [deg] is really Albers (x,y) [meters] */
    if (0 != proj_longlat_to_albers (lon_bounds, lat_bounds, 4*num_pixels))
     goto free_and_return;
   albers_x_bounds = lon_bounds;
   albers_y_bounds = lat_bounds;

   /* Zero-length polygons will indicate lines of sight that
    * have invalid lon-lat coordinates -- usually because
    * they don't intersect the earth. */
   num_sides = 0;

   if ((NULL == (eqarea = Pixel_list_new (num_pixels, num_sides)))
       || (-1 == Pixel_list_use_src_index (eqarea)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initializing pixel list",
                     __func__);
        goto free_and_return;
     }

   /* Assume all granules have the same xtrack size so that xtrack_dimlen = gr->num_xtrack */
   if (-1 == Pixel_list_pack (eqarea, albers_x_bounds, albers_y_bounds,
                              num_pixels, 4, gr->steps, gr->xtrack, gr->num_xtrack,
                              gr->num_xtrack))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: packing pixel list", __func__);
        goto free_and_return;
     }

   status = 0;
free_and_return:
   FREE(bounds);
   if (status)
     {
        Pixel_list_free (eqarea);
        eqarea = NULL;
     }
   return eqarea;
}

static int find_granule_overlaps (Pixel_Regrid_Type *r,
                                  const Granule_Type *gr)
{
   Pixel_List_Type *gr_eqarea = NULL;
   Pixel_List_Type *gr_lonlat = NULL;
   int status = -1;

   if ((NULL == (gr_lonlat = make_lonlat_pixel_list (gr)))
       || (NULL == (gr_eqarea = make_eqarea_pixel_list (gr))))
     {
        goto free_and_return;
     }

   if (Pixel_find_overlaps (r, gr_eqarea, gr_lonlat) < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unexpected error",
                     __func__);
        goto free_and_return;
     }

   status = 0;
free_and_return:
   Pixel_list_free (gr_eqarea);
   Pixel_list_free (gr_lonlat);
   return status;
}

Pixel_Regrid_Type *
scan_init_regrid (const Scan_Type *st, const Pixel_Grid_Param_Type *mesh)
{
   Pixel_Regrid_Type *r = NULL;
   Pixel_List_Type *st_mesh_eqarea = NULL;
   int i, status = -1;

   if (NULL == (st_mesh_eqarea = proj_pixel_list (mesh, proj_longlat_to_albers)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: failed generating equal area coordinates", __func__);
        return NULL;
     }

   if (NULL == (r = Pixel_open_regrid (mesh, st_mesh_eqarea)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: failed initializing regrid operation", __func__);
        goto free_and_return;
     }
   Pixel_regrid_grow_srcdims (r, st->max_step+1, st->max_xtrack+1);

   tell_vlog (TELL_MSGTYPE_INFO, 1, "map granule pixels to working mesh:");

   for (i = 0; i < st->num_granules; i++)
     {
        if (0 != find_granule_overlaps (r, st->granules[i]))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: unexpected error", __func__);
             goto free_and_return;
          }
        tell_vlog (TELL_MSGTYPE_INFO, 1, "finished: %s", st->granules[i]->file);
     }

   status = 0;
free_and_return:
   Pixel_list_free (st_mesh_eqarea);
   if (status)
     {
        Pixel_close_regrid (r);
        r = NULL;
     }

   return r;
}

void scan_vars_free (Scan_Vars_Type *sv)
{
   if (sv == NULL)
     return;
   FREE(sv->slant_column);
   FREE(sv->amf_trop);
   FREE(sv->amf_strat);
   FREE(sv->vert_strat);
   FREE(sv->data_quality_flag);
   FREE(sv);
}

Scan_Vars_Type *scan_vars_alloc (int num_steps, int num_xtrack)
{
   Scan_Vars_Type *sv = NULL;
   int num_pixels = num_steps * num_xtrack;

   if (NULL == (sv = (Scan_Vars_Type *)MALLOC (sizeof *sv)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sv, 0, sizeof *sv);

   if ((NULL == (sv->slant_column = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (sv->amf_strat = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (sv->amf_trop = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (sv->vert_strat = (double *)MALLOC (num_pixels * sizeof(double))))
       || (NULL == (sv->data_quality_flag = (int *)MALLOC (num_pixels * sizeof(int))))
      )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        scan_vars_free (sv);
        return NULL;
     }
   memset ((char *)sv->slant_column, 0, num_pixels * sizeof(double));
   memset ((char *)sv->amf_strat, 0, num_pixels * sizeof(double));
   memset ((char *)sv->amf_trop, 0, num_pixels * sizeof(double));
   memset ((char *)sv->vert_strat, 0, num_pixels * sizeof(double));
   memset ((char *)sv->data_quality_flag, 0, num_pixels * sizeof(int));

   sv->num_steps = num_steps;
   sv->num_xtrack = num_xtrack;

   return sv;
}

static void copy_dbl_field (double *dest, int num_xtrack_dest,
                            const double *src, const Granule_Type *gr)
{
   size_t len = gr->num_xtrack * sizeof(double);
   int j;

   for (j = 0; j < gr->num_steps; j++)
     {
        const double *src_j = src + j * gr->num_xtrack;
        double *dest_j = dest + gr->steps[j] * num_xtrack_dest;
        memcpy ((char *)dest_j, (char *)src_j, len);
     }
}

static void copy_int_field (int *dest, int num_xtrack_dest,
                            const int *src, const Granule_Type *gr)
{
   size_t len = gr->num_xtrack * sizeof(int);
   int j;

   for (j = 0; j < gr->num_steps; j++)
     {
        const int *src_j = src + j * gr->num_xtrack;
        int *dest_j = dest + gr->steps[j] * num_xtrack_dest;
        memcpy ((char *)dest_j, (char *)src_j, len);
     }
}

int scan_vars_pack (const Scan_Type *st, Scan_Vars_Type *sv)
{
   int i, num_pixels;
   int *dqf;

   for (i = 0; i < st->num_granules; i++)
     {
        Granule_Type *gr = st->granules[i];
        copy_dbl_field (sv->slant_column, sv->num_xtrack, gr->slant_column, gr);
        copy_dbl_field (sv->amf_trop, sv->num_xtrack, gr->amf_trop, gr);
        copy_dbl_field (sv->amf_strat, sv->num_xtrack, gr->amf_strat, gr);
        copy_dbl_field (sv->vert_strat, sv->num_xtrack, gr->vert_strat, gr);
        copy_int_field (sv->data_quality_flag, sv->num_xtrack, gr->data_quality_flag, gr);
     }

   num_pixels = sv->num_steps * sv->num_xtrack;
   dqf = sv->data_quality_flag;

   /* FIXME: this seems crude -- change it? */
   for (i = 0; i < num_pixels; i++)
     {
        dqf[i] = (dqf[i] < 0) ? 1 : 0;
     }

   return 0;
}

typedef struct
{
   const char *name;
   const char *long_name;
   int varid;
}
Name_Type;

static int write_vertical_column_attributes (int grp, double fill_value, const char *units,
                                             const Name_Type *vcol)
{
   const char coordinate_string[] = "longitude latitude";

   if (0 != TIO_def_var_fill (grp, vcol->varid, 0, &fill_value))
     return -1;

   if ((0 != TIO_put_att (grp, vcol->varid, "coordinates", TIO_CHAR, sizeof(coordinate_string), coordinate_string))
       || (0 != TIO_put_att (grp, vcol->varid, "long_name", TIO_CHAR, 1+strlen(vcol->long_name), vcol->long_name))
       || (0 != TIO_put_att (grp, vcol->varid, "units", TIO_CHAR, 1 + strlen(units), units)))
     return -1;

   return 0;
}

static int write_granule_vars (const Granule_Type *gr, double fill_value,
                               const double *vtrop_gr, const double *vstrat_gr)
{
   TIO_Var_Info_Type vi;
   Name_Type vtrop =
     {
        .name = "vertical_column_troposphere",
        .long_name = "troposphere nitrogen dioxide vertical column",
        .varid = 0
     };
   Name_Type vstrat =
     {
        .name = "vertical_column_stratosphere",
        .long_name = "stratosphere nitrogen dioxide vertical column",
        .varid = 0
     };
   char units_slant_col[256] = {0};
   int start[2], count[2];
   int i, num_steps, num_xtrack, varid_slant_col;
   int ncid, grp_support, grp_product, status = -1;

   if (0 != TIO_open (gr->file, NC_WRITE, &ncid))
     return -1;

   if (0 != tio_history_append_cmdline (ncid))
     goto close_and_return;

   if (0 != TIO_inq_grp (ncid, "support_data", &grp_support))
     goto close_and_return;

   if ((0 != tio_inq_varid (grp_support, "fitted_slant_column", &varid_slant_col))
       || (0 != TIO_get_att (grp_support, varid_slant_col, "units", TIO_CHAR, units_slant_col)))
     goto close_and_return;

   if (0 != TIO_inq_grp (ncid, "product", &grp_product))
     goto close_and_return;

   if (-1 == TIO_inq_var (grp_product, "vertical_column_total", &vi))
     goto close_and_return;

   if (0 != tio_inq_varid (grp_product, vtrop.name, &vtrop.varid))
     {
        if (0 != TIO_def_var (grp_product, vtrop.name, TIO_DOUBLE, vi.ndims, vi.dimids, &vtrop.varid))
          goto close_and_return;
        if (0 != write_vertical_column_attributes (grp_product, fill_value, units_slant_col, &vtrop))
          goto close_and_return;
     }

   if (0 != tio_inq_varid (grp_product, vstrat.name, &vstrat.varid))
     {
        if (0 != TIO_def_var (grp_product, vstrat.name, TIO_DOUBLE, vi.ndims, vi.dimids, &vstrat.varid))
          goto close_and_return;
        if (0 != write_vertical_column_attributes (grp_product, fill_value, units_slant_col, &vstrat))
          goto close_and_return;
     }

   num_steps = vi.dimlens[0];
   num_xtrack = vi.dimlens[1];

   for (i = 0; i < num_steps; i++)
     {
        const double *vtrop_i = vtrop_gr + gr->steps[i] * num_xtrack;
        const double *vstrat_i = vstrat_gr + gr->steps[i] * num_xtrack;

        start[0] = i;
        start[1] = 0;
        count[0] = 1;
        count[1] = num_xtrack;

        if ((0 != TIO_put_var_section (grp_product, vtrop.name, start, count, TIO_DOUBLE, vtrop_i))
            || (0 != TIO_put_var_section (grp_product, vstrat.name, start, count, TIO_DOUBLE, vstrat_i)))
          {
             goto close_and_return;
          }
     }

   status = 0;
close_and_return:

   if (0 != TIO_close (ncid))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: closing %s", __func__, gr->file);
        return -1;
     }

   return status;
}

int scan_write_split (const Scan_Type *st, double fill_value,
                      const double *vtrop, const double *vstrat)
{
   int i;

   for (i = 0; i < st->num_granules; i++)
     {
        Granule_Type *gr = st->granules[i];
        if (0 != write_granule_vars (gr, fill_value, vtrop, vstrat))
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: writing strat/trop values for granule %s",
                          __func__, gr->file);
             return -1;
          }
        tell_vlog (TELL_MSGTYPE_INFO, 1, "wrote strat/trop split: %s", gr->file);
     }

   return 0;
}
