#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tell.h>
#include <netcdf.h>
#include <tio.h>
#include <tio_template.h>

#include "proj.h"
#include "scan.h"

/* Valid longitude, latitude values */
#define INVALID_LONGITUDE(b) ((0 == isfinite(b)) || ((b) < -180.0) || (360.0 < (b)))
#define INVALID_LATITUDE(b)  ((0 == isfinite(b)) || ((b) <  -90.0) || ( 90.0 < (b)))

typedef struct Granule_Type Granule_Type;

struct Granule_Type
{
   char *file;
   double *lon_bounds;
   double *lat_bounds;
   int *steps;
   int num_xtrack;
   int num_steps;
};

typedef struct  /* FIXME */
{
   char *filename;
   Pixel_Regrid_Type *regrid_obj;
}
Testdata_Info_Type;

struct Scan_Type
{
   Granule_Type **granules;
   int num_granules;
   int min_step;
   int max_step;
   int max_num_xtrack;
   Testdata_Info_Type __t;  /* FIXME */
};

static void granule_free (Granule_Type *gr)
{
   if (NULL == gr)
     return;

   FREE(gr->file);
   FREE(gr->steps);
   FREE(gr->lon_bounds);
   FREE(gr->lat_bounds);
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

static int read_pixel_vertices (Granule_Type *gr, int ncid)
{
   TIO_Var_Info_Type vi;
   int i, grp, num_pixels, start[3], count[3];
   size_t len_bounds;
   double nan_value = nan("");
   double *lonb, *latb;

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_XTRACK, &vi))
     return -1;
   gr->num_xtrack = vi.dimlens[0];

   if (-1 == TIO_inq_var (ncid, TEMPO_DIM_STEP, &vi))
     return -1;
   gr->num_steps = vi.dimlens[0];

   if (NULL == (gr->steps = (int *)MALLOC (gr->num_steps * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start[0] = 0;
   count[0] = gr->num_steps;
   if (-1 == TIO_get_var_section (ncid, TEMPO_DIM_STEP,
                                  start, count, TIO_INT, gr->steps))
     return -1;

   num_pixels = gr->num_steps * gr->num_xtrack;
   len_bounds = 4 * num_pixels * sizeof(double);

   if ((NULL == (gr->lon_bounds = (double *) MALLOC (len_bounds)))
       || (NULL == (gr->lat_bounds = (double *) MALLOC (len_bounds))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

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

static Granule_Type *granule_init (const char *file)
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

   if (0 != read_pixel_vertices (gr, ncid))
     {
        (void) TIO_close (ncid);
        granule_free (gr);
        return NULL;
     }
   (void) TIO_close (ncid);

   return gr;
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
   FREE(st->__t.filename);
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

Scan_Type *scan_read_grids (int num_files, char **files)
{
   Scan_Type *st = NULL;
   int i;

   if (NULL == (st = new_scan (num_files)))
     return NULL;

   st->min_step = 0;
   st->max_step = 0;
   st->max_num_xtrack = 0;

   for (i = 0; i < st->num_granules; i++)
     {
        Granule_Type *gr = granule_init (files[i]);
        int j;

        if (gr == NULL)
          {
             scan_free (st);
             return NULL;
          }
        st->granules[i]  = gr;

        /* record scan dimensions */
        if (gr->num_xtrack > st->max_num_xtrack)
          st->max_num_xtrack = gr->num_xtrack;
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

   if (-1 == Pixel_list_pack (lonlat, gr->lon_bounds, gr->lat_bounds,
                              num_pixels, 4, gr->steps, gr->num_xtrack))
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

   if (-1 == Pixel_list_pack (eqarea, albers_x_bounds, albers_y_bounds,
                              num_pixels, 4, gr->steps, gr->num_xtrack))
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
   Pixel_regrid_grow_srcdims (r, st->max_step, st->max_num_xtrack-1);

   for (i = 0; i < st->num_granules; i++)
     {
        if (0 != find_granule_overlaps (r, st->granules[i]))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: unexpected error", __func__);
             goto free_and_return;
          }
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

int __scan_enable_testdata (Scan_Type *st, const char *filename,
                            Pixel_Regrid_Type *r)
{
   st->__t.regrid_obj = r;

   if (NULL == (st->__t.filename = strdup (filename)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   return 0;
}

static double *read_testdata_value (const char *file, const char *var_name)
{
   TIO_Var_Info_Type vi;
   double *var = NULL;
   size_t num_pixels;
   int ncid, start[2], count[2];

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   if (-1 == TIO_inq_var (ncid, var_name, &vi))
     goto free_and_return;

   num_pixels = vi.dimlens[0] * vi.dimlens[1];
   if (NULL == (var = (double *)MALLOC (num_pixels * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = vi.dimlens[0];
   count[1] = vi.dimlens[1];
   if (0 != TIO_get_var_section (ncid, var_name, start, count, TIO_DOUBLE, var))
     goto free_and_return;

   (void) TIO_close (ncid);
   return var;

free_and_return:
   (void) TIO_close (ncid);
   FREE(var);
   return NULL;
}

/* If pvar != NULL, use the space it points to.  Otherwise, allocate space
 * and return a pointer to the allocated storage */
static double *read_var (const Scan_Type *st, const char *var_name, double *pvar)
{
   int num_xtrack = st->max_num_xtrack;
   int num_step = st->max_step + 1;
   int num_pixels = num_xtrack * num_step;
   double *mesh_value = NULL;
   double *var = NULL;
   int *mesh_mask = NULL;
   double nan_value = nan("");
   double fill_value = nan_value;
   int status = -1;

   if (pvar == NULL)
     {
        if (NULL == (var = (double *) MALLOC (num_pixels * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return NULL;
          }
     }
   else var = pvar;

   /* If we're getting variable values by regridding test data,
    * we'll assume that the test data mesh parameters are the same
    * as the the target mesh parameters (so we can use the existing
    * regrid object) */
   if (NULL == (mesh_value = read_testdata_value (st->__t.filename, var_name)))
        goto free_and_return;

   /* FIXME: prototype performs this transformation */
   if (0 == strcmp (var_name, "slant_column"))
     {
        int i;
        for (i = 0; i < num_pixels; i++)
          {
             if (mesh_value[i] == 0.0)
               mesh_value[i] = nan_value;
          }
     }

   /* Note that regridding the test data means that we lose some of the
    * bad pixels that appeared on the input grid -- simply because there
    * isn't a one-to-one correspondence between the original test data grid
    * and the scan grid.
    */
   if (0 != Pixel_regrid_from_mesh (st->__t.regrid_obj, mesh_mask, fill_value,
                                    mesh_value, var))
     goto free_and_return;

   status = 0;
free_and_return:
   FREE(mesh_value);

   if (status != 0)
     {
        if (pvar == NULL) FREE(var);
        var = NULL;
     }

   return var;
}

void scan_vars_free (Scan_Vars_Type *sv)
{
   if (sv == NULL)
     return;
   FREE(sv->slant_column);
   FREE(sv->amf_trop);
   FREE(sv->amf_strat);
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
       || (NULL == (sv->data_quality_flag = (int *)MALLOC (num_pixels * sizeof(int))))
      )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        scan_vars_free (sv);
        return NULL;
     }

   sv->num_steps = num_steps;
   sv->num_xtrack = num_xtrack;

   return sv;
}

int scan_vars_read (const Scan_Type *st, Scan_Vars_Type *sv)
{
   int i, num_pixels;
   double *slant_column;
   int *dqf;

   if ((NULL == read_var (st, "slant_column", sv->slant_column))
       || (NULL == read_var (st, "amf_trop", sv->amf_trop))
       || (NULL == read_var (st, "amf_strat", sv->amf_strat)))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: reading slant_column, amf_trop, amf_strat",
                     __func__);
        return -1;
     }

   dqf = sv->data_quality_flag;
   slant_column = sv->slant_column;
   num_pixels = sv->num_steps * sv->num_xtrack;

   for (i = 0; i < num_pixels; i++)
     {
        dqf[i] = (0 == isfinite (slant_column[i]));
     }

   return 0;
}

static int write_granule_vars (const Granule_Type *gr, double fill_value,
                               const double *vtrop_gr, const double *vstrat_gr)
{
   TIO_Var_Info_Type vi;
   const char vtrop_name[] = "trop_column_amount";
   const char vstrat_name[] = "strat_column_amount";
   int start[2], count[2], varid_trop, varid_strat;
   int i, num_steps, num_xtrack;
   int ncid, grp, status = -1;

   if (0 != TIO_open (gr->file, NC_WRITE, &ncid))
     return -1;

   if (0 != TIO_inq_grp (ncid, "product", &grp))
     goto close_and_return;

   if (-1 == TIO_inq_var (grp, "column_amount", &vi))
     goto close_and_return;

   if ((0 != TIO_def_var (grp, vtrop_name, TIO_DOUBLE, vi.ndims, vi.dimids, &varid_trop))
       || (0 != TIO_def_var (grp, vstrat_name, TIO_DOUBLE, vi.ndims, vi.dimids, &varid_strat)))
     {
        goto close_and_return;
     }

   if ((0 != TIO_def_var_fill (grp, varid_trop, 0, &fill_value))
       || (0 != TIO_def_var_fill (grp, varid_strat, 0, &fill_value)))
     {
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

        if ((0 != TIO_put_var_section (grp, vtrop_name, start, count, TIO_DOUBLE, vtrop_i))
            || (0 != TIO_put_var_section (grp, vstrat_name, start, count, TIO_DOUBLE, vstrat_i)))
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
     }

   return 0;
}
