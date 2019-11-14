#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <eccodes.h>
#include <proj_api.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include <tio.h>
#include <tell.h>

#define DEGTORAD (M_PI/180.0)
#define BUFSIZE 256

typedef struct
{
   double *v;
   size_t n;
}
Array_Type;

/* Lambert Conformal Conic Projection.
 * For now, this supports only the variant with a single standard parallel.
 */
typedef struct
{
   Array_Type x;     /* projection x coordinate grid [meters] */
   Array_Type y;     /* projection y coordinate grid [meters] */

   double xlim[2];   /* x range spanned by pixel boundaries [meters] */
   double ylim[2];   /* y range spanned by pixel boundaries [meters] */

   /* Proj.4 data structures for the relevant map projections */
   projPJ lcc;
   projPJ longlat;

   /* Angles in radians: */
   double latitude_of_projection_origin;
   double longitude_of_central_meridian;
   double standard_parallel;
   double latitude_of_first_grid_point;
   double longitude_of_first_grid_point;

   /* Lengths in meters: */
   double dx;
   double dy;
   double earth_radius;
}
LCC_Grid_Type;

typedef struct
{
   double *pressure_surface;      /* 2D array: pressure at ground or water surface [hPa] */
   double *pressure_tropopause;   /* 2D array: pressure at tropopause [hPa] */
   double *temperature_on_isobar; /* 3D array: 2D temperature on each of N isobaric surfaces [K] */
}
MET_Data_Type;

#define MFT_PRIVATE_DATA \
   LCC_Grid_Type *grid; \
   Array_Type isobars; \
   MET_Data_Type vars; \
   Array_Type var_on_isobar; \
   unsigned int flags;
#include "met.h"

static void array_free (Array_Type *at)
{
   if (at == NULL)
     return;
   FREE(at->v);
}

static int array_alloc (Array_Type *at, size_t n)
{
   if (at == NULL)
     return -1;

   if (NULL == (at->v = (double *)MALLOC (n * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   at->n = n;

   return 0;
}

static void lcc_grid_free (LCC_Grid_Type *grid)
{
   if (grid == NULL)
     return;
   array_free (&grid->x);
   array_free (&grid->y);
   pj_free (grid->lcc);
   pj_free (grid->longlat);
   FREE(grid);
}

/* wrapper to facilitate reading an int from a grib2 file, instead of a long  */
static int p_codes_get_int (codes_handle *h, const char *name, int *n)
{
   long nl;
   int err = codes_get_long (h, name, &nl);
   *n = nl;
   return err;
}

static LCC_Grid_Type *lcc_grid_new (int nx, int ny)
{
   LCC_Grid_Type *grid = NULL;

   if (NULL == (grid = (LCC_Grid_Type *)MALLOC (sizeof *grid)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)grid, 0, sizeof *grid);

   if ((0 != array_alloc (&grid->x, nx))
       || (0 != array_alloc (&grid->y, ny)))
     {
        lcc_grid_free (grid);
        return NULL;
     }

   return grid;
}

static void mdt_free (MET_Data_Type *mdt)
{
   if (mdt == NULL)
     return;
   FREE(mdt->pressure_surface);
   FREE(mdt->pressure_tropopause);
   FREE(mdt->temperature_on_isobar);
}

static void mft_close (Met_File_Type *mft)
{
   if (mft == NULL)
     return;

   lcc_grid_free (mft->grid);
   array_free (&mft->isobars);
   array_free (&mft->var_on_isobar);
   mdt_free (&mft->vars);
   FREE(mft);
}

static int get_earth_radius (codes_handle *h, double *earth_radius)
{
   int codes_err, earth_shape_code;

   /* GRIB2 files use an integer code system to select the shape of the earth
    * from a hard-coded list of options. Full support for this would be overkill
    * for our purpose.  Instead, we'll implement what we need and warn when
    * there's a mismatch.
    */
   if ((codes_err = p_codes_get_int (h, "shapeOfTheEarth", &earth_shape_code)) != 0)
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading shapeOfTheEarth from grib2 file", __func__);
        return -1;
     }

   if (earth_shape_code != 6)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: grib2 earth shape code = %d (expected 6)",
                     __func__, earth_shape_code);
        return -1;
     }

   /* This is the hard-coded earth radius value from options 6 in the GRIB2 tables.
    * Don't "fix" it unless you're sure you know what you're doing.
    */
   *earth_radius = 6371229.0;

   return 0;
}

static int init_projection (LCC_Grid_Type *g)
{
   /* a,b= ellipsoid axes [m]
    * lon_0, lat_0 = projection origin [deg]
    * lat_1 = latitude of first standard parallel [deg]
    * lat_2 = latitude of second standard parallel [deg]
    */
   const char lcc_fmt[] = "+proj=lcc +units=m +a=%f +b=%f +lon_0=%f +lat_0=%f +lat_1=%f +lat_2=%f";
   const char longlat_fmt[] = "+proj=longlat +ellps=sphere +a=%f";
   char ctl_lcc[BUFSIZE];
   char ctl_longlat[BUFSIZE];
   int len;

   memset (ctl_lcc, 0, BUFSIZE);
   len = snprintf (ctl_lcc, BUFSIZE, lcc_fmt, g->earth_radius, g->earth_radius,
                   g->longitude_of_central_meridian /DEGTORAD,
                   g->latitude_of_projection_origin /DEGTORAD,
                   g->standard_parallel /DEGTORAD,
                   g->standard_parallel /DEGTORAD);
   if (len >= BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: proj4 arg buffer too small", __func__);
        return -1;
     }

   if (NULL == (g->lcc = pj_init_plus (ctl_lcc)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: pj_init_plus(lcc) failed", __func__);
        return -1;
     }

   memset (ctl_longlat, 0, BUFSIZE);
   len = snprintf (ctl_longlat, BUFSIZE, longlat_fmt, g->earth_radius);
   if (len >= BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: proj4 arg buffer too small", __func__);
        return -1;
     }

   if (NULL == (g->longlat = pj_init_plus (ctl_longlat)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: pj_init_plus(longlat) failed", __func__);
        return -1;
     }

   return 0;
}

static LCC_Grid_Type *init_spatial_grid (const char *path)
{
   FILE *fp = NULL;
   LCC_Grid_Type *g = NULL;
   codes_handle *h = NULL;
   int codes_err, err_pj, nx, ny, num_pts, i;
   double x0, y0;
   int status = -1;

   if (NULL == (fp = fopen (path, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: open for reading: %s", __func__, path);
        return NULL;
     }

   if (NULL == (h = codes_handle_new_from_file (0, fp, PRODUCT_GRIB, &codes_err)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: codes_err=%d: unable to create handle from file: %s",
                     __func__, codes_err, path);
        return NULL;
     }

   if (((codes_err = p_codes_get_int (h, "Nx", &nx)) != 0)
       || ((codes_err = p_codes_get_int (h, "Ny", &ny)) != 0))
     {
        tell_verror(TELL_IO_READ_ERROR, "%s: codes_err=%d: reading grid dimensions",
                    __func__, codes_err);
        return NULL;
     }

   if (NULL == (g = lcc_grid_new (nx, ny)))
     return NULL;

   if (((codes_err = codes_get_double (h, "LaDInDegrees", &g->latitude_of_projection_origin)) != 0)
       || ((codes_err = codes_get_double (h, "LoVInDegrees", &g->longitude_of_central_meridian)) != 0)
       || ((codes_err = codes_get_double (h, "Latin1InDegrees", &g->standard_parallel)) != 0))
     {
        tell_verror(TELL_IO_READ_ERROR, "%s: codes_err=%d: reading projection parameters",
                    __func__, codes_err);
        goto return_error;
     }

   g->latitude_of_projection_origin *= DEGTORAD;
   g->longitude_of_central_meridian *= DEGTORAD;
   g->standard_parallel *= DEGTORAD;

   if (0 != get_earth_radius (h, &g->earth_radius))
     goto return_error;

   if (((codes_err = codes_get_double (h, "latitudeOfFirstGridPointInDegrees", &g->latitude_of_first_grid_point)) != 0)
       || ((codes_err = codes_get_double (h, "longitudeOfFirstGridPointInDegrees", &g->longitude_of_first_grid_point)) != 0)
       || ((codes_err = codes_get_double (h, "DxInMetres", &g->dx)) != 0)
       || ((codes_err = codes_get_double (h, "DyInMetres", &g->dy))) != 0)
     {
        tell_verror(TELL_IO_READ_ERROR, "%s: codes_err=%d: reading grid parameters",
                    __func__, codes_err);
        goto return_error;
     }

   g->latitude_of_first_grid_point *= DEGTORAD;
   g->longitude_of_first_grid_point *= DEGTORAD;

   if (0 != init_projection (g))
     goto return_error;

   x0 = g->longitude_of_first_grid_point;
   y0 = g->latitude_of_first_grid_point;
   num_pts = 1;

   if ((err_pj = pj_transform (g->longlat, g->lcc, num_pts, 1, &x0, &y0, NULL) != 0))
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(err_pj));
        goto return_error;
     }

   /* Lambert conformal conic projection grid is Cartesian, with fixed grid spacing.
    * It's assumed that each point is a pixel center, so the pixel boundaries extend
    * to +/- [dx/2, dy/2] from each point.
    */
   for (i = 0; i < nx; i++)
     {
        g->x.v[i] = x0 + i * g->dx;
     }
   g->xlim[0] = g->x.v[0]    - g->dx/2;
   g->xlim[1] = g->x.v[nx-1] + g->dx/2;

   for (i = 0; i < ny; i++)
     {
        g->y.v[i] = y0 + i * g->dy;
     }
   g->ylim[0] = g->y.v[0]    - g->dy/2;
   g->ylim[1] = g->y.v[ny-1] + g->dy/2;

   status = 0;
return_error:
   codes_handle_delete (h);
   fclose (fp);
   if (status)
     {
        lcc_grid_free (g);
        g = NULL;
     }
   return g;
}

static void free_string_array (char **a, size_t n)
{
   size_t i;
   if (a == NULL)
     return;
   for (i = 0; i < n; i++)
     {
        FREE(a[i]);
     }
   FREE(a);
}

static char **get_unique_string_values (codes_index *index, const char *name, size_t *num)
{
   char **s = NULL;
   int codes_err;

   if ((codes_err = codes_index_get_size (index, name, num)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err=%d)", __func__, codes_err);
        return NULL;
     }

   if (NULL == (s = (char **)MALLOC (*num * sizeof(char *))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if ((codes_err = codes_index_get_string (index, name, s, num)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err=%d)", __func__, codes_err);
        FREE(s);
        return NULL;
     }

   return s;
}

static long *get_unique_long_values (codes_index *index, const char *name, size_t *num)
{
   long *s = NULL;
   int codes_err;

   if ((codes_err = codes_index_get_size (index, name, num)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err=%d)", __func__, codes_err);
        return NULL;
     }

   if (NULL == (s = (long *)MALLOC (*num * sizeof(long))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   if ((codes_err = codes_index_get_long (index, name, s, num)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err=%d)", __func__, codes_err);
        FREE(s);
        return NULL;
     }

   return s;
}

static int dbl_compar (const void *va, const void *vb)
{
   double a = *(const double *)va;
   double b = *(const double *)vb;
   if (a < b)
     return -1;
   else if (a > b)
     return 1;
   else return 0;
}

/* This is ridiculously over-complicated for what it does, but
 * I see no other way to accomplish this with the eccodes interface.
 * If you know a better/simpler way, feel free to rewrite it!
 */
static int get_isobars (codes_index *index, Array_Type *isobars)
{
   char **shortnames = NULL;
   char **typeoflevel = NULL;
   long *levels = NULL;
   size_t sn, num_shortnames;
   size_t lev, num_levels;
   size_t ty, num_types;
   size_t num_isobars;
   int max_num_isobars;
   int codes_err=0, status = -1;

   if ((NULL == (shortnames = get_unique_string_values (index, "shortName", &num_shortnames)))
       || (NULL == (typeoflevel = get_unique_string_values (index, "typeOfLevel", &num_types)))
       || (NULL == (levels = get_unique_long_values (index, "level", &num_levels))))
     goto free_and_return;

   max_num_isobars = num_shortnames * num_types * num_levels;

   if (0 != array_alloc (isobars, max_num_isobars))
     goto free_and_return;

   num_isobars = 0;

   for (sn = 0; sn < num_shortnames; sn++)
     {
        if (0 != strcmp (shortnames[sn], "t"))
          continue;
        if ((codes_err = codes_index_select_string (index, "shortName", shortnames[sn])) != 0)
          goto free_and_return;

        for (ty = 0; ty < num_types; ty++)
          {
             if (0 != strcmp (typeoflevel[ty], "isobaricInhPa"))
               continue;
             if ((codes_err = codes_index_select_string (index, "typeOfLevel", typeoflevel[ty])) != 0)
               goto free_and_return;

             for (lev = 0; lev < num_levels; lev++)
               {
                  codes_handle *h;

                  codes_err = codes_index_select_long (index, "level", levels[lev]);
                  if ((codes_err != 0) && (codes_err != GRIB_END_OF_INDEX))
                    goto free_and_return;

                  if (NULL != (h = codes_handle_new_from_index (index, &codes_err)))
                    {
                       isobars->v[num_isobars++] = levels[lev];
                       codes_handle_delete (h);
                    }
               }
          }
     }

   isobars->n = num_isobars;

   /* This will ensure that we later read the temperature slabs in isobar order */
   qsort (isobars->v, num_isobars, sizeof(double), dbl_compar);

   status = 0;
free_and_return:
   free_string_array (shortnames, num_shortnames);
   free_string_array (typeoflevel, num_types);
   FREE(levels);

   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err=%d)", __func__, codes_err);
     }

   return status;
}

static int read_temperature_on_isobars (Met_File_Type *mft, codes_index *index)
{
   LCC_Grid_Type *g = mft->grid;
   codes_handle *h_iso = NULL;
   double *t_iso = NULL;
   size_t i, num_grid, num_isobars;
   int codes_err, status = -1;

   if (0 != get_isobars (index, &mft->isobars))
     goto return_status;

   num_isobars = mft->isobars.n;
   num_grid = g->x.n * g->y.n;

   if (0 != array_alloc (&mft->var_on_isobar, num_isobars))
     goto return_status;

   if (NULL == (t_iso = (double *)MALLOC (num_isobars * num_grid * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }
   mft->vars.temperature_on_isobar = t_iso;

   codes_err = codes_index_select_string (index, "shortName", "t");
   codes_err = codes_index_select_string (index, "typeOfLevel", "isobaricInhPa");

   for (i = 0; i < num_isobars; i++)
     {
        long lev = mft->isobars.v[i];

        if ((codes_err = codes_index_select_long (index, "level", lev)) != 0)
          goto return_status;

        h_iso = codes_handle_new_from_index (index, &codes_err);
        if ((codes_err != 0) || (h_iso == NULL))
          goto return_status;

        t_iso = mft->vars.temperature_on_isobar + i * num_grid;

        codes_err = codes_get_double_array (h_iso, "values", t_iso, &num_grid);
        codes_handle_delete (h_iso);
        h_iso = NULL;
        if (codes_err != 0)
          goto return_status;
     }

   status = 0;
return_status:
   codes_handle_delete (h_iso);
   return status;
}

static double *read_slab (codes_index *index, char *shortname, char *typeoflevel,
                          size_t num_grid)
{
   codes_handle *h = NULL;
   double *values = NULL;
   int codes_err = -1;

   if ((0 != (codes_err = codes_index_select_string (index, "shortName", shortname)))
       || (0 != (codes_err = codes_index_select_string (index, "typeOfLevel", typeoflevel)))
       || (0 != (codes_err = codes_index_select_long (index, "level", 0L))))
     goto return_status;

   if (NULL == (h = codes_handle_new_from_index (index, &codes_err)))
     goto return_status;

   if (NULL == (values = (double *)MALLOC (num_grid * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   codes_err = codes_get_double_array (h, "values", values, &num_grid);

return_status:
   codes_handle_delete (h);
   if (codes_err)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err = %d)", __func__, codes_err);
        FREE(values);
        values = NULL;
     }

   return values;
}

static int read_forecast_vars (Met_File_Type *mft, const char *path)
{
   LCC_Grid_Type *g = mft->grid;
   codes_index *index = NULL;
   double *p = NULL;
   size_t i, num_grid;
   int codes_err, status = -1;

   num_grid = g->x.n * g->y.n;

   index = codes_index_new (0, "shortName,level,typeOfLevel", &codes_err);
   codes_err = codes_index_add_file (index, path);

   if (mft->flags & MET_READ_PRESSURE_SURFACE)
     {
        if (NULL == (p = read_slab (index, "sp", "surface", num_grid)))
          goto return_status;
        /* convert Pa to hPa */
        for (i = 0; i < num_grid; i++)
          p[i] /= 100.0;
        mft->vars.pressure_surface = p;
     }

   if (mft->flags & MET_READ_PRESSURE_TROPOPAUSE)
     {
        if (NULL == (p = read_slab (index, "pres", "tropopause", num_grid)))
          goto return_status;
        /* convert Pa to hPa */
        for (i = 0; i < num_grid; i++)
          p[i] /= 100.0;
        mft->vars.pressure_tropopause = p;
     }

   if (mft->flags & MET_READ_TEMPERATURE_ON_ISOBARS)
     {
        if (0 != read_temperature_on_isobars (mft, index))
          goto return_status;
     }

   status = 0;
return_status:
   codes_index_delete (index);

   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed (codes_err=%d)",
                     __func__, codes_err);
     }

   return status;
}

static int read_grib2_file (Met_File_Type *mft, const char *path)
{
   if (NULL == (mft->grid = init_spatial_grid (path)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error initializing spatial grid from file: %s",
                     __func__, path);
        return -1;
     }

   if (0 != read_forecast_vars (mft, path))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error reading forecast variables from file: %s",
                     __func__, path);
        return -1;
     }

   return 0;
}

static int mft_interp_temp_on_isobars (Met_File_Type *mft, Met_Value_Type *mvt, int offset)
{
   LCC_Grid_Type *g = mft->grid;
   MET_Data_Type *vars = &mft->vars;
   Array_Type *isobars = &mft->isobars;
   Array_Type *var_on_isobar = &mft->var_on_isobar;
   size_t num_grid = g->x.n * g->y.n;
   size_t i, ni = isobars->n;
   const double *pi = isobars->v;
   double *ti = var_on_isobar->v;
   float *t = mvt->temperature_on_isobar;
   gsl_interp *interp = NULL;
   gsl_interp_accel *acc = NULL;
   gsl_error_handler_t *old_error_handler;
   int k, status = -1;

   for (i = 0; i < ni; i++)
     {
        const double *t_slab = vars->temperature_on_isobar + i * num_grid;
        ti[i] = t_slab[offset];
     }

   old_error_handler = gsl_set_error_handler_off();

   if ((NULL == (acc = gsl_interp_accel_alloc ()))
       || (NULL == (interp = gsl_interp_alloc (gsl_interp_linear, ni)))
       || (0 != gsl_interp_init (interp, pi, ti, ni)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initializing interpolation", __func__);
        goto return_status;
     }

   for (k = 0; k < mvt->num_isobars; k++)
     {
        double p = mvt->isobars[k];
        t[k] = gsl_interp_eval (interp, pi, ti, p, acc);
     }

   status = 0;
return_status:
   gsl_set_error_handler (old_error_handler);
   gsl_interp_free (interp);
   gsl_interp_accel_free (acc);

   return status;
}

static int mft_interp (Met_File_Type *mft, float lon_f, float lat_f, Met_Value_Type *mvt)
{
   LCC_Grid_Type *g = mft->grid;
   MET_Data_Type *vars = &mft->vars;
   double x = lon_f * DEGTORAD;
   double y = lat_f * DEGTORAD;
   int nx, ny, ix, iy, o;
   int status, num_pts = 1;

   if ((status = pj_transform (g->longlat, g->lcc, num_pts, 1, &x, &y, NULL)) != 0)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return MFT_INTERP_FAIL;
     }

   if ((x < g->xlim[0] || x > g->xlim[1])
       || (y < g->ylim[0] || y > g->ylim[1]))
     {
        return MFT_INTERP_DOMAIN_ERROR;
     }

   ix = 0.5 + (x - g->x.v[0]) / g->dx;
   iy = 0.5 + (y - g->y.v[0]) / g->dy;

   nx = g->x.n;
   ny = g->y.n;

   if (ix >= nx) ix = nx-1;
   else if (ix < 0) ix = 0;

   if (iy >= ny) iy = ny-1;
   else if (iy < 0) iy = 0;

   /* In grib2 files, the longitude index is the most rapidly varying */
   o = ix + iy * nx;

   if (vars->pressure_surface != NULL)
     {
        mvt->pressure_surface = vars->pressure_surface[o];
     }

   if (vars->pressure_tropopause != NULL)
     {
        mvt->pressure_tropopause = vars->pressure_tropopause[o];
     }

   if ((vars->temperature_on_isobar != NULL) && (mvt->isobars != NULL))
     {
        if (0 != mft_interp_temp_on_isobars (mft, mvt, o))
          return MFT_INTERP_FAIL;
     }

   return MFT_INTERP_SUCCESS;
}

Met_File_Type *met_open_file_grib2 (const char *path, int flags)
{
   Met_File_Type *mft = NULL;

   if (NULL == (mft = (Met_File_Type *)MALLOC (sizeof *mft)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)mft, 0, sizeof *mft);

   mft->flags = flags;
   mft->mft_close = mft_close;
   mft->mft_interp = mft_interp;

   if (0 != read_grib2_file (mft, path))
     {
        mft_close (mft);
        return NULL;
     }

   return mft;
}

int met_linear_interp (double *x0, double *y0, int n0,
                       int n, double *x, double *y)
{
   gsl_interp *interp = NULL;
   gsl_interp_accel *acc = NULL;
   gsl_error_handler_t *old_error_handler;
   int i, status = -1;

   old_error_handler = gsl_set_error_handler_off();

   if ((NULL == (acc = gsl_interp_accel_alloc ()))
       || (NULL == (interp = gsl_interp_alloc (gsl_interp_linear, n0)))
       || (0 != gsl_interp_init (interp, x0, y0, n0)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: initializing interpolation", __func__);
        goto return_status;
     }

   for (i = 0; i < n; i++)
     {
        y[i] = gsl_interp_eval (interp, x0, y0, x[i], acc);
     }

   status = 0;
return_status:
   gsl_set_error_handler (old_error_handler);
   gsl_interp_free (interp);
   gsl_interp_accel_free (acc);

   return status;
}
