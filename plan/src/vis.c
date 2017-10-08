/** @file vis.c
 *  @brief Map solar zenith angle vs position to help visualize plan
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <proj_api.h>

#include <tell.h>
#include <tio.h>

#include "scan.h"
#include "solar.h"
#include "vis.h"

#define DEGTORAD                (M_PI/180.0)
#define PROJ_ARGS_BUFSIZE       80

#define GEO_ALTITUDE  35785831.0   /* meters */

struct Vis_Type
{
   Solar_Geom_Type *solar_geom;

   double tilt;
   double azi;
   double xsize;
   double ysize;
   double x0, y0;

   double *x;
   double *y;
   int nx;
   int ny;

   double *lon;
   double *lat;
   int num_lon;
   int num_lat;

   int dimids_lon_lat[2];
};

void vis_free (Vis_Type *v)
{
   if (v == NULL)
     return;
   FREE(v->x);
   FREE(v);
}

static int vis_alloc_grid (Vis_Type *v, int nx, int ny)
{
   size_t len = nx * ny;

   if (NULL == (v->x = (double *)MALLOC (4 * len * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   v->nx = nx;
   v->ny = ny;
   v->num_lon = nx;
   v->num_lat = ny;

   v->y = v->x + len;
   v->lon = v->y + len;
   v->lat = v->lon + len;

   return 0;
}

static int vis_init_xygrid (Vis_Type *v)
{
   double *v_x = v->x;
   double *v_y = v->y;
   double xmin, xmax, dx;
   double ymin, ymax, dy;
   double x0 = v->x0;
   double y0 = v->y0;
   double cos_a = cos(-v->azi);
   double sin_a = sin(-v->azi);
   int ix, nx = v->nx;
   int iy, ny = v->ny;
   int k, nk = nx * ny;

   /* generate a uniform grid centered on (0,0) */

   xmin = -0.5 * v->xsize;
   xmax = +0.5 * v->xsize;
   ymin = -0.5 * v->ysize;
   ymax = +0.5 * v->ysize;

   dx = (xmax - xmin) / (nx-1);
   dy = (ymax - ymin) / (ny-1);

   for (iy = 0; iy < ny; iy++)
     {
        for (ix = 0; ix < nx; ix++)
          {
             k = ix + iy * nx;
             v_x[k] = xmin + ix * dx;
             v_y[k] = ymin + iy * dy;
          }
     }

   /* rotate grid, shift center to (x0,y0)
    */

   for (k = 0; k < nk; k++)
     {
        double x_k = x0 + v_x[k] * cos_a + v_y[k] * sin_a;
        double y_k = y0 - v_x[k] * sin_a + v_y[k] * cos_a;
        v_x[k] = x_k;
        v_y[k] = y_k;
     }

   return 0;
}

static int vis_xy_to_lonlat (Vis_Type *v, double sat_lon)
{
   projPJ tpers = NULL;
   projPJ longlat = NULL;
   char ctl_tpers[PROJ_ARGS_BUFSIZE];
   const char tpers_fmt[] =
     "+proj=tpers +lat_0=0 +lon_0=%0.3g +h=%0.1f +tilt=%0.3g +azi=%0.3g";
   const char ctl_longlat[] =
     "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
   double *v_x = v->x, *v_y = v->y;
   double *v_lon = v->lon, *v_lat = v->lat;
   int len, status = -1;
   long i, n;

   if (NULL == (longlat = pj_init_plus (ctl_longlat)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus(longlat) failed", __func__);
        goto return_status;
     }

   memset (ctl_tpers, 0, PROJ_ARGS_BUFSIZE);
   len = snprintf (ctl_tpers, PROJ_ARGS_BUFSIZE, tpers_fmt,
                   sat_lon/DEGTORAD, GEO_ALTITUDE,
                   v->tilt/DEGTORAD,
                   v->azi/DEGTORAD);
   if (len >= PROJ_ARGS_BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: proj4 arg buffer too small", __func__);
        goto return_status;
     }

   if (NULL == (tpers = pj_init_plus (ctl_tpers)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus(tpers) failed", __func__);
        goto return_status;
     }

   n = v->nx * v->ny;
   for (i = 0; i < n; i++)
     {
        v_lon[i] = v_x[i] * GEO_ALTITUDE;
        v_lat[i] = v_y[i] * GEO_ALTITUDE;
     }

   if ((status = pj_transform (tpers, longlat, n, 1, v_lon, v_lat, NULL)) != 0)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        goto return_status;
     }

   for (i = 0; i < n; i++)
     {
        if (isfinite (v_lon[i]))
          v_lon[i] /= DEGTORAD;
        else v_lon[i] = TIO_FILL_FLOAT;

        if (isfinite (v_lat[i]))
          v_lat[i] /= DEGTORAD;
        else v_lat[i] = TIO_FILL_FLOAT;
     }

   status = 0;
return_status:
   pj_free (longlat);
   pj_free (tpers);
   return status ? -1 : 0;
}

static int read_vis_params (Vis_Type *v, config_t *cfg, int *img_size)
{
   config_setting_t *s;
   int num_xtrack, num_mirror_steps;
   double size_mirror_step, size_xtrack_pixel;

   if (NULL == (s = config_lookup (cfg, "output_sza_map_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing output_sza_map_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "img_size", img_size))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num_xtrack", &num_xtrack))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num_mirror_steps", &num_mirror_steps))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "size_mirror_step", &size_mirror_step))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "size_xtrack_pixel", &size_xtrack_pixel))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "tilt", &v->tilt))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "azi", &v->azi))
       )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading output_sza_map_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   /* convert microradians to radians */
   size_mirror_step *= 1.e-6;
   size_xtrack_pixel *= 1.e-6;

   /* convert degrees to radians */
   v->tilt *= DEGTORAD;
   v->azi *= DEGTORAD;

   /* FOR size [radians] */
   v->xsize = num_mirror_steps * size_mirror_step;
   v->ysize = num_xtrack * size_xtrack_pixel * cos (v->tilt);

   /* FOR center coordinates [radians] */
   v->x0 = 0.0;
   v->y0 = v->tilt;

   return 0;
}

Vis_Type *vis_init (config_t *cfg, Solar_Geom_Type *solar_geom)
{
   Vis_Type *v = NULL;
   double sat_lon;
   int img_size;

   if (NULL == (v = (Vis_Type *) MALLOC (sizeof *v)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)v, 0, sizeof (*v));

   v->solar_geom = solar_geom;

   if (0 != read_vis_params (v, cfg, &img_size))
     goto return_error;

   if (0 != vis_alloc_grid (v, img_size, img_size))
     goto return_error;

   if (0 != vis_init_xygrid (v))
     goto return_error;

   if (0 != solar_geom->sgt_geosat_longitude (solar_geom, &sat_lon))
     goto return_error;

   if (0 != vis_xy_to_lonlat (v, sat_lon))
     goto return_error;

   return v;

return_error:
   vis_free (v);
   return NULL;
}

double *vis_sza (const Vis_Type *v, double jd_utc, double *psza)
{
   double *sza = NULL;
   double *lon = v->lon;
   double *lat = v->lat;
   Solar_Geom_Type *sgt = v->solar_geom;
   int i, n = v->num_lon * v->num_lat;

   if (psza == NULL)
     {
        if (NULL == (sza = (double *)MALLOC (n * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return NULL;
          }
     }
   else sza = psza;

   for (i = 0; i < n; i++)
     {
        if ((lat[i] == TIO_FILL_FLOAT) || (lon[i] == TIO_FILL_FLOAT))
          {
             sza[i] = TIO_FILL_FLOAT;
             continue;
          }
        if (0 != sgt->sgt_solar_zenith_angle (sgt, jd_utc, lon[i], lat[i], &sza[i]))
          {
             if (psza == NULL) FREE(sza);
             return NULL;
          }
     }

   return sza;
}

int vis_write_grid (Vis_Type *v, int ncid)
{
   const char name_lon[] = "longitude";
   const char name_lat[] = "latitude";
   float float_missing = TIO_FILL_FLOAT;
   int varid_lon, varid_lat;
   int dimid_lon, dimid_lat;
   int start[2], count[2];
   int status = -1;

   if ((0 != TIO_def_dim (ncid, name_lon, v->num_lon, &dimid_lon))
       ||(0 != TIO_def_dim (ncid, name_lat, v->num_lat, &dimid_lat)))
     goto return_status;

   v->dimids_lon_lat[0] = dimid_lat;
   v->dimids_lon_lat[1] = dimid_lon;
   if ((0 != TIO_def_var (ncid, name_lon, NC_FLOAT, 2, v->dimids_lon_lat, &varid_lon))
       || (0 != TIO_def_var (ncid, name_lat, NC_FLOAT, 2, v->dimids_lon_lat, &varid_lat)))
     goto return_status;

   if ((0 != TIO_put_att (ncid, varid_lon, "missing_value", NC_FLOAT, 1, &float_missing))
       || (0 != TIO_def_var_fill (ncid, varid_lon, 0, &float_missing)))
     goto return_status;
   if ((0 != TIO_put_att (ncid, varid_lat, "missing_value", NC_FLOAT, 1, &float_missing))
       || (0 != TIO_def_var_fill (ncid, varid_lat, 0, &float_missing)))
     goto return_status;

   start[0] = 0;
   start[1] = 0;
   count[0] = v->num_lat;
   count[1] = v->num_lon;

   if ((0 != TIO_put_var_section (ncid, name_lon, start, count, NC_DOUBLE, v->lon))
       ||(0 != TIO_put_var_section (ncid, name_lat, start, count, NC_DOUBLE, v->lat)))
     goto return_status;

   status = 0;
return_status:
   return status;
}

int vis_write_value (const Vis_Type *v, int ncid, double jd_utc,
                     const char *name, const double *value)
{
   int varid, start[2], count[2];
   char buf[32];
   float float_missing = TIO_FILL_FLOAT;
   int status = -1;

   if (0 != TIO_def_var (ncid, name, NC_FLOAT, 2, v->dimids_lon_lat, &varid))
     goto return_status;

   if (0 != mkjdtimestr (jd_utc, buf, sizeof(buf)))
     goto return_status;

   if ((0 != TIO_put_att (ncid, varid, "julian_date", NC_DOUBLE, 1, &jd_utc))
       ||(0 != TIO_put_att (ncid, varid, "julian_date_str", NC_CHAR, strlen(buf)+1, buf)))
     goto return_status;

   if ((0 != TIO_put_att (ncid, varid, "missing_value", NC_FLOAT, 1, &float_missing))
       || (0 != TIO_def_var_fill (ncid, varid, 0, &float_missing)))
     goto return_status;

   start[0] = 0;
   start[1] = 0;
   count[0] = v->num_lat;
   count[1] = v->num_lon;

   if (0 != TIO_put_var_section (ncid, name, start, count, NC_DOUBLE, value))
     goto return_status;

   status = 0;
return_status:
   return status;
}
