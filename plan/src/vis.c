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

#include <tell.h>
#include <tio.h>

#include "scan.h"
#include "scan_methods.h"
#include "solar.h"
#include "plan_list.h"
#include "vis.h"

#define DEGTORAD                (M_PI/180.0)
#define PROJ_ARGS_BUFSIZE       80

#define GEO_ALTITUDE  35785831.0   /* meters */

#define MICRORADIAN (1.e-6)

struct Vis_Type
{
   Solar_Geom_Type *solar_geom;
   const char *plan_id;
   double sat_lon;

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

typedef struct Point_Type Point_Type;
struct Point_Type
{
   double lon;
   double lat;
};

typedef struct Box_Type Box_Type;
struct Box_Type
{
   Point_Type min;
   Point_Type max;
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

/* Map mirror angles (x,y) to (lon,lat) */

/* lon, lat are of size (2*nx + 2*ny) */
static int vis_scan_box_to_lonlat (const Vis_Type *v, const Plan_List_Type *entry,
                                   double step_size, int nx, int ny,
                                   double *lon, double *lat)
{
   double xmax = entry->xstart + entry->num_steps * step_size;
   double xmin = entry->xstart;
   double ymax = +0.5 * v->ysize / MICRORADIAN;
   double ymin = -0.5 * v->ysize / MICRORADIAN;
   double dx, dy;
   double *x = NULL;
   double *y = NULL;
   int status = -1;
   int i, k, len;

   len = 2*nx + 2*ny;

   if (NULL == (x = (double *)MALLOC (2 * len * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   y = x + len;

   dx = (xmax - xmin) / (nx-1);
   dy = (ymax - ymin) / (ny-1);

   k = 0;
   for (i = 0; i < nx; i++, k++)
     {
        x[k] = xmin + i * dx;
        y[k] = ymin;
     }
   for (i = 0; i < ny; i++, k++)
     {
        x[k] = xmax;
        y[k] = ymin + i * dy;
     }
   for (i = nx-1; i >= 0; i--, k++)
     {
        x[k] = xmin + i * dx;
        y[k] = ymax;
     }
   for (i = ny-1; i >= 0; i--, k++)
     {
        x[k] = xmin;
        y[k] = ymin + i * dy;
     }

   if (0 != scan_xy_to_lonlat (x, y, len, lon, lat, v->sat_lon))
     goto return_status;

   status = 0;
return_status:
   FREE(x);
   return status;
}

static int vis_azel_to_lonlat (const Vis_Type *v, double az_urad, double el_urad,
                               double *plon, double *plat)
{
   return scan_xy_to_lonlat (&az_urad, &el_urad, 1, plon, plat, v->sat_lon);
}

static int vis_define_mesh (Vis_Type *v, const Box_Type *bbox, double center_lon)
{
   double x0, y0, dx, dy, cos_phi1;
   int i, j;

   /* Plate Carree projection is an Equidistant Cylindrical projection
    * with the equator as standard parallel, so cos(phi1) = 1 */
   cos_phi1 = 1.0;

   /* rectangular grid for Equidistant Cylindrical projection */
   dx = (bbox->max.lon - bbox->min.lon) * cos_phi1 / (v->num_lon-1);
   dy = (bbox->max.lat - bbox->min.lat) / (v->num_lat - 1);

   x0 = (bbox->min.lon - center_lon) * cos_phi1;
   y0 = bbox->min.lat;

   for (j = 0; j < v->num_lat; j++)
     {
        for (i = 0; i < v->num_lon; i++)
          {
             int k = i + j * v->num_lon;
             /* Equidistant Cylindrical projection */
             v->x[k] = x0 + i * dx;
             v->y[k] = y0 + j * dy;
             /* longitude, latitude [deg] */
             v->lat[k] = v->y[k];
             v->lon[k] = center_lon + v->x[k] / cos_phi1;
          }
     }

   return 0;
}

static int read_vis_params (Vis_Type *v, config_t *cfg, int *img_size, Box_Type *box, double *center_lon)
{
   config_setting_t *s;
   config_setting_t *sub;
   int num_xtrack;
   double size_xtrack_pixel, yoffset;

   if (NULL == (s = config_lookup (cfg, "output_sza_map_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing output_sza_map_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "img_size", img_size))
       || (CONFIG_TRUE != config_setting_lookup_int (s, "num_xtrack", &num_xtrack))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "size_xtrack_pixel", &size_xtrack_pixel))
       )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading output_sza_map_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (sub = config_setting_get_member (s, "bounding_box")))
       ||(CONFIG_TRUE != config_setting_lookup_float (sub, "lon_min", &box->min.lon))
       ||(CONFIG_TRUE != config_setting_lookup_float (sub, "lat_min", &box->min.lat))
       ||(CONFIG_TRUE != config_setting_lookup_float (sub, "lon_max", &box->max.lon))
       ||(CONFIG_TRUE != config_setting_lookup_float (sub, "lat_max", &box->max.lat)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading output_sza_map_config.bounding_box: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "yoffset", &yoffset))
     yoffset = 0.0;

   scan_set_lonlat_bounding_box (box->min.lon, box->max.lon,
                                 box->min.lat, box->max.lat, yoffset);

   if (CONFIG_TRUE != config_setting_lookup_float (s, "center_lon", center_lon))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading output_sza_map_config.center_lon: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   /* convert microradians to radians */
   size_xtrack_pixel *= 1.e-6;

   /* FOR size [radians] */
   v->ysize = num_xtrack * size_xtrack_pixel * cos (v->tilt);

   return 0;
}

Vis_Type *vis_init (config_t *cfg, Solar_Geom_Type *solar_geom, const char *plan_id)
{
   Vis_Type *v = NULL;
   Box_Type box = {0};
   double center_lon;
   double sat_lon, unused_azimuth_angle_about_z_axis;
   int img_size;

   if (NULL == (v = (Vis_Type *) MALLOC (sizeof *v)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)v, 0, sizeof (*v));

   v->solar_geom = solar_geom;
   v->plan_id = plan_id;

   if (0 != read_vis_params (v, cfg, &img_size, &box, &center_lon))
     goto return_error;

   /* tilt>0 [deg] is northward tilt of instrument boresight,
    *              about spacecraft roll axis
    * azi>0 [deg] is rotation, eastward from north (CW),
    *             about instrument boresight axis
    * These angles specify the instrument pointing direction
    * for the tilted perspective projection ('tpers') from
    * the Proj4 library.
    *
    * Note that if the sgt_boresight_angles returned an azimuth angle,
    * it would correspond to a rotation about an axis parallel to the
    * Earth's spin axis, which is NOT the same as the azi angle that Proj4 wants.
    * So, if this gets generalized to support the case where the boresight
    * differs from the satellite longitude, fixing these angles will take
    * a bit of care.
    *
    * FIXME!!  Add support for boresight not at satellite longitude
    */
   if (0 != solar_geom->sgt_boresight_angles (solar_geom, &v->tilt, &unused_azimuth_angle_about_z_axis))
     goto return_error;

   v->tilt *= DEGTORAD;
   v->azi = 0.0;

   /* FOR center coordinates [radians] */
   v->x0 = 0.0;
   v->y0 = v->tilt;

   if (0 != solar_geom->sgt_geosat_longitude (solar_geom, &sat_lon))
     goto return_error;
   v->sat_lon = sat_lon;

   if (0 != vis_alloc_grid (v, img_size, img_size))
     goto return_error;

   /* Define a regular grid in the Plate Carree projection,
    * which is convenient for plotting (more convenient than (lon,lat)),
    * yet simply related to (lon,lat)
    * The output file will have both (lon,lat) and (x,y) = Plate Carree
    */
   if (0 != vis_define_mesh (v, &box, center_lon))
     goto return_error;

   return v;

return_error:
   vis_free (v);
   return NULL;
}

static void vec_norm (const double *a, double *norm)
{
   double r = sqrt (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
   norm[0] = a[0]/r;
   norm[1] = a[1]/r;
   norm[2] = a[2]/r;
}

static double vec_dot (const double *a, const double *b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void vec_unit (double lon_deg, double lat_deg, double *p)
{
   double theta = (90.0 - lat_deg) * DEGTORAD;
   double phi = lon_deg * DEGTORAD;
   double sin_t = sin(theta);
   p[0] = sin_t * cos(phi);
   p[1] = sin_t * sin(phi);
   p[2] = cos(theta);
}

double *vis_sza (const Vis_Type *v, double jd_utc, double *psza)
{
   int compare_with_actual_sza = 0; /* non-zero turns on slow comparison */
   double sun_itrs[3], p_sun[3], p[3];
   double max_diff, max_pos[2];
   double *sza = NULL;
   double *lon = v->lon;
   double *lat = v->lat;
   Solar_Geom_Type *sgt = v->solar_geom;
   int i, n = v->num_lon * v->num_lat;

   if (0 != sgt->sgt_solar_xyz (sgt, jd_utc, sun_itrs))
     return NULL;
   vec_norm (sun_itrs, p_sun);

   if (psza == NULL)
     {
        if (NULL == (sza = (double *)MALLOC (n * sizeof(double))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return NULL;
          }
     }
   else sza = psza;

   /* Because the Earth radius is small compared to the Earth-Sun
    * distance, the SZA at a surface point, P, is about the same as
    * the geocenter angle (sun)-(geocenter)-(P).
    * The angle at the geocenter is faster to compute, and the
    * accuracy should be good enough for visualization.
    */

   max_diff = 0.0;
   max_pos[0] = 0.0;
   max_pos[1] = 0.0;

   for (i = 0; i < n; i++)
     {
        if ((lat[i] == TIO_FILL_FLOAT) || (lon[i] == TIO_FILL_FLOAT))
          {
             sza[i] = TIO_FILL_FLOAT;
             continue;
          }

        /* Approximate the SZA using the angle at the geocenter. */
        vec_unit (lon[i], lat[i], p);
        sza[i] = acos(vec_dot (p_sun, p)) / DEGTORAD;

        if (compare_with_actual_sza)
          {
             double sza_novas, diff;
             if (0 != sgt->sgt_solar_zenith_angle (sgt, jd_utc, lon[i], lat[i], &sza_novas))
               return NULL;
             diff = sza_novas - sza[i];
             if (fabs(diff) > max_diff)
               {
                  max_diff = diff;
                  max_pos[0] = lon[i];
                  max_pos[1] = lat[i];
               }
             /* use NOVAS value */
             sza[i] = sza_novas;
          }
     }

   if (compare_with_actual_sza)
     {
        fprintf (stderr, "visualization: max SZA error = %10.6f deg at lon=%10.6f lat=%10.6f\n",
                 max_diff, max_pos[0], max_pos[1]);
     }

   return sza;
}

int vis_write_grid (Vis_Type *v, int ncid, double *control_points)
{
   const char name_lon[] = "longitude";
   const char name_lat[] = "latitude";
   float float_missing = TIO_FILL_FLOAT;
   int varid_lon, varid_lat, varid_x, varid_y;
   int dimid_lon, dimid_lat;
   int start[2], count[2];
   int i, status = -1;

   if (v->plan_id != NULL)
     {
        int len = strlen(v->plan_id) + 1;
        if (0 != TIO_put_att (ncid, NC_GLOBAL, "plan_id", NC_CHAR, len, v->plan_id))
          goto return_status;
     }

   if ((0 != TIO_put_att (ncid, NC_GLOBAL, "day_begin_ctrl_point", NC_DOUBLE, 2, control_points))
       || (0 != TIO_put_att (ncid, NC_GLOBAL, "day_end_ctrl_point", NC_DOUBLE, 2, control_points+2)))
     goto return_status;

   if ((0 != TIO_def_dim (ncid, name_lon, v->num_lon, &dimid_lon))
       ||(0 != TIO_def_dim (ncid, name_lat, v->num_lat, &dimid_lat)))
     goto return_status;

   if ((0 != TIO_def_var (ncid, name_lon, NC_FLOAT, 1, &dimid_lon, &varid_lon))
       || (0 != TIO_def_var (ncid, name_lat, NC_FLOAT, 1, &dimid_lat, &varid_lat)))
     goto return_status;

   start[0] = 0;

   if (0 != TIO_put_var_section (ncid, name_lon, start, &v->num_lon, NC_DOUBLE, v->lon))
     goto return_status;

   start[1] = 0;
   count[0] = 1;
   count[1] = 1;
   for (i = 0; i < v->num_lat; i++)
     {
        double *v_lat = v->lat + i * v->num_lon;
        start[0] = i;
        if (0 != TIO_put_var_section (ncid, name_lat, start, count, NC_DOUBLE, v_lat))
          goto return_status;
     }

   v->dimids_lon_lat[0] = dimid_lat;
   v->dimids_lon_lat[1] = dimid_lon;
   if ((0 != TIO_def_var (ncid, "x", NC_FLOAT, 2, v->dimids_lon_lat, &varid_x))
       || (0 != TIO_def_var (ncid, "y", NC_FLOAT, 2, v->dimids_lon_lat, &varid_y)))
     goto return_status;

   if ((0 != TIO_put_att (ncid, varid_x, "missing_value", NC_FLOAT, 1, &float_missing))
       || (0 != TIO_def_var_fill (ncid, varid_x, 0, &float_missing)))
     goto return_status;
   if ((0 != TIO_put_att (ncid, varid_y, "missing_value", NC_FLOAT, 1, &float_missing))
       || (0 != TIO_def_var_fill (ncid, varid_y, 0, &float_missing)))
     goto return_status;

   start[0] = 0;
   start[1] = 0;
   count[0] = v->num_lat;
   count[1] = v->num_lon;

   if ((0 != TIO_put_var_section (ncid, "x", start, count, NC_DOUBLE, v->x))
       ||(0 != TIO_put_var_section (ncid, "y", start, count, NC_DOUBLE, v->y)))
     goto return_status;

   status = 0;
return_status:
   return status;
}

#define NBOX_LON 100
#define NBOX_LAT 100
#define NBOX (2*NBOX_LON + 2*NBOX_LAT)

int vis_write_value (const Vis_Type *v, int ncid, double jd_utc,
                     const char *name, const double *value,
                     double step_size, const Plan_List_Type *entry)
{
   const char coord_attr[] = "longitude latitude";
   int varid, start[2], count[2];
   char buf[32];
   float float_missing = TIO_FILL_FLOAT;
   double pos[2], scan_angle, scan_duration, solar_boresight_angle;
   double box_lon[NBOX], box_lat[NBOX];
   double box_lon_filtered[NBOX], box_lat_filtered[NBOX];
   int nx = NBOX_LON;
   int ny = NBOX_LAT;
   int n = NBOX;
   int num_scans, i, k;
   int status = -1;

   if ((0 != TIO_def_var (ncid, name, NC_FLOAT, 2, v->dimids_lon_lat, &varid))
       || (0 != TIO_put_att (ncid, varid, "coordinates", NC_CHAR, strlen(coord_attr), &coord_attr)))
     goto return_status;

   if (0 != mkjdtimestr (jd_utc, buf, sizeof(buf)))
     goto return_status;

   if (0 != v->solar_geom->sgt_sat_sun_position (v->solar_geom, jd_utc, &solar_boresight_angle, NULL, NULL))
     goto return_status;

   if (0 != vis_azel_to_lonlat (v, entry->xstart, entry->ystart, &pos[0], &pos[1]))
     goto return_status;

   scan_angle = entry->num_steps * step_size * 1.e-6;  /* [rad] */

   num_scans = (entry->num_repeats_cbm > 0) ? entry->num_repeats_cbm : 1;
   scan_duration = entry->scan_duration / num_scans;

   if ((0 != TIO_put_att (ncid, varid, "julian_date", NC_DOUBLE, 1, &jd_utc))
       ||(0 != TIO_put_att (ncid, varid, "julian_date_str", NC_CHAR, strlen(buf)+1, buf))
       ||(0 != TIO_put_att (ncid, varid, "solar_boresight_angle", NC_DOUBLE, 1, &solar_boresight_angle))
       ||(0 != TIO_put_att (ncid, varid, "scan_duration", NC_DOUBLE, 1, &scan_duration))
       ||(0 != TIO_put_att (ncid, varid, "num_repeats", NC_INT, 1, &entry->num_repeats))
       ||(0 != TIO_put_att (ncid, varid, "num_repeats_cbm", NC_INT, 1, &entry->num_repeats_cbm))
       ||(0 != TIO_put_att (ncid, varid, "maneuver_loss", NC_DOUBLE, 1, &entry->maneuver_loss))
       ||(0 != TIO_put_att (ncid, varid, "start_pos", NC_DOUBLE, 2, pos))
       ||(0 != TIO_put_att (ncid, varid, "scan_angle_rad", NC_DOUBLE, 1, &scan_angle)))
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

   if (0 != vis_scan_box_to_lonlat (v, entry, step_size, nx, ny,
                                    box_lon, box_lat))
     goto return_status;

   /* Filter NaNs from box boundary. */
   k = 0;
   for (i = 0; i < n; i++)
     {
        if ((0 == isnan(box_lon[i])) && (0 == isnan(box_lat[i])))
          {
             box_lon_filtered[k] = box_lon[i];
             box_lat_filtered[k] = box_lat[i];
             k++;
          }
     }
   if (k < n-1)
     {
        /* close polygon */
        box_lon_filtered[k] = box_lon_filtered[0];
        box_lat_filtered[k] = box_lat_filtered[0];
        k++;
     }

   if ((0 != TIO_put_att (ncid, varid, "box_lon", NC_DOUBLE, k, box_lon_filtered))
       || (0 != TIO_put_att (ncid, varid, "box_lat", NC_DOUBLE, k, box_lat_filtered)))
     goto return_status;

   status = 0;
return_status:
   return status;
}
