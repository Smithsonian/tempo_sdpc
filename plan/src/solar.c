/** @file solar.c
 *  @brief Interface for solar illumination geometry
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include <libconfig.h>
#include <libnovas.h>

#include <tell.h>
#include <tio.h>

#ifndef SEC_PER_DAY
# define SEC_PER_DAY 86400.0
#endif

typedef struct
{
   double jd_utc;
   double jd_tt;
   double jd_ut1;
   double delta_t;
}
Times_Type;

#define SGT_PRIVATE_DATA \
   Times_Type *times; \
   Novas_object_t sun; \
   Novas_observer_t geocenter; \
   Novas_observer_t boresight_surface; \
   double sat_longitude; /* radians */ \
   double bs_longitude; /* radians */ \
   double bs_latitude; /* radians */ \
   double bs_elevation_angle; /* radians */ \
   double bs_azimuth_angle; /* radians */ \
   double sat_pos[3]; \
   double xpole; \
   double ypole; \
   double ut1_utc; \
   double unix_epoch_jd; \
   short int accuracy;
#include "solar.h"

#define KM_PER_AU      149597871.0
#define GEO_SAT_RADIUS     42163.968  /* km */
#define EARTH_MEAN_RADIUS   6371.0088 /* km */
#define DEGTORAD       (M_PI/180.0)

#define DEFAULT_HEIGHT        0.0   /* height [meters] */
#define DEFAULT_TEMPERATURE  15.0   /* temperature [Celsius] */
#define DEFAULT_PRESSURE   1010.0   /* pressure [millibars] */

static int times_eval (Times_Type *tt, double jd_utc,
                       int leap_secs, double ut1_utc)
{
#define TT_OFFSET_SECS  32.184

   tt->jd_utc  = jd_utc;
   tt->jd_tt   = jd_utc + (leap_secs + TT_OFFSET_SECS)/SEC_PER_DAY;
   tt->jd_ut1  = jd_utc + ut1_utc/SEC_PER_DAY;
   tt->delta_t = TT_OFFSET_SECS + leap_secs - ut1_utc;

   return 0;
}

static int leap_seconds (const Solar_Geom_Type *sgt, double jd_utc,
                         int *leap_secs)
{
   double utc_secs, tempo, tai_secs;

   utc_secs = (jd_utc - sgt->unix_epoch_jd) * SEC_PER_DAY;

   if ((0 != tio_time_utc_to_taix (utc_secs, &tempo))
       || (0 != tio_time_taix_to_tai (tempo, &tai_secs)))
     return -1;

   *leap_secs = tai_secs - utc_secs;

   return 0;
}

static int sgt_solar_zenith_angle (const Solar_Geom_Type *sgt,
                                   double jd_utc, double lon, double lat,
                                   double *psza)
{
   Times_Type tt;
   Novas_observer_t surf_obs;
   Novas_sky_pos_t sun_place;
   Novas_object_t sun = sgt->sun;  /* struct copy */
   short int error;
   short int coord_sys = 0;     /* 0 means GCRS coordinates */
   short int refr_option = 0;   /* 0 means no refraction */
   double rar, decr;            /* ra,dec adjusted for refraction */
   double zd, az;
   int leap_secs;

   if (0 != leap_seconds (sgt, jd_utc, &leap_secs))
     return -1;

   if (0 != times_eval (&tt, jd_utc, leap_secs, sgt->ut1_utc))
     return -1;

   novas_make_observer_on_surface (lat, lon, DEFAULT_HEIGHT,
                                   DEFAULT_TEMPERATURE, DEFAULT_PRESSURE,
                                   &surf_obs);

   if ((error = novas_place (tt.jd_tt, &sun, &surf_obs, tt.delta_t,
                             coord_sys, sgt->accuracy, &sun_place)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_place",
                     __func__, error);
        return -1;
     }

   novas_equ2hor (tt.jd_ut1, tt.delta_t, sgt->accuracy,
                  sgt->xpole, sgt->ypole, &surf_obs.on_surf,
                  sun_place.ra, sun_place.dec, refr_option,
                  &zd, &az, &rar, &decr);

   *psza = zd;
   return 0;
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

static void vec_cross (const double *a, const double *b, double *c)
{
   c[0] =   a[1] * b[2] - a[2] * b[1];
   c[1] = - a[0] * b[2] + a[2] * b[0];
   c[2] =   a[0] * b[1] - a[1] * b[0];
}

static void vec_scale (double c, const double *a, double *b)
{
   b[0] = a[0] * c;
   b[1] = a[1] * c;
   b[2] = a[2] * c;
}

/* Rotate 3D vector v about axis k by an angle theta, yielding vector vrot,
 * where theta > 0 is defined according to the right hand rule.
 * Rodrigues' rotation formula:
 *   vrot = v cos(t) + (k \cross v) sin(t) + k (k \dot v) (1 - cos(t))
 */
static void vec_rotate (const double *v, const double *k,
                        double theta, double *vrot)
{
   double cross[3], n[3], p[3];
   double dot, cos_t = cos(theta);

   vec_cross (k, v, cross);
   vec_scale (sin(theta), cross, n);

   dot = vec_dot(k, v);
   vec_scale (dot * (1 - cos_t), k, p);

   vec_scale (cos_t, v, vrot);
   vrot[0] += n[0] + p[0];
   vrot[1] += n[1] + p[1];
   vrot[2] += n[2] + p[2];
}

static double vec_angle (double *pa, double *pb)
{
   double a[3], b[3];
   vec_norm (pa, a);
   vec_norm (pb, b);
   return acos (vec_dot(a, b));
}

static int sgt_sat_sun_angles (Solar_Geom_Type *sgt, double jd_utc, double *ptheta, double *pphi)
{
   Times_Type tt;
   Novas_sky_pos_t sun_place;
   double sat_gcrs[3];
   double bs_gcrs[3], bs_gcrs_vel[3];
   double bs_sat[3], sun_sat[3];
   double tilt_angle = sgt->bs_elevation_angle; /* radians */
   double r_sun;
   short int error;
   short int coord_sys = 0;   /* 0 means GCRS coordinates */
   int leap_secs;
   int method = 1; /* 1 = equinox-based method */
   int option = 0; /* 0 = output vector referred to GCRS axes */
   int i;

   if (0 != leap_seconds (sgt, jd_utc, &leap_secs))
     return -1;

   if (0 != times_eval (&tt, jd_utc, leap_secs, sgt->ut1_utc))
     return -1;

   /* convert sat position from ITRS to GCRS system
    * returns sat_gcrs [km] */
   if ((error = novas_ter2cel (tt.jd_ut1, 0.0, tt.delta_t,
                               method, sgt->accuracy, option,
                               sgt->xpole, sgt->ypole, sgt->sat_pos,
                               sat_gcrs)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_ter2cel",
                     __func__, error);
        return -1;
     }

   /* returns bs_gcrs [AU], bs_gcrs_vel [AU/day] */
   if ((error = novas_geo_posvel (tt.jd_tt, tt.delta_t, sgt->accuracy,
                                  &sgt->boresight_surface, bs_gcrs, bs_gcrs_vel)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_geo_posvel",
                     __func__, error);
        return -1;
     }
   for (i = 0; i < 3; i++)
     {
        bs_gcrs[i] *= KM_PER_AU;
     }

   /* returns sun_place.dis [AU] */
   if ((error = novas_place (tt.jd_tt, &sgt->sun, &sgt->geocenter,
                             tt.delta_t, coord_sys, sgt->accuracy,
                             &sun_place)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_place",
                     __func__, error);
        return -1;
     }

   r_sun = sun_place.dis * KM_PER_AU;
   for (i = 0; i < 3; i++)
     {
        double sun_gcrs_i = r_sun * sun_place.r_hat[i];
        double sat_gcrs_i = sat_gcrs[i];
        bs_sat[i] = bs_gcrs[i] - sat_gcrs_i;
        sun_sat[i] = sun_gcrs_i - sat_gcrs_i;
     }

   *ptheta = vec_angle (bs_sat, sun_sat) / DEGTORAD;

   if (pphi)
     {
        double dot, hat_bs_sat[3], hat_sat[3], u[3], hat_u[3];
        double hat_z[3], hat_vel[3], hat_slit[3];
        double cos_phi;

        /* \H = host satellite position vector
         * \S = sun pos. vec.
         * \B = boresight pos. vec.
         * \z = earth rotation axis unit vector
         * \l = unit vector along slit
         * \v = unit tangent vector in the direction of the orbital velocity
         *   = \z x \h,  where \h = \H/norm(\H)
         * R(\v,\k,theta) = rotation of vector \v about axis \k by angle theta,
         *                  where theta>0 is defined according to the right-hand rule.
         *
         * define: \P = \S - \H = vector toward the sun from the satellite
         * define: \n = (\B-\H) / norm(\B-\H) = unit normal to diffuser plate
         * define: \U = \P - (\P dot \n)\n
         *            = vector component of \P perpendicular to the boresight
         * define: \u = \U / norm(\U)
         *
         * \l = R(\z, \v, tilt) = z unit vector rotated about an axis
         *                        along the orbital velocity vector.
         *
         * cos(phi) = \u dot \l
         */

        vec_norm (bs_sat, hat_bs_sat);                      /* \n */
        dot = vec_dot (sun_sat, hat_bs_sat);                /* \P dot \n */
        u[0] = sun_sat[0] - dot * hat_bs_sat[0];
        u[1] = sun_sat[1] - dot * hat_bs_sat[1];
        u[2] = sun_sat[2] - dot * hat_bs_sat[2];            /* \U */
        vec_norm (u, hat_u);                                /* \u */
        hat_z[0] = 0.0;
        hat_z[1] = 0.0;
        hat_z[2] = 1.0;                                     /* \z */
        vec_norm (sat_gcrs, hat_sat);                       /* \h */
        vec_cross (hat_z, hat_sat, hat_vel);                /* \v */
        vec_rotate (hat_z, hat_vel, tilt_angle, hat_slit);  /* \l */
        cos_phi = vec_dot (hat_u, hat_slit);                /* \u dot \l */
        *pphi = acos(cos_phi) / DEGTORAD;
     }

   return 0;
}

static int sgt_geosat_longitude (const Solar_Geom_Type *sgt, double *lon)
{
   if (lon) *lon = sgt->sat_longitude;
   return 0;
}

static int sgt_boresight_angles (const Solar_Geom_Type *sgt, double *elev_deg, double *azi_deg)
{
   if (elev_deg) *elev_deg = sgt->bs_elevation_angle / DEGTORAD;
   if (azi_deg) *azi_deg = sgt->bs_azimuth_angle / DEGTORAD;
   return 0;
}

static int sgt_print_params (const Solar_Geom_Type *sgt, const char *pprefix, FILE *fp)
{
   const char *prefix = pprefix ? pprefix : "";
   (void) fprintf (fp, "%s (%8.3f,%7.3f) = instrument boresight point\n",
                   prefix,
                   sgt->boresight_surface.on_surf.longitude /DEGTORAD,
                   sgt->boresight_surface.on_surf.latitude /DEGTORAD);
   (void) fprintf (fp, "%s (%8.3f,%7.3f) = satellite nadir point \n",
                   prefix, sgt->sat_longitude /DEGTORAD, 0.0);
   return 0;
}

static int read_sat_config (config_t *cfg, double *sat_lon,
                            double *bs_lon, double *bs_lat)
{
   config_setting_t *s, *sub;

   if (NULL == (s = config_lookup (cfg, "sat_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing sat_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "sat_lon", sat_lon))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading sat_lon: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (sub = config_setting_get_member (s, "boresight")))
     return -1;

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "lon", bs_lon))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading boresight longitude: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (sub, "lat", bs_lat))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading boresight latitude: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   *sat_lon *= DEGTORAD;
   *bs_lon *= DEGTORAD;
   *bs_lat *= DEGTORAD;

   return 0;
}

static int read_iers_params (Solar_Geom_Type *sgt, config_t *cfg)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "iers_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing iers_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "xpole", &sgt->xpole))
       ||(CONFIG_TRUE != config_setting_lookup_float (s, "ypole", &sgt->ypole)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading iers_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "ut1_utc", &sgt->ut1_utc))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading iers_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int sgt_initialize (Solar_Geom_Type *sgt, config_t *cfg)
{
   Novas_cat_entry_t dummy_star;
   double a = GEO_SAT_RADIUS / EARTH_MEAN_RADIUS;
   double tan_elev;
   short int error;

   /* We'll need the positions of:
    *  - the sun
    *  - observer at the geocenter
    *  - geostationary satellite
    *  - observer on the surface, at the satellite's boresight coordinates
    */

   if ((error = novas_make_object (0, 10, "sun", &dummy_star, &sgt->sun)) != 0)
   {
      tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from make_object (Sun)",
                   __func__, error);
      return -1;
   }

   novas_make_observer_at_geocenter (&sgt->geocenter);

   if (0 != read_sat_config (cfg, &sgt->sat_longitude, &sgt->bs_longitude, &sgt->bs_latitude))
     return -1;

   /* Elevation and azimuth angles to point boresight at specified (lon,lon). */
   if (fabs (sgt->sat_longitude - sgt->bs_longitude) > DBL_EPSILON)
     {
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR,
                     "%s: The code assumes that the satellite and the boresight are at the same longitude", __func__);
        return -1;
     }

   tan_elev = sin(sgt->bs_latitude) / (a - cos(sgt->bs_latitude));
   sgt->bs_elevation_angle = atan(tan_elev); /* radians */
   sgt->bs_azimuth_angle = 0.0;

   /* WGS84 coordinates of geostationary satellite */
   sgt->sat_pos[0] = GEO_SAT_RADIUS * cos(sgt->sat_longitude);  /* X */
   sgt->sat_pos[1] = GEO_SAT_RADIUS * sin(sgt->sat_longitude);  /* Y */
   sgt->sat_pos[2] = 0.0;                            /* Z */

   novas_make_observer_on_surface (sgt->bs_latitude, sgt->bs_longitude, DEFAULT_HEIGHT,
                                   DEFAULT_TEMPERATURE, DEFAULT_PRESSURE,
                                   &sgt->boresight_surface);

   /* NOVAS accuracy parameter (0=full) */
   sgt->accuracy = 0;

   return read_iers_params (sgt, cfg);
}

static void sgt_delete (Solar_Geom_Type *sgt)
{
   if (NULL == sgt)
     return;
   FREE(sgt);
}

Solar_Geom_Type *solar_geom_init (config_t *cfg)
{
   Solar_Geom_Type *sgt = NULL;

   if (NULL == (sgt = (Solar_Geom_Type *)MALLOC (sizeof *sgt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sgt, 0, sizeof *sgt);

   sgt->unix_epoch_jd = novas_julian_date (1970,1,1,0.0);

   sgt->sgt_delete = sgt_delete;
   sgt->sgt_solar_zenith_angle = sgt_solar_zenith_angle;
   sgt->sgt_sat_sun_angles = sgt_sat_sun_angles;
   sgt->sgt_geosat_longitude = sgt_geosat_longitude;
   sgt->sgt_boresight_angles = sgt_boresight_angles;
   sgt->sgt_print_params = sgt_print_params;

   if (0 != sgt_initialize (sgt, cfg))
     {
        sgt_delete (sgt);
        return NULL;
     }

   return sgt;
}
