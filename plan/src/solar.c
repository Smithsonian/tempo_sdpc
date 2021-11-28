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
   double sat_longitude; /* radians */ \
   double bs_longitude; /* radians */ \
   double bs_latitude; /* radians */ \
   double bs_elevation_angle; /* radians */ \
   double bs_azimuth_angle; /* radians */ \
   double sat_pos[3]; \
   double bs_pos[3]; \
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

static int Compute_Boresight_Sun_Angle_In_ITRS_Frame = 1;
/* The angle between the boresight direction and the sun can be computed
 * in either the ITRS or GCRS frames, and the answer should be essentially
 * the same.  The ITRS frame calculation is simplest because two of
 * three relevant points (sat position and boresight point) are already
 * known in that frame.  The ITRS frame is used by default because it's
 * simplest.  The code for the GCRS calculation is available for comparison.
 */

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

static int sgt_solar_xyz (Solar_Geom_Type *sgt, double jd_utc, double sun_itrs[3])
{
   Times_Type tt;
   Novas_sky_pos_t sun_place;
   Novas_object_t sun = sgt->sun;  /* struct copy */
   short int error;
   short int coord_sys = 0;     /* 0 means GCRS coordinates */
   double r_sun, sun_gcrs[3];
   int method = 1; /* 1 = equinox-based method */
   int option = 0; /* 0 = output vector referred to GCRS axes */
   int i, leap_secs;

   if (0 != leap_seconds (sgt, jd_utc, &leap_secs))
     return -1;

   if (0 != times_eval (&tt, jd_utc, leap_secs, sgt->ut1_utc))
     return -1;

   if ((error = novas_place (tt.jd_tt, &sun, &sgt->geocenter, tt.delta_t,
                             coord_sys, sgt->accuracy, &sun_place)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_place",
                     __func__, error);
        return -1;
     }

   r_sun = sun_place.dis * KM_PER_AU;
   for (i = 0; i < 3; i++)
     {
        sun_gcrs[i] = r_sun * sun_place.r_hat[i];
     }

   /* convert sun position from GCRS to ITRS system
    * returns sun_itrs [km] */
   if ((error = novas_cel2ter (tt.jd_ut1, 0.0, tt.delta_t,
                               method, sgt->accuracy, option,
                               sgt->xpole, sgt->ypole, sun_gcrs,
                               sun_itrs)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_cel2ter",
                     __func__, error);
        return -1;
     }

   return 0;
}

static int sgt_sat_sun_position (Solar_Geom_Type *sgt, double jd_utc, double *ptheta, double *pphi,
                                 double *earth_sun_distance)
{
   Times_Type tt;
   Novas_sky_pos_t sun_place;
   double sun_gcrs[3], sun_itrs[3];
   double bs_gcrs[3], sat_gcrs[3];
   double sat_pos[3];
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

   if (Compute_Boresight_Sun_Angle_In_ITRS_Frame)
     {
        /* The simplest way to compute the solar-boresight angle is to compute
         * the solar position in GCRS and convert that to ITRS which is equivalent to
         * ECEF coordinates (e.g. WGS84, earth-centered, earth-fixed).
         * The other two points of interest are already known in the ECEF coordinates.
         */
        /* returns GCRS sun_place.dis [AU] */
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
             sun_gcrs[i] = r_sun * sun_place.r_hat[i];
          }

        /* convert sun position from GCRS to ITRS system
         * returns sun_itrs [km] */
        if ((error = novas_cel2ter (tt.jd_ut1, 0.0, tt.delta_t,
                                    method, sgt->accuracy, option,
                                    sgt->xpole, sgt->ypole, sun_gcrs,
                                    sun_itrs)) != 0)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_cel2ter",
                          __func__, error);
             return -1;
          }

        for (i = 0; i < 3; i++)
          {
             bs_sat[i] = sgt->bs_pos[i] - sgt->sat_pos[i];
             sun_sat[i] = sun_itrs[i] - sgt->sat_pos[i];
             /* sat_pos is used in the phi calculation below */
             sat_pos[i] = sgt->sat_pos[i];
          }
     }
   else
     {
        /* The solar boresight angle can also be done in GCRS coordinates,
         * but it's more complicated than in ITRS coordinates because the
         * boresight point and satellite point must be rotated from ITRS to GCRS.
         * The end result should be similar enough that it makes no difference.
         */

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

        /* convert boresight position from ITRS to GCRS system
         * returns bs_gcrs [km] */
        if ((error = novas_ter2cel (tt.jd_ut1, 0.0, tt.delta_t,
                                    method, sgt->accuracy, option,
                                    sgt->xpole, sgt->ypole, sgt->bs_pos,
                                    bs_gcrs)) != 0)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_ter2cel",
                          __func__, error);
             return -1;
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
             /* Copy into sat_pos for use below. */
             sat_pos[i] = sat_gcrs[i];
          }
     } /* if (Compute_Boresight_Sun_Angle_In_ITRS_Frame) */

   *ptheta = vec_angle (bs_sat, sun_sat) / DEGTORAD;

   if (earth_sun_distance)
     {
        *earth_sun_distance = r_sun;
     }

   if (pphi)
     {
        double dot, hat_bs_sat[3], hat_sat[3], u[3], hat_u[3];
        double hat_z[3], hat_vel[3], hat_slit[3], cross[3];
        double cos_phi;

        /* \H = host satellite position vector
         * \S = sun pos. vec.
         * \B = boresight pos. vec.
         * \z = earth rotation axis unit vector
         * \l = unit vector northward along slit
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
        vec_norm (sat_pos, hat_sat);                        /* \h */
        vec_cross (hat_z, hat_sat, hat_vel);                /* \v */
        vec_rotate (hat_z, hat_vel, tilt_angle, hat_slit);  /* \l */
        cos_phi = vec_dot (hat_u, hat_slit);                /* \u dot \l */
        *pphi = acos(cos_phi) / DEGTORAD;

        /* acos() returns a value in the range [0,pi].
         * Determine the sign of the angle from the z component of cross(\l, \u) */
        vec_cross (hat_slit, hat_u, cross);
        if (cross[2] < 0) *pphi *= -1;
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
                   sgt->bs_longitude /DEGTORAD,
                   sgt->bs_latitude /DEGTORAD);
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

static int sgt_initialize (Solar_Geom_Type *sgt)
{
   Novas_cat_entry_t dummy_star;
   double a = GEO_SAT_RADIUS / EARTH_MEAN_RADIUS;
   double tan_elev, tan_azi, delta_lon;
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

   /* Elevation and azimuth angles to point boresight at specified (lat,lon). */
   tan_elev = sin(sgt->bs_latitude) / (a - cos(sgt->bs_latitude));
   sgt->bs_elevation_angle = atan(tan_elev); /* radians */

   delta_lon = sgt->bs_longitude - sgt->sat_longitude;
   tan_azi = sin(delta_lon) / (a - cos(delta_lon));
   sgt->bs_azimuth_angle = atan(tan_azi);   /* radians */

   /* WGS84 coordinates of geostationary satellite */
   sgt->sat_pos[0] = GEO_SAT_RADIUS * cos(sgt->sat_longitude);  /* X */
   sgt->sat_pos[1] = GEO_SAT_RADIUS * sin(sgt->sat_longitude);  /* Y */
   sgt->sat_pos[2] = 0.0;                                       /* Z */

   /* WGS84 coordinates of boresight point */
   sgt->bs_pos[0] = EARTH_MEAN_RADIUS * cos(sgt->bs_longitude) * cos(sgt->bs_latitude);  /* X */
   sgt->bs_pos[1] = EARTH_MEAN_RADIUS * sin(sgt->bs_longitude) * cos(sgt->bs_latitude);  /* Y */
   sgt->bs_pos[2] = EARTH_MEAN_RADIUS * sin(sgt->bs_latitude);                           /* Z */

   /* NOVAS accuracy parameter (0=full) */
   sgt->accuracy = 0;

   return 0;
}

static void sgt_delete (Solar_Geom_Type *sgt)
{
   if (NULL == sgt)
     return;
   FREE(sgt);
}

static int init_method_cfg (Solar_Geom_Type *sgt, void *pv)
{
   config_t *cfg = (config_t *)pv;

   if ((0 != read_sat_config (cfg, &sgt->sat_longitude, &sgt->bs_longitude, &sgt->bs_latitude))
        || (0 != read_iers_params (sgt, cfg)))
     return -1;

   return 0;
}

static Solar_Geom_Type *solar_geom_init_using_method (int (*init_method)(Solar_Geom_Type *, void *),
                                                      void *pv)
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
   sgt->sgt_solar_xyz = sgt_solar_xyz;
   sgt->sgt_sat_sun_position = sgt_sat_sun_position;
   sgt->sgt_geosat_longitude = sgt_geosat_longitude;
   sgt->sgt_boresight_angles = sgt_boresight_angles;
   sgt->sgt_print_params = sgt_print_params;

   if (0 != init_method (sgt, pv))
     {
        sgt_delete (sgt);
        return NULL;
     }

   if (0 != sgt_initialize (sgt))
     {
        sgt_delete (sgt);
        return NULL;
     }

   return sgt;
}

Solar_Geom_Type *solar_geom_init (config_t *cfg)
{
   return solar_geom_init_using_method (init_method_cfg, cfg);
}

#ifdef UNIT_TEST

typedef struct
{
   double sat_longitude;
   double bs_longitude;
   double bs_latitude;
   double xpole;
   double ypole;
   double ut1_utc;
   double hour;
   short int year;
   short int month;
   short int day;
}
Sat_Config_Type;

static int init_method_test (Solar_Geom_Type *sgt, void *pv)
{
   Sat_Config_Type *sct = (Sat_Config_Type *)pv;

   sgt->sat_longitude = sct->sat_longitude * DEGTORAD;
   sgt->bs_longitude = sct->bs_longitude * DEGTORAD;
   sgt->bs_latitude = sct->bs_latitude * DEGTORAD;
   sgt->xpole = sct->xpole;
   sgt->ypole = sct->ypole;
   sgt->ut1_utc =sct->ut1_utc;

   return 0;
}

static Sat_Config_Type Config_Table[] =
{  /* sat_lon, bs_lon,bs_lat, xpole,ypole, ut1_utc, hr,year,month,day */
   {  0.0,    0.0,0.0,  0.0,0.0, 0.0,  0.0,2019,3,21},
   {  0.0,    0.0,0.0,  0.0,0.0, 0.0, 12.0,2019,3,21},
   {-90.0,  -90.0,0.0,  0.0,0.0, 0.0,  0.0,2019,3,21},
   {  0.0,    0.0,0.0,  0.0,0.0, 0.0,  0.0,2019,6,21},
   {-90.0,  -90.0,0.0,  0.0,0.0, 0.0,  0.0,2019,6,21},
   {  0.0,    0.0,0.0,  0.0,0.0, 0.0,  0.0,2019,7,15},
   {-90.0,  -90.0,0.0,  0.0,0.0, 0.0,  0.0,2019,7,15},
   {-90.0,  -90.0,0.0,  0.0,0.0, 0.0,  5.0,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,  0.0,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,  1.0,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,  2.0,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,  3.0,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,  4.0,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,    5.3361,2019,7,15},
   {-100.0, -100.0,0.0, 0.0,0.0, 0.0,    5.3361,2013,7,15},  /* Xiong's IRR simulation? e.g. omitted boresight tilt? */
   {-100.0, -94.44,0.0, 0.0,0.0, 0.0,    5.3361,2019,7,15},
   {-100.0, -94.44,34.095, 0.0,0.0, 0.0, 5.3361,2019,7,15},
   {-100.0, -94.44,34.095, 0.0,0.0, 0.0, 5.0715,2019,7,15},
   {-100.0, -100.0,34.095, 0.0,0.0, 0.0, 5.0185,2019,7,15},
};

struct sun_pos_type
{
   double theta, phi, dist;
};

static int fprint_sun_pos (FILE *fp, const Solar_Geom_Type *sgt, const Sat_Config_Type *sct,
                           const char *method_label, const struct sun_pos_type *pos)
{
   return fprintf (fp, "M-D-H: %02d-%02d-%7.3f  sat=%8.4f  bs:%8.4f,%8.4f  %s=> theta=%8.4f  phi=%8.4f, dist=%9.5e\n",
                   sct->month, sct->day, sct->hour,
                   sgt->sat_longitude/DEGTORAD,
                   sgt->bs_longitude/DEGTORAD, sgt->bs_latitude/DEGTORAD,
                   method_label, pos->theta, pos->phi, pos->dist);
}

int main (void)
{
   Solar_Geom_Type *sgt = NULL;
   double jd_begin, jd_end;
   short int de_number, error;
   int i, n;
   struct sun_pos_type itrs_method, gcrs_method;
   n = sizeof(Config_Table)/sizeof(Config_Table[0]);

   tio_time_set_taix_epoch ("2000-01-01T12:00:00Z");

   if ((error = novas_ephem_open ("/soft/tempo/sdpc/install/v1_gnu/ots/share/libnovas/JPLEPH",
                                  &jd_begin, &jd_end,
                                  &de_number)) != 0)
     {
        fprintf (stderr, "*** Error %d while opening ephemeris", error);
        return 1;
     }

   for (i = 0; i < n; i++)
     {
        double jd_utc;
        Sat_Config_Type *sct = &Config_Table[i];

        if (NULL == (sgt = solar_geom_init_using_method (init_method_test, sct)))
          {
             fprintf (stderr, "*** initialization failed (case = %d)\n", i);
             return 1;
          }
        jd_utc = novas_julian_date (sct->year, sct->month, sct->day, sct->hour);

        Compute_Boresight_Sun_Angle_In_ITRS_Frame = 1;
        if (0 != sgt->sgt_sat_sun_position (sgt, jd_utc, &itrs_method.theta, &itrs_method.phi, &itrs_method.dist))
          {
             fprintf (stderr, "*** sgt_sat_sun_position failed (case = %d)\n", i);
             return 1;
          }
        Compute_Boresight_Sun_Angle_In_ITRS_Frame = 0;
        if (0 != sgt->sgt_sat_sun_position (sgt, jd_utc, &gcrs_method.theta, &gcrs_method.phi, &gcrs_method.dist))
          {
             fprintf (stderr, "*** sgt_sat_sun_position failed (case = %d)\n", i);
             return 1;
          }
        (void) fprint_sun_pos (stdout, sgt, sct, "ITRS", &itrs_method);
        (void) fprint_sun_pos (stdout, sgt, sct, "GCRS", &gcrs_method);
        fprintf (stdout, "(GCRS-ITRS) delta:  theta = %f arcmin\n",
                 (gcrs_method.theta - itrs_method.theta) * 60.0);
        fprintf (stdout, "(GCRS-ITRS) delta:    phi = %f arcmin\n",
                 (gcrs_method.phi - itrs_method.phi) * 60.0);
        sgt->sgt_delete (sgt);
     }

   (void) novas_ephem_close();

   return 0;
}
#endif
