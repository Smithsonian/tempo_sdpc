/** @file solar.c
 *  @brief Interface for solar illumination geometry
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <libnovas.h>
#include <tempo_geo.h>

#include <tell.h>

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
   double sat_pos[3]; \
   double xpole; \
   double ypole; \
   double ut1_utc; \
   int leap_secs; \
   short int accuracy;
#include "solar.h"

#define KM_PER_AU      149597871.0
#define GEO_SAT_RADIUS     42163.968  /* km */
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

   if (0 != times_eval (&tt, jd_utc, sgt->leap_secs, sgt->ut1_utc))
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

static void vec_norm (double *a)
{
   double r = sqrt (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
   a[0] /= r;
   a[1] /= r;
   a[2] /= r;
}

static double angle_between_vectors (double *pa, double *pb)
{
   double cos_theta, a[3], b[3], c[3], c_hypot;
   int i;

   /* Assuming A and B are unit vectors, we have:
    *   AxB = sin(theta)
    *   A.B = cos(theta)
    * When A and B are nearly colinear, theta is small
    * and the most accurate answer is obtained from either
    * the sin or tangent:
    *    theta = atan2 ( |AxB|, A.B )
    * Otherwise, we can use the cosine:
    *    theta = acos (A.B)
    * Since the cosine is used in both cases, we compute it first.
    */

   memcpy ((char *)a, (char *)pa, 3 * sizeof(double));
   memcpy ((char *)b, (char *)pb, 3 * sizeof(double));

   vec_norm (a);
   vec_norm (b);

   /* dot product */
   cos_theta = 0.0;
   for (i = 0; i < 3; i++)
     {
        cos_theta += a[i] * b[i];
     }

   if ((-0.7 < cos_theta) && (cos_theta < 0.7))
     return acos (cos_theta);

   /* cross product */
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];

   c_hypot = sqrt (c[0]*c[0] + c[1]*c[1] * c[2]*c[2]);

   return atan2 (c_hypot, cos_theta);
}

static int sgt_sat_sun_angle (Solar_Geom_Type *sgt, double jd_utc,
                              double *psun_angle)
{
   Times_Type tt;
   Novas_sky_pos_t sun_place;
   double sat_gcrs[3];
   double bs_gcrs[3], bs_gcrs_vel[3];
   double bs_sat[3], sun_sat[3];
   double r_sun;
   short int error;
   short int coord_sys = 0;   /* 0 means GCRS coordinates */
   int method = 1; /* 1 = equinox-based method */
   int option = 0; /* 0 = output vector referred to GCRS axes */
   int i;

   if (0 != times_eval (&tt, jd_utc, sgt->leap_secs, sgt->ut1_utc))
     return -1;

   /* convert sat position from ITRS to GCRS system */
   if ((error = novas_ter2cel (tt.jd_ut1, 0.0, tt.delta_t,
                               method, sgt->accuracy, option,
                               sgt->xpole, sgt->ypole, sgt->sat_pos,
                               sat_gcrs)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_ter2cel",
                     __func__, error);
        return -1;
     }

   if ((error = novas_geo_posvel (tt.jd_tt, tt.delta_t, sgt->accuracy,
                                  &sgt->boresight_surface, bs_gcrs, bs_gcrs_vel)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: Error %d from novas_geo_posvel",
                     __func__, error);
        return -1;
     }

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

   *psun_angle = angle_between_vectors (bs_sat, sun_sat) / DEGTORAD;

   return 0;
}

static int sgt_geosat_longitude (const Solar_Geom_Type *sgt, double *lon)
{
   if (lon) *lon = sgt->sat_longitude;
   return 0;
}

static int sgt_print_params (Solar_Geom_Type *sgt, const char *pprefix, FILE *fp)
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

   if ((CONFIG_TRUE != config_setting_lookup_int (s, "leap_secs", &sgt->leap_secs))
       ||(CONFIG_TRUE != config_setting_lookup_float (s, "ut1_utc", &sgt->ut1_utc)))
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
   double bs_lon, bs_lat;
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

   if (0 != read_sat_config (cfg, &sgt->sat_longitude, &bs_lon, &bs_lat))
     return -1;

   /* WGS84 coordinates of geostationary satellite */
   sgt->sat_pos[0] = GEO_SAT_RADIUS * cos(sgt->sat_longitude);  /* X */
   sgt->sat_pos[1] = GEO_SAT_RADIUS * sin(sgt->sat_longitude);  /* Y */
   sgt->sat_pos[2] = 0.0;                            /* Z */

   novas_make_observer_on_surface (bs_lat, bs_lon, DEFAULT_HEIGHT,
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

   sgt->sgt_delete = sgt_delete;
   sgt->sgt_solar_zenith_angle = sgt_solar_zenith_angle;
   sgt->sgt_sat_sun_angle = sgt_sat_sun_angle;
   sgt->sgt_geosat_longitude = sgt_geosat_longitude;
   sgt->sgt_print_params = sgt_print_params;

   if (0 != sgt_initialize (sgt, cfg))
     {
        sgt_delete (sgt);
        return NULL;
     }

   return sgt;
}
