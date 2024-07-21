/** @file granule.c
 *  @brief Load granule geolocation information,
 *         set geolocation-dependent variables
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <wordexp.h>
#include <proj_api.h>

#include <tempo_geo.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"

#define BIT_CLEAR(uc,n)  ((uc) &= ~(1 << (n)))
#define BIT_SET(uc,n)    ((uc) |=  (1 << (n)))

#define BITMASK_CLEAR(f,mask)  ((f) &= (~(mask)))
#define BITMASK_SET(f,mask)    ((f) |=   (mask) )

#define BITMASK_GPQF_BITS_USED   0xffffffff

#define PROJ_ARGS_BUFSIZE       80
#define GEO_ALTITUDE  35785831.0   /* meters */

#define ATTR_STRING_TERRAIN_REFERENCED "terrain_referenced_coordinates"

/* To compute the parallax shift distance and write it to the output file: */
/* #define OUTPUT_PARALLAX_SHIFT 1 */
#undef OUTPUT_PARALLAX_SHIFT

typedef struct
{
   double *lon;
   double *lat;
   double *lon_cnr;
   double *lat_cnr;
   double *shift_km;
   int *inr_quality_flag;
   unsigned int num_mirror_step;
   unsigned int num_xtrack;
   unsigned int num_corner;
   int dimids[3];
   int group;
}
Geoloc_Type;

/* ECEF = Earth-centered Earth-fixed coordinate system.
 * See WGS84 reference frame (X,Y,Z).
 * Origin is Earth center of mass.
 * X,Y axes in the equatorial plane.
 * +X axis lies along the prime meridian, longitude=0
 * +Z axis lies along the rotation axis, +Z toward the North.
 * +Y axis completes the right-handed coordinate system.
 */
typedef struct
{
   double *X;
   double *Y;
   double *Z;
   unsigned int num;
}
ECEF_Position_Type;

typedef struct
{
   double *zenith_angle;
   double *azimuth_angle;
   unsigned int num_mirror_step;
   unsigned int num_xtrack;
}
Angles_Type;

#define NUM_BANDS 2

#define GRANULE_PRIVATE_DATA \
   int ncid; \
   Geoloc_Type geoloc[NUM_BANDS]; \
   Angles_Type *sun_angles[NUM_BANDS]; \
   Angles_Type *sat_angles[NUM_BANDS]; \
   ECEF_Position_Type sun; \
   ECEF_Position_Type moon; \
   ECEF_Position_Type sat;
#include "granule.h"

typedef struct
{
   AltitudeData geoid_height;
   AltitudeData dem;
}
Geoid_Data_Type;

static __inline__ int invalid_lonlat (double lon, double lat)
{
   return ((0 == isfinite(lon)) || (0 == isfinite(lat)) || (fabs(lon) > 360.0) || (fabs(lat) > 90.0));
}
static __inline__ int invalid_point (double x, double y)
{
   return ((0 == isfinite(x)) || (0 == isfinite(y)));
}

static void free_geoloc_fields (Geoloc_Type *geoloc)
{
   if (geoloc == NULL)
     return;
   FREE(geoloc->lon);
   FREE(geoloc->lat);
   FREE(geoloc->inr_quality_flag);
   FREE(geoloc->shift_km);
}

static double *alloc_geoloc_coordinate (Geoloc_Type *geoloc)
{
   size_t num_pixels = geoloc->num_mirror_step * geoloc->num_xtrack;
   size_t alloc_len = num_pixels * (1 + geoloc->num_corner);
   double *x;

   if (NULL == (x = (double *)MALLOC (alloc_len * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)x, 0, alloc_len * sizeof(double));
   return x;
}

static int alloc_geoloc_fields (size_t *dimlens, Geoloc_Type *geoloc)
{
   unsigned int num_pixels;

   geoloc->num_mirror_step = dimlens[0];
   geoloc->num_xtrack = dimlens[1];
   geoloc->num_corner = dimlens[2];

   num_pixels = geoloc->num_mirror_step * geoloc->num_xtrack;

   if ((NULL == (geoloc->lon = alloc_geoloc_coordinate (geoloc)))
       || (NULL == (geoloc->lat = alloc_geoloc_coordinate (geoloc)))
       || (NULL == (geoloc->inr_quality_flag = (int *)MALLOC (num_pixels * sizeof(int))))
      )
     {
        free_geoloc_fields (geoloc);
        return -1;
     }
   memset ((char *)geoloc->inr_quality_flag, 0, num_pixels * sizeof(int));

#ifdef OUTPUT_PARALLAX_SHIFT
   if (NULL == (geoloc->shift_km = alloc_geoloc_coordinate (geoloc)))
     {
        free_geoloc_fields (geoloc);
        return -1;
     }
#endif

   num_pixels = geoloc->num_mirror_step * geoloc->num_xtrack;
   geoloc->lon_cnr = geoloc->lon + num_pixels;
   geoloc->lat_cnr = geoloc->lat + num_pixels;

   return 0;
}

static void free_ecef_position (ECEF_Position_Type *p)
{
   if (p == NULL)
     return;
   FREE(p->X);
}

static int alloc_ecef_position (unsigned int num, ECEF_Position_Type *p)
{
   double *x = NULL;

   if (NULL == (x = (double *)MALLOC (3 * num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   memset ((char *)x, 0, 3 * num * sizeof(double));

   p->X = x;
   p->Y = x + num;
   p->Z = x + num * 2;
   p->num = num;

   return 0;
}

static int def_elevation_vars (Geoloc_Type *geoloc)
{
   TIO_Var_Info_Type info;
   TIO_Attr_Text_Type hgt_text_attrs[] =
     {
        {"units", "m"},
        {"long_name", TEMPO_VAR_TERR_HEIGHT},
        {"comment", "Area-weighted mean terrain height within pixel boundary"},
        {"coordinates", "time longitude latitude"},
        {NULL, NULL}
     };
   TIO_Attr_Text_Type dev_text_attrs[] =
     {
        {"units", "m"},
        {"long_name", TEMPO_VAR_TERR_HEIGHT_STDDEV},
        {"comment", "Area-weighted terrain height standard deviation within pixel boundary"},
        {"coordinates", "time longitude latitude"},
        {NULL, NULL}
     };
   short fill_value = TIO_FILL_SHORT;
   int varid, status;

   tell_push_queue();
   status = TIO_inq_var (geoloc->group, TEMPO_VAR_TERR_HEIGHT, &info);
   tell_pop_queue(1);
   if (status != 0)
     {
        if ((0 != TIO_def_var (geoloc->group, TEMPO_VAR_TERR_HEIGHT, NC_SHORT, 2, geoloc->dimids, &varid))
            || (0 != TIO_def_var_fill (geoloc->group, varid, 0, &fill_value))
            || (0 != TIO_put_text_attrs (geoloc->group, varid, hgt_text_attrs)))
          return -1;
     }

   tell_push_queue();
   status = TIO_inq_var (geoloc->group, TEMPO_VAR_TERR_HEIGHT_STDDEV, &info);
   tell_pop_queue(1);
   if (status != 0)
     {
        if ((0 != TIO_def_var (geoloc->group, TEMPO_VAR_TERR_HEIGHT_STDDEV, NC_SHORT, 2, geoloc->dimids, &varid))
            || (0 != TIO_def_var_fill (geoloc->group, varid, 0, &fill_value))
            || (0 != TIO_put_text_attrs (geoloc->group, varid, dev_text_attrs)))
          return -1;
     }

   return 0;
}

static int set_elevation (Granule_Type *gt, const Elevation_Type *et)
{
   Geoloc_Type *geoloc;
   short *shrt_alloc = NULL;
   short *height = NULL;
   short *stddev = NULL;
   size_t num_pixels, len;
   int start[3], count[3];
   int i, status;

   geoloc = &gt->geoloc[0];
   num_pixels = geoloc->num_mirror_step * geoloc->num_xtrack;
   len = num_pixels * geoloc->num_corner;

   if (NULL == (shrt_alloc = (short *) MALLOC (2 * len * sizeof(short))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)shrt_alloc, 0, 2 * len * sizeof(short));

   height = shrt_alloc;
   stddev = shrt_alloc + len;

   for (i = 0; i < NUM_BANDS; i++)
     {
        geoloc = &gt->geoloc[i];

        if (0 != et->et_regrid (et, num_pixels, geoloc->lon_cnr, geoloc->lat_cnr,
                                height, stddev))
          goto return_error;

        if (0 != def_elevation_vars (geoloc))
          goto return_error;

        start[0] = 0;
        start[1] = 0;
        count[0] = geoloc->num_mirror_step;
        count[1] = geoloc->num_xtrack;

        if (0 != TIO_put_var_section (geoloc->group, TEMPO_VAR_TERR_HEIGHT,
                                       start, count, TIO_SHORT, height))
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: setting pixel center terrain height", __func__);
             goto return_error;
          }
        if (0 != TIO_put_var_section (geoloc->group, TEMPO_VAR_TERR_HEIGHT_STDDEV,
                                       start, count, TIO_SHORT, stddev))
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: setting pixel center terrain height stddev", __func__);
             goto return_error;
          }
     }

   status = 0;

return_error:
   FREE(shrt_alloc);
   return status;
}

static int def_snow_ice_fraction_var (Geoloc_Type *geoloc)
{
   TIO_Var_Info_Type info;
   TIO_Attr_Text_Type text_attrs[] =
     {
        {"long_name", TEMPO_VAR_SNOWICE_FRACTION},
        {"comment", "Fraction of pixel area covered by snow and/or ice"},
        {"coordinates", "time longitude latitude"},
        {NULL, NULL}
     };
   float fill_value = TIO_FILL_FLOAT;
   int varid, status;

   tell_push_queue();
   status = TIO_inq_var (geoloc->group, TEMPO_VAR_SNOWICE_FRACTION, &info);
   tell_pop_queue(1);
   if (status != 0)
     {
        if ((0 != TIO_def_var (geoloc->group, TEMPO_VAR_SNOWICE_FRACTION, NC_FLOAT, 2, geoloc->dimids, &varid))
            || (0 != TIO_def_var_fill (geoloc->group, varid, 0, &fill_value))
            || (0 != TIO_put_text_attrs (geoloc->group, varid, text_attrs)))
          return -1;
     }

   return 0;
}

static int set_snow_ice_fraction (Granule_Type *gt, Snow_Type *sn)
{
   Geoloc_Type *geoloc;
   float *snow_ice_fraction = NULL;
   size_t num_pixels;
   int start[3], count[3];
   int i, status;

   geoloc = &gt->geoloc[0];
   num_pixels = geoloc->num_mirror_step * geoloc->num_xtrack;

   if (NULL == (snow_ice_fraction = (float *) MALLOC (num_pixels * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < NUM_BANDS; i++)
     {
        geoloc = &gt->geoloc[i];

        if (0 != sn->sn_regrid (sn, num_pixels, geoloc->lon_cnr, geoloc->lat_cnr,
                                snow_ice_fraction))
          goto return_error;

        if (0 != def_snow_ice_fraction_var (geoloc))
          goto return_error;

        start[0] = 0;
        start[1] = 0;
        count[0] = geoloc->num_mirror_step;
        count[1] = geoloc->num_xtrack;

        if (0 != TIO_put_var_section (geoloc->group, TEMPO_VAR_SNOWICE_FRACTION,
                                       start, count, TIO_FLOAT, snow_ice_fraction))
          {
             tell_verror (TELL_RUNTIME_ERROR,
                          "%s: setting pixel snow/ice fraction", __func__);
             goto return_error;
          }
     }

   status = 0;

return_error:
   FREE(snow_ice_fraction);
   return status;
}

/* Derivation of the solar glint angle:
 * Notation:
 * Polar unit vector:  u(theta,phi) = Xhat * sin(theta)*cos(phi)
 *                                  + Yhat * sin(theta)*sin(phi)
 *                                  + Zhat * cos(theta)
 * where Xhat, Yhat, Zhat are unit vectors along coordinate axes
 * with origin at earth point, P.
 * In the following, abbreviate theta=t, phi=p.
 * Define s(t0,p0) = unit vector from earth point, P, toward sun
 *        v(t,p)   = unit vector from earth point, P, toward satellite
 * Incident ray from the sun travels along (-s).
 * From the law of reflection, the reflected solar ray, r, travels
 * from earth point, P, along:
 *          r(t0,p0) = (-s(t0,p0)) + 2*cos(t0)*Zhat
 *
 * The dot product, r \dot v, yields the cosine, mu, of the angle
 * between r and v:
 *
 *   mu = r \dot v
 *      = (-sx0)*vx + (-sy0)*vy + (-sz0 + 2*cos(t0))*vz
 *      = ( (-sin(t0)*cos(p0))*(sin(t)*cos(p))
 *        + (-sin(t0)*sin(p0))*(sin(t)*sin(p))
 *        + (+cos(t0)        )*(cos(t)))
 *      = - sin(t0)*sin(t)*( cos(p0)*cos(p) + sin(p0)*sin(p) )
 *        + cos(t0)*cos(t)
 *   mu = cos(t0)*cos(t) - sin(t0)*sin(t)*cos(p-p0)
 *
 * To check this, consider the case where p-p0=PI, t=t0:
 *   mu = cos(t)^2 + sin(t)^2 = 1  => glint_angle = 0,
 * e.g. the reflected solar ray is aimed at the satellite.
 * Also, when t=t0=0, mu=1 independent of p,p0, as expected.
 *
 * Note that for small zenith angles, this is an ill-conditioned expression
 * for determining the glint angle, g.  A better expression can
 * be obtained by using the identity cos(g) = 1 - 2 sin^2(g/2) to obtain:
 *   sin^2(g/2) = 1/2[1 - cos(t0)cos(t) + sin(t0)sin(t)cos(p-p0)].
 * However, in this application, we only need to determine
 * whether or not g<g_threshold, where g_threshold is some large angle
 * like 40 degrees.  The simpler expression is good enough for that.
 */
static double cos_sun_glint_angle (double sza, double saa, double vza, double vaa)
{
   double deg2rad = M_PI/180.0;
   double rel_azimuth;

   /* zenith angles deg -> radians */
   sza *= deg2rad;
   vza *= deg2rad;

   /* ensure azimuth angles span 0-2*PI */
   saa = fmod (saa, 360.0);
   if (saa < 0) saa += 360.0;
   vaa = fmod (vaa, 360.0);
   if (vaa < 0) vaa += 360.0;

   /* relative azimuth angle */
   rel_azimuth = (vaa - saa) * deg2rad;

   return cos(sza)*cos(vza) - sin(sza)*sin(vza)*cos(rel_azimuth);
}

static int set_sun_glint_bit (const Angles_Type *sun_angles,
                               const Angles_Type *sat_angles,
                               double max_glint_angle,
                               unsigned char *illum_flags)
{
   double *sza = sun_angles->zenith_angle;
   double *saa = sun_angles->azimuth_angle;
   double *vza = sat_angles->zenith_angle;
   double *vaa = sat_angles->azimuth_angle;
   double cos_max_glint = cos (max_glint_angle*M_PI/180);
   size_t i, num_pixels = sun_angles->num_mirror_step * sun_angles->num_xtrack;

   for (i = 0; i < num_pixels; i++)
     {
        double mu_glint = cos_sun_glint_angle (sza[i], saa[i], vza[i], vaa[i]);
        if (mu_glint > cos_max_glint)
          {
             BIT_SET(illum_flags[i], 0);
          }
        else BIT_CLEAR(illum_flags[i], 0);
     }

   return 0;
}

static void unit_vector (double dx, double dy, double dz, double v[3])
{
   double len = sqrt(dx*dx + dy*dy + dz*dz);
   v[0] = dx/len;
   v[1] = dy/len;
   v[2] = dz/len;
};

static int set_solar_eclipse_bit (const Geoloc_Type *geoloc,
                                  const ECEF_Position_Type *sun,
                                  const ECEF_Position_Type *moon,
                                  double max_eclipse_angle,
                                  unsigned char *illum_flags)
{
   double cos_max_eclipse = cos (max_eclipse_angle * M_PI/180.0);
   size_t step, i;

   for (step = 0; step < geoloc->num_mirror_step; step++)
     {
        ECEF_Vector s, m, vec;
        unsigned char *illum_flags_row = illum_flags + step * geoloc->num_xtrack;
        double *lon_row = geoloc->lon + step * geoloc->num_xtrack;
        double *lat_row = geoloc->lat + step * geoloc->num_xtrack;

        s.theX =  sun->X[step]; s.theY =  sun->Y[step]; s.theZ =  sun->Z[step];
        m.theX = moon->X[step]; m.theY = moon->Y[step]; m.theZ = moon->Z[step];

        for (i = 0; i < geoloc->num_xtrack; i++)
          {
             TempoGeoErr err;
             EarthPoint pt;
             double us[3], um[3];
             double mu_eclipse;

             BIT_CLEAR(illum_flags_row[i], 1);

             /* input may include fill values */
             if (invalid_lonlat (lon_row[i], lat_row[i]))
               continue;

             pt.theLon = lon_row[i];
             pt.theLat = lat_row[i];
             pt.theAlt = 0.0;

             if (0 != (err = computeSiteVector (&pt, &vec)))
               {
                  const char *err_string = tempoGeoErrorString(err);
                  if (err_string == NULL) err_string = "unknown error";
                  tell_verror (TELL_RUNTIME_ERROR, "%s: computeSiteVector failed (%s)",
                               __func__, err_string);
                  return -1;
               }

             unit_vector (s.theX - vec.theX,
                          s.theY - vec.theY,
                          s.theZ - vec.theZ, us);
             unit_vector (m.theX - vec.theX,
                          m.theY - vec.theY,
                          m.theZ - vec.theZ, um);

             mu_eclipse = us[0]*um[0] + us[1]*um[1] + us[2]*um[2];

             if (mu_eclipse > cos_max_eclipse)
               {
                  BIT_SET(illum_flags_row[i], 1);
               }
          }
     }

   return 0;
}

static int set_ground_pixel_flags (Granule_Type *gt,
                                   double max_glint_angle,
                                   double max_eclipse_angle,
                                   const Land_Cover_Type *lc)
{
   Geoloc_Type *geoloc;
   unsigned int *ground_flags = NULL;
   unsigned char *ubytes = NULL;
   unsigned char *lc_type1, *lc_typeqc, *illum_flags;
   int *inr_quality_flag = NULL;
   size_t num_pixels, ubytes_size;
   int start[2], count[2];
   int i, status = -1;

   geoloc = &gt->geoloc[0];
   num_pixels = geoloc->num_mirror_step * geoloc->num_xtrack;

   ubytes_size = 4 * num_pixels * sizeof(unsigned char);
   if ((NULL == (ground_flags = (unsigned int *) MALLOC (num_pixels * sizeof(unsigned int))))
       || (NULL == (inr_quality_flag = (int *) MALLOC (num_pixels * sizeof(int))))
       || (NULL == (ubytes = (unsigned char *) MALLOC (ubytes_size))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_error;
     }
   memset ((char *)ubytes, 0, ubytes_size);

   lc_type1    = ubytes + num_pixels;
   lc_typeqc   = ubytes + num_pixels * 2;
   illum_flags = ubytes + num_pixels * 3;

   start[0] = 0;
   start[1] = 0;
   count[0] = geoloc->num_mirror_step;
   count[1] = geoloc->num_xtrack;

   for (i = 0; i < NUM_BANDS; i++)
     {
        unsigned int j;
        int varid_gpqf;

        geoloc = &gt->geoloc[i];

        if (0 != tio_inq_varid (geoloc->group, TEMPO_VAR_GROUND_PIXEL_QF, &varid_gpqf))
          {
             if (0 != tio_def_var_ground_pixel_quality_flag (geoloc->group))
               goto return_error;
          }

        if (0 != TIO_get_var_section (geoloc->group, TEMPO_VAR_GROUND_PIXEL_QF,
                                      start, count, TIO_UINT, ground_flags))
          goto return_error;

        if (0 != TIO_get_var_section (geoloc->group, TEMPO_VAR_INRQF,
                                      start, count, TIO_INT, inr_quality_flag))
          goto return_error;

        if ((0 != lc->lc_lookup_type1 (lc, num_pixels, geoloc->lon, geoloc->lat, lc_type1))
            || (0 != lc->lc_lookup_typeqc (lc, num_pixels, geoloc->lon, geoloc->lat, lc_typeqc)))
          goto return_error;

        if ((0 != set_sun_glint_bit (gt->sun_angles[i], gt->sat_angles[i], max_glint_angle, illum_flags))
            || (0 != set_solar_eclipse_bit (geoloc, &gt->sun, &gt->moon, max_eclipse_angle, illum_flags)))
          goto return_error;

        /* Clear bits to be set, preserving the other bits as-is.
         * bit 0 is the least significant bit.
         * bit 0-3 = MODIS land/water mask.
         * bit 4-5 are illumination flags (bit 4=sun glint possibility,
         *                                 bit 5=solar eclipse possibility).
         * bit 6 is the INR quality flag
         * bits 7-15 are currently unused.
         * bits 16-23 = 8-bit MODIS yearly land cover flags, MCD12Q1, IGBP Type 1.
         */

        for (j = 0; j < num_pixels; j++)
          {
             unsigned int flags = ground_flags[j];
             unsigned int inrqf = (inr_quality_flag[j] == 0) ? 0 : 1;

             BITMASK_CLEAR(flags, BITMASK_GPQF_BITS_USED);

             BITMASK_SET(flags,   lc_typeqc[j] >> 4);
             BITMASK_SET(flags, illum_flags[j] << 4);
             BITMASK_SET(flags,          inrqf << 6);
             BITMASK_SET(flags,    lc_type1[j] << 16);

             ground_flags[j] = flags;
          }

        if (0 != TIO_put_var_section (geoloc->group, TEMPO_VAR_GROUND_PIXEL_QF,
                                      start, count, TIO_UINT, ground_flags))
          goto return_error;
     }

   status = 0;
return_error:
   FREE(ubytes);
   FREE(ground_flags);
   FREE(inr_quality_flag);

   return status;
}

static int set_earth_sun_distance (Granule_Type *gt)
{
   const ECEF_Position_Type *sun = &gt->sun;
   double x = sun->X[0], y = sun->Y[0], z = sun->Z[0];
   double dist_km = sqrt (x*x + y*y + z*z);
   return tio_set_earth_sun_distance (gt->ncid, dist_km * 1.e3);
}

static void free_angles_type (Angles_Type *at)
{
   if (at == NULL)
     return;
   FREE(at->zenith_angle);
   FREE(at);
}

static Angles_Type *new_angles_type (int num_mirror_step, int num_xtrack)
{
   Angles_Type *at = NULL;
   double *a = NULL;
   size_t num = num_mirror_step * num_xtrack;

   if ((NULL == (at = (Angles_Type *)MALLOC (sizeof *at)))
       || (NULL == (a = (double *)MALLOC (2 * num * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)at, 0, sizeof (*at));
   memset ((char *)a, 0, 2 * num * sizeof(double));

   at->zenith_angle = a;
   at->azimuth_angle = a + num;
   at->num_mirror_step = num_mirror_step;
   at->num_xtrack = num_xtrack;

   return at;
}

static int object_angles (const ECEF_Vector *object,
                          int num, const double *lon, const double *lat,
                          double *zenith_angle, double *azimuth_angle)
{
   int i;

   for (i = 0; i < num; i++)
     {
        TempoGeoErr err;
        EarthPoint pt;
        AzZen angles;

        zenith_angle[i] = TIO_FILL_FLOAT;
        azimuth_angle[i] = TIO_FILL_FLOAT;

        /* input may include fill values */
        if (invalid_lonlat (lon[i], lat[i]))
          continue;

        pt.theLon = lon[i];
        pt.theLat = lat[i];
        pt.theAlt = 0.0;

        if (0 != (err = computeViewingAngles (&pt, object, &angles)))
          {
             const char *err_string = tempoGeoErrorString(err);
             if (err_string == NULL) err_string = "unknown error";
             tell_verror (TELL_RUNTIME_ERROR, "%s: computeViewingAngles failed (%s)",
                          __func__, err_string);
             return -1;
          }

        zenith_angle[i] = angles.theZen;
        azimuth_angle[i] = angles.theAz;
     }

   return 0;
}

static Angles_Type *map_object_angles (const Geoloc_Type *geoloc,
                                       const ECEF_Position_Type *object)
{
   Angles_Type *at = NULL;
   int step, num_step, num_xtrack;

   num_step = geoloc->num_mirror_step;
   num_xtrack = geoloc->num_xtrack;

   if (NULL == (at = new_angles_type (num_step, num_xtrack)))
     return NULL;

   for (step = 0; step < num_step; step++)
     {
        ECEF_Vector obj;
        double *lon, *lat;
        double *za, *aa;

        lon = geoloc->lon + step * num_xtrack;
        lat = geoloc->lat + step * num_xtrack;
        za = at->zenith_angle + step * num_xtrack;
        aa = at->azimuth_angle + step * num_xtrack;

        obj.theX = object->X[step];
        obj.theY = object->Y[step];
        obj.theZ = object->Z[step];

        if (0 != object_angles (&obj, num_xtrack, lon, lat, za, aa))
          {
             free_angles_type (at);
             return NULL;
          }
     }

   return at;
}

static int write_object_angles (int grp, const Angles_Type *at, const char **varnames)
{
   int start[2], count[2];

   start[0] = 0;
   start[1] = 0;
   count[0] = at->num_mirror_step;
   count[1] = at->num_xtrack;

   if (0 != TIO_put_var_section (grp, varnames[0], start, count, NC_DOUBLE, at->zenith_angle))
     return -1;
   if (0 != TIO_put_var_section (grp, varnames[1], start, count, NC_DOUBLE, at->azimuth_angle))
     return -1;

   return 0;
}

static int is_twilight_granule (Granule_Type *gt, int *is_radt)
{
#define MAX_PRODUCT_TYPE_LEN 16
   char product_type[MAX_PRODUCT_TYPE_LEN];
   size_t len;
   int xtype;

   *is_radt = 0;

   if ((NC_NOERR != nc_inq_att (gt->ncid, NC_GLOBAL, "product_type", &xtype, &len))
       && (len >= MAX_PRODUCT_TYPE_LEN))
     return 0;

   if (0 != TIO_get_att (gt->ncid, NC_GLOBAL, "product_type", NC_CHAR, product_type))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading global attribute: product_type", __func__);
        return -1;
     }

   if (0 != strcmp (product_type, TEMPO_PROD_TYPE_RAD_TWI))
     return 0;

   *is_radt = 1;
   return 0;
}

static int set_exposure_time_valid_max (Granule_Type *gt, float valid_max_exposure_time)
{
   int varid;

   if ((0 != tio_inq_varid (gt->ncid, TEMPO_VAR_EXPOSURE_TIME, &varid))
       || (0 != TIO_put_att (gt->ncid, varid, "valid_max", NC_FLOAT, 1, &valid_max_exposure_time)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: setting %s variable attribute: valid_max=%f",
                     __func__, TEMPO_VAR_EXPOSURE_TIME, valid_max_exposure_time);
        return -1;
     }

   return 0;
}

static int set_object_angles (Granule_Type *gt, int is_radt)
{
   const char *sun_angle_names[] = {TEMPO_VAR_SZ_ANGLE, TEMPO_VAR_SA_ANGLE};
   const char *sat_angle_names[] = {TEMPO_VAR_VZ_ANGLE, TEMPO_VAR_VA_ANGLE};
   int i;

   for (i = 0; i < NUM_BANDS; i++)
     {
        Geoloc_Type *geoloc = &gt->geoloc[i];
        int varid_sza;

        if (NULL == (gt->sun_angles[i] = map_object_angles (geoloc, &gt->sun)))
          return -1;
        if (NULL == (gt->sat_angles[i] = map_object_angles (geoloc, &gt->sat)))
          return -1;

        /* If SZA doesn't exist, assume that means we need to define all the angle
         * variables.  If SZA does exist, assume that means all the angle variables
         * have been correctly defined already. */
        if (0 != tio_inq_varid (geoloc->group, TEMPO_VAR_SZ_ANGLE, &varid_sza))
          {
             if (0 != tio_def_l1_radiance_angle_vars (geoloc->group))
               return -1;
          }

        if (0 != write_object_angles (geoloc->group, gt->sun_angles[i], sun_angle_names))
          return -1;
        if (0 != write_object_angles (geoloc->group, gt->sat_angles[i], sat_angle_names))
          return -1;

        if (is_radt)
          {
             float valid_max_sza = 180.0;
             if ((0 != tio_inq_varid (geoloc->group, TEMPO_VAR_SZ_ANGLE, &varid_sza))
                 || (0 != TIO_put_att (geoloc->group, varid_sza, "valid_max", NC_FLOAT, 1, &valid_max_sza)))
               {
                  tell_verror (TELL_RUNTIME_ERROR, "%s: setting %s variable attribute: valid_max=%f",
                               __func__, TEMPO_VAR_SZ_ANGLE, valid_max_sza);
                  return -1;
               }
          }
     }

   return 0;
}

static void delete_granule_type (Granule_Type *gt)
{
   int i;

   if (gt == NULL)
     return;

   TIO_close (gt->ncid);

   free_ecef_position (&gt->sun);
   free_ecef_position (&gt->moon);
   free_ecef_position (&gt->sat);

   for (i = 0; i < NUM_BANDS; i++)
     {
        free_geoloc_fields (&gt->geoloc[i]);
        free_angles_type (gt->sun_angles[i]);
        free_angles_type (gt->sat_angles[i]);
     }

   FREE(gt);
}

static Granule_Type *new_granule_type (void)
{
   Granule_Type *gt = NULL;

   if (NULL == (gt = (Granule_Type *)MALLOC (sizeof *gt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)gt, 0, sizeof *gt);

   gt->gt_close = delete_granule_type;
   gt->gt_set_elevation = set_elevation;
   gt->gt_set_snow_ice_fraction = set_snow_ice_fraction;
   gt->gt_set_object_angles = set_object_angles;
   gt->gt_set_ground_pixel_flags = set_ground_pixel_flags;
   gt->gt_set_earth_sun_distance = set_earth_sun_distance;
   gt->gt_is_twilight_granule = is_twilight_granule;
   gt->gt_set_exposure_time_valid_max = set_exposure_time_valid_max;

   return gt;
}

#ifdef OUTPUT_PARALLAX_SHIFT
static int write_shift (const Granule_Type *gt, const char *band_name,
                        const Geoloc_Type *geoloc)
{
   const char *varname = "parallax_shift";
   static TIO_Attr_Text_Type attrs[] =
     {
        {"coordinates", "time longitude latitude"},
        {NULL,NULL}
     };
   int grp, varid, start[2], count[2];

   if (0 != TIO_inq_grp (gt->ncid, band_name, &grp))
     return -1;
   if (0 != TIO_def_var (grp, varname, NC_DOUBLE, 2, geoloc->dimids, &varid))
     return -1;

   start[0] = 0;
   start[1] = 0;
   count[0] = geoloc->num_mirror_step;
   count[1] = geoloc->num_xtrack;

   if (0 != TIO_put_var_section (grp, varname, start, count, TIO_DOUBLE, geoloc->shift_km))
     return -1;

   if (0 != TIO_put_text_attrs (grp, varid, attrs))
     return -1;

   return 0;
}
#endif

static int write_band_geolocation (const Granule_Type *gt, const char *band_name,
                                   const Geoloc_Type *geoloc, int terrain_referenced)
{
   const char *attname = ATTR_STRING_TERRAIN_REFERENCED;
   const char *yes = "yes";
   int grp, start[3], count[3];

   if (0 != TIO_inq_grp (gt->ncid, band_name, &grp))
     return -1;

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = geoloc->num_mirror_step;
   count[1] = geoloc->num_xtrack;
   count[2] = geoloc->num_corner;

   if ((0 != TIO_put_var_section (grp, TEMPO_VAR_LONGITUDE, start, count, TIO_DOUBLE, geoloc->lon))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_LATITUDE, start, count, TIO_DOUBLE, geoloc->lat))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS, start, count, TIO_DOUBLE, geoloc->lon_cnr))
       || (0 != TIO_put_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS, start, count, TIO_DOUBLE, geoloc->lat_cnr)))
     {
        return -1;
     }

#ifdef OUTPUT_PARALLAX_SHIFT
   if (0 != write_shift (gt, band_name, geoloc))
     return -1;
#endif

   if (terrain_referenced)
     {
        if (0 != TIO_put_att (grp, NC_GLOBAL, attname, NC_CHAR, 1+strlen(yes), yes))
          return -1;
     }

   return 0;
}

static void free_geoid_data (Geoid_Data_Type *gdt)
{
   if (gdt == NULL)
     return;
   freeAltitudeData (&gdt->geoid_height);
   freeAltitudeData (&gdt->dem);
}

static int read_geoid_data (Geoid_Data_Type *gdt, const char *geoid_dem_path)
{
   if (1 != readGeoidHeightFromFile (geoid_dem_path, &gdt->geoid_height))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading geoid height from file: %s",
                     __func__, geoid_dem_path);
        return -1;
     }

   if (1 != readDEMFromFile (geoid_dem_path, &gdt->dem))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading digital elevation model from file: %s",
                     __func__, geoid_dem_path);
        return -1;
     }

   return 0;
}

static __inline__ double compute_shift (double lat0, double lon0, double lat1, double lon1)
{
   double deg2rad = M_PI/180.0;
   double phi0 = lat0 * deg2rad;
   double th0  = lon0 * deg2rad;
   double phi1 = lat1 * deg2rad;
   double th1  = lon1 * deg2rad;
   /* haversine forumula */
   double sphi = sin((phi1-phi0)*0.5);
   double sth  = sin(( th1- th0)*0.5);
   double delta = 2 * asin (sqrt(sphi*sphi + cos(phi0)*cos(phi1)*sth*sth));
   return delta * 6371.0088;
}

/* This code will only ever process polygons with 4 vertices,
 * so use static temporary arrays for the sort. */
static struct
{
   double x[4];
   double y[4];
   double angles[4];
   int index[4];
}
Angle_Sort;

static int compare_angle_indices (const void *va, const void *vb)
{
   int ia = *(int *)va;
   int ib = *(int *)vb;
   double angle_a = Angle_Sort.angles[ia];
   double angle_b = Angle_Sort.angles[ib];
   if (angle_a < angle_b) return -1;
   if (angle_a > angle_b) return +1;
   return 0;
}

static int polygon_sort_ccw (double *x, double *y, int n, int **pindices)
{
   double xc, yc;
   int i;

   if (n != 4)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: expected polygon with 4 vertices!! num_vertices=%d", __func__, n);
        return -1;
     }

   xc = 0.0;
   yc = 0.0;

   for (i = 0; i < n; i++)
     {
        xc += x[i];
        yc += y[i];
     }

   xc /= n;
   yc /= n;

   /* define angles so that the sort yields points ordered NE,NW,SW,SE */
   for (i = 0; i < n; i++)
     {
        double theta = atan2 (-(y[i]-yc), -(x[i]-xc)) + M_PI;
        Angle_Sort.angles[i] = theta;
        Angle_Sort.index[i] = i;
        Angle_Sort.x[i] = x[i];
        Angle_Sort.y[i] = y[i];
     }

   qsort (Angle_Sort.index, (size_t)n, sizeof(int), compare_angle_indices);

   if (pindices)
     {
        *pindices = Angle_Sort.index;
     }

   for (i = 0; i < n; i++)
     {
        int k = Angle_Sort.index[i];
        x[i] = Angle_Sort.x[k];
        y[i] = Angle_Sort.y[k];
        /* fprintf (stderr, "%f %f %f\n", x[i], y[i], Angle_Sort.angles[k] * 180.0/M_PI); */
     }

   return 0;
}

/* Choose the cyclic permutation of the 4 (x,y) points that maximizes
 * the sum of edge-vector dot-products.  The goal is to make the sequence
 * of points in the final polygon match the sequence of point in the original
 * polygon.  In the limit that the vertices are only slightly perturbed, this
 * should work robustly, but in the general case, this problem may be ill-posed.
 * But let's try, because inconsistent vertex ordering definitely causes problems
 * later on when trying bin up these pixels.
 */
static int select_permutation (double *x, double *y,
                               const double *x0, const double *y0, int n)
{
   double max_sum, sum_p, xs[4], ys[4];
   int i, in, k, kn, p, pmax;
   /* The input points are in CCW order, and all cyclic permutations
    * preserve this order so, by construction, the output is
    * guaranteed to be in CCW order */

   if (n != 4)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: expected polygon with 4 vertices!! num_vertices=%d", __func__, n);
        return -1;
     }

   pmax = 0;
   max_sum = 0.0;

   for (p = 0; p < 4; p++)
     {
        sum_p = 0.0;
        k = p;
        for (i = 0; i < 4; i++)
          {
             in = (i+1) % 4;
             kn = (k+1) % 4;
             sum_p += ((x0[in] - x0[i]) * (x[kn] - x[k])
                       + (y0[in] - y0[i]) * (y[kn] - y[k]));
             k = kn;
          }

        if (sum_p > max_sum)
          {
             pmax = p;
             max_sum = sum_p;
          }
     }

   /* Permutation that maximizes the sum of edge-vector dot-products */
   k = pmax;

   for (i = 0; i < 4; i++)
     {
        xs[i] = x[k];
        ys[i] = y[k];
        k = (k+1) % 4;
     }
   for (i = 0; i < 4; i++)
     {
        x[i] = xs[i];
        y[i] = ys[i];
     }

   return 0;
}

static int apply_parallax_adj (const Geoid_Data_Type *gdt, const ECEF_Vector sat,
                               double *lon, double *lat, int n, double *shift)
{
   /* Convergence criterion on corrected height [km] */
   double height_tol_km = 10.0e-3;
   /* Maximum number of iterations to converge */
   int maxiter = 50;
   int i;

   for (i = 0; i < n; i++)
     {
        double lat_tmp, lon_tmp;

        if (shift)
          {
             shift[i] = 0.0;
          }

        if (invalid_lonlat (lon[i], lat[i]))
          continue;

        if (TEMPO_GEO_NO_ERR != parallaxAdj (lat[i], lon[i],
                                             sat, &gdt->geoid_height, &gdt->dem, height_tol_km, maxiter,
                                             &lat_tmp, &lon_tmp))
          {
             tell_vlog (TELL_MSGTYPE_INFO, 2, "%s: parallaxAdj failed: lat=%f lon=%f",
                        __func__, lat[i], lon[i]);
             continue;
          }

        if (shift)
          {
             shift[i] = compute_shift (lat[i], lon[i], lat_tmp, lon_tmp);
          }

        lat[i] = lat_tmp;
        lon[i] = lon_tmp;
     }

   return 0;
}

static int correct_band_geolocation_for_parallax (const Geoid_Data_Type *gdt,
                                                  const ECEF_Position_Type *satloc,
                                                  Geoloc_Type *geoloc)
{
   unsigned int s, num_xtrack_corners;
   double *lon_step_cnr0 = NULL;
   double *lat_step_cnr0 = NULL;
   int status = -1;

   num_xtrack_corners = geoloc->num_xtrack * geoloc->num_corner;

   if ((NULL == (lon_step_cnr0 = (double *) MALLOC (num_xtrack_corners * sizeof(double))))
       || (NULL == (lat_step_cnr0 = (double *) MALLOC (num_xtrack_corners * sizeof(double)))))
     goto return_status;

   for (s = 0; s < geoloc->num_mirror_step; s++)
     {
        ECEF_Vector sat;
        double *lon_step = geoloc->lon + s * geoloc->num_xtrack;
        double *lat_step = geoloc->lat + s * geoloc->num_xtrack;
        double *lon_step_cnr = geoloc->lon_cnr + s * num_xtrack_corners;
        double *lat_step_cnr = geoloc->lat_cnr + s * num_xtrack_corners;
        double *shift_km_step = NULL;
        double *lat_cnr = NULL;
        double *lon_cnr = NULL;
        unsigned int y, c, num_valid_points;

        /* keep a copy of the uncorrected corners */
        memcpy ((char *)lon_step_cnr0, (char *)lon_step_cnr, num_xtrack_corners * sizeof(double));
        memcpy ((char *)lat_step_cnr0, (char *)lat_step_cnr, num_xtrack_corners * sizeof(double));

        sat.theX = satloc->X[s];
        sat.theY = satloc->Y[s];
        sat.theZ = satloc->Z[s];

#ifdef OUTPUT_PARALLAX_SHIFT
        shift_km_step = geoloc->shift_km + s * geoloc->num_xtrack;
#endif

        /* Adjust pixel centers */
        if (0 != apply_parallax_adj (gdt, sat, lon_step, lat_step, geoloc->num_xtrack, shift_km_step))
          return -1;

        /* Adjust pixel corners */
        if (0 != apply_parallax_adj (gdt, sat, lon_step_cnr, lat_step_cnr, num_xtrack_corners, NULL))
          return -1;

        /* After parallax adjustment, some pixel polygons may be non-convex.
         * Sort the vertices to obtain a convex polygon, and if necessary/possible,
         * permute the vertices to restore the original vertex sequence.  The
         * intent is to keep the relative location of pixel 0 consistent across all
         * pixels.  This simplifies the process of determining the corners of binned
         * pixels later on.
         */
        for (y = 0; y < geoloc->num_xtrack; y++)
          {
             double *lat_cnr0 = lat_step_cnr0 + y * geoloc->num_corner;
             double *lon_cnr0 = lon_step_cnr0 + y * geoloc->num_corner;
             lat_cnr = lat_step_cnr + y * geoloc->num_corner;
             lon_cnr = lon_step_cnr + y * geoloc->num_corner;

             num_valid_points = 0;
             for (c = 0; c < geoloc->num_corner; c++)
               {
                  if (invalid_lonlat (lon_cnr[c], lat_cnr[c]))
                    break;
                  num_valid_points++;
               }

             if (num_valid_points == geoloc->num_corner)
               {
                  (void) polygon_sort_ccw (lon_cnr, lat_cnr, geoloc->num_corner, NULL);
                  (void) select_permutation (lon_cnr, lat_cnr, lon_cnr0, lat_cnr0, geoloc->num_corner);
               }
          }
     }

   status = 0;
return_status:
   FREE(lon_step_cnr0);
   FREE(lat_step_cnr0);

   return status;
}

static char *expand_path (const char *path)
{
   wordexp_t we = {0};
   char *s = NULL;

   if ((0 != wordexp (path, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, path);
        goto return_status;
     }

   if (NULL == (s = strdup (we.we_wordv[0])))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

return_status:
   wordfree (&we);
   return s;
}

static int have_terrain_referenced_coordinates (const Granule_Type *gt, const char *band_name)
{
   const char *attname = ATTR_STRING_TERRAIN_REFERENCED;
   int idp, grp;

   if (0 != TIO_inq_grp (gt->ncid, band_name, &grp))
     return -1;

   /* If the attribute is present, the file coordinates are terrain-referenced */
   if (NC_NOERR == nc_inq_attid (grp, NC_GLOBAL, attname, &idp))
     return 1;

   return 0;
}

/* INR lon-lat coordinates are referenced to the WGS84 ellipsoid. This means that features
 * that sit at a significant vertical offset from the ellipsoid (e.g. mountain-tops, cloud-tops)
 * are assigned the coordinates of the point where the line of sight from the spacecraft
 * through the feature intersects the ellipsoid.  To derive the actual terrain-referenced
 * coordinates of the feature, it is necessary to correct for the offset associated with
 * projection from the satellite viewpoint onto the ellipsoid.  The INR software refers
 * to this as a "parallax" correction. The function parallaxAdj is provided to correct for it.
 */
static int correct_geolocation_for_parallax (Granule_Type *gt, const char **band_names,
                                             TIO_Meta_Type *meta, config_t *cfg)
{
   config_setting_t *s;
   Geoid_Data_Type gdt = {0};
   const char *geoid_dem_setting;
   char *geoid_dem_path = NULL;
   char *geoid_dem_basename;
   int num_to_correct, have_corrected[NUM_BANDS];
   int i, status = -1;

   num_to_correct = NUM_BANDS;
   for (i = 0; i < NUM_BANDS; i++)
     {
        if ((have_corrected[i] = have_terrain_referenced_coordinates (gt, band_names[i])) < 0)
          return -1;
        if (have_corrected[i]) num_to_correct--;
     }

   if (num_to_correct == 0)
     return 0;

   if ((NULL == (s = config_lookup (cfg, "parallax_correction")))
       || (CONFIG_TRUE != config_setting_lookup_string (s, "altitude_file", &geoid_dem_setting)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading parallax_correction.altitude_file", __func__);
        return -1;
     }

   if (NULL == (geoid_dem_path = expand_path (geoid_dem_setting)))
     return -1;

   if (0 != read_geoid_data (&gdt, geoid_dem_path))
     goto free_and_return;

   for (i = 0; i < NUM_BANDS; i++)
     {
        if (0 == have_corrected[i])
          {
             if (0 != correct_band_geolocation_for_parallax (&gdt, &gt->sat, &gt->geoloc[i]))
               goto free_and_return;
             if (0 != write_band_geolocation (gt, band_names[i], &gt->geoloc[i], 1))
               goto free_and_return;
          }
     }

   if (NULL != (geoid_dem_basename = strrchr (geoid_dem_path, '/')))
     {
        geoid_dem_basename++;
     }
   else geoid_dem_basename = geoid_dem_path;
   if (0) tio_meta_append_string (meta, "input_files", geoid_dem_basename);

   status = 0;
free_and_return:
   free_geoid_data (&gdt);
   FREE(geoid_dem_path);
   return status;
}

typedef struct
{
   projPJ geos;
   projPJ longlat;
}
Proj_Type;

static void free_proj (Proj_Type *p)
{
   if (NULL == p)
     return;
   pj_free (p->geos);
   pj_free (p->longlat);
}

static int init_proj (Proj_Type *p, double sat_lon)
{
   char ctl_geos[PROJ_ARGS_BUFSIZE];
   const char geos_fmt[] = "+proj=geos +lon_0=%0.3g +h=%0.1f";
   const char ctl_longlat[] = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
   int len;

   if (NULL == (p->longlat = pj_init_plus (ctl_longlat)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus(longlat) failed", __func__);
        return -1;
     }

   memset (ctl_geos, 0, PROJ_ARGS_BUFSIZE);
   len = snprintf (ctl_geos, PROJ_ARGS_BUFSIZE, geos_fmt, sat_lon, GEO_ALTITUDE);
   if (len >= PROJ_ARGS_BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: proj4 arg buffer too small", __func__);
        return -1;
     }

   if (NULL == (p->geos = pj_init_plus (ctl_geos)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus(geos) failed", __func__);
        return -1;
     }

   return 0;
}

static int impose_band_pixel_corner_sequence (Proj_Type *proj, Geoloc_Type *geoloc)
{
   unsigned int s, c, i, j, k, num_xtrack_corners, num_valid_points;
   double lon_cpy[4], lat_cpy[4];
   double *xx = NULL;
   double *yy = NULL;
   double degtorad = M_PI/180.0;
   int num_points_processed = 0;
   int status = -1;

   num_xtrack_corners = geoloc->num_xtrack * geoloc->num_corner;

   if ((NULL == (xx = (double *) MALLOC (num_xtrack_corners * sizeof(double))))
       || (NULL == (yy = (double *) MALLOC (num_xtrack_corners * sizeof(double)))))
     goto return_status;

   for (s = 0; s < geoloc->num_mirror_step; s++)
     {
        double *lon_step_cnr = geoloc->lon_cnr + s * num_xtrack_corners;
        double *lat_step_cnr = geoloc->lat_cnr + s * num_xtrack_corners;

        /* make a working copy of the corner coordinates */
        memcpy ((char *)xx, (char *)lon_step_cnr, num_xtrack_corners*sizeof(double));
        memcpy ((char *)yy, (char *)lat_step_cnr, num_xtrack_corners*sizeof(double));

        /* Sort the pixel corners in a geostationary projection where
         * the grid is nearly Cartesian with very little curvature.
         * An exact Cartesian grid would be best, but this code would then
         * depend on the exact boresight aim point, and I don't think achieving
         * a perfect projection is that critical here. Even if we the projection
         * were perfect, there would still be some randomness in pixel orientations
         * because the spacecraft attitude and stepping are imperfect.
         * These small differences should not be enough to rotate or distort
         * any of the pixels enough to matter.
         */
        for (i = 0; i < num_xtrack_corners; i++)
          {
             xx[i] *= degtorad;
             yy[i] *= degtorad;
          }
        if ((status = pj_transform (proj->longlat, proj->geos, num_xtrack_corners, 1, xx, yy, NULL)) != 0)
          {
             tell_verror (TELL_APPLICATION_ERROR,
                          "%s: pj_transform failed, status = %d (%s)",
                          __func__, status, pj_strerrno(status));
             status = -1;
             goto return_status;
          }

        for (k = 0; k < geoloc->num_xtrack; k++)
          {
             double *lon_cnr = lon_step_cnr + k * geoloc->num_corner;
             double *lat_cnr = lat_step_cnr + k * geoloc->num_corner;
             double *x = xx + k * geoloc->num_corner;
             double *y = yy + k * geoloc->num_corner;

             num_valid_points = 0;
             for (c = 0; c < geoloc->num_corner; c++)
               {
                  if ((0 != invalid_lonlat (lon_cnr[c], lat_cnr[c]))
                      || (0 != invalid_point (x[c], y[c])))
                    break;
                  num_valid_points++;
               }

             /* Impose a uniform CCW ordering by sorting the geostationary
              * projection of each pixel's corners, then use those sort indices
              * to re-order the lon,lat corner coordinates.
              */
             if (num_valid_points == geoloc->num_corner)
               {
                  int *indices;
                  (void) polygon_sort_ccw (x, y, geoloc->num_corner, &indices);
                  for (j = 0; j < geoloc->num_corner; j++)
                    {
                       int p = indices[j];
                       lon_cpy[j] = lon_cnr[p];
                       lat_cpy[j] = lat_cnr[p];
                    }
                  for (j = 0; j < geoloc->num_corner; j++)
                    {
                       lon_cnr[j] = lon_cpy[j];
                       lat_cnr[j] = lat_cpy[j];
                    }
                  num_points_processed++;
               }
          }
     }

   status = 0;
return_status:
   if (num_points_processed == 0)
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0, "%s: no valid pixel corner points?", __func__);
     }
   FREE(xx);
   FREE(yy);
   return status;
}

static int impose_pixel_corner_sequence (Granule_Type *gt, const char **band_names)
{
   Proj_Type proj = {0};
   const ECEF_Position_Type *sat = &gt->sat;
   double sat_lon, X_sum, Y_sum;
   unsigned int i;
   int status = -1;

   /* Satellite longitude is required to define the geostationary projection */
   X_sum = 0.0;
   Y_sum = 0.0;
   for (i = 0; i < sat->num; i++)
     {
        X_sum += sat->X[i];
        Y_sum += sat->Y[i];
     }
   sat_lon = atan2 (Y_sum/sat->num, X_sum/sat->num) * 180.0/M_PI;

   if (0 != init_proj (&proj, sat_lon))
     goto free_and_return;

   for (i = 0; i < NUM_BANDS; i++)
     {
        if (0 != impose_band_pixel_corner_sequence (&proj, &gt->geoloc[i]))
          goto free_and_return;
        if (0 != write_band_geolocation (gt, band_names[i], &gt->geoloc[i], 0))
          goto free_and_return;
     }

   status = 0;
free_and_return:
   free_proj (&proj);
   return status;
}

static int read_band_geolocation (Granule_Type *gt, const char *band_name,
                                  Geoloc_Type *geoloc, int *pnum_loc)
{
   TIO_Var_Info_Type info;
   int i, grp, start[3], count[3];
   int num_pixels, num_loc;

   if (0 != TIO_inq_grp (gt->ncid, band_name, &grp))
     return -1;

   if (0 != TIO_inq_var (grp, TEMPO_VAR_LONGITUDE_BOUNDS, &info))
     return -1;

   num_pixels = info.dimlens[0] * info.dimlens[1];

   if (0 != alloc_geoloc_fields (info.dimlens, geoloc))
     return -1;

   geoloc->group = grp;

   for (i = 0; i < 3; i++)
     {
        geoloc->dimids[i] = info.dimids[i];
        start[i] = 0;
        count[i] = info.dimlens[i];
     }

   if ((0 != TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE, start, count, TIO_DOUBLE, geoloc->lon))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_LATITUDE, start, count, TIO_DOUBLE, geoloc->lat))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_INRQF, start, count, TIO_INT, geoloc->inr_quality_flag))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS, start, count, TIO_DOUBLE, geoloc->lon_cnr))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS, start, count, TIO_DOUBLE, geoloc->lat_cnr)))
     {
        goto return_error;
     }

   num_loc = 0;
   for (i = 0; i < num_pixels; i++)
     {
        if (geoloc->inr_quality_flag[i] == 0)
          num_loc++;
     }
   tell_vlog (TELL_MSGTYPE_INFO, 1, "band %s: got %d geolocated pixels", band_name, num_loc);
   if (pnum_loc) *pnum_loc = num_loc;

   return 0;

return_error:
   free_geoloc_fields (geoloc);
   return -1;
}

static int read_geolocation (Granule_Type *gt, const char **band_names)
{
   int i, num_loc_tot = 0;

   for (i = 0; i < NUM_BANDS; i++)
     {
        int num_loc = 0;
        if (0 != read_band_geolocation (gt, band_names[i], &gt->geoloc[i], &num_loc))
          return -1;
        num_loc_tot += num_loc;
     }

   if (num_loc_tot == 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: no pixels have valid geolocation", __func__);
        return -1;
     }

   return 0;
}

static int read_ecef_position (int grp, size_t num, const char **varnames,
                               ECEF_Position_Type *p)
{
   int start, count;

   if (0 != alloc_ecef_position (num, p))
     return -1;

   start = 0;
   count = num;

   if ((0 != TIO_get_var_section (grp, varnames[0], &start, &count, NC_DOUBLE, p->X))
       || (0 != TIO_get_var_section (grp, varnames[1], &start, &count, NC_DOUBLE, p->Y))
       || (0 != TIO_get_var_section (grp, varnames[2], &start, &count, NC_DOUBLE, p->Z)))
     {
        free_ecef_position (p);
        return -1;
     }

   return 0;
}

static int read_ecef_geometry (Granule_Type *gt)
{
   TIO_Var_Info_Type info;
   const char *sun_vars[] =
     {TEMPO_VAR_SUN_X, TEMPO_VAR_SUN_Y, TEMPO_VAR_SUN_Z};
   const char *moon_vars[] =
     {TEMPO_VAR_MOON_X, TEMPO_VAR_MOON_Y, TEMPO_VAR_MOON_Z};
   const char *sat_vars[] =
     {TEMPO_VAR_SAT_X, TEMPO_VAR_SAT_Y, TEMPO_VAR_SAT_Z};
   int grp;

   if (0 != TIO_inq_grp (gt->ncid, TEMPO_GRP_GEOMETRY, &grp))
     return -1;
   if (0 != TIO_inq_var (grp, TEMPO_VAR_SUN_X, &info))
     return -1;
   if (0 != read_ecef_position (grp, info.dimlens[0], sun_vars, &gt->sun))
     return -1;
   if (0 != read_ecef_position (grp, info.dimlens[0], moon_vars, &gt->moon))
     return -1;
   if (0 != read_ecef_position (grp, info.dimlens[0], sat_vars, &gt->sat))
     return -1;

   return 0;
}

Granule_Type *granule_open (const char *file, int correct_parallax,
                            TIO_Meta_Type *meta, config_t *cfg)
{
   Granule_Type *gt = NULL;
   const char *band_names[] = {TEMPO_BAND_NAME_UV, TEMPO_BAND_NAME_VIS};

   if (NULL == (gt = new_granule_type ()))
     return NULL;

   if (0 != TIO_open (file, NC_WRITE, &gt->ncid))
     {
        gt->gt_close (gt);
        return NULL;
     }

   if (0 != tio_history_append_cmdline (gt->ncid))
     {
        gt->gt_close (gt);
        return NULL;
     }

   if (0 != tio_meta_ncinit (meta, gt->ncid, "input_files", TIO_META_TYPE_STRING))
     {
        gt->gt_close (gt);
        return NULL;
     }

   if (0 != read_geolocation (gt, band_names))
     {
        gt->gt_close (gt);
        return NULL;
     }

   if (0 != read_ecef_geometry (gt))
     {
        gt->gt_close (gt);
        return NULL;
     }

   /* Ideally, INR should deliver pixel corners with the correct
    * sort order and a uniform sequence, but we sort them anyway.
    */
   if (0 != impose_pixel_corner_sequence (gt, band_names))
     {
        gt->gt_close (gt);
        return NULL;
     }

   /* If we correct lon-lat coordinates for local height above/below
    * the WGS84 ellipsoid, we must do that before doing any other
    * calculations that use the lon-lat coordinates.
    */
   if (correct_parallax)
     {
        if (0 != correct_geolocation_for_parallax (gt, band_names, meta, cfg))
          {
             gt->gt_close (gt);
             return NULL;
          }
     }
   else tell_vlog (TELL_MSGTYPE_INFO, 0, "Skipping parallax correction");

   return gt;
}
