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

/* To compute the parallax shift distance and write it to the output file: */
/* #define OUTPUT_PARALLAX_SHIFT 1 */
#undef OUTPUT_PARALLAX_SHIFT

#define ALWAYS_SORT_POLYGONS 1

typedef struct
{
   double *lon;
   double *lat;
   double *lon_cnr;
   double *lat_cnr;
   double *shift_km;
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
   return ((fabs(lon) > 360.0) || (fabs(lat) > 90.0));
}

static void free_geoloc_fields (Geoloc_Type *geoloc)
{
   if (geoloc == NULL)
     return;
   FREE(geoloc->lon);
   FREE(geoloc->lat);
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

   if ((NULL == (geoloc->lon = alloc_geoloc_coordinate (geoloc)))
       || (NULL == (geoloc->lat = alloc_geoloc_coordinate (geoloc))))
     {
        free_geoloc_fields (geoloc);
        return -1;
     }

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
        {"coordinates", "longitude latitude"},
        {NULL, NULL}
     };
   TIO_Attr_Text_Type dev_text_attrs[] =
     {
        {"units", "m"},
        {"long_name", TEMPO_VAR_TERR_HEIGHT_STDDEV},
        {"comment", "Area-weighted terrain height standard deviation within pixel boundary"},
        {"coordinates", "longitude latitude"},
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
        {"coordinates", "longitude latitude"},
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
                                       start, count, TIO_FILL_FLOAT, snow_ice_fraction))
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

             BITMASK_CLEAR(flags, BITMASK_GPQF_BITS_USED);

             BITMASK_SET(flags,   lc_typeqc[j] >> 4);
             BITMASK_SET(flags, illum_flags[j] << 4);
             BITMASK_SET(flags, (inr_quality_flag[j] & 0x01) << 6);
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
   double dist_km = sqrt (x*x + y*y + z+z);
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

static int set_object_angles (Granule_Type *gt)
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

   return gt;
}

#ifdef OUTPUT_PARALLAX_SHIFT
static int write_shift (const Granule_Type *gt, const char *band_name,
                        const Geoloc_Type *geoloc)
{
   const char *varname = "parallax_shift";
   static TIO_Attr_Text_Type attrs[] =
     {
        {"coordinates", "longitude latitude"},
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
                                   const Geoloc_Type *geoloc)
{
   const char *attname = "terrain_referenced_coordinates";
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

   if (0 != TIO_put_att (grp, NC_GLOBAL, attname, NC_CHAR, 1+strlen(yes), yes))
     return -1;

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

#ifndef ALWAYS_SORT_POLYGONS
static int isconvex_polygon (const double *x, const double *y, int n)
{
   double wsign=0;
   int xsign=0, xsignfirst=0, xflips=0;
   int ysign=0, ysignfirst=0, yflips=0;
   int i, prev, curr, next;

   if (n < 3)
     return 0;

   curr = n-2;
   next = n-1;

   for (i = 0; i < n; i++)
     {
        double w, ax, ay, bx, by;

        prev = curr;
        curr = next;
        next = i;

        /* previous edge vector ("before") */
        bx = x[curr] - x[prev];
        by = y[curr] - y[prev];

        /* next edge vector ("after") */
        ax = x[next] - x[curr];
        ay = y[next] - y[curr];

        /* calculate sign flips using next edge vector ("after"),
         * recording the first sign */
        if (ax > 0)
          {
             if (xsign == 0)
               xsignfirst = +1;
             else if (xsign < 0)
               xflips++;
             xsign = +1;
          }
        else if (ax < 0)
          {
             if (xsign == 0)
               xsignfirst = -1;
             else if (xsign > 0)
               xflips++;
             xsign = -1;
          }

        if (xflips > 2)
          return 0;

        if (ay > 0)
          {
             if (ysign == 0)
               ysignfirst = +1;
             else if (ysign < 0)
               yflips++;
             ysign = +1;
          }
        else if (ay < 0)
          {
             if (ysign == 0)
               ysignfirst = -1;
             else if (ysign > 0)
               yflips++;
             ysign = -1;
          }

        if (yflips > 2)
          return 0;

        /* Find out the orientation of this pair of edges,
         * and ensure it does not differ from previous pairs. */
        w = bx * ay - ax * by;
        if ((wsign == 0) && (w != 0))
          wsign = w;
        else if ((wsign > 0) && (w < 0))
          return 0;
        else if ((wsign < 0) && (w > 0))
          return 0;
     }

   /* final/wraparound sign flips */
   if ((xsign != 0) && (xsignfirst != 0) && (xsign != xsignfirst))
     xflips++;
   if ((ysign != 0) && (ysignfirst != 0) && (ysign != ysignfirst))
     yflips++;

   /* convex polygons have two sign flips along each axis */
   if ((xflips != 2) || (yflips != 2))
     return 0;

   /* This is a convex polygon */
   return 1;
}
#endif

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

static int polygon_sort_ccw (double *x, double *y, int n)
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

   qsort (Angle_Sort.index, sizeof(int), (size_t)n, compare_angle_indices);

   for (i = 0; i < n; i++)
     {
        int k = Angle_Sort.index[i];
        x[i] = Angle_Sort.x[k];
        y[i] = Angle_Sort.y[k];
        /* fprintf (stderr, "%f %f %f\n", x[i], y[i], Angle_Sort.angles[k] * 180.0/M_PI); */
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

   num_xtrack_corners = geoloc->num_xtrack * geoloc->num_corner;

   for (s = 0; s < geoloc->num_mirror_step; s++)
     {
        ECEF_Vector sat;
        double *lon_step = geoloc->lon + s * geoloc->num_xtrack;
        double *lat_step = geoloc->lat + s * geoloc->num_xtrack;
        double *lon_step_cnr = geoloc->lon_cnr + s * num_xtrack_corners;
        double *lat_step_cnr = geoloc->lat_cnr + s * num_xtrack_corners;
        double *shift_km_step = NULL;
        unsigned int y, c;

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
         * When necessary, re-order the vertices to obtain a convex polygon.
         */
        for (y = 0; y < geoloc->num_xtrack; y++)
          {
             double *lat_cnr = lat_step_cnr + y * geoloc->num_corner;
             double *lon_cnr = lon_step_cnr + y * geoloc->num_corner;
             unsigned int num_valid_points = 0;

             for (c = 0; c < geoloc->num_corner; c++)
               {
                  if (invalid_lonlat (lon_cnr[c], lat_cnr[c]))
                    break;
                  num_valid_points++;
               }

             if (num_valid_points == geoloc->num_corner)
               {
#ifdef ALWAYS_SORT_POLYGONS
                  /* Just always sort CCW.  It's slower, but catches the corner case
                   * when the points end up in CW order.
                   */
                  (void) polygon_sort_ccw (lon_cnr, lat_cnr, geoloc->num_corner);
#else
                  if (0 == isconvex_polygon (lon_cnr, lat_cnr, geoloc->num_corner))
                    {
                       tell_vlog (TELL_MSGTYPE_INFO, 2, "%s: non-convex polygon: mirror_step=%d xtrack=%d",
                             __func__, s, y);
                       (void) polygon_sort_ccw (lon_cnr, lat_cnr, geoloc->num_corner);
                    }
#endif
               }
          }
     }

   return 0;
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
   const char *attname = "terrain_referenced_coordinates";
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
static int correct_geolocation_for_parallax (Granule_Type *gt, TIO_Meta_Type *meta, config_t *cfg)
{
   config_setting_t *s;
   Geoid_Data_Type gdt = {0};
   const char *band_names[] = {TEMPO_BAND_NAME_UV, TEMPO_BAND_NAME_VIS};
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
             if (0 != write_band_geolocation (gt, band_names[i], &gt->geoloc[i]))
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

static int read_band_geolocation (Granule_Type *gt, const char *band_name,
                                  Geoloc_Type *geoloc)
{
   TIO_Var_Info_Type info;
   int i, grp, start[3], count[3];

   if (0 != TIO_inq_grp (gt->ncid, band_name, &grp))
     return -1;

   if (0 != TIO_inq_var (grp, TEMPO_VAR_LONGITUDE_BOUNDS, &info))
     return -1;

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
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_LONGITUDE_BOUNDS, start, count, TIO_DOUBLE, geoloc->lon_cnr))
       || (0 != TIO_get_var_section (grp, TEMPO_VAR_LATITUDE_BOUNDS, start, count, TIO_DOUBLE, geoloc->lat_cnr)))
     {
        goto return_error;
     }

   return 0;

return_error:
   free_geoloc_fields (geoloc);
   return -1;
}

static int read_geolocation (Granule_Type *gt)
{
   const char *band_names[] = {TEMPO_BAND_NAME_UV, TEMPO_BAND_NAME_VIS};
   int i;

   for (i = 0; i < NUM_BANDS; i++)
     {
        if (0 != read_band_geolocation (gt, band_names[i], &gt->geoloc[i]))
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

   if (0 != read_geolocation (gt))
     {
        gt->gt_close (gt);
        return NULL;
     }

   if (0 != read_ecef_geometry (gt))
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
        if (0 != correct_geolocation_for_parallax (gt, meta, cfg))
          {
             gt->gt_close (gt);
             return NULL;
          }
     }
   else tell_vlog (TELL_MSGTYPE_INFO, 0, "Skipping parallax correction");

   return gt;
}
