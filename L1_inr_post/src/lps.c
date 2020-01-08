/** @file lps.c
 *  @brief Interface for linear polarization sensitivity
 *         lookup table
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <wordexp.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include <libconfig.h>
#include <proj_api.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "lps.h"

#define DEGTORAD      (M_PI/180.0)
#define EARTH_RADIUS   6378137.0   /* meters */
#define GEO_ALTITUDE  35785831.0   /* meters */
#define GEO_RADIUS    (EARTH_RADIUS + GEO_ALTITUDE)

#define PROJ_ARGS_BUFSIZE       80

typedef struct
{
   double x;
   double y;
   double z;
}
Vector_Type;

typedef struct
{
   double sat_lon;
   /**< GEO satellite orbital station [radians] */

   double tilt;
   /**< tilt > 0 is northward tilt of instrument boresight
    *   about spacecraft roll axis [radians] */

   double azi;
   /**< azi > 0 is rotation, eastward from north (CW),
    *   about instrument boresight axis [radians] */

   Vector_Type spacecraft_u;
   /**< unit vector pointing from Earth's center toward spacecraft */

   Vector_Type irp_u;
   /**< unit vector normal to instrument reference plane (IRP) */
}
Instr_Geom_Type;

typedef struct
{
   double *lpsens;   /**< linear polarization sensivity [num_xtrack, num_wave] */
   double *angmax;   /**< angle of maximum transmission [num_xtrack, num_wave] [rad] */
   double *wave;     /**< wavelength grid [num_wave] */
   int num_xtrack;
   int num_wave;
}
Lps_Table_Type;

struct Lps_Type
{
   Instr_Geom_Type geom;
   double *mirror_x;   /**< azimuth angle from boresight, EAST is positive [rad] */
   int num_mirror_x;
   Lps_Table_Type *table_uv;    /**< UV band [num_mirror_x] */
   Lps_Table_Type *table_vis;   /**< VIS band [num_mirror_x] */

   projPJ tpers;
   projPJ longlat;

   gsl_interp *interp;
   gsl_interp_accel *acc;
   double *tmp;                  /* reusable scratch space */
   int num_tmp;
   size_t size_tmp;
};

static void vec_unit (double theta, double phi, Vector_Type *v)
{
   double sin_t = sin(theta);
   v->x = sin_t * cos(phi);
   v->y = sin_t * sin(phi);
   v->z = cos(theta);
}

static void vec_scale (double c, const Vector_Type *a, Vector_Type *b)
{
   b->x = a->x * c;
   b->y = a->y * c;
   b->z = a->z * c;
}

static void vec_cross (const Vector_Type *a, const Vector_Type *b, Vector_Type *c)
{
   c->x =   a->y * b->z - a->z * b->y;
   c->y = - a->x * b->z + a->z * b->x;
   c->z =   a->x * b->y - a->y * b->x;
}

static void vec_diff (const Vector_Type *a, const Vector_Type *b, Vector_Type *c)
{
   c->x = a->x - b->x;
   c->y = a->y - b->y;
   c->z = a->z - b->z;
}

static double vec_dot (const Vector_Type *a, const Vector_Type *b)
{
   return a->x * b->x + a->y * b->y + a->z * b->z;
}

static double vec_length (const Vector_Type *a)
{
   return sqrt(vec_dot(a,a));
}

static int vec_norm (Vector_Type *p)
{
   double len = vec_length (p);
   if (len == 0.0) return -1;
   p->x /= len;
   p->y /= len;
   p->z /= len;
   return 0;
}

/* Rotate 3D vector v about axis k by an angle theta, yielding vector vrot,
 * where theta > 0 is defined according to the right hand rule.
 * Rodrigues' rotation formula:
 *   vrot = v cos(t) + (k \cross v) sin(t) + k (k \dot v) (1 - cos(t))
 */
static void vec_rotate (const Vector_Type *v, const Vector_Type *k,
                        double theta, Vector_Type *vrot)
{
   Vector_Type cross, n, p;
   double dot, cos_t = cos(theta);

   vec_cross (k, v, &cross);
   vec_scale (sin(theta), &cross, &n);

   dot = vec_dot(k, v);
   vec_scale (dot * (1 - cos_t), k, &p);

   vec_scale (cos_t, v, vrot);
   vrot->x += n.x + p.x;
   vrot->y += n.y + p.y;
   vrot->z += n.z + p.z;
}

static int read_geom (config_t *cfg, Instr_Geom_Type *geom)
{
   config_setting_t *s;

   if (NULL == (s = config_lookup (cfg, "polcorr_instr_geom")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing instr_geom in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "tilt", &geom->tilt))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "azi", &geom->azi)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing %s in param file: %s",
                     __func__, config_setting_name(s),
                     config_setting_source_file (s));
        return -1;
     }

   geom->tilt *= DEGTORAD;
   geom->azi *= DEGTORAD;

   return 0;
}

static void lps_proj_close (Lps_Type *lps)
{
   if (lps == NULL)
     return;
   pj_free (lps->longlat);
   pj_free (lps->tpers);
}

static int lps_proj_open (Lps_Type *lps)
{
   char ctl_tpers[PROJ_ARGS_BUFSIZE];
   const char tpers_fmt[] =
     "+proj=tpers +lat_0=0 +lon_0=%0.3g +h=%0.1f +tilt=%0.3g +azi=%0.3g";
   const char ctl_longlat[] =
     "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
   Instr_Geom_Type *geom;
   int len;

   if (lps == NULL)
     return -1;

   geom = &lps->geom;

   if (NULL == (lps->longlat = pj_init_plus (ctl_longlat)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus(longlat) failed", __func__);
        return -1;
     }

   memset (ctl_tpers, 0, PROJ_ARGS_BUFSIZE);
   len = snprintf (ctl_tpers, PROJ_ARGS_BUFSIZE, tpers_fmt,
                   geom->sat_lon / DEGTORAD, GEO_ALTITUDE,
                   geom->tilt / DEGTORAD,
                   geom->azi / DEGTORAD);
   if (len >= PROJ_ARGS_BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: proj4 arg buffer too small", __func__);
        return -1;
     }

   if (NULL == (lps->tpers = pj_init_plus (ctl_tpers)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: pj_init_plus(tpers) failed", __func__);
        return -1;
     }

   return 0;
}

static int lonlat_to_mirror_xy (Lps_Type *lps, double lon, double lat,
                                double *mirror_x, double *mirror_y)
{
   double x, y, m_x, m_y, cos_a, sin_a;
   int status, n=1;

   /* input (lon,lat) in radians */
   x = lon;
   y = lat;

   if ((status = pj_transform (lps->longlat, lps->tpers, n, 1, &x, &y, NULL)) != 0)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return -1;
     }

   /* convert tilted perspective (x,y) meters to radians */
   x /= GEO_ALTITUDE;
   y /= GEO_ALTITUDE;

   /* (x,y) are angular coordinates in the tilted perspective
    * plane (Proj4 'tpers')
    * x [radian] = azimuth angle
    * y [radian] = altitude angle
    * (x,y) origin is at GEO satellite foot point on equator,
    * with +x east, +y north
    */

   /* To convert angular tpers coordinates to mirror coordinates,
    * offset to FOR center (the boresight), and account for azi rotation
    */
   y -= lps->geom.tilt;
   cos_a = cos(-lps->geom.azi);
   sin_a = sin(-lps->geom.azi);
   m_x = x * cos_a - y * sin_a;
   m_y = x * sin_a + y * cos_a;

   /* SMA mirror coordinates have +X "east" and +Y "SOUTH",
    * so that:
    *        mirror_x = m_x
    *        mirror_y = -m_y
    */
   m_y *= -1;

   if (0) fprintf (stderr, "lon=%g lat=%g => mirror_x = %g\n",
                   lon/DEGTORAD, lat/DEGTORAD, m_x);

   if (mirror_x) *mirror_x = m_x;
   if (mirror_y) *mirror_y = m_y;

   return 0;
}

static int mirror_xy_to_lonlat (Lps_Type *lps, double *lon, double *lat,
                                double mirror_x, double mirror_y)
{
   double x, y, m_x, m_y, cos_a, sin_a;
   int status, n=1;

   /* mirror_x,mirror_y coordinates are offset from tilted perspective
    * angular coordinates by a Y sign flip, a rotation, and an offset:
    */

   m_x = mirror_x;
   m_y = -mirror_y;

   cos_a = cos(-lps->geom.azi);
   sin_a = sin(-lps->geom.azi);
   x =  m_x * cos_a + m_y * sin_a;
   y = -m_x * sin_a + m_y * cos_a;

   y += lps->geom.tilt;

   /* tilted perspective coordinates in meters */
   x *= GEO_ALTITUDE;
   y *= GEO_ALTITUDE;

   /* transform to (lon,lat) */
   if ((status = pj_transform (lps->tpers, lps->longlat, n, 1, &x, &y, NULL)) != 0)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return -1;
     }

   /* return (lon,lat) in radians */
   if (lon) *lon = x;
   if (lat) *lat = y;

   return 0;
}

#if 0
/* this routine assumes a spherical earth */
static void compute_boresight_lon_lat (const Instr_Geom_Type *geom,
                                       double *lon, double *lat)
{
   double a = GEO_RADIUS / EARTH_RADIUS;
   double alpha, ss, cs, delta;

   alpha = tan (geom->tilt);
   ss = sin(geom->tilt) * cos(geom->tilt) * (a - sqrt(1 - (a*a-1.0) * alpha * alpha));
   cs = sqrt (1.0 - ss * ss);

   /* azi>0 is eastward rotation (CW from north) */
   delta = atan(sin(geom->azi) * ss/cs);

   *lat = acos (cs/cos(delta));
   *lon = delta + geom->sat_lon;
}
#endif

static int slit_vector (const Instr_Geom_Type *geom,
                        const Vector_Type *spacecraft_u,
                        const Vector_Type *anti_boresight,
                        Vector_Type *slit)
{
   Vector_Type s0, s1;

   /* In 3 steps, we'll construct a unit vector along the instrument slit.
    * Start with the instrument boresight pointed at the equator at longitude,
    * geom->sat_lon, so the initial slit vector is a unit vector parallel
    * to the Z coordinate axis (e.g. parallel to the Earth's rotation axis,
    * and pointed northward, toward +Z).
    *   1) Tilt the instrument boresight northward by an angle 'geom->tilt'
    *      radians (this is a rotation about spacecraft roll axis).  At this
    *      point, s0 is just the radius vector for theta=tilt, phi=sat_lon.
    *   2) Now, using the satellite's radius vector as the rotation axis,
    *      rotate the instrument by an angle 'geom->azi'.
    *      According to the right-hand rule, azi<0, will shift the
    *      boresight eastward.
    *   3) Now, using the anti-boresight vector as a rotation axis,
    *      rotate the instrument by angle 'geom->azi'.
    *      According to the right-hand rule, azi>0 will rotate the
    *      slit direction northward, toward the Earth's rotation axis.
    */
   vec_unit (geom->tilt, geom->sat_lon, &s0);

   /* rotate about the spacecraft's radius vector direction: */
   vec_rotate (&s0, spacecraft_u, -geom->azi, &s1);

   /* rotate about the anti-boresight direction */
   vec_rotate (&s1, anti_boresight, geom->azi, slit);

   /* Finally, normalize to ensure that we have a unit vector: */
   return vec_norm (slit);
}

static int init_geom_vectors (Lps_Type *lps)
{
   Instr_Geom_Type *geom = &lps->geom;
   Vector_Type sc, rb_u, rb, anti_bs, slit_u;
   double bs_lon, bs_lat;

   /* FIXME:
    * Ball hasn't yet specified how the instrument reference plane
    * (IRP) is defined so, to make progress, I'm arbitrarily defining
    * the IRP in terms of a vector along the instrument slit.
    * When Ball provides a definition, this may have to change.
    *
    * Assume that the instrument reference plane (IRP) contains the
    * boresight direction and the instrument slit.
    * The unit vector, N(irp), normal to the IRP, is then:
    *     N(irp) = \v x \slit
    * where:
    *     \v = (G\S - E\b)/d
    * is a unit vector along the anti-boresight direction, and where:
    *      G = geostationary orbit radius
    *     \S = unit vector from Earth center toward spacecraft
    *      E = earth radius (assumed spherical)
    *     \b = unit vector from Earth center toward boresight surface point
    *      d = distance from spacecraft to boresight surface point
    *
    *  \slit = unit vector along instrument slit (toward north)
    */

   /* Compute (lon,lat) of boresight surface point [radians]
    * by projecting from mirror coordinates: (mirror_x,mirror_y) = (0,0)
    */
   if (0 != mirror_xy_to_lonlat (lps, &bs_lon, &bs_lat, 0.0, 0.0))
     return -1;

   if (0) fprintf (stderr, "boresight: lon=%g, lat=%g\n",
                   bs_lon/DEGTORAD, bs_lat/DEGTORAD);

   /* vector from the Earth's center to the boresight surface point */
   vec_unit (M_PI_2 - bs_lat, bs_lon, &rb_u);
   vec_scale (EARTH_RADIUS, &rb_u, &rb);

   /* vector from Earth's center toward spacecraft */
   vec_unit (M_PI_2, geom->sat_lon, &geom->spacecraft_u);
   vec_scale (GEO_RADIUS, &geom->spacecraft_u, &sc);

   /* unit vector from the boresight surface point toward spacecraft */
   vec_diff (&sc, &rb, &anti_bs);
   if (0 != vec_norm (&anti_bs))
     return -1;

   /* slit = unit vector "northward" along the instrument slit */
   if (0 != slit_vector (geom, &geom->spacecraft_u, &anti_bs, &slit_u))
     return -1;

   /* Westward IRP normal is cross product normalized to unit length. */
   vec_cross (&anti_bs, &slit_u, &geom->irp_u);
   if (0 != vec_norm (&geom->irp_u))
     return -1;

   if (0) fprintf (stderr, "IRP unit vector: (%g, %g, %g)  len=%g\n",
                   geom->irp_u.x, geom->irp_u.y, geom->irp_u.z,
                   vec_length(&geom->irp_u));
   return 0;
}

static int irp_lmp_angle (const Lps_Type *lps, double lon, double lat, double *angle)
{
   const Instr_Geom_Type *geom = &lps->geom;
   Vector_Type zenith_u, lmp_u;
   double theta, phi;

   /* M_PI_2 = pi/2, input (lon,lat) in radians */
   theta = M_PI_2 - lat;
   phi = lon;

   /* unit vector from Earth center toward surface point */
   vec_unit (theta, phi, &zenith_u);

   /* lmp_u = (westward) unit vector normal to local meridian plane */
   vec_cross (&geom->spacecraft_u, &zenith_u, &lmp_u);
   if (0 != vec_norm (&lmp_u))
     return -1;

   /* dot product of IRP, LMP unit (westward) normals
    * gives the cosine of the desired angle [radians]
    */
   *angle = acos(vec_dot(&geom->irp_u, &lmp_u));

   return 0;
}

static int bsearch_d (double t, double *x, int n)
{
   int n0, n1, n2;
   double xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

static int need_temp_space (Lps_Type *lps, size_t num_needed)
{
   double *tmp;

   if (num_needed < lps->size_tmp)
     return 0;

   if (NULL == (tmp = REALLOC (lps->tmp, 2*num_needed * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: realloc failed", __func__);
        return -1;
     }

   lps->size_tmp = num_needed;
   lps->tmp = tmp;
   lps->num_tmp = 0;

   return 0;
}

int lps_eval (Lps_Type *lps, int band_index, int xtrack,
              double lon, double lat, int num_wave, const double *wave,
              double *lpsens, double *angmax, double *lmp_irp_angle)
{
   Lps_Table_Type *tb_array = NULL;
   Lps_Table_Type *tb0 = NULL;
   Lps_Table_Type *tb1 = NULL;
   int num_mirror_x = lps->num_mirror_x;
   double min_mirror_x = lps->mirror_x[0];
   double max_mirror_x = lps->mirror_x[num_mirror_x-1];
   double *tmp_lpsens, *lpsens0, *lpsens1;
   double *tmp_angmax, *angmax0, *angmax1;
   double mirror_x, wt_m;
   int j, m, xtrack_offset, err, err_count;
   size_t i, num_wave_tb0;

   switch (band_index)
     {
      case TEMPO_BAND_UV:  tb_array = lps->table_uv;
        break;
      case TEMPO_BAND_VIS:  tb_array = lps->table_vis;
        break;
      default:
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: unsupported: band_index=%d", __func__, band_index);
        return -1;
     }

   /* convert (lon,lat) to radians */
   lon *= DEGTORAD;
   lat *= DEGTORAD;

   if (0 != irp_lmp_angle (lps, lon, lat, lmp_irp_angle))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: computing IRU-LMP angle for lon=%g lat=%g",
                     __func__, lon, lat);
        return -1;
     }

   if (0 != lonlat_to_mirror_xy (lps, lon, lat, &mirror_x, NULL))
     return -1;

   if ((mirror_x < min_mirror_x) || (max_mirror_x < mirror_x))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: out of range: mirror_x = %g", __func__, mirror_x);
        return -1;
     }

   /* table indices (m,m+1) bracket mirror_x */
   m = bsearch_d (mirror_x, lps->mirror_x, num_mirror_x);

   /* weight for linear interpolation in mirror_x */
   wt_m = ((lps->mirror_x[m+1] - mirror_x)
           / (lps->mirror_x[m+1] - lps->mirror_x[m]));

   tb0 = &tb_array[m];
   tb1 = &tb_array[m+1];

   xtrack_offset = xtrack * num_wave;

   lpsens0 = tb0->lpsens + xtrack_offset;
   lpsens1 = tb1->lpsens + xtrack_offset;

   angmax0 = tb0->angmax + xtrack_offset;
   angmax1 = tb1->angmax + xtrack_offset;

   num_wave_tb0 = tb0->num_wave;  /* size_t */
   if (0 != need_temp_space (lps, num_wave_tb0))
     return -1;

   tmp_lpsens = lps->tmp;
   tmp_angmax = lps->tmp + num_wave_tb0;

   /* interpolate along mirror_x */
   for (i = 0; i < num_wave_tb0; i++)
     {
        tmp_lpsens[i] = wt_m * lpsens0[i] + (1.0 - wt_m) * lpsens1[i];
     }

   for (i = 0; i < num_wave_tb0; i++)
     {
        tmp_angmax[i] = wt_m * angmax0[i] + (1.0 - wt_m) * angmax1[i];
     }

   /* map tabulated wavelengths to data wavelength grid */
   if (lps->interp == NULL)
     {
        if ((NULL == (lps->interp = gsl_interp_alloc (gsl_interp_linear, num_wave_tb0)))
            || (NULL == (lps->acc = gsl_interp_accel_alloc ())))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: gsl_interp alloc failed", __func__);
             return -1;
          }
     }

   err_count = 0;

   gsl_interp_accel_reset (lps->acc);
   if (gsl_interp_init (lps->interp, tb0->wave, tmp_lpsens, num_wave_tb0))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gsl_interp init failed", __func__);
        return -1;
     }

   for (j = 0; j < num_wave; j++)
     {
        double lp_sens = 0.0;
        err = gsl_interp_eval_e (lps->interp, tb0->wave, tmp_lpsens, wave[j], lps->acc, &lp_sens);
        if (err)
          {
             err_count++;
             lp_sens = 0.0;
          }
        lpsens[j] = lp_sens;
     }

   if (err_count)
     {
        tell_vwarn (1, "%s: lpsens interpolation error at %d/%d wavelength points",
                    __func__, err_count, num_wave);
        err_count = 0;
     }

   gsl_interp_accel_reset (lps->acc);
   if (gsl_interp_init (lps->interp, tb0->wave, tmp_angmax, num_wave_tb0))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: gsl_interp init failed", __func__);
        return -1;
     }

   for (j = 0; j < num_wave; j++)
     {
        double ang_max = 0.0;
        err = gsl_interp_eval_e (lps->interp, tb0->wave, tmp_angmax, wave[j], lps->acc, &ang_max);
        if (err)
          {
             err_count++;
             ang_max = 0.0;
          }
        angmax[j] = ang_max;
     }

   if (err_count)
     {
        tell_vwarn (1, "%s: angmax interpolation errors at %d/%d wavelength points",
                    __func__, err_count, num_wave);
     }

   return 0;
}

static void free_table_type (Lps_Table_Type *tbl)
{
   if (tbl == NULL)
     return;
   FREE(tbl->lpsens);
   FREE(tbl->angmax);
   FREE(tbl->wave);
}

static int alloc_table_type (Lps_Table_Type *tbl, int num_xtrack, int num_wave)
{
   size_t len = num_xtrack * num_wave * sizeof(double);
   size_t len_wave = num_wave * sizeof(double);

   if ((NULL == (tbl->lpsens = (double *)MALLOC (len)))
       || (NULL == (tbl->angmax = (double *)MALLOC (len)))
       || (NULL == (tbl->wave = (double *)MALLOC (len_wave))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   tbl->num_xtrack = num_xtrack;
   tbl->num_wave = num_wave;

   return 0;
}

static void free_table_array (Lps_Table_Type *tbl, int num_mirror_x)
{
   int i;
   if (tbl == NULL)
     return;
   for (i = 0; i < num_mirror_x; i++)
     {
        free_table_type (&tbl[i]);
     }
   FREE(tbl);
}

static Lps_Table_Type *alloc_table_array (int num_mirror_x, int num_xtrack, int num_wave)
{
   Lps_Table_Type *tbl = NULL;
   size_t len = num_mirror_x * sizeof (Lps_Table_Type);
   int i;

   if (NULL == (tbl = (Lps_Table_Type *)MALLOC (len)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset((char *)tbl, 0, len);

   for (i = 0; i < num_mirror_x; i++)
     {
        if (0 != alloc_table_type (&tbl[i], num_xtrack, num_wave))
          {
             free_table_array (tbl, num_mirror_x);
             return NULL;
          }
     }

   return tbl;
}

static int read_lps_mirror_x_grid (int ncid, double **xp, int *nxp)
{
   size_t i, dimlen;
   int dimid, start, count;
   double *x;

   if (0 != TIO_inq_dim (ncid, "azimuth_lps", &dimid, &dimlen))
     return -1;

   if (NULL == (x = (double *)MALLOC (dimlen * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   start = 0;
   count = dimlen;

   if (0 != TIO_get_var_section (ncid, "azimuth_lps", &start, &count, NC_DOUBLE, x))
     {
        FREE(x);
        return -1;
     }

   /* convert deg -> radians */
   for (i = 0; i < dimlen; i++)
     {
        x[i] *= DEGTORAD;
     }

   *xp = x;
   *nxp = dimlen;

   return 0;
}

static int read_lps_table_array1 (int grp, Lps_Table_Type **tblp)
{
   Lps_Table_Type *tbl = NULL;
   TIO_Var_Info_Type info = {0};
   const char lps_var[] = "linear_polarization_sensitivity";
   const char angmax_var[] = "angle_of_maximum_transmission";
   const char wave_var[] = "wavelength";
   int i, num_mirror_x, num_xtrack, num_wave;
   int start[3], count[3];

   if (0 != TIO_inq_var (grp, lps_var, &info))
     return -1;

   /* dimlen[0] is slowest varying,
    * dimlen[ndims-1] is fastest varying */
   num_mirror_x = info.dimlens[0];
   num_xtrack = info.dimlens[1];
   num_wave = info.dimlens[2];

   if (NULL == (tbl = alloc_table_array (num_mirror_x, num_xtrack, num_wave)))
     return -1;

   for (i = 0; i < num_mirror_x; i++)
     {
        Lps_Table_Type *ti = &tbl[i];
        int k;

        ti->num_xtrack = num_xtrack;
        ti->num_wave = num_wave;

        start[0] = 0;
        count[0] = num_wave;

        if (0 != TIO_get_var_section (grp, wave_var, start, count, NC_DOUBLE, ti->wave))
          goto return_error;

        start[0] = i;
        start[1] = 0;
        start[2] = 0;

        count[0] = 1;
        count[1] = num_xtrack;
        count[2] = num_wave;

        if ((0 != TIO_get_var_section (grp, lps_var, start, count, NC_DOUBLE, ti->lpsens))
            ||(0 != TIO_get_var_section (grp, angmax_var, start, count, NC_DOUBLE, ti->angmax)))
          goto return_error;

        /* convert deg -> radians */
        for (k = 0; k < num_xtrack * num_wave; k++)
          {
             ti->angmax[k] *= DEGTORAD;
          }
     }

   *tblp = tbl;
   return 0;
return_error:
   free_table_array (tbl, num_mirror_x);
   return -1;
}

static int read_lps_tables (Lps_Type *lps, int ncid)
{
   int grp;

   if (0 != read_lps_mirror_x_grid (ncid, &lps->mirror_x, &lps->num_mirror_x))
     return -1;

   if (0 != TIO_inq_grp (ncid, TEMPO_BAND_NAME_UV, &grp))
     return -1;

   if (0 != read_lps_table_array1 (grp, &lps->table_uv))
     return -1;

   if (0 != TIO_inq_grp (ncid, TEMPO_BAND_NAME_VIS, &grp))
     return -1;

   if (0 != read_lps_table_array1 (grp, &lps->table_vis))
     return -1;

   return 0;
}

static void free_lps (Lps_Type *lps)
{
   if (NULL == lps)
     return;
   lps_proj_close (lps);
   FREE(lps->mirror_x);
   free_table_array (lps->table_uv, lps->num_mirror_x);
   free_table_array (lps->table_vis, lps->num_mirror_x);
   FREE(lps->tmp);
   gsl_interp_free (lps->interp);
   gsl_interp_accel_free (lps->acc);
   FREE(lps);
}

void lps_close (Lps_Type *lps)
{
   free_lps (lps);
}

Lps_Type *lps_open (config_t *cfg, double sat_lon_radians)
{
   Lps_Type *lps = NULL;
   wordexp_t we;
   config_setting_t *s;
   const char *lps_file;
   int ncid, status;

   memset ((char *)&we, 0, sizeof(wordexp_t));

   if (NULL == (s = config_lookup (cfg, "polcorr_tables")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing polcorr in param file: %s",
                     __func__, config_error_file (cfg));
        return NULL;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "lps_lut", &lps_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing %s in param file: %s",
                     __func__, config_setting_name(s),
                     config_setting_source_file (s));
        return NULL;
     }

   if (NULL == (lps = (Lps_Type *)MALLOC (sizeof *lps)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)lps, 0, sizeof *lps);

   lps->geom.sat_lon = sat_lon_radians;

   if (0 != read_geom (cfg, &lps->geom))
     goto return_error;

   if (0 != lps_proj_open (lps))
     goto return_error;

   if (0 != init_geom_vectors (lps))
     goto return_error;

   if ((0 != wordexp (lps_file, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, lps_file);
        goto return_error;
     }

   if (0 != TIO_open (we.we_wordv[0], NC_NOWRITE, &ncid))
     goto return_error;

   status = read_lps_tables (lps, ncid);

   (void) TIO_close (ncid);
   wordfree (&we);

   if (status) goto return_error;

   return lps;

return_error:
   wordfree (&we);
   free_lps (lps);
   return NULL;
}
