#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

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
#define GEO_ALTITUDE  35785831.0   /* meters */

#define PROJ_ARGS_BUFSIZE       80

/* FIXME - these indices should be provided by tio_template.h */
enum
{
   TEMPO_BAND_UV = 0,
   TEMPO_BAND_VIS = 1
};

typedef struct
{
   double sat_lon;  /**< GEO satellite orbital station [deg] */

   double tilt;
   /**< tilt > 0 is northward tilt of instrument boresight
    *   about spacecraft roll axis [deg] */

   double azi;
   /**< azi > 0 is rotation, eastward from north (CW),
    *   about instrument boresight axis [deg] */
}
Instr_Geom_Type;

typedef struct
{
   double *lpsens;   /**< linear polarization sensivity [num_xtrack, num_wave] */
   double *angmax;   /**< angle of maximum transmission [num_xtrack, num_wave] */
   double *wave;     /**< wavelength grid [num_wave] */
   int num_xtrack;
   int num_wave;
}
Lps_Table_Type;

struct Lps_Type
{
   Instr_Geom_Type geom;
   double *mirror_x;             /**< azimuth angle from boresight, EAST is positive [rad] */
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

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "sat_lon", &geom->sat_lon))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "tilt", &geom->tilt))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "azi", &geom->azi)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing %s in param file: %s",
                     __func__, config_setting_name(s),
                     config_setting_source_file (s));
        return -1;
     }

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
                   geom->sat_lon, GEO_ALTITUDE,
                   geom->tilt,
                   geom->azi);
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

/* input (lon,lat) in degrees; output (x,y) in radians */
static int lonlat_to_tpers_xy (Lps_Type *lps,
                               int n, const double *lon, const double *lat,
                               double *x, double *y)
{
   int i, status;

   /* convert (lon,lat) to radians */
   for (i = 0; i < n; i++)
     {
        x[i] = lon[i] * DEGTORAD;
        y[i] = lat[i] * DEGTORAD;
     }

   if ((status = pj_transform (lps->longlat, lps->tpers, n, 1, x, y, NULL)) != 0)
     {
        tell_verror (TELL_APPLICATION_ERROR,
                     "%s: pj_transform failed, status = %d (%s)",
                     __func__, status, pj_strerrno(status));
        return -1;
     }

   /* convert (x,y) to radians */
   for (i = 0; i < n; i++)
     {
        x[i] /= GEO_ALTITUDE;
        y[i] /= GEO_ALTITUDE;
     }

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
              double *lpsens, double *angmax)
{
   Lps_Table_Type *tb_array = NULL;
   Lps_Table_Type *tb0 = NULL;
   Lps_Table_Type *tb1 = NULL;
   int num_mirror_x = lps->num_mirror_x;
   double min_mirror_x = lps->mirror_x[0];
   double max_mirror_x = lps->mirror_x[num_mirror_x-1];
   double *tmp_lpsens, *lpsens0, *lpsens1;
   double *tmp_angmax, *angmax0, *angmax1;
   double mirror_x, y, wt_m;
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

   /* (x,y) are angular coordinates in the tilted perspective plane (Proj4 'tpers')
    * x [radian] = azimuth angle from boresight
    * y [radian] = altitude angle from boresight
    *
    * SMA mirror coordinates have +X "east" and +Y "SOUTH", so that:
    *        mirror_x = x
    *        mirror_y = -y
    * But since we don't need mirror_y, I'll not bother with the sign flip.
    */
   if (0 != lonlat_to_tpers_xy (lps, 1, &lon, &lat, &mirror_x, &y))
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
        tell_vwarn (TELL_MSGTYPE_ERROR, "%s: %d lpsens interpolation errors",
                    __func__, err_count);
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
        tell_vwarn (TELL_MSGTYPE_ERROR, "%s: %d angmax interpolation errors",
                    __func__, err_count);
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
   size_t dimlen;
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

Lps_Type *lps_open (config_t *cfg)
{
   Lps_Type *lps = NULL;
   config_setting_t *s;
   const char *lps_file;
   int ncid, status;

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

   if (0 != read_geom (cfg, &lps->geom))
     goto return_error;

   if (0 != lps_proj_open (lps))
     goto return_error;

   if (0 != TIO_open (lps_file, NC_NOWRITE, &ncid))
     goto return_error;
   status = read_lps_tables (lps, ncid);
   (void) TIO_close (ncid);
   if (status) goto return_error;

   return lps;

return_error:
   free_lps (lps);
   return NULL;
}
