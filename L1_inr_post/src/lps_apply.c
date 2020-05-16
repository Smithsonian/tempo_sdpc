/** @file lps_apply.c
 *  @brief Adjust synthetic radiances using linear polarization sensitivity
 *         lookup table
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "lps_apply.h"

typedef struct
{
   double *lon;
   double *lat;
   size_t num_xtrack;
   size_t num_step;
}
Geoloc_Type;

typedef struct
{
   float *wave;
   float *rad;
   size_t num_waves;
}
Radiance_Type;

typedef struct
{
   float *q;
   float *u;
   double *wave;
   double *lpsens;
   double *angmax;
   size_t num_waves;
}
QU_Type;

static void free_geoloc_type (Geoloc_Type *g)
{
   if (g == NULL)
     return;
   FREE(g->lon);
   FREE(g->lat);
   FREE(g);
}

static Geoloc_Type *alloc_geoloc_type (size_t num_step, size_t num_xtrack)
{
   Geoloc_Type *g = NULL;
   size_t len = num_xtrack * num_step * sizeof(double);

   if (NULL == (g = (Geoloc_Type *)MALLOC (sizeof *g)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)g, 0, sizeof (*g));

   if ((NULL == (g->lon = (double *)MALLOC (len)))
       || (NULL == (g->lat = (double *)MALLOC (len))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_geoloc_type (g);
        return NULL;
     }

   g->num_xtrack = num_xtrack;
   g->num_step = num_step;

   return g;
}

static void free_radiance_type (Radiance_Type *rt)
{
   if (rt == NULL)
     return;
   FREE(rt->wave);
   FREE(rt->rad);
   FREE(rt);
}

static Radiance_Type *alloc_radiance_type (size_t num_xtrack, size_t num_waves)
{
   Radiance_Type *rt = NULL;
   size_t len = num_xtrack * num_waves * sizeof(float);

   if (NULL == (rt = (Radiance_Type *)MALLOC (sizeof *rt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)rt, 0, sizeof (*rt));

   if ((NULL == (rt->rad = (float *)MALLOC (len)))
       || (NULL == (rt->wave = (float *)MALLOC (len))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_radiance_type (rt);
        return NULL;
     }

   rt->num_waves = num_waves;

   return rt;
}

static void free_qu_type (QU_Type *qu)
{
   if (qu == NULL)
     return;
   FREE(qu->q);
   FREE(qu->u);
   FREE(qu->wave);
   FREE(qu->lpsens);
   FREE(qu->angmax);
   FREE(qu);
}

static QU_Type *alloc_qu_type (size_t num_xtrack, size_t num_waves)
{
   QU_Type *qu = NULL;
   size_t len_f = num_waves * sizeof(float);
   size_t len_d = num_waves * sizeof(double);

   if (NULL == (qu = (QU_Type *)MALLOC (sizeof *qu)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)qu, 0, sizeof (*qu));

   if ((NULL == (qu->q = (float *)MALLOC (num_xtrack * len_f)))
       || (NULL == (qu->u = (float *)MALLOC (num_xtrack * len_f)))
       || (NULL == (qu->wave = (double *)MALLOC (len_d)))
       || (NULL == (qu->lpsens = (double *)MALLOC (len_d)))
       || (NULL == (qu->angmax = (double *)MALLOC (len_d))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free_qu_type (qu);
        return NULL;
     }

   qu->num_waves = num_waves;

   return qu;
}

static int lps_apply_band (Lps_Type *lps, int ncid, int band_index)
{
   TIO_Var_Info_Type info = {0};
   const char *band_name = NULL;
   Geoloc_Type *g = NULL;
   Radiance_Type *rt = NULL;
   QU_Type *qu = NULL;
   size_t num_xtrack, num_step, num_waves;
   size_t i, xtrack, step;
   int start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   int grp, status = -1;

   switch (band_index)
     {
      case TEMPO_BAND_UV:
        band_name = TEMPO_BAND_NAME_UV;
        break;
      case TEMPO_BAND_VIS:
        band_name = TEMPO_BAND_NAME_VIS;
        break;
      default:
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: invalid value band_index=%d",
                     __func__, band_index);
        return -1;
     }

   if (0 != TIO_inq_grp (ncid, band_name, &grp))
     return -1;

   if (0 != TIO_inq_var (grp, TEMPO_VAR_RADIANCE, &info))
     return -1;

   /* right-most index varies fastest */
   num_step = info.dimlens[0];
   num_xtrack = info.dimlens[1];
   num_waves = info.dimlens[2];

   if ((NULL == (g = alloc_geoloc_type (num_step, num_xtrack)))
       || (NULL == (rt = alloc_radiance_type (num_xtrack, num_waves)))
       || (NULL == (qu = alloc_qu_type (num_xtrack, num_waves))))
     {
        goto return_status;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = num_step;
   count[1] = num_xtrack;

   /* I'm using truth_lon,truth_lat because this code will only
    * ever apply to synthetic data.  Using the 'truth' (lon,lat)
    * coordinates also allows us to modify the radiances prior to INR.
    */
   if ((0 != TIO_get_var_section (grp, "truth_lat", start, count, NC_DOUBLE, g->lat))
       || (0 != TIO_get_var_section (grp, "truth_lon", start, count, NC_DOUBLE, g->lon)))
     goto return_status;

   for (step = 0; step < num_step; step++)
     {
        start[0] = step;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = num_xtrack;
        count[2] = num_waves;

        if ((0 != TIO_get_var_section (grp, TEMPO_VAR_WAVELENGTH, start, count, NC_FLOAT, rt->wave))
            || (0 != TIO_get_var_section (grp, TEMPO_VAR_RADIANCE, start, count, NC_FLOAT, rt->rad))
            || (0 != TIO_get_var_section (grp, "q_truth", start, count, NC_FLOAT, qu->q))
            || (0 != TIO_get_var_section (grp, "u_truth", start, count, NC_FLOAT, qu->u)))
          goto return_status;

        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             float *wave_x, *rad_x, *q_x, *u_x;
             double lonlat_fill_value = 0.0;
             double delta_irp;
             int k;

             k = xtrack + step * num_xtrack;

             if ((g->lat[k] == lonlat_fill_value)
                 || (g->lon[k] == lonlat_fill_value))
               continue;

             wave_x = rt->wave + xtrack * num_waves;
             for (i = 0; i < num_waves; i++)
               {
                  qu->wave[i] = (double) wave_x[i];
               }

             if (0 != lps_eval (lps, band_index, xtrack, g->lon[k], g->lat[k],
                                num_waves, qu->wave, qu->lpsens, qu->angmax,
                                &delta_irp))
               {
                  tell_verror (TELL_RUNTIME_ERROR,
                               "%s: lps_eval failed step=%ld xtrack=%ld",
                               __func__, step, xtrack);
                  goto return_status;
               }

             q_x = qu->q + xtrack * num_waves;
             u_x = qu->u + xtrack * num_waves;
             rad_x = rt->rad + xtrack * num_waves;

             for (i = 0; i < num_waves; i++)
               {
                  double dolp, pa, lpserr, q, u;

                  /* skip invalid radiances
                   * (assumes radiance fill value < 0) */
                  if (rad_x[i] <= 0.0)
                    continue;

                  /* File contains
                   *   radiance = I = I(0) + I(90)
                   *   q_truth  = Q = I(0) - I(90)
                   *   u_truth  = U = I(45) - I(135)
                   * but we want q = Q/I, u=Q/I,
                   * so divide by radiance:
                   */
                  q = q_x[i] / rad_x[i];
                  u = u_x[i] / rad_x[i];

                  /* degree of linear polarization */
                  dolp = sqrt (q*q + u*u);

                  /* position angle of linear polarization
                   * relative to instrument reference plane */
                  pa = 0.5 * atan2 (u, q);
                  if (pa < 0.0) pa += M_PI;
                  pa += delta_irp;

                  /* radiance error due to linear polarization sensitivity */
                  lpserr = (2.0 * qu->lpsens[i] * dolp
                            * cos (2.0 * (pa - qu->angmax[i])));

                  /* (measured radiance) = (true radiance)*(1 + lpserr) */
                  rad_x[i] *= 1.0 + lpserr;
               }
          }

        /* write measured radiance spectrum */
        if (0 != TIO_put_var_section (grp, TEMPO_VAR_RADIANCE, start, count, NC_FLOAT, rt->rad))
          goto return_status;
     }

   status = 0;
return_status:
   free_geoloc_type (g);
   free_radiance_type (rt);
   free_qu_type (qu);

   return status;
}

int lps_apply (Lps_Type *lps, const char *rad_file)
{
   static int bands[] = {TEMPO_BAND_UV, TEMPO_BAND_VIS};
   unsigned int i, num_bands = sizeof(bands)/sizeof(*bands);
   int ncid, status = -1;

   if (0 != TIO_open (rad_file, NC_WRITE, &ncid))
     return -1;

   if (0 != tio_history_append_cmdline (ncid))
     goto return_status;

   for (i = 0; i < num_bands; i++)
     {
        if (0 != lps_apply_band (lps, ncid, bands[i]))
          goto return_status;
     }

   status = 0;
return_status:
   (void) TIO_close (ncid);

   return status;
}
