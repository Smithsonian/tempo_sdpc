#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "ccd.h"
#include "granule.h"
#include "util.h"

#define NUM_OCTANTS (8)

/* FIXME - better if these are not hard-coded? */

#define NUM_READOUT_SIGNALS (1<<14)   /* assuming a 14-bit pixel at read-out */

#define SATURATION_FUDGE_FACTOR      (0.9)
#define READOUT_SATURATION_THRESHOLD ((NUM_READOUT_SIGNALS-1)*SATURATION_FUDGE_FACTOR)

#define NUM_PIXEL_SAMPLES (10000)

typedef struct
{
   double *octant_means[NUM_OCTANTS];
   double *nonlin_frac[NUM_OCTANTS];
   double *weights;
   double *exposure_time_per_frame;
   size_t num_times;

   double c0[NUM_OCTANTS];
   double c1[NUM_OCTANTS];
}
Trend_Type;

typedef struct
{
   int *signal_input;
   double *signal_output;
   double *nonlin_frac;
   int num_signals;
}
Lin_Corr_Type;

static void free_trend (Trend_Type *tt)
{
   int i;

   if (tt == NULL)
     return;

   FREE(tt->exposure_time_per_frame);
   FREE(tt->weights);

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        FREE(tt->octant_means[i]);
        FREE(tt->nonlin_frac[i]);
     }
}

static int alloc_trend (int num_times, Trend_Type *tt)
{
   int i;

   if ((NULL == (tt->exposure_time_per_frame = (double *)MALLOC (num_times * sizeof(double))))
       || (NULL == (tt->weights = (double *)MALLOC (num_times * sizeof(double)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   for (i = 0; i < NUM_OCTANTS; i++)
     {
        if ((NULL == (tt->octant_means[i] = (double *)MALLOC (num_times * sizeof(double))))
            ||(NULL == (tt->nonlin_frac[i] = (double *)MALLOC (num_times * sizeof(double)))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             free_trend (tt);
             return -1;
          }
     }

   tt->num_times = num_times;

   return 0;
}

static void write_indexed_img (int i, const char *prefix, const Image_Type *img)
{
   char buf[256];
   if (i != 15) return;
   (void) snprintf (buf, sizeof(buf), "%s.%d", prefix, i);
   image_write_raw (img, buf);
}

static double *_pSort_Doubles;
static int index_sort_doubles (const void *va, const void *vb)
{
   int ia = *(const int *)va;
   int ib = *(const int *)vb;
   double a = _pSort_Doubles[ia];
   double b = _pSort_Doubles[ib];
   if (a < b) return -1;
   else if (a > b) return +1;
   else return 0;
}

static int *determine_sweep_order (const Granule_Type *gr, int num_exprecs)
{
   double *exposure_per_frame = NULL;
   int *sweep_order = NULL;
   int i, status = -1;

   if ((NULL == (exposure_per_frame = (double *)MALLOC (num_exprecs * sizeof(double))))
       || (NULL == (sweep_order = (int *)MALLOC (num_exprecs * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto free_and_return;
     }

   if (0 != gr->granule_get_exposure_per_frame (gr, exposure_per_frame))
     goto free_and_return;

   for (i = 0; i < num_exprecs; i++)
     {
        sweep_order[i] = i;
     }

   _pSort_Doubles = exposure_per_frame;
   qsort ((void *)sweep_order, (size_t) num_exprecs, sizeof(int), index_sort_doubles);
   _pSort_Doubles = NULL;

   status = 0;
free_and_return:
   FREE(exposure_per_frame);
   if (status)
     {
        FREE(sweep_order);
        sweep_order = NULL;
     }

   return sweep_order;
}

static int measure_trends (Granule_Type *gr, CCD_Linearity_Type *clt, Trend_Type *tt)
{
   Granule_Exprec_Type *exprec = NULL;
   CCD_Select_Type *sel = NULL;
   size_t num_sample = NUM_PIXEL_SAMPLES;
   int *sweep_order = NULL;
   int i, o, num_exprecs;

   num_exprecs = gr->granule_num_exprecs(gr);

   if (NULL == (sweep_order = determine_sweep_order (gr, num_exprecs)))
     goto return_status;

   if (NULL == (sel = clt_select_alloc (clt, num_sample)))
     goto return_status;

   for (i = 0; i < num_exprecs; i++)
     {
        double octant_means[NUM_OCTANTS];

        if (NULL == gr->granule_read_exprec_by_index (gr, sweep_order[i], &exprec))
          goto return_status;

        if (0 != clt->clt_correct_coadd (clt, exprec->num_coadds, exprec->img))
          goto return_status;
        tt->exposure_time_per_frame[i] = exprec->exposure_time / exprec->num_coadds;

        fprintf (stderr, "exposure record %3d/%d: %8.3f msec/frame\n",
                 i, num_exprecs, 1.e3 * tt->exposure_time_per_frame[i]);

        if (0) write_indexed_img (i, "uncorrected", exprec->img);

        if (0 != clt->clt_correct_offset (clt, exprec->img))
          goto return_status;

        if (0) write_indexed_img (i, "corrected", exprec->img);

        if (0 != clt->clt_trimmed_sample_mean (clt, exprec->img, sel, octant_means))
          goto return_status;

        for (o = 0; o < NUM_OCTANTS; o++)
          {
             tt->octant_means[o][i] = octant_means[o];
          }
     }

return_status:
   FREE(sweep_order);
   gr->granule_free_exprec (exprec);
   clt_select_free (sel);

   return 0;

}

static int write_fit (int octant, size_t num_times, double *exposure_time,
                      double *means_oct, double *nonlin_frac, double c0, double c1)
{
   FILE *fp = NULL;
   char buf[1024];
   size_t k;

   snprintf (buf, sizeof(buf), "octant_%d.dat", octant);

   if (NULL == (fp = fopen (buf, "w")))
     {
        fprintf (stderr, "*** Error: opening %s for writing\n", buf);
        return -1;
     }

   fprintf (fp, "# c0 = %14.6f  c1 = %14.6f\n", c0, c1);
   for (k = 0; k < num_times; k++)
     {
        double linear_signal_k = c0 + c1 * exposure_time[k];
        fprintf (fp, "%10.3f %10.3f %12.4e %10.3f\n",
                 exposure_time[k] * 1.e3, means_oct[k], nonlin_frac[k], linear_signal_k);
     }

   return fclose (fp);
}

static int fit_trends (Trend_Type *tt)
{
   double *weights = tt->weights;
   int o;

   for (o = 0; o < NUM_OCTANTS; o++)
     {
        double *means_oct = tt->octant_means[o];
        double *nonlin_frac = tt->nonlin_frac[o];
        double c0, c1, cov00, cov01, cov11, sumsq;
        int gsl_status;
        size_t k;

        for (k = 0; k < tt->num_times; k++)
          {
             if (means_oct[k] < READOUT_SATURATION_THRESHOLD)
               weights[k] = 1.0;
             else
               weights[k] = 0.0;
          }

        gsl_status = gsl_fit_wlinear (tt->exposure_time_per_frame, 1, weights, 1, means_oct, 1, tt->num_times,
                                      &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
        if (gsl_status)
          {
             fprintf (stderr, "*** Error: gsl_fit_linear: %s\n", gsl_strerror (gsl_status));
             return -1;
          }

        tt->c0[o] = c0;
        tt->c1[o] = c1;

        for (k = 0; k < tt->num_times; k++)
          {
             double linear_signal_k = c0 + c1 * tt->exposure_time_per_frame[k];
             nonlin_frac[k] = (linear_signal_k - means_oct[k]) / READOUT_SATURATION_THRESHOLD;
          }

        if (0) write_fit (o, tt->num_times, tt->exposure_time_per_frame, means_oct, nonlin_frac, c0, c1);
     }

   return 0;
}

static void free_lin_corr (Lin_Corr_Type *lin)
{
   if (lin == NULL)
     return;
   FREE(lin->signal_input);
   FREE(lin->signal_output);
   FREE(lin->nonlin_frac);
}

static int alloc_lin_corr (int num_signals, Lin_Corr_Type *lin)
{
   if ((NULL == (lin->signal_input = (int *)MALLOC (num_signals * sizeof(int))))
       ||(NULL == (lin->signal_output = (double *)MALLOC (num_signals * sizeof(double))))
       ||(NULL == (lin->nonlin_frac = (double *)MALLOC (num_signals * sizeof(double)))))
     {
        free_lin_corr (lin);
        return -1;
     }

   lin->num_signals = num_signals;

   return 0;
}

static int create_output_file (const char *output_file, size_t num_signals,
                               TIO_Meta_Type *meta, int *pncid)
{
   const char *comment_string_array[] =
     {"dimensions: ADC=(0,1) -> (odd,even); quad=(0,1,2,3) -> (A,B,C,D)"};
   int num_dn = num_signals;
   int num_quad = 4;
   int num_adc = 2;
   int dimid_dn, dimid_adc, dimid_quad, varid_lut, varid_frac;
   int ncid, dimids_lut[3];

   tell_vlog (TELL_MSGTYPE_INFO, 1, "creating %s", output_file);

   if (0 != TIO_create (output_file, NC_NETCDF4, &ncid))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: creating %s", __func__, output_file);
        return -1;
     }

   /* Choose dimension names and variable names to match the existing
    * calibration key data file
    */
   if ((0 != TIO_def_dim (ncid, "DN", num_dn, &dimid_dn))
       ||(0 != TIO_def_dim (ncid, "ADC", num_adc, &dimid_adc))
       ||(0 != TIO_def_dim (ncid, "quad", num_quad, &dimid_quad)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dimensions in :%s", __func__, output_file);
        (void) TIO_close (ncid);
        return -1;
     }

   if (0 != TIO_put_att (ncid, NC_GLOBAL, "comment", NC_STRING, 1, comment_string_array))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing comment attribute: %s", __func__, output_file);
        (void) TIO_close (ncid);
        return -1;
     }

   dimids_lut[0] = dimid_adc;
   dimids_lut[1] = dimid_quad;
   dimids_lut[2] = dimid_dn;    /* right-most varies fastest */

   if (0 != TIO_def_var (ncid, "nonlinearity_LUT", NC_FLOAT, 3, dimids_lut, &varid_lut))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining variables in :%s", __func__, output_file);
        (void) TIO_close (ncid);
        return -1;
     }

   if (0 != TIO_def_var (ncid, "nonlinear_fraction", NC_FLOAT, 3, dimids_lut, &varid_frac))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining variables in :%s", __func__, output_file);
        (void) TIO_close (ncid);
        return -1;
     }

   if (0 != tio_meta_write_ncattr (meta, ncid))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing metadata in :%s", __func__, output_file);
        (void) TIO_close (ncid);
        return -1;
     }

   *pncid = ncid;

   return 0;
}

static int write_lut (int ncid, int octant, const Lin_Corr_Type *lin)
{
   const char *varname_lut = "nonlinearity_LUT";
   const char *varname_frac = "nonlinear_fraction";
   int start[3], count[3];
   int adc, quad;

   /* 'octant' is the octant index where octants are ordered like so:
    * (Ao,Bo,Co,Do, Ae,Be,Ce,De).
    * The cal file uses indexing (adc,quad)
    * where adc=(0,1) -> (odd,even)
    * and where quad=(0,1,2,3) -> (A,B,C,D).
    * So, map 'octant' to (adc,quad):
    */

   adc  = octant / 4;
   quad = octant % 4;

   start[0] = adc;
   start[1] = quad;
   start[2] = 0;

   count[0] = 1;
   count[1] = 1;
   count[2] = lin->num_signals;

   if ((0 != TIO_put_var_section (ncid, varname_lut, start, count, TIO_DOUBLE, lin->signal_output))
       || (0 != TIO_put_var_section (ncid, varname_frac, start, count, TIO_DOUBLE, lin->nonlin_frac)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing to output file", __func__);
        return -1;
     }

   return 0;
}

static int create_lookup_tables (const Trend_Type *tt, TIO_Meta_Type *meta,
                                 const char *output_file)
{
   Lin_Corr_Type lin = {0};
   gsl_interp *interp = NULL;
   gsl_interp_accel *acc = NULL;
   size_t num_signals = NUM_READOUT_SIGNALS;
   int status = -1;
   int ncid = -1;
   int o;

   if (0 != alloc_lin_corr (num_signals, &lin))
     return -1;

   if (NULL == (acc = gsl_interp_accel_alloc ()))
     goto return_status;

   if (0 != create_output_file (output_file, num_signals, meta, &ncid))
     goto return_status;

   for (o = 0; o < NUM_OCTANTS; o++)
     {
        double *octant_means = tt->octant_means[o];
        double *nonlin_frac = tt->nonlin_frac[o];
        size_t k, i, num_times;

        for (k = tt->num_times-1; k > 0; k--)
          {
             if (octant_means[k] < READOUT_SATURATION_THRESHOLD)
               break;
          }
        num_times = k+1;

        gsl_interp_free (interp);
        if (NULL == (interp = gsl_interp_alloc (gsl_interp_linear, num_times)))
          goto return_status;

        if (0 != gsl_interp_init (interp, octant_means, nonlin_frac, num_times))
          goto return_status;

        for (i = 0; i < num_signals; i++)
          {
             int sig = i;
             double f;
             if (sig < octant_means[0])
               {
                  f = nonlin_frac[0];
               }
             else if (sig > octant_means[num_times-1])
               {
                  f = nonlin_frac[num_times-1];
               }
             else
               {
                  double d_sig = sig;
                  f = gsl_interp_eval (interp, octant_means, nonlin_frac, d_sig, acc);
               }
             lin.nonlin_frac[i] = f;
             lin.signal_input[i] = sig;
             lin.signal_output[i] = sig * (1.0 + f);
          }

        if (0 != gsl_interp_accel_reset (acc))
          goto return_status;

        if (0 != write_lut (ncid, o, &lin))
          goto return_status;
     }

   status = 0;
return_status:
   if (status)
     {
        fprintf (stderr, "%s: error while generating lookup tables\n", __func__);
     }
   free_lin_corr (&lin);
   gsl_interp_free (interp);
   gsl_interp_accel_free (acc);
   if (ncid > 0) (void) TIO_close (ncid);

   return status;
}

static int process_granule (Granule_Type *gr, CCD_Linearity_Type *clt,
                            TIO_Meta_Type *meta, const char *output_file)
{
   Trend_Type tt = {0};
   int num_exprecs, status = -1;

   num_exprecs = gr->granule_num_exprecs (gr);

   if (0 != alloc_trend (num_exprecs, &tt))
     return -1;

   if (0 != measure_trends (gr, clt, &tt))
     goto return_status;

   if (0 != fit_trends (&tt))
     goto return_status;

   if (0 != create_lookup_tables (&tt, meta, output_file))
     goto return_status;

   status = 0;
return_status:
   free_trend (&tt);

   return status;
}

int derive_linearity (const char *input_file, const char *output_file)
{
   Granule_Type *gr = NULL;
   CCD_Linearity_Type *clt = NULL;
   TIO_Meta_Type *meta = NULL;
   int exposure_type, status = -1;

   if (NULL == (gr = granule_open (input_file)))
     return -1;

   if (0 != gr->granule_type (gr, &exposure_type))
     goto return_status;

   if (exposure_type != EXPREC_TYPE_LIN_IRR)
     {
        tell_verror (TELL_INVALID_DATA_ERROR, "%s: invalid exposure type: %d", __func__, exposure_type);
        goto return_status;
     }

   if (NULL == (clt = ccd_linearity_init ()))
     goto return_status;

   if (NULL == (meta = tio_meta_open ()))
     goto return_status;

   if (0 != tio_meta_set_datetime_production (meta))
     goto return_status;

   if (0 != meta_record_basename (meta, input_file))
     goto return_status;

   if (0 != process_granule (gr, clt, meta, output_file))
     goto return_status;

   status = 0;
return_status:
   if (gr) gr->granule_close (gr);
   if (clt) clt->clt_delete (clt);
   tio_meta_close (meta);

   return status;
}
