#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "ccd.h"
#include "pixelqf.h"
#include "bpix.h"
#include "granule.h"
#include "instr.h"
#include "dark.h"
#include "sensorcal.h"
#include "output.h"

#include "control.h"
#include "process.h"

typedef struct
{
   Granule_Exprec_Type *exprec;
   Image_Type *img_err;
   float storage_region_dark[4];
   float ccd_temp;
   int index;
}
Exprec_Meta_Type;

typedef struct
{
   int saturated_neighbor_hw_serial;
   int saturated_neighbor_hw_parallel;
   double bpix_update_thresh;
   int bpix_update_num_exprecs_needed;
}
Process_Control_Type;

static int get_control_params (config_t *cfg, Process_Control_Type *pct)
{
   config_setting_t *s, *sub;

   if (NULL == (s = config_lookup (cfg, "pixel_quality_flag_params")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing pqf_params in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (sub = config_setting_get_member (s, "saturation")))
       || (CONFIG_TRUE != config_setting_lookup_int (sub, "neighbor_hw_serial", &pct->saturated_neighbor_hw_serial))
       || (CONFIG_TRUE != config_setting_lookup_int (sub, "neighbor_hw_parallel", &pct->saturated_neighbor_hw_parallel))
       )
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading saturation flag parameters",
                     __func__);
        return -1;
     }

   if ((NULL == (sub = config_setting_get_member (s, "badpix_update")))
        || (CONFIG_TRUE != config_setting_lookup_float (sub, "threshold", &pct->bpix_update_thresh))
        || (CONFIG_TRUE != config_setting_lookup_int (sub, "num_exprecs_needed", &pct->bpix_update_num_exprecs_needed))
       )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading badpix_update params in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static void free_exprec_array (Exprec_Meta_Type *a, int num_exprecs,
                               Granule_Type *gr)
{
   int i;

   if ((gr == NULL) || (a == NULL))
     return;

   for (i = 0; i < num_exprecs; i++)
     {
        Exprec_Meta_Type *a_i = &a[i];
        gr->granule_free_exprec (a_i->exprec);
        image_free (a_i->img_err);
     }
   FREE(a);
}

static Exprec_Meta_Type *alloc_exprec_array (int num_exprecs)
{
   Exprec_Meta_Type *a = NULL;
   size_t exprec_array_size;

   exprec_array_size = num_exprecs * sizeof(*a);
   if (NULL == (a = (Exprec_Meta_Type *) MALLOC (exprec_array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)a, 0, exprec_array_size);

   return a;
}

static int validate_exposure_type (int exposure_type0, int exposure_type)
{
   if (exposure_type == EXPREC_TYPE_UNKNOWN)
     {
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains exposures of unknown type",
                     __func__);
        return -1;
     }

   if (exposure_type != exposure_type0)
     {
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains multiple exposure types",
                     __func__);
        return -1;
     }

   return 0;
}

static Dark_Array_Type *create_dark_array (Exprec_Meta_Type *exprec_array,
                                           int num_exprecs)
{
   Dark_Array_Type *dark_array = NULL;
   int i, k;

   if (NULL == (dark_array = dark_array_alloc (num_exprecs)))
     return NULL;

   for (i = 0; i < num_exprecs; i++)
     {
        Exprec_Meta_Type *xr = &exprec_array[i];
        Granule_Exprec_Type *exprec = xr->exprec;
        double sdc;
        sdc = 0.0;
        for (k = 0; k < 4; k++)
          {
             sdc += xr->storage_region_dark[k];
          }
        sdc /= 4;
        if (0 != dark_array_elem_set (dark_array, i, exprec->img, sdc, xr->ccd_temp,
                                      exprec->exposure_time))
          {
             dark_array_free (dark_array);
             return NULL;
          }
     }

   return dark_array;
}

static int make_dark_table (config_t *cfg, int is_linearity,
                            Exprec_Meta_Type *exprec_array, int num_exprecs,
                            const char *output_file)
{
   Dark_Table_Type *dtt = NULL;
   Dark_Config_Type *dcfg = NULL;
   Dark_Array_Type *dark_array = NULL;

   if (NULL == (dcfg = dark_table_config (cfg, is_linearity)))
     return -1;

   if (NULL == (dark_array = create_dark_array (exprec_array, num_exprecs)))
     goto return_error;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "creating dark current table: %s",
              output_file);
   if (NULL == (dtt = dark_table_create (dcfg, dark_array)))
     goto return_error;

   if (0 != dtt->dtt_write (dtt, output_file))
     goto return_error;

   dark_table_config_free (dcfg);
   dark_array_free (dark_array);
   dtt->dtt_delete (dtt);

   return 0;
return_error:
   dark_table_config_free (dcfg);
   dark_array_free (dark_array);
   if (dtt) dtt->dtt_delete (dtt);
   return -1;
}

static int compute_current_and_trim (const CCD_Type *ccd,
                                     const Instr_Type *instr,
                                     const Pixelqf_Type *pt,
                                     const Process_Control_Type *pct,
                                     Exprec_Meta_Type *xr)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   float *mean_sdc = xr->storage_region_dark;
   Image_Type *aimg = NULL;
   double smear_fraction;
   double exposure_time_per_frame;
   int i;

   /* FIXME: is this the correct temp? */
   if (0 != instr->instr_ccd_temp1 (instr, exprec->start_time, &xr->ccd_temp))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "%s: ccd_temp lookup failed, start_time=%15.12e",
                   __func__, exprec->start_time);
        /* drop */
     }

   if (-1 == ccd->ccd_correct_coadd (ccd, exprec->num_coadds, exprec->img))
     return -1;
   exposure_time_per_frame = exprec->exposure_time / exprec->num_coadds;

   if (-1 == ccd->ccd_correct_offset (ccd, exprec->img))
     return -1;

   if (0 == EXPREC_TYPE_IS_LINEARITY(exprec->exposure_type))
     {
        if (-1 == ccd->ccd_correct_nonlinearity (ccd, exprec->img))
          return -1;
     }

   if (-1 == ccd->ccd_correct_gain (ccd, exprec->img))
     return -1;

   smear_fraction = (exprec->frame_transfer_time
                     /(exprec->frame_transfer_time + exposure_time_per_frame));

   if (-1 == ccd->ccd_correct_smear (ccd, &smear_fraction, exprec->img))
     return -1;

   if (-1 == ccd->ccd_mean_storage_region_dark (ccd, exprec->img, mean_sdc))
     return -1;

   if (0) fprintf (stderr, "mean sdc:  %7.1f %7.1f %7.1f %7.1f\n",
                   mean_sdc[0], mean_sdc[1], mean_sdc[2], mean_sdc[3]);

   if (NULL == (aimg = ccd->ccd_copy_active_pixels (ccd, exprec->img)))
     return -1;
   image_free (exprec->img);
   exprec->img = aimg;

   if (EXPREC_TYPE_IS_DARK(exprec->exposure_type))
     {
        if (-1 == pt->pqf_flag_hotcold (pt, exprec->img))
          return -1;
     }

   if (-1 == pt->pqf_flag_neighbor (pt, exprec->img,
                                    pct->saturated_neighbor_hw_serial,
                                    pct->saturated_neighbor_hw_parallel,
                                    IMAGE_PQF_SATURATED, IMAGE_PQF_SATURATED))
     return -1;

   /* FIXME: why not compute an uncertainty for everything,
    * include dark images? */
   if ((exprec->exposure_type == EXPREC_TYPE_RADIANCE)
       || (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type)))
     {
        Image_Type *noisesq = NULL;

        if (NULL == (noisesq = image_dup (exprec->img)))
          return -1;

        if (-1 == ccd->ccd_update_noisesq (ccd, mean_sdc, noisesq))
          {
             image_free (noisesq);
             return -1;
          }

        image_sqrt (noisesq);
        xr->img_err = noisesq;

        image_scale (xr->img_err, 1.0/exposure_time_per_frame);
     }

   image_scale (exprec->img, 1.0/exposure_time_per_frame);

   for (i = 0; i < 4; i++)
     {
        xr->storage_region_dark[i] /= exprec->readout_time;
     }

   return 0;
}

static int process_dark (config_t *cfg, const Control_Type *ctrl,
                         Process_Control_Type *pct,
                         Granule_Type *gr)
{
   CCD_Type *ccd = NULL;
   Granule_Exprec_Type *exprec = NULL;
   Instr_Type *instr = NULL;
   Exprec_Meta_Type *exprec_array = NULL;
   Pixelqf_Type *pt = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Badpix_Map_Occur_Type *bpix_occur = NULL;
   Badpix_Bitmap_Type bpix_occur_mask;
   int ixr, num_exprecs, exposure_type;
   int bpix_occur_threshold;
   int is_linearity, status = -1;

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   if (NULL == (ccd = ccd_init (cfg)))
     goto return_status;

   if (NULL == (pt = pixelqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file, ctrl->instr_glob,
                                    gr->granule_tstart(gr), gr->granule_tend(gr))))
     goto return_status;

   num_exprecs = gr->granule_num_exprecs(gr);
   if (ctrl->limit_num_granules < num_exprecs)
     num_exprecs = ctrl->limit_num_granules;

   if (NULL == (exprec_array = alloc_exprec_array (num_exprecs)))
     goto return_status;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;

   bpix_occur_mask = IMAGE_PQF_HOT_PIXEL | IMAGE_PQF_COLD_PIXEL;
   bpix_occur = bpix_occur_open (bpixmap->num_rows, bpixmap->num_cols,
                                 bpix_occur_mask);
   if (NULL == bpix_occur)
     goto return_status;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "Converting DN to e-/s:");
   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        Exprec_Meta_Type *xr = &exprec_array[ixr];

        tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec %3d/%d", ixr, num_exprecs);

        if (NULL == (exprec = gr->granule_read_exprec_by_index (gr, ixr, NULL)))
          goto return_status;

        xr->exprec = exprec;
        xr->index = ixr;

        if (-1 == validate_exposure_type (exposure_type, exprec->exposure_type))
          goto return_status;

        if (0 != ccd->ccd_apply_pixel_quality_flags (ccd, exprec->img,
                                                     bpixmap->bits, bpixmap->num_rows, bpixmap->num_cols))
          goto return_status;

        if (-1 == compute_current_and_trim (ccd, instr, pt, pct, xr))
          goto return_status;

        if (-1 == bpix_occur_incr (bpix_occur,
                                   exprec->img->pixel_quality_flags))
          goto return_status;
     }

   bpix_occur_threshold = pct->bpix_update_thresh * num_exprecs;
   if (num_exprecs > pct->bpix_update_num_exprecs_needed)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "updating bpix map");
        if (0 != bpix_occur_set (bpix_occur, bpix_occur_threshold,
                                 bpix_occur_mask, bpixmap->bits))
          goto return_status;
     }
   /* FIXME: write out updated badpix map */

   is_linearity = (exposure_type == EXPREC_TYPE_LIN_DARK);
   if (0 != make_dark_table (cfg, is_linearity, exprec_array, num_exprecs, ctrl->output_file))
     goto return_status;

   status = 0;
return_status:

   bpix_free (bpixmap);
   bpix_occur_close (bpix_occur);
   free_exprec_array (exprec_array, num_exprecs, gr);

   if (ccd) ccd->ccd_delete (ccd);
   if (instr) instr->instr_delete (instr);
   if (pt) pt->pqf_delete (pt);

   return status;
}

static int subtract_dark_current_img (Image_Type *img, const Image_Type *dc)
{
   Image_Pixel_Type *img_pixels = img->pixels;
   Image_Pixel_Type *dc_pixels = dc->pixels;
   Image_Pqf_Bitmap_Type *img_pqf = img->pixel_quality_flags;
   Image_Pqf_Bitmap_Type *dc_pqf = dc->pixel_quality_flags;
   int i, n = img->num_rows * img->num_cols;
   Image_Pixel_Type dc_pixels_i;

   for (i = 0; i < n; i++)
     {
        if (img_pixels[i] == IMAGE_PIXEL_FILL_VALUE)
          continue;
        dc_pixels_i = dc_pixels[i];
        if ((dc_pqf[i] == 0) && (dc_pixels_i > 0))
          {
             img_pixels[i] -= dc_pixels_i;
          }
        else img_pqf[i] |= IMAGE_PQF_DARK_CORR_ERROR;
      }

   return 0;
}

static int subtract_dark_current (Granule_Exprec_Type *exprec,
                                  const Dark_Table_Type *dtt,
                                  const float *storage_region_dark,
                                  float ccd_temp, Image_Type *tmp_img)
{
   int k, ordered_by = dtt->dtt_ordering (dtt);
   double key;

   switch (ordered_by)
     {
      case DARK_TABLE_ORDERED_BY_SDC:
        key = 0.0;
        for (k = 0; k < 4; k++)
          {
             key += storage_region_dark[k];
          }
        key /= 4;
        break;

      case DARK_TABLE_ORDERED_BY_TEMP:
        key = ccd_temp;
        break;

      case DARK_TABLE_ORDERED_BY_EXPTIME:
        key = exprec->exposure_time;
        break;
     }

   /* FIXME - may not need to interpolate a new dark image
    * for _every_ exposure record */
   if (0 != dtt->dtt_interp (dtt, key, tmp_img))
     return -1;

   return subtract_dark_current_img (exprec->img, tmp_img);
}

static int validate_dark_table_type (const Dark_Table_Type *dtt,
                                     int exposure_type)
{
   int ordering = dtt->dtt_ordering (dtt);

   if ((exposure_type == EXPREC_TYPE_LIN_IRR)
       && (ordering != DARK_TABLE_ORDERED_BY_EXPTIME))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "Unexpected dark interpolation method -- linearity should use 'exptime'");
     }

   return 0;
}

static int radiometric_correction (const Calibration_Type *cal, const Dark_Table_Type *dtt,
                                   Exprec_Meta_Type *xr, Image_Type *tmp_img)
{
   Granule_Exprec_Type *exprec = xr->exprec;

   /* FIXME? propagate dark-current subtraction error to
    * uncertainty estimate? */
   if (0 != subtract_dark_current (exprec, dtt, xr->storage_region_dark,
                                   xr->ccd_temp, tmp_img))
     return -1;

   /* FIXME: update uncertainty? */
   if (0 != cal->cal_apply_prnu (cal, exprec->img))
     return -1;

   /* >>> Stray light correction goes here <<< */

   /* Multiplicative factor converts e/s to photons/s */
   if ((0 != cal->cal_apply_rcoeffs (cal, exprec->img))
       || (0 != cal->cal_apply_rcoeffs (cal, xr->img_err)))
     return -1;

   if (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type))
     {
        double solar_phi=0.0, solar_theta=0.0;
        /* FIXME: placeholders for angles to be computed based on
         * the observation date/time and orbital ephemeris */
        if ((0 != cal->cal_apply_btdf (cal, solar_phi, solar_theta, exprec->img))
            || (0 != cal->cal_apply_btdf (cal, solar_phi, solar_theta, xr->img_err))) /* FIXME: ok? */
          return -1;
     }

   return 0;
}

static Spectral_Data_Type *
finalize_band (const Calibration_Type *cal,
               const Exprec_Meta_Type *xr, int band_id)
{
   Image_Type *img = xr->exprec->img;
   Image_Type *img_err = xr->img_err;
   Spectral_Data_Type *sdt = NULL;
   int status;

   /* FIXME: It might be slightly more efficient to allocate two
    * Spectral_Data_Type objects at a high level and then re-use them */
   if (NULL == (sdt = sdt_extract_band (cal, band_id, img, img_err)))
     return NULL;

   if (0 != cal->cal_wavecal (cal, band_id, sdt->num_xtrack,
                              sdt->img, sdt->img_err, sdt->wave))
     goto return_status;

   status = 0;
return_status:
   if (status)
     {
        sdt_free (sdt);
        return NULL;
     }

   return sdt;
}

static int apply_cal_then_output (Output_Type *out, Calibration_Type *cal,
                                  config_t *cfg, Dark_Table_Type *dtt,
                                  Exprec_Meta_Type *xr,
                                  Image_Type *tmp_img)
{
   Output_Exprec_Type outrec = {0};
   int status = -1;

   (void) cfg;

   if (0 != radiometric_correction (cal, dtt, xr, tmp_img))
     return -1;

   if (NULL == (outrec.uv = finalize_band (cal, xr, TEMPO_BAND_UV)))
     return -1;

   if (NULL == (outrec.vis = finalize_band (cal, xr, TEMPO_BAND_VIS)))
     goto return_status;

   outrec.meta.start_time = xr->exprec->start_time;
   outrec.meta.exposure_time = xr->exprec->exposure_time;
   outrec.meta.mirror_step = xr->exprec->curr_mirror_step;

   if (0 != out->out_write_rec (out, xr->index, &outrec))
     goto return_status;

   status = 0;
return_status:
   sdt_free (outrec.uv);
   sdt_free (outrec.vis);

   return status;
}

static Exprec_Meta_Type *new_exprec_meta_type (void)
{
   return alloc_exprec_array (1);
}

static void free_exprec_meta_type (Exprec_Meta_Type *xr, Granule_Type *gr)
{
   free_exprec_array (xr, 1, gr);
}

#define QUEUE_DEPTH  3
typedef struct
{
   Exprec_Meta_Type *items[QUEUE_DEPTH];
   int num_queued;
}
Queue_Type;

static void queue_init (Queue_Type *q)
{
   memset ((char *)q->items, 0, sizeof (q->items));
   q->num_queued = 0;
}

static Exprec_Meta_Type *queue_shift (Queue_Type *q)
{
   Exprec_Meta_Type *oldest = q->items[0];
   int i;

   for (i = 1; i < QUEUE_DEPTH; i++)
     {
        q->items[i-1] = q->items[i];
     }

   return oldest;
}

static Exprec_Meta_Type *queue_push (Queue_Type *q,
                                     Exprec_Meta_Type *xr)
{
   Exprec_Meta_Type *oldest = queue_shift (q);
   q->items[QUEUE_DEPTH-1] = xr;
   q->num_queued += 1;
   if (q->num_queued > QUEUE_DEPTH)
     q->num_queued = QUEUE_DEPTH;
   return oldest;
}

static void queue_empty (Queue_Type *q, Granule_Type *gr)
{
   int i;
   if ((q == NULL) || (gr == NULL))
     return;

   for (i = 0; i < q->num_queued; i++)
     {
        free_exprec_meta_type (q->items[i], gr);
     }
}

static void make_transient_ref_img (const Image_Type *prev,
                                    const Image_Type *next,
                                    Image_Type *img_ref)
{
   int num_rows = img_ref->num_rows;
   int num_cols = img_ref->num_cols;
   int s, p;

   for (p = 0; p < num_rows; p++)
     {
        Image_Pixel_Type *prev_pixels = prev->pixels + p * num_cols;
        Image_Pixel_Type *next_pixels = next->pixels + p * num_cols;
        Image_Pixel_Type *ref_pixels = img_ref->pixels + p * num_cols;
        for (s = 0; s < num_cols; s++)
          {
             if ((prev_pixels[s] != IMAGE_PIXEL_FILL_VALUE)
                 && (next_pixels[s] != IMAGE_PIXEL_FILL_VALUE))
               {
                  Image_Pixel_Type pxp = prev_pixels[s] * next_pixels[s];
                  ref_pixels[s] = (pxp >= 0) ? sqrt(pxp) : IMAGE_PIXEL_FILL_VALUE;
               }
             else ref_pixels[s] = IMAGE_PIXEL_FILL_VALUE;
          }
     }
}

static int flag_transients1 (const Pixelqf_Type *pt,
                             const Badpix_Map_Type *bpixmap,
                             int num_exprecs, Queue_Type *q, int exprec_index,
                             Image_Type *img_ref)
{
   Exprec_Meta_Type *prev, *xr, *next;
   int status = 0;

   if (img_ref == NULL)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: img_ref = NULL!", __func__);
        return -1;
     }

   /* The queue lets us lag one record behind so that
    * when exprec_index=i, we process exprec_index=i-1.
    * On the final pass, we process two records, to
    * ensure that every record is processed in the end.
    */

   if (exprec_index == 0)
     {
        /* queue contains {NULL,NULL,img(0)},
         * nothing to do yet */
        return 0;
     }

   /* When we don't need the scratch space at img_ref,
    * we over-write the local copy of the img_Ref pointer so
    * that it points to the pre-existing object that's needed.
    */

   prev = q->items[0];
   xr   = q->items[1];
   next = q->items[2];

   if (prev == NULL)
     {
        /* e.g. exprec_index == 1 for num_exprecs >= 2 */
        img_ref = next->exprec->img;
     }
   else
     {
        /* here's where we actually need the space img_ref points to */
        make_transient_ref_img (prev->exprec->img, next->exprec->img,
                                img_ref);
     }

   status = pt->pqf_flag_transients (pt, bpixmap->bits, img_ref,
                                     xr->exprec->img);
   if (status != 0)
     return status;

   /* When exprec_index==num_exprecs-1, process one more to finish: */
   if (exprec_index == num_exprecs-1)
     {
        /* queue contains {..., img(n-2), img(n-1)} */
        prev = q->items[1];
        img_ref = prev->exprec->img;
        xr = q->items[2];
        status = pt->pqf_flag_transients (pt, bpixmap->bits, img_ref,
                                          xr->exprec->img);
     }

   return status;
}

static int process_exposure (config_t *cfg, const Control_Type *ctrl,
                             Process_Control_Type *pct,
                             Granule_Type *gr)
{
   Queue_Type exprec_queue;
   CCD_Type *ccd = NULL;
   Granule_Exprec_Type *exprec = NULL;
   Instr_Type *instr = NULL;
   Calibration_Type *cal = NULL;
   Exprec_Meta_Type *xr = NULL;
   Exprec_Meta_Type *xr_ready = NULL;
   Pixelqf_Type *pt = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Dark_Table_Type *dtt = NULL;
   Output_Type *out = NULL;
   Image_Type *tmp_img = NULL;
   int num_serial_active_full, num_parallel_active_full;
   int ixr, num_exprecs, exposure_type;
   int status = -1;

   queue_init (&exprec_queue);

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   if (NULL == (ccd = ccd_init (cfg)))
     goto return_status;

   if (NULL == (cal = sensorcal_init (cfg)))
     return -1;

   if (NULL == (pt = pixelqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file, ctrl->instr_glob,
                                    gr->granule_tstart(gr), gr->granule_tend(gr))))
     goto return_status;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;

   if (NULL == (dtt = dark_table_read (ctrl->dark_file)))
     goto return_status;
   if (0 != validate_dark_table_type (dtt, exposure_type))
     goto return_status;

   num_exprecs = gr->granule_num_exprecs(gr);
   if (ctrl->limit_num_granules < num_exprecs)
     num_exprecs = ctrl->limit_num_granules;

   /* Allocate reusable scratch space to hold a full, trimmed CCD image */
   ccd->ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);
   if (NULL == (tmp_img = image_new (num_parallel_active_full, num_serial_active_full)))
     goto return_status;

   /* Dimension the output file to hold num_exprecs trimmed frames,
    * split into two wavelength bands, with wavelength the fastest
    * varying dimension. */
   if ((NULL == (out = output_alloc (cfg, exposure_type)))
       || (0 != out->out_set_file (out, ctrl->output_file)))
     goto return_status;
   out->out_set_dims (out, num_exprecs,
                      num_serial_active_full, num_parallel_active_full/2);
   if (0 != out->out_create (out))
     goto return_status;

   if (exposure_type == EXPREC_TYPE_RADIANCE)
     {
        int ncid_from = gr->granule_ncid (gr);
        int ncid_to = out->out_ncid (out);
        if (0 != TIO_copy_granule_ident (ncid_from, ncid_to))
          goto return_status;
     }

   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        Exprec_Meta_Type *xr_to_delete;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec %3d/%d", ixr, num_exprecs);

        if (NULL == (exprec = gr->granule_read_exprec_by_index (gr, ixr, NULL)))
          goto return_status;
        if (-1 == validate_exposure_type (exposure_type, exprec->exposure_type))
          goto return_status;
        if (0 != ccd->ccd_apply_pixel_quality_flags (ccd, exprec->img,
                                                     bpixmap->bits, bpixmap->num_rows, bpixmap->num_cols))
          goto return_status;

        if (NULL == (xr = new_exprec_meta_type ()))
          goto return_status;
        xr->exprec = exprec;
        xr->index = ixr;

        if (-1 == compute_current_and_trim (ccd, instr, pt, pct, xr))
          goto return_status;

        /* We need at least 2 frames queued to look for transients */

        xr_to_delete = queue_push (&exprec_queue, xr);
        free_exprec_meta_type (xr_to_delete, gr);
        if (exprec_queue.num_queued < 2)
          continue;

        if (0 != flag_transients1 (pt, bpixmap, num_exprecs,
                                   &exprec_queue, ixr, tmp_img))
          goto return_status;

        /* Frame ixr-1 is now ready to continue processing */
        xr_ready = exprec_queue.items[1];
        if (0 != apply_cal_then_output (out, cal, cfg, dtt, xr_ready, tmp_img))
          goto return_status;
     }

   /* Process the last entry in the queue, exprec[num_exprecs-1] */
   xr_ready = exprec_queue.items[2];
   if (0 != apply_cal_then_output (out, cal, cfg, dtt, xr_ready, tmp_img))
     goto return_status;

   status = 0;
return_status:

   queue_empty (&exprec_queue, gr);
   image_free (tmp_img);
   bpix_free (bpixmap);

   if (ccd) ccd->ccd_delete (ccd);
   if (instr) instr->instr_delete (instr);
   if (pt) pt->pqf_delete (pt);
   if (cal) cal->cal_delete (cal);
   if (dtt) dtt->dtt_delete (dtt);
   if (out)
     {
        (void) out->out_close (out);
        out->out_free (out);
     }

   return status;
}

int process_inputs (config_t *cfg, const Control_Type *ctrl)
{
   Granule_Type *gr = NULL;
   Process_Control_Type pct = {0};
   int exposure_type, status = -1;

   if (0 != get_control_params (cfg, &pct))
     goto return_status;

   if (NULL == (gr = granule_open (ctrl->input_file)))
     goto return_status;

   if (0 != gr->granule_type (gr, &exposure_type))
     goto return_status;

   switch (exposure_type)
     {
      case EXPREC_TYPE_DARK:
      case EXPREC_TYPE_LIN_DARK:
        status = process_dark (cfg, ctrl, &pct, gr);
        break;

      case EXPREC_TYPE_RADIANCE:
      case EXPREC_TYPE_IRRADIANCE:
      case EXPREC_TYPE_LIN_IRR:
        status = process_exposure (cfg, ctrl, &pct, gr);
        break;

      default:
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains exposures of unknown or unexpected type",
                     __func__);
        break;
     }

return_status:

   if (gr) gr->granule_close (gr);

   return status;
}
