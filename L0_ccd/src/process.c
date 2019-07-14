#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <wordexp.h>

#include <libconfig.h>
#include <libnovas.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

/* from plan/src */
#include <solar.h>

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
#include "util.h"

typedef struct
{
   Granule_Exprec_Type *exprec;
   Image_Type *img_err;
   float storage_region_dark[4];
   float fpa_temp;
   float fpe_temp;
   /* ConOps 3.3: instrument command parameters: NUM_TG_ROWS, NUM_DG_ROWS */
   int num_dg_rows;   /* index of first row included in storage region dark summation */
   int num_tg_rows;   /* number of rows included in storage region dark summation */
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

typedef struct
{
   char *ephem_name;
   double jd_begin;
   double jd_end;
   short int de_number;
}
Ephem_Type;

typedef struct
{
   short int year;
   short int month;
   short int day;
   double hour;
}
Cal_Date_Type;

static int Write_Nominal_Wavelength_Grid;

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

static void free_exprec_meta (Exprec_Meta_Type *xr, Granule_Type *gr)
{
   if (xr == NULL)
     return;
   gr->granule_free_exprec (xr->exprec);
   image_free (xr->img_err);
   memset ((char *)xr, 0, sizeof (*xr));
   FREE(xr);
}

static Exprec_Meta_Type *alloc_exprec_meta (void)
{
   Exprec_Meta_Type *xr = NULL;

   if (NULL == (xr = (Exprec_Meta_Type *) MALLOC (sizeof *xr)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)xr, 0, sizeof *xr);
   return xr;
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

static int compute_noisesq_for_active_pixels (const CCD_Type *ccd, const Image_Type *img,
                                              Exprec_Meta_Type *xr)
{
   Image_Type *noisesq = NULL;
   Image_Type *aimg = NULL;
   int status = -1;

   if (NULL == (aimg = ccd->ccd_copy_active_pixels (ccd, img)))
     return -1;
   if (NULL == (noisesq = image_dup (aimg)))
     goto return_status;

   if (-1 == ccd->ccd_update_noisesq (ccd, xr->storage_region_dark, noisesq))
     goto return_status;

   xr->img_err = noisesq;

   status = 0;
return_status:
   image_free(aimg);
   if (status)
     {
        image_free (noisesq);
     }

   return status;
}

static int compute_current_and_trim (CCD_Type *ccd,
                                     const Instr_Type *instr,
                                     const Pixelqf_Type *pqft,
                                     const Process_Control_Type *pct,
                                     Exprec_Meta_Type *xr)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   Image_Type *aimg = NULL;
   double smear_fraction;
   double exposure_time_per_frame;
   int i;

   if (-1 == ccd->ccd_correct_coadd (ccd, exprec->num_coadds, exprec->img))
     return -1;
   exposure_time_per_frame = exprec->exposure_time / exprec->num_coadds;

   if (0 != ccd->ccd_configure_using_octant_phase (ccd, exprec->img))
     return -1;

   if (0 != ccd->ccd_correct_offset (ccd, exprec->img))
     return -1;

   if (0) (void) image_write_raw (exprec->img, "offset");

   if (0 == EXPREC_TYPE_IS_LINEARITY(exprec->exposure_type))
     {
        if (-1 == ccd->ccd_correct_nonlinearity (ccd, exprec->img))
          return -1;
     }

   if (0 != ccd->ccd_correct_crosstalk (ccd, exprec->img))
     return -1;

   if (0 != instr->instr_temps (instr, exprec->start_time, &xr->fpa_temp, &xr->fpe_temp))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "%s: temperature lookup failed, start_time=%15.12e",
                   __func__, exprec->start_time);
        /* drop */
     }

   /* convert DN to electrons */
   if (-1 == ccd->ccd_correct_gain (ccd, exprec->img, xr->fpa_temp, xr->fpe_temp))
     return -1;

   if (0) (void) image_write_raw (exprec->img, "gain");

   if (-1 == ccd->ccd_mean_storage_region_dark (ccd, exprec->img,
                                                xr->num_dg_rows, xr->num_tg_rows,
                                                xr->storage_region_dark))
     {
        return -1;
     }

   if (0) fprintf (stderr, "mean sdc:  %7.1f %7.1f %7.1f %7.1f\n",
                   xr->storage_region_dark[0],
                   xr->storage_region_dark[1],
                   xr->storage_region_dark[2],
                   xr->storage_region_dark[3]);

   /* Compute noisesq before smear correction */
   if ((exprec->exposure_type == EXPREC_TYPE_RAD)
       || (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type)))
     {
        if (0 != compute_noisesq_for_active_pixels (ccd, exprec->img, xr))
          return -1;
     }

   smear_fraction = (exprec->frame_transfer_time
                     /(exprec->frame_transfer_time + exposure_time_per_frame));

   if (-1 == ccd->ccd_correct_smear (ccd, &smear_fraction, exprec->img))
     return -1;

   /* trim parallel overclocks, and serial leading and trailing */
   if (NULL == (aimg = ccd->ccd_copy_active_pixels (ccd, exprec->img)))
     return -1;
   image_free (exprec->img);
   exprec->img = aimg;

   if (-1 == pqft->pqf_flag_neighbor (pqft, exprec->img,
                                      pct->saturated_neighbor_hw_serial,
                                      pct->saturated_neighbor_hw_parallel,
                                      IMAGE_PQF_SATURATED, IMAGE_PQF_SATURATED))
     {
        return -1;
     }

   /* Compute pixel current: electrons/sec */
   image_scale (exprec->img, 1.0/exposure_time_per_frame);
   if (xr->img_err)
     {
        image_scale (xr->img_err, 1.0/exposure_time_per_frame);
     }
   for (i = 0; i < 4; i++)
     {
        xr->storage_region_dark[i] /= exprec->readout_time;
     }

   if (0 != ccd->ccd_correct_prnu (ccd, exprec->img))
     return -1;

   if (EXPREC_TYPE_IS_DARK(exprec->exposure_type))
     {
        if (-1 == pqft->pqf_flag_hotcold (pqft, exprec->img))
          return -1;
     }

   if (0) (void) image_write_raw (exprec->img, "prnu");

   return 0;
}

static int write_dark_exprec (int ncid, const Exprec_Meta_Type *xr)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   Image_Type *img = exprec->img;
   int start[3], count[3];

   start[0] = xr->index;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;

   if ((0 != TIO_put_var_section (ncid, "image_start_time", start, count, TIO_DOUBLE, &exprec->start_time))
       || (0 != TIO_put_var_section (ncid, "fpa_temp", start, count, TIO_FLOAT, &xr->fpa_temp))
       || (0 != TIO_put_var_section (ncid, "fpe_temp", start, count, TIO_FLOAT, &xr->fpe_temp))
       || (0 != TIO_put_var_section (ncid, TEMPO_VAR_EXPOSURE_TIME, start, count, TIO_DOUBLE, &exprec->exposure_time)))
     return -1;

   count[0] = 1;
   count[1] = 4;
   if (0 != TIO_put_var_section (ncid, "mean_sdc", start, count, TIO_FLOAT, xr->storage_region_dark))
     return -1;

   count[0] = 1;
   count[1] = img->num_rows;
   count[2] = img->num_cols;

   if (0 != TIO_put_var_section (ncid, "image", start, count, TIO_FLOAT, img->pixels))
     return -1;
   if (0 != TIO_put_var_section (ncid, TEMPO_VAR_PQF, start, count, TIO_USHORT, img->pixel_quality_flags))
     return -1;

   return 0;
}

static int process_dark (config_t *cfg, const Control_Type *ctrl,
                         Process_Control_Type *pct,
                         Granule_Type *gr, TIO_Meta_Type *meta)
{
   CCD_Type *ccd = NULL;
   Instr_Type *instr = NULL;
   Exprec_Meta_Type *xr = NULL;
   Pixelqf_Type *pqft = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Badpix_Map_Occur_Type *bpix_occur = NULL;
   Badpix_Bitmap_Type bpix_occur_mask;
   int ixr, num_exprecs, exposure_type;
   int num_parallel_active_full, num_serial_active_full, drk_ncid = 0;
   int bpix_occur_threshold;
   int status = -1;

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   if (NULL == (ccd = ccd_init (cfg, meta)))
     goto return_status;

   if (NULL == (pqft = pixelqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file, ctrl->instr_glob,
                                    gr->granule_tstart(gr), gr->granule_tend(gr), meta)))
     goto return_status;

   num_exprecs = gr->granule_num_exprecs(gr);
   if (ctrl->limit_num_granules < num_exprecs)
     num_exprecs = ctrl->limit_num_granules;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;
   if (0 != meta_record_basename (meta, ctrl->bpix_file))
     goto return_status;

   bpix_occur_mask = IMAGE_PQF_HOT_PIXEL | IMAGE_PQF_COLD_PIXEL;
   bpix_occur = bpix_occur_open (bpixmap->num_rows, bpixmap->num_cols,
                                 bpix_occur_mask);
   if (NULL == bpix_occur)
     goto return_status;

   /* Open the output file */
   tell_vlog (TELL_MSGTYPE_INFO, 1, "Opening output file: %s", ctrl->output_file);
   if (0 != TIO_create (ctrl->output_file, NC_NETCDF4, &drk_ncid))
     goto return_status;
   ccd->ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);
   if (0 != drk_create_file (drk_ncid, num_exprecs, num_parallel_active_full, num_serial_active_full))
     goto return_status;
   if (0 != tio_write_epoch_timestamp (drk_ncid, NC_GLOBAL))
     goto return_status;
   if (0 != TIO_copy_granule_ident (gr->granule_ncid(gr), drk_ncid))
     goto return_status;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "Converting DN to e-/s:");
   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec %3d/%d", ixr, num_exprecs);

        if (NULL == (xr = alloc_exprec_meta ()))
          goto return_status;

        if (NULL == (xr->exprec = gr->granule_read_exprec_by_index (gr, ixr, NULL)))
          goto return_status;

        xr->index = ixr;

        if (-1 == validate_exposure_type (exposure_type, xr->exprec->exposure_type))
          goto return_status;

        if (0 != ccd->ccd_apply_pixel_quality_flags (ccd, xr->exprec->img,
                                                     bpixmap->bits, bpixmap->num_rows, bpixmap->num_cols))
          goto return_status;

        if (-1 == compute_current_and_trim (ccd, instr, pqft, pct, xr))
          goto return_status;

        if (0 != write_dark_exprec (drk_ncid, xr))
          goto return_status;

        if (-1 == bpix_occur_incr (bpix_occur, xr->exprec->img->pixel_quality_flags))
          goto return_status;

        free_exprec_meta (xr, gr);
        xr = NULL;
     }

   bpix_occur_threshold = pct->bpix_update_thresh * num_exprecs;
   if (num_exprecs > pct->bpix_update_num_exprecs_needed)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "updating internal bpix map");
        if (0 != bpix_occur_set (bpix_occur, bpix_occur_threshold,
                                 bpix_occur_mask, bpixmap->bits))
          goto return_status;
     }

   /* FIXME: eventually, we may write out an updated badpix map */

   status = 0;
return_status:

   bpix_free (bpixmap);
   bpix_occur_close (bpix_occur);
   free_exprec_meta (xr, gr);
   xr = NULL;

   if (drk_ncid)
     {
        if (0 != TIO_close (drk_ncid))
          {
             tell_verror (TELL_IO_ERROR, "Closing output file: %s", ctrl->output_file);
             if (status == 0) status = -1;
          }
        else tell_vlog (TELL_MSGTYPE_INFO, 1, "Closed output file: %s", ctrl->output_file);
     }
   if (ccd) ccd->ccd_delete (ccd);
   if (instr) instr->instr_delete (instr);
   if (pqft) pqft->pqf_delete (pqft);

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
        if ((dc_pqf[i] == 0) && (0 <= dc_pixels_i) && (dc_pixels_i <= img_pixels[i]))
          {
             img_pixels[i] -= dc_pixels_i;
          }
        else img_pqf[i] |= IMAGE_PQF_DARK_CORR_ERROR;
      }

   return 0;
}

static int julian_date_from_taix (double taix, double *jd_utc)
{
   int year, month, day;
   double hour;

   if (0 != tio_time_taix_to_utc_caldate (taix, &year, &month, &day, &hour))
     return -1;

   *jd_utc = novas_julian_date ((short int)year, (short int)month, (short int) day, hour);

   return 0;
}

static int radiometric_correction (const Calibration_Type *cal, Solar_Geom_Type *sgt,
                                   const Dark_Type *drk,
                                   Exprec_Meta_Type *xr, Image_Type *tmp_img)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   Dark_Lookup_Type dlt =
     {
        .fpa_temp = xr->fpa_temp,
        .exposure_time = exprec->exposure_time,
        .storage_region_dark = xr->storage_region_dark,
     };

   /* copy the appropriate dark current image into tmp_img */
   if (0 != drk->drk_get_image (drk, &dlt, tmp_img))
     return -1;

   /* subtract the dark current image, leaving the result in place */
   if (0 != subtract_dark_current_img (exprec->img, tmp_img))
     return -1;

   /* >>> Stray light correction goes here <<< */

   /* Multiplicative factor converts e/s to photons/s */
   if ((0 != cal->cal_apply_radcal_coeffs (cal, exprec->img))
       || (0 != cal->cal_apply_radcal_coeffs (cal, xr->img_err)))
     return -1;

   if (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type))
     {
        int use_reference_diffuser = (exprec->exposure_type == EXPREC_TYPE_IRR_REF);
        double jd_utc, solar_theta, solar_phi;

        if ((0 != julian_date_from_taix (exprec->start_time, &jd_utc))
            || (0 != sgt->sgt_sat_sun_angles (sgt, jd_utc, &solar_theta, &solar_phi)))
          return -1;

        if ((0 != cal->cal_apply_btdf (cal, use_reference_diffuser, solar_phi, solar_theta, exprec->img))
            || (0 != cal->cal_apply_btdf (cal, use_reference_diffuser, solar_phi, solar_theta, xr->img_err)))
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
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: extracting band, band_id=%d", __func__, band_id);
        goto return_status;
     }

   if (Write_Nominal_Wavelength_Grid)
     {
        /* only do this once per band in each granule */
        if (0 != cal->cal_nominal_wavelength_grid (cal, band_id, sdt->wave))
          goto return_status;
     }

   status = 0;
return_status:
   if (status)
     {
        sdt_free (sdt);
        return NULL;
     }

   return sdt;
}

static int apply_cal_and_output (Output_Type *out, Calibration_Type *cal, Solar_Geom_Type *sgt,
                                 Dark_Type *drk, Exprec_Meta_Type *xr, Image_Type *tmp_img)
{
   Output_Exprec_Type outrec = {0};
   int num_negative, status = -1;

   if (0 != radiometric_correction (cal, sgt, drk, xr, tmp_img))
     return -1;

   if (0) (void) image_write_raw (xr->exprec->img, "final");

   if ((num_negative = image_check_negative_pixels (xr->exprec->img, 1)) < 0)
     return -1;
   if (num_negative > 0)
     {
	tell_vwarn (0, "%s: set processing error bit in %d pixels with ((value<0) && (pqf==0))",
		    __func__, num_negative);
     }

   if ((NULL == (outrec.uv = finalize_band (cal, xr, TEMPO_BAND_UV)))
       || (NULL == (outrec.vis = finalize_band (cal, xr, TEMPO_BAND_VIS))))
     goto return_status;

   outrec.meta.start_time = xr->exprec->start_time;
   outrec.meta.exposure_time = xr->exprec->exposure_time;
   outrec.meta.mirror_step = xr->exprec->curr_mirror_step;

   outrec.write_nominal_wavelength_grid = Write_Nominal_Wavelength_Grid;

   if (0 != out->out_write_rec (out, xr->index, &outrec))
     goto return_status;

   Write_Nominal_Wavelength_Grid = 0;

   status = 0;
return_status:
   sdt_free (outrec.uv);
   sdt_free (outrec.vis);

   return status;
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
        free_exprec_meta (q->items[i], gr);
        q->items[i] = NULL;
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

static int flag_transients1 (const Pixelqf_Type *pqft,
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

   status = pqft->pqf_flag_transients (pqft, bpixmap->bits, img_ref, xr->exprec->img);
   if (status != 0)
     return status;

   /* When exprec_index==num_exprecs-1, process one more to finish: */
   if (exprec_index == num_exprecs-1)
     {
        /* queue contains {..., img(n-2), img(n-1)} */
        prev = q->items[1];
        img_ref = prev->exprec->img;
        xr = q->items[2];
        status = pqft->pqf_flag_transients (pqft, bpixmap->bits, img_ref, xr->exprec->img);
     }

   return status;
}

static int process_exposure (config_t *cfg, const Control_Type *ctrl,
                             Process_Control_Type *pct,
                             Granule_Type *gr, TIO_Meta_Type *meta)
{
   Queue_Type exprec_queue = {0};
   CCD_Type *ccd = NULL;
   Instr_Type *instr = NULL;
   Calibration_Type *cal = NULL;
   Exprec_Meta_Type *xr = NULL;
   Exprec_Meta_Type *xr_ready = NULL;
   Pixelqf_Type *pqft = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Dark_Type *drk = NULL;
   Output_Type *out = NULL;
   Image_Type *tmp_img = NULL;
   Solar_Geom_Type *sgt = NULL;
   int num_serial_active_full, num_parallel_active_full;
   int ixr, num_exprecs, exposure_type, ncid_from, ncid_to;
   int do_flag_transients;
   int status = -1;

   queue_init (&exprec_queue);

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   if (NULL == (ccd = ccd_init (cfg, meta)))
     goto return_status;

   if (NULL == (cal = sensorcal_init (cfg, meta)))
     return -1;

   if (NULL == (pqft = pixelqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file, ctrl->instr_glob,
                                    gr->granule_tstart(gr), gr->granule_tend(gr), meta)))
     goto return_status;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;
   if (0 != meta_record_basename (meta, ctrl->bpix_file))
     goto return_status;

   if (NULL == (drk = drk_open (ctrl->dark_file)))
     goto return_status;
   if (0 != meta_record_basename (meta, ctrl->dark_file))
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

   ncid_from = gr->granule_ncid (gr);

   if (exposure_type == EXPREC_TYPE_RAD)
     {
        int scan_type;
        ncid_to = out->out_ncid (out);
        if (0 != TIO_copy_granule_ident (ncid_from, ncid_to))
          goto return_status;
        /* Copy the scan_type attribute to the Level 1 file.
         * It's not in the granule ident struct for now, because putting it there would
         * be too much trouble.
         */
        if ((0 != TIO_get_att (ncid_from, NC_GLOBAL, "scan_type", NC_INT, &scan_type))
            ||(0 != TIO_put_att (ncid_to, NC_GLOBAL, "scan_type", NC_INT, 1, &scan_type)))
          goto return_status;
     }
   else if (EXPREC_TYPE_IS_IRRADIANCE(exposure_type))
     {
        if (NULL == (sgt = solar_geom_init (cfg)))
          goto return_status;
     }

   /* FIXME? */
   do_flag_transients = 0;
   /* We need at least 3 frames to look for transients */
   if (num_exprecs < 3) do_flag_transients = 0;

   Write_Nominal_Wavelength_Grid = 1;

   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        Exprec_Meta_Type *xr_to_delete;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec %3d/%d", ixr, num_exprecs);

        if (NULL == (xr = alloc_exprec_meta ()))
          goto return_status;

        xr->index = ixr;

        if (NULL == (xr->exprec = gr->granule_read_exprec_by_index (gr, ixr, NULL)))
          goto return_status;
        if (-1 == validate_exposure_type (exposure_type, xr->exprec->exposure_type))
          goto return_status;
        if (0 != ccd->ccd_apply_pixel_quality_flags (ccd, xr->exprec->img,
                                                     bpixmap->bits, bpixmap->num_rows, bpixmap->num_cols))
          goto return_status;

        if (-1 == compute_current_and_trim (ccd, instr, pqft, pct, xr))
          goto return_status;

        if (do_flag_transients == 0)
          {
             if (0 != apply_cal_and_output (out, cal, sgt, drk, xr, tmp_img))
               goto return_status;
             free_exprec_meta (xr, gr);
             xr = NULL;
          }
        else
          {
             /* We need at least 2 frames queued to look for transients */
             xr_to_delete = queue_push (&exprec_queue, xr);
             free_exprec_meta (xr_to_delete, gr);
             xr_to_delete = NULL;
             if (exprec_queue.num_queued < 2)
               continue;
             if (0 != flag_transients1 (pqft, bpixmap, num_exprecs, &exprec_queue, ixr, tmp_img))
               goto return_status;
             /* Frame ixr-1 is now ready to continue processing */
             xr_ready = exprec_queue.items[1];
             if (0 != apply_cal_and_output (out, cal, sgt, drk, xr_ready, tmp_img))
               goto return_status;
          }
     }

   if (do_flag_transients)
     {
        /* Process the last entry in the queue, exprec[num_exprecs-1] */
        xr_ready = exprec_queue.items[2];
        if (0 != apply_cal_and_output (out, cal, sgt, drk, xr_ready, tmp_img))
          goto return_status;
     }

   if (0 != out->out_std_metadata (out, meta, ncid_from))
     goto return_status;

   status = 0;
return_status:

   queue_empty (&exprec_queue, gr);
   image_free (tmp_img);
   bpix_free (bpixmap);

   if (ccd) ccd->ccd_delete (ccd);
   if (instr) instr->instr_delete (instr);
   if (pqft) pqft->pqf_delete (pqft);
   if (cal) cal->cal_delete (cal);
   if (drk) drk->drk_close (drk);
   if (sgt) sgt->sgt_delete (sgt);
   if (out)
     {
        (void) out->out_close (out);
        out->out_free (out);
     }

   return status;
}

static int ephem_close (Ephem_Type *eph)
{
   short int error = 0;

   FREE(eph->ephem_name);
   memset ((char *)eph, 0, sizeof (*eph));

   if ((error = novas_ephem_close ()) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: error %d while closing ephemeris",
                     __func__, error);
     }

   return error ? -1 : 0;
}

static char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
     }

   return s_exp;
}

static int ephem_open (config_t *cfg, Ephem_Type *eph)
{
   config_setting_t *s;
   const char *ephem_name;
   short int error;

   if (NULL == (s = config_lookup (cfg, "novas_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing novas_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "ephem_name", &ephem_name))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading ephem_name: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (eph->ephem_name = expand_string (ephem_name)))
     return -1;

   if ((error = novas_ephem_open (eph->ephem_name,
                                  &eph->jd_begin, &eph->jd_end,
                                  &eph->de_number)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: error %d while opening ephemeris %s",
                     __func__, error, eph->ephem_name);
        return -1;
     }

   return 0;
}

static int init_solsys_ephem (config_t *cfg, double tbeg, double tend, Ephem_Type *eph)
{
   double jd_utc0, jd_utc1;

   if (0 != ephem_open (cfg, eph))
     return -1;

   if ((0 != julian_date_from_taix (tbeg, &jd_utc0))
       || (0 != julian_date_from_taix (tend, &jd_utc1)))
     return -1;

   if ((jd_utc0 < eph->jd_begin) || (eph->jd_end < jd_utc1))
     {
        fprintf (stderr, "*** Ephemeris JD=%f-%f doesn't cover JD=%f-%f \n",
                 eph->jd_begin, eph->jd_end,
                 jd_utc0, jd_utc1);
        return -1;
     }

   return 0;
}

int process_inputs (config_t *cfg, const Control_Type *ctrl)
{
   Granule_Type *gr = NULL;
   TIO_Meta_Type *meta = NULL;
   Process_Control_Type pct = {0};
   Ephem_Type eph = {0};
   double tbeg, tend;
   int exposure_type, status = -1;

   if (NULL == (meta = tio_meta_open ()))
     return -1;

   if (0 != tio_meta_set_datetime_production (meta))
     goto return_status;

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
        status = process_dark (cfg, ctrl, &pct, gr, meta);
        break;

      case EXPREC_TYPE_IRR_WRK:
      case EXPREC_TYPE_IRR_REF:
      case EXPREC_TYPE_LIN_IRR:
        /* For irradiances, we'll need to compute the solar illumination geometry */
        tbeg = gr->granule_tstart (gr);
        tend = gr->granule_tend (gr);
        if (0 != init_solsys_ephem (cfg, tbeg, tend, &eph))
          goto return_status;
        /* drop */
      case EXPREC_TYPE_RAD:
        status = process_exposure (cfg, ctrl, &pct, gr, meta);
        break;

      default:
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains exposures of unknown or unexpected type",
                     __func__);
        break;
     }

return_status:

   if (gr) gr->granule_close (gr);
   tio_meta_close (meta);
   (void) ephem_close (&eph);

   return status;
}
