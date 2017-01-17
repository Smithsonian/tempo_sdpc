#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>

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
}
Exprec_Meta_Type;

static int parse_config_file (config_t *cfg, double *bpix_update_thresh,
                              int *bpix_update_num_exprecs_needed)
{
   config_setting_t *setting, *sub;

   if (NULL == (setting = config_lookup (cfg, "ccd_calibration")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_calibration in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((NULL == (sub = config_setting_get_member (setting, "badpix_update")))
        || (CONFIG_TRUE != config_setting_lookup_float (sub, "threshold", bpix_update_thresh))
        || (CONFIG_TRUE != config_setting_lookup_int (sub, "num_exprecs_needed", bpix_update_num_exprecs_needed))
       )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading badpix_update params in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return 0;
}

static int compute_current (CCD_Type *ccd, Pixelqf_Type *pt,
                            Exprec_Meta_Type *xr)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   float *mean_sdc = xr->storage_region_dark;
   Image_Type *aimg = NULL;
   double smear_fraction;
   int i;

   if (-1 == ccd->ccd_correct_coadd (ccd, exprec->num_coadds, exprec->img))
     return -1;

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
                     /(exprec->frame_transfer_time + exprec->exposure_time));

   if (-1 == ccd->ccd_correct_smear (ccd, &smear_fraction, exprec->img))
     return -1;

   if (-1 == ccd->ccd_mean_storage_region_dark (ccd, exprec->img, mean_sdc))
     return -1;

   if (0) fprintf (stderr, "mean sdc:  %7.1f %7.1f %7.1f %7.1f\n",
                   mean_sdc[0], mean_sdc[1], mean_sdc[2], mean_sdc[3]);

   if (NULL == (aimg = ccd->ccd_select_active_pixels (ccd, exprec->img)))
     return -1;
   image_free (exprec->img);
   exprec->img = aimg;

   if (EXPREC_TYPE_IS_DARK(exprec->exposure_type))
     {
        if (-1 == pt->pqf_flag_hotcold (pt, exprec->img))
          return -1;
     }

   if (-1 == pt->pqf_flag_neighbor (pt, exprec->img,
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

        image_scale (xr->img_err, 1.0/exprec->exposure_time);
     }

   image_scale (exprec->img, 1.0/exprec->exposure_time);

   for (i = 0; i < 4; i++)
     {
        xr->storage_region_dark[i] /= exprec->readout_time;
     }

   return 0;
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
             Image_Pixel_Type pxp = prev_pixels[s] * next_pixels[s];
             ref_pixels[s] = (pxp >= 0) ? sqrt(pxp) : IMAGE_PIXEL_FILL_VALUE;
          }
     }
}

static int flag_transients1 (const Pixelqf_Type *pt,
                             Exprec_Meta_Type *exprec_array,
                             int num_exprecs, int exprec_index,
                             Badpix_Map_Type *bpixmap, Image_Type *img_ref)
{
   if (img_ref == NULL)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: img_ref = NULL!", __func__);
        return -1;
     }

   /* While img_ref is assumed to point to allocated storage,
    * we don't always need to write to that storage. When the
    * necessary reference image already exists, it's more efficient
    * to over-write the local copy of the img_ref pointer with the
    * address of the pre-existing reference image.
    */

   if (exprec_index == 0)
     {
        /* over-write local copy of img_ref pointer */
        img_ref = exprec_array[1].exprec->img;
     }
   else if (exprec_index == num_exprecs-1)
     {
        /* over-write local copy of img_ref pointer */
        img_ref = exprec_array[num_exprecs-2].exprec->img;
     }
   else
     {
        /* here's where we actually need the space img_ref points to */
        make_transient_ref_img (exprec_array[exprec_index-1].exprec->img,
                                exprec_array[exprec_index+1].exprec->img,
                                img_ref);
     }

   return pt->pqf_flag_transients (pt, bpixmap->bits, img_ref,
                                   exprec_array[exprec_index].exprec->img);
}

static int flag_transients (const Pixelqf_Type *pt,
                            Exprec_Meta_Type *exprec_array,
                            int num_exprecs, Badpix_Map_Type *bpixmap)
{
   Image_Type *tmp_img = NULL;
   int ixr;

   if (NULL == (tmp_img = image_dup (exprec_array[0].exprec->img)))
     return -1;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "flagging transients");
   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec %3d/%d", ixr, num_exprecs);
        if (-1 == flag_transients1 (pt, exprec_array, num_exprecs, ixr,
                                    bpixmap, tmp_img))
          goto return_error;
     }

   image_free (tmp_img);
   return 0;
return_error:
   image_free (tmp_img);
   return -1;
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
        if (0 != dark_array_elem_init (dark_array, i, exprec->img, sdc, xr->ccd_temp,
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

static int subtract_dark_current_img (Image_Type *img, const Image_Type *dc)
{
   Image_Pixel_Type *img_pixels = img->pixels;
   Image_Pixel_Type *dc_pixels = dc->pixels;
   Image_Pqf_Bitmap_Type *img_pqf = img->pixel_quality_flags;
   Image_Pqf_Bitmap_Type *dc_pqf = dc->pixel_quality_flags;
   int i, n = img->num_rows * img->num_cols;

   for (i = 0; i < n; i++)
     {
        Image_Pixel_Type dc_pixels_i = dc_pixels[i];
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

static int radiometric_correction (config_t *cfg, int exposure_type,
                                   Exprec_Meta_Type *exprec_array, int num_exprecs,
                                   const Control_Type *ctrl)
{
   Calibration_Type *cal = NULL;
   Dark_Table_Type *dtt = NULL;
   Image_Type *tmp_img = NULL;
   int ixr, status = -1;

   if (NULL == (cal = sensorcal_init (cfg)))
     return -1;

   if (NULL == (dtt = dark_table_read (ctrl->dark_file)))
     goto return_status;

   if (0 != validate_dark_table_type (dtt, exposure_type))
     goto return_status;

   tmp_img = exprec_array[0].exprec->img;
   if (NULL == (tmp_img = image_new (tmp_img->num_rows, tmp_img->num_cols)))
     goto return_status;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "Radiometric correction:");

   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        Exprec_Meta_Type *xr = &exprec_array[ixr];
        Granule_Exprec_Type *exprec = xr->exprec;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec %3d/%d", ixr, num_exprecs);

        /* FIXME? propagate dark-current subtraction error to
         * uncertainty estimate? */
        if (0 != subtract_dark_current (exprec, dtt, xr->storage_region_dark,
                                        xr->ccd_temp, tmp_img))
          {
             goto return_status;
          }

        /* FIXME: update uncertainty? */
        if (0 != cal->cal_apply_prnu (cal, exprec->img))
          goto return_status;

        /* >>> Stray light correction goes here <<< */

        /* Multiplicative factor converts e/s to photons/s */
        if ((0 != cal->cal_apply_rcoeffs (cal, exprec->img))
            || (0 != cal->cal_apply_rcoeffs (cal, xr->img_err)))
          goto return_status;

        if (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type))
          {
             double solar_phi=0.0, solar_theta=0.0;  /* FIXME: placeholders */
             if ((0 != cal->cal_apply_btdf (cal, solar_phi, solar_theta, exprec->img))
                 || (0 != cal->cal_apply_btdf (cal, solar_phi, solar_theta, xr->img_err))) /* FIXME: ok? */
               goto return_status;
          }

        /* Wavelength calibration could go here or, to reduce
         * RAM usage, we could write out the current data, free
         * the exprec array, and then do wavelength calibration
         * one frame at a time, re-reading data as needed from the
         * "output" file.
         * FIXME: are we going to store explicit wavelength arrays,
         * or are we going to implement something clever?
         */
     }

   status = 0;
return_status:
   if (cal) cal->cal_delete (cal);
   if (dtt) dtt->dtt_delete (dtt);
   image_free (tmp_img);
   return status;
}

static int validate_exposure_type (int *exposure_type0, int exposure_type)
{
   if (exposure_type == EXPREC_TYPE_UNKNOWN)
     {
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains exposures of unknown type",
                     __func__);
        return -1;
     }

   if (*exposure_type0 == EXPREC_TYPE_UNKNOWN)
     {
        *exposure_type0 = exposure_type;
        return 0;
     }

   if (exposure_type != *exposure_type0)
     {
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains multiple exposure types",
                     __func__);
        return -1;
     }

   return 0;
}

static int write_radiometric (config_t *cfg,
                              const Exprec_Meta_Type *exprec_array,
                              int num_exprecs, const char *output_file)
{
   Output_Type *out = NULL;
   Granule_Exprec_Type *exprec;
   Image_Type *img;
   int k, status = -1;

   if (num_exprecs == 0)
     return 0;

   exprec = exprec_array[0].exprec;
   if (NULL == (out = output_alloc (cfg, exprec->exposure_type)))
     return -1;

   if (0 != out->out_set_file (out, output_file))
     goto return_status;

   img = exprec->img;
   out->out_set_dims (out, num_exprecs, img->num_cols, img->num_rows/2);

   if (0 != out->out_create (out))
     goto return_status;

   for (k = 0; k < num_exprecs; k++)
     {
        Output_Exprec_Type rec;
        const Exprec_Meta_Type *xr = &exprec_array[k];
        rec.exprec = xr->exprec;
        rec.img_err = xr->img_err;
        if (0 != out->out_write_rec (out, k, &rec))
          goto return_status;
     }

   status = 0;
return_status:
   if (out)
     {
        (void) out->out_close (out);
        out->out_free (out);
     }

   return status;
}

int process_inputs (config_t *cfg, const Control_Type *ctrl)
{
   CCD_Type *ccd = NULL;
   Granule_Type *gr = NULL;
   Granule_Exprec_Type *exprec = NULL;
   Instr_Type *instr = NULL;
   Exprec_Meta_Type *exprec_array = NULL;
   Pixelqf_Type *pt = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Badpix_Map_Occur_Type *bpix_occur = NULL;
   Badpix_Bitmap_Type bpix_occur_mask = IMAGE_PQF_HOT_PIXEL | IMAGE_PQF_COLD_PIXEL;
   double bpix_update_thresh;
   int ixr, num_exprecs, exposure_type = EXPREC_TYPE_UNKNOWN;
   int bpix_occur_threshold, bpix_update_num_exprecs_needed;
   int status = -1;

   if (NULL == (ccd = ccd_init (cfg)))
     goto return_status;

   if (NULL == (pt = pqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file)))
     goto return_status;

   if (NULL == (gr = granule_open (ctrl->input_file)))
     goto return_status;

   num_exprecs = gr->granule_num_exprecs(gr);
   if (ctrl->limit_num_granules < num_exprecs)
     num_exprecs = ctrl->limit_num_granules;

   if (NULL == (exprec_array = alloc_exprec_array (num_exprecs)))
     goto return_status;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;

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

        if (-1 == validate_exposure_type (&exposure_type, exprec->exposure_type))
          goto return_status;

        /* FIXME: correct temp? */
        if (0 != instr->instr_ccd_temp1 (instr, exprec->start_time, &xr->ccd_temp))
          {
             tell_vlog (TELL_MSGTYPE_WARN, 0,
                        "%s: ccd_temp lookup failed, start_time=%15.12e",
                        __func__, exprec->start_time);
             /* drop */
          }

        if (-1 == compute_current (ccd, pt, xr))
          goto return_status;

        if (-1 == bpix_occur_incr (bpix_occur,
                                   exprec->img->pixel_quality_flags))
          goto return_status;
     }

   if (-1 == parse_config_file (cfg, &bpix_update_thresh, &bpix_update_num_exprecs_needed))
     goto return_status;

   bpix_occur_threshold = bpix_update_thresh * num_exprecs;
   if (num_exprecs > bpix_update_num_exprecs_needed)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "updating bpix map");
        if (0 != bpix_occur_set (bpix_occur, bpix_occur_threshold,
                                 bpix_occur_mask, bpixmap->bits))
          goto return_status;
     }

   /* FIXME - always do this? */
   if (0 != flag_transients (pt, exprec_array, num_exprecs, bpixmap))
     goto return_status;

   switch (exposure_type)
     {
      default:
        /* Given the validation above, this should never happen */
        tell_verror (TELL_RUNTIME_ERROR, "%s: unknown exposure type: %d",
                     __func__, exposure_type);
        goto return_status;

      case EXPREC_TYPE_DARK:
      case EXPREC_TYPE_LIN_DARK: {
         int is_linearity = (exposure_type == EXPREC_TYPE_LIN_DARK);
         if (0 != make_dark_table (cfg, is_linearity, exprec_array, num_exprecs, ctrl->output_file))
           goto return_status;
      }
        break;

      case EXPREC_TYPE_RADIANCE:
      case EXPREC_TYPE_IRRADIANCE:
      case EXPREC_TYPE_LIN_IRR:
        if (0 != radiometric_correction (cfg, exposure_type, exprec_array, num_exprecs, ctrl))
          goto return_status;
        if (0 != write_radiometric (cfg, exprec_array, num_exprecs, ctrl->output_file))
          goto return_status;
        break;
     }

   status = 0;
return_status:

   bpix_free (bpixmap);
   bpix_occur_close (bpix_occur);
   free_exprec_array (exprec_array, num_exprecs, gr);

   if (ccd)
     {
        ccd->ccd_delete (ccd);
        ccd = NULL;
     }

   if (gr)
     {
        gr->granule_close (gr);
        gr = NULL;
     }

   if (instr)
     {
        instr->instr_delete (instr);
        instr = NULL;
     }

   if (pt)
     {
        pt->pqf_delete (pt);
        pt = NULL;
     }

   return status;
}
