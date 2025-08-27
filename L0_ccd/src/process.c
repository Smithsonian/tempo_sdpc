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
#include <iocsdpc.h>
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
#include "_process.h"
#include "current.h"
#include "util.h"

typedef struct
{
   int saturated_neighbor_hw_serial;
   int saturated_neighbor_hw_parallel;
   double bpix_update_thresh;
   int bpix_update_num_exprecs_needed;
   double exposure_time_per_dark_int_frame;
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
static int _pProcessing_Version = 1;

void process_set_version (int version)
{
   _pProcessing_Version = version;
}

int process_get_version (void)
{
   return _pProcessing_Version;
}

/*{{{ Diagnostic image output */

static struct
{
   int index;
   int ncid;
}
Diagnostic_Controls;

static int set_diagnostic_controls_index (int index)
{
   Diagnostic_Controls.index = index;
   return 0;
}

/* Setting environment variable SDPC_DIAGNOSTIC_INDEX to
 * an integer >= 0 turns on diagnostic output for the specified
 * image index.
 */
static int want_diagnostic_output (int index)
{
   return (index == Diagnostic_Controls.index);
}

static int close_diagnostic_file (void)
{
   if (Diagnostic_Controls.index < 0)
     return 0;
   return TIO_close (Diagnostic_Controls.ncid);
}

static int create_diagnostic_file (const Control_Type *ctrl, int num_parallel_active_full, int num_serial_active_full,
                                   int enabled_straylight_correction, int is_irradiance)
{
   const char *suffix = "_diag.nc";
   const char *output_file = ctrl->output_file;
   int img_dimids[2];
   char *diag_file = NULL;
   char *dot;
   size_t len;
   int varid, ncid;
   int status = -1;

   if (ctrl->diagnostic_index < 0)
     return 0;

   len = strlen(output_file) + strlen(suffix);
   if (NULL == (diag_file = (char *)MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   strncpy (diag_file, output_file, len);
   if (NULL != (dot = strrchr (diag_file, '.')))
     {
        sprintf (dot, "%s", suffix);
     }
   else strncat (diag_file, suffix, len);

   if (0 != TIO_create (diag_file, NC_NETCDF4, &ncid))
     goto return_status;

   if ((0 != TIO_def_dim (ncid, "row", num_parallel_active_full, &img_dimids[0]))
       || (0 != TIO_def_dim (ncid, "col", num_serial_active_full, &img_dimids[1])))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining file array dimensions", __func__);
        goto return_status;
     }

   if ((0 != TIO_def_var (ncid, "dark", TIO_FLOAT, 2, img_dimids, &varid))
       ||(0 != TIO_def_var (ncid, "img_before_dark_subtract", TIO_FLOAT, 2, img_dimids, &varid))
       ||(0 != TIO_def_var (ncid, "img_after_dark_subtract", TIO_FLOAT, 2, img_dimids, &varid))
      )
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark variables", __func__);
        goto return_status;
     }

   if ((0 != TIO_def_var (ncid, "img_before_radiometric_correction", TIO_FLOAT, 2, img_dimids, &varid))
       ||(0 != TIO_def_var (ncid, "img_after_radiometric_correction", TIO_FLOAT, 2, img_dimids, &varid)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining radiometric variables", __func__);
        goto return_status;
     }

   if (enabled_straylight_correction)
     {
        if ((0 != TIO_def_var (ncid, "img_before_straylight_correction", TIO_FLOAT, 2, img_dimids, &varid))
            ||(0 != TIO_def_var (ncid, "img_after_straylight_correction", TIO_FLOAT, 2, img_dimids, &varid)))
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: defining straylight variables", __func__);
             goto return_status;
          }
     }

   if (is_irradiance)
     {
        if ((0 != TIO_def_var (ncid, "img_before_btdf_correction", TIO_FLOAT, 2, img_dimids, &varid))
            ||(0 != TIO_def_var (ncid, "btdf", TIO_FLOAT, 2, img_dimids, &varid))
            ||(0 != TIO_def_var (ncid, "diffuser_polcorr", TIO_FLOAT, 2, img_dimids, &varid))
            ||(0 != TIO_def_var (ncid, "img_after_btdf_correction", TIO_FLOAT, 2, img_dimids, &varid)))
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: defining irradiance variables", __func__);
             goto return_status;
          }
     }

   Diagnostic_Controls.ncid = ncid;

   status = 0;
return_status:
   FREE(diag_file);

   return status;
}

static int write_diagnostic_image (const Image_Type *img, const char *varname)
{
   int ncid = Diagnostic_Controls.ncid;
   int start[2], count[2];

   if (img == NULL)
     return 0;

   start[0] = 0;
   start[1] = 0;
   count[0] = img->num_rows;
   count[1] = img->num_cols;

   if (0 != TIO_put_var_section (ncid, varname, start, count, TIO_FLOAT, img->pixels))
     return -1;

   return 0;
}

/*}}}*/

static int write_std_metadata (int ncid, int exposure_type, const Control_Type *ctrl,
                               TIO_Meta_Type *meta, int ncid_from)
{
#define SHORTNAME_BUFSIZE 32
   char shortname[SHORTNAME_BUFSIZE];
   const char *prod_name = NULL;
   const char *template_basename = NULL;
   char *template_path = NULL;
   int n, status = -1;

   switch (exposure_type)
     {
      case EXPREC_TYPE_DARK:
        prod_name = TEMPO_PROD_TYPE_DRK;
        template_basename = "dark.met.template";
        break;

      case EXPREC_TYPE_IRR_WRK:
        prod_name = TEMPO_PROD_TYPE_IRR;
        template_basename = "irradiance.met.template";
        break;

      case EXPREC_TYPE_IRR_REF:
        prod_name = TEMPO_PROD_TYPE_IRR_REF;
        template_basename = "irradiance.met.template";
        break;

      case EXPREC_TYPE_RAD:
        prod_name = TEMPO_PROD_TYPE_RAD;
        template_basename = "radiance.met.template";
        break;

      case EXPREC_TYPE_RAD_TWI:
        prod_name = TEMPO_PROD_TYPE_RAD_TWI;
        template_basename = "radiance.met.template";
        break;

      default:
        tell_vwarn (0, "%s: no metadata template expansion support for exposure records of type %d",
                    __func__, exposure_type);
        break;
     }

   n = snprintf (shortname, SHORTNAME_BUFSIZE, "TEMPO_%s_L1", prod_name);
   if (n < 0 || n >= SHORTNAME_BUFSIZE)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error generating shortname for %s", __func__, prod_name);
        goto return_status;
     }

   if (0 != tio_meta_set_standard (meta, ctrl->output_file, shortname, process_get_version(), ctrl->pge_version_string))
     goto return_status;

   if (0 != tio_meta_set_datetime_range (meta, ncid_from))
     goto return_status;

   if (0 != tio_meta_write_ncattr (meta, ncid))
     goto return_status;

   if ((exposure_type == EXPREC_TYPE_RAD)
       || (exposure_type == EXPREC_TYPE_RAD_TWI)
       || (exposure_type == EXPREC_TYPE_IRR_WRK)
       || (exposure_type == EXPREC_TYPE_IRR_REF))
     {
        /* For radiance files, input_pointer gets expanded only in the
         * last processing step of Level 0-1, e.g. post-INR
         * For irradiance files, input_pointer gets expanded after
         * wavelength calibration.
         */
        tio_meta_set_noexpand (meta, "input_files", 1);
     }

   if ((ctrl->metadata_template_dir != NULL)
       && (template_basename != NULL))
     {
        if (NULL == (template_path = path_concat (ctrl->metadata_template_dir, template_basename)))
          goto return_status;
        tell_vlog (TELL_MSGTYPE_INFO, 1, "Expanding metadata template: %s", template_path);
        if (0 != tio_meta_expand_file (meta, template_path, ctrl->output_file))
          goto return_status;
     }

   status = 0;
return_status:
   FREE(template_path);
   return status;
}

static int get_control_params (config_t *cfg, Process_Control_Type *pct)
{
   config_setting_t *s, *sub;

   /* This parameter is optional */
   (void) config_lookup_float (cfg, "calibration.exposure_time_per_dark_int_frame",
                               &pct->exposure_time_per_dark_int_frame);

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
   if ((xr == NULL) || (gr == NULL))
     return;
   gr->granule_free_exprec (xr->exprec);
   image_free (xr->img_err);
   trend_collect_free_record (xr->tr);
   memset ((char *)xr, 0, sizeof (*xr));
   FREE(xr);
}

static Exprec_Meta_Type *alloc_exprec_meta (Trend_File_Type *tft)
{
   Exprec_Meta_Type *xr = NULL;

   if (NULL == (xr = (Exprec_Meta_Type *) MALLOC (sizeof *xr)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)xr, 0, sizeof *xr);

   if (NULL == (xr->tr = trend_collect_new_record (tft)))
     {
        FREE(xr);
        return NULL;
     }

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
   float fpa_temp, fpe_temp, spec_temp, tele_temp, bench_temp, fpa_sum, fpe_sum, spec_sum, tele_sum, bench_sum;
   double smear_fraction, redmine204_shift;
   double exposure_time_per_frame, exposure_time_offset, coadd_period, integration_period, sampling_period;
   int i, k, n_sample, n_fpa, n_fpe, n_spec, n_tele, n_bench;
   unsigned int j;

   if (-1 == ccd->ccd_correct_coadd (ccd, exprec->num_coadds, exprec->img))
     return -1;
   exposure_time_per_frame = exprec->exposure_time / exprec->num_coadds;

   if ((exposure_time_per_frame <= 0.0) || (0 == isfinite (exposure_time_per_frame)))
     {
        tell_vwarn (0, "%s: exposure time per frame = %f  (ccd_int_type=%d)", __func__,
                    exposure_time_per_frame, exprec->ccd_int_type);
        if (exprec->ccd_int_type == 3)
          {
             /* DARK_INT mode:
              * Because of a FSW bug, we always have FPE_DARK_LINES=0 in telemetry,
              * so DARK_INT exposure records always have exposure_time=0, which is wrong.
              * As a workaround, we assign a default exposure_time value to enable
              * trending.  The default value can be changed by setting the parameter
              *     calibration.exposure_time_per_dark_int_frame
              * in the config file.
              */
             double etpdif = pct->exposure_time_per_dark_int_frame;
             exposure_time_per_frame = (etpdif > 0.0) ? etpdif : 0.071988;  /* sec */
             tell_vwarn (0, "%s: assuming exposure time per frame = %f  (ccd_int_type=%d)", __func__,
                         exposure_time_per_frame, exprec->ccd_int_type);
          }
     }

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "coadd");

   if (0 != ccd->ccd_configure_using_octant_phase (ccd, exprec->img))
     return -1;

   if (EXPREC_TYPE_IS_DARK(exprec->exposure_type))
     {
        if (0 != ccd->ccd_compute_readnoise (ccd, exprec->num_coadds, exprec->img))
          return -1;
     }

   if (0 != ccd->ccd_correct_offset (ccd, exprec->img))
     return -1;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "offset");

   /* It's not expected that data from a linearity sweep will be processed
    * to Level 1, but if it is, then we'll treat it exactly the same
    * as any other dark or irradiance measurement.
    */

   if (-1 == ccd->ccd_correct_nonlinearity (ccd, exprec->img))
     return -1;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "nonlinearity");

   if (0 != ccd->ccd_correct_crosstalk (ccd, exprec->img))
     return -1;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "crosstalk");

   fpa_sum   = 0.0;
   fpe_sum   = 0.0;
   spec_sum  = 0.0;
   tele_sum  = 0.0;
   bench_sum = 0.0;

   n_fpa   = 0;
   n_fpe   = 0;
   n_spec  = 0;
   n_tele  = 0;
   n_bench = 0;

   /* FIXME: The timing patterns below have been implemented on an empirical basis.
    * The nominal integration mode seems to work reasonably, but there are rooms
    * for improvement for the other modes.
    * It should be revisited once we gain a precise understanding of
    * how it works for different integration modes.
    *
    * September 25, 2024: Code revised to un-tweak image time.
    * This is a temporary solution. New dark correction will be implemented.
    */

   double secs_per_clock      = 8.2987551867219914e-08;
   double frame_transfer_time =  100416 * secs_per_clock;
   double full_flush_time     =  199104 * secs_per_clock;
   double storage_read_time   = 1204992 * secs_per_clock;

   redmine204_shift = 0.0;
   (void) iocsdpc_tweak_image_time_per_redmine_204 (exprec->ccd_int_type, exposure_time_per_frame, &redmine204_shift);

   switch (exprec->ccd_int_type)
     {
      case 3:
      /* DARK_INT */
      case 1:
      /* SHORT_INT */
        integration_period = exposure_time_per_frame;
        coadd_period = full_flush_time + integration_period + frame_transfer_time + storage_read_time;
        exposure_time_offset = (exprec->start_time - redmine204_shift) - coadd_period;
        n_sample = 1;
        break;
      case 0:
      /* NOMINAL */
        integration_period = exposure_time_per_frame;
        coadd_period = integration_period + frame_transfer_time;
        exposure_time_offset = (exprec->start_time - redmine204_shift);
        n_sample = 1;
        break;
      case 2:
      /* LONG_INT */
        integration_period = exposure_time_per_frame - storage_read_time;
        coadd_period = frame_transfer_time + storage_read_time + integration_period;
        exposure_time_offset = (exprec->start_time - redmine204_shift);
        n_sample = 30;
        break;
      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: invalid ccd_int_type=%d", __func__, exprec->ccd_int_type);
        return -1;
     }
   if (exprec->num_coadds == 1)
     {
        exposure_time_offset -= exposure_time_per_frame;
     }
   sampling_period = integration_period / (n_sample + 1);

   for (j = 0; j < exprec->num_coadds; j++)
     for (k = 1; k < (n_sample+1); k++)
       {
         if (0 != instr->instr_fpa_temp (instr, exposure_time_offset + (coadd_period * j) + (sampling_period * k), &fpa_temp))
           {
               return -1;
           }
         else
           {
               fpa_sum += fpa_temp;
               n_fpa += 1;
           }
         if (0 != instr->instr_fpe_temp (instr, exposure_time_offset + (coadd_period * j) + (sampling_period * k), &fpe_temp))
           {
               return -1;
           }
         else
           {
               fpe_sum += fpe_temp;
               n_fpe += 1;
           }
         if (0 != instr->instr_spec_temp (instr, exposure_time_offset + (coadd_period * j) + (sampling_period * k), &spec_temp))
           {
               return -1;
           }
         else
           {
               spec_sum += spec_temp;
               n_spec += 1;
           }
         if (0 != instr->instr_tele_temp (instr, exposure_time_offset + (coadd_period * j) + (sampling_period * k), &tele_temp))
           {
               return -1;
           }
         else
           {
               tele_sum += tele_temp;
               n_tele += 1;
           }
         if (0 != instr->instr_bench_temp (instr, exposure_time_offset + (coadd_period * j) + (sampling_period * k), &bench_temp))
           {
               return -1;
           }
         else
           {
               bench_sum += bench_temp;
               n_bench += 1;
           }
       }
   xr->fpa_temp = fpa_sum / n_fpa;
   tell_vlog (TELL_MSGTYPE_INFO, 2, "FPA temp lookup: time: %0.3f => FPA temp: %0.2f C", exprec->start_time, xr->fpa_temp);

   xr->fpe_temp = fpe_sum / n_fpe;
   tell_vlog (TELL_MSGTYPE_INFO, 2, "FPE temp lookup: time: %0.3f => FPE temp: %.2f C", exprec->start_time, xr->fpe_temp);

   spec_temp = spec_sum / n_spec;
   tell_vlog (TELL_MSGTYPE_INFO, 2, "Spectrometer temp lookup: time %0.3f => Spectrometer temp: %0.2f C", exprec->start_time, spec_temp);

   tele_temp = tele_sum / n_tele;
   tell_vlog (TELL_MSGTYPE_INFO, 2, "Telescope temp lookup: time %0.3f => Telescope temp: %0.2f C", exprec->start_time, tele_temp);

   bench_temp = bench_sum / n_bench;
   tell_vlog (TELL_MSGTYPE_INFO, 2, "Bench temp lookup: time %0.3f => Bench temp: %0.2f C", exprec->start_time, bench_temp);

   if (0 != trend_collect_temp (spec_temp, tele_temp, bench_temp))
     return -1;

   /* convert DN to electrons */
   if (-1 == ccd->ccd_correct_gain (ccd, exprec->img, xr->fpa_temp, xr->fpe_temp))
     return -1;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "gain");

   if (-1 == ccd->ccd_mean_storage_region_dark (ccd, exprec->img,
                                                exprec->num_dg_rows, exprec->num_tg_rows,
                                                xr->storage_region_dark))
     {
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "mean sdc [offset=%d, num=%d]:  %f %f %f %f",
              exprec->num_dg_rows, exprec->num_tg_rows,
              xr->storage_region_dark[0],
              xr->storage_region_dark[1],
              xr->storage_region_dark[2],
              xr->storage_region_dark[3]);

   /* Compute noisesq before smear correction */
   if ((exprec->exposure_type == EXPREC_TYPE_RAD) || (exprec->exposure_type == EXPREC_TYPE_RAD_TWI)
       || (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type)))
     {
        if (0 != compute_noisesq_for_active_pixels (ccd, exprec->img, xr))
          return -1;
     }

   smear_fraction = (exprec->frame_transfer_time
                     /(exprec->frame_transfer_time + exposure_time_per_frame));

   if (-1 == ccd->ccd_correct_smear (ccd, &smear_fraction, exprec->img))
     return -1;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "smear");

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
   if (exposure_time_per_frame > 0.0)
     {
        image_scale (exprec->img, 1.0/exposure_time_per_frame);
        if (xr->img_err)
          {
             image_scale (xr->img_err, 1.0/exprec->num_coadds);
             image_sqrt (xr->img_err);
             image_scale (xr->img_err, 1.0/exposure_time_per_frame);
          }
     }
   for (i = 0; i < 4; i++)
     {
        xr->storage_region_dark[i] /= exprec->readout_time;
     }

   if (0 != trend_collect_sdc (exprec->num_dg_rows, exprec->num_tg_rows, xr->storage_region_dark))
     return -1;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "integration_time");

   if (0 != ccd->ccd_correct_prnu (ccd, exprec->img))
     return -1;

   if (EXPREC_TYPE_IS_DARK(exprec->exposure_type))
     {
        if (-1 == pqft->pqf_flag_hotcold (pqft, exprec->img, &xr->dark_trend))
          return -1;
     }

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (exprec->img, "prnu");

   return 0;
}

static int derive_current (config_t *cfg, const Control_Type *ctrl, Process_Control_Type *pct,
                           Granule_Type *gr, TIO_Meta_Type *meta, Trend_File_Type *tft)
{
   CCD_Type *ccd = NULL;
   Instr_Type *instr = NULL;
   Exprec_Meta_Type *xr = NULL;
   Pixelqf_Type *pqft = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Badpix_Map_Occur_Type *bpix_occur = NULL;
   Badpix_Bitmap_Type bpix_occur_mask;
   int ixr, num_exprecs, exposure_type, is_dark, grp, ncid_from;
   int num_parallel_active_full, num_serial_active_full, ncid = 0;
   int bpix_occur_threshold;
   int status = -1;

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   is_dark = EXPREC_TYPE_IS_DARK(exposure_type);

   num_exprecs = gr->granule_num_exprecs(gr);
   if (ctrl->limit_num_granules < num_exprecs)
     num_exprecs = ctrl->limit_num_granules;

   if (NULL == (ccd = ccd_init (cfg, meta)))
     goto return_status;

   if (NULL == (pqft = pixelqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file, ctrl->instr_glob,
                                    gr->granule_tstart(gr), gr->granule_tend(gr), meta)))
     goto return_status;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;
   if (0)
     {
        if (0 != meta_record_basename (meta, ctrl->bpix_file))
          goto return_status;
     }

   if (is_dark)
     {
        bpix_occur_mask = IMAGE_PQF_HOT_PIXEL | IMAGE_PQF_COLD_PIXEL;
        bpix_occur = bpix_occur_open (bpixmap->num_rows, bpixmap->num_cols,
                                      bpix_occur_mask);
        if (NULL == bpix_occur)
          goto return_status;
     }

   /* Open the output file */
   ccd->ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);
   if (0 != current_create_file_of_type (gr, ctrl->output_file, num_exprecs,
                                         num_parallel_active_full, num_serial_active_full,
                                         &ncid, &grp))
     {
        goto return_status;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "Converting DN to e-/s:");
   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "exposure record %3d/%d", (ixr+1), num_exprecs);

        if (NULL == (xr = alloc_exprec_meta (tft)))
          goto return_status;

        xr->index = ixr;

        if (NULL == (xr->exprec = gr->granule_read_exprec_by_index (gr, ixr, NULL)))
          goto return_status;

        if (xr->exprec->num_coadds == TIO_FILL_UINT)
          {
             tell_vwarn (0, "%s: exposure record %d/%d (missing data)", __func__, (ixr+1), num_exprecs);
             free_exprec_meta (xr, gr);
             xr = NULL;
             continue;
          }

        trend_collect_set_active_record (xr->tr);

        if (0 != trend_collect_time (xr->exprec->start_time, ixr))
          goto return_status;

        if (-1 == validate_exposure_type (exposure_type, xr->exprec->exposure_type))
          goto return_status;
        if (0 != ccd->ccd_apply_pixel_quality_flags (ccd, xr->exprec->img,
                                                     bpixmap->bits, bpixmap->num_rows, bpixmap->num_cols))
          goto return_status;

        if (-1 == compute_current_and_trim (ccd, instr, pqft, pct, xr))
          goto return_status;

        if (0 != current_write_exprec (grp, xr))
          goto return_status;

        if (0 != trend_collect_pqf (xr->exprec->img))
          goto return_status;
        if (0 != trend_collect_write_record (xr->tr))
          goto return_status;

        if (is_dark)
          {
             if (-1 == bpix_occur_incr (bpix_occur, xr->exprec->img->pixel_quality_flags))
               goto return_status;
          }

        free_exprec_meta (xr, gr);
        xr = NULL;
     }

   if (is_dark)
     {
        bpix_occur_threshold = pct->bpix_update_thresh * num_exprecs;
        if (num_exprecs > pct->bpix_update_num_exprecs_needed)
          {
             tell_vlog (TELL_MSGTYPE_INFO, 1, "updating internal bpix map");
             if (0 != bpix_occur_set (bpix_occur, bpix_occur_threshold,
                                      bpix_occur_mask, bpixmap->bits))
               goto return_status;
          }
     }

   ncid_from = gr->granule_ncid (gr);

   if (0 != write_std_metadata (ncid, exposure_type, ctrl, meta, ncid_from))
     goto return_status;

   /* Generate average only for "true" dark measurements,
    * but NOT for linearity data */
   if (exposure_type == EXPREC_TYPE_DARK)
     {
        if (0 != current_write_mean_dark_current (ncid))
          goto return_status;
     }

   /* FIXME: eventually, we may write out an updated badpix map */

   status = 0;
return_status:

   bpix_free (bpixmap);
   bpix_occur_close (bpix_occur);
   free_exprec_meta (xr, gr);
   xr = NULL;

   if (ncid)
     {
        if (0 != TIO_close (ncid))
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
        if (dc_pixels_i > img_pixels[i])
          {
             img_pqf[i] |= IMAGE_PQF_DARK_CORR_ERROR;
          }
        if ((dc_pqf[i] == 0) && (0 <= dc_pixels_i))
          {
             img_pixels[i] -= dc_pixels_i;
          }
        else img_pqf[i] |= IMAGE_PQF_DARK_CORR_ERROR;
      }

   return 0;
}

static int dark_subtract (const Dark_Type *drk, Exprec_Meta_Type *xr, Image_Type *tmp_img)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   const char *method = enable_state_query_enum (ENABLE_DARK);

   if (0 == strcmp (method, "none"))
     return 0;

   /* copy the appropriate dark current image into tmp_img */
   if (0 != drk->drk_image (drk, tmp_img))
     return -1;

   if (0 == strcmp (method, "mean_tfpa_adj"))
     {
        if (0 != drk->drk_image_Tfpa_adj (drk, xr->fpa_temp, tmp_img))
          return -1;
     }
   else if (0 == strcmp (method, "mean_sdc_adj"))
     {
        if (0 != drk->drk_image_sdc_adj (drk, xr->storage_region_dark, tmp_img))
          return -1;
     }

   if (want_diagnostic_output (xr->index))
     {
        (void) write_diagnostic_image (tmp_img, "dark");
        (void) write_diagnostic_image (exprec->img, "img_before_dark_subtract");
     }

   /* subtract the dark current image, leaving the result in place */
   if (0 != subtract_dark_current_img (exprec->img, tmp_img))
     return -1;

   if (want_diagnostic_output (xr->index))
     {
        (void) write_diagnostic_image (exprec->img, "img_after_dark_subtract");
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
                                   Exprec_Meta_Type *xr)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   double jd_utc;

   if (cal->cal_straylight_correction)
     {
        if (want_diagnostic_output(xr->index))
          (void) write_diagnostic_image (exprec->img, "img_before_straylight_correction");

        if (0 != cal->cal_straylight_correction (cal, exprec->img))
          return -1;

        if (want_diagnostic_output(xr->index))
          (void) write_diagnostic_image (exprec->img, "img_after_straylight_correction");
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "radiometric correction");

   if (want_diagnostic_output(xr->index))
     (void) write_diagnostic_image (exprec->img, "img_before_radiometric_correction");

   /* Multiplicative factor converts e/s to photons/s */
   if ((0 != cal->cal_apply_radcal_coeffs (cal, exprec->img))
       || (0 != cal->cal_apply_radcal_coeffs (cal, xr->img_err)))
     return -1;

   if (want_diagnostic_output(xr->index))
     (void) write_diagnostic_image (exprec->img, "img_after_radiometric_correction");

   /* For the irradiance, we need the angular position of the sun
    * relative to the boresight, and the earth-sun distance.
    * For the radiance, we eventually need only the earth-sun distance,
    * and that gets done elsewhere anyway, but to keep the code simple,
    * we compute the solar position for both.  It's cheap, so no worries.
    */
   if ((0 != julian_date_from_taix (exprec->start_time, &jd_utc))
       || (0 != sgt->sgt_sat_sun_position (sgt, jd_utc, &xr->solar_theta, &xr->solar_phi, &xr->earth_sun_distance)))
     return -1;

   if (EXPREC_TYPE_IS_IRRADIANCE(exprec->exposure_type))
     {
        int use_reference_diffuser = (exprec->exposure_type == EXPREC_TYPE_IRR_REF);
        Image_Type *img_btdf = NULL;
        Image_Type *img_polcorr = NULL;

        tell_vlog (TELL_MSGTYPE_INFO, 1, "BTDF correction");

        if (0 != trend_collect_solar_angles (xr->solar_theta, xr->solar_phi, use_reference_diffuser))
          return -1;

        if (want_diagnostic_output(xr->index))
          {
             (void) write_diagnostic_image (exprec->img, "img_before_btdf_correction");
             /* try to allocate space for diagnostic images, but ignore NULLs */
             img_btdf = image_new (exprec->img->num_rows, exprec->img->num_cols);
             img_polcorr = image_new (exprec->img->num_rows, exprec->img->num_cols);
          }

        if ((0 != cal->cal_apply_btdf (cal, use_reference_diffuser, xr->solar_phi, xr->solar_theta, exprec->img, img_btdf))
            || (0 != cal->cal_apply_btdf (cal, use_reference_diffuser, xr->solar_phi, xr->solar_theta, xr->img_err, NULL)))
          {
             image_free (img_btdf);
             image_free (img_polcorr);
             return -1;
          }

        if ((0 != cal->cal_apply_diffuser_polcorr (cal, xr->solar_phi, xr->solar_theta, exprec->img, img_polcorr))
            || (0 != cal->cal_apply_diffuser_polcorr (cal, xr->solar_phi, xr->solar_theta, xr->img_err, NULL)))
          {
             image_free (img_btdf);
             image_free (img_polcorr);
             return -1;
          }

        if (want_diagnostic_output(xr->index))
          {
             (void) write_diagnostic_image (img_btdf, "btdf");
             (void) write_diagnostic_image (img_polcorr, "diffuser_polcorr");
             (void) write_diagnostic_image (exprec->img, "img_after_btdf_correction");
             image_free (img_btdf);
             image_free (img_polcorr);
             img_btdf = NULL;
             img_polcorr = NULL;
          }
     }

   return 0;
}

static Spectral_Data_Type *finalize_band (const Calibration_Type *cal,
                                          const Exprec_Meta_Type *xr, int band_id)
{
   Image_Type *img = xr->exprec->img;
   Image_Type *img_err = xr->img_err;
   Spectral_Data_Type *sdt = NULL;
   int status;

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

static int radcal_and_output (Output_Type *out, Calibration_Type *cal, Solar_Geom_Type *sgt,
                              Exprec_Meta_Type *xr)
{
   Output_Exprec_Type outrec = {0};
   Trend_Record_Type *tr_old = NULL;
   int num_negative, status = -1;

   tr_old = trend_collect_set_active_record (xr->tr);

   if (0 != radiometric_correction (cal, sgt, xr))
     goto return_status;

   if (want_diagnostic_output(xr->index)) (void) image_write_raw (xr->exprec->img, "final");

   if ((num_negative = image_check_negative_pixels (xr->exprec->img, 1)) < 0)
     goto return_status;
   if (num_negative > 0)
     {
	tell_vwarn (0, "%s: set processing error bit in %d pixels with ((value<0) && (pqf==0))",
		    __func__, num_negative);
     }

   if ((NULL == (outrec.uv = finalize_band (cal, xr, TEMPO_BAND_UV)))
       || (NULL == (outrec.vis = finalize_band (cal, xr, TEMPO_BAND_VIS))))
     goto return_status;

   if ((0 != trend_collect_pqf_uv (outrec.uv->pqf, outrec.uv->num_xtrack, outrec.uv->num_channels))
       || (0 != trend_collect_pqf_vis (outrec.vis->pqf, outrec.vis->num_xtrack, outrec.vis->num_channels)))
     goto return_status;

   outrec.meta.start_time = xr->exprec->start_time;
   outrec.meta.exposure_time = xr->exprec->exposure_time;
   outrec.meta.mirror_step = xr->exprec->curr_mirror_step;
   outrec.meta.earth_sun_distance = xr->earth_sun_distance;
   outrec.meta.solar_phi = xr->solar_phi;
   outrec.meta.solar_theta = xr->solar_theta;

   outrec.write_nominal_wavelength_grid = Write_Nominal_Wavelength_Grid;

   if (0 != out->out_write_rec (out, xr->index, &outrec))
     goto return_status;

   if (0 != trend_collect_write_record (xr->tr))
     goto return_status;

   Write_Nominal_Wavelength_Grid = 0;

   status = 0;
return_status:
   (void) trend_collect_set_active_record (tr_old);
   sdt_free (outrec.uv);
   sdt_free (outrec.vis);

   return status;
}

#define QUEUE_DEPTH  3

typedef struct
{
   Image_Type *img;
   Exprec_Meta_Type *xr;
}
Queue_Entry_Type;
typedef struct
{
   Queue_Entry_Type *items[QUEUE_DEPTH];
   int num_queued;
}
Queue_Type;

static void free_queue_entry (Queue_Entry_Type *entry)
{
   if (entry == NULL)
     return;
   image_free (entry->img);
   entry->img = NULL;
   entry->xr = NULL;
   FREE(entry);
}

static Queue_Entry_Type *new_queue_entry (Exprec_Meta_Type *xr)
{
   Queue_Entry_Type *entry = NULL;
   if (NULL == (entry = (Queue_Entry_Type *)MALLOC (sizeof *entry)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   entry->xr = xr;

   if (NULL == (entry->img = image_dup (xr->exprec->img)))
     {
        free_queue_entry (entry);
        return NULL;
     }

   return entry;
}

static void queue_init (Queue_Type *q)
{
   memset ((char *)q->items, 0, sizeof (q->items));
   q->num_queued = 0;
}

static Queue_Entry_Type *queue_shift (Queue_Type *q)
{
   Queue_Entry_Type *oldest = q->items[0];
   int i;

   for (i = 1; i < QUEUE_DEPTH; i++)
     {
        q->items[i-1] = q->items[i];
     }

   return oldest;
}

static int queue_push (Queue_Type *q, Exprec_Meta_Type *xr, Granule_Type *gr)
{
   Queue_Entry_Type *qnew = NULL;
   Queue_Entry_Type *oldest;

   if (NULL == (qnew = new_queue_entry (xr)))
     return -1;

   oldest = queue_shift (q);
   if (oldest)
     {
        free_exprec_meta (oldest->xr, gr);
        free_queue_entry (oldest);
     }

   q->items[QUEUE_DEPTH-1] = qnew;
   q->num_queued += 1;
   if (q->num_queued > QUEUE_DEPTH)
     q->num_queued = QUEUE_DEPTH;

   return 0;
}

static void queue_empty (Queue_Type *q, Granule_Type *gr)
{
   int i;

   if (q == NULL)
     return;

   for (i = 0; i < q->num_queued; i++)
     {
        Queue_Entry_Type *entry = q->items[i];
        if (entry) free_exprec_meta (entry->xr, gr);
        free_queue_entry (entry);
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
   Queue_Entry_Type *prev, *xr, *next;

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
        img_ref = next->img;
     }
   else
     {
        /* here's where we actually need the space img_ref points to */
        make_transient_ref_img (prev->img, next->img, img_ref);
     }

   /* since xr->img is a duplicate, we must transfer any flag bits that get set */
   if ((0 != pqft->pqf_flag_transients (pqft, bpixmap->bits, img_ref, xr->img))
       || (0 != image_transfer_pqf (xr->img, xr->xr->exprec->img)))
     return -1;

   /* When exprec_index==num_exprecs-1, process one more to finish: */
   if (exprec_index == num_exprecs-1)
     {
        /* queue contains {..., img(n-2), img(n-1)} */
        prev = q->items[1];
        img_ref = prev->img;
        xr = q->items[2];
        /* since xr->img is a duplicate, we must transfer any flag bits that get set */
        if ((0 != pqft->pqf_flag_transients (pqft, bpixmap->bits, img_ref, xr->img))
            || (0 != image_transfer_pqf (xr->img, xr->xr->exprec->img)))
          return -1;
     }

   return 0;
}

static int derive_photons (config_t *cfg, const Control_Type *ctrl, Process_Control_Type *pct,
                           Granule_Type *gr, TIO_Meta_Type *meta, Trend_File_Type *tft)
{
   Queue_Type exprec_queue = {0};
   CCD_Type *ccd = NULL;
   Instr_Type *instr = NULL;
   Calibration_Type *cal = NULL;
   Exprec_Meta_Type *xr = NULL;
   Pixelqf_Type *pqft = NULL;
   Badpix_Map_Type *bpixmap = NULL;
   Dark_Type *drk = NULL;
   Output_Type *out = NULL;
   Image_Type *tmp_img = NULL;
   Solar_Geom_Type *sgt = NULL;
   int num_serial_active_full, num_parallel_active_full, flag_transients;
   int ixr, num_exprecs, exposure_type, scan_type, ncid_from, ncid_to;
   int processing_version;
   int status = -1;

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   num_exprecs = gr->granule_num_exprecs(gr);
   if (ctrl->limit_num_granules < num_exprecs)
     num_exprecs = ctrl->limit_num_granules;

   if (NULL == (ccd = ccd_init (cfg, meta)))
     goto return_status;

   if (NULL == (pqft = pixelqf_init (cfg)))
     goto return_status;

   if (NULL == (instr = instr_open (ctrl->instr_status_file, ctrl->instr_glob,
                                    gr->granule_tstart(gr), gr->granule_tend(gr), meta)))
     goto return_status;

   if (NULL == (bpixmap = bpix_read (ctrl->bpix_file)))
     goto return_status;

   if (0)
     {
        if (0 != meta_record_basename (meta, ctrl->bpix_file))
          goto return_status;
     }

   if (NULL == (cal = sensorcal_init (cfg, meta, ctrl->irr_file, exposure_type)))
     goto return_status;
   if (0 != meta_record_basename (meta, ctrl->irr_file))
     goto return_status;

   if (NULL == (drk = drk_init (cfg)))
     goto return_status;

   if (0 != drk->drk_open (drk, ctrl->dark_file))
     goto return_status;
   if (0 != meta_record_basename (meta, ctrl->dark_file))
     goto return_status;

   /* Allocate reusable scratch space to hold a full, trimmed CCD image */
   ccd->ccd_active_image_dims (ccd, &num_parallel_active_full, &num_serial_active_full);
   if (NULL == (tmp_img = image_new (num_parallel_active_full, num_serial_active_full)))
     goto return_status;

   if (0 != create_diagnostic_file (ctrl, num_parallel_active_full, num_serial_active_full,
                                    (cal->cal_straylight_correction != NULL),
                                    EXPREC_TYPE_IS_IRRADIANCE(exposure_type)))
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
   ncid_to = out->out_root_ncid (out);

   processing_version = process_get_version();
   if (0 != TIO_put_att (ncid_to, NC_GLOBAL, "processing_version", NC_INT, 1, &processing_version))
     goto return_status;

   switch (exposure_type)
     {
      case EXPREC_TYPE_RAD:
      case EXPREC_TYPE_RAD_TWI:
        if (0 != TIO_copy_granule_ident (ncid_from, ncid_to))
          goto return_status;
        if (0 != tio_copy_granule_flag_var (ncid_from, ncid_to))
          goto return_status;
        /* Copy the scan_type attribute to the Level 1 file.
         * It's not in the granule ident struct for now, because putting it there would
         * be too much trouble.
         */
        if ((0 != TIO_get_att (ncid_from, NC_GLOBAL, "scan_type", NC_INT, &scan_type))
            ||(0 != TIO_put_att (ncid_to, NC_GLOBAL, "scan_type", NC_INT, 1, &scan_type)))
          goto return_status;
        break;

      case EXPREC_TYPE_IRR_WRK:
        /* drop */
      case EXPREC_TYPE_IRR_REF:
        /* nothing yet */
        break;

      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported exposure record type = %d", __func__, exposure_type);
        goto return_status;
     }

   /* we use this for both radiance and irradiance */
   if (NULL == (sgt = solar_geom_init (cfg)))
     goto return_status;

   if ((flag_transients = enable_state_query_bool (ENABLE_TRANSIENTS)) < 0)
     goto return_status;
   if ((flag_transients != 0) && (num_exprecs < 3))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "%s: flag_transients is ON, but we processing too few exposure records (>= 3 are required)", __func__);
        flag_transients = 0;
     }

   Write_Nominal_Wavelength_Grid = 1;
   queue_init (&exprec_queue);

   for (ixr = 0; ixr < num_exprecs; ixr++)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 1, "exposure record %3d/%d", (ixr+1), num_exprecs);

        if (NULL == (xr = alloc_exprec_meta (tft)))
          goto return_status;

        xr->index = ixr;

        if (NULL == (xr->exprec = gr->granule_read_exprec_by_index (gr, ixr, NULL)))
          goto return_status;

        if (xr->exprec->num_coadds == TIO_FILL_UINT)
          {
             tell_vwarn (0, "%s: exposure record %d/%d (missing data)", __func__, (ixr+1), num_exprecs);
             free_exprec_meta (xr, gr);
             xr = NULL;
             continue;
          }

        trend_collect_set_active_record (xr->tr);

        if (0 != trend_collect_time (xr->exprec->start_time, ixr))
          goto return_status;

        if (-1 == validate_exposure_type (exposure_type, xr->exprec->exposure_type))
          goto return_status;
        if (0 != ccd->ccd_apply_pixel_quality_flags (ccd, xr->exprec->img,
                                                     bpixmap->bits, bpixmap->num_rows, bpixmap->num_cols))
          goto return_status;

        if (-1 == compute_current_and_trim (ccd, instr, pqft, pct, xr))
          goto return_status;

        if (0 != dark_subtract (drk, xr, tmp_img))
          return -1;

        if (flag_transients == 0)
          {
             if (0 != radcal_and_output (out, cal, sgt, xr))
               goto return_status;
             free_exprec_meta (xr, gr);
             xr = NULL;
          }
        else
          {
             /* To look for transients, we need at least 2 pixel current (e-/sec)
              * images queued. To queue these images, we need to make a duplicate copy
              * of the pixel current image because the radiometric correction
              * will occur in-place in the next step.  Along with the duplicate
              * pixel current image, we'll also need the Exprec_Meta_Type pointer
              * for each entry when it comes out of the queue.  queue_push allocates
              * a new queue entry containing both of these objects, and deletes the
              * old queue entry that's no longer needed.  Sure, it's a little complicated,
              * but we only need 3 frames in memory at once, instead of reading the
              * entire granule, which may contain >100 frames.
              */
             if (0 != queue_push (&exprec_queue, xr, gr))
               goto return_status;
             if (exprec_queue.num_queued < 2)
               continue;
             if (0 != flag_transients1 (pqft, bpixmap, num_exprecs, &exprec_queue, ixr, tmp_img))
               goto return_status;
             /* The previous frame, ixr-1, is now ready to continue processing */
             if (0 != radcal_and_output (out, cal, sgt, exprec_queue.items[1]->xr))
               goto return_status;
          }
     }

   if (flag_transients)
     {
        /* Process the last entry in the queue, exprec[num_exprecs-1] */
        if (0 != radcal_and_output (out, cal, sgt, exprec_queue.items[2]->xr))
          goto return_status;
     }

   if (0 != write_std_metadata (ncid_to, exposure_type, ctrl, meta, ncid_from))
     goto return_status;

   if (0 != out->out_finalize (out))
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
   close_diagnostic_file ();

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
   Trend_File_Type *tft = NULL;
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
   if (0 != meta_record_basename (meta, ctrl->input_file))
     goto return_status;

   if (0 != gr->granule_type (gr, &exposure_type))
     goto return_status;

   if (ctrl->trend_file != NULL)
     {
        if (NULL == (tft = trend_collect_open (ctrl->trend_file, exposure_type)))
          goto return_status;
     }

   (void) set_diagnostic_controls_index (ctrl->diagnostic_index);

   /* It's not expected that data from a linearity sweep will be processed
    * to Level 1, but if it is, then we'll treat it exactly the same
    * as any other dark or irradiance measurement */

   switch (exposure_type)
     {
      case EXPREC_TYPE_DARK:
      case EXPREC_TYPE_LIN_DARK:
        status = derive_current (cfg, ctrl, &pct, gr, meta, tft);
        break;

      case EXPREC_TYPE_IRR_WRK:
      case EXPREC_TYPE_IRR_REF:
      case EXPREC_TYPE_LIN_IRR:
      case EXPREC_TYPE_RAD:
      case EXPREC_TYPE_RAD_TWI:
        /* For irradiances, we'll need to compute the solar illumination geometry.
         * For radiances, we'll compute the earth-sun distance.
         */
        tbeg = gr->granule_tstart (gr);
        tend = gr->granule_tend (gr);
        if (0 != init_solsys_ephem (cfg, tbeg, tend, &eph))
          goto return_status;
        status = derive_photons (cfg, ctrl, &pct, gr, meta, tft);
        break;

      default:
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: input granule contains exposures of unknown or unexpected type (type=%d)",
                     __func__, exposure_type);
        break;
     }

return_status:

   if (gr) gr->granule_close (gr);
   (void) trend_collect_close (tft);
   tio_meta_close (meta);
   (void) ephem_close (&eph);

   return status;
}
