#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "granule.h"
#include "instr.h"
#include "dark.h"
#include "process.h"
#include "_process.h"
#include "sensorcal.h"
#include "current.h"

#include "util.h"

/* CCD integration types */
enum
{
   INT_TYPE_NOMINAL = 0,
   INT_TYPE_SHORT = 1,
   INT_TYPE_LONG = 2,
   INT_TYPE_DARK = 3
};

typedef struct
{
   const char *name;
   const char *text;
}
Text_Attr_Type;

static int define_text_attrs (int grp, int varid, const Text_Attr_Type *attrs)
{
   const Text_Attr_Type *a;

   for (a = attrs; a->name != NULL; a++)
     {
        size_t len = strlen(a->text) + 1;
        if (0 != TIO_put_att (grp, varid, a->name, TIO_CHAR, len, a->text))
          return -1;
     }

   return 0;
}

static int create_current_vars (int ncid, int num_times, int num_rows, int num_cols,
                                int exposure_type)
{
   int dimid_time, dimid_row, dimid_col, dimid_quad;
   int varid_time, varid_img, varid_pqf, varid_fpa_temp;
   int varid_fpe_temp, varid_exptime, varid_sdc;
   int varid_num_hot, varid_num_cold, varid_mean_dark;
   int varid_exptime_per_coadd, varid_ccd_int_type;
   int sdc_dimids[2];
   int shuffle = 1;
   int deflate = 1;
   int deflate_level = 1;
   int storage = NC_CHUNKED;
   size_t chunksizes[3];
   int img_dimids[3], start[3], count[3];
   int k, len_pqf_slab;
   unsigned short *pqf_missing = NULL;
   const Text_Attr_Type time_attrs[] =
     {
        {"long_name", "exposure start time"},
        {NULL, NULL}
     };
   const Text_Attr_Type img_attrs[] =
     {
        {"units", "electrons/sec"},
        {"long_name", "pixel current"},
        {NULL, NULL}
     };
   const Text_Attr_Type pqf_attrs[] =
     {
        {"long_name", "pixel quality flag"},
        {NULL, NULL}
     };
   const Text_Attr_Type fpa_temp_attrs[] =
     {
        {"units", "C"},
        {"long_name", "focal plane array (FPA) temperature"},
        {NULL, NULL}
     };
   const Text_Attr_Type fpe_temp_attrs[] =
     {
        {"units", "C"},
        {"long_name", "focal plane electronics (FPE) temperature"},
        {NULL, NULL}
     };
   const Text_Attr_Type exptime_attrs[] =
     {
        {"units", "seconds"},
        {"long_name", "exposure time duration"},
        {NULL, NULL}
     };
   const Text_Attr_Type ccd_int_type_attrs[] =
     {
        {"long_name", "CCD integration type"},
        {"comment",  "CCD integration type: NOMINAL=0, SHORT_INT=1, LONG_INT=2, DARK_INT=3"},
        {NULL, NULL}
     };
   const Text_Attr_Type exptime_per_coadd_attrs[] =
     {
        {"units", "seconds"},
        {"long_name", "exposure time duration per coadd"},
        {NULL, NULL}
     };
   const Text_Attr_Type sdc_attrs[] =
     {
        {"units", "electrons/sec"},
        {"long_name", "mean storage region dark current"},
        {"comment", "mean storage region dark current in each quadrant; A,B,C,D"},
        {NULL, NULL}
     };
   const Text_Attr_Type num_hot_attrs[] =
     {
        {"long_name", "number of hot pixels"},
        {"comment", "number of hot pixels in each quadrant; A,B,C,D"},
        {NULL, NULL}
     };
   const Text_Attr_Type num_cold_attrs[] =
     {
        {"long_name", "number of cold pixels"},
        {"comment", "number of cold pixels in each quadrant; A,B,C,D"},
        {NULL, NULL}
     };
   const Text_Attr_Type mean_dark_attrs[] =
     {
        {"units", "electrons/sec"},
        {"long_name", "mean dark current"},
        {"comment", "mean dark current in each quadrant; A,B,C,D (hot/cold pixels excluded)"},
        {NULL, NULL}
     };

   if ((0 != TIO_def_dim (ncid, "time", num_times, &dimid_time))
       || (0 != TIO_def_dim (ncid, "row", num_rows, &dimid_row))
       || (0 != TIO_def_dim (ncid, "col", num_cols, &dimid_col))
       || (0 != TIO_def_dim (ncid, "quad", 4, &dimid_quad)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining file array dimensions", __func__);
        return -1;
     }

   if ((0 != TIO_def_var (ncid, "image_start_time", TIO_DOUBLE, 1, &dimid_time, &varid_time))
       || (0 != define_text_attrs (ncid, varid_time, time_attrs))
       || (0 != tio_write_timestamp_unit_string (ncid, "image_start_time")))
     return -1;

   if ((0 != TIO_def_var (ncid, "fpa_temp", TIO_FLOAT, 1, &dimid_time, &varid_fpa_temp))
       || (0 != define_text_attrs (ncid, varid_fpa_temp, fpa_temp_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, "fpe_temp", TIO_FLOAT, 1, &dimid_time, &varid_fpe_temp))
       || (0 != define_text_attrs (ncid, varid_fpe_temp, fpe_temp_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, "ccd_int_type", TIO_INT, 1, &dimid_time, &varid_ccd_int_type))
       || (0 != define_text_attrs (ncid, varid_ccd_int_type, ccd_int_type_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, TEMPO_VAR_EXPOSURE_TIME, TIO_FLOAT, 1, &dimid_time, &varid_exptime))
       || (0 != define_text_attrs (ncid, varid_exptime, exptime_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, "exposure_time_per_coadd", TIO_FLOAT, 1, &dimid_time, &varid_exptime_per_coadd))
       || (0 != define_text_attrs (ncid, varid_exptime_per_coadd, exptime_per_coadd_attrs)))
     return -1;

   sdc_dimids[0] = dimid_time;
   sdc_dimids[1] = dimid_quad;

   if ((0 != TIO_def_var (ncid, "mean_sdc", TIO_FLOAT, 2, sdc_dimids, &varid_sdc))
       || (0 != define_text_attrs (ncid, varid_sdc, sdc_attrs)))
     return -1;

   if (EXPREC_TYPE_IS_DARK(exposure_type))
     {
        if ((0 != TIO_def_var (ncid, "mean_dark_current", TIO_FLOAT, 2, sdc_dimids, &varid_mean_dark))
            || (0 != define_text_attrs (ncid, varid_mean_dark, mean_dark_attrs)))
          return -1;

        if ((0 != TIO_def_var (ncid, "num_hot_pixels", TIO_INT, 2, sdc_dimids, &varid_num_hot))
            || (0 != define_text_attrs (ncid, varid_num_hot, num_hot_attrs)))
          return -1;

        if ((0 != TIO_def_var (ncid, "num_cold_pixels", TIO_INT, 2, sdc_dimids, &varid_num_cold))
            || (0 != define_text_attrs (ncid, varid_num_cold, num_cold_attrs)))
          return -1;
     }

   img_dimids[0] = dimid_time;
   img_dimids[1] = dimid_row;
   img_dimids[2] = dimid_col;

   chunksizes[0] = 1;
   chunksizes[1] = 128;
   chunksizes[2] = num_cols / 2;

   if ((0 != TIO_def_var (ncid, "image", TIO_FLOAT, 3, img_dimids, &varid_img))
       || (0 != TIO_def_var_chunking (ncid, varid_img, storage, chunksizes))
       || (0 != TIO_def_var_deflate (ncid, varid_img, shuffle, deflate, deflate_level))
       || (0 != define_text_attrs (ncid, varid_img, img_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, TEMPO_VAR_PQF, TIO_UINT, 3, img_dimids, &varid_pqf))
       || (0 != TIO_def_var_chunking (ncid, varid_pqf, storage, chunksizes))
       || (0 != TIO_def_var_deflate (ncid, varid_pqf, shuffle, deflate, deflate_level))
       || (0 != define_text_attrs (ncid, varid_pqf, pqf_attrs)))
     return -1;

   /* To be certain that all missing data pixels are always labeled as such,
    * we initialize all the pixel_quality_flag values to indicate missing data.
    */

   count[0] = 1;
   count[1] = num_rows;
   count[2] = num_cols;

   start[1] = 0;
   start[2] = 0;

   len_pqf_slab = num_rows * num_cols;
   if (NULL == (pqf_missing = (unsigned short *)MALLOC (len_pqf_slab * sizeof(unsigned short))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   for (k = 0; k < len_pqf_slab; k++)
     {
        pqf_missing[k] = IMAGE_PQF_MISSING_DATA;
     }

   for (k = 0; k < num_times; k++)
     {
        start[0] = k;
        if (0 != TIO_put_var_section (ncid, TEMPO_VAR_PQF, start, count, TIO_USHORT, pqf_missing))
          return -1;
     }

   FREE(pqf_missing);

   return 0;
}

static int write_image_at_index (int ncid, int index, const Image_Type *img)
{
   int start[3], count[3];

   start[0] = index;
   start[1] = 0;
   start[2] = 0;

   count[0] = 1;
   count[1] = img->num_rows;
   count[2] = img->num_cols;

   if (0 != TIO_put_var_section (ncid, "image", start, count, TIO_FLOAT, img->pixels))
     return -1;
   if (0 != TIO_put_var_section (ncid, TEMPO_VAR_PQF, start, count, TIO_USHORT, img->pixel_quality_flags))
     return -1;

   return 0;
}

int current_write_exprec (int ncid, const Exprec_Meta_Type *xr)
{
   Granule_Exprec_Type *exprec = xr->exprec;
   double exposure_time_per_frame;
   int start[3], count[3];

   start[0] = xr->index;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;

   exposure_time_per_frame = exprec->exposure_time / exprec->num_coadds;

   if ((0 != TIO_put_var_section (ncid, "image_start_time", start, count, TIO_DOUBLE, &exprec->start_time))
       || (0 != TIO_put_var_section (ncid, "fpa_temp", start, count, TIO_FLOAT, &xr->fpa_temp))
       || (0 != TIO_put_var_section (ncid, "fpe_temp", start, count, TIO_FLOAT, &xr->fpe_temp))
       || (0 != TIO_put_var_section (ncid, "ccd_int_type", start, count, TIO_INT, &exprec->ccd_int_type))
       || (0 != TIO_put_var_section (ncid, TEMPO_VAR_EXPOSURE_TIME, start, count, TIO_DOUBLE, &exprec->exposure_time))
       || (0 != TIO_put_var_section (ncid, "exposure_time_per_coadd", start, count, TIO_DOUBLE, &exposure_time_per_frame)))
     return -1;

   count[0] = 1;
   count[1] = 4;
   if (0 != TIO_put_var_section (ncid, "mean_sdc", start, count, TIO_FLOAT, xr->storage_region_dark))
     return -1;

   if (EXPREC_TYPE_IS_DARK(exprec->exposure_type))
     {
        const Dark_Trend_Type *dtr = &xr->dark_trend;
        if (0 != TIO_put_var_section (ncid, "num_hot_pixels", start, count, TIO_INT, dtr->num_hot_pixels))
          return -1;
        if (0 != TIO_put_var_section (ncid, "num_cold_pixels", start, count, TIO_INT, dtr->num_cold_pixels))
          return -1;
        if (0 != TIO_put_var_section (ncid, "mean_dark_current", start, count, TIO_FLOAT, dtr->mean_dark_current))
          return -1;
     }

   return write_image_at_index (ncid, xr->index, exprec->img);
}

/* (optional) copy dark_type attribute */
static int copy_dark_type (int ncid_in, int ncid_out)
{
   char dark_type[TIO_MAX_SHORT_NAME_LEN];
   int have_dark_type;

   tell_push_queue();
   have_dark_type = (0 == TIO_get_att (ncid_in, NC_GLOBAL, "dark_type", NC_CHAR, dark_type));
   tell_pop_queue(1);

   if (have_dark_type)
     {
        if (0 != TIO_put_att (ncid_out, NC_GLOBAL, "dark_type", NC_CHAR, 1+strlen(dark_type), dark_type))
          return -1;
     }

   return 0;
}

int current_create_file_of_type (Granule_Type *gr, const char *output_file,
                                 int num_times, int num_rows, int num_cols,
                                 int *pncid, int *pgrp)
{
   const char *product_type;
   int exposure_type, ncid, grp, want_average = 0;
   int status = -1;

   if (0 != gr->granule_type (gr, &exposure_type))
     return -1;

   switch (exposure_type)
     {
      case EXPREC_TYPE_DARK:
        product_type = TEMPO_PROD_TYPE_DRK;
        want_average = 1;
        break;
      case EXPREC_TYPE_LIN_DARK:
        product_type = TEMPO_PROD_TYPE_DRK_LIN;
        break;
      case EXPREC_TYPE_LIN_IRR:
        product_type = TEMPO_PROD_TYPE_IRR_LIN;
        break;
      default:
        tell_verror (TELL_RUNTIME_ERROR, "%s: unsupported exposure record type = %d", __func__, exposure_type);
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "Opening output file: %s", output_file);
   if (0 != TIO_create (output_file, NC_NETCDF4, &ncid))
     return -1;
   if (0 != tio_history_append_cmdline (ncid))
     goto return_status;
   if (0 != tio_write_epoch_timestamp (ncid, NC_GLOBAL))
     goto return_status;
   if (0 != TIO_copy_granule_ident (gr->granule_ncid(gr), ncid))
     goto return_status;
   if (0 != TIO_label_product (ncid, product_type, 1, process_get_version()))
     goto return_status;

   if (exposure_type == EXPREC_TYPE_DARK)
     {
        if (0 != copy_dark_type (gr->granule_ncid(gr), ncid))
          goto return_status;
     }

   if (0 == want_average)
     {
        if (0 != create_current_vars (ncid, num_times, num_rows, num_cols, exposure_type))
          goto return_status;
        *pgrp = ncid;
     }
   else
     {
        /* top-level group will hold the average over non-DARK_INT frames */
        if (0 != create_current_vars (ncid, 1, num_rows, num_cols, exposure_type))
          goto return_status;
        if (0 != TIO_def_grp (ncid, "frames", &grp))
          goto return_status;
        if (0 != create_current_vars (grp, num_times, num_rows, num_cols, exposure_type))
          goto return_status;
        *pgrp = grp;
     }

   *pncid = ncid;
   status = 0;
return_status:

   if (status)
     {
        TIO_close (ncid);
     }

   return status;
}

static int image_add_weighted (Image_Type *img, Image_Type *weights, double weight, const Image_Type *tmp)
{
   int p, pb, pe, s, sb, se; // n = img->num_rows * img->num_cols;

   pb = 0;
   pe = img->num_rows;
   sb = 0;
   se = img->num_cols;

   for (p = pb; p < pe; p++)
     {
        const Image_Pqf_Bitmap_Type *pqf = tmp->pixel_quality_flags + p * img->num_cols;
        const Image_Pixel_Type *pix = tmp->pixels + p * img->num_cols;
        Image_Pixel_Type *wt = weights->pixels + p * img->num_cols;
        Image_Pixel_Type *sum = img->pixels + p * img->num_cols;
        Image_Pqf_Bitmap_Type *sum_pqf = img->pixel_quality_flags + p * img->num_cols;

        for (s = sb; s < se; s++)
          {
             if ((s < 5) || (s >= (se-7)))
               {
                  if ((pqf[s] == 1) && (pix[s] != IMAGE_PIXEL_FILL_VALUE))
                    {
                      if (sum[s] == IMAGE_PIXEL_FILL_VALUE)
                        {  /* first contribution to the total at this pixel */
                            sum[s] = pix[s] * weight;
                            wt[s] = weight;
                            sum_pqf[s] = 0;
                        }
                      else
                        {
                            sum[s] += pix[s] * weight;
                            wt[s] += weight;
                        }
                    }
               }
             else
               {
                  if (pqf[s] == 0)
                    {
                      if (sum[s] == IMAGE_PIXEL_FILL_VALUE)
                        {  /* first contribution to the total at this pixel */
                            sum[s] = pix[s] * weight;
                            wt[s] = weight;
                            sum_pqf[s] = 0;
                        }
                      else
                        {
                            sum[s] += pix[s] * weight;
                            wt[s] += weight;
                        }
                    }
               }
          }
     }

   return 0;
}

static int image_divide (Image_Type *img, const Image_Type *denom)
{
   const Image_Pixel_Type *den = denom->pixels;
   Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags;
   Image_Pixel_Type *pix = img->pixels;
   size_t i, n = img->num_rows * img->num_cols;

   for (i = 0; i < n; i++)
     {
        if (pqf[i] == 0)
          {
             if (isfinite(den[i]) && (den[i] != 0.0))
               {
                  pix[i] /= den[i];
               }
             else
               {
                  pix[i] = IMAGE_PIXEL_FILL_VALUE;
                  pqf[i] |= IMAGE_PQF_BAD_PIXEL;
               }
          }
     }

   return 0;
}

static int get_averaging_weights (int ncid, int num, double *wt)
{
   int k, start, count, num_zero_exposures;

   /* Assume all exposure times are >= 0.
    * If all exposure_times are zero, then weight equally when computing the mean.
    * If any exposure_times are non-zero, then weight by exposure time so that the
    *                               terms with exposure_times=0 contribute nothing.
    * If the data are packaged so that all images in a single file have the same exposure
    * time, then this should do the right thing whether the exposure time is zero or not.
    * If exposure times vary within a single file, then this also seems the right approach.
    */

   start = 0;
   count = num;
   if (0 != TIO_get_var_section (ncid, "exposure_time", &start, &count, TIO_DOUBLE, wt))
     return -1;

   num_zero_exposures = 0;
   for (k = 0; k < num; k++)
     {
        if (wt[k] == 0.0) num_zero_exposures++;
     }

   if (num == num_zero_exposures)
     {
        for (k = 0; k < num; k++)
          {
             wt[k] = 1.0;
          }
     }

   return 0;
}

int current_image_read (int ncid, int k, Image_Type *img)
{
   int start[3], count[3];

   start[0] = k;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = img->num_rows;
   count[2] = img->num_cols;

   if ((0 != TIO_get_var_section (ncid, "image", start, count, TIO_FLOAT, img->pixels))
       || (0 != TIO_get_var_section (ncid, "pixel_quality_flag", start, count, TIO_USHORT, img->pixel_quality_flags)))
     {
        return -1;
     }

   img->image_type = IMAGE_TYPE_ACTIVE;

   return 0;
}

static Image_Type *compute_image_mean (int ncid, const int *filter, unsigned int num_filter)
{
   TIO_Var_Info_Type info = {0};
   Image_Type *img = NULL;
   Image_Type *tmp = NULL;
   Image_Type *weights = NULL;
   double *wt = NULL;
   int k, num, status = -1;

   if (0 != TIO_inq_var (ncid, "image", &info))
     return NULL;

   if (num_filter != info.dimlens[0])
     goto return_status;

   if ((NULL == (img = image_new (info.dimlens[1], info.dimlens[2])))
       || (NULL == (weights = image_dup (img)))
       || (NULL == (tmp = image_dup (img))))
     goto return_status;

   num = info.dimlens[0];

   tell_vlog (TELL_MSGTYPE_INFO, 1, "averaging %d dark frames...", num);

   if (NULL == (wt = (double *)MALLOC (num * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   if (0 != get_averaging_weights (ncid, num, wt))
     goto return_status;

   image_set (img, IMAGE_PIXEL_FILL_VALUE, IMAGE_PQF_MISSING_DATA);

   for (k = 0; k < num; k++)
     {
        if (filter[k] == 0)
          continue;
        if (0 != current_image_read (ncid, k, tmp))
          goto return_status;
        if (0 != image_add_weighted (img, weights, wt[k], tmp))
          goto return_status;
     }

   if (0 != image_divide (img, weights))
     goto return_status;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "done");

   status = 0;

return_status:
   image_free (tmp);
   image_free (weights);
   if (status)
     {
        image_free (img);
        img = NULL;
     }
   FREE(wt);

   return img;
}

static double d_mean_good (int n, const int *filter, double *x, double bad)
{
   double avg = 0.0;
   int i, k;

   if (n <= 0)
     return 0.0;

   k = 0;

   for (i = 0; i < n; i++)
     {
        if ((filter != NULL) && (filter[i] == 0))
          continue;
        if (x[i] != bad)
          {
             avg += x[i];
             k += 1;
          }
     }

   if (k > 0)
     {
        avg /= k;
     }
   else avg = bad;

   return avg;
}

static float f_mean_good (int n, const int *filter, float *x, float bad)
{
   double avg = 0.0;
   int i, k;

   if (n <= 0)
     return 0.0;

   k = 0;

   for (i = 0; i < n; i++)
     {
        if ((filter != NULL) && (filter[i] == 0))
          continue;
        if (x[i] != bad)
          {
             avg += x[i];
             k += 1;
          }
     }

   if (k > 0)
     {
        avg /= k;
     }
   else avg = bad;

   return (float) avg;
}

int current_write_mean_dark_current (int ncid)
{
   TIO_Var_Info_Type info = {0};
   Dark_Trend_Type avg_dtr = {0};
   float avg_sdc[4];
   const char *dbl_var_names[] = {"image_start_time", TEMPO_VAR_EXPOSURE_TIME, "exposure_time_per_coadd", NULL};
   const char *flt_var_names[] = {"fpa_temp", "fpe_temp", NULL};
   const char **var_name;
   Image_Type *img = NULL;
   int *filter = NULL;
   double *d_tmp = NULL;
   double d_mean_value;
   float *f_tmp = NULL;
   float f_mean_value;
   int i, j, grp, num_frames, num_active_frames, start[2], count[2], one=1;
   int status = -1;

   if (0 != TIO_inq_grp (ncid, "frames", &grp))
     return -1;

   /* Average over active-region frames only.
    * Ignore DARK_INT frames, which are images of the storage region */
   if (0 != TIO_inq_var (grp, "ccd_int_type", &info))
     goto return_status;
   num_frames = info.dimlens[0];
   if (NULL == (filter = (int *)MALLOC (num_frames * sizeof(int))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }
   start[0] = 0;
   count[0] = num_frames;
   if (0 != TIO_get_var_section (grp, "ccd_int_type", start, count, TIO_INT, filter))
     goto return_status;

   /* The non-DARK_INT frames should all be the same integration time,
    * but this is not explicitly checked here */
   num_active_frames = 0;
   for (i = 0; i < num_frames; i++)
     {
        if (filter[i] != INT_TYPE_DARK)
          {
             filter[i] = 1;
             num_active_frames++;
          }
        else filter[i] = 0;
     }

   if (NULL == (img = compute_image_mean (grp, filter, num_frames)))
     goto return_status;

   if (0 != write_image_at_index (ncid, 0, img))
     goto return_status;

   if ((NULL == (d_tmp = (double *)MALLOC (num_frames * sizeof (double))))
       || (NULL == (f_tmp = (float *)MALLOC (num_frames * sizeof(float)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_status;
     }

   start[0] = 0;
   count[0] = num_frames;

   for (var_name = dbl_var_names; (var_name != NULL) && (*var_name != NULL); var_name++)
     {
        if (0 != TIO_get_var_section (grp, *var_name, start, count, TIO_DOUBLE, d_tmp))
          goto return_status;
        d_mean_value = d_mean_good (num_frames, filter, d_tmp, TIO_FILL_DOUBLE);
        if (0 != TIO_put_var_section (ncid, *var_name, start, &one, TIO_DOUBLE, &d_mean_value))
          goto return_status;
     }

   for (var_name = flt_var_names; (var_name != NULL) && (*var_name != NULL); var_name++)
     {
        if (0 != TIO_get_var_section (grp, *var_name, start, count, TIO_FLOAT, f_tmp))
          goto return_status;
        f_mean_value = f_mean_good (num_frames, filter, f_tmp, TIO_FILL_FLOAT);
        if (0 != TIO_put_var_section (ncid, *var_name, start, &one, TIO_FLOAT, &f_mean_value))
          goto return_status;
     }

   memset ((char *)avg_dtr.num_hot_pixels, 0, 4 * sizeof(int));
   memset ((char *)avg_dtr.num_cold_pixels, 0, 4 * sizeof(int));
   memset ((char *)avg_dtr.mean_dark_current, 0, 4 * sizeof(float));
   memset ((char *)avg_sdc, 0, 4 * sizeof(float));

   for (i = 0; i < num_frames; i++)
     {
        Dark_Trend_Type dtr;
        float sdc[4];
        if (filter[i] == 0)
          continue;
        start[0] = i;
        start[1] = 0;
        count[0] = 1;
        count[1] = 4;
        if (0 != TIO_get_var_section (grp, "num_hot_pixels", start, count, TIO_INT, dtr.num_hot_pixels))
          goto return_status;
        if (0 != TIO_get_var_section (grp, "num_cold_pixels", start, count, TIO_INT, dtr.num_cold_pixels))
          goto return_status;
        if (0 != TIO_get_var_section (grp, "mean_dark_current", start, count, TIO_FLOAT, dtr.mean_dark_current))
          goto return_status;
        if (0 != TIO_get_var_section (grp, "mean_sdc", start, count, TIO_FLOAT, sdc))
          goto return_status;
        for (j = 0; j < 4; j++)
          {
             avg_dtr.num_hot_pixels[j] += dtr.num_hot_pixels[j];
             avg_dtr.num_cold_pixels[j] += dtr.num_cold_pixels[j];
             avg_dtr.mean_dark_current[j] += dtr.mean_dark_current[j];
             avg_sdc[j] += sdc[j];
          }
     }

   for (j = 0; j < 4; j++)
     {
        avg_dtr.num_hot_pixels[j] /= num_active_frames;
        avg_dtr.num_cold_pixels[j] /= num_active_frames;
        avg_dtr.mean_dark_current[j] /= num_active_frames;
        avg_sdc[j] /= num_active_frames;
     }

   start[0] = 0;
   start[1] = 0;
   count[0] = 1;
   count[1] = 4;
   if (0 != TIO_put_var_section (ncid, "num_hot_pixels", start, count, TIO_INT, avg_dtr.num_hot_pixels))
     goto return_status;
   if (0 != TIO_put_var_section (ncid, "num_cold_pixels", start, count, TIO_INT, avg_dtr.num_cold_pixels))
     goto return_status;
   if (0 != TIO_put_var_section (ncid, "mean_dark_current", start, count, TIO_FLOAT, avg_dtr.mean_dark_current))
     goto return_status;
   if (0 != TIO_put_var_section (ncid, "mean_sdc", start, count, TIO_FLOAT, avg_sdc))
     goto return_status;

   status = 0;
return_status:
   FREE(filter);
   FREE(d_tmp);
   FREE(f_tmp);
   image_free (img);

   return status;
}
