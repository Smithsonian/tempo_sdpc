#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "image.h"

#define DARK_PRIVATE_DATA \
   char *dark_file; \
   Image_Type *image; \
   int ncid;
#include "dark.h"

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

int drk_create_file (int ncid, int num_times, int num_rows, int num_cols)
{
   const char *product_type = TEMPO_PROD_TYPE_DRK;
   int dimid_time, dimid_row, dimid_col, dimid_quad;
   int varid_time, varid_img, varid_pqf, varid_fpa_temp;
   int varid_fpe_temp, varid_exptime, varid_sdc;
   int sdc_dimids[2];
   int shuffle = 1;
   int deflate = 1;
   int deflate_level = 1;
   int storage = NC_CHUNKED;
   size_t chunksizes[3];
   int img_dimids[3];
   const Text_Attr_Type time_attrs[] =
     {
        {"units", "sec since TEMPO epoch"},
        {"comment", "dark collect start time"},
        {NULL, NULL}
     };
   const Text_Attr_Type img_attrs[] =
     {
        {"units", "electrons/sec"},
        {"comment", "dark current"},
        {NULL, NULL}
     };
   const Text_Attr_Type pqf_attrs[] =
     {
        {"comment", "pixel quality flag"},
        {NULL, NULL}
     };
   const Text_Attr_Type fpa_temp_attrs[] =
     {
        {"units", "C"},
        {"comment", "FPA temperature"},
        {NULL, NULL}
     };
   const Text_Attr_Type fpe_temp_attrs[] =
     {
        {"units", "C"},
        {"comment", "FPE temperature"},
        {NULL, NULL}
     };
   const Text_Attr_Type exptime_attrs[] =
     {
        {"units", "seconds"},
        {"comment", "Exposure time"},
        {NULL, NULL}
     };
   const Text_Attr_Type sdc_attrs[] =
     {
        {"units", "electrons/sec"},
        {"comment", "Mean storage region dark current in each quadrant; A,B,C,D"},
        {NULL, NULL}
     };

   if ((0 != TIO_def_dim (ncid, "time", num_times, &dimid_time))
       || (0 != TIO_def_dim (ncid, "row", num_rows, &dimid_row))
       || (0 != TIO_def_dim (ncid, "col", num_cols, &dimid_col))
       || (0 != TIO_def_dim (ncid, "quad", 4, &dimid_quad)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark array dimensions", __func__);
        return -1;
     }

   if (0 != TIO_put_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, strlen(product_type), product_type))
     return -1;

   if ((0 != TIO_def_var (ncid, "image_start_time", TIO_DOUBLE, 1, &dimid_time, &varid_time))
       || (0 != define_text_attrs (ncid, varid_time, time_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, "fpa_temp", TIO_FLOAT, 1, &dimid_time, &varid_fpa_temp))
       || (0 != define_text_attrs (ncid, varid_fpa_temp, fpa_temp_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, "fpe_temp", TIO_FLOAT, 1, &dimid_time, &varid_fpe_temp))
       || (0 != define_text_attrs (ncid, varid_fpe_temp, fpe_temp_attrs)))
     return -1;

   if ((0 != TIO_def_var (ncid, TEMPO_VAR_EXPOSURE_TIME, TIO_FLOAT, 1, &dimid_time, &varid_exptime))
       || (0 != define_text_attrs (ncid, varid_exptime, exptime_attrs)))
     return -1;

   sdc_dimids[0] = dimid_time;
   sdc_dimids[1] = dimid_quad;

   if ((0 != TIO_def_var (ncid, "mean_sdc", TIO_FLOAT, 2, sdc_dimids, &varid_sdc))
       || (0 != define_text_attrs (ncid, varid_sdc, sdc_attrs)))
     return -1;

   img_dimids[0] = dimid_time;
   img_dimids[1] = dimid_row;
   img_dimids[2] = dimid_col;

   chunksizes[0] = 1;
   chunksizes[1] = num_rows / 2;
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

   return 0;
}

static void drk_close (Dark_Type *drk)
{
   if (drk == NULL)
     return;
   FREE(drk->dark_file);
   FREE(drk);
}

static int drk_image_copy (const Dark_Type *drk, const Dark_Lookup_Type *dlt, Image_Type *img)
{
   (void) drk; (void) dlt;
   return image_copy (drk->image, img);
}

static int drk_lookup_exptime (const Dark_Type *drk, const Dark_Lookup_Type *dlt, Image_Type *img)
{
   (void) drk; (void) dlt; (void) img;
   fprintf (stderr, "%s: not implemented\n", __func__);
   return -1;
}

static int drk_lookup_sdc (const Dark_Type *drk, const Dark_Lookup_Type *dlt, Image_Type *img)
{
   (void) drk; (void) dlt; (void) img;
   fprintf (stderr, "%s: not implemented\n", __func__);
   return -1;
}

static int drk_lookup_fptemp (const Dark_Type *drk, const Dark_Lookup_Type *dlt, Image_Type *img)
{
   (void) drk; (void) dlt; (void) img;
   fprintf (stderr, "%s: not implemented\n", __func__);
   return -1;
}

static Image_Type *read_dark_image (int ncid)
{
   TIO_Var_Info_Type info = {0};
   Image_Type *img;
   int start[2], count[2];

   if (0 != TIO_inq_var (ncid, "image", &info))
     return NULL;

   if (NULL == (img = image_new (info.dimlens[0], info.dimlens[1])))
     return NULL;

   start[0] = 0;
   start[1] = 0;
   count[0] = info.dimlens[0];
   count[1] = info.dimlens[1];

   if (0 != TIO_get_var_section (ncid, "image", start, count, TIO_FLOAT, img->pixels))
     {
        image_free(img);
        return NULL;
     }

   return img;
}

Dark_Type *drk_open (const char *path)
{
   Dark_Type *drk = NULL;
   char product_type[TIO_MAX_SHORT_NAME_LEN];
   int method, ncid = 0;
   int status = -1;

   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s for reading",
                     __func__, path);
        return NULL;
     }

   if (0 != TIO_get_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, product_type))
     goto cleanup_and_return;

   if (0 == strcmp (product_type, TEMPO_PROD_TYPE_DRK))
     {
        method = DARK_METHOD_FILE;
     }
   else
     {
        fprintf (stderr, "%s: no support for dark data with product_type=%s\n",
                 __func__, product_type);
        goto cleanup_and_return;
     }

   if (NULL == (drk = (Dark_Type *)MALLOC (sizeof (*drk))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_return;
     }
   memset ((char *)drk, 0, sizeof (*drk));

   if (NULL == (drk->dark_file = strdup (path)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto cleanup_and_return;
     }

   drk->ncid = ncid;
   drk->drk_close = drk_close;

   switch (method)
     {
      case DARK_METHOD_FILE:
        if (NULL == (drk->image = read_dark_image (ncid)))
          goto cleanup_and_return;
        (void) TIO_close (ncid);
        drk->ncid = 0;
        drk->drk_get_image = drk_image_copy;
        break;

      case DARK_METHOD_LOOKUP_EXPTIME:
        drk->drk_get_image = drk_lookup_exptime;
        break;

      case DARK_METHOD_LOOKUP_SDC:
        drk->drk_get_image = drk_lookup_sdc;
        break;

      case DARK_METHOD_LOOKUP_FPTEMP:
        drk->drk_get_image = drk_lookup_fptemp;
        break;

      default:
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: unsupported dark method = %d",
                     __func__, method);
        goto cleanup_and_return;
     }

   status = 0;
cleanup_and_return:
   if (status)
     {
        if (ncid) (void) TIO_close(ncid);
        drk_close (drk);
        drk = NULL;
     }

   return drk;
}
