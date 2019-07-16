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

static void drk_close (Dark_Type *drk)
{
   if (drk == NULL)
     return;
   image_free(drk->image);
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
   Image_Type *img = NULL;
   int start[3], count[3];

   if (0 != TIO_inq_var (ncid, "image", &info))
     return NULL;

   if (NULL == (img = image_new (info.dimlens[1], info.dimlens[2])))
     return NULL;

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = info.dimlens[1];
   count[2] = info.dimlens[2];

   if ((0 != TIO_get_var_section (ncid, "image", start, count, TIO_FLOAT, img->pixels))
       || (0 != TIO_get_var_section (ncid, TEMPO_VAR_PQF, start, count, TIO_USHORT, img->pixel_quality_flags)))
     {
        image_free(img);
        return NULL;
     }

   img->image_type = IMAGE_TYPE_ACTIVE;

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
