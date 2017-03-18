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

#define NUM_BANDS 2

typedef struct
{
   const char *name;
   int num_xtrack;
   int num_channels;
   int ybeg;
   int yend;
}
Band_Info_Type;

#define OUTPUT_PRIVATE_DATA \
   Band_Info_Type bands[NUM_BANDS]; \
   Image_Type *img_outbuf; \
   char *file; \
   int exposure_type; \
   int have_dims; \
   int num_recs; \
   int num_xtrack; \
   int num_waves; \
   int formatted_file_exists; \
   int ncid;
#include "output.h"

static int out_set_file (Output_Type *out, const char *file)
{
   if (NULL == (out->file = strdup (file)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   return 0;
}

static int out_set_dims (Output_Type *out, int num_recs,
                         int num_xtrack, int num_waves)
{
   if ((num_recs < 1) || (num_xtrack < 1) || (num_waves < 1))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: invalid dimensions num_recs=%d, num_xtrack=%d, num_waves=%d",
                     __func__, num_recs, num_xtrack, num_waves);
        return -1;
     }

   out->num_recs = num_recs;
   out->num_xtrack = num_xtrack;
   out->num_waves = num_waves;
   out->have_dims = 1;

   return 0;
}

static int out_close (Output_Type *out)
{
   int status = TIO_close (out->ncid);
   out->ncid = 0;
   return status;
}

static int create_file (Output_Type *out)
{
   if (out->file == NULL)
     {
        tell_verror (TELL_USAGE_ERROR,
                     "%s: output file not specified", __func__);
        return -1;
     }

   if (out->have_dims == 0)
     {
        tell_verror (TELL_USAGE_ERROR,
                     "%s: output granule dimensions not specified", __func__);
        return -1;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "creating %s", out->file);

   if (0 != TIO_create (out->file, NC_NETCDF4, &out->ncid))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: creating %s",
                     __func__, out->file);
        return -1;
     }

   return 0;
}

static void define_bands (Output_Type *out, TIO_Scan_Group_Type bands[NUM_BANDS])
{
   Band_Info_Type *out_bands = out->bands;

   out_bands[0].name = TEMPO_BAND_NAME_UV;
   out_bands[0].num_xtrack = out->num_xtrack;
   out_bands[0].num_channels = out->num_waves;
   out_bands[0].ybeg = 0;
   out_bands[0].yend = out->num_waves;

   out_bands[1].name = TEMPO_BAND_NAME_VIS;
   out_bands[1].num_xtrack = out->num_xtrack;
   out_bands[1].num_channels = out->num_waves;
   out_bands[1].ybeg = out->num_waves;
   out_bands[1].yend = out->num_waves*2;

   bands[0].name = TEMPO_BAND_NAME_UV;
   bands[0].num_xtrack = out->num_xtrack;
   bands[0].num_channels = out->num_waves;

   bands[1].name = TEMPO_BAND_NAME_VIS;
   bands[1].num_xtrack = out->num_xtrack;
   bands[1].num_channels = out->num_waves;
}

static int create_file_of_type (Output_Type *out,
                               int (*tmpl_method)(int, size_t, int, TIO_Scan_Group_Type *))
{
   TIO_Scan_Group_Type bands[NUM_BANDS];

   /* Data arrays in the output Level 1 radiance and irradiance files
    * must have wavelength as the fastest varying dimension,
    * so we'll need an output buffer with that shape.
    */
   if (NULL == (out->img_outbuf = image_new (out->num_xtrack, out->num_waves)))
     return -1;

   if (0 != create_file (out))
     return -1;

   define_bands (out, bands);

   if (0 != (*tmpl_method) (out->ncid, out->num_recs, NUM_BANDS, bands))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: creating %s",
                     __func__, out->file);
        out->out_close (out);
        return -1;
     }

   out->formatted_file_exists = 1;

   return 0;
}

static int create_irr_file (Output_Type *out)
{
   return create_file_of_type (out, TIO_l1_irradiance_template);
}

static int create_rad_file (Output_Type *out)
{
   return create_file_of_type (out, TIO_l1_radiance_template);
}

static void out_free (Output_Type *out)
{
   if (out == NULL)
     return;
   image_free (out->img_outbuf);
   FREE(out->file);
   FREE(out);
}

static int out_file_exists (const Output_Type *out)
{
   return out->formatted_file_exists ? 1 : 0;
}

static void image_pixels_to_outbuf (const Image_Type *img,
                                    int ybeg, int yend,
                                    Image_Type *outbuf)
{
   Image_Pixel_Type *pixels_img = img->pixels;
   int x, nx = img->num_cols;
   int y, ny = yend - ybeg;

   /* Copy wavelength range [ybeg,yend) from img -> outbuf.
    * In img, x varies fastest. In outbuf, y varies fastest. */

   for (x = 0; x < nx; x++)
     {
        Image_Pixel_Type *pixels_out = outbuf->pixels + x * ny;
        for (y = ybeg; y < yend; y++)
          {
             pixels_out[y-ybeg] = pixels_img[x + y * nx];
          }
     }
}

static void image_pqf_to_outbuf (const Image_Type *img,
                                 int ybeg, int yend,
                                 Image_Type *outbuf)
{
   Image_Pqf_Bitmap_Type *pqf_img = img->pixel_quality_flags;
   int x, nx = img->num_cols;
   int y, ny = yend - ybeg;

   /* Copy wavelength range [ybeg,yend) from img -> outbuf.
    * In img, x varies fastest. In outbuf, y varies fastest. */

   for (x = 0; x < nx; x++)
     {
        Image_Pqf_Bitmap_Type *pqf_out = outbuf->pixel_quality_flags + x * ny;
        for (y = ybeg; y < yend; y++)
          {
             pqf_out[y-ybeg] = pqf_img[x + y * nx];
          }
     }
}

static int write_rec_band1 (Output_Type *out, const char *name_var,
                            const char *name_var_err,
                            int band_index, int index,
                            const Output_Exprec_Type *rec)
{
   const Granule_Exprec_Type *exprec = rec->exprec;
   const Image_Type *img_err = rec->img_err;
   const Image_Type *img_waves = rec->img_waves;
   Image_Type *img_outbuf = out->img_outbuf;
   Band_Info_Type *band = &out->bands[band_index];
   int ybeg = band->ybeg;
   int yend = band->yend;
   int grp, start[3], count[3];

   if (0 != TIO_inq_grp (out->ncid, band->name, &grp))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: accessing group %s in %s",
                     __func__, band->name, out->file);
        return -1;
     }

   start[0] = index;
   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = out->num_xtrack;
   count[2] = out->num_waves;

   image_pixels_to_outbuf (exprec->img, ybeg, yend, img_outbuf);
   if (0 != TIO_put_var_section (grp, name_var, start, count, TIO_FLOAT,
                                 img_outbuf->pixels))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, name_var, out->file);
        return -1;
     }

   image_pqf_to_outbuf (exprec->img, ybeg, yend, img_outbuf);
   if (0 != TIO_put_var_section (grp, TEMPO_VAR_PQF, start, count, TIO_USHORT,
                                 img_outbuf->pixel_quality_flags))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_PQF, out->file);
        return -1;
     }

   image_pixels_to_outbuf (img_err, ybeg, yend, img_outbuf);
   if (0 != TIO_put_var_section (grp, name_var_err, start, count, TIO_FLOAT,
                                 img_outbuf->pixels))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, name_var_err, out->file);
        return -1;
     }

   image_pixels_to_outbuf (img_waves, ybeg, yend, img_outbuf);
   if (0 != TIO_put_var_section (grp, TEMPO_VAR_WAVELENGTH, start, count, TIO_FLOAT,
                                 img_outbuf->pixels))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing %s to %s",
                     __func__, TEMPO_VAR_WAVELENGTH, out->file);
        return -1;
     }

   return 0;
}

static int write_rec_bands (Output_Type *out, const char *name_var,
                            const char *name_var_err, int index,
                            const Output_Exprec_Type *rec)
{
   int b;

   if (0 == out_file_exists (out))
     {
        tell_verror (TELL_USAGE_ERROR, "%s: no open output file",
                     __func__);
        return -1;
     }

   for (b = 0; b < NUM_BANDS; b++)
     {
        if (0 != write_rec_band1 (out, name_var, name_var_err, b, index, rec))
          return -1;
     }

   return 0;
}

static int write_irr_rec (Output_Type *out, int index,
                          const Output_Exprec_Type *rec)
{
   return write_rec_bands (out, TEMPO_VAR_IRRADIANCE, TEMPO_VAR_IRRADIANCE_ERROR,
                           index, rec);
}

static int write_rad_rec (Output_Type *out, int index,
                          const Output_Exprec_Type *rec)
{
   return write_rec_bands (out, TEMPO_VAR_RADIANCE, TEMPO_VAR_RADIANCE_ERROR,
                           index, rec);
}

Output_Type *output_alloc (config_t *cfg, int exposure_type)
{
   Output_Type *out = NULL;

   (void) cfg;

   if (NULL == (out = (Output_Type *) MALLOC (sizeof *out)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)out, 0, sizeof *out);

   out->out_free = out_free;
   out->out_set_file = out_set_file;
   out->out_set_dims = out_set_dims;
   out->out_close = out_close;
   out->out_file_exists = out_file_exists;

   switch (exposure_type)
     {
      default:
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR,
                     "%s: no support for exposure type %d",
                     __func__, exposure_type);
        FREE(out);
        return NULL;

      case EXPREC_TYPE_RADIANCE:
        out->out_create = create_rad_file;
        out->out_write_rec = write_rad_rec;
        break;

      case EXPREC_TYPE_LIN_IRR:
        /* drop */
      case EXPREC_TYPE_IRRADIANCE:
        out->out_create = create_irr_file;
        out->out_write_rec = write_irr_rec;
        break;
     }

   return out;
}
