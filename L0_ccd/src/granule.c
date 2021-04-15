#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libconfig.h>
#include <netcdf.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "image.h"

#define GRANULE_TYPE_PRIVATE_DATA \
   int *pixel_buffer; \
   int exposure_type; \
   int ncid; \
   int num_exprecs; \
   int num_rows; \
   int num_cols; \
   double tstart; \
   double tend;
#include "granule.h"

static void granule_free_exprec (Granule_Exprec_Type *exprec)
{
   if (exprec == NULL)
     return;
   image_free (exprec->img);
   FREE(exprec);
}

static Granule_Exprec_Type *new_exprec (int num_rows, int num_cols)
{
   Granule_Exprec_Type *exprec;

   if (NULL == (exprec = (Granule_Exprec_Type *) MALLOC (sizeof *exprec)))
     return NULL;

   memset ((char *)exprec, 0, sizeof *exprec);

   if (NULL == (exprec->img = image_new (num_rows, num_cols)))
     {
        granule_free_exprec (exprec);
        return NULL;
     }

   return exprec;
}

static int unpack_pixel_buffer (int *pixel_buffer, Image_Type *img)
{
   Image_Pixel_Type *pixels = img->pixels;
   Image_Pqf_Bitmap_Type *pixel_quality_flags = img->pixel_quality_flags;
   int i, num_pixels = img->num_rows * img->num_cols;
   int num_invalid_quality_byte = 0;

#define PIXEL_QUALITY_BYTE(v)  (((v) >> 24) & 0xff)
#define PIXEL_VALUE(v)         ((v) & ~(0xff << 24))

   for (i = 0; i < num_pixels; i++)
     {
        int bits = pixel_buffer[i];
        Image_Pixel_Type pixel_value;
        Image_Pqf_Bitmap_Type quality_flag;
        unsigned char quality_byte = PIXEL_QUALITY_BYTE(bits);

        switch (quality_byte)
          {
           case 0x00:
             /* Pixel value is good */
             quality_flag = 0;
             pixel_value = PIXEL_VALUE(bits);
             break;

           case 0x81:
             /* Packet containing the pixel has a bad Cyclic Redundancy Check (CRC) */
             quality_flag = IMAGE_PQF_BAD_PIXEL;
             pixel_value = PIXEL_VALUE(bits);
             break;

           case 0xff:
             /* Pixel value is missing (e.g. from a missing packet) */
             quality_flag = IMAGE_PQF_MISSING_DATA;
             pixel_value = IMAGE_PIXEL_FILL_VALUE;
             break;

           default:
             /* Quality flag is invalid, but we'll preserve the bits just in case */
             num_invalid_quality_byte += 1;
             quality_flag = IMAGE_PQF_BAD_PIXEL;
             pixel_value = PIXEL_VALUE(bits);
             break;
          }

        pixel_quality_flags[i] = quality_flag;
        pixels[i] = pixel_value;
     }

   if (num_invalid_quality_byte)
     {
        tell_vwarn (0, "%s: %d pixels had invalid quality byte values (flagged as bad pixels)",
                    __func__, num_invalid_quality_byte);
     }

   return 0;
}

/* Caller should allocate storage space before the call */
static int granule_get_exposure_per_frame (const Granule_Type *g, double *exposure_per_frame)
{
   double nan_value = nan("");
   int i;

   if (exposure_per_frame == NULL)
     return -1;

   for (i = 0; i < g->num_exprecs; i++)
     {
        double exposure_time;
        unsigned int num_coadds;
        int count=1;

        if ((0 != TIO_get_var_section (g->ncid, "exposure_time", &i, &count, TIO_DOUBLE, &exposure_time))
            || (0 != TIO_get_var_section (g->ncid, "num_coadds", &i, &count, TIO_UINT, &num_coadds)))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading exposure_time, num_coadds (i=%d)", __func__, i);
             return -1;
          }
        if ((num_coadds > 0) && (exposure_time != TIO_FILL_DOUBLE))
          exposure_per_frame[i] = exposure_time / num_coadds;
        else
          exposure_per_frame[i] = nan_value;
     }

   return 0;
}

static int read_ccd_int_type_enum (int ncid, int *start, int *count, int *ccd_int_type)
{
   int varid_ccd_int_type, nc_status, enum_type;
   size_t sstart = start[0];
   size_t scount = count[0];

   if ((NC_NOERR != (nc_status = nc_inq_varid (ncid, "ccd_int_type", &varid_ccd_int_type)))
       || (NC_NOERR != (nc_status = nc_inq_vartype (ncid, varid_ccd_int_type, &enum_type)))
       || (NC_NOERR != (nc_status = nc_get_vara (ncid, varid_ccd_int_type, &sstart, &scount, ccd_int_type))))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: reading ccd_int_type: (%s)",
                     __func__, nc_strerror(nc_status));
        return -1;
     }

   return 0;
}

/* If an exposure record is provided, its contents are overwritten,
 * and a pointer to it is returned.
 * Otherwise, an exposure record is allocated and returned.
 */
static Granule_Exprec_Type *
granule_read_exprec_by_index (const Granule_Type *g, int ith,
                              Granule_Exprec_Type **pexprec)
{
   Granule_Exprec_Type *exprec = NULL;
   int allocated_exprec = 0;
   int start[3], count[3];

   if ((pexprec != NULL) && (*pexprec != NULL))
     {
        exprec = *pexprec;
     }
   else
     {
        if (NULL == (exprec = new_exprec (g->num_rows, g->num_cols)))
          return NULL;
        allocated_exprec = 1;
        if (pexprec) *pexprec = exprec;
     }

   start[0] = ith;
   start[1] = 0;
   start[2] = 0;

   count[0] = 1;
   count[1] = g->num_rows;
   count[2] = g->num_cols;

   exprec->exposure_type = g->exposure_type;

   if ((0 != TIO_get_var_section (g->ncid, "image_start_time", start, count, TIO_DOUBLE,
                                  &exprec->start_time))
       ||(0 != TIO_get_var_section (g->ncid, "exposure_time", start, count, TIO_DOUBLE,
                                    &exprec->exposure_time))
       ||(0 != TIO_get_var_section (g->ncid, "frame_transfer_time", start, count, TIO_DOUBLE,
                                    &exprec->frame_transfer_time))
       ||(0 != TIO_get_var_section (g->ncid, "readout_time", start, count, TIO_DOUBLE,
                                    &exprec->readout_time))
       ||(0 != TIO_get_var_section (g->ncid, "num_coadds", start, count, TIO_UINT,
                                    &exprec->num_coadds))
       ||(0 != TIO_get_var_section (g->ncid, "curr_mirror_step", start, count, TIO_UINT,
                                    &exprec->curr_mirror_step))
       ||(0 != TIO_get_var_section (g->ncid, "num_dg_rows", start, count, TIO_UINT,
                                    &exprec->num_dg_rows))
       ||(0 != TIO_get_var_section (g->ncid, "num_tg_rows", start, count, TIO_UINT,
                                    &exprec->num_tg_rows))
       ||(0 != read_ccd_int_type_enum (g->ncid, start, count, &exprec->ccd_int_type))
      )
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: reading exposure record %d",
                     __func__, ith);
        goto error_return;
     }

   if (0 != TIO_get_var_section (g->ncid, "image", start, count, TIO_INT,
                                 g->pixel_buffer))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: reading exposure record pixels %d",
                     __func__, ith);
        goto error_return;
     }

   if (-1 == unpack_pixel_buffer (g->pixel_buffer, exprec->img))
     goto error_return;

   image_set_type (exprec->img, IMAGE_TYPE_PADDED);

   return exprec;
error_return:
   if (allocated_exprec) granule_free_exprec (exprec);
   return NULL;
}

static void granule_close (Granule_Type *g)
{
   if (g == NULL)
     return;
   if (g->ncid) TIO_close (g->ncid);
   FREE(g->pixel_buffer);
   FREE(g);
}

static int granule_num_exprecs (const Granule_Type *g)
{
   if (g == NULL)
     return -1;
   return g->num_exprecs;
}

static int granule_type (const Granule_Type *g, int *exposure_type)
{
   if (g == NULL)
     return -1;
   if (exposure_type) *exposure_type = g->exposure_type;
   return 0;
}

static int granule_get_ncid (const Granule_Type *g)
{
   return g->ncid;
}

static double granule_tstart (const Granule_Type *g)
{
   return g->tstart;
}
static double granule_tend (const Granule_Type *g)
{
   return g->tend;
}

static Granule_Type *new_granule (void)
{
   Granule_Type *g;

   if (NULL == (g = (Granule_Type *)MALLOC (sizeof *g)))
     return NULL;
   memset ((char *)g, 0, sizeof (*g));

   g->granule_close = granule_close;
   g->granule_num_exprecs = granule_num_exprecs;
   g->granule_read_exprec_by_index = granule_read_exprec_by_index;
   g->granule_get_exposure_per_frame = granule_get_exposure_per_frame;
   g->granule_free_exprec = granule_free_exprec;
   g->granule_type = granule_type;
   g->granule_ncid = granule_get_ncid;
   g->granule_tstart = granule_tstart;
   g->granule_tend = granule_tend;

   return g;
}

typedef struct
{
   char *type;
   int id;
}
Exprec_Type;
static Exprec_Type Exprec_Type_Table[] =
{
   {TEMPO_PROD_TYPESTR_DRK, EXPREC_TYPE_DARK},
   {TEMPO_PROD_TYPESTR_RAD, EXPREC_TYPE_RAD},
   {TEMPO_PROD_TYPESTR_IRR,     EXPREC_TYPE_IRR_WRK},
   {TEMPO_PROD_TYPESTR_IRR_REF, EXPREC_TYPE_IRR_REF},
   {TEMPO_PROD_TYPESTR_DRK_LIN, EXPREC_TYPE_LIN_DARK},
   {TEMPO_PROD_TYPESTR_IRR_LIN, EXPREC_TYPE_LIN_IRR},
   {NULL, EXPREC_TYPE_UNKNOWN}
};

static int identify_exprec_type (const char *exprec_type)
{
   Exprec_Type *t;

   for (t = Exprec_Type_Table; t->type != NULL; t++)
     {
        if (0 == strcmp (t->type, exprec_type))
          return t->id;
     }

   return EXPREC_TYPE_UNKNOWN;
}

/* Level 0 granule file may be dimensioned to hold a
 * fixed number of exposures, but some slots may be empty
 */
static int count_exprecs (Granule_Type *g, int time_dim)
{
   double *start_times;
   int start, count, i, num_exprecs;

   if (NULL == (start_times = (double *) MALLOC (time_dim * sizeof(double))))
     return -1;

   start = 0;
   count = time_dim;
   if (0 != TIO_get_var_section (g->ncid, "image_start_time", &start, &count, TIO_DOUBLE, start_times))
     {
        FREE(start_times);
        return -1;
     }

   num_exprecs = 0;
   for (i = 0; i < count; i++)
     {
        if (start_times[i] != TIO_FILL_DOUBLE)
          num_exprecs++;
     }

   FREE(start_times);

   return num_exprecs;
}

static int get_granule_dims (Granule_Type *g)
{
   TIO_Var_Info_Type info;

   if (0 != TIO_inq_var (g->ncid, "image", &info))
     return -1;

   if (info.ndims != 3)
     {
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: image array has %d dimensions, expected 3",
                     __func__, info.ndims);
        return -1;
     }

   if ((g->num_exprecs = count_exprecs (g, info.dimlens[0])) < 0)
     return -1;

   g->num_rows = info.dimlens[1];
   g->num_cols = info.dimlens[2];

   return 0;
}

Granule_Type *granule_open (const char *file)
{
   Granule_Type *g;
   char exprec_type[TIO_MAX_SHORT_NAME_LEN];
   size_t img_size;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: starting", __func__);

   if (NULL == (g = new_granule ()))
     return NULL;

   if (0 != TIO_open (file, NC_NOWRITE, &g->ncid))
     goto error_return;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != tio_use_file_epoch (g->ncid))
     goto error_return;

   if (0 != get_granule_dims (g))
     goto error_return;

   img_size = g->num_rows * g->num_cols * sizeof(int);
   if (NULL == (g->pixel_buffer = (int *) MALLOC (img_size)))
     goto error_return;

   memset ((char *)exprec_type, 0, sizeof(exprec_type));
   if (0 != TIO_get_att (g->ncid, NC_GLOBAL, "exprec_type", NC_CHAR, exprec_type))
     goto error_return;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "exprec_type: %s", exprec_type);
   g->exposure_type = identify_exprec_type (exprec_type);

   if ((-1 == TIO_get_att (g->ncid, NC_GLOBAL, "time_coverage_start_since_epoch", NC_DOUBLE, &g->tstart))
       || (-1 == TIO_get_att (g->ncid, NC_GLOBAL, "time_coverage_end_since_epoch", NC_DOUBLE, &g->tend)))
     {
        goto error_return;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "%s: succeeded", __func__);

   return g;

error_return:
   granule_close (g);
   return NULL;
}

