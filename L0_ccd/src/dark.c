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

#define DARK_TABLE_PRIVATE_DATA \
   Dark_Array_Type *dark_array; \
   double *key_lo; \
   double *key_hi; \
   double max_delta_temp_per_bin; \
   double max_delta_exptime_per_bin; \
   int ordered_by;
#include "dark.h"

typedef struct Ordering_Type Ordering_Type;
struct Dark_Config_Type
{
   double max_delta_temp_per_bin;
   double max_delta_exptime_per_bin;
   const Ordering_Type *ordering;
};

typedef struct Dark_Type Dark_Type;

struct Dark_Type
{
   Image_Type *img;
   double exposure_time;
   double fp_temp;            /* focal plane temperature */
   double sdc;                /* storage region dark current */
   double key;                /* value of user-selected sort key */
};

struct Dark_Array_Type
{
   Dark_Type *darks;
   int num_darks;
   int internal_alloc;
   /* internal_alloc is non-zero when darks[*].img
    * are allocated within this module */
};

static double *new_darray (int n)
{
   double *a;
   size_t array_size = n * sizeof (*a);

   if (NULL == (a = (double *)MALLOC (array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)a, 0, array_size);

   return a;
}

static int *new_iarray (int n)
{
   int *a;
   size_t array_size = n * sizeof (*a);

   if (NULL == (a = (int *)MALLOC (array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)a, 0, array_size);

   return a;
}

static void free_dark_images (Dark_Array_Type *dark_array)
{
   int i, n;

   if (dark_array->darks == NULL)
     return;

   n = dark_array->num_darks;
   for (i = 0; i < n; i++)
     {
        Dark_Type *dt = &dark_array->darks[i];
        image_free (dt->img);
     }
}

void dark_array_free (Dark_Array_Type *dark_array)
{
   if (dark_array == NULL)
     return;

   if (dark_array->internal_alloc)
     free_dark_images (dark_array);

   FREE(dark_array->darks);
   FREE(dark_array);
}

Dark_Array_Type *dark_array_alloc (int num_darks)
{
   Dark_Array_Type *dark_array = NULL;
   size_t array_size = num_darks * sizeof(Dark_Type);

   if ((NULL == (dark_array = (Dark_Array_Type *)MALLOC (sizeof *dark_array)))
       || (NULL == (dark_array->darks = (Dark_Type *)MALLOC (array_size))))
     {
        FREE(dark_array);
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)dark_array->darks, 0, array_size);
   dark_array->num_darks = num_darks;

   /* By default, assume that the Image_Type storage
    * is managed elsewhere */
   dark_array->internal_alloc = 0;

   return dark_array;
}

int dark_array_elem_set (Dark_Array_Type *dark_array, int i,
                         Image_Type *img, double sdc, double fp_temp,
                         double exposure_time)
{
   Dark_Type *dt;
   if ((dark_array == NULL) || (img == NULL))
     return -1;
   dt = &dark_array->darks[i];
   dt->img = img;
   dt->exposure_time = exposure_time;
   dt->sdc = sdc;
   dt->fp_temp = fp_temp;
   return 0;
}

static int find_bin (const Dark_Table_Type *dtt, double key)
{
   double *key_lo = dtt->key_lo;
   double *key_hi = dtt->key_hi;
   int k, n = dtt->dark_array->num_darks;

   if (key < key_lo[0] || key_hi[n-1] < key)
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "%s: key=%g falls outside tabulated range [%g,%g)",
                   __func__, key, key_lo[0], key_hi[n-1]);
        /* FIXME? */
        if (key < key_lo[0])
          return 0;
        else return n-1;
     }

   for (k = 0; k < n; k++)
     {
        if (key_lo[k] <= key && key < key_hi[k])
          return k;
     }
   if (key == key_hi[n-1]) return n-1;  /* ugly */

   return -1;
}

static int dtt_interp_fptemp (const Dark_Table_Type *dtt,
                              double fptemp, Image_Type *img)
{
   Dark_Array_Type *dark_array = dtt->dark_array;
   Image_Pixel_Type *pixels;
   Image_Pqf_Bitmap_Type *pqf;
   int k, i, n;
   double key;

   /* table provides log(dc) vs 1/T */
   key = 1.0/fptemp;

   if (-1 == (k = find_bin (dtt, key)))
     return -1;

   if (0 != image_copy (dark_array->darks[k].img, img))
     return -1;

   /* FIXME: interpolate on dark table */

   pixels = img->pixels;
   pqf = img->pixel_quality_flags;
   n = img->num_rows * img->num_cols;
   for (i = 0; i < n; i++)
     {
        /* table provides log(dc) vs 1/T */
        if (pqf[i] == 0)
          pixels[i] = exp(pixels[i]);
     }

   return 0;
}

static void dtt_domain_fptemp (const Dark_Table_Type *dtt,
                               double *min_fptemp, double *max_ftpemp)
{
   Dark_Array_Type *dark_array = dtt->dark_array;

   /* key value is 1/T, but user interface value is T */
   if (min_fptemp) *min_fptemp = 1.0/dtt->key_hi[dark_array->num_darks-1];
   if (max_ftpemp) *max_ftpemp = 1.0/dtt->key_lo[0];
}

static int dtt_interp_default (const Dark_Table_Type *dtt,
                               double key, Image_Type *img)
{
   Dark_Array_Type *dark_array = dtt->dark_array;
   int k;

   /* table provides dc vs <key> */
   if (-1 == (k = find_bin (dtt, key)))
     return -1;

   /* FIXME: interpolate on dark table */

   return image_copy (dark_array->darks[k].img, img);
}

static void dtt_domain_default (const Dark_Table_Type *dtt,
                                double *min_key, double *max_key)
{
   Dark_Array_Type *dark_array = dtt->dark_array;

   if (min_key) *min_key = dtt->key_lo[0];
   if (max_key) *max_key = dtt->key_hi[dark_array->num_darks-1];
}

typedef int Interp_Method_Type (const Dark_Table_Type *, double , Image_Type *);
typedef void Domain_Method_Type (const Dark_Table_Type *, double *, double *);

struct Ordering_Type
{
   const char *name;
   int enum_value;
   Interp_Method_Type *interp_method;
   Domain_Method_Type *domain_method;
};
#define ORDERING_TABLE_END {NULL,-1,NULL,NULL}

static Ordering_Type Ordering_Methods[] =
{
   {"sdc",   DARK_TABLE_ORDERED_BY_SDC,  dtt_interp_default, dtt_domain_default},
   {"fptemp", DARK_TABLE_ORDERED_BY_TEMP, dtt_interp_fptemp, dtt_domain_fptemp},
   {"exptime", DARK_TABLE_ORDERED_BY_EXPTIME, dtt_interp_default, dtt_domain_default},
   ORDERING_TABLE_END
};

static const Ordering_Type *find_ordering_method (const char *name)
{
   const Ordering_Type *om;
   for (om = Ordering_Methods; om->name != NULL; om++)
     {
        if (0 == strcmp (name, om->name))
          return om;
     }

   return NULL;
}

static const Ordering_Type *find_ordering_method_by_enum (int enum_value)
{
   const Ordering_Type *om;
   for (om = Ordering_Methods; om->name != NULL; om++)
     {
        if (om->enum_value == enum_value)
          return om;
     }

   return NULL;
}

static int dtt_ordering (const Dark_Table_Type *dtt)
{
   if (dtt == NULL)
     return -1;
   return dtt->ordered_by;
}

static int dark_table_write (const Dark_Table_Type *dtt, int ncid)
{
   Dark_Array_Type *dark_array = dtt->dark_array;
   Image_Type *img = dark_array->darks[0].img;
   const Ordering_Type *ordering;
   const char varname_key[] = "key";
   const char varname_key_lo[] = "key_low";
   const char varname_key_hi[] = "key_high";
   const char varname_sdc[] = "storage_region_dark_current";
   const char varname_fptemp[] = "fp_temp";
   const char varname_exptime[] = "exposure_time";
   const char varname_pqf[] = TEMPO_VAR_PQF;
   const char *order_label, *attname_max_delta;
   const char *varname_img, *varname_img_units, *varname_img_comment;
   const char *varname_key_units, *varname_key_comment;
   const char *varname_key_lo_comment, *varname_key_hi_comment;
   int dimid_key, dimid_channel, dimid_xtrack;
   size_t dimlen_key, dimlen_channel, dimlen_xtrack;
   int dimids_img[3], start[3], count[3];
   int storage = NC_CHUNKED, deflate=1, deflate_level=1;
   int shuffle = deflate;
   size_t chunksizes[3];
   int varid_unused, varid_img, varid_pqf;
   int varid_key, varid_key_lo, varid_key_hi;
   double max_delta_per_bin;
   int k;

   if (NULL == (ordering = find_ordering_method_by_enum (dtt->ordered_by)))
     {
        tell_verror (TELL_INTERNAL_ERROR,
                     "%s: invalid table ordering = %d (this should never happen)",
                     __func__, dtt->ordered_by);
        return -1;
     }

   switch (dtt->ordered_by)
     {
      case DARK_TABLE_ORDERED_BY_SDC:
        order_label = ordering->name;
        varname_img = "dark_current";
        varname_img_units = "electrons/s";
        varname_img_comment = "corrected dark current";
        varname_key_units = "electrons/s";
        varname_key_comment = "Mean per-bin quadrant-averaged storage region dark current";
        varname_key_lo_comment = "Minimum per-bin quadrant-averaged storage region dark current";
        varname_key_hi_comment = "Maximum per-bin quadrant-averaged storage region dark current";
        attname_max_delta = "max_delta_temp_per_bin";
        max_delta_per_bin = dtt->max_delta_temp_per_bin;
        break;
      case DARK_TABLE_ORDERED_BY_TEMP:
        order_label = ordering->name;
        varname_img = "log_of_dark_current";
        varname_img_units = "";
        varname_img_comment = "natural log of corrected dark current";
        varname_key_units = "1/Kelvin";
        varname_key_comment = "Mean per-bin natural log of quadrant-averaged 1/temperature";
        varname_key_lo_comment = "Minimum per-bin quadrant-averaged 1/temperature";
        varname_key_hi_comment = "Maximum per-bin quadrant-averaged 1/temperature";
        attname_max_delta = "max_delta_temp_per_bin";
        max_delta_per_bin = dtt->max_delta_temp_per_bin;
        break;
      case DARK_TABLE_ORDERED_BY_EXPTIME:
        order_label = ordering->name;
        varname_img = "dark_current";
        varname_img_units = "electrons/s";
        varname_img_comment = "corrected dark current";
        varname_key_units = "seconds";
        varname_key_comment = "exposure time";
        varname_key_lo_comment = "minimum applicable exposure time";
        varname_key_hi_comment = "maximum applicable exposure time";
        attname_max_delta = "max_delta_exptime_per_bin";
        max_delta_per_bin = dtt->max_delta_exptime_per_bin;
        break;
      default:
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: invalid table ordering = %d (this should never happen!)",
                     __func__, dtt->ordered_by);
        return -1;
        break;
     }

   dimlen_key = dark_array->num_darks;
   dimlen_channel = img->num_rows;
   dimlen_xtrack = img->num_cols;

   if ((0 != TIO_def_dim (ncid, varname_key, dimlen_key, &dimid_key))
       ||(0 != TIO_def_dim (ncid, TEMPO_DIM_CHANNEL, dimlen_channel, &dimid_channel))
       ||(0 != TIO_def_dim (ncid, TEMPO_DIM_XTRACK, dimlen_xtrack, &dimid_xtrack)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark array dimensions", __func__);
        return -1;
     }

   if ((0 != TIO_def_var (ncid, varname_key, TIO_DOUBLE, 1, &dimid_key, &varid_key))
       || (0 != TIO_def_var (ncid, varname_sdc, TIO_DOUBLE, 1, &dimid_key, &varid_unused))
       || (0 != TIO_def_var (ncid, varname_fptemp, TIO_DOUBLE, 1, &dimid_key, &varid_unused))
       || (0 != TIO_def_var (ncid, varname_exptime, TIO_DOUBLE, 1, &dimid_key, &varid_unused))
      )
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark array grid variables", __func__);
        return -1;
     }

   if ((0 != TIO_def_var (ncid, varname_key_lo, TIO_DOUBLE, 1, &dimid_key, &varid_key_lo))
       || (0 != TIO_def_var (ncid, varname_key_hi, TIO_DOUBLE, 1, &dimid_key, &varid_key_hi)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark array bin edges", __func__);
        return -1;
     }

   /* This is still a level 0 output file, so spatial index varies fastest. */
   dimids_img[0] = dimid_key;
   dimids_img[1] = dimid_channel;
   dimids_img[2] = dimid_xtrack;

   chunksizes[0] = 1;
   chunksizes[1] = dimlen_channel;
   chunksizes[2] = dimlen_xtrack;
   if ((0 != TIO_def_var (ncid, varname_img, TIO_FLOAT, 3, dimids_img, &varid_img))
       || (0 != TIO_def_var_deflate (ncid, varid_img, shuffle, deflate, deflate_level))
       || (0 != TIO_def_var_chunking (ncid, varid_img, storage, chunksizes)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark array image variable", __func__);
        return -1;
     }

   if ((0 != TIO_def_var (ncid, varname_pqf, TIO_USHORT, 3, dimids_img, &varid_pqf))
       || (0 != TIO_def_var_deflate (ncid, varid_pqf, shuffle, deflate, deflate_level))
       || (0 != TIO_def_var_chunking (ncid, varid_pqf, storage, chunksizes)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dark array quality flag variables", __func__);
        return -1;
     }

   if ((0 != TIO_put_att (ncid, NC_GLOBAL, "order_label", TIO_CHAR,
                          strlen(order_label), order_label))
       || (0 != TIO_put_att (ncid, NC_GLOBAL, attname_max_delta, TIO_DOUBLE,
                             1, &max_delta_per_bin)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing dark table attributes", __func__);
        return -1;
     }

   if ((0 != TIO_put_att (ncid, varid_img, "units", TIO_CHAR, strlen(varname_img_units), varname_img_units))
       || (0 != TIO_put_att (ncid, varid_img, "comment", TIO_CHAR, strlen(varname_img_comment), varname_img_comment))
       || (0 != TIO_put_att (ncid, varid_key, "units", TIO_CHAR, strlen(varname_key_units), varname_key_units))
       || (0 != TIO_put_att (ncid, varid_key, "comment", TIO_CHAR, strlen(varname_key_comment), varname_key_comment))
       )
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing dark array attributes", __func__);
        return -1;
     }

   if ((0 != TIO_put_att (ncid, varid_key_lo, "units", TIO_CHAR, strlen(varname_key_units), varname_key_units))
       || (0 != TIO_put_att (ncid, varid_key_lo, "comment", TIO_CHAR, strlen(varname_key_lo_comment), varname_key_lo_comment))
       || (0 != TIO_put_att (ncid, varid_key_hi, "units", TIO_CHAR, strlen(varname_key_units), varname_key_units))
       || (0 != TIO_put_att (ncid, varid_key_hi, "comment", TIO_CHAR, strlen(varname_key_hi_comment), varname_key_hi_comment)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing dark array bin attributes", __func__);
        return -1;
     }

   start[0] = 0;
   count[0] = dimlen_key;
   if ((0 != TIO_put_var_section (ncid, varname_key_lo, start, count, TIO_DOUBLE, dtt->key_lo))
       || (0 != TIO_put_var_section (ncid, varname_key_hi, start, count, TIO_DOUBLE, dtt->key_hi)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing dark array grids", __func__);
        return -1;
     }

   start[1] = 0;
   start[2] = 0;
   count[0] = 1;
   count[1] = img->num_rows;
   count[2] = img->num_cols;

   for (k = 0; k < dark_array->num_darks; k++)
     {
        Dark_Type *dt = &dark_array->darks[k];
        if (dt->img == NULL)
          continue;
        img = dt->img;
        start[0] = k;
        if ((0 != TIO_put_var_section (ncid, varname_img, start, count, TIO_FLOAT, img->pixels))
            || (0 != TIO_put_var_section (ncid, varname_pqf, start, count, TIO_USHORT, img->pixel_quality_flags))
            || (0 != TIO_put_var_section (ncid, varname_key, start, count, TIO_DOUBLE, &dt->key))
            || (0 != TIO_put_var_section (ncid, varname_sdc, start, count, TIO_DOUBLE, &dt->sdc))
            || (0 != TIO_put_var_section (ncid, varname_fptemp, start, count, TIO_DOUBLE, &dt->fp_temp))
            || (0 != TIO_put_var_section (ncid, varname_exptime, start, count, TIO_DOUBLE, &dt->exposure_time))
           )
          {
             tell_verror (TELL_IO_WRITE_ERROR, "%s: writing dark array", __func__);
             return -1;
          }
     }

   return 0;
}

static int dtt_write (const Dark_Table_Type *dtt, const char *file)
{
   int status, ncid;
   const char *product_type = "drk";

   if (0 != TIO_create (file, NC_NETCDF4, &ncid))
     return -1;

   if (0 != TIO_put_att (ncid, NC_GLOBAL, "product_type", NC_CHAR, strlen(product_type), product_type))
     return -1;

   status = dark_table_write (dtt, ncid);

   (void) TIO_close (ncid);
   return status;
}

static void free_dark_table (Dark_Table_Type *dtt)
{
   if (dtt == NULL)
     return;
   FREE(dtt->key_lo);
   FREE(dtt->key_hi);
   dark_array_free (dtt->dark_array);
}

static int alloc_dark_table (Dark_Table_Type *dtt, int num_darks)
{
   free_dark_table (dtt);

   if (NULL == (dtt->dark_array = dark_array_alloc (num_darks)))
     return -1;

   /* We'll be allocating our own Image_Type storage */
   dtt->dark_array->internal_alloc = 1;

   if ((NULL == (dtt->key_lo = new_darray (num_darks)))
       || (NULL == (dtt->key_hi = new_darray (num_darks))))
     return -1;

   return 0;
}

static void dtt_delete (Dark_Table_Type *dtt)
{
   if (dtt == NULL)
     return;

   free_dark_table (dtt);
   FREE(dtt);
}

static Dark_Table_Type *dtt_alloc (void)
{
   Dark_Table_Type *dtt;

   if (NULL == (dtt = (Dark_Table_Type *) MALLOC (sizeof *dtt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)dtt, 0, sizeof *dtt);
   dtt->dtt_delete = dtt_delete;
   dtt->dtt_ordering = dtt_ordering;
   dtt->dtt_write = dtt_write;

   return dtt;
}

Dark_Table_Type *dark_table_read (const char *file)
{
   Dark_Table_Type *dtt = NULL;
   Dark_Array_Type *dark_array;
   const Ordering_Type *ordering;
   TIO_Var_Info_Type info;
   char order_label[32];
   const char *dark_array_variable;
   int k, num_darks, num_xtrack, num_channel;
   int start[3], count[3];
   int ncid;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;

   memset ((char *)order_label, 0, sizeof(order_label));
   if (0 != TIO_get_att (ncid, NC_GLOBAL, "order_label", NC_CHAR, order_label))
     return NULL;

   if (NULL == (ordering = find_ordering_method (order_label)))
     return NULL;

   tell_vlog (TELL_MSGTYPE_INFO, 1, "dark interpolation method: %s",
              order_label);

   if (NULL == (dtt = dtt_alloc ()))
     goto return_error;

   dtt->ordered_by = ordering->enum_value;
   dtt->dtt_interp = ordering->interp_method;
   dtt->dtt_domain = ordering->domain_method;

   switch (dtt->ordered_by)
     {
        case DARK_TABLE_ORDERED_BY_EXPTIME:
        case DARK_TABLE_ORDERED_BY_SDC:
          dark_array_variable = "dark_current";
          break;
        case DARK_TABLE_ORDERED_BY_TEMP:
          dark_array_variable = "log_of_dark_current";
          break;
     }

   /* these parameters should not be needed,
    * so don't bother reading them */
   dtt->max_delta_exptime_per_bin = DBL_MAX;
   dtt->max_delta_temp_per_bin = DBL_MAX;

   if (0 != TIO_inq_var (ncid, dark_array_variable, &info))
     goto return_error;

   if (info.ndims != 3)
     {
        tell_verror (TELL_INVALID_DATA_ERROR, "%s: Expecting 3D array %s",
                     __func__, dark_array_variable);
        goto return_error;
     }

   num_darks = info.dimlens[0];
   num_channel = info.dimlens[1];
   num_xtrack = info.dimlens[2];

   if (-1 == alloc_dark_table (dtt, num_darks))
     goto return_error;

   start[0] = 0;
   count[0] = num_darks;
   if ((0 != TIO_get_var_section (ncid, "key_low", start, count, TIO_DOUBLE, dtt->key_lo))
       || (0 != TIO_get_var_section (ncid, "key_high", start, count, TIO_DOUBLE, dtt->key_hi)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading dark array grids", __func__);
        goto return_error;
     }

   start[1] = 0;
   start[2] = 0;
   count[1] = num_channel;
   count[2] = num_xtrack;

   dark_array = dtt->dark_array;

   start[0] = 0;
   count[0] = 1;
   for (k = 0; k < num_darks; k++)
     {
        Image_Type *img;
        Dark_Type *dt = &dark_array->darks[k];
        if (NULL == (dt->img = image_new (num_channel, num_xtrack)))
          goto return_error;
        img = dt->img;
        start[0] = k;
        if ((0 != TIO_get_var_section (ncid, dark_array_variable, start, count, NC_FLOAT, img->pixels))
            || (0 != TIO_get_var_section (ncid, TEMPO_VAR_PQF, start, count, NC_USHORT, img->pixel_quality_flags))
            || (0 != TIO_get_var_section (ncid, "key", start, count, NC_DOUBLE, &dt->key))
            || (0 != TIO_get_var_section (ncid, "storage_region_dark_current", start, count, NC_DOUBLE, &dt->sdc))
            || (0 != TIO_get_var_section (ncid, "fp_temp", start, count, NC_DOUBLE, &dt->fp_temp))
            || (0 != TIO_get_var_section (ncid, "exposure_time", start, count, NC_DOUBLE, &dt->exposure_time))
           )
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading dark array", __func__);
             goto return_error;
          }
     }

   return dtt;
return_error:
   dtt_delete (dtt);
   if (ncid) TIO_close (ncid);
   return NULL;
}

typedef struct
{
   double min;
   double max;
   double avg;
   double rms;
}
Stats_Type;

static void darray_stats (const double *x, int n, Stats_Type *s)
{
   double sum, sum_sqdev, avg, rms, min, max;
   int i;

   min = x[0];
   max = min;
   sum = 0.0;
   for (i = 0; i < n; i++)
     {
        double x_i = x[i];
        sum += x_i;
        if (x_i < min) min = x_i;
        else if (x_i > max) max = x_i;
     }
   avg = sum / n;

   sum_sqdev = 0.0;
   for (i = 0; i < n; i++)
     {
        double dev = x[i] - avg;
        sum_sqdev += dev * dev;
     }
   rms = sqrt (sum_sqdev / n);

   s->min = min;
   s->max = max;
   s->avg = avg;
   s->rms = rms;
}

static void copy_fptemp (Dark_Type *darks, int num_darks,
                         double *fptemp)
{
   int i;
   for (i = 0; i < num_darks; i++)
     {
        fptemp[i] = darks[i].fp_temp;
     }
}

static void copy_sdc (Dark_Type *darks, int num_darks,
                      double *sdc)
{
   int i;
   for (i = 0; i < num_darks; i++)
     {
        sdc[i] = darks[i].sdc;
     }
}

static void copy_exptime (Dark_Type *darks, int num_darks,
                          double *exptime)
{
   int i;
   for (i = 0; i < num_darks; i++)
     {
        exptime[i] = darks[i].exposure_time;
     }
}

static void free_image_counters (short *img_count)
{
   FREE(img_count);
}

static short *new_image_counters (int num_darks, int num_rows, int num_cols)
{
   short *img_count = NULL;
   size_t img_count_array_size;

   img_count_array_size = num_darks * num_rows * num_cols * sizeof(short);
   if (NULL == (img_count = (short *) MALLOC (img_count_array_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)img_count, 0, img_count_array_size);

   return img_count;
}

static void dark_image_add (Image_Type *img, int use_log,
                            const Image_Type *add, short *img_count)
{
   Image_Pixel_Type *pixels, *add_pixels;
   Image_Pqf_Bitmap_Type *add_pqf;
   int p, num_rows = img->num_rows;
   int s, num_cols = img->num_cols;
   short *count;

   /* Assume that initially, img->pixel_quality_flag[*] = 0, and
    * since we don't update it, there's no reason to test it.
    * By assuming all pixels in img->pixels are good, we ensure
    * that all good-pixel contributions from every instance of
    * add->pixels will be included in the final sum.
    */

   for (p = 0 ; p < num_rows; p++)
     {
        pixels = img->pixels + p * num_cols;
        add_pixels = add->pixels + p * num_cols;
        add_pqf = add->pixel_quality_flags + p * num_cols;
        count = img_count + p * num_cols;
        for (s = 0; s < num_cols; s++)
          {
             if (add_pqf[s] == 0)
               {
                  Image_Pixel_Type add_pixels_s = add_pixels[s];
                  if (use_log == 0)
                    {
                       pixels[s] += add_pixels_s;
                       count[s] += 1;
                    }
                  else if (add_pixels_s > 0)
                    {
                       pixels[s] += log(add_pixels_s);
                       count[s] += 1;
                    }
               }
          }
     }
}

static void dark_image_mean (Image_Type *img, const short *img_count)
{
   Image_Pixel_Type *pixels;
   Image_Pqf_Bitmap_Type *pqf;
   int p, num_rows = img->num_rows;
   int s, num_cols = img->num_cols;
   const short *count;

   for (p = 0 ; p < num_rows; p++)
     {
        pixels = img->pixels + p * num_cols;
        count = img_count + p * num_cols;
        pqf = img->pixel_quality_flags + p * num_cols;
        for (s = 0; s < num_cols; s++)
          {
             if (count[s] > 0)
               {
                  pixels[s] /= count[s];
                  if (pixels[s] < 0) pqf[s] |= IMAGE_PQF_DARK_CORR_ERROR;
               }
             else
               {
                  pixels[s] = IMAGE_PIXEL_FILL_VALUE;
                  pqf[s] |= IMAGE_PQF_MISSING_DATA;
               }
          }
     }
}

static int dark_incr (Dark_Type *dt, int use_log,
                      const Dark_Type *add, short *img_count)
{
   if (dt->img == NULL)
     {
        if (NULL == (dt->img = image_new (add->img->num_rows, add->img->num_cols)))
          return -1;
     }

   dark_image_add (dt->img, use_log, add->img, img_count);
   dt->key += add->key;
   dt->sdc += add->sdc;
   dt->fp_temp += add->fp_temp;
   dt->exposure_time += add->exposure_time;

   return 0;
}

static void dark_divide_by_counters (Dark_Type *dt,
                                     int bin_count, const short *img_count)
{
   if (bin_count == 0)
     return;
   dark_image_mean (dt->img, img_count);
   dt->key /= bin_count;
   dt->sdc /= bin_count;
   dt->fp_temp /= bin_count;
   dt->exposure_time /= bin_count;
}

static double *Qsort_Keys;
static int *Qsort_Indices;
static int Qsort_Index_Compare (const void *a, const void *b)
{
   double va = Qsort_Keys[Qsort_Indices[*(int *)a]];
   double vb = Qsort_Keys[Qsort_Indices[*(int *)b]];
   if (va < vb) return -1;
   else if (va > vb) return +1;
   else return 0;
}

static int make_dark_table_binned (Dark_Table_Type *dtt,
                                   const Dark_Array_Type *dark_array)
{
   Dark_Type *darks = dark_array->darks;
   int num_darks = dark_array->num_darks;
   Dark_Array_Type *dtt_dark_array = NULL;
   Image_Type *img = NULL;
   int *sort_index = NULL;
   int *bin_count = NULL;
   short *img_count = NULL;
   short *img_count_k = NULL;
   size_t img_size;
   double *keys = NULL;
   double *key_lo = NULL;
   double *key_hi = NULL;
   double max_delta_per_bin, key_min, key_max, delta_key;
   Stats_Type stats;
   size_t num_keys;
   int i, k, num_bins, use_log;

   if (NULL == (keys = new_darray (num_darks)))
     return -1;

   switch (dtt->ordered_by)
     {
      case DARK_TABLE_ORDERED_BY_SDC:
      case DARK_TABLE_ORDERED_BY_TEMP:
        copy_fptemp (darks, num_darks, keys);
        darray_stats (keys, num_darks, &stats);
        max_delta_per_bin = dtt->max_delta_temp_per_bin;
        break;

      case DARK_TABLE_ORDERED_BY_EXPTIME:
        copy_exptime (darks, num_darks, keys);
        darray_stats (keys, num_darks, &stats);
        max_delta_per_bin = dtt->max_delta_exptime_per_bin;
        break;
     }

   if (stats.rms > max_delta_per_bin)
     {
        num_bins = ceil ((stats.max - stats.min) / max_delta_per_bin);
     }
   else num_bins = 1;

   if (num_bins > num_darks)
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "%s: implied number of bins (%d) exceeds number of dark frames (%d)",
                         __func__, num_bins, num_darks);
        num_bins = 1;
     }

   switch (dtt->ordered_by)
     {
      case DARK_TABLE_ORDERED_BY_TEMP:
        /* tabulate log(dc) vs 1/T, which we hope is linear */
        use_log = 1;
        copy_fptemp (darks, num_darks, keys);
        for (i = 0; i < num_darks; i++)
          {
             keys[i] = 1.0/keys[i];
          }
        break;
      case DARK_TABLE_ORDERED_BY_SDC:
        /* tabulate dc vs. sdc, which we hope is linear */
        use_log = 0;
        copy_sdc (darks, num_darks, keys);
        break;
      case DARK_TABLE_ORDERED_BY_EXPTIME:
        /* tabulate dc vs. exposure time for linearity correction */
        use_log = 0;
        copy_exptime (darks, num_darks, keys);
        break;
      default:
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: invalid table ordering = %d (this should never happen!)",
                     __func__, dtt->ordered_by);
        goto return_error;
     }

   if (NULL == (sort_index = new_iarray (num_darks)))
     goto return_error;

   for (i = 0; i < num_darks; i++)
     sort_index[i] = i;

   Qsort_Keys = keys;
   Qsort_Indices = sort_index;
   num_keys = num_darks;
   qsort (sort_index, num_keys, sizeof(int), Qsort_Index_Compare);

   if (-1 == alloc_dark_table (dtt, num_bins))
     goto return_error;

   key_lo = dtt->key_lo;
   key_hi = dtt->key_hi;

   key_min = keys[sort_index[0]];
   key_max = keys[sort_index[num_keys-1]];
   delta_key = (key_max - key_min) / num_bins;

   key_lo[0] = key_min;
   for (i = 1; i < num_bins; i++)
     {
        key_lo[i] = key_min + i * delta_key;
        key_hi[i-1] = key_lo[i];
     }
   key_hi[num_bins-1] = key_max;

   img = darks[0].img;
   if ((NULL == (img_count = new_image_counters (num_bins, img->num_rows, img->num_cols)))
       || (NULL == (bin_count = new_iarray (num_bins))))
     goto return_error;
   img_size = img->num_rows * img->num_cols;

   dtt_dark_array = dtt->dark_array;

   k = 0;
   for (i = 0; i < num_darks; i++)
     {
        int is = sort_index[i];
        Dark_Type *dt = &darks[is];
        dt->key = keys[is];
        while (dt->key > key_hi[k])
          k++;
        if (k == num_bins)
          break;
        bin_count[k] += 1;
        img_count_k = img_count + k * img_size;
        if (0 != dark_incr (&dtt_dark_array->darks[k], use_log, dt, img_count_k))
          goto return_error;
     }

   /* Note that bin_count[k]==0 is possible. */
   for (k = 0; k < num_bins; k++)
     {
        img_count_k = img_count + k * img_size;
        dark_divide_by_counters (&dtt_dark_array->darks[k],
                                 bin_count[k], img_count_k);
     }

   FREE(keys);
   FREE(sort_index);
   FREE(bin_count);
   free_image_counters (img_count);

   return 0;
return_error:
   FREE(keys);
   FREE(sort_index);
   FREE(bin_count);
   free_image_counters (img_count);
   return -1;
}

extern Dark_Table_Type *dark_table_create (const Dark_Config_Type *dcfg,
                                           const Dark_Array_Type *dark_array)
{
   Dark_Table_Type *dtt;
   const Ordering_Type *ordering;

   if (dcfg == NULL)
     return NULL;

   ordering = dcfg->ordering;

   if (NULL == (dtt = dtt_alloc()))
     return NULL;

   dtt->max_delta_temp_per_bin = dcfg->max_delta_temp_per_bin;

   dtt->ordered_by = ordering->enum_value;
   dtt->dtt_interp = ordering->interp_method;
   dtt->dtt_domain = ordering->domain_method;
   dtt->dtt_write = dtt_write;

   if (0 != make_dark_table_binned (dtt, dark_array))
     {
        dtt_delete (dtt);
        return NULL;
     }

   return dtt;
}

void dark_table_config_free (Dark_Config_Type *dcfg)
{
   if (dcfg == NULL)
     return;
   FREE(dcfg);
}

Dark_Config_Type *dark_table_config (config_t *cfg, int is_linearity)
{
   config_setting_t *s;
   const Ordering_Type *ordering = NULL;
   Dark_Config_Type *dcfg = NULL;
   double max_delta_temp_per_bin;
   double max_delta_exptime_per_bin;
   const char *table_ordering;

   if (NULL == (s = config_lookup (cfg, "dark_table_creation")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing ccd_methods in param file: %s",
                     __func__, config_error_file (cfg));
        return NULL;
     }

   if ((CONFIG_TRUE != config_setting_lookup_string (s, "ordered_by", &table_ordering))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "max_delta_temp", &max_delta_temp_per_bin))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "max_delta_exptime", &max_delta_exptime_per_bin))
      )
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading dark table parameters",
                     __func__);
        return NULL;
     }

   if (NULL == (ordering = find_ordering_method (table_ordering)))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: unsupported ordering method %s",
                     __func__, table_ordering);
        return NULL;
     }

   tell_vlog (TELL_MSGTYPE_INFO, 1, "dark interpolation method: %s",
              table_ordering);

   /* FIXME? */
   if ((is_linearity != 0)
       && (ordering->enum_value != DARK_TABLE_ORDERED_BY_EXPTIME))
     {
        tell_vlog (TELL_MSGTYPE_WARN, 0,
                   "Using dark interpolation method %s to process"
                   " data of type dark-linear (expected method 'exptime')",
                   table_ordering);
     }

   if (NULL == (dcfg = (Dark_Config_Type *) MALLOC (sizeof *dcfg)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   dcfg->max_delta_temp_per_bin = fabs(max_delta_temp_per_bin);
   dcfg->max_delta_exptime_per_bin = fabs(max_delta_exptime_per_bin);
   dcfg->ordering = ordering;

   return dcfg;
}
