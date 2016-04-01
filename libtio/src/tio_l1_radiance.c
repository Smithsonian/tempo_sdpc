/** @file
 *  @brief TEMPO Level 1 radiance file template generation
 */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include <netcdf.h>
#include <tell.h>

#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

#define TIO_CHUNKSIZE_STEP 1
#define DO_CHUNKING        1

#define COMMENT_WGS84 \
 "Earth-centered WGS84 Cartesian coordinates (z = North Pole, xy=equator, x = prime meridian)"
#define COORDINATE_AT_EXPOSURE_START "coordinate at exposure start"

#define RADIANCE_RELERR_LOG10_MIN  (-4.0)
#define RADIANCE_RELERR_LOG10_MAX  (+2.0)

#define RADIANCE_RELERR_COMMENT \
   "Relative uncertainty, stored as a 16-bit signed integer " \
     "defined as s = 32767*[1 - 2 * (b - log10(f)) /(b-a)], " \
     "where s=short, a<b, f=float, such that 10^a <= f <= 10^b. " \
     "The attributes " RELERR_MIN_LOG10 " and " RELERR_MAX_LOG10 \
     " give the values of a and b. Use attribute " RELERR_MISSING \
     " as the value of elements stored as (short) _FillValue."

/* An instance of a _pDim_Table_Type struct is used as a lookup table
 * for all the dimensions that are defined anywhere in the associated
 * netCDF file.
 */
struct _pDim_Table_Type
{
   _pDim_Type channel;           /* dispersion direction */
   _pDim_Type xtrack;            /* pixel north-south spatial coordinate */
   _pDim_Type step;              /* mirror step position */
   _pDim_Type corner;            /* pixel corner indices */
   _pDim_Type cov;               /* unique elements of a 2x2 symmetric matrix */
   _pDim_Type gyro_axis;         /* will use either 3 or 4 depending on Host */
   _pDim_Type bias_axis;         /* 3 axes */

   _pDim_Type time_ephemeris;    /* ephemeris data point times */
   _pDim_Type time_maneuvers;    /* maneuver times */
   _pDim_Type time_gyroscope;    /* gyroscope sample times */
   _pDim_Type time_bias;         /* gyroscope bias sample times */
   _pDim_Type time_sma;          /* SMA DIT (differential impedance transducer) sample times */
};

static int define_global_dims (int grp, _pDim_Table_Type *dim_table)
{
   static _pDim_Offsets_Type dim_offsets[] =
    {
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_STEP,step),
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_CORNER,corner),
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_COV,cov),
       _pDIM_OFFSETS_END
    };

   return _pTIO_define_dims_using_offsets (grp, dim_offsets, dim_table);
}

static int define_radiance_group_dims (int grp, _pDim_Table_Type *dim_table)
{
   static _pDim_Offsets_Type dim_offsets[] =
    {
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_XTRACK,xtrack),
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_CHANNEL,channel),
       _pDIM_OFFSETS_END
    };

   return _pTIO_define_dims_using_offsets (grp, dim_offsets, dim_table);
}

static int define_global_vars (int grp, const _pDim_Table_Type *dim_table)
{
   int status, varid, dims[TIO_MAX_VAR_DIMS];

   /* coordinate variables */
   dims[0] = dim_table->step.id;
   status = nc_def_var (grp, TEMPO_DIM_STEP, NC_INT, 1, dims, NULL);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining coordinate variable %s (%s)",
                     __func__, TEMPO_DIM_STEP, nc_strerror(status));
        return -1;
     }

   /* time */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             {"comment", "Exposure start time"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME, NC_DOUBLE, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* exposure_time */
     {
        static _pText_Attr_Type exposure_time_attrs[] =
          {
             {"units", "s"},
             {"comment", "Exposure duration"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type exposure_time_float_attrs[] =
          {
             {"valid_min",  0.0},
             {"valid_max", 10.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_EXPOSURE_TIME, NC_FLOAT, 1, dims, exposure_time_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, exposure_time_float_attrs))
          return -1;
     }

   /* granule_flag */
     {
        static _pText_Attr_Type granule_flag_attrs[] =
          {
             {"flag_masks", "0x01, 0x02"},
             {"flag_meanings", "is_first_granule_of_scan, is_last_granule_of_scan"},
             _pTEXT_ATTRS_END
          };
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_GRANULE_FLAG, NC_INT, 0, NULL, granule_flag_attrs, NULL))
          return -1;
     }

   /* earth_sun_distance */
     {
        static _pText_Attr_Type earth_sun_distance_attrs[] =
          {
             {"units", "m"},
             {"comment", "Earth sun distance"},
             _pTEXT_ATTRS_END
          };
        float earth_sun_distance = _pTIO_EARTH_SUN_DISTANCE;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_EARTH_SUN_DISTANCE, NC_FLOAT, 0, NULL, earth_sun_distance_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &earth_sun_distance)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing earth sun distance (%s)",
                          __func__, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

static int define_inr_status (int grp, int inr_status)
{
   int status, enum_typeid;
   static _pEnum_Type enum_table[] =
     {
        {"none", TIO_INR_NONE},
        {"initial", TIO_INR_INITIAL},
        {"final", TIO_INR_FINAL},
        _pENUM_TABLE_END
     };

   if (-1 == _pTIO_define_enum (grp, "inr_status_enum", enum_table, &enum_typeid))
     return -1;
   status = nc_put_att (grp, NC_GLOBAL, "inr_status", enum_typeid, 1,
                        &inr_status);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining inr_status attribute (%s)",
                     __func__, nc_strerror(status));
        return -1;
     }

   return 0;
}

static int define_global_attrs (int grp)
{
   static _pText_Attr_Type text_attrs[] =
     {
        {"Conventions", TIO_CF_CONVENTION_VERSION},
        {"format_version", TIO_L1_FORMAT_VERSION},
        {"time_reference", TIO_TIME_REFERENCE_STRING},
        {"time_coverage_start", _pTIO_TIME_COVERAGE_START},
        {"time_coverage_end", _pTIO_TIME_COVERAGE_END},
        _pTEXT_ATTRS_END
     };
   static _pInt_Attr_Type int_attrs[] =
     {
        MAKE_INT_ATTR1("processing_version", 0),
        MAKE_INT_ATTR1("granule_seq_num", 0),
        MAKE_INT_ATTR1("scan_seq_num", 0),
        MAKE_INT_ATTR1("granule_num", 0),
        _pINT_ATTRS_END
     };

   if (-1 == _pTIO_define_text_attrs (grp, NC_GLOBAL, text_attrs))
     return -1;

   if (-1 == _pTIO_define_int_attrs (grp, NC_GLOBAL, int_attrs))
     return -1;

   if (-1 == _pTIO_define_processing_level (grp, TIO_PROC_LEVEL_1A))
     return -1;

   if (-1 == define_inr_status (grp, TIO_INR_NONE))
     return -1;

   return 0;
}

static _pName_Int_Pair_Type PQF_Pairs[] =
{
   {1<< 0, "missing_data"},
   {1<< 1, "bad_pixel"},
   {1<< 2, "processing_error"},
   {1<< 3, "transient_pixel"},
   {1<< 4, "rts_pixel"},
   {1<< 5, "saturated"},
   {1<< 6, "noise_underflow"},
   {1<< 7, "dark_corr_error"},
   {1<< 8, "offset_corr_error"},
   {1<< 9, "smear_corr_error"},
   {1<<10, "straylight_corr_error"},
   {1<<11, "nonlinear_range"},
   _pNAME_INT_LIST_END
};

int _pEmit_Var_Pixel_Quality_Flag (int grp, _pDim_Table_Type *dim_table)
{
   static _pText_Attr_Type pqf_attrs[] =
     {
        {"comment", "Pixel quality flag"},
        _pTEXT_ATTRS_END
     };
   int status, varid, dims[3], num_masks, len;
   int *flag_masks = NULL;
   char *flag_meanings = NULL;
   int error_status = -1;

   dims[0] = dim_table->step.id;
   dims[1] = dim_table->xtrack.id;
   dims[2] = dim_table->channel.id;
   if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_PQF, NC_SHORT, 3, dims, pqf_attrs, &varid))
     return error_status;
   if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_SHORT))
     return error_status;

   if (-1 == _pTIOMake_Name_Int_Arrays (PQF_Pairs, &num_masks, &flag_meanings,
                                        &flag_masks))
     {
        Tell_verror (TELL_INTERNAL_ERROR,
                     "%s: creating pixel_quality_flag mask arrays",
                     __func__);
        return error_status;
     }

   status = nc_put_att_int (grp, varid, "flag_masks", NC_SHORT,
                            num_masks, flag_masks);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining int attribute %s (%s)",
                     __func__, "flag_masks", nc_strerror(status));
        goto cleanup_and_return;
     }

   len = strlen (flag_meanings) + 1;
   status = nc_put_att_text (grp, varid, "flag_meanings", len, flag_meanings);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining int attribute %s (%s)",
                     __func__, "flag_meanings", nc_strerror(status));
        goto cleanup_and_return;
     }

   error_status = 0;
cleanup_and_return:
   TIO_FREE(flag_meanings);
   TIO_FREE(flag_masks);

   return error_status;
}

static _pName_Int_Pair_Type GPQF_Pairs[] =
{  /* bit 0-3 = MODIS land/water mask (?) */
  { 0,     "shallow_ocean"},
  { 1,     "land"},
  { 2,     "shallow_inland_water"},
  { 3,     "shoreline"},
  { 4,     "intermittent_water"},
  { 5,     "deep_inland_water"},
  { 6,     "continental_shelf_ocean"},
  { 7,     "deep_ocean"},
  {15,     "land_water_error"},
   /* bits 4-7 = misc */
  {1<<4,   "sun_glint_possibility"},
  {1<<5,   "solar_eclipse_possibility"},
   /* bits 8-15 = 8-bit NISE SSM/I-SSMIS EASE-grid snow/ice flags */
  { 0,     "snow_free_land"},
  /* (1-100) = sea ice concentration (%) */
  {101<<8, "permanent_ice"},
  {103<<8, "snow"},
  {252<<8, "mixed_pixels_at_coastline"},
  {253<<8, "suspect_ice_value"},
  {254<<8, "corners_undefined"},
  {255<<8, "ocean"},
   /* bits 16-23 = 8-bit MODIS yearly land cover flags, MCD12Q1, IGBP Type 1 */
  {  0,       "water"},
  {  1<<16,   "evergreen_needleleaf_forest"},
  {  2<<16,   "evergreen_broadleaf_forest"},
  {  3<<16,   "deciduous_needleleaf_forest"},
  {  4<<16,   "deciduous_broadleaf_forest"},
  {  5<<16,   "mixed_forest"},
  {  6<<16,   "closed_shrublands"},
  {  7<<16,   "open_shrublands"},
  {  8<<16,   "woody_savannas"},
  {  9<<16,   "savannas"},
  { 10<<16,   "grasslands"},
  { 11<<16,   "permanent_wetlands"},
  { 12<<16,   "croplands"},
  { 13<<16,   "urban_and_built_up"},
  { 14<<16,   "cropland_natural_vegetation_mosaic"},
  { 15<<16,   "snow_and_ice"},
  { 16<<16,   "barren_or_sparsely_vegetated"},
  {254<<16,   "unclassified"},
  {255<<16,   "fill_value"},
  _pNAME_INT_LIST_END
};

static int emit_var_ground_pixel_quality_flag (int grp, _pDim_Table_Type *dim_table)
{
   static _pText_Attr_Type gpqf_attrs[] =
     {
        {"comment", "Ground pixel quality flag"},
        {"coordinates", "longitude latitude"},
        _pTEXT_ATTRS_END
     };
   int status, varid, dims[2], num_values, len;
   int *flag_values = NULL;
   char *flag_meanings = NULL;
   int error_status = -1;

   dims[0] = dim_table->step.id;
   dims[1] = dim_table->xtrack.id;
   if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_GROUND_PIXEL_QF, NC_INT, 2, dims, gpqf_attrs, &varid))
     return error_status;
   if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_INT))
     return error_status;

   if (-1 == _pTIOMake_Name_Int_Arrays (GPQF_Pairs, &num_values, &flag_meanings,
                                        &flag_values))
     {
        Tell_verror (TELL_INTERNAL_ERROR,
                     "%s: creating ground_pixel_quality_flag value arrays",
                     __func__);
        return error_status;
     }

   status = nc_put_att_int (grp, varid, "flag_values", NC_INT, num_values, flag_values);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining int attribute %s (%s)",
                     __func__, "flag_values", nc_strerror(status));
        goto cleanup_and_return;
     }

   len = strlen (flag_meanings) + 1;
   status = nc_put_att_text (grp, varid, "flag_meanings", len, flag_meanings);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining int attribute %s (%s)",
                     __func__, "flag_meanings", nc_strerror(status));
        goto cleanup_and_return;
     }

   error_status = 0;
cleanup_and_return:
   TIO_FREE(flag_meanings);
   TIO_FREE(flag_values);

   return error_status;
}

static int define_radiance_group (int parent_grp, TIO_Scan_Group_Type *sg,
                                  _pDim_Table_Type *dim_table, int *grp_id)
{
   static _pFloat_Attr_Type lon_float_attrs[] =
     {
        {"valid_min", -180.0},
        {"valid_max", +180.0},
        {_FillValue, TIO_FILL_FLOAT},
        _pFLOAT_ATTRS_END
     };
   static _pFloat_Attr_Type lat_float_attrs[] =
     {
        {"valid_min", -90.0},
        {"valid_max", +90.0},
        {_FillValue, TIO_FILL_FLOAT},
        _pFLOAT_ATTRS_END
     };
   static _pFloat_Attr_Type ell_alt_float_attrs[] =
     {
        {"valid_min", -1.0e2},
        {"valid_max", +1.0e4},
        {_FillValue, TIO_FILL_FLOAT},
        _pFLOAT_ATTRS_END
     };
   static _pFloat_Attr_Type terr_hgt_float_attrs[] =
     {
        {"valid_min", -1.0e2},
        {"valid_max", +1.0e4},
        {_FillValue, TIO_FILL_FLOAT},
        _pFLOAT_ATTRS_END
     };
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];
   int shuffle, deflate=1, deflate_level=1;
#ifdef DO_CHUNKING
   int storage = NC_CHUNKED;
   size_t chunksizes[TIO_MAX_VAR_DIMS];
#endif

   shuffle = deflate;

   if (sg->name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, sg->name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, sg->name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   dim_table->xtrack.len = sg->num_xtrack;
   dim_table->channel.len = sg->num_channels;

   if (-1 == define_radiance_group_dims (grp, dim_table))
     return -1;

   /* group-local coordinate variables */
   dims[0] = dim_table->xtrack.id;
   status = nc_def_var (grp, TEMPO_DIM_XTRACK, NC_INT, 1, dims, NULL);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining coordinate variable %s (%s)",
                     __func__, TEMPO_DIM_XTRACK, nc_strerror(status));
        return -1;
     }

   /* pixel_scale_row */
     {
        static _pText_Attr_Type pixel_scale_row_attrs[] =
          {
             {"units", "nm"},
             {"comment", "Nominal change in dispersed wavelength across one spectral pixel."},
             _pTEXT_ATTRS_END
          };
        float pixel_scale_row = _pTIO_PIXEL_SCALE_ROW;

        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_PIXEL_SCALE_ROW, NC_FLOAT, 0, NULL, pixel_scale_row_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &pixel_scale_row)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing pixel scale (%s)",
                          __func__, nc_strerror(status));
             return -1;
          }
     }

   /* pixel_scale_column */
     {
        static _pText_Attr_Type pixel_scale_column_attrs[] =
          {
             {"units", "microradians"},
             {"comment", "Nominal angular size of one spatial image pixel."},
             _pTEXT_ATTRS_END
          };
        float pixel_scale_column = _pTIO_PIXEL_SCALE_COLUMN;

        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_PIXEL_SCALE_COLUMN, NC_FLOAT, 0, NULL, pixel_scale_column_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &pixel_scale_column)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing pixel scale (%s)",
                          __func__, nc_strerror(status));
             return -1;
          }
     }

   /* mirror_step_size */
     {
        static _pText_Attr_Type mirror_step_size_attrs[] =
          {
             {"units", "microradians"},
             {"comment", "Nominal size of a mirror step from one scan position to the next."},
             _pTEXT_ATTRS_END
          };
        float mirror_step_size = _pTIO_MIRROR_STEP_SIZE;

        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_MIRROR_STEP_SIZE, NC_FLOAT, 0, NULL, mirror_step_size_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &mirror_step_size)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing mirror step size (%s)",
                          __func__, nc_strerror(status));
             return -1;
          }
     }

   /* radiance */
     {
        static _pText_Attr_Type radiance_attrs[] =
          {
             {"units", "TBD"},
             {"coordinates", "longitude latitude spectral_channel"},
             {"ancillary_variables", TEMPO_VAR_RADIANCE_ERROR},
             _pTEXT_ATTRS_END
          };
        float radiance_fill = TIO_FILL_RADIANCE;

        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_RADIANCE, NC_FLOAT, 3, dims, radiance_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_att (grp, varid, _FillValue, NC_FLOAT, 1, &radiance_fill)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "writing %s fill value to grp=%d (%s)",
                          TEMPO_VAR_RADIANCE, grp, nc_strerror(status));
             return -1;
          }
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "defining %s compression parameters (%s)",
                          TEMPO_VAR_RADIANCE, nc_strerror(status));
             return -1;
          }
#ifdef DO_CHUNKING
        /* FIXME */
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = dim_table->xtrack.len;
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (NC_NOERR != (status = nc_def_var_chunking (grp, varid, storage, chunksizes))))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "defining %s chunking parameters (%s)",
                          TEMPO_VAR_RADIANCE, nc_strerror(status));
             return -1;
          }
#endif
     }

   /* radiance error */
     {
        static _pText_Attr_Type radiance_error_attrs[] =
          {
             {"units", ""},  /* dimensionless */
             {"coordinates", "longitude latitude spectral_channel"},
             {"comment", RADIANCE_RELERR_COMMENT},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type radiance_error_float_attrs[] =
          {
             {RELERR_MIN_LOG10, TIO_FILL_FLOAT},  /* min must be first */
             {RELERR_MAX_LOG10, TIO_FILL_FLOAT},  /* max must be second */
             {RELERR_MISSING, TIO_FILL_FLOAT},    /* u_missing must be third */
             _pFLOAT_ATTRS_END
          };
        short radiance_error_fill = -32768;

        radiance_error_float_attrs[0].value = RADIANCE_RELERR_LOG10_MIN;
        radiance_error_float_attrs[1].value = RADIANCE_RELERR_LOG10_MAX;
        radiance_error_float_attrs[2].value = TIO_FILL_RADIANCE_ERROR;

        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_RADIANCE_ERROR, NC_SHORT, 3, dims, radiance_error_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, radiance_error_float_attrs))
          return -1;
        if (NC_NOERR != (status = nc_put_att (grp, varid, _FillValue, NC_SHORT, 1, &radiance_error_fill)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "writing %s fill value to grp=%d (%s)",
                          TEMPO_VAR_RADIANCE_ERROR, grp, nc_strerror(status));
             return -1;
          }
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "defining %s compression parameters (%s)",
                          TEMPO_VAR_RADIANCE_ERROR, nc_strerror(status));
             return -1;
          }
#ifdef DO_CHUNKING
        /* FIXME */
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = dim_table->xtrack.len;
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (NC_NOERR != (status = nc_def_var_chunking (grp, varid, storage, chunksizes))))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "defining %s chunking parameters (%s)",
                          TEMPO_VAR_RADIANCE_ERROR, nc_strerror(status));
             return -1;
          }
#endif
     }

   /* wavelength */
     {
        static _pText_Attr_Type wavelength_attrs[] =
          {
             {"units", "nm"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type wavelength_float_attrs[] =
          {
             {"valid_min", TIO_FILL_FLOAT},
             {"valid_max", TIO_FILL_FLOAT},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_WAVELENGTH, NC_FLOAT, 3, dims, wavelength_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, wavelength_float_attrs))
          return -1;

        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "defining %s compression parameters (%s)",
                          TEMPO_VAR_WAVELENGTH, nc_strerror(status));
             return -1;
          }
#ifdef DO_CHUNKING
        /* FIXME */
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = dim_table->xtrack.len;
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (NC_NOERR != (status = nc_def_var_chunking (grp, varid, storage, chunksizes))))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "defining %s chunking parameters (%s)",
                          TEMPO_VAR_WAVELENGTH, nc_strerror(status));
             return -1;
          }
#endif
     }

   /* longitude */
     {
        static _pText_Attr_Type lon_text_attrs[] =
          {
             {"units", "degrees_east"},
             {"long_name", TEMPO_VAR_LONGITUDE},
             {"comment", "Longitude at pixel center"},
             {"bounds", TEMPO_VAR_LONGITUDE_BOUNDS},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_LONGITUDE, NC_FLOAT, 2, dims, lon_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, lon_float_attrs))
          return -1;
     }

   /* latitude */
     {
        static _pText_Attr_Type lat_text_attrs[] =
          {
             {"units", "degrees_north"},
             {"long_name", TEMPO_VAR_LATITUDE},
             {"comment", "Latitude at pixel center"},
             {"bounds", TEMPO_VAR_LATITUDE_BOUNDS},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_LATITUDE, NC_FLOAT, 2, dims, lat_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, lat_float_attrs))
          return -1;
     }

   /* ellipsoid altitude */
     {
        static _pText_Attr_Type ell_alt_attrs[] =
          {
             {"units", "m"},
             {"long_name", TEMPO_VAR_ELL_ALTITUDE},
             {"comment", "Ellipsoid altitude at pixel center"},
             {"bounds", TEMPO_VAR_ELL_ALTITUDE_BOUNDS},
             {"coordinates", "longitude latitude"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_ELL_ALTITUDE, NC_FLOAT, 2, dims, ell_alt_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, ell_alt_float_attrs))
          return -1;
     }

   /* terrain height */
     {
        static _pText_Attr_Type terr_hgt_attrs[] =
          {
             {"units", "m"},
             {"long_name", TEMPO_VAR_TERR_HEIGHT},
             {"comment", "Terrain height at pixel center"},
             {"bounds", TEMPO_VAR_TERR_HEIGHT_BOUNDS},
             {"coordinates", "longitude latitude"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TERR_HEIGHT, NC_FLOAT, 2, dims, terr_hgt_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, terr_hgt_float_attrs))
          return -1;
     }

   /* longitude bounds */
     {
        static _pText_Attr_Type lon_bnds_attrs[] =
          {
             {"units", "degrees_east"},
             {"long_name", "longitude bounds (NE,NW,SW,SE)"},
             {"comment", "Longitude at pixel corners"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_LONGITUDE_BOUNDS, NC_FLOAT, 3, dims, lon_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, lon_float_attrs))
          return -1;
     }

   /* latitude bounds */
     {
        static _pText_Attr_Type lat_bnds_attrs[] =
          {
             {"units", "degrees_north"},
             {"long_name", "latitude bounds (NE,NW,SW,SE)"},
             {"comment", "Latitude at pixel corners"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_LATITUDE_BOUNDS, NC_FLOAT, 3, dims, lat_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, lat_float_attrs))
          return -1;
     }

   /* ellipsoid altitude bounds */
     {
        static _pText_Attr_Type ell_alt_bnds_attrs[] =
          {
             {"units", "m"},
             {"long_name", "ellipsoid altitude at bounds (NE,NW,SW,SE)"},
             {"comment", "Ellipsoid altitude at pixel corners"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_ELL_ALTITUDE_BOUNDS, NC_FLOAT, 3, dims, ell_alt_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, ell_alt_float_attrs))
          return -1;
     }

   /* terrain height bounds */
     {
        static _pText_Attr_Type terr_hgt_bnds_attrs[] =
          {
             {"units", "m"},
             {"long_name", "terrain height at bounds (NE,NW,SW,SE)"},
             {"comment", "Terrain height at pixel corners"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TERR_HEIGHT_BOUNDS, NC_FLOAT, 3, dims, terr_hgt_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, terr_hgt_float_attrs))
          return -1;
     }

   /* solar zenith angle */
     {
        static _pText_Attr_Type sza_text_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", TEMPO_VAR_SZ_ANGLE},
             {"comment", "solar zenith angle at pixel center"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type sza_float_attrs[] =
          {
             {"valid_min",   0.0},
             {"valid_max", +90.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SZ_ANGLE, NC_FLOAT, 2, dims, sza_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, sza_float_attrs))
          return -1;
     }

   /* solar azimuth angle */
     {
        static _pText_Attr_Type saa_text_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", TEMPO_VAR_SA_ANGLE},
             {"comment", "solar azimuth angle at pixel center"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type saa_float_attrs[] =
          {
             {"valid_min", -180.0},
             {"valid_max", +180.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SA_ANGLE, NC_FLOAT, 2, dims, saa_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, saa_float_attrs))
          return -1;
     }

   /* viewing zenith angle */
     {
        static _pText_Attr_Type vza_text_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", TEMPO_VAR_VZ_ANGLE},
             {"comment", "viewing zenith angle at pixel center"},
             {"bounds", TEMPO_VAR_VZ_ANGLE_BOUNDS},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type vza_float_attrs[] =
          {
             {"valid_min",   0.0},
             {"valid_max", +90.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_VZ_ANGLE, NC_FLOAT, 2, dims, vza_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, vza_float_attrs))
          return -1;
     }

   /* viewing zenith angle bounds */
     {
        static _pText_Attr_Type vza_bnds_text_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", "viewing zenith angle at bounds (NE,NW,SW,SE)"},
             {"comment", "viewing zenith angle at pixel corners"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type vza_bnds_float_attrs[] =
          {
             {"valid_min",   0.0},
             {"valid_max", +90.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_VZ_ANGLE_BOUNDS, NC_FLOAT, 3, dims, vza_bnds_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, vza_bnds_float_attrs))
          return -1;
     }

   /* viewing azimuth angle */
     {
        static _pText_Attr_Type vaa_text_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", TEMPO_VAR_VA_ANGLE},
             {"comment", "viewing azimuth angle at pixel center"},
             {"bounds", TEMPO_VAR_VA_ANGLE_BOUNDS},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type vaa_float_attrs[] =
          {
             {"valid_min", -180.0},
             {"valid_max", +180.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_VA_ANGLE, NC_FLOAT, 2, dims, vaa_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, vaa_float_attrs))
          return -1;
     }

   /* viewing azimuth angle bounds */
     {
        static _pText_Attr_Type vaa_bnds_text_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", "viewing azimuth angle at bounds (NE,NW,SW,SE)"},
             {"comment", "viewing azimuth angle at pixel corners"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type vaa_bnds_float_attrs[] =
          {
             {"valid_min", -180.0},
             {"valid_max", +180.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_VA_ANGLE_BOUNDS, NC_FLOAT, 3, dims, vaa_bnds_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, vaa_bnds_float_attrs))
          return -1;
     }

   /* inr flags */
     {
        static _pText_Attr_Type inrqf_attrs[] =
          {
             {"comment", "INR quality flag"},
             {"coordinates", "longitude latitude"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_INRQF, NC_INT, 2, dims, inrqf_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_INT))
          return -1;
     }

   /* covariance */
     {
        static _pText_Attr_Type cov_attrs[] =
          {
             {"units", "km^2"},
             {"comment", "Unique elements of 2x2 symmetric covariance matrix, cov(0,0), cov(0,1), cov(1,1)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->cov.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_INR_COVARIANCE, NC_FLOAT, 3, dims, cov_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* pixel quality flag */
   if (-1 == _pEmit_Var_Pixel_Quality_Flag (grp, dim_table))
     return -1;

   /* ground pixel quality flag */
   if (-1 == emit_var_ground_pixel_quality_flag (grp, dim_table))
     return -1;

   /* cloud top height */
     {
        static _pText_Attr_Type cloud_top_height_attrs[] =
          {
             {"units", "m"},
             {"coordinates", "longitude latitude"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_CLOUDTOPHEIGHT, NC_FLOAT, 2, dims, cloud_top_height_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_geometry_group (int parent_grp, const char *grp_name,
                                  const _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* satellite position */
     {
        static _pText_Attr_Type satpos_x_attrs[] =
          {
             {"units", "km"},
             {"long_name", "satellite X " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type satpos_y_attrs[] =
          {
             {"units", "km"},
             {"long_name", "satellite Y " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type satpos_z_attrs[] =
          {
             {"units", "km"},
             {"long_name", "satellite Z " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_X, NC_DOUBLE, 1, dims, satpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_Y, NC_DOUBLE, 1, dims, satpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_Z, NC_DOUBLE, 1, dims, satpos_z_attrs, &varid))
          return -1;
     }

   /* sun position */
     {
        static _pText_Attr_Type sunpos_x_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Sun X " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type sunpos_y_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Sun Y " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type sunpos_z_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Sun Z " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SUN_X, NC_DOUBLE, 1, dims, sunpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SUN_Y, NC_DOUBLE, 1, dims, sunpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SUN_Z, NC_DOUBLE, 1, dims, sunpos_z_attrs, &varid))
          return -1;
     }

   /* moon position */
     {
        static _pText_Attr_Type moonpos_x_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Moon X " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type moonpos_y_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Moon Y " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type moonpos_z_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Moon Z " COORDINATE_AT_EXPOSURE_START},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_MOON_X, NC_DOUBLE, 1, dims, moonpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_MOON_Y, NC_DOUBLE, 1, dims, moonpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_MOON_Z, NC_DOUBLE, 1, dims, moonpos_z_attrs, &varid))
          return -1;
     }

   /* Sun angle */
     {
        static _pText_Attr_Type sun_angle_attrs[] =
          {
             {"units", "radians"},
             {"long_name", "Sun cone angle from boresight"},
             {"comment", "Angle between the sun and the boresight of the instrument at the exposure start"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SUN_ANGLE, NC_DOUBLE, 1, dims, sun_angle_attrs, &varid))
          return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_ephemeris_group (int parent_grp, const char *grp_name,
                                   _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TEMPO_VAR_TIME_EPHEM, dim_table->time_ephemeris.len, &dim_table->time_ephemeris.id)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining dimension 'time' in group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_ephemeris.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME_EPHEM, NC_DOUBLE, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* Effective (cross-sectional) area to mass ratio,
    * defined as: (1 + Cr)*A/m
    *      where Cr = a reflection coefficient [dimensionless, 0 <= Cr <= 1]
    *             A = satellite cross-sectional area [m^2]
    *             m = satellite mass [kg]
    */
     {
        static _pText_Attr_Type amr_attrs[] =
          {
             {"units", "m^2/kg"},
             _pTEXT_ATTRS_END
          };
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_AMR, NC_FLOAT, 0, NULL, amr_attrs, &varid))
          return -1;
     }

   /* satellite position */
     {
        static _pText_Attr_Type satpos_x_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Satellite X coordinate"},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type satpos_y_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Satellite Y coordinate"},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type satpos_z_attrs[] =
          {
             {"units", "km"},
             {"long_name", "Satellite Z coordinate"},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_ephemeris.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_X, NC_DOUBLE, 1, dims, satpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_Y, NC_DOUBLE, 1, dims, satpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_Z, NC_DOUBLE, 1, dims, satpos_z_attrs, &varid))
          return -1;
     }

   /* satellite velocity */
     {
        static _pText_Attr_Type satvel_vx_attrs[] =
          {
             {"units", "km/s"},
             {"long_name", "Satellite X velocity"},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type satvel_vy_attrs[] =
          {
             {"units", "km/s"},
             {"long_name", "Satellite Y velocity"},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type satvel_vz_attrs[] =
          {
             {"units", "km/s"},
             {"long_name", "Satellite Z velocity"},
             {"comment", COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_ephemeris.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_VX, NC_DOUBLE, 1, dims, satvel_vx_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_VY, NC_DOUBLE, 1, dims, satvel_vy_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SAT_VZ, NC_DOUBLE, 1, dims, satvel_vz_attrs, &varid))
          return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_maneuvers_group (int parent_grp, const char *grp_name,
                                   _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TEMPO_VAR_TIME_MANEUVER, dim_table->time_maneuvers.len, &dim_table->time_maneuvers.id)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining dimension 'time' in group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_maneuvers.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME_MANEUVER, NC_DOUBLE, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* delta_v */
     {
        static _pText_Attr_Type deltav_x_attrs[] =
          {
             {"units", "m/s"},
             {"long_name", "satellite X delta-v"},
             {"comment", "Velocity change in coordinates defined by satellite body axes (roll, pitch, yaw)"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type deltav_y_attrs[] =
          {
             {"units", "m/s"},
             {"long_name", "satellite Y delta-v"},
             {"comment", "Velocity change in coordinates defined by satellite body axes (roll, pitch, yaw)"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type deltav_z_attrs[] =
          {
             {"units", "m/s"},
             {"long_name", "satellite Z delta-v"},
             {"comment", "Velocity change in coordinates defined by satellite body axes (roll, pitch, yaw)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_maneuvers.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_DELTAV_X, NC_FLOAT, 1, dims, deltav_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_DELTAV_Y, NC_FLOAT, 1, dims, deltav_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_DELTAV_Z, NC_FLOAT, 1, dims, deltav_z_attrs, &varid))
          return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_gyroscope_group (int parent_grp, const char *grp_name,
                                   _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if ((NC_NOERR != (status = nc_def_dim (grp, TEMPO_VAR_TIME_GYRO, dim_table->time_gyroscope.len, &dim_table->time_gyroscope.id)))
       ||(NC_NOERR != (status = nc_def_dim (grp, TEMPO_VAR_TIME_GYRO_BIAS, dim_table->time_bias.len, &dim_table->time_bias.id)))
       ||(NC_NOERR != (status = nc_def_dim (grp, TEMPO_DIM_GYRO_AXIS, dim_table->gyro_axis.len, &dim_table->gyro_axis.id)))
       ||(NC_NOERR != (status = nc_def_dim (grp, TEMPO_DIM_BIAS_AXIS, dim_table->bias_axis.len, &dim_table->bias_axis.id)))
       )
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining dimension 'time' in group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* gyro time coordinate */
     {
        static _pText_Attr_Type gyro_time_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_gyroscope.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME_GYRO, NC_DOUBLE, 1, dims, gyro_time_attrs, NULL))
          return -1;
     }

   /* gyro output */
     {
        static _pText_Attr_Type output_attrs[] =
          {
             {"units", "counts"},
             {"comment", "Raw gyro output in telemetry counts"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_gyroscope.id;
        dims[1] = dim_table->gyro_axis.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_GYRO_OUTPUT, NC_INT, 2, dims, output_attrs, NULL))
          return -1;
     }

   /* gyro bias time coordinate */
     {
        static _pText_Attr_Type bias_time_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_bias.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME_GYRO_BIAS, NC_DOUBLE, 1, dims, bias_time_attrs, NULL))
          return -1;
     }

   /* gyro bias */
     {
        static _pText_Attr_Type bias_attrs[] =
          {
             {"units", "microradian/s"},
             {"comment", "Gyro bias estimate resolved along spacecraft body axes"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_bias.id;
        dims[1] = dim_table->bias_axis.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_GYRO_BIAS, NC_FLOAT, 2, dims, bias_attrs, NULL))
          return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_mirror_group (int parent_grp, const char *grp_name,
                                _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TEMPO_VAR_TIME_SMA, dim_table->time_sma.len, &dim_table->time_sma.id)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining dimension 'time' in group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_sma.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME_SMA, NC_DOUBLE, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* dit processed */
     {
        static _pText_Attr_Type dit_scanx_attrs[] =
          {
             {"units", "microradians"},
             {"comment", "Filtered estimated control axis-x scan position"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type dit_scany_attrs[] =
          {
             {"units", "microradians"},
             {"comment", "Filtered estimated control axis-y scan position"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_sma.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SMADIT_SCANX, NC_FLOAT, 1, dims, dit_scanx_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SMADIT_SCANY, NC_FLOAT, 1, dims, dit_scany_attrs, &varid))
          return -1;
     }

   /* dit raw */
     {
        static _pText_Attr_Type dit_rawx_attrs[] =
          {
             {"units", "counts"},
             {"comment", "Filtered raw DIT readout for sensor axis-x DIT pair"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type dit_rawy_attrs[] =
          {
             {"units", "counts"},
             {"comment", "Filtered raw DIT readout for sensor axis-y DIT pair"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_sma.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SMADIT_RAWX, NC_FLOAT, 1, dims, dit_rawx_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SMADIT_RAWY, NC_FLOAT, 1, dims, dit_rawy_attrs, &varid))
          return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_telemetry_group (int parent_grp, const char *grp_name,
                                   _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp;

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   if ((-1 == define_gyroscope_group (grp, TEMPO_GRP_GYROSCOPE, dim_table, NULL))
       || (-1 == define_mirror_group (grp, TEMPO_GRP_MIRROR, dim_table, NULL)))
     return -1;

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

static int define_inr_input_group (int parent_grp, const char *grp_name,
                                   _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp;

   if (grp_name == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   if ((-1 == define_ephemeris_group (grp, TEMPO_GRP_EPHEMERIS, dim_table, NULL))
       || (-1 == define_maneuvers_group (grp, TEMPO_GRP_MANEUVERS, dim_table, NULL))
       || (-1 == define_telemetry_group (grp, TEMPO_GRP_TELEMETRY, dim_table, NULL)))
     {
        return -1;
     }

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

int TIO_l1_radiance_template (int ncid, size_t num_steps, int num_sgrps,
                              TIO_Scan_Group_Type *sgrps)
{
   _pDim_Table_Type dim_table;
   int i;

   memset ((char *)&dim_table, 0, sizeof (dim_table));

   /* Initialize the dimension sizes that are known at this point.
    * Other dimensions are group-specific and are initialized only
    * when those groups are being defined.
    */
   dim_table.step.len = num_steps;
   dim_table.corner.len = 4;
   dim_table.cov.len = 3;
   dim_table.gyro_axis.len = 4;
   dim_table.bias_axis.len = 3;
   dim_table.time_ephemeris.len = NC_UNLIMITED;
   dim_table.time_maneuvers.len = NC_UNLIMITED;
   dim_table.time_gyroscope.len = NC_UNLIMITED;
   dim_table.time_bias.len = NC_UNLIMITED;
   dim_table.time_sma.len = NC_UNLIMITED;

   if ((-1 == define_global_attrs (ncid))
       || (-1 == define_global_dims (ncid, &dim_table))
       || (-1 == define_global_vars (ncid, &dim_table)))
     {
        Tell_verror (TELL_UNKNOWN_ERROR, "%s failed", __func__);
        return -1;
     }

   for (i = 0; i < num_sgrps; i++)
     {
        if (-1 == define_radiance_group (ncid, &sgrps[i], &dim_table, NULL))
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s failed defining radiance group %d", __func__, i);
             return -1;
          }
     }

   if ((-1 == define_geometry_group (ncid, TEMPO_GRP_GEOMETRY, &dim_table, NULL))
       || (-1 == define_inr_input_group (ncid, TEMPO_GRP_INRINPUT, &dim_table, NULL)))
     {
        Tell_verror (TELL_UNKNOWN_ERROR, "%s failed", __func__);
        return -1;
     }

   return 0;
}
