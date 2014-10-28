#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "netcdf.h"
#include "terr.h"
#include "tio.h"
#include "_tio.h"

#define TIO_CHUNKSIZE_STEP 256

#define COMMENT_WGS84 \
 "Earth-centered WGS84 Cartesian coordinates (z = North Pole, xy=equator, x = prime meridian)"
#define COORDINATE_AT_EXPOSURE_START "coordinate at exposure start"

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

   _pDim_Type time_ephemeris;    /* ephemeris data point times */
   _pDim_Type time_maneuvers;    /* maneuver times */
   _pDim_Type time_gyroscope;    /* gyroscope sample times */
   _pDim_Type time_sma;          /* SMA DIT (differential impedance transducer) sample times */
};

static int define_global_dims (int grp, _pDim_Table_Type *dim_table)
{
   static _pDim_Offsets_Type dim_offsets[] =
    {
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_STEP,step),
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_CORNER,corner),
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_COV,cov),
       _pDIM_OFFSETS_END
    };

   return _pTIO_define_dims_using_offsets (grp, dim_offsets, dim_table);
}

static int define_radiance_group_dims (int grp, _pDim_Table_Type *dim_table)
{
   static _pDim_Offsets_Type dim_offsets[] =
    {
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_XTRACK,xtrack),
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_CHANNEL,channel),
       _pDIM_OFFSETS_END
    };

   return _pTIO_define_dims_using_offsets (grp, dim_offsets, dim_table);
}

static int define_global_vars (int grp, const _pDim_Table_Type *dim_table)
{
   int status, varid, dims[TIO_MAX_VAR_DIMS];

   /* coordinate variables */
   dims[0] = dim_table->step.id;
   status = nc_def_var (grp, TIO_DIM_NAME_STEP, NC_INT, 1, dims, NULL);
   if (NC_NOERR != status)
     {
        Terr_verror (TERR_IO_WRITE_ERROR,
                     "%s: defining coordinate variable %s (%s)",
                     __func__, TIO_DIM_NAME_STEP, nc_strerror(status));
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME, NC_DOUBLE, 1, dims, time_attrs, NULL))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_EXPOSURE_TIME, NC_FLOAT, 1, dims, exposure_time_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GRANULE_FLAG, NC_INT, 0, NULL, granule_flag_attrs, NULL))
          return -1;
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
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining inr_status attribute (%s)",
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
        {"time_coverage_start", "2018-09-01T12:00:00 UTC"},
        {"time_coverage_end", "2018-09-01T13:00:00 UTC"},
        _pTEXT_ATTRS_END
     };
   static _pInt_Attr_Type int_attrs[] =
     {
        {"processing_version", 0},
        {"granule_seq_num", 0},
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
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];
   int shuffle, deflate=1, deflate_level=1;
#if 0
   int storage = NC_CHUNKED;
   size_t chunksizes[TIO_MAX_VAR_DIMS];
#endif

   shuffle = deflate;

   if (sg->name == NULL)
     {
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, sg->name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
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
   status = nc_def_var (grp, TIO_DIM_NAME_XTRACK, NC_INT, 1, dims, NULL);
   if (NC_NOERR != status)
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining coordinate variable %s (%s)",
                     __func__, TIO_DIM_NAME_XTRACK, nc_strerror(status));
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

        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_PIXEL_SCALE_ROW, NC_FLOAT, 0, NULL, pixel_scale_row_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &pixel_scale_row)))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "%s: writing pixel scale (%s)",
                          __func__, nc_strerror(status));
             return -1;
          }
     }

   /* pixel_scale_column */
     {
        static _pText_Attr_Type pixel_scale_column_attrs[] =
          {
             {"units", "microradian"},
             {"comment", "Nominal angular size of one spatial image pixel."},
             _pTEXT_ATTRS_END
          };
        float pixel_scale_column = _pTIO_PIXEL_SCALE_COLUMN;

        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_PIXEL_SCALE_COLUMN, NC_FLOAT, 0, NULL, pixel_scale_column_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &pixel_scale_column)))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "%s: writing pixel scale (%s)",
                          __func__, nc_strerror(status));
             return -1;
          }
     }

   /* mirror_step_size */
     {
        static _pText_Attr_Type mirror_step_size_attrs[] =
          {
             {"units", "microradian"},
             {"comment", "Nominal size of a mirror step from one scan position to the next."},
             _pTEXT_ATTRS_END
          };
        float mirror_step_size = _pTIO_MIRROR_STEP_SIZE;

        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_MIRROR_STEP_SIZE, NC_FLOAT, 0, NULL, mirror_step_size_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &mirror_step_size)))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "%s: writing mirror step size (%s)",
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
             _pTEXT_ATTRS_END
          };
        float radiance_fill = TIO_FILL_RADIANCE;

        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_RADIANCE, NC_FLOAT, 3, dims, radiance_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_att (grp, varid, _FillValue, NC_FLOAT, 1, &radiance_fill)))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "writing %s fill value to grp=%d (%s)",
                          TIO_VAR_NAME_RADIANCE, grp, nc_strerror(status));
             return -1;
          }
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "defining %s compression parameters (%s)",
                          TIO_VAR_NAME_RADIANCE, nc_strerror(status));
             return -1;
          }
#if 0
        /* FIXME */
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = dim_table->xtrack.len;
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (NC_NOERR != (status = nc_def_var_chunking (grp, varid, storage, chunksizes))))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "defining %s chunking parameters (%s)",
                          TIO_VAR_NAME_RADIANCE, nc_strerror(status));
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_WAVELENGTH, NC_FLOAT, 3, dims, wavelength_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, wavelength_float_attrs))
          return -1;

        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "defining %s compression parameters (%s)",
                          TIO_VAR_NAME_WAVELENGTH, nc_strerror(status));
             return -1;
          }
#if 0
        /* FIXME */
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = dim_table->xtrack.len;
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (NC_NOERR != (status = nc_def_var_chunking (grp, varid, storage, chunksizes))))
          {
             Terr_verror (TERR_IO_WRITE_ERROR, "defining %s chunking parameters (%s)",
                          TIO_VAR_NAME_WAVELENGTH, nc_strerror(status));
             return -1;
          }
#endif
     }

   /* longitude */
     {
        static _pText_Attr_Type lon_text_attrs[] =
          {
             {"units", "degrees_east"},
             {"long_name", "longitude"},
             {"comment", "Longitude at pixel center"},
             {"bounds", "longitude_bounds"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "longitude", NC_FLOAT, 2, dims, lon_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, lon_float_attrs))
          return -1;
     }

   /* latitude */
     {
        static _pText_Attr_Type lat_text_attrs[] =
          {
             {"units", "degrees_north"},
             {"long_name", "latitude"},
             {"comment", "Latitude at pixel center"},
             {"bounds", "latitude_bounds"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "latitude", NC_FLOAT, 2, dims, lat_text_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, lat_float_attrs))
          return -1;
     }

   /* ellipsoid altitude */
     {
        static _pText_Attr_Type ell_alt_attrs[] =
          {
             {"units", "m"},
             {"long_name", "ellipsoid_altitude"},
             {"comment", "Ellipsoid altitude at pixel center"},
             {"bounds", "ellipsoid_altitude_bounds"},
             {"coordinates", "longitude latitude"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "ellipsoid_altitude", NC_FLOAT, 2, dims, ell_alt_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, ell_alt_float_attrs))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "longitude_bounds", NC_FLOAT, 3, dims, lon_bnds_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "latitude_bounds", NC_FLOAT, 3, dims, lat_bnds_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "ellipsoid_altitude_bounds", NC_FLOAT, 3, dims, ell_alt_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, ell_alt_float_attrs))
          return -1;
     }

   /* inr flags */
     {
        static _pText_Attr_Type inrqf_attrs[] =
          {
             {"comment", "INR quality flag"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_INRQF, NC_INT, 2, dims, inrqf_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_COVARIANCE, NC_FLOAT, 3, dims, cov_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* data quality flag */
     {
        static _pText_Attr_Type dqf_attrs[] =
          {
             {"comment", "Data quality flag"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_DQF, NC_INT, 2, dims, dqf_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_INT))
          return -1;
     }

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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_CLOUDTOPHEIGHT, NC_FLOAT, 2, dims, cloud_top_height_attrs, &varid))
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
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_X, NC_DOUBLE, 1, dims, satpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_Y, NC_DOUBLE, 1, dims, satpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_Z, NC_DOUBLE, 1, dims, satpos_z_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SUN_X, NC_DOUBLE, 1, dims, sunpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SUN_Y, NC_DOUBLE, 1, dims, sunpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SUN_Z, NC_DOUBLE, 1, dims, sunpos_z_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_MOON_X, NC_DOUBLE, 1, dims, moonpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_MOON_Y, NC_DOUBLE, 1, dims, moonpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_MOON_Z, NC_DOUBLE, 1, dims, moonpos_z_attrs, &varid))
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
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_EPHEM, dim_table->time_ephemeris.len, &dim_table->time_ephemeris.id)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR,
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_EPHEM, NC_DOUBLE, 1, dims, time_attrs, NULL))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_AMR, NC_FLOAT, 0, NULL, amr_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_X, NC_DOUBLE, 1, dims, satpos_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_Y, NC_DOUBLE, 1, dims, satpos_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_Z, NC_DOUBLE, 1, dims, satpos_z_attrs, &varid))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_VX, NC_DOUBLE, 1, dims, satvel_vx_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_VY, NC_DOUBLE, 1, dims, satvel_vy_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SAT_VZ, NC_DOUBLE, 1, dims, satvel_vz_attrs, &varid))
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
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_MANEUVER, dim_table->time_maneuvers.len, &dim_table->time_maneuvers.id)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR,
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_MANEUVER, NC_DOUBLE, 1, dims, time_attrs, NULL))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_DELTAV_X, NC_FLOAT, 1, dims, deltav_x_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_DELTAV_Y, NC_FLOAT, 1, dims, deltav_y_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_DELTAV_Z, NC_FLOAT, 1, dims, deltav_z_attrs, &varid))
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
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_GYRO, dim_table->time_gyroscope.len, &dim_table->time_gyroscope.id)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR,
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
        dims[0] = dim_table->time_gyroscope.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_GYRO, NC_DOUBLE, 1, dims, time_attrs, NULL))
          return -1;
     }

     {
        static _pText_Attr_Type roll_attrs[] =
          {
             {"units", "radians/s"},
             {"comment", "Roll rate"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type pitch_attrs[] =
          {
             {"units", "radians/s"},
             {"comment", "Pitch rate"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type yaw_attrs[] =
          {
             {"units", "radians/s"},
             {"comment", "Yaw rate"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type gyro_attr[] =
          {
             {"bias", 0.0},
             {"scale", 1.0},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->time_gyroscope.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GYRO_ROLL, NC_FLOAT, 1, dims, roll_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, gyro_attr))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GYRO_PITCH, NC_FLOAT, 1, dims, pitch_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, gyro_attr))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GYRO_YAW, NC_FLOAT, 1, dims, yaw_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, gyro_attr))
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
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_SMA, dim_table->time_sma.len, &dim_table->time_sma.id)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR,
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_SMA, NC_DOUBLE, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* dit */
     {
        static _pText_Attr_Type dit_east_attrs[] =
          {
             {"units", "radians"},
             {"comment", "Eastward angular coordinate of scan mirror pointing direction"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type dit_north_attrs[] =
          {
             {"units", "radians"},
             {"comment", "Northward angular coordinate of scan mirror pointing direction"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_sma.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SMADIT_EAST, NC_FLOAT, 1, dims, dit_east_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SMADIT_NORTH, NC_FLOAT, 1, dims, dit_north_attrs, &varid))
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
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   if ((-1 == define_gyroscope_group (grp, TIO_GRP_NAME_GYROSCOPE, dim_table, NULL))
       || (-1 == define_mirror_group (grp, TIO_GRP_NAME_MIRROR, dim_table, NULL)))
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
        Terr_verror (TERR_INVALID_PARM, "%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        Terr_verror (TERR_IO_WRITE_ERROR, "%s: defining group %s (%s)",
                     __func__, grp_name, nc_strerror(status));
        return -1;
     }

   if ((-1 == define_ephemeris_group (grp, TIO_GRP_NAME_EPHEMERIS, dim_table, NULL))
       || (-1 == define_maneuvers_group (grp, TIO_GRP_NAME_MANEUVERS, dim_table, NULL))
       || (-1 == define_telemetry_group (grp, TIO_GRP_NAME_TELEMETRY, dim_table, NULL)))
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

   /* Initialize the dimension sizes that are known at this point.
    * Other dimensions are group-specific and are initialized only
    * when those groups are being defined.
    */
   dim_table.step.len = num_steps;
   dim_table.corner.len = 4;
   dim_table.cov.len = 3;
   dim_table.time_ephemeris.len = NC_UNLIMITED;
   dim_table.time_maneuvers.len = NC_UNLIMITED;
   dim_table.time_gyroscope.len = NC_UNLIMITED;
   dim_table.time_sma.len = NC_UNLIMITED;

   if ((-1 == define_global_attrs (ncid))
       || (-1 == define_global_dims (ncid, &dim_table))
       || (-1 == define_global_vars (ncid, &dim_table)))
     {
        Terr_verror (TERR_UNKNOWN_ERROR, "%s failed", __func__);
        return -1;
     }

   for (i = 0; i < num_sgrps; i++)
     {
        if (-1 == define_radiance_group (ncid, &sgrps[i], &dim_table, NULL))
          {
             Terr_verror (TERR_IO_WRITE_ERROR,
                          "%s failed defining radiance group %d", __func__, i);
             return -1;
          }
     }

   if ((-1 == define_geometry_group (ncid, TIO_GRP_NAME_GEOMETRY, &dim_table, NULL))
       || (-1 == define_inr_input_group (ncid, TIO_GRP_NAME_INRINPUT, &dim_table, NULL)))
     {
        Terr_verror (TERR_UNKNOWN_ERROR, "%s failed", __func__);
        return -1;
     }

   return 0;
}
