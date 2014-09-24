#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "netcdf.h"
#include "tio.h"
#include "_tio.h"

#define TIO_CHUNKSIZE_STEP 256

#define COMMENT_WGS84 "Earth-centered WGS84 Cartesian coordinates (z = North Pole, xy=equator, x = prime meridian)"

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
   _pDim_Type xy_det;            /* detector coordinates */
   _pDim_Type xy_sma;            /* scan mirror assembly (SMA) pointing direction */
   _pDim_Type xyz_sat;           /* satellite body axis coordinates */
   _pDim_Type xyz;               /* WGS84 Cartesian spatial coordinates */
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
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_XYZ,xyz),
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_XYSMA,xy_sma),
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_XYDET,xy_det),
       _pDIM_OFFSET_ENTRY(TIO_DIM_NAME_XYZSAT,xyz_sat),
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
   int status, varid;
   int dims[TIO_MAX_VAR_DIMS];

   /* coordinate variables */
   dims[0] = dim_table->step.id;
   status = nc_def_var (grp, TIO_DIM_NAME_STEP, NC_INT, 1, dims, NULL);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   dims[0] = dim_table->corner.id;
   status = nc_def_var (grp, TIO_DIM_NAME_CORNER, NC_INT, 1, dims, NULL);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   /* time */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME, NC_FLOAT, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* exptime */
     {
        static _pText_Attr_Type exptime_attrs[] =
          {
             {"units", "s"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_EXPTIME, NC_FLOAT, 1, dims, exptime_attrs, NULL))
          return -1;
     }

   /* pixel_size */
     {
        static _pText_Attr_Type pixel_size_attrs[] =
          {
             {"units", "micron"},
             {"comment", "pixel_size[i], applies to dimension i, where i=0 varies slowest"},
             _pTEXT_ATTRS_END
          };
        static float pixel_size[] = {_pTIO_PIXEL_YSIZE, _pTIO_PIXEL_XSIZE};

        dims[0] = dim_table->xy_det.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_PIXELSIZE, NC_FLOAT, 1, dims, pixel_size_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, pixel_size)))
          {
             _pTIO_err_verror_nc (status, "%s: writing pixel size", __func__);
             return -1;
          }
     }

   /* pixel_scale */
     {
        static _pText_Attr_Type pixel_scale_attrs[] =
          {
             {"units", "microradian"},
             {"comment", "pixel_scale[i], applies to dimension i, where i=0 varies slowest"},
             _pTEXT_ATTRS_END
          };
        static float pixel_scale[] = {_pTIO_PIXEL_YSCALE, _pTIO_PIXEL_XSCALE};

        dims[0] = dim_table->xy_det.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_PIXELSCALE, NC_FLOAT, 1, dims, pixel_scale_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, pixel_scale)))
          {
             _pTIO_err_verror_nc (status, "%s: writing pixel scale", __func__);
             return -1;
          }
     }

   return 0;
}

static int define_inr_status (int grp, enum TIO_INR_Status inr_status)
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
        _pTIO_err_verror_nc (status, "%s: defining inr_status attribute", __func__);
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
        {"last_granule_of_scan", 0},
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

static int define_radiance_group (int parent_grp, TIO_Radiance_Group_Type *rg,
                                  _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];
   int shuffle=1, deflate_level=2, deflate;
   size_t total_num;
#if 0
   int storage = NC_CHUNKED;
   size_t chunksizes[TIO_MAX_VAR_DIMS];
#endif

   if (rg->name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, rg->name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, rg->name);
        return -1;
     }

   /* group-local dimensions */
   dim_table->xtrack.len = rg->num_xtrack;
   dim_table->channel.len = rg->num_channels;

   if (-1 == define_radiance_group_dims (grp, dim_table))
     return -1;

   /* group-local coordinate variables */
   dims[0] = dim_table->xtrack.id;
   status = nc_def_var (grp, TIO_DIM_NAME_XTRACK, NC_INT, 1, dims, NULL);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   /* FIXME */
   total_num = (dim_table->channel.len
                * dim_table->xtrack.len
                * dim_table->step.len);
   deflate = (total_num > 1000000);

   /* radiance */
     {
        static _pText_Attr_Type radiance_attrs[] =
          {
             {"units", "TBD"},
             {"coordinates", "longitude latitude spectral_channel"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_RADIANCE, NC_FLOAT, 3, dims, radiance_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             _pTIO_err_verror_nc (status, "defining %s compression parameters", TIO_VAR_NAME_RADIANCE);
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
             _pTIO_err_verror_nc (status, "defining %s chunking parameters", TIO_VAR_NAME_RADIANCE);
             return -1;
          }
#endif
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* wavelength */
     {
        static _pText_Attr_Type wavelength_attrs[] =
          {
             {"units", "nm"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_WAVELENGTH, NC_DOUBLE, 3, dims, wavelength_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             _pTIO_err_verror_nc (status, "defining %s compression parameters", TIO_VAR_NAME_WAVELENGTH);
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
             _pTIO_err_verror_nc (status, "defining %s chunking parameters", TIO_VAR_NAME_WAVELENGTH);
             return -1;
          }
#endif
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_DOUBLE))
          return -1;
     }

   /* longitude */
     {
        static _pText_Attr_Type lon_attrs[] =
          {
             {"units", "degrees_east"},
             {"long_name", "longitude"},
             {"bounds", "longitude_bounds"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "longitude", NC_FLOAT, 2, dims, lon_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* latitude */
     {
        static _pText_Attr_Type lat_attrs[] =
          {
             {"units", "degrees_north"},
             {"long_name", "latitude"},
             {"bounds", "latitude_bounds"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "latitude", NC_FLOAT, 2, dims, lat_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* ellipsoid altitude */
     {
        static _pText_Attr_Type ell_alt_attrs[] =
          {
             {"units", "m"},
             {"long_name", "ellipsoid_altitude"},
             {"bounds", "ellipsoid_altitude_bounds"},
             {"coordinates", "longitude latitude"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "ellipsoid_altitude", NC_FLOAT, 2, dims, ell_alt_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* longitude bounds */
     {
        static _pText_Attr_Type lon_bnds_attrs[] =
          {
             {"units", "degrees_east"},
             {"long_name", "longitude bounds (SW,SE,NE,NW)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "longitude_bounds", NC_FLOAT, 3, dims, lon_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* latitude bounds */
     {
        static _pText_Attr_Type lat_bnds_attrs[] =
          {
             {"units", "degrees_north"},
             {"long_name", "latitude bounds (SW,SE,NE,NW)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "latitude_bounds", NC_FLOAT, 3, dims, lat_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* ellipsoid altitude bounds */
     {
        static _pText_Attr_Type ell_alt_bnds_attrs[] =
          {
             {"units", "m"},
             {"long_name", "ellipsoid altitude at bounds (SW,SE,NE,NW)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "ellipsoid_altitude_bounds", NC_FLOAT, 3, dims, ell_alt_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_INRQF, NC_UINT, 2, dims, inrqf_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_UINT))
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_DQF, NC_UINT, 2, dims, dqf_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_UINT))
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
   int status, grp;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
        return -1;
     }

   /* satellite position */
     {
        static _pText_Attr_Type satpos_attrs[] =
          {
             {"units", "km"},
             {"comment", "Satellite position in " COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xyz.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SATPOS, NC_FLOAT, 2, dims, satpos_attrs, NULL))
          return -1;
     }

   /* sun position */
     {
        static _pText_Attr_Type sunpos_attrs[] =
          {
             {"units", "km"},
             {"comment", "Sun position in " COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xyz.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SUNPOS, NC_FLOAT, 2, dims, sunpos_attrs, NULL))
          return -1;
     }

   /* moon position */
     {
        static _pText_Attr_Type moonpos_attrs[] =
          {
             {"units", "km"},
             {"comment", "Moon position in " COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xyz.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_MOONPOS, NC_FLOAT, 2, dims, moonpos_attrs, NULL))
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
   int status, grp;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_EPHEM, dim_table->time_ephemeris.len, &dim_table->time_ephemeris.id)))
     {
        _pTIO_err_verror_nc (status, "%s: defining dimension 'time' in group %s", __func__, grp_name);
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_EPHEM, NC_FLOAT, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* solar radiation pressure scalar */
     {
        static _pText_Attr_Type srp_attrs[] =
          {
             {"units", "microPascal"},
             _pTEXT_ATTRS_END
          };
        float solar_radiation_pressure = 9.08;  /* perfect reflectance, normal to surface */
        int varid;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SRP, NC_FLOAT, 0, NULL, srp_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, &solar_radiation_pressure)))
          {
             _pTIO_err_verror_nc (status, "%s: writing solar radiation pressure", __func__);
             return -1;
          }
     }

   /* satellite position */
     {
        static _pText_Attr_Type satpos_attrs[] =
          {
             {"units", "km"},
             {"comment", "Satellite position in " COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_ephemeris.id;
        dims[1] = dim_table->xyz.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SATPOS, NC_FLOAT, 2, dims, satpos_attrs, NULL))
          return -1;
     }

   /* satellite velocity */
     {
        static _pText_Attr_Type satvel_attrs[] =
          {
             {"units", "km/s"},
             {"comment", "Satellite velocity in " COMMENT_WGS84},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_ephemeris.id;
        dims[1] = dim_table->xyz.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SATVEL, NC_FLOAT, 2, dims, satvel_attrs, NULL))
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
   int status, grp;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_MANEUVER, dim_table->time_maneuvers.len, &dim_table->time_maneuvers.id)))
     {
        _pTIO_err_verror_nc (status, "%s: defining dimension 'time' in group %s", __func__, grp_name);
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_MANEUVER, NC_FLOAT, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* delta_v */
     {
        static _pText_Attr_Type deltav_attrs[] =
          {
             {"units", "m/s"},
             {"comment", "Velocity change in coordinates defined by satellite body axes (roll, pitch, yaw)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_maneuvers.id;
        dims[1] = dim_table->xyz_sat.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_DELTAV, NC_FLOAT, 2, dims, deltav_attrs, NULL))
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
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_GYRO, dim_table->time_gyroscope.len, &dim_table->time_gyroscope.id)))
     {
        _pTIO_err_verror_nc (status, "%s: defining dimension 'time' in group %s", __func__, grp_name);
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_GYRO, NC_FLOAT, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* gyro_raw */
     {
        static _pText_Attr_Type gyro_raw_attrs[] =
          {
             {"units", "radians/s"},
             {"comment", "Roll, pitch, yaw rate"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_gyroscope.id;
        dims[1] = dim_table->xyz_sat.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GYRORAW, NC_FLOAT, 2, dims, gyro_raw_attrs, NULL))
          return -1;
     }

   /* gyro_bias */
     {
        static _pText_Attr_Type gyro_bias_attrs[] =
          {
             {"units", "radians/s"},
             {"comment", "Roll, pitch, yaw rate bias"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->xyz_sat.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GYROBIAS, NC_FLOAT, 1, dims, gyro_bias_attrs, NULL))
          return -1;
     }

   /* gyro_scale */
     {
        static _pText_Attr_Type gyro_scale_attrs[] =
          {
             {"units", "radians/s"},
             {"comment", "Roll, pitch, yaw rate scale"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->xyz_sat.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_GYROSCALE, NC_FLOAT, 1, dims, gyro_scale_attrs, NULL))
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
   int status, grp;
   int dims[TIO_MAX_VAR_DIMS];

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
        return -1;
     }

   /* group-local dimensions */
   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_SMA, dim_table->time_sma.len, &dim_table->time_sma.id)))
     {
        _pTIO_err_verror_nc (status, "%s: defining dimension 'time' in group %s", __func__, grp_name);
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
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME_SMA, NC_FLOAT, 1, dims, time_attrs, NULL))
          return -1;
     }

   /* dit */
     {
        static _pText_Attr_Type dit_attrs[] =
          {
             {"units", "radians"},
             {"comment", "Scan mirror pointing direction (East, North)"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->time_sma.id;
        dims[1] = dim_table->xy_sma.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_SMADIT, NC_FLOAT, 2, dims, dit_attrs, NULL))
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
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
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
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "%s: defining group %s", __func__, grp_name);
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

int TIO_create_l1_template (int ncid, size_t num_steps, int num_rgrps,
                            TIO_Radiance_Group_Type *rgrps)
{
   _pDim_Table_Type dim_table;
   int i;

   /* Initialize the dimension sizes that are known at this point.
    * Other dimensions are group-specific and are initialized only
    * when those groups are being defined.
    */
   dim_table.step.len = num_steps;
   dim_table.corner.len = 4;
   dim_table.xy_det.len = 2;
   dim_table.xy_sma.len = 2;
   dim_table.xyz_sat.len = 3;
   dim_table.xyz.len = 3;
   dim_table.cov.len = 3;
   dim_table.time_ephemeris.len = NC_UNLIMITED;
   dim_table.time_maneuvers.len = NC_UNLIMITED;
   dim_table.time_gyroscope.len = NC_UNLIMITED;
   dim_table.time_sma.len = NC_UNLIMITED;

   if ((-1 == define_global_attrs (ncid))
       || (-1 == define_global_dims (ncid, &dim_table))
       || (-1 == define_global_vars (ncid, &dim_table)))
     {
        _pTIO_err_verror ("%s failed", __func__);
        return -1;
     }

   for (i = 0; i < num_rgrps; i++)
     {
        if (-1 == define_radiance_group (ncid, &rgrps[i], &dim_table, NULL))
          {
             _pTIO_err_verror ("%s failed defining radiance group %d", __func__, i);
             return -1;
          }
     }

   if ((-1 == define_geometry_group (ncid, TIO_GRP_NAME_GEOMETRY, &dim_table, NULL))
       || (-1 == define_inr_input_group (ncid, TIO_GRP_NAME_INRINPUT, &dim_table, NULL)))
     {
        _pTIO_err_verror ("%s failed", __func__);
        return -1;
     }

   return 0;
}
