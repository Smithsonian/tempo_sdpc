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

typedef struct
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

   _pDim_Type time_ephemeris;  /* ephemeris data point times */
   _pDim_Type time_maneuvers;  /* maneuver times */
   _pDim_Type time_gyroscope;  /* gyroscope sample times */
   _pDim_Type time_sma;        /* SMA DIT (differential impedance transducer) sample times */
}
Dim_Table_Type;

static int define_globals (int grp, Dim_Table_Type *dim_table)
{
   int status;
   int dims[TIO_MAX_VAR_DIMS];

   /* Define global dimensions */
   status = nc_def_dim (grp, TIO_DIM_NAME_XTRACK, dim_table->xtrack.len, &dim_table->xtrack.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_CHANNEL, dim_table->channel.len, &dim_table->channel.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_STEP, dim_table->step.len, &dim_table->step.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_CORNER, dim_table->corner.len, &dim_table->corner.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_XYZ, dim_table->xyz.len, &dim_table->xyz.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_XYSMA, dim_table->xy_sma.len, &dim_table->xy_sma.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_XYDET, dim_table->xy_det.len, &dim_table->xy_det.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_XYZSAT, dim_table->xyz_sat.len, &dim_table->xyz_sat.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   status = nc_def_dim (grp, TIO_DIM_NAME_COV, dim_table->cov.len, &dim_table->cov.id);
   if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

   /* Define global coordinate variables */
     {
        int ignore_id;

        dims[0] = dim_table->step.id;
        status = nc_def_var (grp, TIO_DIM_NAME_STEP, NC_INT, 1, dims, &ignore_id);
        if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

        dims[0] = dim_table->xtrack.id;
        status = nc_def_var (grp, TIO_DIM_NAME_XTRACK, NC_INT, 1, dims, &ignore_id);
        if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

        dims[0] = dim_table->corner.id;
        status = nc_def_var (grp, TIO_DIM_NAME_CORNER, NC_INT, 1, dims, &ignore_id);
        if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;
     }

   /* Define global variables */
     {
        /* time */
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             TEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_TIME, NC_FLOAT, 1, dims, time_attrs, NULL))
          return -1;
     }

     {
        /* exptime */
        static _pText_Attr_Type exptime_attrs[] =
          {
             {"units", "s"},
             TEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_EXPTIME, NC_FLOAT, 1, dims, exptime_attrs, NULL))
          return -1;
     }

     {
        /* pixel_size */
        static _pText_Attr_Type pixel_size_attrs[] =
          {
             {"units", "micrometer"},
             {"comment", "pixel_size[i], applies to dimension i, where i=0 varies slowest"},
             TEXT_ATTRS_END
          };
        static float pixel_size[] = {_pTIO_PIXEL_YSIZE, _pTIO_PIXEL_XSIZE};
        int varid;

        dims[0] = dim_table->xy_det.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_PIXELSIZE, NC_FLOAT, 1, dims, pixel_size_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, pixel_size)))
          {
             _pTIO_err_verror_nc (status, "%s: writing pixel size", __func__);
             return -1;
          }
     }

     {
        /* pixel_scale */
        static _pText_Attr_Type pixel_scale_attrs[] =
          {
             {"units", "microradian"},
             {"comment", "pixel_scale[i], applies to dimension i, where i=0 varies slowest"},
             TEXT_ATTRS_END
          };
        static float pixel_scale[] = {_pTIO_PIXEL_YSCALE, _pTIO_PIXEL_XSCALE};
        int varid;

        dims[0] = dim_table->xy_det.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TIO_VAR_NAME_PIXELSCALE, NC_FLOAT, 1, dims, pixel_scale_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_var_float (grp, varid, pixel_scale)))
          {
             _pTIO_err_verror_nc (status, "%s: writing pixel scale", __func__);
             return -1;
          }
     }

   /* Define global attributes */
     {
        static _pText_Attr_Type text_attrs[] =
          {
             {"Conventions", TIO_CF_CONVENTION_VERSION},
             {"time_reference", TIO_TIME_REFERENCE_STRING},
             TEXT_ATTRS_END
          };
        int zero_i = 0;

        if (-1 == _pTIO_define_text_attrs (grp, NC_GLOBAL, text_attrs))
          return -1;

        status = nc_put_att_int (grp, NC_GLOBAL, "granule_seq_num", NC_INT, 1, &zero_i);
        if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

        status = nc_put_att_int (grp, NC_GLOBAL, "processing_version", NC_INT, 1, &zero_i);
        if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;

        status = nc_put_att_int (grp, NC_GLOBAL, "last_granule_of_scan", NC_INT, 1, &zero_i);
        if (_pTIO_check_verror_nc (status, __LINE__, __FILE__)) return -1;
     }

   return 0;
}

static int define_band_group (int parent_grp, const char *grp_name,
                              Dim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];
   int shuffle=1, deflate_level=2, deflate;
   size_t total_num;
#if 0
   int storage = NC_CHUNKED;
   size_t chunksizes[TIO_MAX_VAR_DIMS];
#endif

   /* FIXME */
   total_num = (dim_table->channel.len
                * dim_table->xtrack.len
                * dim_table->step.len);
   deflate = (total_num > 1000000);

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
        return -1;
     }

   /* radiance */
     {
        static _pText_Attr_Type radiance_attrs[] =
          {
             {"units", "TBD"},
             {"coordinates", "longitude latitude spectral_channel"},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->corner.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, "ellipsoid_altitude_bounds", NC_FLOAT, 3, dims, ell_alt_bnds_attrs, &varid))
          return -1;
        if (-1 == _pTIO_put_fillvalue_attr (grp, varid, NC_FLOAT))
          return -1;
     }

   /* covariance */
     {
        static _pText_Attr_Type cov_attrs[] =
          {
             {"units", "km^2"},
             {"comment", "Unique elements of 2x2 symmetric covariance matrix, cov(0,0), cov(0,1), cov(1,1)"},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
                                  Dim_Table_Type *dim_table, int *grp_id)
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
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
        return -1;
     }

   /* satellite position */
     {
        static _pText_Attr_Type satpos_attrs[] =
          {
             {"units", "km"},
             {"comment", "Satellite position in " COMMENT_WGS84},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
                                   Dim_Table_Type *dim_table, int *grp_id)
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
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_EPHEM, dim_table->time_ephemeris.len, &dim_table->time_ephemeris.id)))
     {
        _pTIO_err_verror_nc (status, "defining dimension 'time' in group %s", __func__, grp_name);
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
                                   Dim_Table_Type *dim_table, int *grp_id)
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
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_MANEUVER, dim_table->time_maneuvers.len, &dim_table->time_maneuvers.id)))
     {
        _pTIO_err_verror_nc (status, "defining dimension 'time' in group %s", __func__, grp_name);
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
                                   Dim_Table_Type *dim_table, int *grp_id)
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
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_GYRO, dim_table->time_gyroscope.len, &dim_table->time_gyroscope.id)))
     {
        _pTIO_err_verror_nc (status, "defining dimension 'time' in group %s", __func__, grp_name);
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
                                Dim_Table_Type *dim_table, int *grp_id)
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
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_dim (grp, TIO_VAR_NAME_TIME_SMA, dim_table->time_sma.len, &dim_table->time_sma.id)))
     {
        _pTIO_err_verror_nc (status, "defining dimension 'time' in group %s", __func__, grp_name);
        return -1;
     }

   /* time coordinate variable */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"units", "s"},
             TEXT_ATTRS_END
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
             TEXT_ATTRS_END
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
                                   Dim_Table_Type *dim_table, int *grp_id)
{
   int status, grp;

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
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
                                   Dim_Table_Type *dim_table, int *grp_id)
{
   int status, grp;

   if (grp_name == NULL)
     {
        _pTIO_err_verror ("%s:  got NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_grp (parent_grp, grp_name, &grp)))
     {
        _pTIO_err_verror_nc (status, "defining group %s", __func__, grp_name);
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

int TIO_create_l1b_template (int ncid, size_t num_steps, size_t num_xtrack,
                             size_t num_channels)
{
   Dim_Table_Type dim_table;

   dim_table.channel.len = num_channels;
   dim_table.xtrack.len = num_xtrack;
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

   if ((-1 == define_globals (ncid, &dim_table))
       || (-1 == define_band_group (ncid, TIO_GRP_NAME_BAND1, &dim_table, NULL))
       || (-1 == define_band_group (ncid, TIO_GRP_NAME_BAND2, &dim_table, NULL))
       || (-1 == define_geometry_group (ncid, TIO_GRP_NAME_GEOMETRY, &dim_table, NULL))
       || (-1 == define_inr_input_group (ncid, TIO_GRP_NAME_INRINPUT, &dim_table, NULL)))
     {
        _pTIO_err_verror ("%s failed", __func__);
        return -1;
     }

   return 0;
}
