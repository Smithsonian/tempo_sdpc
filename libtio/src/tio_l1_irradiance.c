/** @file
 *  @brief TEMPO Level 1 irradiance template generation
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

#define TIO_CHUNKSIZE_XTRACK 128
#define TIO_CHUNKSIZE_STEP 1
#define DO_CHUNKING        1

/* An instance of a _pDim_Table_Type struct is used as a lookup table
 * for all the dimensions that are defined anywhere in the associated
 * netCDF file.
 */
struct _pDim_Table_Type
{
   _pDim_Type channel;           /* dispersion direction */
   _pDim_Type xtrack;            /* pixel north-south spatial coordinate */
   _pDim_Type step;              /* mirror step position */
};

static int define_global_dims (int grp, _pDim_Table_Type *dim_table)
{
   static _pDim_Offsets_Type dim_offsets[] =
    {
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_STEP,step),
       _pDIM_OFFSETS_END
    };

   return _pTIO_define_dims_using_offsets (grp, dim_offsets, dim_table);
}

static int define_irradiance_group_dims (int grp, _pDim_Table_Type *dim_table)
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
   int varid, dims[TIO_MAX_VAR_DIMS];

   /* coordinate variables */
   if (0 != tio_define_dim_step_var (grp, dim_table->step.id))
     return -1;

   /* time */
     {
        static _pText_Attr_Type time_attrs[] =
          {
             {"standard_name", "time"},
             {"long_name", "exposure start time"},
             {"calendar", "gregorian"},
             _pTEXT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if ((-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_TIME, NC_DOUBLE, 1, dims, time_attrs, NULL)
             || (0 != tio_write_timestamp_unit_string (grp, TEMPO_VAR_TIME))))
          return -1;
     }

   /* exposure_time */
     {
        static _pText_Attr_Type exposure_time_attrs[] =
          {
             {"units", "seconds"},
             {"long_name", "exposure duration"},
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

   /* solar boresight angles (phi, theta) */
     {
        static _pText_Attr_Type solar_phi_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", "solar boresight azimuthal angle"},
             _pTEXT_ATTRS_END
          };
        static _pText_Attr_Type solar_theta_attrs[] =
          {
             {"units", "degrees"},
             {"long_name", "solar boresight polar angle"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type solar_phi_float_attrs[] =
          {
             {"valid_min", -180.0},
             {"valid_max",  180.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        static _pFloat_Attr_Type solar_theta_float_attrs[] =
          {
             {"valid_min",    0.0},
             {"valid_max",  180.0},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        dims[0] = dim_table->step.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SOLAR_BORESIGHT_PHI, NC_FLOAT, 1, dims, solar_phi_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, solar_phi_float_attrs))
          return -1;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_SOLAR_BORESIGHT_THETA, NC_FLOAT, 1, dims, solar_theta_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, solar_theta_float_attrs))
          return -1;
     }

   /* earth_sun_distance */
   if (0 != tio_set_earth_sun_distance (grp, _pTIO_EARTH_SUN_DISTANCE))
     return -1;

   return 0;
}

static int set_default_header_timestamps (int grp)
{
   char buf[MAX_ISOTIME_LEN];
   int len, status;

   if (0 != TIO_mktimestamp_str (0.0, 1, buf, sizeof(buf)))
     return -1;

   len = strlen(buf)+1;

   if (NC_NOERR != (status = nc_put_att_text (grp, NC_GLOBAL, "time_coverage_start", len, buf)))
     goto report_error;
   if (NC_NOERR != (status = nc_put_att_text (grp, NC_GLOBAL, "time_coverage_end", len, buf)))
     goto report_error;

   return 0;
report_error:
   tell_verror (TELL_IO_WRITE_ERROR,
                "%s: writing timestamp attributes (%s)", __func__, nc_strerror(status));
   return -1;
}

static int define_global_attrs (int grp, const char *product_type)
{
   static _pText_Attr_Type text_attrs[] =
     {
        {"Conventions", TIO_FORMAT_CONVENTIONS},
        _pTEXT_ATTRS_END
     };
   static _pInt_Attr_Type int_attrs[] =
     {
        MAKE_INT_ATTR1("format_version", TIO_L1_FORMAT_VERSION),
        _pINT_ATTRS_END
     };

   if ((-1 == tio_write_epoch_timestamp (grp, NC_GLOBAL))
       || (-1 == _pTIO_define_text_attrs (grp, NC_GLOBAL, text_attrs)))
     return -1;

   if (0 != set_default_header_timestamps (grp))
     return -1;

   if (-1 == _pTIO_define_int_attrs (grp, NC_GLOBAL, int_attrs))
     return -1;

   if (0 != TIO_label_product (grp, product_type, 1, 0))
     return -1;

   return 0;
}

static int define_irradiance_group (int parent_grp, TIO_Scan_Group_Type *sg,
                                    _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];
   int shuffle, deflate, deflate_level;
#ifdef DO_CHUNKING
   int storage = NC_CHUNKED;
   size_t chunksizes[TIO_MAX_VAR_DIMS];
#endif

   _pTIO_get_level1_compression (&deflate, &deflate_level, &shuffle);

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

   if (-1 == define_irradiance_group_dims (grp, dim_table))
     return -1;

   /* group-local coordinate variables */
   dims[0] = dim_table->xtrack.id;
   status = nc_def_var (grp, TEMPO_DIM_XTRACK, NC_INT, 1, dims, &varid);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining coordinate variable %s (%s)",
                     __func__, TEMPO_DIM_XTRACK, nc_strerror(status));
        return -1;
     }
   if (0 != _pTIO_emit_xtrack_indices (grp, varid, dim_table->xtrack.len))
     return -1;

   /* irradiance */
     {
        static _pText_Attr_Type irradiance_attrs[] =
          {
             {"units", _pTIO_IRRADIANCE_UNITS},
             {"ancillary_variables", TEMPO_VAR_IRRADIANCE_ERROR},
             _pTEXT_ATTRS_END
          };
        float irradiance_fill = TIO_FILL_IRRADIANCE;

        /* It's convenient to make the irradiance a 3D object */
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_IRRADIANCE, NC_FLOAT, 3, dims, irradiance_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_att (grp, varid, _FillValue, NC_FLOAT, 1, &irradiance_fill)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "writing %s fill value to grp=%d (%s)",
                          TEMPO_VAR_IRRADIANCE, grp, nc_strerror(status));
             return -1;
          }
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "defining %s compression parameters for grp = %d (%s)",
                          TEMPO_VAR_IRRADIANCE, grp, nc_strerror(status));
             return -1;
          }
#ifdef DO_CHUNKING
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = ((dim_table->xtrack.len < TIO_CHUNKSIZE_XTRACK) ?
                         dim_table->xtrack.len : TIO_CHUNKSIZE_XTRACK);
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (0 != TIO_def_var_chunking (grp, varid, storage, chunksizes)))
          return -1;
#endif
     }

   /* irradiance error */
     {
        static _pText_Attr_Type irradiance_error_attrs[] =
          {
             {"units", _pTIO_IRRADIANCE_UNITS},
             {"long_name", "irradiance error"},
             _pTEXT_ATTRS_END
          };
        float irradiance_error_fill = TIO_FILL_IRRADIANCE_ERROR;

        /* It's convenient to make the irradiance a 3D object */
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_IRRADIANCE_ERROR, NC_FLOAT, 3, dims, irradiance_error_attrs, &varid))
          return -1;
        if (NC_NOERR != (status = nc_put_att (grp, varid, _FillValue, NC_FLOAT, 1, &irradiance_error_fill)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "writing %s fill value to grp=%d (%s)",
                          TEMPO_VAR_IRRADIANCE_ERROR, grp, nc_strerror(status));
             return -1;
          }
        if (NC_NOERR != (status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level)))
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "defining %s compression parameters for grp = %d (%s)",
                          TEMPO_VAR_IRRADIANCE_ERROR, grp, nc_strerror(status));
             return -1;
          }
#ifdef DO_CHUNKING
        chunksizes[0] = TIO_CHUNKSIZE_STEP;
        chunksizes[1] = ((dim_table->xtrack.len < TIO_CHUNKSIZE_XTRACK) ?
                         dim_table->xtrack.len : TIO_CHUNKSIZE_XTRACK);
        chunksizes[2] = dim_table->channel.len;
        if ((storage == NC_CHUNKED)
            && (0 != TIO_def_var_chunking (grp, varid, storage, chunksizes)))
          return -1;
#endif
     }

   /* nominal_wavelength */
     {
        static _pText_Attr_Type wavelength_attrs[] =
          {
             {"units", "nm"},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type wavelength_float_attrs[] =
          {
             {"valid_min", -TIO_FILL_FLOAT},
             {"valid_max", TIO_FILL_FLOAT},
             {_FillValue, TIO_FILL_FLOAT},
             _pFLOAT_ATTRS_END
          };
        int num_dims = TIO_NOMINAL_WAVELEN_NUM_DIMS;
        switch (num_dims)
          {
           case 1:
             dims[0] = dim_table->channel.id;
             break;
           default:
             num_dims = 2;
             dims[0] = dim_table->xtrack.id;
             dims[1] = dim_table->channel.id;
             break;
          }
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_WAVELEN_NOMINAL, NC_FLOAT, num_dims, dims, wavelength_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, wavelength_float_attrs))
          return -1;
     }

   /* pixel quality flags */
   if (-1 == _pEmit_Var_Pixel_Quality_Flag (grp, dim_table))
     return -1;

   if (grp_id != NULL)
     {
        *grp_id = grp;
     }

   return 0;
}

int TIO_l1_irradiance_template (int ncid, const char *product_type, size_t num_steps,
                                int num_sgrps, TIO_Scan_Group_Type *sgrps)
{
   _pDim_Table_Type dim_table;
   int i;

   memset ((char *)&dim_table, 0, sizeof (dim_table));

   /* Initialize the dimension sizes that are known at this point.
    * Other dimensions are group-specific and are initialized only
    * when those groups are being defined.
    */
   dim_table.step.len = num_steps;

   if ((-1 == define_global_attrs (ncid, product_type))
       || (-1 == define_global_dims (ncid, &dim_table))
       || (-1 == define_global_vars (ncid, &dim_table)))
     {
        Tell_verror (TELL_UNKNOWN_ERROR, "%s failed", __func__);
        return -1;
     }

   for (i = 0; i < num_sgrps; i++)
     {
        if (-1 == define_irradiance_group (ncid, &sgrps[i], &dim_table, NULL))
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s failed defining irradiance group %d",
                          __func__, i);
             return -1;
          }
     }

   return 0;
}
