/** @file
 *  @brief TEMPO Level 1 irradiance template generation
 */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "netcdf.h"
#include "tell.h"
#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

#define IRRADIANCE_RELERR_LOG10_MIN  (-4.0)
#define IRRADIANCE_RELERR_LOG10_MAX  (+2.0)

#define IRRADIANCE_RELERR_COMMENT \
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
       _pDIM_OFFSET_ENTRY(TEMPO_DIM_STEP,step),
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
        _pINT_ATTRS_END
     };

   if (-1 == _pTIO_define_text_attrs (grp, NC_GLOBAL, text_attrs))
     return -1;

   if (-1 == _pTIO_define_int_attrs (grp, NC_GLOBAL, int_attrs))
     return -1;

   if (-1 == _pTIO_define_processing_level (grp, TIO_PROC_LEVEL_1A))
     return -1;

   return 0;
}

static int define_irradiance_group (int parent_grp, TIO_Scan_Group_Type *sg,
                                    _pDim_Table_Type *dim_table, int *grp_id)
{
   int status, grp, varid;
   int dims[TIO_MAX_VAR_DIMS];
   int shuffle, deflate=1, deflate_level=1;

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

   if (-1 == define_irradiance_group_dims (grp, dim_table))
     return -1;

   /* group-local coordinate variables */
   dims[0] = dim_table->xtrack.id;
   status = nc_def_var (grp, TEMPO_DIM_XTRACK, NC_INT, 1, dims, NULL);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining coordinate variable %s (%s)",
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
             {"units", "microradian"},
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

   /* irradiance */
     {
        static _pText_Attr_Type irradiance_attrs[] =
          {
             {"units", "TBD"},
             {"ancillary_variables", TEMPO_VAR_IRRADIANCE_ERROR},
             _pTEXT_ATTRS_END
          };
        float irradiance_fill = TIO_FILL_IRRADIANCE;

        /* Motivation for making the irradiance a 3D object:
         * In practice, we expect the irradiance measurement will occur at only
         * a single mirror step position, so the irradiance is naturally a 2D object.
         * However, by making the radiance and irradiance have the same dimensionality,
         * we can enable re-use of some radiance code for processing irradiances.
         */
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
     }

   /* irradiance error */
     {
        static _pText_Attr_Type irradiance_error_attrs[] =
          {
             {"units", "TBD"},
             {"comment", IRRADIANCE_RELERR_COMMENT},
             _pTEXT_ATTRS_END
          };
        static _pFloat_Attr_Type irradiance_error_float_attrs[] =
          {
             {RELERR_MIN_LOG10, TIO_FILL_FLOAT},  /* min must be first */
             {RELERR_MAX_LOG10, TIO_FILL_FLOAT},  /* max must be second */
             {RELERR_MISSING, TIO_FILL_FLOAT},    /* u_missing must be third */
             _pFLOAT_ATTRS_END
          };
        short irradiance_error_fill = -32768;

        irradiance_error_float_attrs[0].value = IRRADIANCE_RELERR_LOG10_MIN;
        irradiance_error_float_attrs[1].value = IRRADIANCE_RELERR_LOG10_MAX;
        irradiance_error_float_attrs[2].value = TIO_FILL_IRRADIANCE_ERROR;

        /* Motivation for making the irradiance a 3D object:
         * In practice, we expect the irradiance measurement will occur at only
         * a single mirror step position, so the irradiance is naturally a 2D object.
         * However, by making the radiance and irradiance have the same dimensionality,
         * we can enable re-use of some radiance code for processing irradiances.
         */
        dims[0] = dim_table->step.id;
        dims[1] = dim_table->xtrack.id;
        dims[2] = dim_table->channel.id;
        if (-1 == _pTIO_define_var_with_text_attrs (grp, TEMPO_VAR_IRRADIANCE_ERROR, NC_SHORT, 3, dims, irradiance_error_attrs, &varid))
          return -1;
        if (-1 == _pTIO_define_float_attrs (grp, varid, irradiance_error_float_attrs))
          return -1;
        if (NC_NOERR != (status = nc_put_att (grp, varid, _FillValue, NC_SHORT, 1, &irradiance_error_fill)))
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
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "defining %s compression parameters for grp %d (%s)",
                          TEMPO_VAR_WAVELENGTH, grp, nc_strerror(status));
             return -1;
          }
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

int TIO_l1_irradiance_template (int ncid, size_t num_steps, int num_sgrps,
                                TIO_Scan_Group_Type *sgrps)
{
   _pDim_Table_Type dim_table;
   int i;

   /* Initialize the dimension sizes that are known at this point.
    * Other dimensions are group-specific and are initialized only
    * when those groups are being defined.
    */
   dim_table.step.len = num_steps;

   if ((-1 == define_global_attrs (ncid))
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
