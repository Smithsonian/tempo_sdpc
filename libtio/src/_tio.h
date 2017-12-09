/*! @file
 *  @brief Internal interfaces (private)
 */
#ifndef __TIO_INTERNAL_INCLUDE__
#define __TIO_INTERNAL_INCLUDE__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <time.h>
#include <netcdf.h>

/* placeholder values -- FIXME */
#define _pTIO_PIXEL_SCALE_ROW      0.195  /* nm */
#define _pTIO_PIXEL_SCALE_COLUMN   55.0   /* microradian */
#define _pTIO_MIRROR_STEP_SIZE    114.0   /* microradian */

/* default values for template files */
#define _pTIO_TIME_COVERAGE_START   "2017-01-01T12:00:00Z"
#define _pTIO_TIME_COVERAGE_END     "2017-01-01T13:00:00Z"
#define _pTIO_EARTH_SUN_DISTANCE    1.47975e+11
#define _pTIO_PHOTON_UNITS "photons/s/cm^2/nm/sr"

#define TIO_REALLOC realloc
#define TIO_MALLOC malloc
#define TIO_FREE free

typedef struct _pName_UInt_Pair_Type _pName_UInt_Pair_Type;
struct _pName_UInt_Pair_Type
{
   unsigned int value;
   char name[TIO_MAX_NAME_LEN];
};
#define _pNAME_UINT_LIST_END  {0,""}

typedef struct
{
   char *name;
   char *text;
}
_pText_Attr_Type;
#define _pTEXT_ATTRS_END  {NULL,NULL}

typedef struct
{
   char *name;
   int num_values;
#define MAX_NUM_INT_ATTRS  8
   int value[MAX_NUM_INT_ATTRS];
}
_pInt_Attr_Type;
#define _pINT_ATTRS_END   {NULL,0, {0,0,0,0,0,0,0,0}}
#define MAKE_INT_ATTR1(n,v0) {n, 1, {v0,0,0,0,0, 0,0,0}}

typedef struct
{
   char *name;
   short value;
}
_pShort_Attr_Type;
#define _pSHORT_ATTRS_END  {NULL,0}

typedef struct
{
   char *name;
   float value;
}
_pFloat_Attr_Type;
#define _pFLOAT_ATTRS_END  {NULL,0}

typedef struct
{
   int id;      /**< index assigned by nc_def_dim */
   size_t len;  /**< fixed dimension size or NC_UNLIMITED */
}
_pDim_Type;

typedef struct _pDim_Table_Type _pDim_Table_Type;

typedef struct
{
   char *name;
   size_t len_offset;
   size_t id_offset;
}
_pDim_Offsets_Type;
#define _pDIM_OFFSETS_END {NULL,0,0}
#define _pDIM_OFFSET_ENTRY(name,field) \
   {name, \
        (offsetof(_pDim_Table_Type,field) + offsetof(_pDim_Type,len)), \
        (offsetof(_pDim_Table_Type,field) + offsetof(_pDim_Type,id))}

extern int _pTIO_define_dims_using_offsets (int grp,
                                            const _pDim_Offsets_Type *offsets,
                                            _pDim_Table_Type *dim_table);

extern int _pTIO_define_short_attrs (int grp, int varid,
                                     const _pShort_Attr_Type *attrs);
extern int _pTIO_define_int_attrs (int grp, int varid,
                                   const _pInt_Attr_Type *attrs);
extern int _pTIO_define_float_attrs (int grp, int varid,
                                     const _pFloat_Attr_Type *attrs);
extern int _pTIO_define_text_attrs (int grp, int varid,
                                    const _pText_Attr_Type *attrs);

extern int _pTIO_define_var_with_text_attrs (int grp, const char *var_name, nc_type xtype,
                                             int num_dims, const int *dimids,
                                             const _pText_Attr_Type *text_attrs,
                                             int *pvarid);

extern int _pTIO_put_fillvalue_attr (int grp, int varid, nc_type xtype);

extern int _pTIO_define_processing_level (int grp, int level);

extern int _pTIOMake_Name_UInt_Arrays (_pName_UInt_Pair_Type *array,
                                       int *pnum_values, char **pnames,
                                       unsigned int **pvalues);
extern int _pEmit_Var_Pixel_Quality_Flag (int grp,
                                          _pDim_Table_Type *dim_table);

#define MAX_ISOTIME_LEN 32
/*      MAX_ISOTIME_LEN must hold:  yyyy-mm-ddThh:mm:ss.sssZ
 * (Although floating point seconds are currently unsupported
 *  by strftime/strptime and may never be used)
 */

typedef struct _pTIO_Granule_Ident_Type _pTIO_Granule_Ident_Type;
struct _pTIO_Granule_Ident_Type
{
   _pTIO_Granule_Ident_Type *next;
   int scan_seq_num;
   int granule_seq_num;
   int granule_num;
   char tstart_str[MAX_ISOTIME_LEN];   /**< start time, UTC, ISO 8601 time string */
   char tend_str[MAX_ISOTIME_LEN];     /**< end time, UTC, ISO 8601 time string */
   double tstart;  /**< start time, TAI seconds since TEMPO epoch */
   double tend;    /**< end time, TAI seconds since TEMPO epoch */
};

extern int _pTIO_read_granule_ident (int ncid, _pTIO_Granule_Ident_Type *gid);
extern int _pTIO_write_granule_ident (int ncid, const _pTIO_Granule_Ident_Type *gid);
extern int _pTIO_parse_timestr (const char *timestr, struct tm *ptm);
extern int _pTIO_tempo_time_from_utc_timestr (const char *timestr, double *tai_sec);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
