#ifndef __TIO_INTERNAL_INCLUDE__
#define __TIO_INTERNAL_INCLUDE__ 1

/* placeholder values -- FIXME */
#define _pTIO_PIXEL_XSIZE    18.0  /* micron */
#define _pTIO_PIXEL_YSIZE    18.0
#define _pTIO_PIXEL_XSCALE   55.0  /* microradian */
#define _pTIO_PIXEL_YSCALE   55.0

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
   int value;
}
_pInt_Attr_Type;
#define _pINT_ATTRS_END  {NULL,0}

typedef struct
{
   char *name;
   unsigned int value;
}
_pEnum_Type;
#define _pENUM_TABLE_END {NULL,0}

typedef struct
{
   int id;      /* dim id from nc_def_dim */
   size_t len;  /* fixed length of this dimension or NC_UNLIMITED */
}
_pDim_Type;

extern void _pTIO_err_verror (const char *fmt, ...);
extern void _pTIO_err_verror_nc (int status, const char *fmt, ...);
extern int _pTIO_check_verror_nc (int status, int line, const char *file);

extern int _pTIO_define_enum (int grp, const char *name,
                              const _pEnum_Type *enum_table, int *enum_typeid);

extern int _pTIO_define_int_attrs (int grp, int varid, const _pInt_Attr_Type *attrs);
extern int _pTIO_define_text_attrs (int grp, int varid, const _pText_Attr_Type *attrs);
extern int _pTIO_define_var_with_text_attrs (int grp, const char *var_name, nc_type xtype,
                                             int num_dims, const int *dimids,
                                             const _pText_Attr_Type *text_attrs,
                                             int *pvarid);

extern int _pTIO_put_fillvalue_attr (int grp, int varid, nc_type xtype);
extern int _pTIO_define_processing_level (int grp, enum TIO_Processing_Level level);
#endif
