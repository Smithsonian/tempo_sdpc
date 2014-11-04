#ifndef __TIO_INTERNAL_INCLUDE__
#define __TIO_INTERNAL_INCLUDE__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

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
   float value;
}
_pFloat_Attr_Type;
#define _pFLOAT_ATTRS_END  {NULL,0}

typedef struct
{
   char *name;
   int value;
}
_pEnum_Type;
#define _pENUM_TABLE_END {NULL,0}

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

extern int _pTIO_define_enum (int grp, const char *name,
                              const _pEnum_Type *enum_table, int *enum_typeid);

extern int _pTIO_define_dims_using_offsets (int grp,
                                            const _pDim_Offsets_Type *offsets,
                                            _pDim_Table_Type *dim_table);

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

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
