#ifndef __SLIT_FUNCTION_TABLE_INCLUDE__
#define __SLIT_FUNCTION_TABLE_INCLUDE__

typedef struct SF_Table_Type SF_Table_Type;

struct SF_Table_Type
{
   int (*stt_close)(SF_Table_Type *);
   int (*stt_size) (const SF_Table_Type *, int *num_xtrack, int *num_waves, int *num_params);
   int (*stt_get_params) (const SF_Table_Type *, int xtrack, double wave0, double *params);

#ifdef SLIT_FUNCTION_TABLE_PRIVATE_DATA
   SLIT_FUNCTION_TABLE_PRIVATE_DATA
#endif
};

extern SF_Table_Type *sf_table_open (const char *sf_file, const char *band_name);

#endif
