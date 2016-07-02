#ifndef __REGRID_VAR_H
#define __REGRID_VAR_H 1

typedef struct Var_Value_Buffer_Type Var_Value_Buffer_Type;

extern Var_Value_Buffer_Type *
Var_new_value_buffer (int dest_nx, int dest_ny, int src_num_step, int src_num_xtrack);
extern void Var_free_value_buffer (Var_Value_Buffer_Type *vb);

extern int Var_write_lonlat_grid (int ncid, const Pixel_Grid_Param_Type *dest);
extern int Var_write_values (int ncid, const Var_Value_Buffer_Type *vb,
                             int ncid_infile, const char *var_name);
extern int Var_apply_regrid (const Pixel_Regrid_Type *r, Var_Value_Buffer_Type *vb,
                             const char *var_name, const char **files, int num_files);

#endif 
