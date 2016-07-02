#ifndef __REGRID_H__
#define __REGRID_H__

extern Pixel_Regrid_Type *
Regrid_open (const Pixel_Grid_Param_Type *dest,
             const char **files, int num_files,
             int *src_num_steps, int *src_num_xtrack);

extern void Regrid_close (Pixel_Regrid_Type *r);

#endif
