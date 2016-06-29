#ifndef __REGRID_H__
#define __REGRID_H__

extern Pixel_Regrid_Type *
Regrid_open (const char **files, int num_files, 
             const Pixel_Grid_Param_Type *dest, int *src_dims);

extern void Regrid_close (Pixel_Regrid_Type *r);

#endif
