#ifndef __REGRID_H__
#define __REGRID_H__

extern Pixel_Regrid_Type *
Regrid_open (const Pixel_Grid_Param_Type *dest,
             char **files, int num_files, const char *lonlat_grp,
             int *src_num_steps, int *src_num_xtrack);

extern void Regrid_close (Pixel_Regrid_Type *r);

#endif
