#ifndef __REGRID_H__
#define __REGRID_H__

extern int map_strings (const char **str, int n,
                        int (*do_task)(const char *, void *),
                        void *client_data);

extern Pixel_Regrid_Type *
Regrid_open (const char **files, int num_files, int *src_dims,
             const Pixel_Grid_Param_Type *dest);

extern void Regrid_close (Pixel_Regrid_Type *r);

#endif
