#ifndef __REGRID_PIXEL_H__
#define __REGRID_PIXEL_H__ 1

typedef struct
{
   double xmin, xmax;
   double ymin, ymax;
   int nx;
   int ny;
}
Pixel_Grid_Param_Type;

typedef struct
{
   int num_overlaps;
   double min, max;
}
Pixel_Overlap_Info_Type;

typedef struct Pixel_List_Type Pixel_List_Type;
typedef struct Pixel_Regrid_Type Pixel_Regrid_Type;

extern Pixel_List_Type *Pixel_list_new (int num_polys, int num_sides);
extern void Pixel_list_free (Pixel_List_Type *g);
extern int Pixel_list_use_src_index (Pixel_List_Type *p);

extern int
Pixel_list_set_src_index (Pixel_List_Type *lst, int i, int src_index);
extern int
Pixel_list_set_vertices (Pixel_List_Type *lst, int pix, int n,
                         const double *x, const double *y);

extern int
Pixel_grid_arrays (const Pixel_Grid_Param_Type *g,
                   double **x_corners, double **y_corners);

extern Pixel_Regrid_Type *
Pixel_open_regrid (const Pixel_Grid_Param_Type *dest,
                   const Pixel_List_Type *dest_area);

extern int
Pixel_find_overlaps (Pixel_Regrid_Type *r,
                     const Pixel_List_Type *src_area,
                     const Pixel_List_Type *src_lookup);

extern int
Pixel_regrid (const Pixel_Regrid_Type *r, double *src, int *src_mask,
              double *dest, Pixel_Overlap_Info_Type *info);

extern void Pixel_close_regrid (Pixel_Regrid_Type *r);

#endif
