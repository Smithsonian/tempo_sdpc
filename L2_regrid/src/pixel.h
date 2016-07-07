#ifndef __REGRID_PIXEL_H__
#define __REGRID_PIXEL_H__ 1

/** @file pixel.h
 *  @brief Manipulate pixel lists
 */

enum
{
   VALUE_IS_DOUBLE = 0,
   VALUE_IS_BYTE,
   VALUE_IS_UBYTE,
   VALUE_IS_SHORT,
   VALUE_IS_USHORT,
   VALUE_IS_INT,
   VALUE_IS_UINT,
   VALUE_IS_INT64,
   VALUE_IS_UINT64
};

typedef struct
{
   double xmin, xmax;
   double ymin, ymax;
   int nx;
   int ny;
}
Pixel_Grid_Param_Type;

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
Pixel_regrid (const Pixel_Regrid_Type *r, const int *src_mask,
              double fill_value, const double *src, double *dest);
extern int
Pixel_regrid_bytes (const Pixel_Regrid_Type *r, const int *src_mask,
                    int value_type, const void *fill_value,
                    const void *src, void *dest);

extern void Pixel_close_regrid (Pixel_Regrid_Type *r);

#endif
