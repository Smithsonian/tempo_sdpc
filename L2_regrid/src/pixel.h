#ifndef __REGRID_PIXEL_H__
#define __REGRID_PIXEL_H__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file pixel.h
 * @brief Create and manipulate pixel lists.
 *
 * A pixel list is more general than a grid specification
 * because a pixel list treats each pixel as an individual
 * polygon, allowing for gaps and for arbitrary overlaps.
 */

#include <float.h>

#define PIXEL_INIT_NUM_SAMPLES  (-1)
#define PIXEL_INIT_MIN_SAMPLE   (DBL_MAX)
#define PIXEL_INIT_MAX_SAMPLE   (-DBL_MAX)

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
   double xmin; /**< minimum pixel edge X coordinate */
   double xmax; /**< maximum pixel edge X coordinate */
   double ymin; /**< minimum pixel edge Y coordinate */
   double ymax; /**< maximum pixel edge Y coordinate */
   int nx;      /**< Number of pixels in the X dimension */
   int ny;      /**< Number of pixels in the Y dimension */
   int num_extra_xpoints;
   int num_extra_ypoints;
   /**< Number of extra points per side for each dimension
    * The total number of points defining the perimeter
    * of each pixel is then 4 + 2*num_extra_xpoints + 2*num_extra_ypoints
    */
}
Pixel_Grid_Param_Type;

typedef struct
{
   double *min;    /**< minimum input value contributing to regridded pixel result */
   double *max;    /**< maximum input value contributing to regridded pixel result */
   int *num;       /**< number of input values contributing to regridded pixel result */
}
Pixel_Regrid_Stats_Type;

typedef struct Pixel_List_Type Pixel_List_Type;
typedef struct Pixel_Regrid_Type Pixel_Regrid_Type;

/** Allocate a Pixel_List_Type structure
 * @param[in]  num_polys   The number of pixels (polygons)
 *                         to be stored.
 * @param[in]  num_sides   The number of edges associated with each pixel.
 * @return  A Pixel_List_Type pointer on success, NULL on error
 *
 * When no longer needed, the returned structure should be freed
 * by a call to \ref Pixel_list_free.
*/
extern Pixel_List_Type *Pixel_list_new (int num_polys, int num_sides);

/** Free resources associated with a Pixel_List_Type structure
 * @param[in]  g    The Pixel_List_Type structure
 */
extern void Pixel_list_free (Pixel_List_Type *g);

/** Declare intent to use the src_index field of a Pixel_List_Type struct.
 * @param[in] p  The relevant Pixel_List_Type struct
 * @return 0 on success, -1 on error
 */
extern int Pixel_list_use_src_index (Pixel_List_Type *p);

/** Set the src_index field of a specific pixel in a Pixel_List_Type
 *  struct.
 * @param[in]  lst      A Pixel_List_Type struct
 * @param[in]  i        The index of a pixel
 * @param[in] src_index The src_index value
 * @return 0 on success, -1 on error
 *
 * The src_index field is used to support an additional level of
 * indirection when referring to a pixel grid that may contain many
 * more pixels than the relevant Pixel_List_Type struct.
 *
 * It is necessary to call \ref Pixel_list_use_src_index before
 * attempting to set a src_index value with this routine.
 */
extern int
Pixel_list_set_src_index (Pixel_List_Type *lst, int i, int src_index);

/** Define the X,Y vertices of a specific pixel in a Pixel_List_Type
 *  struct.
 * @param[in]  lst      A Pixel_List_Type struct
 * @param[in]  pix      The index of a pixel.
 * @param[in]  n        The number of X,Y vertices
 * @param[in]  x        The vertex X coordinates
 * @param[in]  y        The vertex Y coordinates
 * @return 0 on success, -1 on failure
 *
 * The coordinate vertices are normally provided in
 * counter-clockwise order with no duplicates.
 */
extern int
Pixel_list_set_vertices (Pixel_List_Type *lst, int pix, int n,
                         const double *x, const double *y);

extern int
Pixel_list_pack (Pixel_List_Type *pixel_list,
                 double *xs, double *ys, int num_pixels,
                 int num_pixel_vertices,
                 int *step, int num_xtrack);

/** Use grid parameters to generate X,Y pixel corner arrays
 * @param[in]   g          A grid parameter structure
 * @param[out]  x_corners  Output array of pixel corner X coordinates
 * @param[out]  y_corners  Output array of pixel corner Y coordinates
 * @return 0 on success, -1 on error
 *
 * The pixels are ordered so that the X-coordinate varies fastest.
 * For each pixel, the coordinates of the four corners are
 * packed in counter-clockwise order.
 */
extern int
Pixel_grid_arrays (const Pixel_Grid_Param_Type *g,
                   double **x_corners, double **y_corners);

/** Prepare to perform regridding onto a specified destination grid.
 * @param[in]  dest       The destination grid parameters.
 * @param[in]  dest_area  The destination grid pixel list with
 *                        pixel vertices in coordinates suitable
 *                        for computing areas (NULL if not needed).
 * @return A Pixel_Regrid_Type structure on success, NULL on error.
 *
 * If the destination grid parameters are Cartesian coordinates
 * that may be used to compute polygon areas, the \c dest_area
 * parameter is not needed and may be NULL.
 *
 * Note that the memory for the input parameters is managed by
 * the calling routine.
 * When no longer needed, the memory for the Pixel_Regrid_Type
 * structure should be freed by \ref Pixel_close_regrid.
 */
extern Pixel_Regrid_Type *
Pixel_open_regrid (const Pixel_Grid_Param_Type *dest,
                   const Pixel_List_Type *dest_area);

/** Free resources allocated by \ref Pixel_open_regrid
 * @param[in]  r      Pixel_Regrid_Type structure allocated
 *                    by \ref Pixel_open_regrid
 */
extern void Pixel_close_regrid (Pixel_Regrid_Type *r);

/** Track the largest source coordinates yet referenced
 * @param   r               Pixel_Regrid_Type structure
 * @param   src_max_step    Max mirror_step so far
 * @param   src_max_xtrack  Max xtrack pixel so far
 *
 * When the input source grid is read in pieces, the code
 * infers the source grid's full size by tracking the largest
 * coordinates yet referenced.
 */
extern void Pixel_regrid_grow_srcdims (Pixel_Regrid_Type *r,
                                       int src_max_step, int src_max_xtrack);

/** Query the inferred minimum source grid dimensions
 * @param   r               Pixel_Regrid_Type structure
 * @param   num_src_step    Minimum inferred mirror step dimension
 * @param   num_src_xtrack  Minimum inferred xtrack dimension
 */
extern void Pixel_regrid_get_srcdims (const Pixel_Regrid_Type *r,
                                      int *num_src_step, int *num_src_xtrack);

/** Generate a polygon overlap structure containing, for each
 *  destination pixel, a list of all overlapping source pixels
 *  and the corresponding overlap area.
 *
 * @param[in]  r      Pixel_Regrid_Type structure allocated
 *                    by Pixel_open_regrid
 * @param[in]  src_area  Source pixel list with polygon vertices
 *                    in coordinates suitable for computing areas
 *                    [e.g. Albers projected coordinates]
 * @param[in]  src_lookup  Source pixel list with polygon vertices
 *                    in coordinates that are a linear function of the
 *                    destination grid pixel list array index
 *                    [e.g. longitude-latitude] (NULL if not needed).
 *
 * @return 0 on success, -1 on error
 *
 * The polygon overlap structure is stored in a field private
 * to the Pixel_Regrid_Type structure.
 *
 * For each input source pixel, the source pixel bounding box is
 * used to compute the array indices of destination pixels that
 * might overlap with it.  If the vertex coordinates in the
 * \c src_area pixel list can be used in this way, then the
 * \c src_lookup pixel list is not needed, and the src_lookup
 * parameter may be NULL.  Otherwise, the \c src_lookup list
 * is used for bounding-box coordinate lookup, and the \c src_area
 * list is used only for computing pixel areas.
 */
extern int
Pixel_find_overlaps (Pixel_Regrid_Type *r,
                     const Pixel_List_Type *src_area,
                     const Pixel_List_Type *src_lookup);

/** Perform regridding of variable values using the pixel overlap
 *  structure computed by \ref Pixel_find_overlaps.
 *
 * @param[in]  r      Pixel_Regrid_Type structure allocated
 *                    by \ref Pixel_open_regrid, and with (private)
 *                    pixel overlap structure initialized by
 *                    \ref Pixel_find_overlaps.
 * @param[in]  src_mask  An mask array with one element per
 *                    pixel in the input spatial grid.
 *                    Non-zero elements indicate that the
 *                    corresponding variable value should not
 *                    be used.
 * @param[in]  fill_value  Value used to initialize each pixel
 *                    in the detination array.
 * @param[in]  src    Array of variable values with one element
 *                    per pixel in the input spatial grid.
 * @param[out] dest   Array of regridded variable values, with
 *                    one element per pixel in the destination
 *                    spatial grid.
 * @return 0 on success, -1 on failure.
 *
 * Destination grid pixels that do not overlap the input source
 * grid are set to the value provided in \c fill_value.
 */
extern int
Pixel_regrid (const Pixel_Regrid_Type *r, const int *src_mask,
              double fill_value, const double *src, double *dest,
              Pixel_Regrid_Stats_Type *rs);

/** Free memory associated with a \ref Pixel_Regrid_Stats_Type object.
 * @param[in]  rs    Pointer to a \ref Pixel_Regrid_Stats_Type object,
 *                   allocated by \ref Pixel_alloc_regrid_stats.
 */
extern void Pixel_free_regrid_stats (Pixel_Regrid_Stats_Type *rs);

/** Allocate a \ref Pixel_Regrid_Stats_Type object.
 * @param[in]  num_pixels             Number of spatial pixels.
 * @param[in]  num_values_per_pixel   Number of values per spatial pixel.
 * @return Pointer to an object of type \ref Pixel_Regrid_Stats_Type on success,
 *         NULL on error
 *
 * The second argument is provided to facilitate regridding multidimensional
 * objects.  For example, consider regridding a 4-dimensional array with
 * dimensions [ny, nx, na, nb].  The spatial dimensions are [ny, nx],
 * so that num_pixels = nx * ny, and num_values_per_pixel = na * nb.
 * It is assumed that the spatial coordinates vary slowest.
 */
extern Pixel_Regrid_Stats_Type *
Pixel_alloc_regrid_stats (int num_pixels, int num_values_per_pixel);

/** Perform regridding of bitfield values using the pixel overlap
 *  structure computed by Pixel_find_overlaps.
 *
 * @param[in]  r      Pixel_Regrid_Type structure allocated
 *                    by \ref Pixel_open_regrid, and with (private)
 *                    pixel overlap structure initialized by
 *                    \ref Pixel_find_overlaps.
 * @param[in]  src_mask  An mask array with one element per
 *                    pixel in the input spatial grid.
 *                    Non-zero elements indicate that the
 *                    corresponding variable value should not
 *                    be used.
 * @param[in]  value_type  Data type of the input bitfield
 * @param[in]  fill_value  Value used to initialize each pixel
 *                    in the detination array.
 * @param[in]  src    Array of bitfield values with one element
 *                    per pixel in the input spatial grid.
 * @param[out] dest   Array of regridded bitfield values, with
 *                    one element per pixel in the destination
 *                    spatial grid.
 * @return 0 on success, -1 on failure.
 *
 * Destination grid pixels that do not overlap the input source
 * grid are set to the value provided in \c fill_value.
 */
extern int
Pixel_regrid_bytes (const Pixel_Regrid_Type *r, const int *src_mask,
                    int value_type, const void *fill_value,
                    const void *src, void *dest);

extern int
Pixel_regrid_from_mesh (const Pixel_Regrid_Type *r, const int *mesh_mask,
                        double fill_value, const double *mesh_values,
                        double *values);

extern int *Pixel_regrid_overlap_map (const Pixel_Regrid_Type *r);

/* Debugging tools */

extern int __Pixel_print_overlap (const Pixel_Regrid_Type *r, int dest_pixel_index);
extern void __Pixel_verbose_output (int);
extern void __Pixel_verbose_output_window (double x, double y, double dx, double dy);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
