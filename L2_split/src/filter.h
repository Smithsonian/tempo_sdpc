#ifndef __L2_SPLIT_FILTER_H__
#define __L2_SPLIT_FILTER_H__ 1

#include <pixel.h>
#include <libconfig.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file filter.h
 * @brief Filter the vertical stratospheric column estimate
 *
 * The initial estimate of the vertical stratospheric column is filtered
 * by first smoothing a coarsely gridded image to fill in masked areas,
 * clipping hot-spots, and then interpolating the smooth, coarse-gridded
 * image back onto the original grid.
 */

typedef struct Filter_Type Filter_Type;

struct Filter_Type
{
   void (*filter_delete)(Filter_Type *);
   /**<  Free resources associated with Filter_Type object
    * @param[in]  ft  pointer to Filter_Type object from filter_open
    */

   int (*filter_apply)(const Filter_Type *, const Pixel_Grid_Param_Type *,
                       double *);
   /**<  Apply filter to specified array
    * @param[in]  ft    pointer to Filter_Type object from filter_open
    * @param[in]  mesh  pointer to Pixel_Grid_Param_Type grid definition
    * @param[in/out] mesh_vert_strat  pointer to an array containing the
    *                 estimated vertical stratospheric column.  It is
    *                 interpreted as a 2D array of size mesh->nx*mesh->ny,
    *                 where the X coordinate is the fastest varying index.
    *                 On output, the array is overwritten by the filtered
    *                 result.
    * @return 0 on success, -1 on error.
    */

#ifdef FILTER_PRIVATE_DATA
   FILTER_PRIVATE_DATA
#endif
};

/** Allocate and initialize an object of type Filter_Type
 * @param[in] cfg   pointer to \c config_t object associated with
 *                  \c L2_split parameter file
 * @return pointer to initialized Filter_Type on success, NULL on error
 */
extern Filter_Type *filter_open (config_t *cfg);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
