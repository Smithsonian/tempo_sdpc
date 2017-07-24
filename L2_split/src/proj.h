#ifndef __L2_SPLIT_PROJ_H__
#define __L2_SPLIT_PROJ_H__ 1

#include <pixel.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file proj.h
 * @brief Interface for performing pixel coordinate projections
 */

/** Project longitude-latitude coordinates to Albers equal-area coordinates
 * @param[in/out]   lon    On input, longitude [deg], on output, Albers X-coordinate [meters]
 * @param[in/out]   lat    On input, latitude [deg], on output, Albers Y-coordinate [meters]
 * @param[in]       n      Number of (longitude,latitude) points
 * @return 0 on success, -1 on error
 */
extern int proj_longlat_to_albers (double *lon, double *lat, int n);

typedef int Proj_cvt_type (double *, double *, int);

/** Generate a pixel list for specified mesh parameters
 * @param[in]  mesh   Pointer to a \c Pixel_Grid_Param_Type structure containing
 *                    the mesh specification
 * @param[in]  cvt    Pointer to a coordinate projection to be applied to all
 *                    coordinates in the pixel list.
 * @return  On success, a pointer to a \c Pixel_List_Type structure,
 *          On error, a NULL pointer.
 */
extern Pixel_List_Type *
proj_pixel_list (const Pixel_Grid_Param_Type *mesh, Proj_cvt_type *cvt);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
