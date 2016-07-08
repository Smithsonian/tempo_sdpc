#ifndef __REGRID_H__
#define __REGRID_H__

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file regrid.h
 *  @brief Prepare to regrid a set of Level 2 input files
 *    onto a specified Level 3 destination grid.
 */

/** Prepare to regrid a set of netCDF files onto a specified
 *  destination grid.
 * @param[in]  dest   The destination grid specification.
 * @param[in]  files  An array of input netCDF file names
 * @param[in]  num_files  The number of input netCDF files.
 * @param[in]  lonlat_grp The name of the file group in the
 *                        input files that contains the
 *                        pixel corner coordinates
 * @return A pointer to a Pixel_Regrid_Type structure on success,
 *         NULL on error.
 *
 * This function reads pixel vertices from each of the input files
 * in turn, using \ref Pixel_find_overlaps to populate the pixel overlap
 * field of the \ref Pixel_Regrid_Type data structure.
 */
extern Pixel_Regrid_Type *
Regrid_open (const Pixel_Grid_Param_Type *dest,
             char **files, int num_files, const char *lonlat_grp);

/** Free resources allocated by Regrid_open.
 * @param[in]  r   A \ref Pixel_Regrid_Type structure allocated
 *                 by Regrid_open.
 */
extern void Regrid_close (Pixel_Regrid_Type *r);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
