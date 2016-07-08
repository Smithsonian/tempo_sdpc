#ifndef __REGRID_H__
#define __REGRID_H__

/** @file regrid.h
 *  @brief Construct pixel lists for all granules in a scan.
 *         Compute overlap between input pixels and target grid pixels.
 */

extern Pixel_Regrid_Type *
Regrid_open (const Pixel_Grid_Param_Type *dest,
             char **files, int num_files, const char *lonlat_grp);

extern void Regrid_close (Pixel_Regrid_Type *r);

#endif
