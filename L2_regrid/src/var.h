#ifndef __REGRID_VAR_H
#define __REGRID_VAR_H 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file var.h
 *  @brief Manage the process of regridding selected variables.
 */

typedef struct Var_Value_Buffer_Type Var_Value_Buffer_Type;

/** Allocate a buffer to be faciliate regridding of variable values.
 * @param[in]  dest_nx  Number of X pixels (longitude) in the destination grid
 * @param[in]  dest_ny  Number of Y pixels (latitude) in the destination grid
 * @param[in]  src_num_step   Number of mirror steps in the input source grid
 * @param[in]  src_num_xtrack  Number of cross-track pixels in the input source
 *                            grid.
 * @return A pointer to a Var_Value_Buffer_Type structure on succeess,
 *       NULL on error.
 */
extern Var_Value_Buffer_Type *
Var_new_value_buffer (int dest_nx, int dest_ny, int src_num_step,
                      int src_num_xtrack);

/** Free resources allocated by \ref Var_new_value_buffer
 * @param[in]  vb  Pointer to structure previously allocated by
 *                 \ref Var_new_value_buffer.
 */
extern void Var_free_value_buffer (Var_Value_Buffer_Type *vb);

/** Write longitude-latitude coordinate variables to a specified
 * netCDF output file.
 * @param[in]  ncid   Integer file identifier of an open
 *                    netCDF output file.
 * @param[in]  lonlat_grp  The full path to the group in the
 *                         output file where the coordinate
 *                         variables should be written.
 * @param[in]  dest   The structure containing the output grid
 *                    parameters.
 * @return 0 on success, -1 on error.
 */
extern int Var_write_lonlat_grid (int ncid, const char *lonlat_grp,
                                  const Pixel_Grid_Param_Type *dest);

/** Write regridded variable values to a netCDF file.
 * @param[in]  ncid   Integer file identifier of an open
 *                    netCDF output file.
 * @param[in]  vb     Var_Value_Buffer_Type structure containing
 *                    the regridded variable values to write out.
 * @param[in]  out_var_path  The full path to the output variable
 *                           in the output netCDF file.
 * @param[in] ncid_infile  Integer file identifier of an open,
 *                         read-only netCDF input file.
 * @param[in] in_var_path   The full path to the input variable
 *                           in the input netCDF files (used to
 *                           enable copying of variable attributes).
 */
extern int Var_write_values (int ncid, const Var_Value_Buffer_Type *vb,
                             const char *out_var_path,
                             int ncid_infile, const char *in_var_path);

/** Regrid a specified variable using a specified polygon overlap
 *  structure.
 * @param[in]  r   Pixel_Regrid_Type structure containing a
 *                 polygon overlap information produced by
 *                 \ref Regrid_open.
 * @param[in]  vb  A \ref Var_Value_Buffer_Type structure allocated
 *                 by \ref Var_new_value_buffer
 * @param[in] value_type  The type of the variable to be regridded.
 * @param[in] var_path  The full path of the variable to be regridded,
 *                      as it appears in the netCDF input files.
 * @param[in] files     An array of input file names.
 * @param[in] num_files  The number of input files.
 * @return 0 on success, -1 on error
 *
 * The regridded result is stored in a field of the Var_Value_Buffer_Type
 * structure.
 */
extern int Var_apply_regrid (const Pixel_Regrid_Type *r,
                             Var_Value_Buffer_Type *vb,
                             int value_type, const char *var_path,
                             char **files, int num_files);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
