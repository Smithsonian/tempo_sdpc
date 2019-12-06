#ifndef __L2_SPLIT_GRANULE_H__
#define __L2_SPLIT_GRANULE_H__ 1

#include <pixel.h>
#include <libconfig.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file scan.h
 * @brief Manage I/O to Level 2 data product files for granules
 *        from a single TEMPO scan
 */

typedef struct Scan_Type Scan_Type;

typedef struct
{
   double *slant_column; /**< Slant column [num_steps, num_xtrack] */
   double *amf_trop;     /**< Tropospheric air mass factor [num_steps, num_xtrack] */
   double *amf_strat;    /**< Stratospheric air mass factor [num_steps, num_xtrack] */
   double *vert_strat;   /**< Stratospheric vertical column */
   int *data_quality_flag;  /**< Data quality flag [num_steps, num_xtrack] */
   int num_steps;        /**< number of mirror steps = dimension of the slowest varying index */
   int num_xtrack;       /**< number of cross-track pixels = dimension of the fastest varying index */
}
Scan_Vars_Type;

/** Read granule files comprising a single scan
 * @param[in] num_files   Number of granule files
 * @param[in] files       Pointer to an array of granule file names
 * @param[in] cfg         Initialized pointer to config_t structure
 * @return  On success, an opaque pointer to a \c Scan_Type object.
 *          On error, a NULL pointer
 */
extern Scan_Type *scan_read_granules (int num_files, char **files, config_t *cfg);

/** Free resources associated with a Scan_Type object
 * @param[in] st  Pointer to a Scan_Type object created by \c scan_read_grids.
 */
extern void scan_free (Scan_Type *st);

/** Retrieve the scan start/end times
 * @param[in] st  Pointer to a Scan_Type object created by \c scan_read_grids.
 * @param[out] tstart  Scan start time [seconds since TEMPO epoch]
 * @param[out] tend    Scan end time [seconds since TEMPO epoch]
 * @return 0 on success, -1 on error
 */
extern int scan_time_interval (const Scan_Type *st, double *ptstart, double *ptend);

/** Prepare to regrid variables from the scan grid to a uniform mesh.
 * @param[in]   st   Pointer to a \c Scan_Type object initialized by
 *                  \c scan_read_grids
 * @param[in]  mesh  Pointer to a \c Pixel_Grid_Param_Type structure
 *                   defining the target uniform mesh.
 * @return On success, a pointer to a \c Pixel_Regrid_Type object,
 *         on error, a NULL pointer
 *
 * For each pixel in the destination mesh grid, the \c Pixel_Regrid_Type
 * object contains a list of overlapping scan grid (source) pixels,
 * and the corresponding overlap area.  Given the \c Pixel_Regrid_Type
 * structure, regridding can be performed using the \c Pixel_regrid
 * and \c Pixel_regrid_from_mesh functions.
 */
extern Pixel_Regrid_Type *
scan_init_regrid (const Scan_Type *st, const Pixel_Grid_Param_Type *mesh);

/** Write tropospheric and stratospheric columns to the corresponding granule files
 * @param[in] st   Pointer to \c Scan_Type object initialized by \c scan_read_grids
 * @param[in] fill_value   Fill value to use on output for replacing array values
 *                         containing either NaN or Inf.
 * @param[in] vtrop        Pointer to an array of tropospheric vertical column
 *                         values matching the scan grid dimensions
 * @param[in] vstrat       Pointer to an array of stratospheric vertical column
 *                         values matching the scan grid dimensions
 * @return 0 on success, -1 on error
 */
extern int scan_write_split (const Scan_Type *st, double fill_value,
                             const double *vtrop, const double *vstrat);

/** Allocate storage for required input scan-gridded product variables.
 * @param[in] num_steps   Number of mirror steps in the full scan
 * @param[in] num_xtrack  Number of cross-track pixels in the scan
 * @return On success, a pointer to an allocated \c Scan_Vars_Type object.
 */
extern Scan_Vars_Type *scan_vars_alloc (int num_steps, int num_xtrack);

/** Free resources associated with a \c Scan_Vars_Type object
 * @param[in]  sv   Pointer to an object of type \c Scan_Vars_Type allocated by
 *                  \c scan_vars_alloc
 */
extern void scan_vars_free (Scan_Vars_Type *sv);

/** Pack selected scan-gridded product variables into pre-allocated storage
 * @param[in]  st   Pointer to an object of type \c Scan_Type initialized by
 *                  \c scan_read_grids.
 * @param[in]  sv   Pointer to an object of type \c Scan_Vars_Type allocated
 *                  by \c scan_vars_alloc
 * @return 0 on success, -1 on error
 */
extern int scan_vars_pack (const Scan_Type *st, Scan_Vars_Type *sv);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
