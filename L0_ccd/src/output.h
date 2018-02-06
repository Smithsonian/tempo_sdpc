#ifndef __TEMPO_OUTPUT_H__
#define __TEMPO_OUTPUT_H__ 1
/** @file output.h
 *  @brief Interface to facilitate Level 1 granule output
 */

#include <libconfig.h>
#include "granule.h"
#include "sensorcal.h"

/** @brief Struct to organize the main output objects
 *
 * Along with each calibrated radiance or irradiance exposure record,
 * the output includes the corresponding radiance/irradiance uncertainty,
 * and a calibrated wavelength scale.
 */
typedef struct
{
   Spectral_Data_Type *uv;
   Spectral_Data_Type *vis;
}
Output_Exprec_Type;

typedef struct Output_Type Output_Type;

/** @brief Struct providing functions for selected output operations.
 *
 * This object supports relatively low-level output operations to
 * facilitate interleaving output with processing.  This supports
 * processing radiance granules as a stream of exposure records, without
 * requiring that the entire granule be loaded into memory at once.
 */
struct Output_Type
{
   /** set the output file path
    * @param out  non-NULL pointer to an Output_Type object
    * @param file  path to the output file
    * @return 0 on success, -1 on error
    */
   int (*out_set_file)(Output_Type *, const char *);

   /** set the spatial and wavelength dimensions in the output file
    * @param out  non-NULL pointer to an Output_Type object
    * @param num_recs   number of exposure records
    * @param num_xtrack  number of cross-track (spatial) pixels
    * @param num_waves   number of wavelengths
    * @return 0 on success, -1 on error
    */
   int (*out_set_dims)(Output_Type *, int, int, int);

   /** create the output file
    * @param out  non-NULL pointer to an Output_Type object
    * @return 0 on success, -1 on error
    *
    * The output file path and the associated granule dimensions
    * must be specified before creating the output file.
    */
   int (*out_create)(Output_Type *);

   /** query whether or not the formatted output file exists
    * @param out  non-NULL pointer to an Output_Type object
    * @return non-zero if the file exists, 0 otherwise.
    */
   int (*out_file_exists)(const Output_Type *);

   /** write an exposure record to a specific place in the output file
    * @param out  non-NULL pointer to an Output_Type object
    * @param index integer index indicating which exposure record is being written
    * @param rec   non-NULL pointer to an Output_Exprec_Type object
    * @return 0 on success, non-zero on error
    */
   int (*out_write_rec)(Output_Type *, int, const Output_Exprec_Type *);

   int (*out_copy_metadata)(Output_Type *, int);

   /** close the file associated with an Output_Type object
    * @param out   non-NULL pointer to an Output_Type object
    * @return 0 on success, non-zero on error
    */
   int (*out_close)(Output_Type *);

   /** free resources associated with an Output_Type object
    * @param out   non-NULL pointer to an Output_Type object
    */
   void (*out_free)(Output_Type *);

#ifdef OUTPUT_PRIVATE_DATA
   OUTPUT_PRIVATE_DATA
#endif
};

/** Allocate an Output_Type object
 * @param cfg  non-NULL pointer to a config_t object associated with an open configuration file
 * @param exposure_type  \ref granule_exprec_types "exposure record type" to be written out
 */
extern Output_Type *output_alloc (config_t *cfg, int exposure_type);

#endif
