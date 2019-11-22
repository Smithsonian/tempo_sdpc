#ifndef __L2_SPLIT_PROCESS_H__
#define __L2_SPLIT_PROCESS_H__ 1

#include <libconfig.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** @file process.h
 * @brief Process Level 2 data product files for granules
 *        from a single TEMPO scan
 */

/** Process Level 2 data product files for granules from a single TEMPO scan
 * @param[in] cfg   pointer to \c config_t object associated with
 *                  \c L2_split parameter file
 * @param[in] num_files  Number of granule files
 * @param[in] files      Pointer to an array of granule file names
 * @return 0 on success, -1 on error
 */
extern int process_files (config_t *cfg, int num_files, char **files);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
