#ifndef __INCLUDE_LPS_APPLY_H__
#define __INCLUDE_LPS_APPLY_H__  1

#include "lps.h"

/** Use linear polarization sensitivity to convert synthetic (true) radiance to measured radiance
 * @param[in]  lps       Pointer to \a Lps_Type object allocated by \a lps_open
 * @param[in]  rad_file  Path to radiance file
 * @return 0 on success, -1 on error
 *
 * Radiances in the target file are modified in place.
 */
extern int lps_apply (Lps_Type *lps, const char *rad_file);

#endif
