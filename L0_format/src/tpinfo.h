/** @file tpinfo.h
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Interface for retrieving telemetry point definitions
 */

#ifndef __TPINFO_INCLUDE_H__
#define __TPINFO_INCLUDE_H__ 1

typedef struct TPInfo_Type TPInfo_Type;

typedef struct
{
   char *mnemonic;
   char *description;
   char *units;
   char *enumlist;
}
TPFields_Type;

/** Parse CSV file containing telemetry point definitions
 *
 * @param[in] file  Path to the telemetry point CSV file
 * @return Pointer to a @c TPInfo_Type struct, or @c NULL in case of error.
 *
 * When no longer needed, the @c TPInfo_Type struct should be freed
 * using @c tpinfo_free
 *
 * @see @c tpinfo_free @c tpinfo_find
 */
extern TPInfo_Type *tpinfo_init (const char *file);

/** Free memory allocated by @c tpinfo_init
 *
 * @param[in] tp  Pointer to @c TPInfo_Type struct allocated by @c tpinfo_init
 */
extern void tpinfo_free (TPInfo_Type *tp);

/** Find information on a telemetry point specified by name
 *
 * @param[in] tp  Pointer to @c TPInfo_Type struct allocated by tpinfo_init
 * @param[in] mnemonic  Name of the telemetry point
 * @param[in/out] fields  Pointer to @c TPFields_Type struct containing pointers
 *                     into the @c TPInfo_Type struct.  Do not modify any of
 *                     the objects made available by the @c TPFields_Type
 *                     struct fields.
 * @return A @c TPFields_Type pointer on success, or @c NULL on error.
 *
 * @remark When the @p fields parameter is non-NULL, that @c TPFields_Type
 * struct is used to store the results, and a pointer to that struct is returned.
 * When the @p fields parameter is NULL, a @c TPFields_Type struct is allocated
 * and returned.  When no longer needed, an allocated @c TPFields_Type struct should
 * be freed by calling @c tpinfo_free_tpfields
 *
 * @see @c tpinfo_free_tpfields
 */
extern TPFields_Type *tpinfo_find (const TPInfo_Type *tp, const char *mnemonic, TPFields_Type *fields);

/** Free memory allocated by @c tpinfo_find
 *
 * @param[in] tp  Pointer to @c TPFields_Type struct allocated by @c tpinfo_find
 */
extern void tpinfo_free_tpfields (TPFields_Type *info);

#endif
