/** @file enum.h
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Interface for managing enumerated data types
 */

#ifndef __ENUM_LOOKUP_H__
#define __ENUM_LOOKUP_H__ 1

typedef struct Enum_Lookup_Type Enum_Lookup_Type;

/** Allocate an Enum_Lookup_Type struct
 * @param[in] size  The size of the lookup table to construct
 *                  [should be at least 20% larger than the number
 *                   of enum definitions it is expected to hold]
 * @return  An Enum_Lookup_Type struct on success, or NULL on error
 *
 * The memory allocated by @c elt_open should be freed by @c elt_close
 *
 * @see elt_close
 */
extern Enum_Lookup_Type *elt_open (size_t size);

/** Free an Enum_Lookup_Type struct allocated by elt_open
 * @param[in] elt Enum_Lookup_Type struct allocated by elt_open
 */
extern void elt_close (Enum_Lookup_Type *elt);

/** Determine whether an @c enumlist defines a valid enum
 * @param[in] elt  Enum_Lookup_Type struct allocated by elt_open
 * @param[in] enumlist  String of the form "name0:integer0, name1:integer1, ..."
 *                  providing a list of valid name-integer pairs defining
 *                  an enum.
 * A valid enum must have at least two allowed values.
 * @return non-zero for a valid enum, 0 otherwise.
 */
extern int elt_is_valid (Enum_Lookup_Type *elt, const char *enumlist);

/** Define an enum as a datatype in a netCDF file.
 * @param[in] elt Enum_Lookup_Type struct allocated by @c elt_open
 * @param[in] ncid  netCDF file index of a netCDF file open for writing
 * @param[in] mnemonic  Name of the telemetry point mnemonic that has
 *                      the associated enum type.
 * @param[in] enumlist  String of the form "name0:integer0, name1:integer1, ..."
 *                      defining the allowed values of the enum data type.
 * @param[in] base_type  netCDF type index of the enum's base integer type.
 * @param[out] type_id   The type index of the enum type in the netCDF file.
 *
 * An enum is defined by its enumlist string.
 * When the netCDF file already contains an equivalent enum, the @p type_id
 * of the existing enum is returned and no new enum is defined.
 * Otherwise, a new enum datatype is written to the file, and the @p type_id
 * of that enum datatype is returned.
 */
extern int elt_define (Enum_Lookup_Type *elt, int ncid, const char *mnemonic,
                       const char *enumlist, int base_type, int *type_id);

#endif

