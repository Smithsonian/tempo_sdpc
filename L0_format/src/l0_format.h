/** @file l0_format.h
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Interface shared by modules within L0_format
 */

#ifndef _L0_METHODS_H_
#define _L0_METHODS_H_ 1

#include <libconfig.h>
#include <iocsdpc.h>
#include "tpinfo.h"
#include "enum.h"

/* IOC data types */

#define IOCDB_DTYPE_I8  1
#define IOCDB_DTYPE_U8  2
#define IOCDB_DTYPE_I16 3
#define IOCDB_DTYPE_U16 4
#define IOCDB_DTYPE_I32 5
#define IOCDB_DTYPE_U32 6
#define IOCDB_DTYPE_I64 7
#define IOCDB_DTYPE_U64 8
#define IOCDB_DTYPE_F32 9
#define IOCDB_DTYPE_F64 10
#define IOCDB_DTYPE_ENUM 11            /* AKA "STATE", uint32 in FS */

#ifndef REALLOC
#define REALLOC realloc
#endif

#ifndef MALLOC
#define MALLOC malloc
#endif

#ifndef FREE
#define FREE free
#endif

#define MAX_BASENAME_SIZE 72
#define MAX_ISOTIME_LEN 32

#define L0_DEFAULT_OUTFILE_DURATION_SEC 300.0

/** Write a global time-stamp attribute to an open netCDF file
 * @param[in] ncid  Index of a netCDF file, opened for writing.
 * @param[in] tstamp_name  Name of time-stamp variable, normally one either
 *              @c time_coverage_start or @c time_coverage_end
 * @param[in] tstamp_value  The time-stamp value, expressed as the number of
 *                  seconds elapsed since the TEMPO epoch.
 * @return 0 on success, -1 on error
 */
extern int write_attr_global_timestamp (int ncid, const char *tstamp_name,
                                        double tstamp_value);

/** Generate a Level 0 file basename given the necessary file name components
 * @param[in]  sec_since_epoch   A time stamp value, expressed as the number
 *                               of seconds elapsed since the TEMPO epoch.
 * @param[in] processing_version Integer processing version number.
 * @param[in] suffix      File name suffix indicating the file type.
 * @param[out]  buf     Pointer to a buffer to hold the generated file name.
 * @param[in]  bufsize   Size of the file name buffer.
 *
 * @return 0 on success, -1 on error
 */
extern int make_level0_basename (double sec_since_epoch, int processing_version,
                                 const char *suffix, char *buf, int bufsize);

/** Create a netCDF file with a hidden name.
 *  @param[in] dirname   Path to the target directory that will contain the file.
 *  @param[in] basename  The file basename that will be prefixed by a "." to hide it.
 *  @param[out] ncid     netCDF file index of the opened file
 *  @return 0 on success, -1 on error
 *
 * @remark Standard metadata (e.g. @c time_reference) may also be written to the file.
 */
extern int create_hidden (const char *dirname, const char *basename, int *ncid);

/** Close a hidden file, removing a "." prefix from the name.
 *  @param[in] ncid      netCDF file index of the hidden file.
 *  @param[in] dirname   Path to the target directory that contains the file.
 *  @param[in] basename  The file basename that will be prefixed by a "." to construct
 *                       the hidden file name.
 *  @return 0 on success, -1 on error
 *
 * On return, a "." prefix will be removed from the name of the file on disk.
 */
extern int close_hidden (int ncid, const char *dirname, const char *basename);

/** Write comment and units attributes to a specified netCDF file variable */
extern int annotate_var (int grp, int varid, const char *descr, const char *units);

typedef struct Process_Method_Type Process_Method_Type;
struct Process_Method_Type
{
   int (*process)(Process_Method_Type *, const TPInfo_Type *, const char *);
   void (*delete)(Process_Method_Type *);
   int (*flush_cache)(Process_Method_Type *, const TPInfo_Type *);

#ifdef PROCESS_METHOD_PRIVATE_DATA
   PROCESS_METHOD_PRIVATE_DATA
#endif
};

/** Initialize the module that handles exposure record files
 * @param[in] cfg  Pointer to fully initialized @c config_t struct
 * @return A @c Process_Method_Type structure on success, or @c NULL on error
 */
extern Process_Method_Type *init_exprec_method (config_t *cfg);

/** Initialize the module that handles telemetry point section files
 * @param[in] cfg  Pointer to fully initialized @c config_t struct
 * @return A @c Process_Method_Type structure on success, or @c NULL on error
 */
extern Process_Method_Type *init_tpsec_method (config_t *cfg);

/** Initialize the module that handles scan mechanism controller files
 * @param[in] cfg  Pointer to fully initialized @c config_t struct
 * @return A @c Process_Method_Type structure on success, or @c NULL on error
 */
extern Process_Method_Type *init_smc_method (config_t *cfg);

/** Initialize the module that handles inertial reference unit files
 * @param[in] cfg  Pointer to fully initialized @c config_t struct
 * @return A @c Process_Method_Type structure on success, or @c NULL on error
 */
extern Process_Method_Type *init_iru_method (config_t *cfg);

#endif
