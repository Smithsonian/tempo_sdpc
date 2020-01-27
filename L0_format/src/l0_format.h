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
#define MAX_PATHLEN 1024

#define L0_DEFAULT_OUTFILE_DURATION_SEC 300.0

/** Remove a file
 * @param[in]  dirname   Directory path
 * @param[in]  basename  File name
 * @return 0 on success, -1 on error
 */
extern int remove_file (const char *dirname, const char *basename);

/** Perform shell variable expansion on a string
 * @param[in]  s  Input string
 * @return on success, an allocated string, NULL on error.
 *
 * Use \a free to free the allocated string.
 */
extern char *expand_string (const char *s);

/** Retrieve the processing version number */
extern int get_processing_version (void);

/** Set the TEMPO epoch, and make sure all files use the same epoch
 * @param[in] epoch   TEMPO epoch specification, expressed as the number
 *                    of elapsed seconds (UTC) since the Unix epoch
 * @return 0 on success, -1 on error
 */
extern int verify_epoch (time_t epoch);

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

/** Write standard global metadata to an open netCDF file
 * @param[in] ncid  Index of a netCDF file, opened for writing.
 * @param[in] chdr  Pointer to IOCSDPC common header
 * @return 0 on success, -1 on error
 */
extern int write_std_global_metadata (int ncid, const IOCSDPC_Common_Header_Type *chdr);

typedef struct
{
   int scan_num;
   int scan_type;
   int granule_num;
   int granule_flag;
}
Radiance_Ident_Type;

/** Generate a Level 0 file basename given the necessary file name components
 * @param[out] buf       Pointer to a buffer to hold the generated file name.
 * @param[in]  bufsize   Size of the file name buffer.
 * @param[in]  sec_since_epoch   A time stamp value, expressed as the number
 *                               of seconds elapsed since the TEMPO epoch.
 * @param[in] processing_version Integer processing version number.
 * @param[in] suffix      File name suffix indicating the file type.
 *
 * @return 0 on success, -1 on error
 */
extern int make_level0_basename (char *buf, int bufsize,
                                 double sec_since_epoch, int processing_version,
                                 const char *suffix, const Radiance_Ident_Type *identp);

/** Generate the archive directory path for a Level 0 file
 * @param[inout] archdir_path      The full archive directory path (malloced)
 * @param[in]  sec_since_epoch   A time stamp value, expressed as the number
 *                               of seconds elapsed since the TEMPO epoch.
 * @param[in]  scan_num          Optional scan number (<0 means do not use this in the path)
 * @param[in] suffix      File name suffix indicating the file type.
 *
 * @return 0 on success, -1 on error
 *
 * If \a archdir_path is \a NULL, this function does nothing and returns 0.
 */
extern int make_level0_archdir_path (char **archdir_path, double sec_since_epoch, int scan_num,
                                     const char *suffix);

/** Create a netCDF file with a hidden name.
 *  @param[in] dirname   Path to the target directory that will contain the file.
 *  @param[in] basename  The file basename that will be prefixed by a "." to hide it.
 *  @param[out] ncid     netCDF file index of the opened file
 *  @return 0 on success, -1 on error
 *
 * @remark Standard metadata (e.g. @c time_reference) may also be written to the file.
 */
extern int create_hidden (const char *dirname, const char *basename, int *ncid);

/** Close a hidden file, removing a "." prefix from the name, and optionally put a copy into a specified directory
 *  @param[in] ncid      netCDF file index of the hidden file.
 *  @param[in] dirname   Path to the target directory that contains the file.
 *  @param[in] basename  The file basename that will be prefixed by a "." to construct
 *                       the hidden file name.
 *  @param[in] copydir   Directory to receive a non-hidden copy of the file before the rename occurs.
 *                       If \a copydir is NULL, no copy is created.
 *  @return 0 on success, -1 on error
 *
 * On return, a "." prefix will be removed from the name of the file on disk.
 */
extern int close_hidden (int ncid, const char *dirname, const char *basename, const char *copydir);

/** Write comment and units attributes to a specified netCDF file variable */
extern int annotate_var (int grp, int varid, const char *descr, const char *units);

typedef struct Process_Method_Type Process_Method_Type;
typedef int Process_Method_Callback_Function (Process_Method_Type *, void *);

struct Process_Method_Type
{
   int (*pmt_process)(Process_Method_Type *, const TPInfo_Type *, const char *, void *);
   void (*pmt_delete)(Process_Method_Type *);
   int (*pmt_flush_cache)(Process_Method_Type *, const TPInfo_Type *);
   int (*pmt_query_latest_timestamp)(Process_Method_Type *, int, double *);
   int (*pmt_query_when_last_rec_cached)(const Process_Method_Type *, time_t *);

   Process_Method_Callback_Function *pmt_post_process_callback;

#ifdef PROCESS_METHOD_PRIVATE_DATA
   PROCESS_METHOD_PRIVATE_DATA
#endif
};

/** Initialize the module that handles exposure record files
 * @param[in] cfg  Pointer to fully initialized @c config_t struct
 * @return A @c Process_Method_Type structure on success, or @c NULL on error
 *
 * Note: The caching method must be set before the \a Process_Method_Type
 * is initialized.
 */
extern Process_Method_Type *init_exprec_method (config_t *cfg);

/** Set the exposure record caching method
 * @param[in]  method  Integer method identifier
 */
extern void set_exprec_cache_method (int method);
enum
{
   EXPREC_CACHE_DISK = 0,
   EXPREC_CACHE_MEM  = 1
};

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
