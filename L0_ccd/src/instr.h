#ifndef __INSTR_INCLUDE_H__
#define __INSTR_INCLUDE_H__ 1
/** @file instr.h
 *  @brief Interface to instrument status data
 */

typedef struct Instr_Type Instr_Type;

/** @brief Struct providing functions to access instrument status parameters
 */
struct Instr_Type
{
   /** Free an Instr_Type object
    * @param instr  pointer to an Instr_Type object
    */
   void (*instr_delete)(Instr_Type *);

   /** Look up the value of CCD_TEMP1 at a specific time
    * @param instr   non-NULL pointer to an Instr_Type object
    * @param time    time stamp value in seconds elapsed since the TEMPO epoch
    * @param ccd_temp1   Pointer to a float scalar that will receive the value of CCD_TEMP1
    */
   int (*instr_ccd_temp1)(const Instr_Type *, double, float *);

   /** Look up the value of CCD_TEMP2 at a specific time
    * @param instr   non-NULL pointer to an Instr_Type object
    * @param time    time stamp value in seconds elapsed since the TEMPO epoch
    * @param ccd_temp2   Pointer to a float scalar that will receive the value of CCD_TEMP2
    */
   int (*instr_ccd_temp2)(const Instr_Type *, double, float *);

#ifdef INSTR_PRIVATE_DATA
   INSTR_PRIVATE_DATA
#endif
};

/** Create an Instr_Type object for a specified file
 * @param  path         Path to a file specifier
 * @param  glob_basename  [optional] globbing pattern for instrument file basename
 * @param  tstart       [optional] desired coverage start time
 * @param  tend         [optional] desired coverage end time
 * @return non-NULL pointer to an Instr_Type object on success, NULL on error
 *
 * The file specifier may be one of the following:
 *   1) The path to a single netCDF file of the correct format.
 *   2) The path to an ascii file containing one or more netCDF file paths.
 *      When an ascii file is provided, the first character of the \a path
 *      string should be the '@' symbol.
 *   3) The path to a directory containing netCDF files matching the
 *      provided basename globbing pattern.  In this case, the
 *      provide coverage start/end times are used to select the relevant
 *      files from the directory.
 */
extern Instr_Type *instr_open (const char *file, const char *glob_basename,
                               double tstart, double tend);

#endif
