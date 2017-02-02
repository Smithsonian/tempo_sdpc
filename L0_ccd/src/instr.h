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
 * @param  file   Path to a file
 * @return non-NULL pointer to an Instr_Type object on success, NULL on error
 *
 * The specified file may be either a single netCDF file of the correct format,
 * or may be an ascii file containing one or more netCDF file paths.  When an
 * ascii file is provided, the first character of the \a file string should
 * be the '@' symbol.
 */
extern Instr_Type *instr_open (const char *file);

#endif
