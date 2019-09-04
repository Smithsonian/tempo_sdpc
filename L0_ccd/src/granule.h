#ifndef __GRANULE_INCLUDE__
#define __GRANULE_INCLUDE__ 1
/** @file granule.h
 *  @brief Interface to exposure record granules
 *
 * This module provides functions to access exposure record data
 * in dark, irradiance and radiance granules.
 */

#include "image.h"

/** @brief Exposure record types
 *  @anchor granule_exprec_types
 *
 * These enum values should match the corresponding IOCSDPC definitions
 */
enum
{
   EXPREC_TYPE_RAD,        /**< radiance */
   EXPREC_TYPE_DARK,       /**< dark */
   EXPREC_TYPE_IRR_WRK,    /**< irradiance [working diffuser] */
   EXPREC_TYPE_IRR_REF,    /**< irradiance [reference diffuser] */
   EXPREC_TYPE_LIN_IRR,    /**< irradiance - linearity */
   EXPREC_TYPE_LIN_DARK,   /**< dark - linearity */
   EXPREC_TYPE_UNKNOWN = -1   /**< unknown */
};

#define EXPREC_TYPE_IS_LINEARITY(type) \
   (((type) == EXPREC_TYPE_LIN_DARK) || ((type) == EXPREC_TYPE_LIN_IRR))

#define EXPREC_TYPE_IS_DARK(type) \
   (((type) == EXPREC_TYPE_DARK) || ((type) == EXPREC_TYPE_LIN_DARK))

#define EXPREC_TYPE_IS_IRRADIANCE(type) \
   (((type) == EXPREC_TYPE_IRR_WRK) || ((type) == EXPREC_TYPE_IRR_REF) || ((type) == EXPREC_TYPE_LIN_IRR))

/** @brief Granule exposure record object
 *
 * The content of this struct closely reflects the content of the exposure
 * records provided by the IOC.
 */
typedef struct
{
   double start_time;           /**< exposure start time in sec elapsed since the TEMPO epoch */
   double exposure_time;        /**< total exposure duration [sec] */
   double frame_transfer_time;  /**< frame transfer time [sec] */
   double readout_time;         /**< storage region readout time [sec] */
   int exposure_type;           /**< \ref granule_exprec_types "exposure record type" */
   int num_coadds;              /**< number of co-adds */
   int curr_mirror_step;        /**< current mirror step */
   /* ConOps 3.3: instrument command parameters: NUM_TG_ROWS, NUM_DG_ROWS */
   int num_dg_rows;             /**< index of first row included in storage region dark current sum [rows numbered 1..N] */
   int num_tg_rows;             /**< number of rows included in storage region dark current sum */
   Image_Type *img;             /**< Pointer to image data */
}
Granule_Exprec_Type;

typedef struct Granule_Type Granule_Type;

/** @brief Struct providing functions to access exposure record granules */
struct Granule_Type
{
   /** Free a Granule_Type object
    * @param  Granule_Type * non-NULL pointer to a Granule_Type object
    */
   void (*granule_close) (Granule_Type *);

   /** Query the number of exposure records in a granule
    * @param  Granule_Type * non-NULL pointer to a Granule_Type object
    * @return the number of exposure records on success, -1 on error
    */
   int (*granule_num_exprecs)(const Granule_Type *);

   /** Query the type of exposure records in a granule
    * @param  Granule_Type * non-NULL pointer to a Granule_Type object
    * @param  int *          pointer to an integer that will hold the exposure record type
    * @return 0 on success, -1 on error
    */
   int (*granule_type)(const Granule_Type *, int *);

   /** Read a specified exposure record from a granule
    * @param  Granule_Type * non-NULL pointer to a Granule_Type object
    * @param  int            integer index of the exposure record, 0 <= i < N
    * @param  Granule_Exprec_Type **  (optional) pointer to an existing Granule_Exprec_Type object
    * @return non-NULL pointer to a Granule_Exprec_Type object on success, NULL on error
    *
    * When the Granule_Exprec_Type ** argument is non-NULL, the object it points
    * to is overwritten by the exposure record data that is read from disk.
    * When the Granule_Exprec_Type ** argument is NULL, the return value points
    * to allocated storage which must be subsequently freed by calling the
    * granule_free_exprec method. For example, the following code re-uses a single
    * Granule_Exprec_Type object allocated when i=0 and does not leak memory.
    * @code
    *     Granule_Exprec_Type *x = NULL;
    *     int i, n = gr->granule_num_exprecs (gr);
    *
    *     for (i = 0; i < n; i++)
    *       (void) gr->granule_read_exprec_by_index (gr, i, &x);
    *
    *     gr->granule_free_exprec (x);
    * @endcode
    */
   Granule_Exprec_Type *(*granule_read_exprec_by_index) (const Granule_Type *, int, Granule_Exprec_Type **);

   /** Return the netCDF file descriptor for a Granule_Type object
    * @param  Granule_Type *  non-NULL pointer to a Granule_Type object
    * @return the file descriptor on success, or -1 on error
    */
   int (*granule_ncid)(const Granule_Type *);

   /** Return the coverage start time for a Granule_Type object
    * @param  Granule_Type *  non-NULL pointer to a Granule_Type object
    * @return the coverage start time, in TAI seconds since the TEMPO epoch
    */
   double (*granule_tstart)(const Granule_Type *);

   /** Return the coverage end time for a Granule_Type object
    * @param  Granule_Type *  non-NULL pointer to a Granule_Type object
    * @return the coverage end time, in TAI seconds since the TEMPO epoch
    */
   double (*granule_tend)(const Granule_Type *);

   /** Free a Granule_Exprec_Type object
    * @param  Granule_Exprec_Type *  non-NULL pointer to a Granule_Exprec_Type object
    */
   void (*granule_free_exprec) (Granule_Exprec_Type *);

#ifdef GRANULE_TYPE_PRIVATE_DATA
   GRANULE_TYPE_PRIVATE_DATA
#endif
};

/** Open a granule file
 * @param file  Path to the granule file
 * @return non-NULL Granule_Type struct on success, NULL on error
 */
extern Granule_Type *granule_open (const char *file);

#endif
