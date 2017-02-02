#ifndef __TEMPO_DARK_INCLUDE_H__
#define __TEMPO_DARK_INCLUDE_H__ 1
/** @file dark.h
 *  @brief Facilitate creating and using tables of dark current images.
 */

#include <libconfig.h>
#include "image.h"

/** Opaque pointer data type to manage arrays of dark current images */
typedef struct Dark_Array_Type Dark_Array_Type;

/** Free a Dark_Array_Type object
 * @param da   Pointer to a Dark_Array_Type object allocated by dark_array_alloc
 */
extern void dark_array_free (Dark_Array_Type *da);

/** Allocate a Dark_Array_Type object
 * @param num_darks  Number of dark current arrays to be stored
 * @return non-NULL pointer to a Dark_Array_Type object on success, NULL on error
 */
extern Dark_Array_Type *dark_array_alloc (int num_darks);

/** Set an element in a Dark_Array_Type object
 * @param  da  Pointer to a Dark_Array_Type object
 * @param  i   Integer index specifying which array element to set
 * @param  img  Non-NULL pointer to an Image_Type object
 * @param  sdc  Storage region dark current value from the dark image
 * @param  fp_temp  Focal plane temperature of the dark image
 * @param  exposure_time  Exposure time of the dark image
 * @return 0 on success, non-zero on error
 *
 * Each element in a Dark_Array_Type object contains a dark current image
 * plus additional metadata.
 */
extern int dark_array_elem_set (Dark_Array_Type *da, int i,
                                Image_Type *img, double sdc, double fp_temp,
                                double exposure_time);

/** Write a Dark_Array_Type to a file
 * @param dark_array  non-NULL pointer to a Dark_Array_Type object.
 * @param file        Path to output file
 */
extern int dark_array_write (const Dark_Array_Type *dark_array,
                             const char *file);

typedef struct Dark_Table_Type Dark_Table_Type;
typedef struct Dark_Config_Type Dark_Config_Type;

/** @brief Dark table ordering options */
enum
{
   DARK_TABLE_ORDERED_BY_TEMP = 1,      /**< sort by focal plane temperature */
   DARK_TABLE_ORDERED_BY_SDC = 2,       /**< sort by storage region dark current */
   DARK_TABLE_ORDERED_BY_EXPTIME = 3    /**< sort by exposure time */
};

/** @brief Dark table object
 *
 * Struct member functions facilitate basic operations on dark a current lookup table.
 * A dark current lookup table contains an array of images sorted in increasing order
 * of a particular sort key such as focal plane temperature or mean storage region dark
 * current.  Given a value of the sort key, a dark current image for that sort key value
 * can be constructed by interpolation.
 */
struct Dark_Table_Type
{
   /** delete a Dark_Table_Type object
    * @param Dark_Table_Type  non-NULL pointer to the Dark_Table_Type object
    */
   void (*dtt_delete)(Dark_Table_Type *);

   /** Query the dark table object sort order
    * @param Dark_Table_Type  non-NULL pointer to the Dark_Table_Type object
    */
   int (*dtt_ordering)(const Dark_Table_Type *);

   /** Write the dark table object to a file
    * @param Dark_Table_Type  non-NULL pointer to the Dark_Table_Type object
    * @param const char *     File path string.
    */
   int (*dtt_write)(const Dark_Table_Type *, const char *);

   /** Construct a dark current image by 1-D interpolation on the sort key value
    * @param Dark_Table_Type  non-NULL pointer to the Dark_Table_Type object
    * @param double           sort key value of interest
    * @param Image_Type       Pointer to existing Image_Type object with dimensions
    *                         appropriate to hold the interpolated dark current image
    */
   int (*dtt_interp)(const Dark_Table_Type *, double, Image_Type *);

   /** Query the range of valid sort key values for the dark current table */
   void (*dtt_domain)(const Dark_Table_Type *, double *, double *);

#ifdef DARK_TABLE_PRIVATE_DATA
   DARK_TABLE_PRIVATE_DATA
#endif
};

/** Read a dark current table from a file
 * @param file   Path to the dark current table
 * @return a non-NULL pointer to a Dark_Table_Type object on success, NULL on error
 */
extern Dark_Table_Type *dark_table_read (const char *file);

/** Free a Dark_Config_Type object
 * @param dcfg   Pointer to a Dark_Config_Type object created by dark_table_config
 */
extern void dark_table_config_free (Dark_Config_Type *dcfg);

/** Initialize a Dark_Config_Type object
 * @param cfg   Pointer to a config_t struct associated with an open configuration file
 * @param is_linearity   Flag indicating whether or not the table is to be used for
 *                       linearity calibration (non-zero indicates yes).
 * @return non-NULL pointer to a Dark_Config_Type object on success, NULL on error
 *
 * The configuration file provides the control parameters needed to construct
 * a Dark_Table_Type object.
 */
extern Dark_Config_Type *dark_table_config (config_t *cfg, int is_linearity);

/** Create a Dark_Table_Type object
 * @param dcfg  Pointer to a Dark_Config_Type object initialized by dark_table_config
 * @param dark_array  Pointer to a Dark_Array_Type object
 * @return non-NULL pointer to a Dark_Table_Type object on success, NULL on error
 *
 * The dark images in \a dark_array are sorted in the specified order and then
 * binned and averaged to create a lookup table that captures the dependence of the
 * dark current image on the specified sort key with fewer entries than
 * the input Dark_Array_Type object.
 */
extern Dark_Table_Type *dark_table_create (const Dark_Config_Type *dcfg,
                                           const Dark_Array_Type *dark_array);
#endif
