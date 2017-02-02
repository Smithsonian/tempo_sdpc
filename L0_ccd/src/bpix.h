#ifndef __BPIX_INCLUDE__
#define __BPIX_INCLUDE__ 1
/** @file bpix.h
 *  @brief Support I/O and modification of bad pixel maps.
 *
 * A bad pixel map is an array of integer bit fields.
 * Each array element represents the corresponding pixel in
 * a data image. Each bit indicates whether or not a given
 * \ref pixel_quality_flag_bits "condition" is true for that pixel.
 */

/** bit field data type */
typedef unsigned short Badpix_Bitmap_Type;

/** bit field data type used for netCDF file storage */
#define BADPIX_BITMAP_TIO_TYPE TIO_USHORT

/** @brief Struct containing a bad pixel map */
typedef struct
{
   /** Pointer to \ref pixel_quality_flag_bits "bitfield" array. */
   Badpix_Bitmap_Type *bits;
   int num_rows;                /**< size of slowest varying array dimension */
   int num_cols;                /**< size of fastest varying array dimension */
}
Badpix_Map_Type;

/** Free a Badpix_Map_Type struct
 * @param a  a pointer to a Badpix_Map_Type struct
 */
extern void bpix_free (Badpix_Map_Type *a);

/** Allocate a new Badpix_Map_Type struct
 * @param num_rows   Size of the slowest varying array dimension
 * @param num_cols   Size of the fastest varying array dimension
 * @return a non-NULL pointer to a new Badpix_Map_Type struct on success,
 * or NULL on error.
 */

extern Badpix_Map_Type *bpix_new (int num_rows, int num_cols);

/** Read a bad pixel map from a file
 * @param  file   Path to a file containing a bad pixel map.
 * @return a non-NULL pointer to a new Badpix_Map_Type struct on success,
 * or NULL on error.
 */
extern Badpix_Map_Type *bpix_read (const char *file);

/** Compute the logical OR of a bad pixel map and a bitfield array
 * @param a      non-NULL pointer to a bad pixel map.
 * @param bits   non-NULL pointer to a bitfield array
 * @return 0 on success, non-zero on error
 *
 * The bitfield array is assumed to have the same dimensions as the
 * Badpix_Map_Type object.  Each pixel value of the Badpix_Map_Type
 * object is replaced by the result of a logical OR with the corresponding
 * pixel of the bitfield array.
 */
extern int bpix_logical_or (Badpix_Map_Type *a,
                            const Badpix_Bitmap_Type *bits);

/** Opaque pointer data type used to count occurrences of selected types of bad pixels.
 *
 * A new bad pixel may be added to the bad pixel map if selected criteria
 * are observed to occur in that pixel with sufficient frequency.
 */
typedef struct Badpix_Map_Occur_Type Badpix_Map_Occur_Type;

/** Free resources allocated by bpix_occur_open
 * @param ot  Pointer to Badpix_Map_Occur_Type allocated by bpix_occur_open
 */
extern void bpix_occur_close (Badpix_Map_Occur_Type *ot);

/** Create Badpix_Map_Occur_Type object
 * @param  num_rows   Size of slowest varying array dimension
 * @param  num_cols   Size of fastest varying array dimension
 * @param  mask       Mask indicating which bit occurrences are to be counted
 * @return a non-NULL Badpix_Map_Occur_Type pointer on success, NULL on error.
 *
 * Create a Badpix_Map_Occur_Type object that counts occurrences of specified
 * bits in images with the specified array dimensions.  The non-zero bits of
 * \a mask indicate which bit occurrences are to be counted.
 */
extern Badpix_Map_Occur_Type *bpix_occur_open (int num_rows, int num_cols,
                                               Badpix_Bitmap_Type mask);

/** Increment Badpix_Map_Occur_Type counters
 * @param ot   Pointer to Badpix_Map_Occur_Type allocated by bpix_occur_open
 * @param bits Pointer to a bitfield array.
 * @return 0 on success, non-zero on error
 *
 * The bitfield array, \a bits, is assumed to have array dimensions
 * that match the Badpix_Map_Occur_Type object.  For each array
 * element, \a B, in \a bits, the corresponding counter in the
 * Badpix_Map_Occur_Type object is incremented when B|mask
 * is non-zero, where \a mask is the value specified when the
 * Badpix_Map_Occur_Type was created.
 */
extern int bpix_occur_incr (Badpix_Map_Occur_Type *ot,
                            const Badpix_Bitmap_Type *bits);

/** Set bits in a bad pixel map based on frequency of occurrence
 * @param  ot  Pointer to Badpix_Map_Occur_Type allocated by bpix_occur_open
 * @param num_threshold  Threshold number of occurrences required to set a bit
 *                       in a pixel of a bad pixel map.
 * @param mask  Mask indicating which bit occurrences may be set.
 * @param bits  Target bad pixel map in which bits may be set in qualifying pixels.
 * @return 0 on success, non-zero on error
 *
 * Whenever a bit in the subset specified by \a mask has been observed to occur
 * in a single pixel more than \a num_threshold times, then the corresponding
 * bit of that pixel may be set in the bad pixel map, \a bits.
 */
extern int bpix_occur_set (const Badpix_Map_Occur_Type *ot,
                           int num_threshold, Badpix_Bitmap_Type mask,
                           Badpix_Bitmap_Type *bits);
#endif
