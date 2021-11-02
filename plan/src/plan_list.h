#ifndef __PLAN_PLAN_LIST_H__
#define __PLAN_PLAN_LIST_H__ 1

/** @file plan_list.h
 *  @brief Manage a list of planned instrument scans
 */

#include <stdio.h>
#include <stdint.h>

typedef struct Plan_List_Type Plan_List_Type;
struct Plan_List_Type
{
   Plan_List_Type *next;
   double tstart;            /**< scan start time, UTC [days] */
   double xstart;            /**< scan start mirror coordinate [microradian] */
   double ystart;
   double scan_duration;     /**< scan duration [sec] */
   double integration_time;  /**< integration time for a single exposure in a co-add [sec] */
   double jd_utc_beg_safe;   /**< Earliest time when the aperture may safely open on this day, UTC [days] */
   double jd_utc_end_safe;   /**< Latest time when the aperture may safely open on this day, UTC [days] */
   int num_steps;            /**< number of mirror steps in the scan */
   int num_repeats;          /**< number of scans on this day */
   uint16_t scan_type;       /**< scan type value */
};

/** Allocate a \ref Plan_List_Type structure
 * @param[in]  \ref scan_type   scan type label
 * @return  A \ref Plan_List_Type pointer on success, NULL on error
 *
 * When no longer needed, the returned structure should be freed
 * by a call to \ref plan_list_entry_free
*/
extern Plan_List_Type *plan_list_entry_alloc (uint16_t scan_type);

/** Free resources associated with a \ref Plan_List_Type structure
 * @param[in]  ple  A \ref Plan_List_Type pointer allocated by
 *                     \ref plan_list_entry_alloc
*/
extern void plan_list_entry_free (Plan_List_Type *ple);

/** Append a new plan entry to a Plan_List_Type object
 * @param[in]  head  The head of a plan list.
 * @param[in]  ple   The plan entry to be appended
 * @return 0 on success, -1 on error
*/
extern int plan_list_append (Plan_List_Type **head,
                             Plan_List_Type *ple);

/** Free a list of Plan_List_Type objects
 * @param[in]  head  The head of a plan list.
*/
extern void plan_list_free (Plan_List_Type *head);

/** Write plan list parameters to an ASCII file.
 * @param[in]  fp     Initialized FILE pointer for the destination file
 * @param[in]  mirror_tilt  Pointer to optional function to compute mirror
 *                          tilt angle (NULL is ok)
 * @param[in]  head  The head of a plan list.
 * @param[in]  plan_id  The planID string given in the spacecraft maneuver table.
 * @return 0 on success, -1 on error.
*/
extern int plan_list_write (FILE *fp, double (*mirror_tilt)(double),
                            const Plan_List_Type *head, const char *plan_id);

#endif
