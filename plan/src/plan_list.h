#ifndef __PLAN_PLAN_LIST_H__
#define __PLAN_PLAN_LIST_H__ 1

/** @file plan_list.h
 *  @brief Manage a list of planned instrument scans
 */

#include <stdio.h>

typedef struct Plan_List_Type Plan_List_Type;
struct Plan_List_Type
{
   Plan_List_Type *next;
   double tstart;            /**< scan start time, UTC [days] */
   double xstart;            /**< scan start mirror coordinate [microradian] */
   double scan_duration;     /**< scan duration [sec] */
   double integration_time;  /**< integration time for a single exposure in a co-add [sec] */
   int num_steps;            /**< number of mirror steps in the scan */
   int num_repeats;          /**< number of scans on this day */
};

/** Allocate a \ref Plan_List_Type structure
 * @return  A \ref Plan_List_Type pointer on success, NULL on error
 *
 * When no longer needed, the returned structure should be freed
 * by a call to \ref plan_list_entry_free
*/
extern Plan_List_Type *plan_list_entry_alloc (void);

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
 * @param[in]  fp    Initialized FILE pointer for the destination file
 * @param[in]  head  The head of a plan list.
 * @return 0 on success, -1 on error.
*/
extern int plan_list_write (FILE *fp, const Plan_List_Type *head);

#endif
