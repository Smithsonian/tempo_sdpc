#ifndef __PLAN_PLAN_LIST_H__
#define __PLAN_PLAN_LIST_H__ 1

/** @file plan_list.h
 *  @brief Manage a list of instrument scan parameters
 */

#include <stdio.h>

typedef struct Plan_List_Type Plan_List_Type;
struct Plan_List_Type
{
   Plan_List_Type *next;
   double tstart;         /* UTC [days] */
   double xstart;         /* microradian */
   double scan_duration;  /* sec */
   double step_exposure;  /* sec */
   int num_steps;
   int num_repeats;
};

extern Plan_List_Type *plan_list_entry_alloc (void);
extern void plan_list_entry_free (Plan_List_Type *ple);
extern int plan_list_append (Plan_List_Type **head,
                             Plan_List_Type *ple);
extern void plan_list_free (Plan_List_Type *head);

extern int plan_list_write (FILE *fp, const Plan_List_Type *head);

#endif
