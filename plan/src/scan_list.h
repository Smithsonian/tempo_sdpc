#ifndef __PLAN_SCAN_LIST_H__
#define __PLAN_SCAN_LIST_H__ 1

/** @file scan_list.h
 *  @brief Manage a list of instrument scan parameters
 */

#include <stdio.h>

typedef struct Scan_List_Entry Scan_List_Entry;
struct Scan_List_Entry
{
   Scan_List_Entry *next;
   double tstart;         /* UTC [days] */
   double xstart;         /* microradian */
   double scan_duration;  /* sec */
   double step_exposure;  /* sec */
   int num_steps;
   int num_repeats;
};

extern Scan_List_Entry *scan_list_entry_alloc (void);
extern void scan_list_entry_free (Scan_List_Entry *stt);
extern int scan_list_append (Scan_List_Entry **head,
                             Scan_List_Entry *stt);
extern void scan_list_free (Scan_List_Entry *head);

extern int scan_list_write (FILE *fp, const Scan_List_Entry *head);

#endif
