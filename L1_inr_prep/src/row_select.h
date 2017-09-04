#ifndef __L0_INR_PREP_ROW_SELECT__
#define __L0_INR_PREP_ROW_SELECT__ 1

typedef struct Row_Select_Type Row_Select_Type;

struct Row_Select_Type
{
   Row_Select_Type *next;
   char *file;
   double *times;
   int num_times;
   int start;
   int count;
};

extern void row_select_free (Row_Select_Type *);

extern int row_select_scan (double time_beg, double time_end, int num_pad, 
                            const char *file_glob_pattern,
                            Row_Select_Type **rstp);

#endif
