#ifndef __L0_INR_PREP_ROW_SELECT__
#define __L0_INR_PREP_ROW_SELECT__ 1

typedef struct Row_Select_Type Row_Select_Type;

struct Row_Select_Type
{
   Row_Select_Type *next;    /**< pointer to the next Row_Select_Type struct */
   char *file;               /**< name of file containing time series data */
   double *times;            /**< array of time samples [sec since TEMPO epoch] */
   int num_times;            /**< number of time samples */
   int start;                /**< array index of earliest time sample within specified time interval */
   int count;                /**< number of time samples within specified time interval */
};

/** Free memory allocated by @c row_select_scan
 *
 * @param[in] rst   Pointer to @c Row_Select_Type struct
 */
extern void row_select_free (Row_Select_Type *rst);

/** Find time series data files that cover a specified time interval.
 *
 * @param[in] time_beg    Beginning of time interval, expressed in seconds
 *                        elapsed since the TEMPO epoch
 * @param[in] time_end    End of time interval, expressed in seconds
 *                        elapsed since the TEMPO epoch
 * @param[in] num_pad     Approximate number of additional time samples to
 *                        include as padding before @c time_beg and after @c time_end
 * @param[in] num_files   Number of files to be examined.
 * @param[in] file_list   Pointer to array of strings containing file paths.
 * @param[in] group_path  Pointer to file group containing time series (ignored if NULL)
 * @param[out] rstp      Pointer to a linked list of @c Row_Select_Type objects
 *                       that cover the specified time interval with the
 *                       required padding.
 * @return 0 on success, -1 on error
 */
extern int row_select_scan (double time_beg, double time_end, int num_pad,
                            int num_files, char **file_list,
                            const char *group_path, Row_Select_Type **rstp);

#endif
