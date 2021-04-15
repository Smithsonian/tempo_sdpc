#ifndef __CURRENT_INCLUDE_H__
#define __CURRENT_INCLUDE_H__ 1

#include "_process.h"

extern int current_create_file_of_type (Granule_Type *gr, const char *output_file,
                                        int num_times, int num_rows, int num_cols,
                                        int *pncid, int *pgrp);

extern int current_write_exprec (int ncid, const Exprec_Meta_Type *xr);

extern int current_image_read (int ncid, int k, Image_Type *img);

extern int current_write_mean_dark_current (int ncid);

#endif
