#ifndef __CONTROL_INCLUDE_H__
#define __CONTROL_INCLUDE_H__ 1

typedef struct
{
   const char *input_file;
   const char *output_file;
   const char *bpix_file;
   const char *dark_file;
   const char *instr_status_file;
   const char *instr_glob;
   int limit_num_granules;
}
Control_Type;

#endif
