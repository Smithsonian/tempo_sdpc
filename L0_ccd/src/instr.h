#ifndef __INSTR_INCLUDE_H__
#define __INSTR_INCLUDE_H__ 1

typedef struct Instr_Type Instr_Type;

struct Instr_Type
{
   void (*instr_delete)(Instr_Type *);
   int (*instr_ccd_temp1)(const Instr_Type *, double, float *);
   int (*instr_ccd_temp2)(const Instr_Type *, double, float *);

#ifdef INSTR_PRIVATE_DATA
   INSTR_PRIVATE_DATA
#endif   
};

extern Instr_Type *instr_open (const char *file);

#endif
