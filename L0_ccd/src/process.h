#ifndef __PROCESS_INCLUDE_H__
#define __PROCESS_INCLUDE_H__ 1

#include <libconfig.h>

#include "control.h"

extern int process_inputs (config_t *cfg, const Control_Type *ctrl);
extern void process_set_version (int version);
extern int process_get_version (void);

#endif
