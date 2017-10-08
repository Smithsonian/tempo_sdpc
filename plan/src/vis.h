#ifndef __PLAN_VIS_H__
#define __PLAN_VIS_H__ 1

/** @file vis.c
 *  @brief Map solar zenith angle vs position to help visualize plan
 */

#include <libconfig.h>

typedef struct Vis_Type Vis_Type;

extern void vis_free(Vis_Type *v);
extern Vis_Type *vis_init (config_t *cfg, Solar_Geom_Type *sgt);
extern double *vis_sza (const Vis_Type *v, double jd_utc, double *);

extern int vis_write_grid (Vis_Type *v, int ncid);
extern int vis_write_value (const Vis_Type *v, int ncid, double jd_utc, 
                            const char *name, const double *value);

#endif
