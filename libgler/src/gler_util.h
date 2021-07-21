#ifndef __GLER_UTIL_INCLUDE_H__
#define __GLER_UTIL_INCLUDE_H__ 1

typedef struct GLER_Type GLER_Type;

extern GLER_Type *gler_open (int iwave, const char *config_file);
extern void gler_close (GLER_Type *gt);

extern int gler_land_lookup (const GLER_Type *gt, double taix, int *a, int *b, double *awt);
extern int gler_land_file (const GLER_Type *gt, int k, char *path, int pathlen);

extern int gler_ocean_lookup (const GLER_Type *gt, double taix, int *a, int *b, double *awt);
extern int gler_ocean_file (const GLER_Type *gt, int k, char *path, int pathlen);

#endif
