#ifndef __PLAN_CONFIG_H__
#define __PLAN_CONFIG_H__ 1

#define MALLOC malloc
#define FREE free
#define REALLOC realloc

#define SEC_PER_DAY        86400.0

extern int mkjdtimestr (double jd_utc, char *buf, int bufsize);

#endif
