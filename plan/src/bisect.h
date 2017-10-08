#ifndef __BISECT_H__
#define __BISECT_H__ 1

extern int bisection (int (*func)(double, double *, void *),
                      double a, double b, void *cd, double *xp);

#endif
