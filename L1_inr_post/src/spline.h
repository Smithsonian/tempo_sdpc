#ifndef __INR_POST_SPLINE_H__
#define __INR_POST_SPLINE_H__ 1

/** @file spline.h
 *  @brief Cubic spline interface
 */

extern int spline (const double *x, const double *y, int nx,
                   const double *xs, double *ys, int nxs);

#endif
