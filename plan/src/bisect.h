#ifndef __BISECT_H__
#define __BISECT_H__ 1

/** @file bisect.h
 *  @brief Bisection method root solver.
 *
 * Use bisection method to find a zero of a function, f(x),
 * with a <= x <= b.
 */

/** Find a zero of a function, f(x), with a <= x <= b.
 * @param[in]  func  A pointer to the function
 * @param[in]  a     The lower limit of the X coordinate, a <= x
 * @param[in]  b     The upper limit of the X coordinate, x <= b
 * @param[in]  cd    A pointer to client data required by the function
 *                   of interest (NULL is ok)
 * @param[out] xp    The coordinate for which f(x) == 0
 * @return 0 on success, -1 on error
 *
 * The target function interface is as follows:
 * \verbatim
 *     status = func (x, &f, client_data)
 * \endverbatim
 * where:
 *   \li \c x is the independent coordinate
 *   \li \c f is the returned function value
 *   \li \c client_data is a pointer to additional data needed by the function
 *   \li \c status is 0 on success, or -1 to indicate that an error
 *          occurred during the function evaluation.
*/
extern int bisection (int (*func)(double, double *, void *),
                      double a, double b, void *cd, double *xp);

#endif
