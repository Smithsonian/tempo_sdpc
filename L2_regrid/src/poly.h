#ifndef __REGRID_GEOM_H__
#define __REGRID_GEOM_H__ 1

/** @file poly.h
 *  @brief Manipulate polygons, compute polygon areas,
 *       clip one polygon with another.
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct Polygon_Type Polygon_Type;
typedef struct Polygon_Clip_Type Polygon_Clip_Type;

/** Free memory associated with a Polygon_Type struct.
 * @param[in]   p   The polygon.
 *
 * p is NULL is ok.
 */
extern void Polygon_free (Polygon_Type *p);

/** Allocate a Polygon_Type struct.
 * @param[in]   n    The initial number of vertices.
 * @return A Polygon_Type struct pointer on success, NULL on error
 */
extern Polygon_Type *Polygon_new (int n);

/** Define a polygon by an ordered set of (X,Y) vertex coordinates.
 * @param[in]   p    A Polygon_Type structure obtained
 *                   by calling Polygon_new.
 * @param[in]   n    The number of vertices
 * @param[in]   x,y  The vertex coordinates.
 *
 * @return 0 on success, -1 on error
 *
 * Normally the vertices are provided in counter-clockwise order,
 * with no duplicates.
 * If necessary, the Polygon_Type structure will reallocate enough
 * space to hold the vertices provided.
 */
extern int Polygon_set (Polygon_Type *p, int n,
                        const double *x, const double *y);

/** Append a vertex point to a polygon.
 * @param[in]   p    A Polygon_Type structure obtained
 *                   by calling Polygon_new.
 * @param[in]   x,y  The vertex coordinates.
 *
 * @return 0 on success, -1 on error
 */
extern int Polygon_add (Polygon_Type *p, double x, double y);

/** Return the number of vertices in a Polygon_Type struct.
 * @param[in]   p    The polygon.
 *
 * @return The (non-negative) number of vertices on success,
 *         -1 on error
 */
extern int Polygon_length (const Polygon_Type *p);

/** Return the coordinates of vertex i.
 * @param[in]   p    The polygon.
 * @param[in]   i    The index of a vertex in the range 0<=i<n,
 *                   where n is the number of vertices.
 * @param[out]  x,y  The coordinates of the specified vertex.
 *
 * @return 0 on success, -1 on error
 */
extern int Polygon_vertex (const Polygon_Type *p, int i,
                           double *x, double *y);

/** Compute the bounding box enclosing a polygon.
 * @param[in]    p          The polygon
 * @param[out]   xmin,xmax  The bounding box X interval
 * @param[out]   ymin,ymax  The bounding box Y interval
 */
extern void Polygon_bbox (const Polygon_Type *p,
                          double *xmin, double *xmax,
                          double *ymin, double *ymax);

/** Compute the area enclosed by a polygon.
 * @param[in]    p          The polygon
 * @return The computed area.
 *
 * When the polygon vertices are provided in counter-clockwise
 * order, the computed area is positive.
 *
 * The polygon area is computed via Green's theorem and may
 * be either convex or concave, but should not be self-intersecting.
 */
extern double Polygon_area (const Polygon_Type *p);

/** Allocate working storage to be used for polygon clipping.
 * @return A Polygon_Clip_Type pointer on success, NULL on error.
 */
extern Polygon_Clip_Type *Polygon_open_clip (void);

/** Free working storage allocated by Polygon_open_clip
 * @param[in]  cl     Polygon_Clip_Type structure obtained by
 *                    calling Polygon_open_clip.
 */
extern void Polygon_close_clip (Polygon_Clip_Type *cl);

/** Use the intersection of two polygons to construct a new polygon.
 * @param[in]  cl     Polygon_Clip_Type structure obtained by
 *                    calling Polygon_open_clip.
 * @param[in]  sub    A convex "subject" polygon.
 * @param[in]  clip   A convex "clipper" polygon.
 * @return A Polygon_Type pointer on success, NULL on error.
 *
 * The vertices defining the overlap polygon are returned
 * in a Polygon_Type structure.  When the input polygons
 * do not overlap, the returned structure has length zero.
 * When the returned structure is no longer needed, use
 * Polygon_free to free the associated memory.
 *
 * Polygon clipping is done using the Sutherland-Hodgman
 * algorithm, which works only for convex polygons.
 * Note that the output polygon need not be convex.
 */
extern Polygon_Type *Polygon_clip (Polygon_Clip_Type *cl,
                                   const Polygon_Type *sub,
                                   const Polygon_Type *clip);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
