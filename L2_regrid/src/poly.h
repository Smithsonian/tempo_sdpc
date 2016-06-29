#ifndef __REGRID_GEOM_H__
#define __REGRID_GEOM_H__ 1

typedef struct Polygon_Type Polygon_Type;
typedef struct Polygon_Clip_Type Polygon_Clip_Type;

extern void Polygon_free (Polygon_Type *p);
extern Polygon_Type *Polygon_new (int n);

extern int Polygon_set (Polygon_Type *p, int n, const double *x, const double *y);
extern int Polygon_add (Polygon_Type *p, double x, double y);
extern int Polygon_length (const Polygon_Type *p);
extern int Polygon_vertex (const Polygon_Type *p, int i,
                           double *x, double *y);
extern void Polygon_bbox (const Polygon_Type *p,
                          double *xmin, double *xmax, double *ymin, double *ymax);
extern double Polygon_area (const Polygon_Type *p);

extern Polygon_Clip_Type *Polygon_open_clip (void);
extern void Polygon_close_clip (Polygon_Clip_Type *cl);
extern Polygon_Type *Polygon_clip (Polygon_Clip_Type *cl,
                                   const Polygon_Type *sub,
                                   const Polygon_Type *clip);
#endif
