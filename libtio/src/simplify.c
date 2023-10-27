#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef FREE
# define FREE free
#endif

#ifndef MALLOC
# define MALLOC malloc
#endif

#ifndef REALLOC
# define REALLOC realloc
#endif

#define EARTH_RADIUS_KM (6371.007181)
#define DEGTORAD        (M_PI / 180.0)

#define INITIAL_STACK_SIZE 16

typedef struct
{
   int num;
   int num_alloc;
   int *indices;
}
Stack_Type;

static void stack_free (Stack_Type *s)
{
   if (s == NULL)
     return;
   FREE(s->indices);
   memset ((char *)s, 0, sizeof(*s));
}

static int stack_push (Stack_Type *s, int beg, int end)
{
   if (s->num_alloc == 0)
     {
        if (NULL == (s->indices = (int *)MALLOC (INITIAL_STACK_SIZE * 2 * sizeof(int))))
          {
             fprintf (stderr, "*** %s: malloc failed\n", __func__);
             return -1;
          }
        s->num_alloc = INITIAL_STACK_SIZE;
        s->num = 0;
     }
   else if (s->num == s->num_alloc)
     {
        int new_num = s->num_alloc * 2;
        int *new_indices = NULL;
        if (NULL == (new_indices = (int *)REALLOC (s->indices, new_num * 2 * sizeof(int))))
          {
             fprintf (stderr, "*** %s: realloc failed\n", __func__);
             return -1;
          }
        s->num_alloc = new_num;
        s->indices = new_indices;
     }

   s->indices[s->num  ] = beg;
   s->indices[s->num+1] = end;
   s->num += 2;

   return 0;
}

static int stack_pop (Stack_Type *s, int *beg, int *end)
{
   if (s->num < 2)
     return -1;

   s->num -= 2;
   *end = s->indices[s->num+1];
   *beg = s->indices[s->num  ];

   return 0;
}

typedef struct
{
   const float *lon;
   const float *lat;
   int num;
}
Point_Type;

typedef struct
{
   double dx;
   double dy;
   double len_sqr;
}
Delta_Type;

static void point_sep_sqr (Point_Type *pt, int a, int b, Delta_Type *dt)
{
   double dx = pt->lon[b] - pt->lon[a];
   double dy = pt->lat[b] - pt->lat[a];
   double lat0 = 0.5 * (pt->lat[b] + pt->lat[a]) * DEGTORAD;
   if (fabs(dx) > 180.0) dx = 360.0 - fabs(dx);
   dx *= cos (lat0);
   dt->dx = dx;
   dt->dy = dy;
   dt->len_sqr = dx*dx + dy*dy;
}

/* Douglas-Peucker polyline simplification.
 * This is a stack-based (non-recursive) implementation.
 * The algorithm is described in:
 *   Douglas, D. H., and T. K. Peucker, Algorithms for the reduction
 *   of the number of points required to represent a digitized line
 *   of its caricature, Can. Cartogr., 10, 112-122, 1973.
 * and this implementation was based on an earlier implementation by:
 *   Dr. Gary J. Robinson, Environmental Systems Science Centre,
 *   University of Reading, Reading, UK
 *
 * The polygon is assumed to be specified in geospatial longitude-latitude
 * coordinates in degrees.  The simplification eliminates details within
 * a band of width band_km kilometers.  Vertices of the simplified polygon
 * are returned in an index array:
 *
 *  num_kept = simplify_dp (lon, lat, num, band_km, &indices)
 */
int simplify_dp (const float *lon_deg, const float *lat_deg, int num,
                 float band_km, int **pindex)
{
   Stack_Type s = {0};
   Point_Type pt =
     {
        .lon = lon_deg,
        .lat = lat_deg,
        .num = num
     };
   int *index = NULL;
   double band_sqr, max_dev_sqr;
   int num_kept, beg, end, split, status = -1;

   if (num < 3)
     return -1;

   *pindex = NULL;

   if (NULL == (index = (int *)MALLOC (num * sizeof(int))))
     {
        fprintf (stderr, "%s:  malloc failed\n", __func__);
        return -1;
     }

   band_sqr = (band_km / EARTH_RADIUS_KM) / DEGTORAD;
   band_sqr *= band_sqr;  /* deg^2 */

   num_kept = 0;

   if (0 != stack_push (&s, 0, num-1))
     goto return_error;

   beg = 0;
   end = num-1;

   while (s.num > 0)
     {
        Delta_Type d12, d13, d23;
        int i;

        stack_pop (&s, &beg, &end);

        if (end - beg <= 1)
          {
             /* no intermediate points, so keep the current start point */
             index[num_kept++] = beg;
             continue;
          }

        max_dev_sqr = -1.0;
        split = beg;

        point_sep_sqr (&pt, beg, end, &d12);

        /* Within these endpoints, find the largest deviation from a straight line */
        for (i = beg+1; i < end; i++)
          {
             double dev_sqr;

             point_sep_sqr (&pt, beg, i, &d13);
             point_sep_sqr (&pt, i, end, &d23);

             if (d13.len_sqr >= d12.len_sqr + d23.len_sqr)
               {
                  dev_sqr = d23.len_sqr;
               }
             else if (d23.len_sqr >= d12.len_sqr + d13.len_sqr)
               {
                  dev_sqr = d13.len_sqr;
               }
             else
               {
                  double delta = d13.dx * d12.dy - d13.dy * d12.dx;
                  dev_sqr = delta * delta / d12.len_sqr;
               }

             if (dev_sqr > max_dev_sqr)
               {
                  split = i;
                  max_dev_sqr = dev_sqr;
               }
          }

        /* Was there a significant intermediate point? */
        if (max_dev_sqr < band_sqr)
          {
             /* no, so keep the current start point */
             index[num_kept++] = beg;
          }
        else
          {
             /* yes, so push two segments on the stack for further processing */
             if (0 != stack_push (&s, split, end))
               goto return_error;
             if (0 != stack_push (&s, beg, split))
               goto return_error;
          }
     }

   /* keep the last point */
   index[num_kept++] = num-1;

   status = 0;
return_error:
   stack_free (&s);

   if (status)
     {
        FREE(index);
        index = NULL;
     }
   *pindex = index;

   return status ? -1 : num_kept;
}
