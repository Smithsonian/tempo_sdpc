#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <poly.h>
#include <netcdf.h>
#include <proj_api.h>
#include <geodesic.h>

#include "poly.h"

#define BUFSIZE 1024

/* WGS84 ellipsoid definition */
#define WGS84_SEMIMAJOR_AXIS     6378137.0        /* meters */
#define WGS84_FLATTENING_FACTOR  1/298.257223563

#define GEO_ALTITUDE   35785831.0   /* meters */

#define SAT_LONGITUDE  -100.0       /* deg */
#define AIM_TILT          5.53      /* deg, away from nadir, >0 is northward */
#define AIM_AZIMUTH       4.50      /* deg, >0 is CW rotation, eastward from north */
#define AIM_LONGITUDE   -96.9657    /* deg */
#define AIM_LATITUDE     33.9239    /* deg */

#define CCD_PIXEL_SIZE    41.79e-6  /* radians */
#define MIRROR_STEP_SIZE  114.0e-6  /* radians */
#define SLIT_WIDTH        120.0e-6  /* radians */

/* Nominal TEMPO grid parameters */
#define NUM_STEPS    1280
#define NUM_PIXELS   2048

typedef struct
{
   projPJ tpers;
   projPJ latlong;
   projPJ albers;
   double x_bs;
   double y_bs;
   double step_size;
   double slit_width;
   double pixel_size;
   int num_steps;
   int num_pixels;
}
Obs_Type;

static void close_obs (Obs_Type *o)
{
   pj_free(o->tpers);
   pj_free(o->latlong);
   pj_free(o->albers);
}

static int init_proj (Obs_Type *o)
{
   char buf[BUFSIZE];
   double h    = GEO_ALTITUDE;
   double lon0 = SAT_LONGITUDE;
   double tilt = AIM_TILT;
   double azi  = AIM_AZIMUTH;
   const char tpers_string[] =
     "+proj=tpers +lat_0=0 +lon_0=%0.2f +h=%9.2f +tilt=%0.4f +azi=%0.4f +ellps=WGS84";
   const char albers_string[] =
     "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
   int n;

   n = snprintf (buf, BUFSIZE, tpers_string, lon0, h, tilt, azi);
   if (n >= BUFSIZE)
     {
        fprintf (stderr, "***ERROR:  proj4 format string truncated!!\n");
        return -1;
     }

   if ((NULL == (o->tpers = pj_init_plus (buf)))
       || (NULL == (o->latlong = pj_latlong_from_proj (o->tpers))))
     {
        fprintf (stderr, "***ERROR: proj4 initialization failed:\n  %s\n", buf);
        return -1;
     }

   n = snprintf (buf, BUFSIZE, tpers_string, lon0, h, tilt, azi);
   if (n >= BUFSIZE)
     {
        fprintf (stderr, "***ERROR:  proj4 format string truncated!!\n");
        return -1;
     }
   if (NULL == (o->albers = pj_init_plus (albers_string)))
     {
        fprintf (stderr, "***ERROR: proj4 initialization failed:\n  %s\n", buf);
        return -1;
     }

   return 0;
}

static int open_obs (Obs_Type *o)
{
   double x = AIM_LONGITUDE;
   double y = AIM_LATITUDE;
   int status;

   if (-1 == init_proj (o))
     return -1;

   x *= DEG_TO_RAD;
   y *= DEG_TO_RAD;
   status = pj_transform (o->latlong, o->tpers, 1, 1, &x, &y, NULL);
   if (status)
     {
        fprintf (stderr, "*** ERROR: pj_transform: status=%d (%s)\n",
                 status, pj_strerrno(status));
        return -1;
     }

   o->step_size = MIRROR_STEP_SIZE;   /* mirror step size [radians] */
   o->slit_width = SLIT_WIDTH;        /* spectrometer slit width [radians] */
   o->pixel_size = CCD_PIXEL_SIZE;    /* N/S pixel size [radians] */
   o->x_bs = x / GEO_ALTITUDE;   /* tpers angular X coordinate of aim point [radians] */
   o->y_bs = y / GEO_ALTITUDE;   /* tpers angular Y coordinate of aim point [radians] */

   o->num_steps = NUM_STEPS;
   o->num_pixels = NUM_PIXELS;

   return 0;
}

static int define_latlong_region (Obs_Type *o, double x0, double y0,
                                  double size, int num_side,
                                  int *num, double **px, double **py)
{
   int i, k, status, num_points = 4 * num_side;
   double *x = NULL, *y = NULL;
   double ds;

   if ((NULL == (x = (double *)malloc (num_points * sizeof(double))))
       ||(NULL == (y = (double *)malloc (num_points * sizeof(double)))))
     {
        free(y);
        return -1;
     }

   /* Basic region is a rectangular pixel in the virtual
    * detector plane.
    */

   /* tpers coordinates in meters */
   x0   *= GEO_ALTITUDE;
   y0   *= GEO_ALTITUDE;
   size *= GEO_ALTITUDE;

   ds = size / num_side;

   i = 0;
   for (k = 0; k < num_side; k++)
     {
        x[i] = x0 + k * ds;
        y[i] = y0;
        i++;
     }
   for (k = 0; k < num_side; k++)
     {
        x[i] = x0 + size;
        y[i] = y0 + k * ds;
        i++;
     }
   for (k = 0; k < num_side; k++)
     {
        x[i] = (x0 + size) - k * ds;
        y[i] = y0 + size;
        i++;
     }
   for (k = 0; k < num_side; k++)
     {
        x[i] = x0;
        y[i] = (y0 + size) - k * ds;
        i++;
     }

   status = pj_transform (o->tpers, o->latlong, num_points, 1, x, y, NULL);
   if (status)
     {
        fprintf (stderr, "*** Error: pj_transform failed, status = %d (%s)",
                 status, pj_strerrno(status));
        free(x);
        free(y);
        return -1;
     }

   for (i = 0; i < num_points; i++)
     {
        x[i] /= DEG_TO_RAD;
        y[i] /= DEG_TO_RAD;
     }

   *num = num_points;
   *px = x;
   *py = y;

   return 0;
}

static int compare_areas (Obs_Type *o, double x0, double y0,
                          double tpers_box_size, int num_extra_per_side,
                          double *geodesic_area, double *albers_area)
{
   Polygon_Type *p = NULL;
   struct geod_geodesic g;
   double *lons=NULL, *lats=NULL;
   double geodesic_perimeter;
   int i, n, status = -1;

   /* The "exact" calculation uses a densified polygon with
    * 4*num_extra_per_side additional vertices connected by
    * geodesics */
   if (-1 == define_latlong_region (o, x0, y0, tpers_box_size,
                                    1 + num_extra_per_side,
                                    &n, &lons, &lats))
     return -1;

   geod_init(&g, WGS84_SEMIMAJOR_AXIS, WGS84_FLATTENING_FACTOR);
   geod_polygonarea (&g, lats, lons, n, geodesic_area, &geodesic_perimeter);

   free(lons);
   free(lats);

   /* The "Albers" calculation just uses the 4 corners of the region */
   if (-1 == define_latlong_region (o, x0, y0, tpers_box_size,
                                    1, &n, &lons, &lats))
     return -1;

   for (i = 0; i < n; i++)
     {
        lons[i] *= DEG_TO_RAD;
        lats[i] *= DEG_TO_RAD;
     }
   status = pj_transform (o->latlong, o->albers, n, 1, lons, lats, NULL);
   if (status)
     {
        fprintf (stderr, "*** Error: pj_transform failed, status = %d (%s)",
                 status, pj_strerrno(status));
        free(lons);
        free(lats);
        return -1;
     }
   if ((NULL == (p = Polygon_new (n)))
       || (-1 == Polygon_set (p, n, lons, lats)))
     {
        free(lons);
        free(lats);
        Polygon_free (p);
        return -1;
     }
   *albers_area = Polygon_area (p);

   free(lons);
   free(lats);
   Polygon_free (p);

   return 0;
}

static int print_area_errors (FILE *fp, Obs_Type *o,
                              double x0, double y0, int num_extra_per_side)
{
   double lon0, lat0, bin_factor;
   double geodesic_area, albers_area, tpers_box_size;
   int status;

   lon0 = x0 * GEO_ALTITUDE;
   lat0 = y0 * GEO_ALTITUDE;
   status = pj_transform (o->tpers, o->latlong, 1, 1, &lon0, &lat0, NULL);
   if (status)
     {
        fprintf (stderr, "*** Error: pj_transform failed, status = %d (%s)",
                 status, pj_strerrno(status));
        return -1;
     }

   fprintf (fp, "# Region center: lon0=%0.4f, lat0=%0.4f deg\n",
            lon0/DEG_TO_RAD, lat0/DEG_TO_RAD);
   fprintf (fp,
            "# binfac, box_size [urad], Albers area [km^2], \"Exact\" area [km^2], sqrt(Exact) [km], Frac. error\n");

   for (bin_factor = 1.0; bin_factor < 100.0; bin_factor *= 1.05)
     {
        tpers_box_size = o->slit_width * bin_factor;
        if (-1 == compare_areas (o, x0, y0, tpers_box_size, num_extra_per_side,
                                 &geodesic_area, &albers_area))
          return -1;

        /* If the box gets so big that some vertices miss the earth,
         * those vertices project to HUGE_VAL and the computed area
         * ends up as NaN. */
        if (isnan(geodesic_area) || isnan(albers_area))
          break;

        fprintf (fp, "%6.1f %9.4e %15.8e %15.8e %9.4f %15.8e\n",
                 bin_factor, 1.e6 * tpers_box_size,
                 albers_area/1.e6,
                 geodesic_area/1.e6,
                 sqrt(geodesic_area)/1.e3,
                 1.0 - albers_area/geodesic_area);
     }

   return 0;
}

int main (int argc, char **argv)
{
   Obs_Type obs;
   FILE *fp;
   const char file0[] = "area_check_bs.dat";
   const char file1[] = "area_check_ne.dat";
   double x0, y0;
   int num_extra_per_side = 128;
   int status = 1;

   if (argc > 1)
     {
        if (1 != sscanf (argv[1], "%d", &num_extra_per_side))
          {
             fprintf (stderr, "*** Error parsing arg=%s\n", argv[1]);
             return 1;
          }
     }

   if (-1 == open_obs (&obs))
     return 1;

   fprintf (stdout, "\"Exact\" area uses num_extra_per_side = %d\n",
            num_extra_per_side);

   if (NULL == (fp = fopen (file0, "w")))
     {
        fprintf (stderr, "*** Error opening output file %s\n", file0);
        goto cleanup;
     }

   /* radians */
   x0 = obs.x_bs;
   y0 = obs.y_bs;
   status = print_area_errors (fp, &obs, x0, y0, num_extra_per_side);
   fclose (fp);
   if (status < 0)
     goto cleanup;

   if (NULL == (fp = fopen (file1, "w")))
     {
        fprintf (stderr, "*** Error opening output file %s\n", file1);
        goto cleanup;
     }

   /* Larger errors are expected near northern corners of the FOV
    * where the projected regions are most distorted.
    */
   x0 += (obs.step_size * obs.num_steps/2) * 0.77;
   y0 += (obs.pixel_size * obs.num_pixels/2) * 0.77;
   status = print_area_errors (fp, &obs, x0, y0, num_extra_per_side);
   fclose(fp);
   if (status < 0)
     goto cleanup;

   status = 0;
cleanup:
   close_obs (&obs);
   return status;
}
