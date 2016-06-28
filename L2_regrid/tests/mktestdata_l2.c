#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <poly.h>
#include <netcdf.h>
#include <proj_api.h>

#define MALLOC  malloc
#define FREE    free

#define BUFSIZE 1024

#define GEO_ALTITUDE   35785831.0   /* meters */

#define SAT_LONGITUDE  -100.0       /* deg */
#define AIM_TILT          5.53      /* deg, away from nadir, >0 is northward */
#define AIM_AZIMUTH       4.50      /* deg, >0 is CW rotation, eastward from north */
#define AIM_LONGITUDE   -96.9657    /* deg */
#define AIM_LATITUDE     33.9239    /* deg */

#define CCD_PIXEL_SIZE    41.79e-6  /* radians */
#define MIRROR_STEP_SIZE  114.0e-6  /* radians */
#define SLIT_WIDTH        120.0e-6  /* radians */

#if 1
#define NUM_STEPS  1280
#define NUM_PIXELS 2048
#else
#define NUM_STEPS  40
#define NUM_PIXELS 25
#endif

#define NUM_GRANULES  10

#define NC_CHECK_STATUS(s) \
   do {if (NC_NOERR != (s)) goto cleanup_and_exit; } while (0);

typedef struct
{
   projPJ tpers;
   projPJ latlong;
   double x_bs;
   double y_bs;
   double step_size;
   double slit_width;
   double pixel_size;
   int num_steps;
   int num_pixels;
}
Obs_Type;

typedef struct
{
   double *x0_bounds;  /* xll0,xlr0,xur0,xul0, xll1,xlr1,xur1,xul1, ... */
   double *x1_bounds;
   double *x0;
   double *x1;
   int num_pixels;
   int slit_pos;
}
Slit_Pixel_List_Type;

static void close_obs (Obs_Type *o)
{
   pj_free(o->tpers);
   pj_free(o->latlong);
}

static int init_proj (Obs_Type *o)
{
   char buf[BUFSIZE];
   int n;
   double h    = GEO_ALTITUDE;
   double lon0 = SAT_LONGITUDE;
   double tilt = AIM_TILT;
   double azi  = AIM_AZIMUTH;
   const char ctl_string[] =
     "+proj=tpers +lat_0=0 +lon_0=%0.2f +h=%9.2f +tilt=%0.4f +azi=%0.4f +ellps=WGS84";

   n = snprintf (buf, BUFSIZE, ctl_string, lon0, h, tilt, azi);
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

   o->num_steps = NUM_STEPS;
   o->num_pixels = NUM_PIXELS;
   o->step_size = MIRROR_STEP_SIZE;   /* mirror step size [radians] */
   o->slit_width = SLIT_WIDTH;        /* spectrometer slit width [radians] */
   o->pixel_size = CCD_PIXEL_SIZE;    /* N/S pixel size [radians] */
   o->x_bs = x / GEO_ALTITUDE;   /* tpers angular X coordinate of aim point [radians] */
   o->y_bs = y / GEO_ALTITUDE;   /* tpers angular Y coordinate of aim point [radians] */

   return 0;
}

static void free_slit_grid (Slit_Pixel_List_Type *g)
{
   if (g == NULL)
     return;
   FREE(g->x0_bounds);
   FREE(g->x1_bounds);
   FREE(g);
}

static Slit_Pixel_List_Type *new_slit_grid (int n)
{
   Slit_Pixel_List_Type *g = NULL;
   int len = n*sizeof(double);

   if (NULL == (g = (Slit_Pixel_List_Type *)MALLOC(sizeof *g)))
     return NULL;

   if ( (NULL == (g->x0_bounds = (double *)MALLOC (5*len)))
       || (NULL == (g->x1_bounds = (double *)MALLOC (5*len))))
     {
        free_slit_grid (g);
        return NULL;
     }
   g->x0 = g->x0_bounds + 4*n;
   g->x1 = g->x1_bounds + 4*n;
   g->num_pixels = n;
   return g;
}

static int make_slit_grid (Obs_Type *o, Slit_Pixel_List_Type *g, int slit_pos)
{
   double xoff, yoff, x0, y0, y1, dy, hdx, hdy, cs, sn;
   double *g_x0, *g_x1, *x, *y;
   int i, k, n, status;

   x = g->x0_bounds;
   y = g->x1_bounds;
   g_x0 = g->x0;
   g_x1 = g->x1;

   /* Make detector grid of the corners of rectangular "pixels"
    * for a single mirror position.
    *
    * Note that pixels in the field of regard (FOR) will overlap
    * in X direction because
    *        pixel spacing (= mirror step size)
    * is smaller than
    *        pixel width (= slit width).
    * In the Y direction, the pixel size is just the
    * angular size of a physical CCD pixel, with no overlaps.
    */
   dy = o->pixel_size;

   /* Angular offset of southern-most end of the slit grid from the
    * boresight location. */
   xoff = o->step_size * (o->num_steps/2 - slit_pos);
   yoff = - dy * o->num_pixels/2;

   /* Account for the rotation of the instrument mounting
    * about the boresight by AIM_AZIMUTH. */
   cs = cos(-AIM_AZIMUTH * DEG_TO_RAD);
   sn = sin(-AIM_AZIMUTH * DEG_TO_RAD);
   xoff = xoff * cs + yoff * sn;
   yoff = xoff * (-sn) + yoff * cs;

   /* Shift to get absolute tpers coordinates of the slit grid
    * southern-most end (tpers = tilted perspective projection). */
   x0 = xoff + o->x_bs;
   y0 = yoff + o->y_bs;

   /* pixel half-width, half-height */
   hdx = 0.5 * o->slit_width;
   hdy = 0.5 * dy;

   k = 0;
   for (i = 0; i < o->num_pixels; i++)
     {
        y1 = y0 + i*dy;

        g_x0[i] = x0;
        g_x1[i] = y1;

        x[k]   = x0 - hdx;
        y[k++] = y1 - hdy;

        x[k]   = x0 + hdx;
        y[k++] = y1 - hdy;

        x[k]   = x0 + hdx;
        y[k++] = y1 + hdy;

        x[k]   = x0 - hdx;
        y[k++] = y1 + hdx;
     }

   g->num_pixels = o->num_pixels;
   g->slit_pos = slit_pos;

   /* Project FOR pixel corners onto Earth's surface
    * (Notice that I'm exploiting the fact that both
    *  pixel centers and corners are packed into the same array.)
    */

   n = 5 * o->num_pixels;
   for (i = 0; i < n; i++)
     {
        x[i] *= GEO_ALTITUDE;
        y[i] *= GEO_ALTITUDE;
     }
   status = pj_transform (o->tpers, o->latlong, n, 1, x, y, NULL);
   if (status)
     {
        fprintf (stderr, "*** ERROR: make_slit_grid: pj_transform returned status=%d (%s)\n",
                 status, pj_strerrno(status));
        return -1;
     }

#define MAKE_STORABLE_FLOAT(a) \
   do {if (isinf(a) || isnan(a) || (a == HUGE_VAL)) a = FLT_MAX; } while (0)

   for (i = 0; i < n; i++)
     {
        x[i] /= DEG_TO_RAD;
        y[i] /= DEG_TO_RAD;

        MAKE_STORABLE_FLOAT (x[i]);
        MAKE_STORABLE_FLOAT (y[i]);
     }

   return 0;
}

static int use_compression (int ncid, int varid)
{
   int shuffle = 1;
   int deflate = 1;
   int deflate_level = 1;
   return nc_def_var_deflate (ncid, varid, shuffle, deflate, deflate_level);
}

int main (void)
{
   Obs_Type o = {0};
   Slit_Pixel_List_Type *g = NULL;
   const char coord_lonlat[] = "longitude latitude";
   const char bounds_lonlat[] = "longitude_bounds latitude_bounds";
   const char units_lon[] = "degrees_east";
   const char units_lat[] = "degrees_north";
   int coord_type = NC_FLOAT;
   double *column = NULL;
   int *xtrack = NULL;
   int dims[3];
   size_t i_sizet, start[3], count[3];
   int ncid, id_step, id_xtrack;
   int id_lon_bounds, id_lat_bounds;
   int id_lon, id_lat, id_column;
   int num_steps_per_granule = NUM_STEPS / NUM_GRANULES;
   int granule, num_granules = NUM_GRANULES;
   int i, step, status = 1;

   if (-1 == open_obs (&o))
     return 1;

   if (NULL == (g = new_slit_grid (o.num_pixels)))
     return 1;

   if ((NULL == (column = (double *) MALLOC (o.num_pixels * sizeof(double))))
       || (NULL == (xtrack = (int *) MALLOC (o.num_pixels * sizeof(int)))))
     return 1;

   for (i = 0; i < o.num_pixels; i++)
     {
        xtrack[i] = i;
     }

   step = 0;
   for (granule = 0; granule < num_granules; granule++)
     {
        char outfile[BUFSIZE];
        int n;

        n = sprintf (outfile, "/tmp/test_l2l3_g%02d_grid.nc", granule);
        if (n >= BUFSIZE)
          {
             fprintf (stderr, "**** sprintf failed!!\n");
             goto cleanup_and_exit;
          }
        fprintf (stderr, "Writing %s\n", outfile);

        status = nc_create (outfile, NC_NETCDF4, &ncid);
        NC_CHECK_STATUS(status);

        status = nc_def_dim (ncid, "corner", 4, &dims[2]);
        NC_CHECK_STATUS(status);
        status = nc_def_dim (ncid, "xtrack", o.num_pixels, &dims[1]);
        NC_CHECK_STATUS(status);
        status = nc_def_dim (ncid, "mirror_step", num_steps_per_granule, &dims[0]);
        NC_CHECK_STATUS(status);
        status = nc_def_var (ncid, "mirror_step", NC_INT, 1, &dims[0], &id_step);
        NC_CHECK_STATUS(status);
        status = nc_def_var (ncid, "xtrack", NC_INT, 1, &dims[1], &id_xtrack);
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "column", NC_FLOAT, 2, dims, &id_column);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_column);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_column, "coordinates", strlen(coord_lonlat), coord_lonlat);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_column, "bounds", strlen(bounds_lonlat), bounds_lonlat);
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "longitude", coord_type, 2, dims, &id_lon);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_lon);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_lon, "units", strlen(units_lon), units_lon);
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "latitude", coord_type, 2, dims, &id_lat);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_lat);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_lat, "units", strlen(units_lat), units_lat);
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "longitude_bounds", coord_type, 3, dims, &id_lon_bounds);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_lon_bounds);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_lon_bounds, "units", strlen(units_lon), units_lon);
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "latitude_bounds", coord_type, 3, dims, &id_lat_bounds);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_lat_bounds);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_lat_bounds, "units", strlen(units_lat), units_lat);
        NC_CHECK_STATUS(status);

        start[0] = 0;
        count[0] = o.num_pixels;
        status = nc_put_vara_int (ncid, id_xtrack, &start[0], &count[0], xtrack);
        NC_CHECK_STATUS(status);

        for (i = 0; i < num_steps_per_granule; i++)
          {
             int j;

             if (step == o.num_steps)
               break;

             if (-1 == make_slit_grid (&o, g, step))
               goto cleanup_and_exit;

             for (j = 0; j < o.num_pixels; j++)
               {
#if 0
                  double r = hypot (j - o.num_pixels*0.5,
                                    step - o.num_steps*0.5);
                  column[j] = 100.0 * exp(-r/(o.num_steps*0.25));
#else
                  column[j] = 1.0 * ((step/10) % 10);
#endif
               }

             start[0] = i;    /* step */
             start[1] = 0;    /* xtrack */
             start[2] = 0;    /* corner */
             count[0] = 1;
             count[1] = o.num_pixels;
             count[2] = 4;
             i_sizet = i;
             status = nc_put_vara_int (ncid, id_step, &i_sizet, &count[0], &step);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lon_bounds, start, count, g->x0_bounds);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lat_bounds, start, count, g->x1_bounds);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lon, start, count, g->x0);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lat, start, count, g->x1);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_column, start, count, column);
             NC_CHECK_STATUS(status);

             step++;
          }

        if (NC_NOERR != (status = nc_close (ncid)))
          {
             fprintf (stderr, "*** ERROR: closing file %s\n", outfile);
          }

        ncid = -1;
     }

   status = 0;
cleanup_and_exit:
   if (status)
     {
        fprintf (stderr, "*** ERROR: %s\n", nc_strerror(status));
     }

   FREE(column);
   FREE(xtrack);
   free_slit_grid (g);
   close_obs (&o);

   return status ? 1 : 0;
}
