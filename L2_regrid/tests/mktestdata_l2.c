#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include <poly.h>
#include <netcdf.h>
#include <proj_api.h>

#include <tio.h>
#include <tio_template.h>

#define MALLOC  malloc
#define FREE    free

#define BUFSIZE      1024
#define __MAX_VAR_DIMS 7

#define LOCAL_DIM_SIZE1  3
#define LOCAL_DIM_SIZE2  2

#define GEO_ALTITUDE   35785831.0   /* meters */

#define SAT_LONGITUDE  -100.0       /* deg */
#define AIM_TILT          5.53      /* deg, away from nadir, >0 is northward */
#define AIM_AZIMUTH       4.50      /* deg, >0 is CW rotation, eastward from north */
#define AIM_LONGITUDE   -96.9657    /* deg */
#define AIM_LATITUDE     33.9239    /* deg */

#define CCD_PIXEL_SIZE    41.79e-6  /* radians */
#define MIRROR_STEP_SIZE  114.0e-6  /* radians */
#define SLIT_WIDTH        120.0e-6  /* radians */

#define NUM_STEPS    1280
#define NUM_PIXELS   2048
#define NUM_GRANULES   10

#define OUTPUT_DIR    "."

#define NC_CHECK_STATUS(s) \
   do { \
      if (NC_NOERR != (s)) \
        { \
           fprintf (stderr, "*** NC_ERROR = %d on line %d\n", s, __LINE__); \
           goto cleanup_and_exit; \
        } \
   } while (0);

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
   int num_granules;
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

static int open_obs (Obs_Type *o,
                     int num_steps, int num_pixels, int num_granules)
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

   o->num_granules = (num_steps > num_granules) ? num_granules : 1;
   o->num_steps = num_steps;
   o->num_pixels = num_pixels;

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
        y1 = y0 + (o->num_pixels-1-i)*dy;

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

static int put_lonlat_atts (int ncid, int id_var)
{
   const char coord_lonlat[] = "longitude latitude";
   const char bounds_lonlat[] = "longitude_bounds latitude_bounds";
   int status;
   status = nc_put_att_text (ncid, id_var, "coordinates", strlen(coord_lonlat), coord_lonlat);
   if (NC_NOERR != status) return status;
   status = nc_put_att_text (ncid, id_var, "bounds", strlen(bounds_lonlat), bounds_lonlat);
   return status;
}

static int put_att_text (int ncid, int id_var,
                         const char *name, const char *text)
{
   return nc_put_att_text (ncid, id_var, name, strlen(text), text);
}

static void update_bbox (Slit_Pixel_List_Type *g, double *bbox)
{
   int i, n = g->num_pixels;
   for (i = 0; i < n; i++)
     {
        double *lon = g->x0_bounds + 4*i;
        double *lat = g->x1_bounds + 4*i;
        if ((fabs(lon[0]) < 360.0) && (lon[0] < bbox[0])) bbox[0] = lon[0];
        if ((fabs(lon[1]) < 360.0) && (lon[1] > bbox[1])) bbox[1] = lon[1];
        if ((fabs(lat[0]) <  90.0) && (lat[0] < bbox[2])) bbox[2] = lat[0];
        if ((fabs(lat[2]) <  90.0) && (lat[2] > bbox[3])) bbox[3] = lat[2];
     }
}

static void usage (int argc, char **argv)
{
   (void) argc;
   fprintf (stderr, "Usage: %s [args]\n", argv[0]);
   fprintf (stderr, "Options:\n");
   fprintf (stderr, "  -o <out_dir>         Path to output directory [%s]\n", OUTPUT_DIR);
   fprintf (stderr, "  -g <num_granules>    Number of granules [%d]\n", NUM_GRANULES);
   fprintf (stderr, "  -m <num_steps>       Number of mirror steps [%d]\n", NUM_STEPS);
   fprintf (stderr, "  -p <num_pixels>      Number of pixels along slit [%d]\n", NUM_PIXELS);
}

int main (int argc, char **argv)
{
   Obs_Type o = {0};
   Slit_Pixel_List_Type *g = NULL;
   const char units_lon[] = "degrees_east";
   const char units_lat[] = "degrees_north";
   int coord_type = NC_FLOAT;
   double *var0 = NULL, *var1 = NULL, *var2 = NULL, *var3 = NULL;
   short *var4 = NULL;
   const char *bitfield_names[] = {
      "ul_bitfield","ui_bitfield", "us_bitfield", "uc_bitfield",
       "l_bitfield", "i_bitfield",  "s_bitfield",  "c_bitfield"
     };
   int bitfield_types[] = {
      NC_UINT64, NC_UINT, NC_USHORT, NC_UBYTE,
       NC_INT64,  NC_INT,  NC_SHORT,  NC_BYTE
   };
   unsigned long long *ul_bitfield = NULL;
   unsigned int       *ui_bitfield = NULL;
   unsigned short     *us_bitfield = NULL;
   unsigned char      *uc_bitfield = NULL;
   long long          * l_bitfield = NULL;
   int                * i_bitfield = NULL;
   short              * s_bitfield = NULL;
   signed char        * c_bitfield = NULL;
   int id_bitfield[8];
   int *xtrack = NULL;
   int dims[__MAX_VAR_DIMS], dimid_corner;
   size_t step_start, start[__MAX_VAR_DIMS], count[__MAX_VAR_DIMS];
   int ncid, grp, id_step, id_xtrack;
   int id_lon_bounds, id_lat_bounds;
   int id_lon, id_lat, id_var0, id_var1, id_var2, id_var3, id_var4;
   int num_steps_per_granule, granule, i, step;
   int dimid_local1, num_local1 = LOCAL_DIM_SIZE1;
   int dimid_local2, num_local2 = LOCAL_DIM_SIZE2;
   int num_steps = NUM_STEPS;
   int num_pixels = NUM_PIXELS;
   int num_granules = NUM_GRANULES;
   int b, c, num_var2, num_var3, pattern_scale;
   double bbox[4] = {DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX};
   double tstart, tend, delta_step = 3.0;
   float float_fill = NC_FILL_FLOAT, float_valid_min = 0.0, float_valid_max;
   const char *out_dir = OUTPUT_DIR;
   int status = 1;

   while ((c = getopt (argc, argv, "g:m:o:p:")) != -1)
     {
        switch (c)
          {
           case 'g':
             if ((1 != sscanf (optarg, "%d", &num_granules))
                 || (num_granules < 1 || NUM_STEPS < num_granules))
               {
                  fprintf (stderr, "*** invalid num_granules=%d\n", num_granules);
                  return 1;
               }
             break;
           case 'm':
             if ((1 != sscanf (optarg, "%d", &num_steps))
                 || (num_steps < 1 || NUM_STEPS < num_steps))
               {
                  fprintf (stderr, "*** invalid num_steps=%d\n", num_steps);
                  return 1;
               }
             break;
           case 'o':
             out_dir = optarg;
             break;
           case 'p':
             if ((1 != sscanf (optarg, "%d", &num_pixels))
                 || (num_pixels < 1) || (NUM_PIXELS < num_pixels))
               {
                  fprintf (stderr, "*** invalid num_pixels=%d\n", num_pixels);
                  return 1;
               }
             break;

           case '?':
             fprintf (stderr, "Unknown option -%c'.\n", optopt);
             usage (argc, argv);
             return 1;

           default:
             usage (argc, argv);
             return 0;
          }
     }

   if (-1 == open_obs (&o, num_steps, num_pixels, num_granules))
     return 1;

   if (NULL == (g = new_slit_grid (o.num_pixels)))
     return 1;

   num_steps_per_granule = o.num_steps / o.num_granules;

   num_var2 = o.num_pixels * num_local1;
   num_var3 = o.num_pixels * num_local1 * num_local2;
   if ((NULL == (xtrack = (int *) MALLOC (o.num_pixels * sizeof(int))))
       || (NULL == (var0 = (double *) MALLOC (o.num_pixels * sizeof(double))))
       || (NULL == (var1 = (double *) MALLOC (o.num_pixels * sizeof(double))))
       || (NULL == (var2 = (double *) MALLOC (num_var2 * sizeof(double))))
       || (NULL == (var3 = (double *) MALLOC (num_var3 * sizeof(double))))
       || (NULL == (var4 = (short *) MALLOC (o.num_pixels * sizeof(short))))
       || (NULL == (ul_bitfield = (unsigned long long *) MALLOC (o.num_pixels * sizeof(long long))))
       || (NULL == (ui_bitfield = (unsigned int *) MALLOC (o.num_pixels * sizeof(int))))
       || (NULL == (us_bitfield = (unsigned short *) MALLOC (o.num_pixels * sizeof(short))))
       || (NULL == (uc_bitfield = (unsigned char *) MALLOC (o.num_pixels * sizeof(char))))
       || (NULL == ( l_bitfield = (long long *) MALLOC (o.num_pixels * sizeof(long long))))
       || (NULL == ( i_bitfield = (int *) MALLOC (o.num_pixels * sizeof(int))))
       || (NULL == ( s_bitfield = (short *) MALLOC (o.num_pixels * sizeof(short))))
       || (NULL == ( c_bitfield = (signed char *) MALLOC (o.num_pixels * sizeof(char))))
      )
     return 1;

   for (i = 0; i < o.num_pixels; i++)
     {
        xtrack[i] = i;
     }

   pattern_scale = o.num_steps / 32;
   if (pattern_scale == 0) pattern_scale = 1;

   tstart = 0.0;

   step = 0;
   for (granule = 0; granule < o.num_granules; granule++)
     {
        char outfile[BUFSIZE];
        int processing_version = 1;
        int n, num_steps_left, num_steps_this_granule;

        num_steps_left = o.num_steps - granule * num_steps_per_granule;
        if (num_steps_left < num_steps_per_granule)
          num_steps_this_granule = num_steps_left;
        else
          num_steps_this_granule = num_steps_per_granule;

        n = sprintf (outfile, "%s/test_l2l3_g%02d_grid.nc", out_dir, granule);
        if (n >= BUFSIZE)
          {
             fprintf (stderr, "**** sprintf failed!!\n");
             goto cleanup_and_exit;
          }
        fprintf (stderr, "Writing %s\n", outfile);

        status = nc_create (outfile, NC_NETCDF4, &ncid);
        NC_CHECK_STATUS(status);

        tend = tstart + num_steps_this_granule * delta_step;

        if ((0 != tio_time_set_taix_epoch ("2000-01-01T00:00:00Z"))
#if 0
            || (0 != tio_write_epoch_timestamp (ncid, NC_GLOBAL))
#else
            || (0 != tio_write_granule_ident_times (ncid, tstart, tend))
            || (0 != tio_write_granule_ident_indices (ncid, 1, granule+1))
            || (0 != TIO_label_product (ncid, "TEST", 1, 1))
#endif
           )
          goto cleanup_and_exit;

        status = nc_put_att_int (ncid, NC_GLOBAL, "processing_version", NC_INT, 1, &processing_version);
        NC_CHECK_STATUS(status);

        status = nc_def_grp (ncid, "test_group", &grp);
        NC_CHECK_STATUS(status);

        status = nc_def_dim (ncid, "local1", num_local1, &dimid_local1);
        NC_CHECK_STATUS(status);
        status = nc_def_dim (ncid, "local2", num_local2, &dimid_local2);
        NC_CHECK_STATUS(status);
        status = nc_def_dim (ncid, "corner", 4, &dimid_corner);
        NC_CHECK_STATUS(status);
        status = nc_def_dim (ncid, "xtrack", o.num_pixels, &dims[1]);
        NC_CHECK_STATUS(status);
        status = nc_def_dim (ncid, "mirror_step", num_steps_this_granule, &dims[0]);
        NC_CHECK_STATUS(status);
        status = nc_def_var (ncid, "mirror_step", NC_INT, 1, &dims[0], &id_step);
        NC_CHECK_STATUS(status);
        status = nc_def_var (ncid, "xtrack", NC_INT, 1, &dims[1], &id_xtrack);
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

        dims[2] = dimid_corner;
        status = nc_def_var (ncid, "longitude_bounds", coord_type, 3, dims, &id_lon_bounds);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_lon_bounds);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_lon_bounds, "units", strlen(units_lon), units_lon);
        NC_CHECK_STATUS(status);

        dims[2] = dimid_corner;
        status = nc_def_var (ncid, "latitude_bounds", coord_type, 3, dims, &id_lat_bounds);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_lat_bounds);
        NC_CHECK_STATUS(status);
        status = nc_put_att_text (ncid, id_lat_bounds, "units", strlen(units_lat), units_lat);
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "var0", NC_FLOAT, 2, dims, &id_var0);
        NC_CHECK_STATUS(status);
        status = nc_put_att_float (ncid, id_var0, "valid_min", NC_FLOAT, 1, &float_valid_min);
        NC_CHECK_STATUS(status);
        float_valid_max = 1280.0;
        status = nc_put_att_float (ncid, id_var0, "valid_max", NC_FLOAT, 1, &float_valid_max);
        NC_CHECK_STATUS(status);
        status = nc_def_var_fill (ncid, id_var0, 0, &float_fill);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_var0);
        NC_CHECK_STATUS(status);
        status = put_lonlat_atts (ncid, id_var0);
        NC_CHECK_STATUS(status);
        status = put_att_text (ncid, id_var0, "comment", "This is var0");
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "var1", NC_FLOAT, 2, dims, &id_var1);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_var1);
        NC_CHECK_STATUS(status);
        status = put_lonlat_atts (ncid, id_var1);
        NC_CHECK_STATUS(status);
        status = put_att_text (ncid, id_var1, "comment", "This is var1");
        NC_CHECK_STATUS(status);

        status = nc_def_var (ncid, "var4", NC_SHORT, 2, dims, &id_var4);
        NC_CHECK_STATUS(status);
        status = use_compression (ncid, id_var4);
        NC_CHECK_STATUS(status);
        status = put_lonlat_atts (ncid, id_var4);
        NC_CHECK_STATUS(status);
        status = put_att_text (ncid, id_var4, "comment", "This is var4");
        NC_CHECK_STATUS(status);

        for (b = 0; b < 8; b++)
          {
             status = nc_def_var (ncid, bitfield_names[b], bitfield_types[b], 2, dims, &id_bitfield[b]);
             NC_CHECK_STATUS(status);
             status = use_compression (ncid, id_bitfield[b]);
             NC_CHECK_STATUS(status);
             status = put_lonlat_atts (ncid, id_bitfield[b]);
             NC_CHECK_STATUS(status);
             status = put_att_text (ncid, id_bitfield[b], "comment", "This is a bitfield");
             NC_CHECK_STATUS(status);
          }

        dims[2] = dimid_local1;
        status = nc_def_var (grp, "var2", NC_FLOAT, 3, dims, &id_var2);
        NC_CHECK_STATUS(status);
        status = use_compression (grp, id_var2);
        NC_CHECK_STATUS(status);
        status = put_lonlat_atts (grp, id_var2);  /* FIXME: 3rd dimension? */
        NC_CHECK_STATUS(status);
        status = put_att_text (grp, id_var2, "comment", "This is var2");
        NC_CHECK_STATUS(status);

        dims[2] = dimid_local1;
        dims[3] = dimid_local2;
        status = nc_def_var (grp, "var3", NC_FLOAT, 4, dims, &id_var3);
        NC_CHECK_STATUS(status);
        status = use_compression (grp, id_var3);
        NC_CHECK_STATUS(status);
        status = put_lonlat_atts (grp, id_var3);  /* FIXME: other dimensions? */
        NC_CHECK_STATUS(status);
        status = put_att_text (grp, id_var3, "comment", "This is var3");
        NC_CHECK_STATUS(status);

        /* ------- end of definitions ------- */

        start[0] = 0;
        count[0] = o.num_pixels;
        status = nc_put_vara_int (ncid, id_xtrack, &start[0], &count[0], xtrack);
        NC_CHECK_STATUS(status);

        for (i = 0; i < num_steps_this_granule; i++)
          {
             int j;

             if (-1 == make_slit_grid (&o, g, step))
               goto cleanup_and_exit;

             update_bbox (g, bbox);

             for (j = 0; j < o.num_pixels; j++)
               {
                  int k, m, nm;
                  int use_fill_value = ((abs(j-0.25*o.num_pixels) < 0.1*o.num_pixels)
                                        && (abs(step - 0.75*o.num_steps) < 0.1*o.num_steps));

                  if (use_fill_value)
                    {
                       var0[j] = NC_FILL_FLOAT;
                       var1[j] = NC_FILL_FLOAT;
                       var4[j] = NC_FILL_SHORT;
                       ul_bitfield[j] = NC_FILL_UINT64;
                       ui_bitfield[j] = NC_FILL_UINT;
                       us_bitfield[j] = NC_FILL_USHORT;
                       uc_bitfield[j] = NC_FILL_UBYTE;
                       l_bitfield[j]  = NC_FILL_INT64;
                       i_bitfield[j]  = NC_FILL_INT;
                       s_bitfield[j]  = NC_FILL_SHORT;
                       c_bitfield[j]  = NC_FILL_BYTE;
                    }
                  else
                    {
                       unsigned long long ullb;
                       var0[j] = 1.0 * (step + j/128);
                       var1[j] = 1.0 * ((j/pattern_scale) % 16);
                       var4[j] = (step/pattern_scale) % 16;

                       ullb = ((step % 32) < 16) ? 0x03 : 0x0c;
                       ul_bitfield[j] = ullb;
                       ui_bitfield[j] = ullb;
                       us_bitfield[j] = ullb;
                       uc_bitfield[j] = ullb;
                       l_bitfield[j]  = ullb;
                       i_bitfield[j]  = ullb;
                       s_bitfield[j]  = ullb;
                       c_bitfield[j]  = ullb;
                    }

                  for (k = 0; k < num_local1; k++)
                    {
                       var2[k + j * num_local1] =
                         1.0 * (((k*128 + j + step)/pattern_scale) % 16);
                    }

                  nm = num_local1 * num_local2;
                  for (m = 0; m < nm; m++)
                    {
                       var3[m + j*nm] =
                         1.0 * ((((o.num_pixels-j) + m*step)/pattern_scale) % 16);
                    }
               }

             step_start = i;

             start[0] = step_start;    /* step */
             start[1] = 0;             /* xtrack */
             start[2] = 0;             /* corner */
             count[0] = 1;
             count[1] = o.num_pixels;
             count[2] = 4;
             status = nc_put_vara_int (ncid, id_step, &step_start, &count[0], &step);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lon_bounds, start, count, g->x0_bounds);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lat_bounds, start, count, g->x1_bounds);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lon, start, count, g->x0);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_lat, start, count, g->x1);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_var0, start, count, var0);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_double (ncid, id_var1, start, count, var1);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_short (ncid, id_var4, start, count, var4);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_ulonglong (ncid, id_bitfield[0], start, count, ul_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_uint (ncid, id_bitfield[1], start, count, ui_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_ushort (ncid, id_bitfield[2], start, count, us_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_ubyte (ncid, id_bitfield[3], start, count, uc_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_longlong (ncid, id_bitfield[4], start, count, l_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_int (ncid, id_bitfield[5], start, count, i_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_short (ncid, id_bitfield[6], start, count, s_bitfield);
             NC_CHECK_STATUS(status);
             status = nc_put_vara_schar (ncid, id_bitfield[7], start, count, c_bitfield);
             NC_CHECK_STATUS(status);
             start[2] = 0;
             count[2] = num_local1;
             status = nc_put_vara_double (grp, id_var2, start, count, var2);
             NC_CHECK_STATUS(status);
             start[3] = 0;
             count[3] = num_local2;
             status = nc_put_vara_double (grp, id_var3, start, count, var3);
             NC_CHECK_STATUS(status);

             step++;
          }

        if (NC_NOERR != (status = nc_close (ncid)))
          {
             fprintf (stderr, "*** ERROR: closing file %s\n", outfile);
          }

        tstart = tend;
        ncid = -1;
     }

   fprintf (stdout,
            "Bounding box:\nlon: %9.4f, %9.4f => %2d (0.05 deg steps)\nlat: %9.4f, %9.4f => %2d\n",
            bbox[0], bbox[1], (int) ceil((bbox[1]-bbox[0])/0.05),
            bbox[2], bbox[3], (int) ceil((bbox[3]-bbox[2])/0.05));

   status = 0;
cleanup_and_exit:
   if (status)
     {
        fprintf (stderr, "*** ERROR: %s\n", nc_strerror(status));
     }

   FREE(var0);
   FREE(var1);
   FREE(var2);
   FREE(var3);
   FREE(var4);
   FREE(xtrack);
   FREE(ul_bitfield);
   FREE(ui_bitfield);
   FREE(us_bitfield);
   FREE(uc_bitfield);
   FREE( l_bitfield);
   FREE( i_bitfield);
   FREE( s_bitfield);
   FREE( c_bitfield);
   free_slit_grid (g);
   close_obs (&o);

   return status ? 1 : 0;
}
