/** @file plan.c
 *  @brief Main program
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>

#include <libnovas.h>

#include "solar.h"
#include "scan.h"
#include "plan_list.h"
#include "scan_methods.h"
#include "vis.h"

#define DEFAULT_SCAN_METHOD_NAME "std"
#define DEFAULT_NUM_PLAN_DAYS    14

#define MAX_ALLOWED_AZIMUTH_URAD  1.1e5
/* The maximum SMA coordinate allowed by the flight software is
 * 55000 microradians.  This refers to the mirror tilt angle,
 * and corresponds to a 1.1e5 urad azimuthal offset of the line
 * of sight in the field of regard.  This is the maximum angle
 * allowed by the hardware.
 */

typedef struct
{
   char *ephem_name;
   double jd_begin;
   double jd_end;
   short int de_number;
}
Ephem_Type;

typedef struct
{
   short int year;
   short int month;
   short int day;
   double hour;
}
Cal_Date_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: plan [options]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -d | --date YYYY-MM-DD   Plan start date\n");
   fprintf (stderr, "   -n | --ndays NUM         Number of days to plan [default=14]\n");
   fprintf (stderr, "   -o | --output FILE       Output file [default=stdout]\n");
   fprintf (stderr, "   -s | --scan METHOD       METHOD = std | opt1 [default=std]\n");
   fprintf (stderr, "   -c | --config FILE       Configuration file\n");
   fprintf (stderr, "   -m | --master            Generate master scan table, then exit\n");
   fprintf (stderr, "   -z | --szaout FILE       Generate netCDF output to help visualize\n");
   fprintf (stderr, "                            the solar zenith angle at selected times\n");
   exit (EXIT_SUCCESS);
}

static int read_config_file (const char *config_file, config_t *cfg)
{
   if (0 == config_read_file (cfg, config_file))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: Reading %s:%d - %s",
                     __func__, config_error_file(cfg),
                     config_error_line(cfg), config_error_text(cfg));
        return -1;
     }

   return 0;
}

int mkjdtimestr (double jd_utc, char *buf, int bufsize)
{
   Cal_Date_Type t0;
   double sec;
   int hour, min, status;

   novas_cal_date (jd_utc, &t0.year, &t0.month, &t0.day, &t0.hour);
   hour = t0.hour;
   min  = 60*(t0.hour - hour);
   sec  = 3600*(t0.hour - hour - min/60.0);

   if (((status = snprintf (buf, bufsize, "%4d-%02d-%02dT%02d:%02d:%06.3fZ",
                            t0.year, t0.month, t0.day, hour, min, sec)) < 0)
       || (status >= bufsize))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed (status=%d)",
                     __func__, status);
        return -1;
     }

   return 0;
}

static int ephem_close (Ephem_Type *eph)
{
   short int error = 0;

   FREE(eph->ephem_name);
   memset ((char *)eph, 0, sizeof (*eph));

   if ((error = novas_ephem_close ()) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: error %d while closing ephemeris",
                     __func__, error);
     }

   return error ? -1 : 0;
}

static int ephem_open (config_t *cfg, Ephem_Type *eph)
{
   config_setting_t *s;
   const char *ephem_name;
   short int error;

   if (NULL == (s = config_lookup (cfg, "novas_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing novas_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "ephem_name", &ephem_name))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading ephem_name: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (eph->ephem_name = strdup (ephem_name)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
        return -1;
     }

   if ((error = novas_ephem_open (eph->ephem_name,
                                  &eph->jd_begin, &eph->jd_end,
                                  &eph->de_number)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: error %d while opening ephemeris %s",
                     __func__, error, eph->ephem_name);
        return -1;
     }

   return 0;
}

static int
read_master_scan_table_params
(config_t *cfg, double *pstep_size, double *proll_angle)
{
   config_setting_t *s;
   double step_size, roll_angle;

   if (NULL == (s = config_lookup (cfg, "master_table_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing master_table_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if ((CONFIG_TRUE != config_setting_lookup_float (s, "step_size", &step_size))
       || (CONFIG_TRUE != config_setting_lookup_float (s, "roll_angle", &roll_angle))
      )
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading master_table_config: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (pstep_size) *pstep_size = step_size;
   if (proll_angle) *proll_angle = roll_angle;

   return 0;
}

static double mirror_tilt (double azimuth)
{
   /* The FOR coordinate system refers to the azimuth and elevation angular
    * coordinates in the field of regard indicating the line of sight
    * from which we want to collect photons entering through the slit.
    * To command the instrument, the flight software wants the mirror tilt
    * angle needed to access that line of sight.  Using the law of reflection,
    * the mirror tilt angle is half the azimuthal angle coordinate in the
    * field of regard.
    *
    * The azimuthal angle coordinate in the field of regard increases toward
    * the east (+X axis in a right-handed coordinate system).  The elevation
    * coordinate increases toward the south (+Y axis).
    *
    * The C&THB documentation for the SMA_MOVE command says:
    *
    * "Neglecting alignment tolerances, a motion of the scan mirror in
    * the positive X-direction moves the line of sight in the Spacecraft
    * +X direction (East) . A motion of the scan mirror in the positive
    * Y-direction moves the line of sight in the Spacecraft +Y direction
    * (South)."
    *
    * Therefore, the mirror tilt angle +X coordinate has the same sign
    * as the +X azimuthal angle in the field of regard.
    *
    * Both angles are in microradians.
    */

   return 0.5 * azimuth;
}

static int generate_master_scan_table (config_t *cfg, FILE *fp)
{
   const char fmt[] = "%0.3f,%0.3f,%0.3f\n";
   double mirror_x0, mirror_x1, delta_x, mirror_y, dx_r;
   double xstart, step_size, roll_angle;
   double cos_phi, sin_phi, x;

   if (0 != read_master_scan_table_params (cfg, &step_size, &roll_angle))
     return -1;

   (void) fprintf (fp, "# roll = %0.3f microradian\n", roll_angle);

   /* convert from microradians to radians */
   roll_angle *= 1.e-6;

   xstart = MAX_ALLOWED_AZIMUTH_URAD;
   mirror_x0 = mirror_tilt (xstart);
   mirror_x1 = mirror_tilt (-xstart);
   delta_x = mirror_tilt (step_size);

   (void) fprintf (fp, "mirror_x,delta_x,mirror_y\n");

   if (roll_angle == 0.0)
     {
        mirror_y = 0.0;

        (void) fprintf (fp, fmt, mirror_x1, delta_x, mirror_y);
        (void) fprintf (fp, fmt, mirror_x0, delta_x, mirror_y);
        return 0;
     }

   /* +X is eastward, +Y is southward,  +Z is along the boresight,
    * at (X,Y) = (0,0).
    * +roll_angle is clockwise rotation of a vector
    * about the +Z axis of the right-handed coordinate system.
    * For example, consider the vector (1,0).  Applying roll_angle>0
    * will rotate the vector CW, making the Y coordinate > 0, and
    * making the X coordinate smaller.
    * Applying roll_angle=pi/2 will rotate the vector into (0,1).
    *
    * Vector rotation:
    *  xp = cos(phi) * x - sin(phi) * y
    *  yp = sin(phi) * x + cos(phi) * y
    */

   cos_phi = cos (roll_angle);
   sin_phi = sin (roll_angle);

   /* initially, y=0 and  delta_y=0, so rotation yields
    * x, y components as follows: */

   dx_r = cos_phi * delta_x;

   x = mirror_x1;
   (void) fprintf (fp, fmt, cos_phi*x, dx_r, sin_phi*x);
   x = mirror_x0;
   (void) fprintf (fp, fmt, cos_phi*x, dx_r, sin_phi*x);

   return 0;
}

static int read_sat_time_zone (config_t *cfg, double *hour)
{
   config_setting_t *s;
   double sat_lon;

   if (NULL == (s = config_lookup (cfg, "sat_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing sat_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_float (s, "sat_lon", &sat_lon))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading sat_lon: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   /* 15 degrees per hour */
   *hour = (-sat_lon)/15.0;

   return 0;
}

static int generate_vis (config_t *cfg, const char *filename,
                         Solar_Geom_Type *solar_geom,
                         const Plan_List_Type *plan_list,
                         const Scan_Method_Type *sm)
{
   Vis_Type *v = NULL;
   int ncid, tio_status, status = -1;

   if (filename == NULL)
     return 0;

   if (NULL == (v = vis_init (cfg, solar_geom)))
     return -1;

   if (0 != TIO_create (filename, NC_NETCDF4, &ncid))
     goto return_status;

   if (0 != sm->sm_vis (v, plan_list, ncid))
     goto return_status;

   tio_status = TIO_close (ncid);

   status = 0;
return_status:
   vis_free (v);
   return ((tio_status == 0) && (status == 0)) ? 0 : -1;
}

static int generate_plan (config_t *cfg, const Cal_Date_Type *t0,
                          int num_plan_days, const char *scan_method,
                          const char *vis_output_file, FILE *fp)
{
   Ephem_Type eph = {0};
   Scan_Type *scan = NULL;
   Plan_List_Type *plan_list = NULL;
   const Scan_Method_Type *sm = NULL;
   Solar_Geom_Type *solar_geom = NULL;
   double jd_utc, jd_utc0, jd_utc1;
   int status = -1;

   jd_utc0 = novas_julian_date (t0->year, t0->month, t0->day, t0->hour);
   jd_utc1 = jd_utc0 + num_plan_days;

   if (0 != ephem_open (cfg, &eph))
     return -1;

   if ((jd_utc0 < eph.jd_begin) || (eph.jd_end < jd_utc1))
     {
        fprintf (stderr, "*** Ephemeris JD=%0.2f-%0.2f doesn't cover JD=%0.2f + %d days \n",
                 eph.jd_begin, eph.jd_end,
                 jd_utc0, num_plan_days);
        goto return_status;
     }

   if (NULL == (sm = find_scan_method (scan_method)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unrecognized scan method '%s'",
                     __func__, scan_method ? scan_method : "(null)");
        goto return_status;
     }

   if (NULL == (solar_geom = solar_geom_init (cfg)))
     goto return_status;

   if (NULL == (scan = scan_open (cfg)))
     goto return_status;

   for (jd_utc = jd_utc0; jd_utc < jd_utc1; jd_utc += 1.0)
     {
        Scan_Limit_Times_Type limit_times = {0};
        Plan_List_Type *entry = NULL;

        if (0 != scan_limit_times (scan, jd_utc, solar_geom, &limit_times))
          goto return_status;

        if (NULL == (entry = sm->sm_plan (scan, solar_geom, &limit_times)))
          goto return_status;

        if (0 != plan_list_append (&plan_list, entry))
          goto return_status;
     }

   (void) fprintf (fp, "# %s = scan method\n", scan_method);
   if (0 != scan->st_print_params (scan, "#", fp))
     goto return_status;
   if (0 != solar_geom->sgt_print_params (solar_geom, "#", fp))
     goto return_status;
   (void) fprintf (fp, "# NOVAS ephemeris: %s\n", eph.ephem_name);
   (void) fprintf (fp, "#\n");

   if (0 != plan_list_write (fp, mirror_tilt, plan_list))
     goto return_status;

   if (0 != generate_vis (cfg, vis_output_file, solar_geom, plan_list, sm))
     goto return_status;

   status = 0;
return_status:
   (void) ephem_close (&eph);
   if (solar_geom)
     {
        solar_geom->sgt_delete (solar_geom);
        solar_geom = NULL;
     }
   if (scan)
     {
        scan->st_delete (scan);
        scan = NULL;
     }
   plan_list_free (plan_list);

   return status;
}

int main (int argc, char **argv)
{
   const char appname[] = "plan";
   char *config_file = "plan.cfg";
   char *outfile = NULL;
   FILE *fp = stdout;
   char *scan_method = DEFAULT_SCAN_METHOD_NAME;
   int num_plan_days = DEFAULT_NUM_PLAN_DAYS;
   int status = EXIT_FAILURE;
   int do_master_scan_table = 0;
   int have_date = 0;
   const char *vis_output_file = NULL;
   Cal_Date_Type t0 = {0};
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"date",    required_argument, 0, 'd'},
        {"config",  required_argument, 0, 'c'},
        {"ndays",   required_argument, 0, 'n'},
        {"scan",    required_argument, 0, 's'},
        {"output",  required_argument, 0, 'o'},
        {"szaout",  required_argument, 0, 'z'},
        {"master",  no_argument, 0, 'm'},
        {0,0,0,0}
     };
   config_t cfg;

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);

   config_init (&cfg);

   /* Try reading the default config file, but if it doesn't exist,
    * keep going in case there's a config file on the command line */
   if (0 == access (config_file, F_OK | R_OK))
     {
        if (-1 == read_config_file (config_file, &cfg))
          goto return_status;
     }

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hmc:d:n:o:s:z:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: getopt returned character %d??",
                          __func__, c);
             goto return_status;
             break;
           case 'c':
             /* This config file will override the default one
              * that might have been read previously.
              * Subsequent command-line args will override
              * any corresponding config file values */
             config_file = optarg;
             if (-1 == read_config_file (config_file, &cfg))
               goto return_status;
             break;
           case 'd':
             if (3 != sscanf (optarg, "%hd%*c%hd%*c%hd", &t0.year, &t0.month, &t0.day))
               {
                  fprintf (stderr, "*** error reading date option %s\n", optarg);
                  usage();
               }
             have_date++;
             break;
           case 'h':
             usage();
             break;
           case 'm':
             do_master_scan_table++;
             break;
           case 'n':
             if (1 != sscanf (optarg, "%d", &num_plan_days))
               usage ();
             break;
           case 'o':
             outfile = optarg;
             if (NULL == (fp = fopen (outfile, "w")))
               {
                  fprintf (stderr, "*** error opening file %s for writing\n", outfile);
                  goto return_status;
               }
             break;
           case 's':
             scan_method = optarg;
             break;
           case 'z':
             vis_output_file = optarg;
             break;
          }
     }

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   {
      time_t now_tt = time(NULL);
      if (fprintf (fp, "# Created: %s", ctime(&now_tt)) < 0)
        goto return_status;
   }

   if (do_master_scan_table)
     {
        if (0 != generate_master_scan_table (&cfg, fp))
          goto return_status;
     }
   else
     {
        if (have_date == 0)
          {
             fprintf (stderr, "Usage error: plan start date not specified (--date option missing)\n");
             goto return_status;
          }

        /* satellite orbital station determines effective time zone */
        if (0 != read_sat_time_zone (&cfg, &t0.hour))
          goto return_status;

        if (0 != generate_plan (&cfg, &t0, num_plan_days, scan_method,
                                vis_output_file, fp))
          goto return_status;
     }

   status = EXIT_SUCCESS;
return_status:
   if ((fp != NULL) && (outfile != NULL))
     {
        if (0 != fclose (fp))
          {
             tell_verror (TELL_IO_WRITE_ERROR, "closing file %s\n", outfile);
          }
     }
   config_destroy (&cfg);
   tell_close ();

   return status;
}
