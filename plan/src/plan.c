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
#include <wordexp.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include <libnovas.h>

#include "solar.h"
#include "scan.h"
#include "plan_list.h"
#include "scan_methods.h"
#include "vis.h"

#define DEFAULT_SCAN_METHOD_NAME "std"
#define DEFAULT_NUM_PLAN_DAYS    14
#define DEFAULT_NUM_SZA_DAYS     1

#define SMA_MAX_CALIBRATED_MIRROR_X  49600.0
#define SMA_MAX_CALIBRATED_MIRROR_Y   4400.0
/* The maximum range over which the scan mirror should be commanded is:
 *    |X| <= SMA_MAX_CALIBRATED_MIRROR_X,
 *    |Y| <= SMA_MAX_CALIBRATED_MIRROR_Y.
 * These coordinates refer to the scan mirror tilt angle in microradians.
 * (From TEMPO ConOps, Ball doc 2418231, Rev E, section 12.3.2, page 75,
 *  10/12/2017)
 */

#define SMA_MAX_SCAN_TABLE_STEP   100.0
/* The maximum mirror step size within a scan table [microradians]
 * (From TEMPO ConOps, Ball doc 2418231, Rev E, section 12.3.2, page 75,
 *  10/12/2017)
 */

#define IRRADIANCE_SUN_ANGLE_DEG 30.0

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

typedef struct
{
   const char *scan_tailoring_file;
   const char *vis_output_file;
   int num_sza_days;
}
Optional_Output_Type;

static void usage (void)
{
   fprintf (stderr, "Usage: plan [options]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -d | --date DATE         Plan start date:\n");
   fprintf (stderr, "                                DATE = (YYYY-MM-DD | DDDD days since the epoch)\n");
   fprintf (stderr, "   -n | --ndays N[,M]       N=number of days to plan [default=14]\n");
   fprintf (stderr, "                            M=number of days for SZA map output [default=1]\n");
   fprintf (stderr, "   -s | --scan METHOD       METHOD = std | opt1 [default=std]\n");
   fprintf (stderr, "   -t | --type SCAN_TYPE    Scan type [default=%d (TEMPO_SCAN_TYPE_STANDARD)]\n", TEMPO_SCAN_TYPE_STANDARD);
   fprintf (stderr, "   -N | --night             Enable night-lights scans\n");
   fprintf (stderr, "   -o | --output FILE       Radiance scan output file [default=stdout]\n");
   fprintf (stderr, "   -I | --irr FILE          Generate irradiance plan output file\n");
   fprintf (stderr, "   -m | --master FILE       Generate master scan table\n");
   fprintf (stderr, "   -z | --szaout FILE       Generate netCDF SZA map output to visualize\n");
   fprintf (stderr, "                            the solar illumination at the start of each scan\n\n");
   fprintf (stderr, "   -c | --config FILE       Configuration file\n");
   fprintf (stderr, "  For testing:\n");
   fprintf (stderr, "   -T | --tailor FILE       Output a 'nominal' scan tailoring file [default=none]\n");
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

static char *expand_string (const char *s)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        tell_verror (TELL_UNKNOWN_ERROR,
                     "%s: expanding path: %s", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
     }

   return s_exp;
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

   if (NULL == (eph->ephem_name = expand_string (ephem_name)))
     return -1;

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

static int print_standard_scan_table (FILE *fp, double roll_angle, double delta_x,
                                      double mirror_x0, double mirror_x1)
{
   const char fmt[] = "%0.3f,%0.7g,%0.4g\n";
   double cos_phi, sin_phi, x, dx_r;

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

static int timestamp_created (FILE *fp)
{
   time_t now_tt = time(NULL);
   return fprintf (fp, "# Created: %s", ctime(&now_tt));
}

typedef int Print_Method_Type(FILE *, double, double, double, double);

static int generate_scan_table (config_t *cfg, FILE *fp,
                               Print_Method_Type *print_method)
{
   double mirror_x0, mirror_x1, delta_x;
   double step_size, roll_angle, max_roll_angle;

   timestamp_created (fp);

   if (0 != read_master_scan_table_params (cfg, &step_size, &roll_angle))
     return -1;

   /* convert from microradians to radians */
   roll_angle *= 1.e-6;

   max_roll_angle = atan2 (SMA_MAX_CALIBRATED_MIRROR_Y,
                           SMA_MAX_CALIBRATED_MIRROR_X);

   if (fabs(roll_angle) > max_roll_angle)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: roll_angle = %0.6g exceeds allowed maximum (+/-%0.6g microradians)",
                     __func__, roll_angle*1.e6, max_roll_angle*1.e6);
        return -1;
     }

   mirror_x0 =  SMA_MAX_CALIBRATED_MIRROR_X;
   mirror_x1 = -SMA_MAX_CALIBRATED_MIRROR_X;
   delta_x = mirror_tilt (step_size);

   if (fabs(delta_x) > SMA_MAX_SCAN_TABLE_STEP)
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: step_size = %0.6g exceeds allowed maximum (+/-%0.6g microradians)",
                     __func__, delta_x, SMA_MAX_SCAN_TABLE_STEP);
        return -1;
     }

   (void) fprintf (fp, "# roll = %0.6g microradian\n", roll_angle * 1.e6);
   (void) fprintf (fp, "mirror_x,delta_x,mirror_y\n");

   return print_method (fp, roll_angle, delta_x, mirror_x0, mirror_x1);
}

static int generate_master_scan_table (config_t *cfg, const char *type, FILE *fp)
{
   (void) type;
   return generate_scan_table (cfg, fp, print_standard_scan_table);
}

static int read_epoch (config_t *cfg)
{
   config_setting_t *s;
   const char *epoch;

   if (NULL == (s = config_lookup (cfg, "sat_config")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing sat_config in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "epoch", &epoch))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,"%s: reading : %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   return tio_time_set_taix_epoch (epoch);
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

static int generate_scan_vis (config_t *cfg, const char *filename, int num_days,
                              Solar_Geom_Type *solar_geom,
                              const Plan_List_Type *plan_list,
                              double step_size,
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

   if (0 != sm->sm_vis (v, plan_list, step_size, num_days, ncid))
     goto return_status;

   tio_status = TIO_close (ncid);

   status = 0;
return_status:
   vis_free (v);
   return ((tio_status == 0) && (status == 0)) ? 0 : -1;
}

static Plan_List_Type *
attach_nightlights_scans (const Scan_Type *scan, Solar_Geom_Type *solar_geom,
                          Scan_Limit_Times_Type *limit_times,
                          Plan_List_Type *entry)
{
   const Scan_Method_Type *nl_dawn = find_scan_method ("nl_dawn");
   const Scan_Method_Type *nl_dusk = find_scan_method ("nl_dusk");
   Plan_List_Type *dawn = NULL;
   Plan_List_Type *dusk = NULL;
   Plan_List_Type *last = NULL;
   double radiance_scans_end_time;

   /* Because radiance scans are allowed to start a little early,
    * it may be necessary to end dawn scans earlier.
    */
   if (entry->tstart < limit_times->jd_utc_beg)
     {
        limit_times->jd_utc_beg = entry->tstart;
     }
   dawn = nl_dawn->sm_plan (scan, solar_geom, limit_times);
   if (dawn && (dawn->num_repeats > 0))
     {
        if (0 != plan_list_append (&dawn, entry))
          return NULL;
        entry = dawn;
     }

   /* Because radiance scans are allowed to run a little bit late,
    * it may be necessary to start dusk scans later.
    */
   for (last = entry; last != NULL; last = last->next)
     {
        if (last->next == NULL)
          break;
     }
   radiance_scans_end_time = last->tstart + last->num_repeats * last->scan_duration / SEC_PER_DAY;
   if (radiance_scans_end_time > limit_times->jd_utc_end)
     {
        limit_times->jd_utc_end = radiance_scans_end_time;
     }

   dusk = nl_dusk->sm_plan (scan, solar_geom, limit_times);
   if (dusk && (dusk->num_repeats > 0))
     {
        if (0 != plan_list_append (&entry, dusk))
          return NULL;
     }

   return entry;
}

static int print_scan_tailoring_file (Solar_Geom_Type *sgt, double jd_utc0, double jd_utc1,
                                      const char *filename)
{
   double delta = 15.0 * 60.0 / SEC_PER_DAY;
   double unix_epoch_jd;
   int i, n;
   FILE *fp;

   if (filename == NULL)
     return 0;

   jd_utc0 = floor (jd_utc0) + 0.5;
   jd_utc1 = ceil (jd_utc1) + 0.5;

   n = (jd_utc1 - jd_utc0) / delta;

   if (n <= 0)
     return 0;

   unix_epoch_jd = novas_julian_date (1970,1,1,0.0);

   if (NULL == (fp = fopen (filename, "w")))
     {
        fprintf (stderr, "*** Error: opening %s for writing\n", filename);
        return -1;
     }

   timestamp_created (fp);

   fprintf (fp, "# time: Time since the TEMPO epoch [TAI hours]\n");
   fprintf (fp, "# xoff: Mirror pointing offset in the East/West direction [degrees]\n");
   fprintf (fp, "# yoff: Mirror pointing offset in the North/South direction [degrees]\n");
   fprintf (fp, "# solar_boresight_angle: Angle between the sun and the TEMPO boresight [degrees]\n");
   fprintf (fp, "time,xoff,yoff,solar_boresight_angle\n");

   for (i = 0; i < n; i++)
     {
        double jd_utc = jd_utc0 + i * delta;
        double angle, t_utc, taix;

        if (0 != sgt->sgt_sat_sun_angles (sgt, jd_utc, &angle, NULL))
          return -1;

        t_utc = (jd_utc - unix_epoch_jd) * SEC_PER_DAY;

        if (0 != tio_time_utc_to_taix (t_utc, &taix))
          return -1;

        if (fprintf (fp, "%0.8f,%f,%f,%0.6f\n", taix/3600.0, 0.0, 0.0, angle) < 0)
          {
             fprintf (stderr, "*** Error: writing to %s\n", filename);
             break;
          }
     }

   return fclose (fp);
}

static int write_scan_plan (FILE *fp, const Ephem_Type *eph, const Solar_Geom_Type *solar_geom, const Scan_Type *scan,
                            const char *scan_method, const Plan_List_Type *plan_list)
{
   char epoch_str[32];

   /* Write out scan plan */
   if (0 != TIO_mktimestamp_str (0.0, 1, epoch_str, sizeof(epoch_str)))
     return -1;

   timestamp_created (fp);
   (void) fprintf (fp, "# %s = scan method\n", scan_method);
   if (0 != scan->st_print_params (scan, "#", fp))
     return -1;
   if (0 != solar_geom->sgt_print_params (solar_geom, "#", fp))
     return -1;
   (void) fprintf (fp, "# NOVAS ephemeris: %s\n", eph->ephem_name);
   (void) fprintf (fp, "# TEMPO epoch: %s\n", epoch_str);
   (void) fprintf (fp, "#\n");

   return plan_list_write (fp, mirror_tilt, plan_list);
}

static int write_irradiance_plan (FILE *fp, Solar_Geom_Type *solar_geom, const Cal_Date_Type *t0, int num_days)
{
   const char header[] = "time,solar_theta,solar_phi,timestamp\n";
   double unix_epoch_jd;
   double jd_utc0, jd_utc1, jd_utc;
   double irr_angle = IRRADIANCE_SUN_ANGLE_DEG;

   unix_epoch_jd = novas_julian_date (1970,1,1,0.0);

   jd_utc0 = novas_julian_date (t0->year, t0->month, t0->day, t0->hour);
   jd_utc1 = jd_utc0 + num_days;

   timestamp_created (fp);
   if (fprintf (fp, header) < 0)
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: fprintf failed", __func__);
        return -1;
     }

   for (jd_utc = jd_utc0; jd_utc < jd_utc1; jd_utc += 1.0)
     {
        char buf[32];
        double jd_utc_irr, solar_theta, solar_phi, tirr_utc, tirr_tai;

        if (0 != scan_irradiance_time (solar_geom, irr_angle, jd_utc, &jd_utc_irr))
          return -1;
        if (0 != solar_geom->sgt_sat_sun_angles (solar_geom, jd_utc_irr, &solar_theta, &solar_phi))
          return -1;

        if (0 != mkjdtimestr (jd_utc_irr, buf, sizeof(buf)))
          return -1;

        tirr_utc = (jd_utc_irr - unix_epoch_jd) * SEC_PER_DAY;

        if (0 != tio_time_utc_to_taix (tirr_utc, &tirr_tai))
          return -1;

        (void) fprintf (fp, "%0.3f,%0.3f,%0.3f,\"%s\"\n",
                        tirr_tai, solar_theta, solar_phi, buf);
     }

   return 0;
}

static Plan_List_Type *generate_scan_plan (const Ephem_Type *eph, Solar_Geom_Type *solar_geom,
                                           const Scan_Type *scan, const Scan_Method_Type *sm,
                                           const Cal_Date_Type *t0, int num_plan_days, int enable_night_scan,
                                           const Optional_Output_Type *oot)
{
   Plan_List_Type *plan_list = NULL;
   double jd_utc, jd_utc0, jd_utc1;
   int status = -1;

   jd_utc0 = novas_julian_date (t0->year, t0->month, t0->day, t0->hour);
   jd_utc1 = jd_utc0 + num_plan_days;

   if ((jd_utc0 < eph->jd_begin) || (eph->jd_end < jd_utc1))
     {
        fprintf (stderr, "*** Ephemeris JD=%0.2f-%0.2f doesn't cover JD=%0.2f + %d days \n",
                 eph->jd_begin, eph->jd_end,
                 jd_utc0, num_plan_days);
        return NULL;
     }

   if (0 != print_scan_tailoring_file (solar_geom, jd_utc0, jd_utc1, oot->scan_tailoring_file))
     return NULL;

   for (jd_utc = jd_utc0; jd_utc < jd_utc1; jd_utc += 1.0)
     {
        Scan_Limit_Times_Type limit_times = {0};
        Plan_List_Type *entry = NULL;

        if (0 != scan_limit_times (scan, jd_utc, solar_geom, &limit_times))
          goto return_status;

        if (NULL == (entry = sm->sm_plan (scan, solar_geom, &limit_times)))
          goto return_status;

        if (enable_night_scan)
          {
             if (NULL == (entry = attach_nightlights_scans (scan, solar_geom, &limit_times, entry)))
               goto return_status;
          }

        if (0 != plan_list_append (&plan_list, entry))
          goto return_status;
     }

   status = 0;
return_status:
   if (status)
     {
        plan_list_free (plan_list);
        plan_list = NULL;
     }

   return plan_list;
}

static void close_outfile (FILE *fp, const char *filename)
{
   if ((fp == NULL) || (fp == stdout) || (fp == stderr))
     return;

   if (0 != fclose (fp))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "closing file: %s\n",
                     filename ? filename : "<null>");
     }
}

static FILE *handle_outfile_arg (const char *arg, const char *mode)
{
   if (arg == NULL)
     return NULL;
   else if ((strlen(arg) == 1) && (arg[0] == '-'))
     return stdout;
   else
     return fopen (arg, mode);
}

int main (int argc, char **argv)
{
   const char appname[] = "plan";
   char *config_file = "plan.cfg";
   char *scan_outfile = NULL;
   char *irr_outfile = NULL;
   char *master_outfile = NULL;
   FILE *fp_scan = stdout;
   FILE *fp_master = NULL;
   FILE *fp_irr = NULL;
   char *scan_method = DEFAULT_SCAN_METHOD_NAME;
   uint16_t scan_type = TEMPO_SCAN_TYPE_STANDARD;
   int num_plan_days = DEFAULT_NUM_PLAN_DAYS;
   int status = EXIT_FAILURE;
   int enable_night_scan = 0;
   int have_date = 0;
   Cal_Date_Type t0 = {0};
   int ndays_since_epoch = 0;
   Optional_Output_Type oot =
     {
        .scan_tailoring_file = NULL,
        .vis_output_file = NULL,
        .num_sza_days = DEFAULT_NUM_SZA_DAYS,
     };
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"date",    required_argument, 0, 'd'},
        {"config",  required_argument, 0, 'c'},
        {"ndays",   required_argument, 0, 'n'},
        {"night",   no_argument,       0, 'N'},
        {"scan",    required_argument, 0, 's'},
        {"type",    required_argument, 0, 't'},
        {"output",  required_argument, 0, 'o'},
        {"irr",     required_argument, 0, 'I'},
        {"szaout",  required_argument, 0, 'z'},
        {"tailor",  required_argument, 0, 'T'},
        {"master",  no_argument, 0, 'm'},
        {0,0,0,0}
     };
   config_t cfg = {0};
   Ephem_Type eph = {0};
   Scan_Type *scan = NULL;
   Plan_List_Type *plan_list = NULL;
   Solar_Geom_Type *solar_geom = NULL;
   const Scan_Method_Type *sm = NULL;

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
        int c = getopt_long (argc, argv, "hNc:d:I:m:n:o:s:t:T:z:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: getopt returned character %d??",
                          __func__, c);
             usage();
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
             if (NULL == strchr(optarg, '-'))
               {
                  if (1 != sscanf (optarg, "%d", &ndays_since_epoch))
                    {
                       fprintf (stderr, "*** error reading date option %s\n", optarg);
                       usage();
                    }
               }
             else if (3 != sscanf (optarg, "%hd%*c%hd%*c%hd", &t0.year, &t0.month, &t0.day))
               {
                  fprintf (stderr, "*** error reading date option %s\n", optarg);
                  usage();
               }
             have_date++;
             break;
           case 'h':
             usage();
             break;
           case 'N':
             enable_night_scan++;
             break;
           case 'm':
             master_outfile = optarg;
             if (NULL == (fp_master = handle_outfile_arg (master_outfile, "w")))
               {
                  fprintf (stderr, "*** error opening file %s for writing\n", master_outfile);
                  goto return_status;
               }
             break;
           case 'I':
             irr_outfile = optarg;
             if (NULL == (fp_irr = handle_outfile_arg (irr_outfile, "w")))
               {
                  fprintf (stderr, "*** error opening file %s for writing\n", irr_outfile);
                  goto return_status;
               }
             break;
           case 'n':
             if (NULL != strchr (optarg, ','))
               {
                  if (2 != sscanf (optarg, "%d,%d", &num_plan_days, &oot.num_sza_days))
                    usage ();
               }
             else
               {
                  if (1 != sscanf (optarg, "%d", &num_plan_days))
                    usage ();
               }
             break;
           case 'o':
             scan_outfile = optarg;
             if (NULL == (fp_scan = handle_outfile_arg (scan_outfile, "w")))
               {
                  fprintf (stderr, "*** error opening file %s for writing\n", scan_outfile);
                  goto return_status;
               }
             break;
           case 's':
             scan_method = optarg;
             break;
           case 't':
             if (1 != sscanf (optarg, "%hd", &scan_type))
               {
                  fprintf (stderr, "*** error reading scan_type: %s\n", optarg);
                  goto return_status;
               }
             if (scan_type == TEMPO_SCAN_TYPE_NIGHTLIGHTS)
               {
                  fprintf (stderr, "*** error: scan_type=%d is not user-selectable\n", scan_type);
                  goto return_status;
               }
             break;
           case 'T':
             oot.scan_tailoring_file = optarg;
             break;
           case 'z':
             oot.vis_output_file = optarg;
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

   if (fp_master)
     {
        if (0 != generate_master_scan_table (&cfg, NULL, fp_master))
          goto return_status;
     }

   if (have_date == 0)
     {
        fprintf (stderr, "Usage error: plan start date not specified (--date option missing)\n");
        goto return_status;
     }

   if (0 != read_epoch (&cfg))
     goto return_status;

   if (ndays_since_epoch > 0)
     {
        double taix = SEC_PER_DAY * ((double) ndays_since_epoch);
        int year, month, day;
        if (0 != tio_time_taix_to_utc_caldate (taix, &year, &month, &day, &t0.hour))
          goto return_status;
        t0.year = year;
        t0.month = month;
        t0.day = day;
     }

   /* satellite orbital station determines effective time zone */
   if (0 != read_sat_time_zone (&cfg, &t0.hour))
     goto return_status;

   if ((0 != ephem_open (&cfg, &eph))
       || (NULL == (solar_geom = solar_geom_init (&cfg))))
     goto return_status;

   if (fp_irr)
     {
        if (0 != write_irradiance_plan (fp_irr, solar_geom, &t0, num_plan_days))
          goto return_status;
     }

   if (NULL == (scan = scan_open (&cfg, scan_type)))
     goto return_status;

   if (NULL == (sm = find_scan_method (scan_method)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: unrecognized scan method '%s'",
                     __func__, scan_method ? scan_method : "(null)");
        goto return_status;
     }

   if (NULL == (plan_list = generate_scan_plan (&eph, solar_geom, scan, sm, &t0, num_plan_days, enable_night_scan, &oot)))
     goto return_status;

   if (0 != write_scan_plan (fp_scan, &eph, solar_geom, scan, scan_method, plan_list))
     goto return_status;

   /* Optionally, generate some plots */
   if (0 != generate_scan_vis (&cfg, oot.vis_output_file, oot.num_sza_days, solar_geom, plan_list,
                               scan->st_step_size(scan), sm))
     goto return_status;

   status = EXIT_SUCCESS;
return_status:
   if (solar_geom) solar_geom->sgt_delete (solar_geom);
   if (scan) scan->st_delete (scan);
   (void) ephem_close (&eph);
   plan_list_free (plan_list);
   close_outfile (fp_scan, scan_outfile);
   close_outfile (fp_master, master_outfile);
   close_outfile (fp_irr, irr_outfile);
   config_destroy (&cfg);
   tell_close ();

   return status;
}
