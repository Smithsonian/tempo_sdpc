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
#include <ioclib.h>

#include <tio.h>
#include <tio_template.h>

#include <libnovas.h>

#include "solar.h"
#include "scan.h"
#include "plan_list.h"
#include "scan_methods.h"
#include "vis.h"

int Plan_Verbose;

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

#define MIN_SCAN_DURATION_SEC (2 * 60.0)

static double __Unix_Epoch_JD;

enum
{
   PARTIAL_SCAN_END = 0,
   PARTIAL_SCAN_BEGIN = 1
};

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
   const char *vis_output_file;
   int num_sza_days;
}
Optional_Output_Type;

static double get_unix_epoch_jd (void)
{
   return __Unix_Epoch_JD;
}

static void set_unix_epoch_jd (void)
{
   __Unix_Epoch_JD = novas_julian_date (1970,1,1,0.0);
}

static void usage (void)
{
   fprintf (stderr, "Usage: plan [options]\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -h | --help              Print this usage message\n");
   fprintf (stderr, "   -d | --date DATE         Plan start date:\n");
   fprintf (stderr, "                                DATE = (YYYY-MM-DD | DDDD days since the epoch)\n");
   fprintf (stderr, "   -n | --ndays N[,M]       N=number of days to plan [default=14]\n");
   fprintf (stderr, "                            M=number of days for SZA map output [default=1]\n");
   fprintf (stderr, "   -s | --scan METHOD       METHOD = std | opt1 | split-METHOD-NAME [default=std]\n");
   fprintf (stderr, "                                e.g. split-opt1-CA, where CA is a setting in the config file\n");
   fprintf (stderr, "   -t | --type SCAN_TYPE    Scan type [default=%d (TEMPO_SCAN_TYPE_STANDARD)]\n", TEMPO_SCAN_TYPE_STANDARD);
   fprintf (stderr, "   -N | --nightlights       Enable night-lights scans\n");
   fprintf (stderr, "   -M | --maneuver SOURCE   Read maneuver windows from SOURCE.\n");
   fprintf (stderr, "                                SOURCE may be a file or a directory.\n");
   fprintf (stderr, "                                SOURCE=@FILE indicates FILE contains a list of files.\n");
   fprintf (stderr, "   -o | --output FILE       Radiance scan output file [default=stdout]\n");
   fprintf (stderr, "   -I | --irr FILE          Generate irradiance plan output file\n");
   fprintf (stderr, "   -m | --master FILE       Generate master scan table\n");
   fprintf (stderr, "   -z | --szaout FILE       Generate netCDF SZA map output to visualize\n");
   fprintf (stderr, "                            the solar illumination at the start of each scan\n\n");
   fprintf (stderr, "   -c | --config FILE       Configuration file\n");
   fprintf (stderr, "   -v | --verbose           Increase verbosity\n");
   fprintf (stderr, "  For testing:\n");
   fprintf (stderr, "   -T | --tailor FILE       Output ONLY a 'nominal' scan tailoring file\n");
   fprintf (stderr, "   -Z | --Zenith STR        Check solar zenith angle calculation\n");
   fprintf (stderr, "                            where STR=\"LON,LAT,YYYY-MM-DDTHH:MM:SSZ\"\n");
   fprintf (stderr, "                            and LON,LAT are in degrees\n");
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
   double fhour, sec;
   int hour, min, status;

   /* Try to avoid printing silly things like: 16:59:60Z or 16:60:00Z */

   novas_cal_date (jd_utc, &t0.year, &t0.month, &t0.day, &t0.hour);
   /* round at msec resolution */
   fhour = round (t0.hour * 3600.0e3)/3600.0e3;
   hour =   fhour;
   min  =  (fhour - hour)*60;
   sec  = ((fhour - hour)*60 - min)*60;
   /* truncate at msec resolution */
   sec = round (sec * 1.e3) / 1.e3;

   /* Seconds should be printed with no more precision than the above calculation,
    * e.g. compute to msec, and print msec */
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
attach_twilight_scans (const Scan_Type *scan, Solar_Geom_Type *solar_geom,
                       Scan_Limit_Times_Type *limit_times,
                       Twilight_Scan_Type *twilight_scan,
                       Plan_List_Type *entry)
{
   const Scan_Method_Type *twl_dawn = find_scan_method ("twilight_dawn");
   const Scan_Method_Type *twl_dusk = find_scan_method ("twilight_dusk");
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
   dawn = twl_dawn->sm_plan (scan, solar_geom, limit_times, twilight_scan);
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

   dusk = twl_dusk->sm_plan (scan, solar_geom, limit_times, twilight_scan);
   if (dusk && (dusk->num_repeats > 0))
     {
        if (0 != plan_list_append (&entry, dusk))
          return NULL;
     }

   return entry;
}

static int print_scan_tailoring_file (Solar_Geom_Type *sgt, const Cal_Date_Type *t0,
                                      int num_plan_days, const char *filename)
{
   double delta = 15.0 * 60.0 / SEC_PER_DAY;
   double unix_epoch_jd = get_unix_epoch_jd();
   double jd_utc0, jd_utc1;
   char epoch[32];
   int i, n;
   FILE *fp;

   if (filename == NULL)
     return -1;

   jd_utc0 = novas_julian_date (t0->year, t0->month, t0->day, t0->hour);
   jd_utc1 = jd_utc0 + num_plan_days;

   jd_utc0 = floor (jd_utc0) + 0.5;
   jd_utc1 = ceil (jd_utc1) + 0.5;

   n = (jd_utc1 - jd_utc0) / delta;

   if (n <= 0)
     return 0;

   if (0 != TIO_mktimestamp_str (0.0, 1, epoch, sizeof(epoch)))
     return -1;

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
   fprintf (fp, "# Epoch = %s\n", epoch);
   fprintf (fp, "time,xoff,yoff,solar_boresight_angle\n");

   for (i = 0; i < n; i++)
     {
        double jd_utc = jd_utc0 + i * delta;
        double angle, t_utc, taix;

        if (0 != sgt->sgt_sat_sun_position (sgt, jd_utc, &angle, NULL, NULL))
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

static int violates_sun_angle_constraint (const Scan_Type *st, Solar_Geom_Type *solar_geom,
                                          double jd_utc)
{
   double unix_epoch_jd = get_unix_epoch_jd();
   double min_sun_angle = st->st_min_sun_angle (st);
   double angle, t_utc, taix;

   if (0 != solar_geom->sgt_sat_sun_position (solar_geom, jd_utc, &angle, NULL, NULL))
     return -1;

   t_utc = (jd_utc - unix_epoch_jd) * SEC_PER_DAY;
   if (0 != tio_time_utc_to_taix (t_utc, &taix))
     return -1;

   if (angle < min_sun_angle)
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: scan violates sun-angle safety constraint: jd_utc = %f days (taix = %f sec since epoch) angle = %f deg (min_sun_angle = %f deg)",
                     __func__, jd_utc, taix, angle, min_sun_angle);
        return 1;
     }

   return 0;
}

static int verify_safety_constraints (Solar_Geom_Type *solar_geom, const Scan_Type *st,
                                      const Plan_List_Type *plan_list)
{
   const Plan_List_Type *entry;

   for (entry = plan_list; entry != NULL; entry = entry->next)
     {
        double jd_utc_start = entry->tstart;
        double jd_utc_end = entry->tstart + entry->num_repeats * entry->scan_duration / SEC_PER_DAY;

        if ((0 != violates_sun_angle_constraint (st, solar_geom, jd_utc_start))
            || (0 != violates_sun_angle_constraint (st, solar_geom, jd_utc_end)))
          return -1;
     }

   return 0;
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
   double unix_epoch_jd = get_unix_epoch_jd();
   double jd_utc0, jd_utc1, jd_utc;
   double irr_angle = IRRADIANCE_SUN_ANGLE_DEG;

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
        if (0 != solar_geom->sgt_sat_sun_position (solar_geom, jd_utc_irr, &solar_theta, &solar_phi, NULL))
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
                                           const Cal_Date_Type *t0, int num_plan_days,
                                           Twilight_Scan_Type *twilight_scan,
                                           Split_Scan_Type *split_scan)
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

   for (jd_utc = jd_utc0; jd_utc < jd_utc1; jd_utc += 1.0)
     {
        Scan_Limit_Times_Type limit_times = {0};
        Plan_List_Type *entry = NULL;

        if (0 != scan_limit_times (scan, jd_utc, solar_geom, &limit_times))
          goto return_status;

        if (NULL == (entry = sm->sm_plan (scan, solar_geom, &limit_times, split_scan)))
          goto return_status;

        if (twilight_scan)
          {
             if (NULL == (entry = attach_twilight_scans (scan, solar_geom, &limit_times, twilight_scan, entry)))
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

static int partial_scan (const Plan_List_Type *entry, int is_start, double t_bound, Plan_List_Type **pnew_entry)
{
   Plan_List_Type *new_entry = NULL;
   double scan_duration_days = entry->scan_duration / SEC_PER_DAY;
   double entry_end = entry->tstart + entry->num_repeats * scan_duration_days;
   double fnum = (entry_end - t_bound) / scan_duration_days;
   int num = (int) fnum;
   double frac, duration_sec, tstart;

   *pnew_entry = NULL;

   if (is_start)
     {
        frac = (fnum-num);
        duration_sec = floor (frac * entry->scan_duration);
        tstart = t_bound;
     }
   else
     {
        frac = 1.0 - (fnum-num);
        duration_sec = floor (frac * entry->scan_duration);
        tstart = t_bound - duration_sec / SEC_PER_DAY;
     }

   if (duration_sec < MIN_SCAN_DURATION_SEC)
     return 0;

   if (NULL == (new_entry = plan_list_entry_alloc (entry->scan_type)))
     return -1;

   *new_entry = *entry;

   /* Note that the partial scan always begins from xstart,ystart */
   new_entry->tstart = tstart;
   new_entry->scan_duration = duration_sec;
   new_entry->num_steps = frac * entry->num_steps;
   new_entry->num_repeats = 1;

   *pnew_entry = new_entry;

   return 1;
}

static int insert_maneuver_gap (Plan_List_Type *plan_list, double mnv_beg, double mnv_end)
{
   Plan_List_Type *entry = NULL;
   Plan_List_Type *parent_entry = NULL;
   Plan_List_Type *pre_gap_partial = NULL;
   Plan_List_Type *post_gap_partial = NULL;

   for (entry = plan_list; entry != NULL; parent_entry = entry, entry = entry->next)
     {
        double scan_duration_days = entry->scan_duration / SEC_PER_DAY;
        double entry_beg = entry->tstart;
        double entry_end = entry->tstart + entry->num_repeats * scan_duration_days;
        int orig_num_repeats, num_remaining, prev_num, need_partial_scan;

        /* Skip entries that precede the maneuver interval */
        if (entry_end < mnv_beg)
          continue;
        /* When remaining plan entries follow the maneuver, we're done */
        if (mnv_end < entry_beg)
          return 0;

        if (Plan_Verbose)
          {
             fprintf (stderr, "maneuver window=(%f,%f) overlaps radiance scan=(%f,%f)\n",
                      mnv_beg, mnv_end, entry_beg, entry_end);
          }

        if ((mnv_beg < entry_beg) && (entry_end < mnv_end))
          {
             /* Entire plan entry is lost */
             entry->num_repeats = 0;
          }
        else if ((entry_beg < mnv_beg) && (mnv_end < entry_end))
          {
             Plan_List_Type *curr = entry;
             Plan_List_Type *save_next = entry->next;
             Plan_List_Type *post_gap;

             /* Original plan entry may be split into 4 entries:
              *  1. pre-gap full scans
              *  2. pre-gap partial scan
              *   >>>> the maneuver <<<<
              *  3. post-gap partial scan
              *  4. post-gap full scans
              */
             orig_num_repeats = entry->num_repeats;

             /* pre-gap partial scan */
             if ((need_partial_scan = partial_scan (entry, PARTIAL_SCAN_END, mnv_beg, &pre_gap_partial)) < 0)
               return -1;
             if (need_partial_scan)
               {
                  curr->next = pre_gap_partial;
                  curr = pre_gap_partial;
               }

             /* post-gap partial scan */
             if ((need_partial_scan = partial_scan (entry, PARTIAL_SCAN_BEGIN, mnv_end, &post_gap_partial)) < 0)
               return -1;
             if (need_partial_scan)
               {
                  curr->next = post_gap_partial;
                  curr = post_gap_partial;
               }

             /* pre-gap full scans (zero is ok) */
             entry->num_repeats = floor((mnv_beg - entry_beg) / scan_duration_days);

             /* post-gap full scans (zero is ok) */
             num_remaining = floor((entry_end - mnv_end) / scan_duration_days);
             if (num_remaining)
               {
                  if (NULL == (post_gap = plan_list_entry_alloc (entry->scan_type)))
                    return -1;

                  curr->next = post_gap;
                  curr = post_gap;

                  *post_gap = *entry;

                  prev_num = orig_num_repeats - num_remaining;
                  post_gap->tstart = entry->tstart + prev_num * scan_duration_days;
                  post_gap->num_repeats = num_remaining;
               }

             curr->next = save_next;
             entry = curr;
          }
        else if (entry_beg < mnv_beg)
          {
             /* pre-gap partial scan */
             if ((need_partial_scan = partial_scan (entry, PARTIAL_SCAN_END, mnv_beg, &pre_gap_partial)) < 0)
               return -1;
             if (need_partial_scan)
               {
                  pre_gap_partial->next = entry->next;
                  entry->next = pre_gap_partial;
               }

             /* pre-gap full scans */
             entry->num_repeats = floor((mnv_beg - entry_beg) / scan_duration_days);
          }
        else
          {
             /* post-gap partial scan */
             if (parent_entry)
               {
                  if ((need_partial_scan = partial_scan (entry, PARTIAL_SCAN_BEGIN, mnv_end, &post_gap_partial)) < 0)
                    return -1;

                  if (need_partial_scan)
                    {
                       parent_entry->next = post_gap_partial;
                       post_gap_partial->next = entry;
                    }
               }

             /* post-gap full scans */
             num_remaining = floor ((entry_end - mnv_end) / scan_duration_days);
             prev_num = entry->num_repeats - num_remaining;
             entry->tstart = entry->tstart + prev_num * scan_duration_days;
             entry->num_repeats = num_remaining;
          }
     }

   return 0;
}

static int utcstr_to_jd_utc (const char *utc_str, double *jd_utc)
{
   Cal_Date_Type t0 = {0};
   double taix_sec, hour;
   int year, month, day;

   if ((0 != tio_time_utcstr_to_taix (utc_str, &taix_sec))
       || (0 != tio_time_taix_to_utc_caldate (taix_sec, &year, &month, &day, &hour)))
     return -1;

   t0.year = year;
   t0.month = month;
   t0.day = day;
   t0.hour = hour;
   *jd_utc = novas_julian_date (t0.year, t0.month, t0.day, t0.hour);

   return 0;
}

static int utc_to_jd_utc (double utc_time, double *jd_utc)
{
   Cal_Date_Type t0 = {0};
   double taix_sec, hour;
   int year, month, day;

   if ((0 != tio_time_utc_to_taix (utc_time, &taix_sec))
       || (0 != tio_time_taix_to_utc_caldate (taix_sec, &year, &month, &day, &hour)))
     return -1;

   t0.year = year;
   t0.month = month;
   t0.day = day;
   t0.hour = hour;
   *jd_utc = novas_julian_date (t0.year, t0.month, t0.day, t0.hour);

   return 0;
}

typedef struct Maneuver_Window_Type Maneuver_Window_Type;
struct Maneuver_Window_Type
{
   Maneuver_Window_Type *next;
   double beg_timet;
   double end_timet;
};

static void free_maneuver_window (Maneuver_Window_Type *mw)
{
   if (mw == NULL)
     return;
   FREE (mw);
}

static void free_maneuver_window_list (Maneuver_Window_Type *lst)
{
   while (lst)
     {
        Maneuver_Window_Type *next = lst->next;
        free_maneuver_window (lst);
        lst = next;
     }
}

static void append_maneuver_window_list (Maneuver_Window_Type *lst, Maneuver_Window_Type *x)
{
   Maneuver_Window_Type *w;

   for (w = lst; w != NULL; w = w->next)
     {
        if (w->next == NULL)
          {
             w->next = x;
             break;
          }
     }
}

static Maneuver_Window_Type *new_maneuver_window1 (double beg_timet, double end_timet)
{
   Maneuver_Window_Type *mw = NULL;

   if (NULL == (mw = (Maneuver_Window_Type *) MALLOC (sizeof *mw)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   mw->beg_timet = beg_timet;
   mw->end_timet = end_timet;
   mw->next = NULL;

   return mw;
}

typedef struct
{
   Maneuver_Window_Type *lst;
   double table_beg_timet;
   double table_end_timet;
}
Maneuver_Table_Type;
#define MANEUVER_TABLE_DEFAULT_INIT {NULL, -1.0, -1.0}

static void free_maneuver_table (Maneuver_Table_Type *mt)
{
   if (mt == NULL)
     return;
   free_maneuver_window_list (mt->lst);
}

/* The newest maneuver plan always supercedes prior maneuver plans,
 * therefore we process the maneuver files from newest to oldest,
 * halting either when we run out of files, or when we've
 * covered the target observing plan time interval.
 * It follows that each additional (prior) maneuver file may only
 * prepend or append maneuver windows.
 */
static int merge_prior_maneuver_file (Maneuver_Table_Type *mt, const char *maneuver_file,
                                      double plan_beg_timet, double plan_end_timet)
{
   IOCLib_String_Table_Type *st = NULL;
   IOCLib_KV_Table_Type *kv = NULL;
   Maneuver_Window_Type *head = NULL;
   Maneuver_Window_Type *win = NULL;
   Maneuver_Window_Type *prefix = NULL;
   Maneuver_Window_Type *suffix = NULL;
   double table_beg_timet, table_end_timet, beg_timet, end_timet;
   double mt_table_beg_timet, mt_table_end_timet;
   double last_end = -1.0;
   const char *column_names[] = {"window_start_time", "window_end_time"};
   unsigned int i, num_columns = sizeof(column_names) / sizeof (*column_names);
   unsigned int beg_col = 0, end_col = 1;
   int status = -1;

   if (NULL == (st = ioclib_csv_read_string_table (maneuver_file, column_names, num_columns)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading maneuver times from: %s", __func__, maneuver_file);
        return -1;
     }

   if (NULL == (kv = ioclib_extract_csv_metadata (st)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: extracting metadata from: %s", __func__, maneuver_file);
        goto return_status;
     }

   if ((0 != ioclib_kv_table_get_double (kv, "table_begin_time", &table_beg_timet))
       || (0 != ioclib_kv_table_get_double (kv, "table_end_time", &table_end_timet)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: extracting metadata from: %s", __func__, maneuver_file);
        goto return_status;
     }

   /* Does this maneuver table overlap the target plan interval?
    * If not, do nothing.
    */
   if ((table_end_timet < plan_beg_timet)
       || (plan_end_timet < table_beg_timet))
     {
        status = 0;
        goto return_status;
     }

   if (Plan_Verbose) fprintf (stderr, "using maneuver table: %s\n", maneuver_file);

   /* Validate maneuver table */
   for (i = 0; i < st->num_rows; i++)
     {
        if (*st->data[0][i] == ':')
          continue;
        if ((0 != ioclib_string_to_double (st->data[beg_col][i], &beg_timet))
            || (0 != ioclib_string_to_double (st->data[end_col][i], &end_timet)))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: error parsing maneuver window start/end times: row %d: %s",
                          __func__, i, maneuver_file);
             goto return_status;
          }
        if ((0 == isfinite(beg_timet))
            || (0 == isfinite(end_timet))
            || (beg_timet < table_beg_timet)
            || (table_end_timet < end_timet)
            || (end_timet < beg_timet))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: invalid maneuver table: %s", __func__, maneuver_file);
             goto return_status;
          }
        if (last_end < 0)
          {
             last_end = end_timet;
          }
        else if (beg_timet < last_end)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: invalid maneuver table: %s", __func__, maneuver_file);
             goto return_status;
          }
     }

   /* This time range has been finalized.  The new table may only append or prepend to it. */
   mt_table_beg_timet = (mt->table_beg_timet < 0) ? 0.0 : mt->table_beg_timet;
   mt_table_end_timet = (mt->table_end_timet < 0) ? 0.0 : mt->table_end_timet;

   head = mt->lst;

   /* Process the new table */
   for (i = 0; i < st->num_rows; i++)
     {
        if (*st->data[0][i] == ':')
          continue;

        if ((0 != ioclib_string_to_double (st->data[beg_col][i], &beg_timet))
            || (0 != ioclib_string_to_double (st->data[end_col][i], &end_timet)))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: error parsing maneuver window start/end times: row %d: %s",
                          __func__, i, maneuver_file);
             goto return_status;
          }

        /* Skip maneuver windows outside the target plan interval */
        if (end_timet < plan_beg_timet)
          continue;
        if (plan_end_timet < beg_timet)
          break;

        /* When the very first table is processed, we'll have
         * head=NULL, and 0=mt_table_beg_timet=mt_table_end_timet,
         * so that the first window will define 'head' and
         * the remaining windows will define a suffix list.
         */
        if (head == NULL)
          {
             if (NULL == (head = new_maneuver_window1 (beg_timet, end_timet)))
               goto return_status;
             continue;
          }

        /* Accumulate prefix and suffix lists, skipping maneuver windows
         * that overlap the already finalized time interval.
         */
        if (end_timet < mt_table_beg_timet)
          {
             if (NULL == (win = new_maneuver_window1 (beg_timet, end_timet)))
               goto return_status;
             if (prefix)
               {
                  append_maneuver_window_list (prefix, win);
               }
             else prefix = win;
          }
        else if (beg_timet > mt_table_end_timet)
          {
             if (NULL == (win = new_maneuver_window1 (beg_timet, end_timet)))
               goto return_status;
             if (suffix)
               {
                  append_maneuver_window_list (suffix, win);
               }
             else suffix = win;
          }
     }

   /* Connect prefix and suffix lists */
   append_maneuver_window_list (head, suffix);
   if (prefix)
     {
        append_maneuver_window_list (prefix, head);
        mt->lst = prefix;
     }
   else mt->lst = head;

   /* Update finalized time range */
   if ((mt->table_beg_timet < 0) || (table_beg_timet < mt->table_beg_timet))
     {
        mt->table_beg_timet = table_beg_timet;
     }
   if ((mt->table_end_timet < 0) || (mt->table_end_timet < table_end_timet))
     {
        mt->table_end_timet = table_end_timet;
     }

   status = 0;
return_status:
   ioclib_free_string_table (st);
   ioclib_kv_table_free (kv);
   if (status)
     {
        free_maneuver_window_list (prefix);
        free_maneuver_window_list (suffix);
     }

   return status;
}

typedef struct Line_Type Line_Type;
struct Line_Type
{
   Line_Type *next;
   char *text;
};

#define MT_BUFSIZE 1024

static void free_line1 (Line_Type *lt)
{
   if (lt == NULL)
     return;
   FREE(lt->text);
   FREE(lt);
}

static void free_line_list (Line_Type *lt)
{
   if (lt == NULL)
     return;

   while (lt)
     {
        Line_Type *next = lt->next;
        free_line1 (lt);
        lt = next;
     }
}

static Line_Type *new_line (const char *s)
{
   Line_Type *lt = NULL;
   char buf[MT_BUFSIZE];

   /* use sscanf to filter whitespace */
   if (1 != sscanf (s, "%s", buf))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: error parsing text: %s", __func__, s);
        return NULL;
     }

   if (NULL == (lt = (Line_Type *)MALLOC (sizeof *lt)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   lt->next = NULL;

   if (NULL == (lt->text = strdup (buf)))
     {
        FREE(lt);
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   return lt;
}

static Line_Type *read_lines_and_reverse_order (const char *file)
{
   FILE *fp = NULL;
   Line_Type *rlines = NULL;
   Line_Type *lt;
   char buf[MT_BUFSIZE];

   if (NULL == (fp = fopen (file, "r")))
     {
        tell_verror (TELL_IO_OPEN_ERROR, "%s: opening file for reading: %s", __func__, file);
        return NULL;
     }

   for (;;)
     {
        if (NULL == fgets (buf, sizeof(buf), fp))
          break;

        if (NULL == (lt = new_line (buf)))
          {
             free_line_list (rlines);
             rlines = NULL;
             break;
          }

        if (rlines == NULL)
          {
             rlines = lt;
          }
        else
          {
             /* final line list will contain the lines in reverse order */
             lt->next = rlines;
             rlines = lt;
          }
     }

   fclose (fp);
   return rlines;
}

static int process_maneuver_file_dir (Maneuver_Table_Type *mt, const char *maneuver_file_dir,
                                      double plan_beg_timet, double plan_end_timet)
{
   IOCLib_Listdir_Type *lst = NULL;
   int status = -1;
   size_t i;

   if (NULL == (lst = ioclib_listdir (maneuver_file_dir, IOCLIB_LISTDIR_SORT)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: listdir failed: %s", __func__, maneuver_file_dir);
        return -1;
     }

   /* Files are in ascending order. The maneuver file naming scheme
    * ensures that this corresponds to time order, so the last file
    * is the newest.  Process the files from newest to oldest:
    */
   i = lst->num_files;
   while (i-- > 0)
     {
        char buf[MT_BUFSIZE];
        int n;
        n = snprintf (buf, MT_BUFSIZE, "%s/%s", maneuver_file_dir, lst->files[i]);
        if ((n < 0) || (n >= MT_BUFSIZE))
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: file path exceeds buffer size", __func__);
             goto return_status;
          }
        if (0 != merge_prior_maneuver_file (mt, buf, plan_beg_timet, plan_end_timet))
          goto return_status;
     }

   status = 0;
return_status:
   ioclib_listdir_free (lst);
   return status;
}

static int process_maneuver_file (Maneuver_Table_Type *mt, const char *maneuver_file,
                                 double plan_beg_timet, double plan_end_timet)
{
   /* '@' indicates the maneuver_file contains a list of maneuver files,
    * assumed to be sorted in time order, oldest first.
    */
   if (*maneuver_file == '@')
     {
        Line_Type *rlines = NULL;
        Line_Type *lt;
        maneuver_file++;
        if (NULL == (rlines = read_lines_and_reverse_order (maneuver_file)))
          return -1;
        for (lt = rlines; lt != NULL; lt = lt->next)
          {
             if (0 != merge_prior_maneuver_file (mt, lt->text, plan_beg_timet, plan_end_timet))
               {
                  free_line_list (rlines);
                  return -1;
               }
          }
        free_line_list (rlines);
        return 0;
     }

   if (1 == ioclib_isdir (maneuver_file, "r"))
     {
        return process_maneuver_file_dir (mt, maneuver_file, plan_beg_timet, plan_end_timet);
     }

   return merge_prior_maneuver_file (mt, maneuver_file, plan_beg_timet, plan_end_timet);
}

static int include_maneuvers (Plan_List_Type *plan_list, const char *maneuver_file,
                              const Cal_Date_Type *t0, int num_plan_days)
{
   Maneuver_Table_Type mt = MANEUVER_TABLE_DEFAULT_INIT;
   Maneuver_Window_Type *win;
   double jd_utc0, jd_utc1, plan_beg_timet, plan_end_timet;
   double unix_epoch_jd = get_unix_epoch_jd();
   int status = -1;

   if (maneuver_file == NULL)
     return 0;

   /* Set target plan time interval */

   jd_utc0 = novas_julian_date (t0->year, t0->month, t0->day, t0->hour);
   jd_utc1 = jd_utc0 + num_plan_days;

   plan_beg_timet = (jd_utc0 - unix_epoch_jd) * SEC_PER_DAY;
   plan_end_timet = (jd_utc1 - unix_epoch_jd) * SEC_PER_DAY;

   if (0 != process_maneuver_file (&mt, maneuver_file, plan_beg_timet, plan_end_timet))
     return -1;

   for (win = mt.lst; win != NULL; win = win->next)
     {
        double mnv_beg, mnv_end;

        if (0 != utc_to_jd_utc (win->beg_timet, &mnv_beg))
          goto return_status;
        if (0 != utc_to_jd_utc (win->end_timet, &mnv_end))
          goto return_status;

        if (0 != insert_maneuver_gap (plan_list, mnv_beg, mnv_end))
          goto return_status;
     }

   status = 0;
return_status:
   free_maneuver_table (&mt);
   return status;
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

static int perform_sza_check (Solar_Geom_Type *solar_geom,
                              const char *sza_check_string)
{
   double jd_utc, lon, lat, sza;
   char *utc_str = NULL;
   char *p = NULL;

   if (2 != sscanf (sza_check_string, "%le,%le,", &lon, &lat))
     {
        fprintf (stderr, "*** Error: parsing string: %s\n",
                 sza_check_string ? sza_check_string : "(null)");
        return -1;
     }

   if (NULL == (p = strrchr (sza_check_string, ',')))
     {
        fprintf (stderr, "*** Error: parsing string: %s\n",
                 sza_check_string ? sza_check_string : "(null)");
        return -1;
     }

   utc_str = p + 1;

   if (0 != utcstr_to_jd_utc (utc_str, &jd_utc))
     return -1;

   if (0 != solar_geom->sgt_solar_zenith_angle (solar_geom, jd_utc, lon, lat, &sza))
     return -1;

   (void) fprintf (stdout, "%f\n", sza);

   return 0;
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
   int enable_twilight_scan = 0;
   int have_date = 0;
   Cal_Date_Type t0 = {0};
   int ndays_since_epoch = 0;
   const char *maneuver_file = NULL;
   const char *scan_tailoring_file = NULL;
   const char *sza_check_string = NULL;
   Optional_Output_Type oot =
     {
        .vis_output_file = NULL,
        .num_sza_days = DEFAULT_NUM_SZA_DAYS,
     };
   static struct option long_options[] =
     {
        {"help",         no_argument,       0, 'h'},
        {"date",         required_argument, 0, 'd'},
        {"config",       required_argument, 0, 'c'},
        {"ndays",        required_argument, 0, 'n'},
        {"nightlights",  no_argument,       0, 'N'},
        {"scan",         required_argument, 0, 's'},
        {"type",         required_argument, 0, 't'},
        {"output",       required_argument, 0, 'o'},
        {"irr",          required_argument, 0, 'I'},
        {"szaout",       required_argument, 0, 'z'},
        {"tailor",       required_argument, 0, 'T'},
        {"maneuver",     required_argument, 0, 'M'},
        {"master",       no_argument,       0, 'm'},
        {"verbose",      no_argument,       0, 'v'},
        {"Zenith",       required_argument, 0, 'Z'},
        {0,0,0,0}
     };
   config_t cfg = {0};
   Ephem_Type eph = {0};
   Scan_Type *scan = NULL;
   Twilight_Scan_Type *twilight_scan = NULL;
   Split_Scan_Type *split_scan = NULL;
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
        int c = getopt_long (argc, argv, "hNvZ:M:c:d:I:m:n:o:s:t:T:z:", long_options, &option_index);
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
           case 'Z':
             sza_check_string = optarg;
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
           case 'M':
             maneuver_file = optarg;
             break;
           case 'N':
             enable_twilight_scan++;
             break;
           case 'v':
             Plan_Verbose++;
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
             scan_tailoring_file = optarg;
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

   set_unix_epoch_jd ();

   if (sza_check_string)
     {
        if ((0 != read_epoch (&cfg))
            || (0 != ephem_open (&cfg, &eph))
            || (NULL == (solar_geom = solar_geom_init (&cfg))))
          goto return_status;
        if (0 == perform_sza_check (solar_geom, sza_check_string))
          status = EXIT_SUCCESS;
        goto return_status;
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

   if (scan_tailoring_file)
     {
        if (0 == print_scan_tailoring_file (solar_geom, &t0, num_plan_days, scan_tailoring_file))
          status = EXIT_SUCCESS;
        goto return_status;
     }

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

   if (enable_twilight_scan)
     {
        if (NULL == (twilight_scan = twilight_scan_open (&cfg)))
          goto return_status;
     }

   if (0 == strncmp (scan_method, "split", 5))
     {
        if (NULL == (split_scan = split_scan_open (&cfg, scan_method)))
          goto return_status;
     }

   plan_list = generate_scan_plan (&eph, solar_geom, scan, sm, &t0, num_plan_days,
                                   twilight_scan, split_scan);
   if (NULL == plan_list)
     goto return_status;

   if (0 != include_maneuvers (plan_list, maneuver_file, &t0, num_plan_days))
     goto return_status;

   if (0 != verify_safety_constraints (solar_geom, scan, plan_list))
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
   if (twilight_scan) twilight_scan->tst_delete (twilight_scan);
   if (split_scan) split_scan->sst_delete (split_scan);
   (void) ephem_close (&eph);
   plan_list_free (plan_list);
   close_outfile (fp_scan, scan_outfile);
   close_outfile (fp_master, master_outfile);
   close_outfile (fp_irr, irr_outfile);
   config_destroy (&cfg);
   tell_close ();

   return status;
}
