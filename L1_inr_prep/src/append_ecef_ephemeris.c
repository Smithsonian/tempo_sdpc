#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>

#include <libnovas.h>
#include <J2K2ECEF.h>
#include <readIERSBulletinA.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#define TT_OFFSET_SECS  32.184
#define SEC_PER_DAY     86400

#define TIMESTAMP_UNIT_STRING_SIZE 128
#define MAX_ISOTIME_LEN 32
/*      MAX_ISOTIME_LEN must hold:  yyyy-mm-ddThh:mm:ss.sssZ */

typedef struct
{
   double *t;
   double *x;
   double *y;
   double *z;
   int *quality;
   size_t n;
}
Vector_Array_Type;

typedef struct J2K_Type J2K_Type;
struct J2K_Type
{
   int (*taix_to_tdb)(const J2K_Type *j2k, double taix, double *tdb);
   int (*j2k_to_ecef)(const J2K_Type *j2k, double tdb, double *vec_j2k, double *vec_ecef);
   double last_jd_utc;
   double ut1_minus_utc;
   double lon;
   int leap_secs;
};

static int taix_to_jd (double taix, double *jd_utc)
{
   struct cal_date_type {
      short int year;
      short int month;
      short int day;
      double hour;
   } t0;
   int year, month, day;

   if (0 != tio_time_taix_to_utc_caldate (taix, &year, &month, &day, &t0.hour))
     return -1;
   t0.year = year;
   t0.month = month;
   t0.day = day;

   *jd_utc = novas_julian_date (t0.year, t0.month, t0.day, t0.hour);

#ifdef DEBUG
   fprintf (stdout, "jd_utc = %10.4f   y/m/d/h = %d-%02d-%02d %7.4f\n",
            *jd_utc, t0.year, t0.month, t0.day, t0.hour);
#endif

   return 0;
}

static int taix_to_tdb (const J2K_Type *j2k, double taix, double *tdb)
{
   double jd_utc;

   if (0 != taix_to_jd (taix, &jd_utc))
     return -1;

   *tdb = jd_utc + (j2k->leap_secs + TT_OFFSET_SECS)/SEC_PER_DAY;

   return 0;
}

static int j2k_to_ecef (const J2K_Type *j2k, double tdb, double *vec_j2k, double *vec_ecef)
{
   TempoECEFErr err;

   err = j2KtoECEF (tdb, j2k->lon, j2k->leap_secs, j2k->ut1_minus_utc, vec_j2k, vec_ecef);
   if (err != TEMPO_ECEF_NO_ERR)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: j2ktoECEF failed (err=%d)", __func__, err);
        return -1;
     }

   return 0;
}

static int j2k_update (J2K_Type *j2k, double taix, char *iers_bulletin)
{
   TempoBulletinErr err;
   double jd_utc;

   if (j2k == NULL)
     return -1;

   if (0 != taix_to_jd (taix, &jd_utc))
     return -1;

   /* Don't re-read the IERS bulletin unless the date has changed */
   if ((j2k->last_jd_utc > 0.0)
       && (fabs (jd_utc - j2k->last_jd_utc) < 1.0))
     return 0;

   j2k->last_jd_utc = jd_utc;
   j2k->taix_to_tdb = taix_to_tdb;
   j2k->j2k_to_ecef = j2k_to_ecef;

   /* Greenwich longitude for standard ECEF */
   j2k->lon = 0.0;

   if ((err = readIERSBulletinA (iers_bulletin, jd_utc, &j2k->leap_secs, &j2k->ut1_minus_utc)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: readIERSBulletinA failed (err=%d): %s",
                     __func__, err, iers_bulletin);
        return -1;
     }

   return 0;
}

static int map_j2k_to_ecef (const J2K_Type *j2k, int n, double *taix, double *X, double *Y, double *Z)
{
   double tshift_days, tdb, vec_ecef[3], vec_j2k[3];
   int i;

   if (0 != j2k->taix_to_tdb (j2k, taix[0], &tdb))
     return -1;

   /* Assuming no leap seconds are inserted during taix[0:n-1],
    * the conversion from TAI to TDB is just a linear shift, so
    * we can shift all sample times by the same amount (in days).
    */
   tshift_days = tdb - taix[0] / SEC_PER_DAY;

   for (i = 0; i < n; i++)
     {
        vec_j2k[0] = X[i];
        vec_j2k[1] = Y[i];
        vec_j2k[2] = Z[i];
        tdb = taix[i] / SEC_PER_DAY + tshift_days;
        if (0 != j2k->j2k_to_ecef (j2k, tdb, vec_j2k, vec_ecef))
          return -1;
        X[i] = vec_ecef[0];
        Y[i] = vec_ecef[1];
        Z[i] = vec_ecef[2];
     }

   return 0;
}

static Vector_Array_Type *alloc_vecarray (size_t n)
{
   Vector_Array_Type *v = NULL;

   if ((NULL == (v = (Vector_Array_Type *)malloc (sizeof *v)))
       || (NULL == (v->t = (double *)malloc (4 * n * sizeof(double))))
       || (NULL == (v->quality = (int *)malloc (n * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        free(v);
        return NULL;
     }
   v->x = v->t + n;
   v->y = v->x + n;
   v->z = v->y + n;
   v->n = n;

   return v;
}

static void free_vecarray (Vector_Array_Type *v)
{
   if (v == NULL)
     return;
   free(v->t);
   free(v->quality);
   free(v);
}

static int convert_j2k_to_ecef (const J2K_Type *j2k, Vector_Array_Type *v)
{
   return map_j2k_to_ecef (j2k, v->n, v->t, v->x, v->y, v->z);
}

static int read_vec (int grp, Vector_Array_Type *v,
                     const char *xname, const char *yname, const char *zname)
{
   int start = 0;
   int count = v->n;

   if ((0 != TIO_get_var_section (grp, xname, &start, &count, NC_DOUBLE, v->x))
       ||(0 != TIO_get_var_section (grp, yname, &start, &count, NC_DOUBLE, v->y))
       ||(0 != TIO_get_var_section (grp, zname, &start, &count, NC_DOUBLE, v->z)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: reading ephemeris vector", __func__);
        return -1;
     }

   return 0;
}

static int write_vec (int grp, int start, Vector_Array_Type *v,
                      const char *xname, const char *yname, const char *zname)
{
   int count = v->n;

   if ((0 != TIO_put_var_section (grp, xname, &start, &count, NC_DOUBLE, v->x))
       ||(0 != TIO_put_var_section (grp, yname, &start, &count, NC_DOUBLE, v->y))
       ||(0 != TIO_put_var_section (grp, zname, &start, &count, NC_DOUBLE, v->z)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: writing ephemeris vector", __func__);
        return -1;
     }

   return 0;
}

/* returns 1: success, 0: no data, -1: error */
static int read_ephem (int grp, const char *time_var, const char *quality_var,
                       const char **pos_vars, const char **vel_vars,
                       Vector_Array_Type **posp, Vector_Array_Type **velp)
{
   Vector_Array_Type *pos = NULL;
   Vector_Array_Type *vel = NULL;
   size_t time_dimlen;
   int time_dimid, varid, start, count;
   int status = -1;

   if (0 != tio_inq_varid (grp, time_var, &varid))
     {
        tell_vlog (TELL_MSGTYPE_INFO, 2, "%s variable is not present", time_var);
        return 0;
     }

   if (0 != TIO_inq_dim (grp, "time", &time_dimid, &time_dimlen))
     return -1;

   if (time_dimlen == 0)
     {
        tell_vlog (TELL_MSGTYPE_INFO, 2, "%s variable has 0 samples", time_var);
        return 0;
     }

   if ((NULL == (pos = alloc_vecarray (time_dimlen)))
       || (NULL == (vel = alloc_vecarray (time_dimlen))))
     goto return_status;

   start = 0;
   count = time_dimlen;

   if (0 != TIO_get_var_section (grp, time_var, &start, &count, NC_DOUBLE, pos->t))
     goto return_status;

   if (quality_var)
     {
        if (0 != TIO_get_var_section (grp, quality_var, &start, &count, NC_INT, pos->quality))
          goto return_status;
        memcpy ((char *)vel->quality, (char *)pos->quality, time_dimlen * sizeof(int));
     }

   if (0 != read_vec (grp, pos, pos_vars[0], pos_vars[1], pos_vars[2]))
     goto return_status;
   if (0 != read_vec (grp, vel, vel_vars[0], vel_vars[1], vel_vars[2]))
     goto return_status;

   memcpy ((char *)vel->t, (char *)pos->t, time_dimlen * sizeof(*pos->t));

   *posp = pos;
   *velp = vel;

   status = 1;
return_status:
   if (status < 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading ephemeris", __func__);
        free_vecarray (pos);
        free_vecarray (vel);
     }
   return status;
}

static int process_dop_ephem (const J2K_Type *j2k, Vector_Array_Type *pos, Vector_Array_Type *vel)
{
   if ((0 != convert_j2k_to_ecef (j2k, pos))
       || (0 != convert_j2k_to_ecef (j2k, vel)))
     return -1;

   return 0;
}

static void __add_orbital_velocity_vector (const double satx, const double saty,
                                           double *satvx, double *satvy)
{
   double v_geo = 3.074666284127684;  /* km/sec */
   double phi = atan2 (saty, satx);
   /* Add geostationary orbital velocity component:
    * x = R cos(phi)  -> dx/dt = - R sin(phi) dphi/dt = - v_geo * sin(phi)
    * y = R sin(phi)  -> dy/dt =   R cos(phi) dphi/dt =   v_geo * cos(phi)
    */
   *satvx += - v_geo * sin(phi);
   *satvy +=   v_geo * cos(phi);
}

static int process_gpsr_ephem (const Vector_Array_Type *pos, Vector_Array_Type *vel)
{
   double convert_mm_per_sec_to_km_per_sec = 1.e-6;
   size_t i;

   for (i = 0; i < pos->n; i++)
     {
        vel->x[i] *= convert_mm_per_sec_to_km_per_sec;
        vel->y[i] *= convert_mm_per_sec_to_km_per_sec;
        vel->z[i] *= convert_mm_per_sec_to_km_per_sec;
        __add_orbital_velocity_vector (pos->x[i], pos->y[i], &vel->x[i], &vel->y[i]);
     }

   return 0;
}

static int process_file (int dest_ncid, J2K_Type *j2k, char *iers_bulletin, const char *file)
{
   typedef struct
     {
        Vector_Array_Type *pos;
        Vector_Array_Type *vel;
     }
   ephem_type;
   ephem_type dop = {0}, gpsr = {0};
   typedef struct
     {
        const char *time_var;
        const char *pos_vars[3];
        const char *vel_vars[3];
     }
   ephem_var_type;
   ephem_var_type dop_n =
     {
        "anc_gps_time",
        {"anc_satx", "anc_saty", "anc_satz"},
        {"anc_satvx", "anc_satvy", "anc_satvz"}
     };
   ephem_var_type gpsr_n =
     {
        "anc_gpsr_gps_time",
        {"anc_gpsr_satx", "anc_gpsr_saty", "anc_gpsr_satz"},
        {"anc_gpsr_satvx", "anc_gpsr_satvy", "anc_gpsr_satvz"}
     };
   struct
     {
        int grp, start, count;
     }
   dest;
   int ncid, grp, dimid;
   int have_dop, have_gpsr;
   size_t dimlen;
   int status = -1;

   if (0 != TIO_inq_grp (dest_ncid, "ephemeris", &dest.grp))
     goto return_status;

   if (0 != TIO_inq_dim (dest.grp, "time", &dimid, &dimlen))
     return -1;
   dest.start = dimlen;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   if (0 != TIO_inq_grp (ncid, "anc_gps", &grp))
     goto return_status;

   if ((have_dop = read_ephem (grp, dop_n.time_var, NULL, dop_n.pos_vars, dop_n.vel_vars,
                               &dop.pos, &dop.vel)) < 0)
     {
        goto return_status;
     }

   if (have_dop)
     {
        if (0 != j2k_update (j2k, dop.pos->t[0], iers_bulletin))
          return -1;
        if (0 != process_dop_ephem (j2k, dop.pos, dop.vel))
          goto return_status;
        dest.count = dop.pos->n;
        if (0 != TIO_put_var_section (dest.grp, "dop_time", &dest.start, &dest.count, NC_DOUBLE, dop.pos->t))
          goto return_status;
        if (0 != write_vec (dest.grp, dest.start, dop.pos, "dop_x", "dop_y", "dop_z"))
          goto return_status;
        if (0 != write_vec (dest.grp, dest.start, dop.vel, "dop_vx", "dop_vy", "dop_vz"))
          goto return_status;
     }

   if ((have_gpsr = read_ephem (grp, gpsr_n.time_var, "anc_gpsr_navsoln_mth_raw",
                                gpsr_n.pos_vars, gpsr_n.vel_vars, &gpsr.pos, &gpsr.vel)) < 0)
     {
        goto return_status;
     }

   if (have_gpsr)
     {
        if (0 != process_gpsr_ephem (gpsr.pos, gpsr.vel))
          goto return_status;
        dest.count = gpsr.pos->n;
        if (0 != TIO_put_var_section (dest.grp, "gpsr_time", &dest.start, &dest.count, NC_DOUBLE, gpsr.pos->t))
          goto return_status;
        if (0 != write_vec (dest.grp, dest.start, gpsr.pos, "gpsr_x", "gpsr_y", "gpsr_z"))
          goto return_status;
        if (0 != write_vec (dest.grp, dest.start, gpsr.vel, "gpsr_vx", "gpsr_vy", "gpsr_vz"))
          goto return_status;
        if (0 != TIO_put_var_section (dest.grp, "gpsr_navsoln_mth", &dest.start, &dest.count, NC_INT, gpsr.pos->quality))
          goto return_status;
     }

   status = 0;
return_status:
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: processing: %s", __func__, file ? file : "<null>");
     }
   (void) TIO_close (ncid);
   free_vecarray (dop.pos);
   free_vecarray (dop.vel);
   free_vecarray (gpsr.pos);
   free_vecarray (gpsr.vel);
   return status;
}

typedef struct
{
   const char *names[3];
   const char *units[3];
   const char *descr[3];
}
Vector_Info_Type;

static int enable_compression (int grp, int varid)
{
#ifndef ENABLE_COMPRESSION
   (void) grp; (void) varid;
#else
   int shuffle=1, deflate=1, deflate_level=1;
   size_t chunksizes = 256;
   if (-1 == TIO_def_var_deflate (grp, varid, shuffle, deflate, deflate_level))
     return -1;
   if (-1 == TIO_def_var_chunking (grp, varid, NC_CHUNKED, &chunksizes))
     return -1;
#endif
   return 0;
}

static int create_int_time_series (int grp, int dimid, const char *name, const char *descr)
{
   int len, varid;

   if (0 != TIO_def_var (grp, name, NC_INT, 1, &dimid, &varid))
     return -1;
   len = strlen(descr) + 1;
   if (0 != TIO_put_att (grp, varid, "comment", NC_CHAR, len, descr))
     return -1;

   return enable_compression (grp, varid);
}

static int create_dbl_time_series (int grp, int dimid, const char *name, const char *units, const char *descr)
{
   int len, varid;

   if (0 != TIO_def_var (grp, name, NC_DOUBLE, 1, &dimid, &varid))
     return -1;
   len = strlen(descr) + 1;
   if (0 != TIO_put_att (grp, varid, "comment", NC_CHAR, len, descr))
     return -1;
   len = strlen(units) + 1;
   if (0 != TIO_put_att (grp, varid, "units", NC_CHAR, len, units))
     return -1;

   return enable_compression (grp, varid);
}

static int create_vector_vs_time (int grp, int dimid, const Vector_Info_Type *vi)
{
   int i;

   for (i = 0; i < 3; i++)
     {
        if (0 != create_dbl_time_series (grp, dimid, vi->names[i], vi->units[i], vi->descr[i]))
          return -1;
     }

   return 0;
}

static int create_output_file (const char *path)
{
   Vector_Info_Type dop_pos =
     {
        {"dop_x", "dop_y", "dop_z"}, {"km", "km", "km"},
        {"satellite X coordinate derived from DOP ephemeris",
         "satellite Y coordinate derived from DOP ephemeris",
         "satellite Z coordinate derived from DOP ephemeris"}
     };
   Vector_Info_Type dop_vel =
     {
        {"dop_vx", "dop_vy", "dop_vz"}, {"km/s", "km/s", "km/s"},
        {"satellite X velocity component derived from DOP ephemeris",
         "satellite Y velocity component derived from DOP ephemeris",
         "satellite Z velocity component derived from DOP ephemeris"}
     };
   Vector_Info_Type gpsr_pos =
     {
        {"gpsr_x", "gpsr_y", "gpsr_z"}, {"km", "km", "km"},
        {"satellite X coordinate derived from GPSR ephemeris",
         "satellite Y coordinate derived from GPSR ephemeris",
         "satellite Z coordinate derived from GPSR ephemeris"}
     };
   Vector_Info_Type gpsr_vel =
     {
        {"gpsr_vx", "gpsr_vy", "gpsr_vz"}, {"km/s", "km/s", "km/s"},
        {"satellite X velocity component derived from GPSR ephemeris",
         "satellite Y velocity component derived from GPSR ephemeris",
         "satellite Z velocity component derived from GPSR ephemeris"}
     };
   char unit_string[TIMESTAMP_UNIT_STRING_SIZE];
   char epoch[MAX_ISOTIME_LEN];
   int ncid, grp, dimid_time;
   int n, len;

   if (0 != TIO_mktimestamp_str (0.0, 1, epoch, sizeof(epoch)))
     return -1;

   len = sizeof(unit_string);

   memset ((char *)unit_string, 0, len);
   if (((n = snprintf (unit_string, len, "seconds since %s", epoch)) < 0)
       || (n >= (int) len))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: snprintf failed: buffer size=%d  return value=%d",
                     __func__, len, n);
        return -1;
     }

   if (0 != TIO_create (path, NC_NETCDF4, &ncid))
     return -1;

   if (0 != tio_write_epoch_timestamp (ncid, NC_GLOBAL))
     return -1;

   if (0 != TIO_def_grp (ncid, "ephemeris", &grp))
     return -1;

   if (0 != TIO_def_dim (grp, "time", NC_UNLIMITED, &dimid_time))
     return -1;

   if (0 != create_dbl_time_series (grp, dimid_time, "dop_time", unit_string, "DOP ephemeris timestamp"))
     return -1;
   if ((0 != create_vector_vs_time (grp, dimid_time, &dop_pos))
       || (0 != create_vector_vs_time (grp, dimid_time, &dop_vel)))
     return -1;

   if (0 != create_dbl_time_series (grp, dimid_time, "gpsr_time", unit_string, "GPSR ephemeris timestamp"))
     return -1;
   if ((0 != create_vector_vs_time (grp, dimid_time, &gpsr_pos))
       || (0 != create_vector_vs_time (grp, dimid_time, &gpsr_vel)))
     return -1;
   if (0 != create_int_time_series (grp, dimid_time, "gpsr_navsoln_mth", "GPSR solution method"))
     return -1;

   return ncid;
}

static int use_data_epoch (const char *path)
{
   int ncid, status;
   if (0 != TIO_open (path, NC_NOWRITE, &ncid))
     return -1;
   status = tio_use_file_epoch (ncid);
   TIO_close (ncid);
   return status;
}

static int open_outfile (const char *path)
{
   if (0 == access (path, F_OK | W_OK))
     {
        int ncid;
        if (0 != TIO_open (path, NC_WRITE, &ncid))
          return -1;
        return ncid;
     }
   return create_output_file (path);
}

static void usage (void)
{
   fprintf (stderr, "Usage: append_ecef_ephemeris [options] FILE\n");
   fprintf (stderr, "  Options:\n");
   fprintf (stderr, "   -h | --help          Print this message\n");
   fprintf (stderr, "   -i | --iers FILE     IERS bulletin A file\n");
   fprintf (stderr, "   -o | --output FILE   Destination file\n");
   fprintf (stderr, "   -v | --verbose       Verbose level [-vvv is more verbose]\n");
   exit (EXIT_SUCCESS);
}

int main (int argc, char **argv)
{
   const char appname[] = "append_ecef_ephemeris";
   int status = EXIT_FAILURE;
   static struct option long_options[] =
     {
        {"help",    no_argument,       0, 'h'},
        {"iers",    required_argument, 0, 'i'},
        {"output",  required_argument, 0, 'o'},
        {"verbose", no_argument,       0, 'v'},
        {0,0,0,0}
     };
   const char *dest_file = NULL;
   char *iers_bulletin = NULL;
   int log_level = 0;
   int dest_ncid;
   J2K_Type j2k = {0};

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hi:o:v", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "%s: getopt returned character %d??\n", __func__, c);
             goto return_status;
             break;
           case 'h':
             usage();
             break;
           case 'i':
             iers_bulletin = optarg;
             break;
           case 'o':
             dest_file = optarg;
             break;
           case 'v':
             log_level++;
             break;
          }
     }

   if ((dest_file == NULL) || (iers_bulletin == NULL)
       || (optind == 0))
     {
        usage();
     }

   (void) tell_set_log_level (TELL_MSGTYPE_INFO, log_level);

   if (0 != use_data_epoch (argv[optind]))
     goto return_status;

   if ((dest_ncid = open_outfile (dest_file)) < 0)
     goto return_status;

   for ( ; optind < argc; optind++)
     {
        const char *path = argv[optind];
        tell_vlog (TELL_MSGTYPE_INFO, 1, "processing: %s", path);
        if (process_file (dest_ncid, &j2k, iers_bulletin, path) < 0)
          goto return_status;
     }

   if (0 != TIO_close(dest_ncid))
     goto return_status;

   status = EXIT_SUCCESS;
return_status:
   return status;
}
