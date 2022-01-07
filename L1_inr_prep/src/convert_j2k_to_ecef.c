#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <libnovas.h>
#include <J2K2ECEF.h>
#include <readIERSBulletinA.h>

#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "convert_j2k_to_ecef.h"

#define TT_OFFSET_SECS  32.184
#define SEC_PER_DAY     86400

typedef struct
{
   double *t;
   double *x;
   double *y;
   double *z;
   size_t n;
}
Vector_Array_Type;

typedef struct J2K_Type J2K_Type;
struct J2K_Type
{
   int (*taix_to_tdb)(const J2K_Type *j2k, double taix, double *tdb);
   int (*j2k_to_ecef)(const J2K_Type *j2k, double tdb, double *vec_j2k, double *vec_ecef);
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
        tell_verror (TELL_RUNTIME_ERROR, "%s: j2ktoECEF failed: err=%d\n", __func__, err);
        return -1;
     }

   return 0;
}

static void j2k_close (J2K_Type *j2k)
{
   (void) j2k;
}

static int j2k_open (J2K_Type *j2k, double taix, char *iers_bulletin)
{
   TempoBulletinErr err;
   double jd_utc;

   if (j2k == NULL)
     return -1;

   j2k->taix_to_tdb = taix_to_tdb;
   j2k->j2k_to_ecef = j2k_to_ecef;

   /* Greenwich longitude for standard ECEF */
   j2k->lon = 0.0;

   if (0 != taix_to_jd (taix, &jd_utc))
     return -1;

   if ((err = readIERSBulletinA (iers_bulletin, jd_utc, &j2k->leap_secs, &j2k->ut1_minus_utc)) != 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: readIERSBulletinA failed (err=%d)\n", __func__, err);
        return -1;
     }
   tell_vinfo (0, "Using IERS bulletin %s", iers_bulletin);
   tell_vinfo (0, "TEMPO time=%f is JD=%f: leap_secs=%d UT1-UTC=%f",
               taix, jd_utc, j2k->leap_secs, j2k->ut1_minus_utc);

   return 0;
}

static Vector_Array_Type *alloc_vecarray (size_t n)
{
   Vector_Array_Type *v = NULL;

   if ((NULL == (v = (Vector_Array_Type *)malloc (sizeof *v)))
       || (NULL == (v->t = (double *)malloc (4 * n * sizeof(double)))))
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
   free(v);
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

static int convert_vec_j2k_to_ecef (int grp, const J2K_Type *j2k, Vector_Array_Type *v,
                                    const char *xname, const char *yname, const char *zname)
{
   int start = 0;
   int count = v->n;

   if ((0 != TIO_get_var_section (grp, xname, &start, &count, NC_DOUBLE, v->x))
       ||(0 != TIO_get_var_section (grp, yname, &start, &count, NC_DOUBLE, v->y))
       ||(0 != TIO_get_var_section (grp, zname, &start, &count, NC_DOUBLE, v->z)))
     return -1;

   if (0 != map_j2k_to_ecef (j2k, v->n, v->t, v->x, v->y, v->z))
     return -1;

   if ((0 != TIO_put_var_section (grp, xname, &start, &count, NC_DOUBLE, v->x))
       ||(0 != TIO_put_var_section (grp, yname, &start, &count, NC_DOUBLE, v->y))
       ||(0 != TIO_put_var_section (grp, zname, &start, &count, NC_DOUBLE, v->z)))
     return -1;

   return 0;
}

int convert_j2k_to_ecef (char *iers_bulletin, const char *file)
{
   Vector_Array_Type *v = NULL;
   J2K_Type j2k = {0};
   int ncid, grp, time_dimid, start, count;
   size_t time_dimlen;
   int status = -1;

   if (0 != TIO_open (file, NC_WRITE, &ncid))
     return -1;

   if (0 != tio_use_file_epoch (ncid))
     goto return_status;

   if (0 != TIO_inq_grp (ncid, "/inr_input/ephemeris", &grp))
     goto return_status;

   if (0 != TIO_inq_dim (grp, "ephemeris_time", &time_dimid, &time_dimlen))
     goto return_status;

   if (NULL == (v = alloc_vecarray (time_dimlen)))
     goto return_status;

   start = 0;
   count = time_dimlen;

   if (0 != TIO_get_var_section (grp, "ephemeris_time", &start, &count, NC_DOUBLE, v->t))
     goto return_status;

   if (0 != j2k_open (&j2k, v->t[0], iers_bulletin))
     goto return_status;

   if (0 != convert_vec_j2k_to_ecef (grp, &j2k, v, "satellite_X", "satellite_Y", "satellite_Z"))
     goto return_status;

   if (0 != convert_vec_j2k_to_ecef (grp, &j2k, v, "satellite_velocity_X", "satellite_velocity_Y", "satellite_velocity_Z"))
     goto return_status;

   status = 0;
return_status:
   if (status)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: failed", __func__);
     }
   j2k_close (&j2k);
   free_vecarray (v);
   (void) TIO_close (ncid);
   return status;
}

#ifdef __TEST__

int main (int argc, char **argv)
{
   char *iers_bulletin;

   if (argc < 3)
     {
        fprintf (stderr, "\nUsage: j2k_to_ecef <iers-bulletin-A> FILE [FILE ...]\n\n");
        fprintf (stderr, "  Transform Level 1 radiance file ephemeris vectors:\n");
        fprintf (stderr, "    /inr_input/ephemeris/satellite_[XYZ]\n");
        fprintf (stderr, "    /inr_input/ephemeris/satellite_velocity_[XYZ]\n");
        fprintf (stderr, "  from J2000 coordinate frame to ECEF coordinate frame.\n");
        fprintf (stderr, "\n  *** WARNING: Files are modified in place!\n\n");
        return 1;
     }

   argc--;
   argv++;

   iers_bulletin = *argv;
   argc--;
   argv++;

   do
     {
        fprintf (stdout, "processing %s\n", *argv);
        if (0 != convert_j2k_to_ecef (iers_bulletin, *argv))
          return 1;
        argv++;
        argc--;
     }
   while (argc > 0);

   fprintf (stdout, "done\n");

   return 0;
}
#endif
