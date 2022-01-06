#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <libnovas.h>
#include <J2K2ECEF.h>
#include <readIERSBulletinA.h>

#include <tio.h>
#include <tio_template.h>

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
   int (*ecef_to_j2k)(const J2K_Type *j2k, double tdb, const double *vec_ecef, double *vec_j2k);
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
        fprintf (stderr, "*** Error: j2ktoECEF failed: err=%d\n", err);
        return -1;
     }

   return 0;
}

static int ecef_to_j2k (const J2K_Type *j2k, double tdb, const double *vec_ecef, double *vec_j2k)
{
   double x[3] = {1.0, 0.0, 0.0};
   double y[3] = {0.0, 1.0, 0.0};
   double z[3] = {0.0, 0.0, 1.0};
   double Rt[9];
   double X, Y, Z;

   /*  Linear mapping J2K -> ECEF corresponds to a matrix multiplication:
    *    vec_ecef = R * vec_j2k
    * where the first column of R is the transformation of a
    * unit vector [1, 0, 0]; the second [0, 1, 0]; the third [0, 0, 1]:
    *
    *   R(:,0) = J2K2ECEF(t,[1, 0, 0]);
    *   R(:,1) = J2K2ECEF(t,[0, 1, 0]);
    *   R(:,2) = J2K2ECEF(t,[0, 0, 1]);
    *
    * R is an orthonormal rotation matrix, so the transpose, Rt,
    * gives the inverse transformation, ECEF -> J2K:
    *   vec_j2k = Rt * vec_ecef
    */

   if (   (0 != j2k->j2k_to_ecef (j2k, tdb, x, Rt+0))
       || (0 != j2k->j2k_to_ecef (j2k, tdb, y, Rt+3))
       || (0 != j2k->j2k_to_ecef (j2k, tdb, z, Rt+6)))
     {
        return -1;
     }

   X = vec_ecef[0];
   Y = vec_ecef[1];
   Z = vec_ecef[2];

   vec_j2k[0] = Rt[0] * X + Rt[1] * Y + Rt[2] * Z;
   vec_j2k[1] = Rt[3] * X + Rt[4] * Y + Rt[5] * Z;
   vec_j2k[2] = Rt[6] * X + Rt[7] * Y + Rt[8] * Z;

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
   j2k->ecef_to_j2k = ecef_to_j2k;

   /* Greenwich longitude for standard ECEF */
   j2k->lon = 0.0;

   if (0 != taix_to_jd (taix, &jd_utc))
     return -1;

   if ((err = readIERSBulletinA (iers_bulletin, jd_utc, &j2k->leap_secs, &j2k->ut1_minus_utc)) != 0)
     {
        fprintf (stderr, "*** Error: readIERSBulletinA failed (err=%d)\n", err);
        return -1;
     }

   return 0;
}

#ifdef DEBUG
static double diff_length (double *a, double *b, int n)
{
   double s = 0.0;
   int i;
   if (b)
     {
        for (i = 0; i < n; i++)
          {
             double d = a[i] - b[i];
             s += d*d;
          }
     }
   else
     {
        for (i = 0; i < n; i++)
          {
             s += a[i] * a[i];
          }
     }
   return sqrt(s);
}

static int check_conversion (const J2K_Type *j2k, double tdb, double *vec_ecef, double *vec_j2k)
{
   double delta, norm, vec[3];
   int i;

   if (0 != j2k->j2k_to_ecef (j2k, tdb, vec_j2k, vec))
     return -1;

   delta = diff_length (vec_ecef, vec, 3);
   norm = diff_length (vec_ecef, NULL, 3);
   if (delta < 1.e-6 * norm)
     return 0;

   fprintf (stderr, "*** Check failed:\nexpected got\n");
   for (i = 0; i < 3; i++)
     {
        fprintf (stderr, "%2d %15.7f %15.7f\n", i, vec_ecef[i], vec[i]);
     }

   return 1;
}
#endif

static int map_ecef_to_j2k (const J2K_Type *j2k, int n, double *taix, double *X, double *Y, double *Z)
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
        vec_ecef[0] = X[i];
        vec_ecef[1] = Y[i];
        vec_ecef[2] = Z[i];
        tdb = taix[i] / SEC_PER_DAY + tshift_days;
        if (0 != j2k->ecef_to_j2k (j2k, tdb, vec_ecef, vec_j2k))
          return -1;
        X[i] = vec_j2k[0];
        Y[i] = vec_j2k[1];
        Z[i] = vec_j2k[2];
#ifdef DEBUG
        if (0 != check_conversion (j2k, tdb, vec_ecef, vec_j2k))
          return -1;
#endif
     }

   return 0;
}

static Vector_Array_Type *alloc_vecarray (size_t n)
{
   Vector_Array_Type *v = NULL;

   if ((NULL == (v = (Vector_Array_Type *)malloc (sizeof *v)))
       || (NULL == (v->t = (double *)malloc (4 * n * sizeof(double)))))
     {
        fprintf (stderr, "***Error: malloc failed\n");
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

static int convert_ecef_to_j2k (int grp, const J2K_Type *j2k, Vector_Array_Type *v,
                                const char *xname, const char *yname, const char *zname)
{
   int start = 0;
   int count = v->n;

   if ((0 != TIO_get_var_section (grp, xname, &start, &count, NC_DOUBLE, v->x))
       ||(0 != TIO_get_var_section (grp, yname, &start, &count, NC_DOUBLE, v->y))
       ||(0 != TIO_get_var_section (grp, zname, &start, &count, NC_DOUBLE, v->z)))
     return -1;

   if (0 != map_ecef_to_j2k (j2k, v->n, v->t, v->x, v->y, v->z))
     return -1;

   if ((0 != TIO_put_var_section (grp, xname, &start, &count, NC_DOUBLE, v->x))
       ||(0 != TIO_put_var_section (grp, yname, &start, &count, NC_DOUBLE, v->y))
       ||(0 != TIO_put_var_section (grp, zname, &start, &count, NC_DOUBLE, v->z)))
     return -1;

   return 0;
}

static int process_file (char *iers_bulletin, const char *file)
{
   Vector_Array_Type *v = NULL;
   J2K_Type j2k = {0};
   int ncid, grp, time_dimid, start, count;
   size_t time_dimlen;
   int status = -1;

   if (0 != TIO_open (file, NC_WRITE, &ncid))
     return -1;

   if (0 != tio_use_file_epoch (ncid))
     return -1;

   if (0 != TIO_inq_grp (ncid, "anc_sec_102", &grp))
     goto return_status;

   if (0 != TIO_inq_dim (grp, "time", &time_dimid, &time_dimlen))
     goto return_status;

   if (NULL == (v = alloc_vecarray (time_dimlen)))
     goto return_status;

   start = 0;
   count = time_dimlen;

   if (0 != TIO_get_var_section (grp, "time", &start, &count, NC_DOUBLE, v->t))
     goto return_status;

   if (0 != j2k_open (&j2k, v->t[0], iers_bulletin))
     goto return_status;

   if (0 != convert_ecef_to_j2k (grp, &j2k, v, "anc_satx", "anc_saty", "anc_satz"))
     goto return_status;

   if (0 != convert_ecef_to_j2k (grp, &j2k, v, "anc_satvx", "anc_satvy", "anc_satvz"))
     goto return_status;

   status = 0;
return_status:
   j2k_close (&j2k);
   free_vecarray (v);
   (void) TIO_close (ncid);
   return status;
}

int main (int argc, char **argv)
{
   char *iers_bulletin;

   if (argc < 3)
     {
        fprintf (stderr, "\nUsage: ecef_to_j2k <iers-bulletin-A> FILE [FILE ...]\n\n");
        fprintf (stderr, "  Transform ephemeris vectors:\n");
        fprintf (stderr, "    anc_sec_102/anc_sat[xyz]\n");
        fprintf (stderr, "    anc_sec_102/anc_satv[xyz]\n");
        fprintf (stderr, "  from ECEF coordinate frame to J2000 coordinate frame.\n");
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
        if (0 != process_file (iers_bulletin, *argv))
          return 1;
        argv++;
        argc--;
     }
   while (argc > 0);

   fprintf (stdout, "done\n");

   return 0;
}

