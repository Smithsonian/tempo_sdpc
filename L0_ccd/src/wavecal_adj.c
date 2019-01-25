/* -*- mode: C; mode: fold -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include "config.h"
#include "util.h"
#include "wavecal_adj.h"

typedef struct
{
   Wadj_Cheb_Series_Type attr;
   double *coeff;
}
Table_Type;

struct Wadj_Type
{
   Table_Type full;      /* Chebyshev series coefficients for full wavelength band */
   Table_Type narrow;    /* Chebyshev series coefficients for narrow wavelength band being fitted */
   Table_Type adjust;    /* Chebyshev series coefficients for narrow-band wavelength shift adjustment */
   size_t num_xtrack;    /* number of pixels along the slit */
};

static void table_free (Table_Type *tbl)
{
   FREE(tbl->coeff);
   tbl->coeff = NULL;
}

static int table_alloc (Table_Type *tbl, size_t num_xtrack)
{
   double *coeff = NULL;
   size_t len = num_xtrack * tbl->attr.num_series_coeff;

   if (NULL == (coeff = (double *)MALLOC (len * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   memset ((char *)coeff, 0, len * sizeof(double));
   tbl->coeff = coeff;

   return 0;
}

static int table_read (int grp, const char *name, Table_Type *tbl)
{
   TIO_Var_Info_Type info = {0};
   Wadj_Cheb_Series_Type *attr = &tbl->attr;
   int start[2], count[2];
   size_t num_xtrack;

   if (0 != TIO_inq_var (grp, name, &info))
     return -1;

   if ((0 != TIO_get_att (grp, info.varid, "num_params", NC_INT, &attr->num_series_coeff))
       ||(0 != TIO_get_att (grp, info.varid, "pix_min", NC_INT, &attr->pix_min))
       ||(0 != TIO_get_att (grp, info.varid, "pix_max", NC_INT, &attr->pix_max)))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: error reading attributes of variable: %s", __func__, name);
        return -1;
     }

   if (attr->num_series_coeff == 0)
     return 0;

   num_xtrack = info.dimlens[0];

   if (0 != table_alloc (tbl, num_xtrack))
     return -1;

   start[0] = 0;
   start[1] = 0;

   count[0] = num_xtrack;
   count[1] = attr->num_series_coeff;

   if (0 != TIO_get_var_section (grp, name, start, count, NC_DOUBLE, tbl->coeff))
     {
        table_free (tbl);
        return -1;
     }

   return 0;
}

void wadj_close (Wadj_Type *wadj)
{
   if (wadj == NULL)
     return;
   table_free (&wadj->full);
   table_free (&wadj->narrow);
   table_free (&wadj->adjust);
   FREE(wadj);
}

static int read_table_file (config_t *cfg, char **file)
{
   config_setting_t *s;
   const char *path = NULL;

   if (NULL == (s = config_lookup (cfg, "wavecal_control")))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: accessing wavecal_control in param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (CONFIG_TRUE != config_setting_lookup_string (s, "shift_adjust_table", &path))
     {
        tell_verror (TELL_INVALID_PARM_ERROR,
                     "%s: reading shift_adjust_table from param file: %s",
                     __func__, config_error_file (cfg));
        return -1;
     }

   if (NULL == (*file = expand_path (path)))
     return -1;

   return 0;
}

static int invalid_tables (const Table_Type *full, const Table_Type *adj)
{
   const Wadj_Cheb_Series_Type *full_attr = &full->attr;
   const Wadj_Cheb_Series_Type *adj_attr = &adj->attr;

   if ((full_attr->pix_min != adj_attr->pix_min)
       || (full_attr->pix_max != adj_attr->pix_max))
     {
        tell_verror (TELL_RUNTIME_ERROR,
                     "%s: inconsistent window parameters in shift adjust table (mismatched adjustment coefficients)",
                     __func__);
        return -1;
     }

   return 0;
}

Wadj_Type *wadj_open (config_t *cfg, const char *group)
{
   Wadj_Type *wadj = NULL;
   char *file = NULL;
   int ncid, grp, dimid_xtrack;
   int status = -1;

   if (0 != read_table_file (cfg, &file))
     return NULL;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     goto return_error;

   if (NULL == (wadj = (Wadj_Type *) MALLOC (sizeof *wadj)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto return_error;
     }
   memset ((char *)wadj, 0, sizeof (*wadj));

   if (0 != TIO_inq_dim (ncid, "xtrack", &dimid_xtrack, &wadj->num_xtrack))
     goto return_error;

   if (0 != TIO_inq_grp (ncid, group, &grp))
     goto return_error;

   if ((0 != table_read (grp, "full_band", &wadj->full))
       || (0 != table_read (grp, "narrow_band", &wadj->narrow))
       || (0 != table_read (grp, "full_band_adjust", &wadj->adjust)))
     goto return_error;

   if (invalid_tables (&wadj->full, &wadj->adjust))
     goto return_error;

   status = 0;
return_error:
   if (status)
     {
        wadj_close (wadj);
        wadj = NULL;
     }
   TIO_close (ncid);
   FREE(file);
   return wadj;
}

int wadj_narrow_band_get_attr (const Wadj_Type *wadj, Wadj_Cheb_Series_Type *cheb)
{
   if ((wadj == NULL) || (cheb == NULL))
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: received NULL pointer", __func__);
        return -1;
     }

   *cheb = wadj->narrow.attr;
   return 0;
}

int wadj_full_band_get_attr (const Wadj_Type *wadj, Wadj_Cheb_Series_Type *cheb)
{
   if ((wadj == NULL) || (cheb == NULL))
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: received NULL pointer", __func__);
        return -1;
     }

   *cheb = wadj->full.attr;
   return 0;
}

const double *wadj_narrow_band_coeff (const Wadj_Type *wadj, size_t xtrack)
{
   const Table_Type *n;
   const double *coeff;

   if (wadj == NULL)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: received NULL pointer", __func__);
        return NULL;
     }

   if (xtrack >= wadj->num_xtrack)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: invalid xtrack=%ld (expected range = [0,%ld])",
                     __func__, xtrack, wadj->num_xtrack);
        return NULL;
     }

   n = &wadj->narrow;

   coeff = n->coeff + n->attr.num_series_coeff * xtrack;

   return coeff;
}

int wadj_num_final_coeff (const Wadj_Type *wadj, int *num_coeff)
{
   int num_full, num_adj;

   if ((wadj == NULL) || (num_coeff == NULL))
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: received NULL pointer", __func__);
        return -1;
     }

   num_full = wadj->full.attr.num_series_coeff;
   num_adj = wadj->adjust.attr.num_series_coeff;

   *num_coeff = (num_full > num_adj) ? num_full : num_adj;

   return 0;
}

/* code assumes that the coeff array is large enough (use wadj_num_adj_coeff to compute the size!!) */
int wadj_final_coeff (const Wadj_Type *wadj, size_t xtrack, double shift, double *coeff)
{
   const Table_Type *full;
   const Table_Type *adj;
   double *full_coeff_x;
   double *adj_coeff_x;
   size_t xtrack_offset;
   int i;

   if ((wadj == NULL) || (coeff == NULL))
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: received NULL pointer", __func__);
        return -1;
     }

   if (xtrack >= wadj->num_xtrack)
     {
        tell_verror (TELL_INTERNAL_ERROR, "%s: invalid xtrack=%ld (expected range = [0,%ld])",
                     __func__, xtrack, wadj->num_xtrack);
        return -1;
     }

   full = &wadj->full;
   xtrack_offset = xtrack * full->attr.num_series_coeff;
   full_coeff_x = full->coeff + xtrack_offset;

   for (i = 0; i < full->attr.num_series_coeff; i++)
     {
        coeff[i] = full_coeff_x[i];
     }

   adj = &wadj->adjust;

   if (adj->attr.num_series_coeff == 0)
     return 0;

   /* Here, we're assuming that the adjustment Chebyshev polynomial terms have
    * the same 'x' coordinate as the full-band series Chebyshev polynomial terms.
    * This is true when both series are derived for the same pixel interval.
    * The number of terms in each the series need not be the same.
    * The consistency of these two Chebyshev series expansions is assumed to
    * have been checked when the lookup table was loaded, so there's no need to
    * repeat that check here.
    */

   adj_coeff_x = adj->coeff + xtrack_offset;

   for (i = 0; i < adj->attr.num_series_coeff; i++)
     {
        coeff[i] += shift * adj_coeff_x[i];
     }

   return 0;
}
