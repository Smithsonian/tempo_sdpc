#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <tell.h>
#include <tio.h>

#include "config.h"

#define SNOW_TYPE_PRIVATE_DATA \
   unsigned char *map; \
   size_t num_rows; \
   size_t num_cols;
#include "snow.h"
#include "ezlhconv.h"

#define NISE_MASK_NODATA 255

#define NISE_PROJECTION_NAME  "Nl"
/* NISE projection names have the form [NSM][lh]
 * where N=northern hemisphere
 *       S=southern hemisphere
 *       M=global
 *       l = "low" = 25km resolution
 *       h = "high" = 12.5km resolution
 */

static void free_snow_type (Snow_Type *sn)
{
   if (sn == NULL)
     return;
   FREE(sn->map);
   FREE(sn);
}

static int snow_lookup (const Snow_Type *sn, unsigned int num,
                        const float *lon, const float *lat,
                        unsigned char *mask)
{
   char grid[] = NISE_PROJECTION_NAME;
   unsigned int i;

   for (i = 0; i < num; i++)
     {
        double col, row;
        if (0 == ezlh_convert (grid, lat[i], lon[i], &col, &row))
          {
             int k = (int) col + ((int) row) * sn->num_cols;
             mask[i] = sn->map[k];
          }
        else mask[i] = NISE_MASK_NODATA;
     }

   return 0;
}

static Snow_Type *new_snow_type (void)
{
   Snow_Type *sn = NULL;

   if (NULL == (sn = (Snow_Type *)MALLOC (sizeof *sn)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }
   memset ((char *)sn, 0, sizeof *sn);

   sn->sn_delete = free_snow_type;
   sn->sn_lookup = snow_lookup;

   return sn;
}

static int alloc_map (Snow_Type *sn, size_t *dimlens)
{
   size_t map_size;

   sn->num_rows = dimlens[0];
   sn->num_cols = dimlens[1];

   map_size = dimlens[0] * dimlens[1] * sizeof (unsigned char);

   if (NULL == (sn->map = (unsigned char *)MALLOC (map_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   return 0;
}

static int read_snow_file (Snow_Type *sn, const char *file)
{
   int dimids[TIO_MAX_VAR_DIMS];
   size_t dimlens[2];
   int i, ncid, err, varid, status = -1;

   /* FIXME:  netcdf bug? 24 byte memory leak when opening,
    * and then closing an HDF4 file */
   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   /* Here, we're using the netcdf interface to read an HDF4 file.
    *
    * I'd prefer to use libtio for this, but the libtio interface assumes
    * that variables within a given group have unique names, and that isn't
    * true for these snow and ice coverage maps(!).  Because the variables
    * in the files don't have unique names, it's necessary to access each
    * variable in the file using its numerical id.
    *
    * Maybe it's ok/expected to assume that each variable will always
    * have the same numerical id, but I'd rather not assume that.
    * Since the names aren't unique and I don't want to assume
    * a fixed numerical id, I'll look at the variables in sequence
    * until I find the one that "looks right".  (sigh)
    */

   varid = 0;

   for (;;)
     {
        char name[NC_MAX_NAME];
        int vartype, num_dims, num_attrs;

        if (NC_NOERR != (err = nc_inq_var (ncid, varid, name, &vartype, &num_dims, dimids, &num_attrs)))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: accessing varid=%d, file=%s (%s)",
                          __func__, varid, file, nc_strerror (err));
             goto close_and_return;
          }

        /* Look for a variable with the expected name, data type, and dimensionality */
        if (! ((num_dims == 2)
               && (vartype == NC_UBYTE)
               && (0 == strcasecmp (name, "extent"))))
          continue;

        /* Check the coordinate variable names to make sure this is the right hemisphere */
        if (NC_NOERR != (err = nc_inq_dimname (ncid, dimids[0], name)))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: accessing dimid=%d (%s)",
                          __func__, dimids[0], nc_strerror (err));
             goto close_and_return;
          }

        /* found it? */
        if (NULL != strstr (name, "Northern"))
          break;

        varid++;
     }

   /* Now that we have (varid, dimids) for the Northern hemisphere map,
    * we can allocate space and read it: */

   for (i = 0; i < 2; i++)
     {
        if (NC_NOERR != (err = nc_inq_dimlen (ncid, dimids[i], &dimlens[i])))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading dimension lengths", __func__);
             goto close_and_return;
          }
     }

   if (0 != alloc_map (sn, dimlens))
     goto close_and_return;

   if (NC_NOERR != (err = nc_get_var_ubyte (ncid, varid, sn->map)))
     {
        tell_verror (TELL_IO_READ_ERROR,
                     "%s: reading snow and ice coverage from %s (%s)",
                     __func__, file, nc_strerror (err));
        goto close_and_return;
     }

   status = 0;
close_and_return:
   TIO_close (ncid);

   return status;
}

Snow_Type *snow_init (const char *file)
{
   Snow_Type *sn = NULL;

   if (NULL == (sn = new_snow_type ()))
     return NULL;

   if (0 != read_snow_file (sn, file))
     {
        free_snow_type (sn);
        return NULL;
     }

   return sn;
}
