#include <stdio.h>
#include <stdlib.h>
#include "met.h"

#define NUM_ISOBARS 5
#define BUFSIZE     128

static int perform_test (Met_List_Type *met, Met_Value_Type *mvt,
                         int argc, char **argv)
{
   FILE *fp = NULL;
   const char *points_file;
   float isobars[NUM_ISOBARS];
   float temperature_on_isobar[NUM_ISOBARS];
   int k, j, status = -1;

   if (argc < 3)
     {
        fprintf (stderr, "Usage: %s POINTS_FILE GRIB_FILE [GRIB_FILE ...]\n", argv[0]);
        return 0;
     }

   argv++; argc--;
   points_file = argv[0];
   argv++; argc--;

   for (j = 0; j < argc; j++)
     {
        if (0 != met_list_add_file (met, argv[j]))
          {
             return -1;
          }
     }

   if (NULL == (fp = fopen (points_file, "r")))
     {
        fprintf (stderr, "*** Error: opening file %s for reading", points_file ? points_file : "(null)");
        return -1;
     }

   while (0 == feof(fp))
     {
        float lon, lat, p0, dp;
        char buf[BUFSIZE];

        /* This isn't bulletproof, but it doesn't need to be, it's only for testing */
        if (NULL == fgets (buf, sizeof(buf), fp))
          break;

        /* allow comment lines */
        if (buf[0] == '#')
          continue;

        if (2 != sscanf (buf, "%f %f", &lon, &lat))
          {
             fprintf (stderr, "*** Error: reading (lon,lat) coordinates from file: %s\n", points_file);
             goto cleanup_and_return;
          }

        mvt->num_isobars = 0;
        mvt->isobars = NULL;
        mvt->temperature_on_isobar = NULL;

        status = met_list_interp (met, lon, lat, mvt);

        if (status < 0)
          {
             fprintf (stderr, "*** Error: met_list_interp failed (status=%d)\n", status);
             goto cleanup_and_return;
          }
        else if (status == MFT_INTERP_DOMAIN_ERROR)
          {
             fprintf (stdout, "==> mft_interp: domain error\n");
             continue;
          }

        mvt->num_isobars = NUM_ISOBARS;
        mvt->isobars = isobars;
        mvt->temperature_on_isobar = temperature_on_isobar;

        /* the forecast isobar grid does not extend all the way to the surface */
        p0 = 0.95 * mvt->pressure_surface;
        dp = (mvt->pressure_tropopause - mvt->pressure_surface)/(mvt->num_isobars - 1);

        for (k = 0; k < mvt->num_isobars; k++)
          {
             mvt->isobars[k] = p0 + k * dp;
             mvt->temperature_on_isobar[k] = 0.0;
          }

        status = met_list_interp (met, lon, lat, mvt);

        if (status < 0)
          {
             fprintf (stderr, "*** Error: met_list_interp failed (status=%d)\n", status);
             goto cleanup_and_return;
          }
        else if (status == MFT_INTERP_DOMAIN_ERROR)
          {
             fprintf (stdout, "==> met_list_interp: domain error\n");
             goto cleanup_and_return;
          }

        fprintf (stdout, "%0.1f\n", mvt->pressure_tropopause);
        fprintf (stdout, "%0.1f\n", mvt->pressure_surface);
        for (k = 0; k < mvt->num_isobars; k++)
          {
             fprintf (stdout, "%6.1f\n", mvt->temperature_on_isobar[k]);
          }
        fflush (stdout);
     }

   status = 0;
cleanup_and_return:
   (void) fclose (fp);

   return status;
}

int main (int argc, char **argv)
{
   Met_List_Type *met = NULL;

   Met_Value_Type mvt =
     {
        .pressure_surface = 0.0,
        .pressure_tropopause = 0.0,
        .isobars = NULL,
        .num_isobars = 0,
        .temperature_on_isobar = NULL
     };
   unsigned int flags;
   int status;

   flags = 0;
   flags |= MET_READ_PRESSURE_SURFACE;
   flags |= MET_READ_PRESSURE_TROPOPAUSE;
   flags |= MET_READ_TEMPERATURE_ON_ISOBARS;

   if (NULL == (met = met_list_new (flags)))
     return 1;

   status = perform_test (met, &mvt, argc, argv);

   met_list_free (met);

   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
