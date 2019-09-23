#include <stdio.h>
#include "met.h"

int main (int argc, char **argv)
{
   Met_List_Type *met = NULL;

   float isobars[] = {90000.0, 9000.0, 900.0};
#define NLEV sizeof(isobars)/sizeof(*isobars)

   float temperature_on_isobar[NLEV];
   Met_Value_Type mvt =
     {
        .pressure_surface = 0.0,
        .pressure_tropopause = 0.0,
        .isobars = isobars,
        .num_isobars = NLEV,
        .temperature_on_isobar = temperature_on_isobar
     };

   struct points_type
     {
        float lon;
        float lat;
     }
   points[] =
     {
        {214.5, 1.0}            // lower left corner of NAM 221 grid
        , {-90.0, 36.0}         // Roughly near the center of the TEMPO FOR
        , {-52.712, 47.5605}    // St. John's, Newfoundland
        , {-124.4049, 40.4370}  // Capetown, CA
     };

   unsigned int flags;
   size_t i, n;
   int k, j, status;

   if (argc == 1)
     {
        fprintf (stderr, "Usage: %s FILE [FILE ...]\n", argv[0]);
        return 1;
     }

   flags = 0;
   flags |= MET_READ_PRESSURE_SURFACE;
   flags |= MET_READ_PRESSURE_TROPOPAUSE;
   flags |= MET_READ_TEMPERATURE_ON_ISOBARS;

   if (NULL == (met = met_list_new (flags)))
     return 1;

   for (j = 1; j < argc; j++)
     {
        if (0 != met_list_add_file (met, argv[j]))
          {
             met_list_free (met);
             return 1;
          }
     }

   n = sizeof(points)/sizeof(points[0]);
   fprintf (stdout, "Will test n = %ld points\n", n);

   for (i = 0; i < n; i++)
     {
        struct points_type *p = &points[i];

        fprintf (stdout, "point[%ld] lon=%f lat=%f\n", i, p->lon, p->lat);

        status = met_list_interp (met, p->lon, p->lat, &mvt);

        if (status < 0)
          {
             fprintf (stderr, "*** Error: mft_interp failed (status=%d)\n", status);
             return 1;
          }
        else if (status == MFT_INTERP_DOMAIN_ERROR)
          {
             fprintf (stdout, "==> mft_interp: domain error\n");
             continue;
          }

        fprintf (stdout, "P(surf)=%f hPa  P(tropopause)=%f hPa\n",
                 mvt.pressure_surface, mvt.pressure_tropopause);
        for (k = 0; k < mvt.num_isobars; k++)
          {
             fprintf (stdout, "%d: p_isobar: %f   T_isobar: %f\n",
                      k, mvt.isobars[k], mvt.temperature_on_isobar[k]);
          }
     }

   met_list_free (met);

   return 0;
}
