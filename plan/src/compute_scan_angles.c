/** @file compute_scan_angles.c
 *  @brief Interface to INR library function calculateScanFG
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <libconfig.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

#include <tempo_geo.h>

#include "compute_scan_angles.h"

#define DEGTORAD       (M_PI/180.0)
#define DEGTOMICRORAD  (1.e6*DEGTORAD)

/* #define USE_PREFLIGHT_CONFIG 1 */
#undef USE_PREFLIGHT_CONFIG

/* The following parameter values should be set to the flight defaults,
 * and should match the indicated parameter values in the [satellite]
 * section of the operational INR config file. These parameters should
 * be consistent with a nominal flight boresight point of:
 *   boresight: {lon = -89.2170, lat = 33.5231};
 * Boresight coordinates from 2023-10-30 email to JCH from
 * Jim Carr <jcarr@carrastro.com>
 */
static Geometry_Param_Type Geometry_Params =
{
   .ewbias_rad    = 3.599225,          /* 'ewbias' config file parameter */
   .nsbias_rad    = 0.09479,           /* 'nwbias' config file parameter */
   .clockbias_rad = 110.0e-6,          /* 'clockingbias' config file parameter */
   .theta0_rad    = 13.0 * DEGTORAD    /* 'telescopeOffset' config file parameter, converted to radians */
};

/* This function enables overriding the hard-coded defaults,
 * mostly so that pre-flight regression tests run correctly.
 */
int __set_geometry_params (double ewbias, double nsbias, double clockingbias, double telescopeOffset)
{
   Geometry_Param_Type *g = &Geometry_Params;
   g->ewbias_rad    = ewbias;
   g->nsbias_rad    = nsbias;
   g->clockbias_rad = clockingbias;
   g->theta0_rad    = telescopeOffset * DEGTORAD;
   return 0;
}

static int geometry_params (double sat_lon_deg, Geometry_Param_Type *p)
{
#ifndef USE_PREFLIGHT_CONFIG
   (void) sat_lon_deg;
   *p = Geometry_Params;  /* struct copy */
#else
   double ratio = 6.610702780451408;  /* (GEO orbit radius)/(Earth equatorial radius) */
   double dlon = (-57.32 - 0.6288 * sat_lon_deg) * DEGTORAD;

   /* From TEMPO-SER-4008_TEMPO_Instrument_Geometric_Model_for_INR.pdf
    * theta0 = mirror incidence angle
    *        = 13.0 deg  (nominal)
    * Define:  nsbias = 5.53 rad   (nominal)
    *          ewbias = \pi + 2*theta0 + C
    *                 = 3.595378259108319  (nominal, for C=0)
    *
    *   lon_target = -57.32 + 0.3712 * lon_sat [deg]
    *   C = clock angle
    *     = arctan (Re * sin(lon_target - lon_sat) / (a0 - Re*cos(lon_target - lon_sat)))
    *   lon_sat = nominal satellite longitude [deg]
    *   Re = equatorial radius of Earth
    *   a0 = nominal satellite semi-major axis
    *
    * Simplifying, we get:
    *    C = arctan (sin(dlon) / (ratio - cos(dlon)))
    * where
    *    dlon = (lon_target - lon_sat)
    *         = (-57.32 - 0.6288 * lon_sat)
    *    ratio = a0/Re = (42163.968 / 6378.1370) = 6.610702780451408
    */

   /* These values below are derived from the nominal equations,
    * but for operations we may have different (off-nominal) values.
    * The operational values can be derived from the following
    * INRSW config file parameters:
    * [satellite]
    *     ewbias
    *     nsbias
    *     telescopeOffset
    * There's also a clock angle parameter in the INRSW config file,
    * but I don't know the param name.
    */
   p->theta0_rad = 13.0 * DEGTORAD;
   p->nsbias_rad = 5.53 * DEGTORAD;
   p->clockbias_rad = atan2 (sin(dlon), ratio - cos(dlon));
   p->ewbias_rad = M_PI + 2*p->theta0_rad + p->clockbias_rad;
#endif

   return 0;
}

int __compute_scan_angles (const EarthPoint *c_pt, double sat_lon, AziElev_Type *apt)
{
   EarthPoint pt = *c_pt;  /* struct copy */
   EarthPoint *fg_pts = NULL;
   EarthPolygon polygon =
     {
        .thePointCount = 1,
        .theAllocatedPoints = 1,
        .thePoints = &pt
     };
   EarthPolygon fg_polygon = {0};
   Geometry_Param_Type g = {0};

   /* need sat_lon in degrees */
   sat_lon /= DEGTORAD;

   geometry_params (sat_lon, &g);

   /* WARNING: This subroutine call modifies polygon.thePoints values!! */
   calculateScanFG (&polygon, sat_lon,
                    g.ewbias_rad, g.nsbias_rad, g.clockbias_rad, g.theta0_rad,
                    &fg_polygon);

   /* fg_polygon values are in degrees,
    * we want values in microradians */
   fg_pts = fg_polygon.thePoints;
#define SCANFG_HAS_COORDINATES_SWAPPED 1
#ifndef SCANFG_HAS_COORDINATES_SWAPPED
   apt->azimuth = fg_pts->theLon * DEGTOMICRORAD;
   apt->elevation = fg_pts->theLat * DEGTOMICRORAD;
#else
   apt->azimuth = fg_pts->theLat * DEGTOMICRORAD;
   apt->elevation = fg_pts->theLon * DEGTOMICRORAD;
#endif

   free (fg_polygon.thePoints);

   return 0;
}
