#include "cfortHdf.h"

/* static float           fill_float32 = -0X1P+100;
   static double          fill_float64 = -0X1P+100; 
*/
static float           fill_float32 = -0X1P100;
static double          fill_float64 = -0X1P100;

int 
r4Fill( float *r4 )
{
  *r4 = fill_float32;
  return 0;
}

int 
r8Fill( double *r8 )
{
  *r8 = fill_float64;
  return 0;
}

FCALLSCFUN1(INT, r4Fill, R4FILL, r4fill, PFLOAT )
FCALLSCFUN1(INT, r8Fill, R8FILL, r8fill, PDOUBLE )
