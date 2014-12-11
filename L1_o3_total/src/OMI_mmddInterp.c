#include <cfortran.h>

static char daytab[2][13] = { 
      { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
      { 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31} };

static short bindoy[2][13] = { 
      { 0, 15, 46, 74,105,135,166,196,227,258,288,319,349},
      { 0, 15, 46, 75,106,136,167,197,228,259,289,320,350} };


/* day of year: set day of year from month and day */
int 
day_of_year( int year, 
             int month,
             int day )
{
   int i, leap;

   leap = (year%4 == 0 && year%100 != 0 || year%400 == 0);

   for( i =1; i < month; i++ )
       day += daytab[leap][i];
   return day;
}


float 
OMI_mmddInterp( int  year, 
                int  month,
                int  day,
                int *dday,
                int *mm_cur,
                int *mm_pre )
{  
   int   doy, leap;
   float frac;

   doy = day_of_year( year, month, day );
   leap = (year%4 == 0 && year%100 != 0 || year%400 == 0);
   if( doy <= bindoy[leap][month] )
   { *mm_cur = month; 
     *mm_pre = month-1; 
      if( *mm_pre < 1 ) *mm_pre = 12;
   }
   else
   { *mm_cur = month+1; 
      if( *mm_cur > 12 ) *mm_cur = 1;
     *mm_pre = month; 
   }

   *dday = bindoy[leap][*mm_cur] - doy; 
   
   if( *mm_pre == 12 && doy > bindoy[leap][12] && doy <= bindoy[leap][12]+16) 
      *dday = 31 - (doy - bindoy[leap][12]);

   frac = (float) *dday/ 30.5;
   return frac;
} 

/* FORTRAN bindings */

FCALLSCFUN3( INT, day_of_year, DAY_OF_YEAR, day_of_year, INT, INT, INT )
