MODULE time_module

  USE parameters_module
  USE error_module,     ONLY : ErrorType, CheckError

  IMPLICIT NONE
!   PRIVATE

  PUBLIC :: JULDAY
  PUBLIC :: CALDATE
  PUBLIC :: GET_DAY_OF_YEAR
  PUBLIC :: MINT
  
  TYPE TimeType
    REAL(KIND=8) :: Tau
    REAL(KIND=8) :: Tau0 ! ??
    REAL(KIND=8) :: JulianDay
    INTEGER      :: NYMD
    INTEGER      :: NHMS
    INTEGER      :: Year
    INTEGER      :: Month
    INTEGER      :: Day
    INTEGER      :: Hour
    INTEGER      :: Minute
    INTEGER      :: Second
    REAL(KIND=8) :: Hour_Decimal
    INTEGER      :: Season
    INTEGER      :: DayOfYear
    INTEGER      :: DayOfWeek
  ENDTYPE TimeType
  
  ! Astronomical Julian Date at 0 GMT, 1 Jan 1985
  REAL(KIND=8), PARAMETER :: JD85 = 2446066.5d0
  
  
CONTAINS

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_current_time
!
! !DESCRIPTION: Subroutine SET\_CURRENT\_TIME takes in the elapsed time in 
!  minutes since the start of a GEOS-Chem simulation and sets the GEOS-Chem 
!  time variables accordingly.  NOTE: All time variables are returned
!  w/r/t Greenwich Mean Time (aka Universal Time).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_TIME_FROM_TAU(Tau,Time,Error)
   
! !REMARKS:
!  The GEOS met fields are assimilated data, and therefore contain data on
!  the leap-year days.  However, the GCAP met fields are climatological GCM
!  output, and do not have data on the leap-year days.  SET_CURRENT_TIME
!  computes the days according to the Astronomical Julian Date algorithms
!  (in "julday_mod.f"), which contain leap-year days.  For GCAP, whenever a
!  February 29th is encountered, we shall just skip ahead a day to March 1st
!  and return the corresponding time & date values.
! 
! !REVISION HISTORY: 
!  05 Feb 2006 - R. Yantosca - Initial Version
!  (1 ) GCAP/GISS fields don't have leap years, so if JULDAY says it's 
!        Feb 29th, reset MONTH, DAY, JD1 to Mar 1st. (swu, bmy, 8/29/05)
!  (2 ) Now references "define.h".  Now add special handling to skip from
!        Feb 28th to Mar 1st for GCAP model. (swu, bmy, 4/24/06)
!  (3 ) Fix bug in case of GCAP fields for runs that start during leap year
!       and after February 29 (phs, 9/27/06)  
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  29 Jul 2011 - R. Yantosca - Bug fix: For GCAP, we need to skip over the
!                              # of leap-year-days that have already occurred
!                              when going from Julian date to Y/M/D date
!  14 Jun 2013 - R. Yantosca - Now move the day of week computation here
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! --------------------
      ! Subroutine Arguments
      ! --------------------
      REAL(KIND=8),    INTENT(IN)    :: Tau
      TYPE(TimeType),  INTENT(INOUT) :: Time
      TYPE(ErrorType), INTENT(INOUT) :: Error
      
      ! ---------------
      ! Local variables
      ! ---------------
      REAL(KIND=8)  :: A, B,  JD0, THISDAY

      !=================================================================
      ! SET_CURRENT_TIME starts here
      !=================================================================
      
      ! Check error status before computation
      IF(CheckError(Error)) RETURN

      ! Copy Tau
      Time%Tau = Tau
      
      ! Compute Julian Day
      Time%JulianDay = JD85 + Time%Tau / 24.0d0
      
      ! Floor it
      JD0 = REAL(FLOOR(Time%JulianDay),KIND=8)
      
      ! Compute Y/M/D
      CALL CALDATE( Time%JulianDay, Time%NYMD, Time%NHMS )
      
      
      ! Extract current year, month, day from NYMD
      CALL YMD_EXTRACT( Time%NYMD, Time%YEAR, Time%MONTH, Time%DAY )

      ! Extract current hour, minute, second from NHMS
      CALL YMD_EXTRACT( Time%NHMS, Time%HOUR, Time%MINUTE, Time%SECOND )
      
      ! Current Greenwich Mean Time
      Time%Hour_Decimal = ( DBLE( Time%HOUR )            ) &
                        + ( DBLE( Time%MINUTE ) / 60d0   ) &
                        + ( DBLE( Time%SECOND ) / 3600d0 )

      ! Season index (1=DJF, 2=MAM, 3=JJA, 4=SON)
      SELECT CASE ( Time%MONTH )
         CASE ( 12, 1, 2 )
            Time%SEASON = 1
         CASE ( 3, 4, 5 )
            Time%SEASON = 2
         CASE ( 6, 7, 8 )
            Time%SEASON = 3
         CASE ( 9, 10, 11 )
            Time%SEASON = 4
      END SELECT
      
      ! Days elapsed in this year (0-366)
      Time%DayOfYear = GET_DAY_OF_YEAR(Time%Month,Time%Day,Time%Year,&
                                       Time%Hour,Time%Minute,REAL(Time%Second,KIND=8))
      
      ! Get fractional GMT day
      THISDAY     = Time%DAY + ( Time%Hour_Decimal / 24d0 )

      ! Get current Julian date 
      !JD          = JULDAY( Time%YEAR, Time%MONTH, THISDAY )

      ! Add 1.5 to JD and divide by 7
      A           = ( Time%JulianDay + 1.5d0 ) / 7d0

      ! Take fractional part and multiply by 7
      B           = ( A - INT( A ) ) * 7d0
      B           = INT( NINT( B*1d5 + SIGN( 5d0, B ) ) / 10d0 ) / 1d4

      ! Discard the fractional part of B
      Time%DayOfWeek = INT( B )
      
      

  END SUBROUTINE SET_TIME_FROM_TAU

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymd_extract
!
! !DESCRIPTION: Subroutine YMD\_EXTRACT extracts the year, month, and date 
!  from an integer variable in YYYYMMDD format.  It can also extract the 
!  hours, minutes, and seconds from a variable in HHMMSS format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE YMD_EXTRACT( NYMD, Y, M, D )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: NYMD      ! YYYY/MM/DD format date
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: Y, M, D   ! Separated YYYY, MM, DD values
! 
! !REVISION HISTORY: 
!  21 Nov 2001 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(KIND=8) :: REM

      ! Extract YYYY from YYYYMMDD 
      Y = INT( DBLE( NYMD ) / 1d4 )

      ! Extract MM from YYYYMMDD
      REM = DBLE( NYMD ) - ( DBLE( Y ) * 1d4 )
      M   = INT( REM / 1d2 )

      ! Extract DD from YYYYMMDD
      REM = REM - ( DBLE( M ) * 1d2 )
      D   = INT( REM )

      ! Return to calling program
      END SUBROUTINE YMD_EXTRACT

! NAME:
!	JULDAY
!
! PURPOSE:
!	Calculate the Julian Day Number for a given month, day, and year.
!	This is the inverse of the library function CALDAT.
!	See also caldat, the inverse of this function.
!
! INPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year. Year parameters must be valid
!              values from the civil calendar. Years B.C.E. are represented
!              as negative integers. Years in the common era are represented
!              as positive integers. In particular, note that there is no
!              year 0 in the civil calendar. 1 B.C.E. (-1) is followed by
!              1 C.E. (1).
!
!	HOUR:	Number of the hour of the day.
!
!	MINUTE:	Number of the minute of the hour.
!
!	SECOND:	Number of the second of the minute (double precision number).
!
! OUTPUTS:
!	JULDAY returns the Julian Day Number (which begins at noon) of the
!	specified calendar date.
!
! RESTRICTIONS:
!	Accuracy using IEEE double precision numbers is approximately
!      1/10000th of a second, with higher accuracy for smaller (earlier)
!      Julian dates.
!
! MODIFICATION HISTORY:
!	Translated from "Numerical Recipies in C", by William H. Press,
!	Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
!	Cambridge University Press, 1988 (second printing).
!
!	AB, September, 1988
!	DMS, April, 1995, Added time of day.
!
!      Translated into Fortran90 by Vijay Natraj, JPL, February 24 2012

  function JULDAY(MONTH, DAY, YEAR, Hour, Minute, Second)
    
    implicit none
    
    !  Inputs
    
    integer(kind=4), intent(in)  :: MONTH
    integer(kind=4), intent(in)  :: DAY
    integer(kind=4), intent(in)  :: YEAR
    integer(kind=4), intent(in)  :: HOUR
    integer(kind=4), intent(in)  :: MINUTE
    real(kind=8),    intent(in)  :: SECOND
    
    !  Outputs
    
    real(kind=8)                 :: JULDAY
    
    !  Local variables
    
    integer(kind=4)              :: GREG
    integer(kind=4)              :: min_calendar
    integer(kind=4)              :: max_calendar
    integer(kind=4)              :: bc
    integer(kind=4)              :: L_YEAR
    integer(kind=4)              :: inJanFeb
    integer(kind=4)              :: JY
    integer(kind=4)              :: JM
    integer(kind=4)              :: JA
    integer(kind=4)              :: JUL
    real(kind=8)                 :: eps
    
    ! Gregorian Calendar was adopted on Oct. 15, 1582
    ! skipping from Oct. 4, 1582 to Oct. 15, 1582
    
    GREG = 2299171  ! incorrect Julian day for Oct. 25, 1582
    
    ! check if date is within allowed range
    
    min_calendar = -4716
    max_calendar = 5000000
    IF ((YEAR .LT. min_calendar) .OR. (YEAR .GT. max_calendar)) &
         write(0,*) 'Value of Julian date is out of allowed range'
    
    IF (YEAR .LT. 0) THEN
       bc = 1
    ELSE
       bc = 0
    ENDIF
    L_YEAR = YEAR + bc
    
    IF (MONTH .LE. 2) THEN
       inJanFeb = 1
    ELSE
       inJanFeb = 0
    ENDIF
    
    JY = YEAR - inJanFeb
    JM = MONTH + 1 + 12*inJanFeb
    
    JUL = FLOOR(365.25d0 * JY) + FLOOR(30.6001d0 * JM) + DAY + 1720995
    
    ! Test whether to change to Gregorian Calendar.
    
    IF (JUL .GE. GREG) THEN
       JA = FLOOR(0.01d0 * JY)
       JUL = JUL + 2 - JA + FLOOR(0.25d0 * JA)
    ENDIF
    
    ! Add a small offset so we get the hours, minutes, & seconds back correctly
    ! if we convert the Julian dates back. This offset is proportional to the
    ! Julian date, so small dates (a long, long time ago) will be "more" accurate.
    
    ! eps = (MACHAR(/DOUBLE)).eps
    
    eps = 2.2204460d-16 ! For Ganesha, calculated from IDL output of above statement
    
    IF (ABS(JUL) .GE. 1) THEN
       eps = eps*ABS(JUL)
    ENDIF
    
    ! For Hours, divide by 24, then subtract 0.5, in case we have unsigned integers.
    
    JULDAY = real(JUL,kind=8) + real(Hour,kind=8)/24.d0 - &
         0.5d0 + real(Minute,kind=8)/1440.d0 + &
         Second/86400.d0 + eps
    
    RETURN
  END FUNCTION JULDAY
  
  FUNCTION NYMD2TAU( NYMD ) RESULT ( TAU )
    !
    ! !INPUT PARAMETERS: 
    !
    INTEGER, INTENT(IN) :: NYMD
    !
    ! !RETURN VALUE:
    !
    REAL(KIND=8) :: TAU
    
    ! Local variables
    INTEGER :: Y, M, D
    REAL(KIND=8) :: JDAY0, JDAY1
    
    ! ==============================================================
    
    ! Extract date
    CALL YMD_EXTRACT( NYMD, Y, M, D )
    
    ! Compute Julday for 1985
    JDAY0 =  JULDAY(1, 1, 1985, 0, 0, 0.0d0)
    
    ! Compute current julian day
    JDAY1 = JULDAY(M, D, Y, 0, 0, 0.0d0)
    
    ! Compute tau time
    TAU = (JDAY1-JDAY0)*24.0d0
    
    RETURN
  END FUNCTION NYMD2TAU
  
  ! !DESCRIPTION: Function MINT is the modified integer function.
  FUNCTION MINT( X ) RESULT ( VALUE )
    !
    ! !INPUT PARAMETERS: 
    !
    REAL(KIND=8), INTENT(IN) :: X
    !
    ! !RETURN VALUE:
    !
    REAL(KIND=8) :: VALUE
    !
    ! !REMARKS:
    !  The modified integer function is defined as follows:
    !
    !            { -INT( ABS( X ) )   for X < 0
    !     MINT = {
    !            {  INT( ABS( X ) )   for X >= 0
    !
    ! !REVISION HISTORY: 
    !  20 Nov 2001 - R. Yantosca - Initial version
    !  20 Nov 2009 - R. Yantosca - Added ProTeX headers
    IF ( X < 0d0 ) THEN 
       VALUE = -INT( ABS( X ) )        
    ELSE
       VALUE =  INT( ABS( X ) )        
    ENDIF
    
  END FUNCTION MINT
  
  ! !DESCRIPTION: Subroutine CALDATE converts an astronomical Julian day to 
  !  the YYYYMMDD and HHMMSS format.
  SUBROUTINE CALDATE( JULIANDAY, YYYYMMDD, HHMMSS )
    !
    ! !INPUT PARAMETERS: 
    !
    ! Arguments
    REAL(KIND=8), INTENT(IN) :: JULIANDAY  ! Astronomical Julian Date 
    !
    ! !OUTPUT PARAMETERS: 
    !
    INTEGER, INTENT(OUT) :: YYYYMMDD   ! Date in YYYY/MM/DD format
    INTEGER, INTENT(OUT) :: HHMMSS     ! Time in hh:mm:ss format
    !
    ! !REMARKS:
    !   Algorithm taken from "Practical Astronomy With Your Calculator",
    !   Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
    !
    ! !REVISION HISTORY: 
    !  (1 ) Now compute HHMMSS correctly.  Also use REAL*4 variables HH, MM, SS
    !        to avoid roundoff errors. (bmy, 11/21/01)
    !  (2 ) Renamed NYMD to YYYYMMDD and NHMS to HHMMSS for documentation
    !        purposes (bmy, 6/26/02)
    !  20 Nov 2009 - R. Yantosca - Added ProTeX header
    !
    ! !LOCAL VARIABLES:
    !   
    REAL(KIND=4) :: HH, MM, SS
    REAL(KIND=8) :: A, B, C, D, DAY, E, F 
    REAL(KIND=8) :: FDAY, G, I, JD, M, Y
    
    !=================================================================
    ! CALDATE begins here!
    ! See "Practical astronomy with your calculator", Peter Duffett-
    ! Smith 1992, for an explanation of the following algorithm.
    !=================================================================
    JD = JULIANDAY + 0.5d0
    I  = INT( JD )
    F  = JD - INT( I )
    
    IF ( I > 2299160d0 ) THEN
       A = INT( ( I - 1867216.25d0 ) / 36524.25d0 )
       B = I + 1 + A - INT( A / 4 )
    ELSE
       B = I
    ENDIF
    
    C = B + 1524d0
    
    D = INT( ( C - 122.1d0 ) / 365.25d0 )
    
    E = INT( 365.25d0 * D )
    
    G = INT( ( C - E ) / 30.6001d0 )
    
    ! DAY is the day number
    DAY  = C - E + F - INT( 30.6001d0 * G ) 
    
    ! FDAY is the fractional day number
    FDAY = DAY - INT( DAY )
    
    ! M is the month number
    IF ( G < 13.5d0 ) THEN
       M = G - 1d0
    ELSE
       M = G - 13d0
    ENDIF
    
    ! Y is the year number
    IF ( M > 2.5d0 ) THEN
       Y = D - 4716d0
    ELSE
       Y = D - 4715d0
    ENDIF
    
    ! Year-month-day value
    YYYYMMDD = ( INT( Y ) * 10000 ) + ( INT( M ) * 100 ) + INT( DAY )
    
    ! Hour-minute-second value
    ! NOTE: HH, MM, SS are REAL*4 to avoid numerical roundoff errors
    HH     = REAL(FDAY * 24d0,KIND=4)
    MM     = REAL(( HH - INT( HH ) ) * 60d0,KIND=4)
    SS     = REAL(( MM - INT( MM ) ) * 60d0,KIND=4)
    HHMMSS = ( INT( HH ) * 10000 ) + ( INT( MM ) * 100 ) + INT( SS )
    
  END SUBROUTINE CALDATE
  
  FUNCTION GET_DAY_OF_YEAR(MONTH,DAY,YEAR,HOUR,MINUTE,SECOND) RESULT( DOY )
    
    integer(kind=4), intent(in)  :: MONTH
    integer(kind=4), intent(in)  :: DAY
    integer(kind=4), intent(in)  :: YEAR
    integer(kind=4), intent(in)  :: HOUR
    integer(kind=4), intent(in)  :: MINUTE
    real(kind=8),    intent(in)  :: SECOND

    !
    ! !RETURN VALUE:
    ! 
    INTEGER :: DOY   ! Astronomical Julian Date       
        
    ! Compute DOY using difference in julian day
    DOY = INT(JULDAY(MONTH,DAY,YEAR,HOUR,MINUTE,SECOND) - JULDAY(1,0,YEAR,0,0,0d0) )
    
    
  END FUNCTION GET_DAY_OF_YEAR
  
  
  
   


END MODULE time_module