MODULE UTIL_tools_class
   IMPLICIT NONE
   CHARACTER (LEN=26), PARAMETER :: ALPHABETu = &
                       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
   CHARACTER (LEN=26), PARAMETER :: alphabetL = &
                       "abcdefghijklmnopqrstuvwxyz"

!   PUBLIC :: slit
!   PUBLIC :: Itrapez
   PUBLIC :: hunt
   PUBLIC :: adjustDEG
   PUBLIC :: toLowerString
!   PUBLIC :: toUpperString
   PUBLIC :: deQuote

!   INTERFACE  locate
!     MODULE PROCEDURE locateIK4
!     MODULE PROCEDURE locateRK4
!     MODULE PROCEDURE locateRK8
!   END INTERFACE  

!   INTERFACE  Lintext
!     MODULE PROCEDURE LintextRK4p
!     MODULE PROCEDURE LintextRK4a
!     MODULE PROCEDURE LintextRK8p
!     MODULE PROCEDURE LintextRK8a
!   END INTERFACE  

   CONTAINS

!   FUNCTION locateRK8( xx, x ) RESULT(locate)
!     IMPLICIT NONE
!     REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: xx
!     REAL (KIND = 8), INTENT(IN) :: x
!     INTEGER (KIND = 4) :: locate
!     !!
!      ! Given an array xx(1:N), and given a value x, return a
!      ! value j, such that x is between xx(j) and xx(j+1). xx
!      ! must be montonic, either increasing or decreasing. j = 0
!      ! or j = N is returned to indicate that x is out of range.
!     !!
!     INTEGER (KIND = 4) :: n, jl, jm, ju
!     LOGICAL :: ascnd
!
!     n = size(xx)
!     ascnd = ( xx(n) >= xx(1) ) ! True is ascending order of table
!                                ! False otherwise
!     jl = 0                     ! initialize lower
!     ju = n + 1                 ! and upper limits 
!
!     DO
!        IF( ju - jl <= 1 ) EXIT ! Repeat until this condition is satisfied
!        jm = (ju + jl)/2
!        IF( ascnd .EQV. (x >= xx(jm)) ) THEN
!           jl = jm 
!        ELSE
!           ju = jm
!        ENDIF
!     END DO
!
!     IF( x == xx(1) ) THEN      ! set the output, becareful with 
!        locate = 1              ! the end points
!     ELSE IF ( x == xx(n) ) THEN 
!        locate = n - 1
!     ELSE
!        locate = jl
!     ENDIF
!   END FUNCTION locateRK8
!
!   FUNCTION locateRK4( xx, x ) RESULT(locate)
!     IMPLICIT NONE
!     REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: xx
!     REAL (KIND = 4), INTENT(IN) :: x
!     INTEGER (KIND = 4) :: locate
!     !!
!      ! Given an array xx(1:N), and given a value x, return a
!      ! value j, such that x is between xx(j) and xx(j+1). xx
!      ! must be montonic, either increasing or decreasing. j = 0
!      ! or j = N is returned to indicate that x is out of range.
!     !!
!     INTEGER (KIND = 4) :: n, jl, jm, ju
!     LOGICAL :: ascnd
!
!     n = size(xx)
!     ascnd = ( xx(n) >= xx(1) ) ! True is ascending order of table
!                                ! False otherwise
!     jl = 0                     ! initialize lower
!     ju = n + 1                 ! and upper limits 
!
!     DO
!        IF( ju - jl <= 1 ) EXIT ! Repeat until this condition is satisfied
!        jm = (ju + jl)/2
!        IF( ascnd .EQV. (x >= xx(jm)) ) THEN
!           jl = jm 
!        ELSE
!           ju = jm
!        ENDIF
!     END DO
!
!     IF( x == xx(1) ) THEN      ! set the output, becareful with 
!        locate = 1              ! the end points
!     ELSE IF ( x == xx(n) ) THEN 
!        locate = n - 1
!     ELSE
!        locate = jl
!     ENDIF
!   END FUNCTION locateRK4
!
!   FUNCTION locateIK4( xx, x ) RESULT(locate)
!     IMPLICIT NONE
!     INTEGER (KIND = 4), DIMENSION(:), INTENT(IN) :: xx
!     INTEGER (KIND = 4), INTENT(IN) :: x
!     INTEGER (KIND = 4) :: locate
!     !!
!      ! Given an array xx(1:N), and given a value x, return a
!      ! value j, such that x is between xx(j) and xx(j+1). xx
!      ! must be montonic, either increasing or decreasing. j = 0
!      ! or j = N is returned to indicate that x is out of range.
!     !!
!     INTEGER (KIND = 4) :: n, jl, jm, ju
!     LOGICAL :: ascnd
!
!     n = size(xx)
!     ascnd = ( xx(n) >= xx(1) ) ! True is ascending order of table
!                                ! False otherwise
!     jl = 0                     ! initialize lower
!     ju = n + 1                 ! and upper limits 
!
!     DO
!        IF( ju - jl <= 1 ) EXIT ! Repeat until this condition is satisfied
!        jm = (ju + jl)/2
!        IF( ascnd .EQV. (x >= xx(jm)) ) THEN
!           jl = jm 
!        ELSE
!           ju = jm
!        ENDIF
!     END DO
!
!     IF( x == xx(1) ) THEN      ! set the output, becareful with 
!        locate = 1              ! the end points
!     ELSE IF ( x == xx(n) ) THEN 
!        locate = n - 1
!     ELSE
!        locate = jl
!     ENDIF
!   END FUNCTION locateIK4
   
   FUNCTION hunt(xx,x,jlo) RESULT( locate )
     IMPLICIT NONE
     INTEGER(KIND=4), INTENT(INOUT) :: jlo
     REAL(KIND=4), INTENT(IN) :: x
     REAL(KIND=4), DIMENSION(:), INTENT(IN) :: xx
     INTEGER(KIND=4) :: locate
     !! Given an array xx(1:N), and given a value x, returns a value jlo 
     !! such that x is between xx(jlo) and xx(jlo+1). xx must be monotonic, 
     !! either increasing or decreasing. jlo = 0 or jlo = N is returned 
     !! to indicate that x is out of range. jlo on input is taken as the
     !! initial guess for jlo on output.
     INTEGER(KIND=4) :: n,inc,jhi,jm
     LOGICAL :: ascnd
     n = SIZE(xx)
     ascnd = (xx(n) >= xx(1))  !! True if ascending order of table, 
                               !! false otherwise.
     IF( jlo <= 0 .OR. jlo > n ) THEN !! Input guess not useful. 
        jlo = 0                       !! Go immediately to bisection.
        jhi = n+1
     ELSE
        inc = 1                              !! Set the hunting increment.
        IF( x >= xx(jlo) .EQV. ascnd ) THEN  !! Hunt up:
          DO
             jhi = jlo+inc
             IF( jhi > n ) THEN              !! Done hunting, since off 
                jhi = n+1                    !! end of table.
                EXIT
             ELSE
                IF( x < xx(jhi) .EQV. ascnd ) EXIT
                jlo = jhi                    !! Not done hunting,
                inc = inc+inc                 ! so double the increment
             ENDIF
          ENDDO                              !! and try again.
        ELSE                                 !! Hunt down:
          jhi = jlo
          DO
            jlo = jhi-inc
            IF( jlo < 1 ) THEN               !! Done hunting, since off 
              jlo = 0                         ! end of table.
              EXIT
            ELSE
              IF( x >= xx(jlo) .EQV. ascnd ) EXIT
              jhi = jlo                      !! Not done hunting,
              inc = inc+inc                  !! so double the increment
            ENDIF
          ENDDO                              !! and try again.
        ENDIF
     ENDIF                                   !! Done hunting, value bracketed.
     DO                                      !!Hunt is done, so begin the
        IF( jhi-jlo <= 1 ) THEN               !final bisection phase:
          IF( abs( x-xx(n) ) < 1.0D-4 ) jlo = n-1
          IF( abs( x-xx(1) ) < 1.0D-4 ) jlo = 1
          EXIT
        ELSE
          jm = (jhi+jlo)/2
          IF( x >= xx(jm) .EQV. ascnd ) THEN
            jlo = jm
          ELSE
            jhi = jm
          ENDIF
        ENDIF
     ENDDO
     locate = jlo
   END FUNCTION hunt

!   FUNCTION LintextRK8p( xa, ya, x ) RESULT( fun )
!     IMPLICIT NONE
!     REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: xa, ya
!     REAL (KIND = 8), INTENT(IN) :: x
!     REAL (KIND = 8) :: fun
!     REAL (KIND = 8) :: h, fra
!     INTEGER (KIND = 4) :: klo, khi, n
!     n = SIZE( xa ) 
!     IF( n .NE. SIZE( ya ) ) THEN
!        WRITE(*, '( "inconsistent input data points in Lintext")' )
!        STOP
!     ENDIF
!  
!     klo = MAX( MIN( locate( xa, x ), n-1 ), 1 )
!     khi = klo + 1
!     h   = xa(khi) - xa(klo)
!     IF( h == 0.0D0 ) THEN
!        WRITE(*,*) "bad xa input in Lintext"
!        STOP
!     ENDIF
!
!     fra = (ya(khi) - ya(klo))/h
!     fun = ya(klo) + (x - xa(klo))*fra
!   END FUNCTION LintextRK8p
!
!   FUNCTION LintextRK8a( xa, ya, x ) RESULT( fun )
!     IMPLICIT NONE
!     REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: xa, ya
!     REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: x
!     REAL (KIND = 8), DIMENSION(SIZE(x)) :: fun
!     REAL (KIND = 8) :: h, fra
!     INTEGER (KIND = 4) :: klo, khi, n, ii
!     n = SIZE( xa ) 
!     IF( n .NE. SIZE( ya ) ) THEN
!        WRITE(*, '( "inconsistent input data points in Lintext")' )
!        STOP
!     ENDIF
!
!     DO ii = 1, SIZE(x)
!       klo = MAX( MIN( locate( xa, x(ii) ), n-1 ), 1 )
!       khi = klo + 1
!       h   = xa(khi) - xa(klo)
!       IF( h == 0.0D0 ) THEN
!          WRITE(*,*) "bad xa input in Lintext"
!          STOP
!       ENDIF
! 
!       fra     = (ya(khi) - ya(klo))/h
!       fun(ii) = ya(klo) + (x(ii) - xa(klo))*fra
!     ENDDO
!   END FUNCTION LintextRK8a
!
!   FUNCTION LintextRK4a( xa, ya, x ) RESULT( fun )
!     IMPLICIT NONE
!     REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: xa, ya
!     REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: x
!     REAL (KIND = 4), DIMENSION(SIZE(x)) :: fun
!     REAL (KIND = 4) :: h, fra
!     INTEGER (KIND = 4) :: klo, khi, n, ii
!     n = SIZE( xa )
!     IF( n .NE. SIZE( ya ) ) THEN
!        WRITE(*, '( "inconsistent input data points in Lintext")' )
!        STOP
!     ENDIF
!   
!     DO ii = 1, SIZE(x)
!       klo = MAX( MIN( locate( xa, x(ii) ), n-1 ), 1 )
!       khi = klo + 1
!       h   = xa(khi) - xa(klo)
!       IF( h == 0.0D0 ) THEN
!          WRITE(*,*) "bad xa input in Lintext"
!          STOP
!       ENDIF
! 
!       fra     = (ya(khi) - ya(klo))/h
!       fun(ii) = ya(klo) + (x(ii) - xa(klo))*fra
!     ENDDO
!   END FUNCTION LintextRK4a
!
!   FUNCTION LintextRK4p( xa, ya, x ) RESULT( fun )
!     IMPLICIT NONE
!     REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: xa, ya
!     REAL (KIND = 4), INTENT(IN) :: x
!     REAL (KIND = 4) :: fun
!     REAL (KIND = 4) :: h, fra
!     INTEGER (KIND = 4) :: klo, khi, n
!     n = SIZE( xa )
!     IF( n .NE. SIZE( ya ) ) THEN
!        WRITE(*, '( "inconsistent input data points in Lintext")' )
!        STOP
!     ENDIF
! 
!     klo = MAX( MIN( locate( xa, x ), n-1 ), 1 )
!     khi = klo + 1
!     h   = xa(khi) - xa(klo)
!     IF( h == 0.0D0 ) THEN
!        WRITE(*,*) "bad xa input in Lintext"
!        STOP
!     ENDIF
!
!     fra = (ya(khi) - ya(klo))/h
!     fun = ya(klo) + (x - xa(klo))*fra
!   END FUNCTION LintextRK4p

!   FUNCTION slit( x_center, x_fwhm, x, type )
!     IMPLICIT NONE
!     REAL (KIND = 8), INTENT(IN) :: x_center, x_fwhm, x
!     INTEGER (KIND = 4), INTENT(IN), OPTIONAL :: type
!     INTEGER (KIND = 4) :: ltype
!     REAL (KIND = 8) :: x_dist
!     REAL (KIND = 8) :: slit
!  
!     IF( .not. PRESENT( type ) ) THEN
!        ltype = 1
!     ELSE
!        ltype = type
!     ENDIF
!  
!     x_dist = ABS( x - X_center )
!     IF( ltype == 1 ) THEN            ! triangular
!       IF( x_dist >= x_fwhm ) THEN 
!          slit = 0.0D0
!       ELSE
!          slit = (1.0D0 - x_dist / x_fwhm)/x_fwhm
!       ENDIF
!     ELSE IF( ltype == 0 )  THEN      ! square
!       IF( x_dist <= 0.5D0*x_fwhm ) THEN
!          slit = 1.0D0/x_fwhm
!       ELSE
!          slit = 0.0D0
!       ENDIF
!     ELSE
!       WRITE(*,'("unknown filter type =", I3, " in FUNCTION slit" )' ) ltype
!       STOP
!     ENDIF
!   END FUNCTION slit

!   FUNCTION Itrapez( xa, ya, x_center, x_fwhm, type ) RESULT (Isum)
!     REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: xa, ya
!     REAL (KIND = 8), DIMENSION(:), ALLOCATABLE :: xla, yla, wt
!     REAL (KIND = 8), INTENT(IN) :: x_center, x_fwhm
!     INTEGER (KIND = 4), INTENT(IN), OPTIONAL :: type
!     INTEGER (KIND = 4) :: ltype
!     REAL (KIND = 8) :: y_center, xmin, ymin, xmax, ymax, Isum
!     INTEGER (KIND = 4) :: i, imin, imax, ictr, n
!     INTEGER (KIND = 4) :: ierr
!
!     IF( .not. PRESENT( type ) ) THEN
!        ltype = 1
!     ELSE
!        ltype = type
!     ENDIF
!  
!     IF( x_fwhm < EPSILON(1.0D0) ) THEN
!        WRITE(*,'("Wrong input: x_fwhm = ", F12.6 )' ) x_fwhm
!        STOP
!     ENDIF
!
!     IF( ltype == 1 ) THEN    !triangular
!       xmin = x_center - x_fwhm
!       xmax = x_center + x_fwhm
!     ELSE IF( ltype == 0 ) THEN !square
!       xmin = x_center - (0.5D0-EPSILON(0.5D0))*x_fwhm  
!       xmax = x_center + (0.5D0-EPSILON(0.5D0))*x_fwhm 
!     ELSE
!       WRITE(*,'("unknown filter type =", I3, " in FUNCTION Itrapez" )' ) ltype
!       STOP
!     ENDIF
!     
!     imin = locate( xa, xmin )
!     imax = locate( xa, xmax )
!     ictr = locate( xa, x_center )
!
!     n = imax - imin
!     n = n + 3
!    
!     ALLOCATE( xla(n), yla(n), wt(n), STAT = ierr )
!     IF( ierr .NE. 0 ) THEN
!       WRITE(*, '("Allocate xla, yla, wt denied in Itrapez" )' )
!       STOP
!     ENDIF
!
!     ymin     = Lintext( xa, ya, xmin )
!     ymax     = Lintext( xa, ya, xmax )
!     y_center = Lintext( xa, ya, x_center ) 
!
!     xla = (/ xmin, xa(imin+1:ictr), x_center, xa(ictr+1:imax), xmax /)
!     yla = (/ ymin, ya(imin+1:ictr), y_center, ya(ictr+1:imax), ymax /)
!     wt  = (/ (slit( x_center, x_fwhm, xla(i), ltype ), i = 1,n ) /)
!     Isum = SUM( (/ (0.5D0*(xla(i+1)-xla(i))*(wt(i+1)*yla(i+1)+wt(i)*yla(i)), &
!                    i=1,n-1 ) /) )
!     DEALLOCATE( xla, yla, wt, STAT = ierr )
!     IF( ierr .NE. 0 ) THEN
!       WRITE(*, '("Deallocate xla, yla, wt denied in Itrapez" )' )
!       STOP
!     ENDIF
!
!   END FUNCTION Itrapez

   FUNCTION adjustDEG( phi ) RESULT( adjphi )
     REAL( KIND = 4 ), INTENT(IN) :: phi
     REAL( KIND = 4 ) :: adjphi

     IF( phi > 180 ) THEN
        adjphi = phi - 360.0
     ELSE IF( phi < -180 ) THEN
        adjphi = phi + 360.0
     ELSE
        adjphi = phi 
     ENDIF
     RETURN
   END FUNCTION adjustDEG 

   FUNCTION toLowerString( strIN ) RESULT( str_low )
     CHARACTER(LEN=*), INTENT(IN) :: strIN
     CHARACTER(LEN=LEN(strIN)) :: str_low
    
     INTEGER (KIND=4) :: i, letter
     
     DO i = 1, LEN_TRIM( strIN )
       letter = INDEX( ALPHABETu, strIN(i:i) ) 
       IF( letter > 0 ) THEN
         str_low(i:i) = alphabetL(letter:letter)
       ELSE
         str_low(i:i) = strIN(i:i)
       ENDIF
     ENDDO
     RETURN
   END FUNCTION toLowerString

!   FUNCTION toUpperString( strIN ) RESULT( str_up )
!     CHARACTER(LEN=*), INTENT(IN) :: strIN
!     CHARACTER(LEN=LEN(strIN)) :: str_up
!    
!     INTEGER (KIND=4) :: i, letter
!     
!     DO i = 1, LEN_TRIM( strIN )
!       letter = INDEX( alphabetL, strIN(i:i) ) 
!       IF( letter > 0 ) THEN
!         str_up(i:i) = ALPHABETu(letter:letter)
!       ELSE
!         str_up(i:i) = strIN(i:i)
!       ENDIF
!     ENDDO
!     RETURN
!   END FUNCTION toUpperString

   FUNCTION deQuote( StringValue ) RESULT( nc )
     CHARACTER(LEN=*), INTENT(INOUT) :: StringValue
     CHARACTER(LEN=LEN(StringValue)) :: strTemp
     INTEGER :: nc, di

     di = INDEX( StringValue, '"' )
     IF( di > 0 ) THEN
        strTemp = StringValue(di+1:)
        di = INDEX( strTemp, '"', BACK = .TRUE. )
        IF( di > 0 ) THEN
           StringValue = strTemp(1:di-1)
        ELSE
               StringValue = strTemp(1:)
        ENDIF
     ENDIF
     nc = LEN( TRIM(StringValue) )
     RETURN
   END FUNCTION deQuote

END MODULE UTIL_tools_class
