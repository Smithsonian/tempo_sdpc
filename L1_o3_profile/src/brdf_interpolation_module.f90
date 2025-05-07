MODULE interpolation_module
  
  USE parameters_module
  USE level1_def,       ONLY : GeolocationType
  USE error_module,     ONLY : ErrorType, CheckError
  
  IMPLICIT NONE
  
  PUBLIC :: BSPLINE, &
            SPLINE1, &
            SPLINT1, &
            NearestNeighbourSampler


  TYPE XYGridType
    INTEGER                   :: imx
    INTEGER                   :: jmx
    INTEGER                   :: tmx
    REAL(KIND=8), ALLOCATABLE :: LongitudeMid(:)
    REAL(KIND=8), ALLOCATABLE :: LongitudeEdge(:)
    REAL(KIND=8), ALLOCATABLE :: LatitudeMid(:)
    REAL(KIND=8), ALLOCATABLE :: LatitudeEdge(:)
    REAL(KIND=8), ALLOCATABLE :: TauTime(:)
  ENDTYPE XYGridType

  TYPE XYGridWtType
    INTEGER                   :: SampleTypeIndex
    INTEGER                   :: imin
    INTEGER                   :: imax
    INTEGER                   :: jmin
    INTEGER                   :: jmax
    INTEGER                   :: tmin
    INTEGER                   :: tmax
    INTEGER                   :: ni
    INTEGER                   :: nj
    INTEGER                   :: nt
    REAL(KIND=8), ALLOCATABLE :: Weight(:,:,:)
    LOGICAL                   :: PeriodicOverlap ! Flag set if pixel overlaps +/- 180
    INTEGER                   :: imin2
    INTEGER                   :: imax2
    INTEGER                   :: ni2
    REAL(KIND=8), ALLOCATABLE :: Weight2(:,:,:)

  ENDTYPE XYGridWtType

  CONTAINS
  
  SUBROUTINE PIECELIN_INTERP(X,Y,U,V,NX,NU,CallingSubroutine)
    
    ! Piecewise linear interpolation
    ! Assumes x is already sorted in ascending order
    
    ! Input arguments
    INTEGER,       INTENT(IN)                :: NX, NU
    REAL(KIND=8),  INTENT(IN), DIMENSION(NX) :: X,Y
    REAL(KIND=8),  INTENT(IN), DIMENSION(0:NU) :: U
    
    ! Output argument
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NU) :: V
    
    ! Optional Arguments
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CallingSubroutine

    ! Local variables
    REAL(KIND=8), DIMENSION(NX-1) :: delta
    INTEGER,      DIMENSION(0:NU)   :: k
    INTEGER                       :: i, j
    REAL(KIND=8)                  :: s
    
    
    ! =================================================
    ! SUBROUTINE PIECELIN_INTERP
    ! =================================================
    
    ! Compute gradient
    DO i=1,nx-1
      delta(i) = (y(i+1)-y(i))/(x(i+1)-x(i))
    ENDDO
    
    ! Find subinterval indices
    k(:) = 1
    
    ! Get subintervals
    DO i=2,nx-1
      
      DO j=0,nu
        IF( x(i) .LE. u(j) ) k(j) = i
      ENDDO
      
    ENDDO
    
    ! Evaluate interpolant
    DO j=0,nu
      s = u(j) - x(k(j))
      v(j) = y(k(j)) + s*delta(k(j))
    ENDDO
    
    RETURN
    
  END SUBROUTINE PIECELIN_INTERP
  
  INTEGER FUNCTION RegularGridNearestIndex(Wvl0,dWvl,nWvl,Wvl,FlagOutsideRange)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL(KIND=8),      INTENT(IN) :: Wvl0
    REAL(KIND=8),      INTENT(IN) :: dWvl
    INTEGER,           INTENT(IN) :: nWvl
    REAL(KIND=8),      INTENT(IN) :: Wvl
    LOGICAL, OPTIONAL, INTENT(IN) :: FlagOutsideRange

    ! ---------------
    ! local variables
    ! ---------------

    ! =====================================================================
    ! RegularGridNearestIndex Starts here
    ! =====================================================================

    
    RegularGridNearestIndex = NINT((Wvl-Wvl0)/dWvl)+1

    ! Set to boundary if out of range
    IF(RegularGridNearestIndex .LT. 1) THEN
      RegularGridNearestIndex = 1
      IF(PRESENT(FlagOutsideRange)) THEN
        IF(FlagOutsideRange) RegularGridNearestIndex = -1
      ENDIF
    ELSEIF(RegularGridNearestIndex .GT. nWvl) THEN
      RegularGridNearestIndex = nWvl
      IF(PRESENT(FlagOutsideRange)) THEN
        IF(FlagOutsideRange) RegularGridNearestIndex = -1
      ENDIF
    ENDIF

  END FUNCTION RegularGridNearestIndex



  ! This subroutine combines spline and splint function
  ! in Numerical Recipes by Press et al., 1997.
  SUBROUTINE BSPLINE(xa, ya, n, x, y, np, errstat)

    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    
    INTEGER, INTENT (IN)                      :: n, np
    INTEGER, INTENT (OUT)                     :: errstat
    REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
    REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y
    
    REAL (KIND=dp), DIMENSION(n)              :: y2a, xacpy
    REAL (KIND=dp), DIMENSION(np)             :: xcpy
    REAL (KIND=dp), DIMENSION(n-1)            :: diff
    REAL (KIND=dp)                            :: xmin, xmax, xamin, xamax
    
    errstat = 0
    IF (n < 3) THEN
      errstat = - 1; RETURN
    ENDIF
    
    diff = xa(2:n) - xa(1:n-1)
    IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
      errstat =  -2; RETURN
    ENDIF

    xmin = MINVAL(x); xmax = MAXVAL(x)
    xamin = MINVAL(xa); xamax = MAXVAL(xa)
    IF (xmin < xamin .OR. xmax > xamax) THEN
      errstat =  -3; RETURN
    ENDIF
    
    IF (xa(1) < xa(n)) THEN
      CALL SPLINE1(xa, ya, n, y2a)
      CALL SPLINT1(xa, ya, y2a, n, x, y, np)
    ELSE
      xacpy = -xa; xcpy = -x
      CALL SPLINE1(xacpy, ya, n, y2a)
      CALL SPLINT1(xacpy, ya, y2a, n, xcpy, y, np)
    ENDIF

    RETURN
  END SUBROUTINE BSPLINE
  
  ! Not efficient - use for small things
  SUBROUTINE LinearInt_Edgefill(xa, ya, n, x_in, y, np,errstat)
  
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,      INTENT(IN)  :: n, np
    REAL(KIND=8), INTENT(IN)  :: xa(n), ya(n)
    REAL(KIND=8), INTENT(IN)  :: x_in(np)
    REAL(KIND=8), INTENT(OUT) :: y(np)
    INTEGER,      INTENT(OUT) :: errstat
    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER      :: idx0, i
    REAL(KIND=8) :: xwt0,xwt1

    ! =====================================================================
    ! LinearInt_Edgefill starts here
    ! =====================================================================

    errstat = 0

    ! Compute interpolation weights
    DO i=1,np
      
      ! Compute Linear interpolation weights
      CALL LinearInterpolationWeights(n,xa,x_in(n),idx0,xwt0,xwt1)

      ! Compute output value
      y(i) = ya(idx0)*xwt0 + ya(idx0+1)*xwt1

    ENDDO

  END SUBROUTINE LinearInt_Edgefill


  ! This subroutine combines spline and splint function
  ! in Numerical Recipes by Press et al., 1997.
  ! Points outside range will be set to the value at end of interval
  SUBROUTINE BSPLINE_EdgeFill(xa, ya, n, x_in, y, np, errstat,oob_low,oob_hi,fill_value)

    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    
    INTEGER, INTENT (IN)                       :: n, np
    INTEGER, INTENT (OUT)                      :: errstat
    REAL (KIND=dp), DIMENSION(n),  INTENT(IN)  :: xa, ya
    REAL (KIND=dp), DIMENSION(np), INTENT(IN)  :: x_in
    REAL (KIND=dp), DIMENSION(np), INTENT(OUT) :: y
    LOGICAL,        OPTIONAL,      INTENT(OUT) :: oob_low, oob_hi
    REAL(KIND=8),   OPTIONAL,      INTENT(IN)  :: fill_value


    REAL (KIND=dp), DIMENSION(n)               :: y2a, xacpy
    REAL (KIND=dp), DIMENSION(np)              :: xcpy, x
    REAL (KIND=dp), DIMENSION(n-1)             :: diff
    REAL (KIND=dp)                             :: xmin, xmax, xamin, xamax
    INTEGER                                    :: i
    LOGICAL                                    :: oob_low_tmp, oob_hi_tmp

    errstat = 0
    IF (n < 3) THEN
      errstat = - 1; RETURN
    ENDIF
    
    ! Zero output
    y(:) = 0.0d0

    ! Initialize Out of bounds flags
    oob_low_tmp =.FALSE. ; oob_hi_tmp = .FALSE.

    ! Initialize x
    x = x_in

    diff = xa(2:n) - xa(1:n-1)
    IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
      errstat =  -2; RETURN
    ENDIF
    
    xamin = MINVAL(xa); xamax = MAXVAL(xa)

    ! Fill values outside range
    DO i=1,np
      IF(x(i) .LT. xamin) THEN
        x(i) = xamin+TINY(0.0d0) ; oob_low_tmp = .TRUE.
      ELSEIF(x(i) .GT. xamax) THEN
        x(i) = xamax-TINY(0.0d0) ; oob_hi_tmp = .TRUE.
      ENDIF
    ENDDO

    IF (xa(1) < xa(n)) THEN
      CALL SPLINE1(xa, ya, n, y2a)
      CALL SPLINT1(xa, ya, y2a, n, x, y, np)
    ELSE
      xacpy = -xa; xcpy = -x
      CALL SPLINE1(xacpy, ya, n, y2a)
      CALL SPLINT1(xacpy, ya, y2a, n, xcpy, y, np)
    ENDIF

    IF(PRESENT(fill_value)) THEN
      DO i=1,np
        IF(x_in(i) .LT. xamin .OR. x_in(i) .GT. xamax) THEN
          y(i) = fill_value
        ENDIF
      ENDDO
    ENDIF

    IF(PRESENT(oob_low)) oob_low = oob_low_tmp
    IF(PRESENT(oob_hi )) oob_hi  = oob_hi_tmp

    RETURN

  END SUBROUTINE BSPLINE_EdgeFill

  ! modified to always use "natural" boundary conditions
  SUBROUTINE SPLINE1 (x, y, n, y2, CallingSubroutine)
    
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: x, y
    REAL (KIND=dp), DIMENSION(n), INTENT(OUT) :: y2
    
    ! Optional Arguments
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CallingSubroutine

    REAL (KIND=dp), DIMENSION(n) :: u
    INTEGER       :: i, k
    REAL(KIND=dp) :: sig, p, qn, un
    
    y2 (1) = 0.0
    u (1) = 0.0
    
    DO i = 2, n - 1
      sig = (x (i) - x (i - 1)) / (x (i + 1) -x (i - 1))
      p = sig * y2 (i - 1) + 2.D0
      y2 (i) = (sig - 1.) / p
      u (i) = (6._dp * ((y (i + 1) - y (i)) / (x (i + 1) - x (i)) -  & 
            (y (i) - y (i - 1)) / (x (i) - x (i - 1))) / (x (i + 1) - &
            x (i - 1)) - sig * u (i - 1)) / p
    ENDDO
    
    qn = 0.0
    un = 0.0
    y2 (n) = (un - qn * u (n - 1)) / (qn * y2 (n - 1) + 1.0)
    DO k = n - 1, 1, -1
      y2 (k) = y2 (k) * y2 (k + 1) + u (k)
    ENDDO
    
    RETURN
    
  END SUBROUTINE SPLINE1

  ! This code could be optimized if x is in increasing/decreasing order
  SUBROUTINE SPLINT1 (xa, ya, y2a, n, x, y, m, CallingSubroutine)
    
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    INTEGER, INTENT(IN) :: n, m
    REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: xa, ya, y2a
    REAL (KIND=dp), DIMENSION(m), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(m), INTENT(OUT):: y
    
    ! Optional Arguments
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CallingSubroutine

    INTEGER        :: ii, klo, khi, k 
    REAL (KIND=dp) :: h, a, b

    !klo = 1; khi = n
    DO ii = 1, m
      klo = 1; khi = n
      
      DO WHILE (khi - klo > 1)
          k = (khi + klo) / 2
          IF (xa (k) > x(ii)) THEN
            khi = k
          ELSE
            klo = k
          ENDIF
      ENDDO
      
      h = xa (khi) - xa (klo)
      IF (h == 0.0) STOP 'Bad xa input in: splint!!!'
      a = (xa (khi) - x(ii)) / h
      b = (x(ii) - xa (klo)) / h
      
      y(ii) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + &
            (b**3 - b) * y2a (khi)) * (h**2) / 6.0
    ENDDO
    
    RETURN
  END SUBROUTINE SPLINT1


  SUBROUTINE NearestNeighbourSampler(XYGrid, Geolocation, GridWt, Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(XYGridType),      INTENT(IN)    :: XYGrid
    TYPE(GeolocationType), INTENT(IN)    :: Geolocation
    TYPE(XYGridWtType),    INTENT(INOUT) :: GridWt
    TYPE(ErrorType),       INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local Variables
    ! ---------------
    

    ! =====================================================================
    ! NearestNeighbourSampler starts here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Periodic overlap is not needed in this case
    GridWt%PeriodicOverlap = .FALSE.

    ! Check Allocatation of weighting array
    iF(ALLOCATED(GridWt%Weight)) DEALLOCATE(GridWt%Weight)
    ALLOCATE(GridWt%Weight(1,1,1)) ; GridWt%ni = 1 ; GridWt%nj = 1 ; GridWt%nt = 1
    
    ! Set Grid Weight
    GridWt%Weight(1,1,1) = 1.0d0

    ! Find Nearest longitude index
    GridWt%imin = MINLOC(ABS(XYGrid%LongitudeMid(:)-Geolocation%Longitude),DIM=1) 
    GridWt%imax = GridWt%imin
    
    ! Find nearest latitude index
    GridWt%jmin = MINLOC(ABS(XYGrid%LatitudeMid(:)-Geolocation%Latitude),DIM=1) 
    GridWt%jmax = GridWt%jmin

    ! Find nearest tau time index
    GridWt%tmin = MINLOC(ABS(XYGrid%TauTime(:)-Geolocation%Time%Tau),DIM=1) 
    GridWt%tmax = GridWt%tmin
    
  END SUBROUTINE NearestNeighbourSampler
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: PointInPolygonSampler
  ! 
  ! DESCRIPTION: Computes XY regrid weights - Every point whose centroid 
  !              is within the pixel boundaries is given equal weight
  !              Observations are linearly interpolated in time

  SUBROUTINE PointInPolygonSampler(XYGrid, Geolocation, GridWt, Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(XYGridType),      INTENT(IN)    :: XYGrid
    TYPE(GeolocationType), INTENT(IN)    :: Geolocation
    TYPE(XYGridWtType),    INTENT(INOUT) :: GridWt
    TYPE(ErrorType),       INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local Variables
    ! ---------------
    REAL(KIND=8)              :: elons(2), elats(2), elons2(2)
    REAL(KIND=8)              :: twt0, twt1, npt, npt2
    INTEGER                   :: I,J, II, JJ
    
    ! =====================================================================
    ! PointInPolygonSampler starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Get the longitude/latitude bounds of the pixel
    elons(1) = minval(Geolocation%CornerLongitudes)
    elons(2) = maxval(Geolocation%CornerLongitudes)
    elats(1) = minval(Geolocation%CornerLatitudes)
    elats(2) = maxval(Geolocation%CornerLatitudes)

    ! Check for periodic overlap
    GridWt%PeriodicOverlap = .FALSE.
    IF(elons(2) - elons(1) .GT. 180.0d0 ) THEN
      
      ! Set flag
      GridWt%PeriodicOverlap = .TRUE.

      ! Create two new polygons
      elons2(1) = elons(1) ; elons2(2) = 180.0 
      elons(2) = elons(1)  ; elons(1) = -180.0
      
    ENDIF

    ! Get Indices

    ! Longitude
    GridWt%imin = MINVAL( MINLOC( ABS(XYGrid%LongitudeMid-elons(1)) ) ) - 1
    IF(GridWt%imin .LE. 0) GridWt%imin = 1
    GridWt%imax = MAXVAL( MINLOC( ABS(XYGrid%LongitudeMid-elons(2)) ) ) + 1
    IF(GridWt%imax .GT. XYGrid%imx) GridWt%imax = XYGrid%imx
    GridWt%ni = GridWt%imax - GridWt%imin + 1

    ! Latitude
    GridWt%jmin = MINVAL( MINLOC( ABS(XYGrid%LatitudeMid-elats(1)) ) ) - 1
    IF(GridWt%jmin .LE. 0) GridWt%jmin = 1
    GridWt%jmax = MAXVAL( MINLOC( ABS(XYGrid%LatitudeMid-elats(2)) ) ) + 1
    IF(GridWt%jmax .GT. XYGrid%jmx) GridWt%jmax = XYGrid%jmx
    GridWt%nj = GridWt%jmax - GridWt%jmin + 1

    ! Second set of longitudes if required
    IF(GridWt%PeriodicOverlap) THEN
      GridWt%imin2 = MINVAL( MINLOC( ABS(XYGrid%LongitudeMid-elons2(1)) ) ) - 1
      IF(GridWt%imin2 .LE. 0) GridWt%imin2 = 1
      GridWt%imax = MAXVAL( MINLOC( ABS(XYGrid%LongitudeMid-elons2(2)) ) ) + 1
      IF(GridWt%imax2 .GT. XYGrid%imx) GridWt%imax2 = XYGrid%imx
      GridWt%ni2 = GridWt%imax2 - GridWt%imin2 + 1
    ENDIF

    ! Get indices and time-weighted fractions
    IF(Geolocation%Time%Tau .GE. MAXVAL(XYGrid%TauTime) ) THEN
      GridWt%tmin = XYGrid%tmx ; GridWt%nt = 1 ; twt0 = 1.0d0 ; twt1 = 0.0d0
    ELSE
      GridWt%tmin = MINVAL(MAXLOC(XYGrid%TauTime, MASK=(XYGrid%TauTime .LT. Geolocation%Time%Tau)))
      GridWt%tmax = GridWt%tmin + 1
      twt0 = (Geolocation%Time%Tau-XYGrid%TauTime(GridWt%tmin)) &
           / (XYGrid%TauTime(GridWt%tmax)-XYGrid%TauTime(GridWt%tmin))
      !twt0 = (Geolocation%Time%Tau-XYGrid%TauTime(idt0))/(XYGrid%TauTime(idt1)-XYGrid%TauTime(idt0))
      twt1 = 1.0d0 - twt0
      GridWt%nt = 2
    ENDIF

    ! Allocate Grid weighting arrays
    iF(ALLOCATED(GridWt%Weight)) DEALLOCATE(GridWt%Weight)
    iF(ALLOCATED(GridWt%Weight2)) DEALLOCATE(GridWt%Weight2)
    ALLOCATE(GridWt%Weight(GridWt%ni,GridWt%nj,GridWt%nt))
    IF(GridWt%PeriodicOverlap) ALLOCATE(GridWt%Weight2(GridWt%ni2,GridWt%nj,GridWt%nt))

    ! Zero arrays
    GridWt%Weight = 0.0d0 ; IF(GridWt%PeriodicOverlap) GridWt%Weight2 = 0.0d0
    
    ! Initialize Point counters
    npt =  0.0d0 ; npt2 = 0.0d0

    ! Do the PNPOLY Tests for the first subpixel
    DO J=1,GridWt%nj
      JJ = GridWt%jmin + J - 1
      DO I=1,GridWt%ni
        II = GridWt%imin + I - 1

        IF( PNPOLY_r8(4,Geolocation%CornerLongitudes,&
                        Geolocation%CornerLatitudes, &
                        XYGrid%LongitudeMid(II),     &
                        XYGrid%LatitudeMid(JJ) )     ) THEN

          IF(GridWt%nt .EQ. 1) THEN
            GridWt%Weight(I,J,1) = 1.0d0
          ELSE
            GridWt%Weight(I,J,1) = twt0
            GridWt%Weight(I,J,2) = twt1
          ENDIF

          npt = npt + 1.0d0

        ENDIF
      ENDDO
    ENDDO

    ! Do the PNPOLY Tests for the second if required
    IF(GridWt%PeriodicOverlap) THEN
      
      DO J=1,GridWt%nj
        JJ = GridWt%jmin + J - 1
        DO I=1,GridWt%ni2
          II = GridWt%imin2 + I - 1

          IF( PNPOLY_r8(4,Geolocation%CornerLongitudes,&
                          Geolocation%CornerLatitudes, &
                          XYGrid%LongitudeMid(II),     &
                          XYGrid%LatitudeMid(JJ))      ) THEN

            IF(GridWt%nt .EQ. 1) THEN
              GridWt%Weight(I,J,1) = 1.0d0
            ELSE
              GridWt%Weight(I,J,1) = twt0
              GridWt%Weight(I,J,2) = twt1
            ENDIF

            npt2 = npt2 + 1.0d0

          ENDIF
        ENDDO
      ENDDO
    
    ENDIF

    ! If there were no points inside the polygon do nearest neighbour sampling
    IF(npt .LT. TINY(0.0d0) .AND. npt2 .LT. TINY(0.0d0)) THEN
      CALL NearestNeighbourSampler(XYGrid, Geolocation, GridWt, Error)
    
    ELSE

      ! Renormalize 
      IF(npt .GT. TINY(0.0D0)) THEN
        GridWt%Weight = GridWt%Weight/(npt+npt2)
      ENDIF
      IF(npt2 .GT. TINY(0.0d0)) THEN
        GridWt%Weight2 = GridWt%Weight2 / (npt+npt2)
      ENDIF

    ENDIF
    
  END SUBROUTINE PointInPolygonSampler
  
  LOGICAL FUNCTION PNPOLY_r8( nvert, vertx, verty, testx, testy)
    
    IMPLICIT NONE
    
    ! Return variable
    ! ---------------
    !LOGICAL :: PNPOLY ! True if point is in polygon
    
    ! Input
    INTEGER :: nvert ! # Vertices in polygon
    REAL(KIND=8), DIMENSION(nvert) :: vertx ! X coordinates of the polygon vertices
    REAL(KIND=8), DIMENSION(nvert) :: verty ! Y coordinates of the polygon vertices
    REAL(KIND=8)                   :: testx ! X coordinate of test point
    REAL(KIND=8)                   :: testy ! Y coordinate of test point
    
    ! Local variables
    INTEGER :: i, j
    
    ! =====================================================================
    ! PNPOLY_r8 starts here
    ! =====================================================================

    ! Initialize PNPOLY
    PNPOLY_r8 = .FALSE.
    
    ! Initialize j
    j = nvert
    
    ! Perform crossings test (Jordan curve thm.) to see if point is in polygon
    DO i=1,nvert
      
      IF ( ((verty(i) .gt. testy) .neqv. (verty(j) .gt. testy)) .and. &
          (testx .lt. (vertx(j)-vertx(i)) * &
          (testy-verty(i)) / (verty(j)-verty(i)) + vertx(i)) ) THEN
         
         PNPOLY_r8 = .NOT. PNPOLY_r8
         
      ENDIF
      
      j = i
      
    ENDDO
    
    RETURN
    
  END FUNCTION PNPOLY_r8

  SUBROUTINE LinearInterpolationWeights(nx,x,x_int,idx0,xwt0,xwt1)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,      INTENT(IN)  :: nx
    REAL(KIND=8), INTENT(IN)  :: x(nx)
    REAL(KIND=8), INTENT(IN)  :: x_int
    INTEGER,      INTENT(OUT) :: idx0
    REAL(KIND=8), INTENT(OUT) :: xwt0
    REAL(KIND=8), INTENT(OUT) :: xwt1

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: i
    
    ! =====================================================================
    ! LinearInterpolationWeights starts here
    ! =====================================================================
    
    IF(x_int .LE. x(1)) THEN
      idx0 = 1
      xwt0 = 1.0d0
      xwt1 = 0.0d0
    ELSEIF(x_int .GE. x(nx)) THEN
      idx0 = nx-1
      xwt0 = 0.0d0
      xwt1 = 1.0d0
    ELSE
      idx0 = -1
      i    =  1
      DO WHILE(idx0 .LT. 0)
        ! IF(i .GT. nx) THEN
        !   print*,x
        !   print*,'--->',x_int
        ! ENDIF
        IF(x(i) .GT. x_int) THEN
          idx0 = i-1
          xwt0 = (x(idx0+1)-x_int)/(x(idx0+1)-x(idx0))
          xwt1 = 1.0d0 - xwt0
        ENDIF
        i = i+1
      ENDDO

    ENDIF

  END SUBROUTINE LinearInterpolationWeights
  
  SUBROUTINE LinearInterpolationWeights_InterpArray(nx,x,n_int,x_int,idx0,xwt0,xwt1)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,      INTENT(IN)  :: nx
    REAL(KIND=8), INTENT(IN)  :: x(nx)
    INTEGER,      INTENT(IN)  :: n_int
    REAL(KIND=8), INTENT(IN)  :: x_int(n_int) ! Must be sorted in ascending order
    INTEGER,      INTENT(OUT) :: idx0(n_int)
    REAL(KIND=8), INTENT(OUT) :: xwt0(n_int)
    REAL(KIND=8), OPTIONAL, INTENT(OUT) :: xwt1(n_int)

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: i, n
    LOGICAL :: below_min, in_range
    ! =====================================================================
    ! LinearInterpolationWeights starts here
    ! =====================================================================
    
    ! Initialize Output
    idx0(:) = -1 ; xwt0(:) = 0.0d0

    IF(n_int .LE. 0) RETURN
    
    ! First set weights for interpolation values below limit
    n = 1 ; below_min = .TRUE.
    DO WHILE(below_min)
      
      IF(x_int(n) .LE. x(1)) THEN

        ! Update weights for points below
        idx0(n) = 1
        xwt0(n) = 1.0d0
        n = n+1

      ELSE

        ! Set flag to exit
        below_min = .FALSE.

      ENDIF

      ! Also Check if we have finished interpolating
      IF(n .GT. n_int) below_min = .FALSE.

    ENDDO
    
    ! Values in range
    IF(idx0(n_int) .LT. 0) THEN

      i = 1 ; in_range = .TRUE.
      DO WHILE(in_range)
        
        ! Find if we have reached the first index > x_int
        IF(x(i) .GT. x_int(n)) THEN

          ! Update interpolation for point i
          idx0(n) = i-1
          xwt0(n) = (x(idx0(n)+1)-x_int(n))/(x(idx0(n)+1)-x(idx0(n)))

          ! Move on to next search point
          n = n+1

        ELSE
          
          ! Increment Data array point
          i = i+1

        ENDIF
        
        ! Check for exit
        IF(i .GT. nx .OR. n .GT. n_int) in_range = .FALSE. 

      ENDDO
    
    ENDIF
    
    ! The remaining points are above the maximum limit
    IF(idx0(n_int) .LT. 0) THEN
      idx0(n:n_int) = nx-1
      xwt0(n:n_int) = 0.0d0
    ENDIF

    IF(PRESENT(xwt1)) xwt1 = 1.0d0 - xwt0

  END SUBROUTINE LinearInterpolationWeights_InterpArray

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: VerticalRegridWeights
  ! 
  ! DESCRIPTION: Computes vertical regrid weights between two pressure 
  !              grids. 
  !              Code adapted from the GAMAP IDL regrid_column procedure
  
  SUBROUTINE VerticalRegridWeights( n_in, pedge_in, n_out, pedge_out, weight )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,       INTENT(IN)  :: n_in
    REAL(KIND=8),  INTENT(IN)  :: pedge_in(n_in+1)
    INTEGER,       INTENT(IN)  :: n_out
    REAL(KIND=8),  INTENT(IN)  :: pedge_out(n_out+1)
    REAL(KIND=8),  INTENT(OUT) :: weight(n_in,n_out)
    
    ! ---------------
    ! local variables
    ! ---------------
    LOGICAL :: valid
    INTEGER :: l,ll, k, kk
    
    REAL(KIND=8) :: TotWt
    
    ! =====================================================================
    ! VerticalRegridWeights starts here
    ! =====================================================================
    
    ! Zero weights
    weight(:,:) = 0.0d0
    
    DO l=1,n_in
      ll = n_in+2-l
    ENDDO
    
    DO l=1,n_in
      
      ! Reset valid flag
      valid = .FALSE.
      
      ! Indices need to go from BOA -> TOA
      ll = n_in+2-l
      
      ! If the thickness of this pressure level is zero, then this 
      ! means that this pressure level lies below the surface 
      ! pressure (due to topography), as set up in the calling
      ! program.  Therefore, skip to the next INPUT level.
      ! This also helps avoid divide by zero errors. (bmy, 8/6/01)
      IF( (pedge_in(ll) - pedge_in(ll-1)) .GT. 1e-5 ) THEN
        
        ! Loop over output layers
        DO k=1,n_out
          
          kk = n_out+2-k
          
          !==============================================================
          ! No contribution if:
          ! -------------------
          ! Bottom of OUTPUT layer above Top    of INPUT layer  OR
          ! Top    of OUTPUT layer below Bottom of INPUT layer
          !==============================================================
          IF( .NOT. ( pedge_out(kk  ) .LT. pedge_in(ll-1) .OR. &
                      pedge_out(kk-1) .GT. pedge_in(ll)      ) ) THEN
            
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Entire INPUT layer in OUTPUT layer
            !==============================================================
            IF( pedge_out(kk)   .GE. pedge_in(ll) .AND. &
                pedge_out(kk-1) .LE. pedge_in(ll-1)     ) THEN
              
              weight(ll-1,kk-1) = 1.0d0
              valid = .TRUE.
              
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Top of OUTPUT layer in INPUT layer
            !==============================================================
            ELSEIF( pedge_out(kk-1) .LE. pedge_in(ll) .AND. &
                    pedge_out(kk  ) .GE. pedge_in(ll)       ) THEN
              
              weight(ll-1,kk-1) = (pedge_in(ll) - pedge_out(kk-1)) &
                                / (pedge_in(ll) - pedge_in(ll-1) )
              valid = .TRUE.
              
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Entire OUTPUT layer in INPUT layer
            !==============================================================
            ELSEIF( pedge_out(kk  ) .LE. pedge_in(ll  ) .AND. &
                    pedge_out(kk-1) .GE. pedge_in(ll-1)       ) THEN
              
              weight(ll-1,kk-1) = (pedge_out(kk) - pedge_out(kk-1)) &
                                / (pedge_in(ll) - pedge_in(ll-1))
              valid = .TRUE.
              
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Bottom of OUTPUT layer in INPUT layer
            !==============================================================
            ELSEIF(pedge_out(kk  ) .GE. pedge_in(ll-1) .AND. &
                   pedge_out(kk-1) .LE. pedge_in(ll-1)       ) THEN
                   
              weight(ll-1,kk-1) = (pedge_out(kk)-pedge_in(ll-1)) &
                                / (pedge_in(ll) -pedge_in(ll-1))
              valid = .TRUE.
              
            ENDIF
            
          ENDIF
          
        ENDDO
        
      ENDIF
      
      
      
      TotWt = 0.0D0
      DO k=1,n_out
        TotWt = TotWt + weight(ll-1,k)
      ENDDO
      
      
    ENDDO
    
  END SUBROUTINE VerticalRegridWeights
  
  
END MODULE interpolation_module