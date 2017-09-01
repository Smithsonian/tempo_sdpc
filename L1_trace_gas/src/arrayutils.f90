module arrayutils

  private
  public array_smooth, array_sort_r8, array_locate_r4, array_locate_r8, &
    array_find_bounding_indices_r8, array_find_inner_bounding_indices_r8!, &
    !array_find_bounding_indices_r4, array_find_inner_bounding_indices_r4

contains
  subroutine array_smooth (x, n, errstat)

    use OMSAO_precision_module, only: i4, r8
    implicit none
    integer (kind=i4), intent(in) :: n
    real (kind=r8), intent(inout), dimension(n) :: x
    integer, intent(inout) :: errstat
    !
    real (kind=r8), dimension(n) :: tmp

    if (errstat /= 0) return

    ! smoothing with kernel (1/16,1/4,3/8,1/4,1/16)
    if (n < 5) return
    tmp(1:n) = x(1:n)
    x(3:n-2) = 0.375_r8 * tmp (3:n-2) &
      + 0.25_r8 * (tmp (4:n-1) + tmp (2:n-3)) &
      + 0.0625_r8 * (tmp (5:n) + tmp (1:n-4))

  end subroutine array_smooth

  ! This is not a very good general purpose sorting routine.  It
  ! may be fast for nearly sorted arrays, for which it is mainly used.
  ! This is something that needs to be investigated.  --JED
  SUBROUTINE array_sort_r8 (n, x, y )

    USE OMSAO_precision_module, ONLY: r8
    IMPLICIT NONE

    ! --------------
    ! Input variable
    ! --------------
    INTEGER, INTENT (IN) :: n

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (n), INTENT (INOUT) :: x, y

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: xfirst, xtemp, yfirst, ytemp
    INTEGER        :: index, i, j

    ! -----------------------------------
    ! Loop to sort N-1 entries into order
    ! -----------------------------------
    DO i = 1, n-1

      ! ----------------------------------------------------------------
      ! Initialize earliest so far to be in the first place in this pass
      ! ----------------------------------------------------------------
      xfirst = x(i) ; yfirst = y(i)
      index = i

      ! ----------------------------------------------
      ! If X(I+1) > X(I) we can skip to the next index
      ! ----------------------------------------------
      IF ( x(i+1) > x(i) ) CYCLE

      ! --------------------------------------------------
      ! Search remaining (unsorted) items for earliest one
      ! --------------------------------------------------
      DO j = i+1, n
        IF ( x(j) < xfirst ) THEN
          xfirst = x(j) ; yfirst = y(j)
          index = j
        END IF
      END DO

      IF ( index /= i ) THEN
        xtemp    = x(i)     ; ytemp    = y(i)
        x(i)     = x(index) ; y(i)     = y(index)
        x(index) = xtemp    ; y(index) = ytemp
      END IF

    END DO

    RETURN
  END SUBROUTINE array_sort_r8

  SUBROUTINE array_locate_r4 ( n, x, x0, psel, ipos )

    ! -------------------------------------------------------------
    ! Given an ordered array X of size N and a point X0, return the
    ! position IPOS of X0 in X based on the selection PSEL.
    ! -------------------------------------------------------------

    USE OMSAO_precision_module, ONLY: i4, r4
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                INTENT (IN) :: n
    REAL    (KIND=r4), DIMENSION (n), INTENT (IN) :: x
    REAL    (KIND=r4),                INTENT (IN) :: x0
    CHARACTER (LEN=2),                INTENT (IN) :: psel

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: ipos

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: ilow, iupp, imid

    ! --------------------
    ! Initialize variables
    ! --------------------
    ipos = -1
    ilow = 0
    iupp = n+1
    imid = 0

    ! -------------------------------------
    ! The actual bisection. Short and sweet
    ! -------------------------------------
    DO WHILE ( iupp-ilow > 1 )
      imid = ( iupp + ilow ) / 2
      IF ( x(imid) >= x0 ) THEN
        iupp = imid
      ELSE
        ilow = imid
      END IF
    END DO

    ! -------------------------------------------------------
    ! The final value of IPOS depends on what we are actually
    ! looking for. Depending on our selection, we may have to
    ! return IPOS=-1 for cases that don't fit the selection.
    ! -------------------------------------------------------
    SELECT CASE ( psel )
    CASE ( 'LT' )
      ! X(IPOS)  < X0
      ipos = ilow
      IF ( x(ipos) == x0 ) ipos = ipos - 1
    CASE ( 'LE' )
      ! X(IPOS) <= X0
      ipos = ilow
      IF ( x(iupp) == x0 ) ipos = iupp
    CASE ( 'GE' )
      ! X(IPOS) >= X0
      ipos = iupp
    CASE ( 'GT' )
      ! X(IPOS)  > X0
      ipos = iupp
      IF ( x(ipos) == x0 ) ipos = ipos + 1
    CASE DEFAULT
      ! Return IPOS=IUP
      ipos = iupp
    END SELECT

    IF ( ipos <= 0 .OR. ipos > n ) ipos = -1

    RETURN
  END SUBROUTINE array_locate_r4

  SUBROUTINE array_locate_r8 ( n, x, x0, psel, ipos )

    ! -------------------------------------------------------------
    ! Given an ordered array X of size N and a point X0, return the
    ! position IPOS of X0 in X based on the selection PSEL.
    ! -------------------------------------------------------------

    USE OMSAO_precision_module, ONLY: i4, r8
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                INTENT (IN) :: n
    REAL    (KIND=r8), DIMENSION (n), INTENT (IN) :: x
    REAL    (KIND=r8),                INTENT (IN) :: x0
    CHARACTER (LEN=2),                INTENT (IN) :: psel

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: ipos

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: ilow, iupp, imid

    ! --------------------
    ! Initialize variables
    ! --------------------
    ipos = -1
    ilow = 0
    iupp = n+1
    imid = 0

    ! -------------------------------------
    ! The actual bisection. Short and sweet
    ! -------------------------------------
    DO WHILE ( iupp-ilow > 1 )
      imid = ( iupp + ilow ) / 2
      IF ( x(imid) >= x0 ) THEN
        iupp = imid
      ELSE
        ilow = imid
      END IF
    END DO

    ! -------------------------------------------------------------
    ! Force the Upper and Lower indices into the range of the array
    ! -------------------------------------------------------------
    ilow = MAX ( ilow, 1 )
    iupp = MIN ( iupp, n )

    ! -------------------------------------------------------
    ! The final value of IPOS depends on what we are actually
    ! looking for. Depending on our selection, we may have to
    ! return IPOS=-1 for cases that don't fit the selection.
    ! -------------------------------------------------------
    SELECT CASE ( psel )
    CASE ( 'LT' )
      ! X(IPOS)  < X0
      ipos = ilow
      IF ( x(ipos) == x0 ) ipos = ipos - 1
    CASE ( 'LE' )
      ! X(IPOS) <= X0
      ipos = ilow
      IF ( x(iupp) == x0 ) ipos = iupp
    CASE ( 'GE' )
      ! X(IPOS) >= X0
      ipos = iupp
    CASE ( 'GT' )
      ! X(IPOS)  > X0
      ipos = iupp
      IF ( x(ipos) == x0 ) ipos = ipos + 1
    CASE DEFAULT
      ! Return IPOS=IUP
      ipos = iupp
    END SELECT

    IF ( ipos <= 0 .OR. ipos > n ) ipos = -1

    RETURN
  END SUBROUTINE array_locate_r8

!  subroutine array_find_bounding_indices_r4 (n, xs, xmin, xmax, iminp, imaxp)
!
!    USE OMSAO_precision_module, ONLY: i4, r4
!    implicit none
!    integer (kind=i4), intent(in) :: n
!    real (kind=r4), intent(in), dimension(n) :: xs
!    real (kind=r4), intent(in) :: xmin, xmax
!    integer (kind=i4), intent(out) :: iminp, imaxp
!
!    ! xs[iminp] <= xmin
!    CALL array_locate_r4 (n, xs, xmin, 'LE', iminp)
!    ! xmax <= xs[imaxp]
!    CALL array_locate_r4 (n, xs, xmax, 'GE', imaxp)
!
!  end subroutine
!
!  subroutine array_find_inner_bounding_indices_r4 (n, xs, xmin, xmax, iminp, imaxp)
!
!    USE OMSAO_precision_module, ONLY: i4, r4
!    implicit none
!    integer (kind=i4), intent(in) :: n
!    real (kind=r4), intent(in), dimension(n) :: xs
!    real (kind=r4), intent(in) :: xmin, xmax
!    integer (kind=i4), intent(out) :: iminp, imaxp
!
!    ! xs[iminp] <= xmin
!    CALL array_locate_r4 (n, xs, xmin, 'GE', iminp)
!    ! xmax <= xs[imaxp]
!    CALL array_locate_r4 (n, xs, xmax, 'LE', imaxp)
!
!  end subroutine

  subroutine array_find_bounding_indices_r8 (n, xs, xmin, xmax, iminp, imaxp)

    USE OMSAO_precision_module, ONLY: i4, r8
    implicit none
    integer (kind=i4), intent(in) :: n
    real (kind=r8), intent(in), dimension(n) :: xs
    real (kind=r8), intent(in) :: xmin, xmax
    integer (kind=i4), intent(out) :: iminp, imaxp

    ! xs[iminp] <= xmin
    call array_locate_r8 (n, xs, xmin, 'LE', iminp)
    ! xmax <= xs[imaxp]
    call array_locate_r8 (n, xs, xmax, 'GE', imaxp)

  end subroutine

  subroutine array_find_inner_bounding_indices_r8 (n, xs, xmin, xmax, iminp, imaxp)

    USE OMSAO_precision_module, ONLY: i4, r8
    implicit none
    integer (kind=i4), intent(in) :: n
    real (kind=r8), intent(in), dimension(n) :: xs
    real (kind=r8), intent(in) :: xmin, xmax
    integer (kind=i4), intent(out) :: iminp, imaxp

    ! xs[iminp] <= xmin
    call array_locate_r8 (n, xs, xmin, 'GE', iminp)
    ! xmax <= xs[imaxp]
    call array_locate_r8 (n, xs, xmax, 'LE', imaxp)

  end subroutine

end module arrayutils
