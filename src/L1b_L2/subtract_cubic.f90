MODULE subtract_cubic
  USE OMSAO_precision_module, ONLY: r8, i4
  use optimizer_interface_module
  use errormodule

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: cubic_x, cubic_y, cubic_w

  PRIVATE cubic_x, cubic_y, cubic_w, r8, i4

CONTAINS
  subroutine eval_cubic (a, x, y)
    implicit none
    real (kind=r8), dimension(:), intent(in) :: a
    real (kind=r8), dimension(:), intent(in) :: x
    real (kind=r8), dimension(:), intent(out) :: y
    ! do it this way to reproduce the previous answer  ! FIXME!!
    y = a(1) + a(2)*x + a(3)*x*x + a(4)*x*x*x
    !y = a(1) + x*(a(2) + x*(a(3) + x*a(4)))      ! a better way
  end subroutine eval_cubic

  SUBROUTINE cubic_objective (this, a, na, y, m, return_status)
    implicit none
    type(optimizer_type) :: this
    real (kind=r8), dimension(:), intent(in) :: a
    real (kind=r8), dimension(:), intent(out) :: y
    integer (kind=i4), intent(in) :: na, m
    integer (kind=i4), intent(out) :: return_status
    ! local
    real (kind=r8), dimension(m) :: y0

    call eval_cubic (a, cubic_x(1:m), y0)
    y = (y0 - cubic_y(1:m)) / cubic_w(1:m)
    return_status = 0

  end subroutine cubic_objective

  SUBROUTINE cubic_subtract ( locwvl, npts, ll_rad, lu_rad, errstat )

    USE OMSAO_indices_module,        ONLY : max_rs_idx, solar_idx
    USE OMSAO_parameters_module,     ONLY : max_spec_pts, doas_npol
    USE OMSAO_variables_module,      ONLY : database

    IMPLICIT NONE

    INTEGER (KIND=i4),                   INTENT (IN) :: npts, ll_rad, lu_rad
    REAL    (KIND=r8), DIMENSION (npts), INTENT (IN) :: locwvl
    integer, intent(inout) :: errstat

    INTEGER (KIND=i4)                   :: i, nlower, nupper, nfitted
    REAL    (KIND=r8)                   :: locavg, chisq
    REAL    (KIND=r8), DIMENSION (npts) :: x, ptemp, sig

    ! optimization variables
    INTEGER (KIND=i4)                                     :: exval
    REAL    (KIND=r8), DIMENSION (doas_npol)              :: blow, bupp
    REAL    (KIND=r8), DIMENSION (max_spec_pts)           :: f
    REAL    (KIND=r8), DIMENSION (doas_npol)              :: par
    type(optimizer_type) :: opt
    integer (kind=i4) :: return_status

    if (errstat < 0) return

    ! ======================
    ! Assign fitting weights
    ! ======================
    sig = 1.0_r8

    ! -------------------------------------------------------------------------
    ! Find limits for polynomial fitting, with ~1 nm overlap
    ! -------------------------------------------------------------------------
    ! ARRAY_LOCATE is the preferred way to calculate the bounds; however, this
    ! particular  section of the code needs still to be tested. tpk 09 Feb 2007
    ! -------------------------------------------------------------------------
    nlower = MINVAL ( MINLOC ( locwvl(1:npts), MASK=(locwvl(1:npts) >= locwvl(ll_rad)-1.0_r8) ) )
    nupper = MAXVAL ( MAXLOC ( locwvl(1:npts), MASK=(locwvl(1:npts) <= locwvl(lu_rad)+1.0_r8) ) )
    !CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl(ll_rad)-1.0_r8, 'GE', nlower )
    !CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl(lu_rad)+1.0_r8, 'LE', nupper )

    nfitted = nupper - nlower + 1

    ALLOCATE (cubic_x(nfitted))
    ALLOCATE (cubic_y(nfitted))
    ALLOCATE (cubic_w(nfitted))

    ! Find average position over fitted region
    locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=r8 )

    ! Load temporary position file: re-define positions in order to fit
    ! about mean position
    DO i = 1, nfitted
      ptemp(i) = locwvl(i+nlower-1) - locavg
    END DO

    call optimizer_open (opt, cubic_objective, doas_npol, return_status, &
                         mode=opt_unbounded, max_num_iterations=5)
    if (return_status < 0) then
      call err_message_error ("cubic_subtract:  optimizer_open failed", errstat)
      return
    endif

    !     Load and fit database spectra nos. 2-11
    DO i = 1, max_rs_idx
      IF ( i /= solar_idx ) THEN
        exval = 0
        blow(1:doas_npol) = 0.0_r8  ;  bupp(1:doas_npol) = 0.0_r8

        cubic_x(1:nfitted) = ptemp(1:nfitted)
        cubic_y(1:nfitted) = database(1+nlower-1:nfitted+nlower-1,i)
        cubic_w(1:nfitted) = sig(1:nfitted)

        par = 0.0_r8 ; f = 0.0_r8
        call opt%optimize (opt, par, doas_npol, f(1:nfitted), nfitted, exval)
        chisq = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines
        x(1:npts) = locwvl(1:npts) - locavg

        cubic_x(1:npts) = x(1:npts)
        call eval_cubic (par(1:doas_npol), cubic_x(1:npts), cubic_y(1:npts))
        database(1:npts,i) = database(1:npts,i) - cubic_y(1:npts)
      END IF
    END DO

    call optimizer_close (opt, return_status)
    if (return_status < 0) then
      call err_message_error ("cubic_subtract:  optimizer_close failed", errstat)
      return
    endif

    DEALLOCATE (cubic_x)
    DEALLOCATE (cubic_y)
    DEALLOCATE (cubic_w)

    RETURN

  END SUBROUTINE cubic_subtract

  SUBROUTINE cubic_subtract_meas ( locwvl, npoints, locspec, ll_rad, lu_rad, errstat )

    USE OMSAO_parameters_module,     ONLY : doas_npol
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: npoints, ll_rad, lu_rad
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: locwvl

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (npoints), INTENT (INOUT) :: locspec
    integer, intent(inout) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    REAL    (KIND=r8), DIMENSION (npoints)             :: tmp, ptmp, sig
    REAL    (KIND=r8), DIMENSION (doas_npol)           :: r
    INTEGER (KIND=i4)                                  :: i, nlower, nupper, nfitted
    REAL    (KIND=r8)                                  :: locavg
    REAL    (KIND=r8), DIMENSION (npoints)             :: x

    ! optimization variables
    INTEGER (KIND=i4)                                :: exval
    REAL    (KIND=r8)                                :: chisq2
    REAL    (KIND=r8), DIMENSION (doas_npol)         :: blow, bupp
    REAL    (KIND=r8), DIMENSION (npoints)           :: f
    type(optimizer_type) :: opt
    integer (kind=i4) :: return_status

    if (errstat < 0) return

    ! ======================
    ! Assign fitting weights
    ! ======================
    sig = 1.0_r8

    ! -------------------------------------------------------------------------
    !     Find limits for polynomial fitting, with ~1 nm overlap
    ! -------------------------------------------------------------------------
    ! ARRAY_LOCATE is the preferred way to calculate the bounds; however, this
    ! particular  section of the code needs still to be tested. tpk 09 Feb 2007
    ! -------------------------------------------------------------------------
    nlower = MINVAL(MINLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) >= locwvl(ll_rad)-1.0_r8) ))
    nupper = MAXVAL(MAXLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) <= locwvl(lu_rad)+1.0_r8) ))
    !CALL array_locate_r8 ( npoints, locwvl(1:npoints), locwvl(ll_rad)-1.0_r8, 'GE', nlower )
    !CALL array_locate_r8 ( npoints, locwvl(1:npoints), locwvl(lu_rad)+1.0_r8, 'LE', nupper )
    nfitted = nupper - nlower + 1

    !     Find average position over fitted region
    locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=r8 )

    !     Load temporary position file: re-define positions in order to fit
    !     about mean position
    DO i = 1, nfitted
      ptmp(i) = locwvl(i+nlower-1) - locavg
    END DO

    call optimizer_open (opt, cubic_objective, doas_npol, return_status, &
                         mode=opt_unbounded, max_num_iterations=5)
    if (return_status < 0) then
      call err_message_error ("cubic_subtract_meas:  optimizer_open failed", errstat)
      return
    endif

    !     Load and fit spectrum
    tmp(1:nfitted) = locspec(1+nlower-1:nfitted+nlower-1)

    blow(1:doas_npol) = 0.0_r8  ;  bupp(1:doas_npol) = 0.0_r8

    ALLOCATE (cubic_x(nfitted))
    ALLOCATE (cubic_y(nfitted))
    ALLOCATE (cubic_w(nfitted))

    cubic_x(1:nfitted) = ptmp(1:nfitted)
    cubic_y(1:nfitted) =  tmp(1:nfitted)
    cubic_w(1:nfitted) =  sig(1:nfitted)

    exval = 0
    r = 0.0_r8 ; f = 0.0_r8

    call opt%optimize (opt, r, doas_npol, f(1:nfitted), nfitted, exval)
    chisq2 = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines

    call optimizer_close (opt, return_status)
    if (return_status < 0) then
      call err_message_error ("cubic_subtract_meas:  optimizer_close failed", errstat)
      return
    endif

    ! Re-load spec with high-pass filtered data, over whole  spectral region
    x(1:npoints) = locwvl(1:npoints) - locavg

    cubic_x(1:npoints) = x(1:npoints)
    call eval_cubic (r(1:doas_npol), cubic_x(1:npoints), cubic_y(1:npoints))
    locspec(1:npoints) = locspec(1:npoints) - cubic_y(1:npoints)

    DEALLOCATE(cubic_x)
    DEALLOCATE(cubic_y)
    DEALLOCATE(cubic_w)

    RETURN

  END SUBROUTINE cubic_subtract_meas

END MODULE

