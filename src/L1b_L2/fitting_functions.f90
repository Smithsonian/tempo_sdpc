MODULE fitting_functions

  USE OMSAO_precision_module, ONLY: i4, r8
  ! ----------------------------------------------------------------------------
  ! Number of calls to fitting function. This is counted in the fitting function
  ! itself, in an attempt to catch infinite loops of the ELSUNC routine, which
  ! occur on occasion for reasons not yet known.
  ! ----------------------------------------------------------------------------
  INTEGER (KIND=i4) :: num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls

  INTEGER (KIND=i4), PARAMETER, PRIVATE :: forever = HUGE(1_i4)

CONTAINS
  SUBROUTINE specfit_func_sol ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

    !
    ! Calculates the Solar spectrum and its derivatives for ELSUNC
    !
    ! NOTE: the variable DYDA as required here is the transpose of that
    !       rquired for the Numerical Recipes
    !

    USE OMSAO_elsunc_fitting_module, ONLY : &
      ELSUNC_PARSOOB_EVAL, ELSUNC_INFLOOP_EVAL
    USE OMSAO_indices_module, ONLY: MAX_CALFIT_IDX
    USE OMSAO_variables_module,  ONLY : &
      fitwavs, fitweights, currspec, sol_wav_avg, &
      lobnd, upbnd, mask_fitvar_cal, fitvar_cal
    USE spectra, ONLY: spectrum_solar

    IMPLICIT NONE

    INTEGER (KIND=i4),                              INTENT (IN)    :: nfitvar, npoints, mdy
    INTEGER (KIND=i4),                              INTENT (INOUT) :: ctrl
    REAL    (KIND=r8), DIMENSION (nfitvar),         INTENT (INOUT)    :: fitvar
    REAL    (KIND=r8), DIMENSION (npoints),         INTENT (INOUT) :: ymod
    REAL    (KIND=r8), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

    !REAL    (KIND=r8), DIMENSION (npoints) :: locwvl
    INTEGER (KIND=i4) :: abs_ctrl, i, idx

    abs_ctrl = ABS(ctrl)

    select case (abs_ctrl)
    case (1, 3)

      if (abs_ctrl == 1) then
        num_fitfunc_calls = num_fitfunc_calls + 1
        IF ( num_fitfunc_calls > max_fitfunc_calls ) THEN
          ctrl = ELSUNC_INFLOOP_EVAL
          RETURN
        END IF
        IF ( ANY(fitvar(1:nfitvar) < lobnd(1:nfitvar)) .OR. &
            ANY(fitvar(1:nfitvar) > upbnd(1:nfitvar)) ) THEN
          ctrl = ELSUNC_PARSOOB_EVAL
          RETURN
        END IF
      endif

      ! -----------------------------------------------------------------------------------
      ! First, we have to undo the compression of the FITVAR_SOL array. This compression
      ! is performed in the SOLAR_FIT and RADIANCE_WAVCAL subroutines and accelerates the
      ! fitting process, because ELSUNC has to handle less indices. But here we require the
      ! original layout, otherwise the index assingment is screwed.
      ! -----------------------------------------------------------------------------------
      DO i = 1, nfitvar
        idx = mask_fitvar_cal(i)
        fitvar_cal(idx) = fitvar(i)
      END DO

      !locwvl(1:npoints) = fitwavs(1:npoints)
      CALL spectrum_solar ( &
        npoints, sol_wav_avg, fitwavs(1:npoints), ymod(1:npoints), &
        fitvar_cal, MAX_CALFIT_IDX)

      ! -----------------------------------------------------------------------
      ! Calculate the weighted difference between fitted and measured spectrum.
      ! -----------------------------------------------------------------------
      if (abs_ctrl == 1)  ymod(1:npoints) &
        = fitweights(1:npoints) * (ymod(1:npoints) - currspec(1:npoints))
      RETURN

    case (2)
      ! -------------------------------------------------
      ! Count the number of calls to the fitting function
      ! with request for the Jacobian and terminate if we
      ! exceed the allowed maximum (just to be safe!).
      ! -------------------------------------------------
      num_fitfunc_jacobi = num_fitfunc_jacobi + 1
      IF ( num_fitfunc_jacobi > max_fitfunc_calls ) THEN
        ctrl = ELSUNC_INFLOOP_EVAL
      ELSE
      ! ---------------------------------------------------------------------
      ! The following sets up ELSUNC for numerical computation of the fitting
      ! function derivative. It is faster and more flexible than the original
      ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
      ! ---------------------------------------------------------------------
        ctrl = 0
      END IF

    end select

  END SUBROUTINE specfit_func_sol

  SUBROUTINE specfit_func ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

    !
    !     Calculates the spectrum and its derivatives for ELSUNC
    !
    ! NOTE: the variable DYDA as required here is the transpose of that
    !       rquired for the Numerical Recipes
    !

    USE OMSAO_elsunc_fitting_module, ONLY : &
      ELSUNC_INFLOOP_EVAL
    USE OMSAO_variables_module, ONLY : &
      yn_doas, fitwavs, fitweights, currspec, rad_wav_avg, &
      lobnd, upbnd !, num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls
    USE spectra, ONLY: spectrum_earthshine

    IMPLICIT NONE

    INTEGER (KIND=i4),                              INTENT (IN)    :: nfitvar, npoints, mdy
    INTEGER (KIND=i4),                              INTENT (INOUT) :: ctrl
    REAL    (KIND=r8), DIMENSION (nfitvar),         INTENT (INOUT) :: fitvar
    REAL    (KIND=r8), DIMENSION (npoints),         INTENT (INOUT) :: ymod
    REAL    (KIND=r8), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

    REAL    (KIND=r8), DIMENSION (npoints)       :: locwvl

    locwvl(1:npoints) = fitwavs(1:npoints)

    SELECT CASE ( ABS ( ctrl ) )
    CASE ( 1 )
      ! -------------------------------------------------
      ! Count the number of calls to the fitting function
      ! and terminate if we exceed the allowed maximum.
      ! -------------------------------------------------
      num_fitfunc_calls = num_fitfunc_calls + 1
      IF ( num_fitfunc_calls > max_fitfunc_calls ) THEN
        ctrl = INT(ELSUNC_INFLOOP_EVAL, KIND=i4)
        RETURN
      END IF
      ! -------------------------------------------------------------
      ! Check whether any of the fitting variables are out of bounds.
      ! If so, indicate uncomputability and return. "Uncomputability"
      ! is indicated by setting  CTRL to < -10.
      ! (NOTE that this is slightly different from the description in
      ! the ELSUNC fitting module).
      ! -------------------------------------------------------------
      IF ( ANY(fitvar(1:nfitvar) < lobnd(1:nfitvar)) .OR. &
          ANY(fitvar(1:nfitvar) > upbnd(1:nfitvar)) ) THEN
        ! DON'T USE THIS! ctrl = INT(ELSUNC_PARSOOB_EVAL, KIND=i4)
        ctrl = -1
        RETURN
      END IF

      ! -----------------------------------------------------------------------
      ! Calculate the weighted difference between fitted and measured spectrum.
      ! -----------------------------------------------------------------------
      CALL spectrum_earthshine ( &
        npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), &
        ymod(1:npoints), fitvar(1:nfitvar), & !database,
        yn_doas )
      !IF ( .NOT. yn_reference_fit ) THEN
      !   WRITE (99,'(3I6)') num_fitfunc_calls, nfitvar, npoints
      !   WRITE (99,'(1P100(E15.5:))') fitvar(1:nfitvar)
      !   DO i = 1, npoints
      !      WRITE (99,'(0PF12.4, 1P3E15.5)') locwvl(i), currspec(i), ymod(i), fitweights(i)
      !   END DO
      !END IF

      ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) * fitweights(1:npoints)

    CASE ( 2 )
      ! -------------------------------------------------
      ! Count the number of calls to the fitting function
      ! with request for the Jacobian and terminate if we
      ! exceed the allowed maximum (just to be safe!).
      ! -------------------------------------------------
      num_fitfunc_jacobi = num_fitfunc_jacobi + 1
      IF ( num_fitfunc_jacobi > max_fitfunc_calls ) THEN
        ctrl = INT(ELSUNC_INFLOOP_EVAL, KIND=i4)
        RETURN
      END IF
      ! ---------------------------------------------------------------------
      ! The following sets up ELSUNC for numerical computation of the fitting
      ! function derivative. It is faster and more flexible than the original
      ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
      ! ---------------------------------------------------------------------
      ctrl = 0; RETURN

    CASE ( 3 )
      ! Calculate the spectrum, without weighting
      CALL spectrum_earthshine ( &
        npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), &
        ymod(1:npoints), fitvar(1:nfitvar), & !database,
        yn_doas )

    CASE DEFAULT
      !WRITE (*,'(A,I4)') &
      !     "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE specfit_func

  SUBROUTINE specfit_func_o3exp ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

    !
    !     Calculates the spectrum and its derivatives for ELSUNC
    !
    ! NOTE: the variable DYDA as required here is the transpose of that
    !       rquired for the Numerical Recipes
    !

    USE OMSAO_elsunc_fitting_module, ONLY: ELSUNC_INFLOOP_EVAL
    USE OMSAO_variables_module, ONLY : &
      yn_doas, rad_wav_avg, fitwavs, fitweights, currspec, &
      lobnd, upbnd !, num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls
    USE spectra, ONLY: spectrum_earthshine_o3exp

    IMPLICIT NONE

    INTEGER (KIND=i4),                              INTENT (IN)    :: nfitvar, npoints, mdy
    INTEGER (KIND=i4),                              INTENT (INOUT) :: ctrl
    REAL    (KIND=r8), DIMENSION (nfitvar),         INTENT (INOUT) :: fitvar
    REAL    (KIND=r8), DIMENSION (npoints),         INTENT (INOUT) :: ymod
    REAL    (KIND=r8), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

    REAL    (KIND=r8), DIMENSION (npoints)       :: locwvl

    locwvl(1:npoints) = fitwavs(1:npoints)

    SELECT CASE ( ABS ( ctrl ) )
    CASE ( 1 )
      ! -------------------------------------------------
      ! Count the number of calls to the fitting function
      ! and terminate if we exceed the allowed maximum.
      ! -------------------------------------------------
      num_fitfunc_calls = num_fitfunc_calls + 1
      IF ( num_fitfunc_calls > max_fitfunc_calls ) THEN
        ctrl = INT(ELSUNC_INFLOOP_EVAL, KIND=i4)
        RETURN
      END IF
      ! -------------------------------------------------------------
      ! Check whether any of the fitting variables are out of bounds.
      ! If so, indicate uncomputability and return. "Uncomputability"
      ! is indicated by setting  CTRL to < -10.
      ! (NOTE that this is slightly different from the description in
      ! the ELSUNC fitting module).
      ! -------------------------------------------------------------
      IF ( ANY(fitvar(1:nfitvar) < lobnd(1:nfitvar)) .OR. &
          ANY(fitvar(1:nfitvar) > upbnd(1:nfitvar)) ) THEN
        ! DON'T USE THIS! ctrl = INT(ELSUNC_PARSOOB_EVAL, KIND=i4)
        ctrl = -1
        RETURN
      END IF

      ! -----------------------------------------------------------------------
      ! Calculate the weighted difference between fitted and measured spectrum.
      ! -----------------------------------------------------------------------
      CALL spectrum_earthshine_o3exp ( &
        npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), &
        ymod(1:npoints), fitvar(1:nfitvar), & !database,
        yn_doas )

      !WRITE (*,'(1P100(E12.4:))') fitvar(1:nfitvar)
      !IF ( .NOT. yn_reference_fit ) THEN
      !   WRITE (99,'(3I6)') num_fitfunc_calls, nfitvar, npoints
      !   WRITE (99,'(1P100(E15.5:))') fitvar(1:nfitvar)
      !   DO i = 1, npoints
      !      WRITE (99,'(0PF12.4, 1P3E15.5)') locwvl(i), currspec(i), ymod(i), fitweights(i)
      !   END DO
      !END IF

      ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) * fitweights(1:npoints)
    CASE ( 2 )
      ! -------------------------------------------------
      ! Count the number of calls to the fitting function
      ! with request for the Jacobian and terminate if we
      ! exceed the allowed maximum (just to be safe!).
      ! -------------------------------------------------
      num_fitfunc_jacobi = num_fitfunc_jacobi + 1
      IF ( num_fitfunc_jacobi > max_fitfunc_calls ) THEN
        ctrl = INT(ELSUNC_INFLOOP_EVAL, KIND=i4)
        RETURN
      END IF
      ! ---------------------------------------------------------------------
      ! The following sets up ELSUNC for numerical computation of the fitting
      ! function derivative. It is faster and more flexible than the original
      ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
      ! ---------------------------------------------------------------------
      ctrl = 0; RETURN

    CASE ( 3 )
      ! Calculate the spectrum, without weighting
      CALL spectrum_earthshine_o3exp ( &
        npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), &
        ymod(1:npoints), fitvar(1:nfitvar), & !database,
        yn_doas )

    CASE DEFAULT
      !WRITE (*,'(A,I4)') &
      !     "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE specfit_func_o3exp

  SUBROUTINE cubic_func ( x, afunc, ma )

    ! ***************************************************
    !
    !   Computes a third-order polynomial. Used in LFIT
    !
    ! ***************************************************

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER (KIND=i4), INTENT (IN) :: ma
    REAL    (KIND=r8), INTENT (IN) :: x

    ! ================
    ! Output variables
    ! ================
    REAL (KIND=r8), DIMENSION (ma), INTENT (OUT) :: afunc

    ! ===============
    ! Local variables
    ! ===============
    INTEGER (KIND=i4) :: i

    afunc(1) = 1.0_r8
    DO i = 2, ma
      afunc(i) = afunc(i-1) * x
    END DO

    RETURN
  END SUBROUTINE cubic_func

!---------------------------------------------------------------------------

  SUBROUTINE spec_fit (                                       &
      nfitvar, fitvar, nspecpts, lowbnd, uppbnd, max_itnum, &
      covar, fitspec, fitres, exval, itnum, fitfunc           )

    USE OMSAO_indices_module,        ONLY: elsunc_userdef
    USE OMSAO_elsunc_fitting_module, ONLY: ELSUNC_NP, ELSUNC_NW
    USE OMSAO_variables_module,      ONLY: &
      tol, epsrel, epsabs, epsx !, num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls
    USE OMSAO_elsunc_fitting_module, ONLY: elsunc

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER (KIND=i4),                       INTENT (IN) :: nfitvar, nspecpts, max_itnum
    REAL    (KIND=r8), DIMENSION (nfitvar),  INTENT (IN) :: lowbnd, uppbnd

    ! ==================
    ! Modified variables
    ! ==================
    REAL (KIND=r8), DIMENSION (nfitvar),  INTENT (INOUT) :: fitvar

    ! ================
    ! Output variables
    ! ================
    INTEGER (KIND=i4),                              INTENT (OUT) :: exval, itnum
    REAL    (KIND=r8), DIMENSION (nfitvar,nfitvar), INTENT (OUT) :: covar
    REAL    (KIND=r8), DIMENSION (nspecpts),        INTENT (OUT) :: fitres, fitspec

    ! ===============
    ! Local variables
    ! ===============
    INTEGER (KIND=i4)                               :: elbnd
    INTEGER (KIND=i4), DIMENSION (ELSUNC_NP)        :: p
    REAL    (KIND=r8), DIMENSION (ELSUNC_NW)        :: w
    REAL    (KIND=r8), DIMENSION (nfitvar)          :: blow, bupp
    REAL    (KIND=r8), DIMENSION (nspecpts)         :: f
    REAL    (KIND=r8), DIMENSION (nspecpts,nfitvar) :: dfda
    INTEGER :: exit_val_local

    EXTERNAL fitfunc

    ! ===============================================================
    ! ELBND: 0 = unconstrained
    !        1 = all variables have same lower bound
    !        else: lower and upper bounds must be supplied by the use
    ! ===============================================================
    elbnd = elsunc_userdef

    exit_val_local = 0 ; itnum = 0

    p   = -1    ;  p(1)   = 0  ;  p(3) = max_itnum
    w   = -1.0  ;  w(1:4) = (/ tol,  epsrel,  epsabs,  epsx /)

    blow(1:nfitvar) = lowbnd(1:nfitvar)
    bupp(1:nfitvar) = uppbnd(1:nfitvar)

    ! ---------------------------------------------------------------------------------
    ! Reset to ZERO the number of calls to the fitting function, NUM_FITFUNC_CALLS.
    ! This variable is incremented with each call to the fitting function and
    ! checked agains MAX_FITFUNC_CALLS. Exceeded function calls lead to uncomputability
    ! followed by termination of the iteration process.
    ! ---------------------------------------------------------------------------------
    num_fitfunc_calls = 0 ; num_fitfunc_jacobi = 0
    max_fitfunc_calls = max_itnum * nfitvar * nfitvar ! Empiric value
    IF ( max_itnum < 0 ) max_fitfunc_calls = forever

    ! Use a local variable for the exit value since elsunc uses a integer, which
    ! may be a different size from i4.
    CALL elsunc ( &
      fitvar(1:nfitvar), nfitvar, nspecpts, nspecpts, fitfunc, elbnd, blow(1:nfitvar), &
      bupp(1:nfitvar), p, w, exit_val_local, &
      f(1:nspecpts), dfda(1:nspecpts,1:nfitvar) )

    exval = exit_val_local

    ! -----------------------------
    ! Save the number of iterations
    ! -----------------------------
    itnum = p(6)

    ! ------------------------------------------------------------------
    ! Call to ELSUNC fitting function to obtain complete fitted spectrum
    ! ------------------------------------------------------------------
    CALL fitfunc ( fitvar(1:nfitvar), nfitvar, fitspec(1:nspecpts), nspecpts, 3, dfda, 0 )

    ! ---------------------------------------------------------------
    ! Compute fitting residual
    ! FITRES is the negative of the returned function F = Model-Data.
    ! ---------------------------------------------------------------
    fitres(1:nspecpts) = -f(1:nspecpts)

    ! Covariance matrix.
    ! ------------------
    covar(1:nfitvar,1:nfitvar) = dfda(1:nfitvar,1:nfitvar)

    RETURN
  END SUBROUTINE spec_fit

END MODULE
