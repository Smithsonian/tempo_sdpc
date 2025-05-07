!> Two-step retrieval routines
module m_twostep

  public step2_specfit, twostep_inversion
  private

contains

  SUBROUTINE twostep_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, &
       epsrel, last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, &
       lowbnd, upbnd, sa, varname, ffidx, flidx, delta_x, covar, ncovar, &
       conv, avg_kernel,  contri, ozdfs, ozinfo, lchisq)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,  ONLY: elsunc_np, elsunc_nw!, &
    !maxlay, max_spec_pts
    USE ozprof_data_module,       ONLY: fgasidxs, ngas, do_bothstep!, gasidxs
    USE OMSAO_variables_module,   ONLY: tol, epsabs, epsx, step2_y, &
         step2_dyda!, max_itnum_rad, mask_fitvar_rad, fitvar_rad, fitwavs
    !USE OMSAO_indices_module,     ONLY: elsunc_userdef
    USE bounded_nonlin_LS,        ONLY: elsunc
    USE OMSAO_errstat_module
    use m_oe_inversion


    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT (IN)                             :: ns, nf, ozwrtint_unit,&
         num_iter, ffidx, flidx
    REAL (KIND=dp), INTENT (IN)                      :: epsrel
    REAL (KIND=dp), INTENT (OUT)                     :: ozdfs, ozinfo, lchisq
    REAL (KIND=dp), DIMENSION(nf), INTENT (IN)       :: lowbnd, upbnd, xap, &
         xold
    REAL (KIND=dp), DIMENSION(nf), INTENT (OUT)      :: delta_x
    REAL (KIND=dp), DIMENSION(ns), INTENT (IN)       :: gspec, sig
    CHARACTER (LEN=6), DIMENSION(nf), INTENT (IN)    :: varname
    REAL (KIND=dp), DIMENSION(ns, nf), INTENT (INOUT):: dyda
    REAL (KIND=dp), DIMENSION(nf, ns), INTENT (OUT)  :: contri
    REAL (KIND=dp), DIMENSION(nf, nf), INTENT(INOUT) :: sa
    REAL (KIND=dp), DIMENSION(nf, nf), INTENT(OUT)   :: covar, ncovar, &
         avg_kernel
    LOGICAL, INTENT(IN)                              :: do_sa_diagonal, &
         ozwrtint, last_iter
    LOGICAL, INTENT(OUT)                             :: conv

    ! =======================
    ! Local variables
    ! =======================
    INTEGER                                          :: i, j, idx, n2f, &
         elbnd, exval
    INTEGER, DIMENSION(nf)                           :: step2idxs
    REAL (KIND=dp), DIMENSION(nf, nf)                :: sa_sav
    REAL (KIND=dp), DIMENSION(nf)                    :: step2_fitvar
    CHARACTER (LEN=6), DIMENSION(nf)                 :: step2_varname
    REAL (KIND=dp), DIMENSION(ns)                    :: gspec1

    INTEGER,        DIMENSION (elsunc_np)            :: p
    REAL (KIND=dp), DIMENSION (elsunc_nw)            :: w
    REAL (KIND=dp), DIMENSION (nf)                   :: blow, bupp !, stderr
    REAL (KIND=dp), DIMENSION (ns)                   :: step2_fitres
    REAL (KIND=dp), DIMENSION (ns, nf)               :: dfda
    !REAL (KIND=dp), DIMENSION (nf, nf)               :: correl

    !EXTERNAL step2_specfit

    ! Freeze trace gas variables other than ozone by either
    ! 1. Setting weighting function to zero
    ! 2. Setting the a priori covariance matrix to zero
    sa_sav = sa

    n2f = 0
    DO i = 1, ngas
      IF (fgasidxs(i) > 0) THEN
        idx = fgasidxs(i)
        IF (.NOT. do_bothstep) sa(idx, idx) = 0.0
        n2f = n2f + 1
        step2idxs(n2f) = idx
      ENDIF
    ENDDO

    !DO i = flidx + 1, nf
    !   IF (varname(i)(3:4) == 'a1' .OR. varname(i)(3:4) == 'a2' .OR. varname(i)(3:4) == 'a3') THEN
    !      IF (.NOT. do_bothstep) sa(i, i) = 0.0
    !      n2f = n2f + 1
    !      step2idxs(n2f) = i
    !   ENDIF
    !ENDDO

    IF (n2f > 0) THEN
      step2_fitvar(1:n2f)  = xold(step2idxs(1:n2f))
      step2_dyda(1:ns, 1:n2f) = dyda(:, step2idxs(1:n2f))
      step2_varname(1:n2f) = varname(step2idxs(1:n2f))
    ENDIF

    delta_x = 0.0

    ! Call optimal estimation routine for those un-frozen parameters
    CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel, &
         last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa, &
         varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), &
         ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf), &
         contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec1)
    !WRITE(www_lun, '(I3,1X,A6,2d14.6)') ((i, varname(i), xold(i), &
    !     delta_x(i)), i=1, nf)

    IF (n2f <= 0 ) RETURN

    ! Second step: use NLLS (ELSUNC, but linear here, i.e., 1 iterations) to
    ! fit those freezed trace gases
    elbnd = 0   ! Unconstrained
    exval = 0

    p  = -1   ;          p(1)   = 0;     p(3) = 1   !max_itnum_rad
    w  = -1.0 ;          w(1:4) = (/ tol,  epsrel,  epsabs,  epsx /)
    blow(1:n2f) = -1.0D99;   bupp(1:n2f) = 1.0D99
    blow(1:n2f) = -1.0D99;   bupp(1:n2f) = 1.0D99
    !blow(2) = 0.0; bupp(2) = 1.0
    !blow(5) = 0.0; bupp(5) = 1.0

    step2_fitres = 0.0;  dfda(1:ns, 1:n2f) = 0.0
    step2_y(1:ns) = gspec1(1:ns)

    CALL elsunc ( step2_fitvar(1:n2f), n2f, ns, ns, step2_specfit,          &
         elbnd, blow(1:n2f),  bupp(1:n2f), p, w, exval, step2_fitres(1:ns), &
         dfda(1:ns, 1:n2f) )

    ! ---------------------------------------------------------------
    ! Compute fitting residual
    ! FITRES is the negative of the returned function F = Model-Data.
    ! -------------------------------------opl--------------------------
    step2_fitres(1:ns) = -step2_fitres(1:ns)

    !WRITE(90, '(3D14.6)') ((gspec(i), gspec1(i), step2_fitres(i)), i=1, ns)

    ! Fitting RMS and CHI**2
    ! ----------------------
    lchisq = SQRT (SUM (step2_fitres(1:ns)**2 ) / ns )

    !! compute standard deviation for each variable
    !DO i = 1, n2f
    !   stderr(i) =  SQRT(dfda(i, i) * ns / (ns - n2f)) !* lchisq
    !END DO

    delta_x(step2idxs(1:n2f)) = delta_x(step2idxs(1:n2f)) + step2_fitvar(1:n2f)

    DO i = 1, n2f
      DO j = 1, n2f
        covar(step2idxs(i), step2idxs(j)) = covar(step2idxs(i), &
             step2idxs(j)) + dfda(i, j)
        ncovar(step2idxs(i), step2idxs(j)) = ncovar(step2idxs(i), &
             step2idxs(j)) + dfda(i, j)
      ENDDO
    ENDDO

    DO i = 1, ngas
      IF (fgasidxs(i) > 0) THEN
        idx = fgasidxs(i)
        sa(idx, idx) = sa_sav(idx, idx)
        !IF (ABS(xold(idx) + delta_x(idx)) > SQRT(sa(idx, idx))) &
        !     sa(idx, idx)= (xold(idx) + delta_x(idx))**2
      ENDIF
    ENDDO

    !WRITE(www_lun, '(I3,1X,A6,2d14.6)') ((i, varname(i), xold(i), delta_x(i)), i=1, nf)
    !print *, lchisq

    RETURN
  END SUBROUTINE twostep_inversion

  SUBROUTINE step2_specfit ( a, na, y, m, ctrl, dyda, mdy )

    USE OMSAO_precision_module
    !USE OMSAO_parameters_module, ONLY : max_spec_pts
    USE OMSAO_variables_module,  ONLY : step2_y, step2_dyda
    USe OMSAO_errstat_module

    IMPLICIT NONE

    ! Input parameters
    ! ================
    INTEGER,                         INTENT (IN)  :: na, m, mdy
    REAL (KIND=dp), DIMENSION (na),  INTENT (IN)  :: a

    ! Modified parameters
    ! ===================
    INTEGER, INTENT (INOUT) :: ctrl

    ! Output parameters
    ! =================
    REAL (KIND=dp), DIMENSION (m),     INTENT (OUT) :: y
    REAL (KIND=dp), DIMENSION (m, na), INTENT (OUT) :: dyda

    ! Local variables
    ! ===============
    REAL (KIND=dp), DIMENSION (m) :: y0
    INTEGER                       :: i

    y0 = 0.0_dp
    DO i = 1, na
      y0 =  y0 + a(i) * (step2_dyda(1:m, i))
    END DO

    SELECT CASE ( ABS(ctrl) )
    CASE ( 1 )
      y  =  y0 - step2_y(1:m)
    CASE ( 2 )
      dyda(1:m, 1:na) = step2_dyda(1:m, 1:na)
    CASE ( 3 )
      ! This CASE is included to get the complete fitted spectrum
      y  = y0
    CASE DEFAULT
      WRITE(www_lun, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE step2_specfit


end module m_twostep
