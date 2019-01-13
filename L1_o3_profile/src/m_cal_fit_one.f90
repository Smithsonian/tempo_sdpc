!
module m_cal_fit_one
  USE bounded_nonlin_LS, ONLY: elsunc  
  USE m_spectra, ONLY:spectrum_solar
  USE m_fitting_util, ONLY: poly_fit
  public cal_fit_one
  private specfit, specfit_func_sol
contains

  ! ********************************************************************
  ! Note: make sure currspec, fitwavs, fitweights, slitcal are set properly
  ! ********************************************************************
  SUBROUTINE cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, &
       wrt_to_file, slitcal, slit_unit, avgwav, varstd, solfit_exval)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,     ONLY : max_calfit_idx, solar_idx, & 
         hwe_idx, asy_idx, shi_idx, squ_idx, wvl_idx, spc_idx, sin_idx, &
         vgl_idx, vgr_idx, hwl_idx, hwr_idx, spk_idx, &
         bl0_idx, bl1_idx, bl2_idx, bl3_idx, bl4_idx, bl5_idx, bl6_idx, bl7_idx, &
         sc0_idx, sc1_idx, sc2_idx, sc3_idx, sc4_idx, sc5_idx, sc6_idx, sc7_idx, &
         wr0_idx, wr1_idx, wr2_idx, wr3_idx, wr4_idx, wr5_idx, wr6_idx, wr7_idx

    USE OMSAO_variables_module,   ONLY: weight_sun,  &
        fitvar_sol, fitvar_sol_saved,chisq,  sol_wav_avg, which_slit, &
        mask_fitvar_sol, rmask_fitvar_sol, fitwavs, fitweights, currspec,& 
        lo_sunbnd, up_sunbnd,max_itnum_sol, poly_order, &
        n_refspec_pts, refspec_orig_data,tmp_rad
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ================
    ! Input variables
    ! ================
    INTEGER, INTENT(IN) :: slit_unit, n_fit_pts, n_fitvar_sol
    LOGICAL, INTENT(IN) :: wrt_to_screen, wrt_to_file, slitcal

    ! ================
    ! Output variables
    ! ================
    INTEGER,            INTENT (OUT)              :: solfit_exval
    REAL (KIND = dp),   INTENT (OUT)              :: avgwav  
    REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: varstd

    ! ===============
    ! Local variables
    ! ===============
    REAL (KIND=dp)  :: asum, ssum, rms,  niter
    REAL (KIND=dp), DIMENSION (n_fit_pts)                 :: fitres, fitspec, polyx
    REAL (KIND=dp), DIMENSION (n_fitvar_sol, n_fitvar_sol):: covar
    REAL (KIND=dp)  :: hw1e, e_asym, vgl, vgr, hwl, hwr, sin, shi, squ, spk
    REAL (KIND=dp), DIMENSION (n_fitvar_sol) :: fitvar, lobnd, upbnd, stderr
    REAL (KIND=dp), DIMENSION (8) :: polycoeffs
    INTEGER         :: i, j, ref_fidx, ref_lidx, ref_all_pts, ll, lu

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=11), PARAMETER :: modulename = 'cal_fit_one'

    ! Initialization for wavelength registration
    IF (ANY(rmask_fitvar_sol(wr0_idx:wr7_idx) > 0)) THEN ! updated v2
       DO i = 1, n_fit_pts
        polyx(i) = 1.0d0 * i
       ENDDO
        polyx = (polyx - n_fit_pts / 2.0) / n_fit_pts

       DO i = wr0_idx, wr7_idx
        IF ( fitvar_sol(i) > lo_sunbnd(i) .and. fitvar_sol(i) < up_sunbnd(i) ) THEN
           poly_order = i - wr0_idx + 1
        ENDIF
       ENDDO
       ll = 1; lu = n_fit_pts
       CALL poly_fit(polyx, n_fit_pts, fitwavs(1:n_fit_pts), ll, lu, polycoeffs(1:poly_order))
       j = 1
       DO i = wr0_idx, wr7_idx
          fitvar_sol(i) = polycoeffs(j)
          lo_sunbnd(i) = -1.0D+99
          up_sunbnd(i) =  1.0D+99
          j = j + 1
       ENDDO
    ENDIF

    ! -------------------------------------------------------
    ! Calculate SOL_WAV_AVG of measured solar spectra here,
    ! for use in calculated spectra.
    ! -------------------------------------------------------
    IF ( weight_sun ) THEN
        asum = SUM ( fitwavs(1:n_fit_pts) / (fitweights(1:n_fit_pts) &
           * fitweights(1:n_fit_pts)))
        ssum = SUM (1.0/ (fitweights(1:n_fit_pts)*fitweights(1:n_fit_pts)))
        sol_wav_avg = asum / ssum
     ELSE
        sol_wav_avg = ( fitwavs(1) + fitwavs(n_fit_pts)) / 2.0
     END IF

     ! --------------------------------------------------------------
     ! ELSUNC FIT: Calculate and iterate on the irradiance spectrum.
     ! --------------------------------------------------------------     
     ! initialize the sin variables for faster convergence
     ref_all_pts = n_refspec_pts(solar_idx)   
     ref_fidx=MINVAL(MAXLOC(refspec_orig_data(solar_idx, 1:ref_all_pts, wvl_idx),&
          MASK=(refspec_orig_data(solar_idx, 1:ref_all_pts, wvl_idx) <= fitwavs(1))))
  
     ref_lidx=MINVAL(MINLOC(refspec_orig_data(solar_idx, 1:ref_all_pts, wvl_idx),&
          MASK=(refspec_orig_data(solar_idx,1:ref_all_pts,wvl_idx)>=fitwavs(n_fit_pts))))
 
     IF (ref_fidx > ref_all_pts .OR. ref_fidx <= 0 .OR. ref_lidx > ref_all_pts &
          .OR. ref_lidx <= 0 .OR. ref_fidx > ref_lidx) THEN
        WRITE(www_lun, *) ref_fidx, ref_lidx, n_fit_pts
        WRITE(www_lun, *) fitwavs(1), fitwavs(n_fit_pts)
        WRITE(www_lun, *) refspec_orig_data(solar_idx, ref_fidx, wvl_idx), refspec_orig_data(solar_idx, ref_lidx, wvl_idx)
        WRITE(www_lun, *) modulename,' : Solar spectra not cover fitting window!!!'; STOP
     END IF
     
     fitvar_sol(sin_idx) = SUM(currspec(1:n_fit_pts)) * (ref_lidx-ref_fidx + 1.0) / &
          SUM(refspec_orig_data(solar_idx, ref_fidx:ref_lidx, spc_idx)) /n_fit_pts

     fitvar = 0.0; lobnd = 0.0; upbnd = 0.0
     fitvar(1:n_fitvar_sol) = fitvar_sol(mask_fitvar_sol(1:n_fitvar_sol))
     lobnd(1:n_fitvar_sol) = lo_sunbnd(mask_fitvar_sol(1:n_fitvar_sol)) 
     upbnd(1:n_fitvar_sol) = up_sunbnd(mask_fitvar_sol(1:n_fitvar_sol))  
  
  
!  CALL specfit ( &
!       n_fitvar_sol, fitvar(1:n_fitvar_sol), n_fitvar_sol, n_fit_pts, &
!       lobnd(1:n_fitvar_sol), upbnd(1:n_fitvar_sol), max_itnum_sol,   &
!       rms, chisq, covar(1:n_fitvar_sol,1:n_fitvar_sol), fitspec(1:n_fit_pts), &
!       fitres(1:n_fit_pts), stderr, solfit_exval, specfit_func_sol)
   CALL specfit ( n_fitvar_sol, fitvar(1:n_fitvar_sol), &
         n_fit_pts, lobnd(1:n_fitvar_sol), upbnd(1:n_fitvar_sol), &
         max_itnum_sol, rms, chisq, covar(1:n_fitvar_sol,1:n_fitvar_sol), &
         fitspec(1:n_fit_pts), fitres(1:n_fit_pts), stderr, solfit_exval, &
         specfit_func_sol)


     tmp_rad(1:n_fit_pts) = fitspec(1:n_fit_pts)
     !DO i = 1, n_fit_pts
     !   WRITE(90, '(f10.4, 3D14.6)') fitwavs(i), currspec(i), fitspec(i), fitres(i)
     !ENDDO
     fitvar_sol_saved = fitvar_sol

     !standard variables variables and standard error
     varstd(1:max_calfit_idx, 1:2) = 0.0
     avgwav = ( fitwavs(1) + fitwavs(n_fit_pts)) / 2.0
  
     DO i = 1, max_calfit_idx
        varstd(i, 1) = fitvar_sol(i)
     ENDDO

     DO i = 1, n_fitvar_sol
        varstd(mask_fitvar_sol(i), 2) = stderr(i)
     END DO

     IF (wrt_to_screen .OR. (wrt_to_file .AND. solfit_exval > 0)) THEN
        IF (which_slit == 3) THEN ! triangle (symmetric)
           hw1e = fitvar_sol(hwe_idx); e_asym = 0.0
        ELSE IF (which_slit == 5) THEN ! OMI-preflight_slit
           hw1e = 0.0; e_asym = 0.0; spk = 0.0
        ELSE IF (which_slit == 2) THEN ! Voigt
           vgl  = fitvar_sol(vgl_idx);  vgr    = fitvar_sol(vgr_idx)
           hwl  = fitvar_sol(hwl_idx);  hwr    = fitvar_sol(hwr_idx)
        ELSE
           hw1e = fitvar_sol(hwe_idx); e_asym =  fitvar_sol(asy_idx)
           spk =  fitvar_sol(spk_idx)
        ENDIF
          shi = fitvar_sol(shi_idx); squ = fitvar_sol(squ_idx); sin = fitvar_sol(sin_idx)
     ENDIF
     IF (wrt_to_screen) THEN
        IF ( which_slit /= 2) THEN
           WRITE(*, '(5(A, 1pd14.6), 2(A, I6))') 'wav=',avgwav ,'hwle = ',  hw1e,  &
                ' e_asym =', e_asym, ' spk = ', spk, ' rms = ', rms, ' exval = ', solfit_exval, ' niter= ', niter
        ELSE
           WRITE(www_lun, '(5(A,1pd14.6),A,I6)') 'vgl = ',  vgl, ' vgr = ',  vgr,  &
                ' hwl = ', hwl, ' hwr = ', hwr,' rms = ',rms,' exval = ', solfit_exval
        END IF
  
        WRITE(*, '(3(A, 1pd14.6))') 'shi = ', shi, ' squ =', squ, ' sin = ', sin
        WRITE(*, '(4(A, 1pd14.6))') 'b10 = ', fitvar_sol(bl0_idx), ' b11 = ', &
             fitvar_sol(bl1_idx),' b12 = ', fitvar_sol(bl2_idx), &
             ' b13 = ', fitvar_sol(bl3_idx)
        WRITE(*, '(4(A, 1pd14.6))') 'sc0 = ', fitvar_sol(sc0_idx), ' sc1 = ', &
             fitvar_sol(sc1_idx), ' sc2 = ', fitvar_sol(sc2_idx), &
             ' sc3 = ', fitvar_sol(sc3_idx)
     END IF

     IF (wrt_to_file .AND. solfit_exval > 0) THEN
        IF ( slitcal ) THEN   ! for solar slit width calibration
           IF (which_slit == 2) THEN
              WRITE(slit_unit, '(f8.3,1p8d11.3,I6)') avgwav, vgl, vgr, hwl, hwr, &
                   shi, squ, sin, rms, solfit_exval
           ELSE
              WRITE(slit_unit, '(f8.3,1p7d11.3,I6,1pd11.3)') avgwav, hw1e, e_asym,&
                   spk, shi, squ, sin, rms, solfit_exval, &
                   varstd(hwe_idx, 2)
           END IF
        ELSE                ! for solar/rad wavelength calibration
           WRITE(slit_unit, '(f8.3,1p4d11.3, I6, 1pd11.3)') avgwav, shi, squ, sin,&
                rms, solfit_exval, varstd(shi_idx, 2)
        END IF
     END IF

     RETURN

   END SUBROUTINE cal_fit_one


   SUBROUTINE specfit ( nfitvar, fitvar, nspecpts, lowbnd, &
       uppbnd, max_itnum, rms, chisq, covar, fitspec, fitres, stderr, &
       exval, fitfunc )

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: elsunc_userdef
    USE OMSAO_parameters_module, ONLY: elsunc_np, elsunc_nw!, max_spec_pts
    USE OMSAO_variables_module, ONLY: tol,  epsrel,  epsabs,  epsx

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER,                              INTENT (IN) :: &
                                           nfitvar, nspecpts, max_itnum
    REAL (KIND=dp), DIMENSION (nfitvar),  INTENT (IN) :: lowbnd, uppbnd

    ! ==================
    ! Modified variables
    ! ==================
    REAL (KIND=dp), DIMENSION (nfitvar),  INTENT (INOUT) :: fitvar

    ! ================
    ! Output variables
    ! ================
    INTEGER,                                     INTENT (OUT) :: exval
    REAL (KIND=dp), DIMENSION (nfitvar, nfitvar),INTENT (OUT) :: covar
    REAL (KIND=dp), DIMENSION (nspecpts),  INTENT (OUT) :: fitres, fitspec
    REAL (KIND=dp),                              INTENT (OUT) :: rms, chisq
    REAL (KIND=dp), DIMENSION (nfitvar),         INTENT (OUT) :: stderr

    ! ===============
    ! Local variables
    ! ===============
    INTEGER                                      :: elbnd, i!, j
    INTEGER,        DIMENSION (elsunc_np)        :: p
    REAL (KIND=dp), DIMENSION (elsunc_nw)        :: w
    REAL (KIND=dp), DIMENSION (nfitvar)          :: blow, bupp
    REAL (KIND=dp), DIMENSION (nspecpts)         :: f
    REAL (KIND=dp), DIMENSION (nspecpts,nfitvar) :: dfda
    !REAL (KIND=dp), DIMENSION (nfitvar, nfitvar) :: correl

    EXTERNAL fitfunc

    ! ===============================================================
    ! ELBND: 0 = unconstrained
    !        1 = all variables have same lower bound
    !        else: lower and upper bounds must be supplied by the use
    ! ===============================================================
    elbnd = elsunc_userdef

    exval = 0

    p   = -1    ;  p(1)   = 0  ;    p(3) = max_itnum
    w   = -1.0  ;  w(1:4) = (/ tol,  epsrel,  epsabs,  epsx /)

    blow(1:nfitvar) = lowbnd(1:nfitvar)
    bupp(1:nfitvar) = uppbnd(1:nfitvar)
CALL elsunc ( fitvar(1:nfitvar), nfitvar, nspecpts, nspecpts,&
          fitfunc,elbnd, blow(1:nfitvar),  bupp(1:nfitvar), &
          p, w, exval, f(1:nspecpts),dfda(1:nspecpts,1:nfitvar) )
    ! ------------------------------------------------------------------
    ! Call to ELSUNC fitting function to obtain complete fitted spectrum
    ! ------------------------------------------------------------------
    CALL fitfunc ( fitvar, nfitvar, fitspec, nspecpts, 3, &
         dfda(1:nspecpts,1:nfitvar), 0)

    ! ---------------------------------------------------------------
    ! Compute fitting residual
    ! FITRES is the negative of the returned function F = Model-Data.
    ! ---------------------------------------------------------------
    fitres(1:nspecpts) = -f(1:nspecpts)

    ! Covariance matrix.
    ! ------------------
    covar(1:nfitvar,1:nfitvar) = dfda(1:nfitvar,1:nfitvar)

    ! Fitting RMS and CHI**2
    ! ----------------------
    rms  = SQRT(SUM(fitres(1:nspecpts)**2 ) / nspecpts)
    ! This gives the same CHI**2 as the NR routines
    chisq = SUM  (fitres(1:nspecpts)**2 )

    ! compute standard deviation for each variable
    DO i = 1, nfitvar
      stderr(i) = rms * SQRT(covar(i, i) * nspecpts / (nspecpts - nfitvar))
    END DO
    RETURN
  END SUBROUTINE specfit

  SUBROUTINE specfit_func_sol ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )
    !
    ! Calculates the Solar spectrum and its derivatives for ELSUNC
    ! NOTE: the variable DYDA as required here is the transpose of that
    !       rquired for the Numerical Recipes

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY : fitwavs, fitweights, &
                                       currspec,sol_wav_avg
    USE OMSAO_errstat_module

    IMPLICIT NONE

    INTEGER, INTENT (IN)    :: mdy
    INTEGER, INTENT (INOUT) :: ctrl, nfitvar, npoints
    REAL (KIND=dp), DIMENSION (nfitvar),         INTENT (INOUT) :: fitvar
    REAL (KIND=dp), DIMENSION (npoints),         INTENT (INOUT) :: ymod
    REAL (KIND=dp), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

    REAL (KIND=dp), DIMENSION (npoints) :: locwvl

    locwvl(1:npoints) = fitwavs(1:npoints)

    SELECT CASE ( ABS ( ctrl ) )
    CASE ( 1 )
      ! -----------------------------------------------------------------------
      ! Calculate the weighted difference between fitted and measured spectrum.
      ! -----------------------------------------------------------------------

   CALL spectrum_solar ( npoints, nfitvar, sol_wav_avg, &
           locwvl(1:npoints), ymod(1:npoints), fitvar(1:nfitvar) )
      ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) /fitweights(1:npoints)
    CASE ( 2 )
      ! ---------------------------------------------------------------------
      ! The following sets up ELSUNC for numerical computation of the fitting
      ! function derivative. It is faster and more flexible than the original
      ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
      ! ---------------------------------------------------------------------
      ctrl = 0
      RETURN

    CASE ( 3 )
      ! Calculate the spectrum, without weighting
      CALL spectrum_solar ( npoints, nfitvar, sol_wav_avg, &
           locwvl(1:npoints), ymod(1:npoints), fitvar(1:nfitvar))

    CASE DEFAULT
      WRITE (www_lun,'(A,I4)') &
     "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL= ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE specfit_func_sol

!  SUBROUTINE specfit_func ( vars, npars, ymod, npoints, ctrl, dyda, mdy )
!
!    !
!    !     Calculates the spectrum and its derivatives for ELSUNC
!    !
!    ! NOTE: the variable DYDA as required here is the transpose of that
!    !       rquired for the Numerical Recipes
!    !
!
!    USE OMSAO_precision_module
!    !USE OMSAO_parameters_module, ONLY : max_spec_pts
!    USE OMSAO_variables_module, ONLY : database, yn_doas, yn_smooth, &
!         rad_wav_avg, fitwavs, fitweights, currspec
!    USE OMSAO_errstat_module
!    use spectra, only: spectrum_earthshine
!
!    IMPLICIT NONE
!
!    INTEGER, INTENT (IN)    :: npars, mdy
!    INTEGER, INTENT (INOUT) :: ctrl, npoints
!    REAL (KIND=dp), DIMENSION (npars),         INTENT (IN)    :: vars
!    REAL (KIND=dp), DIMENSION (npoints),       INTENT (INOUT) :: ymod
!    REAL (KIND=dp), DIMENSION (npoints,npars), INTENT (INOUT) :: dyda
!
!    !INTEGER                                   :: i
!    !REAL (KIND=dp), DIMENSION (npars)         :: vartmp
!    REAL (KIND=dp), DIMENSION (npoints)       :: locwvl!, dyplus, dyminus
!
!    locwvl(1:npoints) = fitwavs(1:npoints)
!
!    SELECT CASE ( ABS ( ctrl ) )
!    CASE ( 1 )
!      ! Calculate the weighted difference between fitted and measured spectrum.
!      CALL spectrum_earthshine ( &
!           npoints, npars, yn_smooth, rad_wav_avg, locwvl(1:npoints), &
!           ymod(1:npoints), vars(1:npars),database, yn_doas )
!
!      ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) )&
!           / fitweights(1:npoints)
!
!    CASE ( 2 )
!      ! ---------------------------------------------------------------------
!      ! The following sets up ELSUNC for numerical computation of the fitting
!      ! function derivative. It is faster and more flexible than the original
!      ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
!      ! ---------------------------------------------------------------------
!      ctrl = 0; RETURN
!
!    CASE ( 3 )
!      ! Calculate the spectrum, without weighting
!      CALL spectrum_earthshine ( npoints, npars, yn_smooth, rad_wav_avg, &
!           locwvl(1:npoints), ymod(1:npoints), vars(1:npars), database, &
!           yn_doas )
!
!    CASE DEFAULT
!      WRITE (www_lun,'(A,I4)') &
!           "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle
!           CTRL = ", ctrl
!    END SELECT
!
!    RETURN
!  END SUBROUTINE specfit_func

END MODULE m_cal_fit_one
