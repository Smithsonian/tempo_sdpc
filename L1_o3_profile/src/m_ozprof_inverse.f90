!
module m_ozprof_inverse

  public ozprof_inverse
  private get_caloz, negativeo3_inversion

contains

  ! ***************************************************************************
  ! Author: Xiong Liu
  ! Date:   July 24, 2003
  ! Purpose: Ozone profile retrieval using PTR with GSVD-Lcurve/GCV
  !		     or Optimal Estimation
  !
  ! Modification History
  ! 1. xliu, jan 6, 2004, interface with optimal estimation
  !    a. Remove convergence determination here (done in gsvd or oem)
  !       The convergence criterion is determined just using relative
  !       chisq change.
  !    b. Remove writing fitting results for each iteration
  !    c. Modify setting up return value
  !    d. Deal with negative values in retrieved ozone (initialization
  !       is too difficult from the actual)
  ! Exit status
  ! exval = 0, not converge
  !       = 1, using absolute radiance converge in OE/PTR (delchi < epsrel)
  !       = 2, converge in LIDORT (delchi < epsrel):
  !       = 4, ozone parameters change less than epsrel!
  !       =+10, used modified a  priori (e.g., ozone hole condition)
  !       =+20, used UV2 retrievals (maybe modified a priori)
  !       = +100, negative values occur for the last iteration
  !       = 3, 1 + 2
  !       = 5, 1 + 4
  !       = 6, 2 + 4
  !       = 7, 4 + 2 + 1
  !       = -1, failure due to set_cldalb (specfit_ozprof.f90)
  !       = -2, failure due to regular matrix with Tikhonov
  !             (specfit_ozprof.f90) or making atmosphere
  !       = -3, failure due to out of bounds
  !       = -4, failure in get_raman
  !       = -5, failure in pseudo_model.f90
  !       = -6, failure due to NaN values
  ! ***************************************************************************

  SUBROUTINE ozprof_inverse (nf, varname, fitvar, fitvarap, lowbnd, upbnd,  &
       ns, np, sa, bb, nchisq, fitspec, fitres, exval)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,  ONLY: maxlay
    USE ozprof_data_module,       ONLY: ffidx=>ozfit_start_index, &
         flidx=>ozfit_end_index, ozwrtint_unit, ozwrtint, num_iter, &
         avg_kernel, contri, covar, ncovar, ozdfs, ozinfo, use_oe, nlay, &
         pfidx=>ozprof_start_index, plidx=>ozprof_end_index, ring_on_line, &
         albfidx, nfalb, use_logstate, fgasidxs, gasidxs, ngas, tracegas, &
         ozwrtfavgk, favg_kernel, radcalwrt, use_lograd, which_caloz, &
         start_layer, end_layer, atmosprof, ozwrtwf, weight_function, &
         do_simu, wrtring, wfcfidx, nfwfc, ecfrind, ecfrfind, so2zfind, &
         so2valts, use_large_so2_aperr, use_effcrs, &
         do_simu_rmring!, tf_fidx, tf_lidx, nt_fit, saa_flag, lcurve_unit,
    !do_bothstep, do_twostep,
    USE OMSAO_variables_module,   ONLY: fitvar_rad, mask_fitvar_rad, epsrel, &
         fitwavs, fitweights, maxit=>max_itnum_rad, clmspec_rad, nradpix, &
         numwin, currpix, currline, currloop, the_surfalt, band_selectors!, &
    !the_sza_atm, the_vza_atm, scnwrt, npix_fitted, fitvar_rad_init, &
    !fitvar_rad_apriori
    USE OMSAO_indices_module,     ONLY: no2_t1_idx, so2_idx, bro_idx, &
         hcho_idx, so2v_idx, o2o2_idx!, bro2_idx
    USE OMSAO_pixelcorner_module, ONLY: omi_allNSPC
    USE OMSAO_errstat_module
    use m_make_atm, only: adjust_so2vplumez
    use m_get_raman
    use m_gsvd_lcurve_gcv, only: gsvd_lcurve_gcv
    use m_oe_inversion
    use m_pseudo_model

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER, INTENT (IN)                              :: ns, nf, np
    REAL (KIND=dp), DIMENSION(nf), INTENT (IN)        :: lowbnd, upbnd
    CHARACTER (LEN=6), DIMENSION(nf), INTENT (IN)     :: varname
    REAL (KIND=dp), DIMENSION(nf), INTENT (INOUT)     :: fitvarap
    REAL (KIND=dp), DIMENSION(nf, nf), INTENT(INOUT)  :: sa
    REAL (KIND=dp), DIMENSION(nf, nf), INTENT(IN)     :: bb

    ! ==================
    ! Modified variables
    ! ==================
    REAL (KIND=dp), DIMENSION (nf), INTENT (INOUT) :: fitvar

    ! ================
    ! Output variables
    ! ================
    INTEGER,        INTENT (OUT)                   :: exval
    REAL (KIND=dp), INTENT (OUT)                   :: nchisq
    REAL (KIND=dp), DIMENSION (ns), INTENT (OUT)   :: fitspec, fitres

    ! ================
    ! Local Variables
    ! ================
    LOGICAL :: refl_only, proceed, do_sa_diagonal, conv, varconv, negval, &
         last_iter, so2aperr_update, decorrelate, use_uv2init, correct_merr
    INTEGER :: i, j, k, errstat, uv2fy, uv2ly, nuv2, fidx, lidx, &
         cmerr_niter, uv12_retflg!, so2zfidx
    REAL (KIND=dp) :: ochisq, lchisq, delchi, nradrms, oradrms, &
         readout_noise!, avgres, sumdfs, sumo, sumn
    REAL (KIND=dp), DIMENSION(nf) :: delta_x, xold, xap, std, nstd, aperr, &
         fitvarap0
    REAL (KIND=dp), DIMENSION(ns) :: gspec, sig, fitspec1, fitres1, &
         gspec_new, simrad!, selidx
    REAL (KIND=dp), DIMENSION(ns, nf)   :: dyda
    !REAL (KIND=dp), DIMENSION(ns+1)     :: sig1
    !REAL (KIND=dp), DIMENSION(ns+1, nf) :: dyda1
    !REAL (KIND=dp), DIMENSION(nf, ns+1) :: contri1
    !REAL (KIND=dp), DIMENSION(nf, nf)   :: correl
    REAL (KIND=dp), DIMENSION(maxlay)   :: tozprof
    real (kind=dp), dimension(nf,nf) :: tmp_covar, tmp_ncovar, tmp_avg_kernel
    real (kind=dp), dimension(nf,ns) :: tmp_contri

    LOGICAL, SAVE                       :: first = .TRUE.
    INTEGER, SAVE                       :: no2fidx, so2vfidx, so2fidx, &
         brofidx, hchofidx, o4fidx

    ! Variables for adjusting measurement errors based on fitting residuals
    LOGICAL                           :: adjust_merr
    INTEGER, PARAMETER                :: nreg = 5
    INTEGER, DIMENSION(nreg)          :: regfidxs, reglidxs, regnpts
    REAL (KIND=dp), DIMENSION(nreg)   :: reg_rms, reg_res
    REAL (KIND=dp), DIMENSION(0:nreg) :: reg_waves = &
         (/260.0, 290.0, 300.0, 310.0, 325.0, 350./)

    ! For detect a NaN
    INTEGER, PARAMETER :: DBPRECISION = SELECTED_INT_KIND(PRECISION(1.0d0))
    INTEGER (DBPRECISION), PARAMETER :: NAN = Z"7FF8000000000000"
    INTEGER, PARAMETER :: DPSB = BIT_SIZE(NAN) - 1

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=14), PARAMETER :: modulename = 'ozprof_inverse'

    IF (first) THEN
      DO i = 1, ngas
        IF (gasidxs(i) == no2_t1_idx) no2fidx  = fgasidxs(i)
        IF (gasidxs(i) == so2v_idx)   so2vfidx = fgasidxs(i)
        IF (gasidxs(i) == so2_idx)    so2fidx  = fgasidxs(i)
        IF (gasidxs(i) == bro_idx)    brofidx  = fgasidxs(i)
        IF (gasidxs(i) == hcho_idx)   hchofidx = fgasidxs(i)
        IF (gasidxs(i) == o2o2_idx)   o4fidx   = fgasidxs(i)
      ENDDO

      first = .FALSE.
    ENDIF

    use_uv2init = .FALSE.
    IF (numwin >= 2 ) THEN
      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        IF (band_selectors(i) == 2) THEN
          uv2fy = fidx; uv2ly = lidx; nuv2 = nradpix(i)
          use_uv2init = .TRUE.; EXIT
        ENDIF
        fidx = lidx + 1
      ENDDO
    ENDIF
    uv12_retflg = 0 ! 0: uv1+uv2 1: uv1+uv2+modified a priori 2: uv2 only


    num_iter  = 0
    refl_only = .FALSE.       ! need both radiances and weighting function
    proceed   = .TRUE.;    conv      = .FALSE.; varconv = .FALSE.
    do_sa_diagonal = .FALSE.; negval = .FALSE.
    sig = 1.0 ! measurement error included in dyda and gspec (normalized to 1)
    exval = 0
    covar = 0.0
    ncovar = 0.0 ! covariance matrix
    fitvarap0 = fitvarap
    decorrelate = .FALSE.
    correct_merr = .FALSE.
    cmerr_niter = 0 !maxit + 1
    readout_noise=1.0

    ! After retrievals are done with currently assumed measurement errors,
    ! perform retrievals again with adjusted measurement errors based on
    ! fitting residuals in several different spectral regions
    adjust_merr = .FALSE.
    IF (adjust_merr) THEN
      correct_merr = .FALSE.
      fidx = 1;     i = fidx
      DO j = 1, nreg
        DO WHILE (fitwavs(i) < reg_waves(j) .AND. i <= ns)
          i = i + 1
        ENDDO
        lidx = i - 1
        regfidxs(j) = fidx; reglidxs(j) = lidx; regnpts(j) = lidx - fidx + 1
        fidx = lidx + 1
      ENDDO
    ENDIF

    ! Calculate ring spectrum
    IF (ring_on_line .AND. .NOT. (do_simu .AND. radcalwrt .AND. &
         .NOT. wrtring .AND. .NOT. do_simu_rmring) ) THEN
      CALL GET_RAMAN(nlay, fitvar_rad(pfidx:plidx), errstat)
      IF (errstat == pge_errstat_error) THEN
        exval = -4; RETURN
      ENDIF
    ENDIF

    IF (ozwrtint) WRITE(ozwrtint_unit, '(A,I5,A10,I5, A10, I5)')  'Line = ', &
         currline, ' XPix = ', currpix, ' Loop = ', currloop

    ochisq = 10.0D20;  oradrms = 100.0
    inverse: DO WHILE (proceed)

      ! compute radiance, spectrum, and weighting function
      ! xliu: 0/1/28/1010: add use_effcrs here because use_hres always
      !  needs weighting function for correction
      IF (radcalwrt .AND. do_simu .AND. .NOT. do_simu_rmring .AND. &
           use_effcrs) THEN
        refl_only = .TRUE.
      ENDIF
      IF (num_iter == cmerr_niter .AND. correct_merr) fitweights(1:ns) = &
           fitweights(1:ns) / SQRT(1.0d0 * omi_allNSPC(currline)) * &
           readout_noise

      CALL pseudo_model(num_iter, refl_only, ns, nf, fitvar, fitvarap, &
           dyda, gspec, fitres, fitspec, nchisq, nradrms, errstat)

      IF (errstat == pge_errstat_error) THEN
        proceed = .FALSE.
        exval = -5
        CYCLE
      ENDIF
      IF ( radcalwrt .AND. do_simu .AND. .NOT. do_simu_rmring) THEN
        exval = 0
        RETURN
      ENDIF

      xold = fitvar
      IF (use_logstate) THEN
        xold(ffidx:flidx) = LOG(xold(ffidx:flidx))
        DO i = ffidx, flidx
          dyda(:, i) = dyda(:, i) * fitvar(i)
        ENDDO
      ENDIF
      xap = fitvarap
      IF (use_logstate) xap(ffidx:flidx) = LOG(xap(ffidx:flidx))

      delchi = ABS(nradrms - oradrms) / SQRT(oradrms)
      ! check for NAN
      IF (IEOR(IBCLR(TRANSFER(delchi, NAN), DPSB), NAN) == 0) THEN
        proceed = .FALSE.; exval = -6; CYCLE
      ENDIF

      ochisq = nchisq; oradrms = nradrms
      IF (delchi < epsrel) THEN     ! converge, exit
        exval = 1
        !proceed = .FALSE.; CYCLE
      ENDIF

      num_iter = num_iter + 1
      so2aperr_update = .FALSE.

      DO
        IF (use_oe) THEN  ! use optimal estimation
          last_iter = .TRUE.
          !          IF (.NOT. do_twostep) THEN
          ! FIXME - masking array temporaries
          CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, &
               epsrel, last_iter, num_iter, ns, nf, gspec, sig, dyda, &
               xap, xold, sa, varname, ffidx, flidx, delta_x, &
                      ! covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, &
                      ! avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, &
               tmp_covar(:,:), tmp_ncovar(:,:), conv, &
               tmp_avg_kernel(:,:), tmp_contri(:,:), ozdfs, &
               ozinfo, lchisq, gspec_new)
          covar(1:nf, 1:nf)=tmp_covar
          ncovar(1:nf, 1:nf)=tmp_ncovar
          avg_kernel(1:nf, 1:nf)=tmp_avg_kernel
          contri(1:nf, 1:ns)=tmp_contri

          ! ELSE
          ! CALL twostep_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, &
          !      epsrel, last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, &
          !      xold, lowbnd, upbnd, sa, varname, ffidx, flidx, delta_x, &
          !      covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, &
          !      avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, &
          !      ozinfo, lchisq)
          ! ENDIF
        ELSE    ! use Phillips-Tikhonov regularization
          CALL gsvd_lcurve_gcv (ozwrtint, ozwrtint_unit, epsrel, num_iter, &
               ns, nf, np, bb(1:np, 1:nf), gspec, dyda, fitvarap, xold, &
               varname, delta_x, covar(1:nf, 1:nf), conv, &
               avg_kernel(1:nf,1:nf), contri(1:nf,1:ns), ozdfs, ozinfo, lchisq)
        ENDIF

        IF ( radcalwrt .AND. do_simu .AND. do_simu_rmring) THEN
          ! Remove Ring effect and I/F wavelength shift
          simrad(1:ns) = fitspec(1:ns) - fitres(1:ns)
          !DO i = 1, ns
          !  WRITE(90, '(4D14.6)') fitwavs(i), fitspec(i), simrad(i), fitres(i)
          !ENDDO

          ! Ring effect, assume ring effect last variable, one single band
          gspec(1:ns) = gspec(1:ns) - dyda(1:ns, nf) * delta_x(nf)

          ! I/F wavelength shift (should not since wavelength scale is not
          ! changed)
          !gspec(1:ns) = gspec(1:ns) - dyda(1:ns, nf-1) * delta_x(nf-1)
          fitres(1:ns) = gspec(1:ns) * fitweights(1:ns)
          fitspec(1:ns) = fitres(1:ns) + simrad(1:ns)

          exval = 0; RETURN
        ENDIF

        fitvar = delta_x + xold
        ! check for NAN
        IF (IEOR(IBCLR(TRANSFER(lchisq, NAN), DPSB), NAN) == 0) THEN
          proceed = .FALSE.; exval = -6; EXIT
        ENDIF

        !DO i = 1, nf
        !   DO j = 1, nf
        !      correl(i, j) = covar(i, j) / SQRT(covar(i, i) * covar(j, j))
        !   ENDDO
        !ENDDO
        !WRITE(www_lun, '(10F10.3)') correl(so2zfind, :)
        !WRITE(www_lun, *)

        ! Special treatment for SO2
        !IF (do_twostep .OR. (so2vfidx <= 0 .AND. so2fidx <= 0) .OR. &
        !use_large_so2_aperr) EXIT
        IF ((so2vfidx <= 0 .AND. so2fidx <= 0) .OR. use_large_so2_aperr) EXIT
        so2aperr_update = .FALSE.
        IF ( so2vfidx > 0) THEN
          IF ( ABS(fitvar(so2vfidx)) > SQRT(sa(so2vfidx, so2vfidx)) * 0.5) THEN
            sa(so2vfidx, so2vfidx) = 4.0 * fitvar(so2vfidx) ** 2.0
            so2aperr_update  = .TRUE.
          ENDIF
        ENDIF

        IF ( so2fidx > 0 ) THEN
          IF ( ABS(fitvar(so2fidx)) > SQRT(sa(so2fidx, so2fidx)) * 0.5) THEN
            sa(so2fidx, so2fidx) = 4.0 * fitvar(so2fidx) ** 2.0
            so2aperr_update = .TRUE.
          ENDIF
        ENDIF
        IF (.NOT. so2aperr_update ) EXIT
      ENDDO

      IF (use_logstate) THEN
        fitvar(ffidx:flidx) = EXP(fitvar(ffidx:flidx))
        xold(ffidx:flidx) = EXP(xold(ffidx:flidx))
        delta_x(ffidx:flidx) = fitvar(ffidx:flidx) - xold(ffidx:flidx)
      ENDIF

      !WRITE(www_lun, '(I3,1X,A6,3d14.6)') ((i, varname(i), xold(i), &
      !     delta_x(i), fitvar(i)), i=1, nf)
      !WRITE(www_lun, *) SUM(xold(ffidx:flidx)), SUM(fitvar(ffidx:flidx))

      uv12_retflg = 0
      IF (ANY (fitvar(ffidx:flidx) <= lowbnd(ffidx:flidx)) .OR. &
           ANY (fitvar(ffidx:flidx) >= upbnd(ffidx:flidx))) THEN

        IF (use_uv2init) THEN
          CALL negativeo3_inversion (uv2fy, uv2ly, nuv2, do_sa_diagonal,  &
               ozwrtint, ozwrtint_unit, epsrel, last_iter, num_iter, ns, nf, &
               gspec, sig, dyda, xap, xold, lowbnd, upbnd, sa, varname, &
               ffidx, flidx, delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), &
               conv, avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, &
               ozinfo, lchisq, gspec_new, uv12_retflg)

          fitvar = delta_x + xold
          fitvarap(ffidx:flidx) = xap(ffidx:flidx)
        ENDIF
      END IF

      IF (ANY (fitvar(ffidx:flidx) <= lowbnd(ffidx:flidx)) .OR. &
           ANY (fitvar(ffidx:flidx) >= upbnd(ffidx:flidx))) THEN

        IF (so2aperr_update) THEN
          ! Redo the retrieval using original state with updated SO2 fields
          IF (so2vfidx > 0) THEN
            xold(so2vfidx) = fitvar(so2vfidx); fitvar = xold
            fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar
            fitvarap(so2vfidx) = fitvar(so2vfidx)
          ENDIF
        ELSE
          WRITE(www_lun, *) modulename, &
               ': Retrieved ozone values out of bounds!!!'
          proceed = .FALSE.;     exval = -3; CYCLE
        ENDIF
      ELSE

        ! Unnecessary here, it is the a priori error that matters
        !IF (so2aperr_update) THEN
        !   ! Update the a priori for the next iteration
        !   fitvarap(so2vfidx) = fitvar(so2vfidx)
        !ENDIF

        ! Update a priori value for ozone at every iterations
        !IF (num_iter >= 2) xap(ffidx:flidx) = fitvar(ffidx:flidx)

        ! update the uncondensed fitting variables
        fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar
      END IF

      IF (ALL(ABS(delta_x(ffidx:flidx) / fitvar(ffidx:flidx)) <= epsrel)) THEN
        varconv = .TRUE.
      ENDIF

      IF (num_iter >= maxit .AND. .NOT. conv)  THEN
        proceed = .FALSE.;      exval = 0; CYCLE
      END IF

      IF ((conv .OR. varconv)) THEN
        IF (conv) exval = exval + 2
        IF (varconv) exval = exval + 4
        proceed = .FALSE.; CYCLE
      ENDIF

      IF (ring_on_line .AND. ABS(SUM(delta_x(ffidx:flidx))) > 20.0D0) THEN
        CALL GET_RAMAN(nlay, fitvar_rad(pfidx:plidx), errstat)
        IF (errstat == pge_errstat_error) THEN
          exval = -4; proceed = .FALSE.
        ENDIF
      ENDIF

      ! Check if need to update the SO2V profile shape
      IF (so2zfind > 0 ) THEN
        IF (fitvar(so2zfind) /= xold(so2zfind)) THEN
          IF (fitvar(so2zfind) > the_surfalt) THEN
            so2valts(0) = fitvar(so2zfind)
            so2valts(-1) = so2valts(0) - 1.0
            IF (so2valts(-1) < the_surfalt) so2valts(-1) = &
                 (the_surfalt + so2valts(0)) / 2.0
            so2valts(1)  = so2valts(0) + 1.0
            ! Update the a priori value for SO2 plume height for next iteration
            ! (since the apriori value is an arbitrary guess)
            fitvarap(so2zfind) = fitvar(so2zfind)
          ELSE
            fitvar(so2zfind) = xold(so2zfind)
            fitvar_rad(mask_fitvar_rad(so2zfind)) = fitvar(so2zfind)
          ENDIF
          CALL ADJUST_SO2VPLUMEZ(errstat)
        ENDIF
      ENDIF
    END DO inverse

    DO i = 1, nf
      aperr(i) = SQRT(sa(i, i))
    ENDDO

    DO k = 1, ngas
      i = fgasidxs(k)
      IF (i > 0) THEN
        tracegas(k, 9) = avg_kernel(i, i)
        ! consider interference from others
        tracegas(k, 10) = SUM(avg_kernel(i, 1:nf) * aperr(1:nf) / aperr(i) )
      ENDIF
    ENDDO
    IF (ozwrtfavgk) favg_kernel(1:nf, 1:nf) = avg_kernel(1:nf, 1:nf)

    IF (ozwrtwf .AND. exval >= 0) THEN
      DO i = 1, nf
        weight_function(1:ns, i) = dyda(1:ns, i) * fitweights(1:ns)
      ENDDO
    ENDIF

    ! xliu: 03/19/2010
    ! Averaging kernels have already been calculated for each iteration
    ! No need to call oe_inversion unless:
    ! a. Decorrelate ozone with other varaibles
    ! b. Recalculate averaging kernels with corrected measurement errors
    ! 05/26/2010
    ! c. Adjust measurement error to reflect actual fitting residuals
    IF ( exval >= 0 .AND. (decorrelate .OR. &
         (correct_merr .AND. cmerr_niter == maxit + 1) .OR. adjust_merr) ) THEN

      IF (correct_merr .AND. cmerr_niter == maxit + 1) THEN
        gspec(1:ns) = gspec(1:ns) * &
             SQRT(1.0d0 * omi_allNSPC(currline)) / readout_noise
        DO i = 1, nf
          dyda(1:ns, i) = dyda(1:ns, i) * &
               SQRT(1.0d0 * omi_allNSPC(currline)) / readout_noise
        ENDDO
        fitweights(1:ns) = fitweights(1:ns) / &
             SQRT(1.0d0 * omi_allNSPC(currline)) * readout_noise
        nchisq = nchisq  * omi_allNSPC(currline) / &
             readout_noise / readout_noise
      ENDIF

      IF (adjust_merr) THEN
        DO i = 1, nreg
          fidx = regfidxs(i); lidx = reglidxs(i)
          reg_rms(i) = SQRT(SUM((fitres(fidx:lidx)/fitweights(fidx:lidx))**2) &
               / regnpts(i))
          reg_res(i) = SQRT(SUM(fitres(fidx:lidx)**2) / regnpts(i))

          fitweights(fidx:lidx) = fitweights(fidx:lidx) * reg_rms(i)
          gspec(fidx:lidx) = gspec(fidx:lidx) / reg_rms(i)
          dyda(fidx:lidx, 1:nf) = dyda(fidx:lidx, 1:nf) / reg_rms(i)
          nchisq  = SUM((fitres(1:ns) / fitweights(1:ns))**2.0)
          !WRITE(*, '(2F10.4,3I5,2D14.5)') reg_waves(i-1), reg_waves(i), &
          !     fidx, lidx, regnpts(i), reg_rms(i), reg_res(i)
        ENDDO
      ENDIF

      ! save standard deviations for those variables
      DO i = 1, nf
        std(i) = covar(i, i); nstd(i) = ncovar(i, i)
      ENDDO

      IF ( decorrelate )  THEN
        DO i = 1, nf
          IF (i < ffidx .OR. i > flidx ) THEN
            dyda(:, i)  = 0.0D0
            fitvarap(i) = fitvar(i)
          ENDIF
        ENDDO
      ENDIF
      !xold = fitvar
      xap = fitvarap

      IF (use_logstate) THEN
        xold(ffidx:flidx) = LOG(xold(ffidx:flidx))
        xap(ffidx:flidx)  = LOG(xap(ffidx:flidx))
      ENDIF

      IF (use_oe) THEN  ! use optimal estimation
        last_iter = .TRUE.
        !IF (.NOT. do_twostep) THEN
        CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel, &
             last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa, &
             varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), &
             ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf), &
             contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec_new)
        !ELSE
        !  CALL twostep_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, &
        !       epsrel, last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, &
        !       xold, lowbnd, upbnd, sa, varname, ffidx, flidx, delta_x, &
        !       covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, &
        !       avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, ozinfo, &
        !       lchisq)
        !ENDIF

        DO i = 1, nf
          contri(i, 1:ns) = contri(i, 1:ns) / fitweights(1:ns)
        ENDDO

      ELSE    ! use Phillips-Tikhonov regularization
        CALL gsvd_lcurve_gcv (ozwrtint, ozwrtint_unit, epsrel, num_iter, &
             ns, nf, np, bb(1:np, 1:nf), gspec, dyda, fitvarap, xold, varname, &
             delta_x, covar(1:nf, 1:nf), conv, avg_kernel(1:nf,1:nf), &
             contri(1:nf,1:ns), ozdfs, ozinfo, lchisq)
      END IF
      fitvar = xold  + delta_x

      !WRITE(www_lun, '(I3,1X,A6,d14.6,f10.4,d14.6)') &
      ! ((i, varname(i), xold(i), delta_x(i) / xold(i)*100.0, fitvar(i)), &
      ! i=1, nf)

      IF (use_logstate) THEN
        fitvar(ffidx:flidx) = EXP(fitvar(ffidx:flidx))
        xold(ffidx:flidx) = EXP(xold(ffidx:flidx))
        delta_x(ffidx:flidx) = fitvar(ffidx:flidx) - xold(ffidx:flidx)
      ENDIF

      IF (ANY (fitvar(ffidx:flidx) <= lowbnd(ffidx:flidx)) .OR. &
           ANY (fitvar(ffidx:flidx) >= upbnd(ffidx:flidx))) THEN
        WRITE(www_lun, *) modulename, &
             ': Retrieved ozone values out of bounds!!!'
        exval = -3
      ELSE
        ! update the uncondensed fitting variables
        fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar
      END IF

      DO i = 1, nf
        IF (i < ffidx .OR. i > flidx) THEN
          covar(i, i) = std(i); ncovar(i, i) = nstd(i)
        ENDIF
      ENDDO
    ENDIF

    ! This is incorrect, and need to be changed
    IF  (use_logstate) THEN
      DO i = ffidx, flidx
        covar(i, i) = covar(i, i) * fitvar(i) ** 2.0
        ncovar(i, i) = ncovar(i, i) * fitvar(i) ** 2.0
        DO j = i+1, flidx
          covar(i, j) =  covar(i, j) * fitvar(i) * fitvar(j)
          covar(j, i) = covar(i,j)
          ncovar(i, j) = ncovar(i, j) * fitvar(i) * fitvar(j)
          ncovar(j, i) = ncovar(i,j)
        ENDDO
      ENDDO
    ENDIF

    ! For exval  > 1, ideally still need to calculate final spectrum
    ! but this step is unnecssary, since retrieval is already done.
    ! On the other hand, it can save computation by one iteration (significant)
    ! for exval == 1, final spectrum already calculated
    ! IF (exval > 1) THEN
    !   refl_only = .TRUE.
    !   CALL pseudo_model(num_iter, refl_only, ns, nf, fitvar, fitvarap, &
    !        dyda, gspec, fitres, fitspec, nchisq, nradrms, errstat)
    !   delchi = ABS(nradrms - oradrms) / oradrms   ! converge, exit
    !   IF (delchi < epsrel)  exval = exval + 1
    !END IF

    IF (exval >= 0 .AND. radcalwrt) THEN
      refl_only = .TRUE.
      xold = fitvar

      !Save retrievals
      !Use a priori albedo and ozone, the other are the same
      IF (which_caloz == 1) THEN
        fitvar(ffidx:flidx) = fitvarap(ffidx:flidx)
      ELSE
        tozprof(1:nlay)                = fitvar_rad(pfidx:plidx)
        tozprof(start_layer:end_layer) = fitvarap(ffidx:flidx)
        CALL get_caloz(nlay, atmosprof(1, 0:nlay), tozprof(1:nlay))
        fitvar(ffidx:flidx) = tozprof(start_layer:end_layer)
      ENDIF

      IF (nfalb > 0) fitvar(albfidx:albfidx+nfalb-1) = &
           fitvarap0(albfidx:albfidx+nfalb-1)
      IF (nfwfc > 0) fitvar(wfcfidx:wfcfidx+nfwfc-1) = &
           fitvarap0(wfcfidx:wfcfidx+nfwfc-1)
      IF (ecfrfind > 0) fitvar(ecfrind) = fitvarap0(ecfrind)
      fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar

      CALL pseudo_model(num_iter, refl_only, ns, nf, fitvar, fitvarap, &
           dyda, gspec, fitres1, fitspec1, nchisq, nradrms, errstat)

      ! Restore the retrieved variables
      fitvar = xold
      fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar

      IF (.NOT. use_lograd) THEN
        clmspec_rad(1:ns)  = fitspec1(1:ns) - fitres1(1:ns)
      ELSE
        clmspec_rad(1:ns)  = EXP(fitspec1(1:ns) - fitres1(1:ns))
      ENDIF
    ENDIF
    IF (exval > 0) THEN
      IF (uv12_retflg == 1) exval = exval + 10
      IF (uv12_retflg == 2) exval = exval + 20
      IF ( ANY(fitvar(ffidx:flidx) <= 0.0)) negval = .TRUE.
      IF (negval) exval = exval + 100
    ENDIF

    RETURN

  END SUBROUTINE ozprof_inverse


  SUBROUTINE get_caloz (nl, pres, ozprof)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY : mflay, rearth
    USE ozprof_data_module,      ONLY : caloz_fname, profunit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline, reverse

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER, INTENT (IN)                            :: nl
    REAL (KIND=dp), DIMENSION(0:nl), INTENT (INOUT) :: pres
    REAL (KIND=dp), DIMENSION(nl),   INTENT (INOUT) :: ozprof

    ! ==================
    ! Logical variables
    ! ==================
    INTEGER :: i, nz, oztyp, sidx, eidx, nlay, errstat
    REAL (KIND=dp)                                  :: alt, gcorr
    REAL (KIND=dp), DIMENSION(0:mflay)              :: zs
    REAL (KIND=dp), DIMENSION(0:nl)                 :: cozprof

    LOGICAL                                         :: first = .true.
    REAL (KIND=dp), DIMENSION(0:mflay), SAVE        :: ps, ozs, cozs


    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=9), PARAMETER :: modulename = 'get_caloz'

    IF ( first ) THEN
      OPEN(profunit, FILE=TRIM(ADJUSTL(caloz_fname)), STATUS='old')
      READ(profunit, *) nz, oztyp
      IF (nz > mflay) THEN
        WRITE(www_lun, *) modulename, ': Need to increase mflay!!!'; STOP
      ENDIF
      READ (profunit, *)
      IF (oztyp <= 2) THEN
        READ (profunit, *) ozs(1:nz)
      ELSE
        READ (profunit, *) ozs(0:nz)
      ENDIF

      IF ( oztyp == 1 .AND. nz /= nl ) THEN
        WRITE(www_lun, *) modulename, ': Number of layers are inconsistent!!!'
        STOP 1
      ELSE
        READ (profunit, *)
        READ (profunit, *) ps(0:nz)  ! mb
      ENDIF
      CLOSE (profunit)

      IF (oztyp >= 2) THEN
        ! make profiles top-down
        IF ( ps(0) > ps(1) ) THEN
          CALl REVERSE(ozs(1:nz), nz)
          CALL REVERSE(ps(0:nz), nz+1)
        ENDIF

        ! Convert ozone from mixing ratio to partial column ozone
        IF (oztyp == 2) THEN        ! Get cumulative ozone profile
          cozs(0) = 0.0
          DO i = 1, nz
            cozs(i) = cozs(i-1) + ozs(i)
          ENDDO
        ELSE IF (oztyp == 3) THEN
          ! Integrate from ppbv to DU (with gravity correction)
          zs = - 16.0 * LOG10( ps / 1013.25)
          cozs(0) = 0.0
          DO i = 1, nz             ! 2533.12 = 1.25 / 0.5 * 1013.25
            alt     = ( zs(i-1) + zs(i) ) / 2.0
            gcorr   = ( rearth / (rearth + alt) ) ** 2.0 * 2533.125
                                               !* 1.25 / 0.5 * 1013.25
            cozs(i) = cozs(i-1) + ( ozs(i-1) + ozs(i) ) * &
                 (ps(i) - ps(i-1)) / gcorr
          ENDDO
        ENDIF
      ENDIF

      first = .FALSE.
    ENDIF

    IF (oztyp == 1 ) THEN
      ozprof(1:nl) = ozs(1:nz)
    ELSE
      ! Interpolate to the retrieval grid
      sidx = MINVAL(MINLOC(pres, MASK=(pres >= ps(0 )))) - 1
      eidx = MINVAL(MAXLOC(pres, MASK=(pres <= ps(nz)))) - 1
      nlay = eidx - sidx

      ps = LOG(ps); pres = LOG(pres)
      CALL BSPLINE(ps, cozs, nz+1, pres(sidx:eidx), cozprof(sidx:eidx), &
           nlay+1, errstat)
      ozprof(sidx+1:eidx) = cozprof(sidx+1:eidx) - cozprof(sidx:eidx-1)
      ps = EXP(ps); pres = EXP(pres)

    ENDIF

    RETURN

  END SUBROUTINE get_caloz




  SUBROUTINE negativeo3_inversion (uv2fy, uv2ly, nuv2, do_sa_diagonal,  &
       ozwrtint, ozwrtint_unit, epsrel, last_iter, num_iter, ns, nf, gspec, &
       sig, dyda, xap, xold, lowbnd, upbnd, sa, varname, ffidx, flidx,      &
       delta_x, covar, ncovar, conv, avg_kernel, contri, ozdfs, ozinfo,     &
       lchisq, gspec_new, uv12_retflg)

    USE OMSAO_precision_module
    !USE OMSAO_parameters_module,  ONLY: maxlay, max_spec_pts
    USE ozprof_data_module,       ONLY: atmosprof, nlay
    USE OMSAO_variables_module,   ONLY: the_lat, the_month, the_day!, &
    !mask_fitvar_rad, fitwavs, fitvar_rad
    USE OMSAO_errstat_module
    use prepare_atmosphere, only: get_tomsv8_clima
    use m_oe_inversion

    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT (IN)                             :: uv2fy, uv2ly, nuv2, &
         ns, nf, ozwrtint_unit, num_iter, ffidx, flidx
    INTEGER, INTENT (OUT)                            :: uv12_retflg
    REAL (KIND=dp), INTENT (IN)                      :: epsrel
    REAL (KIND=dp), INTENT (OUT)                     :: ozdfs, ozinfo, lchisq
    REAL (KIND=dp), DIMENSION(nf), INTENT (IN)       :: lowbnd, upbnd, xold
    REAL (KIND=dp), DIMENSION(nf), INTENT (OUT)      :: delta_x
    REAL (KIND=dp), DIMENSION(ns), INTENT (IN)       :: gspec, sig
    CHARACTER (LEN=6), DIMENSION(nf), INTENT (IN)    :: varname
    REAL (KIND=dp), DIMENSION(ns, nf), INTENT (INOUT):: dyda
    REAL (KIND=dp), DIMENSION(nf, ns), INTENT (OUT)  :: contri
    REAL (KIND=dp), DIMENSION(nf, nf), INTENT(INOUT) :: sa
    REAL (KIND=dp), DIMENSION(nf, nf), INTENT(OUT)   :: covar, ncovar, &
         avg_kernel
    REAL (KIND=dp), DIMENSION(ns), INTENT (OUT)      :: gspec_new
    REAL (KIND=dp), DIMENSION(nf), INTENT (OUT)      :: xap
    LOGICAL, INTENT(IN)                              :: do_sa_diagonal, &
         ozwrtint, last_iter
    LOGICAL, INTENT(OUT)                             :: conv

    ! =======================
    ! Local variables
    ! =======================
    INTEGER                         :: errstat
    REAL (KIND=dp)                  :: newtoz, oldtoz, aptoz
    REAL (KIND=dp), DIMENSION(nlay) :: newapoz, apoz
    REAL (KIND=dp), DIMENSION(nf) :: delta_x1, tmp_fitvar

    ! First try inversion with UV2 only
    uv12_retflg = 2
    CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel, &
         last_iter, num_iter, nuv2, nf, gspec(uv2fy:uv2ly), sig(uv2fy:uv2ly), &
         dyda(uv2fy:uv2ly, :), xap, xold, sa, varname, ffidx, flidx, delta_x, &
         covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf), &
         contri(1:nf, uv2fy:uv2ly), ozdfs, ozinfo, lchisq, &
         gspec_new(uv2fy:uv2ly))

    aptoz  = SUM(xap(ffidx:flidx))
    apoz   = xap(ffidx:flidx)
    oldtoz = SUM(xold(ffidx:flidx))
    newtoz = SUM(delta_x(ffidx:flidx)) + oldtoz

    !print *, 'Inside'
    !WRITE(*, '(4F8.3)') aptoz, oldtoz, newtoz, SUM(delta_x(ffidx:flidx))

    ! If UV2 total ozone is larger than a priori by 50 DU, then
    ! use total ozone dependent TOMS V8 climatology to replace xap
    IF (ABS(newtoz - aptoz) >= 50.0 .AND. newtoz >= 100. .AND. &
         newtoz <= 600.) THEN
      CALL get_tomsv8_clima(the_month, the_day, the_lat, newtoz, &
           nlay, atmosprof(1, 0:nlay), apoz(1:nlay), newapoz(1:nlay), errstat)

      IF (errstat /= pge_errstat_error) THEN
        xap(ffidx:flidx) = newapoz(1:nlay)

        ! Try rertievals with updated a priori and with both retrievals,
        ! accept results if it is succesful, otherwise, use UV2 retrievals
        CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel, &
             last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa,  &
             varname, ffidx, flidx, delta_x1, covar(1:nf, 1:nf), &
             ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf),  &
             contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec_new)

        tmp_fitvar = xold + delta_x1

        IF (ALL (tmp_fitvar(ffidx:flidx) > lowbnd(ffidx:flidx)) .AND. &
             ALL (tmp_fitvar(ffidx:flidx) < upbnd(ffidx:flidx))) THEN
          delta_x = delta_x1
          uv12_retflg = 1
          !print *, 'Use modified UV1 + UV2'
          !WRITE(*, '(4F8.3)') aptoz, oldtoz, &
          !  oldtoz+SUM(delta_x(ffidx:flidx)), SUM(delta_x(ffidx:flidx)
        ELSE
          CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel, &
               last_iter, num_iter, nuv2, nf, gspec(uv2fy:uv2ly), &
               sig(uv2fy:uv2ly), dyda(uv2fy:uv2ly, :), xap, xold, sa, &
               varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), &
               ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf), &
               contri(1:nf, uv2fy:uv2ly), ozdfs, ozinfo, lchisq, &
               gspec_new(uv2fy:uv2ly))
          !print *, 'Still use UV2'
        ENDIF
      ENDIF
    ENDIF

    RETURN
  END SUBROUTINE negativeo3_inversion

end module m_ozprof_inverse
