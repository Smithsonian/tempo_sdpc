MODULE OMSAO_solar_wavcal_module
  use optimizer_interface_module
  use errormodule

  IMPLICIT NONE

  PRIVATE solar_fit

CONTAINS

  SUBROUTINE xtrack_solar_calibration_loop ( first_pix, last_pix, errstat )

    USE OMSAO_precision_module
    USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
    USE omi_pge_fitting_aux, ONLY: omi_adjust_irradiance_data
    USE OMSAO_omidata_module, ONLY: &
      omi_cross_track_skippix, &
      omi_solcal_chisq, omi_solcal_pars, omi_solcal_xflag, &
      omi_solcal_itnum, omi_irradiance_wght, omi_irradiance_ccdpix
    USE OMSAO_indices_module, ONLY: wvl_idx, sig_idx, spc_idx, ccd_idx, &
      max_calfit_idx, shi_idx, squ_idx, solcal_idx
    USE OMSAO_parameters_module, ONLY: r8_missval, i2_missval, i4_missval, MAX_STR_LEN
    USE OMSAO_variables_module,  ONLY: verb_thresh_lev, Slit_Half_Width_1e, &
      Slit_Asym_Factor, curr_sol_spec, fitvar_cal, fitvar_cal_saved,  &
      fitvar_sol_init, ctrl_n_fitres_loop, ctrl_fitres_range, &
      curr_xtrack_pixnum
    USE OMSAO_errstat_module
    use irradiance_data, only : Irr_Data

    IMPLICIT NONE
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i2)              :: solcal_itnum
    INTEGER   (KIND=i4)              :: locerrstat, ipix, solcal_exval, n_irradwvl
    CHARACTER (LEN=MAX_STR_LEN)         :: addmsg
    REAL      (KIND=r8)              :: chisquav, curr_sol_wav_avg
    LOGICAL                          :: yn_skip_pix, is_bad_pixel

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=29), PARAMETER :: modulename = 'xtrack_solar_calibration_loop'

    if (errstat < 0) return

    omi_solcal_chisq = r8_missval

    fitvar_cal_saved(1:max_calfit_idx) = fitvar_sol_init(1:max_calfit_idx)

    ! ---------------------------------------------------------------
    ! Loop for solar wavelength calibration and slit function fitting
    ! ---------------------------------------------------------------

    XtrackSolCal: DO ipix = first_pix, last_pix

      locerrstat = pge_errstat_ok

      curr_xtrack_pixnum = ipix

      n_irradwvl = Irr_Data%nwaves(ipix)

      IF ( n_irradwvl <= 0 ) CYCLE

      saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

      ! Set up generic fitting arrays
      CALL omi_adjust_irradiance_data ( &
        Irr_Data, ipix, &
        omi_irradiance_ccdpix(1:n_irradwvl,ipix), &
        curr_sol_spec(wvl_idx:ccd_idx,1:n_irradwvl), &
        curr_sol_wav_avg, &
        yn_skip_pix, locerrstat )

      IF ( yn_skip_pix .OR. locerrstat >= pge_errstat_error ) THEN
        errstat = MAX ( errstat, locerrstat )
        omi_cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
          modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
          locerrstat )
        CYCLE
      END IF

      is_bad_pixel   = .FALSE.
      CALL solar_fit ( &   ! Solar wavelength calibration
        ctrl_n_fitres_loop(solcal_idx), ctrl_fitres_range(solcal_idx), n_irradwvl, &
        curr_sol_wav_avg, &
        curr_sol_spec(wvl_idx:ccd_idx,1:n_irradwvl), Slit_Half_Width_1e, &
        Slit_Asym_Factor, solcal_exval, solcal_itnum, chisquav, &
        is_bad_pixel, locerrstat )
      ! solar_fit modifies the following variables:
      !   curr_sol_spec, Slit_Half_Width_1e, Slit_Asym_Factor, solcal_exval,
      !   solcal_itnum, chisquav, is_bad_pixel, locerrstat

      IF ( is_bad_pixel .OR. locerrstat >= pge_errstat_error ) THEN
        errstat = MAX ( errstat, locerrstat )
        omi_cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
          modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
          locerrstat )
        CYCLE
      END IF

      ! -----------------------------------------------------------------------
      ! Save crucial variables for across-track reference in Earthshine fitting
      ! -----------------------------------------------------------------------
      if (solcal_exval == i4_missval) solcal_exval = i2_missval
      omi_solcal_chisq(ipix)                     = chisquav
      omi_solcal_pars (1:max_calfit_idx,ipix)    = fitvar_cal(1:max_calfit_idx)
      omi_solcal_xflag(ipix)                     = INT (solcal_exval, KIND=i2)
      omi_solcal_itnum(ipix)                     = INT (solcal_itnum, KIND=i2)

      ! ------------------------------------------------------------------------
      ! Save the processed solar spectrum in its original array. Note that the
      ! spectrum is now normalized, has bad pixels set to -1, and that the
      ! wavelength array is calibrated.
      ! ------------------------------------------------------------------------
      Irr_Data%wavelengths(1:n_irradwvl,ipix) = curr_sol_spec(wvl_idx,1:n_irradwvl)
      Irr_Data%spectrum(1:n_irradwvl, ipix) = curr_sol_spec(spc_idx,1:n_irradwvl)
      Irr_Data%avg_wavelengths(ipix) = curr_sol_wav_avg
      omi_irradiance_wght(1:n_irradwvl,ipix) = curr_sol_spec(sig_idx,1:n_irradwvl)

      addmsg = ''
      WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'SOLAR FIT          #', ipix, &
        ': hw 1/e = ', Slit_Half_Width_1e, '; e_asy = ', Slit_Asym_Factor, '; shift = ', &
        fitvar_cal(shi_idx), '; squeeze = ', fitvar_cal(squ_idx), '; exit val = ', &
        solcal_exval, '; iter num = ', solcal_itnum
      CALL error_check ( &
        0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), &
        vb_lev_omidebug, errstat )
      IF ( verb_thresh_lev >= vb_lev_screen  ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

    END DO XtrackSolCal
    errstat = MAX ( errstat, locerrstat )

    ! ----------------------------------------------------------------------------
    ! After the successful wavelength calibration of all cross-track solar spectra,
    ! we compute(d) the average spectrum over all fitted spectra. This is (was) an
    ! attempt to remove the stripes from the data, at least the part that is due
    ! to spikes in the solar irradiances.
    !
    ! HOWEVER, the performance severerly degraded fitting performance (10 times
    ! larger fitting residuals in the radiance fit), so it is not used any more.
    ! It is left here, commented out, for possible future improvement and use.
    ! ----------------------------------------------------------------------------
    ! !! CALL solar_xtrack_average ( &
    ! !!     first_pix, last_pix, MAXVAL(omi_nwav_irrad(first_pix:last_pix)),  errstat )

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_solar_calibration_loop

  SUBROUTINE solar_fit ( &
      n_fitres_loop, fitres_range, n_irradwvl, avg_sol_wav, &
      sol_spec, hw1e, e_asym, solcal_exval, solcal_itnum, chisquav, &
      is_bad_pixel, errstat )

    ! ***************************************************************
    !
    !   Perform solar wavelength calibration and slit width fitting
    !
    ! ***************************************************************

    USE OMSAO_precision_module, ONLY: i2, i4, r8
    USE OMSAO_parameters_module, ONLY: r8_missval, &
      i2_missval, i4_missval, downweight
    USE OMSAO_variables_module,  ONLY: yn_newshift, fitwavs, fitweights, &
      currspec, fitvar_cal, n_fitvar_cal, lobnd, upbnd, fitvar_cal_saved, &
      mask_fitvar_cal, fitvar_sol_init, sol_wav_avg, &
      max_itnum_sol, up_sunbnd, lo_sunbnd, tol, epsrel, epsabs, epsx
    USE OMSAO_indices_module, ONLY: wvl_idx, ccd_idx, asy_idx, hwe_idx, &
      shi_idx, sig_idx, squ_idx, spc_idx, max_calfit_idx
    USE OMSAO_errstat_module, ONLY: pge_errstat_ok
    use radiance_wavcal, only: solar_residuals
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: n_fitres_loop, n_irradwvl, fitres_range
    real (kind=r8), intent(in) :: avg_sol_wav

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), INTENT (OUT)   :: hw1e, e_asym, chisquav
    INTEGER (KIND=i4), INTENT (OUT)   :: solcal_exval
    INTEGER (KIND=i2), INTENT (OUT)   :: solcal_itnum

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,                                                   INTENT (OUT)   :: is_bad_pixel
    INTEGER (KIND=i4),                                         INTENT (INOUT) :: errstat
    REAL    (KIND=r8), DIMENSION(wvl_idx:ccd_idx,1:n_irradwvl), INTENT (INOUT) :: sol_spec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: locerrstat, i, j, locitnum, n_nozero_wgt
    REAL    (KIND=r8)  :: mean, sdev, loclim
    REAL    (KIND=r8), DIMENSION (n_irradwvl)         :: fitres
    REAL    (KIND=r8), DIMENSION (MAX_CALFIT_IDX)    :: fitvar

    type(optimizer_type) :: opt
    integer (kind=i4) :: return_status

    ! ----------------------------------------------------------------
    ! Initialize local error status variable; note that error handling
    ! is rudimentary in this subroutine - no error is reported.
    ! ----------------------------------------------------------------
    locerrstat = pge_errstat_ok

    solcal_exval = i4_missval
    solcal_itnum = i2_missval

    chisquav = r8_missval

    is_bad_pixel = .FALSE.

    ! --------------------------------------------------------------
    ! Calculate and iterate on the irradiance spectrum.
    ! --------------------------------------------------------------
    fitwavs   (1:n_irradwvl) = sol_spec(wvl_idx,1:n_irradwvl)
    fitweights(1:n_irradwvl) = sol_spec(sig_idx,1:n_irradwvl)
    currspec  (1:n_irradwvl) = sol_spec(spc_idx,1:n_irradwvl)

    ! -------------------------------------------------------------
    ! Initialize the fitting variables. FITVAR_CAL_SAVED has been
    ! set to the initial values in the calling routine. outside the
    ! pixel loop. Here we use FITVAR_CAL_SAVED, which will be
    ! updated with current values from the previous fit if that fit
    ! has gone well.
    ! -------------------------------------------------------------
    fitvar_cal(1:max_calfit_idx) = fitvar_cal_saved(1:max_calfit_idx)

    ! ---------------------------------------------------------
    ! Assign varied fitting variables to array passed to optimizer
    ! ---------------------------------------------------------
    fitvar = 0.0_r8 ; lobnd = 0.0_r8 ; upbnd = 0.0_r8
    n_fitvar_cal = 0
    DO i = 1, max_calfit_idx
      IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
        n_fitvar_cal  = n_fitvar_cal + 1
        mask_fitvar_cal(n_fitvar_cal) = i
        fitvar(n_fitvar_cal) = fitvar_cal(i)
        lobnd (n_fitvar_cal) = lo_sunbnd(i)
        upbnd (n_fitvar_cal) = up_sunbnd(i)
      END IF
    END DO

    ! --------------------------------------------------------------------
    ! Check whether we enough spectral points to carry out the fitting. If
    ! not, call it a bad pixel and return.
    ! --------------------------------------------------------------------
    IF ( n_fitvar_cal >= n_irradwvl ) THEN
      is_bad_pixel = .TRUE.  ;  RETURN
    END IF

    ! ---------------------------------------------------------------------
    ! Attempt to standardize the re-iteration with spectral points excluded
    ! that have fitting residuals larger than a pre-set window. Needs more
    ! thinking before it can replace a simple window determined empirically
    ! from fitting lots of spectra.
    ! ---------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! Refit if any part of the fitting residual computed above is larger than
    ! the pre-set window given by FITRES_RANGE. N_FITRES_LOOP must be set > 0
    ! since it determines the maximum number of re-iterations (we don't want
    ! to fit forever!).
    ! -----------------------------------------------------------------------

    loclim = 0.0_r8
    solcal_itnum = 0
    j = 0
    sol_wav_avg = avg_sol_wav

    call optimizer_open (opt, solar_residuals, n_fitvar_cal, return_status, &
                         mode=opt_bounded, tol=tol, epsabs=epsabs, epsrel=epsrel, epsx=epsx, &
                         param_min = lobnd(1:n_fitvar_cal), &
                         param_max = upbnd(1:n_fitvar_cal), &
                         param_mask = mask_fitvar_cal(1:n_fitvar_cal), &
                         max_num_iterations = max_itnum_sol)
    if (return_status < 0) then
      call err_message_error ("solar_fit: optimizer_open failed", errstat)
      return
    endif
    
    fit_loop: do
      call opt%optimize (opt, fitvar(1:n_fitvar_cal), n_fitvar_cal, &
                         fitres(1:n_irradwvl), n_irradwvl, return_status)
      locitnum = opt%num_iterations
      solcal_exval = return_status

      solcal_itnum = solcal_itnum + INT ( locitnum, KIND=i2 )
      j = j + 1

      n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_irradwvl)) ) )

      IF ( n_nozero_wgt > 0 ) THEN
        chisquav = SUM (fitres(1:n_irradwvl)**2)
      ELSE
        chisquav = r8_missval
      END IF

      if (1 < j .and. j <= n_fitres_loop) then
        if ( solcal_exval > 0 ) then
          fitvar_cal_saved(1:max_calfit_idx) = fitvar_cal(1:max_calfit_idx)
        else
          fitvar_cal_saved(1:max_calfit_idx) = fitvar_sol_init(1:max_calfit_idx)
        end if
        IF ( MAXVAL(ABS(fitres(1:n_irradwvl))) <= loclim ) exit fit_loop
      else if (j == 1) then
        mean = SUM  ( fitres(1:n_irradwvl) )                 / REAL(n_nozero_wgt,   KIND=r8)
        sdev = SQRT ( SUM ( (fitres(1:n_irradwvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
        loclim = REAL (fitres_range, KIND=r8)*sdev
        if (.not.((n_fitres_loop > 0) &
                  .and.(loclim > 0.0_r8) &
                  .and.(MAXVAL(ABS(fitres(1:n_irradwvl))) >= loclim) &
                  .and.(n_nozero_wgt > n_fitvar_cal))) exit fit_loop
      else
        exit fit_loop  ! (j > n_fitres_loop)
      endif

      WHERE ( ABS(fitres(1:n_irradwvl)) > loclim )
        fitweights(1:n_irradwvl) = downweight
      END WHERE

    enddo fit_loop

    call optimizer_close (opt, return_status)
    if (return_status < 0) then
      call err_message_error ("solar_fit: optimizer_close failed", errstat)
      return
    endif

    ! ---------------------------------------------------------------
    ! The following assignment makes sense only because FITVAR_CAL is
    ! updated with FITVAR (using the proper mask) in SPECTRUM_SOLAR.
    ! ---------------------------------------------------------------
    IF ( solcal_exval == opt_convergence_good ) THEN
      fitvar_cal_saved(1:max_calfit_idx) = fitvar_cal(1:max_calfit_idx)
    ELSE
      fitvar_cal_saved(1:max_calfit_idx) = fitvar_sol_init(1:max_calfit_idx)
    END IF

    ! ---------------------------------------------------------------
    ! Save shifted&squeezed wavelength array, and the fitting weights
    ! ---------------------------------------------------------------
    IF (yn_newshift .EQV. .true.) THEN !gga
      sol_spec(wvl_idx,1:n_irradwvl) = &
        (fitwavs (1:n_irradwvl) - fitvar_cal_saved(shi_idx) + &
        sol_wav_avg * fitvar_cal_saved(squ_idx)) /          &
        (1.0_r8 + fitvar_cal_saved(squ_idx))
    ELSE !gga
      sol_spec(wvl_idx,1:n_irradwvl) = fitwavs (1:n_irradwvl)
    END IF
    sol_spec(sig_idx,1:n_irradwvl) = fitweights (1:n_irradwvl)

    ! ------------------------------------------------
    !  Save the slit function parameters for later use
    ! in the undersampling correction.
    ! ------------------------------------------------
    hw1e   = fitvar_cal(hwe_idx)  ;  e_asym = fitvar_cal(asy_idx)

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE solar_fit

END MODULE OMSAO_solar_wavcal_module
