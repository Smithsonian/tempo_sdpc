MODULE OMSAO_solar_wavcal_module
  use optimizer_interface_module
  use errormodule

  IMPLICIT NONE

  private
  public xtrack_solar_calibration_loop

CONTAINS

  SUBROUTINE adjust_irradiance_data (irr, xtpix, irrad_ccd, &
                                     sol_spec, avg_sol_wav, &
                                     yn_skip_pix, errstat )

    !USE sao_pge_utils, ONLY: print_array
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: downweight, normweight, r4_missval
    USE OMSAO_indices_module,         ONLY: &
      wvl_idx, spc_idx, sig_idx, ccd_idx, &
      qflg_mis_idx, qflg_bad_idx, qflg_err_idx
    USE OMSAO_variables_module, ONLY: yn_spectrum_norm
    USE OMSAO_errstat_module
    USE ezspline_interpolation, ONLY: ezspline_1d_interpolation
    USE strutils, ONLY: convert_2bytes_to_16bits
    USE irradiance_data, only: Irradiance_Data_Type

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    type (Irradiance_Data_Type), intent(in) :: irr
    integer (kind=i4), intent(in) :: xtpix

    ! ----------------
    ! Output variables
    ! ----------------
    LOGICAL,                                               INTENT (OUT) :: yn_skip_pix
    INTEGER (KIND=i4), DIMENSION (irr%nwaves(xtpix)),         INTENT (OUT) :: irrad_ccd
    REAL    (KIND=r8), DIMENSION (ccd_idx,irr%nwaves(xtpix)), INTENT (OUT) :: sol_spec
    real (kind=r8), intent(out) :: avg_sol_wav

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i2), PARAMETER                            :: nbits = 16
    INTEGER (KIND=i4)                                       :: &
      i, j, locerrstat, imin1, imax1, imin2, imax2, j1, j2
    LOGICAL                                                 :: have_good_window
    INTEGER (KIND=i2), DIMENSION (irr%nwaves(xtpix),0:nbits-1) :: irrad_qflg_bit
    INTEGER (KIND=i2), DIMENSION (irr%nwaves(xtpix))           :: irrad_qflg_mask
    REAL    (KIND=r8), DIMENSION (irr%nwaves(xtpix))           :: weightsum
    REAL    (KIND=r8)                                       :: sol_spec_avg, asum, ssum
    INTEGER (KIND=i4) :: num_irr_wvl

    ! ----------------------------------------------
    ! Variables for separating the good from the bad
    ! ----------------------------------------------
    INTEGER (KIND=i4) :: ngood, nbad
    INTEGER (KIND=i4), DIMENSION (irr%nwaves(xtpix)) :: bad_idx
    REAL    (KIND=r8), DIMENSION (irr%nwaves(xtpix)) :: wvl_good, wvl_bad, spc_good, spc_bad

    locerrstat  = pge_errstat_ok
    yn_skip_pix = .FALSE.

    ! The total window
    imin1 = irr%ccdpix_selection (1,xtpix)
    imax1 = irr%ccdpix_selection (4,xtpix)
    ! The fitting window
    imin2 = irr%ccdpix_selection (2,xtpix)
    imax2 = irr%ccdpix_selection (3,xtpix)

    ! ---------------------------------------------------------------
    ! Assign irradiance spectrum to generic variables that are passed
    ! through the fitting routines down to the spectrum function.
    ! ---------------------------------------------------------------
    num_irr_wvl                          = irr%nwaves(xtpix)
    sol_spec(wvl_idx,1:num_irr_wvl) = irr%wavelengths(1:num_irr_wvl, xtpix)
    sol_spec(spc_idx,1:num_irr_wvl) = irr%spectrum(1:num_irr_wvl, xtpix)
    irrad_ccd(        1:num_irr_wvl) = (/ (i, i = imin1, imax1) /)

    ! ---------------------------------------------------------
    ! Compute the weights. This is a bit tedious, as we have to
    ! check for a number of things that can go wrong. We start
    ! out assuming "all is well" and exclude/modify only those
    ! entries that are expected to give us trouble.
    ! ---------------------------------------------------------

    sol_spec(sig_idx,1:num_irr_wvl) = normweight

    ! -----------------------------------
    ! Make sure wavelengths are ascending
    ! -----------------------------------
    DO i = 2, num_irr_wvl
      IF ( sol_spec(wvl_idx,i) <= sol_spec(wvl_idx,i-1) ) THEN
        sol_spec(wvl_idx,i) = sol_spec(wvl_idx,i-1) + 0.001_r8
        sol_spec(sig_idx,i) = downweight
      END IF
    END DO

    ! -----------------------------
    ! No missing values in spectrum
    ! -----------------------------
    !call print_array (sol_spec (spc_idx, 1:num_irr_wvl), num_irr_wvl)

    WHERE ( sol_spec(spc_idx,1:num_irr_wvl) <= REAL( r4_missval, KIND=r8 ) )
      sol_spec(sig_idx,1:num_irr_wvl) = downweight
      sol_spec(spc_idx,1:num_irr_wvl) = 0.0_r8
    END WHERE

    ! ----------------------------------------------------------------------
    ! Find the pixel quality flags (not assigned correctly in the L1 product
    ! as of 13 September 2004 and thus not used yet; tpk note to himself)
    ! ----------------------------------------------------------------------
    ! -------------------------------------------------------------------
    ! CAREFUL: Only 15 flags/positions (0:14) can be returned or else the
    !          conversion will result in a numeric overflow.
    ! -------------------------------------------------------------------
    CALL convert_2bytes_to_16bits (nbits-1, num_irr_wvl, &
                                   irr%qflags(1:num_irr_wvl, xtpix), &
                                   irrad_qflg_bit(1:num_irr_wvl,0:nbits-2) )

    ! --------------------------------------------------------------------
    ! Add contributions from various quality flags. Any CCD pixel that has
    ! a cumulative quality flag > 0 will be excluded form the fitting.
    !
    ! Choice of flags is based on the recommendations of the L1b README.
    ! --------------------------------------------------------------------
    irrad_qflg_mask(1:num_irr_wvl) = 0_i2
    irrad_qflg_mask(1:num_irr_wvl) =                  &
      irrad_qflg_bit(1:num_irr_wvl,qflg_mis_idx) + &   ! Missing pixel
      irrad_qflg_bit(1:num_irr_wvl,qflg_bad_idx) + &   ! Bad pixel
      irrad_qflg_bit(1:num_irr_wvl,qflg_err_idx) !+ &   ! Processing error
    !irrad_qflg_bit(1:num_irr_wvl,qflg_rts_idx)       ! RTS

    WHERE ( irrad_qflg_mask(1:num_irr_wvl) > 0_i2 )
      sol_spec(sig_idx,1:num_irr_wvl) = downweight
      sol_spec(spc_idx,1:num_irr_wvl) = 0.0_r8
    END WHERE

    ! -------------------------------------------------------------------------------
    ! Translate window limit wavelenghts into indices; making sure that Shift&Squeeze
    ! doesn't shift the wavelength array off the limits we set at the beginning. Do
    ! this iteratively until we have found good window limits. In the best case, all
    ! is well the first time around, but we might have to adjust the window margins
    ! if we fail to read all the data.
    ! -------------------------------------------------------------------------------
    have_good_window = .TRUE.

    ! -----------------------------------------------
    ! Compute normalization factor for solar spectrum
    ! -----------------------------------------------
    weightsum = 0.0_r8
    WHERE ( sol_spec(sig_idx,1:num_irr_wvl) /= downweight )
      weightsum = 1.0_r8
    END WHERE
    sol_spec_avg = SUM ( sol_spec(spc_idx,1:num_irr_wvl)*weightsum(1:num_irr_wvl) ) / &
      MAX(1.0_r8, SUM(weightsum(1:num_irr_wvl)))
    IF ( sol_spec_avg == 0.0_r8 ) sol_spec_avg = 1.0_r8

    ! -------------------------------------------------------------------------
    ! So far we have only taken care of/excluded any negative values in the
    ! spectrum, but there may abnormally high or low positive values also.
    ! Now we check for any values exceeding 100 times the average, which should
    ! be a large enough window to keep anything sensible and reject the real
    ! outliers.
    ! -------------------------------------------------------------------------
    WHERE ( weightsum(1:num_irr_wvl) /= 0.0_r8 .AND. &
           ABS(sol_spec(spc_idx,1:num_irr_wvl)) >= 100.0_r8 * sol_spec_avg )
      weightsum(1:num_irr_wvl) = 0.0_r8
      sol_spec(sig_idx,1:num_irr_wvl) = downweight
      sol_spec(spc_idx,1:num_irr_wvl) = 0.0_r8
    ENDWHERE

    ! ------------------------------------------------------------------
    ! Recompute the solar spectrum average, because it may have changed.
    ! ------------------------------------------------------------------
    sol_spec_avg = SUM ( sol_spec(spc_idx,1:num_irr_wvl)*weightsum(1:num_irr_wvl) ) / &
      MAX(1.0_r8, SUM(weightsum(1:num_irr_wvl)))
    IF ( sol_spec_avg <= 0.0_r8 ) THEN
      yn_skip_pix = .TRUE.
      sol_spec_avg = 1.0_r8
    END IF

    ! --------------------------------------
    ! Finally, normalize the solar spectrum.
    ! --------------------------------------
    IF ( yn_spectrum_norm ) &
      sol_spec(spc_idx,1:num_irr_wvl) = sol_spec(spc_idx,1:num_irr_wvl) / sol_spec_avg

    ! ---------------------------------------------
    ! Calculate SOL_WAV_AVG of measured solar spectra here,
    ! for use in calculated spectra.
    ! ---------------------------------------------
    asum = SUM ( sol_spec(wvl_idx,1:num_irr_wvl) * &
                ( sol_spec(sig_idx,1:num_irr_wvl)*sol_spec(sig_idx,1:num_irr_wvl) ) )
    ssum = SUM ( 1.0_r8 * &
                ( sol_spec(sig_idx,1:num_irr_wvl)*sol_spec(sig_idx,1:num_irr_wvl) ) )
    avg_sol_wav = asum / ssum

    ! ------------------------------------------------------------
    ! Count the good and the bad, and interpolate one to the other
    ! ------------------------------------------------------------
    ngood = 0 ; nbad = 0 ; bad_idx = 0
    wvl_good = 0.0_r8 ; wvl_bad = 0.0_r8 ; spc_good = 0.0_r8 ; spc_bad = 0.0_r8
    DO i = 1, num_irr_wvl
      IF ( sol_spec(spc_idx,i) <= 0.0_r8 .OR. sol_spec(spc_idx,i) == downweight ) THEN
        nbad          = nbad + 1
        bad_idx(nbad) = i
        wvl_bad(nbad) = sol_spec(wvl_idx,i)
      ELSE
        ngood           = ngood + 1
        spc_good(ngood) = sol_spec(spc_idx,i)
        wvl_good(ngood) = sol_spec(wvl_idx,i)
      END IF
    END DO
    IF ( nbad > 0 ) THEN
      CALL ezspline_1d_interpolation (                      &
        ngood, wvl_good(1:ngood), spc_good(1:ngood),     &
        nbad, wvl_bad(1:nbad), spc_bad(1:nbad), locerrstat )
      DO i = 1, nbad
        j = bad_idx(i)
        sol_spec(spc_idx,j) = spc_bad(i)
        sol_spec(sig_idx,j) = downweight
      END DO
    END IF

    ! ------------------------------------------------------------
    ! Anything outside the fitting window will receive Zero weight
    ! ------------------------------------------------------------
    ! (the CCD indices are absolute positions, i.e., unlikely to be "1:num_irr_wvl")
    ! ----------------------------------------------------------------------------
    IF ( imin2 > imin1 ) sol_spec(sig_idx,1:imin2-imin1+1)         = downweight
    IF ( imax2 < imax1 ) sol_spec(sig_idx,imax2-imin1+1:num_irr_wvl) = downweight

    ! ------------------------------------------------------------------------
    ! Also any window excluded by the user (specified in fitting control file)
    ! ------------------------------------------------------------------------
    j1 = irr%ccdpix_exclusion(1, xtpix) ; j2 = irr%ccdpix_exclusion(2, xtpix)
    j1 = j1 - imin1 + 1    ; j2 = j2 - imin1 + 1
    IF ( j1 >= 1 .AND. j2 <= num_irr_wvl ) sol_spec(sig_idx,j1:j2) = downweight

    IF ( locerrstat /= pge_errstat_ok ) errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE adjust_irradiance_data

  SUBROUTINE xtrack_solar_calibration_loop ( first_pix, last_pix, errstat )

    USE OMSAO_precision_module
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_omidata_module, ONLY: &
      omi_cross_track_skippix, &
      omi_solcal_chisq, omi_solcal_pars, omi_solcal_xflag, &
      omi_irradiance_wght, omi_irradiance_ccdpix
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
    omi_solcal_xflag = i2_missval

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
      CALL adjust_irradiance_data ( &
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
      currspec, fitvar_cal, n_fitvar_cal, fitvar_cal_saved, &
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
    REAL    (KIND=r8), DIMENSION (MAX_CALFIT_IDX)    :: fitvar, lobnd, upbnd

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
