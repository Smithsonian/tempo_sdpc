MODULE OMSAO_solar_wavcal_module

  use tell_module
  implicit none

  private
  public xtrack_solar_calibration_loop

CONTAINS

  SUBROUTINE adjust_irradiance_data (irr, xtpix, irrad_ccd, &
                                     adj_wvl, adj_spec, adj_wgts, &
                                     avg_sol_wav, &
                                     do_skip_pix, errstat )

    !USE sao_pge_utils, ONLY: print_array
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: downweight, normweight, r4_missval
    USE OMSAO_indices_module,         ONLY: &
      qual_flag_mis, qual_flag_bad, qual_flag_err, qual_flag_sat
    use ctrlvars, only: yn_spectrum_norm
    !USE OMSAO_errstat_module
    USE ezspline_interpolation, ONLY: ezspline_1d_interpolation
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
    LOGICAL,                                               INTENT (OUT) :: do_skip_pix
    INTEGER (KIND=i4), DIMENSION (irr%nwaves(xtpix)),         INTENT (OUT) :: irrad_ccd
    real (kind=r8), dimension (irr%nwaves(xtpix)), intent(out) :: adj_wvl, adj_spec, adj_wgts
    real (kind=r8), intent(out) :: avg_sol_wav

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                       :: &
      i, j, imin1, imax1, imin2, imax2, j1, j2 ! locerrstat,
    LOGICAL                                                 :: have_good_window
    REAL    (KIND=r8), DIMENSION (irr%nwaves(xtpix))           :: weightsum
    REAL    (KIND=r8)                                       :: sol_spec_avg, asum, ssum
    INTEGER (KIND=i4) :: num_irr_wvl
    integer (kind=i4) :: bad_qflg_mask

    ! ----------------------------------------------
    ! Variables for separating the good from the bad
    ! ----------------------------------------------
    INTEGER (KIND=i4) :: ngood, nbad
    INTEGER (KIND=i4), DIMENSION (irr%nwaves(xtpix)) :: bad_idx
    REAL    (KIND=r8), DIMENSION (irr%nwaves(xtpix)) :: wvl_good, wvl_bad, spc_good, spc_bad

    if (errstat /= 0) return

    !locerrstat  = pge_errstat_ok
    do_skip_pix = .FALSE.

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
    adj_wvl(1:num_irr_wvl) = irr%wavelengths(1:num_irr_wvl, xtpix)
    adj_spec(1:num_irr_wvl) = irr%spectrum(1:num_irr_wvl, xtpix)
    irrad_ccd(        1:num_irr_wvl) = (/ (i, i = imin1, imax1) /)

    ! ---------------------------------------------------------
    ! Compute the weights. This is a bit tedious, as we have to
    ! check for a number of things that can go wrong. We start
    ! out assuming "all is well" and exclude/modify only those
    ! entries that are expected to give us trouble.
    ! ---------------------------------------------------------

    adj_wgts(1:num_irr_wvl) = normweight

    ! -----------------------------------
    ! Make sure wavelengths are ascending
    ! -----------------------------------
    DO i = 2, num_irr_wvl
      IF ( adj_wvl(i) <= adj_wvl(i-1) ) THEN
        adj_wvl(i) = adj_wvl(i-1) + 0.001_r8
        adj_wgts(i) = downweight
      END IF
    END DO

    ! -----------------------------
    ! No missing values in spectrum
    ! -----------------------------
    !call print_array (sol_spec (spc_idx, 1:num_irr_wvl), num_irr_wvl)

    WHERE ( adj_spec(1:num_irr_wvl) <= REAL( r4_missval, KIND=r8 ) )
      adj_wgts(1:num_irr_wvl) = downweight
      adj_spec(1:num_irr_wvl) = 0.0_r8
    END WHERE

    ! ----------------------------------------------------------------------
    ! Find the pixel quality flags (not assigned correctly in the L1 product
    ! as of 13 September 2004 and thus not used yet; tpk note to himself)
    ! Choice of flags is based on the recommendations of the L1b README.
    ! --------------------------------------------------------------------
    bad_qflg_mask = 0
    bad_qflg_mask = ior(bad_qflg_mask, qual_flag_mis)
    bad_qflg_mask = ior(bad_qflg_mask, qual_flag_bad)
    bad_qflg_mask = ior(bad_qflg_mask, qual_flag_err)
    bad_qflg_mask = ior(bad_qflg_mask, qual_flag_sat)

    where (iand (irr%qflags(1:num_irr_wvl, xtpix), bad_qflg_mask) /= 0)
      adj_wgts(1:num_irr_wvl) = downweight
      adj_spec(1:num_irr_wvl) = 0.0_r8
    end where

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
    WHERE ( adj_wgts(1:num_irr_wvl) /= downweight )
      weightsum = 1.0_r8
    END WHERE
    sol_spec_avg = SUM ( adj_spec(1:num_irr_wvl)*weightsum(1:num_irr_wvl) ) / &
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
           ABS(adj_spec(1:num_irr_wvl)) >= 100.0_r8 * sol_spec_avg )
      weightsum(1:num_irr_wvl) = 0.0_r8
      adj_wgts(1:num_irr_wvl) = downweight
      adj_spec(1:num_irr_wvl) = 0.0_r8
    ENDWHERE

    ! ------------------------------------------------------------------
    ! Recompute the solar spectrum average, because it may have changed.
    ! ------------------------------------------------------------------
    sol_spec_avg = SUM ( adj_spec(1:num_irr_wvl)*weightsum(1:num_irr_wvl) ) / &
      MAX(1.0_r8, SUM(weightsum(1:num_irr_wvl)))
    IF ( sol_spec_avg <= 0.0_r8 ) THEN
      do_skip_pix = .TRUE.
      sol_spec_avg = 1.0_r8
    END IF

    ! --------------------------------------
    ! Finally, normalize the solar spectrum.
    ! --------------------------------------
    IF ( yn_spectrum_norm ) &
      adj_spec(1:num_irr_wvl) = adj_spec(1:num_irr_wvl) / sol_spec_avg

    ! ---------------------------------------------
    ! Calculate SOL_WAV_AVG of measured solar spectra here,
    ! for use in calculated spectra.
    ! ---------------------------------------------
    asum = SUM ( adj_wvl(1:num_irr_wvl) * &
                ( adj_wgts(1:num_irr_wvl)*adj_wgts(1:num_irr_wvl) ) )
    ssum = SUM ( 1.0_r8 * &
                ( adj_wgts(1:num_irr_wvl)*adj_wgts(1:num_irr_wvl) ) )
    avg_sol_wav = asum / ssum

    ! ------------------------------------------------------------
    ! Count the good and the bad, and interpolate one to the other
    ! ------------------------------------------------------------
    ngood = 0 ; nbad = 0 ; bad_idx = 0
    wvl_good = 0.0_r8 ; wvl_bad = 0.0_r8 ; spc_good = 0.0_r8 ; spc_bad = 0.0_r8
    DO i = 1, num_irr_wvl
      IF ( adj_spec(i) <= 0.0_r8 .OR. adj_wgts(i) == downweight ) THEN
        nbad          = nbad + 1
        bad_idx(nbad) = i
        wvl_bad(nbad) = adj_wvl(i)
      ELSE
        ngood           = ngood + 1
        spc_good(ngood) = adj_spec(i)
        wvl_good(ngood) = adj_wvl(i)
      END IF
    END DO
    if (ngood == 0) THEN
      do_skip_pix = .TRUE.
    else  IF ( nbad > 0 ) THEN
      CALL ezspline_1d_interpolation (                      &
        ngood, wvl_good(1:ngood), spc_good(1:ngood),     &
        nbad, wvl_bad(1:nbad), spc_bad(1:nbad), errstat) !locerrstat )
      if (errstat /= 0) return
      DO i = 1, nbad
        j = bad_idx(i)
        adj_spec(j) = spc_bad(i)
        adj_wgts(j) = downweight
      END DO
    END IF

    ! ------------------------------------------------------------
    ! Anything outside the fitting window will receive Zero weight
    ! ------------------------------------------------------------
    ! (the CCD indices are absolute positions, i.e., unlikely to be "1:num_irr_wvl")
    ! ----------------------------------------------------------------------------
    IF ( imin2 > imin1 ) adj_wgts(1:imin2-imin1+1)         = downweight
    IF ( imax2 < imax1 ) adj_wgts(imax2-imin1+1:num_irr_wvl) = downweight

    ! ------------------------------------------------------------------------
    ! Also any window excluded by the user (specified in fitting control file)
    ! ------------------------------------------------------------------------
    j1 = irr%ccdpix_exclusion(1, xtpix) ; j2 = irr%ccdpix_exclusion(2, xtpix)
    j1 = j1 - imin1 + 1    ; j2 = j2 - imin1 + 1
    IF ( j1 >= 1 .AND. j2 <= num_irr_wvl ) adj_wgts(j1:j2) = downweight

    !IF ( locerrstat /= pge_errstat_ok ) errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE adjust_irradiance_data

  SUBROUTINE solar_fit ( &
      n_fitres_loop, fitres_range, n_irradwvl, avg_sol_wav, &
      sol_wvl, sol_spec, sol_wgts, sol_resid, &
      hw1e, e_asym, k, solcal_exval, solcal_itnum, chisquav, &
      is_bad_pixel, errstat )

    ! ***************************************************************
    !
    !   Perform solar wavelength calibration and slit width fitting
    !
    ! ***************************************************************

    USE OMSAO_precision_module, ONLY: i2, i4, r8
    USE OMSAO_parameters_module, ONLY: r8_missval, &
      i2_missval, i4_missval
    USE OMSAO_variables_module,  ONLY: fitvar_cal, fitvar_cal_saved, &
      fitvar_sol_init, & !sol_wav_avg,
      max_itnum_sol, up_sunbnd, lo_sunbnd
    use ctrlvars, only: yn_newshift
    USE OMSAO_indices_module, ONLY: asy_idx, hwe_idx, sgk_idx, &
      shi_idx, squ_idx, max_calfit_idx
    !USE OMSAO_errstat_module, ONLY: pge_errstat_ok
    use optimizer_interface_module, only: opt_convergence_good
    use wavecal

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: n_fitres_loop, n_irradwvl, fitres_range
    real (kind=r8), intent(in) :: avg_sol_wav

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), INTENT (OUT)   :: hw1e, e_asym, k, chisquav
    INTEGER (KIND=i4), INTENT (OUT)   :: solcal_exval
    INTEGER (KIND=i2), INTENT (OUT)   :: solcal_itnum

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,                                                   INTENT (OUT)   :: is_bad_pixel
    INTEGER (KIND=i4),                                         INTENT (INOUT) :: errstat
    real (kind=r8), dimension(:), intent(inout) :: sol_wvl, sol_spec, sol_wgts, sol_resid

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: locitnum !locerrstat,

    ! ----------------------------------------------------------------
    ! Initialize local error status variable; note that error handling
    ! is rudimentary in this subroutine - no error is reported.
    ! ----------------------------------------------------------------
    !locerrstat = pge_errstat_ok

    solcal_exval = i4_missval
    solcal_itnum = i2_missval

    chisquav = r8_missval

    is_bad_pixel = .FALSE.

    ! --------------------------------------------------------------
    ! Calculate and iterate on the irradiance spectrum.
    ! --------------------------------------------------------------

    ! -------------------------------------------------------------
    ! Initialize the fitting variables. FITVAR_CAL_SAVED has been
    ! set to the initial values in the calling routine. outside the
    ! pixel loop. Here we use FITVAR_CAL_SAVED, which will be
    ! updated with current values from the previous fit if that fit
    ! has gone well.
    ! -------------------------------------------------------------
    fitvar_cal(1:max_calfit_idx) = fitvar_cal_saved(1:max_calfit_idx)

    ! upon input loclim is the max per fit, on output it will be the total number
    locitnum = max_itnum_sol
    call wavecal_fit (sol_wvl, sol_spec, sol_wgts, sol_resid, n_irradwvl, avg_sol_wav, &
                      fitvar_cal, lo_sunbnd, up_sunbnd, max_calfit_idx, &
                      n_fitres_loop, real(fitres_range, kind=r8), &
                      is_bad_pixel, locitnum, chisquav, solcal_exval, errstat)
    if (errstat /= 0) return
    solcal_itnum = int (locitnum, kind=i2)

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
    if (yn_newshift) then !gga
      ! JCH: <comment-start>
      !      This is the assignment from the old code.
      !      I think 'sol_wav_avg' was supposed to be 'avg_sol_wav'.  When
      !      I made that replacement, the computed numbers didn't change, so I
      !      got rid of the reference to 'sol_wav_avg' from OMSAO_variables_module.
      !sol_wvl(1:n_irradwvl) = &
      !  (sol_wvl(1:n_irradwvl) - fitvar_cal_saved(shi_idx) &
      !   + sol_wav_avg * fitvar_cal_saved(squ_idx)) &
      !  / (1.0_r8 + fitvar_cal_saved(squ_idx))
      ! JCH: <comment-end>
      sol_wvl(1:n_irradwvl) = &
        (sol_wvl(1:n_irradwvl) - fitvar_cal_saved(shi_idx) &
         + avg_sol_wav * fitvar_cal_saved(squ_idx)) &
        / (1.0_r8 + fitvar_cal_saved(squ_idx))
    endif
    ! ------------------------------------------------
    !  Save the slit function parameters for later use
    ! in the undersampling correction.
    ! ------------------------------------------------
    hw1e   = fitvar_cal(hwe_idx)
    e_asym = fitvar_cal(asy_idx)
    k      = fitvar_cal(sgk_idx)

    RETURN

  END SUBROUTINE solar_fit

  SUBROUTINE xtrack_solar_calibration_loop (first_pix, last_pix, &
                                            save_wvl, save_resid, errstat)
    USE OMSAO_precision_module
    use ctrlvars, only : yn_diagnostic_run
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_omidata_module, ONLY: &
      omi_cross_track_skippix, &
      omi_solcal_chisq, omi_solcal_pars, omi_solcal_xflag, &
      omi_irradiance_wght, omi_irradiance_ccdpix
    USE OMSAO_indices_module, ONLY: &
      max_calfit_idx, shi_idx, squ_idx, solcal_idx
    USE OMSAO_parameters_module, ONLY: r8_missval, i2_missval, i4_missval, MAX_STR_LEN, &
      nwavel_max, nxtrack_max
    USE OMSAO_variables_module,  ONLY: Slit_Half_Width_1e, & ! verb_thresh_lev,
      Slit_Asym_Factor, Slit_Shape_Factor, fitvar_cal, fitvar_cal_saved,  &
      fitvar_sol_init, ctrl_n_fitres_loop, ctrl_fitres_range, &
      curr_xtrack_pixnum
    !USE OMSAO_errstat_module
    use irradiance_data, only : Irr_Data

    IMPLICIT NONE
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix

    ! -----------------
    ! Modified variable
    ! -----------------
    real (kind=r8), dimension(:,:), allocatable, intent(inout) :: save_wvl, save_resid
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i2)              :: solcal_itnum
    INTEGER   (KIND=i4)              :: locerrstat, ipix, solcal_exval, n_irradwvl
    CHARACTER (LEN=MAX_STR_LEN)         :: addmsg
    REAL      (KIND=r8)              :: chisquav, curr_sol_wav_avg
    LOGICAL                          :: do_skip_pix, is_bad_pixel
    real (kind=r8), dimension(:), allocatable :: adj_wvl, adj_spec, adj_wgts, adj_resid
    integer (kind=i4) :: adj_len
    integer locerr

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=29), PARAMETER :: modulename = 'xtrack_solar_calibration_loop'

    if (errstat /= 0) return

    omi_solcal_chisq = r8_missval
    omi_solcal_xflag = i2_missval

    fitvar_cal_saved(1:max_calfit_idx) = fitvar_sol_init(1:max_calfit_idx)

    ! ---------------------------------------------------------------
    ! Loop for solar wavelength calibration and slit function fitting
    ! ---------------------------------------------------------------

    if (yn_diagnostic_run) then
      allocate (save_wvl(nwavel_max, nxtrack_max), &
                save_resid(nwavel_max, nxtrack_max), stat=locerrstat)
      if (locerrstat /= 0) then
        call tell_error (tell_malloc_error, "xtrack_solar_calibration_loop: allocate failed", &
                         errstat)
        return
      endif
      save_wvl(:,:) = r8_missval
      save_resid(:,:) = r8_missval
    endif

    adj_len = 0
    XtrackSolCal: DO ipix = first_pix, last_pix

      !locerrstat = pge_errstat_ok

      curr_xtrack_pixnum = ipix

      n_irradwvl = Irr_Data%nwaves(ipix)

      IF ( n_irradwvl <= 0 ) CYCLE
      if (n_irradwvl > adj_len) then
        if (adj_len > 0) then
          deallocate (adj_wvl, adj_spec, adj_wgts, adj_resid, stat=locerr)
          if (locerr /= 0) then
            call tell_error (tell_malloc_error, &
                 "xtrack_solar_calibration_loop: deallocate failed", errstat)
            return
          endif
        endif
        allocate (adj_wvl(n_irradwvl), adj_spec(n_irradwvl), &
             adj_wgts(n_irradwvl), adj_resid(n_irradwvl), stat=locerr)
        if (locerr /= 0) then
          call tell_error (tell_malloc_error, &
               "xtrack_solar_calibration_loop: allocate failed", errstat)
          return
        endif
        adj_len = n_irradwvl
      endif

      saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

      ! Set up generic fitting arrays
      CALL adjust_irradiance_data ( &
        Irr_Data, ipix, &
        omi_irradiance_ccdpix(1:n_irradwvl,ipix), &
        adj_wvl, adj_spec, adj_wgts, &
        curr_sol_wav_avg, &
        do_skip_pix, errstat) !locerrstat )

      IF ( do_skip_pix .OR. errstat /= 0) then !locerrstat >= pge_errstat_error ) THEN
        !errstat = MAX ( errstat, locerrstat )
        omi_cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I5)') 'xtrack_solar_calibration_loop: SKIPPING cross track pixel #', ipix
        call tell_log (0, addmsg)
        !CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
        !                  modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
        !                  locerrstat )
        CYCLE
      END IF

      is_bad_pixel   = .FALSE.
      CALL solar_fit ( &   ! Solar wavelength calibration
        ctrl_n_fitres_loop(solcal_idx), ctrl_fitres_range(solcal_idx), n_irradwvl, &
        curr_sol_wav_avg, &
        adj_wvl, adj_spec, adj_wgts, adj_resid, &
        Slit_Half_Width_1e, Slit_Asym_Factor, Slit_Shape_Factor, &
        solcal_exval, solcal_itnum, chisquav, &
        is_bad_pixel, errstat) ! locerrstat )
      ! solar_fit modifies the following variables:
      !   adj_wvl, adj_spec, adj_wgts, Slit_Half_Width_1e, Slit_Asym_Factor, solcal_exval,
      !   solcal_itnum, chisquav, is_bad_pixel, locerrstat

      if (yn_diagnostic_run) then
        save_wvl(1:n_irradwvl,ipix) = adj_wvl(1:n_irradwvl)
        save_resid(1:n_irradwvl,ipix) = adj_resid(1:n_irradwvl)
      endif

      IF ( is_bad_pixel .OR. errstat /= 0) then ! locerrstat >= pge_errstat_error ) THEN
        !errstat = MAX ( errstat, locerrstat )
        omi_cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I5)') 'xtrack_solar_calibration_loop: SKIPPING cross track pixel #', ipix
        call tell_log (0, addmsg)
        !CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
        !  modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
        !  locerrstat )
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
      Irr_Data%wavelengths(1:n_irradwvl,ipix) = adj_wvl(1:n_irradwvl)
      Irr_Data%spectrum(1:n_irradwvl, ipix) = adj_spec(1:n_irradwvl)
      Irr_Data%avg_wavelengths(ipix) = curr_sol_wav_avg
      omi_irradiance_wght(1:n_irradwvl,ipix) = adj_wgts(1:n_irradwvl)

      addmsg = ''
      WRITE (addmsg, '(A,I4,6(A,1PE10.3),2(A,I9))') &
           'SOLAR FIT          #', ipix, &
           ': hw 1/e = ', Slit_Half_Width_1e, &
           '; e_asy = ', Slit_Asym_Factor, &
           '; k = ', Slit_Shape_Factor, &
           '; shift = ', fitvar_cal(shi_idx),&
           '; squeeze = ', fitvar_cal(squ_idx), &
           '; rms = ', sqrt(sum(adj_resid(1:n_irradwvl)**2)/real(n_irradwvl, kind=8)), &
           '; exit val = ', solcal_exval, '; iter num = ', solcal_itnum
      call tell_log (1, addmsg)
      !CALL error_check ( &
      !  0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), &
      !  vb_lev_omidebug, errstat )
      !IF ( verb_thresh_lev >= vb_lev_screen  ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

    END DO XtrackSolCal
    !errstat = MAX ( errstat, locerrstat )

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

    !errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_solar_calibration_loop

END MODULE OMSAO_solar_wavcal_module

