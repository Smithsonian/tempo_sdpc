MODULE fitting_loops

  use tell_module

  private
  public xtrack_radiance_wvl_calibration, xtrack_radiance_fitting_loop

CONTAINS
  SUBROUTINE xtrack_radiance_wvl_calibration (first_pix, last_pix, &
                                              nxtrack, n_max_rspec, &
                                              save_wvl, save_resid, &
                                              errstat)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: &
      max_calfit_idx, max_rs_idx, hwe_idx, asy_idx, sgk_idx, &
      shi_idx, squ_idx, solar_idx,    &
      radcal_idx
    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, downweight, &
      NWAVEL_MAX, NXTRACK_MAX, r8_missval
    USE OMSAO_variables_module,  ONLY:  &
      Slit_Half_Width_1e, Slit_Asym_Factor, Slit_Shape_Factor, &
      database, fitvar_cal, fitvar_cal_saved, &  ! sol_wav_avg,
      fitvar_rad_init, ctrl_n_fitres_loop, ctrl_fitres_range, &
      curr_xtrack_pixnum, refspecs_original, radwavcal_freq
    use ctrlvars, only: yn_radiance_reference, yn_diagnostic_run, yn_do_he5_output
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_radiance_ref_module, ONLY: omi_adjust_radiance_data
    USE OMSAO_omidata_module, ONLY: omi_nwav_radref, omi_nwav_rad, &
      omi_irradiance_wght, omi_radiance_wavl, &
      omi_radiance_spec, omi_radiance_qflg, n_omi_radwvl, omi_radref_wavl, &
      omi_radref_spec, omi_radref_qflg, rad_ccdpix_selection, &
      rad_ccdpix_exclusion, omi_radiance_ccdpix, omi_cross_track_skippix, &
      omi_radcal_pars, omi_radcal_xflag, &
      omi_radcal_chisq, &
      omi_radref_wght, omi_database, n_omi_database_wvl, &
      omi_database_wvl, & ! omi_radref_wav_avg,
      omi_solcal_pars
    USE prepare_databases, ONLY: prep_databases
    !USE OMSAO_errstat_module
    USE he5_output_tools, ONLY: he5_write_omi_database
    use output_tools, only : write_refspec_database
    USE sao_pge_utils, ONLY: interpolation
    USE radiance_wavcal, ONLY: radiance_wavecal
    USE irradiance_data, ONLY: Irr_Data
    !USE EZspline_obj
    !USE EZspline

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix, nxtrack, n_max_rspec

    ! ---------------
    ! Output variable
    ! ---------------
    !INTEGER (KIND=i4), INTENT (OUT) :: n_comm_wvl_out

    ! -----------------
    ! Modified variable
    ! -----------------
    real (kind=r8), dimension(:,:), allocatable, intent(inout) :: save_wvl, save_resid
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i2)      :: radcal_itnum
    INTEGER   (KIND=i4)      :: locerrstat, ipix, radcal_exval, i, imax, n_ref_wvl !, nxtloc, xtr_add
    REAL      (KIND=r8)      :: chisquav, rad_spec_avg
    LOGICAL                  :: do_skip_pix, is_bad_pixel, did_full_range
    CHARACTER (LEN=MAX_STR_LEN) :: addmsg
    INTEGER (KIND=i4), DIMENSION (4)            :: select_idx
    INTEGER (KIND=i4), DIMENSION (2)            :: exclud_idx
    REAL    (KIND=r8), DIMENSION (n_max_rspec) :: ref_wvl, ref_spc, ref_wgt
    real (kind=r8), dimension(:), allocatable :: adj_wvls, adj_spec, adj_wgts, adj_resid
    integer (kind=i4) :: adj_num, adj_num_allocated, adj_num_max
    integer (kind=i4) :: n_irradwvl
    integer :: locerr

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=31), PARAMETER :: modulename = 'xtrack_radiance_wvl_calibration'

    if (errstat /= 0) return

    !locerrstat = 0 !pge_errstat_ok

    fitvar_cal_saved(1:max_calfit_idx) = fitvar_rad_init(1:max_calfit_idx)

    ! -------------------------------------------------
    ! Set the number of wavelengths for the common mode
    ! -------------------------------------------------
    !
    !JCH: The value of n_comm_wvl should have been correctly determined
    !     already, in a way that accounts for the fitting window size,
    !     so there's no need to re-define it in this way.
    !
    !n_comm_wvl_out = MAXVAL ( omi_nwav_radref(first_pix:last_pix) )
    !IF ( MAXVAL(omi_nwav_rad(first_pix:last_pix,0)) > n_comm_wvl_out ) &
    !  n_comm_wvl_out = MAXVAL(omi_nwav_rad(first_pix:last_pix,0))

    if (yn_diagnostic_run) then
      allocate (save_wvl(nwavel_max,nxtrack), &
                save_resid(nwavel_max,nxtrack), &
                stat=locerrstat)
      if (locerrstat /= 0) then
        call tell_error (tell_malloc_error, &
                         "xtrack_radiance_wvl_calibration: allocate failed", &
                         errstat)
        return
      endif
      save_wvl(:,:) = r8_missval
      save_resid(:,:) = r8_missval
    endif

    adj_num_allocated = 0
    adj_num_max = 0
    ! --------------------------------
    ! Loop over cross-track positions.
    ! --------------------------------
    XTrackWavCal: DO ipix = first_pix, last_pix

      !locerrstat = 0 ! pge_errstat_ok

      curr_xtrack_pixnum = ipix

      ! ---------------------------------------------------------------------
      ! If we already determined that this cross track pixel position carries
      ! an error, we don't even have to start processing.
      ! ---------------------------------------------------------------------
      IF ( omi_cross_track_skippix(ipix) ) CYCLE

      ! ---------------------------------------------------------------------------
      ! For each cross-track position we have to initialize the saved Shift&Squeeze
      ! ---------------------------------------------------------------------------
      saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

      ! ----------------------------------------------------
      ! Assign number of radiance and irradiance wavelengths
      ! ----------------------------------------------------
      n_irradwvl = Irr_Data % nwaves(ipix)
      n_omi_radwvl   = omi_nwav_rad  (ipix,0)
      adj_num = n_omi_radwvl
      if (adj_num > adj_num_max) adj_num_max = adj_num

      ! -----------------------------------------------------------------
      ! tpk: Should the following be "> n_fitvar_rad"??? No, because that
      !      value is set only inside OMI_ADJUST_RADIANCE_DATA!!!
      ! -----------------------------------------------------------------
      IF ( n_irradwvl <= 0 .OR. adj_num <= 0 ) CYCLE

      ! ---------------------------------------------------------------
      ! Restore solar fitting variables for across-track reference in
      ! Earthshine fitting. Use the Radiance References if appropriate.
      ! ---------------------------------------------------------------
      !sol_wav_avg = omi_radref_wav_avg(ipix)    ! JCH: no need to set this here
      Slit_Half_Width_1e = omi_solcal_pars(hwe_idx,ipix)
      Slit_Asym_Factor = omi_solcal_pars(asy_idx,ipix)
      Slit_Shape_Factor = omi_solcal_pars(sgk_idx,ipix)

      ! -----------------------------------------------------
      ! Assign (hopefully predetermined) "reference" weights.
      ! -----------------------------------------------------
      ref_wgt(1:n_irradwvl) = omi_irradiance_wght(1:n_irradwvl,ipix)

      ! -----------------------------------------------------
      ! Catch the possibility that adj_num > N_IRRADWVL
      ! -----------------------------------------------------
      IF ( adj_num > n_irradwvl ) THEN
        i = adj_num - n_irradwvl
        ref_wgt(n_irradwvl+1:adj_num) = downweight
        ! n_irradwvl = adj_num !! BAD BAD  --JED
      END IF

      ! ---------------------------------------------------------------
      ! If a Radiance Reference is being used, then it must be calibrated
      ! rather than the swath line that has been read.
      ! ---------------------------------------------------------------
      IF ( yn_radiance_reference ) THEN
        n_omi_radwvl = omi_nwav_radref(ipix)   ! JCH let's get the array length right
        adj_num = n_omi_radwvl
        if (adj_num > adj_num_max) adj_num_max = adj_num

        omi_radiance_wavl(1:adj_num,ipix,0) = omi_radref_wavl(1:adj_num,ipix)
        omi_radiance_spec(1:adj_num,ipix,0) = omi_radref_spec(1:adj_num,ipix)
        omi_radiance_qflg(1:adj_num,ipix,0) = omi_radref_qflg(1:adj_num,ipix)
      END IF

      ! ---------------------------------------------------------------------------
      ! Set up generic fitting arrays. Remember that OMI_RADIANCE_XXX arrays are
      ! 3-dim with the last dimension being the scan line numbers. For the radiance
      ! wavelength calibration we only have one scan line at index "0".
      ! ---------------------------------------------------------------------------

      ! reallocate buffers if needed
      if (adj_num > adj_num_allocated) then
        if (adj_num_allocated > 0) then
          deallocate (adj_wvls, adj_spec, adj_wgts, adj_resid, stat=locerr)
          if (locerr /= 0) then
            call tell_error (tell_malloc_error, &
                "xtrack_radiance_wvl_calibration: deallocate failed", errstat)
            return
          endif
        endif
        allocate (adj_wvls(adj_num), adj_spec(adj_num), adj_wgts(adj_num), &
                  adj_resid(adj_num), stat=locerr)
        if (locerr /= 0) then
          call tell_error (tell_malloc_error, &
              "xtrack_radiance_wvl_calibration: allocate failed", errstat)
          return
        endif
        adj_wvls(1:adj_num) = r8_missval
        adj_spec(1:adj_num) = r8_missval
        adj_wvls(1:adj_num) = r8_missval
        adj_num_allocated = adj_num
      endif
      adj_wvls(1:adj_num) = omi_radiance_wavl (1:adj_num, ipix, 0)
      adj_spec(1:adj_num) = omi_radiance_spec (1:adj_num, ipix, 0)
      adj_wgts(1:adj_num) = ref_wgt (1:adj_num)

      select_idx(1:4) = rad_ccdpix_selection(ipix,1:4)
      exclud_idx(1:2) = rad_ccdpix_exclusion(ipix,1:2)

      CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
        select_idx(1:4), exclud_idx(1:2),            &
        adj_num,                                &
        adj_wvls(1:adj_num), adj_spec(1:adj_num), adj_wgts(1:adj_num), &
        omi_radiance_qflg  (1:adj_num,ipix,0),  &
        omi_radiance_ccdpix(1:adj_num,ipix,0),  &
        rad_spec_avg, do_skip_pix )

      ! ------------------------------------------------------------------------------------
      IF (do_skip_pix) then
        omi_cross_track_skippix (ipix) = .FALSE.
        addmsg = ''
        WRITE (addmsg, '(A,I5)') 'xtrack_radiance_wvl_calibration: SKIPPING cross track pixel #', ipix
        call tell_log (0, addmsg)
      END IF

      ! ------------------------------------------------
      ! Logincal to skip radiance wavelength calibration
      ! ------------------------------------------------
      IF (radwavcal_freq > 0) THEN
        is_bad_pixel = .FALSE.

        CALL radiance_wavecal ( & ! Radiance wavelength calibration
          ipix, adj_num, &
          adj_wvls(1:adj_num), adj_spec(1:adj_num), &
          adj_wgts(1:adj_num), adj_resid(1:adj_num), &
          ctrl_n_fitres_loop(radcal_idx), ctrl_fitres_range(radcal_idx), &
          radcal_exval, radcal_itnum, chisquav, is_bad_pixel, errstat)

        if (yn_diagnostic_run) then
          save_wvl(1:adj_num,ipix) = adj_wvls(1:adj_num)
          save_resid(1:adj_num,ipix) = adj_resid(1:adj_num)
        endif

        IF ( is_bad_pixel .OR. errstat /= 0) then
          omi_cross_track_skippix (ipix) = .FALSE.
          addmsg = ''
          WRITE (addmsg, '(A,I5)') 'xtrack_radiance_wvl_calibration: SKIPPING cross track pixel #', ipix
          call tell_log (0, addmsg)
          call tell_set_error (0)
          errstat = 0
        END IF

        addmsg = ''
        WRITE (addmsg, '(A,I4,6(A,1PE10.3),2(A,I9))') 'RADIANCE Wavcal    #', ipix, &
          ': hw 1/e = ', Slit_Half_Width_1e, '; e_asy = ', Slit_Asym_Factor, &
          '; k = ', Slit_Shape_Factor, &
          '; shift = ', fitvar_cal(shi_idx), '; squeeze = ', fitvar_cal(squ_idx), &
          '; rms = ', sqrt(sum(adj_resid(1:adj_num)**2)/real(adj_num, kind=8)), &
          '; exit val = ', radcal_exval, '; iter num = ', radcal_itnum
        call tell_log (1, addmsg)

        ! ---------------------------------
        ! Save crucial variables for output
        ! ---------------------------------
        omi_radcal_pars (1:max_calfit_idx,ipix) = fitvar_cal(1:max_calfit_idx)
        omi_radcal_xflag(ipix)                  = INT (radcal_exval, KIND=i2)
        omi_radcal_chisq(ipix)                  = chisquav
        ! -----------------------------------------------------------------------
      ENDIF

      IF (yn_radiance_reference) THEN
        n_ref_wvl = adj_num
        ref_wvl(1:adj_num) = adj_wvls(1:adj_num)
        ref_spc(1:adj_num) = adj_spec(1:adj_num)
        ref_wgt(1:adj_num) = adj_wgts(1:adj_num)

        omi_nwav_radref(ipix)           = adj_num
        omi_radref_wavl(1:adj_num,ipix) = adj_wvls(1:adj_num)
        omi_radref_spec(1:adj_num,ipix) = adj_spec(1:adj_num)
        omi_radref_wght(1:adj_num,ipix) = adj_wgts(1:adj_num)
      ELSE
        n_ref_wvl = n_irradwvl
        ref_wvl(1:n_ref_wvl) = Irr_Data%wavelengths(1:n_ref_wvl,ipix)
        ref_spc(1:n_ref_wvl) = Irr_Data%spectrum(1:n_ref_wvl,ipix)
        ref_wgt(1:n_ref_wvl) = omi_irradiance_wght(1:n_ref_wvl,ipix)
      ENDIF

      ! ----------------------------------------------------
      ! Spline reference spectra to current wavelength grid.
      ! ----------------------------------------------------
      Call prep_databases ( &
        ipix, n_ref_wvl, ref_wvl(1:n_ref_wvl), ref_spc(1:n_ref_wvl), &
        adj_num, adj_wvls(1:adj_num), n_max_rspec, errstat) ! locerrstat )
      ! --------------------------------------------------------------------------------
      if (errstat /= 0) then !exit XTrackWavCal
        addmsg = ''
        write (addmsg, '(A,I5)') 'xtrack_radiance_wvl_calibration: prep_databases failed: SKIPPING cross track pixel #', ipix
        call tell_log (0, addmsg)
        call tell_set_error (0)
        errstat = 0
        cycle
      endif

      ! ---------------------------------------------------------
      ! Save DATABASE in OMI_DATABASE for radiance fitting loops.
      ! ---------------------------------------------------------
      omi_database (1:adj_num,ipix,1:max_rs_idx) = database (1:adj_num,1:max_rs_idx)
      n_omi_database_wvl(ipix) = adj_num
      omi_database_wvl(1:adj_num, ipix) = adj_wvls(1:adj_num)
      if (adj_num < adj_num_max) then
        omi_database (adj_num+1:adj_num_max,ipix,1:max_rs_idx) = r8_missval
        omi_database_wvl (adj_num+1:adj_num_max, ipix) = r8_missval
      endif

      ! ----------------------------------------------------------------------
      ! Update the radiance reference with the wavelength calibrated values.
      ! ----------------------------------------------------------------------
      IF ( yn_radiance_reference ) THEN
        !omi_radref_wavl(1:adj_num,ipix) = adj_wvls(1:adj_num)
        !omi_radref_spec(1:adj_num,ipix) = adj_spec(1:adj_num)
        omi_radref_wght(adj_num+1:nwavel_max,ipix) = downweight

        ! --------------------------------------------------------
        ! Update the solar spectrum entry in OMI_DATABASE. First
        ! re-sample the solar reference spectrum to the OMI grid
        ! then assign to data base.
        !
        ! We need to keep the irradiance spectrum because we still
        ! have to fit the radiance reference, and we can't really
        ! do that against itself. In a later module the irradiance
        ! is replaced by the radiance reference.
        ! --------------------------------------------------------

        ! ------------------------------------------------------------------
        ! Prevent failure of interpolation by finding the maximum wavelength
        ! of the irradiance wavelength array.
        ! ------------------------------------------------------------------
        imax = MAXVAL ( MAXLOC (Irr_Data%wavelengths(1:n_irradwvl,ipix) ) )
        !imin = MINVAL ( MINLOC ( omi_irradiance_wavl(1:imax,          ipix) ) )

        CALL interpolation ( &
          imax, Irr_Data%wavelengths(1:imax,ipix),                     &
          Irr_Data%spectrum(1:imax,ipix),                           &
          adj_num, omi_database_wvl(1:adj_num,ipix),              &
          omi_database(1:adj_num,ipix,solar_idx),                   &
          'endpoints', 0.0_r8, did_full_range, errstat) ! locerrstat )

        IF (errstat /= 0) then ! locerrstat >= pge_errstat_error ) THEN
          !errstat = MAX ( errstat, locerrstat )
          omi_cross_track_skippix (ipix) = .TRUE.
          addmsg = ''
          WRITE (addmsg, '(A,I5)') 'xtrack_radiance_wvl_calibration: SKIPPING cross track pixel #', ipix
          call tell_log (0, addmsg)
          call tell_set_error (0)
          errstat = 0
          !CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
          !                  modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
          !                  locerrstat )
          CYCLE
        END IF

      END IF

    END DO XTrackWavCal

    ! CCM Write splined/convolved databases if necessary
    IF(yn_diagnostic_run) THEN
      if (yn_do_he5_output) then
      ! omi_database maybe omi_database_wvl?
        CALL he5_write_omi_database(omi_database(1:adj_num,1:nxtrack_max,1:max_rs_idx), &
                                  omi_database_wvl(1:adj_num, 1:nxtrack_max), &
                                  max_rs_idx, adj_num, nxtrack_max, errstat)
      endif
      ! JCH:  I don't think adj_num should be used to define the subarrays that
      ! get written out because it's value might have changed with each
      ! pass through the above loop.  I added adj_num_max and will use that.
      ! (Note that this makes the he5 and netcdf output files different.)
      call write_refspec_database (omi_database(1:adj_num_max,1:nxtrack_max,1:max_rs_idx), &
                                   omi_database_wvl(1:adj_num_max,1:nxtrack_max), &
                                   refspecs_original(1:max_rs_idx), &
                                   max_rs_idx, adj_num_max, nxtrack, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_write_error, &
                         'xtrack_radiance_wvl_calibration: error writing refspec database', &
                         errstat)
        return
      endif
    ENDIF

    !errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_radiance_wvl_calibration

  SUBROUTINE check_wavelength_overlap ( &
      n_fitvar_rad, n_sol_wvl, irradiance_wvl, n_rad_wvl, radiance_wvl, &
      do_cycle_this_pix )

    USE OMSAO_precision_module,  ONLY: i4, r8

    IMPLICIT NONE

    ! Input variables
    INTEGER (KIND=i4),                        INTENT (IN) :: n_sol_wvl, n_rad_wvl, n_fitvar_rad
    REAL    (KIND=r8), DIMENSION (n_sol_wvl), INTENT (IN) :: irradiance_wvl
    REAL    (KIND=r8), DIMENSION (n_rad_wvl), INTENT (IN) :: radiance_wvl

    ! Output variable
    LOGICAL, INTENT (OUT) :: do_cycle_this_pix

    ! Local variables
    INTEGER (KIND=i4) :: j, n_overlap1, n_overlap2

    do_cycle_this_pix = .FALSE.

    IF ( radiance_wvl(1)         >= irradiance_wvl(n_sol_wvl) .OR. &
        radiance_wvl(n_rad_wvl) <= irradiance_wvl(1)                 ) THEN
      do_cycle_this_pix = .TRUE.
      RETURN
    END IF

    n_overlap1 = 0
    DO j = 1, n_rad_wvl
      IF ( radiance_wvl(j) >= irradiance_wvl(1)         .AND. &
          radiance_wvl(j) <= irradiance_wvl(n_sol_wvl)         ) n_overlap1 = n_overlap1 + 1
    END DO
    IF ( n_overlap1 < n_fitvar_rad ) THEN
      do_cycle_this_pix = .TRUE.
      RETURN
    END IF

    n_overlap2 = 0
    DO j = 1, n_sol_wvl
      IF ( irradiance_wvl(j) >= radiance_wvl(1)         .AND. &
          irradiance_wvl(j) <= radiance_wvl(n_rad_wvl)         ) n_overlap2 = n_overlap2 + 1
    END DO
    IF ( n_overlap2 < n_fitvar_rad ) THEN
      do_cycle_this_pix = .TRUE.
      RETURN
    END IF

    RETURN
  END SUBROUTINE check_wavelength_overlap

  SUBROUTINE xtrack_radiance_fitting_loop (pge_idx, &
      n_max_rspec, first_pix, last_pix, iloop,                &
      n_fitvar_rad, allfit_cols, allfit_errs, corr_matrix, &
      target_var, errstat, fitspc_out, fitspc_out_dim0                 )

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: &
      wvl_idx, spc_idx, & ! shi_idx,
      o3_t1_idx, o3_t3_idx, hwe_idx, asy_idx, sgk_idx, &
      pge_o3_idx, & !pge_hcho_idx, &
      solar_idx, radfit_idx, & !pge_gly_idx, &
      max_rs_idx! , mxs_idx ,max_calfit_idx, lbe_idx, hcho_idx,
    USE OMSAO_parameters_module, ONLY: &
      i2_missval, r8_missval, nxtrack_max
    USE OMSAO_variables_module,  ONLY:  &
      database, curr_sol_spec, n_rad_wvl, & ! sol_wav_avg,
      Slit_Half_Width_1e, Slit_Asym_Factor, Slit_Shape_Factor,    &
      n_database_wvl, ctrl_n_fitres_loop, ctrl_fitres_range,     &
      szamax, n_fincol_idx, curr_xtrack_pixnum!, fitvar_rad
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_prefitcol_module, ONLY: prefit_type, copy_prefit_values
    USE OMSAO_omidata_module, ONLY: omi_database_wvl, omi_radiance_wavl, &
      omi_database, rad_ccdpix_selection, rad_ccdpix_exclusion, &
      omi_fitconv_flag, omi_itnum_flag, omi_radfit_chisq, &
      omi_fit_rms, omi_radiance_spec, omi_column_amount, omi_column_uncert, &
      omi_o3_amount, omi_o3_uncert, n_omi_radwvl, &
      omi_szenith, n_omi_database_wvl, omi_nwav_rad, &
      omi_radiance_qflg, omi_cross_track_skippix, & ! omi_radref_wav_avg,
      omi_solcal_pars, omi_radiance_ccdpix, omi_radref_wght
    USE OMSAO_radiance_ref_module, ONLY: omi_adjust_radiance_data
    !USE OMSAO_errstat_module
    USE radiance_fit, ONLY: fit_radiance

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: &
      pge_idx, iloop, first_pix, last_pix, n_max_rspec, n_fitvar_rad, &
      fitspc_out_dim0

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    REAL    (KIND=r8), INTENT (OUT  ), DIMENSION (n_fitvar_rad,first_pix:last_pix) :: &
      allfit_cols, allfit_errs, corr_matrix

    ! ---------------------------------------------------------
    ! Optional output variable (fitted variable for target gas)
    ! ---------------------------------------------------------
    REAL (KIND=r8), DIMENSION(n_fincol_idx,first_pix:last_pix), INTENT (OUT) :: target_var

    ! CCM Output fit spectra
    !REAL (KIND=r8), DIMENSION(n_comm_wvl,nxtrack_max,4), INTENT (OUT) :: fitspc_out
    REAL (KIND=r8), DIMENSION(fitspc_out_dim0,nxtrack_max,4), INTENT (OUT) :: fitspc_out

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ipix, radfit_exval, radfit_itnum
    REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg
    REAL    (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3fit_cols, o3fit_dcols
    LOGICAL                                     :: do_skip_pix, do_cycle_this_pix
    LOGICAL                                     :: is_bad_pixel, yn_reference_fit
    INTEGER (KIND=i4), DIMENSION (4)            :: select_idx
    INTEGER (KIND=i4), DIMENSION (2)            :: exclud_idx
    INTEGER (KIND=i4)                           :: n_solar_pts
    REAL    (KIND=r8), DIMENSION (n_max_rspec)  :: solar_wvl
    real (kind=r8), dimension(:), allocatable :: adj_wvls, adj_spec, adj_wgts
    integer (kind=i4) :: adj_num, adj_num_allocated
    integer :: locerr

    ! CCM Array for holding fitted spectra
    REAL    (KIND=r8), DIMENSION (fitspc_out_dim0)   :: fitspc

    type (prefit_type) :: prefit

    !CHARACTER (LEN=28), PARAMETER :: modulename = 'xtrack_radiance_fitting_loop'

    if (errstat /= 0) return
    locerrstat = 0 ! pge_errstat_ok

    !!!fitvar_rad_saved = fitvar_rad_init

    adj_num_allocated = 0
    XTrackPix: DO ipix = first_pix, last_pix

      curr_xtrack_pixnum = ipix

      ! ---------------------------------------------------------------------
      ! If we already determined that this cross track pixel position carries
      ! an error, we don't even have to start processing.
      ! ---------------------------------------------------------------------
      IF ( omi_cross_track_skippix(ipix) .OR. szamax < omi_szenith(ipix,iloop) ) CYCLE

      !locerrstat = pge_errstat_ok

      n_database_wvl = n_omi_database_wvl(ipix)
      n_omi_radwvl = omi_nwav_rad (ipix,iloop)
      adj_num = n_omi_radwvl

      ! ---------------------------------------------------------------------------
      ! For each cross-track position we have to initialize the saved Shift&Squeeze
      ! ---------------------------------------------------------------------------
      saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

      ! ----------------------------------------------------------------------------
      ! Assign the solar wavelengths. Those should be current in the DATABASE array
      ! and can be taken from there no matter which case - YN_SOLAR_COMP and/or
      ! YN_RADIANCE_REFRENCE we are processing.
      ! ----------------------------------------------------------------------------
      n_solar_pts              = n_omi_database_wvl(ipix)
      if (n_solar_pts < 1) cycle  ! JED fix

      solar_wvl(1:n_solar_pts) = omi_database_wvl  (1:n_solar_pts, ipix)

      CALL check_wavelength_overlap ( &
        n_fitvar_rad,                                                &
        n_solar_pts,          solar_wvl (1:n_solar_pts),             &
        n_omi_radwvl, omi_radiance_wavl (1:n_omi_radwvl,ipix,iloop), &
        do_cycle_this_pix )

      IF (do_cycle_this_pix &
          .or. (n_database_wvl <= 0) &
          .or. (n_omi_radwvl <= 0) ) cycle

      ! ----------------------------------------------
      ! Restore DATABASE from OMI_DATABASE (see above)
      ! ----------------------------------------------
      database (1:n_database_wvl,1:max_rs_idx) = omi_database (1:n_database_wvl,ipix,1:max_rs_idx)

      ! ---------------------------------------------------------------------------------
      ! Restore solar fitting variables for across-track reference in Earthshine fitting.
      ! Note that, for the YN_SOLAR_COMP case, some variables have been assigned already
      ! in the XTRACK_RADIANCE_WAVCAL loop.
      ! ---------------------------------------------------------------------------------
      !sol_wav_avg                             = omi_radref_wav_avg(ipix)  ! JCH: no need to set this here
      Slit_Half_Width_1e                     = omi_solcal_pars(hwe_idx,ipix)
      Slit_Asym_Factor                       = omi_solcal_pars(asy_idx,ipix)
      Slit_Shape_Factor                       = omi_solcal_pars(sgk_idx,ipix)
      curr_sol_spec(1:n_database_wvl,wvl_idx) = omi_database_wvl(1:n_database_wvl,ipix)
      curr_sol_spec(1:n_database_wvl,spc_idx) = omi_database    (1:n_database_wvl,ipix,solar_idx)
      ! --------------------------------------------------------------------------------

      ! reallocate buffers if needed
      if (adj_num > adj_num_allocated) then
        if (adj_num_allocated > 0) then
          deallocate (adj_wvls, adj_spec, adj_wgts, stat=locerr)
          if (locerr /= 0) then
            call tell_error (tell_malloc_error, &
                 "xtrack_radiance_fitting_loop: deallocate failed", errstat)
            return
          endif
        endif
        allocate (adj_wvls(adj_num), adj_spec(adj_num), adj_wgts(adj_num), &
                  stat=locerr)
        if (locerr /= 0) then
          call tell_error (tell_malloc_error, &
               "xtrack_radiance_fitting_loop: allocate failed", errstat)
          return
        endif
        adj_num_allocated = adj_num
      endif

      adj_wvls(1:adj_num) = omi_radiance_wavl (1:adj_num, ipix, iloop)
      adj_spec(1:adj_num) = omi_radiance_spec (1:adj_num, ipix, iloop)
      adj_wgts(1:adj_num) = omi_radref_wght(1:adj_num,ipix)

      select_idx(1:4) = rad_ccdpix_selection(ipix,1:4)
      exclud_idx(1:2) = rad_ccdpix_exclusion(ipix,1:2)

      ! Set up generic fitting arrays
      CALL omi_adjust_radiance_data ( &
        select_idx(1:4), exclud_idx(1:2),                        &
        adj_num,                                            &
        adj_wvls(1:adj_num), adj_spec(1:adj_num), adj_wgts(1:adj_num), &
        omi_radiance_qflg  (1:n_omi_radwvl,ipix,iloop),          &
        omi_radiance_ccdpix(1:n_omi_radwvl,ipix,iloop),          &
        rad_spec_avg, do_skip_pix )

      n_rad_wvl = adj_num  !! FIXME: Gid rid of this global

      call copy_prefit_values (prefit, pge_idx, ipix, iloop)

      ! --------------------
      ! The radiance fitting
      ! --------------------
      fitcol       = r8_missval
      dfitcol      = r8_missval
      radfit_exval = INT(i2_missval, KIND=i4)
      radfit_itnum = INT(i2_missval, KIND=i4)
      rms          = r8_missval

      IF ((MAXVAL(adj_spec(1:adj_num)) > 0.0_r8) &
          .and. (adj_num > n_fitvar_rad) &
          .and. (.not. do_skip_pix)) then

        ! FIXME(?) JCH: in the original code, yn_reference_fit=.false.
        ! occurred _after_ the call to fit_radiance.  In other words,
        ! yn_reference_fit was uninitialized at the time it was passed
        ! to fit_radiance.  I assume the intent was to call fit_radiance
        ! with yn_reference_fit=.false.
        yn_reference_fit = .false.

        call tell_log (2, 'xtrack_radiance_fitting_loop: call fit_radiance')
        is_bad_pixel = .FALSE.
        CALL fit_radiance ( &
          pge_idx, ipix, ctrl_n_fitres_loop(radfit_idx), &
          ctrl_fitres_range(radfit_idx), &
          adj_num, adj_wvls, adj_spec, adj_wgts, &
          fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav, &
          prefit, o3fit_cols, o3fit_dcols, &
          target_var(1:n_fincol_idx,ipix), &
          allfit_cols(1:n_fitvar_rad,ipix), allfit_errs(1:n_fitvar_rad,ipix), &
          corr_matrix(1:n_fitvar_rad,ipix), is_bad_pixel, fitspc(1:adj_num), &
          yn_reference_fit, errstat)

        IF ( is_bad_pixel ) CYCLE

      END IF

      ! -----------------------------------
      ! Assign pixel values to final arrays
      ! -----------------------------------
      omi_fitconv_flag (ipix,iloop) = INT (radfit_exval, KIND=i2)
      omi_itnum_flag   (ipix,iloop) = INT (radfit_itnum, KIND=i2)
      omi_radfit_chisq (ipix,iloop) = chisquav
      omi_fit_rms      (ipix,iloop) = rms
      omi_column_amount(ipix,iloop) = fitcol
      omi_column_uncert(ipix,iloop) = dfitcol

      ! On occasion, the fit returns an exit value of 0 or 1 
      ! (suspicious or good) but a fitting uncertainty of zero. 
      ! In these cases, the SCD fit is almost always unphysical. 
      ! If this happens, force the fit convergence flag to "failed".
      IF ( (ABS(dfitcol) < 1e-15) .AND. (radfit_exval > -1) ) THEN
         omi_fitconv_flag (ipix,iloop) = INT (-2, KIND=i2)
      END IF

      ! CCM assign fit residual
      fitspc_out(1:adj_num,ipix,1) = fitspc(1:adj_num)
      fitspc_out(1:adj_num,ipix,2) = adj_spec(1:adj_num)
      fitspc_out(1:adj_num,ipix,3) = adj_wvls(1:adj_num)
      fitspc_out(1:adj_num,ipix,4) = adj_wgts(1:adj_num)

      IF ( pge_idx == pge_o3_idx ) THEN
        omi_o3_amount(o3_t1_idx:o3_t3_idx,ipix,iloop) = o3fit_cols (o3_t1_idx:o3_t3_idx)
        omi_o3_uncert(o3_t1_idx:o3_t3_idx,ipix,iloop) = o3fit_dcols(o3_t1_idx:o3_t3_idx)
      END IF

    END DO XTrackPix

    !errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_radiance_fitting_loop
END MODULE
