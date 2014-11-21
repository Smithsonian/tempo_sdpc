MODULE OMSAO_radiance_ref_module

  USE OMSAO_precision_module, ONLY: i2, i4, r4, r8
  USE OMSAO_parameters_module,   ONLY: MAX_STR_LEN
  USE OMSAO_variables_module, ONLY: n_rad_wvl_max
  use errormodule
  use terr_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC omi_adjust_radiance_data, remove_target_from_radiance, &
    omi_get_radiance_reference, xtrack_radiance_reference_loop

CONTAINS

  ! This function is called from omi_fitting with l1bfile=l1b_radref_filename.
  ! Among other things, this function sets
  !     rad_ccdpix_selection(*,1)=1 and  rad_ccdpix_selection(*,4)=nwrr.
  ! This is important to note since omi_read_radiance_lines uses those
  ! values to set omi_nwav_rad (producing omi_nwav_rad(*,*)=nwrr)
  ! It is this value that is used to set the number of wavelengths.
  !
  SUBROUTINE omi_get_radiance_reference (rpt_rr, &
                                         xtrange, &
                                         radwcal_lines, errstat )

    USE OMSAO_parameters_module, ONLY: &
      r4_missval, downweight, normweight, nlines_max
    USE OMSAO_indices_module,    ONLY: &
      qual_flag_mis, qual_flag_bad, qual_flag_err
    USE OMSAO_variables_module,  ONLY: ctrl_fit_winwav_lim, &
      ctrl_fit_winexc_lim, radiance_reference_lnums, radref_latrange, &
      Radiance_Paras_Type
    USE OMSAO_omidata_module,    ONLY: &
      rad_ccdpix_selection, omi_radiance_qflg, omi_radiance_spec, omi_radiance_wavl, &
      omi_szenith, omi_vzenith, omi_nwav_radref, omi_radref_spec, omi_radref_wavl,   &
      omi_radref_qflg, omi_radref_sza, omi_radref_vza, omi_radref_wght, &
      rad_ccdpix_exclusion, n_comm_wvl, &
      omi_nwav_rad, omi_radref_wav_avg
    USE OMSAO_errstat_module
    USE omi_pge_fitting_aux, ONLY: find_swathline_by_latitude, read_latitude
    USE omi_read_l1b_data, ONLY: omi_read_radiance_lines
    USE arrayutils, only: array_locate_r8
    USE errormodule
    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=*), PARAMETER :: modulename = 'omi_get_radiance_reference' ! JED fixed

    ! ---------------
    ! Input variables
    ! ---------------
    TYPE(Radiance_Paras_Type), INTENT(IN) :: rpt_rr
    INTEGER (KIND=i4), DIMENSION (0:rpt_rr%ntimes-1,2), INTENT(in) :: xtrange

    ! -----------------------------
    ! Output and Modified variables
    ! -----------------------------
    INTEGER (KIND=i4), DIMENSION (2), INTENT (OUT)   :: radwcal_lines
    INTEGER (KIND=i4),                INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL                      :: have_scanline
    LOGICAL, DIMENSION (2)       :: have_limits

    INTEGER (KIND=i4) :: fpix, lpix, midpt_line
    INTEGER (KIND=i4) :: nloop, j1, iline, ix, iloop, imin, imax, icnt
    REAL    (KIND=r4) :: lat_midpt
    REAL    (KIND=r8) :: specsum
    real    (kind=r4), dimension (:,:), allocatable :: latr4

    !INTEGER (KIND=i1), DIMENSION (0:rpt_rr%ntimes-1)       :: binfac
    !LOGICAL,           DIMENSION (0:rpt_rr%ntimes-1)       :: ynzoom

    REAL    (KIND=r4), DIMENSION (rpt_rr%nxtrack)           :: szacount
    REAL    (KIND=r8), DIMENSION (rpt_rr%nxtrack, rpt_rr%nwavel_ccd)     :: radref_spec, radref_wavl
    REAL    (KIND=r8), DIMENSION (rpt_rr%nwavel_ccd)           :: radref_wavl_ix
    REAL    (KIND=r8), DIMENSION (rpt_rr%nxtrack, rpt_rr%nwavel_ccd)     :: allcount, dumcount
    REAL    (KIND=r8), DIMENSION (rpt_rr%nwavel_ccd      )     :: cntr8
    integer :: locerrstat
    INTEGER (KIND=i4) :: ntrr, nxrr, nwrr
    integer (kind=i2) :: bad_qflg_mask
    real (kind=r8) :: sum_cntr8

    if (errstat < 0) return

    ! ------------------------------
    ! Initialize some some variables
    ! ------------------------------
    radiance_reference_lnums = -1  ! This will be written to file, hence needs a value
    lat_midpt = SUM ( radref_latrange ) / 2.0_r4
    ntrr = rpt_rr % ntimes
    nxrr = rpt_rr % nxtrack
    nwrr = rpt_rr % nwavel_ccd

    ALLOCATE (latr4(1:nxrr,0:ntrr-1), STAT=locerrstat)
    if (locerrstat /= 0) then
      errstat = -1
      call err_message_error ("omi_get_radiance_reference: allocate failed", errstat)
      return
    endif

    CALL read_latitude (rpt_rr%l1bfilename, rpt_rr%swathname, &
                        0, ntrr, latr4, errstat)
    if (errstat < 0) &
      return

    ! ----------------------------------------------------------------------
    ! Locate the swath line numbers corresponding the center of the latitude
    ! range to average into radiance reference spectrum.
    ! ----------------------------------------------------------------------
    call terr_trace (1, 'omi_get_radiance_reference:  calling '// &
                     'find_swathline_by_latitude (midpt_line)')
    CALL find_swathline_by_latitude ( &
      nxrr, 0, ntrr-1, latr4(1:nxrr,0:ntrr-1), lat_midpt, &
      xtrange(0:ntrr-1,1:2), midpt_line, have_scanline )

    ! --------------------------------------------------------------------
    ! If lower and upper bounds of the radiance reference block to average
    ! are identical, then we keep the midpoint line number as the only
    ! reference. Else locate the corresponding swath line numbers.
    ! --------------------------------------------------------------------
    IF ( radref_latrange(1) == radref_latrange(2) ) THEN
      radiance_reference_lnums(1:2) = midpt_line
      have_limits(1:2)           = .TRUE.
    ELSE
      call terr_trace (1, 'omi_get_radiance_reference:  calling '// &
                       'find_swathline_by_latitude (have_limits(1))')
      CALL find_swathline_by_latitude ( &
        nxrr, 0, midpt_line, latr4(1:nxrr,0:midpt_line), radref_latrange(1), &
        xtrange, radiance_reference_lnums(1), have_limits(1)   )
        !xtrange(0:midpt_line,1:2), radiance_reference_lnums(1), have_limits(1)   )
      call terr_trace (1, 'omi_get_radiance_reference:  calling '// &
                       'find_swathline_by_latitude (have_limits(2))')
      CALL find_swathline_by_latitude ( &
        nxrr, midpt_line, ntrr-1, latr4(1:nxrr,midpt_line:ntrr-1), radref_latrange(2), &
        xtrange, radiance_reference_lnums(2), have_limits(2) )
        !xtrange(midpt_line:ntrr-1,1:2), radiance_reference_lnums(2), have_limits(2) )
    END IF

    deallocate (latr4)

    ! -----------------------------------------------------
    ! If we don't find a working scan line, we have to fold
    ! -----------------------------------------------------
    IF ( ( .NOT. have_scanline )               .OR. &
        ( ANY ( .NOT. have_limits(1:2) ) )    .OR. &
        ( midpt_line < 0 )                       .OR. &
        ( ANY ( radiance_reference_lnums < 0 ) )        ) THEN
      CALL error_check ( 1, 0, pge_errstat_fatal, OMSAO_E_READ_L1B_FILE, &
                        modulename//f_sep//"Failed to find working radiance spectrum.", vb_lev_default, errstat )
      RETURN
    END IF

    ! ------------------------------------------------------------------
    ! A kludge for now, or maybe not. Assign the radiance reference
    ! swath line numbers to the radiance wavelength calibration swath
    ! line numbers.
    ! ------------------------------------------------------------------
    radwcal_lines(1:2) = radiance_reference_lnums(1:2)

    ! ---------------------------------------------------------------------------
    ! Fudge the selection of the complete spectrum; this is required because the
    ! routine that does the radiance line reading expects it to be available from
    ! after the Irradiance reading (which we are skipping). We will be setting up
    ! the proper numbers below.
    ! ---------------------------------------------------------------------------
    rad_ccdpix_selection(1:nxrr,1) = 1 ; rad_ccdpix_selection(1:nxrr,4) = nwrr

    ! --------------------------------------------------------------------
    ! Now we can average the spectra and the wavelength arrays. Loop over
    ! the block of swath lines in multiples of NLINES_MAX (100 by default)
    ! --------------------------------------------------------------------
    allcount    = 0.0_r8  ;  dumcount    = 0.0_r8  ;  szacount = 0.0_r4
    radref_wavl = 0.0_r8  ;  radref_spec = 0.0_r8
    omi_radref_sza = 0.0_r4 ; omi_radref_vza = 0.0_r4

    bad_qflg_mask = ior(qual_flag_mis, ior (qual_flag_bad, qual_flag_err))

    DO iline = radiance_reference_lnums(1), radiance_reference_lnums(2), nlines_max

      ! --------------------------------------------------------
      ! Check if loop ends before n_times_loop max is exhausted,
      ! or if we are outside the FIRST_LINE -> LAST_LINE range.
      ! --------------------------------------------------------
      nloop = MIN( nlines_max, radiance_reference_lnums(2)-radiance_reference_lnums(1)+1 )

      IF ( (iline+nloop) > ntrr  )  nloop = ntrr - iline

      ! ------------------------------
      ! Get NTIMES_LOOP radiance lines
      ! ------------------------------
      ! omi_read_radiance_lines also sets omi_nwav_rad
      write(*,*)'omi_get_radiance_reference calling omi_read_radiance_lines, iline=',iline
      CALL omi_read_radiance_lines (              &
        rpt_rr%l1bfilename, iline, nxrr, nloop, nwrr, errstat )

      ! Global used to set the dimension of fitspc
      n_rad_wvl_max = MAXVAL(omi_nwav_rad(:,0))

      DO iloop = 0, nloop-1

        ! ------------------------------------------------------
        ! Skip this cross-track position if there isn't any data
        ! ------------------------------------------------------
        fpix = xtrange(iline+iloop,1)
        lpix = xtrange(iline+iloop,2)

        DO ix = fpix, lpix

          cntr8(1:nwrr) = 1.0_r8

          where (iand(omi_radiance_qflg(1:nwrr,ix,iloop), bad_qflg_mask) /= 0)
            cntr8(1:nwrr) = 0.0_r8
          end where

          sum_cntr8 = sum (cntr8)
          if (sum_cntr8 > 0.0) then
          ! ------------------------------------
          ! Only proceed if we have a good value
          ! ------------------------------------

            omi_radiance_spec(1:nwrr,ix,iloop) = &
              omi_radiance_spec(1:nwrr,ix,iloop)*cntr8(1:nwrr) * cntr8(1:nwrr)

            !specsum = SUM ( omi_radiance_spec(1:nwrr,ix,iloop) ) / sum_cntr8
            !IF ( specsum == 0.0_r8 ) specsum = 1.0_r8
            ! Why do the above if specsum is set to 1.0 below?? --JED
            specsum = 1.0_r8

            radref_spec(ix,1:nwrr) = &
              radref_spec(ix,1:nwrr) + omi_radiance_spec(1:nwrr,ix,iloop)/specsum
            radref_wavl(ix,1:nwrr) = &
              radref_wavl(ix,1:nwrr) + omi_radiance_wavl(1:nwrr,ix,iloop)
            allcount(ix,1:nwrr) = allcount(ix,1:nwrr) + cntr8(1:nwrr)
            dumcount(ix,1:nwrr) = dumcount(ix,1:nwrr) + 1.0_r8

            IF ( omi_szenith(ix,iloop) /= r4_missval .AND. &
                omi_vzenith(ix,iloop) /= r4_missval         ) THEN
              omi_radref_sza(ix) = omi_radref_sza(ix) + omi_szenith(ix,iloop)
              omi_radref_vza(ix) = omi_radref_vza(ix) + omi_vzenith(ix,iloop)
              szacount      (ix) = szacount  (ix) + 1.0_r4
            END IF

          END IF

        END DO

      END DO

    END DO

    ! -----------------------------------------------------------
    ! Now for the actual averaging and assignment of final arrays
    ! -----------------------------------------------------------
    n_comm_wvl = 0
    DO ix = 1, nxrr

      ! -----------------------------------
      ! Average the wavelengths and spectra
      ! -----------------------------------
      WHERE ( allcount(ix,1:nwrr) /= 0.0_r8 )
        radref_spec(ix,1:nwrr) = radref_spec(ix,1:nwrr) / allcount(ix,1:nwrr)
      END WHERE
      WHERE ( dumcount(ix,1:nwrr) /= 0.0_r8 )
        radref_wavl(ix,1:nwrr) = radref_wavl(ix,1:nwrr) / dumcount(ix,1:nwrr)
      END WHERE

      ! -------------------------------------------
      ! Average the Solar and Viewing Zenith Angles
      ! -------------------------------------------
      IF ( szacount(ix) > 0.0_r4 ) THEN
        omi_radref_sza(ix) = omi_radref_sza(ix) / szacount(ix)
        omi_radref_vza(ix) = omi_radref_vza(ix) / szacount(ix)
      ELSE
        omi_radref_sza(ix) = r4_missval
        omi_radref_vza(ix) = r4_missval
      END IF

      ! -------------------------------------------------------------------------------
      ! Determine the CCD pixel numbers based on the selected wavelength fitting window
      ! -------------------------------------------------------------------------------

      radref_wavl_ix = radref_wavl (ix, 1:nwrr)
      DO j1 = 1, 3, 2
        CALL array_locate_r8 ( &
          nwrr, radref_wavl_ix, ctrl_fit_winwav_lim(j1  ), 'LE', &
          rad_ccdpix_selection(ix,j1  ) )
        CALL array_locate_r8 ( &
          nwrr, radref_wavl_ix, ctrl_fit_winwav_lim(j1+1), 'GE', &
          rad_ccdpix_selection(ix,j1+1) )
      END DO

      imin = rad_ccdpix_selection(ix,1)
      imax = rad_ccdpix_selection(ix,4)

      icnt = imax - imin + 1
      omi_nwav_radref(       ix) = icnt
      omi_radref_spec(1:icnt,ix) = radref_spec(ix,imin:imax)
      omi_radref_wavl(1:icnt,ix) = radref_wavl_ix(imin:imax)
      omi_radref_qflg(1:icnt,ix) = 0_i2
      omi_radref_wght(1:icnt,ix) = normweight

      ! -----------------------------------------------------------------
      ! Re-assign the average solar wavelength variable, sinfe from here
      ! on we are concerned with radiances.
      ! -----------------------------------------------------------------
      omi_radref_wav_avg(ix) = &
        SUM( omi_radref_wavl(1:icnt,ix) ) / REAL(icnt, KIND=r8)

      ! ------------------------------------------------------------------
      ! Set weights and quality flags to "bad" for missing spectral points
      ! ------------------------------------------------------------------
      allcount(ix,1:icnt) = allcount(ix,imin:imax)
      WHERE ( allcount(ix,1:icnt) == 0.0_r8 )
        omi_radref_qflg(1:icnt,ix) = 7_i2
        omi_radref_wght(1:icnt,ix) = downweight
      END WHERE

      ! ------------------------------------------------------------------------------
      ! If any window is excluded, find the corresponding indices. This has to be done
      ! after the array assignements above because we need to know which indices to
      ! exclude from the final arrays, not the complete ones read from the HE4 file.
      ! ------------------------------------------------------------------------------
      rad_ccdpix_exclusion(ix,1:2) = -1
      IF ( MINVAL(ctrl_fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r8 ( &
          nwrr, radref_wavl_ix, ctrl_fit_winexc_lim(1), 'GE', &
          rad_ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
          nwrr, radref_wavl_ix, ctrl_fit_winexc_lim(2), 'LE', &
          rad_ccdpix_exclusion(ix,2) )
      END IF

      ! ----------------------------------------
      ! Update the maximum number of wavelengths
      ! ----------------------------------------
      n_comm_wvl = MAX ( n_comm_wvl, icnt )

    END DO

    RETURN
  END SUBROUTINE omi_get_radiance_reference

  SUBROUTINE xtrack_radiance_reference_loop (&
      do_remove_target, nx, nw, &
      fpix, lpix, pge_idx, errstat )

    USE OMSAO_indices_module,    ONLY: &
      spc_idx, wvl_idx, &
      o3_t1_idx, o3_t3_idx, hwe_idx, asy_idx, shi_idx, &
      squ_idx, solar_idx, radref_idx, max_rs_idx, max_calfit_idx
    USE OMSAO_parameters_module, ONLY:  &
      i2_missval, r8_missval, downweight, normweight
    USE OMSAO_variables_module,  ONLY:  &
      database, curr_sol_spec, sol_wav_avg,                  &
      Slit_Half_Width_1e, Slit_Asym_Factor, n_fitvar_rad, verb_thresh_lev,  &
      n_database_wvl, fitvar_rad, n_fincol_idx, fincol_idx,                            &
      ctrl_n_fitres_loop, ctrl_fitres_range, xtrack_fitres_limit, &
      n_rad_wvl_max, target_npol, &
      curr_xtrack_pixnum, fitvar_rad_saved, fitvar_rad_init
    USE OMSAO_prefitcol_module, ONLY:  prefit_type, copy_prefit_values
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_omidata_module, ONLY: &
      omi_nwav_rad, n_omi_database_wvl, &
      omi_cross_track_skippix, n_omi_radwvl, &
      omi_database, omi_database_wvl, omi_solcal_pars,      &
      omi_radiance_wavl, omi_radref_wavl, omi_radiance_spec, omi_radref_spec,&
      omi_radiance_qflg, omi_radref_qflg, omi_radiance_spec,                 &
      rad_ccdpix_selection, omi_radiance_ccdpix, rad_ccdpix_exclusion,       &
      omi_radref_wght, omi_radref_pars,    &
      omi_radref_xflag, omi_radref_chisq, omi_radref_col,  &
      omi_radref_rms, omi_radref_dcol, omi_radref_xtrcol, omi_radref_wav_avg
    USE OMSAO_errstat_module
    USE radiance_fit, ONLY: fit_radiance
    use irradiance_data, only: Irr_Data
    use ctrlvars, only: yn_radiance_reference

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nx, nw, fpix, lpix
    LOGICAL,           INTENT (IN) :: do_remove_target

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ipix, radfit_exval, radfit_itnum
    REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg
    REAL    (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3fit_cols, o3fit_dcols
    REAL    (KIND=r8), DIMENSION (n_fitvar_rad)        :: corr_matrix_tmp, allfit_cols_tmp, allfit_errs_tmp
    LOGICAL                  :: do_skip_pix
    CHARACTER (LEN=MAX_STR_LEN) :: addmsg
    LOGICAL                                          :: is_bad_pixel
    INTEGER (KIND=i4), DIMENSION (4)                 :: select_idx
    INTEGER (KIND=i4), DIMENSION (2)                 :: exclud_idx
    INTEGER (KIND=i4)                                :: n_solar_pts
    REAL    (KIND=r8), DIMENSION (1:nw)              :: solar_wgt
    REAL    (KIND=r8), DIMENSION (n_fincol_idx,1:nx) :: target_var
    integer (kind=i4) :: nwav_irrad
    real (kind=r8), dimension(:), allocatable :: adj_wvls, adj_spec, adj_wgts
    integer (kind=i4) :: num_adj_allocated, n_rad_wvl_loc
    integer locerr

    !CHARACTER (LEN=30), PARAMETER :: modulename = 'xtrack_radiance_reference_loop'

    ! CCM fitted spectrum now returned from radiance_fit.f90
    !REAL    (KIND=r8), DIMENSION (n_rad_wvl)         :: fitspctmp
    ! The above will not work since size varies with the loop index --JED
    REAL    (KIND=r8), DIMENSION (n_rad_wvl_max) :: fitspctmp

    type (prefit_type) :: prefit

    if (errstat < 0) return

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    locerrstat          = pge_errstat_ok
    xtrack_fitres_limit = 0.0_r8
    target_var          = r8_missval
    fitvar_rad_saved    = fitvar_rad_init

    ! ---------------------------------------------------
    ! Note that this initialization will overwrite valid
    ! results on any second call to this subroutine. This
    ! happens, for example, when yn_radiance_reference
    ! and do_remove_target are selected simultaneously.
    ! In that case, however, we write the results to file
    ! before the second call.
    ! ---------------------------------------------------
    omi_radref_pars  (1:max_calfit_idx,1:nx) = r8_missval
    omi_radref_xflag (1:nx)                  = i2_missval
    omi_radref_chisq (1:nx)                  = r8_missval
    omi_radref_col   (1:nx)                  = r8_missval
    omi_radref_dcol  (1:nx)                  = r8_missval
    omi_radref_rms   (1:nx)                  = r8_missval
    omi_radref_xtrcol(1:nx)                  = r8_missval

    num_adj_allocated = 0

    XTrackPix: DO ipix = fpix, lpix

      locerrstat = pge_errstat_ok

      ! -----------------------------------------------------------
      ! The current cross-track pixel number is required further on
      ! when the slit function is computed: The CCD position based
      ! hyper-parameterization requires the knowledge of the row#.
      ! -----------------------------------------------------------
      curr_xtrack_pixnum = ipix

      ! ---------------------------------------------------------------------
      ! If we already determined that this cross track pixel position carries
      ! an error, we don't even have to start processing.
      ! ---------------------------------------------------------------------
      IF ( omi_cross_track_skippix(ipix) ) CYCLE

      n_database_wvl = n_omi_database_wvl(ipix)
      n_omi_radwvl   = omi_nwav_rad      (ipix,0)
      n_rad_wvl_loc = n_omi_radwvl

      ! ---------------------------------------------------------------------------
      ! For each cross-track position we have to initialize the saved Shift&Squeeze
      ! ---------------------------------------------------------------------------
      saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

      nwav_irrad = Irr_Data%nwaves(ipix)

      ! ------------------------------------------------------------------
      ! Assign number of irradiance wavelengths and the fitting weights
      ! from the solar wavelength calibration. Why? gga
      ! ------------------------------------------------------------------
      n_solar_pts              = nwav_irrad

      !       IF (yn_solar_comp) THEN !gga
      !          solar_wgt(1:n_solar_pts) = 1.0_r8
      !       ELSE
      !          solar_wgt(1:n_solar_pts) = omi_irradiance_wght(1:n_solar_pts,ipix)
      !       END IF !gga
      solar_wgt(1:n_solar_pts) = normweight

      ! -----------------------------------------------------
      ! Catch the possibility that N_OMI_RADWVL > N_SOLAR_PTS
      ! -----------------------------------------------------
      IF ( n_omi_radwvl > n_solar_pts ) THEN
        solar_wgt(n_solar_pts+1:n_omi_radwvl) = downweight
      END IF

      ! tpk: Should this be "> n_fitvar_rad"?
      if ((n_database_wvl <= 0) .or. (n_rad_wvl_loc <= 0)) cycle

      ! ----------------------------------------------
      ! Restore DATABASE from OMI_DATABASE (see above)
      ! ----------------------------------------------
      !
      ! Note: In xtrack_radiance_wvl_calibration, the database array was
      !       computed and then assigned to omi_database.  Here, it is
      !       restored to for use in radiance fitting.  --JED
      database (1:n_database_wvl,1:max_rs_idx) = omi_database (1:n_database_wvl,ipix,1:max_rs_idx)

      ! -----------------------------------------------------------------------
      ! Restore solar fitting variables for across-track reference in Earthshine fitting
      ! --------------------------------------------------------------------------------
      sol_wav_avg = omi_radref_wav_avg (ipix)
      Slit_Half_Width_1e = omi_solcal_pars(hwe_idx,ipix)
      Slit_Asym_Factor = omi_solcal_pars(asy_idx,ipix)
      curr_sol_spec(1:n_database_wvl,wvl_idx) = omi_database_wvl(1:n_database_wvl,ipix)
      curr_sol_spec(1:n_database_wvl,spc_idx) = omi_database    (1:n_database_wvl,ipix,solar_idx)
      ! --------------------------------------------------------------------------------

      ! ---------------------------------------------------------------
      ! If a Radiance Reference is being used, then it must be calibrated
      ! rather than the swath line that has been read.
      ! ---------------------------------------------------------------
      IF ( yn_radiance_reference ) THEN
        omi_radiance_wavl(1:n_omi_radwvl,ipix,0) = omi_radref_wavl(1:n_omi_radwvl,ipix)
        omi_radiance_spec(1:n_omi_radwvl,ipix,0) = omi_radref_spec(1:n_omi_radwvl,ipix)
        omi_radiance_qflg(1:n_omi_radwvl,ipix,0) = omi_radref_qflg(1:n_omi_radwvl,ipix)
      END IF

      ! -------------------------------------------------------------------------

      ! reallocate buffers if needed
      if (n_rad_wvl_loc > num_adj_allocated) then
        if (num_adj_allocated > 0) then
          deallocate (adj_wvls, adj_spec, adj_wgts)
        endif
        allocate (adj_wvls(n_rad_wvl_loc), adj_spec(n_rad_wvl_loc), adj_wgts(n_rad_wvl_loc), &
                  stat=locerr)
        if (locerr /= 0) then
          call err_message_error ("xtrack_solar_calibration_loop: allocate failed", errstat)
          return
        endif
        num_adj_allocated = n_rad_wvl_loc
      endif

      adj_wvls(1:n_rad_wvl_loc) = omi_radiance_wavl (1:n_rad_wvl_loc, ipix, 0)
      adj_spec(1:n_rad_wvl_loc) = omi_radiance_spec (1:n_rad_wvl_loc, ipix, 0)
      adj_wgts(1:n_rad_wvl_loc) = solar_wgt(1:n_rad_wvl_loc)

      select_idx(1:4) = rad_ccdpix_selection(ipix,1:4)
      exclud_idx(1:2) = rad_ccdpix_exclusion(ipix,1:2)

      ! Set up generic fitting arrays
      CALL omi_adjust_radiance_data ( & 
        select_idx(1:4), exclud_idx(1:2),            &
        n_rad_wvl_loc, &
        adj_wvls(1:n_rad_wvl_loc), adj_spec(1:n_rad_wvl_loc), adj_wgts(1:n_rad_wvl_loc), &
        omi_radiance_qflg   (1:n_omi_radwvl,ipix,0), &
        omi_radiance_ccdpix (1:n_omi_radwvl,ipix,0), &
        rad_spec_avg, do_skip_pix )

      ! --------------------------------------------------------------------
      ! Update the weights for the Reference/Wavelength Calibration Radiance
      ! --------------------------------------------------------------------
      omi_radref_wght(1:n_rad_wvl_loc,ipix) = adj_wgts(1:n_rad_wvl_loc)

      IF (.NOT. yn_radiance_reference) THEN
        call copy_prefit_values (prefit, pge_idx, ipix, 0)
      ELSE
        prefit%o3_col = 0.0_r8
        prefit%o3_dcol = 0.0_r8
      END IF

      ! --------------------
      ! The radiance fitting
      ! --------------------
      fitcol       = r8_missval
      dfitcol      = r8_missval
      radfit_exval = INT(i2_missval, KIND=i4)
      radfit_itnum = INT(i2_missval, KIND=i4)
      rms          = r8_missval

      addmsg = ''
      IF ((MAXVAL(adj_spec(1:n_rad_wvl_loc)) > 0.0_r8) &
          .AND. (n_rad_wvl_loc > n_fitvar_rad) &
          .AND. (.NOT. do_skip_pix)) THEN
        is_bad_pixel     = .FALSE.

        call terr_trace (2, 'OMSAO_radiance_ref_module: call fit_radiance')
        CALL fit_radiance ( &
          pge_idx, ipix, ctrl_n_fitres_loop(radref_idx), &
          ctrl_fitres_range(radref_idx), &
          n_rad_wvl_loc, adj_wvls, adj_spec, adj_wgts, &
          fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,               &
          prefit, o3fit_cols, o3fit_dcols,                                          &
          target_var(1:n_fincol_idx,ipix),                                          &
          allfit_cols_tmp(1:n_fitvar_rad), allfit_errs_tmp(1:n_fitvar_rad),         &
          corr_matrix_tmp(1:n_fitvar_rad), is_bad_pixel, fitspctmp, &
          errstat)

        IF ( is_bad_pixel ) CYCLE

        WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'RADIANCE Reference #', ipix, &
          ': hw 1/e = ', Slit_Half_Width_1e, '; e_asy = ', Slit_Asym_Factor, '; shift = ', &
          fitvar_rad(shi_idx), '; squeeze = ', fitvar_rad(squ_idx),&
          '; exit val = ', radfit_exval, '; iter num = ', radfit_itnum
      ELSE
        WRITE (addmsg, '(A,I2,A)') 'RADIANCE Reference #', ipix, ': Skipped!'
      END IF

      ! ------------------
      ! Report on progress
      ! ------------------
      CALL error_check ( &
        0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), vb_lev_omidebug, errstat )
      IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

      ! -----------------------------------
      ! Assign pixel values to final arrays
      ! -----------------------------------
      omi_radref_pars (1:max_calfit_idx,ipix) = fitvar_rad(1:max_calfit_idx)
      omi_radref_xflag(ipix)                  = INT (radfit_exval, KIND=i2)
      omi_radref_chisq(ipix)                  = chisquav
      omi_radref_col  (ipix)                  = fitcol
      omi_radref_dcol (ipix)                  = dfitcol
      omi_radref_rms  (ipix)                  = rms

      ! -------------------------------------------------------------------------
      ! Remember weights for the reference radiance, to be used as starting point
      ! in the regular radiance fitting
      ! -------------------------------------------------------------------------
      omi_radref_wght(1:n_rad_wvl_loc,ipix) = adj_wgts(1:n_rad_wvl_loc)

      ! -----------------------------------------------
      ! Update the solar spectrum entry in OMI_DATABASE
      ! -----------------------------------------------
      IF ( yn_radiance_reference ) &
        omi_database (1:n_rad_wvl_loc,ipix,solar_idx) = omi_radref_spec(1:n_rad_wvl_loc,ipix)

    END DO XTrackPix

    ! -----------------------------------------
    ! Remove target gas from radiance reference
    ! -----------------------------------------

    IF ( yn_radiance_reference .AND. do_remove_target ) THEN
      ! ----------------------------------------------------------------
      ! Removing the target gas from the radiance reference will alter
      ! OMI_RADREF_SPEC (1:NWVL,FPIX:LPIX). This is being passed to the
      ! subroutine via MODULE use rather than through the argument list.
      ! ----------------------------------------------------------------
      CALL remove_target_from_radiance (                                  &
        nw, fpix, lpix, n_fincol_idx, fincol_idx(1:2,1:n_fincol_idx),  &
        target_npol, target_var(1:n_fincol_idx,fpix:lpix), omi_radref_xtrcol(fpix:lpix) )

    END IF

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_radiance_reference_loop

  SUBROUTINE remove_target_from_radiance (       &
      nw, ipix, jpix, n_loc_fincol_idx, loc_fincol_idx, &
      target_npol, target_var, target_fit         )

    ! USE OMSAO_indices_module,   ONLY: solar_idx
    USE OMSAO_parameters_module,ONLY: downweight, r8_missval
    USE OMSAO_omidata_module,   ONLY: omi_radref_spec, omi_database, omi_nwav_radref
    USE OMSAO_variables_module, ONLY: refspecs_original
    USE OMSAO_median_module,    ONLY: median
    USE SLATEC_davint, ONLY: dpolft

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4),                                     INTENT (IN) :: &
      nw, ipix, jpix, n_loc_fincol_idx, target_npol
    INTEGER (KIND=i4), DIMENSION (2,n_loc_fincol_idx),         INTENT (IN) :: loc_fincol_idx

    ! ------------------
    ! Modified Variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (n_loc_fincol_idx,ipix:jpix), INTENT (INOUT) :: target_var
    REAL (KIND=r8), DIMENSION              (ipix:jpix), INTENT (OUT)   :: target_fit

    ! ------------------------------
    ! Local Variables and Parameters
    ! ------------------------------
    INTEGER (KIND=i4)                 :: i, j, k, l, nwvl
    REAL    (KIND=r8)                 :: yfloc
    REAL    (KIND=r8), DIMENSION (nw) :: tmpexp

    !CHARACTER (LEN=27), PARAMETER :: modulename = 'remove_target_from_radiance'

    ! ----------------
    ! DPOLFt variables
    ! ----------------
    INTEGER (KIND=i4)            :: ndeg, ierr, nx, npol, nfit
    REAL    (KIND=r8)            :: eps
    REAL    (KIND=r8), DIMENSION (jpix-ipix+1) :: x, y, yf, w
    REAL    (KIND=r8), DIMENSION (3*((jpix-ipix+1)+target_npol+1)) :: a

    target_fit = 0.0_r8

    npol = target_npol
    nx   = jpix-ipix+1
    DO j = 1, n_loc_fincol_idx
      k = loc_fincol_idx(2,j)

      ! ----------------------------------------------------------------------
      ! If we can/have to, fit a cross-track polynomial to the fitted columns,
      ! we do this individually for each loc_fincol_idx and remove the smoothed
      ! column loading rather than the originally fitted one. In any case, YF
      ! will contain the column values to be removed. Hence the outer loop
      ! over N_loc_fincol_idx rather than cross-track position.
      ! ----------------------------------------------------------------------
      nfit = 0
      DO i = ipix, jpix
        IF ( target_var(j,i) > r8_missval ) THEN
          nfit    = nfit + 1
          y(nfit) = target_var(j,i)
          w(nfit) = 1.0_r8
        END IF
      END DO

      IF ( nfit /= (jpix-ipix+1) .AND. nfit > 0 ) THEN
        WHERE ( target_var(j,ipix:jpix) <= r8_missval )
          target_var(j,ipix:jpix) = SUM(y(1:nfit))/REAL(nfit,KIND=r8)
        ENDWHERE
      END IF

      ! ----------------------------------------------------
      ! We either fit a polynomial or simply use the Median.
      ! The distinction is made depending on
      !
      ! (a) the order of the cross-track polynomial, and
      ! (b) the number of cross-track points we can fit.
      ! ----------------------------------------------------
      IF ( npol > 0 .AND. jpix-ipix+1 > npol ) THEN

        !IF ( npol >=0 .AND. nfit > npol ) THEN
        !eps  = -0.1_r8  ! Chose the best-fitting order

        eps =  0.0_r8  ! Fit the complete NPOL polynomial
        ndeg = npol
        x(1:nx) = (/ ( REAL(i-nx/2, KIND=r8), i = 1, nx ) /) / REAL(nx/2, KIND=r8)
        y(1:nx) = target_var(j,ipix:jpix)
        WHERE ( y(1:nx) > r8_missval )
          w(1:nx) = 1.0_r8
        ELSEWHERE
          w(1:nx) = downweight
        END WHERE
        CALL dpolft (&
          nx, x(1:nx), y(1:nx), w(1:nx), npol, ndeg, eps, yf(1:nx), ierr, a )

        !CALL dpolft (&
        !     nfit, x(1:fit), y(1:nfit), w(1:nfit), npol, ndeg, eps, yf(1:nfit), ierr, a )

      ELSE
        ! -----------------------------------------------------------------
        ! The Median is a better choice than the Mean, since the former
        ! is less sensitive to outliers. The Mean may be skewed towards
        ! abnormally high values at the edges of the swath.
        ! -----------------------------------------------------------------
        nfit = nx
        yf(1:nx) = median(nx, target_var(j,ipix:jpix))
      END IF

      DO i = ipix, jpix

        l = i - ipix + 1

        nwvl = omi_nwav_radref(i)

        IF ( yf(l) > r8_missval ) THEN
          yfloc = yf(l)
        ELSE
          yfloc = 0.0_r8
        END IF

        tmpexp(1:nwvl) = yfloc * omi_database(1:nwvl,i,k)
        WHERE ( tmpexp >= MAXEXPONENT(1.0_r8) )
          tmpexp = MAXEXPONENT(1.0_r8) - 1.0_r8
        ENDWHERE
        WHERE ( tmpexp <= MINEXPONENT(1.0_r8) )
          tmpexp = MINEXPONENT(1.0_r8) + 1.0_r8
        ENDWHERE

        omi_radref_spec(1:nwvl,i) = omi_radref_spec(1:nwvl,i) * &
          EXP(+tmpexp(1:nwvl))

        target_fit(i) = target_fit(i) + yfloc/refspecs_original(k)%NormFactor

      END DO
    END DO

    RETURN
  END SUBROUTINE remove_target_from_radiance

  SUBROUTINE omi_adjust_radiance_data ( &
      omi_ccdpix_idx, omi_ccdpix_exc, n_adj, &
      adj_wvls, adj_spec, adj_wgts, &
      omi_rad_qflg, &
      omi_rad_ccd, &
      rad_spec_avg, do_skip_pix )

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: &
      qual_flag_mis, qual_flag_bad, qual_flag_err
    USE OMSAO_parameters_module,    ONLY: downweight, r4_missval
    use ctrlvars, only: yn_radiance_reference, yn_spectrum_norm, yn_solar_comp
    USE OMSAO_solcomp_module,       ONLY: solarcomp_pars
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                             INTENT (IN) :: n_adj
    INTEGER (KIND=i4), DIMENSION (2),              INTENT (IN) :: omi_ccdpix_exc
    INTEGER (KIND=i4), DIMENSION (4),              INTENT (IN) :: omi_ccdpix_idx
    INTEGER (KIND=i2), DIMENSION (n_adj),   INTENT (IN) :: omi_rad_qflg

    ! ----------------
    ! Output variables
    ! ----------------
    LOGICAL,                                                INTENT (OUT) :: do_skip_pix
    REAL    (KIND=r8),                                      INTENT (OUT) :: rad_spec_avg
    INTEGER (KIND=i4),  DIMENSION (n_adj),           INTENT (OUT) :: omi_rad_ccd
    !REAL    (KIND=r8),  DIMENSION (ccd_idx,1:n_adj), INTENT (OUT) :: curr_rad_spec
    real(kind=r8), dimension(n_adj), intent(inout) :: adj_wvls, adj_spec, adj_wgts

    ! ------------------
    ! Modified variables
    ! ------------------
    !INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                     :: &
      i, locerrstat, imin1, imax1, imin2, imax2, j1, j2
    LOGICAL                                               :: have_good_window
    REAL    (KIND=r8), DIMENSION (n_adj)           :: weightsum
    integer (kind=i2) :: bad_qflg_mask

    locerrstat  = pge_errstat_ok
    do_skip_pix = .FALSE.

    imin1 = omi_ccdpix_idx(1) ; imax1 = omi_ccdpix_idx(4)  ! The total window
    imin2 = omi_ccdpix_idx(2) ; imax2 = omi_ccdpix_idx(3)  ! The fitting window
    omi_rad_ccd (1:n_adj) = (/ (i, i = imin1, imax1) /)

    ! ---------------------------------------------------------
    ! Compute the weights. This is a bit tedious, as we have to
    ! check for a number of things that can go wrong. We start
    ! out assuming "all is well" and exclude/modify only those
    ! entries that are expected to give us trouble.
    ! ---------------------------------------------------------
    ! EXCEPT, of course, that we also want to exclude any CCD
    ! pixel already excluded from the solar fit. So instead of
    ! starting out with "all is well" we start out with
    ! "anything that is well in the solar spectrum".
    !
    ! We assume in this that the number of spectral points is
    ! equal in both spectra. A reminder to ourselves: Check
    ! that this indeed the case!
    ! ---------------------------------------------------------

    ! ---------------------------------------
    ! (1) Anything outside the fitting window
    ! ---------------------------------------
    ! (the CCD indices are absolute positions, i.e., unlikely to be "1:n_sol_wvl")
    ! ----------------------------------------------------------------------------
    IF ( imin2 > imin1 ) adj_wgts(1:imin2-imin1+1) = downweight
    IF ( imax2 < imax1 ) adj_wgts(imax2-imin1+1:n_adj) = downweight

    ! ----------------------------------------------------------------------
    ! (2) Any window excluded by the user (specified in fitting control file
    ! ----------------------------------------------------------------------
    IF ( ALL ( omi_ccdpix_exc(1:2) > 0 ) ) THEN
      j1 = omi_ccdpix_exc(1) - imin1 + 1 ; j2 = omi_ccdpix_exc(2) - imin1 + 1
      IF ( j1 >= 1 .AND. j2 <= n_adj ) &
        adj_wgts(j1:j2) = downweight
    END IF

    ! ----------------------------------
    ! (3) Wavelengths in ascending order
    ! ----------------------------------
    DO i = 2, n_adj
      IF ( adj_wvls(i) <= adj_wvls(i-1) ) THEN
        adj_wvls(i) = adj_wvls(i-1) + 0.001_r8
        adj_wgts(i) = downweight
      END IF
    END DO

    ! ---------------------------------
    ! (4) No missing values in spectrum
    ! ---------------------------------
    WHERE ( adj_spec(1:n_adj) <= REAL( r4_missval, KIND=r8 ) )
      adj_wgts(1:n_adj) = downweight
      adj_spec(1:n_adj) = 0.0_r8
    END WHERE

    ! -------------------------------------------------------------------------------
    ! Translate window limit wavelenghts into indices; making sure that Shift&Squeeze
    ! doesn't shift the wavelength array off the limits we set at the beginning. Do
    ! this iteratively until we have found good window limits. In the best case, all
    ! is well the first time around, but we might have to adjust the window margins
    ! if we fail to read all the data.
    ! -------------------------------------------------------------------------------
    have_good_window = .TRUE.

    ! --------------------------------------------------------------------
    ! Find the pixel quality flags.
    ! Choice of flags is based on the recommendations of the L1b README.
    ! --------------------------------------------------------------------
    bad_qflg_mask = ior(qual_flag_mis, ior (qual_flag_bad, qual_flag_err))
    WHERE (iand (omi_rad_qflg(1:n_adj), bad_qflg_mask) /= 0)
      adj_wgts(1:n_adj) = downweight
      adj_spec(1:n_adj) = 0.0_r8
    END WHERE

    ! --------------------------------------------------
    ! Compute normalization factor for radiance spectrum
    ! --------------------------------------------------
    weightsum = 0.0_r8
    WHERE ( adj_wgts(1:n_adj) /= downweight )
      weightsum = 1.0_r8
    END WHERE

    rad_spec_avg = SUM ( &
      ABS(adj_spec(1:n_adj))*weightsum(1:n_adj) ) / &
      MAX(1.0_r8, SUM(weightsum(1:n_adj)))

    IF ( rad_spec_avg == 0.0_r8 ) rad_spec_avg = 1.0_r8

    ! -------------------------------------------------------------------------
    ! So far we have only taken care of/excluded any negative values in the
    ! spectrum, but there may abnormally high or low positive values also.
    ! Now we check for any values exceeding 100 times the average, which should
    ! be a large enough window to keep anything sensible and reject the real
    ! outliers.
    ! -------------------------------------------------------------------------
    WHERE ( weightsum(1:n_adj) /= 0.0_r8 .AND. &
           ABS(adj_spec(1:n_adj)) >= 100.0_r8 * rad_spec_avg )
      weightsum(1:n_adj) = 0.0_r8
      adj_wgts(1:n_adj) = downweight
      adj_spec(1:n_adj) = 0.0_r8
    ENDWHERE

    ! ---------------------------------------------------------------------
    ! Recompute the radiance spectrum average, because it may have changed.
    ! ---------------------------------------------------------------------
    rad_spec_avg = SUM ( &
      adj_spec(1:n_adj)*weightsum(1:n_adj) ) / &
      MAX(1.0_r8, SUM(weightsum(1:n_adj)))
    IF ( rad_spec_avg <= 0.0_r8 ) THEN
      do_skip_pix = .TRUE.
      rad_spec_avg = 1.0_r8
    ELSE
      ! -----------------------------------------
      ! Finally, normalize the radiance spectrum.
      ! -----------------------------------------
      ! There are three possibilities:
      ! (1) If YN_SPECTRUM_NORM = .TRUE. the spectrum will be normalized to 1.
      ! (2) If YN_SPECTRUM_NORM = .FALSE. but YN_SOLAR_COMP = .TRUE. and
      !     YN_RADIANCE_REFERENCE = .FALSE. then the radiance spectrum will be
      !      normalized with the composite solar norm.
      ! (3) If both of the above are .FALSE. do nothing.
      !
      ! "(2)" assures that radiance and irradiance retain their relative
      ! magnitudes, which means that the Solar Intensity parameter will be more
      ! closely associated with (but not necessarily identical to) the scene albedo.
      ! -----------------------------------------------------------------------------
      IF ( .NOT. yn_spectrum_norm ) THEN                               ! branch for "(2)" or "(3)"
        IF ( yn_solar_comp .AND. (.NOT. yn_radiance_reference) ) THEN ! "(2)"
          rad_spec_avg = solarcomp_pars%SolarNorm
        ELSE                                                          ! "(3)"
          rad_spec_avg = 1.0_r8
        END IF
      END IF
      adj_spec(1:n_adj) = adj_spec(1:n_adj) / rad_spec_avg
    END IF

    RETURN
  END SUBROUTINE omi_adjust_radiance_data

END MODULE OMSAO_radiance_ref_module
