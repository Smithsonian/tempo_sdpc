MODULE omi_pge_fitting_process

  use errormodule
  use terr_module
  use l1bread_utils
  private
  public omi_pge_fitting

CONTAINS
  
SUBROUTINE omi_pge_fitting ( pge_idx, n_max_rspec, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,      ONLY: pge_errstat_ok, pge_errstat_error, pge_errstat_fatal
  USE OMSAO_he5_module,          ONLY: NrofScanLines, NrofCrossTrackPixels
  USE OMSAO_variables_module,    ONLY: l1b_rad_filename, Radiance_Paras_Type, &
    l1b_radref_filename, l1b_channel
  use ctrlvars, only: yn_radiance_reference
  USE OMSAO_omidata_module,      ONLY: omi_radiance_swathname, EarthSunDistance
  USE omi_pge_fitting_aux, ONLY: omi_set_fitting_parameters
  USE omi_read_l1b_data, ONLY: L1Bga_EarthSunDistance
  !use l1bread, only: l1bread_radiance_info
  USE OMSAO_errstat_module
  USE OMSAO_solcomp_module, ONLY: soco_pars_deallocate
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: errstat
  TYPE(Radiance_Paras_Type) :: rpt_rad, rpt_rr

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_pge_fitting_process'

  pge_error_status = pge_errstat_ok

  ! -------------------------------------------------------------------------------------
  ! Set the swath name of various ESDTs
  ! -------------------------------------------------------------------------------------
  CALL omi_set_fitting_parameters ( pge_idx, errstat )
  ! -------------------------------------------------------------------------------------
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error ) GO TO 666

  ! -----------------------------------------------------------------------------------
  ! Get dimensions the L1B radiance granule
  ! -----------------------------------------------------------------------------------
  errstat = pge_errstat_ok
  !CALL l1bread_radiance_info (l1b_rad_filename, l1b_channel, rpt_rad, errstat)
  call read_l1_radiance_info (l1b_rad_filename, l1b_channel, rpt_rad, errstat)
  if (errstat < 0) goto 666

  EarthSunDistance = L1Bga_EarthSunDistance( l1b_rad_filename, rpt_rad%swathname )
  omi_radiance_swathname = rpt_rad%swathname

  NrofScanLines        = rpt_rad%ntimes
  NrofCrossTrackPixels = rpt_rad%nxtrack

  ! ---------------------------------------------------------------
  ! Dimensions for Radiance Reference granule
  ! ---------------------------------------------------------------
  IF ( .NOT. yn_radiance_reference ) THEN
    l1b_radref_filename = l1b_rad_filename
    rpt_rr%l1bfilename = rpt_rad%l1bfilename
    rpt_rr%ntimes = rpt_rad%ntimes
    rpt_rr%nxtrack = rpt_rad%nxtrack
    rpt_rr%nwavel_ccd = rpt_rad%nwavel_ccd
    rpt_rr%swathname = rpt_rad%swathname
  ELSE
    !CALL l1bread_radiance_info (l1b_radref_filename, l1b_channel, rpt_rr, errstat)
    call read_l1_radiance_info (l1b_radref_filename, l1b_channel, rpt_rr, errstat)
    if (errstat < 0) goto 666
  ENDIF

  ! ----------------------------------------------------------------
  ! Number of cross-track positions must be the same; fold otherwise
  ! ----------------------------------------------------------------
  IF ( rpt_rad%nxtrack /= rpt_rr%nxtrack ) THEN
    CALL error_check ( 0, 1, pge_errstat_fatal, OMSAO_F_XTRMISRAD, &
                      modulename, vb_lev_default, errstat )
    GO TO 666
  END IF

  ! -----------------------------------------------------------------------------------
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ! -----------------------------------------------------------------
  CALL omi_fitting (pge_idx, rpt_rad, rpt_rr, n_max_rspec, pge_error_status)

  IF ( pge_error_status >= pge_errstat_fatal ) GO TO 666

  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
 666 CONTINUE

  ! -------------------------------------------------
  ! Deallocation of some potentially allocated memory
  ! -------------------------------------------------
  CALL soco_pars_deallocate (errstat)

  IF ( pge_error_status >= pge_errstat_fatal ) RETURN

  RETURN
END SUBROUTINE omi_pge_fitting

SUBROUTINE omi_fitting (pge_idx, rpt_rad, rpt_rr, &
                        n_max_rspec, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: i2_missval, MAX_STR_LEN
  USE OMSAO_indices_module,    ONLY: sao_molecule_names
  USE OMSAO_variables_module,  ONLY: &
    l1b_rad_filename, &
    l2_filename, pixnum_lim,    &
    radfit_latrange,                &
    common_latrange,    &
    Radiance_Paras_Type, &
    radiance_reference_lnums, l1b_radref_filename
  use ctrlvars, only: yn_radiance_reference, yn_common_iter, &
    yn_diagnostic_run, yn_remove_target
  USE OMSAO_he5_module,       ONLY:  pge_swath_name
  USE OMSAO_solar_wavcal_module, ONLY: xtrack_solar_calibration_loop
  USE OMSAO_radiance_ref_module, ONLY: omi_get_radiance_reference, &
    xtrack_radiance_reference_loop
  USE OMSAO_prefitcol_module, ONLY: read_prefit_columns, init_prefit_files
  USE OMSAO_errstat_module
  USE OMSAO_wfamf_module, ONLY: omi_read_climatology, CmETA
  USE he5_output_tools, ONLY: he5_init_swath, he5_define_fields, &
    he5_close_output_file, he5_set_field_attributes, &
    he5_write_global_attributes, he5_write_swath_attributes, &
    he5_write_wavcal_output, he5_write_common_mode !, he5_open_readwrite
  USE omi_read_l1b_data, ONLY: omi_read_binning_factor, &
    omi_read_radiance_lines, omi_read_radiance_lines
  USE omi_pge_fitting_aux, ONLY: omi_set_xtrpix_range, &
    read_latitude, find_swathline_range
  use commonmode, only: finalize_common_mode
  USE fitting_loops, ONLY: xtrack_radiance_wvl_calibration
  USE metadata_tools, ONLY: check_metadata_consistency, set_l2_metadata
  USE omi_pge_postprocessing, ONLY: omi_pge_postprocess
  USE swathline_loop, ONLY: swathline_loops
  use datafields, only: he5_initialize_datafields
  USE OMSAO_omidata_module, ONLY: n_comm_wvl, ntimes_loop, &
    omi_cross_track_skippix, omi_radcal_xflag, &
    omi_radiance_swathname
  USE irradiance_data, only: irradiance_data_init

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec
  TYPE(Radiance_Paras_Type), INTENT(IN) :: rpt_rad, rpt_rr

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! -------------------------
  ! Local variables (for now)
  ! -------------------------
  INTEGER   (KIND=i4) ::                                                 &
    iline, first_line, last_line, locerrstat, first_wc_pix, last_wc_pix, &
    first_pix, last_pix
  INTEGER (kind=i4) :: ntimes_rad, nxtrack_rad, nwavel_rad
  INTEGER (kind=i4) :: ntimes_rr, nxtrack_rr, nwavel_rr
  INTEGER (KIND=i4), DIMENSION (2) :: radiance_wavcal_lnums

  ! ----------------------------------------------------------------------
  ! Swath dimensions and variables that aren't passed from calling routine
  ! ----------------------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: molname

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! and Common Mode spectrum
  ! ----------------------------------------------------------
  INTEGER (KIND=i1), DIMENSION (0:rpt_rad%ntimes-1)   :: omi_binfac
  INTEGER (KIND=i4), DIMENSION (0:rpt_rad%ntimes-1,2) :: omi_xtrpix_range
  LOGICAL,           DIMENSION (0:rpt_rad%ntimes-1)   :: &
    omi_is_szoom, is_common_range, do_radfit_range

  INTEGER (KIND=i1), DIMENSION (0:rpt_rr%ntimes-1)   :: omi_binfac_rr
  INTEGER (KIND=i4), DIMENSION (0:rpt_rr%ntimes-1,2) :: omi_xtrpix_range_rr
  LOGICAL,           DIMENSION (0:rpt_rr%ntimes-1)   :: omi_is_szoom_rr

  ! ----------------------------------------------------------
  ! OMI L1b latitudes
  ! ----------------------------------------------------------
  REAL (KIND=r4), DIMENSION (:,:), allocatable :: l1b_rad_latitudes

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'omi_fitting'

  ! ------------------
  ! External functions
  ! ------------------
  !INTEGER (KIND=i4), EXTERNAL :: &
  !  he5_init_swath, he5_define_fields, he5_close_output_file, &
  !  he5_set_field_attributes, he5_write_global_attributes,    &
  !  he5_write_swath_attributes, he5_open_readwrite

  if (errstat < 0) return

  ntimes_rad = rpt_rad%ntimes
  nxtrack_rad = rpt_rad%nxtrack
  nwavel_rad = rpt_rad%nwavel_ccd

  ntimes_rr = rpt_rr%ntimes
  nxtrack_rr = rpt_rr%nxtrack
  nwavel_rr = rpt_rr%nwavel_ccd

  ! ---------------------------------------------------------------
  ! Some initializations that will save us headaches in cases where
  ! a proper set-up of those variables failes or is bypassed.
  ! ---------------------------------------------------------------
  first_pix = 1 ; last_pix = 1

  ! --------------------------------
  ! Name of the main output molecule
  ! --------------------------------
  molname = sao_molecule_names(pge_idx)

  ! -------------------------------------------------------------------
  ! Range of cross-track pixels to fit. This is based on the selection
  ! in the fitting control file and whether the granule being processed
  ! is in global or spatial zoom mode, or even a mixture thereof.
  !
  ! NOTE that we set OMI_XTRPIX_RANGE for all swath lines because the
  ! choice of swath lines to process may not contain the radiance
  ! reference/calibration line.
  ! -------------------------------------------------------------------
  CALL omi_read_binning_factor ( &
    TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
    ntimes_rad, omi_binfac, omi_is_szoom, errstat )
  if (errstat < 0) return

  CALL omi_set_xtrpix_range ( &
    ntimes_rad, nxtrack_rad, pixnum_lim(3:4),                         &
    omi_binfac, omi_xtrpix_range, &
    first_wc_pix, last_wc_pix, errstat )

  if (errstat < 0) return

  ! --------------------------------------------------------------------
  ! If the radiance reference is obtained from the same L1b file, we can
  ! simply copy the variables we have just read to the corresponding
  ! "rr" ones (in this case, the dimensions are the same). Otherwise we
  ! have to read them from the radiance reference granule.
  ! --------------------------------------------------------------------
  IF ( TRIM(ADJUSTL(l1b_radref_filename)) /= TRIM(ADJUSTL(l1b_rad_filename)) ) THEN
    CALL omi_read_binning_factor ( &
      TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
      ntimes_rr, omi_binfac_rr, omi_is_szoom_rr, &
      errstat )
    CALL omi_set_xtrpix_range ( &
      ntimes_rr, nxtrack_rad, pixnum_lim(3:4),                                 &
      omi_binfac_rr, omi_xtrpix_range_rr, &
      first_wc_pix, last_wc_pix, errstat )
    if (errstat < 0) return
  ELSE
    omi_binfac_rr      (0:ntimes_rad-1)     = omi_binfac      (0:ntimes_rad-1)
    omi_is_szoom_rr    (0:ntimes_rad-1)     = omi_is_szoom    (0:ntimes_rad-1)
    omi_xtrpix_range_rr(0:ntimes_rad-1,1:2) = omi_xtrpix_range(0:ntimes_rad-1,1:2)
  END IF

  call irradiance_data_init (rpt_rad, errstat);
  if (errstat < 0) return

  ! ---------------------------------------------------------------
  ! Solar wavelength calibration, done even when we use a composite
  ! solar spectrum to avoid un-initialized variables. However, no
  ! actual fitting is performed in the latter case.
  ! ---------------------------------------------------------------
  CALL xtrack_solar_calibration_loop ( first_wc_pix, last_wc_pix, errstat )
  if (errstat < 0) return

  ! ---------------------------------------------------------------
  ! No matter what, we need a swath line for radiance wavelength
  ! calibration. This may be a single line or it may be the average
  ! over a block of lines. However, if we are not using a radiance
  ! reference, then we are still doing a radiance calibration and
  ! need to make sure that we are using a radiance from the current
  ! granule.
  ! ---------------------------------------------------------------

  ! Should radiance and radiance-ref swathnames be equal?  The original
  ! code did not have this restriction.  --JED
  if (trim(omi_radiance_swathname) /= trim(rpt_rr%swathname)) then
    write (*,*) "swathnames are not equal: ", trim(omi_radiance_swathname), &
      " /= ", trim(rpt_rr%swathname)
    write (*,*) "modify omi_get_radiance_reference to use omi_radiance_swathname"
    stop
  endif
  CALL omi_get_radiance_reference (rpt_rr, &
                                   omi_xtrpix_range_rr, &
                                   radiance_wavcal_lnums, errstat)
  if (errstat < 0) return

  ! ---------------------------------------------------------------
  ! The Climatology is going to be read here and kept in memory. If
  ! this has a bad impact in the efficiency of the application then
  ! I will find a different way. We are doing this to be able to in
  ! itialize the output he5 with the correct number of levels for
  ! the Scattering weights and Gas_profiles output. We are going to
  ! use the number of levels in the climatology as the number of le
  ! vels of the reported scattering weights.
  ! ---------------------------------------------------------------
  CALL omi_read_climatology (pge_idx, errstat )

  ! ----------------------------------------
  ! Initialization of HE5 output data fields
  ! ----------------------------------------

  ! FIXME: error handling needs worked here
  errstat = HE5_Init_Swath ( l2_filename, pge_swath_name, ntimes_rad, nxtrack_rad, CmETA )
  if (errstat < 0) return

  CALL he5_initialize_datafields ( )
  errstat = HE5_Define_Fields ( pge_idx, pge_swath_name, ntimes_rad, nxtrack_rad, CmETA )
  if (errstat < 0) return

  ! -----------------------------------------------------------------------------------
  ! If we are NOT using a radiance reference, then we need to read the
  ! swath line for radiance wavelength calibration. In this case, the
  ! value of  RADIANCE_WAVCAL_LNUMS is that of the first line of the radiance reference
  ! selected
  ! -----------------------------------------------------------------------------------
  IF ( .NOT. yn_radiance_reference ) THEN
    ntimes_loop = 1 ! The number of scan lines to read
    iline = radiance_wavcal_lnums(1)
    ! Get NTIMES_LOOP radiance lines.
    ! Note: omi_read_radiance_lines sets the global omi_nwav_rad
    CALL omi_read_radiance_lines (&
      l1b_rad_filename, iline, nxtrack_rad, ntimes_loop, &
      nwavel_rad, errstat )
    if (errstat < 0) return
  END IF

  ! -----------------------------------------------------
  ! Across-track loop for radiance wavelength calibration
  ! -----------------------------------------------------
  omi_radcal_xflag = i2_missval
  CALL xtrack_radiance_wvl_calibration (                          &
    first_wc_pix, last_wc_pix, n_max_rspec, n_comm_wvl, errstat )
  if (errstat < 0) return

  ! --------------------------------------------------------------
  ! Terminate on not having any cross-track pixels left to process
  ! --------------------------------------------------------------
  IF ( ALL ( omi_cross_track_skippix ) ) THEN
    CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_NOPIXEL, &
      modulename, vb_lev_default, errstat )
    GO TO 400
  END IF

  if (.not.yn_radiance_reference) then
    call init_prefit_files (pge_idx, ntimes_rad, nxtrack_rad, errstat)
    if ( errstat >= pge_errstat_error ) return
  endif

  ! ---------------------------------------------------------------------
  ! If we are using a radiance reference AND want to remove the target
  ! gas from it (important for BrO, for example), we have to run through
  ! all spectra that go into the reference, compute the average column,
  ! and then remove that from the averaged radiance reference spectrum
  ! (owing to the fact that the average of hundreds of OMI spectra still
  !  doesn't produce a decently fitted column).
  ! ---------------------------------------------------------------------
  IF ( yn_radiance_reference .AND. yn_remove_target ) THEN

    iline = SUM ( radiance_reference_lnums(1:2) ) / 2
    IF ( iline < 0 .OR. iline > ntimes_rr ) iline = ntimes_rr / 2

    first_pix = omi_xtrpix_range_rr(iline,1)
    last_pix  = omi_xtrpix_range_rr(iline,2)

    CALL xtrack_radiance_reference_loop ( &
      yn_remove_target, & ! note: yn_remove_target=TRUE here
      nxtrack_rr, nwavel_rr, first_pix, last_pix, pge_idx, errstat )
    if (errstat < 0) return

    ! -------------------------------------------------------------
    ! Write the output from solar/earthshine wavelength calibration
    ! and radiance reference to file. The latter results will be
    ! overwritten in the call to XTRACK_RADIANCE_REFERENCE_LOOP
    ! below, hence we need to write them out here.
    ! -------------------------------------------------------------
    CALL he5_write_wavcal_output ( nxtrack_rad, first_pix, last_pix, errstat )

  END IF

  ! -----------------------------------------------------------------
  ! Before we go any further we need to read the L1b latitude values,
  ! since we base our screening of which swath lines to process on
  ! those values. Both common mode, if used, and the radiance fit
  ! uses the same arrays, so we read this only ones.
  !
  ! We could shave off some fractional minute from the run time by
  ! not reading the latitudes in cases where no radiance reference
  ! is used, i.e., where both radiance granule and radiance reference
  ! granule are the same. The down-side is an increase in virtual
  ! memory program uses, plus some more logic to find out whether to
  ! read the latitudes or not. For now we are going with a second
  ! read, particularly since the current algorithm settings would
  ! require it anyway.
  ! -----------------------------------------------------------------

  allocate (l1b_rad_latitudes (1:nxtrack_rad, 0:ntimes_rad-1), stat=locerrstat)
  if (locerrstat /= 0) then
    call err_message_error (modulename // ": allocate failed", errstat)
    return
  endif

  CALL read_latitude ( &
    TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
    0, ntimes_rad, l1b_rad_latitudes, errstat)

  if (errstat < 0) return

  ! -----------------------------------------------------------------
  ! Now we enter the on-line computation of the common mode spectrum.
  ! -----------------------------------------------------------------
  IF ( yn_common_iter ) THEN

    ! ----------------------------------------------------------
    ! Set the logical YN array that determines which swath lines
    ! will be used in the common mode
    ! ----------------------------------------------------------
    is_common_range = .FALSE.
    CALL find_swathline_range ( &
      TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
      ntimes_rad, nxtrack_rad, l1b_rad_latitudes,       &
      common_latrange(1:2), is_common_range, errstat             )

    ! -------------------------------------------------------------
    ! First and last swath line number will be overwritten. Hence
    ! we save them for further reference.
    ! -------------------------------------------------------------
    !first_line_save = first_line
    !last_line_save  = last_line

    ! ------------------------------------------------------------------
    ! There is a time saver catch built into the assignment below,
    ! which we probably want to rethink. If we are only doing a few
    ! lines but the common mode extends over a wide range of a
    ! latitudes, then we would be processing a lot of lines to derive
    ! the common mode. This is currently being excluded, but the better
    ! way may be to remember that we can control everything through
    ! the fitting control file.
    ! ------------------------------------------------------------------
    !IF ( MINVAL(common_latlines(1:2)) >= 0 .AND. &
    !     MAXVAL(common_latlines(1:2)) <= ntimes_rad ) THEN
    !   first_line = common_latlines(1) !MAX(first_line, common_latlines(1))
    !   last_line  = common_latlines(2) !MIN(last_line,  common_latlines(2))
    !END IF

    ! ------------------------------------------
    ! Interface to the loop over all swath lines
    ! ------------------------------------------
    write(*,*)' omi_fitting: calling swathline_loops (common mode)'
    CALL swathline_loops ( &
      pge_idx, rpt_rad, n_max_rspec, &
      is_common_range, &
      omi_xtrpix_range, &
      .FALSE., -1, &
      .TRUE., errstat )

    if (errstat < 0) return

    ! -----------------------------------------------------------------
    ! Reset first and last swath line number to non-common mode values.
    ! -----------------------------------------------------------------
    !first_line = first_line_save
    !last_line  = last_line_save

    ! ---------------------------------------------------
    ! Set the index value of the Common Mode spectrum and
    ! assign values to the fitting parameter arrays
    ! ---------------------------------------------------
    CALL finalize_common_mode (nxtrack_rad)

    ! -------------------------------------------
    ! Write the just computed common mode to file
    ! -------------------------------------------
    IF ( yn_diagnostic_run ) then
      CALL he5_write_common_mode ( nxtrack_rad, n_comm_wvl, errstat )
      if (errstat < 0) return
    endif

  END IF

  ! ----------------------------------------------------------
  ! Now into the proper fitting, with or without common mode.
  ! ----------------------------------------------------------

  ! -------------------------------------------------------------------
  ! Radiance Reference Fit: Only if we have not selected to remove the
  ! target gas from the radiance reference (in which case we have done
  ! this fit already), or N_COMM_ITER==0 (which means that we can refit
  ! with a common mode spectrum)
  ! -------------------------------------------------------------------
  IF ( radiance_wavcal_lnums(1) >= 0 ) THEN

    ! --------------------------------
    ! The number of scan lines to read
    ! --------------------------------
    ntimes_loop = 1
    iline = radiance_wavcal_lnums(1)

    ! -------------------------------------------------
    ! Only use prefit if not using a Radiance Reference
    ! -------------------------------------------------
    IF (.NOT. yn_radiance_reference ) THEN
      CALL read_prefit_columns ( pge_idx, nxtrack_rad, ntimes_loop, iline, errstat )
      if (errstat < 0) return
    END IF

    ! ------------------------------
    ! Get NTIMES_LOOP radiance lines
    ! ------------------------------
    CALL omi_read_radiance_lines ( &
      l1b_rad_filename, iline, nxtrack_rad, ntimes_loop, nwavel_rad, errstat )
    if (errstat < 0) return

    ! -----------------------------------------------
    ! Radiance Reference Fit (or WavCal Radiance Fit)
    ! -----------------------------------------------
    first_pix = omi_xtrpix_range(iline,1)
    last_pix  = omi_xtrpix_range(iline,2)

    ! -------------------------------------
    ! Initialize saved fitting variables
    ! -------------------------------------
    ! fitvar_rad_saved(1:n_max_fitpars ) = fitvar_rad_init(1:n_max_fitpars)
    ! Not needed --- xtrack_radiance_reference_loop does this.  --JED
    CALL xtrack_radiance_reference_loop ( &
      .FALSE., nxtrack_rr, nwavel_rr, first_pix, last_pix, pge_idx, errstat )

  END IF

  ! -------------------------------------------------------------------
  ! Output of fit results for solar and radiance wavelength calibration
  ! but ONLY if we haven't done it already (see above). The Radiance
  ! Reference values under YN_REMOVE_TARGET settings carry valuable
  ! information only BEFORE the target has been removed.
  ! -------------------------------------------------------------------
  IF ( .NOT. (yn_radiance_reference .AND. yn_remove_target) ) THEN
    CALL he5_write_wavcal_output ( nxtrack_rr, first_pix, last_pix, errstat )
    if (errstat < 0) return
  END IF

  ! ----------------------------------------------------------
  ! Set the logical YN array that determines which swath lines
  ! will be processed. Unless we have constrained either swath
  ! line numbers or the latitude range, this can be set to
  ! .TRUE. universically.
  ! ----------------------------------------------------------
  ! First, set the range of swath lines to process
  ! ----------------------------------------------
  first_line = 0  ;  last_line = ntimes_rad-1
  IF ( pixnum_lim(1) > 0 ) first_line = MIN(pixnum_lim(1), last_line)
  IF ( pixnum_lim(2) > 0 ) last_line  = MAX( MIN(pixnum_lim(2), last_line), first_line )

  do_radfit_range = .FALSE.
  IF ( first_line         > 0           .OR. &
    last_line          < ntimes_rad-1 .OR. &
    radfit_latrange(1) > -90.0_r4    .OR. &
    radfit_latrange(2) < +90.0_r4           ) THEN

    IF ( radfit_latrange(1) > -90.0_r4    .OR. &
      radfit_latrange(2) < +90.0_r4           ) THEN
      write(*,*)' omi_fitting:  setting do_radfit_range:  call find_swathline_range (radiances)'
      CALL find_swathline_range ( &
        TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
        ntimes_rad, nxtrack_rad, l1b_rad_latitudes,       &
        radfit_latrange(1:2), do_radfit_range, errstat             )
    ELSE
      write(*,*)' omi_fitting: setting do_radfit_range:  applying first_line/last_line mask'
      do_radfit_range = .TRUE.
      IF ( first_line > 0           ) do_radfit_range(0:first_line-1)          = .FALSE.
      IF ( last_line  < ntimes_rad-1 ) do_radfit_range(last_line+1:ntimes_rad-1) = .FALSE.
    END IF
  ELSE
    write(*,*)' omi_fitting: setting do_radfit_range:  all TRUE'
    do_radfit_range = .TRUE.
  END IF

  deallocate (l1b_rad_latitudes)

  ! ------------------------------------------
  ! Interface to the loop over all swath lines
  ! ------------------------------------------
  write(*,*)' omi_fitting: calling swathline_loops (radiances)'
  CALL swathline_loops ( &
    pge_idx, rpt_rad, n_max_rspec,     &
    do_radfit_range,                           &
    omi_xtrpix_range,                      &
    .FALSE., -1,                       &
    .FALSE., errstat)
  if (errstat < 0) return

  ! FIXME: Instead of jumping here, the following should be put
  ! into a separate function, and that function called here and
  ! where the goto happended.  --JED
 400 CONTINUE

  ! ----------------------------------------------
  ! This subroutine completes the following tasks:
  !    (1) Compute pixel geololcation corners
  !    (2) Compute AMFs
  !    (3) Apply cross-track destriping correction
  ! ---------------------------------------
  CALL omi_pge_postprocess ( &
    l1b_rad_filename, pge_idx, ntimes_rad, nxtrack_rad,                    &
    do_radfit_range, omi_xtrpix_range, &
    omi_is_szoom, n_max_rspec, errstat                 )

  ! ---------------------
  ! Write some attributes
  ! ---------------------
  errstat = he5_write_global_attributes( )
  if (errstat < 0) then
    call err_message_error (modulename//f_sep// &
                            "he5_write_global_attributes failed", &
                            errstat)
    return
  endif

  errstat = he5_write_swath_attributes ( pge_idx )
  if (errstat < 0) then
    call err_message_error (modulename//f_sep// &
                            "he5_write_swath_attributes failed", &
                            errstat)
    return
  endif

  errstat = he5_set_field_attributes   ( pge_idx )
  if (errstat < 0) then
    call err_message_error (modulename//f_sep// &
                            "he5_set_field_attributes failed", &
                            errstat)
    return
  endif
  ! -----------------
  ! Close output file
  ! -----------------
  errstat = he5_close_output_file ( pge_idx)
  if (errstat < 0) then
    call err_message_error (modulename//f_sep// &
                            "he5_close_output_file failed", &
                            errstat)
    return
  endif

  ! ----------------------------------------
  ! Get Metadata from MCF, PCF, and L1B file
  ! ----------------------------------------
  CALL check_metadata_consistency ( errstat )
  if (errstat < 0) then
      call err_message_error (modulename//f_sep// &
                            "check_metadata_consistency failed", &
                            errstat)
      return
  endif

  CALL set_l2_metadata ( pge_idx, errstat )
  if (errstat < 0) then
      call err_message_error (modulename//f_sep// &
                            "set_l2_metadata failed", &
                            errstat)
      return
  endif

END SUBROUTINE omi_fitting

END MODULE
