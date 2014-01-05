MODULE omi_pge_fitting_process
CONTAINS
SUBROUTINE omi_pge_fitting ( pge_idx, n_max_rspec, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,      ONLY: pge_errstat_ok, pge_errstat_error, pge_errstat_fatal
  USE OMSAO_he5_module,          ONLY: NrofScanLines, NrofCrossTrackPixels
  USE OMSAO_variables_module,    ONLY: l1b_rad_filename, Radiance_Paras_Type, &
    yn_radiance_reference, l1b_radref_filename
  USE OMSAO_omidata_module,      ONLY: omi_radiance_swathname
  USE omi_pge_fitting_aux, ONLY: omi_set_fitting_parameters
  USE omi_read_l1b_data, ONLY: omi_read_radiance_paras
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
  CALL omi_read_radiance_paras (l1b_rad_filename, rpt_rad, errstat)
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
    CALL omi_read_radiance_paras (l1b_radref_filename, rpt_rr, errstat)
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
                        n_max_rspec, pge_error_status)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: i2_missval, MAX_STR_LEN
  USE OMSAO_indices_module,    ONLY: sao_molecule_names, pge_hcho_idx, pge_gly_idx
  USE OMSAO_variables_module,  ONLY: &
    l1b_rad_filename, &
    l2_filename, pixnum_lim,    &
    radfit_latrange,                &
    yn_solar_comp, yn_diagnostic_run,              &
    yn_common_iter, common_latrange,    &
    radiance_wavcal_lnums, yn_solmonthave, Radiance_Paras_Type, &
    radiance_reference_lnums, l1b_radref_filename, yn_radiance_reference, &
    yn_remove_target !, fitvar_rad_init, fitvar_rad_saved
  USE OMSAO_he5_module,       ONLY:  pge_swath_name
  USE OMSAO_solar_wavcal_module, ONLY: xtrack_solar_calibration_loop
  USE OMSAO_radiance_ref_module, ONLY: omi_get_radiance_reference, &
    xtrack_radiance_reference_loop
  USE OMSAO_prefitcol_module, ONLY: read_prefit_columns, init_prefit_files, &
    bro_prefit_col, bro_prefit_dcol, lqh2o_prefit_col, lqh2o_prefit_dcol, &
    o3_prefit_col, o3_prefit_dcol, &
    yn_o3_prefit, yn_bro_prefit, yn_lqh2o_prefit
  USE OMSAO_errstat_module
  USE OMSAO_solmonthave_module, ONLY: omi_read_monthly_average_irradiance
  USE OMSAO_wfamf_module, ONLY: omi_read_climatology, CmETA
  USE he5_output_tools, ONLY: he5_init_swath, he5_define_fields, &
    he5_close_output_file, he5_set_field_attributes, &
    he5_write_global_attributes, he5_write_swath_attributes, &
    he5_write_wavcal_output, he5_write_common_mode !, he5_open_readwrite
  USE omi_read_l1b_data, ONLY: omi_read_binning_factor, omi_read_irradiance_data, &
    omi_read_radiance_lines, omi_read_radiance_lines
  USE omi_pge_fitting_aux, ONLY: omi_set_xtrpix_range, omi_create_solcomp_irradiance, &
    read_latitude, find_swathline_range, finalize_common_mode
  USE fitting_loops, ONLY: xtrack_radiance_wvl_calibration
  USE metadata_tools, ONLY: check_metadata_consistency, set_l2_metadata
  USE omi_pge_postprocessing, ONLY: omi_pge_postprocess
  USE omi_pge_swathline_loop, ONLY: omi_pge_swathline_loops
  use datafields, only: he5_initialize_datafields
  USE OMSAO_omidata_module, ONLY: n_comm_wvl, ntimes_loop, &
    omi_cross_track_skippix, omi_radcal_itnum, omi_radcal_xflag, &
    omi_radiance_swathname, omi_solcal_itnum, omi_solcal_xflag

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec
  TYPE(Radiance_Paras_Type), INTENT(IN) :: rpt_rad, rpt_rr

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! -------------------------
  ! Local variables (for now)
  ! -------------------------
  INTEGER   (KIND=i4) ::                                                 &
    iline, first_line, last_line, errstat, first_wc_pix, last_wc_pix, &
    first_pix, last_pix
  INTEGER (kind=i4) :: ntimes_rad, nxtrack_rad, nwavel_rad
  INTEGER (kind=i4) :: ntimes_rr, nxtrack_rr, nwavel_rr

  ! ----------------------------------------------------------------------
  ! Swath dimensions and variables that aren't passed from calling routine
  ! ----------------------------------------------------------------------
  INTEGER   (KIND=i4)      :: nTimesIrr, nXtrackIrr
  CHARACTER (LEN=MAX_STR_LEN) :: molname

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! and Common Mode spectrum
  ! ----------------------------------------------------------
  INTEGER (KIND=i1), DIMENSION (0:rpt_rad%ntimes-1)   :: omi_binfac
  INTEGER (KIND=i4), DIMENSION (0:rpt_rad%ntimes-1,2) :: omi_xtrpix_range
  LOGICAL,           DIMENSION (0:rpt_rad%ntimes-1)   :: &
    omi_yn_szoom, yn_common_range, yn_radfit_range

  INTEGER (KIND=i1), DIMENSION (0:rpt_rr%ntimes-1)   :: omi_binfac_rr
  INTEGER (KIND=i4), DIMENSION (0:rpt_rr%ntimes-1,2) :: omi_xtrpix_range_rr
  LOGICAL,           DIMENSION (0:rpt_rr%ntimes-1)   :: omi_yn_szoom_rr

  ! ----------------------------------------------------------
  ! OMI L1b latitudes
  ! ----------------------------------------------------------
  REAL (KIND=r4), DIMENSION (1:rpt_rad%nxtrack, 0:rpt_rad%ntimes-1) :: l1b_latitudes
  REAL (KIND=r4), DIMENSION(:,:), ALLOCATABLE :: latitudes_rr

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

  errstat = pge_errstat_ok

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
    ntimes_rad, omi_binfac, omi_yn_szoom, errstat )
  CALL omi_set_xtrpix_range ( &
    ntimes_rad, nxtrack_rad, pixnum_lim(3:4),                         &
    omi_binfac, omi_xtrpix_range, &
    first_wc_pix, last_wc_pix, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! --------------------------------------------------------------------
  ! If the radiance reference is obtained from the same L1b file, we can
  ! simply copy the variables we have just read to the corresponding
  ! "rr" ones (in this case, the dimensions are the same). Otherwise we
  ! have to read them from the radiance reference granule.
  ! --------------------------------------------------------------------
  IF ( TRIM(ADJUSTL(l1b_radref_filename)) /= TRIM(ADJUSTL(l1b_rad_filename)) ) THEN
    CALL omi_read_binning_factor ( &
      TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
      ntimes_rr, omi_binfac_rr, omi_yn_szoom_rr, &
      errstat )
    CALL omi_set_xtrpix_range ( &
      ntimes_rr, nxtrack_rad, pixnum_lim(3:4),                                 &
      omi_binfac_rr, omi_xtrpix_range_rr, &
      first_wc_pix, last_wc_pix, errstat )
    pge_error_status = MAX ( pge_error_status, errstat )
    IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ELSE
    omi_binfac_rr      (0:ntimes_rad-1)     = omi_binfac      (0:ntimes_rad-1)
    omi_yn_szoom_rr    (0:ntimes_rad-1)     = omi_yn_szoom    (0:ntimes_rad-1)
    omi_xtrpix_range_rr(0:ntimes_rad-1,1:2) = omi_xtrpix_range(0:ntimes_rad-1,1:2)
  END IF

  ! --------------------------------------------------------------------
  ! Solar Irradiance Processing: If we don't do a solar composite, we can
  ! use a solar monthly average, if not we have to read the irradiance
  ! data.
  ! Otherwise we need to compute them from the solar composite
  ! parameterization on a equidistant grid.
  ! -------------------------------------------------------------------
  omi_solcal_itnum = i2_missval ; omi_solcal_xflag = i2_missval
  ! --------------------------------------------------------------------------
  ! Check than only one or non of yn_solar_comp are yn_solmonthva are set True
  ! --------------------------------------------------------------------------
  IF ( yn_solar_comp .AND. yn_solmonthave ) THEN
    CALL error_check ( 1, 0, pge_errstat_fatal, OMSAO_F_SOLCOM_VS_SOLAVE, &
      modulename, vb_lev_gt1mb, errstat )
    GO TO 666
  END IF

  IF ( yn_solar_comp ) THEN
    ! -----------------------------------
    ! Compute composite solar irradiances
    ! -----------------------------------
    CALL omi_create_solcomp_irradiance ( nxtrack_rad )
  ELSE IF (yn_solmonthave) THEN
    ! -----------------------------------
    ! Read solar monthly mean irradiance
    ! -----------------------------------
    CALL omi_read_monthly_average_irradiance (nTimesIrr, nXtrackIrr, errstat)
    pge_error_status = MAX ( pge_error_status, errstat )
    IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ELSE
    ! --------------------
    ! Read OMI irradiances
    ! --------------------
    CALL omi_read_irradiance_data ( nTimesIrr, nXtrackIrr, errstat )
    pge_error_status = MAX ( pge_error_status, errstat )
    IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  END IF

  ! ---------------------------------------------------------------
  ! Solar wavelength calibration, done even when we use a composite
  ! solar spectrum to avoid un-initialized variables. However, no
  ! actual fitting is performed in the latter case.
  ! ---------------------------------------------------------------
  CALL xtrack_solar_calibration_loop ( first_wc_pix, last_wc_pix, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! ---------------------------------------------------------------
  ! No matter what, we need a swath line for radiance wavelength
  ! calibration. This may be a single line or it may be the average
  ! over a block of lines. However, if we are not using a radiance
  ! reference, then we are still doing a radiance calibration and
  ! need to make sure that we are using a radiance from the current
  ! granule.
  ! ---------------------------------------------------------------
  ALLOCATE (latitudes_rr(1:nxtrack_rr,0:ntimes_rr-1), STAT=errstat)
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  CALL read_latitude (l1b_radref_filename, omi_radiance_swathname, &
                      ntimes_rr, nxtrack_rr, latitudes_rr)
  CALL omi_get_radiance_reference (l1b_radref_filename, &
                                   ntimes_rr, nxtrack_rr, nwavel_rr, &
                                   omi_xtrpix_range_rr, latitudes_rr, &
                                   radiance_wavcal_lnums, errstat)
  DEALLOCATE(latitudes_rr)
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

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
  errstat = HE5_Init_Swath ( l2_filename, pge_swath_name, ntimes_rad, nxtrack_rad, CmETA )
  CALL he5_initialize_datafields ( )
  errstat = HE5_Define_Fields ( pge_idx, pge_swath_name, ntimes_rad, nxtrack_rad, CmETA )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

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
    pge_error_status = MAX ( pge_error_status, errstat )
    IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  END IF

  ! -----------------------------------------------------
  ! Across-track loop for radiance wavelength calibration
  ! -----------------------------------------------------
  omi_radcal_itnum = i2_missval ; omi_radcal_xflag = i2_missval
  CALL xtrack_radiance_wvl_calibration (                          &
    first_wc_pix, last_wc_pix, n_max_rspec, n_comm_wvl, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! --------------------------------------------------------------
  ! Terminate on not having any cross-track pixels left to process
  ! --------------------------------------------------------------
  IF ( ALL ( omi_cross_track_skippix ) ) THEN
    CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_NOPIXEL, &
      modulename, vb_lev_default, errstat )
    GO TO 666
  END IF

  ! ----------------------------------
  ! CCM - Add lqH2O prefit
  ! ----------------------------------
  SELECT CASE(pge_idx )

  CASE( pge_hcho_idx )

    ! -------------------------------------------------------------------
    ! First access to pre-fitted O3 and BrO columns. At this time we have
    ! already set up some of the fitting arrays, so any error would lead
    ! to some major headaches to undo things. Hence we return if access
    ! fails.
    ! -------------------------------------------------------------------
    IF ( ( .NOT. yn_radiance_reference ) .AND. &
      ( pge_idx == pge_hcho_idx ) .AND. &
      ANY((/yn_o3_prefit(1),yn_bro_prefit(1)/)) ) THEN
      CALL init_prefit_files ( pge_idx, ntimes_rad, nxtrack_rad, errstat )
      IF ( errstat >= pge_errstat_error ) RETURN
    END IF

  CASE( pge_gly_idx )

    ! -------------------------------------------------------------------
    ! First access to pre-fitted Liquid Water columns. At this time we
    ! have already set up some of the fitting arrays, so any error would
    ! lead to some major headaches to undo things. Hence we return if
    ! access fails.
    ! -------------------------------------------------------------------

    IF ( (.NOT. yn_radiance_reference) .AND. yn_lqh2o_prefit(1) ) THEN
      CALL init_prefit_files ( pge_idx, ntimes_rad, nxtrack_rad, errstat )
      IF ( errstat >= pge_errstat_error ) RETURN
    END IF

  CASE DEFAULT
    ! Nothing here yet
  END SELECT

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

    CALL xtrack_radiance_reference_loop (                                   &
      yn_radiance_reference, yn_remove_target,                           &
      nxtrack_rr, nwavel_rr, first_pix, last_pix, pge_idx, pge_error_status )

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
  CALL read_latitude ( &
    TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
    ntimes_rad, nxtrack_rad, l1b_latitudes)

  ! -----------------------------------------------------------------
  ! Now we enter the on-line computation of the common mode spectrum.
  ! -----------------------------------------------------------------
  IF ( yn_common_iter ) THEN

    ! ----------------------------------------------------------
    ! Set the logical YN array that determines which swath lines
    ! will be used in the common mode
    ! ----------------------------------------------------------
    yn_common_range = .FALSE.
    CALL find_swathline_range ( &
      TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
      ntimes_rad, nxtrack_rad, l1b_latitudes,       &
      common_latrange(1:2), yn_common_range, errstat             )

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
    CALL omi_pge_swathline_loops ( &
      pge_idx, rpt_rad, n_max_rspec, &
      yn_common_range, &
      omi_xtrpix_range, &
      yn_radiance_reference, .FALSE., -1, &
      yn_common_iter, pge_error_status )

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
    IF ( yn_diagnostic_run ) &
      CALL he5_write_common_mode ( nxtrack_rad, n_comm_wvl, pge_error_status )

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
    ! CCM modify to include lqH2O
    ! -------------------------------------------------
    IF ( ( .NOT. yn_radiance_reference )                             .AND. &
      ( (pge_idx == pge_hcho_idx) .OR. (pge_idx == pge_gly_idx) ) .AND. &
      ANY((/yn_o3_prefit(1),yn_bro_prefit(1),yn_lqh2o_prefit(1)/)) ) THEN
      CALL read_prefit_columns ( pge_idx, nxtrack_rad, ntimes_loop, iline, errstat )
      pge_error_status = MAX ( pge_error_status, errstat )
      IF ( errstat >= pge_errstat_error ) GO TO 666
    ELSE
      o3_prefit_col    = 0.0_r8 ; o3_prefit_dcol    = 0.0_r8
      bro_prefit_col   = 0.0_r8 ; bro_prefit_dcol   = 0.0_r8
      lqh2o_prefit_col = 0.0_r8 ; lqh2o_prefit_dcol = 0.0_r8
    END IF

    ! ------------------------------
    ! Get NTIMES_LOOP radiance lines
    ! ------------------------------
    CALL omi_read_radiance_lines ( &
      l1b_rad_filename, iline, nxtrack_rad, ntimes_loop, nwavel_rad, errstat )
    pge_error_status = MAX ( pge_error_status, errstat )
    IF ( pge_error_status >= pge_errstat_error )  GO TO 666

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
      yn_radiance_reference, .FALSE.,  &
      nxtrack_rr, nwavel_rr, first_pix, last_pix, pge_idx, errstat )

  END IF

  ! -------------------------------------------------------------------
  ! Output of fit results for solar and radiance wavelength calibration
  ! but ONLY if we haven't done it already (see above). The Radiance
  ! Reference values under YN_REMOVE_TARGET settings carry valuable
  ! information only BEFORE the target has been removed.
  ! -------------------------------------------------------------------
  IF ( .NOT. (yn_radiance_reference .AND. yn_remove_target) ) THEN
    CALL he5_write_wavcal_output ( nxtrack_rr, first_pix, last_pix, errstat )
    pge_error_status = MAX ( pge_error_status, errstat )
    IF ( pge_error_status >= pge_errstat_error )  GO TO 666
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

  yn_radfit_range = .FALSE.
  IF ( first_line         > 0           .OR. &
    last_line          < ntimes_rad-1 .OR. &
    radfit_latrange(1) > -90.0_r4    .OR. &
    radfit_latrange(2) < +90.0_r4           ) THEN

    IF ( radfit_latrange(1) > -90.0_r4    .OR. &
      radfit_latrange(2) < +90.0_r4           ) THEN
      CALL find_swathline_range ( &
        TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
        ntimes_rad, nxtrack_rad, l1b_latitudes,       &
        radfit_latrange(1:2), yn_radfit_range, errstat             )
    ELSE
      yn_radfit_range = .TRUE.
      IF ( first_line > 0           ) yn_radfit_range(0:first_line-1)          = .FALSE.
      IF ( last_line  < ntimes_rad-1 ) yn_radfit_range(last_line+1:ntimes_rad-1) = .FALSE.
    END IF
  ELSE
    yn_radfit_range = .TRUE.
  END IF

  ! ------------------------------------------
  ! Interface to the loop over all swath lines
  ! ------------------------------------------
  CALL omi_pge_swathline_loops ( &
    pge_idx, rpt_rad, n_max_rspec,     &
    yn_radfit_range,                           &
    omi_xtrpix_range,                      &
    yn_radiance_reference, .FALSE., -1,                       &
    .FALSE., pge_error_status )

  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
 666 CONTINUE
  IF ( pge_error_status >= pge_errstat_fatal ) RETURN

  ! ----------------------------------------------
  ! This subroutine completes the following tasks:
  !    (1) Compute pixel geololcation corners
  !    (2) Compute AMFs
  !    (3) Apply cross-track destriping correction
  ! ---------------------------------------
  CALL omi_pge_postprocess ( &
    l1b_rad_filename, pge_idx, ntimes_rad, nxtrack_rad,                    &
    yn_radfit_range, omi_xtrpix_range, &
    omi_yn_szoom, n_max_rspec, errstat                 )

  ! ---------------------
  ! Write some attributes
  ! ---------------------
  errstat = pge_errstat_ok
  errstat = he5_write_global_attributes( )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    modulename//f_sep//"HE5_WRITE_GLOBAL_ATTRIBUTES", vb_lev_default, pge_error_status )

  errstat = pge_errstat_ok
  errstat = he5_write_swath_attributes ( pge_idx )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    modulename//f_sep//"HE5_WRITE_SWATH_ATTRIBUTES", vb_lev_default, pge_error_status )

  errstat = pge_errstat_ok
  errstat = he5_set_field_attributes   ( pge_idx )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    modulename//f_sep//"HE5_SET_FIELD_ATTRIBUTES", vb_lev_default, pge_error_status )
  ! -----------------
  ! Close output file
  ! -----------------
  errstat = pge_errstat_ok
  errstat = he5_close_output_file ( pge_idx)
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    modulename//f_sep//"HE5_CLOSE_OUTPUT_FILE", vb_lev_default, pge_error_status )

  ! ----------------------------------------
  ! Get Metadata from MCF, PCF, and L1B file
  ! ----------------------------------------
  errstat = pge_errstat_ok
  CALL check_metadata_consistency ( errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    modulename//f_sep//"CHECK_METADATA_CONSISTENCY.", vb_lev_default, pge_error_status )
  errstat = pge_errstat_ok
  CALL set_l2_metadata ( pge_idx, errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    modulename//f_sep//"SET_L2_METADATA.", vb_lev_default, pge_error_status )

  RETURN
END SUBROUTINE omi_fitting

END MODULE
