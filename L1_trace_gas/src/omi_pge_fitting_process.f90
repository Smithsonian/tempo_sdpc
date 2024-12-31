MODULE omi_pge_fitting_process

  use tell_module
  use l1bread_utils
  private
  public omi_pge_fitting

CONTAINS

SUBROUTINE omi_pge_fitting ( pge_idx, n_max_rspec, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module,    ONLY: l1b_rad_filename, Radiance_Paras_Type, &
    l1b_radref_filename, l1b_channel, solcal_cache_mode
  use ctrlvars, only: yn_radiance_reference, yn_gems
  USE OMSAO_omidata_module,      ONLY: omi_radiance_swathname, EarthSunDistance
  USE omi_pge_fitting_aux, ONLY: omi_set_fitting_parameters
  USE omi_read_l1b_data, ONLY: read_earth_sun_distance 
  USE OMSAO_solcomp_module, ONLY: soco_pars_deallocate
  use m_read_gems, only: gems_read_l1_rad_info, gems_read_latitude, &
       gems_read_earth_sun_distance
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE(Radiance_Paras_Type) :: rpt_rad, rpt_rr

  if (errstat /= 0) return

  ! -------------------------------------------------------------------------------------
  ! Set the swath name of various ESDTs
  ! -------------------------------------------------------------------------------------
  CALL omi_set_fitting_parameters ( pge_idx, errstat )
  if (errstat /= 0) goto 666

  if (solcal_cache_mode /= 'save') then
    ! -----------------------------------------------------------------------------------
    ! Get dimensions the L1B radiance granule
    ! -----------------------------------------------------------------------------------
    if (.not. yn_gems) then
      call read_l1_radiance_info(l1b_rad_filename, l1b_channel, rpt_rad, errstat)
    else ! GEMS
      call gems_read_l1_rad_info (l1b_rad_filename, rpt_rad, errstat)
    endif
    if (errstat /= 0) goto 666

    if (.not. yn_gems) then
      call read_earth_sun_distance (l1b_rad_filename, EarthSunDistance, errstat)
    else
      call gems_read_earth_sun_distance(l1b_rad_filename, EarthSunDistance, &
          errstat)
    endif
    if (errstat /= 0) goto 666
    omi_radiance_swathname = rpt_rad%swathname

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
      if (.not. yn_gems) then !TEMPO
        call read_l1_radiance_info (l1b_radref_filename, l1b_channel, rpt_rr, &
            errstat)
      else !GEMS
        call gems_read_l1_rad_info (l1b_radref_filename, rpt_rr, errstat)
      endif
      if (errstat /= 0) goto 666
    ENDIF

    ! ----------------------------------------------------------------
    ! Number of cross-track positions must be the same; fold otherwise
    ! ----------------------------------------------------------------
    IF ( rpt_rad%nxtrack /= rpt_rr%nxtrack ) THEN
      call tell_error (tell_runtime_error, &
                      "Cross-track position different between L1b radiance granules", &
                      errstat)
      GO TO 666
    END IF
  end if

  CALL omi_fitting (pge_idx, rpt_rad, rpt_rr, n_max_rspec, errstat)
  if (errstat /= 0) goto 666

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

  RETURN
END SUBROUTINE omi_pge_fitting

SUBROUTINE omi_fitting (pge_idx, rpt_rad, rpt_rr, &
     n_max_rspec, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval, MAX_STR_LEN, &
       nwavel_max, nxtrack_max
  USE OMSAO_indices_module,    ONLY: sao_molecule_names, max_rs_idx
  USE OMSAO_variables_module,  ONLY: &
       l1b_rad_filename, l1b_irrad_filename, ecs_version_id, &
       l2_filename, pixnum_lim, n_fitvar_rad,   &
       radfit_latrange,                &
       common_latrange,    &
       Radiance_Paras_Type, &
       radiance_reference_lnums, l1b_radref_filename, common_mode_spec, &
       solcal_cache_mode, solcal_filename
  use ctrlvars, only: yn_radiance_reference, yn_common_iter, &
       yn_diagnostic_run, yn_remove_target, yn_disable_omi_features, &
       yn_do_he5_output, yn_wrt_odl, yn_gems, &
       yn_do_solar_cal, yn_write_solar_cal, yn_read_solar_cal, yn_exit_post_solar_cal
  USE OMSAO_he5_module,       ONLY:  pge_swath_name, n_lun_inp, lun_input
  USE OMSAO_solar_wavcal_module, ONLY: xtrack_solar_calibration_loop
  USE OMSAO_radiance_ref_module, ONLY: omi_get_radiance_reference, &
       xtrack_radiance_reference_loop
  USE OMSAO_prefitcol_module, ONLY: read_prefit_columns, init_prefit_files
  USE OMSAO_wfamf_module, ONLY: read_profiles_dimensions, CmETA, amf_wvl
  use output_tools, only : create_output_file, close_output_file, &
       write_fitting_statistics, write_common_mode, write_wavcal_output, &
       write_solar_wavecal_diagnostics, write_radiance_wavecal_diagnostics, &
       label_output_file
  USE he5_output_tools, ONLY: he5_init_swath, he5_define_fields, &
       he5_close_output_file, he5_set_field_attributes, &
       he5_write_global_attributes, he5_write_swath_attributes, &
       he5_write_wavcal_output, he5_write_common_mode
  USE omi_read_l1b_data, ONLY: omi_read_binning_factor, &
       omi_read_radiance_lines, omi_read_radiance_lines
  USE omi_pge_fitting_aux, ONLY: omi_set_xtrpix_range, &
       read_latitude, find_swathline_range, fitting_statistics_type
  use commonmode, only: finalize_common_mode
  USE fitting_loops, ONLY: xtrack_radiance_wvl_calibration
  USE metadata_tools, ONLY: check_metadata_consistency, set_l2_metadata
  USE omi_pge_postprocessing, ONLY: omi_pge_postprocess
  USE swathline_loop, ONLY: swathline_loops
  use datafields, only: he5_initialize_datafields
  USE OMSAO_omidata_module, ONLY: n_comm_wvl, ntimes_loop, correlation_names, &
       omi_cross_track_skippix, omi_radcal_xflag, &
       omi_radiance_swathname, result_vars
  USE irradiance_data, only: irradiance_data_init, Irr_Data
  use m_write_odl_metadata
  use OMSAO_casestring_module, only : upper_case
  use slitfunction_tempo, only : solarcal_write_file, solarcal_read_file
  use m_read_gems, only: gems_read_latitude

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
  INTEGER (kind=i4) :: ntimes_rr, nxtrack_rr, nwavel_rr, extension_dot
  INTEGER (KIND=i4), DIMENSION (2) :: radiance_wavcal_lnums
  real (kind=r8), dimension(:,:), allocatable :: save_solcal_wvl, &
       save_solcal_spec, save_solcal_resid, save_radcal_wvl, save_radcal_resid

  ! ----------------------------------------------------------------------
  ! Swath dimensions and variables that aren't passed from calling routine
  ! ----------------------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: molname
  character (len=1024) :: l2_filename_netcdf

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
  type (fitting_statistics_type) :: fit_stats
  character (len=256) :: logmsg


  if (errstat /= 0) return

  ! ---------------------------------------------------------------
  ! Some initializations that will save us headaches in cases where
  ! a proper set-up of those variables failes or is bypassed.
  ! ---------------------------------------------------------------
  first_pix = 1 ; last_pix = 1

  ! --------------------------------
  ! Name of the main output molecule
  ! --------------------------------
  molname = sao_molecule_names(pge_idx)

  if (solcal_cache_mode /= 'save') then
  
    ntimes_rad = rpt_rad%ntimes
    nxtrack_rad = rpt_rad%nxtrack
    nwavel_rad = rpt_rad%nwavel_ccd

    ntimes_rr = rpt_rr%ntimes
    nxtrack_rr = rpt_rr%nxtrack
    nwavel_rr = rpt_rr%nwavel_ccd

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
    if (errstat /= 0) return

    CALL omi_set_xtrpix_range ( &
        ntimes_rad, nxtrack_rad, pixnum_lim(3:4),                         &
        omi_binfac, omi_xtrpix_range, &
        first_wc_pix, last_wc_pix, errstat )

    if (errstat /= 0) return

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
      if (errstat /= 0) return
    ELSE
      omi_binfac_rr      (0:ntimes_rad-1)     = omi_binfac      (0:ntimes_rad-1)
      omi_is_szoom_rr    (0:ntimes_rad-1)     = omi_is_szoom    (0:ntimes_rad-1)
      omi_xtrpix_range_rr(0:ntimes_rad-1,1:2) = omi_xtrpix_range(0:ntimes_rad-1,1:2)
    END IF
  end if

  call irradiance_data_init (errstat);
  if (errstat /= 0) return

  ! ---------------------------------------------------------------
  ! Solar wavelength calibration, done even when we use a composite
  ! solar spectrum to avoid un-initialized variables. However, no
  ! actual fitting is performed in the latter case.
  ! ---------------------------------------------------------------

  ! Set logicals to decide how to proceed with solar calibration
  ! 'none': do solar calibration then proceed with fitting
  ! 'save': save solar calibration results then exit retrieval
  ! 'read': skip solar calibration and read results from file
  SELECT CASE (solcal_cache_mode)
     CASE ('none')
        yn_do_solar_cal = .TRUE.
        yn_exit_post_solar_cal = .FALSE.
        yn_write_solar_cal = .FALSE.
        yn_read_solar_cal = .FALSE.
     CASE ('save')
        ! ------------------------------------------
        ! Find the range of XT pixels to process
        ! ------------------------------------------
        first_wc_pix = MAX (                                        1, pixnum_lim(3) )
        last_wc_pix  = MAX ( MIN ( Irr_Data%nxtrack,  pixnum_lim(4) ), first_wc_pix     )
        yn_do_solar_cal = .TRUE.
        yn_exit_post_solar_cal = .TRUE.
        yn_write_solar_cal = .TRUE.
        yn_read_solar_cal = .FALSE.
     CASE ('read')
        yn_do_solar_cal = .FALSE.
        yn_exit_post_solar_cal = .FALSE. 
        yn_write_solar_cal = .FALSE.
        yn_read_solar_cal = .TRUE.
     CASE DEFAULT
        CALL tell_log(1, 'Solar calibration cache mode = ' // solcal_cache_mode)
        CALL tell_error (tell_runtime_error, &
           "omi_pge_fitting_process: solcal_cache_mode input is invalid",  errstat)
        RETURN
     END SELECT

  CALL tell_log(1, 'Solar calibration cache mode = ' // solcal_cache_mode)

  ALLOCATE (save_solcal_wvl(nwavel_max, nxtrack_max), &
            save_solcal_spec(nwavel_max, nxtrack_max), &
            save_solcal_resid(nwavel_max, nxtrack_max))
  save_solcal_wvl(:,:) = r8_missval
  save_solcal_spec(:,:) = r8_missval
  save_solcal_resid(:,:) = r8_missval

  
  if (yn_do_solar_cal) then
    CALL xtrack_solar_calibration_loop (first_wc_pix, last_wc_pix, &
       save_solcal_wvl, save_solcal_spec, save_solcal_resid, errstat)
    if (errstat /= 0) return

    if (yn_write_solar_cal) THEN
      ! Save the solar calibration results
      CALL solarcal_write_file(l1b_irrad_filename, molname, solcal_filename, &
        save_solcal_wvl, save_solcal_spec, save_solcal_resid, errstat)
      if (errstat /= 0) return 
    endif
  endif

  if (yn_read_solar_cal) then
      ! real previously computed solar calibration data
      CALL solarcal_read_file(molname, solcal_filename, &
           save_solcal_wvl, save_solcal_spec, save_solcal_resid, errstat)
      if (errstat /= 0) return      
  endif

  if (yn_exit_post_solar_cal) then 
     deallocate (save_solcal_wvl, save_solcal_spec, save_solcal_resid, stat = errstat)
     if (errstat /= 0) then
        call tell_error(tell_malloc_error, "omi_fitting: deallocate failed", &
                        errstat)
     endif
     return
  endif

  ! ---------------------------------------------------------------
  ! No matter what, we need a swath line for radiance wavelength
  ! calibration. This may be a single line or it may be the average
  ! over a block of lines. However, if we are not using a radiance
  ! reference, then we are still doing a radiance calibration and
  ! need to make sure that we are using a radiance from the current
  ! granule.
  ! JCH:  The radiance reference also yields a value for n_comm_wvl
  !       that accounts for the fit window definition.
  ! ---------------------------------------------------------------
  ! Should radiance and radiance-ref swathnames be equal?  The original
  ! code did not have this restriction.  --JED
  if (.not.yn_gems .and. trim(omi_radiance_swathname) /= trim(rpt_rr%swathname)) then
    write (*,*) "swathnames are not equal: ", trim(omi_radiance_swathname), &
         " /= ", trim(rpt_rr%swathname)
    write (*,*) "modify omi_get_radiance_reference to use omi_radiance_swathname"
    stop
  endif
  call tell_log (1, 'omi_fitting: calling omi_get_radiance_reference')
  CALL omi_get_radiance_reference (rpt_rr, &
       omi_xtrpix_range_rr, &
       radiance_wavcal_lnums, errstat)
  if (errstat /= 0) return

  ! -------------------------------------------------
  ! Obtain vertical dimension of apriori gas profiles
  ! -------------------------------------------------
  call tell_log (1, 'omi_fitting: calling read_profiles_dimensions')
  CALL read_profiles_dimensions ( errstat )


  ! FIXME: for now, we'll define the netcdf output file name by
  ! changing the file extension.
  extension_dot = index(l2_filename,".",.true.) ! find rightmost '.'
  if (extension_dot > 1) then
    l2_filename_netcdf = l2_filename(1:extension_dot-1)//".nc"
  else
    l2_filename_netcdf = "l2_output.nc"
  endif

  write(logmsg,'(a,i4,a)')'omi_fitting: n_comm_wvl=',n_comm_wvl, &
       ', calling create_output_file: '//trim(l2_filename_netcdf)
  call tell_log (1, logmsg)
  call create_output_file (l2_filename_netcdf, pge_idx, amf_wvl, &
                           ntimes_rad, nxtrack_rad, CmETA, &
                           n_comm_wvl, nwavel_max, max_rs_idx, &
                           n_fitvar_rad, errstat)
  if (errstat /= 0) return

  call label_output_file (upper_case(molname), ecs_version_id, errstat)
  if (errstat /= 0) return

  ! FIXME: he5 output stuff to be removed once netcdf conversion is complete.
  !        netcdf output file creation occurs a bit later after some output dimensions
  !        have been determined
  if (yn_do_he5_output) then
    errstat = HE5_Init_Swath ( l2_filename, pge_swath_name, ntimes_rad, nxtrack_rad, CmETA )
    if (errstat /= 0) return

    CALL he5_initialize_datafields ( )
    errstat = HE5_Define_Fields ( pge_idx, pge_swath_name, ntimes_rad, nxtrack_rad, CmETA )
    if (errstat /= 0) return
  endif

  ! -----------------------------------------------------------------------------------
  ! If we are NOT using a radiance reference, then we need to read the
  ! swath line for radiance wavelength calibration. In this case, the
  ! value of RADIANCE_WAVCAL_LNUMS is the selected radiance reference line.
  ! For OMI, the first line worked best.
  ! For TEMPO, the first line can be clipped by the limb of the earth,
  !            but the last line may be (is?) fortuitously unaffected by that.
  ! Obviously, a more general algorithm would look at all the radiances first and
  ! pick the very best line, but since we're trying to avoid reading all the data
  ! at once, we try to get by with a simple choice.
  ! -----------------------------------------------------------------------------------
  IF ( .NOT. yn_radiance_reference ) THEN
    ntimes_loop = 1 ! The number of scan lines to read
    iline = radiance_wavcal_lnums(2)
    ! Get NTIMES_LOOP radiance lines.
    ! Note: omi_read_radiance_lines sets the global omi_nwav_rad
    CALL omi_read_radiance_lines (&
         l1b_rad_filename, iline, nxtrack_rad, ntimes_loop, &
         nwavel_rad, errstat )
    if (errstat /= 0) return
  END IF
  
  ! -----------------------------------------------------
  ! Across-track loop for radiance wavelength calibration
  ! -----------------------------------------------------

  call tell_log (1, 'omi_fitting: calling xtrack_radiance_wvl_calibration')
  omi_radcal_xflag = i2_missval
  CALL xtrack_radiance_wvl_calibration (first_wc_pix, last_wc_pix, &
       nxtrack_rad, n_max_rspec, &
       save_radcal_wvl, save_radcal_resid, &
       errstat )
  if (errstat /= 0) return
  
  if (yn_diagnostic_run) then
    call write_solar_wavecal_diagnostics (nwavel_max, nxtrack_rad, &
         save_solcal_wvl, save_solcal_resid, &
         errstat)
    call write_radiance_wavecal_diagnostics (nwavel_max, nxtrack_rad, &
         save_radcal_wvl, save_radcal_resid, &
         errstat)
    if (errstat /= 0) return
    ! FIXME: in diagnostic mode, there's a memory leak
    ! if we don't make it to this deallocate statement.
    deallocate (save_solcal_wvl, save_solcal_spec, save_solcal_resid, &
         save_radcal_wvl, save_radcal_resid, stat=errstat)
    if (errstat /= 0) then
       call tell_error(tell_malloc_error, "omi_fitting: deallocate failed", &
                       errstat)
       return
     endif
    if (errstat /= 0) return
  endif

  ! --------------------------------------------------------------
  ! Terminate on not having any cross-track pixels left to process
  ! --------------------------------------------------------------
  IF ( ALL ( omi_cross_track_skippix ) ) THEN
    call tell_log (0, "omi_fitting: no valid cross-track positions to process")
    !CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_NOPIXEL, &
    !  modulename, vb_lev_default, errstat )
  else

    if (.not.yn_radiance_reference) then
      call init_prefit_files (pge_idx, ntimes_rad, nxtrack_rad, errstat)
      if (errstat /= 0) return
    endif
    
    ! ---------------------------------------------------------------------
    ! If we are using a radiance reference AND want to remove the target
    ! gas from it (important for BrO, for example), we have to run through
    ! all spectra that go into the reference, compute the average column,
    ! and then remove that from the averaged radiance reference spectrum
    ! (owing to the fact that the average of hundreds of OMI spectra still
    !  doesn't produce a decently fitted column).
    ! ---------------------------------------------------------------------
    IF ( yn_radiance_reference ) THEN

      iline = SUM ( radiance_reference_lnums(1:2) ) / 2
      IF ( iline < 0 .OR. iline > ntimes_rr ) iline = ntimes_rr / 2

      first_pix = omi_xtrpix_range_rr(iline,1)
      last_pix  = omi_xtrpix_range_rr(iline,2)

      call tell_log (1, 'omi_fitting:  calling xtrack_radiance_reference_loop (1)')
      CALL xtrack_radiance_reference_loop ( &
           yn_remove_target, & ! note: yn_remove_target=TRUE here
           nxtrack_rr, nwavel_rr, first_pix, last_pix, pge_idx, errstat )
      if (errstat /= 0) return

      ! -------------------------------------------------------------
      ! Write the output from solar/earthshine wavelength calibration
      ! and radiance reference to file. The latter results will be
      ! overwritten in the call to XTRACK_RADIANCE_REFERENCE_LOOP
      ! below, hence we need to write them out here.
      ! -------------------------------------------------------------
      if (yn_do_he5_output) then
        CALL he5_write_wavcal_output ( nxtrack_rad, first_pix, last_pix, &
             errstat)
      endif
      call write_wavcal_output (result_vars, nxtrack_rad, errstat)
      if (errstat /= 0) return

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
    allocate (l1b_rad_latitudes (1:nxtrack_rad, 0:ntimes_rad-1), &
         stat=locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "omi_fitting: allocate failed", &
           errstat)
      return
    endif

    if (.not. yn_gems) then !TEMPO
      CALL read_latitude ( &
         TRIM(ADJUSTL(l1b_rad_filename)), &
         TRIM(ADJUSTL(omi_radiance_swathname)), &
         0, ntimes_rad, l1b_rad_latitudes, errstat)
    else !GEMS
       call gems_read_latitude (trim(adjustl(l1b_rad_filename)), &
           0, ntimes_rad, l1b_rad_latitudes, errstat)
    endif
    if (errstat /= 0) return

    ! -----------------------------------------------------------------
    ! Now we enter the on-line computation of the common mode spectrum.
    ! -----------------------------------------------------------------
    IF ( yn_common_iter ) THEN

      ! ----------------------------------------------------------
      ! Set the logical YN array that determines which swath lines
      ! will be used in the common mode
      ! ----------------------------------------------------------
      if (yn_disable_omi_features) then
        is_common_range = .true.
      else
        is_common_range = .FALSE.
        CALL find_swathline_range ( &
             TRIM(ADJUSTL(l1b_rad_filename)), &
             TRIM(ADJUSTL(omi_radiance_swathname)), &
             ntimes_rad, nxtrack_rad, l1b_rad_latitudes, &
             common_latrange(1:2), is_common_range, errstat )
      endif

      ! ------------------------------------------
      ! Interface to the loop over all swath lines
      ! ------------------------------------------
      call tell_log (1, 'omi_fitting: calling swathline_loops (common mode)-------------------------------')
      CALL swathline_loops ( &
           pge_idx, rpt_rad, n_max_rspec, &
           is_common_range, &
           omi_xtrpix_range, &
           .FALSE., -1, &
           .TRUE., errstat )

      if (errstat /= 0) return

      ! ---------------------------------------------------
      ! Set the index value of the Common Mode spectrum and
      ! assign values to the fitting parameter arrays
      ! ---------------------------------------------------
      CALL finalize_common_mode (nxtrack_rad)

      ! -------------------------------------------
      ! Write the just computed common mode to file
      ! -------------------------------------------
      IF ( yn_diagnostic_run ) then
        write(logmsg,'(a,i4)')'omi_fitting: writing out common mode, n_comm_wvl=',n_comm_wvl
        call tell_log (1, logmsg)
        if (yn_do_he5_output) then
          CALL he5_write_common_mode ( nxtrack_rad, n_comm_wvl, errstat )
        endif
        call write_common_mode (nxtrack_rad, n_comm_wvl, common_mode_spec, &
             errstat)
        if (errstat /= 0) return
      endif

    END IF

    ! ----------------------------------------------------------
    ! Now into the proper fitting, with or without common mode.
    ! ----------------------------------------------------------

    ! --------------------------------------------------------------
    ! Radiance Reference Fit: Only if we have selected to remove the
    ! target gas from the radiance reference
    ! -------------------------------------------------------------------
    !IF ( radiance_wavcal_lnums(1) >= 0 ) THEN
    IF (yn_radiance_reference .and. yn_remove_target) then

      ! --------------------------------
      ! The number of scan lines to read
      ! --------------------------------
      ntimes_loop = 1
      iline = radiance_wavcal_lnums(1)

      ! ------------------------------
      ! Get NTIMES_LOOP radiance lines
      ! ------------------------------
      CALL omi_read_radiance_lines ( &
           l1b_rad_filename, iline, nxtrack_rad, ntimes_loop, nwavel_rad, &
           errstat )
      if (errstat /= 0) return

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
      call tell_log (1, 'omi_fitting:  calling xtrack_radiance_reference_loop (2)')
      CALL xtrack_radiance_reference_loop ( &
           .FALSE., nxtrack_rr, nwavel_rr, first_pix, last_pix, pge_idx, &
           errstat )

    END IF

    ! -------------------------------------------------------------------
    ! Output of fit results for solar and radiance wavelength calibration
    ! but ONLY if we haven't done it already (see above). The Radiance
    ! Reference values under YN_REMOVE_TARGET settings carry valuable
    ! information only BEFORE the target has been removed.
    ! -------------------------------------------------------------------
    IF ( .NOT. (yn_radiance_reference .AND. yn_remove_target) ) THEN
      if (yn_do_he5_output) then
        CALL he5_write_wavcal_output ( nxtrack_rr, first_pix, last_pix, &
             errstat )
        if (errstat /= 0) return
      endif
      call write_wavcal_output (result_vars, nxtrack_rr, errstat)
      if (errstat /= 0) return

    END IF

    ! ----------------------------------------------------------
    ! Set the logical YN array that determines which swath lines
    ! will be processed. Unless we have constrained either swath
    ! line numbers or the latitude range, this can be set to
    ! .TRUE. universically.
    ! ----------------------------------------------------------
    ! First, set the range of swath lines to process
    ! ----------------------------------------------
    if (yn_disable_omi_features) then
      first_line = 0
      last_line = ntimes_rad - 1
      if ((0 <= pixnum_lim(1)) &
           .and. (pixnum_lim(1) < ntimes_rad)) first_line = pixnum_lim(1)
      if ((first_line <= pixnum_lim(2)) &
           .and. (pixnum_lim(2) < ntimes_rad)) last_line = pixnum_lim(2)
      do_radfit_range(:) = .false.
      do_radfit_range(first_line:last_line) = .true.
    else
      first_line = 0  ;  last_line = ntimes_rad-1
      IF ( pixnum_lim(1) > 0 ) first_line = MIN(pixnum_lim(1), last_line)
      IF ( pixnum_lim(2) > 0 ) last_line = MAX( MIN(pixnum_lim(2), last_line),&
           first_line )

      do_radfit_range = .FALSE.
      IF ( first_line         > 0           .OR. &
           last_line          < ntimes_rad-1 .OR. &
           radfit_latrange(1) > -90.0_r4    .OR. &
           radfit_latrange(2) < +90.0_r4           ) THEN

        IF ( radfit_latrange(1) > -90.0_r4    .OR. &
             radfit_latrange(2) < +90.0_r4           ) THEN
          CALL find_swathline_range ( &
               TRIM(ADJUSTL(l1b_rad_filename)), &
               TRIM(ADJUSTL(omi_radiance_swathname)), &
               ntimes_rad, nxtrack_rad, l1b_rad_latitudes, &
               radfit_latrange(1:2), do_radfit_range, errstat )
        ELSE
          do_radfit_range = .TRUE.
          IF ( first_line > 0 ) do_radfit_range(0:first_line-1) = .FALSE.
          IF ( last_line  < ntimes_rad-1 ) &
               do_radfit_range(last_line+1:ntimes_rad-1) = .FALSE.
        END IF
      ELSE
        do_radfit_range = .TRUE.
      END IF
    endif

    deallocate (l1b_rad_latitudes, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "omi_fitting: deallocated l1b_rad_latitudes failed", errstat)
      return
    endif

    ! ------------------------------------------
    ! Interface to the loop over all swath lines
    ! ------------------------------------------
    call tell_log (1, 'omi_fitting: calling swathline_loops (radiances)----------------------------')
    CALL swathline_loops ( &
         pge_idx, rpt_rad, n_max_rspec,     &
         do_radfit_range,                           &
         omi_xtrpix_range,                      &
         .FALSE., -1,                       &
         .FALSE., errstat)
    if (errstat /= 0) return

  endif ! if (all (omi_cross_track_skippix))

  ! ----------------------------------------------
  ! This subroutine completes the following tasks:
  !    (1) Compute pixel geololcation corners
  !    (2) Compute AMFs
  !    (3) Apply cross-track destriping correction
  ! ---------------------------------------
    call tell_log (1, 'omi_fitting:  calling omi_pge_postprocess ----------------------------')
    CALL omi_pge_postprocess ( &
       l1b_rad_filename, omi_radiance_swathname, pge_idx, &
       ntimes_rad, nxtrack_rad, &
       do_radfit_range, omi_xtrpix_range, &
       omi_is_szoom, n_max_rspec, fit_stats, errstat )
    if (errstat /= 0) return

  call tell_log (1, 'omi_fitting:  writing output...')
  ! ---------------------
  ! Write some attributes
  ! ---------------------
  call write_fitting_statistics (fit_stats, correlation_names, n_fitvar_rad, errstat)
  if (errstat /= 0) return
  if (yn_do_he5_output) then
    errstat = he5_write_global_attributes (fit_stats)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "omi_fitting: he5_write_global_attributes failed", &
           errstat)
      return
    endif
  endif

  if (yn_do_he5_output) then
    errstat = he5_write_swath_attributes ( pge_idx )
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "omi_fitting: he5_write_swath_attributes failed", &
           errstat)
      return
    endif
  endif

  if (yn_do_he5_output) then
    errstat = he5_set_field_attributes ( pge_idx )
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "omi_fitting: he5_set_field_attributes failed", &
           errstat)
      return
    endif
  endif
  ! -----------------
  ! Close output file
  ! -----------------
  call close_output_file (errstat)
  if (errstat /= 0) return

  if (yn_do_he5_output) then
    errstat = he5_close_output_file ( pge_idx)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "omi_fitting: he5_close_output_file failed", &
           errstat)
      return
    endif
  endif

  ! --------------------------------------------------------
  ! Write ODL metadata to go with netCDF, if switch is set
  ! --------------------------------------------------------
  ! Note this occurs now because postprocessing has the whole lon, lat grid
  ! available, whereas main retrieval seems to only ever have blocks
  if (yn_wrt_odl) then
    errstat = write_odl_metadata (l1b_rad_filename, l2_filename_netcdf, &
         ecs_version_id, nxtrack_rad-1, ntimes_rad-1, lun_input, n_lun_inp)
    if (errstat /= 0) then
      call tell_error(tell_io_write_error, "failed writing ODL metadata", &
           errstat)
      return
    endif
  endif

  ! ----------------------------------------
  ! Get Metadata from MCF, PCF, and L1B file
  ! ----------------------------------------
  CALL check_metadata_consistency ( errstat )
  if (errstat /= 0) then
    call tell_error (tell_application_error, &
         "omi_fitting: check_metadata_consistency failed", &
         errstat)
    return
  endif

  CALL set_l2_metadata ( pge_idx, errstat )
  if (errstat /= 0) then
    call tell_error (tell_io_write_error, &
         "omi_fitting: set_l2_metadata failed", &
         errstat)
    return
  endif

END SUBROUTINE omi_fitting

END MODULE
