MODULE pcf_file_module

  use tell_module
  implicit none
  character (LEN=*), parameter :: modulename = 'read_pcf_file'
  private modulename

CONTAINS

  subroutine do_pgs_get_reference (lun, lunstr, msgtype, &
       severity, outfile, errstat)

    ! This is a wrapper around the C function PGS_PC_GetReference.
    ! PGS_PC_GetReference has the prototype:
    !
    !   int PGS_PC_GetReference (int lun, int *version, char *str)
    !
    ! Here, version on input specifies the version to get.  On output
    ! it is set to the number of versions remaining.  Upon return,
    ! str will be set to the HD5 locater.
    !
    ! Note that version must be passed by reference.
    use OMSAO_precision_module, ONLY: i4, C_INT
    use OMSAO_parameters_module, ONLY: MAX_STR_LEN
    USE OMSAO_errstat_module, ONLY: pge_errstat_warning, error_check, &
         PGSd_PC_VALUE_LENGTH_MAX, ERR_SEP, vb_lev_default, &
         PGS_SMF_MASK_LEV_S, pge_errstat_error, &
         vb_lev_default

    implicit none

    integer (kind=i4), intent(IN) :: lun, severity, msgtype
    character (len=*), intent(IN) :: lunstr
    character (len=*), intent(OUT) :: outfile
    integer (kind=i4), intent(OUT) :: errstat

    ! locals
    character (len=PGSd_PC_VALUE_LENGTH_MAX) :: file
    integer (kind=C_INT) :: version, locerr
    integer (kind=i4) :: filelen
    character (len=MAX_STR_LEN) :: errstr
    integer (kind=C_INT), external :: PGS_PC_GetReference, &
         pgs_smf_teststatuslevel

    version = 1; file = " "
    locerr = PGS_PC_GetReference (lun, version, file)
    ! Now version is set to the number of versions remaining.

    locerr = pgs_smf_teststatuslevel (locerr)

    write (errstr,'(a,i9,a,a,a,i4,a,a)')"pgs_get_reference: lun=", lun, ", lunstr=", TRIM(lunstr), &
         " version=", version, " file=", trim(file)
    call tell_log (1, errstr)

    file = ADJUSTL(file)
    filelen = LEN_TRIM(file)

    if (( locerr /= PGS_SMF_MASK_LEV_S ).OR.(filelen == 0)) then

      file = ""
      if ((filelen == 0) &
           .and. (locerr == PGS_SMF_MASK_LEV_S)) &
           locerr = locerr + 1

      write (errstr, '(i0)') lun
      errstr = modulename//ERR_SEP//lunstr//" "//TRIM(ADJUSTL(errstr))

      select case (severity)
      case (pge_errstat_warning)
        call error_check (0, 1, severity, msgtype, &
             errstr, vb_lev_default, errstat)

      case (pge_errstat_error)
        call error_check (locerr, PGS_SMF_MASK_LEV_S, severity, &
             msgtype, errstr, vb_lev_default, errstat)

      case default
        call error_check (locerr, PGS_SMF_MASK_LEV_S, severity, &
             msgtype, errstr, vb_lev_default, errstat)
      end select
    else
      outfile = trim(adjustl(file))
    endif

  end subroutine do_pgs_get_reference

  SUBROUTINE read_pcf_file (pge_idx, pge_name, errstat )

    USE OMSAO_precision_module
    USE OMSAO_indices_module
    USE OMSAO_errstat_module
    USE metadata_tools, ONLY: &
         pcf_granule_s_time,  pcf_granule_e_time, &
         n_mdata_int, mdata_integer_fields, mdata_integer_values
    USE OMSAO_parameters_module,     ONLY: zerospec_string, str_missval
    USE OMSAO_he5_module,            ONLY: &
         pge_swath_name, process_level, instrument_name, pge_version
    !USE OMSAO_omidata_module,        ONLY: l1b_radiance_esdt
    USE OMSAO_variables_module,      ONLY: &
         verb_thresh_lev, orbit_number, &
         ecs_version_id, l1b_rad_filename, l1b_irrad_filename, l2_filename,  &
         static_input_fnames, Have_AMF_Table, omi_slitfunc_fname,            &
         solcal_filename, OMSAO_I0_filename, voc_amf_filenames,              &
         refspecs_original, OMSAO_refseccor_filename, OMSAO_OMLER_filename,  &
         OMSAO_refseccor_cld_filename,                                       &
         solcal_cache_mode, solcal_source
    USE OMSAO_wfamf_module,         ONLY: wfamf_table_lun, climatology_lun,  &
         OMSAO_wfamf_table_filename, OMSAO_climatology_filename, &
         num_met_luns, meteorology_lun, OMSAO_meteorology_filename
    USE control_module, ONLY: read_fitting_control_file
    USE sao_pge_utils, ONLY: get_pge_ident
    USE strutils, ONLY: remove_quotes
    !use ctrlvars, only: yn_wrt_odl

    IMPLICIT NONE

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=I4), INTENT(OUT) :: pge_idx
    CHARACTER (LEN=12), INTENT(OUT) :: pge_name
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                      :: i, j

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: pge_error_status, estat

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER   (KIND=i4),      EXTERNAL :: &
         pgs_pc_getnumberoffiles, pgs_pc_getconfigdata, &
         pgs_smf_teststatuslevel

    pge_error_status = pge_errstat_ok

    ! --------------------------------------------------------------------
    ! Read all ConfigData from PCF and set variables associated with them.
    ! Some, like MoleculeID, might lead to PGE termination if in error.
    ! --------------------------------------------------------------------
    DO i = 1, n_config_luns
      estat = &
           PGS_PC_GetConfigData ( config_lun_array(i), config_lun_values(i) )
      estat = PGS_SMF_TestStatusLevel(estat)
      CALL error_check ( estat, PGS_SMF_MASK_LEV_S, pge_errstat_warning, &
           OMSAO_W_GETLUN, &
           modulename//f_sep//TRIM(ADJUSTL(config_lun_strings(i))), &
           vb_lev_default, pge_error_status )
      ! -------------------------------------------------------------
      ! Remove any quotes (") that might have made in into the string
      ! -------------------------------------------------------------
      CALL remove_quotes ( config_lun_values(i) )

      SELECT CASE ( config_lun_array(i) )
        ! --------------------------------------------------------------------
        ! Get the verbosity threshold. The setting of this variable determines
        ! how "chatty" in terms of screen output the PGE will be. If the READ
        ! was in error, proceed with DEBUG level
        ! --------------------------------------------------------------------
      CASE ( verbosity_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "1"
        READ (config_lun_values(i), '(I1)') verb_thresh_lev
        ! ----------------------------------------------------------------------
        ! Get the orbit number. This will be checked against L1B MetaData later.
        ! Note that the source for OrbitNumber is the L1B file, NOT the PCF. But
        ! this only matters in case they are different, which will be checked
        ! later. If the ARE different, then the L1B value will take precidence.
        ! ----------------------------------------------------------------------
      CASE ( orbitnumber_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "00000"
        READ (config_lun_values(i), '(I5)') orbit_number

        ! --------------------------------------------------------------------
        ! Get the solar irradiance calibration mode. This is used to either
        ! 1. "none": run solar irradiance slit and wavelength calibration in trace gas code
        !    and then proceed to spectral fitting;
        ! 2. "save": perform calibration only then save results; or 
        ! 3. "read": previously saved calibration file and then perform spectral fitting.
        ! --------------------------------------------------------------------
      CASE ( solcal_mode_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "none"
        solcal_cache_mode = TRIM(ADJUSTL(config_lun_values(i)))
        ! -----------------------------------------------------------------
        ! Get the source (solar irradiance file or I0 irradiance) used in
        ! solar calibration; options are
        ! 1. solar_irradiance (default)
        ! 2. I0_irradiance
        ! -----------------------------------------------------------------
      CASE ( solcal_source_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "solar_irradiance"
        solcal_source = TRIM(ADJUSTL(config_lun_values(i)))
        ! -----------------------------------------------------------------
        ! Get Granule Start and End time. This is provided by the Scheduler
        ! and has to be checked against the L1B MetaData fields for
        ! RangeBeginningTime and RangeEndingTime.
        ! -----------------------------------------------------------------
      CASE ( granule_s_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) &
             config_lun_values(i) = "T00:00:00.000Z"
        pcf_granule_s_time = TRIM(ADJUSTL(config_lun_values(i)))
      CASE ( granule_e_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) &
             config_lun_values(i) = "T00:00:00.000Z"
        pcf_granule_e_time = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the Process Level
        ! ------------------------------------------------------------------
      CASE ( proclevel_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        process_level = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the PGE Version
        ! ------------------------------------------------------------------
      CASE ( pge_version_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        pge_version = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the Instrument Name
        ! ------------------------------------------------------------------
      CASE ( instrument_name_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        instrument_name = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the HE5 Swath Name; set to MISSING VALUE if it can't be found.
        ! ------------------------------------------------------------------
      CASE ( swathname_lun )
        IF ( estat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        pge_swath_name = TRIM(ADJUSTL(config_lun_values(i)))

        ! ----------------------------------------------------------------------
        ! Get ECS Collection version number and assign it to a Metadata variable
        ! ----------------------------------------------------------------------
      CASE ( versionid_lun )
        READ (config_lun_values(i), '(I1)') ecs_version_id
        getidx: DO j = 1, n_mdata_int
          IF ( mdata_integer_fields(3,j) == "pcf"              .AND. ( &
               INDEX(mdata_integer_fields(1,j), 'VERSIONID') /= 0 .OR. &
               INDEX(mdata_integer_fields(1,j), 'VersionID') /= 0     )  ) THEN
            mdata_integer_values(j) = ecs_version_id
            EXIT getidx
          END IF
        END DO getidx

        ! -------------------------------------------------------------------------
        ! Get the SAO PGE Name string. This string is converted to the PGE
        ! index number (10, 11, 12), which in turn is used to access LUNs and other
        ! elements that are PGE specific. All of this has to be done here, because
        ! we require the PGE index to identify LUNs and other PGE specific items.
        ! >>> We MUST know this, or else we can't execute the PGE <<<
        ! -------------------------------------------------------------------------
      CASE ( pge_molid_lun )
        CALL error_check ( estat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, &
             OMSAO_F_GETLUN, &
             modulename//f_sep//TRIM(ADJUSTL(config_lun_strings(i))), &
             vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_error ) then
          call tell_error(tell_io_read_error, &
               "read_pcf_file: failed to read molecule ID LUN", errstat)
          return
        endif

        ! ---------------------------------------------------------------
        ! Translate molecule string into index (required for LUN look-up)
        ! ---------------------------------------------------------------
        estat = pge_errstat_ok
        CALL get_pge_ident (config_lun_values(i), pge_idx, pge_name, estat)
        CALL error_check ( estat, pge_errstat_ok, pge_errstat_fatal, &
             OMSAO_F_GET_MOLINDEX, modulename, vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_fatal ) then
          call tell_error(tell_io_read_error, &
               "read_pcf_file: failed to read molecule index", errstat)
          RETURN
        endif
        CALL error_check ( 0, 1, pge_errstat_ok, OMSAO_S_GET_MOLINDEX, &
             TRIM(ADJUSTL(pge_name)), vb_lev_stmdebug, estat )
      END SELECT
    END DO

    estat = pge_errstat_ok

    ! ---------------------------------------------------------
    ! Static input file with tabulated OMI slit function values
    ! ---------------------------------------------------------
    call do_pgs_get_reference (omi_slitfunc_lun, "OMI_SLITFUNC_LUN", &
         OMSAO_F_GETLUN, pge_errstat_fatal, &
         omi_slitfunc_fname, pge_error_status)
    if (pge_error_status >= pge_errstat_error) then
      call tell_error(tell_io_read_error, &
           "read_pcf_file: failed to read slitfunction LUN", errstat)
      return
    endif

    ! ---------------------------------------------------------
    ! Get input name for on-orbit solar irradiance calibration file.
    ! This file is generated by trace gas code and is either used as an output
    ! or input filename, depending on the solar_cache_mode in the pipeline.
    ! ---------------------------------------------------------
    call do_pgs_get_reference (solcal_fname_lun, "SOLCAL_FNAME_LUN", &
         OMSAO_F_GETLUN, pge_errstat_fatal, &
         solcal_filename, pge_error_status)
    if (pge_error_status >= pge_errstat_error) then
      call tell_error(tell_io_read_error, &
           "read_pcf_file: failed to read solar irradiance calibration LUN", errstat)
      return
    endif

    ! ------------------------------------------------------------------
    ! Read names of static input files from PCF. The first entry is the
    ! algorithm control file (initial fitting values, etc.), the other
    ! entries are reference spectra.
    ! ------------------------------------------------------------------
    static_input_fnames = zerospec_string
    DO i = icf_idx, max_rs_idx

      ! JED note: before I modified this, the call to error_check was:
      !
      !   error_check ( 0, 1, pge_errstat_fatal, OMSAO_F_GETLUN, &
      !          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
      !          vb_lev_default, pge_error_status )
      !
      ! But, passing error_check was first arg=0 will produce pge_errstat_ok for
      ! pge_error_status.  So, I am going to make this a warning.

      call do_pgs_get_reference (pge_static_input_luns(i), "PGE_STATIC_INPUT_LUN", &
           OMSAO_W_GETLUN, pge_errstat_warning, &
           static_input_fnames(i), pge_error_status)
      if (pge_error_status >= pge_errstat_error) then
        call tell_error(tell_io_read_error, &
             "read_pcf_file: failed to read static input LUNs", errstat)
        return
      endif
      if (0 /= index (trim(static_input_fnames(i)), zerospec_string)) &
           static_input_fnames(i) = zerospec_string
      !write (*,*) "static_input_fnames(", i, ")=", trim(static_input_fnames(i))
    end do

    ! ------------------------------------------------------
    ! Save file names with reference spectra to global array
    ! ------------------------------------------------------
    refspecs_original(1:max_rs_idx)%FileName = static_input_fnames(1:max_rs_idx)

    ! ----------------------------------------------------
    ! Read fitting conrol parameters from input file; this
    ! returns L1B_RADIANCE_ESDT, which determines the
    ! ingestion of the L1b radiance file.
    ! ----------------------------------------------------
    estat = pge_errstat_ok
    CALL read_fitting_control_file ( pge_idx, errstat ) ! l1b_radiance_esdt, errstat )
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, &
           "read_pcf_file: reading fitting control file", errstat)
      return
    endif
    !CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
    !  modulename//f_sep//"READ_FITTING_CONTROL_FILE.", vb_lev_default, pge_error_status )
    !IF ( pge_error_status >= pge_errstat_fatal ) RETURN

    ! -----------------------------
    ! Read Irradiance L1B file name
    ! -----------------------------
    call do_pgs_get_reference (l1b_irradiance_lun, "L1B_IRRADIANCE_LUN", &
         OMSAO_F_GETLUN, pge_errstat_fatal, &
         l1b_irrad_filename, pge_error_status)
    if (pge_error_status >= pge_errstat_error) then
      call tell_error(tell_io_read_error, &
           "read_pcf_file: failed to read L1B irradiance LUN", errstat)
      return
    endif
    ! ---------------------------
    ! Read Radiance L1B file name
    ! ---------------------------
    call do_pgs_get_reference (l1b_radiance_lun, "L1B_RADIANCE_LUN", &
         OMSAO_F_GETLUN, pge_errstat_fatal, &
         l1b_rad_filename, pge_error_status)
    if (pge_error_status >= pge_errstat_error) then
      call tell_error(tell_io_read_error, &
           "read_pcf_file: failed to read L1B radiance LUN", errstat)
      return
    endif

    ! -------------------------------------------------------------------------
    ! Read name of AMF table file(s). Remember that a missing AMF table is
    ! not a fatal problem, since in that case the slant columns will be written
    ! to the output file.
    ! -------------------------------------------------------------------------
    Have_AMF_Table = .TRUE.
     DO i = 1, n_voc_amf_luns
       call do_pgs_get_reference (voc_amf_luns(i), "PGE_STATIC_INPUT_LUN", &
            OMSAO_W_GETLUN, pge_errstat_warning, &
            voc_amf_filenames(i), pge_error_status)
       Have_AMF_Table = Have_AMF_Table.and.(len_trim(voc_amf_filenames(i)) > 0)
     END DO

    ! -------------------------------------------------------------------------
    ! gga. The new implementation of the wf amf is intended to be molecule inde
    ! pendent. Therfore I read the file names outside the CASE stament.
    ! -------------------------------------------------------------------------
    ! wavelength dependent AMF table file
    ! -----------------------------------
    call do_pgs_get_reference (wfamf_table_lun, "PGE_STATIC_INPUT_LUN", &
         OMSAO_W_GETLUN, pge_errstat_warning, &
         OMSAO_wfamf_table_filename, pge_error_status)

    ! ----------------
    ! Climatology
    ! ----------------
    call do_pgs_get_reference (climatology_lun, "PGE_STATIC_INPUT_LUN", &
         OMSAO_W_GETLUN, pge_errstat_warning, &
         OMSAO_climatology_filename, pge_error_status)

    !-----------------
    ! Meteorology
    !-----------------
    do i=1, num_met_luns
      call do_pgs_get_reference (meteorology_lun(i), "PGE_STATIC_INPUT_LUN", &
                                 OMSAO_W_GETLUN, pge_errstat_warning, &
                                 OMSAO_meteorology_filename(i), pge_error_status)
    enddo

    ! --------------------------------------------------------
    ! Read name of file with I0 irradiance-replacement spectra
    ! --------------------------------------------------------
    call do_pgs_get_reference (OMSAO_I0_lun, "PGE_STATIC_INPUT_LUN", &
         OMSAO_W_GETLUN, pge_errstat_warning, &
         OMSAO_I0_filename, pge_error_status)

    ! ------------------------------------------------------------
    ! Read name of file with GEOS-Chem background Reference Sector
    ! concentrations !gga
    ! ------------------------------------------------------------
    call do_pgs_get_reference (OMSAO_refseccor_lun, "PGE_STATIC_INPUT_LUN", &
         OMSAO_W_GETLUN, pge_errstat_warning, &
         OMSAO_refseccor_filename, pge_error_status)

    ! ---------------------------------------------------------
    ! Read name of file with the radiance reference clouds !gga
    ! ---------------------------------------------------------
    call do_pgs_get_reference (OMSAO_refseccor_cld_lun, "PGE_STATIC_INPUT_LUN", &
         OMSAO_W_GETLUN, pge_errstat_warning, &
         OMSAO_refseccor_cld_filename, pge_error_status)

    ! ----------------------------------
    ! Read name of OMLER albedo file gga
    ! ----------------------------------
    call do_pgs_get_reference (OMSAO_OMLER_lun, "PGE_STATIC_INPUT_LUN", &
         OMSAO_W_GETLUN, pge_errstat_warning, &
         OMSAO_OMLER_filename, pge_error_status)

    ! ---------------------------
    ! Read name of L2 output file
    ! ---------------------------
    call do_pgs_get_reference (pge_l2_output_lun, "PGE_L2_OUTPUT_LUN", &
         OMSAO_F_GETLUN, pge_errstat_fatal, &
         l2_filename, pge_error_status)
    IF ( pge_error_status >= pge_errstat_error ) then
      call tell_error(tell_io_read_error, &
           "read_pcf_file: failed to read L2 output LUN", errstat)
      return
    endif

    !-------------------------------------------------------
    ! if writing metadata with netCDF, get namelist filename
    !-------------------------------------------------------
    !if (yn_wrt_odl) then
    !  call do_pgs_get_reference (mdlist_LUN, "PGE_MDLIST_LUN", &
    !       OMSAO_F_GETLUN, pge_errstat_fatal, &
    !       mdlist_filename, pge_error_status)
    !  IF ( pge_error_status >= pge_errstat_error ) then
    !    call tell_error(tell_io_read_error, &
    !         "read_pcf_file: failed to read metadata list LUN", errstat)
    !    return
    !  endif
    !endif

    RETURN
  END SUBROUTINE read_pcf_file

END MODULE pcf_file_module
