MODULE OMSAO_Reference_sector_module

  ! ------------------------------------------------------------------
  ! This module defines variables associated with the Reference Sector
  ! Correction and contains the subroutines needed to apply it to the
  ! HCHO retrieval.
  ! ------------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i2, r8
  use tell_module

  IMPLICIT NONE

  PRIVATE i2, r8

  ! ------------------------------------------------------------------
  ! maxngrid: Parameter to define maximum size of arrays; set to 1000
  ! grid_lat: Latitudes of GEOS-Chem background levels from the Pacif
  !            ic, Hawaii is not included
  ! background_correction: Background level-reference_sector_concentra
  !                        tion
  ! background_level: Total column median obtained from the Radiance
  !                   Reference granule (molecules/cm2)
  ! Reference_sector_concentration: Total column (molecules/cm2) obtai
  !                                 ned using the GEOS-Chem climatolog
  !                                 y by D. Millet
  ! ngridpoints: Number of points in grid_lat and Reference_Sector_con
  !              centration
  ! ------------------------------------------------------------------

  INTEGER (KIND=i2), PRIVATE, PARAMETER               :: maxngrid = 1000
  REAL    (KIND=r8), PRIVATE, DIMENSION (maxngrid)    :: grid_lat!, background_correction, &
    !background_level
  REAL    (KIND=r8), PRIVATE, DIMENSION (maxngrid,12) :: Reference_sector_concentration
  INTEGER (KIND=i2), PRIVATE  :: ngridpoints

CONTAINS

  SUBROUTINE Reference_Sector_correction (ntimes, nxtrack, & !xtrange, lat
      saocol, saodco, saoamf, saomqf, pge_idx, n_max_rspec, &
      errstat)

    !USE l1bread, only: l1bread_radiance_info
    use l1bread_utils, only : read_l1_radiance_info
    USE OMSAO_precision_module, ONLY: i4
    USE OMSAO_variables_module, ONLY: Radiance_Paras_Type, &
      l1b_radref_filename, l1b_channel
    USE OMSAO_omidata_module, ONLY: omi_radiance_swathname
    use OMSAO_indices_module, only : pge_hcho_idx
    use output_tools, only : write_reference_sector_corrected_column
    use ctrlvars, only : yn_do_he5_output
    !USE OMSAO_errstat_module
    ! ---------------------------------------------------------------
    ! This subroutine is a wrapper for the Reference Background corre
    ! ction
    ! ---------------------------------------------------------------
    ! VARIABLE DECLARATION
    ! --------------------
    IMPLICIT NONE
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                   INTENT (IN) :: ntimes, nxtrack
    !INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),       INTENT (IN) :: xtrange
    !REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN) :: lat
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN) :: saoamf
    INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN) :: saomqf
    INTEGER (KIND=i4),                                   INTENT (IN) :: pge_idx, n_max_rspec
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN) :: saocol, saodco

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                                   INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                   :: nTimesRadRR, nXtrackRadRR, nWvlCCDrr, itrack, iline
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: int_saocol, int_saodco
    TYPE(Radiance_Paras_Type) :: rpt_rr
    REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE     :: mem_correction

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: locerrstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=27), PARAMETER :: modulename = 'Reference_sector_correction'

    if (errstat /= 0) return

    !locerrstat = pge_errstat_ok

    int_saocol   = saocol
    int_saodco   = saodco

    ! -------------------------------------------------
    ! Read the concentrations from the Reference Sector
    ! Total column molecules/cm2
    ! -------------------------------------------------
    CALL Read_reference_sector_concentration(errstat)

    ! ---------------------------------------------------
    ! Obtain dimensions of the Radiance Reference granule
    ! ---------------------------------------------------
    !CALL l1bread_radiance_info (l1b_radref_filename, l1b_channel, &
    !                            rpt_rr, errstat)
    call read_l1_radiance_info (l1b_radref_filename, l1b_channel, &
                                rpt_rr, errstat)
    if (errstat /= 0) return

    nTimesRadRR  = rpt_rr%ntimes
    nXtrackRadRR = rpt_rr%nxtrack
    nWvlCCDrr    = rpt_rr%nwavel_ccd
    omi_radiance_swathname = rpt_rr%swathname

    ALLOCATE (mem_correction(1:nXtrackRadRR,0:nTimesRadRR-1), stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "Reference_Sector_correction: allocate failed", errstat)
      return
    endif
    ! -------------------------------------------------------------------------
    ! We need to retrieve the concentrations for the Radiance Reference Granule
    ! to be able to compute the latitude dependent correction due to the Refere
    ! ce sector hypothesis. This is a bit pedestrian, in the case of the Radian
    ! ce granule being the same that the Radiance Reference granule we could sk
    ! ip this step but since time is not an issue. Let us do it anyway.
    ! -------------------------------------------------------------------------
    CALL Reference_Sector_radref_retrieval_and_median &
      (pge_idx, rpt_rr, n_max_rspec, errstat, mem_correction)

    ! -----------------------------------------------------------------------------
    ! Apply the correction to the slant columns and covert back to vertical columns
    ! All based in problems with the instrument in same pixels.
    ! -----------------------------------------------------------------------------
    DO itrack = 1, nxtrack
      DO iline = 0, ntimes-1

        IF (saomqf(itrack,iline) .NE. 0) CYCLE
        int_saocol(itrack,iline) = int_saocol(itrack,iline) * saoamf(itrack,iline)
        int_saocol(itrack,iline) = int_saocol(itrack,iline) - mem_correction(itrack,iline)
        int_saocol(itrack,iline) = int_saocol(itrack,iline) / saoamf(itrack,iline)

      END DO
    END DO

    DEALLOCATE (mem_correction, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "Reference_Sector_correction: deallocate failed", errstat)
      return
    endif

    ! ------------------------------------------
    ! Compute the back ground correction for the
    ! particular radiance reference granule
    ! ------------------------------------------
    ! CALL compute_background_correction(locerrstat)

    ! --------------------------------------------------
    ! Apply the already calculated background_correction
    ! to the saocol
    ! --------------------------------------------------
    ! CALL aply_background_correction(ntimes, nxtrack, int_saocol, & !int_saodco,
    !     saomqf, lat, locerrstat)

    ! ---------------------------------------------------
    ! Final step to output the results, the new reference
    ! sector corrected total columns
    ! ---------------------------------------------------
    if (yn_do_he5_output) then
      CALL he5_write_reference_sector_corrected_column &   ! FIXME <-- (to be removed)
        (pge_idx, &
         ntimes, nxtrack, int_saocol, &!int_saodco,
         locerrstat)
    endif
    if (pge_idx == pge_hcho_idx) then
      call write_reference_sector_corrected_column (nxtrack, ntimes, int_saocol, errstat)
      if (errstat /= 0) return
    endif

    !errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE Reference_Sector_correction

  SUBROUTINE Read_reference_sector_concentration(errstat)

    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN
    USE OMSAO_variables_module, ONLY: OMSAO_refseccor_filename
    USE OMSAO_indices_module,   ONLY: OMSAO_refseccor_lun
    USE OMSAO_precision_module, ONLY: i4
    USE OMSAO_errstat_module, only : pge_errstat_ok, pgs_smf_mask_lev_s, &
      pgsd_io_gen_rseqfrm

    IMPLICIT NONE
    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                                  INTENT (INOUT) :: errstat
    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER (KIND=i4), EXTERNAL :: &
      pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)           :: funit, igrid
    CHARACTER(LEN=1), PARAMETER :: hstr='#'
    LOGICAL                     :: file_header
    CHARACTER (LEN=MAX_STR_LEN)    :: header_line
    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: version, locerrstat, ios
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=35), PARAMETER :: modulename = 'Read_reference_sector_concentration'

    if (errstat /= 0) return

    locerrstat = pge_errstat_ok

    ! --------------------
    ! Initialize variables
    ! --------------------
    ngridpoints = 0_i2
    grid_lat(1:maxngrid) =  0.0_r8
    Reference_sector_concentration(1:maxngrid,1:12) = 0.0_r8

    ! -----------------------------------------
    ! Open Reference Sector concentrations file
    ! -----------------------------------------
    version = 1
    locerrstat = PGS_IO_GEN_OPENF ( OMSAO_refseccor_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
    if (locerrstat > pgs_smf_mask_lev_s) then
      call tell_error (tell_io_open_error, &
                       "opening ref. sector concentrations file: "// &
                       trim(OMSAO_refseccor_filename), errstat)
      return
    endif
    !CALL error_check ( &
    !  locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_REFSECCOR_FILE, &
    !  modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
    !IF (  errstat /= pge_errstat_ok ) RETURN

    ! ------------
    ! Reading file
    ! -----------------------------------------------------------------------------------
    ! Skip header lines. The header is done when the first character of the line is not #
    ! -----------------------------------------------------------------------------------
    file_header = .TRUE.
    skip_header: DO WHILE (file_header .EQV. .TRUE.)
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) header_line
      IF ( ios /= 0 ) THEN
        call tell_error (tell_io_read_error, &
                         "reading ref. sector concentrations file: "// &
                         trim(OMSAO_refseccor_filename), errstat)
        return
        !CALL error_check ( &
        !  ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
        !  modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
        !IF (  errstat /= pge_errstat_ok ) RETURN
      END IF
      IF (header_line(1:1) /= hstr) THEN
        file_header = .FALSE.
      ENDIF
    END DO skip_header

    ! -----------------------------------------------
    ! Read number of grid points. Variable via module
    ! -----------------------------------------------
    READ (UNIT=funit, FMT='(I5)', IOSTAT=ios) ngridpoints
    IF ( ios /= 0 ) THEN
      call tell_error (tell_io_read_error, &
                       "reading ref. sector concentrations file: "// &
                       trim(OMSAO_refseccor_filename), errstat)
      return
      !CALL error_check ( &
      !  ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
      !  modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
      !IF (  errstat /= pge_errstat_ok ) RETURN
    END IF

    ! ---------------------------------------------------------
    ! Read reference sector concentrations. Variable via module
    ! ---------------------------------------------------------
    DO igrid = 1, ngridpoints
      READ (UNIT=funit, FMT='(F6.2,2x,12(1x,E14.7))', IOSTAT=ios) grid_lat(igrid), &
        Reference_sector_concentration(igrid, 1:12)
      IF ( ios /= 0 ) THEN
        call tell_error (tell_io_read_error, &
                         "reading ref. sector concentrations file: "// &
                         trim(OMSAO_refseccor_filename), errstat)
        return
        !CALL error_check ( &
        !  ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
        !  modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
        !IF (  errstat /= pge_errstat_ok ) RETURN
      END IF
    END DO

    ! -----------------------------------------------
    ! Close monthly average file, report SUCCESS read
    ! -----------------------------------------------
    locerrstat = PGS_IO_GEN_CLOSEF ( funit )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
    if (locerrstat > pgs_smf_mask_lev_s) then
      call tell_error (tell_io_error, &
                       "closing ref. sector concentrations file: "// &
                       trim(OMSAO_refseccor_filename), errstat)
      return
    endif
    !CALL error_check ( &
    !  locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_REFSECCOR_FILE, &
    !  modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
    !IF ( errstat >= pge_errstat_error ) RETURN

  END SUBROUTINE Read_reference_sector_concentration

  SUBROUTINE Reference_Sector_radref_retrieval_and_median &
      (pge_idx, rpt_rr, n_max_rspec, errstat, mem_correction)

    USE OMSAO_precision_module, ONLY: i1, i4!, r4
    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, i2_missval, r8_missval, & ! , r4_missval
      MAX_STR_LEN
    USE OMSAO_wfamf_module,     ONLY: amf_calculation_bis
    USE OMSAO_variables_module, ONLY: OMSAO_refseccor_cld_filename, voc_amf_filenames, &
      Radiance_Paras_Type, common_latrange, l1b_rad_filename, l1b_radref_filename
    USE OMSAO_indices_module,   ONLY: voc_omicld_idx
    USE omi_pge_fitting_aux, ONLY: read_latitude, find_swathline_range, &
      compute_fitting_statistics, fitting_statistics_type
    use commonmode, only: finalize_common_mode
    USE omi_read_l1b_data, ONLY: omi_read_glint_ice_flags, omi_read_binning_factor
    USE swathline_loop, ONLY: swathline_loops
    !USE omi_pge_swathline_loop_memory, ONLY: omi_pge_swathline_loops_mem
    USE OMSAO_errstat_module, only : pge_errstat_ok, pge_errstat_error
    USE OMSAO_omidata_module, ONLY: omi_radiance_swathname, &
      retrieval_type, alloc_retrieval_type, dealloc_retrieval_type
    use ctrlvars, only: yn_disable_omi_features

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec
    TYPE (Radiance_Paras_Type), INTENT(IN) :: rpt_rr

    ! Output Variables
    REAL (KIND=r8), DIMENSION (1:rpt_rr%nxtrack, 0:rpt_rr%ntimes-1), &
      INTENT(OUT) :: mem_correction
    integer, intent(inout) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                           :: first_line, last_line
    INTEGER (KIND=i4), DIMENSION (0:rpt_rr%ntimes-1,2)            :: omi_xtrpix_range_rr
    LOGICAL,           DIMENSION (0:rpt_rr%ntimes-1)              :: radfitref_range_ok
    CHARACTER(LEN=MAX_STR_LEN)                                     :: l1b_rad_save_filename
    REAL    (KIND=r8), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: mem_amf
    !REAL    (KIND=r8), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: &
    !  mem_column_amount, mem_column_uncertainty, mem_rms
    !REAL    (KIND=r4), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: mem_latitude, mem_longitude, mem_sza, mem_vza, &
    !  mem_height
    !INTEGER (KIND=i2), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: mem_fit_flag, mem_xtrflg
    INTEGER (KIND=i2), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: mem_snow, mem_glint
    !INTEGER (KIND=i2), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: refmqf
    LOGICAL,           DIMENSION (0:rpt_rr%ntimes-1)              :: yn_szoom_rs, common_range_ok
    INTEGER (KIND=i1), DIMENSION (0:rpt_rr%ntimes-1)              :: binfac_rs
    LOGICAL                                                     :: yn_write
    INTEGER (KIND=i4) :: nTimesRadRR, nXtrackRadRR, nWvlCCDrr
    type (retrieval_type) :: rt
    type (fitting_statistics_type) :: ref_stats

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: locerrstat
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=64), PARAMETER :: modulename = &
    !  'Reference_Sector_radiance_reference_granule_retrieval' !JED fix

    if (errstat /= 0) return

    nTimesRadRR = rpt_rr%ntimes
    nXtrackRadRR = rpt_rr%nxtrack
    nWvlCCDrr = rpt_rr%nwavel_ccd

    locerrstat = pge_errstat_ok

    ! ---------------------------------
    ! A bit of mess with the file names
    ! ---------------------------------
    l1b_rad_save_filename = l1b_rad_filename
    l1b_rad_filename      = l1b_radref_filename

    ! -----------------------
    ! Variable initialization
    ! -----------------------
    call alloc_retrieval_type (rt, nXtrackRadRR, nTimesRadRR, locerrstat)
    if (locerrstat /= 0) return

    allocate (ref_stats % quality_flag(rpt_rr%nxtrack,0:rpt_rr%ntimes-1), &
              stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
                       "Reference_Sector_radiance_reference_granule_retrieval:  allocate failed", &
                       errstat)
      return
    endif

    !mem_column_amount      = r8_missval
    !mem_column_uncertainty = r8_missval
    ! ---------------------------------------------------
    ! mem_correction needs to be initialized to 0.0 since
    ! we want no correction if any of the two reference
    ! sector retrievals does fail.
    ! ---------------------------------------------------
    mem_correction         = 0.0
    ! ---------------------------------------------------
    mem_amf                = r8_missval
    !mem_rms                = r8_missval
    !mem_latitude           = r4_missval
    !mem_longitude          = r4_missval
    !mem_sza                = r4_missval
    !mem_vza                = r4_missval
    !mem_fit_flag           = i2_missval
    ref_stats % quality_flag  = i2_missval
    !mem_xtrflg             = i2_missval
    mem_snow               = i2_missval
    mem_glint              = i2_missval
    !mem_height             = r4_missval

    ! -----------------------------------------------------
    ! I want to perform the retrieval for the whole granule
    ! -----------------------------------------------------
    first_line = 0  ;  last_line = nTimesRadRR-1
    radfitref_range_ok = .TRUE.
    omi_xtrpix_range_rr(0:nTimesRadRR-1,1) = 1
    omi_xtrpix_range_rr(0:nTimesRadRR-1,2) = nXtrackRadRR

    ! ---------------------------------
    ! Read radiance reference latitudes
    ! ---------------------------------
    CALL read_latitude ( &
      TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
      0, nTimesRadRR, rt%latitude(1:nXtrackRadRR,0:nTimesRadRR-1), &
      errstat)
    if (errstat /= 0) &
      return

    if (yn_disable_omi_features) then
      common_range_ok(0:nTimesRadRR-1) = .true.
    else
    ! --------------------------------------------------
    ! Compute the common mode for the Radiance Reference
    ! granule
    ! --------------------------------------------------
    common_range_ok(0:nTimesRadRR-1) = .FALSE.
    CALL find_swathline_range ( &
      TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)),  &
      nTimesRadRR, nXtrackRadRR, rt%latitude(1:nXtrackRadRR,0:nTimesRadRR-1), &
      common_latrange(1:2), common_range_ok(0:nTimesRadRR-1), locerrstat        )
    endif

    ! ----------------------------------------------------------
    ! Interface to the loop over all swath lines for common mode
    ! ----------------------------------------------------------
    call tell_log (1, "Reference_Sector_radiance_reference_granule_retrieval:  "// &
                   "calling swathline_loops (common mode)")
    CALL swathline_loops (                                    &
      pge_idx, rpt_rr, n_max_rspec, &
      common_range_ok(0:nTimesRadRR-1),                           &
      omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2),                   &
      .FALSE., -1,                         &
      .TRUE., locerrstat)

    ! ---------------------------------------------------
    ! Set the index value of the Common Mode spectrum and
    ! assign values to the fitting parameter arrays
    ! ---------------------------------------------------
    CALL finalize_common_mode (nXtrackRadRR)

    ! --------------------------------------
    ! Interface to loop over all swath lines
    ! --------------------------------------
    call tell_log (1, "Reference_Sector_radiance_reference_granule_retrieval:  "// &
                   "calling swathline_loops (rad. reference)")
    CALL swathline_loops (                             &
      pge_idx, rpt_rr, n_max_rspec, &
      radfitref_range_ok(0:nTimesRadRR-1),                        &
      omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2),                   &
      .FALSE., -1,                         &
      .TRUE., locerrstat, retrieval_opt=rt)
      !mem_column_amount(1:nXtrackRadRR,0:nTimesRadRR-1),  &
      !mem_column_uncertainty(1:nXtrackRadRR,0:nTimesRadRR-1),     &
      !mem_rms(1:nXtrackRadRR,0:nTimesRadRR-1),                    &
      !mem_fit_flag(1:nXtrackRadRR,0:nTimesRadRR-1),               &
      !mem_xtrflg(1:nXtrackRadRR,0:nTimesRadRR-1),                 &
      !mem_latitude(1:nXtrackRadRR,0:nTimesRadRR-1),               &
      !mem_longitude(1:nXtrackRadRR,0:nTimesRadRR-1),              &
      !mem_sza(1:nXtrackRadRR,0:nTimesRadRR-1),                    &
      !mem_vza(1:nXtrackRadRR,0:nTimesRadRR-1),                    &
      !mem_height(1:nXtrackRadRR,0:nTimesRadRR-1),                 &
      !locerrstat)

    ! ------------------------------------------------------
    ! mem_column_uncertainty, men_latitude, mem_fit_flag and
    ! mem_column_uncertainty are holding the radiance refere
    ! nce in memory (SCD).
    ! Now I need to work out VCD.
    ! ------------------------------------------------------
    ! Read snow and ice data
    ! ----------------------
    ! ----------------------------------
    ! Read L1b glint and snow/ice flags
    ! ----------------------------------
    CALL omi_read_glint_ice_flags ( &
      l1b_radref_filename, nXtrackRadRR, nTimesRadRR, mem_snow, &
      mem_glint, locerrstat )

    ! --------------------
    ! OMI zoom mode or not
    ! --------------------
    CALL omi_read_binning_factor ( &
      TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
      nTimesRadRR, binfac_rs(0:nTimesRadRR-1), yn_szoom_rs(0:nTimesRadRR-1),     &
      locerrstat )

    ! ------------------------------------------------------
    ! Compute AMF using internal variables to avoud conflict
    ! with main retrieval.
    ! ---------------------------------------------------------
    ! No output of the amf calculation for the reference sector
    ! Height equal to 0 everyehwere. Pacific Ocean
    ! ---------------------------------------------------------
    yn_write = .FALSE.
    ! --------------------------------------------------------
    ! To read clouds corresponding with the radiance reference
    ! granule
    ! --------------------------------------------------------
    voc_amf_filenames(voc_omicld_idx) = TRIM(ADJUSTL(OMSAO_refseccor_cld_filename))

    CALL amf_calculation_bis (                                            &
      pge_idx, nTimesRadRR, nXtrackRadRR, rt%latitude, rt%longitude, &
      rt%sza, rt%vza, rt%saa, rt%vaa, rt%time, mem_snow, mem_glint, &
      omi_xtrpix_range_rr, yn_szoom_rs, rt%column_amount, &
      rt%column_uncertainty, mem_amf, rt%height, yn_write, locerrstat )

    ! --------------------------------------------------------
    ! Compute average fitting statistics and main quality flag
    ! --------------------------------------------------------
    CALL compute_fitting_statistics ( &
      nTimesRadRR, nXtrackRadRR, omi_xtrpix_range_rr,     &
      rt%column_amount, rt%column_uncertainty, rt%rms,          &
      rt%fit_flag, ref_stats, locerrstat )

    ! ------------------------------------------------
    ! Apply the background correction to Slant columns
    ! ------------------------------------------------
    rt%column_amount = rt%column_amount * mem_amf

    ! -----------------------------------------------------
    ! Once we have the radiance reference retrievals we can
    ! compute the background correction
    ! -----------------------------------------------------
    !CALL compute_background_median(nTimesRadRR, nXtrackRadRR,     &
    !     mem_column_amount, &! mem_column_uncertainty, mem_amf,
    !     mem_latitude, mem_longitude, mem_xtrflg, refmqf, locerrstat)

    CALL compute_background_correction_bis(rt%column_amount, mem_correction, &
      rt%latitude, mem_amf,                 &
      nXtrackRadRR, nTimesRadRR, ref_stats % quality_flag, locerrstat)

    call dealloc_retrieval_type (rt, errstat)
    if (errstat /= 0) return

    ! -------------------------------
    ! No more mess with the filenames
    ! -------------------------------
    l1b_rad_filename = l1b_rad_save_filename

    ! -----------
    ! Error check
    ! -----------
    locerrstat = MAX ( locerrstat, errstat )
    IF ( locerrstat >= pge_errstat_error )  RETURN

  END SUBROUTINE Reference_Sector_radref_retrieval_and_median

  SUBROUTINE compute_background_correction_bis(column_amount, correction, latitude, amf, &
      nXtrackRadRR, nTimesRadRR, refmqf, locerrstat)

    USE OMSAO_precision_module, ONLY: r4, i4
    USE ezspline_interpolation, ONLY: ezspline_1d_interpolation
    USE OMSAO_he5_module, ONLY: granule_month
    USE OMSAO_errstat_module, only : pge_errstat_ok
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                           INTENT(IN) :: &
      nTimesRadRR, nXtrackRadRR
    REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(IN) :: &
      column_amount, amf
    INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(IN) :: refmqf
    REAL    (KIND=r4), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(IN) :: latitude

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT(INOUT) :: locerrstat
    REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(INOUT) :: &
      correction

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: iline, itrack
    REAL    (KIND=r8), DIMENSION(1)        :: Ref_column, latitude_r8
    REAL    (KIND=r8), DIMENSION(maxngrid) :: Ref_column_month

    ! -------------------
    ! Routine starts here
    ! -------------------
    locerrstat       = pge_errstat_ok

    Ref_column_month = Reference_sector_concentration(1:maxngrid,granule_month)
    ! -----------------------------------------------------------
    ! Loop pixel by pixel to compute the background correction as
    ! sao_colum - reference_colum
    ! -----------------------------------------------------------
    DO iline = 0,nTimesRadRR-1

      DO itrack = 1,nXtrackRadRR

        ! -----------------------------------------------------
        ! If the main quality flag of the retrieval is not good
        ! then cycle.
        ! -----------------------------------------------------
        IF (refmqf(itrack,iline) .NE. 0) CYCLE
        latitude_r8(1) = REAL(latitude(itrack,iline), KIND=r8)
        ! ---------------------------------------------------
        ! For the rest of the pixels find out the latitude,
        ! interpolate the reference column and substrack from
        ! the retrieved colum
        ! ---------------------------------------------------
        CALL ezspline_1d_interpolation ( INT(ngridpoints, KIND=i4),     &
          grid_lat(1:ngridpoints), Ref_column_month(1:ngridpoints), &
          1, latitude_r8(1:1), Ref_column(1), locerrstat )
        correction(itrack,iline) = column_amount(itrack,iline) - ( Ref_column(1) * amf(itrack,iline) )

      END DO
    END DO

  END SUBROUTINE compute_background_correction_bis

!UNUSED!   SUBROUTINE compute_background_median(nTimesRadRR, nXtrackRadRR, &
!UNUSED!       mem_column_amount, & !mem_column_uncertainty, mem_amf,
!UNUSED!       mem_latitude, mem_longitude, mem_xtrflg, refmqf, locerrstat)
!UNUSED!
!UNUSED!     USE OMSAO_median_module, ONLY: median
!UNUSED!
!UNUSED!     IMPLICIT NONE
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Input variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4),                                           &
!UNUSED!       INTENT(IN) :: nTimesRadRR, nXtrackRadRR
!UNUSED!     REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
!UNUSED!       INTENT(IN) :: mem_column_amount !, mem_column_uncertainty, mem_amf
!UNUSED!     REAL    (KIND=r4), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
!UNUSED!       INTENT(IN) :: mem_latitude, mem_longitude
!UNUSED!     INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
!UNUSED!       INTENT(IN) :: refmqf, mem_xtrflg
!UNUSED!
!UNUSED!     ! ------------------
!UNUSED!     ! Modified variables
!UNUSED!     ! ------------------
!UNUSED!     INTEGER (KIND=i4), INTENT(INOUT) :: locerrstat
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Local variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4)                           :: igrid, itrack, iline
!UNUSED!     REAL    (KIND=r8), DIMENSION(ngridpoints+1) :: grid
!UNUSED!     REAL    (KIND=r8), DIMENSION(1000)          :: into_median
!UNUSED!     REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) &
!UNUSED!       :: column_amounts
!UNUSED!     INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) &
!UNUSED!       :: yn_median
!UNUSED!     INTEGER (KIND=i4)                           :: npixels, ipixel
!UNUSED!
!UNUSED!     ! -------------------
!UNUSED!     ! Routine starts here
!UNUSED!     ! -------------------
!UNUSED!
!UNUSED!     locerrstat = pge_errstat_ok
!UNUSED!
!UNUSED!     background_level = r8_missval
!UNUSED!
!UNUSED!     ! ----------------------------------------------------------------
!UNUSED!     ! Computing the median of all the pixels of the radiance reference
!UNUSED!     ! granule within two consequtive latitudes of the background conce
!UNUSED!     ! trations grid readed from the background_concentrations file.
!UNUSED!     ! Only pixels between 180W and 140W are considered. Indeed because
!UNUSED!     ! of the conditions impossed on the radiance reference granule, as
!UNUSED!     ! close as possible to the 165W, they should be no pixels above
!UNUSED!     ! 180W and below 140W.
!UNUSED!     ! ----------------------------------------------------------------
!UNUSED!     ! Generating grid
!UNUSED!     ! ---------------
!UNUSED!     DO igrid = 1, ngridpoints+1
!UNUSED!       grid(igrid) = (180.0_r8 / (REAL(ngridpoints+1, KIND=r8)) * &
!UNUSED!         (REAL(igrid, KIND=r8))) - 90.0_r8
!UNUSED!     END DO
!UNUSED!
!UNUSED!     DO igrid = 1, ngridpoints
!UNUSED!
!UNUSED!       column_amounts        = 0.0_r8
!UNUSED!       yn_median             = 0_i2
!UNUSED!       into_median           = 0.0_r8
!UNUSED!
!UNUSED!       ! -------------------------------------------------------------------
!UNUSED!       ! Finging pixels in the granule between grid(igrid) and grid(igrid+1)
!UNUSED!       ! Mid point is grid_lat(igrid).
!UNUSED!       ! -------------------------------------------------------------------
!UNUSED!       WHERE (mem_latitude .gt. grid(igrid) .AND. mem_latitude .lt. grid(igrid+1) &
!UNUSED!           .AND. refmqf .EQ. 0 .AND. mem_xtrflg .EQ. 0 .AND. mem_longitude     &
!UNUSED!           .GT. -180 .AND. mem_longitude .LT. -140)
!UNUSED!         column_amounts = mem_column_amount
!UNUSED!         yn_median      = 1
!UNUSED!       END WHERE
!UNUSED!       npixels = 0_i4
!UNUSED!       npixels = SUM(yn_median)
!UNUSED!       ipixel  = 1_i4
!UNUSED!       DO itrack = 1, nXtrackRadRR
!UNUSED!         DO iline = 0, nTimesRadRR-1
!UNUSED!           IF (yn_median(itrack,iline) .EQ. 1) THEN
!UNUSED!             into_median(ipixel) = column_amounts(itrack,iline)
!UNUSED!             ipixel = ipixel + 1
!UNUSED!           END IF
!UNUSED!         END DO
!UNUSED!       END DO
!UNUSED!       IF (npixels .GT. 0) THEN
!UNUSED!         background_level(igrid) = median(npixels, into_median(1:npixels))
!UNUSED!       END IF
!UNUSED!
!UNUSED!     END DO
!UNUSED!
!UNUSED!   END SUBROUTINE COMPUTE_BACKGROUND_median

!UNUSED!   SUBROUTINE compute_background_correction(locerrstat)
!UNUSED!
!UNUSED!     IMPLICIT NONE
!UNUSED!
!UNUSED!     INTEGER (KIND=i4), INTENT (INOUT) :: locerrstat
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Local variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i2) :: igrid
!UNUSED!
!UNUSED!     locerrstat = pge_errstat_ok
!UNUSED!
!UNUSED!     background_correction = 0.0_r8
!UNUSED!
!UNUSED!     ! ----------------------------------------
!UNUSED!     ! Month we are dealing with: granule_month
!UNUSED!     ! Working out the difference between the r
!UNUSED!     ! eference sector median values and the GE
!UNUSED!     ! OS-Chem values
!UNUSED!     ! ----------------------------------------
!UNUSED!     DO igrid = 1, ngridpoints
!UNUSED!       IF (background_level(igrid) .NE. r8_missval) THEN
!UNUSED!         background_correction(igrid) = background_level(igrid) - &
!UNUSED!           Reference_sector_concentration(igrid,granule_month)
!UNUSED!
!UNUSED!       END IF
!UNUSED!     END DO
!UNUSED!
!UNUSED!   END SUBROUTINE compute_background_correction

!UNUSED!   SUBROUTINE aply_background_correction(ntimes, nxtrack, saocol, &! saodco,
!UNUSED!                                         saomqf, latitude, locerrstat)
!UNUSED!
!UNUSED!     USE ezspline_interpolation, ONLY: ezspline_1d_interpolation
!UNUSED!     IMPLICIT NONE
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Input variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
!UNUSED!     REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN) :: latitude
!UNUSED!     INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN) :: saomqf
!UNUSED!
!UNUSED!     ! ------------------
!UNUSED!     ! Modified variables
!UNUSED!     ! ------------------
!UNUSED!     REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (INOUT) :: saocol
!UNUSED!     !REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (INOUT) :: saodco
!UNUSED!     INTEGER (KIND=i4),                                   INTENT (INOUT) :: locerrstat
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Local variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4) :: itimes, itrack
!UNUSED!     REAL    (KIND=r8), DIMENSION(1) :: correction, out_lat
!UNUSED!
!UNUSED!     locerrstat = pge_errstat_ok
!UNUSED!     !npixels = 1
!UNUSED!
!UNUSED!     DO itrack = 1, nxtrack
!UNUSED!       DO itimes = 0, ntimes-1
!UNUSED!
!UNUSED!         correction(1) = 0.0_r8
!UNUSED!
!UNUSED!         ! ----------------------------------------------
!UNUSED!         ! Only to be applied if the main quality flag is
!UNUSED!         ! OK: .EQ. 0
!UNUSED!         ! ----------------------------------------------
!UNUSED!         IF (saomqf(itrack,itimes) .EQ. 0) THEN
!UNUSED!
!UNUSED!           ! ----------------------------------------------
!UNUSED!           ! For each pixel, take the background_correction
!UNUSED!           ! and interpolate it to latitude(itrack, itimes)
!UNUSED!           ! ----------------------------------------------
!UNUSED!           out_lat(1) = REAL(latitude(itrack,itimes), KIND=r8)
!UNUSED!           CALL ezspline_1d_interpolation ( INT(ngridpoints, KIND=i4),  &
!UNUSED!             grid_lat, background_correction,                       &
!UNUSED!             1, out_lat(1), correction(1), &
!UNUSED!             locerrstat )
!UNUSED!
!UNUSED!           ! ---------------------------------------------------------
!UNUSED!           ! And apply the result of the interpolation to the original
!UNUSED!           ! saocol
!UNUSED!           ! ---------------------------------------------------------
!UNUSED!           saocol(itrack,itimes) = saocol(itrack,itimes) - correction(1)
!UNUSED!
!UNUSED!         END IF
!UNUSED!
!UNUSED!       END DO
!UNUSED!     END DO
!UNUSED!
!UNUSED!   END SUBROUTINE aply_background_correction

  SUBROUTINE he5_write_reference_sector_corrected_column &
      (pge_idx, &
       nt, nx, column, & !uncertainty,
       errstat)

    USE OMSAO_he5_module, ONLY: pge_swath_id, he5_swwrfld, &
      he5_start_2d, he5_stride_2d, he5_edge_2d
    USE OMSAO_precision_module, ONLY: i4
    !USE OMSAO_omidata_module,   ONLY: n_roff_dig
    USE OMSAO_indices_module,   ONLY: pge_hcho_idx
    !USE sao_pge_utils, ONLY: roundoff_2darr_r8
    use datafields, only: rscol_field
    USE OMSAO_errstat_module, only : pge_errstat_ok, pge_errstat_error, &
      he5_stat_ok, omsao_e_he5swwrfld, vb_lev_default, error_check

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nt, nx
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: column
    !REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: uncertainty

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=43), PARAMETER :: modulename = &
      'he5_write_reference_sector_corrected_column'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                          :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

    ! -----------------------------------------------
    ! Total column corrected by the reference sector.
    ! -----------------------------------------------
    ! Only implemented for HCHO.
    ! -----------------------------------------------

    locerrstat = pge_errstat_ok

    he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)

    ! ----------------------------------------------------------------------
    ! All PGEs: Output of columns and column uncertainties left as a comment
    ! for the future. So far only for HCHO
    ! ----------------------------------------------------------------------
    IF (pge_idx .EQ. pge_hcho_idx) THEN
      colloc = column
      !CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
      locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(rscol_field)), &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
      errstat = MAX ( errstat, locerrstat )

      !colloc = uncertainty
      !CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
      !locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(rsdcol_field)), &
      !     he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
      !errstat = MAX ( errstat, locerrstat )

    END IF

    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
      modulename, vb_lev_default, errstat )

    RETURN

  END SUBROUTINE he5_write_reference_sector_corrected_column

END MODULE OMSAO_Reference_sector_module
