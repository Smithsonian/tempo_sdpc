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
    use ctrlvars, only : yn_gems
    use m_read_gems, only : gems_read_l1_rad_info
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
    if (.not. yn_gems) then !TEMPO
      call read_l1_radiance_info (l1b_radref_filename, l1b_channel, &
                                rpt_rr, errstat)
    else ! GEMS
      call gems_read_l1_rad_info (l1b_radref_filename, rpt_rr, errstat)
    endif
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

  END SUBROUTINE Read_reference_sector_concentration

  SUBROUTINE Reference_Sector_radref_retrieval_and_median &
      (pge_idx, rpt_rr, n_max_rspec, errstat, mem_correction)

    USE OMSAO_precision_module, ONLY: i1, i4!, r4
    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, i2_missval, r8_missval, & ! , r4_missval
      MAX_STR_LEN
    USE OMSAO_wfamf_module,     ONLY: amf_calculation
    USE OMSAO_variables_module, ONLY: OMSAO_refseccor_cld_filename, voc_amf_filenames, &
      Radiance_Paras_Type, common_latrange, l1b_rad_filename, l1b_radref_filename
    USE OMSAO_indices_module,   ONLY: voc_omicld_idx
    USE omi_pge_fitting_aux, ONLY: read_latitude, find_swathline_range, &
      compute_fitting_statistics, fitting_statistics_type
    use commonmode, only: finalize_common_mode
    USE omi_read_l1b_data, ONLY: omi_read_glint_ice_land_flags, omi_read_binning_factor
    USE swathline_loop, ONLY: swathline_loops
    !USE omi_pge_swathline_loop_memory, ONLY: omi_pge_swathline_loops_mem
    USE OMSAO_errstat_module, only : pge_errstat_ok, pge_errstat_error
    USE OMSAO_omidata_module, ONLY: omi_radiance_swathname, &
      retrieval_type, alloc_retrieval_type, dealloc_retrieval_type
    use ctrlvars, only: yn_disable_omi_features, yn_gems
    use m_read_gems, only: gems_read_latitude, gems_read_ice_glint

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
    INTEGER (KIND=i2), DIMENSION (rpt_rr%nxtrack,0:rpt_rr%ntimes-1) :: mem_snow, mem_glint, mem_land_water, mem_amfdiag
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

    ! ---------------------------------------------------
    ! mem_correction needs to be initialized to 0.0 since
    ! we want no correction if any of the two reference
    ! sector retrievals does fail.
    ! ---------------------------------------------------
    mem_correction         = 0.0
    ! ---------------------------------------------------
    mem_amf                = r8_missval
    ref_stats % quality_flag  = i2_missval
    mem_snow               = i2_missval
    mem_glint              = i2_missval
    mem_land_water         = i2_missval

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
    if (.not. yn_gems) then !TEMPO
      CALL read_latitude ( TRIM(ADJUSTL(l1b_rad_filename)), &
           TRIM(ADJUSTL(omi_radiance_swathname)), &
           0, nTimesRadRR, rt%latitude(1:nXtrackRadRR,0:nTimesRadRR-1), &
           errstat)
    else
      call gems_read_latitude (trim(adjustl(l1b_rad_filename)), 0, &
           nTimesRadRR, rt%latitude(1:nXtrackRadRR,0:nTimesRadRR-1), errstat)
    endif
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
    if (.not. yn_gems) then !TEMPO
      CALL omi_read_glint_ice_land_flags ( &
           l1b_radref_filename, nXtrackRadRR, nTimesRadRR, mem_snow, &
           mem_glint, mem_land_water, locerrstat )
    else !GEMS
      call gems_read_ice_glint(l1b_radref_filename, nXtrackRadRR, nTimesRadRR,&
           mem_snow, mem_glint, locerrstat)
    endif

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

    CALL amf_calculation (                                            &
      pge_idx, nTimesRadRR, nXtrackRadRR, rt%latitude, rt%longitude, &
      rt%sza, rt%vza, rt%saa, rt%vaa, rt%time, mem_snow, mem_glint, &
      mem_land_water, omi_xtrpix_range_rr, rt%column_amount, &
      rt%column_uncertainty, mem_amf, mem_amfdiag, rt%height, yn_write, locerrstat )

    ! --------------------------------------------------------
    ! Compute average fitting statistics and main quality flag
    ! --------------------------------------------------------
    CALL compute_fitting_statistics ( &
      nTimesRadRR, nXtrackRadRR, omi_xtrpix_range_rr,     &
      rt%column_amount, rt%column_uncertainty, rt%rms, rt%fit_flag, mem_amf, &
      mem_amfdiag, rt%sza, rt%vza, ref_stats, locerrstat )

    ! ------------------------------------------------
    ! Apply the background correction to Slant columns
    ! ------------------------------------------------
    rt%column_amount = rt%column_amount * mem_amf

    ! -----------------------------------------------------
    ! Once we have the radiance reference retrievals we can
    ! compute the background correction
    ! -----------------------------------------------------

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

END MODULE OMSAO_Reference_sector_module
