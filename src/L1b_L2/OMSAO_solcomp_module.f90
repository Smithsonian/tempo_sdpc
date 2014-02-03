MODULE OMSAO_solcomp_module

  ! =================================================================
  !
  ! This module defines variables associated with parameterization of
  ! the composite solar spectrum.
  !
  ! =================================================================

  USE OMSAO_precision_module,  ONLY: i4, r8, C_LONG
  USE OMSAO_parameters_module,   ONLY: MAX_STR_LEN

  IMPLICIT NONE

  PRIVATE i4, r8, C_LONG, MAX_STR_LEN
  ! ----------------------------------------------------------------------
  ! Type declaration for variable with composite solar spectrum parameters
  ! ----------------------------------------------------------------------
  TYPE, PUBLIC :: SolarCompositeParameters
    CHARACTER (LEN=MAX_STR_LEN)                            :: Title
    INTEGER   (KIND=I4)                                 :: nBins, nPol, nXtr
    REAL      (KIND=r8)                                 :: SolarNorm
    REAL      (KIND=r8), DIMENSION (2)                  :: FirstLastWav
    REAL      (KIND=r8), DIMENSION (:),     ALLOCATABLE :: SolComParWavs
    REAL      (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: SolComParData
  END TYPE SolarCompositeParameters

  ! -----------------------------------------------------------
  ! Variable that holds the tabulated solar spectrum parameters
  ! -----------------------------------------------------------
  TYPE (SolarCompositeParameters) :: solarcomp_pars

  ! ---------------------------------------------------
  ! The division orbit between OPF versions < and >= 25
  ! ---------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: opf_ge25 = 25

  ! ----------------------------------------------------------------
  ! Strings to aid in the construction of the HE5 field name to read
  ! ----------------------------------------------------------------
  CHARACTER (LEN= 4), DIMENSION (3), PARAMETER :: ch_typ = (/ 'MEA_', 'MED_', 'PC0_' /)
  CHARACTER (LEN= 8), DIMENSION (2), PARAMETER :: ch_opf = (/ 'OPFlt25_', 'OPFge25_' /)
  CHARACTER (LEN= 5), DIMENSION (2), PARAMETER :: ch_orb = (/ '1st7_',    'comp_'    /)
  CHARACTER (LEN= 3), DIMENSION (3), PARAMETER :: ch_cha = (/ 'UV1',  'UV2',   'VIS' /)

  ! ---------------------------
  ! 32bit/64bit C_LONG integers
  ! ---------------------------
  INTEGER (KIND=C_LONG), PARAMETER, PRIVATE :: zerocl = 0, onecl = 1

  PRIVATE soco_get_dims

CONTAINS

  SUBROUTINE soco_pars_read ( &
      soco_fname, typ_idx, l1b_channel, winwav_min, winwav_max, errstat )

    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN
    USE OMSAO_he5_module, ONLY: HE5_SWopen, HE5_SWattach, HE5_SWrdfld, &
      HE5_SWrdlattr, HE5_SWdetach, HE5_SWclose, HE5_SWinqswath, &
      he5f_acc_rdonly
    USE sao_pge_utils, ONLY: array_locate_r8
    USE OMSAO_errstat_module
    IMPLICIT NONE

    !INCLUDE 'hdfeos5.inc'

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    soco_fname ........... file name with Solar Composite spectra
    !    typ_idx .............. whether "MEA", "MED" or "PC0" is requested
    !    l1b_channel .......... OMI channel ("UV1", "UV2", or "VIS")
    !    winwav_min ........... starting wavelength to read from file
    !    winwav_max ........... final wavelength to read from file
    !    errstat .............. error status returned from the subroutine
    ! ---------------------------------------------------------------------

    CHARACTER (LEN=14), PARAMETER :: modulename = 'soco_pars_read'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: typ_idx
    REAL    (KIND=r8), INTENT (IN) :: winwav_min, winwav_max
    CHARACTER (LEN=*), INTENT (IN) :: soco_fname
    CHARACTER (LEN=3), INTENT (IN) :: l1b_channel

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! --------------
    ! Local variable
    ! --------------
    INTEGER (KIND=i4) :: n_xtr, n_wvl, n_pars
    INTEGER (KIND=i4) :: ios, ntmp, j1, j2, swlen
    REAL    (KIND=r8) :: soco_norm
    REAL    (KIND=r8), DIMENSION (:),     ALLOCATABLE :: tmp_wvls
    REAL    (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: tmp_pars

    INTEGER   (KIND=i4)                :: he5stat, swath_file_id, swath_id
    INTEGER   (KIND=C_LONG)               :: he5statcl, he5_start_1, he5_stride_1, he5_edge_1
    INTEGER   (KIND=C_LONG), DIMENSION(3) :: he5_start_3, he5_stride_3, he5_edge_3
    CHARACTER (LEN=MAX_STR_LEN)           :: soco_he5_field, swath_name

    ! -----------------------
    ! Initialize error status
    ! -----------------------
    errstat = pge_errstat_ok

    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(soco_fname)), he5f_acc_rdonly )
    IF ( swath_file_id == he5_stat_fail ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, he5stat )
      errstat = MAX ( errstat, he5stat )
      RETURN
    END IF
    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    swath_name = "" !JED fix
    he5statcl  = HE5_SWinqswath  ( TRIM(ADJUSTL(soco_fname)), swath_name, swlen )
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == he5_stat_fail ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )
      errstat = MAX ( errstat, he5stat )
      RETURN
    END IF

    ! -----------------------------
    ! Read dimensions from HE5 file
    ! -----------------------------
    CALL soco_get_dims ( swath_id, l1b_channel, n_xtr, n_wvl, n_pars )

    ! ------------------------------------
    ! Allocate memory for temporary arrays
    ! ------------------------------------
    ALLOCATE ( tmp_wvls(1:n_wvl), STAT=ios )
    if (ios == 0) ALLOCATE ( tmp_pars(1:n_xtr,1:n_wvl,1:n_pars), STAT=ios )
    if (ios /= 0) then
      call error_check (0, 1, pge_errstat_fatal, OMSAO_F_MEM_ALLOC, &
                        modulename, vb_lev_default, errstat)
      return
    endif

    ! --------------------
    ! Read wavelength bins
    ! --------------------
    he5_start_1 = zerocl ; he5_stride_1 = onecl ; he5_edge_1 = INT(n_wvl,KIND=C_LONG)
    he5stat = HE5_SWrdfld ( swath_id, "SoCo_Wvl_"//l1b_channel,   &
      he5_start_1, he5_stride_1, he5_edge_1, tmp_wvls(1:n_wvl) )
    IF ( he5stat /= pge_errstat_ok ) &
      CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
      modulename//f_sep//'Solar Composite access failed.', vb_lev_default, errstat )

    ! ----------------------
    ! Compose HE5 field name
    ! ----------------------
    soco_he5_field = ' '
    soco_he5_field = ch_typ(typ_idx)//l1b_channel

    ! -------------------------------
    ! Read Solar Composite parameters
    ! -------------------------------
    he5_start_3  = (/               zerocl,               zerocl,                zerocl /)
    he5_stride_3 = (/               onecl,               onecl,                onecl /)
    he5_edge_3   = (/ INT(n_xtr,KIND=C_LONG), INT(n_wvl,KIND=C_LONG), INT(n_pars,KIND=C_LONG) /)
    he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(soco_he5_field)),   &
      he5_start_3, he5_stride_3, he5_edge_3, tmp_pars(1:n_xtr,1:n_wvl,1:n_pars) )
    IF ( he5stat /= pge_errstat_ok ) &
      CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
      modulename//f_sep//'Solar Composite access failed.', vb_lev_default, errstat )

    ! -------------------------
    ! Read Solar Norm attribute
    ! -------------------------
    he5stat = HE5_SWrdlattr ( &
      swath_id, TRIM(ADJUSTL(soco_he5_field)), "SolarSpectrumNorm", soco_norm )

    ! -----------------------------------------------
    ! Detach from HE5 swath and close HE5 output file
    ! -----------------------------------------------
    he5stat = HE5_SWdetach ( swath_id )
    he5stat = HE5_SWclose  ( swath_file_id )

    ! --------------------------------------------
    ! Find the slice to store in the global arrays
    ! --------------------------------------------
    CALL array_locate_r8 ( n_wvl, tmp_wvls(1:n_wvl), winwav_min, 'LE', j1 )
    CALL array_locate_r8 ( n_wvl, tmp_wvls(1:n_wvl), winwav_max, 'GE', j2 )

    j1 = MAX ( j1, 1 )
    IF ( j2 < 1 ) j2 = n_wvl

    ntmp = j2 - j1 + 1

    ! ----------------
    ! Assign the title
    ! ----------------
    solarcomp_pars%Title = 'OMI Solar Composite Spectrum --- '//TRIM(ADJUSTL(soco_he5_field))

    ! --------------------------------
    ! Assign variables to global array
    ! --------------------------------
    solarcomp_pars%nXtr  = ntmp
    solarcomp_pars%nBins = ntmp
    solarcomp_pars%nPol  = n_pars

    ALLOCATE ( solarcomp_pars%SolComParWavs(1:ntmp), STAT=ios )
    if (ios == 0) ALLOCATE ( solarcomp_pars%SolComParData(1:n_xtr,1:ntmp,1:n_pars), STAT=ios )
    if (ios /= 0) then
      call error_check (0, 1, pge_errstat_fatal, OMSAO_F_MEM_ALLOC, &
                        modulename, vb_lev_default, errstat)
      return
    endif

    solarcomp_pars%SolComParWavs(1:ntmp)                  = tmp_wvls(j1:j2)
    solarcomp_pars%SolComParData(1:n_xtr,1:ntmp,1:n_pars) = tmp_pars(1:n_xtr,j1:j2,1:n_pars)

    solarcomp_pars%SolarNorm = soco_norm

    ! ---------------------------------------
    ! De-allocate memory for temporary arrays
    ! ---------------------------------------
    IF ( ALLOCATED(tmp_wvls) ) DEALLOCATE (tmp_wvls)
    IF ( ALLOCATED(tmp_pars) ) DEALLOCATE (tmp_pars)

    RETURN
  END SUBROUTINE soco_pars_read

!UNUSED!   SUBROUTINE soco_pars_read_v002 ( &
!UNUSED!       soco_fname, typ_idx, sel_idx, l1b_channel, l1br_opf_version, winwav_min, winwav_max, errstat )
!UNUSED! 
!UNUSED!     ! -----------------------------------------------------------
!UNUSED!     ! Old version (Collection 2 L1b data) of SOCO_PARS_READ above
!UNUSED!     ! -----------------------------------------------------------
!UNUSED! 
!UNUSED!     USE OMSAO_parameters_module, ONLY: MAX_STR_LEN
!UNUSED!     USE OMSAO_he5_module, ONLY: HE5_SWopen, HE5_SWattach, HE5_SWrdfld, &
!UNUSED!       HE5_SWrdlattr, HE5_SWdetach, HE5_SWclose, HE5_SWinqswath, &
!UNUSED!       he5f_acc_rdonly
!UNUSED!     USE sao_pge_utils, ONLY: array_locate_r8
!UNUSED! 
!UNUSED!     IMPLICIT NONE
!UNUSED! 
!UNUSED!     !INCLUDE 'hdfeos5.inc'
!UNUSED! 
!UNUSED!     ! ---------------------------------------------------------------------
!UNUSED!     ! Explanation of subroutine arguments:
!UNUSED!     !
!UNUSED!     !    soco_fname ........... file name with Solar Composite spectra
!UNUSED!     !    typ_idx .............. whether "MED" or "PC0" is requested
!UNUSED!     !    sel_idx .............. whether "1st7" or "comp" is requested
!UNUSED!     !    l1b_channel .......... OMI channel ("UV1", "UV2", or "VIS")
!UNUSED!     !    l1br_opf_version ..... OPF version of L1B radiance file
!UNUSED!     !    winwav_min ........... starting wavelength to read from file
!UNUSED!     !    winwav_max ........... final wavelength to read from file
!UNUSED!     !    errstat .............. error status returned from the subroutine
!UNUSED!     ! ---------------------------------------------------------------------
!UNUSED! 
!UNUSED!     CHARACTER (LEN=19), PARAMETER :: modulename = 'soco_pars_read_v002'
!UNUSED! 
!UNUSED!     ! ---------------
!UNUSED!     ! Input variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4), INTENT (IN) :: typ_idx, sel_idx, l1br_opf_version
!UNUSED!     REAL    (KIND=r8), INTENT (IN) :: winwav_min, winwav_max
!UNUSED!     CHARACTER (LEN=*), INTENT (IN) :: soco_fname
!UNUSED!     CHARACTER (LEN=3), INTENT (IN) :: l1b_channel
!UNUSED! 
!UNUSED!     ! -----------------
!UNUSED!     ! Modified variable
!UNUSED!     ! -----------------
!UNUSED!     INTEGER (KIND=i4), INTENT (INOUT) :: errstat
!UNUSED! 
!UNUSED!     ! --------------
!UNUSED!     ! Local variable
!UNUSED!     ! --------------
!UNUSED!     INTEGER (KIND=i4) :: n_xtr, n_wvl, n_pars
!UNUSED!     INTEGER (KIND=i4) :: ios, ntmp, j1, j2, opf_idx, swlen
!UNUSED!     REAL    (KIND=r8) :: soco_norm
!UNUSED!     REAL    (KIND=r8), DIMENSION (:),     ALLOCATABLE :: tmp_wvls
!UNUSED!     REAL    (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: tmp_pars
!UNUSED! 
!UNUSED!     INTEGER   (KIND=i4)                :: he5stat, swath_file_id, swath_id
!UNUSED!     INTEGER   (KIND=C_LONG)            :: he5stat_long
!UNUSED!     INTEGER   (KIND=C_LONG)               :: he5_start_1, he5_stride_1, he5_edge_1
!UNUSED!     INTEGER   (KIND=C_LONG), DIMENSION(3) :: he5_start_3, he5_stride_3, he5_edge_3
!UNUSED!     CHARACTER (LEN=MAX_STR_LEN)           :: soco_he5_field, swath_name
!UNUSED! 
!UNUSED!     ! ---------------------------------
!UNUSED!     ! External OMI and Toolkit routines
!UNUSED!     ! ---------------------------------
!UNUSED!     !INTEGER (KIND=i4),  EXTERNAL :: &
!UNUSED!     !  HE5_SWopen, HE5_SWattach, HE5_SWrdfld, HE5_SWrdlattr, HE5_SWdetach, HE5_SWclose
!UNUSED!     !INTEGER (KIND=C_LONG), EXTERNAL :: &
!UNUSED!     !  HE5_SWinqswath
!UNUSED! 
!UNUSED!     ! -----------------------
!UNUSED!     ! Initialize error status
!UNUSED!     ! -----------------------
!UNUSED!     errstat = pge_errstat_ok
!UNUSED! 
!UNUSED!     ! --------------------------------------------------------------------------
!UNUSED!     ! Select either pre- or post-OPFv25 data. The OPF version has been read from
!UNUSED!     ! the L1b radiance file in INIT_METADATA.
!UNUSED!     ! --------------------------------------------------------------------------
!UNUSED!     IF ( l1br_opf_version >= opf_ge25 ) THEN
!UNUSED!       opf_idx = 2
!UNUSED!     ELSE
!UNUSED!       opf_idx = 1
!UNUSED!     END IF
!UNUSED! 
!UNUSED!     ! -----------------------------------------------------------
!UNUSED!     ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
!UNUSED!     ! -----------------------------------------------------------
!UNUSED!     swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(soco_fname)), he5f_acc_rdonly )
!UNUSED!     IF ( swath_file_id == he5_stat_fail ) THEN
!UNUSED!       CALL error_check ( &
!UNUSED!         0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, he5stat )
!UNUSED!       errstat = MAX ( errstat, he5stat )
!UNUSED!       RETURN
!UNUSED!     END IF
!UNUSED! 
!UNUSED!     ! ---------------------------------------------
!UNUSED!     ! Check for existing HE5 swath and attach to it
!UNUSED!     ! ---------------------------------------------
!UNUSED!     he5stat_long = HE5_SWinqswath  ( TRIM(ADJUSTL(soco_fname)), swath_name, swlen )
!UNUSED!     he5stat = INT(he5stat_long, KIND=i4)
!UNUSED!     swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
!UNUSED!     IF ( swath_id == he5_stat_fail ) THEN
!UNUSED!       CALL error_check ( &
!UNUSED!         0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )
!UNUSED!       errstat = MAX ( errstat, he5stat )
!UNUSED!       RETURN
!UNUSED!     END IF
!UNUSED! 
!UNUSED!     ! -----------------------------
!UNUSED!     ! Read dimensions from HE5 file
!UNUSED!     ! -----------------------------
!UNUSED!     CALL soco_get_dims ( swath_id, l1b_channel, n_xtr, n_wvl, n_pars )
!UNUSED!     ! ------------------------------------
!UNUSED!     ! Allocate memory for temporary arrays
!UNUSED!     ! ------------------------------------
!UNUSED!     ALLOCATE ( tmp_wvls(1:n_wvl),                  STAT=ios )
!UNUSED!     if (ios == 0) ALLOCATE ( tmp_pars(1:n_xtr,1:n_wvl,1:n_pars), STAT=ios )
!UNUSED!     if (ios /= 0) then
!UNUSED!       call error_check (0, 1, pge_errstat_fatal, OMSAO_F_MEM_ALLOC, &
!UNUSED!                         modulename, vb_lev_default, errstat)
!UNUSED!       return
!UNUSED!     endif
!UNUSED! 
!UNUSED!     ! --------------------
!UNUSED!     ! Read wavelength bins
!UNUSED!     ! --------------------
!UNUSED!     he5_start_1 = 0 ; he5_stride_1 = 1 ; he5_edge_1 = n_wvl
!UNUSED!     he5stat = HE5_SWrdfld ( swath_id, "SoCo_Wvl_"//l1b_channel,   &
!UNUSED!       he5_start_1, he5_stride_1, he5_edge_1, tmp_wvls(1:n_wvl) )
!UNUSED!     IF ( he5stat /= pge_errstat_ok ) &
!UNUSED!       CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
!UNUSED!       modulename//f_sep//'Solar Composite access failed.', vb_lev_default, errstat )
!UNUSED! 
!UNUSED!     ! ----------------------
!UNUSED!     ! Compose HE5 field name
!UNUSED!     ! ----------------------
!UNUSED!     soco_he5_field = ' '
!UNUSED!     soco_he5_field = ch_typ(typ_idx)//ch_opf(opf_idx)//ch_orb(sel_idx)//l1b_channel
!UNUSED! 
!UNUSED!     ! -------------------------------
!UNUSED!     ! Read Solar Composite parameters
!UNUSED!     ! -------------------------------
!UNUSED!     he5_start_3  = (/     0,     0,      0 /)
!UNUSED!     he5_stride_3 = (/     1,     1,      1 /)
!UNUSED!     he5_edge_3   = (/ n_xtr, n_wvl, n_pars /)
!UNUSED!     he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(soco_he5_field)),   &
!UNUSED!       he5_start_3, he5_stride_3, he5_edge_3, tmp_pars(1:n_xtr,1:n_wvl,1:n_pars) )
!UNUSED!     IF ( he5stat /= pge_errstat_ok ) &
!UNUSED!       CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
!UNUSED!       modulename//f_sep//'Solar Composite access failed.', vb_lev_default, errstat )
!UNUSED! 
!UNUSED!     ! -------------------------
!UNUSED!     ! Read Solar Norm attribute
!UNUSED!     ! -------------------------
!UNUSED!     he5stat = HE5_SWrdlattr ( &
!UNUSED!       swath_id, TRIM(ADJUSTL(soco_he5_field)), "SolarSpectrumNorm", soco_norm )
!UNUSED! 
!UNUSED!     ! -----------------------------------------------
!UNUSED!     ! Detach from HE5 swath and close HE5 output file
!UNUSED!     ! -----------------------------------------------
!UNUSED!     he5stat = HE5_SWdetach ( swath_id )
!UNUSED!     he5stat = HE5_SWclose  ( swath_file_id )
!UNUSED! 
!UNUSED!     ! --------------------------------------------
!UNUSED!     ! Find the slice to store in the global arrays
!UNUSED!     ! --------------------------------------------
!UNUSED!     CALL array_locate_r8 ( n_wvl, tmp_wvls(1:n_wvl), winwav_min, 'LE', j1 )
!UNUSED!     CALL array_locate_r8 ( n_wvl, tmp_wvls(1:n_wvl), winwav_max, 'GE', j2 )
!UNUSED! 
!UNUSED!     j1 = MAX ( j1, 1 )
!UNUSED!     IF ( j2 < 1 ) j2 = n_wvl
!UNUSED! 
!UNUSED!     ntmp = j2 - j1 + 1
!UNUSED! 
!UNUSED!     ! ----------------
!UNUSED!     ! Assign the title
!UNUSED!     ! ----------------
!UNUSED!     solarcomp_pars%Title = 'OMI Solar Composite Spectrum --- '//TRIM(ADJUSTL(soco_he5_field))
!UNUSED! 
!UNUSED!     ! --------------------------------
!UNUSED!     ! Assign variables to global array
!UNUSED!     ! --------------------------------
!UNUSED!     solarcomp_pars%nXtr  = ntmp
!UNUSED!     solarcomp_pars%nBins = ntmp
!UNUSED!     solarcomp_pars%nPol  = n_pars
!UNUSED! 
!UNUSED!     ALLOCATE ( solarcomp_pars%SolComParWavs(1:ntmp), STAT=ios )
!UNUSED!     if (ios == 0) ALLOCATE ( solarcomp_pars%SolComParData(1:n_xtr,1:ntmp,1:n_pars), STAT=ios )
!UNUSED!     if (ios /= 0) then
!UNUSED!       call error_check (0, 1, pge_errstat_fatal, OMSAO_F_MEM_ALLOC, &
!UNUSED!                         modulename, vb_lev_default, errstat)
!UNUSED!       return
!UNUSED!     endif    
!UNUSED! 
!UNUSED!     solarcomp_pars%SolComParWavs(1:ntmp)                  = tmp_wvls(j1:j2)
!UNUSED!     solarcomp_pars%SolComParData(1:n_xtr,1:ntmp,1:n_pars) = tmp_pars(1:n_xtr,j1:j2,1:n_pars)
!UNUSED! 
!UNUSED!     solarcomp_pars%SolarNorm = soco_norm
!UNUSED! 
!UNUSED!     ! ---------------------------------------
!UNUSED!     ! De-allocate memory for temporary arrays
!UNUSED!     ! ---------------------------------------
!UNUSED!     IF ( ALLOCATED(tmp_wvls) ) DEALLOCATE (tmp_wvls)
!UNUSED!     IF ( ALLOCATED(tmp_pars) ) DEALLOCATE (tmp_pars)
!UNUSED! 
!UNUSED!     RETURN
!UNUSED!   END SUBROUTINE soco_pars_read_v002

  SUBROUTINE soco_get_dims ( swath_id, l1b_channel, n_xtr, n_wvl, n_pars )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swath_id
    CHARACTER (LEN=3), INTENT (IN) :: l1b_channel

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: n_xtr, n_wvl, n_pars

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=C_LONG)                  :: ndimcl
    INTEGER (KIND=i4)                   :: ndim, swlen, istart, iend, nsep, i
    INTEGER (KIND=i4),  DIMENSION (0:9) :: dim_array, dim_seps
    INTEGER (KIND=C_LONG), DIMENSION (0:9) :: dim_arraycl
    CHARACTER (LEN=MAX_STR_LEN)            :: dim_chars

    ! ------------------
    ! EXTERNAL functions
    ! ------------------
    INTEGER (KIND=C_LONG), EXTERNAL :: HE5_SWinqdims

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    dim_chars="" ! --JED
    ndimcl = HE5_SWinqdims  ( swath_id, dim_chars, dim_arraycl(0:6) )
    ndim   = INT ( ndimcl, KIND=i4 )
    dim_array(0:6) = INT ( dim_arraycl(0:6), KIND=i4 )

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    n_xtr = 0 ; n_wvl = 0 ; n_pars = 0

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    dim_chars = TRIM(ADJUSTL(dim_chars))
    swlen = LEN_TRIM(ADJUSTL(dim_chars))
    istart = 1  ;  iend = 1 ; nsep = 0 ; dim_seps = 0

    getseps: DO i = 1, swlen
      IF ( dim_chars(i:i) == ',' ) THEN
        nsep = nsep + 1
        dim_seps(nsep) = i
      END IF
    END DO getseps
    nsep = nsep + 1 ; dim_seps(nsep) = swlen+1

    ! --------------------------------------------------------------------
    ! Hangle along the NSEP indices until we have found the two dimensions
    ! we are interested in.
    ! --------------------------------------------------------------------
    getdims:DO i = 0, nsep-1
      istart = dim_seps(i)+1 ; iend = dim_seps(i+1)-1
      IF  ( dim_chars(istart:iend) == "N_Xtr_"//l1b_channel  ) n_xtr  = dim_array(i)
      IF  ( dim_chars(istart:iend) == "N_Wvl_"//l1b_channel  ) n_wvl  = dim_array(i)
      IF  ( dim_chars(istart:iend) == "N_SoCo_Pars"          ) n_pars = dim_array(i)
      IF ( n_xtr > 0 .AND. n_wvl > 0 .AND. n_pars > 0 ) EXIT getdims
    END DO getdims

    RETURN
  END SUBROUTINE soco_get_dims

  SUBROUTINE soco_compute ( do_norm, ixt, nwvl, wvl, spc)

    USE sao_pge_utils, ONLY: array_locate_r8
    IMPLICIT NONE

    ! ------------------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    ixt .................. current cross-track position
    !    nwvl ................. number of wavlengths to be computed
    !    wvl .................. wavlength at which to compute spectrum
    !    spc .................. the solar composite spectrum
    ! ------------------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    LOGICAL,                             INTENT (IN) :: do_norm
    INTEGER (KIND=i4),                   INTENT (IN) :: ixt, nwvl
    REAL    (KIND=r8), DIMENSION (nwvl), INTENT (IN) :: wvl

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (nwvl), INTENT (OUT) :: spc

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: iwvl, iloc, ipow, nbin, npol
    REAL    (KIND=r8) :: binwvl, binpow, snorm

    ! -------------------------------
    ! Name of the function/subroutine
    ! -------------------------------
    !CHARACTER (LEN=12), PARAMETER :: modulename = 'soco_compute'

    nbin  = solarcomp_pars%nBins
    npol  = solarcomp_pars%nPol
    snorm = solarcomp_pars%SolarNorm

    ! ---------------------------------------------------------------------------
    ! Adding the contributions one-by-one and power-by-power is a bit pedestrian,
    ! but anything else involves allocating and deallocating memory, and we don't
    ! really want to go to that trouble.
    ! ---------------------------------------------------------------------------
    DO iwvl = 1, nwvl
      ! ----------------------------------------------------------------------------
      ! CAREFUL: Prevent OOB if "iloc = n_soco_bin" or "wvl = soco_wvl(n_soco_bins)"
      ! ----------------------------------------------------------------------------
      IF ( wvl(iwvl) <= solarcomp_pars%SolComParWavs(1) ) THEN
        iloc = 1 ; binwvl = -1.0_r8
      ELSE IF ( wvl(iwvl) >= solarcomp_pars%SolComParWavs(nbin) ) THEN
        iloc = nbin ; binwvl = 1.0_r8
      ELSE
        CALL array_locate_r8 ( nbin, solarcomp_pars%SolComParWavs(1:nbin), wvl(iwvl), 'LE', iloc )
        binwvl = 2.0_r8*(wvl(iwvl)-solarcomp_pars%SolComParWavs(iloc)) / &
          (solarcomp_pars%SolComParWavs(iloc+1)-solarcomp_pars%SolComParWavs(iloc)) - 1.0_r8
      END IF
      binpow     = 1.0_r8
      spc (iwvl) = solarcomp_pars%SolComParData(ixt,iloc,1)
      DO ipow = 2, npol
        binpow    = binpow * binwvl
        spc(iwvl) = spc(iwvl) + solarcomp_pars%SolComParData(ixt,iloc,ipow) * binpow
      END DO

    END DO

    ! -----------------------------------------------------------------
    ! Multiply by norm of spectrum, unless spectrum is to be normalized
    ! (do_norm = .TRUE.)
    ! -----------------------------------------------------------------
    !!! tpk(11/09/2010) IF ( .NOT. do_norm ) spc(1:nwvl) = spc(1:nwvl) * snorm
    IF ( .NOT. do_norm ) spc(1:nwvl) = spc(1:nwvl) * snorm

    ! -------------------------------------------------------------------
    ! Multiply by norm of spectrum anyway: It is better to take care of
    ! normalization in the routine that computes the final spectrum norm.
    ! -------------------------------------------------------------------
    !!!spc(1:nwvl) = spc(1:nwvl) * snorm

    RETURN
  END SUBROUTINE soco_compute

  SUBROUTINE soco_pars_deallocate ( errstat )

    USE OMSAO_errstat_module, ONLY: pge_errstat_ok
    !USE OMSAO_parameters_module,   ONLY: MAX_STR_LEN
    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    errstat .............. error status returned from the subroutine
    ! ---------------------------------------------------------------------

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    errstat = pge_errstat_ok

    ! ------------------------------------
    ! De-allocate memory for global arrays
    ! ------------------------------------
    IF ( ALLOCATED(solarcomp_pars%SolComParWavs) ) DEALLOCATE (solarcomp_pars%SolComParWavs)
    IF ( ALLOCATED(solarcomp_pars%SolComParData) ) DEALLOCATE (solarcomp_pars%SolComParData)

    RETURN
  END SUBROUTINE soco_pars_deallocate

END MODULE OMSAO_solcomp_module
