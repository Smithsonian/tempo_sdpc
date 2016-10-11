MODULE OMSAO_he5_module

  ! ==============================
  ! Module for HDF-EOS5 parameters
  ! ==============================

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: &
    MAX_STR_LEN, i2_missval, i4_missval, r4_missval, r8_missval, blank13,       &
    valid_max_i2, valid_max_i4, valid_min_r8, valid_max_r8, zero_r8, one_r8, &
    main_qa_min_flag_r8, main_qa_max_flag_r8
  USE OMSAO_indices_module,    ONLY: &
    max_calfit_idx, sao_pge_min_idx, sao_pge_max_idx, &
    o3_t1_idx, o3_t3_idx

  IMPLICIT NONE

  ! ----------------------------------------------
  ! Swath IDs for current HE5 swath and swath file
  ! ----------------------------------------------
  INTEGER (KIND=i4), PUBLIC :: pge_swath_file_id, pge_swath_id

  ! --------------------------------------------------------------
  ! Swath IDs for pre-fitted BrO and O3 HE5 swaths and swath files
  ! --------------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN), PUBLIC :: &
    o3fit_swath_name    = 'undefined', &
    brofit_swath_name   = 'undefined', &
    lqh2ofit_swath_name = 'undefined'
  INTEGER (KIND=i4), PUBlIC :: &
    o3fit_swath_id    = -1, o3fit_swath_file_id    = -1, &
    brofit_swath_id   = -1, brofit_swath_file_id   = -1, &
    lqh2ofit_swath_id = -1, lqh2ofit_swath_file_id = -1

  ! --------------------
  ! Fill Value attribute
  ! --------------------
  CHARACTER (LEN=*), PUBLIC, PARAMETER :: &
    missval_attr = "MissingValue", &
    offset_attr  = "Offset", &
    rstemp_attr  = "ReferenceSpectrumTemperature", &
    title_attr   = "Title", &
    ufd_attr     = "UniqueFieldDefinition", &
    units_attr   = "Units", &
    valids_attr  = "ValidRange", &
    scafac_attr  = "ScaleFactor"

  ! ----------------------------------
  ! Strings for Global File Attributes
  ! ----------------------------------
  CHARACTER (LEN=*), PUBLIC, PARAMETER :: tai_attr = "TAI93At0zOfGranule"

  ! -----------------------------------------------
  ! Variables that will hold Global File Attributes
  ! -----------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN), PUBLIC :: &
    pge_swath_name, process_level, instrument_name, pge_version, l1b_orbitdata
  INTEGER (KIND=i4), PUBLIC :: granule_day, granule_month, granule_year
  REAL (KIND=r8), PUBLIC    :: TAI93At0zOfGranule

  ! ----------------------
  ! Swath Level Attributes
  ! ----------------------
  CHARACTER (LEN=18), PUBLIC, PARAMETER :: vcoordinate_field = "VerticalCoordinate"

  ! ------------------------------------------------------------------
  ! Integer variables for writing to swath data and geolocation fields
  ! ------------------------------------------------------------------
  INTEGER (KIND=C_LONG), PUBLIC, DIMENSION(1) :: he5_start_1d, he5_stride_1d, he5_edge_1d
  INTEGER (KIND=C_LONG), PUBLIC, DIMENSION(2) :: he5_start_2d, he5_stride_2d, he5_edge_2d
  INTEGER (KIND=C_LONG), PUBLIC, DIMENSION(3) :: he5_start_3d, he5_stride_3d, he5_edge_3d
  INTEGER (KIND=C_LONG), PUBLIC, DIMENSION(4) :: he5_start_4d, he5_stride_4d, he5_edge_4d

  ! -----------------------------------
  ! The Vertical Coordinate of the PGEs
  ! -----------------------------------
  ! Since we are doing nadir, this is either "Slant Column" or "Total Column"
  ! -------------------------------------------------------------------------
  CHARACTER (LEN=12), PUBLIC, &
    DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
    vertical_coordinate = (/                                                 &
    "Slant Column", &  ! OClO
    "Total Column", &  ! BrO
    "Total Column", &  ! HCHO
    "Slant Column", &  ! O3
    "Slant Column", &  ! NO2
    "Slant Column", &  ! SO2
    "Total Column", &  ! CHO-CHO
    "Slant Column", &  ! IO
    "Slant Column", &  ! H2O
    "Slant Column", &  ! HONO
    "Slant Column", &  ! O2O2
    "Slant Column", &  ! LqH2O
    "Slant Column"  /) ! NO2(D)

  ! -------------------------------
  ! Base names for L2 output swaths
  ! ------------------------------------------------------------------
  ! These Base Name are appended with the PGE target molecules during
  ! the initialization of the HE5 output swath name. They will be
  ! accessed as SWATH_BASE_NAME(pge_idx), and so have to correspond to
  ! the order of PGE indices defined in OMSAO_indices_module:
  !
  !        [1] OClO  [2] BrO  [3] HCHO  [4] O3
  !
  ! The "unofficial" ones are
  !
  !        [5] NO2   [6] SO2  [7] C2H2O2 (Glyoxal)  [8] IO  [9] NO2
  !
  ! Note that these strings will be used IF AND ONLY IF we cannot find
  ! the Swath Name in the PCF file.
  ! ------------------------------------------------------------------
  CHARACTER (LEN=23), PUBLIC, &
    DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
    swath_base_name = (/ &
    "OMI Slant Column Amount", &  ! OClO
    "OMI Total Column Amount", &  ! BrO
    "OMI Total Column Amount", &  ! HCHO
    "OMI Slant Column Amount", &  ! O3
    "OMI Slant Column Amount", &  ! NO2
    "OMI Slant Column Amount", &  ! SO2
    "OMI Total Column Amount", &  ! CHO-CHO
    "OMI Slant Column Amount", &  ! IO
    "OMI Slant Column Amount", &  ! H2O
    "OMI Slant Column Amount", &  ! HONO
    "OMI Slant Column Amount", &  ! O2O2
    "OMI Slant Column Amount", &  ! LqH2O
    "OMI Slant Column Amount"  /) ! NO2(D)

  ! -------------------------
  ! PGE HDF Global Attributes
  ! ---------------------------------------------------------------
  ! NOTE: These values MIGHT ventually become PSAs, because we want
  !       to search for them in the data set. For now we keep them
  !       as Global Attributes, since this is a much more painless
  !       state than any newly defined PSA.
  ! ---------------------------------------------------------------
  !INTEGER (KIND=i4), PUBLIC:: &
  !  NrofScanLines                   = i4_missval,  &
  !  NrofCrossTrackPixels            = i4_missval,  &
  !  NrofInputSamples                = i4_missval,  &
  !  NrofGoodInputSamples            = i4_missval,  &
  !  NrofGoodOutputSamples           = i4_missval,  &
  !  NrofMissingSamples              = i4_missval,  &
  !  NrofSuspectOutputSamples        = i4_missval,  &
  !  NrofBadOutputSamples            = i4_missval,  &
  !  NrofConvergedSamples            = i4_missval,  &
  !  NrofFailedConvergenceSamples    = i4_missval,  &
  !  NrofExceededIterationsSamples   = i4_missval,  &
  !  NrofOutofBoundsSamples          = i4_missval
  !REAL (KIND=r4), PUBLIC :: &
  !  PercentGoodOutputSamples        = r4_missval,  &
  !  PercentSuspectOutputSamples     = r4_missval,  &
  !  PercentBadOutputSamples         = r4_missval,  &
  !  AbsolutePercentMissingSamples   = r4_missval

  ! --------------------------------------------
  ! Variables for InputPointer and InputVersions
  ! --------------------------------------------
  INTEGER   (KIND=i4), PUBLIC, PARAMETER :: n_lun_inp_max = 20
  INTEGER   (KIND=i4), PUBLIC, DIMENSION (n_lun_inp_max) :: lun_input
  INTEGER   (KIND=i4), PUBLIC :: n_lun_inp
  CHARACTER (LEN=MAX_STR_LEN), PUBLIC            :: input_versions

  ! ---------------------------------------------------------------
  ! Finally some system definitions that come with the HE5 Library.
  ! We include this here so that we don't have to worry about it
  ! inside the subroutines that use it.
  ! ---------------------------------------------------------------
  INCLUDE 'hdfeos5.inc'
  INTEGER (KIND = i4), PUBLIC, EXTERNAL :: &
    he5_ehrdglatt,  he5_ehwrglatt,  he5_swattach,  he5_swclose,    he5_swcreate,  &
    he5_swdefdfld,  he5_swdefdim,   he5_swdefgfld, he5_swdetach,   he5_swopen,    &
    he5_swrdattr,   he5_swrdfld,    he5_swrdgattr, he5_swrdlattr,  he5_swwrattr,  &
    he5_swwrfld,    he5_swwrgattr,  he5_swwrlattr, HE5_SWfldinfo,                 &
    HE5_SWdefcomp,  HE5_SWdefchunk, HE5_SWsetfill, HE5_SWdefcomch,  HE5_GDOPEN,   &
    HE5_GDATTACH,   HE5_GDRDFLD,    HE5_GDCLOSE,   HE5_GDINQFLDS,   HE5_GDDETACH, &
    HE5_GDINQLATTRS,HE5_GDRDLATTR,  HE5_SWinqdflds

  INTEGER (KIND=C_LONG), PUBLIC, EXTERNAL :: HE5_SWinqswath, HE5_SWinqdims

  ! ---------------------------------------------------------------------------
  ! Parameters for HE5 compression. From the limited experimentation performed,
  ! this combination has shown to produce the smallest file sizes for SAO PGEs.
  ! ---------------------------------------------------------------------------
  INTEGER (KIND=i4), PUBLIC,  PARAMETER :: &
    he5_comp_type   = HE5_HDFE_COMP_SHUF_DEFLATE, &
    he5_comp_par    = 9,                          &
    he5_nocomp_type = HE5_HDFE_COMP_NONE,         &
    he5_nocomp_par  = 0

  PRIVATE
  PUBLIC he5_init_input_file, &
    HE5T_NATIVE_CHAR, HE5T_NATIVE_INT, HE5T_NATIVE_FLOAT, HE5T_NATIVE_DOUBLE, &
    HE5T_NATIVE_INT8, HE5T_NATIVE_INT16, &
    he5f_acc_rdonly, he5_hdfe_nomerge, he5f_acc_trunc

CONTAINS

  SUBROUTINE he5_init_input_file ( &
      file_name, swath_name, swath_id, swath_file_id, ntimes_aux, nxtrack_aux, he5stat )

    !------------------------------------------------------------------------------
    ! This subroutine initializes the HE5 input files for reading
    ! prefitted BrO and O3 that are used in OMHCHO
    !
    ! Input:
    !   file_name         - Name of HE5 input file
    !
    ! Output:
    !   swath_name.............Name of existing swath in file
    !   swath_file_id_inp......id number for HE5 input file (required for closing it)
    !   swath_id_inp...........id number for swath (required for reading from swath)
    !   ntimes_aux.............nTimes  as given in product file
    !   nxtrack_aux............nXtrack as given in product file
    !   he5stat................OMI_E_SUCCESS if everything went well
    !
    ! No Swath ID Variables passed through MODULE.
    !
    !------------------------------------------------------------------------------



    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=19), PARAMETER :: modulename = 'he5_init_input_file'

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*), INTENT (IN) :: file_name

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER   (KIND=i4),      INTENT (INOUT) :: he5stat
    INTEGER   (KIND=i4),      INTENT (OUT)   :: swath_id, swath_file_id, nxtrack_aux, ntimes_aux
    CHARACTER (LEN=MAX_STR_LEN), INTENT (OUT)   :: swath_name

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                       :: nsep, ndim
    INTEGER   (KIND=i4)                       :: i, errstat, iend, istart
    integer (kind=C_LONG) :: swlen
    INTEGER   (KIND=C_LONG)                   :: ndimcl, errstatcl
    INTEGER   (KIND=i4),      DIMENSION(0:12) :: dim_array, dim_seps
    INTEGER   (KIND=C_LONG),  DIMENSION(0:12) :: dim_arraycl
    CHARACTER (LEN=MAX_STR_LEN)                  :: dim_chars

    errstat = pge_errstat_ok

    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_rdonly )
    IF ( swath_file_id == he5_stat_fail ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, he5stat )
      IF ( he5stat >= pge_errstat_error ) RETURN
    END IF

    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    swath_name=""
    errstatcl = HE5_SWinqswath  ( TRIM(ADJUSTL(file_name)), swath_name, swlen )
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == he5_stat_fail ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )
      IF ( he5stat >= pge_errstat_error ) RETURN
    END IF

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    dim_chars="" !JED
    ndimcl = HE5_SWinqdims  ( swath_id, dim_chars, dim_arraycl(0:12) )
    ndim   = INT ( ndimcl, KIND=i4 )
    IF ( ndim <= 0 ) THEN
      he5stat = MAX ( he5stat, pge_errstat_error )
      RETURN
    END IF
    dim_array(0:12) = INT ( dim_arraycl(0:12), KIND=i4 )

    ! -----------------------------------------
    ! Extract nTimes and nXtrack from the swath
    ! -----------------------------------------
    nxtrack_aux = 0  ;  ntimes_aux = 0
    dim_chars = TRIM(ADJUSTL(dim_chars))
    swlen = LEN_TRIM(ADJUSTL(dim_chars))
    istart = 1  ;  iend = 1

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    nsep = 0 ; dim_seps = 0 ; dim_seps(0) = 0
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
      IF  ( dim_chars(istart:iend) == "nTimes"  ) ntimes_aux  = dim_array(i)
      IF  ( dim_chars(istart:iend) == "nXtrack" ) nxtrack_aux = dim_array(i)
      IF ( ntimes_aux > 0 .AND. nxtrack_aux > 0 ) EXIT getdims
    END DO getdims

    RETURN
  END SUBROUTINE he5_init_input_file

END MODULE OMSAO_he5_module
