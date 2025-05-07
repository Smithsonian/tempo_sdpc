!
module he5_input_tools

  public he5_init_input_file
  private

contains

  SUBROUTINE he5_init_input_file ( file_name, swath_name, swath_id, &
       swath_file_id, ntimes_aux, nxtrack_aux, he5stat )

    !--------------------------------------------------------------------------
    ! This subroutine initializes the HE5 input files for for
    ! prefitted BrO and O3 that are used in OMHCHO
    !
    ! Input:
    !   file_name         - Name of HE5 input file
    !   swath_name        - Name of existing swath in file
    !   swath_file_id_inp - id number for HE5 input file (required to close it)
    !   swath_id_inp      - id number for swath (required to read from swath)
    !
    ! Output:
    !   ntimes_aux.............nTimes  as given in product file
    !   nxtrack_aux............nXtrack as given in product file
    !   he5stat................OMI_E_SUCCESS if everything went well
    !
    ! No Swath ID Variables passed through MODULE.
    !
    !--------------------------------------------------------------------------

    !USE OMSAO_indices_module,    ONLY: pge_hcho_idx, pge_bro_idx, &
         !max_calfit_idx
    USE OMSAO_parameters_module, ONLY: maxchlen, vb_lev_default!, &
         !r8_missval, r4_missval, i4_missval
    USE OMSAO_he5_module
    USE OMSAO_errstat_module
    USE OMSAO_precision_module,  ONLY: i4, r8
    use ISO_C_BINDING, only: C_LONG

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
    INTEGER   (KIND=i4),      INTENT (OUT)   :: swath_id, swath_file_id, &
         nxtrack_aux, ntimes_aux
    CHARACTER (LEN=maxchlen), INTENT (OUT)   :: swath_name

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)   :: ndim, nsep
    INTEGER   (KIND=i4)   :: errstat!, iend, istart, swlen, i, j, astat
    integer (kind=r8) :: i, istart, iend
    integer (kind=C_LONG) :: swlen
    !    INTEGER   (KIND=i4),      DIMENSION(0:9) :: dim_array, dim_seps
    INTEGER   (KIND=C_LONG),      DIMENSION(0:9) :: dim_array, dim_seps
    CHARACTER (LEN=maxchlen)                 :: dim_chars!, tmp_char

    !REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: o3fit_cols, o3fit_dcols
    !REAL (KIND=r8), DIMENSION (:,:),   ALLOCATABLE :: tmp_col

    errstat = pge_errstat_ok

    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_rdonly )
    IF ( swath_file_id == he5_stat_fail ) THEN
      CALL error_check ( &
           0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, &
           vb_lev_default, he5stat )
      IF ( he5stat >= pge_errstat_error ) RETURN
    END IF

    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    errstat  = HE5_SWinqswath  ( TRIM(ADJUSTL(file_name)), swath_name, swlen )

    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == he5_stat_fail ) THEN
      CALL error_check ( &
           0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, &
           vb_lev_default, he5stat )
      IF ( he5stat >= pge_errstat_error ) RETURN
    END IF

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    errstat  = HE5_SWinqdims  ( swath_id, dim_chars, dim_array(0:9) )


    ! -------------------------------------------------------------
    ! The number of recorded dimension is the number of non-zero
    ! entries - i.a.w., it't the index of the first ZERO minus ONE.
    ! -------------------------------------------------------------
    ndim = MINVAL( MINLOC( dim_array, MASK=(dim_array == 0) ) ) - 1

    ! -----------------------------------------
    ! Extract nTimes and nXtrack from the swath
    ! -----------------------------------------
    nxtrack_aux = 0  
    ntimes_aux = 0
    dim_chars = TRIM(ADJUSTL(dim_chars))
    swlen = LEN_TRIM(ADJUSTL(dim_chars))
    istart = 1  
    iend = 1

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    nsep = 0 
    dim_seps = 0 
    dim_seps(0) = 0
    getseps: DO i = 1, swlen
      IF ( dim_chars(i:i) == ',' ) THEN
        nsep = nsep + 1
        dim_seps(nsep) = i
      END IF
    END DO getseps
    nsep = nsep + 1 
    dim_seps(nsep) = swlen+1


    ! --------------------------------------------------------------------
    ! Hangle along the NSEP indices until we have found the two dimensions
    ! we are interested in.
    ! --------------------------------------------------------------------
    getdims:DO i = 0, nsep-1

      istart = dim_seps(i)+1 
      iend = dim_seps(i+1)-1
      IF  ( dim_chars(istart:iend) == "nTimes"  ) ntimes_aux  = INT(dim_array(i), kind=4)
      IF  ( dim_chars(istart:iend) == "nXtrack" ) nxtrack_aux = INT(dim_array(i), kind=4)
      IF ( ntimes_aux > 0 .AND. nxtrack_aux > 0 ) EXIT getdims
    END DO getdims

    RETURN
  END SUBROUTINE he5_init_input_file


end module he5_input_tools
