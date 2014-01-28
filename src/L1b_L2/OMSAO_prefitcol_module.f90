MODULE OMSAO_prefitcol_module

  ! =================================================================
  ! This module collects all variables and subroutines required for
  ! the ingestion of prefitted columns as may be required by the
  ! OMI OMHCHO PGE. They have been collected in this module to
  ! unclutter the main OMI PGE fitting routine.
  ! =================================================================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, nxtrack_max, nlines_max
  USE OMSAO_indices_module,    ONLY: &
    o3_t1_idx, o3_t2_idx, o3_t3_idx, bro_idx, lqh2o_idx, &
    pge_hcho_idx, pge_gly_idx

  IMPLICIT NONE

  public init_prefit_files, read_prefit_columns
  public copy_prefit_values, apply_prefit_values_and_bounds, apply_prefit_values
  type, public :: prefit_type
    real (kind=r8) :: bro_col, bro_dcol
    real (kind=r8) :: lqh2o_col, lqh2o_dcol
    real (kind=r8), dimension(o3_t1_idx:o3_t3_idx) :: o3_col, o3_dcol
  end type prefit_type

  ! -------------------------------------
  ! Logicals for use of prefitted columns
  ! -------------------------------------
  LOGICAL, public, DIMENSION (2) :: yn_bro_prefit, yn_o3_prefit, yn_lqh2o_prefit

  ! ------------------------------------------
  ! Total number of prefitted column variables
  ! ------------------------------------------
  INTEGER (KIND=i4), public :: n_prefit_vars

  ! --------------------------------------------
  ! BrO prefitted column variables
  ! --------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: bro_prefit_fname
  INTEGER (KIND=i4) :: bro_prefit_fitidx
  REAL (KIND=r8), private :: bro_prefit_var
  REAL (KIND=r8), private, DIMENSION (nxtrack_max,0:nlines_max-1) :: &
    bro_prefit_col, bro_prefit_dcol

  ! --------------------------------------------
  ! O3 prefitted column variables
  ! --------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: o3_prefit_fname
  INTEGER (KIND=i4), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3_prefit_fitidx
  REAL (KIND=r8), private, DIMENSION (o3_t1_idx:o3_t3_idx) :: o3_prefit_var
  REAL (KIND=r8), private, DIMENSION (o3_t1_idx:o3_t3_idx,nxtrack_max,0:nlines_max-1) :: &
    o3_prefit_col, o3_prefit_dcol

  ! --------------------------------------------
  ! LqH2O prefitted column variables CCM
  ! --------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: lqh2o_prefit_fname
  INTEGER (KIND=i4) :: lqh2o_prefit_fitidx
  REAL (KIND=r8), private :: lqh2o_prefit_var
  REAL (KIND=r8), private, DIMENSION (nxtrack_max,0:nlines_max-1) :: &
    lqh2o_prefit_col, lqh2o_prefit_dcol

CONTAINS

  SUBROUTINE init_prefit_files ( pge_idx, ntimes, nxtrack, errstat )

    USE OMSAO_errstat_module
    USE OMSAO_he5_module, ONLY: &
      o3fit_swath_id,    o3fit_swath_file_id,    o3fit_swath_name,  &
      brofit_swath_id,   brofit_swath_file_id,   brofit_swath_name, &
      lqh2ofit_swath_id, lqh2ofit_swath_file_id, lqh2ofit_swath_name, &
      he5_init_input_file

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ntimes, nxtrack

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ntimes_o3, nxtrack_o3, ntimes_bro, nxtrack_bro
    INTEGER (KIND=i4) :: ntimes_lqh2o, nxtrack_lqh2o
    CHARACTER (LEN=17), PARAMETER :: modulename = 'init_prefit_files'

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    o3fit_swath_id  = -1 ; o3fit_swath_file_id  = -1
    brofit_swath_id = -1 ; brofit_swath_file_id = -1
    lqh2ofit_swath_id = -1 ; lqh2ofit_swath_file_id = -1

    ! --------------------
    ! Return if no prefits
    ! --------------------
    IF ( .NOT. ANY((/yn_o3_prefit, yn_bro_prefit,yn_lqh2o_prefit/)) ) RETURN

    ! ----------------------------------
    ! Add prefits for specific retrieval
    ! ----------------------------------
    SELECT CASE( pge_idx )
    CASE ( pge_hcho_idx )

      ! ----------
      ! O3 prefits
      ! ----------
      locerrstat = pge_errstat_ok
      IF ( yn_o3_prefit(1) ) THEN
        CALL he5_init_input_file ( &
          o3_prefit_fname, o3fit_swath_name, o3fit_swath_id, o3fit_swath_file_id, &
          ntimes_o3, nxtrack_o3, errstat )
        IF ( ntimes_o3 /= ntimes .OR. nxtrack_o3 /= nxtrack ) THEN
          locerrstat = pge_errstat_error
          CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_PREFITDIM, &
            modulename//f_sep//"O3 access failed.", vb_lev_default, errstat )
          yn_o3_prefit = .FALSE.
        END IF
      END IF

      ! -----------
      ! BrO prefits
      ! -----------
      locerrstat = pge_errstat_ok
      IF ( yn_bro_prefit(1) ) THEN
        CALL he5_init_input_file ( &
          bro_prefit_fname, brofit_swath_name, brofit_swath_id, &
          brofit_swath_file_id, ntimes_bro, nxtrack_bro, locerrstat )
        IF ( ntimes_bro /= ntimes .OR. nxtrack_bro /= nxtrack ) THEN
          locerrstat = pge_errstat_error
          CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_PREFITDIM, &
            modulename//f_sep//"BrO access failed.", vb_lev_default, errstat )
          yn_bro_prefit = .FALSE.
        END IF
      END IF

    CASE ( pge_gly_idx )

      ! -------------
      ! lqH2O prefits
      ! -------------

      locerrstat = pge_errstat_ok
      IF ( yn_lqh2o_prefit(1) ) THEN
        CALL he5_init_input_file ( &
          lqh2o_prefit_fname, lqh2ofit_swath_name, lqh2ofit_swath_id, &
          lqh2ofit_swath_file_id, ntimes_lqh2o, nxtrack_lqh2o, locerrstat )
        IF ( ntimes_lqh2o /= ntimes .OR. nxtrack_lqh2o /= nxtrack ) THEN
          locerrstat = pge_errstat_error
          CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_PREFITDIM, &
            modulename//f_sep//"lqH2O access failed.", vb_lev_default, errstat )
          yn_lqh2o_prefit = .FALSE.
        END IF
      END IF

    CASE DEFAULT
      RETURN
    END SELECT

    RETURN
  END SUBROUTINE init_prefit_files

  SUBROUTINE read_prefit_columns ( pge_idx, nxtrack, nloop, iline, errstat )

    USE OMSAO_parameters_module, ONLY: r8_missval
    USE OMSAO_variables_module,  ONLY: refspecs_original
    USE datafields, ONLY: DataField_HE5, o3_prefit_he5fields, o3_prefit_uncert_he5fields
    USE OMSAO_he5_module,        ONLY:            &
      o3fit_swath_id,  &
      brofit_swath_id,   &
      lqh2ofit_swath_id
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nxtrack, nloop, iline

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, i!,iloop, it, nxtloc
    LOGICAL           :: yn_read_amf
    CHARACTER (LEN=12), PARAMETER :: col_str  = "ColumnAmount"
    CHARACTER (LEN=17), PARAMETER :: dcol_str = "ColumnUncertainty"
    INTEGER  (KIND=i4), PARAMETER :: lcolstr = LEN(col_str), ldcolstr = LEN(dcol_str)
    TYPE(DataField_HE5), pointer :: ptr_prefit, ptr_prefit_uncert

    !CHARACTER (LEN=19), PARAMETER :: modulename = 'read_prefit_columns'

    ! -----------------------------------
    ! Read prefits for specific retrieval
    ! -----------------------------------
    SELECT CASE ( pge_idx )
    CASE( pge_hcho_idx )

      ! ---------------------------------------------
      ! O3 prefitted columns and column uncertainties
      ! ---------------------------------------------
      yn_read_amf = .FALSE. ; locerrstat = pge_errstat_ok
      IF ( yn_o3_prefit(1) ) THEN
        o3_prefit_col = 0.0_r8  ;  o3_prefit_dcol = 0.0_r8

        ptr_prefit => o3_prefit_he5fields % head
        ptr_prefit_uncert => o3_prefit_uncert_he5fields % head
        i = 1
        do while (associated(ptr_prefit) .and. associated(ptr_prefit_uncert))
          CALL he5_read_prefit_columns ( &
            o3fit_swath_id, nloop, nxtrack, iline, &
            len_trim(ptr_prefit%name), trim(ptr_prefit%name), &
            o3_prefit_col (i,1:nxtrack,0:nloop-1), &
            len_trim(ptr_prefit_uncert%name), trim(ptr_prefit_uncert%name), &
            o3_prefit_dcol(i,1:nxtrack,0:nloop-1), &
            yn_read_amf, locerrstat )
          errstat = MAX( errstat, locerrstat )

          ! ----------------------------------------------------------------------
          ! Multiply O3 columns with normalization factor to return to true values
          ! ----------------------------------------------------------------------
          WHERE ( o3_prefit_col (i,1:nxtrack,0:nloop-1) > r8_missval )
            o3_prefit_col (i,1:nxtrack,0:nloop-1) = &
              o3_prefit_col (i,1:nxtrack,0:nloop-1) * refspecs_original(i)%NormFactor
            o3_prefit_dcol(i,1:nxtrack,0:nloop-1) = &
              o3_prefit_dcol(i,1:nxtrack,0:nloop-1) * refspecs_original(i)%NormFactor
          END WHERE
          ptr_prefit => ptr_prefit % next
          ptr_prefit_uncert => ptr_prefit_uncert % next
          i = i + 1
        end do
      END IF

      ! -----------------------------------------------
      ! BrO prefitted columns and column uncertainties
      ! -----------------------------------------------
      yn_read_amf = .TRUE. ; locerrstat = pge_errstat_ok
      IF ( yn_bro_prefit(1) ) THEN
        CALL he5_read_prefit_columns (                                 &
          brofit_swath_id, nloop, nxtrack, iline,                   &
          lcolstr,   col_str, bro_prefit_col (1:nxtrack,0:nloop-1), &
          ldcolstr, dcol_str, bro_prefit_dcol(1:nxtrack,0:nloop-1), &
          yn_read_amf, locerrstat )
        errstat = MAX( errstat, locerrstat )

        ! -----------------------------------------------------------------------
        ! Multiply BrO columns with normalization factor to return to true values
        ! -----------------------------------------------------------------------
        WHERE ( bro_prefit_col (1:nxtrack,0:nloop-1) > r8_missval )
          bro_prefit_col (1:nxtrack,0:nloop-1) = &
            bro_prefit_col (1:nxtrack,0:nloop-1) * refspecs_original(bro_idx)%NormFactor
          bro_prefit_dcol(1:nxtrack,0:nloop-1) = &
            bro_prefit_dcol(1:nxtrack,0:nloop-1) * refspecs_original(bro_idx)%NormFactor
        END WHERE
      END IF

    CASE ( pge_gly_idx )

      ! ------------------------------------------------
      ! lqH2O prefitted columns and column uncertainties
      ! ------------------------------------------------
      ! ccm - Retrieved "Slant Columns"
      yn_read_amf = .FALSE. ; locerrstat = pge_errstat_ok
      IF ( yn_lqh2o_prefit(1) ) THEN
        CALL he5_read_prefit_columns (                                 &
          lqh2ofit_swath_id, nloop, nxtrack, iline,                   &
          lcolstr,   col_str, lqh2o_prefit_col (1:nxtrack,0:nloop-1), &
          ldcolstr, dcol_str, lqh2o_prefit_dcol(1:nxtrack,0:nloop-1), &
          yn_read_amf, locerrstat )
        errstat = MAX( errstat, locerrstat )

        ! -------------------------------------------------------------------------
        ! Multiply lqH2O columns with normalization factor to return to true values
        ! -------------------------------------------------------------------------
        WHERE ( lqh2o_prefit_col (1:nxtrack,0:nloop-1) > r8_missval )
          lqh2o_prefit_col (1:nxtrack,0:nloop-1) = &
            lqh2o_prefit_col (1:nxtrack,0:nloop-1) * refspecs_original(lqh2o_idx)%NormFactor
          lqh2o_prefit_dcol(1:nxtrack,0:nloop-1) = &
            lqh2o_prefit_dcol(1:nxtrack,0:nloop-1) * refspecs_original(lqh2o_idx)%NormFactor
        END WHERE
      END IF

    END SELECT

    RETURN
  END SUBROUTINE read_prefit_columns

  SUBROUTINE he5_read_prefit_columns (                       &
      swath_id, ntimes_mol, nxtrack_mol, iline,             &
      mcol_len,  molcol_field, col_mol, mdcol_len,          &
      moldcol_field, dcol_mol, yn_read_amf, errstat )

    USE OMSAO_he5_module, ONLY: he5_start_2d, he5_stride_2d, he5_edge_2d, &
      HE5_SWrdfld
    USE datafields, ONLY: amfmol_field
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_read_prefit_columns'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4),       INTENT (IN) :: iline, swath_id, ntimes_mol, nxtrack_mol
    INTEGER   (KIND=i4),       INTENT (IN) :: mcol_len, mdcol_len
    CHARACTER (LEN=mcol_len) , INTENT (IN) :: molcol_field
    CHARACTER (LEN=mdcol_len), INTENT (IN) :: moldcol_field
    LOGICAL,                   INTENT (IN) :: yn_read_amf

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4),                                          INTENT (INOUT) :: errstat
    REAL    (KIND=r8), DIMENSION (nxtrack_mol, 0:ntimes_mol-1), INTENT (OUT)   :: col_mol, dcol_mol

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    REAL    (KIND=r8), DIMENSION (nxtrack_mol, 0:ntimes_mol-1) :: amf

    locerrstat = pge_errstat_ok

    ! ----------------------------------------------------
    ! Read current data block fitting output from HE5 file
    ! ----------------------------------------------------
    he5_start_2d = (/ 0, iline /) ; he5_stride_2d = (/ 1, 1 /)
    he5_edge_2d  = (/ nxtrack_mol, ntimes_mol /)

    ! -----------------------------
    ! Column amount and uncertainty
    ! -----------------------------
    locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(molcol_field)),     &
      he5_start_2d, he5_stride_2d, he5_edge_2d, col_mol(1:nxtrack_mol,0:ntimes_mol-1) )
    locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(moldcol_field)),    &
      he5_start_2d, he5_stride_2d, he5_edge_2d, dcol_mol(1:nxtrack_mol,0:ntimes_mol-1) )

    ! --------------------------------------
    ! Air Mass Factor, but only if requested
    ! (e.g., for BrO, but not O3)
    ! --------------------------------------
    IF ( yn_read_amf ) THEN
      locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(amfmol_field)),     &
        he5_start_2d, he5_stride_2d, he5_edge_2d, amf(1:nxtrack_mol,0:ntimes_mol-1) )
      WHERE ( amf > 0.0_r8 )
        col_mol  =  col_mol * amf
        dcol_mol = dcol_mol * amf
      END WHERE
    END IF

    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, &
      OMSAO_E_HE5SWRDFLD, modulename, vb_lev_default, errstat )

    RETURN
  END SUBROUTINE he5_read_prefit_columns

!UNUSED!   FUNCTION he5_close_prefit_file ( swath_id, swath_file_id ) RESULT ( he5stat )
!UNUSED!
!UNUSED!     !------------------------------------------------------------------------------
!UNUSED!     ! This function detatches from the HE5 swath and closes the HE5 input file.
!UNUSED!     !
!UNUSED!     ! Input:
!UNUSED!     !
!UNUSED!     !   swath_id       - ID of pre-fitted input swath      (for BrO or O3)
!UNUSED!     !   swath_file_id  - ID of pre-fitted input swath file (for BrO or O3)
!UNUSED!     !
!UNUSED!     !------------------------------------------------------------------------------
!UNUSED!
!UNUSED!     USE OMSAO_he5_module
!UNUSED!     USE OMSAO_errstat_module
!UNUSED!
!UNUSED!     IMPLICIT NONE
!UNUSED!
!UNUSED!     ! ---------------------------------------
!UNUSED!     ! Name of this module/subroutine/function
!UNUSED!     ! ---------------------------------------
!UNUSED!     CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_close_prefit_file'
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Input variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4), INTENT (IN) :: swath_id, swath_file_id
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Result variable
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4) :: he5stat
!UNUSED!
!UNUSED!     ! --------------
!UNUSED!     ! Local variable
!UNUSED!     ! --------------
!UNUSED!     INTEGER (KIND=i4) :: locerr
!UNUSED!
!UNUSED!     he5stat = pge_errstat_ok
!UNUSED!     locerr  = pge_errstat_ok
!UNUSED!
!UNUSED!     ! -----------------------------------------------
!UNUSED!     ! Detach from HE5 swath and close HE5 output file
!UNUSED!     ! -----------------------------------------------
!UNUSED!     locerr = HE5_SWdetach ( swath_id )
!UNUSED!     locerr = HE5_SWclose  ( swath_file_id )
!UNUSED!     CALL error_check ( locerr, HE5_STAT_OK, pge_errstat_warning, &
!UNUSED!       OMSAO_W_HE5SWCLOSE, modulename, vb_lev_default, he5stat )
!UNUSED!
!UNUSED!     RETURN
!UNUSED!   END FUNCTION he5_close_prefit_file

  subroutine copy_prefit_values (prefit, pge_idx, ipix, iloop)
    implicit none
    type (prefit_type), intent(inout) :: prefit
    integer (kind=i4), intent(in) :: pge_idx, ipix, iloop

    select case (pge_idx)
    case (pge_hcho_idx)
      prefit%o3_col = o3_prefit_col (o3_t1_idx:o3_t3_idx,ipix,iloop)
      prefit%o3_dcol = o3_prefit_dcol (o3_t1_idx:o3_t3_idx,ipix,iloop)
      prefit%bro_col = bro_prefit_col (ipix,iloop)
      prefit%bro_dcol = bro_prefit_dcol (ipix,iloop)
    case (pge_gly_idx)
      prefit%lqh2o_col = lqh2o_prefit_col (ipix,iloop)
      prefit%lqh2o_dcol = lqh2o_prefit_dcol (ipix,iloop)
    case default
      !nothing
    end select

  end subroutine copy_prefit_values

  subroutine apply_prefit_values (fitvar)
    implicit none
    real (kind=r8), dimension(:), intent(inout) :: fitvar
    integer (kind=i4) :: j

    IF ( yn_lqh2o_prefit(1) .AND. (.NOT. yn_lqh2o_prefit(2)) .AND. lqh2o_prefit_var /= 0.0_r8 ) &
      fitvar(lqh2o_prefit_fitidx) = lqh2o_prefit_var

    IF ( yn_bro_prefit(1) .AND. (.NOT. yn_bro_prefit(2)) .AND. bro_prefit_var /= 0.0_r8 ) &
      fitvar(bro_prefit_fitidx) = bro_prefit_var

    IF ( yn_o3_prefit(1)  .AND. (.NOT. yn_o3_prefit(2))  ) THEN
      DO j = o3_t1_idx, o3_t3_idx
        IF ( o3_prefit_var(j) /= 0.0_r8 ) fitvar(o3_prefit_fitidx(j)) = o3_prefit_var(j)
      END DO
    END IF

  end subroutine apply_prefit_values

  subroutine apply_prefit_values_and_bounds (prefit, pge_idx, lobnd, upbnd, fitvar, errstat)
    use OMSAO_parameters_module, only: r8_missval
    implicit none
    type (prefit_type), intent(in) :: prefit
    integer (kind=i4), intent(in) :: pge_idx
    real (kind=r8), dimension(:), intent(inout) :: lobnd, upbnd, fitvar
    integer, intent(inout) :: errstat
    ! local
    real (kind=r8) :: col_tmp, dcol_tmp
    integer (kind=i4) :: j

    if (errstat < 0) return

    select case (pge_idx)
      case (pge_hcho_idx)
        ! BrO
        IF ( yn_bro_prefit(1) ) THEN
          bro_prefit_var = 0.0_r8
          col_tmp = prefit%bro_col
          dcol_tmp = prefit%bro_dcol
          IF ( col_tmp > r8_missval ) THEN
            fitvar(bro_prefit_fitidx) = col_tmp
            IF ( yn_bro_prefit(2) ) THEN
              lobnd(bro_prefit_fitidx) = col_tmp - 1.0_r8 * dcol_tmp
              upbnd(bro_prefit_fitidx) = col_tmp + 1.0_r8 * dcol_tmp
            ELSE
              lobnd(bro_prefit_fitidx) = col_tmp
              upbnd(bro_prefit_fitidx) = col_tmp
              bro_prefit_var           = col_tmp
            END IF
          END IF
        END IF
        ! O3
        IF ( yn_o3_prefit(1) ) THEN
          o3_prefit_var(o3_t1_idx:o3_t3_idx) = 0.0_r8
          DO j = o3_t1_idx, o3_t3_idx
            col_tmp = prefit%o3_col(j)
            dcol_tmp = prefit%o3_dcol(j)
            IF ( col_tmp > r8_missval ) THEN
              fitvar(o3_prefit_fitidx(j)) = col_tmp
              IF ( yn_o3_prefit(2) ) THEN
                lobnd(o3_prefit_fitidx(j)) = col_tmp - 2.0_r8 * dcol_tmp
                upbnd(o3_prefit_fitidx(j)) = col_tmp + 2.0_r8 * dcol_tmp
              ELSE
                lobnd(o3_prefit_fitidx(j)) = col_tmp
                upbnd(o3_prefit_fitidx(j)) = col_tmp
                o3_prefit_var(j)           = col_tmp
              END IF
            END IF
          END DO
        END IF
      case (pge_gly_idx)
        ! liquid water
        IF ( yn_lqh2o_prefit(1) ) THEN
          lqh2o_prefit_var = 0.0_r8
          col_tmp = prefit%lqh2o_col
          dcol_tmp = prefit%lqh2o_dcol
          IF ( col_tmp > r8_missval ) THEN
            fitvar(lqh2o_prefit_fitidx) = col_tmp
            IF ( yn_lqh2o_prefit(2) ) THEN
              lobnd(lqh2o_prefit_fitidx) = col_tmp - 1.0_r8 * dcol_tmp
              upbnd(lqh2o_prefit_fitidx) = col_tmp + 1.0_r8 * dcol_tmp
            ELSE
              lobnd(lqh2o_prefit_fitidx) = col_tmp
              upbnd(lqh2o_prefit_fitidx) = col_tmp
              lqh2o_prefit_var           = col_tmp
            END IF
          END IF
        END IF
      case default
        ! nothing
    end select

  end subroutine apply_prefit_values_and_bounds

END MODULE OMSAO_prefitcol_module
