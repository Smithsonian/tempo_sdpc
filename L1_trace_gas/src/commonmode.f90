module commonmode
  use tell_module
  implicit none
  private
  public compute_common_mode, finalize_common_mode
contains

  SUBROUTINE finalize_common_mode (xti)

    USE OMSAO_precision_module, ONLY: i4, r8
    USE OMSAO_indices_module,   ONLY: max_calfit_idx, comm_idx, mxs_idx
    USE OMSAO_variables_module, ONLY:                                           &
      common_mode_spec, fitvar_rad_init, lo_radbnd, up_radbnd,               &
      common_fitpos, common_fitvar, refspecs_original
    USE OMSAO_omidata_module,   ONLY:                                           &
      n_omi_database_wvl, omi_database

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                   INTENT (IN) :: xti

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=i4) :: i, j, k
    REAL    (KIND=r8) :: comnorm

    !CHARACTER (LEN=19), PARAMETER :: modulename = 'compute_common_mode'

    ! ---------------------------------------------------
    ! Set the index value of the Common Mode spectrum and
    ! assign values to the fitting parameter arrays
    ! ---------------------------------------------------
    i = max_calfit_idx + (comm_idx-1)*mxs_idx + common_fitpos
    fitvar_rad_init(i) = common_fitvar(1)
    lo_radbnd      (i) = common_fitvar(2)
    up_radbnd      (i) = common_fitvar(3)
    DO i = 1, xti  ! NOTE: "xti == nxtrack" for this call
      j = n_omi_database_wvl(i)

      ! ------------------------------------------
      ! Average the wavelength and spectrum arrays
      ! ------------------------------------------
      k = MAX(1,common_mode_spec%RefSpecCount(i))
      common_mode_spec%RefSpecWavs(1:j,i)  = &
        common_mode_spec%RefSpecWavs(1:j,i) / REAL(k, KIND=r8)
      common_mode_spec%RefSpecData(1:j,i)  = &
        common_mode_spec%RefSpecData(1:j,i) / REAL(k, KIND=r8)

      ! ---------------------------------------
      ! Normalize the Common Mode Spectrum to 1
      ! ---------------------------------------
      ! Skip this for now until we bother with excluding the the
      ! low weights, which otherwise skew the norm.
      ! --------------------------------------------------------
      comnorm = 1.0_r8
      !comnorm = SUM(common_mode_spec%RefSpecData(1:j,i)) / REAL(k, KIND=r8)
      !IF ( comnorm == 0.0_r8 ) comnorm = 1.0_r8
      !common_mode_spec%RefSpecData(1:j,i) = &
      !     common_mode_spec%RefSpecData(1:j,i) / comnorm

      ! -------------------------------------------------
      ! Assign the common mode to the OMI data base array
      ! -------------------------------------------------
      omi_database(1:j,i,comm_idx) = common_mode_spec%RefSpecData(1:j,i)

      ! --------------------------------------------------------------
      ! Now assign a normalization factor to the original data base of
      ! reference spectra. This is needed in the computation of the
      ! columns and uncertainties of all fitting parameters.
      ! --------------------------------------------------------------
      refspecs_original(comm_idx)%NormFactor = 1.0_r8
    END DO

    !WRITE (88,'(A,2I6)'), '#', 36, n_omi_database_wvl(36)
    !DO i = 1, n_omi_database_wvl(36)
    !   WRITE (88,'(0PF10.3, 1PE15.5)') &
    !        common_mode_spec%RefSpecWavs(i,36), common_mode_spec%RefSpecData(i,36)
    !END DO

    ! ------------------------------------
    ! That is all we do for the final call
    ! ------------------------------------
    return

  end subroutine finalize_common_mode

  SUBROUTINE compute_common_mode ( &
      yn_reference_fit, xti, nwvl, fitwvl, fitres, errstat)

    USE OMSAO_precision_module, ONLY: i2, i4, r8
    USE OMSAO_indices_module,   ONLY: max_calfit_idx, comm_idx, mxs_idx
    USE OMSAO_variables_module, ONLY:                                           &
      common_mode_spec, fitvar_rad_init, lo_radbnd, up_radbnd,               &
      common_fitpos, common_latrange
    USE OMSAO_omidata_module,   ONLY:                                           &
      rad_ccdpix_selection, omi_blockline_no, omi_latitude, n_comm_wvl

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    LOGICAL,                             INTENT (IN) :: yn_reference_fit
    INTEGER (KIND=i4),                   INTENT (IN) :: xti, nwvl
    REAL    (KIND=r8), DIMENSION (nwvl), INTENT (IN) :: fitwvl, fitres
    integer, intent(inout) :: errstat

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=i4) :: i
    REAL    (KIND=r8) :: comnorm
    character (len=128) :: errmsg

    !CHARACTER (LEN=19), PARAMETER :: modulename = 'compute_common_mode'

    if (errstat < 0) return

    IF ( yn_reference_fit ) THEN
      ! -------------------------------------------------------------
      ! The Radiance Reference Fit branch saves the wavelength values
      ! and initializes the count and spectrum arrays for the current
      ! cross track position. For goot measure, we also set the
      ! common mode spectrum fitting parameters to Zero.
      ! -------------------------------------------------------------

      ! ---------------------------------------------------
      ! Set the index value of the Common Mode spectrum and
      ! assign values to the fitting parameter arrays
      ! ---------------------------------------------------
      i = max_calfit_idx + (comm_idx-1)*mxs_idx + common_fitpos
      fitvar_rad_init(i) = 0.0_r8
      lo_radbnd      (i) = 0.0_r8
      up_radbnd      (i) = 0.0_r8

      ! -------------------------------------------------------------
      ! Note that we set COMMON_CNT and COMMON_SPC to 0 for _all_
      ! positions _all_ the time. This is a safety measure just in
      ! case one of the x-track positions is skipped during the
      ! reference fit due to non-convergence. The Common Mode can
      ! function without wavelengths (actually, we still have to
      ! create a good rationale for them), but not without both
      ! spectrum and cound arrays starting from Zero
      ! -------------------------------------------------------------
      !common_wvl(xti, 1:nwvl) = fitwvl(1:nwvl)
      !common_cnt              = 0_i4
      !common_spc              = 0.0_r8

      common_mode_spec%nPoints      = n_comm_wvl
      common_mode_spec%RefSpecWavs  = 0.0_r8
      common_mode_spec%RefSpecData  = 0.0_r8
      common_mode_spec%RefSpecCount = 0

      common_mode_spec%CCDPixel(xti,1) = INT(rad_ccdpix_selection(xti,1), KIND=i2)
      common_mode_spec%CCDPixel(xti,2) = INT(rad_ccdpix_selection(xti,4), KIND=i2)
    ELSE

      ! --------------------------------------------------------
      ! The Reguar Fitting branch updates the spectrum and count
      ! --------------------------------------------------------
      IF ( omi_latitude(xti,omi_blockline_no) >= common_latrange(1) .AND. &
          omi_latitude(xti,omi_blockline_no) <= common_latrange(2)         )  THEN

        ! the summed spectra should all have the same number of data points
        if (common_mode_spec % num_wavelengths(xti) == 0) then
          common_mode_spec % num_wavelengths(xti) = nwvl
        else if (nwvl /= common_mode_spec % num_wavelengths(xti)) then
          write (errmsg, '(a,i4,a,i4,a,i4)')'compute_common_mode:  (xti=',xti, &
            ') expected nwvl=', common_mode_spec % num_wavelengths(xti), &
            ', got nwvl=',nwvl
          call tell_error (tell_runtime_error, errmsg, errstat)
          return
        endif

        comnorm = 1.0_r8
        IF ( nwvl > 0 ) THEN
          comnorm = SUM(ABS(fitres(1:nwvl)))/REAL(nwvl, KIND=r8)
          IF ( comnorm == 0.0_r8 ) comnorm = 1.0_r8
        END IF

        !common_cnt(xti)        = common_cnt(xti) + 1
        !common_spc(xti,1:nwvl) = common_spc(xti,1:nwvl) + fitres(1:nwvl)/comnorm

        common_mode_spec%RefSpecWavs(1:nwvl,xti)  = &
          common_mode_spec%RefSpecWavs(1:nwvl,xti) + fitwvl(1:nwvl)
        common_mode_spec%RefSpecData(1:nwvl,xti)  = &
          common_mode_spec%RefSpecData(1:nwvl,xti) + fitres(1:nwvl)/comnorm
        common_mode_spec%RefSpecCount(xti)        = &
          common_mode_spec%RefSpecCount(xti) + 1

      END IF
    END IF

    RETURN
  END SUBROUTINE compute_common_mode
end module
