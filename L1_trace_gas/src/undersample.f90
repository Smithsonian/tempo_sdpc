MODULE undersample
CONTAINS
SUBROUTINE undersample_spectrum ( xtrack_pix, n_sensor_pts, curr_wvl, hw1e, e_asym, k, phase, errstat )

  !     Convolves input spectrum with Gaussian slit function of specified
  !     HW1e, and samples at a particular input phase to give the OMI
  !     undersampling spectrum. This version calculates both phases of the
  !     undersampling spectrum, phase1 - i.e., underspec (i, 1) - being the
  !     more common in OMI spectra.

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: solar_idx, us1_idx, us2_idx
  USE OMSAO_variables_module,  ONLY: &
    refspecs_original, database, have_undersampling
  use slitfunction, only : slitfunction_convolve
  !USE OMSAO_errstat_module
  USE sao_pge_utils, ONLY: interpolation
  use tell_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                           INTENT (IN) :: n_sensor_pts, xtrack_pix
  REAL    (KIND=r8),                           INTENT (IN) :: hw1e, e_asym, k, phase
  REAL    (KIND=r8), DIMENSION (n_sensor_pts), INTENT (IN) :: curr_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  LOGICAL                                     :: did_full_range
  REAL (KIND=r8), DIMENSION (n_sensor_pts,2)  :: underspec
  !REAL (KIND=r8), DIMENSION (max_spec_pts)    :: &
  !  locwvl, locspec, specmod, tmpwav, over, under, resample
  REAL (KIND=r8), DIMENSION (:), allocatable    :: &          ! JCH
    locwvl, locspec, specmod, tmpwav, over, under, resample

  INTEGER (KIND=i4) :: npts, locerrstat
  integer :: num_alloc

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  !CHARACTER (LEN=11), PARAMETER :: modulename = 'undersample'

  if (errstat /= 0) return

  !locerrstat = pge_errstat_ok

  ! ==================================================
  ! Assign solar reference spectrum to local variables
  ! ==================================================
  npts            = refspecs_original(solar_idx)%nPoints

  num_alloc = max(npts, n_sensor_pts)
  allocate (locwvl(num_alloc), locspec(num_alloc), &
            specmod(num_alloc), tmpwav(num_alloc), &
            over(num_alloc), under(num_alloc), &
            resample(num_alloc), stat=locerrstat)
  if (locerrstat /= 0) then
    call tell_error (tell_malloc_error, "undersample: allocate failed", &
                     errstat)
    return
  endif

  locwvl (1:npts) = refspecs_original(solar_idx)%RefSpecWavs(1:npts)
  locspec(1:npts) = refspecs_original(solar_idx)%RefSpecData(1:npts)

  !IF ( yn_use_labslitfunc ) THEN
  !  CALL omi_slitfunc_convolve ( &
  !    xtrack_pix, npts, locwvl(1:npts), locspec(1:npts), specmod(1:npts), locerrstat )
  !ELSE
  !  CALL asymmetric_gaussian_sf ( &
  !    npts, hw1e, e_asym, locwvl(1:npts), locspec(1:npts), specmod(1:npts))
  !END IF
  !CALL error_check ( &
  !  locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
  !  modulename//f_sep//'Convolution', vb_lev_default, errstat )
  !IF ( locerrstat >= pge_errstat_error ) RETURN
  call slitfunction_convolve (npts, locwvl(1:npts), locspec(1:npts), &
       specmod(1:npts), xtrack_pix, [hw1e, e_asym, k], 3, errstat)
  if (errstat /= 0) return

  ! Phase1 calculation: Calculate spline derivatives for KPNO data
  !                     Calculate solar spectrum at OMI positions

  CALL interpolation (                                                                &
    npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, curr_wvl(1:n_sensor_pts), &
    resample(1:n_sensor_pts), 'endpoints', 0.0_r8, did_full_range, locerrstat )
  if (locerrstat /= 0) then
    call tell_error (tell_runtime_error, &
                     "undersample_spectrum (phase 1a): interpolation error", &
                     errstat)
    return
  endif
  !CALL error_check (                                                    &
  !  locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
  !  modulename//f_sep//'Phase 1a', vb_lev_default, errstat )
  !IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  IF ( .NOT. did_full_range ) then
    call tell_log (2, "undersample_spectrum (phase 1a): interpolating on wavelength subset")
    !CALL error_check (           &
    !  0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
    !  modulename//f_sep//'Phase 1a', vb_lev_develop, errstat )
  endif

  ! Calculate solar spectrum at OMI + phase positions, original and resampled.

  ! ------------------------------------------------------------------------------
  ! The original ("modified K.C.) scheme to compute the UNDERSPEC wavelength array
  ! ------------------------------------------------------------------------------
  ! ( assumes ABS(PHASE) < 1.0 )
  ! ----------------------------
  tmpwav(2:n_sensor_pts) = (1.0_r8-phase)*curr_wvl(1:n_sensor_pts-1) + phase*curr_wvl(2:n_sensor_pts)
  tmpwav(1)              = curr_wvl(1)
  tmpwav(n_sensor_pts)   = curr_wvl(n_sensor_pts)
  IF ( tmpwav(2) <= tmpwav(1) ) tmpwav(2) = (tmpwav(1)+tmpwav(3))/2.0_r8

  CALL interpolation (                                                              &
    npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, tmpwav(1:n_sensor_pts), &
    over(1:n_sensor_pts), 'endpoints', 0.0_r8, did_full_range, locerrstat )
  if (locerrstat /= 0) then
    call tell_error (tell_runtime_error, &
                     "undersample_spectrum (phase 1b): interpolation error", &
                     errstat)
    return
  endif
  !CALL error_check (                                                    &
  !  locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
  !  modulename//f_sep//'Phase 1b', vb_lev_default, errstat )
  !IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. did_full_range ) CALL error_check (                 &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE,       &
  !     modulename//f_sep//'Phase 1b', vb_lev_omidebug, errstat )

  CALL interpolation (                                                                 &
    n_sensor_pts, curr_wvl(1:n_sensor_pts), resample(1:n_sensor_pts), n_sensor_pts, &
    tmpwav(1:n_sensor_pts), under(1:n_sensor_pts), 'endpoints', 0.0_r8,             &
    did_full_range, locerrstat )
  if (locerrstat /= 0) then
    call tell_error (tell_runtime_error, &
                     "undersample_spectrum (phase 1c): interpolation error", &
                     errstat)
    return
  endif
  !CALL error_check (                                                    &
  !  locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
  !  modulename//f_sep//'Phase 1c', vb_lev_default, errstat )
  !IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. did_full_range ) CALL error_check (                 &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE,       &
  !     modulename//f_sep//'Phase 1c', vb_lev_omidebug, errstat )

  underspec(1:n_sensor_pts,1) = over(1:n_sensor_pts) - under(1:n_sensor_pts)
  resample (1:n_sensor_pts  ) = over(1:n_sensor_pts)

  ! ------------------------------------------------
  ! Save spectra to final arrays for Undersampling 1
  ! ------------------------------------------------
  refspecs_original(us1_idx)%nPoints                     = n_sensor_pts
  refspecs_original(us1_idx)%NormFactor                  = 1.0E+00_R8
  refspecs_original(us1_idx)%RefSpecWavs(1:n_sensor_pts) = tmpwav(1:n_sensor_pts)
  refspecs_original(us1_idx)%RefSpecData(1:n_sensor_pts) = underspec(1:n_sensor_pts,1)

  database(1:n_sensor_pts,us1_idx) = underspec (1:n_sensor_pts,1)

  ! ---------------------------------------------------------
  ! If we haven't selected Undersampling 2 then we return now
  ! ---------------------------------------------------------
  IF ( .NOT. have_undersampling(us2_idx) ) RETURN

  ! --------------------------------------------------------------------------------------
  ! Phase2 calculation: Calculate solar spectrum at OMI positions, original and resampled.
  ! --------------------------------------------------------------------------------------
  CALL interpolation (                                                                &
    npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, curr_wvl(1:n_sensor_pts), &
    over(1:n_sensor_pts), 'endpoints', 0.0_r8, did_full_range, locerrstat )
  if (locerrstat /= 0) then
    call tell_error (tell_runtime_error, &
                     "undersample_spectrum (phase 2a): interpolation error", &
                     errstat)
    return
  endif
  !CALL error_check ( &
  !  locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
  !  modulename//f_sep//'Phase 2a', vb_lev_default, errstat )
  !IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. did_full_range ) CALL error_check ( &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
  !     modulename//f_sep//'Phase 2a', vb_lev_omidebug, errstat )

  CALL interpolation ( &
    n_sensor_pts, tmpwav(1:n_sensor_pts), resample(1:n_sensor_pts), n_sensor_pts, &
    curr_wvl(1:n_sensor_pts), under(1:n_sensor_pts), 'endpoints', 0.0_r8,         &
    did_full_range, locerrstat )
  if (locerrstat /= 0) then
    call tell_error (tell_runtime_error, &
                     "undersample_spectrum (phase 2b): interpolation error", &
                     errstat)
    return
  endif
  !CALL error_check ( &
  !  locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
  !  modulename//f_sep//'Phase 2b', vb_lev_default, errstat )
  !IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. did_full_range ) CALL error_check ( &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
  !     modulename//f_sep//'Phase 2b', vb_lev_omidebug, errstat )

  ! ------------------------------------
  ! Compute final undersampling spectrum
  ! ------------------------------------
  underspec(1:n_sensor_pts,2) = over(1:n_sensor_pts) - under(1:n_sensor_pts)

  ! ------------------------------------------------
  ! Save spectra to final arrays for Undersampling 2
  ! ------------------------------------------------
  refspecs_original(us2_idx)%nPoints                     = n_sensor_pts
  refspecs_original(us2_idx)%NormFactor                  = 1.0E+00_R8
  refspecs_original(us2_idx)%RefSpecWavs(1:n_sensor_pts) = curr_wvl(1:n_sensor_pts)
  refspecs_original(us2_idx)%RefSpecData(1:n_sensor_pts) = underspec(1:n_sensor_pts,2)

  database(1:n_sensor_pts,us2_idx) = underspec (1:n_sensor_pts,2)

  RETURN
END SUBROUTINE undersample_spectrum

!UNUSED! SUBROUTINE undersample_new ( xtrack_pix, n_sensor_pts, curr_wvl, n_solar_pts, solar_wvl, &
!UNUSED!     hw1e, e_asym, errstat )
!UNUSED!
!UNUSED!   !     Convolves input spectrum with Gaussian slit function of specified
!UNUSED!   !     HW1e, and samples at a particular input phase to give the OMI
!UNUSED!   !     undersampling spectrum. This version calculates both phases of the
!UNUSED!   !     undersampling spectrum, phase1 - i.e., underspec (i, 1) - being the
!UNUSED!   !     more common in OMI spectra.
!UNUSED!
!UNUSED!   USE OMSAO_precision_module
!UNUSED!   USE OMSAO_indices_module,    ONLY: solar_idx, us1_idx, us2_idx
!UNUSED!   USE OMSAO_parameters_module, ONLY: max_spec_pts
!UNUSED!   USE OMSAO_variables_module,  ONLY: &
!UNUSED!     refspecs_original, database
!UNUSED!   USE OMSAO_slitfunction_module
!UNUSED!   USE OMSAO_errstat_module
!UNUSED!   USE sao_pge_utils, ONLY: interpolation
!UNUSED!
!UNUSED!   IMPLICIT NONE
!UNUSED!
!UNUSED!   ! ---------------
!UNUSED!   ! Input variables
!UNUSED!   ! ---------------
!UNUSED!   INTEGER (KIND=i4),                           INTENT (IN) :: n_sensor_pts, xtrack_pix, &
!UNUSED!     n_solar_pts
!UNUSED!   REAL    (KIND=r8),                           INTENT (IN) :: hw1e, e_asym
!UNUSED!   REAL    (KIND=r8), DIMENSION (n_sensor_pts), INTENT (IN) :: curr_wvl
!UNUSED!   REAL    (KIND=r8), DIMENSION (n_solar_pts), INTENT (IN)  :: solar_wvl
!UNUSED!
!UNUSED!   ! ---------------
!UNUSED!   ! Output variable
!UNUSED!   ! ---------------
!UNUSED!   INTEGER (KIND=i4), INTENT (INOUT) :: errstat
!UNUSED!
!UNUSED!   ! ---------------
!UNUSED!   ! Local variables
!UNUSED!   ! ---------------
!UNUSED!   LOGICAL                                     :: did_full_range
!UNUSED!   REAL (KIND=r8), DIMENSION (n_sensor_pts,2)  :: underspec
!UNUSED!   REAL (KIND=r8), DIMENSION (max_spec_pts)    :: &
!UNUSED!     locwvl, locspec, specmod, tmpwav, over, under, resample
!UNUSED!
!UNUSED!   INTEGER (KIND=i4) :: npts, locerrstat
!UNUSED!
!UNUSED!   ! ------------------------------
!UNUSED!   ! Name of this subroutine/module
!UNUSED!   ! ------------------------------
!UNUSED!   CHARACTER (LEN=23), PARAMETER :: modulename = 'undersample_new'
!UNUSED!
!UNUSED!   locerrstat = pge_errstat_ok
!UNUSED!
!UNUSED!   ! ==================================================
!UNUSED!   ! Assign solar reference spectrum to local variables
!UNUSED!   ! ==================================================
!UNUSED!   npts            = refspecs_original(solar_idx)%nPoints
!UNUSED!   locwvl (1:npts) = refspecs_original(solar_idx)%RefSpecWavs(1:npts)
!UNUSED!   locspec(1:npts) = refspecs_original(solar_idx)%RefSpecData(1:npts)
!UNUSED!
!UNUSED!   tmpwav(2:n_solar_pts) = solar_wvl
!UNUSED!   tmpwav(1)             = curr_wvl(1)
!UNUSED!   tmpwav(n_solar_pts)   = curr_wvl(n_sensor_pts)
!UNUSED!   IF ( tmpwav(2) <= tmpwav(1) ) tmpwav(2) = (tmpwav(1)+tmpwav(3))/2.0_r8
!UNUSED!   IF ( tmpwav(n_solar_pts-1) >= tmpwav(n_solar_pts) ) tmpwav(n_solar_pts-1) = &
!UNUSED!     (tmpwav(n_solar_pts)+tmpwav(n_solar_pts-2))/2.0_r8
!UNUSED!
!UNUSED!   IF ( yn_use_labslitfunc ) THEN
!UNUSED!     CALL omi_slitfunc_convolve ( &
!UNUSED!       xtrack_pix, npts, locwvl(1:npts), locspec(1:npts), specmod(1:npts), locerrstat )
!UNUSED!   ELSE
!UNUSED!     CALL asymmetric_gaussian_sf ( &
!UNUSED!       npts, hw1e, e_asym, locwvl(1:npts), locspec(1:npts), specmod(1:npts))
!UNUSED!   END IF
!UNUSED!   CALL error_check ( &
!UNUSED!     locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
!UNUSED!     modulename//f_sep//'Convolution', vb_lev_default, errstat )
!UNUSED!   IF ( locerrstat >= pge_errstat_error ) RETURN
!UNUSED!
!UNUSED!   ! Phase1 calculation: Calculate spline derivatives for KPNO data
!UNUSED!   !                     Calculate solar spectrum at OMI radiance positions
!UNUSED!
!UNUSED!   CALL interpolation (                                                                &
!UNUSED!     npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, curr_wvl(1:n_sensor_pts), &
!UNUSED!     over(1:n_sensor_pts), 'endpoints', 0.0_r8, did_full_range, locerrstat )
!UNUSED!   CALL error_check (                                                    &
!UNUSED!     locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
!UNUSED!     modulename//f_sep//'Phase 1a', vb_lev_default, errstat )
!UNUSED!   IF ( locerrstat >= pge_errstat_error ) RETURN
!UNUSED!
!UNUSED!   ! -------------------------------------------------------------
!UNUSED!   ! Issue a warning if we don't have the full interpolation range
!UNUSED!   ! -------------------------------------------------------------
!UNUSED!   IF ( .NOT. did_full_range ) CALL error_check (           &
!UNUSED!     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
!UNUSED!     modulename//f_sep//'Phase 1a', vb_lev_develop, errstat )
!UNUSED!
!UNUSED!   ! Convolved High resolution solar to OMI solar grid
!UNUSED!   CALL interpolation (                                                               &
!UNUSED!     npts, locwvl(1:npts), specmod(1:npts), n_solar_pts, tmpwav(1:n_solar_pts), &
!UNUSED!     resample(1:n_solar_pts), 'endpoints', 0.0_r8, did_full_range, locerrstat )
!UNUSED!   CALL error_check (                                                    &
!UNUSED!     locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
!UNUSED!     modulename//f_sep//'Phase 1b', vb_lev_default, errstat )
!UNUSED!   IF ( locerrstat >= pge_errstat_error ) RETURN
!UNUSED!
!UNUSED!   ! -------------------------------------------------------------
!UNUSED!   ! Issue a warning if we don't have the full interpolation range
!UNUSED!   ! -------------------------------------------------------------
!UNUSED!   !IF ( .NOT. did_full_range ) CALL error_check (                 &
!UNUSED!   !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE,       &
!UNUSED!   !     modulename//f_sep//'Phase 1b', vb_lev_omidebug, errstat )
!UNUSED!
!UNUSED!   ! Undersample solar to radiance grid
!UNUSED!   CALL interpolation (                                                                 &
!UNUSED!     n_solar_pts, tmpwav(1:n_solar_pts), resample(1:n_solar_pts), n_sensor_pts, &
!UNUSED!     curr_wvl(1:n_sensor_pts), under(1:n_sensor_pts), 'endpoints', 0.0_r8,             &
!UNUSED!     did_full_range, locerrstat )
!UNUSED!   CALL error_check (                                                    &
!UNUSED!     locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
!UNUSED!     modulename//f_sep//'Phase 1c', vb_lev_default, errstat )
!UNUSED!   IF ( locerrstat >= pge_errstat_error ) RETURN
!UNUSED!
!UNUSED!   ! Calculate undersample spectrum
!UNUSED!   underspec(1:n_sensor_pts,1) = over(1:n_sensor_pts) - under(1:n_sensor_pts)
!UNUSED!
!UNUSED!   ! ------------------------------------------------
!UNUSED!   ! Save spectra to final arrays for Undersampling 1
!UNUSED!   ! ------------------------------------------------
!UNUSED!   refspecs_original(us1_idx)%nPoints                     = n_sensor_pts
!UNUSED!   refspecs_original(us1_idx)%NormFactor                  = 1.0E+00_R8
!UNUSED!   refspecs_original(us1_idx)%RefSpecWavs(1:n_sensor_pts) = curr_wvl(1:n_sensor_pts)
!UNUSED!   refspecs_original(us1_idx)%RefSpecData(1:n_sensor_pts) = underspec(1:n_sensor_pts,1)
!UNUSED!   refspecs_original(us2_idx)%nPoints                     = n_sensor_pts
!UNUSED!   refspecs_original(us2_idx)%NormFactor                  = 1.0E+00_R8
!UNUSED!   refspecs_original(us2_idx)%RefSpecWavs(1:n_sensor_pts) = curr_wvl(1:n_sensor_pts)
!UNUSED!   refspecs_original(us2_idx)%RefSpecData(1:n_sensor_pts) = underspec(1:n_sensor_pts,1)
!UNUSED!
!UNUSED!   ! Save undersample spectrum to database, for compliance with previous versions
!UNUSED!   ! it is saved to both under sample spectra but only one of them needs to be used.
!UNUSED!   database(1:n_sensor_pts,us1_idx) = underspec (1:n_sensor_pts,1)
!UNUSED!   database(1:n_sensor_pts,us2_idx) = underspec (1:n_sensor_pts,1)
!UNUSED!
!UNUSED!   RETURN
!UNUSED! END SUBROUTINE undersample_new

END MODULE
