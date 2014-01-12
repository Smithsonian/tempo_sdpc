MODULE radiance_wavcal
  use optimizer_interface_module
  use elsunc_interface_module
  use spectra, only: spectrum_solar

  public radiance_wavecal
  private solar_residuals

CONTAINS
SUBROUTINE radiance_wavecal (                            &
    ipix, n_fitres_loop, fitres_range, n_rad_wvl,        &
    curr_rad_spec, radcal_exval, radcal_itnum, chisquav, &
    is_bad_pixel, errstat )

  USE OMSAO_precision_module,     ONLY: i2, i4, r8
  USE OMSAO_parameters_module,    ONLY: &
    i2_missval, i4_missval, r8_missval, downweight
  USE OMSAO_elsunc_fitting_module, ONLY: ELSUNC_LESS_IS_NOISE
  USE OMSAO_indices_module,       ONLY: &
    max_calfit_idx, shi_idx, squ_idx, wvl_idx, spc_idx, sig_idx, ccd_idx, &
    hwe_idx, asy_idx
  USE OMSAO_variables_module,   ONLY: &
    fitwavs, fitweights, currspec,                                 &
    fitvar_cal, fitvar_rad_init, fitvar_sol_init, fitvar_cal_saved,&
    mask_fitvar_cal, n_fitvar_cal,           &
    lo_radbnd, up_radbnd, lobnd, upbnd,  &
    max_itnum_sol, Slit_Half_Width_1e, Slit_Asym_Factor, yn_newshift, sol_wav_avg
  USE OMSAO_errstat_module
  USE omi_pge_fitting_aux, ONLY: compute_common_mode
  IMPLICIT NONE

  ! ----------------
  ! Input parameters
  ! ----------------
  INTEGER (KIND=i4), INTENT (IN) :: ipix, n_fitres_loop, n_rad_wvl, fitres_range

  ! -------------------
  ! Modified parameters
  ! -------------------
  LOGICAL,                                          INTENT (OUT)   :: is_bad_pixel
  INTEGER (KIND=i2),                                INTENT (OUT)   :: radcal_itnum
  INTEGER (KIND=i4),                                INTENT (OUT)   :: radcal_exval
  REAL    (KIND=r8),                                INTENT (OUT)   :: chisquav
  INTEGER (KIND=i4),                                INTENT (INOUT) :: errstat
  REAL    (KIND=r8), DIMENSION (ccd_idx,n_rad_wvl), INTENT (INOUT) :: curr_rad_spec

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: i, locerrstat, locitnum, n_nozero_wgt
  REAL    (KIND=r8)  :: mean, sdev, loclim
  REAL    (KIND=r8), DIMENSION (n_rad_wvl)         :: fitres, fitspec
  REAL    (KIND=r8), DIMENSION (max_calfit_idx)    :: fitvar

  type(optimizer_type) :: opt
  integer (kind=i4) :: return_status

  is_bad_pixel = .FALSE.

  ! Select and wavelength calibrate radiance spectrum

  locerrstat = pge_errstat_ok

  radcal_exval = i4_missval
  radcal_itnum = i2_missval
  chisquav     = r8_missval

  ! ----------
  ! ELSUNC fit
  ! ----------
  fitwavs   (1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)
  currspec  (1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl)
  fitweights(1:n_rad_wvl) = curr_rad_spec(sig_idx, 1:n_rad_wvl)

  ! ------------------------------------------
  ! Update wavelegths for common mode spectrum
  ! (the .TRUE. in the call below selects the
  !  "wavelength update only" branch)
  ! ------------------------------------------
  CALL compute_common_mode ( &
    .TRUE., ipix, n_rad_wvl, fitwavs(1:n_rad_wvl), currspec(1:n_rad_wvl))

  ! -------------------------------------------------------------
  ! Initialize the fitting variables. FITVAR_CAL_SAVED has been
  ! set to the initial values in the calling routine. outside the
  ! pixel loop. Here we use FITVAR_CAL_SAVED, which will be
  ! updated with current values from the previous fit if that fit
  ! has gone well.
  ! -------------------------------------------------------------
  !fitvar_cal(1:max_calfit_idx) = fitvar_cal_saved(1:max_calfit_idx)
  !fitvar_cal(1:max_calfit_idx) = fitvar_rad_init(1:max_calfit_idx)
  fitvar_cal(1:max_calfit_idx) = fitvar_sol_init(1:max_calfit_idx)

  ! -------------------------------------------------------------------------
  ! Keep the slit function variables from solar fit fixed. Remember to reduce
  ! the number of solar fitting variables if previously varied.
  ! -------------------------------------------------------------------------
  fitvar = 0.0_r8 ; lobnd = 0.0_r8 ; upbnd = 0.0_r8 ; n_fitvar_cal = 0
  DO i = 1, max_calfit_idx
    SELECT CASE ( i )
    CASE ( hwe_idx )
      fitvar_cal(hwe_idx) = Slit_Half_Width_1e
      lo_radbnd (hwe_idx) = Slit_Half_Width_1e
      up_radbnd (hwe_idx) = Slit_Half_Width_1e
    CASE ( asy_idx )
      fitvar_cal(asy_idx) = Slit_Asym_Factor
      lo_radbnd (asy_idx) = Slit_Asym_Factor
      up_radbnd (asy_idx) = Slit_Asym_Factor
    CASE DEFAULT
      IF (lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_cal  = n_fitvar_cal + 1
        mask_fitvar_cal(n_fitvar_cal) = i
        fitvar(n_fitvar_cal) = fitvar_cal(i)
        lobnd (n_fitvar_cal) = lo_radbnd(i)
        upbnd (n_fitvar_cal) = up_radbnd(i)
      END IF
    END SELECT
  END DO

  ! --------------------------------------------------------------------
  ! Check whether we enough spectral points to carry out the fitting. If
  ! not, call it a bad pixel and return.
  ! --------------------------------------------------------------------
  IF ( n_fitvar_cal >= n_rad_wvl ) THEN
    is_bad_pixel = .TRUE.  ;  RETURN
  END IF

  loclim = 0.0_r8
  radcal_itnum = 0
  i = 0

  ! ---------------------------------------------------------------------
  ! Attempt to standardize the re-iteration with spectral points excluded
  ! that have fitting residuals larger than a pre-set window. Needs more
  ! thinking before it can replace a simple window determined empirically
  ! from fitting lots of spectra.
  ! ---------------------------------------------------------------------

  fit_loop: do
    call optimizer_open (opt, elsunc_optimizer, solar_residuals, n_fitvar_cal, return_status, &
                         param_min = lobnd(1:n_fitvar_cal), &
                         param_max = upbnd(1:n_fitvar_cal), &
                         param_mask = mask_fitvar_cal(1:n_fitvar_cal), &
                         max_num_iterations = max_itnum_sol)
    if (return_status < 0) then
      write(*,*)' radiance_wavecal: optimizer_open failed '
      stop  ! FIXME!!!
    endif

    call opt%optimize (opt, fitvar(1:n_fitvar_cal), n_fitvar_cal, &
                       fitres(1:n_rad_wvl), n_rad_wvl, return_status)
    locitnum = opt%num_iterations
    radcal_exval = return_status

    call optimizer_close (opt, return_status)
    if (return_status < 0) then
      write(*,*)'radiance_wavecal: optimizer_close failed'
      stop  ! FIXME!
    endif

    CALL spectrum_solar ( &
      n_rad_wvl, sol_wav_avg, fitwavs(1:n_rad_wvl), fitspec(1:n_rad_wvl), &
      fitvar_cal, max_calfit_idx)

    n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_rad_wvl)) ) )
    IF (n_nozero_wgt > 0) THEN
      chisquav = SUM  ( fitres(1:n_rad_wvl)**2 )
    ELSE
      chisquav = r8_missval
    END IF

    radcal_itnum = radcal_itnum + INT ( locitnum, KIND=i2 )
    i = i + 1

    if (1 < i .and. i <= n_fitres_loop) then
      if (maxval(abs(fitres(1:n_rad_wvl))) <= loclim) exit fit_loop
    else if (i == 1) then
      mean    = SUM  ( fitres(1:n_rad_wvl) )                 / REAL(n_nozero_wgt,   KIND=r8)
      sdev    = SQRT ( SUM ( (fitres(1:n_rad_wvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
      loclim  = REAL (fitres_range, KIND=r8)*sdev
      if (.not. ((n_fitres_loop > 0) &
                 .and. (loclim > 0.0_r8) &
                 .and. (MAXVAL(ABS(fitres(1:n_rad_wvl))) >= loclim) &
                 .and. (n_nozero_wgt > n_fitvar_cal))) exit fit_loop
    else
      exit fit_loop ! i > n_fitres_loop
    endif

    WHERE ( ABS(fitres(1:n_rad_wvl)) > loclim )
      fitweights(1:n_rad_wvl) = downweight
    END WHERE

  enddo fit_loop

  ! ------------------------------------------------------------------
  ! The following assignment makes sense only because FITVAR_CAL is
  ! updated with FITVAR (using the proper mask) in SPECTRUM_SOLAR.
  ! ------------------------------------------------------------------
  IF ( radcal_exval >= INT(ELSUNC_LESS_IS_NOISE, KIND=i4) ) THEN
    fitvar_cal_saved(1:max_calfit_idx) = fitvar_cal(1:max_calfit_idx)
  ELSE
    fitvar_cal_saved(1:max_calfit_idx) = fitvar_rad_init(1:max_calfit_idx)
  END IF

  ! -----------------------------------------------------------
  ! Reality check: set SQUEEZE to 1.0 to avoid division by Zero
  ! -----------------------------------------------------------
  IF ( fitvar_cal(squ_idx) == -1.0_r8 ) fitvar_cal(squ_idx) = 0.0_r8

  ! ---------------------
  ! Perform Shift&Squueze
  ! ---------------------
  ! gga to include Xiong comments
  !  write(*,*) curr_rad_spec(wvl_idx,1:n_rad_wvl), fitvar_cal(shi_idx), fitvar_cal(squ_idx)
  IF (yn_newshift .EQV. .true.) THEN
    curr_rad_spec(wvl_idx,1:n_rad_wvl) = ( &
      curr_rad_spec(wvl_idx,1:n_rad_wvl) - fitvar_cal(shi_idx) + sol_wav_avg * fitvar_cal(squ_idx)) / &
      (1.0_r8 + fitvar_cal(squ_idx))
  ELSE
    curr_rad_spec(wvl_idx,1:n_rad_wvl) = ( &
      curr_rad_spec(wvl_idx,1:n_rad_wvl) - fitvar_cal(shi_idx) ) / (1.0_r8 + fitvar_cal(squ_idx))
  END IF
  !  write(*,*) curr_rad_spec(wvl_idx,1:n_rad_wvl)

  ! -----------------------------------------------------------------
  ! We haven't implemented any error checks in this subroutine, hence
  ! it doesn't make sense yet to update the error variable. We'll do
  ! it anyway to make ourselves feel better.
  ! -----------------------------------------------------------------
  IF ( locerrstat /= pge_errstat_ok ) errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE radiance_wavecal

  subroutine solar_residuals (this_optimizer, params, num_params, residuals, num_residuals, return_status)
    USE OMSAO_indices_module, ONLY: max_calfit_idx
    USE OMSAO_variables_module,  ONLY : &
      fitwavs, fitweights, currspec, sol_wav_avg, fitvar_cal
    implicit none
    type(optimizer_type) :: this_optimizer
    real (kind=r8), dimension (:), intent(in) :: params
    real (kind=r8), dimension (:), intent(out) :: residuals
    integer (kind=i4), intent(in) :: num_params, num_residuals
    integer (kind=i4), intent(out) :: return_status
    ! local variables
    integer (kind=i4) :: i, idx

    ! unpack fit parameters
    DO i = 1, num_params
      idx = this_optimizer%param_mask(i)
      fitvar_cal(idx) = params(i)
    END DO

    CALL spectrum_solar (num_residuals, sol_wav_avg, fitwavs(1:num_residuals), &
                         residuals(1:num_residuals), fitvar_cal, max_calfit_idx)

    residuals = (currspec(1:num_residuals) - residuals(1:num_residuals)) * fitweights(1:num_residuals)
    return_status = 0

  end subroutine solar_residuals

END MODULE

