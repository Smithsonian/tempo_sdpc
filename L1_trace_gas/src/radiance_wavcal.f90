MODULE radiance_wavcal

  private
  public radiance_wavecal

CONTAINS
SUBROUTINE radiance_wavecal ( &
    ipix, n_rad_wvl, adj_wvls, adj_spec, adj_wgts, adj_resid, &
    n_fitres_loop, fitres_range, &
    radcal_exval, radcal_itnum, chisquav, &
    is_bad_pixel, errstat )

  USE OMSAO_precision_module,     ONLY: i2, i4, r8
  USE OMSAO_parameters_module,    ONLY: &
    i2_missval, i4_missval, r8_missval
  USE OMSAO_indices_module,       ONLY: &
    max_calfit_idx, shi_idx, squ_idx, &
    hwe_idx, asy_idx, sgk_idx
  USE OMSAO_variables_module,   ONLY: &
    fitvar_cal, fitvar_rad_init, fitvar_cal_saved, &
    lo_radbnd, up_radbnd, max_itnum_sol, &
    Slit_Half_Width_1e, Slit_Asym_Factor, Slit_Shape_Factor
  USE commonmode, ONLY: compute_common_mode
  use wavecal
  use optimizer_interface_module, only: opt_convergence_good
  IMPLICIT NONE

  ! ----------------
  ! Input parameters
  ! ----------------
  INTEGER (KIND=i4), INTENT (IN) :: ipix, n_fitres_loop, n_rad_wvl, &
       fitres_range

  ! -------------------
  ! Modified parameters
  ! -------------------
  LOGICAL,           INTENT (OUT)   :: is_bad_pixel
  INTEGER (KIND=i2), INTENT (OUT)   :: radcal_itnum
  INTEGER (KIND=i4), INTENT (OUT)   :: radcal_exval
  REAL    (KIND=r8), INTENT (OUT)   :: chisquav
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
  real (kind=r8), dimension(n_rad_wvl), intent(inout) :: adj_wvls, adj_spec, &
       adj_wgts, adj_resid

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: locitnum
  real (kind=r8) :: sol_wav_avg

  integer (kind=i4) :: i, num_fitvar
  real (kind=r8), dimension(max_calfit_idx) :: fitvar_saved, fitvar, lobnd, upbnd

  if (errstat /= 0) return

  is_bad_pixel = .FALSE.

  ! Select and wavelength calibrate radiance spectrum
  radcal_exval = i4_missval
  radcal_itnum = i2_missval
  chisquav     = r8_missval

  ! ------------------------------------------
  ! Update wavelegths for common mode spectrum
  ! (the .TRUE. in the call below selects the
  !  "wavelength update only" branch)
  ! ------------------------------------------
  CALL compute_common_mode ( &
    .TRUE., ipix, n_rad_wvl, adj_wvls(1:n_rad_wvl), adj_spec(1:n_rad_wvl), errstat)
  if (errstat /= 0) return

  ! -------------------------------------------------------------
  ! Initialize the fitting variables. FITVAR_CAL_SAVED has been
  ! set to the initial values in the calling routine. outside the
  ! pixel loop. Here we use FITVAR_CAL_SAVED, which will be
  ! updated with current values from the previous fit if that fit
  ! has gone well.
  ! -------------------------------------------------------------
  fitvar_cal(1:max_calfit_idx) = fitvar_cal_saved(1:max_calfit_idx)

  ! -------------------------------------------------------------------------
  ! Keep the slit function variables from solar fit fixed. Remember to reduce
  ! the number of solar fitting variables if previously varied.
  ! -------------------------------------------------------------------------
  fitvar_cal(hwe_idx) = Slit_Half_Width_1e
  lo_radbnd (hwe_idx) = Slit_Half_Width_1e
  up_radbnd (hwe_idx) = Slit_Half_Width_1e
  fitvar_cal(asy_idx) = Slit_Asym_Factor
  lo_radbnd (asy_idx) = Slit_Asym_Factor
  up_radbnd (asy_idx) = Slit_Asym_Factor
  fitvar_cal(sgk_idx) = Slit_Shape_Factor
  lo_radbnd (sgk_idx) = Slit_Shape_Factor
  up_radbnd (sgk_idx) = Slit_Shape_Factor

  ! -----------------------------------------------------
  ! Assign the solar average wavelength - the wavelength
  ! calibration will not converge without it!
  ! -----------------------------------------------------
  sol_wav_avg = SUM (adj_wvls(1:n_rad_wvl) ) / REAL(n_rad_wvl,KIND=r8)

  ! on input, set locitnum to the max per fit, upon output it will be set to the total
  locitnum = max_itnum_sol
  call wavecal_fit (adj_wvls, adj_spec, adj_wgts, adj_resid, n_rad_wvl, sol_wav_avg, &
                    fitvar_cal, lo_radbnd, up_radbnd, max_calfit_idx, &
                    n_fitres_loop, real(fitres_range, kind=r8), &
                    is_bad_pixel, locitnum, chisquav, radcal_exval, errstat)
  if (errstat /= 0) return
  radcal_itnum = int(locitnum, kind=i2)

  ! pack paramters to update fitvar_cal
  num_fitvar = 0
  do i = 1, max_calfit_idx
    if (lo_radbnd(i) >= up_radbnd(i)) cycle
    num_fitvar = num_fitvar + 1
    fitvar_saved(num_fitvar) = fitvar_cal_saved(i)
    fitvar(num_fitvar) = fitvar_cal(i)
    lobnd (num_fitvar) = lo_radbnd(i)
    upbnd (num_fitvar) = up_radbnd(i)
  enddo

  ! ------------------------------------------------------------------
  ! The following assignment makes sense only because FITVAR_CAL is
  ! updated with FITVAR (using the proper mask) in SPECTRUM_SOLAR.
  ! ------------------------------------------------------------------
  IF ( radcal_exval == opt_convergence_good .and. &
     ( any(fitvar(1:num_fitvar) .ne. fitvar_saved(1:num_fitvar)) .and. &
       all(fitvar(1:num_fitvar) .ne. lobnd(1:num_fitvar)) .and. &
       all(fitvar(1:num_fitvar) .ne. upbnd(1:num_fitvar))        )   ) THEN
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
  adj_wvls(1:n_rad_wvl) = ( &
    adj_wvls(1:n_rad_wvl) - fitvar_cal_saved(shi_idx) &
    + sol_wav_avg * fitvar_cal_saved(squ_idx)) / &
    (1.0_r8 + fitvar_cal_saved(squ_idx))

  RETURN

END SUBROUTINE radiance_wavecal

END MODULE

