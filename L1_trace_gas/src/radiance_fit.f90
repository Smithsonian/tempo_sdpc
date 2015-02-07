MODULE radiance_fit
  USE OMSAO_precision_module, only : i4, r8
  use OMSAO_indices_module, only : max_rs_idx, n_max_fitpars
  use OMSAO_parameters_module, only : nwavel_max
  use optimizer_interface_module
  use tell_module

  private
  public fit_radiance

  interface
    subroutine earthshine_spectrum_interface (npts, avg_wavl, wavelengths, spectrum, params)
      import i4, r8
      implicit none
      integer (kind=i4), intent(in) :: npts
      real (kind=r8), intent(in) :: avg_wavl
      real (kind=r8), dimension(:), intent(inout) :: params
      real (kind=r8), dimension(npts), intent(in) :: wavelengths
      real (kind=r8), dimension(npts), intent(out) :: spectrum
    end subroutine earthshine_spectrum_interface
  end interface

  type spectrum_type
    real (kind=r8), dimension(nwavel_max) :: spec, wavs, weights
    real (kind=r8) :: wav_avg
  end type spectrum_type

  type (spectrum_type) :: Spec
  logical, dimension(n_max_fitpars) :: param_frozen_at_zero
  logical, dimension(max_rs_idx)    :: database_j_is_zero

CONTAINS

  SUBROUTINE fit_radiance ( &
      pge_idx, ipix, num_fitres_loops, fitres_range, &
      n_rad_wvl_loc, adj_wvls, adj_spec, adj_wgts, &
      fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,   &
      prefit, o3fit_cols, o3fit_dcols,                              &
      target_var, allfit, allerr, corrmat, is_bad_pixel, fitspc_out, &
      yn_reference_fit, errstat)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,      ONLY: &
      solar_idx, n_max_fitpars, &
      max_calfit_idx, max_rs_idx, o3_t1_idx, o3_t3_idx,  pge_o3_idx
    USE OMSAO_parameters_module,   ONLY: &
      r8_missval, i2_missval, downweight
    USE OMSAO_variables_module,    ONLY:                                    &
      n_fincol_idx, fincol_idx, database, &
      fitvar_rad, n_fitvar_rad,      &
      lo_radbnd, up_radbnd, &
      fit_winwav_idx, mask_fitvar_rad, max_itnum_rad, refspecs_original, &
      all_radfit_idx, &
      n_rad_wvl_max, fitvar_rad_init, fitvar_rad_saved, &
      tol, epsrel, epsabs, epsx
    use ctrlvars, only: yn_smooth, yn_doas, yn_o3amf_cor
    USE OMSAO_prefitcol_module, ONLY:  prefit_type, apply_prefit_values_and_bounds, n_prefit_vars
    USE commonmode, ONLY: compute_common_mode
    USE subtract_cubic, ONLY: cubic_subtract_meas
    use arrayutils, only: array_smooth
    IMPLICIT NONE

    ! *******************************************************************
    ! CAREFUL: Assumes that radiance and solar wavelength arrays have the
    ! same number of points. That must not be the case if we read in a
    ! general EL1 file. Examine and adjust! (tpk, note to himself)
    ! *******************************************************************

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: &
      pge_idx, ipix, n_rad_wvl_loc, num_fitres_loops, fitres_range

    type (prefit_type), intent(in) :: prefit
    logical, intent(in) :: yn_reference_fit

    ! -----------------------------
    ! (Possibly) Modified Variables
    ! -----------------------------
    real (kind=r8), dimension (n_rad_wvl_loc), intent(inout) :: &
      adj_wvls, adj_spec, adj_wgts
    REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx), INTENT (INOUT) &
      :: o3fit_cols, o3fit_dcols
    integer, intent(inout) :: errstat

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,                                     INTENT (OUT) :: is_bad_pixel
    INTEGER (KIND=i4),                           INTENT (OUT) :: radfit_exval, radfit_itnum
    REAL    (KIND=r8),                           INTENT (OUT) :: fitcol, rms, dfitcol, chisquav
    REAL    (KIND=r8), DIMENSION (n_fitvar_rad), INTENT (OUT) :: allfit, allerr, corrmat
    REAL    (KIND=r8), DIMENSION (n_fincol_idx), INTENT (OUT) :: target_var

    ! CCM Also return fitted spectrum
    REAL (KIND=r8), DIMENSION (n_rad_wvl_max), INTENT (OUT) :: fitspc_out

    ! ---------------
    ! Local variables
    ! ---------------
    procedure(earthshine_spectrum_interface), pointer :: earthshine_spectrum => null()
    INTEGER (KIND=i4) :: i, j, idx, j1, j2, k1, k2, l, ll_rad, lu_rad, indx
    INTEGER (KIND=i4) :: n_fitwav_rad, locitnum, n_nozero_wgt
    REAL    (KIND=r8)                                         :: asum, ssum
    REAL    (KIND=r8)                                         :: mean, sdev, loclim, normfac, mfac
    REAL    (KIND=r8), DIMENSION (n_rad_wvl_loc)              :: fitres, fitspec
    REAL    (KIND=r8), DIMENSION (n_max_fitpars)              :: fitvar, lobnd, upbnd
    REAL    (KIND=r8), DIMENSION (n_rad_wvl_loc, n_fitvar_rad) :: covar_matrix

    REAL    (KIND=r8) :: fitcol_saved, covar_xx, pm_one

    type(optimizer_type) :: opt
    integer (kind=i4) :: return_status
    character (len=128) :: log_msg

    SAVE fitcol_saved

    if (errstat < 0) return

    is_bad_pixel = .FALSE.

    radfit_exval = INT(i2_missval, KIND=i4)
    radfit_itnum = INT(i2_missval, KIND=i4)
    chisquav     = r8_missval

    IF ( pge_idx == pge_o3_idx ) THEN
      o3fit_cols = r8_missval  ;  o3fit_dcols = r8_missval
    END IF

    ! ============================================================
    ! Assign LL_RAD, LU_RAD, and weights for each earthshine radiance
    ! Assign LL_RAD, LU_RAD for earthshine radiance
    ! ============================================================
    ll_rad = fit_winwav_idx(2)  ;  lu_rad = fit_winwav_idx(3)

    Spec%wavs(1:n_rad_wvl_loc) = adj_wvls(1:n_rad_wvl_loc)
    Spec%spec(1:n_rad_wvl_loc) = adj_spec(1:n_rad_wvl_loc)
    Spec%weights(1:n_rad_wvl_loc) = adj_wgts(1:n_rad_wvl_loc)

    ! ---------------------------------------------------------------
    ! High pass filtering for DOAS. First, take log (rad/irrad), then
    ! filter by subtracting a cubic, then re-add the log (irradiance).
    ! This way we are fitting the log (rad) with the proper filtering
    ! having been done. The spectrum subroutine will start by
    ! subtracting the irradiance spectrum, where the proper shifting
    ! and squeezing can take place. In this high-pass filtering, we
    ! are ignoring the small extra wavelength calibration change, as
    ! for Ring effect, above.
    ! ---------------------------------------------------------------
    IF ( yn_doas ) THEN
      Spec%spec(1:n_rad_wvl_loc) = LOG ( Spec%spec(1:n_rad_wvl_loc) / database(1:n_rad_wvl_loc,solar_idx) )
      CALL cubic_subtract_meas (Spec%wavs(1:n_rad_wvl_loc), n_rad_wvl_loc, Spec%spec(1:n_rad_wvl_loc), ll_rad, lu_rad, errstat)
      if (errstat < 0) return
      Spec%spec(1:n_rad_wvl_loc) = Spec%spec(1:n_rad_wvl_loc) + LOG ( database(1:n_rad_wvl_loc,solar_idx) )
    END IF

    ! --------------------------------------------------------------------
    ! Apply smoothing (1/16,1/4,3/8,1/4,1/16); 2/98 uhe/ife recommendation
    ! --------------------------------------------------------------------
    if ( yn_smooth ) THEN
      call array_smooth (spec%spec, n_rad_wvl_loc, errstat)
      if (errstat < 0) return
    endif

    ! ---------------------------------------
    ! Compute average of radiance wavelengths
    ! ---------------------------------------
    asum = SUM ( Spec%wavs(1:n_rad_wvl_loc) * ( Spec%weights(1:n_rad_wvl_loc)*Spec%weights(1:n_rad_wvl_loc) ) )
    ssum = SUM (             1.0_r8   * ( Spec%weights(1:n_rad_wvl_loc)*Spec%weights(1:n_rad_wvl_loc) ) )
    Spec%wav_avg = asum / ssum

    radfit_exval = 0

    ! --------------------------------------------------------------------
    ! Initialize the fitting variables with the initial guess. Rather than
    ! subjecting ourselves to the vagarities of a wrong convergence, we
    ! bite the computationally more expensive bullet of starting from
    ! scratch each and every time.
    ! --------------------------------------------------------------------

    ! -----------------------------------------------------------
    ! Initialize the fitting variables. FITVAR_RAD_SAVED has been
    ! set to the initial values in the calling routine outside the
    ! pixel loop. Here we use FITVAR_RAD_SAVED, which will be
    ! updated with current values from the previous fit if that
    ! fit has gone well.
    ! -----------------------------------------------------------
    fitvar_rad(1:n_max_fitpars) = fitvar_rad_saved(1:n_max_fitpars)

    ! ---------------------------------------------------
    ! For HCHO we check whether we will be doing pre-fits
    ! ---------------------------------------------------------------------------
    ! Note that the conversion to slant columns (BrO only) and the multiplication
    ! with the cross section normalization has already been done at the time when
    ! the pre-fitted columns are read in.
    ! ---------------------------------------------------------------------------

    call apply_prefit_values_and_bounds (prefit, pge_idx, lo_radbnd, up_radbnd, fitvar_rad, errstat)

    ! -----------------------------------------------------------------
    ! Create a condensed array of fitting variables that only contains
    ! the ones that are varied. This considerably reduces the execution
    ! time of the fitting routine.
    ! -----------------------------------------------------------------
    fitvar = 0.0_r8 ; lobnd = 0.0_r8 ; upbnd = 0.0_r8
    DO i = 1, n_fitvar_rad
      idx       = mask_fitvar_rad(i)
      fitvar(i) = fitvar_rad(idx)
      lobnd (i) = lo_radbnd(idx)
      upbnd (i) = up_radbnd(idx)
    END DO

    write (log_msg, *)'fit_radiance:  n_fitvar_rad=',n_fitvar_rad
    call tell_log (2, log_msg)

    ! --------------------------------------------------------------------
    ! Check whether we enough spectral points to carry out the fitting. If
    ! not, call it a bad pixel and return.
    ! --------------------------------------------------------------------
    IF ( (n_fitvar_rad-n_prefit_vars) >= n_rad_wvl_loc ) THEN
      is_bad_pixel = .TRUE.
      RETURN
    END IF

    covar_matrix = r8_missval

    ! ------------------------------------------------------------------------------
    ! For certain PGEs we use a fitting function that corrects the O3 cross sections
    ! with a quadratic polynomial in wavelength, to account for the fact that O3
    ! absorption can vary greatly in magnitude over the fitting window (HH bands).
    ! ------------------------------------------------------------------------------
    if (yn_o3amf_cor) then
      earthshine_spectrum => spectrum_earthshine_o3exp
    else
      earthshine_spectrum => spectrum_earthshine
    endif

    ! To improve computational efficiency, it's useful to avoid adding
    ! spectrum model contributions that are all zeros.  To do that,
    ! it helps to know which params are frozen at zero and which spectrum
    ! model components are all zeros.
    ! Note that defining database_j_is_zero by looking at only the wavelength
    ! sub-range j1:j2 that's actually used seems like a good idea, but
    ! doesn't gain anything in terms of efficiency because the range j1:j2
    ! can change during the fit. The cost of updating database_j_is_zero on
    ! every function evaluation cancels out the gain from avoiding the zeros.
    param_frozen_at_zero = (fitvar_rad == 0.0_r8 .and. lo_radbnd == 0.0_r8 .and. up_radbnd == 0.0_r8)
    database_j_is_zero = .true.
    do j=1, max_rs_idx
      database_j_is_zero(j) = size(database,2).eq.count(database(:,j)==0.0_r8)
    enddo

    radfit_itnum = 0
    j = 0

    call optimizer_open (opt, earthshine_residuals, n_fitvar_rad, errstat, &
                         mode=opt_bounded, tol=tol, epsrel=epsrel, epsabs=epsabs, epsx=epsx, &
                         param_min = lobnd(1:n_fitvar_rad), &
                         param_max = upbnd(1:n_fitvar_rad), &
                         param_mask = mask_fitvar_rad(1:n_fitvar_rad), &
                         max_num_iterations = max_itnum_rad)
    if (errstat < 0) then
      call tell_error (tell_runtime_error, "fit_radiance: optimizer_open failed", errstat)
      return
    endif

    fit_loop: do
      call opt%optimize (opt, fitvar(1:n_fitvar_rad), n_fitvar_rad, &
                         fitres(1:n_rad_wvl_loc), n_rad_wvl_loc, return_status, &
                         cov_matrix=covar_matrix)
      locitnum = opt%num_iterations
      radfit_exval = return_status

      call earthshine_spectrum (n_rad_wvl_loc, Spec%wav_avg, Spec%wavs(1:n_rad_wvl_loc), &
                                fitspec(1:n_rad_wvl_loc), fitvar_rad)

      n_nozero_wgt = MAX ( INT ( ANINT ( SUM(Spec%weights(1:n_rad_wvl_loc)) ) ), 1 )
      mean         = SUM  ( fitres(1:n_rad_wvl_loc) )                 / REAL(n_nozero_wgt, KIND=r8)
      sdev         = SQRT ( SUM ( (fitres(1:n_rad_wvl_loc)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
      loclim       = mean + REAL(fitres_range, KIND=r8)*sdev

      ! ----------------------
      ! Fitting RMS and CHI**2
      ! ----------------------
      IF ( n_nozero_wgt > 0 ) THEN
        rms     = SQRT ( SUM ( fitres(1:n_rad_wvl_loc)**2 ) / REAL(n_nozero_wgt, KIND=r8) )
        ! ---------------------------------------------
        ! This gives the same CHI**2 as the NR routines
        ! ---------------------------------------------
        chisquav = SUM  ( fitres(1:n_rad_wvl_loc)**2 )
      ELSE
        rms      = r8_missval
        chisquav = r8_missval
      END IF

      radfit_itnum = radfit_itnum + locitnum
      j = j + 1

      if (1 < j .and. j <= num_fitres_loops) then
        if (MAXVAL(ABS(fitres(1:n_rad_wvl_loc))) < loclim) exit fit_loop
      else if (j == 1) then
        ! (jch) Original code saved loclim only from the first iteration,
        !       so I kept that behavior, even though it looked odd.
        !       Is this a bug? FIXME?

        ! (jch)  It seems to me that loclim could be negative,
        !        but the original code iterates further only if loclim > 0.
        !        Is this a bug?  FIXME?
        if (.not. ((num_fitres_loops > 0) &
                   .and. (loclim > 0.0_r8) &
                   .and. (n_nozero_wgt >  n_fitvar_rad))) exit fit_loop
      else
        exit fit_loop ! (j > num_fitres_loops)
      endif

      WHERE ( ABS(fitres(1:n_rad_wvl_loc)) >= loclim )
        Spec%weights(1:n_rad_wvl_loc) = downweight
      END WHERE

    enddo fit_loop

    call optimizer_close (opt, errstat)
    if (errstat < 0) then
      call tell_error (tell_runtime_error, "fit_radiance: optimizer_close failed", errstat)
      return
    endif

    ! ---------------------------------------------------------------------
    ! Save correlation information from covariance matrix !gga to real corr
    ! elation. Only correlation of the main variable (retrieved molecule) w
    ! ith the other variables is kept
    ! ---------------------------------------------------------------------
    corrmat = r8_missval
    indx   = fincol_idx(1,1)
    covar_xx = covar_matrix(indx, indx)
    if (covar_xx /= 0.0) then
      DO i = 1, n_fitvar_rad
        IF (covar_matrix (i,i) == 0.0) cycle
        corrmat(i) = (covar_matrix (indx, i ) &
                      / SQRT(covar_xx * covar_matrix (i, i)))
      END DO
    endif

    ! --------------------------------------------------------------------
    ! Save fitting weights for possible use through radiance reference fit
    ! --------------------------------------------------------------------
    adj_wgts(1:n_rad_wvl_loc) = Spec%weights(1:n_rad_wvl_loc)

    ! CCM save fitted spectrum
    fitspc_out(1:n_rad_wvl_loc) = fitspec(1:n_rad_wvl_loc)

    ! --------------------------------------------------------------------
    ! At this point we have to make sure we only report the fitted columns
    ! if EXVAL, the exit variable from the fitting routine, is >=0. All
    ! negative values mean trouble and are likely to have produced strange
    ! uncertainties. In any case, the values for the fitted
    ! columns should not be trusted.
    ! --------------------------------------------------------------------
    SELECT CASE ( radfit_exval )
      ! ======================================================
      ! Anything but EXVAL > 0 means that trouble has occurred,
      ! and the fit most likely is screwed.
      ! ======================================================
    CASE ( :-1 )
      fitcol = r8_missval ; dfitcol = r8_missval
    CASE ( 0: )

      ! ---------------------------
      ! Update common mode spectrum
      ! ---------------------------
      if (.not.yn_reference_fit) then
        CALL compute_common_mode ( &
          yn_reference_fit, ipix, n_rad_wvl_loc, Spec%wavs(1:n_rad_wvl_loc), &
          fitres(1:n_rad_wvl_loc), errstat)
      endif
      if (errstat < 0) return

      ! =====================================================================
      ! Compute the actual number of radiance wavelengths used in the fitting
      ! =====================================================================
      n_fitwav_rad = INT (SUM(1.0_r8 * Spec%weights(1:n_rad_wvl_loc)**2))

      ! --------------------------------------------------------------------------
      ! Assign total column. We have done the preliminary work with FITCOL_IDX.
      !
      ! Reminder about the (admittedly confusing) variable names:
      !
      !   FITCOL:       Fitted column
      !   DFITCOL:      Uncertainty of fitted column
      !   N_FINCOL_IDX: Number of fitting variables that make up the final FITCOL
      !                 (e.g., O3 at more than one temperature).
      !   FINCOL_IDX:   Array of dimension (2,N_FINCOL_IDX*MXS_IDX), where MXS_IDX
      !                 is the maximum fitting sub-index (e.g., AD1_IDX, LBE_IDX).
      !                 FINCOL_IDX(1,*) carries the relative indices of the varied
      !                 final column variables in the array that was passed to
      !                 the fitting routine. FINCOL_IDX(2,*) contains the index
      !                 for the associated reference spectrum; this we need for
      !                 access to the normalization factor. See subroutine
      !                 READ_CONTROL_FILE for assignment of these indices.
      !
      !  TARGET_VAR     Saved fitting parameter values for possible use to remove
      !                 target gas from radiance reference
      ! --------------------------------------------------------------------------

      pm_one = 1.0_r8
      if (yn_doas) pm_one = -1.0_r8

      fitcol = 0.0_r8  ;  dfitcol = 0.0_r8 ; target_var = 1.0_r8
      DO i = 1, n_fincol_idx
        ! --------------------------------------------------
        ! First add the contribution of the diagonal element
        ! --------------------------------------------------
        j1 = fincol_idx(1, i) ; k1 = fincol_idx(2,i)

        fitcol  = fitcol  + pm_one * fitvar(j1) / refspecs_original(k1)%NormFactor
        dfitcol = dfitcol + covar_matrix(j1,j1) / refspecs_original(k1)%NormFactor**2

        IF ( pge_idx == pge_o3_idx ) THEN
          o3fit_cols (k1) = fitvar(j1) / refspecs_original(k1)%NormFactor
          ! -------------------------------------------------------------------------
          ! Just as for the combined column uncertainty, we can't compute O3FIT_DCOLS
          ! if any of the major elements are missing.
          ! -------------------------------------------------------------------------
          IF (((n_fitwav_rad - n_fitvar_rad + n_prefit_vars) <= 0) &
              .OR. (rms == r8_missval) .OR. (n_rad_wvl_loc == 0)) THEN
            o3fit_dcols(k1) = r8_missval
          ELSE
            o3fit_dcols(k1) = &
              rms * SQRT (covar_matrix(j1,j1) &
                          / refspecs_original(k1)%NormFactor**2 &
                          * REAL(n_fitwav_rad, KIND=r8 ) &
                          / REAL(n_fitwav_rad-n_fitvar_rad+n_prefit_vars, &
                                 KIND=r8))
          END IF
        END IF

        ! -------------------------------------------------------------------------
        ! Then add the contributions from off-diagonal elements (correlations)
        ! -------------------------------------------------------------------------
        DO l = i+1, n_fincol_idx
          j2 = fincol_idx(1,l) ; k2 = fincol_idx(2,l)
          dfitcol = dfitcol + 2.0_r8 * covar_matrix(j1,j2) / &
            (refspecs_original(k1)%NormFactor*refspecs_original(k2)%NormFactor)
        END DO

        ! ------------------------------------------
        ! Save fitting parameter value of target gas
        ! ------------------------------------------
        target_var(i) = fitvar(j1)

      END DO

      ! ------------------------------------------------------------
      ! Fitting columns and uncertainties for all fitting parameters
      ! ------------------------------------------------------------
      DO i = 1, n_fitvar_rad
        allfit(i) = fitvar(i)
        allerr(i) = covar_matrix(i,i)
        j = all_radfit_idx(i)
        IF ( j > max_calfit_idx ) THEN
          j = j - max_calfit_idx
          normfac = refspecs_original(j)%NormFactor
          IF ( normfac == 0.0_r8 ) normfac = 1.0_r8
          allfit(i) = allfit(i) / normfac
          allerr(i) = allerr(i) / normfac**2
        END IF
      END DO

      ! ---------------------------------------------------------------------
      ! Computing rationale for DFITCOL:
      !
      ! We are in the CASE branch where iteration has converged to something,
      ! but we don't know whether that something actually makes sense. It is
      ! well possible that some quantities are missing because of numerical
      ! overflows, to name just one. By setting DFITCOL to "missing" we
      ! acknowledge that we have no idea what the fitting uncertainties are
      ! for the current pixel. Setting DFITCOL to a large value (+1E+30 for
      ! example) would be an alternative, but that is just another way of
      ! expressing our ingnorance about what really is going on.
      ! ---------------------------------------------------------------------
      IF (((n_fitwav_rad-n_fitvar_rad+n_prefit_vars) <= 0 ) &
          .OR. (rms == r8_missval) &
          .OR. (n_rad_wvl_loc == 0)) THEN
        dfitcol                = r8_missval
        allerr(1:n_fitvar_rad) = r8_missval
      ELSE
        ! --------------------------------------------------------
        ! Adjust fitting errors according to formulas given in
        ! Numerical Recipes, $15.5: Multiply by SQRT[chi**2/(n-2)]
        ! --------------------------------------------------------
        mfac = SQRT ( chisquav &
                     / REAL(n_fitwav_rad-n_fitvar_rad+n_prefit_vars, KIND=r8))
        dfitcol = SQRT ( dfitcol ) * mfac
        allerr(1:n_fitvar_rad) = SQRT ( allerr(1:n_fitvar_rad) ) * mfac
      END IF

      ! -------------------------------------------------------------------
      ! The following assignment makes sense only because FITVAR_RAD is
      ! updated with FITVAR (using the proper mask) in SPECTRUM_EARTHSHINE.
      !
      ! HOWEVER: The OMI data are rather noisy, and fitting uncertainties
      !          are generally large. Hence this assignment is rather
      !          dangerous. By going back to the initial guess, we may
      !          lose some speed, but we gain predictability.
      ! -------------------------------------------------------------------
      !IF ( (radfit_exval == opt_convergence_good) .AND. &
      !     (fitcol+1.0_r8*dfitcol >= 0.0_r8)                    .AND. &
      !     ( .NOT. yn_reference_fit )                                  )  THEN
      !   fitvar_rad_saved(1:n_max_fitpars) = fitvar_rad(1:n_max_fitpars)
      !ELSE
      fitvar_rad_saved(1:n_max_fitpars) = fitvar_rad_init(1:n_max_fitpars)
      !END IF
    CASE DEFAULT
      ! -------------------------------------------------------------------
      ! We should never reach here, because the above CASE statements cover
      ! all possible values of EXVAL. But better safe than sorry.
      fitcol = r8_missval ; dfitcol = r8_missval
      ! -------------------------------------------------------------------
    END SELECT

    fitcol_saved = fitcol

    RETURN
  END SUBROUTINE fit_radiance

  subroutine earthshine_residuals (this_optimizer, params, num_params, residuals, num_residuals, return_status)
    use OMSAO_variables_module, only: fitvar_rad
    use ctrlvars, only: yn_o3amf_cor
    implicit none
    type(optimizer_type) :: this_optimizer
    real (kind=r8), dimension (:), intent(in) :: params
    real (kind=r8), dimension (:), intent(out) :: residuals
    integer (kind=i4), intent(in) :: num_params, num_residuals
    integer (kind=i4), intent(out) :: return_status
    ! local variables
    procedure(earthshine_spectrum_interface), pointer :: earthshine_spectrum => null()
    integer (kind=i4) :: i, idx

    DO i = 1, num_params
      idx = this_optimizer%param_mask(i)
      fitvar_rad(idx) = params(i)
    END DO

    if (yn_o3amf_cor) then
      earthshine_spectrum => spectrum_earthshine_o3exp
    else
      earthshine_spectrum => spectrum_earthshine
    endif
    call earthshine_spectrum (num_residuals, Spec%wav_avg, &
                              Spec%wavs(1:num_residuals), residuals(1:num_residuals), &
                              fitvar_rad)
    residuals = (Spec%spec(1:num_residuals) - residuals(1:num_residuals)) * Spec%weights(1:num_residuals)

    return_status = 0

  end subroutine earthshine_residuals

  SUBROUTINE spectrum_earthshine (npts, rad_wav_avg, locwvl, fit, rad_fitvar)

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: &
      max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
      lbe_idx, ad2_idx, mxs_idx, wvl_idx, spc_idx,                   &
      bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, sc2_idx, &
      sc3_idx, sin_idx, shi_idx, squ_idx
    USE OMSAO_parameters_module, ONLY: max_spec_pts, downweight
    USE OMSAO_variables_module,  ONLY: &
      n_database_wvl, curr_sol_spec, &
      database, curr_xtrack_pixnum
    use ctrlvars, only: yn_radiance_reference, yn_spectrum_norm, &
      yn_doas, yn_newshift, yn_solar_comp
    USE OMSAO_prefitcol_module,  ONLY:  apply_prefit_values
    USE OMSAO_omidata_module,      ONLY: omi_solcal_pars
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_errstat_module
    USE OMSAO_solcomp_module, ONLY: soco_compute
    USE sao_pge_utils, ONLY: interpolation
    USE arrayutils, only: array_locate_r8, array_sort_r8

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER (KIND=i4),                    INTENT (IN)    :: npts
    REAL    (KIND=r8),                    INTENT (IN)    :: rad_wav_avg
    REAL    (KIND=r8), DIMENSION (:),     INTENT (INOUT) :: rad_fitvar
    REAL    (KIND=r8), DIMENSION (npts),  INTENT (IN)    :: locwvl
    !REAL    (KIND=r8), DIMENSION (max_rs_idx,max_spec_pts), INTENT (IN) :: database

    ! ================
    ! Output variables
    ! ================
    REAL (KIND=r8), DIMENSION (npts), INTENT (OUT) :: fit

    ! ===============
    ! Local variables
    ! ===============
    REAL    (KIND=r8), PARAMETER                  :: expmax = REAL(MAXEXPONENT(1.0_r4), KIND=r8)
    REAL    (KIND=r8), PARAMETER                  :: expmin = REAL(MINEXPONENT(1.0_r4), KIND=r8)
    LOGICAL                                       :: did_full_range, is_solsynth
    INTEGER (KIND=i4)                             :: i, j, errstat, j1, j2, n_sunpos
    REAL    (KIND=r8)                             :: shift, squeeze, soco_shi

    ! Try to save some stack space by reusing some arrays via pointers.  --JED
    REAL (KIND=r8), DIMENSION(npts), TARGET       :: tmpspace
    REAL (KIND=r8), POINTER                       :: del(:), sunspec_ss(:), sumexp(:)
    ! REAL    (KIND=r8), DIMENSION (npts)           :: del, sunspec_ss, sumexp
    REAL    (KIND=r8), DIMENSION (npts)           :: database_j, fit_final_add_on
    REAL    (KIND=r8), DIMENSION (npts)           :: locwvl_shift
    REAL    (KIND=r8), DIMENSION (max_spec_pts)   :: sunpos_ss, sunspec_loc, sunspec_save
    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=19), PARAMETER :: modulename = 'spectrum_earthshine'

    SAVE sunspec_save

    !     Calculate the spectrum:
    !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
    !     1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
    !     OMI data.

    errstat = pge_errstat_ok

    ! ----------------------------------------------------------------------------
    ! Here is a logical to determine whether we need to compute a "sythetic"
    ! solar spectrum from the solar composite. The cases for YES are
    !
    ! (1) We are using the solar composite and are NOT doing a radiance reference
    ! (2) We are using the solar composite and ARE doing a radiance reference, and
    !     this happens to be the radiance reference fit.
    ! ----------------------------------------------------------------------------
    ! ---------------------------------------------------------------------
    ! The solar composite spectrum may have an additional shift, which was
    ! determined during the solar wavelength calibration. This needs to be
    ! taken into account when computing the spectra. But careful: It should
    ! NOT be added to the local wavelength array, since that is related to
    ! the radiance only. The Solar Composite shift must be subtracted from
    ! the wavelength array, hence the negative sign.
    ! ---------------------------------------------------------------------
    !IF (( yn_solar_comp .AND. (.NOT. yn_radiance_reference) ) &
    !    .OR. (yn_solar_comp .AND. yn_radiance_reference &
    !          .AND. yn_reference_fit)) ) THEN
    ! The above test can be simplified to the following: --JED
    IF (yn_solar_comp .and. (.not.yn_radiance_reference)) then
      is_solsynth = .TRUE.
      soco_shi = -omi_solcal_pars(shi_idx,curr_xtrack_pixnum)
    ELSE
      is_solsynth = .FALSE.
      soco_shi = 0.0_r8
    END IF

    shift   = rad_fitvar(shi_idx)
    squeeze = rad_fitvar(squ_idx)

    ! -------------------------------------
    ! Dealing with any pre-fitted variables
    ! -------------------------------------
    call apply_prefit_values (rad_fitvar)

    ! -----------------------------------------------------------------------------------------
    ! Assign current solar spectrum to local arrays. This depends on whether we are using
    ! actual measured solar spectra or solar composites. Since there is no point to interpolate
    ! already interpolated spectra, we use the original solar composites here as base for the
    ! interpolation to the final radiance wavelengths.
    ! -----------------------------------------------------------------------------------------
    n_sunpos                = n_database_wvl
    sunpos_ss  (1:n_sunpos) = curr_sol_spec(1:n_sunpos, wvl_idx)
    sunspec_loc(1:n_sunpos) = curr_sol_spec(1:n_sunpos, spc_idx)

    ! ----------------------------------------------
    ! Sort local arrays - important to pass EZspline
    ! ----------------------------------------------
    ! Most of the time, the array is already sorted.  Try to avoid the function
    ! call overhead.  --JED
    DO i=2, n_sunpos
      IF (sunpos_ss(i-1) < sunpos_ss(i)) CYCLE
      CALL array_sort_r8 ( n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos) )
      EXIT
    ENDDO

    ! ---------------------------------------------------------------------
    ! Apply Shift&Squeeze
    ! Changed to include Xiong comments (gga) if yn_newshift equal .true. :
    ! Lambda = Lambda * (1 + squeeze) + shift - solar_wavel_avg * squeeze
    ! ---------------------------------------------------------------------
    j1 = -1; j2 = -1
    IF ( squeeze == 0.0_r8 .AND. is_solsynth ) THEN
      locwvl_shift(1:npts) = locwvl(1:npts) - shift
      CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(   1), 'GE', j1 )
      CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(npts), 'LE', j2 )
    ELSE IF (yn_newshift) THEN !gga
      sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) +       &
        shift - rad_wav_avg * squeeze
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 ) !gga
    ELSE
      sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) + shift
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 )
    END IF

    sunspec_ss => tmpspace
    ! ---------------------------------------------------------------------
    ! Re-sample the solar reference spectrum to the current radiance grid
    ! ---------------------------------------------------------------------
    !
    ! The endpoints may be problematic due to no-strict ascendence. If that
    ! happens, exclude end-points.
    ! ---------------------------------------------------------------------

    IF ( j1 <= 0 .OR. j2 <= 0 ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
        modulename//f_sep//'Resampling to Radiance Grid -- no solar spectrum!!!', &
        vb_lev_default, errstat )
    ELSE

      IF ( squeeze /= saved_squeeze .OR. shift /= saved_shift ) THEN

        IF ( squeeze == 0.0_r8 .AND. is_solsynth ) THEN
          CALL soco_compute ( &
            yn_spectrum_norm, curr_xtrack_pixnum, npts, &
            locwvl_shift(1:npts)+soco_shi, sunspec_ss(1:npts) )
        ELSE
          CALL interpolation (                                                 &
            n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos),       &
            npts, locwvl(1:npts), sunspec_ss(1:npts), 'endpoints', 0.0_r8,  &
            did_full_range, errstat                                            )
          CALL error_check ( &
            errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL,      &
            modulename//f_sep//'Resampling to Radiance Grid -- interpolation', &
            vb_lev_default, errstat )
        END IF
        sunspec_save(1:npts) = sunspec_ss(1:npts)
        saved_shift          = shift
        saved_squeeze        = squeeze
      ELSE
        sunspec_ss(1:npts) = sunspec_save(1:npts)
      END IF

    END IF

    ! Add up the contributions, with solar intensity as rad_fitvar(sin_idx), trace
    ! species beginning at rad_fitvar(SQU_IDX+1), to include possible linear and
    ! Beer's law forms.  Do these as linear-Beer's-linear. In order to
    ! do DOAS I only need to be careful to include just linear
    ! contributions, since I already high-pass filtered them.

    IF ( j1 > 1    )  Spec%weights(1:j1-1)    = downweight
    IF ( j2 < npts )  Spec%weights(j2+1:npts) = downweight

    fit = 0.0_r8

    ! ==================================================================
    ! For BOAS or any wavelength calibration, we have the following line
    ! ==================================================================

    fit(j1:j2) = rad_fitvar(sin_idx) * sunspec_ss(j1:j2)

    !     DOAS here - the spectrum to be fitted needs to be re-defined:
    IF ( yn_doas ) THEN

      i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

      ! For DOAS, rad_fitvar(SIN_IDX) should == 1., and not be varied
      fit(j1:j2) = &
        rad_fitvar(sin_idx) * LOG ( sunspec_ss(j1:j2) ) &
        + rad_fitvar(i) * (database(j1:j2, ring_idx) / sunspec_ss (j1:j2))

      DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. j /= ring_idx ) THEN
          if (database_j_is_zero(j)) cycle
          i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
          fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j1:j2,j)
        END IF
      END DO

    ELSE
      sumexp => tmpspace
      sumexp(j1:j2) = 0.0_r8
      fit_final_add_on(j1:j2) = 0.0_r8
      DO j = 1, max_rs_idx
        IF ( j.eq.solar_idx ) CYCLE
        if (database_j_is_zero(j)) cycle
        database_j(j1:j2) = database(j1:j2,j)
        ! -----------------------------
        ! Initial add-on contributions.
        ! -----------------------------
        i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
        if (.not. param_frozen_at_zero(i)) then
          fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database_j(j1:j2)
        endif

        ! -----------------------------
        ! Beer's law contributions.
        ! -----------------------------
        ! ---------------------------------------------------------------
        ! We sum over all contributions and take the EXP only at the end.
        ! This should shave a few seconds off the execution time.
        ! ---------------------------------------------------------------
        i = max_calfit_idx + (j-1)*mxs_idx + lbe_idx
        if (.not. param_frozen_at_zero(i)) then
          sumexp(j1:j2) = sumexp(j1:j2) - rad_fitvar(i)*database_j(j1:j2)
        endif

        ! Final add-on contributions.
        i = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
        if (.not. param_frozen_at_zero(i)) then
          fit_final_add_on(j1:j2) = fit_final_add_on(j1:j2) + &
            rad_fitvar(i) * database_j(j1:j2)
        endif
      END DO

      WHERE ( sumexp(j1:j2) >= expmax )
        sumexp(j1:j2) = expmax
      ENDWHERE
      WHERE ( sumexp(j1:j2) <= expmin )
        sumexp(j1:j2) = expmin
      ENDWHERE

      ! FIXME?? (jch):  Note that the common mode contribution is included in
      ! fit_final_add_on(j1:j2), but after that, the model is then multiplied
      ! by a scaling polynomial. As I understand it, the purpose of the common
      ! mode is to take out features that consistently appear in the fit residuals.
      ! For that reason, the common mode presumably represents an additive term
      ! in the final model. Why then is it being multiplied by a scaling polynomial
      ! that's a function of wavelength?
      ! Shouldn't the common mode be added _after_ the scaling??

      fit(j1:j2) = fit(j1:j2) * EXP(sumexp(j1:j2)) + fit_final_add_on(j1:j2)

    ENDIF

    ! Add the scaling.
    del => tmpspace
    ! Use the form: A+BX+CX^2+DX^3 = A + X*(B + X*(C + X*D))
    del(j1:j2) = locwvl(j1:j2) - rad_wav_avg
    fit(j1:j2) = fit(j1:j2) &
      * (rad_fitvar(sc0_idx) + &
         del(j1:j2) * (rad_fitvar(sc1_idx) + &
                       del(j1:j2) * (rad_fitvar(sc2_idx) + &
                                     del(j1:j2) * rad_fitvar(sc3_idx))))

    ! Add baseline parameters.
    fit(j1:j2) = fit(j1:j2)                                        + &
      rad_fitvar(bl0_idx)                                       + &
      rad_fitvar(bl1_idx) * del(j1:j2)                          + &
      rad_fitvar(bl2_idx) * del(j1:j2)*del(j1:j2)               + &
      rad_fitvar(bl3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2)

    ! This form is better than the above, but introduces differences
    ! in the last digits of the output, breaking simple-minded diff-based
    ! regression tests.
    !fit(j1:j2) = fit(j1:j2) + rad_fitvar(bl0_idx) &
    !  + del(j1:j2) * (rad_fitvar(bl1_idx) + &
    !                  del(j1:j2) * (rad_fitvar(bl2_idx) + &
    !                                del(j1:j2) * rad_fitvar(bl3_idx)))

    ! ----------------------------------------------------------------
    ! Final sanity check: If the various multiplications and additions
    ! have lead to NaN values, we set those to ZERO. This is somewhat
    ! experimental, and if we come up with a better way of doing this,
    ! then the logic below should be changed accordingly.
    ! ----------------------------------------------------------------
    !WHERE ( .NOT. ( fit(j1:j2) > -HUGE(1.0_r8) .AND. fit(j1:j2) < HUGE(1.0_r8) ) )!
    !  fit(j1:j2) = 0.0_r8!
    WHERE (.NOT. (abs(fit(j1:j2)) < HUGE(1.0_r8)))
      fit(j1:j2) = 0.0_r8
    END WHERE

    RETURN
  END SUBROUTINE spectrum_earthshine

  SUBROUTINE spectrum_earthshine_o3exp (npts, rad_wav_avg, locwvl, fit, rad_fitvar)

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: &
      max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
      lbe_idx, ad2_idx, mxs_idx, wvl_idx, spc_idx,                   &
      bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, sc2_idx, &
      sc3_idx, sin_idx, shi_idx, squ_idx, &
      o3_t1_idx, o3_t2_idx, o3_t3_idx
    USE OMSAO_parameters_module, ONLY: max_spec_pts, downweight
    USE OMSAO_variables_module,  ONLY: &
      n_database_wvl, curr_sol_spec, &
      database, curr_xtrack_pixnum
    use ctrlvars, only: yn_radiance_reference, yn_spectrum_norm, yn_doas, &
      yn_solar_comp, yn_newshift
    USE OMSAO_prefitcol_module,  ONLY:  apply_prefit_values
    USE OMSAO_omidata_module,      ONLY: omi_solcal_pars
    USE cache_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_errstat_module
    USE OMSAO_solcomp_module, ONLY: soco_compute
    USE sao_pge_utils, ONLY: interpolation
    USE arrayutils, only: array_locate_r8, array_sort_r8

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER (KIND=i4),                    INTENT (IN)    :: npts
    REAL    (KIND=r8),                    INTENT (IN)    :: rad_wav_avg
    REAL    (KIND=r8), DIMENSION (:),     INTENT (INOUT) :: rad_fitvar
    REAL    (KIND=r8), DIMENSION (npts),  INTENT (IN)    :: locwvl

    !REAL    (KIND=r8), DIMENSION (max_rs_idx,max_spec_pts), INTENT (IN) :: database
    ! database was passed as an argument.  However, it came from the OMSAO_variables_module, and
    ! had a very different size.
    ! ================
    ! Output variables
    ! ================
    REAL (KIND=r8), DIMENSION (npts), INTENT (OUT) :: fit

    ! ===============
    ! Local variables
    ! ===============
    REAL    (KIND=r8), PARAMETER                  :: expmax = REAL(MAXEXPONENT(1.0_r4), KIND=r8)
    REAL    (KIND=r8), PARAMETER                  :: expmin = REAL(MINEXPONENT(1.0_r4), KIND=r8)
    LOGICAL                                       :: did_full_range, is_solsynth
    INTEGER (KIND=i4)                             :: i, j, errstat, j1, j2, n_sunpos, k1, k2
    REAL    (KIND=r8)                             :: shift, squeeze, soco_shi
    REAL    (KIND=r8), DIMENSION (npts)           :: del, sunspec_ss, tmpexp, sumexp
    REAL    (KIND=r8), DIMENSION (npts)           :: locwvl_shift
    REAL    (KIND=r8), DIMENSION (max_spec_pts)   :: sunpos_ss, sunspec_loc, sunspec_save

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=25), PARAMETER :: modulename = 'spectrum_earthshine_o3exp'

    SAVE sunspec_save

    !  Calculate the spectrum:
    !  First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
    !  1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
    !  OMI data.

    errstat = pge_errstat_ok

    ! ----------------------------------------------------------------------------
    ! Here is a logical to determine whether we need to compute a "sythetic"
    ! solar spectrum from the solar composite. The cases for YES are
    !
    ! (1) We are using the solar composite and are NOT doing a radiance reference
    ! (2) We are using the solar composite and ARE doing a radiance reference, and
    !     this happens to be the radiance reference fit.
    ! ----------------------------------------------------------------------------
    ! ---------------------------------------------------------------------
    ! The solar composite spectrum may have an additional shift, which was
    ! determined during the solar wavelength calibration. This needs to be
    ! taken into account when computing the spectra. But careful: It should
    ! NOT be added to the local wavelength array, since that is related to
    ! the radiance only. The Solar Composite shift must be subtracted from
    ! the wavelength array, hence the negative sign.
    ! ---------------------------------------------------------------------

    !IF (( yn_solar_comp .AND. (.NOT. yn_radiance_reference) ) &
    !    .OR. (yn_solar_comp .AND. yn_radiance_reference &
    !          .AND. yn_reference_fit)) ) THEN
    ! The above test can be simplified to the following: --JED
    IF (yn_solar_comp .and. (.not.yn_radiance_reference)) then
      is_solsynth = .TRUE.
      soco_shi = -omi_solcal_pars(shi_idx,curr_xtrack_pixnum)
    ELSE
      is_solsynth = .FALSE.
      soco_shi = 0.0_r8
    END IF

    shift   = rad_fitvar(shi_idx)
    squeeze = rad_fitvar(squ_idx)

    call apply_prefit_values (rad_fitvar)

    ! ---------------------------------------------------------------------------------
    ! Assign current solar spectrum to local arrays. This depends on whether we are
    ! using actual measured solar spectra or solar composites. Since there is no point
    ! to interpolate already interpolated spectra, we use the original solar composites
    ! here as base for the interpolation to the final radiance wavelengths.
    ! ---------------------------------------------------------------------------------
    n_sunpos                = n_database_wvl
    sunpos_ss  (1:n_sunpos) = curr_sol_spec(1:n_sunpos,wvl_idx)
    sunspec_loc(1:n_sunpos) = curr_sol_spec(1:n_sunpos,spc_idx)

    ! ----------------------------------------------
    ! Sort local arrays - important to pass EZspline
    ! ----------------------------------------------
    CALL array_sort_r8 ( n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos) )

    ! ---------------------------------------------------------------------
    ! Apply Shift&Squeeze
    ! Changed to include Xiong comments (gga) if yn_newshift equal .true. :
    ! Lambda = Lambda * (1 + squeeze) + shift - solar_wavel_avg * squeeze
    ! ---------------------------------------------------------------------
    j1 = -1; j2 = -1
    IF ( squeeze == 0.0_r8 .AND. is_solsynth ) THEN
      locwvl_shift(1:npts) = locwvl(1:npts) - shift
      CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(   1), 'GE', j1 )
      CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(npts), 'LE', j2 )
    ELSE IF (yn_newshift .EQV. .true.) THEN !gga
      sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) +       &
        shift - rad_wav_avg * squeeze
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 ) !gga
    ELSE
      sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) + shift
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
      CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 )
    END IF

    ! ---------------------------------------------------------------------
    ! Re-sample the solar reference spectrum to the current radiance grid
    ! ---------------------------------------------------------------------
    !
    ! The endpoints may be problematic due to no-strict ascendence. If that
    ! happens, exclude end-points.
    ! ---------------------------------------------------------------------

    IF ( j1 <= 0 .OR. j2 <= 0 ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
        modulename//f_sep//'Resampling to Radiance Grid -- no solar spectrum!!!', &
        vb_lev_default, errstat )
    ELSE

      IF ( squeeze /= saved_squeeze .OR. shift /= saved_shift ) THEN

        IF ( squeeze == 0.0_r8 .AND. is_solsynth ) THEN
          CALL soco_compute ( &
            yn_spectrum_norm, curr_xtrack_pixnum, npts, &
            locwvl_shift(1:npts)+soco_shi, sunspec_ss(1:npts))
        ELSE
          CALL interpolation (                                                         &
            n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos),               &
            npts, locwvl(1:npts), sunspec_ss(1:npts), 'endpoints', 0.0_r8, &
            did_full_range, errstat                                                   )
          CALL error_check ( &
            errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
            modulename//f_sep//'Resampling to Radiance Grid -- interpolation', &
            vb_lev_default, errstat )
        END IF

        sunspec_save(1:npts) = sunspec_ss(1:npts)
        saved_shift   = shift
        saved_squeeze = squeeze
      ELSE
        sunspec_ss(1:npts) = sunspec_save(1:npts)
      END IF

    END IF

    !     Add up the contributions, with solar intensity as rad_fitvar(sin_idx), trace
    !     species beginning at rad_fitvar(SQU_IDX+1), to include possible linear and
    !     Beer's law forms.  Do these as linear-Beer's-linear. In order to
    !     do DOAS I only need to be careful to include just linear
    !     contributions, since I already high-pass filtered them.

    IF ( j1 > 1    )  Spec%weights(1:j1-1)    = downweight
    IF ( j2 < npts )  Spec%weights(j2+1:npts) = downweight

    fit = 0.0_r8

    ! --------------------------------------------------------
    ! Compute abcissae for exponential x-section modification:
    ! Values between -1 and +1 on the fitting wavelength grid.
    ! --------------------------------------------------------
    del(j1:j2) = (locwvl(j1:j2) - locwvl(j1))/(locwvl(j2)-locwvl(j1)) - 0.5_r8

    ! ==================================================================
    ! For BOAS or any wavelength calibration, we have the following line
    ! ==================================================================

    fit(j1:j2) = rad_fitvar(sin_idx) * sunspec_ss(j1:j2)

    !     DOAS here - the spectrum to be fitted needs to be re-defined:
    IF ( yn_doas ) THEN

      i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

      ! For DOAS, rad_fitvar(SIN_IDX) should == 1., and not be varied
      fit(j1:j2) = &
        rad_fitvar(sin_idx) * LOG ( sunspec_ss(j1:j2) ) &
        + rad_fitvar(i) * (database(j1:j2, ring_idx) / sunspec_ss (j1:j2))

      DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. j /= ring_idx  .AND. &
            j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
          if (database_j_is_zero(j)) cycle
          i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
          if (.not. param_frozen_at_zero(i)) then
            fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j1:j2,j)
          endif
        END IF
      END DO

    ELSE
      ! -----------------------------
      ! Initial add-on contributions.
      ! -----------------------------
      DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. &
            j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
          if (database_j_is_zero(j)) cycle
          i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
          if (.not. param_frozen_at_zero(i)) then
            fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j1:j2,j)
          endif
        END IF
      END DO
      ! -----------------------------
      ! Beer's law contributions.
      ! -----------------------------
      ! ---------------------------------------------------------------
      ! We sum over all contributions and take the EXP only at the end.
      ! This should shave a few seconds off the execution time.
      ! ---------------------------------------------------------------
      sumexp(j1:j2) = 0.0_r8
      DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
          if (database_j_is_zero(j)) cycle
          tmpexp = 0.0_r8
          i = max_calfit_idx + (j-1)*mxs_idx + lbe_idx
          if (.not.param_frozen_at_zero(i)) then
            IF ( j == o3_t1_idx .OR. j == o3_t2_idx .OR. j == o3_t3_idx ) THEN
              k1 = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
              k2 = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
              tmpexp(j1:j2) = rad_fitvar(i)*database(j1:j2,j) *  &
                (1.0_r8 + rad_fitvar(k1)*del(j1:j2) + &
                 rad_fitvar(k2)*del(j1:j2)*del(j1:j2))
            ELSE
              tmpexp(j1:j2) = rad_fitvar(i)*database(j1:j2,j)
            END IF
          endif

          WHERE ( tmpexp(j1:j2) >= expmax )
            tmpexp(j1:j2) = expmax
          ENDWHERE
          WHERE ( tmpexp(j1:j2) <= expmin )
            tmpexp(j1:j2) = expmin
          ENDWHERE
          sumexp(j1:j2) = sumexp(j1:j2) - tmpexp(j1:j2)
        END IF
      END DO
      WHERE ( sumexp(j1:j2) >= expmax )
        sumexp(j1:j2) = expmax
      ENDWHERE
      WHERE ( sumexp(j1:j2) <= expmin )
        sumexp(j1:j2) = expmin
      ENDWHERE
      fit(j1:j2) = fit(j1:j2) * EXP(sumexp(j1:j2))

      ! Final add-on contributions.
      DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. &
            j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
          i = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
          if (database_j_is_zero(j)) cycle
          if (.not.param_frozen_at_zero(i)) then
            fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j1:j2,j)
          endif
        END IF
      END DO

    END IF

    ! ----------------------------------------
    ! Compute abcissae for closure polynomials
    ! ----------------------------------------
    del(j1:j2) = locwvl(j1:j2) - rad_wav_avg

    ! Add the scaling.
    fit(j1:j2) = fit(j1:j2) * ( &
      rad_fitvar(sc0_idx)                                       + &
      rad_fitvar(sc1_idx) * del(j1:j2)                          + &
      rad_fitvar(sc2_idx) * del(j1:j2)*del(j1:j2)               + &
      rad_fitvar(sc3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2) )

    ! Add baseline parameters.
    fit(j1:j2) = fit(j1:j2)                                        + &
      rad_fitvar(bl0_idx)                                       + &
      rad_fitvar(bl1_idx) * del(j1:j2)                          + &
      rad_fitvar(bl2_idx) * del(j1:j2)*del(j1:j2)               + &
      rad_fitvar(bl3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2)

    ! ----------------------------------------------------------------
    ! Final sanity check: If the various multiplications and additions
    ! have lead to NaN values, we set those to ZERO. This is somewhat
    ! experimental, and if we come up with a better way of doing this,
    ! then the logic below should be changed accordingly.
    ! ----------------------------------------------------------------
    WHERE ( .NOT. ( fit(j1:j2) > -HUGE(1.0_r8) .AND. fit(j1:j2) < HUGE(1.0_r8) ) )
      fit(j1:j2) = 0.0_r8
    END WHERE

    RETURN
  END SUBROUTINE spectrum_earthshine_o3exp

END MODULE
