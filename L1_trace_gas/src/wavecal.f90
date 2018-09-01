module wavecal

  use OMSAO_precision_module, only: i2, i4, r8
  use optimizer_interface_module
  use tell_module
  !use OMSAO_variables_module, only: sol_wav_avg
  use OMSAO_indices_module, only: MAX_CAL_PARMS

  implicit none
  real (kind=r8), dimension(MAX_CAL_PARMS) :: cal_parms
  real (kind=r8), dimension (:), allocatable :: &
    cal_wavelengths, cal_weights, cal_spectrum

  private
  public wavecal_fit

  real (kind=r8), private :: private_avg_wavelength

contains

  subroutine deallocate_module_variables (errstat)
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return
    if (allocated (cal_wavelengths)) deallocate(cal_wavelengths, stat=errstat)
    if (allocated (cal_spectrum) .and. errstat == 0) &
         deallocate(cal_spectrum, stat=errstat)
    if (allocated (cal_weights) .and. errstat == 0) &
         deallocate(cal_weights, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "deallocate_module_variables failed", errstat)
      return
    endif
  end subroutine

  subroutine allocate_module_variables (wvls, spec, wgts, n, errstat)

    implicit none
    real (kind=r8), dimension(:), intent(in) :: wvls, spec, wgts
    integer (kind=i4), intent(in) :: n
    integer, intent(inout) :: errstat
    !
    integer :: locerr

    if (errstat /= 0) return

    call deallocate_module_variables (errstat)
    if (errstat /= 0) return

    allocate(cal_wavelengths(n), cal_spectrum(n), cal_weights(n), stat=locerr)
    if (locerr /= 0) then
      call tell_error (tell_malloc_error, "allocate_module_variables: allocate failed", errstat)
      return
    endif

    cal_wavelengths(1:n) = wvls(1:n)
    cal_spectrum(1:n) = spec(1:n)
    cal_weights(1:n) = wgts(1:n)

  end subroutine

  subroutine wavecal_residuals (this_optimizer, params, num_params, &
                                residuals, num_wavelengths, return_status)
    implicit none
    type(optimizer_type) :: this_optimizer
    real (kind=r8), dimension (:), intent(in) :: params
    real (kind=r8), dimension (:), intent(out) :: residuals
    integer (kind=i4), intent(in) :: num_params, num_wavelengths
    integer (kind=i4), intent(out) :: return_status
    ! local variables
    integer (kind=i4) :: i, idx, err

    ! unpack fit parameters
    DO i = 1, num_params
      idx = this_optimizer%param_mask(i)
      cal_parms(idx) = params(i)
    END DO

    err = 0

    CALL spectrum_solar (num_wavelengths, private_avg_wavelength, & !sol_wav_avg, &
                         cal_wavelengths(1:num_wavelengths), &
                         residuals(1:num_wavelengths), cal_parms, err)

    if (err /= 0) then
      return_status = err
      return
    endif

    residuals(1:num_wavelengths) = cal_weights(1:num_wavelengths) &
      * (cal_spectrum(1:num_wavelengths) - residuals(1:num_wavelengths))

    return_status = 0

  end subroutine wavecal_residuals

  ! ---------------------------------------------------------------------------

  subroutine spectrum_solar (npoints, solar_wavel_avg, locwvl, fit, loc_cal_parms, err)

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: &
      solar_idx, &
      bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, &
      sc2_idx, sc3_idx, sin_idx, hwe_idx, asy_idx, shi_idx, squ_idx
    USE OMSAO_variables_module,  ONLY: &
      refspecs_original, curr_xtrack_pixnum
    use ctrlvars, only: yn_spectrum_norm, yn_newshift
    use slitfunction, only : slitfunction_convolve
    USE cache_module, ONLY: saved_shift, saved_squeeze, saved_hwe, saved_asy
    !USE OMSAO_errstat_module
    USE sao_pge_utils, ONLY: interpolation
    use OMSAO_parameters_module, only: MAX_SPEC_PTS

    IMPLICIT NONE

    INTEGER (KIND=i4),                      INTENT (IN)    :: npoints
    REAL    (KIND=r8),                      INTENT (IN)    :: solar_wavel_avg
    REAL    (KIND=r8), DIMENSION (:), INTENT (IN) :: loc_cal_parms
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: locwvl
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: fit
    integer (kind=i4), intent(inout) :: err

    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL                                                 :: did_full_range
    REAL    (KIND=r8), DIMENSION (npoints)                  :: del, sunspec_ss
    INTEGER (KIND=i4)                                       :: npts
    ! Shorthands for solar reference spectrum
    REAL    (KIND=r8), DIMENSION (refspecs_original(solar_idx)%nPoints) :: &
      solar_wvls, solar_spec
    real (kind=r8), dimension (MAX_SPEC_PTS), save :: &
      saved_solar_spec_convolved = 0.0_r8

    if (err /= 0) return

    npts               = refspecs_original(solar_idx)%nPoints
    solar_wvls(1:npts) = refspecs_original(solar_idx)%RefSpecWavs(1:npts)
    solar_spec(1:npts) = refspecs_original(solar_idx)%RefSpecData(1:npts)
    IF ( .NOT. yn_spectrum_norm ) &
      solar_spec(1:npts) = solar_spec(1:npts) * refspecs_original(solar_idx)%NormFactor

    ! =========================================================================
    !     Spectrum Calculation for Solar and Radiance Wavelength Calibration
    ! =========================================================================

    !     Calculate the spectrum:
    !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
    !     1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
    !     OMI data.
    !     Now, after Xiong recommendation if yn_newfit equal true then (gga):
    !     Lambda = Lambda * (1 + squeeze) + shift - solar_wavel_avg * squeeze

    IF (yn_newshift) THEN ! gga
      solar_wvls(1:npts) = solar_wvls(1:npts) * (1.0_r8 + loc_cal_parms(squ_idx)) &
        +  loc_cal_parms(shi_idx) - solar_wavel_avg * loc_cal_parms(squ_idx)
    ELSE ! gga
      solar_wvls(1:npts) = solar_wvls(1:npts) * (1.0_r8 + loc_cal_parms(squ_idx)) &
        + loc_cal_parms(shi_idx)
    END IF

    ! ----------------------------------------------
    ! Convolve only if we don't do a solar composite
    ! ----------------------------------------------
    !IF ( yn_use_labslitfunc ) THEN
    !  ! ------------------------------------------------------------------------
    !  ! Only if either SHIFT or SQUEEZE have changed from the last iteration do
    !  ! we need to reconvolve the solar spectrum.
    !  !
    !  ! The choice of OMI lab slit function vs. Gaussian is made in the fitting
    !  ! control file: If the initial value of FITVAR(hwe_idx) is 0.0 then we are
    !  ! using the lab measurements, otherwise the Gaussian.
    !  ! ------------------------------------------------------------------------
    !  IF ( loc_cal_parms(squ_idx) /= saved_squeeze .OR. &
    !    loc_cal_parms(shi_idx) /= saved_shift ) THEN
    !    saved_squeeze = loc_cal_parms(squ_idx)
    !    saved_shift   = loc_cal_parms(shi_idx)
    !    saved_solar_spec_convolved = 0.0_r8
    !    CALL omi_slitfunc_convolve (                                  &
    !      curr_xtrack_pixnum, npts, solar_wvls(1:npts),             &
    !      solar_spec(1:npts), saved_solar_spec_convolved(1:npts), errstat )
    !    CALL error_check ( &
    !      errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
    !      modulename//f_sep//'Convolution', vb_lev_default, errstat )
    !    IF ( errstat >= pge_errstat_error ) RETURN
    !  END IF
    !ELSE
    !  CALL asymmetric_gaussian_sf (                                           &
    !    npts, loc_cal_parms(hwe_idx), loc_cal_parms(asy_idx),                    &
    !    solar_wvls(1:npts), solar_spec(1:npts), saved_solar_spec_convolved(1:npts) )
    !END IF

!    if (loc_cal_parms(squ_idx) /= saved_squeeze &
!        .OR. loc_cal_parms(shi_idx) /= saved_shift) then
    if (loc_cal_parms(hwe_idx) /= saved_hwe &
        .OR. loc_cal_parms(asy_idx) /= saved_asy) then
      ! The slit-function convolved solar spectrum is cached in
      ! saved_solar_spec_convolved and need not be updated unless the
      ! shift/squeeze parameters have changed, modifying the wavelength grid.
      saved_squeeze = loc_cal_parms(squ_idx)
      saved_shift   = loc_cal_parms(shi_idx)
      saved_hwe   = loc_cal_parms(hwe_idx)
      saved_asy   = loc_cal_parms(asy_idx)
      saved_solar_spec_convolved = 0.0_r8
      CALL slitfunction_convolve ( &
        npts, solar_wvls(1:npts), solar_spec(1:npts), &
        saved_solar_spec_convolved(1:npts), &
        curr_xtrack_pixnum, loc_cal_parms ([hwe_idx, asy_idx]), 2, &
        err)
      if (err /= 0) return
    endif

    ! =============================================
    ! Broadening and re-sampling of solar spectrum:
    ! =============================================
    ! Case for wavelength fitting of irradiance and radiance
    ! Broaden the solar reference by the hw1e value
    ! ------------------------------------------------------

    ! ------------------------------------------------------
    ! Re-sample the solar reference spectrum to the OMI grid
    ! ------------------------------------------------------
    CALL interpolation ( &
      npts, solar_wvls(1:npts), saved_solar_spec_convolved(1:npts), &
      npoints, locwvl(1:npoints), sunspec_ss(1:npoints), 'endpoints', 0.0_r8, &
      did_full_range, err )
    if (err /= 0) then
      call tell_log (0, "spectrum_solar: interpolation failed while resampling to solar grid")
      return
    endif

    ! --------------------------------------------------------------------
    ! Add up the contributions, with solar intensity as FITVAR (SIN_IDX),
    ! to include possible linear and Beer's law forms.  Do these as
    ! linear-Beer's-linear. In order to do DOAS we only need to be careful
    ! to include just linear contributions, since I already high-pass
    ! filtered them.
    ! --------------------------------------------------------------------

    ! -----------
    !  Doing BOAS
    ! -----------
    fit(1:npoints) = loc_cal_parms(sin_idx) * sunspec_ss(1:npoints)

    ! ----------------
    ! Add the scaling.
    ! ----------------
    del(1:npoints) = locwvl(1:npoints) - solar_wavel_avg
    fit(1:npoints) = fit(1:npoints) * ( &
      loc_cal_parms(sc0_idx)                                               + &
      loc_cal_parms(sc1_idx) * del(1:npoints)                              + &
      loc_cal_parms(sc2_idx) * del(1:npoints)*del(1:npoints)               + &
      loc_cal_parms(sc3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints) )

    ! ------------------------
    ! Add baseline parameters.
    ! ------------------------
    fit(1:npoints) = fit(1:npoints) + &
      loc_cal_parms(bl0_idx)                                               + &
      loc_cal_parms(bl1_idx) * del(1:npoints)                              + &
      loc_cal_parms(bl2_idx) * del(1:npoints)*del(1:npoints)               + &
      loc_cal_parms(bl3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints)

    RETURN
  END SUBROUTINE spectrum_solar

  !---------------------------------------------------------------------------
  ! In this routine, the parameter num_iterations serves a dual role.
  ! On input, it is that value of the maximum number of iterations per fit.
  ! On output, it will set to the total number of iterations for all fits.
  subroutine wavecal_fit ( &
      wavelengths, spectrum, weights, residuals, &
      num_wavelengths, avg_wavelength, &
      loc_cal_parms, min_cal_parms, max_cal_parms, num_cal_parms, &
      n_refits, sdev_factor, &
      is_bad_pixel, num_iterations, chisqr, exit_value, &
      errstat)

    use OMSAO_parameters_module, only: downweight, r8_missval
    use OMSAO_variables_module, only: tol, epsrel, epsabs, epsx

    implicit none
    integer (kind=i4), intent(in) :: num_wavelengths, num_cal_parms
    real (kind=r8), intent(in), dimension(num_wavelengths) :: &
      wavelengths, spectrum
    real (kind=r8), intent(inout), dimension(num_wavelengths) :: &
      weights, residuals
    real (kind=r8), intent(in) :: avg_wavelength, sdev_factor
    real (kind=r8), dimension(num_cal_parms), intent(inout) ::loc_cal_parms
    real (kind=r8), dimension(num_cal_parms), intent(in) :: &
      min_cal_parms, max_cal_parms
    integer (kind=i4), intent (in) :: n_refits
    integer (kind=i4), intent(inout) :: num_iterations
    logical, intent(out) :: is_bad_pixel
    real (kind=r8), intent(out) :: chisqr
    integer (kind=i4), intent(out) :: exit_value
    integer, intent(inout) :: errstat
    !
    INTEGER (KIND=i4)  :: i, idx, n_nozero_wgt, num_iterations_per_fit, fit_loop_limit
    REAL    (KIND=r8)  :: mean, sdev, loclim
    type(optimizer_type) :: opt
    real (kind=r8), dimension(num_cal_parms) :: fitvar, lobnd, upbnd
    integer (kind=i4), dimension(num_cal_parms) :: param_mask
    real (kind=r8), dimension (num_wavelengths) :: fitres
    integer :: num_fitvar
    character (len=256) :: log_msg

    if (errstat /= 0) return

    ! pack the parameters
    num_fitvar = 0
    do i = 1, num_cal_parms
      if (min_cal_parms(i) >= max_cal_parms(i)) cycle
      num_fitvar = num_fitvar + 1
      param_mask(num_fitvar) = i
      fitvar(num_fitvar) = loc_cal_parms(i)
      lobnd (num_fitvar) = min_cal_parms(i)
      upbnd (num_fitvar) = max_cal_parms(i)
    enddo

    write (log_msg, *)'wavecal_fit:  num_fitvar=',num_fitvar
    call tell_log (2, log_msg)

    num_iterations_per_fit = num_iterations
    num_iterations = 0

    ! --------------------------------------------------------------------
    ! Check whether we enough spectral points to carry out the fitting. If
    ! not, call it a bad pixel and return.
    ! --------------------------------------------------------------------
    if (num_fitvar >= num_wavelengths ) then
      is_bad_pixel = .TRUE.
      chisqr = r8_missval
      return
    endif

    ! Set up module local variables
    call allocate_module_variables (wavelengths, spectrum, weights, &
                                    num_wavelengths, errstat)
    if (errstat /= 0) return
    cal_parms(1:num_cal_parms) = loc_cal_parms(1:num_cal_parms)

    private_avg_wavelength = avg_wavelength   !sol_wav_avg = avg_wavelength

    ! ---------------------------------------------------------------------
    ! Attempt to standardize the re-iteration with spectral points excluded
    ! that have fitting residuals larger than a pre-set window. Needs more
    ! thinking before it can replace a simple window determined empirically
    ! from fitting lots of spectra.
    ! ---------------------------------------------------------------------
    call optimizer_open (opt, wavecal_residuals, num_fitvar, errstat, &
                         mode=opt_bounded, tol=tol, epsabs=epsabs, &
                         epsrel=epsrel, epsx=epsx, &
                         param_min = lobnd(1:num_fitvar), &
                         param_max = upbnd(1:num_fitvar), &
                         param_mask = param_mask(1:num_fitvar), &
                         max_num_iterations = num_iterations_per_fit)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "wavecal_fit: optimizer_open failed", errstat)
      goto 666
    endif

    loclim = 0.0_r8

    fit_loop_limit = MAX(n_refits, 0)
    fit_loop: do i=0, fit_loop_limit

      if (i /= 0) then
        where (abs(fitres(1:num_wavelengths)) > loclim )
          cal_weights(1:num_wavelengths) = downweight
        end where
      endif

      call opt%optimize (opt, fitvar(1:num_fitvar), num_fitvar, &
                         fitres(1:num_wavelengths), num_wavelengths, &
                         exit_value)
      num_iterations = num_iterations + opt%num_iterations

      n_nozero_wgt = int ( anint ( sum(cal_weights(1:num_wavelengths)) ) )
      if (n_nozero_wgt > 0) then
        chisqr = sum  ( fitres(1:num_wavelengths)**2 )
      else
        chisqr = r8_missval
      endif

      if (i == 0) then
        mean = sum (fitres(1:num_wavelengths)) / real(n_nozero_wgt, kind=r8)
        sdev = sqrt (sum ((fitres(1:num_wavelengths)-mean)**2)/real(n_nozero_wgt-1, kind=r8))
        loclim  = sdev_factor*sdev

        if ((loclim <= 0.0_r8) &
            .or. (n_nozero_wgt <= num_fitvar)) &
          exit fit_loop
      endif

      if (maxval(abs(fitres(1:num_wavelengths))) <= loclim) exit fit_loop

    enddo fit_loop

    call optimizer_close (opt, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "wavecal_fit: optimizer_close failed", errstat)
      goto 666
    endif

    ! unpack the parameters
    do i = 1, num_fitvar
      idx = param_mask(i)
      loc_cal_parms(idx) = fitvar(i)
    enddo
    weights(1:num_wavelengths) = cal_weights(1:num_wavelengths)
    residuals(1:num_wavelengths) = fitres(1:num_wavelengths)
    ! drop

 666 continue

    call deallocate_module_variables (errstat)
    return

  end subroutine wavecal_fit

end module
