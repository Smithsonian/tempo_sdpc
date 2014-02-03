MODULE spectra
  use OMSAO_precision_module, only : i4, r8

  use OMSAO_indices_module, only : max_rs_idx, n_max_fitpars
  logical, dimension(n_max_fitpars) :: param_frozen_at_zero
  logical, dimension(max_rs_idx) :: database_j_is_zero

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

CONTAINS

SUBROUTINE spectrum_earthshine (npts, rad_wav_avg, locwvl, fit, rad_fitvar)

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
    max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
    lbe_idx, ad2_idx, mxs_idx, wvl_idx, spc_idx,                   &
    bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, sc2_idx, &
    sc3_idx, sin_idx, shi_idx, squ_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts, downweight
  USE OMSAO_variables_module,  ONLY: &
    n_database_wvl, curr_sol_spec, fitweights, &
    database, curr_xtrack_pixnum
  use ctrlvars, only: yn_radiance_reference, yn_spectrum_norm, &
    yn_doas, yn_newshift, yn_solar_comp
  USE OMSAO_prefitcol_module,  ONLY:  apply_prefit_values
  USE OMSAO_omidata_module,      ONLY: omi_solcal_pars
  USE cache_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_errstat_module
  USE OMSAO_solcomp_module, ONLY: soco_compute
  USE sao_pge_utils, ONLY: array_locate_r8, interpolation, array_sort_r8

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
  sunpos_ss  (1:n_sunpos) = curr_sol_spec(wvl_idx,1:n_sunpos)
  sunspec_loc(1:n_sunpos) = curr_sol_spec(spc_idx,1:n_sunpos)

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

  IF ( j1 > 1    )  fitweights(1:j1-1)    = downweight
  IF ( j2 < npts )  fitweights(j2+1:npts) = downweight

  fit = 0.0_r8

  ! ==================================================================
  ! For BOAS or any wavelength calibration, we have the following line
  ! ==================================================================

  fit(j1:j2) = rad_fitvar(sin_idx) * sunspec_ss(j1:j2)

  !     DOAS here - the spectrum to be fitted needs to be re-defined:
  IF ( yn_doas ) THEN

    i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

    fit(j1:j2) = &
      ! For DOAS, rad_fitvar(SIN_IDX) should == 1., and not be varied
      rad_fitvar(sin_idx) * LOG ( sunspec_ss(j1:j2) ) + &
      ! Ring adjustment
      rad_fitvar(i) * (database(ring_idx,j1:j2) / sunspec_ss (j1:j2))

    DO j = 1, max_rs_idx
      IF ( j /= solar_idx .AND. j /= ring_idx ) THEN
        if (database_j_is_zero(j)) cycle
        i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
        fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j,j1:j2)
      END IF
    END DO

  ELSE
    sumexp => tmpspace
    sumexp(j1:j2) = 0.0_r8
    fit_final_add_on(j1:j2) = 0.0_r8
    DO j = 1, max_rs_idx
      IF ( j.eq.solar_idx ) CYCLE
      if (database_j_is_zero(j)) cycle
      database_j(j1:j2) = database(j, j1:j2)
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
    n_database_wvl, curr_sol_spec, fitweights, &
    database, curr_xtrack_pixnum
  use ctrlvars, only: yn_radiance_reference, yn_spectrum_norm, yn_doas, &
    yn_solar_comp, yn_newshift
  USE OMSAO_prefitcol_module,  ONLY:  apply_prefit_values
  USE OMSAO_omidata_module,      ONLY: omi_solcal_pars
  USE cache_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_errstat_module
  USE OMSAO_solcomp_module, ONLY: soco_compute
  USE sao_pge_utils, ONLY: array_locate_r8, interpolation, array_sort_r8

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
  sunpos_ss  (1:n_sunpos) = curr_sol_spec(wvl_idx,1:n_sunpos)
  sunspec_loc(1:n_sunpos) = curr_sol_spec(spc_idx,1:n_sunpos)

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

  IF ( j1 > 1    )  fitweights(1:j1-1)       = downweight
  IF ( j2 < npts )  fitweights(j2+1:npts) = downweight

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

    fit(j1:j2) = &
      ! For DOAS, rad_fitvar(SIN_IDX) should == 1., and not be varied
      rad_fitvar(sin_idx) * LOG ( sunspec_ss(j1:j2) ) + &
      ! Ring adjustment
      rad_fitvar(i) * (database(ring_idx,j1:j2) / sunspec_ss (j1:j2))

    DO j = 1, max_rs_idx
      IF ( j /= solar_idx .AND. j /= ring_idx  .AND. &
        j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
        if (database_j_is_zero(j)) cycle
        i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
        if (.not. param_frozen_at_zero(i)) then
          fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j,j1:j2)
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
          fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j,j1:j2)
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
            tmpexp(j1:j2) = rad_fitvar(i)*database(j,j1:j2) *  &
              (1.0_r8 + rad_fitvar(k1)*del(j1:j2) + &
               rad_fitvar(k2)*del(j1:j2)*del(j1:j2))
          ELSE
            tmpexp(j1:j2) = rad_fitvar(i)*database(j,j1:j2)
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
          fit(j1:j2) = fit(j1:j2) + rad_fitvar(i) * database(j,j1:j2)
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
