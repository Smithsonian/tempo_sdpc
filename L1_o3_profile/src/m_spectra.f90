!
module m_spectra

  public spectrum_solar, spectrum_earthshine
  private

contains

  ! *********************** Modification History **************
  ! xliu:
  ! 1. Replace call to asym_gauss to voigt_gauss
  ! 2. Add variable vgr_idx, vgl_idx, hwr_idx, hwl_idx
  ! Jbak:
  ! 1. add super_gauss_module
  ! **********************************************************

  SUBROUTINE spectrum_solar ( &
        npoints, nfitvar, sol_wav_avg, locwvl, fit, fitvar)

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, solar_idx, shi_idx, squ_idx, &
        sin_idx, &
        bl0_idx, bl7_idx, sc0_idx, sc7_idx, wr0_idx, wr7_idx
    USE OMSAO_variables_module,  ONLY: n_refspec_pts, refspec_orig_data, &
         fitwavs, fitvar_sol, mask_fitvar_sol,  &
         correct_lambda , rmask_fitvar_sol
    USE OMSAO_errstat_module
    USE OMSAO_slitfunction_module
    USE m_convol, ONLY: simple_convol
    USE m_ezspline_interpolation, only: interpolation

    IMPLICIT NONE


    INTEGER,                            INTENT (INOUT) :: npoints, nfitvar
    REAL (KIND=dp),                     INTENT (INOUT) :: sol_wav_avg
    REAL (KIND=dp), DIMENSION (nfitvar),INTENT (INOUT) :: fitvar
    REAL (KIND=dp), DIMENSION (npoints),INTENT (INOUT) :: locwvl, fit

    REAL (KIND=dp), DIMENSION (npoints) :: del, delx, sunspec_ss , &
        tempspec, tempwave, deli

    ! =======================================
    ! Variable declarations for IMPLICIT NONE
    ! =======================================
    INTEGER :: i, npts, errstat, ref_fidx, ref_lidx

    ! =======================================
    ! Shorthands for solar reference spectrum
    ! =======================================
    REAL (KIND=dp), DIMENSION (n_refspec_pts(solar_idx)) :: kppos, kpspec, kpspec_gauss

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=14), PARAMETER :: modulename = 'spectrum_solar'

    errstat = pge_errstat_ok

    fitvar_sol(mask_fitvar_sol(1:nfitvar)) = fitvar(1:nfitvar)


    IF (ANY(rmask_fitvar_sol(wr0_idx:wr7_idx) > 0)) THEN
      tempwave = 0.0d0
      DO i = 1, npoints
         del(i) = 1.0d0 * i
      ENDDO
      del = (del - npoints / 2.0) / npoints
      deli = 1.0d0

      DO i = wr0_idx, wr7_idx
         tempwave = tempwave + fitvar_sol(i) * deli
         deli = deli * del
      ENDDO
      fitwavs(1:npoints) = tempwave
    ENDIF

    locwvl(1:npoints) = fitwavs(1:npoints)
    sol_wav_avg = ( fitwavs(1) + fitwavs(npoints)) / 2.0

    npts = n_refspec_pts(solar_idx)
    ref_fidx = MINVAL(MINLOC(refspec_orig_data(solar_idx, 1:npts, wvl_idx), &
         MASK = (refspec_orig_data(solar_idx, 1:npts, wvl_idx) >= &
         fitwavs(1) - 2.0)))

    ref_lidx = MINVAL(MAXLOC(refspec_orig_data(solar_idx, 1:npts, wvl_idx), &
         MASK=(refspec_orig_data(solar_idx, 1:npts, wvl_idx) <= &
         fitwavs(npoints)+2.0)))

    npts = ref_lidx - ref_fidx + 1
    kppos (1:npts) = refspec_orig_data(solar_idx, ref_fidx:ref_lidx, wvl_idx)
    kpspec(1:npts) = refspec_orig_data(solar_idx, ref_fidx:ref_lidx, spc_idx)
    !     Spectrum calculation for both fitting and non-fitting cases.

    !     Calculate the spectrum:
    !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
    !     1 + FITVAR(SQU_IDX); do in absolute sense, 
    !     to make it easy to back-convert OMI data.

    IF (correct_lambda == 1) THEN
       kppos(1:npts)  = kppos(1:npts) * (1.0 + fitvar_sol(squ_idx)) + fitvar_sol(shi_idx)
    ELSE
       kppos(1:npts)  = kppos(1:npts) * (1.0 + fitvar_sol(squ_idx)) + fitvar_sol(shi_idx) - sol_wav_avg * fitvar_sol(squ_idx)
    ENDIF
    ! =============================================
    ! Broadening and re-sampling of solar spectrum:
    ! =============================================
   
    CALL simple_convol (kppos(1:npts), kpspec(1:npts), kpspec_gauss(1:npts), npts)
  
    !write(*,'(2f15.7, f8.2)') sum(kpspec(1:npts)), sum(kpspec_gauss(1:npts)), fitvar_sol(hwe_idx), fitvar_sol(spk_idx) !;stop
  

    ! ------------------------------------------------------
    ! Re-sample the solar reference spectrum to the OMI grid
    ! ------------------------------------------------------
    CALL interpolation ( &
         npts, kppos(1:npts), kpspec_gauss(1:npts), &
         npoints, locwvl(1:npoints), sunspec_ss(1:npoints), errstat )
    IF ( errstat > pge_errstat_warning ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
      STOP 1
    END IF

    ! Add up the contributions, with solar intensity as FITVAR_SOL (SIN_IDX),
    ! to include possible linear and Beer's law forms.  Do these as 
    ! linear-Beer's-linear. In order to do DOAS we only need to be careful
    ! to include just linear contributions, since I already high-pass 
    ! filtered them.

    ! -----------
    !  Doing BOAS
    ! -----------
    fit(1:npoints) = fitvar_sol(sin_idx) * sunspec_ss(1:npoints)

    ! ----------------
    ! Add the scaling.
    ! ----------------
    del(1:npoints) = locwvl(1:npoints) - sol_wav_avg
    tempspec = 0.0d0
    delx(1:npoints) = 1.0d0
    DO i = sc0_idx, sc7_idx
     IF (fitvar_sol(i) /= 0.0) THEN
        tempspec = tempspec + fitvar_sol(i) * delx
     ENDIF
     delx = delx * del
    ENDDO
    fit(1:npoints) = fit(1:npoints) * tempspec

    ! ------------------------
    ! Add baseline parameters.
    ! ------------------------
    tempspec = 0.0d0
    delx(1:npoints) = 1.0d0
    DO i = bl0_idx, bl7_idx
     IF (fitvar_sol(i) /= 0.0) THEN
        tempspec = tempspec + fitvar_sol(i) * delx
     ENDIF
     delx = delx * del
    ENDDO
    fit(1:npoints) = fit(1:npoints) + tempspec
     
    RETURN
  END SUBROUTINE spectrum_solar

  SUBROUTINE spectrum_earthshine ( npoints, n_fitvar, rad_wav_avg, &
       locwvl, fit, fitvar, doas )

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: &
         max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
         lbe_idx, ad2_idx, mxs_idx, &
         bl0_idx, sc0_idx, sin_idx, shi_idx, squ_idx, bl7_idx, sc7_idx
    USE OMSAO_variables_module,  ONLY: fitvar_rad, mask_fitvar_rad, &
         n_refwvl, refwvl, database !, curr_rad_spec, curr_sol_spec
    use m_ezspline_interpolation, only: interpolation
    USE OMSAO_errstat_module

    IMPLICIT NONE


    ! ===============
    ! Input variables
    ! ===============
    LOGICAL,  INTENT (IN) :: doas
    INTEGER,  INTENT (IN) :: npoints, n_fitvar
    REAL (KIND=dp),  INTENT (IN) :: rad_wav_avg
    REAL (KIND=dp), DIMENSION (n_fitvar),   INTENT (IN) :: fitvar
    REAL (KIND=dp), DIMENSION (npoints),    INTENT (IN) :: locwvl

    ! ================
    ! Output variables
    ! ================
    REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: fit

    ! ===============
    ! Local variables
    ! ===============
    INTEGER                              :: i, j, idx, errstat
    REAL (KIND=dp), DIMENSION (npoints)  :: del, delx, sunspec_ss, tempspec
    REAL (KIND=dp), DIMENSION (n_refwvl) :: sunpos_ss

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=8), PARAMETER :: modulename = 'spectrum'

    !     Spectrum calculation for both fitting and non-fitting cases.

    !     Calculate the spectrum:
    !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
    !     1 + FITVAR(SQU_IDX); do in absolute sense, 
    !     to make it easy to back-convert OMI data.


    errstat = pge_errstat_ok

    ! -------------------------------------------------------------------------
    ! First, we have to undo the compression of the FITVAR_RAD array. 
    ! This compression is performed in the RADIANCE_FIT subroutine and 
    ! accelerates the fitting process, because ELSUNC has to handle less 
    ! indices. But here we require the original layout,
    ! otherwise the index assingment is screwed.
    ! -------------------------------------------------------------------------
    DO i = 1, n_fitvar
      idx = mask_fitvar_rad(i)
      fitvar_rad(idx) = fitvar(i)
    END DO

    sunpos_ss(1:n_refwvl) = refwvl(1:n_refwvl) * &
         (1.0_dp + fitvar_rad(squ_idx)) + &
         fitvar_rad(shi_idx)

    ! Broadening and re-sampling of solar spectrum:
    !       Re-sample the solar reference spectrum to the radiance grid

    CALL interpolation ( &
         n_refwvl, sunpos_ss(1:n_refwvl), database(solar_idx,1:n_refwvl), &
         npoints, locwvl(1:npoints), sunspec_ss(1:npoints), errstat )

    IF ( errstat > pge_errstat_warning ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
      STOP 1
    END IF

    !     Add up the contributions, with solar intensity as 
    !     FITVAR_RAD(sin_idx), trace species beginning at 
    !     FITVAR_RAD(SQU_IDX+1), to include possible linear and
    !     Beer's law forms.  Do these as linear-Beer's-linear. In order to
    !     do DOAS I only need to be careful to include just linear
    !     contributions, since I already high-pass filtered them.

    ! ==================================================================
    ! For BOAS or any wavelength calibration, we have the following line
    ! ==================================================================
    fit(1:npoints) = fitvar_rad(sin_idx) * sunspec_ss(1:npoints)

    !     DOAS here - the spectrum to be fitted needs to be re-defined:
    IF ( doas ) THEN

      i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

      fit(1:npoints) = &
           ! For DOAS, FITVAR_RAD(SIN_IDX) should == 1., and not be varied
           fitvar_rad(sin_idx) * LOG ( sunspec_ss(1:npoints) ) + &
                                ! Ring adjustment
             fitvar_rad(i) * (database(ring_idx,3:n_refwvl-2) / &
             sunspec_ss (1:npoints))

      DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. j /= ring_idx ) THEN 
          i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
          fit(1:npoints) = fit(1:npoints) + fitvar_rad(i) * &
               database(j,3:n_refwvl-2)
        END IF
      END DO

    ELSE

      ! Initial add-on contributions.
      DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
          i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
          fit(1:npoints) = fit(1:npoints) + fitvar_rad(i) * &
               database(j,3:n_refwvl-2)
        END IF
      END DO
      ! Beer's law contributions.
      DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
          i = max_calfit_idx + (j-1)*mxs_idx + lbe_idx
          fit(1:npoints) = fit(1:npoints) * &
               EXP(-fitvar_rad(i)*database(j,3:n_refwvl-2))

        END IF
      END DO
      ! Final add-on contributions.
      DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
          i = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
          fit(1:npoints) = fit(1:npoints) + fitvar_rad(i) * &
               database(j,3:n_refwvl-2)
        END IF
      END DO
    END IF


   ! Add the scaling.
    del(1:npoints) = locwvl(1:npoints) - rad_wav_avg
    tempspec = 0.0d0
    delx(1:npoints) = 1.0d0
    DO i = sc0_idx, sc7_idx
       IF (fitvar_rad(i) /= 0.0) THEN
          tempspec = tempspec + fitvar_rad(i) * delx
       ENDIF
       delx = delx * del
    ENDDO
    fit(1:npoints) = fit(1:npoints) * tempspec

    ! Add baseline parameters.
    tempspec = 0.0d0
    delx(1:npoints) = 1.0d0
    DO i = bl0_idx, bl7_idx
       IF (fitvar_rad(i) /= 0.0) THEN
          tempspec = tempspec + fitvar_rad(i) * delx
       ENDIF
       delx = delx * del
    ENDDO
    fit(1:npoints) = fit(1:npoints) + tempspec



    RETURN
  END SUBROUTINE spectrum_earthshine


end module m_spectra
