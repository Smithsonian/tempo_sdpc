!
module m_solar_fit

  public solar_fit
  private

contains

  ! *********************** Modification History ********
  ! xliu: 
  ! 1. Add subroutine solar_fit_vary
  ! 2. Print fitting variables at the end of solar_fit
  ! 3. Add indices for voigt function
  ! *****************************************************

  SUBROUTINE solar_fit (error)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,  ONLY: max_fit_pts
    USE OMSAO_indices_module,     ONLY: max_calfit_idx, & 
      shi_idx, squ_idx,&
      wvl_idx, spc_idx, sig_idx, hwr_idx, hwe_idx, vgl_idx, asy_idx, spk_idx,&
      wr0_idx, wr7_idx
    USE OMSAO_variables_module,   ONLY: curr_sol_spec, n_fitvar_sol,     &
         fitvar_sol, mask_fitvar_sol, lo_sunbnd, up_sunbnd, n_irrad_wvl, &
         wincal_wav, solwinfit, fixslitcal, fitwavs, nsolpix,    &
         fitweights, currspec, numwin, which_slit, fitvar_sol_init, &
         scnwrt, lo_sunbnd_init, up_sunbnd_init, wavcal, &
         sol_wav_avg, rmask_fitvar_sol, poly_order, correct_lamda
    USE OMSAO_errstat_module
    use m_cal_fit_one

    IMPLICIT NONE

    ! ================
    ! Output variables
    ! ================
    LOGICAL,            INTENT (OUT)             :: error

    ! ===============
    ! Local variables
    ! ===============
    INTEGER         :: i,j, iwin, fidx, lidx, n_fit_pts,  &
                       solfit_exval, ll, lu
    REAL (KIND=dp), DIMENSION(8) :: polycoeffs
    REAL (KIND=dp), DIMENSION(max_fit_pts) :: polyx
    REAL (KIND=dp), DIMENSION (n_irrad_wvl)      :: allwaves, del
    REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_varstd
    INTEGER, SAVE   :: slit_unit
    LOGICAL, SAVE   :: wrt_to_screen, wrt_to_file, slitcal
    LOGICAL, SAVE   :: first = .TRUE.

    IF (first) THEN
      IF (scnwrt) THEN
        wrt_to_screen = .TRUE.
      ELSE
        wrt_to_screen = .FALSE.
      ENDIF

      fixslitcal = .TRUE.; slitcal = .TRUE.
      wrt_to_file = .FALSE.; slit_unit = 1000
      fitvar_sol = fitvar_sol_init
      lo_sunbnd  = lo_sunbnd_init; up_sunbnd  = up_sunbnd_init

      ! find the locations of actually used fitting variables
      IF (which_slit == 5) THEN
        fitvar_sol(hwe_idx:asy_idx) = 0_dp
        lo_sunbnd(hwe_idx:asy_idx)  = 0_dp
        up_sunbnd(hwe_idx:asy_idx)  = 0_dp

        fitvar_sol(vgl_idx:spk_idx) = 0_dp
        lo_sunbnd(vgl_idx:spk_idx)  = 0_dp
        up_sunbnd(vgl_idx:spk_idx)  = 0_dp
        fixslitcal = .FALSE.; slitcal = .FALSE.
      ENDIF

      n_fitvar_sol = 0
      DO i = 1, max_calfit_idx
        IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
          n_fitvar_sol =  n_fitvar_sol + 1
          mask_fitvar_sol(n_fitvar_sol) = i
        END IF
      END DO
      first = .FALSE.
    ENDIF

    error = .FALSE.
    allwaves = curr_sol_spec(wvl_idx, 1:n_irrad_wvl)
    fidx = 1
    DO iwin = 1, numwin    

      ! get spectra
      n_fit_pts = nsolpix(iwin)
      lidx = fidx + n_fit_pts - 1
      fitwavs   (1:n_fit_pts) = curr_sol_spec(wvl_idx, fidx:lidx)
      currspec  (1:n_fit_pts) = curr_sol_spec(spc_idx, fidx:lidx)
      fitweights(1:n_fit_pts) = curr_sol_spec(sig_idx, fidx:lidx)

     !DO i = 1, n_fit_pts
     ! print * , i, fitwavs(i), currspec(i), fitweights(i)
     !enddo
      ! UV-2, problem in solar reference above 325 nm
      !IF (reduce_resolution) THEN
      !   IF (fitwavs(2) < 325.1 .AND. fitwavs(n_fit_pts) > 325.1) THEN  
      !      n_fit_pts = MAXVAL(MAXLOC(fitwavs(1:n_fit_pts), &
      !      MASK=(fitwavs(1:n_fit_pts) <= 325.1)))
      !   ENDIF
      !ENDIF

      IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, fitwavs(1), &
           fitwavs(n_fit_pts), nsolpix(iwin)


      fitvar_sol = fitvar_sol_init

      ! Initialization for wavelength registration block this to return_v1
      IF (ANY(rmask_fitvar_sol(wr0_idx:wr7_idx) > 0)) THEN
        DO i = 1, n_fit_pts
           polyx(i) = 1.0d0 * i - 1.0
        ENDDO
        polyx(1:n_fit_pts) = (polyx(1:n_fit_pts) - n_fit_pts / 2.0) / n_fit_pts

        DO i = wr0_idx, wr7_idx
           IF ( fitvar_sol(i) > lo_sunbnd(i) .and. fitvar_sol(i) < up_sunbnd(i)) THEN
              poly_order = i - wr0_idx + 1
           ENDIF
        ENDDO
        ll = 1; lu = n_fit_pts
      !  CALL poly_fit(polyx(1:n_fit_pts), n_fit_pts, fitwavs(1:n_fit_pts), ll,lu, polycoeffs(1:poly_order))

        j = 1
        DO i = wr0_idx, wr7_idx
           IF ( fitvar_sol(i) > lo_sunbnd(i) .and. fitvar_sol(i) < up_sunbnd(i)) THEN
              fitvar_sol(i) = polycoeffs(j)
              lo_sunbnd(i) = -1.0D+99
              up_sunbnd(i) =  1.0D+99
              j = j + 1
           ENDIF
        ENDDO
      ENDIF

      CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
           slitcal, slit_unit, wincal_wav(iwin), &
!           solwinfit(iwin,1:max_calfit_idx, 1:2), solfit_exval)
           tmp_varstd, solfit_exval)
      solwinfit(iwin,1:max_calfit_idx, 1:2)=tmp_varstd

      IF (solfit_exval < 0) THEN
        WRITE(www_lun, *) &
             'Solar_fit: solar calibration not converge for window: ', iwin
        error = .TRUE.; RETURN
      END IF

      ! =================================
      ! Shift and squeeze solar spectrum.
      ! =================================
      ! fitvar_sol is updated in solar_fit_one through common module variables
      IF (wavcal) THEN
        IF (correct_lamda == 1) THEN
            allwaves(fidx:lidx) = (allwaves(fidx:lidx) - fitvar_sol(shi_idx)) / (1.0 + fitvar_sol(squ_idx))     
        ELSE
            allwaves(fidx:lidx) = (allwaves(fidx:lidx) - fitvar_sol(shi_idx) + &
                                   sol_wav_avg * fitvar_sol(squ_idx)) / (1.0 +fitvar_sol(squ_idx))
        ENDIF
     ENDIF

      fidx = lidx + 1
    END DO

    !  fidx = 1; sfidx  =1
    !  DO i = 1, numwin
    !     lidx  = fidx + nsolpix(i) - 1
    !     IF (i == numwin) THEN 
    !        slidx = nsol_ring
    !     ELSE
    !        slidx = MINVAL(MAXLOC(sol_spec_ring(1, 1:nsol_ring), &
    !             MASK=(sol_spec_ring(1, 1:nsol_ring) < curr_sol_spec(wvl_idx, lidx+1))))
    !     ENDIF
    !
    !     finter = MINVAL(MAXLOC(sol_spec_ring(1, 1:nsol_ring), &
    !             MASK=(sol_spec_ring(1, 1:nsol_ring) == curr_sol_spec(wvl_idx, fidx))))
    !     linter = finter + nsolpix(i) - 1
    !
    !     sol_spec_ring(1, finter:linter) = allwaves(fidx:lidx)
    !     IF (finter > sfidx) sol_spec_ring(1, sfidx:finter-1) = (sol_spec_ring(1,sfidx:finter-1)  &
    !          - solwinfit(i, shi_idx, 1)) / (1.0 + solwinfit(i, squ_idx, 1))
    !     IF (linter < slidx) sol_spec_ring(1, linter+1:slidx) = (sol_spec_ring(1, linter+1:slidx) &
    !          - solwinfit(i, shi_idx, 1)) / (1.0 + solwinfit(i, squ_idx, 1))
    !    
    !     fidx = lidx+ 1
    !     sfidx= slidx + 1
    !  ENDDO

    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = allwaves

    RETURN

  END SUBROUTINE solar_fit


end module m_solar_fit
