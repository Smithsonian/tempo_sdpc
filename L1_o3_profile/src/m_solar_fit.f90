!
module m_solar_fit

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen
  USE OMSAO_indices_module,    ONLY : max_calfit_idx, &
      shi_idx, squ_idx, wvl_idx, spc_idx, sig_idx, &
      hwe_idx, asy_idx, hwr_idx, hwl_idx, vgr_idx, vgl_idx, spk_idx, &
      wr0_idx, wr7_idx
  USE OMSAO_variables_module, ONLY : scnwrt, numwin, winlim, currpixchar, &
      n_irrad_wvl,nsolpix, curr_sol_spec, nsol_ring,sol_spec_ring, &
      fitwavs, currspec, fitweights,&
      n_fitvar_sol, fitvar_sol, lo_sunbnd, up_sunbnd, &
      mask_fitvar_sol, rmask_fitvar_sol, sol_wav_avg, &
      fitvar_sol_init, lo_sunbnd_init, up_sunbnd_init, fitvar_sol_saved, &
      which_slit, instrument_sidx, wavcal_sol, wavcal, fixslitcal, slit_fname, poly_order, &
      correct_lambda, xbin_decerr, fitspec_rad !, calscn

  USE OMSAO_errstat_module
  USE m_cal_fit_one, ONLY: calfitone , cal_fit_one
  USE m_fitting_util, ONLY: poly_fit
  IMPLICIT NONE
  LOGICAL, PARAMETER :: calscn=.false.
  INTEGER, PARAMETER, PRIVATE :: slit_unit = 1000
  
  public solar_fit, solar_fit_vary
  private

CONTAINS

  ! *********************** Modification History ********
  ! xliu: 
  ! 1. Add subroutine solar_fit_vary
  ! 2. Print fitting variables at the end of solar_fit
  ! 3. Add indices for voigt function
  ! Jbak
  ! 1. solar_fit and solar_fit_vary is merged here
  ! 2. Add indices for super gaussian and wavelength polynominal fitting 
  ! *****************************************************

  SUBROUTINE solar_fit (error)

  USE OMSAO_variables_module,   ONLY: wincal_wav, solwinfit
  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
  LOGICAL, INTENT (OUT)  :: error

  ! ===============
  ! Local variables
  ! ===============
  INTEGER         :: i,j, iwin, fidx, lidx, n_fit_pts,  &
                    solfit_exval, ll, lu
  REAL (KIND=dp), DIMENSION(8)  :: polycoeffs
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: polyx
  REAL (KIND=dp), DIMENSION (n_irrad_wvl)      :: allwaves
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_varstd

  LOGICAL, SAVE   :: wrt_to_screen, wrt_to_file, slitcal
  LOGICAL, SAVE   :: first = .TRUE.

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'solar_fit'
  
  IF (first) THEN
     wrt_to_screen = calscn
     wrt_to_file = .FALSE.
     fixslitcal = .TRUE.; slitcal = .TRUE.
     ! find the locations of actually used fitting variables
     IF (which_slit >= instrument_sidx) THEN
       fitvar_sol_init(hwe_idx:asy_idx) = 0_dp
       lo_sunbnd_init(hwe_idx:asy_idx)  = 0_dp
       up_sunbnd_init(hwe_idx:asy_idx)  = 0_dp

       fitvar_sol_init(vgl_idx:spk_idx) = 0_dp
       lo_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
       up_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
       fixslitcal = .FALSE.; slitcal = .FALSE.
     ENDIF
       fitvar_sol = fitvar_sol_init
       lo_sunbnd  = lo_sunbnd_init; up_sunbnd  = up_sunbnd_init

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
      fitvar_sol = fitvar_sol_init

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

      ! Initialization for wavelength registration block this to return_v1
      IF (ANY(rmask_fitvar_sol(wr0_idx:wr7_idx) > 0)) THEN
        allocate(polyx(n_fit_pts))
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
        CALL poly_fit(polyx(1:n_fit_pts), n_fit_pts, fitwavs(1:n_fit_pts), ll,lu, polycoeffs(1:poly_order))

        j = 1
        DO i = wr0_idx, wr7_idx
           IF ( fitvar_sol(i) > lo_sunbnd(i) .and. fitvar_sol(i) < up_sunbnd(i)) THEN
              fitvar_sol(i) = polycoeffs(j)
              lo_sunbnd(i) = -1.0D+99
              up_sunbnd(i) =  1.0D+99
              j = j + 1
           ENDIF
        ENDDO
        deallocate(polyx)
      ENDIF

      IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, fitwavs(1), &
           fitwavs(n_fit_pts), nsolpix(iwin)

      CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
           slitcal, slit_unit, wincal_wav(iwin), &
           tmp_varstd, solfit_exval)
      solwinfit(iwin,1:max_calfit_idx, 1:2)=tmp_varstd
      fitspec_rad(fidx:lidx) = calfitone(1:n_fit_pts)
      IF (solfit_exval < 0) THEN
        WRITE(www_lun, *) &
             'Solar_fit: solar calibration not converge for window: ', iwin
        error = .TRUE.; RETURN
      END IF
      
      ! =================================
      ! Shift and squeeze solar spectrum.
      ! =================================
      ! fitvar_sol is updated in solar_fit_one through common module variables
      IF (wavcal_sol) THEN
        IF (correct_lambda == 1) THEN
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

  SUBROUTINE solar_fit_vary (error )

    USE OMSAO_variables_module,  ONLY :  &
         slitwav_sol, nslit_sol, slit_fit_pts, n_slit_step, smooth_slit, &
         slit_redo, solslitfit, nslit, slitwav, slitfit
    USE m_ezspline_interpolation, only: interpolation
    USE m_subtract_poly, only: subtract_poly_meas
    IMPLICIT NONE

    ! ================
    ! Input/OUTPUT variables
    ! ================
    LOGICAL, INTENT (OUT)  :: error

    ! ===============
    ! Local variables
    ! ===============
    REAL (KIND=dp)                          :: tmpwave
    REAL (KIND=dp), DIMENSION (n_irrad_wvl) :: allwaves, locshi, locsqu,locspec
    REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: tmp_varstd
    INTEGER :: npoints, i, iwin, fidx, lidx, errstat = pge_errstat_ok, &
         islit, fpos, lpos, fslit, lslit, ios, finter, &
         linter, npoly, solfit_exval
    CHARACTER(LEN=maxchlen)                       :: tmpchar, fname
    LOGICAL :: calfname_exist = .FALSE.

    ! Save variables
    LOGICAL, SAVE :: wrt_to_screen, wrt_to_file, slitcal, first=.true.
    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=14), PARAMETER :: modulename = 'solar_fit_vary'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    IF (first) THEN 
      wrt_to_screen = calscn
      slitcal=.TRUE.
      fixslitcal = .TRUE.
      wrt_to_file = .FALSE.

      ! find the locations of actually used fitting variables
      IF (which_slit >= instrument_sidx ) THEN 
        fitvar_sol_init(hwe_idx:asy_idx) = 0_dp
        lo_sunbnd_init(hwe_idx:asy_idx)  = 0_dp
        up_sunbnd_init(hwe_idx:asy_idx)  = 0_dp

        fitvar_sol_init(vgl_idx:spk_idx) = 0_dp
        lo_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
        up_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
        fixslitcal = .FALSE.; slitcal = .FALSE.
      ENDIF

      fitvar_sol_saved = fitvar_sol_init 
      fitvar_sol = fitvar_sol_init
      lo_sunbnd  = lo_sunbnd_init; up_sunbnd  = up_sunbnd_init
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

    ! Determine if file exists or not
    fname = TRIM(ADJUSTL(slit_fname)) // currpixchar // '.dat'
    INQUIRE (FILE=TRIM(ADJUSTL(fname)), EXIST=calfname_exist)
    !print * , fname, calfname_exist, slit_redo
    slit_redo = .true.
    IF (slit_redo .OR. .NOT. calfname_exist) THEN

      islit  = 0              ! number of sucessful calibrations      
      fidx = 1                ! first pixel
      DO iwin = 1, numwin    
       
        IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') &
             'win = ', iwin, winlim(iwin,1), &
             winlim(iwin,2), nsolpix(iwin)

        lidx = fidx + nsolpix(iwin) - 1   ! ending index in this window

        ! do each fit using points (fpos:lpos)
        fpos = fidx

        DO WHILE (fpos < lidx)

          DO i = fpos+1, fpos + slit_fit_pts - 1
            ! until end of window or sudden gap of > 2.0 um
            IF (i > lidx) EXIT
            IF (allwaves(i) - allwaves(i-1) > 2.0) EXIT
          ENDDO
          lpos = i - 1
          npoints = lpos - fpos + 1

          IF (npoints < slit_fit_pts / 2) THEN
            ! Either until the end with not enough points or gap behind
            fpos = lpos + 1
            CYCLE
          ENDIF

          fitwavs   (1:npoints) = curr_sol_spec(wvl_idx, fpos:lpos)
          currspec  (1:npoints) = curr_sol_spec(spc_idx, fpos:lpos)
          fitweights(1:npoints) = curr_sol_spec(sig_idx, fpos:lpos)
          fitvar_sol =  fitvar_sol_saved
          CALL cal_fit_one (npoints, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
               slitcal, slit_unit, tmpwave, tmp_varstd, solfit_exval)
          fitspec_rad(fpos:lpos) = calfitone(1:npoints)
          IF (solfit_exval > 0) THEN
            islit = islit + 1
            slitwav_sol(islit) = tmpwave
            solslitfit(islit, 1:max_calfit_idx, 1:2) = &
                tmp_varstd(1:max_calfit_idx, 1:2)
          END IF

          fpos = fpos + n_slit_step          
        END DO

        fidx = lidx + 1       
      END DO
      nslit_sol = islit 

    ELSE  
      READ (slit_unit, '(A)') tmpchar
      islit = 0

      DO 
        islit = islit + 1
        IF (which_slit == 2) THEN
          READ (slit_unit, *, IOSTAT = ios) slitwav_sol(islit), &
               solslitfit(islit, vgl_idx, 1), solslitfit(islit, vgr_idx, 1), &
               solslitfit(islit, hwl_idx, 1), solslitfit(islit, hwr_idx, 1), &
               solslitfit(islit, shi_idx, 1), solslitfit(islit, squ_idx, 1)
        ELSE
          READ (slit_unit, *, IOSTAT = ios) slitwav_sol(islit), &
           solslitfit(islit, hwe_idx, 1), solslitfit(islit, asy_idx, 1), &
           solslitfit(islit, shi_idx, 1), solslitfit(islit, squ_idx, 1), solslitfit(islit, spk_idx, 1)
        END IF

        IF (ios < 0) THEN
          islit = islit - 1
          nslit_sol = islit
          EXIT  ! end of file
        ELSE IF (ios > 0) THEN
          WRITE(www_lun, *) modulename, ': read slit file error'
          error = .TRUE. 
          RETURN 
        END IF
      END DO

    END IF
    !CLOSE(slit_unit)

    ! Disable smoothing slit funciton (only smooth shift)
    IF (smooth_slit ) THEN

      fidx = 1
      DO iwin = 1, numwin

        lidx = fidx + nsolpix(iwin) - 1     
        fslit = MINVAL(MINLOC(slitwav_sol(1:nslit_sol), &
             MASK=(slitwav_sol(1:nslit_sol) >= winlim(iwin, 1))))
        lslit = MINVAL(MAXLOC(slitwav_sol(1:nslit_sol), &
             MASK=(slitwav_sol(1:nslit_sol) <= winlim(iwin, 2))))

        IF (lslit > fslit + 3) THEN           ! use a 4-order polynomial
          poly_order = 3

          DO i = hwe_idx, spk_idx
            IF (solslitfit(fslit,i,1)==0. .AND. &
                 solslitfit(lslit, i, 1)==0.0) CYCLE
            IF (wavcal_sol .AND. (i == shi_idx .OR. i == squ_idx)) CYCLE

            npoly = lslit - fslit + 1
            locspec(1:npoly) = solslitfit(fslit:lslit, i, 1)
            CALL subtract_poly_meas (slitwav_sol(fslit:lslit), npoly, &
                 locspec(1:npoly), 1, npoly)
            solslitfit(fslit:lslit, i,1) = solslitfit(fslit:lslit, i, 1)&
                 - locspec(1:npoly) 
          END DO
        ELSE                                 ! use average
          DO i = hwe_idx, spk_idx
            solslitfit(fslit:lslit, i, 1) = SUM(solslitfit(fslit:lslit, i, 1)) / &
                 (lslit - fslit + 1.0)
          END DO
        ENDIF

        fidx = lidx + 1
      END DO
    END IF

    ! to be used in voigt.f90 or gauss.f90
    nslit = nslit_sol
    slitwav = slitwav_sol
    slitfit = solslitfit


    ! squeeze and shift wavelength position
    IF (wavcal_sol .OR. .NOT. wavcal) RETURN     ! done in next iteration

    fslit = 1
    fidx = 1

    DO iwin = 1, numwin

      lidx = fidx + nsolpix(iwin) - 1     
      lslit = MINVAL(MAXLOC(slitwav_sol, MASK=(slitwav_sol(1:nslit_sol)&
           <= winlim(iwin, 2))))

      IF (lslit < fslit + 3) THEN
        locshi(fidx:lidx) = solslitfit(fslit, shi_idx, 1)
        locsqu(fidx:lidx) = solslitfit(fslit, squ_idx, 1)
      ELSE
        finter = MINVAL(MINLOC(allwaves, MASK=(allwaves > slitwav_sol(fslit))))
        linter = MINVAL(MAXLOC(allwaves, MASK=(allwaves < slitwav_sol(lslit))))
        IF (finter == 0) CYCLE

        CALL interpolation (lslit-fslit+1,  slitwav_sol(fslit:lslit), &
             solslitfit(fslit:lslit, shi_idx, 1),   linter-finter+1, &
             allwaves(finter:linter), locshi(finter:linter),errstat )
        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
          STOP 1
        END IF

        CALL interpolation (lslit-fslit+1, slitwav_sol(fslit:lslit),  &
             solslitfit(fslit:lslit, squ_idx, 1),  linter-finter+1,  &
             allwaves(finter:linter), locsqu(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
          STOP 1
        END IF

        IF (finter > fidx ) THEN
          locshi(fidx:finter-1)=solslitfit(fslit, shi_idx, 1)
          locsqu(fidx:finter-1)=solslitfit(fslit, squ_idx, 1)
        END IF

        IF (linter < lidx)  THEN
          locshi(linter+1:lidx)  = solslitfit(lslit, shi_idx, 1)
          locsqu(linter+1:lidx) =  solslitfit(lslit, squ_idx, 1)
        END IF
      END IF
      
      IF (correct_lambda == 1) THEN 
          allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx) ) &
                                 / ( 1.0 + locsqu(fidx:lidx))
        ELSE
          allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx)  + &
                                 sol_wav_avg*locsqu(fidx:lidx)/ ( 1.0 + locsqu(fidx:lidx)) )
      ENDIF

      fidx = lidx + 1
      fslit = lslit + 1     
    END DO

    !fidx = 1 ; sfidx = 1
    !DO i = 1, numwin
    !  lidx  = fidx + nsolpix(i) - 1
    !  IF (i == numwin) THEN 
    !    slidx = nsol_ring
    !  ELSE
    !    slidx = MINVAL(MAXLOC(sol_spec_ring(1, 1:nsol_ring), &
    !         MASK=(sol_spec_ring(1, 1:nsol_ring) < &
    !         curr_sol_spec(wvl_idx, lidx+1))))
    !  ENDIF

    !  finter = MINVAL(MAXLOC(sol_spec_ring(1, 1:nsol_ring), &
    !       MASK=(sol_spec_ring(1, 1:nsol_ring) == &
    !       curr_sol_spec(wvl_idx, fidx))))
    !  linter = finter + nsolpix(i) - 1

    !  sol_spec_ring(1, finter:linter) = allwaves(fidx:lidx)
    !  IF (finter > sfidx) sol_spec_ring(1, sfidx:finter-1) = &
    !       (sol_spec_ring(1,sfidx:finter-1) - &
    !       locshi(fidx)) / (1.0 + locsqu(fidx))
    !  IF (linter < slidx) sol_spec_ring(1, linter+1:slidx) = &
    !       (sol_spec_ring(1, linter+1:slidx) - &
    !       locshi(lidx)) / (1.0 + locsqu(lidx))

    !  fidx = lidx+ 1
    !  sfidx= slidx + 1
    !ENDDO

    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = allwaves

    RETURN
  END SUBROUTINE solar_fit_vary

end module m_solar_fit
