!
module m_radiance_fit

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, max_fit_pts
  USE OMSAO_indices_module,    ONLY : max_calfit_idx, &
      shi_idx, squ_idx, wvl_idx, spc_idx, sig_idx, &
      hwe_idx, asy_idx, hwr_idx, hwl_idx, vgr_idx, vgl_idx, spk_idx, &
      wr0_idx, wr7_idx
  USE OMSAO_variables_module, ONLY : scnwrt, numwin, winlim, currpixchar, &
      n_rad_wvl, nradpix, curr_rad_spec,  &
      fitwavs, currspec, fitweights,&
      n_fitvar_sol, fitvar_sol, lo_sunbnd, up_sunbnd, &
      mask_fitvar_sol, rmask_fitvar_sol, sol_wav_avg, &
      fitvar_sol_init, lo_sunbnd_init, up_sunbnd_init, fitvar_sol_saved, &
      which_slit, wavcal_sol, wavcal, fixslitcal, rslit_fname, poly_order, &
      correct_lambda, xbin_decerr
  USE OMSAO_errstat_module
  USE m_cal_fit_one

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE ::rslit_unit = 1001

  
  public radiance_fit, radiance_fit_vary
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

  SUBROUTINE radiance_fit (error)

  USE OMSAO_variables_module,   ONLY: wincal_wav, radwinfit

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
  REAL (KIND=dp), DIMENSION(8) :: polycoeffs
  REAL (KIND=dp), DIMENSION(max_fit_pts) :: polyx
  REAL (KIND=dp), DIMENSION (n_rad_wvl)      :: allwaves
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_varstd

  LOGICAL, SAVE   :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  !CHARACTER (LEN=*), PARAMETER :: modulename = 'radiance_fit'

  IF (first) THEN
      wrt_to_screen = .FALSE.

      fixslitcal = .TRUE.; slitcal = .TRUE.
      wrt_to_file = .FALSE.
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
    allwaves = curr_rad_spec(wvl_idx, 1:n_rad_wvl)
    fidx = 1
    DO iwin = 1, numwin    

      ! get spectra
      n_fit_pts = nradpix(iwin)
      lidx = fidx + n_fit_pts - 1
      fitwavs   (1:n_fit_pts) = curr_rad_spec(wvl_idx, fidx:lidx)
      currspec  (1:n_fit_pts) = curr_rad_spec(spc_idx, fidx:lidx)
      fitweights(1:n_fit_pts) = curr_rad_spec(sig_idx, fidx:lidx)
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

      IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, fitwavs(1), &
           fitwavs(n_fit_pts), nradpix(iwin)

      CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
           slitcal, rslit_unit, wincal_wav(iwin), &
           tmp_varstd, solfit_exval)
      radwinfit(iwin,1:max_calfit_idx, 1:2)=tmp_varstd

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

    curr_rad_spec(wvl_idx, 1:n_rad_wvl) = allwaves

    RETURN

  END SUBROUTINE radiance_fit

  SUBROUTINE radiance_fit_vary (n_rad_wvl, curr_rad_spec, error )

    USE OMSAO_variables_module,  ONLY : slitwav_rad, nslit_rad,  &
         slit_fit_pts, n_slit_step, smooth_slit,  slit_redo, rslit_fname,&
         nradpix, radslitfit, &        
         nslit, slitwav, slitfit

    use m_ezspline_interpolation, only: interpolation
    use m_subtract_poly, only: subtract_poly_meas


    IMPLICIT NONE

    ! ================
    ! Input/OUTPUT variables
    ! ================
    LOGICAL, INTENT (OUT)            :: error
    INTEGER, INTENT (IN)             :: n_rad_wvl
    REAL (KIND=dp), DIMENSION (sig_idx,n_rad_wvl), INTENT (INOUT) :: &
         curr_rad_spec

    ! ===============
    ! Local variables
    ! ===============
    REAL (KIND=dp)                        :: tmpwave
    REAL (KIND=dp), DIMENSION (n_rad_wvl) :: allwaves, locshi, locsqu, locspec
    REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: tmp_fitvar
    INTEGER :: npoints, i, iwin, fidx, lidx, errstat = pge_errstat_ok, &
         islit, fpos, lpos, fslit, lslit, ios, finter, linter, &
         npoly, solfit_exval
    CHARACTER(LEN=maxchlen)                    :: tmpchar, fname
    LOGICAL        :: calfname_exist
    LOGICAL, SAVE  :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=12), PARAMETER :: modulename = 'rad_fit_vary'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    IF (first) THEN
      slitcal=.TRUE.
      fixslitcal = .TRUE.
      wrt_to_screen = .FALSE.
      wrt_to_file = .FALSE.
      
      fitvar_sol = fitvar_sol_init

      ! need to restore back the initial fitting variables and bounds
      IF (wavcal_sol) THEN
        lo_sunbnd  = lo_sunbnd_init
        up_sunbnd  = up_sunbnd_init
        fitvar_sol_saved = fitvar_sol_init     
      END IF

      ! find the number of actual used fitting variables
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
    allwaves = curr_rad_spec(wvl_idx, 1 : n_rad_wvl) 

    ! Determine if file exists or not
    fname = TRIM(ADJUSTL(rslit_fname)) // currpixchar // '.dat'
    INQUIRE (FILE=TRIM(ADJUSTL(fname)), EXIST=calfname_exist)


    IF (slit_redo .OR. .NOT. calfname_exist) THEN

      islit  = 0               ! number of sucessful calibrations      
      fidx = 1                 ! first pixel
      DO iwin = 1, numwin     

        IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, winlim(iwin,1), &
             winlim(iwin,2), nradpix(iwin)
        lidx = fidx + nradpix(iwin) - 1

        ! do each fit using points (fpos:lpos)
        fpos = fidx

        DO WHILE (fpos < lidx)

          DO i = fpos+1, fpos + slit_fit_pts - 1
            ! until end of window or sudden gap of > 2.0 um
            IF (i > lidx) EXIT
            IF (allwaves(i) - allwaves(i-1) > 2.0) EXIT
          ENDDO
          lpos = i - 1; npoints = lpos - fpos + 1

          IF (npoints < slit_fit_pts / 2) THEN
            ! Either until the end with not enough points or gap behind
            fpos = lpos + 1; CYCLE
          ENDIF

          npoints = lpos - fpos + 1
          fitwavs   (1:npoints) = curr_rad_spec(wvl_idx, fpos:lpos)
          currspec  (1:npoints) = curr_rad_spec(spc_idx, fpos:lpos)
          fitweights(1:npoints) = curr_rad_spec(sig_idx, fpos:lpos)
          fitvar_sol =  fitvar_sol_saved
          CALL cal_fit_one (npoints, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
               slitcal, rslit_unit, tmpwave, tmp_fitvar, solfit_exval)

          IF (solfit_exval > 0) THEN
            islit = islit + 1;   slitwav_rad(islit) = tmpwave
            radslitfit(islit, 1:max_calfit_idx, 1:2) = &
                 tmp_fitvar(1:max_calfit_idx, 1:2)
          END IF

          fpos = fpos + n_slit_step           
        END DO

        fidx = lidx + 1
      END DO

      nslit_rad = islit 
    ELSE

      READ (rslit_unit, '(A)') tmpchar

      islit = 0
      DO 
        islit = islit + 1
        IF (which_slit == 2) THEN
          READ (rslit_unit, *, IOSTAT = ios) slitwav_rad(islit), &
               radslitfit(islit, vgl_idx, 1), radslitfit(islit, vgr_idx, 1), &
               radslitfit(islit, spk_idx, 1), radslitfit(islit, hwr_idx, 1), &
               radslitfit(islit, shi_idx, 1), radslitfit(islit, squ_idx, 1)
        ELSE
          READ (rslit_unit, *, IOSTAT = ios) slitwav_rad(islit), &
               radslitfit(islit, hwe_idx, 1), radslitfit(islit, asy_idx, 1), &
               radslitfit(islit, shi_idx, 1), radslitfit(islit, squ_idx, 1), radslitfit(islit, spk_idx, 1)
        END IF

        IF (ios < 0) THEN
          islit = islit - 1; nslit_rad = islit; EXIT  ! end of file
        ELSE IF (ios > 0) THEN
          WRITE(www_lun, *) modulename, ': read slit file error'
          error = .TRUE. ; RETURN
        END IF

      END DO

    END IF
    !CLOSE(rslit_unit)

    IF (smooth_slit) THEN

      fidx = 1
      DO iwin = 1, numwin

        lidx = fidx + nradpix(iwin) - 1
        fslit = MINVAL(MINLOC(slitwav_rad(1:nslit_rad), &
             MASK=(slitwav_rad(1:nslit_rad) >= winlim(iwin, 1))))
        lslit = MINVAL(MAXLOC(slitwav_rad(1:nslit_rad), &
             MASK=(slitwav_rad(1:nslit_rad) <= winlim(iwin, 2))))

        ! use a 4-order polynomial
        IF (lslit > fslit + 3) THEN
          poly_order = 4
          DO i = hwe_idx, spk_idx
            IF (radslitfit(fslit,i,1)==0. .AND. radslitfit(lslit,i,1)==0.0) CYCLE
            IF (wavcal_sol .AND. (i == shi_idx .OR. i == squ_idx)) CYCLE

            npoly = lslit - fslit + 1
            locspec(1:npoly) = radslitfit(fslit:lslit, i, 1)
            CALL subtract_poly_meas ( slitwav_rad(fslit:lslit), npoly, &
                 locspec(1:npoly), 1, npoly )
            radslitfit(fslit:lslit, i, 1) = radslitfit(fslit:lslit, i, 1) &
                 - locspec(1:npoly) 

            !DO j = 3, npoly - 2
            !  locspec(j) = (radslitfit(fslit+j-1, i, 1)*6_dp + &
            !    radslitfit(fslit+j-2, i, 1)*4_dp &
            !  + radslitfit(fslit+j-3, i, 1) + radslitfit(fslit+j, i, 1)*4_dp &
            !  + radslitfit(fslit+j+1, i, 1)) / 15_dp
            !END DO
            radslitfit(fslit:lslit, i, 1) =  locspec(1:npoly)

          END DO
        ELSE
          DO i = hwe_idx, spk_idx
            radslitfit(fslit:lslit, i, 1) = SUM(radslitfit(fslit:lslit, i, 1)) / &
                 (lslit - fslit + 1.0)
          END DO
        ENDIF

        fidx = lidx + 1
      END DO
    END IF

    ! to be used in voigt.f90 or gauss.f90
    nslit = nslit_rad; slitwav = slitwav_rad; slitfit = radslitfit

    !OPEN (unit=77, file='rslit_after_smooth.dat')
    !WRITE(77, *) nslit_rad
    !DO i = 1, nslit_rad
    !   WRITE(77, '(3d14.6)') slitwav_rad(i), radslitfit(i, hwe_idx, 1), &
    !        radslitfit(i, shi_idx, 1)
    !END DO
    !CLOSE (77)

    ! squeeze and shift wavelength position

    IF (wavcal_sol .OR. .NOT. wavcal) RETURN     ! do wavelength registration in next iteration

    fslit = 1; fidx = 1; 
    DO iwin = 1, numwin

      lidx = fidx + nradpix(iwin) - 1

      lslit = MINVAL(MAXLOC(slitwav_rad, MASK=(slitwav_rad(1:nslit_rad) &
           <= winlim(iwin, 2))))

      IF (lslit < fslit + 3) THEN
        locshi(fidx:lidx) = radslitfit(fslit, shi_idx, 1)
        locsqu(fidx:lidx) = radslitfit(fslit, squ_idx, 1)
      ELSE

        finter = MINVAL(MINLOC(allwaves, MASK = &
             (allwaves > slitwav_rad(fslit))))
        linter = MINVAL(MAXLOC(allwaves, MASK = &
             (allwaves < slitwav_rad(lslit))))
        CALL interpolation (lslit-fslit+1, slitwav_rad(fslit:lslit), &
             radslitfit(fslit:lslit, shi_idx, 1), linter-finter+1, &
             allwaves(finter:linter), locshi(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( lslit-fslit+1, slitwav_rad(fslit:lslit), &
             radslitfit(fslit:lslit, squ_idx, 1), linter-finter+1, &
             allwaves(finter:linter), locsqu(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        IF (finter > fidx ) THEN
          locshi(fidx:finter-1)=radslitfit(fslit, shi_idx, 1)
          locsqu(fidx:finter-1)=radslitfit(fslit, squ_idx, 1)
        END IF

        IF (linter < lidx)  THEN
          locshi(linter+1:lidx)  = radslitfit(lslit, shi_idx, 1)
          locsqu(linter+1:lidx) =  radslitfit(lslit, squ_idx, 1)
        END IF
      END IF

      IF (correct_lambda == 1) THEN
         allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx)) &
               / (1.0 + locsqu(fidx:lidx))
     ELSE
         allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx) +  sol_wav_avg * locsqu(fidx:lidx) ) &
               / (1.0 + locsqu(fidx:lidx))
     ENDIF
      fidx = lidx + 1; fslit = lslit + 1     
    END DO

    curr_rad_spec(wvl_idx, 1:n_rad_wvl) = allwaves

    RETURN
  END SUBROUTINE radiance_fit_vary
   
end module m_radiance_fit
