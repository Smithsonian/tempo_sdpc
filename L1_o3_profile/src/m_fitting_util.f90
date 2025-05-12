!
module m_fitting_util

  public  reduce_irrad_resolution, reduce_rad_resolution, &
          rough_spike_detect, uv1_spike_detect , spike_detect_correct, &
          cubic_specfit, poly_specfit, poly_fit
  !spike_detect_correct1, rough_spike_detect1, specfit_func, cubic_func
CONTAINS

  SUBROUTINE reduce_irrad_resolution (spec, qflag, np, nx, which_slit, slwth, &
       samprate, dwav, retswav, retewav, swav, ewav, np_out, pge_error_status)

    USE OMSAO_precision_module 
    USE OMSAO_parameters_module, ONLY: nmax=>max_spec_pts
    USE OMSAO_indices_module,    ONLY: sig_idx!, wvl_idx, spc_idx
    USE OMSAO_variables_module,  ONLY: use_redfixwav, nredfixwav, redfixwav
    USE ozprof_data_module,      ONLY: pos_alb
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: interpolation

    IMPLICIT NONE

    ! ====================
    ! In/Output variables
    ! ====================
    INTEGER, INTENT(IN)                                 :: np, nx, which_slit
    INTEGER, INTENT(OUT)                                :: np_out, pge_error_status
    INTEGER (KIND=i2), DIMENSION(np, nx), INTENT(IN)    :: qflag                   
    REAL (KIND=dp), INTENT(IN)                          :: samprate, slwth, retswav, retewav
    REAL (KIND=dp), INTENT(INOUT)                       :: dwav
    REAL (KIND=dp), INTENT(INOUT)                       :: swav, ewav
    REAL (KIND=dp), DIMENSION(sig_idx, np, nx), INTENT(INOUT) :: spec

    ! ====================
    ! Local variables
    ! ====================
    INTEGER                :: i, j, ix, mslit, nf, nsamp, nsamp1, nslit, errstat, iwin, idx
    INTEGER, DIMENSION(nx) :: nmod
    INTEGER, DIMENSION( 3) :: sidx, eidx, nstep
    REAL (KIND=dp)         :: dlam0, slitsum, redsnr, dx, fwav, lwav
    REAL (KIND=dp), DIMENSION(nmax)      :: slit
    REAL (KIND=dp), DIMENSION(sig_idx, np, nx) :: specmod
    REAL (KIND=dp), DIMENSION(sig_idx, nmax)   :: fnspec

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=28), PARAMETER :: modulename = 'reduce_irradiance_resolution'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    pge_error_status = pge_errstat_ok

    ! If slwth == 0, not allowed, unless that user provides a fixed wavelength grid
    IF (slwth == 0 .AND. .NOT. use_redfixwav) THEN
      WRITE(www_lun, '(2A)') modulename, ': Zero Slit Width Not allowed!!!'
      pge_error_status = pge_errstat_error
    ENDIF

    ! Filter bad measurements, determine wavelength range for all x-track positions
    dlam0 = spec(1, 2, 1) - spec(1, 1, 1)
    ewav = MAXVAL(spec(1, np, :)) 
    DO ix = 1, nx
      j = 0

      DO i = 1, np 
        IF (spec(2, i, ix) > 0. .AND. spec(2, i, ix) <= 4.0E14 .AND. qflag(i, ix) == 0) THEN
          j = j + 1
          specmod(:, j, ix) = spec(:, i, ix)
        ENDIF
      ENDDO
      nmod(ix) = j
      IF (specmod(1, j, ix) < ewav) ewav = specmod(1, j, ix)
    ENDDO
    swav = MAXVAL(specmod(1, 1, :))

    IF (slwth == 0 .AND. use_redfixwav) THEN

      j = 0
      DO i = 1, nredfixwav
        IF (redfixwav(i) >= swav .AND. redfixwav(i) <=ewav) THEN
          j = j + 1
          spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
      ENDDO

    ELSE

      ! Establish fine wavelength scale (common to all across-track positions)
      nf = INT((ewav - swav) / dwav + 1)
      fnspec(1, 1) = swav
      DO i = 2, nf
        fnspec(1, i) = fnspec(1, i - 1) + dwav
      ENDDO

      IF (dwav > dlam0 .OR. dwav <= 0) THEN
        dwav = dlam0
      ENDIF

      IF (which_slit == 1) THEN               ! Symmetric Gaussian (slwth is hw1e)
        ! hw1e * 1.6551 = FWHM
        mslit = NINT(2.62826 * slwth / dwav) ! slit truncation (<0.1%)
        nsamp = INT(slwth * 1.66511 / samprate / dwav)
      ELSE IF (which_slit == 2) THEN          ! Triangular (slwth is FWHM)
        mslit = NINT(slwth / dwav)
        nsamp = INT(slwth / samprate / dwav)
      ENDIF
      nslit = mslit * 2 + 1
      nsamp1 = INT(dlam0 / dwav)
      IF (nsamp1 < 1.0) nsamp1 = 1

      ! Set up slit function
      IF (which_slit == 1) THEN
        slit(1:nslit) = EXP( -((fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth)**2 )
      ELSE
        slit(1:nslit) = 1.0 - ABS(fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth
        WHERE (slit(1:nslit) < 0.0)
          slit(1:nslit) = 0.0
        ENDWHERE
      ENDIF
      slitsum = SUM(slit(1:nslit))
      slit(1:nslit) = slit(1:nslit) / slitsum       ! Normalization
      redsnr = 1.0 / SQRT(slitsum / dlam0 * dwav)   ! Measurement error/noise reduction

      ! Convolve and Sample
      ! More sampling at the two ends of the selected spectral range
      ! Especially for the first and last 4 positions 
      ! This is to avoid extrapolation while keeping as many measurements as possible
      IF (retswav <= fnspec(1, mslit+1)) THEN
        i = mslit + 1 + 3 * nsamp1
      ELSE
        i = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= retswav + dlam0 * 3)))
      ENDIF
      sidx(1) = mslit+1
      eidx(1) = i
      nstep(1) = nsamp1

      IF (retewav >= fnspec(1, nf - mslit)) THEN
        i = nf - mslit - nsamp1 * 3
      ELSE
        i = MAXVAL(MAXLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) <= retewav - dlam0 * 3)))
      ENDIF
      sidx(2) = eidx(1) + nsamp1
      eidx(2) = i-1
      nstep(2) = nsamp
      sidx(3) = i + nsamp1
      eidx(3) = nf - mslit
      nstep(3) = nsamp1

      j = 0
      DO i = 1, nredfixwav
        IF (redfixwav(i) > fnspec(1, mslit + 1) .AND. redfixwav(i) < fnspec(1, nf - mslit)) THEN
          j = j + 1
          spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
      ENDDO
    ENDIF

    IF (use_redfixwav) THEN
      ! Add 3 * 2 extra wavelengths for irradiance and 2 * 2 extra wavelengths for radiance
      ! Add 3 wavelengths at the beginning of a spectra region
      DO i = 1, j
        fwav = MAX(swav, retswav)
        IF (spec(1, i, 1) > fwav) THEN
          IF (i == 1) THEN
            dx = (spec(1, i, 1) - fwav) / 3.
          ELSE
            dx = (spec(1, i, 1) - MAX(spec(1, i-1, 1), fwav)) / 3.
          ENDIF
          IF (dx > dlam0) dx = dlam0
          spec(1, i+3:j+3, 1:nx) = spec(1, i:j, 1:nx)
          spec(1, i, 1:nx)   = spec(1, i+3, 1:nx) - dx * 3.0
          spec(1, i+1, 1:nx)   = spec(1, i+3, 1:nx) - dx * 2.0
          spec(1, i+2, 1:nx)   = spec(1, i+3, 1:nx) - dx * 1.0
          j = j + 3
          EXIT
        ENDIF
      ENDDO

      ! Add 3 wavelengths at the end of a spectra region
      DO i = j, 1, -1
        lwav = MIN(ewav, retewav)
        IF (spec(1, i, 1) < lwav) THEN
          IF (i == j) THEN
            dx = (lwav - spec(1, i, 1)) / 3.
          ELSE
            dx = ( MAX(spec(1, i+1, 1), lwav) - spec(1, i, 1)) / 3.
          ENDIF
          IF (dx > dlam0) dx = dlam0
          spec(1, i+4:j+3, 1:nx) = spec(1, i+1:j, 1:nx)
          spec(1, i+1, 1:nx)   = spec(1, i, 1:nx) + dx * 1.0
          spec(1, i+2, 1:nx)   = spec(1, i, 1:nx) + dx * 2.0
          spec(1, i+3, 1:nx)   = spec(1, i, 1:nx) + dx * 3.0
          j = j + 3
          EXIT
        ENDIF
      ENDDO

      ! Add a wavelength at the wavelength to derive initial cloud fraction
      IF (swav < pos_alb .AND. ewav > pos_alb) THEN
        DO i = j, 1, -1 
          IF (spec(1, i, 1) < pos_alb) THEN
            spec(1, i+1:j, 1:nx) = spec(1, i+2:j+1, 1:nx)
            spec(1, i+1, 1:nx) = pos_alb
            j = j + 1
            EXIT
          ENDIF
        ENDDO
      ENDIF

    ENDIF
    np_out = j

    ! Perform direct interpolation
    IF (slwth == 0 .AND. use_redfixwav) THEN

      DO ix = 1, nx
        DO i = 2, 3
          CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
               np_out, spec(1, 1:np_out, ix),  spec(i, 1:np_out, ix), pge_error_status )
          IF ( pge_error_status > pge_errstat_warning ) THEN
            errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
            RETURN
          END IF
        ENDDO
      ENDDO
    ELSE

      DO ix = 1, nx 
        ! Pre-interpolation
        DO i = 2, 3
          CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
               nf, fnspec(1, 1:nf), fnspec(i, 1:nf), pge_error_status )
          IF ( pge_error_status > pge_errstat_warning ) THEN
            print *, ix, 'interpolation problem'
            errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
            RETURN
          END IF
        ENDDO

        IF (.NOT. use_redfixwav) THEN
          j = 0
          DO iwin = 1, 3
            DO i = sidx(iwin), eidx(iwin), nstep(iwin)
              j = j + 1
              spec(1, j, ix) = fnspec(1, i)
              spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, i-mslit:i+mslit))
              spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, i-mslit:i+mslit)) * redsnr ! Reduce noise
            ENDDO
          ENDDO
        ELSE
          DO j = 1,  np_out
            idx = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= spec(1, j, ix))))
            IF ( ABS(fnspec(1, idx-1) - spec(1, j, ix)) < ABS(fnspec(1, idx) - spec(1, j, ix))) idx = idx - 1
            spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, idx-mslit:idx+mslit))
            spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, idx-mslit:idx+mslit)) * redsnr   ! Reduce noise              
          ENDDO
        ENDIF
      ENDDO

      np_out = j
    ENDIF
    !print *, nx
    !IF (nx == 30) THEN
    !   WRITE(www_lun, '(2F10.4)') spec(1, np_out, 15)
    !   WRITE(www_lun, '(2D14.6)') spec(2, np_out, 15)
    !ELSE
    !   WRITE(www_lun, '(2F10.4)') spec(1, np_out, 29:30)
    !   WRITE(www_lun, '(2D14.6)') SUM(spec(2, np_out, 29:30)) / 2.
    !ENDIF
    !IF (np_out > np) THEN
    !   WRITE(www_lun, '(2A)') modulename, ': Improper sampling rate or slit width!!!'
    !   pge_error_status = pge_errstat_error
    !ENDIF

    RETURN
  END SUBROUTINE reduce_irrad_resolution

  SUBROUTINE reduce_rad_resolution (spec, qflag, np, nx, which_slit, slwth, &
       samprate, dwav, retswav, retewav, swav, ewav, np_out, pge_error_status)

    USE OMSAO_precision_module 
    USE OMSAO_parameters_module, ONLY: nmax=>max_spec_pts
    USE OMSAO_indices_module,    ONLY: sig_idx!, wvl_idx, spc_idx
    USE OMSAO_variables_module,  ONLY: use_redfixwav, nredfixwav, redfixwav
    USE ozprof_data_module,      ONLY: pos_alb
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: interpolation

    IMPLICIT NONE

    ! ====================
    ! In/Output variables
    ! ====================
    INTEGER, INTENT(IN)                              :: np, nx, which_slit
    INTEGER, INTENT(OUT)                             :: np_out, pge_error_status
    INTEGER (KIND=i2), DIMENSION(np, nx), INTENT(IN) :: qflag                   
    REAL (KIND=dp), INTENT(IN)                       :: samprate, slwth, dwav, swav, ewav, retswav, retewav
    REAL (KIND=dp), DIMENSION(sig_idx, np, nx), INTENT(INOUT) :: spec

    ! ====================
    ! Local variables
    ! ====================
    INTEGER                :: i, j, ix, mslit, nf, nsamp, nsamp1, nslit, errstat, fidx, lidx, iwin, idx
    INTEGER, DIMENSION(nx) :: nmod
    INTEGER, DIMENSION( 3) :: sidx, eidx, nstep
    REAL (KIND=dp)         :: dlam0, slitsum, redsnr, dx, fwav, lwav
    REAL (KIND=dp), DIMENSION(nmax)      :: slit
    REAL (KIND=dp), DIMENSION(sig_idx, np, nx) :: specmod
    REAL (KIND=dp), DIMENSION(sig_idx, nmax)   :: fnspec

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=28), PARAMETER :: modulename = 'reduce_irradiance_resolution'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    pge_error_status = pge_errstat_ok

    ! If slwth == 0, not allowed, unless that user provides a fixed wavelength grid
    IF (slwth == 0 .AND. .NOT. use_redfixwav) THEN
      WRITE(www_lun, '(2A)') modulename, ': Zero Slit Width Not allowed!!!'
      pge_error_status = pge_errstat_error
    ENDIF

    dlam0 = spec(1, 2, 1) - spec(1, 1, 1)  

    ! Filter bad measurements
    DO ix = 1, nx
      j = 0

      DO i = 1, np 
        IF (spec(2, i, ix) > 0. .AND. spec(2, i, ix) <= 4.0E14 .AND. qflag(i, ix) == 0) THEN
          j = j + 1
          specmod(:, j, ix) = spec(:, i, ix)
        ENDIF
      ENDDO
      nmod(ix) = j
    ENDDO

    IF (slwth == 0 .AND. use_redfixwav) THEN

      j = 0
      DO i = 1, nredfixwav
        IF (redfixwav(i) >= swav .AND. redfixwav(i) <=ewav) THEN
          j = j + 1
          spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
      ENDDO

    ELSE

      ! Establish fine wavelength scale (common to all across-track positions)
      nf = INT((ewav - swav) / dwav + 1)
      fnspec(1, 1) = swav
      DO i = 2, nf
        fnspec(1, i) = fnspec(1, i - 1) + dwav
      ENDDO

      IF (which_slit == 1) THEN               ! Symmetric Gaussian (slwth is hw1e)
        ! hw1e * 1.6551 = FWHM
        mslit = NINT(2.62826 * slwth / dwav) ! slit truncation (<0.1%)
        nsamp = INT(slwth * 1.66511 / samprate / dwav)
      ELSE IF (which_slit == 2) THEN          ! Triangular (slwth is FWHM)
        mslit = NINT(slwth / dwav)
        nsamp = INT(slwth / samprate / dwav)
      ENDIF
      nslit  = mslit * 2 + 1
      nsamp1 = INT(dlam0 / dwav)
      IF (nsamp1 < 1.0) nsamp1 = 1

      ! Set up slit function
      IF (which_slit == 1) THEN
        slit(1:nslit) = EXP( -((fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth)**2 )
      ELSE
        slit(1:nslit) = 1.0 - ABS(fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth
        WHERE (slit(1:nslit) < 0.0)
          slit(1:nslit) = 0.0
        ENDWHERE
      ENDIF
      slitsum = SUM(slit(1:nslit))
      slit(1:nslit) = slit(1:nslit) / slitsum ! Normalization
      redsnr  = 1.0 / SQRT(slitsum / dlam0 * dwav)                          ! Measurement error/noise reduction

      ! Convolve and Sample
      ! More sampling at the two ends of the selected spectral range
      ! Especially for the first and last 4 positions 
      ! This is to avoid extrapolation while keeping as many measurements as possible
      IF (retswav <= fnspec(1, mslit+1)) THEN
        i = mslit + 1 + 3 * nsamp1
      ELSE
        i = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= retswav + dlam0 * 3)))
      ENDIF
      sidx(1) = mslit+1
      eidx(1) = i
      nstep(1) = nsamp1

      IF (retewav >= fnspec(1, nf - mslit)) THEN
        i = nf - mslit - nsamp1 * 3
      ELSE
        i = MAXVAL(MAXLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) <= retewav - dlam0 * 3)))
      ENDIF
      sidx(2) = eidx(1) + nsamp1
      eidx(2) = i-1
      nstep(2) = nsamp
      sidx(3) = i + nsamp1
      eidx(3) = nf - mslit
      nstep(3) = nsamp1

      j = 0
      DO i = 1, nredfixwav
        IF (redfixwav(i) > fnspec(1, mslit + 1) .AND. redfixwav(i) < fnspec(1, nf - mslit)) THEN
          j = j + 1
          spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
      ENDDO

    ENDIF

    IF (use_redfixwav) THEN
      ! Add 3 * 2 extra wavelengths for irradiance and 2 * 2 extra wavelengths for radiance
      ! Add 3 wavelengths at the beginning of a spectra region
      DO i = 1, j
        fwav = MAX(swav, retswav)
        IF (spec(1, i, 1) > fwav) THEN
          IF (i == 1) THEN
            dx = (spec(1, i, 1) - fwav) / 3.
          ELSE
            dx = (spec(1, i, 1) - MAX(spec(1, i-1, 1), fwav)) / 3.
          ENDIF
          IF (dx > dlam0) dx = dlam0
          spec(1, i+3:j+3, 1:nx) = spec(1, i:j, 1:nx)
          spec(1, i, 1:nx)   = spec(1, i+3, 1:nx) - dx * 3.0
          spec(1, i+1, 1:nx)   = spec(1, i+3, 1:nx) - dx * 2.0
          spec(1, i+2, 1:nx)   = spec(1, i+3, 1:nx) - dx * 1.0
          j = j + 3
          EXIT
        ENDIF
      ENDDO

      ! Add 3 wavelengths at the end of a spectra region
      DO i = j, 1, -1
        lwav = MIN(ewav, retewav)
        IF (spec(1, i, 1) < lwav) THEN
          IF (i == j) THEN
            dx = (lwav - spec(1, i, 1)) / 3.
          ELSE
            dx = ( MAX(spec(1, i+1, 1), lwav) - spec(1, i, 1)) / 3.
          ENDIF
          IF (dx > dlam0) dx = dlam0
          spec(1, i+4:j+3, 1:nx) = spec(1, i+1:j, 1:nx)
          spec(1, i+1, 1:nx)   = spec(1, i, 1:nx) + dx * 1.0
          spec(1, i+2, 1:nx)   = spec(1, i, 1:nx) + dx * 2.0
          spec(1, i+3, 1:nx)   = spec(1, i, 1:nx) + dx * 3.0
          j = j + 3
          EXIT
        ENDIF
      ENDDO

      ! Add a wavelength at the wavelength to derive initial cloud fraction
      IF (swav < pos_alb .AND. ewav > pos_alb) THEN
        DO i = j, 1, -1 
          IF (spec(1, i, 1) < pos_alb) THEN
            spec(1, i+1:j, 1:nx) = spec(1, i+2:j+1, 1:nx)
            spec(1, i+1, 1:nx) = pos_alb
            j = j + 1
            EXIT
          ENDIF
        ENDDO
      ENDIF

    ENDIF
    np_out = j

    IF (slwth == 0 .AND. use_redfixwav) THEN

      DO ix = 1, nx
        DO i = 2, 3
          CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
               np_out, spec(1, 1:np_out, ix),  spec(i, 1:np_out, ix), pge_error_status )
          IF ( pge_error_status > pge_errstat_warning ) THEN
            errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
            RETURN
          END IF
        ENDDO
      ENDDO

    ELSE
      DO ix = 1, nx 
        ! Pre-interpolation
        fidx = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= specmod(1, 1, ix))))
        lidx = MAXVAL(MAXLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) <= specmod(1, nmod(ix), ix))))
        fnspec(2:3, 1:fidx-1) = 0.0
        fnspec(2:3, lidx+1:nf) = 0.0

        IF (nmod(ix) > np * 0.75) THEN
          DO i = 2, 3           
            CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
                 lidx-fidx+1, fnspec(1, fidx:lidx), fnspec(i, fidx:lidx), pge_error_status )
            IF ( pge_error_status > pge_errstat_warning ) THEN
              errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 
              RETURN
            END IF
          ENDDO
        ELSE
          fnspec(2:3, 1:nf) = 0.0
        ENDIF

        IF (.NOT. use_redfixwav) THEN
          ! Convolve and Sample
          j = 0
          DO iwin = 1, 3
            DO i = sidx(iwin), eidx(iwin), nstep(iwin)
              j = j + 1

              spec(1, j, ix) = fnspec(1, i)
              spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, i-mslit:i+mslit))
              spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, i-mslit:i+mslit)) * redsnr ! Reduce noise

              IF ((i - mslit < fidx) .OR. (i + mslit > lidx)) THEN
                spec(2:3, j, ix) = 0.0
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO j = 1,  np_out
            idx = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= spec(1, j, ix))))
            IF ( ABS(fnspec(1, idx-1) - spec(1, j, ix)) < ABS(fnspec(1, idx) - spec(1, j, ix))) idx = idx - 1
            spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, idx-mslit:idx+mslit))
            spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, idx-mslit:idx+mslit)) * redsnr   ! Reduce noise              
          ENDDO
        ENDIF
      ENDDO

      np_out = j
    ENDIF


    RETURN
  END SUBROUTINE reduce_rad_resolution

  !  Unused
  !
  !  SUBROUTINE SPIKE_DETECT_CORRECT1(ns, fitspec, simrad)
  !
  !    USE OMSAO_precision_module
  !    USE OMSAO_variables_module, ONLY : fitweights, currspec!, &
  !         !poly_order, fitwavs
  !    USE ozprof_data_module,     ONLY : use_lograd
  !    USE OMSAO_errstat_module,   ONLY : www_lun
  !    IMPLICIT NONE
  !
  !    ! =======================
  !    ! Input/Output variables
  !    ! =======================
  !    INTEGER, INTENT(IN)                             :: ns
  !    REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: simrad
  !    REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: fitspec
  !
  !    ! =======================
  !    ! local variables
  !    ! =======================
  !    INTEGER                        :: i, nspike, ncorr, approach, iter!, j
  !    REAL (KIND=dp), DIMENSION(ns)  :: dfthresh, diff, mratio, sratio!, &
  !         !reldf, reldv, relavg
  !    REAL (KIND=dp)                 :: rms, thresh
  !
  !    IF (use_lograd) THEN
  !      simrad = EXP(simrad)
  !      fitspec = EXP(fitspec)
  !    ENDIF
  !    !WRITE(90, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)
  !
  !    sratio(1:ns-1) = simrad(1:ns-1)  /  simrad(2:ns)
  !    sratio(ns) = sratio(ns-1)
  !    dfthresh(1:ns-1) = &
  !         2.0 * SQRT(fitweights(1:ns-1)**2.0 + fitweights(2:ns)**2.0)
  !    rms = SUM(dfthresh(1:ns-1)) / (ns-1.0D0)
  !
  !    nspike = 0
  !    ncorr = 0
  !    approach = 1
  !
  !    DO iter = 1, 1
  !      mratio(1:ns-1) = fitspec(1:ns-1) / fitspec(2:ns)
  !      mratio(ns) = mratio(ns-1)
  !      diff(1:ns-1)   = mratio(1:ns-1) - sratio(1:ns-1)
  !
  !      DO i = ns - 2, 1, - 1
  !        thresh = MAX(dfthresh(i), rms)
  !        IF (ABS(diff(i)) > thresh ) THEN
  !          nspike = nspike + 1
  !        ENDIF
  !
  !        IF (diff(i) < -thresh ) THEN
  !          ncorr = ncorr + 1
  !          fitspec (i+1)  = fitspec(i+1)  * mratio(i) / sratio(i)
  !          currspec(i+1)  = currspec(i+1) * mratio(i) / sratio(i)
  !          mratio(i) = fitspec(i)/fitspec(i+1)
  !          mratio(i+1) = fitspec(i+1)/fitspec(i+2)
  !          diff(i+1) = mratio(i+1) - sratio(i+1)
  !        ELSE IF  (diff(i) > thresh) THEN
  !          ncorr = ncorr + 1
  !          fitspec(i)  = fitspec(i)  * sratio(i) / mratio(i)
  !          currspec(i) = currspec(i) * sratio(i) / mratio(i)
  !          mratio(i-1) = fitspec(i-1) / fitspec(i)
  !          diff(i-1)   = mratio(i-1) - sratio(i-1)
  !        ENDIF
  !      ENDDO
  !
  !      WRITE(www_lun, *) 'Number of spikes = ', nspike, thresh
  !      WRITE(www_lun, *) 'Number of corrections = ', ncorr
  !    ENDDO
  !
  !    !WRITE(91, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)
  !    IF (use_lograd) THEN
  !      simrad = LOG(simrad)
  !      fitspec = LOG(fitspec)
  !    ENDIF
  !
  !    RETURN
  !
  !  END SUBROUTINE SPIKE_DETECT_CORRECT1




  SUBROUTINE ROUGH_SPIKE_DETECT(ns, waves, rad, sol, nspike)

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN)                             :: ns
    INTEGER, INTENT(OUT)                            :: nspike
    REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: rad
    REAL (KIND=dp), INTENT(IN), DIMENSION (ns)      :: sol, waves

    ! =======================
    ! local variables
    ! =======================
    INTEGER                       :: i!, j
    REAL (KIND=dp), PARAMETER     :: thresh = 0.20
    REAL (KIND=dp)                :: diff, oldratio, newratio
    REAL (KIND=dp), DIMENSION(ns) :: normrad

    !WRITE(90, '(f8.3, 2d14.6)') ((waves(i), rad(i), sol(i)), i=1, ns)

    normrad = rad / sol
    nspike  = 0
    DO i = ns - 3, 1, - 1      ! only detect spike below 305 nm

      IF (waves(i) > 305.0D0) CYCLE

      oldratio = normrad(i) / normrad(i+1)
      diff = oldratio-1.0
      IF (diff > thresh ) THEN   !  This pixel got large error
        nspike = nspike + 1
        newratio = normrad(i+1) / normrad(i+2)

        IF (ABS(newratio - 1.0) < 0.5 * thresh) THEN
          normrad(i) = normrad(i+1) * newratio
          rad(i) = rad(i) * newratio / oldratio
        ELSE
          rad(i+1) =  rad(i+1) * normrad(i+2) / normrad(i+3) / newratio
          newratio =  normrad(i+2) / normrad(i+3)
          rad(i) = rad(i) * newratio / oldratio
        ENDIF
      ENDIF
    ENDDO
    !WRITE(www_lun, *) 'Number of spikes = ', nspike

    RETURN
  END SUBROUTINE ROUGH_SPIKE_DETECT


    SUBROUTINE SPIKE_DETECT_CORRECT(ns, fitspec, simrad)
  
     USE OMSAO_precision_module
      USE OMSAO_variables_module, ONLY : currspec!, fitwavs, fitweights
      USE ozprof_data_module,     ONLY : use_lograd
      IMPLICIT NONE
  
      ! =======================
      ! Input/Output variables
      ! =======================
      INTEGER, INTENT(IN)                             :: ns
      REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: simrad
      REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: fitspec
  
      ! =======================
      ! local variables
      ! =======================
      INTEGER                        :: i, approach
      REAL (KIND=dp), DIMENSION(ns)  :: mratio, sratio
      REAL (KIND=dp)                 :: fitavg, simavg, fitavg3, &
           simavg3, fitavg2, simavg2
  
      IF (use_lograd) THEN
        simrad = EXP(simrad)
        fitspec = EXP(fitspec)
      ENDIF
      !WRITE(90, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)
  
     ! use approach 1 is better
      approach = 1
  
    IF (approach == 1) THEN
  
       ! Assume the average of last 10 pixels are without errors
        ! Also assume the ratio of adjacent pixels is less dependent on
        ! ozone profile
        sratio(1:ns-1) = simrad(1:ns-1)  /  simrad(2:ns)
        sratio(ns) = sratio(ns-1)
        mratio(1:ns-1) = fitspec(1:ns-1) / fitspec(2:ns)
        mratio(ns) = mratio(ns-1)
  
        fitavg = SUM(fitspec(ns-9:ns)) / 10.0D0
        simavg = SUM(simrad(ns-9:ns))  / 10.0D0
  
       sratio(ns-10)   = simrad(ns-10)  / simavg
        mratio(ns-10)   = fitspec(ns-10) / fitavg
        fitspec(ns-10)  = fitspec(ns-10)  * sratio(ns-10) / mratio(ns-10)
        currspec(ns-10) = currspec(ns-10) * sratio(ns-10) / mratio(ns-10)
  
        DO i = ns - 11, 1, - 1
          mratio(i)  = fitspec(i)  / fitspec(i + 1)
          fitspec(i) = fitspec(i)  * sratio(i) / mratio(i)
          currspec(i)= currspec(i) * sratio(i) / mratio(i)
        ENDDO
      ELSE
  
       ! Assume the average of last 20 pixels are without errors
        ! Also assume the ratio of adjacent pixels (first * third / second^2) &
        ! is less dependent on ozone profile
  
       sratio(1:ns-2) = simrad(1:ns-2)  * simrad(3:ns)  / (simrad(2:ns-1)**2.0)
        mratio(1:ns-2) = fitspec(1:ns-2) * fitspec(3:ns) / (fitspec(2:ns-1)**2.0)
  
        fitavg3= SUM(fitspec(ns-9:ns)) / 10.0D0
        simavg3= SUM(simrad(ns-9:ns))  / 10.0D0
        fitavg2= SUM(fitspec(ns-19:ns-10)) / 10.0D0
        simavg2= SUM(simrad(ns-19:ns-10))  / 10.0D0
  
        sratio(ns-20)  = simrad(ns-20)  * simavg3 / (simavg2**2.0)
        mratio(ns-20)  = fitspec(ns-20) * fitavg3 / (fitavg2**2.0)
        fitspec (ns-20)= fitspec(ns-20) * sratio(ns-20) / mratio(ns-20)
        currspec(ns-20)= currspec(ns-20)* sratio(ns-20) / mratio(ns-20)
  
        sratio(ns-21)  = simrad(ns-21)  * simavg2 / (simrad(ns-20)**2.0)
        mratio(ns-21)  = fitspec(ns-21) * fitavg2 / (fitspec(ns-20)**2.0)
        fitspec (ns-21)= fitspec(ns-21) * sratio(ns-21) / mratio(ns-21)
       currspec(ns-21)= currspec(ns-21)* sratio(ns-21) / mratio(ns-21)
  
        DO i = ns - 22, 1, - 1
          mratio(i)  = fitspec(i)  * fitspec(i+2) / (fitspec(i + 1)**2.0)
          fitspec(i) = fitspec(i)  * sratio(i)    /  mratio(i)
         currspec(i)= currspec(i) * sratio(i)    /  mratio(i)
       ENDDO
     ENDIF
  
      !WRITE(91, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)
      IF (use_lograd) THEN
        simrad = LOG(simrad)
        fitspec = LOG(fitspec)
      ENDIF
 
     RETURN
  
    END SUBROUTINE SPIKE_DETECT_CORRECT


  !  Unused
  !
  !  SUBROUTINE ROUGH_SPIKE_DETECT1(ns, waves, rad, sol)
  !
  !    USE OMSAO_precision_module
  !    USE OMSAO_errstat_module
  !    IMPLICIT NONE
  !
  !    ! =======================
  !    ! Input/Output variables
  !    ! =======================
  !    INTEGER, INTENT(IN)                             :: ns
  !    REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: rad
  !    REAL (KIND=dp), INTENT(IN), DIMENSION (ns)      :: sol, waves
  !
  !    ! =======================
  !    ! local variables
  !    ! =======================
  !    INTEGER                       :: i, nspike!, j
  !    REAL (KIND=dp), PARAMETER     :: thresh = 0.15
  !    REAL (KIND=dp)                :: diff, oldratio, newratio, newratio2
  !    REAL (KIND=dp), DIMENSION(ns) :: normrad
  !
  !    !WRITE(90, '(f8.3, 2d14.6)') ((waves(i), rad(i), sol(i)), i=1, ns)
  !
  !    normrad = rad / sol
  !    nspike  = 0
  !    DO i = ns - 5, 1, - 1      ! only detect spike below 305 nm
  !
  !      IF (waves(i) > 305.0D0) CYCLE
  !
  !      oldratio = normrad(i) * normrad(i+2) / (normrad(i+1)**2.0)
  !      diff = oldratio-1.0
  !
  !      IF (diff > thresh ) THEN   !  This pixel got large error
  !        nspike = nspike + 1
  !        newratio = normrad(i+1) * normrad(i+3) / (normrad(i+2)**2.0)
  !
  !        IF (ABS(newratio - 1.0) < 0.5 * thresh) THEN
  !          !write(*, *) i, oldratio, normrad(i)
  !          normrad(i) = normrad(i+1)**2.0 * newratio / normrad(i+2)
  !          rad(i) = rad(i) * newratio / oldratio
  !          !write(*, *) i, newratio, normrad(i)
  !        ELSE
  !          !write(*, *) i, oldratio, normrad(i), normrad(i+1)
  !          newratio2 = normrad(i+2) * normrad(i+4) / (normrad(i+3)**2.0)
  !          rad(i+1) =  rad(i+1) * newratio2 / newratio
  !          rad(i) = rad(i) * newratio2 / oldratio
  !          !write(*, *) i, newratio, normrad(i), normrad(i+1)
  !        ENDIF
  !      ENDIF
  !    ENDDO
  !    WRITE(www_lun, *) 'Number of spikes = ', nspike
  !
  !    RETURN
  !  END SUBROUTINE ROUGH_SPIKE_DETECT1



  ! xliu, 05/23/2010
  ! Identify spikes (e.g. due to emission lines, protons) in UV-1
  ! And modify measurement error
  SUBROUTINE UV1_SPIKE_DETECT(ns, fitspec, simrad, nspike)

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY : fitweights, nradpix!, fitwavs, currspec
    USE ozprof_data_module,     ONLY : use_lograd

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN)                             :: ns
    REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: simrad
    REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: fitspec
    INTEGER, INTENT(OUT)                            :: nspike

    ! =======================
    ! local variables
    ! =======================
    INTEGER                        :: iw, j
    REAL (KIND=dp)                 :: rms, thresh, diff, fratio, sratio

    IF (use_lograd) THEN
      simrad = EXP(simrad)
      fitspec = EXP(fitspec)
    ENDIF

    nspike = 0
    j = nradpix(1)
    DO iw =  nradpix(1) - 1, 1, -1
      fratio = fitspec(iw) / fitspec(j)
      sratio = simrad(iw) / simrad(j)

      rms = SQRT(fitweights(iw)**2 + fitweights(j)**2) * 3.0 * &
           SQRT(( 1.0 * j - iw))
      thresh = MAX(rms, 0.06)
      diff   = ABS(fratio - sratio)
      IF ( diff > thresh) THEN
        nspike = nspike + 1
        fitweights(iw) = diff / 3.0
        !WRITE(*, '(I5, 4F10.4)') nspike, fitwavs(iw), rms, thresh, diff
      ELSE
        j = iw
      ENDIF

    ENDDO
    !WRITE(*, *) 'Number of spikes = ', nspike

    IF (use_lograd) THEN
      simrad = LOG(simrad)
      fitspec = LOG(fitspec)
    ENDIF

    !DO iw = 1, ns
    !   WRITE(90, *) fitwavs(iw), simrad(iw), fitspec(iw)
    !ENDDO
    !STOP 1

    RETURN
  END SUBROUTINE UV1_SPIKE_DETECT

 SUBROUTINE cubic_specfit ( a, na, y, m, ctrl, dyda, mdy )

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY  : cubic_x, cubic_y, cubic_w
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! Input parameters
    ! ================
    INTEGER,                         INTENT (IN)  :: na, m, mdy
    REAL (KIND=dp), DIMENSION (na),  INTENT (IN)  :: a

    ! Modified parameters
    ! ===================
    INTEGER, INTENT (INOUT) :: ctrl

    ! Output parameters
    ! =================
    REAL (KIND=dp), DIMENSION (m),    INTENT (OUT)  :: y
    REAL (KIND=dp), DIMENSION (m,na), INTENT (OUT)  :: dyda

    ! Local variables
    ! ===============
    REAL (KIND=dp), DIMENSION (m) :: x, y0
    x  = cubic_x(1:m)
    y0 = a(1) + a(2)*x + a(3)*x*x + a(4)*x*x*x


    SELECT CASE ( ABS(ctrl) )
    CASE ( 1 )
      y  = ( y0 - cubic_y(1:m) ) / cubic_w(1:m)
    CASE ( 2 )
      dyda = 0.0_dp
      dyda(1:m,1) = 1.0_dp
      dyda(1:m,2) = x(1:m)
      dyda(1:m,3) = x(1:m)*x(1:m)
      dyda(1:m,4) = x(1:m)*x(1:m)*x(1:m)
    CASE ( 3 )
      ! This CASE is included to get the complete fitted spectrum
      y  = y0
    CASE DEFAULT
      WRITE (www_lun, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE cubic_specfit
 SUBROUTINE poly_specfit ( a, na, y, m, ctrl, dyda, mdy )

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY  : poly_x, poly_y, poly_w
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! Input parameters
    ! ================
    INTEGER,                         INTENT (IN)  :: na, m, mdy
    REAL (KIND=dp), DIMENSION (na),  INTENT (IN)  :: a

    ! Modified parameters
    ! ===================
    INTEGER, INTENT (INOUT) :: ctrl

    ! Output parameters
    ! =================
    REAL (KIND=dp), DIMENSION (m),    INTENT (OUT)  :: y
    REAL (KIND=dp), DIMENSION (m,na), INTENT (OUT)  :: dyda

    ! Local variables
    ! ===============
    REAL (KIND=dp), DIMENSION (m) :: x, y0
    INTEGER                       :: i
    x  = poly_x(1:m)

    y0 = 0.0_dp
    DO i = 1, na
      y0 = y0 + a(i) * (x ** (i-1))
    END DO

    SELECT CASE ( ABS(ctrl) )
    CASE ( 1 )
      y  = ( y0 - poly_y(1:m) ) / poly_w(1:m)
    CASE ( 2 )
      dyda = 0.0_dp
      DO i = 1, na
        dyda(1:m, i) = x(1:m) ** (i-1)
      END DO
    CASE ( 3 )
      ! This CASE is included to get the complete fitted spectrum
      y  = y0
    CASE DEFAULT
      WRITE (www_lun, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE poly_specfit

  SUBROUTINE poly_fit (locwvl, npoints, locspec, ll_rad, lu_rad, r )

  USE OMSAO_precision_module
  !USE OMSAO_parameters_module, ONLY : elsunc_np, elsunc_nw
  USE OMSAO_variables_module,  ONLY : poly_x, poly_y, poly_w, poly_order
  USE m_bounded_nonlin_LS,       ONLY : elsunc

  IMPLICIT NONE

  INTEGER,                             INTENT (IN) :: npoints, ll_rad, lu_rad
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: locwvl

  REAL (KIND=dp), DIMENSION (npoints),    INTENT (INOUT) :: locspec
  REAL (KIND=dp), DIMENSION (poly_order), INTENT (OUT)   :: r

  REAL (KIND=dp), DIMENSION (npoints)              :: tmp, ptmp, sig
  !REAL (KIND=dp), DIMENSION (poly_order)           :: par, polfunc
  !REAL (KIND=dp), DIMENSION (poly_order,poly_order):: covar

  INTEGER                             :: i, nlower, nupper, nfitted
  REAL (KIND=dp)                      :: locavg
  REAL (KIND=dp), DIMENSION (npoints) :: x

  ! ================
  ! ELSUNC variables
  ! ================
  INTEGER                                       :: exval
  !REAL (KIND=dp), DIMENSION (npoints)           :: fitres
  REAL (KIND=dp)                                :: chisq2 !, rms
  INTEGER                                       :: elbnd
  INTEGER,        DIMENSION (11)                :: p
  REAL (KIND=dp), DIMENSION (6)                 :: w
  REAL (KIND=dp), DIMENSION (poly_order)         :: blow, bupp
  REAL (KIND=dp), DIMENSION (npoints)           :: f !, fitspec
  REAL (KIND=dp), DIMENSION (npoints,poly_order) :: dfda


  ! ======================
  ! Assign fitting weights
  ! ======================
  sig = 1.0_dp

  !     Find limits for polynomial fitting, with ~1 nm overlap
  nlower = MINVAL(MINLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) >= locwvl(ll_rad)-1.0) ))
  nupper = MAXVAL(MAXLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) <= locwvl(lu_rad)+1.0) ))
  nfitted = nupper - nlower + 1

  !     Find average position over fitted region
  locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=dp )

  !     Load temporary position file: re-define positions in order to fit
  !     about mean position
  DO i = 1, nfitted
     ptmp(i) = locwvl(i+nlower-1) - locavg
  END DO

  !     Load and fit spectrum
  tmp(1:nfitted) = locspec(1+nlower-1:nfitted+nlower-1)

  ! ===============================================================
  ! ELBND: 0 = unconstrained
  !        1 = all variables have same lower bound
  !        else: lower and upper bounds must be supplied by the use
  ! ===============================================================  
  elbnd = 0  ;  exval = 0
  p   = -1 ; p(1) = 0  ;  p(3) = 5  ; w = -1.0
  blow(1:poly_order) = 0.0  ;  bupp(1:poly_order) = 0.0
  
  poly_x(1:nfitted) = ptmp(1:nfitted)
  poly_y(1:nfitted) =  tmp(1:nfitted)
  poly_w(1:nfitted) =  sig(1:nfitted)

  r = 0.0 ; f = 0.0 ; dfda = 0.0

  CALL elsunc ( &
       r, poly_order, nfitted, nfitted, poly_specfit, elbnd, blow(1:poly_order), &
       bupp(1:poly_order), p, w, exval, f(1:nfitted), dfda(1:nfitted,1:poly_order) )
  chisq2 = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines

  ! Re-load spec with high-pass filtered data, over whole  spectral region
  x(1:npoints) = locwvl(1:npoints) - locavg

  poly_x(1:npoints) = x(1:npoints)
  exval = 3
  CALL poly_specfit ( &
       r(1:poly_order), poly_order, poly_y(1:npoints), npoints, exval, dfda, 0 )
  locspec(1:npoints) = poly_y(1:npoints)

  RETURN
  END SUBROUTINE poly_fit

end module m_fitting_util
