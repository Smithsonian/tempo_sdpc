!
module m_cross_calibrate
  USE m_solar_fit, ONLY:solar_fit, solar_fit_vary
  USE m_solar_wavcal_vary
  USE m_radiance_fit, ONLY: radiance_fit_vary
  USE m_radiance_wavcal, ONLY: radiance_wavcal, radiance_wavcal_vary
  USE m_ezspline_interpolation, only: interpol
  USE OMSAO_variables_module, ONLY:calwrt, calunit,the_pix, the_line, currpix, &
            ring_group, irrad_group, cali_group,cali_group, rad_group
  USE OMSAO_indices_module
  INTEGER, PARAMETER, PRIVATE :: n_fitvar = 6
  INTEGER, PARAMETER, DIMENSION(1:n_fitvar), PRIVATE :: & 
     fit_idx=(/shi_idx, squ_idx, sin_idx, hwe_idx, asy_idx, spk_idx/)
  public calibrate_irrad_cross, calibrate_rad_cross, & 
         solwavcal_coadd, radwavcal_coadd

  private

contains

  SUBROUTINE calibrate_irrad_cross (irrad, ring, cali, first_pix, last_pix, pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: normweight
    USE OMSAO_variables_module, ONLY: numwin, nxbin, scnwrt,band_selectors, & 
         use_meas_sig, yn_varyslit, which_slit, wavcal, wavcal_sol, &
         reduce_resolution, redslw, &
         curr_sol_spec, n_irrad_wvl, nsolpix, &
         wincal_wav, solwinfit, sol_spec_ring, nsol_ring, slitwav_sol, &
         sswav_sol, nslit_sol, nwavcal_sol, solslitfit,slit_fname  
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ================
    ! Output variables
    ! ================
    TYPE(irrad_group), INTENT(INOUT) :: irrad
    TYPE(ring_group), INTENT(INOUT) :: ring
    TYPE(cali_group), INTENT(INOUT) :: cali
    INTEGER, INTENT (IN)          :: first_pix, last_pix
    INTEGER, INTENT (OUT)         :: pge_error_status

    ! ===============
    ! Local variables
    ! ===============
    INTEGER                     :: i, ix, fidx, lidx
    REAL (KIND=dp)              :: dwvl
    LOGICAL                     :: error
    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=25), PARAMETER :: modulename = 'm_irrad_cross_calibrate'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg
    INTEGER :: errstat

    IF (calwrt) THEN
       OPEN(calunit, FILE=TRIM(ADJUSTL(slit_fname)), STATUS='unknown') 
       WRITE(calunit, '(5i5)') which_slit, n_fitvar + 1
       WRITE(calunit, '(5i5)') last_pix-first_pix + 1, first_pix, last_pix
       WRITE(calunit,'(10(a5))') 'wav', calfit_strings(fit_idx)
    ENDIF
    pge_error_status = pge_errstat_ok
    IF (scnwrt) WRITE(*, *) 'yn_varyslit=',yn_varyslit,'which_slit=',which_slit
    
    DO ix = first_pix, last_pix
      the_pix = (ix-1)*nxbin +1
      currpix = ix
      ! Bad solar, no calibration
      IF (irrad%errstat(ix) == pge_errstat_error) CYCLE
      error   = .FALSE.

      ! Spectrum
      n_irrad_wvl = irrad%nwav(ix)
      curr_sol_spec(wvl_idx,1:n_irrad_wvl) = irrad%wavl(1:n_irrad_wvl, ix)
      curr_sol_spec(spc_idx,1:n_irrad_wvl) = irrad%spec(1:n_irrad_wvl, ix)

      IF (use_meas_sig) THEN
        curr_sol_spec(sig_idx, 1:n_irrad_wvl) = irrad%prec(1:n_irrad_wvl, ix)
      ELSE
        curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
      ENDIF
      nsolpix(1:numwin) = irrad%npix(1:numwin, ix)


      ! Ring Spectrum
      nsol_ring = ring%winpix(ix,2) - ring%winpix(ix,1) + 1
      sol_spec_ring(1, 1:nsol_ring) = ring%wavl(ring%winpix(ix,1):ring%winpix(ix,2), ix)

      ! Perform calibration
      IF (scnwrt) WRITE(*, '(A,I5)') &
           'Performing solar wavelength calibration: ', (ix - 1) * nxbin + 1

      IF (yn_varyslit) THEN
        IF (which_slit < 5 ) THEN
          CALL solar_fit_vary (error )
          cali%nslit_sol(ix) = nslit_sol
          cali%slitfit_sol(1:nslit_sol, 1:max_calfit_idx, :, ix) = &
                          solslitfit(1:nslit_sol, 1:max_calfit_idx, :)
          cali%slitwav_sol(1:nslit_sol, ix) = slitwav_sol(1:nslit_sol)
        ELSE IF (which_slit == 5) THEN 
          wavcal_sol = .true.
        ENDIF

        IF (wavcal .AND. wavcal_sol) THEN
          CALL solar_wavcal_vary(calunit, error)
          cali%nwavcal_sol(ix) = nwavcal_sol
          cali%sswav_sol(1:nwavcal_sol, ix) = sswav_sol(1:nwavcal_sol)
        ENDIF
      ELSE
        IF (wavcal .OR. which_slit /= 5) THEN
          CALL solar_fit (error )
          cali%wincal_wav(1:numwin, ix) = wincal_wav(1:numwin)
          IF (reduce_resolution) THEN
            solwinfit(1:numwin, hwe_idx, 1) = &
                 SQRT(solwinfit(1:numwin, hwe_idx, 1) ** 2. &
                 + redslw(band_selectors(1:numwin))**2)
          ENDIF
            cali%solwinfit(1:numwin, 1:max_calfit_idx, ix) = &
               solwinfit(1:numwin, 1:max_calfit_idx, 1)
        ENDIF
      END IF
      
      IF(calwrt) CALL write_cali_cross (calunit,ix,cali,.false.)
      ! Check wavelength shifts, if it is too large, do not apply
      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nsolpix(i) - 1
        dwvl =  curr_sol_spec(wvl_idx, fidx+1) - curr_sol_spec(wvl_idx, fidx)

        IF ( ANY(ABS(curr_sol_spec(wvl_idx, fidx:lidx) - &
             irrad%wavl(fidx:lidx, ix)) > 1.0) .OR. & !dwvl*2.0) .OR. &
             ANY(curr_sol_spec(wvl_idx, fidx+1:lidx) - &
             curr_sol_spec(wvl_idx, fidx:lidx-1) < 0.0) ) THEN
             WRITE(*,*) 'solshi is ', maxval(abs(curr_sol_spec(wvl_idx,fidx:lidx)-irrad%wavl(fidx:lidx, ix)  ))
             error = .TRUE.
        ELSE
             irrad%wavl(fidx:lidx, ix) = &
               real(curr_sol_spec(wvl_idx, fidx:lidx) , kind=r4)
        ENDIF
        fidx = lidx + 1
      ENDDO

      IF (error) THEN
        irrad%errstat(ix) = pge_errstat_warning; CYCLE
      ENDIF

      ring%wavl(ring%winpix(ix,1):ring%winpix(ix,2), ix) = &
           real(sol_spec_ring(1, 1:nsol_ring) , kind=r4)
    ENDDO

    errstat = OMI_SMF_setmsg(OMI_S_SUCCESS, "done irrad cross calibrate.", &
         modulename, 0)

    RETURN
  END SUBROUTINE calibrate_irrad_cross


  SUBROUTINE calibrate_rad_cross (rad,cali, first_pix, last_pix, first_line, last_line,pge_error_status)
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: normweight, max_fit_pts!, downweight
    USE OMSAO_variables_module,  ONLY: curr_rad_spec, n_rad_wvl, &
         use_meas_sig, yn_varyslit, wavcal_sol, slit_rad, nradpix, &
         wincal_wav, solwinfit, numwin, slitwav_rad, nslit_rad, radslitfit, &
         radwinfit, nslit, slitwav, slitfit, nwavcal_rad,  sswav_rad, &
         which_slit, scnwrt, wavcal, nxbin, nybin,ntimes_loop,offset_line,  rslit_fname       
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ================
    ! Output variables
    ! ================
    INTEGER, INTENT (IN)        :: first_pix, last_pix, first_line, last_line
    INTEGER, INTENT (OUT)       :: pge_error_status
    TYPE (rad_group), INTENT(INOUT) :: rad
    TYPE (cali_group), INTENT(INOUT) :: cali

    ! ===============
    ! Local variables
    ! ===============
    REAL (KIND=dp), DIMENSION(max_fit_pts) :: diff, corr, wav
    REAL (KIND=dp)                         :: dwvl
    INTEGER  :: i, ix, iline, mline, inter_errstat, np, fidx, lidx!, j
    LOGICAL  :: error
    LOGICAL, SAVE :: first=.true.

    IF (calwrt) THEN
      IF (first) THEN 
       OPEN(calunit, FILE=TRIM(ADJUSTL(rslit_fname)), STATUS='unknown') 
       WRITE(calunit, '(5i5)') which_slit, n_fitvar + 1
       WRITE(calunit, '(5i5)') last_pix-first_pix + 1, first_pix, last_pix
       WRITE(calunit,'(10(a5))') 'wav', calfit_strings(fit_idx)
       first = .false.
      ENDIF
    ENDIF

    pge_error_status = pge_errstat_ok; error   = .FALSE.
    mline = NINT(ntimes_loop / 2.0) - 1
    
    DO ix = first_pix, last_pix
      currpix = ix
      the_pix = (ix-1)*nxbin + 1
      np = MAXVAL(rad%nwav(ix, 0:ntimes_loop-1))

      IF (ALL( rad%pix_errstat(ix, 0:ntimes_loop-1) == pge_errstat_error)) &
           CYCLE

      ! Only process a scan line in the middle
      DO i = mline, 0, -1  ! first half
        IF (rad%pix_errstat(ix, i) == pge_errstat_ok .AND. &
             rad%nwav(ix, i) == np) EXIT
      ENDDO

      IF ( i < 0 ) THEN    ! second half
        DO i = mline + 1, ntimes_loop - 1
          IF (rad%pix_errstat(ix, i) == pge_errstat_ok .AND. &
               rad%nwav(ix, i) == np) EXIT
        ENDDO
      ENDIF

      IF (i > ntimes_loop - 1) CYCLE ! no valid pixels for this xtrack position
      iline = i

      ! Spectrum

      n_rad_wvl = rad%nwav(ix, iline)
      curr_rad_spec(wvl_idx, 1:n_rad_wvl) = rad%wavl(1:n_rad_wvl, ix, iline)
      curr_rad_spec(spc_idx, 1:n_rad_wvl) = rad%spec(1:n_rad_wvl, ix, iline)

      IF (use_meas_sig) THEN
        curr_rad_spec(sig_idx, 1:n_rad_wvl) =rad%prec(1:n_rad_wvl, ix, iline)
      ELSE
        curr_rad_spec(sig_idx, 1:n_rad_wvl) = normweight
      ENDIF
      nradpix(1:numwin) = rad%npix(1:numwin, ix, iline)
      
      ! Perform calibration
      the_line = iline*nybin + offset_line + 1
      IF (scnwrt) WRITE(*, '(A, 3I5)') &
           '-Performing radiance wavelength calibration: ', &
           (ix - 1) * nxbin + 1, iline * nybin + offset_line + 1, ntimes_loop
      IF (which_slit < 5) THEN
        IF (yn_varyslit .AND. slit_rad ) THEN
          CALL radiance_fit_vary(n_rad_wvl, &
               curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), error)
          cali%nslit_rad(ix) = nslit_rad
          cali%slitfit_rad(1:nslit_rad, 1:max_calfit_idx, :, ix) = &
               radslitfit(1:nslit_rad, 1:max_calfit_idx, :)
          cali%slitwav_rad(1:nslit_rad, ix) = slitwav_rad(1:nslit_rad)
        ELSE IF (yn_varyslit ) THEN
          nslit = cali%nslit_sol(ix)
          slitwav(1:nslit) = cali%slitwav_sol(1:nslit, ix)
          slitfit = cali%slitfit_sol(:, :, :, ix)
        ELSE
          solwinfit(1:numwin,:,1 ) = cali%solwinfit(1:numwin, :,  ix)
        ENDIF
      ENDIF

      IF (wavcal) THEN
        IF (yn_varyslit .AND. (wavcal_sol .OR. which_slit == 5) ) THEN
          CALL radiance_wavcal_vary ( n_rad_wvl, &
               curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), error)
          cali%nwavcal_rad(ix) = nwavcal_rad
          cali%sswav_rad(1:nwavcal_rad, ix) = sswav_rad(1:nwavcal_rad)
        ELSE
          CALL radiance_wavcal (n_rad_wvl, &
               curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), error)
          cali%wincal_wav(1:numwin, ix) = wincal_wav(1:numwin)
          cali%radwinfit(1:numwin, 1:max_calfit_idx,  ix) = &
               radwinfit(1:numwin, 1:max_calfit_idx, 1)
        END IF
      ENDIF

      IF (calwrt) CALL write_cali_cross (calunit,ix, cali,.true.)
      IF (error) THEN
        pge_error_status = pge_errstat_warning
      ENDIF

      ! Save calibrated wavelengths
      diff(1:n_rad_wvl)  = curr_rad_spec(wvl_idx, 1:n_rad_wvl) - &
                           rad%wavl(1:n_rad_wvl, ix, iline)
          
      ! Check wavelength shifts, if it is too large, do not apply
      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        dwvl =  curr_rad_spec(wvl_idx, fidx+1) - curr_rad_spec(wvl_idx, fidx)
        IF (ANY(ABS(diff(fidx:lidx)) > dwvl)) diff(fidx:lidx) = 0.0
        fidx = lidx + 1
      ENDDO
      IF (MAXVAL(ABS(diff(1:n_rad_wvl))) == 0.0) RETURN

      wav (1:n_rad_wvl)  = rad%wavl(1:n_rad_wvl, ix, iline)
      DO i = 0, ntimes_loop - 1
        IF (rad%pix_errstat(ix, i) == pge_errstat_ok) THEN
          np = rad%nwav(ix, i)

          IF (np == n_rad_wvl) THEN
            rad%wavl(1:np, ix, i) = &
                 real(rad%wavl(1:np, ix, i) + diff(1:np) , kind=r4)
            IF ( ANY(rad%wavl(2:np, ix, i) - &
                 rad%wavl(1:np-1, ix, i) < 0.0) ) THEN
                 rad%wavl(1:np, ix, i) = &
                   real(rad%wavl(1:np, ix, i) - diff(1:np) , kind=r4)
              rad%pix_errstat(ix, i) = pge_errstat_warning
            ENDIF
          ELSE  ! Perform interpolation
            CALL INTERPOL(wav(1:n_rad_wvl), diff(1:n_rad_wvl), n_rad_wvl,  &
                 real(rad%wavl(1:np, ix, i), kind=8), corr(1:np), &
                 np, inter_errstat)
            IF (inter_errstat < 0)  THEN
              ! wavelength not calibrated
              rad%pix_errstat(ix, i) = pge_errstat_warning
            ELSE
              rad%wavl(1:np, ix, i) = &
                   real(rad%wavl(1:np, ix, i) + corr(1:np) , kind=r4)
              IF ( ANY(rad%wavl(2:np, ix, i) - &
                   rad%wavl(1:np-1, ix, i) < 0.0) ) THEN
                   rad%wavl(1:np, ix, i) = &
                   real(rad%wavl(1:np, ix, i) - corr(1:np), kind=r4)
                   rad%pix_errstat(ix, i) = pge_errstat_warning
              ENDIF
            ENDIF
          ENDIF
        ENDIF

      ENDDO
    ENDDO !first, last+pix:w

    RETURN
  END SUBROUTINE calibrate_rad_cross

  !SUBROUTINE read_cali_cross (cunit, ix, cali, radcal)
  !USE OMSAO_variables_module, ONLY: cali_group, slit_fname
  !IMPLICIT NONE
  !! INPUT/OUTPUT
  !INTEGER, INTENT(IN) :: cunit, ix
  !LOGICAL, INTENT(IN) :: radcal
  !TYPE(cali_group) :: cali
  !! LOGICAL variables
  !INTEGER :: npix
  ! OPEN(cunit, FILE=TRIM(ADJUSTL(slit_fname)), STATUS='unknown') 
  ! READ(cunit, '(5i5)') !which_slit, n_fitvar + 1
  ! READ(cunit, '(5i5)') npix !last_pix-first_pix + 1, first_pix, last_pix
  ! READ(cunit,'(10(a5))') !'wav', calfit_strings(fit_idx)
  !
  !END SUBROUTINE  

  SUBROUTINE write_cali_cross (cunit,ix, cali, radcal)
   USE OMSAO_variables_module, ONLY: yn_varyslit, numwin,&
                                  cali_group, fitspec_rad, the_pix, the_line, &
      n_rad_wvl, n_irrad_wvl, curr_rad_spec, curr_sol_spec
   IMPLICIT NONE
   ! INPUT/OUTPUT
   INTEGER, INTENT(IN) :: cunit, ix
   LOGICAL, INTENT(IN) :: radcal
   TYPE(cali_group) :: cali
   ! Local variables
   INTEGER :: iw, nw
   CHARACTER (len=20) :: form1, form2, form3
   !assign format
   WRITE(form1,'(a,i2,a)')  '(f8.2,',n_fitvar,'e15.5)'
   IF (radcal) THEN 
      nw = n_rad_wvl
   ELSE
      nw = n_irrad_wvl
   ENDIF
     WRITE(form2,'(a,i5,a)')  '(',700,'f8.3)'
     WRITE(form3,'(a,i5,a)')  '(',700,'e15.5)'
   !ELSE
    ! WRITE(form2,'(a,i5,a)')  '(',n_irrad_wvl,'f8.3)'
    ! WRITE(form3,'(a,i5,a)')  '(',n_irrad_wvl,'e15.5)'
   !ENDIF
   !write slit variables
   IF (.not. yn_varyslit) THEN 
       write(cunit, *) the_pix, the_line,  numwin, nw
       DO iw = 1, numwin
       IF (radcal) THEN 
           WRITE(cunit,ADJUSTL(TRIM(form1))) cali%wincal_wav(iw,ix), cali%radwinfit(iw, fit_idx(1:n_fitvar), ix)
       ELSE 
           WRITE(cunit,ADJUSTL(TRIM(form1))) cali%wincal_wav(iw, ix), cali%solwinfit(iw, fit_idx(1:n_fitvar), ix)
       ENDIF
       ENDDO
   ENDIF
   IF (.not. radcal) THEN
     WRITE(cunit,ADJUSTL(TRIM(form2))) curr_sol_spec(1,1:n_irrad_wvl)
     WRITE(cunit,ADJUSTL(TRIM(form3))) fitspec_rad(1:n_irrad_wvl)
     WRITE(cunit,ADJUSTL(TRIM(form3))) curr_sol_spec(2,1:n_irrad_wvl)
   ELSE
     WRITE(cunit,ADJUSTL(TRIM(form2))) curr_rad_spec(1,1:n_rad_wvl)
     WRITE(cunit,ADJUSTL(TRIM(form3))) fitspec_rad(1:n_rad_wvl)
     WRITE(cunit,ADJUSTL(TRIM(form3))) curr_rad_spec(2,1:n_rad_wvl)
   ENDIF
   RETURN
  RETURN
  END SUBROUTINE
end module m_cross_calibrate
