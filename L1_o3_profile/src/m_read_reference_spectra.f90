!
module m_read_reference_spectra

  public read_reference_spectra
  private read_one_refspec

contains

  SUBROUTINE read_reference_spectra ( specunit, pge_error_status )

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: max_rs_idx, wvl_idx, spc_idx, &
         solar_idx!, pge_static_input_luns, so2_idx
    USE OMSAO_parameters_module, ONLY: max_spec_pts, &
         zerospec_string, vb_lev_develop!, maxchlen
    USE OMSAO_variables_module,  ONLY: &
         refspec_fname, refspec_firstlast_wav, refspec_norm, n_refspec_pts, &
         refspec_orig_data, verb_thresh_lev, winwav_min, winwav_max, scnwrt, &
         solar_refspec
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=27), PARAMETER :: modulename = 'read_reference_spectra'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER, INTENT (IN) :: specunit

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER, INTENT (OUT) :: pge_error_status

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER         :: i
    real (kind=dp), dimension(2) :: tmp_specwav
    real (kind=dp), dimension(1:max_spec_pts,wvl_idx:spc_idx) :: tmp_onespec

    ! =================================
    ! External OMI and Toolkit routines
    ! =================================
    INTEGER :: OMI_SMF_setmsg

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    n_refspec_pts         = 0
    refspec_norm          = 0.0
    refspec_orig_data     = 0.0
    refspec_firstlast_wav = 0.0
    pge_error_status = pge_errstat_ok

    ! -----------------------------------------------------------
    ! Read spectra one by one. Skip if name of file is ZEROSPEC
    ! -----------------------------------------------------------
    DO i = 1, max_rs_idx
      IF( scnwrt ) WRITE(*,*) i, TRIM(ADJUSTL(refspec_fname(i))) !! Kai
      IF ( INDEX (TRIM(ADJUSTL(refspec_fname(i))), &
           zerospec_string ) == 0 ) THEN
!  FIXME
!  Modified to mask array temporaries
        CALL read_one_refspec1 ( &
             specunit, &
             refspec_fname(i), winwav_min, winwav_max, n_refspec_pts(i), &
!             refspec_firstlast_wav(i,wvl_idx:spc_idx), refspec_norm(i), &
!             refspec_orig_data(i,1:max_spec_pts,wvl_idx:spc_idx), &
             tmp_specwav, refspec_norm(i), &
             tmp_onespec, &
             pge_error_status)
        refspec_firstlast_wav(i,wvl_idx:spc_idx) = tmp_specwav
        refspec_orig_data(i,1:max_spec_pts,wvl_idx:spc_idx) = tmp_onespec

        IF( scnwrt ) THEN !! start Kai add
          WRITE(www_lun, '(A80,/,2F10.2,I5,2F10.2, D14.6)') refspec_fname(i), winwav_min, &
               winwav_max, n_refspec_pts(i), refspec_orig_data(i, 1, wvl_idx), &
               refspec_orig_data(i, n_refspec_pts(i), wvl_idx), refspec_norm(i)
        ENDIF !! END Kai add               
        IF ( pge_error_status >= pge_errstat_error ) RETURN

        !IF (i == so2_idx) refspec_orig_data(i,1:max_spec_pts,wvl_idx) = refspec_orig_data(i,1:max_spec_pts,wvl_idx)-0.269
      END IF
    END DO

    solar_refspec(1:n_refspec_pts(solar_idx)) = refspec_orig_data(solar_idx, 1:n_refspec_pts(solar_idx), spc_idx)

    ! ------------------------------------
    ! Report successful reading of spectra
    ! ------------------------------------
    IF ( verb_thresh_lev >= vb_lev_develop ) &
         errstat = OMI_SMF_setmsg (omsao_s_read_refspec_file, '', modulename, 0)

    RETURN
  END SUBROUTINE read_reference_spectra

  SUBROUTINE read_one_refspec ( specunit, specname, winwav_min, winwav_max,&
       nspec, specwav, specnorm, onespec, pge_error_status )

    USE OMSAO_precision_module,   ONLY: dp
    USE OMSAO_indices_module,     ONLY: wvl_idx, spc_idx
    USE OMSAO_parameters_module,  ONLY: maxchlen, max_spec_pts, &
         lm_start_of_table!, vb_lev_omidebug
    !USE OMSAO_variables_module,   ONLY: verb_thresh_lev
    USE OMSAO_errstat_module
    use utilities, only: skip_to_filemark

    IMPLICIT NONE

    ! ----------------
    ! Input Parameters
    ! ----------------
    INTEGER,             INTENT (IN) :: specunit
    CHARACTER (LEN=*),   INTENT (IN) :: specname
    REAL      (KIND=dp), INTENT (IN) :: winwav_min, winwav_max

    ! -----------------
    ! Output Parameters
    ! -----------------
    INTEGER,       INTENT (OUT) :: pge_error_status
    INTEGER,       INTENT (OUT) :: nspec
    REAL(KIND=dp), INTENT (OUT) :: specnorm
    REAL(KIND=dp), DIMENSION (2), INTENT (OUT) :: specwav
    REAL(KIND=dp), DIMENSION (max_spec_pts,wvl_idx:spc_idx), INTENT (OUT) :: onespec

    ! ----------------
    ! Local Variables
    ! ----------------
    INTEGER  :: i, file_read_stat, j1, j2, imin, imax!, ios
    INTEGER,        DIMENSION (max_spec_pts) :: irev
    REAL (KIND=dp), DIMENSION (max_spec_pts) :: x, y
    REAL (KIND=dp)                           :: xdum, xmin, xmax
    CHARACTER (LEN=maxchlen)                 :: lastline
    LOGICAL                                  :: isinc

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'read_one_refspec'

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER :: errstat!, version

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER :: OMI_SMF_setmsg

    nspec = 0 
    specnorm = 0.0 
    specwav = 0.0 
    onespec = 0.0
    pge_error_status = pge_errstat_ok

    x = 0.0 
    y = 0.0


    ! -----------------------
    ! Open reference spectrum
    ! -----------------------
    OPEN ( UNIT=specunit, FILE=TRIM(ADJUSTL(specname)), STATUS='OLD', IOSTAT=errstat)
    IF ( errstat /= pge_errstat_ok ) THEN
      errstat = OMI_SMF_setmsg ( &
           omsao_e_open_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! --------------------------------------
    ! Skip comments header to start of table
    ! --------------------------------------
    CALL skip_to_filemark ( specunit, lm_start_of_table, lastline, file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg ( &
           omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! -----------------------------------------------
    ! Read dimension, start&end, and norm of spectrum
    ! -----------------------------------------------
    READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) nspec, specwav, specnorm
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg ( &
           omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF


    ! ---------------------------------------------------------------------
    ! Find first and last index to read, based on WINWAV_MIN and WINWAV_MAX
    ! ---------------------------------------------------------------------
    IF (specwav(1) < specwav(2)) THEN
      xmin = specwav(1)
      xmax = specwav(2)
      imin = 1
      imax = nspec
      isinc = .TRUE.
    ELSE
      xmin = specwav(2)
      xmax = specwav(1)
      imin = nspec
      imax = 1
      isinc = .FALSE.
    ENDIF

    IF (xmin > winwav_max .OR. xmax < winwav_min) THEN
      WRITE(www_lun, *) TRIM(ADJUSTL(specname)), ' does not cover any wavelength range!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ELSE IF (xmin > winwav_min .AND. xmax < winwav_max) THEN    ! partly covered
      j1 = imin
      j2 = imax
    ELSE IF (xmin > winwav_min) THEN    ! partly covered
      j1 = imin
      j2 = 0
    ELSE IF (xmax < winwav_max) THEN    ! partly covered
      j2 = imax
      j1 = 0
    ELSE
      j1 = 0
      j2 = 0                    ! fully covered
    ENDIF

    IF (j1 < 1 .OR. j2 < 1) THEN
      getidx: DO i = 1, nspec
        READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) xdum
        IF ( i /= nspec .AND. file_read_stat /= file_read_ok ) THEN
          errstat = OMI_SMF_setmsg ( &
               omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        END IF

        IF (j1 == 0) THEN
          IF (isinc) THEN
            IF ( xdum > winwav_min ) j1 = i - 1
          ELSE
            IF ( xdum < winwav_min ) j1 = i 
          ENDIF
        ENDIF

        IF (j2 == 0) THEN
          IF (isinc) THEN
            IF ( xdum > winwav_max ) j2 = i 
          ELSE
            IF ( xdum < winwav_max ) j2 = i - 1
          ENDIF
        ENDIF
        IF ( j1 > 0 .AND. j2 > 0 ) EXIT getidx

      END DO getidx

      ! --------------------------------------------------------------------------
      ! Rewind file, skip to start of table line, and re-read first entry as dummy
      ! --------------------------------------------------------------------------
      REWIND ( specunit )
      ! --------------------------------------
      ! Skip comments header to start of table
      ! --------------------------------------
      CALL skip_to_filemark ( specunit, lm_start_of_table, lastline, file_read_stat )
      IF ( file_read_stat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg ( &
             omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF

      ! -----------------------------------------------
      ! Read dimension, start&end, and norm of spectrum
      ! -----------------------------------------------
      READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) nspec, specwav, specnorm
      IF ( file_read_stat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg ( &
             omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF
    ENDIF

    IF (j2 -j1 + 1 > max_spec_pts) THEN
      WRITE(www_lun, *) TRIM(ADJUSTL(specname)), ', increase max_spec_pts ', max_spec_pts, &
           ' to at least ', j2 - j1 + 1
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! -------------------------------------------------
    ! Skip first J1-1 lines, then read J1 to J2 entries
    ! -------------------------------------------------
    DO i = 1, j2
      IF (i < j1) THEN
        READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) xdum  ! skip these lines
      ELSE
        READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) x(i-j1+1), y(i-j1+1)
      ENDIF

      IF ( i /= nspec .AND. file_read_stat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg ( &
             omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF
    END DO

    ! -----------------------------------------------
    ! Close fitting control file, report SUCCESS read
    ! -----------------------------------------------
    CLOSE ( UNIT=specunit )

    ! ------------------------------------------------------------
    ! Reassign number of spectral points and first/last wavelength
    ! ------------------------------------------------------------
    nspec = j2 - j1 + 1

    ! ---------------------------------------------------
    ! Reorder spectrum so that wavelengths are increasing
    ! ---------------------------------------------------
    IF ( .NOT. isinc ) THEN
      irev = (/ (i, i = nspec, 1, -1) /)
      x(1:nspec) = x(irev(1:nspec))
      y(1:nspec) = y(irev(1:nspec))
    END IF

    ! ------------------------
    ! Assign output quantities
    ! ------------------------
    onespec(1:nspec,wvl_idx) = x(1:nspec)
    onespec(1:nspec,spc_idx) = y(1:nspec)
    specwav(1:2)             = (/ x(1), x(nspec) /)

    !WRITE(www_lun, *) TRIM(ADJUSTL(specname))
    !WRITE(www_lun, *) x(1), x(nspec), winwav_min, winwav_max

    RETURN
  END SUBROUTINE read_one_refspec

  !xliu: 10/25/2011, replace the routine above with the following one
  !      This subroutine can allow data gap in the reference if that portion of
  !      the orbit is not used so that we could use a smaller max_spec_pts
  SUBROUTINE read_one_refspec1 ( specunit, specname, winwav_min, winwav_max,&
     nspec, specwav, specnorm, onespec, pge_error_status )

    USE OMSAO_precision_module,   ONLY: dp
    USE OMSAO_indices_module,     ONLY: wvl_idx, spc_idx
    USE OMSAO_parameters_module,  ONLY: maxchlen, max_spec_pts, lm_start_of_table,&
        vb_lev_omidebug, maxwin
    USE OMSAO_variables_module,   ONLY: verb_thresh_lev, winlim, numwin
    USE OMSAO_errstat_module
    use utilities, only: skip_to_filemark
    IMPLICIT NONE

    ! ----------------
    ! Input Parameters
    ! ----------------
    INTEGER,             INTENT (IN) :: specunit
    CHARACTER (LEN=*),   INTENT (IN) :: specname
    REAL      (KIND=dp), INTENT (IN) :: winwav_min, winwav_max
    ! -----------------
    ! Output Parameters
    ! -----------------
    INTEGER,                                                  INTENT (OUT) :: pge_error_status
    INTEGER,                                                  INTENT (OUT) :: nspec
    REAL (KIND=dp),                                           INTENT (OUT) :: specnorm
    REAL (KIND=dp), DIMENSION (2),                            INTENT (OUT) :: specwav
    REAL (KIND=dp), DIMENSION (max_spec_pts,wvl_idx:spc_idx), INTENT (OUT) :: onespec

    ! ----------------
    ! Local Variables
    ! ----------------
    INTEGER  :: i, ip, ios, file_read_stat, imin, imax, nwin, iwin, fwin, lwin,wstep
    INTEGER,        DIMENSION (max_spec_pts) :: irev
    REAL (KIND=dp), DIMENSION (max_spec_pts) :: x, y
    REAL (KIND=dp), DIMENSION (maxwin, 2)    :: wlim
    REAL (KIND=dp)                           :: xdum, xmin, xmax
    CHARACTER (LEN=maxchlen)                 :: lastline
    LOGICAL                                  :: isinc

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'read_one_refspec'

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER :: errstat, version

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER :: OMI_SMF_setmsg

    IF (numwin == 1) THEN
       wlim(1, 1) = winlim(1, 1) - 5.0d0
       wlim(1, 2) = winlim(1, 2) + 5.0d0
       nwin = 1
    ELSE
       iwin = 1
       wlim(1, 1) = winlim(1, 1) - 5.0d0
       DO i = 2, numwin
          IF (winlim(i, 1) - winlim(i-1, 2) > 20.d0) THEN
             wlim(iwin, 2) = winlim(i-1, 2) + 5.d0
             iwin = iwin + 1
             wlim(iwin, 1) = winlim(i, 1) - 5.d0
          ENDIF
       ENDDO
       wlim(iwin, 2) = winlim(numwin, 2) + 5.0d0
       nwin = iwin
    ENDIF

    nspec = 0 ; specnorm = 0.0 ; specwav = 0.0 ; onespec = 0.0
    pge_error_status = pge_errstat_ok
    x = 0.0 ; y = 0.0
    ! -----------------------
    ! Open reference spectrum
    ! -----------------------
    OPEN ( UNIT=specunit, FILE=TRIM(ADJUSTL(specname)), STATUS='OLD', IOSTAT=errstat)
    IF ( errstat /= pge_errstat_ok ) THEN
       errstat = OMI_SMF_setmsg ( &
          omsao_e_open_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
         print * , 'no file', specname
       pge_error_status = pge_errstat_error; RETURN
    END IF

    ! --------------------------------------
    ! Skip comments header to start of table
    ! --------------------------------------
    CALL skip_to_filemark ( specunit, lm_start_of_table, lastline, file_read_stat)
    IF ( file_read_stat /= file_read_ok ) THEN
          errstat = OMI_SMF_setmsg ( &
          omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
          pge_error_status = pge_errstat_error; RETURN
    END IF

    ! -----------------------------------------------
    ! Read dimension, start&end, and norm of spectrum
    ! -----------------------------------------------
    READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) nspec, specwav, specnorm
    IF ( file_read_stat /= file_read_ok ) THEN
           errstat = OMI_SMF_setmsg ( &
          omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
          pge_error_status = pge_errstat_error; RETURN
    END IF

    ! ---------------------------------------------------------------------
    ! Find first and last index to read, based on WINWAV_MIN and WINWAV_MAX
    ! ---------------------------------------------------------------------
    IF (specwav(1) < specwav(2)) THEN
       xmin = specwav(1); xmax = specwav(2)
       imin = 1; imax = nspec; isinc = .TRUE.
    ELSE
       xmin = specwav(2); xmax = specwav(1)
       imin = nspec; imax = 1; isinc = .FALSE.
    ENDIF

    IF (xmin > winwav_max .OR. xmax < winwav_min) THEN
       WRITE(*, *) TRIM(ADJUSTL(specname)), ' does not cover any wavelength range!!!'
       pge_error_status = pge_errstat_warning; RETURN
    ENDIF

    ip = 1; i = 1
    IF (isinc) THEN
       fwin = 1; lwin = nwin; wstep=1
    ELSE
       fwin = nwin; lwin = 1; wstep=-1
    ENDIF

    DO iwin = fwin, lwin, wstep
       xmin = wlim(iwin, 1); xmax = wlim(iwin, 2)

       DO WHILE (i <= nspec )
          READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) x(ip), y(ip)
          i = i + 1
          IF ( file_read_stat /= file_read_ok ) THEN
                errstat = OMI_SMF_setmsg ( &
                omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
           pge_error_status = pge_errstat_error; RETURN
          ENDIF
          IF (x(ip) >= xmin .AND. x(ip) <= xmax) THEN
             IF (ip > max_spec_pts) THEN
                WRITE(*, *) TRIM(ADJUSTL(specname)), ', increase max_spec_pts ', max_spec_pts
                pge_error_status = pge_errstat_error; RETURN
             ENDIF
             ip = ip + 1
          ELSE IF (isinc .AND. x(ip) > xmax) THEN
             IF (iwin < lwin .AND. x(ip) > wlim(iwin+1, 1)) ip = ip + 1
             EXIT
          ELSE IF (.NOT. isinc .AND. x(ip) < xmin) THEN
             IF (iwin < lwin .AND. x(ip) < wlim(iwin+1, 2)) ip = ip + 1
             EXIT
          ENDIF
       ENDDO
    ENDDO

    ! -----------------------------------------------
    ! Close fitting control file, report SUCCESS read
    ! -----------------------------------------------
    CLOSE ( UNIT=specunit )

    ! ------------------------------------------------------------
    ! Reassign number of spectral points and first/last wavelength
    ! ------------------------------------------------------------
    nspec = ip - 1

    ! ---------------------------------------------------
    ! Reorder spectrum so that wavelengths are increasing
    ! ---------------------------------------------------
    IF ( .NOT. isinc ) THEN
       irev = (/ (i, i = nspec, 1, -1) /)
       x(1:nspec) = x(irev(1:nspec));  y(1:nspec) = y(irev(1:nspec))
    END IF

    ! ------------------------
    ! Assign output quantities
    ! ------------------------
    onespec(1:nspec,wvl_idx) = x(1:nspec)
    onespec(1:nspec,spc_idx) = y(1:nspec)
    specwav(1:2)             = (/ x(1), x(nspec) /)

    RETURN
  END SUBROUTINE read_one_refspec1
end module m_read_reference_spectra
