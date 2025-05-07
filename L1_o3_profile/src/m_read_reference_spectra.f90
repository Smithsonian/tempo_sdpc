!
module m_read_reference_spectra

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  IMPLICIT NONE

  PUBLIC read_reference_spectra, read_one_refspec1
  PRIVATE read_nc_refspec

CONTAINS  

  SUBROUTINE read_reference_spectra ( specunit, pge_error_status )
    USE OMSAO_indices_module,    ONLY: max_rs_idx, wvl_idx, spc_idx, &
         solar_idx
    USE OMSAO_parameters_module, ONLY: &
         zerospec_string, vb_lev_develop, maxwin
    USE OMSAO_variables_module,  ONLY: numwin, winlim, &
         refspec_fname, refspec_norm, n_refspec_pts, &
         refspec_orig_data, verb_thresh_lev !, scnwrt

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=27), PARAMETER :: modulename = 'read_reference_spectra:'
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
    INTEGER, PARAMETER :: maxline=354068
    INTEGER         :: i, max_pts
    REAL (kind=dp)  :: buffer
    REAL (kind=dp), DIMENSION(:,:),   ALLOCATABLE :: onespec
    REAL (kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: database
    TYPE ref_ctr_type
      INTEGER            :: nref
      REAL (KIND=dp), DIMENSION (maxwin, 2) :: winlim ! winwav of reference spectrum
    END TYPE ref_ctr_type
    TYPE (ref_ctr_type) :: ref_ctr
    ! =================================
    ! External OMI and Toolkit routines
    ! =================================
    INTEGER :: OMI_SMF_setmsg
    
    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    n_refspec_pts         = 0
    refspec_norm          = 1.0
    pge_error_status = pge_errstat_ok
    allocate(database (max_rs_idx, maxline, spc_idx))
    !-----------------------------
    ! set the boundary of refspec 
    !-----------------------------
    ref_ctr%nref  = 1
    ref_ctr%winlim(1, 1) = winlim(1, 1) 
    DO i = 2, numwin
      IF (winlim(i, 1) - winlim(i-1, 2) > 20.d0) THEN ! to skip the huge jump btw channels
       ref_ctr%winlim(ref_ctr%nref, 2) = winlim(i-1, 2) 
       ref_ctr%nref = ref_ctr%nref + 1
       ref_ctr%winlim(ref_ctr%nref, 1) = winlim(i, 1) 
      ENDIF
    ENDDO
    ref_ctr%winlim(ref_ctr%nref, 2) = winlim(numwin, 2) 
    ! -----------------------------------------------------------
    ! Read spectra one by one. Skip if name of file is ZEROSPEC
    ! -----------------------------------------------------------
    DO i = 1, max_rs_idx
      pge_error_status = pge_errstat_ok
      IF ( INDEX (TRIM(ADJUSTL(refspec_fname(i))),zerospec_string ) == 0 ) THEN
        buffer = 5.0
        IF ( i == solar_idx) THEN 
          buffer = 6.0
        ENDIF

        IF ( INDEX (TRIM(ADJUSTL(refspec_fname(i))), '.nc' ) /= 0 ) THEN
         CALL read_nc_refspec (refspec_fname(i), & 
                               ref_ctr%nref, ref_ctr%winlim(1:ref_ctr%nref, :), buffer, &  
                               n_refspec_pts(i),refspec_norm(i), onespec, pge_error_status)
        ELSE
         CALL read_one_refspec1 (specunit, refspec_fname(i), &
              ref_ctr%nref, ref_ctr%winlim(1:ref_ctr%nref, :), buffer, &  
              n_refspec_pts(i), refspec_norm(i), onespec, pge_error_status)
        ENDIF

        IF ( pge_error_status >= pge_errstat_error ) THEN 
          WRITE(*, *) TRIM(ADJUSTL(modulename))//TRIM(ADJUSTL(www_message))
          RETURN
        ENDIF
        IF ( n_refspec_pts(i) > maxline ) THEN 
          WRITE(*, *) TRIM(ADJUSTL(modulename))//':increase maxline>',n_refspec_pts(i),'!!!' ; STOP 1
        ENDIF
        WRITE(www_lun, '(A80,/,D14.6, I5,2F10.2)') refspec_fname(i), refspec_norm(i), & 
         n_refspec_pts(i), onespec(1, wvl_idx), onespec( n_refspec_pts(i), wvl_idx)
        IF (n_refspec_pts(i) == 0 ) cycle
        database(i,1:n_refspec_pts(i),1:spc_idx) = onespec(1:n_refspec_pts(i), 1:spc_idx)
        deallocate(onespec)
      ENDIF
    ENDDO

    ! ------------------------------------
    ! Report successful reading of spectra
    ! ------------------------------------
    IF ( verb_thresh_lev >= vb_lev_develop ) &
         errstat = OMI_SMF_setmsg (omsao_s_read_refspec_file, '', modulename, 0)
    !-----------------------------------------------
    ! end with assigning output variables and deallocating local variables
    !-----------------------------------------------
    max_pts = maxval(n_refspec_pts)
    allocate (refspec_orig_data(max_rs_idx, max_pts, 3))
    refspec_orig_data = 0.0D0
    refspec_orig_data (:,:,1:2) = database(:, 1:max_pts, 1:2)
    deallocate(database)
    RETURN
  END SUBROUTINE read_reference_spectra

  SUBROUTINE read_one_refspec1 ( specunit, specname, & 
     nwin, wlim, buffer, &
     nspec,specnorm, onespec, pge_error_status )

    USE OMSAO_precision_module,   ONLY: dp
    USE OMSAO_indices_module,     ONLY: wvl_idx, spc_idx
    USE OMSAO_parameters_module,  ONLY: maxchlen, lm_start_of_table
    USE OMSAO_errstat_module
    use m_utilities, only: skip_to_filemark
    IMPLICIT NONE

    ! ----------------
    ! Input Parameters
    ! ----------------
    INTEGER,             INTENT (IN) :: specunit
    CHARACTER (LEN=*),   INTENT (IN) :: specname
    INTEGER, INTENT(IN) :: nwin
    REAL (KIND=dp), DIMENSION (nwin, 2) :: wlim
    REAL (KIND=dp) :: buffer
    ! -----------------
    ! Output Parameters
    ! -----------------
    INTEGER,                         INTENT (INOUT) :: pge_error_status
    INTEGER,                         INTENT (OUT) :: nspec
    REAL (KIND=dp),                  INTENT (OUT) :: specnorm
    REAL (KIND=dp), DIMENSION (:,:), INTENT (OUT), ALLOCATABLE :: onespec

    ! ----------------
    ! Local Variables
    ! ----------------
    INTEGER  :: i, ip, file_read_stat, imin, imax, iwin, fwin, lwin,wstep
    INTEGER,        DIMENSION (:), ALLOCATABLE   :: irev
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE   :: x, y
    REAL (KIND=dp), DIMENSION (2)            :: specwav
    REAL (KIND=dp)                           :: xmin, xmax
    CHARACTER (LEN=maxchlen)                 :: lastline
    LOGICAL                                  :: isinc

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'read_one_refspec'

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER :: errstat

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER :: OMI_SMF_setmsg

    nspec = 0 ; specnorm = 0.0 ; specwav = 0.0

    ! -----------------------
    ! Open reference spectrum
    ! -----------------------
    OPEN ( UNIT=specunit, FILE=TRIM(ADJUSTL(specname)), STATUS='OLD', IOSTAT=errstat)
    IF ( errstat /= pge_errstat_ok ) THEN
       errstat = OMI_SMF_setmsg ( &
        omsao_e_open_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
        www_message='Not exist:'//TRIM(ADJUSTL(specname))  
       pge_error_status = pge_errstat_error; RETURN
    END IF

    ! --------------------------------------
    ! Skip comments header to start of table
    ! --------------------------------------
    CALL skip_to_filemark ( specunit, lm_start_of_table, lastline, file_read_stat)
    IF ( file_read_stat /= file_read_ok ) THEN
          errstat = OMI_SMF_setmsg ( &
          omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
          www_message='Erros in reading :'//TRIM(ADJUSTL(specname))  
          pge_error_status = pge_errstat_error; RETURN
    END IF

    ! -----------------------------------------------
    ! Read dimension, start&end, and norm of spectrum
    ! -----------------------------------------------
    READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) nspec, specwav, specnorm
    IF ( file_read_stat /= file_read_ok ) THEN
           errstat = OMI_SMF_setmsg ( &
          omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
          www_message='Erros in reading dimension :'//TRIM(ADJUSTL(specname))  
          pge_error_status = pge_errstat_error; RETURN
    END IF

    allocate(x(nspec), y(nspec), irev(nspec))
    x = 0.0 ; y = 0.0 ; irev = 0
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

    IF (xmin > wlim(nwin,2)+buffer .OR. xmax < wlim(1,1)-buffer) THEN
       www_message= 'does not cover any wavelength range!!!'
       pge_error_status = pge_errstat_warning; RETURN
    ENDIF

    ip = 1; i = 1
    IF (isinc) THEN
       fwin = 1; lwin = nwin; wstep=1
    ELSE
       fwin = nwin; lwin = 1; wstep=-1
    ENDIF

    DO iwin = fwin, lwin, wstep
       xmin = wlim(iwin, 1)-buffer; xmax = wlim(iwin, 2)+buffer
       DO WHILE (i <= nspec )
          READ (UNIT=specunit, FMT=*, IOSTAT=file_read_stat) x(ip), y(ip)
          i = i + 1
          IF ( file_read_stat /= file_read_ok ) THEN
                errstat = OMI_SMF_setmsg ( &
                omsao_e_read_refspec_file, TRIM(ADJUSTL(specname)), modulename, 0)
          www_message='Erros in reading dimension :'//TRIM(ADJUSTL(specname))  
           pge_error_status = pge_errstat_error; RETURN
          ENDIF
          IF (x(ip) >= xmin .AND. x(ip) <= xmax) THEN
             ip = ip + 1
          ELSE IF (isinc .AND. x(ip) > xmax) THEN
             IF (iwin < lwin ) THEN 
                IF (x(ip) > wlim(iwin+1, 1)) ip = ip + 1
             ENDIF
             EXIT
          ELSE IF (.NOT. isinc .AND. x(ip) < xmin) THEN
             IF (iwin < lwin ) THEN 
                IF ( x(ip) < wlim(iwin+1, 2)) ip = ip + 1
             ENDIF
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
    allocate (onespec(nspec, spc_idx))
    onespec(1:nspec,wvl_idx) = x(1:nspec)
    onespec(1:nspec,spc_idx) = y(1:nspec)
    specwav(1:2)             = (/ x(1), x(nspec) /)

    !-------------------------------
    ! end with deallocate
    !-------------------------------
    deallocate(x, y, irev)
    RETURN
  END SUBROUTINE read_one_refspec1

  SUBROUTINE read_nc_refspec (specname, nwin, wlim, buffer, nspec, specnorm, onespec, pge_error_status) 
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    ! ----------------
    ! Input Parameters
    ! ----------------
    CHARACTER (LEN=*),   INTENT (IN) :: specname
    INTEGER, INTENT(IN) :: nwin
    REAL (KIND=dp), DIMENSION (nwin, 2) :: wlim
    REAL (KIND=dp) :: buffer
    ! -----------------
    ! Output Parameters
    ! -----------------
    INTEGER,                         INTENT (INOUT) :: pge_error_status
    INTEGER,                         INTENT (OUT)   :: nspec
    REAL (KIND=dp),                  INTENT (OUT)   :: specnorm
    REAL (KIND=dp), DIMENSION (:,:), INTENT (OUT), ALLOCATABLE :: onespec
    
    ! -----------------
    ! Local variables
    ! -----------------
    LOGICAL :: isinc
    INTEGER :: ncid, vid, rcode, wmx, fidx, lidx,iwin, ntmp, imin, imax
    INTEGER, DIMENSION (:,:), ALLOCATABLE :: winpix
    REAL (KIND=dp) :: xmin, xmax
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: specwav
    CHARACTER (LEN=14)    :: tmpchar

    ! Open this file
    ncid = ncopn(trim(adjustl(specname)), nf_Nowrite, rcode)
    IF (rcode .eq. -1 ) THEN 
      www_message=TRIM(ADJUSTL(specname))//' is not exist'
      pge_error_status = pge_errstat_error ; RETURN
    ENDIF

    ! Read specnorm
    rcode = nf_inq_varid(ncid, 'Specnorm', vid)
    rcode = nf_get_var_double(ncid, vid, specnorm)

    ! Read Wavelengths
    rcode = nf_inq_varid(ncid, 'Wavelength', vid)
    rcode = nf_inq_dim(ncid, 1, tmpchar, wmx)
    allocate (specwav(wmx))
    rcode = nf_get_var_double(ncid, vid, specwav)

    IF (specwav(1) < specwav(2)) THEN
       xmin = specwav(1); xmax = specwav(wmx)
       imin = 1; imax = wmx; isinc = .TRUE.
    ELSE
       xmin = specwav(wmx); xmax = specwav(1)
       imin = wmx; imax = 1; isinc = .FALSE.
    ENDIF
     
    IF (.NOT. isinc) THEN 
       www_message= 'should be increasing order !!!'
       pge_error_status = pge_errstat_error; RETURN
    ENDIF

    IF (xmin > wlim(nwin,2)+buffer .OR. xmax < wlim(1,1)-buffer) THEN
       www_message= 'does not cover any wavelength range!!!'
       pge_error_status = pge_errstat_warning; RETURN
    ENDIF

    ! Find position of selected wavrange
    allocate (winpix(nwin, 2))  
    winpix = 0
    nspec = 0
    DO iwin = 1, nwin
      fidx = MINVAL(MINLOC(specwav, MASK=(specwav >= wlim(iwin,1)-buffer )))
      lidx = MINVAL(MAXLOC(specwav, MASK=(specwav <= wlim(iwin,2)+buffer )))
      !if (wlim(iwin, 2) > specwav(wmx)) cycle !xl,12/28/2021 move to end of loop
      IF (fidx <= 0 .OR. lidx <= 0) CYCLE
      winpix(iwin, 1:2) = (/fidx, lidx/)
      nspec = nspec + lidx - fidx + 1
      IF (wlim(iwin, 2) > specwav(wmx)) CYCLE
    ENDDO
    allocate (onespec(nspec, 2))

    ! Read spectrum 
    fidx = 1
    DO iwin = 1, nwin 
      ntmp = winpix(iwin, 2) - winpix(iwin,1) + 1
      if (ntmp <= 1) cycle
      lidx = fidx + ntmp -1
      onespec(fidx:lidx, 1) = specwav(winpix(iwin,1):winpix(iwin,2))
      rcode = nf_inq_varid(ncid, 'Spectrum', vid)
      rcode = nf_get_vara_double (ncid, vid, (/winpix(iwin,1)/), (/ntmp/), onespec(fidx:lidx,2))
      fidx = lidx + 1
    ENDDO
    !-------------------------------
    ! end with deallocate
    !-------------------------------
    deallocate (specwav, winpix)
  END SUBROUTINE
end module m_read_reference_spectra
