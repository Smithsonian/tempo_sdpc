!
module tmpo_read_l1b_data
  USE tell_module
  USE tio_module
  USE m_read_l1_tio, ONLY: read_L1_dims_tio,  open_L1_tio, close_l1_tio, &
      read_l1_rad_line_tio
  USE m_convert_coadd, ONLY: convert_2bytes_to_16bits, coadd_2bytes_qflgs,prespec_align, solwavcal_coadd, radwavcal_coadd
  USE m_fitting_util, ONLY: reduce_rad_resolution,reduce_irrad_resolution
  USE OMSAO_indices_module, ONLY: spc_idx
  USE OMSAO_parameters_module, ONLY: max_ring_pts
  USE OMSAO_precision_module
  USE OMSAO_tmpodata_module, ONLY: nwavel_max, nwavel_ccd,nxtrack_max,ntimes_max,&
      irrad_swathname, rad_swathname,tmpo_rad, tmpo_irrad, tmpo_ring, tmpo_refl,tmpo_geo=>tmpo_geo1


  REAL (KIND=r4), DIMENSION (spc_idx, max_ring_pts, nxtrack_max), PRIVATE ::tmpo_solspec_ring

  
  public tmpo_read_irradiance, tmpo_read_radiance_lines, & 
         tmpo_set_parameters
  private
  CONTAINS
  
  SUBROUTINE tmpo_set_parameters (errstat)
  USE OMSAO_variables_module, ONLY: l1b_rad_filename, &
      band_selectors, inschs, nswath, nxtrack, ntimes, &
       GranuleJDay, GranuleYear, GranuleMonth, GranuleDay
  USE OMSAO_tmpodata_module, ONLY: rad_swathname
  USE m_read_metadata_tio, ONLY:read_date_tio

  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: errstat
  INTEGER ::  nwavel
  

  ! Set up dimensions
  CALL read_l1_dims_tio (l1b_rad_filename, rad_swathname(1), ntimes, nxtrack, nwavel, errstat)
  IF (errstat /= 0) RETURN
 
  ! Set up instrument channels
  IF (nswath == 2) THEN 
    inschs(1) = 1
    inschs(2) = 2
  ELSE IF (nswath == 1) THEN 
    IF (band_selectors(1) == 1) THEN
        inschs(1) = 1
    ELSE
        inschs(1) = 2
    ENDIF
  ELSE
    WRITE(*,'(A)') 'Check nswath !!!'
    stop 1
  ENDIF 

  ! Set up dates from filename
  CALL read_date_tio (l1b_rad_filename, GranuleYear, GranuleMOnth, GranuleDay, GranuleJDay, errstat)
  IF (errstat /= 0) RETURN

  RETURN
  END SUBROUTINE

  SUBROUTINE tmpo_read_irradiance ( first_pix, last_pix, pge_error_status)
   USE OMSAO_indices_module, ONLY : sig_idx, wvl_idx, spc_idx
   USE OMSAO_parameters_module, ONLY: maxchlen, maxwin, max_fit_pts, mswath,max_ring_pts
   USE OMSAO_precision_module
   USE OMSAO_errstat_module
   USE OMSAO_variables_module, ONLY: nswath, nxtrack, numwin, currpix, &
    winlim, lower_spec, upper_spec, inschs, band_selectors, &
    retlbnd, retubnd, reduce_lbnd, reduce_ubnd,redlam, redsampr,redslw, & 
    reduce_slit, reduce_resolution, &
    l1b_irrad_filename, nxbin, &
    use_backup, refdbdir, calunit,&
    scnwrt, use_redfixwav, which_slit, wcal_bef_coadd, dwavmax,GranuleJDay
   USE ozprof_data_module, ONLY:toms_fwhm, pos_alb,nrefl, mrefl
   USE m_gauss, ONLY:gauss_uneven
   USE m_triangle, ONLY:triangle_uneven
   IMPLICIT NONE
   
   !----------------
   ! Input/output variables
   !----------------
   INTEGER, INTENT (IN) :: first_pix, last_pix
   INTEGER, INTENT (INOUT) ::  pge_error_status
   ! type variables
   TYPE (tiof_file_type) :: tio_l1obj
   ! local variables
   LOGICAL :: read_irrad=.true.
   LOGICAL :: problems
   INTEGER :: n_ring_off = 12
   INTEGER :: i, j,fidx,lidx, is, iline,ch,iw,ix,iix, irefl, errstat, idx, ic
   INTEGER :: nx, nt,nw, nsub, nbin, nring, nbad, nwavel, noff1, noff2, thedoy
   ! variable used to read TEMPO data
   INTEGER (KIND=i2) :: irrad_mflg
   INTEGER (KIND=i4), DIMENSION(mswath)  :: nwls 
   INTEGER (KIND=i4), DIMENSION(mswath)  :: spos, epos  ! pos btw UV and VIS
   INTEGER (KIND=i4), DIMENSION(:), ALLOCATABLE :: idxs
   INTEGER (KIND=i2), DIMENSION(:, :), ALLOCATABLE :: irrad_qflg
   REAL (kind=4),DIMENSION (:, :), ALLOCATABLE :: irrad_spec, irrad_prec, irrad_wavl
   ! variables used for reduced resultion
   INTEGER :: npos, np
   REAL (KIND=dp) :: tmpsampr, retswav, retewav
   REAL (KIND=dp),   DIMENSION(:,:,:), ALLOCATABLE :: tmpspec
   INTEGER(KIND=i2), DIMENSION(:, :), ALLOCATABLE :: tmpqflg
   ! Subset variables
   INTEGER, DIMENSION(maxwin,2) :: winpix ! pos btw winlim
   INTEGER (KIND=i4), DIMENSION (maxwin) :: nwbin
   INTEGER, PARAMETER :: nbits = 16
   INTEGER (KIND=i2), DIMENSION(:), ALLOCATABLE :: flgmsks
   INTEGER (kind=i2), DIMENSION(:, :,:), ALLOCATABLE :: flgbits
   REAL (KIND=dp) :: wcenter, normsc
   REAL (KIND=dp), DIMENSION(nxbin) :: wshis, wsqus
   REAL (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: subspec
   REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: subring
   ! Error Message
   CHARACTER(LEN=maxchlen) :: message, bkfname
   CHARACTER(LEN=*), PARAMETER :: modulename='tmpo_read_irradiance_data'

   
   !--------------------------------------------------------------------------
   ! Starting with allocating local variables
   !--------------------------------------------------------------------------
   allocate (irrad_qflg (nwavel_max, nxtrack_max))
   allocate (irrad_prec (nwavel_max, nxtrack_max))
   allocate (irrad_spec (nwavel_max, nxtrack_max))
   allocate (irrad_wavl (nwavel_max, nxtrack_max))
   allocate (flgmsks(nwavel_max), idxs(nwavel_ccd))
   allocate (flgbits(nxbin, nwavel_max, 0:nbits-1))
   allocate (tmpspec (sig_idx, nwavel_max, nxtrack_max))
   allocate (tmpqflg(nwavel_max, nxtrack_max))
   allocate (subspec(nxbin, sig_idx, max_fit_pts))
   allocate (subring(sig_idx, max_ring_pts))
   ! ----------------------------
   ! Initialize irradiance arrays
   ! ----------------------------
   tmpo_irrad%errstat(:) = pge_errstat_ok
   tmpo_irrad%nwav (:) = 0
   tmpo_irrad%npix (:,:) = 0
   tmpo_irrad%prec (:,:) = 0.0
   tmpo_irrad%spec (:,:) = 0.0
   tmpo_irrad%wavl (:,:) = 0.0
   tmpo_irrad%qflg (:,:) = 0
   ! ----------------------------------------------------
   ! read irradiance : arracy is merged with UV and VIS
   ! ----------------------------------------------------
   fidx = 1
   iline = 0
   nwavel = 0
   errstat = 0
   IF (.NOT. use_backup) THEN 
     DO is = 1, nswath
        ch = inschs(is)
        call read_L1_dims_tio (l1b_irrad_filename, irrad_swathname(ch),&
             nt, nx,nw, errstat = errstat)
        call open_L1_tio (l1b_irrad_filename, tio_l1obj, errstat)
       
        if (errstat /= 0) then
           message = ADJUSTL(TRIM(modulename))//": failed to open irradiance file"
           go to 123
        endif
        lidx = fidx + nw -1
        call read_L1_rad_line_tio (tio_l1obj, irrad_swathname(ch),iline, &
             irrad_spec(fidx:lidx,1:nx), irrad_prec(fidx:lidx,1:nx), &
             irrad_qflg(fidx:lidx,1:nx), irrad_wavl(fidx:lidx,1:nx), &
             irrad_mflg, nwls(ch), read_irrad, errstat)
        IF (errstat /= 0) THEN 
           message = ADJUSTL(TRIM(modulename))//": failed to read irradiance data"
           go to 123
        endif

        call close_L1_tio (tio_l1obj, errstat)
        IF (errstat /= 0) THEN
           message = ADJUSTL(TRIM(modulename))//": failed to close irradiance data"
           go to 123
        ENDIF
        IF (irrad_wavl(fidx, 1) > irrad_wavl(lidx, 1)) THEN
          idxs(1:nw) = (/ (i, i = lidx, fidx , -1) /)
          irrad_wavl(fidx:lidx, :) = irrad_wavl(idxs(1:nw), :)
          irrad_spec(fidx:lidx, :) = irrad_spec(idxs(1:nw), :)
          irrad_prec(fidx:lidx, :) = irrad_prec(idxs(1:nw), :)
          irrad_qflg(fidx:lidx, :) = irrad_qflg(idxs(1:nw), :)   
        ENDIF
        nwavel = nwavel + nwls(ch)
        spos(ch) = fidx
        epos(ch) = lidx 
        fidx = lidx + 1     
     ENDDO
   ELSE
     thedoy = GranuleJDay
     IF (thedoy == 366) thedoy = 265
     OPEN (UNIT=calunit, FILE= ADJUSTL(TRIM(refdbdir)) // 'solar-distance.dat', &
           STATUS='UNKNOWN', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Cannot open Sun-Earth Distancedatafile!!!'
        pge_error_status = pge_errstat_error
        RETURN
     END IF

     DO i = 1, 12
        READ(calunit, *)
      ENDDO
     DO i = 1, thedoy
        READ(calunit, *) normsc, normsc
     ENDDO
     !   earthsundistance = normsc
     CLOSE(calunit)
     normsc = 1.0 / normsc ** 2  ! solar energy is inversely proportional tosquare distance

     ! Determine backup filename
     bkfname = ADJUSTL(TRIM(refdbdir)) //'TMPO/tmposol_v003_avg_nshi_backup.dat'

      IF( scnwrt ) WRITE(*,*) 'use_backup=(T):'//ADJUSTL(TRIM( bkfname ))
      OPEN (UNIT=calunit, FILE=TRIM(ADJUSTL(bkfname)), STATUS='UNKNOWN',IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Cannot open solar backup file!!!'
        pge_error_status = pge_errstat_error
        RETURN
      END IF     
     ! Determine sun-earth distance correction
     WRITE(*,*) "Please fill up when it is performed with backup = T"
     !IF (scnwrt) WRITE(*,*) 'Backup(T):'//ADJUSTL(TRIM(bkfname))'
   ENDIF  

   !-----------------------------------------------------
   ! Degrade spectral resolution if necessary
   !-----------------------------------------------------
   ! Do not coadd wavelengths with a gap (e.g., filter Mg absorption lines), need
   ! to determine
   ! delta-lamda in UV-1
   ! Note in OMI delta-lamda varies with wavelength (largest for the first two
   ! pixels in each channel)
   dwavmax = ( irrad_wavl(2, 1) - irrad_wavl(1, 1) ) * 1.1

   IF (reduce_resolution) THEN 
     nwavel = 0 ; j = 1
     DO is = 1, nswath
       ch = inschs(is)
       npos = nxtrack
       retswav = retlbnd(ch) ; retewav = retubnd(ch)
       tmpsampr = redsampr
       np = nwls(ch)
       tmpspec(wvl_idx, 1:np, 1:npos) = irrad_wavl(spos(ch):epos(ch), 1:npos)
       tmpspec(spc_idx, 1:np, 1:npos) = irrad_spec(spos(ch):epos(ch), 1:npos)
       tmpspec(sig_idx, 1:np, 1:npos) = irrad_prec(spos(ch):epos(ch), 1:npos)
       tmpqflg(1:np, 1:npos) = irrad_qflg(spos(ch):epos(ch), 1:npos)
       DO ix = 1, npos
         CALL convert_2bytes_to_16bits (nbits, np, tmpqflg(1:np, ix), &
         flgbits(1, 1:np, 0:nbits-1))
         tmpqflg(1:np, ix) = flgbits(1, 1:np, 0) & !Missing
                           + flgbits(1, 1:np, 1) & !Missing
                           + flgbits(1, 1:np, 2)  !Missing
       ENDDO
       CALL reduce_irrad_resolution(tmpspec(:, 1:np, 1:npos), & 
            tmpqflg(1:np,1:npos), np, npos, reduce_slit, redslw(is), & 
         tmpsampr, redlam, retswav, retewav, reduce_lbnd(ch),reduce_ubnd(ch), &
         nwls(ch), pge_error_status)
       IF (pge_error_status == pge_errstat_error) RETURN
       nwavel = nwavel + nwls(ch)
       spos(ch) = j
       j = nwavel + 1
       epos(ch) = nwavel
       irrad_wavl(spos(ch):epos(ch), 1:npos) = real (tmpspec(wvl_idx, 1:nwls(ch),1:npos), kind=r4)
       irrad_spec(spos(ch):epos(ch), 1:npos) = real (tmpspec(spc_idx, 1:nwls(ch),1:npos), kind=r4)
       irrad_prec(spos(ch):epos(ch), 1:npos) = real (tmpspec(sig_idx, 1:nwls(ch),1:npos), kind=r4)
       irrad_qflg(spos(ch):epos(ch), 1:npos) = 0 ! All are good (pre-filtered)
     ENDDO
   ENDIF
   
   ! ----------------------------------------------------
   ! check N of wavelengths for ozone and albedo fitting spectra
   ! ----------------------------------------------------
   IF (nwavel > nwavel_max) THEN
     WRITE(*, '(A)') "Need to increase nwavel_max!!!"
     stop 1
     !pge_error_status = pge_errstat_error
     RETURN
   ENDIF

   ! Determine number of wavelengths to be read for deteriming cloud fraction
   fidx = MAXVAL ( MINLOC ( irrad_wavl(1:nwavel, 1), MASK = &
        (irrad_wavl(1:nwavel, 1) > pos_alb - toms_fwhm * 1.4) ))
   lidx = MAXVAL ( MAXLOC ( irrad_wavl(1:nwavel, 1), MASK = &
        (irrad_wavl(1:nwavel, 1) < pos_alb + toms_fwhm * 1.4) ))
   IF (fidx <1 .OR. lidx > nwavel) THEN
     message= modulename//': Need to change pos_alb/toms_fwhm!!!'
     go to 123
   ENDIF

   nrefl = lidx - fidx + 1
   IF (nrefl > mrefl ) THEN
     message= modulename//': Need to increase mrefl!!!'
     go to 123
   ENDIF

   ! ----------------------------------------------------
   ! subset for valid spectra : 1) winlim, 2) lower/upper spec 3) flgmsks
   ! cadded irradiance spectra
   ! ----------------------------------------------------
   ! determine number of binning for different fitting windows
   nwbin(1:numwin) = nxbin ! TEMPO has same dimension between UV and VIS, 
                           ! So it is more simple than OMI
   DO ix = first_pix, last_pix
     ! Get quality flags bits, coadd flags
     currpix = ix
     flgmsks (:) = 0
     DO is = 1, nswath 
       ch = inschs(is)
       nbin = nxbin
       iix = (ix -1 ) *nbin
       IF (.NOT. reduce_resolution) THEN 
         ! properly align cross track position to be coadded (should be
         ! within one pixel)
         IF (nbin > 2)  THEN 
           CALL prespec_align(nwls(ch), nbin, &
                irrad_wavl(spos(ch):epos(ch), iix+1:iix+nbin), &
                irrad_spec(spos(ch):epos(ch), iix+1:iix+nbin), &
                irrad_prec(spos(ch):epos(ch), iix+1:iix+nbin), &
                irrad_qflg(spos(ch):epos(ch), iix+1:iix+nbin))
         ENDIF
         DO ic = 1, nbin
           CALL convert_2bytes_to_16bits ( nbits, nwls(ch),irrad_qflg(spos(ch):epos(ch), iix + ic ), &
                  flgbits(ic, spos(ch):epos(ch), 0:nbits-1))

           flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) &
                 + flgbits(ic, spos(ch):epos(ch), 0)                &   !Missing
                 + flgbits(ic, spos(ch):epos(ch), 1)                &   !Bad
                 + flgbits(ic, spos(ch):epos(ch), 2)                &   !Processing error
!                 + flgbits(ic, spos(ch):epos(ch), 3)                &   !transient_pixel
!                 + flgbits(ic, spos(ch):epos(ch), 4)                &   !RTS_Pixel_Warning Flag
!                 + flgbits(ic, spos(ch):epos(ch), 5)                &   !Saturation Possibility Flag
!                 + flgbits(ic, spos(ch):epos(ch), 7)                &   !Dark Current Warning Flag
                 + flgbits(ic, spos(ch):epos(ch), 8)                &   ! offset correction error
                 + flgbits(ic, spos(ch):epos(ch), 9)                &   ! smear correction error
                 + flgbits(ic, spos(ch):epos(ch), 10)               &   ! stray light correction error
                 + flgbits(ic, spos(ch):epos(ch), 11)                  ! nonlinear range error

         ENDDO
       ELSE
         ! Already aligned because of using common wavelength scale
         DO ic = 1, nbin 
           flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) + & 
                                        irrad_qflg(spos(ch):epos(ch), iix + ic)
         ENDDO
       ENDIF        
     ENDDO ! loop of nswath

     !1) subset solar spectrum for ozone fitting spectra 
     nsub = 0 ; subspec=0.0 
     DO iw = 1, numwin
       ch = band_selectors(iw)
       nbin = nwbin(iw)
       iix = (ix-1)*nbin

       winpix(iw, 1) = MINVAL ( MINLOC ( irrad_wavl(spos(ch):epos(ch),iix + 1), &
             MASK = irrad_wavl(spos(ch):epos(ch),iix + 1) >= &
             winlim(iw, 1)) ) + spos(ch) - 1
       winpix(iw, 2) = MAXVAL ( MAXLOC ( irrad_wavl(spos(ch):epos(ch),iix + 1),&
             MASK = irrad_wavl(spos(ch):epos(ch),iix + 1) <= &
             winlim(iw, 2)) ) + spos(ch) - 1
       IF (winpix(iw,1) .le. 0 )  winpix(iw, 1) = 1
       IF (winpix(iw,2) .gt. nwavel )  winpix(iw, 2) = nwavel

       tmpo_irrad%winpix(iw, ix, 1:2) =  0
       tmpo_irrad%npix  (iw, ix) = nsub
       !IF (rm_mgline .AND. winlim(iw, 1) < 286.0 .AND. winlim(iw, 2) > 286) THEN 
       DO i = winpix(iw,1), winpix(iw, 2)
         IF (ALL(irrad_spec(i, iix+1:iix+nbin)  > lower_spec)   .AND. &
             ALL(irrad_spec(i, iix+1:iix+nbin)  < upper_spec) .AND. flgmsks(i) == 0 ) THEN
           nsub = nsub + 1
           subspec(1:nbin, wvl_idx, nsub) = irrad_wavl(i,iix+1:iix+nbin)
           subspec(1:nbin, spc_idx, nsub) = irrad_spec(i,iix+1:iix+nbin)
           subspec(1:nbin, sig_idx, nsub) = irrad_prec(i,iix+1:iix+nbin)
           IF (tmpo_irrad%winpix(iw, ix, 1) == 0) tmpo_irrad%winpix(iw,ix, 1) = i
           tmpo_irrad%winpix(iw, ix, 2) = i
           tmpo_irrad%wind(nsub, ix) = i
         END IF 
         !  print * , i, nsub,irrad_wavl(i, iix+1), flgmsks(i)  ,irrad_spec(i, iix+1), flgmsks(i)
       ENDDO
       tmpo_irrad%npix(iw, ix) = nsub - tmpo_irrad%npix(iw, ix)
     ENDDO
     tmpo_irrad%nwav(ix) = nsub

     IF (nsub == 0) tmpo_irrad%errstat(ix) = pge_errstat_error
     IF (tmpo_irrad%errstat(ix) == pge_errstat_error)    CYCLE

     !------------------------------------
     ! Perform coadding when necessary
     !-------------------------------------
     fidx = 1
     DO iw = 1, numwin
       nbin = nwbin(iw)
       lidx = fidx + tmpo_irrad%npix(iw, ix) - 1
       IF (nbin > 1) THEN
         CALL solwavcal_coadd(wcal_bef_coadd, tmpo_irrad%npix(iw, ix), nbin, &
              subspec(1:nbin, :, fidx:lidx), wshis(1:nbin), wsqus(1:nbin), problems)
         IF (problems) THEN
           WRITE(*, '(A)') 'No solar wavelength calibration beforecoadding!!!'
           tmpo_irrad%errstat(ix) = pge_errstat_warning
         ENDIF
       ENDIF
       fidx = lidx + 1
     ENDDO
     !----------------------------------------
     ! Subset solar spectrum for Ring effect
     !----------------------------------------
     IF (.NOT. reduce_resolution .OR. (reduce_resolution .AND. .NOT. use_redfixwav)) THEN
        ! ozone fitting spectrum of first window (coadded) is transfered to ring spectrum 
       subring = 0.0
       noff1 = n_ring_off
       nring = tmpo_irrad%npix(1, ix) + noff1
       subring(1:spc_idx, noff1+1 : nring) = subspec(1, 1:spc_idx,1:tmpo_irrad%npix(1, ix))
       tmpo_ring%ndiv(ix) = 0
       ! add extra spectra before first window (uncoadded)
       ! if unavailable, needed to ammened with solar reference spectrum
       noff1 = noff1 + 1 
       nbin = nwbin(1)
       iix = (ix -1)*nbin
       nbad  = 0
       DO i = winpix(1, 1) - 1, 1, -1
         !print * , '(a)', irrad_wavl(i, iix+1:iix+nbin)
         IF (ALL(irrad_spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
             ALL(irrad_spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
           noff1 = noff1 - 1
           subring(wvl_idx, noff1) = SUM(irrad_wavl(i, iix+1:iix+nbin))/ nbin
           subring(spc_idx, noff1) = SUM(irrad_spec(i, iix+1:iix+nbin))/ nbin
           IF (noff1 == 1) EXIT
         ELSE !jbak
           nbad = nbad + 1
           IF (nbad == 5) EXIT
         ENDIF
       ENDDO
       ! add the ring spectrum after the first ozone fitting window
       DO iw = 2, numwin 
         ! the btw the ozone fitting window
         ch = band_selectors(iw)
         nbin = nwbin(iw -1)
         iix = (ix - 1)*nbin
         IF (ch == band_selectors(iw -1)) THEN ! the current fitting window is belong to the same channel
           DO i = winpix(iw-1, 2)+1, winpix(iw, 1)-1
                !print * , '(b)', irrad_wavl(i, iix+1:iix+nbin)
             IF (ALL(irrad_spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
                 ALL(irrad_spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
               nring = nring + 1
               subring(wvl_idx, nring) = SUM(irrad_wavl(i,iix+1:iix+nbin)) / nbin
               subring(spc_idx, nring) = SUM(irrad_spec(i,iix+1:iix+nbin)) / nbin
             ENDIF
           ENDDO
         ELSE ! the current fittign window is belong to the different channel
           wcenter = (winlim(iw-1, 2) + winlim(iw, 1)) / 2.0
           idx = MAXVAL ( MAXLOC ( irrad_wavl(spos(ch-1):epos(ch-1),iix+1), &
                   MASK = irrad_wavl(spos(ch-1):epos(ch-1), iix+1) < wcenter ) )+ spos(ch-1) - 1
           j = 0        
           DO i = winpix(iw-1, 2)+1, idx
               ! print * ,'(c)', irrad_wavl(i, iix+1:iix+nbin)
             IF (ALL(irrad_spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
                   ALL(irrad_spec(i, iix+1:iix+nbin) < upper_spec) .AND.flgmsks(i) == 0 .AND. &
                   ALL(irrad_wavl(i, iix+1:iix+nbin) > subring(wvl_idx,nring)) ) THEN
               nring = nring + 1
               subring(wvl_idx, nring) = SUM(irrad_wavl(i,iix+1:iix+nbin)) / nbin
               subring(spc_idx, nring) = SUM(irrad_spec(i,iix+1:iix+nbin)) / nbin            
               j = j + 1
               IF (j == n_ring_off ) exit
                 !print * , nring, subring(1, nring), '2'
             ENDIF
           ENDDO
           tmpo_ring%ndiv(ix) = nring         
           nbin = nwbin(iw)
           idx = MAXVAL ( MINLOC ( irrad_wavl(spos(ch):epos(ch), iix+1),   &
                   MASK = irrad_wavl(spos(ch):epos(ch), iix+1) > wcenter ) ) + spos(ch) - 1        
           j = 0
           DO i = idx, winpix(iw, 1) - 1 
               ! print * ,'(d)', irrad_wavl(i, iix+1:iix+nbin)
             IF (ALL(irrad_spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
                 ALL(irrad_spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0 .AND. &
                 ALL(irrad_wavl(i, iix+1:iix+nbin) > subring(wvl_idx, nring)) ) THEN
               nring = nring + 1
               subring(wvl_idx, nring) = SUM(irrad_wavl(i,iix+1:iix+nbin)) / nbin
               subring(spc_idx, nring) = SUM(irrad_spec(i,iix+1:iix+nbin)) / nbin
               j = j + 1
               if (j == n_ring_off) exit
             ENDIF
           ENDDO
         ENDIF
                   
         idx = SUM(tmpo_irrad%npix(1:iw-1, ix))
         ! print *,'(e)', subspec(1, wvl_idx, idx+1), subspec(1, wvl_idx, idx+tmpo_irrad%npix(iw,ix))
         subring(wvl_idx:spc_idx, nring+1:nring+tmpo_irrad%npix(iw, ix)) = &
                subspec(1, wvl_idx:spc_idx,idx+1:idx+tmpo_irrad%npix(iw,ix))
         nring = nring + tmpo_irrad%npix(iw, ix)
       ENDDO
       ! Add extra spectra after fitting window
       noff2 = nring
       nring = nring + n_ring_off
       nbad = 0
       DO i =  winpix(numwin, 2) + 1, nwavel
         IF (ALL(irrad_spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
             ALL(irrad_spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
           noff2 = noff2 + 1
           subring(wvl_idx, noff2) = SUM(irrad_wavl(i, iix+1:iix+nbin))/ nbin
           subring(spc_idx, noff2) = SUM(irrad_spec(i, iix+1:iix+nbin))/ nbin 
         ENDIF
         ! print* , i, noff2,subring(1,noff2)
         IF (noff2 == nring) EXIT
       ENDDO
       tmpo_ring%nsol(ix) = nring
       tmpo_ring%winpix(ix,1) = noff1
       tmpo_ring%winpix(ix,2) = noff2
     ELSE
       ! Need to convole the saved solar spectra with additional slit width
       nring = tmpo_ring%nsol(ix)
       subring(1,1:nring) = tmpo_solspec_ring(1, 1:nring, ix)
       subring(2,1:nring) = tmpo_solspec_ring(2, 1:nring, ix)

       IF (which_slit == 0) THEN
          CALL gauss_uneven(subring(1, 1:nring),subring(2, 1:nring), nring, &
               nswath, redslw(inschs(1:nswath)), retlbnd(inschs(1:nswath)),retubnd(inschs(1:nswath)))
       ELSE IF (which_slit == 3) THEN
          CALL triangle_uneven(subring(1, 1:nring), subring(2, 1:nring), nring,&
               nswath, redslw(inschs(1:nswath)), retlbnd(inschs(1:nswath)),retubnd(inschs(1:nswath)))
       ELSE
         WRITE(*, *) 'This type of slit convolution is not implemented!!!'
         stop 1
            !pge_error_status = pge_errstat_error
       ENDIF  
       WRITE(*, *) 'This type is  not implemented!!!'
       stop 1
     ENDIF
     !===============================================================================
     ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 15 pixels
     !=============================================================================
     irefl = 0; tmpo_refl%winpix(ix, 1:2) = 0
     nbin = nxbin
     idx = MAXVAL ( MINLOC ( irrad_wavl(1:nwavel, iix+1), MASK = &
            (irrad_wavl(1:nwavel, iix+1) > pos_alb - toms_fwhm * 1.4) ))
     DO i  = idx, nwavel
       IF (ALL(irrad_spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
           ALL(irrad_spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
         irefl = irefl + 1
         tmpo_refl%solwavl( irefl, ix) = SUM(irrad_wavl(i,iix+1:iix+nbin)) / nbin
         tmpo_refl%solspec( irefl, ix) = SUM(irrad_spec(i,iix+1:iix+nbin)) / nbin
         IF (tmpo_refl%winpix(ix, 1) == 0) tmpo_refl%winpix(ix, 1) = i
         tmpo_refl%winpix(ix, 2) = i
       ENDIF
       IF (irefl == nrefl) EXIT
     ENDDO
     IF (irefl /= nrefl) THEN
       WRITE(www_lun, '(i5,A,2i4)') ix,'Could not get enough irrad for cf!!!',irefl, nrefl
       tmpo_irrad%errstat(ix) = pge_errstat_error
     ENDIF

     !IF(scnwrt) WRITE(*,'(A, I4, A,I3, A, 2I4)') & 
     !    '(irrad) ix=',ix,' irefl = ',irefl, ' Nw=',tmpo_irrad%npix(1:numwin, ix)

     DO i = 1, numwin
         IF (tmpo_irrad%npix(i, ix) < 10) THEN
           WRITE(www_lun, '(A,f8.3,A3,f8.3)') ' Not enough irrad (>=10) in win:', winlim(i,1), ' - ', winlim(i,2)
           tmpo_irrad%errstat(ix) = pge_errstat_error 
         ENDIF
     ENDDO

     !fidx = tmpo_irrad%npix(1,ix)+1
     tmpo_irrad%norm(ix) = SUM ( subspec(1, spc_idx, 1:nsub) ) / nsub
     !tmpo_irrad%norm(ix) = SUM ( subspec(1, spc_idx, fidx:nsub) ) / (nsub-fidx +1)
     IF ( tmpo_irrad%norm(ix) <= 0.0 ) THEN
        tmpo_irrad%errstat(ix) = pge_errstat_error; CYCLE
     ENDIF

     tmpo_irrad%wavl(1:nsub, ix)  = subspec(1, wvl_idx, 1:nsub)
     tmpo_irrad%spec(1:nsub, ix)  = subspec(1, spc_idx, 1:nsub) /tmpo_irrad%norm(ix)
     tmpo_irrad%prec(1:nsub, ix)  = subspec(1, sig_idx, 1:nsub) /tmpo_irrad%norm(ix)
     tmpo_ring%wavl(1:nring, ix)  = subring(1, 1:nring)
     tmpo_ring%spec(1:nring, ix)  = subring(2, 1:nring) / tmpo_irrad%norm(ix)
     dwavmax = (tmpo_irrad%wavl(2,1) - tmpo_irrad%wavl(1,1))*1.1
   ENDDO ! end cross track loop
   !--------------------------------------------------------------------------
   ! Ending  with deallocating local variables
   !--------------------------------------------------------------------------
   deallocate (irrad_qflg , irrad_prec, irrad_spec, irrad_wavl)
   deallocate (idxs, flgmsks, flgbits)
   deallocate (tmpspec, tmpqflg)
   deallocate (subspec, subring)
  RETURN
  123 continue
   pge_error_status = pge_errstat_error
    WRITE(*,*) message
    call tell_error (tell_io_error, message, pge_error_status)
  RETURN
  END SUBROUTINE tmpo_read_irradiance

  SUBROUTINE tmpo_read_radiance_lines (iline, first_pix, last_pix, sline, eline,  pge_error_status)
   USE OMSAO_parameters_module, ONLY:max_fit_pts, maxwin, mswath
   USE OMSAO_indices_module, ONLY:sig_idx, spc_idx, wvl_idx
   USE OMSAO_variables_module, ONLY: nswath, nxtrack, inschs,&
       l1b_rad_filename, nxbin, nybin, inschs, numwin, reduce_resolution, &
       upper_spec, lower_spec, &
       wcal_bef_coadd, ybin_decerr, szamax
   USE OMSAO_precision_module
   USE OMSAO_errstat_module
   USE ozprof_data_module, ONLY: nrefl
   IMPLICIT NONE
   ! ---------------
   ! Input variables
   ! ---------------
   INTEGER,  INTENT (IN) :: iline, first_pix, last_pix, sline, eline

   ! ----------------
   ! Output variables
   ! ----------------
   INTEGER, INTENT (OUT) :: pge_error_status

   !--------------------
   ! local variables
   !--------------------
   TYPE (tiof_file_type) :: tio_l1obj
   !----------------------------------------------------
   !final spectrum input after subset and coadding
   !----------------------------------------------------
   LOGICAL :: read_irrad=.false.
   LOGICAL :: problems
   INTEGER (KIND=i2) :: tmp_mflg
   INTEGER :: i, j, fidx, lidx, is, ch,  errstat, iloop, ix, iix, nbin, iw, &
              ic, nsub, blockline, nx, nt, nw, nwavel, irefl, ii, nl
   ! variables used to read TEMPO_data 
   INTEGER, DIMENSION (mswath) :: nwls, epos, spos
   INTEGER (KIND=i4), DIMENSION (:), ALLOCATABLE :: idxs
   INTEGER (kind=2), ALLOCATABLE, DIMENSION (:,:,:) :: rad_qflg
   REAL (kind=i4), ALLOCATABLE, DIMENSION(:,:,:) ::  rad_spec,rad_prec,rad_wavl
   !REAL (kind=i4), DIMENSION (nwavel_ccd,nxtrack_max) ::  ccd_spec,ccd_prec,ccd_wavl
   !INTEGER (kind=2), dimension(nwavel_ccd,nxtrack_max)::  ccd_qflg
   REAL (kind=i4), DIMENSION (:,:), ALLOCATABLE :: ccd_spec,ccd_prec,ccd_wavl
   INTEGER (kind=2), dimension(:,:),ALLOCATABLE :: ccd_qflg
   ! variables used for reduced resolution
   ! Subset variables
   INTEGER, PARAMETER :: nbits = 16
   INTEGER (KIND=i4), DIMENSION (maxwin)     :: nwbin
   INTEGER (KIND=i2), DIMENSION (:), ALLOCATABLE :: flgmsks
   INTEGER (kind=i2), DIMENSION (:,:,:), ALLOCATABLE :: flgbits
   LOGICAL, DIMENSION (maxwin, nxtrack_max)  :: wavcals
   REAL (KIND=dp), DIMENSION (maxwin,nxbin) :: wshis, wsqus
   ! subset
   REAL (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: subspec
   CHARACTER (LEN=100) :: message
   ! ------------------------------
   ! Name of this module/subroutine
   ! ------------------------------
   CHARACTER (LEN=*), PARAMETER :: modulename = 'tmpo_read_radiance_lines'

   !--------------------------------------------------------------------------
   ! Starting with allocating local variables
   !--------------------------------------------------------------------------
   nl = (eline - sline + 1)/nybin
   allocate (rad_qflg (nwavel_max, nxtrack_max, 0:nl-1))
   allocate (rad_prec (nwavel_max, nxtrack_max, 0:nl-1))
   allocate (rad_wavl (nwavel_max, nxtrack_max, 0:nl-1))
   allocate (rad_spec (nwavel_max, nxtrack_max, 0:nl-1))
   allocate (flgmsks(nwavel_max),flgbits (nxbin, nwavel_max, 0:nbits-1))
   allocate (idxs(nwavel_ccd))
   allocate (ccd_prec(nwavel_ccd, nxtrack_max))
   allocate (ccd_spec(nwavel_ccd, nxtrack_max))
   allocate (ccd_wavl(nwavel_ccd, nxtrack_max))
   allocate (ccd_qflg(nwavel_ccd, nxtrack_max))
   allocate (subspec(nxbin, sig_idx, max_fit_pts))
   !--------------------------------------------------------------------------
   ! Initializing
   !--------------------------------------------------------------------------
   errstat = 0
   pge_error_status = pge_errstat_ok
   tmpo_rad%errstat(0:nl-1)          = pge_errstat_ok
   tmpo_rad%pix_errstat(1:nxtrack, 0:nl-1) = pge_errstat_ok
   tmpo_rad%npix (:,:,:) = 0
   tmpo_rad%spec (:,:,:) = 0.0
   tmpo_rad%prec (:,:,:) = 0.0
   tmpo_rad%qflg (:,:,:) = 0
   tmpo_rad%wavl (:,:,:) = 0.0

   rad_spec = 0.0
   rad_wavl = 0.0
   rad_prec = 0.0  
   rad_qflg = 0

   fidx = 1
   nwavel = 0
   DO is = 1, nswath
     ch = inschs(is)
     errstat=0 ! FIXME remove when errstat properly defined
     call read_L1_dims_tio (l1b_rad_filename, rad_swathname(ch),&
            nt, nx,nw, errstat = errstat)
     call open_L1_tio (l1b_rad_filename, tio_l1obj, errstat)
     IF (errstat /= 0) then
       message = ADJUSTL(TRIM(modulename))//": failed to open radiance file"
       go to 123
     ENDIF
     lidx = fidx + nw -1
     DO iloop = 0,  nl-1
       DO i = 0 , nybin -1
         blockline = (sline-1) + iloop*nybin + i ! this suroutine read data from 0 to nl-1
          IF ( mod(nl, 50) .eq. 0) WRITE(*,*) '==>readl1b:', is, blockline
         !print * , blockline, nx, nw, nwavel_ccd, rad_swathname(is)
         call read_L1_rad_line_tio (tio_l1obj, rad_swathname(ch), &
              blockline, &
              radiance           = ccd_spec(1:nw, 1:nx), &
              rad_precision      = ccd_prec(1:nw, 1:nx), &
              pixel_quality_flag = ccd_qflg(1:nw, 1:nx), &
              wavelengths        = ccd_wavl(1:nw, 1:nx), &
              meas_qual_flag     = tmp_mflg, &
              num_wavelengths    = nwls(ch),&
              read_irrad         = read_irrad, &
              errstat = errstat)
         IF (errstat /= 0) then
           message=ADJUSTL(TRIM(modulename))//": failed to read from radiancefile"
           go to  123
         ENDIF
         rad_spec(fidx:lidx,1:nx, iloop) = rad_spec(fidx:lidx,1:nx, iloop) + ccd_spec(1:nw,1:nx)
         rad_wavl(fidx:lidx,1:nx, iloop) = rad_wavl(fidx:lidx,1:nx, iloop) + ccd_wavl(1:nw,1:nx)
         rad_prec(fidx:lidx,1:nx, iloop) = rad_prec(fidx:lidx,1:nx, iloop) + ccd_prec(1:nw,1:nx)
         DO ix = 1, nx
            CALL coadd_2bytes_qflgs(nbits, nwls(ch), rad_qflg(fidx:lidx,ix, iloop), ccd_qflg(1:nw, ix))
         ENDDO
       ENDDO ! end binloop
   
       rad_spec(fidx:lidx,1:nx, iloop) = rad_spec(fidx:lidx,1:nx, iloop) / nybin
       rad_wavl(fidx:lidx,1:nx, iloop) = rad_wavl(fidx:lidx,1:nx, iloop) / nybin
       IF (ybin_decerr) THEN
         rad_prec(fidx:lidx,1:nx, iloop) = real(rad_prec(fidx:lidx,1:nx,iloop) / nybin / &
                                   SQRT(1.0D0 * nybin) , kind=r4)
       ELSE
         rad_prec(fidx:lidx,1:nx, iloop) = rad_prec(fidx:lidx,1:nx, iloop) / nybin
       ENDIF
       IF ( iloop == 0) THEN 
        ! Sort data in wavelength increasing order
          IF (ccd_wavl(1, 1) > ccd_wavl(nw, 1))THEN
          idxs(1:nw) = (/ (j, j = lidx, fidx , -1) /)
          ELSE
          idxs(1:nw) = (/ (j, j = fidx, lidx ) /)
          ENDIF
        ENDIF
        rad_wavl(fidx:lidx,1:nx, iloop) = rad_wavl(idxs(1:nw),1:nx, iloop)
        rad_spec(fidx:lidx,1:nx, iloop) = rad_spec(idxs(1:nw),1:nx, iloop)
        rad_prec(fidx:lidx,1:nx, iloop) = rad_prec(idxs(1:nw),1:nx, iloop)
        rad_qflg(fidx:lidx,1:nx, iloop) = rad_qflg(idxs(1:nw),1:nx, iloop)
      END DO     ! end iloop

      call close_L1_tio (tio_l1obj, errstat)
      IF (errstat /= 0) then
        message=ADJUSTL(TRIM(modulename))//": failed to close radiance file"
        go to 123
      ENDIF
      nwavel  = nwavel + nwls(ch)
      spos(ch) = fidx
      epos(ch) = lidx
      fidx = lidx + 1 
    ENDDO ! end swath loop


    !deallocate (ccd_spec,ccd_prec, ccd_wavl, ccd_qflg)

    IF (reduce_resolution) THEN
      message=ADJUSTL(TRIM(modulename))//": reduce_resolution = .true."
      goto 123
    ENDIF

    IF (nwavel > nwavel_max) THEN
      message=ADJUSTL(TRIM(modulename))//": Need to increase nwavel_max!!!"
      go to 123
      RETURN
    ENDIF
    nwbin(1:numwin) = nxbin

    DO iloop = 0, nl-1
      ! Subset and coadd radiance spectrum
      IF (ALL(tmpo_rad%pix_errstat(first_pix:last_pix, iloop) == pge_errstat_error)) CYCLE
      
      DO ix = first_pix, last_pix
        IF (tmpo_geo%sza (ix, iloop+iline) > szamax .OR. tmpo_geo%sza (ix,iloop+iline) < 0 ) THEN
           tmpo_rad%pix_errstat(ix, iloop) = pge_errstat_error
          CYCLE
        ENDIF
        flgmsks = 0
        DO is = 1, nswath 
          ch = inschs(is)
          nbin = nxbin
          iix = (ix -1)*nbin
          IF (.NOT. reduce_resolution) THEN 
            IF (nbin > 2) CALL prespec_align(nwls(ch), nbin,rad_wavl(spos(ch):epos(ch),&
                 iix+1:iix+nbin, iloop), rad_spec(spos(ch):epos(ch),iix+1:iix+nbin, iloop), &
                 rad_prec(spos(ch):epos(ch), iix+1:iix+nbin, iloop), &
                 rad_qflg(spos(ch):epos(ch), iix+1:iix+nbin, iloop))
            DO ic = 1, nbin
               CALL convert_2bytes_to_16bits ( nbits, nwls(ch),rad_qflg(spos(ch):epos(ch), &
                   iix + ic, iloop), flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
               flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) &
                   + flgbits(ic, spos(ch):epos(ch), 0)                &   !Missing
                   + flgbits(ic, spos(ch):epos(ch), 1)                &   !Bad
                   + flgbits(ic, spos(ch):epos(ch), 2)                &   !Processing error
!                  + flgbits(ic, spos(ch):epos(ch), 3)                &   !transient pixel 
!                  + flgbits(ic, spos(ch):epos(ch), 4)                &   !RTS_Pixel_Warning Flag
!                  + flgbits(ic, spos(ch):epos(ch), 5)                &   !Saturation Possibility Flag
!                  + flgbits(ic, spos(ch):epos(ch), 7)                &   !Dark Current Warning Flag
                   + flgbits(ic, spos(ch):epos(ch), 8)                &   ! offset correction error
                   + flgbits(ic, spos(ch):epos(ch), 9)                &   ! smear correction error
                   + flgbits(ic, spos(ch):epos(ch), 10)               &   ! stray light correction error
                   + flgbits(ic, spos(ch):epos(ch), 11)                  ! nonlinear range error
        !print * , flgbits(ic, tmpo_refl%winpix(ix,1),0:11) 
            ENDDO
          ELSE
            ! Already aligned because of using common wavelength scale
            DO ic = 1, nbin
             flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) + &
             rad_qflg(spos(ch):epos(ch), iix + ic, iloop)
            ENDDO
          ENDIF
        ENDDO ! end xtrack
     
        ! Subset valid data
        nsub = 0
        subspec = 0.0
        fidx = 1

        DO iw = 1, numwin
          tmpo_rad%npix(iw, ix, iloop) = nsub
          nbin = nwbin(iw)
          iix = (ix - 1) * nbin
          lidx = fidx + tmpo_irrad%npix(iw, ix) - 1

          DO ii = fidx, lidx
            i = tmpo_irrad%wind(ii, ix)
            IF (ALL(rad_spec(i, iix+1:iix+nbin, iloop) > lower_spec) .AND. &
                 ALL(rad_spec(i, iix+1:iix+nbin, iloop) < upper_spec) .AND. flgmsks(i) == 0 ) THEN
              nsub = nsub + 1
              subspec(1:nbin, wvl_idx, nsub) = rad_wavl(i,iix+1:iix+nbin, iloop)
              subspec(1:nbin, spc_idx, nsub) = rad_spec(i,iix+1:iix+nbin, iloop)
              subspec(1:nbin, sig_idx, nsub) = rad_prec(i,iix+1:iix+nbin, iloop)
              tmpo_rad%wind(nsub, ix, iloop) = int(ii, kind=i2)

            ENDIF
            !WRITE(*,'(4i4,f8.2,2e17.5)')  ii, i,nsub,flgmsks(i), &
            ! rad_wavl(i, iix+1, iloop),rad_spec(i,iix+1:iix+nbin, iloop), rad_prec(i, iix+1, iloop)
          ENDDO
          fidx = lidx + 1
          tmpo_rad%npix(iw, ix, iloop) = nsub - tmpo_rad%npix(iw, ix, iloop)
          ! processing this pixel
          IF (tmpo_rad%npix(iw, ix, iloop) <= tmpo_irrad%npix(iw, ix) * 0.80 ) THEN
            WRITE(*, '(3I5, A, 2I5, F8.2)') ix, iloop,iw, ': Too fewer #of rad=',  &
            tmpo_rad%npix(iw, ix, iloop),tmpo_irrad%npix(iw, ix), tmpo_geo%sza(ix, iloop+iline)
            tmpo_rad%pix_errstat(ix, iloop) = pge_errstat_error
          ENDIF
        ENDDO ! end numwin
        tmpo_rad%nwav(ix, iloop) = nsub
        IF (tmpo_rad%pix_errstat(ix, iloop) == pge_errstat_error) CYCLE   ! This pixel will not be processed.

        !-------------------------------------------
        ! Perform coadding
        !-------------------------------------------
        fidx = 1
        DO iw = 1, numwin
          nbin = nwbin(iw)
          lidx = fidx + tmpo_rad%npix(iw, ix, iloop) - 1
          IF (nbin > 1) THEN
            CALL radwavcal_coadd(wcal_bef_coadd, wavcals(iw, ix), & !iw, ix, &
              tmpo_rad%npix(iw, ix, iloop), nbin, &
              subspec(1:nbin, :, fidx:lidx), wshis(iw, 1:nbin), &
              wsqus(iw,  1:nbin), problems)
              wavcals(iw, ix) = .FALSE.
            IF (problems) THEN
              WRITE(*, '(A)') 'No radiance wavelength calibration before coadding!!!'
              pge_error_status = pge_errstat_warning
            ENDIF
          ENDIF
          fidx = lidx + 1
        ENDDO
        !--------------------------------------------------------------------
        ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 20 pixels
        irefl = 0
        nbin = nxbin
        iix = (ix - 1) * nbin
        DO i  = tmpo_refl%winpix(ix, 1), tmpo_refl%winpix(ix, 2) + 10
          IF ( ALL(rad_spec(i, iix+1:iix+nbin, iloop) > lower_spec) .AND. & 
            ALL(rad_spec(i, iix+1:iix+nbin, iloop) < upper_spec) .AND. flgmsks(i) == 0) THEN
            irefl = irefl + 1
            tmpo_refl%radwavl(irefl, ix, iloop) = SUM(rad_wavl(i,iix+1:iix+nbin, iloop)) / nbin
            tmpo_refl%radspec(irefl, ix, iloop) = SUM(rad_spec(i,iix+1:iix+nbin, iloop)) / nbin
          ENDIF
          !print * , i, irefl, flgmsks(i) , rad_wavl(i, iix, iloop) 
          IF (irefl == nrefl) EXIT
        ENDDO
        IF (irefl /= nrefl) THEN
          WRITE(*, '( 2I5, A,i2.2)') ix, iloop,  'Number of rad/sol (cf) do not match: ',irefl,nrefl 
              !print * , tmpo_refl%solwavl(1:nrefl, ix), tmpo_refl%winpix(ix, :)
             !print * ,rad_spec(tmpo_refl%winpix(ix,1):tmpo_refl%winpix(ix, 2), ix, iloop)
             !print * , flgmsks(tmpo_refl%winpix(ix,1):tmpo_refl%winpix(ix, 2))
             tmpo_rad%pix_errstat(ix, iloop) = pge_errstat_error
        ENDIF

        !if (nsub > 0) then
          tmpo_rad%norm(ix, iloop) = SUM ( subspec(1, spc_idx, 1:nsub) ) / nsub
        !else
        !  tmpo_rad%norm(ix, iloop) = 0.0
        !endif

        IF ( tmpo_rad%norm(ix, iloop) <= 0.0 ) THEN
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF

        tmpo_rad%wavl(1:nsub, ix, iloop) = real(subspec(1, wvl_idx, 1:nsub) , kind=r4)
        tmpo_rad%spec(1:nsub, ix, iloop) = real(subspec(1, spc_idx, 1:nsub) / & 
                                           tmpo_rad%norm(ix, iloop) , kind=r4)
        tmpo_rad%prec(1:nsub, ix, iloop) =  real(subspec(1, sig_idx, 1:nsub) / &
                                           tmpo_rad%norm(ix, iloop) , kind=r4)
      ENDDO  ! end xtrack 
    ENDDO ! end iloop
    !--------------------------------------------------------------------------
    ! Finishing with deallocating local variables
    !--------------------------------------------------------------------------
    deallocate (rad_qflg, rad_prec, rad_wavl, rad_spec)
    deallocate (flgmsks,flgbits, idxs, subspec)
    deallocate (ccd_spec, ccd_prec, ccd_wavl, ccd_qflg)

    RETURN
  123 continue
     pge_error_status = pge_errstat_error
     WRITE(*,*) message
     call tell_error (tell_io_error, message, pge_error_status)
  END SUBROUTINE tmpo_read_radiance_lines

end module tmpo_read_l1b_data
