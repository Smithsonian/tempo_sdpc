!! We don't consider Zoom-in mode
module omi_read_l1b_data
  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: spc_idx
  USE OMSAO_parameters_module, ONLY: max_ring_pts
  use m_utilities, ONLY: get_doy
  use m_convert_coadd, only: coadd_2bytes_qflgs, convert_2bytes_to_16bits, &
       prespec_align, solwavcal_coadd
  USE m_ezspline_interpolation, only: interpolation
  USE m_fitting_util, ONLY: reduce_rad_resolution, reduce_irrad_resolution
  USE OMSAO_omidata_module,    ONLY:  omi_radiance_swathname,&
    nfxtrack,nwavel_max, nxtrack_max, ntimes_max,nlines_max, & 
    omi_geo, omi_irrad, omi_rad, omi_ring, omi_refl

  CHARACTER(len=9),PRIVATE   :: omiraddate
  INTEGER (KIND=i2), PRIVATE :: omi_mflg
  INTEGER (KIND=i1), DIMENSION (:), PRIVATE, POINTER :: omi_saa_flag
  REAL (KIND=r4), DIMENSION (:,:,:),PRIVATE, POINTER :: omi_solspec_ring

  public omi_read_radiance_paras, find_scan_line_range, &
         omi_read_irradiance_data, omi_read_radiance_lines, &
         omi_set_parameters, replace_solar_irradiance
  private 

contains

  SUBROUTINE omi_set_parameters (pge_error_status )
    USE OMSAO_parameters_module, ONLY: vb_lev_omidebug!, maxchlen
    USE OMSAO_variables_module, ONLY:ntimes, nxtrack, inschs, band_selectors, ncoadd,&
        l1b_rad_filename, GranuleYear,  GranuleMonth,GranuleDay, GranuleJDay, nswath,verb_thresh_lev
    USE OMSAO_omidata_module , ONLY: zoom_mode, omi_radiance_swathname,nfxtrack
    USE OMSAO_errstat_module
    USE hdfeos4_parameters
    USE L1B_Reader_class
    IMPLICIT NONE

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER, INTENT (OUT) :: pge_error_status

    ! ---------------
    ! Local variables
    ! ---------------
    TYPE (l1b_block_type) :: omi_data_block
    INTEGER               :: i
    INTEGER (KIND=i4)     :: errstat, iline
    REAL (KIND=r4), DIMENSION (1:nxtrack_max) :: tmp_vza  
    ! Exteranl functions
    INTEGER               :: estat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_paras'

    ! --------------------------
    ! Initialize OUTPUT variable
    ! --------------------------
    pge_error_status = pge_errstat_ok
    errstat          = omi_s_success
    ntimes = 0
    nxtrack = 0
    zoom_mode = .FALSE.

    allocate (omi_saa_flag(ntimes_max))
    allocate (omi_solspec_ring(spc_idx, max_ring_pts, nxtrack_max))

    ! ----------------------------------------------------------------
    ! Name of solar and earthshine swaths (normally obtained from PCF)
    ! ----------------------------------------------------------------
    IF (nswath == 2) THEN
      inschs(1) = 1
      inschs(2) = 2
    ELSE IF (nswath == 1) THEN
      ! Ozone profile retrieval with channel 1 only (impossible due to always
      ! using uv2 for fc)
      IF (band_selectors(1) == 1) THEN
        inschs(1) = 1
      ELSE
        ! Ozone profile retrieval with channel 2 only (total ozone retrieval)
        inschs(1) = 2
      ENDIF
    ELSE
      WRITE(*, '(A)') 'Need and only need UV swathes for ozone profileretrieval!!!'
      pge_error_status = pge_errstat_error
    ENDIF

    i = INDEX(l1b_rad_filename, '-o') -14
    omiraddate = l1b_rad_filename(i : i + 8)
    READ(omiraddate, '(I4,1X,2I2)') GranuleYear,  GranuleMonth,GranuleDay
    CALL GET_DOY(GranuleYear,  GranuleMonth,GranuleDay, GranuleJDay)

    ! Determine if UV2 data are observed by zoom mode
    IF (nswath == 2 .OR. inschs(1) == 2)  THEN

      ! -----------------------------------------------------------------------
      ! Open data block called 'omi_data_block' with default size of 100 lines
      ! -----------------------------------------------------------------------
      errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(nswath))
      IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_open_l1b_file, 'L1Br_open failed.', modulename, 0)
        STOP 1
      END IF

      ! -------------------------------
      ! Get dimensions of current Swath
      ! -------------------------------
      errstat = L1Br_getSWdims ( omi_data_block, NumTimes_k=ntimes, nXtrack_k=nfxtrack)
      IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getSWdims failed.', modulename, 0)
        STOP 1
      END IF

      iline = ntimes / 2
      errstat = L1Br_getGEOline ( omi_data_block, iline,    &
           ViewingZenithAngle_k    = tmp_vza (1:nfxtrack))

      IF (tmp_vza(nfxtrack) == tmp_vza(nfxtrack-1) .AND. tmp_vza(nfxtrack-1) &
           == tmp_vza(nfxtrack-2) ) zoom_mode = .TRUE.

      ! --------------------------
      ! Close data block structure
      ! --------------------------
      errstat = L1Br_CLOSE ( omi_data_block )
      IF ( errstat /= omi_s_success .AND. verb_thresh_lev >= vb_lev_omidebug ) THEN
        estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, 'L1Br_CLOSE failed.', modulename, 0)
        STOP 1    
      END IF
    ENDIF

    IF (nswath ==2 .OR. inschs(1) == 1 ) THEN
      ! -----------------------------------------------------------------------
      ! Open data block called 'omi_data_block' with default size of 100 lines
      ! -----------------------------------------------------------------------
      errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(1))
      IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_open_l1b_file, 'L1Br_open failed.', modulename, 0)
        STOP 1
      END IF

      ! -------------------------------
      ! Get dimensions of current Swath
      ! -------------------------------
      errstat = L1Br_getSWdims ( omi_data_block, NumTimes_k=ntimes, nXtrack_k=nfxtrack)
      IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getSWdims failed.', modulename, 0)
        STOP 1
      END IF

      ! --------------------------
      ! Close data block structure
      ! --------------------------
      errstat = L1Br_CLOSE ( omi_data_block )
      IF ( errstat /= omi_s_success .AND. verb_thresh_lev >= vb_lev_omidebug ) THEN
        estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, 'L1Br_CLOSE failed.', modulename, 0)
        STOP 1    
      END IF
    ENDIF

    ! Note that nfxtrack means the number of pixels for UV-1 if it is selected
    nxtrack = nfxtrack * ncoadd

    IF (ntimes > ntimes_max) THEN
      pge_error_status = pge_errstat_error
      WRITE(www_lun, '(A,I5)') 'Need to increase ntimes_max >= ', ntimes 
    ENDIF

    IF (nxtrack > nxtrack_max) THEN
      pge_error_status = pge_errstat_error
      WRITE(www_lun, '(A)') 'Need to increase nxtrack_max!!!'
    ENDIF

    RETURN
  END SUBROUTINE omi_set_parameters

  ! Find the scan line/x track position based on inut lat/lon range
  ! Avoid including descending orbits
  SUBROUTINE find_scan_line_range ( slat, elat, slon, elon, sline, eline, &
       spix, epix, pge_error_status )

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY: l1b_rad_filename, szamax, ntimes
    USE OMSAO_omidata_module,   ONLY: omi_radiance_swathname, &
         nxtrack_max, ntimes_max,  nfxtrack!, nxtrack, nswath
    USE L1B_Reader_class
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! -----------------------
    ! Input/Output variables
    ! -----------------------
    REAL (KIND=dp), INTENT(IN) :: slat, elat, slon, elon
    INTEGER,      INTENT (OUT) :: pge_error_status, sline, eline, spix, epix

    REAL (KIND=r4), DIMENSION (:,:), POINTER :: lons, lats, szas
    REAL (KIND=r4), DIMENSION (:),   POINTER :: latdf
    LOGICAL, DIMENSION(:,:), POINTER         :: any_pixels
    REAL (KIND=dp)                           :: tmplon, dlat
    INTEGER, DIMENSION (:,:), POINTER        :: okline
    INTEGER, DIMENSION (:), POINTER          :: lines
    INTEGER                                  :: iline, errstat, ix

    ! Exteranl functions
    INTEGER                    :: estat

    TYPE (L1B_block_type)      :: omi_data_block   

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=20), PARAMETER :: modulename = 'find_scan_line_range'

    ! --------------------------
    ! Initialize OUTPUT variable
    ! --------------------------
    allocate (lons(nxtrack_max, 0:ntimes_max), & 
              lats(nxtrack_max, 0:ntimes_max), & 
              szas(nxtrack_max, 0:ntimes_max))
    allocate (any_pixels (nxtrack_max, 0:ntimes_max))
    allocate (latdf (ntimes_max))
    allocate (okline(1:nxtrack_max, 2), lines(1:ntimes_max))
    pge_error_status = pge_errstat_ok
 

    CALL omi_set_parameters ( pge_error_status )

    IF ( pge_error_status >= pge_errstat_error ) RETURN

    ! -----------------------------------------------------------
    ! Open data block (UV-1, if both are selected) 
    ! called 'omi_data_block' with nTimes lines
    ! -----------------------------------------------------------
    errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(1))
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", modulename, 0 )
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! ---------------------------------
    ! Read all Latitudes and Longitudes
    ! ---------------------------------
    sline = -5
    eline = -5
    spix =-5
    epix =-5
    DO iline = 0, ntimes - 1    
      errstat = L1Br_getGEOline ( omi_data_block, iline,  &
           Latitude_k              = lats (1:nfxtrack, iline), & 
           Longitude_k             = lons (1:nfxtrack, iline), &
           SolarZenithANgle_k      = szas (1:nfxtrack, iline))  

      IF( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
             "L1Br_getGEOline failed.", modulename, 0 )
        pge_error_status = pge_errstat_error
        RETURN
      END IF
    ENDDO

    lines(1:ntimes-1) = (/(iline, iline=1, ntimes-1)/)
    DO ix = 1, nfxtrack 
      latdf(1:ntimes-1) = lats(ix, 1:ntimes-1) - lats(ix, 0:ntimes-2)
      dlat  = latdf(ntimes / 2)

      okline(ix, 1) = MINVAL(MINLOC( lines(1:ntimes-1), MASK = (latdf(1:ntimes-1) > dlat/2.)))
      okline(ix, 2) = MINVAL(MAXLOC( lines(1:ntimes-1), MASK = (latdf(1:ntimes-1) > dlat/2.)))
      !print *, ix, okline(ix, 1), okline(ix, 2), ntimes, dlat
    ENDDO

    DO iline = 0, ntimes - 1      
      IF (ANY(lats(1:nfxtrack, iline) >= slat) .AND. ANY(lats(1:nfxtrack, iline) <= elat)) THEN
        any_pixels(1:nfxtrack, iline) = .FALSE.
        DO ix = 1, nfxtrack
          tmplon = lons(ix, iline)
          IF (elon > 180 .AND. tmplon < 0) tmplon = tmplon + 360.
          IF (tmplon >= slon .AND. tmplon <= elon .AND. szas(ix, iline) <= szamax &
               .AND. iline >= okline(ix, 1) .AND. iline <= okline(ix, 2)) any_pixels(ix, iline) = .TRUE.
          !IF (tmplon >= slon .AND. tmplon <= elon .AND. szas(ix, iline) <= szamax) any_pixels(ix, iline) = .TRUE.
        ENDDO
        IF ( ANY(any_pixels(1:nfxtrack, iline)) ) THEN
          IF (sline < 0) sline = iline
          !WRITE(www_lun, '(I5,60L2)') iline, any_pixels(1:nfxtrack, iline)
          !WRITE(www_lun, '(I5,60F6.1)') iline, szas(1:nfxtrack, iline)
          eline = iline
        ENDIF
      ENDIF
    ENDDO

    ! --------------------------
    ! Close data block structure
    ! --------------------------
    errstat = L1Br_close ( omi_data_block )
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
           "L1Br_close failed.", modulename, 0 )
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    IF (sline >=0 .AND. eline >= sline) THEN
      DO ix = 1, nfxtrack
        IF (ANY(any_pixels(ix, sline:eline))) THEN
          IF (spix < 0) spix = ix
          epix = ix
        ENDIF
      ENDDO
    ENDIF

    deallocate( lons, lats, szas, any_pixels, latdf, okline, lines)
    !WRITE(www_lun, '(4F8.2)') slat, elat, slon, elon
    !WRITE(www_lun, '(4I8)')  sline, eline, spix, epix
    !STOP

    RETURN
  END SUBROUTINE find_scan_line_range

  SUBROUTINE omi_read_irradiance_data (lun, nxcoadd, first_pix, last_pix, &
       pge_error_status ) 

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, maxwin
    USE OMSAO_parameters_module, ONLY: maxchlen, maxwin, max_ring_pts, &
         mrefl, vb_lev_omidebug, mswath, max_fit_pts
    USE OMSAO_variables_module,  ONLY: verb_thresh_lev, l1b_irrad_filename, &
         wcal_bef_coadd, currpix, numwin, coadd_uv2, band_selectors, winpix, &
         winlim, scnwrt, refdbdir, use_backup, &
         reduce_resolution, redlam, redsampr, reduce_slit, rm_mgline, &
         dwavmax, use_redfixwav, which_slit, avgsol_allorb, &
         nxbin,nswath, GranuleJDay, inschs, nxtrack,ncoadd,&
         reduce_ubnd, reduce_lbnd, retlbnd, retubnd, redslw, earthsundistance
    USE OMSAO_omidata_module,  ONLY: nwavel_max,  omi_irradiance_swathname, &
                                      nfxtrack, zoom_p1
    USE ozprof_data_module, ONLY: pos_alb, toms_fwhm, nrefl
    USE hdfeos4_parameters
    USE L1B_Reader_class
    USE OMSAO_errstat_module
    use m_gauss, only: gauss_uneven
    use m_triangle, only: triangle_uneven
    use tell_module


    IMPLICIT NONE
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER, INTENT (IN)        :: nxcoadd, first_pix, last_pix, lun
    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER , INTENT (OUT)      :: pge_error_status
    ! ---------------
    ! Local variables
    ! ---------------
    TYPE (L1B_block_type)                   :: omi_data_block
    INTEGER, PARAMETER :: noff_uv1=12, noff_uv2=25
    INTEGER  :: nwavel, is, ix, i, j, iix, nomi, fidx, lidx, ch, idum,  &
         iw, ic, idx, noff1, noff2, nring, irefl, nbin,  &
         thedoy, nsolbin, nbad
    ! variables used to read original spectra
    INTEGER (KIND=i2)                     :: mflg
    INTEGER (KIND=i4), DIMENSION(mswath)  :: nwls
    INTEGER   (KIND=i2), DIMENSION (mswath) :: spos, epos
    INTEGER, DIMENSION (:), POINTER         :: idxs
    INTEGER (KIND=i2), DIMENSION(:,:), POINTER :: tmpnavg
    INTEGER (KIND=i2), DIMENSION(:, :), POINTER :: irrad_qflg
    REAL (kind=4),DIMENSION (:, :), POINTER     :: irrad_spec, irrad_prec, irrad_wavl
    INTEGER (KIND=i4)                     :: nx, nt
    ! variables used for reduce resoltuion
    INTEGER :: npos, np
    REAL (KIND = dp)  :: tmpsampr, retswav, retewav!, fdum
    INTEGER (KIND=i2), DIMENSION (:,:), POINTER :: tmpqflg
    REAL (KIND = dp), DIMENSION (:,:,:),POINTER :: tmpspec
    INTEGER (kind=i2), dimension(1) :: temp_mflg
    ! Subset variables 
    INTEGER (KIND=i4), DIMENSION(maxwin)   :: nwbin
    INTEGER, PARAMETER                      :: nbits = 16
    INTEGER (KIND=i2), DIMENSION(0:nbits-1) :: mflgbits
    INTEGER (KIND=i2), DIMENSION(:,:,:), POINTER :: flgbits
    INTEGER (KIND=i2), DIMENSION(nwavel_max)     :: flgmsks
    REAL (KIND = dp)  :: wcenter, normsc
    REAL (KIND = dp), DIMENSION (maxwin, nxcoadd) :: wshis, wsqus
    REAL (KIND = dp), DIMENSION (:,:,:),POINTER :: subspec
    REAL (KIND = dp), DIMENSION (:,:), POINTER    :: subring 

    INTEGER (KIND=i4) :: errstat, iline
    LOGICAL                                 :: error
    CHARACTER (LEN=maxchlen)                :: bkfname
    ! Exteranl functions
    INTEGER                                 :: estat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=24), PARAMETER :: modulename = 'omi_read_irradiance_data' 

    !--------------------------------------------------------------------------
    ! Starting with allocating local variables
    !--------------------------------------------------------------------------
    allocate (irrad_qflg (nwavel_max, nxtrack_max))
    allocate (irrad_prec (nwavel_max, nxtrack_max))
    allocate (irrad_spec (nwavel_max, nxtrack_max))
    allocate (irrad_wavl (nwavel_max, nxtrack_max))
    allocate (idxs (nwavel_max))
    allocate (flgbits(nxcoadd, nwavel_max, 0:nbits-1))
    allocate (tmpspec (sig_idx, nwavel_max, nxtrack_max))
    allocate (tmpqflg(nwavel_max, nxtrack_max))
    allocate (tmpnavg(nwavel_max, nxtrack_max))
    allocate (subspec(nxcoadd*2, sig_idx, max_fit_pts))
    allocate (subring(sig_idx, max_ring_pts))

    j = 1
    nwavel = 0
    iline = 0
    pge_error_status = pge_errstat_ok
    errstat = omi_s_success

    ! For zoom-in global products (once every 32 days), solar irradiance (not in radiance)
    ! in UV1 are provided at 60 Xtrack positions. Under this condition, solar irradiance needs
    ! to be rebinned additionally to the normal.
    ! In UV-2, also provided at 60 across-track position, but corresponding to positions 16-45
    ! in normal mode after coadd every two spectra
    nsolbin = 1  

    ! ----------------------------
    ! Initialize irradiance arrays
    ! ----------------------------  
    omi_irrad%errstat(1:nxtrack) = pge_errstat_ok
    omi_irrad%nwav (1:nxtrack) = 0
    omi_irrad%npix (1:numwin, 1:nxtrack) = 0
    omi_irrad%prec (:,1:nxtrack) = 0.0
    omi_irrad%spec (:,1:nxtrack) = 0.0
    omi_irrad%wavl (:,1:nxtrack) = 0.0
    omi_irrad%qflg (:,1:nxtrack) = 0

    IF (.NOT. use_backup) THEN 
      DO is = 1, nswath
        ch = inschs(is)
        ! ------------------------------------------------------
        ! Open data block structure with default size of 1 lines
        ! ------------------------------------------------------
         errstat = L1Br_open ( omi_data_block, l1b_irrad_filename, &
                              TRIM(ADJUSTL(omi_irradiance_swathname(ch))) )
         IF( errstat /= omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", modulename, 0)
         STOP 1

        ! ----------------------------------
        ! Obtain irradiance swath dimensions
        ! ----------------------------------
          errstat = L1Br_getSWdims ( omi_data_block, NumTimes_k=nt, nXtrack_k=nx)!, nWavel_k=nw, nWavelCoef_k=nwc)
          IF( errstat /=  omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "IL1Br_getSWdims failed.", modulename, 0)
            STOP 1
          END IF
        endif

        ! ----------------------------------------------------------
        ! Obtain time, geolocation, and angular information on block
        ! ----------------------------------------------------------
          errstat = L1Br_getDATA ( omi_data_block, iline, &
               MeasurementQualityFlags_k = mflg)
          IF( errstat /= omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, &
                 "L1Br_getDATA failed.", modulename, 0)
            STOP 1

          errstat = L1Br_getSIGline ( omi_data_block, iline,     &
               Signal_k            = irrad_spec(j:, :), &
               SignalPrecision_k   = irrad_prec(j:, :), &
               PixelQualityFlags_k = irrad_qflg(j:, :), &
               Wavelength_k        = irrad_wavl(j:, :), &
               !NumberSmallPixelColumns_k = tmpNinteg,           &
               Nwl_k  = nwls(ch) )
          IF( errstat /= omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, &
                 "L1Br_getSIGline failed.", modulename, 0)
            STOP 1
          END IF
        endif

        temp_mflg=mflg
        CALL convert_2bytes_to_16bits ( nbits, 1, temp_mflg, mflgbits(0:nbits-1))

        IF (mflgbits(0) == 1 .OR. mflgbits(1) == 1 .OR. mflgbits(3) == 1 .OR. mflgbits(12) == 1) THEN
          WRITE(www_lun, *) 'All irradiances could not be used, use backup irradiances: ', &
               TRIM(ADJUSTL(omi_irradiance_swathname(ch)))
          omi_irrad%errstat(1:nxtrack) = pge_errstat_error
          pge_error_status = pge_errstat_error
          RETURN
        ELSE IF (ANY(mflgbits == 1)) THEN
          WRITE(www_lun, *) 'Warning set on all irradiances: ', TRIM(ADJUSTL(omi_irradiance_swathname(ch)))
          omi_irrad%errstat(1:nxtrack) = pge_errstat_error  
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF

        nwavel = nwavel + nwls(ch)
        spos(ch) = int(j, kind=i2)
        j = nwavel + 1
        epos(ch) = int(nwavel, kind=i2)

        ! --------------------------
        ! Close data block structure
        ! --------------------------
          errstat = L1Br_CLOSE ( omi_data_block )
          IF( errstat /= omi_s_success .AND. verb_thresh_lev >= &
               vb_lev_omidebug ) THEN
            estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, &
                 "L1Br_CLOSE failed.", modulename, 0)
            STOP 1
          END IF


        ! Need to sort the data in increasing wavelength
        IF (irrad_wavl(spos(ch), 1) > irrad_wavl(epos(ch), 1)) THEN
          idxs(spos(ch):epos(ch)) = (/ (i, i = epos(ch), spos(ch), -1) /)
          irrad_wavl(spos(ch):epos(ch), :) = irrad_wavl(idxs(spos(ch):epos(ch)), :)
          irrad_spec(spos(ch):epos(ch), :) = irrad_spec(idxs(spos(ch):epos(ch)), :)
          irrad_prec(spos(ch):epos(ch), :) = irrad_prec(idxs(spos(ch):epos(ch)), :)
          irrad_qflg(spos(ch):epos(ch), :) = irrad_qflg(idxs(spos(ch):epos(ch)), :)     
        ENDIF

        !IF ( ch == 2 ) CALL corruv2wav(nwls(ch), nx, irrad_wavl(spos(ch):epos(ch), 1:nx)) 

        !OPEN(unit=90, FILE='/data/dumbo/xliu/OMIHCLD/OMIL1BBIRR-o05168_vis.dat', STATUS='old')
        !WRITE(90, *) nx, nw
        !DO ix = 1, nx 
        !   WRITE(90, *) ix
        !   DO iw = spos(ch), epos(ch)
        !      !CALL convert_2bytes_to_16bits ( nbits, 1, irrad_qflg(iw, ix), mflgbits(0:nbits-1))
        !      WRITE(90, '(F10.4,D14.6,1X)') irrad_wavl(iw, ix), irrad_spec(iw, ix) !, &
        !      !mflgbits(0:nbits-1)
        !   ENDDO
        !ENDDO
      ENDDO

    ! Use backup solar spectrum
    ELSE
      ! Determine sun-earth distance correction
      thedoy=GranuleJDay
      IF (thedoy == 366) thedoy = 365
     
      OPEN (UNIT=lun, FILE= ADJUSTL(TRIM(refdbdir)) // 'solar-distance.dat', &
           STATUS='UNKNOWN', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Cannot open Sun-Earth Distance datafile!!!'
        pge_error_status = pge_errstat_error
        RETURN
      END IF
      DO i = 1, 12
        READ(LUN, *)
      ENDDO
      DO i = 1, thedoy
        READ(LUN, *) normsc, normsc
      ENDDO
      earthsundistance = normsc
      CLOSE(LUN)
      normsc = 1.0 / normsc ** 2  ! solar energy is inversely proportional to square distance
      IF (avgsol_allorb) THEN
        bkfname = ADJUSTL(TRIM(refdbdir))//'OMI/omisol_v003_avg_nshi_backup.dat'
      ELSE
        bkfname = ADJUSTL(TRIM(refdbdir)) // 'omisolrunavg/omisol_'//omiraddate // '_v3_31runavg_backup.dat'
      ENDIF
      IF( scnwrt ) WRITE(*,*) 'use_backup=(T):'//ADJUSTL(TRIM( bkfname ))
      OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(bkfname)), STATUS='UNKNOWN', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Cannot open solar backup file!!!'
        pge_error_status = pge_errstat_error
        RETURN
      END IF

      ! Open the file and read all the data
      nwavel = 0
      DO is = 1, 2
        READ(lun, *) nx, nwls(is)
        spos(is) = int(nwavel + 1, kind=i2)
        epos(is) = int(nwavel + nwls(is) , kind=i2)
        DO i = 1, nx
          READ(lun, *) 
          IF (avgsol_allorb) THEN
            DO j = 1, nwls(is)
              READ(lun, *) irrad_wavl(nwavel + j, i),  irrad_spec(nwavel + j, i),&
                           irrad_prec(nwavel + j, i), idum, idum
              IF (idum > 0) irrad_prec(nwavel + j, i) = &
                             real(irrad_prec(nwavel + j, i) &
                             / SQRT( REAL(idum, KIND=dp) ) , kind=r4)
            ENDDO
          ELSE
            DO j = 1, nwls(is)
              READ(lun, *) irrad_wavl(nwavel + j, i), irrad_spec(nwavel + j, i), &
                           irrad_prec(nwavel + j, i), idum, tmpnavg(nwavel +j,i)
            ENDDO
            idum = MAXVAL(tmpnavg(spos(is):epos(is), i))
            DO j = 1, nwls(is)
               IF (tmpnavg(nwavel + j, i) > 10) THEN
                   irrad_prec(nwavel + j, i) = real( irrad_prec(nwavel+j, i) &
                   / SQRT( REAL(tmpnavg(nwavel + j, i), KIND=dp) ), kind=r4)
               ELSE
                 irrad_spec(nwavel + j, i) = 0.0
                 irrad_prec(nwavel + j, i) = 0.0
               ENDIF
            ENDDO
          ENDIF
        ENDDO
        irrad_spec(spos(is):epos(is), 1:nx) = real(irrad_spec(spos(is):epos(is), 1:nx) *normsc , kind=r4)
        irrad_prec(spos(is):epos(is), 1:nx) = real(irrad_prec(spos(is):epos(is), 1:nx) *normsc , kind=r4)
        nwavel = epos(is)
        !IF ( is == 2 ) CALL corruv2wav(nwls(is), nx, irrad_wavl(spos(is):epos(is), 1:nx)) 
      ENDDO
      IF (nswath == 1) THEN
        IF (band_selectors(1) == 1) THEN
          nwavel = nwls(1)         
        ELSE
          nwavel = nwls(2)
          irrad_wavl(1:nwavel, :) = &
               irrad_wavl(spos(2):epos(2), :)
          irrad_spec(1:nwavel, :) = &
               irrad_spec(spos(2):epos(2), :) 
          irrad_prec(1:nwavel, :) = &
               irrad_prec(spos(2):epos(2), :)
          spos(2) = int(spos(2) - nwls(1) , kind=i2)
          epos(2) = int(epos(2) - nwls(1) , kind=i2)
        ENDIF
      ENDIF
      CLOSE(LUN)
    ENDIF
    !  ! xliu: Feb/19/2008, read straylight spectra
    !  straylight_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omi_irrad_sl_v3.dat' 
    !  OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(straylight_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
    !  IF ( errstat /= pge_errstat_ok ) THEN
    !     WRITE(www_lun, '(2A)') modulename, ': Cannot open irradiance straylight file!!!'
    !     pge_error_status = pge_errstat_error; RETURN
    !  END IF
    !  
    !  ! Open the file and read all the data
    !  nwavel = 0
    !  DO is = 1, 2
    !     READ(lun, *) nx, nwls(is)
    !     spos(is) = nwavel + 1; epos(is) = nwavel + nwls(is)
    !     DO i = 1, nx
    !        READ(lun, *) 
    !        DO j = 1, nwls(is)
    !           READ(lun, *) fdum, omi_irrad_stray(nwavel + j, i)
    !        ENDDO
    !     ENDDO
    !     nwavel = epos(is)
    !  ENDDO
    !
    !  IF (nswath == 1) THEN
    !     IF (band_selectors(1) == 1) THEN
    !        nwavel = nwls(1)         
    !     ELSE
    !        nwavel = nwls(2)
    !        omi_irrad_stray(1:nwavel, :) = omi_irrad_stray(spos(2):epos(2), :) 
    !        spos(2) = spos(2) - nwls(1); epos(2) = epos(2) - nwls(1) 
    !     ENDIF
    !  ENDIF
    !  CLOSE(LUN)
    !
    !  straylight_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omi_rad_sl_v3.dat' 
    !  OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(straylight_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
    !  IF ( errstat /= pge_errstat_ok ) THEN
    !     WRITE(www_lun, '(2A)') modulename, ': Cannot open radiance straylight file!!!'
    !     pge_error_status = pge_errstat_error; RETURN
    !  END IF
    !  
    !  ! Open the file and read all the data
    !  nwavel = 0
    !  DO is = 1, 2
    !     READ(lun, *) nx, nwls(is)
    !     spos(is) = nwavel + 1; epos(is) = nwavel + nwls(is)
    !     DO i = 1, nx
    !        READ(lun, *) 
    !        DO j = 1, nwls(is)
    !           READ(lun, *) fdum, omi_rad_stray(nwavel + j, i)
    !        ENDDO
    !     ENDDO
    !     nwavel = epos(is)
    !  ENDDO
    !
    !  IF (nswath == 1) THEN
    !     IF (band_selectors(1) == 1) THEN
    !        nwavel = nwls(1)         
    !     ELSE
    !        nwavel = nwls(2)
    !        omi_irrad_stray(1:nwavel, :) = omi_rad_stray(spos(2):epos(2), :) 
    !        spos(2) = spos(2) - nwls(1); epos(2) = epos(2) - nwls(1) 
    !     ENDIF
    !  ENDIF
    !  CLOSE(LUN) 

    ! Do not coadd wavelengths with a gap (e.g., filter Mg absorption lines), need to determine
    ! delta-lamda in UV-1
    ! Note in OMI delta-lamda varies with wavelength (largest for the first two pixels in each channel)
    dwavmax = ( irrad_wavl(2, 1) - irrad_wavl(1, 1) ) * 1.1

    ! Degrade spectral resolution if necessary
    IF (reduce_resolution) THEN
      nwavel = 0
      j = 1
      DO is = 1, nswath
        ch = inschs(is)
        IF (coadd_uv2 .AND. is == 1) THEN
          npos = nfxtrack * nsolbin
        ELSE
          npos = nxtrack
        ENDIF
        !IF (nswath == 2) THEN
        !   retswav = retlbnd(ch); retewav = retubnd(ch)
        !ELSE
        !   retswav = retlbnd(band_selectors(1)); retewav = retubnd(band_selectors(1))
        !ENDIF
        retswav = retlbnd(ch)
        retewav = retubnd(ch)

        tmpsampr = redsampr 
        IF (is == 1 .AND. band_selectors(1) == 1) tmpsampr = redsampr / 3.0
        np = nwls(ch)

        tmpspec(wvl_idx, 1:np, 1:npos) = irrad_wavl(spos(ch):epos(ch), 1:npos)
        tmpspec(spc_idx, 1:np, 1:npos) = irrad_spec(spos(ch):epos(ch), 1:npos)
        tmpspec(sig_idx, 1:np, 1:npos) = irrad_prec(spos(ch):epos(ch), 1:npos)
        tmpqflg(   1:np, 1:npos) = irrad_qflg(spos(ch):epos(ch), 1:npos)       

        DO ix = 1, npos
          CALL convert_2bytes_to_16bits ( nbits, np, tmpqflg(1:np, ix), flgbits(1, 1:np, 0:nbits-1))
          tmpqflg(1:np, ix) = flgbits(1, 1:np, 0) &   ! Missing
               + flgbits(1, 1:np, 1) &   ! Bad
               + flgbits(1, 1:np, 2)     ! Processing error
        ENDDO

        CALL reduce_irrad_resolution (tmpspec(:, 1:np, 1:npos), &
             tmpqflg(1:np, 1:npos), np, npos, reduce_slit, redslw(is), &
             tmpsampr, redlam, retswav, retewav, reduce_lbnd(ch), &
             reduce_ubnd(ch), nwls(ch), pge_error_status)

        IF (pge_error_status == pge_errstat_error) RETURN

        nwavel = nwavel + nwls(ch)
        spos(ch) = int(j, kind=i2)
        j = nwavel + 1
        epos(ch) = int(nwavel, kind=i2)
        irrad_wavl(spos(ch):epos(ch), 1:npos) = &
             real(tmpspec(wvl_idx, 1:nwls(ch), 1:npos), kind=r4)
        irrad_spec(spos(ch):epos(ch), 1:npos) = &
             real(tmpspec(spc_idx, 1:nwls(ch), 1:npos), kind=r4)
        irrad_prec(spos(ch):epos(ch), 1:npos) = &
             real(tmpspec(sig_idx, 1:nwls(ch), 1:npos), kind=r4)
        irrad_qflg(spos(ch):epos(ch), 1:npos) = 0   ! All data are good  (pre filtered) 
      ENDDO
    ENDIF

    IF (nwavel > nwavel_max) THEN
      WRITE(www_lun, '(A)') "Need to increase nwavel_max!!!"
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF

    ! Determine number of wavelengths to be read for deteriming cloud fraction
    fidx = MAXVAL ( MINLOC ( irrad_wavl(1:nwavel, 1), MASK = &
         (irrad_wavl(1:nwavel, 1) > pos_alb - toms_fwhm * 1.4) ))
    lidx = MAXVAL ( MAXLOC ( irrad_wavl(1:nwavel, 1), MASK = &
         (irrad_wavl(1:nwavel, 1) < pos_alb + toms_fwhm * 1.4) ))
    IF (fidx <1 .OR. lidx > nwavel) THEN
      WRITE(www_lun, '(2A)') modulename, ': Need to change pos_alb/toms_fwhm!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF
    nrefl = lidx - fidx + 1 
    IF (nrefl > mrefl ) THEN
      WRITE(www_lun, '(2A)') modulename, ': Need to increase mrefl!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF

    ! Determine number of binning for different fitting windows
    DO iw = 1, numwin
      ch = band_selectors(iw)
      IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
        nwbin(iw) = nxbin 
        IF (ch == 1) nwbin(iw) = nxbin 
      ELSE
        nwbin(iw) = nxbin * ncoadd 
      ENDIF
      nwbin(iw) = nwbin(iw) * nsolbin
    ENDDO

    ! Subset and coadd irradiance spectrum
    DO ix = first_pix, last_pix
      currpix = ix          ! Global variable, to be used in wavlength calibraiton with omi slit before coadding

      ! Indices for UV-1 are from 1 to nfxtrack
      ! Indices for UV-2 are from 1 to nxtrack (nfxtrack * 2)
      ! UV2 indices corresponding to UV-1 pixel ix are ix * 2 -1 & ix * 2 (or iix+1, iix+2), respectively
      ! If additional across track coadding (e.g., nxbin) is performed, then for a particular ix
      ! UV1: nbin = nxbin; UV-2: nbin = nxbin * 2
      ! Coadded original across track pixels are: iix + 1 : iix + nbin (iix = (ix-1) * nbin)

      ! Get quality flag bits, coadd flags if necessary to avoid coadding inconsistent # of pixels 
      flgmsks = 0
      DO is = 1, nswath
        ch = inschs(is)

        !Do not use nwbin
        IF (is == 1) THEN
          nbin = nxbin 
        ELSE
          nbin = nxbin * ncoadd
        ENDIF
        nbin = nbin * nsolbin  ! Zoom mode UV1 solar irradiance got 60 positions
        iix = (ix - 1) * nbin

        ! Shift the position by 15 
        IF (ch == 2 .AND. nsolbin == 2) THEN
          iix = iix - (zoom_p1 - 1) * nsolbin
        ENDIF

        IF (.NOT. reduce_resolution) THEN
          ! properly align cross track positions to be coadded (should be within one pixel)
          IF (nbin / nsolbin > 2) CALL prespec_align(nwls(ch), nbin, &
               irrad_wavl(spos(ch):epos(ch), iix+1:iix+nbin), &
               irrad_spec(spos(ch):epos(ch), iix+1:iix+nbin), &
                                !omi_irrad_stray(spos(ch):epos(ch), iix+1:iix+nbin), &
                                !omi_rad_stray(spos(ch):epos(ch), iix+1:iix+nbin), &
               irrad_prec(spos(ch):epos(ch), iix+1:iix+nbin), &       
               irrad_qflg(spos(ch):epos(ch), iix+1:iix+nbin))

          DO ic = 1, nbin
            CALL convert_2bytes_to_16bits ( nbits, nwls(ch), irrad_qflg(spos(ch):epos(ch), iix + ic ), &
                 flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
            flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) &
                 + flgbits(ic, spos(ch):epos(ch), 0)                &   ! Missing
                 + flgbits(ic, spos(ch):epos(ch), 1)                &   ! Bad 
                 + flgbits(ic, spos(ch):epos(ch), 2)                    ! Processing error
          ENDDO

          !DO i = spos(is), epos(is)
          !   WRITE(*, '(F10.4, D14.6, 16I2)') irrad_wavl(i, iix+1), &
          !        irrad_spec(i, iix+1), flgbits(1, i, 0:nbits-1)
          !ENDDO

        ELSE
          ! Already aligned because of using common wavelength scale
          DO ic = 1, nbin
            flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) + &
                 irrad_qflg(spos(ch):epos(ch), iix + ic )
          ENDDO
        ENDIF
      ENDDO
      ! Subset valid data
      nomi = 0
      subspec = 0.0
      ! strayspec = 0.0
      DO iw = 1, numwin
        ch = band_selectors(iw)
        nbin = nwbin(iw)
        iix = (ix - 1) * nbin

        ! Shift the position by 15 
        IF (ch == 2 .AND. nsolbin == 2) THEN
          iix = iix - (zoom_p1 - 1) * nsolbin
        ENDIF

        winpix(iw, 1) = MINVAL ( MINLOC ( irrad_wavl(spos(ch):epos(ch),iix + 1), &
             MASK = irrad_wavl(spos(ch):epos(ch),iix + 1) >= &
             winlim(iw, 1)) ) + spos(ch) - 1
        winpix(iw, 2) = MAXVAL ( MAXLOC ( irrad_wavl(spos(ch):epos(ch),iix + 1), &
             MASK = irrad_wavl(spos(ch):epos(ch),iix + 1) <= &
             winlim(iw, 2)) ) + spos(ch) - 1

        IF (winpix(iw, 1) .le. 0) winpix(iw,1) = 1
        IF (winpix(iw, 2) .gt. nwavel) winpix(iw, 2) = nwavel

        omi_irrad%winpix(iw, ix, 1:2) = 0       
        omi_irrad%npix  (iw, ix) = nomi
        fidx = winpix(iw, 1) 
        lidx = winpix(iw, 2) 

        IF (rm_mgline .AND. winlim(iw, 1) < 286.0 .AND. &
             winlim(iw, 2) > 286.0) THEN
          DO i = fidx, lidx
            IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. &
                flgmsks(i) == 0 .AND. &
                 !(ALL(irrad_wavl(i, iix+1:iix+nbin) < 273.8)   .OR. &
                !ALL(irrad_wavl(i, iix+1:iix+nbin)  > 275.2))  .AND. &
                 (ALL(irrad_wavl(i, iix+1:iix+nbin) < 278.8)   .OR. &
                 ALL(irrad_wavl(i, iix+1:iix+nbin)  > 281.0))  .AND. &
                 (ALL(irrad_wavl(i, iix+1:iix+nbin) < 284.7)   .OR. &
                 ALL(irrad_wavl(i, iix+1:iix+nbin)  > 285.7))) THEN
              !(ALL(irrad_wavl(i, iix+1:iix+nbin) < 278.0)   .OR. &
              !ALL(irrad_wavl(i, iix+1:iix+nbin)  > 282.0))  .AND. &
              !(ALL(irrad_wavl(i, iix+1:iix+nbin) < 284.0)   .OR. &
              !ALL(irrad_wavl(i, iix+1:iix+nbin)  > 286.0))) THEN
              nomi = nomi + 1
              subspec(1:nbin, wvl_idx, nomi) =irrad_wavl(i, iix+1:iix+nbin)
              subspec(1:nbin, spc_idx, nomi) =irrad_spec(i, iix+1:iix+nbin)
              subspec(1:nbin, sig_idx, nomi) =irrad_prec(i, iix+1:iix+nbin)
           !strayspec(1:nbin, 1, nomi)     = omi_irrad_stray(i, iix+1:iix+nbin)
           !strayspec(1:nbin, 2, nomi)     = omi_rad_stray(i, iix+1:iix+nbin)
              IF (omi_irrad%winpix(iw, ix, 1) == 0) omi_irrad%winpix(iw, ix, 1) = i
              omi_irrad%winpix(iw, ix, 2) = i
              omi_irrad%wind(nomi, ix) = int(i, kind=i2)
            ENDIF
          ENDDO
        ELSE
         DO i = fidx, lidx

            IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                 ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. &
                 flgmsks(i) == 0 ) THEN
              nomi = nomi + 1
              subspec(1:nbin, wvl_idx, nomi) =irrad_wavl(i, iix+1:iix+nbin)
              subspec(1:nbin, spc_idx, nomi) =irrad_spec(i, iix+1:iix+nbin)
              subspec(1:nbin, sig_idx, nomi) =irrad_prec(i, iix+1:iix+nbin)
           !strayspec(1:nbin, 1, nomi)     = omi_irrad_stray(i, iix+1:iix+nbin)
           !strayspec(1:nbin, 2, nomi)     = omi_rad_stray(i, iix+1:iix+nbin)
              IF (omi_irrad%winpix(iw, ix, 1) == 0) omi_irrad%winpix(iw, ix, 1) = i
              omi_irrad%winpix(iw, ix, 2) = i
              omi_irrad%wind(nomi, ix) = int(i, kind=i2)
            ENDIF
          ENDDO
        ENDIF

        omi_irrad%npix(iw, ix) = nomi - omi_irrad%npix(iw, ix)

        !WRITE(www_lun, '(2I5, 2F8.3, 2I5, 2F8.3, 6I5)') ix, iw, irrad_wavl(spos(ch), iix+1), &
        !     irrad_wavl(epos(ch), iix+1), spos(ch), epos(ch), winlim(iw, 1), winlim(iw, 2), &
        !     winpix(iw, 1), winpix(iw, 2), fidx, lidx, lidx - fidx + 1, omi_irrad%npix(iw, ix)
      ENDDO
      omi_irrad%nwav(ix) = nomi

      ! Perform coadding when necessary
 
      fidx = 1
      DO iw = 1, numwin       
        ch = band_selectors(iw)
        nbin = nwbin(iw)
        lidx = fidx + omi_irrad%npix(iw, ix) - 1 
        IF (nbin > 1) THEN
          CALL solwavcal_coadd(wcal_bef_coadd, omi_irrad%npix(iw, ix), nbin, &
               subspec(1:nbin, :, fidx:lidx), wshis(iw, 1:nbin), wsqus(iw, 1:nbin), error)
          IF (error) THEN
            WRITE(www_lun, '(A)') 'No solar wavelength calibration before coadding!!!'
            omi_irrad%errstat(ix) = pge_errstat_warning
          ENDIF
        ENDIF
        fidx = lidx + 1    
      ENDDO
      ! Subset solar spectrum for Ring effect
      IF (.NOT. reduce_resolution .OR. (reduce_resolution .AND. .NOT. use_redfixwav)) THEN
        subring = 0.0
        ch = band_selectors(1)
        IF (band_selectors(1) == 1) THEN
          noff1 = noff_uv1
        ELSE
          noff1 = noff_uv2
        ENDIF
        nring = omi_irrad%npix(1, ix) + noff1
        subring(1:spc_idx, noff1+1 : nring) = subspec(1, 1:spc_idx, 1:omi_irrad%npix(1, ix)) 
        omi_ring%ndiv(ix) = 0

        ! add extra spectra before first window (uncoadded)
        ! if unavailable, needed to ammened with solar reference spectrum
        noff1 = noff1 + 1 
        nbin = nwbin(1)
        iix = (ix - 1) * nbin
        
        IF (ch == 2 .AND. nsolbin == 2) THEN  ! Shift the position by 15 
          iix = iix - (zoom_p1 - 1) * nsolbin
        ENDIF
        ! jbak 2017-12-25, there are very huge bad pixels for near upper boundary of
        ! shorther wavelengths
        ! wavelengths ==> causing fail because reference solar spectrum does not
        ! cover in apending ring
        ! therefore there is 5 continuous bad pixels then exit
        nbad = 0
        DO i = winpix(1, 1) - 1, 1, -1
          IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
               ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
            noff1 = noff1 - 1
            subring(wvl_idx, noff1) = SUM(irrad_wavl(i, iix+1:iix+nbin)) / nbin
            subring(spc_idx, noff1) = SUM(irrad_spec(i, iix+1:iix+nbin)) / nbin
           ! print * , i, subring(wvl_idx, noff1)
            IF (noff1 == 1) EXIT
          ELSE ! JBAK
            nbad = nbad + 1
            IF (nbad == 5) EXIT
          ENDIF
        ENDDO
        DO iw = 2, numwin
          ch = band_selectors(iw)
          nbin = nwbin(iw - 1) 
          iix = (ix - 1) * nbin
          IF (ch == 2 .AND. nsolbin == 2) THEN  ! Shift the position by 15 
            iix = iix - (zoom_p1 - 1) * nsolbin
          ENDIF
          IF (ch == band_selectors(iw - 1) ) THEN   
            DO i = winpix(iw-1, 2)+1, winpix(iw, 1)-1 
              IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                   ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
                nring = nring + 1
                subring(wvl_idx, nring) = SUM(irrad_wavl(i, iix+1:iix+nbin)) / nbin
                subring(spc_idx, nring) = SUM(irrad_spec(i, iix+1:iix+nbin)) / nbin
                !print * , i,nring, subring(wvl_idx, nring)
              ENDIF
            ENDDO
          ELSE  ! first channel 1 and second channel 2
            wcenter = (winlim(iw-1, 2) + winlim(iw, 1)) / 2.0
            idx = MAXVAL ( MAXLOC ( irrad_wavl(spos(ch-1):epos(ch-1),iix+1), &
                   MASK = irrad_wavl(spos(ch-1):epos(ch-1), iix+1) < wcenter ) ) + spos(ch-1) - 1
            DO i = winpix(iw-1, 2)+1, idx 
              IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                   ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
                   ALL(irrad_wavl(i, iix+1:iix+nbin) > subring(wvl_idx, nring)) ) THEN
                nring = nring + 1
                subring(wvl_idx, nring) = SUM(irrad_wavl(i, iix+1:iix+nbin)) / nbin
                subring(spc_idx, nring) = SUM(irrad_spec(i, iix+1:iix+nbin)) / nbin
             ! print * , 'a',nring, subring(wvl_idx, nring) !, subring(wvl_idx, nring)-subring(wvl_idx,nring-1)
              ENDIF
            ENDDO
            omi_ring%ndiv(ix) = nring               ! 1:nring is from the same channel
            nbin = nwbin(iw) 
            iix = (ix - 1) * nbin
            IF (ch == 2 .AND. nsolbin == 2) THEN  ! Shift the position by 15 
              iix = iix - (zoom_p1 - 1) * nsolbin
            ENDIF

            idx = MAXVAL ( MINLOC ( irrad_wavl(spos(ch):epos(ch), iix+1),   &
                   MASK = irrad_wavl(spos(ch):epos(ch), iix+1) > wcenter ) ) + spos(ch) - 1            
            DO i = idx, winpix(iw, 1) - 1 
              IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                   ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
                   ALL(irrad_wavl(i, iix+1:iix+nbin) > subring(wvl_idx, nring)) ) THEN
                nring = nring + 1
                subring(wvl_idx, nring) = SUM(irrad_wavl(i, iix+1:iix+nbin)) / nbin
                subring(spc_idx, nring) = SUM(irrad_spec(i, iix+1:iix+nbin)) / nbin
            ! print * , 'b',nring, subring(wvl_idx, nring), subring(wvl_idx,nring)-subring(wvl_idx,nring-1)
              ENDIF
            ENDDO
          ENDIF
          idx = SUM(omi_irrad%npix(1:iw-1, ix))
          subring(wvl_idx:spc_idx, nring+1:nring+omi_irrad%npix(iw, ix)) = subspec(1, wvl_idx:spc_idx, &
               idx+1:idx+omi_irrad%npix(iw, ix))
          !print * , subring(1, nring+1), subring(1, nring+omi_irrad%npix(iw,ix))
          nring = nring + omi_irrad%npix(iw, ix)
        ENDDO
        ! Add extra spectra after fitting window
        noff2 = nring
        IF (ch == 2) THEN
          nring = nring + noff_uv2
        ELSE
          nring = nring + noff_uv1
        ENDIF
        nbad = 0
        DO i = winpix(numwin, 2) + 1, nwavel
          IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
            nbin = nxbin
          ELSE
            nbin = nxbin * ncoadd
          ENDIF
          iix = (ix - 1) * nbin

          IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
               ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
            noff2 = noff2 + 1
            subring(wvl_idx, noff2) = SUM(irrad_wavl(i, iix+1:iix+nbin)) / nbin
            subring(spc_idx, noff2) = SUM(irrad_spec(i, iix+1:iix+nbin)) / nbin
            !print * , noff2, subring(wvl_idx, noff2)
          ELSE
             nbad = nbad + 1
          ENDIF
          IF (noff2 == nring) EXIT
        ENDDO
    
        omi_ring%nsol(ix) = nring
        omi_ring%winpix(ix, 1) = noff1
        omi_ring%winpix(ix, 2) = noff2
        !print * ,ix, omi_ring%winpix(15,:), nring, 'noff1, noff2'
        !DO j = 1, nring
        !  print * , j, subring(1, j), subring(1, j)-subring(1, j-1)
        !ENDDO
        !STOP 1
      ELSE
        ! Need to convole the saved solar spectra with additional slit width
        nring = omi_ring%nsol(ix)
        subring(1, 1:nring) = omi_solspec_ring(1, 1:nring, ix)
        subring(2, 1:nring) = omi_solspec_ring(2, 1:nring, ix)

        IF (which_slit == 0) THEN
          CALL gauss_uneven(subring(1, 1:nring), subring(2, 1:nring), nring, &
               nswath, redslw(inschs(1:nswath)), retlbnd(inschs(1:nswath)), retubnd(inschs(1:nswath)))
        ELSE IF (which_slit == 3) THEN
          CALL triangle_uneven(subring(1, 1:nring), subring(2, 1:nring), nring, &
               nswath, redslw(inschs(1:nswath)), retlbnd(inschs(1:nswath)), retubnd(inschs(1:nswath)))
        ELSE
          WRITE(www_lun, *) 'This type of slit convolution is not implemented!!!'
          pge_error_status = pge_errstat_error
        ENDIF
      ENDIF

      ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 15 pixels
      irefl = 0
      omi_refl%winpix(ix, 1:2) = 0
      IF (.NOT. coadd_uv2) THEN
        nbin = nxbin
      ELSE
        nbin = nxbin * ncoadd
      ENDIF
      nbin = nbin * nsolbin
      iix = (ix - 1) * nbin
      IF ( nsolbin == 2 ) THEN  ! Shift the position by 15 
        iix = iix - (zoom_p1 - 1) * nsolbin
      ENDIF

      idx = MAXVAL ( MINLOC ( irrad_wavl(1:nwavel, iix+1), MASK = &
           (irrad_wavl(1:nwavel, iix+1) > pos_alb - toms_fwhm * 1.4) ))
      DO i  = idx, nwavel
        IF (ALL(irrad_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
             ALL(irrad_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
          irefl = irefl + 1
          omi_refl%solwavl(irefl, ix) = SUM(irrad_wavl(i, iix+1:iix+nbin)) / nbin
          omi_refl%solspec(irefl, ix) = SUM(irrad_spec(i, iix+1:iix+nbin)) / nbin
          IF (omi_refl%winpix(ix, 1) == 0) omi_refl%winpix(ix, 1) = i
          omi_refl%winpix(ix, 2) = i
        ENDIF
        IF (irefl == nrefl) EXIT
      ENDDO

      IF (irefl /= nrefl) THEN
        WRITE(www_lun, *) &
             'Could not get enough irradiance points for cloud fraction!!!'
        omi_irrad%errstat(ix) = pge_errstat_error
        CYCLE
      ENDIF

      IF (scnwrt) THEN
        WRITE(www_lun, *) 'End Of Reading Irradiance Spectrum: ', ix
        DO i = 1, numwin
          WRITE(www_lun,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), &
               winlim(i,2), omi_irrad%npix(i, ix)
          IF (omi_irrad%npix(i, ix) < 4) THEN
            WRITE(www_lun, '(A,f8.3,A3,f8.3)') &
                 ' Not enough points (>=4)  in window: ', winlim(i,1), &
                 ' - ', winlim(i,2)
            pge_error_status = pge_errstat_error
          ENDIF
        ENDDO
      ENDIF

      omi_irrad%norm(ix) = SUM ( subspec(1, spc_idx, 1:nomi) ) / nomi

      IF ( omi_irrad%norm(ix) <= 0.0 ) THEN 
        omi_irrad%errstat(ix) = pge_errstat_error
        CYCLE
      ENDIF

      omi_irrad%wavl(1:nomi, ix) = &
           real(subspec(1, wvl_idx, 1:nomi) , kind=r4)
      omi_irrad%spec(1:nomi, ix) = &
           real(subspec(1, spc_idx, 1:nomi) / omi_irrad%norm(ix) , kind=r4)
      omi_irrad%prec(1:nomi, ix) = &
           real(subspec(1, sig_idx, 1:nomi) / omi_irrad%norm(ix) , kind=r4)
      !omi_irrad_stray(1:nomi, ix) = strayspec(1, 1, 1:nomi) / omi_irrad%norm(ix)
      !omi_rad_stray  (1:nomi, ix) = strayspec(1, 2, 1:nomi) / omi_irrad%norm(ix)

      omi_ring%wavl(1:nring, ix) = real(subring(1, 1:nring) , kind=r4)
      omi_ring%spec(1:nring, ix) = real(subring(2, 1:nring) / &
           omi_irrad%norm(ix) , kind=r4)
      dwavmax = (omi_irrad%wavl(2,1) - omi_irrad%wavl(1,1))*1.1
      !DO i = 1, nomi
      !   WRITE(90, '(F10.4, 2D16.7)') irrad_wavl(i, ix), &
      !        irrad_spec(i, ix) * omi_irrad%norm(ix) , irrad_prec(i, ix) * omi_irrad%norm(ix)
      !ENDDO
      !STOP 1

      ! Back up solar irradiance from OMTO3
      ! IF (orbnumsol == 99999) irrad_prec(1:nomi, ix)  = 1.0
    ENDDO
   !--------------------------------------------------------------------------
   ! Ending  with deallocating local variables
   !--------------------------------------------------------------------------
   deallocate (irrad_qflg, irrad_prec, irrad_spec, irrad_wavl, idxs)
   deallocate (flgbits)
   deallocate (tmpspec, tmpqflg, tmpnavg)
   deallocate (subspec, subring)
   RETURN
  END SUBROUTINE omi_read_irradiance_data


  SUBROUTINE omi_read_radiance_lines ( iline, ny, first_line, &
       first_pix, last_pix, nxcoadd, pge_error_status )

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
    USE OMSAO_parameters_module, ONLY: maxwin, vb_lev_omidebug, mswath
    USE OMSAO_variables_module,  ONLY: verb_thresh_lev, l1b_rad_filename, &
         wcal_bef_coadd, numwin, coadd_uv2, band_selectors, &
         szamax, currpix, reduce_resolution, redlam, redsampr, &
         reduce_slit, correct_merr, ybin_decerr, & 
         nxbin, nybin, ncoadd,redslw, nswath, inschs, &
         nxtrack, ntimes, reduce_lbnd, reduce_ubnd, retlbnd, retubnd
    USE ozprof_data_module,       ONLY:  nrefl
    USE OMSAO_pixelcorner_module, ONLY: rowanomaly_flg
    USE OMSAO_errstat_module
    USE hdfeos4_parameters
    USE L1B_Reader_class
    use m_convert_coadd, only: prespec_align, radwavcal_coadd
    use tell_module
  

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,  INTENT (IN) :: iline, ny, first_pix, last_pix, first_line, &
         nxcoadd

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER, INTENT (OUT) :: pge_error_status


    ! ---------------
    ! Local variables
    ! ---------------
    TYPE (l1b_block_type) :: omi_data_block
    INTEGER (KIND=i4)     :: blockline, errstat
    INTEGER, PARAMETER    :: nbits = 16
    INTEGER (KIND=i2), DIMENSION(0:nbits-1) :: mflgbits , tmp_mflgbits
    INTEGER (KIND=i2), DIMENSION(:,:,:), POINTER :: flgbits
    INTEGER (KIND=i2), DIMENSION(:), POINTER     :: flgmsks
    INTEGER (KIND=i4), DIMENSION(maxwin)         :: nwbin
    integer (kind=i2), dimension(1) :: temp_omi_mflg
    REAL (KIND=r4)  :: tmp_ExposureTime

    REAL (KIND= dp)                                        :: tmpNinteg

    INTEGER   (KIND=4), DIMENSION (mswath)    :: nwls
    INTEGER                                   :: nwavel, is, iloop, nwl, &
         i, j, ix, iix, ii, nomi, fidx, lidx, ch, iw, ic, irefl, nx, &
         nbin, fpix, lpix, np, npos
   ! variables used to read original spectrum
    LOGICAL                                   :: error
    INTEGER   (KIND=i2), DIMENSION(mswath)    :: spos, epos
    INTEGER, DIMENSION (:), POINTER           :: idxs
    INTEGER (kind=2), POINTER, DIMENSION (:,:,:) :: rad_qflg
    REAL (kind=i4), POINTER, DIMENSION(:,:,:) ::  rad_spec,rad_prec,rad_wavl
    REAL (KIND=I4), POINTER, DIMENSION(:,:) :: ccd_spec, ccd_prec, ccd_wavl
    INTEGER (KIND=2),POINTER, DIMENSION(:,:) :: ccd_qflg
    ! variables used for reduced resolution
    REAL (KIND = dp), DIMENSION (:,:,:), POINTER :: tmpspec
    INTEGER (KIND=2), DIMENSION (:,:), POINTER :: tmpqflg
    REAL (KIND = dp)                          :: tmpsampr, retswav, retewav
    REAL (KIND = dp), DIMENSION (:,:,:), POINTER :: subspec
    LOGICAL, DIMENSION (numwin)                  :: wavcals 
    REAL (KIND = dp), DIMENSION (numwin,nxcoadd) :: wshis, wsqus

    ! Exteranl functions
    INTEGER                                   :: estat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_lines'

    ! JCH: try to catch nwl being used without being properly initialized
    nwl = 2147483647

    nx = nfxtrack / nxbin
    ! allocation
    allocate (rad_qflg (nwavel_max, nxtrack_max, 0:ny-1))
    allocate (rad_prec (nwavel_max, nxtrack_max, 0:ny-1))
    allocate (rad_wavl (nwavel_max, nxtrack_max, 0:ny-1))
    allocate (rad_spec (nwavel_max, nxtrack_max, 0:ny-1))
    allocate (ccd_prec (nwavel_max, nxtrack_max))
    allocate (ccd_spec (nwavel_max, nxtrack_max))
    allocate (ccd_wavl (nwavel_max, nxtrack_max))
    allocate (ccd_qflg (nwavel_max, nxtrack_max))
    allocate (idxs(nwavel_max), flgmsks(nwavel_max))
    allocate (flgbits(nxcoadd, nwavel_max, 0:nbits-1))
    allocate (tmpspec(sig_idx, nwavel_max, nxtrack_max), &
              tmpqflg(nwavel_max,nxtrack_max))
    allocate (subspec(nxcoadd, sig_idx, nwavel_max) )

    ! Initialize all local data arrays
    rad_spec (1:nwavel_max, 1:nxtrack_max, 0:ny-1) = 0.0
    rad_prec (1:nwavel_max, 1:nxtrack_max, 0:ny-1) = 0.0
    rad_qflg (1:nwavel_max, 1:nxtrack_max, 0:ny-1) = 0
    rad_wavl (1:nwavel_max, 1:nxtrack_max, 0:ny-1) = 0.0

    ccd_spec (1:nwavel_max, 1:nxtrack_max) = 0.0
    ccd_prec (1:nwavel_max, 1:nxtrack_max) = 0.0
    ccd_qflg (1:nwavel_max, 1:nxtrack_max) = 0
    ccd_wavl (1:nwavel_max, 1:nxtrack_max) = 0.0


    errstat = omi_s_success
    pge_error_status = pge_errstat_ok
    omi_rad%errstat(0:ny-1)          = pge_errstat_ok
    omi_rad%pix_errstat(1:nxtrack, 0:ny-1) = pge_errstat_ok
    omi_rad%npix = 0
    omi_rad%nwav = 0
    omi_rad%spec (:, 1:nxtrack, 0:ny-1) = 0.0
    omi_rad%prec (:, 1:nxtrack, 0:ny-1) = 0.0
    omi_rad%qflg (:, 1:nxtrack, 0:ny-1) = 0
    omi_rad%wavl (:, 1:nxtrack, 0:ny-1) = 0.0

    j = 1
    nwavel = 0
    wavcals = .TRUE.

    DO is = 1, nswath
      ch = inschs(is)

      ! Open data block called 'omi_data_block' with default size of 100 lines
!      errstat = L1Br_OPEN ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(is), ny )
! FIXME - for some reason insists on having block size = full swath length 
!      errstat = L1Br_OPEN ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(is), 1644 )
        errstat = L1Br_OPEN ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(ch), ntimes) !1643)
        IF ( errstat /= omi_s_success ) THEN
          estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, 'L1Br_OPEN failed.', modulename, 0 ) 
          STOP 1
        END IF


      DO iloop = 0, ny - 1
        ! The current scan line number we are reading
        mflgbits = 0

        DO i = 0, nybin - 1
          blockline = first_line + iloop * nybin + i


          !Measurement flag for the line
            errstat = L1Br_getDATA ( omi_data_block, blockline, &
                 MeasurementQualityFlags_k = omi_mflg)
            IF( errstat /= omi_s_success ) THEN
              estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, &
                   "L1Br_getDATA failed.", modulename, 0 )
              STOP 1
            END IF

          ! Get radiances associated with wavelength range
            errstat = L1Br_getSIGline ( omi_data_block, blockline,   &
                 Signal_k            = ccd_spec(j:, :),  &
                 SignalPrecision_k   = ccd_prec(j:, :),  &
                 PixelQualityFlags_k = ccd_qflg(j:, :),  &
                 Wavelength_k        = ccd_wavl(j:, :),  &
                 Nwl_k  = nwls(ch) )
            IF( errstat /= omi_s_success ) THEN
              estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, &
                   "L1Br_getSIGline failed.", modulename, 0)
              STOP 1
            END IF
            
          temp_omi_mflg=omi_mflg
          CALL convert_2bytes_to_16bits ( nbits, 1, temp_omi_mflg, &
               tmp_mflgbits(0:nbits-1))
          mflgbits = mflgbits + tmp_mflgbits(0:nbits-1)

          IF (correct_merr) THEN
              errstat = L1Br_getDATA ( omi_data_block, blockline,   &
                   ExposureTime_k      = tmp_ExposureTime)
              IF( errstat /= omi_s_success ) THEN
                estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, &
                     "L1Br_getDATA failed.", modulename, 0)
                STOP 1
              END IF
            tmpNinteg = 2.d0 / tmp_ExposureTime
          ENDIF

          nwl = j + nwls(ch) - 1

          !print *, j, nwl, nwls(ch), ny, iloop
          !print *, ccd_wavl(j, 6), ccd_wavl(nwl, 6)
          !print *, rad_wavl(j, 6, iloop), rad_wavl(nwl, 6, iloop)

          rad_spec(j:nwl, :, iloop) = &
               rad_spec(j:nwl, :, iloop) + ccd_spec(j:nwl, :)
          IF (correct_merr) THEN
            rad_prec(j:nwl, :, iloop) = &
                 real(rad_prec(j:nwl, :, iloop) + &
                 ccd_prec(j:nwl, :) / SQRT( tmpNinteg ) , kind=r4)
          ELSE
            rad_prec(j:nwl, :, iloop) = &
                 rad_prec(j:nwl, :, iloop) + ccd_prec(j:nwl, :)
          ENDIF
          rad_wavl(j:nwl, :, iloop) = &
               rad_wavl(j:nwl, :, iloop) + ccd_wavl(j:nwl, :)

          DO ix = 1, nxtrack
            CALL coadd_2bytes_qflgs(nbits, nwls(ch), rad_qflg(j:nwl, ix, iloop), ccd_qflg(j:nwl, ix))
          ENDDO

!          IF( errstat /= omi_s_success ) THEN
!            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, &
!                 "L1Br_getSIGline failed.", modulename, 0)
!            STOP 1
!          END IF
        ENDDO



        IF (mflgbits(0) >= 1 .OR. mflgbits(1) >= 1 .OR. &
             mflgbits(3) >= 1 .OR. mflgbits(12) >= 1) THEN
          WRITE(www_lun, *) 'All radiances could not be used: line ', &
               blockline, ' Swath ', is
          omi_rad%errstat(iloop) = pge_errstat_error      
        ELSE IF (ANY(mflgbits >= 1)) THEN
          !WRITE(www_lun, *) 'Warning set on all radiances: line',  blockline, ' Swath ', is
          IF (omi_rad%errstat(iloop) /= pge_errstat_error) &
               omi_rad%errstat(iloop) = pge_errstat_warning      
          ! Over SAA region
          IF (mflgbits(10) >= 1) omi_saa_flag(iloop) = 1
        ENDIF

        rad_spec(j:nwl, :, iloop) = rad_spec(j:nwl, :, iloop) / nybin
        IF (ybin_decerr) THEN 
        rad_prec(j:nwl, :, iloop) =  real(rad_prec(j:nwl, :, iloop) / nybin / &
             SQRT(1.0D0 * nybin) , kind=r4)
        ELSE
        rad_prec(j:nwl, :, iloop) =   rad_prec(j:nwl, :, iloop) / nybin 
        ENDIF
        rad_wavl(j:nwl, :, iloop) =  rad_wavl(j:nwl, :, iloop) / nybin 
      END DO     ! end iloop

      nwavel = nwavel + nwls(ch)
      spos(ch) = int(j, kind=i2)
      j = nwavel + 1
      epos(ch) = int(nwavel, kind=i2)

      ! Close data block structure
        errstat = L1Br_CLOSE ( omi_data_block )
        IF ( errstat /= omi_s_success .AND. verb_thresh_lev >= &
             vb_lev_omidebug ) THEN
          estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, &
               'L1Br_CLOSE failed.', modulename, 0 )
        END IF

      ! Sort data in wavelength increasing order   
      IF (rad_wavl(spos(ch), 1, 0) > rad_wavl(epos(ch), 1, 0)) THEN
        idxs(spos(ch):epos(ch)) = (/ (i, i = epos(ch), spos(ch), -1) /)
        rad_wavl(spos(ch):epos(ch), :, :) = rad_wavl(idxs(spos(ch):epos(ch)), :, :)
        rad_spec(spos(ch):epos(ch), :, :) = rad_spec(idxs(spos(ch):epos(ch)), :, :)
        rad_prec(spos(ch):epos(ch), :, :) = rad_prec(idxs(spos(ch):epos(ch)), :, :)
        rad_qflg(spos(ch):epos(ch), :, :) = rad_qflg(idxs(spos(ch):epos(ch)), :, :)     
      ENDIF
    ENDDO ! end swath loop

    ! Degrade spectral resolution if necessary
    IF (reduce_resolution) THEN
      nwavel = 0
      j = 1
      DO is = 1, nswath
        ch = inschs(is)

        IF (coadd_uv2 .AND. is == 1) THEN
          npos = nfxtrack
        ELSE
          npos = nxtrack
        ENDIF
        fpix = first_pix
        lpix = last_pix
        IF (is == 2) THEN
          fpix = first_pix * 2 -1
          lpix = last_pix * 2
        ENDIF
        !IF (nswath == 2) THEN
        !   retswav = retlbnd(ch); retewav = retubnd(ch)
        !ELSE
        !   retswav = retlbnd(band_selectors(1)); retewav = retubnd(band_selectors(1))
        !ENDIF
        retswav = retlbnd(ch)
        retewav = retubnd(ch)
        tmpsampr = redsampr
        IF (is == 1 .AND. band_selectors(1) == 1) tmpsampr = redsampr / 3.0
        np = nwls(ch)

        DO iloop = 0, ny - 1
          tmpspec(wvl_idx, 1:np, fpix:lpix) = rad_wavl(spos(ch):epos(ch), fpix:lpix, iloop)
          tmpspec(spc_idx, 1:np, fpix:lpix) = rad_spec(spos(ch):epos(ch), fpix:lpix, iloop)
          tmpspec(sig_idx, 1:np, fpix:lpix) = rad_prec(spos(ch):epos(ch), fpix:lpix, iloop)
          tmpqflg(       1:np, fpix:lpix) = rad_qflg(spos(ch):epos(ch), fpix:lpix, iloop)       

          DO ix = fpix, lpix
            CALL convert_2bytes_to_16bits ( nbits, np, tmpqflg(1:np, ix), flgbits(1, 1:np, 0:nbits-1))
            tmpqflg(1:np, ix) = flgbits(1, 1:np, 0) &   ! Missing
                 + flgbits(1, 1:np, 1)                &   ! Bad
                 + flgbits(1, 1:np, 2)                &   ! Processing error
           !      + flgbits(1, 1:np, 4)                &   ! RTS_Pixel_Warning Flag
                 + flgbits(1, 1:np, 5)                &   ! Saturation Possibility Flag
                 + flgbits(1, 1:np, 7)                    ! Dark Current Warning Flag
          ENDDO

          CALL reduce_rad_resolution (tmpspec(:, 1:np, fpix:lpix), tmpqflg(1:np, fpix:lpix),   &
               np, lpix-fpix+1, reduce_slit, redslw(is), tmpsampr, redlam, retswav, retewav, reduce_lbnd(ch), &
               reduce_ubnd(ch), nwls(ch), pge_error_status)
          IF (pge_error_status == pge_errstat_error) RETURN

          nwavel = j + nwls(ch) - 1
          rad_wavl(j:nwavel, fpix:lpix, iloop) = &
               real(tmpspec(wvl_idx, 1:nwls(ch), fpix:lpix), kind=r4)
          rad_spec(j:nwavel, fpix:lpix, iloop) = &
               real(tmpspec(spc_idx, 1:nwls(ch), fpix:lpix) , kind=r4)
          rad_prec(j:nwavel, fpix:lpix, iloop) = &
               real(tmpspec(sig_idx, 1:nwls(ch), fpix:lpix) , kind=r4)
          rad_qflg(j:nwavel, fpix:lpix, iloop) = 0   ! All data are good  (pre filtered)  
        ENDDO
        spos(ch) = int(j, kind=i2)
        epos(ch) = int(nwavel, kind=i2)
        j = nwavel + 1
      ENDDO
    ENDIF

    IF (nwavel > nwavel_max) THEN
      WRITE(www_lun, *) "Need to increase nwavel_max!!!"
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF


    ! Determine number of binning for different fitting windows
    DO iw = 1, numwin
      ch = band_selectors(iw)
      IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
        nwbin(iw) = nxbin
      ELSE
        nwbin(iw) = nxbin * ncoadd
      ENDIF
    ENDDO

    ! Subset and coadd radiance spectrum
    DO iloop = 0, ny - 1 
      IF (omi_rad%errstat(iloop) == pge_errstat_error) THEN
        omi_rad%pix_errstat(first_pix:last_pix, iloop) = pge_errstat_error
        CYCLE
      ENDIF

      DO ix = first_pix, last_pix
        currpix = ix
!        IF (omi_szenith(ix, iloop) > szamax .OR. omi_szenith(ix, iloop) < 0 ) THEN
        IF (omi_geo%sza (ix, iline+iloop) > szamax .OR. omi_geo%sza (ix, iline+iloop) < 0 ) THEN
          omi_rad%pix_errstat(ix, iloop) = pge_errstat_error
          CYCLE
        ENDIF

        IF (rowanomaly_flg(ix, iline+iloop) == 1 .AND. ix /= 24) THEN
          omi_rad%pix_errstat(ix, iloop) = pge_errstat_error
          CYCLE
        ENDIF

        ! Get quality flag bits
        ! Coadd uv-2 flags if necessary to avoid coadding inconsistent # of pixels 
        flgmsks = 0
        DO is = 1, nswath
          ch = inschs(is)
          IF (is == 1) THEN
            nbin = nxbin
          ELSE
            nbin = nxbin * ncoadd
          ENDIF
          iix = (ix - 1) * nbin 

          ! properly align cross track positions to be coadded (should be within one pixel)
          IF (.NOT. reduce_resolution) THEN
            IF (nbin > 2) CALL prespec_align(nwls(ch), nbin, rad_wavl(spos(ch):epos(ch),&
                 iix+1:iix+nbin, iloop), rad_spec(spos(ch):epos(ch), iix+1:iix+nbin, iloop), &
                 rad_prec(spos(ch):epos(ch), iix+1:iix+nbin, iloop), &       
                 rad_qflg(spos(ch):epos(ch), iix+1:iix+nbin, iloop))
            DO ic = 1, nbin
              CALL convert_2bytes_to_16bits ( nbits, nwls(ch), rad_qflg(spos(ch):epos(ch), &
                   iix + ic, iloop), flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
              flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) &
                   + flgbits(ic, spos(ch):epos(ch), 0)                &   ! Missing
                   + flgbits(ic, spos(ch):epos(ch), 1)                &   ! Bad 
                   + flgbits(ic, spos(ch):epos(ch), 2)                &   ! Processing error
           !        + flgbits(ic, spos(ch):epos(ch), 4)                &   ! RTS_Pixel_Warning Flag
                   + flgbits(ic, spos(ch):epos(ch), 5)                &   ! Saturation Possibility Flag
                   + flgbits(ic, spos(ch):epos(ch), 7)                     ! Dark Current Warning Flag
             !print *, ic,flgbits(ic, 40, 0:7), rad_qflg(40, iix+ic, iloop), iloop, iline
             !if (is == 1) flgmsks(40) = 1
            ENDDO
            
            !DO i = spos(ch), epos(ch)
            !   WRITE(91, '(F10.4, D14.6, 16I2)') rad_wavl(i, iix+1, iloop), &
            !        rad_spec(i, iix+1, iloop), flgbits(1, i, 0:nbits-1)
            !ENDDO

          ELSE
            ! Already aligned because of using common wavelength scale
            DO ic = 1, nbin
              flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) + &
                   rad_qflg(spos(ch):epos(ch), iix + ic, iloop)
            ENDDO
          ENDIF
        ENDDO

        ! Subset valid data
        nomi = 0
        subspec = 0.0
        fidx = 1
        DO iw = 1, numwin
          ch = band_selectors(iw)                  
          omi_rad%npix(iw, ix, iloop) = nomi
          nbin = nwbin(iw)
          iix = (ix - 1) * nbin

          !fidx = omi_irrad%winpix(iw, ix, 1) 
          !lidx = omi_irrad%winpix(iw, ix, 2) 
          lidx = fidx + omi_irrad%npix(iw, ix) - 1

          DO ii = fidx, lidx
            i = omi_irrad%wind(ii, ix)
            IF (ALL(rad_spec(i, iix+1:iix+nbin, iloop) > 0.0) .AND. &
                 ALL(rad_spec(i, iix+1:iix+nbin, iloop) < 4.0E14) .AND. flgmsks(i) == 0 ) THEN
              nomi = nomi + 1
              subspec(1:nbin, wvl_idx, nomi) = rad_wavl(i, iix+1:iix+nbin, iloop)
              subspec(1:nbin, spc_idx, nomi) = rad_spec(i, iix+1:iix+nbin, iloop)
              subspec(1:nbin, sig_idx, nomi) = rad_prec(i, iix+1:iix+nbin, iloop)
              omi_rad%wind(nomi, ix, iloop) = int(ii, kind=i2)
            ENDIF
            !print * , ii, i, nomi, flgmsks(i), rad_wavl(i, iix+1, iloop)
          ENDDO
          fidx = lidx + 1
          omi_rad%npix(iw, ix, iloop) = nomi - omi_rad%npix(iw, ix, iloop)
          ! If the # of wavelengths is <= 75% of the # of irradiances, stop processing this pixel
          IF (omi_rad%npix(iw, ix, iloop) <= omi_irrad%npix(iw, ix) * 0.9 ) THEN
            WRITE(*, '(A,5I5,F9.2)') 'Too fewer radiance points: ', ix, iloop, iw, &
                 omi_rad%npix(iw, ix, iloop), omi_irrad%npix(iw, ix)!, omi_szenith(ix, iloop)
                 omi_rad%pix_errstat(ix, iloop) = pge_errstat_error
            STOP 1
            EXIT
          ENDIF

          !WRITE(www_lun, '(2I5, 2F8.3, 2I5, 2F8.3, 6I5)') ix, iw, rad_wavl(spos(ch), iix+1, iloop), &
          !     rad_wavl(epos(ch), iix+1, iloop), spos(ch), epos(ch),  &
          !     fidx, lidx, lidx - fidx + 1, omi_rad%npix(iw, ix, iloop)
        ENDDO

        IF (omi_rad%pix_errstat(ix, iloop) == pge_errstat_error) CYCLE   ! This pixel will not be processed.     
        omi_rad%nwav(ix, iloop) = nomi

        ! Perform coadding if UV-2 is selected with UV-1
        fidx = 1
        DO iw = 1, numwin
          ch = band_selectors(iw)
          nbin = nwbin(iw)
          lidx = fidx + omi_rad%npix(iw, ix, iloop) - 1 

          IF (nbin > 1) THEN
            CALL radwavcal_coadd(wcal_bef_coadd, wavcals(iw), & !iw, ix, &
                 omi_rad%npix(iw, ix, iloop), nbin, &
                 subspec(1:nbin, :, fidx:lidx), wshis(iw,1:nbin), &
                 wsqus(iw, 1:nbin), error)
            wavcals(iw) = .FALSE.
            IF (error) THEN
              WRITE(www_lun, '(A)') 'No radiance wavelength calibration before coadding!!!'
              pge_error_status = pge_errstat_warning
            ENDIF
          ENDIF
          fidx = lidx + 1
        ENDDO

        ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 20 pixels
        irefl = 0
        fidx = int(omi_refl%winpix(ix, 1) , kind=i4)
        IF (.NOT. coadd_uv2) THEN
          nbin = nxbin
        ELSE
          nbin = nxbin * ncoadd
        ENDIF
        iix = (ix - 1) * nbin


        DO i  = fidx, nwavel
          IF ( ALL(rad_spec(i, iix+1:iix+nbin, iloop) > 0.0) .AND. &
               ALL(rad_spec(i, iix+1:iix+nbin, iloop) < 4.0E14) .AND. flgmsks(i) == 0) THEN
            irefl = irefl + 1
            omi_refl%radwavl( irefl, ix, iloop) = SUM(rad_wavl(i, iix+1:iix+nbin, iloop)) / nbin
            omi_refl%radspec( irefl, ix, iloop) = SUM(rad_spec(i, iix+1:iix+nbin, iloop)) / nbin
          ENDIF
          IF (irefl == nrefl) EXIT
        ENDDO

 
        IF (irefl /= nrefl) THEN
          WRITE(www_lun, '(A, 2I5, F9.2)') 'Number of rad/sol points (cloud fraction) do not match: ', &
               ix, iloop !, omi_szenith(ix, iloop)
          omi_rad%pix_errstat(ix, iloop) = pge_errstat_error  
        ENDIF

        !IF (scnwrt) THEN
        !   WRITE(www_lun, *) 'End Of Reading Radiance Spectrum: ', ix, iloop + iline
        !   DO i = 1, numwin
        !      WRITE(www_lun,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i, 1), winlim(i, 2), omi_rad%npix(i, ix, iloop)
        !   
        !      IF (omi_rad%npix(i, ix, iloop) <= 20) THEN
        !         WRITE(www_lun, '(A,f8.3,A3,f8.3)') ' Not enough points in window: ', winlim(i, 1), ' - ', winlim(i, 2)
        !         pge_error_status = pge_errstat_error
        !        
        !      ENDIF
        !   ENDDO
        !   print *, subspec(1, 1, 138:139)
        !ENDIF

        omi_rad%norm(ix, iloop) = SUM ( subspec(1, spc_idx, 1:nomi) ) / nomi
        !omi_rad%norm(ix, iloop) = 1.0E11 
        IF ( omi_rad%norm(ix, iloop) <= 0.0 ) THEN 
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
        omi_rad%wavl(1:nomi, ix, iloop) = real(subspec(1, wvl_idx, 1:nomi) , kind=r4)
        omi_rad%spec(1:nomi, ix, iloop) = real(subspec(1, spc_idx, 1:nomi) /omi_rad%norm(ix, iloop) , kind=r4)
        omi_rad%prec(1:nomi, ix, iloop) = real(subspec(1, sig_idx, 1:nomi) /omi_rad%norm(ix, iloop) , kind=r4)
      ENDDO
    ENDDO

   !-------------------------------------------------
   ! finishing with deallocation
   !-------------------------------------------------
   deallocate(rad_qflg, rad_prec, rad_wavl, rad_spec)
   deallocate(ccd_qflg, ccd_prec, ccd_wavl, ccd_spec)
   deallocate(idxs, flgmsks, flgbits)
   deallocate(tmpspec, tmpqflg, subspec)
    !DO ix = 1, nfxtrack
    !   WRITE(90, '(10I5)') ix, omi_irrad%npix(1:numwin, ix), omi_irrad%nwav(ix)
    !   DO i = 0, nloop - 1
    !      WRITE(90, '(I5, F10.3, 3I5)') i, omi_szenith(ix, i), omi_rad%npix(1:numwin, ix, i), %omi_rad%nwav(ix, i)
    !   ENDDO
    !ENDDO

    RETURN
  END SUBROUTINE omi_read_radiance_lines

  ! Replace Solar Composite with original OMI solar irradiance
  SUBROUTINE replace_solar_irradiance (lun, first_pix, last_pix, &
       pge_error_status )
   
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: mswath
    !USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
    USE OMSAO_parameters_module, ONLY: maxchlen, maxwin!, mrefl
    USE OMSAO_variables_module,  ONLY: avg_solcomp, avgsol_allorb, numwin, &
         coadd_uv2, band_selectors, refdbdir, nxbin, nswath, nxtrack, ncoadd, orbnum
    USE OMSAO_omidata_module,    ONLY: nxtrack_max, nwavel_max
    !USE ozprof_data_module,      ONLY: nrefl
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,  INTENT (IN) :: first_pix, last_pix, lun

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER, INTENT (OUT) :: pge_error_status

    ! ----------------
    ! Loacal variables
    ! ----------------
    INTEGER, PARAMETER :: ntype = 2, mscpt = 7000!, norbtype = 2
    CHARACTER(LEN=3), DIMENSION(ntype) :: comp_types = (/'med', 'pc0'/)  ! average and 1st principal compoment
    CHARACTER(LEN=3), DIMENSION(mswath):: channels   = (/'uv1', 'uv2'/)  ! average and 1st principal compoment
    CHARACTER(LEN=4), DIMENSION(ntype) :: orb_types  = (/'comp','1st7'/) ! all orbits and first 7 days after OPF change

    ! For the solar composite data
    INTEGER,        DIMENSION(mswath)                      :: nscpts   
    REAL (KIND=dp), DIMENSION(mswath)                      :: snorms                
    REAL (KIND=dp), DIMENSION(:,:,:), POINTER  :: solcomp
    REAL (KIND=dp), DIMENSION(:,:), POINTER    :: solcomp_wvl

    INTEGER :: i, j, ch, nx, errstat, iw
    INTEGER (KIND=i4), DIMENSION(maxwin)  :: nwbin
    REAL (KIND=dp)                        :: swav, ewav
    CHARACTER(LEN=maxchlen)               :: scfname
    CHARACTER(LEN=3)                      :: chc, typec
    CHARACTER(LEN=4)                      :: orbtypec, opfc
    REAL (KIND=dp), DIMENSION(:), POINTER :: tmpspec, tmpwvl
    
    ! ------------------
    ! External functions
    ! ------------------
    !INTEGER :: OMI_SMF_setmsg

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=24), PARAMETER :: modulename = 'replace_solar_irradiance'

    allocate (solcomp(mswath, nxtrack_max, mscpt))
    allocate (solcomp_wvl (mswath, mscpt))
    allocate (tmpspec(nwavel_max), tmpwvl(nwavel_max))

    pge_error_status = pge_errstat_ok

    ! Selectively read the solar composite data
    DO i = 1, nswath
      IF (nswath == 1) THEN
        ch = band_selectors(1)
      ELSE
        ch = i
      ENDIF

      chc = channels(ch)
      IF (avg_solcomp) THEN 
        typec = comp_types(1)
      ELSE
        typec = comp_types(2)
      ENDIF

      IF (orbnum >= 6551) THEN
        opfc = 'ge25'
      ELSE
        opfc = 'lt25'
      ENDIF

      IF (avgsol_allorb) THEN 
        orbtypec = orb_types(1)
      ELSE
        orbtypec = orb_types(2)
      ENDIF

      scfname = ADJUSTL(TRIM(refdbdir)) // 'OMI/SolarComposite/' // chc // '_' // typec // '_' &
           // opfc // '_' // orbtypec // '.dat'
      !print *, i, ch, ADJUSTL(TRIM(scfname))

      OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(scfname)), STATUS='UNKNOWN', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Cannot open solar composite file!!!'
        pge_error_status = pge_errstat_error
        RETURN
      END IF

      READ(lun, *) nx, nscpts(ch), swav, ewav, snorms(ch)
      IF (  (ch == 1 .AND. nx /= nxtrack / ncoadd) .OR. (ch == 2 .AND. nx /= nxtrack) ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Solar composite does not cover all xtrack positions!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
      IF (nscpts(ch) > mscpt) THEN
        WRITE(www_lun, '(2A)') modulename, ': Increase the dimension mscpt for solar composite!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF

      DO j = 1, nscpts(ch)
        READ(lun, *) solcomp_wvl(ch, j), solcomp(ch, 1:nx, j)
      ENDDO
      CLOSE(lun)         
    ENDDO

    ! Determine number of binning for different fitting windows
    DO iw = 1, numwin
      ch = band_selectors(iw)
      IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
        nwbin(iw) = nxbin
      ELSE
        nwbin(iw) = nxbin * ncoadd
      ENDIF
    ENDDO

    deallocate (solcomp, solcomp_wvl, tmpspec, tmpwvl)
    RETURN
  END SUBROUTINE replace_solar_irradiance

! Correction for wavelength registration at 1:67 and 498:557
  !SUBROUTINE corruv2wav(nw, nx, waves)
  !  USE OMSAO_precision_module
  !
  !  ! ---------------
  !  ! Input variables
  !  ! ---------------
  !  INTEGER (KIND=i4),                  INTENT (IN)    :: nw, nx
  !  REAL (KIND=r4), DIMENSION (nw, nx), INTENT (INOUT) :: waves
  !
  !  ! ----------------
  !  ! Local variables
  !  ! ---------------
  !  INTEGER        :: ix
  !  REAL (KIND=r4) :: del, ndel, delp1, delp2, delm1, delm2, ndelp1, ndelm1, sh1, sh2
  !
  !  !RETURN
  !
  !  !print *, nw, nx
  !  DO ix = 1, nx   
  !    ! At position 67
  !    del    = waves(68, ix) - waves(67, ix)
  !    delm1  = waves(67, ix) - waves(66, ix)
  !    delp1  = waves(69, ix) - waves(68, ix)
  !    delp2  = waves(70, ix) - waves(69, ix)
  !    ndel   = delp1 * 2 - delp2
  !    ndelm1 = ndel * 2  - delp1
  !
  !    ! Shifts
  !    sh1 = (ndel - del)
  !    sh2 = sh1 + (ndelm1 - delm1)
  !
  !    !print *, 67, sh1, sh2
  !    waves(67, ix)   = waves(67, ix) - sh1
  !    waves(1:66, ix) = waves(1:66, ix) - sh2
  !
  !    ! At position 498
  !    delm2 = waves(496, ix) - waves(495, ix)
  !    delm1 = waves(497, ix) - waves(496, ix)
  !    del   = waves(498, ix) - waves(497, ix)
  !    delp1 = waves(499, ix) - waves(498, ix)
  !    ndel = delm1 * 2 - delm2
  !    ndelp1 = ndel * 2 - delm1
  !
  !    sh1 = (ndel - del)
  !    sh2 = sh1 + (ndelp1 - delp1)
  !    !print *, 498, sh1, sh2
  !
  !    waves(498, ix) = waves(498, ix) + sh1
  !    waves(499:nw, ix) = waves(499:nw, ix) + sh2
  !  ENDDO
  !
  !  RETURN
  !END SUBROUTINE corruv2wav

end module omi_read_l1b_data


