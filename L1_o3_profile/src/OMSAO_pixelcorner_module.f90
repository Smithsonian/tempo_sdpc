MODULE OMSAO_pixelcorner_module

  ! ========================================================================= !
  ! 1. Compute geolocation of OMI pixel corners from lat and longitude fields !
  ! 2. Compute one single effective viewing geometry for each pixel since     !
  !    current viewing geometry is only provided at pixel center              !
  ! 3. Deal with binning along and across the track                           !
  ! ========================================================================= !

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE L1B_Reader_class
  USE OMSAO_variables_module, ONLY: l1b_rad_filename, verb_thresh_lev, &
       the_utc, TAI93At0ZOfGranule, TAI93StartOfGranule, GranuleYear, &
       GranuleMonth, GranuleDay, GranuleHour, GranuleMinute, GranuleSecond, &
       GranuleJDay, nxbin, nybin, nswath, offset_line
  USE OMSAO_omidata_module,   ONLY: omi_radiance_swathname, nxtrack_max, &
       ntimes_max,zoom_mode, zoom_p1       
  USE m_angle_sat2toa, only: omi_angle_sat2toa
  USE m_convert_coadd, only: coadd_byte_qflgs, convert_2bytes_to_16bits,convert_byte_to_8bits

  IMPLICIT NONE

  ! Public parameters
  REAL (KIND=r4), DIMENSION (0:ntimes_max-1) :: omi_SpcftLat, &
       omi_SpcftLon, omi_SpcftAlt, omi_SecondsInDay
  INTEGER (KIND=i2), DIMENSION(0:ntimes_max-1)              :: omi_Mflg
  INTEGER (KIND=i1), DIMENSION(0:ntimes_max-1)              :: omi_NSPC
  INTEGER (KIND=i1),DIMENSION  (nxtrack_max, 0:ntimes_max-1) :: rowanomaly_flg
  INTEGER (KIND=i2), DIMENSION (nxtrack_max, 0:ntimes_max-1)  :: omi_land_water_flg, omi_glint_flg, omi_snow_ice_flg ! (nxtrack_max,0:ntimes_max-1) :: 
  ! ----------------
  ! Local Parameters
  ! ---------------------------------------------------------------------
  ! * Values for Pi (rad, deg) and Conversions between Degree and Radians
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi         = 3.14159265358979_r8  ! 2*ASIN(1.0_r8)
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf     = 0.5_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi      = 2.0_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi_deg     = 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf_deg =  90.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi_deg  = 360.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: deg2rad    = pi / 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: rad2deg    = 180.0_r8 / pi
  ! ---------------------------------------------------------------------
  ! * Precison for DEG <-> RAD conversion - anything less than EPS
  !   is effectively ZERO.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: eps = 1.0E-10_r8
  ! ---------------------------------------------------------------------
  REAL (KIND=r4), PARAMETER, PRIVATE :: rearth0 = 6378  ! equatorial radius
  REAL (KIND=r4), PARAMETER, PRIVATE :: minza = 0.0, maxza=90.0, &
       minaza = -360., maxaza = 360.0

  public compute_pixel_corners
  private get_sphgeoview_corners, sphergeom_intermediate, angle_minus_twopi, &
       circle_rdis,  convert_gpqualflag_info,convert_xtrackqflag_info
    !, sphergeom_baseline_comp, lonlat_to_pi, circle_dis


CONTAINS

  SUBROUTINE compute_pixel_corners ( geo, ntimes, nxtrack, nl, pge_error_status )

    USE m_utilities, only: day_of_year
    USE OMSAO_variables_module, ONLY:geo_group

    ! =======================================================
    ! Computes OMI pixel corner coordinates, start to finish:
    !
    !  * Reads L1b geolocation data
    !  * Computes the corners
    !  * Writes output to file
    ! =======================================================

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
    INTEGER, INTENT (IN)           :: nl
    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: pge_error_status
    TYPE(geo_group), INTENT(OUT) :: geo
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: errstat, iline, i, j, nline, nx, sline, &
         eline, ix, nxtrack1, xoff, ysidx, yeidx, ymidx, xsidx, xeidx, &
         xmidx, iy, estat, nbits, ndim!, ny
    ! * origianl dataset at each line
    REAL (KIND=r4),    DIMENSION (:), POINTER :: ccd_lat, ccd_lon, ccd_sza, &
                                                 ccd_vza, ccd_saza, ccd_vaza
    INTEGER (KIND=i2), DIMENSION (:), POINTER :: ccd_height, ccd_geoflg
    INTEGER (KIND=i1), DIMENSION (nxtrack_max):: ccd_xqflg, ccd_xqflg1
    INTEGER (KIND=i2)                         :: ccd_mflg
    INTEGER (KIND=i1)                         :: ccd_NSPC
    REAL (KIND=r4)                            :: ccd_auralon, ccd_auralat, &
                                                 ccd_auraalt, ccd_SecondsInDay  
    REAL (KIND=r8)                            :: ccd_time
    TYPE (L1B_block_type)                     :: omi_data_block
    LOGICAL :: ascendQ 

    INTEGER (KIND=i1), DIMENSION (nxtrack_max) :: &
                  waveshift_flg, blockage_flg, straysun_flg, strayearth_flg
    INTEGER (KIND=i2), DIMENSION (nxtrack_max) ::land_water_flg, glint_flg, snow_ice_flg

    ! * processing variables
    REAL (KIND=r8), DIMENSION (:,:),POINTER :: sza, vza, saza, vaza
    REAL (KIND=r8),    DIMENSION (:,:), POINTER :: & 
           omi_lat,   omi_lon, omi_sza, omi_vza, omi_aza, omi_sca
    REAL (KIND=r8),    DIMENSION (:,:), POINTER :: omi_elat, omi_elon
    REAL (KIND=r8),    DIMENSION (:,:), POINTER :: omi_clat, omi_clon
    INTEGER (KIND=i2), DIMENSION (:,:), POINTER :: omi_GeoFlg
    INTEGER (KIND=i2), DIMENSION (:,:), POINTER :: omi_Height
    INTEGER (KIND=i1), DIMENSION (:,:), POINTER :: omi_XTrackQFlg
    REAL (KIND=r8),    DIMENSION (:),   POINTER :: omi_time
    ! Exteranl functions
    INTEGER                       :: OMI_SMF_setmsg
    INTEGER (KIND=i4), EXTERNAL   :: PGS_TD_TAItoUTC, PGS_TD_UTCtoTAI !, day_of_year

    CHARACTER(LEN = 28 ) :: GranuleDAY0Z
    CHARACTER(LEN = 21)  :: modulename = 'compute_pixel_corners'

    pge_error_status = pge_errstat_ok
    errstat          = omi_s_success
    !---------------------------------------------------
    ! allocate
    !---------------------------------------------------
    ! * ccd read variables
    allocate (ccd_sza(nxtrack_max) , ccd_vza(nxtrack_max), &
              ccd_saza(nxtrack_max), ccd_vaza(nxtrack_max))
    allocate (ccd_lon(nxtrack_max),  ccd_lat(nxtrack_max),&
              ccd_height(nxtrack_max), ccd_geoflg(nxtrack_max))
   ! * processing variables
    allocate(sza(nxtrack_max, 0:ntimes_max-1), saza(nxtrack_max,0:ntimes_max-1))
    allocate(vza(nxtrack_max, 0:ntimes_max-1), vaza(nxtrack_max,0:ntimes_max-1))
    allocate(omi_time(0:ntimes_max-1))
    allocate(omi_lat(nxtrack_max, 0:ntimes_max-1), omi_lon(nxtrack_max,0:ntimes_max-1))
    allocate(omi_sza(nxtrack_max, 0:ntimes_max-1), omi_vza(nxtrack_max,0:ntimes_max-1))
    allocate(omi_aza(nxtrack_max, 0:ntimes_max-1), omi_sca(nxtrack_max,0:ntimes_max-1))
    allocate(omi_elat(0:nxtrack_max, 0:ntimes_max-1),omi_elon(0:nxtrack_max, 0:ntimes_max-1))
    allocate(omi_clat(0:nxtrack_max, 0:ntimes_max), omi_clon(0:nxtrack_max, 0:ntimes_max)) 
    allocate(omi_GeoFlg(nxtrack_max, 0:ntimes_max-1), omi_Height(nxtrack_max,0:ntimes_max-1))
    allocate(omi_XTrackQFlg(nxtrack_max, 0:ntimes_max-1))

   !  omi_Mflg(0:ntimes_max-1), omi_NSPC(0:ntimes_max-1), omi_time(0:ntimes_max-1), &
   !  omi_rowanomaly_flg(nxtrack_max,0:ntimes_max-1), omi_land_water_flg(nxtrack_max,0:ntimes_max-1), &
   !  omi_glint_flg (nxtrack_max,0:ntimes_max-1), omi_snow_ice_flg(nxtrack_max,0:ntimes_max-1), &
   !  omi_SpcftLat(0:ntimes_max-1), omi_SpcftLon(0:ntimes_max-1), omi_SpcftAlt(0:ntimes_max-1), omi_SecondsInDay(0:ntimes_max-1))

    ! -----------------------------------------------------------
    ! Open data block (UV-1, if both are selected) 
    ! called 'omi_data_block' with nTimes lines
    ! -----------------------------------------------------------
    ! ** all geolocation variables is arranged in all domain (not subset) 
    !    to be saved to L2
    errstat = L1Br_open ( omi_data_block, l1b_rad_filename, &
         omi_radiance_swathname(1))
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
           "L1Br_open failed.", modulename, 0 )
      STOP 1
    END IF

    sline = offset_line
    eline = offset_line + nl * nybin - 1
    IF (eline == sline) eline = sline + 1
    nline = eline - sline + 1

    IF (zoom_mode .AND. nswath == 1) THEN
      nxtrack1 = nxtrack / 2
      xoff = (zoom_p1 - 1) / nxbin
    ELSE
      nxtrack1 = nxtrack
      xoff = 0
    ENDIF

    errstat = L1Br_getGEOline ( omi_data_block, 0, Time_k = TAI93StartOfGranule)
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
           "L1Br_getGEOline failed.", modulename, 0 )
      STOP 1
    ENDIF
    estat = PGS_TD_TAItoUTC(TAI93StartOfGranule, the_utc)
    READ (the_utc, '(I4, 1x, I2, 1x, I2, 1x, I2, 1x, I2, 1x, F9.6)') GranuleYear, &
         GranuleMonth, GranuleDay, GranuleHour, GranuleMinute, GranuleSecond
    WRITE( GranuleDAY0Z, FMT = '(I4.4,A1,I2.2,A1,I2.2,A)' ) &
         GranuleYear, '-', GranuleMonth, '-', GranuleDay, 'T00:00:00.000Z'
    estat = PGS_TD_UTCtoTAI( GranuleDAY0Z, TAI93At0zOfGranule )
    GranuleJDay = day_of_year(GranuleYear, GranuleMonth, GranuleDay)

    ! ---------------------------------
    ! Read all Latitudes and Longitudes
    ! ---------------------------------
    DO iline = sline, eline

      errstat = L1Br_getGEOline ( omi_data_block, iline,           &
           Time_k                      = ccd_time,          &
           SecondsInDay_k              = ccd_SecondsInDay,         &
           SpacecraftAltitude_k        = ccd_auraalt,              &
           SpacecraftLongitude_k       = ccd_auralon,              &
           SpacecraftLatitude_k        = ccd_auralat,              &
           TerrainHeight_k             = ccd_height(1:nxtrack),    &
           GroundPixelQualityFlags_k   = ccd_geoflg(1:nxtrack),    &
           XTrackQualityFlags_k        = ccd_xqflg(1:nxtrack),&
           Latitude_k                  = ccd_lat (1:nxtrack),& 
           Longitude_k                 = ccd_lon (1:nxtrack),&             
           SolarZenithANgle_k          = ccd_sza (1:nxtrack),      & 
           SolarAzimuthAngle_k         = ccd_saza(1:nxtrack),      & 
           ViewingZenithAngle_k        = ccd_vza (1:nxtrack),      & 
           ViewingAzimuthAngle_k       = ccd_vaza(1:nxtrack))  
      
      omi_time(iline)               = ccd_time
      omi_SecondsInDay(iline)       = ccd_SecondsInDay    
      omi_SpcftAlt(iline)           = ccd_auraalt
      omi_SpcftLon(iline)           = ccd_auralon
      omi_SpcftLat(iline)           = ccd_auralat
      omi_Height(1:nxtrack1, iline) = ccd_height(1:nxtrack1)  
      omi_GeoFlg(1:nxtrack1, iline) = ccd_geoflg(1:nxtrack1) 
      omi_XTrackQFlg(1:nxtrack1, iline) = ccd_xqflg(1:nxtrack1)   
      omi_lat(1:nxtrack1, iline)    = ccd_lat(1:nxtrack1)   
      omi_lon(1:nxtrack1, iline)    = ccd_lon(1:nxtrack1)   
      sza (1:nxtrack1, iline)       = ccd_sza(1:nxtrack1)   
      saza(1:nxtrack1, iline)       = ccd_saza(1:nxtrack1)  
      vza (1:nxtrack1, iline)       = ccd_vza(1:nxtrack1)   
      vaza(1:nxtrack1, iline)       = ccd_vaza(1:nxtrack1) 

      IF( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
             "L1Br_getGEOline failed.", modulename, 0 )
        STOP 1
      ENDIF
    ENDDO

    ! --------------------------
    ! Close data block structure
    ! --------------------------
    errstat = L1Br_close ( omi_data_block )
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
           "L1Br_close failed.", modulename, 0 )
      STOP 1
    ENDIF

    !----------------------------------------------------------------
    ! measurement quality flags and NumberSmallPixelColumns 
    !  ==> only in UV2 swath
    ! XtrackQuality flags
    !  ==> only in UV2 because UV1 and UV2 are different
    !      and more pixels are filtered in UV2
    !errstat = L1Br_open ( omi_data_block, l1b_rad_filename, 'Earth UV-2 Swath')
    ! FIXME
    ! for some reason code fails if we don't set block size to full swath length,
    ! seems to be unable to access any line beyond block limit of 100
    !    errstat = L1Br_open ( omi_data_block, l1b_rad_filename, 'Earth UV-2 Swath', 1644)
    errstat = L1Br_open ( omi_data_block, l1b_rad_filename, 'Earth UV-2 Swath', 100) !1643
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", &
           modulename, 0 )
      STOP 1
    ENDIF

    ascendQ = .TRUE.
    IF( eline > sline ) THEN
      IF( omi_SpcftLat(sline) <  omi_SpcftLat(sline+1) ) THEN
        ascendQ = .TRUE.
      ELSE
        ascendQ = .FALSE.
      ENDIF
    ENDIF

    nbits = 8
    ndim = 1
    DO iline = sline, eline
      ccd_mflg = 0
      errstat = L1Br_getDATA ( omi_data_block, iline, &
           MeasurementQualityFlags_k = ccd_mflg)
      IF( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
             "L1Br_getData failed.", modulename, 0 )
        STOP 1
      END IF
      IF( iline > sline ) THEN
        IF( omi_SpcftLat(iline-1) <  omi_SpcftLat(iline) ) THEN
          ascendQ = .TRUE.
        ELSE
          ascendQ = .FALSE.
        ENDIF
      ENDIF
      IF( .NOT. ascendQ  ) THEN
        omi_Mflg(iline) = IBSET( ccd_mflg, 7  ) !! set bit 7 to 1 for decending pixel
      ELSE
        omi_Mflg(iline) = ccd_mflg
      ENDIF

      errstat = L1Br_getSIGline ( omi_data_block, iline, &
           NumberSmallPixelColumns_k = ccd_NSPC)
      IF( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
             "L1Br_getSIGline failed.", modulename, 0 )
        STOP 1
      END IF
      omi_NSPC(iline) = ccd_NSPC

      IF (nswath == 2) THEN
        errstat = L1Br_getGEOline ( omi_data_block, iline,           &
             XTrackQualityFlags_k = ccd_xqflg1(1:nxtrack_max))

        ! Coadd the flags
        DO ix = 1, nxtrack1
          i = ix * 2 - 1
          j = i + 1
          !CALL coadd_byte_qflgs(nbits, ndim, ccd_xqflg1(i))!, & changed by someson TEMPO team
          !     !ccd_xqflg1(j))
          CALL coadd_byte_qflgs(nbits, ndim, ccd_xqflg1(i), ccd_xqflg1(j)) ! returned by jbak
          omi_XTrackQFlg(ix, iline) = ccd_xqflg1(i)
        ENDDO
      ENDIF

    ENDDO
    errstat = L1Br_close ( omi_data_block )
    IF( errstat /= omi_s_success ) THEN
      estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
           "L1Br_close failed.", modulename, 0 )
      STOP 1
    ENDIF

    WHERE (omi_XTrackQFlg(1:nxtrack1, sline:eline) == -127)
      omi_XTrackQFlg(1:nxtrack1, sline:eline) = 0
    ENDWHERE


    !*******************************************************
    ! co-adding, convert flag
    !*****************************************************
    ! Compute the corner coordinates/viewing geometry for
    !    the spatially coadded pixels
    ! -----------------------------------------------------
    CALL get_sphgeoview_corners (nxtrack1, nline,  &
         omi_lon(1:nxtrack1, sline:eline), omi_lat(1:nxtrack1, sline:eline),       &
         sza(1:nxtrack1, sline:eline), saza(1:nxtrack1, sline:eline),                    &
         vza(1:nxtrack1, sline:eline), vaza(1:nxtrack1, sline:eline), &
         omi_clon(0:nxtrack1, sline:eline+1), omi_clat(0:nxtrack1, sline:eline+1), &
         omi_elon(0:nxtrack1, sline:eline),   omi_elat(0:nxtrack1, sline:eline),   &
         omi_sza(1:nxtrack1, sline:eline),    omi_vza(1:nxtrack1, sline:eline),    &
         omi_aza(1:nxtrack1, sline:eline),    omi_sca(1:nxtrack1, sline:eline))

    ! Re-store the data
    nx = nxtrack1 / nxbin 
    i = sline
    j = i + nl - 1

    omi_lon (1+xoff:nx+xoff, 0:nl-1)  = omi_lon (1:nx, i:j)
    omi_lat (1+xoff:nx+xoff, 0:nl-1)  = omi_lat (1:nx, i:j)
    !sza        (1+xoff:nx+xoff, 0:nl-1)  = sza        (1:nx, i:j)
    !saza       (1+xoff:nx+xoff, 0:nl-1)  = saza       (1:nx, i:j)
    !vza        (1+xoff:nx+xoff, 0:nl-1)  = vza        (1:nx, i:j)
    !vaza       (1+xoff:nx+xoff, 0:nl-1)  = vaza       (1:nx, i:j)
    omi_clon(0+xoff:nx+xoff, 0:nl)    = omi_clon(0:nx, i:j+1)
    omi_clat(0+xoff:nx+xoff, 0:nl)    = omi_clat(0:nx, i:j+1)
    omi_elon(0+xoff:nx+xoff, 0:nl-1)  = omi_elon(0:nx, i:j)
    omi_elat(0+xoff:nx+xoff, 0:nl-1)  = omi_elat(0:nx, i:j)
    omi_sza (1+xoff:nx+xoff, 0:nl-1)  = omi_sza (1:nx, i:j)
    omi_vza (1+xoff:nx+xoff, 0:nl-1)  = omi_vza (1:nx, i:j)
    omi_aza (1+xoff:nx+xoff, 0:nl-1)  = omi_aza (1:nx, i:j)
    omi_sca (1+xoff:nx+xoff, 0:nl-1)  = omi_sca (1:nx, i:j)

    ! Need to get Time, SecondsInDay, Spacecfraft altitude/latitude/longitude, 
    ! GroundPixelQualityFlags, and Terrain Height for the spatially coadded pixels
    ! GroundPixelQualityFlags, XTrackQualityFlags, and Terrain Height for the spatially coadded pixels
    ! Derive Row Anomaly Related Flags



    DO iy = 0, nl - 1 
      ysidx = sline + iy * nybin 
      yeidx = ysidx + nybin - 1
      ymidx = ysidx + nybin / 2
      omi_time(iy)         = SUM(omi_time(ysidx:yeidx))      / nybin
      omi_SecondsInDay(iy) = SUM(omi_SecondsInDay(ysidx:yeidx)) / nybin
      omi_SpcftAlt(iy)     = SUM(omi_SpcftAlt(ysidx:yeidx))  / nybin

      ! Use those from the middle point (avoid dealing with polar, dateline regions)
      omi_SpcftLat(iy)     = omi_SpcftLat(ymidx)
      omi_SpcftLon(iy)     = omi_SpcftLon(ymidx)
      omi_Mflg(iy)         = omi_Mflg(ymidx)
      omi_NSPC(iy)         = omi_NSPC(ymidx)       

      DO ix = 1, nx 
        xsidx = (ix - 1) * nxbin + 1
        xeidx = xsidx + nxbin - 1
        xmidx = xsidx + nxbin / 2
        omi_Height(ix, iy)  = INT( &
             SUM(1.0 * omi_Height(xsidx:xeidx, ysidx:yeidx)) &
             / (1.0 * nxbin * nybin), kind=i2)
        omi_GeoFlg(ix, iy) = omi_GeoFlg(xmidx, ymidx)
        omi_XTrackQFlg(ix, iy) = omi_XTrackQFlg(xmidx, ymidx)
      ENDDO

      CALL convert_xtrackqflag_info (nx, omi_XTrackQFlg(1:nx, iy), &
           rowanomaly_flg(1:nx, iy), waveshift_flg(1:nx), &
           blockage_flg(1:nx), straysun_flg(1:nx), &
           strayearth_flg(1:nx) )  
      ! get separate land/water, glint, snow/ice flags
      CALL convert_gpqualflag_info (nx, omi_GeoFlg(1:nx, iy), &
           omi_land_water_flg(1:nx, iy), omi_glint_flg(1:nx, iy), &
           omi_snow_ice_flg(1:nx, iy)) 
    ENDDO

    IF (xoff > 0) THEN
      omi_GeoFlg(1+xoff:nx+xoff, 0:nl-1) = omi_GeoFlg(1:nx, 0:nl-1)
      omi_XTrackQFlg(1+xoff:nx+xoff, 0:nl-1) = omi_XTrackQFlg(1:nx, 0:nl-1)
      omi_Height(1+xoff:nx+xoff, 0:nl-1) = omi_Height(1:nx, 0:nl-1)
    ENDIF

    i = 0
    j = nl -1
    ! MOVE TO type variables
    geo%time(i:j) = omi_time(i:j)
    geo%height(1:nx, i:j) = omi_height(1:nx,i:j)
    ! angles
    geo%sza(1:nx, i:j) = omi_sza(1:nx, i:j)
    geo%vza(1:nx, i:j) = omi_vza(1:nx, i:j)
    geo%aza(1:nx, i:j) = omi_aza(1:nx, i:j)
    geo%sca(1:nx, i:j) = omi_sca(1:nx, i:j)
    ! lon/lat
    geo%lon(1:nx, i:j) = omi_lon(1:nx, i:j)
    geo%lat(1:nx, i:j) = omi_lat(1:nx, i:j)
    geo%elon(0:nx,i:j) = omi_elon(0:nx, i:j)
    geo%elat(0:nx,i:j) = omi_elat(0:nx, i:j)
    geo%clon(1, 1:nx, i:j) = omi_clon(0:nx-1,i:j)
    geo%clon(2, 1:nx, i:j) = omi_clon(0:nx-1,i+1:j+1)
    geo%clon(3, 1:nx, i:j) = omi_clon(1:nx,i+1:j+1)
    geo%clon(4, 1:nx, i:j) = omi_clon(1:nx,i:j)
    geo%clat(1, 1:nx, i:j) = omi_clat(0:nx-1,i:j)
    geo%clat(2, 1:nx, i:j) = omi_clat(0:nx-1,i+1:j+1)
    geo%clat(3, 1:nx, i:j) = omi_clat(1:nx,i+1:j+1)
    geo%clat(4, 1:nx, i:j) = omi_clat(1:nx,i:j)


    ! flag
    geo%xflg(1:nx, i:j) = omi_XTrackQFlg(1:nx,i:j)
    geo%gflg(1:nx, i:j) = omi_geoflg(1:nx,i:j)
    geo%glint_flg(1:nx, i:j) = omi_glint_flg(1:nx,i:j)
    geo%snow_ice_flg(1:nx,i:j) = omi_snow_ice_flg(1:nx, i:j)
    geo%land_water_flg(1:nx,i:j) = omi_land_water_flg(1:nx,i:j)
    !---------------------------------------------------------------
    ! de-allocate
    !--------------------------------------------------------------
    deallocate (ccd_sza, ccd_vza, ccd_saza, ccd_vaza)
    deallocate (ccd_lon, ccd_lat, ccd_height, ccd_geoflg)
    !deallocate (ccd_xqflg, ccd_xqflg2))
    deallocate (sza, vza, saza, vaza)
    deallocate (omi_time,omi_lat, omi_lon)
    deallocate (omi_sza, omi_vza, omi_aza, omi_sca)
    deallocate (omi_elat, omi_elon, omi_clat, omi_clon)
    deallocate (omi_GeoFlg, omi_Height)
    deallocate (omi_XTrackQFlg)
     !omi_Mflg, omi_NSPC
     !omi_rowanomaly_flg, & 
     !omi_land_water_flg, &
     !omi_glint_flg , omi_snow_ice_flg, &
     !omi_SpcftLat, omi_SpcftLon, omi_SpcftAlt, omi_SecondsInDay)
    RETURN
  END SUBROUTINE compute_pixel_corners

  SUBROUTINE convert_gpqualflag_info ( &
       nxtrack, omi_geoflg, land_water_flg, glint_flg, snow_ice_flg )

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (IN) :: omi_geoflg
!    INTEGER (KIND=i2), DIMENSION (:), INTENT (IN) :: omi_geoflg

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg
!    INTEGER (KIND=i2), DIMENSION (:), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4),                PARAMETER      :: nbyte = 16
    INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = int((/ 1, 2, 4, 8, 16, 32, 64 /), kind=i2)
    INTEGER (KIND=i4)                                :: i
    INTEGER (KIND=i2), DIMENSION (nxtrack)           :: tmp_flg
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    land_water_flg = 0 
    glint_flg = 0 
    snow_ice_flg = 0

    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    tmp_flg(1:nxtrack) = int(omi_geoflg(1:nxtrack), kind=2)  ;  tmp_bytes = 0
    ! FIXME - TEMPO ground_pixel_flag is now int4, but all subroutines
    ! in this module assume int2

    CALL convert_2bytes_to_16bits ( &
         nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1) )

    ! ------------------------------
    ! The Glint flag is easy: Byte 4
    ! ------------------------------
    glint_flg(1:nxtrack) = tmp_bytes(1:nxtrack,4)

    ! ------------------------------------------------------------------
    ! Land/Water and Ice require a bit more work. The BIT slices must be
    ! multiplied with the corresponding powers of 2. The sum over this
    ! product is the information we seek.
    ! ------------------------------------------------------------------
    DO i = 1, nxtrack
      land_water_flg(i) = SUM(tmp_bytes(i,0:3 )*seven_byte(1:4))
      snow_ice_flg  (i) = SUM(tmp_bytes(i,8:14)*seven_byte(1:7))
    END DO

    RETURN
  END SUBROUTINE convert_gpqualflag_info

  SUBROUTINE convert_xtrackqflag_info ( nxtrack, omi_xtrackqflg, &
       rowanomaly_flg, waveshift_flg, blockage_flg, straysun_flg, strayearth_flg )

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
    INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (IN) :: omi_xtrackqflg

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (OUT) :: rowanomaly_flg, &
         waveshift_flg, blockage_flg, straysun_flg, strayearth_flg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4),                PARAMETER      :: nbyte = 8
    INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = int((/ 1, 2, 4, 8, 16, 32, 64 /), kind=i2)
    INTEGER (KIND=i4)                                :: i
    INTEGER (KIND=i1), DIMENSION (nxtrack)           :: tmp_flg
    INTEGER (KIND=i1), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    rowanomaly_flg = 0; waveshift_flg = 0; blockage_flg = 0
    straysun_flg = 0; strayearth_flg = 0

    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    tmp_flg(1:nxtrack) = omi_xtrackqflg(1:nxtrack)  ;  tmp_bytes = 0

    CALL convert_byte_to_8bits (nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1))

    waveshift_flg(1:nxtrack)  = tmp_bytes(1:nxtrack,4)
    blockage_flg(1:nxtrack)   = tmp_bytes(1:nxtrack,5)
    straysun_flg(1:nxtrack)   = tmp_bytes(1:nxtrack,6)
    strayearth_flg(1:nxtrack) = tmp_bytes(1:nxtrack,7)

    ! ------------------------------------------------------------------
    ! Row anomaly require a bit more work. The BIT slices must be
    ! multiplied with the corresponding powers of 2. The sum over this
    ! product is the information we seek.
    ! ------------------------------------------------------------------
    DO i = 1, nxtrack
      rowanomaly_flg(i) = int(SUM(tmp_bytes(i,0:2)*seven_byte(1:3)), kind=i1)
    ENDDO

    RETURN
  END SUBROUTINE convert_xtrackqflag_info

  !   Unused?
  !
  !  SUBROUTINE sphergeom_baseline_comp ( a0, b0, gam0, c0 )
  !    ! -------------------------------------------------------
  !    ! Finds the lengh of the baseline of a spherical triangle
  !    ! -------------------------------------------------------
  !    IMPLICIT NONE
  !
  !    ! ---------------
  !    ! Input variables
  !    ! ---------------
  !    REAL (KIND=r8), INTENT (IN) :: a0, b0, gam0
  !
  !    ! ---------------
  !    ! Output variable
  !    ! ---------------
  !    REAL (KIND=r8), INTENT (OUT) :: c0
  !
  !    ! --------------
  !    ! Local variable
  !    ! --------------
  !    REAL (KIND=r8) :: tmp
  !
  !    ! --------------------------------------------------
  !    ! Initialize output variables to keep compiler happy
  !    ! --------------------------------------------------
  !    c0 = 0.0_r8
  !
  !    ! -------------------------------------------------
  !    ! Compute length of baseline between the two points
  !    ! -------------------------------------------------
  !    tmp = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
  !    IF ( ABS(tmp) < eps ) THEN
  !      c0 = pihalf
  !    ELSE
  !      c0 = ACOS(tmp)
  !    END IF
  !
  !    RETURN
  !  END SUBROUTINE sphergeom_baseline_comp




  !   Unused?
  !
  !  SUBROUTINE lonlat_to_pi ( lon, lat )
  !
  !    IMPLICIT NONE
  !
  !    ! ------------------
  !    ! Modified variables
  !    ! ------------------
  !    REAL (KIND=r8), INTENT (INOUT) :: lon, lat
  !
  !    ! ------------------------------------
  !    ! Adjust longitude values to [-pi,+pi]
  !    ! ------------------------------------
  !    IF ( ABS(lon) > twopi ) lon = MOD(lon, twopi)
  !    IF ( lon >  pi ) lon = lon - twopi
  !    IF ( lon < -pi ) lon = lon + twopi
  !    ! ---------------------------------------
  !    ! Adjust latitude values to [-pi/2,+pi/2]
  !    ! ---------------------------------------
  !    IF ( ABS(lat) > pihalf ) lat = MOD(lat, pihalf)
  !    IF ( lat >  pihalf ) lat =   pi - lat
  !    IF ( lat < -pihalf ) lat = -(pi + lat)
  !
  !    RETURN
  !  END SUBROUTINE lonlat_to_pi

  REAL (KIND=r8) FUNCTION angle_minus_twopi ( gamma0, pival ) RESULT ( gamma )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: gamma0, pival

    IF ( gamma0 > pival ) THEN
      gamma = gamma0 - 2.0_r8 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    ELSE IF ( gamma0 < -pival ) THEN
      gamma = gamma0 + 2.0_r8 * pival 
    ELSE
      gamma = gamma0
    END IF

    RETURN
  END FUNCTION angle_minus_twopi

  SUBROUTINE get_sphgeoview_corners (nxtrack, ntimes, lon, lat, sza, saza, &
       vza, vaza, clon, clat, elon, elat, esza, evza, eaza, esca)

    use m_ezspline_interpolation, only: interpol

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                    INTENT(IN)    :: nxtrack, ntimes
    REAL    (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT(INOUT) :: lon, lat, sza, saza, vza, vaza
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes),   INTENT(OUT)   :: clon, clat
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes-1), INTENT(OUT)   :: elon, elat
    REAL    (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT(OUT)   :: esza, evza, eaza, esca     

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                   :: i, j, jj, ix, mpix, nx, ny
    INTEGER                                             :: errstat
    !REAL    (KIND=r8)                                  :: lat1, lat2, phi1, phi2, dis
    !REAL    (KIND=r8)                                  :: a0, b0, c0, gam0, alp0, a, gam
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: omixsize
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes-1):: edsza, edsazm, edvza, edvazm
    REAL    (KIND=r8), DIMENSION (1:nxtrack)            :: tmpdisx, xsize, tmpxmid
    REAL    (KIND=r8), DIMENSION (0:nxtrack)            :: tmpx, tmpsza, tmpvza, tmpsaza, tmpvaza
    REAL    (KIND=r8), DIMENSION (0:ntimes-1)           :: tmpdisy
    REAL    (KIND=dp), DIMENSION (3)                    :: zen0, zen, sazm, vazm, relaza

    ! ------------------------------------------------------------------
    ! Convert geolocation to radians; do everything in R8 rather than R4
    ! ------------------------------------------------------------------
    lon = lon * deg2rad
    lat = lat * deg2rad

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    clon   = -999.9d0
    clat = -999.9d0
    elon   = -999.9d0 
    elat = -999.9d0
    esza   = -999.9d0 
    evza = -999.9d0
    eaza = -999.9d0
    esca = -999.9d0

    ! Perform interpolation across the track
    DO i = 0, ntimes - 1
      ! Compute the distances between two pixels: (x1 + x2) / 2.
      DO ix = 1, nxtrack - 1 
        tmpdisx(ix) = circle_rdis(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i))
      ENDDO

      ! Compute the pixel size across the track
      ! Assume the center two pixels have equal pixel size (which causes about < 0.1 km error for UV-2)
      mpix = nxtrack / 2
      xsize(mpix) = tmpdisx(mpix) / 2.0
      xsize(mpix + 1) = xsize(mpix)

      DO ix = mpix -1, 1, -1
        xsize(ix) = tmpdisx(ix) - xsize(ix + 1)
      ENDDO
      DO ix = mpix + 2, nxtrack
        xsize(ix) = tmpdisx(ix - 1) - xsize(ix - 1)
      ENDDO
      omixsize(:, i) = xsize

      !!  This is to test SUBROUTINE sphergeom_intermediate
      !!  Works for both interpolation and extrapolation (with certain limitation) 
      !lat1 = -88.0 * deg2rad; lat2 = -88.0 * deg2rad
      !phi1 = 0.0  * deg2rad; phi2 = 180.0 * deg2rad 
      !dis  = circle_rdis(lat1, phi1, lat2, phi2)     
      !print *, dis
      !CALL sphergeom_intermediate ( lat1, phi1, lat2, phi2, dis, dis*0.8, a, gam )
      !WRITE(www_lun, '(6D14.6)')  lat1, phi1, a, gam, lat2, phi2
      !WRITE(www_lun, *) circle_rdis(lat1, phi1, a, gam)/dis
      !WRITE(www_lun, *) circle_rdis(a, gam, lat2, phi2)/dis

      ! Perform interpolation 
      DO ix = 1, nxtrack - 1          
        CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i), &
             tmpdisx(ix), xsize(ix), elat(ix, i), elon(ix, i))
      ENDDO
      ix = 1
      CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i),   &
           tmpdisx(ix), -xsize(ix), elat(ix-1, i), elon(ix-1, i))
      ix = nxtrack
      CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix-1, i), lon(ix-1, i),   &
           tmpdisx(ix-1), -xsize(ix), elat(ix, i), elon(ix, i))   

      ! Compute viewing geometry for west and east edge
      ! Performal interpolation/extrapolation (2 points) along the spherical lines (good enough)
      tmpx = 0.0
      DO ix = 1, nxtrack
        tmpx(ix) = tmpx(ix-1) + xsize(ix)
      ENDDO
      tmpxmid(1:nxtrack) = (tmpx(0:nxtrack-1) + tmpx(1:nxtrack)) / 2.0
      CALL interpol(tmpxmid(1:nxtrack), sza(1:nxtrack, i),  nxtrack, tmpx(0:nxtrack), &
           tmpsza(0:nxtrack),  nxtrack+1, errstat)      
      CALL interpol(tmpxmid(1:nxtrack), saza(1:nxtrack, i), nxtrack, tmpx(0:nxtrack), &
           tmpsaza(0:nxtrack), nxtrack+1, errstat)
      CALL interpol(tmpxmid(1:nxtrack), vza(1:nxtrack, i),  nxtrack, tmpx(0:nxtrack), &
           tmpvza(0:nxtrack),  nxtrack+1, errstat)
      CALL interpol(tmpxmid(1:nxtrack), vaza(1:nxtrack, i), nxtrack, tmpx(0:nxtrack), &
           tmpvaza(0:nxtrack), nxtrack+1, errstat)

      ! Check center pixel
      DO ix = mpix-1, mpix + 1
        IF (tmpvaza(ix) < 0) THEN
          tmpvaza(ix) = -tmpvaza(ix-1)
          EXIT
        ENDIF
      ENDDO

      edsza(0:nxtrack, i) = tmpsza(0:nxtrack)
      edsazm(0:nxtrack, i) = tmpsaza(0:nxtrack)
      edvza(0:nxtrack, i) = tmpvza(0:nxtrack)
      edvazm(0:nxtrack, i) = tmpvaza(0:nxtrack)
    ENDDO

    !PRINT *
    !WRITE(www_lun, *) 'Solar Angle'
    !WRITE(www_lun, '(10F9.2)') esza(1:nxtrack, 0)
    !WRITE(www_lun, *) 'View Angle'
    !WRITE(www_lun, '(10F9.2)') evza(1:nxtrack, 0)
    !WRITE(www_lun, *) 'Relative Azmimuthal Angle'
    !WRITE(www_lun, '(10F9.2)') eaza(1:nxtrack, 0)
    !WRITE(www_lun, *) 'Scattering Angle'
    !WRITE(www_lun, '(10F9.2)') esca(1:nxtrack, 0)

    ! Perform interpolation along the track with a simipler but simpler (center) approach
    ! since the pixel size along the track does not vary much
    IF (ntimes == 1) tmpdisy(0) = 0.00212031  ! ~ 13.5 km
    DO ix = 0, nxtrack
      DO i = 0, ntimes - 2 
        tmpdisy(i) = circle_rdis(elat(ix, i), elon(ix, i), elat(ix, i+1), elon(ix, i+1))
        CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1), &
             elon(ix, i+1), tmpdisy(i), tmpdisy(i)*0.5, clat(ix, i+1), clon(ix, i+1))  
        !print *, i+1, clat(ix, i+1), clon(ix, i+1)
      ENDDO

      i = 0
      CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1),    &
           elon(ix, i+1), tmpdisy(i), -tmpdisy(i)*0.5, clat(ix, i), clon(ix, i))
      !print *, i, clat(ix, i), clon(ix, i)

      i = ntimes - 1       
      CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i-1),    &
           elon(ix, i-1), tmpdisy(i-1), -tmpdisy(i-1)*0.5, clat(ix, i+1), clon(ix, i+1))
      !print *, i+1, clat(ix, i+1), clon(ix, i+1)
    ENDDO

    ! Perform coadding
    IF (nxbin > 1 .OR. nybin > 1) THEN
      nx = nxtrack / nxbin    
      ny = ntimes  / nybin

      ! cornor coordinates (only need sampling)
      j = 0
      DO ix = 0, nxtrack, nxbin
        clon(j, :) = clon(ix, :)
        clat(j, :) = clat(ix, :)
        j = j + 1
      ENDDO

      j = 0
      DO i = 0, ntimes, nybin
        clon(0:nx, j) = clon(0:nx, i)
        clat(0:nx, j) = clat(0:nx, i)
        j = j + 1
      ENDDO

      ! edge coordinates (easy to be re-computed from cornor coordinates)
      DO ix = 0, nx
        DO i = 0, ny - 1 
          tmpdisy(i) = circle_rdis(clat(ix, i), clon(ix, i), clat(ix, i+1), clon(ix, i+1))
          CALL sphergeom_intermediate(clat(ix, i), clon(ix, i), clat(ix, i+1), clon(ix, i+1), &
               tmpdisy(i), tmpdisy(i)*0.5, elat(ix, i), elon(ix, i))  
          !print *, tmpdisy(i), elat(ix, i)*rad2deg, elon(ix, i) * rad2deg
          !print *, clat(ix, i)*rad2deg, clon(ix, i)*rad2deg, clat(ix, i+1)*rad2deg, clon(ix, i+1)*rad2deg
        ENDDO
      ENDDO

      ! Center coordinates (computed from edge coordinates)
      DO ix = 1, nx
        DO i = 0, ny - 1 
          tmpdisx(ix) = circle_rdis(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i))
          CALL sphergeom_intermediate(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i), &
               tmpdisx(ix), tmpdisx(ix)*0.5, lat(ix, i), lon(ix, i))
          !print *, elat(ix-1, i)*rad2deg, elon(ix-1, i)*rad2deg, elat(ix, i)*rad2deg, elon(ix, i)*rad2deg
          !print *, lat(ix, i) *rad2deg, lon(ix, i)*rad2deg
        ENDDO
      ENDDO

      ! Average edge viewing geometries along the track, sample along the track
      i = 0
      DO ix = 0, nxtrack, nxbin
        jj = 0
        DO j = 0, ntimes-1, nybin
          edsza (i, jj) = SUM(edsza (ix, j:j+nybin-1)) / nybin
          edsazm(i, jj) = SUM(edsazm(ix, j:j+nybin-1)) / nybin
          edvza (i, jj) = SUM(edvza (ix, j:j+nybin-1)) / nybin
          edvazm(i, jj) = SUM(edvazm(ix, j:j+nybin-1)) / nybin
          jj = jj + 1
        ENDDO
        i = i + 1
      ENDDO

      ! Compute center viewing geometries (interpolate across the track)
      DO i = 0, ny-1
        tmpx = 0.0
        DO ix = 1, nx
          tmpx(ix) = tmpx(ix-1) + circle_rdis(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i))
        ENDDO
        tmpxmid(1:nx) = (tmpx(0:nx-1) + tmpx(1:nx)) / 2.0

        CALL interpol(tmpx(0:nx), edsza (0:nx, i),  nx+1, tmpxmid(1:nx), sza (1:nx, i),  nx, errstat) 
        CALL interpol(tmpx(0:nx), edsazm(0:nx, i),  nx+1, tmpxmid(1:nx), saza(1:nx, i),  nx, errstat)  

        CALL interpol(tmpx(0:nx), edvza (0:nx, i),  nx+1, tmpxmid(1:nx), vza (1:nx, i),  nx, errstat) 
        CALL interpol(tmpx(0:nx), edvazm(0:nx, i),  nx+1, tmpxmid(1:nx), vaza(1:nx, i),  nx, errstat)      

        ! Check center pixel
        mpix = nx / 2
        DO ix = mpix - 1, mpix + 1
          IF (vaza(ix, i) < 0) THEN
            vaza(ix, i) = -vaza(ix-1, i)
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ELSE
      nx = nxtrack
      ny = ntimes
    ENDIF

    ! Now compute effective viewing geometry
    ! Compute effective viewing geometry for each pixel (at a certain atmosphere) 
    DO i = 0, ny-1
      DO ix = 1, nx
        IF ( sza(ix, i)  >= minza  .AND. sza(ix, i)  < maxza  .AND. &
             vza(ix, i)  >= minza  .AND. vza(ix, i)  < maxza  .AND. &
             saza(ix, i) >= minaza .AND. saza(ix, i) < maxaza .AND. &
             vaza(ix, i) >= minaza .AND. vaza(ix, i) < maxaza) THEN
          zen0(1) = edsza(ix-1, i) 
          zen0(2) = sza(ix, i) 
          zen0(3) = edsza(ix, i)
          zen(1)  = edvza(ix-1, i) 
          zen(2)  = vza(ix, i) 
          zen(3)  = edvza(ix, i)
          sazm(1) = edsazm(ix-1, i)
          sazm(2) = saza(ix, i)
          sazm(3) = edsazm(ix, i)
          vazm(1) = edvazm(ix-1, i)
          vazm(2) = vaza(ix, i)
          vazm(3) = edvazm(ix, i)
          relaza  = ABS(vazm - sazm)

          WHERE (vazm > 0) 
            zen = - zen
          ENDWHERE

          ! -180 < relaza < 180          
          WHERE (relaza > 180.0)
            relaza = 360.0 - relaza
          ENDWHERE

          CALL omi_angle_sat2toa (3, zen0, zen, relaza, esza(ix, i), evza(ix, i), eaza(ix, i), esca(ix, i) )
          !print *, zen0 * rad2deg
          !print *, zen  * rad2deg
          !print *, 180.0 - relaza * rad2deg
          !print *, esza(ix, i), evza(ix, i), eaza(ix, i), esca(ix, i)
          !STOP             
        ENDIF
      ENDDO
    ENDDO

    clon(0:nx, 0:ny)   = clon(0:nx, 0:ny)   * rad2deg    
    clat(0:nx, 0:ny)   = clat(0:nx, 0:ny)   * rad2deg
    lon(1:nx, 0:ny-1)  = lon(1:nx, 0:ny-1)  * rad2deg    
    lat(1:nx, 0:ny-1)  = lat(1:nx, 0:ny-1)  * rad2deg
    elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) * rad2deg
    elat(0:nx, 0:ny-1) = elat(0:nx, 0:ny-1) * rad2deg

    WHERE (clon(0:nx, 0:ny) > 180.0)
      clon(0:nx, 0:ny) = clon(0:nx, 0:ny) - 360.0
    ENDWHERE
    WHERE (clon(0:nx, 0:ny) < -180.0)
      clon(0:nx, 0:ny) = clon(0:nx, 0:ny) + 360.0
    ENDWHERE

    WHERE (lon(1:nx, 0:ny-1) > 180.0)
      lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) - 360.0
    ENDWHERE
    WHERE (lon(1:nx, 0:ny-1) < -180.0)
      lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) + 360.0
    ENDWHERE

    WHERE (elon(0:nx, 0:ny-1) > 180.0)
      elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) - 360.0
    ENDWHERE
    WHERE (elon(0:nx, 0:ny-1) < -180.0)
      elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) + 360.0
    ENDWHERE


    RETURN

  END SUBROUTINE get_sphgeoview_corners



  !  Unused
  !
  !  ! Compute spherical distance between two points
  !  ! Haversine Formula (from R.W. Sinnott, "virtue of the Haversine", 
  !  ! Sky and Telescope V68 (2), 1984, p159  
  !  ! lat1, lon1, lat2, lon2 in radians
  !  FUNCTION circle_dis(lat1, lon1, lat2, lon2) RESULT(dis)
  !
  !    IMPLICIT NONE
  !
  !    ! ----------------------
  !    ! Input/output variables
  !    ! -----------------------
  !    REAL (KIND=r8), INTENT (IN) :: lat1, lon1, lat2, lon2
  !    REAL (KIND=r8)              :: dis
  !
  !    ! Local variable
  !    REAL (KIND=r8)              :: dlon, dlat, a, rdis, mlat
  !
  !    dlat = lat2 - lat1;    dlon = lon2 - lon1; mlat = (lat1 + lat2) / 2.0
  !
  !    a = MIN(1.0, SQRT( SIN(dlat/2.0)**2.0 + COS(lat1) * COS(lat2) * SIN(dlon/2.0)**2.0 )  )
  !    rdis = 2.0 * ASIN(a)                                    ! relative distaince in radiances
  !    dis = rdis * (rearth0 - 21.0 * SIN(mlat))               ! in the same unit of rearth0
  !
  !    RETURN
  !
  !  END FUNCTION circle_dis


  FUNCTION circle_rdis(lat1, lon1, lat2, lon2) RESULT(rdis)

    IMPLICIT NONE

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    REAL (KIND=r8), INTENT (IN) :: lat1, lon1, lat2, lon2
    REAL (KIND=r8)              :: rdis

    ! Local variable
    REAL (KIND=r8)              :: dlon, dlat, a

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = MIN(1.0, SQRT( SIN(dlat/2.0)**2.0 + COS(lat1) * COS(lat2) * SIN(dlon/2.0)**2.0 )  )
    rdis = 2.0 * ASIN(a)                                    ! relative distaince in radiances

    RETURN

  END FUNCTION circle_rdis




  ! This one works, gam is the longitude difference (-pi < gam0 < pi) between 2 and 1 (phi2 - phi1)
  SUBROUTINE sphergeom_intermediate ( lat1, lon1, lat2, lon2, c0, c, lat, lon )

    ! -----------------------------------------------------------------
    ! Finds the co-ordinates of C the baseline extended from two
    ! lon/lat points (A, B) on a sphere given the hypotenuse C_IN.
    ! ----------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: lat1, lat2, lon1, lon2, c0, c

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8),  INTENT (OUT)  :: lat, lon

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8)  :: x, y, z, tmp1, tmp2, frc, gamsign, theta, gam0, gam

    lat = 0.0_r8 
    lon = 0.0_r8    
    gam0 = angle_minus_twopi ( lon2 - lon1, pi )
    gamsign = ABS(gam0) / gam0
    gam0 = ABS(gam0)

    ! Get straight line (AB) segment fraction frc intercepted by the line from center to C
    ! If frc < 0, extrapolation, but it is limited to |c| < (180-c0)/2.0
    tmp1 = SIN(c)
    frc = tmp1 / (SIN(c0 - c) + tmp1)

    ! Work in Cartesian Coordinate 
    tmp1 = frc * COS(lat2)
    tmp2 = 1.0 - frc
    x = tmp2 *   COS(lat1) + tmp1 * COS(gam0)
    y = tmp1 *   SIN(gam0)
    z = tmp2 *   SIN(lat1) + frc * SIN(lat2)

    gam = ATAN(y/x)                          ! -90 < gam < 90
    IF (frc >= 0) THEN
      IF (gam < 0) gam = gam + pi           ! 0 <= gam <= 180
    ELSE
      IF (gam > 0) gam = gam - pi           ! -180 <= gam <= 0
    ENDIF
    gam = gamsign * gam                      ! Get correct sign
    lon = gam + lon1

    theta = ATAN (SQRT(x**2 + y**2) / z)     ! -90 < theta < 90
    IF (theta < 0) theta = theta + pi        ! 0 <= theta <= 180
    lat = pihalf - theta                     ! Convert to latitude

    RETURN
  END SUBROUTINE sphergeom_intermediate


END MODULE OMSAO_pixelcorner_module
