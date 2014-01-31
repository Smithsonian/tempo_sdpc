MODULE omi_read_l1b_data
  INCLUDE 'hdf.f90'
CONTAINS

  SUBROUTINE omi_read_binning_factor ( &
      l1bfile, l1bswath, ntimes, binfac, yn_szoom, errstat )

    USE OMSAO_precision_module
    USE OMSAO_omidata_module,    ONLY : global_mode, szoom_mode
    use l1bread

    IMPLICIT NONE

    ! --------------
    ! Input Variable
    ! --------------
    CHARACTER (LEN=*), INTENT (IN) :: l1bfile, l1bswath
    INTEGER (KIND=i4), INTENT (IN) :: ntimes

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat
    INTEGER (KIND=i1), DIMENSION (0:ntimes-1), INTENT (OUT)   :: binfac
    LOGICAL,           DIMENSION (0:ntimes-1), INTENT (OUT)   :: yn_szoom
    !
    type (L1B_Object_Type) :: l1bobj

    if (errstat < 0) return

    write (*,*) "DEBUG: In omi_read_binning_factor, l1bfile=", &
      TRIM(l1bfile), ", l1bswath=", TRIM(l1bswath), ", ntimes=", ntimes

    ! Allow error to flow
    call l1bread_open_swath (l1bfile, l1bswath, l1bobj, errstat)
    call l1bread_get1d_i1 (l1bobj, "ImageBinningFactor", 0, ntimes, binfac, errstat)
    call l1bread_close (l1bobj)
    if (errstat < 0) return

    ! ----------------------------------------------------------------------
    ! Check whether we have a Spatial Zoom granule, in which case we need to
    ! process 60 cross-track positions rather than 30.
    ! ----------------------------------------------------------------------
    IF ( ( INDEX(l1bfile, 'OML1BRUZ') > 0 ) .OR. &
        ( INDEX(l1bfile, 'OML1BRVZ') > 0 ) )    &
      binfac(0:ntimes-1) = global_mode

    ! ------------------------------------------------------------------------------
    ! Check for GLOBAL and SPATIAL ZOOM mode and set up arrays for index adjustment.
    ! ------------------------------------------------------------------------------
    WHERE ( binfac(0:ntimes-1) == szoom_mode )
      yn_szoom (0:nTimes-1) = .TRUE.
    ELSEWHERE
      yn_szoom (0:nTimes-1) = .FALSE.
    END WHERE

    RETURN
  END SUBROUTINE omi_read_binning_factor

  FUNCTION L1Bga_EarthSunDistance( he4filename, swathname ) RESULT( ESdistance )

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    CHARACTER (LEN=*), INTENT(IN) :: he4filename, swathname

    ! ---------------
    ! Result Variable
    ! ---------------
    REAL (KIND=r4) :: ESdistance

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=i4) :: swfid, swid, status

    ! ------------------
    ! External Functions
    ! ------------------
    INTEGER (KIND=i4) :: swopen, swattach, swrdattr, swdetach, swclose

    ESdistance = -1.0_r4

    swfid = swopen( he4filename, DFACC_READ )

    IF( swfid /=  -1) THEN
      swid = swattach( swfid, swathname )
      IF( swid /= -1 ) THEN
        status = swrdattr( swid, "EarthSunDistance", ESdistance )
        !IF( status == -1 ) THEN
        !   WRITE(msg,'(A)') "Get swath attribute EarthSunDistance"// &
        !        "failed from "//TRIM(swathname)//","  // &
        !        TRIM( he4filename )
        !   ierr = OMI_SMF_setmsg( OZT_E_INPUT,  msg, &
        !        "L1Bga_EarthSunDistance", zero )
        !ENDIF
        status = swdetach(swid)
      ENDIF
      status = swclose(swfid)
    ENDIF

  END FUNCTION L1Bga_EarthSunDistance

  SUBROUTINE omi_read_radiance_lines ( &
      l1bfile, iline, nxtrack, nloop, nwavel_ccd, errstat )

    USE OMSAO_precision_module
    !USE OMSAO_indices_module, ONLY: pge_o3_idx
    USE OMSAO_parameters_module, ONLY: &
      i2_missval, r4_missval, &
      min_zenith, min_azimuth, max_azimuth, &  ! "non-inclusive"
      max_latitude, max_longitude, &
      earth_radius_avg
    USE OMSAO_variables_module,  ONLY: zatmos
    USE OMSAO_omidata_module,  ONLY: &
      omi_radiance_swathname, omi_radiance_spec,  &
      omi_radiance_wavl, omi_radiance_qflg, omi_height, omi_geoflg, omi_latitude,             &
      omi_longitude, omi_szenith, omi_sazimuth, omi_vzenith, omi_vazimuth,                    &
      omi_razimuth, omi_auraalt, omi_time, omi_nwav_rad, omi_radiance_errstat,                &
      rad_ccdpix_selection,                                                                   &
      omi_xtrflg_l1b, omi_xtrflg
    USE OMSAO_errstat_module
    USE angle_sat2toa, ONLY: gnome_angle_sat2toa
    use l1bread

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*), INTENT (IN) :: l1bfile
    INTEGER (KIND=i4), INTENT (IN) :: iline, nloop, nxtrack, nwavel_ccd

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                           :: &
      iloop, imin, imax, ix, icnt
    REAL      (KIND=r4), DIMENSION (nxtrack)      :: tmp_sazm, tmp_vazm
    REAL      (KIND=r4), DIMENSION (nwavel_ccd,nxtrack,0:nloop-1) :: tmp_wvl, tmp_spc
    INTEGER   (KIND=i2), DIMENSION (nwavel_ccd,nxtrack,0:nloop-1) :: tmp_flg

    type (L1B_Object_Type) :: l1bobj

    omi_radiance_errstat = pge_errstat_ok

    ! let errstat flow
    call l1bread_open_swath (l1bfile, omi_radiance_swathname, l1bobj, errstat)

    call l1bread_get1d_r8 (l1bobj, "Time", iline, nloop, omi_time, errstat)
    call l1bread_get1d_r4 (l1bobj, "SpacecraftAltitude", iline, nloop, omi_auraalt, errstat)
    call l1bread_get2d_r4 (l1bobj, "Latitude", iline, nloop, omi_latitude, errstat)
    call l1bread_get2d_r4 (l1bobj, "Longitude", iline, nloop, omi_longitude, errstat)
    call l1bread_get2d_r4 (l1bobj, "SolarZenithAngle", iline, nloop, omi_szenith, errstat)
    call l1bread_get2d_r4 (l1bobj, "SolarAzimuthAngle", iline, nloop, omi_sazimuth, errstat)
    call l1bread_get2d_r4 (l1bobj, "ViewingZenithAngle", iline, nloop, omi_vzenith, errstat)
    call l1bread_get2d_r4 (l1bobj, "ViewingAzimuthAngle", iline, nloop, omi_vazimuth, errstat)
    call l1bread_get2d_i2 (l1bobj, "TerrainHeight", iline, nloop, omi_height, errstat)
    call l1bread_get2d_i2 (l1bobj, "GroundPixelQualityFlags", iline, nloop, omi_geoflg, errstat)
    call l1bread_get2d_i1 (l1bobj, "XTrackQualityFlags", iline, nloop, omi_xtrflg_l1b, errstat)

    call l1bread_get3d_r4 (l1bobj, "Radiance", iline, nloop, tmp_spc, errstat)
    call l1bread_get3d_i2 (l1bobj, "PixelQualityFlags", iline, nloop, tmp_flg, errstat)
    call l1bread_get3d_r4 (l1bobj, "Wavelength", iline, nloop, tmp_wvl, errstat)

    call l1bread_close (l1bobj)
    if (errstat < 0) return

    do iloop = 0, nloop-1
      ! --------------------------------------------------------
      ! Expand XTR Quality Flags to something easily parse-able.
      ! Assign Missing Values where necessary. gga
      ! --------------------------------------------------------
      CALL convert_xtqualflag_info ( &
        nxtrack, omi_xtrflg_l1b(1:nxtrack,iloop), omi_xtrflg(1:nxtrack,iloop))

      ! --------------------------------------------------------------
      ! Check for missing data and reinitialize to PGEs local MissVals
      ! --------------------------------------------------------------
      WHERE (omi_height(1:nxtrack,iloop) <= i2_missval )
        omi_height(1:nxtrack,iloop) = i2_missval
      ENDWHERE
      WHERE ( ABS(omi_latitude(1:nxtrack,iloop)) > max_latitude )
        omi_latitude(1:nxtrack,iloop) = r4_missval
      ENDWHERE
      WHERE ( ABS(omi_longitude(1:nxtrack,iloop)) > max_longitude )
        omi_longitude(1:nxtrack,iloop) = r4_missval
      ENDWHERE

      ! ---------------------------------------------------------------------
      ! For the Zenith Angles we only correct those values < 0, since we want
      ! to maintain the information of the value of SZA even if it is out of
      ! the bounds required for the computation of AMFs
      ! ---------------------------------------------------------------------
      WHERE (omi_szenith(1:nxtrack,iloop) < min_zenith )
        omi_szenith(1:nxtrack,iloop) = r4_missval
      ENDWHERE
      WHERE (omi_vzenith(1:nxtrack,iloop) < min_zenith )
        omi_vzenith(1:nxtrack,iloop) = r4_missval
      ENDWHERE
      WHERE ((omi_sazimuth(1:nxtrack,iloop) < min_azimuth) &
             .or. (omi_sazimuth(1:nxtrack,iloop) > max_azimuth))
        omi_sazimuth(1:nxtrack,iloop) = r4_missval
      ENDWHERE
      WHERE ((omi_vazimuth(1:nxtrack,iloop) < min_azimuth) &
             .or. (omi_vazimuth(1:nxtrack,iloop) > max_azimuth))
        omi_vazimuth(1:nxtrack,iloop) = r4_missval
      ENDWHERE

      ! ----------------------------------------------------------------
      ! Relative azimuth angles (requires some checks and adjustments).
      ! It also has become a Geolocation Field rather than a Data Field.
      ! ----------------------------------------------------------------
      ! (1) Map [-180, +180] to [0, 360]
      tmp_sazm(1:nxtrack) = omi_sazimuth(1:nxtrack,iloop)
      tmp_vazm(1:nxtrack) = omi_vazimuth(1:nxtrack,iloop)
      WHERE ( tmp_sazm(1:nxtrack) /= r4_missval .AND. &
             tmp_sazm(1:nxtrack) < 0.0_r4 .AND. tmp_sazm(1:nxtrack) >= -180.0_r4 )
        tmp_sazm(1:nxtrack) = 360.0_r4 + tmp_sazm(1:nxtrack)
      ENDWHERE
      WHERE ( tmp_vazm(1:nxtrack) /= r4_missval .AND. &
             tmp_vazm(1:nxtrack) < 0.0_r4 .AND. tmp_vazm(1:nxtrack) >= -180.0_r4 )
        tmp_vazm(1:nxtrack) = 360.0_r4 + tmp_vazm(1:nxtrack)
      ENDWHERE
      ! (2) Compute relative azimuth angle (RELATIVE means absolute value)
      WHERE ( tmp_sazm(1:nxtrack) >= 0.0_r4 .AND. tmp_vazm(1:nxtrack) >= 0.0_r4 )
        omi_razimuth(1:nxtrack,iloop) = ABS( tmp_sazm(1:nxtrack) - tmp_vazm(1:nxtrack) )
      ENDWHERE
      WHERE ( omi_razimuth(1:nxtrack,iloop) > 180.0_r4 )
        omi_razimuth(1:nxtrack,iloop) = 360.0_r4 - omi_razimuth(1:nxtrack,iloop)
      ENDWHERE

      ! -------------------------------------------------------------------------
      ! Compute zenith and relative azimuth angles at PGE TOA (from control file)
      ! -------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------
      ! NOTE: In the GOME data products all angles are given at the Spacecraft, while
      ! in the OMI data product, angles seem to be given at the surface. This difference
      ! requires a change of the adjustment.
      ! --------------------------------------------------------------------------------
      ! OMI
      IF ( zatmos > 0.0_r8 ) &
        CALL gnome_angle_sat2toa ( &
        earth_radius_avg, 0.0_r4, zatmos, nxtrack, &
        omi_szenith(1:nxtrack,iloop), omi_vzenith(1:nxtrack,iloop), omi_razimuth(1:nxtrack,iloop))
      !! GOME
      !CALL angle_sat2toa ( &
      !     earth_radius_avg, eos_aura_avgalt, zatmos, nxtrack, &
      !     omi_szenith(1:nxtrack,iloop), omi_vzenith(1:nxtrack,iloop), omi_razimuth(1:nxtrack,iloop))

      ! ----------------------------------------------
      ! Get radiances associated with wavelength range
      ! ----------------------------------------------
      DO ix = 1, nxtrack
        imin = rad_ccdpix_selection(ix,1)
        imax = rad_ccdpix_selection(ix,4)
        icnt = imax - imin + 1
        omi_radiance_wavl(1:icnt,ix,iloop) = REAL ( tmp_wvl(imin:imax,ix, iloop), KIND=r8 )
        omi_radiance_spec(1:icnt,ix,iloop) = REAL ( tmp_spc(imin:imax,ix, iloop), KIND=r8 )
        omi_radiance_qflg(1:icnt,ix,iloop) =        tmp_flg(imin:imax,ix, iloop)
        omi_nwav_rad     (       ix,iloop) = icnt
      END DO

    END DO
    RETURN

  END SUBROUTINE omi_read_radiance_lines

  SUBROUTINE omi_read_glint_ice_flags ( l1bfile, nx, nt, snow_ice_flg, glint_flg, errstat )

    USE OMSAO_precision_module
    USE OMSAO_omidata_module,    ONLY: omi_radiance_swathname
    use l1bread

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*), INTENT (IN) :: l1bfile
    INTEGER (KIND=i4), INTENT (IN) :: nx, nt

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nx,0:nt-1), INTENT (OUT)   :: glint_flg, snow_ice_flg
    INTEGER (KIND=i4),                        INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)               :: iline
    INTEGER (KIND=i2), DIMENSION (nx,0:nt-1) :: geoflg
    INTEGER (KIND=i2), DIMENSION (nx) :: land_water_flg
    type (L1B_Object_Type) :: l1bobj

    ! let errstat flow
    call l1bread_open_swath (l1bfile, omi_radiance_swathname, l1bobj, errstat)
    call l1bread_get2d_i2 (l1bobj, "GroundPixelQualityFlags", 0, nt, geoflg, errstat)
    call l1bread_close (l1bobj)
    if (errstat < 0) return

    ! ---------------------------------------------
    ! convert to snow/glint flags
    ! ------------------------------------------------------------------------
    ! NOTE: GLINT_FLG and SNOW_ICE_FLG are defined on the whole swath because
    !       they are being used in the AMF computation routine. Hence we need
    !       to save the full array.
    ! ------------------------------------------------------------------------
    DO iline = 0, nt-1

      CALL convert_gpqualflag_info (   &
        nx,                         &
        geoflg        (1:nx, iline),       &
        land_water_flg(1:nx),       &
        glint_flg     (1:nx,iline), &
        snow_ice_flg  (1:nx,iline)    )

    END DO

  END SUBROUTINE omi_read_glint_ice_flags

  SUBROUTINE convert_gpqualflag_info ( &
      nxtrack, omi_geoflg, land_water_flg, glint_flg, snow_ice_flg )

    USE OMSAO_precision_module
    USE strutils
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (IN) :: omi_geoflg

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i2),                PARAMETER      :: nbyte = 16
    INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = INT((/ 1, 2, 4, 8, 16, 32, 64 /),KIND=i2)
    INTEGER (KIND=i4)                                :: i
    INTEGER (KIND=i2), DIMENSION (nxtrack)           :: tmp_flg
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    land_water_flg = 0 ; glint_flg = 0 ; snow_ice_flg = 0

    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    tmp_flg(1:nxtrack) = omi_geoflg(1:nxtrack)  ;  tmp_bytes = 0

    ! -------------------------------------------------------------------
    ! CAREFUL: Only 15 flags/positions (0:14) can be returned or else the
    !          conversion will result in a numeric overflow.
    ! -------------------------------------------------------------------
    CALL convert_2bytes_to_16bits ( &
      nbyte-1, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-2) )

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

  SUBROUTINE convert_xtqualflag_info ( nxtrack, omi_xtrflg_l1b, omi_xtrflg )

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: i1_missval, i2_missval
    USE strutils

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (INOUT) :: omi_xtrflg_l1b

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: omi_xtrflg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4),                       PARAMETER    :: nbit = 16
    INTEGER (KIND=i1), DIMENSION (0:2),      PARAMETER    :: three_bit = INT((/ 1, 2, 4 /),KIND=i1)
    INTEGER (KIND=i2), DIMENSION (3:nbit-1), PARAMETER    :: add_value = INT((/ 10, 30, 100, 1000, 10000, &
                                                                              0,  0,   0,    0,     0, &
                                                                              0,  0,   0 /), KIND=i2)
    INTEGER (KIND=i4)                                     :: i
    INTEGER (KIND=i2), DIMENSION (nxtrack)                :: tmp_flg
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbit-1)       :: tmp_bits

    ! -----------------------------------------------------------------------
    ! Initialize output quantities. We can't initialize to "I2_MISSVAL" since
    ! we will be recursively adding values to OMI_XTRFLG and hence have to
    ! start out from Zero.
    ! -----------------------------------------------------------------------
    omi_xtrflg = 0_i2

    ! --------------------------------------------------------
    ! Save input variable in TMP_FLG for modification; in that
    ! process, perform a INT8 --> UINT8 conversion.
    ! --------------------------------------------------------
    tmp_bits = 0
    DO i = 1, nXtrack

      IF ( omi_xtrflg_l1b(i) > -127_i1 .AND. omi_xtrflg_l1b(i) < 0_i1 ) THEN
        tmp_flg(i) = INT ( omi_xtrflg_l1b(i), KIND=i2 ) + 256_i2
      ELSE
        tmp_flg(i) = INT ( omi_xtrflg_l1b(i), KIND=i2 )
      END IF

      ! -----------------------------------------------------------
      ! The code below fails on 32 bit platforms due to "255" being
      ! outside the range of 8bit Integers. On 64 bit platforms it
      ! works perfectly fine.
      ! -----------------------------------------------------------
      !IF ( omi_xtrflg_l1b(i) > -127_i1 .AND. omi_xtrflg_l1b(i) < 0_i1 ) THEN
      !   tmp_flg(i) = INT ( IAND(omi_xtrflg_l1b(i),255), KIND=i2 )
      !ELSE
      !   tmp_flg(i) = INT ( omi_xtrflg_l1b(i), KIND=i2 )
      !END IF

    END DO

    ! gga
    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    !tmp_flg(1:nxtrack) = omi_xtrflg_l1b(1:nxtrack)  ;  tmp_bits = 0
    ! gga
    CALL convert_2bytes_to_16bits ( &
      nbit, nxtrack, tmp_flg(1:nxtrack), tmp_bits(1:nxtrack,0:nbit-1) )
    ! ------------------------------------------------------------------
    ! The Row Anomaly Flags, Bits 0-2
    ! ------------------------------------------------------------------
    DO i = 1, nxtrack
      omi_xtrflg(i) = INT ( SUM(tmp_bits(i,0:2 )*three_bit(0:2)), KIND=i2 )
    END DO

    ! ----------------------------------------------
    ! Add the other bit flags:
    !
    !  Bit  Effect                      Added Value
    !   3   Reserved for future use        10
    !   4   wavelength-shift               30
    !   5   blockage                      100
    !   6   stray sunlight               1000
    !   7   stray earth radiance        10000
    ! ----------------------------------------------
    DO i = 1, nxtrack
      IF ( omi_xtrflg(i) < 0_i2 ) THEN
        omi_xtrflg(i)     = i2_missval
        omi_xtrflg_l1b(i) = i1_missval
      ELSE
        omi_xtrflg(i) = omi_xtrflg(i) + SUM(INT(tmp_bits(i,3:nbit-1),KIND=i2) * add_value(3:nbit-1))
      END IF
    END DO

    RETURN
  END SUBROUTINE convert_xtqualflag_info

  SUBROUTINE omi_xtract_swathname ( l1bfile, l1bchan, omiswath )

    USE OMSAO_precision_module

    USE OMSAO_parameters_module, ONLY : MAX_STR_LEN
    USE hdfeos4_parameters
    USE strutils

    IMPLICIT NONE

    ! --------------
    ! Input Variable
    ! --------------
    CHARACTER (LEN=*), INTENT (IN) :: l1bfile
    CHARACTER (LEN=3), INTENT (IN) :: l1bchan

    ! ----------------
    ! Output variables
    ! ----------------
    CHARACTER (LEN=*), INTENT (OUT) :: omiswath

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)      :: is, ie
    INTEGER   (KIND=i4)      :: swfid, nswath, strbufsize, xswath
    CHARACTER (LEN=MAX_STR_LEN) :: swathlist

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=20), PARAMETER :: modulename = 'omi_xtract_swathname'

    ! --------------------------
    ! Initialize OUTPUT variable
    ! --------------------------
    omiswath = '?'

    ! ---------------------------------------------------------
    ! Inquire about the swaths in the current L1b radiance file
    ! ---------------------------------------------------------
    swfid  = SWOpen     ( l1bfile, DFACC_READ )
    swathlist="" !JED
    nswath = SWInqswath ( l1bfile, swathlist, strbufsize )
    xswath = SWClose    ( swfid )

    ! --------------------------------------------------------------------
    ! Extract the swath name we need. Either there is one one swath in the
    ! file (VIS) or there are two (UV-1, UV-2)
    ! --------------------------------------------------------------------
    is = 1
    SELECT CASE ( nswath )
    CASE ( 1 )
      is = 1 ; ie = strbufsize
    CASE ( 2 )
      CALL find_endstring ( strbufsize, swathlist, 1, ie )
      SELECT CASE ( l1bchan )
      CASE ( 'UV1' )
        is = 1 ; ie = ie-1
      CASE ( 'UV2' )
        is = ie+1 ; ie = strbufsize
      CASE DEFAULT
        ! Nothing to do here except to fold
      END SELECT
    CASE DEFAULT
      ! Nothing to do here except to fold.
    END SELECT

    omiswath = TRIM(ADJUSTL(swathlist(is:ie)))

    RETURN
  END SUBROUTINE omi_xtract_swathname
END MODULE
