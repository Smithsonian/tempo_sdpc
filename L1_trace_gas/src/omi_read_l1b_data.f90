MODULE omi_read_l1b_data
  use tell_module
  INCLUDE 'hdf.f90'

CONTAINS

  SUBROUTINE omi_read_binning_factor ( &
      l1bfile, l1bswath, ntimes, binfac, is_szoom, errstat )

    USE OMSAO_precision_module
    USE OMSAO_omidata_module,    ONLY : global_mode, szoom_mode
    !use l1bread
    use tio_module

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
    LOGICAL,           DIMENSION (0:ntimes-1), INTENT (OUT)   :: is_szoom
    !
    !type (L1B_Object_Type) :: l1bobj
    type (tiof_l1_object_type) :: tio_l1obj

    if (errstat < 0) return

    write (*,*) "DEBUG: In omi_read_binning_factor, l1bfile=", &
      TRIM(l1bfile), ", l1bswath=", TRIM(l1bswath), ", ntimes=", ntimes

    ! Allow error to flow
    !call l1bread_open_swath (l1bfile, l1bswath, l1bobj, errstat)
    !call l1bread_get1d_i1 (l1bobj, "ImageBinningFactor", 0, ntimes, binfac, errstat)
    !call l1bread_close (l1bobj)
    call tiof_open (l1bfile, tio_l1obj, errstat)
    call tiof_inq_group (tio_l1obj, l1bswath, errstat)
    call tiof_get1d_i1 (tio_l1obj, "ImageBinningFactor", 0, ntimes, binfac, errstat)
    call tiof_close (tio_l1obj, errstat)
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
      is_szoom (0:nTimes-1) = .TRUE.
    ELSEWHERE
      is_szoom (0:nTimes-1) = .FALSE.
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
    !use l1bread
    use tio_module

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

    !type (L1B_Object_Type) :: l1bobj
    type (tiof_l1_object_type) :: tio_l1obj

    omi_radiance_errstat = pge_errstat_ok

    write(*,*)'omi_read_radiance_lines: reading '//trim(omi_radiance_swathname)
    ! let errstat flow
    !call l1bread_open_swath (l1bfile, omi_radiance_swathname, l1bobj, errstat)
    !call l1bread_get1d_r8 (l1bobj, "Time", iline, nloop, omi_time, errstat)
    !call l1bread_get1d_r4 (l1bobj, "SpacecraftAltitude", iline, nloop, omi_auraalt, errstat)
    !call l1bread_get2d_r4 (l1bobj, "Latitude", iline, nloop, omi_latitude, errstat)
    !call l1bread_get2d_r4 (l1bobj, "Longitude", iline, nloop, omi_longitude, errstat)
    !call l1bread_get2d_r4 (l1bobj, "SolarZenithAngle", iline, nloop, omi_szenith, errstat)
    !call l1bread_get2d_r4 (l1bobj, "SolarAzimuthAngle", iline, nloop, omi_sazimuth, errstat)
    !call l1bread_get2d_r4 (l1bobj, "ViewingZenithAngle", iline, nloop, omi_vzenith, errstat)
    !call l1bread_get2d_r4 (l1bobj, "ViewingAzimuthAngle", iline, nloop, omi_vazimuth, errstat)
    !call l1bread_get2d_i2 (l1bobj, "TerrainHeight", iline, nloop, omi_height, errstat)
    !call l1bread_get2d_i2 (l1bobj, "GroundPixelQualityFlags", iline, nloop, omi_geoflg, errstat)
    !call l1bread_get2d_i1 (l1bobj, "XTrackQualityFlags", iline, nloop, omi_xtrflg_l1b, errstat)
    !call l1bread_get3d_r4 (l1bobj, "Radiance", iline, nloop, tmp_spc, errstat)
    !call l1bread_get3d_i2 (l1bobj, "PixelQualityFlags", iline, nloop, tmp_flg, errstat)
    !call l1bread_get3d_r4 (l1bobj, "Wavelength", iline, nloop, tmp_wvl, errstat)
    !call l1bread_close (l1bobj)
    call tiof_open (l1bfile, tio_l1obj, errstat)
    call tiof_get1d_r8 (tio_l1obj, "time", iline, nloop, omi_time, errstat)
    write(*,*)' tiof_inq_group: opening swath='//trim(omi_radiance_swathname)
    call tiof_inq_group (tio_l1obj, omi_radiance_swathname, errstat)
    call tiof_get1d_r4 (tio_l1obj, "SpacecraftAltitude", iline, nloop, omi_auraalt, errstat)
    call tiof_get2d_r4 (tio_l1obj, "latitude", iline, nloop, omi_latitude, errstat)
    call tiof_get2d_r4 (tio_l1obj, "longitude", iline, nloop, omi_longitude, errstat)
    call tiof_get2d_r4 (tio_l1obj, "solar_zenith_angle", iline, nloop, omi_szenith, errstat)
    call tiof_get2d_r4 (tio_l1obj, "solar_azimuth_angle", iline, nloop, omi_sazimuth, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_zenith_angle", iline, nloop, omi_vzenith, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_azimuth_angle", iline, nloop, omi_vazimuth, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ellipsoid_altitude", iline, nloop, omi_height, errstat)
    call tiof_get2d_i2 (tio_l1obj, "GroundPixelQualityFlags", iline, nloop, omi_geoflg, errstat)
    call tiof_get2d_i1 (tio_l1obj, "XTrackQualityFlags", iline, nloop, omi_xtrflg_l1b, errstat)
    call tiof_get3d_r4 (tio_l1obj, "radiance", iline, nloop, tmp_spc, errstat)
    call tiof_get3d_i2 (tio_l1obj, "data_quality_flag", iline, nloop, tmp_flg, errstat)
    call tiof_get3d_r4 (tio_l1obj, "wavelength", iline, nloop, tmp_wvl, errstat)
    call tiof_close (tio_l1obj, errstat)
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
        ! FIXME (JCH) this can fail if icnt > declared size, nwavel_max
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
    !use l1bread
    use tio_module

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
    INTEGER (KIND=i2), DIMENSION (nx,0:nt-1) :: geoflg
    !type (L1B_Object_Type) :: l1bobj
    type (tiof_l1_object_type) :: tio_l1obj

    ! let errstat flow
    !call l1bread_open_swath (l1bfile, omi_radiance_swathname, l1bobj, errstat)
    !call l1bread_get2d_i2 (l1bobj, "GroundPixelQualityFlags", 0, nt, geoflg, errstat)
    !call l1bread_close (l1bobj)
    call tiof_open (l1bfile, tio_l1obj, errstat)
    call tiof_inq_group (tio_l1obj, omi_radiance_swathname, errstat)
    call tiof_get2d_i2 (tio_l1obj, "GroundPixelQualityFlags", 0, nt, geoflg, errstat)
    call tiof_close (tio_l1obj, errstat)
    if (errstat < 0) return

    ! ---------------------------------------------
    ! convert to snow/glint flags
    ! ------------------------------------------------------------------------
    ! NOTE: GLINT_FLG and SNOW_ICE_FLG are defined on the whole swath because
    !       they are being used in the AMF computation routine. Hence we need
    !       to save the full array.
    ! ------------------------------------------------------------------------

    ! Bits 0-3 are land/water -- not used here
    ! land_water_flg = iand (geoflg, 15_i2)

    ! Bit 4 is glint
    glint_flg = iand (ishft(geoflg, -4), 1_i2)

    ! Bits 8-14 are snow/ice
    snow_ice_flg  = iand (ishft(geoflg, -8), 127_i2)

  END SUBROUTINE omi_read_glint_ice_flags

  SUBROUTINE convert_xtqualflag_info ( nxtrack, loc_xtrflg_l1b, loc_xtrflg )

    USE OMSAO_precision_module
    USE strutils

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (INOUT) :: loc_xtrflg_l1b

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: loc_xtrflg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i2), DIMENSION (3:7), PARAMETER :: &
      add_value = INT((/ 10, 30, 100, 1000, 10000 /), kind=i2)
    INTEGER (KIND=i4) :: i

    ! The Row Anomaly Flags, Bits 0-2, add 1,2,4
    loc_xtrflg(1:nxtrack) = iand (loc_xtrflg_l1b(1:nxtrack), 7_i1)
    ! others:
    !  Bit  Effect                      Added Value
    !   3   Reserved for future use        10
    !   4   wavelength-shift               30
    !   5   blockage                      100
    !   6   stray sunlight               1000
    !   7   stray earth radiance        10000
    ! ----------------------------------------------
    do i = 3, 7
      where (0 /= iand (loc_xtrflg_l1b(1:nxtrack), ishft (1_i1, i)))
        loc_xtrflg(1:nxtrack) = loc_xtrflg(1:nxtrack) + add_value(i)
      end where
    end do
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
