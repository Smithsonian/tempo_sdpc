MODULE omi_read_l1b_data
  use tell_module
!  INCLUDE 'hdf.f90'

  private
  public omi_read_binning_factor, read_earth_sun_distance, &
    omi_read_radiance_lines, omi_read_glint_ice_flags

CONTAINS

  SUBROUTINE omi_read_binning_factor ( &
      l1bfile, l1bswath, ntimes, binfac, is_szoom, errstat )

    USE OMSAO_precision_module
    USE OMSAO_omidata_module,    ONLY : global_mode, szoom_mode
    !use l1bread
    use tio_module
    use tg_names_module
    use netcdf, only : nf90_nowrite
    use ctrlvars, only : yn_omi_data

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
    type (tiof_file_type) :: tio_l1obj
    character (len=256) :: logmsg

    if (errstat /= 0) return

    write (logmsg, '(a, i5)')"DEBUG: In omi_read_binning_factor, l1bfile="// &
      TRIM(l1bfile)//", l1bswath="//TRIM(l1bswath)//", ntimes=", ntimes
    call tell_log (1, logmsg)

    ! Allow error to flow
    !call l1bread_open_swath (l1bfile, l1bswath, l1bobj, errstat)
    !call l1bread_get1d_i1 (l1bobj, "ImageBinningFactor", 0, ntimes, binfac, errstat)
    !call l1bread_close (l1bobj)
    if (yn_omi_data) then
      call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
      call tiof_inq_group (tio_l1obj, l1bswath, errstat)
      call tiof_get1d_i1 (tio_l1obj, "ImageBinningFactor", [0], [ntimes], &
           binfac, errstat)
      call tiof_close (tio_l1obj, errstat)
      if (errstat /= 0) return

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

    else ! TEMPO data - no zoom
      binfac = global_mode
      is_szoom = .false.
    endif


    RETURN
  END SUBROUTINE omi_read_binning_factor

  subroutine read_earth_sun_distance (filename, dist, errstat)
    use tio_module
    use tg_names_module
    use netcdf, only : nf90_nowrite
    use OMSAO_parameters_module, only : r4
    implicit none
    character (len=*), intent(in) :: filename
    real (kind=r4), intent(out) :: dist
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj

    if (errstat /= 0) return

    call tiof_open (filename, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "Error opening file: "//trim(filename), errstat)
      return
    endif

    call tiof_get_r4 (obj, tg_var_earth_sun_distance, dist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "Error reading earth-sun distance from file"//trim(filename), &
                       errstat)
    endif

    call tiof_close (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, "Error closing file: "//trim(filename), errstat)
      return
    endif

  end subroutine read_earth_sun_distance

!unused  FUNCTION L1Bga_EarthSunDistance( he4filename, swathname ) RESULT( ESdistance )
!unused
!unused     USE OMSAO_precision_module
!unused     IMPLICIT NONE
!unused
!unused     ! ---------------
!unused     ! Input Variables
!unused     ! ---------------
!unused     CHARACTER (LEN=*), INTENT(IN) :: he4filename, swathname
!unused
!unused     ! ---------------
!unused     ! Result Variable
!unused     ! ---------------
!unused     REAL (KIND=r4) :: ESdistance
!unused
!unused     ! ---------------
!unused     ! Local Variables
!unused     ! ---------------
!unused     INTEGER (KIND=i4) :: swfid, swid, status
!unused
!unused     ! ------------------
!unused     ! External Functions
!unused     ! ------------------
!unused     INTEGER (KIND=i4) :: swopen, swattach, swrdattr, swdetach, swclose
!unused
!unused     ESdistance = -1.0_r4
!unused
!unused     swfid = swopen( he4filename, DFACC_READ )
!unused
!unused     IF( swfid /=  -1) THEN
!unused       swid = swattach( swfid, swathname )
!unused       IF( swid /= -1 ) THEN
!unused         status = swrdattr( swid, "EarthSunDistance", ESdistance )
!unused         !IF( status == -1 ) THEN
!unused         !   WRITE(msg,'(A)') "Get swath attribute EarthSunDistance"// &
!unused         !        "failed from "//TRIM(swathname)//","  // &
!unused         !        TRIM( he4filename )
!unused         !   ierr = OMI_SMF_setmsg( OZT_E_INPUT,  msg, &
!unused         !        "L1Bga_EarthSunDistance", zero )
!unused         !ENDIF
!unused         status = swdetach(swid)
!unused       ENDIF
!unused       status = swclose(swfid)
!unused     ENDIF
!unused
!unused   END FUNCTION L1Bga_EarthSunDistance

  SUBROUTINE omi_read_radiance_lines ( &
      l1bfile, iline, nxtrack, nloop, nwavel_ccd, errstat )

    USE OMSAO_precision_module
    !USE OMSAO_indices_module, ONLY: pge_o3_idx
    USE OMSAO_parameters_module, ONLY: &
      i2_missval, r4_missval, &
      min_zenith, min_azimuth, max_azimuth, &  ! "non-inclusive"
      max_legal_zenith, &
      max_latitude, max_longitude, &
      earth_radius_avg
    USE OMSAO_variables_module,  ONLY: zatmos
    USE OMSAO_omidata_module,  ONLY: &
      omi_radiance_swathname, omi_radiance_spec,  &
      omi_radiance_wavl, omi_radiance_qflg, omi_height, omi_geoflg, omi_latitude,             &
      omi_longitude, omi_szenith, omi_sazimuth, omi_vzenith, omi_vazimuth,                    &
      omi_razimuth, omi_auraalt, omi_time, omi_nwav_rad, & !omi_radiance_errstat,                &
      rad_ccdpix_selection,                                                                   &
      omi_xtrflg_l1b, omi_xtrflg
    !USE OMSAO_errstat_module
    USE angle_sat2toa, ONLY: gnome_angle_sat2toa
    !use l1bread
    use tio_module
    use tg_names_module
    use netcdf, only : nf90_nowrite
    use ctrlvars, only : yn_omi_data

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
    character (len=256) :: logmsg

    !type (L1B_Object_Type) :: l1bobj
    type (tiof_file_type) :: tio_l1obj

    !omi_radiance_errstat(:) = pge_errstat_ok

    write (logmsg, '(a,i4,a)')'omi_read_radiance_lines: iline=',iline, &
      ' reading swathname='//trim(omi_radiance_swathname)//' file='//trim(l1bfile)
    call tell_log (2, logmsg)

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
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, tg_var_time, [iline], [nloop], omi_time, errstat)
    call tiof_inq_group (tio_l1obj, omi_radiance_swathname, errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_latitude, [iline,0], [nloop,nxtrack], &
                        omi_latitude(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_longitude, [iline,0], [nloop,nxtrack], &
                        omi_longitude(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_sz_angle, [iline,0], [nloop,nxtrack], &
                        omi_szenith(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_sa_angle, [iline,0], [nloop,nxtrack], &
                        omi_sazimuth(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_vz_angle, [iline,0], [nloop,nxtrack], &
                        omi_vzenith(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_va_angle, [iline,0], [nloop,nxtrack], &
                        omi_vazimuth(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_i2 (tio_l1obj, tg_var_terrain_height, [iline,0], [nloop,nxtrack], &
                        omi_height(1:nxtrack,0:nloop-1), errstat)
    call tiof_get2d_ui4 (tio_l1obj, tg_var_gpqf, [iline,0], [nloop,nxtrack], &
                        omi_geoflg(1:nxtrack,0:nloop-1), errstat)
    call tiof_get3d_r4 (tio_l1obj, tg_var_radiance, [iline,0,0], [nloop,nxtrack,nwavel_ccd], &
                        tmp_spc(:,1:nxtrack,0:nloop-1), errstat)
    call tiof_get3d_i2 (tio_l1obj, tg_var_pqf, [iline,0,0], [nloop,nxtrack,nwavel_ccd], &
                        tmp_flg(:,1:nxtrack,0:nloop-1), errstat)
    call tiof_get3d_r4 (tio_l1obj, tg_var_wavelength, [iline,0,0], [nloop,nxtrack,nwavel_ccd], &
                        tmp_wvl(:,1:nxtrack,0:nloop-1), errstat)
    if (yn_omi_data) then
      call tiof_get1d_r4 (tio_l1obj, "SpacecraftAltitude", [iline], &
           [nloop], omi_auraalt, errstat)
      call tiof_get2d_i1 (tio_l1obj, "XTrackQualityFlags", [iline,0], &
           [nloop,nxtrack], omi_xtrflg_l1b(1:nxtrack,0:nloop-1), errstat)
    else
      !No need to set omi_auraalt as it is only used when writing HE5
      omi_xtrflg_l1b(1:nxtrack,0:nloop-1) = 0
    endif
    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) return

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
      ! EJOS - TEMPO zenith angles have high fill values which must be 
      ! excluded, so we require a check against a maximum legal value
      ! ---------------------------------------------------------------------
      WHERE ((omi_szenith(1:nxtrack,iloop) < min_zenith ) &
           .or. (omi_szenith(1:nxtrack,iloop) > max_legal_zenith ))
        omi_szenith(1:nxtrack,iloop) = r4_missval
      ENDWHERE
      WHERE ((omi_vzenith(1:nxtrack,iloop) < min_zenith ) &
           .or. (omi_vzenith(1:nxtrack,iloop) > max_legal_zenith ))
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
    use tg_names_module
    use netcdf, only : nf90_nowrite

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
    INTEGER (KIND=i4), DIMENSION (nx,0:nt-1) :: geoflg
    !type (L1B_Object_Type) :: l1bobj
    type (tiof_file_type) :: tio_l1obj

    ! let errstat flow
    !call l1bread_open_swath (l1bfile, omi_radiance_swathname, l1bobj, errstat)
    !call l1bread_get2d_i2 (l1bobj, "GroundPixelQualityFlags", 0, nt, geoflg, errstat)
    !call l1bread_close (l1bobj)
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_l1obj, omi_radiance_swathname, errstat)
    call tiof_get2d_ui4  (tio_l1obj, tg_var_gpqf, [0,0], [nt,nx], geoflg(1:nx,0:nt-1), errstat)
    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) return

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
    glint_flg = int (iand (ishft(geoflg, -4), 1_i4), kind=i2)

    ! Bits 8-14 are snow/ice
    snow_ice_flg  = int (iand (ishft(geoflg, -8), 127_i4), kind=i2)

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

!unused   SUBROUTINE omi_xtract_swathname ( l1bfile, l1bchan, omiswath )
!unused 
!unused     USE OMSAO_precision_module
!unused 
!unused     USE OMSAO_parameters_module, ONLY : MAX_STR_LEN
!unused     USE hdfeos4_parameters
!unused     USE strutils
!unused 
!unused     IMPLICIT NONE
!unused 
!unused     ! --------------
!unused     ! Input Variable
!unused     ! --------------
!unused     CHARACTER (LEN=*), INTENT (IN) :: l1bfile
!unused     CHARACTER (LEN=3), INTENT (IN) :: l1bchan
!unused 
!unused     ! ----------------
!unused     ! Output variables
!unused     ! ----------------
!unused     CHARACTER (LEN=*), INTENT (OUT) :: omiswath
!unused 
!unused     ! ---------------
!unused     ! Local variables
!unused     ! ---------------
!unused     INTEGER   (KIND=i4)      :: is, ie
!unused     INTEGER   (KIND=i4)      :: swfid, nswath, strbufsize, xswath
!unused     CHARACTER (LEN=MAX_STR_LEN) :: swathlist
!unused 
!unused     ! ------------------------------
!unused     ! Name of this module/subroutine
!unused     ! ------------------------------
!unused     !CHARACTER (LEN=20), PARAMETER :: modulename = 'omi_xtract_swathname'
!unused 
!unused     ! --------------------------
!unused     ! Initialize OUTPUT variable
!unused     ! --------------------------
!unused     omiswath = '?'
!unused 
!unused     ! ---------------------------------------------------------
!unused     ! Inquire about the swaths in the current L1b radiance file
!unused     ! ---------------------------------------------------------
!unused     swfid  = SWOpen     ( l1bfile, DFACC_READ )
!unused     swathlist="" !JED
!unused     nswath = SWInqswath ( l1bfile, swathlist, strbufsize )
!unused     xswath = SWClose    ( swfid )
!unused 
!unused     ! --------------------------------------------------------------------
!unused     ! Extract the swath name we need. Either there is one one swath in the
!unused     ! file (VIS) or there are two (UV-1, UV-2)
!unused     ! --------------------------------------------------------------------
!unused     is = 1
!unused     SELECT CASE ( nswath )
!unused     CASE ( 1 )
!unused       is = 1 ; ie = strbufsize
!unused     CASE ( 2 )
!unused       CALL find_endstring ( strbufsize, swathlist, 1, ie )
!unused       SELECT CASE ( l1bchan )
!unused       CASE ( 'UV1' )
!unused         is = 1 ; ie = ie-1
!unused       CASE ( 'UV2' )
!unused         is = ie+1 ; ie = strbufsize
!unused       CASE DEFAULT
!unused         ! Nothing to do here except to fold
!unused       END SELECT
!unused     CASE DEFAULT
!unused       ! Nothing to do here except to fold.
!unused     END SELECT
!unused 
!unused     omiswath = TRIM(ADJUSTL(swathlist(is:ie)))
!unused 
!unused     RETURN
!unused   END SUBROUTINE omi_xtract_swathname

END MODULE
