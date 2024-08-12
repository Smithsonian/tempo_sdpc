!> Read input from TEMPO L1 and L2 netCDF4 files
module m_read_input_tio
  use tio_module
  use tell_module

  public open_tio, close_tio, read_irr_tio, read_rad_tio
  private read_rad_dims, quick_lin_interpol, allocate_rad_vars
  private quick_irr_interpol

contains

  !> Open a netCDF file
  !-----------------------------------------------------------------------
  !> @param[in]   l1file     filename
  !> @param[out]  tio_l1obj  libtio file object
  !> @param       errstat    error handling integer, non-zero = problem
  !
  !> @author E. O'Sullivan April 2021
  !-----------------------------------------------------------------------
  subroutine open_tio(l1file, tio_l1obj, errstat)

    use netcdf, only: nf90_nowrite

    implicit none
    !input variables
    character (len=*), intent(in) :: l1file
    !output variables
    type (tiof_file_type) :: tio_l1obj
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    call tiof_open(l1file, tio_l1obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "open_tio: failed to open "//trim(adjustl(l1file)), errstat)
      return
    endif

  end subroutine open_tio

  !> Get 1 global attribute hqw & gga
  !-----------------------------------------------------------------------
  subroutine get_tio_global_attr(l1file, attrname, attrval, errstat)

    use netcdf, only: nf90_global, nf90_get_att, nf90_nowrite

    implicit none

    character (len=*), intent(in) :: l1file, attrname
    integer (kind=4), intent(inout) :: errstat
    character (len=*), intent(out) :: attrval

    integer (kind=4) :: ncerr
    type (tiof_file_type) :: tio_l1obj

    call tiof_open(l1file, tio_l1obj, nf90_nowrite, errstat)
    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, attrname, attrval)
    call tiof_close(tio_l1obj, errstat)

  end subroutine get_tio_global_attr

  !> Get TEMPO L1 RAD attributes hqw&gga
  !----------------------------------------------------------------------
  !could use get_tio_global_attr, but this one aviods multiple open & close

  subroutine get_tio_l1rad_glbattr(l1file, errstat)
    use m_vars, only: gmetadata

    use netcdf, only: nf90_global, nf90_get_att, nf90_nowrite
    implicit none

    include 'GetConfig.inc'

    character (len=*), intent(in) :: l1file
    integer (kind=4), intent(inout) :: errstat

    character(len=CFG_VAL_LEN) :: attrval
    integer (kind=4) :: ncerr, tmpint

    real (kind=4) :: tmpreal

    type(tiof_file_type) :: tio_l1obj

    call tiof_open(l1file, tio_l1obj, nf90_nowrite, errstat)

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'time_coverage_start', attrval)
    gmetadata%startdate=attrval(01:11)
    gmetadata%starttime=attrval(12:19)

    read(attrval(01:04),*) tmpint
    gmetadata%granule_year = tmpint
    read(attrval(06:07),*) tmpint
    gmetadata%granule_month = tmpint
    read(attrval(09:10),*) tmpint
    gmetadata%granule_day = tmpint
    read(attrval(12:13),*) tmpint
    gmetadata%granule_hour_start = tmpint
    read(attrval(15:16),*) tmpint
    gmetadata%granule_minute_start = tmpint
    read(attrval(18:19),*) tmpint
    gmetadata%granule_seconds_start = tmpint

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'time_coverage_end', attrval)
    gmetadata%enddate = attrval(01:11)
    gmetadata%endtime = attrval(12:19)
    read(attrval(12:13), *) tmpint
    gmetadata%granule_hour_end = tmpint
    read(attrval(15:16), *) tmpint
    gmetadata%granule_minute_end = tmpint
    read(attrval(18:19), *) tmpint
    gmetadata%granule_seconds_end = tmpint

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'scan_num',tmpint)
    gmetadata%scan_num = tmpint

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'granule_num',tmpint)
    gmetadata%granule_num = tmpint

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'geospatial_lon_min',tmpreal)
    gmetadata%geospatial_lon_min = tmpreal

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'geospatial_lon_max',tmpreal)
    gmetadata%geospatial_lon_max = tmpreal

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'geospatial_lat_min',tmpreal)
    gmetadata%geospatial_lat_min = tmpreal

    ncerr = nf90_get_att(tio_l1obj%fileid, nf90_global, 'geospatial_lat_max',tmpreal)
    gmetadata%geospatial_lat_max = tmpreal

    call tiof_close(tio_l1obj, errstat)

  end subroutine get_tio_l1rad_glbattr

  !> Close a netCDF file
  !-----------------------------------------------------------------------
  !> @param[in]  tio_l1obj  libtio file object
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !> @author E. O'Sullivan April 2021
  !-----------------------------------------------------------------------
  subroutine close_tio (tio_l1obj, errstat)

    implicit none
    !input variables
    type (tiof_file_type) :: tio_l1obj
    !output variables
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    call tiof_close(tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "close_tio: failed", errstat)
      return
    endif

  end subroutine close_tio

  !> Read the dimensions of a L1 radiance or irradiance file
  !-----------------------------------------------------------------------
  !
  !> @param[in]  swathname  name of group from which to read dimensions
  !> @param[in]  tio_l1obj  libtio file object
  !> @param[in]  nxtrack    cross-track dimension
  !> @param[in]  ntimes     along-track dimension
  !> @param[in]  nwavel     wavelength dimension
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !> @author E. O'Sullivan April 2021
  !-----------------------------------------------------------------------
  subroutine read_rad_dims(tio_l1obj, swathname, nxtrack, ntimes, nwavel, &
       errstat)

    implicit none
    !input variables
    character (len=*), intent(in) :: swathname
    type (tiof_file_type) :: tio_l1obj
    !output variables
    integer (kind=4), intent(out) :: nxtrack, ntimes, nwavel
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    call tiof_push_group (tio_l1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, "xtrack", nXtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "mirror_step", nTimes, errstat)
    call tiof_inq_dimlen (tio_l1obj, "spectral_channel", nWavel, errstat)
    call tiof_pop_group (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_rad_dims: failed", errstat)
      return
    endif

  end subroutine read_rad_dims

  !> Read L1 irradiance file and build 440, 466, 477nm irradiance arrays
  !-----------------------------------------------------------------------
  !
  !> @param[in]  l1_file    TEMPO-format irradiance filename
  !> @param[in]  swathname  name of group from which to read dimensions
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !> @author E. O'Sullivan April 2021
  !-----------------------------------------------------------------------
  subroutine read_irr_tio (l1_file, swathname, errstat)

    use m_vars, only: irr_out_irradiance_440nm, irr_out_irradiance_466nm, &
         irr_out_irradiance_477nm, irr_EarthSunDist, w440, w466, w477

    use m_vars, only: irr_NumTimes, irr_nXtrack, irr_nWavel

    use m_vars, only: rad_NumTimes, rad_nXtrack, out_ProcessingQualityFlags

    use m_vars, only: irr_waveshift, option_apply_solshift

    use m_vars, only: negative999

    use m_vars, only: ixdebug, run_mode, lun_debug_irr

    implicit none

    !input variables
    character (len=*), intent(in) :: l1_file, swathname
    !output variables
    integer (kind=4), intent(inout) :: errstat
    !local variables
    integer (kind=4) :: nxtrack, ntimes, nwavel, ix, it, iw
    real (kind=4), dimension(:,:,:), allocatable :: tio_irr, tio_wvl
    real (kind=4), dimension(:), allocatable :: tmpwvl
    integer (kind=2), dimension(:,:,:), allocatable :: tio_pqf
    real (kind=4) :: thisirr440, thisirr466, thisirr477

    character(len=80) :: logmsg
    integer (kind=4) :: lerr

    type(tiof_file_type) :: tio_l1obj

    ! only execute if no previous error
    if (errstat /= 0) return

    ! open file, get dimensions, allocate arrays
    call open_tio (l1_file, tio_l1obj, errstat)
    call read_rad_dims (tio_l1obj, swathname, nxtrack, ntimes, nwavel, errstat)
    if (errstat /= 0) return

    irr_NumTimes = ntimes
    irr_nXtrack = nxtrack
    irr_nWavel = nwavel
    write(*,*)'irr: ntimes,nxtrack,nwavel=',irr_NumTimes,&
         irr_nXtrack,irr_nWavel

    ! sanity check, OMCDO2N reads RAD before IRR
    if (irr_nXtrack .NE. rad_nXtrack) then
         write(*,*) 'irr_nXtrack and rad_nXtrack differ'
         errstat = -1
    endif

    allocate (irr_out_irradiance_440nm(nxtrack), &
         irr_out_irradiance_466nm(nxtrack), &
         irr_out_irradiance_477nm(nxtrack), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_irr_tio: allocation failed",&
           errstat)
      return
    endif

    ! initialize 
    irr_out_irradiance_440nm = negative999
    irr_out_irradiance_466nm = negative999
    irr_out_irradiance_477nm = negative999

    allocate (tio_irr(nwavel, nxtrack, 1), &
              tio_wvl(nwavel, nxtrack, 1), &
              tio_pqf(nwavel, nxtrack, 1), &
              stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_irr_tio: allocation failed",&
           errstat)
      return
    endif

    !read earth_sun_distance
    call tiof_get_r4 (tio_l1obj, "earth_sun_distance", &
         irr_EarthSunDist, errstat)

    !read wavelength, irradiance, quality flag
    call tiof_push_group (tio_l1obj, swathname, errstat)
    call tiof_get3d_r4 (tio_l1obj, "irradiance", [0,0,0], &
         [1,nXtrack,nWavel], tio_irr, errstat)
    call tiof_get3d_r4 (tio_l1obj, "wavelength", [0,0,0], &
         [1,nXtrack,nWavel], tio_wvl, errstat)
    call tiof_get3d_i2 (tio_l1obj, "pixel_quality_flag", [0,0,0], &
         [1,nXtrack,nWavel], tio_pqf, errstat)
    call close_tio (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_irr_tio: failed to read irr", errstat)
      return
    endif

    ! apply solar wavelength shift
    if (option_apply_solshift == 1) then
       write(*,*) 'applying solar wavelength shift'
       allocate (tmpwvl(nwavel),stat=errstat) 
       do ix = 1, nxtrack
          tmpwvl(:) = tio_wvl(:,ix,1) - irr_waveshift(ix)
          tio_wvl(:,ix,1) = tmpwvl(:)
       enddo       
       ! tio_wvl now agrees with solar reference spectrum in fitting
       deallocate(tmpwvl)
    endif

    ! interpolate values at 440, 466, 477nm
    ! problem with pixel_quality_flag, set to -9999.0 in quick_lin_interpol
     
    do ix = 1, nxtrack ! loop over pixels
      ! lerr is local error for each case, donot affect global errstat
      ! failed interpol will have negative999 returned
      thisirr440 = negative999 
      call quick_lin_interpol (tio_wvl(:,ix,1), w440, tio_irr(:,ix,1), &
           thisirr440, tio_pqf(:,ix,1),nwavel, negative999, lerr)
      ! 440nm does not affect ecf
      if (lerr /= 0) then
        write (logmsg,'(A,I4)') "440nm irr interpol failed for cross-track ",ix
        call tell_log (1, logmsg)
      endif
      if (thisirr440 .GT. 0.) then
          irr_out_irradiance_440nm(ix) = thisirr440
      else
          irr_out_irradiance_440nm(ix) = negative999
          do it = 1, rad_NumTimes
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
          end do 
      endif

      thisirr466 = negative999
      call quick_lin_interpol (tio_wvl(:,ix,1), w466, tio_irr(:,ix,1), &
           thisirr466, tio_pqf(:,ix,1),nwavel,negative999,lerr)
      ! the following is for testing
      !call quick_irr_interpol(tio_wvl(:,ix,1), w466, tio_irr(:,ix,1), &
      !      thisirr466, tio_pqf(:,ix,1), negative999,lerr)
      if (lerr /= 0) then
        write (logmsg,'(A,I4)') "466nm interpol failed for cross-track ",ix
        call tell_log (1, logmsg)
      endif
      if (thisirr466 .GT. 0.) then
           irr_out_irradiance_466nm(ix) = thisirr466
      else
           irr_out_irradiance_466nm(ix) = negative999
           do it = 1, rad_NumTimes
              out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it), 8)
           end do
      endif

      ! irr_out_irradiance_477nm does not affect out_ProcessingQualityFlags
      thisirr477 = negative999
      call quick_lin_interpol (tio_wvl(:,ix,1), w477, tio_irr(:,ix,1), &
           thisirr477, tio_pqf(:,ix,1),nwavel,negative999,lerr)
      if (lerr /= 0) then
        write (logmsg,'(A,I4)') "477nm interpol failed for cross-track ",ix
        call tell_log (1, logmsg)
      endif
      if (thisirr477 .GT. 0.) then
           irr_out_irradiance_477nm(ix) = thisirr477
      else
           irr_out_irradiance_477nm(ix) = negative999
      endif
      
    enddo ! end of ix loop

    ! debug
    if ((trim(run_mode) .eq. 'development').and.(ixdebug .ge. 0)) then
       write(*,*) 'writing debug_irr.txt'
       open(unit=lun_debug_irr, file='debug_irr.txt')
       write(lun_debug_irr,*)'irr_EarthSunDist=',irr_EarthSunDist
       write(lun_debug_irr,*)'ixdebug=',ixdebug
       write(lun_debug_irr,*)'iw   wavelength    irr'
       do iw = 1, nwavel
          write(lun_debug_irr,*) iw, tio_wvl(iw,ixdebug,1),tio_irr(iw,ixdebug,1)
       enddo
       
       write(lun_debug_irr,*) 'ix   IRR440   IRR466   IRR477'
       ix = ixdebug
       !do ix = 1, nxtrack
       write(lun_debug_irr,*) ix, irr_out_irradiance_440nm(ix),&
          irr_out_irradiance_466nm(ix),irr_out_irradiance_477nm(ix)
       !enddo
    endif

    close(lun_debug_irr)

    ! deallocate large local arrays after use
    deallocate (tio_wvl, tio_irr, tio_pqf, stat=errstat)

  end subroutine read_irr_tio

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_rad_tio (l1_file, swathname, errstat)

    use m_vars, only: rad_Time, rad_Latitude, rad_Longitude, &
         rad_SolarZenithAngle, rad_ViewingZenithAngle, &
         rad_ViewingAzimuthAngle, rad_SolarAzimuthAngle, & 
         out_TerrainHeight, scddes_hour,&
         rad_GroundPixelQualityFlags, &
         out_ProcessingQualityFlags, &
         w440, w466, w477, fFillValue,&
         rad_440nm,rad_466nm,rad_477nm, rad_EarthSunDist, &
         rad_NumTimes, rad_nXtrack, rad_nWavel, rad_SnowIceFraction

    use m_vars, only: rad_waveshift, option_apply_radshift

    use m_vars, only: ixdebug, itdebug, run_mode, lun_debug_rad

    use m_vars, only: negative999, dFillValue

    !hqw added rad_440nm,rad_466nm, rad_477nm in m_vars
    implicit none

    !input variables
    character (len=*), intent(in) :: l1_file, swathname
    !output variables
    integer (kind=4), intent(inout) :: errstat

    !local variables
    type(tiof_file_type) :: tio_l1obj

    integer (kind=4) :: ntimes, nxtrack, nwavel, ix, it
    integer :: iw1, iw2, iw, nw
    real(kind=4) :: rad466, rad477, rad440
    real(kind=4) :: ww1, ww2, www, yy1, yy2
    real(kind=4) :: min_rad, delta_w
  
    ! moved rad_Radiance, rad_Wavelength from m_vars.f90 here
    ! local variables save memory as they are deallocated after use
    real(kind=4), dimension(:,:,:), allocatable:: rad_Radiance
    real(kind=4), dimension(:,:,:), allocatable:: rad_Wavelength
    integer(kind=2), dimension(:,:,:), allocatable:: rad_PixelQualityFlags
    integer(kind=2), dimension(:,:), allocatable:: rad_TerrainHeight

    real(kind=4), dimension(:), allocatable:: temp_wav, temp_rad
    integer(kind=2) :: temp_radflags

    ! only execute if there is no previous error
    if (errstat /= 0) return

    ! min radiance, anything smaller is treated as invalid
    min_rad = 0. ! 1.e9
    ! maximum wavelength difference around w440, w466, w477
    ! to prevent drifting too far from the target due to unusable wavelengths
    delta_w = 0.5 ! nm

    ! Open file, get dimensions
    call open_tio (l1_file, tio_l1obj, errstat)
    call read_rad_dims (tio_l1obj, swathname, nxtrack, ntimes, nwavel, errstat)
    if (errstat /= 0) return

    rad_NumTimes = ntimes
    rad_nXtrack = nxtrack
    rad_nWavel = nwavel
    write(*,*) 'read_rad_tio:ntimes,nxtrack,nwavel=',&
          rad_NumTimes,rad_nXtrack,rad_nWavel

    !allocate m_vars arrays
    call allocate_rad_vars (ntimes, nxtrack, errstat)
    call allocate_init_desvars (nxtrack, errstat)
    if (errstat /= 0) return

    !allocate local arrays & initilize
    allocate(rad_Radiance(nwavel, nxtrack, ntimes), stat=errstat)
    allocate(rad_Wavelength(nwavel, nxtrack, ntimes), stat=errstat)
    allocate(rad_PixelQualityFlags(nwavel, nxtrack, ntimes), stat=errstat)
    allocate(rad_TerrainHeight(nxtrack, ntimes), stat=errstat)
    rad_Radiance = fFillValue
    rad_Wavelength = fFillValue
    rad_PixelQualityFlags = 0
    rad_TerrainHeight = 0

    ! read the arrays
    call tiof_use_file_epoch (tio_l1obj, errstat)
    call tiof_get_r4 (tio_l1obj, "earth_sun_distance", rad_EarthSunDist, &
         errstat)
    call tiof_get1d_r8 (tio_l1obj, "time", [0], [ntimes], rad_Time, errstat, &
         replace_fill=dFillValue)
    call tiof_push_group (tio_l1obj, swathname, errstat)
    call tiof_get2d_r4 (tio_l1obj, "longitude", [0,0], [ntimes, nxtrack], &
         rad_Longitude, errstat)
    call tiof_get2d_r4 (tio_l1obj, "latitude", [0,0], [ntimes, nxtrack], &
         rad_Latitude, errstat)
    call tiof_get2d_r4 (tio_l1obj, "solar_zenith_angle", [0,0], &
         [ntimes, nxtrack], rad_SolarZenithAngle, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_zenith_angle", [0,0], &
         [ntimes, nxtrack], rad_ViewingZenithAngle, errstat)
   ! what needed is relative azimuth angle read in by read_cldo4_tio
   ! SAA and VAA are read as original code uses them to calculate RAA
    call tiof_get2d_r4 (tio_l1obj, "solar_azimuth_angle", [0,0], &
         [ntimes, nxtrack], rad_SolarAzimuthAngle, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_azimuth_angle", [0,0], &
         [ntimes, nxtrack], rad_ViewingAzimuthAngle, errstat)

    ! added rad_SnowIceFraction & initialize to zero
    rad_SnowIceFraction=0.
    call tiof_get2d_r4 (tio_l1obj, "snow_ice_fraction", [0,0], &
         [ntimes, nxtrack], rad_SnowIceFraction, errstat)

    ! read L1b terrain_height as local variable rad_TerrainHeight
    ! tranfer int rad_TerrainHeight to real out_TerrainHeight
    call tiof_get2d_i2 (tio_l1obj, "terrain_height", [0,0], &
         [ntimes,nxtrack], rad_TerrainHeight, errstat)
    out_TerrainHeight = real(rad_TerrainHeight, kind=4)
    where (rad_TerrainHeight < -1000)
        out_TerrainHeight = fFillValue
    endwhere

    call tiof_get2d_ui4 (tio_l1obj, "ground_pixel_quality_flag", [0,0], &
         [ntimes,nxtrack], rad_GroundPixelQualityFlags, errstat)
    call tiof_get3d_ui2 (tio_l1obj, "pixel_quality_flag", [0,0,0], &
         [ntimes,nxtrack,nwavel], rad_PixelQualityFlags, errstat)
    call tiof_get3d_r4 (tio_l1obj, "radiance", [0,0,0], &
         [ntimes,nxtrack,nwavel], rad_Radiance, errstat)
    call tiof_get3d_r4 (tio_l1obj, "wavelength", [0,0,0], &
         [ntimes,nxtrack,nwavel], rad_Wavelength, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_rad_tio: failed", errstat)
      return
    endif

    call close_tio (tio_l1obj, errstat)

   ! moved interpolation to rad440, rad466 and rad440 here
   ! from cal_ecf.f90 so that big arrays
   ! rad_Wavelength, radRadiance, rad_PixelQualityFlags can be deallocated

    allocate(temp_wav(nwavel), temp_rad(nwavel), stat = errstat)

   ! apply rad wavelength shift
   if (option_apply_radshift .eq. 1) then 
      write(*,*) 'applying rad wavelength shift'
      do it = 1, ntimes
         do ix = 1, nxtrack
            temp_wav(:) = rad_Wavelength(:,ix,it)-rad_waveshift(ix,it)
            rad_Wavelength(:,ix,it) = temp_wav(:)
         enddo
      enddo 
      ! rad_Wavelength now agrees with solar reference used in fitting
   endif

   ! debug
   if ((trim(run_mode) .eq. 'development').and. &
       (ixdebug .ge. 0) .and. (itdebug .ge. 0)) then
       write(*,*) 'writing debug_rad.txt'
       open(unit=lun_debug_rad, file='debug_rad.txt')
       write(lun_debug_rad,*)'rad_EarthSunDist=',rad_EarthSunDist
       write(lun_debug_rad,*) 'iw   lamda     rad '
   endif

   nw = nwavel

   ! loop through pixels
   do it = 1, ntimes
      do ix = 1, nxtrack
        ! get local spectrum, if PixelQualiyFlags indicate spectral issue
        ! set the corresponding temp_rad to negative value -9999.
        ! negative temp_rad will not be used in subsequent interpolation
        do iw = 1, nw
         temp_wav(iw) = rad_Wavelength(iw,ix,it)
         temp_radflags = rad_PixelQualityFlags(iw,ix,it)
         ! all bit = 0 is fast but may be too restrictive, 
         !if (rad_PixelQualityFlags(iw,ix,it) .NE. 0) then
         !  switch back to bit test if desired 
          if (btest(temp_radflags,0) .or. & !missing
              btest(temp_radflags,1) .or. & !bad
              btest(temp_radflags,2) .or. & !processing_error
              btest(temp_radflags,5) .or. & !saturated
         !     btest(temp_radflags,6) .or. & !noise_underflow
              btest(temp_radflags,7) .or. & !dark_corr_error
              btest(temp_radflags,8) .or. & !offset_corr_error
              btest(temp_radflags,9) .or. & !smear_corr_error
         !     btest(temp_radflags,11) .or. &  !nonlinear_range
              btest(temp_radflags,10)) then  !straylight_corr_error
          temp_rad(iw)=negative999 
         else
          temp_rad(iw)=rad_Radiance(iw,ix,it)
         end if
        enddo !iw

        ! debug
        if ((trim(run_mode).eq.'development').and. &
            (ix.eq.ixdebug).and.(it.eq.itdebug)) then 
            do iw = 1, nw
               write(lun_debug_rad,*) iw, rad_Wavelength(iw,ix,it), rad_Radiance(iw,ix,it)
            enddo
        endif

      ! interpolation could be done using function call, however,
      ! explicit code is faster as this is done for each pixel
      ! ----------------------------
      ! calculate radiance at 466 nm
      ! ----------------------------
      ! find the wavelength index nearest 466 nm
      rad466=negative999 
      iw1=-9
      iw2=-9
      do iw=1,nw
        if((w466 .lt. temp_wav(iw)) .and. (temp_rad(iw) .gt. min_rad)) then
          iw2=iw
          exit
        endif
      end do

      do iw=nw,1,-1
        if((w466 .ge. temp_wav(iw)) .and. (temp_rad(iw) .gt. min_rad)) then
          iw1=iw
          exit
        endif
      end do

      !when iw1 & iw2 >=0, yy1 and yy2 should be valid due to above
      if((iw1 .ge. 0) .and. (iw2 .ge. iw1)) then
        yy1=temp_rad(iw1)
        yy2=temp_rad(iw2)
        ww1=w466-temp_wav(iw1)
        ww2=temp_wav(iw2)-w466
        www = ww1 + ww2
        if ((iw2 .gt. iw1) .and. (www .lt. delta_w) .and. (www .gt. 0.)) then 
       ! should be the case
           rad466=(ww1*yy2+ww2*yy1)/www
        else ! just for safeguard
           rad466=yy1
        endif
      else
        rad466=negative999 
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),8)
      endif

      rad_466nm(ix,it) = rad466
      ! ----------------------------
      ! calculate radiance at 477 nm
      ! ----------------------------
      rad477=negative999 
      iw1=-9
      iw2=-9
      do iw=1,nw
        if((w477 .lt. temp_wav(iw)) .and. (temp_rad(iw) .gt. min_rad)) then
          iw2=iw
          exit
        endif
      end do

      do iw=nw,1,-1
        if((w477 .ge. temp_wav(iw)) .and. (temp_rad(iw) .gt. min_rad)) then
          iw1=iw
          exit
        endif
      end do

      if((iw1 .ge. 0) .and. (iw2 .ge. iw1)) then
        yy1=temp_rad(iw1)
        yy2=temp_rad(iw2)
        ww1=w477-temp_wav(iw1)
        ww2=temp_wav(iw2)-w477
        www = ww1 + ww2
        if ((iw2 .gt. iw1) .and. (www .lt. delta_w) .and. (www .gt. 0.)) then
           rad477=(ww1*yy2+ww2*yy1)/www
        else
           rad477=yy1
        endif
      else
        ! rad477 is not actually used for calculation, only a diagnostic
        ! thus did not assign ProcessingQualityFlag for it
        rad477=negative999 
      endif
      rad_477nm(ix,it) = rad477

      ! ----------------------------
      ! calculate radiance at 440 nm
      ! ----------------------------
      rad440=negative999 
      iw1=-9
      iw2=-9
      do iw=1,nw
        if((w440 .lt. temp_wav(iw)) .and. (temp_rad(iw) .gt. min_rad)) then
          iw2=iw
          exit
        endif
      end do

      do iw=nw,1,-1
        if((w440 .ge. temp_wav(iw)) .and. (temp_rad(iw) .gt. min_rad)) then
          iw1=iw
          exit
        endif
      end do

      if((iw1 .ge. 0) .and. (iw2 .ge. iw1)) then
        yy1=temp_rad(iw1)
        yy2=temp_rad(iw2)
        ww1=w440-temp_wav(iw1)
        ww2=temp_wav(iw2)-w440
        www = ww1 + ww2
        if ((iw2 .gt. iw1) .and. (www .lt. delta_w) .and. (www .gt. 0.)) then 
           rad440=(ww1*yy2+ww2*yy1)/www
        else
           rad440=yy1
        endif
      else
        rad440=negative999 
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
      endif
      rad_440nm(ix,it) = rad440

      enddo !ix
    enddo !it

    ! debug 
    if ((trim(run_mode) .eq. 'development').and. &
        (itdebug .ge. 0) .and. (ixdebug .ge. 0)) then
       write(lun_debug_rad,*)'ixdebug,itdebug=',ixdebug,itdebug
       write(lun_debug_rad,*)'rad440, rad466, rad477=',rad_440nm(ixdebug,itdebug),&
          rad_466nm(ixdebug,itdebug),rad_477nm(ixdebug,itdebug)
    endif

    close(lun_debug_rad)

    ! deallocate large temporary rad arrays
    deallocate(rad_Radiance, rad_Wavelength, rad_PixelQualityFlags)
    deallocate(rad_TerrainHeight)
    deallocate(temp_rad, temp_wav)

  end subroutine read_rad_tio

!-------------------------------------------------------------------
  subroutine read_cldo4_tio (l2_file, errstat)
     use m_vars, only: nasa_SlantColumnAmountO2O2,fFillValue
     use m_vars, only: nasa_NumTimes, nasa_nXtrack, run_mode
     use m_vars, only: scd_mdqfl,nasa_scdrms,nasa_scduncertainty
     use m_vars, only: rad_RelativeAzimuthAngle, out_RelativeAzimuthAngle
     use m_vars, only: max_o2o2_scd, max_o2o2_uncertainty, max_o2o2_relerr
     use m_vars, only: o2o2_norm
     use m_vars, only: fFillValue, dFillValue, option_scdfullfilter

     implicit none

     !input variables
     character (len=*), intent(in) :: l2_file

     !output variables
     integer (kind=4), intent(inout) :: errstat

     !local variables
     type(tiof_file_type) :: tio_l2obj
     integer (kind=4) :: ntimes, nxtrack, ix, it
     real (kind=4), allocatable, dimension(:,:) :: scd_relerr
     real (kind=8), allocatable, dimension(:,:) :: tmp_dbl
     real (kind=8):: dspecial

     real(kind=4):: tmp_raa, temp_raa, fspecial

     integer :: errstat1 

     if (errstat /= 0) return

     fspecial = fFillValue
     dspecial = dFillValue !-9999.d0
     errstat1 = 0

     !Open file, get dimensions
     call open_tio (l2_file, tio_l2obj, errstat)

     ! get information from support_data group
     call read_cldo4_dims(tio_l2obj, 'support_data', nXtrack, nTimes, errstat)
     if (errstat /= 0) then
          call tell_error (tell_runtime_error, "read_cldo4_dims: failed", errstat)
          return
     endif
     nasa_nXtrack = nXtrack
     nasa_numTimes = nTimes
     write(*,*) 'read_cldo4_tio:nXtrack,nTimes=',nXtrack,nTimes

     ! Allocate cldo4 variables
     call allocate_cldo4_vars(nTimes,nXtrack,errstat)
     if (errstat /= 0) then
          call tell_error (tell_runtime_error, "allocate_cldo4_vars: failed", errstat)
          return
     endif

     ! read main_data_quality_flag or SCD_MainDataQualityFlags, 
     ! errstat1 suggest which one to read
     ! try product group to read main_data_quality_flag
     call tiof_push_group(tio_l2obj, "product", errstat)

     call tiof_get2d_i2(tio_l2obj,"main_data_quality_flag", [0,0], [ntimes, nxtrack],&
            scd_mdqfl, errstat1)

     if (errstat1 /=0) then
          write(*,*)'-->no main_data_quality_flag found in product group'
          ! get out of product group
          call tiof_pop_group (tio_l2obj, errstat)
          ! try support_data group to read SCD_MainDataQualityFlags
          call tiof_push_group(tio_l2obj, "support_data", errstat) 
          call tiof_get2d_i2(tio_l2obj,"SCD_MainDataQualityFlags",[0,0], &
                   [ntimes,nxtrack],scd_mdqfl, errstat)
         if (errstat /=0) then
          call tell_error(tell_runtime_error, "read_cldo4_tio: scd_mdqfl fail",errstat)
          return
         endif  
         write(*,*)' -->SCD_MainDataQualityFlags read from support_data group'
         ! get out of support_data group
         call tiof_pop_group(tio_l2obj, errstat) 
     else
        write(*,*)'-->main_data_quality_flag read from product group.'
        ! scd_mdqfl read, get out of product group
        call tiof_pop_group (tio_l2obj, errstat)
     endif

     ! allocate tmp_dbl array 
     allocate(tmp_dbl(nXtrack, nTimes), stat = errstat)
     if (errstat /= 0) then
          call tell_error (tell_runtime_error, &
                "allocate tmp_dbl: failed in read_cldo4_tio", errstat)
          return
     endif

     ! open support_data group to read SCD
     call tiof_push_group (tio_l2obj,"support_data", errstat)
     if (errstat /= 0) then
          call tell_error (tell_io_read_error, "read_cldo4_tio: pushing support_data group failed", errstat)
          return
     endif

     ! read fitted_slant_column
     call tiof_get2d_r8 (tio_l2obj, "fitted_slant_column", [0,0], [ntimes, nxtrack],&
           tmp_dbl, errstat)
     if (errstat /= 0) then
          call tell_error (tell_runtime_error, "read_cldo4_tio: fitted_slant_column failed", errstat)
          return
     endif

     ! normalize scd 
     tmp_dbl = tmp_dbl/o2o2_norm

     ! filter out unphysical values
     ! note: mdqfl=0 (normal), mdqfl=1 (suspicious), mdqfl=2 (bad) 
     ! for o2o2, suspicious scds are extremely large, they should not be used
     ! any scd<0. will be skipped for ocp and pscene calculation,
     where ((tmp_dbl < 0.) .or. (tmp_dbl > max_o2o2_scd) .or. &
            (scd_mdqfl .ne. 0)) 
           tmp_dbl = dspecial ! negative value
     end where

     ! assign nasa_SlantColumnAmountO2O2 
     ! de-stripe scd, nasa_SlantColumnAmountO2O2 is real (kind=4)
     call destripe_o2o2scd(tmp_dbl,nXtrack,nTimes,nasa_SlantColumnAmountO2O2)

    ! set filtered out scd to fFillValue for output
    where (nasa_SlantColumnAmountO2O2 < 0.)
          nasa_SlantColumnAmountO2O2 = fFillValue
    endwhere

     ! ------
     ! read fitted SCD uncertainty
     call tiof_get2d_r8(tio_l2obj, "fitted_slant_column_uncertainty",[0,0],&
              [ntimes, nxtrack], tmp_dbl, errstat)
     if (errstat /=0) then
         call tell_error(tell_runtime_error,"read_cldo4_tio: failed scduncertainty", errstat)
         return
     endif
     ! normalize scd uncertainty and assign nasa_scduncertainty
     ! filter out large uncertainty, use max_o2o2_scd below to be permissive
     tmp_dbl = tmp_dbl/o2o2_norm
     where((tmp_dbl < 0.).or.(tmp_dbl > max_o2o2_scd).or. &
           (scd_mdqfl .ne. 0))
           tmp_dbl = dspecial
     endwhere

     ! assign nasa_scduncertainty 
     nasa_scduncertainty = real(tmp_dbl,kind=4)
     where(tmp_dbl < 0.)
           nasa_scduncertainty = fFillValue
     endwhere

     ! deallocate tmp_dbl array
     deallocate(tmp_dbl)

     ! Get out of support_data group
     call tiof_pop_group(tio_l2obj, errstat)
     if (errstat /= 0) then
          call tell_error (tell_io_read_error, "read_cldo4_tio: surport_data failed", errstat)
          return
     endif

    !-------------
    ! stricter filtering according to relative & abs uncertainty
    if (option_scdfullfilter .eq. 1) then

     ! calculate relative scd error
     allocate(scd_relerr(nXtrack, nTimes), stat = errstat)
     if (errstat /= 0) then
          call tell_error (tell_runtime_error, "allocate scd_relerr: failed in read_cldo4_tio", errstat)
          return
     endif

     ! nasa_scduncertainty & nasa_SlantColumnAmountO2O2 are both normalized
     scd_relerr = nasa_scduncertainty / nasa_SlantColumnAmountO2O2
     ! assign huge positive value for negative scd or scduncertainty
     ! so that they will be filtered out 
     where (nasa_scduncertainty < 0.) 
           scd_relerr = 9999. !make it a large POSITIVE value
     endwhere
     where (nasa_SlantColumnAmountO2O2 < 0.)
           scd_relerr = 9999.
     endwhere
    ! bad scd or bad scduncertainty will satisfy the following
    ! and therefore be filtered out
    where (scd_relerr > max_o2o2_relerr)
          nasa_SlantColumnAmountO2O2 = fFillValue
    endwhere

    ! filter out large abs uncertainty
    where (nasa_scduncertainty > max_o2o2_uncertainty)
          nasa_SlantColumnAmountO2O2 = fFillValue
    endwhere

    ! scd_relerr no longer needed
    deallocate(scd_relerr)

    endif ! option_scdfullfilter

    !-------------------------------------------------
    ! now read rad_RelativeAzimuthAngle here
    ! open geolocation group to read raa
    ! other angles are read from L1B file
    call tiof_push_group(tio_l2obj,"geolocation", errstat)

    call tiof_get2d_r4 (tio_l2obj, "relative_azimuth_angle", [0,0], &
         [ntimes, nxtrack], rad_RelativeAzimuthAngle, errstat)

    call tiof_pop_group(tio_l2obj, errstat)

     if (run_mode .NE. 'production') then
        ! open qa_statistics group to read fit_rms_residual
        call tiof_push_group(tio_l2obj, "qa_statistics", errstat)
        call tiof_get2d_r4(tio_l2obj, "fit_rms_residual", [0,0], [ntimes, nXtrack], &
              nasa_scdrms, errstat)
         if (errstat /= 0) then 
              call tell_error(tell_runtime_error,"read_cldo4_tio: failed scdrms", errstat) 
              return
         endif
     endif ! run_mode

     ! Close level 2 file
     call close_tio (tio_l2obj, errstat)

     if (errstat /= 0) then
       call tell_error (tell_io_read_error, "read_cldo4_tio: failed", errstat)
       return
     endif
 
    ! rad_RelativeAzimuthAngle -> out_RelativeAzimuthAngle 
    ! so that RAA is within [0.,180) range for valid pixels
    ! to be consistent with LUT RAA range
    do it = 1, ntimes
      do ix = 1, nxtrack
         tmp_raa = rad_RelativeAzimuthAngle(ix,it)
         ! bad tmp_raa should be a large negative fill value
         if ((tmp_raa .lt. -360.)) then
             out_RelativeAzimuthAngle(ix,it) = fspecial
         else
             ! correct to [0., 360.) range first
             temp_raa = tmp_raa
             do while(temp_raa .lt. 0.)
                temp_raa = temp_raa + 360.
             end do
             do while(temp_raa .ge. 360.)
                temp_raa = temp_raa - 360.
             end do
             
             ! change to LUT RAA range of [0.,180.], using symmetry w.r.t. 0
             if (temp_raa .gt. 180.) temp_raa = 360.-temp_raa
    
             out_RelativeAzimuthAngle(ix,it) = temp_raa
         endif
      end do
    end do
    
   end subroutine read_cldo4_tio

  !> Use simple linear interpolation to find irrad at target wavelength
  !-----------------------------------------------------------------------
  !
  !> @param[in]  w_array    1D wavelength array
  !> @param[in]  w_target   target wavelength
  !> @param[in]  irr_array  1D irradiance array
  !> @param[out] irr_out    irradiance at target wavelength
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !> @author E. O'Sullivan April 2021
  !-----------------------------------------------------------------------
  subroutine quick_irr_interpol (w_array, w_target, irr_array, irr_out, &
       q_array, badout, errstat)
  ! alternative to quick_lin_interpol, not used  
    implicit none

    !input variables
    real (kind=4), dimension(:), intent(in) :: w_array, irr_array
    integer (kind=2), dimension(:), intent(in) :: q_array
    real (kind=4), intent(in) :: w_target
    real (kind=4) :: badout ! value for bad output

    !output variables
    real (kind=4), intent(out) :: irr_out
    integer (kind=4), intent(out) :: errstat ! change inout to out

    !local variables
    integer (kind=4), dimension(1) :: iw1, iw2
    real (kind=4) :: yy1, yy2, ww1, ww2

    errstat = 0 ! initialize

    iw1=maxloc(w_array-w_target, mask=w_array-w_target.lt.0)
    iw2=minloc(w_array-w_target, mask=w_array-w_target.gt.0)
    if (iw2(1) >= iw1(1)) then
      yy1=irr_array(iw1(1))
      yy2=irr_array(iw2(1))
      if ((yy1 .GT. 0.) .and. (yy2 .GT. 0.) .and. &
         (q_array(iw1(1)) .EQ. 0) .and. (q_array(iw2(1)) .EQ. 0)) then 
         ww1=w_target-w_array(iw1(1))
         ww2=w_array(iw2(1))-w_target
         irr_out=(ww1*yy2+ww2*yy1)/(ww1+ww2)
      else
         irr_out = badout !-9999.
      endif
    else
      errstat = -1
    endif

  end subroutine quick_irr_interpol

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Use simple linear interpolation to find irrad at target wavelength
  !-----------------------------------------------------------------------
  !
  !> @param[in]  w_array    1D wavelength array
  !> @param[in]  w_target   target wavelength
  !> @param[in]  irr_array  1D irradiance array
  !> @param[out] irr_out    irradiance at target wavelength
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !-----------------------------------------------------------------------
  subroutine quick_lin_interpol (w_array, w_target, irr_array, irr_out, &
       q_array, nw, badout, errstat)
    implicit none

    !input variables
    real (kind=4), dimension(:), intent(in) :: w_array, irr_array
    integer (kind=2), dimension(:), intent(in) :: q_array
    real (kind=4), intent(in) :: w_target
    real (kind=4) :: badout ! value for bad output

    !output variables
    real (kind=4), intent(out) :: irr_out
    integer (kind=4), intent(in) :: nw
    integer (kind=4), intent(out) :: errstat ! changed inout to out

    !local variables
    integer (kind=4) :: iw1, iw2, iw
    real (kind=4) :: yy1, yy2, ww1, ww2, www, min_val

    ! initial
    min_val = 0.
    irr_out = badout !-9999.

    errstat = 0 ! initialize 
    
    iw1 = -9
    iw2 = -9

    do iw = 1, nw
       if ((w_target .lt. w_array(iw)) .and. (irr_array(iw) .gt. min_val) &
          .and. (q_array(iw) .eq. 0)) then
          iw2=iw
          exit
       endif
    enddo

    do iw = 1, nw
       if ((w_target .ge. w_array(iw)) .and. (irr_array(iw) .gt. min_val) &
          .and. (q_array(iw) .eq. 0)) then
          iw1 = iw
          exit
       endif
    enddo
        
    if ((iw2 .ge. iw1) .and. (iw1 .ge. 0)) then
      yy1=irr_array(iw1)
      yy2=irr_array(iw2)
      ww1=w_target-w_array(iw1)
      ww2=w_array(iw2)-w_target
      www = ww1 + ww2
      if (www .gt. 0.) then
         irr_out=(ww1*yy2+ww2*yy1)/www
      else
         irr_out=badout !-9999.
      endif
    else
      irr_out = badout !-9999.
      errstat = -1
    endif

  end subroutine quick_lin_interpol


  !> Allocate variables associated with radiance
  !-----------------------------------------------------------------------
  !
  !> @param[in]  ntimes     along-track dimesnion size
  !> @param[in]  nxtrack    cross-track dimension size
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !-----------------------------------------------------------------------
  subroutine allocate_rad_vars (ntimes, nxtrack, errstat)

    use m_vars, only: rad_Time, rad_Latitude, rad_Longitude, &
         rad_SolarZenithAngle, rad_ViewingZenithAngle, &
         rad_ViewingAzimuthAngle, rad_SolarAzimuthAngle,&
         out_TerrainHeight, rad_SnowIceFraction,&
         rad_440nm, rad_466nm, rad_477nm, &
         out_ProcessingQualityFlags, rad_GroundPixelQualityFlags
! the following arrays are very large, thus move them out
!         rad_PixelQualityFlags , &
!         rad_Radiance, rad_Wavelength

    use m_vars, only: fFillValue, prev_processingflags
         

    implicit none

    !input variables
    integer (kind=4), intent(in) :: ntimes, nxtrack
    !output variables
    integer (kind=4), intent(inout) :: errstat
    real(kind=4) :: fspecial

    if (errstat /= 0) return

    fspecial = fFillValue 

    allocate (rad_Time(ntimes), &
         rad_Latitude(nxtrack, ntimes), &
         rad_Longitude(nxtrack, ntimes), &
         rad_SolarZenithAngle(nxtrack, ntimes), &
         rad_ViewingZenithAngle(nxtrack, ntimes), &
         rad_SolarAzimuthAngle(nxtrack, ntimes), &
         rad_ViewingAzimuthAngle(nxtrack, ntimes), &
         out_TerrainHeight(nxtrack, ntimes), &
         rad_SnowIceFraction(nxtrack, ntimes), &
         rad_GroundPixelQualityFlags(nxtrack, ntimes), &
         rad_440nm(nxtrack, ntimes), &
         rad_466nm(nxtrack, ntimes), &
         rad_477nm(nxtrack, ntimes), &
! the following large rad arrays are moved to subroutines to save memory
!         rad_Radiance(nwavel, nxtrack, ntimes), &
!         rad_Wavelength(nwavel, nxtrack, ntimes), &
!         rad_PixelQualityFlags(nwavel, nxtrack, ntimes), &
         out_ProcessingQualityFlags(nxtrack, ntimes), &
         prev_processingflags(nxtrack, ntimes), &
         stat = errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "allocated_rad_vars: failed", &
           errstat)
      return
    endif

    ! initialize allocated m_vars variables
    rad_Latitude = fspecial
    rad_Longitude = fspecial
    rad_SolarZenithAngle = fspecial
    rad_ViewingZenithAngle = fspecial
    rad_SolarAzimuthAngle = fspecial
    rad_ViewingAzimuthAngle = fspecial
    rad_SnowIceFraction = fspecial
    rad_GroundPixelQualityFlags = -1
    out_TerrainHeight = fspecial
    rad_440nm = fspecial
    rad_466nm = fspecial
    rad_477nm = fspecial

    ! init out_ProcessingQualityFalgs to Zero
    ! equivalent to all bits cleared
    out_ProcessingQualityFlags = 0

    ! init prev_processingflags 
    prev_processingflags = 0

  end subroutine allocate_rad_vars

  !> Read the dimensions of a L2 CLDO4 file
  !-----------------------------------------------------------------------
  !
  !> @param[in]  swathname  name of group from which to read dimensions
  !> @param[in]  tio_l2obj  libtio file object
  !> @param[in]  nxtrack    cross-track dimension
  !> @param[in]  ntimes     along-track dimension
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !> @author G. Gonzalez Abad May 2021
  !-----------------------------------------------------------------------
  subroutine read_cldo4_dims(tio_l2obj, swathname, nXtrack, nTimes, &
     errstat)

  implicit none
  !input variables
  character (len=*), intent(in) :: swathname
  type (tiof_file_type) :: tio_l2obj
  !output variables
  integer (kind=4), intent(out) :: nxtrack, ntimes
  integer (kind=4), intent(inout) :: errstat

  if (errstat /= 0) return

  call tiof_push_group (tio_l2obj, swathname, errstat)
  call tiof_inq_dimlen (tio_l2obj, "xtrack", nXtrack, errstat)
  call tiof_inq_dimlen (tio_l2obj, "mirror_step", nTimes, errstat)
  call tiof_pop_group (tio_l2obj, errstat)

  if (errstat /= 0) then
    call tell_error (tell_io_read_error, "read_cldo4_dims: failed", errstat)
    return
  endif

end subroutine read_cldo4_dims

  !> Allocate variables associated with CLDO4
  !-----------------------------------------------------------------------
  !
  !> @param[in]  ntimes     along-track dimesnion size
  !> @param[in]  nxtrack    cross-track dimension size
  !> @param      errstat    error handling integer, non-zero = problem
  !
  !> @author G. Gonzalez Abad May 2021
  !-----------------------------------------------------------------------
  subroutine allocate_cldo4_vars (ntimes, nxtrack, errstat)

     use m_vars, only: nasa_SlantColumnAmountO2O2, l2_TerrainPressure,&
           scd_mdqfl, fFillValue,&
           nasa_scdrms, nasa_scduncertainty,&
           rad_RelativeAzimuthAngle,out_RelativeAzimuthAngle

     implicit none

     !input variables
     integer (kind=4), intent(in) :: ntimes, nxtrack
     !output variables
     integer (kind=4), intent(inout) :: errstat

     if (errstat /= 0) return

     allocate (nasa_SlantColumnAmountO2O2(nxtrack, ntimes), &
          stat = errstat)
     nasa_SlantColumnAmountO2O2 = fFillValue

     allocate (l2_TerrainPressure(nxtrack, ntimes), &
          stat = errstat)
     l2_TerrainPressure = fFillValue
 
     allocate (rad_RelativeAzimuthAngle(nxtrack, ntimes),stat=errstat)
     rad_RelativeAzimuthAngle = fFillValue

     allocate (out_RelativeAzimuthAngle(nxtrack, ntimes),stat=errstat)
     out_RelativeAzimuthAngle = fFillValue

     allocate (scd_mdqfl(nxtrack, ntimes), stat=errstat)

     if (errstat /= 0) then
       call tell_error (tell_malloc_error, "allocated_cldo4_vars: failed at scd_msqfl", &
            errstat)
       return
     endif

     !if (run_mode .NE. 'production') then
         allocate (nasa_scdrms(nxtrack, ntimes),stat=errstat)
         nasa_scdrms = fFillValue
         allocate (nasa_scduncertainty(nxtrack, ntimes),stat=errstat)
         nasa_scduncertainty = fFillValue
     !endif

   end subroutine allocate_cldo4_vars

!--------------------------------------------------------------------
   subroutine allocate_diaglog_vars(ntimes,nxtrack,errstat)

   use m_vars, only: rad_waveshift, irr_waveshift

   implicit none
   integer, intent(in):: nxtrack, ntimes
   integer (kind=4), intent(inout) :: errstat

   allocate(rad_waveshift(nxtrack,ntimes),stat=errstat)
   allocate(irr_waveshift(nxtrack), stat=errstat)
 
   ! initialize all shifts to zero
   rad_waveshift = 0.
   irr_waveshift = 0.

   end subroutine allocate_diaglog_vars
!-------------------------------------------------------------------
   subroutine read_diaglog_vars(diaglogfnm,errstat)

   use m_vars, only: rad_waveshift, irr_waveshift
   use m_vars, only: maxradshift,  maxirrshift
   use m_vars, only: max_ns_pixels, max_ew_pixels
   use m_vars, only: option_apply_solshift, option_apply_radshift
   use m_vars, only: lun_debug_shift, ixdebug, itdebug, run_mode
   use netcdf, only: nf90_nowrite

   implicit none
   character(len=*),intent(in):: diaglogfnm
   integer(kind=4),intent(inout):: errstat

   !local variables
   integer(kind=4) :: ix, nx, nt
   integer(kind=4) :: errhere, err1
   integer(kind=2),dimension(:,:), allocatable:: local_mdqfl 

   type (tiof_file_type) :: diaglog_obj

   errhere = 0
   err1 = 0 ! for loca_mdqfl
   
   ! initial nx,nt large enough, they will be replaced if successful
   nx = max_ns_pixels
   nt = max_ew_pixels 

   call tiof_open(trim(diaglogfnm),diaglog_obj,nf90_nowrite,errhere)
   if (errhere /= 0) then
      write(*,*) 'read_diaglog_vars: failed to open diaglog file.'
      write(*,*) 'continue with assumption of shifts = 0.'
      write(*,*) 'force option_apply_solshift=0, option_apply_radshift=0'
      ! force option_apply_XXXshift to 0
      option_apply_solshift = 0
      option_apply_radshift = 0
      ! allocate & initialize to 0 wavelength shift
      write(*,*) '      diaglog: fake nx,nt=',nx,nt
      call allocate_diaglog_vars(nt,nx,errstat)
      ! if allocation fails, errstat tells the calling program 
      return
   else 
      call tiof_inq_dimlen(diaglog_obj, "xtrack", nx, errhere)
      call tiof_inq_dimlen(diaglog_obj, "mirror_step", nt, errhere)
      write(*,*) '      diaglog: nx,nt=',nx,nt
      ! allocate & initialize to 0 wavelength shift 
      call allocate_diaglog_vars(nt,nx,errstat)

      allocate(local_mdqfl(nx,nt),stat=errstat)
      local_mdqfl = 0

      ! read wavelength shift 
      call tiof_get2d_r4(diaglog_obj,"radshi", [0,0], [nt,nx], &
        rad_waveshift, errhere)

      call tiof_get1d_r4(diaglog_obj,"solshi", [0], [nx], &
        irr_waveshift, errhere)

      call tiof_get2d_i2(diaglog_obj,"mdqfl",[0,0], [nt,nx], &
        local_mdqfl, err1)

      call tiof_close(diaglog_obj, errhere)

      if (errhere /=0) then
          rad_waveshift = 0.
          irr_waveshift = 0.
          write(*,*) "warning: read_diag_nvars failed to read, assume 0 shifts"
          return
      else 
          if (err1 /=0) then !diaglog does not contain mdqfl
             write(*,*) "warning: rad_diag_nvars assume local_mdqfl=0 for radshi"
          else !diaglog contains mdqfl from fitting
             ! do not use radshi for failed or suspect fitting
             write(*,*) "read_diag_nvars: set rad_waveshift=0 for bad/suspect fit"
             where (local_mdqfl /= 0) ! bad or suspect fit
                 rad_waveshift = 0.
             endwhere
             ! local_mdqfl no longer needed
             deallocate(local_mdqfl)
          endif ! err1
      endif !errhere
   endif

   ! set large shift to zero
   where ((rad_waveshift > maxradshift).or.(rad_waveshift < -maxradshift)) 
       rad_waveshift = 0.
   endwhere
   where((irr_waveshift > maxirrshift).or.(irr_waveshift < -maxirrshift))
       irr_waveshift = 0.
   endwhere

   ! option for wavelength shift
   write(*,*) 'option_apply_solshift=',option_apply_solshift
   if (option_apply_solshift == 0) then
       irr_waveshift = 0.
   endif
   write(*,*)'option_apply_radshift=',option_apply_radshift
   if (option_apply_radshift == 0) then
       rad_waveshift = 0.
   endif

   ! debug
   if ((trim(run_mode) .eq. 'development').and.(ixdebug .ge. 0)) then
      write(*,*) 'writing debug_shift.txt'
      open(unit=lun_debug_shift,file='debug_shift.txt')
      write(lun_debug_shift,*) 'option_apply_solshift=',option_apply_solshift
      write(lun_debug_shift,*) 'option_apply_radshift=',option_apply_radshift
      write(lun_debug_shift,*) 'ixdebug,itdebug=',ixdebug,itdebug
      write(lun_debug_shift,*)'ix    irr_waveshift         rad_waveshift'
      do ix = 1, nx
         write(lun_debug_shift,*)ix,irr_waveshift(ix),rad_waveshift(ix,itdebug)
      enddo
      close(lun_debug_shift)
   endif

   end subroutine read_diaglog_vars

!--------------------------------------------------------------------
   subroutine allocate_init_desvars(nxtrack, errstat)

   use m_vars, only: scddes_hour
   use m_vars, only: lun_desfac_fnm
   use m_vars, only: name_desfac_fnm
   use m_vars, only: name_desfac_dir,option_destripe_scd

   implicit none
   integer, intent(in):: nxtrack
   integer, intent(inout):: errstat
   character(len=255):: filename_deshour
   real :: tmpvalue
   integer :: ix

   filename_deshour = trim(name_desfac_dir)//trim(name_desfac_fnm)

     ! allocate array
     allocate (scddes_hour(nxtrack), stat=errstat)
     if (errstat /= 0) return
     ! initialize scddes
     scddes_hour = 1.

     if (option_destripe_scd .eq. 1) then 
     ! temporary fix below: read scddes from filename_deshour
        write(*,*) '   reading destripe factor '//trim(filename_deshour)
        open(unit=lun_desfac_fnm,file=trim(filename_deshour))
        do ix = 1, nxtrack
           read(lun_desfac_fnm,*) tmpvalue
           scddes_hour(ix) = tmpvalue
        enddo
        close(lun_desfac_fnm)
     else
        write(*,*) '  desfac=1.0 for all pixels'
     endif 
        
   end subroutine allocate_init_desvars
!-----------------------------------------------------------

   subroutine destripe_o2o2scd(arrin,nXtrack,nTimes,arrout)

   use m_vars, only: scddes_hour, option_destripe_scd
   use m_vars, only: fFillValue

   implicit none
   
   real(kind=8),dimension(nXtrack,nTimes),intent(in)::arrin
   integer, intent(in):: nXtrack,nTimes
   real(kind=4),dimension(nXtrack,nTimes),intent(out)::arrout
   
   real:: fspecial = fFillValue 
   real:: desfac, thisarr
   ! local variable
   real, dimension(nXtrack):: scddes

   integer:: ix, it

   arrout = fspecial

   scddes = scddes_hour

   if (option_destripe_scd .eq. 1) then 
   write(*,*) '    destriping scd'
   do ix = 1, nXtrack
      desfac = scddes(ix)
      if (desfac .gt. 0.) then
        do it = 1, nTimes
          thisarr = real(arrin(ix,it),kind=4)
          if (thisarr < 0.) then
             arrout(ix,it) = fFillValue
          else
             arrout(ix,it) = thisarr/desfac
          endif
        enddo
      endif
   enddo 
   else
       write(*,*) '     scd is not destriped.'
       arrout = real(arrin,kind=4)
   endif

   end subroutine destripe_o2o2scd

  !-----------------------------------------------------------
end module m_read_input_tio
