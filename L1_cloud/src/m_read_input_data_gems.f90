!> Read radiance and geolocation data from L1C GEMS files
module m_read_input_data_gems
  use tio_module
  use tell_module
  use netcdf, only: nf90_nowrite, nf90_noerr, nf90_global, nf90_get_att

  public read_input_data_gems, read_solar_data_gems, bad_rad_lambda_gems
  private read_cld_dimensions_gems, read_cld_geo_data_gems, &
       read_cld_rad_data_gems, read_date_gems, read_earth_sun_dist_gems, &
       read_sol_data_gems, bad_irrad_lambda_gems

contains

  !> Top-level subroutine to read in an L1C GEMS radiance file
  !---------------------------------------------------------------------
  !
  !> @param[in] l1bfile  filename for L1B netCDF radiance file
  !> @param     errstat  error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2020
  !---------------------------------------------------------------------
  subroutine read_input_data_gems(l1bfile, errstat)

    use m_vars, only: input_data_path, iLine, wrt_solar, wmin, wmax, &
         wmin2, wmax2, set_wmin, set_wmax, wave_long, wave_short, nLines, &
         start_line, max_lines, nXtrack, nTimes, nWavel, &
         n_input, n_missing, nwl, qc, nwave, min_wl, nWavel, &
         ws, fs, nsolwave, Year, Month, Day, nTimes, nXtrack
    use m_read_input_data_tio, only: alloc_scan
    use m_read_solar_data_tio, only: write_solar_tio
    use m_strpos

    implicit none

    !input variables
    character (len=*), intent (in) :: l1bfile
    !output variables
    integer (kind=4), intent (inout) :: errstat

    !local variables
    character (len = 200) :: filenamepath
    character (len = 128) :: logmsg
    integer (kind = 4) :: i

    if (errstat /= 0) return

    !set file path+name, swath name
    filenamepath=trim(input_data_path)//l1bfile

    if (iLine == 0) then
      !-----------------------------------------------------------------
      !If on first line, perform some setup

      write(logmsg,"(A30, 3A)")'read_input_data_gems: filename ', &
           trim(filenamepath)
      call tell_log(1,logmsg)

      !set wavelength bounds
      if (wrt_solar) then
        wmin=355.0d0
        wmax=500.0d0
        wmin2=310.0d0
        wmax2=375.0d0
      endif
      wmin2 = 330.0d0
      wmax2 = 367.0d0
      if (.not. set_wmin) wmin = 345.5d0
      if (.not. set_wmax) wmax = 354.5d0
      wave_long=362.5d0
      wave_short=345.4d0

      !open file, read variable dimension sizes
      call read_cld_dimensions_gems(l1bfile, nTimes, nXtrack, nWavel, errstat)
      if(errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "read_cld_dimensions_gems: failed", &
             errstat)
        return
      endif

      iLine=start_line
      if (max_lines > 0) then
        write(logmsg,"(A40, I6)") 'read_input_data_gems: changing nTimes to ',&
             max_lines
        call tell_log(1,logmsg)
        nTimes=max_lines+start_line
      endif
      nLines=nTimes

      !Allocate arrays for variables to be read in
      call tell_log(1,'read_input_data_gems: calling alloc_scan')
      call alloc_scan(errstat)
      if(errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "alloc_scan: failed", &
             errstat)
        return
      endif

      call read_cld_geo_data_gems(l1bfile, errstat)
      if(errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "read_cld_geo_data_gems: failed", &
             errstat)
        return
      endif
      call tell_log(2,'read_cld_geo_data_gems: success')

      !Need the month in order to read in correct calibration climatologies
      call read_date_gems (l1bfile, Year, Month, Day, errstat)
      if (errstat /= 0) then
        call tell_log (1, "read_input_data_gems: failed to read date")
        Month=1
      else
        write (logmsg, *) "Date is: ", Year, Month, Day
        call tell_log (1, logmsg)
      endif

      write (logmsg,"(A10,2F7.2)") 'wmin2 wmax2 ',wmin2,wmax2
      call tell_log(3,logmsg)


    endif !iLine == 0


    !----------------------------------------------------------------
    ! For each line of the input file

    !Keep track of number of value input
    n_input = n_input + nXtrack

    call read_cld_rad_data_gems(l1bfile, errstat)
    if(errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_rad_data_gems: failed", &
           errstat)
      return
    endif
    call tell_log(2,'read_cld_rad_data_gems: success')

    !Print check of wavelengths
    if (iLine == start_line) then
      call tell_log(1,'nwl, iLine, wmin, wmax')
      write(logmsg,"(2I6,2F7.3)") nwl, iLine, wmin, wmax
      call tell_log(1,logmsg)
    endif

    nwave=nwl

    if (iLine /= start_line) then
      !check for missing data
      if (nwl > nWavel .or. nwl < min_wl) then
        qc(:,iLine) = ibset(qc(:,iLine),14)
        n_missing = n_missing + nXtrack
        write (logmsg,"(A13, 4I6)") 'missing line ',iLine, nwl, nWavel, min_wl
        call tell_log(1,logmsg)
      endif ! missing wavelength data
    endif ! start_line

    ! read solar flux
    !===================
    if (iLine == start_line) then
      call read_solar_data_gems(errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "read_solar_data_gems: failed", &
             errstat)
        return
      endif
      call tell_log(3,'irradiance')
      do i=0,nsolwave-1
        write(logmsg,'(i4,2e12.4)') i,ws(i,0),fs(i,0)
        call tell_log(3,logmsg)
      enddo
    endif ! iLine==start_line

    ! option to write out solar data and quit
    if (iLine == start_line) then
      if (wrt_solar) then
        call write_solar_tio(errstat)
        if (errstat /= 0) then
          call tell_error (tell_io_read_error, &
               "write_solar_tio: failed", &
               errstat)
          return
        endif
      endif
    endif

  end subroutine read_input_data_gems


  !>Open GEMS L1C radiance file and get dimensions
  !---------------------------------------------------------------------
  !
  !> @param[in]  l1bfile    filename for L1C netCDF radiance file
  !> @param[out] nTimes     Along-track dimension size
  !> @param[out] nXtrack    Across-track dimension size
  !> @param[out] nWavel     Spectral dimension size
  !> @param      errstat    error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2020
  !---------------------------------------------------------------------
  subroutine read_cld_dimensions_gems(l1bfile, nTimes, nXtrack, nWavel, &
       errstat)

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile
    !output variables
    integer (kind=4), intent(out) :: nTimes, nXtrack, nWavel
    integer (kind=4), intent (inout) :: errstat

    type (tiof_file_type) :: tio_l1obj

    !local variables
    character (len=128) :: logmsg

    if (errstat /= 0) return
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_time_set_taix_epoch ("2000-01-01T00:00:00Z", errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_y", nXtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_x", nTimes, errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_band", nWavel, errstat)

    call tell_log(1,'read_input_data_gems: nTimes, nXtrack, nWavel')
    write(logmsg,"(3I6)")nTimes,nXtrack,nWavel
    call tell_log(1,logmsg)

    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_cld_dimensions_gems: failed to open L1B file", &
           errstat)
      return
    endif

  end subroutine read_cld_dimensions_gems


  !>Read all geolocation data from GEMS L1C radiance file in one pass
  !---------------------------------------------------------------------
  !
  !> @param[in] l1bfile    filename for L1C netCDF radiance file
  !> @param     errstat    error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2020
  !---------------------------------------------------------------------
  subroutine read_cld_geo_data_gems(l1bfile, errstat)

    use m_vars, only: time, lat, lon, sza, sazimuth, sat_zen, &
         vazimuth, terr_height, geoflg, anomflg, mflg, nLines, nXtrack, &
         azimuth, fill_value, gems_snow_index

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=4), dimension(nLines,nXtrack) :: tmp_lat, tmp_lon, tmp_sza, &
         tmp_saa, tmp_vza, tmp_vaa
    integer (kind=2), dimension(nLines,nXtrack) :: tmp_hgt, tmp_gflg, tmp_snow
    integer (kind=1), dimension(nLines,nXtrack) :: tmp_xtqf
    integer, dimension(1:2), parameter :: flip = (/2,1/)
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, "image_acquisition_time", [0], [nLines], &
         time, errstat)
    call tiof_get2d_r4 (tio_l1obj, "pixel_latitude", [0,0], &
         [nXtrack,nLines], tmp_lat, errstat)
    call tiof_get2d_r4 (tio_l1obj, "pixel_longitude", [0,0], &
         [nXtrack,nLines], tmp_lon, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sun_zenith_angle", [0,0], &
         [nXtrack,nLines], tmp_sza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sun_azimuth_angle", [0,0], &
         [nXtrack,nLines], tmp_saa, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sc_zenith_angle", [0,0], &
         [nXtrack,nLines], tmp_vza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sc_azimuth_angle", [0,0], &
         [nXtrack,nLines], tmp_vaa, errstat)
    call tiof_get2d_i2 (tio_l1obj, "terrain_height", [0,0], &
         [nXtrack,nLines], tmp_hgt, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ground_pixel_quality_flag", [0,0], &
         [nXtrack,nLines], tmp_gflg, errstat)
    call tiof_get2d_i2 (tio_l1obj, "snow_index", [0,0], &
         [nXtrack,nLines], tmp_snow, errstat)
    call tiof_get2d_i1 (tio_l1obj, "xtrack_quality_flag", [0,0], &
         [nXtrack,nLines], tmp_xtqf, errstat)
    !No main data quality flag
    mflg=0
    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_geo_data_gems: failed to read geolocation data", &
           errstat)
      return
    endif

    ! flip 2D array axes
    lat = reshape(tmp_lat,(/nXtrack,nLines/),order=flip)
    lon = reshape(tmp_lon,(/nXtrack,nLines/),order=flip)
    sza = reshape(tmp_sza,(/nXtrack,nLines/),order=flip)
    sazimuth = reshape(tmp_saa,(/nXtrack,nLines/),order=flip)
    sat_zen = reshape(tmp_vza,(/nXtrack,nLines/),order=flip)
    sazimuth = reshape(tmp_vaa,(/nXtrack,nLines/),order=flip)
    terr_height = reshape(tmp_hgt,(/nXtrack,nLines/),order=flip)
    ! FIXME - flag translation!!!!
    anomflg = reshape(tmp_xtqf,(/nXtrack,nLines/),order=flip)
    gems_snow_index = reshape(tmp_snow,(/nXtrack,nLines/),order=flip)
    geoflg = int(reshape(tmp_gflg,(/nXtrack,nLines/),order=flip), kind=4)


    !Calculate relative azimuth angle
    azimuth(:,:)=sazimuth(:,:)+180.0-vazimuth(:,:)
    where(azimuth(:,:) < -180.) azimuth(:,:)=azimuth(:,:)+360.
    where(azimuth(:,:) > 180.) azimuth(:,:)=azimuth(:,:)-360.
    azimuth(:,:)=abs(azimuth(:,:))
    where(azimuth(:,:) > 360.0) azimuth(:,:)=fill_value

  end subroutine read_cld_geo_data_gems


  !>Read one line of radiance data from GEMS L1C file
  !---------------------------------------------------------------------
  !
  !> @param[in]  l1bfile    filename for L1C netCDF radiance file
  !> @param      errstat    error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2020
  !---------------------------------------------------------------------
  subroutine read_cld_rad_data_gems(l1bfile, errstat)

    use m_vars, only: iLine, wmin2, wmax2, nwl, f12d, w12d, quality_flagL, &
         nWavel, nXtrack
    use m_read_solar_data_tio, only: calc_wl_line

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=4), dimension(nWavel, nXtrack, 1) :: tio_rad
    real (kind=4), dimension(1, nXtrack, nWavel) :: in_rad
    real (kind=4), dimension(nWavel, nXtrack) :: wl_local
    real (kind=4), dimension(nXtrack, nWavel) :: in_wvl
    integer (kind=2), dimension(nWavel, nXtrack, 1) :: tio_flg
    integer (kind=2), dimension(1, nXtrack, nWavel) :: in_flg
    integer (kind=4) :: ih, il
    integer, dimension(1:2), parameter :: flip = (/2,1/)
    integer, dimension(1:3), parameter :: flip3d = (/3,2,1/)

    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    !open file, read wavelength, radiance, flag values
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get3d_r4 (tio_l1obj, "image_pixel_values", [0,0,iLine-1], &
         [nWavel,nXtrack,1], in_rad, errstat)
    call tiof_get2d_r4 (tio_l1obj, "wavelength", [0,0], &
         [nWavel,nXtrack], in_wvl, errstat)
    call tiof_get3d_i2 (tio_l1obj,"bad_pixel_mask", [0,0,iLine-1], &
         [nWavel,nXtrack,1], in_flg, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_rad_data_gems: failed to read radiance data", &
           errstat)
      return
    endif

    ! Flip arrays
    tio_rad = reshape(in_rad,(/nWavel,nXtrack,1/),order=flip3d)
    tio_flg = reshape(in_flg,(/nWavel,nXtrack,1/),order=flip3d)
    wl_local = reshape(in_wvl,(/nWavel,nXtrack/),order=flip)

    !Determine limits of wavelength window
    errstat=calc_wl_line(nXtrack, nWavel, wmin2, wmax2, wl_local, &
         il, ih, nwl)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_rad_data_gems: calc_wl_line: failed", &
           errstat)
      return
    endif

    !Copy radiance, wavelength and quality in window to main arrays
    f12d(0:nwl-1,0:nXtrack-1)=tio_rad(il:ih,1:nXtrack,1)
    w12d(0:nwl-1,0:nXtrack-1)=wl_local(il:ih,1:nXtrack)
    quality_flagL(1:nwl,1:nXtrack)=tio_flg(il:ih,1:nXtrack,1)

    w12d(nwl:nWavel-1,:)=0.d0
    f12d(nwl:nWavel-1,:)=0.d0

  end subroutine read_cld_rad_data_gems


  !> Read numerical month for use by calibration files
  !-----------------------------------------------------------------------
  !
  !> @param[in]  l1bfile    L1C filename
  !> @param[out] year       integer year
  !> @param[out] month      integer month
  !> @param[out] day        integer day
  !> @param[out] errstat    error handling integer, non-zero = failure
  !
  !> @author E. O'Sullivan October 2020
  !-----------------------------------------------------------------------
  subroutine read_date_gems (l1bfile, year, month, day, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile

    !output variables
    integer (kind=4), intent(out) :: year, month, day
    integer (kind=4), intent(inout) :: errstat

    !local variables
    character (len=15):: datestamp
    integer :: ncerr
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    ncerr = nf90_get_att (tio_l1obj%fileid, nf90_global, &
         "mission_reference_time", datestamp)
    call tiof_close (tio_l1obj, errstat)
    if (ncerr /= nf90_noerr .or. errstat /=0) then
      errstat = ncerr
      call tell_error (tell_io_read_error, &
           "read_date_gems: failed to read date", errstat)
      return
    endif

    read(datestamp, '(i4,i2,i2,7x)') year, month, day

  end subroutine read_date_gems


  !> Read Earth-Sun distance from GEMS radiance file
  !-----------------------------------------------------------------------
  !
  !> @param[in]  l1bfile    L1C filename
  !> @param[out] dist       Earth-Sun distance in m
  !> @param[out] errstat    error handling integer, non-zero = failure
  !
  ! NB: If file attribute is empty, use default earth-sun dist
  !
  !> @author E. O'Sullivan October 2020
  !-----------------------------------------------------------------------
  subroutine read_earth_sun_dist_gems (l1bfile, dist, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile

    !output variables
    real (kind=4), intent(out) :: dist
    integer (kind=4), intent(inout) :: errstat

    !local variables
    integer :: ncerr
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    ncerr = nf90_get_att (tio_l1obj%fileid, nf90_global, &
         "earth_sun_distance", dist)
    call tiof_close (tio_l1obj, errstat)
    if (ncerr /= nf90_noerr .or. errstat /=0) then
      errstat = ncerr
      call tell_error (tell_io_read_error, &
           "read_earth_sun_dist_gems: failed to read dist", errstat)
      return
    endif

    ! Check whether value is sensible
    if (dist < 100.0) then
      dist = 149597870000.0  ! default
      call tell_log (1, "WARNING: Default earth-sun distance used")
    endif

  end subroutine read_earth_sun_dist_gems


  !>Top-level subroutine to read in a GEMS irradiance file
  !---------------------------------------------------------------------
  !
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2020
  !---------------------------------------------------------------------
  subroutine read_solar_data_gems(errstat)

    use m_vars, only: fs, nsolwave, dist_rad, dist_irrad, &
         irrad_filename_nc, filename_in_nc
    use m_lambda_qual
    use m_LUN_set
    use m_pgs_include

    implicit none

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !internal variables
    integer (kind = 4) :: nTimes, nXtrack, nWavel
    character (len = 200) ::  logmsg


    if (errstat /= 0) return

    !open IRR1B file
    write(logmsg,*) 'read_solar_data_gems: opening ',trim(irrad_filename_nc)
    call tell_log(2,logmsg)

    call read_cld_dimensions_gems(filename_in_nc, nTimes, nXtrack, nWavel, &
         errstat)
    if(errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read solar dimensions: failed", &
           errstat)
      return
    endif

    !Read in solar data
    call read_sol_data_gems(irrad_filename_nc, nXtrack, nWavel, errstat)
    if(errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_sol_data_gems: failed", &
           errstat)
      return
    endif

    !set processing quality_flags
    call bad_irrad_lambda_gems(nXtrack, errstat)

    ! FIXME - AS YET GEMS IRRAD DOES NOT CONTAIN E-S DISTANCE.
    !correction for earth-sun distance
    !call read_earth_sun_dist_gems(irrad_filename_nc,dist_irrad,errstat)
    !if (errstat /= 0) then
    !  call tell_error (tell_io_error, &
    !       "read_earth_sun_dist_gems: failed for irradiance file", &
    !       errstat)
    !  return
    !endif
    dist_irrad=149597870000.0

    call read_earth_sun_dist_gems(filename_in_nc,dist_rad,errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_earth_sun_dist_gems: failed for radiance file", &
           errstat)
      return
    endif

    fs(0:nsolwave-1,:)=fs(0:nsolwave-1,:)*(dist_irrad/dist_rad)**2


  end subroutine read_solar_data_gems


  !>Read irradiance data from GEMS netCDF file
  !---------------------------------------------------------------------
  !
  !> @param[in]    filename_sol_nc filename for L1B netCDF irradiance file
  !> @param[in]    nXtrack         size of dimension across direction of scan
  !> @param[in]    nWavel          size of spectral dimension
  !> @param[inout] errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2020
  !---------------------------------------------------------------------
  subroutine read_sol_data_gems(filename_sol_nc, nXtrack, nWavel, errstat)

    use m_vars, only: wmin2, wmax2, ws, fs, nsolwave, ierr, irr_quality_flagL
    use m_read_solar_data_tio, only: calc_wl_line

    implicit none
    !input variables
    character (len=*), intent (in) :: filename_sol_nc
    integer (kind=4), intent (in) :: nXtrack, nWavel
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=4), dimension(nWavel, nXtrack) :: tio_wvl, tio_rad
    integer (kind=2), dimension(nWavel,nXtrack) :: tio_flg
    real (kind=4), dimension(nXtrack, nWavel) :: tmp_wvl, tmp_irr
    integer (kind=2), dimension(nXtrack, nWavel) :: tmp_pqf
    integer, dimension(1:2), parameter :: flip=(/2, 1/)
    integer (kind=4) :: ih, il
    character (len=128) :: logmsg

    type (tiof_file_type) :: tio_irrobj

    if (errstat /= 0) return

    !open file, read wavelength, radiance, flag values
    call tiof_open (filename_sol_nc, tio_irrobj, nf90_nowrite, errstat)
    call tiof_get2d_r4 (tio_irrobj, "wavelength", [0,0], &
         [nWavel, nXtrack], tmp_wvl, errstat)
    call tiof_get2d_r4 (tio_irrobj, "image_pixel_values", [0,0], &
         [nWavel, nXtrack], tmp_irr, errstat)
    call tiof_get2d_i2 (tio_irrobj, "bad_pixel_mask", [0,0], &
         [nWavel, nXtrack], tmp_pqf, errstat)
    call tiof_close (tio_irrobj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_sol_data_gems: failed to read irradiance data", &
           errstat)
      return
    endif

    !Reshape arrays
    tio_wvl = reshape (tmp_wvl,(/nwavel,nxtrack/),order=flip)
    tio_rad = reshape (tmp_irr,(/nwavel,nxtrack/),order=flip)
    tio_flg = reshape (tmp_pqf,(/nwavel,nxtrack/),order=flip)


    !Determine limits of wavelength window
    errstat=calc_wl_line(nXtrack, nWavel, wmin2, wmax2, tio_wvl, &
         il, ih, nsolwave)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_sol_data_gems: calc_wl_line: failed", &
           errstat)
      return
    endif

    write(logmsg,"(A9,I6)") 'nsolwave ',nsolwave
    call tell_log(2,logmsg)

    !allocate solar data arrays
    if (allocated(ws)) deallocate(ws)
    if (allocated(fs)) deallocate(fs)
    if (allocated (irr_quality_flagL)) deallocate( irr_quality_flagL)
    allocate( ws(0:nsolwave-1,0:nXtrack-1), &
         fs(0:nsolwave-1,0:nXtrack-1), &
         irr_quality_flagL(nWavel,nXtrack), stat=ierr )
    if (ierr .ne. 0) then
      call tell_error (tell_malloc_error, &
           "read_sol_data_gems: failed to allocate memory", &
           errstat)
      return
    endif
    irr_quality_flagL(1:nWavel,1:nXtrack) = -1

    !Copy radiance, wavelength and quality in window to main arrays
    fs(0:nsolwave-1,0:nXtrack-1)=tio_rad(il:ih,1:nXtrack)
    ws(0:nsolwave-1,0:nXtrack-1)=tio_wvl(il:ih,1:nXtrack)
    irr_quality_flagL(1:nsolwave,1:nXtrack)=tio_flg(il:ih,1:nXtrack)


  end subroutine read_sol_data_gems


  !>Set irradiance quality flags
  !-------------------------------------------------------------------------
  !
  !> @param[in]  nXtrack  Cross-track dimension size
  !> @param      errstat  error-tracking integer
  !
  ! @author  E. O'Sullivan  October 2020
  !-------------------------------------------------------------------------
  subroutine bad_irrad_lambda_gems(nXtrack, errstat)

    use m_vars, ONLY: irr_quality_flagL, qc, nsolwave, ws, fs, wmin, wmax
    use m_find

    implicit none

    integer, intent(in) :: nXtrack, errstat
    integer :: iw, ip
    integer :: iw_start, iw_end
    logical, dimension(nsolwave,nXtrack) :: pxl_error
    character (len=128) :: logmsg


    if (errstat /= 0) return

    pxl_error = .false.

    do ip=0,nXtrack-1

      !find starting and ending wavelengths
      !====================================
      iw_end=maxval(find2(ws(:,ip) <= wmax,count(ws(:,ip) <= wmax)))-1
      iw_start=minval(find2(ws(:,ip) >= wmin,count(ws(:,ip) >= wmin)))-1

      ! check the image quality flags
      !==============================
      do iw=iw_start,iw_end
        if (irr_quality_flagL(iw+1,ip+1) /= 0) pxl_error(iw+1,ip+1) = .true.
        if(pxl_error(iw+1,ip+1)) then
          qc(ip,:)=IBSET(qc(ip,:),11)
          fs(iw,ip)=0.
          write(logmsg,"(A29,I6,2X,F10.6)") 'bad irradiance scan position ', &
               ip, ws(iw,ip)
          call tell_log(1,logmsg)
        endif
      enddo !iw
    enddo !ip

  end subroutine bad_irrad_lambda_gems


  !>Set radiance quality flags
  !-------------------------------------------------------------------------
  !
  !> @param[in]  ip       Cross-track pixel of interest
  !> @param[in]  iLine    Along track step of interest
  !> @param      errstat  error-tracking integer
  !
  ! @author  E. O'Sullivan  October 2020
  !-------------------------------------------------------------------------
  subroutine bad_rad_lambda_gems(ip, iLine, errstat)

    use m_vars, ONLY: quality_flagL, qc, wmin, wmax
    use m_cloud_pres_mod, ONLY: f1p, w1p
    use m_find

    implicit none

    integer, intent(in) :: ip, iLine, errstat
    integer :: iw
    integer :: iw_start, iw_end
    logical :: pxl_error!, pxl_warning

    if (errstat /= 0) return

    pxl_error=.false.

    !find starting and ending wavelengths
    !====================================
    if (count(w1p(:)<= wmax .and. w1p(:) >= wmin) > 0) then
      iw_end=maxval(find2(w1p(:) <= wmax .and. w1p(:) >= wmin, &
           count(w1p(:) <= wmax .and. w1p(:) >= wmin)))-1
    else
      iw_end=-1
    endif
    if (count(w1p(:) >= wmin) > 1) then
      iw_start=minval(find2(w1p(:) >= wmin,count(w1p(:) >= wmin)))-1
    else
      iw_start=-1
    endif

    ! check the image quality flags
    !================================
    if (iw_start >= 0 .and. iw_end > 0) then
      do iw=iw_start, iw_end
        if (quality_flagL(iw+1,ip+1) /= 0 .or. f1p(iw+1) .le. 0.0) then
          pxl_error = .true.
        endif
        if(pxl_error) then
          f1p(iw+1)=0.
          qc(ip,iLine)=IBSET(qc(ip,iLine),9)
        endif ! if pxl_error
      enddo ! loop over
    else
      qc(ip,iLine)=IBSET(qc(ip,iLine),9)
    endif

  end subroutine bad_rad_lambda_gems


end module m_read_input_data_gems
