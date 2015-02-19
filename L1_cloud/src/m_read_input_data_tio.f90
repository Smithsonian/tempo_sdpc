!read input data from L1B netCDF file
module m_read_input_data_tio
  use cld_names_module
  use tio_module
  use tell_module
  use netcdf, only : nf90_nowrite

  public read_input_data_tio
  private read_cld_dimensions, read_cld_geo_data, read_cld_rad_data, &
       alloc_scan

contains

  subroutine read_input_data_tio(l1bfile, errstat)
    !read in a netCDF radiance file
    use m_vars, only: iprt, input_data_path, iLine, wrt_solar, wmin, wmax, &
         wmin2, wmax2, set_wmin, set_wmax, wave_long, wave_short, nLines, &
         start_line, max_lines, nXtrack, nTimes, nWavel, meas_qual_flg, &
         mflg, n_input, n_missing, nwl, w12d, qc, nwave, min_wl, &
         ws, fs, nsolwave, read_he4
    use m_strpos
    use m_read_solar_data_tio

    implicit none

    !input variables
    character (len=*), intent (in) :: l1bfile 

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !local variables
    character (len=200) :: filenamepath, swathname
    integer (kind = 4) :: PGS_TD_TAItoUTC, i
    logical :: uvswath

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    !set file path+name, swath name
    filenamepath=trim(input_data_path)//l1bfile
    uvswath = strpos (l1bfile, 'BRUG') > 0
    if (uvswath) then
      swathname = "band_540_740_nm"
    else
      call tell_error (tell_invalid_parm, &
           "read_input_data_tio: input file is not OMI L1 UV swath", &
           errstat)
      return
    endif
    if (iLine == 0) then
      !-----------------------------------------------------------------
      !If on first line, perform some setup

      if (iprt > 0) print *,'read_input_data_tio: filename ', &
           trim(filenamepath), '   ', trim(swathname)
      !MAY ALREADY DONE IN he4 READ
      if (.not. read_he4) then
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
      endif

      !open file, read variable dimension sizes
      call read_cld_dimensions(l1bfile, tio_l1obj, swathname, errstat)
      if(errstat < 0) then 
        call tell_error (tell_io_read_error, &
             "read_cld_dimensions: failed", &
             errstat)
        return
      endif

      iLine=start_line
      if (max_lines > 0 .and. iprt > 0) then
        print *,'read_input_data_tio: changing nTimes to ',max_lines
        nTimes=max_lines+start_line
      endif
      nLines=nTimes

      !MAY ALREADY BEDONE DURING he4 READ
      if (.not. read_he4) then
        !Allocate arrays for variables to be read in
        if (iprt >= 1) print *,'read_input_data_tio: calling alloc_scan'
        call alloc_scan(errstat)
        if(errstat < 0) then 
          call tell_error (tell_malloc_error, &
               "alloc_scan: failed", &
               errstat)
          return
        endif
      endif

    !Read in all geolocation data at once
!    call read_cld_geo_data2(l1bfile, tio_l1obj, swathname, errstat)
!    if(errstat < 0) then 
!      call tell_error (tell_io_read_error, &
!           "read_cld_geo_data2: failed", &
!           errstat)
!      return
!    endif
!    if(iprt >=2) print *,'read_cld_geo_data2: success'

    endif !iLine == 0


    !----------------------------------------------------------------
    ! For each line of the input file

    !Read in geolocation data line by line
    call read_cld_geo_data(l1bfile, tio_l1obj, swathname, errstat)
    if(errstat < 0) then 
      call tell_error (tell_io_read_error, &
           "read_cld_geo_data: failed", &
           errstat)
      return
    endif
    if(iprt >=2) print *,'read_cld_geo_data: success'

    !Keep track of number of value input
    n_input = n_input + nXtrack



    !Check data quality flags, if data missing set radiance error flag
    if(btest(mflg(iLine),0) .or. btest(mflg(iLine),1) &
         .or. btest(mflg(iLine),3) .or. btest(mflg(iLine),12)) then
      meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),1)
      n_missing = n_missing + nXtrack 
      if (iprt >= 1) then
        print *,'missing line ',iLine, btest(mflg(iLine),0), &
             btest(mflg(iLine),1), btest(mflg(iLine),3), btest(mflg(iLine),12)
      endif
    endif
    !Set internal measurement quality flags
    if(btest(mflg(iLine),2) .or. btest(mflg(iLine),4) &
         .or. btest(mflg(iLine),5) .or. btest(mflg(iLine),6) &
         .or. btest(mflg(iLine),8) .or. btest(mflg(iLine),9) &
         .or. btest(mflg(iLine),11)) &
         meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),2)

    if(btest(mflg(iLine),7)) meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),3)
    if(btest(mflg(iLine),10)) meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),4)

    !Read in radiance data 
    call read_cld_rad_data(l1bfile, tio_l1obj, swathname, errstat)
    if(errstat < 0) then 
      call tell_error (tell_io_read_error, &
           "read_cld_rad_data: failed", &
           errstat)
      return
    endif
    if(iprt > 1) print *,'read_cld_rad_data: success'

    !Print check of wavelengths
    if (iprt >= 1 .and. iLine == start_line) then
      print *, 'nwl, iLine, wmin, wmax'
      print *, nwl, iLine, wmin, wmax
    endif

    nwave=nwl

    if (iLine == start_line) then
      if (iprt >= 3) then
        write(6,"(6f12.2)") w12d(0:nwl-1,0)
      endif
    else 
      !check for missing data
      if (nwl > nWavel .or. nwl < min_wl) then
        qc(:,iLine) = IBSET(qc(:,iLine),14)
        n_missing = n_missing + nXtrack 
        errstat=-1
        if (iprt >= 1) print *,'missing line ',iLine, nwl, nWavel, min_wl
        if (iprt >= 3) then
          write(6,"(6f12.2)") w12d(0:nwl-1,0)
        endif
      endif ! missing wavelength data
    endif ! start_line

    ! read solar flux
    !===================
    if (iLine == start_line) then
      call read_solar_data_tio(errstat)
      if (errstat < 0) then
        call tell_error (tell_io_read_error, &
             "read_solar_data_tio: failed", &
             errstat)
        return
      endif
      if (iprt > 1) then
        print *,'irradiance'
        do i=0,nsolwave-1
          write(*,'(i4,2e12.4)') i,ws(i,0),fs(i,0)
        enddo
      endif ! iprt > 1
    endif ! iLine==start_line

    ! option to write out solar data and quit
    if (iLine == start_line) then
      if (wrt_solar) then
        call write_solar_tio(errstat)
        if (errstat < 0) then
          call tell_error (tell_io_read_error, &
               "write_solar_tio: failed", &
               errstat)
          return
        endif
      endif
    endif

  end subroutine read_input_data_tio


  subroutine read_cld_dimensions(l1bfile, tio_l1obj, swathname, errstat)
    !open netCDF radiance file and get dimensions
    use m_vars, only: nXtrack, nTimes, nWavel, iprt

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=*), intent (in) :: swathname

    !output variables
    integer (kind=4), intent (inout) :: errstat

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, "xtrack", nXtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "mirror_step", nTimes, errstat)
    call tiof_inq_dimlen (tio_l1obj, "spectral_channel", nWavel, errstat)
    if(iprt > 0) then
      print *,'read_input_data_tio: nTimes, nXtrack, nWavel'
      print *, nTimes,nXtrack,nWavel
    endif
    call tiof_close (tio_l1obj, errstat)
    
    if (errstat < 0) then
      call tell_error (tell_io_open_error, &
           "read_cld_dimensions: failed to open L1B file", &
           errstat)
      return
    endif

  end subroutine read_cld_dimensions
 

  subroutine read_cld_geo_data(l1bfile, tio_l1obj, swathname, errstat)
    !read geolocation data from netCDF input file
    use m_vars, only: iLine, time, lat, lon, sza, sazimuth, sat_zen, &
         vazimuth, terr_height, geoflg, anomflg, mflg, nLines, nXtrack, &
         azimuth, fill_value, read_he4

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=*), intent (in) :: swathname
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=8), dimension(1) :: tio_time
    real (kind=4), dimension(nXtrack,1) :: tio_lat, tio_lon, tio_sza, &
         tio_sat_zen, tio_vazimuth, tio_sazimuth
    integer (kind=2), dimension(nXtrack,1) :: tio_terr_height, tio_geoflg
    integer (kind=1), dimension(nXtrack,1) :: tio_anomflg 
    integer (kind=2), dimension(1) :: tio_mflg

    integer (kind=4) :: i

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, cld_var_time, [iLine-1], [1], &
         tio_time, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_latitude, [iLine-1,0], [1,-1], &
         tio_lat, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_longitude, [iLine-1,0], [1,-1], &
         tio_lon, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_sz_angle, [iLine-1,0], [1,-1], &
         tio_sza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "solar_azimuth_angle", [iLine-1,0], &
         [1,-1], tio_sazimuth, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_vz_angle, [iLine-1,0], [1,-1], &
         tio_sat_zen, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_azimuth_angle", [iLine-1,0], &
         [1,-1], tio_vazimuth, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ellipsoid_altitude", [iLine-1,0], &
         [1,-1], tio_terr_height, errstat)
    call tiof_get2d_i2 (tio_l1obj, "GroundPixelQualityFlags", [iLine-1,0], &
         [1,-1], tio_geoflg, errstat)
    call tiof_get2d_i1 (tio_l1obj, "XTrackQualityFlags", [iLine-1,0], &
         [1,-1], tio_anomflg, errstat)
    call tiof_get1d_i2 (tio_l1obj, "MeasurementQualityFlags", [iLine-1], &
         [1], tio_mflg, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_geo_data: failed to read geolocation data", &
           errstat)
      return
    endif

    !Test array values match he4 version
    if (read_he4) then
      if(time(iLine).ne.tio_time(1)) print *,'mismatch:time'
      if(mflg(iLine).ne.tio_mflg(1)) print *,'mismatch:mflg'
      do i=1,nXtrack
        if(lat(i,iLine).ne.tio_lat(i,1)) print *,'mismatch:lat'
        if(lon(i,iLine).ne.tio_lon(i,1)) print *,'mismatch:lon'
        if(sza(i-1,iLine).ne.tio_sza(i,1)) print *,'mismatch:sza'
        if(sazimuth(i,iLine).ne.tio_sazimuth(i,1)) print *,'mismatch:sazimuth'
        if(vazimuth(i,iLine).ne.tio_vazimuth(i,1)) print *,'mismatch:vazimuth'
        if(sat_zen(i-1,iLine).ne.tio_sat_zen(i,1)) print *,'mismatch:sat_zen'
        if(terr_height(i,iLine).ne.tio_terr_height(i,1)) &
             print *,'mismatch:terr_height'
        if(geoflg(i,iLine).ne.tio_geoflg(i,1)) print *,'mismatch:geoflg'
        if(anomflg(i,iLine).ne.tio_anomflg(i,1)) print *,'mismatch:anomflg'
      enddo
    endif

    !Copy from temporary arrays to appropriate line of main arrays
    time(iLine)=tio_time(1)
    mflg(iLine)=tio_mflg(1)
    lat(:,iLine)=tio_lat(:,1)
    lon(:,iLine)=tio_lon(:,1)
    sza(:,iLine)=tio_sza(:,1)
    sazimuth(:,iLine)=tio_sazimuth(:,1)
    sat_zen(:,iLine)=tio_sat_zen(:,1)
    vazimuth(:,iLine)=tio_vazimuth(:,1)
    terr_height(:,iLine)=tio_terr_height(:,1)
    geoflg(:,iLine)=tio_geoflg(:,1)
    anomflg(:,iLine)=tio_anomflg(:,1)
    !Calculate relative azimuth angle
    azimuth(:,iLine)=sazimuth(:,iLine)+180.0-vazimuth(:,iLine)
    where(azimuth(:,iLine) < -180.) azimuth(:,iLine)=azimuth(:,iLine)+360.
    where(azimuth(:,iLine) > 180.) azimuth(:,iLine)=azimuth(:,iLine)-360.
    azimuth(:,iLine)=abs(azimuth(:,iLine))
    where(azimuth(:,iLine) > 360.0) azimuth(:,iLine)=fill_value


  end subroutine read_cld_geo_data


  subroutine read_cld_geo_data2(l1bfile, tio_l1obj, swathname, errstat)
    !read geolocation data from netCDF input file
    use m_vars, only: iLine, time, lat, lon, sza, sazimuth, sat_zen, &
         vazimuth, terr_height, geoflg, anomflg, mflg, nLines, nXtrack, &
         azimuth, fill_value, read_he4

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=*), intent (in) :: swathname
    !output variables
    integer (kind=4), intent (inout) :: errstat

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, cld_var_time, [0], [-1], &
         time, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_latitude, [0,0], [-1,-1], &
         lat, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_longitude, [0,0], [-1,-1], &
         lon, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_sz_angle, [0,0], [-1,-1], &
         sza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "solar_azimuth_angle", [0,0], &
         [-1,-1], sazimuth, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_vz_angle, [0,0], [-1,-1], &
         sat_zen, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_azimuth_angle", [0,0], &
         [-1,-1], vazimuth, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ellipsoid_altitude", [0,0], &
         [-1,-1], terr_height, errstat)
    call tiof_get2d_i2 (tio_l1obj, "GroundPixelQualityFlags", [0,0], &
         [-1,-1], geoflg, errstat)
    call tiof_get2d_i1 (tio_l1obj, "XTrackQualityFlags", [0,0], &
         [-1,-1], anomflg, errstat)
    call tiof_get1d_i2 (tio_l1obj, "MeasurementQualityFlags", [0], &
         [-1], mflg, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_geo_data: failed to read geolocation data", &
           errstat)
      return
    endif

    !Calculate relative azimuth angle
    azimuth(:,:)=sazimuth(:,:)+180.0-vazimuth(:,:)
    where(azimuth(:,:) < -180.) azimuth(:,:)=azimuth(:,:)+360.
    where(azimuth(:,:) > 180.) azimuth(:,:)=azimuth(:,:)-360.
    azimuth(:,:)=abs(azimuth(:,:))
    where(azimuth(:,:) > 360.0) azimuth(:,:)=fill_value


  end subroutine read_cld_geo_data2



  subroutine read_cld_rad_data(l1bfile, tio_l1obj, swathname, errstat)
    !read radiance data from netCDF input file
    use m_vars, only: iLine, wmin2, wmax2, nwl, f12d, w12d, quality_flagL, &
         nWavel, nXtrack, read_he4
    use m_read_solar_data_tio, only: calc_wl_line

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=*), intent (in) :: swathname
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=4), dimension(nWavel, nXtrack, 1) :: tio_wvl, tio_rad
    real (kind=4), dimension(nWavel, nXtrack) :: wl_local
    integer (kind=2), dimension(nWavel,nXtrack, 1) :: tio_flg
    integer (kind=4) :: ih, il, i, j

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    !open file, read wavelength, radiance, flag values
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_get3d_r4 (tio_l1obj, cld_var_radiance, [iLine-1,0,0], &
         [1,-1,-1], tio_rad, errstat)
    call tiof_get3d_r4 (tio_l1obj, cld_var_wavelength, [iLine-1,0,0], &
         [1,-1,-1], tio_wvl, errstat)
    call tiof_get3d_i2 (tio_l1obj, cld_var_dqf, [iLine-1,0,0], &
         [1,-1,-1], tio_flg, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_rad_data: failed to read radiance data", &
           errstat)
      return
    endif

    !Determine limits of wavelength window
    wl_local(:,:)=tio_wvl(:,:,1)
    errstat=calc_wl_line(iLine, nXtrack, nWavel, wmin2, wmax2, wl_local, &
         il, ih, nwl)
    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
           "read_cld_rad_data: calc_wl_line: failed", &
           errstat)
      return
    endif

    if (read_he4) then
      !Test array values match he4 version
      do i=1,nXtrack
        do j=1,nwl
          if(f12d(j-1,i-1).ne.tio_rad(il+j-1,i,1)) print *,'mismatch: rad'
          !NB match to precision of wavelength values
          if(w12d(j-1,i-1)-tio_wvl(il+j-1,i,1).ge.4e-5) then
            print *,'mismatch: wvl',w12d(j-1,i-1)-tio_wvl(il+j-1,i,1)
          endif
          if(quality_flagL(j,i).ne.tio_flg(il+j-1,i,1)) print *,'mismatch: flg'
        enddo
      enddo
    endif

    !Copy radiance, wavelength and quality in window to main arrays
    f12d(0:nwl-1,0:nXtrack-1)=tio_rad(il:ih,1:nXtrack,1)
    w12d(0:nwl-1,0:nXtrack-1)=tio_wvl(il:ih,1:nXtrack,1)
    quality_flagL(1:nwl,1:nXtrack)=tio_flg(il:ih,1:nXtrack,1)

    w12d(nwl:nWavel-1,:)=0.d0
    f12d(nwl:nWavel-1,:)=0.d0

  end subroutine read_cld_rad_data






  subroutine alloc_scan(errstat)
  !allocate memory for variables 
    use m_vars, only: lat, lon, sza, sat_zen, sazimuth, vazimuth, azimuth, &
         terr_height, geoflg, anomflg, mflg, quality_flagL, w12d, f12d, &
         time, meas_qual_flg, cloud_pres, refl, dIdR, ps, ref_clr, &
         reflect_cld, eff_cld_frac, eff_cld_frac2, rad_cld_frac, cld_pres2, & 
         chlorophyll, biases, biases2, stds, stds2, chi_sqr, chi_sqr2, &
         land_flg, chlcl, qc, qc2, fill, shifts, shifts2, squeezes, &
         nXtrack, nLines, nWavel, fill_value

    implicit none
    integer (kind=4), intent (inout) :: errstat

    if (errstat < 0) return

    if (allocated(lat)) deallocate (lat)   
    allocate( lat(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: lat", &
           errstat)
      return
    endif

    if (allocated(lon)) deallocate (lon)   
    allocate( lon(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: lon", &
           errstat)
      return
    endif

    if (allocated(sza)) deallocate (sza)
    allocate( sza(0:nXtrack-1,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: sza", &
           errstat)
      return
    endif

    if (allocated(sat_zen)) deallocate (sat_zen)
    allocate( sat_zen(0:nXtrack-1,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: sat_zen", &
           errstat)
      return
    endif

    if (allocated(sazimuth)) deallocate (sazimuth)
    allocate( sazimuth(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: sazimuth", &
           errstat)
      return
    endif

    if (allocated(vazimuth)) deallocate (vazimuth)
    allocate( vazimuth(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: vazimuth", &
           errstat)
      return
    endif

    if (associated(terr_height)) nullify (terr_height)
    allocate( terr_height(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: terr_height", &
           errstat)
      return
    endif

    if (allocated(geoflg)) deallocate (geoflg)
    allocate( geoflg(nXtrack,nLines), STAT=errstat ) ; geoflg=0 
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: geoflg", &
           errstat)
      return
    endif

    if (allocated(anomflg)) deallocate (anomflg)
    allocate( anomflg(nXtrack,nLines), STAT=errstat ) ; anomflg=0
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: anomflg", &
           errstat)
      return
    endif

    if (allocated(mflg)) deallocate (mflg)
    allocate( mflg(nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: mflg", &
           errstat)
      return
    endif

    if (allocated(quality_flagL)) deallocate (quality_flagL)
    allocate( quality_flagL(nWavel,nXtrack), STAT=errstat ) ; quality_flagL=0
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: quality_flagL", &
           errstat)
      return
    endif

    if (allocated(w12d)) deallocate (w12d)
    allocate( w12d(0:nWavel-1,0:nXtrack-1), STAT=errstat )   ; w12d = 0.0
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: w12d", &
           errstat)
      return
    endif

    if (allocated(f12d)) deallocate (f12d)
    allocate( f12d(0:nWavel-1,0:nXtrack-1), STAT=errstat ) ; f12d=0.
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: f12d", &
           errstat)
      return
    endif

    if (associated(time)) nullify (time)
    allocate( time(nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: time", &
           errstat)
      return
    endif

    if (allocated(meas_qual_flg)) deallocate (meas_qual_flg)
    allocate  (meas_qual_flg(nLines)) ; meas_qual_flg = 0
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: meas_qual_flg", &
           errstat)
      return
    endif

    if (allocated(cloud_pres)) deallocate (cloud_pres)
    allocate (cloud_pres (0:nXtrack - 1,nLines)) ; cloud_pres=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: cloud_pres", &
           errstat)
      return
    endif

    if (allocated(azimuth)) deallocate (azimuth)
    allocate (azimuth (0:nXtrack-1,nLines))      ; azimuth=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: azimuth", &
           errstat)
      return
    endif

    if (allocated(refl)) deallocate (refl)  
    allocate (refl    (0:nXtrack-1,nLines))      ; refl=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: refl", &
           errstat)
      return
    endif

    if (allocated(dIdR)) deallocate (dIdR)  
    allocate (dIdR    (0:nXtrack-1,nLines))      ; dIdR=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: dIdR", &
           errstat)
      return
    endif

    if (allocated(ps)) deallocate (ps)     
    allocate (ps      (0:nXtrack-1,nLines))      ; ps=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: ps", &
           errstat)
      return
    endif

    if (allocated(ref_clr)) deallocate (ref_clr)     
    allocate (ref_clr      (0:nXtrack-1,nLines))      ; ref_clr=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: ref_clr", &
           errstat)
      return
    endif

    if (allocated(reflect_cld)) deallocate (reflect_cld)      
    allocate (reflect_cld      (0:nXtrack-1,nLines)) ; reflect_cld=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: reflect_cld", &
           errstat)
      return
    endif

    if (allocated(rad_cld_frac)) deallocate (rad_cld_frac)
    allocate (rad_cld_frac(0:nXtrack-1,nLines))  ; rad_cld_frac=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: rad_cld_frac", &
           errstat)
      return
    endif

    if (allocated(eff_cld_frac)) deallocate (eff_cld_frac)
    allocate (eff_cld_frac(0:nXtrack-1,nLines))  ; eff_cld_frac=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: eff_cld_frac", &
           errstat)
      return
    endif

    if (allocated(eff_cld_frac2)) deallocate (eff_cld_frac2)
    allocate (eff_cld_frac2(0:nXtrack-1,nLines))  ; eff_cld_frac2=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: eff_cld_frac2", &
           errstat)
      return
    endif

    if (allocated(cld_pres2)) deallocate (cld_pres2)
    allocate (cld_pres2(0:nXtrack-1,nLines))  ; cld_pres2=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: cld_pres2", &
           errstat)
      return
    endif

    if (allocated(chlorophyll)) deallocate (chlorophyll)
    allocate (chlorophyll(0:nXtrack - 1,nLines)) ; chlorophyll=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: chlorophyll", &
           errstat)
      return
    endif

    if (allocated(biases)) deallocate (biases)  
    allocate (biases  (0:nXtrack-1,nLines))      ; biases=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: biases", &
           errstat)
      return
    endif

    if (allocated(biases2)) deallocate (biases2)  
    allocate (biases2  (0:nXtrack-1,nLines))      ; biases2=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: biases2", &
           errstat)
      return
    endif

    if (allocated(stds)) deallocate (stds) 
    allocate (stds    (0:nXtrack-1,nLines))      ; stds=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: stds", &
           errstat)
      return
    endif

    if (allocated(stds2)) deallocate (stds2) 
    allocate (stds2    (0:nXtrack-1,nLines))      ; stds2=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: stds2", &
           errstat)
      return
    endif

    if (allocated(chi_sqr)) deallocate (chi_sqr)
    allocate (chi_sqr (0:nXtrack-1,nLines))      ; chi_sqr=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: chi_sqr", &
           errstat)
      return
    endif

    if (allocated(chi_sqr2)) deallocate (chi_sqr2)
    allocate (chi_sqr2 (0:nXtrack-1,nLines))      ; chi_sqr2=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: chi_qsr2", &
           errstat)
      return
    endif

    if (allocated(land_flg)) deallocate (land_flg)
    allocate (land_flg(0:nXtrack-1)) ; land_flg=.FALSE.
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: land_flg", &
           errstat)
      return
    endif

    if (allocated(chlcl)) deallocate (chlcl)
    allocate (chlcl   (0:nXtrack-1)) ; chlcl=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: chlcl", &
           errstat)
      return
    endif

    if (allocated(qc)) deallocate (qc)      
    allocate (qc      (0:nXtrack-1,nLines)) ; qc=0
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: qc", &
           errstat)
      return
    endif

    if (allocated(qc2)) deallocate (qc2)      
    allocate (qc2      (0:nXtrack-1,nLines)) ; qc2=0
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: qc2", &
           errstat)
      return
    endif

    if (allocated(fill)) deallocate (fill)      
    allocate (fill    (0:nXtrack-1,nLines)) ; fill=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: fill", &
           errstat)
      return
    endif

    if (allocated(shifts)) deallocate (shifts)      
    allocate (shifts  (0:nXtrack-1,nLines)) ; shifts=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: shifts", &
           errstat)
      return
    endif

    if (allocated(shifts2)) deallocate (shifts2)      
    allocate (shifts2  (0:nXtrack-1,nLines)) ; shifts2=fill_value
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: shifts2", &
           errstat)
      return
    endif

    if (allocated(squeezes)) deallocate (squeezes)      
    allocate (squeezes(0:nXtrack-1,nLines)) ; squeezes=1
    if (errstat < 0) then
            call tell_error (tell_malloc_error, &
           "alloc_scan: allocation failure: squeezes", &
           errstat)
      return
    endif

  end subroutine alloc_scan





end module m_read_input_data_tio
