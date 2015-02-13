!read input data from L1B netCDF file
module m_read_input_data_tio
  use cld_names_module
  use tio_module
  use tell_module
  use netcdf, only : nf90_nowrite

  public read_input_data_tio
  private read_cld_dimensions, read_cld_geo_data, read_cld_rad_data, &
       alloc_scan, calc_wl_line

contains

  subroutine read_input_data_tio(l1bfile, errstat)
    !read in a netCDF radiance file
    use m_vars, only: iprt, input_data_path, iLine, wrt_solar, wmin, wmax, &
         wmin2, wmax2, set_wmin, set_wmax, wave_long, wave_short, nLines, &
         start_line, max_lines, nXtrack, nTimes, nWavel, meas_qual_flg, &
         mflg, n_input, n_missing, nwl, w12d, qc, nwave, min_wl
    use m_strpos

    implicit none

    !input variables
    character (len=*), intent (in) :: l1bfile 

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !local variables
    character (len=200) :: filenamepath, swathname
    integer (kind = 4) :: PGS_TD_TAItoUTC 
    logical :: uvswath

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    if (iLine == 0) then
      !-----------------------------------------------------------------
      !If on first line, perform some setup
      !set file path+name, swath name
      filenamepath=trim(input_data_path)//l1bfile
      uvswath = strpos (l1bfile, 'BRUG') > 0
      if (iprt > 0) print *,'read_input_data_tio: filename ', &
           trim(filenamepath), '   ', trim(swathname)
      if (uvswath) then
        swathname = "band_540_740_nm"
      else
        call tell_error (tell_io_write_error, &
             "read_input_data_tio: input file is not OMI L1 UV swath", &
             errstat)
        return
      endif
      !ALREADY DONE IN he4 READ
      !set wavelength bounds
      !   if (wrt_solar) then
      !     wmin=355.0d0
      !     wmax=500.0d0
      !     wmin2=310.0d0
      !     wmax2=375.0d0
      !   endif
      !   wmin2 = 330.0d0
      !   wmax2 = 367.0d0
      !   if (.not. set_wmin) wmin = 345.5d0
      !   if (.not. set_wmax) wmax = 354.5d0
      !   wave_long=362.5d0 
      !   wave_short=345.4d0

      !open file, read variable dimension sizes
      call read_cld_dimensions(l1bfile, tio_l1obj, swathname, errstat)
      if(errstat < 0) then 
        call tell_error (tell_io_write_error, &
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

      !COMMENTED OUT AS ALREADY DONE DURING he4 READ
      !    !Allocate arrays for variables to be read in
      !    call alloc_scan(errstat)
      !    if(errstat < 0) then 
      !      call tell_error (tell_io_write_error, &
      !           "alloc_scan: failed", &
      !           errstat)
      !      return
      !    endif


    endif !iLine == 0


    !----------------------------------------------------------------
    ! For each line of the input file

    !Read in geolocation data 
    call read_cld_geo_data(l1bfile, tio_l1obj, swathname, errstat)
    if(errstat < 0) then 
      call tell_error (tell_io_write_error, &
           "read_cld_geo_data: failed", &
           errstat)
      return
    endif
    if(iprt >=2) print *,'read_cld_geo_data: success'

    !Keep track of number of value input
    n_input = n_input + nXtrack


    ! COMMENTED AS UNNECESSARY UNTIL WE SORT OUT OUR FLAGGING SYSTEM
    !  !Check data quality flags, if data missing set radiance error flag
    !  if(btest(mflg(iLine),0) .or. btest(mflg(iLine),1) &
    !       .or. btest(mflg(iLine),3) .or. btest(mflg(iLine),12)) then
    !    meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),1)
    !    n_missing = n_missing + nXtrack 
    !    if (iprt >= 1) then
    !      print *,'missing line ',iLine, btest(mflg(iLine),0), &
    !           btest(mflg(iLine),1), btest(mflg(iLine),3), btest(mflg(iLine),12)
    !    endif
    !  endif
    !  !Set internal measurement quality flags
    !  if(btest(mflg(iLine),2) .or. btest(mflg(iLine),4) &
    !       .or. btest(mflg(iLine),5) .or. btest(mflg(iLine),6) &
    !       .or. btest(mflg(iLine),8) .or. btest(mflg(iLine),9) &
    !       .or. btest(mflg(iLine),11)) &
    !       meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),2)
    !
    !  if(btest(mflg(iLine),7)) meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),3)
    !  if(btest(mflg(iLine),10)) meas_qual_flg(iLine)=ibset(meas_qual_flg(iLine),4)

    !Read in radiance data 
    call read_cld_rad_data(l1bfile, tio_l1obj, swathname, errstat)
    if(errstat < 0) then 
      call tell_error (tell_io_write_error, &
           "read_cld_rad_data: failed", &
           errstat)
      return
    endif
    if(iprt >=1) print *,'read_cld_rad_data: success'

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


  end subroutine read_input_data_tio


  subroutine read_cld_dimensions(l1bfile, tio_l1obj, swathname, errstat)
    !open netCDF radiance file and get dimensions
    use m_vars, only: nXtrack, nTimes, nWavel, iprt

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=200), intent (in) :: swathname

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
      call tell_error (tell_io_write_error, &
           "read_cld_dimensions: failed to open L1B file", &
           errstat)
      return
    endif

  end subroutine read_cld_dimensions
 

  subroutine read_cld_geo_data(l1bfile, tio_l1obj, swathname, errstat)
    !read geolocation data from netCDF input file
    use m_vars, only: iLine, time, lat, lon, sza, sazimuth, sat_zen, &
         vazimuth, terr_height, geoflg, anomflg, mflg, nLines, nXtrack, &
         azimuth, fill_value

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=200), intent (in) :: swathname
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=8), dimension(1) :: tio_time
    real (kind=4), dimension(nXtrack,1) :: tio_lat, tio_lon, tio_sza, &
         tio_sat_zen, tio_vazimuth, tio_sazimuth
    integer (kind=2), dimension(nXtrack,1) :: tio_terr_height, tio_geoflg
    integer (kind=1), dimension(nXtrack,1) :: tio_anomflg 

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
    !    call tiof_get1d_i2 (tio_l1obj, "MeasurementQualityFlags", [0,0], &
    !         [1,-1], mflg, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "read_cld_geo_data: failed to read geolocation data", &
           errstat)
      return
    endif

    !Test array values match he4 version
    if(time(iLine).ne.tio_time(1)) print *,'mismatch:time'
    do i=1,nXtrack
      if(lat(i,iLine).ne.tio_lat(i,1)) print *,'mismatch:lat'
      if(lon(i,iLine).ne.tio_lon(i,1)) print *,'mismatch:lon'
      if(sza(i-1,iLine).ne.tio_sza(i,1)) print *,'mismatch:sza'
      if(sazimuth(i,iLine).ne.tio_sazimuth(i,1)) print *,'mismatch:sazimuth'
      if(vazimuth(i,iLine).ne.tio_vazimuth(i,1)) print *,'mismatch:vazimuth'
      if(sat_zen(i-1,iLine).ne.tio_sat_zen(i,1)) print *,'mismatch:sat_zen'
      if(terr_height(i,iLine).ne.tio_terr_height(i,1)) print *,'mismatch:terr_height'
      if(geoflg(i,iLine).ne.tio_geoflg(i,1)) print *,'mismatch:geoflg'
      if(anomflg(i,iLine).ne.tio_anomflg(i,1)) print *,'mismatch:anomflg'
    enddo

    !Copy from temporary arrays to appropriate line of main arrays
    time(iLine)=tio_time(1)
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


  subroutine read_cld_rad_data(l1bfile, tio_l1obj, swathname, errstat)
    !read radiance data from netCDF input file
    use m_vars, only: iLine, wmin2, wmax2, nwl, f12d, w12d, quality_flagL, &
         nWavel, nXtrack

    implicit none
    !input variables
    character (len=*), intent (in) :: l1bfile 
    character (len=200), intent (in) :: swathname
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
      call tell_error (tell_io_write_error, &
           "read_cld_rad_data: failed to read radiance data", &
           errstat)
      return
    endif

    !Determine limits of wavelength window
    wl_local(:,:)=tio_wvl(:,:,1)
    errstat=calc_wl_line(iLine, nXtrack, nWavel, wmin2, wmax2, wl_local, &
         il, ih, nwl)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "read_cld_rad_data: calc_wl_line: failed", &
           errstat)
      return
    endif

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
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: lat", &
           errstat)
      return
    endif

    if (allocated(lon)) deallocate (lon)   
    allocate( lon(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: lon", &
           errstat)
      return
    endif

    if (allocated(sza)) deallocate (sza)
    allocate( sza(0:nXtrack-1,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: sza", &
           errstat)
      return
    endif

    if (allocated(sat_zen)) deallocate (sat_zen)
    allocate( sat_zen(0:nXtrack-1,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: sat_zen", &
           errstat)
      return
    endif

    if (allocated(sazimuth)) deallocate (sazimuth)
    allocate( sazimuth(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: sazimuth", &
           errstat)
      return
    endif

    if (allocated(vazimuth)) deallocate (vazimuth)
    allocate( vazimuth(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: vazimuth", &
           errstat)
      return
    endif

    if (associated(terr_height)) nullify (terr_height)
    allocate( terr_height(nXtrack,nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: terr_height", &
           errstat)
      return
    endif

    if (allocated(geoflg)) deallocate (geoflg)
    allocate( geoflg(nXtrack,nLines), STAT=errstat ) ; geoflg=0 
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: geoflg", &
           errstat)
      return
    endif

    if (allocated(anomflg)) deallocate (anomflg)
    allocate( anomflg(nXtrack,nLines), STAT=errstat ) ; anomflg=0
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: anomflg", &
           errstat)
      return
    endif

  !  if (allocated(mflg)) deallocate (mflg)
  !  allocate( mflg(nLines), STAT=errstat )
  !  if (errstat < 0) then
  !          call tell_error (tell_io_write_error, &
  !         "alloc_scan: allocation failure: mflg", &
  !         errstat)
  !    return
  !  endif

    if (allocated(quality_flagL)) deallocate (quality_flagL)
    allocate( quality_flagL(nWavel,nXtrack), STAT=errstat ) ; quality_flagL=0
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: quality_flagL", &
           errstat)
      return
    endif

    if (allocated(w12d)) deallocate (w12d)
    allocate( w12d(0:nWavel-1,0:nXtrack-1), STAT=errstat )   ; w12d = 0.0
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: w12d", &
           errstat)
      return
    endif

    if (allocated(f12d)) deallocate (f12d)
    allocate( f12d(0:nWavel-1,0:nXtrack-1), STAT=errstat ) ; f12d=0.
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: f12d", &
           errstat)
      return
    endif

    if (associated(time)) nullify (time)
    allocate( time(nLines), STAT=errstat )
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: time", &
           errstat)
      return
    endif

    if (allocated(meas_qual_flg)) deallocate (meas_qual_flg)
    allocate  (meas_qual_flg(nLines)) ; meas_qual_flg = 0
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: meas_qual_flg", &
           errstat)
      return
    endif

    if (allocated(cloud_pres)) deallocate (cloud_pres)
    allocate (cloud_pres (0:nXtrack - 1,nLines)) ; cloud_pres=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: cloud_pres", &
           errstat)
      return
    endif

    if (allocated(azimuth)) deallocate (azimuth)
    allocate (azimuth (0:nXtrack-1,nLines))      ; azimuth=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: azimuth", &
           errstat)
      return
    endif

    if (allocated(refl)) deallocate (refl)  
    allocate (refl    (0:nXtrack-1,nLines))      ; refl=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: refl", &
           errstat)
      return
    endif

    if (allocated(dIdR)) deallocate (dIdR)  
    allocate (dIdR    (0:nXtrack-1,nLines))      ; dIdR=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: dIdR", &
           errstat)
      return
    endif

    if (allocated(ps)) deallocate (ps)     
    allocate (ps      (0:nXtrack-1,nLines))      ; ps=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: ps", &
           errstat)
      return
    endif

    if (allocated(ref_clr)) deallocate (ref_clr)     
    allocate (ref_clr      (0:nXtrack-1,nLines))      ; ref_clr=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: ref_clr", &
           errstat)
      return
    endif

    if (allocated(reflect_cld)) deallocate (reflect_cld)      
    allocate (reflect_cld      (0:nXtrack-1,nLines)) ; reflect_cld=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: reflect_cld", &
           errstat)
      return
    endif

    if (allocated(rad_cld_frac)) deallocate (rad_cld_frac)
    allocate (rad_cld_frac(0:nXtrack-1,nLines))  ; rad_cld_frac=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: rad_cld_frac", &
           errstat)
      return
    endif

    if (allocated(eff_cld_frac)) deallocate (eff_cld_frac)
    allocate (eff_cld_frac(0:nXtrack-1,nLines))  ; eff_cld_frac=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: eff_cld_frac", &
           errstat)
      return
    endif

    if (allocated(eff_cld_frac2)) deallocate (eff_cld_frac2)
    allocate (eff_cld_frac2(0:nXtrack-1,nLines))  ; eff_cld_frac2=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: eff_cld_frac2", &
           errstat)
      return
    endif

    if (allocated(cld_pres2)) deallocate (cld_pres2)
    allocate (cld_pres2(0:nXtrack-1,nLines))  ; cld_pres2=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: cld_pres2", &
           errstat)
      return
    endif

    if (allocated(chlorophyll)) deallocate (chlorophyll)
    allocate (chlorophyll(0:nXtrack - 1,nLines)) ; chlorophyll=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: chlorophyll", &
           errstat)
      return
    endif

    if (allocated(biases)) deallocate (biases)  
    allocate (biases  (0:nXtrack-1,nLines))      ; biases=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: biases", &
           errstat)
      return
    endif

    if (allocated(biases2)) deallocate (biases2)  
    allocate (biases2  (0:nXtrack-1,nLines))      ; biases2=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: biases2", &
           errstat)
      return
    endif

    if (allocated(stds)) deallocate (stds) 
    allocate (stds    (0:nXtrack-1,nLines))      ; stds=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: stds", &
           errstat)
      return
    endif

    if (allocated(stds2)) deallocate (stds2) 
    allocate (stds2    (0:nXtrack-1,nLines))      ; stds2=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: stds2", &
           errstat)
      return
    endif

    if (allocated(chi_sqr)) deallocate (chi_sqr)
    allocate (chi_sqr (0:nXtrack-1,nLines))      ; chi_sqr=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: chi_sqr", &
           errstat)
      return
    endif

    if (allocated(chi_sqr2)) deallocate (chi_sqr2)
    allocate (chi_sqr2 (0:nXtrack-1,nLines))      ; chi_sqr2=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: chi_qsr2", &
           errstat)
      return
    endif

    if (allocated(land_flg)) deallocate (land_flg)
    allocate (land_flg(0:nXtrack-1)) ; land_flg=.FALSE.
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: land_flg", &
           errstat)
      return
    endif

    if (allocated(chlcl)) deallocate (chlcl)
    allocate (chlcl   (0:nXtrack-1)) ; chlcl=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: chlcl", &
           errstat)
      return
    endif

    if (allocated(qc)) deallocate (qc)      
    allocate (qc      (0:nXtrack-1,nLines)) ; qc=0
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: qc", &
           errstat)
      return
    endif

    if (allocated(qc2)) deallocate (qc2)      
    allocate (qc2      (0:nXtrack-1,nLines)) ; qc2=0
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: qc2", &
           errstat)
      return
    endif

    if (allocated(fill)) deallocate (fill)      
    allocate (fill    (0:nXtrack-1,nLines)) ; fill=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: fill", &
           errstat)
      return
    endif

    if (allocated(shifts)) deallocate (shifts)      
    allocate (shifts  (0:nXtrack-1,nLines)) ; shifts=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: shifts", &
           errstat)
      return
    endif

    if (allocated(shifts2)) deallocate (shifts2)      
    allocate (shifts2  (0:nXtrack-1,nLines)) ; shifts2=fill_value
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: shifts2", &
           errstat)
      return
    endif

    if (allocated(squeezes)) deallocate (squeezes)      
    allocate (squeezes(0:nXtrack-1,nLines)) ; squeezes=1
    if (errstat < 0) then
            call tell_error (tell_io_write_error, &
           "alloc_scan: allocation failure: squeezes", &
           errstat)
      return
    endif

  end subroutine alloc_scan


!! Private function: calc_wl_line
 !
 ! Functionality:
 !
 !   This function calculates the wavelengths values for a specific line
 !   and returns those values, plus the limits of the specified wavelength
 !   range.
 !
 ! Calling Arguments:
 !
 ! Inputs:
 !
 !    i         Index into L1B structure, line number
 !    this      L1B block
 !    minwl     minimum wavelength requested
 !    maxwl     maximum wavelength requested
 !
 ! Outputs:
 !
 !    wl_local  wavelength values calculated from L1B coefficients
 !    il        index of lower bound of wavelength range (in wl_local)
 !    ih        index of upper bound of wavelength range (in wl_local)
 !    Nwl_l     number of wavelengths in range
 !
 !    status    the return PGS_SMF status value
 !
 ! Change History:
 !
 !    Date            Author          Modifications
 !    ====            ======          =============
 !    January 2005    Jeremy Warner   Original Source
 !    February 2014   E. O'Sullivan   Adapted for use with TEMPO
 !
!!
  function calc_wl_line(j, nXtrack, nWavel, minwl, maxwl, wl_local, &
       il, ih, Nwl_l) result (errstat)
    integer (kind = 4), intent(in) :: j
    real (kind = 4), intent(inout) :: minwl, maxwl
    integer (kind = 4), intent(out) :: il, ih, Nwl_l
    real (kind = 4), dimension(:,:), intent(inout) :: wl_local

    integer (kind = 4) :: errstat
    integer (kind = 4) :: fflag, i, k, q        
    integer (kind = 4) :: fac


    ! First check wavelengths for fill values and set fflag if found
    ! FIX ME: ought to use the actual fill value, but 0-10000 is quicker
    fflag = 0
    do i = 1, nXtrack
      do k = 1, nWavel
        if (wl_local(k,i) .lt. 0.0 .or. wl_local(k,i) .gt. 10000.0) then
          fflag = 1
       endif
      enddo
    enddo

    ! Now find the limits of the wavelength range.

    if (fflag .eq. 0) then  ! We found no fill values, so do it the fast way.

      ! Check to see if minwl and maxwl are consistent (i.e., minwl < maxwl)

      if(minwl .gt. 0.0 .and. maxwl .gt. 0.0 .and. minwl .gt. maxwl) then
        if((minwl-maxwl) > 0.0) then
          errstat=-1
          call tell_error (tell_io_write_error, &
               "calc_wl_line: wlmin > wlmax", &
               errstat)
          return
        endif
      endif

      ! Check to see if minwl is compatable with the wavelength range
      ! If minwl was not given, set the lower bound to the first element.

      if(minwl .gt. 0.0) then
        if(minwl > maxval(wl_local))then
          errstat=-1
          call tell_error (tell_io_write_error, &
               "calc_wl_line: wlmin out of bound", &
               errstat)
          return
        endif

        ! Set the lower bound based on wl_local and minwl

        if(minwl < minval(wl_local)) then
          il = 1
        else
          il = 0
          do i = 1, nXtrack
            do k = 1, nWavel
              if (wl_local(k,i) > minwl) then
                if (il .eq. 0 .or. il .gt. k) il = k
                exit
              endif
            enddo
          enddo
          if (il > 1) then 
            il = il - 1
          else 
            il = 1
          endif
        endif
      else
        il = 1
      endif

      ! Check to see if maxwl is compatable with the wavelength range
      ! If maxwl was not given, set the upper bound to the last element.

      if(maxwl .gt. 0.0) then
        if(maxwl < minval(wl_local)) then
          errstat=-1
          call tell_error (tell_io_write_error, &
               "calc_wl_line: wlmax out of bound", &
               errstat)
          return
        endif

        ! Set the upper bound based on wl_local and maxwl

        if(maxwl > maxval(wl_local)) then
          ih = nWavel
        else
          ih = nWavel+1
          do i = 1, nXtrack
            do k = 1, nWavel
              if (wl_local(k,i) > maxwl) then
                if (ih .eq. nWavel+1 .or. ih .lt. k) ih = k
                exit
              endif
            enddo
          enddo
          if (ih < nWavel) then 
            ih = ih + 1
          else 
            ih = nWavel
          endif
        endif
      else
        ih = nWavel
      endif

    else  ! => we had a fill value

      ! Initialize il and ih

      ih = 0
      il = 0

      ! Loop through all wavelengths, and defind bounds based on minwl and maxwl

      if (minwl < 0.0) minwl = 0.0
      if (maxwl < 0.0) maxwl = 1000.0
      do k = 1, nWavel
        do i = 1, nXtrack
          if (wl_local(k,i) >= minwl .and. wl_local(k,i) <= maxwl) then
            if (il .eq. 0) il = k
            if (ih .eq. 0) ih = k
            if (ih .lt. k) ih = k
          endif
        enddo
      enddo
      if (il .eq. 0 .and. ih .eq. 0) then
      else
        if (il > 1) il = il - 1
        if (ih < nWavel) ih = ih + 1
      endif
    endif

    ! Set the number of wavelengths in the span to Nwl_l

    if (il .eq. 0 .and. ih .eq. 0) then
      Nwl_l = 0
    else
      Nwl_l = ih-il+1
    endif
    return
  end function calc_wl_line



end module m_read_input_data_tio
