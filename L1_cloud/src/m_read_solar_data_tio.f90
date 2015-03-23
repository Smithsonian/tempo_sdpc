!read solar data from L1B netCDF file
module m_read_solar_data_tio
  use cld_names_module
  use tio_module
  use tell_module
  use netcdf, only : nf90_nowrite

  public read_solar_data_tio, write_solar_tio, calc_wl_line
  private read_sol_dimensions, read_sol_data, read_earth_sun_dist

contains

  subroutine read_solar_data_tio(errstat)
    !read in a netCDF irradiance file
    use m_vars, only: fs, nsolwave, iprt, status, dist_rad, dist_irrad, &
         filename, nc_swathname
    use m_lambda_qual
    use m_LUN_set
    use m_pgs_include

    implicit none

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !internal variables
    integer (kind = 4) :: version
    integer (kind = 4) :: pgs_pc_getreference, nTimes, nXtrack, nWavel, &
         ext_index
    character (len = 200) :: filename_sol, filename_sol_nc, swathname, &
         filename_rad_nc

    type (tiof_file_type) :: tio_irrl1obj


    if (errstat /= 0) return

    ! obtain name of IRR1B data file
    version = 1
    status = pgs_pc_getreference( IRR1B_FILE, version, filename_sol )
    if( status .ne. PGS_S_SUCCESS ) then
      errstat = -1
      call tell_error (tell_io_read_error, &
           "read_solar_data_tio: failed to obtain irradiance L1B filename", &
           errstat)
      return
    endif
    ext_index=index(filename_sol,'.he4')
    filename_sol_nc=filename_sol(1:ext_index-1)//'.nc'

    !open IRR1B file
    swathname=trim(nc_swathname)

    if (iprt >= 2) then
      print *,'read_solar_data_tio: opening ',trim(filename_sol_nc),' ', &
           trim(swathname)
    endif
    call read_sol_dimensions(filename_sol_nc, tio_irrl1obj, swathname, &
         nTimes, nXtrack, nWavel, errstat)
    if(errstat /= 0) then 
      call tell_error (tell_io_read_error, &
           "read_sol_dimensions: failed", &
           errstat)
      return
    endif

    ! Read and test measurement quality flags 
    call read_sol_mflg(filename_sol_nc, tio_irrl1obj, swathname, errstat)
    if(errstat /= 0) then 
      call tell_error (tell_io_read_error, &
           "read_sol_mflg: ended with error", &
           errstat)
      return
    endif

    !Read in solar data
    call read_sol_data(filename_sol_nc, tio_irrl1obj, swathname, &
         nXtrack, nWavel, errstat)
    if(errstat /= 0) then 
      call tell_error (tell_io_read_error, &
           "read_sol_data: failed", &
           errstat)
      return
    endif


    !set processing quality_flags
    call bad_irrad_lambda(nXtrack)

    !correction for earth-sun distance
    call read_earth_sun_distance(filename_sol_nc,dist_irrad,errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_earth_sun_distance: failed for irradiance file", &
           errstat)
      return
    endif

    ext_index=index(filename,'.he4')
    filename_rad_nc=filename(1:ext_index-1)//'.nc'
    call read_earth_sun_distance(filename_rad_nc,dist_rad,errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_earth_sun_distance: failed for radiance file", &
           errstat)
      return
    endif

    fs(0:nsolwave-1,:)=fs(0:nsolwave-1,:)*(dist_irrad/dist_rad)**2


  end subroutine read_solar_data_tio



  subroutine read_sol_dimensions(filename_sol_nc, tio_irrl1obj, swathname, &
       nTimes, nXtrack, nWavel, errstat)
    !open netCDF irradiance file and get dimensions
    use m_vars, only: iprt

    implicit none
    !input variables
    character (len=*), intent (in) :: filename_sol_nc
    character (len=200), intent (in) :: swathname

    !output variables
    integer (kind=4), intent (inout) :: errstat
    integer (kind=4), intent (out) :: nTimes, nXtrack, nWavel

    type (tiof_file_type) :: tio_irrl1obj

    if (errstat /= 0) return

    call tiof_open (filename_sol_nc, tio_irrl1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_irrl1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_irrl1obj, cld_dim_xtrack, nXtrack, errstat)
    call tiof_inq_dimlen (tio_irrl1obj, cld_dim_step, nTimes, errstat)
    call tiof_inq_dimlen (tio_irrl1obj, cld_dim_channel, nWavel, errstat)
    if(iprt > 0) then
      print *,'read_solar_data_tio: nTimes, nXtrack, nWavel'
      print *, nTimes,nXtrack,nWavel
    endif
    call tiof_close (tio_irrl1obj, errstat)
    
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_sol_dimensions: failed to open L1B file", &
           errstat)
      return
    endif

  end subroutine read_sol_dimensions


  subroutine read_sol_data(filename_sol_nc, tio_irrl1obj, swathname, &
       nXtrack, nWavel, errstat)
    !read irradiance data from netCDF solar file
    use m_vars, only: wmin2, wmax2, ws, fs, nsolwave, iprt, ierr, &
         dist_rad, dist_irrad, irr_quality_flagL, read_he4

    implicit none
    !input variables
    character (len=*), intent (in) :: filename_sol_nc 
    character (len=200), intent (in) :: swathname
    integer (kind=4), intent (in) :: nXtrack, nWavel
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    real (kind=4), dimension(nWavel, nXtrack, 1) :: tio_wvl, tio_rad
    real (kind=4), dimension(nWavel, nXtrack) :: wl_local, tio_rad2
    integer (kind=2), dimension(nWavel,nXtrack, 1) :: tio_flg
    integer (kind=4) :: ih, il, i, j

    type (tiof_file_type) :: tio_irrl1obj

    if (errstat /= 0) return

    !open file, read wavelength, radiance, flag values
    call tiof_open (filename_sol_nc, tio_irrl1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_irrl1obj, swathname, errstat)
    call tiof_get3d_r4 (tio_irrl1obj, cld_var_irradiance, [0,0,0], &
         [1,nXtrack,nWavel], tio_rad, errstat)
    call tiof_get3d_r4 (tio_irrl1obj, cld_var_wavelength, [0,0,0], &
         [1,nXtrack, nWavel], tio_wvl, errstat)
    call tiof_get3d_i2 (tio_irrl1obj, cld_var_pixelqf, [0,0,0], &
         [1,nXtrack, nWavel], tio_flg, errstat)
    call tiof_close (tio_irrl1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_sol_data: failed to read irradiance data", &
           errstat)
      return
    endif

    !Determine limits of wavelength window
    wl_local(:,:)=tio_wvl(:,:,1)
    errstat=calc_wl_line(nXtrack, nWavel, wmin2, wmax2, wl_local, &
         il, ih, nsolwave)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_sol_data: calc_wl_line: failed", &
           errstat)
      return
    endif

    if (iprt .ge. 2) print *,'nsolwave ',nsolwave

    !MAY ALREADY BE DONE IN he4 READ
    if (.not. read_he4) then
      !allocate solar data arrays
      if (allocated(ws)) deallocate(ws)
      allocate( ws(0:nsolwave-1,0:nXtrack-1), stat=ierr )
      if (ierr .ne. 0) then
        errstat=-1
        call tell_error (tell_malloc_error, &
             "read_sol_data: failed to allocate memory: ws", &
             errstat)
        return
      endif

      if (allocated(fs)) deallocate(fs)
      allocate( fs(0:nsolwave-1,0:nXtrack-1), stat=ierr )
      if (ierr .ne. 0) then
        errstat=-1
        call tell_error (tell_malloc_error, &
             "read_sol_data: failed to allocate memory: fs", &
             errstat)
        return
      endif

      if (allocated (irr_quality_flagL)) deallocate( irr_quality_flagL)
      allocate( irr_quality_flagL(nWavel,nXtrack), stat=ierr )
      if(ierr .ne. 0) then
        errstat=-1
        call tell_error (tell_malloc_error, &
             "read_solar_data_tio: failed to allocate irr_quality_flagL", &
             errstat)
        return
      end if
      irr_quality_flagL(1:nWavel,1:nXtrack) = -1

    endif


    if (read_he4) then
      !Test array values match he4 version
      tio_rad2(il:ih,:)=tio_rad(il:ih,:,1)*(dist_irrad/dist_rad)**2
      do i=1,nXtrack
        do j=1,nsolwave
          !NB match to precision of corrected irradiances
          if(fs(j-1,i-1)-tio_rad2(il+j-1,i).ge.1e8) then
            if (iprt >= 1) print *,'mismatch: irrad',i,j,fs(j-1,i-1)-tio_rad2(il+j-1,i)
          endif
          !NB match to precision of wavelength values
          if(ws(j-1,i-1)-tio_wvl(il+j-1,i,1).ge.4e-5) then
            print *,'mismatch: solar wvl',ws(j-1,i-1)-tio_wvl(il+j-1,i,1)
          endif
          if(irr_quality_flagL(j,i).ne.tio_flg(il+j-1,i,1)) &
               print *,'mismatch: irrad flag'
        enddo
      enddo
    endif

    !Copy radiance, wavelength and quality in window to main arrays
    fs(0:nsolwave-1,0:nXtrack-1)=tio_rad(il:ih,1:nXtrack,1)
    ws(0:nsolwave-1,0:nXtrack-1)=tio_wvl(il:ih,1:nXtrack,1)
    irr_quality_flagL(1:nsolwave,1:nXtrack)=tio_flg(il:ih,1:nXtrack,1)


  end subroutine read_sol_data



  subroutine read_sol_mflg(filename_sol_nc, tio_irrl1obj, swathname, &
       errstat)
    !read measurement quality flag from netCDF solar file
    use m_vars, only: meas_qual_flg

    implicit none
    !input variables
    character (len=*), intent (in) :: filename_sol_nc 
    character (len=200), intent (in) :: swathname
    !output variables
    integer (kind=4), intent (inout) :: errstat
    !local variables
    integer (kind=2), dimension(1) :: mflg

    type (tiof_file_type) :: tio_irrl1obj

    if (errstat /= 0) return

    call tiof_open (filename_sol_nc, tio_irrl1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_irrl1obj, swathname, errstat)
    call tiof_get1d_i2 (tio_irrl1obj, "MeasurementQualityFlags", [0], &
         [1], mflg, errstat)
    call tiof_close (tio_irrl1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_sol_mflg: unable to read Measurement Quality Flags", &
           errstat)
      return
    endif

    ! testing MeasurementQualityFlags
    if(btest(mflg(1),0) .or. btest(mflg(1),1) .or. btest(mflg(1),3) &
         .or. btest(mflg(1),12) .or. btest(mflg(1),11)) then
      errstat = -1
      call tell_error (tell_runtime_error, &
           "read_solar_data_tio: solar data flagged as bad, aborting", &
           errstat)
      return
    endif

    if(btest(mflg(1),2) .or. btest(mflg(1),4) .or. btest(mflg(1),5) &
         .or. btest(mflg(1),6) .or. btest(mflg(1),7) .or. btest(mflg(1),8) &
         .or. btest(mflg(1),9) .or. btest(mflg(1),10)) &
         meas_qual_flg(:)=ibset(meas_qual_flg(:),5)

  end subroutine read_sol_mflg



  !FIXME: VIOLATES RULE THAT SUBROUTINES SHOULD NOT END MAIN PROGRAM
  subroutine write_solar_tio(errstat)
  !write out solar data to a file and quit
    use m_vars, only: solar_path, sfile, nsolwave, ws, fs

    !local variables
    integer (kind=4) :: ios, i, errstat
    character(len=80) :: sfile2

    if (errstat /= 0) return

    sfile2=trim(solar_path)//trim(sfile)
    open(1,file=sfile2,form='formatted',status='unknown',action='write',iostat=ios)
    if (ios .ne. 0) then
      errstat=-1
      call tell_error (tell_io_open_error, &
           "write_solar_tio: failed to open solar output file", &
           errstat)
      return
    else   
      print *,'opened solar file ', sfile2
      write(1,*) nsolwave
      do i=0,59
        write(1,*) ws(0:nsolwave-1,i)
        write(1,*) fs(0:nsolwave-1,i)
      enddo
      print *,'wrote solar file'
      close (1) 

      stop
      
      return
    endif

  end subroutine write_solar_tio


  subroutine read_earth_sun_distance (filename_nc, dist, errstat)
  ! read earth-sun distance from a netCDF radiance or irradiance file
    implicit none
    character (len=*), intent(in) :: filename_nc
    real (kind=4), intent(out) :: dist
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj

    if (errstat /= 0) return

    call tiof_open (filename_nc, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "Error opening file: "//trim(filename_nc), errstat)
      return
    endif

    call tiof_get_r4 (obj, cld_var_earth_sun_distance, dist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "Error reading earth-sun distance from file"//trim(filename_nc), &
                       errstat)
    endif

    call tiof_close (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, "Error closing file: "//trim(filename_nc), errstat)
      return
    endif

  end subroutine read_earth_sun_distance


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
  function calc_wl_line(nXtrack, nWavel, minwl, maxwl, wl_local, &
       il, ih, Nwl_l) result (errstat)
    real (kind = 4), intent(inout) :: minwl, maxwl
    integer (kind = 4), intent(out) :: il, ih, Nwl_l
    real (kind = 4), dimension(:,:), intent(inout) :: wl_local

    integer (kind = 4) :: errstat
    integer (kind = 4) :: fflag, i, k     

    errstat=0

    ! First check wavelengths for fill values and set fflag if found
    ! FIXME: ought to use the actual fill value, but 0-10000 is quicker
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
          call tell_error (tell_invalid_parm, &
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
          call tell_error (tell_invalid_parm, &
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
          call tell_error (tell_invalid_parm, &
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



end module m_read_solar_data_tio
