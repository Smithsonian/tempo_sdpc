!> Read meteorology parameters from netCDF simulated NAM meteorology file
module met_module
  use tell_module
  use tio_module

  public read_met_data, read_synth_met_data
  private open_synth_troppres, close_synth_troppres

contains

  !> determine whether file is netCDF or GRIB2, call appropriate subroutine
  !--------------------------------------------------------------------------
  !
  !> @param[in]  metfile  filename of netcdf synthetic meteorology file
  !> @param[in]  lat      latitude of target pixel
  !> @param[in]  lon      longitude of target pixel
  !> @param[out] troppres tropopause pressure value (hPa)
  !> @param[out] surfpres surface pressure value (hPa) [OPTIONAL]
  !> @param[out] tprof    teperature profile vector [OPTIONAL]
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Jan 2018
  !--------------------------------------------------------------------------
  subroutine read_met_data (metfile, lat, lon, troppres, surfpres, &
       tprof, errstat)

    implicit none

    character (len=*), intent(in) :: metfile
    real (kind=4), intent(in) :: lat, lon
    real (kind=4), intent(out) :: troppres
    real (kind=4), intent(out), optional :: surfpres
    real (kind=4), dimension(:), intent(out), optional :: tprof
    integer (kind=4), intent(inout) :: errstat

    integer :: ext
    logical :: read_nc, read_grib2

    if (errstat /= 0) return

    read_nc = .false.
    read_grib2 = .false.

    ext = index(metfile, '.nc')
    if (ext > 0) then
      read_nc = .true.
    else
      ext = index(metfile, '.netcdf')
      if (ext > 0) then
        read_nc = .true.
      else
        ext = index(metfile, '.NC')
        if (ext > 0) then
          read_nc = .true.
        endif
      endif
    endif

    ext = index(metfile, '.grib2')
    if (ext > 0) then
      read_grib2 = .true.
    else
      ext = index(metfile, '.grb2')
      if (ext > 0) then
        read_grib2 = .true.
      else
      ext = index(metfile, '.GRIB2')
      if (ext > 0) then
        read_grib2 = .true.
      else
        ext = index(metfile, '.GRB2')
        if (ext > 0) then
          read_grib2 = .true.
        endif
      endif
    endif
  endif

  if (read_nc .and. .not. read_grib2) then
    call read_synth_met_data (metfile, lat, lon, troppres, surfpres, &
       tprof, errstat)
  else if (read_grib2 .and. .not. read_nc) then
    print *, "GRIB2 meteorology reading not yet implemented"
    stop 1
  else
    print *, &
         "Unable to determine meteorology file type from filename extension"
    stop 1
  endif


  end subroutine read_met_data


  !> read parameters from synthetic meteorology file
  !--------------------------------------------------------------------------
  !
  !> @param[in]  metfile  filename of netcdf synthetic meteorology file
  !> @param[in]  lat      latitude of target pixel
  !> @param[in]  lon      longitude of target pixel
  !> @param[out] troppres tropopause pressure value (hPa)
  !> @param[out] surfpres surface pressure value (hPa) [OPTIONAL]
  !> @param[out] tprof    teperature profile vector [OPTIONAL]
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Jan 2018
  !--------------------------------------------------------------------------
  subroutine read_synth_met_data (metfile, lat, lon, troppres, surfpres, &
       tprof, errstat)

    implicit none

    ! number of pressure levels (fixed)
    integer(kind=4) , parameter :: nlev=72

    !Input variables
    character (len=*), intent(in) :: metfile
    real (kind=4), intent(in) :: lon, lat

    !Output variables
    real(kind=4), intent(out) :: troppres
    real(kind=4), intent(out), optional :: surfpres
    real(kind=4), intent(out), dimension(nlev), optional :: tprof
    integer(kind=4), intent(inout) :: errstat

    !local variables
    integer(kind=4) :: nlon, nlat
    integer(kind=4), dimension(2) :: lonlatidx
    real(kind=4), dimension(:,:), allocatable :: longrid, latgrid
    real(kind=4), dimension(1,1) :: tmp_surfpres, tmp_troppres
    real(kind=4), dimension(1,1,nlev) :: tmp_tprof

    type (tiof_file_type) :: metobj

    if (errstat /= 0) return


    !open the file
    call open_synth_troppres(metfile, metobj, errstat)
    if (errstat /= 0) return

    !determine dimensions of meteorology grid
    call tiof_inq_dimlen (metobj, "x", nlon, errstat)
    call tiof_inq_dimlen (metobj, "y", nlat, errstat)

    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_troppres: failed to read dimensions", errstat)
      return
    endif

    !allocate arrays where necessary
    allocate(longrid(0:nlon-1,0:nlat-1), latgrid(0:nlon-1,0:nlat-1), &
         stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "read_synth_troppres: failed to allocate arrays", errstat)
      return
    endif

    !read in lat, lon values
    call tiof_get2d_r4(metobj, "lon", [0,0], [nlat,nlon], &
         longrid(0:nlon-1,0:nlat-1), errstat)
    call tiof_get2d_r4(metobj, "lat", [0,0], [nlat,nlon], &
         latgrid(0:nlon-1,0:nlat-1), errstat)
    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_troppres: failed to read lon, lat positions", errstat)
      return
    endif

    !determine the nearest location to the target lon, lat
    lonlatidx = minloc(abs(latgrid-lat)+abs(longrid-lon))
    !correct each index by -1 since we are using 0-indexed arrays
    lonlatidx = lonlatidx-1

    !read in values to be output, convert pressures from Pa to hPa
    call tiof_get2d_r4(metobj, "TROPPB", [lonlatidx(2),lonlatidx(1)], [1,1], &
         tmp_troppres, errstat)
    troppres = real(tmp_troppres(1,1)/100.0d0, kind=4)
    if (present(surfpres)) then
      call tiof_get2d_r4(metobj, "PS", [lonlatidx(2),lonlatidx(1)], [1,1], &
           tmp_surfpres, errstat)
      surfpres = real(tmp_surfpres(1,1)/100.0d0, kind=4)
    endif
    if (present(tprof)) then
      call tiof_get3d_r4(metobj, "T", [0,lonlatidx(2),lonlatidx(1)], &
           [nlev,1,1], tmp_tprof, errstat)
      tprof = tmp_tprof(1,1,:)
    endif
    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_troppres: failed to read pressure, temp", errstat)
      return
    endif

    !close the file
    call close_synth_troppres(metobj, errstat)
    if (errstat /= 0) return

    !deallocate arrays
    if (allocated(longrid)) deallocate(longrid, stat=errstat)
    if (allocated(latgrid)) deallocate(latgrid, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "read_synth_troppres: failed to deallocate arrays", errstat)
      return
    endif

  end subroutine read_synth_met_data




  !> open netCDF synthetic meteorology file
  !--------------------------------------------------------------------------
  !
  !> @param[in]  metfile  filename of netcdf synthetic meteorology file
  !> @param[out] metobj
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Jan 2018
  !--------------------------------------------------------------------------
  subroutine open_synth_troppres (metfile, metobj, errstat)

    use netcdf, only: nf90_nowrite

    implicit none

    !Input variables
    character(len=*), intent(in) :: metfile
    !Output variables
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: metobj

    if (errstat /= 0) return

    call tiof_open (metfile, metobj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "open_synth_troppres: failed to open synth meteorology file", &
           errstat)
      return
    endif

  end subroutine open_synth_troppres


  !> close netCDF synthetic meteorology file
  !--------------------------------------------------------------------------
  !
  !> @param[in]  metobj
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Jan 2018
  !--------------------------------------------------------------------------
  subroutine close_synth_troppres (metobj, errstat)

    implicit none

    !Output variables
    integer, intent(inout) :: errstat
    !Input variables

    type (tiof_file_type) :: metobj

    if (errstat /= 0) return

    call tiof_close (metobj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "close_synth_troppres: failed to close synth meteorology file", &
           errstat)
      return
    endif

  end subroutine close_synth_troppres




end module met_module
