!> Read meteorology parameters from real or simulated NAM meteorology files
module met_module
  use tell_module
  use tio_module
  use eccodes

  public read_met_data
  private open_synth_troppres, close_synth_troppres, make_grib2_index, &
       close_grib2_index, read_synth_met_data, read_grib2_met_data

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
    call read_grib2_met_data (metfile, lat, lon, troppres, surfpres, &
       tprof, errstat)
  else
    call tell_error (tell_runtime_error, &
         "Unable to determine meteorology file type from filename extension",&
         errstat)
    return
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



  !> Subroutine to create an index of a grib2 file using eccodes
  !-------------------------------------------------------------------------
  !
  !> @param[in]  gribfile filename of GRIB2 meteorology file
  !> @param[out] grib_idx index id for GRIB2 meteorology file
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Nov 2018
  !-------------------------------------------------------------------------
  subroutine make_grib2_index (gribfile, grib_idx, errstat)

    implicit none

    ! input variables
    character (len=*), intent(in) :: gribfile

    ! output variables
    integer (kind=4), intent(out) :: grib_idx
    integer (kind=4), intent(inout) :: errstat

    ! local variables
    integer (kind=4) :: status
    !integer (kind=4) ::shortnamesize, levelsize
    !character(len=20), dimension(:), allocatable :: shortname, level

    if (errstat /= 0) return

    ! index on keywords needed for surface & tropo pressure, t profile
    call codes_index_create(grib_idx, gribfile, 'shortName,level',status)

    if (status /= CODES_SUCCESS) then
      call tell_error (tell_io_error, &
           "make_grib2_index: failed to index GRIB2 file", errstat)
      return
    endif

    !uncomment to print shortname and level keywords
    !call codes_index_get_size(grib_idx,'shortName',shortnamesize)
    !call codes_index_get_size(grib_idx,'level',levelsize)
    !allocate(shortname(shortnamesize), level(levelsize))
    !call codes_index_get(grib_idx,'shortName',shortname)
    !call codes_index_get(grib_idx,'level',level)
    !print *, shortname, level

  end subroutine make_grib2_index



  !> Subroutine to close the index of a grib2 file after use
  !-------------------------------------------------------------------------
  !
  !> @param[in]  grib_idx index id for GRIB2 meteorology file
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Nov 2018
  !-------------------------------------------------------------------------
  subroutine close_grib2_index (grib_idx, errstat)

    implicit none

    ! input variables
    integer (kind=4), intent(out) :: grib_idx

    ! output variables
    integer (kind=4), intent(inout) :: errstat

    ! local variables
    integer (kind=4) :: status

    if (errstat /= 0) return

    ! index on keywords needed for surface & tropo pressure, t profile
    call codes_index_release(grib_idx, status)

    if (status /= CODES_SUCCESS) then
      call tell_error (tell_io_error, &
           "close_grib2_index: failed to close GRIB2 index", errstat)
      return
    endif

  end subroutine close_grib2_index


  !> Subroutine to read from an indexed GRIB2 meteorology file using eccodes
  !-------------------------------------------------------------------------
  !
  !> @param[in]  gribfile filename for GRIB2 meteorology file
  !> @param[in]  lat      latitude  of target pixel
  !> @param[in]  lon      longitude of target pixel
  !> @param[out] troppres tropopause pressure value
  !> @param[out] surfpres surface pressure value [OPTIONAL]
  !> @param[out] tprof    temperature profile vector [OPTIONAL]
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Nov 2018
  !-------------------------------------------------------------------------
  subroutine read_grib2_met_data (gribfile, lat, lon, troppres, &
       surfpres, tprof, errstat)

    implicit none

    !FIXME - pressure profile number of levels and pressure values
    !        hard coded. Tough to see how to avoid this for GRIB2.
    !        Need to figure out how many / which levels to read
    !number of pressure levels in temperature profile
    integer (kind=4), parameter :: npres=10

    !input variables
    character (len=*), intent(in) :: gribfile
    real (kind=4), intent(in) :: lon, lat

    !output variables
    real (kind=4), intent(out), optional :: surfpres
    real (kind=4), intent(out) :: troppres
    real (kind=4), dimension(npres), intent(out), optional :: tprof
    integer (kind=4), intent(inout) :: errstat

    !local variables
    integer (kind=4) :: grib_idx
    integer (kind=4) :: igrib, status, i
    real (kind=8) :: inlon, inlat
    real (kind=8), dimension(4) :: outlats, outlons, distance, weights
    real (kind=8), dimension(4) :: tmp_spres, tmp_tpres, tmp_t
    integer (kind=4), dimension(4) :: kindex
    integer (kind=4), dimension (npres), parameter :: presvals = (/1000, &
         900, 800, 700, 600, 500, 400, 300, 200, 100/)
    real (kind=8) :: weightsum

    if (errstat /= 0) return

    call make_grib2_index(gribfile, grib_idx, errstat)
    if (errstat /= 0) return

    inlat = real(lat, kind=8)
    if (lon < 0.0d0) inlon = real(lon, kind=8)+360.0d0

    ! tropopause pressure in 4 pixels nearest lon, lat position
    call codes_index_select(grib_idx,'shortName','pres',status)
    call codes_index_select(grib_idx,'level','0',status)
    call codes_new_from_index(grib_idx, igrib, status)
    call codes_grib_find_nearest_four_single(igrib, .false., &
         inlat, inlon, outlats, outlons, tmp_tpres, distance, kindex, status)
    if (status /= CODES_SUCCESS) then
      call tell_error (tell_io_error, &
           "read_grib2_met_data: failed to read tropopause pressure", errstat)
      return
    endif
    ! weight by 1/distance^2
    weightsum = sum(1/distance**2)
    weights = (1/distance**2)/weightsum
    troppres = real(sum(tmp_tpres*weights)/100.0d0, kind=4)

    ! surface pressure
    if (present(surfpres)) then
      call codes_release(igrib)
      call codes_index_select(grib_idx,'shortName','sp',status)
      call codes_index_select(grib_idx,'level','0',status)
      call codes_new_from_index(grib_idx, igrib, status)
      call codes_get_element(igrib,"values", kindex, tmp_spres)
      surfpres = real(sum(tmp_spres*weights)/100.0d0, kind=4)
      if (status /= CODES_SUCCESS) then
        call tell_error (tell_io_error, &
             "read_grib2_met_data: failed to read surface error", errstat)
        return
      endif
    endif

    ! loop over pressure levels to get temperature profile
    if (present(tprof)) then
      do i=1, npres
        call codes_release(igrib)
        call codes_index_select(grib_idx,'shortName','t',status)
        call codes_index_select(grib_idx,'level',presvals(i),status)
        call codes_new_from_index(grib_idx, igrib, status)
        call codes_get_element(igrib,"values", kindex, tmp_t)
        tprof(i) = real(sum(tmp_t*weights), kind=4)
        if (status /= CODES_SUCCESS) then
          call tell_error (tell_io_error, &
               "read_grib2_met_data: failed to read temperature", errstat)
          return
        endif
      enddo
    endif

    !clean up
    call codes_release(igrib, status)

    call close_grib2_index(grib_idx, errstat)
    if (errstat /= 0) return

  end subroutine read_grib2_met_data




end module met_module
