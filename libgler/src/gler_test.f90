program gler_test
  use, intrinsic :: iso_c_binding, only : c_char, c_null_char
  use gler_module
  use tio_module
  use tell_module

  character (kind=c_char, len=1024) :: land_glob, ocean_glob
  integer :: iwavelen, errstat, nlon, nlat, i,j
  type(gler_type) :: glt
  real (kind=8) :: taix
  real (kind=4), allocatable, dimension(:) :: lon, lat
  real (kind=4), allocatable, dimension(:,:) :: alb
  real (kind=4) :: wind_speed, lon_range(2), lat_range(2)

  errstat = 0
  iwavelen = 466

  call tell_open ("gler_test", 0)
  call tell_set_log_level (0)

  write (land_glob, '(''$SDPC_REFDATA_DIR/gler/climatology_'',i3,''nm_fixup/modis_land_d???.nc'')') iwavelen
  write (ocean_glob, '(''$SDPC_REFDATA_DIR/gler/climatology_'',i3,''nm_fixup/modis_ocean_d???.nc'')') iwavelen

  taix = 425963157.845324
  call tiof_time_set_taix_epoch ('2000-01-01T12:00:00Z', errstat)

  call gler_open (glt, trim(land_glob)//c_null_char, trim(ocean_glob)//c_null_char, errstat)
  if (errstat /= 0) then
    write(*,*)'gler_open failed'
    stop 1
  endif

  call gler_interp_time (glt, taix, errstat)
  if (errstat /= 0) then
    write(*,*)'gler_interp_time failed'
    stop 1
  endif

  nlon = 1500
  nlat = 1000
  lon_range(:) = (/-135.0, -50.0/)
  lat_range(:) = (/  15.0,  60.0/)
  allocate (lon(nlon), lat(nlat), alb (nlon, nlat))

  lon(1:nlon) = lon_range(1) + (/(ilon, ilon=0,nlon-1)/) * (lon_range(2)-lon_range(1))/(nlon - 1)
  lat(1:nlat) = lat_range(1) + (/(ilat, ilat=0,nlat-1)/) * (lat_range(2)-lat_range(1))/(nlat - 1)
  wind_speed = 10.0

  do j = 1, nlat
    do i = 1,nlon
      call gler_albedo (glt, lon(i), lat(j), wind_speed, alb(i,j), errstat)
      if (errstat /= 0) then
        write(*,*)'gler_albedo failed:  lon=',lon(i),' lat=',lat(j)
        call tell_set_error (0)
        errstat = 0
        continue
      endif
    enddo
  enddo

  call gler_close (glt)

  call write_output_file (lon, lat, alb, errstat)
  deallocate (lon, lat, alb)

  if (errstat /= 0) then
    write (*,*)'*** Error writing output file'
    stop 1
  endif

  write(*,*)'Success'

contains

  subroutine write_output_file (lon, lat, alb, errstat)
    use netcdf
    use tio_module
    implicit none
    real (kind=4), intent(in), dimension(:) :: lon, lat
    real (kind=4), intent(in), dimension(:,:) :: alb
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj
    type (tiof_dimlist_type) :: dimlist
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: attlist
    character (len=3), dimension(2) :: dimnames
    integer, dimension(2) :: dimids, start, edge
    integer :: nlon, nlat

    if (errstat /= 0) return

    nlon = size(lon)
    nlat = size(lat)

    call tiof_create (obj, 'gler_test.nc', nf90_clobber, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'creating file', errstat)
      return
    endif

    call tiof_dimlist_append (dimlist, 'lon', nlon, errstat)
    call tiof_dimlist_append (dimlist, 'lat', nlat, errstat)
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'creating dimensions', errstat)
      return
    endif

    dimnames(1) = 'lon'
    dimnames(2) = 'lat'
    call tiof_dimlist_lookup (dimlist, dimnames, dimids, errstat)

    call tiof_attlist_append (attlist, errstat, "coordinates", &
                              att_text="lat lon")

    call tiof_varlist_append (varlist, errstat, 'lon', nf90_float, &
                              dimids=dimids(1:1))
    call tiof_varlist_append (varlist, errstat, 'lat', nf90_float, &
                              dimids=dimids(2:2))
    call tiof_varlist_append (varlist, errstat, 'alb', nf90_float, &
                              dimids=dimids, attlist=attlist)
    call tiof_def_vars (obj, varlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'creating variables', errstat)
      return
    endif

    call tiof_attlist_free (attlist)
    call tiof_varlist_free (varlist)
    call tiof_dimlist_free (dimlist)

    start(1:2) = 0
    edge(1) = nlat
    edge(2) = nlon

    call tiof_put1d_r4 (obj, 'lat', start(1:1), edge(1:1), lat, errstat)
    call tiof_put1d_r4 (obj, 'lon', start(2:2), edge(2:2), lon, errstat)
    call tiof_put2d_r4 (obj, 'alb', start, edge, alb, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'writing output', errstat)
      return
    endif

    call tiof_close (obj, errstat)

  end subroutine

end program
