program test_nc

  use netcdf
  use tio_module
  use tell_module
  use met_module

  implicit none

  integer(kind=4), parameter :: nx=1, ny=1, nlev=72
  real (kind=4), parameter :: target_lon=-1.0d0, target_lat=1.0d0
  real (kind=4), parameter :: correct_psurf=101100.0d0
  real (kind=4), parameter :: correct_ptrop=11800.0d0
  real (kind=4), dimension(nlev), parameter :: correct_tprof= &
       (/188.1d0, &
       188.2d0, 188.3d0, 188.4d0, 188.5d0, 188.6d0, 188.7d0, 188.8d0, &
       188.9d0, 189.0d0, 189.1d0, 189.2d0, 189.3d0, 189.4d0, 189.5d0, &
       189.6d0, 189.7d0, 189.8d0, 189.9d0, 190.0d0, 190.1d0, 190.2d0, &
       190.3d0, 190.4d0, 190.5d0, 190.6d0, 190.7d0, 190.8d0, 190.9d0, &
       191.0d0, 191.1d0, 191.2d0, 191.3d0, 191.4d0, 191.5d0, 191.6d0, &
       191.7d0, 191.8d0, 191.9d0, 192.0d0, 192.1d0, 192.2d0, 192.3d0, &
       192.4d0, 192.5d0, 192.6d0, 192.7d0, 192.8d0, 192.9d0, 193.0d0, &
       193.1d0, 193.2d0, 193.3d0, 193.4d0, 193.5d0, 193.6d0, 193.7d0, &
       193.8d0, 193.9d0, 194.0d0, 194.1d0, 194.2d0, 194.3d0, 194.4d0, &
       194.5d0, 194.6d0, 194.7d0, 194.8d0, 194.9d0, 195.0d0, 195.1d0, &
       195.2d0 /)
  real (kind=4), dimension(nx,ny) :: lat_in, lon_in, ps_in, pt_in
  real (kind=4), dimension(nx,ny,nlev) :: tprof_in

  character (len=11), parameter :: filename="test_met.nc"
  real (kind=4), parameter :: tolerance=1e-8

  integer(kind=4), dimension(2) :: dimid_2d
  integer(kind=4), dimension(3) :: dimid_3d

  real (kind=4) :: psurf, ptrop
  real (kind=4), dimension(nlev) :: tprof

  integer (kind=4) :: errstat

  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist
  type (tiof_file_type) :: obj

  errstat = 0

  !Create a netCDF file
  call tiof_dimlist_append (dimlist, "x", nx, errstat)
  call tiof_dimlist_append (dimlist, "y", ny, errstat)
  call tiof_dimlist_append (dimlist, "z", nlev, errstat)
  if (errstat /= 0) then
    write (*,*)'*** test_nc: tiof_dimlist_append failed'
    stop 1
  endif

  call tiof_create (obj, filename, nf90_clobber, errstat)
  if (errstat /= 0) then
    write(*,*)'*** test_nc: tiof_create failed'
    stop 1
  endif

  call tiof_def_dims (obj, dimlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** test_nc: tiof_def_dims failed'
    stop 1
  endif

  call tiof_dimlist_lookup (dimlist, ["x", "y"], dimid_2d, errstat)
  call tiof_dimlist_lookup (dimlist, ["x", "y", "z"], dimid_3d, errstat)
  if (errstat < 0) then
    write (*,*)'*** test_nc: tiof_dimlist_lookup failed'
    stop 1
  endif

  call tiof_varlist_append (varlist, errstat, "x", nf90_int, &
       dimids=[dimid_2d(2)])
  call tiof_varlist_append (varlist, errstat, "y", nf90_int, &
       dimids=[dimid_2d(1)])
  call tiof_varlist_append (varlist, errstat, "z", nf90_int, &
       dimids=[dimid_3d(1)])
  call tiof_def_vars (obj, varlist, errstat)
  call tiof_varlist_free (varlist)
  if (errstat < 0) then
    write (*,*)'*** test_nc: failed to define coordinate vars'
    stop 1
  endif

  call tiof_varlist_append (varlist, errstat, "lon", nf90_float, &
       dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, "lat", nf90_float, &
       dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, "PS", nf90_float, &
       dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, "TROPPB", nf90_float, &
       dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, "T", nf90_float, &
       dimids=dimid_3d)
  call tiof_def_vars (obj, varlist, errstat)
  call tiof_varlist_free (varlist)
  if (errstat < 0) then
    write (*,*)'*** test_nc: failed to define data vars'
    stop 1
  endif

  lon_in = target_lon
  lat_in = target_lat
  ps_in(nx,ny) = correct_psurf
  pt_in(nx,ny) = correct_ptrop
  tprof_in(nx,ny,1:nlev) = correct_tprof
  call tiof_put2d_r4 (obj, "lon", [0,0], [1,1], lon_in, errstat)
  call tiof_put2d_r4 (obj, "lat", [0,0], [1,1], lat_in, errstat)
  call tiof_put2d_r4 (obj, "PS", [0,0], [1,1], ps_in, errstat)
  call tiof_put2d_r4 (obj, "TROPPB", [0,0], [1,1], pt_in, errstat)
  call tiof_put3d_r4 (obj, "T", [0,0,0], [nlev,1,1], tprof_in, errstat)
  if (errstat < 0) then
    write (*,*)'*** test_nc: failed to write data vars'
    stop 1
  endif

  call tiof_close(obj, errstat)
  if (errstat < 0) then
    write (*,*)'*** test_nc: failed to close test file'
    stop 1
  endif

  call tiof_dimlist_free (dimlist)


  ! Read and check content
  call read_met_data (filename, target_lat, target_lon, &
       ptrop, psurf, tprof, errstat)
  if (errstat /= 0) then
    print *, "*** test_nc: failed to read netCDF file"
    stop 1
  endif

  if (ptrop-correct_ptrop .gt. tolerance) then
    print *, "*** test_nc: failed: troposphere pressure incorrect"
    stop 1
  endif

  if (psurf-correct_psurf .gt. tolerance) then
    print *, "*** test_nc: failed: surface pressure incorrect"
    stop 1
  endif

  if (any(abs(tprof-correct_tprof) .gt. tolerance)) then
    print *, "*** test_nc: failed: temperature profile incorrect"
    stop 1
  endif


end program test_nc
