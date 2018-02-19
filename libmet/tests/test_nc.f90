program test_nc

  use netcdf
  use tio_module
  use tell_module
  use met_module

  implicit none

  integer(kind=4), parameter :: nx=1, ny=1, nlev=72
  real (kind=4), parameter :: target_lon=-1.0, target_lat=1.0
  real (kind=4), parameter :: correct_psurf=101100.0
  real (kind=4), parameter :: correct_ptrop=11800.0
  real (kind=4), dimension(nlev), parameter :: correct_tprof= &
       (/188.1, &
       188.2, 188.3, 188.4, 188.5, 188.6, 188.7, 188.8, &
       188.9, 189.0, 189.1, 189.2, 189.3, 189.4, 189.5, &
       189.6, 189.7, 189.8, 189.9, 190.0, 190.1, 190.2, &
       190.3, 190.4, 190.5, 190.6, 190.7, 190.8, 190.9, &
       191.0, 191.1, 191.2, 191.3, 191.4, 191.5, 191.6, &
       191.7, 191.8, 191.9, 192.0, 192.1, 192.2, 192.3, &
       192.4, 192.5, 192.6, 192.7, 192.8, 192.9, 193.0, &
       193.1, 193.2, 193.3, 193.4, 193.5, 193.6, 193.7, &
       193.8, 193.9, 194.0, 194.1, 194.2, 194.3, 194.4, &
       194.5, 194.6, 194.7, 194.8, 194.9, 195.0, 195.1, &
       195.2 /)
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
