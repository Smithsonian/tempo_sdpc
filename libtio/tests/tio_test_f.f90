program tio_test
  use netcdf
  use tio_module
  implicit none

  type (tiof_object_type) :: obj
  character (len=*), parameter :: filename = "delete_radiance.nc"
  character (len=*), parameter :: groupname = "band_290_490_nm"
  integer :: step0, num_steps, num_wavelengths, num_xtrack
  integer :: foo_varid
  integer, dimension(2) :: dimids
  real (kind=4), allocatable, dimension(:,:,:) :: radiance
  integer :: errstat

  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist

  errstat = 0

  call tiof_dimlist_append (dimlist, "dim1", 10000, errstat)
  call tiof_dimlist_append (dimlist, "dim2", 20000, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_dimlist_append failed'
    stop 1
  endif

  call tiof_open (filename, obj, nf90_write, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_open failed:  file='//filename
    stop 2
  endif

  call tiof_def_dims (obj, dimlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_def_dims failed'
    stop 2
  endif

  call tiof_dimlist_lookup (dimlist, 2, ["dim2", "dim1"], dimids, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_dimlist_lookup failed'
    stop 2
  endif

  call tiof_varlist_append (varlist, "foo", nf90_float, dimids, errstat, &
                            shuffle=.true., deflate_level=5, &
                            contiguous=.false., chunksizes=[500,500])
  if (errstat < 0) then
    write (*,*)'*** tiof_varlist_append failed'
    stop 2
  endif

  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_def_vars failed'
    stop 2
  endif

  call tiof_varlist_lookup (varlist, "foo", foo_varid, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_varlist_lookup failed'
    stop 2
  endif

  call tiof_inq_group (obj, groupname, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_group failed:  group='//groupname
    stop 2
  endif

  call tiof_inq_dimlen (obj, "spectral_channel", num_wavelengths, errstat);
  call tiof_inq_dimlen (obj, "xtrack", num_xtrack, errstat);
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_dimlen failed:  group='//groupname
    stop 2
  endif

  step0 = 5
  num_steps = 2
  allocate (radiance(num_wavelengths, num_xtrack, num_steps))

  call tiof_get3d_r4 (obj, "radiance", step0, num_steps, radiance, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed'
    stop 3
  endif

  !write(*,*)radiance
  deallocate (radiance)

  call tiof_close (obj, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_close failed'
    stop 4
  endif

end program tio_test
