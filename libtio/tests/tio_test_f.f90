program tio_test
  use netcdf
  use tio_module
  implicit none

  type (L1B_Object_Type) :: l1bobj
  character (len=*), parameter :: filename = "delete_radiance.nc"
  character (len=*), parameter :: groupname = "band_290_490_nm"
  integer :: step0, num_steps
  real (kind=4), allocatable, dimension(:,:,:) :: radiance
  integer :: errstat

  errstat = 0

  call tiof_open (filename, l1bobj, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_open failed:  file='//filename
    stop
  endif

  call tiof_inq_group (l1bobj, groupname, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_group failed:  group='//groupname
    stop
  endif

  step0 = 5
  num_steps = 2
  allocate (radiance(l1bobj%num_wavelengths, l1bobj%num_xtrack, num_steps))

  call tiof_get3d_r4 (l1bobj, "radiance", step0, num_steps, radiance, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed'
    stop
  endif

  !write(*,*)radiance
  deallocate (radiance)

  call tiof_close (l1bobj, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_close failed'
    stop
  endif

end program tio_test
