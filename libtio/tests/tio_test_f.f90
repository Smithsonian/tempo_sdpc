program tio_test
  use netcdf
  use tio_module
  implicit none

  type (tiof_object_type) :: obj
  character (len=*), parameter :: filename = "delete_radiance.nc"
  character (len=*), parameter :: groupname = "band_290_490_nm"
  integer :: num_wavelengths, num_xtrack, num_steps
  integer, dimension(2) :: dimids
  integer, dimension(3) :: start, edge
  real (kind=4), allocatable, dimension(:,:,:) :: radiance, rx
  integer :: errstat, i,j,k, iwave_start
  real (kind=4) :: n

  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist
  type (tiof_attlist_type) :: foo_attlist

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

  call tiof_dimlist_lookup (dimlist, ["dim2", "dim1"], dimids, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_dimlist_lookup failed'
    stop 2
  endif
  if (any(dimids < 0)) then
    write (*,*)'*** tiof_dimlist_lookup:  dimids = ',dimids
    stop 2
  endif

  call tiof_varlist_append (varlist, errstat, "scalar_int", nf90_int)
  call tiof_varlist_append (varlist, errstat, "dim1", nf90_int, &
                            dimids=[dimids(2)], &
                            comment="dim1 coordinate variable")
  call tiof_varlist_append (varlist, errstat, "dim2", nf90_int, &
                            dimids=[dimids(1)], &
                            comment="dim2 coordinate variable")
  if (errstat < 0) then
    write (*,*)'*** tiof_varlist_append failed'
    stop 2
  endif

  ! global variables
  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_def_vars failed'
    stop 2
  endif

  call tiof_inq_group (obj, groupname, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_group failed:  group='//groupname
    stop 3
  endif

  call tiof_attlist_append (foo_attlist, errstat, &
                            "i4_attr", att_i4=[1,2,3,4])
  call tiof_attlist_append (foo_attlist, errstat, &
                            "r8_attr", att_r8=[1.234_8, 3.1415_8])
  call tiof_attlist_append (foo_attlist, errstat, "text_attr", &
                            att_text='This is a text attribute'// &
                            ' that can wrap around if necessary')

  call tiof_varlist_append (varlist, errstat, "foo", nf90_float, &
                            dimids = dimids, &
                            shuffle=.true., deflate_level=5, &
                            contiguous=.false., chunksizes=[500,500], &
                            comment = "This is foo", &
                            units = "foo-units", valid_range = [-10.0_8, 10.0_8], &
                            attlist = foo_attlist)
  if (errstat < 0) then
    write (*,*)'*** tiof_varlist_append failed'
    stop 3
  endif

  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_def_vars failed'
    stop 3
  endif

  call tiof_inq_group (obj, groupname, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_group failed:  group='//groupname
    stop 4
  endif

  call tiof_inq_dimlen (obj, tempo_dim_channel, num_wavelengths, errstat);
  call tiof_inq_dimlen (obj, tempo_dim_xtrack, num_xtrack, errstat);
  call tiof_inq_dimlen (obj, tempo_dim_step, num_steps, errstat);
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_dimlen failed:  group='//groupname
    stop 4
  endif

  allocate (radiance(num_wavelengths, num_xtrack, num_steps), &
            rx(num_wavelengths,num_xtrack,num_steps), &
           stat=errstat)
  if (errstat /= 0) then
    write(*,*)'*** allocate failed'
    stop 4
  endif

  ! indices in C order, so first is slowest, last changes quickest
  n = 0.0
  do k=0,num_steps-1
    do j=0,num_xtrack-1
      do i=0,num_wavelengths-1
        rx(i+1,j+1,k+1) = n
        n = n + 1.0
      enddo
    enddo
  enddo

  ! Test reading sub-arrays
  ! read test
  start = [2, 0, 0]
  edge  = [2, num_xtrack, num_wavelengths]
  radiance(:,:,:) = 0.0
  call tiof_get3d_r4 (obj, tempo_var_radiance, start, edge, radiance, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed '
    stop 5
  endif
  call compare_arrays

  ! read test
  start = [0, 0, 0]
  edge  = [num_steps, num_xtrack/2, num_wavelengths]
  radiance(:,:,:) = 0.0
  call tiof_get3d_r4 (obj, tempo_var_radiance, start, edge, &
                      radiance(:,1:num_xtrack/2,:), errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed '
    stop 5
  endif
  call compare_arrays

  ! read test
  start = [0, 0, num_wavelengths/2]
  edge  = [num_steps, num_xtrack, num_wavelengths/2]
  radiance(:,:,:) = 0.0
  call tiof_get3d_r4 (obj, tempo_var_radiance, start, edge, &
                      radiance(1:num_wavelengths/2,:,:), errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed '
    stop 5
  endif
  call compare_arrays

  ! Test writing sub-arrays

  ! write/read test
  start = [0, 0, 0]
  edge  = [num_steps, num_xtrack/2, num_wavelengths]
  rx = -1 * rx;
  call tiof_put3d_r4 (obj, tempo_var_radiance, start, edge, &
                      rx(:,1:num_xtrack/2,:), errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_put3d_r4 failed '
    stop 6
  endif
  radiance(:,:,:) = 0.0
  call tiof_get3d_r4 (obj, tempo_var_radiance, start, edge, &
                      radiance(:,1:num_xtrack/2,:), errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed '
    stop 6
  endif
  call compare_arrays

  ! write/read test
  iwave_start = 2   ! starting wavelength index to write out (fortran, 1-based)
  start = [0, 0, iwave_start-1]
  edge  = [num_steps, num_xtrack, num_wavelengths-iwave_start+1]
  rx = -1 * rx;
  call tiof_put3d_r4 (obj, tempo_var_radiance, start, edge, &
                      rx(iwave_start:num_wavelengths,:,:), errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_put3d_r4 failed '
    stop 6
  endif
  radiance(:,:,:) = 0.0
  call tiof_get3d_r4 (obj, tempo_var_radiance, start, edge, &
                      radiance(1:num_wavelengths-iwave_start+1,:,:), errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_get3d_r4 failed '
    stop 6
  endif
  call compare_arrays

  deallocate (radiance,rx)

  call tiof_close (obj, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_close failed'
    stop 4
  endif

contains

  subroutine write_arrays ()
    integer :: i,j,k=1
    do k=0,edge(1)-1
      do j=0,edge(2)-1
        do i=0,edge(3)-1
          write(*,'(i3,i3,i3,a,f6.1,a,f6.1)')i,j,k, &
            ": rx=", rx(1+start(3)+i,1+start(2)+j,1+start(1)+k), &
            " r=", radiance(1+i,1+j,1+k)
        enddo
      enddo
    enddo
    write(*,*)'radiance:'
    do k=1,num_steps
      do j=1,num_xtrack
        do i=1,num_wavelengths
          write(*,'(i3,i3,i3,2x,f6.1)')i,j,k,radiance(i,j,k)
        enddo
      enddo
    enddo
  end subroutine write_arrays

  subroutine compare_arrays ()
    implicit none
    if (any (rx(1+start(3):start(3)+edge(3), &
                1+start(2):start(2)+edge(2), &
                1+start(1):start(1)+edge(1)) &
             /= radiance(1:edge(3),1:edge(2),1:edge(1)))) then
      write(*,*)'*** unexpected input radiance values'
      call write_arrays
      stop 5
    endif
  end subroutine compare_arrays

end program tio_test
