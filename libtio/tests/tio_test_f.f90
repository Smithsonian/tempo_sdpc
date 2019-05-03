program tio_test
  use netcdf
  use tio_module
  use iso_c_binding, only : c_null_char
  implicit none
  integer, parameter :: &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4), &
    r4 = kind(1.0), &
    r8 = selected_real_kind (2*precision(1.0_r4))

  integer, parameter :: dim_strlen_size = 32, dim_name_size = 20
  integer, parameter :: dim1_size = 5, dim2_size = 7
  integer, parameter :: half_dim1_size = int(dim1_size/2.0)
  integer, parameter :: half_dim2_size = int(dim2_size/2.0)
  character (len=dim_strlen_size), parameter :: dim_strlen = "strlen"
  character (len=dim_strlen_size), parameter :: dim_name = "name"
  integer, dimension(2) :: list_of_names_dimids
  character (len=dim_strlen_size), dimension(dim_name_size) :: list_of_names, &
    input_names

  type (tiof_file_type) :: obj
  character (len=*), parameter :: filename = "delete_radiance.nc"
  character (len=*), parameter :: groupname = "band_290_490_nm"
  integer :: num_wavelengths, num_xtrack, num_steps
  integer, dimension(2) :: dimids
  integer, dimension(3) :: start, edge
  real (kind=r4), allocatable, dimension(:,:,:) :: radiance, rx
  real (kind=r4), dimension(dim1_size,dim2_size) :: fv
  real (kind=r4), dimension(dim1_size) :: iv
  integer (kind=i2), dimension(dim1_size) :: iv_in_i2
  integer (kind=i4), dimension(dim1_size) :: iv_in_i4
  integer (kind=i8), dimension(dim1_size) :: iv_in_i8
  real (kind=r4), dimension(dim1_size) :: rv_in_r4
  real (kind=r8), dimension(dim1_size) :: rv_in_r8
  real (kind=r8), parameter :: fill_value_r4 = 1.25e30, fill_value_r8 = 2.55e30, &
    fill_value_i2 = -9999.0, fill_value_i4 = -99999.0, fill_value_i8 = -1.0e12
  real (kind=r4), parameter :: replace_fill_r4=1.0e20, iv_value = 1234.0
  integer (kind=i2), parameter ::  replace_fill_i2=-15000
  integer (kind=i4), parameter ::  replace_fill_i4=-55000
  integer (kind=i8), parameter ::  replace_fill_i8=-95000
  integer :: errstat, i,j,k, iwave_start
  integer :: scalar_int=-12345, scalar_int_read=0
  real (kind=r4) :: n

  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist
  type (tiof_attlist_type) :: fv_attlist

  errstat = 0

  call tiof_dimlist_append (dimlist, "dim1", dim1_size, errstat)
  call tiof_dimlist_append (dimlist, "dim2", dim2_size, errstat)
  call tiof_dimlist_append (dimlist, dim_strlen, dim_strlen_size, errstat)
  call tiof_dimlist_append (dimlist, dim_name, dim_name_size, errstat)
  if (errstat /= 0) then
    write (*,*)'*** tiof_dimlist_append failed'
    stop 1
  endif

  call tiof_open (filename, obj, nf90_write, errstat)
  if (errstat /= 0) then
    write(*,*)'*** tiof_open failed:  file='//filename
    stop 2
  endif

  call tiof_use_file_epoch (obj, errstat)
  if (errstat /= 0) then
    stop 2
  endif

  call test_granule_ident (obj, errstat)
  if (errstat /= 0) then
    stop 2
  endif

  call tiof_def_dims (obj, dimlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_def_dims failed'
    stop 2
  endif

  call tiof_dimlist_lookup (dimlist, ["dim1", "dim2"], dimids, errstat)
  call tiof_dimlist_lookup (dimlist, [dim_strlen, dim_name], &
                            list_of_names_dimids, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_dimlist_lookup failed'
    stop 2
  endif
  if (any(dimids < 0)) then
    write (*,*)'*** tiof_dimlist_lookup:  dimids = ',dimids
    stop 2
  endif

  call tiof_varlist_append (varlist, errstat, "scalar_int", nf90_int)
  call tiof_varlist_append (varlist, errstat, "scalar_uint", nf90_uint)
  call tiof_varlist_append (varlist, errstat, "dim1", nf90_int, &
                            dimids=[dimids(2)], &
                            comment="dim1 coordinate variable")
  call tiof_varlist_append (varlist, errstat, "list_of_names", nf90_char, &
                            dimids=list_of_names_dimids, &
                            comment="A list of names")
  call tiof_varlist_append (varlist, errstat, "list_of_strings", nf90_string, &
                            dimids=[list_of_names_dimids(2)], &
                            comment="A list of names")
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

  call tiof_put_git_commit_hash (obj, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_put_git_commit_hash failed'
    stop 2
  endif

  call tiof_put_i4 (obj, "scalar_int", scalar_int, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_put_i4 failed'
    stop 2
  endif
  call tiof_get_i4 (obj, "scalar_int", scalar_int_read, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_get_i4 failed'
    stop 2
  endif
  if (scalar_int_read /= scalar_int) then
    write (*,*)'*** I/O of scalar variables failed'
    stop 2
  endif

  call tiof_put_ui4 (obj, "scalar_uint", scalar_int, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_put_ui4 failed'
    stop 2
  endif
  call tiof_get_ui4 (obj, "scalar_uint", scalar_int_read, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_get_ui4 failed'
    stop 2
  endif
  if (scalar_int_read /= scalar_int) then
    write (*,*)'*** I/O of scalar variables failed'
    stop 2
  endif

  list_of_names(:)(:) = ' '
  input_names(:)(:) = ' '
  list_of_names(1) = "Fred"
  list_of_names(2) = "Barney"
  list_of_names(3) = "Aunt Bea"
  ! write names as text
  call tiof_put1d_text (obj, "list_of_names", 0, 3, list_of_names, errstat)
  if (errstat /= 0) then
    write(*,*)'*** tiof_put1d_text failed'
    stop 3
  endif
  call tiof_get1d_text (obj, "list_of_names", 0, 3, input_names, errstat)
  if (errstat /= 0) then
    write(*,*)'*** tiof_get1d_text failed'
    stop 3
  endif
  do i=1,3
    if (trim(list_of_names(i)) /= trim(input_names(i))) then
      write(*,*)'*** string I/O mismatch:'
      write(*,*)' wrote:(', trim(list_of_names(i)),') i=',i
      write(*,*)' read:(', trim(input_names(i)),')'
      stop 3
    endif
  enddo
  ! write names as strings
  input_names(:)(:) = ' '
  do i=1,3
    list_of_names(i) = trim(list_of_names(i))//c_null_char
  enddo
  call tiof_put1d_string (obj, "list_of_strings", 0, 3, list_of_names, errstat)
  if (errstat /= 0) then
    write(*,*)'*** tiof_put1d_string failed'
    stop 3
  endif
  call tiof_get1d_string (obj, "list_of_strings", 0, 3, input_names, errstat)
  if (errstat /= 0) then
    write(*,*)'*** tiof_get1d_string failed'
    stop 3
  endif
  do i=1,3
    if (trim(list_of_names(i)) /= trim(input_names(i))) then
      write(*,*)'*** string I/O mismatch:'
      write(*,*)' wrote:(', trim(list_of_names(i)),') i=',i
      write(*,*)' read:(', trim(input_names(i)),')'
      stop 3
    endif
  enddo

  call tiof_inq_group (obj, groupname, errstat)
  if (errstat < 0) then
    write(*,*)'*** tiof_inq_group failed:  group='//groupname
    stop 3
  endif

  call tiof_attlist_append (fv_attlist, errstat, &
                            "i4_attr", att_i4=[1,2,3,4])
  call tiof_attlist_append (fv_attlist, errstat, &
                            "r8_attr", att_r8=[1.234_8, 3.1415_8])
  call tiof_attlist_append (fv_attlist, errstat, "text_attr", &
                            att_text='This is a text attribute'// &
                            ' that can wrap around if necessary')

  call tiof_varlist_append (varlist, errstat, "fv", nf90_float, &
                            dimids = dimids, &
                            fillvalue = fill_value_r4, &
                            shuffle=.true., deflate_level=5, &
                            contiguous=.false., chunksizes=[1,dim1_size], &
                            comment = "This is fv", &
                            long_name = "ffffvvvv", &
                            units = "fv-units", valid_range = [-10.0_8, 10.0_8], &
                            attlist = fv_attlist)
  call tiof_varlist_append (varlist, errstat, "iv_i2", nf90_short, &
                            dimids=[dimids(1)], &
                            fillvalue = fill_value_i2)
  call tiof_varlist_append (varlist, errstat, "iv_i4", nf90_int, &
                            dimids=[dimids(1)], &
                            fillvalue = fill_value_i4)
  call tiof_varlist_append (varlist, errstat, "iv_i8", nf90_int64, &
                            dimids=[dimids(1)], &
                            fillvalue = fill_value_i8)
  call tiof_varlist_append (varlist, errstat, "rv_r4", nf90_float, &
                            dimids=[dimids(1)], &
                            fillvalue = fill_value_r4)
  call tiof_varlist_append (varlist, errstat, "rv_r8", nf90_double, &
                            dimids=[dimids(1)], &
                            fillvalue = fill_value_r8)
  if (errstat < 0) then
    write (*,*)'*** tiof_varlist_append failed'
    stop 3
  endif

  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_def_vars failed'
    stop 3
  endif

  ! First, test integer I/O and replacing fill values:

  iv = iv_value
  call tiof_put1d_i2 (obj, "iv_i2", [0], [half_dim1_size], int(iv,kind=i2), errstat)
  call tiof_put1d_i4 (obj, "iv_i4", [0], [half_dim1_size], int(iv,kind=i4), errstat)
  call tiof_put1d_i8 (obj, "iv_i8", [0], [half_dim1_size], int(iv,kind=i8), errstat)
  if (errstat < 0) then
    write (*,*)'*** failed: writing iv arrays'
    stop 3
  endif

  ! First try reading the same type that was written to the file
  call tiof_get1d_i2 (obj, "iv_i2", [0], [dim1_size], iv_in_i2, errstat, &
                      replace_fill=replace_fill_i2)
  if (errstat < 0) then
    write(*,*)'*** error reading iv_i2'
    stop
  endif
  if (any(iv_in_i2(1:half_dim1_size) /= int(iv_value,kind=i2)) &
      .or. any(iv_in_i2(1+half_dim1_size:) /= int(replace_fill_i2,kind=i2))) then
    write(*,*)'*** unexpected iv_i2 values'
    stop 3
  endif
  call tiof_get1d_i4 (obj, "iv_i4", [0], [dim1_size], iv_in_i4, errstat, &
                      replace_fill=replace_fill_i4)
  if (errstat < 0) then
    write(*,*)'*** error reading iv_i4'
    stop
  endif
  if (any(iv_in_i4(1:half_dim1_size) /= int(iv_value,kind=i4)) &
      .or. any(iv_in_i4(1+half_dim1_size:) /= int(replace_fill_i4,kind=i4))) then
    write(*,*)'*** unexpected iv_i4 values'
    stop 3
  endif
  call tiof_get1d_i8 (obj, "iv_i8", [0], [dim1_size], iv_in_i8, errstat, &
                      replace_fill=replace_fill_i8)
  if (errstat < 0) then
    write(*,*)'*** error reading iv_i8'
    stop
  endif
  if (any(iv_in_i8(1:half_dim1_size) /= int(iv_value,kind=i8)) &
      .or. any(iv_in_i8(1+half_dim1_size:) /= int(replace_fill_i8,kind=i8))) then
    write(*,*)'*** unexpected iv_i8 values'
    stop 3
  endif

  ! First try reading a larger type than was written to the file
  call tiof_get1d_i4 (obj, "iv_i2", [0], [dim1_size], iv_in_i4, errstat, &
                      replace_fill=replace_fill_i4)
  if (errstat < 0) then
    write(*,*)'*** error reading iv_i2'
    stop
  endif
  if (any(iv_in_i4(1:half_dim1_size) /= int(iv_value,kind=i4)) &
      .or. any(iv_in_i4(1+half_dim1_size:) /= int(replace_fill_i4,kind=i4))) then
    write(*,*)'*** unexpected iv_i2 values'
    !write(*,*)iv_in_i4
    stop 3
  endif

  ! Now, test floating point I/O and replacing fill values:
  iv = iv_value
  call tiof_put1d_r4 (obj, "rv_r4", [0], [half_dim1_size], real(iv,kind=r4), errstat)
  call tiof_put1d_r8 (obj, "rv_r8", [0], [half_dim1_size], real(iv,kind=r8), errstat)
  if (errstat < 0) then
    write (*,*)'*** failed: writing iv arrays'
    stop 3
  endif

  ! First try reading the same type that was written to the file
  call tiof_get1d_r4 (obj, "rv_r4", [0], [dim1_size], rv_in_r4, errstat, &
                      replace_fill=real(replace_fill_i4,kind=r4))
  if (errstat < 0) then
    write(*,*)'*** error reading rv_r4'
    stop
  endif
  if (any(rv_in_r4(1:half_dim1_size) /= real(iv_value,kind=r4)) &
      .or. any(rv_in_r4(1+half_dim1_size:) /= real(replace_fill_i4,kind=r4))) then
    write(*,*)'*** unexpected rv_r4 values'
    stop 3
  endif
  call tiof_get1d_r8 (obj, "rv_r8", [0], [dim1_size], rv_in_r8, errstat, &
                      replace_fill=real(replace_fill_i4,kind=r8))
  if (errstat < 0) then
    write(*,*)'*** error reading rv_r8'
    stop
  endif
  if (any(rv_in_r8(1:half_dim1_size) /= real(iv_value,kind=r8)) &
      .or. any(rv_in_r8(1+half_dim1_size:) /= real(replace_fill_i4,kind=r8))) then
    write(*,*)'*** unexpected rv_r8 values'
    stop 3
  endif

  ! First try reading a larger type than was written to the file
  call tiof_get1d_r8 (obj, "rv_r4", [0], [dim1_size], rv_in_r8, errstat, &
                      replace_fill=real(replace_fill_i4,kind=r8))
  if (errstat < 0) then
    write(*,*)'*** error reading rv_r4'
    stop
  endif
  if (any(rv_in_r8(1:half_dim1_size) /= real(iv_value,kind=r8)) &
      .or. any(rv_in_r8(1+half_dim1_size:) /= real(replace_fill_i4,kind=r8))) then
    write(*,*)'*** unexpected rv_r4 values'
    !write(*,*)rv_in_r8
    stop 3
  endif

  fv = 0.0
  fv(1:dim1_size,1:half_dim2_size) = 1.0e3
  !write(*,'(7(1pe9.3,1x))')transpose(fv)
  call tiof_put2d_r4 (obj, "fv", [0,0], [half_dim2_size, dim1_size], fv, errstat)
  if (errstat < 0) then
    write (*,*)'*** tiof_put2d_r4 failed: writing fv'
    stop 4
  endif

  call tiof_get2d_r4 (obj, "fv", &
                      [0,0], [dim2_size, dim1_size], &
                      fv, errstat, replace_fill=replace_fill_r4)
  !write(*,'(7(1pe9.3,1x))')transpose(fv)
  if (errstat < 0) then
    write (*,*)'*** tiof_get2d_r4 failed: reading fv'
    stop 4
  endif
  if (any (fv(1:dim1_size,1:half_dim2_size) /= 1.0e3)) then
    write(*,*)'*** tiof_get2d_r4: non-fill values were modified!'
    !write(*,'(7(1pe9.3,1x))')transpose(fv)
    stop 4
  endif
  if (any (fv(1:dim1_size,1+half_dim2_size:dim2_size) /= replace_fill_r4)) then
    write (*,*)'*** tiof_get2d_r4 failed: fv has unexpected fill values'
    !write(*,'(7(1pe9.3,1x))')transpose(fv)
    stop 4
  endif

  call tiof_inq_group (obj, "/"//groupname, errstat)
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

  call tiof_dimlist_free (dimlist)
  call tiof_varlist_free (varlist)
  call tiof_attlist_free (fv_attlist)

contains

  subroutine test_granule_ident (obj, errstat)
    use iso_c_binding, only : c_null_char
    implicit none
    type (tiof_file_type), intent(in) :: obj
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj_to
    character (len=1024) :: namebuf
    real (kind=r8) :: tempo_time, hour, hour_expected
    integer :: year, month, day

    if (errstat /= 0) return

    namebuf = repeat('X',len(namebuf))
    call tiof_filename_from_granule (obj, "test"//c_null_char, 1, 1, &
                                     namebuf, errstat)
    if (errstat /= 0) then
      write(*,*)'*** Error: generating filename'
      return
    endif

    tempo_time = 5.42937157308202d+08
    hour_expected = 11.8756411672300768d0;
    call tiof_taix_time_to_utc_caldate (tempo_time, year, &
                                         month, day, hour, errstat)
    if (errstat /= 0) return
    if (year /= 2017 .or. month /= 3 .or. day /= 16 &
        .or. abs(hour - hour_expected)*1.e6 > 1) then
      write(*,*)'*** Error: time conversion yielded unexpected time'
      write(*,*)'year=',year,' month=',month,' day=',day
      write(*,*)'hour=',hour,' expected hour=',hour_expected
      errstat = -1;
      return
    endif

    call tiof_create (obj_to, namebuf, nf90_clobber, errstat)
    if (errstat /= 0) then
      write(*,*)'*** Error: creating output file, errstat=',errstat
      return
    endif

    call tiof_copy_granule_ident (obj, obj_to, errstat)
    if (errstat /= 0) then
      write(*,*)'*** Error: copying granule label'
      return
    endif

    if (tiof_same_granule_ident (obj, obj_to) <= 0) then
      write(*,*)'*** Error: verifying granule label'
      return
    endif

    call tiof_label_product (obj_to, "XXX", 2, errstat)
    if (errstat /= 0) then
      write(*,*)'*** Error: labeling granule'
      return
    endif

    call tiof_close (obj_to, errstat)

    call execute_command_line ('/bin/rm '//namebuf)

  end subroutine test_granule_ident

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
