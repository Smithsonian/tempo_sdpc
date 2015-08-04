!> \example example1.f90
!  This is a simple of example of how to use the
!  libtio Fortran interface

program example1
  use netcdf
  use tio_module
  implicit none

  integer, parameter :: n=3, m=5
  integer, dimension(2) :: dimids
  real (kind=8), dimension(n,m) :: twod_out
  real (kind=4), dimension(n,m) :: twod_in

  character (len=*), parameter :: file = "example1.nc"
  type (tiof_file_type) :: obj
  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist
  type (tiof_attlist_type) :: global_attlist, var_attlist
  integer :: errstat = 0

  ! Create an empty HDF5 file:
  call tiof_create (obj, file, nf90_clobber, errstat)
  if (errstat < 0) then
    write(*,*)'*** error creating file='//trim(file)
    stop
  endif

  ! Create a dimlist structure containing a list of dimensions
  ! to be used in the file:
  call tiof_dimlist_append (dimlist, "n", n, errstat)
  call tiof_dimlist_append (dimlist, "m", m, errstat)

  ! Write out the dimension list:
  ! The dimlist structure will save the netcdf id number
  ! assigned to each dimension.
  call tiof_def_dims (obj, dimlist, errstat)
  if (errstat < 0) then
    write(*,*)'*** error defining dimensions'
    stop
  endif

  ! To define a multidimensional variable, we'll need an array of
  ! the associated dimension id numbers (dimids).
  ! Look up the dimid values in the dimlist structure:
  call tiof_dimlist_lookup (dimlist, ["n", "m"], dimids, errstat)

  ! Create a list of global attributes:
  call tiof_attlist_append (global_attlist, errstat, &
                            "informative_text", &
                            att_text = "Standardized attribute"// &
                            " names are preferred.")
  call tiof_attlist_append (global_attlist, errstat, &
                            "global_integer_array", &
                            att_i4 = [2,4,6,8])

  ! Write out the global attribute list:
  call tiof_def_atts (obj, global_attlist, nf90_global, errstat)

  ! Create a list of attributes to attach to a variable:
  call tiof_attlist_append (var_attlist, errstat, &
                            "e", att_r8 = [exp(1.0_8)])

  ! Create a list of variables, with attribute lists if necessary:
  call tiof_varlist_append (varlist, errstat, &
                            "pi", nf90_double, &
                            comment = "A scalar double");
  call tiof_varlist_append (varlist, errstat, &
                            "int1d", nf90_int, &
                            dimids = [dimids(1)], &
                            comment = "An integer array");
  call tiof_varlist_append (varlist, errstat, &
                            "twod", nf90_float, &
                            dimids = dimids, &
                            comment = "A 2D array of floats", &
                            units = "m/s", &
                            valid_range = [-10.0_8, 500.0_8], &
                            fillvalue = -1.0e30_8, &
                            shuffle = .true., &
                            deflate_level = 5, &
                            attlist = var_attlist)

  ! Write out the variable list to the file:
  ! The varlist structure will save the netcdf id
  ! assigned to each variable.
  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    write(*,*)'*** error defining variables'
    stop
  endif

  ! Write data to the file:
  ! Note: 1) in this example, the internal and external types differ,
  !          so type conversion will take place when the values are
  !          written out.
  !       2) the datablock start and size arrays refer to the
  !          array indices in C order, and not in Fortran order.
  !       2) the datablock size array could have been [-1,-1]
  !          instead of [m,n] (a negative size means "use the
  !          full size defined in the file)
  call init_array (twod_out)
  call tiof_put_r8 (obj, "pi", 4.0_8 * atan(1.0_8), errstat)
  call tiof_put1d_i4 (obj, "int1d", [0], [n], [3,6,9], errstat)
  call tiof_put2d_r8 (obj, "twod", [0,0], [m,n], twod_out, errstat)
  call tiof_close (obj, errstat)
  if (errstat < 0) then
    write(*,*)'*** error writing hdf5 file'
  endif

  ! Open the file read-only
  call tiof_open (file, obj, nf90_nowrite, errstat)
  if (errstat < 0) then
    write(*,*)'*** error opening file='//trim(file)
    stop
  endif

  ! Read a subset of the 2D array, and verify that the
  ! input sub-array has the correct values:
  ! (note the correspondence between array elements
  ! of twod_out and twod_in)
  twod_in(:,:) = 0.0_4
  call tiof_get2d_r4 (obj, "twod", [1,0], [3, -1], twod_in, errstat)
  if (any (real(twod_out(:,2:4), kind=4) &
              /= twod_in(:,1:3))) then
    call write_arrays
    stop
  endif

  ! Read a different subset of the 2D array, and verify that the
  ! input sub-array has the correct values:
  ! (note the correspondence between array elements
  ! of twod_out and twod_in)
  twod_in(:,:) = 0.0_4
  call tiof_get2d_r4 (obj, "twod", [0,1], [m,2], twod_in, errstat)
  if (any (real(twod_out(2:3,1:m), kind=4) &
           /= reshape(twod_in, [2,m]))) then
    call write_arrays
    stop
  endif

  ! Close the file
  call tiof_close (obj, errstat)
  if (errstat < 0) then
    write(*,*)'*** error reading from file='//trim(file)
  endif

  write(*,*)'ok'

contains

  subroutine init_array (a)
    implicit none
    real (kind=8), dimension(:,:), intent(inout) :: a
    integer :: i
    a(:,:) = reshape ([(i*1.0_8, i=1,size(a))], shape(a))
  end subroutine init_array

  subroutine write_arrays
    implicit none
    integer :: i,j
    write(*,*)'*** error: array value mismatch'
    write(*,*)'twod_out='
    write(*,'(5(1x,f4.1))')((twod_out(i,j),j=1,m),i=1,n)
    write(*,*)'twod_in ='
    write(*,'(5(1x,f4.1))')((twod_in (i,j),j=1,m),i=1,n)
  end subroutine write_arrays

end program
