program ftest1
  use iso_c_binding, only : c_null_char
  use tell_module
  use ftest1_module
  implicit none
  integer, parameter :: &
    r4 = kind(1.0), &
    r8 = selected_real_kind (2*precision(1.0_r4))

  character (len=*), parameter :: filename = "ftest1.nc"
  integer, parameter :: max_dims = 6
  integer, dimension(max_dims), parameter :: dimlens = [3,4,5,6,7,8]

  ! this choice of len and dimension simplifies the test
  character (len=dimlens(max_dims)), dimension(dimlens(1)) :: strings

  integer, dimension(max_dims) :: start, edge
  real (kind=r8), dimension(:), allocatable :: values
  integer :: errstat, num_values, j
  real (kind=r8), parameter :: &
    fillvalue_out = 9.2345, fillvalue_in = -1.5432

  errstat = 0

  num_values = product(dimlens)
  allocate (values(num_values), stat=errstat)
  if (errstat /= 0) then
    call tell_error (tell_malloc_error, "allocate failed", errstat)
    stop 1
  endif

  values = real([(j, j=0,num_values-1)], kind=r8)
  !ensure values do not exceed size that can fit into single byte
  values = 100.0_r8*values/maxval(values)
  values(2::3) = fillvalue_out

  strings = ["I"//c_null_char//"      ", "dislike"//c_null_char, &
            "Fortran"//c_null_char]

  call check_create (filename, values, max_dims, dimlens, &
                     fillvalue_out, strings, errstat)
  if (errstat < 0) stop 1

  values(2::3) = fillvalue_in

  start(:) = 0
  edge(:) = dimlens(:)
  call check_read (filename, values, max_dims, start, edge, &
                   fillvalue_in, strings, errstat)
  if (errstat < 0) stop 1

  deallocate (values)

end program
