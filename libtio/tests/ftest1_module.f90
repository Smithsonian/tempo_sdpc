module ftest1_module

  integer, private, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4), &
    r4 = kind(1.0), &
    r8 = selected_real_kind (2*precision(1.0_r4))

contains

subroutine expected_fail (errstat)
  implicit none
  integer, intent(inout) :: errstat

  if (errstat < 0) then
    errstat = 0
    return
  endif
  write (*,*) "*** expected failure did not occur!!"
  stop 1

end subroutine

subroutine create_group (obj, errstat)
  use tell_module
  use tio_module
  use netcdf
  implicit none
  type (tiof_file_type), intent(inout) :: obj
  integer, intent(inout) :: errstat

  character (len=*), parameter :: grpname = "a_group"
  type (tiof_varlist_type) :: varlist
  type (tiof_var_type), pointer :: var_info
  type (tiof_attlist_type) :: attlist
  integer :: grp

  if (errstat < 0) return

  call tiof_def_group (obj, grpname, errstat, grp)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, "tiof_def_group failed", &
                     errstat);
    return
  endif

  call tiof_push_group (obj, grpname, errstat)
  if (errstat < 0) then
    call tell_error (tell_runtime_error, "tiof_push_group failed", &
                     errstat);
    return
  endif

  call tiof_attlist_append (attlist, errstat, "i4_att", att_i4=[111_i4])
  call tiof_attlist_append (attlist, errstat, "i8_att", att_i8=[222_i8])
  call tiof_attlist_append (attlist, errstat, "r4_att", att_r4=[111.0_r4])
  call tiof_attlist_append (attlist, errstat, "r8_att", att_r8=[222.0_r8])
  call tiof_attlist_append (attlist, errstat, "txt_att", &
                            att_text="attribute text")

  call tiof_varlist_append (varlist, errstat, "var", nf90_int, &
                            comment="comment: a group variable", &
                            attlist=attlist)
  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, &
                     "failed defining group variable", errstat);
    return
  endif

  call tiof_put_i4 (obj, "var", 9913, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, &
                     "failed writing group variable", errstat);
    return
  endif

  call tiof_pop_group (obj, errstat)

  call tiof_varlist_lookup (varlist, "var", var_info, errstat)
  if (errstat < 0) then
    call tell_error (tell_runtime_error, &
                     "varlist lookup failed", errstat);
    return
  endif

  call tiof_def_atts (obj, attlist, nf90_global, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, &
                     "failed writing global attributes", errstat);
    return
  endif

  call tiof_varlist_free (varlist)
  call tiof_attlist_free (attlist)

end subroutine create_group

subroutine read_from_group (obj, errstat)
  use tell_module
  use tio_module
  implicit none
  type (tiof_file_type), intent(inout) :: obj
  integer, intent(inout) :: errstat

  integer (kind=i4) :: var

  if (errstat < 0) return

  call tiof_inq_group (obj, "a_group", errstat)
  if (errstat < 0) then
    call tell_error (tell_io_read_error, "failed accessing group", &
                     errstat)
    return
  endif

  call tiof_get_i4 (obj, "var", var, errstat)
  if (var /= 9913) then
    call tell_error (tell_io_read_error, &
                     "read wrong value for i4 group variable", &
                     errstat)
    return
  endif

  call tiof_inq_group (obj, "/", errstat)

end subroutine read_from_group

subroutine test_put_errors (obj, max_dims, errstat)
  use tell_module
  use tio_module
  implicit none
  type (tiof_file_type), intent(inout) :: obj
  integer, intent(in) :: max_dims
  integer, intent(inout) :: errstat

  integer, dimension(max_dims) :: start, edge
  real (kind=r8), dimension(max_dims) :: values
  character (len=1), dimension(max_dims) :: strings
  character (len=*), parameter :: absent = "absent"

  if (errstat < 0) return

  call tell_push_queue ()

  start(:) = 0
  edge(:) = 1
  values(:) = 0.0
  strings(:) = " "

  call tiof_put_i1 (obj, absent, int(values(1),kind=i1), errstat)
  call expected_fail (errstat)
  call tiof_put_i2 (obj, absent, int(values(1),kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put_i4 (obj, absent, int(values(1),kind=i4), errstat)
  call expected_fail (errstat)
  call tiof_put_i8 (obj, absent, int(values(1),kind=i8), errstat)
  call expected_fail (errstat)
  call tiof_put_ui1 (obj, absent, int(values(1),kind=i1), errstat)
  call expected_fail (errstat)
  call tiof_put_ui2 (obj, absent, int(values(1),kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put_ui4 (obj, absent, int(values(1),kind=i4), errstat)
  call expected_fail (errstat)
  call tiof_put_ui8 (obj, absent, int(values(1),kind=i8), errstat)
  call expected_fail (errstat)
  call tiof_put_r4 (obj, absent, real(values(1),kind=r4), errstat)
  call expected_fail (errstat)
  call tiof_put_r8 (obj, absent, real(values(1),kind=r8), errstat)
  call expected_fail (errstat)

  call tiof_put1d_text (obj, absent, start(1), edge(1), &
                        strings, errstat)
  call expected_fail (errstat)
  call tiof_put1d_string (obj, absent, start(1), edge(1), &
                          strings, errstat)
  call expected_fail (errstat)

  ! integer (kind=i1)
  call tiof_put1d_i1 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i1), errstat)
  call expected_fail (errstat)
  call tiof_put2d_i1 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i1), errstat)
  call expected_fail (errstat)
  call tiof_put3d_i1 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i1), errstat)
  call expected_fail (errstat)

  call tiof_put1d_ui1 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i1), errstat)
  call expected_fail (errstat)
  call tiof_put2d_ui1 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i1), errstat)
  call expected_fail (errstat)
  call tiof_put3d_ui1 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i1), errstat)
  call expected_fail (errstat)

  ! integer (kind=i2)
  call tiof_put1d_i2 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put2d_i2 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put3d_i2 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put4d_i2 (obj, absent, start(1:4), edge(1:4), &
                      int(reshape(values, edge(1:4)), kind=i2), errstat)
  call expected_fail (errstat)

  call tiof_put1d_ui2 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put2d_ui2 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i2), errstat)
  call expected_fail (errstat)
  call tiof_put3d_ui2 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i2), errstat)
  call expected_fail (errstat)

  ! integer (kind=i4)
  call tiof_put1d_i4 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i4), errstat)
  call expected_fail (errstat)
  call tiof_put2d_i4 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i4), errstat)
  call expected_fail (errstat)
  call tiof_put3d_i4 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i4), errstat)
  call expected_fail (errstat)

  call tiof_put1d_ui4 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i4), errstat)
  call expected_fail (errstat)
  call tiof_put2d_ui4 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i4), errstat)
  call expected_fail (errstat)
  call tiof_put3d_ui4 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i4), errstat)
  call expected_fail (errstat)

  ! integer (kind=i8)
  call tiof_put1d_i8 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i8), errstat)
  call expected_fail (errstat)
  call tiof_put2d_i8 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i8), errstat)
  call expected_fail (errstat)
  call tiof_put3d_i8 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i8), errstat)
  call expected_fail (errstat)

  call tiof_put1d_ui8 (obj, absent, start(1:1), edge(1:1), &
                      int(values,kind=i8), errstat)
  call expected_fail (errstat)
  call tiof_put2d_ui8 (obj, absent, start(1:2), edge(1:2), &
                      int(reshape(values, edge(1:2)), kind=i8), errstat)
  call expected_fail (errstat)
  call tiof_put3d_ui8 (obj, absent, start(1:3), edge(1:3), &
                      int(reshape(values, edge(1:3)), kind=i8), errstat)
  call expected_fail (errstat)

  ! real (kind=r4)
  call tiof_put1d_r4 (obj, absent, start(1:1), edge(1:1), &
                      real(values,kind=r4), errstat)
  call expected_fail (errstat)
  call tiof_put2d_r4 (obj, absent, start(1:2), edge(1:2), &
                      real(reshape(values, edge(1:2)), kind=r4), errstat)
  call expected_fail (errstat)
  call tiof_put3d_r4 (obj, absent, start(1:3), edge(1:3), &
                      real(reshape(values, edge(1:3)), kind=r4), errstat)
  call expected_fail (errstat)
  call tiof_put4d_r4 (obj, absent, start(1:4), edge(1:4), &
                      real(reshape(values, edge(1:4)), kind=r4), errstat)
  call expected_fail (errstat)
  call tiof_put5d_r4 (obj, absent, start(1:5), edge(1:5), &
                      real(reshape(values, edge(1:5)), kind=r4), errstat)
  call expected_fail (errstat)
  call tiof_put6d_r4 (obj, absent, start(1:6), edge(1:6), &
                      real(reshape(values, edge(1:6)), kind=r4), errstat)
  call expected_fail (errstat)

  ! real (kind=r8)
  call tiof_put1d_r8 (obj, absent, start(1:1), edge(1:1), &
                      real(values,kind=r8), errstat)
  call expected_fail (errstat)
  call tiof_put2d_r8 (obj, absent, start(1:2), edge(1:2), &
                      real(reshape(values, edge(1:2)), kind=r8), errstat)
  call expected_fail (errstat)
  call tiof_put3d_r8 (obj, absent, start(1:3), edge(1:3), &
                      real(reshape(values, edge(1:3)), kind=r8), errstat)
  call expected_fail (errstat)
  call tiof_put4d_r8 (obj, absent, start(1:4), edge(1:4), &
                      real(reshape(values, edge(1:4)), kind=r8), errstat)
  call expected_fail (errstat)
  call tiof_put5d_r8 (obj, absent, start(1:5), edge(1:5), &
                      real(reshape(values, edge(1:5)), kind=r8), errstat)
  call expected_fail (errstat)
  call tiof_put6d_r8 (obj, absent, start(1:6), edge(1:6), &
                      real(reshape(values, edge(1:6)), kind=r8), errstat)
  call expected_fail (errstat)

  call tell_pop_queue (1)

end subroutine test_put_errors

subroutine test_get_errors (obj, max_dims, errstat)
  use tell_module
  use tio_module
  implicit none
  type (tiof_file_type), intent(inout) :: obj
  integer, intent(in) :: max_dims
  integer, intent(inout) :: errstat

  character (len=*), parameter :: absent = "absent"
  integer, dimension(max_dims) :: start, edge
  integer (kind=i1) :: i1s
  integer (kind=i2) :: i2s
  integer (kind=i4) :: i4s
  integer (kind=i8) :: i8s
  real (kind=r4) :: r4s
  real (kind=r8) :: r8s
  character (len=1), dimension(1) :: strings
  integer (kind=i1) :: i1_1d(1), i1_2d(1,1), i1_3d(1,1,1)
  integer (kind=i2) :: i2_1d(1), i2_2d(1,1), i2_3d(1,1,1), i2_4d(1,1,1,1)
  integer (kind=i4) :: i4_1d(1), i4_2d(1,1), i4_3d(1,1,1)
  integer (kind=i8) :: i8_1d(1), i8_2d(1,1), i8_3d(1,1,1)
  real (kind=r4) :: &
    r4_1d(1), r4_2d(1,1), r4_3d(1,1,1), &
    r4_4d(1,1,1,1), r4_5d(1,1,1,1,1), r4_6d(1,1,1,1,1,1)
  real (kind=r8) :: &
    r8_1d(1), r8_2d(1,1), r8_3d(1,1,1), &
    r8_4d(1,1,1,1), r8_5d(1,1,1,1,1), r8_6d(1,1,1,1,1,1)

  if (errstat < 0) return

  call tell_push_queue ()

  start(:) = 0
  edge(:) = 1

  call tiof_get_i1 (obj, absent, i1s, errstat)
  call expected_fail (errstat)
  call tiof_get_i2 (obj, absent, i2s, errstat)
  call expected_fail (errstat)
  call tiof_get_i4 (obj, absent, i4s, errstat)
  call expected_fail (errstat)
  call tiof_get_i8 (obj, absent, i8s, errstat)
  call expected_fail (errstat)
  call tiof_get_ui1 (obj, absent, i1s, errstat)
  call expected_fail (errstat)
  call tiof_get_ui2 (obj, absent, i2s, errstat)
  call expected_fail (errstat)
  call tiof_get_ui4 (obj, absent, i4s, errstat)
  call expected_fail (errstat)
  call tiof_get_ui8 (obj, absent, i8s, errstat)
  call expected_fail (errstat)
  call tiof_get_r4 (obj, absent, r4s, errstat)
  call expected_fail (errstat)
  call tiof_get_r8 (obj, absent, r8s, errstat)
  call expected_fail (errstat)

  call tiof_get1d_text (obj, absent, start(1), edge(1), &
                        strings, errstat)
  call expected_fail (errstat)
  call tiof_get1d_string (obj, absent, start(1), edge(1), &
                          strings, errstat)
  call expected_fail (errstat)

  ! integer (kind=i1)
  call tiof_get1d_i1 (obj, absent, start(1:1), edge(1:1), &
                      i1_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_i1 (obj, absent, start(1:2), edge(1:2), &
                      i1_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_i1 (obj, absent, start(1:3), edge(1:3), &
                      i1_3d, errstat)
  call expected_fail (errstat)

  call tiof_get1d_ui1 (obj, absent, start(1:1), edge(1:1), &
                      i1_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_ui1 (obj, absent, start(1:2), edge(1:2), &
                      i1_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_ui1 (obj, absent, start(1:3), edge(1:3), &
                      i1_3d, errstat)
  call expected_fail (errstat)

  ! integer (kind=i2)
  call tiof_get1d_i2 (obj, absent, start(1:1), edge(1:1), &
                      i2_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_i2 (obj, absent, start(1:2), edge(1:2), &
                      i2_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_i2 (obj, absent, start(1:3), edge(1:3), &
                      i2_3d, errstat)
  call expected_fail (errstat)
  call tiof_get4d_i2 (obj, absent, start(1:4), edge(1:4), &
                      i2_4d, errstat)
  call expected_fail (errstat)

  call tiof_get1d_ui2 (obj, absent, start(1:1), edge(1:1), &
                      i2_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_ui2 (obj, absent, start(1:2), edge(1:2), &
                      i2_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_ui2 (obj, absent, start(1:3), edge(1:3), &
                      i2_3d, errstat)
  call expected_fail (errstat)

  ! integer (kind=i4)
  call tiof_get1d_i4 (obj, absent, start(1:1), edge(1:1), &
                      i4_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_i4 (obj, absent, start(1:2), edge(1:2), &
                      i4_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_i4 (obj, absent, start(1:3), edge(1:3), &
                      i4_3d, errstat)
  call expected_fail (errstat)

  call tiof_get1d_ui4 (obj, absent, start(1:1), edge(1:1), &
                      i4_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_ui4 (obj, absent, start(1:2), edge(1:2), &
                      i4_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_ui4 (obj, absent, start(1:3), edge(1:3), &
                      i4_3d, errstat)
  call expected_fail (errstat)

  ! integer (kind=i8)
  call tiof_get1d_i8 (obj, absent, start(1:1), edge(1:1), &
                      i8_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_i8 (obj, absent, start(1:2), edge(1:2), &
                      i8_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_i8 (obj, absent, start(1:3), edge(1:3), &
                      i8_3d, errstat)
  call expected_fail (errstat)

  call tiof_get1d_ui8 (obj, absent, start(1:1), edge(1:1), &
                      i8_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_ui8 (obj, absent, start(1:2), edge(1:2), &
                      i8_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_ui8 (obj, absent, start(1:3), edge(1:3), &
                      i8_3d, errstat)
  call expected_fail (errstat)

  ! real (kind=r4)
  call tiof_get1d_r4 (obj, absent, start(1:1), edge(1:1), &
                      r4_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_r4 (obj, absent, start(1:2), edge(1:2), &
                      r4_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_r4 (obj, absent, start(1:3), edge(1:3), &
                      r4_3d, errstat)
  call expected_fail (errstat)
  call tiof_get4d_r4 (obj, absent, start(1:4), edge(1:4), &
                      r4_4d, errstat)
  call expected_fail (errstat)
  call tiof_get5d_r4 (obj, absent, start(1:5), edge(1:5), &
                      r4_5d, errstat)
  call expected_fail (errstat)
  call tiof_get6d_r4 (obj, absent, start(1:6), edge(1:6), &
                      r4_6d, errstat)
  call expected_fail (errstat)

  ! real (kind=r8)
  call tiof_get1d_r8 (obj, absent, start(1:1), edge(1:1), &
                      r8_1d, errstat)
  call expected_fail (errstat)
  call tiof_get2d_r8 (obj, absent, start(1:2), edge(1:2), &
                      r8_2d, errstat)
  call expected_fail (errstat)
  call tiof_get3d_r8 (obj, absent, start(1:3), edge(1:3), &
                      r8_3d, errstat)
  call expected_fail (errstat)
  call tiof_get4d_r8 (obj, absent, start(1:4), edge(1:4), &
                      r8_4d, errstat)
  call expected_fail (errstat)
  call tiof_get5d_r8 (obj, absent, start(1:5), edge(1:5), &
                      r8_5d, errstat)
  call expected_fail (errstat)
  call tiof_get6d_r8 (obj, absent, start(1:6), edge(1:6), &
                      r8_6d, errstat)
  call expected_fail (errstat)

  call tell_pop_queue (1)

end subroutine test_get_errors

subroutine test_open_errors (errstat)
  use tell_module
  use tio_module
  use netcdf
  implicit none
  integer, intent(inout) :: errstat

  type (tiof_file_type) :: obj

  call tell_push_queue ()

  ! FIXME: As a test, I tried creating /dev/full (the directory,
  !         not a file _in_ the directory).
  !        tiof_create succeeded (!) but I got a segv when the
  !        program exited.  Since netcdf didnt catch this,
  !        the library probably should.

  ! FIXME: netcdf doesn't catch this either:
  !call tiof_create (obj, "/dev/full/xxx", nf90_clobber, errstat)
  !call expected_fail (errstat)

  call tiof_open ("/tmp/absent", obj, nf90_nowrite, errstat)
  call expected_fail (errstat)

  call tell_pop_queue (1)

end subroutine test_open_errors

subroutine test_usage_errors (obj, errstat)
  use tio_module
  implicit none
  type (tiof_file_type), intent(inout) :: obj
  integer, intent(inout) :: errstat

  integer :: dimlen

  if (errstat < 0) return

  call tell_push_queue ()

  ! this doesn't cause an error, but it does exercise
  ! some special case code
  call tiof_pop_group (obj, errstat)

  call tiof_inq_dimlen (obj, "absent", dimlen, errstat)
  call expected_fail (errstat)

  call tell_pop_queue (1)

end subroutine test_usage_errors

subroutine check_create (filename, values, max_dims, dimlens, &
                         fillvalue, strings, errstat)
  use tell_module
  use tio_module
  use netcdf
  implicit none

  character (len=*), intent(in) :: filename
  real (kind=r8), dimension(:), intent(in) :: values
  real (kind=r8), intent(in) :: fillvalue
  character (len=*), dimension(:), intent(in) :: strings
  integer, intent(in) :: max_dims
  integer, dimension(:), intent(in) :: dimlens
  integer, intent(inout) :: errstat

  type (tiof_file_type) :: obj
  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist

  integer, dimension(max_dims) :: dimids, dimsizes
  character (len=2), dimension(max_dims) :: dimnames
  character (len=2), dimension(max_dims) :: anames

  integer, parameter :: max_types = 10
  character (len=2), dimension(max_types), parameter :: type_prefixes = &
    ["sb","ss","si","sj", &  !   signed byte, short, int, int64
     "ub","us","ui","uj", &  ! unsigned byte, short, int, int64
     "rf","rd"]              ! float, double
  integer :: num_dims

  integer, dimension(max_dims) :: start, edge
  integer, dimension(max_types) :: types
  character (len=4) :: typed_aname
  character (len=3) :: typed_sname
  integer :: i, j

  if (max_dims>9) then
    call tell_error (tell_internal_error, &
                     "*** name strings are sized for <= 9 dimensions", &
                     errstat)
    return
  endif

  do i=1,max_dims
    write(dimnames(i), "('n',i1)")i
    write(anames(i),"('A',i1)")i
  enddo

  do i=1,max_dims
    call tiof_dimlist_append (dimlist, dimnames(i), dimlens(i), errstat)
    if (errstat < 0) then
      call tell_error (tell_runtime_error, "tiof_dimlist_append failed", errstat);
      return
    endif
  enddo

  call test_open_errors (errstat)

  call tiof_create (obj, filename, nf90_clobber, errstat)
  call tiof_put_git_commit_hash (obj, errstat)
  call tiof_put_git_commit_hash (obj, errstat, "")  ! exercise option

  call test_usage_errors (obj, errstat)

  call tiof_def_dims (obj, dimlist, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, "file init failed", errstat);
    return
  endif

  ! don't need sizes, but exercise the code anyway
  call tiof_dimlist_lookup (dimlist, dimnames, dimids, errstat, &
                            dimsizes)
  if (errstat < 0) then
    call tell_error (tell_runtime_error, "dimlist lookup failed", errstat);
    return
  endif
  if (any(dimsizes /= dimlens)) then
    call tell_error (tell_runtime_error, &
                     "dimlist lookup returned incorrect dimsizes", &
                     errstat)
    return
  endif

  call tiof_varlist_append (varlist, errstat, "txt1d", nf90_char, &
                            dimids=[dimids(max_dims), dimids(1)], &
                            comment="comment:  txt1d")
  call tiof_varlist_append (varlist, errstat, "str1d", nf90_string, &
                            dimids=dimids(1:1), &
                            comment="comment:  str1d")

  types = [ nf90_byte,  nf90_short,  nf90_int,  nf90_int64, &
           nf90_ubyte, nf90_ushort, nf90_uint, nf90_uint64, &
           nf90_float, nf90_double]

  do j=1,max_types
    if (types(j) == nf90_float .or. types(j) == nf90_double) then
      num_dims = max_dims
    else
      num_dims = 3
    endif

    typed_sname = type_prefixes(j)//'S'
    call tiof_varlist_append (varlist, errstat, typed_sname, types(j), &
                              fillvalue=fillvalue, &
                              comment="comment: "//typed_sname)

    do i=1,num_dims
      typed_aname = type_prefixes(j)//anames(i)
      call tiof_varlist_append (varlist, errstat, typed_aname, types(j), &
                                dimids=dimids(1:i), &
                                fillvalue=fillvalue, &
                                comment="comment: "//typed_aname )
      if (errstat < 0) then
        call tell_error (tell_runtime_error, "varlist append failed", errstat);
        return
      endif
    enddo
  enddo

  ! rather than support all 4D arrays, just add this one explicitly
  call tiof_varlist_append (varlist, errstat, "ssA4", nf90_short, dimids=dimids(1:4), &
                            fillvalue=fillvalue, comment="comment: ssA4")
  if (errstat < 0) then
    call tell_error (tell_runtime_error, "varlist append failed", errstat);
    return
  endif

  ! exercise variable attributes
  call tiof_varlist_append (varlist, errstat, "many_attributes", &
                            nf90_int, &
                            dimids=[dimids(max_dims)], &
                            fillvalue=fillvalue, &
                            contiguous=.false., &
                            chunksizes=[2], &
                            shuffle=.false., &
                            deflate_level=1, &
                            units="whatever", &
                            valid_range=[0.0_r8,999.0_r8], &
                            no_fill=1, &
                            comment="comment: xxx")

  call tiof_def_vars (obj, varlist, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, "def_vars failed", errstat);
    return
  endif

  call tiof_put_i1 (obj, "sbS", int(values(1),kind=i1), errstat)
  call tiof_put_i2 (obj, "ssS", int(values(1),kind=i2), errstat)
  call tiof_put_i4 (obj, "siS", int(values(1),kind=i4), errstat)
  call tiof_put_i8 (obj, "sjS", int(values(1),kind=i8), errstat)
  call tiof_put_ui1 (obj, "ubS", int(values(1),kind=i1), errstat)
  call tiof_put_ui2 (obj, "usS", int(values(1),kind=i2), errstat)
  call tiof_put_ui4 (obj, "uiS", int(values(1),kind=i4), errstat)
  call tiof_put_ui8 (obj, "ujS", int(values(1),kind=i8), errstat)
  call tiof_put_r4 (obj, "rfS", real(values(1),kind=r4), errstat)
  call tiof_put_r8 (obj, "rdS", real(values(1),kind=r8), errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, "scalar write failed", errstat);
    return
  endif

  start(:) = 0
  edge(:) = dimlens(:)

  ! text array
  call tiof_put1d_text (obj, "txt1d", start(1), edge(1), &
                        strings, errstat)
  ! string array
  call tiof_put1d_string (obj, "str1d", start(1), edge(1), &
                          strings, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_write_error, "string array write failed", errstat);
    return
  endif

  ! integer (kind=i1)
  call tiof_put1d_i1 (obj, "sbA1", start(1:1), edge(1:1), &
                      int(values,kind=i1), errstat)
  call tiof_put2d_i1 (obj, "sbA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i1), errstat)
  call tiof_put3d_i1 (obj, "sbA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i1), errstat)

  call tiof_put1d_ui1 (obj, "ubA1", start(1:1), edge(1:1), &
                      int(values,kind=i1), errstat)
  call tiof_put2d_ui1 (obj, "ubA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i1), errstat)
  call tiof_put3d_ui1 (obj, "ubA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i1), errstat)

  if (errstat < 0) then
    call tell_error (tell_io_write_error, "i1 write failed", errstat);
    return
  endif

  ! integer (kind=i2)
  call tiof_put1d_i2 (obj, "ssA1", start(1:1), edge(1:1), &
                      int(values,kind=i2), errstat)
  call tiof_put2d_i2 (obj, "ssA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i2), errstat)
  call tiof_put3d_i2 (obj, "ssA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i2), errstat)
  call tiof_put4d_i2 (obj, "ssA4", start(4:1:-1), edge(4:1:-1), &
                      int(reshape(values(1:product(edge(1:4))), edge(1:4)),kind=i2), errstat)

  call tiof_put1d_ui2 (obj, "usA1", start(1:1), edge(1:1), &
                      int(values,kind=i2), errstat)
  call tiof_put2d_ui2 (obj, "usA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i2), errstat)
  call tiof_put3d_ui2 (obj, "usA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i2), errstat)

  if (errstat < 0) then
    call tell_error (tell_io_write_error, "i2 write failed", errstat);
    return
  endif

  ! integer (kind=i4)
  call tiof_put1d_i4 (obj, "siA1", start(1:1), edge(1:1), &
                      int(values,kind=i4), errstat)
  call tiof_put2d_i4 (obj, "siA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i4), errstat)
  call tiof_put3d_i4 (obj, "siA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i4), errstat)

  call tiof_put1d_ui4 (obj, "uiA1", start(1:1), edge(1:1), &
                      int(values,kind=i4), errstat)
  call tiof_put2d_ui4 (obj, "uiA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i4), errstat)
  call tiof_put3d_ui4 (obj, "uiA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i4), errstat)

  if (errstat < 0) then
    call tell_error (tell_io_write_error, "i4 write failed", errstat);
    return
  endif

  ! integer (kind=i8)
  call tiof_put1d_i8 (obj, "sjA1", start(1:1), edge(1:1), &
                      int(values,kind=i8), errstat)
  call tiof_put2d_i8 (obj, "sjA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i8), errstat)
  call tiof_put3d_i8 (obj, "sjA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i8), errstat)

  call tiof_put1d_ui8 (obj, "ujA1", start(1:1), edge(1:1), &
                      int(values,kind=i8), errstat)
  call tiof_put2d_ui8 (obj, "ujA2", start(2:1:-1), edge(2:1:-1), &
                      int(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=i8), errstat)
  call tiof_put3d_ui8 (obj, "ujA3", start(3:1:-1), edge(3:1:-1), &
                      int(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=i8), errstat)

  if (errstat < 0) then
    call tell_error (tell_io_write_error, "i8 write failed", errstat);
    return
  endif

  ! real (kind=r4)
  call tiof_put1d_r4 (obj, "rfA1", start(1:1), edge(1:1), &
                      real(values,kind=r4), errstat)
  call tiof_put2d_r4 (obj, "rfA2", start(2:1:-1), edge(2:1:-1), &
                      real(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=r4), errstat)
  call tiof_put3d_r4 (obj, "rfA3", start(3:1:-1), edge(3:1:-1), &
                      real(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=r4), errstat)
  call tiof_put4d_r4 (obj, "rfA4", start(4:1:-1), edge(4:1:-1), &
                      real(reshape(values(1:product(edge(1:4))), edge(1:4)),kind=r4), errstat)
  call tiof_put5d_r4 (obj, "rfA5", start(5:1:-1), edge(5:1:-1), &
                      real(reshape(values(1:product(edge(1:5))), edge(1:5)),kind=r4), errstat)
  call tiof_put6d_r4 (obj, "rfA6", start(6:1:-1), edge(6:1:-1), &
                      real(reshape(values(1:product(edge(1:6))), edge(1:6)),kind=r4), errstat)

  if (errstat < 0) then
    call tell_error (tell_io_write_error, "r4 write failed", errstat);
    return
  endif

  ! real (kind=r8)
  call tiof_put1d_r8 (obj, "rdA1", start(1:1), edge(1:1), &
                      real(values,kind=r8), errstat)
  call tiof_put2d_r8 (obj, "rdA2", start(2:1:-1), edge(2:1:-1), &
                      real(reshape(values(1:product(edge(1:2))), edge(1:2)),kind=r8), errstat)
  call tiof_put3d_r8 (obj, "rdA3", start(3:1:-1), edge(3:1:-1), &
                      real(reshape(values(1:product(edge(1:3))), edge(1:3)),kind=r8), errstat)
  call tiof_put4d_r8 (obj, "rdA4", start(4:1:-1), edge(4:1:-1), &
                      real(reshape(values(1:product(edge(1:4))), edge(1:4)),kind=r8), errstat)
  call tiof_put5d_r8 (obj, "rdA5", start(5:1:-1), edge(5:1:-1), &
                      real(reshape(values(1:product(edge(1:5))), edge(1:5)),kind=r8), errstat)
  call tiof_put6d_r8 (obj, "rdA6", start(6:1:-1), edge(6:1:-1), &
                      real(reshape(values(1:product(edge(1:6))), edge(1:6)),kind=r8), errstat)

  if (errstat < 0) then
    call tell_error (tell_io_write_error, "r8 write failed", errstat);
    return
  endif

  call test_put_errors (obj, max_dims, errstat)

  call create_group (obj, errstat)
  if (errstat < 0) then
    call tell_error (tell_runtime_error, "failed creating group", errstat)
    return
  endif

  call tiof_close (obj, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_error, "close failed", errstat);
    return
  endif

  call tiof_dimlist_free (dimlist)
  call tiof_varlist_free (varlist)

end subroutine

subroutine check_read (filename, values, max_dims, start, edge, &
                       fillvalue, strings, errstat)
  use tell_module
  use tio_module
  use netcdf
  implicit none

  character (len=*), intent(in) :: filename
  real (kind=r8), dimension(:), intent(in) :: values
  integer, intent(in) :: max_dims
  integer, dimension(:), intent(in) :: start, edge
  real (kind=r8), intent(in) :: fillvalue
  character (len=*), dimension(:), intent(in) :: strings
  integer, intent(inout) :: errstat

  type (tiof_file_type) :: obj

  integer (kind=i1) :: i1s
  integer (kind=i2) :: i2s
  integer (kind=i4) :: i4s
  integer (kind=i8) :: i8s
  real (kind=r4) :: r4s
  real (kind=r8) :: r8s

  character (len=8), dimension(3) :: c1d
  !character (len=:), dimension(:), allocatable :: c1d
  integer :: len_text, dim_text

  integer (kind=i1), allocatable :: &
    i1_1d(:), i1_2d(:,:), i1_3d(:,:,:)
  integer (kind=i2), allocatable :: &
    i2_1d(:), i2_2d(:,:), i2_3d(:,:,:), i2_4d(:,:,:,:)
  integer (kind=i4), allocatable :: &
    i4_1d(:), i4_2d(:,:), i4_3d(:,:,:)
  integer (kind=i8), allocatable :: &
    i8_1d(:), i8_2d(:,:), i8_3d(:,:,:)
  real (kind=r4), allocatable :: &
    r4_1d(:), r4_2d(:,:), r4_3d(:,:,:), &
    r4_4d(:,:,:,:), r4_5d(:,:,:,:,:), r4_6d(:,:,:,:,:,:)
  real (kind=r8), allocatable :: &
    r8_1d(:), r8_2d(:,:), r8_3d(:,:,:), &
    r8_4d(:,:,:,:), r8_5d(:,:,:,:,:), r8_6d(:,:,:,:,:,:)

  real (kind=r4), parameter :: r4_huge_1 = huge(1.0_4)

  call tiof_open (filename, obj, nf90_nowrite, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_error, "open failed", errstat);
    return
  endif

  i1s = -1
  call tiof_get_i1 (obj, "sbS", i1s, errstat)
  if (int(values(1),kind=i1) /= i1s) then
    call tell_error (tell_runtime_error, "i1s value mismatch", errstat)
    return
  endif

  i2s = -1
  call tiof_get_i2 (obj, "ssS", i2s, errstat)
  if (int(values(1),kind=i2) /= i2s) then
    call tell_error (tell_runtime_error, "i2s value mismatch", errstat)
    return
  endif

  i4s = -1
  call tiof_get_i4 (obj, "siS", i4s, errstat)
  if (int(values(1),kind=i4) /= i4s) then
    call tell_error (tell_runtime_error, "i4s value mismatch", errstat)
    return
  endif

  i8s = -1
  call tiof_get_i8 (obj, "sjS", i8s, errstat)
  if (int(values(1),kind=i8) /= i8s) then
    call tell_error (tell_runtime_error, "i8s value mismatch", errstat)
    return
  endif

  i1s = -1
  call tiof_get_ui1 (obj, "ubS", i1s, errstat)
  if (int(values(1),kind=i1) /= i1s) then
    call tell_error (tell_runtime_error, "ui1s value mismatch", errstat)
    return
  endif

  i2s = -1
  call tiof_get_ui2 (obj, "usS", i2s, errstat)
  if (int(values(1),kind=i2) /= i2s) then
    call tell_error (tell_runtime_error, "ui2s value mismatch", errstat)
    return
  endif

  i4s = -1
  call tiof_get_ui4 (obj, "uiS", i4s, errstat)
  if (int(values(1),kind=i4) /= i4s) then
    call tell_error (tell_runtime_error, "ui4s value mismatch", errstat)
    return
  endif

  i8s = -1
  call tiof_get_ui8 (obj, "ujS", i8s, errstat)
  if (int(values(1),kind=i8) /= i8s) then
    call tell_error (tell_runtime_error, "ui8s value mismatch", errstat)
    return
  endif

  r4s = -1
  call tiof_get_r4 (obj, "rfS", r4s, errstat)
  if (real(values(1),kind=r4) /= r4s) then
    call tell_error (tell_runtime_error, "r4s value mismatch", errstat)
    return
  endif

  r8s = -1
  call tiof_get_r8 (obj, "rdS", r8s, errstat)
  if (real(values(1),kind=r8) /= r8s) then
    call tell_error (tell_runtime_error, "r8s value mismatch", errstat)
    return
  endif

  call tiof_inq_dimlen (obj, "n6", len_text, errstat)
  call tiof_inq_dimlen (obj, "n1", dim_text, errstat)
  if (errstat < 0) then
    call tell_error (tell_runtime_error, &
                     "failed reading text array dimensions", errstat);
    return
  endif

  !allocate (character (len=len_text) :: c1d(dim_text), stat=errstat)
  !if (errstat /= 0) then
  !  call tell_error (tell_malloc_error, "allocate failed", errstat);
  !  return
  !endif
  ! FIXME: compiler bug?? gfortran prints len_text=8, len(cld(1))=0 !!
  !write(*,*)len_text, len(c1d(1))

  allocate (i1_1d(edge(1)), &
            i1_2d(edge(1),edge(2)), &
            i1_3d(edge(1),edge(2),edge(3)), &
            i2_1d(edge(1)), &
            i2_2d(edge(1),edge(2)), &
            i2_3d(edge(1),edge(2),edge(3)), &
            i2_4d(edge(1),edge(2),edge(3),edge(4)), &
            i4_1d(edge(1)), &
            i4_2d(edge(1),edge(2)), &
            i4_3d(edge(1),edge(2),edge(3)), &
            i8_1d(edge(1)), &
            i8_2d(edge(1),edge(2)), &
            i8_3d(edge(1),edge(2),edge(3)), &
            r4_1d(edge(1)), &
            r4_2d(edge(1),edge(2)), &
            r4_3d(edge(1),edge(2),edge(3)), &
            r4_4d(edge(1),edge(2),edge(3),edge(4)), &
            r4_5d(edge(1),edge(2),edge(3),edge(4),edge(5)), &
            r4_6d(edge(1),edge(2),edge(3),edge(4),edge(5),edge(6)), &
            r8_1d(edge(1)), &
            r8_2d(edge(1),edge(2)), &
            r8_3d(edge(1),edge(2),edge(3)), &
            r8_4d(edge(1),edge(2),edge(3),edge(4)), &
            r8_5d(edge(1),edge(2),edge(3),edge(4),edge(5)), &
            r8_6d(edge(1),edge(2),edge(3),edge(4),edge(5),edge(6)), &
            stat=errstat)
  if (errstat /= 0) then
    call tell_error (tell_malloc_error, "allocate failed", errstat);
    return
  endif

  ! text array
  c1d(:) = "X"
  call tiof_get1d_text (obj, "txt1d", start(1), edge(1), &
                        c1d, errstat)
  if (any (strings /= c1d)) then
    call tell_error (tell_runtime_error, "txt1d value mismatch", errstat)
    return
  endif

  ! string array
  c1d(:) = "X"
  call tiof_get1d_string (obj, "str1d", start(1), edge(1), &
                          c1d, errstat)
  if (any (strings /= c1d)) then
    call tell_error (tell_runtime_error, "str1d value mismatch", errstat)
    return
  endif

  ! integer (kind=i1)
  i1_1d(:) = -1
  call tiof_get1d_i1 (obj, "sbA1", start(1:1), edge(1:1), i1_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i1))
  if (any (int(values(1:edge(1)),kind=i1) /= i1_1d)) then
    call tell_error (tell_runtime_error, "i1_1d value mismatch", errstat)
    return
  endif
  i1_1d(:) = -1
  call tiof_get1d_ui1 (obj, "ubA1", start(1:1), edge(1:1), i1_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i1))
  if (any (int(values(1:edge(1)),kind=i1) /= i1_1d)) then
    call tell_error (tell_runtime_error, "i1_1d u-value mismatch", errstat)
    return
  endif

  i1_2d(:,:) = -1
  call tiof_get2d_i1 (obj, "sbA2", start(2:1:-1), edge(2:1:-1), i1_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i1))
  if (any (int(values(1:product(edge(1:2))),kind=i1) &
           /= reshape(i1_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i1_2d value mismatch", errstat)
    return
  endif
  i1_2d(:,:) = -1
  call tiof_get2d_ui1 (obj, "ubA2", start(2:1:-1), edge(2:1:-1), i1_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i1))
  if (any (int(values(1:product(edge(1:2))),kind=i1) &
           /= reshape(i1_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i1_2d u-value mismatch", errstat)
    return
  endif

  i1_3d(:,:,:) = -1
  call tiof_get3d_i1 (obj, "sbA3", start(3:1:-1), edge(3:1:-1), i1_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i1))
  if (any (int(values(1:product(edge(1:3))),kind=i1) &
           /= reshape(i1_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i1_3d value mismatch", errstat)
    return
  endif
  i1_3d(:,:,:) = -1
  call tiof_get3d_ui1 (obj, "ubA3", start(3:1:-1), edge(3:1:-1), i1_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i1))
  if (any (int(values(1:product(edge(1:3))),kind=i1) &
           /= reshape(i1_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i1_3d u-value mismatch", errstat)
    return
  endif

  ! integer (kind=i2)
  i2_1d(:) = -1
  call tiof_get1d_i2 (obj, "ssA1", start(1:1), edge(1:1), i2_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:edge(1)),kind=i2) /= i2_1d)) then
    call tell_error (tell_runtime_error, "i2_1d value mismatch", errstat)
    return
  endif
  i2_1d(:) = -1
  call tiof_get1d_ui2 (obj, "usA1", start(1:1), edge(1:1), i2_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:edge(1)),kind=i2) /= i2_1d)) then
    call tell_error (tell_runtime_error, "i2_1d u-value mismatch", errstat)
    return
  endif

  i2_2d(:,:) = -1
  call tiof_get2d_i2 (obj, "ssA2", start(2:1:-1), edge(2:1:-1), i2_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:product(edge(1:2))),kind=i2) &
           /= reshape(i2_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i2_2d value mismatch", errstat)
    return
  endif
  i2_2d(:,:) = -1
  call tiof_get2d_ui2 (obj, "usA2", start(2:1:-1), edge(2:1:-1), i2_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:product(edge(1:2))),kind=i2) &
           /= reshape(i2_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i2_2d u-value mismatch", errstat)
    return
  endif

  i2_3d(:,:,:) = -1
  call tiof_get3d_i2 (obj, "ssA3", start(3:1:-1), edge(3:1:-1), i2_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:product(edge(1:3))),kind=i2) &
           /= reshape(i2_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i2_3d value mismatch", errstat)
    return
  endif

  i2_4d(:,:,:,:) = -1
  call tiof_get4d_i2 (obj, "ssA4", start(4:1:-1), edge(4:1:-1), i2_4d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:product(edge(1:4))),kind=i2) &
           /= reshape(i2_4d, [product(edge(1:4))]))) then
    call tell_error (tell_runtime_error, "i2_4d value mismatch", errstat)
    return
  endif

  i2_3d(:,:,:) = -1
  call tiof_get3d_ui2 (obj, "usA3", start(3:1:-1), edge(3:1:-1), i2_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i2))
  if (any (int(values(1:product(edge(1:3))),kind=i2) &
           /= reshape(i2_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i2_3d u-value mismatch", errstat)
    return
  endif

  ! integer (kind=i4)
  i4_1d(:) = huge(1)
  call tiof_get1d_i4 (obj, "siA1", start(1:1), edge(1:1), i4_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i4))
  if (any (int(values(1:edge(1)),kind=i4) /= i4_1d)) then
    call tell_error (tell_runtime_error, "i4_1d value mismatch", errstat)
    return
  endif
  i4_1d(:) = huge(1)
  call tiof_get1d_ui4 (obj, "uiA1", start(1:1), edge(1:1), i4_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i4))
  if (any (int(values(1:edge(1)),kind=i4) /= i4_1d)) then
    call tell_error (tell_runtime_error, "i4_1d u-value mismatch", errstat)
    return
  endif

  i4_2d(:,:) = huge(1)
  call tiof_get2d_i4 (obj, "siA2", start(2:1:-1), edge(2:1:-1), i4_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i4))
  if (any (int(values(1:product(edge(1:2))),kind=i4) &
           /= reshape(i4_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i4_2d value mismatch", errstat)
    return
  endif
  i4_2d(:,:) = huge(1)
  call tiof_get2d_ui4 (obj, "uiA2", start(2:1:-1), edge(2:1:-1), i4_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i4))
  if (any (int(values(1:product(edge(1:2))),kind=i4) &
           /= reshape(i4_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i4_2d u-value mismatch", errstat)
    return
  endif

  i4_3d(:,:,:) = huge(1)
  call tiof_get3d_i4 (obj, "siA3", start(3:1:-1), edge(3:1:-1), i4_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i4))
  if (any (int(values(1:product(edge(1:3))),kind=i4) &
           /= reshape(i4_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i4_3d value mismatch", errstat)
    return
  endif
  i4_3d(:,:,:) = huge(1)
  call tiof_get3d_ui4 (obj, "uiA3", start(3:1:-1), edge(3:1:-1), i4_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i4))
  if (any (int(values(1:product(edge(1:3))),kind=i4) &
           /= reshape(i4_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i4_3d u-value mismatch", errstat)
    return
  endif

  ! integer (kind=i8)
  i8_1d(:) = huge(1)
  call tiof_get1d_i8 (obj, "sjA1", start(1:1), edge(1:1), i8_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i8))
  if (any (int(values(1:edge(1)),kind=i8) /= i8_1d)) then
    call tell_error (tell_runtime_error, "i8_1d value mismatch", errstat)
    return
  endif
  i8_1d(:) = huge(1)
  call tiof_get1d_ui8 (obj, "ujA1", start(1:1), edge(1:1), i8_1d, errstat, &
                     replace_fill=int(fillvalue,kind=i8))
  if (any (int(values(1:edge(1)),kind=i8) /= i8_1d)) then
    call tell_error (tell_runtime_error, "i8_1d u-value mismatch", errstat)
    return
  endif

  i8_2d(:,:) = huge(1)
  call tiof_get2d_i8 (obj, "sjA2", start(2:1:-1), edge(2:1:-1), i8_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i8))
  if (any (int(values(1:product(edge(1:2))),kind=i8) &
           /= reshape(i8_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i8_2d value mismatch", errstat)
    return
  endif
  i8_2d(:,:) = huge(1)
  call tiof_get2d_ui8 (obj, "ujA2", start(2:1:-1), edge(2:1:-1), i8_2d, errstat, &
                     replace_fill=int(fillvalue,kind=i8))
  if (any (int(values(1:product(edge(1:2))),kind=i8) &
           /= reshape(i8_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "i8_2d u-value mismatch", errstat)
    return
  endif

  i8_3d(:,:,:) = huge(1)
  call tiof_get3d_i8 (obj, "sjA3", start(3:1:-1), edge(3:1:-1), i8_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i8))
  if (any (int(values(1:product(edge(1:3))),kind=i8) &
           /= reshape(i8_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i8_3d value mismatch", errstat)
    return
  endif
  i8_3d(:,:,:) = huge(1)
  call tiof_get3d_ui8 (obj, "ujA3", start(3:1:-1), edge(3:1:-1), i8_3d, errstat, &
                     replace_fill=int(fillvalue,kind=i8))
  if (any (int(values(1:product(edge(1:3))),kind=i8) &
           /= reshape(i8_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "i8_3d u-value mismatch", errstat)
    return
  endif

  ! real (kind=r4)
  r4_1d(:) = r4_huge_1
  call tiof_get1d_r4 (obj, "rfA1", start(1:1), edge(1:1), r4_1d, errstat, &
                     replace_fill=real(fillvalue,kind=r4))
  if (any (real(values(1:edge(1)),kind=r4) /= r4_1d)) then
    call tell_error (tell_runtime_error, "r4_1d value mismatch", errstat)
    return
  endif

  r4_2d(:,:) = r4_huge_1
  call tiof_get2d_r4 (obj, "rfA2", start(2:1:-1), edge(2:1:-1), r4_2d, errstat, &
                     replace_fill=real(fillvalue,kind=r4))
  if (any (real(values(1:product(edge(1:2))),kind=r4) &
           /= reshape(r4_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "r4_2d value mismatch", errstat)
    return
  endif

  r4_3d(:,:,:) = r4_huge_1
  call tiof_get3d_r4 (obj, "rfA3", start(3:1:-1), edge(3:1:-1), r4_3d, errstat, &
                     replace_fill=real(fillvalue,kind=r4))
  if (any (real(values(1:product(edge(1:3))),kind=r4) &
           /= reshape(r4_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "r4_3d value mismatch", errstat)
    return
  endif

  r4_4d(:,:,:,:) = r4_huge_1
  call tiof_get4d_r4 (obj, "rfA4", start(4:1:-1), edge(4:1:-1), r4_4d, errstat, &
                     replace_fill=real(fillvalue,kind=r4))
  if (any (real(values(1:product(edge(1:4))),kind=r4) &
           /= reshape(r4_4d, [product(edge(1:4))]))) then
    call tell_error (tell_runtime_error, "r4_4d value mismatch", errstat)
    return
  endif

  r4_5d(:,:,:,:,:) = r4_huge_1
  call tiof_get5d_r4 (obj, "rfA5", start(5:1:-1), edge(5:1:-1), r4_5d, errstat, &
                     replace_fill=real(fillvalue,kind=r4))
  if (any (real(values(1:product(edge(1:5))),kind=r4) &
           /= reshape(r4_5d, [product(edge(1:5))]))) then
    call tell_error (tell_runtime_error, "r4_5d value mismatch", errstat)
    return
  endif

  r4_6d(:,:,:,:,:,:) = r4_huge_1
  call tiof_get6d_r4 (obj, "rfA6", start(6:1:-1), edge(6:1:-1), r4_6d, errstat, &
                     replace_fill=real(fillvalue,kind=r4))
  if (any (real(values(1:product(edge(1:6))),kind=r4) &
           /= reshape(r4_6d, [product(edge(1:6))]))) then
    call tell_error (tell_runtime_error, "r4_6d value mismatch", errstat)
    return
  endif

  ! real (kind=r8)
  r8_1d(:) = huge(1)
  call tiof_get1d_r8 (obj, "rdA1", start(1:1), edge(1:1), r8_1d, errstat, &
                     replace_fill=real(fillvalue,kind=r8))
  if (any (real(values(1:edge(1)),kind=r8) /= r8_1d)) then
    call tell_error (tell_runtime_error, "r8_1d value mismatch", errstat)
    return
  endif

  r8_2d(:,:) = huge(1)
  call tiof_get2d_r8 (obj, "rdA2", start(2:1:-1), edge(2:1:-1), r8_2d, errstat, &
                     replace_fill=real(fillvalue,kind=r8))
  if (any (real(values(1:product(edge(1:2))),kind=r8) &
           /= reshape(r8_2d, [product(edge(1:2))]))) then
    call tell_error (tell_runtime_error, "r8_2d value mismatch", errstat)
    return
  endif

  r8_3d(:,:,:) = huge(1)
  call tiof_get3d_r8 (obj, "rdA3", start(3:1:-1), edge(3:1:-1), r8_3d, errstat, &
                     replace_fill=real(fillvalue,kind=r8))
  if (any (real(values(1:product(edge(1:3))),kind=r8) &
           /= reshape(r8_3d, [product(edge(1:3))]))) then
    call tell_error (tell_runtime_error, "r8_3d value mismatch", errstat)
    return
  endif

  r8_4d(:,:,:,:) = huge(1)
  call tiof_get4d_r8 (obj, "rdA4", start(4:1:-1), edge(4:1:-1), r8_4d, errstat, &
                     replace_fill=real(fillvalue,kind=r8))
  if (any (real(values(1:product(edge(1:4))),kind=r8) &
           /= reshape(r8_4d, [product(edge(1:4))]))) then
    call tell_error (tell_runtime_error, "r8_4d value mismatch", errstat)
    return
  endif

  r8_5d(:,:,:,:,:) = huge(1)
  call tiof_get5d_r8 (obj, "rdA5", start(5:1:-1), edge(5:1:-1), r8_5d, errstat, &
                     replace_fill=real(fillvalue,kind=r8))
  if (any (real(values(1:product(edge(1:5))),kind=r8) &
           /= reshape(r8_5d, [product(edge(1:5))]))) then
    call tell_error (tell_runtime_error, "r8_5d value mismatch", errstat)
    return
  endif

  r8_6d(:,:,:,:,:,:) = huge(1)
  call tiof_get6d_r8 (obj, "rdA6", start(6:1:-1), edge(6:1:-1), r8_6d, errstat, &
                     replace_fill=real(fillvalue,kind=r8))
  if (any (real(values(1:product(edge(1:6))),kind=r8) &
           /= reshape(r8_6d, [product(edge(1:6))]))) then
    call tell_error (tell_runtime_error, "r8_6d value mismatch", errstat)
    return
  endif

  call test_get_errors (obj, max_dims, errstat)

  call read_from_group (obj, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_read_error, "failed reading from group", errstat)
    return
  endif

  call tiof_close (obj, errstat)
  if (errstat < 0) then
    call tell_error (tell_io_error, "close failed", errstat)
    return
  endif

end subroutine check_read

end module ftest1_module
