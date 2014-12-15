module tio_module
  use netcdf
  use tell_module
  implicit none

  integer, private, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4)

  integer :: tiof_get_l1bvar, tiof_put_att1
  external   tiof_get_l1bvar, tiof_put_att1

  type, public :: tiof_l1_object_type
    integer :: fileid = -1
    integer :: groupid = -1
    integer :: num_times=0, num_xtrack=0, num_wavelengths=0
  end type tiof_l1_object_type

  private

  public tiof_open, tiof_close, tiof_inq_group, &
    tiof_get1d_r8, &
    tiof_get1d_r4, tiof_get2d_r4, tiof_get3d_r4, &
    tiof_get2d_i2, tiof_get3d_i2, &
    tiof_get1d_i1, tiof_get2d_i1

contains

  subroutine get_dimlen (grp, name, dimlen, errstat)
    implicit none
    integer, intent(in) :: grp
    character (len=*), intent(in) :: name
    integer, intent(out) :: dimlen
    integer, intent(inout) :: errstat

    integer :: status, dimid

    if (errstat < 0) return

    status = nf90_inq_dimid (grp, name, dimid)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    status = nf90_inquire_dimension (grp, dimid, len=dimlen)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif
  end subroutine get_dimlen

  subroutine tiof_open (file, l1bobj, errstat)
    implicit none
    character (len=*), intent(in) :: file
    type (tiof_l1_object_type), intent(out) :: l1bobj
    integer, intent(inout) :: errstat

    integer :: fileid, status

    if (errstat < 0) return

    status = nf90_open (file, nf90_nowrite, fileid)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_open_error, "opening file "//file//" ("//trim(nf90_strerror(status))//")", errstat)
      l1bobj % fileid = -1
      return
    endif

    l1bobj % fileid = fileid
    l1bobj % groupid = fileid

    call get_dimlen (fileid, "mirror_step", l1bobj % num_times, errstat)
    if (errstat < 0) return
  end subroutine tiof_open

  subroutine tiof_close (l1bobj, errstat)
    implicit none
    type (tiof_l1_object_type), intent(inout) :: l1bobj
    integer, intent(inout) :: errstat

    integer :: status

    if (errstat < 0) return

    if (l1bobj%fileid >= 0) then
      status = nf90_close (l1bobj % fileid)
      if (status /= nf90_noerr) then
        call tell_error (tell_io_error, "closing file ("//trim(nf90_strerror(status))//")", errstat)
      endif
      l1bobj % fileid = -1
      l1bobj % groupid = -1
    endif
  end subroutine tiof_close

  subroutine tiof_inq_group (l1bobj, grpname, errstat)
    implicit none
    type (tiof_l1_object_type), intent(inout) :: l1bobj
    character (len=*), intent(in) :: grpname
    integer, intent(inout) :: errstat

    integer :: status, grp

    if (errstat < 0) return

    status = nf90_inq_ncid (l1bobj % fileid, grpname, grp)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing group "//grpname//" ("//trim(nf90_strerror(status))//")", errstat)
    endif
    l1bobj % groupid = grp

    call get_dimlen (grp, "xtrack", l1bobj % num_xtrack, errstat)
    if (errstat < 0) return

    call get_dimlen (grp, "spectral_channel", l1bobj % num_wavelengths, errstat)
    if (errstat < 0) return
  end subroutine tiof_inq_group

  subroutine tiof_get1d_r8 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=8), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_double, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get1d_r8

  subroutine tiof_get3d_r4 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get3d_r4

  subroutine tiof_get2d_r4 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get2d_r4

  subroutine tiof_get1d_r4 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get1d_r4

  subroutine tiof_get3d_i2 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i2), dimension (:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get3d_i2

  subroutine tiof_get2d_i2 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i2), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get2d_i2

  subroutine tiof_get2d_i1 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i1), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_byte, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get2d_i1

  subroutine tiof_get1d_i1 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_l1_object_type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i1), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_byte, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get1d_i1

end module tio_module
