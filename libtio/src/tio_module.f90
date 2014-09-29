module tio_module
  use netcdf
  implicit none

  integer (kind=4) :: tiof_get_l1bvar
  external            tiof_get_l1bvar

  type, public :: L1B_Object_Type
    integer (kind=4) :: fileid = -1
    integer (kind=4) :: groupid = -1
    integer (kind=4) :: num_times=0, num_xtrack=0, num_wavelengths=0
  end type L1B_Object_Type

  private

  public tiof_open, tiof_close, tiof_inq_group, &
    tiof_get1d_r8, &
    tiof_get1d_r4, tiof_get2d_r4, tiof_get3d_r4, &
    tiof_get2d_i2, tiof_get3d_i2, &
    tiof_get1d_i1, tiof_get2d_i1

contains

    subroutine err_message_error (msg, errcode)

    implicit none
    character (len=*), intent(in) :: msg
    integer, intent(inout) :: errcode

    write (*,*) "ERROR: ", trim(msg)
    if (errcode >= 0) errcode = -1
    return
  end subroutine err_message_error

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
      call err_message_error ("accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    status = nf90_inquire_dimension (grp, dimid, len=dimlen)
    if (status /= nf90_noerr) then
      call err_message_error ("accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif
  end subroutine get_dimlen

  subroutine tiof_open (file, l1bobj, errstat)
    implicit none
    character (len=*), intent(in) :: file
    type (L1B_Object_Type), intent(out) :: l1bobj
    integer, intent(inout) :: errstat

    integer (kind=4) :: fileid
    integer :: status

    if (errstat < 0) return

    status = nf90_open (file, nf90_nowrite, fileid)
    if (status /= nf90_noerr) then
      call err_message_error ("opening file "//file//" ("//trim(nf90_strerror(status))//")", errstat)
      l1bobj % fileid = -1
      return
    endif

    l1bobj % fileid = fileid

    call get_dimlen (fileid, "mirror_step", l1bobj % num_times, errstat)
    if (errstat < 0) return
  end subroutine tiof_open

  subroutine tiof_close (l1bobj, errstat)
    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    integer, intent(inout) :: errstat

    integer :: status

    if (l1bobj%fileid >= 0) then
      status = nf90_close (l1bobj % fileid)
      if (status /= nf90_noerr) then
        call err_message_error ("closing file ("//trim(nf90_strerror(status))//")", errstat)
      endif
      l1bobj % fileid = -1
      l1bobj % groupid = -1
    endif
  end subroutine tiof_close

  subroutine tiof_inq_group (l1bobj, grpname, errstat)
    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent(in) :: grpname
    integer, intent(inout) :: errstat

    integer :: status, grp

    if (errstat < 0) return

    status = nf90_inq_ncid (l1bobj % fileid, grpname, grp)
    if (status /= nf90_noerr) then
      call err_message_error ("accessing group "//grpname//" ("//trim(nf90_strerror(status))//")", errstat)
    endif
    l1bobj % groupid = grp

    call get_dimlen (grp, "xtrack", l1bobj % num_xtrack, errstat)
    if (errstat < 0) return

    call get_dimlen (grp, "spectral_channel", l1bobj % num_wavelengths, errstat)
    if (errstat < 0) return
  end subroutine tiof_inq_group

  subroutine tiof_get1d_r8 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=8), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_double, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get1d_r8

  subroutine tiof_get3d_r4 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get3d_r4

  subroutine tiof_get2d_r4 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get2d_r4

  subroutine tiof_get1d_r4 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get1d_r4

  subroutine tiof_get3d_i2 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=2), dimension (:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get3d_i2

  subroutine tiof_get2d_i2 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=2), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get2d_i2

  subroutine tiof_get2d_i1 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=1), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get2d_i1

  subroutine tiof_get1d_i1 (l1bobj, name, step0, numsteps, array, errstat)
    implicit none
    type (L1B_Object_Type), intent(in) :: l1bobj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=1), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_l1bvar (l1bobj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif
  end subroutine tiof_get1d_i1

end module tio_module
