module l1bread

  use errormodule
  use hdfeos4_parameters ! , only: DFNT_INT16, DFNT_FLOAT32, DFNT_FLOAT64
  implicit none

  integer (KIND = 4) :: newgetl1bblk
  external              newgetl1bblk

  type, public :: L1B_Object_Type
    integer (kind=4) :: fileid = -1
    integer (kind=4) :: swathid = -1
    integer (kind=4) :: num_times=0, num_xtrack=0, num_wavelengths=0
    integer (kind=4) :: start_line = -1, end_line = -1
  end type L1B_Object_Type

  private
  public l1bread_open, l1bread_open_swath, l1bread_attach_swatch, &
    l1bread_close, l1bread_radiance_info, l1bread_swathname, &
    l1bread_get1d_r8, &
    l1bread_get1d_r4, l1bread_get2d_r4, l1bread_get3d_r4, &
    l1bread_get2d_i2, l1bread_get3d_i2, &
    l1bread_get1d_i1, l1bread_get2d_i1

contains

  subroutine do_swrdattr_i4 (swathid, name, value, errstat)

    implicit none
    integer (kind=4), intent(in) :: swathid
    character (len=*), intent(in) :: name
    integer (kind=4), intent(out) :: value
    integer, intent(inout) :: errstat
    integer :: err

    if (errstat < 0) return

    value = swdiminfo (swathid, name)
    if (value < 0) then
      err = swrdattr (swathid, name, value)
      if (err /= 0) then
        errstat = -1
        call err_message_error ("swrdattr/swdiminfo " // name // " failed", errstat)
      endif
    endif

    end subroutine

  subroutine l1bread_open (file, l1bobj, errstat)

    implicit none
    character (len=*), intent(in) :: file
    type (L1B_Object_Type), intent(out) :: l1bobj
    integer, intent(inout) :: errstat

    integer (kind=4) :: fileid

    if (errstat < 0) return

    fileid = swopen (file, DFACC_READ)
    if (fileid < 0) then
      call err_message_error ("swopen "//file//" failed", errstat)
      fileid = -1
      ! drop
    endif

    l1bobj%fileid = fileid
    l1bobj%swathid = -1
    return

  end subroutine

  subroutine l1bread_attach_swatch (l1bobj, swathname, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent(in) :: swathname
    integer, intent(inout) :: errstat

    integer (kind=4) :: swathid
    integer :: err

    if (l1bobj%swathid /= -1) then
      err = swclose (l1bobj%swathid)
      l1bobj%swathid = -1
    endif

    swathid = swattach (l1bobj%fileid, swathname)
    if (swathid < 0) then
      call err_message_error ("swattach " // swathname // " failed", errstat)
      return
    endif

    call do_swrdattr_i4 (swathid, "NumTimes", l1bobj%num_times, errstat)
    call do_swrdattr_i4 (swathid, "nXtrack", l1bobj%num_xtrack, errstat)
    call do_swrdattr_i4 (swathid, "nWavel", l1bobj%num_wavelengths, errstat)
    if (errstat < 0) goto 666

    l1bobj%swathid = swathid
    return

 666 continue
    if (swathid >= 0) err = swdetach (swathid)

  end subroutine

  subroutine l1bread_open_swath (file, swathname, l1bobj, errstat)

    use hdfeos4_parameters
    use errormodule

    implicit none
    character (len=*), intent(in) :: file, swathname
    type (L1B_Object_Type), intent(out) :: l1bobj
    integer, intent(inout) :: errstat

    write (*,*) "l1bread_open_swath ", trim(file), trim(swathname)

    ! allow errors to flow through
    call l1bread_open (file, l1bobj, errstat);
    call l1bread_attach_swatch (l1bobj, swathname, errstat)

  end subroutine

  subroutine l1bread_close (l1bobj)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    integer err

    if (l1bobj % swathid >= 0) then
      err = swdetach (l1bobj % swathid)
      l1bobj % swathid = -1
    endif
    if (l1bobj % fileid >= 0) then
      err = swclose (l1bobj % fileid)
      l1bobj % fileid = -1
    endif
  end subroutine

  !subroutine check_array_size (l1bobj, name, array, errstat)
  !
  !  implicit none
  !  type (L1B_Object_Type), intent(inout) :: l1bobj
  !  character (len=*), intent (in) :: name
  !  integer (kind=4), dimension(:), intent(in) :: array
  !  integer, intent(inout) :: errstat
  !
  !  integer (kind=8) :: num, min_num
  !
  !  if (errstat < 0) return
  !
  !  ! FIXME: This routine should check all dimensions
  !  num = size(array, kind=8)
  !
  !  if (num < min_num) &
  !    call err_message_error ("check_array_size: " // name // &
  !                            " array size is too small", errstat)
  !
  !end subroutine

  subroutine l1bread_get1d_r8 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    real (kind=8), dimension(:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_FLOAT64, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get3d_r4 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    real (kind=4), dimension(:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_FLOAT32, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get2d_r4 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    real (kind=4), dimension(:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_FLOAT32, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get1d_r4 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    real (kind=4), dimension(:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_FLOAT32, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get3d_i2 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    integer (kind=2), dimension(:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_INT16, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get2d_i2 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    integer (kind=2), dimension(:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_INT16, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get2d_i1 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    integer (kind=1), dimension(:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_INT8, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_get1d_i1 (l1bobj, name, tstep, numsteps, array, errstat)

    implicit none
    type (L1B_Object_Type), intent(inout) :: l1bobj
    character (len=*), intent (in) :: name
    integer (kind=4), intent (in) :: tstep, numsteps
    integer (kind=1), dimension(:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer err

    if (errstat < 0) return

    err = newgetl1bblk (l1bobj%swathid, name, tstep, numsteps, &
                        DFNT_INT8, array)

    if (err < 0) then
      call err_message_error ("Unable to read " // name // " from L1b file", errstat)
      return
    endif

  end subroutine

  subroutine l1bread_swathname (l1bfile, l1bchan, swathname, errstat)

    implicit none
    character (len=*), intent (in) :: l1bfile
    character (len=*), intent (in) :: l1bchan
    character (len=*), intent (out) :: swathname
    integer, intent(inout) :: errstat
    !
    integer (kind=4) :: strbufsize, nswath
    integer :: i0, i, i1
    character (len=1024) :: swathlist
    character (len=32) local_l1bchan

    if (errstat < 0) return

    swathlist=""
    nswath = SWInqswath ( l1bfile, swathlist, strbufsize )
    if (nswath < 0) then
      call err_message_error ("SWInqswath failed for file "//l1bfile, errstat);
      return
    endif

    if (nswath == 0) then
      call err_message_error ("SWInqswath returned no swaths for "//l1bfile, &
                         errstat)
      return
    endif

    swathlist = swathlist(1:strbufsize)

    ! special omi hack
    local_l1bchan = l1bchan
    if (l1bchan == "UV1") local_l1bchan = "UV-1"
    if (l1bchan == "UV2") local_l1bchan = "UV-2"

    i0 = 1
    i1 = strbufsize
    do while (i0 < i1)
      i = index (swathlist(i0:i1), ",")
      if (i == 0) i = i1 - i0 + 2
      i = i0 + i - 1
      swathname = swathlist (i0:i-1)
      ! How robust is this?  Here we just look for, e.g., UV-1 in the name
      ! write (*,*) "swathname=", trim(swathname)
      if (0 /= index (trim(swathname), trim(local_l1bchan))) then
        return
      endif
      i0 = i+1
    end do

    call err_message_warn ("No suitable swathname found for l1bchan="//l1bchan);
    swathname="?"
    return
  end subroutine l1bread_swathname

  subroutine l1bread_radiance_info (l1bfile, l1bchannel, rpt, errstat)

    USE OMSAO_precision_module
    !USE OMSAO_parameters_module, ONLY : MAX_STR_LEN
    USE OMSAO_variables_module,  ONLY : Radiance_Paras_Type
    USE OMSAO_errstat_module
    USE hdfeos4_parameters

    implicit none
    character (len=*), intent (in) :: l1bfile, l1bchannel
    type(Radiance_Paras_Type), INTENT(out) :: rpt
    integer (kind=i4), intent (inout) :: errstat

    type (L1B_Object_Type) :: l1bobj

    if (errstat < 0) return

    rpt%ntimes = 0 ; rpt%nxtrack = 0 ; rpt%nwavel_ccd = 0
    rpt%l1bfilename = l1bfile
    rpt%l1bchannel = l1bchannel

    ! allow error to flow through
    call l1bread_swathname (l1bfile, l1bchannel, rpt%swathname, errstat)
    call l1bread_open_swath (l1bfile, rpt%swathname, l1bobj, errstat)
    if (errstat < 0) goto 666

    rpt%ntimes = l1bobj%num_times
    rpt%nxtrack = l1bobj%num_xtrack
    rpt%nwavel_ccd = l1bobj%num_wavelengths

    write (*,*) "DEBUG: In l1bread_radiance_info, l1bfile=", &
      trim(l1bfile), ", l1bswath=", trim(rpt%swathname)

    ! drop

 666 call l1bread_close (l1bobj)

  end subroutine l1bread_radiance_info

end module
