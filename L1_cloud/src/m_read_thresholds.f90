!>Read thresholds reference file needed for cloud mask computation
module m_read_thresholds

  private
  public read_thresholds

contains

  !-------------------------------------------------------------------------
  !
  ! !USES:
  !  read_thresholds reads precomputed thresholds needed for cloud
  !               mask
  !
  ! !OUTPUT PARAMETERS:
  !
  !> @param errstat Error return code, non-zero indicates problem
  ! !DESCRIPTION: 
  !
  ! !REVISION HISTORY:
  !
  !> @author  25aug04   Joiner      original fortran 90
  !> @author  26mar15   O'Sullivan  updating for TEMPO
  !
  !-------------------------------------------------------------------------
  subroutine read_thresholds (errstat)

    use m_vars, ONLY: npixs, nscanpos, thresholds, npixels, ex, &
         stddev_thresh 
    use m_LUN_set
    use m_pgs_include
    use tell_module
    Implicit NONE

    ! !INPUT PARAMETERS:
    !
    ! !OUTPUT PARAMETERS:
    !
    integer, intent(inout)         :: errstat        ! Error return code:
    !  0   all is well
    !  -1   files not found
    !-------------------------------------------------------------------------
    integer :: i, j, l, m
    character(len=100) :: text, logmsg
    integer :: lun=10
    integer :: pgs_io_gen_openf, pgs_io_gen_closef
    integer :: status, version

    if (errstat /= 0) return

    !    if (ex) then
    version = 1
    status = pgs_io_gen_openf ( thresh_id, PGSd_IO_Gen_RSeqFrm, &
         0,lun, version)
    call tell_log(1,'read_thresholds: trying to open threshold file')
    call tell_log(1,'status, lun')
    write(logmsg,"(I6,I5)") status, lun
    call tell_log(1,logmsg)
    if(status.ne.0) then
      call tell_error(tell_io_open_error, &
           "read_thresholds: error opening threshold table file", errstat)
      return
    endif
    !    else
    !      filename=trim(input_data_path)//trim(thresh_file)
    !      print *,'default output file name = ',filename
    !      open(lun,file=filename,status='old', form="formatted", iostat=rc)
    !      if (rc /= 0) then
    !        if (iprt > 0) print *,'read_thresholds: error reading ',filename
    !        stop
    !      endif
    !    endif
    read(lun, *, err=100) text
    read(lun, *, err=100) npixs, nscanpos
    call tell_log(1,text)
    call tell_log(1,'read_thresholds: npixs,nscanpos')
    write(logmsg,"(2I4)") npixs,nscanpos
    call tell_log(1,logmsg)

    !allocate the threshold table, but fill with default
    !constant value in case file not found
    !===================================================
    if (.not. allocated(thresholds)) then
      allocate(thresholds(npixs,nscanpos), npixels(npixs), stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "read_thresholds: failed to allocate memory", &
             errstat)
        return
      endif
      thresholds=stddev_thresh
    endif
    read(lun, *, err=100) npixels
    do i=1, npixs
      do j=1, nscanpos
        read(lun, *, err=100) l, m, thresholds(i, j) 
!        if (iprt > 3) print *, l, m, i, j, thresholds(i, j)
      enddo
    enddo

    if (ex) then
      status = pgs_io_gen_closef (lun)
    else
      close(lun)
    endif

    return

100 status = 1
    call tell_error(tell_io_read_error, &
         'read_thresholds: error reading file', errstat)
    return

  end subroutine read_thresholds

end module m_read_thresholds
