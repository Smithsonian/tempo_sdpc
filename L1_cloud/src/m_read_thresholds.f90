module m_read_thresholds

contains
  !-------------------------------------------------------------------------
  !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:  read_thresholds
  !                
  !
  ! !INTERFACE:
  !
  subroutine read_thresholds (errstat)

    ! !USES:
    !  read_thresholds reads precomputed thresholds needed for cloud
    !               mask
    !
    use m_vars, ONLY: npixs, nscanpos, thresholds, npixels, ex, iprt, &
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

    ! !DESCRIPTION: 
    !
    ! !REVISION HISTORY:
    !
    !  25aug04   Joiner      original fortran 90
    !  26mar15   O'Sullivan  updating for TEMPO
    !
    !EOP
    !-------------------------------------------------------------------------
    integer :: i, j, l, m
    character(len=100) :: text
    integer :: lun=10
    integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
    integer :: status, version, ierr

    if (errstat /= 0) return

    !    if (ex) then
    version = 1
    status = pgs_io_gen_openf ( thresh_id, PGSd_IO_Gen_RSeqFrm, &
         0,lun, version)
    if (iprt > 0) then
      print *,'read_thresholds: trying to open threshold file ',status, lun
    endif
    if(status.ne.0) then
      ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening threshold table file', &
           'read_thresholds, module m_read_thresholds',1)
      ierr = OMI_SMF_setmsg( status, &
           "PGE aborting", "read_thresholds", 1 ) 
      errstat = -1
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
    if (iprt >= 1) then
      print *,text
      print *,'read_thresholds: npixs,nscanpos'
      print *,npixs,nscanpos
    endif

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
        if (iprt > 3) print *, l, m, i, j, thresholds(i, j)
      enddo
    enddo

    if (ex) then
      status = pgs_io_gen_closef (lun)
    else
      close(lun)
    endif

    return

100 status = 1
    if (iprt > 0) print *,'read_thresholds: error reading file'
    ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
         "Error reading threshold table, PGE aborting", "read_thresholds", 1 )
    errstat = -1
    return

  end subroutine read_thresholds

end module m_read_thresholds
