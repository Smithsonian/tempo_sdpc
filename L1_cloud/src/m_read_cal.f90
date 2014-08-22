module m_read_cal

contains
  !-------------------------------------------------------------------------
  !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:  read_cals
  !                
  !
  ! !INTERFACE:
  !
  subroutine read_cals (rc)

    ! !USES:
    !  read_cals reads precomputed resids spectra needed for cloud
    !               pressures
    !
    use m_vars, ONLY: nscanpos, cal_fact, iprt, using_cal
    use m_LUN_set
    use m_pgs_include
    Implicit NONE

    ! !INPUT PARAMETERS:
    !

    ! !OUTPUT PARAMETERS:
    !

    integer, intent(out)         :: rc        ! Error return code:
    !  0   all is well
    !  1   files not found

    ! !DESCRIPTION: 
    !
    ! !REVISION HISTORY:
    !
    !  25aug04   Joiner     original fortran 90
    !
    !EOP
    !-------------------------------------------------------------------------
    !integer :: i, ii
    character(len=100) :: text

    !integer, parameter :: lun=11
    integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
    integer :: status, version, ierr, lun
    !include 'PGS_IO.f'
    !include 'PGS_IO_1.f'
    !include 'PGS_OMI_1900.f'
    !include 'PGS_SMF.f'

    version = 1
    rc=0
    status = pgs_io_gen_openf ( cal_id, PGSd_IO_Gen_RSeqFrm, &
         0,lun, version)
    if (iprt > 0) then
      print *,'read_cals: trying to open cal file ',status, lun
    endif
    if(status.ne.0) then
      ! ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening cal table file', &
      ! 'read_cals, module m_read_cals',1)
      ierr = OMI_SMF_setmsg( status, &
           "PGE aborting, exit code = 1", "read_cals", 1 ) 
      rc=1
      call exit(1)
      using_cal=.false.
    else
      read(lun, *, err=100) text
      read(lun, *, err=100) nscanpos
      if (iprt >= 1) then
        print *,text
        print *,'read_cals: nscanpos'
        print *,nscanpos
      endif

      !allocate the cal table, but fill with default
      !constant value in case file not found
      !===================================================
      if (.not. allocated(cal_fact)) then
        allocate(cal_fact(nscanpos))
      endif
      read(lun, *, err=100) cal_fact(1:nscanpos) 
      if (iprt >= 2) print *, 'read_cals: ', cal_fact(1:nscanpos)

      status = pgs_io_gen_closef (lun)
      if (iprt > 2) print *,'read_cals: done reading file'
    endif

    return

100 status = 1
    if (iprt > 0) print *,'read_cals: error reading file'
    !   ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
    !     "Error reading reflectivity calibration, PGE aborting, exit code = 1", "read_cals", 1 )
    rc=1
    call exit(1)

  end subroutine read_cals

end module m_read_cal
