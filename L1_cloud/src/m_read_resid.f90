module m_read_resid

contains

  subroutine read_resids (errstat) 

    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !IROUTINE:  read_resids
    !                
    !
    ! !INTERFACE:
    !
    ! !USES: read_resids reads precomputed resids spectra needed for cloud
    !               pressures
    !
    use m_vars, ONLY: nwav, nscanpos, resid_spec, iprt, read_he4
    use m_LUN_set
    use m_pgs_include
    use tell_module
    Implicit NONE

    ! !REVISION HISTORY:
    !
    !  25aug04   Joiner      original fortran 90
    !  26mar15   O'Sullivan  updated for TEMPO
    !
    !EOP
    !-------------------------------------------------------------------------
    integer, intent(inout) :: errstat

    !local variables
    integer :: i                       
    character(len=100) :: text
    integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
    integer :: status, version, ierr, lun
    integer :: pgs_met_getPCAttr_i, pgs_pc_getconfigdata
    integer :: OrbitNumber, ThreshOrbitNumber
    character(len=200) :: buf

    if (errstat /= 0) return

    version = 1
    !If reading from he4 input, get orbit number
    if (read_he4) then
      status = pgs_met_getPCAttr_i(L1B_LUN, version , "CoreMetadata.0", &
           "OrbitNumber.1",OrbitNumber)
      IF(status /= 0 ) THEN
        ierr = OMI_SMF_setmsg( status, &
             "Warning: Could not get orbit number", "read_resids", 1 )
        OrbitNumber=1
      ENDIF
    else !use the orbit number in the PCF file
      status = pgs_pc_getconfigdata(OrbNum_LUN,buf)
      read(buf,*) OrbitNumber
      if (status /= 0) then
        ierr = OMI_SMF_setmsg( status, &
             "Warning: Could not get orbit number from PCF", "read_resids", 1 )
      endif
      if (iprt >= 1) then
        print *,'read_resids: using orbit number from PCF ',OrbitNumber
      endif
    endif




    status = pgs_pc_getconfigdata(ThreshOrbNum_LUN,buf)
    IF(status /= 0 ) THEN
      ierr = OMI_SMF_setmsg( status, &
           "Error getting ThreshOrbNum", "read_resids", 1 )
      ThreshOrbitNumber=999999
    ELSE
      read(buf,*) ThreshOrbitNumber
    ENDIF

    version = 1
    if(OrbitNumber .le. ThreshOrbitNumber) then 
      status = pgs_io_gen_openf (resid_id_early, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version) 
    else 
      status = pgs_io_gen_openf (resid_id_late, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
    endif

    if (iprt > 0) then
      print *,'read_resids: trying to open resid file ',status, lun
    endif

    if(status.ne.0) then
      ierr = OMI_SMF_setmsg( status, &
           "PGE aborting", "read_resids", 1 ) 
      errstat = -1
      return
    else
      read(lun, *, err=100) text
      read(lun, *, err=100) nscanpos, nwav
      if (iprt >= 1) then
        print *,text
        print *,'read_resids: nscanpos, nwav'
        print *,nscanpos, nwav
      endif

      !allocate the resid table, but fill with default
      !constant value in case file not found
      !===================================================
      if (.not. allocated(resid_spec)) then
        allocate(resid_spec(nwav,nscanpos), stat=errstat)
        if (errstat /= 0) then
          call tell_error (tell_malloc_error, &
               "read_resids: failed to allocate memory", &
               errstat)
          return
        endif
      endif
      do i=1, nwav
        read(lun, *, err=100) resid_spec(i,1:nscanpos) 
        if (iprt >= 2) print *, 'read_resids: ',i, resid_spec(i,1:nscanpos)
      enddo

      status = pgs_io_gen_closef (lun)
      if (iprt > 2) print *,'read_resids: done reading file'
    endif

    return

100 status = 1
    if (iprt > 0) print *,'read_resids: error reading file'
    ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
         "Error reading resid table, PGE aborting", "read_resids", 1 )
    errstat = -1
    return

  end subroutine read_resids

  subroutine read_o3 (errstat)

    use m_vars, ONLY: wave_o3, xsect_o3, iprt
    use m_LUN_set
    use m_pgs_include
    use tell_module
    Implicit NONE

    ! !REVISION HISTORY:
    !
    !  25aug04   Joiner      original fortran 90
    !  26mar15   O'Sullivan  update for TEMPO
    !
    !EOP
    !-------------------------------------------------------------------------
    integer, intent(inout) :: errstat
    !local variables
    integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
    integer :: status, version, ierr, lun, nwav_o3

    if (errstat /= 0) return

    version = 1
    status = pgs_io_gen_openf ( o3_id, PGSd_IO_Gen_RSeqFrm, &
         0,lun, version)
    if (iprt > 0) then
      print *,'read_o3: trying to open o3 file ',status, lun
    endif
    if(status.ne.0) then
      ierr = OMI_SMF_setmsg( status, &
           "PGE aborting", "read_o3", 1 ) 
      errstat = -1
      return
    else
      read(lun, *, err=100) nwav_o3
      if (iprt >= 1) then
        print *,'read_o3: nwav_o3'
        print *, nwav_o3
      endif

      !allocate the resid table, but fill with default
      !constant value in case file not found
      !===================================================
      if (.not. allocated(wave_o3)) then
        allocate(wave_o3(nwav_o3), xsect_o3(nwav_o3), stat=errstat)
        if (errstat /= 0) then
          call tell_error (tell_malloc_error, &
               "read_resids: failed to allocate memory", &
               errstat)
          return
        endif
      endif

      read(lun, *, err=100) wave_o3
      read(lun, *, err=100) xsect_o3
      if (iprt >= 2) print *, 'read_o3: ', wave_o3
      if (iprt >= 2) print *, 'read_o3: ', xsect_o3

      status = pgs_io_gen_closef (lun)
      if (iprt > 2) print *,'read_o3: done reading file'
    endif

    return

100 status = 1
    if (iprt > 0) print *,'read_o3: error reading file'
    errstat = -1
    call tell_error (tell_io_read_error, &
         "read_o3: failed to read O3 cross-section file", &
         errstat)
    return


  end subroutine read_o3

end module m_read_resid
