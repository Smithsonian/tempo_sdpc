module m_write_HDFEOS_attr

contains

  subroutine write_HDFEOS_attr(outfile, swathname, ierr)
    !***********************************
    use m_vars, ONLY: status,version
    use m_swathnames
    use m_HDFEOS_attr
    use m_LUN_set
    use m_strpos
    use m_pgs_include
    use tell_module
    !**********************************
    implicit none
    !
    ! !ROUTINE:  write_HDFEOS_attr
    ! 
    ! !DESCRIPTION: writes HDF-EOS attributes to output file
    !
    ! !CALLING SEQUENCE: 
    !
    !        call write_HDFEOS_attr
    !     
    ! !INPUT PARAMETERS:   
    !
    ! !OUTPUT PARAMETERS:  
    !
    ! !REVISION HISTORY: 
    !
    !  28Apr03   Vasilkov     original fortran 90
    !  09Jun03   Vasilkov     excluded writing time
    !  09Sep07   Joiner       mods for zoom data
    !-------------------------------------------------------------------------
    !
    !inputs
    !------
    character(len=*), intent(in) :: outfile, swathname

    !ouputs
    !------
    integer, intent(out) :: ierr

    !include
    ! INCLUDE 'PGS_SMF.f'
    ! INCLUDE 'PGS_MET_13.f'
    ! INCLUDE 'PGS_OMI_1900.f'
    ! INCLUDE 'PGS_OMCLDRR_52251.f'

    ! Declare the HDF-EOS file and swath identification numbers, and
    ! the status of the HDF-EOS functions calls.
    !-----------------------------------------------------------------
    integer (kind = 4) swfid, swid

    ! Declare the HDF-EOS functions.
    !-------------------------------
    integer, parameter :: HE5_ACC_RDWR=100
    !  integer, parameter :: HE5_ACC_TRUNC=102
    integer (kind = 4) :: he5_swopen, he5_swattach, he5_swdetach, he5_swclose
    integer (kind = 4) :: pgs_pc_getreference, pgs_met_getPCAttr_s
    !  integer (kind = 4) :: pgs_pc_getuniversalref
    integer (kind = 4) :: OMI_SMF_setmsg 
    integer (kind = 4) :: year, month, day
    !  integer (kind = 4) :: NumTimes, NumTimesSmallPix
    !Local variables
    !---------------
    character(len=255) :: rad_flnm, str_value, swn, logmsg
    character(len=1) :: buf
    !**************************************************************

    ! open the file
    !----------------------
    swfid = he5_swopen (outfile, HE5_ACC_RDWR)

    !Get L1B radiance filename
    version = 1
    status = pgs_pc_getreference(L1B_LUN,version,rad_flnm)
    IF(status /= 0 ) THEN
      status = OMI_SMF_setmsg(OMCLDRR_W_MET, "get L1B radiance filename failed", &
           "write_HDFEOS_attr", 0 )
      ierr=status
    ENDIF

    !Write the global attribute
    !-------------------------------------------------
    version = 1
    status = pgs_met_getPCAttr_s(L1B_LUN, version , "CoreMetadata.0", &
         "RANGEBEGINNINGDATE",str_value)
    IF(status /= 0 ) THEN
      status = OMI_SMF_setmsg(OMCLDRR_F_FAILURE, "get RANGEBEGINNINGDATE failed", &
           "write_HDFEOS_attr", 0 )
      ierr=status
    ENDIF
    read(str_value, '(i4,a1,i2,a1,i2)') year,buf,month,buf,day

    status = CLDRR_WriteGlobalAttr(swfid, year, month, day)  

    ! Attach to the swath
    !----------------------
    swid = he5_swattach(swfid, swathname)

    !Write swath attribute
    !-------------------------------------------------
    if (visz) then
      swn = visswathz
    else if (uvsz) then
      swn = uv2swathz
    else if (vis) then
      swn = visswath
    else
      swn = uv2swath
    endif
    status = CLDRR_WriteSwathAttr(swid, rad_flnm, swn)

    ! Detach from the swath interface.
    !-------------------------------------------------
    status = he5_swdetach (swid)
    write(logmsg,"(A34, I12)") 'write_HDFEOS_attr: detached swath ',status
    call tell_log(2,logmsg)

    ! Close the OMI Level2  HDF-EOS output file.
    !----------------------------------------------
    status = he5_swclose (swfid)
    write(logmsg,"(A32, I12)") 'write_HDFEOS_attr: closed swath ',status
    call tell_log(2,logmsg)
  end subroutine write_HDFEOS_attr

end module m_write_HDFEOS_attr
