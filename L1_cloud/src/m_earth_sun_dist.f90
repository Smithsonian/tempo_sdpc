module m_earth_sun_dist

contains

  subroutine EarthSunDist(fn,dist1,dist2) 

    use m_LUN_set
    use m_swathnames
    use m_strpos
    use m_pgs_include
    use tell_module

    implicit none

    ! INCLUDE 'PGS_TD_3.f'
    ! INCLUDE 'PGS_MET.f'
    ! INCLUDE 'PGS_SMF.f'
    ! INCLUDE 'PGS_MET_13.f'
    ! INCLUDE 'PGS_OMI_1900.f'
    ! INCLUDE 'PGS_OMCLDRR_52251.f'

    character(len=*), intent(in) :: fn
    real(kind=4), intent(out) :: dist1, dist2
    character(len=200)  :: swn, filenm, logmsg
    integer (kind=4), parameter :: zero = 0
    integer (kind=4) :: OMI_SMF_setmsg, swopen, swattach, swrdattr, &
         swdetach, swclose
    integer, parameter :: DFACC_READ = 1
    integer (kind=4) :: pgs_pc_getreference 
    integer (kind=4) :: swfid, swid, ierr
    integer (kind=4) :: status, version

    version = 1
    status = pgs_pc_getreference( L1B_LUN, version, filenm )
    IF( status .NE. PGS_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_W_MET, &
           "get L1B name failed ", "EarthSunDist", zero )
    END IF

    !    vis  = strpos (filenm, 'BRVG') > 0
    !    visz = strpos (filenm, 'BRVZ') > 0
    uvsz = strpos (filenm, 'BRUZ') > 0

    ! if (visz) then
    !   swn = sunvisswathz
    ! else if (uvsz) then
    !   swn = sunuv2swathz
    ! else if (vis) then

    !    if (vis .or. visz) then
    !      swn = sunvisswath
    !    else
    swn = sunuv2swath
    !    endif

    write(logmsg,*) 'earth_sun_dist: ',trim(fn),' ',trim(swn)
    call tell_log(2,logmsg)

    !! open the  swath file
    swfid = swopen( fn, DFACC_READ )
    IF( swfid < zero ) THEN
      status = OMI_E_FAILURE
      ierr = OMI_SMF_setmsg( OMI_E_FILE_OPEN, fn, "EarthSunDist", zero )
      RETURN
    ENDIF

    !! attach to the swath
    swid = swattach( swfid, swn )
    IF( swid < zero ) THEN
      status = OMI_E_FAILURE
      ierr = OMI_SMF_setmsg( OMI_E_SWATH_ATTACH, swn, "EarthSunDist", zero )
      RETURN
    ENDIF

    status = swrdattr( swid, "EarthSunDistance", dist1 )
    If( status < zero ) THEN
      status = OMI_E_FAILURE
      ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
           "get EarthSundistance failed", "EarthSunDist", zero )
      RETURN
    ENDIF
    ierr = swdetach( swid )
    ierr = swclose( swfid )

    !    if (visz) then
    !      swn = visswathz
    !    else if (uvsz) then
    !      swn = uv2swathz
    !    else if (vis) then
    !      swn = visswath
    !    else
    !      swn = uv2swath
    !    endif

    if (uvsz) then
      swn = uv2swathz
    else
      swn = uv2swath
    endif

    write(logmsg,*) 'earth_sun_dist: ',trim(filenm),' ',trim(swn)
    call tell_log(2,logmsg)

    !! open the  swath file
    swfid = swopen( filenm, DFACC_READ )
    IF( swfid < zero ) THEN
      status = OMI_E_FAILURE
      ierr = OMI_SMF_setmsg( OMI_E_FILE_OPEN, filenm, "EarthSunDist", zero )
      RETURN
    ENDIF

    !! attach to the swath
    swid = swattach( swfid, swn )
    IF( swid < zero ) THEN
      status = OMI_E_FAILURE
      ierr = OMI_SMF_setmsg( OMI_E_SWATH_ATTACH, swn, "EarthSunDist", zero )
      RETURN
    ENDIF

    status = swrdattr( swid, "EarthSunDistance", dist2 )
    If( status < zero ) THEN
      status = OMI_E_FAILURE
      ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
           "get EarthSundistance failed", "EarthSunDist", zero )
    ENDIF
    ierr = swdetach( swid )
    ierr = swclose( swfid )
    write(logmsg,*) 'earth_sun_dist: distances ', dist1, dist2
    call tell_log(2,logmsg)
    write(logmsg,*) trim(fn)
    call tell_log(2,logmsg)
    write(logmsg,*) trim(filenm)
    call tell_log(2,logmsg)
    write(logmsg,*) trim(swn)
    call tell_log(2,logmsg)

  end subroutine EarthSunDist

end module m_earth_sun_dist
