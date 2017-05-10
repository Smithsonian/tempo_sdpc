!>Read terrain pressure climatology reference file
!
!-------------------------------------------------------------------------
!
! !ROUTINE: rd_terr
!
!> DESCRIPTION: get terrain pressure (0.5X0.5 deg resolution)
!
! !CALLING SEQUENCE: call rd_terr (lat, lon, terr_pres)
!
! !INPUT PARAMETERS: 
!> @param lat[in]      latitude
!> @param lon[in]      longitude
!
! !OUTPUT PARAMETERS:  
!> @param terr_pres[out]       terrain pressure
!
! !SEE ALSO: 
!
! !REVISION HISTORY:
!
!> @author  01Oct96   Joiner       Original code 
!> @author  19Mar02   Vasilkov     To read filenames from PCF!
!> @author  26Mar15   O'Sullivan   Update for TEMPO
!-------------------------------------------------------------------------
module m_rd_terr

  private
  public rd_terr

contains

  subroutine rd_terr (errstat)   

    use m_vars, ONLY: p_terr, done_read_terr, lat, lon, ps, &
         iLine, nXtrack
    use m_LUN_set
    use m_pgs_include
    use tell_module
    implicit NONE          

    integer, intent (inout) :: errstat

    !local variables
    !================
    integer :: lun=2 
    integer :: pgs_io_gen_openf, pgs_io_gen_closef
    integer :: status, version
    integer :: ipts, i, j
    real (KIND=8) :: lont, latt
    character (len=128) :: logmsg

    version = 1

    if (errstat /= 0) return

    !=======================
    !read terrain data set
    !=======================
    if (.not. done_read_terr) then
      status = pgs_io_gen_openf ( terr_prs_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        call tell_error(tell_io_open_error, &
             "rd_terr: error opening terrain pressure file", errstat)
        return
      endif
      write(logmsg,"(A39, I4)") 'rd_terr: opening terrain file, status :',&
           status
      call tell_log(1,logmsg)
      read (lun,*,err=200)  p_terr
      status = pgs_io_gen_closef (lun)
      write(logmsg,"(A39, I4)") 'rd_terr: closing terrain file, status :',&
           status
      call tell_log(1,logmsg)
      done_read_terr=.true.
    endif

    do ipts=1, nXtrack
      ! Check lat/lon within bounds, ps=0 if not
      if (lat(ipts,iLine) > 90.0 .or. lat(ipts,iLine) < -90.0 .or. &
           lon(ipts,iLine) > 180.0 .or. lon(ipts,iLine) < -180.0) then
        ps(ipts-1,iLine)=0.0d0
      else
        latt=-lat(ipts,iLine)
        i=int(anint((latt+90.0)/0.5),kind(i))+1
        if (i == 361) i=360
        if(lon(ipts,iLine) < 0) then 
          lont=lon(ipts,iLine)+360 
        else 
          lont=lon(ipts,iLine)   
        endif
        j=int(anint((lont)/0.5),kind(j))+1
        if(j == 721) j=1   
        ps(ipts-1,iLine)=p_terr(i,j)/1013.0   
      endif
    enddo   ! ipts

    return

    !FIXME - should this error cause the code to abort?
200 call tell_error(tell_io_read_error, &
         'rd_terr: error reading terrain file', errstat)
    p_terr=0.d0

  end subroutine rd_terr

end module m_rd_terr
