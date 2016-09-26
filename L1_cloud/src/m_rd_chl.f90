!>Read chlorophyll concentration climatology reference file
!
!--------------------------------------------------------------------------
!> @author Joanna Joiner - original code, 
!> @author modified by Alexander P. Vasilkov
!> @author updated for TEMPO by E. O'Sullivan
!
!         FUNCTION: get chlorophyll concentration (0.5X0.5 deg resolution)
!
!         CALLING SEQUENCE: call rd_chl (lat, lon, chl_out)
!
!         INPUT:         
!> @param lat[in]            latitude (-90.0,90.0)
!> @param lon[in]     longitude (-180-180) -> (0.0,360.0)
!
!         OUTPUT:  
!> @param chl_out[out]        chlorophyll concentration
!
!        HISTORY: Developed Nov. 15, 2001
!
!--------------------------------------------------------------------------
module m_rd_chl

  private
  public rd_chl

contains

  subroutine rd_chl (errstat)   

    use m_vars, ONLY: done_read_chl, chl2d, lat, lon, chlcl, iLine, &
         nXtrack
    use m_LUN_set
    use m_pgs_include
    use tell_module

    implicit none          

    integer, intent(inout) :: errstat
    !local variables
    integer                    :: lun=3
    integer                    :: ip, i, j
    integer :: pgs_io_gen_openf, pgs_io_gen_closef
    integer :: status, ierr,version
    character (len=128) :: logmsg

    version = 1

    if (errstat /= 0) return

    !read in chlorophyll data set
    !============================
    if (.not. done_read_chl) then
      status = pgs_io_gen_openf ( chl_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        errstat = -1
        call tell_error (tell_io_open_error, &
             "rd_chl: error opening clorophyll file", errstat)
        return
      endif

      write(logmsg,"(A33,I4)") 'rd_chl: opening chl file, status ',status
      call tell_log(1,logmsg)
      read (lun,*,err=200)  chl2d
      status = pgs_io_gen_closef (lun)
      write(logmsg,"(A33,I4)") 'rd_chl: closing chl file, status ',status
      call tell_log(1,logmsg)
      done_read_chl=.true.
    endif

    do ip=1, nXtrack !size(lat,dim=2)
      i=int(anint((lat(ip,iLine)+90.0)/0.5, kind=4),kind(i))+1
      j=int(anint((lon(ip,iLine)+180.0)/0.5, kind=4),kind(j))+1
      if(j <= 0) j=1
      if(j == 721) j=1   
      chlcl(ip-1)=chl2d(j,i)   
      ! print *, ip, lat(ip), lon(ip), i, j, chlcl(ip)
    enddo   ! ip

    return

    !FIXME - should this failure cause the code to abort?
200 call tell_error(tell_io_read_error, &
         'rd_chl: error reading chlorophyll file', errstat)
    chl2d=0.d0

  end subroutine rd_chl

end module m_rd_chl
