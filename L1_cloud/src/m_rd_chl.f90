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

    use m_vars, ONLY: done_read_chl, chl2d, iprt, lat, lon, chlcl, iLine, &
         nXtrack
    use m_LUN_set
    use m_pgs_include

    implicit none          

    integer, intent(inout) :: errstat
    !local variables
    integer                    :: lun=3
    integer                    :: ip, i, j
    integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
    integer :: status, ierr,version

    version = 1

    if (errstat /= 0) return

    !read in chlorophyll data set
    !============================
    if (.not. done_read_chl) then
      status = pgs_io_gen_openf ( chl_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening clorophyll file', &
             'rd_chl, module m_rd_chl',2)
        errstat = -1
        return
      endif

      if (iprt > 0) print *,'rd_chl: opening chl file, status ',status
      read (lun,*,err=200)  chl2d
      status = pgs_io_gen_closef (lun)
      if (iprt > 0) print *,'rd_chl: closing chl file, status ',status
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
200 print *, 'rd_chl: error reading chlorophyll file'
    chl2d=0.

  end subroutine rd_chl

end module m_rd_chl
