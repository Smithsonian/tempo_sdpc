module m_rd_chl

private
public rd_chl

contains

subroutine rd_chl( )   

use m_vars, ONLY: done_read_chl, chl2d, iprt, lat, lon, chlcl, iLine, nXtrack, form
use m_LUN_set

!*************************************************************************
!        AUTHOR:  Joanna Joiner - original code, 
!                modified by Alexander P. Vasilkov
!
!         FUNCTION: get chlorophyll concentration (0.5X0.5 deg resolution)
!
!         CALLING SEQUENCE: call rd_chl (lat, lon, chl_out)
!
!         INPUT:         lat            latitude (-90.0,90.0)
!                      lon     longitude (-180-180) -> (0.0,360.0)
!
!         OUTPUT:  chl_out        chlorophyll concentration
!
!        HISTORY: Developed Nov. 15, 2001
!
!***************************************************************************
implicit none          
!real (KIND=8), dimension(:), intent(in)  :: lat, lon
!real (KIND=8), dimension(:), intent(out) :: chl_out

!local variables

!integer, parameter :: lun=3
! lun cannot be a parameter as it's an output of pgs_io_gen_openf
integer :: lun=3
integer                    :: ip, i, j
integer                    :: iret

!integer, parameter :: chl_id = 510004
integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
integer :: status, ierr,version=1
include 'PGS_IO.f'
include 'PGS_IO_1.f'
include 'PGS_OMI_1900.f'
include 'PGS_SMF.f'

!*****************************************************************
   
!read in chlorophyll data set
!============================
if (.not. done_read_chl) then
  status = pgs_io_gen_openf ( chl_id, PGSd_IO_Gen_RSeqFrm, &
        0,lun, version)
!  IF( status .NE. OMI_S_SUCCESS ) THEN
if(status.ne.0) then
  ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening clorophyll file', &
  'rd_chl, module m_rd_chl',2)
  call exit(1)
endif

  if (iprt > 0) print *,'rd_chl: opening chl file, status ',status
  read (lun,*,err=200)  chl2d
  status = pgs_io_gen_closef (lun)
  if (iprt > 0) print *,'rd_chl: closing chl file, status ',status
  done_read_chl=.true.
endif
   
do ip=1, nXtrack !size(lat,dim=2)
 i=anint((lat(ip,iLine)+90.0)/0.5)+1  
 if(form==5) j=anint((lon(ip,iLine)+180.0)/0.5)+1
 if(form==2) j=anint((lon(ip,iLine)+0.0)/0.5)+1
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
