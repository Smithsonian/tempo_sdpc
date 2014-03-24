module m_rd_terr

private
public rd_terr

contains

subroutine rd_terr( )   

use m_vars, ONLY: p_terr, done_read_terr, iprt, lat, lon, ps, iLine, nXtrack
use m_LUN_set
implicit NONE          

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!!BOP
!
! !ROUTINE: rd_terr
!
! !DESCRIPTION: get terrain pressure (0.5X0.5 deg resolution)
!
! !CALLING SEQUENCE: call rd_terr (lat, lon, terr_pres)
!
! !INPUT PARAMETERS: 
!real (KIND=8), dimension(:), intent(in)  :: lat, lon
!                     lat     : latitude
!                     lon     : longitude
!
! !OUTPUT PARAMETERS:  
!real (KIND=8), dimension(:), intent(out) :: terr_pres
!                     terr_pres      : terrain pressure
!
! !SEE ALSO: 
!
! !REVISION HISTORY:
!
!  01Oct96   Joiner     Original code 
!  19Mar02   Vasilkov   To read filenames from PCF!
!!EOP
!-------------------------------------------------------------------------

!local variables
!================
integer, parameter :: lun=2 

!integer, parameter :: terr_prs_id = 510003
integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
integer :: status,ierr, version=1
include 'PGS_IO.f'
include 'PGS_IO_1.f'
include 'PGS_OMI_1900.f'
include 'PGS_SMF.f'

integer                    :: ipts, i, j
real (KIND=8)                       :: lont, latt
integer                    :: iret
   
!=======================
!read terrain data set
!=======================
if (.not. done_read_terr) then
  status = pgs_io_gen_openf ( terr_prs_id, PGSd_IO_Gen_RSeqFrm, &
        0,lun, version)
if(status.ne.0) then
  ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening terrain pressure file', &
  'rd_terr, module m_rd_terr',2)
  call exit(1)
endif
  if (iprt > 0) print *,'rd_terr: opening terrain file, status :',status
  read (lun,*,err=200)  p_terr
  status = pgs_io_gen_closef (lun)
  if (iprt > 0) print *,'rd_terr: closing terrain file, status :',status
  done_read_terr=.true.
endif
   
do ipts=1, nXtrack
 latt=-lat(ipts,iLine)
 i=anint((latt+90.0)/0.5)+1   
 if (i == 361) i=360
  if(lon(ipts,iLine) < 0) then 
    lont=lon(ipts,iLine)+360 
  else 
    lont=lon(ipts,iLine)   
  endif
 j=anint((lont)/0.5)+1   
 if(j == 721) j=1   
!print *, 'rd_terr : ',ipts, lat(ipts,iLine), lon(ipts,iLine), &
!   i, j!, ps(ipts)
 ps(ipts-1,iLine)=p_terr(i,j)/1013.0   
enddo   ! ipts
   
return
200 print *, 'rd_terr: error reading terrain file'
p_terr=0.

end subroutine rd_terr   

end module m_rd_terr
