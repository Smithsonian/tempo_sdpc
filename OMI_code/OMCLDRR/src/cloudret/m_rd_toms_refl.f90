module m_rd_toms_refl

private
public rd_toms_refl

contains

subroutine rd_toms_refl( )   

use m_vars, ONLY: done_read_refl, iprt, lat, lon, toms_refl, ref_nmon, &
  iLine, nXtrack, ref_nlat, ref_nlon, ref_clr, ref_lats, ref_lons, month, &
  ler_sz, ler_th, ler_ph, ler354
use m_LUN_set
implicit NONE          

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!!BOP
!
! !ROUTINE: rd_toms_refl
!
! !DESCRIPTION: get TOMS reflectivity climatology
!
! !CALLING SEQUENCE: call rd_toms_refl (lat, lon, terr_pres)
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
!  01Aug07   Joiner     Original code 
!!EOP
!-------------------------------------------------------------------------

!local variables
!================
!integer, parameter :: lun=2 
!lun cannot be a parameter since it's an output of pgs_io_gen_openf
integer :: lun=2 

integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
integer :: status,ierr, version=1

include 'PGS_IO.f'
include 'PGS_IO_1.f'
include 'PGS_OMI_1900.f'
include 'PGS_SMF.f'

integer                    :: ipts, i, j
real (KIND=8)                       :: lont, latt
real (KIND=8)                       :: deltlat, deltlon
real (KIND=8)                       :: startlat, startlon
!integer                    :: iret
character(len=100)         :: txt
   
!=======================
!read terrain data set
!=======================
if (.not. done_read_refl) then
  status = pgs_io_gen_openf ( refl_id, PGSd_IO_Gen_RSeqFrm, &
        0,lun, version)
  if(status.ne.0) then
    ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening reflectivity file', &
    'rd_toms_refl, module m_rd_toms_refl',2)
    call exit(1)
  endif
  if (iprt > 0) print *,'rd_toms_refl: opening reflectivity file, status :',status
  read (lun,*,err=200)  txt
  if (iprt > 0) print *, txt
  read (lun,*,err=200)  txt
  if (iprt > 0) print *, txt
  read (lun,*,err=200)  ref_nlon, ref_nlat, ref_nmon
  if (iprt > 0) print *, ref_nlon, ref_nlat, ref_nmon
  allocate(ref_lats(ref_nlat))
  allocate(ref_lons(ref_nlon))
  allocate(toms_refl(ref_nmon,ref_nlon,ref_nlat))
  read (lun,*,err=200)  ref_lons
  read (lun,*,err=200)  ref_lats
  if (iprt > 1) print *, ref_lons  
  if (iprt > 1) print *, ref_lats
  read (lun,*,err=200)  toms_refl
  status = pgs_io_gen_closef (lun)
  if (iprt > 0) print *,'rd_toms_refl: closing reflectivity file, status :',status
  status = pgs_io_gen_openf ( ler354_id, PGSd_IO_Gen_RSeqFrm, &
        0,lun, version)
  if(status.ne.0) then
    ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening ler354_cox_munk file', &
    'rd_toms_refl, module m_rd_toms_refl',2)
    call exit(1)
  endif
  read (lun,*,err=201) ler_sz
  read (lun,*,err=201) ler_th
  read (lun,*,err=201) ler_ph
  read (lun,*,err=201) ler354
  status = pgs_io_gen_closef (lun)
  if (iprt > 0) print *,'rd_toms_refl: closing ler354_cox_munk file, status :',status
  done_read_refl=.true.
  deltlat=ref_nlat/180.
  deltlon=ref_nlon/360.
  startlat=ref_lats(1)
  startlon=ref_lons(1)
endif ! not done reading

!if (iprt > 0) print *,'deltlat, lon ',deltlat, deltlon
   
do ipts=1, nXtrack
 latt=lat(ipts,iLine)
 i=anint((startlat-latt)/deltlat)+1   
 if (i <= 0) i=1
 if (i >= ref_nlat+1) i=ref_nlat
 lont=lon(ipts,iLine)   
 if(lont > 180.) then 
   lont=lont-360 
 endif
 j=anint((lont-startlon)/deltlon)+1   
 if(j <= 0) j=1   
 if(j >= ref_nlon+1) j=1   
 !print *, 'rd_toms_refl : ',ipts, lat(ipts,iLine), lon(ipts,iLine), &
 !   i, j
 !print *, toms_refl(imon,j,i)
 ref_clr(ipts-1,iLine)=toms_refl(month,j,i)/100.
enddo   ! ipts
   
return
200 print *, 'rd_toms_refl: error reading reflectivity file'
toms_refl=0.
201 print *, 'rd_toms_refl: error reading Cox-Munk LER file'
ler354=0.

end subroutine rd_toms_refl   

end module m_rd_toms_refl
