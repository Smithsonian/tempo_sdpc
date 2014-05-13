module m_read_ocean_table

private
public read_ocean_table

contains

subroutine read_ocean_table(rc)

use m_vars, ONLY: nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl, &
 oc_perms, nwave2, &
 wgrid_out_oc, wgrid_oc, oc_table, w_grid, &
 theta_oc,scan_oc,phi_oc,ocrefl,chl,iprt
use m_interpol
use m_LUN_set
use m_pgs_include

implicit none

      integer, intent(out)         :: rc        ! Error return code:

!character(len=255) :: fn
!integer            :: i,j,k,l,m
integer            :: i,j,m
real (KIND=8), allocatable, dimension(:) :: oc_perms2
! ocean Raman correction coefficient
real (KIND=8), parameter :: coef = 1.5

!integer, parameter :: lun=10
!lun is an output of pgs_io_gen_openf so cannot be parameter!
integer :: lun=10
integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
integer :: status, version, err_code
!include 'PGS_IO.f'
!include 'PGS_IO_1.f'
!include 'PGS_OMI_1900.f'
!include 'PGS_SMF.f'

version = 1
rc=0
status = pgs_io_gen_openf ( oc_ram_id, PGSd_IO_Gen_RSeqUnf, &
        0,lun, version)
if(status.ne.0) then
  err_code=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening ocean table file', &
  'read_ocean_table, module m_read_ocean_table',2)
  rc=1
  call exit(1)
endif

read(lun) nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl
if (iprt >= 1) then
  print *, 'read_ocean_table: nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl'
  print *, nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl
endif
allocate(oc_perms(nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl))
allocate(oc_perms2(nwave_oc))
allocate(oc_table(nthet_oc,nscan_oc,nchl,nwave2))
allocate(wgrid_out_oc(nwave_oc),theta_oc(nthet_oc), &
   scan_oc(nscan_oc),phi_oc(nphi_oc),wgrid_oc(nwave_oc))
allocate(ocrefl(nocrefl),chl(nchl))

read(lun) oc_perms,wgrid_out_oc,theta_oc,scan_oc,phi_oc,ocrefl,chl
wgrid_oc=wgrid_out_oc
 do i=1,nthet_oc
  do j=1,nscan_oc
     do m=1,nchl
       ! use 10% reflectivity
       oc_perms2(:)=oc_perms(:,i,j,1,2,m)
       oc_table(i,j,m,:)= coef * interpol(oc_perms2(:),wgrid_oc,w_grid)
     enddo
  enddo
 enddo

if (iprt >= 6) then
  print *,'ocean_table wavelengths'
  write(6,100) wgrid_out_oc
  print *,'ocean_table_new (%)'
  write(6,100) oc_table(1,1,1,:)*100
endif
deallocate(oc_perms2)
deallocate(oc_perms)
deallocate(wgrid_out_oc)
status = pgs_io_gen_closef (lun)

100 format (6f12.3) 
end subroutine read_ocean_table   
end module m_read_ocean_table   

