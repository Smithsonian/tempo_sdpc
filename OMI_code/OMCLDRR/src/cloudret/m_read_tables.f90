module m_read_tables

contains
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  read_tables
!                
!
! !INTERFACE:
!
subroutine read_tables (rc)

! !USES:
!  read_tables reads precomputed tables needed for cloud
!               parameter retrievals                
!
use m_vars, ONLY: w_grid, nwave2, nwave, & 
  wmin, wmax, &
  ntheta, nscan, nphi, iprt, theta, scan, phi, &
  sflx, wgrid_out, npres, pres,  &
  k1bar, sba, nba, i01a, i0a, tra, nia, nra, out_path, ring_file_pre, &
  ring_file_suf, filename, ex, z1, z2

use m_find
use m_LUN_set

      Implicit NONE

! !INPUT PARAMETERS:
!

! !OUTPUT PARAMETERS:
!

      integer, intent(out)         :: rc        ! Error return code:
                                                !  0   all is well
                                                !  1   files not found

! !DESCRIPTION: uses piecewise parabolic
!
! !REVISION HISTORY:
!
!  05Jan01   Joiner     original fortran 90
!  12Aug02   Vasilkov   read filenames from PCF
!
!EOP
!-------------------------------------------------------------------------
integer :: i, j, k, l, m
logical :: first

integer, parameter :: lun=10
integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
integer :: status, version, ierr, OMCLDRR_F_FAILURE
include 'PGS_IO.f'
include 'PGS_IO_1.f'
include 'PGS_OMI_1900.f'
include 'PGS_SMF.f'

if (ex) then
version = 1
status = pgs_io_gen_openf ( ring_id, PGSd_IO_Gen_RSeqUnf, &
        0,lun, version)
if(status.ne.0) then
  ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening Ring table file', &
  'read_tables, module m_read_tables',1)
  ierr = OMI_SMF_setmsg( status, &
                            "PGE aborting, exit code = 1", "read_tables", 1 ) 
 call exit(1)
endif
else
 filename=trim(out_path)//trim(ring_file_pre)//'_all'//trim(ring_file_suf)
  print *,'default output file name = ',filename
 open(lun,file=filename,status='old', form="unformatted", iostat=rc)
 if (rc /= 0) then
   if (iprt > 0) print *,'read_tables: error reading ',filename
   stop
 endif
endif
 read(lun, err=100) nwave, ntheta, nscan, nphi, npres !nrefl
 if (iprt >= 1) then
   print *,'read_tables: nwave,ntheta,nscan,nphi,npres'
   print *,nwave,ntheta,nscan,nphi,npres
 endif
 if (.not. allocated(wgrid_out)) then
  allocate(sflx(nwave))
  allocate(wgrid_out(nwave))
  allocate(theta(ntheta))
  allocate(scan(nscan))
  allocate(phi(nphi))
  allocate(pres(npres))
  allocate(k1bar(nwave))
  allocate(sba(npres,nwave))
  allocate(nba(npres,nwave))
  allocate(i01a(npres,nphi,nscan,ntheta,nwave))
  allocate(i0a(npres,nscan,ntheta,nwave))
  allocate(z1(npres,nscan,ntheta,nwave))
  allocate(z2(npres,nscan,ntheta,nwave))
  allocate(tra(npres,nscan,ntheta,nwave))
  allocate(nia(npres,nphi,nscan,ntheta,nwave))
  allocate(nra(npres,nscan,ntheta,nwave))
 endif
 read(lun, err=100) theta, scan, phi, &
  pres, wgrid_out, sflx
do i=1, npres
 first=.true.
 do j=1, nphi
  do k=1, nscan
   do l=1, ntheta
    do m=1, nwave
     if (first) then
       read(lun, err=100) k1bar(m), nba(i,m), sba(i,m) ! hen using solar flux, only depends on wavelength
     endif
     read(lun, err=100) i01a(i,j,k,l,m),nia(i,j,k,l,m),nra(i,k,l,m),tra(i,k,l,m),i0a(i,k,l,m), z1(i,k,l,m), z2(i,k,l,m)
    enddo
    first=.false.
   enddo
  enddo
 enddo
enddo

i0a=log(i0a)

if (ex) then
  status = pgs_io_gen_closef (lun)
else
  close(lun)
endif

nwave2=nwave
allocate(w_grid(nwave))
w_grid=wgrid_out
! ierr = OMI_SMF_setmsg( OMCLDRR_S_SUCCESS, &
!      "Ring Table", "read_table", 1 )
return

100 status = 1
   ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
     "Error reading Ring table, PGE aborting, exit code = 1", "read_tables", 1 )
      call exit(1)


end subroutine read_tables

end module m_read_tables
