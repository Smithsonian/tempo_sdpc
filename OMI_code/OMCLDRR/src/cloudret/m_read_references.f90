module m_read_references

contains
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  read_references
!                
!
! !INTERFACE:
!
subroutine read_references (rc)

! !USES:
!  read_references reads precomputed references spectra needed for cloud
!               pressures
!
use m_vars, ONLY: nwav, nscanpos, reference_spec, reference_wave, &
                  input_data_path, refl_Ref, psurf_ref, sza_ref, satz_ref, &
                  az_Ref, lat_ref, lon_ref, thresh_file, filename, ex, iprt, &
                  reference_ring, reference_rad, ntheta, theta, nscan, scan, &
                  nphi, phi, npres, pres
use m_LUN_set
use m_interp_ring_rad
use m_interp_pres
      Implicit NONE

! !INPUT PARAMETERS:
!

! !OUTPUT PARAMETERS:
!

      integer, intent(out)         :: rc        ! Error return code:
                                                !  0   all is well
                                                !  1   files not found

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  25aug04   Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------
integer :: i, ii
character(len=100) :: text

integer, parameter :: lun=11
integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
integer :: status, version, ierr, OMCLDRR_F_FAILURE
include 'PGS_IO.f'
include 'PGS_IO_1.f'
include 'PGS_OMI_1900.f'
include 'PGS_SMF.f'

if (ex) then
version = 1
status = pgs_io_gen_openf ( ref_id, PGSd_IO_Gen_RSeqFrm, &
        0,lun, version)
if (iprt > 0) then
 print *,'read_references: trying to open reference file ',status, lun
endif
if(status.ne.0) then
  ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening reference table file', &
  'read_references, module m_read_references',1)
  ierr = OMI_SMF_setmsg( status, &
                            "PGE aborting, exit code = 1", "read_references", 1 ) 
 call exit(1)
endif
else
 filename=trim(input_data_path)//trim(thresh_file)
  print *,'default output file name = ',filename
 open(lun,file=filename,status='old', form="formatted", iostat=rc)
 if (rc /= 0) then
   if (iprt > 0) print *,'read_references: error reading ',filename
   stop
 endif
endif
 read(lun, *, err=100) text
 read(lun, *, err=100) nscanpos, nwav
 if (iprt >= 1) then
   print *,text
   print *,'read_references: nscanpos, nwav'
   print *,nscanpos, nwav
 endif

!allocate the reference table, but fill with default
!constant value in case file not found
!===================================================
 if (.not. allocated(reference_spec)) then
  allocate(reference_wave(nwav,nscanpos))
  allocate(reference_spec(nwav,nscanpos))
  allocate(refl_ref(nscanpos))
  allocate(psurf_ref(nscanpos))
  allocate(sza_ref(nscanpos))
  allocate(satz_ref(nscanpos))
  allocate(az_ref(nscanpos))
  allocate(lat_ref(nscanpos))
  allocate(lon_ref(nscanpos))
 endif
 do i=1, nscanpos
    read(lun, *,err=100) ii, refl_ref(i), psurf_ref(i), sza_ref(i), satz_ref(i), az_ref(i), lat_ref(i), lon_ref(i)
    if (iprt > 2) &
    print *, i, refl_ref(i), psurf_ref(i), sza_ref(i), satz_ref(i), az_ref(i), lat_ref(i), lon_ref(i)
 
    read(lun, *, err=100) reference_wave(1:nwav,i) 
    read(lun, *, err=100) reference_spec(1:nwav,i) 
    if (iprt >= 3) print *, 'read_references: ',i, reference_wave(1:nwav, i)
    if (iprt >= 3) print *, 'read_references: ',i, reference_spec(1:nwav, i)
 enddo

if (ex) then
  status = pgs_io_gen_closef (lun)
else
  close(lun)
endif
if (iprt > 2) print *,'read_references: done reading file'

return

100 status = 1
   if (iprt > 0) print *,'read_references: error reading file'
   ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
     "Error reading reference table, PGE aborting, exit code = 1", "read_references", 1 )
      call exit(1)

end subroutine read_references

end module m_read_references
