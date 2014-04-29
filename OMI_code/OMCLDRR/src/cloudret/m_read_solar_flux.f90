module m_read_solar_flux

contains

subroutine read_solar_flux()

   use m_vars, ONLY: ws, fs, wmin2, wmax2, nsolwave, meas_qual_flg, &
            irr_quality_flagL, gomi, iprt, dist_rad, dist_irrad, config_irr 
   USE hdfeos4_parameters
   USE L1B_Reader_class
   USE m_LUN_set
   USE m_lambda_qual
   USE m_earth_sun_dist
   USE m_swathnames

   IMPLICIT NONE
 
   INCLUDE 'PGS_PC.f'
   INCLUDE 'PGS_PC_9.f'
   INCLUDE 'PGS_SMF.f'
   INCLUDE 'PGS_OMI_1900.f'
   INCLUDE 'PGS_OMCLDRR_52251.f'

   INTEGER (KIND = 4) :: version, status, pgs_pc_getreference, ierr, &
                         nTimes, nXtrack, nWavel, nWavelCoef, &
                         nwl        !,iLine
   CHARACTER (LEN = 200) :: filename, swathname
   TYPE (L1B_block_type) :: blk
   INTEGER (KIND = 4), PARAMETER :: zero = 0
   INTEGER (KIND = 2) :: mflg
 
   REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: irradianceL, &
                                                   wavelengthL
   REAL (KIND = 4) :: Wl_vis_beg, Wl_vis_end

 Wl_vis_beg = wmin2 
 Wl_vis_end = wmax2
  
! obtain name of IRR1B data file
   version = 1
   status = pgs_pc_getreference( IRR1B_FILE, version, filename )
   IF( status .NE. PGS_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, & 
      "get IRR1B_FILE name failed, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF
 
! open data block structure with default size of 1 lines
   if (vis .or. visz) then
     swathname = sunvisswath
!  if (vis) then
!  elseif (visz) then
!    swathname = sunvisswathz
!  elseif (uvsz) then
!    swathname = sunuv2swathz
   elseif (gomi) then
     swathname = "Sun Volume UV-2 Swath"  
   else
     swathname = sunuv2swath 
   endif
   if (iprt >= 2) then
    print *,'opening ',filename,' ',swathname
   endif
print *,'solar: ',filename,swathname
   status = L1Br_open( blk, filename, swathname )
   IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, & 
      "L1Br_open failed, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF
 
! obtain sizes of dimensions defined in swath
   status = L1Br_getSWdims( blk, NumTimes_k=nTimes,  nXtrack_k=nXtrack, &
     nWavel_k=nWavel, nWavelCoef_k=nWavelCoef )
   IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
      "IRR1Br_getSWdims failed, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   else
     if (iprt .ge. 2) print *,'read_solar_flux: nwavel, nwavelcoef, nXtrack ', &
        nWavel, nWavelCoef, nXtrack, nTimes
   END IF
 
! allocate memory for arrays

   ALLOCATE( irradianceL(nWavel,nXtrack), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, & 
      "irradianceL allocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )      
      call exit(1)
    END IF

   IF (ALLOCATED (irr_quality_flagL)) then
     DEALLOCATE( irr_quality_flagL, STAT=ierr )
     IF( ierr .NE. zero ) THEN
       ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
       "irr_quality_flagL deallocation failure, PGE aborting, exit code = 1", &
       "read_solar_flux", 1 )
       call exit(1)
     END IF
   END IF
   ALLOCATE( irr_quality_flagL(nWavel,nXtrack), STAT=ierr )
   irr_quality_flagL=0
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "irr_quality_flagL allocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
       call exit(1)
    END IF

   ALLOCATE( wavelengthL(nWavel,nXtrack), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "wavelengthL allocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF
 
!  Initialize all local data arrays

   irradianceL(1:nWavel,1:nXtrack) = -1.0
   irr_quality_flagL(1:nWavel,1:nXtrack) = -1.0
   wavelengthL(1:nWavel,1:nXtrack) = -1.0
   
   status = L1Br_getDATA( blk, 0, &
        MeasurementQualityFlags_k=mflg,InstrumentConfigurationId_k=config_irr )!
   IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( status, &
      "L1Brd_getDATA failed, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF

! testing MeasurementQualityFlags
if(btest(mflg,0) .or. btest(mflg,1) .or. btest(mflg,3) .or. btest(mflg,12) &
  .or. btest(mflg,11)) then
      ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, & 
      "MeasurementQualityFlags: error, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
endif

!JJ to Sasha check it, checking bit 10 twice, should check 11?
if(btest(mflg,2) .or. btest(mflg,4) .or. btest(mflg,5) .or. btest(mflg,6) &
.or. btest(mflg,7) .or. btest(mflg,8) .or. btest(mflg,9) .or. btest(mflg,10) &
) meas_qual_flg(:)=IBSET(meas_qual_flg(:),5)

   status = L1Br_getSIGline( blk, 0, Wlmin_k=Wl_vis_beg, Wlmax_k=Wl_vis_end, &
                             Signal_k=irradianceL, &
                             PixelQualityFlags_k=irr_quality_flagL, &
                             Wavelength_k=wavelengthL, Nwl_k=nwl )
   IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, &
      "L1Brd_getSIGline failed, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF

nsolwave=nwl
if (iprt .ge. 2) print *,'nsolwave ',nsolwave

   ALLOCATE( ws(0:nsolwave-1,0:nXtrack-1), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "solar wave allocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF

   ALLOCATE( fs(0:nsolwave-1,0:nXtrack-1), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "solar flux allocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF

ws(0:nsolwave-1,:)=wavelengthL(1:nsolwave,:)
fs(0:nsolwave-1,:)=irradianceL(1:nsolwave,:)
 
!set processing qulity flags
call bad_irrad_lambda(nXtrack)

! deallocate memory
 
   DEALLOCATE( irradianceL, STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "irradianceL deallocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF
 
 
   DEALLOCATE( wavelengthL, STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "wavelengthL deallocation failure, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF

! close data block structure
   status = L1Br_close( blk )
   IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( status, &
      "L1Br_close failed, PGE aborting, exit code = 1", &
      "read_solar_flux", 1 )
      call exit(1)
   END IF

  call EarthSunDist(filename,dist_irrad, dist_rad)
  fs(0:nsolwave-1,:)=fs(0:nsolwave-1,:)*(dist_irrad/dist_rad)**2

!      ierr = OMI_SMF_setmsg( PGS_S_SUCCESS, &
      ierr = OMI_SMF_setmsg( PGS_S_SUCCESS, &
      "Solar Irradiance", "read_solar_flux", 1 )

end subroutine read_solar_flux

end module m_read_solar_flux
