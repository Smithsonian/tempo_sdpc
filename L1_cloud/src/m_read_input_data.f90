module m_read_input_data

private
public read_input_data

contains

subroutine read_input_data(blk, rc)

   use m_vars
   use L1B_Reader_class
   use m_read_solar_flux
   use m_swathnames
   use m_strpos
   use m_instr_config

   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  read_input_data
! 
! !DESCRIPTION: read_input_data reads level 1b data including
!		solar irradiance and observed radiance
!
! !CALLING SEQUENCE: 
!
!        call read_input_data
!     
! !INPUT PARAMETERS:   
!
! !OUTPUT PARAMETERS:  
      integer, intent(out)         :: rc        ! Error return code:
                                                !  0   all is well
                                                !  1   files not found
!
! !SEE ALSO:  
!
! !REVISION HISTORY: 
!
!  05Jan01   Joiner     original fortran 90
!  18Mar02   Vasilkov   changes to read simulated OMI Level 1B radiance
!			data and ASCII solar flux 
!  03Jan03   Vasilkov   updates for V1.2 L1B Reader
!  21Jan03   Vasilkov   updates for V1.3 L1B Reader
!
!EOP
!-------------------------------------------------------------------------
!

INTEGER :: i, izoom
!INTEGER :: il, lun=7
!character(len=1) :: buff
!character(len=255) :: buff1
!real (KIND=8) :: bufn,bufn1,bufn2,bufn3
!logical :: old
INTEGER (KIND = 4) :: PGS_TD_TAItoUTC !pgs_pc_getreference
INTEGER (KIND = 1) :: imbin
!integer, parameter :: sz=16
!integer, parameter :: ez=45
!integer, parameter :: nz=ez-sz+1

! declaration of variables used in both examples
   CHARACTER (LEN = 200) :: swathname,filenamen
   CHARACTER (LEN = 30) :: DateTime
   CHARACTER (LEN = 70) :: msg

   TYPE (L1B_block_type), intent(inout) :: blk


!***********************************************************************
include 'PGS_IO.f'
include 'PGS_IO_1.f'
INCLUDE 'PGS_SMF.f'
include 'PGS_TD_3.f'
include 'PGS_OMI_1900.f'
include 'PGS_OMCLDRR_52251.f'
!***********************************************************************

rc=0
if (form == 2) then
  filename_gome=trim(input_data_path)//filename
print *, trim(input_data_path)
print *, 'filename ',filename
call rdgome( filename_gome )

nXtrack=size(info,dim=1)
n_sol_spec=1
nwave=size(rad2,dim=2)
nsolwave=size(sol2,dim=1)
nWaveL=nwave
nTimes=1
nLines=1
call alloc_scan()

!info = (/(/pnum/),(/date/), (/utc/),(/glat/),(/glon/),(/sza/), &
!        (/scan/),(/glint/),(/los/),(/sat_ht/),(/re/)/)

date=info(1,2)
sza(:,1)=info(:,6)
lat(:,1)=info(:,4)
lon(:,1)=info(:,5)
where (lon > 180) lon = lon - 360.
sat_zen(:,1)=info(:,7)
ps=1.
azimuth=0.
iLine=1
if (allocated(ws)) deallocate (ws)
if (allocated(fs)) deallocate (fs)
ALLOCATE( ws(0:nsolwave-1,0:nXtrack-1), STAT=ierr )
ALLOCATE( fs(0:nsolwave-1,0:nXtrack-1), STAT=ierr )

do i=0, nXtrack-1
  ws(:,i)=sol2(:,1)
  fs(:,i)=sol2(:,2)
enddo
w12d(:,:)=transpose(rad2(:,:,1))
f12d(:,:)=transpose(rad2(:,:,2))
print *, 'wavelength range'
write(6,102) w12d(0,0), w12d(nwave-1,0)
print *, 'wavelength range solar'
write(6,102) ws(0,0), ws(nsolwave-1,0)
if (max_lines > 0) nXtrack=max_lines

endif

if (form == 5) then ! use new l1b reader

 if (iLine == 0) then
!**********************************************************************
  filenamen=trim(input_data_path)//filename
  vis  = strpos (filename, 'BRVG') > 0
  visz = strpos (filename, 'BRVZ') > 0
  uvsz = strpos (filename, 'BRUZ') > 0
  if (visz) then
     swathname = visswathz
  else if (uvsz) then
     swathname = uv2swathz
  else if (vis) then
     swathname = visswath
  else
     swathname = uv2swath
  endif
  gomi = index(filename, 'GOMI') > 0
  if (wrt_solar) then
    wmin=355
    wmax=500
  endif
  if (cloud_clear) then
    wmin=383!360
    wmax=385!405
  endif
  if (.not. vis .and. .not. visz) then
   if (.not. gomi) then
    wmin2 = 330.! 356.
    wmax2 = 367.
    if (.not. set_wmin) & 
      !wmin = 345.5 
      wmin = 358. !355.5 ! 356.
    if (.not. set_wmax) & 
      !wmax = 354.5 
      wmax = 363. !365. !362.5
    if (wrt_solar) then
      wmin2=310
      wmax2=375
    endif ! write_solar
    wave_long=362.5 !353.4
    wave_short=345.4
   else
    wmin = 355.! 356.
    wmax = 365.
    wave_long=364.0
    wave_short=340.4
    shift=.true.
    squeeze=.true. 
   endif ! gomi
  endif ! not vis
 if (iprt > 0) print *,'read_input_data: filename ',filenamen, swathname
!**********************************************************************
   status = L1Br_open( blk, filenamen, swathname )
 if (iprt > 0)   print *,'read_input_data: opening l1b status ',status
   IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( status, & 
      "PGE aborting, exit code = 1", "read_input_data", 1 )
      call exit(1)
   ELSE
!     ierr = OMI_SMF_setmsg( OMCLDRR_S_SUCCESS, &
     ierr = OMI_SMF_setmsg( PGS_S_SUCCESS, &
      "Opened Earth Radiance", "read_input_data", 1 )
  END IF

! obtain sizes of dimensions defined in swath
   status = L1Br_getSWdims( blk, NumTimes_k=nTimes, nXtrack_k=nXtrack, &
      nWavel_k=nWavel, nWavelCoef_k=nWavelCoef )
 if(iprt > 0) print *,'read_input_data: nTimes, nXtrack, nWavel, nWavelCoef '
 if(iprt > 0) print *, nTimes,nXtrack,nWavel,nWavelCoef
 iLine=start_line
 if (max_lines > 0 .and. iprt > 0) then
   print *,'read_input_data: changing nTimes to ',max_lines
   nTimes=max_lines+start_line
 endif
 IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( status, & 
      "PGE aborting, exit code = 1", "read_input_data", 1 )
      call exit(1)
 ELSE
!      ierr = OMI_SMF_setmsg( OMCLDRR_S_SUCCESS, &
      ierr = OMI_SMF_setmsg( PGS_S_SUCCESS, &
      "Read Earth Radiance Dims", "read_input_data", 1 )
 END IF

!*********************************************************************
 nLines=nTimes
!*********************************************************************

 if (do_alloc2) then
    call alloc_scan()
 endif ! alloc

endif ! if iLine == 0

!call pzeitbeg('rd_geo')

status = L1Br_getDATA ( blk, iLine-1, MeasurementQualityFlags_k=mflg(iLine), &
         MeasurementClass_k=meas_class(iLine), &
         InstrumentConfigurationID_k=config_rad(iLine), &
         ImageBinningFactor_k=imbin )
IF( status .NE. OMI_S_SUCCESS ) THEN
         ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
           "L1Brd_getDATA failed", "read_input_data", 1 )
         ierr = OMI_SMF_setmsg( status, & 
           "PGE aborting, exit code = 1", "read_input_data", 1 )
         call exit(1)
END IF

status = L1Br_getGEOline( blk, iLine-1, Time_k=time(iLine), & 
         Latitude_k=lat(:,iLine), Longitude_k=lon(:,iLine), &
         SolarZenithAngle_k=sza(:,iLine), SolarAzimuthAngle_k=sazimuth(:,iLine), &
         ViewingZenithAngle_k=sat_zen(:,iLine), ViewingAzimuthAngle_k=vazimuth(:,iLine), &
         TerrainHeight_k=terr_height(:,iLine), GroundPixelQualityFlags_k=geoflg(:,iLine), &
         XTrackQualityFlags_k=anomflg(:,iLine)) 
          !Geoflag_k=geoflg(:,iLine), MeasFlag_k=mflg(iLine), &
         !MeasClass_k=meas_class(iLine), Config_k=config_rad(iLine) , ImgBinFact_k=imbin )
IF( status .NE. OMI_S_SUCCESS ) THEN
         ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
           "L1Brd_getGEOline failed", "read_input_data", 1 )
         ierr = OMI_SMF_setmsg( status, & 
           "PGE aborting, exit code = 1", "read_input_data", 1 )
         call exit(1)
END IF


!call pzeitend

call instr_config(ierr,izoom)
if(iprt>1) print *,'instrum config compatability code',ierr,'izoom',izoom

!check for skipping the zoom mode
if(izoom==1 .and. .not. do_zoom) then
  ierr = OMI_SMF_setmsg( OMCLDRR_W_ZOOM, &
         "zoom mode measurements skipped ", "read_input_data", 1 )
  call exit(0)  
endif

!JJ - is the 180 correct?
azimuth(:,iLine)=sazimuth(:,iLine)+180.0-vazimuth(:,iLine)
where(azimuth(:,iLine) < -180.) azimuth(:,iLine)=azimuth(:,iLine)+360.
where(azimuth(:,iLine) > 180.) azimuth(:,iLine)=azimuth(:,iLine)-360.
azimuth(:,iLine)=abs(azimuth(:,iLine))
where(azimuth(:,iLine) > 360.0) azimuth(:,iLine)=fill_value

!==============================================================================
! Find Day, Month and Year that the Input data was collected
!==============================================================================
if (iLine == start_line) then
   status = PGS_TD_TAItoUTC(time(1),DateTime)

   IF(status .NE. PGS_S_SUCCESS) THEN
     ierr = OMI_SMF_setmsg(OMI_E_FAILURE, "TAI time conversion failed", &
        "read_input_data", 0)
     month=1
   ELSE
  10 FORMAT (I4,1X,I2,1x,I2,17X)
     READ  (DateTime,10) Year, Month, Day
     WRITE( msg,* ) "Date is: ", Year, Month, Day
     status = OMI_SMF_setmsg(PGS_S_SUCCESS, msg, "read_input_data", 1)
      if (iprt >= 1) print *, msg
   ENDIF
   if (iprt >= 3) print *,'wmin2 wmax2 ',wmin2,wmax2
endif ! get month

n_input = n_input + nXtrack
! check MeasurementQualityFlags for missing data
!=====================================================
if(btest(mflg(iLine),0) .or. btest(mflg(iLine),1) &
  .or. btest(mflg(iLine),3) .or. btest(mflg(iLine),12)) then

  ! if missing data, skip the line, set radiance measurement error flag
  ! ===================================================================  
  meas_qual_flg(iLine)=IBSET(meas_qual_flg(iLine),1)
  n_missing = n_missing + nXtrack 
  rc=1
  if (iprt >= 1) print *,'missing line ',iLine, btest(mflg(iLine),0), btest(mflg(iLine),1), &
                btest(mflg(iLine),3), btest(mflg(iLine),12)
!  goto 999
endif

! check other MeasurementQualityFlags and set our measurement quality flag
!=========================================================================
if(btest(mflg(iLine),2) .or. btest(mflg(iLine),4) .or. btest(mflg(iLine),5) .or. btest(mflg(iLine),6) &
.or. btest(mflg(iLine),8) .or. btest(mflg(iLine),9) .or. btest(mflg(iLine),11)) &
meas_qual_flg(iLine)=IBSET(meas_qual_flg(iLine),2)

if(btest(mflg(iLine),7)) meas_qual_flg(iLine)=IBSET(meas_qual_flg(iLine),3)
if(btest(mflg(iLine),10)) meas_qual_flg(iLine)=IBSET(meas_qual_flg(iLine),4)

!read a line of data
!====================
!call pzeitbeg('rd_rad')
status = L1Br_getSIGline( blk, iLine-1, Wlmin_k=wmin2, Wlmax_k=wmax2, &
                          Signal_k=f12d, &! rad_precisionL, &
                          PixelQualityFlags_k=quality_flagL, &
                          Wavelength_k=w12d, Nwl_k=nwl )

!check if it was good and abort if bad
!=====================================
IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
         "L1Br_getSIGline failed", "read_input_data", 1 )
      ierr = OMI_SMF_setmsg( status, & 
         "PGE aborting, exit code = 1", "read_input_data", 1 )
      call exit(1)
END IF

!call pzeitend

if(iLine==start_line)  ll=1
!Next two lines fail under gfortran
!presuming goal is to copy over part of w12d and f12d, restrict input ranges
!w1(0:nwl-1,:,ll)=w12d
!f1(0:nwl-1,:,ll)=f12d
w1(0:nwl-1,:,ll)=w12d(0:nwl-1,:)
f1(0:nwl-1,:,ll)=f12d(0:nwl-1,:)
w12d(nwl:nWavel-1,:)=0.
f12d(nwl:nWavel-1,:)=0.

!print check of wavelengths
!==================================
if (iprt >= 1 .and. iLine == start_line) then
  print *, 'nwl, iLine, Wmin, Wmax'
  print *, nwl, iLine, Wmin, Wmax
endif

!if (iprt >= 5) then
! print *,'wavelength ',w1(:,:,iLine)
! print *,'radiance ', f1(:,:,iLine)
!endif

nwave=nwl
if (iLine == start_line) then
   if (iprt >= 3) then
     write(6,102) w12d(0:nwl-1,0)
   endif
else 
 !check for missing data
 !=======================
 if (nwl > nWavel .or. nwl < min_wl) then
   qc(:,iLine) = IBSET(qc(:,iLine),14)
   n_missing = n_missing + nXtrack 
   rc=2
   if (iprt >= 1) print *,'missing line ',iLine, nwl, nWavel, min_wl
   if (iprt >= 3) then
     write(6,102) w12d(0:nwl-1,0)
   endif
 endif ! missing wavelength data
endif ! start_line

!999 continue

! read solar flux
!===================
if (iLine == start_line) then
!  call pzeitbeg('rd_sol')
  call read_solar_flux()
!  call pzeitend
  if (iprt > 1) then
    print *,'irradiance'
    do i=0,nsolwave-1
     write(*,'(i4,2e12.4)') i,ws(i,0),fs(i,0)
    enddo
  endif ! iprt > 1
endif ! iLine==start_line

endif ! different formats
!************************************************************************

if (iLine == start_line) then
  if (wrt_solar) call write_solar()
endif

!100 format(i6,6f10.1)
!101 format(i6,f10.1,e12.3)
102 format(6f12.2) 

end subroutine read_input_data

subroutine alloc_scan()

use m_vars
use L1B_Reader_class

implicit none

include 'PGS_OMI_1900.f'
include 'PGS_OMCLDRR_52251.f'

! allocate memory for arrays
   if (allocated(lat)) deallocate (lat)   
   ALLOCATE( lat(nXtrack,nlines), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "latitude allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   if (allocated(lon)) deallocate (lon)   
   ALLOCATE( lon(nXtrack,nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "longitude allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   if (allocated(sza)) deallocate (sza)
   ALLOCATE( sza(0:nXtrack-1,nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, & 
      "szenith allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   if (allocated(sat_zen)) deallocate (sat_zen)
   ALLOCATE( sat_zen(0:nXtrack-1,nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, & 
      "vzenith allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( sazimuth(nXtrack,nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "sazimuth allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( vazimuth(nXtrack,nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "vazimuth allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( terr_height(nXtrack,nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "terrain height allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( geoflg(nXtrack,nLines), STAT=ierr )
   geoflg=0 
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "geoflg allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( anomflg(nXtrack,nLines), STAT=ierr )
   anomflg=0
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "anomflg allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( mflg(nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN 
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "measflg allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( rad_precisionL(nWavel,nXtrack), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "rad_precisionL allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( quality_flagL(nWavel,nXtrack), STAT=ierr )
   quality_flagL=0
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "quality_flagL allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

  if (allocated(w1)) deallocate (w1)
   ALLOCATE( w1(0:nWavel-1,0:nXtrack-1,ny), STAT=ierr ) ; w1=0.
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "wavelengthL allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF
   
   if (allocated(w12d)) deallocate (w12d)
   ALLOCATE( w12d(0:nWavel-1,0:nXtrack-1), STAT=ierr )   ; w12d = 0.0
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "wavelengthL allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   if (allocated(f1)) deallocate (f1)
   ALLOCATE( f1(0:nWavel-1,0:nXtrack-1,ny), STAT=ierr ) ; f1=0.
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "radianceL allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   if (allocated(f12d)) deallocate (f12d)
   ALLOCATE( f12d(0:nWavel-1,0:nXtrack-1), STAT=ierr ) ; f12d=0.
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "radianceL allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF
   
   ALLOCATE( time(nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "Time allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

   ALLOCATE( meas_class(nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "MeasurementClass allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF
   ALLOCATE( config_rad(nLines), STAT=ierr )
   IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
      "InstrumentConfigurationID allocation failure, PGE aborting, exit code = 1", &
      "read_input_data", 1 )
      call exit(1)
   END IF

if (allocated(meas_qual_flg)) deallocate (meas_qual_flg)
allocate  (meas_qual_flg(nLines)) ; meas_qual_flg = 0
if (allocated(cloud_pres)) deallocate (cloud_pres)
allocate (cloud_pres (0:nXtrack - 1,nLines)) ; cloud_pres=fill_value
if (allocated(azimuth)) deallocate (azimuth)
allocate (azimuth (0:nXtrack-1,nLines))      ; azimuth=fill_value
if (allocated(refl)) deallocate (refl)  
allocate (refl    (0:nXtrack-1,nLines))      ; refl=fill_value
if (allocated(dIdR)) deallocate (dIdR)  
allocate (dIdR    (0:nXtrack-1,nLines))      ; dIdR=fill_value
if (allocated(ps)) deallocate (ps)     
allocate (ps      (0:nXtrack-1,nLines))      ; ps=fill_value
if (allocated(ref_clr)) deallocate (ref_clr)     
allocate (ref_clr      (0:nXtrack-1,nLines))      ; ref_clr=fill_value
if (allocated(ai)) deallocate (ai)      
allocate (ai      (0:nXtrack-1,nLines))      ; ai=fill_value
if (allocated(reflect_cld)) deallocate (reflect_cld)      
allocate (reflect_cld      (0:nXtrack-1,nLines)) ; reflect_cld=fill_value
if (allocated(rad_cld_frac)) deallocate (rad_cld_frac)
allocate (rad_cld_frac(0:nXtrack-1,nLines))  ; rad_cld_frac=fill_value
if (allocated(eff_cld_frac)) deallocate (eff_cld_frac)
allocate (eff_cld_frac(0:nXtrack-1,nLines))  ; eff_cld_frac=fill_value
if (allocated(eff_cld_frac2)) deallocate (eff_cld_frac2)
allocate (eff_cld_frac2(0:nXtrack-1,nLines))  ; eff_cld_frac2=fill_value
if (allocated(cld_pres2)) deallocate (cld_pres2)
allocate (cld_pres2(0:nXtrack-1,nLines))  ; cld_pres2=fill_value
 if (allocated(chlorophyll)) deallocate (chlorophyll)
allocate (chlorophyll(0:nXtrack - 1,nLines)) ; chlorophyll=fill_value
if (allocated(eta)) deallocate (eta)  
allocate (eta(0:nXtrack - 1,nLines))         ; eta=fill_value
if (allocated(biases)) deallocate (biases)  
allocate (biases  (0:nXtrack-1,nLines))      ; biases=fill_value
if (allocated(biases2)) deallocate (biases2)  
allocate (biases2  (0:nXtrack-1,nLines))      ; biases2=fill_value
if (allocated(stds)) deallocate (stds) 
allocate (stds    (0:nXtrack-1,nLines))      ; stds=fill_value
if (allocated(stds2)) deallocate (stds2) 
allocate (stds2    (0:nXtrack-1,nLines))      ; stds2=fill_value
 if (allocated(rms)) deallocate (rms)   
allocate (rms     (0:nXtrack-1,nLines))      ; rms=fill_value
 if (allocated(chi_sqr)) deallocate (chi_sqr)
allocate (chi_sqr (0:nXtrack-1,nLines))      ; chi_sqr=fill_value
 if (allocated(chi_sqr2)) deallocate (chi_sqr2)
allocate (chi_sqr2 (0:nXtrack-1,nLines))      ; chi_sqr2=fill_value
 if (allocated(land_flg)) deallocate (land_flg)
!allocate (land_flg(0:nXtrack-1)) ; land_flg=-1
!!!land_flg is a logical. presumably -1 indicates false...
allocate (land_flg(0:nXtrack-1)) ; land_flg=.FALSE.
if (allocated(chlcl)) deallocate (chlcl)
allocate (chlcl   (0:nXtrack-1)) ; chlcl=fill_value
if (allocated(qc)) deallocate (qc)      
allocate (qc      (0:nXtrack-1,nLines)) ; qc=0
if (allocated(qc2)) deallocate (qc2)      
allocate (qc2      (0:nXtrack-1,nLines)) ; qc2=0
if (allocated(fill)) deallocate (fill)      
allocate (fill    (0:nXtrack-1,nLines)) ; fill=fill_value
if (allocated(shifts)) deallocate (shifts)      
allocate (shifts  (0:nXtrack-1,nLines)) ; shifts=fill_value
if (allocated(shifts2)) deallocate (shifts2)      
allocate (shifts2  (0:nXtrack-1,nLines)) ; shifts2=fill_value
if (allocated(squeezes)) deallocate (squeezes)      
allocate (squeezes(0:nXtrack-1,nLines)) ; squeezes=1

end subroutine alloc_scan

!;.........................................................................
!;Reads spectral info from a GOME orbit
!;Includes correction of solar spectrum for Doppler shift
!;
!;FILENAME =scalar string for complete directory/filename.
!;Reads files of type .BIN only.    These are found in directory
!;  /misc/jfg19/gome/small_swath/
!;If files are ftp'd to local pc, you may have to un-comment
!; the SWAP_ENDIAN statements.
!;
!;Output Keywords:
!;INFO(J,K) is information about spectrum number J
!;  K=0:  PNUM =  pixel number
!;  K=1:  DATE =  date (YYYYMMDD)
!;  K=2:  UTC  =  Coordinated UT (seconds)
!;  K=3:  GLAT =  latitude of central pixel (of 5)
!;  K=4:  GLON =  longitude of central pixel (of 5)
!;  K=5:  SZA  =  solar zenith angle at center of pixel
!;  K=6:  SCAN =  pixel number within a scan (west,centr,east)
!;  K=7:  GLINT=  =1 if sun glints, =0 otherwise
!;  K=8:  LOS  =  instrument line of sight w/resp North
!;  K=9:  SAT_HT= satellite height above earth
!;  K=10: RE   =  radius of earth under satellite 
!;
!;.....................................................................
subroutine getspec( s,r,wlmin,wlmax, sol,rad )

use m_find
use mathcons
implicit NONE     

real (KIND=8), intent(in) :: wlmin, wlmax
real(kind=4), dimension(:,:,:), intent(in) :: r
real (KIND=8), dimension(:,:), intent(in) :: s
real (KIND=8), dimension(:,:,:), pointer :: rad
real (KIND=8), dimension(:,:), pointer :: sol
  
integer :: imins, imaxs, imine, imaxe, n
!real (KIND=8) :: dum
  

n=size(r,dim=1)
imins=find1(abs(s(1,:)-wlmin) == minval(abs(s(1,:)-wlmin)))
imaxs=find1(abs(s(1,:)-wlmax) == minval(abs(s(1,:)-wlmax)))
imine=find1(abs(r(n/2+1,1,:)-wlmin) == minval(abs(r(n/2+1,1,:)-wlmin)))
imaxe=imine+imaxs-imins  
allocate(rad(n,abs(imaxe-imine)+1,2))
allocate(sol(imaxs-imins+1,2))
sol(:,1) = s(1,imins:imaxs)
sol(:,2) = s(2,imins:imaxs)
rad(:,:,1) = (r(:,1,imine:imaxe))
rad(:,:,2) = (r(:,2,imine:imaxe))
end subroutine getspec

subroutine rdgome( filename )

       use m_vars, ONLY: sol2,rad2,sol3,rad3, s2, rad2b, s3, rad3b, info, &
        s2b, s3b
       use m_strpos
       implicit NONE          

!.....................................................................
! Set I/O parameters
real (KIND=8), parameter :: wlmin2=313., wlmax2=405., wlmin3=410., wlmax3=500., &
 doppler=1.000019     !doppler WL shift (v=7.5km/sec, theta=40 deg )
!.....................................................................
character(len=*), intent(inout) :: filename
real (KIND=8), dimension(:), pointer :: sza, glat, glon

!real(kind=4) :: test
!integer :: recl
integer :: ios, i
integer :: lun1=1
real(kind=4), allocatable, dimension(:) :: glint, pnum, utc, scan, sat_ht
real(kind=4), allocatable, dimension(:) :: los, re, date
!real(kind=4), allocatable, dimension(:) :: yearday
logical :: eof
real(kind=4), dimension(10) :: geoloc
real(kind=4), dimension(6) :: losn, szan, szasc, lossc
real(kind=4), dimension(3,16) :: pmd
integer :: temp, filesize, npix, pixsiz=15204
!integer :: j,k
real(kind=4) :: npixr

!check for filesize appended onto end of filename
!==================================================
  filesize=0
  temp=strpos(filename,':')
  if (temp > 0) then
    read(filename(temp+1:len(filename)),*) filesize
    filename=filename(1:temp-1)
  endif

!gmpx = {pix, geoloc:fltarr(10),   glint:0.0,           losn:fltarr(6),   &^M   
!             lossc:fltarr(6),     pixnum:0.0,          pmd:fltarr(3,16), &^M   
!             rad2b:fltarr(2,832), rad3:fltarr(2,1024), re:0.0,           &^M   
!             sat_ht:0.0,          scan:0.0,            szan:fltarr(6),   &^M   
!             szasc:fltarr(6),     yearday:0.0,         utc:0.0  }^M   
                                                                                allocate(s2b(2,841))   ;    allocate(s3b(2,1024))
allocate(s2(2,841))   ;    allocate(s3(2,1024))
   
print *,'rd_gome: opening ',filename
open (unit= lun1, file=  filename, action='READ', iostat=ios, &
  status='old',form='unformatted' )
if (ios /= 0) then
  print *,'rd_gome: error opening file '
  stop
endif

npix =(filesize -((1682+2048)*4))/pixsiz
if (filesize /= 0) print *,'filesize ',filesize, 'pixsiz ',pixsiz
print *,'leftover ',modulo((filesize -((1682+2048)*4)),pixsiz)
print *,'reading ',filename
read (lun1,err=100) npixr
npix = npixr     
print *,'npix ',npix
read (lun1,err=100) s2b      
s2=s2b
read (lun1,err=100) s3b      
s3=s3b
allocate(glint(npix))
allocate(pnum(npix))
allocate(re(npix))
allocate(sat_ht(npix))
allocate(scan(npix))
allocate(date(npix)) 
allocate(utc(npix))
allocate(glat(npix))
allocate(glon(npix))
allocate(los(npix))
allocate(sza(npix))
allocate(rad2b(npix,2,832))
allocate(rad3b(npix,2,1024))
allocate(info(npix,11))
eof = .false.
do i=1, npix
read (lun1, err=100) geoloc, glint(i), losn, lossc, pnum(i), &
  pmd, rad2b(i,:,:), rad3b(i,:,:), re(i), sat_ht(i), scan(i), &
  szan, szasc, date(i), utc(i)
los(i)=losn(3)
sza(i)=szan(3)
glon(i)=geoloc(10)
glat(i)=geoloc(9)
enddo
close (lun1)
print *,'rd_gome: closing ',filename
   
info(:,1)=pnum
info(:,2)=date
info(:,3)=utc
info(:,4)=glat
info(:,5)=glon
info(:,6)=sza
info(:,7)=scan
info(:,8)=glint
info(:,9)=los
info(:,10)=sat_ht
info(:,11)=re

call getspec (s2, rad2b, wlmin2, wlmax2, sol2, rad2)
call getspec (s3, rad3b, wlmin3, wlmax3, sol3, rad3)

sol2(:,1)=sol2(:,1)*doppler
sol3(:,1)=sol3(:,1)*doppler

deallocate(s2)
deallocate(s3)
deallocate(glint)
deallocate(pnum)
deallocate(re)
deallocate(sat_ht)
deallocate(scan)
deallocate(date)
deallocate(utc)
deallocate(glat)
deallocate(glon)
deallocate(los)
deallocate(sza)

return

100 print *,'rd_gome: error reading file'
    stop

end subroutine rdgome

subroutine write_solar

use m_vars
integer :: ios 
character(len=80) :: sfile2

sfile2=trim(solar_path)//trim(sfile)
open(1,file=sfile2,form='formatted',status='unknown',action='WRITE',iostat=ios,err=100)
print *,'opened solar file ', sfile2
write(1,*) nsolwave
do i=0,59
write(1,*) ws(0:nsolwave-1,i)
write(1,*) fs(0:nsolwave-1,i)
enddo
print *,'wrote solar file'
close (1) 
stop

return
100 print *,'error writing solar file'
    stop

end subroutine write_solar

end module m_read_input_data
