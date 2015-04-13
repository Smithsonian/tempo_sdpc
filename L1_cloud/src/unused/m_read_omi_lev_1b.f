       module m_read_omi_lev_1b

       public read_omi_lev_1b
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Program read_omi_lev_1b.f:
c   This program reads in the OMI UV-2 Swath data fields from a Level 1B
c   HDF-EOS file.  Such a file might be called "omi_lev_1b.hdf", and
c   might have been produced by a program called "write_omi_lev_1b.f".
c
c Geolocation Data Arrays:
c   Time(nTimes)                - time of exposure, Toolkit time (r*8)
c   SecInDay(nTimes)            - time of exposure, seconds in day (r*4)
c   SpaCraLat(nTimes)           - spacecraft latitude (r*4)
c   SpaCraLon(nTimes)           - spacecraft longitude (r*4)
c   SpaCraAlt(nTimes)           - spacecraft altitude (r*4)
c   Lat(nXtrack, nTimes)        - ground pixel latitude (r*4)
c   Lon(nXtrack, nTimes)        - ground pixel longitude (r*4)
c   SolZenAng(nXtrack, nTimes)  - ground pixel solar zenith angle (r*4)
c   SolAziAng(nXtrack, nTimes)  - ground pixel solar azimuth angle (r*4)
c   ViewZenAng(nXtrack, nTimes) - ground pixel viewing zenith angle (r*4)
c   ViewAziAng(nXtrack, nTimes) - ground pixel viewing azimuth angle (r*4)
c   TerHei(nXtrack, nTimes)     - ground pixel terrain height (i*2)
c   GeoFlags(nXtrack, nTimes)   - ground pixel geolocation flags (i*2)
c
c Radiance Data Arrays:
c   RadMant(nWavel, nXtrack, nTimes)           - radiance mantissa (i*2)
c   RadPrecMant(nWavel, nXtrack, nTimes)       - radiance precision mantissa (i*2)
c   RadExpo(nWavel, nXtrack, nTimes)           - radiance exponent (i*1)
c   WavelCoef(nWavelCoef, nXtrack, nTimes)     - wavelength coefficient (r*4)
c   WavelCoefPrec(nWavelCoef, nXtrack, nTimes) - wavelength coefficient precision (r*4)
c   UVsmallPixRad(nTimesSmallPixel, nXtrack)   - UV small pixel radiance (r*4)
c   UVsmallPixWavel(nTimesSmallPixel, nXtrack) - UV small pixel wavelength (r*4)
c   AutoFlags(nWavel, nXtrack, nTimes)         - automated flags (i*2)
c   ReadOutDisc(nTimes)                        - read out discipline (i*1)
c
c Author:
c   Peter J.T. Leonard, SSAI
c   Joanna Joiner, GSFC, Modified to be more generic
c
c Date:
c   June  16, 2001
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
c Declare the geolocation data arrays.
       !real (kind=8), dimension(:), pointer :: Time
       real (kind=4), dimension(:), pointer :: SecInDay
       real (kind=4), dimension(:), pointer :: SpaCraLat
       real (kind=4), dimension(:), pointer :: SpaCraLon
       real (kind=4), dimension(:), pointer :: SpaCraAlt
!      real (kind=4), dimension(:,:), pointer :: Latitude
!      real (kind=4), dimension(:,:), pointer :: Longitude 
!      real (kind=4), dimension(:,:), pointer :: SolZenAng
!      real (kind=4), dimension(:,:), pointer :: SolAziAng
!      real (kind=4), dimension(:,:), pointer :: ViewZenAng
!      real (kind=4), dimension(:,:), pointer :: ViewAziAng
       !integer*2 TerHei(nXtrack, nTimes)
       !integer*2 GeoFlags(nXtrack, nTimes)
 
c Declare the radiance data arrays.
       integer (kind=2), dimension(:,:,:), pointer :: RadMant
       integer (kind=2), dimension(:,:,:), pointer :: RadPrecMant
       integer (kind=1), dimension(:,:,:), pointer :: RadExpo
       !real (kind=4), dimension(:,:,:), pointer :: wave
       real (kind=4), dimension(:,:,:), pointer :: WavelCoef
       real (kind=4), dimension(:,:,:), pointer :: WavelCoefPrec
       !real*4 UVsmallPixRad(nTimesSmallPixel, nXtrack)
       !real*4 UVsmallPixWavel(nTimesSmallPixel, nXtrack)
       !integer*2 AutoFlags(nWavel, nXtrack, nTimes)
       !integer*1 ReadOutDisc(nTimes)
 
c Declare the HDF-EOS file and swath identification numbers, and
c the status of the HDF-EOS functions calls.
       integer (kind=4) swfid, swid
 
c Declare the HDF-EOS functions.
!      include "HdfEosDef.h" 
 
c Declare and set some parameters to be used by the HDF-EOS routines.
       integer (kind=4) DFACC_READ
       parameter (DFACC_READ = 1)
c DFACC_READ = 1 means that the named HDF-EOS file is read only.

       contains

       subroutine L1Br_getGEOsw (infile, swathname, ierr, 
!    .  latitude, longitude, SolZenAng, SolAziAng, ViewZenAng,
!    .  ViewAziAng, 
     .  iprt_in)!, solar, old)

       use m_vars, ONLY: lat, lon, SZA, SAzimuth, sat_Zen, 
     .  VAzimuth, geoflg, mflg, time, terr_height, config_rad, anomflg
       use m_read_swath_field
       implicit none
!      real (kind=4), dimension(:,:), intent(out) :: Latitude
!      real (kind=4), dimension(:,:), intent(out) :: Longitude 
!      real (kind=4), dimension(:,:), intent(out) :: SolZenAng
!      real (kind=4), dimension(:,:), intent(out) :: SolAziAng
!      real (kind=4), dimension(:,:), intent(out) :: ViewZenAng
!      real (kind=4), dimension(:,:), intent(out) :: ViewAziAng
       integer, optional :: iprt_in

       integer (kind=4) swopen, swattach, swdetach, swclose
       integer (kind=4) status
       integer :: iprt
!inputs
       character(len=*), intent(in) :: infile, swathname
!      logical, optional :: solar, old
!ouputs
       integer, intent(out) :: ierr
 
c Declare the counters for the do loops.
       integer (kind=4) :: i, k

!      logical :: sol ! flag for reading solar file
!      logical :: old_format ! flag for old Veefkind format 
!      character(len=3) :: pre
!      character(len=32) :: fieldname

c Assign values to the start, stride and edge arrays that correspond
c to the five types of data arrays.
!      sol=.false.
!      if (present(solar)) then
!        if (solar) sol=.true.
!      endif
       iprt=0
       if (present(iprt_in)) then
         iprt=iprt_in
       endif
 
c Open the OMI Level 1B HDF-EOS input file.
       if (iprt > 0) write(6,*) 'opening ',infile
       swfid = swopen (infile, DFACC_READ)
       if (iprt > 0) write (6, *) 'swfid ', swfid
 
c Attach to the UV-2 swath.
       swid = swattach (swfid, swathname)
       if (iprt > 0) write (6, *) 'swid ', swid, swathname

c Initialize error code to 0
           ierr=0
 
c Read in the geolocation data fields.
 
c Read in the time of exposure (seconds in day) data field.
!      status = get_data (swid, "SecondsInDay",
!    *     SecInDay)
 
c Read in the spacecraft latitude data field.
!      status = get_data (swid, "SpacecraftLatitude",
!    *     SpaCraLat)
 
c Read in the spacecraft longitude data field.
!      status = get_data(swid, "SpacecraftLongitude",
!    *     SpaCraLon)
 
c Read in the spacecraft altitude data field.
!      status = get_data(swid, "SpacecraftAltitude",
!    *     SpaCraAlt)
 
c Read in the ground pixel latitude data field.
       status = get_data(swid,"Latitude",Lat)
 
c Read in the ground pixel longitude data field.
       status = get_data(swid,"Longitude",Lon)
 
c Read in the ground pixel solar zenith angle data field.
       status = get_data(swid,"SolarZenithAngle",SZA)
 
c Read in the ground pixel solar azimuth angle data field.
       status = get_data(swid,"SolarAzimuthAngle",SAzimuth)
 
c Read in the ground pixel viewing zenith angle data field.
       status = get_data(swid,"ViewingZenithAngle",sat_Zen)
 
c Read in the ground pixel viewing azimuth angle data field.
       status = get_data(swid,"ViewingAzimuthAngle",VAzimuth)

c Read in the ground pixel quality flag.
       status = get_data(swid,"GroundPixelQualityFlags",geoflg)

c Read in the ground pixel quality flag.
       status = get_data(swid,"XTrackQualityFlags",anomflg)

c Read in the measurement quality flag.
       status = get_data(swid,"MeasurementQualityFlags",mflg)

c Read in the time field.
       status = get_data(swid,"Time",time)

c Read in the ground pixel terrain height data field.
      status = get_data (swid, "TerrainHeight", terr_height)

c Read in the InstrumentConfigurationId field.
       status = get_data(swid,"InstrumentConfigurationId",config_rad)

c Detach from the swath interface.
       status = swdetach (swid)
       if (iprt > 0) write (6, *) 'L1Br_getGEOsw: detaching swath ', status
 
c Close the OMI Level 1B HDF-EOS input file.
       status = swclose (swfid)
       if (iprt > 0) write (6, *) 'L1Br_getGEOsw: closing file ',status
       end subroutine L1Br_getGEOsw 
 
c Read in the ground pixel geolocation flags data field.
       !status = get_data (swid, "GeoFlags",
!      status = get_data (swid, "QFlags", GeoFlags)
!     endif ! not sol
       subroutine L1Br_getRADsw (infile, swathname, ierr, 
     .  wmin, wmax, nwl, nWavel, iprt_in)!
 
       use m_vars, ONLY: f1, w1
       use m_read_swath_field
       use m_find
       implicit none

!inputs
       character(len=*), intent(in) :: infile, swathname
       integer, intent(out) :: nwl
       integer, intent(in) :: nWavel
       real (kind=4), intent(in) :: wmin, wmax
!      real (kind=4), dimension(:,:,:), intent(out) :: rad, wave
       integer, intent(out) :: ierr
       integer, intent(in), optional :: iprt_in

       integer (kind=4) swopen, swattach, swdetach, swclose
       integer (kind=4) :: status
       integer :: ngoodwave
       integer, parameter :: maxWavel = 800
       integer, dimension(maxWavel) :: ind 
       real, dimension(:), allocatable :: wave1
       integer :: i, j, k, n
       integer :: iprt

       iprt=0
       if (present(iprt_in)) then
         iprt=iprt_in
       endif
 
c Open the OMI Level 1B HDF-EOS input file.
       if (iprt > 0) write(6,*) 'opening ',infile
       swfid = swopen (infile, DFACC_READ)
       if (iprt > 0) write (6, *) 'swfid ', swfid
 
c Attach to the UV-2 swath.
       swid = swattach (swfid, swathname)
       if (iprt > 0) write (6, *) 'swid ', swid, swathname

c Initialize error code to 0
           ierr=0
 
c Read in the radiance mantissa data field.
       call pzeitbeg('rdrads')
       status = get_data_3Di2(swid,"RadianceMantissa")
       if (iprt > 0) print *,'read RadMant ',status
       RadMant => dummy_3Di2
 
c Read in the radiance precision mantissa data field.
!      status = get_data (swid, "RadiancePrecisionMantissa",
!    *     RadPrecMant)
 
c Read in the radiance exponent data field.
       status = get_data_3Di1 (swid, "RadianceExponent")
!    *     ,RadExpo)
       RadExpo => dummy_3Di1
       if (iprt > 0) print *,'read RadExpo ',status

       !rad=RadMant * 10.0**RadExpo

c Read in the wavelength coefficient data field.
         status = get_data_3Dr4 (swid, "WavelengthCoefficient")!,
!    *     WavelCoef)
       WavelCoef => dummy_3Dr4
       if (iprt > 0) print *,'read WavelCoef ',status
       call pzeitend
       !print *, size(WavelCoef,dim=1), size(WavelCoef,dim=2), size(WavelCoef,dim=3)

       call pzeitbeg('initrd')
       f1=0.
       w1=0.
       call pzeitend
       call pzeitbeg('comprd')
       allocate(wave1(nWavel))
       DO n = 1, size(WavelCoef,dim=3)
        DO i = 1, size(WavelCoef,dim=2)
         DO k = 1, nWavel
          wave1(k) = WavelCoef(1,i,n)
          DO j = 2, size(WavelCoef,dim=1) !nWavelCoef
            !ave(k,i) = WavelCoef(1,j) + WavelCoef(2,j)*k &
            !  + WavelCoef(3,j) + (WavelCoef(4,j) + WavelCoef(5,j)*i)*i
            wave1(k) = wave1(k) + (k-1)**(j-1)*WavelCoef(j,i,n)
            !wave(k,i) = wave(1,j) + (k-1)**(q-1)*WavelCoef(q,i,j)
          ENDDO ! nWavelCoef
         ENDDO ! nWavel
         !print *,'wave1 ',wave1(1:nWavel)
         ngoodwave = count(wave1(:) >= wmin .and. wave1(:) <= wmax)
         !print *,'ngoodwave ',ngoodwave
         if (ngoodwave > 0) then
           ind(1:ngoodwave)=find2(wave1(:) >= wmin .and. wave1(:) <= wmax,ngoodwave)
           w1(0:ngoodwave-1,i-1,n)=wave1(ind(1:ngoodwave))
           f1(0:ngoodwave-1,i-1,n)=RadMant(ind(1:ngoodwave),i,n) * 
     .             10.0**RadExpo(ind(1:ngoodwave),i,n)
         else
           ierr=1 ! no good wavelengths
           !if (iprt > 1) print *,'no good wavelengths '
         endif
        ENDDO ! nXtrack
       ENDDO ! nLines
       call pzeitend

       ! just check the first pixel (close enough for now)
       deallocate(wave1)
       deallocate(WavelCoef)
       deallocate(RadMant)
       deallocate(RadExpo)

c Read in the wavelength coefficient precision data field.
       !status = get_data (swid, "WavelengthCoefficientPrecision",
!      status = get_data (swid, "WavelengthCoeficientsPrecision",
!    *     WavelCoefPrec)
 
c Read in the UV small pixel radiance data field.
!      status = get_data (swid, "UVsmallPixelRadiance",
!    *     UVsmallPixRad)
 
c Read in the UV small pixel wavelength data field.
!      status = get_data (swid, "UVsmallPixelWavelength",
!    *     UVsmallPixWavel)
 
c Read in the automated flags data field.
!      status = get_data (swid, "AutomatedFlags",
!    *     AutoFlags)
 
c Read in the read out discipline data field.
!      status = get_data (swid, "ReadOutDiscipline",
!    *     ReadOutDisc)
c Detach from the swath interface.
       status = swdetach (swid)
       if (iprt > 0) write (6, *) 'L1Br_getGEOsw: detaching swath ', status
 
c Close the OMI Level 1B HDF-EOS input file.
       status = swclose (swfid)
       if (iprt > 0) write (6, *) 'L1Br_getGEOsw: closing file ',status
 
       end subroutine L1Br_getRADsw

       end module m_read_omi_lev_1b
 
