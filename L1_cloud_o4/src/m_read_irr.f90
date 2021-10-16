!*****************
module m_read_irr
  !*****************

contains

  !1111111111111111111
  subroutine read_irr
    !1111111111111111111

    use hdfeos4_parameters
    use he5_swreader
    use m_vars

    implicit none

    !---------------------------------------------------------------------72
    ! ROUTINE: m_read_irr
    ! 
    ! DESCRIPTION: This program reads the OMI SUN Volume VIS Swath data 
    !   from a Level 1B HDF-EOS file
    !
    ! Flags
    !   - PixelQualityFlag (nw,nx,nt=1)
    !   - MeasurementQualityFlag = 0 on 12/22/2004
    !   - No xTrackQualityFlags!
    !   - No GroundPixelQualityFlags!
    !   
    !   --> avoid 1 in bit0/bit1/bit2 of PixelQualityFlag
    !
    ! Outputs:
    !   - irr_EarthSunDist
    !   - irr_out_irradiance_466nm (nx)
    !   - irr_out_irradiance_477nm (nx): option
    !
    ! REVISION HISTORY: 
    !
    !  04/23/15 Yang original fortran 90
    !---------------------------------------------------------------------72

    ! Data Type                                Element Name    L1B Dataset Name
    ! ==============================           ============    ================
    ! Geolocation fields
    !  real(kind=8)                          ::Time          ! Time
    !  real(kind=4)                          ::SecInDay      ! SecondsInDay
    !  real(kind=4),   dimension(:,:),pointer::Lon           ! Longitude
    !  real(kind=4),   dimension(:,:),pointer::Lat           ! Latitude
    !  real(kind=4),   dimension(:,:),pointer::SolZenAng     ! SolarZenithAngle
    !  real(kind=4),   dimension(:,:),pointer::SolAziAng     ! SolarAzimuthAngle
    !  real(kind=4),   dimension(:,:),pointer::ViewZenAng    ! ViewingZenithAngle
    !  real(kind=4),   dimension(:,:),pointer::ViewAziAng    ! ViewingAzimuthAngle
    !  integer(kind=2),dimension(:,:),pointer::TerHgt        ! TerrainHeight
    !  integer(kind=2),dimension(:,:),pointer::GPQFlag       ! GroundPixelQualityFlags
    !  integer(kind=1),dimension(:,:),pointer::XTQFlag       ! XTrackQualityFlags
    !
    ! Data fields
    !  integer(kind=2),dimension(:,:),pointer::RadMantissa   ! RadianceMantissa or IrradianceMantissa or SignalMantissa
    !  integer(kind=2),dimension(:,:),pointer::RadPrecision  ! RadiancePrecision or IrradiancePrecision or SignalPrecision
    !  integer(kind=1),dimension(:,:),pointer::RadExponent   ! RadianceExponent or IrradianceExponent or SignalExponent
    !  integer(kind=2),dimension(:,:),pointer::PQFlag        ! PixelQualityFlags
    !  real(kind=4),   dimension(:,:),pointer::WavelenCoef   ! WavelengthCoefficient
    !  real(kind=4),   dimension(:,:),pointer::WavelenPrec   ! WavelengthCoefficientPrecision
    !  integer(kind=2)                       ::WavelenRefCol ! WavelengthReferenceColumn
    !  integer(kind=2)                       ::MQFlag        ! MeasurementQualityFlags

    character(len=255)::filename
    !integer::ndim,ntmp
    integer(kind=4)::swfid,swid,status,ierr
    integer(kind=4)::swrdfld!,swfldinfo
    integer(kind=4)::nt,nx,nw,nwc
    integer(kind=4)::ix,iw,iw1,iw2!,it,iwc
    real(kind=4)::yy1,yy2,ww1,ww2
    real(kind=4),dimension(:,:),allocatable::wtmp,itmp

    integer(kind=4)::start1,stride1,edge1
    !integer(kind=4),dimension(2)::start2,stride2,edge2
    integer(kind=4),dimension(3)::start3,stride3,edge3

    !---------------
    ! read OMLIBIRR
    !---------------
    filename=trim(name_irr_dir)//name_irr_file
    !filename=trim(name_irr_file)
    swfid=swopen(filename,DFACC_READ)   
    swid=swattach(swfid,name_irr_swath)   
if (swfid < 0 .or. swid < 0) then
  print *, "read_irr fail"
  stop 1
endif

    status=swrdattr(swid,"NumTimes",irr_NumTimes)
    status=swrdattr(swid,"EarthSunDistance",irr_EarthSunDist)

    irr_nXtrack=swdiminfo(swid,"nXtrack")
    irr_nWavel=swdiminfo(swid,"nWavel")
    irr_nWavelCoef=swdiminfo(swid,"nWavelCoef")

    nt=irr_NumTimes
    nx=irr_nXtrack
    nw=irr_nWavel
    nwc=irr_nWavelCoef

    !------------------------------------------------------------
    ! Array(nt=1): Time, SecondsInDay, WavelengthReferenceColumn
    !              MeasurementQualityFlags=0
    !------------------------------------------------------------
    start1=0
    stride1=1
    edge1=1

    status=swrdfld(swid,"Time",start1,stride1,edge1,irr_Time)
    status=swrdfld(swid,"SecondsInDay",start1,stride1,edge1,irr_SecondsInDay)
    status=swrdfld(swid,"MeasurementQualityFlags",start1,stride1,edge1,irr_MeasurementQualityFlags)
    status=swrdfld(swid,"WavelengthReferenceColumn",start1,stride1,edge1,irr_WavelengthReferenceColumn)

    !---------------------------------------------------------------------------
    ! Array(nwc,nx,nt=1): WavelengthCoefficient, WavelengthCoefficientPrecision
    !---------------------------------------------------------------------------
    !  start2(1)=0
    !  stride2(1)=1
    !  edge2(1)=nwc
    !  start2(2)=0
    !  stride2(2)=1
    !  edge2(2)=nx

    start3=0 ! all start at index 0
    stride3=1 
    edge3=(/nwc, nx, 1/)

    allocate(irr_WavelengthCoefficient(nwc,nx),stat=ierr)
    status=swrdfld(swid,"WavelengthCoefficient",start3,stride3,edge3,irr_WavelengthCoefficient)

    allocate(irr_WavelengthCoefficientPrecision(nwc,nx),stat=ierr)
    status=swrdfld(swid,"WavelengthCoefficientPrecision",start3,stride3,edge3,irr_WavelengthCoefficientPrecision)

    !-------------------------------------------------------------
    ! Array(nw,nx,nt=1): IrradianceMantissa, IrradiancePrecision, 
    !                    IrradianceExponent, PixelQualityFlags
    !-------------------------------------------------------------
    !  start2(1)=0
    !  stride2(1)=1
    !  edge2(1)=nw
    !  start2(2)=0
    !  stride2(2)=1
    !  edge2(2)=nx

    start3=0
    stride3=1
    edge3=(/nw, nx, 1/)

    allocate(irr_IrradianceMantissa(nw,nx),stat=ierr)
    status=swrdfld(swid,"IrradianceMantissa",start3,stride3,edge3,irr_IrradianceMantissa)
    !  allocate(irr_IrradiancePrecision(nw,nx),stat=ierr)
    !  status=swrdfld(swid,"IrradiancePrecision",start3,stride3,edge3,irr_IrradiancePrecision)
    !  status=swrdfld(swid,"IrradiancePrecisionMantissa",start3,stride3,edge3,irr_IrradiancePrecision)
    allocate(irr_IrradianceExponent(nw,nx),stat=ierr)
    status=swrdfld(swid,"IrradianceExponent",start3,stride3,edge3,irr_IrradianceExponent)
    allocate(irr_PixelQualityFlags(nw,nx),stat=ierr)
    status=swrdfld(swid,"PixelQualityFlags",start3,stride3,edge3,irr_PixelQualityFlags)

    !-------
    ! close
    !-------
    ierr = swdetach(swid)
    ierr = swclose(swfid)

    !---------
    ! Outputs
    !   **index of the spectral pixel in nWavel direction, starting at 0
    !---------

    allocate(wtmp(nw,nx),stat=ierr)
    allocate(itmp(nw,nx),stat=ierr)
    allocate(irr_out_irradiance_440nm(nx),stat=ierr)
    allocate(irr_out_irradiance_466nm(nx),stat=ierr)
    allocate(irr_out_irradiance_477nm(nx),stat=ierr)

    do ix=1,nx
      do iw=1,nw
        wtmp(iw,ix)=irr_WavelengthCoefficient(1,ix) &
             +(real(iw-1)-irr_WavelengthReferenceColumn)**1*irr_WavelengthCoefficient(2,ix) &
             +(real(iw-1)-irr_WavelengthReferenceColumn)**2*irr_WavelengthCoefficient(3,ix) &
             +(real(iw-1)-irr_WavelengthReferenceColumn)**3*irr_WavelengthCoefficient(4,ix) &
             +(real(iw-1)-irr_WavelengthReferenceColumn)**4*irr_WavelengthCoefficient(5,ix)

        !    if((irr_IrradianceMantissa(iw,ix) .eq. -32767) .or. (irr_IrradianceExponent(iw,ix) .eq. -127)) then
        if(btest(irr_PixelQualityFlags(iw,ix),0) .or. &
             btest(irr_PixelQualityFlags(iw,ix),1) .or. &
             btest(irr_PixelQualityFlags(iw,ix),2)) then 
          itmp(iw,ix)=-9999.
        else
          itmp(iw,ix)=irr_IrradianceMantissa(iw,ix)*10.**irr_IrradianceExponent(iw,ix)
        end if
      enddo
    enddo

    do ix=1,nx
      iw1=-9
      iw2=-9
      do iw=1,nw
        if((w466 .le. wtmp(iw,ix)) .and. (itmp(iw,ix) .gt. 0.0)) then
          iw2=iw
          exit
        end if
      enddo

      do iw=nw,1,-1
        if((w466 .ge. wtmp(iw,ix)) .and. (itmp(iw,ix) .gt. 0.0)) then
          iw1=iw
          exit
        end if
      enddo

      yy1=itmp(iw1,ix)
      yy2=itmp(iw2,ix)
      ww1=w466-wtmp(iw1,ix)
      ww2=wtmp(iw2,ix)-w466
      irr_out_irradiance_466nm(ix)=(ww1*yy2+ww2*yy1)/(ww1+ww2)
      !    write(4,115) ix,iw1,iw2,ww1,ww2,yy1,yy2,irr_out_irradiance_466nm(ix)
    enddo

    !111 format(5(1pe15.5))
    !113 format(5i5)
    !115 format(3i5,5(1pe15.5))
    !117 format(a45,3i5,2(1pe15.5))

    write(*,*)

    ! ----------
    ! add 477 nm
    ! ----------
    do ix=1,nx
      iw1=-9
      iw2=-9
      do iw=1,nw
        if((w477 .le. wtmp(iw,ix)) .and. (itmp(iw,ix) .gt. 0.0)) then
          iw2=iw
          exit
        end if
      enddo

      do iw=nw,1,-1
        if((w477 .ge. wtmp(iw,ix)) .and. (itmp(iw,ix) .gt. 0.0)) then
          iw1=iw
          exit
        end if
      enddo

      yy1=itmp(iw1,ix)
      yy2=itmp(iw2,ix)
      ww1=w477-wtmp(iw1,ix)
      ww2=wtmp(iw2,ix)-w477
      irr_out_irradiance_477nm(ix)=(ww1*yy2+ww2*yy1)/(ww1+ww2)
    enddo

    ! ----------
    ! add 440 nm
    ! ----------
    do ix=1,nx
      iw1=-9
      iw2=-9
      do iw=1,nw
        if((w440 .le. wtmp(iw,ix)) .and. (itmp(iw,ix) .gt. 0.0)) then
          iw2=iw
          exit
        end if
      enddo

      do iw=nw,1,-1
        if((w440 .ge. wtmp(iw,ix)) .and. (itmp(iw,ix) .gt. 0.0)) then
          iw1=iw
          exit
        end if
      enddo

      yy1=itmp(iw1,ix)
      yy2=itmp(iw2,ix)
      ww1=w440-wtmp(iw1,ix)
      ww2=wtmp(iw2,ix)-w440
      irr_out_irradiance_440nm(ix)=(ww1*yy2+ww2*yy1)/(ww1+ww2)
    enddo

    !11111111111111111111111
  end subroutine read_irr
  !11111111111111111111111

  !*********************
end module m_read_irr
!*********************
