module m_cal_ocp
  public cal_ocp
contains
!******************
subroutine cal_ocp
  !******************
  ! -----------------------------
  ! define ProcessingQualityFlags
  ! -----------------------------
  ! bit??  meaning:WhereSet
  ! bit00  (Error) invalid geolocation: m_cal_ecf.f90
  ! bit01  (Error) SZA,VZA,RAA out of LUT range: m_cal_ecf.f90, m_cal_ocp.f90 
  ! bit02  (Warning) ecf < minECF (0.05): m_cal_ocp.f90 
  ! bit03  (ERROR) input surface pressure or albedo error: m_cal_ecf.f90
  ! bit04  (Warning) snow_ice_fraction > min_snowice
  ! bit05  (Warning) SCD correction max_scd_iter reached in ocp : m_cal_ocp.f90
  ! bit06  (Error) SCD < 0, : m_cal_ocp.f90
  ! bit07  (Warning) 440nm radiance or irradiance error: m_read_input_tio.f90
  ! bit08  (ERROR) 466nm radiance or irradiance error: m_read_input_tio.f90
  ! bit09  (Warning) ecf out of normal range, clipped: m_cal_ecf.g90
  ! bit10  SceneAlbedoAtTerrain.eq.'yes' // 'both' skipped, 
  !        SCD correction proble or max_scd_iter reached
  ! bit11  SceneAlbedoAtTerrain.eq.'np' // 'both' skipped,
  !        SCD correction problem or max_scd_iter reeached
  ! bit12  (ERROR) skipped cloud ecf calculation 
  !        due to any problem during processing: m_cal_ecf.f90
  ! bit13  (ERROR) skipped cloud ocp calculation
  !        due to any problem during processing: m_cal_ocp.f90
  ! bit14  (Warning) ocp out of normal range, clipped: m_cal_ocp.f90
  ! bit15  (Warning) skipped pscene calculation during processing

  use m_vars
  use m_read_GMI
  !use m_read_DEM
  use m_read_hdf5
  use m_scd_adjust

  implicit none

  real:: cal_ecf,cal_crf
  integer::ialb, isza, ivza, iraa, ipsfc, ipcld
  integer::ialb1,isza1,ivza1,iraa1,ipsfc1
  integer::ialb2,isza2,ivza2,iraa2,ipsfc2
  real::   walb1,wsza1,wvza1,wraa1,wpsfc1,vpsfc1,apsfc1
  real::   walb2,wsza2,wvza2,wraa2,wpsfc2,vpsfc2,apsfc2
  real::   alb440

  real::yy1,yy2,ww1,ww2
  real(kind=8)::cpp, t8p
  real :: temp_cpp, temp_t8p, delta_temp
  integer::ipsfc0,ipm0,ipm1,ipm2
  integer(kind=4)::iflag

  integer(kind=4)::ierr
  integer(kind=4)::nt,nx
  integer(kind=4)::it,ix, iternum

  integer(kind=4)::ip

  integer(kind=4)::gmi_ix1,gmi_ix2,gmi_iy1,gmi_iy2
  real::gmi_wx1,gmi_wx2,gmi_wy1,gmi_wy2
  real::pp11,pp12,pp21,pp22,pp1,pp2
  real::tt11,tt12,tt21,tt22,tt1,tt2
!hqw changed tt,pp to allocatable and removed pp_geos and tt_geos
!  real(kind=4),dimension(gmi_np)::tt
!  real(kind=4),dimension(gmi_np+1)::pp !include Psfc
!  real(kind=4),dimension(geos_np)::tt_geos
!  real(kind=4),dimension(geos_np+1)::pp_geos !include Psfc
  real(kind=4), dimension(:), allocatable :: tt, pp

  integer(kind=4)::kleipool_ix,kleipool_iy

  integer ::isnowice

!hqw move pflag00, pflag01 from m_vars.f90 to local variable
  integer(kind=4):: pflag00, pflag01

  real::a1111,a1112,a1121,a1122,a1211,a1212,a1221,a1222,a2111,a2112,a2121,a2122,a2211,a2212,a2221,a2222
  real::a111,a112,a121,a122,a211,a212,a221,a222
  real::a11,a12,a21,a22
  real::a1,a2

  real::vpsfc0,apsfc0, scdm, scdadj, scdmorg
  real,dimension(npcld)::amfvcd_int,amfvcd_ext

  real::diff,diff_save
  real::x0,x1,x2,xx
  real::y0,y1,y2,yy
  integer(kind=4)::ipp
  integer::option_psfc_clear

  real::pi,dtor

  ! hqw comments out DEM and BDEM related stuff
  ! ------
  ! local useful variables 
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.

  nt=rad_NumTimes
  nx=rad_nXtrack

! allocate m_vars variables
!hqw disabled STD arrays, they are not actually calculated
! allocate dimensions for outputs
!  allocate(out_CloudPressureSTD(nx,nt),stat=ierr)
!  out_CloudPressureSTD=int(iFillValue, kind=2)
!  allocate(out_TerrainPressureStdDev(nx,nt),stat=ierr)
!  out_TerrainPressureStdDev=fFillValue
!hqw out_TerrainHeight now moved to m_read_input_tio.f90
!  allocate(out_TerrainHeight(nx,nt),stat=ierr)
!  out_TerrainHeight=fFillValue
!  allocate(out_TerrainHeightStdDev(nx,nt),stat=ierr)
!  out_TerrainHeightStdDev=fFillValue
!hqw LandAreaFraction is not actually used
!  allocate(out_LandAreaFraction(nx,nt),stat=ierr)
!  out_LandAreaFraction=int(iFillValue, kind=2)

  allocate(out_CloudPressure(nx,nt),stat=ierr)
  allocate(out_CloudPressureNotClipped(nx,nt),stat=ierr)
  allocate(out_TerrainPressure(nx,nt),stat=ierr)
  allocate(out_SurfaceReflectivity440(nx,nt),stat=ierr)
  allocate(out_SurfaceReflectivity466(nx,nt),stat=ierr)
  allocate(out_SlantColumnAmountO2O2(nx,nt),stat=ierr)
  allocate(out_O2O2CloudTemperature(nx,nt),stat=ierr)

!hqw initialize to (negative) fill value
  out_CloudPressure=int(iFillValue, kind=2)
  out_CloudPressureNotClipped=int(iFillValue, kind=2)
  out_TerrainPressure=fFillValue ! surface pressure used in ocp calculation 
  out_SurfaceReflectivity440=fFillValue
  out_SurfaceReflectivity466=fFillValue
  out_SlantColumnAmountO2O2=fFillValue
  out_O2O2CloudTemperature=fFillValue

  !hqw allocate and initialize local array
  allocate(tt(nlayers),stat=ierr)
  allocate(pp(nlayers+1),stat=ierr)

!hqw debug
  !write(*,*) 'writing debug_scd_adjust.txt'
  !open(unit=19,file='debug_scd_adjust.txt')
  !write(19,*)'    ix   scdmorg    scdm     scdadj      temp_t8p     t8p    temp_cpp'

  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========
      ! initialize cloud pressure to fFillValue
      cpp=fFillValue

      ! initialize for next iteration
      scdadj = fFillValue
      temp_t8p = fFillValue
      scdm = fFillValue
      t8p = fFillValue
      temp_cpp = fFillValue
      scdmorg = fFillValue

      ! initialize iflag to -1
      iflag=-1

      !---------------------
      ! geolocation check
      !--------------------
      pflag00=0
      if(abs(rad_Latitude(ix,it)) .gt. 90.) pflag00=pflag00+1
      if(abs(rad_longitude(ix,it)) .gt. 180.) pflag00=pflag00+1

      ! --------------------
      ! input angles check
      ! --------------------
      pflag01=0
      sza0=rad_SolarZenithAngle(ix,it)
      vza0=rad_ViewingZenithAngle(ix,it)
      raa0=out_RelativeAzimuthAngle(ix,it)
      !hqw changed 86. to max_SZA
      if((sza0 .lt. 0.) .or. (sza0 .gt. max_SZA)) pflag01=pflag01+1
      !hqw changed 72. to max_VZA
      if((vza0 .lt. 0.) .or. (vza0 .gt. max_VZA)) pflag01=pflag01+1
      !hqw skip invalid RAA
      !valid RAA is converted to [0.,180) range in m_read_input_tio.f90
      !invalid RAA is set to fspecial < 0. there
      if (raa0 .lt. 0.) pflag01=pflag01+1
      ! skip ocp calculation if any angle is invalid 
      if((pflag00 .ge. 1) .or. (pflag01 .ge. 1)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         ! bit 1 is also set in m_cal_ecf.f90, but only for SZA&VZA
         ! thus the following provides a safeguard for RAA
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),1)
         go to 990
      endif

      ! ----------------------------
      ! trim cloud fraction, fc & fr
      ! ----------------------------
      ! EffectiveCloudFraction and CloudRadianceFraction are clipped within
      ! [0.,1.), negative values signal bad data
      cal_ecf=out_EffectiveCloudFraction(ix,it)
      cal_crf=out_CloudRadianceFraction466(ix,it)

      !hqw skip calculation if cal_ecf or cal_crf are bad
      if ((cal_ecf .lt. 0.) .or. (cal_crf .lt. 0.)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)  
         go to 990
      endif

      !hqw temporaryalb0 and psfc0 which will be replaced by climatology
      alb0 = 0.05
      psfc0 = 1000.

      ! initialize local array
      tt = fFillValue
      pp = fFillValue
      vvcd = fFillValue

      ! ----------------------------------------------
      ! option for TemperaturePressure/SurfacePressure
      ! ----------------------------------------------
      if((name_option_TemperaturePressure.eq.'GMI')) then !.or. &
      !     (name_option_TemperaturePressure.eq.'DEM').or. &
      !     (name_option_TemperaturePressure.eq.'BDEM')) then
        gmi_ix1=floor((rad_Longitude(ix,it)+180.0)/1.25)+1
        gmi_ix2=gmi_ix1+1

        gmi_iy1=floor(rad_Latitude(ix,it)+90.)+1
        gmi_iy2=gmi_iy1+1

        if(gmi_ix1.lt.1) gmi_ix1=1
        if(gmi_ix1.gt.gmi_nx) gmi_ix1=gmi_nx
        if(gmi_ix2.lt.1) gmi_ix2=1
        if(gmi_ix2.gt.gmi_nx) gmi_ix2=gmi_nx
        if(gmi_iy1.lt.1) gmi_iy1=1
        if(gmi_iy1.gt.gmi_ny) gmi_iy1=gmi_ny
        if(gmi_iy2.lt.1) gmi_iy2=1
        if(gmi_iy2.gt.gmi_ny) gmi_iy2=gmi_ny

        gmi_wx1=rad_Longitude(ix,it)-gmi_lon(gmi_ix1)
        gmi_wx2=gmi_lon(gmi_ix2)-rad_Longitude(ix,it)

        gmi_wy1=rad_Latitude(ix,it)-gmi_lat(gmi_iy1)
        gmi_wy2=gmi_lat(gmi_iy2)-rad_Latitude(ix,it)
      endif

      if(name_option_TemperaturePressure.eq.'GMI') then
        pp11=gmi_TerrainPressure(gmi_ix1,gmi_iy1)
        pp12=gmi_TerrainPressure(gmi_ix1,gmi_iy2)
        pp21=gmi_TerrainPressure(gmi_ix2,gmi_iy1)
        pp22=gmi_TerrainPressure(gmi_ix2,gmi_iy2)
        pp1=(gmi_wy2*pp11+gmi_wy1*pp12)/(gmi_wy1+gmi_wy2)
        pp2=(gmi_wy2*pp21+gmi_wy1*pp22)/(gmi_wy1+gmi_wy2)
        gmi_psfc=(gmi_wx2*pp1+gmi_wx1*pp2)/(gmi_wx1+gmi_wx2)
        psfc0=gmi_psfc
        if (nlayers .NE. gmi_np) then
           write(*,*)'gmi_np,nlayers incompatible',gmi_np,nlayers
           call exit(-1)
        endif
        do ip=1,gmi_np
          pp11=gmi_Pressure(gmi_ix1,gmi_iy1,ip)
          pp12=gmi_Pressure(gmi_ix1,gmi_iy2,ip)
          pp21=gmi_Pressure(gmi_ix2,gmi_iy1,ip)
          pp22=gmi_Pressure(gmi_ix2,gmi_iy2,ip)
          pp1=(gmi_wy2*pp11+gmi_wy1*pp12)/(gmi_wy1+gmi_wy2)
          pp2=(gmi_wy2*pp21+gmi_wy1*pp22)/(gmi_wy1+gmi_wy2)
          pp(ip)=(gmi_wx2*pp1+gmi_wx1*pp2)/(gmi_wx1+gmi_wx2)
          tt11=gmi_temperature(gmi_ix1,gmi_iy1,ip)
          tt12=gmi_temperature(gmi_ix1,gmi_iy2,ip)
          tt21=gmi_temperature(gmi_ix2,gmi_iy1,ip)
          tt22=gmi_temperature(gmi_ix2,gmi_iy2,ip)
          tt1=(gmi_wy2*tt11+gmi_wy1*tt12)/(gmi_wy1+gmi_wy2)
          tt2=(gmi_wy2*tt21+gmi_wy1*tt22)/(gmi_wy1+gmi_wy2)
          tt(ip)=(gmi_wx2*tt1+gmi_wx1*tt2)/(gmi_wx1+gmi_wx2)
        end do
        pp(gmi_np+1)=gmi_psfc
        call read_GMI_VCD(pp,tt)
        vvcd=gmi_vcd
      endif

     !This option should not be used for TEMPO without development
     ! if(name_option_TemperaturePressure.eq.'DEM') then
     !!   psfc0=real(inp_TerrainPressure(ix,it))
     !   do ip=1,gmi_np
     !     tt11=gmi_temperature(gmi_ix1,gmi_iy1,ip)
     !     tt12=gmi_temperature(gmi_ix1,gmi_iy2,ip)
     !     tt21=gmi_temperature(gmi_ix2,gmi_iy1,ip)
     !     tt22=gmi_temperature(gmi_ix2,gmi_iy2,ip)
     !     tt1=(gmi_wy2*tt11+gmi_wy1*tt12)/(gmi_wy1+gmi_wy2)
     !     tt2=(gmi_wy2*tt21+gmi_wy1*tt22)/(gmi_wy1+gmi_wy2)
     !     tt(ip)=(gmi_wx2*tt1+gmi_wx1*tt2)/(gmi_wx1+gmi_wx2)
     !   end do
     !   call read_DEM_VCD(psfc0,tt,pp)
     !   vvcd=dem_vcd
     ! endif

     ! if(name_option_TemperaturePressure.eq.'BDEM') then
     !   psfc0=BDEM_TerrainPressure(ix,it)
     !   do ip=1,gmi_np
     !     tt11=gmi_temperature(gmi_ix1,gmi_iy1,ip)
     !     tt12=gmi_temperature(gmi_ix1,gmi_iy2,ip)
     !     tt21=gmi_temperature(gmi_ix2,gmi_iy1,ip)
     !     tt22=gmi_temperature(gmi_ix2,gmi_iy2,ip)
     !     tt1=(gmi_wy2*tt11+gmi_wy1*tt12)/(gmi_wy1+gmi_wy2)
     !     tt2=(gmi_wy2*tt21+gmi_wy1*tt22)/(gmi_wy1+gmi_wy2)
     !     tt(ip)=(gmi_wx2*tt1+gmi_wx1*tt2)/(gmi_wx1+gmi_wx2)
     !   end do
     !   call read_DEM_VCD(psfc0,tt,pp)
     !   vvcd=dem_vcd
     ! endif

      if(name_option_TemperaturePressure.eq.'GEOS5') then
        psfc0=l2_TerrainPressure(ix,it)
        if (nlayers .NE. geos_np) then
           write(*,*)'ERROR: geos_np nlayers mismatch',geos_np,nlayers
           call exit(-1)
        endif
        do ip=1,geos_np
          pp(ip)=geos_Pressure(ix,it,ip)
          tt(ip)=geos_temperature(ix,it,ip)
        end do
        pp(geos_np+1) = psfc0
        call read_GEOS5_VCD(pp,tt)
        vvcd=geos_vcd
      endif

!hqw debug
      if ((it .eq. itdebug).AND. (ix .eq. ixdebug)) then
         write(*,*) 'writing debug_tp.txt'
         open(unit=20,file='debug_tp.txt')
         write(20,*) trim(name_option_TemperaturePressure)
         write(20,*) 'Level, Pressure(hPa), Temperature(K)'
         write(20,*)'ix, it=',ix,it
         write(20,*)'latitude=',rad_latitude(ix,it)
         write(20,*)'longitude=',rad_longitude(ix,it)
         write(20,*) 'psfc=',pp(nlayers+1)
         do ip = 1, nlayers
            write(20,*)ip, pp(ip), tt(ip)
         end do
         close(20)
      endif

      ! ------------------------------
      ! option for SurfaceReflectivity
      ! ------------------------------
      if(name_option_SurfaceReflectivity.eq.'Kleipool') then
        kleipool_ix=nint((rad_Longitude(ix,it)+180.0)/0.5)
        kleipool_iy=nint((rad_Latitude(ix,it)+90.0)/0.5)
        if(kleipool_ix.lt.1) kleipool_ix=1
        if(kleipool_ix.gt.kleipool_nx) kleipool_ix=kleipool_nx
        if(kleipool_iy.lt.1) kleipool_iy=1
        if(kleipool_iy.gt.kleipool_ny) kleipool_iy=kleipool_ny
        alb0=kleipool_SurfaceReflectivity466(kleipool_ix,kleipool_iy)
        alb440=kleipool_SurfaceReflectivity440(kleipool_ix,kleipool_iy)
      endif

      if(name_option_SurfaceReflectivity.eq.'BRDF') then
        alb0=BRDF_SurfaceReflectivity466(ix,it)
        !hqw current BRDF tables does not have 440, thus
        !  all BRDF_SurfaceReflectivity440 are filled with fspecial
        alb440=BRDF_SurfaceReflectivity440(ix,it)
      endif

      ! clip and assign out_SurfaceReflectivity which contains 
      ! input LER or GLER depending on name_option
      if((alb0 .ge. -2.0) .and. (alb0 .lt. 0.0)) alb0=0.0
      if((alb0 .gt.  1.0) .and. (alb0 .le. 2.0)) alb0=1.0
      out_SurfaceReflectivity466(ix,it)=alb0

      if((alb440 .ge. -2.0) .and. (alb440 .lt. 0.0)) alb440=0.0
      if((alb440 .gt.  1.0) .and. (alb440 .le. 2.0)) alb440=1.0
      out_SurfaceReflectivity440(ix,it)=alb440

      ! -----------------------------------------------------
      ! option for SnowIce:
      !   - 'Pscene': calculate Pscene
      !   - 'Pcld': calculate Pcld over snow/ice
      ! -----------------------------------------------------
      isnowice=0

      if (rad_SnowIceFraction(ix,it) .gt. min_snowice) isnowice = 1

      ! -------------------------------------------
      ! set ProcessingQualityFlags:
      ! -------------------------------------------

      if(isnowice .gt. 0) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),4)
      end if

      if((cal_ecf .gt. 0.).and.(cal_ecf .lt. min_ecf)) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),2)
      endif

      if(nasa_SlantColumnAmountO2O2(ix,it).lt.0.0) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),6)
      endif

      ! skip calcultion if these bits are set
      if(btest(out_ProcessingQualityFlags(ix,it),0).or. & ! geolocation error
           btest(out_ProcessingQualityFlags(ix,it),1).or. & ! SZA//VZA//RAA error
           !hqw do not skip ocp calculation when 0<ecf<min_ecf  
           !btest(out_ProcessingQualityFlags(ix,it),2).or. & ! < min_ecf
           !hqw do not skip snow/ice scene for ocp calculation
           !btest(out_ProcessingQualityFlags(ix,it),4).or. & !snowice
           btest(out_ProcessingQualityFlags(ix,it),6).or. & ! scd<0
           !hqw skip ocp calculation if effective cloud fraction is invalid (<0.)
           btest(out_ProcessingQualityFlags(ix,it),12)) then
            out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13) 
            go to 990 
      endif

      ! skip calculation if cal_ecf are out of range or ecf=0.0
      ! NOTE: cloud pressure is skipped when ecf=0.0
      ! when there is no cloud there is no need to calculate cloud pressure 
      if((cal_ecf.le.0.0).or.(cal_ecf.gt.1.0)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13) 
         go to 990
      endif

      !hqw moved this section from before to after check ProcessingQualityFlags
      ! as these calculation is not needed if we deice to skip calculation
      ! -----------------
      ! set nodes for LUT
      ! -----------------
      ialb1=-9; ialb2=-9
      walb1=0.; walb2=0.
      do ialb=1,nalb-1
        if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
          ialb1=ialb
          ialb2=ialb+1
          walb1=alb0-lut_alb(ialb)
          walb2=lut_alb(ialb+1)-alb0
        endif
      end do
      if(ialb1 .lt. 0) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      isza1=-9; isza2=-9
      wsza1=0.; wsza2=0.
      do isza=1,nsza-1
        if((sza0 .ge. lut_sza(isza)) .and. (sza0 .le. lut_sza(isza+1))) then
          isza1=isza
          isza2=isza+1
          wsza1=sza0-lut_sza(isza)
          wsza2=lut_sza(isza+1)-sza0
        endif
      end do
      if(isza1 .lt. 0) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      ivza1=-9; ivza2=-9
      wvza1=0.; wvza2=0.
      do ivza=1,nvza-1
        if((vza0 .ge. lut_vza(ivza)) .and. (vza0 .le. lut_vza(ivza+1))) then
          ivza1=ivza
          ivza2=ivza+1
          wvza1=vza0-lut_vza(ivza)
          wvza2=lut_vza(ivza+1)-vza0
        endif
      end do
      if(ivza1 .lt. 0) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      iraa1=-9; iraa2=-9
      wraa1=0.; wraa2=0.
      do iraa=1,nraa-1
        if((raa0 .ge. lut_raa(iraa)) .and. (raa0 .le. lut_raa(iraa+1))) then
          iraa1=iraa
          iraa2=iraa+1
          wraa1=raa0-lut_raa(iraa)
          wraa2=lut_raa(iraa+1)-raa0
        endif
      end do
      if(iraa1 .lt. 0) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
        go to 990
      endif

      ! -------------------------
      ! AMF at surface: clear sky
      ! -------------------------
      do ipsfc=1,npsfc
        a1111=lut_amf_clr(ialb1,isza1,ivza1,iraa1,ipsfc)
        a1112=lut_amf_clr(ialb1,isza1,ivza1,iraa2,ipsfc)
        a1121=lut_amf_clr(ialb1,isza1,ivza2,iraa1,ipsfc)
        a1122=lut_amf_clr(ialb1,isza1,ivza2,iraa2,ipsfc)
        a1211=lut_amf_clr(ialb1,isza2,ivza1,iraa1,ipsfc)
        a1212=lut_amf_clr(ialb1,isza2,ivza1,iraa2,ipsfc)
        a1221=lut_amf_clr(ialb1,isza2,ivza2,iraa1,ipsfc)
        a1222=lut_amf_clr(ialb1,isza2,ivza2,iraa2,ipsfc)
        a2111=lut_amf_clr(ialb2,isza1,ivza1,iraa1,ipsfc)
        a2112=lut_amf_clr(ialb2,isza1,ivza1,iraa2,ipsfc)
        a2121=lut_amf_clr(ialb2,isza1,ivza2,iraa1,ipsfc)
        a2122=lut_amf_clr(ialb2,isza1,ivza2,iraa2,ipsfc)
        a2211=lut_amf_clr(ialb2,isza2,ivza1,iraa1,ipsfc)
        a2212=lut_amf_clr(ialb2,isza2,ivza1,iraa2,ipsfc)
        a2221=lut_amf_clr(ialb2,isza2,ivza2,iraa1,ipsfc)
        a2222=lut_amf_clr(ialb2,isza2,ivza2,iraa2,ipsfc)

        a111=(wraa2*a1111+wraa1*a1112)/(wraa1+wraa2)
        a112=(wraa2*a1121+wraa1*a1122)/(wraa1+wraa2)
        a121=(wraa2*a1211+wraa1*a1212)/(wraa1+wraa2)
        a122=(wraa2*a1221+wraa1*a1222)/(wraa1+wraa2)
        a211=(wraa2*a2111+wraa1*a2112)/(wraa1+wraa2)
        a212=(wraa2*a2121+wraa1*a2122)/(wraa1+wraa2)
        a221=(wraa2*a2211+wraa1*a2212)/(wraa1+wraa2)
        a222=(wraa2*a2221+wraa1*a2222)/(wraa1+wraa2)

        a11=(wvza2*a111+wvza1*a112)/(wvza1+wvza2)
        a12=(wvza2*a121+wvza1*a122)/(wvza1+wvza2)
        a21=(wvza2*a211+wvza1*a212)/(wvza1+wvza2)
        a22=(wvza2*a221+wvza1*a222)/(wvza1+wvza2)

        a1=(wsza2*a11+wsza1*a12)/(wsza1+wsza2)
        a2=(wsza2*a21+wsza1*a22)/(wsza1+wsza2)

        cal_amf_clr(ipsfc)=(walb2*a1+walb1*a2)/(walb1+walb2)
      end do

      ! -----------
      ! AMF at Pcld
      ! -----------
      do ipcld=1,npcld
        ipsfc=ipcld
        a111=lut_amf_cld(isza1,ivza1,iraa1,ipcld,ipsfc)
        a112=lut_amf_cld(isza1,ivza1,iraa2,ipcld,ipsfc)
        a121=lut_amf_cld(isza1,ivza2,iraa1,ipcld,ipsfc)
        a122=lut_amf_cld(isza1,ivza2,iraa2,ipcld,ipsfc)
        a211=lut_amf_cld(isza2,ivza1,iraa1,ipcld,ipsfc)
        a212=lut_amf_cld(isza2,ivza1,iraa2,ipcld,ipsfc)
        a221=lut_amf_cld(isza2,ivza2,iraa1,ipcld,ipsfc)
        a222=lut_amf_cld(isza2,ivza2,iraa2,ipcld,ipsfc)

        a11=(wraa2*a111+wraa1*a112)/(wraa1+wraa2)
        a12=(wraa2*a121+wraa1*a122)/(wraa1+wraa2)
        a21=(wraa2*a211+wraa1*a212)/(wraa1+wraa2)
        a22=(wraa2*a221+wraa1*a222)/(wraa1+wraa2)

        a1=(wvza2*a11+wvza1*a12)/(wvza1+wvza2)
        a2=(wvza2*a21+wvza1*a22)/(wvza1+wvza2)

        cal_amf_cld(ipcld)=(wsza2*a1+wsza1*a2)/(wsza1+wsza2)
      end do

      ! -----------
      ! check psfc0
      ! -----------
      ! make sure psfc0 is <= lut_psfc(npsfc)
      if(psfc0 .gt. lut_psfc(npsfc)) psfc0=lut_psfc(npsfc)

      ipsfc0=-9
      ipsfc1=-9; ipsfc2=-9
      do ipsfc=1,npsfc-1
        if((psfc0 .gt. lut_psfc(ipsfc)) .and. (psfc0 .le. lut_psfc(ipsfc+1))) then
          ipsfc1=ipsfc
          ipsfc2=ipsfc+1
          wpsfc1=psfc0-lut_psfc(ipsfc)
          wpsfc2=lut_psfc(ipsfc+1)-psfc0
          vpsfc1=vvcd(ipsfc)
          vpsfc2=vvcd(ipsfc+1)
          apsfc1=real(cal_amf_clr(ipsfc), kind=4)
          apsfc2=real(cal_amf_clr(ipsfc+1), kind=4)
        endif
        if((psfc0 .ge. lut_psfc(ipsfc)) .and. (psfc0 .lt. lut_psfc(ipsfc+1))) then
          ipsfc0=ipsfc
        endif
      end do

      if(ipsfc1 .lt. 0) then
        write(*,*) " *** Pcld: Check Surface Pressure *** ",ix,it,psfc0
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
        go to 990
      else
        vpsfc0=(wpsfc1*vpsfc2+wpsfc2*vpsfc1)/(wpsfc1+wpsfc2)
        apsfc0=(wpsfc1*apsfc2+wpsfc2*apsfc1)/(wpsfc1+wpsfc2)
      endif

      !hqw added initizlization to ipm?
      !ipm? are used for extrapolation at the large pressure end
      ipm0=-9; ipm1=-9; ipm2=-9
      if(ipsfc0 .gt. 0) then
        ipm0=ipsfc0-0
        ipm1=ipsfc0-1
        ipm2=ipsfc0-2
      else
        ipm0=npsfc-0
        ipm1=npsfc-1
        ipm2=npsfc-2
      endif

      ! -----------------
      ! calculate AMF*VCD
      ! -----------------

      ! hqw changed to NASA_SlantColumnAmountO2O2 and moved up here
      ! Temperature correction will be applied through iteration
      !scdm=out_SlantColumnAmountO2O2(ix,it)
      scdmorg = nasa_SlantColumnAmountO2O2(ix,it)

      !hqw initial iteration use 273K reference
      iternum = 0
      t8p = 273. ! initial reference temperature for SCD retrieval
      temp_t8p = 273.
      scdm = scdmorg
      scdadj = scdmorg

!hqw moved 777 downward, because amfvcds do not depend on scdm
!777   continue

      do ipcld=1,npcld
        !ipsfc=ipcld !hqw seems unnecessary
        amfvcd_int(ipcld)=vvcd(ipcld)*cal_crf*real(cal_amf_cld(ipcld),kind=4) &
             +vpsfc0*(1.0-cal_crf)*apsfc0
        amfvcd_ext(ipcld)=vvcd(ipcld)*cal_crf*real(cal_amf_cld(ipcld), kind=4)&
             +vvcd(ipcld)*(1.0-cal_crf)*real(cal_amf_clr(ipcld), kind=4)
      end do

      !hqw looks like this hardcode Pclr=Psfc when Pcld>Psfc
      option_psfc_clear=0
      ! find pressure for AMF*VCD
      !    1: Pclr = Pcld if Pcld > Psfc
      !    0: Pclr = Psfc (fixed) & Pcld > Psfc
      ! ????????????????????????????????????

!hqw temperature iteration comes back here to 777
777   continue ! iteration

      iflag=-1

      if(scdm.le.amfvcd_int(5)) then
        iflag=0
        !    write(3,311) it,ix,cal_ecf,cal_crf,scdm,amfvcd_int(1),amfvcd_int(2),amfvcd_int(3),amfvcd_int(4),amfvcd_int(5)
        ! 311 format(2i3,3f10.4,5f10.4)
      endif

      yy1=0. ; yy2=0.
      ww1=0. ; ww2=0.
      do ipcld=1,npcld-1
        if((scdm.gt.amfvcd_int(ipcld)).and.(scdm.le.amfvcd_int(ipcld+1))) then
          iflag=1
          yy1=lut_pcld(ipcld)
          yy2=lut_pcld(ipcld+1)
          ww1=scdm-amfvcd_int(ipcld)
          ww2=amfvcd_int(ipcld+1)-scdm
        endif
      end do

      if(iflag .eq. 1) then ! normal interpolation
        cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)
      ! the choice below is for low pressure cases
      else if(iflag .eq. 0) then ! scdm<= amfvcd_int(5)
        x0=0.0
        x1=lut_pcld(5)
        x2=lut_pcld(6)
        y0=0.0
        y1=amfvcd_int(5)
        y2=amfvcd_int(6)

        xx=x1
        !hqw 1st & 3rd term below is always 0, yy=y1,is it needed?
        !but it makes the formula consistent with the one inside ipp loop
        yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
             +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
             +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
        diff_save=abs(scdm-yy)

        !hqw increment ipp by 1Pa at a time until min diff found
        ! 150 Pa > lut_pcld[5:6], thus a safe choice
        ! may need work if LUT is changed
        do ipp=1,150
          xx=x1-real(ipp)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff=abs(scdm-yy)
          if(diff.ge.diff_save) then
            go to 970
          else
            diff_save=diff
            if(ipp.le.1) then
              xx=-9999.
            endif
          endif
        end do
970     continue
        cpp=real(nint(xx))
        !    write(3,312) '0',it,ix,cal_ecf,cal_crf,scdm,x0,x1,x2,y0,y1,y2,cpp
        ! 312 format(a1,2x,2i3,3f10.4,7f10.4)
      else ! very large scdm case, ipm0 should tend to npsfc
        iflag=2
        !hqw added safeguard logic, though the program is expected to
        !actually bypass it because all ipms should be valid
        if ((ipm2 .lt. 1) .or. (ipm1 .lt. 1) .or. (ipm0 .lt. 1)) then
           write(*,*) 'ipm < 0 for large scdm, which should not happen.'
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
           go to 990
        endif
        x0=lut_pcld(ipm0)
        x1=lut_pcld(ipm1)
        x2=lut_pcld(ipm2)

        if(option_psfc_clear.eq.0) then !hqw current hardcoded choice
          y0=amfvcd_int(ipm0)
          y1=amfvcd_int(ipm1)
          y2=amfvcd_int(ipm2)
        endif
        if(option_psfc_clear.eq.1) then
          y0=amfvcd_ext(ipm0)
          y1=amfvcd_ext(ipm1)
          y2=amfvcd_ext(ipm2)
        endif

        xx=psfc0
        yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
             +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
             +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
        diff_save=abs(scdm-yy)

        !hqw seems this tries ipp one at a time until mininal difference is found
        ! 5000 is a large number which safely covers psfc, 2000 should be enough
        ! but should make no difference to the computer
        do ipp=1,5000
          xx=psfc0+real(ipp)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff=abs(scdm-yy)
          if(diff.ge.diff_save) then
            go to 980
          else
            diff_save=diff
            if(ipp.ge.5000) then !this should not happen
              xx=-9999. ! if it ever does, set to invalid value
            endif
          endif
        end do

980     continue
        cpp=real(nint(xx))
      endif

      !hqw adjust scd according to T at cpp
      ! use the temperature at temp_cpp when in range
      ! try using T at half of the pressure
      temp_cpp = real (cpp * 0.5, kind=4)
      if (temp_cpp .gt. 50. .and. temp_cpp .lt. 1200.) then
        if (name_option_TemperaturePressure .eq. 'GMI') then
          call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
        else if (name_option_TemperaturePressure .eq. 'GEOS5') then
          call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
        else
          temp_t8p = real(t8p, kind=4)
        endif
      else
         temp_t8p = real(t8p, kind=4) !this will terminate iteration below
      endif
      iternum = iternum + 1

      !hqw debug
      !if ((it .eq. itdebug) .and. (ix .eq. itdebug)) then
      !   write(19,*) iternum, scdm, scdadj, temp_cpp, temp_t8p
      !endif

      !hqw test if terminate temperature iteration
      delta_temp = real(abs(t8p - temp_t8p), kind=4)
      if (delta_temp .lt. dt_threshold) then
         goto 990 ! exit iteration
      endif

      if (iternum .lt. max_scd_iter) then 
         t8p = temp_t8p  !update t8p from previous step
         scdm = scdadj !update scdm from previous step
         goto 777  ! goto iteration start
      endif

990   continue

      !hqw debug
      !if (it .eq. 102) then
      !   write(19,*) ix, scdmorg, scdm, scdadj, temp_t8p, t8p, temp_cpp
      !endif

      !hqw scdm & t8p is the step right before final iteration
      !   scdadj & temp_t8p is the step right after final iterateion
      if (scdm .gt. 0.) then
         out_SlantColumnAmountO2O2(ix,it) = scdm ! scdadj
         out_O2O2CloudTemperature(ix,it) = real(t8p, kind=4) ! temp_t8p
      else !hqw skipped pixels will end up here
         out_SlantColumnAmountO2O2(ix,it) = nasa_SlantColumnAmountO2O2(ix,it)
         out_O2O2CloudTemperature(ix,it) = 273.
      endif

      !set out_ProcessingQualityFlags digit 5 for max_scd_iter
      if (iternum .eq. max_scd_iter) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),5)
      endif

      out_CloudPressure(ix,it)=nint(cpp, kind=2)
      out_CloudPressureNotClipped(ix,it)=nint(cpp, kind=2)
      if(cpp.lt. 0.) then
        out_CloudPressure(ix,it)=int(iFillValue, kind=2)
        out_CloudPressureNotClipped(ix,it)=int(iFillValue, kind=2)
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
      endif
      ! clip slightly out-of-range values and set out_ProcessingQualityFlags
      if((cpp.gt.0.).and.(cpp.le.100.)) then
          out_CloudPressure(ix,it)=100
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),14)
      endif
      if((cpp.gt.psfc0).and.(cpp.lt.1200.)) then
          out_CloudPressure(ix,it)=nint(psfc0, kind=2)
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),14)
      endif

      !hqw skip out_CloudPressureSTD as inp_CloudPressureSTD is not read
      !out_CloudPressureSTD(ix,it)=inp_CloudPressureSTD(ix,it)
      !if(out_CloudPressureSTD(ix,it).le.iFillValue) &
      !     out_CloudPressureSTD(ix,it)=int(iFillValue, kind=2)

      !hqw moved out_TerrainPressure here, as psfc0 is the actual surfp used
      out_TerrainPressure(ix,it) = psfc0

      !hqw if BDEM replace with out_ with BDEM_
      !if(name_option_TemperaturePressure.eq.'BDEM') then
      !  out_TerrainPressure(ix,it)=BDEM_TerrainPressure(ix,it)
      !  !out_TerrainPressureStdDev(ix,it)=BDEM_TerrainPressureStdDev(ix,it)
      !  out_TerrainHeight(ix,it)=BDEM_TerrainHeight(ix,it)
      !  !out_TerrainHeightStdDev(ix,it)=BDEM_TerrainHeightStdDev(ix,it)
      !  !hqw LandAreaFraction is not used,  comment out
      !  !out_LandAreaFraction(ix,it)=int(BDEM_LandAreaFraction(ix,it), kind=2)
      !endif

      !  if(iflag.eq.0) out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),2)
      !  if(iflag.eq.2) out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)

      !=====
!888   continue
    end do
  end do
  !=====

  !hqw debug
  !close(19)

  ! deallocate allocated local variables
  deallocate(pp, tt)

  !**********************
end subroutine cal_ocp
!**********************
end module m_cal_ocp
