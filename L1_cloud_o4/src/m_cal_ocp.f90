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
  ! bit00  (Error) invalid lat/lon/SZA/VZA/RAA: m_cal_ecf.f90
  ! bit01  N/A reserved
  ! bit02  (Warning) pcld replaced by pscene because ecf<minECF: m_cal_pscene.f90 
  ! bit03  (ERROR) input surface pressure or albedo error: m_cal_ecf.f90
  ! bit04  (Warning) pcld replaced by pscene as snow_ice_fraction>min_snowice: m_cal_pscene.f90
  ! bit05  (Warning) SCD correction problem or max_scd_iter reached in ocp : m_cal_ocp.f90
  ! bit06  (Error) SCD < 0 or SCD quality issue, : m_cal_ocp.f90
  ! bit07  (Warning) 440nm radiance or irradiance error: m_read_input_tio.f90
  ! bit08  (ERROR) 466nm radiance or irradiance error: m_read_input_tio.f90
  ! bit09  (ERROR) ecf out of normal range, clipped: m_cal_ecf.g90
  ! bit10  SceneAlbedoAtTerrain.eq.'yes' // 'both' skipped, 
  !        SCD correction problem or max_scd_iter reached
  ! bit11  SceneAlbedoAtTerrain.eq.'np' // 'both' skipped,
  !        SCD correction problem or max_scd_iter reeached
  ! bit12  (ERROR) skipped cloud ecf calculation 
  !        due to any problem during processing: m_cal_ecf.f90
  ! bit13  (ERROR) skipped cloud ocp calculation due to
  !        any problem during processing,or invalid ocp: m_cal_ocp.f90
  ! bit14  (ERROR) ocp out of normal range and clipped: m_cal_ocp.f90
  ! bit15  (Warning) skipped pscene calculation during processing

  use m_vars
  use m_read_GMI
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
  real(kind=4), dimension(:), allocatable :: tt, pp

  real::a1111,a1112,a1121,a1122,a1211,a1212,a1221,a1222,a2111,a2112,a2121,a2122,a2211,a2212,a2221,a2222
  real::a111,a112,a121,a122,a211,a212,a221,a222
  real::a11,a12,a21,a22
  real::a1,a2

  real::vpsfc0,apsfc0, scdm, scdadj, scdmorg
  real,dimension(npcld)::amfvcd_int,amfvcd_ext

  real::diff,diff_save, maxpress
  real::x0,x1,x2,xx
  real::y0,y1,y2,yy
  integer(kind=4)::ipp
!  integer::option_psfc_clear !hqw moved to m_vars

  real::pi,dtor

!hqw add local variables
  real:: lat0, lon0
  real:: fFillValue9

  ! ------
  ! local useful variables 
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.

  fFillValue9 = -9999.

  nt=rad_NumTimes
  nx=rad_nXtrack

  maxpress = 2000 !Pa

! allocate m_vars variables
!hqw disabled STD arrays, they are not actually calculated
! allocate dimensions for outputs
!  allocate(out_CloudPressureSTD(nx,nt),stat=ierr)
!  out_CloudPressureSTD=int(iFillValue, kind=2)
!  allocate(out_TerrainPressureStdDev(nx,nt),stat=ierr)
!  out_TerrainPressureStdDev=fFillValue9
!  allocate(out_TerrainHeightStdDev(nx,nt),stat=ierr)
!  out_TerrainHeightStdDev=fFillValue9

  allocate(out_CloudPressure(nx,nt),stat=ierr)
  allocate(out_CloudPressureNotClipped(nx,nt),stat=ierr)
!hqw moved out_TerrainPressure to m_cal_pscene
! repurpose to hold calculated cpp for name_optio_SceneAlbedoAtTerrain='yes' 
!  allocate(out_TerrainPressure(nx,nt),stat=ierr)
  allocate(out_SlantColumnAmountO2O2(nx,nt),stat=ierr)
  allocate(out_O2O2CloudTemperature(nx,nt),stat=ierr)

!hqw initialize to (negative) fill value
  out_CloudPressure=int(iFillValue, kind=2)
  out_CloudPressureNotClipped=int(iFillValue, kind=2)
!  out_TerrainPressure=fFillValue9 
  out_SlantColumnAmountO2O2=fFillValue9
  out_O2O2CloudTemperature=fFillValue9

  !hqw allocate and initialize local array
  allocate(tt(nlayers),stat=ierr)
  allocate(pp(nlayers+1),stat=ierr)
  tt = fFillValue9
  pp = fFillValue9

!hqw debug
  !write(*,*) 'writing debug_scd_adjust.txt'
  !open(unit=19,file='debug_scd_adjust.txt')
  !write(19,*)'    ix   scdmorg    scdm     scdadj      temp_t8p     t8p    temp_cpp'

  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========
      ! initialize local variables
      ! this is also the value if calculation skipped
      alb0 = fFillValue9
      alb440 = fFillValue9
      psfc0 = fFillValue9
      cpp=fFillValue9

      ! initialize local variable for T correction
      scdadj = fFillValue9
      temp_t8p = fFillValue9
      scdm = fFillValue9
      t8p = fFillValue9
      temp_cpp = fFillValue9
      scdmorg = fFillValue9

      !initialize local array
      cal_amf_clr = fFillValue
      cal_amf_cld = fFillValue

      ! get local angles and geolocation
      sza0=rad_SolarZenithAngle(ix,it)
      vza0=rad_ViewingZenithAngle(ix,it)
      raa0=out_RelativeAzimuthAngle(ix,it)
      lat0=rad_Latitude(ix,it)
      lon0=rad_Longitude(ix,it)

      ! initialize iflag, iternum 
      iflag=-1
      iternum = 0

!hqw geolocation and angles were checked in m_cal_ecf
!invalid values triggers bit0 to be set
!thus, check bit0 to decide whether to skip calculation
      if (btest(out_ProcessingQualityFlags(ix,it),0)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      ! ----------------------------
      ! cloud fraction
      ! ----------------------------
      ! ecf and crf are already clipped within [0.,1.) in m_cal_ecf
      ! negative values signal bad data
      cal_ecf=out_EffectiveCloudFraction(ix,it)
      cal_crf=out_CloudRadianceFraction466(ix,it)

      !hqw skip ocp if cal_ecf or cal_crf are bad or ZERO
      !Note when ecf//crf=0, there is no need to calculate ocp
      if ((cal_ecf .le. 0.) .or. (cal_crf .le. 0.) .or. (cal_ecf .gt. 1.)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)  
         go to 990
      endif

      !hqw skip if bit3 is set (psfc//rsfc error)
      if (btest(out_ProcessingQualityFlags(ix,it),3)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      !----------------------------------------------
      ! T-P profile and vvcd
      !----------------------------------------------
      ! initialize local array
      tt = fFillValue9
      pp = fFillValue9
      vvcd = fFillValue9

      ! ----------------------------------------------
      ! option for TemperaturePressure/SurfacePressure
      ! ----------------------------------------------
      ! bad psfc should have been skipped before
      if((name_option_TemperaturePressure.eq.'GMI')) then 
        gmi_ix1=floor((lon0+180.0)/1.25)+1
        gmi_ix2=gmi_ix1+1
        gmi_iy1=floor(lat0+90.)+1
        gmi_iy2=gmi_iy1+1

        if(gmi_ix1.lt.1) gmi_ix1=1
        if(gmi_ix1.gt.gmi_nx) gmi_ix1=gmi_nx
        if(gmi_ix2.lt.1) gmi_ix2=1
        if(gmi_ix2.gt.gmi_nx) gmi_ix2=gmi_nx
        if(gmi_iy1.lt.1) gmi_iy1=1
        if(gmi_iy1.gt.gmi_ny) gmi_iy1=gmi_ny
        if(gmi_iy2.lt.1) gmi_iy2=1
        if(gmi_iy2.gt.gmi_ny) gmi_iy2=gmi_ny

        gmi_wx1=lon0-gmi_lon(gmi_ix1)
        gmi_wx2=gmi_lon(gmi_ix2)-lon0

        gmi_wy1=lat0-gmi_lat(gmi_iy1)
        gmi_wy2=gmi_lat(gmi_iy2)-lat0
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
        !hqw adds safeguard
        if ((psfc0 .lt. lut_psfc(1)).or.(psfc0.gt.lut_psfc(npsfc))) then 
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
          go to 990
        endif
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

      if(name_option_TemperaturePressure.eq.'GEOS5') then
        psfc0=l2_TerrainPressure(ix,it)
        !hqw safeguard
        if ((psfc0.lt.lut_psfc(1)).or.(psfc0.gt.lut_psfc(npsfc))) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
          go to 990
        endif
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

     !---------------------------------
     ! surface reflectivity
     !---------------------------------
!hqw now directly use out_SurfaceReflectivity assigned in ecf
!   instead of repeating calculation

      ! skip if bit 3 (psfc or rsfc error) is set
      if (btest(out_ProcessingQualityFlags(ix,it),3)) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
        go to 990
      endif

      alb0=out_SurfaceReflectivity466(ix,it)
      alb440=out_SurfaceReflectivity440(ix,it)

      ! skip if alb0 is negative, safeguard
      if (alb0 .lt. 0.) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      !hqw in read_cldo4_tio, negative or bad SCD are set to fspecial=-9999.,
      if(nasa_SlantColumnAmountO2O2(ix,it).lt.0.0) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),6)
        go to 990
      endif

      ! skip calcultion if these bits are set
      !if(btest(out_ProcessingQualityFlags(ix,it),0).or. & ! already checked
           !hqw do not skip ocp calculation when 0<ecf<min_ecf  
           !btest(out_ProcessingQualityFlags(ix,it),2).or. & ! ecf< min_ecf
           ! bit2 (min_ecf)is now handled in pscene
           !hqw do not skip snow/ice scene for ocp calculation
           !btest(out_ProcessingQualityFlags(ix,it),4).or. & ! snowice
           ! snowice is now handled in pscene
           !btest(out_ProcessingQualityFlags(ix,it),6).or. & ! checked before

      !hqw skip ocp if effective cloud fraction is skipped 
      if (btest(out_ProcessingQualityFlags(ix,it),12)) then
            out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13) 
            go to 990 
      endif

      !hqw moved this section from before 
      ! as these are not needed if we decide to skip 
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
      !hqw in LUT when ipcld>ipsfc, lut_amf_cld<0.
      ! but this would not happen below, as ipsfc=ipcld 
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
      ! make sure psfc0 is <= lut_psfc(npsfc),safeguard
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

      !hqw added initizlization to ipms
      !ipms are used for extrapolation at high pressure end
      ipm0=-9; ipm1=-9; ipm2=-9
      if(ipsfc0 .gt. 0) then ! psfc0 < lut_psfc(npsfc)
        ipm0=ipsfc0-0
        ipm1=ipsfc0-1
        ipm2=ipsfc0-2
      else ! psfc0>=lut_psfc(npsfc)
        ipm0=npsfc-0
        ipm1=npsfc-1
        ipm2=npsfc-2
      endif

      ! -----------------
      ! calculate AMF*VCD
      ! -----------------

      ! Temperature correction will be applied through iteration
      scdmorg = nasa_SlantColumnAmountO2O2(ix,it)

      !hqw initial iteration use 273K reference
      iternum = 0
      t8p = 273. ! initial reference temperature for SCD retrieval
      temp_t8p = 273.
      scdm = scdmorg
      scdadj = scdmorg

      !hqw amfvcd_int uses psfc0, amfvcd_ext uses Pcld
      !hqw amfvcd_ext > amfvcd_int when pcld > psfc0
      do ipcld=1,npcld
        !ipsfc=ipcld !hqw seems unnecessary
        amfvcd_int(ipcld)=vvcd(ipcld)*cal_crf*real(cal_amf_cld(ipcld),kind=4) &
             +vpsfc0*(1.0-cal_crf)*apsfc0
        amfvcd_ext(ipcld)=vvcd(ipcld)*cal_crf*real(cal_amf_cld(ipcld), kind=4)&
             +vvcd(ipcld)*(1.0-cal_crf)*real(cal_amf_clr(ipcld), kind=4)
      end do

      !????????????????????????????????????????
      ! hqw move option_psfc_clear to m_vars.f90
      !hqw looks like this hardcodes Pclr=Psfc
      !option_psfc_clear=0
      ! find pressure for AMF*VCD
      !    0: Pclr = Psfc (fixed) & Pcld > Psfc
      !    1: Pclr = Pcld if Pcld > Psfc
      ! ???????????????????????????????????????

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!hqw SCD iteration comes back here to 777
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
777   continue ! iteration

      iflag=-1

!hqw pcld(5)=104Pa, pcld(6)=121Pa
      if(scdm.le.amfvcd_int(5)) then ! low pressure end
        iflag=0
        !    write(3,311) it,ix,cal_ecf,cal_crf,scdm,amfvcd_int(1),amfvcd_int(2),amfvcd_int(3),amfvcd_int(4),amfvcd_int(5)
        ! 311 format(2i3,3f10.4,5f10.4)
      endif

      yy1=-9. ; yy2=-9.
      ww1=0. ; ww2=0.
      do ipcld=1,npcld-1
        if((scdm.gt.amfvcd_int(ipcld)).and.(scdm.le.amfvcd_int(ipcld+1))) then
          iflag=1 ! node found
          yy1=lut_pcld(ipcld)
          yy2=lut_pcld(ipcld+1)
          ww1=scdm-amfvcd_int(ipcld)
          ww2=amfvcd_int(ipcld+1)-scdm
        endif
      end do

      if(iflag .eq. 1) then ! normal interpolation
        cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)

      ! the choice below is for low pressure end
      else if(iflag .eq. 0) then ! scdm<= amfvcd_int(5)
        x0=0.0
        x1=lut_pcld(5)
        x2=lut_pcld(6)
        y0=0.0
        y1=amfvcd_int(5)
        y2=amfvcd_int(6)

        xx=x1
        !hqw 1st & 3rd term below is always 0, yy=y1,it is verbose, but
        !it makes the formula consistent with the one inside ipp loop later
        yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
             +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
             +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
        diff_save=abs(scdm-yy)

        !hqw increase ipp (dec xx) 1Pa at a time until min diff found
        ! 150 Pa > lut_pcld[5:6], thus it is safe 
        ! may need change if LUT is changed
        do ipp=1,150
          xx=x1-real(ipp)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff=abs(scdm-yy)
          if(diff.ge.diff_save) then
            go to 970 !inflection point found
          else
            diff_save=diff
            !hqw logic change
            !if(ipp.le.1) then
            if (xx .le. 0.)then
              xx=-9999.
              go to 970 !hqw addition, exit when no solution found
            endif
          endif
        end do
970     continue
        cpp=real(nint(xx))
        !    write(3,312) '0',it,ix,cal_ecf,cal_crf,scdm,x0,x1,x2,y0,y1,y2,cpp
        ! 312 format(a1,2x,2i3,3f10.4,7f10.4)

      else ! large scdm case: scdm>amfvcd_int(npcld) 
        iflag=2 
        !hqw added safeguard, though the program is expected to
        !bypass it because all ipms should be valid
        if ((ipm2 .lt. 1) .or. (ipm1 .lt. 1) .or. (ipm0 .lt. 1)) then
           write(*,*) 'ipm <= 0 for large scdm, this should not happen.'
           cpp = fFillValue9
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
           go to 990
        endif
        x0=lut_pcld(ipm0)
        x1=lut_pcld(ipm1)
        x2=lut_pcld(ipm2)

        if(option_psfc_clear.eq.0) then !original hardcoded choice
        !hqw use pclr=psfc,ipm0=ipsfc0
          y0=amfvcd_int(ipm0)
          y1=amfvcd_int(ipm1)
          y2=amfvcd_int(ipm2)
        endif
        if(option_psfc_clear.eq.1) then
        !hqw use pcld=pcld,ipm0=npsfc
          y0=amfvcd_ext(ipm0)
          y1=amfvcd_ext(ipm1)
          y2=amfvcd_ext(ipm2)
        endif

        xx=psfc0
        yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
             +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
             +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
        diff_save=abs(scdm-yy)

        !hqw this increase xx 1Pa at a time from psfc0
        !until mininal difference is found
        ! 5000 is large and safe, 2000 should be enough
        ! but should make no difference to the computer
        do ipp=1,5000
          xx=psfc0+real(ipp)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff=abs(scdm-yy)
          if(diff.ge.diff_save) then
            go to 980 !inflection point found 
          else
            diff_save=diff
            if (xx .gt. maxpress) then !no solution found
              xx=-9999. ! if it ever does, set to invalid value
              go to 980
            endif
          endif
        end do
980     continue
        cpp=real(nint(xx))
      endif

      !skip if cpp <0.
      if (cpp .lt. 0.) then 
          cpp = fFillValue9
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
          go to 990
      endif

      !hqw add clip cpp to within LUT range, safeguard
      if (cpp .gt. lut_psfc(npsfc)) then
          cpp = lut_psfc(npsfc)
      endif

      !hqw adjust scd according to T at temp_cpp
      ! use the temperature at temp_cpp when in range
      ! try using T at half of the pressure
      temp_cpp = real (cpp * 0.5, kind=4)
      if ((temp_cpp .ge. pp(1)) .and. (temp_cpp .le. psfc0)) then
        if (name_option_TemperaturePressure .eq. 'GMI') then
          call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
        else if (name_option_TemperaturePressure .eq. 'GEOS5') then
      ! scdmorg<0. should already been skipped
      ! returned scdadj always > 0., because if negative
      ! scdadj = scdmorg and temp_t8p=273K
          call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
        else
          temp_t8p = real(t8p, kind=4)
        endif
      else
         temp_t8p = real(t8p, kind=4) !this will terminate iteration below
         ! signal iteration problem
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),5)
      endif

      ! increment iternum
      iternum = iternum + 1

      !hqw debug
      !if ((it .eq. itdebug) .and. (ix .eq. itdebug)) then
      !   write(19,*) iternum, scdm, scdadj, temp_cpp, temp_t8p
      !endif

      !hqw test if terminate temperature iteration
      !technically, temp_t8p can be -999. when scdmorg<0.
      !but it won't happen as they should have been skipped
      !if negative scdadj occurs within scd_adjust_geos,
      ! scdadj=scdmorg and temp_t8p=273. 
      delta_temp = real(abs(t8p - temp_t8p), kind=4)
      if ((delta_temp .lt. dt_threshold).or.(temp_t8p .eq. 273.)) then
         goto 990 ! exit iteration
      endif

      if (iternum .lt. max_scd_iter) then 
         t8p = temp_t8p  !update t8p from previous step
         scdm = scdadj !update scdm from previous step
         goto 777  ! goto iteration start
      endif

! skipped calculation will end up here
990   continue

      !set out_ProcessingQualityFlags bit 5 for max_scd_iter
      if (iternum .eq. max_scd_iter) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),5)
      endif

      !hqw debug
      !if (it .eq. 102) then
      !   write(19,*) ix, scdmorg, scdm, scdadj, temp_t8p, t8p, temp_cpp
      !endif

      !hqw scdm & t8p is the step right before final iteration
      !   scdadj & temp_t8p is the step right after final iterateion
      if (scdm .gt. 0.) then
         out_SlantColumnAmountO2O2(ix,it) = scdm ! scdadj
         out_O2O2CloudTemperature(ix,it) = real(t8p, kind=4) ! temp_t8p
      else !hqw skipped pixels will satisfy this condition
         out_SlantColumnAmountO2O2(ix,it) = nasa_SlantColumnAmountO2O2(ix,it)
         out_O2O2CloudTemperature(ix,it) = 273.
      endif

      out_CloudPressure(ix,it)=nint(cpp, kind=2)
      out_CloudPressureNotClipped(ix,it)=nint(cpp, kind=2)
      if((cpp.le. 0.).or.(cpp.ge.lut_psfc(npsfc))) then
        out_CloudPressure(ix,it)=int(iFillValue, kind=2)
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
      endif

       if (option_clip_pcld .eq. 'yes') then
          if ((cpp .gt. psfc0).and.(cpp.le.lut_psfc(npsfc))) then
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),14)
             out_CloudPressure(ix,it) = nint(psfc0)
          endif
          if ((cpp.gt.0.).and.(cpp .lt. lut_pcld(1))) then
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),14)
             out_CloudPressure(ix,it) = nint(lut_pcld(1))
          endif
       endif

      !hqw skip out_CloudPressureSTD as inp_CloudPressureSTD is not read
      !out_CloudPressureSTD(ix,it)=inp_CloudPressureSTD(ix,it)
      !if(out_CloudPressureSTD(ix,it).le.iFillValue) &
      !     out_CloudPressureSTD(ix,it)=int(iFillValue, kind=2)

      !hqw moved out_TerrainPressure to m_cal_pscene and repurposed it
      !out_TerrainPressure(ix,it) = psfc0

      !=====
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
