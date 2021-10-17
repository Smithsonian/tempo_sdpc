module m_cal_pscene
  public cal_pscene
contains
!*********************
subroutine cal_pscene
  !*********************
  use m_vars
  use m_read_GMI
  use m_read_DEM
  use m_read_hdf5
  use m_scd_adjust

  implicit none

  !real,dimension(nsza)::lut_rsza
  integer::ialb, isza, ivza, iraa, ipsfc, ipcld!, irsfc
  integer::ialb1,isza1,ivza1,iraa1,ipsfc1!,irsfc1
  integer::ialb2,isza2,ivza2,iraa2,ipsfc2!,irsfc2
  real::   walb1,wsza1,wvza1,wraa1,wpsfc1,vpsfc1!,apsfc1,vcd1,tvcd1,wrsfc1
  real::   walb2,wsza2,wvza2,wraa2,wpsfc2,vpsfc2!,apsfc2,vcd2,tvcd2,wrsfc2
  real::yy1,yy2,ww1,ww2,rr1,rr2,wr1,wr2!,xxx,yyy,xx1,xx2,
  real(kind=8)::cpp, aaa
  integer(kind=4)::iflag
  integer(kind=4)::ierr!,status
  integer(kind=4)::nt,nx!nw,nwc
  integer(kind=4)::it,ix!,iw,iwc,iw1,iw2,dww

  integer(kind=4)::gmi_ix1,gmi_ix2,gmi_iy1,gmi_iy2, iternum
  real::gmi_wx1,gmi_wx2,gmi_wy1,gmi_wy2
  real::pp11,pp12,pp21,pp22,pp1,pp2
  real::tt11,tt12,tt21,tt22,tt1,tt2
  real (kind=4), dimension(:), allocatable:: tt, pp
  !real(kind=4)::sum1_vcd,avg_tvcd
  integer(kind=4)::ip!,gmi_ix,gmi_iy

  !real::a11111,a11112,a11121,a11122,a11211,a11212,a11221,a11222,a12111,a12112,a12121,a12122,a12211,a12212,a12221,a12222
  !real::a21111,a21112,a21121,a21122,a21211,a21212,a21221,a21222,a22111,a22112,a22121,a22122,a22211,a22212,a22221,a22222
  real::a1111,a1112,a1121,a1122,a1211,a1212,a1221,a1222,a2111,a2112,a2121,a2122,a2211,a2212,a2221,a2222
  real::a111,a112,a121,a122,a211,a212,a221,a222
  real::a11,a12,a21,a22
  real::a1,a2
  real::rad0,rad1,rad2,rrr0,rrr1,rrr2
  real::sbar,tran!,ler
  real::TerrainLER440,TerrainLER466
  real::SceneLER440,SceneLER466
  real::SceneCPP
  real::scdm!,omi_amf,vsfc0,cal_vcd
  real::scdmorg, scdadj, temp_cpp, t8p, temp_t8p, delta_temp !hqw addition
  real(kind=8),dimension(nalb)::temp_ler_alb466,temp_ler_alb440
  real(kind=8),dimension(npsfc)::lev_ler_alb466,lev_ler_alb440
  real(kind=8),dimension(npsfc)::lev_ler_amf
  real::ler440,ler466

  real,dimension(npcld)::amfvcd

  real::diff,diff_save,pdiff
  real::x0,x1,x2,xx
  real::y0,y1,y2,yy
  integer(kind=4)::ipp!,ipsfc_save
  !integer::iyes

  real::pi,dtor

  real::vpsfc0
  integer::ipsfc0
  !real::cal_ecf,cal_crf

  ! ------
  ! refine
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.

  nt=rad_NumTimes
  nx=rad_nXtrack
  !nw=rad_nWavel
  !nwc=rad_nWavelCoef

  ! allocate dimensions for outputs
  allocate(out_SurfaceLER440(nx,nt),stat=ierr)
  allocate(out_SurfaceLER466(nx,nt),stat=ierr)
  allocate(out_SceneLER440(nx,nt),stat=ierr)
  allocate(out_SceneLER466(nx,nt),stat=ierr)
  allocate(out_ScenePressure(nx,nt),stat=ierr)
  allocate(out_SlantColumnSceneO2O2(nx,nt),stat=ierr)
  allocate(out_SlantColumnTerrainO2O2(nx,nt),stat=ierr)
  allocate(out_O2O2SceneTemperature(nx,nt),stat=ierr)
  allocate(out_O2O2TerrainTemperature(nx,nt),stat=ierr)

  out_SurfaceLER440=fFillValue
  out_SurfaceLER466=fFillValue
  out_SceneLER440=fFillValue
  out_SceneLER466=fFillValue
  out_ScenePressure=fFillValue
  out_SlantColumnSceneO2O2=fFillValue
  out_SlantColumnTerrainO2O2=fFillValue
  out_O2O2SceneTemperature=fFillValue
  out_O2O2TerrainTemperature=fFillValue

  ! allocate and initialize local arrays
  allocate(tt(nlayers), pp(nlayers+1))
  tt(1:nlayers) = -999.0
  pp(1:nlayers+1) = -999.0

  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========

      pflag00=0
      pflag01=0
      if(rad_Latitude(ix,it) .lt. -90.) pflag00=pflag00+1
      !hqw changed 86. to max_SZA, 72. to max_VZA
      sza0 = rad_SolarZenithAngle(ix,it)
      vza0 = rad_ViewingZenithAngle(ix,it)
      if((sza0 .lt. 0.) .or. (sza0 .gt. max_SZA)) pflag01=pflag01+1
      if((vza0 .lt. 0.) .or. (vza0 .gt. max_VZA)) pflag01=pflag01+1
      if((pflag00 .ge. 1) .or. (pflag01 .ge. 1)) go to 990 ! skip all, start next pixel

      !the following has been commented out, as they are not used
      !cal_ecf=out_EffectiveCloudFraction(ix,it)
      !cal_crf=out_CloudRadianceFraction466(ix,it)
      !if((cal_ecf .gt. -990.0).and.(cal_ecf .lt. 0.0)) cal_ecf=0.0
      !if((cal_ecf .gt. 1.0).and.(cal_ecf .lt. 990.0)) cal_ecf=1.0
      !if((cal_crf .gt. -990.0).and.(cal_crf .lt. 0.0)) cal_crf=0.0
      !if((cal_crf .gt. 1.0).and.(cal_crf .lt. 990.0)) cal_crf=1.0

      raa0=out_RelativeAzimuthAngle(ix,it)
      !hqw out_RelativeAzimuthAbgle is calculated in cal_ecf
      ! it should have already be within [0.,180.] range
      ! it does not hurt to check again here, though seems unnecessary
      ! pixels with invalid angles should have been skipped vis pflag01
      if(raa0 .gt. 180.) raa0=360.-raa0

      !hqw inp_TerrainPressure used to come from OMCLDO2
      !psfc0 will be replaced by climatology
      psfc0 = -999. ! set to a temporary value here

      ! -----------------------------
      ! option for SlantColumnDensity
      ! -----------------------------
      !hqw out_SlantColumnAmount(ix,it) was now assigned in cal_ocp.f90
      !   thus comment out this section
      !if(name_option_SlantColumnDensity.eq.'NASA') then
      !  out_SlantColumnAmountO2O2(ix,it)=nasa_SlantColumnAmountO2O2(ix,it)
      !endif
      !hqw NASA product and the LUTs assumes O2O2 in unit of 10^43 molec^2/cm^5
      !  SCD < -9. are set to fFillValue -1.2676506E30 defined in m_vars.f90
      !this is now taken care of in m_read_input_tio.f90, thus comment it out
      !if(out_SlantColumnAmountO2O2(ix,it).le.-9.) out_SlantColumnAmountO2O2(ix,it)=fFillValue

      !hqw add scdmorg, skip calculation if <0., start next pixel
      scdmorg = nasa_SlantColumnAmountO2O2(ix,it)
      if (scdmorg .lt. 0.) go to 990

      ! ----------------------------------------------
      ! option for TemperaturePressure/SurfacePressure
      ! ----------------------------------------------
      if((name_option_TemperaturePressure.eq.'GMI').or. &
           (name_option_TemperaturePressure.eq.'DEM').or. &
           (name_option_TemperaturePressure.eq.'BDEM')) then
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
        if (psfc0.gt.0.0) then
          call read_GMI_VCD(pp,tt)
          vvcd=gmi_vcd
        else
          vvcd(1:npcld) = -999.0
        endif
      endif

      !hqw DEM is not currently supported, DO NOT USE
      !comment out for now
      !if(name_option_TemperaturePressure.eq.'DEM') then
      !  !hqw TEMPO does not read inp_TerrainPressure
      !  !psfc0=real(inp_TerrainPressure(ix,it))
      !  psfc0 = l2_TerrainPressure(ix,it)
      !  do ip=1,gmi_np
      !    tt11=gmi_temperature(gmi_ix1,gmi_iy1,ip)
      !    tt12=gmi_temperature(gmi_ix1,gmi_iy2,ip)
      !    tt21=gmi_temperature(gmi_ix2,gmi_iy1,ip)
      !    tt22=gmi_temperature(gmi_ix2,gmi_iy2,ip)
      !    tt1=(gmi_wy2*tt11+gmi_wy1*tt12)/(gmi_wy1+gmi_wy2)
      !    tt2=(gmi_wy2*tt21+gmi_wy1*tt22)/(gmi_wy1+gmi_wy2)
      !    tt(ip)=(gmi_wx2*tt1+gmi_wx1*tt2)/(gmi_wx1+gmi_wx2)
      !  end do
      !  call read_DEM_VCD(psfc0,tt,pp)
      !  vvcd=dem_vcd
      !endif

      !hqw BDEM is not currently supported, DO NOT USE
      ! comment out for now
      !if(name_option_TemperaturePressure.eq.'BDEM') then
      !  psfc0=BDEM_TerrainPressure(ix,it)
      !  do ip=1,gmi_np
      !    tt11=gmi_temperature(gmi_ix1,gmi_iy1,ip)
      !    tt12=gmi_temperature(gmi_ix1,gmi_iy2,ip)
      !    tt21=gmi_temperature(gmi_ix2,gmi_iy1,ip)
      !    tt22=gmi_temperature(gmi_ix2,gmi_iy2,ip)
      !    tt1=(gmi_wy2*tt11+gmi_wy1*tt12)/(gmi_wy1+gmi_wy2)
      !    tt2=(gmi_wy2*tt21+gmi_wy1*tt22)/(gmi_wy1+gmi_wy2)
      !    tt(ip)=(gmi_wx2*tt1+gmi_wx1*tt2)/(gmi_wx1+gmi_wx2)
      !  end do
      !  call read_DEM_VCD(psfc0,tt,pp)
      !  vvcd=dem_vcd
      !endif

      if(name_option_TemperaturePressure.eq.'GEOS5') then
        !geos_Pressure & geos_temperature are assigned in read_geoscf
        !they are TOA->BOA, the opposite of the GEOS-CF original order
        psfc0=geos_Pressure(ix,it,geos_np+1)
        do ip=1,geos_np
          pp(ip)=geos_Pressure(ix,it,ip)
          tt(ip)=geos_temperature(ix,it,ip)
        end do
        pp(geos_np+1)=psfc0
        if (psfc0 .gt. 0.0) then
          call read_GEOS5_VCD(pp,tt)
          vvcd=geos_vcd
        else
          vvcd(1:npcld) = -999.0
        endif
      endif

      !hqw skip if vvcd<0. due to improper psfc0
      if (vvcd(1) .lt. 0.) go to  990

      ! -----------------
      ! set nodes for LUT
      ! -----------------

      isza1=-9; isza2=-9
      do isza=1,nsza-1
        if((sza0 .ge. lut_sza(isza)) .and. (sza0 .le. lut_sza(isza+1))) then
          isza1=isza
          isza2=isza+1
          wsza1=sza0-lut_sza(isza)
          wsza2=lut_sza(isza+1)-sza0
        endif
      end do
      if(isza1 .lt. 0) go to 990

      ivza1=-9; ivza2=-9
      do ivza=1,nvza-1
        if((vza0 .ge. lut_vza(ivza)) .and. (vza0 .le. lut_vza(ivza+1))) then
          ivza1=ivza
          ivza2=ivza+1
          wvza1=vza0-lut_vza(ivza)
          wvza2=lut_vza(ivza+1)-vza0
        endif
      end do
      if(ivza1 .lt. 0) go to 990

      iraa1=-9; iraa2=-9
      do iraa=1,nraa-1
        if((raa0 .ge. lut_raa(iraa)) .and. (raa0 .le. lut_raa(iraa+1))) then
          iraa1=iraa
          iraa2=iraa+1
          wraa1=raa0-lut_raa(iraa)
          wraa2=lut_raa(iraa+1)-raa0
        endif
      end do
      if(iraa1 .lt. 0) go to 990

      !the following line ensures psfc0 to be within LUT range
      !note psfc0 < 0. should have already been skipped before
      if(psfc0 .gt. lut_psfc(npsfc)) psfc0=lut_psfc(npsfc)
      ipsfc1=-9; ipsfc2=-9
      do ipsfc=1,npsfc-1
        if((psfc0 .gt. lut_psfc(ipsfc)) .and. (psfc0 .le. lut_psfc(ipsfc+1))) then
          ipsfc1=ipsfc
          ipsfc2=ipsfc+1
          wpsfc1=psfc0-lut_psfc(ipsfc)
          wpsfc2=lut_psfc(ipsfc+1)-psfc0
        endif
      end do
      if(ipsfc1 .lt. 0) then
        !    write(*,*) " *** Pscene: Check Surface Pressure *** ",psfc0
        go to 990
      endif

      !------------------------------------------------------------
      ! name_option_SceneAlbedoAtTerrain
      !   yes (Ascene at Psfc); no (Ascene at each P level)
      !------------------------------------------------------------

      !+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
      if((name_option_SceneAlbedoAtTerrain.eq.'yes') .or. &
           (name_option_SceneAlbedoAtTerrain.eq.'both')) then
      !+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0

         !hqw add local variable initialization for pixel
         scdm = fFillValue
         scdadj = fFillValue
         t8p = 273.
         temp_t8p = 273.
         delta_temp = fFillValue
         temp_cpp = fFillValue
         cpp = fFillValue

        !------------------------------------------
        ! calculate LER at terrain pressure (Psfc)
        !------------------------------------------
        !------------------------
        !1. SurfaceLER at 466 nm
        !------------------------

        do ialb=1,nalb
          a1111=lut_rad_clr(ialb,isza1,ivza1,iraa1,ipsfc1)
          a1112=lut_rad_clr(ialb,isza1,ivza1,iraa1,ipsfc2)
          a1121=lut_rad_clr(ialb,isza1,ivza1,iraa2,ipsfc1)
          a1122=lut_rad_clr(ialb,isza1,ivza1,iraa2,ipsfc2)
          a1211=lut_rad_clr(ialb,isza1,ivza2,iraa1,ipsfc1)
          a1212=lut_rad_clr(ialb,isza1,ivza2,iraa1,ipsfc2)
          a1221=lut_rad_clr(ialb,isza1,ivza2,iraa2,ipsfc1)
          a1222=lut_rad_clr(ialb,isza1,ivza2,iraa2,ipsfc2)
          a2111=lut_rad_clr(ialb,isza2,ivza1,iraa1,ipsfc1)
          a2112=lut_rad_clr(ialb,isza2,ivza1,iraa1,ipsfc2)
          a2121=lut_rad_clr(ialb,isza2,ivza1,iraa2,ipsfc1)
          a2122=lut_rad_clr(ialb,isza2,ivza1,iraa2,ipsfc2)
          a2211=lut_rad_clr(ialb,isza2,ivza2,iraa1,ipsfc1)
          a2212=lut_rad_clr(ialb,isza2,ivza2,iraa1,ipsfc2)
          a2221=lut_rad_clr(ialb,isza2,ivza2,iraa2,ipsfc1)
          a2222=lut_rad_clr(ialb,isza2,ivza2,iraa2,ipsfc2)

          a111=(wpsfc2*a1111+wpsfc1*a1112)/(wpsfc1+wpsfc2)
          a112=(wpsfc2*a1121+wpsfc1*a1122)/(wpsfc1+wpsfc2)
          a121=(wpsfc2*a1211+wpsfc1*a1212)/(wpsfc1+wpsfc2)
          a122=(wpsfc2*a1221+wpsfc1*a1222)/(wpsfc1+wpsfc2)
          a211=(wpsfc2*a2111+wpsfc1*a2112)/(wpsfc1+wpsfc2)
          a212=(wpsfc2*a2121+wpsfc1*a2122)/(wpsfc1+wpsfc2)
          a221=(wpsfc2*a2211+wpsfc1*a2212)/(wpsfc1+wpsfc2)
          a222=(wpsfc2*a2221+wpsfc1*a2222)/(wpsfc1+wpsfc2)

          a11=(wraa2*a111+wraa1*a112)/(wraa1+wraa2)
          a12=(wraa2*a121+wraa1*a122)/(wraa1+wraa2)
          a21=(wraa2*a211+wraa1*a212)/(wraa1+wraa2)
          a22=(wraa2*a221+wraa1*a222)/(wraa1+wraa2)

          a1=(wvza2*a11+wvza1*a12)/(wvza1+wvza2)
          a2=(wvza2*a21+wvza1*a22)/(wvza1+wvza2)

          cal_ler_r466(ialb)=(wsza2*a1+wsza1*a2)/(wsza1+wsza2)
        end do

        ! calculate transmittance and sbar at R=0.0(1),0.1(7), and 0.2(12)
        rad0=real(cal_ler_r466(1), kind=4)
        rad1=real(cal_ler_r466(7), kind=4)
        rad2=real(cal_ler_r466(12), kind=4)
        rrr0=lut_alb(1)
        rrr1=lut_alb(7)
        rrr2=lut_alb(12)
        tran=(1./rrr1-1./rrr2)/(1./(rad1-rad0)-1./(rad2-rad0))
        sbar=1./rrr1-tran/(rad1-rad0)
        ler466=(rad_of_irr466(ix,it)-rad0)/(tran+sbar*(rad_of_irr466(ix,it)-rad0))
        if(ler466.lt.0.0) ler466=0.0
        if(ler466.gt.1.0) ler466=1.0

        alb0=ler466
        !  alb0=inp_SceneAlbedo(ix,it) ! this is from KNMI

        !------------------------
        !2. SurfaceLER at 440 nm
        !------------------------

        do ialb=1,nalb
          a1111=lut_rad_clr440(ialb,isza1,ivza1,iraa1,ipsfc1)
          a1112=lut_rad_clr440(ialb,isza1,ivza1,iraa1,ipsfc2)
          a1121=lut_rad_clr440(ialb,isza1,ivza1,iraa2,ipsfc1)
          a1122=lut_rad_clr440(ialb,isza1,ivza1,iraa2,ipsfc2)
          a1211=lut_rad_clr440(ialb,isza1,ivza2,iraa1,ipsfc1)
          a1212=lut_rad_clr440(ialb,isza1,ivza2,iraa1,ipsfc2)
          a1221=lut_rad_clr440(ialb,isza1,ivza2,iraa2,ipsfc1)
          a1222=lut_rad_clr440(ialb,isza1,ivza2,iraa2,ipsfc2)
          a2111=lut_rad_clr440(ialb,isza2,ivza1,iraa1,ipsfc1)
          a2112=lut_rad_clr440(ialb,isza2,ivza1,iraa1,ipsfc2)
          a2121=lut_rad_clr440(ialb,isza2,ivza1,iraa2,ipsfc1)
          a2122=lut_rad_clr440(ialb,isza2,ivza1,iraa2,ipsfc2)
          a2211=lut_rad_clr440(ialb,isza2,ivza2,iraa1,ipsfc1)
          a2212=lut_rad_clr440(ialb,isza2,ivza2,iraa1,ipsfc2)
          a2221=lut_rad_clr440(ialb,isza2,ivza2,iraa2,ipsfc1)
          a2222=lut_rad_clr440(ialb,isza2,ivza2,iraa2,ipsfc2)

          a111=(wpsfc2*a1111+wpsfc1*a1112)/(wpsfc1+wpsfc2)
          a112=(wpsfc2*a1121+wpsfc1*a1122)/(wpsfc1+wpsfc2)
          a121=(wpsfc2*a1211+wpsfc1*a1212)/(wpsfc1+wpsfc2)
          a122=(wpsfc2*a1221+wpsfc1*a1222)/(wpsfc1+wpsfc2)
          a211=(wpsfc2*a2111+wpsfc1*a2112)/(wpsfc1+wpsfc2)
          a212=(wpsfc2*a2121+wpsfc1*a2122)/(wpsfc1+wpsfc2)
          a221=(wpsfc2*a2211+wpsfc1*a2212)/(wpsfc1+wpsfc2)
          a222=(wpsfc2*a2221+wpsfc1*a2222)/(wpsfc1+wpsfc2)

          a11=(wraa2*a111+wraa1*a112)/(wraa1+wraa2)
          a12=(wraa2*a121+wraa1*a122)/(wraa1+wraa2)
          a21=(wraa2*a211+wraa1*a212)/(wraa1+wraa2)
          a22=(wraa2*a221+wraa1*a222)/(wraa1+wraa2)

          a1=(wvza2*a11+wvza1*a12)/(wvza1+wvza2)
          a2=(wvza2*a21+wvza1*a22)/(wvza1+wvza2)

          cal_ler_r440(ialb)=(wsza2*a1+wsza1*a2)/(wsza1+wsza2)
        end do

        ! calculate transmittance and sbar at R=0.0(1),0.1(7), and 0.2(12)
        rad0=real(cal_ler_r440(1), kind=4)
        rad1=real(cal_ler_r440(7), kind=4)
        rad2=real(cal_ler_r440(12), kind=4)
        rrr0=lut_alb(1)
        rrr1=lut_alb(7)
        rrr2=lut_alb(12)
        tran=(1./rrr1-1./rrr2)/(1./(rad1-rad0)-1./(rad2-rad0))
        sbar=1./rrr1-tran/(rad1-rad0)
        ler440=(rad_of_irr440(ix,it)-rad0)/(tran+sbar*(rad_of_irr440(ix,it)-rad0))
        if(ler440.lt.0.0) ler440=0.0
        if(ler440.gt.1.0) ler440=1.0

        !-----
        ! find alb node for lut_amf_ler
        ! alb0 was assigned ler466 above

        ialb1=-9; ialb2=-9
        do ialb=1,nalb-1
          if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
            ialb1=ialb
            ialb2=ialb+1
            walb1=alb0-lut_alb(ialb)
            walb2=lut_alb(ialb+1)-alb0
          endif
        end do
        if(ialb1 .lt. 0) go to 990

        ! -----------------------
        ! AMF at each cloud level
        ! -----------------------
        do ipcld=1,npcld
          ipsfc=ipcld
          a1111=lut_amf_ler(ialb1,isza1,ivza1,iraa1,ipsfc,ipcld)
          a1112=lut_amf_ler(ialb1,isza1,ivza1,iraa2,ipsfc,ipcld)
          a1121=lut_amf_ler(ialb1,isza1,ivza2,iraa1,ipsfc,ipcld)
          a1122=lut_amf_ler(ialb1,isza1,ivza2,iraa2,ipsfc,ipcld)
          a1211=lut_amf_ler(ialb1,isza2,ivza1,iraa1,ipsfc,ipcld)
          a1212=lut_amf_ler(ialb1,isza2,ivza1,iraa2,ipsfc,ipcld)
          a1221=lut_amf_ler(ialb1,isza2,ivza2,iraa1,ipsfc,ipcld)
          a1222=lut_amf_ler(ialb1,isza2,ivza2,iraa2,ipsfc,ipcld)
          a2111=lut_amf_ler(ialb2,isza1,ivza1,iraa1,ipsfc,ipcld)
          a2112=lut_amf_ler(ialb2,isza1,ivza1,iraa2,ipsfc,ipcld)
          a2121=lut_amf_ler(ialb2,isza1,ivza2,iraa1,ipsfc,ipcld)
          a2122=lut_amf_ler(ialb2,isza1,ivza2,iraa2,ipsfc,ipcld)
          a2211=lut_amf_ler(ialb2,isza2,ivza1,iraa1,ipsfc,ipcld)
          a2212=lut_amf_ler(ialb2,isza2,ivza1,iraa2,ipsfc,ipcld)
          a2221=lut_amf_ler(ialb2,isza2,ivza2,iraa1,ipsfc,ipcld)
          a2222=lut_amf_ler(ialb2,isza2,ivza2,iraa2,ipsfc,ipcld)

          ! -----------------
          ! in case wsfc2=0.0
          ! -----------------
          if((a1111.lt.0.0) .or. (a1112.lt.0.0)) then
            cal_ler_amf(ipcld)=-999.
            go to 898 ! goto next cloud level
          endif

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

          cal_ler_amf(ipcld)=(walb2*a1+walb1*a2)/(walb1+walb2)
898       continue
        end do !ipcld

        ! -----------
        ! check psfc0
        ! -----------
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
          endif
          if((psfc0 .ge. lut_psfc(ipsfc)) .and. (psfc0 .lt. lut_psfc(ipsfc+1))) then
            ipsfc0=ipsfc
          endif
        end do

        if(ipsfc1 .lt. 0) then
          write(*,*) " *** Pcld2: Check Surface Pressure *** ",psfc0
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),0)
        else
          vpsfc0=(wpsfc1*vpsfc2+wpsfc2*vpsfc1)/(wpsfc1+wpsfc2)
        endif

        ! -----------------
        ! calculate AMF*VCD
        ! -----------------
        !hqw: but cal_ler_amf may be -999., how is that handled?
        do ipcld=1,npcld
          aaa = real(cal_ler_amf(ipcld),kind=4)
          !hqw added check for aaa>0.
          !vvcd should always > 0., otherwise it should have been skipped
          if (aaa .gt. 0.) then
             amfvcd(ipcld)=real(aaa*vvcd(ipcld), kind=4)
          else
             ! skipped to next pixel
             amfvcd(ipcld) = -999.
             go to 990
          endif
        end do

        !hqw initialize local vairable before scd T-correction iteration
        scdm = scdmorg
        scdadj = scdmorg
        t8p = 273.
        temp_t8p = 273.

        iternum = 0
777     continue ! T iteration come back here

        iflag=-1

        if(scdm.le.amfvcd(1)) then
          iflag=0
        endif

        do ipcld=1,npcld-1
          if((scdm.gt.amfvcd(ipcld)).and.(scdm.le.amfvcd(ipcld+1))) then
            iflag=1
            yy1=lut_pcld(ipcld)
            yy2=lut_pcld(ipcld+1)
            ww1=scdm-amfvcd(ipcld)
            ww2=amfvcd(ipcld+1)-scdm
          endif
        end do

        if(iflag .ge. 1) then
          cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)
        else if(iflag .eq. 0) then
          x0=0.0
          x1=lut_pcld(1)
          x2=lut_pcld(2)
          y0=0.0
          y1=amfvcd(1)
          y2=amfvcd(2)

          xx=x1
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff_save=abs(scdm-yy)

          do ipp=1,100
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
970       continue
          cpp=xx
        else
          x0=lut_pcld(npcld-0)
          x1=lut_pcld(npcld-1)
          x2=lut_pcld(npcld-2)
          y0=amfvcd(npcld-0)
          y1=amfvcd(npcld-1)
          y2=amfvcd(npcld-2)

          xx=lut_pcld(npcld-0)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff_save=abs(scdm-yy)

          do ipp=1,1000
            xx=x0+real(ipp)
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 980
            else
              diff_save=diff
              if(ipp.ge.5000) then
                xx=9999.
              endif
            endif
          end do

980       continue
          cpp=xx
        endif

        !hqw scd temperature correction
        !Oct2021 temporarily use T at half cpp
        !       need to tune with synthetic profile
        temp_cpp = cpp * 0.5
        if (temp_cpp .gt. 50. .and. temp_cpp .lt. 1100.) then
          if (name_option_TemperaturePressure .eq. 'GMI') then
            call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else if (name_option_TemperaturePressure .eq. 'GEOS5') then
            call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else
            temp_t8p = t8p
          endif
        else
           temp_t8p = t8p
        endif
        iternum = iternum + 1

        !hqw test if terminate iteration
        delta_temp = abs(t8p - temp_t8p)
        if (delta_temp .ge. dt_threshold .and. iternum .lt. max_scd_iter) then
            t8p = temp_t8p
            scdm = scdadj
            go to 777 ! do another iteration
        endif

        !hqw assign output
        if (scdm .gt. 0.) then
           out_SlantColumnTerrainO2O2(ix,it) = scdm
           out_O2O2TerrainTemperature(ix,it) = t8p
        else
           out_SlantColumnTerrainO2O2(ix,it)=nasa_SlantColumnAmountO2O2(ix,it)
           out_O2O2TerrainTemperature(ix,it) = 273.
        endif

      !+0+0+0+0+0+0+0+0+0+0
      endif !name_option_SceneAlbedoAtTerrain=='yes'//'both'
      !+0+0+0+0+0+0+0+0+0+0

      TerrainLER466=ler466
      TerrainLER440=ler440

      !+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1
      if((name_option_SceneAlbedoAtTerrain.eq.'no') .or. &
           (name_option_SceneAlbedoAtTerrain.eq.'both')) then
        !+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1

         !hqw initialize local variables for name_option
         scdm = fFillValue
         scdadj = fFillValue
         temp_cpp = fFillValue
         delta_temp = fFillValue
         t8p = 273.
         temp_t8p = 273.
         cpp = fFillValue

        !--------------------------------------
        ! calculate LER at each pressure level
        !--------------------------------------
        ! real(kind=8),dimension(nalb)::temp_ler_alb
        ! real(kind=8),dimension(npsfc)::lev_ler_alb
        ! real(kind=8),dimension(npsfc)::lev_ler_amf

        !----------------------
        !1. SceneLER at 466 nm
        !----------------------
        do ipsfc=1,npsfc
          ipcld=ipsfc
          do ialb=1,nalb
            a111=lut_rad_clr(ialb,isza1,ivza1,iraa1,ipsfc)
            a112=lut_rad_clr(ialb,isza1,ivza1,iraa2,ipsfc)
            a121=lut_rad_clr(ialb,isza1,ivza2,iraa1,ipsfc)
            a122=lut_rad_clr(ialb,isza1,ivza2,iraa2,ipsfc)
            a211=lut_rad_clr(ialb,isza2,ivza1,iraa1,ipsfc)
            a212=lut_rad_clr(ialb,isza2,ivza1,iraa2,ipsfc)
            a221=lut_rad_clr(ialb,isza2,ivza2,iraa1,ipsfc)
            a222=lut_rad_clr(ialb,isza2,ivza2,iraa2,ipsfc)
            a11=(wraa2*a111+wraa1*a112)/(wraa1+wraa2)
            a12=(wraa2*a121+wraa1*a122)/(wraa1+wraa2)
            a21=(wraa2*a211+wraa1*a212)/(wraa1+wraa2)
            a22=(wraa2*a221+wraa1*a222)/(wraa1+wraa2)
            a1=(wvza2*a11+wvza1*a12)/(wvza1+wvza2)
            a2=(wvza2*a21+wvza1*a22)/(wvza1+wvza2)
            temp_ler_alb466(ialb)=(wsza2*a1+wsza1*a2)/(wsza1+wsza2)
          end do

          ! calculate transmittanc and sbar at R=0.0(1),0.1(7), and 0.2(12)
          rad0=real(temp_ler_alb466(1), kind=4)
          rad1=real(temp_ler_alb466(7), kind=4)
          rad2=real(temp_ler_alb466(12), kind=4)
          rrr0=lut_alb(1)
          rrr1=lut_alb(7)
          rrr2=lut_alb(12)
          tran=(1./rrr1-1./rrr2)/(1./(rad1-rad0)-1./(rad2-rad0))
          sbar=1./rrr1-tran/(rad1-rad0)
          ler466=(rad_of_irr466(ix,it)-rad0)/(tran+sbar*(rad_of_irr466(ix,it)-rad0))
          if(ler466.lt.0.0) ler466=0.0
          if(ler466.gt.1.0) ler466=1.0
          lev_ler_alb466(ipsfc)=ler466
        end do

        !----------------------
        !2. SceneLER at 440 nm
        !----------------------
        do ipsfc=1,npsfc
          ipcld=ipsfc
          do ialb=1,nalb
            a111=lut_rad_clr440(ialb,isza1,ivza1,iraa1,ipsfc)
            a112=lut_rad_clr440(ialb,isza1,ivza1,iraa2,ipsfc)
            a121=lut_rad_clr440(ialb,isza1,ivza2,iraa1,ipsfc)
            a122=lut_rad_clr440(ialb,isza1,ivza2,iraa2,ipsfc)
            a211=lut_rad_clr440(ialb,isza2,ivza1,iraa1,ipsfc)
            a212=lut_rad_clr440(ialb,isza2,ivza1,iraa2,ipsfc)
            a221=lut_rad_clr440(ialb,isza2,ivza2,iraa1,ipsfc)
            a222=lut_rad_clr440(ialb,isza2,ivza2,iraa2,ipsfc)
            a11=(wraa2*a111+wraa1*a112)/(wraa1+wraa2)
            a12=(wraa2*a121+wraa1*a122)/(wraa1+wraa2)
            a21=(wraa2*a211+wraa1*a212)/(wraa1+wraa2)
            a22=(wraa2*a221+wraa1*a222)/(wraa1+wraa2)
            a1=(wvza2*a11+wvza1*a12)/(wvza1+wvza2)
            a2=(wvza2*a21+wvza1*a22)/(wvza1+wvza2)
            temp_ler_alb440(ialb)=(wsza2*a1+wsza1*a2)/(wsza1+wsza2)
          end do

          ! calculate transmittance and sbar at R=0.0(1),0.1(7), and 0.2(12)
          rad0=real(temp_ler_alb440(1), kind=4)
          rad1=real(temp_ler_alb440(7), kind=4)
          rad2=real(temp_ler_alb440(12), kind=4)
          rrr0=lut_alb(1)
          rrr1=lut_alb(7)
          rrr2=lut_alb(12)
          tran=(1./rrr1-1./rrr2)/(1./(rad1-rad0)-1./(rad2-rad0))
          sbar=1./rrr1-tran/(rad1-rad0)
          ler440=(rad_of_irr440(ix,it)-rad0)/(tran+sbar*(rad_of_irr440(ix,it)-rad0))
          if(ler440.lt.0.0) ler440=0.0
          if(ler440.gt.1.0) ler440=1.0
          lev_ler_alb440(ipsfc)=ler440
        end do

        ! -----------------------
        ! AMF at each cloud level
        ! -----------------------
        do ipsfc=1,npsfc
          ipcld=ipsfc
          alb0=real(lev_ler_alb466(ipsfc), kind=4)
          ialb1=-9; ialb2=-9
          do ialb=1,nalb-1
            if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
              ialb1=ialb
              ialb2=ialb+1
              walb1=alb0-lut_alb(ialb)
              walb2=lut_alb(ialb+1)-alb0
            endif
          end do

          a1111=lut_amf_ler(ialb1,isza1,ivza1,iraa1,ipsfc,ipcld)
          a1112=lut_amf_ler(ialb1,isza1,ivza1,iraa2,ipsfc,ipcld)
          a1121=lut_amf_ler(ialb1,isza1,ivza2,iraa1,ipsfc,ipcld)
          a1122=lut_amf_ler(ialb1,isza1,ivza2,iraa2,ipsfc,ipcld)
          a1211=lut_amf_ler(ialb1,isza2,ivza1,iraa1,ipsfc,ipcld)
          a1212=lut_amf_ler(ialb1,isza2,ivza1,iraa2,ipsfc,ipcld)
          a1221=lut_amf_ler(ialb1,isza2,ivza2,iraa1,ipsfc,ipcld)
          a1222=lut_amf_ler(ialb1,isza2,ivza2,iraa2,ipsfc,ipcld)
          a2111=lut_amf_ler(ialb2,isza1,ivza1,iraa1,ipsfc,ipcld)
          a2112=lut_amf_ler(ialb2,isza1,ivza1,iraa2,ipsfc,ipcld)
          a2121=lut_amf_ler(ialb2,isza1,ivza2,iraa1,ipsfc,ipcld)
          a2122=lut_amf_ler(ialb2,isza1,ivza2,iraa2,ipsfc,ipcld)
          a2211=lut_amf_ler(ialb2,isza2,ivza1,iraa1,ipsfc,ipcld)
          a2212=lut_amf_ler(ialb2,isza2,ivza1,iraa2,ipsfc,ipcld)
          a2221=lut_amf_ler(ialb2,isza2,ivza2,iraa1,ipsfc,ipcld)
          a2222=lut_amf_ler(ialb2,isza2,ivza2,iraa2,ipsfc,ipcld)
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
          lev_ler_amf(ipsfc)=(walb2*a1+walb1*a2)/(walb1+walb2)
        end do

        ! -----------------
        ! calculate AMF*VCD
        ! -----------------
        !hqw re-assign local vairables
        scdm = scdmorg
        scdadj = scdm
        t8p = 273.
        temp_t8p = t8p

        !hqw comment out the following
        !if((scdm.gt.-10.0).and.(scdm.lt.0.0)) scdm=0.1
        !the following would not happen, as scdmorg<0. has been skipped
        if(scdm.lt.0.0) then
          !out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),6) in cal_ocp
          cpp=fFillValue
         ! ler=fFillValue
          go to 988
        endif

        iternum = 0
776     continue !hqw iteration come back here

        do ipcld=1,npcld
          ipsfc=ipcld
          amfvcd(ipcld)=real(lev_ler_amf(ipsfc), kind=4)*vvcd(ipcld)
        end do

        iflag=-1

        if(scdm.le.amfvcd(1)) then
          iflag=0
        endif

        do ipcld=1,npcld-1
          if((scdm.gt.amfvcd(ipcld)).and.(scdm.le.amfvcd(ipcld+1))) then
            iflag=1
            yy1=lut_pcld(ipcld)
            yy2=lut_pcld(ipcld+1)
            ww1=scdm-amfvcd(ipcld)
            ww2=amfvcd(ipcld+1)-scdm
          endif
        end do

        if(iflag .ge. 1) then
          cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)
        else if(iflag .eq. 0) then
          x0=0.0
          x1=lut_pcld(1)
          x2=lut_pcld(2)
          y0=0.0
          y1=amfvcd(1)
          y2=amfvcd(2)

          xx=x1
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff_save=abs(scdm-yy)

          do ipp=1,100
            xx=x1-real(ipp)
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 972
            else
              diff_save=diff
              if(ipp.le.1) then
                xx=-9999.
              endif
            endif
          end do
972       continue
          cpp=xx
        else
          x0=lut_pcld(npcld-0)
          x1=lut_pcld(npcld-1)
          x2=lut_pcld(npcld-2)
          y0=amfvcd(npcld-0)
          y1=amfvcd(npcld-1)
          y2=amfvcd(npcld-2)

          xx=lut_pcld(npcld-0)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff_save=abs(scdm-yy)

          do ipp=1,1000
            xx=x0+real(ipp)
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 982
            else
              diff_save=diff
              if(ipp.ge.5000) then
                xx=9999.
              endif
            endif
          end do

982       continue
          cpp=xx
        endif

        !hqw scd T-correction
        temp_cpp = cpp * 0.5
        if (temp_cpp .gt. 50. .and. temp_cpp .lt. 1100.) then
          if (name_option_TemperaturePressure .eq. 'GMI') then
            call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else if (name_option_TemperaturePressure .eq. 'GEOS5') then
            call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else
            temp_t8p = t8p
          endif
        else
           temp_t8p = t8p
        endif
        iternum = iternum + 1

        !hqw iteration termination
        delta_temp = abs(t8p - temp_t8p)
        if (delta_temp .gt. dt_threshold .and. iternum .lt. max_scd_iter) then
           scdm = scdadj
           t8p = temp_t8p
           go to 776
        endif

        !hqw assign output
        if (scdm .gt. 0.) then
            out_SlantColumnSceneO2O2(ix,it) = scdm
            out_O2O2SceneTemperature(ix,it) = t8p
        else
            out_SlantColumnSceneO2O2(ix,it)=nasa_SlantColumnAMountO2O2(ix,it)
            out_O2O2SceneTemperature(ix,it) = 273.
        endif

        !------------------------
        !calculate ler466 at cpp
        !------------------------
        iflag=-1
        do ipcld=1,npcld-1
          if((cpp.gt.lut_pcld(ipcld)).and.(cpp.le.lut_pcld(ipcld+1))) then
            iflag=1
            rr1=real(lev_ler_alb466(ipcld), kind=4)
            rr2=real(lev_ler_alb466(ipcld+1), kind=4)
            wr1=real(cpp, kind=4)-lut_pcld(ipcld)
            wr2=lut_pcld(ipcld+1)-real(cpp, kind=4)
          endif
        end do

        if(iflag .ge. 1) then
          ler466=(wr1*rr2+wr2*rr1)/(wr1+wr2)
        else
          x0=lut_pcld(npcld-0)
          x1=lut_pcld(npcld-1)
          x2=lut_pcld(npcld-2)
          y0=real(lev_ler_alb466(npcld-0), kind=4)
          y1=real(lev_ler_alb466(npcld-1), kind=4)
          y2=real(lev_ler_alb466(npcld-2), kind=4)
          xx=real(cpp, kind=4)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          ler466=yy
        endif
        if(ler466.lt.0.0) ler466=0.0
        if(ler466.gt.1.0) ler466=1.0

        !-------------------------------
        !calculate ler440 at cpp level
        !-------------------------------
        iflag=-1
        do ipcld=1,npcld-1
          if((cpp.gt.lut_pcld(ipcld)).and.(cpp.le.lut_pcld(ipcld+1))) then
            iflag=1
            rr1=real(lev_ler_alb440(ipcld), kind=4)
            rr2=real(lev_ler_alb440(ipcld+1), kind=4)
            wr1=real(cpp, kind=4)-lut_pcld(ipcld)
            wr2=lut_pcld(ipcld+1)-real(cpp, kind=4)
          endif
        end do

        if(iflag .ge. 1) then
          ler440=(wr1*rr2+wr2*rr1)/(wr1+wr2)
        else
          x0=lut_pcld(npcld-0)
          x1=lut_pcld(npcld-1)
          x2=lut_pcld(npcld-2)
          y0=real(lev_ler_alb440(npcld-0), kind=4)
          y1=real(lev_ler_alb440(npcld-1), kind=4)
          y2=real(lev_ler_alb440(npcld-2), kind=4)
          xx=real(cpp, kind=4)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          ler440=yy
        endif
        if(ler440.lt.0.0) ler440=0.0
        if(ler440.gt.1.0) ler440=1.0

988     continue

        !+1+1+1+1+1+1+1+1
      endif
      !+1+1+1+1+1+1+1+1

      SceneLER466=ler466
      SceneLER440=ler440
      SceneCPP=real(cpp, kind=4)

      if((TerrainLER466 .lt. -990.) .or. (TerrainLER466 .gt. 990.)) TerrainLER466=fFillValue
      if((TerrainLER466 .ge. -990.) .and. (TerrainLER466 .lt. 0.0)) TerrainLER466=0.0
      if((TerrainLER466 .gt.  1.0) .and. (TerrainLER466 .le. 990.)) TerrainLER466=1.0
      if((TerrainLER440 .lt. -990.) .or. (TerrainLER440 .gt. 990.)) TerrainLER440=fFillValue
      if((TerrainLER440 .ge. -990.) .and. (TerrainLER440 .lt. 0.0)) TerrainLER440=0.0
      if((TerrainLER440 .gt.  1.0) .and. (TerrainLER440 .le. 990.)) TerrainLER440=1.0
      if((SceneLER466 .lt. -990.) .or. (SceneLER466 .gt. 990.)) SceneLER466=fFillValue
      if((SceneLER466 .ge. -990.) .and. (SceneLER466 .lt. 0.0)) SceneLER466=0.0
      if((SceneLER466 .gt.  1.0) .and. (SceneLER466 .le. 990.)) SceneLER466=1.0
      if((SceneLER440 .lt. -990.) .or. (SceneLER440 .gt. 990.)) SceneLER440=fFillValue
      if((SceneLER440 .ge. -990.) .and. (SceneLER440 .lt. 0.0)) SceneLER440=0.0
      if((SceneLER440 .gt.  1.0) .and. (SceneLER440 .le. 990.)) SceneLER440=1.0
      if(SceneCPP .gt. 9999.) SceneCPP=fFillValue
      if(SceneCPP .lt. -9999.) SceneCPP=fFillValue
      if((SceneCPP.gt.psfc0).and.(SceneCPP.lt.9999.)) SceneCPP=psfc0

      out_ScenePressure(ix,it)=SceneCPP
      out_SceneLER466(ix,it)=SceneLER466
      out_SceneLER440(ix,it)=SceneLER440
      out_SurfaceLER466(ix,it)=TerrainLER466
      out_SurfaceLER440(ix,it)=TerrainLER440

      ! ---------------------------------------------------------
      ! Use scene pressure
      !   1. Over snow/ice
      !   2. ECF < min_ecf
      ! ---------------------------------------------------------
      if(name_option_SnowIce.eq.'Pscene') then
        if(btest(out_ProcessingQualityFlags(ix,it),4)) then
          pdiff=psfc0-SceneCPP
          if((SceneLER466.ge.0.2).and.(pdiff.lt.100.)) then
            out_CloudPressure(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
            out_CloudPressureNotClipped(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
            out_EffectiveCloudFraction(ix,it)=0.
            out_EffectiveCloudFractionNotClipped(ix,it)=0.
            out_CloudRadianceFraction466(ix,it)=0.
            out_CloudRadianceFractionNotClipped466(ix,it)=0.
            out_CloudRadianceFraction440(ix,it)=0.
            out_CloudRadianceFractionNotClipped440(ix,it)=0.
          endif
          if((SceneLER466.ge.0.2).and.(pdiff.ge.100.)) then
            out_CloudPressure(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
            out_CloudPressureNotClipped(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
            !hqw changed 1000 to 1. after out_EffectiveCloudFraction changed
            ! from [0,1000] to [0.,1.] range
            out_EffectiveCloudFraction(ix,it)=1.
            out_EffectiveCloudFractionNotClipped(ix,it)=1.
            out_CloudRadianceFraction466(ix,it)=1.
            out_CloudRadianceFractionNotClipped466(ix,it)=1.
            out_CloudRadianceFraction440(ix,it)=1.
            out_CloudRadianceFractionNotClipped440(ix,it)=1.
          endif
        endif
      endif

      if(name_option_MinECF.eq.'yes') then
        if(btest(out_ProcessingQualityFlags(ix,it),2)) then
          out_CloudPressure(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
          out_CloudPressureNotClipped(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
        endif
        if((out_CloudPressure(ix,it).gt.iFillValue).and.(out_CloudPressure(ix,it).le.100)) &
             out_CloudPressure(ix,it)=100
        if((out_CloudPressure(ix,it).ge.nint(psfc0)).and.(out_CloudPressure(ix,it).lt.5000)) &
             out_CloudPressure(ix,it)=nint(psfc0, kind=2)
      endif

      !  if(btest(out_ProcessingQualityFlags(ix,it),4)) then
      !    if(SceneLER466.lt.0.2) then
      !      write(48,*) it,ix,out_EffectiveCloudFraction(ix,it),SceneLER466,out_CloudPressure(ix,it),out_ScenePressure(ix,it),psfc0
      !    endif
      !  endif

990   continue

      !=====
    end do
  end do
  !=====

  !************************
end subroutine cal_pscene
!************************
end module m_cal_pscene
