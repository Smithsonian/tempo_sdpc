module m_cal_pscene
  public cal_pscene

contains
!*********************
subroutine cal_pscene
  !*********************
  use m_vars
  use m_read_GMI
!  use m_read_DEM
  use m_read_hdf5
  use m_scd_adjust

  implicit none

  !real,dimension(nsza)::lut_rsza
  integer::ialb, isza, ivza, iraa, ipsfc, ipcld
  integer::ialb1,isza1,ivza1,iraa1,ipsfc1
  integer::ialb2,isza2,ivza2,iraa2,ipsfc2
  real::   walb1,wsza1,wvza1,wraa1,wpsfc1,vpsfc1
  real::   walb2,wsza2,wvza2,wraa2,wpsfc2,vpsfc2
  real::yy1,yy2,ww1,ww2,rr1,rr2,wr1,wr2
  real(kind=4)::cpp, aaa, temp_cpp, cpp1st
  integer(kind=4)::iflag
  integer(kind=4)::ierr
  integer(kind=4)::nt,nx
  integer(kind=4)::it,ix

  integer(kind=4)::gmi_ix1,gmi_ix2,gmi_iy1,gmi_iy2, iternum
  real::gmi_wx1,gmi_wx2,gmi_wy1,gmi_wy2
  real::pp11,pp12,pp21,pp22,pp1,pp2
  real::tt11,tt12,tt21,tt22,tt1,tt2
  real (kind=4), dimension(:), allocatable:: tt, pp
  !real(kind=4)::sum1_vcd,avg_tvcd
  integer(kind=4)::ip

  integer(kind=4):: pflag00, pflag01

  real::a1111,a1112,a1121,a1122,a1211,a1212,a1221,a1222,a2111,a2112,a2121,a2122,a2211,a2212,a2221,a2222
  real::a111,a112,a121,a122,a211,a212,a221,a222
  real::a11,a12,a21,a22
  real::a1,a2
  real::rad0,rad1,rad2,rrr0,rrr1,rrr2
  real::sbar,tran
  real::TerrainLER440,TerrainLER466
  real::SceneLER440,SceneLER466
  real::SceneCPP
  real::scdm
  real::scdmorg, scdadj, t8p, temp_t8p, delta_temp !hqw addition
  real(kind=8),dimension(nalb)::temp_ler_alb466,temp_ler_alb440
  real(kind=8),dimension(npsfc)::lev_ler_alb466,lev_ler_alb440
  real(kind=8),dimension(npsfc)::lev_ler_amf
  real::ler440,ler466
  real:: fFillValue9, maxpress

  real,dimension(npcld)::amfvcd

  real::diff,diff_save,pdiff
  real::x0,x1,x2,xx
  real::y0,y1,y2,yy
  integer(kind=4)::ipp

  real::pi,dtor

  real::vpsfc0
  integer::ipsfc0

  ! ------
  ! refine
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.

  nt=rad_NumTimes
  nx=rad_nXtrack

  fFillValue9 = -9999.
  maxpress = 1500. !Pa for pscene

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

  !initialize array
  out_SurfaceLER440=fFillValue9
  out_SurfaceLER466=fFillValue9
  out_SceneLER440=fFillValue9
  out_SceneLER466=fFillValue9
  out_ScenePressure=fFillValue9
  out_SlantColumnSceneO2O2=fFillValue9
  out_SlantColumnTerrainO2O2=fFillValue9
  out_O2O2SceneTemperature=fFillValue9
  out_O2O2TerrainTemperature=fFillValue9

  ! allocate and initialize local arrays
  allocate(tt(nlayers), pp(nlayers+1))
  tt(1:nlayers) = fFillValue9
  pp(1:nlayers+1) = fFillValue9

  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========
! m_cal_ocp is called before this module, out_ProcessingQualityFlags
! bit 0 and 1 can be used to skip pixels with invalid geo or angles
      if (btest(out_ProcessingQualityFlags(ix,it),0) .or. &
          btest(out_ProcessingQualityFlags(ix,it),1)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15)
          go to 990
      endif
      
      sza0 = rad_SolarZenithAngle(ix,it)
      vza0 = rad_ViewingZenithAngle(ix,it)
      raa0 = out_RelativeAzimuthAngle(ix,it)

      !psfc0 will be replaced by climatology
      psfc0 = fFillValue9

      ! -----------------------------
      ! option for SlantColumnDensity
      ! -----------------------------
      !hqw add scdmorg, skip calculation if <0., start next pixel
      scdmorg = nasa_SlantColumnAmountO2O2(ix,it)
      if (scdmorg .lt. 0.) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15) 
         go to 990
      endif

      ! ----------------------------------------------
      ! option for TemperaturePressure/SurfacePressure
      ! ----------------------------------------------
      if((name_option_TemperaturePressure.eq.'GMI')) then !.or. &
!           (name_option_TemperaturePressure.eq.'DEM').or. &
!           (name_option_TemperaturePressure.eq.'BDEM')) then
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
          vvcd(1:npcld) = -9999.
        endif
      endif

      !hqw DEM is not currently supported, DO NOT USE
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
          vvcd(1:npcld) = -9999.
        endif
      endif

      !hqw skip if vvcd<0. due to improper psfc0
      if ((vvcd(1) .lt. 0.) .and. (vvcd(npcld) .lt. 0.)) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15)
          go to 990
      endif

      ! -----------------
      ! set nodes for LUT
      ! -----------------

      isza1=-9; isza2=-9; wsza1=0. ; wsza2=0.
      do isza=1,nsza-1
        if((sza0 .ge. lut_sza(isza)) .and. (sza0 .le. lut_sza(isza+1))) then
          isza1=isza
          isza2=isza+1
          wsza1=sza0-lut_sza(isza)
          wsza2=lut_sza(isza+1)-sza0
        endif
      end do
      if(isza1 .lt. 0) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15)
         go to 990
      endif

      ivza1=-9; ivza2=-9; wvza1=0.; wvza2=0.
      do ivza=1,nvza-1
        if((vza0 .ge. lut_vza(ivza)) .and. (vza0 .le. lut_vza(ivza+1))) then
          ivza1=ivza
          ivza2=ivza+1
          wvza1=vza0-lut_vza(ivza)
          wvza2=lut_vza(ivza+1)-vza0
        endif
      end do
      if(ivza1 .lt. 0) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15)
         go to 990
      endif

      iraa1=-9; iraa2=-9; wraa1=0.; wraa2=0.
      do iraa=1,nraa-1
        if((raa0 .ge. lut_raa(iraa)) .and. (raa0 .le. lut_raa(iraa+1))) then
          iraa1=iraa
          iraa2=iraa+1
          wraa1=raa0-lut_raa(iraa)
          wraa2=lut_raa(iraa+1)-raa0
        endif
      end do
      if(iraa1 .lt. 0) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15)
         go to 990
      endif

      !note psfc0 < 0. should have already been skipped 
      !the following ensures psfc0 to be within LUT range
      if(psfc0 .gt. lut_psfc(npsfc)) psfc0=lut_psfc(npsfc)

      ipsfc1=-9; ipsfc2=-9; wpsfc1=0.; wpsfc2=0.
      do ipsfc=1,npsfc-1
        if((psfc0 .gt. lut_psfc(ipsfc)) .and. (psfc0 .le. lut_psfc(ipsfc+1))) then
          ipsfc1=ipsfc
          ipsfc2=ipsfc+1
          wpsfc1=psfc0-lut_psfc(ipsfc)
          wpsfc2=lut_psfc(ipsfc+1)-psfc0
        endif
      end do
      if(ipsfc1 .lt. 0) then
         write(*,*) " Pscene: Check Surface Pressure *** ",psfc0,ix,it
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),15)
        go to 990
      endif

      !------------------------------------------------------------
      ! name_option_SceneAlbedoAtTerrain
      !   yes (Ascene at Psfc); no (Ascene at each P level)
      !------------------------------------------------------------
      if (name_option_SceneAlbedoAtTerrain .ne. 'both') then
        if (name_option_SceneAlbedoAtTerrain .eq. 'yes') then
          ! 'no' will be skipped, set bit 11
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11) 
        else ! .eq. 'no'
          ! 'yes' will be skipped, set bit 10
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
        endif
      endif
      !------------------------------------------------------------

      !+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
      if((name_option_SceneAlbedoAtTerrain.eq.'yes') .or. &
         (name_option_SceneAlbedoAtTerrain.eq.'both')) then
      !+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0

         !hqw local variable initialization for pixel
         scdm = fFillValue9
         scdadj = fFillValue9
         t8p = 273.
         temp_t8p = 273.
         delta_temp = fFillValue9
         temp_cpp = fFillValue9
         cpp = fFillValue9

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
        !hqw added logic
        if (rad_of_irr466(ix,it) .gt. 0.) then
           ler466=(rad_of_irr466(ix,it)-rad0)/(tran+sbar*(rad_of_irr466(ix,it)-rad0))
           if(ler466.lt.0.0) ler466=0.0
           if(ler466.gt.1.0) ler466=1.0
           alb0=ler466
        else
           ler466 = fFillValue9
           alb0 = ler466
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
           go to 333
        endif

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

        ! calculate tran and sbar at R=0.0(1),0.1(7), and 0.2(12)
        rad0=real(cal_ler_r440(1), kind=4)
        rad1=real(cal_ler_r440(7), kind=4)
        rad2=real(cal_ler_r440(12), kind=4)
        rrr0=lut_alb(1)
        rrr1=lut_alb(7)
        rrr2=lut_alb(12)
        tran=(1./rrr1-1./rrr2)/(1./(rad1-rad0)-1./(rad2-rad0))
        sbar=1./rrr1-tran/(rad1-rad0)
        !hqw added logic
        if (rad_of_irr440(ix,it) .gt. 0.) then
           ler440=(rad_of_irr440(ix,it)-rad0)/(tran+sbar*(rad_of_irr440(ix,it)-rad0))
           if(ler440.lt.0.0) ler440=0.0
           if(ler440.gt.1.0) ler440=1.0
        else
           ler440 = fFillValue9
        endif 

        !-----
        ! find alb node for lut_amf_ler using the calculated ler466
        ! alb0 was assigned ler466 before
        !-----

        ialb1=-9; ialb2=-9; walb1=0.; walb2=0.
        do ialb=1,nalb-1
          if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
            ialb1=ialb
            ialb2=ialb+1
            walb1=alb0-lut_alb(ialb)
            walb2=lut_alb(ialb+1)-alb0
          endif
        end do
        !if node not found, calculation will be skipped
        if(ialb1 .lt. 0) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
          go to 333
        endif

        ! -----------------------
        ! AMF at each cloud level
        ! ----------------------- 
        ! use iab1 & iab2 found for the calculated ler466 
        ! to interpolate AMF at each cloud level: cal_ler_amf(ipcld)
        !hqw lut_amf_ler(alb,sza,vza,raa,pcld,psfc),thus,
        !it is better to switch ipsfc and ipcld below
        !however, it does not matter for running, as ipsfc=ipcld
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
          ! in case wpsfc2=0.0
          ! -----------------
          !in AMF LUT, when Pcld>Psfc, entries are set to -999.
          !when clouds are below surface,the condition will be true
          !however, ipsfc=ipcld ensures this will not happen
          !as it assumes clouds are at the 'surface' for each cloud level
          if((a1111.lt.0.0) .or. (a1112.lt.0.0)) then
            cal_ler_amf(ipcld)=-9999. 
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
        ! calculate VCD (vpsfc0) at psfc0
        ! -----------
        ipsfc0=-9
        ipsfc1=-9; ipsfc2=-9; wpsfc1=0.; wpsfc2=0.
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
          write(*,*) " *** Pcld2: Surface Pressure error *** ",psfc0,ix,it
          vpsfc0 = -9999.
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
          go to 333 !goto next name_option_SceneAlbedoAtTerrain 
        else
          vpsfc0=(wpsfc1*vpsfc2+wpsfc2*vpsfc1)/(wpsfc1+wpsfc2)
        endif

        ! -----------------
        ! calculate AMF*VCD for each cloud level using cal_ler_amf
        ! -----------------
        !hqw: note cal_ler_amf may be -9999., they are skipped below
        do ipcld=1,npcld
          aaa = real(cal_ler_amf(ipcld),kind=4)
          !hqw added check for aaa>0.
          !vvcd should always > 0., otherwise it would have been skipped
          if (aaa .gt. 0.) then
             amfvcd(ipcld)=real(aaa*vvcd(ipcld), kind=4)
          else
             ! skipped to next pixel
             amfvcd(ipcld) = -9999.
             ! when this happens, no need to go further, skip calculation
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
             go to 333
          endif
        end do
        ! all amfvcds should >0. from here on

        !hqw initialize local vairable before scd T-correction iteration
        scdm = scdmorg
        scdadj = scdmorg
        t8p = 273.
        temp_t8p = 273.

        iternum = 0

777     continue ! SCD T iteration comes back here

        iflag=-1

        if (scdm .lt. 0.) then ! skip invalid scdm
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
          go to 333
        endif

        if (scdm.le.amfvcd(1)) then !low pressure end
          iflag=0
        endif

        ! find node and weight for scdm, set iflag=1 if found
        ! if not found, iflag=-1 and scdm > amfvcd(npcld): high pressure end
        yy1=-999. ; yy2=-999; ww1=0.; ww2=0.
        do ipcld=1,npcld-1
          if((scdm.gt.amfvcd(ipcld)).and.(scdm.le.amfvcd(ipcld+1))) then
            iflag=1
            yy1=lut_pcld(ipcld)
            yy2=lut_pcld(ipcld+1)
            ww1=scdm-amfvcd(ipcld)
            ww2=amfvcd(ipcld+1)-scdm
          endif
        end do

        ! calculate cpp
        cpp = fFillValue9 ! initialization
        cpp1st = fFillValue9

        if(iflag .ge. 1) then !found node, normal interpolation
          cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)
        else if(iflag .eq. 0) then !low pressure end
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

          ! decrease pressure by 1Pa at a time until minimum diff found
          do ipp=1,100 !safe as pcld(1)=55 < 100
            xx=x1-real(ipp)
!hqw adds condition to exit loop when xx becomes negative
            if (xx .lt. 0.) then
               xx = -9999.
               go to 970 
            endif
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 970
            else
              diff_save=diff
!hqw corrected the following logic
!              if(ipp.le.1) then
               if(ipp .ge. 100) then
                xx=-9999. !no solution found, set to -9999.
              endif
            endif
          end do
970       continue
          cpp=xx 

        else ! iflag .eq. -1: high pressure end
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

          ! increment 1Pa at a time until minimum diff found
          do ipp=1,1000
            xx=x0+real(ipp)
!hqw adds condition to exit loop when xx becomes too large
            if (xx .gt. maxpress) then
               xx = -9999.
               go to 980
            endif
                
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 980
            else
              diff_save=diff
              if(ipp.ge.1000) then
                xx=-9999. !no solution found, set to -9999.
              endif
            endif
          end do
980       continue
          cpp=xx
        endif

        !hqw if calculated cpp < 0., skip calculation
        if (cpp .lt. 0.) then
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
           go to 333
        endif

        if (iternum .eq. 0) then ! iternum has not been incremented yet
           cpp1st=cpp
        endif 

        !hqw scd temperature correction
        !Oct2021 temporarily use T at half cpp, need to tune later
        temp_cpp = cpp * 0.5
        if (temp_cpp .gt. 0. .and. temp_cpp .lt. maxpress) then !reasonable
          ! adjustment always start from scdmorg
          if (name_option_TemperaturePressure .eq. 'GMI') then
            call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else if (name_option_TemperaturePressure .eq. 'GEOS5') then
            call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else ! the following should not happen
            temp_t8p = t8p
          endif
        else ! temp_cpp unreasonable   
           ! set temp_t8p to t8p to terminate iteration below
           temp_t8p = t8p
           ! use cpp1st  & scdmorg
           scdm = scdmorg
           cpp = cpp1st
           ! signal SCD correction problem during processing
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
        endif

        ! increment iternum
        iternum = iternum + 1

        !hqw test if terminate iteration
        delta_temp = abs(t8p - temp_t8p)
        if (delta_temp .ge. dt_threshold .and. iternum .lt. max_scd_iter) then
            t8p = temp_t8p
            scdm = scdadj
            go to 777 ! do another iteration
        endif

        ! signal max_scd_iter reached
        if (iternum .ge. max_scd_iter) then
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
        endif

        !hqw assign output
        if (scdm .gt. 0.) then
           !scdm is right before iteration terminates, scdadj right afterwards
           out_SlantColumnTerrainO2O2(ix,it) = scdm
           out_O2O2TerrainTemperature(ix,it) = t8p
        else
           out_SlantColumnTerrainO2O2(ix,it) = fFillValue9
           out_O2O2TerrainTemperature(ix,it) = fFillValue9
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10) 
        endif

      !+0+0+0+0+0+0+0+0+0+0
      endif !name_option_SceneAlbedoAtTerrain=='yes'//'both'
      !+0+0+0+0+0+0+0+0+0+0

333   continue !skip here if something goes wrong within this name option

      !****************************************************************

      ! assign result
      TerrainLER466=ler466
      TerrainLER440=ler440

      if ((TerrainLER466 .lt. -0.2).or.(TerrainLER466 .gt. 1.2)) then
         TerrainLER466=fFillValue9
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),10)
      endif
      if((TerrainLER466 .ge. -0.2) .and. (TerrainLER466 .lt. 0.0)) TerrainLER466=0.0
      if((TerrainLER466 .gt. 1.0) .and. (TerrainLER466 .le. 1.2)) TerrainLER466=1.0

      if((TerrainLER440 .lt. -0.2) .or. (TerrainLER440 .gt. 1.2)) TerrainLER440=fFillValue9
      if((TerrainLER440 .ge. -0.2) .and. (TerrainLER440 .lt. 0.0)) TerrainLER440=0.0
      if((TerrainLER440 .gt.  1.0) .and. (TerrainLER440 .le. 1.2)) TerrainLER440=1.0

      
      ! re-init ler466, ler440, cpp
      ler466 = fFillValue9
      ler440 = fFillValue9
      cpp = fFillValue9
      cpp1st = fFillValue9
 
      !**************************************************************
      !+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1
      if((name_option_SceneAlbedoAtTerrain.eq.'no') .or. &
           (name_option_SceneAlbedoAtTerrain.eq.'both')) then
      !+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1

         !hqw initialize local variables for name_option
         scdm = fFillValue9
         scdadj = fFillValue9
         temp_cpp = fFillValue9
         delta_temp = fFillValue9
         t8p = 273.
         temp_t8p = 273.

        !--------------------------------------
        ! calculate LER at each pressure level
        !--------------------------------------

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

          ! calculate transmitance and sbar at R=0.0(1),0.1(7), and 0.2(12)
          rad0=real(temp_ler_alb466(1), kind=4)
          rad1=real(temp_ler_alb466(7), kind=4)
          rad2=real(temp_ler_alb466(12), kind=4)
          rrr0=lut_alb(1)
          rrr1=lut_alb(7)
          rrr2=lut_alb(12)
          tran=(1./rrr1-1./rrr2)/(1./(rad1-rad0)-1./(rad2-rad0))
          sbar=1./rrr1-tran/(rad1-rad0)
          !hqw add logic for rad_of_irr
          if (rad_of_irr466(ix,it) .gt. 0.) then
             ler466=(rad_of_irr466(ix,it)-rad0)/(tran+sbar*(rad_of_irr466(ix,it)-rad0))
             if(ler466.lt.0.0) ler466=0.0
             if(ler466.gt.1.0) ler466=1.0
          else
             ler466 = fFillValue9
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
             go to 444
          endif 
          lev_ler_alb466(ipsfc)=ler466
        end do ! ipsfc

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
          if (rad_of_irr440(ix,it) .gt. 0.) then
             ler440=(rad_of_irr440(ix,it)-rad0)/(tran+sbar*(rad_of_irr440(ix,it)-rad0))
             if(ler440.lt.0.0) ler440=0.0
             if(ler440.gt.1.0) ler440=1.0
          else
             ler440 = fFillValue9
          endif
          lev_ler_alb440(ipsfc)=ler440
        end do

        ! -----------------------
        ! AMF at each cloud level using lev_ler_alb466 at each cloud level
        ! -----------------------
        do ipsfc=1,npsfc
          ipcld=ipsfc
          alb0=real(lev_ler_alb466(ipsfc), kind=4)

          ! negative alb0 should have already been skipped before
          ialb1=-9; ialb2=-9 ; walb1=0.; walb2=0.
          do ialb=1,nalb-1
            if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
              ialb1=ialb
              ialb2=ialb+1
              walb1=alb0-lut_alb(ialb)
              walb2=lut_alb(ialb+1)-alb0
            endif
          end do !iab

          !hqw adds logic for iab1
          if (ialb1 .gt. 0) then
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
          else ! ialb1 <= 0
             lev_ler_amf(ipsfc) = -9999.
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
             go to 444
          endif ! ialb1

        end do !ipsfc

        ! -----------------
        ! calculate AMF*VCD
        ! -----------------
        do ipcld=1,npcld
          ipsfc=ipcld
          amfvcd(ipcld)=real(lev_ler_amf(ipsfc), kind=4)*vvcd(ipcld)
        end do

        !hqw re-assign local vairables
        scdm = scdmorg
        scdadj = scdm
        t8p = 273.
        temp_t8p = t8p

        cpp = -9999.
        cpp1st = -9999.

        iternum = 0

776     continue !SCD iteration comes back here

        if (scdm .lt. 0.) then 
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
          go to 444
        endif

        iflag=-1

        if(scdm.le.amfvcd(1)) then ! low pressure end
          iflag=0
        endif

        yy1=-9999.; yy2=-9999.; ww1=0.; ww2=0.
        do ipcld=1,npcld-1
          if((scdm.gt.amfvcd(ipcld)).and.(scdm.le.amfvcd(ipcld+1))) then
            iflag=1
            yy1=lut_pcld(ipcld)
            yy2=lut_pcld(ipcld+1)
            ww1=scdm-amfvcd(ipcld)
            ww2=amfvcd(ipcld+1)-scdm
          endif
        end do

        if(iflag .ge. 1) then !normal interpolation
          cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)
        else if(iflag .eq. 0) then !low pressure end
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
!hqw adds condition to terminate loop when xx<0.
            if (xx .lt. 0.) then
                xx = -9999.
                go to 972
            endif
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 972
            else
              diff_save=diff
!hqw corrected the following logic, so that
!when solution not found within the loop, set xx=-9999.
!              if(ipp.le.1) then
              if (ipp .ge.100) then
                xx=-9999.
              endif
            endif
          end do

972       continue

          cpp=xx

        else ! high pressure end
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
            xx=x0+real(ipp) !increment xx 1Pa at a time
!hqw adds condition to terminate loop when xx too large
            if (xx .gt. maxpress) then
                xx = -9999.
                go to 982
            endif
            yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
                 +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
                 +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
            diff=abs(scdm-yy)
            if(diff.ge.diff_save) then
              go to 982
            else
              diff_save=diff
              if(ipp.ge.1000) then
                xx=-9999.
              endif
            endif
          end do

982       continue

          cpp=xx

        endif

        if (cpp .lt. 0.) then 
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
          go to 444
        endif

        if (iternum .eq. 0) then !iternum has not been incremented yet
           cpp1st = cpp
        endif

        !hqw scd T-correction
        temp_cpp = cpp * 0.5
        if ((temp_cpp .gt. 0.).and.(temp_cpp.lt.maxpress)) then !reasonable
          if (name_option_TemperaturePressure .eq. 'GMI') then
            call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else if (name_option_TemperaturePressure .eq. 'GEOS5') then
            call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
          else ! this should not happen
            temp_t8p = t8p
          endif
        else ! unresonable
           temp_t8p = t8p ! to terminate iteration below
           cpp = cpp1st ! use values before iteration
           scdm = scdmorg
           ! indicate SCD correction problem 
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
        endif

        ! increment iternum
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
            out_SlantColumnSceneO2O2(ix,it)= fFillValue9
            out_O2O2SceneTemperature(ix,it) = fFillValue9
            out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
            go to 444
        endif

        if (iternum .ge. max_scd_iter) then
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
        endif

        !-----------------------
        ! Assign SceneCPP
        !-----------------------
        SceneCPP=real(cpp, kind=4)

        if ((SceneCPP .lt. 0.)) then
         SceneCPP=fFillValue9
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
         go to 444
        endif

        !clip SceneCPP
        if((SceneCPP.gt.psfc0).and.(SceneCPP.le.maxpress)) SceneCPP=psfc0
 
        !-----------------------------
        !calculate ler466 at SceneCPP
        !-----------------------------
        cpp = SceneCPP

        ler466 = fFillValue9

        iflag=-1

        ! cpp < 0. should have been skipped before
        if (cpp .lt. lut_pcld(1)) then ! low pressure end 
           iflag = 0
        endif

        rr1=-9. ; rr2=-9.; wr1=0. ; wr2=0.
        do ipcld=1,npcld-1
          if((cpp.gt.lut_pcld(ipcld)).and.(cpp.le.lut_pcld(ipcld+1))) then
            iflag=1
            rr1=real(lev_ler_alb466(ipcld), kind=4)
            rr2=real(lev_ler_alb466(ipcld+1), kind=4)
            wr1=real(cpp, kind=4)-lut_pcld(ipcld)
            wr2=lut_pcld(ipcld+1)-real(cpp, kind=4)
          endif
        end do

        if(iflag .ge. 1) then ! normal interpolation
          ler466=(wr1*rr2+wr2*rr1)/(wr1+wr2)
        else if (iflag .eq. -1) then ! high pressure end
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
        else if (iflag .eq. 0) then !hqw add low pressure end
          x0=lut_pcld(1)
          x1=lut_pcld(2)
          x2=lut_pcld(3)
          y0=real(lev_ler_alb466(1),kind=4)
          y1=real(lev_ler_alb466(2),kind=4)
          y2=real(lev_ler_alb466(3),kind=4)
          xx=real(cpp,kind=4)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          ler466=yy
       endif
          
        !-------------------------------
        !calculate ler440 at cpp level
        !-------------------------------
        ler440 = fFillValue9

        iflag=-1

        if (cpp .le. lut_pcld(1)) then ! low pressure end
           iflag = 0
        endif

        rr1=-9.; rr2=-9; wr1=0.; wr2=0. 
        do ipcld=1,npcld-1
          if((cpp.gt.lut_pcld(ipcld)).and.(cpp.le.lut_pcld(ipcld+1))) then
            iflag=1
            rr1=real(lev_ler_alb440(ipcld), kind=4)
            rr2=real(lev_ler_alb440(ipcld+1), kind=4)
            wr1=real(cpp, kind=4)-lut_pcld(ipcld)
            wr2=lut_pcld(ipcld+1)-real(cpp, kind=4)
          endif
        end do

        if(iflag .ge. 1) then !normal interpolation
          ler440=(wr1*rr2+wr2*rr1)/(wr1+wr2)
        else if (iflag .eq. -1) then ! high pressure end
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
        else if (iflag .eq. 0) then !hqw add low pressure end
          x0=lut_pcld(1)
          x1=lut_pcld(2)
          x2=lut_pcld(3)
          y0=real(lev_ler_alb440(1),kind=4)
          y1=real(lev_ler_alb440(2),kind=4)
          y2=real(lev_ler_alb440(3),kind=4)
          xx=real(cpp,kind=4)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          ler466=yy
       endif

      !+1+1+1+1+1+1+1+1
      endif !name_option_SceneAlbedoAtTerrain .eq. 'no' // 'both'
      !+1+1+1+1+1+1+1+1

444   continue ! skip here when things go wrong within the name option

      SceneLER466=ler466
      SceneLER440=ler440

      if((SceneLER466 .lt. -0.2) .or. (SceneLER466 .gt. 1.2)) then
         SceneLER466=fFillValue9
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),11)
      endif
      if((SceneLER466 .ge. -0.2) .and. (SceneLER466 .lt. 0.0)) SceneLER466=0.0
      if((SceneLER466 .gt.  1.0) .and. (SceneLER466 .le. 1.2)) SceneLER466=1.0

      if((SceneLER440 .lt. -0.2) .or. (SceneLER440 .gt. 0.2)) SceneLER440=fFillValue9
      if((SceneLER440 .ge. -0.2) .and. (SceneLER440 .lt. 0.0)) SceneLER440=0.0
      if((SceneLER440 .gt.  1.0) .and. (SceneLER440 .le. 0.2)) SceneLER440=1.0

      !--------------
      ! Assign result to array
      !--------------
      out_SurfaceLER466(ix,it)=TerrainLER466
      out_SurfaceLER440(ix,it)=TerrainLER440

      out_ScenePressure(ix,it)=SceneCPP

      out_SceneLER466(ix,it)=SceneLER466
      out_SceneLER440(ix,it)=SceneLER440

      ! ---------------------------------------------------------
      ! Use scene pressure for
      !   1. snow/ice > min_snowice
      ! ---------------------------------------------------------
      if(name_option_SnowIce.eq.'Pscene') then
        !if(btest(out_ProcessingQualityFlags(ix,it),4)) then
        !hqw replce condition
         if((rad_SnowIceFraction(ix,it) .gt. min_snowice).and. &
            (rad_SnowIceFraction(ix,it) .le. 1.0)) then

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
            out_EffectiveCloudFraction(ix,it)=1.
            out_EffectiveCloudFractionNotClipped(ix,it)=1.
            out_CloudRadianceFraction466(ix,it)=1.
            out_CloudRadianceFractionNotClipped466(ix,it)=1.
            out_CloudRadianceFraction440(ix,it)=1.
            out_CloudRadianceFractionNotClipped440(ix,it)=1.
          endif

          !signal pcld replacement by pscene
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),4)
        endif ! min_snowice
      endif !name_option_SnowIce

      !--------------------------------------------------------
      ! Use scene pressure for 
      !   2. ECF < min_ecf
      !--------------------------------------------------------
      ! original default choice was 'yes', now changed to 'no'
      ! so that user has the choice of if/how to clip
      if(name_option_MinECF.eq.'yes') then
!        if(btest(out_ProcessingQualityFlags(ix,it),2)) then 
        if ((out_EffectiveCloudFraction(ix,it).gt.0.).and. &
            (out_EffectiveCloudFraction(ix,it).lt. min_ecf)) then 
          out_CloudPressure(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
          out_CloudPressureNotClipped(ix,it)=nint(out_ScenePressure(ix,it), kind=2)
          !clip
          if((out_CloudPressure(ix,it).gt.0).and.(out_CloudPressure(ix,it).le.100)) &
             out_CloudPressure(ix,it)=100
          if((out_CloudPressure(ix,it).ge.nint(psfc0)).and.(out_CloudPressure(ix,it).lt.maxpress)) &
             out_CloudPressure(ix,it)=nint(psfc0, kind=2)
          !signal pcld replacement by pscene
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),2)
        endif ! min_ecf
      endif !name_option_MinECF

990   continue

      !=====
    end do
  end do
  !=====

  !************************
end subroutine cal_pscene
!************************
end module m_cal_pscene
