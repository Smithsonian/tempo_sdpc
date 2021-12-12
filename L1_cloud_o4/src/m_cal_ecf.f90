module m_cal_ecf
  public cal_ecf
contains
!******************
subroutine cal_ecf
  !******************
  use m_vars
  use m_read_GMI
  use m_read_DEM
  use m_read_hdf5

  implicit none

  real::temp_raa,temp_psfc,temp_rsfc
  real::rad466,rad477,rad440
  real(kind=4)::rout_ecf,rout_crf440,rout_crf466
  integer::ialb, isza, ivza, iraa, ipsfc
  integer::ialb1,isza1,ivza1,iraa1,ipsfc1
  integer::ialb2,isza2,ivza2,iraa2,ipsfc2
  real::   walb1,wsza1,wvza1,wraa1,wpsfc1
  real::   walb2,wsza2,wvza2,wraa2,wpsfc2
  integer(kind=4)::ierr
  integer(kind=4)::nt,nx
  integer(kind=4)::it,ix

  integer(kind=4)::gmi_ix1,gmi_ix2,gmi_iy1,gmi_iy2
  real::gmi_wx1,gmi_wx2,gmi_wy1,gmi_wy2
  real::pp11,pp12,pp21,pp22,pp1,pp2

  integer(kind=4)::kleipool_ix,kleipool_iy

  real::r11111,r11112,r11121,r11122,r11211,r11212,r11221,r11222,r12111,r12112,r12121,r12122,r12211,r12212,r12221,r12222
  real::r21111,r21112,r21121,r21122,r21211,r21212,r21221,r21222,r22111,r22112,r22121,r22122,r22211,r22212,r22221,r22222
  real::r1111,r1112,r1121,r1122,r1211,r1212,r1221,r1222,r2111,r2112,r2121,r2122,r2211,r2212,r2221,r2222
  real::r111,r112,r121,r122,r211,r212,r221,r222
  real::r11,r12,r21,r22
  real::r1,r2

  real::pi,dtor
  real(kind=4), parameter:: fspecial = -9999. !make this a large negative value

  ! ------
  ! initialization
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.

  nt=rad_NumTimes
  nx=rad_nXtrack


  ! allocate arrays & fill values with fspecial

  allocate(cal_rad_clr(nx,nt),stat=ierr)
  allocate(cal_rad_cld(nx,nt),stat=ierr)
  allocate(cal_rad_cld440(nx,nt),stat=ierr)

  allocate(rad_of_irr440(nx,nt),stat=ierr)
  allocate(rad_of_irr466(nx,nt),stat=ierr)
  allocate(rad_of_irr477(nx,nt),stat=ierr)

  cal_rad_clr=fspecial
  cal_rad_cld=fspecial
  cal_rad_cld440=fspecial
  rad_of_irr440=fspecial
  rad_of_irr466=fspecial
  rad_of_irr477=fspecial

!  hqw STDs are not calculated in OMCDO2N, thus disabled
  allocate(out_EffectiveCloudFraction(nx,nt),stat=ierr)
  allocate(out_EffectiveCloudFractionNotClipped(nx,nt),stat=ierr)
!  allocate(out_EffectiveCloudFractionSTD(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFraction440(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFractionNotClipped440(nx,nt),stat=ierr)
!  allocate(out_CloudRadianceFractionSTD440(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFraction466(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFractionNotClipped466(nx,nt),stat=ierr)
!  allocate(out_CloudRadianceFractionSTD466(nx,nt),stat=ierr)
  allocate(out_RelativeAzimuthAngle(nx,nt),stat=ierr)
  allocate(out_ReflectanceFactor(nx,nt),stat=ierr)

  out_EffectiveCloudFraction=fspecial
  out_EffectiveCloudFractionNotClipped=fspecial
!  out_EffectiveCloudFractionSTD=int(iFillValue, kind=2)
  out_CloudRadianceFraction440=fspecial
  out_CloudRadianceFractionNotClipped440=fspecial
!  out_CloudRadianceFractionSTD440=int(iFillValue, kind=2)
  out_CloudRadianceFraction466=fspecial
  out_CloudRadianceFractionNotClipped466=fspecial
!  out_CloudRadianceFractionSTD466=int(iFillValue, kind=2)

  out_RelativeAzimuthAngle=fspecial
  out_ReflectanceFactor=fspecial

  !hqw moved GMI lat/lon to read_GMI_TMP and read_DEM_GMI_TMP
  !-----------------
  ! read GMI lat/lon
  !-----------------
   !hqw debug
   if ((ixdebug .gt. 0) .and. (ixdebug .le. nx) .and. (itdebug .gt. 0)) then
      close(19)
      write(*,*) 'writing debug_ecf.txt'
      open(unit=19, file='debug_ecf.txt')
      write(19,*)'ix, alb0,  psfc0,  rad_of_irr466, cal_rad_clr, cal_rad_cld,cldfrac'
   endif

  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========

      !initialize to special negative values
      rout_ecf =fspecial
      rout_crf440 =fspecial
      rout_crf466 =fspecial

      pflag00=0
      pflag01=0
      if((rad_Latitude(ix,it) .lt. -90.) .or. (rad_Latitude(ix,it) .gt. 90.) .or. &
        (rad_Longitude(ix,it) .lt. -180.) .or. (rad_Longitude(ix,it) .gt. 180.)) then
        pflag00=pflag00+1
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),0)
      endif
      !hqw changed 86. to max_SZA
      if((rad_SolarZenithAngle(ix,it) .lt. 0.) .or. (rad_SolarZenithAngle(ix,it) .gt. max_SZA)) then
        pflag01=pflag01+1
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),1)
      endif
      !hqw changed 72. to max_VZA
      if((rad_ViewingZenithAngle(ix,it) .lt. 0.) .or. (rad_ViewingZenithAngle(ix,it) .gt. max_VZA)) then
        pflag01=pflag01+1
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),1)
      endif

      if((rad_SolarAzimuthAngle(ix,it) .ge. -360.) .and. (rad_SolarAzimuthAngle(ix,it) .le. 360.) .and. &
           (rad_ViewingAzimuthAngle(ix,it) .ge. -360.) .and. (rad_ViewingAzimuthAngle(ix,it) .le. 360.)) then
      !hqw RAA = SAA - VAA + PI, Why +PI?
      !xliu: this is related to how the SAA and VAA are fined
      !      RAA of forward scattering = 0, RAA of backward scattering = 180.
      !also see Eun-su Yang email slide for explanation
        temp_raa=rad_SolarAzimuthAngle(ix,it)+180.0-rad_ViewingAzimuthAngle(ix,it)
        ! ensure temp_raa is within [0., 360.) range
        do while((temp_raa .lt. 0.0) .or. (temp_raa .ge. 360.0))
          if(temp_raa .ge. 360.0) temp_raa=temp_raa-360.
          if(temp_raa .lt. 0.0) temp_raa=temp_raa+360.
        end do
      else
        temp_raa= fspecial
        pflag01 = pflag01 + 1
        ! temp_raa = -9999. will be skipped in calculation
      endif

      !hqw skip calculation if lat//angles are invalid
      if((pflag00 .ge. 1) .or. (pflag01 .ge. 1)) go to 990

      ! get local radiances
      rad466=rad_466nm(ix,it)
      rad477=rad_477nm(ix,it)
      rad440=rad_440nm(ix,it)

      ! ------------------------
      ! calculate cloud fraction
      ! ------------------------
      ! hqw alb0 & psfc0 were originally initialized to OMCLDO2
      ! but will be overwritten later in this subroutine by climatology
      ! to get rid of dependency on OMCLDO2 , set to temporary values
       temp_psfc = 1000.
       temp_rsfc = 0.03

      ! radiance unavailable at 440 or 466, set out_ProcessingQualityFlags bit 7
      if((rad440 .le. 0.) .or. (rad466 .le. 0.)) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
        go to 990
      endif

      !hqw izoom=0 for TEMPO in m_vars.f90 as a parameter
      if (irr_out_irradiance_466nm(ix) .gt. 0.) then
         rad_of_irr466(ix,it)=rad466/irr_out_irradiance_466nm(ix)*(rad_EarthSunDist/irr_EarthSunDist)**2
      else
         rad_of_irr466(ix,it) = fspecial
      endif
      !hqw skip calculation if rad_of_irr < 0.
      if (rad_of_irr466(ix,it) .lt. 0.) go to 990

      if (irr_out_irradiance_440nm(ix) .gt. 0.) then
         rad_of_irr440(ix,it)=rad440/irr_out_irradiance_440nm(ix)*(rad_EarthSunDist/irr_EarthSunDist)**2
      else
         rad_of_irr440(ix,it) = fspecial
      endif
      if (rad_of_irr440(ix,it) .lt. 0.) go to 990

      !hqw 477nm is currently not used in cloud fraction calculation
      ! thus even if rad_of_irr477< 0., still do the rest of the calculation
      if (rad477 .gt. 0.) then
         rad_of_irr477(ix,it)=rad477/irr_out_irradiance_477nm(ix)*(rad_EarthSunDist/irr_EarthSunDist)**2
      else
         rad_of_irr477(ix,it) = fspecial
      endif

      sza0=rad_SolarZenithAngle(ix,it)
      vza0=rad_ViewingZenithAngle(ix,it)
      raa0=temp_raa
      !hqw Why 360.-raa?
      !xliu: +raa has the same effect as -raa,
      !      and RAA needs to be within [0.,180] for use with LUT
      !also consult E.Yang's slides for definition of RAA
      if (raa0 .lt. 0.) raa0 = - raa0
      !hqw added the line above to make sure raa0>0. because of LUT RAA range
      !though this may be redundant as temp_raa should be within [0..,360.)
      if (raa0 .gt. 180.) raa0=360.-raa0
      out_RelativeAzimuthAngle(ix,it)=raa0

      !initialize alb0 and psfc0 to temporary value
      alb0=temp_rsfc
      psfc0=temp_psfc

      !----------------
      !get actual psfc0
      !out-of-range rad_Longitude should have been skipped
      if(name_option_TemperaturePressure.eq.'GMI') then
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

        pp11=gmi_TerrainPressure(gmi_ix1,gmi_iy1)
        pp12=gmi_TerrainPressure(gmi_ix1,gmi_iy2)
        pp21=gmi_TerrainPressure(gmi_ix2,gmi_iy1)
        pp22=gmi_TerrainPressure(gmi_ix2,gmi_iy2)
        pp1=(gmi_wy2*pp11+gmi_wy1*pp12)/(gmi_wy1+gmi_wy2)
        pp2=(gmi_wy2*pp21+gmi_wy1*pp22)/(gmi_wy1+gmi_wy2)
        gmi_psfc=(gmi_wx2*pp1+gmi_wx1*pp2)/(gmi_wx1+gmi_wx2)

        psfc0=gmi_psfc
        !hqw assign l2_TerrainPressure for GMI here
        l2_TerrainPressure(ix,it) = psfc0
      endif

      !hqw do not use BDEM without further development
      !if(name_option_TemperaturePressure.eq.'BDEM') then
      !  psfc0=BDEM_TerrainPressure(ix,it)
      !  l2_TerrainPressure(ix,it) = psfc0
      !endif

      if(name_option_TemperaturePressure.eq.'GEOS5') then
        !hqw l2_TerrainPressure already assigned in read_geoscf
        ! which has been called before this subroutine, use it as psfc0
        psfc0=l2_TerrainPressure(ix,it)
      endif

      !hqw if psfc0[hPa] is out of range, skip calculation
      if ((psfc0 .lt. 0.) .or. (psfc0 .gt. 1200.)) goto 990

      !---------------
      ! get actual alb0
      if(name_option_SurfaceReflectivity.eq.'Kleipool') then
        kleipool_ix=nint((rad_Longitude(ix,it)+180.0)/0.5)
        kleipool_iy=nint((rad_Latitude(ix,it)+90.0)/0.5)
        if(kleipool_ix.lt.1) kleipool_ix=1
        if(kleipool_ix.gt.kleipool_nx) kleipool_ix=kleipool_nx
        if(kleipool_iy.lt.1) kleipool_iy=1
        if(kleipool_iy.gt.kleipool_ny) kleipool_iy=kleipool_ny
        alb0=kleipool_SurfaceReflectivity466(kleipool_ix,kleipool_iy)
      endif

      if(name_option_SurfaceReflectivity.eq.'BRDF') then
        alb0=BRDF_SurfaceReflectivity466(ix,it)
      endif

      ! bound alb0 within [0.,1.] if it is in reasonable range
      ! otherwise skip the calculation
      if ((alb0 .lt. -0.2) .or. (alb0 .gt. 1.2)) go to 990
      if((alb0 .ge. -0.2) .and. (alb0 .lt. 0.0)) alb0=0.0
      if((alb0 .gt.  1.0) .and. (alb0 .le. 1.2)) alb0=1.0
 
      !-------------------
      ! interpolation prep for alb/sza/vza/raa/psfc

      ! initialize local weight for interpolation
      walb1 = fspecial
      walb2 = fspecial
      wsza1 = fspecial
      wsza2 = fspecial
      wvza1 = fspecial
      wvza2 = fspecial
      wraa1 = fspecial
      wraa2 = fspecial
      wpsfc1 = fspecial
      wpsfc2 = fspecial
 
      ! if interpolation nodes cannot be found skip calculation
      ! find nodes for alb0
      ialb1=-9; ialb2=-9
      do ialb=1,nalb-1
        if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
          ialb1=ialb
          ialb2=ialb+1
          walb1=alb0-lut_alb(ialb)
          walb2=lut_alb(ialb+1)-alb0
        endif
      end do
      if ((ialb1 .lt. 0) .or. (ialb2 .lt. 0) .or. (walb1 .lt. 0.) &
          .or. (walb2 .lt. 0.)) go to 990

      ! find nodes for sza0
      isza1=-9; isza2=-9
      do isza=1,nsza-1
        if((sza0 .ge. lut_sza(isza)) .and. (sza0 .le. lut_sza(isza+1))) then
          isza1=isza
          isza2=isza+1
          wsza1=sza0-lut_sza(isza)
          wsza2=lut_sza(isza+1)-sza0
        endif
      end do
      if ((isza1 .lt. 0) .or. (isza2 .lt. 0) .or. (wsza1 .lt. 0.) &
          .or. (wsza2 .lt. 0)) go to 990

      ! find nodes for vza0
      ivza1=-9; ivza2=-9
      do ivza=1,nvza-1
        if((vza0 .ge. lut_vza(ivza)) .and. (vza0 .le. lut_vza(ivza+1))) then
          ivza1=ivza
          ivza2=ivza+1
          wvza1=vza0-lut_vza(ivza)
          wvza2=lut_vza(ivza+1)-vza0
        endif
      end do
      if ((ivza1 .lt. 0) .or. (ivza2 .lt. 0) .or. (wvza1 .lt. 0.) &
          .or. (wvza2 .lt. 0.)) go to 990

      ! find nodes for raa0
      iraa1=-9; iraa2=-9
      do iraa=1,nraa-1
        if((raa0 .ge. lut_raa(iraa)) .and. (raa0 .le. lut_raa(iraa+1))) then
          iraa1=iraa
          iraa2=iraa+1
          wraa1=raa0-lut_raa(iraa)
          wraa2=lut_raa(iraa+1)-raa0
        endif
      end do
      if ((iraa1 .lt. 0) .or. (iraa2 .lt. 0) .or. (wraa1 .lt. 0.) &
           .or. (wraa2 .lt. 0.)) go to 990

      !find nodes for psfc0
      !bound psfc by the largest lut_psfc
      if (psfc0 .gt. lut_psfc(npsfc)) psfc0=lut_psfc(npsfc)
      ipsfc1=-9; ipsfc2=-9
      do ipsfc=1,npsfc-1
        if((psfc0 .ge. lut_psfc(ipsfc)) .and. (psfc0 .le. lut_psfc(ipsfc+1))) then
          ipsfc1=ipsfc
          ipsfc2=ipsfc+1
          wpsfc1=psfc0-lut_psfc(ipsfc)
          wpsfc2=lut_psfc(ipsfc+1)-psfc0
        endif
      end do
      if ((ipsfc1 .lt. 0) .or. (ipsfc2 .lt. 0)) then
        write(*,*) " *** Surface Pressure too small *** ",ix,it
        go to 990
      endif
      if ((wpsfc1 .lt. 0.) .or. (wpsfc2 .lt. 0.)) go to 990

      !-------------
      ! get LUT values at interpolation nodes
      !-------------
      !------------------------------
      !radiance at surface: clear sky
      !------------------------------
      r11111=lut_rad_clr(ialb1,isza1,ivza1,iraa1,ipsfc1)
      r11112=lut_rad_clr(ialb1,isza1,ivza1,iraa1,ipsfc2)
      r11121=lut_rad_clr(ialb1,isza1,ivza1,iraa2,ipsfc1)
      r11122=lut_rad_clr(ialb1,isza1,ivza1,iraa2,ipsfc2)
      r11211=lut_rad_clr(ialb1,isza1,ivza2,iraa1,ipsfc1)
      r11212=lut_rad_clr(ialb1,isza1,ivza2,iraa1,ipsfc2)
      r11221=lut_rad_clr(ialb1,isza1,ivza2,iraa2,ipsfc1)
      r11222=lut_rad_clr(ialb1,isza1,ivza2,iraa2,ipsfc2)
      r12111=lut_rad_clr(ialb1,isza2,ivza1,iraa1,ipsfc1)
      r12112=lut_rad_clr(ialb1,isza2,ivza1,iraa1,ipsfc2)
      r12121=lut_rad_clr(ialb1,isza2,ivza1,iraa2,ipsfc1)
      r12122=lut_rad_clr(ialb1,isza2,ivza1,iraa2,ipsfc2)
      r12211=lut_rad_clr(ialb1,isza2,ivza2,iraa1,ipsfc1)
      r12212=lut_rad_clr(ialb1,isza2,ivza2,iraa1,ipsfc2)
      r12221=lut_rad_clr(ialb1,isza2,ivza2,iraa2,ipsfc1)
      r12222=lut_rad_clr(ialb1,isza2,ivza2,iraa2,ipsfc2)
      r21111=lut_rad_clr(ialb2,isza1,ivza1,iraa1,ipsfc1)
      r21112=lut_rad_clr(ialb2,isza1,ivza1,iraa1,ipsfc2)
      r21121=lut_rad_clr(ialb2,isza1,ivza1,iraa2,ipsfc1)
      r21122=lut_rad_clr(ialb2,isza1,ivza1,iraa2,ipsfc2)
      r21211=lut_rad_clr(ialb2,isza1,ivza2,iraa1,ipsfc1)
      r21212=lut_rad_clr(ialb2,isza1,ivza2,iraa1,ipsfc2)
      r21221=lut_rad_clr(ialb2,isza1,ivza2,iraa2,ipsfc1)
      r21222=lut_rad_clr(ialb2,isza1,ivza2,iraa2,ipsfc2)
      r22111=lut_rad_clr(ialb2,isza2,ivza1,iraa1,ipsfc1)
      r22112=lut_rad_clr(ialb2,isza2,ivza1,iraa1,ipsfc2)
      r22121=lut_rad_clr(ialb2,isza2,ivza1,iraa2,ipsfc1)
      r22122=lut_rad_clr(ialb2,isza2,ivza1,iraa2,ipsfc2)
      r22211=lut_rad_clr(ialb2,isza2,ivza2,iraa1,ipsfc1)
      r22212=lut_rad_clr(ialb2,isza2,ivza2,iraa1,ipsfc2)
      r22221=lut_rad_clr(ialb2,isza2,ivza2,iraa2,ipsfc1)
      r22222=lut_rad_clr(ialb2,isza2,ivza2,iraa2,ipsfc2)

      ! -----------------
      ! in case wsfc2=0.0
      ! -----------------
      if((r11111.lt.0.0) .or. (r11112.lt.0.0)) then
        cal_rad_clr(ix,it)= fspecial
        go to 897
      endif

      r1111=(wpsfc2*r11111+wpsfc1*r11112)/(wpsfc1+wpsfc2)
      r1112=(wpsfc2*r11121+wpsfc1*r11122)/(wpsfc1+wpsfc2)
      r1121=(wpsfc2*r11211+wpsfc1*r11212)/(wpsfc1+wpsfc2)
      r1122=(wpsfc2*r11221+wpsfc1*r11222)/(wpsfc1+wpsfc2)
      r1211=(wpsfc2*r12111+wpsfc1*r12112)/(wpsfc1+wpsfc2)
      r1212=(wpsfc2*r12121+wpsfc1*r12122)/(wpsfc1+wpsfc2)
      r1221=(wpsfc2*r12211+wpsfc1*r12212)/(wpsfc1+wpsfc2)
      r1222=(wpsfc2*r12221+wpsfc1*r12222)/(wpsfc1+wpsfc2)
      r2111=(wpsfc2*r21111+wpsfc1*r21112)/(wpsfc1+wpsfc2)
      r2112=(wpsfc2*r21121+wpsfc1*r21122)/(wpsfc1+wpsfc2)
      r2121=(wpsfc2*r21211+wpsfc1*r21212)/(wpsfc1+wpsfc2)
      r2122=(wpsfc2*r21221+wpsfc1*r21222)/(wpsfc1+wpsfc2)
      r2211=(wpsfc2*r22111+wpsfc1*r22112)/(wpsfc1+wpsfc2)
      r2212=(wpsfc2*r22121+wpsfc1*r22122)/(wpsfc1+wpsfc2)
      r2221=(wpsfc2*r22211+wpsfc1*r22212)/(wpsfc1+wpsfc2)
      r2222=(wpsfc2*r22221+wpsfc1*r22222)/(wpsfc1+wpsfc2)

      r111=(wraa2*r1111+wraa1*r1112)/(wraa1+wraa2)
      r112=(wraa2*r1121+wraa1*r1122)/(wraa1+wraa2)
      r121=(wraa2*r1211+wraa1*r1212)/(wraa1+wraa2)
      r122=(wraa2*r1221+wraa1*r1222)/(wraa1+wraa2)
      r211=(wraa2*r2111+wraa1*r2112)/(wraa1+wraa2)
      r212=(wraa2*r2121+wraa1*r2122)/(wraa1+wraa2)
      r221=(wraa2*r2211+wraa1*r2212)/(wraa1+wraa2)
      r222=(wraa2*r2221+wraa1*r2222)/(wraa1+wraa2)

      r11=(wvza2*r111+wvza1*r112)/(wvza1+wvza2)
      r12=(wvza2*r121+wvza1*r122)/(wvza1+wvza2)
      r21=(wvza2*r211+wvza1*r212)/(wvza1+wvza2)
      r22=(wvza2*r221+wvza1*r222)/(wvza1+wvza2)

      r1=(wsza2*r11+wsza1*r12)/(wsza1+wsza2)
      r2=(wsza2*r21+wsza1*r22)/(wsza1+wsza2)

      cal_rad_clr(ix,it)=(walb2*r1+walb1*r2)/(walb1+walb2)

897   continue

      !--------------------------------
      !radiance at 700 hPa: cloudy sky
      !hqw in LUT_4660_RAD.h5, ALB(18)=0.8, Psfc(18)=701hPa
      !--------------------------------
      ialb= LUT466rad_cloud_albid !hqw remove hardcoded index 18
      ipsfc= LUT466rad_cloud_psfcid !18

      r111=lut_rad_clr(ialb,isza1,ivza1,iraa1,ipsfc)
      r112=lut_rad_clr(ialb,isza1,ivza1,iraa2,ipsfc)
      r121=lut_rad_clr(ialb,isza1,ivza2,iraa1,ipsfc)
      r122=lut_rad_clr(ialb,isza1,ivza2,iraa2,ipsfc)
      r211=lut_rad_clr(ialb,isza2,ivza1,iraa1,ipsfc)
      r212=lut_rad_clr(ialb,isza2,ivza1,iraa2,ipsfc)
      r221=lut_rad_clr(ialb,isza2,ivza2,iraa1,ipsfc)
      r222=lut_rad_clr(ialb,isza2,ivza2,iraa2,ipsfc)
      r11=(wraa2*r111+wraa1*r112)/(wraa1+wraa2)
      r12=(wraa2*r121+wraa1*r122)/(wraa1+wraa2)
      r21=(wraa2*r211+wraa1*r212)/(wraa1+wraa2)
      r22=(wraa2*r221+wraa1*r222)/(wraa1+wraa2)
      r1=(wvza2*r11+wvza1*r12)/(wvza1+wvza2)
      r2=(wvza2*r21+wvza1*r22)/(wvza1+wvza2)
      cal_rad_cld(ix,it)=(wsza2*r1+wsza1*r2)/(wsza1+wsza2)

      !--------------------------------
      !radiance at 700 hPa: cloudy sky
      !in LUT_4400_RAD.h5, ALB(18)=0.8, Psfc(18)=701hPa
      !--------------------------------
      ialb=LUT440rad_cloud_albid ! 18
      ipsfc=LUT440rad_cloud_psfcid !18

      r111=lut_rad_clr440(ialb,isza1,ivza1,iraa1,ipsfc)
      r112=lut_rad_clr440(ialb,isza1,ivza1,iraa2,ipsfc)
      r121=lut_rad_clr440(ialb,isza1,ivza2,iraa1,ipsfc)
      r122=lut_rad_clr440(ialb,isza1,ivza2,iraa2,ipsfc)
      r211=lut_rad_clr440(ialb,isza2,ivza1,iraa1,ipsfc)
      r212=lut_rad_clr440(ialb,isza2,ivza1,iraa2,ipsfc)
      r221=lut_rad_clr440(ialb,isza2,ivza2,iraa1,ipsfc)
      r222=lut_rad_clr440(ialb,isza2,ivza2,iraa2,ipsfc)
      r11=(wraa2*r111+wraa1*r112)/(wraa1+wraa2)
      r12=(wraa2*r121+wraa1*r122)/(wraa1+wraa2)
      r21=(wraa2*r211+wraa1*r212)/(wraa1+wraa2)
      r22=(wraa2*r221+wraa1*r222)/(wraa1+wraa2)
      r1=(wvza2*r11+wvza1*r12)/(wvza1+wvza2)
      r2=(wvza2*r21+wvza1*r22)/(wvza1+wvza2)
      cal_rad_cld440(ix,it)=(wsza2*r1+wsza1*r2)/(wsza1+wsza2)

      !-----------------------------------
      !calculate cloud fraction, ecf & crf
      !-----------------------------------
      !hqw added condition
      if ((cal_rad_clr(ix,it) .gt. 0.) .and. (cal_rad_cld(ix,it) .gt. 0.)) then
         rout_ecf=(rad_of_irr466(ix,it)-cal_rad_clr(ix,it))/(cal_rad_cld(ix,it)-cal_rad_clr(ix,it))
      else
         rout_ecf = fspecial
      endif

      ! caculate cloud radiance fraction at 440
      if ((cal_rad_cld440(ix,it) .gt. 0.) .and. (rad_of_irr440(ix,it).gt. 0.)) then
         rout_crf440=rout_ecf*cal_rad_cld440(ix,it)/rad_of_irr440(ix,it)
      else
        rout_crf440 = fspecial
      endif

      ! calculate cloud radiance fraction at 466
      if ((cal_rad_cld(ix,it) .gt. 0.).and.(rad_of_irr466(ix,it).gt. 0.)) then
         rout_crf466=rout_ecf*cal_rad_cld(ix,it)/rad_of_irr466(ix,it)
      else
         rout_crf466 = fspecial
      endif

      out_EffectiveCloudFractionNotClipped(ix,it)= rout_ecf
      out_CloudRadianceFractionNotClipped440(ix,it)= rout_crf440
      out_CloudRadianceFractionNotClipped466(ix,it)= rout_crf466

      ! clip cloud fraction
      !hqw added logic to differentiate skipped or bad calculation
      if((rout_ecf .lt. 0.) .and. (rout_ecf .ge. -0.2)) rout_ecf=0.
      if((rout_ecf .gt. 1.) .and. (rout_ecf .le. 1.2)) rout_ecf=1.
      if((rout_ecf .lt. -0.2) .or. (rout_ecf .gt. 1.2)) rout_ecf=fspecial 
      if((rout_crf440 .lt. 0.).and.(rout_crf440 .ge. -0.2)) rout_crf440=0.
      if((rout_crf440 .gt. 1.).and.(rout_crf440 .le. 1.2)) rout_crf440=1.
      if((rout_crf440 .lt. -0.2) .or. (rout_crf440 .gt. 1.2)) rout_crf440=fspecial
      if((rout_crf466 .lt. 0.).and.(rout_crf466 .ge. -0.2)) rout_crf466=0.
      if((rout_crf466 .gt. 1.).and.(rout_crf466 .le. 1.2)) rout_crf466=1.
      if((rout_crf466 .lt. -0.2) .or. (rout_crf466 .gt. 1.2)) rout_crf466=fspecial

      out_EffectiveCloudFraction(ix,it)=rout_ecf
      out_CloudRadianceFraction440(ix,it)=rout_crf440
      out_CloudRadianceFraction466(ix,it)=rout_crf466

      !hqw out_Reflectance is equivalent to observed Lambertian reflectance at 466
      out_ReflectanceFactor(ix,it)=pi*rad_of_irr466(ix,it)/cos(dtor*sza0)
      if(rad_of_irr466(ix,it).le.0.0) out_ReflectanceFactor(ix,it)=fFillValue

990   continue

      !hqw debug
      if (it .eq. itdebug) then
         write(19,*)ix,alb0,psfc0,rad_of_irr466(ix,it),cal_rad_clr(ix,it),cal_rad_cld(ix,it),out_EffectiveCloudFraction(ix,it)
      endif

      !=====
    end do
  end do
  !=====

  close(19)

!**********************
end subroutine cal_ecf
!**********************
end module m_cal_ecf
