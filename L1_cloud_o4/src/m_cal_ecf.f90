module m_cal_ecf
  public cal_ecf

contains
!11111111111111111111
!******************
subroutine cal_ecf(ecfocp_iternum)
  !******************
!11111111111111111111
  use m_vars
  use m_read_GMI
  use m_read_hdf5
  use m_read_input_kleipool

  implicit none

! input
  integer, intent(in):: ecfocp_iternum

  real(kind=4):: rad466,rad440,rad477
  real(kind=4):: irr466,irr440,irr477

  real:: xradofirr466, yradofirr466

  real(kind=4):: rout_ecf,rout_crf440,rout_crf466,rout_crf477
  integer::ialb, isza, ivza, iraa, ipsfc
  integer::ialb1,isza1,ivza1,iraa1,ipsfc1
  integer::ialb2,isza2,ivza2,iraa2,ipsfc2
  real::   walb1,wsza1,wvza1,wraa1,wpsfc1
  real::   walb2,wsza2,wvza2,wraa2,wpsfc2

  integer(kind=4)::nt,nx
  integer(kind=4)::it,ix

! perturbation polynomial order counter for rad_of_irr466
  integer(kind=4):: iord, iord1

  real:: gmi_psfc

  real(kind=4):: earthsunfactor2
 
  integer(kind=4):: pflag00, pflag01

! local value moved from m_vars
  real:: alb0, sza0, vza0, raa0, psfc0
  
! add local variable
  real::lat0, lon0
  real:: alb440
  real:: kleipool466, kleipool440, kleipool477
  real:: thisecf, thatecf, ecf_change, thisscd

  real::r11111,r11112,r11121,r11122,r11211,r11212,r11221,r11222,r12111,r12112,r12121,r12122,r12211,r12212,r12221,r12222
  real::r21111,r21112,r21121,r21122,r21211,r21212,r21221,r21222,r22111,r22112,r22121,r22122,r22211,r22212,r22221,r22222
  real::r1111,r1112,r1121,r1122,r1211,r1212,r1221,r1222,r2111,r2112,r2121,r2122,r2211,r2212,r2221,r2222
  real::r111,r112,r121,r122,r211,r212,r221,r222
  real::r11,r12,r21,r22
  real::r1,r2

  real:: this_cal_radcld, that_cal_radcld, thiswp, thatwp, thisocp
  integer:: ipsnext

  real::pi,dtor
  real(kind=4) :: fspecial, fspecial9 

  ! ------
  ! initialization
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.
  fspecial = fFillValue ! large negative value in m_vars
  fspecial9 = -999.

  nt=rad_NumTimes
  nx=rad_nXtrack

   if ((trim(run_mode) .eq. 'development').and.(itdebug .ge. 0)) then
      write(*,*) 'writing debug_ecf.txt'
      open(unit=lun_debug_ecf, file='debug_ecf.txt')
      write(lun_debug_ecf,*)'ix, alb0,  psfc0,  rad_of_irr466, cal_rad_clr, cal_rad_cld, cldfrac'
   endif

  ! earthsunfactor2 accounts for earth-sun distance between irr and rad
   earthsunfactor2 = (irr_EarthSunDist/rad_EarthSunDist)**2


  ! loop through each ground pixel
  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========
      ! initial local variables
      ipsnext = -9
      thisocp = fspecial9
      this_cal_radcld = fspecial9
      that_cal_radcld = fspecial9
      thiswp = fspecial9
      thatwp = fspecial9

      ! thisscd has been normalized and filtered in read_cldo4_tio
      ! if thisscd is invalid, only one pass is needed
      ! because ocp will be missing 
      thisscd = nasa_SlantColumnAmountO2O2(ix,it)
      if ((thisscd .lt. 0.).and.(ecfocp_iternum .gt. 1)) then
         ! first iteration will not end up here
         go to 3455 ! skip calculation
      endif   

      ! if ecf_change is below threshold, no need to recalculate
      ! previous iteration values and quality flags are still valid
      ! simply skip to next ground pixel

      if (ecfocp_iternum .gt. 2) then ! check only after 2 passes
         thisecf = out_EffectiveCloudFractionNotClipped(ix,it)
         thatecf = prev_ecf_notclipped(ix,it)
         ecf_change = abs(thisecf - thatecf)
         if ((thisecf .gt. -1.) .and. (thatecf .gt. -1.).and. &
             (ecf_change .lt. delta_ecf)) then
             ! no ProcessingQualityFlags change here
             ! keep previous values and flags, skip calculation
             go to 3455
          else
             ! assign current ecf to previous ecf for next ecfocp iter
             ! current ecf will be re-calculated
             ! and compared with previous ecf next time around
             prev_ecf_notclipped(ix,it)=out_EffectiveCloudFractionNotClipped(ix,it)
             prev_ecf(ix,it) = out_EffectiveCloudFraction(ix,it)
             ! update ecf_niter, skipped calculation will not end up here
             ecf_niter(ix,it) = ecfocp_iternum
          endif
       else
          ! 1st & 2nd iterations will end up here
          ! update ecf_niter for the 1st couple of iter
          ecf_niter(ix,it) = ecfocp_iternum
       endif

      ! the first pass will always go through the whole thing

      !---------------
      ! initialize out_ProcessingQualityFalgs relavent bits to zero
      ! bit7 & bit8 were set before in m_read_input_tio, however,
      ! rad_of_irr466 & rad_of_irr440 are checked again here, 
      ! for ecfocp iteration, it is easier to clear them here as well. 
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),0)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),1)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),3)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),7)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),8)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),9)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),12)
      
      !initialize local variable
      rout_ecf =fspecial
      rout_crf440 =fspecial
      rout_crf466 =fspecial
      rout_crf477 =fspecial

      psfc0=fspecial
      alb0=fspecial ! 466nm
      alb440=fspecial ! 440nm

      ! get local location and angles
      lat0=rad_Latitude(ix,it)
      lon0=rad_Longitude(ix,it)
      sza0=rad_SolarZenithAngle(ix,it)
      vza0=rad_ViewingZenithAngle(ix,it)
      raa0 = out_RelativeAzimuthAngle(ix,it)
 
!  now use the out_RelativeAzimuthAngle from m_read_input_tio
!  ! +raa has the same effect as -raa due to symmetry,
!  !    and RAA needs to be within [0.,180] for use with LUT
!  this is taken care of in m_read_input_tio

      ! flags for deciding whether to skip calculation 
      pflag00=0 ! for lat/lon
      pflag01=0 ! for sza/vza/raa

      ! bit 0 of out_ProcessingQaulityFlag is for geolocation & geometry
      if((rad_Latitude(ix,it) .lt. -90.) .or. (rad_Latitude(ix,it) .gt. 90.) .or. &
        (rad_Longitude(ix,it) .lt. -180.) .or. (rad_Longitude(ix,it) .gt. 180.)) then
        pflag00=pflag00+1
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),0)
      endif

      if((sza0 .lt. 0.) .or. (sza0 .gt. max_SZA)) then
        pflag01=pflag01+1
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),0)
      endif

      if((vza0 .lt. 0.) .or. (vza0 .gt. max_VZA)) then
        pflag01=pflag01+1
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),0)
      endif

!~~~~~~~~~~
! out_RelativeAzimuthAngle is now taken care of within m_read_input_tio
!      invalid raa is set to a large negative value there
! the following is no longer needed, kept here as a clarification for RAA definition
!      if((rad_SolarAzimuthAngle(ix,it) .ge. -360.) .and. (rad_SolarAzimuthAngle(ix,it) .le. 360.) .and. &
!           (rad_ViewingAzimuthAngle(ix,it) .ge. -360.) .and. (rad_ViewingAzimuthAngle(ix,it) .le. 360.)) then
!      ! RAA = SAA - VAA + PI, Why +PI?
!      !xliu: this is related to how the SAA and VAA are defined
!      !    RAA of forward scattering=0, RAA of backward scattering=180.
!      !also see Eun-Su Yang slide for explanation
!        temp_raa=rad_SolarAzimuthAngle(ix,it)+180.0-rad_ViewingAzimuthAngle(ix,it)
!        ! ensure temp_raa is within [0., 360.) range
!        do while((temp_raa .lt. 0.0) .or. (temp_raa .ge. 360.0))
!          if(temp_raa .ge. 360.0) temp_raa=temp_raa-360.
!          if(temp_raa .lt. 0.0) temp_raa=temp_raa+360.
!        end do
!      else
!        temp_raa= fspecial
!        pflag01 = pflag01 + 1
!        ! temp_raa = -9999. will be skipped in calculation
!      endif
!~~~~~~~~~~~

       if (raa0 .lt. -360.) then
          pflag01=pflag01+1
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),0)
       endif

      ! skip calculation if location/angle are invalid
      if((pflag00 .ge. 1) .or. (pflag01 .ge. 1)) then
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
           go to 990
      endif

      ! get local radiances
      rad466=rad_466nm(ix,it)
      rad440=rad_440nm(ix,it)
      rad477=rad_477nm(ix,it)
      irr466=irr_out_irradiance_466nm(ix)
      irr440=irr_out_irradiance_440nm(ix)
      irr477=irr_out_irradiance_477nm(ix)

      ! ------------------------
      ! calculate cloud fraction
      ! ------------------------
      ! 466nm
      ! bit12 is for skipped ecf
      if ((rad466 .gt. 0.).and.(irr_out_irradiance_466nm(ix) .gt. 0.)) then
      ! earthsunfactor2 = (irr_EarthSunDist/rad_EathSunDist)**2 defined above
      ! calculate the irr expected at time of rad observation
      ! as earthsunfactor is very close to one for TEMPO,this is not important
         xradofirr466=rad466/(irr_out_irradiance_466nm(ix)*earthsunfactor2)
         if (PerturbRadOfIrr466) then
            yradofirr466 = 0.
            do iord = 0, nord_RoI466pert
              !fortran index starts from 1
              iord1 = iord + 1
               yradofirr466=yradofirr466+&
                         RoI466PertCoef(iord1)*(xradofirr466**iord)
            enddo
            ! avoid unphysical values from perturbation
            if (yradofirr466 .gt. 0.) then
               rad_of_irr466(ix,it) = yradofirr466
            else
               rad_of_irr466(ix,it) = xradofirr466
            endif ! yradofirr466
         else ! PerturbRadOfIrr466
            rad_of_irr466(ix,it)=xradofirr466
         endif !PerturbRadOfIrr466
      else
         rad_of_irr466(ix,it) = fspecial
         ! bit8 was set in m_read_input_tio, set again here 
         ! but is also useful later when we implement iteration
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),8)
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12) 
         go to 990
      endif 

      ! 440nm
      ! bit7 was set in m_read_input_tio, set again here
      if ((rad440 .gt. 0.).and.(irr_out_irradiance_440nm(ix) .gt. 0.)) then
         rad_of_irr440(ix,it)=rad440/(irr_out_irradiance_440nm(ix)*earthsunfactor2)
      else
         rad_of_irr440(ix,it) = fspecial
         ! bit7 was set in m_read_input_tio, thus a safeguard below
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
         ! 440 is not used for ECF, do not skip calculation
      endif

      ! 477nm
      if ((rad477 .gt. 0.).and.(irr_out_irradiance_477nm(ix) .gt. 0.)) then
         rad_of_irr477(ix,it)=rad477/(irr_out_irradiance_477nm(ix)*earthsunfactor2)
      else
         rad_of_irr477(ix,it) = fspecial
         ! no bit in processing quality flag for 477nm 
         ! 477 is not used for ECF, do not skip calculation
      endif

      !----------------
      !get psfc0 from model
      !----------------
      !out-of-range rad_Longitude should have been skipped already
      if(name_option_TemperaturePressure.eq.'GMI') then ! for test
        call get_GMIpsfc_lonlat(lon0, lat0, gmi_psfc)

        ! gmi_psfc is NOT corrected for topography
        ! assign l2_TerrainPressure for GMI here
        l2_TerrainPressure(ix,it) = gmi_psfc

        ! psfc0 is used later
        psfc0=gmi_psfc
      endif

      if(name_option_TemperaturePressure.eq.'GEOS5') then ! for TEMPO
        ! l2_TerrainPressure was assigned in read_geoscf
        ! which is called before this subroutine
        ! for this option, it is adjusted for topography
        psfc0=l2_TerrainPressure(ix,it)
      endif

      ! if psfc0[hPa] is out of LUT range, skip calculation
      ! bit 3 is for surface pressure or surface albedo problem
      if ((psfc0 .lt. 0.) .or. (psfc0 .gt. lut_psfc(npsfc))) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
         goto 990
      endif

      !--------------------------
      ! get actual alb0 & alb440
      !--------------------------
      if(name_option_SurfaceReflectivity.eq.'Kleipool') then ! testing
        kleipool466 = fspecial
        kleipool440 = fspecial
        kleipool477 = fspecial
        call get_kleipool_lonlat(lon0, lat0, kleipool466, kleipool440, kleipool477)
        alb0 = kleipool466
        alb440 = kleipool440
      endif

      if(name_option_SurfaceReflectivity.eq.'BRDF') then ! TEMPO
        alb0=BRDF_SurfaceReflectivity466(ix,it)
        alb440=BRDF_SurfaceReflectivity440(ix,it)
      endif

      ! moved out_SurfaceReflectivity assignment from m_cal_ocp here
      ! bound alb0 within [0.,1.] if it is in reasonable range,+/-0.2
      if((alb0 .ge. -0.2) .and. (alb0 .lt. 0.0)) alb0=0.0
      if((alb0 .gt.  1.0) .and. (alb0 .le. 1.2)) alb0=1.0
      out_SurfaceReflectivity466(ix,it)=alb0
 
      ! otherwise set processing quality flag and skip the calculation
      ! note bad alb0 is a large negative number
      if ((alb0 .lt. -0.2) .or. (alb0 .gt. 1.2)) then
          out_SurfaceReflectivity466(ix,it)=fspecial
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
          go to 990
      endif

      ! bound alb440 within [0.,1.] if it is reasonable
      if((alb440 .ge. -0.2) .and. (alb440 .lt. 0.0)) alb440=0.0
      if((alb440 .gt.  1.0) .and. (alb440 .le. 1.2)) alb440=1.0
      if ((alb440 .lt. -0.2) .or. (alb440 .gt. 1.2)) alb440=fspecial
      out_SurfaceReflectivity440(ix,it)=alb440
 
      !-------------------
      ! interpolation prep for alb/sza/vza/raa/psfc
      !-------------------

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

      ! added exit within do loops for alb,sza,vza,raa,psfc
      ! to terminate loop once nodes are found 
      ! if interpolation nodes cannot be found skip calculation

      ! find nodes for alb0
      ialb1=-9; ialb2=-9
      do ialb=1,nalb-1
        if((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
          ialb1=ialb
          ialb2=ialb+1
          walb1=alb0-lut_alb(ialb)
          walb2=lut_alb(ialb+1)-alb0
          exit
        endif
      end do
      if ((ialb1 .lt. 0) .or. (ialb2 .lt. 0) .or. (walb1 .lt. 0.) &
          .or. (walb2 .lt. 0.)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
         go to 990
      endif

      ! find nodes for sza0
      isza1=-9; isza2=-9
      do isza=1,nsza-1
        if((sza0 .ge. lut_sza(isza)) .and. (sza0 .le. lut_sza(isza+1))) then
          isza1=isza
          isza2=isza+1
          wsza1=sza0-lut_sza(isza)
          wsza2=lut_sza(isza+1)-sza0
          exit
        endif
      end do
      if ((isza1 .lt. 0) .or. (isza2 .lt. 0) .or. (wsza1 .lt. 0.) &
          .or. (wsza2 .lt. 0)) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
          go to 990
      endif

      ! find nodes for vza0
      ivza1=-9; ivza2=-9
      do ivza=1,nvza-1
        if((vza0 .ge. lut_vza(ivza)) .and. (vza0 .le. lut_vza(ivza+1))) then
          ivza1=ivza
          ivza2=ivza+1
          wvza1=vza0-lut_vza(ivza)
          wvza2=lut_vza(ivza+1)-vza0
          exit
        endif
      end do
      if ((ivza1 .lt. 0) .or. (ivza2 .lt. 0) .or. (wvza1 .lt. 0.) &
          .or. (wvza2 .lt. 0.)) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
          go to 990
      endif

      ! find nodes for raa0
      iraa1=-9; iraa2=-9
      do iraa=1,nraa-1
        if((raa0 .ge. lut_raa(iraa)) .and. (raa0 .le. lut_raa(iraa+1))) then
          iraa1=iraa
          iraa2=iraa+1
          wraa1=raa0-lut_raa(iraa)
          wraa2=lut_raa(iraa+1)-raa0
          exit
        endif
      end do
      if ((iraa1 .lt. 0) .or. (iraa2 .lt. 0) .or. (wraa1 .lt. 0.) &
           .or. (wraa2 .lt. 0.)) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
          go to 990
      endif

      !----------------------------------
      ! bound psfc by the largest lut_psfc
      ! Note lut_psfc fully covers the expected psfc range
      ! bad psfc0 should have been skipped before
      ! this is only for safeguard 
      if (psfc0 .gt. lut_psfc(npsfc)) psfc0=lut_psfc(npsfc)

      ! find nodes for psfc0
      ipsfc1=-9; ipsfc2=-9
      do ipsfc=1,npsfc-1
        if((psfc0 .ge. lut_psfc(ipsfc)) .and. (psfc0 .le. lut_psfc(ipsfc+1))) then
          ipsfc1=ipsfc
          ipsfc2=ipsfc+1
          wpsfc1=psfc0-lut_psfc(ipsfc)
          wpsfc2=lut_psfc(ipsfc+1)-psfc0
          exit
        endif
      end do
      if ((ipsfc1 .lt. 0) .or. (ipsfc2 .lt. 0) .or. &
          (wpsfc1 .lt. 0.) .or. (wpsfc2 .lt. 0)) then
        write(*,*) "Error *** Surface Pressure too small *** ",ix,it
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
        go to 990
      endif

      !-------------
      ! get LUT values at interpolation node
      ! all node indices should have been found now
      ! may use a function in future
      !-------------

      !------------------------------
      !radiance at surface: clear sky 466nm
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
      ! safeguard
      ! -----------------
      ! entries in lut_rad_clr are always >0.
      ! the following is a safeguard
      if((r11111 .lt. 0.0) .or. (r11112 .lt. 0.0)) then
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
      ! original OMCDO2N uses 466nm radiance at 700 hPa for cloudy sky
      ! in LUT_4660_RAD.h5, ALB(18)=0.8, Psfc(18)=701hPa, 1-indexed
      ! NOTE: this ASSUMES Acloud=0.8 & Pcloud=701hPa for ecf
      ! for cloud albedo rationale, refer to Stammes et al. [2008]
      ! the linear (1) ecf (2) ocp process is also inherited from OMCDO2N
      !
      ! As cal_rad_cld depends on pcld, it makes sense to do ECFOCP iteration
      ! with first pass uses pcld=700hPa as a start
      ! However, TOA rad for overcast is insensitive to pcld 
      ! error associated with the pcld assumption for ecf
      ! is quite small, as cal_rad_cld >> cal_rad_clr
      ! thus, in most cases, iteration is not required
      ! Nonetheless, some geometries still trigger a few iterations
      !--------------------------------
      ialb= LUTrad_cloud_albid ! ALB(18)=0.8 in LUT
 
      ! change from fixed Psfc(18)=701hPa to iteration with ocp
      ipsfc = lut_pcld_indarr(ix,it) ! LUTrad_cloud_psfcid=18
      ! lut_pcld_indarr is the level where lut_pcld just<= ocp
      ! lut_pcld_indarr+1 have lut_pcld>ocp
      ipsnext = ipsfc + 1

      ! on 1st pass, lut_pcld_indarr is initialized to 18
      ! on subsequent passes, some elements may be negative
      if (ipsfc .lt. 0.) then ! failed ocp retrieval
          ! revert back to 701hPa
          ipsfc = LUT_pcld_700hPa 
          ipsnext = LUT_pcld_700hPa
      else if (ipsfc .eq. npcld) then
          ipsnext = npcld
      endif

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

      this_cal_radcld =(wsza2*r1+wsza1*r2)/(wsza1+wsza2)

      if ((ipsfc .eq. ipsnext).or.(ecfocp_iternum .eq. 1)) then
          ! boundary cases end up here
          cal_rad_cld(ix,it) = this_cal_radcld
      else ! vertical interplation
         r111=lut_rad_clr(ialb,isza1,ivza1,iraa1,ipsnext)
         r112=lut_rad_clr(ialb,isza1,ivza1,iraa2,ipsnext)
         r121=lut_rad_clr(ialb,isza1,ivza2,iraa1,ipsnext)
         r122=lut_rad_clr(ialb,isza1,ivza2,iraa2,ipsnext)
         r211=lut_rad_clr(ialb,isza2,ivza1,iraa1,ipsnext)
         r212=lut_rad_clr(ialb,isza2,ivza1,iraa2,ipsnext)
         r221=lut_rad_clr(ialb,isza2,ivza2,iraa1,ipsnext)
         r222=lut_rad_clr(ialb,isza2,ivza2,iraa2,ipsnext)
         r11=(wraa2*r111+wraa1*r112)/(wraa1+wraa2)
         r12=(wraa2*r121+wraa1*r122)/(wraa1+wraa2)
         r21=(wraa2*r211+wraa1*r212)/(wraa1+wraa2)
         r22=(wraa2*r221+wraa1*r222)/(wraa1+wraa2)
         r1=(wvza2*r11+wvza1*r12)/(wvza1+wvza2)
         r2=(wvza2*r21+wvza1*r22)/(wvza1+wvza2)

         that_cal_radcld =(wsza2*r1+wsza1*r2)/(wsza1+wsza2)

         thisocp = out_CloudPressure(ix,it) ! from previous ecfocp iteration 
         ! thisocp should be valid and between the 2 levels
         ! otherwise, ipsc=ipsnext would have been taken care of before
         ! however, as a safeguard, add condition just to make sure
         thiswp = thisocp - lut_pcld(ipsfc)
         thatwp = lut_pcld(ipsnext) - thisocp
         if ((thiswp .ge. 0.) .and. (thatwp .gt. 0.)) then 
            cal_rad_cld(ix,it)=(this_cal_radcld*thatwp+&
                   that_cal_radcld*thiswp)/(thiswp + thatwp)
         else
            cal_rad_cld(ix,it) = this_cal_radcld
         endif
      endif

      !--------------------------------
      ! 440nm radiance at 700 hPa: cloudy sky
      ! in LUT_4400_RAD.h5, ALB(18)=0.8, Psfc(18)=701hPa
      ! pcld=701hPa is an used in 1st pass of ecf calculation
      !--------------------------------
      ! use the same ialb & ipsfc as 466 above 
      ! 440nm does not affect ECF
      ! cal_rad_cld440 is used for crf440

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

      !--------------------------------
      ! 477nm radiance at 700 hPa: cloudy sky
      ! in LUT_4770_CLEAR.h5, ALB(18)=0.8, Psfc(18)=701hPa
      ! pcld=701hPa is an used in 1st pass of ecf calculation
      !--------------------------------
      ! use the same ialb & ipsfc as 466 above 
      ! 477nm does not affect ECF
      ! cal_rad_cld477 is used for crf477

      r111=lut_rad_ler(ialb,isza1,ivza1,iraa1,ipsfc)
      r112=lut_rad_ler(ialb,isza1,ivza1,iraa2,ipsfc)
      r121=lut_rad_ler(ialb,isza1,ivza2,iraa1,ipsfc)
      r122=lut_rad_ler(ialb,isza1,ivza2,iraa2,ipsfc)
      r211=lut_rad_ler(ialb,isza2,ivza1,iraa1,ipsfc)
      r212=lut_rad_ler(ialb,isza2,ivza1,iraa2,ipsfc)
      r221=lut_rad_ler(ialb,isza2,ivza2,iraa1,ipsfc)
      r222=lut_rad_ler(ialb,isza2,ivza2,iraa2,ipsfc)
      r11=(wraa2*r111+wraa1*r112)/(wraa1+wraa2)
      r12=(wraa2*r121+wraa1*r122)/(wraa1+wraa2)
      r21=(wraa2*r211+wraa1*r212)/(wraa1+wraa2)
      r22=(wraa2*r221+wraa1*r222)/(wraa1+wraa2)
      r1=(wvza2*r11+wvza1*r12)/(wvza1+wvza2)
      r2=(wvza2*r21+wvza1*r22)/(wvza1+wvza2)

      cal_rad_cld477(ix,it)=(wsza2*r1+wsza1*r2)/(wsza1+wsza2)

      !-----------------------------------
      !calculate effective cloud fraction ecf and cloud radiance fraction crf
      !-----------------------------------
      ! added conditional safeguard
      if ((cal_rad_clr(ix,it) .gt. 0.) .and. (cal_rad_cld(ix,it) .gt. 0.)) then
         rout_ecf=(rad_of_irr466(ix,it)-cal_rad_clr(ix,it))/(cal_rad_cld(ix,it)-cal_rad_clr(ix,it))
      else
         rout_ecf = fspecial
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
      endif

      ! perturb ecf if requested 
      if ((PerturbECF) .and. (rout_ecf .ge. 0.)) then 
          rout_ecf = ECFPertCoef(1) + ECFPertCoef(2) * rout_ecf
      endif

      ! calculate cloud radiance fraction at 466
      ! crf definition follows Vasilkov et al.
      if ((cal_rad_cld(ix,it) .gt. 0.).and.(rad_of_irr466(ix,it).gt. 0.)) then
         rout_crf466=rout_ecf*cal_rad_cld(ix,it)/rad_of_irr466(ix,it)
      else
         rout_crf466 = fspecial
      endif

      ! caculate cloud radiance fraction at 440
      if ((cal_rad_cld440(ix,it) .gt. 0.) .and. (rad_of_irr440(ix,it).gt. 0.)) then
         rout_crf440=rout_ecf*cal_rad_cld440(ix,it)/rad_of_irr440(ix,it)
      else
        rout_crf440 = fspecial
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
      endif

      ! calculate cloud raiance fraction at 477
      if ((cal_rad_cld477(ix,it) .gt. 0.) .and. (rad_of_irr477(ix,it).gt. 0.)) then
         rout_crf477=rout_ecf*cal_rad_cld477(ix,it)/rad_of_irr477(ix,it)
      else
        rout_crf440 = fspecial
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
      endif

      ! assign non-clipped ecf & crf to array
      ! due to assumption of Acld=0.8 (some pcld=701hPa) and IPA
      ! ecf and crf can be negative or above 1.0 within a reasonable range 
      out_EffectiveCloudFractionNotClipped(ix,it)= rout_ecf
      out_CloudRadianceFractionNotClipped440(ix,it)= rout_crf440
      out_CloudRadianceFractionNotClipped466(ix,it)= rout_crf466
      out_CloudRadianceFractionNotClipped477(ix,it)= rout_crf477

      !------------------------------------------
      ! clip ecf & crf to [0.0, 1.0]
      ! ecf_lowclip & ecf_highclip are used for both ECF and CRF
      ! added logic to differentiate skipped or bad calculation
      ! out_ProcessingQualityFlag bit 9 for out-of-range,clipped (WARNING)
      !     these are set to 0.0 or 1.0, and may still be usable    
      ! out_ProcessingQualityFlag bit12 for unreasonable values (ERROR)
      !     these should not be used

      if ((rout_ecf .lt. 0.) .and. (rout_ecf .ge. ecf_lowclip)) then
         rout_ecf=0.
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),9)
      endif
      if ((rout_ecf .gt. 1.) .and. (rout_ecf .le. ecf_highclip)) then
         rout_ecf=1.
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),9)
      endif
      if((rout_ecf .lt. ecf_lowclip) .or. (rout_ecf .gt. ecf_highclip)) then
         rout_ecf=fspecial 
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),12)
      endif

      if ((rout_crf466 .lt. 0.).and.(rout_crf466 .ge. ecf_lowclip)) rout_crf466=0.
      if ((rout_crf466 .gt. 1.).and.(rout_crf466 .le. ecf_highclip)) rout_crf466=1.
      if ((rout_crf466 .lt. ecf_lowclip) .or. (rout_crf466 .gt. ecf_highclip)) then
        rout_crf466=fspecial
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),1)
      endif

      if ((rout_crf440 .lt. 0.).and.(rout_crf440 .ge. ecf_lowclip)) rout_crf440=0.
      if ((rout_crf440 .gt. 1.).and.(rout_crf440 .le. ecf_highclip)) rout_crf440=1.
      if ((rout_crf440 .lt. ecf_lowclip) .or. (rout_crf440 .gt. ecf_highclip)) then
        rout_crf440=fspecial
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),7)
      endif

      if ((rout_crf477 .lt. 0.).and.(rout_crf477 .ge. ecf_lowclip)) rout_crf477=0.
      if ((rout_crf477 .gt. 1.).and.(rout_crf477 .le. ecf_highclip)) rout_crf477=1.
      if ((rout_crf477 .lt. ecf_lowclip) .or. (rout_crf477 .gt. ecf_highclip)) then
        rout_crf477=fspecial
        ! no bit in processing quality flag for crf477
      endif

      ! assign clipped ecf & crf to array
      out_EffectiveCloudFraction(ix,it)=rout_ecf
      out_CloudRadianceFraction440(ix,it)=rout_crf440
      out_CloudRadianceFraction466(ix,it)=rout_crf466
      out_CloudRadianceFraction477(ix,it)=rout_crf477

      ! out_ReflectanceFactor is equivalent to observed Lambertian reflectance at 466
      if(rad_of_irr466(ix,it) .gt. 0.0) then
         out_ReflectanceFactor(ix,it)=pi*rad_of_irr466(ix,it)/cos(dtor*sza0)
      else
         out_ReflectanceFactor(ix,it)=fspecial
      endif

! skip to here when something goes wrong during processing
990   continue

      ! debug
      if ((trim(run_mode).eq.'development').and.(it.eq.itdebug)) then
         write(lun_debug_ecf,*)ix,alb0,psfc0,rad_of_irr466(ix,it),&
         cal_rad_clr(ix,it),cal_rad_cld(ix,it),&
         out_EffectiveCloudFractionNotClipped(ix,it)
      endif

      ! safeguard ecf 
      if (btest(out_ProcessingQualityFlags(ix,it),12)) then
         out_EffectiveCloudFraction(ix,it) = fspecial
         out_CloudRadianceFraction440(ix,it) = fspecial
         out_CloudRadianceFraction466(ix,it) = fspecial
         out_EffectiveCloudFractionNotClipped(ix,it) = fspecial
         out_CloudRadianceFractionNotClipped440(ix,it) = fspecial
         out_CloudRadianceFractionNotClipped466(ix,it) = fspecial
         out_CloudRadianceFractionNotClipped477(ix,it) = fspecial
      endif   

 3455  continue ! skip to here when ecf does not need re-calculation

      !=====
    end do ! ix
  end do ! it
  !=====

  !close debug file unit
  close(lun_debug_ecf)

!**********************
end subroutine cal_ecf
!**********************

!22222222222222222222222
!**********************
subroutine allocate_ecf_arrays(nx,nt,fspecial,ierr)
!**********************
!222222222222222222222222
   use m_vars, only: cal_rad_clr, cal_rad_cld, cal_rad_cld440, cal_rad_cld477,&
       rad_of_irr440, rad_of_irr466, rad_of_irr477, out_ReflectanceFactor,&
       out_SurfaceReflectivity466,out_SurfaceReflectivity440,&
       out_EffectiveCloudFraction,out_EffectiveCloudFractionNotClipped,&
       out_CloudRadianceFraction440, out_CloudRadianceFractionNotClipped440,&
       out_CloudRadianceFraction466, out_CloudRadianceFractionNotClipped466,&
       out_CloudRadianceFraction477, out_CloudRadianceFractionNotClipped477,&
       lut_pcld_indarr, LUT_pcld_700hPa, ecf_niter,ocp_niter
   use m_vars, only: prev_ecf, prev_ecf_notclipped

   implicit none

   real, intent(in):: fspecial
   integer, intent(in):: nx,nt
   integer, intent(inout):: ierr

  allocate(cal_rad_clr(nx,nt),stat=ierr)
  allocate(cal_rad_cld(nx,nt),stat=ierr)
  allocate(cal_rad_cld440(nx,nt),stat=ierr)
  allocate(cal_rad_cld477(nx,nt),stat=ierr)

  allocate(rad_of_irr440(nx,nt),stat=ierr)
  allocate(rad_of_irr466(nx,nt),stat=ierr)
  allocate(rad_of_irr477(nx,nt),stat=ierr)

  allocate(out_SurfaceReflectivity466(nx,nt),stat=ierr)
  allocate(out_SurfaceReflectivity440(nx,nt),stat=ierr)

  cal_rad_clr=fspecial
  cal_rad_cld=fspecial
  cal_rad_cld440=fspecial
  cal_rad_cld477=fspecial
  rad_of_irr440=fspecial
  rad_of_irr466=fspecial
  rad_of_irr477=fspecial

  out_SurfaceReflectivity466=fspecial
  out_SurfaceReflectivity440=fspecial

  allocate(out_EffectiveCloudFraction(nx,nt),stat=ierr)
  allocate(out_EffectiveCloudFractionNotClipped(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFraction440(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFractionNotClipped440(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFraction466(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFractionNotClipped466(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFraction477(nx,nt),stat=ierr)
  allocate(out_CloudRadianceFractionNotClipped477(nx,nt),stat=ierr)
  allocate(out_ReflectanceFactor(nx,nt),stat=ierr)

  out_EffectiveCloudFraction=fspecial
  out_EffectiveCloudFractionNotClipped=fspecial
  out_CloudRadianceFraction440=fspecial
  out_CloudRadianceFractionNotClipped440=fspecial
  out_CloudRadianceFraction466=fspecial
  out_CloudRadianceFractionNotClipped466=fspecial
  out_CloudRadianceFraction477=fspecial
  out_CloudRadianceFractionNotClipped477=fspecial
  out_ReflectanceFactor=fspecial

  allocate(prev_ecf(nx,nt),stat=ierr)
  allocate(prev_ecf_notclipped(nx,nt),stat=ierr)
  prev_ecf=fspecial 
  prev_ecf_notclipped=fspecial

  allocate(lut_pcld_indarr(nx,nt),stat=ierr)
  ! initialize to 700hPa which is index 18 in lut_pcld
  ! this will be used on 1st pass through ecfocp
  lut_pcld_indarr= LUT_pcld_700hPa

  ! ecfocp number of iteration 
  allocate(ecf_niter(nx,nt),stat=ierr)
  ecf_niter = 0 
  allocate(ocp_niter(nx,nt),stat=ierr)
  ocp_niter = 0

!**********************
end subroutine allocate_ecf_arrays
!**********************

end module m_cal_ecf
