module m_cal_ocp
  public cal_ocp

contains
!1111111111111111111
!******************
subroutine cal_ocp(ecfocp_iternum)
!******************
!11111111111111111111
  ! -----------------------------
  ! define ProcessingQualityFlags
  ! -----------------------------
  ! bit??  meaning:WhereSet
  !------------------------------
  ! bit00  (Error) invalid lat/lon/SZA/VZA/RAA: m_cal_ecf.f90
  ! bit01  (Warning) invalid 466nm cloud radiance fraction (crf): m_cal_ecf.f90
  ! bit02  (Warning) ecf<minECF pcld replaced by pscene: m_cal_pscene.f90 
  ! bit03  (ERROR) input surface pressure or albedo error: m_cal_ecf.f90
  ! bit04  (Warning) snow_ice_fraction>min_snowice pcld replaced by pscene: m_cal_pscene.f90
  ! bit05  (Warning) SCD correction problem or max_scd_iter reached in ocp : m_cal_ocp.f90
  ! bit06  (Error) SCD < 0 or SCD quality issue, : m_cal_ocp.f90
  ! bit07  (Warning) invalid 440nm irr,rad or crf: m_read_input_tio.f90,m_cal_ecf.f90
  ! bit08  (ERROR) 466nm radiance or irradiance error: m_read_input_tio.f90
  ! bit09  (Warning) ecf out of normal range, clipped: m_cal_ecf.g90
  ! bit10  (Info) derived SurfaceLER or TerrainPressure error in pscene
  ! bit11  (Info) derived ScenePressure or SceneLER error in pscene
  ! bit12  (ERROR) skipped cloud ecf calculation 
  !        due to any problem during processing: m_cal_ecf.f90
  ! bit13  (ERROR) skipped cloud ocp calculation due to
  !        any problem during processing,or invalid ocp: m_cal_ocp.f90
  ! bit14  (Warning) ocp out of normal range and clipped: m_cal_ocp.f90
  ! bit15  (Info) skipped pscene calculation during processing
  ! bit 9 & bit 14 have changed from Error to Warning

  use m_vars
  use m_read_GMI
  use m_read_hdf5
  use m_scd_adjust

  implicit none

  ! input variable
  integer, intent(in):: ecfocp_iternum

  !local variable moved from m_vars
  real:: alb0, sza0, vza0, raa0, psfc0

  !local variables
  real:: cal_ecf,cal_crf
  integer::ialb, isza, ivza, iraa, ipsfc, ipcld
  integer::ialb1,isza1,ivza1,iraa1,ipsfc1
  integer::ialb2,isza2,ivza2,iraa2,ipsfc2
  real::   walb1,wsza1,wvza1,wraa1,wpsfc1,vpsfc1,apsfc1
  real::   walb2,wsza2,wvza2,wraa2,wpsfc2,vpsfc2,apsfc2
  real::   alb440

  real::yy1,yy2,ww1,ww2
  real(kind=4)::cpp, t8p
  real :: temp_cpp, temp_t8p, delta_temp
  integer::ipsfc0,ipm0,ipm1,ipm2
  integer(kind=4)::iflag

  integer(kind=4)::ierr
  integer(kind=4)::nt,nx
  integer(kind=4)::it,ix, iternum

  integer(kind=4)::ip

  integer :: index_pcld_lut

  real:: gmi_psfc
  real(kind=4), dimension(:), allocatable :: tt, pp, qq, ppdry

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

  real::pi,dtor

  ! add local variables
  real:: lat0, lon0
  real:: thisocp, thatocp, ocp_change
  real:: fFillValue9

  ! ------
  ! local initialization
  ! ------
  pi=4.*atan(1.)
  dtor=pi/180.

  fFillValue9 = -9999.

  nt=rad_NumTimes
  nx=rad_nXtrack

  maxpress = 1200 !hPa

  ! allocate and initialize local array
  allocate(tt(nlayers),stat=ierr)
  allocate(qq(nlayers),stat=ierr)
  allocate(pp(nlayers+1),stat=ierr)
  allocate(ppdry(nlayers+1),stat=ierr)
  tt = fFillValue9
  pp = fFillValue9
  ppdry = fFillValue9
  qq = 0.

! debug
  if ((trim(run_mode) .eq. 'development').and.(ecfocp_iternum .eq. 1) & 
     .and. (ixdebug .ge. 0) .and. (itdebug .ge. 0)) then
  write(*,*) 'writing debug_scd_adjust.txt'
  open(unit=lun_debug_scdadj,file='debug_scd_adjust.txt')
  write(lun_debug_scdadj,*) &
    '    ix   scdmorg    scdm     scdadj      temp_t8p     t8p    temp_cpp'
  endif

  ! ==========
  do it=1,nt
    do ix=1,nx
      ! ==========
      ! if ocp_change is small enough, no need to recalculate 
      ! previous iteration values and quality flags are valid here
      ! simply skip to next ground pixels

      if (ecfocp_iternum .gt. 2) then ! check only after 2 passes
         thisocp = out_CloudPressureNotClipped(ix,it)
         thatocp = prev_ocp_notclipped(ix,it)
         ocp_change = abs(thisocp - thatocp)
         if ((thisocp .gt. 0.) .and. (thatocp .gt. 0.).and. &
             (ocp_change .lt. delta_ocp)) then
             ! no ProcessingQualityFlags change here
             ! keep previous flags, skip calculation
             go to 3456
         else
             ! assign current ocp to previous ocp for next ecfocp iteration
             ! current ocp will be re-calculated 
             ! and compared with previous ocp next time around
             prev_ocp_notclipped(ix,it) = out_CloudPressureNotClipped(ix,it)
             prev_ocp(ix,it) = out_CloudPressure(ix,it)
             ! update ocp_niter, skipped calculation will not end up here
             ocp_niter(ix,it) = ecfocp_iternum
         endif
      else ! for 1st couple of iterations
         ! update ocp_niter
         ocp_niter(ix,it) = ecfocp_iternum       
      endif

      ! -----------
      ! clear relevant out_ProcessingQualityFlags bits
      ! Note, not all bits used here need clearing
      ! as some(e.g. bit0, bit3)  were cleared in cal_ecf called before
      ! and are used here
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),5)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),6)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),13)
      out_ProcessingQualityFlags(ix,it)=ibclr(out_ProcessingQualityFlags(ix,it),14)

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

! geolocation and angles were checked in m_cal_ecf
! invalid values triggers bit0 to be set
! thus, check bit0 to decide whether to skip calculation
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

      ! skip ocp if cal_ecf or cal_crf are bad or ZERO
      ! when ecf//crf=0, there is no need to calculate ocp, remain fill value
      ! bit13 indicates ocp is skipped
      if ((cal_ecf .le. 0.) .or. (cal_crf .le. 0.) .or. (cal_ecf .gt. 1.)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)  
         go to 990
      endif

      ! skip if bit3 (psfc//rsfc error) is set
      if (btest(out_ProcessingQualityFlags(ix,it),3)) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
      endif

      !----------------------------------------------
      ! T-P profile and vvcd
      !----------------------------------------------
      ! initialize local array
      tt(:) = fFillValue9
      pp(:) = fFillValue9
      ppdry(:) = fFillValue9
      qq(:) = 0.
      vvcd(:) = fFillValue9

      ! ----------------------------------------------
      ! option for TemperaturePressure & SurfacePressure
      ! ----------------------------------------------
      ! bad psfc should have been skipped before
      ! OMI option as a backup
      if((name_option_TemperaturePressure.eq.'GMI')) then 
        call get_GMItmp_lonlat(lon0,lat0,tt,pp,gmi_psfc)

        psfc0=gmi_psfc
        ! adds safeguard
        if ((psfc0 .lt. lut_psfc(1)).or.(psfc0.gt.lut_psfc(npsfc))) then 
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
          go to 990
        endif
        if (nlayers .NE. gmi_np) then
           write(*,*)'gmi_np,nlayers incompatible',gmi_np,nlayers
           call exit(-1)
        endif

        ppdry = pp ! GMI is not adjusted to ppdry
 
        call read_GMI_VCD(pp,tt)
        vvcd=gmi_vcd
      endif ! GMI

      !---TEMPO option for GEOS-CF
      if(name_option_TemperaturePressure.eq.'GEOS5') then
        psfc0=l2_TerrainPressure(ix,it)
        ! safeguard
        if ((psfc0 .lt. lut_psfc(1)).or.(psfc0 .gt. lut_psfc(npsfc))) then
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),3)
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
          go to 990
        endif
        if (nlayers .NE. geos_np) then
           write(*,*)'ERROR: geos_np nlayers mismatch',geos_np,nlayers
           call exit(-1)
        endif

        ! geos_Pressure,temperature,Q is from TOA to BOA
        do ip=1,geos_np
          pp(ip)=geos_Pressure(ix,it,ip)
          tt(ip)=geos_Temperature(ix,it,ip)
          qq(ip)=geos_Q(ix,it,ip)
        end do
        pp(geos_np+1) = psfc0

       ! pp, tt, qq are on GEOS-CF grid
       ! vvcd = geos_vcd is on LUT pcld grid
        call read_GEOS5_VCD(pp,tt,qq,ppdry)
        vvcd=geos_vcd
      endif ! GEOS5

! debug
      if ((trim(run_mode).eq.'development').and.(ecfocp_iternum .eq. 1) &
          .and.(it .eq. itdebug).and. (ix .eq. ixdebug)) then
         write(*,*) '  writing debug_tpocp.txt'
         open(unit=lun_debug_ocp,file='debug_tpocp.txt')
         write(lun_debug_ocp,*) 'name_option_TemperaturePressure=', &
               trim(name_option_TemperaturePressure)
         write(lun_debug_ocp,*)'ix, it=',ix,it
         write(lun_debug_ocp,*)'latitude=',rad_latitude(ix,it)
         write(lun_debug_ocp,*)'longitude=',rad_longitude(ix,it)
         write(lun_debug_ocp,*) 'psfc=',pp(nlayers+1)
         write(lun_debug_ocp,*) 'nlayers=',nlayers
         write(lun_debug_ocp,*) 'Level, Pressure(hPa), DryPressure,  Temperature(K)'
         do ip = 1, nlayers
            write(lun_debug_ocp,*)ip, pp(ip),ppdry(ip),tt(ip)
         end do

         write(lun_debug_ocp,*) 'LUT_Level, Pressure (hPa), vvcd'
         do ip = 1, npcld
            write(lun_debug_ocp,*) ip, lut_pcld(ip), vvcd(ip)
         enddo
         close(lun_debug_ocp)
      endif

     !---------------------------------
     ! surface reflectivity
     !---------------------------------
     ! directly use out_SurfaceReflectivity assigned in ecf
     ! instead of repeating calculation

      ! skip if bit3 (psfc or rsfc error) is set
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

      ! in read_cldo4_tio, negative or bad SCD are set to fspecial
      if (nasa_SlantColumnAmountO2O2(ix,it) .lt. 0.0) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),6)
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
        go to 990
      endif

      ! whether to skip ocp when 0<ecf<min_ecf
      if (btest(out_ProcessingQualityFlags(ix,it),2)) then ! ecf< min_ecf
        if (name_option_skipECFminocp .EQ. 1) then ! requested skip ocp for ecf<min_ecf
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
        endif
      endif

      ! whether to skip snow/ice for ocp
      ! snowice is handled in pscene, flag was set in m_read_gler
      if (btest(out_ProcessingQualityFlags(ix,it),4)) then ! snowice surface
        if (name_option_skipSnowocp .eq. 1) then ! requested skip ocp for snowice
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
         go to 990
        endif
      endif

      ! skip ocp if ecf is skipped 
      if (btest(out_ProcessingQualityFlags(ix,it),12)) then
            out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13) 
            go to 990 
      endif

      ! -----------------
      ! set nodes for LUT
      ! -----------------
      ! now exit loop as soon as node is found
      ! skip calculation if not found
      ialb1=-9; ialb2=-9
      walb1=0.; walb2=0.
      do ialb=1,nalb-1
        if ((alb0 .ge. lut_alb(ialb)) .and. (alb0 .le. lut_alb(ialb+1))) then
          ialb1=ialb
          ialb2=ialb+1
          walb1=alb0-lut_alb(ialb)
          walb2=lut_alb(ialb+1)-alb0
          exit
        endif
      end do
      if ((ialb1 .lt. 0) .or. (walb1 .lt. 0.) .or. (walb2 .lt. 0.)) then
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
          exit
        endif
      end do
      if ((isza1 .lt. 0) .or. (wsza1 .lt. 0.) .or. (wsza2 .lt. 0.)) then
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
          exit
        endif
      end do
      if ((ivza1 .lt. 0) .or. (wvza1 .lt. 0.) .or. (wvza2 .lt. 0.)) then
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
          exit
        endif
      end do
      if ((iraa1 .lt. 0) .or. (wraa1 .lt. 0.) .or. (wraa2 .lt. 0.)) then
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
      ! in LUT when ipcld>ipsfc, lut_amf_cld<0.
      ! but this would not happen below, as ipsfc=ipcld 
      do ipcld=1,npcld
        ipsfc=ipcld ! as if cloud at surface
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
      ! lut_psfc currently covers until 1100 hPa
      ! reasoable psfc0 should not trigger this warning
      if (psfc0 .gt. lut_psfc(npsfc)) then
         write(*,*) '*** WARNING: psfc0 in cal_ocp is too large ***'
         write(*,*) 'ix,it,psfc0=',ix,it,psfc0
         write(*,*) '*** limit psfc0 to max(lut_pcld) ***'
         psfc0 = lut_psfc(npsfc)
      endif

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
          exit
        endif
      end do
     
      if (ipsfc1 .lt. 0) then ! this should not happen, safeguard
        write(*,*) " *** Pcld: Check Surface Pressure *** ",ix,it,psfc0
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
        go to 990
      else
        vpsfc0=(wpsfc1*vpsfc2+wpsfc2*vpsfc1)/(wpsfc1+wpsfc2)
        apsfc0=(wpsfc1*apsfc2+wpsfc2*apsfc1)/(wpsfc1+wpsfc2)
      endif

      ipsfc0=-9
      do ipsfc=1,npsfc-1
        if ((psfc0 .ge. lut_psfc(ipsfc)) .and. (psfc0 .lt. lut_psfc(ipsfc+1))) then
          ipsfc0=ipsfc
          exit
        endif
      end do
      ! if found: ipsfc0>0 is the lut layer where psfc0 resides
      ! if not found: ipsfc0 remains negative

      ! ipms are used for extrapolation at high pressure end
      ipm0=-9; ipm1=-9; ipm2=-9
      if (ipsfc0 .gt. 0) then ! psfc0 < lut_psfc(npsfc)
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
      ! skip calculation if scdmorg is invalid
      if (scdmorg .lt. 0.) then
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
        go to 990
      endif 

      ! initial iteration use TrefO4 reference
      iternum = 0
      t8p = TrefO4 ! initial reference temperature for SCD retrieval
      temp_t8p = TrefO4
      scdm = scdmorg
      scdadj = scdmorg

      ! amfvcd_int uses psfc0, amfvcd_ext uses Pcld
      ! amfvcd_ext > amfvcd_int when pcld > psfc0
      do ipcld=1,npcld
        amfvcd_int(ipcld)=vvcd(ipcld)*cal_crf*real(cal_amf_cld(ipcld),kind=4) &
             +vpsfc0*(1.0-cal_crf)*apsfc0
        amfvcd_ext(ipcld)=vvcd(ipcld)*cal_crf*real(cal_amf_cld(ipcld), kind=4)&
             +vvcd(ipcld)*(1.0-cal_crf)*real(cal_amf_clr(ipcld), kind=4)
      end do

      !????????????????????????????????????????
      ! move option_psfc_clear to m_vars.f90
      ! where option_psfc_clear=0 is specified
      !       looks like this hardcodes Pclr=Psfc
      ! find pressure for AMF*VCD
      !    0: Pclr = Psfc (fixed) & Pcld > Psfc
      !    1: Pclr = Pcld if Pcld > Psfc
      ! ???????????????????????????????????????

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SCD iteration comes back here to 777
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
777   continue ! iteration

      iflag=-1 

      ! lut_pcld is in increasing order (hpa)
      ! pcld(1)=55, pcld(3) = 76, pcld(5)=104
      if (scdm .le. amfvcd_int(5)) then ! low pressure end
        iflag=0
      endif

      yy1=-9. ; yy2=-9.
      ww1=0. ; ww2=0.
      do ipcld=1,npcld-1
        if ((scdm.gt.amfvcd_int(ipcld)).and.(scdm.le.amfvcd_int(ipcld+1))) then
          iflag=1 ! node found
          yy1=lut_pcld(ipcld)
          yy2=lut_pcld(ipcld+1)
          ww1=scdm-amfvcd_int(ipcld)
          ww2=amfvcd_int(ipcld+1)-scdm
        endif
      end do
      ! for scdm > amfvcd_int(1), iflag change to 1
      ! for scdm <= amfvcd_int(1), iflag remains 0

      if (iflag .eq. 1) then ! normal interpolation
        cpp=(ww1*yy2+ww2*yy1)/(ww1+ww2)

      ! the choice below is for low pressure end
      else if(iflag .eq. 0) then ! scdm<= amfvcd_int(1)
        x0=0.0
        x1=lut_pcld(5)
        x2=lut_pcld(6)
        y0=0.0
        y1=amfvcd_int(5)
        y2=amfvcd_int(6)

        xx=x1 ! the target pressure
        !as 1st & 3rd term below is always 0, yy=y1,it is verbose, however,
        ! the formula is consistent with the one inside ipp loop below
        yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
             +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
             +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
        diff_save=abs(scdm-yy)

        ! increase ipp (dec xx) 1 hPa at a time until min diff found
        ! as 150 hPa > lut_pcld[5:6], it is safe  to use 150
        ! may need change if LUT is changed
        do ipp=1,150
          xx=x1-real(ipp)
          yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
               +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
               +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
          diff=abs(scdm-yy)
          if (diff .ge. diff_save) then
            go to 970 !inflection point found
          else
            diff_save=diff
            if (xx .le. 0.)then
              xx=-9999.
              go to 970 ! exit when no solution found
            endif
          endif
        end do
970     continue
        cpp=real(nint(xx))

      else ! large scdm case: scdm>amfvcd_int(npcld) 
        iflag=2 
        ! the program is expected to bypass this safeguard
        ! because all ipms should be valid
        if ((ipm2 .lt. 1) .or. (ipm1 .lt. 1) .or. (ipm0 .lt. 1)) then
           write(*,*) 'ipm <= 0 for large scdm, this should not happen.'
           cpp = fFillValue9
           out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
           go to 990
        endif

        x0=lut_pcld(ipm0)
        x1=lut_pcld(ipm1)
        x2=lut_pcld(ipm2)

        ! defult in m_vars option_psfc_clear=0 
        if(option_psfc_clear.eq.0) then !original OMCDO2N choice 
        ! use pclr=psfc, ipm0=ipsfc0
          y0=amfvcd_int(ipm0)
          y1=amfvcd_int(ipm1)
          y2=amfvcd_int(ipm2)
        endif

        if(option_psfc_clear.eq.1) then
        ! use pclr=pcld, ipm0=npsfc
          y0=amfvcd_ext(ipm0)
          y1=amfvcd_ext(ipm1)
          y2=amfvcd_ext(ipm2)
        endif

        xx=psfc0
        yy=(xx-x1)*(xx-x2)/(x0-x1)/(x0-x2)*y0 &
             +(xx-x0)*(xx-x2)/(x1-x0)/(x1-x2)*y1 &
             +(xx-x0)*(xx-x1)/(x2-x0)/(x2-x1)*y2
        diff_save=abs(scdm-yy)

        !increase xx 1 hPa at a time from psfc0
        ! until mininal difference is found
        ! 5000 is large and safe, 2000 should be enough
        ! but should make no difference to computer
        do ipp=1,2000
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

      ! skip if cpp <0.
      if (cpp .lt. 0.) then 
          cpp = fFillValue9
          out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
          go to 990
      endif

      ! add clip cpp to within LUT range, safeguard
      if (cpp .gt. lut_psfc(npsfc)) then
          cpp = lut_psfc(npsfc)
      endif

      ! adjust scd according to T at temp_cpp
      ! use the temperature at temp_cpp when in range
      ! initially used T at 0.5*cpp as it is in the middle of pressure
      ! temp_cpp = real (cpp * 0.5, kind=4)
      ! now changed to 0.7937*cpp as it is in the middle of o2o2 column
      ! set frac4cpp = 0.7937 in m_vars
      temp_cpp = real (cpp * frac4cpp, kind=4)
      if ((temp_cpp .ge. pp(1)) .and. (temp_cpp .le. psfc0)) then
        if (name_option_TemperaturePressure .eq. 'GMI') then
          call scd_adjust_gmi(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
        else if (name_option_TemperaturePressure .eq. 'GEOS5') then
      ! scdmorg<0. should have already been skipped
      ! returned scdadj always > 0., because if negative
      ! scdadj = scdmorg and temp_t8p=TrefO4
          call scd_adjust_geos(pp,tt,temp_cpp,scdmorg,scdadj,temp_t8p)
        else
          temp_t8p = real(t8p, kind=4)
        endif
      else
         temp_t8p = real(t8p, kind=4) !this will terminate further iteration 
         ! signal iteration problem
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),5)
      endif

      ! increment iternum
      iternum = iternum + 1

      ! debug
      if ((trim(run_mode).eq.'development').and. &
          (it .eq. itdebug) .and. (ix .eq. ixdebug)) then
         write(*,*) iternum, scdm, scdadj, temp_cpp, temp_t8p
      endif

      ! test if terminate temperature iteration
      !technically, temp_t8p can be -999. when scdmorg<0.
      !but it won't happen as they should have been skipped
      !if negative scdadj occurs within scd_adjust_geos,
      ! scdadj=scdmorg and temp_t8p=TrefO4
      delta_temp = real(abs(t8p - temp_t8p), kind=4)
      if ((delta_temp .lt. dt_threshold).or.(temp_t8p .eq. TrefO4)) then
         goto 990 ! exit iteration
      endif

      if (iternum .lt. max_scd_iter) then 
         t8p = temp_t8p  !update t8p from previous step
         scdm = scdadj !update scdm from previous step
         goto 777  ! goto iteration start
      endif

! skipped calculation will end up here
990   continue

      ! set out_ProcessingQualityFlags bit 5 for max_scd_iter
      if (iternum .ge. max_scd_iter) then
         out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),5)
      endif

      ! debug
      if ((trim(run_mode).eq.'development').and.(it .eq. itdebug) &
          .and.(ecfocp_iternum .eq. 1)) then
         write(lun_debug_scdadj,*) ix, scdmorg, scdm, scdadj, temp_t8p, t8p, temp_cpp
      endif

      ! scdm & t8p is the step right before final iteration
      ! scdadj & temp_t8p is the step right after final iterateion
      if (scdm .gt. 0.) then
         out_SlantColumnAmountO2O2(ix,it) = scdadj !scdm // scdadj
         out_O2O2CloudTemperature(ix,it) = real(temp_t8p, kind=4) ! t8p //temp_t8p
      else ! skipped pixels will satisfy this condition
         out_SlantColumnAmountO2O2(ix,it) = nasa_SlantColumnAmountO2O2(ix,it)
         out_O2O2CloudTemperature(ix,it) = TrefO4
      endif

      out_CloudPressureNotClipped(ix,it)= cpp 

      !-----------
      ! calculate lut_pcld index from cpp
      ! bad cpp will return negative index_pcld_lut 
      ! lut_pcld_indarr will be used for ecf in next iteration
      ! negative elements will revert back to pcld=700hPa assumption there
      !-----------
      call find_pcld_lutind(cpp,index_pcld_lut)
      lut_pcld_indarr(ix,it) = index_pcld_lut

      out_CloudPressure(ix,it)= cpp 
      if ((cpp .le. 0.).or.(cpp .ge. lut_psfc(npsfc))) then
        out_CloudPressure(ix,it)= fFillValue 
        out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),13)
      endif

      if (option_clip_pcld .eq. 'yes') then
         if ((cpp .gt. psfc0).and.(cpp .le. lut_psfc(npsfc))) then
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),14)
             out_CloudPressure(ix,it) = psfc0 
         endif
         if ((cpp .gt. 0.).and.(cpp .lt. lut_pcld(1))) then
             out_ProcessingQualityFlags(ix,it)=ibset(out_ProcessingQualityFlags(ix,it),14)
             out_CloudPressure(ix,it) = lut_pcld(1) 
         endif
      endif

      ! safeguard ocp for anything with processing problem (bit 13)  
      if (btest(out_ProcessingQualityFlags(ix,it),13)) then
         out_CloudPressure(ix,it) = fFillValue
      endif

 3456 continue ! skip to here when ocp does not need re-calculation

      !=====
    end do
  end do
  !=====

  ! deallocate allocated local variables
  deallocate(pp, tt)

  ! debug
  close(lun_debug_scdadj)

  if ((trim(run_mode).eq.'development') &
      .and.(itdebug.ge. 0).and.(ixdebug.ge.0)) then
     write(*,*) '  writing debug_pcldind.txt'
     open(unit=lun_debug_pcldind,file='debug_pcldind.txt')
     write(lun_debug_pcldind,*)'it=',ix,itdebug
     write(lun_debug_pcldind,*) 'OCP(hPa), LUT_pcld_ind'
     do ix = 1, nx
        cpp = out_CloudPressure(ix,itdebug)
        index_pcld_lut = lut_pcld_indarr(ix,itdebug) 
        write(lun_debug_pcldind,*)cpp, index_pcld_lut
     enddo
  close(lun_debug_pcldind)
  endif

!**********************
end subroutine cal_ocp
!**********************

!222222222222222222222222
!**********************
subroutine allocate_ocp_arrays(nx,nt,fFillValue9,ierr)
!**********************
!222222222222222222222222
   use m_vars, only: out_CloudPressure,out_CloudPressureNotClipped,&
        out_SlantColumnAmountO2O2,out_O2O2CloudTemperature
   use m_vars, only: prev_ocp, prev_ocp_notclipped

   implicit none
   integer, intent(in):: nx,nt
   integer, intent(inout)::ierr
   real, intent(in):: fFillValue9

  allocate(out_CloudPressure(nx,nt),stat=ierr)
  allocate(out_CloudPressureNotClipped(nx,nt),stat=ierr)
  allocate(out_SlantColumnAmountO2O2(nx,nt),stat=ierr)
  allocate(out_O2O2CloudTemperature(nx,nt),stat=ierr)

! initialize 
  out_CloudPressure= fFillValue9 
  out_CloudPressureNotClipped= fFillValue9 
  out_SlantColumnAmountO2O2=fFillValue9
  out_O2O2CloudTemperature=fFillValue9

   allocate(prev_ocp(nx,nt),stat=ierr)
   allocate(prev_ocp_notclipped(nx,nt),stat=ierr)
   prev_ocp = fFillValue9
   prev_ocp_notclipped = fFillValue9

end subroutine allocate_ocp_arrays
!**********************

end module m_cal_ocp
