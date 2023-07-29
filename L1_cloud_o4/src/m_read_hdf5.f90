!******************
module m_read_hdf5
!******************
! this module contains misc subroutines

contains

!111111111111111111111
!subroutine read_GEOS5
!111111111111111111111

! read_GEOS5 is replaced with read_geoscf, thus remove it

!1111111111111111111111111
!end subroutine read_GEOS5
!1111111111111111111111111

!2222222222222222222222222222222
subroutine read_GEOS5_VCD(ptot,tt,qq,ppdry)
!2222222222222222222222222222222
! for individual pixel 
  use m_vars, only : vcd_convfac, geos_np, geos_vcd, &
                    npcld, lut_pcld
  use m_vars, only : name_option_adjdry

  implicit none

  real(kind=4),dimension(:),intent(in)::tt  ! tt(geos_np)
  real(kind=4),dimension(:),intent(in)::qq  ! qq(geos_np)
  real(kind=4),dimension(:),intent(in)::ptot ! pp(geos_np+1) include Psfc
  real(kind=4),dimension(size(ptot)),intent(out)::ppdry 

  ! local variable
  real(kind=4),dimension(size(ptot))::tmp_vcd, pp
  real::sum_vcd
  real::xx1,xx2,yy1,yy2,xxx,yyy

  integer(kind=4)::iflag,ipcld
  integer(kind=4)::ip, nlayer, nlevel

! ptot is total pressure (in cal_ocp and cal_pscene)
! ---------------------------
  if (name_option_adjdry .EQ. 1) then 
    nlayer = geos_np
    nlevel = nlayer + 1
    call totp_to_dryp(ptot,qq,nlevel,nlayer,ppdry)
    pp = ppdry
  else ! not adjust to dry pressure
     ppdry = ptot
     pp = ppdry
  endif

! pp here is a local variable
! different than the pp in cal_ocp & cal_pscene
! there pp refers to total pressure

! ---------------------------
! GEOS-5 pressure coordinate
! ---------------------------
  tmp_vcd(1) = vcd_convfac/2.0/tt(1)*(pp(1)**2)
  sum_vcd=tmp_vcd(1)
  do ip=1,geos_np
    ! assumes ppdry on levels, tt between levels
    sum_vcd=sum_vcd+vcd_convfac/2.0/tt(ip)*(pp(ip+1)**2-pp(ip)**2)
    tmp_vcd(ip+1)=sum_vcd
  end do

! LUT is based on US-standard air which is dry
! pp is best provided as dry pressure instead of total pressure
! the following convert tmp_vcd on GEOS levels to geos_vcd on LUT levels
! geos_vcd is converted to SCD through LUT AMF in cal_ocp and cal_pscene
! and compare with retrieved O4 SCD
! -----------------------
! LUT pressure coordinate
! -----------------------
! current min(pcld) in LUT is 55hPa
! GEOS-CF TOA is <0.01 hPa
  do ipcld=1,npcld
    xxx=lut_pcld(ipcld)

    iflag=0
    do ip=1,geos_np
      if ((xxx.gt.pp(ip)) .and. (xxx.le.pp(ip+1))) then
        xx1=pp(ip)
        xx2=pp(ip+1)
        yy1=tmp_vcd(ip)
        yy2=tmp_vcd(ip+1)
        yyy=(yy1-yy2)/(xx1-xx2)*xxx+(xx1*yy2-xx2*yy1)/(xx1-xx2)
        iflag=iflag+1
        exit
      endif
    end do

    if (iflag.ge.1) then ! node found above
!!! these have lut_pcld(ipcld) .le. psfc, use the interpolated
      geos_vcd(ipcld)=yyy
    else ! node not found
!!! should be because lut_pcld(ipcld) .gt. Psfc
!!! this could happen for high topography
!!! use T at BOA to extropolate and accumulate upon previous ipcld 
!!! i.e. add ghost layers vcd beneath the GEOS psfc and max LUT pressure
      sum_vcd=vcd_convfac/2.0/tt(geos_np)*(lut_pcld(ipcld)**2-lut_pcld(ipcld-1)**2)
      geos_vcd(ipcld)=geos_vcd(ipcld-1)+sum_vcd
    endif
  end do

!22222222222222222222222222222
end subroutine read_GEOS5_VCD
!22222222222222222222222222222

!33333333333333333333333333333
subroutine totp_to_dryp(pp,qq,nlevel,nlayer,ppdry)
!33333333333333333333333333333
   ! nlevel is dimension of pp, nlayer for qq
   ! pp is from TOA to BOA, so is qq [kg H2O/ kg dry air]
   ! pp includes surface, qq does not
   ! qq is for each layer between 2 levels

   integer, intent(in):: nlevel, nlayer
   real(kind=4), dimension(nlevel), intent(in):: pp
   real(kind=4), dimension(nlayer), intent(in):: qq
   ! ppdry is in the same unit as pp
   real(kind=4), dimension(nlevel), intent(out):: ppdry

   integer :: iz
   real(kind=4) :: detptot, detpdry, sumdry

   ! initial
   do iz = 1, geos_np+1
      ppdry(iz) = pp(iz) 
   enddo

   if (pp(2) .gt. pp(3)) then ! should not happen, safeguard
      write(*,*) 'totp_to_dryp: pp should increase, not decrease'
      write(*,*) '****set ppdry = pp, nothing is done!!!'
      return
   endif 

   ! ptot = pdry + pwet = pdry + pdry * qq
   ppdry(1) = pp(1) !TOA
   sumdry = ppdry(1)
   do iz = 1, nlayer 
      detptot = pp(iz+1) - pp(iz)
      detpdry = detptot / (1. + qq(iz)) 
      sumdry = sumdry + detpdry ! accumulate
      ppdry(iz+1) = sumdry
   enddo

!3333333333333333333333333333
end subroutine totp_to_dryp
!3333333333333333333333333333

!5555555555555555555555555555
!subroutine read_BRDF_Rsfc_h5
!5555555555555555555555555555

! GLER is read in m_read_input_gler.f90
!the following is no longer needed, thus removed
 
!55555555555555555555555555555555
!end subroutine read_BRDF_Rsfc_h5
!55555555555555555555555555555555

!**********************
end module m_read_hdf5
!**********************
