!******************
module m_read_hdf5
!******************
! this module contains misc subroutines

contains

!111111111111111111111
!subroutine read_GEOS5
!111111111111111111111

!hqw read_GEOS5 is replaced with read_geoscf, thus remove it

!1111111111111111111111111
!end subroutine read_GEOS5
!1111111111111111111111111

!2222222222222222222222222222222
subroutine read_GEOS5_VCD(pp,tt)
!2222222222222222222222222222222

! for an individual pixel 

  use m_vars, only : vcd_convfac, geos_np, geos_vcd, &
                    npcld, lut_pcld

  implicit none

  real(kind=4),dimension(geos_np)::tt
  real(kind=4),dimension(geos_np+1)::pp !include Psfc
  real(kind=4),dimension(geos_np+1)::tmp_vcd
  real::sum_vcd
  real::xx1,xx2,yy1,yy2,xxx,yyy
  integer(kind=4)::iflag,ipcld
  integer(kind=4)::ip

! ---------------------------
! GEOS-5 pressure coordinate
! ---------------------------
  !hqw tmp_vcd(1)=(6.765e-4)/2.0/tt(1)*(pp(1)**2)
  tmp_vcd(1) = vcd_convfac/2.0/tt(1)*(pp(1)**2)
  sum_vcd=tmp_vcd(1)
  do ip=1,geos_np
    !hqw assumes pp on levels, tt between levels
    !sum_vcd=sum_vcd+(6.765e-4)/2.0/tt(ip)*(pp(ip+1)**2-pp(ip)**2)
    sum_vcd=sum_vcd+vcd_convfac/2.0/tt(ip)*(pp(ip+1)**2-pp(ip)**2)
    tmp_vcd(ip+1)=sum_vcd
  end do

! -----------------------
! LUT pressure coordinate
! -----------------------
  do ipcld=1,npcld
    xxx=lut_pcld(ipcld)

    iflag=0
    do ip=1,geos_np
      if((xxx.gt.pp(ip)).and.(xxx.le.pp(ip+1))) then
        xx1=pp(ip)
        xx2=pp(ip+1)
        yy1=tmp_vcd(ip)
        yy2=tmp_vcd(ip+1)
        yyy=(yy1-yy2)/(xx1-xx2)*xxx+(xx1*yy2-xx2*yy1)/(xx1-xx2)
        iflag=iflag+1
      endif
    end do

    if(iflag.ge.1) then
      geos_vcd(ipcld)=yyy
    else
!!!!  if lut_pcld(ipcld) > Psfc, use T at the bottom layer
    !  sum_vcd=(6.765e-4)/2.0/tt(geos_np)*(lut_pcld(ipcld)**2-lut_pcld(ipcld-1)**2)
      sum_vcd=vcd_convfac/2.0/tt(geos_np)*(lut_pcld(ipcld)**2-lut_pcld(ipcld-1)**2)
      geos_vcd(ipcld)=geos_vcd(ipcld-1)+sum_vcd
    endif
  end do

!22222222222222222222222222222
end subroutine read_GEOS5_VCD
!22222222222222222222222222222

!hqw GLER is read in m_read_input_gler.f90
!the following is no longer needed, thus removed
!5555555555555555555555555555
!subroutine read_BRDF_Rsfc_h5
!5555555555555555555555555555
 
!55555555555555555555555555555555
!end subroutine read_BRDF_Rsfc_h5
!55555555555555555555555555555555

!hqw GLER now handled by m_read_input_gler.f90
!the following is no longer used, thus removed
!6666666666666666666666666666
!subroutine read_BDEM_Psfc_h5
!6666666666666666666666666666

!66666666666666666666666666666666
!end subroutine read_BDEM_Psfc_h5
!66666666666666666666666666666666

!**********************
end module m_read_hdf5
!**********************
