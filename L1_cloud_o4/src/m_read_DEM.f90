
module m_read_DEM
!****************
!hqw this module is not actually needed by TEMPO
use m_vars

contains

!1111111111111111111111111111111111111111
subroutine read_DEM_GMI_TMP(name_gmi_tmp)
!1111111111111111111111111111111111111111

  implicit none

  character(len=255)::name_gmi_tmp

  character(len=255)::head
  integer(kind=4)::ix,iy,ip, ierr
  real(kind=4)::t01,t02,t03,t04,t05,t06
!---------------------
!hqw this option is not tested
!---------------------
   if (.not. allocated(gmi_lon)) then
      allocate(gmi_lon(gmi_nx), gmi_lat(gmi_ny),stat=ierr)
      allocate(gmi_Temperature(gmi_nx,gmi_ny,gmi_np),stat=ierr)
   endif
!---------------------
!hqw define GMI lat/lon
!--------------------
if (gmi_lon(0) .lt. -180.) then
    do ix=1,gmi_nx
      gmi_lon(ix)=-180.0+1.25*real(ix-1)
    end do

    gmi_lat(1)=-89.75
    do iy=2,gmi_ny-1
      gmi_lat(iy)=-89.0+1.0*real(iy-2)
    end do
    gmi_lat(gmi_ny)=89.75
endif
!---------------------
! read GMI Temperature
!---------------------
! temperature defined from top (1) to near surface (72)
!   because VCD will be calculated from TOA to BOA
  open(unit=93,file=trim(name_gmi_tmp),form='formatted',status='old',action='read')

  read(93,*) head
  do ip=1,12
    read(93,113) t01,t02,t03,t04,t05,t06
  end do
  read(93,*) head
  do ip=1,12
    read(93,113) t01,t02,t03,t04,t05,t06
  end do
  read(93,*) head

  do iy=1,gmi_ny
    do ix=1,gmi_nx
    do ip=1,12
      read(93,115) t01,t02,t03,t04,t05,t06
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-1)=t01
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-2)=t02
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-3)=t03
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-4)=t04
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-5)=t05
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-6)=t06
    end do
    end do
  end do
  close(unit=93)

113 format(6f13.7)
115 format(6f13.3)

!111111111111111111111111111111
end subroutine read_DEM_GMI_TMP
!111111111111111111111111111111


!2222222222222222222222222222222
!hqw added pp to the arg list to facilitate scd T-corr
subroutine read_DEM_VCD(psfc,tt,pp)
!2222222222222222222222222222222

! for an individual pixel 

  implicit none

  real(kind=4),dimension(gmi_np)::tt
  real(kind=4),dimension(gmi_np+1)::pp !include Psfc
  real(kind=4),dimension(gmi_np+1)::tmp_vcd
  real::psfc,sum_vcd
  real::xx1,xx2,yy1,yy2,xxx,yyy
  !real::x1,x2,x3,y1,y2,y3
  integer(kind=4)::iflag,ipcld!,iline

  integer(kind=4)::ip!,ix,iy
  real,dimension(gmi_np):: &
    gmi_am=(/0.0240240, 3.32090,  9.86428, 16.3740,  22.8526,  29.3314, &
            35.7764,   42.1605,  48.5176,  54.8747,  61.2079,  67.4923, &
            74.7281,   83.9671,  94.2326, 104.273,  114.070,  123.827, &
           135.803,   149.585,  162.934,  175.614,  187.358,  198.178, &
           207.704,   215.463,  221.337,  224.131,  220.614,  209.029, &
           189.061,   163.662,  139.115,  118.250,  100.514,   85.4390, &
            72.5579,   61.4957,  52.0159,  43.9097,  36.9927,  31.0889, &
            26.0491,   21.7610,  18.1244,  15.0503,  12.4602,  10.2849, &
             8.45639,   6.91832,  5.63180,  4.56169,  3.67650,  2.94832, &
             2.35259,   1.86788,  1.47565,  1.15998,  0.907287, 0.705957, &
             0.546293,  0.420424, 0.321783, 0.244938, 0.185422, 0.139599, &
             0.104524,  0.0776725,0.0567925,0.0401425,0.0263500,0.0150000/)
  real,dimension(gmi_np):: &
    gmi_bm=(/0.992500,0.974200,0.952650,0.931150,0.909650,0.888150, &
             0.866700,0.845350,0.824000,0.802600,0.781250,0.760000, &
             0.735300,0.703550,0.668250,0.633200,0.598400,0.563650, &
             0.520450,0.469150,0.418300,0.368150,0.318900,0.270550, &
             0.223550,0.178300,0.134900,0.088650,0.045850,0.017500, &
             0.003500,0.000000,0.000000,0.000000,0.000000,0.000000, &
             0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
             0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
             0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
             0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
             0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
             0.000000,0.000000,0.000000,0.000000,0.000000,0.000000/)

!-----------------------
! calculate GMI Pressure
!-----------------------
! Pressure defined from top (1) to near surface (72)
! The terrain pressure is defined at pressure level 73 

  do ip=1,gmi_np
    pp(gmi_np+1-ip)=gmi_am(ip)+gmi_bm(ip)*psfc
  end do
  pp(gmi_np+1)=psfc

! -----------------------
! GMI pressure coordinate
! -----------------------
  !tmp_vcd(1)=(6.765e-4)/2.0/tt(1)*(pp(1)**2)
  tmp_vcd(1)=vcd_convfac/2.0/tt(1)*(pp(1)**2)
  sum_vcd=tmp_vcd(1)
  do ip=1,gmi_np
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
    do ip=1,gmi_np
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
      dem_vcd(ipcld)=yyy
!      write(3,*) ipcld,xxx,dem_vcd(ipcld)
    else
!!!!  if lut_pcld(ipcld) > Psfc, use T at the bottom layer
! FIXME: There was sloppy code here, treating the bottom layer as if there
!        was a layer below it. I am assuming that layer was supposed to have
!        a zero value, but who knows?
      if (ipcld > 1) then
        !sum_vcd=(6.765e-4)/2.0/tt(gmi_np)*(lut_pcld(ipcld)**2-lut_pcld(ipcld-1)**2)
        sum_vcd=vcd_convfac/2.0/tt(gmi_np)*(lut_pcld(ipcld)**2-lut_pcld(ipcld-1)**2)
        dem_vcd(ipcld)=dem_vcd(ipcld-1)+sum_vcd
      else
        dem_vcd=0.0
      endif
    endif
  end do


!22222222222222222222222222
end subroutine read_DEM_VCD
!22222222222222222222222222
    
 
!********************
end module m_read_DEM
!********************
