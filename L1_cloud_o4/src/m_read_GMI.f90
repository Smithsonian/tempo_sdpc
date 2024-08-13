!****************
module m_read_GMI
!****************
use m_vars, only: gmi_nx,gmi_ny,gmi_np, gmi_vcd
use m_vars, only: gmi_lon, gmi_lat, gmi_Temperature
use m_vars, only: gmi_Pressure, gmi_TerrainPressure
use m_vars, only: npcld, vcd_convfac, nlayers, lut_pcld
use m_vars, only: gmi_np

contains

!11111111111111111111111111111111111111111111111111
subroutine read_GMI_TMP(gmonth,name_gmi_dir,ierr)
!11111111111111111111111111111111111111111111111111

  use m_vars, only: ilun_gmi_psfc, ilun_gmi_tmp
  
  implicit none
     
  integer, intent(in):: gmonth 
  character(len=255), intent(in):: name_gmi_dir
  integer(kind=4), intent(inout):: ierr

  character(len=255)::name_gmi_psfc
  character(len=255)::name_gmi_tmp 
  character(len=2):: gmonthstr

  character(len=255)::head
  integer(kind=4)::ix,iy,ip
  real(kind=4)::temp,t01,t02,t03,t04,t05,t06

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

!-------------
! filename for gmi climatology
!-------------
    if (gmonth < 10) then
       write(gmonthstr,'(I1)') gmonth
       name_gmi_psfc = trim(name_gmi_dir)//'GMI_Psfc_0'//trim(gmonthstr)//'.txt'
       name_gmi_tmp = trim(name_gmi_dir)//'GMI_TMP_0'//trim(gmonthstr)//'.txt'
    else
       write(gmonthstr,'(I2)') gmonth
       name_gmi_psfc = trim(name_gmi_dir)//'GMI_Psfc_'//trim(gmonthstr)//'.txt'
       name_gmi_tmp = trim(name_gmi_dir)//'GMI_TMP_'//trim(gmonthstr)//'.txt'
    endif

    write(*,*) ' read GMI climatology for psfc & temperature from: '
    write(*,*) '   '//trim(name_gmi_psfc)
    write(*,*) '   '//trim(name_gmi_tmp)

!-------------
! allocate GMI variables
!-------------
    allocate(gmi_lon(gmi_nx), gmi_lat(gmi_ny),stat=ierr)
    allocate(gmi_Temperature(gmi_nx,gmi_ny,gmi_np),stat=ierr)
    allocate(gmi_Pressure(gmi_nx,gmi_ny,gmi_np+1),stat=ierr)
    allocate(gmi_TerrainPressure(gmi_nx,gmi_ny),stat=ierr)

    if (ierr .ne. 0) then
        write(*,*) ' cannot allocate arrays for GMI.'
        return
    endif

    nlayers = gmi_np
!-------------
! define GMI lon/lat
!------------
    do ix=1,gmi_nx
      gmi_lon(ix)=-180.0+1.25*real(ix-1)
    end do

    gmi_lat(1)=-89.75
    do iy=2,gmi_ny-1
      gmi_lat(iy)=-89.0+1.0*real(iy-2)
    end do
    gmi_lat(gmi_ny)=89.75

!--------------
! read GMI Psfc
!--------------
    open(unit=ilun_gmi_psfc,file=trim(name_gmi_psfc),err=121,&
      form='formatted',status='old',action='read',iostat=ierr)

    do ix=1,gmi_nx
    do iy=1,gmi_ny
      read(ilun_gmi_psfc,111,iostat=ierr) temp 
      gmi_TerrainPressure(ix,iy)=temp
    end do
    end do
121 close(unit=ilun_gmi_psfc)

    if (ierr .ne. 0) then
       write(*,*) 'error reading ',trim(name_gmi_psfc)
       return
    endif 
!---------------------
! read GMI Temperature
!---------------------
! temperature defined from top (1) to near surface (72)

    open(unit=ilun_gmi_tmp,file=trim(name_gmi_tmp),err=122,&
      form='formatted',status='old',action='read',iostat=ierr)

    read(ilun_gmi_tmp,*,iostat=ierr,err=122) head
    do ip=1,12
      read(ilun_gmi_tmp,113,iostat=ierr,err=122) t01,t02,t03,t04,t05,t06
    end do
    read(ilun_gmi_tmp,*,iostat=ierr,err=122) head
    do ip=1,12
      read(ilun_gmi_tmp,113,iostat=ierr,err=122) t01,t02,t03,t04,t05,t06
    end do
    read(ilun_gmi_tmp,*,iostat=ierr,err=122) head

    do iy=1,gmi_ny
      do ix=1,gmi_nx
      do ip=1,12
      read(ilun_gmi_tmp,115,iostat=ierr,err=122) t01,t02,t03,t04,t05,t06
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-1)=t01
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-2)=t02
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-3)=t03
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-4)=t04
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-5)=t05
      gmi_Temperature(ix,iy,gmi_np+1-6*(ip-1)-6)=t06
      end do
      end do
    end do
122  close(unit=ilun_gmi_tmp,iostat=ierr)

    if (ierr .ne. 0) then 
       write(*,*)'error reading ',trim(name_gmi_tmp)
       return
    endif
!-----------------------
! calculate GMI Pressure
!-----------------------
! Pressure defined from top (1) to near surface (72)
! The terrain pressure is defined at pressure level 73 

  do ix=1,gmi_nx
  do iy=1,gmi_ny
    do ip=1,gmi_np
      gmi_Pressure(ix,iy,gmi_np+1-ip)=gmi_am(ip)+gmi_bm(ip)*gmi_TerrainPressure(ix,iy)
    end do
    gmi_Pressure(ix,iy,gmi_np+1)=gmi_TerrainPressure(ix,iy)
  end do
  end do

111 format(10x,f10.3)
113 format(6f13.7)
115 format(6f13.3)

!11111111111111111111111111
end subroutine read_GMI_TMP
!11111111111111111111111111


!22222222222222222222222222222
subroutine read_GMI_VCD(pp,tt)
!22222222222222222222222222222

  implicit none

  real(kind=4),dimension(gmi_np)::tt
  real(kind=4),dimension(gmi_np+1)::pp !include Psfc
  real(kind=4),dimension(gmi_np+1)::tmp_vcd
  real::sum_vcd
  real::xx1,xx2,yy1,yy2,xxx,yyy
  integer(kind=4)::iflag,ip,ipcld

! -----------------------
! GMI pressure coordinate
! -----------------------
  tmp_vcd(1)=vcd_convfac/2.0/tt(1)*(pp(1)**2)
  sum_vcd=tmp_vcd(1)
  do ip=1,gmi_np
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
      gmi_vcd(ipcld)=yyy
    else
      sum_vcd=vcd_convfac/2.0/tt(gmi_np)*(lut_pcld(ipcld)**2-lut_pcld(ipcld-1)**2)
      gmi_vcd(ipcld)=gmi_vcd(ipcld-1)+sum_vcd
    endif
  end do

!22222222222222222222222222
end subroutine read_GMI_VCD
!22222222222222222222222222
    
!33333333333333333333333333
subroutine get_GMIpsfc_lonlat(lon0, lat0, psfcout)
!33333333333333333333333333
! get GMI psfc at TEMPO pixel (ix,it)

    implicit none

    real, intent(IN):: lon0, lat0
    real, intent(OUT):: psfcout

    integer :: gmi_ix1,gmi_ix2,gmi_iy1,gmi_iy2
    real:: gmi_wx1, gmi_wx2, gmi_wy1, gmi_wy2
    real::pp11,pp12,pp21,pp22,pp1,pp2

        gmi_wx1 = 0.
        gmi_wx2 = 0.
        gmi_wy1 = 0.
        gmi_wy2 = 0.

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

        pp11=gmi_TerrainPressure(gmi_ix1,gmi_iy1)
        pp12=gmi_TerrainPressure(gmi_ix1,gmi_iy2)
        pp21=gmi_TerrainPressure(gmi_ix2,gmi_iy1)
        pp22=gmi_TerrainPressure(gmi_ix2,gmi_iy2)
        pp1=(gmi_wy2*pp11+gmi_wy1*pp12)/(gmi_wy1+gmi_wy2)
        pp2=(gmi_wy2*pp21+gmi_wy1*pp22)/(gmi_wy1+gmi_wy2)

        psfcout=(gmi_wx2*pp1+gmi_wx1*pp2)/(gmi_wx1+gmi_wx2)

!33333333333333333333333333
end subroutine get_GMIpsfc_lonlat
!33333333333333333333333333

!44444444444444444444444444
subroutine get_GMItmp_lonlat(lon0,lat0,tt,pp,gmi_psfc)
!44444444444444444444444444
! get GMI T-P profile at lon0 & lat0

    implicit none

    real, intent(IN):: lon0, lat0
    real, intent(OUT):: gmi_psfc
    real, dimension(gmi_np), intent(out):: tt
    real, dimension(gmi_np+1), intent(out)::pp

    integer :: gmi_ix1,gmi_ix2,gmi_iy1,gmi_iy2
    real:: gmi_wx1, gmi_wx2, gmi_wy1, gmi_wy2
    real::pp11,pp12,pp21,pp22,pp1,pp2
    real::tt11,tt12,tt21,tt22,tt1,tt2

    integer :: ip

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

        pp11=gmi_TerrainPressure(gmi_ix1,gmi_iy1)
        pp12=gmi_TerrainPressure(gmi_ix1,gmi_iy2)
        pp21=gmi_TerrainPressure(gmi_ix2,gmi_iy1)
        pp22=gmi_TerrainPressure(gmi_ix2,gmi_iy2)
        pp1=(gmi_wy2*pp11+gmi_wy1*pp12)/(gmi_wy1+gmi_wy2)
        pp2=(gmi_wy2*pp21+gmi_wy1*pp22)/(gmi_wy1+gmi_wy2)

        gmi_psfc=(gmi_wx2*pp1+gmi_wx1*pp2)/(gmi_wx1+gmi_wx2)

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

!444444444444444444444
end subroutine get_GMItmp_lonlat
!444444444444444444444

!********************
end module m_read_GMI
!********************
