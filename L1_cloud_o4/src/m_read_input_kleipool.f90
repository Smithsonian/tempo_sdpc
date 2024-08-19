!*******************
module m_read_input_kleipool
!*******************

contains

  !111111111111111111111111111111111111111111
  subroutine read_Kleipool_Rsfc(filename,m12,ierr)
  !111111111111111111111111111111111111111111

    use hdfeos4_parameters
    use he5_swreader
    use m_vars, only: kleipool_lon, kleipool_lat
    use m_vars, only: kleipool_SurfaceReflectivity466
    use m_vars, only: kleipool_SurfaceReflectivity440
    use m_vars, only: kleipool_SurfaceReflectivity477
    use m_vars, only: kleipool_nx, kleipool_ny
    use m_vars, only: PerturbAlb466, nord_Alb466pert, Alb466PertCoef

    implicit none

    !---------------------------------------------------------------------72
    ! ROUTINE: read_Kleipool_Rsfc
    !
    ! DESCRIPTION: This program reads Earth surface reflectance climatology
    !
    ! REVISION HISTORY:
    !  10/07/15 Yang original fortran 90
    !---------------------------------------------------------------------72

    ! Data Type                                Element Name
    ! ==============================           ============
    ! Data fields
    !  real(kind=4),   dimension(:,:),pointer::Latitude
    !  real(kind=4),   dimension(:,:),pointer::Longitude
    !  real(kind=4),   dimension(:),  pointer::Wavelength
    !  integer(kind=2),dimension(:,:),pointer::MonthlySurfaceReflectance
    !  integer(kind=1),dimension(:,:),pointer::MonthlySurfaceReflectanceFlag
    !  integer(kind=2),dimension(:,:),pointer::MonthlyMinimumSurfaceReflectance

    integer, intent(inout):: ierr

    integer(kind=4)::he5_gdopen,he5_gdattach,he5_gdrdfld,&
         he5_gddetach,he5_gdclose

    character(len=255), intent(in)::filename
    integer (kind=4), intent(in)::m12
    integer::ix,iy,iord
    integer::nx,ny,nw,n12
    integer(kind=4)::gdfid,gdid,status
    integer(kind=2),dimension(:,:,:,:),pointer::rtemp
    real(kind=4)::r440,r463,r466,r471,r477
    real(kind=4)::xalb,yalb

    integer(kind=8)::start1,stride1,edge1
    !integer(kind=8),dimension(2)::start2,stride2,edge2
    integer(kind=8),dimension(4)::start4,stride4,edge4

    nx=kleipool_nx
    ny=kleipool_ny
    nw=23
    n12=12

    !----------------------------------
    !  allocate kleipool variables
    !  477 is not actually used
    !  currently 440 is only used to assign out_SurfaceReflectivity440
    !  only 466 is used in calculation
    !----------------------------------
    allocate(kleipool_lon(nx), stat=ierr)
    allocate(kleipool_lat(nx), stat=ierr)
    allocate(kleipool_SurfaceReflectivity466(nx,ny), stat=ierr)
    allocate(kleipool_SurfaceReflectivity440(nx,ny), stat=ierr)
    allocate(kleipool_SurfaceReflectivity477(nx,ny), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*) '   Error allocating kleipool arrays.'
       return
    endif

    write(*,*) '   reading Kleipool rsfc '//trim(filename)
    !-----------------------
    ! read OMI-Aura_L3-OMLER
    !-----------------------
    gdfid=he5_gdopen(filename,HE5F_ACC_RDONLY)
    gdid=he5_gdattach(gdfid,'EarthSurfaceReflectanceClimatology')
    if (gdfid < 0 .or. gdid < 0) write(*,*) "   read_Kleipool_Rsfc failed"

    !----------
    ! Array(nx)
    !----------
    start1=0
    stride1=1
    edge1=nx
    status=he5_gdrdfld(gdid,"Longitude",start1,stride1,edge1,kleipool_lon)

    !----------
    ! Array(ny)
    !----------
    start1=0
    stride1=1
    edge1=ny
    status=he5_gdrdfld(gdid,"Latitude",start1,stride1,edge1,kleipool_lat)

    !--------------------
    ! Array(nx,ny,nw,n12)
    !--------------------
    start4(1)=0
    stride4(1)=1
    edge4(1)=nx
    start4(2)=0
    stride4(2)=1
    edge4(2)=ny
    start4(3)=0
    stride4(3)=1
    edge4(3)=nw
    start4(4)=0
    stride4(4)=1
    edge4(4)=n12

    ! allocate temporary array
    allocate(rtemp(nx,ny,nw,n12),stat=ierr)

    status=he5_gdrdfld(gdid,"MonthlySurfaceReflectance",start4,stride4,edge4,rtemp)

    !-------
    ! close
    !-------
    status=he5_gddetach(gdid)
    status=he5_gdclose(gdfid)

    !----------------------------------
    ! LER for the month 
    !    w15=440 w18=463 w19=471 w20=477
    !----------------------------------
    do ix=1,nx
      do iy=1,ny
        r477=real(rtemp(ix,iy,20,m12))/1000.
        r471=real(rtemp(ix,iy,19,m12))/1000.
        r463=real(rtemp(ix,iy,18,m12))/1000.
        r440=real(rtemp(ix,iy,15,m12))/1000.
        if(r477 .gt. 1.0) r477=1.0
        if(r471 .gt. 1.0) r471=1.0
        if(r463 .gt. 1.0) r463=1.0
        if(r440 .gt. 1.0) r440=1.0
        kleipool_SurfaceReflectivity477(ix,iy)=r477
        kleipool_SurfaceReflectivity440(ix,iy)=r440
        r466=(r463*5.+r471*3.)/8. ! interpolate to 466
        kleipool_SurfaceReflectivity466(ix,iy)=r466
      end do
    end do

  ! deallocate rtemp for memory
  deallocate(rtemp)

  ! Perturb r466
  if (PerturbAlb466) then
     write(*,*) '   Apply Alb466PertCoef to Kleipool 466'
     do ix = 1, nx
        do iy = 1, ny
           xalb = kleipool_SurfaceReflectivity466(ix,iy)
           if ((xalb .ge. 0.) .and. (xalb .le. 1.)) then
              yalb = 0. 
              do iord = 0, nord_Alb466Pert
                 yalb = yalb + Alb466PertCoef(iord+1)*(xalb**iord)
              enddo
              if (yalb .gt. 1.) yalb = 1.0
              if (yalb .lt. 0.) yalb = 0.0
           else if (xalb .lt. 0.) then
              yalb = 0.0 
           else
              yalb = 1.0
           endif
           kleipool_SurfaceReflectivity466(ix,iy) = yalb
         enddo !iy
      enddo !ix
   endif
         
  !111111111111111111111111111111111
  end subroutine read_Kleipool_Rsfc
  !111111111111111111111111111111111

  !222222222222222222222222222222222
  subroutine get_kleipool_lonlat(lon0,lat0,rsfc466out,&
             rsfc440out,rsfc477out)
  !222222222222222222222222222222222
  ! get Kleipool rsfc at TEMPO pixel (ix,it)

  use m_vars, only: kleipool_SurfaceReflectivity466,kleipool_lon
  use m_vars, only: kleipool_SurfaceReflectivity440,kleipool_lat
  use m_vars, only: kleipool_SUrfaceReflectivity477
  use m_vars, only: kleipool_nx, kleipool_ny
 
  implicit none
  real, intent(in):: lon0, lat0
  real, intent(out):: rsfc466out, rsfc440out, rsfc477out

  real:: kle_wx1,kle_wx2,kle_wy1,kle_wy2
  real:: rsfc11,rsfc12,rsfc21,rsfc22,rsfc1,rsfc2
  integer:: kle_ix1,kle_ix2,kle_iy1,kle_iy2

  kle_wx1 = 0. 
  kle_wx2 = 0.
  kle_wy1 = 0.
  kle_wy2 = 0.

  kle_ix1 = floor((lon0+180.)/0.5)+1
  kle_ix2 = kle_ix1+1
  kle_iy1 = floor((lat0+90.)/0.5)+1
  kle_iy2 = kle_iy1+1

  if (kle_ix1 .lt. 1) kle_ix1 = 1
  if (kle_ix1 .gt. kleipool_nx) kle_ix1 = kleipool_nx
  if (kle_ix2 .lt. 1) kle_ix2 = 1
  if (kle_ix2 .gt. kleipool_nx) kle_ix2 = kleipool_nx
  if (kle_iy1 .lt. 1) kle_iy1 = 1
  if (kle_iy1 .gt. kleipool_ny) kle_iy1 = kleipool_ny
  if (kle_iy2 .lt. 1) kle_iy2 = 1
  if (kle_iy2 .gt. kleipool_ny) kle_iy2 = kleipool_ny

  kle_wx1 = lon0 - kleipool_lon(kle_ix1)
  kle_wx2 = kleipool_lon(kle_ix2) - lon0
  kle_wy1 = lat0 - kleipool_lat(kle_iy1)
  kle_wy2 = kleipool_lat(kle_iy2) - lat0

  rsfc11 = kleipool_SurfaceReflectivity466(kle_ix1,kle_iy1)
  rsfc12 = kleipool_SurfaceReflectivity466(kle_ix1,kle_iy2)
  rsfc21 = kleipool_SurfaceReflectivity466(kle_ix2,kle_iy1)
  rsfc22 = kleipool_SurfaceReflectivity466(kle_ix2,kle_iy2)
  rsfc1 = (kle_wy2*rsfc11+kle_wy1*rsfc12)/(kle_wy1+kle_wy2)
  rsfc2 = (kle_wy2*rsfc21+kle_wx1*rsfc22)/(kle_wy1+kle_wy2)
  rsfc466out=(kle_wx2*rsfc1+kle_wx1*rsfc2)/(kle_wx1+kle_wx2)
  
  rsfc11 = kleipool_SurfaceReflectivity440(kle_ix1,kle_iy1)
  rsfc12 = kleipool_SurfaceReflectivity440(kle_ix1,kle_iy2)
  rsfc21 = kleipool_SurfaceReflectivity440(kle_ix2,kle_iy1)
  rsfc22 = kleipool_SurfaceReflectivity440(kle_ix2,kle_iy2)
  rsfc1 = (kle_wy2*rsfc11+kle_wy1*rsfc12)/(kle_wy1+kle_wy2)
  rsfc2 = (kle_wy2*rsfc21+kle_wx1*rsfc22)/(kle_wy1+kle_wy2)
  rsfc440out =(kle_wx2*rsfc1+kle_wx1*rsfc2)/(kle_wx1+kle_wx2)
  
  rsfc11 = kleipool_SurfaceReflectivity477(kle_ix1,kle_iy1)
  rsfc12 = kleipool_SurfaceReflectivity477(kle_ix1,kle_iy2)
  rsfc21 = kleipool_SurfaceReflectivity477(kle_ix2,kle_iy1)
  rsfc22 = kleipool_SurfaceReflectivity477(kle_ix2,kle_iy2)
  rsfc1 = (kle_wy2*rsfc11+kle_wy1*rsfc12)/(kle_wy1+kle_wy2)
  rsfc2 = (kle_wy2*rsfc21+kle_wx1*rsfc22)/(kle_wy1+kle_wy2)
  rsfc477out =(kle_wx2*rsfc1+kle_wx1*rsfc2)/(kle_wx1+kle_wx2)

  !2222222222222222222222222
  end subroutine get_kleipool_lonlat
  !2222222222222222222222222

!***********************
end module m_read_input_kleipool
!***********************
