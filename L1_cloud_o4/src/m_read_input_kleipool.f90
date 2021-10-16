!*******************
module m_read_input_kleipool
  !*******************

  use hdfeos4_parameters
  use he5_swreader
  use m_vars

contains

!hqw read_input is no longer needed, thus deleted

!  111111111111111111111
!  subroutine read_input
!  111111111111111111111


    !1111111111111111111111111
!  end subroutine read_input
  !1111111111111111111111111


!hqw read_OMO4SCD is no longer needed, thus deleted
  !22222222222222222222222
!  subroutine read_OMO4SCD
    !22222222222222222222222

    !222222222222222222222222222
!  end subroutine read_OMO4SCD
  !222222222222222222222222222


  !333333333333333333333333333333333333333333
  subroutine read_Kleipool_Rsfc(filename,m12)
    !333333333333333333333333333333333333333333

    use hdfeos4_parameters
    use he5_swreader
    use m_vars

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

    integer(kind=4)::he5_gdopen,he5_gdattach,he5_gdrdfld,&
         he5_gddetach,he5_gdclose

    character(len=255), intent(in)::filename
    integer (kind=4), intent(in)::m12
    integer::ix,iy, ierr
    integer::nx,ny,nw,n12
    integer(kind=4)::gdfid,gdid,status
    integer(kind=2),dimension(:,:,:,:),pointer::rtemp
    real(kind=4)::r440,r463,r466,r471,r477

    integer(kind=8)::start1,stride1,edge1
    !integer(kind=8),dimension(2)::start2,stride2,edge2
    integer(kind=8),dimension(4)::start4,stride4,edge4

    nx=kleipool_nx
    ny=kleipool_ny
    nw=23
    n12=12

    !----------------------------------
    ! hqw allocate kleipool variables
    ! in fact, 477 is not actually used
    !          440 is only used to assign out_SurfaceReflectivity440
    !     only 466 is used in calculation 
    !----------------------------------
    allocate(kleipool_lon(nx), stat=ierr)
    allocate(kleipool_lat(nx), stat=ierr)
    allocate(kleipool_SurfaceReflectivity466(nx,ny), stat=ierr)
    allocate(kleipool_SurfaceReflectivity440(nx,ny), stat=ierr)
    !hqw as 477 is not actually used, comment out to save memory
    !allocate(kleipool_SurfaceReflectivity477(nx,ny), stat=ierr)

    write(*,*) '   reading Kleipool rsfc '//trim(filename)
    !-----------------------
    ! read OMI-Aura_L3-OMLER
    !-----------------------
    gdfid=he5_gdopen(filename,HE5F_ACC_RDONLY)
    gdid=he5_gdattach(gdfid,'EarthSurfaceReflectanceClimatology')
    if (gdfid < 0 .or. gdid < 0) print *, "read_Kleipool_Rsfc fail"

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
    allocate(rtemp(nx,ny,nw,n12),stat=ierr)
    status=he5_gdrdfld(gdid,"MonthlySurfaceReflectance",start4,stride4,edge4,rtemp)

    !-------
    ! close
    !-------
    status=he5_gddetach(gdid)
    status=he5_gdclose(gdfid)

    !----------------------------------
    ! interpolate LER for the month i12
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
        !hqw commented 477 out 
        !kleipool_SurfaceReflectivity477(ix,iy)=r477
        kleipool_SurfaceReflectivity440(ix,iy)=r440
        r466=(r463*5.+r471*3.)/8.
        kleipool_SurfaceReflectivity466(ix,iy)=r466
      end do
    end do

    !333333333333333333333333333333333
  end subroutine read_Kleipool_Rsfc
  !333333333333333333333333333333333

!hqw read_BRDF_Rsfc will be replaced by read_gler
  !4444444444444444444444444
  subroutine read_BRDF_Rsfc
  !4444444444444444444444444

    implicit none

    ! DESCRIPTION: This program reads BRDF Rsfc data in txt form

    integer(kind=4)::ix,it,ierr
    integer(kind=4)::nx,nt
    character(len=255)::filename

    !hqw replace
    !nt=inp_NumTimes
    !nx=inp_nXtrack
    nt = rad_NumTimes
    nx = rad_nXtrack

    ! 440 was used in cal_ocp to assign out_SurfaceReflectivity
    ! in theory it should be allocated and assigned here 
    !if (.not. allocated(BRDF_SurfaceReflectivity466)) then 
       allocate(BRDF_SurfaceReflectivity466(nx,nt),stat=ierr)
    !endif

    filename=trim(name_brdf_dir)//trim(name_brdf_file)

    open(unit=94,file=filename,status='old',action='read',iostat=ierr)
    if (ierr /= 0) print *, "read_BRDF_Rsfc fail"

    do it=1,nt
      do ix=1,nx
        read(94,114) BRDF_SurfaceReflectivity466(ix,it)
        if(BRDF_SurfaceReflectivity466(ix,it) .gt. 1.0) BRDF_SurfaceReflectivity466(ix,it)=1.0
      end do
    end do
    close(unit=94)

114 format(f12.4)

  !44444444444444444444444444444
  end subroutine read_BRDF_Rsfc
  !44444444444444444444444444444

!***********************
end module m_read_input_kleipool
!***********************
