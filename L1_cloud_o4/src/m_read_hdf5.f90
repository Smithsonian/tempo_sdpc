!******************
module m_read_hdf5
!******************
! this module contains misc subroutines
use HDF5
use m_vars

contains

!hqw subroutines used inp_Numtimes & inp_nXtrack originally from OMCLDO2
!    changed to rad_NumTimes * rad_Xtrack

!111111111111111111111
subroutine read_GEOS5
!111111111111111111111

!hqw read_GEOS5 is replaced with read_geoscf, thus remove it

!1111111111111111111111111
end subroutine read_GEOS5
!1111111111111111111111111

!2222222222222222222222222222222
subroutine read_GEOS5_VCD(pp,tt)
!2222222222222222222222222222222

! for an individual pixel 

  use m_vars, only : vcd_convfac

  implicit none

  real(kind=4),dimension(geos_np)::tt
  real(kind=4),dimension(geos_np+1)::pp !include Psfc
  real(kind=4),dimension(geos_np+1)::tmp_vcd
  real::sum_vcd!,psfc
  real::xx1,xx2,yy1,yy2,xxx,yyy
  !real::x1,x2,x3,y1,y2,y3
  integer(kind=4)::iflag,ipcld!,iline
  integer(kind=4)::ip!,ix,iy

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

!5555555555555555555555555555
subroutine read_BRDF_Rsfc_h5
!5555555555555555555555555555

  IMPLICIT NONE

  CHARACTER(LEN=255)::filename
  CHARACTER(LEN=255)::dsetname

! HID_T type integers.
  INTEGER(HID_T)::file_id,dset_id,datatype_id

! Regular four-byte integer.
  INTEGER(KIND=4)::hdf_err,ierr
  real(kind=4),dimension(:,:,:),allocatable::BRDF_SurfaceReflectivity

! HSIZE_T type integer.
  !INTEGER(HSIZE_T),DIMENSION(1)::dims1
  !INTEGER(HSIZE_T),DIMENSION(2)::dims2
  INTEGER(HSIZE_T),DIMENSION(3)::dims3

  integer(kind=4)::nt,nx,nw,it,ix

!-------
! start
!-------
  filename=trim(name_brdf_dir)//trim(name_brdf_file)

  !hqw replace
  !nt=inp_NumTimes
  !nx=inp_nXtrack
  nt = rad_NumTimes
  nx = rad_nXtrack
  nw=2
  
! allocate dimensions & fill values for outputs
!hqw BRDF_SurfaceReflectivity is local, will be deallocated at end
  allocate(BRDF_SurfaceReflectivity(nw,nx,nt),stat=ierr)
  BRDF_SurfaceReflectivity=-9.9
  allocate(BRDF_SurfaceReflectivity466(nx,nt),stat=ierr)
  BRDF_SurfaceReflectivity466=-9.9
  allocate(BRDF_SurfaceReflectivity440(nx,nt),stat=ierr)
  BRDF_SurfaceReflectivity440=-9.9

! Initialize FORTRAN interface.
  CALL h5open_f(hdf_err)

! Open an existing file.
  CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_err)
if (hdf_err /= 0) print *, "read_BRDF_Rsfc_h5 fail"

! Open an existing dataset.
  dsetname="/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/GLER"
  CALL h5dopen_f(file_id,dsetname,dset_id,hdf_err)
  CALL h5dget_type_f(dset_id,datatype_id,hdf_err)
  dims3=SHAPE(BRDF_SurfaceReflectivity)
  CALL H5Dread_f(dset_id,datatype_id,BRDF_SurfaceReflectivity,dims3,hdf_err)
  CALL h5dclose_f(dset_id,hdf_err)
    
! Close the file.
  CALL h5fclose_f(file_id,hdf_err)

! Close FORTRAN interface.
  CALL h5close_f(hdf_err)

  do it=1,nt
  do ix=1,nx
    BRDF_SurfaceReflectivity466(ix,it)=BRDF_SurfaceReflectivity(2,ix,it)
    BRDF_SurfaceReflectivity440(ix,it)=BRDF_SurfaceReflectivity(1,ix,it)
  end do
  end do

  !hqw deallocate helpers
  deallocate(BRDF_SurfaceReflectivity)
 
!55555555555555555555555555555555
end subroutine read_BRDF_Rsfc_h5
!55555555555555555555555555555555

!hqw GLER now handled by m_read_input_gler.f90
!the following is no longer used 
!6666666666666666666666666666
subroutine read_BDEM_Psfc_h5
!6666666666666666666666666666

  IMPLICIT NONE

  CHARACTER(LEN=255)::filename
  CHARACTER(LEN=255)::dsetname

! HID_T type integers.
  INTEGER(HID_T)::file_id,dset_id,datatype_id

! Regular four-byte integer.
  INTEGER(KIND=4)::hdf_err,ierr

! HSIZE_T type integer.
  !INTEGER(HSIZE_T),DIMENSION(1)::dims1
  INTEGER(HSIZE_T),DIMENSION(2)::dims2
  !INTEGER(HSIZE_T),DIMENSION(3)::dims3

  integer(kind=4)::nt,nx!,it,ix

!-------
! start
!-------
  filename=trim(name_brdf_dir)//trim(name_brdf_file)

  !hqw replace nt, nx as inp_ were from OMCLDO2
  !nt=inp_NumTimes
  !nx=inp_nXtrack
  nt = rad_NumTimes
  nx = rad_nXtrack
  
! allocate dimensions & fill values for outputs
! hqw STDev & LandAreaFraction are not actually used in calculation,
! thus, comment out to save memory
  allocate(BDEM_TerrainPressure(nx,nt),stat=ierr)
!  allocate(BDEM_TerrainPressureStdDev(nx,nt),stat=ierr)
  allocate(BDEM_TerrainHeight(nx,nt),stat=ierr)
!  allocate(BDEM_TerrainHeightStdDev(nx,nt),stat=ierr)
!  allocate(BDEM_LandAreaFraction(nx,nt),stat=ierr)

! Initialize FORTRAN interface.
  CALL h5open_f(hdf_err)

! Open an existing file.
  CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_err)
if (hdf_err /= 0) print *, "read_BDEM_Psfc_h5 fail"

! Open an existing dataset.
  dsetname="/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/TerrainPressure"
  CALL h5dopen_f(file_id,dsetname,dset_id,hdf_err)
  CALL h5dget_type_f(dset_id,datatype_id,hdf_err)
  dims2=SHAPE(BDEM_TerrainPressure)
  CALL H5Dread_f(dset_id,datatype_id,BDEM_TerrainPressure,dims2,hdf_err)
  CALL h5dclose_f(dset_id,hdf_err)

!  dsetname="/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/TerrainPressureStdDev"
!  CALL h5dopen_f(file_id,dsetname,dset_id,hdf_err)
!  CALL h5dget_type_f(dset_id,datatype_id,hdf_err)
!  dims2=SHAPE(BDEM_TerrainPressureStdDev)
!  CALL H5Dread_f(dset_id,datatype_id,BDEM_TerrainPressureStdDev,dims2,hdf_err)
!  CALL h5dclose_f(dset_id,hdf_err)

  dsetname="/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/TerrainHeight"
  CALL h5dopen_f(file_id,dsetname,dset_id,hdf_err)
  CALL h5dget_type_f(dset_id,datatype_id,hdf_err)
  dims2=SHAPE(BDEM_TerrainHeight)
  CALL H5Dread_f(dset_id,datatype_id,BDEM_TerrainHeight,dims2,hdf_err)
  CALL h5dclose_f(dset_id,hdf_err)

!  dsetname="/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/TerrainHeightStdDev"
!  CALL h5dopen_f(file_id,dsetname,dset_id,hdf_err)
!  CALL h5dget_type_f(dset_id,datatype_id,hdf_err)
!  dims2=SHAPE(BDEM_TerrainHeightStdDev)
!  CALL H5Dread_f(dset_id,datatype_id,BDEM_TerrainHeightStdDev,dims2,hdf_err)
!  CALL h5dclose_f(dset_id,hdf_err)

!  dsetname="/HDFEOS/SWATHS/Geometry Dependent Surface LER/Data Fields/LandAreaFraction"
!  CALL h5dopen_f(file_id,dsetname,dset_id,hdf_err)
!  CALL h5dget_type_f(dset_id,datatype_id,hdf_err)
!  dims2=SHAPE(BDEM_LandAreaFraction)
!  CALL H5Dread_f(dset_id,datatype_id,BDEM_LandAreaFraction,dims2,hdf_err)
!  CALL h5dclose_f(dset_id,hdf_err)

! Close the file.
  CALL h5fclose_f(file_id,hdf_err)

! Close FORTRAN interface.
  CALL h5close_f(hdf_err)

!66666666666666666666666666666666
end subroutine read_BDEM_Psfc_h5
!66666666666666666666666666666666

!**********************
end module m_read_hdf5
!**********************
