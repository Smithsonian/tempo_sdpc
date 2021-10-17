!*****************
module m_read_lut
  !*****************

  use tell_module

contains

  !11111111111111111111111
  subroutine read_lut_rad (errstat)
    !11111111111111111111111

    use hdfeos4_parameters
    use he5_swreader
    use hdf5
    use H5LT
    use H5DS
    use m_vars

    implicit none

    !---------------------------------------------------------------------72
    ! ROUTINE: m_read_lut
    ! 
    ! DESCRIPTION: This program reads lookup table for radiance and AMF
    !
    ! REVISION HISTORY: 
    !
    !  04/23/15 Yang original fortran 90
    !---------------------------------------------------------------------72

    !  real(kind=4),dimension(:,:,:,:,:),pointer:lut_rad_clr
    !  real(kind=4),dimension(:,:,:,:,:),pointer:lut_amf_clr
    !  real(kind=4),dimension(:,:,:,:,:),pointer:lut_amf_cld
    !  real(kind=4),dimension(:),pointer:lut_sza
    !  real(kind=4),dimension(:),pointer:lut_vza
    !  real(kind=4),dimension(:),pointer:lut_raa
    !  real(kind=4),dimension(:),pointer:lut_alb
    !  real(kind=4),dimension(:),pointer:lut_psfc
    !  real(kind=4),dimension(:),pointer:lut_pcld

    integer (kind=4), intent(inout) :: errstat

    character(len=255)::filename
    !integer::ndim,ntmp
    integer(kind=4)::ierr!,status
    integer(hid_t)::h5fid,h5id
    integer(hsize_t),dimension(1)::dims1
    integer(hsize_t),dimension(5)::dims5

    if (errstat /= 0) return

    !-------------------
    ! read lookup table
    !-------------------
    filename=trim(name_lut_dir)//'/'//trim(name_lut_rad)
    write(*,*) '   reading '//trim(filename)

    call h5open_f(ierr)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,h5fid,ierr)
    if (ierr .lt. 0) then
      call tell_error (tell_io_read_error, "read_lut_rad: failed open", &
           errstat)
      return
    endif


    allocate(lut_alb(1:nalb), lut_sza(1:nsza), &
         lut_vza(1:nvza), lut_raa(1:nraa), &
         lut_psfc(1:npsfc),lut_rad_clr(1:nalb,1:nsza,1:nvza,1:nraa,1:npsfc), &
         stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_lut_rad: allocate fail", &
           errstat)
    endif
    !----------------------
    ! ALB,SZA,VZA,RAA,PSFC
    !----------------------

    dims1(1)=nalb
    call h5dopen_f(h5fid,"ALB",h5id,ierr)
    call h5dread_f(h5id,H5T_NATIVE_REAL,lut_alb,dims1,ierr)
    call h5dclose_f(h5id,ierr)

    if (ierr .ge. 0) then
      dims1(1)=nsza
      call h5dopen_f(h5fid,"SZA",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_sza,dims1,ierr)
      call h5dclose_f(h5id,ierr)
    endif

    if (ierr .ge. 0) then
      dims1(1)=nvza
      call h5dopen_f(h5fid,"VZA",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_vza,dims1,ierr)
      call h5dclose_f(h5id,ierr)
    endif

    if (ierr .ge. 0) then
      dims1(1)=nraa
      call h5dopen_f(h5fid,"RAA",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_raa,dims1,ierr)
      call h5dclose_f(h5id,ierr)
    endif

    if (ierr .ge. 0) then
      dims1(1)=npsfc
      call h5dopen_f(h5fid,"Psfc",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_psfc,dims1,ierr)
      call h5dclose_f(h5id,ierr)
    endif
    !---------------------
    ! Radiance: clear sky
    !---------------------
    if (ierr .ge. 0) then
      dims5(1)=nalb
      dims5(2)=nsza
      dims5(3)=nvza
      dims5(4)=nraa
      dims5(5)=npsfc
      call h5dopen_f(h5fid,"RAD",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_rad_clr,dims5,ierr)
      call h5dclose_f(h5id,ierr)
    endif

    if (ierr .lt. 0) then
      call tell_error (tell_io_read_error, &
           "read_lut_rad: failed to read data", errstat)
      return
    endif

    !------------------------------------------------------------
    ! Radiance: cloud --> lut_rad_clr(18,1:nsza,1:nvza,1:nraa,18)
    !------------------------------------------------------------

    !-------
    ! close
    !-------
    call h5fclose_f(h5fid,ierr)
    call h5close_f(ierr)

  !===========================
  end subroutine read_lut_rad
  !===========================

  !440440440440440440440440440
  subroutine read_lut_rad440 (errstat)
    !440440440440440440440440440

    use hdfeos4_parameters
    use he5_swreader
    use hdf5
    use H5LT
    use H5DS
    use m_vars

    implicit none

    integer (kind=4), intent(inout) :: errstat

    character(len=255)::filename
    !integer::ndim,ntmp
    integer(kind=4)::ierr!,status
    integer(hid_t)::h5fid,h5id
    !integer(hsize_t),dimension(1)::dims1
    integer(hsize_t),dimension(5)::dims5


    if (errstat /= 0) return

    !-------------------
    ! read lookup table
    !-------------------
    filename=trim(name_lut_dir)//'/'//trim(name_lut_rad440)
    write(*,*)'   reading '//trim(filename)

    call h5open_f(ierr)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,h5fid,ierr)
    if (ierr .lt. 0) then
      call tell_error (tell_io_open_error, "read_lut_rad440: failed open", &
           errstat)
      return
    endif

    !----------------------
    ! ALB,SZA,VZA,RAA,PSFC
    !----------------------
    ! EOS - these are all identical to those already read in read_lut_rad!!!
    !       SKIPPING
    !dims1(1)=nalb
    !allocate(lut_alb(1:nalb),stat=ierr)
    !call h5dopen_f(h5fid,"ALB",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_alb,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nsza
    !allocate(lut_sza(1:nsza),stat=ierr)
    !call h5dopen_f(h5fid,"SZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_sza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nvza
    !allocate(lut_vza(1:nvza),stat=ierr)
    !call h5dopen_f(h5fid,"VZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_vza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nraa
    !allocate(lut_raa(1:nraa),stat=ierr)
    !call h5dopen_f(h5fid,"RAA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_raa,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=npsfc
    !allocate(lut_psfc(1:npsfc),stat=ierr)
    !call h5dopen_f(h5fid,"Psfc",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_psfc,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !---------------------
    ! Radiance: clear sky
    !---------------------
    allocate(lut_rad_clr440(1:nalb,1:nsza,1:nvza,1:nraa,1:npsfc),stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_lut_rad440: allocate fail", &
           errstat)
      return
    endif
    dims5(1)=nalb
    dims5(2)=nsza
    dims5(3)=nvza
    dims5(4)=nraa
    dims5(5)=npsfc
    call h5dopen_f(h5fid,"RAD",h5id,ierr)
    call h5dread_f(h5id,H5T_NATIVE_REAL,lut_rad_clr440,dims5,ierr)
    call h5dclose_f(h5id,ierr)
    if (ierr .lt. 0) then
      call tell_error (tell_io_read_error, &
           "read_lut_rad440: failed to read data", errstat)
      return
    endif

    !------------------------------------------------------------
    ! Radiance: cloud --> lut_rad_clr(18,1:nsza,1:nvza,1:nraa,18)
    !------------------------------------------------------------

    !-------
    ! close
    !-------
    call h5fclose_f(h5fid,ierr)
    call h5close_f(ierr)

    !440440440440440440440440440440
  end subroutine read_lut_rad440
  !440440440440440440440440440440


  !222222222222222222222222222
  subroutine read_lut_amf_clr (errstat)
    !222222222222222222222222222

    use hdfeos4_parameters
    use he5_swreader
    use hdf5
    use H5LT
    use H5DS
    use m_vars

    implicit none

    integer (kind=4), intent(inout) :: errstat

    character(len=255)::filename
    !integer::ndim,ntmp
    integer(kind=4)::ierr!,status
    integer(hid_t)::h5fid,h5id
    !integer(hsize_t), dimension(1) :: dims1
    integer(hsize_t), dimension(5) :: dims5

    !-------------------
    ! read lookup table
    !-------------------
    filename=trim(name_lut_dir)//'/'//trim(name_lut_amf_clr)
    write(*,*)'   reading '//trim(filename)

    call h5open_f(ierr)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,h5fid,ierr)
    if (ierr .lt. 0) then
      call tell_error (tell_io_open_error, "read_lut_amf_clr: failed open", &
           errstat)
      return
    endif

    !----------------------
    ! ALB,SZA,VZA,RAA,PSFC
    !----------------------
    ! EOS - Once again, all these have already been read in and are identical
    !       in file, SKIP
    !dims1(1)=nalb
    !allocate(lut_alb(1:nalb),stat=ierr)
    !call h5dopen_f(h5fid,"ALB",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_alb,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nsza
    !allocate(lut_sza(1:nsza),stat=ierr)
    !call h5dopen_f(h5fid,"SZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_sza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nvza
    !allocate(lut_vza(1:nvza),stat=ierr)
    !call h5dopen_f(h5fid,"VZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_vza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nraa
    !allocate(lut_raa(1:nraa),stat=ierr)
    !call h5dopen_f(h5fid,"RAA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_raa,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=npsfc
    !allocate(lut_psfc(1:npsfc),stat=ierr)
    !call h5dopen_f(h5fid,"Psfc",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_psfc,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !-----------------
    ! Air Mass Factor
    !-----------------
    dims5(1)=nalb
    dims5(2)=nsza
    dims5(3)=nvza
    dims5(4)=nraa
    dims5(5)=npsfc
    allocate(lut_amf_clr(1:nalb,1:nsza,1:nvza,1:nraa,1:npsfc),&
         lut_rad_ler(1:nalb,1:nsza,1:nvza,1:nraa,1:npsfc),stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_lut_amf_clr: allocate fail", &
           errstat)
      return
    endif

    call h5dopen_f(h5fid,"AMF",h5id,ierr)
    call h5dread_f(h5id,H5T_NATIVE_REAL,lut_amf_clr,dims5,ierr)
    call h5dclose_f(h5id,ierr)

    if (ierr .ge. 0) then
      call h5dopen_f(h5fid,"RAD",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_rad_ler,dims5,ierr)
      call h5dclose_f(h5id,ierr)
    endif

    if (ierr .lt. 0) then
      call tell_error (tell_io_read_error, &
           "read_lut_amf_clr: failed to read data", errstat)
      return
    endif

    !-------
    ! close
    !-------
    call h5fclose_f(h5fid,ierr)
    call h5close_f(ierr)

    !===============================
  end subroutine read_lut_amf_clr
  !===============================


  !333333333333333333333333333
  subroutine read_lut_amf_cld (errstat)
    !333333333333333333333333333

    use hdfeos4_parameters
    use he5_swreader
    use hdf5
    use H5LT
    use H5DS
    use m_vars

    implicit none

    integer(kind=4), intent(inout) :: errstat

    character(len=255)::filename
    !integer::ndim,ntmp
    integer(kind=4)::ierr!,status
    integer(hid_t)::h5fid,h5id
    integer(hsize_t), dimension(1) :: dims1
    integer(hsize_t), dimension(5) :: dims5

    !-------------------
    ! read lookup table
    !-------------------
    filename=trim(name_lut_dir)//'/'//trim(name_lut_amf_cld)
    !filename=trim(name_lut_amf_cld)
    write(*,*)'   reading '//trim(filename)

    call h5open_f(ierr)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,h5fid,ierr)
    if (ierr .lt. 0) then
      call tell_error (tell_io_open_error, "read_lut_amf_cld: failed open", &
           errstat)
      return
    endif

    !---------------------------
    ! SZA,VZA,RAA,PSFC,PCLD
    !---------------------------
    !EOS - SKIP the first four which have already been read
    !dims1(1)=nsza
    !allocate(lut_sza(1:nsza),stat=ierr)
    !call h5dopen_f(h5fid,"SZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_sza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nvza
    !allocate(lut_vza(1:nvza),stat=ierr)
    !call h5dopen_f(h5fid,"VZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_vza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nraa
    !allocate(lut_raa(1:nraa),stat=ierr)
    !call h5dopen_f(h5fid,"RAA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_raa,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=npsfc
    !allocate(lut_psfc(1:npsfc),stat=ierr)
    !call h5dopen_f(h5fid,"Psfc",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_psfc,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    allocate(lut_pcld(1:npcld), &
         lut_amf_cld(1:nsza,1:nvza,1:nraa,1:npcld,1:npsfc),stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_lut_amf_cld: allocate error", &
           errstat)
      return
    endif

    dims1(1)=npcld
    call h5dopen_f(h5fid,"Pcld",h5id,ierr)
    call h5dread_f(h5id,H5T_NATIVE_REAL,lut_pcld,dims1,ierr)
    call h5dclose_f(h5id,ierr)

    !-----------------
    ! Air Mass Factor
    !-----------------
    dims5(1)=nsza
    dims5(2)=nvza
    dims5(3)=nraa
    dims5(4)=npcld
    dims5(5)=npsfc
    if (ierr .ge. 0) then
      call h5dopen_f(h5fid,"AMF",h5id,ierr)
      call h5dread_f(h5id,H5T_NATIVE_REAL,lut_amf_cld,dims5,ierr)
      call h5dclose_f(h5id,ierr)
    endif

    if (ierr .lt. 0) then
      call tell_error (tell_io_read_error, &
           "read_lut_amf_cld: failed to read data", errstat)
      return
    endif

    !-------
    ! close
    !-------
    call h5fclose_f(h5fid,ierr)
    call h5close_f(ierr)

    !===============================
  end subroutine read_lut_amf_cld
  !===============================

  !444444444444444444444444444
  subroutine read_lut_amf_ler (errstat)
    !444444444444444444444444444

    use hdfeos4_parameters
    use he5_swreader
    use hdf5
    use H5LT
    use H5DS
    use m_vars

    implicit none

    integer (kind=4), intent(inout) :: errstat

    character(len=255)::filename
    !integer::ndim,ntmp
    integer(kind=4)::ierr!,status
    integer(hid_t)::h5fid,h5id
    !integer(hsize_t),dimension(1)::dims1
    integer(hsize_t),dimension(6)::dims6


    if (errstat /= 0) return

    !-------------------
    ! read lookup table
    !-------------------
    filename=trim(name_lut_dir)//'/'//trim(name_lut_ler)
    !filename=trim(name_lut_ler)
    write(*,*)'   reading '//trim(filename)

    call h5open_f(ierr)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,h5fid,ierr)
    if (ierr .lt. 0) then
      call tell_error (tell_io_open_error, "read_lut_amf_ler: failed open", &
           errstat)
      return
    endif

    !----------------------
    ! ALB,SZA,VZA,RAA,PSFC
    !----------------------
    ! EOS - SKIPPING the first few again, already read in and identical 
    !       between files
    !dims1(1)=nalb
    !allocate(lut_alb(1:nalb),stat=ierr)
    !call h5dopen_f(h5fid,"ALB",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_alb,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nsza
    !allocate(lut_sza(1:nsza),stat=ierr)
    !call h5dopen_f(h5fid,"SZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_sza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nvza
    !allocate(lut_vza(1:nvza),stat=ierr)
    !call h5dopen_f(h5fid,"VZA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_vza,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=nraa
    !allocate(lut_raa(1:nraa),stat=ierr)
    !call h5dopen_f(h5fid,"RAA",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_raa,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=npsfc
    !allocate(lut_psfc(1:npsfc),stat=ierr)
    !call h5dopen_f(h5fid,"Psfc",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_psfc,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !dims1(1)=npcld
    !allocate(lut_pcld(1:npcld),stat=ierr)
    !call h5dopen_f(h5fid,"Pcld",h5id,ierr)
    !call h5dread_f(h5id,H5T_NATIVE_REAL,lut_pcld,dims1,ierr)
    !call h5dclose_f(h5id,ierr)

    !------------------
    ! AMF and radiance
    !------------------
    dims6(1)=nalb
    dims6(2)=nsza
    dims6(3)=nvza
    dims6(4)=nraa
    dims6(5)=npcld
    dims6(6)=npsfc

    allocate(lut_amf_ler(1:nalb,1:nsza,1:nvza,1:nraa,1:npcld,1:npsfc),&
         stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_rad_amf_ler: allocate fail", &
           errstat)
      return
    endif

    call h5dopen_f(h5fid,"AMF",h5id,ierr)
    call h5dread_f(h5id,H5T_NATIVE_REAL,lut_amf_ler,dims6,ierr)
    call h5dclose_f(h5id,ierr)

    if (ierr .lt. 0) then
      call tell_error ( tell_io_read_error, &
           "read_rad_amf_ler: failed to read data", errstat)
      return
    endif
    !-------
    ! close
    !-------
    call h5fclose_f(h5fid,ierr)
    call h5close_f(ierr)

    !===============================
  end subroutine read_lut_amf_ler
  !===============================

  !*********************
end module m_read_lut
!*********************
