module irradiance_data

  use OMSAO_precision_module, only: i2, i4, r4, r8
  use tell_module
  USE sao_pge_utils, ONLY: print_array
  implicit none

  type, public :: Irradiance_Data_Type
    integer (kind=i4) :: max_nwave, nxtrack
    ! Each 2d array A has size A[max_nwave, nxtrack], where only
    ! A[1:nwaves[j], j], j=1:nxtrack are meaningful.
    ! The 1d arrays have size [nxtrack].
    integer(kind=i2), dimension(:,:), allocatable :: qflags
    real (kind=r8), dimension(:,:), allocatable :: &
      wavelengths, spectrum
    integer (kind=i4), dimension(:), allocatable :: nwaves
    real (kind=r8), dimension(:), allocatable :: avg_wavelengths
    integer (kind=i4), dimension(:,:), allocatable :: ccdpix_selection
    integer (kind=i4), dimension(:,:), allocatable :: ccdpix_exclusion
  end type

  type (Irradiance_Data_Type) :: Irr_Data

  private
  public irradiance_data_init, Irr_Data, deallocate_irr_data_type

contains

  subroutine irradiance_data_init (source, errstat)

    use OMSAO_variables_module, only: OMSAO_I0_filename, l1b_irrad_filename
    use netcdf, only : nf90_nowrite
    use tio_module

    implicit none
    character (len=*), intent (in) :: source
    integer, intent(inout) :: errstat
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    ! --------------------------------------------------------------------
    ! The definition of TEMPO's epoch is always needed. Read it from
    ! solar irradiance file prior to start the solar calibration
    ! --------------------------------------------------------------------
    call tell_log (1, 'reading epoch from '//trim(l1b_irrad_filename))
    call tiof_open (l1b_irrad_filename, tio_l1obj, nf90_nowrite, errstat)
    call tiof_use_file_epoch (tio_l1obj, errstat)
    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, &
           "irradiance_data_init:  failed reading epoch",  errstat)
      return
    endif
    
    select case ( source )
      case ('I0_irradiance')
        call read_I0_irradiance(OMSAO_I0_filename, errstat)
      case ('solar_irradiance')
        call read_irradiance_data (errstat)
      case default
        call tell_error (tell_runtime_error, &
            "irradiance_data_init: not a valid source: "//source, errstat)
        return
    end select

  end subroutine irradiance_data_init

  ! =========================================================================

  subroutine allocate_irr_data_type (idt, max_nwave, nxtrack, errstat)

    implicit none
    type (Irradiance_Data_Type), intent(inout) :: idt
    integer (kind=i4), intent(in) :: max_nwave, nxtrack
    integer, intent(inout) :: errstat

    integer :: locerr

    if (errstat /= 0) return

    allocate (idt%qflags(max_nwave, nxtrack), &
              idt%wavelengths(max_nwave, nxtrack), &
              idt%spectrum(max_nwave, nxtrack), &
              idt%nwaves(nxtrack), &
              idt%avg_wavelengths(nxtrack), &
              idt%ccdpix_selection(4, nxtrack), &
              idt%ccdpix_exclusion(2, nxtrack), &
              stat=locerr)
    if (locerr /= 0) then
      call tell_error (tell_malloc_error, &
                       "allocate_irr_data_type: allocate failed", &
                       errstat)
      return
    endif

    idt%max_nwave = max_nwave
    idt%nxtrack = nxtrack

  end subroutine allocate_irr_data_type

  ! =========================================================================

  subroutine deallocate_irr_data_type (idt, errstat)

    implicit none

    type (Irradiance_Data_Type), intent(inout) :: idt
    integer, intent(inout) :: errstat

    if (allocated(idt%qflags)) deallocate (idt%qflags, stat=errstat)
    if (allocated(idt%wavelengths)) deallocate (idt%wavelengths, stat=errstat)
    if (allocated(idt%spectrum)) deallocate (idt%spectrum, stat=errstat)
    if (allocated(idt%nwaves)) deallocate (idt%nwaves, stat=errstat)
    if (allocated(idt%avg_wavelengths)) deallocate (idt%avg_wavelengths, &
         stat=errstat)
    if (allocated(idt%ccdpix_selection)) deallocate (idt%ccdpix_selection, &
         stat=errstat)
    if (allocated(idt%ccdpix_exclusion)) deallocate (idt%ccdpix_exclusion, &
         stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "deallocate_irr_data_type: failed", &
           errstat)
      return
    endif


  end subroutine deallocate_irr_data_type

  ! =========================================================================

  subroutine package_irradiance_data (nwl, nxtrack, &
                                      wavelengths, spectrum, qflags, &
                                      errstat)
    ! In this routine, only the first nwl wavelengths are meaningful.
    use arrayutils, only: array_find_inner_bounding_indices_r8, &
      array_find_bounding_indices_r8
    use OMSAO_variables_module, only: &
      ctrl_fit_winexc_lim, ctrl_fit_winwav_lim

    implicit none
    integer (kind=i4), intent(in) :: nwl, nxtrack
    real (kind=r8), dimension(:, :), intent(in) :: wavelengths, spectrum
    integer (kind=i2), dimension(:,:), intent(in) :: qflags
    integer, intent(inout) :: errstat

    integer (kind=i4), dimension (:,:), allocatable :: ccdpix_sel
    integer (kind=i4) :: ix, icnt, imin, imax, max_nwavel
    integer :: locerrstat

    if (errstat /= 0) return

    allocate (ccdpix_sel (4, nxtrack), stat = locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "allocate_irr_data_type: allocate failed", errstat)
      return
    endif
    ccdpix_sel = -1

    ! ----------------------------------------------------
    ! Limit irradiance arrays to fitting window. Check for
    ! strictly ascending wavelengths in the process.
    ! ----------------------------------------------------
    DO ix = 1, nxtrack

      ! Restrict wavelengths to selected wavelength fitting window
      call array_find_bounding_indices_r8 (nwl, wavelengths(1:nwl, ix), &
                                           ctrl_fit_winwav_lim(1), &
                                           ctrl_fit_winwav_lim(2), &
                                           ccdpix_sel (1,ix), ccdpix_sel(2,ix))

      call array_find_bounding_indices_r8 (nwl, wavelengths(1:nwl, ix), &
                                           ctrl_fit_winwav_lim(3), &
                                           ctrl_fit_winwav_lim(4), &
                                           ccdpix_sel (3,ix), ccdpix_sel(4,ix))
    end do

    max_nwavel = 1 + MAXVAL (ccdpix_sel(4, :) - ccdpix_sel(1,:))

    call allocate_irr_data_type (Irr_Data, max_nwavel, nxtrack, errstat)
    if (errstat /= 0) return
    
    Irr_Data%ccdpix_selection = ccdpix_sel

    do ix=1, nxtrack

      imin = ccdpix_sel(1,ix)
      imax = ccdpix_sel(4,ix)
      icnt = imax - imin + 1

      Irr_Data%wavelengths(1:icnt,ix) = wavelengths(imin:imax,ix)

      Irr_Data%spectrum(1:icnt, ix) = spectrum(imin:imax,ix)

      Irr_Data%qflags(1:icnt, ix) = qflags(imin:imax,ix)

      Irr_Data%avg_wavelengths(ix) &
        = SUM(wavelengths(imin:imax, ix))/REAL(icnt, kind=r8)

      Irr_Data%nwaves(ix) = icnt

      ! ------------------------------------------------------------------------------
      ! If any window is excluded, find the corresponding indices. This has to be done
      ! after the array assignements above because we need to know which indices to
      ! exclude from the final arrays, not the complete ones read from the HE4 file.
      ! ------------------------------------------------------------------------------
      Irr_Data%ccdpix_exclusion(1:2, ix) = -1

      if ( minval(ctrl_fit_winexc_lim(1:2)) > 0.0_r8 ) then

        call array_find_inner_bounding_indices_r8 ( &
          nwl, wavelengths(1:nwl, ix), &
          ctrl_fit_winexc_lim(1), &
          ctrl_fit_winexc_lim(2), &
          Irr_Data%ccdpix_exclusion(1, ix), &
          Irr_Data%ccdpix_exclusion(2, ix))
      end if

    end do

  end subroutine

  ! =========================================================================

  SUBROUTINE read_irradiance_data (errstat)

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  ONLY: &
      l1b_irrad_filename, l1b_channel
    USE arrayutils, only: array_locate_r4
    use ctrlvars, only: yn_disable_omi_features, yn_gems
    use m_read_gems, only: gems_read_irrad_data

    use l1bread_utils
    use tio_module
    use netcdf, only : nf90_nowrite

    implicit none
    integer (kind=i4), intent (inout) :: errstat
    !
    integer :: locerrstat
    integer (kind=i4) :: nwavel, ix, nxtrack
    real (kind=r4), dimension(:,:,:), allocatable :: &
      tmp_wavelengths, tmp_spectrum
    real (kind=r8), dimension(:,:), allocatable :: &
      wavelengths, spectrum
    integer (kind=i2), dimension (:,:,:), allocatable :: tmp_qflags
    integer (kind=i2), dimension (:,:), allocatable :: tmp_qflags_2d
    character (len=64) :: swathname
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    ! Allow errstat to flow
    if (yn_gems) then
      call gems_read_irrad_data (nwavel, nxtrack, wavelengths, spectrum, &
           tmp_qflags_2d, errstat)
    else !TEMPO
      call tell_log (1, 'reading irradiances = '//trim(l1b_irrad_filename))
      call tiof_open (l1b_irrad_filename, tio_l1obj, nf90_nowrite, errstat)
      call lookup_swathname (l1b_channel, swathname, errstat)
      call tiof_inq_group (tio_l1obj, swathname, errstat)
      call tiof_inq_dimlen (tio_l1obj, "xtrack", nxtrack, errstat)
      call tiof_inq_dimlen (tio_l1obj, "spectral_channel", nwavel, errstat)
      if (errstat /= 0) return

      allocate (tmp_wavelengths(nwavel, nxtrack,1), &
              tmp_spectrum(nwavel, nxtrack,1), &
              tmp_qflags (nwavel, nxtrack,1), &
              wavelengths (nwavel, nxtrack), &
              spectrum (nwavel, nxtrack), &
              stat=locerrstat)
      if (locerrstat /= 0) then
        call tell_error (tell_malloc_error, &
             "read_irradiance_data: allocate failed", errstat)
        return
      endif

      call tiof_get3d_r4 (tio_l1obj, "irradiance", [0,0,0], &
           [1,nxtrack,nwavel], tmp_spectrum(:,1:nxtrack,1:1), errstat)
      call tiof_get3d_i2 (tio_l1obj, "pixel_quality_flag", [0,0,0], &
           [1,nxtrack,-1], tmp_qflags(:,1:nxtrack,1:1), errstat)
      call tiof_get3d_r4 (tio_l1obj, "wavelength", [0,0,0], &
           [1,nxtrack,nwavel], tmp_wavelengths(:,1:nxtrack,1:1), errstat)
      call tiof_close (tio_l1obj, errstat)
      if (errstat /= 0) then
        call tell_error (tell_runtime_error, &
             "read_irradiance_data:  failed reading irradiance data",  errstat)
        return
      endif

      ! -------------------------------
      ! Reverse arrays for UV-1 channel.
      !   Why?  Are they stored in descending order?  --JED
      ! -------------------------------
      if (.not.yn_disable_omi_features) then
        IF ( l1b_channel == 'UV1' ) THEN
          DO ix = 1, nxtrack
            tmp_wavelengths(nwavel:1:-1, ix,1) = tmp_wavelengths(1:nwavel,ix,1)
            tmp_spectrum(nwavel:1:-1, ix,1) = tmp_spectrum(1:nwavel,ix,1)
            tmp_qflags(nwavel:1:-1, ix,1) = tmp_qflags(1:nwavel,ix,1)
          END DO
        END IF
      endif

      wavelengths = real (tmp_wavelengths(:,:,1), kind=r8)
      spectrum = real (tmp_spectrum(:,:,1), kind=r8)
      deallocate (tmp_wavelengths, stat=errstat)
      if (errstat == 0) deallocate (tmp_spectrum, stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "read_irradiance_data: deallocate failed", errstat)
        return
      endif
      tmp_qflags_2d = tmp_qflags(:,:,1)
    endif ! TEMPO/GEMS

    call package_irradiance_data (nwavel, nxtrack, &
                                  wavelengths, spectrum, tmp_qflags_2d, &
                                  errstat)

    return
  end subroutine read_irradiance_data

  !--------------------------------------------------------------------
  !> Read I0 irradiance-replacement spectra from file designated in PCF
  !--------------------------------------------------------------------
  !
  ! @param[in]  filename   Name of I0 spectra file
  ! @param      errstat    Error tracking integer, non-zero = problem
  !
  ! @author  E. O'Sullivan   November 2020
  !-------------------------------------------------------------------------
  subroutine read_I0_irradiance (filename, errstat)

    use ctrlvars, only: yn_spectrum_norm
    USE OMSAO_parameters_module, ONLY: downweight, normweight
    USE OMSAO_omidata_module, ONLY: omi_irradiance_wght
    implicit none

    !input variables
    character (len=*), intent(in) :: filename
    integer (kind=4), intent(inout) :: errstat
    !local variables
    integer (kind=4) :: nxtrack, nwavel, i
    real (kind=8), dimension(:,:), allocatable :: spectra, wavelengths
    real (kind=8), dimension(:), allocatable :: weightsum
    integer (kind=2), dimension(:,:), allocatable :: qflags
    integer (kind=2), parameter :: bad_pixel = 1

    if (errstat /= 0) return

    call tell_log (1, 'reading I0 irradiances = '//trim(filename))

    call read_I0_dims (filename, nxtrack, nwavel, errstat)
    if (errstat /= 0) return

    allocate ( spectra(nwavel,nxtrack), wavelengths(nwavel, nxtrack), &
         qflags(nwavel, nxtrack), weightsum(nwavel), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_I0_irradiance: allocation error", errstat)
      return
    endif

    call read_I0_spectra (filename, nxtrack, nwavel, spectra, wavelengths, &
         qflags, errstat)
    if (errstat /= 0) return

    call package_irradiance_data (nwavel, nxtrack, wavelengths, spectra, &
         qflags, errstat)

    ! Set omi_irradiance_wght values. This is important if reading
    ! slift calibration parameters since the spectra will not be fitted
    ! and bad pixels identified. If the I0_irradiance is used for
    ! slift calibration then omi_irradiance_wght is set later on
    do i = 1, Irr_Data%nxtrack
      nwavel = Irr_Data%nwaves(i)
      where (Irr_Data%qflags(1:nwavel,i) == bad_pixel)
        omi_irradiance_wght(1:nwavel,i) = downweight
      endwhere
    enddo

    ! It is important to normalize the I0_irradiance 
    if (yn_spectrum_norm) then
      do i = 1, Irr_Data%nxtrack
        weightsum = 1.0
        nwavel = Irr_Data%nwaves(i)
        where (Irr_Data%qflags(1:nwavel,i) == bad_pixel)
          weightsum = 0.0
          ! Irr_Data%spectrum(1:nwavel,i) = 0.0
        endwhere
        Irr_Data%spectrum(1:nwavel,i) = &
            Irr_Data%spectrum(1:Irr_Data%nwaves(i),i) / &
            (SUM(Irr_Data%spectrum(1:nwavel,i)*weightsum(1:nwavel)) / SUM(weightsum(1:nwavel)))
      end do
    endif

    deallocate (spectra, wavelengths, qflags, weightsum, stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_I0_irradiance: deallocation error", errstat)
    endif

  end subroutine read_I0_irradiance

  !> Subroutine to read precalculated irradiance-like I0 spectra
  !--------------------------------------------------------------------------
  !
  ! @param[in]  filename      Name of I0 spectra netCDF4 file
  ! @param[out]  nxtrack       size of xtrack dimension
  ! @param[out]  nwavel        size of I0 wavelength dimension
  ! @param      errstat       error tracking integer (non-zero = problem)
  !
  ! @ author  E. O'Sullivan  November 2020
  !--------------------------------------------------------------------------
  subroutine read_I0_dims (filename, nxtrack, nwavel, errstat)
    use tio_module
    implicit none

    !input variables
    character (len=*), intent(in) :: filename
    !output variables
    integer (kind=4), intent(out) :: nxtrack
    integer (kind=4), intent(out) :: nwavel
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: tio_obj

    if (errstat /= 0) return

    call open_nc (trim(adjustl(filename)), tio_obj, errstat)
    call tiof_inq_dimlen (tio_obj, "xtrack", nxtrack, errstat)
    call tiof_inq_dimlen (tio_obj, "spectral_channel", nwavel, errstat)
    call close_nc (tio_obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_IO_dims: failed", errstat)
      return
    endif

  end subroutine read_I0_dims

  !> Subroutine to read precalculated irradiance-like I0 spectra
  !--------------------------------------------------------------------------
  !
  ! @param[in]  filename      Name of I0 spectra netCDF4 file
  ! @param[in]  nxtrack       size of xtrack dimension
  ! @param[in]  nwavel        size of I0 wavelength dimension
  ! @param[out] spectra       2D array of spectra
  ! @param[out] wavelength    2D array of wavelengths
  ! @param[out] quality_flag  2D array of quality flags (non-zero = bad)
  ! @param      errstat       error tracking integer (non-zero = problem)
  !
  ! @ author  E. O'Sullivan  November 2020
  !--------------------------------------------------------------------------
  subroutine read_I0_spectra (filename, nxtrack, nwavel, spectra, wavelength, &
    quality_flag, errstat)
    use tio_module
    implicit none

    !input variables
    character (len=*), intent(in) :: filename
    integer (kind=4), intent(in) :: nxtrack
    integer (kind=4), intent(in) :: nwavel
    !output variables
    real (kind=8), dimension(:,:), intent(out) :: spectra, &
          wavelength
    integer (kind=2), dimension(:,:), intent(out) :: quality_flag
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: tio_obj

    if (errstat /= 0) return

    call open_nc (trim(adjustl(filename)), tio_obj, errstat)
    call tiof_get2d_r8 (tio_obj, "spectrum", [0,0], [nxtrack, nwavel], &
          spectra, errstat)
    call tiof_get2d_r8 (tio_obj, "wavelength", [0,0], [nxtrack, nwavel], &
          wavelength, errstat)
    call tiof_get2d_i2 (tio_obj, "quality_flag", [0,0], [nxtrack, nwavel], &
          quality_flag, errstat)
    call close_nc (tio_obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_IO_spectra: failed", errstat)
      return
    endif

  end subroutine read_I0_spectra

  !---------------------------------------------------------------------
  !
  !> @param l1file filename for L1 netCDF file
  !> @param tio_l1obj file object
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2018
  !---------------------------------------------------------------------
  subroutine open_nc (l1file, tio_l1obj, errstat)
    use tio_module
    use netcdf, only: nf90_nowrite
    implicit none

    !input variables
    character (len=*), intent (in) :: l1file

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !local variables
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    call tiof_open (l1file, tio_l1obj, nf90_nowrite, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "open_nc: failed to open L1 file", &
           errstat)
      return
    endif

  end subroutine open_nc
 
  !> Subroutine to close a netCDF file
  !---------------------------------------------------------------------
  !
  !> @param[in] tio_l1obj file object
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan October 2018
  !---------------------------------------------------------------------
  subroutine close_nc (tio_l1obj, errstat)
    use tio_module
    implicit none

    !input variables

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !local variables
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "close_nc: failed to close L1 file", &
           errstat)
      return
    endif

  end subroutine close_nc
 
end module
