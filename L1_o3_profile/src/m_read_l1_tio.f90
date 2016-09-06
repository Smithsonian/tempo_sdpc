!> Read radiance and geolocation data from L1B netCDF file
module m_read_l1_tio
  use o3p_names_module
  use tio_module
  use tell_module
  use netcdf, only: nf90_nowrite
  use m_convert_coadd, only: convert_2bytes_to_16bits
  use OMSAO_omidata_module, only: ncoadd, nxtrack_max, ntimes_max 
! nxtrack, ntimes, nfxtrack, 
  use OMSAO_variables_module, only: scnwrt


  public read_L1_dims_tio, read_L1_rad_line_tio, &
       open_L1_tio, close_L1_tio
  private

contains

  !> Subroutine to read in dimensions of an L1 netCDF file
  !---------------------------------------------------------------------
  !
  !> @param[in] l1file filename for L1 netCDF file
  !> @param[in] swathname swath name in L1 netCDF file
  !> @param[out] nstep size of along-track dimension
  !> @param[out] nxtrack size of xtrack dimension
  !> @param[out] OPTIONAL size of xtrack dimension * ncoadd
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan June 2016
  !---------------------------------------------------------------------
  subroutine read_L1_dims_tio (l1file, swathname, nstep, nxtrack, &
       nxtrack_coadd, errstat)

    implicit none
    
    !input variables
    character (len=*), intent (in) :: l1file, swathname

    !output variables
    integer (kind=4), intent (out) :: nstep, nxtrack
    integer (kind=4), optional, intent (out) :: nxtrack_coadd
    integer (kind=4), intent (inout) :: errstat

    !local variables
    type (tiof_file_type) :: tio_l1obj


    if (errstat /= 0) return

    call tiof_open (l1file, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (tio_l1obj, o3p_dim_step, nstep, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, o3p_dim_xtrack, nxtrack, errstat)
    call tiof_close (tio_l1obj, errstat)
    
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_L1_dims_tio: failed to open L1 file", &
           errstat)
      return
    endif

    ! Binning factor
    if (present(nxtrack_coadd)) then
      nxtrack_coadd = nxtrack * ncoadd
    endif

    ! check dimension sizes are within bounds
    if (nstep > ntimes_max) then
      errstat = -1
      call tell_error (tell_invalid_parm, &
           "read_L1_dims_tio: ntimes greater than ntimes_max", errstat)
      return
    else if (nxtrack > nxtrack_max) then
      errstat = -1
      call tell_error (tell_invalid_parm, &
           "read_L1_dims_tio: nxtrack greater than nxtrack_max", errstat)
      return
    endif
    if (present(nxtrack_coadd)) then
      if (nxtrack > nxtrack_max) then
        errstat = -1
        call tell_error (tell_invalid_parm, &
             "read_L1_dims_tio: nxtrack_coadd greater than nxtrack_max", &
             errstat)
        return
      endif
    endif

  end subroutine read_L1_dims_tio 



  !> Subroutine to read one line of (ir)radiances from an L1 netCDF file
  !---------------------------------------------------------------------
  !
  !> @param[in] tio_l1obj (ir)radiance file object
  !> @param[in] swathname swath to read data from
  !> @param[in] step along-track step to read data from
  !> @param[out] radiance array of radiance (or irradiance) values
  !> @param[out] rad_precision array of radiance precision values
  !> @param[out] pixel_quality_flag array of pixel quality flag values
  !> @param[out] wavelengths array of wavelength values
  !> @param[out] meas_qual_flag measurement quality flag for the step
  !> @param[out] num_wavelengths size of wavelength axis
  !> @param[in] read_irrad set to true if reading from irradiance file
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan July 2016
  !---------------------------------------------------------------------
  subroutine read_L1_rad_line_tio (tio_l1obj, swathname, step, &
       radiance, rad_precision, pixel_quality_flag, wavelengths, &
       meas_qual_flag, num_wavelengths, read_irrad, errstat)

    implicit none
    
    !input variables
    character (len=*), intent (in) :: swathname
    integer, intent (in) :: step
    logical, intent (in) :: read_irrad
    
    !output variables
    real (kind=4), dimension(:,:), intent(out) :: &
         radiance, rad_precision, wavelengths
    integer (kind=2), dimension (:,:), intent(out) :: &
         pixel_quality_flag
    integer (kind=4), intent(out) :: num_wavelengths
    integer (kind=2), intent(out) :: meas_qual_flag
    integer (kind=4), intent(inout) :: errstat

    !local_variables
    real (kind=4), dimension(:,:,:), allocatable :: tio_rad, &
         tio_prec, tio_wvl
    integer (kind=2), dimension(:,:,:), allocatable :: tio_pqf
    integer (kind=2), dimension(1) :: tio_mqf
    character (len=25) :: rad_param_name, prec_param_name
    integer :: nx,nw

    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    !Setup to read either radiance or irradiance
    if (read_irrad) then
      rad_param_name = o3p_var_irradiance
      !prec_param_name = o3p_var_irradiance_prec
      !prec_param_name = o3p_var_irradiance_relerr
      prec_param_name = o3p_var_irradiance_error
    else
      rad_param_name = o3p_var_radiance
      !prec_param_name = o3p_var_radiance_prec
      !prec_param_name = o3p_var_radiance_relerr
      prec_param_name = o3p_var_radiance_error
    endif

    call tiof_push_group (tio_l1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, o3p_dim_channel, nw, errstat)
    call tiof_inq_dimlen (tio_l1obj, o3p_dim_xtrack, nx, errstat)

    allocate(tio_rad(nw, nx, 1), tio_prec(nw, nx, 1), tio_wvl(nw, nx, 1), &
         tio_pqf(nw, nx, 1), stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_L1_rad_line_tio: failed to allocate memory", &
           errstat)
      return
    endif
 
!    call tiof_get3d_r4 (tio_l1obj, o3p_var_radiance, [step,0,0], &
!         [1,nx,nw], tio_rad, errstat)
!    call tiof_get3d_r4 (tio_l1obj, o3p_var_radiance_prec, [step,0,0], &
!         [1,nx,nw], tio_prec, errstat)
    call tiof_get3d_r4 (tio_l1obj, rad_param_name, [step,0,0], &
         [1,nx,nw], tio_rad, errstat)
    call tiof_get3d_r4 (tio_l1obj, prec_param_name, [step,0,0], &
         [1,nx,nw], tio_prec, errstat)
    call tiof_get3d_r4 (tio_l1obj, o3p_var_wavelength, [step,0,0], &
         [1,nx,nw], tio_wvl, errstat)
    call tiof_get3d_i2 (tio_l1obj, o3p_var_pqf, [step,0,0], &
         [1,nx,nw], tio_pqf, errstat)
    call tiof_get1d_i2 (tio_l1obj, o3p_var_mqf, [step], &
         [1], tio_mqf, errstat)
    call tiof_pop_group (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_L1_rad_line_tio: failed to read data from L1 file", &
           errstat)
      return
    endif

    ! if using relerr convert relative error to actual error
    !tio_prec(1:nw,1:nx,1) = tio_rad(1:nw,1:nx,1)*tio_prec(1:nw,1:nx,1)

    ! shift data into output arrays
    radiance(1:nw,1:nx) = tio_rad(1:nw,1:nx,1)
    rad_precision(1:nw,1:nx) = tio_prec(1:nw,1:nx,1)
    wavelengths(1:nw,1:nx) = tio_wvl(1:nw,1:nx,1)
    pixel_quality_flag(1:nw,1:nx) = tio_pqf(1:nw,1:nx,1)
    meas_qual_flag = tio_mqf(1)
    num_wavelengths = nw

    if(allocated(tio_rad)) deallocate(tio_rad, stat=errstat)
    if(allocated(tio_prec)) deallocate(tio_prec, stat=errstat)
    if(allocated(tio_wvl)) deallocate(tio_wvl, stat=errstat)
    if(allocated(tio_pqf)) deallocate(tio_pqf, stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_L1_rad_line_tio: failed to deallocate memory", &
           errstat)
      return
    endif


  end subroutine read_L1_rad_line_tio






  !> Subroutine to open an L1 netCDF file
  !---------------------------------------------------------------------
  !
  !> @param l1file filename for L1 netCDF file
  !> @param tio_l1obj file object
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan July 2016
  !---------------------------------------------------------------------
  subroutine open_L1_tio (l1file, tio_l1obj, errstat)

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
           "open_L1_tio: failed to open L1 file", &
           errstat)
      return
    endif

  end subroutine open_L1_tio


  !> Subroutine to close an L1 netCDF file
  !---------------------------------------------------------------------
  !
  !> @param[in] tio_l1obj file object
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan July 2016
  !---------------------------------------------------------------------
  subroutine close_L1_tio (tio_l1obj, errstat)

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
           "close_L1_tio: failed to close L1 file", &
           errstat)
      return
    endif

  end subroutine close_L1_tio




end module m_read_l1_tio
