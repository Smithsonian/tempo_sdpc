! read and convolve with TEMPO preflight or irradiance file slit function
MODULE slitfunction_tempo
  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: nxtrack_max
  IMPLICIT NONE

  public tempo_slitfunc_read, tempo_slitfunc_convolve

  PRIVATE

  ! Mean supergaussian slit function parameters
  real (kind=8) :: mean_asym, mean_hw1e, mean_shape, mean_wl

CONTAINS

  !------------------------------------------------------------------------
  !
  ! @param  errstat  error-tracking variable, non-zero = failure
  !
  ! @author  E.O'Sullivan July 2020
  !------------------------------------------------------------------------
  SUBROUTINE tempo_slitfunc_read ( errstat )

    USE OMSAO_variables_module,  ONLY: omi_slitfunc_fname, l1b_channel, &
         ctrl_fit_winwav_lim
    use omsao_parameters_module, only: N_FIT_WINWAV
    use l1bread_utils, only: lookup_swathname
    use tio_module
    use tell_module
    use netcdf, only: nf90_nowrite, nf90_noerr, nf90_global, nf90_enotatt, &
         nf90_get_att, nf90_inquire_attribute

    IMPLICIT NONE

    !input variables
    integer (kind=4), intent(inout) :: errstat

    ! --------------
    ! Local variable
    ! --------------
    character (len=15) :: swathname
    character (len=3) :: prod_str
    integer (kind=4) :: ncerr, n
    integer (kind=4), parameter :: nxtrack = 2048, nwl = 1024
    logical :: preflight
    real (kind=4), dimension(nwl, nxtrack) :: sf_asym, sf_hw1e, sf_shape, &
         sf_wavelength
    real (kind=4), dimension(nwl, nxtrack,1) :: tmp_asym, tmp_hw1e, tmp_shape
    real (kind=4), dimension(nwl) :: tmp_wl
    real (kind=4) :: minwl, maxwl
    logical, dimension(nwl, nxtrack) :: mask
    
    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    ! determine which band we're using
    call lookup_swathname (l1b_channel, swathname, errstat)

    ! determine whether we're reading from irradiance or pre-flight slit func.
    call tiof_open (omi_slitfunc_fname, tio_l1obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "tempo_slitfunc_read: error opening slit function file", errstat)
      return
    endif
    prod_str=""
    ncerr = nf90_inquire_attribute (tio_l1obj%fileid, nf90_global, &
         "product_type")
    if (ncerr == nf90_enotatt) then
      preflight = .true.
    else if (ncerr == nf90_noerr) then
      ncerr = nf90_get_att (tio_l1obj%fileid, nf90_global, "product_type", &
           prod_str)
      if (ncerr .ne. nf90_noerr .or. prod_str .ne. "IRR") then
        call tell_error (tell_io_open_error, &
             "tempo_slitfunc_read: undetermined slit function file type", &
             errstat)
        return
      else
        preflight = .false.
      endif
    endif

    if (preflight) then ! read from pre flight slit function
      call tiof_push_group (tio_l1obj, swathname, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_asym", [0,0], [nxtrack,nwl], &
           sf_asym, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_hw1e", [0,0], [nxtrack,nwl], &
           sf_hw1e, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_shape", [0,0], [nxtrack,nwl], &
           sf_shape, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_wavelength", [0,0], [nxtrack,nwl], &
           sf_wavelength, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "tempo_slitfunc_read: failed to read pre-flight file", &
             errstat)
        return
      endif
      call tell_log (0, "read pre-flight slit function")
    else ! read from irradiance file
      call tiof_push_group (tio_l1obj, swathname, errstat)
      call tiof_get3d_r4 (tio_l1obj, "sf_asym", [0,0,0], [1,nxtrack,nwl], &
           tmp_asym, errstat)
      call tiof_get3d_r4 (tio_l1obj, "sf_hw1e", [0,0,0], [1,nxtrack,nwl], &
           tmp_hw1e, errstat)
      call tiof_get3d_r4 (tio_l1obj, "sf_shape", [0,0,0], [1,nxtrack,nwl], &
           tmp_shape, errstat)
      call tiof_get1d_r4 (tio_l1obj, "nominal_wavelength", [0], [nwl], &
           tmp_wl, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "tempo_slitfunc_read: failed to read pre-flight file", &
             errstat)
        return
      endif
      call tell_log (0, "read slit function from irradiance file")
      sf_asym = tmp_asym(:,:,1)
      sf_hw1e = tmp_hw1e(:,:,1)
      sf_shape = tmp_shape(:,:,1)
      do n=1,nxtrack
        sf_wavelength(:,n) = tmp_wl(:)
      enddo
    endif

    ! determine the mean variable values in the wavelength window
    mask = .true.
    minwl = real(ctrl_fit_winwav_lim(2), kind=4)
    maxwl = real(ctrl_fit_winwav_lim(N_FIT_WINWAV-1), kind=4)
    where (sf_wavelength .lt. minwl .or. sf_wavelength .gt. maxwl)
      mask = .false.
    end where
    !also mask out any bad pixels in the slit function files
    where (sf_asym .gt. 1000.0 .or. sf_asym .lt. 0.0)
      mask = .false.
    end where
    mean_asym = real(sum(sf_asym,mask=mask)/count(mask), kind=8)
    mean_hw1e = real(sum(sf_hw1e,mask=mask)/count(mask), kind=8)
    mean_shape = real(sum(sf_shape,mask=mask)/count(mask), kind=8)
    mean_wl = real(sum(sf_wavelength,mask=mask)/count(mask), kind=8)
print *, mean_asym, mean_hw1e, mean_shape, mean_wl

  END SUBROUTINE tempo_slitfunc_read


  !--------------------------------------------------------------------------
  !
  ! @param[in]  xtrack_pix  current cross-track pixel number
  ! @param[in]  nwl         number of wavelengths in current spectrum
  ! @param[in]  wvl         wavelength array of current spectrum
  ! @param[in]  spec        current spectrum
  ! @param[out] spec_conv   convolved spectrum
  ! @param      errstat     error handling integer, non-zero = failure
  !
  ! @author  E. O'Sullivan, July 2020
  !--------------------------------------------------------------------------
  SUBROUTINE tempo_slitfunc_convolve ( xtrack_pix, nwvl, wvl, spec, &
       spec_conv, errstat )

    use slitfunction_super_gaussian, only: super_gaussian_sf
    use tell_module

    implicit none

    !input variables
    integer (kind=4), intent(in) :: xtrack_pix, nwvl
    real (kind=8), dimension(:), intent(in) :: wvl, spec

    !output variables
    real (kind=8), dimension(:), intent(out) :: spec_conv
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    call super_gaussian_sf (nwvl, mean_hw1e, mean_asym, mean_shape, 0.0_r8, &
         wvl, spec, spec_conv)

  END SUBROUTINE tempo_slitfunc_convolve

END MODULE slitfunction_tempo
