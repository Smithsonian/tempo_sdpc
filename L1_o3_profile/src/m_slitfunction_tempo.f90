!> Read TEMPO slit function from preflight slit file or irradiance file
module m_slitfunction_tempo

  implicit none

  public tempo_slitfunc_read
  private

contains

  !-----------------------------------------------------------------------
  ! Read parameters of a supergaussian from either the irradiance file
  ! or a pre-flight slit function file, and detectermine the mean value
  ! of those parameters in each fitting window, as defined by wavelength
  ! limits. Note that no cross-track limits are imposed.
  !
  ! @param[out] mean_hw1e  1D array of mean half-widths in each window
  ! @param[out] mean_asym  1D array of mean asymmetry parameters
  ! @param[out] mean_shape 1D array of mean shape parameters
  ! @param      errstat    error tracking integer, non-zero mean failure
  !
  ! @author  E.O'Sullivan  April 2021
  !-----------------------------------------------------------------------
  subroutine tempo_slitfunc_read (mean_hw1e, mean_asym, mean_shape, &
       errstat)

    use OMSAO_variables_module, only: winlim, numwin
    use OMSAO_tmpodata_module, only: irrad_swathname
    use OMSAO_parameters_module, only: maxchlen
    use OMI_LUN_set, only: slitfunc_lun
    use PGS_PC_class
    use tio_module
    use tell_module
    use netcdf, only: nf90_nowrite, nf90_noerr, nf90_global, nf90_enotatt, &
         nf90_get_att, nf90_inquire_attribute

    implicit none

    !output_variables
    real (kind=4), dimension(numwin), intent(out) :: mean_hw1e, mean_asym, &
         mean_shape
    integer (kind=4), intent(inout) :: errstat
    integer :: version = 1

    !local variables
    character (len=maxchlen) :: filename
    character (len=3) :: prod_str
    character (len=128) :: msg_str
    integer (kind=4) :: ncerr, nw
    integer (kind=4), parameter :: nxtrack = 2048, nwl = 1024
    logical :: preflight
    real (kind=8), dimension(nwl, nxtrack) :: sf_asym, sf_hw1e, sf_shape, &
         sf_wavelength
    real (kind=4), dimension(nwl, nxtrack,1) :: tmp_asym, tmp_hw1e, &
         tmp_shape, tmp_wl
    real (kind=4), dimension(nwl,nxtrack) :: r4_asym, r4_hw1e, r4_shape, r4_wavelength
    real (kind=4), dimension(numwin) :: mean_wl
    real (kind=4) :: minwl, maxwl
    logical, dimension(nwl, nxtrack) :: mask

    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    ! Get filename from PCF file
    errstat = PGS_PC_getreference( slitfunc_lun, version, filename )
    if( errstat /= 0 ) then
      call tell_error (tell_io_read_error, &
           "m_slitfunction_tempo: slitfunction file not found in PCF", errstat)
      return
    end if

    ! determine whether we're reading from irradiance or pre-flight slit func.
    call tiof_open (trim(adjustl(filename)), tio_l1obj, nf90_nowrite, errstat)
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

    !Each window may reference any swath, so for each in turn, read the
    !slit parameter arrays and determine the mean values
    do nw=1,numwin
      sf_asym=0.0d0
      sf_hw1e=0.0d0
      sf_shape=0.0d0
      sf_wavelength=0.0d0
      if (preflight) then ! read from pre flight slit function
        call tiof_push_group (tio_l1obj, trim(adjustl(irrad_swathname(nw))), &
             errstat)
        call tiof_get2d_r4 (tio_l1obj, "sf_asym", [0,0], [nxtrack,nwl], &
             r4_asym, errstat)
        call tiof_get2d_r4 (tio_l1obj, "sf_hw1e", [0,0], [nxtrack,nwl], &
             r4_hw1e, errstat)
        call tiof_get2d_r4 (tio_l1obj, "sf_shape", [0,0], [nxtrack,nwl], &
             r4_shape, errstat)
        call tiof_get2d_r4 (tio_l1obj, "sf_wavelength", [0,0], [nxtrack,nwl], &
             r4_wavelength, errstat)
        if (errstat /= 0) then
          call tell_error (tell_io_read_error, &
               "tempo_slitfunc_read: failed to read pre-flight file", &
               errstat)
          return
        endif
        ! note we do maths in double precision to avoid rounding errors
        sf_asym=real(r4_asym, kind=8)
        sf_hw1e=real(r4_hw1e, kind=8)
        sf_shape=real(r4_shape, kind=8)
        sf_wavelength=real(r4_wavelength, kind=8)
        write(msg_str,'(A,I2)') "read pre-flight slit function for window ",nw
        call tell_log (0, trim(adjustl(msg_str)))
        call tiof_pop_group (tio_l1obj, errstat)
      else ! read from irradiance file
        call tiof_push_group (tio_l1obj, trim(adjustl(irrad_swathname(nw))), &
             errstat)
        call tiof_get3d_r4 (tio_l1obj, "sf_asym", [0,0,0], [1,nxtrack,nwl], &
             tmp_asym, errstat)
        call tiof_get3d_r4 (tio_l1obj, "sf_hw1e", [0,0,0], [1,nxtrack,nwl], &
             tmp_hw1e, errstat)
        call tiof_get3d_r4 (tio_l1obj, "sf_shape", [0,0,0], [1,nxtrack,nwl], &
             tmp_shape, errstat)
        call tiof_get3d_r4 (tio_l1obj, "wavelength", [0,0,0], [1,nxtrack,nwl],&
             tmp_wl, errstat)
        if (errstat /= 0) then
          call tell_error (tell_io_read_error, &
               "tempo_slitfunc_read: failed to read pre-flight file", &
               errstat)
          return
        endif
        write(msg_str,'(A,I2)') "read irradiance slit function for window ",nw
        call tell_log (0, trim(adjustl(msg_str)))
        call tiof_pop_group (tio_l1obj, errstat)
        sf_asym = real(tmp_asym(:,:,1), kind=8)
        sf_hw1e = real(tmp_hw1e(:,:,1), kind=8)
        sf_shape = real(tmp_shape(:,:,1), kind=8)
        sf_wavelength = real(tmp_wl(:,:,1), kind=8)
      endif
      ! get mean values
      mask = .true.
      minwl = real(winlim(nw,1), kind=4)
      maxwl = real(winlim(nw,2), kind=4)
      where (sf_wavelength .lt. minwl .or. sf_wavelength .gt. maxwl)
        mask = .false.
      end where
      !also mask out any bad pixels in the slit function files
      where (sf_asym .gt. 1000.0 .or. sf_asym .lt. 0.0)
        mask = .false.
      end where
      mean_asym(nw) = real(sum(sf_asym,mask=mask)/count(mask), kind=4)
      mean_hw1e(nw) = real(sum(sf_hw1e,mask=mask)/count(mask), kind=4)
      mean_shape(nw) = real(sum(sf_shape,mask=mask)/count(mask), kind=4)
      mean_wl(nw) = real(sum(sf_wavelength,mask=mask)/count(mask), kind=4)
      write(msg_str,'(I2,4(1x,f9.5))') nw, mean_asym(nw), mean_hw1e(nw), &
           mean_shape(nw), mean_wl(nw)
      call tell_log(0, msg_str)
    enddo

    call tiof_close(tio_l1obj,errstat)

  end subroutine tempo_slitfunc_read

end module m_slitfunction_tempo
