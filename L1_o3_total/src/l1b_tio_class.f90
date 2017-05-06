!> Level 1 data input functions
!! @file
!! @note Ideally, \a use statements in this module will provide
!!       access only to subroutines, type definitions, and compile-time
!!       constants. All variables should be passed via the subroutine
!!       parameter argument lists.
module l1b_tio_class
  use netcdf
  use tio_module
  use tell_module
  use o3t_names_module
  implicit none
  private

  integer, parameter :: nblock = 10 ! default number of lines in a block

  !> Level 1 file object
  !! @note This is an opaque data type.
  type, public :: l1b_tio_type
    private
    type (tiof_file_type) :: ft
    character (len=1024) :: filename
    integer :: num_steps
    integer :: num_xtrack
    integer :: num_wavelengths
    real (kind=4) :: distance   ! earth-sun distance
  end type

  !> Level 1 metadata object
  type, public :: l1b_metadata_type
    integer :: year, month, day, jday
    integer :: orbit_number
    character (len=32) :: shortname
  end type

  !> geolocated radiance spectra
  type, public :: l1b_radgeo_type
    real (kind=4), dimension(:,:), allocatable :: lon, lat, sza, saz, vza, vaz
    real (kind=4), dimension(:,:,:), allocatable :: lon_bounds, lat_bounds
    real (kind=4), dimension(:,:,:), allocatable :: radiance, wavelength
    integer (kind=4), dimension(:), allocatable :: step_indices
    integer (kind=2), dimension(:,:,:), allocatable :: qa_flags
    integer (kind=2), dimension(:,:), allocatable :: geoflg, hgt
    integer (kind=1), dimension(:,:), allocatable :: anomflg
    integer (kind=2), dimension(:), allocatable :: measurement_quality_flags
    integer (kind=1), dimension(:), allocatable :: instid
    integer :: beg_line, end_line, num_lines
  end type

  public l1b_tio_open, l1b_tio_close, l1b_tio_getdims, l1b_tio_get_irr, &
    l1b_tio_earthsun_distance, l1b_tio_init_rad, l1b_tio_getrad, &
    l1b_tio_get_etc, l1b_tio_getgeo, l1b_tio_get_metadata, l1b_file_object

contains

  function l1b_file_object (this) result(obj)
    implicit none
    type (l1b_tio_type), intent(in) :: this
    type (tiof_file_type) :: obj

    obj = this % ft

  end function l1b_file_object

  !> Open Level 1 data product file
  !! @param[inout]  this  Level 1 file object
  !! @param[in] filename  Level 1 data product file name
  !! @param[in] groupname  Name of group to open
  !! @param[inout] errstat  Error status variable
  subroutine l1b_tio_open (this, filename, groupname, errstat)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    character (len=*), intent(in) :: filename, groupname
    integer, intent(inout) :: errstat

    if (errstat < 0) return

    call tiof_open (filename, this % ft, nf90_nowrite, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_open_error, &
                       "l1b_tio_open:  opening "//trim(filename), errstat)
      return
    endif

    this % filename = filename

    call tiof_get_r4 (this % ft, o3t_var_earth_sun_distance, this % distance, errstat)
    call tiof_inq_dimlen (this % ft, o3t_dim_step, this % num_steps, errstat)
    call tiof_inq_group (this % ft, groupname, errstat)
    call tiof_inq_dimlen (this % ft, o3t_dim_xtrack, this % num_xtrack, errstat)
    call tiof_inq_dimlen (this % ft, o3t_dim_channel, this % num_wavelengths, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
                       "l1b_tio_open:  reading "//trim(filename), errstat)
      return
    endif
  end subroutine

  !> Close Level 1 data product file
  !! @param[inout]  this  Level 1 file object
  !! @param[inout] errstat  Error status variable
  subroutine l1b_tio_close (this, errstat)
    implicit none
    type(l1b_tio_type), intent(inout) :: this
    integer, intent(inout) :: errstat
    if (errstat < 0) continue  ! no-op
    this % num_steps = -1
    this % num_xtrack = -1
    this % num_wavelengths = -1
    this % filename = ''
  end subroutine

  !> Get dimension sizes from an open \a l1b_tio_type object
  !! @param[inout]  this  Level 1 file object
  !! @param[out]  num_steps  Number of scan steps
  !! @param[out]  num_xtrack  Number of cross-track pixels
  !! @param[out] num_wavelengths Number of wavelengths
  !! @param[inout] errstat  Error status variable
  subroutine l1b_tio_getdims (this, num_steps, num_xtrack, num_wavelengths, errstat)
    implicit none
    type(l1b_tio_type), intent(inout) :: this
    integer, intent(out) :: num_steps, num_xtrack, num_wavelengths
    integer, intent(inout) :: errstat

    if (errstat < 0) return

    num_steps = this % num_steps
    num_xtrack = this % num_xtrack
    num_wavelengths = this % num_wavelengths
  end subroutine

  !> Get Earth-Sun distance from an open \a l1b_tio_type object
  !! @param[inout]  this  Level 1 file object
  !! @param[out] distance  Earth-Sun distance
  !! @param[inout] errstat  Error status variable
  subroutine l1b_tio_earthsun_distance (this, distance, errstat)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    real (kind=4), intent(out) :: distance
    integer, intent(inout) :: errstat
    if (errstat < 0) return
    distance = this % distance
  end subroutine

  !> Read irradiance spectra from an open \a l1b_tio_type object
  !! @param[inout]  this  Level 1 file object
  !! @param[out] irr  Array of irradiance spectra, one for each cross-track pixel.
  !! @param[out] wavelength  Array of irradiance wavelength grids, one for each cross-track pixel.
  !! @param[out] qaflags  Array of quality assurance flags, one for each spectral bin.
  !! @param[out] measurement_quality_flag  Measurement quality flag
  !! @param[out] instid  Instrument id
  !! @param[inout] errstat  Error status variable
  subroutine l1b_tio_get_irr (this, irr, wavelength, qaflags, &
                              measurement_quality_flag, instid, errstat)
    use netcdf, only : nf90_noerr, nf90_inq_varid ! FIXME <---
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    real (kind=4), dimension(:,:), intent(out) :: irr, wavelength
    integer (kind=2), dimension(:,:), intent(out) :: qaflags
    integer (kind=2), intent(out) :: measurement_quality_flag
    integer (kind=1), intent(out) :: instid
    integer, intent(inout) :: errstat

    integer :: nx, nw, status, instid_varid
    real (kind=4), dimension(:,:,:), allocatable :: &
      tmp_wavelengths, tmp_spectrum
    integer (kind=2), dimension(:,:,:), allocatable :: tmp_qflags
    integer (kind=2), dimension(1) :: tmp_mqf
    integer (kind=1), dimension(1) :: tmp_instid

    if (errstat < 0) return

    nx = this % num_xtrack
    nw = this % num_wavelengths

    if ((size(irr,1) < nw .or. size(irr,2) < nx) &
        .or. (size(wavelength,1) < nw .or. size(wavelength,2) < nx) &
        .or. (size(qaflags,1) < nw &
              .or. size(qaflags,2) < nx)) then
      call tell_error (tell_runtime_error, &
                       "l1b_tio_get_irr: array too small", errstat)
      return
    endif

    allocate (tmp_spectrum(nw, nx,1), &
              tmp_wavelengths(nw, nx,1), &
              tmp_qflags (nw, nx,1), &
              stat=status)
    if (status /= 0) then
      call tell_error (tell_malloc_error, &
                       'l1b_tio_get_irr:  allocate failed', errstat)
      return
    endif

    ! FIXME:  This instid stuff is here only to support tests with OMI data.
    !         Maybe it will eventually go away?
    status = nf90_inq_varid (this%ft%groupid, o3t_var_instid, instid_varid)
    if (status == nf90_noerr) then
      call tiof_get1d_ui1 (this % ft, o3t_var_instid, [0], [1], &
                           tmp_instid(1:1), errstat)
      instid = tmp_instid(1)
    else
      instid = 0_1
    endif

    call tiof_get3d_r4 (this % ft, o3t_var_wavelength, [0,0,0], [1,nx,nw], &
                        tmp_wavelengths(:,1:nx,1:1), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_irradiance, [0,0,0], [1,nx,nw], &
                        tmp_spectrum(:,1:nx,1:1), errstat)
    call tiof_get3d_ui2 (this % ft, o3t_var_pqf, [0,0,0], [1,nx,nw], &
                        tmp_qflags(:,1:nx,1:1), errstat)
    call tiof_get1d_ui2 (this % ft, o3t_var_mqf, [0], [1], tmp_mqf(1:1), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
                       'l1b_tio_get_irr: failed reading irradiance: '//trim(this % filename), &
                       errstat)
      return
    endif

    irr(1:nw, 1:nx) = tmp_spectrum(1:nw, 1:nx, 1)
    wavelength(1:nw, 1:nx) = tmp_wavelengths(1:nw, 1:nx, 1)
    qaflags(1:nw, 1:nx) = tmp_qflags (1:nw, 1:nx, 1)
    measurement_quality_flag = tmp_mqf(1)

  end subroutine

  !> Initialize an \a l1b_radgeo_type object using an open \a l1b_tio_type
  !! @param[inout] this  Level 1 file object
  !! @param[inout] rg  The \a l1b_radgeo_type object to initialize
  !! @param[inout] errstat  Error status variable
  !! @details
  !! For efficiency, radiance spectra will be loaded in blocks.
  !! The \a l1b_radgeo_type object is initialized by allocating space
  !! for all variables, defining the number of radiance spectra per cached
  !! block (for processing efficiency), and loading the geolocation variables.
  subroutine l1b_tio_init_rad (this, rg, errstat)
    implicit none
    type(l1b_tio_type), intent(inout), target :: this
    type(l1b_radgeo_type), intent(inout) :: rg
    integer, intent(inout) :: errstat

    integer :: i, nx, ns, nw, nl, ierr
    integer (kind=4), allocatable :: tmp_geoflg(:,:)

    if (errstat < 0) return

    rg % beg_line = -1
    rg % end_line = -1
    rg % num_lines = min(nblock, this % num_steps)

    nl = rg % num_lines
    nx = this % num_xtrack
    nw = this % num_wavelengths
    ns = this % num_steps

    allocate (rg % lon (nx, ns), &
              rg % lat (nx, ns), &
              rg % lon_bounds (4, nx, ns), &
              rg % lat_bounds (4, nx, ns), &
              rg % sza (nx, ns), &
              rg % saz (nx, ns), &
              rg % vza (nx, ns), &
              rg % vaz (nx, ns), &
              rg % hgt (nx, ns), &
              rg % geoflg (nx, ns), &
              rg % anomflg (nx, ns), &
              rg % radiance (nw, nx, nl), &
              rg % wavelength (nw, nx, nl), &
              rg % step_indices (ns), &
              rg % qa_flags (nw, nx, nl), &
              rg % measurement_quality_flags (nl), &
              rg % instid (nl), &
              tmp_geoflg (nx, ns), &
              stat = ierr)
    if (ierr /= 0) then
      call tell_error (tell_malloc_error, "l1b_tio_init_rad: allocate failed", errstat)
      return
    endif

    ! If present, read mirror step indices from input radiance file
    call tell_push_queue
    call tiof_push_group (this % ft, "/", errstat)
    call tiof_get1d_i4 (this % ft, o3t_dim_step, [0], [ns], &
                        rg % step_indices(1:ns), errstat)
    if (errstat == 0) then
      call tell_pop_queue(0)
    else
      rg % step_indices(1:ns) = [(i, i=0,ns-1)]
      call tell_pop_queue (1)
    endif
    call tiof_pop_group (this % ft, errstat)

    call tiof_get2d_r4 (this % ft, o3t_var_longitude, [0,0], [ns,nx], &
                        rg % lon(1:nx,1:ns), errstat)
    call tiof_get2d_r4 (this % ft, o3t_var_latitude, [0,0], [ns,nx], &
                        rg % lat(1:nx,1:ns), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_longitude_bounds, [0,0,0], [ns,nx,4], &
                        rg % lon_bounds(1:4,1:nx,1:ns), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_latitude_bounds, [0,0,0], [ns,nx,4], &
                        rg % lat_bounds(1:4,1:nx,1:ns), errstat)
    call tiof_get2d_r4 (this % ft, o3t_var_sz_angle, [0,0], [ns,nx], &
                        rg % sza(1:nx,1:ns), errstat)
    call tiof_get2d_r4 (this % ft, o3t_var_sa_angle, [0,0], [ns,nx], &
                        rg % saz(1:nx,1:ns), errstat)
    call tiof_get2d_r4 (this % ft, o3t_var_vz_angle, [0,0], [ns,nx], &
                        rg % vza(1:nx,1:ns), errstat)
    call tiof_get2d_r4 (this % ft, o3t_var_va_angle, [0,0], [ns,nx], &
                        rg % vaz(1:nx,1:ns), errstat)
    call tiof_get2d_i2 (this % ft, o3t_var_terrain_height, [0,0], [ns,nx], &
                        rg % hgt(1:nx,1:ns), errstat)
    ! call tiof_get2d_ui2 (this % ft, o3t_var_geoflg, [0,0], [ns,nx], &
    !                      rg % geoflg(1:nx,1:ns), errstat)
    call tiof_get2d_ui4 (this % ft, o3t_var_geoflg, [0,0], [ns,nx], &
                         tmp_geoflg(1:nx,1:ns), errstat)
    call tiof_get2d_ui1 (this % ft, o3t_var_anomflg, [0,0], [ns,nx], &
                        rg % anomflg(1:nx,1:ns), errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "l1b_tio_init_rad: reading geo variables", errstat)
      return
    endif

    ! FIXME?  For TEMPO, the input ground_pixel_quality_flag has 32-bits,
    !         but this code uses only the low order 16-bits,
    !         including them in the output.  For now, we're discarding
    !         the unused high-order 16 bits on input.  If we want to
    !         carry the full 32-bits through to the output, more extensive
    !         changes will have to be made.
    rg % geoflg(1:nx, 1:ns) = tmp_geoflg(1:nx, 1:ns)
    deallocate (tmp_geoflg)

  end subroutine

  !> Load a block of radiance spectra into an initialized \a l1b_radgeo_type object
  !! @param[inout] this  Level 1 file object
  !! @param[inout] rg  The \a l1b_radgeo_type object to receive the spectra
  !! @param[in]  iline  Index of first scan line to read.
  !! @param[inout] errstat  Error status variable
  !! @details
  !! Note that as each block is loaded into the cache, the
  !! \a l1b_radgeo_type object records the beginning and ending scan lines
  !! that were loaded.
  subroutine load_radiance_block (this, rg, iline, errstat)
    use netcdf, only : nf90_noerr, nf90_inq_varid ! FIXME <---
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    integer, intent(inout) :: errstat

    integer :: num_read, nx, nw, status, instid_varid

    if (errstat < 0) return

    num_read = rg % num_lines
    if ((iline + num_read) > this % num_steps) then
      num_read = this % num_steps - iline
    endif

    rg % beg_line = iline
    rg % end_line = iline + num_read - 1

    nx = this % num_xtrack
    nw = this % num_wavelengths

    ! FIXME:  This instid stuff is here only to support tests with OMI data.
    !         Maybe it will eventually go away?
    status = nf90_inq_varid (this%ft%groupid, o3t_var_instid, instid_varid)
    if (status == nf90_noerr) then
      call tiof_get1d_ui1 (this % ft, o3t_var_instid, [iline], [num_read], &
                           rg % instid (1:num_read), errstat)
    else
      rg % instid (1:num_read) = 0_1
    endif

    call tiof_get1d_ui2 (this % ft, o3t_var_mqf, [iline], [num_read], &
                        rg % measurement_quality_flags (1:num_read), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_wavelength, [iline,0,0], [num_read,nx,nw], &
                        rg % wavelength(1:nw,1:nx,1:num_read), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_radiance, [iline,0,0], [num_read,nx,nw], &
                        rg % radiance(1:nw,1:nx,1:num_read), errstat)
    call tiof_get3d_ui2 (this % ft, o3t_var_pqf, [iline,0,0], [num_read,nx,nw], &
                        rg % qa_flags(1:nw,1:nx,1:num_read), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_read_error, &
                       'load_radiance_block: reading radiance block', errstat)
      return
    endif

  end subroutine

  subroutine maybe_load_new_block (this, rg, iline, errstat)
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    integer, intent(inout) :: errstat

    character (len=256) :: msg

    if (errstat < 0) return

    if (iline < 0 .or. iline >= this % num_steps) then
      write(msg, '(a,i9,a)')'maybe_load_new_block: iline=',iline,' out of range'
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    if (iline < rg % beg_line .or. iline > rg % end_line) then
      call load_radiance_block (this, rg, iline, errstat)
      if (errstat < 0) then
        write(msg, '(a,i9)')'maybe_load_new_block: loading radiance block for iline=',iline
        call tell_error (tell_runtime_error, msg, errstat)
        return
      endif
    endif

  end subroutine

  !> Copy radiance spectra for a scan line from an \a l1b_radgeo_type object
  !! @param[inout]  this  Level 1 file object
  !! @param[inout] rg  The \a l1b_radgeo_type object
  !! @param[in]  iline  The scan line to copy.
  !! @param[inout] errstat  Error status variable
  !! @param[out] radiance  Array of radiance spectra, one for each cross-track pixel.
  !! @param[out] wavelength  Array of radiance wavelength grids, one for each cross-track pixel.
  !! @param[out] qa_flags  Array of quality assurance flags, one for each spectral bin.
  !! @details
  !! If \a iline isn't in the block cached
  !! by the \a l1b_radgeo_type object, a block of spectra at
  !! \a iline will be loaded from the open Level 1 file (\a this).
  subroutine l1b_tio_getrad (this, rg, iline, errstat, &
                             radiance, wavelength, qa_flags)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    integer, intent(inout) :: errstat
    real (kind=4), dimension(:,:), intent(out) :: radiance, wavelength
    integer (kind=2), dimension (:,:), intent(out) :: qa_flags

    integer :: nx, nw, j

    call maybe_load_new_block (this, rg, iline, errstat)
    if (errstat < 0) return

    nx = this % num_xtrack
    nw = this % num_wavelengths
    j  = iline - rg % beg_line + 1

    wavelength(1:nw, 1:nx) = rg % wavelength (1:nw, 1:nx, j)
    radiance(1:nw, 1:nx) = rg % radiance (1:nw, 1:nx, j)
    qa_flags(1:nw, 1:nx) = rg % qa_flags (1:nw, 1:nx, j)

  end subroutine

  !> Copy the measurement quality flag and instrument id for a scan line
  !! from an \a l1b_radgeo_type object
  !! @param[inout]  this  Level 1 file object
  !! @param[inout] rg  The \a l1b_radgeo_type object
  !! @param[in]  iline  The scan line to copy.
  !! @param[inout] errstat  Error status variable
  !! @param[out] instid  (Optional) Instrument ID
  !! @param[out] mqf  (Optional) Measurement quality flag
  !! @details
  !! If \a iline isn't in the block cached
  !! by the \a l1b_radgeo_type object, a block of spectra at
  !! \a iline will be loaded from the open Level 1 file (\a this).
  subroutine l1b_tio_get_etc (this, rg, iline, errstat, instid, mqf)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    integer (kind=1), optional, intent(out) :: instid
    integer (kind=2), optional, intent(out) :: mqf
    integer, intent(inout) :: errstat

    integer :: j

    call maybe_load_new_block (this, rg, iline, errstat)
    if (errstat < 0) return

    j  = iline - rg % beg_line + 1
    if (present(instid)) instid = rg % instid (j)
    if (present(mqf)) mqf = rg % measurement_quality_flags (j)

  end subroutine

  !> Copy geolocation data for a scan line from an \a l1b_radgeo_type object
  !! @param[inout]  this  Level 1 file object
  !! @param[inout] rg  The \a l1b_radgeo_type object
  !! @param[in]  iline  The scan line to copy.
  !! @param[out] lat  Latitude for each cross-track pixel
  !! @param[out] lon  Longitude for each cross-track pixel
  !! @param[out] lat_bounds  Latitude corners for each cross-track pixel
  !! @param[out] lon_bounds  Longitude corners for each cross-track pixel
  !! @param[out] step_index  Mirror step index
  !! @param[out] sza  Solar zenith angle for each cross-track pixel
  !! @param[out] saz  Solar azimuth angle for each cross-track pixel
  !! @param[out] vza  Viewing zenith angle for each cross-track pixel
  !! @param[out] vaz  Viewing azimuth angle for each cross-track pixel
  !! @param[out] height  Terrain height for each cross-track pixel
  !! @param[out] geoflg  Geolocation flag for each cross-track pixel
  !! @param[inout] errstat  Error status variable
  !! @param[out] anomflg  (Optional) Anomaly flag for each cross-track pixel
  !! @details
  !! If \a iline isn't in the block cached
  !! by the \a l1b_radgeo_type object, a block of spectra at
  !! \a iline will be loaded from the open Level 1 file (\a this).
  subroutine l1b_tio_getgeo (this, rg, iline, lat, lon, &
                             lat_bounds, lon_bounds, step_index, &
                             sza, saz, vza, vaz, height, geoflg, errstat, &
                             anomflg)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    real (kind=4), dimension(:), intent(out) :: lat, lon, sza, saz, vza, vaz
    real (kind=4), dimension(:,:), intent(out) :: lat_bounds, lon_bounds
    integer (kind=4) :: step_index
    integer (kind=2), dimension(:), intent(out) :: height, geoflg
    integer (kind=1), dimension(:), intent(out), optional :: anomflg
    integer, intent(inout) :: errstat

    character (len=256) :: msg
    integer :: nx, i

    if (errstat < 0) return

    if (iline < 0 .or. iline >= this % num_steps) then
      write(msg, '(a,i9,a)')'l1b_tio_getgeo: iline=',iline,' out of range'
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    nx = this % num_xtrack
    i  = iline + 1

    lat(1:nx) = rg % lat(1:nx, i)
    lon(1:nx) = rg % lon(1:nx, i)
    lat_bounds(1:4,1:nx) = rg % lat_bounds(1:4, 1:nx, i)
    lon_bounds(1:4,1:nx) = rg % lon_bounds(1:4, 1:nx, i)
    step_index = rg % step_indices (i)
    sza(1:nx) = rg % sza(1:nx, i)
    saz(1:nx) = rg % saz(1:nx, i)
    vza(1:nx) = rg % vza(1:nx, i)
    vaz(1:nx) = rg % vaz(1:nx, i)
    height(1:nx) = rg % hgt (1:nx, i)
    geoflg(1:nx) = rg % geoflg(1:nx, i)

    if (present(anomflg)) then
      anomflg(1:nx) = rg % anomflg(1:nx, i)
    endif

  end subroutine

  !> Read selected metadata attributes from a Level 1 product file
  !! @param[inout] obj   Level 1 file object
  !! @param[out]  m     Metadata object
  !! @param[inout] errstat  Error status variable
  subroutine l1b_tio_get_metadata (obj, m, errstat)
    use netcdf
    use iso_c_binding, only : c_null_char
    implicit none
    type (l1b_tio_type), intent(inout) :: obj
    type (l1b_metadata_type), intent(out) :: m
    integer, intent(inout) :: errstat

    character (len=64) :: tmpstr
    integer :: ncid, ierr, status, i_null
    integer (kind=4), external :: day_of_year

    if (errstat < 0) return

    ncid = obj % ft % fileid

    status = nf90_inquire_attribute (ncid, nf90_global, "omi_shortname")
    if (status == nf90_noerr) then
      status = nf90_get_att (ncid, nf90_global, "omi_shortname", tmpstr)
      if (status /= nf90_noerr) then
        call tell_error (tell_io_read_error, &
                         "reading 'omi_shortname' from file: "//trim(obj % filename), &
                         errstat)
        return
      endif
      i_null = index (tmpstr, c_null_char)
      if (i_null > 0) tmpstr = tmpstr(1:i_null-1)
      m % shortname = tmpstr(1:len(m%shortname))
    else
      m % shortname = ''
    endif

    status = nf90_inquire_attribute (ncid, nf90_global, "omi_orbit_number")
    if (status == nf90_noerr) then
      status = nf90_get_att (ncid, nf90_global, "omi_orbit_number", m % orbit_number)
      if (status /= nf90_noerr) then
        call tell_error (tell_io_read_error, &
                         "reading 'omi_orbit_number' from file: "//trim(obj % filename), &
                       errstat)
        return
      endif
    else
      m % orbit_number = -1
    endif

    status = nf90_get_att (ncid, nf90_global, "time_coverage_start", tmpstr)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
                       "reading 'time_coverage_start' from file: "//trim(obj % filename), &
                       errstat)
      return
    endif

    read (tmpstr, '(i4,1x,i2,1x,i2)', iostat=ierr) m % year, m % month, m % day
    if (ierr /= 0) then
      call tell_error (tell_runtime_error, &
                       "parsing time_coverage_start='"//trim(tmpstr)//"'", &
                       errstat)
      return
    endif
    m % jday = day_of_year (m % year, m % month, m % day)

  end subroutine

end module
