module l1b_tio_class
  use netcdf
  use tio_module
  use tell_module
  use o3t_names_module
  implicit none
  private

  integer, parameter :: nblock = 10 ! default number of lines in a block

  ! opaque type
  type, public :: l1b_tio_type
    private
    type (tiof_file_type) :: ft
    character (len=1024) :: filename
    character (len=tiof_max_name_len) :: groupname
    integer :: num_steps
    integer :: num_xtrack
    integer :: num_wavelengths
    real (kind=4) :: distance   ! earth-sun distance
  end type

  type, public :: l1b_radgeo_type
    real (kind=4), dimension(:,:), allocatable :: lon, lat, sza, saz, vza, vaz
    real (kind=4), dimension(:,:,:), allocatable :: radiance, wavelength
    integer (kind=2), dimension(:,:,:), allocatable :: qa_flags
    integer (kind=2), dimension(:,:), allocatable :: geoflg, hgt
    integer (kind=1), dimension(:,:), allocatable :: anomflg
    integer (kind=2), dimension(:), allocatable :: measurement_quality_flags
    integer (kind=1), dimension(:), allocatable :: instid
    integer :: beg_line, end_line, num_lines
  end type

  public l1b_tio_open, l1b_tio_close, l1b_tio_getdims, l1b_tio_get_irr, &
    l1b_tio_earthsun_distance, l1b_tio_init_rad, l1b_tio_getrad, &
    l1b_tio_get_etc, l1b_tio_getgeo

contains

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

  subroutine l1b_tio_close (this, errstat)
    implicit none
    type(l1b_tio_type), intent(inout) :: this
    integer, intent(inout) :: errstat
    if (errstat < 0) continue  ! no-op
    this % num_steps = -1
    this % num_xtrack = -1
    this % num_wavelengths = -1
    this % filename = ''
    this % groupname = ''
  end subroutine

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

  subroutine l1b_tio_earthsun_distance (this, distance, errstat)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    real (kind=4), intent(out) :: distance
    integer, intent(inout) :: errstat
    if (errstat < 0) return
    distance = this % distance
  end subroutine

  subroutine l1b_tio_get_irr (this, irr, wavelength, qaflags, &
                              measurement_quality_flag, instid, errstat)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    real (kind=4), dimension(:,:), intent(out) :: irr, wavelength
    integer (kind=2), dimension(:,:), intent(out) :: qaflags
    integer (kind=2), intent(out) :: measurement_quality_flag
    integer (kind=1), intent(out) :: instid
    integer, intent(inout) :: errstat

    integer :: nx, nw, status
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

    call tiof_get3d_r4 (this % ft, o3t_var_wavelength, [0,0,0], [1,nx,nw], &
                        tmp_wavelengths(:,1:nx,1:1), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_irradiance, [0,0,0], [1,nx,nw], &
                        tmp_spectrum(:,1:nx,1:1), errstat)
    call tiof_get3d_ui2 (this % ft, o3t_var_pqf, [0,0,0], [1,nx,nw], &
                        tmp_qflags(:,1:nx,1:1), errstat)
    call tiof_get1d_ui2 (this % ft, o3t_var_mqf, [0], [1], tmp_mqf(1:1), errstat)
    call tiof_get1d_ui1 (this % ft, o3t_var_instid, [0], [1], tmp_instid(1:1), errstat)
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
    instid = tmp_instid(1)

  end subroutine

  subroutine l1b_tio_init_rad (this, rg, errstat)
    implicit none
    type(l1b_tio_type), intent(inout), target :: this
    type(l1b_radgeo_type), intent(inout) :: rg
    integer, intent(inout) :: errstat

    integer :: nx, ns, nw, nl, ierr

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
              rg % sza (nx, ns), &
              rg % saz (nx, ns), &
              rg % vza (nx, ns), &
              rg % vaz (nx, ns), &
              rg % hgt (nx, ns), &
              rg % geoflg (nx, ns), &
              rg % anomflg (nx, ns), &
              rg % radiance (nw, nx, nl), &
              rg % wavelength (nw, nx, nl), &
              rg % qa_flags (nw, nx, nl), &
              rg % measurement_quality_flags (nl), &
              rg % instid (nl), stat = ierr)
    if (ierr /= 0) then
      call tell_error (tell_malloc_error, "l1b_tio_init_rad: allocate failed", errstat)
      return
    endif

    call tiof_get2d_r4 (this % ft, o3t_var_longitude, [0,0], [ns,nx], &
                        rg % lon(1:nx,1:ns), errstat)
    call tiof_get2d_r4 (this % ft, o3t_var_latitude, [0,0], [ns,nx], &
                        rg % lat(1:nx,1:ns), errstat)
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
    call tiof_get2d_ui2 (this % ft, o3t_var_geoflg, [0,0], [ns,nx], &
                         rg % geoflg(1:nx,1:ns), errstat)
    call tiof_get2d_ui1 (this % ft, o3t_var_anomflg, [0,0], [ns,nx], &
                        rg % anomflg(1:nx,1:ns), errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "l1b_tio_init_rad: reading geo variables", errstat)
      return
    endif

  end subroutine

  subroutine load_radiance_block (this, rg, iline, errstat)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    integer, intent(inout) :: errstat

    integer :: num_read, nx, nw

    if (errstat < 0) return

    num_read = rg % num_lines
    if ((iline + num_read) > this % num_steps) then
      num_read = this % num_steps - iline
    endif

    rg % beg_line = iline
    rg % end_line = iline + num_read - 1

    nx = this % num_xtrack
    nw = this % num_wavelengths

    call tiof_get1d_ui2 (this % ft, o3t_var_mqf, [iline], [num_read], &
                        rg % measurement_quality_flags (1:num_read), errstat)
    call tiof_get1d_ui1 (this % ft, o3t_var_instid, [iline], [num_read], &
                        rg % instid (1:num_read), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_wavelength, [iline,0,0], [num_read,nx,nw], &
                        rg % wavelength(1:nw,1:nx,1:num_read), errstat)
    call tiof_get3d_r4 (this % ft, o3t_var_radiance, [iline,0,0], [num_read,nx,nw], &
                        rg % radiance(1:nw,1:nx,1:num_read), errstat)
    call tiof_get3d_ui2 (this % ft, o3t_var_dqf, [iline,0,0], [num_read,nx,nw], &
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
  
  subroutine l1b_tio_getgeo (this, rg, iline, lat, lon, &
                             sza, saz, vza, vaz, height, geoflg, errstat, &
                             anomflg)
    implicit none
    type (l1b_tio_type), intent(inout) :: this
    type (l1b_radgeo_type), intent(inout) :: rg
    integer, intent(in) :: iline
    real (kind=4), dimension(:), intent(out) :: lat, lon, sza, saz, vza, vaz
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

end module
