module gler_module
  use netcdf, only : nf90_nowrite
  use tell_module
  use tio_module
  use ndinterp_module
  use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_double, &
    c_loc, c_ptr, c_null_ptr, c_associated
  implicit none

  private

  public gler_open, gler_close, gler_interp_time, gler_albedo

  integer, private, parameter :: path_bufsize = 1024, msg_bufsize = 128
  real (kind=4), private, parameter :: albedo_fill_value = -1.e30

  ! Mask values
  integer (kind=2), private, parameter :: m_land=0, m_water=1, m_shadow=2

  ! NOTE: The ocean (lon,lat) grid is much coarser than on land.

  type :: gler_land_type
    type (ndi_dim), dimension(2) :: dims   ! dims = [lat, lon]
    real (kind=4), allocatable, dimension(:,:) :: albedo
    integer (kind=2), allocatable, dimension(:,:) :: mask
  end type

  type :: gler_ocean_type
    type (ndi_dim), dimension(3) :: dims   ! dims = [windspeed, lat, lon]
    real (kind=4), allocatable, dimension(:,:,:) :: albedo
  end type

  type, public :: gler_type
    private
    type(c_ptr) :: c_gler_type
    type(gler_land_type) :: land
    type(gler_ocean_type) :: ocean
  end type

  interface
    function c_gler_open (iwave, cfg_file_ptr) bind (c, name='gler_open')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int
      implicit none
      integer (kind=c_int), value :: iwave
      type (c_ptr), value :: cfg_file_ptr
      type (c_ptr) :: c_gler_open
    end function
  end interface

  interface
    subroutine c_gler_close (ptr) bind (c, name='gler_close')
      use, intrinsic :: iso_c_binding, only: c_ptr
      implicit none
      type (c_ptr), value :: ptr
    end subroutine
  end interface

  interface
    function c_gler_land_lookup (ptr, taix, a, b, awt) &
        bind (c, name='gler_land_lookup')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int
      implicit none
      type (c_ptr), value :: ptr
      real (kind=c_double), value :: taix
      integer (c_int), intent(out) :: a, b
      real (kind=c_double), intent(out) :: awt
      integer (c_int) :: c_gler_land_lookup
    end function
  end interface

  interface
    function c_gler_land_file (ptr, k, path, pathlen) &
        bind (c, name='gler_land_file')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
      implicit none
      type (c_ptr), value :: ptr
      integer (c_int), value, intent(in) :: k, pathlen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_gler_land_file
    end function
  end interface

  interface
    function c_gler_ocean_lookup (ptr, taix, a, b, awt) &
        bind (c, name='gler_ocean_lookup')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int
      implicit none
      type (c_ptr), value :: ptr
      real (kind=c_double), value :: taix
      integer (c_int), intent(out) :: a, b
      real (kind=c_double), intent(out) :: awt
      integer (c_int) :: c_gler_ocean_lookup
    end function
  end interface

  interface
    function c_gler_ocean_file (ptr, k, path, pathlen) &
        bind (c, name='gler_ocean_file')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
      implicit none
      type (c_ptr), value :: ptr
      integer (c_int), value, intent(in) :: k, pathlen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_gler_ocean_file
    end function
  end interface

contains

  ! read land GLER file and interpolate albedo for a specified time of day
  subroutine gler_read_land (gl, hourf_in, path, errstat)
    implicit none
    type (gler_land_type), intent(inout) :: gl
    real (kind=8), intent(in) :: hourf_in
    character (kind=c_char,len=*), intent(in) :: path
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj
    integer :: dim_times, dim_lon, dim_lat, i0, lon, lat
    integer (kind=2) :: m1, m2, mm
    real (kind=4) :: a1, a2
    real (kind=8) :: wt0, hourf, alb
    real (kind=8), allocatable, dimension(:) :: times
    real (kind=4), allocatable, dimension(:,:,:) :: albedo
    integer (kind=2), allocatable, dimension(:,:,:) :: mask

    if (errstat /= 0) return

    hourf = mod(hourf_in, 24.0)

    call tiof_open (path, obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (obj, 't', dim_times, errstat)
    call tiof_inq_dimlen (obj, 'x', dim_lon, errstat)
    call tiof_inq_dimlen (obj, 'y', dim_lat, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading file: '//trim(path), errstat)
      return
    endif

    allocate (times(dim_times))

    call tiof_get1d_r8 (obj, 'hour', (/0/), (/dim_times/), times, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading times from file: '//trim(path), errstat)
      return
    endif

    ! Find indices that bracket the specified time-of-day, hourf
    if (hourf < times(1) .or. times(dim_times) < hourf) hourf = times(1)
    i0 = binary_search (times, hourf)
    wt0 = (times(i0+1) - hourf) / (times(i0+1) - times(i0))

    allocate (gl % albedo (dim_lat, dim_lon), &
              gl % mask (dim_lat, dim_lon))
    allocate (albedo (dim_lat,dim_lon,2), &
              mask (dim_lat, dim_lon, 2))

    call ndi_dims_alloc (gl % dims, (/dim_lat, dim_lon/), errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, 'allocating dimension structures', errstat)
      return
    endif

    call tiof_get1d_r8 (obj, 'lat', (/0/), (/dim_lat/), gl % dims(1) % x, errstat)
    call tiof_get1d_r8 (obj, 'lon', (/0/), (/dim_lon/), gl % dims(2) % x, errstat)
    call tiof_get3d_i2 (obj, 'qf', (/i0-1,0,0/), (/2,dim_lon,dim_lat/), mask, errstat)

    ! Read albedo tables only for the relevant times
    call tiof_get3d_r4 (obj, 'alb', (/i0-1,0,0/), (/2,dim_lon,dim_lat/), albedo, errstat, &
                       replace_fill=albedo_fill_value)
    call tiof_close (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading file: '//trim(path), errstat)
      return
    endif

    ! Interpolation to specified time-of-day
    do lon = 1, gl % dims(2) % dimlen
      do lat = 1, gl % dims(1) % dimlen

        a1 = albedo (lat,lon,1)
        a2 = albedo (lat,lon,2)
        m1 = mask (lat,lon,1)
        m2 = mask (lat,lon,2)

        ! assume land, and change below when necessary
        mm = m_land

        ! The land/water pattern is identical for all times.
        ! Only the shadow pattern changes with time.
        if (m1 == m_water .or. m2 == m_water) then
          alb = albedo_fill_value
          mm = m_water
        else if (m1 == m_land) then
          if (m2 == m_land) then
            alb = wt0 * a1 + (1.0 - wt0) * a2
          else  ! m2=water is not possible, so must be m2=shadow
            alb = a1
          endif
        else if (m2 == m_land) then  ! m1 must be shadow
          alb = a2
        else      ! case m2=water already checked, so m2 must be shadow
          alb = 0 ! both in shadow
          mm = m_shadow
        endif
        gl % albedo(lat,lon) = real (alb, kind=4)
        gl % mask(lat,lon) = mm
      enddo
    enddo

  end subroutine gler_read_land

  ! read ocean GLER file and interpolate albedo for a specified time of day
  subroutine gler_read_ocean (gl, hourf_in, path, errstat)
    implicit none
    type (gler_ocean_type), intent(inout) :: gl
    real (kind=8), intent(in) :: hourf_in
    character (kind=c_char,len=*), intent(in) :: path
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj
    integer :: dim_times, dim_ws, dim_lon, dim_lat, i0
    integer :: k, lon, lat
    real (kind=8) :: wt0, hourf
    real (kind=4) :: wt0f, alb, a1, a2
    real (kind=8), allocatable, dimension(:) :: times
    real (kind=4), allocatable, dimension(:,:,:,:) :: albedo

    if (errstat /= 0) return

    hourf = mod(hourf_in, 24.0)

    call tiof_open (path, obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (obj, 't', dim_times, errstat)
    call tiof_inq_dimlen (obj, 'x', dim_lon, errstat)
    call tiof_inq_dimlen (obj, 'y', dim_lat, errstat)
    call tiof_inq_dimlen (obj, 'ws', dim_ws, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading file: '//trim(path), errstat)
      return
    endif

    allocate (times(dim_times))

    call tiof_get1d_r8 (obj, 'hour', (/0/), (/dim_times/), times, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading times from file: '//trim(path), errstat)
      return
    endif

    ! Find indices that bracket the specified time-of-day, hourf
    if (hourf < times(1) .or. times(dim_times) < hourf) hourf = times(1)
    i0 = binary_search (times, hourf)
    wt0 = (times(i0+1) - hourf) / (times(i0+1) - times(i0))
    wt0f = real(wt0,kind=4)

    allocate (gl % albedo (dim_ws, dim_lat, dim_lon))
    allocate (albedo (dim_ws,dim_lat,dim_lon,2))

    call ndi_dims_alloc (gl % dims, (/dim_ws, dim_lat, dim_lon/), errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, 'allocating dimension structures', errstat)
      return
    endif

    call tiof_get1d_r8 (obj, 'wspd', (/0/), (/dim_ws/), gl % dims(1) % x, errstat)
    call tiof_get1d_r8 (obj, 'lat', (/0/), (/dim_lat/), gl % dims(2) % x, errstat)
    call tiof_get1d_r8 (obj, 'lon', (/0/), (/dim_lon/), gl % dims(3) % x, errstat)

    ! Read albedo tables only for the relevant times
    ! Arrays start,count use C index ordering;  Declare arrays using Fortran index ordering (obviously)
    call tiof_get4d_r4 (obj, 'alb', (/i0-1,0,0,0/), (/2,dim_lon,dim_lat,dim_ws/), albedo, errstat, &
                       replace_fill=albedo_fill_value)
    call tiof_close (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading file: '//trim(path), errstat)
      return
    endif

    ! The "ocean" albedo table apparently has albedo values for locations that
    ! are not on the ocean, and also fill values for pixels that may or may not
    ! be on the ocean.

    ! Interpolation to specified time-of-day
    do lon = 1, gl % dims(3) % dimlen
      do lat = 1, gl % dims(2) % dimlen
        do k = 1, gl % dims(1) % dimlen
          a1 = albedo(k,lat,lon,1)
          a2 = albedo(k,lat,lon,2)
          if (a1 /= albedo_fill_value .and. a2 /= albedo_fill_value) then
            alb = wt0f * a1 + (1.0 - wt0f) * a2
          else if (a1 == albedo_fill_value) then
            alb = a2
          else if (a2 == albedo_fill_value) then
            alb = a1
          else
            alb = albedo_fill_value
          endif
          gl % albedo (k,lat,lon) = alb
        enddo
      enddo
    enddo

  end subroutine gler_read_ocean

  subroutine gler_close (glt)
    implicit none
    type (gler_type), intent(inout) :: glt
    call c_gler_close (glt % c_gler_type)
  end subroutine

  subroutine gler_open (glt, iwave, errstat, config_file)
    implicit none
    type (gler_type), intent(inout) :: glt
    integer, intent(in) :: iwave
    character (len=*), optional, intent(in) :: config_file
    integer, intent(inout) :: errstat

    character (kind=c_char, len=:), allocatable, target:: c_config_file

    if (present (config_file)) then
      allocate (character (len=len(config_file)) :: c_config_file)
      c_config_file = trim(config_file) // c_null_char
      glt % c_gler_type = c_gler_open (iwave, c_loc(c_config_file))
      deallocate (c_config_file)
    else
      glt % c_gler_type = c_gler_open (iwave, c_null_ptr)
    endif
    if (.not.c_associated(glt % c_gler_type)) then
      call tell_error (tell_runtime_error, 'gler_open: failed', errstat)
      return
    endif
  end subroutine gler_open

  subroutine gler_interp_time (glt, taix, errstat)
    implicit none
    type (gler_type), intent(inout) :: glt
    real (kind=8), intent(in) :: taix
    integer, intent(inout) :: errstat

    integer :: i0, i1, k, year, month, day, status
    character (kind=c_char, len=path_bufsize) :: path0, path1
    real (kind=8) :: hourf, wt0_yday
    real (kind=4) :: wt0f
    type (gler_land_type) :: tmp_lt
    type (gler_ocean_type) :: tmp_ot

    if (errstat /= 0) return

    call tiof_taix_time_to_utc_caldate (taix, year, month, day, hourf, errstat)

    ! Land GLER =============================================

    ! Which files bracket the specified day-of-year? */
    status = c_gler_land_lookup (glt % c_gler_type, taix, i0, i1, wt0_yday)
    if (status /= 0) then
      call tell_set_error (tell_runtime_error)
      return
    endif
    wt0f = real(wt0_yday,kind=4)

    ! Retrieve the path to file i0
    status = c_gler_land_file (glt % c_gler_type, i0, path0, path_bufsize)
    if (status /= 0) then
      call tell_set_error (tell_runtime_error)
      return
    endif

    ! Retrieve the path to file i1
    status = c_gler_land_file (glt % c_gler_type, i1, path1, path_bufsize)
    if (status /= 0) then
      call tell_set_error (tell_runtime_error)
      return
    endif

    ! Read the two files, interpolating to the relevant measurement time-of-day, hourf
    call gler_read_land (glt % land, hourf, path0, errstat)
    call gler_read_land (tmp_lt, hourf, path1, errstat)
    if (errstat /= 0) return

    ! Perform day-of-year interpolation
    where (glt % land % mask /= m_water)
      glt % land % albedo(:,:) = (wt0f * glt % land % albedo(:,:) + &
                                  (1.0-wt0f) * tmp_lt % albedo(:,:))
    endwhere

    ! Ocean GLER =============================================

    ! Retrieve the path to file i0
    status = c_gler_ocean_file (glt % c_gler_type, i0, path0, path_bufsize)
    if (status /= 0) then
      call tell_set_error (tell_runtime_error)
      return
    endif

    ! Retrieve the path to file i1
    status = c_gler_ocean_file (glt % c_gler_type, i1, path1, path_bufsize)
    if (status /= 0) then
      call tell_set_error (tell_runtime_error)
      return
    endif

    ! Read the two files, interpolating to the relevant measurement time-of-day, hourf
    call gler_read_ocean (glt % ocean, hourf, path0, errstat)
    call gler_read_ocean (tmp_ot, hourf, path1, errstat)
    if (errstat /= 0) return

    ! Perform day-of-year interpolation
    do k = 1, glt % ocean % dims(1) % dimlen
      where (glt % ocean % albedo(k,:,:) /= albedo_fill_value)
        glt % ocean % albedo(k,:,:) = (wt0f * glt % ocean % albedo(k,:,:) + &
                                       (1.0-wt0f) * tmp_ot % albedo(k,:,:))
      endwhere
    enddo

  end subroutine

  integer function binary_search (xa, t)
    implicit none
    real (kind=8), dimension(:), intent(in) :: xa
    real (kind=8), intent(in) :: t

    integer :: n, n0, n1, n2
    real (kind=8) :: xt

    n = size(xa)

    n0 = 1
    n1 = n+1

    ! Require xa in ascending order
    if (xa(1) >= xa(2)) then
      write (*,*)'*** Error: binary_search requires xa in ascending order'
      binary_search = -1
      return
    endif

    ! Don't extrapolate
    if (t < xa(1)) then
      binary_search = -1
      return
    else if (xa(n) < t) then
      binary_search = -n
      return
    endif

    do while (n1 > n0+1)
      n2 = (n0 + n1) / 2
      xt = xa(n2)
      if (t <= xt) then
        if (xt == t) then
          binary_search = n2
          return
        endif
        n1 = n2
      else
        n0 = n2
      endif
    enddo

    binary_search = n0

  end function

  subroutine gler_albedo (glt, lon, lat, wind_speed, alb, errstat)
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    implicit none
    type (gler_type), intent(in) :: glt
    real (kind=4), intent(in) :: lon, lat, wind_speed
    real (kind=4), intent(out) :: alb
    integer, intent(inout) :: errstat

    integer :: ilat, ilon, nw
    real (kind=8), dimension(2) :: x_land
    real (kind=8), dimension(3) :: x_ocean
    real (kind=8) :: wmin, wmax
    real (kind=4) :: wind_speed_trim

    if (errstat /= 0) return

    ilat = binary_search (glt % land % dims(1) % x, real(lat,kind=8))
    ilon = binary_search (glt % land % dims(2) % x, real(lon,kind=8))

    ! don't extrapolate
    if (ilat < 0 .or. ilon < 0) then
      alb = ieee_value (alb, ieee_quiet_nan)
      return
    endif

    if (glt % land % mask(ilat, ilon) == m_land) then
      x_land(:) = real ((/lat, lon/), kind=8)
      call ndi_table_interp (glt % land % dims, &
                             glt % land % albedo, &
                             x_land, alb, errstat, fill_value=albedo_fill_value)
    else if (glt % land % mask(ilat, ilon) == m_water) then
      ! Force the wind speed coordinate to be on the tabulated grid
      nw = size (glt % ocean % dims(1) % x)
      wmax = glt % ocean % dims(1) % x(nw)
      wmin = glt % ocean % dims(1) % x(1)
      wind_speed_trim = wind_speed
      if (wind_speed < wmin) then
        wind_speed_trim = real (wmin, kind=4)
      else if (wind_speed > wmax) then
        wind_speed_trim = real (wmax, kind=4)
      endif

      x_ocean(:) = real ((/wind_speed_trim, lat, lon/), kind=8)
      call ndi_table_interp (glt % ocean % dims, &
                             glt % ocean % albedo, &
                             x_ocean, alb, errstat, fill_value=albedo_fill_value)
    else
      alb = 0.0
    endif

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "gler_albedo: interpolation failed", errstat)
      return
    endif

  end subroutine

end module
