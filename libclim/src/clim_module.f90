!> Climatology interface  module
!! @file
!! @sa  clim_util.c   (Utility functions implemented in C)
!!
!! This module provides an interface to a trace gas climatology database.
!!
!! The trace gas climatology files provide the volume mixing ratio (VMR) of
!! selected trace gases.  For each trace gas species, there is one climatology
!! file for each hour of the day, for each month, for a total of 24 * 12 = 288
!! climatology files.  Each file provides VMR(z,lat,lon).
!! The z coordinate is defined in terms of pressure.
!! Each file also contains the surface pressure, PS (lat,lon,hour)
!! The pressure vs height, P(z,lat,lon,hour), is obtained using the
!! EtaA(z) and EtaB(z) variables from a separate file:
!!   P(z,lat,lon,hour) = EtaA(z) + EtaB(z) * PS(lat,lon,hour)
!!
module clim_module
  use netcdf, only : nf90_nowrite
  use tell_module
  use tio_module
  implicit none

  private

  public clim_query_nz
  public clim_pres, clim_pres_init, clim_pres_eta
  public clim_val_interp, clim_val_init
  public clim_cloud, clim_cloud_init
  public clim_partial_column

  interface clim_pres_init
    module procedure clim_pres_init_args
    module procedure clim_pres_init_struct
  end interface

  real (kind=4), private, parameter :: mean_days_per_month = 365.25/12

  integer, private, parameter :: path_bufsize = 1024, msg_bufsize = 128

  ! Vertical grid has Num_Layers with Num_Layers+1 edges
  ! At each layer edge:     EtaA, EtaB, i=1,Num_Layers+1
  ! At each layer midpoint: have VMR(k), k=1,Num_Layers
  real (kind=4), private, dimension(:), allocatable :: EtaA, EtaB
  integer, private :: Num_Layers
  logical, private :: Have_Forecast

  type :: dim_subset_type
    integer :: dimlen
    real (kind=4) :: min, max    ! domain limits [min, max)
    integer ::      imin, imax   ! file variable array index limits implied by [min, max)
    integer :: num_values        ! = imax-imin+1

    real (kind=4), dimension(:), allocatable :: values  ! (num_values)
  end type

  type :: clim_pres_slab_type
    real (kind=4), dimension(:,:,:), allocatable :: p_surf    ! (nlat,nlon,nhour)
    real (kind=4), dimension(:,:,:), allocatable :: p_trop    ! (nlat,nlon,nhour)
  end type

  type, public :: clim_pres_type
    private
    integer :: year, month, day
    integer :: month0, month1
    real (kind=4) :: month0_weight
    type(dim_subset_type) :: lon_subset  ! nlon, grid box centers
    type(dim_subset_type) :: lat_subset  ! nlat, grid box centers
    integer :: nhours
    real (kind=4), dimension(:), allocatable :: hours         ! (nhour)
    real (kind=4), dimension(:,:,:), allocatable :: p_surf    ! (nlat,nlon,nhour)
    real (kind=4), dimension(:,:,:), allocatable :: p_trop    ! (nlat,nlon,nhour)
  end type

  type, public :: clim_pres_bounds_type
    real (kind=4) :: hour_beg, hour_end !< UTC begin/end of time interval of interest
    real (kind=4) :: lon_min, lon_max   !< longitude range of interest [deg]
    real (kind=4) :: lat_min, lat_max   !< latitude range of interest [deg]
  end type

  type :: clim_type
    integer :: nz, nlat, nlon
    real (kind=4), dimension(:,:,:), allocatable :: values   ! (nz, nlat, nlon)
  end type

  type, public :: clim_val_type
    private
    integer :: nhours
    real (kind=4), dimension (:), allocatable :: hours   ! (nhours)
    type(clim_type), dimension(:), allocatable :: clim   ! (nhours)
  end type

  type, public :: clim_cloud_type
    private
    integer :: nlon, nlat, nmon
    real (kind=4), dimension(:,:,:), allocatable :: cloud_pressure  ! (nlon,nlat,nmon)
    real (kind=4), dimension(:), allocatable :: lon   ! (nlon)
    real (kind=4), dimension(:), allocatable :: lat   ! (nlat)
  end type

  interface
    function c_make_climatology_path (month, hour, path, path_buflen) &
        bind (c, name='make_climatology_path')
      use, intrinsic :: iso_c_binding, only : c_char, c_int
      implicit none
      integer (c_int), value, intent(in) :: month, hour, path_buflen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_make_climatology_path
    end function
  end interface

  interface
    function c_make_forecast_path (timet, path, path_buflen) &
        bind (c, name='make_forecast_path')
      use, intrinsic :: iso_c_binding, only : c_char, c_int, c_long
      implicit none
      integer (c_long), value, intent(in) :: timet
      integer (c_int), value, intent(in) :: path_buflen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_make_forecast_path
    end function
  end interface

  interface
    function c_make_pressure_eta_path (path, path_buflen) &
        bind (c, name='make_pressure_eta_path')
      use, intrinsic :: iso_c_binding, only : c_char, c_int
      implicit none
      integer (c_int), value, intent(in) :: path_buflen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_make_pressure_eta_path
    end function
  end interface

  interface
    function c_make_cloud_climatology_path (path, path_buflen) &
        bind (c, name='make_cloud_climatology_path')
      use, intrinsic :: iso_c_binding, only : c_char, c_int
      implicit none
      integer (c_int), value, intent(in) :: path_buflen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_make_cloud_climatology_path
    end function
  end interface

  interface
    function c_make_timet (year, month, day, hour) &
        bind (c, name='make_timet')
      use, intrinsic :: iso_c_binding, only : c_int, c_long
      implicit none
      integer (c_int), value, intent(in) :: year, month, day, hour
      integer (c_long) :: c_make_timet
    end function
  end interface

  interface
    function c_have_forecast_files (tt, num_hours) &
        bind (c, name='have_forecast_files')
      use, intrinsic :: iso_c_binding, only : c_int, c_long
      implicit none
      integer (c_long), value, intent(in) :: tt
      integer (c_int), value, intent(in) :: num_hours
      integer (c_int) :: c_have_forecast_files
    end function
  end interface

contains

  subroutine define_month_interp (cpt, month, day)
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(in) :: month, day

    real (kind=4) :: f

    cpt % month = month
    cpt % day = day

    ! At mid-month, month0 is weighted 100%.
    ! At each month0/month1 transition, both are weighted equally,
    ! with linear interpolation on all other days.

    f = (day - 1)/mean_days_per_month

    cpt % month0 = month

    if (f < 0.5) then
      cpt % month0_weight = 0.5 + f
      cpt % month1 = month - 1
      if (cpt % month1 < 1) cpt % month1 = 12
    else
      cpt % month0_weight = 1.5 - f
      cpt % month1 = month + 1
      if (cpt % month1 > 12) cpt % month1 = 1
    endif

  end subroutine

  subroutine make_pressure_eta_path (file, errstat)
    use, intrinsic :: iso_c_binding, only : c_char, c_int
    implicit none
    character (kind=c_char, len=*), target, intent(inout) :: file
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    err = c_make_pressure_eta_path (file, len(file))
    if (err /= 0) then
      call tell_error (tell_runtime_error, "make_pressure_eta_path: failed", errstat)
      return
    endif
  end subroutine

  subroutine read_pressure_eta (errstat)
    use, intrinsic :: iso_c_binding, only : c_char
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj
    character (kind=c_char, len=path_bufsize) :: file_eta
    integer :: num_edges, err

    if (errstat /= 0) return
    if (allocated (EtaA)) return

    call make_pressure_eta_path (file_eta, errstat)
    if (errstat /= 0) return

    call tiof_open (file_eta, obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (obj, 'z', num_edges, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading file: '//trim(file_eta), errstat)
      return
    endif

    ! Define global
    Num_Layers = num_edges-1

    allocate (EtaA(num_edges), &
              EtaB(num_edges), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'read_pressure_eta: allocate failed', errstat)
      return
    endif

    call tiof_get1d_r4 (obj, 'Ap', (/0/), (/num_edges/), EtaA, errstat)
    call tiof_get1d_r4 (obj, 'Bp', (/0/), (/num_edges/), EtaB, errstat)
    call tiof_close (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading pressure eta tables', errstat)
      return
    endif

  end subroutine

  subroutine subset_dim (obj, var_name, dim_name, dim_min, dim_max, &
                         dim_subset, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    character (len=*), intent(in) :: var_name, dim_name
    real (kind=4), intent(in) :: dim_min, dim_max
    type(dim_subset_type), intent(inout) :: dim_subset
    integer, intent(inout) :: errstat

    integer :: n, lo, hi, err
    real (kind=4), dimension(:), allocatable :: val

    if (errstat /= 0) return

    call tiof_inq_dimlen (obj, dim_name, n, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'subset_dim: read failed: '//trim(dim_name), errstat)
      return
    endif

    allocate (val(n), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'subset_dim: allocate failed', errstat)
      return
    endif

    call tiof_get1d_r4 (obj, var_name, (/0/), (/n/), val, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'subset_dim: read failed: '//trim(var_name), errstat)
      return
    endif

    do lo = 1, n
      if (dim_min < val(lo)) exit
    enddo
    lo = max(lo-1,1)

    do hi = lo, n
      if (val(hi) > dim_max) exit
    enddo

    dim_subset % dimlen = n
    dim_subset % min = dim_min
    dim_subset % max = dim_max
    dim_subset % imin = lo - 1
    dim_subset % imax = hi - 1
    dim_subset % num_values = hi - lo + 1
    dim_subset % values = val (lo:hi)

  end subroutine

  subroutine read_pressure_subset (obj, cpt, cps, errstat)
    implicit none
    type (tiof_file_type) :: obj
    type (clim_pres_type), intent(in) :: cpt
    type (clim_pres_slab_type), intent(inout) :: cps
    integer, intent(inout) :: errstat

    integer :: ilat0, ilon0, nlat, nlon
    integer :: istart(3), icount(3)

    if (errstat /= 0) return

    ilon0 = cpt % lon_subset % imin - 1
    nlon = cpt % lon_subset % num_values

    ilat0 = cpt % lat_subset % imin - 1
    nlat = cpt % lat_subset % num_values

    istart(1) = 0
    istart(2) = ilon0
    istart(3) = ilat0

    icount(1) = 1
    icount(2) = nlon
    icount(3) = nlat

    if (.not.allocated (cps % p_surf)) then
      allocate (cps % p_surf (nlat, nlon, 1))
    endif
    if (.not.allocated (cps % p_trop)) then
      allocate (cps % p_trop (nlat, nlon, 1))
    endif

    call tiof_get3d_r4 (obj, 'PS', istart, icount, cps % p_surf, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading surface pressure (PS)', errstat)
      return
    endif

    call tiof_get3d_r4 (obj, 'TROPPB', istart, icount, cps % p_trop, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading tropopause pressure (TROPPB)', errstat)
      return
    endif

    ! GEOS-CF files have pressure in Pa, but this interface returns hPa:
    cps % p_surf = cps % p_surf / 100.0
    cps % p_trop = cps % p_trop / 100.0

  end subroutine

  subroutine read_pressure_slab (cpt, cps, clim_file, errstat)
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    type (clim_pres_slab_type), intent(inout) :: cps
    character (len=*), intent(in) :: clim_file
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj

    if (errstat /= 0) return

    call tell_log (1, 'clim_module: reading: '//trim(clim_file))

    call tiof_open (clim_file, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       'read_pressure_slab: error opening file: '//trim(clim_file), &
                       errstat)
      return
    endif

    call read_pressure_subset (obj, cpt, cps, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'read_pressure_slab: failed', errstat)
    endif

    call tiof_close (obj, errstat)

  end subroutine

  subroutine maybe_alloc_subset (cpt, source_file, lon_min, lon_max, lat_min, lat_max, &
                                 nhours, errstat)
    use, intrinsic :: iso_c_binding, only : c_char
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    character (kind=c_char, len=path_bufsize), intent(in) :: source_file
    real (kind=4), intent(in) :: lon_min, lon_max, lat_min, lat_max
    integer, intent(in) :: nhours
    integer, intent(inout) :: errstat

    integer :: nlon, nlat, err
    type (tiof_file_type) :: obj

    if (errstat /= 0) return
    if (allocated (cpt % lon_subset % values)) return

    call tiof_open (source_file, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       'maybe_alloc_subset: error opening file: '//trim(source_file), &
                       errstat)
      return
    endif

    call subset_dim (obj, 'lon', 'lon', lon_min, lon_max, cpt % lon_subset, errstat)
    call subset_dim (obj, 'lat', 'lat', lat_min, lat_max, cpt % lat_subset, errstat)
    call tiof_close (obj, errstat)
    if (errstat /= 0) return

    nlon = cpt % lon_subset % num_values
    nlat = cpt % lat_subset % num_values
    allocate (cpt % p_surf(nlat,nlon,nhours), &
              cpt % p_trop(nlat,nlon,nhours), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'maybe_alloc_subset: allocate failed', errstat)
      return
    endif

  end subroutine

  !> @brief
  !> Initialize pressure climatology (struct args)
  !> @param[out] cpt   Instance of opaque @a type(clim_pres_type) to hold
  !>                   the selected pressure climatology data
  !> @param[in] year   Integer year
  !> @param[in] month  Integer month [1,12] of climatology data to read
  !> @param[in] day    Integer day [1,30] of data to read
  !> @param[in] b      Instance of @a type(clim_pres_bounds_type) providing
  !>                   the hour, longitude, latitude range of interest
  !> @param[inout] errstat        Error status code (0 on success)
  !> @param[in,optional] use_fcast   If present, and .false., the forecast
  !>                                 files are ignored.
  subroutine clim_pres_init_struct (cpt, year, month, day, b, errstat, use_fcast)
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(in) :: year, month, day
    type (clim_pres_bounds_type), intent(in) :: b
    integer, intent(inout) :: errstat
    logical, optional, intent(in) :: use_fcast

    call clim_pres_init_args (cpt, year, month, day, &
                              b % hour_beg, b % hour_end, &
                              b % lon_min, b % lon_max, &
                              b % lat_min, b % lat_max, errstat, &
                              use_fcast)
  end subroutine

  !> @brief
  !> Initialize pressure climatology (explicit args)
  !> @param[out] cpt   Instance of opaque @a type(clim_pres_type) to hold
  !>                   the selected pressure climatology data
  !> @param[in] year   Integer year
  !> @param[in] month  Integer month [1,12] of climatology data to read
  !> @param[in] day    Integer day [1,30] of data to read
  !> @param[in] hour_beg, hour_end  Selected UTC interval [hour]
  !> @param[in] lon_min, lon_max    Selected longitude range [deg]
  !> @param[in] lat_min, lat_max    Selected latitude range
  !> @param[inout] errstat        Error status code (0 on success)
  !> @param[in,optional] use_fcast   If present, and .false., the forecast
  !>                                 files are ignored.
  subroutine clim_pres_init_args (cpt, year, month, day, &
                                  hour_beg, hour_end, &
                                  lon_min, lon_max, &
                                  lat_min, lat_max, errstat, use_fcast)
    use, intrinsic :: iso_c_binding, only : c_char
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(in) :: year, month, day
    real (kind=4), intent(in) :: hour_beg, hour_end
    real (kind=4), intent(in) :: lon_min, lon_max, lat_min, lat_max
    integer, intent(inout) :: errstat
    logical, optional, intent(in) :: use_fcast

    character (kind=c_char, len=path_bufsize) :: file_month0, file_month1
    integer :: i, nhours, have_forecast_files
    integer (kind=8) :: timet
    real (kind=4) :: wt0
    type (clim_pres_slab_type) :: cps0, cps1

    if (errstat /= 0) return

    call read_pressure_eta (errstat)
    if (errstat /= 0) return

    call define_month_interp (cpt, month, day)

    nhours = ceiling(hour_end) - floor(hour_beg) + 1

    allocate (cpt % hours(nhours))
    cpt % nhours = nhours
    cpt % hours(:) = real((/(i, i=0,nhours-1)/) + floor(hour_beg), 4)
    cpt % year = year

    timet = c_make_timet (year, month, day, int(hour_beg))
    have_forecast_files = c_have_forecast_files (timet, nhours)
    Have_Forecast = (have_forecast_files == 1)

    if (present(use_fcast)) then
      if (.not.use_fcast) Have_Forecast = .false.
    endif

    if (Have_Forecast) then

      do i = 1, nhours
        timet = c_make_timet (year, month, day, int(cpt % hours(i)))
        call make_forecast_path (timet, file_month0, errstat)
        call maybe_alloc_subset (cpt, file_month0, lon_min, lon_max, lat_min, lat_max, nhours, errstat)
        if (errstat /= 0) return

        call read_pressure_slab (cpt, cps0, file_month0, errstat)
        if (errstat /= 0) return

        cpt % p_surf(:,:,i) = cps0 % p_surf (:,:,1)
        cpt % p_trop(:,:,i) = cps0 % p_trop (:,:,1)
      enddo

    else

      wt0 = cpt % month0_weight

      do i = 1, nhours
        call make_climatology_path (cpt % month0, int(cpt % hours(i)), file_month0, errstat)
        call make_climatology_path (cpt % month1, int(cpt % hours(i)), file_month1, errstat)
        call maybe_alloc_subset (cpt, file_month0, lon_min, lon_max, lat_min, lat_max, nhours, errstat)
        if (errstat /= 0) return

        call read_pressure_slab (cpt, cps0, file_month0, errstat)
        call read_pressure_slab (cpt, cps1, file_month1, errstat)
        if (errstat /= 0) return

        cpt % p_surf(:,:,i) = wt0 * cps0 % p_surf(:,:,1) + (1.0 - wt0) * cps1 % p_surf(:,:,1)
        cpt % p_trop(:,:,i) = wt0 * cps0 % p_trop(:,:,1) + (1.0 - wt0) * cps1 % p_trop(:,:,1)
      enddo

    endif

  end subroutine

  subroutine clim_query_nz (nz, errstat)
    implicit none
    integer, intent(out) :: nz
    integer, intent(inout) :: errstat

    call read_pressure_eta (errstat)
    if (errstat /= 0) return

    nz = Num_Layers

  end subroutine

  subroutine lonlat_lookup (cpt, lon, lat, ilon0, ilat0, errstat)
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: lon, lat
    integer, intent(out) :: ilon0, ilat0
    integer, intent(inout) :: errstat

    logical :: domain_error
    real (kind=4) :: dlon, lon0, lon1
    real (kind=4) :: dlat, lat0, lat1
    integer :: nlon, nlat
    character (len=msg_bufsize) :: msg

    if (errstat /= 0) return

    ilon0 = -1
    ilat0 = -1

    nlon = cpt % lon_subset % num_values
    nlat = cpt % lat_subset % num_values

    ! Assume that the lon/lat grids have fixed spacing;
    ! (lon,lat) give coordinates of grid box center.

    lon0 = cpt % lon_subset % values(1)
    lon1 = cpt % lon_subset % values(nlon)
    dlon = cpt % lon_subset % values(2) - lon0
    lon0 = lon0 - dlon/2
    lon1 = lon1 + dlon/2

    lat0 = cpt % lat_subset % values(1)
    lat1 = cpt % lat_subset % values(nlat)
    dlat = cpt % lat_subset % values(2) - lat0
    lat0 = lat0 - dlat/2
    lat1 = lat1 + dlat/2

    domain_error = ((lon < lon0) .or. (lon1 < lon) &
                    .or. (lat < lat0) .or. (lat1 < lat))

    if (domain_error) then
      write (msg, '(a,f10.1,a,f10.1)')'lonlat_lookup: domain error: lon=',lon,' lat=',lat
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    ilon0 = ceiling ((lon - lon0)/dlon)
    if (ilon0 <= 0) ilon0 = 1

    ilat0 = ceiling ((lat - lat0)/dlat)
    if (ilat0 <= 0) ilat0 = 1

  end subroutine

  !> @brief
  !> Interpolate pressure vs height
  !> @param[out] eta_a    Output Eta_A array
  !> @param[out] eta_b    Output Eta_B array
  !>
  !> \a Eta_A and \a Eta_B are the pressure parameterization arrays such that
  !> the pressure vs height is defined in terms of \a eta_a, \a eta_b, and
  !> the surface pressure as:
  !> @v+
  !> p(z) = eta_a(z) + eta_b(z) * p_surf
  !> @v-
  subroutine clim_pres_eta (eta_a, eta_b, errstat)
    implicit none
    real (kind=4), dimension(Num_Layers+1), intent(out) :: eta_a, eta_b
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    call read_pressure_eta (errstat)
    if (errstat /= 0) return

    eta_a(:) = EtaA(:)
    eta_b(:) = EtaB(:)

  end subroutine

  !> @brief
  !> Interpolate pressure vs height
  !> @param[in] cpt       Initialized instance of opaque @a type(clim_pres_type)
  !> @param[in] hour_utc  UTC hour of interest [hours]
  !> @param[in] lon, lat  Longitude, latitude coordinates of interest [deg]
  !> @param[out] pres_z   Output pressure [hPa] vs height [nz+1]
  !> @param[inout] errstat        Error status code (0 on success)
  !> @param[out] p_surf    Optional output surface pressure [hPa]
  !> @param[out] p_trop    Optional output tropopause pressure [hPa]
  subroutine clim_pres (cpt, hour_utc, lon, lat, pres_z, errstat, &
                        p_surf, p_trop)
    implicit none
    type(clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: hour_utc, lon, lat
    real (kind=4), dimension(Num_Layers+1), intent(out) :: pres_z
    integer, intent(inout) :: errstat
    real (kind=4), optional, intent(out) :: p_surf, p_trop

    integer :: ilon0, ilat0, ihr0
    real (kind=4) :: hr0, wt0, psurf, hr_min, hr_max
    character (len=msg_bufsize) :: msg

    if (errstat /= 0) return

    call lonlat_lookup (cpt, lon, lat, ilon0, ilat0, errstat)
    if (errstat /= 0) return

    hr_min = cpt % hours(1)
    hr_max = cpt % hours(cpt % nhours)

    if ((hour_utc < hr_min) .or. (hr_max < hour_utc)) then
      write (msg, '(a,f10.1)')'clim_pres: domain error: hour_utc=', hour_utc
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    hr0  = hr_min
    ihr0 = ceiling (hour_utc - hr0)
    if (ihr0 <= 0) ihr0 = 1

    if (ihr0 + 1 > cpt % nhours) then
      write (msg, '(a,f10.1)')'clim_pres: domain error: hour_utc=', hour_utc
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    ! This assumes the diurnal grid spacing is 1 hour
    wt0 = 1.0 - (hour_utc - cpt % hours(ihr0))

    psurf = (wt0 * cpt % p_surf (ilat0, ilon0, ihr0) &
             + (1.0 - wt0) * cpt % p_surf (ilat0, ilon0, ihr0+1))

    pres_z(1:Num_Layers+1) = EtaA(1:Num_Layers+1) + EtaB (1:Num_Layers+1) * psurf

    if (present(p_surf)) then
      p_surf = psurf
    endif

    if (present(p_trop)) then
      p_trop = (wt0 * cpt % p_trop (ilat0, ilon0, ihr0) &
                + (1.0 - wt0) * cpt % p_trop (ilat0, ilon0, ihr0+1))
    endif

  end subroutine

  subroutine allocate_clim_type (cpt, ct, errstat)
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    type (clim_type), intent(inout) :: ct
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    ct % nlon = cpt % lon_subset % num_values
    ct % nlat = cpt % lat_subset % num_values
    ct % nz = Num_Layers

    allocate (ct % values (ct % nz, ct % nlat, ct % nlon), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'allocate_clim_type: allocate failed', errstat)
      return
    endif

  end subroutine

  subroutine read_climatology (cpt, file, name, ct, errstat)
    use, intrinsic :: iso_c_binding, only : c_null_char
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    character (len=*), intent(in) :: file
    character (len=*), intent(in) :: name
    type (clim_type), intent(inout) :: ct
    integer, intent(inout) :: errstat

    integer :: istart(0:3), icount(0:3)
    type (tiof_file_type) :: obj

    real (kind=4), dimension (:,:,:,:), allocatable :: values
    integer :: nlon, nlat, nz

    if (errstat /= 0) return

    nlon = cpt % lon_subset % num_values
    nlat = cpt % lat_subset % num_values
    nz = ct % nz

    istart(0) = 0
    istart(1) = cpt % lon_subset % imin
    istart(2) = cpt % lat_subset % imin
    istart(3) = 0

    icount(0) = 1
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = nz

    call tell_log (1, 'clim_module: reading: '//trim(file))

    call tiof_open (file, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       'read_climatology: error opening file: '//trim(file), &
                       errstat)
      return
    endif

    allocate (values(nz,nlat,nlon,1))
    call tiof_get4d_r4 (obj, trim(name)//c_null_char, &
                        istart, icount, values, errstat)
    ct % values = reshape (values, (/nz,nlat,nlon/))

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       'read_climatology: error reading file: '//trim(file), &
                       errstat)
      return
    endif

    call tiof_close (obj, errstat)

  end subroutine

  subroutine make_climatology_path (imonth, ihour, file, errstat)
    use, intrinsic :: iso_c_binding, only : c_null_char, c_char, c_int
    implicit none
    integer (kind=c_int), intent(in) :: imonth, ihour
    character (kind=c_char, len=*), intent(inout) :: file
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    err = c_make_climatology_path (imonth, ihour, file, len(file))
    if (err /= 0) then
      call tell_error (tell_runtime_error, "make_climatology_path: failed", errstat)
      return
    endif

  end subroutine

  subroutine make_forecast_path (ttime, file, errstat)
    use, intrinsic :: iso_c_binding, only : c_null_char, c_char, c_int, c_long
    implicit none
    integer (kind=c_long), intent(in) :: ttime
    character (kind=c_char, len=*), intent(inout) :: file
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    err = c_make_forecast_path (ttime, file, len(file))
    if (err /= 0) then
      call tell_error (tell_runtime_error, "make_forecast_path: failed", errstat)
      return
    endif

  end subroutine

  subroutine interp_months (wt0, ct0, ct1, ct_avg, errstat)
    implicit none
    real (kind=4), intent(in) :: wt0
    type (clim_type), intent(in) :: ct0, ct1
    type (clim_type), intent(inout) :: ct_avg
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    ct_avg % values(:,:,:) = &
      (wt0 * ct0 % values(:,:,:) + (1.0 - wt0) * ct1 % values(:,:,:))

  end subroutine

  !> @brief
  !> Initialize tables for named variable
  !> @param[out] cst  Instance of opaque @a type(clim_val_type) to hold
  !>                  the selected values
  !> @param[in] cpt   Initialized instance of @a type(clim_pres_type)
  !> @param[in] name  String name of the value in the database files
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_val_init (cst, cpt, name, errstat)
    implicit none
    type (clim_val_type), intent(out) :: cst
    type (clim_pres_type), intent(in) :: cpt
    character (len=*), intent(in) :: name
    integer, intent(inout) :: errstat

    type (clim_type) :: ct_month0, ct_month1
    integer :: i, nhours, err
    integer (kind=8) :: timet
    real (kind=4) :: hour_beg, hour_end
    character (len=path_bufsize) :: file_month0, file_month1

    if (errstat /= 0) return

    hour_beg = cpt % hours(1)
    hour_end = cpt % hours(cpt % nhours)

    if (hour_end < hour_beg) then
      call tell_error (tell_invalid_parm_error, 'invalid hour range', errstat)
      return
    endif

    ! allocate storage for the required number of hours
    ! (usually either 2 or 3 hours)

    nhours = ceiling(hour_end) - floor(hour_beg) + 1

    allocate (cst % hours(nhours), &
              cst % clim(nhours), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'clim_val_init: allocate failed', errstat)
      return
    endif

    ! This assumes the diurnal grid spacing is 1 hour

    cst % nhours = nhours
    cst % hours(:) = real((/(i, i=0,nhours-1)/) + floor(hour_beg), 4)
    do i = 1,nhours
      call allocate_clim_type (cpt, cst % clim(i), errstat)
      if (errstat /= 0) return
    enddo

    if (Have_Forecast) then
      ! read composition forecast

      do i = 1, nhours
        timet = c_make_timet (cpt % year, cpt % month, cpt % day, int(cpt % hours(i)))
        call make_forecast_path (timet, file_month0, errstat)
        if (errstat /= 0) return
        call read_climatology (cpt, file_month0, name, cst % clim(i), errstat)
      enddo

    else
      ! read composition climatologies and perform month interpolation

      call allocate_clim_type (cpt, ct_month0, errstat)
      call allocate_clim_type (cpt, ct_month1, errstat)
      if (errstat /= 0) return

      do i = 1,nhours
        call make_climatology_path (cpt % month0, int(cst % hours(i)), file_month0, errstat)
        call make_climatology_path (cpt % month1, int(cst % hours(i)), file_month1, errstat)
        if (errstat /= 0) return

        call read_climatology (cpt, file_month0, name, ct_month0, errstat)
        call read_climatology (cpt, file_month1, name, ct_month1, errstat)
        if (errstat /= 0) return

        call interp_months (cpt % month0_weight, ct_month0, ct_month1, cst % clim(i), errstat)
        if (errstat /= 0) return
      enddo
    endif

  end subroutine

  subroutine hrlonlat_lookup (cst, cpt, hour_utc, lon, lat, &
                              ihr0, ilon0, ilat0, errstat)
    implicit none
    type (clim_val_type), intent(in) :: cst
    type (clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: hour_utc, lon, lat
    integer, intent(out) :: ihr0, ilon0, ilat0
    integer, intent(inout) :: errstat

    character (len=msg_bufsize) :: msg
    real (kind=4) :: hr0, hr_min, hr_max

    if (errstat /= 0) return

    call lonlat_lookup (cpt, lon, lat, ilon0, ilat0, errstat)
    if (errstat /= 0) return

    hr_min = cpt % hours(1)
    hr_max = cpt % hours(cpt % nhours)

    if ((hour_utc < hr_min) .or. (hr_max < hour_utc)) then
      write (msg, '(a,f10.1)')'hrlonlat_lookup: domain error: hour_utc=', hour_utc
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    ! This assumes the diurnal grid spacing is 1 hour
    hr0 = cst % hours(1)
    ihr0 = ceiling (hour_utc - hr0)
    if (ihr0 <= 0) ihr0 = 1

    if (ihr0 + 1 > cst % nhours) then
      write (msg, '(a,f10.1)')'hrlonlat_lookup: domain error: hour_utc=', hour_utc
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

  end subroutine

  !> @brief
  !> Interpolate table values vs height
  !> @param[out] cst  Initialized instance of opaque @a type(clim_val_type)
  !> @param[in] cpt   Initialized instance of @a type(clim_pres_type)
  !> @param[in] hour_utc  UTC hour of interest
  !> @param[in] lon, lat  Longitude, latitude coordinates of interest [deg]
  !> @param[out] values_z   Output value vs height.
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_val_interp (cst, cpt, hour_utc, lon, lat, values_z, &
                              errstat)
    implicit none
    type(clim_val_type), intent(in) :: cst
    type(clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: hour_utc, lon, lat
    real (kind=4), intent(out), dimension(Num_Layers) :: values_z
    integer, intent(inout) :: errstat

    integer :: ihr0, ilon0, ilat0
    real (kind=4) :: wt0

    if (errstat /= 0) return

    call hrlonlat_lookup (cst, cpt, hour_utc, lon, lat, ihr0, ilon0, ilat0, errstat)
    if (errstat /= 0) return

    ! This assumes the diurnal grid spacing is 1 hour
    wt0 = 1.0 - (hour_utc - cst % hours(ihr0))

    values_z(:) = (wt0 * cst % clim(ihr0) % values(:,ilat0,ilon0) &
                + (1.0 - wt0) * cst % clim(ihr0+1) % values(:,ilat0,ilon0))

  end subroutine

  !> @brief
  !> Compute the partial column from the volume mixing ratio and pressure
  !> @param[in] pres_z  Atmospheric pressure in @a nz layers [hPa]
  !> @param[in] vmr_z   Volume mixing ratio of a trace constituent in
  !>                    @a [nz-1] layers.
  !> @param[out] col_z  Partial column [cm^-2] in each layer [nz-1]
  !> @param[inout]  errstat  Error status code (0 on success)
  !> @param[in]  vmr_h2o  Optional volume mixing ratio of water in @a [nz-1] layers
  !>
  !> The pressure array is assumed to provide the pressure at both boundaries
  !> of each layer, therefore the pressure array dimension is one larger than
  !> the VMR and partial column arrays.  The partial column includes the
  !> variation of the gravitational acceleration with height.
  !> When the water VMR is provided, the partial column value will include
  !> a correction for the presence of water.  Otherwise, the partial column
  !> is computed for dry air.
  subroutine clim_partial_column (pres_z, vmr_z, col_z, errstat, vmr_h2o)
    implicit none
    real (kind=4), intent(in),  dimension(:)              :: pres_z
    real (kind=4), intent(in),  dimension(size(pres_z)-1) :: vmr_z
    real (kind=4), intent(out), dimension(size(pres_z)-1) :: col_z
    integer, intent(inout) :: errstat
    real (kind=4), intent(in),  dimension(size(pres_z)-1), optional :: vmr_h2o

    real (kind=4), parameter :: cm  = 1.e-2                        ! m per cm
    real (kind=4), parameter :: km  = 1.e3                         ! m per km
    real (kind=4), parameter :: hpa = 1.e2                         ! Pa per hPa
    real (kind=4), parameter :: gm_constant = 3.986004418e14       ! G * M_earth [m^3 s^-2]
    real (kind=4), parameter :: r_earth = 6371.0088e3              ! Earth mean radius [m]
    real (kind=4), parameter :: avogadros_number = 6.02214076e23   ! [mol^-1]
    real (kind=4), parameter :: mu_dry_air = 28.97e-3              ! [kg/mol]
    real (kind=4), parameter :: mu_water = 18.01528e-3             ! [kg/mol]
    real (kind=4), parameter :: pres_std_sea = 1013.25             ! sea-level pressure [hPa]

    real (kind=4), parameter :: coef = (avogadros_number / mu_dry_air) * (hpa * cm**2)
    real (kind=4), parameter :: g0 = gm_constant / r_earth**2
    real (kind=4), parameter :: f = mu_water / mu_dry_air - 1.0

    real (kind=4) :: z_km, f_earth, grav_accel, pmid
    integer :: i, nz

    if (errstat/= 0) return

    nz = size(pres_z)

    do i = 1, nz-1
      pmid = 0.5 * (pres_z(i) + pres_z(i+1))       ! layer midpoint pressure

      z_km = -16.0 * log (pmid / pres_std_sea)     ! P [hPa] -> altitude [km]
      f_earth = 1.0 + km * z_km / r_earth
      grav_accel = g0 / (f_earth * f_earth)        ! grav [m/s^2]

      col_z(i) = (coef / grav_accel) * vmr_z(i) * (pres_z(i) - pres_z(i+1))
    enddo

    if (present (vmr_h2o)) then
      col_z(:) = col_z(:) / (1.0 + f * vmr_h2o(:))
    endif

  end subroutine

  subroutine make_cloud_climatology_path (file, errstat)
    use, intrinsic :: iso_c_binding, only : c_char, c_int
    implicit none
    character (kind=c_char, len=*), target, intent(inout) :: file
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    err = c_make_cloud_climatology_path (file, len(file))
    if (err /= 0) then
      call tell_error (tell_runtime_error, "make_cloud_climatology_path: failed", errstat)
      return
    endif
  end subroutine

  !> @brief
  !> Initialize cloud climatology
  !> @param[out] cct   Instance of opaque @a type(clim_cloud_type) to hold
  !>                   the cloud climatology data
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_cloud_init (cct, errstat)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    implicit none
    type (clim_cloud_type), intent(out) :: cct
    integer, intent(inout) :: errstat

    character (kind=c_char, len=path_bufsize) :: file
    type (tiof_file_type) :: obj
    integer, dimension(3) :: istart, icount
    integer :: err

    if (errstat /= 0) return

    call make_cloud_climatology_path (file, errstat)
    if (errstat /= 0) return

    call tiof_open (file, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       'clim_cloud_init: error opening file: '//trim(file), &
                       errstat)
      return
    endif

    call tiof_inq_dimlen (obj, 'nlon', cct % nlon, errstat)
    call tiof_inq_dimlen (obj, 'nlat', cct % nlat, errstat)
    call tiof_inq_dimlen (obj, 'nmon', cct % nmon, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading pressure file dimensions', errstat)
      return
    endif

    allocate (cct % lon(cct % nlon), &
              cct % lat(cct % nlat), &
              cct % cloud_pressure (cct % nlon, cct % nlat, cct % nmon), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'clim_cloud_init: allocate failed', errstat)
      return
    endif

    istart(1) = 0
    icount(1) = cct % nlon

    call tiof_get1d_r4 (obj, 'lon'//c_null_char, &
                        istart, icount, cct % lon, errstat)

    istart(1) = 0
    icount(1) = cct % nlat

    call tiof_get1d_r4 (obj, 'lat'//c_null_char, &
                        istart, icount, cct % lat, errstat)

    istart(1) = 0
    istart(2) = 0
    istart(3) = 0
    icount(1) = cct % nmon
    icount(2) = cct % nlat
    icount(3) = cct % nlon

    call tiof_get3d_r4 (obj, 'omcldrr'//c_null_char, &
                        istart, icount, cct % cloud_pressure, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       'clim_cloud_init: error reading file: '//trim(file), &
                       errstat)
      return
    endif

    call tiof_close (obj, errstat)

  end subroutine

  !> @brief
  !> Interpolate cloud pressure
  !> @param[in] cct    Initialized instance of opaque @a type(clim_cloud_type)
  !> @param[in] month  Integer month [1,12] of climatology data to read
  !> @param[in] day    Integer day [1,30] of data to read
  !> @param[in] lon, lat  Longitude, latitude coordinates of interest [deg]
  !> @param[out] cloud_pressure   Output cloud pressure [hPa]
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_cloud (cct, month, day, lon, lat, cloud_pressure, errstat)
    implicit none
    type (clim_cloud_type), intent(in) :: cct
    integer, intent(in) :: month, day
    real (kind=4), intent(in) :: lon, lat
    real (kind=4), intent(out) :: cloud_pressure
    integer, intent(inout) :: errstat

    integer :: ilon0, ilat0, month0, month1
    real (kind=4) :: lon0, lat0, dlon, dlat, f, wt0

    if (errstat /= 0) return

    f = (day - 1)/ mean_days_per_month

    month0 = month

    if (f < 0.5) then
      wt0 = 0.5 + f
      month1 = month - 1
      if (month1 < 1) month1 = 12
    else
      wt0 = 1.5 - f
      month1 = month + 1
      if (month1 > 12) month1 = 1
    endif

    lon0 = cct % lon(1)
    dlon = cct % lon(2) - lon0
    ilon0 = ceiling ((lon - lon0)/dlon)
    if (ilon0 == 0) ilon0 = 1

    lat0 = cct % lat(1)
    dlat = cct % lat(2) - lat0
    ilat0 = ceiling ((lat - lat0)/dlat)
    if (ilat0 == 0) ilat0 = 1

    cloud_pressure = (wt0 * cct % cloud_pressure (ilon0, ilat0, month0) &
                      + (1.0 - wt0) * cct % cloud_pressure (ilon0, ilat0, month1))

  end subroutine

end module clim_module
