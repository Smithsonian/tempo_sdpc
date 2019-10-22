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
!! The z coordinate is defined in terms of pressure.  The pressure
!! climatology is contained in a set of 12 files, one for each month.
!! Each (monthly) pressure climatology file provides a parameterization
!! of the pressure, P(z,lat,lon,hour).
!!
module clim_module
  use netcdf, only : nf90_nowrite
  use tell_module
  use tio_module
  implicit none

  private

  public clim_pres, clim_pres_init, clim_pres_nz
  public clim_species_vmr, clim_species_init

  interface clim_pres_init
    module procedure clim_pres_init_args
    module procedure clim_pres_init_struct
  end interface

  integer, private, parameter :: path_bufsize = 1024, msg_bufsize = 128

  type :: dim_subset_type
    integer :: dimlen
    real (kind=4) :: min, max    ! domain limits [min, max)
    integer ::      imin, imax   ! file variable array index limits implied by [min, max)
    integer :: num_values        ! = imax-imin+1

    real (kind=4), dimension(:), allocatable :: values  ! (num_values)
  end type

  type, public :: clim_pres_type
    private
    real (kind=4) :: fmonth  ! label for month-interpolated result
    type(dim_subset_type) :: lon_subset  ! nlon, grid box centers
    type(dim_subset_type) :: lat_subset  ! nlat, grid box centers
    type(dim_subset_type) :: hour_subset ! nhour, interval start times
    integer :: nz
    real (kind=4), dimension(:), allocatable :: eta_a, eta_b  ! (nz)
    real (kind=4), dimension(:,:,:), allocatable :: p_surf    ! (nlat,nlon,nhour)
    real (kind=4), dimension(:,:,:), allocatable :: p_trop    ! (nlat,nlon,nhour)
  end type

  type, public :: clim_pres_bounds_type
    real (kind=4) :: hour_beg, hour_end
    real (kind=4) :: lon_min, lon_max
    real (kind=4) :: lat_min, lat_max
  end type

  ! Note: clim_type % nz = clim_pres_type % nz - 1
  type :: clim_type
    integer :: nz, nlat, nlon
    real (kind=4), dimension(:,:,:), allocatable :: vmr ! (nz, nlat, nlon)
  end type

  type, public :: clim_species_type
    private
    integer :: nhours
    real (kind=4), dimension (:), allocatable :: hours   ! (nhours)
    type(clim_type), dimension(:), allocatable :: clim   ! (nhours)
  end type

  interface
    function c_make_climatology_filename (species, month, hour, path, path_buflen) &
        bind (c, name='make_climatology_filename')
      use, intrinsic :: iso_c_binding, only : c_char, c_int
      implicit none
      character (kind=c_char), intent(in) :: species(*)
      integer (c_int), value, intent(in) :: month, hour, path_buflen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_make_climatology_filename
    end function
  end interface

  interface
    function c_make_pressure_filename (month, path, path_buflen) &
        bind (c, name='make_pressure_filename')
      use, intrinsic :: iso_c_binding, only : c_char, c_int
      implicit none
      integer (c_int), value, intent(in) :: month, path_buflen
      character (kind=c_char), intent(out) :: path(*)
      integer (c_int) :: c_make_pressure_filename
    end function
  end interface

contains

  subroutine alloc_pres_type (cpt, nhour, nlon, nlat, nz, errstat)
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(in) :: nhour, nlon, nlat, nz
    integer, intent(inout) :: errstat
    integer :: err

    if (errstat /= 0) return

    allocate (cpt % eta_a(nz), &
              cpt % eta_b(nz), &
              cpt % p_surf(nlat,nlon,nhour), &
              cpt % p_trop(nlat,nlon,nhour), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'alloc_pres_type: allocate failed', errstat)
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

  subroutine read_pressure_subset (obj, cpt, errstat)
    implicit none
    type (tiof_file_type) :: obj
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(inout) :: errstat

    integer :: ilat0, ilon0, ihour0, nlat, nlon, nhour
    integer :: istart(3), icount(3)

    if (errstat /= 0) return

    call tiof_inq_dimlen (obj, 'z', cpt % nz, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading pressure file dimensions', errstat)
      return
    endif

    ilon0 = cpt % lon_subset % imin - 1
    nlon = cpt % lon_subset % num_values

    ilat0 = cpt % lat_subset % imin - 1
    nlat = cpt % lat_subset % num_values

    ihour0 = cpt % hour_subset % imin - 1
    nhour = cpt % hour_subset % num_values

    call alloc_pres_type (cpt, nhour, nlon, nlat, cpt % nz, errstat)
    if (errstat /= 0) return

    call tiof_get1d_r4 (obj, 'EtaA', (/0/), (/cpt % nz/), cpt % eta_a, errstat)
    call tiof_get1d_r4 (obj, 'EtaB', (/0/), (/cpt % nz/), cpt % eta_b, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading pressure file', errstat)
      return
    endif

    istart(1) = ihour0
    istart(2) = ilon0
    istart(3) = ilat0

    icount(1) = nhour
    icount(2) = nlon
    icount(3) = nlat

    call tiof_get3d_r4 (obj, 'SurfacePressure', istart, icount, &
                        cpt % p_surf, errstat)
    call tiof_get3d_r4 (obj, 'TropopausePressure', istart, icount, &
                        cpt % p_trop, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, 'reading pressure file', errstat)
      return
    endif

  end subroutine

  subroutine read_pressure_file (cpt, pressure_file, hour_beg, hour_end, &
                                 lon_min, lon_max, lat_min, lat_max, errstat)
    implicit none
    type (clim_pres_type), intent(out) :: cpt
    character (len=*), intent(in) :: pressure_file
    real (kind=4), intent(in) :: hour_beg, hour_end
    real (kind=4), intent(in) :: lon_min, lon_max, lat_min, lat_max
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj

    if (errstat /= 0) return

    call tiof_open (pressure_file, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       'read_pressure_file: error opening file: '//trim(pressure_file), &
                       errstat)
      return
    endif

    call subset_dim (obj, 'Longitude', 'x', lon_min, lon_max, cpt % lon_subset, errstat)
    call subset_dim (obj, 'Latitude', 'y', lat_min, lat_max, cpt % lat_subset, errstat)
    call subset_dim (obj, 'Hour', 't', hour_beg, hour_end, cpt % hour_subset, errstat)

    call read_pressure_subset (obj, cpt, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'read_pressure_file: failed', errstat)
    endif

    call tiof_close (obj, errstat)

  end subroutine

  subroutine make_pressure_filename (imonth, file, errstat)
    use, intrinsic :: iso_c_binding, only : c_char, c_int
    implicit none
    integer (kind=c_int), intent(in) :: imonth
    character (kind=c_char, len=*), target, intent(inout) :: file
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    err = c_make_pressure_filename (imonth, file, len(file))
    if (err /= 0) then
      call tell_error (tell_runtime_error, "make_pressure_filename: failed", errstat)
      return
    endif
  end subroutine

  !> @brief
  !> Initialize pressure climatology (struct args)
  !> @param[out] cpt   Instance of opaque @a type(clim_pres_type) to hold
  !>                   the selected pressure climatology data
  !> @param[in] month  Integer month [1,12] of climatology data to read
  !> @param[in] day    Integer day [1,30] of data to read
  !> @param[in] b      Instance of @a type(clim_pres_bounds_type) providing
  !>                   the hour, longitude, latitude range of interest
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_pres_init_struct (cpt, month, day, b, errstat)
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(in) :: month, day
    type (clim_pres_bounds_type), intent(in) :: b
    integer, intent(inout) :: errstat

    call clim_pres_init_args (cpt, month, day, &
                              b % hour_beg, b % hour_end, &
                              b % lon_min, b % lon_max, &
                              b % lat_min, b % lat_max, errstat)
  end subroutine

  !> @brief
  !> Initialize pressure climatology (explicit args)
  !> @param[out] cpt   Instance of opaque @a type(clim_pres_type) to hold
  !>                   the selected pressure climatology data
  !> @param[in] month  Integer month [1,12] of climatology data to read
  !> @param[in] day    Integer day [1,30] of data to read
  !> @param[in] hour_beg, hour_end  Selected UTC interval [hour]
  !> @param[in] lon_min, lon_max    Selected longitude range [deg]
  !> @param[in] lat_min, lat_max    Selected latitude range
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_pres_init_args (cpt, month, day, &
                                  hour_beg, hour_end, &
                                  lon_min, lon_max, &
                                  lat_min, lat_max, errstat)
    use, intrinsic :: iso_c_binding, only : c_char
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer, intent(in) :: month, day
    real (kind=4), intent(in) :: hour_beg, hour_end
    real (kind=4), intent(in) :: lon_min, lon_max, lat_min, lat_max
    integer, intent(inout) :: errstat

    real (kind=4), parameter :: mean_days_per_month = 365.25/12
    character (kind=c_char, len=path_bufsize) :: file_month0, file_month1
    integer :: month0, month1
    real (kind=4) :: fmonth, wt0
    type (clim_pres_type) :: cpt1

    if (errstat /= 0) return

    fmonth = month + (day - 1)/mean_days_per_month

    month0 = month
    month1 = modulo (month0 + 1, 12)
    wt0 = 1.0 - (fmonth - month0)

    call make_pressure_filename (month0, file_month0, errstat)
    call read_pressure_file (cpt, file_month0, hour_beg, hour_end, &
                             lon_min, lon_max, lat_min, lat_max, errstat)
    if (errstat /= 0) return

    call make_pressure_filename (month1, file_month1, errstat)
    call read_pressure_file (cpt1, file_month1, hour_beg, hour_end, &
                             lon_min, lon_max, lat_min, lat_max, errstat)
    if (errstat /= 0) return

    cpt % fmonth = fmonth
    cpt % eta_a(:)      = (wt0 * cpt % eta_a(:)      + (1.0 - wt0) * cpt1 % eta_a(:))
    cpt % eta_b(:)      = (wt0 * cpt % eta_b(:)      + (1.0 - wt0) * cpt1 % eta_b(:))
    cpt % p_surf(:,:,:) = (wt0 * cpt % p_surf(:,:,:) + (1.0 - wt0) * cpt1 % p_surf(:,:,:))
    cpt % p_trop(:,:,:) = (wt0 * cpt % p_trop(:,:,:) + (1.0 - wt0) * cpt1 % p_trop(:,:,:))

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

  function clim_pres_nz (cpt)
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    integer :: clim_pres_nz
    clim_pres_nz = cpt % nz
  end function

  !> @brief
  !> Interpolate pressure vs height
  !> @param[in] cpt       Initialized instance of opaque @a type(clim_pres_type)
  !> @param[in] hour_utc  UTC hour of interest
  !> @param[in] lon, lat  Longitude, latitude coordinates of interest
  !> @param[out] pres_z   Output pressure [hPa] vs height
  !> @param[inout] errstat        Error status code (0 on success)
  !> @param[out] eta_a, eta_b    Optional output pressure parameters
  !>                             p(z) = eta_a(z) + eta_b(z) * p_surf
  !> @param[out] p_surf    Optional output surface pressure
  !> @param[out] p_trop    Optional output tropopause pressure
  subroutine clim_pres (cpt, hour_utc, lon, lat, pres_z, errstat, &
                        eta_a, eta_b, p_surf, p_trop)
    implicit none
    type(clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: hour_utc, lon, lat
    real (kind=4), dimension(cpt % nz), intent(out) :: pres_z
    integer, intent(inout) :: errstat
    real (kind=4), dimension(cpt % nz), optional, intent(out) :: eta_a, eta_b
    real (kind=4), optional, intent(out) :: p_surf, p_trop

    integer :: ilon0, ilat0, ihr0
    real (kind=4) :: hr0, wt0, psurf
    character (len=msg_bufsize) :: msg

    if (errstat /= 0) return

    call lonlat_lookup (cpt, lon, lat, ilon0, ilat0, errstat)
    if (errstat /= 0) return

    if ((hour_utc < cpt % hour_subset % min) &
        .or. (cpt % hour_subset % max < hour_utc)) then
      write (msg, '(a,f10.1)')'clim_pres: domain error: hour_utc=', hour_utc
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    hr0  = cpt % hour_subset % values(1)
    ihr0 = ceiling (hour_utc - hr0)
    if (ihr0 <= 0) ihr0 = 1

    if (ihr0 + 1 > cpt % hour_subset % num_values) then
      write (msg, '(a,f10.1)')'clim_pres: domain error: hour_utc=', hour_utc
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    ! This assumes the diurnal grid spacing is 1 hour
    wt0 = 1.0 - (hour_utc - cpt % hour_subset % values(ihr0))

    psurf = (wt0 * cpt % p_surf (ilat0, ilon0, ihr0) &
             + (1.0 - wt0) * cpt % p_surf (ilat0, ilon0, ihr0+1))

    pres_z(:) = cpt % eta_a(:) + cpt % eta_b (:) * psurf

    if (present(eta_a)) then
      eta_a(:) = cpt % eta_a(:)
    endif

    if (present(eta_b)) then
      eta_b(:) = cpt % eta_b(:)
    endif

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
    ct % nz = cpt % nz - 1

    allocate (ct % vmr (ct % nz, ct % nlat, ct % nlon), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, 'allocate_clim_type: allocate failed', errstat)
      return
    endif

  end subroutine

  subroutine read_climatology (cpt, file, species, ct, errstat)
    use, intrinsic :: iso_c_binding, only : c_null_char
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    character (len=*), intent(in) :: file
    character (len=*), intent(in) :: species
    type (clim_type), intent(inout) :: ct
    integer, intent(inout) :: errstat

    integer :: istart(3), icount(3)
    type (tiof_file_type) :: obj

    if (errstat /= 0) return

    istart(1) = cpt % lon_subset % imin
    istart(2) = cpt % lat_subset % imin
    istart(3) = 0

    icount(1) = cpt % lon_subset % num_values
    icount(2) = cpt % lat_subset % num_values
    icount(3) = ct % nz

    call tiof_open (file, obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       'read_climatology: error opening file: '//trim(file), &
                       errstat)
      return
    endif

    call tiof_get3d_r4 (obj, 'TRC'//trim(species)//c_null_char, &
                        istart, icount, ct % vmr, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       'read_climatology: error reading file: '//trim(file), &
                       errstat)
      return
    endif

    call tiof_close (obj, errstat)

  end subroutine

  subroutine make_climatology_filename (species, imonth, ihour, file, errstat)
    use, intrinsic :: iso_c_binding, only : c_null_char, c_char, c_int
    implicit none
    character (kind=c_char, len=*), intent(in) :: species
    integer (kind=c_int), intent(in) :: imonth, ihour
    character (kind=c_char, len=*), intent(inout) :: file
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat /= 0) return

    err = c_make_climatology_filename (trim(species)//c_null_char, imonth, ihour, file, len(file))
    if (err /= 0) then
      call tell_error (tell_runtime_error, "make_climatology_filename: failed", errstat)
      return
    endif

  end subroutine

  subroutine interp_months (cpt, wt0, ct0, ct1, ct_avg, errstat)
    implicit none
    type (clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: wt0
    type (clim_type), intent(in) :: ct0, ct1
    type (clim_type), intent(inout) :: ct_avg
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    ct_avg % vmr(:,:,:) = &
      (wt0 * ct0 % vmr(:,:,:) + (1.0 - wt0) * ct1 % vmr(:,:,:))
  end subroutine

  !> @brief
  !> Initialize trace gas species climatology
  !> @param[out] cst  Instance of opaque @a type(clim_species_type) to hold
  !>                  the selected trace gas species climatology data
  !> @param[in] cpt   Initialized instance of @a type(clim_pres_type)
  !> @param[in] species  String name of the trace gas species (must exactly
  !>                     match a species name in the climatology database)
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_species_init (cst, cpt, species, errstat)
    implicit none
    type (clim_species_type), intent(out) :: cst
    type (clim_pres_type), intent(in) :: cpt
    character (len=*), intent(in) :: species
    integer, intent(inout) :: errstat

    type (clim_type) :: ct_month0, ct_month1
    integer :: i, nhours, month0, month1, err
    real (kind=4) :: hour_beg, hour_end
    real (kind=4) :: weight0
    character (len=path_bufsize) :: file_month0, file_month1

    if (errstat /= 0) return

    hour_beg = cpt % hour_subset % min
    hour_end = cpt % hour_subset % max

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
      call tell_error (tell_malloc_error, 'alloc_pres_type: allocate failed', errstat)
      return
    endif

    ! This assumes the diurnal grid spacing is 1 hour

    cst % nhours = nhours
    cst % hours(:) = real((/(i, i=0,nhours-1)/) + floor(hour_beg), 4)
    do i = 1,nhours
      call allocate_clim_type (cpt, cst % clim(i), errstat)
      if (errstat /= 0) return
    enddo

    ! read climatologies and perform month interpolation

    call allocate_clim_type (cpt, ct_month0, errstat)
    call allocate_clim_type (cpt, ct_month1, errstat)
    if (errstat /= 0) return

    month0 = floor(cpt % fmonth)
    month1 = modulo (month0 + 1, 12)
    weight0 = 1.0 - (cpt % fmonth - month0)

    do i = 1,nhours
      call make_climatology_filename (species, month0, int(cst % hours(i)), file_month0, errstat)
      call read_climatology (cpt, file_month0, species, ct_month0, errstat)
      if (errstat /= 0) return

      call make_climatology_filename (species, month1, int(cst % hours(i)), file_month1, errstat)
      call read_climatology (cpt, file_month1, species, ct_month1, errstat)
      if (errstat /= 0) return

      call interp_months (cpt, weight0, ct_month0, ct_month1, cst % clim(i), errstat)
      if (errstat /= 0) return
    enddo

  end subroutine

  subroutine hrlonlat_lookup (cst, cpt, hour_utc, lon, lat, &
                              ihr0, ilon0, ilat0, errstat)
    implicit none
    type (clim_species_type), intent(in) :: cst
    type (clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: hour_utc, lon, lat
    integer, intent(out) :: ihr0, ilon0, ilat0
    integer, intent(inout) :: errstat

    character (len=msg_bufsize) :: msg
    real (kind=4) :: hr0

    if (errstat /= 0) return

    call lonlat_lookup (cpt, lon, lat, ilon0, ilat0, errstat)
    if (errstat /= 0) return

    if ((hour_utc < cpt % hour_subset % min) &
        .or. (cpt % hour_subset % max < hour_utc)) then
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
  !> Interpolate species volume mixing ratio vs height
  !> @param[out] cst  Initialized instance of opaque @a type(clim_species_type)
  !> @param[in] cpt   Initialized instance of @a type(clim_pres_type)
  !> @param[in] hour_utc  UTC hour of interest
  !> @param[in] lon, lat  Longitude, latitude coordinates of interest
  !> @param[out] vmr_z   Output volume mixing ratio vs height.
  !>                     For (nz) pressure values, the VMR array size is (nz-1).
  !> @param[inout] errstat        Error status code (0 on success)
  subroutine clim_species_vmr (cst, cpt, hour_utc, lon, lat, vmr_z, errstat)
    implicit none
    type(clim_species_type), intent(in) :: cst
    type(clim_pres_type), intent(in) :: cpt
    real (kind=4), intent(in) :: hour_utc, lon, lat
    real (kind=4), intent(out), dimension(cpt % nz - 1) :: vmr_z
    integer, intent(inout) :: errstat

    integer :: ihr0, ilon0, ilat0
    real (kind=4) :: wt0

    if (errstat /= 0) return

    call hrlonlat_lookup (cst, cpt, hour_utc, lon, lat, ihr0, ilon0, ilat0, errstat)
    if (errstat /= 0) return

    ! This assumes the diurnal grid spacing is 1 hour
    wt0 = 1.0 - (hour_utc - cst % hours(ihr0))

    vmr_z(:) = (wt0 * cst % clim(ihr0) % vmr(:,ilat0,ilon0) &
                + (1.0 - wt0) * cst % clim(ihr0+1) % vmr(:,ilat0,ilon0))

  end subroutine

end module clim_module
