!> Read meteorology parameters from real or simulated NAM meteorology files
module met_module
  use, intrinsic :: iso_c_binding, only: c_float, c_int, c_ptr
  use tell_module
  use eccodes
  implicit none

  public met_list_interp_f, read_synth_met_data

  !> Fortran interface for C struct \a Met_Value_Type
  type, bind(c), public :: met_value_type
    real (c_float) :: pressure_surface
    real (c_float) :: pressure_tropopause
    type (c_ptr) :: temperature_on_isobar
    type (c_ptr) :: isobars
    integer (c_int) :: num_isobars
  end type

  interface
    function met_list_new (flags) bind (c, name='met_list_new')
      use, intrinsic :: iso_c_binding, only : c_ptr, c_int
      implicit none
      integer (c_int), value :: flags
      type (c_ptr) :: met_list_new
    end function
  end interface

  interface
    subroutine met_list_free (met_list) bind (c, name='met_list_free')
      use, intrinsic :: iso_c_binding, only : c_ptr
      implicit none
      type (c_ptr), value :: met_list
    end subroutine
  end interface

  interface
    function met_list_add_file (met_list, path) bind (c, name='met_list_add_file')
      use, intrinsic :: iso_c_binding, only : c_ptr, c_char, c_int
      implicit none
      type (c_ptr), value :: met_list
      character (kind=c_char,len=1) :: path
      integer (c_int) :: met_list_add_file
    end function
  end interface

  interface
    function met_list_interp (met_list, lon, lat, mvt) bind (c, name='met_list_interp')
      use, intrinsic :: iso_c_binding, only : c_ptr, c_float, c_int
      implicit none
      type (c_ptr), value :: met_list
      real (c_float), value :: lon, lat
      type (c_ptr), value :: mvt
      integer (c_int) :: met_list_interp
    end function
  end interface

contains

  subroutine met_list_interp_f (met_list, lon, lat, errstat, &
                                psurf, ptrop, isobars, temp_on_isobar)
    use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_loc
    implicit none
    type (c_ptr), value :: met_list
    real (kind=4), intent(in) :: lon, lat
    integer, intent(inout) :: errstat
    real (kind=4), optional, intent(out) :: psurf
    real (kind=4), optional, intent(out) :: ptrop
    real (kind=4), optional, intent(in), dimension(:), target :: isobars
    real (kind=4), optional, intent(inout), dimension(:), target :: temp_on_isobar

    type (met_value_type), target :: mvt
    integer :: status

    if (errstat /= 0) return

    if (present(isobars) .and. present(temp_on_isobar)) then
      if (size(isobars) .ne. size(temp_on_isobar)) then
        call tell_error (tell_runtime_error, &
                         'size(isobars) /= size(temp_on_isobar)', errstat)
        return
      endif
      mvt % num_isobars = size(isobars)
      mvt % isobars = c_loc (isobars)
      mvt % temperature_on_isobar = c_loc(temp_on_isobar)
    else
      mvt % num_isobars = 0
      mvt % isobars = c_null_ptr
      mvt % temperature_on_isobar = c_null_ptr
    endif

    status = met_list_interp (met_list, lon, lat, c_loc(mvt))
    if (status /= 0) then
      call tell_error (tell_runtime_error, &
                       'met_list_interp_f: interpolation failed', &
                       errstat)
      return
    endif

    if (present(psurf)) psurf = mvt % pressure_surface
    if (present(ptrop)) ptrop = mvt % pressure_tropopause

  end subroutine

  !> read parameters from synthetic meteorology file
  !--------------------------------------------------------------------------
  !
  !> @param[in]  metfile  filename of netcdf synthetic meteorology file
  !> @param[in]  lat      latitude of target pixel
  !> @param[in]  lon      longitude of target pixel
  !> @param[out] troppres tropopause pressure value (hPa)
  !> @param[out] surfpres surface pressure value (hPa) [OPTIONAL]
  !> @param[out] tprof    teperature profile vector [OPTIONAL]
  !> @param      errstat  error handling integer, non-zero indicates failure
  !
  !> @author     E. O'Sullivan Jan 2018
  !--------------------------------------------------------------------------
  subroutine read_synth_met_data (metfile, lat, lon, troppres, surfpres, &
       tprof, errstat)
    use netcdf, only : nf90_nowrite
    use tio_module
    implicit none

    ! number of pressure levels (fixed)
    integer(kind=4) , parameter :: nlev=72

    !Input variables
    character (len=*), intent(in) :: metfile
    real (kind=4), intent(in) :: lon, lat

    !Output variables
    real(kind=4), intent(out) :: troppres
    real(kind=4), intent(out), optional :: surfpres
    real(kind=4), intent(out), dimension(nlev), optional :: tprof
    integer(kind=4), intent(inout) :: errstat

    !local variables
    integer(kind=4) :: nlon, nlat
    integer(kind=4), dimension(2) :: lonlatidx
    real(kind=4), dimension(:,:), allocatable :: longrid, latgrid
    real(kind=4), dimension(1,1) :: tmp_surfpres, tmp_troppres
    real(kind=4), dimension(1,1,nlev) :: tmp_tprof

    type (tiof_file_type) :: metobj

    if (errstat /= 0) return

    !open the file
    call tiof_open (metfile, metobj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_synth_met_data: failed to open synth meteorology file", &
           errstat)
      return
    endif

    !determine dimensions of meteorology grid
    call tiof_inq_dimlen (metobj, "x", nlon, errstat)
    call tiof_inq_dimlen (metobj, "y", nlat, errstat)

    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_troppres: failed to read dimensions", errstat)
      return
    endif

    !allocate arrays where necessary
    allocate(longrid(0:nlon-1,0:nlat-1), latgrid(0:nlon-1,0:nlat-1), &
         stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "read_synth_troppres: failed to allocate arrays", errstat)
      return
    endif

    !read in lat, lon values
    call tiof_get2d_r4(metobj, "lon", [0,0], [nlat,nlon], &
         longrid(0:nlon-1,0:nlat-1), errstat)
    call tiof_get2d_r4(metobj, "lat", [0,0], [nlat,nlon], &
         latgrid(0:nlon-1,0:nlat-1), errstat)
    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_troppres: failed to read lon, lat positions", errstat)
      return
    endif

    !determine the nearest location to the target lon, lat
    lonlatidx = minloc(abs(latgrid-lat)+abs(longrid-lon))
    !correct each index by -1 since we are using 0-indexed arrays
    lonlatidx = lonlatidx-1

    !read in values to be output, convert pressures from Pa to hPa
    call tiof_get2d_r4(metobj, "TROPPB", [lonlatidx(2),lonlatidx(1)], [1,1], &
         tmp_troppres, errstat)
    troppres = real(tmp_troppres(1,1)/100.0d0, kind=4)
    if (present(surfpres)) then
      call tiof_get2d_r4(metobj, "PS", [lonlatidx(2),lonlatidx(1)], [1,1], &
           tmp_surfpres, errstat)
      surfpres = real(tmp_surfpres(1,1)/100.0d0, kind=4)
    endif
    if (present(tprof)) then
      call tiof_get3d_r4(metobj, "T", [0,lonlatidx(2),lonlatidx(1)], &
           [nlev,1,1], tmp_tprof, errstat)
      tprof = tmp_tprof(1,1,:)
    endif
    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_troppres: failed to read pressure, temp", errstat)
      return
    endif

    !close the file
    call tiof_close (metobj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_synth_troppres: failed to close synth meteorology file", &
           errstat)
      return
    endif

    !deallocate arrays
    if (allocated(longrid)) deallocate(longrid, stat=errstat)
    if (allocated(latgrid)) deallocate(latgrid, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "read_synth_troppres: failed to deallocate arrays", errstat)
      return
    endif

  end subroutine read_synth_met_data

end module met_module
