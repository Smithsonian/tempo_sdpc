!> Read meteorology parameters from real or simulated NAM meteorology files
module met_module
  use, intrinsic :: iso_c_binding, only: c_float, c_int, c_ptr
  use tell_module
  use tio_module
  use eccodes
  implicit none

  public met_list_interp_f
  public open_synth_met_data, close_synth_met_data, read_synth_met_data

  type, public :: synth_met_type
    private
    real (kind=4), dimension(:,:), allocatable :: longrid, latgrid
    integer (kind=4) :: nlon, nlat
    type (tiof_file_type) :: metobj
  end type

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

  interface
    function met_linear_interp (x0, y0, n0, n, x, y) bind (c, name='met_linear_interp')
      use, intrinsic :: iso_c_binding, only : c_double, c_int
      real (kind=c_double), dimension(*), intent(in) :: x0, y0, x
      real (kind=c_double), dimension(*), intent(out) :: y
      integer (c_int), value :: n0, n
      integer (c_int) :: met_linear_interp
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

  subroutine open_synth_met_data (smt, metfile, errstat)
    use netcdf, only : nf90_nowrite
    character (len=*), intent(in) :: metfile
    type (synth_met_type), intent(inout) :: smt
    integer, intent(inout) :: errstat

    integer :: nlat, nlon

    if (errstat /= 0) return

    !open the file
    call tiof_open (metfile, smt % metobj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_synth_met_data: failed to open synth meteorology file", &
           errstat)
      return
    endif

    !determine dimensions of meteorology grid
    call tiof_inq_dimlen (smt % metobj, "x", nlon, errstat)
    call tiof_inq_dimlen (smt % metobj, "y", nlat, errstat)

    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_met_data: failed to read dimensions", errstat)
      return
    endif

    smt % nlon = nlon
    smt % nlat = nlat

    !allocate arrays where necessary
    allocate(smt % longrid(0:nlon-1,0:nlat-1), &
             smt % latgrid(0:nlon-1,0:nlat-1), &
             stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "read_synth_met_data: failed to allocate arrays", errstat)
      return
    endif

    !read in lat, lon values
    call tiof_get2d_r4(smt % metobj, "lon", [0,0], [nlat,nlon], &
                       smt % longrid(0:nlon-1,0:nlat-1), errstat)
    call tiof_get2d_r4(smt % metobj, "lat", [0,0], [nlat,nlon], &
                       smt % latgrid(0:nlon-1,0:nlat-1), errstat)
    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "read_synth_met_data: failed to read lon, lat positions", errstat)
      return
    endif

  end subroutine

  subroutine close_synth_met_data (smt, errstat)
    type (synth_met_type), intent(inout) :: smt
    integer, intent(inout) :: errstat

    !close the file
    call tiof_close (smt % metobj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_synth_met_data: failed to close synth meteorology file", &
           errstat)
      return
    endif

    !deallocate arrays
    if (allocated(smt % longrid)) deallocate(smt % longrid, stat=errstat)
    if (allocated(smt % latgrid)) deallocate(smt % latgrid, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "read_synth_met_data: failed to deallocate arrays", errstat)
      return
    endif

  end subroutine

  !> read parameters from synthetic meteorology file
  !--------------------------------------------------------------------------
  !
  !> @param[in]  smt      Instance of synth_met_type from open_synth_met_data()
  !> @param[in]  lat      latitude of target pixel
  !> @param[in]  lon      longitude of target pixel
  !> @param[out] troppres tropopause pressure value (hPa)
  !> @param      errstat  error handling integer, non-zero indicates failure
  !> @param[out] surfpres surface pressure value (hPa) [OPTIONAL]
  !> @param[in]  pprof    pressure profile vector (hPa) [OPTIONAL]
  !> @param[out] tprof    temperature profile vector [OPTIONAL]
  !
  ! The optional \a pprof, and \a tprof vectors are ordered so that
  ! array index 1 is at the surface.
  !
  !> @author     E. O'Sullivan Jan 2018
  !--------------------------------------------------------------------------
  subroutine read_synth_met_data (smt, lat, lon, troppres, errstat, &
                                  surfpres, tprof, pprof)
    use tio_module
    implicit none

    type (synth_met_type), intent(inout) :: smt

    ! number of pressure levels (fixed)
    integer(kind=c_int) , parameter :: nlev=72

    !Input variables
    real (kind=4), intent(in) :: lon, lat

    !Output variables
    real(kind=4), intent(out) :: troppres
    real(kind=4), intent(out), optional :: surfpres
    real(kind=8), intent(out), dimension(:), optional, target :: tprof
    real(kind=8), intent(out), dimension(:), optional :: pprof
    integer(kind=4), intent(inout) :: errstat

    !local variables
    integer (kind=c_int) :: nprof
    integer(kind=4) :: err
    integer(kind=4), dimension(2) :: lonlatidx
    real(kind=4), dimension(1,1) :: tmp_surfpres, tmp_troppres
    real(kind=4), dimension(1,1,nlev) :: tmp_tprof, tmp_pprof
    real(kind=8), dimension(nlev) :: tlev_d, plev_d

    if (errstat /= 0) return

    !determine the nearest location to the target lon, lat
    lonlatidx = minloc(abs(smt % latgrid-lat)+abs(smt % longrid-lon))
    !correct each index by -1 since we are using 0-indexed arrays
    lonlatidx = lonlatidx-1

    !read in values to be output, convert pressures from Pa to hPa
    call tiof_get2d_r4(smt % metobj, "TROPPB", [lonlatidx(2),lonlatidx(1)], [1,1], &
         tmp_troppres, errstat)
    troppres = real(tmp_troppres(1,1)/100.0d0, kind=4)
    if (present(surfpres)) then
      call tiof_get2d_r4(smt % metobj, "PS", [lonlatidx(2),lonlatidx(1)], [1,1], &
           tmp_surfpres, errstat)
      surfpres = real(tmp_surfpres(1,1)/100.0d0, kind=4)
    endif

    if (present(tprof)) then
      ! read the temperature profile
      call tiof_get3d_r4(smt % metobj, "T", [0,lonlatidx(2),lonlatidx(1)], &
                         [nlev,1,1], tmp_tprof, errstat)
      if (errstat /= 0) then
        call tell_error(tell_io_error, &
                        "read_synth_met_data: failed to read temp profile", errstat)
        return
      endif

      if (present(pprof)) then
        ! If an array of isobars was provided, we'll interpolate the
        ! temperatures onto that pressure grid.
        ! Read the pressure grid, and interpolate on arrays of doubles
        ! because the interpolation routine wants doubles.
        call tiof_get3d_r4(smt % metobj, "PL", [0,lonlatidx(2),lonlatidx(1)], &
                           [nlev,1,1], tmp_pprof, errstat)
        if (errstat /= 0) then
          call tell_error(tell_io_error, &
                          "read_synth_met_data: failed to read pressure profile", errstat)
          return
        endif
        tlev_d(:) = real(tmp_tprof(1,1,1:nlev), kind=8)
        plev_d(:) = real(tmp_pprof(1,1,1:nlev)/100.0, kind=8) ! convert to hPa
        nprof = size(pprof)
        ! The file stores pressures in increasing order, as wanted by the
        ! interpolation routine, but the calling routine will provide the
        ! pressures in decreasing order, hence the need to reverse array order:
        err = met_linear_interp (plev_d, tlev_d, nlev, nprof, pprof(nprof:1:-1), tprof)
        if (err /= 0) then
          call tell_error(tell_runtime_error, &
                          "read_synth_met_data: interpolation error", errstat)
          return
        endif
        ! since we reversed the pprof() array, we now need to reverse the tprof()
        ! array to provide the expected order on output.
        tprof = tprof(nprof:1:-1)
      else
        ! If no array of isobars was given, we just return the temperature
        ! profile on the native grid
        if (size(tprof).lt.nlev) then
          call tell_error(tell_runtime_error, &
                          "read_synth_met_data: array size mismatch", errstat)
          return
        endif
        tprof(:) = real(tmp_tprof(1,1,nlev:1:-1), kind=8)
      endif
    endif

  end subroutine read_synth_met_data

end module met_module
