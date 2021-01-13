!> Read OMI cloud pressure climatology and interpolate
module m_read_pres_clim

  use tell_module

  public get_tomsv8_ctp
  private get_gridfrac1

contains

  !> Get TOMS V8 monthly mean cloud for one pixel
  !------------------------------------------------------------------------
  !
  ! @param[in]   month    month of observation (1-12)
  ! @param[in]   day      day of observation within month (1-31)
  ! @param[in]   lon      longitude of pixel
  ! @param[in]   lat      latitude of pixel
  ! @param[out]  ctp      cloud top pressure (hPa)
  ! @param       errstat  error tracking integer, non-zero = [roblem
  !
  ! @author  E. O'Sullivan  Jan 2021
  !  taken from O3P code with only cosmetic changes
  !------------------------------------------------------------------------
  subroutine get_tomsv8_ctp(month, day, lon, lat, ctp, errstat)

    use m_LUN_set, only: pres_clim_LUN
    use m_pgs_include

    implicit none

    ! Input variables
    integer (kind=4), intent(in) :: month, day
    real (kind=8), intent(in) :: lon, lat
    ! Output variables
    integer (kind=4), intent(inout) :: errstat
    real (kind=8), intent(out) :: ctp
    ! Local variables
    integer (kind=4), parameter :: nlat=180, nlon=360, nmon=12, &
         atmos_unit=7771777
    real (kind=8), parameter :: longrid = 1.0, latgrid = 1.0, mongrid=1.0, &
         lon0=-180.0, lat0=-90.0, mon0=0.0
    character (len=128) :: isccp_fname
    integer (kind=4), dimension(2) :: latin, lonin, monin
    real (kind=8), dimension(2) :: latfrac, lonfrac, monfrac
    real (kind=8), save, dimension(:,:,:), pointer :: glbctp
    integer (kind=4) :: i, j, k, nblat, nblon, nbmon
    real (kind=8) :: mon
    logical :: file_exist
    logical, save :: first = .true.
    integer :: status, version, pgs_pc_getreference


    if (errstat /= 0) return

    if (first) then
      allocate(glbctp(nmon, nlon, nlat), stat = errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, "get_tomsv8_ctp: allocate failed",&
             errstat)
        return
      endif

      ! get cloud pressure climatology filename
      version = 1
      status = pgs_pc_getreference ( pres_clim_LUN, version, isccp_fname)
      if (status /= 0) then
        call tell_error (tell_io_read_error, &
             "get_tomsv8_ctp: failed to get pressure climatology filename", &
             errstat)
        return
      endif
      !isccp_fname = "../refdata/data/omcldrr_clim_rev.dat"

      ! Determine if file exists or not
      inquire (FILE= isccp_fname, EXIST= file_exist)
      if (.not. file_exist) then
        call tell_error (tell_runtime_error, &
             "get_tomsv8_ctp: pressure climatology not found", errstat)
        return
      endif
      ! file exists, so read it
      open (UNIT = atmos_unit, file=isccp_fname , status = 'unknown')
      do i = 1, nmon
        do j = nlat, 1, -1
          read (atmos_unit, *) glbctp(i, :, j)
        enddo
      enddo
      close (atmos_unit)
      first = .false.
    endif

    mon = (month - 1.0) + day / 31.0

    call get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, &
         lon0, lat0, mon0, lon, lat, mon, nblon, nblat, nbmon, &
         lonfrac, latfrac, monfrac, lonin, latin, monin)
    ctp = 0.0
    do i = 1, nblon
      do j = 1, nblat
        do k = 1, nbmon
          ctp = ctp + glbctp(monin(k), lonin(i), latin(j)) * &
               lonfrac(i) * latfrac(j) * monfrac(k)
        enddo
      enddo
    enddo

    return

  end subroutine get_tomsv8_ctp




  !> Interpolate pressure climatology in lat, lon, time
  subroutine get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, &
       lon0, lat0, mon0, lon, lat, mon, nblon, nblat, nbmon, &
       lonfrac, latfrac, monfrac,lonin, latin, monin)

  implicit none

  ! Input variables
  integer (kind=4), intent(in) :: nlon, nlat, nmon
  real (kind=8), intent(in) :: lon0, lat0, mon0, lat,lon, mon, longrid, &
       latgrid, mongrid
  ! Output variables
  integer (kind=4), intent(out) :: nblon, nblat, nbmon
  integer (kind=4), dimension(2), intent(out) :: latin, lonin, monin
  real (kind=8), dimension(2), intent(out) :: latfrac, lonfrac, monfrac
  ! Local variables
  real (kind=8) :: frac, lat_offset, lon_offset, mon_offset

  lat_offset   = lat0   + latgrid / 2.0
  lon_offset   = lon0   + longrid / 2.0
  mon_offset   = mon0   + mongrid / 2.0

  nblat = 2
  frac = (lat - lat_offset) / latgrid + 1
  latin(1) = int(frac)
  latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac
  latfrac(2) = 1.0 - latfrac(1)
  if (latin(1) == 0)   then
     latin(1) = 1
     latfrac(1) = 1.0
     nblat = 1
  endif

  if (latin(2) > nlat) then
     latin(1) = nlat
     latfrac(1) = 1.0
     nblat = 1
  endif

  ! Circular in longitude direction
  nblon = 2
  frac = (lon - lon_offset) / longrid + 1
  lonin(1) = int(frac)
  lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac
  lonfrac(2) = 1.0 - lonfrac(1)
  if (lonin(1) == 0) lonin(1) = nlon
  if (lonin(2) > nlon) lonin(2) = 1

  ! Circular in year
  nbmon = 2
  frac = (mon - mon_offset) / mongrid + 1
  monin(1) = int(frac)
  monin(2) = monin(1) + 1
  monfrac(1) = monin(2) - frac
  monfrac(2) = 1.0 - monfrac(1)
  if (monin(1) == 0) monin(1) = nmon
  if (monin(2) > nmon) monin(2) = 1

  return

  end subroutine get_gridfrac1



end module m_read_pres_clim
