module polpredict_module
  use netcdf
  use tell_module
  use tio_module
  use types_module
  use ndinterp_module
  implicit none
  private

  public pp_read, pp_dealloc, pp_initialized, pp_get_swav, pp_get_wav
  public pp_interp_ozone_pre, pp_interp_ctp, pp_interp_surface_albedo
  public pp_interp_ozone_profile, pp_ozone_zone_info, pp_derive_scene_pressure
  public pp_derive_ctp, pp_derive_to3, pp_derive_albcld, pp_get_qu

  type, public :: polpredict_type
    private
    type (ndi_dim) :: qu_dims(7)  ! (wav,alb,pre,ozo,sza,vza,raa)
    real (kind=r4), allocatable, dimension(:,:,:, :,:,:, :) :: q
    real (kind=r4), allocatable, dimension(:,:,:, :,:,:, :) :: u
    real (kind=r4), allocatable, dimension(:,:) :: ozcol          ! ozcol(ozo,pre)
    type (ndi_dim) :: si_dims(7)  ! (swav,alb,pre,ozo,sza,vza,raa)
    real (kind=r4), allocatable, dimension(:,:,:, :,:,:, :) :: si
    integer (kind=i2), allocatable, dimension(:) :: salbflg       ! salbflg(swav)
    type (ndi_dim) :: ll_dims(2)  ! (lat,lon)
    real (kind=r4), allocatable, dimension(:,:) :: ps
    real (kind=r4), allocatable, dimension(:,:) :: tom
    real (kind=r4), allocatable, dimension(:,:) :: toz_surf
    real (kind=r4), allocatable, dimension(:,:) :: ctp
    type (ndi_dim) :: ler_dims(3)  ! (albwav,lat,lon)
    real (kind=r4), allocatable, dimension(:,:,:) :: ler
  end type

  ! array indices for polpredict_type % qu_dims(:)
  integer, parameter :: &
    iqu_wav=1, iqu_alb=2, iqu_pre=3, &
    iqu_ozo=4, iqu_sza=5, iqu_vza=6, iqu_raa=7

  ! array indices for polpredict_type % ll_dims(:)
  integer, parameter :: ll_lat=1, ll_lon=2

  ! array indices for polpredict_type % ler_dims(:)
  integer, parameter :: ler_albwav=1, ler_lat=2, ler_lon=3

  integer, parameter :: maxzone = 2

  type, public :: pp_ozone_zone_info_type
    integer :: nzone, beg_index, end_index
    integer, dimension(maxzone) :: nzoneo3
    real (kind=r8), dimension(maxzone) :: zonefracs
  end type

contains

  subroutine pp_interp_ozone_profile (lut_s, pre, errstat)
    implicit none
    type (polpredict_type), target, intent(inout) :: lut_s
    real (kind=r8), intent(in) :: pre
    integer, intent(inout) :: errstat

    integer, parameter :: nd=1, nterms=2**nd
    integer, dimension(nd) :: indices, id_k
    real (kind=r8), dimension(nd) :: weights
    real (kind=r8), pointer, dimension(:) :: oz
    real (kind=r8) :: wt_k
    integer :: k, nz

    if (errstat /= 0) return

    nz = lut_s % qu_dims(iqu_ozo) % dimlen
    oz => lut_s % qu_dims(iqu_ozo) % x

    ! interpolate ozone profile vector, oz(nz)

    call ndi_find_indices (lut_s % qu_dims (iqu_pre:iqu_pre), (/pre/), indices)
    call ndi_calc_weights (lut_s % qu_dims (iqu_pre:iqu_pre), (/pre/), indices, weights)

    oz(:) = 0.0
    do k = 1, nterms
      call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
      oz(:) = oz(:) + wt_k * lut_s % ozcol(:, id_k(1))
    enddo

    ! Update the SI copy of the ozo dimension array.
    ! This duplication is unfortunate and error prone, but
    ! fortran pointers are crippled and I don't see a better solution.
    lut_s % si_dims (iqu_ozo) % x(:) = oz(:)

  end subroutine pp_interp_ozone_profile

  subroutine pp_interp_ozone_pre (lut_s, lat0, lon0, pre, oz0_r8, errstat)
    implicit none
    type (polpredict_type), intent(in) :: lut_s
    real (kind=r8), intent(in) :: lat0, lon0, pre
    real (kind=r8), intent(out) :: oz0_r8
    integer, intent(inout) :: errstat

    real (kind=r8), dimension(2) :: x
    real (kind=r4) :: tom0, pre0, oz0

    if (errstat /= 0) return

    x(:) = (/lat0, lon0/)

    ! ozone value for pre0
    call ndi_table_interp (lut_s % ll_dims, lut_s % toz_surf, x, oz0, errstat)
    call ndi_table_interp (lut_s % ll_dims, lut_s % ps, x, pre0, errstat)
    call ndi_table_interp (lut_s % ll_dims, lut_s % tom, x, tom0, errstat)
    if (errstat /= 0) return

    ! Adjust ozo0 for surface pressure difference, get ozone for pre
    oz0_r8 = oz0 + tom0 * (pre - pre0) / 1266.56d0

  end subroutine pp_interp_ozone_pre

  subroutine pp_ozone_zone_info (lut_s, lat, oz, info, errstat)
    implicit none
    type (polpredict_type), intent(in) :: lut_s
    real (kind=r8), intent(in) :: lat, oz
    type (pp_ozone_zone_info_type), intent(inout) :: info
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    call define_o3_zones (lut_s % qu_dims(iqu_ozo) % x, oz, lat, maxzone, &
                          info % nzone, info % nzoneo3, &
                          info % beg_index, info % end_index, &
                          info % zonefracs)

  end subroutine pp_ozone_zone_info

  subroutine pp_interp_ctp (lut_s, lat0, lon0, ctp_r8, errstat)
    implicit none
    type (polpredict_type), intent(in) :: lut_s
    real (kind=r8), intent(in) :: lat0, lon0
    real (kind=r8), intent(out) :: ctp_r8
    integer, intent(inout) :: errstat

    real (kind=r8), dimension(2) :: x
    real (kind=r4) :: ctp0

    if (errstat /= 0) return

    x(:) = (/lat0, lon0/)

    call ndi_table_interp (lut_s % ll_dims, lut_s % ctp, x, ctp0, errstat)
    if (errstat /= 0) return

    ctp_r8 = real(ctp0, kind=r8)

  end subroutine pp_interp_ctp

  subroutine pp_interp_surface_albedo (lut_s, lat0, lon0, swav, snalbs, errstat)
    implicit none
    type (polpredict_type), intent(in), target :: lut_s
    real (kind=r8), intent(in) :: lat0, lon0
    real (kind=r8), dimension(:), intent(in) :: swav
    real (kind=r8), dimension(:), intent(out) :: snalbs
    integer, intent(inout) :: errstat

    integer, parameter :: nd = 2, nterms = 2**(nd)
    integer, dimension(nd) :: indices, id_k
    real (kind=r8), dimension(nd) :: weights, x
    real (kind=r8), allocatable, dimension(:) :: alb
    real (kind=r8) :: wt_k
    real (kind=r8), pointer, dimension(:) :: albwav_grid
    integer :: err, k, num_ler_waves, nswav, fidx, lidx

    if (errstat /= 0) return

    x(:) = (/lat0, lon0/)

    num_ler_waves = lut_s % ler_dims(ler_albwav) % dimlen

    ! allocate space for interpolated function of wavelength
    allocate (alb(num_ler_waves), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    ! We want to interpolate a vector-valued quantity
    ! so it's more efficient to use the lower level ndi_ routines.
    ! This lets us move the (implied) loop over wavelengths inside
    ! the loop over terms.

    call ndi_find_indices (lut_s % ler_dims(2:3), x, indices)
    call ndi_calc_weights (lut_s % ler_dims(2:3), x, indices, weights)

    alb(:) = 0.0
    do k = 1, nterms
      call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
      alb(:) = alb(:) + wt_k * lut_s % ler(:, id_k(1), id_k(2))
    enddo

    ! Now, spline interpolate onto the SI wavelength grid
    albwav_grid => lut_s % ler_dims(ler_albwav) % x
    fidx = minloc (swav, dim=1, mask=(swav >= albwav_grid(1)))
    lidx = maxloc (swav, dim=1, mask=(swav <= albwav_grid(num_ler_waves)))
    err = spline (albwav_grid, alb, num_ler_waves, &
                  swav(fidx:lidx), snalbs(fidx:lidx), lidx-fidx+1)
    if (err /= 0) then
      call tell_error (tell_runtime_error, &
                       "pp_interp_surface_albedo: spline failed", errstat)
      return
    endif
    nswav = size(swav)
    if (fidx > 1) snalbs(1:fidx-1) = snalbs(fidx)
    if (lidx < nswav) snalbs(lidx+1:nswav) = snalbs(lidx)

  end subroutine pp_interp_surface_albedo

  logical function pp_initialized (lut_s)
    implicit none
    type (polpredict_type), intent(in) :: lut_s
    pp_initialized = allocated (lut_s % qu_dims(1) % x)
  end function pp_initialized

  subroutine pp_dealloc (lut_s, errstat)
    implicit none
    type (polpredict_type), intent(inout) :: lut_s
    integer, intent(inout) :: errstat

    call ndi_dims_dealloc (lut_s % qu_dims)
    call ndi_dims_dealloc (lut_s % ll_dims)
    call ndi_dims_dealloc (lut_s % si_dims)
    call ndi_dims_dealloc (lut_s % ler_dims)

    if (allocated (lut_s % q)) deallocate (lut_s % q)
    if (allocated (lut_s % u)) deallocate (lut_s % u)
    if (allocated (lut_s % ozcol)) deallocate (lut_s % ozcol)
    if (allocated (lut_s % si)) deallocate (lut_s % si)
    if (allocated (lut_s % salbflg)) deallocate (lut_s % salbflg)
    if (allocated (lut_s % ps)) deallocate (lut_s % ps)
    if (allocated (lut_s % tom)) deallocate (lut_s % tom)
    if (allocated (lut_s % toz_surf)) deallocate (lut_s % toz_surf)
    if (allocated (lut_s % ctp)) deallocate (lut_s % ctp)
    if (allocated (lut_s % ler)) deallocate (lut_s % ler)

    if (errstat /= 0) return ! silence compiler warning

  end subroutine pp_dealloc

  subroutine pp_get_wav (lut_s, wav, errstat)
    implicit none
    type (polpredict_type), intent(in) :: lut_s
    real (kind=r8), allocatable, dimension(:), intent(out) :: wav
    integer, intent(inout) :: errstat
    integer :: err
    if (errstat /= 0) return
    allocate (wav(lut_s % qu_dims(iqu_wav) % dimlen), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "pp_get_wav: malloc failed", errstat)
      return
    endif
    wav(:) = lut_s % qu_dims(iqu_wav) % x(:)
  end subroutine pp_get_wav

  subroutine pp_get_swav (lut_s, ibeg, iend, swav, errstat)
    implicit none
    type (polpredict_type), intent(in) :: lut_s
    integer, intent(in) :: ibeg, iend
    real (kind=r8), allocatable, dimension(:), intent(out) :: swav
    integer, intent(inout) :: errstat
    integer :: num_swav, max_num_swav, err

    if (errstat /= 0) return

    max_num_swav = lut_s % si_dims(iqu_wav) % dimlen

    num_swav = iend - ibeg + 1
    if (num_swav <= 0 .or. num_swav > max_num_swav) then
      call tell_error (tell_runtime_error, &
                       "pp_get_swav: invalid wavelength selection", errstat)
      return
    endif

    allocate (swav(num_swav), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "pp_get_swav: malloc failed", errstat)
      return
    endif

    swav(:) = lut_s % si_dims(iqu_wav) % x(ibeg:iend)

  end subroutine pp_get_swav

  subroutine pp_read (obj, subset, fmonth, lut, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (radiance_subset_type), intent(inout) :: subset
    real (kind=r4), intent(in) :: fmonth
    type (polpredict_type), intent(out) :: lut
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    call read_dim (obj, "sza", lut % qu_dims(iqu_sza), errstat, subset % sza)
    call read_dim (obj, "vza", lut % qu_dims(iqu_vza), errstat, subset % vza)
    call read_dim (obj, "raa", lut % qu_dims(iqu_raa), errstat, subset % raa)
    call read_dim (obj, "alb", lut % qu_dims(iqu_alb), errstat)
    call read_dim (obj, "pre", lut % qu_dims(iqu_pre), errstat)
    call read_dim (obj, "wav", lut % qu_dims(iqu_wav), errstat, subset % wav)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "reading QU lookup table dimensions", errstat)
      return
    endif

    ! The 'ozo' coordinate grid depends on the surface pressure
    ! through a separate lookup table: ozcol(ozo,pre)
    call init_ozo_dim (obj, "ozo", lut % qu_dims(iqu_ozo), errstat)
    call read_ozcol_array (obj, &
                           lut % qu_dims(iqu_ozo) % dimlen, &
                           lut % qu_dims(iqu_pre) % dimlen, &
                           lut % ozcol, errstat)
    if (errstat /= 0) return

    ! SI array has wavelength dimension 'swav',
    ! while Q,U arrays have wavelength dimension 'wav'
    call dup_dims (lut % qu_dims, lut % si_dims, errstat)
    deallocate (lut % si_dims(iqu_wav) % x)
    call read_dim (obj, "swav", lut % si_dims(iqu_wav), errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "reading SI lookup table dimension", errstat)
      return
    endif

    call read_dim (obj, "lon", lut % ll_dims(ll_lon), errstat, subset % lon)
    call read_dim (obj, "lat", lut % ll_dims(ll_lat), errstat, subset % lat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "reading QU lookup table dimensions", errstat)
      return
    endif

    call read_dim (obj, "albwav", lut % ler_dims(1), errstat)
    call dup_dims (lut % ll_dims, lut % ler_dims(2:3), errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "reading ler lookup table dimensions", errstat)
      return
    endif

    call read_qu_arrays (obj, subset, lut % qu_dims, &
                         lut % q, &
                         lut % u, &
                         errstat)
    if (errstat /= 0) return

    call read_swav_arrays (obj, subset, lut % si_dims, &
                           lut % si, &
                           lut % salbflg, &
                           errstat)
    if (errstat /= 0) return

    call read_ll_arrays (obj, subset, fmonth, lut % ll_dims, lut % ler_dims, &
                         lut % ps, &
                         lut % toz_surf, &
                         lut % tom, &
                         lut % ctp, &
                         lut % ler, &
                         errstat)
    if (errstat /= 0) return

  end subroutine pp_read

  subroutine read_ll_arrays (obj, subset, fmonth, ll_dims, ler_dims, &
                             ps, toz_surf, tom, ctp, ler, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (radiance_subset_type), intent(in) :: subset
    real (kind=r4), intent(in) :: fmonth
    type (ndi_dim), dimension(:), intent(inout) :: ll_dims, ler_dims
    real (kind=r4), allocatable, dimension(:,:), intent(out) :: &
      ps, toz_surf, tom, ctp
    real (kind=r4), allocatable, dimension(:,:,:), intent(out) :: ler
    integer, intent(inout) :: errstat

    real (kind=r4), allocatable, dimension(:,:,:) :: tmp0_3d, tmp1_3d
    real (kind=r4), allocatable, dimension(:,:) :: tmp0, tmp1
    real (kind=r4) :: weight
    integer, dimension(4) :: start, edge
    integer :: imon0, imon1, err, num_lon, num_lat, num_albwav

    if (errstat /= 0) return

    ! The month coordinate, fmonth, satisfies 0.0 <= fmonth <= 12.0.
    ! The monthly tables give the value of a quantity at mid-month.

    if (1.0 < fmonth .and. fmonth < 12.0) then
      imon0 = int(fmonth)
      imon1 = imon0 + 1
    else
      imon0 = 12
      imon1 = 1
    endif

    weight = fmonth - int(fmonth)

    num_lon = ll_dims(ll_lon) % dimlen
    num_lat = ll_dims(ll_lat) % dimlen
    num_albwav = ler_dims(ler_albwav) % dimlen

    allocate (tmp0(num_lat, num_lon), &
              tmp1(num_lat, num_lon), &
              tmp0_3d(num_albwav, num_lat, num_lon), &
              tmp1_3d(num_albwav, num_lat, num_lon), &
              tom(num_lat, num_lon), &
              ctp(num_lat, num_lon), &
              ps(num_lat, num_lon), &
              toz_surf(num_lat, num_lon), &
              ler(num_albwav, num_lat, num_lon), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    ! In the file, the lat index is fastest varying,
    start(1:2) = (/ &
      subset % lon % imin - 1, &
      subset % lat % imin - 1 &
      /)
    call tiof_get2d_r4 (obj, "ps", start(1:2), (/num_lon,num_lat/), ps, errstat)
    if (errstat /= 0) return

    edge(1:3) = (/1,num_lon,num_lat/)

    start(1:3) = (/ &
      imon0-1, &
      subset % lon % imin - 1, &
      subset % lat % imin - 1 &
      /)
    call tiof_get2d_r4 (obj, "tom", start(1:3), edge(1:3), tmp0, errstat)
    start(1) = imon1-1
    call tiof_get2d_r4 (obj, "tom", start(1:3), edge(1:3), tmp1, errstat)
    if (errstat /= 0) return

    tom(:,:) = (1.0 - weight) * tmp0 + weight * tmp1

    start(1) = imon0-1
    call tiof_get2d_r4 (obj, "toz_surf", start(1:3), edge(1:3), tmp0, errstat)
    start(1) = imon1-1
    call tiof_get2d_r4 (obj, "toz_surf", start(1:3), edge(1:3), tmp1, errstat)
    if (errstat /= 0) return

    toz_surf(:,:) = (1.0 - weight) * tmp0 + weight * tmp1

    start(1) = imon0-1
    call tiof_get2d_r4 (obj, "ctp", start(1:3), edge(1:3), tmp0, errstat)
    start(1) = imon1-1
    call tiof_get2d_r4 (obj, "ctp", start(1:3), edge(1:3), tmp1, errstat)
    if (errstat /= 0) return

    ctp(:,:) = (1.0 - weight) * tmp0 + weight * tmp1

    ! we're done with these now
    deallocate(tmp0, tmp1, stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "dealloc failed", errstat)
      return
    endif

    ! in the file, albwav is fastest varying, (which is good!)
    start(1:4) = (/ &
      imon0-1, &
      subset % lon % imin - 1, &
      subset % lat % imin - 1, &
      0/)
    edge(1:4) = (/1,num_lon,num_lat,num_albwav/)

    call tiof_get3d_r4 (obj, "ler", start(1:4), edge(1:4), tmp0_3d, errstat)
    start(1) = imon1-1
    call tiof_get3d_r4 (obj, "ler", start(1:4), edge(1:4), tmp1_3d, errstat)
    if (errstat /= 0) return

    ler(:,:,:) = (1.0 - weight) * tmp0_3d + weight * tmp1_3d

    deallocate (tmp0_3d, tmp1_3d, stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "dealloc failed", errstat)
      return
    endif

  end subroutine read_ll_arrays

  subroutine read_qu_arrays (obj, subset, qu_dims, q, u, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (radiance_subset_type), intent(in) :: subset
    type (ndi_dim), dimension(:), intent(in) :: qu_dims
    real (kind=r4), allocatable, dimension(:,:,:,:,:,:,:), intent(inout) :: q, u
    integer, intent(inout) :: errstat

    integer, dimension(7) :: start, edge
    integer :: err

    ! Note that these indices follow the C array index ordering
    ! (as does the file itself)
    start(:) = (/ &
      subset % raa % imin - 1, &
      subset % vza % imin - 1, &
      subset % sza % imin - 1, &
      0, 0, 0, &
      subset % wav % imin - 1/)

    edge(:) = (/ &
      qu_dims(iqu_raa) % dimlen, &
      qu_dims(iqu_vza) % dimlen, &
      qu_dims(iqu_sza) % dimlen, &
      qu_dims(iqu_ozo) % dimlen, &
      qu_dims(iqu_pre) % dimlen, &
      qu_dims(iqu_alb) % dimlen, &
      qu_dims(iqu_wav) % dimlen /)

    ! allocate space using fortran index ordering (duh):
    allocate (q(edge(7),edge(6),edge(5),edge(4),edge(3),edge(2),edge(1)), &
              u(edge(7),edge(6),edge(5),edge(4),edge(3),edge(2),edge(1)), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
    endif

    call tiof_get7d_r4 (obj, "q", start, edge, q, errstat)
    call tiof_get7d_r4 (obj, "u", start, edge, u, errstat)
    if (errstat /= 0) return

  end subroutine read_qu_arrays

  subroutine read_swav_arrays (obj, subset, si_dims, si, salbflg, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (radiance_subset_type), intent(in) :: subset
    type (ndi_dim), dimension(:), intent(in) :: si_dims
    real (kind=r4), allocatable, dimension(:,:,:,:,:,:,:), intent(inout) :: si
    integer (kind=i2), allocatable, dimension(:), intent(inout) :: salbflg
    integer, intent(inout) :: errstat

    integer, dimension(7) :: start, edge
    integer :: err

    ! Note that these indices follow the C array index ordering
    ! (as does the file itself)
    start(:) = (/ &
      subset % raa % imin - 1, &
      subset % vza % imin - 1, &
      subset % sza % imin - 1, &
      0, 0, 0, 0/)

    edge(:) = (/ &
      si_dims(iqu_raa) % dimlen, &
      si_dims(iqu_vza) % dimlen, &
      si_dims(iqu_sza) % dimlen, &
      si_dims(iqu_ozo) % dimlen, &
      si_dims(iqu_pre) % dimlen, &
      si_dims(iqu_alb) % dimlen, &
      si_dims(iqu_wav) % dimlen /)

    ! allocate space using fortran index ordering (duh):
    allocate (si(edge(7),edge(6),edge(5),edge(4),edge(3),edge(2),edge(1)), &
              salbflg(edge(7)), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
    endif

    call tiof_get7d_r4 (obj, "SI", start, edge, si, errstat)
    call tiof_get1d_i2 (obj, "salbflg", start(7:7), edge(7:7), salbflg, errstat)
    if (errstat /= 0) return

  end subroutine read_swav_arrays

  subroutine read_ozcol_array (obj, dimlen_ozo, dimlen_pre, ozcol, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    integer, intent(in) :: dimlen_ozo, dimlen_pre
    real (kind=r4), allocatable, dimension(:,:), intent(out) :: ozcol
    integer, intent(inout) :: errstat

    real (kind=r4), allocatable, dimension(:,:) :: tmp_ozcol
    integer, dimension(2) :: edge
    integer :: err

    if (errstat /= 0) return

    ! FIXME: In the netcdf file, the 'pre' dimension is fastest varying
    !        but this is not optimal. We'll read the array as-is,
    !        and then transpose it so that the 'ozo' dimension varies fastest
    allocate (tmp_ozcol (dimlen_pre, dimlen_ozo), &
              ozcol (dimlen_ozo, dimlen_pre), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    edge = (/dimlen_ozo, dimlen_pre/)
    call tiof_get2d_r4 (obj, "ozcol", (/0, 0/), edge, tmp_ozcol, errstat)
    if (errstat /= 0) return

    ozcol = transpose (tmp_ozcol)

    deallocate (tmp_ozcol, stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "dealloc failed", errstat)
      return
    endif

  end subroutine read_ozcol_array

  subroutine init_ozo_dim (obj, dim_name, the_dim, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    character (len=*), intent(in) :: dim_name
    type (ndi_dim), intent(out) :: the_dim
    integer, intent(inout) :: errstat
    integer :: dimlen, err

    if (errstat /= 0) return

    call tiof_inq_dimlen (obj, dim_name, dimlen, errstat)
    if (errstat /= 0) return

    the_dim % dimlen = dimlen
    the_dim % name = dim_name
    allocate (the_dim % x(dimlen), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif
  end subroutine init_ozo_dim

  subroutine read_dim (obj, dim_name, the_dim, errstat, subset)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    character (len=*), intent(in) :: dim_name
    type (ndi_dim), intent(out) :: the_dim
    integer, intent(inout) :: errstat
    type(range_type), optional, intent(inout) :: subset

    type(range_type) :: dim_range
    real (kind=r8), allocatable, dimension(:) :: dim_tmp
    real (kind=r8) :: dim_min, dim_max
    integer :: start(1), edge(1)
    integer :: dimlen, imin, imax, dimlen_subset, err, i

    if (errstat /= 0) return

    call tiof_inq_dimlen (obj, dim_name, dimlen, errstat)
    if (errstat /= 0) return

    allocate (dim_tmp(dimlen), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    start(1) = 0
    edge(1)  = dimlen
    call tiof_get1d_r8 (obj, dim_name, start, edge, dim_tmp, &
                        errstat)
    if (errstat /= 0) return

    if (present(subset)) then
      dim_range = subset
    else
      dim_range % min = -huge(0.0d0)
      dim_range % max =  huge(0.0d0)
      dim_range % imin = 1
      dim_range % imax = dimlen
    endif

    dim_min = dim_range % min
    dim_max = dim_range % max

    imin = 1
    imax = 1
    do i = 1, dimlen
      if (dim_min < dim_tmp(i)) exit
      imin = i
    enddo
    do i = imin, dimlen
      imax = i
      if (dim_max < dim_tmp(i)) exit
    enddo

    dim_range % imin = imin
    dim_range % imax = imax

    dimlen_subset = imax - imin + 1
    allocate (the_dim % x(dimlen_subset), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "allocate failed", errstat)
      return
    endif

    the_dim % x(:) = real (dim_tmp(imin:imax), kind=r4)
    the_dim % dimlen = dimlen_subset
    the_dim % name = dim_name

    if (present (subset)) then
      subset = dim_range
    endif

  end subroutine read_dim

  subroutine dup_dims (from_dims, to_dims, errstat)
    implicit none
    type (ndi_dim), dimension(:), intent(in) :: from_dims
    type (ndi_dim), dimension(:), intent(out) :: to_dims
    integer, intent(inout) :: errstat

    integer, dimension(size(from_dims)) :: dimlens

    integer :: i

    if (errstat /= 0) return

    do i = 1, size(from_dims)
      dimlens(i) = from_dims(i) % dimlen
    enddo

    call ndi_dims_alloc (to_dims, dimlens, errstat)
    if (errstat /= 0) return

    to_dims(:) = from_dims(:)

  end subroutine dup_dims

  subroutine pp_derive_scene_pressure (lut_s, cld_wave_index, srad, &
                                       sza, vza, raa, oz, oz_info, &
                                       snalbs, snps, ctp, errstat)
    implicit none
    type (polpredict_type), target, intent(in) :: lut_s
    integer, intent(in) :: cld_wave_index
    real (kind=r8), dimension(:), intent(in) :: srad
    real (kind=r8), intent(in) :: sza, vza, raa, oz
    type (pp_ozone_zone_info_type), intent(in) :: oz_info
    real (kind=r8), dimension(:,:), intent(in) :: snalbs
    real (kind=r8), dimension(:), intent(in) :: snps
    real (kind=r8), intent(out) :: ctp
    integer, intent(inout) :: errstat

    integer, parameter :: nd = 6, nterms=2**nd
    integer :: is, iz, k, fidx, lidx, ioz
    integer, dimension(nd) :: ids, indices, id_k
    real (kind=r8), dimension(nd) :: x, weights
    real (kind=r8), dimension(2) :: snrad
    real (kind=r8), pointer, dimension(:) :: alb_grid, pre_grid, oz_grid
    real (kind=r8) :: val, radcfrac, cfrac, wt_k

    if (errstat /= 0) return

    ids(:) = (/iqu_alb,iqu_pre,iqu_ozo,iqu_sza,iqu_vza,iqu_raa/)
    x(4:6) = (/sza, vza, raa/)

    call ndi_find_indices (lut_s % si_dims(ids(4:6)), x(4:6), indices(4:6))

    alb_grid => lut_s % si_dims(iqu_alb) % x
    pre_grid => lut_s % si_dims(iqu_pre) % x
    oz_grid => lut_s % si_dims(iqu_ozo) % x

    snrad(:) = 0.0

    do is = 1, 2
      x(1) = snalbs(cld_wave_index, is)
      indices(1) = ndi_find_index (snalbs(cld_wave_index, is), alb_grid)
      x(2) = snps(is)
      indices(2) = ndi_find_index (snps(is), pre_grid)

      fidx = oz_info % beg_index
      do iz = 1, oz_info % nzone
        lidx = fidx + oz_info % nzoneo3(iz) - 1
        ioz = ndi_find_index (oz, oz_grid(fidx:lidx))
        indices(3) = fidx + ioz - 1
        x(3) = oz

        call ndi_calc_weights (lut_s % si_dims (ids), x, indices, weights)

        val = 0.0
        do k = 1, nterms
          call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
          val = val + wt_k * lut_s % si (cld_wave_index,id_k(1),id_k(2),id_k(3),id_k(4),id_k(5),id_k(6))
        enddo
        snrad(is) = snrad(is) + val * oz_info % zonefracs(iz)
        fidx = lidx + 1
      enddo
    enddo

    cfrac = (srad(cld_wave_index) - snrad(1)) / (snrad(2) - snrad(1))
    if (cfrac <= 0.0) then
      ctp = snps(1)
    else if (cfrac > 0.0 .and. cfrac < 1.0) then
      radcfrac = cfrac + snrad(2)/srad(cld_wave_index)
      ctp = snps(1) * (1.0 - radcfrac) + snps(2) * radcfrac
    else
      ctp = snps(2)
    endif

  end subroutine pp_derive_scene_pressure

  subroutine pp_derive_ctp (lut_s, use_mler, srad, &
                            sza, vza, raa, oz, oz_info, wave_indices, &
                            nscene, ctp, cfracs, snalbs, snps, status, errstat)
    implicit none
    integer, parameter :: nw = 3
    type (polpredict_type), target, intent(in) :: lut_s
    logical, intent(in) :: use_mler
    real (kind=r8), dimension(nw), intent(in) :: srad
    real (kind=r8), intent(in) :: sza, vza, raa, oz
    type (pp_ozone_zone_info_type), intent(in) :: oz_info
    integer, dimension(nw), intent(in) :: wave_indices
    integer, intent(in) :: nscene
    real (kind=r8), intent(out) :: ctp
    real (kind=r8), dimension(:), intent(inout) :: cfracs
    real (kind=r8), dimension(:,:), intent(inout) :: snalbs
    real (kind=r8), dimension(:), intent(inout) :: snps
    integer, intent(out) :: status
    integer, intent(inout) :: errstat

    integer, parameter :: nd = 4, nterms_4d=2**4, nterms_2d=2**2
    integer :: iz, ia, ip, is, iw, k, nalb, npre, err, fidx, lidx, ioz
    integer, dimension(nd) :: ids, indices, id_k
    real (kind=r8), dimension(nd) :: x, weights
    real (kind=r8), allocatable, dimension(:,:,:) :: snrad
    real (kind=r8), allocatable, dimension(:,:) :: snrad0
    real (kind=r8), pointer, dimension(:) :: oz_grid
    real (kind=r8) :: wt_k, snrad0_val

    if (errstat /= 0) return

    ! 4 SI dimensions we'll be interpolating on:
    ids(:) = (/iqu_ozo, iqu_sza, iqu_vza, iqu_raa/)

    ! 4-D interpolation point
    x(:) = (/oz, sza, vza, raa/)

    ! locate 3 of the 4 interpolation coordinates (the geometry variables)
    ! in the main SI lookup table
    call ndi_find_indices (lut_s % si_dims (ids(2:4)), x(2:4), indices(2:4))

    ! sizes of dimensions we'll be looping over (plus nw=3)
    nalb = lut_s % si_dims (iqu_alb) % dimlen
    npre = lut_s % si_dims (iqu_pre) % dimlen

    allocate (snrad(nw, nalb, npre), &
              snrad0(nw, nscene), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    ! Initialize the 3-D object we'll be constructing --
    ! scene radiance spectrum as a function of pressure and albedo:
    snrad(:,:,:) = 0.0

    ! first index of the first zone in the ozone profile grid
    fidx = oz_info % beg_index

    oz_grid => lut_s % si_dims(iqu_ozo) % x

    do iz = 1, oz_info % nzone

      ! locate the 4th interpolation coordinate (ozone)
      ! in the derived ozone profile grid
      lidx = fidx + oz_info % nzoneo3(iz) - 1
      ioz = ndi_find_index (oz, oz_grid(fidx:lidx))
      indices(1) = fidx + ioz - 1

      ! compute the 4-D linear interpolation weights
      call ndi_calc_weights (lut_s % si_dims (ids), x, indices, weights)

      ! Using 4-D interpolation weights, linearly interpolate
      ! a 3-D object with dimensions (nw,nalb,npre) and accumulate
      ! the weighted sum over zones, iz
      do k = 1, nterms_4d
        call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
        wt_k = wt_k * oz_info % zonefracs(iz)
        do ip = 1, npre
          do ia = 1, nalb
            snrad(1:nw,ia,ip) = snrad(1:nw,ia,ip) &
            + wt_k * lut_s % si(wave_indices(1:nw), ia, ip, id_k(1), id_k(2), id_k(3), id_k(4))
          enddo
        enddo

      enddo

      fidx = lidx + 1
    enddo

    ! Interpolate on pressure and albedo to get the
    ! initial radiance spectrum for each scene, snrad0(:,:):
    ids(1:2) = (/iqu_alb, iqu_pre/)
    do is = 1, nscene
      x(2) = snps(is)
      indices(2) = ndi_find_index (x(2), lut_s % si_dims(iqu_pre) % x)
      do iw = 1, nw
        x(1) = snalbs (iw, is)
        indices(1) = ndi_find_index (x(1), lut_s % si_dims(iqu_alb) % x)
        call ndi_calc_weights (lut_s % si_dims(ids(1:2)), x(1:2), indices(1:2), weights(1:2))
        snrad0_val = 0.0
        do k = 1, nterms_2d
          call ndi_calc_term_weight (k, indices(1:2), weights(1:2), id_k, wt_k)
          snrad0_val = snrad0_val + wt_k * snrad(iw, id_k(1), id_k(2))
        enddo
        snrad0(iw,is) = snrad0_val
      enddo
    enddo

    if (use_mler) then
      call derive_ctp_mler (lut_s, srad, snrad, snrad0, &
                            snalbs, snps, cfracs, ctp, status, errstat)
    else
      call derive_ctp_ler (lut_s, srad, snrad, snrad0, &
                           snalbs, snps, ctp, status, errstat)
    endif

  end subroutine pp_derive_ctp

  subroutine derive_ctp_ler (lut_s, srad, snrad, snrad0, &
                             snalbs, snps, ctp, status, errstat)
    implicit none
    integer, parameter :: nw=3
    type (polpredict_type), target, intent(in) :: lut_s
    real (kind=r8), dimension(nw), intent(in) :: srad
    real (kind=r8), dimension(:,:,:), intent(in) :: snrad
    real (kind=r8), dimension(:,:), intent(inout) :: snrad0
    real (kind=r8), dimension(:,:), intent(inout) :: snalbs
    real (kind=r8), dimension(:), intent(inout) :: snps
    real (kind=r8), intent(out) :: ctp
    integer, intent(out) :: status
    integer, intent(inout) :: errstat

    integer, parameter :: nd = 2, nterms=2**nd
    real (kind=r8), pointer, dimension(:) :: alb_grid, pre_grid
    real (kind=r8), allocatable, dimension(:,:) :: avgrad, ratio
    real (kind=r8), dimension(nd) :: weights, x
    real (kind=r8) :: savgrad, sratio, snalb, wt_k, wt_p
    real (kind=r8) :: val, rad_a0p, rad_a1p, drad_dalb
    integer, dimension(nd) :: ids, indices, id_k
    integer :: iw, k, npre, nalb, err

    status = 0

    if (errstat /= 0) return

    npre = lut_s % si_dims (iqu_pre) % dimlen
    nalb = lut_s % si_dims (iqu_alb) % dimlen

    allocate (avgrad(nalb,npre), &
              ratio(nalb,npre), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    avgrad(:,:) = 0.5 * (snrad(1,:,:) + snrad(3,:,:))
    ratio(:,:) = snrad(2,:,:) / avgrad(:,:)

    savgrad = 0.5 * (srad(1) + srad(3))
    sratio = srad(2) / savgrad

    alb_grid => lut_s % si_dims(iqu_alb) % x
    pre_grid => lut_s % si_dims(iqu_pre) % x

    call cldret_find_xy (avgrad, ratio, alb_grid, pre_grid, &
                         nalb, npre, savgrad, sratio, snalb, ctp, status)

    if (ctp < pre_grid(1)) then
      ctp = pre_grid(1)
      status = 2
    else if (ctp > pre_grid(npre)) then
      ! surface p is not known, use pre_grid(npre)
      ctp = pre_grid(npre)
      status = 2
    endif

    if (snalb < alb_grid(1)) then
      snalb = alb_grid(1)
    else if (snalb > alb_grid(nalb)) then
      snalb = alb_grid(nalb)
    endif

    snps(1) = ctp
    snalbs(1:nw, 1) = snalb

    ids(:) = (/iqu_alb, iqu_pre/)

    do iw = 1,nw
      x(:) = (/snalbs(iw,1), snps(1)/)
      call ndi_find_indices (lut_s % si_dims(ids), x, indices)
      call ndi_calc_weights (lut_s % si_dims(ids), x, indices, weights)

      val = 0.0
      do k = 1, nterms
        call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
        val = val + wt_k * snrad(iw, id_k(1), id_k(2))
      enddo
      snrad0(iw,1) = val

      wt_p = (pre_grid(indices(2)+1) - snps(1)) / (pre_grid(indices(2)+1) - pre_grid(indices(2)))
      rad_a1p = (           wt_p  * snrad(iw, indices(1)+1, indices(2)) &
                   + (1.0 - wt_p) * snrad(iw, indices(1)+1, indices(2)+1))
      rad_a0p = (           wt_p  * snrad(iw, indices(1)  , indices(2)) &
                   + (1.0 - wt_p) * snrad(iw, indices(1)  , indices(2)+1))
      drad_dalb = (rad_a1p - rad_a0p) / (alb_grid(indices(1)+1) - alb_grid(indices(1)))
      snalbs(iw,1) = snalbs(iw,1) + (srad(iw) - snrad0(iw,1)) / drad_dalb

      if (snalbs(iw,1) > 1.0) then
        snalbs(iw,1) = 1.0
        status = 1
      else if (snalbs(iw,1) < 0.0) then
        snalbs(iw,1) = 0.0
        status = 1
      endif
    enddo

  end subroutine derive_ctp_ler

  subroutine derive_ctp_mler (lut_s, srad, snrad, snrad0, &
                              snalbs, snps, cfracs, ctp, status, errstat)
    implicit none
    integer, parameter :: nw = 3
    type (polpredict_type), target, intent(in) :: lut_s
    real (kind=r8), dimension(nw), intent(in) :: srad
    real (kind=r8), dimension(:,:,:), intent(in) :: snrad
    real (kind=r8), dimension(:,:), intent(inout) :: snrad0
    real (kind=r8), dimension(:,:), intent(inout) :: snalbs
    real (kind=r8), dimension(:), intent(inout) :: snps
    real (kind=r8), dimension(:), intent(inout) :: cfracs
    real (kind=r8), intent(out) :: ctp
    integer, intent(out) :: status
    integer, intent(inout) :: errstat

    integer, parameter :: nfc0 = 11, nd=2, nterms=2**nd
    integer :: ia, ip, iw, nalb, npre, ja, err, k
    real (kind=r8), parameter, dimension(nfc0) :: fc0 = (/(ia*0.1, ia=0,nfc0-1)/)
    real (kind=r8), dimension(nd) :: x, weights
    real (kind=r8), pointer, dimension(:) :: alb_grid, pre_grid
    real (kind=r8), allocatable, dimension(:,:,:) :: mlerad
    real (kind=r8), allocatable, dimension(:,:) :: mavgrad, mratio
    real (kind=r8) :: a_iw, wt_a, cldrad, savgrad, sratio, cfrac, val, wt_k
    real (kind=r8) :: wt_p, rad_a1p, rad_a0p, drad_dalb
    integer, dimension(nd) :: ids, indices, id_k

    status = 0

    if (errstat /= 0) return

    npre = lut_s % si_dims (iqu_pre) % dimlen
    nalb = lut_s % si_dims (iqu_alb) % dimlen

    allocate (mlerad(nfc0,npre,nw), &
              mavgrad(nfc0,npre), &
              mratio(nfc0,npre), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    alb_grid => lut_s % si_dims (iqu_alb) % x

    do iw = 1, nw
      a_iw = snalbs(iw,2)
      ja = ndi_find_index (a_iw, alb_grid)
      wt_a = (alb_grid(ja+1) - a_iw) / (alb_grid(ja+1) - alb_grid(ja))

      do ip = 1, npre
        cldrad = wt_a * snrad(iw,ja,ip) + (1.0 - wt_a) * snrad(iw,ja+1,ip)
        do ia = 1, nfc0
          mlerad(ia, ip, iw) = snrad0(iw, 1) * (1.0 - fc0(ia)) + cldrad * fc0(ia)
        enddo
      enddo
    enddo

    mavgrad(:,:) = 0.5 * (mlerad(:,:,1) + mlerad(:,:,3))
    mratio(:,:) = mlerad(:,:,2) / mavgrad(:,:)

    deallocate(mlerad)

    savgrad = 0.5 * (srad(1) + srad(3))
    sratio = srad(2) / savgrad

    pre_grid => lut_s % si_dims(iqu_pre) % x

    call cldret_find_xy(mavgrad, mratio, fc0, pre_grid, &
                        nfc0, npre, savgrad, sratio, cfrac, ctp, status)

    if (ctp < pre_grid(1)) then
      ctp = pre_grid(1)
      status = 2
    else if (ctp > snps(1)) then
      ! cloud should not be below the surface
      ctp = snps(1)
      status = 2
    endif

    snps(2) = ctp
    cfracs(1:nw) = cfrac

    ids(:) = (/iqu_alb, iqu_pre/)

    do iw = 1, nw

      x(:) = (/snalbs(iw,2), snps(2)/)
      call ndi_find_indices (lut_s % si_dims(ids), x, indices)
      call ndi_calc_weights (lut_s % si_dims(ids), x, indices, weights)

      val = 0.0
      do k = 1, nterms
        call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
        val = val + wt_k * snrad(iw, id_k(1), id_k(2))
      enddo
      snrad0(iw,2) = val

      ! Derive effective cloud fraction using MLER model
      cfracs(iw) = (srad(iw) - snrad0(iw,1)) / (snrad0(iw,2) - snrad0(iw,1))

      if (cfracs(iw) > 1.0) then
        ! If cloud fraction is  >1, re-derive CLOUD albedo
        cfracs(iw) = 1.0
        wt_p = (pre_grid(indices(2)+1) - snps(2)) / (pre_grid(indices(2)+1) - pre_grid(indices(2)))
        rad_a1p = (         wt_p  * snrad(iw, indices(1)+1, indices(2)) &
                   + (1.0 - wt_p) * snrad(iw, indices(1)+1, indices(2)+1))
        rad_a0p = (         wt_p  * snrad(iw, indices(1)  , indices(2)) &
                   + (1.0 - wt_p) * snrad(iw, indices(1)  , indices(2)+1))
        drad_dalb = (rad_a1p - rad_a0p) / (alb_grid(indices(1)+1) - alb_grid(indices(1)))
        snalbs(iw,2) = snalbs(iw,2) + (srad(iw) - snrad0(iw,2)) / drad_dalb
        if (snalbs(iw,2) > 1.0) then
          snalbs(iw,2) = 1.0
          status = 1
        endif
      else if (cfracs(iw) < 0.0) then
        ! If cloud fraction is  <0, re-derive SURFACE albedo
        cfracs(iw) = 0.0
        indices(1) = ndi_find_index (snalbs(iw,1), alb_grid)
        indices(2) = ndi_find_index (snps(1), pre_grid)
        wt_p = (pre_grid(indices(2)+1) - snps(1)) / (pre_grid(indices(2)+1) - pre_grid(indices(2)))
        rad_a1p = (         wt_p  * snrad(iw, indices(1)+1, indices(2)) &
                   + (1.0 - wt_p) * snrad(iw, indices(1)+1, indices(2)+1))
        rad_a0p = (         wt_p  * snrad(iw, indices(1)  , indices(2)) &
                   + (1.0 - wt_p) * snrad(iw, indices(1)  , indices(2)+1))
        drad_dalb = (rad_a1p - rad_a0p) / (alb_grid(indices(1)+1) - alb_grid(indices(1)))
        snalbs(iw,1) = snalbs(iw,1) + (srad(iw) - snrad0(iw,1)) / drad_dalb
        if (snalbs(iw,1) < 0.0) then
          snalbs(iw,1) = 0.0
          status = 1
        endif
      endif

    enddo

  end subroutine derive_ctp_mler

  subroutine pp_derive_to3 (lut_s, use_mler, swav, srad, &
                            sza, vza, raa, oz, oz_info, wave_indices, &
                            nscene, cfracs, snalbs, snps, status, errstat)
    implicit none
    integer, parameter :: nw = 3
    type (polpredict_type), target, intent(in) :: lut_s
    logical, intent(in) :: use_mler
    real (kind=r8), dimension(nw), intent(in) :: swav, srad
    real (kind=r8), intent(in) :: sza, vza, raa
    real (kind=r8), intent(inout) :: oz
    type (pp_ozone_zone_info_type), target, intent(in) :: oz_info
    integer, dimension(nw), intent(in) :: wave_indices
    integer, intent(in) :: nscene
    real (kind=r8), dimension(:), intent(out) :: cfracs
    real (kind=r8), dimension(:,:), intent(inout) :: snalbs
    real (kind=r8), dimension(:), intent(in) :: snps
    integer, intent(out) :: status
    integer, intent(inout) :: errstat

    real (kind=r8), parameter :: min_vis_wave = 520.0

    integer, parameter :: nd = 4, nterms_2d = 2**2, nterms_4d = 2**4
    integer :: is, ia, iw, iz, izb, ize, err, k, nalb, noz
    integer :: nzone, off, fidx, lidx, sidx, eidx, nw1, num_iter
    integer, dimension(nd) :: id_k, indices, ids
    integer, pointer, dimension(:) :: nzoneo3
    real (kind=r8), dimension(nd) :: weights, x
    real (kind=r8), allocatable, dimension(:,:,:,:) :: snalbrad
    real (kind=r8), allocatable, dimension(:,:) :: snrad, snalbwfs, snozwfs
    real (kind=r8), allocatable, dimension(:) :: oz_grid_copy
    real (kind=r8), pointer, dimension(:) :: alb_grid, oz_grid, zonefracs
    real (kind=r8), dimension(nw) :: cfracwfs, mrad, rad, albcldwfs, ozwfs
    real (kind=r8) :: wt_k, ra, ra1, wt, wscl, wt_a, wt_z, val, avg, da, do3, dr
    real (kind=r8) :: dratio, dratioda, dratiodo3, drda, drdo3, mavg, mratio, ratio
    real (kind=r8) :: snozwfs_z1, snozwfs_z, snalbwfs_a1, snalbwfs_a

    logical :: converge
    ! Use log(radiance) in improve the linearity of the O3 retrieval in UV
    ! In VIS, disable it due to the complexity of implemention
    logical :: use_lograd = .true.

    if (errstat /= 0) return

    if (swav(1) > min_vis_wave) use_lograd = .false.

    nalb = lut_s % si_dims(iqu_alb) % dimlen
    noz  = lut_s % si_dims(iqu_ozo) % dimlen

    izb = oz_info % beg_index
    ize = oz_info % end_index

    allocate (snalbrad(nw,noz,nalb,nscene), &
              snrad(nscene,nw), &
              snalbwfs(nscene,nw), &
              snozwfs(nscene,nw), &
              oz_grid_copy(noz), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    status = 0
    snalbrad(:,:,:,:) = 0
    albcldwfs(:) = 0

    ! 4 SI dimensions we'll be interpolating on:
    ids(:) = (/iqu_pre, iqu_sza, iqu_vza, iqu_raa/)

    do is = 1, nscene

      ! 4-D interpolation point
      x(:) = (/snps(is), sza, vza, raa/)

      call ndi_find_indices (lut_s % si_dims(ids), x, indices)
      call ndi_calc_weights (lut_s % si_dims(ids), x, indices, weights)

      do k = 1, nterms_4d
        call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)

        do ia = 1, nalb
          do iw = 1, nw
            do iz = izb, ize
              snalbrad(iw,iz,ia,is) = snalbrad(iw,iz,ia,is) &
                + wt_k * lut_s % si(wave_indices(iw), ia, id_k(1), iz, id_k(2), id_k(3), id_k(4))
              ! radiance at last UV wavelength independent of ozone
              if (swav(1) < min_vis_wave .and. iw == nw) exit
            enddo
          enddo
        enddo !ia
      enddo !k
    enddo !is

    ! Derive scene albedo or cloud fraction from last UV wavelength that
    ! is insensitive to ozone, and apply it to shorter wavelengths
    ! Or use last vis wavelength to initialize cloud fraction/scene albedo
    alb_grid => lut_s % si_dims(iqu_alb) % x
    iw = nw
    if (use_mler) then
      do is = 1, nscene
        ia = ndi_find_index (snalbs(iw,is), alb_grid)
        wt = (alb_grid(ia+1) - snalbs(iw,is)) / (alb_grid(ia+1) - alb_grid(ia))
        ra  = snalbrad(iw, izb, ia,   is)
        ra1 = snalbrad(iw, izb, ia+1, is)
        snrad(is,iw) = wt  * ra + (1.0 - wt) * ra1
        snalbwfs(is,iw) = (ra1 - ra)/(alb_grid(ia+1) - alb_grid(ia))
      enddo

      ! Derive effective cloud fraction using MLER model
      cfracs(iw) = (srad(iw) - snrad(1,iw)) / (snrad(2,iw) - snrad(1,iw))
      cfracwfs(iw) = snrad(2,iw) - snrad(1,iw)
      cfracs(1:iw-1) = cfracs(iw)

      if (cfracs(iw) > 1.0) then
        cfracs(iw) = 1.0
        snalbs(iw,2) = snalbs(iw,2) + (srad(iw) - snrad(2,iw)) / snalbwfs(2, iw)
        snalbs(1:iw-1,2) = snalbs(iw,2)
        cfracs(1:iw-1) = cfracs(iw)
      else if (cfracs(iw) < 0.0) then
        cfracs(iw) = 0.0
        snalbs(iw,1) = snalbs(iw,1) + (srad(iw) - snrad(1, iw)) / snalbwfs(1, iw)
        snalbs(1:iw-1,1) = snalbs(iw,1)
        cfracs(1:iw-1) = cfracs(iw)
      endif
    else
      fidx = maxloc(snalbrad(iw,izb,1:nalb,1), dim=1, &
                    mask=(snalbrad(iw,izb,1:nalb,1) < srad(iw)))
      if (fidx < 1) then
        fidx = 1
      else if (fidx == nalb) then
        fidx = nalb-1
      endif
      lidx = fidx + 1
      snalbwfs(1,iw) = ((snalbrad(iw,izb,lidx,1) - snalbrad(iw,izb,fidx,1)) &
                        / (alb_grid(lidx)- alb_grid(fidx)))
      snalbs(iw,1) = (alb_grid(fidx) &
                      + (srad(iw) - snalbrad(iw,izb,fidx,1))/snalbwfs(1, iw))
      snalbs(1:iw-1,1) = snalbs(iw,1)
    endif
    if (any(snalbs(1:nw,1:nscene) < 0.0) .or. &
        any(snalbs(1:nw,1:nscene) > 1.0) ) status = 1
    where(snalbs(1:nw,1:nscene) < 0.0)
      snalbs(1:nw,1:nscene) = 0.0
    endwhere
    where(snalbs(1:nw,1:nscene) > 1.0)
      snalbs(1:nw,1:nscene) = 1.0
    endwhere

    ! Weight ozone profiles according to latitude bands
    ! if ozone profiles from two latitude bands are selected
    nzoneo3 => oz_info % nzoneo3
    nzone = oz_info % nzone

    oz_grid => lut_s % si_dims(iqu_ozo) % x
    oz_grid_copy(:) = oz_grid(:)

    zonefracs => oz_info % zonefracs

    if (nzone == 2) then
      fidx = ize - nzoneo3(nzone) + 1
      if (nzoneo3(2) == 10) then ! M+H
        off = 2
      else
        off = 0
      endif
      fidx = fidx + off
      lidx = fidx + nzoneo3(1) - 1
      sidx = izb
      eidx = izb + nzoneo3(1) - 1

      oz_grid_copy(fidx:lidx) = (oz_grid(sidx:eidx) * zonefracs(1) &
                                 + oz_grid(fidx:lidx) * zonefracs(2))
      snalbrad(1:nw, fidx:lidx, 1:nalb, 1:nscene) = &
        snalbrad(1:nw, sidx:eidx, 1:nalb, 1:nscene) * zonefracs(1) + &
        snalbrad(1:nw, fidx:lidx, 1:nalb, 1:nscene) * zonefracs(2)

      fidx = fidx - off
      lidx = ize
    else
      fidx = izb
      lidx = ize
    endif

    if (swav(1) < min_vis_wave) then
      nw1 = nw - 1
      wscl = (swav(3) - swav(1)) / (swav(3) - swav(2))
    else
      nw1 = nw
      wscl = 1.0
    endif

    if (use_lograd .and. swav(1) < min_vis_wave) then
      mrad(1:nw1) =  log(srad(1:nw1))
      snalbrad(1:nw1, fidx:lidx, 1:nalb, :) = log(snalbrad(1:nw1, fidx:lidx, 1:nalb, :))
    else
      mrad(1:nw1) = srad(1:nw1)
    endif

    indices(:) = 0
    weights(:) = 0.0

    ! Derive delta_oz based on 317/333 nm ratio and delta_albcld at 333 nm
    converge = .false.
    num_iter = 5
    do while (.not.converge .and. num_iter > 0)

      iz = ndi_find_index (oz, oz_grid_copy(fidx:lidx))
      iz = iz + fidx - 1
      wt_z = (oz_grid_copy(iz+1) - oz) / (oz_grid_copy(iz+1) - oz_grid_copy(iz))
      indices(1) = iz
      weights(1) = wt_z

      do is = 1, nscene
        do iw = 1, nw1
          ia = ndi_find_index (snalbs(iw,is), alb_grid)
          wt_a = (alb_grid(ia+1) - snalbs(iw,is)) / (alb_grid(ia+1) - alb_grid(ia))
          indices(2) = ia
          weights(2) = wt_a

          ! 2D linear interpolation:
          val = 0.0
          do k = 1, nterms_2d
            call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
            val = val + wt_k * snalbrad(iw, id_k(1), id_k(2), is)
          enddo
          snrad(is,iw) = val

          snozwfs_z1 = (         wt_a  * snalbrad(iw, iz+1, ia  , is) &
                        + (1.0 - wt_a) * snalbrad(iw, iz+1, ia+1, is))
          snozwfs_z  = (         wt_a  * snalbrad(iw, iz  , ia  , is) &
                        + (1.0 - wt_a) * snalbrad(iw, iz  , ia+1, is))
          snozwfs(is,iw) = (snozwfs_z1 - snozwfs_z) / (oz_grid_copy(iz+1) - oz_grid_copy(1))

          snalbwfs_a1 = (         wt_z  * snalbrad(iw, iz  , ia+1, is) &
                         + (1.0 - wt_z) * snalbrad(iw, iz+1, ia+1, is))
          snalbwfs_a  = (         wt_z  * snalbrad(iw, iz  , ia  , is) &
                         + (1.0 - wt_z) * snalbrad(iw, iz+1, ia  , is))
          snalbwfs(is,iw) = (snalbwfs_a1 - snalbwfs_a) / (alb_grid(ia+1) - alb_grid(ia))
        enddo ! iw
      enddo ! is

      if (use_mler) then
        do iw = 1, nw1
          if (.not. use_lograd .or. swav(1) > min_vis_wave) then
            rad(iw) = snrad(1, iw) * (1.0 - cfracs(iw)) + snrad(2, iw) * cfracs(iw)
            ozwfs(iw) = snozwfs(1, iw) * (1.0 -cfracs(iw)) + snozwfs(2, iw) * cfracs(iw)
            cfracwfs(1:nw1) = snrad(2, 1:nw1) - snrad(1, 1:nw1)
          else
            rad(iw) = log(exp(snrad(1, iw)) * (1.0 - cfracs(iw)) + exp(snrad(2, iw)) * cfracs(iw))
            ozwfs(iw) = (snozwfs(1, iw) * (1.0 -cfracs(iw)) * exp(snrad(1, iw)) + &
                         snozwfs(2, iw) * cfracs(iw) * exp(snrad(2, iw)) ) / exp(rad(iw))
            cfracwfs(1:nw1) = (exp(snrad(2, 1:nw1)) - exp(snrad(1, 1:nw1))) / exp(rad(iw))
          endif
        enddo

        if (cfracs(2) < 0.0) then
          albcldwfs(1:nw1) = snalbwfs(1, 1:nw1)
        else if (cfracs(2) > 1.0) then
          albcldwfs(1:nw1) = snalbwfs(2, 1:nw1)
        else
          albcldwfs(1:nw1) = cfracwfs(1:nw1)
        endif
      else
        rad(1:nw1) = snrad(1, 1:nw1)
        albcldwfs(1:nw1) = snalbwfs(1, 1:nw1)
        ozwfs(1:nw1) = snozwfs(1, 1:nw1)
      endif

      if (swav(1) < min_vis_wave) then
        ! uv retrievals: use the ratio of first two wavelengths and
        ! radiance of second wavelength
        if (use_lograd) then
          mratio = mrad(1) - mrad(2)
          ratio = rad(1) - rad(2)
          dratio = mratio - ratio
          dratiodo3 = ozwfs(1) - ozwfs(2)
          dratioda = albcldwfs(1) * wscl - albcldwfs(2)
        else
          mratio = mrad(1) / mrad(2)
          ratio = rad(1) / rad(2)
          dratio = mratio - ratio
          dratiodo3 = (rad(2) * ozwfs(1) - rad(1) * ozwfs(2)) / (rad(2)**2)
          dratioda = (rad(2) * albcldwfs(1) * wscl - rad(1) * albcldwfs(2)) / (rad(2)**2)
        endif
        dr    = mrad(2) - rad(2)
        drdo3 = ozwfs(2)
        drda  = albcldwfs(2)
      else
        ! vis retrievals, use the average radiance of 1st and 3rd wavelengths and
        ! the ratio of radiance at second wavelength to average radiance
        mavg  = (mrad(1) + mrad(3))/2.0
        avg   = (rad(1) + rad(3)) /2.0
        dr    = mavg - avg
        drdo3 = (ozwfs(1) + ozwfs(3))/2.0
        drda  = (albcldwfs(1) + albcldwfs(3))/2.0

        mratio = mrad(2) / mavg
        ratio = rad(2) / avg
        dratio = mratio - ratio
        dratiodo3 = (avg * ozwfs(2) - rad(2) * drdo3) / avg**2
        dratioda = (avg * albcldwfs(2) - rad(2) * drda) / avg**2
      endif

      ! two equations, two variables
      ! dratio = dratiodo3 * do3 + dratioda * da
      ! dr = drdo3 * do3 + drda * da
      da = (dratio * drdo3 - dr * dratiodo3) / &
        (dratioda * drdo3 - drda * dratiodo3)
      do3 = (dratio - dratioda * da) / dratiodo3
      oz = oz + do3

      !print *, niter, ' oz = ', oz, 'do3 = ', do3, ' da = ', da

      if (use_mler) then
        if (cfracs(2) < 0.0) then
          snalbs(2, 1) = snalbs(2, 1) + da
          if (swav(1) < min_vis_wave) then
            snalbs(1, 1) = snalbs(1, 1) + da * wscl
          else
            snalbs(1, 1) = snalbs(1, 1) + da
            snalbs(3, 1) = snalbs(3, 1) + da
          endif
        else if (cfracs(2) > 1.0) then
          snalbs(2, 2) = snalbs(2, 2) + da
          if (swav(1) < min_vis_wave) then
            snalbs(1, 2) = snalbs(1, 2) + da * wscl
          else
            snalbs(1, 2) = snalbs(1, 2) + da
            snalbs(3, 2) = snalbs(3, 2) + da
          endif
        else
          cfracs(2) = cfracs(2) + da
          if (swav(1) < min_vis_wave) then
            cfracs(1) = cfracs(1) + da * wscl
          else
            cfracs(1) = cfracs(1) + da
            cfracs(3) = cfracs(3) + da
          endif
        endif
      else
        snalbs(2, 1) = snalbs(2, 1) + da
        if (swav(1) < min_vis_wave) then
          snalbs(1, 1) = snalbs(1, 1) + da * wscl
        else
          snalbs(1, 1) = snalbs(1, 1) + da
          snalbs(3, 1) = snalbs(3, 1) + da
        endif
      endif
      if (any(snalbs(1:nw1, 1:nscene) < 0.0) .or. &
          any(snalbs(1:nw1, 1:nscene) > 1.0) ) status = 1
      where(snalbs(1:nw1, 1:nscene) < 0.0)
        snalbs(1:nw1, 1:nscene) = 0.0
      endwhere
      where(snalbs(1:nw1, 1:nscene) > 1.0)
        snalbs(1:nw1, 1:nscene) = 1.0
      endwhere
      !print *, snalbs(1:nw1, 1:nscene), cfracs(1:nw1)

      ! special conditions, exit loop
      if (abs(do3) < 1.0) converge = .true.
      if (oz < oz_grid_copy(fidx)) then
        oz = oz_grid_copy(fidx)
        status = 2
        exit
      endif
      if (oz > oz_grid_copy(lidx)) then
        oz = oz_grid_copy(lidx)
        status = 2
        exit
      endif

      num_iter = num_iter - 1
      if (num_iter == 0 .and. .not. converge) status = 3
    enddo

  end subroutine pp_derive_to3

  ! Purpose: derive effective cloud fraction (for MLER model) or
  !          scene albedo (for LER model)
  subroutine pp_derive_albcld (lut_s, use_mler, swav, srad, &
                               sza, vza, raa, oz, oz_info, &
                               nscene, cfrac, snalbs, snps, errstat)
    implicit none
    type (polpredict_type), target, intent(in) :: lut_s
    logical, intent(in) :: use_mler
    real (kind=r8), dimension(:), intent(in) :: swav, srad
    real (kind=r8), intent(in) :: sza, vza, raa
    real (kind=r8), intent(inout) :: oz
    type (pp_ozone_zone_info_type), target, intent(in) :: oz_info
    integer, intent(in) :: nscene
    real (kind=r8), dimension(:), intent(out) :: cfrac
    real (kind=r8), dimension(:,:), intent(inout) :: snalbs
    real (kind=r8), dimension(:), intent(in) :: snps
    integer, intent(inout) :: errstat

    integer, parameter :: nd = 5, nterms_5d = 2**5
    real (kind=r8), dimension(nd) :: weights
    integer, dimension(nd) :: id_k, indices
    integer, dimension(3) :: ids
    real (kind=r8), dimension(3) :: xs
    real (kind=r8), allocatable, dimension(:,:) :: snalbrad
    real (kind=r8), dimension(nscene) :: snrad, snalbwfs
    real (kind=r8), pointer, dimension(:) :: pre_grid, oz_grid, alb_grid
    integer :: is, ip, ia, iz, ioz, fidx, lidx, k, nalb, err, iw, nswav
    real (kind=r8) :: wt, wt_k, cfracwf

    if (errstat /= 0) return

    nswav = size(swav)
    nalb = lut_s % si_dims(iqu_alb) % dimlen

    allocate (snalbrad(nalb,nscene), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    pre_grid => lut_s % si_dims(iqu_pre) % x
    oz_grid => lut_s % si_dims(iqu_ozo) % x

    !SI dimensions = (swav,alb, pre,ozo,sza,vza,raa)
    ! indices(1:5) =           (pre,ozo,sza,vza,raa) plus matching weights(:)
    ids(:) = (/iqu_sza, iqu_vza, iqu_raa/)
    xs(:)  = (/    sza,     vza,     raa/)
    call ndi_find_indices (lut_s % si_dims(ids), xs, indices(3:5))
    call ndi_calc_weights (lut_s % si_dims(ids), xs, indices(3:5), weights(3:5))

    do iw = 1, nswav

      if (lut_s % salbflg (iw) < 1) cycle

      snalbrad(:,:) = 0.0

      do is = 1,nscene
        ip = ndi_find_index (snps(is), pre_grid)
        indices(1) = ip
        weights(1) = (pre_grid(ip+1) - snps(is)) / (pre_grid(ip+1) - pre_grid(ip))

        do ia = 1, nalb
          fidx = oz_info % beg_index
          do iz = 1, oz_info % nzone
            lidx = fidx + oz_info % nzoneo3(iz) - 1
            ioz = ndi_find_index (oz, oz_grid(fidx:lidx))
            ioz = fidx + ioz - 1

            indices(2) = ioz
            weights(2) = (oz_grid(ioz+1) - oz) / (oz_grid(ioz+1) - oz_grid(ioz))

            do k = 1, nterms_5d
              call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
              snalbrad(ia, is) = snalbrad(ia, is) + oz_info % zonefracs(iz) &
                * wt_k * lut_s % si(iw, ia, id_k(1), id_k(2), id_k(3), id_k(4), id_k(5))
            enddo

            fidx = lidx + 1
          enddo
        enddo
      enddo

      alb_grid => lut_s % si_dims (iqu_alb) % x

      ! Derive scene albedo or cloud fraction
      if (use_mler) then
        do is = 1, nscene
          ia = ndi_find_index (snalbs(iw, is), alb_grid)
          wt = (alb_grid(ia+1) - snalbs(iw, is)) / (alb_grid(ia+1) - alb_grid(ia))
          snrad(is) = wt * snalbrad(ia,is) + (1.0 - wt) * snalbrad(ia+1, is)
          snalbwfs(is) = (snalbrad(ia+1,is) - snalbrad(ia,is)) / (alb_grid(ia+1) - alb_grid(ia))
        enddo

        ! Derive effective cloud fraciton using MLER model
        cfrac(iw) = (srad(iw) - snrad(1)) / (snrad(2) - snrad(1))
        cfracwf = snrad(2) - snrad(1)
        if (cfrac(iw) > 1.0) then
          cfrac(iw) = 1.0
          snalbs(iw,2) = snalbs(iw,2) + (srad(iw) - snrad(2)) / snalbwfs(2)
        else if (cfrac(iw) < 0.0) then
          cfrac(iw) = 0.0
          snalbs(iw,1) = snalbs(iw,1) + (srad(iw) - snrad(1)) / snalbwfs(1)
        endif
      else
        fidx = maxloc(snalbrad(1:nalb,1), dim=1, mask=(snalbrad(1:nalb,1) < srad(iw)))
        if (fidx < 1) then
          fidx = 1
        else if (fidx == nalb) then
          fidx = nalb - 1
        endif
        lidx = fidx + 1
        snalbwfs(1) = ((snalbrad(lidx, 1) - snalbrad(fidx, 1)) &
                       / (alb_grid(lidx)- alb_grid(fidx)))
        snalbs(iw,1) = alb_grid(fidx) + (srad(iw) - snalbrad(fidx, 1))/snalbwfs(1)
      endif

      !if (any(snalbs(iw,1:nscene) < 0.0) .or. any(snalbs(iw,1:nscene) > 1.0) ) status = 1
      where(snalbs(iw,1:nscene) < 0.0)
        snalbs(iw,1:nscene) = 0.0
      end where
      where(snalbs(iw,1:nscene) > 1.0)
        snalbs(iw,1:nscene) = 1.0
      end where

    enddo ! iw

  end subroutine pp_derive_albcld

  subroutine pp_get_qu (lut_s, use_mler, swav, sza, vza, raa, oz, oz_info, &
                        nscene, cfrac0, snalbs0, snps, qu_range, out_q, out_u, &
                        errstat)
    implicit none
    type (polpredict_type), target, intent(in) :: lut_s
    logical, intent(in) :: use_mler
    real (kind=r8), dimension(:), intent(in) :: swav
    real (kind=r8), intent(in) :: sza, vza, raa, oz
    type (pp_ozone_zone_info_type), intent(in) :: oz_info
    integer, intent(in) :: nscene
    real (kind=r8), dimension(:), intent(in) :: cfrac0, snps
    real (kind=r8), dimension(:,:), intent(in) :: snalbs0
    type (range_type), intent(in) :: qu_range
    real (kind=r8), dimension(:), intent(inout) :: out_q, out_u
    integer, intent(inout) :: errstat

    integer, parameter :: nd=6, nterms_6d=2**nd
    integer, dimension(nd) :: indices, id_k
    real (kind=r8), dimension(nd) :: weights
    real (kind=r8), pointer, dimension(:) :: wav_grid, pre_grid, alb_grid, oz_grid
    real (kind=r8), allocatable, dimension(:,:) :: snalbs
    real (kind=r8), allocatable, dimension(:) :: cfrac
    integer, dimension(size(swav)) :: albcld_indices
    integer, dimension(3) :: ids
    integer :: nalbcld, ia, is, ip, iz, ioz, iw, nswav, nwav, err, lidx, fidx, k
    real (kind=r8) :: wt_k, frac, q_val, u_val

    if (errstat /= 0) return

    ! All wavelengths with derived cloud fractions/scene albedo
    nswav = size(swav)

    nalbcld = 0
    do iw = 1, nswav
      if (lut_s % salbflg(iw) > 0 ) then
        nalbcld = nalbcld + 1
        albcld_indices(nalbcld) = iw
      endif
    enddo

    wav_grid => lut_s % qu_dims(iqu_wav) % x
    nwav = lut_s % qu_dims(iqu_wav) % dimlen
    allocate (snalbs(nwav,maxscene), &
              cfrac(nwav), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, &
                       "pp_get_qu: malloc failed", errstat)
      return
    endif

    snalbs(:,:) = 0.0
    cfrac(:) = 0.0

    ! Interpolate scene albedo, and optionally cfrac, to output wavelength grid
    fidx = minloc (wav_grid, dim=1, mask=(wav_grid >= swav(albcld_indices(1))))
    lidx = maxloc (wav_grid, dim=1, mask=(wav_grid <= swav(albcld_indices(nalbcld))))
    do is = 1, nscene
      err = spline (swav(albcld_indices(1:nalbcld)), &
                    snalbs0(albcld_indices(1:nalbcld),is), nalbcld, &
                    wav_grid(fidx:lidx), snalbs(fidx:lidx,is), lidx-fidx+1)
      if (err /= 0) then
        call tell_error (tell_runtime_error, &
                         "pp_get_qu: spline failed (snalbs)", errstat)
        return
      endif
      snalbs(1:fidx-1, is) = snalbs(fidx, is)
      snalbs(lidx+1:nwav, is) = snalbs(lidx, is)
    enddo
    if (use_mler) then
      err = spline (swav(albcld_indices(1:nalbcld)), &
                    cfrac0(albcld_indices(1:nalbcld)), nalbcld, &
                    wav_grid(fidx:lidx), cfrac(fidx:lidx), lidx-fidx+1)
      if (err /= 0) then
        call tell_error (tell_runtime_error, &
                         "pp_get_qu: spline failed (cfrac)", errstat)
        return
      endif
      cfrac(1:fidx-1) = cfrac(fidx)
      cfrac(lidx+1:nwav) = cfrac(lidx)
    endif

    ! For each wavelength, we'll interpolate in a 6-D table: (alb,pre,ozo,sza,vza,raa)
    ! The geometric variables (sza,vza,raa) are fixed, while the structural variables
    ! (alb,pre,ozo) vary within the loops over scenes and vertical structure

    ids(:) = (/iqu_sza, iqu_vza, iqu_raa/)
    call ndi_find_indices (lut_s % qu_dims(ids), (/sza, vza, raa/), indices(4:6))
    call ndi_calc_weights (lut_s % qu_dims(ids), (/sza, vza, raa/), indices(4:6), weights(4:6))

    alb_grid => lut_s % qu_dims(iqu_alb) % x
    pre_grid => lut_s % qu_dims(iqu_pre) % x
    oz_grid  => lut_s % qu_dims(iqu_ozo) % x

    out_q(:) = 0.0d0
    out_u(:) = 0.0d0

    ! Interpolate Q, U
    do is = 1, nscene
      if (is == 2 .and. all(cfrac <= 0)) exit ! if no clouds

      ! pressure
      ip = ndi_find_index (snps(is), pre_grid);
      indices(2) = ip
      weights(2) = (pre_grid(ip+1) - snps(is))/(pre_grid(ip+1) - pre_grid(ip));

      fidx = oz_info % beg_index
      do iz = 1, oz_info % nzone
        lidx = fidx + oz_info % nzoneo3(iz) - 1

        if (oz_info % zonefracs(iz) > 0.0) then
          ! ozone
          ioz = ndi_find_index (oz, oz_grid(fidx:lidx))
          ioz = ioz + fidx - 1
          indices(3) = ioz;
          weights(3) = (oz_grid(ioz+1) - oz)/(oz_grid(ioz+1) - oz_grid(ioz));

          do iw = qu_range % imin, qu_range % imax  ! 1, nwav

            ! albedo
            ia = ndi_find_index (snalbs(iw,is), alb_grid)
            indices(1) = ia;
            weights(1) = (alb_grid(ia+1) - snalbs(iw,is)) / (alb_grid(ia+1) - alb_grid(ia));

            if (nscene == 1) then
              frac = 1
            else if (is == 1) then
              frac = 1.0 - cfrac(iw)
            else
              frac = cfrac(iw)
            endif

            q_val = 0.0
            u_val = 0.0
            do k = 1, nterms_6d
              call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
              q_val = q_val + wt_k * lut_s % q(iw, id_k(1),id_k(2),id_k(3),id_k(4),id_k(5),id_k(6))
              u_val = u_val + wt_k * lut_s % u(iw, id_k(1),id_k(2),id_k(3),id_k(4),id_k(5),id_k(6))
            enddo ! k

            out_q(iw) = out_q(iw) + q_val * frac * oz_info % zonefracs(iz)
            out_u(iw) = out_u(iw) + u_val * frac * oz_info % zonefracs(iz)

          enddo ! iw
        endif
        fidx = lidx + 1
      enddo !iz
    enddo ! is

  end subroutine pp_get_qu

  ! Find climatological ozone profiles to be used for the look-up table
  ! based on latitude and total ozone amount
  !
  ! WARNING: assumes size(ozs) = lut_s % qu_dims(iqu_ozo) % dimlen
  !
  subroutine define_o3_zones(ozs, oz, lat, maxzone, &
                             nzone, nzoneo3, ozsidx, ozeidx, zonefracs)
    implicit none
    real*8, dimension(:), intent(in) :: ozs
    real*8, intent(in)                 :: lat, oz

    integer, intent(in)                :: maxzone
    integer, intent(out)               :: nzone, ozsidx, ozeidx
    integer, dimension(maxzone), intent(out) :: nzoneo3
    real*8,  dimension(maxzone), intent(out) :: zonefracs

    ! L: Low latitude profiles (4), M: Mid-latitude profiles (8), H: High-latitude profiles (10)
    REAL*8, PARAMETER :: max_low = 20.0     ! L profile for lat <= max_low
    REAL*8, PARAMETER :: max_lowmid = 30.0  ! L/M mixture for max_low < lat < max_lowmid
    REAL*8, PARAMETER :: max_mid = 50.0     ! M profile for max_lowmid <= lat <= max_mid
    REAL*8, PARAMETER :: max_midhigh = 60.0 ! M/H mixture for max_mid < lat < max_midhigh

    ! To force using 1 ozone zone to speed up the interpolation process,
    ! let max_low = max_lowmid, max_mid == max_midhigh
    !REAL*8, PARAMETER :: max_low = 25.0     ! L profile for lat <= max_low
    !REAL*8, PARAMETER :: max_lowmid = 25.0  ! L/M mixture for max_low < lat < max_lowmid
    !REAL*8, PARAMETER :: max_mid = 55.0     ! M profile for max_lowmid <= lat <= max_mid
    !REAL*8, PARAMETER :: max_midhigh = 55.0 ! M/H mixture for max_mid < lat < max_midhigh

    ! Local variable
    REAL*8 :: abslat

    ! Define zones for ozone profile interpolation
    abslat = ABS(lat)

    zonefracs = 0.0
    IF (abslat <= max_low) THEN
      nzone = 1; nzoneo3(1) = 4
      ozsidx = 1; ozeidx = 4
      zonefracs(1) = 1.0
      IF (oz > ozs(ozeidx)) THEN  ! Switch to mid-latitude
        ozsidx = 5; ozeidx = 12; nzoneo3(1) = 8
      ENDIF
    ELSE IF (abslat < max_lowmid) THEN
      IF (oz <= ozs(4)) THEN
        nzone = 2; nzoneo3(1) = 4; nzoneo3(2) = 8
        ozsidx = 1; ozeidx = 12
        zonefracs(2) = (abslat - max_low) / (max_lowmid - max_low)
        zonefracs(1) = 1.0 - zonefracs(2)
      ELSE
        nzone = 1; nzoneo3(1) = 8
        ozsidx = 5; ozeidx = 12
        zonefracs(1) = 1.0
      ENDIF
    ELSE IF (abslat <= max_mid) THEN
      nzone = 1; nzoneo3(1) = 8
      ozsidx = 5; ozeidx = 12
      zonefracs(1) = 1.0
      IF (oz <= ozs(ozsidx)) THEN
        nzoneo3(1) = 10
        ozsidx = 13; ozeidx = 22
      ENDIF
    ELSE IF (abslat < max_midhigh) THEN
      nzone = 2; nzoneo3(1) = 8; nzoneo3(2) = 10
      ozsidx = 5; ozeidx = 22
      zonefracs(2) = (abslat - max_mid) / (max_midhigh - max_mid)
      zonefracs(1) = 1.0 - zonefracs(2)
      IF (oz <= ozs(ozsidx)) THEN
        nzone = 1; nzoneo3(1) = 10
        ozsidx = 13; ozeidx = 22
        zonefracs(1) = 1.0
      ENDIF
    ELSE
      nzone = 1; nzoneo3(1) = 10
      ozsidx = 13; ozeidx = 22
      zonefracs(1) = 1.0
    ENDIF

    RETURN

  END SUBROUTINE define_o3_zones

  ! Given 2-D tables A, B that depend on x, y and a given point A0, B0
  ! Find x0 and y0 that are associated with A0, B0
  subroutine cldret_find_xy (a, b, x, y, nx, ny, a0, b0, x0, y0, status)

    implicit none
    integer, intent(in)                   :: nx, ny
    real*8, dimension(nx, ny), intent(in) :: a, b
    real*8, dimension(nx), intent(in)     :: x
    real*8, dimension(ny), intent(in)     :: y
    real*8, intent(in)                    :: a0, b0
    real*8, intent(out)                   :: x0, y0
    integer, intent(out)                  :: status ! 1: out of bounds

    ! local variables
    integer :: i, j, xidx, yidx, idx1, idx2, idx
    logical :: found
    real*8  :: tmp1, tmp2, tmp3, tmp4, da, db, dadx, dady, dbdx, dbdy, dx, dy

    found = .false.; status = 0
    i = 1 ! search from the beginning

    ! first find initial search position based on average radiance
    ! to speed up the search
    idx1 = minval(maxloc(a(1:nx, 1), mask=(a(1:nx, 1) < a0)))
    idx2 = minval(maxloc(a(1:nx, ny), mask=(a(1:nx, ny) < a0)))
    i = min(idx1, idx2); if (i < 1) i = 1

    xidx = 0; yidx = 0
    do while (.not. found .and. i <= nx -1)
      j = 1
      do while (.not. found .and. j <= ny -1)
        ! determine if point (a0, b0) is left/right side of a line
        ! <0, left side, >0 right side
        tmp1 = (a(i, j+1) - a(i, j)) * (b0 - b(i, j)) - &
          (b(i, j+1) - b(i, j)) * (a0 - a(i, j))
        tmp2 = (a(i+1, j+1) - a(i+1, j)) * (b0 - b(i+1, j)) - &
          (b(i+1, j+1) - b(i+1, j)) * (a0 - a(i+1, j))

        ! determine if point (a0, b0) is above/below a line
        ! <0, above, >0 below
        tmp3 = -(a(i+1, j) - a(i, j)) * (b0 - b(i, j)) + &
          (b(i+1, j) - b(i, j)) * (a0 - a(i, j))
        tmp4 = -(a(i+1, j+1) - a(i, j+1)) * (b0 - b(i, j+1)) + &
          (b(i+1, j+1) - b(i, j+1)) * (a0 - a(i, j+1))

        if (tmp1 * tmp2 < 0.0 .and. tmp3*tmp4 <= 0.0) then ! brack it, point in a
          xidx = i; yidx = j; found = .true.
        else if (tmp1 < 0 .and. tmp3 < 0) then
          xidx = i - 1; yidx = j - 1; found = .true.
        else if (tmp1 < 0 .and. j == ny - 1 .and. tmp3 > 0.0) then
          xidx = i - 1; yidx = j + 1; found = .true.
        else if (tmp1 < 0 .and. tmp3*tmp4 < 0.0) then
          xidx = i - 1; yidx = j; found = .true.
        else if (tmp1 * tmp2 <= 0.0 .and. j == 1 .and. tmp3 < 0.0) then
          xidx = i; yidx = j - 1; found = .true.
        else if (tmp1 * tmp2 <= 0.0 .and. j == ny-1 .and. tmp3 > 0.0) then
          xidx = i; yidx = j + 1; found = .true.
        else if (i == nx - 1 .and. tmp2 > 0.0 .and. j == 1 .and. tmp3 < 0.0) then
          xidx = i + 1; yidx = j - 1; found = .true.
        else if (i == nx - 1 .and. tmp2 > 0.0 .and. j == ny - 1 .and. tmp3 > 0.0) then
          xidx = i + 1; yidx = j + 1 ; found = .true.
        else if (i == nx - 1 .and. tmp2 > 0.0 .and. tmp3 * tmp4 < 0.0) then
          xidx = i + 1; yidx = j; found = .true.
        endif

        j = j + 1
      enddo
      i = i + 1
    enddo

    if ( xidx >= 1 .and. xidx < nx .and. yidx >= 1 .and. yidx < ny) then ! in bounds
      da = a0 - a(xidx, yidx)
      db = b0 - b(xidx, yidx)
      dx = x(xidx + 1) - x(xidx)
      dy = y(yidx + 1) - y(yidx)

      dadx = (a(xidx+1, yidx) - a(xidx, yidx)) / dx
      dady = (a(xidx, yidx+1) - a(xidx, yidx)) / dy
      dbdx = (b(xidx+1, yidx) - b(xidx, yidx)) / dx
      dbdy = (b(xidx, yidx+1) - b(xidx, yidx)) / dy

      ! special case: [a(xidx, yidx), b(xidx, yidx)] and
      ! [a(xidx, yidx+1), b(xidx, yidx+1)] are the same point
      if (dady == 0.0 .and. dbdy == 0.0) then
        dbdy = -1.0e-4 * b(xidx, yidx) / dy
      endif

      dy = (da * dbdx - db * dadx) / (dady * dbdx - dbdy * dadx)
      dx = (da * dbdy - db * dady) / (dadx * dbdy - dbdx * dady)
      x0 = x(xidx) + dx; y0 = y(yidx) + dy
    else  ! retrievals out of bounds
      if (xidx == 0) then
        x0 = x(1)
      else if (xidx == nx) then
        x0 = x(nx)
      else
        if (yidx == 0) then
          idx = 1
        else
          idx = ny
        endif

        da = a0 - a(xidx, idx)
        dadx = (a(xidx+1, idx) - a(xidx, idx)) / (x(xidx+1) -x(xidx))
        x0 = da / dadx + x(xidx)
      endif

      if (yidx == 0) then
        y0 = y(1)
      else if (yidx == ny) then
        y0 = y(ny)
      else
        if (xidx == 0) then
          idx = 1
        else
          idx = nx
        endif
        db = b0 - b(idx, yidx)
        dbdy = (b(idx, yidx+1) - b(idx, yidx)) / (y(yidx+1) -y(yidx))
        y0 = db / dbdy + y(yidx)
      endif

      status = 1
    endif

    return
  end subroutine cldret_find_xy

end module polpredict_module
