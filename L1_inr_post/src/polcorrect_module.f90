module polcorrect_module
  use, intrinsic :: iso_c_binding
  use netcdf
  use tell_module
  use tio_module
  use types_module
  implicit none
  private

  public :: polcorrect

  type, bind(c), public :: polcorrect_type
    type(c_ptr) :: rad_file    !< radiance file name
    type(c_ptr) :: qu_file     !< QU lookup table file name
    type(c_ptr) :: lps         !< opaque pointer to linear polarization sensitivity struct
    integer (c_int) :: uv_beg, uv_end
    integer (c_int) :: vis_beg, vis_end
    integer (c_int) :: use_mler
    integer (c_int) :: step
    integer (c_int) :: xtrack
    real (kind=c_double) :: delta_pa
  end type

  ! FIXME: these parameters should be provided by tio_module
  integer (c_int), private, parameter :: &
    tempo_band_uv = 0, &
    tempo_band_vis = 1
  character (len=*), private, parameter :: &
    tempo_band_name_uv = 'band_290_490_nm', &
    tempo_band_name_vis = 'band_540_740_nm'

  interface
    function lps_eval (lps, band_index, xtrack, lon, lat, &
                       num_wave, wave, lpsens, angmax) &
        bind (c, name='lps_eval')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int
      implicit none
      type (c_ptr), value :: lps
      integer (c_int), value :: band_index, xtrack, num_wave
      real (kind=c_double), value :: lon, lat
      real (kind=c_double), dimension(num_wave), intent(in) :: wave
      real (kind=c_double), dimension(num_wave), intent(out) :: lpsens
      real (kind=c_double), dimension(num_wave), intent(out) :: angmax
      integer (c_int) :: lps_eval
    end function lps_eval
  end interface

  interface
    function tio_get_fill_value (grp, name, nofill, fillvalue) &
        bind (c, name='TIO_get_fill_value')
      use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_float, c_char
      implicit none
      integer(c_int), intent(in) :: grp
      character (len=1, kind=c_char), intent(in) :: name
      integer(c_int), intent(out) :: nofill
      real (c_float), intent(out) :: fillvalue
      integer(c_int) :: tio_get_fill_value
    end function tio_get_fill_value
  end interface

  interface
    function tio_time_tempo_to_utc_caldate (tempo_time,year,month,day,hour) &
      bind (c, name='tio_time_tempo_to_utc_caldate')
      use, intrinsic :: iso_c_binding, only : c_double, c_int
      implicit none
      real (c_double), value :: tempo_time
      integer (c_int), intent(out) :: year, month, day
      real (c_double), intent(out) :: hour
      integer (c_int) :: tio_time_tempo_to_utc_caldate
    end function tio_time_tempo_to_utc_caldate
  end interface

  real (kind=r8), parameter :: r8_fill = nf90_fill_double
  real (kind=r8), parameter :: cldalb0 = 0.8d0  ! cloud albedo

  ! radiance_status bit for polarization correction status
  integer, parameter :: polcorr_status_bit = 0

contains

  integer(c_int) function polcorrect (pt) bind (C, name='polcorrect')
    use polpredict_module, only: pp_dealloc, polpredict_type
    implicit none
    type (polcorrect_type), intent(in) :: pt
    type (polpredict_type) :: lut_s
    type (radiance_type) :: rad_s
    character(:,kind=c_char), pointer :: rad_file
    integer :: errstat

    errstat = 0

    rad_file => c_f_string (pt % rad_file)

    call tiof_open (rad_file, rad_s % obj, nf90_write, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, 'reading radiance file: '//rad_file, errstat)
      polcorrect = errstat
      return
    endif

    call process_group (pt, lut_s, rad_s, tempo_band_uv, errstat)
    call process_group (pt, lut_s, rad_s, tempo_band_vis, errstat)

    call tiof_close (rad_s % obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, 'error closing file:'//rad_file, errstat)
      polcorrect = errstat
      return
    endif

    call pp_dealloc (lut_s, errstat)

    polcorrect = 0
  end function

  subroutine check_pol_correction_status (obj, has_pol_correction, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    logical, intent(out) :: has_pol_correction
    integer, intent(inout) :: errstat
    integer :: status_flag
    if (errstat /= 0) return

    has_pol_correction = .false.

    call tiof_get_i4 (obj, tempo_var_radiance_status, status_flag, errstat)
    if (errstat /= 0) return

    has_pol_correction = btest (status_flag, polcorr_status_bit)
  end subroutine check_pol_correction_status

  subroutine set_polcorr_status_bit (obj, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    integer, intent(inout) :: errstat
    integer :: status_flag
    if (errstat /= 0) return

    call tiof_get_i4 (obj, tempo_var_radiance_status, status_flag, errstat)
    if (errstat /= 0) return

    status_flag = ibset (status_flag, polcorr_status_bit)

    call tiof_put_i4 (obj, tempo_var_radiance_status, status_flag, errstat)
    if (errstat /= 0) return

  end subroutine set_polcorr_status_bit

  subroutine process_group (pt, lut_s, rad_s, band_id, errstat)
    use polpredict_module
    implicit none
    type(polcorrect_type), intent(in) :: pt
    type(polpredict_type), intent(inout) :: lut_s
    type(radiance_type), intent(inout) :: rad_s
    integer(c_int), intent(in) :: band_id
    integer, intent(inout) :: errstat

    ! FIXME: Will these booleans always be true in the operational code?
    logical, parameter :: retoz = .true., retctp = .true.
    integer, parameter :: num_ctp_waves = 3, num_oz_waves = 3

    type (radiance_subset_type) :: subset
    type (pp_ozone_zone_info_type) :: oz_info
    integer :: step, ix, err, beg_step, end_step, beg_xtrack, end_xtrack
    integer :: retoz_errstat, retctp_errstat
    integer :: iter, num_iter, nscene, num_swav, num_wav, cld_wave_index
    integer, dimension(2) :: swav_limits
    integer, dimension(num_ctp_waves) :: ctp_wave_indices
    integer, dimension(num_oz_waves) :: oz_wave_indices
    real (kind=r4) :: fmonth
    real (kind=r8) :: lon, lat, pre, sza, vza, raa
    real (kind=r8) :: oz, oz0, oz_saved
    real (kind=r8) :: ctp, ctp0, ctp_saved
    real (kind=r8), allocatable, dimension(:) :: srad, swav, wav, q, u, lpserr
    real (kind=r8), allocatable, dimension(:) :: cfracs0, tmp_cfracs0
    real (kind=r8), allocatable, dimension(:,:) :: snalbs0, tmp_snalbs0
    real (kind=r8), dimension(maxscene) :: snps
    type (range_type) :: qu_range

    logical :: use_mler

    if (errstat /= 0) return

    use_mler = (pt % use_mler /= 0)

    call read_radiance_fmonth (rad_s, fmonth, errstat)
    if (errstat /= 0) return

    ! Wavelength intervals and array indices to be used for processing this band
    select case (band_id)
      case (tempo_band_uv)
        swav_limits = (/pt % uv_beg, pt % uv_end/)
        cld_wave_index = 5   ! 367 nm
        ctp_wave_indices = (/10, 11, 12/)
        oz_wave_indices = (/1, 2, 5/)
        qu_range % min = 285.0
        qu_range % max = 495.0
      case (tempo_band_vis)
        swav_limits = (/pt % vis_beg, pt % vis_end/)
        cld_wave_index = 7   ! 670 nm
        ctp_wave_indices = (/8, 9, 10/)
        oz_wave_indices = (/3, 4, 6/)
        qu_range % min = 535.0
        qu_range % max = 745.0
    end select

    call read_radiance_geometry (rad_s, band_id, errstat)
    if (errstat /= 0) return

    call calc_relative_azimuth_angle (rad_s, errstat)
    call calc_surface_pressure (rad_s, errstat)
    call define_subset (rad_s, subset, errstat)
    if (errstat /= 0) return

    call init_pp_lut (pt % qu_file, fmonth, subset, lut_s, errstat)
    call pp_get_swav (lut_s, swav_limits(1), swav_limits(2), swav, errstat)
    if (errstat /= 0) return

    num_swav = size(swav)

    call pp_get_wav (lut_s, wav, errstat)
    if (errstat /= 0) return

    num_wav = size(wav)

    ! Computing Q, U is costly, so minimize the number of wavelengths computed:
    qu_range % imin = ndi_find_index (qu_range % min, wav)
    qu_range % imax = ndi_find_index (qu_range % max, wav) + 1

    allocate (snalbs0(num_swav,maxscene), tmp_snalbs0(num_swav,maxscene), &
              cfracs0(num_swav), tmp_cfracs0(num_swav), &
              srad(num_swav), &
              q(num_wav), &
              u(num_wav), &
              lpserr(rad_s % num_wave), &
              stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    if (pt % step == 0) then
      beg_step = 1
      end_step = rad_s % num_step
    else
      beg_step = pt % step
      end_step = pt % step
    endif

    if (pt % xtrack == 0) then
      beg_xtrack = 1
      end_xtrack = rad_s % num_xtrack
    else
      beg_xtrack = pt % xtrack
      end_xtrack = pt % xtrack
    endif

    snalbs0 = 0.0
    cfracs0 = 0.0

    retoz_errstat = 0
    retctp_errstat = 0
    oz_saved = 0.0
    ctp_saved = 0.0

    do step = beg_step, end_step

      call read_radiance_for_mirror_step (rad_s, step, errstat)
      if (errstat /= 0) return

      do ix = beg_xtrack, end_xtrack

        lat = rad_s % lat(ix,step)
        lon = rad_s % lon(ix,step)
        pre = rad_s % pre(ix,step)

        ! get ozone for pre
        call pp_interp_ozone_pre (lut_s, lat, lon, pre, oz0, errstat)
        call pp_interp_ctp (lut_s, lat, lon, ctp0, errstat)
        if (errstat /= 0) cycle

        if (use_mler) then
          nscene = 2
          snps(:) = (/pre, ctp0/)
          call pp_interp_surface_albedo (lut_s, lat, lon, swav, snalbs0(:,1), errstat)
          snalbs0(:,2) = cldalb0
        else
          nscene = 1
          snps(1) = pre
          snalbs0(:,1) = 0.05d0
        endif

        ! If possible, use oz, ctp from previous iteration
        oz = oz0
        ctp = ctp0
        if (ix > beg_xtrack) then
          if (retoz .and. retoz_errstat == 0) oz = oz_saved
          if (retctp .and. retctp_errstat == 0) ctp = ctp_saved
        endif
        retoz_errstat = 0
        retctp_errstat = 0

        ! Define ozone profile down to surface pressure
        call pp_interp_ozone_profile (lut_s, pre, errstat)
        if (errstat /= 0) return

        ! Define o3 zones
        call pp_ozone_zone_info (lut_s, lat, oz, oz_info, errstat)
        if (errstat /= 0) return

        if (.not. (use_mler.or.retctp)) then
          ! If not retcp, need to derive scene pressure here, otherwise
          ! scene pressure will be derived during cloud retrieval
          nscene = 2
          snps(1) = pre
          snps(2) = ctp0
          snalbs0(cld_wave_index, 1) = 0.05
          snalbs0(cld_wave_index, 2) = cldalb0
          call pp_derive_scene_pressure (lut_s, cld_wave_index, srad, sza, vza, raa, oz, &
                                         oz_info, snalbs0, snps, ctp, errstat)
          if (errstat /= 0) ctp = pre
          nscene = 1
          snps(1) = ctp
        endif

        ! Ozone (oz) and cloud-top pressure (ctp) depend on each other.
        ! Iterate unless a previously successful retrieval is available.
        num_iter = 1
        if (retctp .and. retoz) then
          num_iter = 2
          if (ix /= beg_xtrack .and. retoz_errstat == 0 .and. retctp_errstat == 0) num_iter = 1
        endif

        err = spline (rad_s % wave(:,ix), rad_s % radiance (:,ix), rad_s % num_wave, &
                      swav, srad, num_swav)
        if (err /= 0) then
          call tell_error (tell_runtime_error, "spline interpolation failed", errstat)
          return
        endif

        sza = rad_s % sza(ix, step)
        vza = rad_s % vza(ix, step)
        raa = rad_s % raa(ix, step)

        do iter=1,num_iter

          if (retctp) then
            ! Derive cloud-top pressure
            tmp_snalbs0 = snalbs0(ctp_wave_indices,:)
            tmp_cfracs0 = cfracs0(ctp_wave_indices)
            call pp_derive_ctp (lut_s, use_mler, srad(ctp_wave_indices), &
                                sza, vza, raa, oz, oz_info, ctp_wave_indices, &
                                nscene, ctp, tmp_cfracs0, tmp_snalbs0, snps, &
                                retctp_errstat, errstat)
            snalbs0(ctp_wave_indices,:) = tmp_snalbs0
            cfracs0(ctp_wave_indices) = tmp_cfracs0
            if (errstat /= 0) return

            if (isnan(ctp)) retctp_errstat = 4
            if (use_mler .and. ctp > pre) then
              retctp_errstat = 2
              ctp = pre
            endif
            if (retctp_errstat == 0) then
              ctp_saved = ctp
            else
              ctp = ctp0
            endif
          endif

          if (retoz) then
            ! Derive total ozone
            tmp_snalbs0 = snalbs0(oz_wave_indices,:)
            tmp_cfracs0 = cfracs0(oz_wave_indices)
            call pp_derive_to3 (lut_s, use_mler, &
                                swav(oz_wave_indices), srad(oz_wave_indices), &
                                sza, vza, raa, oz, oz_info, oz_wave_indices, &
                                nscene, tmp_cfracs0, tmp_snalbs0, snps, &
                                retoz_errstat, errstat)
            snalbs0(oz_wave_indices,:) = tmp_snalbs0
            cfracs0(oz_wave_indices) = tmp_cfracs0
            if (errstat /= 0) return

            if (isnan(oz)) retoz_errstat = 4
            if (retoz_errstat == 0) then
              oz_saved = oz
            else
              oz = oz0
            endif

            ! Define o3 zones
            call pp_ozone_zone_info (lut_s, lat, oz, oz_info, errstat)
            if (errstat /= 0) return
          endif

        enddo ! iter

        call pp_derive_albcld (lut_s, use_mler, swav, srad, sza, vza, raa, oz, &
                               oz_info, nscene, cfracs0, snalbs0, snps, errstat)
        if (errstat /= 0) return

        call pp_get_qu (lut_s, use_mler, swav, sza, vza, raa, oz, oz_info, &
                        nscene, cfracs0, snalbs0, snps, qu_range, q, u, errstat)
        if (errstat /= 0) return

        call calc_lpserr (pt % lps, pt % delta_pa, wav, q, u, &
                          rad_s, band_id, ix, step, lpserr, errstat)
        if (errstat /= 0) return

        where (rad_s % radiance (:,ix) /= r8_fill)
          rad_s % radiance(:, ix) = rad_s % radiance(:,ix) / (1.0 + lpserr(:))
        end where

      enddo ! ix
      call write_radiance_for_mirror_step (rad_s, step, errstat)
      if (errstat /= 0) return
    enddo ! step

    call dealloc_radiance (rad_s, errstat)

  end subroutine process_group

  subroutine calc_lpserr (lps, delta_pa, wav, q0, u0, rad_s, &
                          band_id, ix, step, lpserr, errstat)
    use, intrinsic :: ieee_arithmetic
    implicit none
    type (c_ptr), value :: lps
    real (kind=r8), intent(in) :: delta_pa
    real (kind=r8), dimension(:), intent(in) :: wav, q0, u0
    type (radiance_type), target, intent(in) :: rad_s
    integer, intent(in) :: band_id, ix, step
    real (kind=r8), dimension(:), intent(inout) :: lpserr
    integer, intent(inout) :: errstat

    real (kind=r8), parameter :: degtorad = 4.0_r8*atan(1.0_r8)/180.0_r8
    real (kind=r8), allocatable, dimension(:) :: q, u, dolps, lpsens, angmax
    real (kind=r8), target, allocatable, dimension(:) :: pa
    real (kind=r8), pointer, dimension(:) :: rad_wave
    real (kind=r8) :: lon, lat
    integer :: err, nwav, num_rad_wave

    if (errstat /= 0) return

    if (rad_s % sza(ix,step) == r8_fill) then
      lpserr(:) = 0.0
      return
    endif

    nwav = size(wav)
    num_rad_wave = rad_s % num_wave
    rad_wave => rad_s % wave (:, ix)

    allocate (q(num_rad_wave), u(num_rad_wave), &
              dolps(num_rad_wave), pa(num_rad_wave), &
              lpsens(num_rad_wave), angmax(num_rad_wave), stat=err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "calc_lpserr: malloc failed", errstat)
      return
    endif

    ! Map Q, U from LUT wav grid to radiance measurement wavelength grid:

    err = spline (wav, q0, nwav, rad_wave, q, num_rad_wave)
    if (err /= 0) then
      call tell_error (tell_runtime_error, "calc_lpserr: Q interpolation failed", errstat)
      return
    endif

    err = spline (wav, u0, nwav, rad_wave, u, num_rad_wave)
    if (err /= 0) then
      call tell_error (tell_runtime_error, "calc_lpserr: U interpolation failed", errstat)
      return
    endif

    !do i=1,num_rad_wave
    !  write(*,'(3f10.4)')rad_wave(i), q(i), u(i)
    !enddo

    ! dolps = degree of linear polarization
    dolps = sqrt (q*q + u*u)

    ! pa = angle between the plane of linear polarization,
    !      and the LUT plane of reference, 0 <= pa <= 180
    pa = 0.5 * atan2(u, q) / degtorad
    where (pa < 0.0)
      pa = pa + 180.0
    end where

    ! At this point, the plane of reference for pa is the local meridian
    ! plane, LMP. However, the instrument reference plane, IRP, is rotated
    ! relative to the LMP by an angle, delta_pa. To apply the instrument
    ! linear polarization sensitivity tables, we need the angle between
    ! the plane of linear polarization and the IRP, which is:

    pa = pa + delta_pa

    ! Use (lon,lat) to derive this spectrum's angular offset from the
    ! instrument boresight, then perform a table lookup to obtain the
    ! associated values of:
    !     lpsens = linear polarization sensitivity
    !     angmax = angle of maximum transmission

    lon = rad_s % lon (ix, step)
    lat = rad_s % lat (ix, step)

    err = lps_eval (lps, band_id, ix, lon, lat, num_rad_wave, rad_wave, lpsens, angmax)
    if (err /= 0) then
      call tell_set_error (errstat)
      return
    endif

    ! Measured radiance, I', true radiance, I:
    ! I' = I * ( 1 + 2.0 * lps * Dolp * cos(2.0 * (pa - maxang))
    ! I = I' / ( 1 + 2.0 * lps * Dolp * cos(2.0 * (pa - maxang))
    ! where pa: phase angle of polarization wrt to instrument reference plane (irp),
    ! and maxang: angle of maximum transmission

    lpserr = 2.0 * lpsens * dolps * cos(2.0 * (pa - angmax) * degtorad)

    where (.not.ieee_is_finite(lpserr))
      lpserr = 0.0
    end where

  end subroutine calc_lpserr

  subroutine init_pp_lut (fileptr, fmonth, subset, lut_s, errstat)
    use polpredict_module, only : pp_read, pp_initialized, polpredict_type
    implicit none
    type (c_ptr), value :: fileptr
    real (kind=r4), intent(in) :: fmonth
    type (radiance_subset_type), intent(inout) :: subset
    type (polpredict_type), intent(inout) ::lut_s
    integer, intent(inout) :: errstat

    character (:,kind=c_char), pointer :: lut_file
    type (tiof_file_type) :: lut_obj

    if (errstat /= 0) return
    if (pp_initialized (lut_s)) return

    lut_file => c_f_string (fileptr)

    call tiof_open (lut_file, lut_obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, 'reading lookup table: '//lut_file, errstat)
      return
    endif

    ! subset lookup tables on input
    call pp_read (lut_obj, subset, fmonth, lut_s, errstat)
    call tiof_close (lut_obj, errstat)

  end subroutine init_pp_lut

  subroutine dealloc_radiance (rad_s, errstat)
    implicit none
    type(radiance_type), intent(inout) :: rad_s
    integer, intent(inout) :: errstat
    integer :: err

    deallocate (rad_s % radiance, &
                rad_s % wave, &
                rad_s % lon, rad_s % lat, &
                rad_s % sza, rad_s % saa, &
                rad_s % vza, rad_s % vaa, &
                rad_s % raa, &
                rad_s % hgt, rad_s % pre, &
                stat = err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "dealloc_radiance: dealloc failed", errstat)
      return
    endif
  end subroutine dealloc_radiance

  subroutine alloc_radiance (rad_s, num_step, num_xtrack, num_wave, errstat)
    implicit none
    type(radiance_type), intent(inout) :: rad_s
    integer, intent(in) :: num_step, num_xtrack, num_wave
    integer, intent(inout) :: errstat
    integer :: err

    rad_s % num_step = num_step
    rad_s % num_xtrack = num_xtrack
    rad_s % num_wave = num_wave

    allocate (rad_s % radiance(num_wave, num_xtrack), &
              rad_s % wave(num_wave, num_xtrack), &
              rad_s % lon(num_xtrack, num_step), &
              rad_s % lat(num_xtrack, num_step), &
              rad_s % sza(num_xtrack, num_step), &
              rad_s % saa(num_xtrack, num_step), &
              rad_s % vza(num_xtrack, num_step), &
              rad_s % vaa(num_xtrack, num_step), &
              rad_s % raa(num_xtrack, num_step), &
              rad_s % hgt(num_xtrack, num_step), &
              rad_s % pre(num_xtrack, num_step), &
              stat = err)
    if (err /= 0) then
      call tell_error (tell_malloc_error, "alloc_radiance: malloc failed", errstat)
      return
    endif

  end subroutine alloc_radiance

  subroutine calc_relative_azimuth_angle (rad_s, errstat)
    implicit none
    type(radiance_type), intent(inout) :: rad_s
    integer, intent(inout) :: errstat
    integer :: i, j, dimlens(2)
    real (kind=r8) :: saa, vaa, raa

    if (errstat /= 0) return

    dimlens = shape(rad_s % vaa)

    do j = 1, dimlens(2)
      do i = 1, dimlens(1)
        vaa = rad_s % vaa(i,j)
        saa = rad_s % saa(i,j)

        if (vaa == r8_fill .or. saa == r8_fill) then
          raa = r8_fill
        else
          raa = 180.d0 + (saa - vaa)
          if (raa > 180.d0) then
            raa = raa - 360.d0
          else if (raa < -180.d0) then
            raa = raa + 360.d0
          endif
        endif

        rad_s % raa(i,j) = raa
      enddo
    enddo

  end subroutine calc_relative_azimuth_angle

  subroutine calc_surface_pressure (rad_s, errstat)
    implicit none
    real (kind=r8), parameter :: &
      pressure0_hpa = 1013.25, &
      pressure_scale_height_meters = 16.0d3
    type(radiance_type), intent(inout) :: rad_s
    integer, intent(inout) :: errstat
    integer :: i, j, dimlens(2)

    real (kind=r8) :: z_meters, pre

    if (errstat /= 0) return

    dimlens = shape(rad_s % hgt)

    do j = 1, dimlens(2)
      do i = 1, dimlens(1)

        ! terrain height [meters]
        z_meters = rad_s % hgt(i,j)

        if (z_meters == r8_fill) then
          pre = r8_fill
        else
          ! surface pressure [hPa]
          pre = pressure0_hpa * 10**(-z_meters/pressure_scale_height_meters)
        endif

        rad_s % pre(i,j) = pre

      enddo
    enddo

  end subroutine calc_surface_pressure

  subroutine read_radiance_fmonth (rad_s, fmonth, errstat)
    implicit none
    type (radiance_type), intent(in) :: rad_s
    real (kind=r4), intent(out) :: fmonth
    integer, intent(inout) :: errstat

    integer :: err, year, month, day
    real (kind=r8) :: tstart, hour

    if (errstat /= 0) return

    err = nf90_get_att (rad_s % obj % fileid, nf90_global, &
                        "time_coverage_start_since_epoch", tstart)
    if (err /= 0) then
      call tell_error (tell_io_read_error, &
                       "process_group: error reading granule start time", &
                       errstat)
      return
    endif

    err = tio_time_tempo_to_utc_caldate (tstart, year, month, day, hour)
    if (err /= 0) then
      call tell_error (tell_io_read_error, &
                       "process_group: error processing granule start time", &
                       errstat)
      return
    endif

    ! The monthly tables give values at mid-month, so fmonth=1.0 means Jan 15.
    fmonth = (month-0.5) + day/30.0
    if (fmonth > 12.0) fmonth = fmonth - 12.0

  end subroutine read_radiance_fmonth

  subroutine read_radiance_geometry (rad_s, band_id, errstat)
    implicit none
    type(radiance_type), intent(inout) :: rad_s
    integer(c_int), intent(in) :: band_id
    integer, intent(inout) :: errstat

    integer :: num_step, num_xtrack, num_wave
    integer, dimension(2) :: start, edge
    logical :: has_pol_correction

    if (errstat /= 0) return

    call tiof_inq_group (rad_s % obj, "/", errstat)
    call tiof_inq_dimlen (rad_s % obj, tempo_dim_step, num_step, errstat)

    select case (band_id)
      case (tempo_band_uv)
        call tiof_inq_group (rad_s % obj, tempo_band_name_uv, errstat)
      case (tempo_band_vis)
        call tiof_inq_group (rad_s % obj, tempo_band_name_vis, errstat)
    end select

    call check_pol_correction_status (rad_s % obj, has_pol_correction, errstat)
    if (errstat /= 0) return
    if (has_pol_correction) then
      call tell_error (tell_runtime_error, &
                       "polarization correction has already been applied", errstat)
      return
    endif

    call tiof_inq_dimlen (rad_s % obj, tempo_dim_xtrack, num_xtrack, errstat)
    call tiof_inq_dimlen (rad_s % obj, tempo_dim_channel, num_wave, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "reading granule dimensions", errstat)
      return
    endif

    call alloc_radiance (rad_s, num_step, num_xtrack, num_wave, errstat)
    if (errstat /= 0) return

    start = (/0,0/)
    edge = (/num_step, num_xtrack/)

    call tiof_get2d_r8 (rad_s % obj, tempo_var_longitude, start, edge, &
                        rad_s % lon, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_latitude, start, edge, &
                        rad_s % lat, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_sz_angle, start, edge, &
                        rad_s % sza, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_sa_angle, start, edge, &
                        rad_s % saa, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_vz_angle, start, edge, &
                        rad_s % vza, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_va_angle, start, edge, &
                        rad_s % vaa, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_terr_height, start, edge, &
                        rad_s % hgt, errstat, replace_fill=r8_fill)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "reading geometry variables", errstat)
      return
    endif

  end subroutine read_radiance_geometry

  subroutine read_radiance_for_mirror_step (rad_s, step, errstat)
    implicit none
    type(radiance_type), intent(inout) :: rad_s
    integer, intent(in) :: step
    integer, intent(inout) :: errstat
    integer, dimension(3) :: start, edge

    character (len=128) :: msg

    if (errstat /= 0) return

    start = (/step-1, 0,0/)
    edge  = (/1, rad_s % num_xtrack, rad_s % num_wave/)

    rad_s % this_step = step

    call tiof_get2d_r8 (rad_s % obj, tempo_var_radiance, start, edge, &
                        rad_s % radiance, errstat, replace_fill=r8_fill)
    call tiof_get2d_r8 (rad_s % obj, tempo_var_wavelength, start, edge, &
                        rad_s % wave, errstat, replace_fill=r8_fill)
    if (errstat /= 0) then
      write(msg,'(a,i0)')'reading radiances: step=',step
      call tell_error (tell_io_read_error, msg, errstat)
      return
    endif

  end subroutine read_radiance_for_mirror_step

  subroutine write_radiance_for_mirror_step (rad_s, step, errstat)
    use, intrinsic :: iso_c_binding, only : c_null_char
    implicit none
    type(radiance_type), intent(inout) :: rad_s
    integer, intent(in) :: step
    integer, intent(inout) :: errstat
    integer, dimension(3) :: start, edge

    character (len=128) :: msg
    integer :: err, nofill, varid
    real (kind=r4) :: fill_value

    if (errstat /= 0) return

    start = (/step-1, 0,0/)
    edge  = (/1, rad_s % num_xtrack, rad_s % num_wave/)

    ! FIXME - a libtio interface is preferable
    err = nf90_inq_varid (rad_s % obj % groupid, tempo_var_radiance, varid)
    err = nf90_inq_var_fill (rad_s % obj % groupid, varid, nofill, fill_value)
    if (err == 0) then
      where (rad_s % radiance == r8_fill)
        rad_s % radiance = real(fill_value, kind=r8)
      end where
    endif

    call tiof_put2d_r8 (rad_s % obj, tempo_var_radiance, start, edge, &
                        rad_s % radiance, errstat)
    if (errstat /= 0) then
      write(msg,'(a,i0)')'writing radiances: step=',step
      call tell_error (tell_io_write_error, msg, errstat)
      return
    endif

    call set_polcorr_status_bit (rad_s % obj, errstat)

  end subroutine write_radiance_for_mirror_step

  subroutine define_subset (rad_s, subset, errstat)
    implicit none
    type(radiance_type), intent(in) :: rad_s
    type(radiance_subset_type), intent(out) :: subset
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    subset % lon % min = minval(rad_s % lon, rad_s % lon /= r8_fill)
    subset % lon % max = maxval(rad_s % lon, rad_s % lon /= r8_fill)

    subset % lat % min = minval(rad_s % lat, rad_s % lat /= r8_fill)
    subset % lat % max = maxval(rad_s % lat, rad_s % lat /= r8_fill)

    subset % sza % min = minval(rad_s % sza, rad_s % sza /= r8_fill)
    subset % sza % max = maxval(rad_s % sza, rad_s % sza /= r8_fill)

    subset % vza % min = minval(rad_s % vza, rad_s % vza /= r8_fill)
    subset % vza % max = maxval(rad_s % vza, rad_s % vza /= r8_fill)

    subset % raa % min = minval(rad_s % raa, rad_s % raa /= r8_fill)
    subset % raa % max = maxval(rad_s % raa, rad_s % raa /= r8_fill)

    ! This is used to subset the lookup table (LUT) on input.
    ! Since we only want to read the LUT *once*, the subset
    ! wavelength range should cover *both* TEMPO wavelength bands.

    subset % wav % min = lut_wav_min
    subset % wav % max = lut_wav_max

  end subroutine define_subset

  function c_f_string(c_str) result(f_str)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char
    type(c_ptr), intent(in) :: c_str
    character(:,kind=c_char), pointer :: f_str
    character(kind=c_char), pointer :: arr(:)
    interface
      function strlen(s) bind (c, name='strlen')
        use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
        implicit none
        type(c_ptr), intent(in), value :: s
        integer (c_size_t) :: strlen
      end function strlen
    end interface
    call c_f_pointer (c_str, arr, [strlen(c_str)])
    call get_scalar_pointer (size(arr), arr, f_str)
  end function

  subroutine get_scalar_pointer (scalar_len, scalar, ptr)
    use, intrinsic :: iso_c_binding, only : c_char
    integer, intent(in) :: scalar_len
    character (kind=c_char, len=scalar_len), intent(in), target :: scalar(1)
    character (:, kind=c_char), intent(out), pointer :: ptr
    ptr => scalar(1)
  end subroutine get_scalar_pointer

end module polcorrect_module
