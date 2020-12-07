!> Reader routines for use with GEMS data
module m_read_gems
  use tell_module
  use tio_module
  use tg_names_module
  use netcdf, only: nf90_nowrite

  implicit none

  private
  public gems_read_l1_rad_info, gems_read_irrad_data, gems_read_latitude, &
       gems_read_ice_glint, gems_read_radiance_lines, gems_read_geofields, &
       gems_read_cld, gems_read_earth_sun_distance

contains

  !> Read dimension info from GEMS L1C radiance file
  !------------------------------------------------------------------------
  !
  ! @param[in]  l1bfile   L1 radiance file name
  ! @param[out] rpt       radiance parameter type strcuture
  ! @param      errstat   error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_l1_rad_info (l1bfile, rpt, errstat)

    USE OMSAO_variables_module,  ONLY : Radiance_Paras_Type

    implicit none
    character (len=*), intent(in) :: l1bfile
    type (Radiance_Paras_Type), intent(out) :: rpt
    integer (kind=4), intent(inout) :: errstat

    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    rpt%ntimes = 0 ; rpt%nxtrack = 0 ; rpt%nwavel_ccd = 0
    rpt%l1bfilename = l1bfile
    rpt%l1bchannel = "UV2"

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    !call tiof_use_file_epoch (tio_l1obj, errstat)
    call tiof_time_set_taix_epoch("2000-01-01T00:00:00Z",errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_x", rpt%ntimes, errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_y", rpt%nxtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_band", rpt%nwavel_ccd, errstat)
    if (errstat /= 0) return

    call tiof_close (tio_l1obj, errstat)

  end subroutine gems_read_l1_rad_info


  !> Read data from GEMS irradiance file
  !------------------------------------------------------------------------
  !
  ! @param[out] nwavel      irradiance wavelength dimension size
  ! @param[out] nxtrack     irradiance xtrack dimension size
  ! @param[out] wavelengths 2D irradiance wavelength array
  ! @param[out] spectrum    2D irradiance value array
  ! @param[out] qflags      2D irradiance pixel quality flag array
  ! @param      errstat     error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_irrad_data (nwavel, nxtrack, wavelengths, spectrum, &
       qflags, errstat)

    use OMSAO_variables_module, only: l1b_irrad_filename

    implicit none

    !output variables
    integer (kind=4), intent(out) :: nwavel, nxtrack
    real (kind=8), dimension(:,:), allocatable, intent(out) :: wavelengths, &
         spectrum
    integer (kind=2), dimension(:,:), allocatable, intent(out) :: qflags
    integer (kind=4), intent(inout) :: errstat
    integer, dimension(1:2), parameter :: flip = (/2,1/)

    !local variables
    integer (kind=2), dimension(:,:), allocatable :: tmp_qflags
    real (kind=4), dimension(:,:), allocatable :: tmp_wavelengths, tmp_spectrum
    type(tiof_file_type) :: tio_irrobj
    integer :: staterr

    if (errstat /= 0) return

    call tell_log (1, 'reading GEMS irradiances = '//trim(l1b_irrad_filename))

    call tiof_open(l1b_irrad_filename, tio_irrobj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (tio_irrobj, "dim_image_y", nxtrack, errstat)
    call tiof_inq_dimlen (tio_irrobj, "dim_image_band", nwavel, errstat)

    allocate (wavelengths (nwavel, nxtrack), &
         spectrum (nwavel, nxtrack), &
         qflags (nwavel, nxtrack), &
         tmp_wavelengths (nxtrack, nwavel), &
         tmp_spectrum (nxtrack, nwavel), &
         tmp_qflags (nxtrack, nwavel), stat=staterr)
    if (staterr /= 0) then
      call tell_error (tell_malloc_error, &
           "gems_read_irrad_data: allocation failed", errstat)
      return
    endif

    call tiof_get2d_r4 (tio_irrobj, "wavelength", [0,0], &
         [nwavel, nxtrack], tmp_wavelengths(1:nxtrack,:), errstat)
    call tiof_get2d_r4 (tio_irrobj, "image_pixel_values", [0,0], &
         [nwavel, nxtrack], tmp_spectrum(1:nxtrack,:), errstat)
    call tiof_get2d_i2 (tio_irrobj, "bad_pixel_mask", [0,0], &
         [nwavel, nxtrack], tmp_qflags(1:nxtrack,:), errstat)
    call tiof_close (tio_irrobj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "gems_read_irrad failed", &
           errstat)
      return
    endif

    wavelengths = real (reshape (tmp_wavelengths,(/nwavel,nxtrack/), &
         order=flip), kind=8)
    spectrum = real (reshape (tmp_spectrum,(/nwavel,nxtrack/),order=flip), &
         kind=8)
    qflags = reshape (tmp_qflags,(/nwavel,nxtrack/),order=flip)

    deallocate (tmp_wavelengths, tmp_spectrum, tmp_qflags, stat=errstat)

  end subroutine gems_read_irrad_data



  !> Read subset of GEMS latitude array
  !------------------------------------------------------------------------
  !
  ! @param[in]  l1bfile   L1 radiance file name
  ! @param[in]  tstart    starting along-track step
  ! @param[in}  ntimes    number of steps to read
  ! @param[out] latitude  2D latitude array
  ! @param      errstat   error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_latitude (l1bfile, tstart, ntimes, latitude, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile
    integer (kind=4), intent(in) :: tstart, ntimes
    !output variables
    real (kind=4), dimension(:,:), intent(out) :: latitude
    integer (kind=4), intent(inout) :: errstat
    !local variables
    real (kind=4), dimension(:,:), allocatable :: tmp_lat
    type (tiof_file_type) :: tio_l1obj
    integer :: nxtrack, lxtrack, lstep
    integer, dimension(1:2), parameter :: flip = (/2,1/)

    if (errstat /= 0) return

    lxtrack=size(latitude,1)
    lstep=size(latitude,2)

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (tio_l1obj, "dim_image_y", nxtrack, errstat)
    if (errstat /= 0 .or. nxtrack /= lxtrack) then
      call tell_error (tell_io_read_error, &
           "gems_read_latitude: xtrack dimension mismatch", errstat)
      return
    endif

    allocate (tmp_lat(ntimes, lxtrack), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "gems_read_latitude: allocate failed", errstat)
      return
    endif

    call tiof_get2d_r4 (tio_l1obj, "pixel_latitude", [0,tstart], &
         [nxtrack,ntimes], tmp_lat, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "gems_read_lat: failed to read latitudes", &
                       errstat)
      return
    endif
    call tiof_close(tio_l1obj, errstat)

    latitude = reshape(tmp_lat, (/nxtrack,ntimes/),order=flip)
    deallocate (tmp_lat, stat=errstat)

  end subroutine gems_read_latitude


  !> Read GEMS glint and ice flags
  !------------------------------------------------------------------------
  !
  ! @param[in]  l1bfile        L1 radiance file name
  ! @param[in]  nxtrack        Cross-track dimension size
  ! @param[in}  ntimes         Along-track dimension size
  ! @param[out] snow_ice_flag  2D snow/ice flag array
  ! @param[out] glint_flag     2D glint flag array
  ! @param      errstat        error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_ice_glint (l1bfile, nxtrack, ntimes, &
       snow_ice_flag, glint_flag, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile
    integer (kind=4), intent(in) :: nxtrack, ntimes
    !output files
    integer (kind=2), dimension(nxtrack, 0:ntimes-1), intent(out) :: &
         snow_ice_flag, glint_flag
    integer (kind=4), intent(inout) :: errstat
    !local variables
    integer (kind=2), dimension(0:ntimes-1, nxtrack) :: tmp_flg
    integer (kind=2), dimension(nxtrack, 0:ntimes-1) :: tmp_flg2
    type(tiof_file_type) :: tio_l1obj
    integer, dimension(1:2), parameter :: flip = (/2,1/)

    if (errstat /= 0) return

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ground_pixel_quality_flag", [0, 0], &
         [nxtrack,ntimes], tmp_flg, errstat)
    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "gems_read_ice_glint: failed", &
                       errstat)
      return
    endif

    snow_ice_flag=0
    glint_flag=0
    tmp_flg2 = reshape(tmp_flg, (/nxtrack, ntimes/), order=flip)
    snow_ice_flag = int(ibits(tmp_flg2,2,1), kind=2)
    glint_flag = int(ibits(tmp_flg2,1,1), kind=2)

  end subroutine gems_read_ice_glint


  !> Read subset of GEMS radiance and geolocation
  !------------------------------------------------------------------------
  !
  ! @param[in]  l1bfile        L1 radiance file name
  ! @param[in]  iline          starting along-track position
  ! @param[in]  nxtrack        Cross-track dimension size
  ! @param[in}  nloop          Along-track block size
  ! @param[in}  nwavel_ccd     spectral dimension size
  ! @param[out] tmp_spc        3D radiance array
  ! @param[out] tmp_wvl        3D wavelength array
  ! @param[out] tmp_flg        3D pixel quality flag array
  ! @param      errstat        error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_radiance_lines (l1bfile, iline, nxtrack, nloop, &
       nwavel_ccd, tmp_spc, tmp_wvl, tmp_flg, errstat)

    USE OMSAO_omidata_module,  ONLY: omi_height, omi_geoflg, omi_latitude,&
      omi_longitude, omi_szenith, omi_sazimuth, omi_vzenith, omi_vazimuth, &
      omi_time, omi_xtrflg_l1b

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile
    integer (kind=4), intent(in) :: iline, nxtrack, nloop, nwavel_ccd
    !output variables
    real (kind=4), dimension(:,:,:), intent(out) :: tmp_spc, tmp_wvl
    integer (kind=2), dimension(:,:,:), intent(out) :: tmp_flg
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: tio_l1obj
    integer, dimension(1:2), parameter :: flip = (/2,1/)
    integer, dimension(1:3), parameter :: flip3d = (/3,2,1/)
    real (kind=4), dimension(nloop,nxtrack) :: tmp_lat, tmp_lon, tmp_sza, &
         tmp_saa, tmp_vza, tmp_vaa
    real (kind=4), dimension(nxtrack,nwavel_ccd) :: tmp_wvl2
    real (kind=4), dimension(:,:,:), allocatable :: tmp_rad
    integer(kind=2), dimension(:,:,:), allocatable :: tmp_pqf
    integer(kind=2), dimension(nloop,nxtrack) :: tmp_hgt, tmp_gpqf
    integer :: n

    if (errstat /= 0) return

    allocate (tmp_rad(nloop, nxtrack, nwavel_ccd), &
         tmp_pqf(nloop, nxtrack, nwavel_ccd), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "gems_read_radiance_lines: allocation failed", errstat)
      return
    endif

    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get2d_r4 (tio_l1obj, "pixel_latitude", [0, iline], &
         [nxtrack, nloop], tmp_lat, errstat)
    call tiof_get2d_r4 (tio_l1obj, "pixel_longitude", [0, iline], &
         [nxtrack, nloop], tmp_lon, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sun_zenith_angle", [0, iline], &
         [nxtrack, nloop], tmp_sza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sun_azimuth_angle", [0, iline], &
         [nxtrack, nloop], tmp_saa, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sc_zenith_angle", [0, iline], &
         [nxtrack, nloop], tmp_vza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sc_azimuth_angle", [0, iline], &
         [nxtrack, nloop], tmp_vaa, errstat)
    call tiof_get2d_i2 (tio_l1obj, "terrain_height", [0, iline], &
         [nxtrack, nloop], tmp_hgt, errstat)
    call tiof_get1d_r8 (tio_l1obj, "image_acquisition_time", [iline], &
         [nloop], omi_time, errstat)
    call tiof_get3d_r4 (tio_l1obj, "image_pixel_values", [0,0,iline], &
         [nwavel_ccd,nxtrack,nloop], tmp_rad, errstat)
    call tiof_get2d_r4 (tio_l1obj, "wavelength", [0,0], &
         [nwavel_ccd,nxtrack], tmp_wvl2, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ground_pixel_quality_flag", [0,iline], &
         [nxtrack,nloop], tmp_gpqf, errstat)
    call tiof_get3d_i2 (tio_l1obj, "bad_pixel_mask", &
         [0, 0,iline], [nwavel_ccd,nxtrack,nloop], tmp_pqf, errstat)

    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "gems_read_radiance_lines: failed", errstat)
      return
    endif

    ! flip array axes
    omi_latitude(1:nxtrack,0:nloop-1) = &
         reshape(tmp_lat,(/nxtrack,nloop/),order=flip)
    omi_longitude(1:nxtrack,0:nloop-1) = &
         reshape(tmp_lon,(/nxtrack,nloop/),order=flip)
    omi_szenith(1:nxtrack,0:nloop-1) = &
         reshape(tmp_sza,(/nxtrack,nloop/),order=flip)
    omi_sazimuth(1:nxtrack,0:nloop-1) = &
         reshape(tmp_saa,(/nxtrack,nloop/),order=flip)
    omi_vzenith(1:nxtrack,0:nloop-1) = &
         reshape(tmp_vza,(/nxtrack,nloop/),order=flip)
    omi_vazimuth(1:nxtrack,0:nloop-1) = &
         reshape(tmp_vaa,(/nxtrack,nloop/),order=flip)
    omi_height(1:nxtrack,0:nloop-1) = &
         reshape(tmp_hgt,(/nxtrack,nloop/),order=flip)
    omi_geoflg(1:nxtrack,0:nloop-1) = &
         reshape(tmp_gpqf,(/nxtrack,nloop/),order=flip)
    ! variables below output via subroutine interface
    tmp_spc(:,1:nxtrack,1:nloop) = &
         reshape(tmp_rad,(/nwavel_ccd,nxtrack,nloop/),order=flip3d)
    tmp_flg(:,1:nxtrack,1:nloop) = &
         reshape(tmp_pqf,(/nwavel_ccd,nxtrack,nloop/),order=flip3d)
    tmp_wvl(:,1:nxtrack,1) = &
         reshape(tmp_wvl2,(/nwavel_ccd,nxtrack/),order=flip)
    !There has to be a better way to do this:
    do n=2,nloop
      tmp_wvl(:,:,n)=tmp_wvl(:,:,1)
    enddo

    !For now set omi_xtrflg_l1b to zero, since GEMS xtrflg is undefined as yet
    omi_xtrflg_l1b(1:nxtrack,0:nloop-1)=0

    deallocate(tmp_rad, tmp_pqf, stat=errstat)

  end subroutine gems_read_radiance_lines



  !> Read entire GEMS geolocation fields
  !------------------------------------------------------------------------
  !
  ! @param[in]  l1bfile        L1 radiance file name
  ! @param[in]  ntimes         Along-track dimension size
  ! @param[in]  nxtrack        Cross-track dimension size
  ! @param[out] lat            2D latitude array
  ! @param[out] lon            2D longitude array
  ! @param[out] sza            2D solar zenith angle array
  ! @param[out] vza            2D viewing zenith angle array
  ! @param[out] saa            2D solar azimutha angle array
  ! @param[out] vaa            2D viewing azimuth angle array
  ! @param[out] tght           2D terrain height array
  ! @param[out] time           1D time array
  ! @param      errstat        error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_geofields (l1bfile, ntimes, nxtrack, lat, lon, &
       sza, vza, saa, vaa, thgt, time, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile
    integer (kind=4), intent(in) :: ntimes, nxtrack
    !output variables
    real (kind=4), dimension(:,:), intent(out) :: lat, lon, sza, vza, saa, &
         vaa, thgt
    real (kind=8), dimension(:), intent(out) :: time
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: tio_l1obj
    integer, dimension(1:2), parameter :: flip = (/2,1/)
    real (kind=4), dimension(0:ntimes-1,1:nxtrack) :: tmp_lat, tmp_lon, &
         tmp_sza, tmp_saa, tmp_vza, tmp_vaa, tmp_hgt

    if (errstat /= 0) return


    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get2d_r4 (tio_l1obj, "pixel_latitude", [0, 0], &
         [nxtrack, ntimes], tmp_lat, errstat)
    call tiof_get2d_r4 (tio_l1obj, "pixel_longitude", [0, 0], &
         [nxtrack, ntimes], tmp_lon, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sun_zenith_angle", [0, 0], &
         [nxtrack, ntimes], tmp_sza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sun_azimuth_angle", [0, 0], &
         [nxtrack, ntimes], tmp_saa, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sc_zenith_angle", [0, 0], &
         [nxtrack, ntimes], tmp_vza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "sc_azimuth_angle", [0, 0], &
         [nxtrack, ntimes], tmp_vaa, errstat)
    call tiof_get2d_r4 (tio_l1obj, "terrain_height", [0, 0], &
         [nxtrack, ntimes], tmp_hgt, errstat)
    call tiof_get1d_r8 (tio_l1obj, "image_acquisition_time", [0], &
         [ntimes], time, errstat)

    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "gems_read_geofields: failed", errstat)
      return
    endif

    ! flip array axes
    lat = reshape(tmp_lat,(/nxtrack,ntimes/),order=flip)
    lon = reshape(tmp_lon,(/nxtrack,ntimes/),order=flip)
    sza = reshape(tmp_sza,(/nxtrack,ntimes/),order=flip)
    saa = reshape(tmp_saa,(/nxtrack,ntimes/),order=flip)
    vza = reshape(tmp_vza,(/nxtrack,ntimes/),order=flip)
    vaa = reshape(tmp_vaa,(/nxtrack,ntimes/),order=flip)
    thgt = reshape(tmp_hgt,(/nxtrack,ntimes/),order=flip)

  end subroutine gems_read_geofields


  !> Read cloud data from L2 GEMS cloud product file
  !----------------------------------------------------------------------
  !
  ! @param[in]   cloud_file      l2 cloud product file
  ! @param[in]   ntimes          Along-track dimension size
  ! @param[in]   nxtrack         Cross-track dimension size
  ! @param[out]  cloud_fraction  Cloud fraction
  ! @param[out]  cloud_pressure  Cloud pressure
  ! @param       errstat         error tracking integer
  !
  ! NB: GEMS cloud product is TBD, so for now this just returns empty
  !     "missing value" arrays. Outline code in place for when product
  !     becomes available.
  !
  ! @author  E. O'Sullivan  October 2020
  !----------------------------------------------------------------------
  subroutine gems_read_cld (cloud_file, ntimes, nxtrack, cloud_fraction, &
       cloud_pressure, errstat)

    use OMSAO_parameters_module, only: r8_missval

    implicit none

    ! input variables
    character (len=*), intent(in) :: cloud_file
    integer (kind=4), intent(in) :: ntimes, nxtrack
    ! output variables
    real (kind=8), dimension(1:nxtrack, 0:ntimes-1), intent(out) :: &
         cloud_fraction, cloud_pressure
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: cld
    character (len=22) :: cf_name, cp_name, group_name
    real(kind=4), dimension(ntimes, nxtrack) :: gems_cp, gems_cf
    integer, dimension(1:2), parameter :: flip=(/2, 1/)

    if (errstat /= 0) return

    cf_name="EffectiveCloudFraction"
    cp_name="CloudCentroidPressure"
    group_name="/Data Fields"

    call tell_log (1, "Reading GEMS cloud file "//trim(cloud_file) )

    ! Outline code for when GEMS cloud product available
    call tiof_open (cloud_file, cld, nf90_nowrite, errstat)
    call tiof_push_group (cld, trim(group_name), errstat)
    call tiof_get2d_r4 (cld, trim(cp_name), [0,0], [nxtrack,ntimes], &
         gems_cp, errstat)
    call tiof_get2d_r4 (cld, trim(cf_name), [0,0], [nxtrack,ntimes], &
         gems_cf, errstat)
    call tiof_close (cld, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "gems_read_cld: failed", errstat)
      return
    endif

    ! Reshape to match TEMPO dimension order, deal with undeclared fill values
    cloud_pressure(1:nxtrack,0:ntimes-1) = &
         reshape (gems_cp,(/nxtrack,ntimes/),order=flip)
    where (cloud_pressure <= -998.0d0)
      cloud_pressure = r8_missval
    endwhere
    cloud_fraction(1:nxtrack,0:ntimes-1) = &
         reshape (gems_cf,(/nxtrack,ntimes/),order=flip)
    where (cloud_fraction <= -998.0d0)
      cloud_fraction = r8_missval
    endwhere

    ! Force cloud params into physical bounds, ignoring missing values
    where (cloud_fraction > r8_missval .and. cloud_fraction < 0.0d0)
      cloud_fraction = 0.0d0
    elsewhere (cloud_fraction > 1.0d0)
      cloud_fraction = 1.0d0
    endwhere
    where (cloud_pressure > r8_missval .and. cloud_pressure < 0.0d0)
      cloud_pressure = 0.0d0
    endwhere

  end subroutine gems_read_cld


  !> Read Earth-sun distance from GEMS L1C radiance file
  !------------------------------------------------------------------------
  !
  ! @param[in]  l1bfile   L1 radiance file name
  ! @param[out] dist      Earth-sun distance in m
  ! @param      errstat   error tracking integer
  !
  ! @author E. O'Sullivan October 2020
  !------------------------------------------------------------------------
  subroutine gems_read_earth_sun_distance (l1bfile, dist, errstat)

    use netcdf, only: nf90_nowrite, nf90_noerr, nf90_global, nf90_get_att

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile
    !output variables
    real (kind=4), intent(out) :: dist
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: tio_l1obj
    integer :: ncerr
    real (kind=8) :: localdist

    if (errstat /= 0) return

    call tiof_open(l1bfile, tio_l1obj, nf90_nowrite, errstat)
    ncerr = nf90_get_att (tio_l1obj%fileid, nf90_global, "earth_sun_distance",&
         localdist)
    call tiof_close(tio_l1obj, errstat)
    if (errstat /= 0 .or. ncerr /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
           "gems_read_erath_sun_distance: failed", errstat)
      return
    endif

    !Looks like value can be empty in early files, so just in case:
    if (localdist < 100.0d0) then
      dist=149957870700.0
    else
      dist=real(localdist,kind=4)
    endif

  end subroutine gems_read_earth_sun_distance



end module m_read_gems
