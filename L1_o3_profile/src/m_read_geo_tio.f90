!> Subroutines to read in geographic data from L1 radiance netCDF file
module m_read_geo_tio
  use netcdf
  use tio_module
  use tell_module
  use o3p_names_module
  use OMSAO_omidata_module, only: nxtrack_max, ntimes_max

  implicit none

  ! Public parameters
  real (kind=4), dimension (4, 0:nxtrack_max, 0:ntimes_max) :: tio_allclat, &
       tio_allclon

  ! Local parameters used in geometry calculations
  real (kind=8), parameter, private :: pi         = 3.14159265358979d0
  real (kind=8), parameter, private :: pihalf     = 0.5d0  * pi
  real (kind=8), parameter, private :: twopi      = 2.0d0  * pi
  real (kind=8), parameter, private :: deg2rad    = pi / 180.0d0
  real (kind=8), parameter, private :: rad2deg    = 180.0d0 / pi
  real (kind=4), parameter, private :: minza = 0.0, maxza=90.0, &
       minaza = -360., maxaza = 360.0


  public read_geo_tio!, read_geo_line_tio
  private get_sphgeoview_corners, sphergeom_intermediate, circle_rdis, &
       angle_minus_twopi


contains

  !> Read geolocation data from L1 radiance netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] l1file  L1 netCDF radiance file name
  !> @param[in] l1swath L1 netCDF swath name
  !> @param[in] nstep   size of dimension in scan direction
  !> @param[in] nxtrack size of dimention across scan direction
  !> @param[in] nl      number of binned lines in scan direction
  !> @param     errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan    July 2016
  !-----------------------------------------------------------------------
  subroutine read_geo_tio ( l1file, l1swath, nstep, nxtrack, nl, errstat)
    use OMSAO_omidata_module, only:  offset_line, nswath, &
         land_water_flg, glint_flg, snow_ice_flg
    use OMSAO_variables_module, only: use_he5_in, coadd_uv2, nxbin, nybin
    ! replace with variables fined in this module once he5 input obsoleted
    use OMSAO_pixelcorner_module, only: omi_alllat, omi_alllon, omi_allsza, &
       omi_allvza, omi_allaza, omi_allsca, omi_alltime, omi_allGeoFlg, &
       omi_allHeight, omi_allXTrackQFlg, omi_allMflg, &
       omi_allelat, omi_allelon, omi_allclat, omi_allclon, &
       omi_allsza, omi_allvza, omi_allaza, omi_allsca
    use m_convert_coadd, only: coadd_byte_qflgs, convert_gpqualflag_info

    implicit none

    ! input variables
    integer, intent (in) :: nstep, nxtrack, nl
    character (len=*), intent(in) :: l1file,  l1swath
    ! output variables
    integer, intent (inout) :: errstat
    ! local variables
    type (tiof_file_type) :: tio_l1obj
    integer :: eline, sline, nline, iline, nx, i, j, ix, iy, nbits, ndim
    integer :: ysidx, yeidx, ymidx, xsidx, xeidx, xmidx
    real (kind=4), dimension (1:nxtrack, 0:nstep-1) :: tio_lat, &
         tio_lon, tio_sza, tio_vza, tio_saza, tio_vaza
    real (kind=8), dimension (1:nxtrack, 0:nstep-1) :: tio_alllat, &
         tio_alllon
    real (kind=8), dimension (1:nxtrack, 0:nstep-1) :: tmp_sza, &
         tmp_vza, tmp_saza, tmp_vaza
    integer (kind=4), dimension (1:nxtrack_max, 0:ntimes_max-1) :: tio_geoflg
    integer (kind=2), dimension (1:nxtrack_max, 0:ntimes_max-1) :: tio_height
    integer (kind=1), dimension (1:nxtrack_max, 0:ntimes_max-1) :: tio_xtrackqflg, tmp_xtrackqflg
!    integer (kind=i1), dimension (1:nxtrack_max)      :: tio_xtrackqflg1
    integer (kind=2), dimension (0:ntimes_max-1) :: tio_mflg
    real (kind=8), dimension (0:ntimes_max-1)    :: tio_time
    real (kind=8), dimension (0:nxtrack_max, 0:ntimes_max-1) :: tio_allelat, &
         tio_allelon
    real (kind=4), dimension (4, 0:nxtrack_max, 0:ntimes_max) :: tmp_allclon, &
         tmp_allclat
    real (kind=8), dimension (0:nxtrack_max, 0:ntimes_max) :: dummy_clat, &
         dummy_clon, dummy_allclat, dummy_allclon
    real (kind=8), dimension (1:nxtrack_max, 0:ntimes_max-1) :: tio_allsza, &
         tio_allvza, tio_allaza, tio_allsca

    if (errstat /= 0) return

    ! limits of region to be read in
    sline = offset_line
    eline  = offset_line + nl * nybin - 1
    if (eline == sline) eline = sline + 1
    nline = eline - sline + 1

    ! read in geolocation data for the chosen subset of step positions
    call tiof_open (l1file, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, o3p_var_time, [sline], [nline], &
         tio_time(sline:eline), errstat)
    call tiof_push_group (tio_l1obj, l1swath, errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_latitude, [sline,0], &
         [nline,nxtrack], tio_lat(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_longitude, [sline,0], &
         [nline,nxtrack], tio_lon(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_sz_angle, [sline,0], &
         [nline,nxtrack], tio_sza(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_sa_angle, [sline,0], &
         [nline,nxtrack], tio_saza(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_vz_angle, [sline,0], &
         [nline,nxtrack], tio_vza(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_va_angle, [sline,0], &
         [nline,nxtrack], tio_vaza(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_i2 (tio_l1obj, o3p_var_terrain_height, [sline,0], &
         [nline,nxtrack], tio_height(1:nxtrack, sline:eline), errstat)
    call tiof_get2d_i4 (tio_l1obj, o3p_var_geoflg, [sline,0], &
         [nline,nxtrack], tio_geoflg(1:nxtrack, sline:eline), errstat)
    call tiof_get1d_i2 (tio_l1obj, o3p_var_mqf, [sline], &
         [nline], tio_mflg(sline:eline), errstat)
    call tiof_get2d_i1 (tio_l1obj, o3p_var_anomflg, [sline,0], &
         [nline,nxtrack], tio_xtrackqflg(1:nxtrack,sline:eline), errstat)
    call tiof_get3d_r4 (tio_l1obj,o3p_var_latitude_bounds, [sline,0,0], &
         [nline,nxtrack,4], tmp_allclat(:,1:nxtrack,sline:eline), errstat)
    call tiof_get3d_r4 (tio_l1obj,o3p_var_longitude_bounds, [sline,0,0], &
         [nline,nxtrack,4], tmp_allclon(:,1:nxtrack,sline:eline), errstat)
    call tiof_pop_group (tio_l1obj, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_geo_data: failed to read geolocation data", &
           errstat)
      return
    endif

    ! coadd and correct any missing xtrack quality flags
    nbits = 8
    ndim = 1
    if (nswath == 2 .AND. coadd_uv2) then
      do iline = sline, eline
        do ix = 1, nxtrack
          i = ix * 2 - 1
          j = i + 1
          tmp_xtrackqflg = tio_xtrackqflg
          !CALL coadd_byte_qflgs(nbits, ndim, tmp_xtrackqflg(i, iline))!, &
               !tmp_xtrackqflg(j, iline)) ! by someone by TEMPO teams
          CALL coadd_byte_qflgs(nbits, ndim, tmp_xtrackqflg(i, iline), tmp_xtrackqflg(j, iline)) ! returned by jbak
          tio_xtrackqflg(ix, iline) = tmp_xtrackqflg(i, iline)
        end do
      end do
    endif
    where (tio_xtrackqflg(1:nxtrack, sline:eline) == -127)
      tio_xtrackqflg(1:nxtrack, sline:eline) = 0
    end where


    ! Resample corners coordinates from L1 file to account for binning
    nx = nxtrack / nxbin
    do i = 0, nx-1
      do j = 0, nl-1
        tio_allclon(1,i,j) = tmp_allclon(1,((i+1)*nxbin),&
             sline+((j+1)*nybin)-1)
        tio_allclon(2,i,j) = tmp_allclon(2,((i+1)*nxbin),sline+(j*nybin))
        tio_allclon(3,i,j) = tmp_allclon(3,1+(i*nxbin),sline+(j*nybin))
        tio_allclon(4,i,j) = tmp_allclon(4,1+(i*nxbin),&
             sline+((j+1)*nybin)-1)
        tio_allclat(1,i,j) = tmp_allclat(1,((i+1)*nxbin),&
             sline+((j+1)*nybin)-1)
        tio_allclat(2,i,j) = tmp_allclat(2,((i+1)*nxbin),sline+(j*nybin))
        tio_allclat(3,i,j) = tmp_allclat(3,1+(i*nxbin),sline+(j*nybin))
        tio_allclat(4,i,j) = tmp_allclat(4,1+(i*nxbin),&
             sline+((j+1)*nybin)-1)
      end do
    end do

    ! FIXME - note that the geometric calculations below probably need
    ! to be re-thought for TEMPO. For now we keep them for testing

    ! Add small offset to SZA if it is zero, to avoid VLIDORT failure
    do i = 1, nxtrack
      do j = sline,eline
        if (tio_sza(i,j) .lt. 1e-5) then
          tio_sza(i,j) = tio_sza(i,j) + 1e-5
        endif
        if (tio_vza(i,j) .lt. 1e-5) then
          tio_vza(i,j) = tio_vza(i,j) + 1e-5
        endif
      enddo
    enddo

    ! compute pixel corner coordinates for the spatially coadded pixels
    ! get_sphgeoview_corners works in r8, so need temporary arrays
    ! dummy corner values, not used except to check calculation against OMI
    tio_alllon(1:nxtrack,sline:eline) = tio_lon(1:nxtrack,sline:eline)
    tio_alllat(1:nxtrack,sline:eline) = tio_lat(1:nxtrack,sline:eline)
    tmp_sza(1:nxtrack,sline:eline) = tio_sza(1:nxtrack,sline:eline)
    tmp_saza(1:nxtrack,sline:eline) = tio_saza(1:nxtrack,sline:eline)
    tmp_vza(1:nxtrack,sline:eline) = tio_vza(1:nxtrack,sline:eline)
    tmp_vaza(1:nxtrack,sline:eline) = tio_vaza(1:nxtrack,sline:eline)
    call get_sphgeoview_corners (nxtrack, nline, &
         tio_alllon(:, sline:eline), &
         tio_alllat(:, sline:eline), &
         tmp_sza(:, sline:eline), &
         tmp_saza(:, sline:eline), &
         tmp_vza(:, sline:eline), &
         tmp_vaza(:, sline:eline), &
         dummy_clon(0:nxtrack, sline:eline+1), &
         dummy_clat(0:nxtrack, sline:eline+1), &
         tio_allelon(0:nxtrack, sline:eline), &
         tio_allelat(0:nxtrack, sline:eline), &
         tio_allsza(1:nxtrack, sline:eline), &
         tio_allvza(1:nxtrack, sline:eline), &
         tio_allaza(1:nxtrack, sline:eline), &
         tio_allsca(1:nxtrack, sline:eline))

    ! move data within arrays
    i = sline
    j = i + nl - 1

    tio_alllon (1:nx, 0:nl-1)  = tio_alllon (1:nx, i:j)
    tio_alllat (1:nx, 0:nl-1)  = tio_alllat (1:nx, i:j)
    dummy_allclon(0:nx, 0:nl)  = dummy_clon(0:nx, i:j+1)
    dummy_allclat(0:nx, 0:nl)  = dummy_clat(0:nx, i:j+1)
    tio_allelon(0:nx, 0:nl-1)  = tio_allelon(0:nx, i:j)
    tio_allelat(0:nx, 0:nl-1)  = tio_allelat(0:nx, i:j)
    tio_allsza (1:nx, 0:nl-1)  = tio_allsza (1:nx, i:j)
    tio_allvza (1:nx, 0:nl-1)  = tio_allvza (1:nx, i:j)
    tio_allaza (1:nx, 0:nl-1)  = tio_allaza (1:nx, i:j)
    tio_allsca (1:nx, 0:nl-1)  = tio_allsca (1:nx, i:j)



    ! determine time and flags for coadded pixels
    do iy = 0, nl - 1
      ysidx = sline + iy * nybin
      yeidx = ysidx + nybin - 1
      ymidx = ysidx + nybin / 2
      tio_time(iy) = sum(tio_time(ysidx:yeidx)) / nybin
      ! Use those from the middle point
      ! (avoid dealing with polar, dateline regions)
      tio_mflg(iy) = tio_mflg(ymidx)
      ! get separate land/water, glint, snow/ice flags
      if (.NOT. use_he5_in) then
        call convert_gpqualflag_info (nxtrack, &
             omi_allGeoFlg(1:nxtrack, ymidx), &
             land_water_flg(1:nxtrack, ymidx), glint_flg(1:nxtrack, ymidx), &
             snow_ice_flg(1:nxtrack, ymidx))
      endif

      do ix = 1, nx
        xsidx = (ix - 1) * nxbin + 1
        xeidx = xsidx + nxbin - 1
        xmidx = xsidx + nxbin / 2

        tio_height(ix, iy)  = int( &
             sum(1.0 * tio_height(xsidx:xeidx, ysidx:yeidx)) &
             / (1.0 * nxbin * nybin), kind=2)

        tio_geoflg(ix, iy) = tio_geoflg(xmidx, ymidx)

        tio_xtrackqflg(ix, iy) = tio_xtrackqflg(xmidx, ymidx)

      enddo
    enddo

    ! if using he5 input, check values are identical
    if (use_he5_in) then
      do iy = 0, nl-1
        if (tio_time(iy).ne.omi_alltime(iy)) print *, 'mismatch: time'
        if (tio_mflg(iy).ne.omi_allMflg(iy)) print *, 'mismatch: mflg'
        do ix = 1, nx
          if (tio_height(ix, iy).ne.omi_allHeight(ix, iy)) &
               print *, 'mismatch: height'
!          if (tio_geoflg(ix, iy).ne.omi_allGeoFlg(ix, iy)) &
!               print *, 'mismatch: geoflg', &
!               ix, iy, tio_geoflg(ix, iy), omi_allGeoFlg(ix, iy)
          if (tio_xtrackqflg(ix, iy).ne.omi_allXtrackQFlg(ix, iy)) &
               print *, 'mismatch: xtrackqflg'
          if (tio_alllat(ix, iy).ne.omi_alllat(ix, iy)) &
               print *, 'mismatch: lat', tio_alllat(ix, iy), omi_alllat(ix, iy)
          if (tio_alllon(ix, iy).ne.omi_alllon(ix, iy)) &
               print *, 'mismatch: lon', tio_alllon(ix, iy), omi_alllon(ix, iy)
          if (tio_allsza(ix, iy).ne.omi_allsza(ix, iy)) &
               print *, 'mismatch: sza'
          if (tio_allaza(ix, iy).ne.omi_allaza(ix, iy)) &
               print *, 'mismatch: aza'
          if (tio_allvza(ix, iy).ne.omi_allvza(ix, iy)) &
               print *, 'mismatch: vza'
          if (tio_allsca(ix, iy).ne.omi_allsca(ix, iy)) &
               print *, 'mismatch: sca'
        enddo
        do ix = 0, nx
          if (dummy_allclat(ix, iy).ne.omi_allclat(ix, iy)) &
               print *, 'mismatch: clat'
          if (dummy_allclon(ix, iy).ne.omi_allclon(ix, iy)) &
               print *, 'mismatch: clon'
          if (tio_allelat(ix, iy).ne.omi_allelat(ix, iy)) &
               print *, 'mismatch: elat'
          if (tio_allelon(ix, iy).ne.omi_allelon(ix, iy)) &
               print *, 'mismatch: elon'
        enddo
      enddo
    endif

    ! move tio values into omi variables
    omi_alltime = tio_time
    omi_allMflg = tio_mflg
    omi_allHeight = tio_height
    omi_allGeoFlg = INT(tio_geoflg, kind=2)
    omi_allXtrackQFlg = tio_xtrackqflg
    omi_alllat(1:nx,0:nl-1) = tio_alllat(1:nx,0:nl-1)
    omi_alllon(1:nx,0:nl-1) = tio_alllon(1:nx,0:nl-1)
    omi_allsza = tio_allsza
    omi_allaza = tio_allaza
    omi_allvza = tio_allvza
    omi_allsca = tio_allsca
    omi_allclat = dummy_allclat
    omi_allclon = dummy_allclon
    omi_allelat = tio_allelat
    omi_allelon = tio_allelon


  end subroutine read_geo_tio


  ! FIXME - geometry probably needs checking for TEMPO
  !
  !> Calculate parameter values at binned pixel corners
  !---------------------------------------------------------------------
  !
  !> @param[in] nxtrack dimension size across track
  !> @param[in] times dimension size along track
  !> @param lon longitude array
  !> @param lat latitude array
  !> @param sza solar zenith angle
  !> @param saza solar azimuth angle
  !> @param vza viewing zenith angle
  !> @param vaza viewing azimuth angle
  !> @param[out] clat corner latitudes
  !> @param[out] clon corner longitudes
  !> @param[out] elat edge latitudes?
  !> @param[out] elon edge longitudes?
  !> @param[out] esza edgs solar zenith angles?
  !> @param[out] eaza edge solar azimuth angles?
  !> @param[out] evza edge viewing zenith angles?
  !> @param[out] esca edge
  !
  !---------------------------------------------------------------------

  subroutine get_sphgeoview_corners (nxtrack, ntimes, lon, lat, sza, saza, &
       vza, vaza, clon, clat, elon, elat, esza, evza, eaza, esca)

    use OMSAO_variables_module, only: nxbin, nybin
    use m_ezspline_interpolation, only: interpol
    use m_angle_sat2toa

    implicit none

    ! ---------------
    ! Input variables
    ! ---------------
    integer (kind=4), intent(in)    :: nxtrack, ntimes
    real (kind=8), dimension (1:nxtrack, 0:ntimes-1), intent(inout) :: lon, &
         lat, sza, saza, vza, vaza
    real (kind=8), dimension (0:nxtrack, 0:ntimes), intent(out) :: clon, clat
    real (kind=8), dimension (0:nxtrack, 0:ntimes-1), intent(out) :: elon, elat
    real (kind=8), dimension (1:nxtrack, 0:ntimes-1), intent(out) :: esza, &
         evza, eaza, esca

    ! ---------------
    ! Local variables
    ! ---------------
    integer (kind=4) :: i, j, jj, ix, mpix, nx, ny
    integer                                             :: errstat
    real (kind=8), dimension (1:nxtrack,0:ntimes-1) :: omixsize
    real (kind=8), dimension (0:nxtrack, 0:ntimes-1):: edsza, edsazm, &
         edvza, edvazm
    real (kind=8), dimension (1:nxtrack)  :: tmpdisx, xsize, tmpxmid
    real (kind=8), dimension (0:nxtrack)  :: tmpx, tmpsza, tmpvza, &
         tmpsaza, tmpvaza
    real (kind=8), dimension (0:ntimes-1) :: tmpdisy
    real (kind=8), dimension (3) :: zen0, zen, sazm, vazm, relaza

    ! ------------------------------------------------------------------
    ! Convert geolocation to radians; do everything in R8 rather than R4
    ! ------------------------------------------------------------------

    lon = lon * deg2rad
    lat = lat * deg2rad

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    clon   = -999.9d0
    clat = -999.9d0
    elon   = -999.9d0
    elat = -999.9d0
    esza   = -999.9d0
    evza = -999.9d0
    eaza = -999.9d0
    esca = -999.9d0

    ! Perform interpolation across the track
    do i = 0, ntimes - 1
      ! Compute the distances between two pixels: (x1 + x2) / 2.
      do ix = 1, nxtrack - 1
        tmpdisx(ix) = circle_rdis(lat(ix, i), lon(ix, i), lat(ix+1, i), &
             lon(ix+1, i))
      enddo

      ! Compute the pixel size across the track
      ! Assume the center two pixels have equal pixel size
      ! (which causes about < 0.1 km error for UV-2)
      mpix = nxtrack / 2
      xsize(mpix) = tmpdisx(mpix) / 2.0
      xsize(mpix + 1) = xsize(mpix)

      do ix = mpix -1, 1, -1
        xsize(ix) = tmpdisx(ix) - xsize(ix + 1)
      enddo
      do ix = mpix + 2, nxtrack
        xsize(ix) = tmpdisx(ix - 1) - xsize(ix - 1)
      enddo
      omixsize(:, i) = xsize

      !!  This is to test SUBROUTINE sphergeom_intermediate
      !!  Works for both interpolation and extrapolation
      !! (with certain limitation)

      ! Perform interpolation
      do ix = 1, nxtrack - 1
        call sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), &
             lon(ix+1, i), tmpdisx(ix), xsize(ix), elat(ix, i), elon(ix, i))
      enddo
      ix = 1
      call sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), &
           lon(ix+1, i), tmpdisx(ix), -xsize(ix), elat(ix-1, i), elon(ix-1, i))
      ix = nxtrack
      call sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix-1, i), &
           lon(ix-1, i), tmpdisx(ix-1), -xsize(ix), elat(ix, i), elon(ix, i))

      ! Compute viewing geometry for west and east edge
      ! Performal interpolation/extrapolation (2 points) along the
      ! spherical lines (good enough)
      tmpx = 0.0
      do ix = 1, nxtrack
        tmpx(ix) = tmpx(ix-1) + xsize(ix)
      enddo
      tmpxmid(1:nxtrack) = (tmpx(0:nxtrack-1) + tmpx(1:nxtrack)) / 2.0
      call interpol(tmpxmid(1:nxtrack), sza(1:nxtrack, i),  nxtrack, &
           tmpx(0:nxtrack), tmpsza(0:nxtrack), nxtrack+1, errstat)
      call interpol(tmpxmid(1:nxtrack), saza(1:nxtrack, i), nxtrack, &
           tmpx(0:nxtrack), tmpsaza(0:nxtrack), nxtrack+1, errstat)
      call interpol(tmpxmid(1:nxtrack), vza(1:nxtrack, i),  nxtrack, &
           tmpx(0:nxtrack), tmpvza(0:nxtrack),  nxtrack+1, errstat)
      call interpol(tmpxmid(1:nxtrack), vaza(1:nxtrack, i), nxtrack, &
           tmpx(0:nxtrack), tmpvaza(0:nxtrack), nxtrack+1, errstat)

      ! Check center pixel
      do ix = mpix-1, mpix + 1
        if (tmpvaza(ix) < 0) then
          tmpvaza(ix) = -tmpvaza(ix-1)
          exit
        endif
      enddo

      edsza(0:nxtrack, i) = tmpsza(0:nxtrack)
      edsazm(0:nxtrack, i) = tmpsaza(0:nxtrack)
      edvza(0:nxtrack, i) = tmpvza(0:nxtrack)
      edvazm(0:nxtrack, i) = tmpvaza(0:nxtrack)
    enddo

    ! Perform interpolation along the track with a simipler (center) approach
    ! since the pixel size along the track does not vary much
    !FIXME - hard coded OMI pixel size
    if (ntimes == 1) tmpdisy(0) = 0.00212031  ! ~ 13.5 km
    do ix = 0, nxtrack
      do i = 0, ntimes - 2
        tmpdisy(i) = circle_rdis(elat(ix, i), elon(ix, i), elat(ix, i+1), &
             elon(ix, i+1))
        call sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1), &
             elon(ix, i+1), tmpdisy(i), tmpdisy(i)*0.5, clat(ix, i+1), &
             clon(ix, i+1))
      enddo
      i = 0
      call sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1), &
           elon(ix, i+1), tmpdisy(i), -tmpdisy(i)*0.5, clat(ix, i), &
           clon(ix, i))

      i = ntimes - 1
      call sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i-1),  &
           elon(ix, i-1), tmpdisy(i-1), -tmpdisy(i-1)*0.5, clat(ix, i+1), &
           clon(ix, i+1))
    enddo

    ! Perform coadding
    if (nxbin > 1 .or. nybin > 1) then
      nx = nxtrack / nxbin
      ny = ntimes  / nybin

      ! cornor coordinates (only need sampling)
      j = 0
      do ix = 0, nxtrack, nxbin
        clon(j, :) = clon(ix, :)
        clat(j, :) = clat(ix, :)
        j = j + 1
      enddo

      j = 0
      do i = 0, ntimes, nybin
        clon(0:nx, j) = clon(0:nx, i)
        clat(0:nx, j) = clat(0:nx, i)
        j = j + 1
      enddo

      ! edge coordinates (easy to be re-computed from corner coordinates)
      do ix = 0, nx
        do i = 0, ny - 1
          tmpdisy(i) = circle_rdis(clat(ix, i), clon(ix, i), clat(ix, i+1), &
               clon(ix, i+1))
          call sphergeom_intermediate(clat(ix, i), clon(ix, i), &
               clat(ix, i+1), clon(ix, i+1), &
               tmpdisy(i), tmpdisy(i)*0.5, elat(ix, i), elon(ix, i))
        enddo
      enddo

      ! Center coordinates (computed from edge coordinates)
      do ix = 1, nx
        do i = 0, ny - 1
          tmpdisx(ix) = circle_rdis(elat(ix-1, i), elon(ix-1, i), &
               elat(ix, i), elon(ix, i))
          call sphergeom_intermediate(elat(ix-1, i), elon(ix-1, i), &
               elat(ix, i), elon(ix, i), &
               tmpdisx(ix), tmpdisx(ix)*0.5, lat(ix, i), lon(ix, i))
        enddo
      enddo

      ! Average edge viewing geometries along the track, sample along the track
      i = 0
      do ix = 0, nxtrack, nxbin
        jj = 0
        do j = 0, ntimes-1, nybin
          edsza (i, jj) = sum(edsza (ix, j:j+nybin-1)) / nybin
          edsazm(i, jj) = sum(edsazm(ix, j:j+nybin-1)) / nybin
          edvza (i, jj) = sum(edvza (ix, j:j+nybin-1)) / nybin
          edvazm(i, jj) = sum(edvazm(ix, j:j+nybin-1)) / nybin
          jj = jj + 1
        enddo
        i = i + 1
      enddo

      ! Compute center viewing geometries (interpolate across the track)
      do i = 0, ny-1
        tmpx = 0.0
        do ix = 1, nx
          tmpx(ix) = tmpx(ix-1) + circle_rdis(elat(ix-1, i), elon(ix-1, i), &
               elat(ix, i), elon(ix, i))
        enddo
        tmpxmid(1:nx) = (tmpx(0:nx-1) + tmpx(1:nx)) / 2.0

        call interpol(tmpx(0:nx), edsza (0:nx, i),  nx+1, tmpxmid(1:nx), &
             sza (1:nx, i),  nx, errstat)
        call interpol(tmpx(0:nx), edsazm(0:nx, i),  nx+1, tmpxmid(1:nx), &
             saza(1:nx, i),  nx, errstat)

        call interpol(tmpx(0:nx), edvza (0:nx, i),  nx+1, tmpxmid(1:nx), &
             vza (1:nx, i),  nx, errstat)
        call interpol(tmpx(0:nx), edvazm(0:nx, i),  nx+1, tmpxmid(1:nx), &
             vaza(1:nx, i),  nx, errstat)

        ! Check center pixel
        mpix = nx / 2
        do ix = mpix - 1, mpix + 1
          if (vaza(ix, i) < 0) then
            vaza(ix, i) = -vaza(ix-1, i)
            exit
          endif
        enddo
      enddo
    else
      nx = nxtrack
      ny = ntimes
    endif

    ! Now compute effective viewing geometry
    ! Compute effective viewing geometry for each pixel
    ! (at a certain atmosphere)
    do i = 0, ny-1
      do ix = 1, nx
        if ( sza(ix, i)  >= minza  .and. sza(ix, i)  < maxza  .and. &
             vza(ix, i)  >= minza  .and. vza(ix, i)  < maxza  .and. &
             saza(ix, i) >= minaza .and. saza(ix, i) < maxaza .and. &
             vaza(ix, i) >= minaza .and. vaza(ix, i) < maxaza) then
          zen0(1) = edsza(ix-1, i)
          zen0(2) = sza(ix, i)
          zen0(3) = edsza(ix, i)
          zen(1)  = edvza(ix-1, i)
          zen(2)  = vza(ix, i)
          zen(3)  = edvza(ix, i)
          sazm(1) = edsazm(ix-1, i)
          sazm(2) = saza(ix, i)
          sazm(3) = edsazm(ix, i)
          vazm(1) = edvazm(ix-1, i)
          vazm(2) = vaza(ix, i)
          vazm(3) = edvazm(ix, i)
          relaza  = abs(vazm - sazm)

          where (vazm > 0)
            zen = - zen
          end where

          ! -180 < relaza < 180
          where (relaza > 180.0)
            relaza = 360.0 - relaza
          end where

          call omi_angle_sat2toa (3, zen0, zen, relaza, esza(ix, i), &
               evza(ix, i), eaza(ix, i), esca(ix, i) )
        endif
      enddo
    enddo

    clon(0:nx, 0:ny)   = clon(0:nx, 0:ny)   * rad2deg
    clat(0:nx, 0:ny)   = clat(0:nx, 0:ny)   * rad2deg
    lon(1:nx, 0:ny-1)  = lon(1:nx, 0:ny-1)  * rad2deg
    lat(1:nx, 0:ny-1)  = lat(1:nx, 0:ny-1)  * rad2deg
    elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) * rad2deg
    elat(0:nx, 0:ny-1) = elat(0:nx, 0:ny-1) * rad2deg

    where (clon(0:nx, 0:ny) > 180.0)
      clon(0:nx, 0:ny) = clon(0:nx, 0:ny) - 360.0
    end where
    where (clon(0:nx, 0:ny) < -180.0)
      clon(0:nx, 0:ny) = clon(0:nx, 0:ny) + 360.0
    end where

    where (lon(1:nx, 0:ny-1) > 180.0)
      lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) - 360.0
    end where
    where (lon(1:nx, 0:ny-1) < -180.0)
      lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) + 360.0
    end where

    where (elon(0:nx, 0:ny-1) > 180.0)
      elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) - 360.0
    end where
    where (elon(0:nx, 0:ny-1) < -180.0)
      elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) + 360.0
    end where


    return

  end subroutine get_sphgeoview_corners



  subroutine sphergeom_intermediate ( lat1, lon1, lat2, lon2, c0, c, lat, lon )

    ! -----------------------------------------------------------------
    ! Finds the co-ordinates of C the baseline extended from two
    ! lon/lat points (A, B) on a sphere given the hypotenuse C_IN.
    ! ----------------------------------------------------------------
    implicit none

    ! ---------------
    ! Input variables
    ! ---------------
    real (kind=8),    intent (IN) :: lat1, lat2, lon1, lon2, c0, c

    ! ----------------
    ! Output variables
    ! ----------------
    real (kind=8),  intent (OUT)  :: lat, lon

    ! ---------------
    ! Local variables
    ! ---------------
    real (kind=8)  :: x, y, z, tmp1, tmp2, frc, gamsign, theta, gam0, gam
    real (kind=8) :: one=1.0d0

    lat = 0.0d0
    lon = 0.0d0
    gam0 = angle_minus_twopi ( lon2 - lon1, pi )
    !gamsign = abs(gam0) / gam0
    ! version above fails for gam0==0
    gamsign = sign(one,gam0)
    gam0 = abs(gam0)


    ! Get straight line (AB) segment fraction frc intercepted by the line from center to C
    ! If frc < 0, extrapolation, but it is limited to |c| < (180-c0)/2.0
    tmp1 = sin(c)
    frc = tmp1 / (sin(c0 - c) + tmp1)

    ! Work in Cartesian Coordinate
    tmp1 = frc * cos(lat2)
    tmp2 = 1.0 - frc
    x = tmp2 *   cos(lat1) + tmp1 * cos(gam0)
    y = tmp1 *   sin(gam0)
    z = tmp2 *   sin(lat1) + frc * sin(lat2)

    gam = atan(y/x)                          ! -90 < gam < 90
    if (frc >= 0) then
      if (gam < 0) gam = gam + pi           ! 0 <= gam <= 180
    else
      if (gam > 0) gam = gam - pi           ! -180 <= gam <= 0
    endif

    gam = gamsign * gam                      ! Get correct sign

    lon = gam + lon1

    theta = atan (sqrt(x**2 + y**2) / z)     ! -90 < theta < 90
    if (theta < 0) theta = theta + pi        ! 0 <= theta <= 180
    lat = pihalf - theta                     ! Convert to latitude

    return
  end subroutine sphergeom_intermediate



  function circle_rdis(lat1, lon1, lat2, lon2) result(rdis)

    implicit none

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    real (kind=8), intent (IN) :: lat1, lon1, lat2, lon2
    real (kind=8)              :: rdis

    ! Local variable
    real (kind=8)              :: dlon, dlat, a

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = min(1.0, sqrt( sin(dlat/2.0)**2.0 + cos(lat1) * cos(lat2) * &
         sin(dlon/2.0)**2.0 )  )
    rdis = 2.0 * asin(a)           ! relative distance in radiances

    return

  end function circle_rdis




  real (kind=8) function angle_minus_twopi ( gamma0, pival ) result ( gamma )

    implicit none

    ! ---------------
    ! Input variables
    ! ---------------
    real (kind=8), intent (IN) :: gamma0, pival

    if ( gamma0 > pival ) then
      gamma = gamma0 - 2.0d0 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    else if ( gamma0 < -pival ) then
      gamma = gamma0 + 2.0d0 * pival
    else
      gamma = gamma0
    end if

    return
  end function angle_minus_twopi



end module m_read_geo_tio
