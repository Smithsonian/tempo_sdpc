!> Subroutines to read in geographic data from L1 radiance netCDF file
module m_read_geo_tio
  
  USE netcdf
  USE tio_module
  USE tell_module
  USE o3p_names_module
  USE OMSAO_tmpodata_module, only: nxtrack_max, ntimes_max

  IMPLICIT NONE

  ! Local PARAMETERs used in geometry calculations
  REAL (kind=8), PARAMETER, PRIVATE :: pi         = 3.14159265358979d0
  REAL (kind=8), PARAMETER, PRIVATE :: pihalf     = 0.5d0  * pi
  REAL (kind=8), PARAMETER, PRIVATE :: twopi      = 2.0d0  * pi
  REAL (kind=8), PARAMETER, PRIVATE :: deg2rad    = pi / 180.0d0
  REAL (kind=8), PARAMETER, PRIVATE :: rad2deg    = 180.0d0 / pi
  REAL (kind=4), PARAMETER, PRIVATE :: minza = 0.0, maxza=90.0, &
       minaza = -360., maxaza = 360.0


  PUBLIC read_geo_tio!, read_geo_line_tio
  PRIVATE get_sphgeoview_corners, sphergeom_intermediate, circle_rdis, &
       angle_minus_twopi, convert_gpqualflag_info


CONTAINS

  !> Read geolocation data from L1 radiance netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] l1file  L1 netCDF radiance file name
  !> @param[in] l1swath L1 netCDF swath name
  !> @param[in] ntimes   size of dimension in scan direction
  !> @param[in] nxtrack size of dimention across scan direction
  !> @param[in] nl      number of binned lines in scan direction
  !> @param     errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan    July 2016
  !-----------------------------------------------------------------------
  SUBROUTINE read_geo_tio (l1swath, geo,ntimes, nxtrack, & 
                           spix, lpix, sline, eline,do_deallo,  errstat)
    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY: geo_group, nxbin, nybin, & 
                                      l1file=>l1b_rad_filename, szamax

    IMPLICIT NONE 

    ! input variables
    LOGICAL, INTENT (in) :: do_deallo
    INTEGER, INTENT (in) :: ntimes, nxtrack
    INTEGER, INTENT(inout) ::  spix, lpix, sline, eline
    CHARACTER (len=*), INTENT(in) :: l1swath
    ! output variables
    TYPE (geo_group), INTENT(OUT) :: geo
    INTEGER, INTENT (out) :: errstat
    ! local variables
    type (tiof_file_type) :: tio_l1obj
    INTEGER :: nline, iline, nx, i, j, ix, iy,iix,iiy, nbits, ndim, nl, sline1, eline1, nbin
    INTEGER :: ysidx, yeidx, ymidx, xsidx, xeidx, xmidx
    !-----------------------------------------------------------
    ! variables for reading original tempo geolocation dataset
    !------------------------------------------------------------
    real (kind=8), dimension (:), POINTER    :: tio_time
    real (kind=4), dimension (:,:), POINTER :: tio_lat, &
                    tio_lon, tio_sza, tio_vza, tio_saza, tio_vaza
    INTEGER (kind=2), dimension (:,:), POINTER :: tio_height
    INTEGER (kind=4), dimension (:,:), POINTER :: tio_geoflg
    real (kind=4), dimension (:,:,:), POINTER:: tio_clon, tio_clat
    !---------------------------------------------------------
    ! variables for binning/arrange tempo geolocation dataset
    !--------------------------------------------------------
    INTEGER (kind=i2), dimension(1:nxtrack_max) ::land_water_flg, glint_flg, snow_ice_flg
    real (kind=8) :: sza, vza,sazm, vazm, relaza
    LOGICAL, save :: first =.true.

    !IF (first) THEN 
       allocate(tio_clon (4, nxtrack, 0:ntimes-1) , tio_clat(4, nxtrack, 0:ntimes-1) )
       allocate(tio_lon ( nxtrack, 0:ntimes-1) , tio_lat(nxtrack, 0:ntimes-1))
       allocate(tio_sza ( nxtrack, 0:ntimes-1) , tio_saza(nxtrack, 0:ntimes-1))
       allocate(tio_vza ( nxtrack, 0:ntimes-1) , tio_vaza(nxtrack, 0:ntimes-1))
       allocate(tio_height ( nxtrack, 0:ntimes-1) , tio_geoflg(nxtrack, 0:ntimes-1))
       allocate(tio_time (0:ntimes-1))
    !   first = .false.
    !ENDIF
   
    !------------------------------------------------------------------------
    ! Initialize
    !------------------------------------------------------------------------

    errstat = 0
    tio_sza = -999
    ! limits of region to be read in
    nl = eline - sline +1
    sline1 = sline
    eline1 = eline
    if (eline == sline) eline1 = sline1 + 1       
    IF (eline1 > ntimes) THEN 
        eline1 = eline1-1 ; sline1=sline1-1
    ENDIF
    sline1 = sline1-1
    eline1 = eline1-1
    nline = eline1 - sline1 + 1
    ! read in geolocation data for the chosen subset of step positions
    call tiof_open (l1file, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, o3p_var_time, [sline1], [nline], &
                        tio_time(sline1:eline1), errstat)
    call tiof_push_group (tio_l1obj, l1swath, errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_latitude, [sline1,0],[nline,nxtrack], &
                        tio_lat (1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_longitude,[sline1,0],[nline,nxtrack], &
                        tio_lon (1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_sz_angle, [sline1,0],[nline,nxtrack], &
                        tio_sza (1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_sa_angle, [sline1,0],[nline,nxtrack], &
                        tio_saza(1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_vz_angle, [sline1,0],[nline,nxtrack], &
                        tio_vza (1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_r4 (tio_l1obj, o3p_var_va_angle, [sline1,0],[nline,nxtrack], &
                        tio_vaza(1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_i2 (tio_l1obj, o3p_var_terrain_height, [sline1,0],  [nline,nxtrack], &
                        tio_height(1:nxtrack,sline1:eline1), errstat)
    call tiof_get2d_i4 (tio_l1obj, o3p_var_geoflg,         [sline1,0],  [nline,nxtrack], &
                        tio_geoflg(1:nxtrack,sline1:eline1), errstat)
    call tiof_get3d_r4 (tio_l1obj,o3p_var_latitude_bounds, [sline1,0,0],[nline,nxtrack,4],&
                        tio_clat(:,1:nxtrack,sline1:eline1), errstat)
    call tiof_get3d_r4 (tio_l1obj,o3p_var_longitude_bounds,[sline1,0,0],[nline,nxtrack,4],&
                        tio_clon(:,1:nxtrack,sline1:eline1), errstat)
    call tiof_pop_group (tio_l1obj, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_geo_data: failed to read geolocation data", &
           errstat)
      return
    endif

    ! FIXME - note that the geometric calculations below probably need
    ! to be re-thought for TEMPO. For now we keep them for testing

    ! Add small offset to SZA if it is zero, to avoid VLIDORT failure
    do i = 1, nxtrack
      do j = sline1,eline1
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

    !------------------------------------------
    ! dimension after coadding
    !-----------------------------------------
    nl = nl/nybin
    nx = nxtrack/nxbin
    nbin = nxbin*nybin

    geo%sza = -999
    geo%gflg (:,:) = 0
    
    Do iy = 0, nl-1
       ysidx = (sline -1)+ iy * nybin 
       yeidx = ysidx + nybin - 1
       ymidx = ysidx + nybin / 2

       geo%time(iy) = sum(tio_time(ysidx:yeidx)) / nybin
       DO ix = 1, nx
         iix = (ix-1)*nxbin
         xsidx = (ix - 1) * nxbin + 1
         xeidx = xsidx + nxbin - 1
         xmidx = xsidx + nxbin / 2
         IF (ALL( (tio_clat(:,xsidx:xeidx,ysidx:yeidx) < 90 .and. tio_clat(:,xsidx:xeidx, ysidx:yeidx) > -90)  )) THEN 

             sza = sum(tio_sza(xsidx:xeidx, ysidx:yeidx))/nbin
             vza = sum(tio_vza(xsidx:xeidx, ysidx:yeidx))/nbin
             IF (sza .lt. 1e-5) sza = sza + 1e-5
             IF (vza .lt. 1e-5) vza = vza + 1e-5
             sazm = sum(tio_saza(xsidx:xeidx, ysidx:yeidx))/nbin
             vazm = sum(tio_vaza(xsidx:xeidx, ysidx:yeidx))/nbin
             relaza = ABS( vazm - sazm) 
             IF ( relaza < 180) relaza = 360.0 - relaza 
           
             geo%sza(ix, iy) = sza
             geo%vza(ix, iy) = vza
             geo%aza(ix, iy) = relaza
             geo%sca(ix, iy) = 0.0 
                
             geo%gflg(ix, iy)  = tio_geoflg(xmidx, ymidx)
             geo%height(ix, iy)= int( sum(1.0 * tio_height(xsidx:xeidx, ysidx:yeidx)) &
                                 / (1.0 * nbin), kind=2)
             geo%lon(ix, iy) = sum(tio_lon(xsidx:xeidx, ysidx:yeidx))/nbin
             geo%lat(ix, iy) = sum(tio_lat(xsidx:xeidx, ysidx:yeidx))/nbin
             geo%clon(1,ix, iy) = tio_clon(1,xeidx, yeidx)
             geo%clon(2,ix, iy) = tio_clon(2,xsidx, yeidx)
             geo%clon(3,ix, iy) = tio_clon(3,xeidx, ysidx)
             geo%clon(4,ix, iy) = tio_clon(4,xsidx, ysidx)

             geo%clat(1,ix, iy) = tio_clat(1,xeidx, yeidx)
             geo%clat(2,ix, iy) = tio_clat(2,xsidx, yeidx)
             geo%clat(3,ix, iy) = tio_clat(3,xeidx, ysidx)
             geo%clat(4,ix, iy) = tio_clat(4,xsidx, ysidx)              
             geo%elon(ix-1:ix, iy) = geo%clon(1:2,ix,iy)
             geo%elat(ix-1:ix, iy) = geo%clat(1:2,ix,iy)
             !print * ,xsidx,ix,geo%sza(ix, iy), geo%vza(ix, iy)
             
         ENDIF
      ENDDO
      call convert_gpqualflag_info (nx,geo%gflg(1:nx, iy), &
               land_water_flg(1:nx), glint_flg(1:nx),snow_ice_flg(1:nx))
      geo%glint_flg(1:nx,iy) = glint_flg(1:nx)
      geo%snow_ice_flg(1:nx,iy) = snow_ice_flg(1:nx)
      geo%land_water_flg(1:nx,iy) = land_water_flg(1:nx)
    ENDDO
      
    !IF (do_deallo) THEN 
    deallocate(tio_clon, tio_clat, tio_lon, tio_lat, tio_sza, tio_saza, & 
               tio_vza,tio_vaza, tio_height, tio_geoflg, tio_time)
    !ENDIF
    RETURN
 
  end subroutine read_geo_tio

 SUBROUTINE convert_gpqualflag_info ( &
       nxtrack, geoflg, land_water_flg, glint_flg, snow_ice_flg )
    
    USE OMSAO_precision_module
    USE m_convert_coadd, ONLY: convert_2bytes_to_32bits
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
    INTEGER (KIND=i4), DIMENSION (nxtrack), INTENT (IN) :: geoflg

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: land_water_flg, glint_flg,snow_ice_flg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4),                PARAMETER      :: nbyte = 32
    INTEGER (KIND=i4), DIMENSION (7), PARAMETER      :: seven_byte = int((/ 1,2, 4, 8, 16, 32, 64 /), kind=i4)
    INTEGER (KIND=i4)                                :: i
    INTEGER (KIND=i4), DIMENSION (nxtrack)           :: tmp_flg
    INTEGER (KIND=i4), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    land_water_flg = 0
    glint_flg = 0
    snow_ice_flg = 0

    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    tmp_flg(1:nxtrack) = geoflg(1:nxtrack) ;  tmp_bytes = 0
    ! FIXME - TEMPO ground_pixel_flag is now int4, but all subroutines
    ! in this module assume int2

    CALL convert_2bytes_to_32bits ( &
         nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1) )
    ! ------------------------------
    ! The Glint flag is easy: Byte 4
    ! ------------------------------
    glint_flg(1:nxtrack) = tmp_bytes(1:nxtrack,4)

    ! ------------------------------------------------------------------
    ! Land/Water and Ice require a bit more work. The BIT slices must be
    ! multiplied with the corresponding powers of 2. The sum over this
    ! product is the information we seek.
    ! ------------------------------------------------------------------
    DO i = 1, nxtrack
      land_water_flg(i) = SUM(tmp_bytes(i,0:3 )*seven_byte(1:4))
      snow_ice_flg  (i) = SUM(tmp_bytes(i,8:14)*seven_byte(1:7))
    END DO

    RETURN
  END SUBROUTINE convert_gpqualflag_info


 ! FIXME - geometry probably needs checking for TEMPO
  !
  !> Calculate PARAMETER values at binned pixel corners
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

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (kind=4), INTENT(in)    :: nxtrack, ntimes
    real (kind=8), dimension (1:nxtrack, 0:ntimes-1), INTENT(inout) :: lon, &
         lat, sza, saza, vza, vaza
    real (kind=8), dimension (0:nxtrack, 0:ntimes), INTENT(out) :: clon, clat
    real (kind=8), dimension (0:nxtrack, 0:ntimes-1), INTENT(out) :: elon, elat
    real (kind=8), dimension (1:nxtrack, 0:ntimes-1), INTENT(out) :: esza, &
         evza, eaza, esca

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (kind=4) :: i, j, jj, ix, mpix, nx, ny
    INTEGER                                             :: errstat
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
      IF (errstat /= 0) THEN 
         print * , 'get_corners: errors in interpol'
         STOP
      ENDIF
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
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    real (kind=8),    INTENT (IN) :: lat1, lat2, lon1, lon2, c0, c

    ! ----------------
    ! Output variables
    ! ----------------
    real (kind=8),  INTENT (OUT)  :: lat, lon

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

    IMPLICIT NONE

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    real (kind=8), INTENT (IN) :: lat1, lon1, lat2, lon2
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

  REAL (kind=8) FUNCTION angle_minus_twopi ( gamma0, pival ) result ( gamma )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (kind=8), INTENT (IN) :: gamma0, pival

    if ( gamma0 > pival ) then
      gamma = gamma0 - 2.0d0 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    else if ( gamma0 < -pival ) then
      gamma = gamma0 + 2.0d0 * pival
    else
      gamma = gamma0
    end if

    return
  END FUNCTION angle_minus_twopi

end module m_read_geo_tio
