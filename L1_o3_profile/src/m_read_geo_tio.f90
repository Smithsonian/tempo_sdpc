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

  PUBLIC read_geo_tio
  PRIVATE convert_gpqualflag_info, scattering_angle_deg, edge_midpoint

CONTAINS

  function scattering_angle_deg (sza_deg, vza_deg, relaza_deg) result (sca_deg)
    implicit none
    real (kind=8), intent(in) :: sza_deg, vza_deg, relaza_deg
    real (kind=8) :: sca, sza, vza, relaza, sca_deg

    sza = sza_deg * deg2rad
    vza = vza_deg * deg2rad
    relaza = relaza_deg * deg2rad

    sca = acos(cos(sza)*cos(vza) + sin(sza)*sin(vza)*cos(relaza))
    sca_deg = 180.d0 - sca * rad2deg

  end function scattering_angle_deg

  subroutine edge_midpoint (clon_f, clat_f, elon_f, elat_f, sgn)
    use OMSAO_pixelcorner_module, ONLY:  sphergeom_intermediate, circle_rdis
    implicit none
    real (kind=4), intent(in) :: clon_f(2), clat_f(2)
    real (kind=4), intent(out) :: elon_f, elat_f
    integer, intent(in) :: sgn
    ! local variables
    real (kind=8) :: clon(2), clat(2), elon, elat, dist_gc

    clon(:) = real(clon_f * deg2rad, kind=8)
    clat(:) = real(clat_f * deg2rad, kind=8)

    ! great circle distance
    dist_gc = circle_rdis (clat(1),clon(1), clat(2),clon(2))

    ! edge mid-point
    call sphergeom_intermediate (clat(1),clon(1), clat(2),clon(2), dist_gc, &
                                 0.5*dist_gc*sign(1,sgn), elat, elon)

    elon_f = real(elon * rad2deg, kind=4)
    elat_f = real(elat * rad2deg, kind=4)

  end subroutine edge_midpoint

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
                           spix, lpix, sline, eline, errstat)
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: r8_missval
    USE OMSAO_variables_module, ONLY: geo_group, nxbin, nybin, &
                                      l1file=>l1b_rad_filename
    use netcdf, only : nf90_noerr, nf90_inq_varid

    IMPLICIT NONE

    ! input variables
    INTEGER, INTENT (in) :: ntimes, nxtrack
    INTEGER, INTENT(inout) ::  spix, lpix, sline, eline
    CHARACTER (len=*), INTENT(in) :: l1swath
    ! output variables
    TYPE (geo_group), INTENT(INOUT) :: geo
    INTEGER, INTENT (out) :: errstat
    ! local variables
    type (tiof_file_type) :: tio_l1obj
    INTEGER :: nline, nx, i, j, ix, iy,iix, nl, sline1, eline1, nbin
    INTEGER :: ysidx, yeidx, ymidx, xsidx, xeidx, xmidx, varid, status, num_sf
    !-----------------------------------------------------------
    ! variables for reading original tempo geolocation dataset
    !------------------------------------------------------------
    real (kind=8), dimension (:), allocatable    :: tio_time
    real (kind=4), dimension (:,:), allocatable :: tio_lat, &
                    tio_lon, tio_sza, tio_vza, tio_saza, tio_vaza
    INTEGER (kind=2), dimension (:,:), allocatable :: tio_height
    INTEGER (kind=4), dimension (:,:), allocatable :: tio_geoflg
    real (kind=4), dimension (:,:), allocatable:: tio_snowice_fraction
    real (kind=4), dimension (:,:,:), allocatable:: tio_clon, tio_clat
    integer (kind=4), dimension(:), allocatable :: step_idx
    !---------------------------------------------------------
    ! variables for binning/arrange tempo geolocation dataset
    !--------------------------------------------------------
    real (kind=8) :: sza, vza,sazm, vazm, relaza
    real (kind=4) :: sf_array(nxbin,nybin)
    logical :: sf_mask (nxbin,nybin)

    allocate(tio_clon (4, nxtrack, 0:ntimes-1) , tio_clat(4, nxtrack, 0:ntimes-1) )
    allocate(tio_lon ( nxtrack, 0:ntimes-1) , tio_lat(nxtrack, 0:ntimes-1))
    allocate(tio_sza ( nxtrack, 0:ntimes-1) , tio_saza(nxtrack, 0:ntimes-1))
    allocate(tio_vza ( nxtrack, 0:ntimes-1) , tio_vaza(nxtrack, 0:ntimes-1))
    allocate(tio_height ( nxtrack, 0:ntimes-1) , tio_geoflg(nxtrack, 0:ntimes-1))
    allocate(tio_time (0:ntimes-1), step_idx(0:ntimes-1))
    geo % have_snowice_fraction = .false.

    !------------------------------------------------------------------------
    ! Initialize
    !------------------------------------------------------------------------

    errstat = 0
    tio_sza = -999
    ! limits of region to be read in
    nl = eline - sline +1
    sline1 = sline
    eline1 = eline
    !if (eline == sline) eline1 = sline1 + 1 !xl, 12/3/2021
    IF (eline1 > ntimes) THEN
        eline1 = eline1-1 ; sline1=sline1-1
    ENDIF
    sline1 = sline1-1
    eline1 = eline1-1
    nline = eline1 - sline1 + 1
    ! read in geolocation data for the chosen subset of step positions
    call tiof_open (l1file, tio_l1obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (tio_l1obj, o3p_var_time, [sline1], [nline], &
                        tio_time(sline1:eline1), errstat, replace_fill=r8_missval)
    call tiof_get1d_i4 (tio_l1obj, tempo_dim_step, [sline1], [nline], &
                        step_idx(sline1:eline1), errstat)
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
    status = nf90_inq_varid (tio_l1obj % groupid, tempo_var_snowice_fraction, varid)
    if (status == nf90_noerr) then
      allocate(tio_snowice_fraction( nxtrack, 0:ntimes-1))
      call tiof_get2d_r4 (tio_l1obj, tempo_var_snowice_fraction, [sline1,0], [nline, nxtrack], &
                          tio_snowice_fraction(1:nxtrack,sline1:eline1), errstat)
      geo % have_snowice_fraction = .true.
    endif
    call tiof_pop_group (tio_l1obj, errstat)
    call tiof_close (tio_l1obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_geo_data: failed to read geolocation data", &
           errstat)
      return
    endif

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
    !print *, nl-1, sline1, eline1, nybin, sline, eline
    geo%step_idx(0:nl-1) = step_idx (sline1:eline1:nybin)

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
             relaza = ABS( vazm - (180.0 + sazm))
             IF ( relaza > 360) relaza = relaza - 360
             IF ( relaza > 180) relaza = 360.0 - relaza

             geo%sza(ix, iy) = real(sza, kind=r4)
             geo%vza(ix, iy) = real(vza, kind=r4)
             geo%aza(ix, iy) = real(relaza, kind=r4)
             geo%sca(ix, iy) = real(scattering_angle_deg (sza, vza, relaza), kind=r4)

             ! FIXME: For now, the center pixel determines the binned pixel ground pixel status flag
             !        How should status flag binning work?
             geo%gflg(ix, iy)  = tio_geoflg(xmidx, ymidx)

             if (geo % have_snowice_fraction) then
               ! JCH: Since snow_ice_flg is mostly used like a measure of fractional snow/ice
               ! cover, I'm assuming it's ok to replace it with something that actual measures
               ! the fractional area covered by snow/ice.
               sf_array = tio_snowice_fraction(xsidx:xeidx, ysidx:yeidx)
               sf_mask = 0.0 .le. sf_array .and. sf_array .le. 1.0   ! exclude missing data
               num_sf = count(sf_mask)
               if (num_sf > 0) then
                 geo%snow_ice_flg(ix,iy) = aint(100.0*sum(sf_array, mask=sf_mask)/num_sf)
               else
                 geo%snow_ice_flg(ix,iy) = 0
               endif
             endif

             geo%height(ix, iy)= int( sum(1.0 * tio_height(xsidx:xeidx, ysidx:yeidx)) &
                                 / (1.0 * nbin), kind=2)
             geo%lon(ix, iy) = sum(tio_lon(xsidx:xeidx, ysidx:yeidx))/nbin
             geo%lat(ix, iy) = sum(tio_lat(xsidx:xeidx, ysidx:yeidx))/nbin

             ! ix increases southward, iy increases westward
             geo%clon(1,ix, iy) = tio_clon(1,xsidx, ysidx) ! NE
             geo%clon(2,ix, iy) = tio_clon(2,xsidx, yeidx) ! NW
             geo%clon(3,ix, iy) = tio_clon(3,xeidx, yeidx) ! SW
             geo%clon(4,ix, iy) = tio_clon(4,xeidx, ysidx) ! SE

             geo%clat(1,ix, iy) = tio_clat(1,xsidx, ysidx)
             geo%clat(2,ix, iy) = tio_clat(2,xsidx, yeidx)
             geo%clat(3,ix, iy) = tio_clat(3,xeidx, yeidx)
             geo%clat(4,ix, iy) = tio_clat(4,xeidx, ysidx)

             ! file stores pixel corners in order: (NE,NW,SW,SE)
             call edge_midpoint (geo%clon(1:2,ix,iy), geo%clat(1:2, ix,iy), &
                                 geo%elon(ix,iy), geo%elat(ix,iy), 1)
             if (ix == 1) then
               call edge_midpoint (geo%clon(1:2,ix,iy), geo%clat(1:2, ix,iy), &
                                   geo%elon(ix-1,iy), geo%elat(ix-1,iy), -1)
             endif
         ENDIF
      ENDDO
      call convert_gpqualflag_info (nx, geo%gflg(1:nx, iy), &
                                    geo%land_water_flg(1:nx, iy), &
                                    geo%glint_flg(1:nx, iy), &
                                    geo%snow_ice_flg(1:nx, iy), &
                                    geo%have_snowice_fraction)
    ENDDO

    deallocate(tio_clon, tio_clat, tio_lon, tio_lat, tio_sza, tio_saza, &
               tio_vza,tio_vaza, tio_height, tio_geoflg, tio_time)
    if (geo%have_snowice_fraction) then
      deallocate (tio_snowice_fraction)
    endif
    RETURN

  end subroutine read_geo_tio

 SUBROUTINE convert_gpqualflag_info (nxtrack, geoflg, land_water_flg, glint_flg, &
       snow_ice_flg, have_snowice_fraction )

    USE OMSAO_precision_module
    IMPLICIT NONE
    INTEGER (KIND=i4), INTENT (IN) :: nxtrack
    INTEGER (KIND=i4), DIMENSION (:), INTENT (IN) :: geoflg
    INTEGER (KIND=i2), DIMENSION (:), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg
    logical, intent(in) :: have_snowice_fraction
    INTEGER (KIND=i4) :: i

    ! bit 0 is least signifcant, n-1 is most signifcant

    do i = 1, nxtrack
      ! MODIS land-water mask is in bits 0-3
      land_water_flg(i) = int(ibits(geoflg(i), 0, 4), kind=i2)
      ! NISE snow-ice mask is in bits 8-15
      if (.not. have_snowice_fraction) then
        snow_ice_flg(i) = int(ibits(geoflg(i), 8, 8), kind=i2)
      endif
      ! glint possibility is bit 4
      glint_flg(i) = int(ibits(geoflg(i), 4, 1), kind=i2)
    enddo

    RETURN
  END SUBROUTINE convert_gpqualflag_info

end module m_read_geo_tio
