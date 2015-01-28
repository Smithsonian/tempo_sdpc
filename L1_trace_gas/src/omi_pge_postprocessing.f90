MODULE omi_pge_postprocessing
  use tell_module
  implicit none
  private
  public omi_pge_postprocess
CONTAINS

SUBROUTINE omi_pge_postprocess ( &
    l1bfile, pge_idx, ntimes, nxtrack, do_process_line, xtrange, is_szoom, n_max_rspec, &
    fit_stats, errstat )

  ! ---------------------------------------------------------
  ! In this subroutine we collect all those computations that
  ! are done "post fitting". These include:
  !
  ! (1) AMF calculation
  ! (2) Fitting statistics
  ! (3) Cross-track destriping
  ! (4) Ground-pixel corner computation
  ! (5) Reference Sector Background Correction for HCHO
  ! ---------------------------------------------------------

  USE OMSAO_precision_module
  USE OMSAO_pixelcorner_module, ONLY: compute_pixel_corners
  USE OMSAO_destriping_module, ONLY: xtrack_destriping
  USE OMSAO_errstat_module
  use ctrlvars, only: yn_radiance_reference, yn_refseccor
  USE OMSAO_indices_module, ONLY: pge_hcho_idx
  USE OMSAO_Reference_sector_module, ONLY: reference_sector_correction
  USE OMSAO_wfamf_module, ONLY: amf_calculation_bis, climatology_allocate, Cmlat, Cmlon, CmETA, CmEp1
  USE he5_output_tools, ONLY: saopge_geofield_read, saopge_columninfo_read, &
    he5_write_fitting_statistics
  use output_tools, only : read_geofields, read_column_results
  USE omi_read_l1b_data, ONLY: omi_read_glint_ice_flags
  USE omi_pge_fitting_aux, ONLY: compute_fitting_statistics, fitting_statistics_type
  USE OMSAO_variables_module, ONLY: max_good_col
  use datafields, only: lat_field, lon_field, sza_field, thgt_field, vza_field
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),                              INTENT (IN) :: l1bfile
  INTEGER (KIND=i4),                              INTENT (IN) :: ntimes, nxtrack, n_max_rspec, pge_idx
  INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),  INTENT (IN) :: xtrange
  LOGICAL,           DIMENSION (0:ntimes-1),      INTENT (IN) :: do_process_line, is_szoom

  ! -----------------
  ! Modified variable
  ! -----------------
  type (fitting_statistics_type), intent(inout) :: fit_stats
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ----------------
  ! Local variables
  ! ----------------
  ! (1) OMI data
  ! ----------------
  REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1) :: lat, lon, sza, vza, thg
  REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: saocol, saodco, saorms, saoamf
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: saofcf !, saomqf
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: glint_flg, snow_ice_flg
  LOGICAL                                             :: do_write

  ! --------------
  ! Error handling
  ! --------------
  INTEGER (KIND=i4) :: locerrstat

  if (errstat < 0) return

  ! -------------------------
  ! Initialize error variable
  ! -------------------------
  locerrstat = pge_errstat_ok

  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  if (.false.) then
    ! FIXME (to be removed)
  CALL  saopge_geofield_read ( ntimes, nxtrack, lat_field,  lat, locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, lon_field,  lon, locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, sza_field,  sza, locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, vza_field,  vza, locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, thgt_field, thg, locerrstat )
  else
    call read_geofields (ntimes, nxtrack, lat, lon, sza, vza, thg, errstat)
    if (errstat < 0) return
  endif

  ! ----------------------------------------------------
  ! Compute ground pixel corner latitudes and longitudes
  ! ----------------------------------------------------
  CALL compute_pixel_corners ( ntimes, nXtrack, lat, lon, is_szoom, locerrstat )

  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  if (.false.) then
  CALL saopge_columninfo_read (                 &  ! FIXME (<--- to be removed)
    ntimes, nxtrack, saocol, saodco, saorms, &
    saoamf, saofcf, locerrstat                 )
  else
    call read_column_results (ntimes, nxtrack, saocol, saodco, saorms, &
                              saoamf, saofcf, errstat)
    if (errstat < 0) return
  endif

  ! ----------------------------------
  ! Read L1b glint and snow/ice flags
  ! ----------------------------------
  CALL omi_read_glint_ice_flags ( &
    l1bfile, nxtrack, ntimes, snow_ice_flg, glint_flg, errstat )
  if (errstat < 0) return

  ! -----------
  ! Compute AMF
  ! -----------
  !!$  CALL amf_calculation (                             &
  !!$       pge_idx, ntimes, nxtrack, lat, lon, sza, vza, &
  !!$       snow_ice_flg, glint_flg, xtrange, is_szoom,   &
  !!$       saocol, saodco, saoamf, locerrstat              )

  ! ----------------
  ! Comnpute AMF bis
  ! ----------------
  do_write = .TRUE.
  call tell_log (1, 'omi_pge_postprocess:  calling amf_calculation_bis')
  CALL amf_calculation_bis (                             &
    pge_idx, ntimes, nxtrack, lat, lon, sza, vza,     &
    snow_ice_flg, glint_flg, xtrange, is_szoom,       &
    saocol, saodco, saoamf, thg, do_write, &
    locerrstat              )

  ! ----------------------------------
  ! Compute average fitting statistics
  ! ----------------------------------
  allocate (fit_stats%quality_flag (nxtrack, 0:ntimes-1), stat=errstat)
  if (errstat /= 0) then
    call tell_error (tell_malloc_error, "omi_pge_postprocess:  allocate failed", errstat)
    return
  endif
  call tell_log (1, 'omi_pge_postprocess:  calling compute_fitting_statistics')
  CALL compute_fitting_statistics ( &
    pge_idx, ntimes, nxtrack, xtrange, &
    saocol, saodco, saorms, saofcf, fit_stats, locerrstat )
  CALL he5_write_fitting_statistics ( &
    pge_idx, max_good_col, nxtrack, ntimes, fit_stats % quality_flag, &
    fit_stats%col_avg, fit_stats%dcol_avg, fit_stats%rms_avg, locerrstat)
  errstat = max(locerrstat, errstat)
  ! ---------------------------------------
  ! Apply cross-track destriping correction
  ! ---------------------------------------
  call tell_log (1, 'omi_pge_postprocess:  calling xtrack_destriping')
  CALL xtrack_destriping (                                    &
    pge_idx, ntimes, nxtrack, do_process_line, xtrange,         &
    lat, saocol, & !saodco, saoamf, saofcf,
    fit_stats % quality_flag, locerrstat )

  ! ---------------------------------------------------------------
  ! Apply Reference Sector Correction; Only for HCHO retrieval !gga
  ! ---------------------------------------------------------------
  IF ((yn_refseccor) .AND. ( pge_idx == pge_hcho_idx ) .AND.  &
    (yn_radiance_reference)) THEN
    call tell_log (1, 'omi_pge_postprocess:  calling Reference_Sector_correction')
    CALL Reference_Sector_correction (ntimes, nxtrack, & !xtrange, lat,
      saocol, saodco, saoamf, fit_stats % quality_flag, pge_idx, n_max_rspec, &
      locerrstat)
  ENDIF

  ! --------------------------------
  ! Deallocate Climatology variables
  ! --------------------------------
  CALL climatology_allocate ( "d", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE omi_pge_postprocess
END MODULE

