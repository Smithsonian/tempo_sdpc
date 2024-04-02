MODULE omi_pge_postprocessing
  use tell_module
  implicit none
  private
  public omi_pge_postprocess
CONTAINS

SUBROUTINE omi_pge_postprocess ( &
    l1bfile, omi_radiance_swathname, pge_idx, &
    ntimes, nxtrack, do_process_line, xtrange, is_szoom, n_max_rspec, &
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
  use ctrlvars, only: yn_radiance_reference, yn_refseccor, yn_do_he5_output, &
       yn_gems
  USE OMSAO_indices_module, ONLY: pge_hcho_idx
  USE OMSAO_Reference_sector_module, ONLY: reference_sector_correction
  USE OMSAO_wfamf_module, ONLY: amf_calculation, &
    wfamf_deallocate
  USE he5_output_tools, ONLY: saopge_geofield_read, saopge_columninfo_read, &
    he5_write_fitting_statistics, saopge_geofieldtime_read
  use output_tools, only : read_geofields, read_column_results, &
    copy_pixel_corners, copy_metadata, copy_gpqf_attributes
  USE omi_read_l1b_data, ONLY: omi_read_glint_ice_land_flags
  USE omi_pge_fitting_aux, ONLY: compute_fitting_statistics, fitting_statistics_type
  USE OMSAO_variables_module, ONLY: max_good_col, l1b_rad_filename
  use m_read_gems, only: gems_read_ice_glint, gems_read_geofields
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),                              INTENT (IN) :: l1bfile, omi_radiance_swathname
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
  REAL    (KIND=r8), DIMENSION (0:ntimes-1) :: time
  REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1) :: lat, lon, sza, vza, &
       saa, vaa, thg
  REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: saocol, saodco, saorms, saoamf
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: saofcf, amfdiag
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: glint_flg, snow_ice_flg, land_water_flg
  LOGICAL                                             :: do_write
  logical :: corners_copied

  ! --------------
  ! Error handling
  ! --------------
  INTEGER (KIND=i4) :: locerrstat

  if (errstat /= 0) return

  ! -------------------------
  ! Initialize error variable
  ! -------------------------
  locerrstat = 0

  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  if (.not. yn_gems) then !TEMPO
    call read_geofields (l1b_rad_filename, ntimes, nxtrack, lat, lon, &
         sza, vza, saa, vaa, thg, time, errstat)
  else !GEMS
    call gems_read_geofields (l1b_rad_filename, ntimes, nxtrack, lat, lon, &
         sza, vza, saa, vaa, thg, time, errstat)
  endif
  if (errstat /= 0) return

  if (.not. yn_gems) then !TEMPO
    call copy_metadata (l1bfile, errstat)
    if (errstat /= 0) return

    call copy_gpqf_attributes (l1bfile, omi_radiance_swathname, errstat)
    if (errstat /= 0) return

    call copy_pixel_corners (l1bfile, omi_radiance_swathname, &
         ntimes, nXtrack, corners_copied, errstat)
    if (errstat /= 0) return
  else!GEMS
    corners_copied = .false.
  endif
  if (.not.corners_copied) then
    CALL compute_pixel_corners ( ntimes, nXtrack, lat, lon, is_szoom, errstat)
    if (errstat /= 0) return
  endif

  call read_column_results (ntimes, nxtrack, saocol, saodco, saorms, &
       saoamf, saofcf, errstat)
  if (errstat /= 0) return

  ! ----------------------------------
  ! Read L1b glint and snow/ice flags
  ! ----------------------------------
  if (.not. yn_gems) then !TEMPO
    CALL omi_read_glint_ice_land_flags ( &
         l1bfile, nxtrack, ntimes, snow_ice_flg, glint_flg, land_water_flg, errstat )
  else !GEMS
    call gems_read_ice_glint (l1bfile, nxtrack, ntimes, snow_ice_flg, &
         glint_flg, errstat)
  endif
  if (errstat /= 0) return

  ! -----------
  ! Compute AMF
  ! -----------
  do_write = .TRUE.
  call tell_log (1, 'omi_pge_postprocess:  calling amf_calculation')
  CALL amf_calculation (                             &
    pge_idx, ntimes, nxtrack, lat, lon, sza, vza, saa, vaa, time,  &
    snow_ice_flg, glint_flg, land_water_flg, xtrange,       &
    saocol, saodco, saoamf, amfdiag, thg, do_write, &
    errstat              )
  if (errstat /= 0) return

  ! ----------------------------------
  ! Compute average fitting statistics
  ! ----------------------------------
  allocate (fit_stats%quality_flag (nxtrack, 0:ntimes-1), stat=errstat)
  if (errstat /= 0) then
    call tell_error (tell_malloc_error, &
         "omi_pge_postprocess:  allocate failed", errstat)
    return
  endif
  call tell_log (1, 'omi_pge_postprocess:  calling compute_fitting_statistics')
  CALL compute_fitting_statistics ( ntimes, nxtrack, xtrange, &
    saocol, saodco, saorms, saofcf, amfdiag, fit_stats, errstat)
  if (errstat /= 0) return
  if (yn_do_he5_output) then
    CALL he5_write_fitting_statistics ( &
      pge_idx, max_good_col, nxtrack, ntimes, fit_stats % quality_flag, &
      fit_stats%col_avg, fit_stats%dcol_avg, fit_stats%rms_avg, locerrstat)
  endif
  ! ---------------------------------------
  ! Apply cross-track destriping correction
  ! ---------------------------------------
  call tell_log (1, 'omi_pge_postprocess:  calling xtrack_destriping')
  CALL xtrack_destriping (ntimes, nxtrack, do_process_line, xtrange, &
    lat, saocol, &
    fit_stats % quality_flag, errstat)

  ! ---------------------------------------------------------------
  ! Apply Reference Sector Correction; Only for HCHO retrieval !gga
  ! ---------------------------------------------------------------
  IF ((yn_refseccor) .AND. ( pge_idx == pge_hcho_idx ) .AND.  &
    (yn_radiance_reference)) THEN
    call tell_log (1, 'omi_pge_postprocess:  calling Reference_Sector_correction')
    CALL Reference_Sector_correction (ntimes, nxtrack, &
      saocol, saodco, saoamf, fit_stats % quality_flag, pge_idx, n_max_rspec, &
      errstat)
  ENDIF

  ! Deallocate AMF variables
  CALL wfamf_deallocate (errstat)
  if (errstat /= 0) then
    call tell_error(tell_runtime_error, "wfamf_deallocate failed", errstat)
    return
  endif

  RETURN
END SUBROUTINE omi_pge_postprocess
END MODULE

