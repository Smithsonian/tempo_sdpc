MODULE OMSAO_omidata_module

  USE OMSAO_precision_module, ONLY: r4, r8, i4, i2, i1
  USE OMSAO_parameters_module, ONLY: &
    MAX_STR_LEN, max_spec_pts, nxtrack_max, nlines_max, nutcdim, nwavel_max
  USE OMSAO_indices_module, ONLY: n_max_fitpars, max_rs_idx, max_calfit_idx, &
    o3_t1_idx, o3_t3_idx
  IMPLICIT NONE

  type retrieval_type
    ! 2d arrays are (nxtrack, 0:ntimes-1)
    real (kind=r8), dimension(:,:), allocatable :: column_amount, column_uncertainty
    real (kind=r8), dimension(:,:), allocatable :: rms
    real (kind=r8), dimension(:), allocatable :: time
    real (kind=r4), dimension(:,:), allocatable :: latitude, longitude, height
    real (kind=r4), dimension(:,:), allocatable :: sza, vza, saa, vaa
    integer (kind=i2), dimension(:,:), allocatable :: fit_flag, xtr_flag
    integer (kind=i4) :: nxtrack, ntimes
  end type retrieval_type

  !> Input variables
  ! FIXME: (JCH)
  ! I've defined instances of input_vars_type and result_vars_type as a
  ! crutch to help gather related module variables into structures without
  ! making these changes everywhere at once.  New code can access needed
  ! data through these structures, while old code can continue to work
  ! unmodified, until it also gets updated to use these structures. Once
  ! these structures are used everywhere, we can greatly reduce the
  ! number of exported module symbols.
  type, public :: input_vars_type
    real (kind=r8), dimension(:), pointer :: time => null()
    real (kind=r4), dimension(:,:), pointer :: latitude => null()
    real (kind=r4), dimension(:,:), pointer :: longitude => null()
    real (kind=r4), dimension(:,:), pointer :: solar_zenith => null()
    real (kind=r4), dimension(:,:), pointer :: solar_azimuth => null()
    real (kind=r4), dimension(:,:), pointer :: viewing_zenith => null()
    real (kind=r4), dimension(:,:), pointer :: viewing_azimuth => null()
    real (kind=r4), dimension(:,:), pointer :: snow_ice_fraction => null()
    integer (kind=i2), dimension(:,:), pointer :: terrain_height => null()
    integer (kind=i4), dimension(:,:), pointer :: ground_pixel_quality_flag => null()
  end type input_vars_type

  !> Radiance fit results
  type, public :: result_vars_type
    real (kind=r8), dimension(:,:), pointer :: column_amount => null()
    real (kind=r8), dimension(:,:), pointer :: column_uncert => null()
    real (kind=r8), dimension(:,:), pointer :: fit_rms_residual => null()
    integer(kind=i2), dimension(:,:), pointer :: fit_convergence_flag => null()
    integer (kind=i2), dimension(:,:), pointer :: fit_iteration_count => null()
    integer (kind=i2), dimension(:), pointer :: solcal_convergence_flag => null()
    real (kind=r8), dimension(:), pointer :: solcal_shift => null()
    integer (kind=i2), dimension(:), pointer :: radcal_convergence_flag => null()
    integer (kind=i2), dimension(:), pointer :: radref_convergence_flag => null()
    real (kind=r8), dimension(:), pointer :: radref_column_amount => null()
    real (kind=r8), dimension(:), pointer :: radref_column_uncert => null()
    real (kind=r8), dimension(:), pointer :: radref_column_xtrfit => null()
    real (kind=r8), dimension(:), pointer :: radref_fit_rms => null()
  end type result_vars_type

  type (input_vars_type), save :: input_vars
  type (result_vars_type), save :: result_vars

  !> Radiance fit diagnostics
  type, public :: radfit_diagnostics_type
    ! these arrays are dimension(n_fitvar_rad,num_xtrack,0:nblock-1):
    real (kind=r8), dimension (:,:,:), pointer :: params => null()
    real (kind=r8), dimension (:,:,:), pointer :: errors => null()
    real (kind=r8), dimension (:,:,:), pointer :: correl => null()
    ! this array is dimension(n_rad_wvl,num_xtrack,4,0:nblock-1):
    real (kind=r8), dimension (:,:,:,:), pointer :: fitspc => null()
  end type radfit_diagnostics_type

  type, public :: amf_correction_type
    ! 2d arrays have dimension (1:nxtrack,0:ntimes-1)
    real (kind=r8), dimension (:,:), pointer :: amf_molecule_specific => null()
    real (kind=r8), dimension (:,:), pointer :: amf_molecule_stratospheric => null()
    real (kind=r8), dimension (:,:), pointer :: amf_molecule_tropospheric => null()
    real (kind=r8), dimension (:,:), pointer :: amf_molecule_specific_clear_sky => null()
    real (kind=r8), dimension (:,:), pointer :: amf_molecule_stratospheric_clear_sky => null()
    real (kind=r8), dimension (:,:), pointer :: amf_molecule_tropospheric_clear_sky => null()
    real (kind=r8), dimension (:,:), pointer :: amf_geometric => null()
    real (kind=r8), dimension (:,:), pointer :: cloud_fraction => null()
    real (kind=r8), dimension (:,:), pointer :: cloud_pressure => null()
    integer (kind=i2), dimension (:,:), pointer :: diagnostic_flag => null()
    real (kind=r4), dimension (:,:), pointer :: surface_pressure => null() 
    real (kind=r4), dimension (:,:), pointer :: tropopause_pressure => null()
    real (kind=r4), dimension (:), pointer :: eta_a => null()
    real (kind=r4), dimension (:), pointer :: eta_b => null()
  end type amf_correction_type

  PRIVATE MAX_STR_LEN, max_spec_pts, nxtrack_max, nlines_max, nutcdim, nwavel_max
  PRIVATE r4, r8, i4, i2, i1
  PRIVATE n_max_fitpars, max_rs_idx, max_calfit_idx, &
    o3_t1_idx, o3_t3_idx

  ! ------------------------------------------------------------
  ! Boundary wavelengths (approximate) for UV-2 and VIS channels
  ! ------------------------------------------------------------
  REAL (KIND=r4), PARAMETER :: uv2_upper_wvl = 385.0_r4, vis_lower_wvl = 350.0_r4

  ! ---------------------------------------
  ! Minimum OMI spectral resolution (in nm)
  ! ---------------------------------------
  REAL (KIND=r8), PARAMETER :: omi_min_specres = 0.5_r8

  ! -------------------------------------------------------
  ! Maximum dimension of OMI data fields in SAO L2 products
  ! -------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_field_maxdim = 3

  ! --------------------------------------------------
  ! Parameters defined by the NISE snow cover approach
  ! --------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
    NISE_snowfree =   0, NISE_allsnow = 100, NISE_permice = 101, NISE_drysnow = 103, &
    NISE_ocean    = 104, NISE_suspect = 125, NISE_error   = 127

  ! --------------------------------------------------------------------
  ! Values to go into the diagnostic array that shows how the AMF was
  ! computed, with values that indicate missing cloud products, glint,
  ! and geometric or no AMF
  ! --------------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
    omi_cld_addmiss = 200, omi_glint_add = 10000, &
    omi_geo_amf = -1, omi_oobview_amf = -2, omi_wfmod_amf = -9, &
    omi_bigsza_amf = 1000, omi_ooblut_amf = 2000, omi_scattfail_amf = -3

  ! -----------------------
  ! Arrays for OMI L1b data
  ! -----------------------
  REAL    (KIND=r4), DIMENSION (0:nlines_max-1)                        :: omi_auraalt
  REAL    (KIND=r8), DIMENSION (0:nlines_max-1), target                :: omi_time
!unused  INTEGER (KIND=i4), DIMENSION (0:nlines_max-1)                        :: omi_radiance_errstat
  INTEGER (KIND=i1), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_xtrflg_l1b
  INTEGER (KIND=i4), DIMENSION (nxtrack_max,0:nlines_max-1), target    :: omi_geoflg
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_xtrflg
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1), target    :: omi_height, land_water_flg
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1), target    :: omi_latitude, omi_longitude
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1), target    :: omi_szenith, omi_sazimuth
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1), target    :: omi_vzenith, omi_vazimuth
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1), target    :: snow_ice_fraction
  !REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_spec
  REAL    (KIND=r8), DIMENSION (:,:,:), allocatable :: omi_radiance_spec
  !REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_wavl
  REAL    (KIND=r8), DIMENSION (:,:,:), allocatable :: omi_radiance_wavl
  !INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_qflg
  INTEGER (KIND=i2), DIMENSION (:,:,:), allocatable :: omi_radiance_qflg
  !INTEGER (KIND=i4), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_ccdpix
  INTEGER (KIND=i4), DIMENSION (:,:,:), allocatable :: omi_radiance_ccdpix
  INTEGER (KIND=i2), DIMENSION (nUTCdim,0:nlines_max-1)                :: omi_time_utc

  ! ---------------------------------------------------------------
  ! Snow/Ice and Glint flags are used in the AMF computation module
  ! outside the "nlines_max" loops and hence need to be defined on
  ! the maximum swath dimensions.
  ! ---------------------------------------------------------------
  !INTEGER (KIND=i4), DIMENSION (nwavel_max,nxtrack_max) :: omi_irradiance_ccdpix
  INTEGER (KIND=i4), DIMENSION (:,:), allocatable :: omi_irradiance_ccdpix
  !INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max) :: omi_radref_qflg
  INTEGER (KIND=i2), DIMENSION (:,:), allocatable :: omi_radref_qflg
  !REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max) :: &
  !  omi_irradiance_wght, omi_radref_spec, omi_radref_wavl, omi_radref_wght
  REAL    (KIND=r8), DIMENSION (:,:), allocatable :: &
    omi_irradiance_wght, omi_radref_spec, omi_radref_wavl, omi_radref_wght
  !REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max) :: &
  !  omi_irradiance_prec, omi_irradiance_wavl, omi_irradiance_spec, &
  !INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max) :: omi_irradiance_qflg
  REAL    (KIND=r4), DIMENSION (nxtrack_max) :: omi_radref_sza, omi_radref_vza

  ! ---------------------------------------------------------------
  ! Scene Albedo, from OMCLDO2, for Wavelength-Modified AMF lookup
  ! ---------------------------------------------------------------
  REAL (KIND=r4), DIMENSION (nxtrack_max) :: omi_scene_albedo

  INTEGER (KIND=i4), DIMENSION (nxtrack_max,4) :: rad_ccdpix_selection
  INTEGER (KIND=i4), DIMENSION (nxtrack_max,2) :: rad_ccdpix_exclusion

  ! ----------------------------------------
  ! Arrays for fitting and/or derived output
  ! ----------------------------------------
  REAL    (KIND=r8), PARAMETER                              :: d_comm_wvl = 0.01_r8
  INTEGER (KIND=i4)                                         :: n_comm_wvl
  !INTEGER (KIND=i4), DIMENSION (nxtrack_max)                :: common_cnt
  !REAL    (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts)   :: common_spc, common_wvl
  REAL    (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1), target  :: &
    omi_column_amount, omi_column_uncert, &
    omi_fit_rms, omi_radfit_chisq
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_razimuth
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1), target :: omi_fitconv_flag
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1), target :: omi_itnum_flag

  ! ----------------------------------------------------------------------------
  ! Correlations with main output product. Due to a bug in the HDF-EOS5 routines
  ! (non-TLCF implementation), STRING fields cannot be written to file directly.
  ! A work-around solution is to convert the CHARACTERs to INTEGERs. Thus the
  ! need for the additional array CORRELATION_NAMES_INT.
  ! ----------------------------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN), DIMENSION (n_max_fitpars) :: correlation_names
  CHARACTER (LEN=n_max_fitpars*MAX_STR_LEN)              :: correlation_names_concat

  ! --------------------------------------------------------
  ! Ozone is a special case: We can have up to 3 temperatues
  ! --------------------------------------------------------
  REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx, nxtrack_max,0:nlines_max-1) :: &
    omi_o3_amount, omi_o3_uncert

  ! ---------------------------------
  ! Dimensions for measurement swaths
  ! ---------------------------------
  !INTEGER (KIND=i4) :: ntimes, ntimessmallpixel, nxtrack, nwavel, ntimes_loop, nwavel_ccd
  INTEGER (KIND=i4) :: ntimessmallpixel, nwavel, ntimes_loop
  !INTEGER (KIND=i4), DIMENSION (nxtrack_max)                  :: omi_nwav_irrad, omi_nwav_radref
  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                  :: omi_nwav_radref
  INTEGER (KIND=i4), DIMENSION (nxtrack_max,0:nlines_max-1)   :: omi_nwav_rad

  ! ---------------------------------------
  ! Swath attributes for measurement swaths
  ! ---------------------------------------
  INTEGER (KIND=i2) :: ImageBinningFactor, BinnedImageRows, StopColumn
  INTEGER (KIND=i4) :: NumTimes, NumTimesSmallPixel

  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                         :: n_omi_database_wvl
  INTEGER (KIND=i2), DIMENSION (nxtrack_max), TARGET                 :: &
    omi_solcal_xflag, omi_solcal_itnum, omi_radcal_xflag, omi_radref_xflag
  REAL    (KIND=r8), DIMENSION (max_calfit_idx, nxtrack_max)         :: &
    omi_solcal_pars,  omi_radcal_pars,  omi_radref_pars
  !REAL    (KIND=r8), DIMENSION (nwavel_max, nxtrack_max, max_rs_idx) :: omi_database
  REAL    (KIND=r8), DIMENSION (:,:,:), allocatable                  :: omi_database
  !REAL    (KIND=r8), DIMENSION (nwavel_max, nxtrack_max            ) :: omi_database_wvl
  REAL    (KIND=r8), DIMENSION (:,:), allocatable                    :: omi_database_wvl
  !JCH -> unused: REAL    (KIND=r8), DIMENSION (nxtrack_max) :: omi_radref_wav_avg
  REAL    (KIND=r8), DIMENSION (nxtrack_max), TARGET :: &
    omi_solcal_chisq, omi_solcal_shift, omi_radcal_chisq, omi_solcal_rms, &
    omi_radref_chisq, omi_radref_col,   omi_radref_dcol,  &
    omi_radref_rms, omi_radref_xtrcol
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max,0:nlines_max-1)        :: omi_wavwin_rad, omi_fitwin_rad
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max)                       :: omi_wavwin_sol, omi_fitwin_sol

  ! ---------------
  ! OMI swath names
  ! ---------------
  CHARACTER (LEN=MAX_STR_LEN) :: omi_radiance_swathname, omi_irradiance_swathname
  !CHARACTER (LEN=MAX_STR_LEN) :: l1b_radiance_esdt

  ! ------------------------------
  ! Distance between Earth and Sun
  ! ------------------------------
  REAL (KIND=r4) :: EarthSunDistance

  ! ---------------------------------------------------------
  ! OMI scan line, block line, and across-track pixel numbers
  ! ---------------------------------------------------------
  !INTEGER (KIND=i4) :: omi_scanline_no, omi_blockline_no, omi_xtrackpix_no
  INTEGER (KIND=i4) :: omi_blockline_no

  ! ---------------------------
  ! OMI L2 output data QA flags
  ! ----------------------------------------------------------------
  ! Per Centages of output column data that are used to classify the
  ! scientific data quality:
  !         Good data >= QA_PERCENT_PASS   : "Passed"
  !         Good data >= QA_PERCENT_SUSPECT: "Suspect"
  !         Good data <  QA_PERCENT_SUSPECT: "Failed"
  ! ----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: qa_percent_passed = 75, qa_percent_suspect = 50

  ! ------------------------------------------------------
  ! Finally some variables that will be initialized in the
  ! course of the processing.
  ! ------------------------------------------------------
  INTEGER (KIND=i4) :: &
    n_omi_radwvl, n_omi_radrefwvl,  &
    nwavelcoef_irrad, nwavelcoef_rad,               &
    ntimes_smapix_irrad, ntimes_smpix_rad, nclenfit

  LOGICAL, DIMENSION (nxtrack_max) ::  omi_cross_track_skippix = .FALSE.

  ! ------------------------------------
  ! Number of significant digits to keep
  ! ------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_roff_dig = 5

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! ----------------------------------------------------------
  CHARACTER (LEN=16), PARAMETER :: uv1_glob_swath = 'Earth UV-1 Swath'
  CHARACTER (LEN=10), PARAMETER ::    &
    uv1_zoom_swath = '(60x159x4)', &
    uv2_zoom_swath = '(60x557x4)', &
    vis_zoom_swath = '(60x751x4)'

  INTEGER (KIND=i4), PARAMETER :: gzoom_spix = 16, gzoom_epix = 45, gzoom_npix = 30
  INTEGER (KIND=i1), PARAMETER :: global_mode = 8_i1, szoom_mode = 4_i1
  INTEGER (KIND=i1)            :: truezoom, fullswath

  INTEGER (KIND=i4) :: omi_l1b_idx

contains

  subroutine initialize_omidata_structs (errstat)
    use tell_module
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    allocate (omi_radiance_spec(nwavel_max,nxtrack_max,0:nlines_max-1), &
              omi_radiance_wavl(nwavel_max,nxtrack_max,0:nlines_max-1), &
              omi_radiance_qflg(nwavel_max,nxtrack_max,0:nlines_max-1), &
              omi_radiance_ccdpix(nwavel_max,nxtrack_max,0:nlines_max-1), &
              omi_irradiance_ccdpix(nwavel_max,nxtrack_max),&
              omi_radref_qflg(nwavel_max, nxtrack_max), &
              omi_irradiance_wght(nwavel_max, nxtrack_max), &
              omi_radref_spec(nwavel_max, nxtrack_max), &
              omi_radref_wavl(nwavel_max, nxtrack_max), &
              omi_radref_wght(nwavel_max, nxtrack_max), &
              omi_database(nwavel_max, nxtrack_max, max_rs_idx), &
              omi_database_wvl(nwavel_max, nxtrack_max), &
              stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "initialize_omidata_structs:  allocate failed", errstat)
      return
    endif
              
    ! FIXME: (JCH)  Eventually, these struct fields will be arrays and
    ! not pointers and they'll probably be initialized elsewhere.  While
    ! pointers are being used, we'll initialize them here.
    ! Note that all the target objects have the 'target' attribute solely
    ! to support these pointers.   If/when these pointers aren't needed
    ! any longer, consider removing the 'target' attribute where necessary.
    input_vars % time => omi_time
    input_vars % latitude => omi_latitude
    input_vars % longitude => omi_longitude
    input_vars % solar_zenith => omi_szenith
    input_vars % solar_azimuth => omi_sazimuth
    input_vars % viewing_zenith => omi_vzenith
    input_vars % viewing_azimuth => omi_vazimuth
    input_vars % snow_ice_fraction => snow_ice_fraction
    input_vars % terrain_height => omi_height
    input_vars % ground_pixel_quality_flag => omi_geoflg
    
    result_vars % column_amount => omi_column_amount
    result_vars % column_uncert => omi_column_uncert
    result_vars % fit_rms_residual => omi_fit_rms
    result_vars % fit_convergence_flag => omi_fitconv_flag
    result_vars % fit_iteration_count => omi_itnum_flag
    result_vars % solcal_convergence_flag => omi_solcal_xflag
    result_vars % solcal_shift => omi_solcal_shift
    result_vars % radcal_convergence_flag => omi_radcal_xflag
    result_vars % radref_convergence_flag => omi_radref_xflag
    result_vars % radref_column_amount => omi_radref_col
    result_vars % radref_column_uncert => omi_radref_dcol
    result_vars % radref_column_xtrfit => omi_radref_xtrcol
    result_vars % radref_fit_rms => omi_radref_rms
  end subroutine initialize_omidata_structs

  subroutine deallocate_omidata_structs (errstat)
    use tell_module
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    if (allocated(omi_radiance_spec)) deallocate(omi_radiance_spec, &
         stat=errstat)
    if (allocated(omi_radiance_wavl)) deallocate(omi_radiance_wavl, &
         stat=errstat)
    if (allocated(omi_radiance_qflg)) deallocate(omi_radiance_qflg, &
         stat=errstat)
    if (allocated(omi_radiance_ccdpix)) deallocate(omi_radiance_ccdpix, &
         stat=errstat)
    if (allocated(omi_irradiance_ccdpix)) deallocate(omi_irradiance_ccdpix, &
         stat=errstat)
    if (allocated(omi_radref_qflg)) deallocate(omi_radref_qflg, &
         stat=errstat)
    if (allocated(omi_irradiance_wght)) deallocate(omi_irradiance_wght, &
         stat=errstat)
    if (allocated(omi_radref_spec)) deallocate(omi_radref_spec, &
         stat=errstat)
    if (allocated(omi_radref_wavl)) deallocate(omi_radref_wavl, &
         stat=errstat)
    if (allocated(omi_radref_wght)) deallocate(omi_radref_wght, &
         stat=errstat)
    if (allocated(omi_database)) deallocate(omi_database, &
         stat=errstat)
    if (allocated(omi_database_wvl)) deallocate(omi_database_wvl, &
              stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "deallocate_omidata_structs:  deallocate failed", errstat)
      return
    endif

!    input_vars % time => omi_time
!    input_vars % latitude => omi_latitude
!    input_vars % longitude => omi_longitude
!    input_vars % solar_zenith => omi_szenith
!    input_vars % solar_azimuth => omi_sazimuth
!    input_vars % viewing_zenith => omi_vzenith
!    input_vars % viewing_azimuth => omi_vazimuth
!    input_vars % terrain_height => omi_height
!
!    result_vars % column_amount => omi_column_amount
!    result_vars % column_uncert => omi_column_uncert
!    result_vars % fit_rms_residual => omi_fit_rms
!    result_vars % fit_convergence_flag => omi_fitconv_flag
!    result_vars % fit_iteration_count => omi_itnum_flag
!    result_vars % solcal_convergence_flag => omi_solcal_xflag
!    result_vars % radcal_convergence_flag => omi_radcal_xflag
!    result_vars % radref_convergence_flag => omi_radref_xflag
!    result_vars % radref_column_amount => omi_radref_col
!    result_vars % radref_column_uncert => omi_radref_dcol
!    result_vars % radref_column_xtrfit => omi_radref_xtrcol
!    result_vars % radref_fit_rms => omi_radref_rms
  end subroutine deallocate_omidata_structs

  subroutine dealloc_retrieval_type (rt, errstat)
    use tell_module
    implicit none
    type (retrieval_type), intent(inout) :: rt
    integer, intent (inout) :: errstat

    if (errstat /= 0) return
    if (allocated (rt%column_amount)) &
         deallocate (rt%column_amount, stat=errstat)
    if (allocated (rt%column_uncertainty) .and. errstat == 0) &
         deallocate (rt%column_uncertainty, stat=errstat)
    if (allocated (rt%rms) .and. errstat == 0) &
         deallocate (rt%rms, stat=errstat)
    if (allocated (rt%latitude) .and. errstat == 0) &
         deallocate (rt%latitude, stat=errstat)
    if (allocated (rt%longitude) .and. errstat == 0) &
         deallocate (rt%longitude, stat=errstat)
    if (allocated (rt%height) .and. errstat == 0) &
         deallocate (rt%height, stat=errstat)
    if (allocated (rt%sza) .and. errstat == 0) &
         deallocate (rt%sza, stat=errstat)
    if (allocated (rt%vza) .and. errstat == 0) &
         deallocate (rt%vza, stat=errstat)
    if (allocated (rt%saa) .and. errstat == 0) &
         deallocate (rt%saa, stat=errstat)
    if (allocated (rt%vaa) .and. errstat == 0) &
         deallocate (rt%vaa, stat=errstat)
    if (allocated (rt%fit_flag) .and. errstat == 0) &
         deallocate (rt%fit_flag, stat=errstat)
    if (allocated (rt%xtr_flag) .and. errstat == 0) &
         deallocate (rt%xtr_flag, stat=errstat)
    if (allocated (rt%time) .and. errstat == 0) &
         deallocate (rt%time, stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "dealloc_retrieval_type: failed", &
           errstat)
      return
    endif
  end subroutine dealloc_retrieval_type

  subroutine alloc_retrieval_type (rt, nxtrack, ntimes, errstat)
    use OMSAO_parameters_module, only: r4_missval, r8_missval, i2_missval
    use tell_module
    implicit none
    type (retrieval_type), intent(inout) :: rt
    integer (kind=i4), intent(in) :: nxtrack, ntimes
    integer, intent(inout) :: errstat
    integer :: locerrstat

    if (errstat /= 0) return

    allocate ( &
      rt%column_amount(nxtrack, 0:ntimes-1), &
      rt%column_uncertainty(nxtrack, 0:ntimes-1), &
      rt%rms(nxtrack, 0:ntimes-1), &
      rt%latitude(nxtrack, 0:ntimes-1), &
      rt%longitude(nxtrack, 0:ntimes-1), &
      rt%height(nxtrack, 0:ntimes-1), &
      rt%sza(nxtrack, 0:ntimes-1), &
      rt%vza(nxtrack, 0:ntimes-1), &
      rt%saa(nxtrack, 0:ntimes-1), &
      rt%vaa(nxtrack, 0:ntimes-1), &
      rt%fit_flag(nxtrack, 0:ntimes-1), &
      rt%xtr_flag(nxtrack, 0:ntimes-1), &
      rt%time(0:ntimes-1), stat=locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "alloc_retrieval_type:  allocation failed", errstat)
      return
    endif
    rt%nxtrack = nxtrack
    rt%ntimes = ntimes

    rt%column_amount = r8_missval
    rt%column_uncertainty = r8_missval
    rt%rms = r8_missval
    rt%latitude = r4_missval
    rt%longitude = r4_missval
    rt%height = r4_missval
    rt%sza = r4_missval
    rt%vza = r4_missval
    rt%saa = r4_missval
    rt%vaa = r4_missval
    rt%fit_flag = i2_missval
    rt%xtr_flag = i2_missval

  end subroutine alloc_retrieval_type

END MODULE OMSAO_omidata_module
