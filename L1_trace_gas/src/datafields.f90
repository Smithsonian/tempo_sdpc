MODULE datafields

  use OMSAO_precision_module
  use OMSAO_parameters_module, ONLY: MAX_STR_LEN
  use iso_c_binding, ONLY: c_ptr

  INTEGER (KIND=i4), PARAMETER, PRIVATE :: maxrank = 3

  TYPE, PUBLIC :: DataField_HE5
    REAL      (KIND=r8) :: FillValue, MissingValue, ScaleFactor, &
      Offset, SpecTemp
    REAL      (KIND=r8), DIMENSION (2) :: ValidRange
    CHARACTER (LEN=MAX_STR_LEN)     :: Name, Dimensions, Units, Title, UniqueFD
    INTEGER   (KIND=i4)                      :: LenUnits, LenTitle, LenUniqueFD
    INTEGER   (KIND=i4)                      :: Swath_ID
    INTEGER   (KIND=i4)                      :: HE5_DataType
    INTEGER   (KIND=i4)                      :: Rank
    INTEGER   (KIND=i4), DIMENSION (maxrank) :: Dims

    type (c_ptr) :: vdata
    logical :: output = .TRUE.
    type (DataField_HE5), pointer :: next
  END TYPE DataField_HE5

  type, public :: Datafield_List_Type
    integer :: num_items = 0
    type (DataField_HE5), pointer :: head => NULL(), tail => NULL()
  end type Datafield_List_Type

  ! Field Name Constants-- use these instead of literals

  ! Swath geolocation fields
  type (Datafield_List_Type), public, save :: Geo_he5fields
  character (len=*), public, parameter :: &
    auraalt_field = "SpacecraftAltitude", &
    time_field    = "Time", &
    utc_field     = "TimeUTC", &
    thgt_field    = "TerrainHeight", &
    scno_field    = "ScanNumber", &
    lat_field     = "Latitude", &
    lon_field     = "Longitude", &
    saa_field     = "SolarAzimuthAngle", &
    sza_field     = "SolarZenithAngle", &
    vaa_field     = "ViewingAzimuthAngle", &
    vza_field     = "ViewingZenithAngle", &
    xtr_field     = "XtrackQualityFlags", &
    extr_field    = "XtrackQualityFlagsExpanded"

  ! Common Mode data fields
  type (Datafield_List_Type), public, save :: comdata_he5fields
  character (len=*), parameter, public :: &
    amfmol_field   = "AirMassFactor", &
    amfdiag_field  = "AirMassFactorDiagnosticFlag", &
    amfgeo_field   = "AirMassFactorGeometric", &
    avgcol_field   = "AverageColumnAmount", &
    avgdcol_field  = "AverageColumnUncertainty", &
    avgrms_field   = "AverageFittingRMS", &
    col_field      = "ColumnAmount", &
    dstrcol_field  = "ColumnAmountDestriped", &
    dcol_field     = "ColumnUncertainty", &
    commspc_field  = "CommonModeSpectrum", &
    commwvl_field  = "CommonModeWavelengths", &
    fitrms_field   = "FittingRMS", &
    fitcon_field   = "FitConvergenceFlag", &
    mainqa_field   = "MainDataQualityFlag", &
    maxcol_field   = "MaximumColumnAmount", &
    pxclat_field   = "PixelCornerLatitudes", &
    pxclon_field   = "PixelCornerLongitudes", &
    pxarea_field   = "PixelArea"

  type (Datafield_List_Type), public, save :: sol_calfit_he5fields
  character (len=*), parameter, public :: &
    swccf_field    = "SolarWavCalConvergenceFlag"

  type (Datafield_List_Type), public, save :: rad_calfit_he5fields
  character (len=*), parameter, public :: &
    rwccf_field    = "RadianceWavCalConvergenceFlag", &
    rwclr_field    = "RadianceWavCalLatitudeRange"

  type (Datafield_List_Type), public, save :: rad_reffit_he5fields
  character (len=*), parameter, public :: &
    rrcf_field     = "RadianceReferenceConvergenceFlag", &
    rrlr_field     = "RadianceReferenceLatitudeRange", &
    rrcol_field    = "RadianceReferenceColumnAmount", &
    rrdcol_field   = "RadianceReferenceColumnUncertainty", &
    rrxcol_field   = "RadianceReferenceColumnXTRFit", &
    rrrms_field    = "RadianceReferenceFittingRMS"

  ! Swath data fields, additional for "diagnostic" runs
  type (Datafield_List_Type), public, save :: diagnostic_he5fields
  character (len=*), parameter, public :: &
    ccdpix_field   = "CCDPixelRange", &
    commcnt_field  = "CommonModeCount", &
    spcfit_field   = "FittedSpectrum", &
    corrcol_field  = "FittingParameterColumns", &
    corr_field     = "FittingParameterCorrelations", &
    correlm_field  = "FittingParameterNames", &
    correrr_field  = "FittingParameterUncertainty", &
    spcobs_field   = "MeasuredSpectrum", &
    posobs_field   = "MeasuredWavelength", &
    spcres_field   = "ResidualSpectrum", &
    fitwt_field    = "SpectralFitWeight", &
    spdata_field   = "DatabaseSpec", &
    spdatw_field   = "DatabaseWavl", &
    spnrmf_field   = "DatabaseNormFactor", &
    spname_field   = "DatabaseSpecNames", &
    itnum_field    = "IterationCount", &
    xtrcor_field   = "CrossTrackStripeCorrection"

  ! Swath data fields unique to OMHCHO and OMCHOCHO
  type (Datafield_List_Type), public, save :: voc_he5fields
  character (len=*), parameter, public :: &
    amfcfr_field   = "AMFCloudFraction", &
    amfctp_field   = "AMFCloudPressure"

  ! Special data fields for wavelength-modified AMF fitting
  type (Datafield_List_Type), public, save :: wmamf_he5fields
  character (len=*), parameter, public :: &
    adalb_field     = "AdjustedSceneAlbedo", &
    scol_field      = "SlantColumnAmount", &
    sdstrcol_field  = "SlantColumnAmountDestriped", &
    sdcol_field     = "SlantColumnUncertainty", &
    sfitrms_field   = "SlantFittingRMS", &
    sfitcon_field   = "SlantFitConvergenceFlag"

  ! Swath data fields for Reference Sector -- GGA
  type (Datafield_List_Type), public, save :: rs_he5fields
  character (len=*), parameter, public :: &
    rscol_field     = "ReferenceSectorCorrectedVerticalColumn", &
    rscod_field     = "ReferenceSectorCorrectedUncertainty"

  ! Swath data field for Scattering Weights, Gas Profile Averaging Kernels and albedo -- GGA
  type (Datafield_List_Type), public, save :: sw_he5fields
  character (len=*), parameter, public :: &
    scaweights_field = "ScatteringWeights", &
    gasprofile_field = "GasProfile", &
    albedo_field     = "Albedo"

  type (Datafield_List_Type), public, save :: o3_prefit_he5fields
  character (len=*), parameter, public :: &
    slantcol_temp1_field = "SlantColumnAmountTemperatureT1", &
    slantcol_temp2_field = "SlantColumnAmountTemperatureT2", &
    slantcol_temp3_field = "SlantColumnAmountTemperatureT3"

  type (Datafield_List_Type), public, save :: o3_prefit_uncert_he5fields
  character (len=*), parameter, public :: &
    slantcol_dtemp1_field = "SlantColumnUncertaintyTemperatureT1", &
    slantcol_dtemp2_field = "SlantColumnUncertaintyTemperatureT2", &
    slantcol_dtemp3_field = "SlantColumnUncertaintyTemperatureT3"

  private

  public he5_initialize_datafields

CONTAINS

  subroutine he5_initialize_datafields ()

    call new_datafield(geo_he5fields, &
                       "Latitude", &                       ! name
                       "Geodetic Latitude", &              ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       -90.0_r8, 90.0_r8, &                ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! r4: omi_latitude (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "Longitude", &                      ! name
                       "Geodetic Longitude", &             ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       -180.0_r8, 180.0_r8, &              ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! r4: omi_longitude (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "SolarAzimuthAngle", &              ! name
                       "Solar Azimuth Angle", &            ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       -180.0_r8, 180.0_r8, &              ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! r4 omi_sazimuth (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "SolarZenithAngle", &               ! name
                       "Solar Zenith Angle", &             ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       0.0_r8, 180.0_r8, &                 ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! r4 omi_szenith (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "SpacecraftAltitude", &             ! name
                       "Altitude of Aura Spacecraft", &    ! title
                       "nTimes", &                         ! dimensions
                       "m", &                              ! units
                       "r4", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! r4 omi_auraalt(0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "TerrainHeight", &                  ! name
                       "Terrain Height", &                 ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "m", &                              ! units
                       "i2", &                             ! datatype
                       -1000.0_r8, 10000.0_r8, &           ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! i2 omi_height (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "Time", &                           ! name
                       "Time in TAI units", &              ! title
                       "nTimes", &                         ! dimensions
                       "s", &                              ! units
                       "r8", &                             ! datatype
                       0.0_r8, 10000000000.0_r8, &         ! validrange
                       "Aura-Shared" &                     ! uniquefd
                      )
    ! r8 omi_auraalt(0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "TimeUTC", &                        ! name
                       "Coordianted Universal Time", &     ! title
                       "nUTCdim,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       0.0_r8, 9999.0_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! i2 omi_time_utc (nUTCdim, 0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "ViewingAzimuthAngle", &            ! name
                       "Viewing Azimuth Angle", &          ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       -180.0_r8, 180.0_r8, &              ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r4 omi_vazimuth (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "ViewingZenithAngle", &             ! name
                       "Viewing Zenith Angle", &           ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       0.0_r8, 180.0_r8, &                 ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r4 omi_vzenith (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "XtrackQualityFlags", &             ! name
                       "Cross-Track Quality Flags", &      ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i1", &                             ! datatype
                       0.0_r8, 127.0_r8, &                 ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! i1 omi_xtrflg_l1b (nxtrack_max,0:nlines_max-1)

    call new_datafield(geo_he5fields, &
                       "XtrackQualityFlagsExpanded", &     ! name
                       "Expanded Cross-Track Quality Flags", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       0.0_r8, 11147.0_r8, &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! i2 omi_xtrflg (nxtrack_max,0:nlines_max-1)

    call new_datafield(sol_calfit_he5fields, &
                       "SolarWavCalConvergenceFlag", &     ! name
                       "Solar Wavelength Calibration Convergence Flag", &! title
                       "nXtrack", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -10.0_r8, 12344.0_r8, &             ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    ! The code in he5_output_tools assumes that these two are sequential and
    ! in this order.
    call new_datafield(rad_calfit_he5fields, &
                       "RadianceWavCalConvergenceFlag", &  ! name
                       "Radiance Wavelength Calibration Convergence Flag", &! title
                       "nXtrack", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -10.0_r8, 12344.0_r8, &             ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(rad_calfit_he5fields, &
                       "RadianceWavCalLatitudeRange", &    ! name
                       "Radiance Wavelength Calibration Latitude Range", &! title
                       "4", &                              ! dimensions
                       "NoUnits", &                        ! units
                       "r4", &                             ! datatype
                       -90.0_r8, 90.0_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(rad_reffit_he5fields, &
                       "RadianceReferenceConvergenceFlag", &! name
                       "Radiance Reference Fit Convergence Flag", &! title
                       "nXtrack", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -10.0_r8, 12344.0_r8, &             ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(rad_reffit_he5fields, &
                       "RadianceReferenceLatitudeRange", & ! name
                       "Radiance Reference Fit Latitude Range", &! title
                       "4", &                              ! dimensions
                       "NoUnits", &                        ! units
                       "r4", &                             ! datatype
                       -90.0_r8, 90.0_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(rad_reffit_he5fields, &
                       "RadianceReferenceColumnAmount", &  ! name
                       "Radiance Reference Fit Column Amount", &! title
                       "nXtrack", &                        ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    !r8 omi_radref_col(nxtrack_max)

    call new_datafield(rad_reffit_he5fields, &
                       "RadianceReferenceColumnUncertainty", &! name
                       "Radiance Reference Fit Column Uncertainty", &! title
                       "nXtrack", &                        ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_radref_dcol(nxtrack_max)

    call new_datafield(rad_reffit_he5fields, &
                       "RadianceReferenceColumnXTRFit", &  ! name
                       "Radiance Reference Fit Column XTR Fit", &! title
                       "nXtrack", &                        ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_radref_xtrcol (nxtrack_max)

    call new_datafield(rad_reffit_he5fields, &
                       "RadianceReferenceFittingRMS", &    ! name
                       "Radiance Reference Fit RMS", &     ! title
                       "nXtrack", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_radref_rms (nxtrack_max)

    call new_datafield(comdata_he5fields, &
                       "AirMassFactor", &                  ! name
                       "Molecule Specific Air Mass Factor (AMF)", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r4", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! derived from amfgeo in amf_calculation

    call new_datafield(comdata_he5fields, &
                       "AirMassFactorDiagnosticFlag", &    ! name
                       "Diagnostic Flag for Molecule Specific AMF", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -2.0_r8, 13127.0_r8, &              ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local array amfdiag from amf_calculation

    call new_datafield(comdata_he5fields, &
                       "AirMassFactorGeometric", &         ! name
                       "Geometric Air Mass Factor (AMF)", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r4", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! derived from amfgeo in amf_calculation

    call new_datafield(comdata_he5fields, &
                       "AverageColumnAmount", &            ! name
                       "Average Column Amount", &          ! title
                       "1", &                              ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local from compute_fitting_statistics

    call new_datafield(comdata_he5fields, &
                       "AverageColumnUncertainty", &       ! name
                       "Average Column Uncertainty", &     ! title
                       "1", &                              ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local from compute_fitting_statistics

    call new_datafield(comdata_he5fields, &
                       "AverageFittingRMS", &              ! name
                       "Average Fitting RMS", &            ! title
                       "1", &                              ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local from compute_fitting_statistics

    call new_datafield(comdata_he5fields, &
                       "ColumnAmount", &                   ! name
                       "Column Amount", &                  ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_column_amount (nxtrack_max,0:nlines_max-1)

    call new_datafield(comdata_he5fields, &
                       "ColumnAmountDestriped", &          ! name
                       "Column Amount with Destriping Correction", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(comdata_he5fields, &
                       "ColumnUncertainty", &              ! name
                       "Column Uncertainty", &             ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_column_uncert (nxtrack_max,0:nlines_max-1)

    call new_datafield(comdata_he5fields, &
                       "FittingRMS", &                     ! name
                       "Fitting RMS", &                    ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_fit_rms (nxtrack_max,0:nlines_max-1)

    call new_datafield(comdata_he5fields, &
                       "FitConvergenceFlag", &             ! name
                       "Fitting Convergence Flag", &       ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -10.0_r8, 12344.0_r8, &             ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! i2 omi_fitconv_flag (nxtrack_max,0:nlines_max-1)

    call new_datafield(comdata_he5fields, &
                       "MainDataQualityFlag", &            ! name
                       "Main Data Quality Flag", &         ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -1.0_r8, 2.0_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local variable array (saomqf)

    call new_datafield(comdata_he5fields, &
                       "MaximumColumnAmount", &            ! name
                       "Maximum Column Amount for QA Flag 'good'", &! title
                       "1", &                              ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local variable max_col from compute_fitting_statistics
    call new_datafield(comdata_he5fields, &
                       "PixelCornerLatitudes", &           ! name
                       "Pixel Corner Latitude Coordinates", &! title
                       "nXtrack+1,nTimes+1", &             ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       -90.0_r8, 90.0_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(comdata_he5fields, &
                       "PixelCornerLongitudes", &          ! name
                       "Pixel Corner Longitude Coordinates", &! title
                       "nXtrack+1,nTimes+1", &             ! dimensions
                       "deg", &                            ! units
                       "r4", &                             ! datatype
                       -180.0_r8, 180.0_r8, &              ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(comdata_he5fields, &
                       "PixelArea", &                      ! name
                       "Pixel Area", &                     ! title
                       "nXtrack", &                        ! dimensions
                       "km^2", &                           ! units
                       "r4", &                             ! datatype
                       0.0_r8, 180.0_r8*360.0_r8*6378.0_r8, &! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(diagnostic_he5fields, &
                       "CCDPixelRange", &                  ! name
                       "First and Last CCD Pixel Number Fitted", &! title
                       "nXtrack,2", &                      ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       1.0_r8, 780.0_r8, &                 ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! derived from common_mode_spec%CCDPixel

    call new_datafield(diagnostic_he5fields, &
                       "CommonModeCount", &                ! name
                       "Common Mode Spectrum Averaging Count", &! title
                       "nXtrack", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "i4", &                             ! datatype
                       0.0_r8, 2147483647.0_r8, &          ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! derived from common_mode_spec%RefSpecCount

    call new_datafield(diagnostic_he5fields, &
                       "CommonModeSpectrum", &             ! name
                       "Common Mode Spectrum", &           ! title
                       "nXtrack,nCommonWavl", &            ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! derived from common_mode_spec%RefSpecData
    call new_datafield(diagnostic_he5fields, &
                       "CommonModeWavelengths", &          ! name
                       "Common Mode Spectrum Wavelengths", &! title
                       "nXtrack,nCommonWavl", &            ! dimensions
                       "nm", &                             ! units
                       "r4", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! derived from common_mode_spec%RefSpecWavs

    call new_datafield(diagnostic_he5fields, &
                       "CrossTrackStripeCorrection", &     ! name
                       "Correction Factor for Cross Track Stripes", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnit", &                         ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(diagnostic_he5fields, &
                       "FittingParameterColumns", &        ! name
                       "Colum Values of all Fitting Parameters", &! title
                       "nFitElements,nXtrack,nTimes", &    ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! all_fitted_columns allocated

    call new_datafield(diagnostic_he5fields, &
                       "FittingParameterCorrelations", &   ! name
                       "Correlations with Main Fitting Parameter", &! title
                       "nFitElements,nXtrack,nTimes", &    ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 correlation columns allocated in code

    call new_datafield(diagnostic_he5fields, &
                       "FittedSpectrum", &                 ! name
                       "Fitted Spectrum", &                ! title
                       "nCommonWavl,nXtrack,nTimes", &     ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! omi_fitspc (:,:,1,:) allocated in code

    call new_datafield(diagnostic_he5fields, &
                       "MeasuredSpectrum", &               ! name
                       "Observed L1B Radiance", &          ! title
                       "nCommonWavl,nXtrack,nTimes", &     ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! omi_fitspc (:, :, 2, :) allocated in code

    call new_datafield(diagnostic_he5fields, &
                       "MeasuredWavelength", &             ! name
                       "Spectral Position of L1B Radiance", &! title
                       "nCommonWavl,nXtrack,nTimes", &     ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! omi_fitspc (:, :, 3, :) allocated in code

    call new_datafield(diagnostic_he5fields, &
                       "SpectralFitWeight", &              ! name
                       "Weighting of L1B Radiance Pixel", &! title
                       "nCommonWavl,nXtrack,nTimes", &     ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! omi_fitspc (:, :, 4, :) allocated in code

    call new_datafield(diagnostic_he5fields, &
                       "FittingParameterNames", &          ! name
                       "Names of all Fitting Parameters", &! title
                       "nCharLenFitElements", &            ! dimensions
                       "NoUnits", &                        ! units
                       "ch", &                             ! datatype
                       0.0_r8, 32767.0_r8, &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! correlation_names_concat (see comment in OMSAO_omidata_module.f90)

    call new_datafield(diagnostic_he5fields, &
                       "FittingParameterUncertainty", &    ! name
                       "Uncertainties of all Fitting Parameters", &! title
                       "nFitElements,nXtrack,nTimes", &    ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! all_fitted_errors allocated

    call new_datafield(diagnostic_he5fields, &
                       "DatabaseSpec", &                   ! name
                       "Reference Spectra used in fitting process", &! title
                       "nRfSpec,nwavel_max,nXtrack", &     ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local array derived from omi_database

    call new_datafield(diagnostic_he5fields, &
                       "DatabaseWavl", &                   ! name
                       "Reference Spectra Wavelengths", &  ! title
                       "nwavel_max,nXtrack", &             ! dimensions
                       "nm", &                             ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local array derived from omi_database_wvl

    call new_datafield(diagnostic_he5fields, &
                       "DatabaseNormFactor", &             ! name
                       "Normalisation factors for DatabaseSpec", &! title
                       "nRfSpec", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! refspecs_original(ii)%NormFactor

    call new_datafield(diagnostic_he5fields, &
                       "DatabaseSpecNames", &              ! name
                       "Species Names for DatabaseSpec", & ! title
                       "nRfSpec", &                        ! dimensions
                       "NoUnits", &                        ! units
                       "ch", &                             ! datatype
                       0.0_r8, 32767.0_r8, &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(diagnostic_he5fields, &
                       "IterationCount", &                 ! name
                       "Radiance Fit Iteration Count", &   ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       0.0_r8, 32767.0_r8, &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! i2 omi_itnum_flag (nxtrack_max,0:nlines_max-1)

    call new_datafield(diagnostic_he5fields, &
                       "ResidualSpectrum", &               ! name
                       "Residual Spectrum", &              ! title
                       "nCommonWavl,nXtrack,nTimes", &     ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! omi_fitspc(:, :, 2, :) - omi_fitspc (:,:,1,:) allocated in code

    call new_datafield(voc_he5fields, &
                       "AMFCloudFraction", &               ! name
                       "Adjusted Cloud Fraction for AMF Computation", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r4", &                             ! datatype
                       0.0_r8, 1.0_r8, &                   ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! local variable amfcfr

    call new_datafield(voc_he5fields, &
                       "AMFCloudPressure", &               ! name
                       "Adjusted Cloud Pressure for AMF Computation", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "hPa", &                            ! units
                       "r4", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    !local variable amfctp

    call new_datafield(wmamf_he5fields, &
                       "AdjustedSceneAlbedo", &            ! name
                       "Adjusted Scene Albedo", &          ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r4", &                             ! datatype
                       0.0_r8, 1.0_r8, &                   ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(wmamf_he5fields, &
                       "SlantColumnAmount", &              ! name
                       "Slant Column Amount", &            ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(wmamf_he5fields, &
                       "SlantColumnAmountDestriped", &     ! name
                       "Slant Column Amount with Destriping Correction", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(wmamf_he5fields, &
                       "SlantColumnUncertainty", &         ! name
                       "Slant Column Uncertainty", &       ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(wmamf_he5fields, &
                       "SlantFittingRMS", &                ! name
                       "Slant Fitting RMS", &              ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       0.0_r8, 1e30_r8, &                  ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(wmamf_he5fields, &
                       "SlantFitConvergenceFlag", &        ! name
                       "Slant Fitting Convergence Flag", & ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "i2", &                             ! datatype
                       -10.0_r8, 12344.0_r8, &             ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(rs_he5fields, &
                       "ReferenceSectorCorrectedVerticalColumn", &! name
                       "Reference Sector Corrected Vertical Column", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(rs_he5fields, &
                       "ReferenceSectorCorrectedUncertainty", &! name
                       "Reference Sector Corrected Uncertainty", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(sw_he5fields, &
                       "ScatteringWeights", &              ! name
                       "Scattering Weights", &             ! title
                       "nXtrack,nTimes,nLevels", &         ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(sw_he5fields, &
                       "GasProfile", &                     ! name
                       "Gas Profile", &                    ! title
                       "nXtrack,nTimes,nLevels", &         ! dimensions
                       "ppb", &                            ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(sw_he5fields, &
                       "Albedo", &                         ! name
                       "Albedo", &                         ! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "NoUnits", &                        ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8, &                ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )

    call new_datafield(o3_prefit_he5fields, &
                       "SlantColumnAmountTemperatureT1", & ! name
                       "Slant Column Amount at Temperature T1", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8 , &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_o3_amount_t1 (nxtrack_max,0:nlines_max-1)

    call new_datafield(o3_prefit_he5fields, &
                       "SlantColumnAmountTemperatureT2", & ! name
                       "Slant Column Amount at Temperature T2", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8 , &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_o3_amount_t2 (nxtrack_max,0:nlines_max-1)

    call new_datafield(o3_prefit_he5fields, &
                       "SlantColumnAmountTemperatureT3", & ! name
                       "Slant Column Amount at Temperature T3", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r8", &                             ! datatype
                       -1e30_r8, 1e30_r8 , &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_o3_amount_t3 (nxtrack_max,0:nlines_max-1)

    call new_datafield(o3_prefit_uncert_he5fields, &
                       "SlantColumnUncertaintyTemperatureT1", &! name
                       "Slant Column Uncertainty at Temperature T1", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r4", &                             ! datatype
                       -1e30_r8, 1e30_r8 , &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_o3_amount_uncert_t1 (nxtrack_max,0:nlines_max-1)

    call new_datafield(o3_prefit_uncert_he5fields, &
                       "SlantColumnUncertaintyTemperatureT2", &! name
                       "Slant Column Uncertainty at Temperature T2", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r4", &                             ! datatype
                       -1e30_r8, 1e30_r8 , &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_o3_amount_uncert_t2 (nxtrack_max,0:nlines_max-1)

    call new_datafield(o3_prefit_uncert_he5fields, &
                       "SlantColumnUncertaintyTemperatureT3", &! name
                       "Slant Column Uncertainty at Temperature T3", &! title
                       "nXtrack,nTimes", &                 ! dimensions
                       "molec/cm2", &                      ! units
                       "r4", &                             ! datatype
                       -1e30_r8, 1e30_r8 , &               ! validrange
                       "OMI-Specific" &                    ! uniquefd
                      )
    ! r8 omi_o3_amount_uncert_t3 (nxtrack_max,0:nlines_max-1)

  end subroutine

  ! ---------------------------------------------------------------------------

  subroutine new_datafield (list, name, title, dimstring, units, &
                            dtype, minv, maxv, uniquefd)

    use OMSAO_he5_module, only: pge_swath_id
    implicit none

    type (Datafield_List_Type), intent(inout) :: list
    character (len=*), intent(in) :: name, title, dimstring, uniquefd, &
      units, dtype
    real (kind=r8), intent(in) :: minv, maxv
                                                           ! local vars
    type (DataField_HE5), pointer :: item
    INTEGER (KIND=i4) :: numtype
    REAL    (KIND=r8) :: missval

    allocate (item)
    item%next => NULL()

    if (.not.associated(list%head)) then
      list % head => item
    else
      list % tail % next => item
    endif
    list % tail => item

    item % name = adjustl(name)
    item % dimensions = adjustl(dimstring)
    item % units = adjustl(units)
    item % lenunits = len_trim(item%units)
    item % title = adjustl(title)
    item % lentitle = len_trim(item%title)
    item % uniquefd = uniquefd
    item % lenuniquefd = len_trim(item%uniquefd)
    item % validrange = (/ minv, maxv /)

    item % dims = -1
    item % rank = 1                                        !FIXME: is this correct??

    CALL find_datatype (dtype, numtype, missval)
    item % fillvalue = missval
    item % missingvalue = missval

    item % he5_datatype = numtype
    item % scalefactor = 1.0_r8
    item % offset = 0.0_r8

    item % swath_id = pge_swath_id

  end subroutine new_datafield

  SUBROUTINE find_datatype ( tyc, numtype, missval )

    USE OMSAO_parameters_module, ONLY: &
      r8_missval, r4_missval, i4_missval, i2_missval, i1_missval
    USE OMSAO_he5_module, ONLY: HE5T_NATIVE_CHAR, HE5T_NATIVE_INT, &
      HE5T_NATIVE_FLOAT, HE5T_NATIVE_DOUBLE, HE5T_NATIVE_INT8, &
      HE5T_NATIVE_INT16
    IMPLICIT NONE

    CHARACTER (LEN=2), INTENT (IN) :: tyc

    INTEGER (KIND=i4), INTENT (OUT) :: numtype
    REAL    (KIND=r8), INTENT (OUT) :: missval

    numtype = 0  ;  missval = 0.0_r8

    SELECT CASE (tyc)
    CASE ("ch")
      numtype = HE5T_NATIVE_CHAR
      missval = REAL ( i1_missval, KIND=r8 )
    CASE ("i1")
      numtype = HE5T_NATIVE_INT8
      missval = REAL ( i1_missval, KIND=r8 )
    CASE ("i2")
      numtype = HE5T_NATIVE_INT16
      missval = REAL ( i2_missval, KIND=r8 )
    CASE ("i4")
      numtype = HE5T_NATIVE_INT
      missval = REAL ( i4_missval, KIND=r8 )
    CASE ("r4")
      numtype = HE5T_NATIVE_FLOAT
      missval = REAL ( r4_missval, KIND=r8 )
    CASE ("r8")
      numtype = HE5T_NATIVE_DOUBLE
      missval = r8_missval
    END SELECT

    RETURN
  END SUBROUTINE find_datatype

end module
