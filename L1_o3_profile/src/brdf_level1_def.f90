MODULE level1_def
  
  USE error_module,      ONLY : ErrorType,CheckError
  USE time_module,       ONLY : TimeType
  USE parameters_module
  USE covariance_module, ONLY : UncertType
  
  TYPE AuxL2FileOptType
    LOGICAL                :: OverwriteProfileClim
    LOGICAL                :: InLevel1File
    CHARACTER(LEN=maxChar) :: FileType
    CHARACTER(LEN=maxChar) :: FileName
    INTEGER                :: ncid
    INTEGER                :: geoid
    INTEGER                :: auxid
  ENDTYPE AuxL2FileOptType
  
  TYPE AuxL2ProfOptType
    LOGICAL                             :: OverwriteProfileClim
    CHARACTER(LEN=maxChar)              :: FileType
    CHARACTER(LEN=maxChar)              :: FileName
    INTEGER                             :: ncid
    INTEGER                             :: prfid
    INTEGER                             :: nGas
    CHARACTER(LEN=maxChar), ALLOCATABLE :: GasName(:)
    INTEGER                             :: nAer
    CHARACTER(LEN=maxChar), ALLOCATABLE :: AerName(:)
    INTEGER                             :: nCld ! # optical prop types
    CHARACTER(LEN=maxChar), ALLOCATABLE :: CldName(:)
  ENDTYPE AuxL2ProfOptType

  TYPE AuxL2CloudOptType
    LOGICAL                             :: OverwriteProfileClim
    LOGICAL                             :: UseL2MetCloud
    INTEGER                             :: nType
    CHARACTER(LEN=maxChar), ALLOCATABLE :: TypeName(:)

    LOGICAL                             :: UseAuxLambertian
    TYPE(AuxL2FileOptType)              :: AuxLambertian ! Alternative 2D Cloud Fr/Press
  ENDTYPE AuxL2CloudOptType 

  TYPE AuxL2SurfOptType
    LOGICAL                             :: OverwriteSurface
    CHARACTER(LEN=maxChar), ALLOCATABLE :: FileType(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: FileName(:)
    INTEGER,                ALLOCATABLE :: FileBandIndex(:)
  ENDTYPE AuxL2SurfOptType

  TYPE GeolocationType
    REAL(KIND=8)   :: SZA
    REAL(KIND=8)   :: VZA
    REAL(KIND=8)   :: AZA
    REAL(KIND=8)   :: VAA
    REAL(KIND=8)   :: SAA
    REAL(KIND=8)   :: Longitude
    REAL(KIND=8)   :: Latitude
    REAL(KIND=8)   :: CornerLongitudes(4)
    REAL(KIND=8)   :: CornerLatitudes(4)
    REAL(KIND=8)   :: SurfaceAltitude
    REAL(KIND=8)   :: ObservationAltitude
    TYPE(TimeType) :: Time
    REAL(KIND=8)   :: OpticalBenchTemperature
  ENDTYPE GeolocationType
  
  ! Auxiliary information packaged with the L1 (e.g. other L2 products)
  TYPE AuxiliaryL1Type
    REAL(KIND=8) :: CloudFraction
    REAL(KIND=8) :: CloudPressure
    REAL(KIND=8) :: WindSpeed
    REAL(KIND=8) :: WindDirection
    REAL(KIND=8) :: Chlorophyll
    REAL(KIND=8) :: OceanSalinity
    REAL(KIND=8) :: SnowFraction
    REAL(KIND=8) :: SeaIceFraction
    REAL(KIND=8) :: SnowDepth
    REAL(KIND=8) :: SnowAge
    REAL(KIND=8) :: LandCoverFraction(17)
  ENDTYPE AuxiliaryL1Type
  
  TYPE AuxProfType
    
    ! Vertical Dimension
    INTEGER                     :: lmx
    INTEGER(KIND=2)             :: GridTypeIndex
    
    ! Met
    REAL(KIND=8), ALLOCATABLE   :: AP(:)
    REAL(KIND=8), ALLOCATABLE   :: BP(:)
    REAL(KIND=8), ALLOCATABLE   :: PressureEdge(:)
    REAL(KIND=8), ALLOCATABLE   :: TemperatureEdge(:)
    REAL(KIND=8), ALLOCATABLE   :: RH(:)
    TYPE(UncertType)            :: TempUncertainty
    TYPE(UncertType)            :: SurfPresUncertainty
    
    ! Gas 
    REAL(KIND=8),   ALLOCATABLE   :: GasMixingRatio(:,:)
    REAL(KIND=8),   ALLOCATABLE   :: GasBackgroundMixingRatio(:,:)
    INTEGER(KIND=2),ALLOCATABLE   :: GasSigmaAdjustType(:)
    TYPE(UncertType), ALLOCATABLE :: GasUncertainty(:)
    
    ! Aerosol
    REAL(KIND=8),           ALLOCATABLE :: AerLayerOpticalDepth(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AerColumnOpticalDepth(:)
    REAL(KIND=8),           ALLOCATABLE :: AerAltMin(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerAltMax(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerAltPeak(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerAltSigma(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerAltExp(:) ![km]
    TYPE(UncertType),       ALLOCATABLE :: AerUncertainty(:)
    INTEGER(KIND=2),        ALLOCATABLE :: AerTypeIndex(:)
    INTEGER,                ALLOCATABLE :: AernProfilePar(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: AerProfileParName(:,:)

    ! Cloud profile 
    LOGICAL                   :: UseL2MetCloud
    LOGICAL                   :: LambertianCloudsPresent
    LOGICAL                   :: ScatteringCloudsPresent
    INTEGER                   :: MaxCldPix
    INTEGER(KIND=2)           :: nCldPix
    REAL(KIND=8), ALLOCATABLE :: CldPixFraction(:)
    REAL(KIND=8), ALLOCATABLE :: CldPixPressure(:)
    LOGICAL                   :: L2LambertianCloudsExist
    LOGICAL                   :: L2ScatteringCloudsExist

    ! Cloud Optics fields
    REAL(KIND=8),           ALLOCATABLE :: CldLayerOpticalDepth(:,:,:)
    REAL(KIND=8),           ALLOCATABLE :: CldColumnOpticalDepth(:,:)
    INTEGER(KIND=2),        ALLOCATABLE :: CldTypeIndex(:)

    ! Currently these are not implemented
    REAL(KIND=8),           ALLOCATABLE :: CldAltMin(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldAltMax(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldAltPeak(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldAltSigma(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldAltExp(:,:) ![km]
    ! TYPE(UncertType),       ALLOCATABLE :: CldUncertainty(:)
    
    INTEGER,                ALLOCATABLE :: CldnProfilePar(:,:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: CldProfileParName(:,:,:)
    
  ENDTYPE AuxProfType
  
  TYPE AuxSurfType
    LOGICAL                      :: OverwriteSurface
    INTEGER                      :: ncid  ! File
    INTEGER                      :: gid   ! Surface Group ID
    INTEGER                      :: kmx ! # Kernels
    INTEGER                      :: wmx ! Wavelength or Chebyshev Coeff dim
    INTEGER                      :: pmx ! Parameter
    INTEGER(KIND=2), ALLOCATABLE :: kern_idx(:)
    REAL(KIND=8),    ALLOCATABLE :: kern_amp(:,:)
    REAL(KIND=8),    ALLOCATABLE :: kern_par(:,:,:)
    LOGICAL                      :: isLambertian
    LOGICAL                      :: isChebyshevParam
    REAL(KIND=8),    ALLOCATABLE :: wvl(:) ! Wavelength
    REAL(KIND=8)                 :: WvlMin ! Chebyshev bounds
    REAL(KIND=8)                 :: WvlMax
  ENDTYPE AuxSurfType

  TYPE SinglePixelType
    TYPE(GeolocationType)               :: Geolocation
    TYPE(AuxiliaryL1Type)               :: Auxiliary
    INTEGER                             :: nBand
    REAL(KIND=8),           ALLOCATABLE :: StartWvl(:)
    REAL(KIND=8),           ALLOCATABLE :: EndWvl(:)
    REAL(KIND=8),           ALLOCATABLE :: dWvl(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: BandName(:)
  ENDTYPE SinglePixelType
  
  TYPE L1OptType
    LOGICAL                             :: DoSinglePixelCalc
    TYPE(SinglePixelType)               :: SinglePix
    INTEGER                             :: nBand
    INTEGER,                ALLOCATABLE :: BandIndex(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: Infile(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: FileType(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: BandName(:)
    TYPE(AuxL2CloudOptType)             :: AuxCloud
    TYPE(AuxL2FileOptType)              :: AuxSurfAlt
    TYPE(AuxL2FileOptType)              :: AuxSurfWind
    TYPE(AuxL2FileOptType)              :: AuxChlorophyll
    TYPE(AuxL2FileOptType)              :: AuxSalinity
    TYPE(AuxL2FileOptType)              :: AuxSnow
    TYPE(AuxL2ProfOptType)              :: AuxProfile
    TYPE(AuxL2SurfOptType)              :: AuxSurface
  ENDTYPE L1OptType
  
  TYPE L1BandType
    REAL(KIND=8),    ALLOCATABLE :: Radiance(:)
    REAL(KIND=8),    ALLOCATABLE :: RadianceUncertainty(:)
    REAL(KIND=8)                 :: RadianceFillValue
    REAL(KIND=8),    ALLOCATABLE :: Wavelength(:)
    INTEGER(KIND=2), ALLOCATABLE :: RadianceFlags(:)
    INTEGER                      :: nwvl
    CHARACTER(LEN=maxChar)       :: ItfName
    CHARACTER(LEN=maxChar)       :: RadianceUnit
    REAL(KIND=8)                 :: SCD ! Used for AMF Mode
    REAL(KIND=8)                 :: SCDUncertainty
    REAL(KIND=8)                 :: AMFWavelength

    LOGICAL                      :: ReplaceGeolocationSZA
    LOGICAL                      :: ReplaceGeolocationVZA
    LOGICAL                      :: ReplaceGeolocationAZA
    REAL(KIND=8)                 :: SZA
    REAL(KIND=8)                 :: VZA
    REAL(KIND=8)                 :: AZA
    
    ! Future Update - Direct Multigeometry input (more efficient for VLIDORT)
    ! INTEGER                      :: nGeometries
    ! REAL(KIND=8),    ALLOCATABLE :: SZA(:)
    ! REAL(KIND=8),    ALLOCATABLE :: VZA(:)
    ! REAL(KIND=8),    ALLOCATABLE :: AZA(:)
    ! REAL(KIND=8),    ALLOCATABLE :: VAA(:)
    ! REAL(KIND=8),    ALLOCATABLE :: SAA(:)
  ENDTYPE L1BandType
  
  TYPE L1Type
    CHARACTER(LEN=maxChar)              :: CalculationMode ! Store calc mode
    TYPE(L1OptType)                     :: Opt         ! Save Options
    INTEGER                             :: nBand       ! For convenience
    INTEGER,                ALLOCATABLE :: ncid(:)     ! NetCDF file index
    INTEGER,                ALLOCATABLE :: spcid(:)    ! Band Group Index
    INTEGER,                ALLOCATABLE :: auxid(:)    ! Support Group Index
    INTEGER,                ALLOCATABLE :: geoid(:)    ! Geolocation Group Index
    INTEGER                             :: ix          ! Current Cross Track Index
    INTEGER                             :: iy          ! Current Time index
    INTEGER                             :: imx         ! Number of cross track pixels
    INTEGER                             :: jmx         ! Number of along track pixels
    TYPE(L1BandType),       ALLOCATABLE :: Spectrum(:) ! Band spectrum
    TYPE(GeolocationType)               :: Geolocation
    TYPE(AuxiliaryL1Type)               :: Auxiliary
    TYPE(AuxProfType)                   :: Profile
    TYPE(AuxSurfType),      ALLOCATABLE :: Surface(:)
    LOGICAL,  ALLOCATABLE               :: BadQualityFlag(:,:)
  ENDTYPE L1Type
  
  TYPE L1DiagOptType
    LOGICAL :: ViewingGeometry
  ENDTYPE L1DiagOptType
  
  CONTAINS

  ! Figure out how to combine L1 into single object later
  TYPE(GeolocationType) FUNCTION GetBandGeolocation(Level1,BandIndex,Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(L1Type),    INTENT(IN)    :: Level1
    INTEGER,         INTENT(IN)    :: BandIndex
    TYPE(ErrorType), INTENT(INOUT) :: Error

    ! ---------------
    ! Local variables
    ! ---------------

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'GetBandGeolocation'

    ! =====================================================================
    ! GetBandGeolocation starts here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Initialize with the base geolocation
    GetBandGeolocation = Level1%Geolocation

    ! Check Viewing Geometry
    IF(Level1%Spectrum(BandIndex)%ReplaceGeolocationSZA) &
      GetBandGeolocation%SZA = Level1%Spectrum(BandIndex)%SZA
    IF(Level1%Spectrum(BandIndex)%ReplaceGeolocationVZA) &
      GetBandGeolocation%VZA = Level1%Spectrum(BandIndex)%VZA
    IF(Level1%Spectrum(BandIndex)%ReplaceGeolocationAZA) &
      GetBandGeolocation%AZA = Level1%Spectrum(BandIndex)%AZA
    
  END FUNCTION GetBandGeolocation

END MODULE level1_def