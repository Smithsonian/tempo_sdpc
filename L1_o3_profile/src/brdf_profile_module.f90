MODULE profile_module
  
  USE parameters_module
  USE covariance_module,    ONLY : UncertType, GetOutputCovarDim, ComputeZCorrCovar
  USE error_module,         ONLY : ErrorType, CheckError, RaiseFatalError, RaiseWarning
  USE level1_def,           ONLY : GeolocationType, L1Type, L1OptType, &
                                   AuxiliaryL1Type, AuxProfType
  USE interpolation_module, ONLY : BSPLINE, XYGridType,SPLINE1,SPLINT1,&
                                   NearestNeighbourSampler, XYGridWtType,&
                                   VerticalRegridWeights, BSPLINE_EdgeFill,&
                                   PointInPolygonSampler,LinearInt_Edgefill
  USE profilepar_module,    ONLY : profiles_uniform, profiles_expone, profiles_gdfone
  USE netcdf_module,        ONLY : NCDimType, match_names_in_dimlist, &
                                   nc_fld_1d, netcdf_handle_error, &
                                   CheckNetCDFErrorStatus, ncdf_var_exists

  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  
  TYPE ProfOptType
    CHARACTER(LEN=maxChar)              :: infile
    
    ! Sampling Option
    INTEGER                             :: SamplingIndex

    ! Gas Profile Options
    INTEGER                             :: nGas
    CHARACTER(LEN=maxChar), ALLOCATABLE :: GasName(:)
    REAL(KIND=8),           ALLOCATABLE :: GasMolecularWeight(:)
    

    ! Species to subtract from total air column
    INTEGER                             :: nWetSpc
    CHARACTER(LEN=maxChar), ALLOCATABLE :: WetSpcName(:)
    INTEGER,                ALLOCATABLE :: WetIdx(:)

    ! For computing the proxy column average mixing ratio
    CHARACTER(LEN=maxCHar)              :: ProxyNormName
    INTEGER                             :: ProxyNormIdx


    ! Aerosol Profile Options
    INTEGER                             :: nAer
    CHARACTER(LEN=maxChar), ALLOCATABLE :: AerName(:)
    REAL(KIND=8),           ALLOCATABLE :: AerMolecularWeight(:)
    REAL(KIND=8),           ALLOCATABLE :: AerDryDensity(:)
    REAL(KIND=8)                        :: AerRefWavelength
    LOGICAL                             :: UseFileAerosol
    REAL(KIND=8),           ALLOCATABLE :: AerColAOD(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: AerProfType(:)
    REAL(KIND=8),           ALLOCATABLE :: AerZmin(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerZmax(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerZpeak(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerZwdth(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AerZexp(:) ![km]
    
    ! Cloud Profile Options
    INTEGER                             :: nCld
    CHARACTER(LEN=maxChar), ALLOCATABLE :: CldName(:)
    LOGICAL                             :: UseFileCloud
    INTEGER                             :: nCldSubPixels
    LOGICAL                             :: DoLambertianCloud
    REAL(KIND=8),           ALLOCATABLE :: CldFraction(:)
    REAL(KIND=8),           ALLOCATABLE :: CldPressure(:) ![hPa]
    REAL(KIND=8)                        :: CldRefWavelength ![nm]
    REAL(KIND=8),           ALLOCATABLE :: CldColOD(:,:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: CldProfType(:)
    REAL(KIND=8),           ALLOCATABLE :: CldZmin(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldZmax(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldZpeak(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldZwdth(:,:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: CldZexp(:,:) ![km]
    
    ! Gravity Options
    REAL(KIND=8)                        :: StandardGravity
    REAL(KIND=8)                        :: PlanetaryRadius

    ! 
    LOGICAL                             :: LinearizeHybridGrid
  ENDTYPE ProfOptType
  
  TYPE GasProfileType
    REAL(KIND=8),           ALLOCATABLE :: PartialColumn(:,:) ![molec/cm2]
    TYPE(UncertType),       ALLOCATABLE :: Uncertainty(:)
    REAL(KIND=8),           ALLOCATABLE :: MixingRatio(:,:) ![v/v]
    REAL(KIND=8),           ALLOCATABLE :: BackgroundMixingRatio(:,:) ![v/v]
    REAL(KIND=8),           ALLOCATABLE :: MolecularWeight(:) ![kg/mol]
    CHARACTER(LEN=maxChar), ALLOCATABLE :: Name(:)
    INTEGER(KIND=2),        ALLOCATABLE :: SigmaAdjustType(:)
    INTEGER                             :: nSpecies
    LOGICAL,                ALLOCATABLE :: UseL2Gas(:)
    INTEGER,                ALLOCATABLE :: L2Idx(:)

    ! New Parameters for Dry Air Mixing Ratios
    INTEGER                             :: nWetAirSpc
    INTEGER,                ALLOCATABLE :: WetAirIdx(:)
    INTEGER                             :: ProxyNormIdx
    REAL(KIND=8)                        :: AprioriProxyNormMixingRatio

  ENDTYPE GasProfileType
  
  TYPE AerProfileType
    REAL(KIND=8),           ALLOCATABLE :: LayerOpticalDepth(:,:)
    REAL(KIND=8),           ALLOCATABLE :: ColumnOpticalDepth(:)
    REAL(KIND=8),           ALLOCATABLE :: ColumnOptDepthDeriv(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltMin(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AltMax(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AltPeak(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AltSigma(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AltExp(:) ![km]
    REAL(KIND=8),           ALLOCATABLE :: AltPeakDeriv(:,:) ![1/km]
    REAL(KIND=8),           ALLOCATABLE :: AltSigmaDeriv(:,:) ![1/km]
    REAL(KIND=8),           ALLOCATABLE :: AltExpDeriv(:,:) ![1/km]
    TYPE(UncertType),       ALLOCATABLE :: Uncertainty(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: Name(:)
    INTEGER(KIND=2),        ALLOCATABLE :: TypeIndex(:)
    INTEGER,                ALLOCATABLE :: nProfilePar(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: ProfileParName(:,:)
    INTEGER                             :: nSpecies
    LOGICAL                             :: UseFileAerosol
    LOGICAL,                ALLOCATABLE :: UseL2Aerosol(:)
    INTEGER,                ALLOCATABLE :: L2Idx(:)
  ENDTYPE AerProfileType
  
  TYPE CldProfileType
    LOGICAL                             :: DoLambertianCloud
    INTEGER                             :: nSubPix
    REAL(KIND=8)                        :: TotalCloudFraction
    REAL(KIND=8),           ALLOCATABLE :: CloudFraction(:)
    REAL(KIND=8),           ALLOCATABLE :: CloudPressure(:)
    INTEGER,                ALLOCATABLE :: LambertianCldLevel(:)
    REAL(KIND=8),           ALLOCATABLE :: LambertianCldLayerFrac(:)
    REAL(KIND=8),           ALLOCATABLE :: LayerOpticalDepth(:,:,:)
    REAL(KIND=8),           ALLOCATABLE :: ColumnOpticalDepth(:,:)
    REAL(KIND=8),           ALLOCATABLE :: ColumnOptDepthDeriv(:,:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltMin(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltMax(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltPeak(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltSigma(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltExp(:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltPeakDeriv(:,:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltSigmaDeriv(:,:,:)
    REAL(KIND=8),           ALLOCATABLE :: AltExpDeriv(:,:,:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: Name(:)
    INTEGER(KIND=2),        ALLOCATABLE :: TypeIndex(:)
    INTEGER,                ALLOCATABLE :: nProfilePar(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: ProfileParName(:,:)
    INTEGER                             :: nSpecies
    LOGICAL                             :: UseFileCloud
    LOGICAL,                ALLOCATABLE :: UseL2Cloud(:)
    INTEGER,                ALLOCATABLE :: L2Idx(:)
  ENDTYPE CldProfileType

  TYPE MetProfileType
    REAL(KIND=8), ALLOCATABLE :: PressureEdge(:)
    REAL(KIND=8), ALLOCATABLE :: PressureMid(:)
    REAL(KIND=8), ALLOCATABLE :: TemperatureEdge(:)
    REAL(KIND=8), ALLOCATABLE :: TemperatureMid(:)
    TYPE(UncertType)          :: TempUncertainty
    TYPE(UncertType)          :: TempShiftUncertainty
    TYPE(UncertType)          :: SurfPresUncertainty
    REAL(KIND=8), ALLOCATABLE :: AltitudeEdge(:)
    REAL(KIND=8), ALLOCATABLE :: AltitudeMid(:)
    REAL(KIND=8), ALLOCATABLE :: AirMolecularWeight(:)
    REAL(KIND=8), ALLOCATABLE :: Gravity(:)
    REAL(KIND=8), ALLOCATABLE :: RH(:)
    REAL(KIND=8), ALLOCATABLE :: AP(:)
    REAL(KIND=8), ALLOCATABLE :: BP(:)
    REAL(KIND=8), ALLOCATABLE :: AirPartialColumn(:)
    REAL(KIND=8), ALLOCATABLE :: DryAirPartialColumn(:)

    ! For CIA
    REAL(KIND=8), ALLOCATABLE :: AirPartialColumnSquared(:)

    ! INTEGER                   :: nFine ! Fine p/T For GasXS 
    ! REAL(KIND=8), ALLOCATABLE :: PressureFine(:)
    ! REAL(KIND=8), ALLOCATABLE :: TemperatureFine(:)

  ENDTYPE MetProfileType
  
  TYPE SurfProfType
    REAL(KIND=8)              :: WindSpeed
    REAL(KIND=8)              :: WindDirection
    REAL(KIND=8)              :: Chlorophyll
    REAL(KIND=8)              :: OceanSalinity
    REAL(KIND=8)              :: SnowFraction
    REAL(KIND=8)              :: SeaIceFraction
    REAL(KIND=8)              :: SnowDepth
    REAL(KIND=8)              :: SnowAge
    REAL(KIND=8)              :: LandCoverFraction(17)
    REAL(KIND=8)              :: SIF_734nm
  ENDTYPE SurfProfType
  
  TYPE ClimVGridType
    INTEGER                   :: lmx
    INTEGER(KIND=2)           :: GridTypeIndex
    REAL(KIND=8), ALLOCATABLE :: RegridWt(:,:)
    REAL(KIND=8), ALLOCATABLE :: RegridWtSigma(:,:) ! Compute regrid weights adjusting pres.
    REAL(KIND=8), ALLOCATABLE :: PressureEdge(:)
    REAL(KIND=8), ALLOCATABLE :: PressureMid(:)
    REAL(KIND=8), ALLOCATABLE :: PressureEdgeSigma(:) ! After converting to same grid
    REAL(KIND=8), ALLOCATABLE :: PressureMidSigma(:)
    REAL(KIND=8), ALLOCATABLE :: AP(:)
    REAL(KIND=8), ALLOCATABLE :: BP(:)

  END TYPE ClimVGridType

  TYPE ProfileType
    CHARACTER(LEN=maxChar) :: CalculationMode
    CHARACTER(LEN=maxChar) :: Infile
    INTEGER                :: SamplingIndex
    INTEGER                :: lmx ! # vertical layers
    INTEGER                :: ncid
    INTEGER(KIND=2)        :: GridTypeIndex
    TYPE(XYGridType)       :: XYGrid
    TYPE(MetProfileType)   :: Met
    TYPE(GasProfileType)   :: Gas
    TYPE(AerProfileType)   :: Aer
    TYPE(CldProfileType)   :: Cld
    TYPE(SurfProfType)     :: Surface
    REAL(KIND=8)           :: PlanetaryRadius
    REAL(KIND=8)           :: StandardGravity
    LOGICAL                :: LinearizeHybridGrid
    
    ! If climatology is being usurped
    LOGICAL                   :: UseL2Met
    TYPE(ClimVgridType)       :: ClimVertGrid
    
  ENDTYPE ProfileType

  ! --------------------------------------------
  ! Diagnostics 
  ! --------------------------------------------

  TYPE ProfDiagOptType
    
    ! Meteorological
    LOGICAL          :: PressureEdge
    LOGICAL          :: PressureMid
    LOGICAL          :: TemperatureEdge
    LOGICAL          :: TemperatureMid
    LOGICAL          :: AltitudeEdge
    LOGICAL          :: AltitudeMid
    LOGICAL          :: AirMolecularWeight
    LOGICAL          :: Gravity
    LOGICAL          :: RH
    LOGICAL          :: AirPartialColumn
    
    ! Gas
    TYPE(DiagSpcOpt) :: GasPartialColumn
    TYPE(DiagSpcOpt) :: GasUncertainty
    TYPE(DiagSpcOpt) :: GasMixingRatio
    TYPE(DiagSpcOpt) :: DryGasMixingRatio
    TYPE(DiagSpcOpt) :: ProxyColumnMixingRatio

    
    ! Aerosol
    LOGICAL          :: TotalAOD
    TYPE(DiagSpcOpt) :: LayerOpticalDepth
    TYPE(DiagSpcOpt) :: ProfilePar
    TYPE(DiagSpcOpt) :: ProfileParDeriv
    
  ENDTYPE ProfDiagOptType
  
  TYPE ProfDiagType
    TYPE(ProfDiagOptType)     :: Options
    INTEGER                   :: lmx
  ENDTYPE ProfDiagType

  PUBLIC :: InitProfile,                    &
            SampleProfile,                  &
            ComputeProfileDerivedQuantities
  
  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'profile_module'
  PRIVATE :: ModuleName

  CONTAINS
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: InitProfile
  ! 
  ! DESCRIPTION: Initializes profile:
  !                Determines vertical dimension
  !                Reads climatology geolocation
  !                Allocates arrays for profile variables
  !                Sets options passed from input
  !                Sets cloud/aersol profiles from input file if 
  !                option is specified
  
  SUBROUTINE InitProfile( ProfOpt, CalculationMode, L1Opt, L2Prof, Profile, Error )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ProfOptType),      INTENT(IN)    :: ProfOpt
    CHARACTER(LEN=maxChar), INTENT(IN)    :: CalculationMode
    TYPE(L1OptType),        INTENT(IN)    :: L1Opt
    TYPE(AuxProfType),      INTENT(IN)    :: L2Prof
    TYPE(ProfileType),      INTENT(INOUT) :: Profile
    TYPE(ErrorType),        INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER                :: rcode, vid, dimid(4), n, c, i, lmx
    INTEGER(KIND=2)        :: vgridtype
    CHARACTER(LEN=maxChar) :: tmpchar

    REAL(KIND=8), ALLOCATABLE :: ap(:), bp(:)
    !REAL(KIND=8)              :: zmin, zmax
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitProfile'

    ! ============================================================
    ! InitProfile starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Store the calculation mode
    Profile%CalculationMode = CalculationMode

    ! Save profile file name
    Profile%infile = ProfOpt%infile
    
    ! Save sampling Option
    Profile%SamplingIndex = ProfOpt%SamplingIndex
    
    ! Save the cloud option
    Profile%Cld%DoLambertianCloud = ProfOpt%DoLambertianCloud
    
    ! Are we using Level-2 Meteorology?
    Profile%UseL2Met = L1Opt%AuxProfile%OverwriteProfileClim
    
    ! Store the "wet" gas indices and proxy
    Profile%Gas%nWetAirSpc = ProfOpt%nWetSpc
    
    Profile%Gas%WetAirIdx = ProfOpt%WetIdx
    
    Profile%Gas%ProxyNormIdx = ProfOpt%ProxyNormIdx
    
    ! Attach file
    rcode = nf_open(TRIM(ADJUSTL(ProfOpt%infile)), nf_Share, Profile%ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

    ! Get the type index
    rcode = nf_inq_varid(Profile%ncid, 'VerticalGridType', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_get_vara_int2(Profile%ncid, vid,      &
                             (/1/), (/1/), vgridtype )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    
    ! Determine the vertical dimension based on profile type
    ! ------------------------------------------------------
    
    ! (1) Pressure Edges are explicitly specified
    IF( vgridtype .EQ. 1 ) THEN
      
      ! Use the pressure field to get the vertical dimension
      rcode = nf_inq_varid(Profile%ncid, 'pedge', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_inq_vardimid(Profile%ncid, vid, dimid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_inq_dim(Profile%ncid, dimid(3),tmpchar,lmx)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      lmx = lmx - 1

    ! (2) Eta Grid Parameterization
    ELSEIF( vgridtype .EQ. 2 ) THEN

      ! Use the Eta A coefficient to get vertical dimension
      rcode = nf_inq_varid(Profile%ncid, 'ap', vid) 
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_inq_vardimid(Profile%ncid, vid, dimid(1))
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_inq_dim(Profile%ncid, dimid(1),tmpchar,lmx) 
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      lmx = lmx - 1

    ELSE

      print*,'Vertical Grid type index in profile is unrecognized',&
             vgridtype
      STOP 1

    ENDIF
    
    ! Set dimensions based on whether Met is coming from climatology or L2 
    IF(Profile%UseL2Met) THEN

      Profile%lmx = L2Prof%lmx
      Profile%GridTypeIndex = L2Prof%GridTypeIndex
      Profile%ClimVertGrid%lmx = lmx
      Profile%ClimVertGrid%GridTypeIndex = vgridtype

      ! Allocate arrays for the Climatology grid
      ALLOCATE(Profile%ClimVertGrid%PressureMid(lmx))
      ALLOCATE(Profile%ClimVertGrid%PressureEdge(lmx+1))
      ALLOCATE(Profile%ClimVertGrid%AP(lmx+1))
      ALLOCATE(Profile%ClimVertGrid%BP(lmx+1))
      ALLOCATE(Profile%ClimVertGrid%RegridWt(lmx,L2Prof%lmx))

      ! On the adjusted grid (using the L2 surf. pres)
      ALLOCATE(Profile%ClimVertGrid%PressureMidSigma(lmx))
      ALLOCATE(Profile%ClimVertGrid%PressureEdgeSigma(lmx+1))
      ALLOCATE(Profile%ClimVertGrid%RegridWtSigma(lmx,L2Prof%lmx))

    ELSE
      Profile%lmx = lmx
      Profile%GridTypeIndex = vgridtype
    ENDIF

    ! Allocate temporary array for reading ap/bp
    ALLOCATE(ap(lmx+1),bp(lmx+1))
    
    ! Get the Lon-Lat-time grid dimensions
    ! ------------------------------------
    
    ! Longitude dimension
    rcode = nf_inq_varid(Profile%ncid, 'xmid', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_vardimid(Profile%ncid, vid, dimid(1))
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_dim(Profile%ncid, dimid(1),tmpchar,Profile%XYGrid%imx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

    ! Latitude dimension
    rcode = nf_inq_varid(Profile%ncid, 'ymid', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_vardimid(Profile%ncid, vid, dimid(1))
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_dim(Profile%ncid, dimid(1),tmpchar,Profile%XYGrid%jmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

    ! Time
    rcode = nf_inq_varid(Profile%ncid, 'time', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_vardimid(Profile%ncid, vid, dimid(1))
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_dim(Profile%ncid, dimid(1),tmpchar,Profile%XYGrid%tmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

    ! Maximum Cloud sub-pixels
    Profile%Cld%UseFileCloud = ProfOpt%UseFileCloud 
    
    ! -------------------------------------------------
    ! Allocate Arrays
    ! -------------------------------------------------
    
    ! Geolocation
    ALLOCATE(Profile%XYGrid%LongitudeMid(Profile%XYGrid%imx))
    ALLOCATE(Profile%XYGrid%LongitudeEdge(Profile%XYGrid%imx+1))
    ALLOCATE(Profile%XYGrid%LatitudeMid(Profile%XYGrid%jmx))
    ALLOCATE(Profile%XYGrid%LatitudeEdge(Profile%XYGrid%jmx+1))
    ALLOCATE(Profile%XYGrid%TauTime(Profile%XYGrid%tmx))

    ! Meteorological 
    ALLOCATE(Profile%Met%PressureEdge(Profile%lmx+1))
    ALLOCATE(Profile%Met%PressureMid(Profile%lmx))
    ALLOCATE(Profile%Met%AltitudeEdge(Profile%lmx+1))
    ALLOCATE(Profile%Met%AltitudeMid(Profile%lmx))
    ALLOCATE(Profile%Met%TemperatureEdge(Profile%lmx+1))
    ALLOCATE(Profile%Met%TemperatureMid(Profile%lmx))
    ALLOCATE(Profile%Met%RH(Profile%lmx))
    ALLOCATE(Profile%Met%AirMolecularWeight(Profile%lmx))
    ALLOCATE(Profile%Met%Gravity(Profile%lmx))
    ALLOCATE(Profile%Met%AirPartialColumn(Profile%lmx))
    ALLOCATE(Profile%Met%AirPartialColumnSquared(Profile%lmx))
    ALLOCATE(Profile%Met%DryAirPartialColumn(Profile%lmx))
    ALLOCATE(Profile%Met%AP(Profile%lmx+1))
    ALLOCATE(Profile%Met%BP(Profile%lmx+1))
    
    ! Gas
    Profile%Gas%nSpecies = ProfOpt%nGas
    ALLOCATE(Profile%Gas%PartialColumn(Profile%lmx,Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%Uncertainty(Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%MixingRatio(Profile%lmx,Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%BackgroundMixingRatio(Profile%lmx,Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%SigmaAdjustType(Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%MolecularWeight(Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%Name(Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%UseL2Gas(Profile%Gas%nSpecies))
    ALLOCATE(Profile%Gas%L2Idx(Profile%Gas%nSpecies))

    ! Aerosol
    Profile%Aer%nSpecies = ProfOpt%nAer
    ALLOCATE(Profile%Aer%LayerOpticalDepth(Profile%lmx,Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%ColumnOpticalDepth(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%ColumnOptDepthDeriv(Profile%lmx,Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%Uncertainty(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltMin(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltMax(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltPeak(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltSigma(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltExp(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%Name(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%TypeIndex(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%nProfilePar(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%ProfileParName(Profile%Aer%nSpecies,MaxProfPar))
    ALLOCATE(Profile%Aer%AltPeakDeriv(Profile%lmx,Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltSigmaDeriv(Profile%lmx,Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%AltExpDeriv(Profile%lmx,Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%UseL2Aerosol(Profile%Aer%nSpecies))
    ALLOCATE(Profile%Aer%L2Idx(Profile%Aer%nSpecies))

    ! Read Horizontal Grid
    rcode = nf_inq_varid(Profile%ncid, 'xmid', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_get_vara_double( Profile%ncid, vid,            &
                                (/1/), (/Profile%XYGrid%imx/),&
                                Profile%XYGrid%LongitudeMid(:))
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_varid(Profile%ncid, 'ymid', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_get_vara_double( Profile%ncid, vid,            &
                                (/1/), (/Profile%XYGrid%jmx/),&
                                Profile%XYGrid%LatitudeMid(:) )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_varid(Profile%ncid, 'xedge', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_get_vara_double( Profile%ncid, vid,              &
                                (/1/), (/Profile%XYGrid%imx+1/),&
                                Profile%XYGrid%LongitudeEdge(:) )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_varid(Profile%ncid, 'yedge', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_get_vara_double( Profile%ncid, vid,              &
                                (/1/), (/Profile%XYGrid%jmx+1/),&
                                Profile%XYGrid%LatitudeEdge(:)  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_inq_varid(Profile%ncid, 'time', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    rcode = nf_get_vara_double( Profile%ncid, vid,            &
                                (/1/), (/Profile%XYGrid%tmx/),&
                                Profile%XYGrid%TauTime(:)     )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    
    ! Read Eta Grid Coefficients if they are defined
    IF( vgridtype .EQ. 2) THEN
      rcode = nf_inq_varid(Profile%ncid, 'ap', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_get_vara_double( Profile%ncid, vid,      &
                                 (/1/), (/lmx+1/), ap)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_inq_varid(Profile%ncid, 'bp', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      rcode = nf_get_vara_double( Profile%ncid, vid,      &
                                 (/1/), (/lmx+1/), bp)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

      ! Set the coefficients
      IF(Profile%UseL2Met) THEN
        Profile%ClimVertGrid%AP = ap
        Profile%ClimVertGrid%AP = bp
      ELSE
        Profile%Met%AP(:) = ap
        Profile%Met%BP(:) = bp
      ENDIF

    ENDIF
    
    ! AP/BP Coefficients in L2 Profile set if needed
    IF(Profile%UseL2Met) THEN
      Profile%Met%AP(:) = L2Prof%AP
      Profile%Met%BP(:) = L2Prof%BP
    ENDIF
    
    ! Initialize flags for using L2
    Profile%Gas%UseL2Gas(:) = .FALSE. ; Profile%Gas%L2Idx(:) = -1

    ! Transfer input option names
    DO n=1,Profile%Gas%nSpecies
      
      ! Get Gas name and molecular weight
      Profile%Gas%Name(n) = ProfOpt%GasName(n)
      Profile%Gas%MolecularWeight(n) = ProfOpt%GasMolecularWeight(n)*1e-3 ! g/mol->kg/mol
      
      ! Read the pressure adjustment type for each gas
      rcode = nf_inq_varid(Profile%ncid,TRIM(ADJUSTL(Profile%Gas%Name(n)))//'_SigmaAdjustType',vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,            &
             'nf_inq_varid:'//TRIM(ADJUSTL(Profile%Gas%Name(n)))//'_SigmaAdjustType')
      rcode = nf_get_vara_int2(Profile%ncid, vid,(/1/), (/1/),&
                               Profile%Gas%SigmaAdjustType(n) )
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,                &
             'nf_get_vara_int2:'//TRIM(ADJUSTL(Profile%Gas%Name(n)))//'_SigmaAdjustType')
      
      ! Check if we are overwriting the field
      DO i=1,L1Opt%AuxProfile%nGas
        IF(TRIM(ADJUSTL(Profile%Gas%Name(n))) .EQ.    &
           TRIM(ADJUSTL(L1Opt%AuxProfile%GasName(i))) ) THEN
          Profile%Gas%UseL2Gas(n) = .TRUE.
          Profile%Gas%L2Idx(n)    = i
          Profile%Gas%SigmaAdjustType = L2Prof%GasSigmaAdjustType(i)
        ENDIF
      ENDDO

    ENDDO

    ! Store 
    Profile%Aer%UseFileAerosol = ProfOpt%UseFileAerosol
    
    ! Initialize the flags for using L2 
    Profile%Aer%UseL2Aerosol(:) = .FALSE. ; Profile%Aer%L2Idx(:) = -1

    ! Transfer aerosol input options
    DO n=1,Profile%Aer%nSpecies
      
      ! Set Aerosol name
      Profile%Aer%Name(n) = ProfOpt%AerName(n)
      
      IF( Profile%Aer%UseFileAerosol ) THEN
        
        ! Check if aerosol field is present in file
        DO i=1,L1Opt%AuxProfile%nAer
          IF(TRIM(ADJUSTL(Profile%Aer%Name(n))) .EQ.    &
            TRIM(ADJUSTL(L1Opt%AuxProfile%AerName(i))) ) THEN

            Profile%Aer%UseL2Aerosol(n) = .TRUE.
            Profile%Aer%L2Idx(n)    = i

          ENDIF
        ENDDO

        ! Set field based on input
        IF(Profile%Aer%UseL2Aerosol(n)) THEN

          ! Copy values from L2 Auxiliary data
          Profile%Aer%TypeIndex(n)   = L2Prof%AerTypeIndex(Profile%Aer%L2Idx(n))
          Profile%Aer%nProfilePar(n) = L2Prof%AernProfilePar(Profile%Aer%L2Idx(n))
          DO i=1,Profile%Aer%nProfilePar(n) 
            Profile%Aer%ProfileParName(n,i) = L2Prof%AerProfileParName(Profile%Aer%L2Idx(n),i)
          ENDDO

        ELSE
          
          ! Load Profile type index
          rcode = nf_inq_varid(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))) // '_StateType', vid)
          CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,     &
              'nf_inq_varid:'//TRIM(ADJUSTL(Profile%Aer%Name(n))) //'_StateType')
          
          ! Read
          rcode = nf_get_vara_int2(Profile%ncid, vid,     &
                                  (/1/), (/1/),           &
                                  Profile%Aer%TypeIndex(n))
          CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,         &
              'nf_get_vara_int2:'//TRIM(ADJUSTL(Profile%Aer%Name(n))) //'_StateType')
          
          ! Set profile parameters based on climatology profile index
          IF( Profile%Aer%TypeIndex(n) .EQ. 1 ) THEN
            ! Linearization parameters
            Profile%Aer%nProfilePar(n) = 1
            Profile%Aer%ProfileParName(n,1) = 'AerODCol'
          
          ELSEIF( Profile%Aer%TypeIndex(n) .EQ. 2 ) THEN

            ! Linearization parameters
            Profile%Aer%nProfilePar(n) = 3
            Profile%Aer%ProfileParName(n,1) = 'AerODCol'
            Profile%Aer%ProfileParName(n,2) = 'PeakHeight'
            Profile%Aer%ProfileParName(n,3) = 'PeakWidth'

          ELSEIF(Profile%Aer%TypeIndex(n) .EQ. 3) THEN

            ! Linearization parameters
            Profile%Aer%nProfilePar(n) = 2
            Profile%Aer%ProfileParName(n,1) = 'AerODCol'
            Profile%Aer%ProfileParName(n,2) = 'RelaxHeight'

          ELSEIF(Profile%Aer%TypeIndex(n) .EQ. 4) THEN
            
            ! Linearization parameters
            Profile%Aer%nProfilePar(n) = 1
            Profile%Aer%ProfileParName(n,1) = 'AerODCol'
          ELSE
            STOP 'Incorrect aerosol type index in climatology file'
          ENDIF

        ENDIF

      ELSE
        
        ! Set profile parameters based on input file
        IF( TRIM(ADJUSTL(ProfOpt%AerProfType(n))) .EQ. 'GDF' ) THEN
          Profile%Aer%TypeIndex(n) = 2

          ! Linearization parameters
          Profile%Aer%nProfilePar(n) = 3
          Profile%Aer%ProfileParName(n,1) = 'AerODCol'
          Profile%Aer%ProfileParName(n,2) = 'PeakHeight'
          Profile%Aer%ProfileParName(n,3) = 'PeakWidth'

        ELSEIF(TRIM(ADJUSTL(ProfOpt%AerProfType(n))) .EQ. 'EXP') THEN
          Profile%Aer%TypeIndex(n) = 3

          ! Linearization parameters
          Profile%Aer%nProfilePar(n) = 2
          Profile%Aer%ProfileParName(n,1) = 'AerODCol'
          Profile%Aer%ProfileParName(n,2) = 'RelaxHeight'

        ELSEIF(TRIM(ADJUSTL(ProfOpt%AerProfType(n))) .EQ. 'BOX') THEN
          Profile%Aer%TypeIndex(n) = 4
          Profile%Aer%nProfilePar(n) = 1
          Profile%Aer%ProfileParName(n,1) = 'AerODCol'
        ELSE
          STOP 'Incorrect profile parameter specified in input file (allowed options GDF,EXP, and BOX)'
        ENDIF
        
        ! Store parameters from input file
        Profile%Aer%ColumnOpticalDepth(n) = ProfOpt%AerColAOD(n)
        Profile%Aer%AltMin(n)             = ProfOpt%AerZmin(n)
        Profile%Aer%AltMax(n)             = ProfOpt%AerZmax(n)
        Profile%Aer%AltPeak(n)            = ProfOpt%AerZpeak(n)
        Profile%Aer%AltSigma(n)           = ProfOpt%AerZwdth(n)
        Profile%Aer%AltExp(n)             = ProfOpt%AerZexp(n)
        
      ENDIF
    ENDDO
    
    ! Cloud - subpixel-independent arrays
    IF(L1Opt%AuxCloud%UseL2MetCloud) THEN
      Profile%Cld%nSpecies = L1Opt%AuxProfile%nCld
    ELSE
      Profile%Cld%nSpecies = ProfOpt%nCld
    ENDIF
    
    ALLOCATE(Profile%Cld%Name(Profile%Cld%nSpecies))
    ALLOCATE(Profile%Cld%TypeIndex(Profile%Cld%nSpecies))
    ALLOCATE(Profile%Cld%nProfilePar(Profile%Cld%nSpecies))
    ALLOCATE(Profile%Cld%ProfileParName(Profile%Cld%nSpecies,MaxProfPar))

    ! Copy cloud options
    Profile%Cld%UseFileCloud = ProfOpt%UseFileCloud

    ! ------------------------------------------------------------------
    ! (1) Lambertian style cloud from auxiliary information
    ! ------------------------------------------------------------------
    IF(L1Opt%AuxCloud%AuxLambertian%OverwriteProfileClim .AND.  &
       Profile%Cld%DoLambertianCloud                            ) THEN
      
      ! Allocate arrays for single cloudy pixel if we are using L2 product
      ALLOCATE(Profile%Cld%CloudFraction(1))
      ALLOCATE(Profile%Cld%CloudPressure(1))
      ALLOCATE(Profile%Cld%LambertianCldLevel(1))
      ALLOCATE(Profile%Cld%LambertianCldLayerFrac(1))
    
    ! ----------------------------------------------------------------------
    ! (2) Option where we are directly replacing clouds from L2 sampled file
    ! ----------------------------------------------------------------------
    ELSEIF(L1Opt%AuxCloud%UseL2MetCloud) THEN
      
      ! Species names are from the level1 options
      DO n=1,L1Opt%AuxProfile%nCld
        Profile%Cld%Name(n) = L1Opt%AuxProfile%CldName(n)
      ENDDO

      ! Currently all cloud inputs are explicit profiles
      Profile%Cld%TypeIndex(:)        = 1
      Profile%Cld%nProfilePar(:)      = 1
      Profile%Cld%ProfileParName(:,1) = 'CldODCol'

    ! ----------------------------------------------------------------------------
    ! (3) Option Where we are taking clouds from the climatology (lambertian opt. doesnt req any linearization)
    ! ----------------------------------------------------------------------------
    ELSEIF(Profile%Cld%UseFileCloud .AND. .NOT. Profile%Cld%DoLambertianCloud) THEN
      
      DO n=1,ProfOpt%nCld

        ! Cloud Names are from the input file profile options
        Profile%Cld%Name(n) = ProfOpt%CldName(n)

        ! Load Profile type index
        rcode = nf_inq_varid(Profile%ncid, TRIM(ADJUSTL(Profile%Cld%Name(n))) // '_StateType', vid)
        CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,    &
             'nf_inq_varid:'//TRIM(ADJUSTL(Profile%Cld%Name(n))) // '_StateType')
          
        ! Read
        rcode = nf_get_vara_int2(Profile%ncid, vid,      &
                                (/1/), (/1/),           &
                                Profile%Cld%TypeIndex(n))
        CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,    &
           'nf_get_vara_int2:'//TRIM(ADJUSTL(Profile%Cld%Name(n))) // '_StateType')
          
        ! Set profile parameters based on climatology profile index
        IF( Profile%Cld%TypeIndex(n) .EQ. 1 ) THEN

          ! Linearization parameters
          Profile%Cld%nProfilePar(n) = 1
          Profile%Cld%ProfileParName(n,1) = 'CldODCol'
          
        ELSEIF( Profile%Cld%TypeIndex(n) .EQ. 2 ) THEN

          ! Linearization parameters
          Profile%Cld%nProfilePar(n) = 3
          Profile%Cld%ProfileParName(n,1) = 'CldODCol'
          Profile%Cld%ProfileParName(n,2) = 'PeakHeight'
          Profile%Cld%ProfileParName(n,3) = 'PeakWidth'

        ELSEIF(Profile%Cld%TypeIndex(n) .EQ. 3) THEN

          ! Linearization parameters
          Profile%Cld%nProfilePar(n) = 2
          Profile%Cld%ProfileParName(n,1) = 'CldODCol'
          Profile%Cld%ProfileParName(n,2) = 'RelaxHeight'

        ELSEIF(Profile%Cld%TypeIndex(n) .EQ. 4) THEN
            
          ! Linearization parameters
          Profile%Cld%nProfilePar(n) = 1
          Profile%Cld%ProfileParName(n,1) = 'CldODCol'
        ELSE
          STOP 'Incorrect cloud type index in climatology file'
        ENDIF

      ENDDO
    
    ! ------------------------------------------------------------------
    ! (4) Option Where we are taking clouds from the input file
    ! ------------------------------------------------------------------
    ELSE
        
      ! Set subpixel number
      Profile%Cld%nSubPix  = ProfOpt%nCldSubPixels
        
      ! Set Cloud Fraction
      ALLOCATE(Profile%Cld%CloudFraction(Profile%Cld%nSubPix))
      Profile%Cld%TotalCloudFraction = 0.0d0
      DO c=1,Profile%Cld%nSubPix
        Profile%Cld%CloudFraction(c) = ProfOpt%CldFraction(c)
        Profile%Cld%TotalCloudFraction = Profile%Cld%TotalCloudFraction &
                                       + ProfOpt%CldFraction(c)
      ENDDO
        
      IF(Profile%Cld%DoLambertianCloud) THEN

        ! Arrays needed for lambertian clouds
        ALLOCATE(Profile%Cld%CloudPressure(Profile%Cld%nSubPix))
        ALLOCATE(Profile%Cld%LambertianCldLevel(Profile%Cld%nSubPix))
        ALLOCATE(Profile%Cld%LambertianCldLayerFrac(Profile%Cld%nSubPix))

        ! Store Cloud pressures
        DO c=1,Profile%Cld%nSubPix
          Profile%Cld%CloudPressure(c) = ProfOpt%CldPressure(c)
        ENDDO

      ELSE
        
        ! Cloud Names are from input file
        DO n=1,ProfOpt%nCld
          Profile%Cld%Name(n) = ProfOpt%CldName(n)
        ENDDO
        
        ! Allocate profile parameter arrays
        ALLOCATE(Profile%Cld%LayerOpticalDepth(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%ColumnOpticalDepth(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%ColumnOptDepthDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltMin(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltMax(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltPeak(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltSigma(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltExp(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltPeakDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltSigmaDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
        ALLOCATE(Profile%Cld%AltExpDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))

        ! Zero Arrays
        Profile%Cld%LayerOpticalDepth(:,:,:) = 0.0d0
        Profile%Cld%ColumnOpticalDepth(:,:) = 0.0d0
        Profile%Cld%ColumnOptDepthDeriv(:,:,:) = 0.0d0
        Profile%Cld%AltMin(:,:) = 0.0d0
        Profile%Cld%AltMax(:,:) = 0.0d0
        Profile%Cld%AltPeak(:,:) = 0.0d0
        Profile%Cld%AltSigma(:,:) = 0.0d0
        Profile%Cld%AltExp(:,:) = 0.0d0
        Profile%Cld%AltPeakDeriv(:,:,:) = 0.0d0
        Profile%Cld%AltSigmaDeriv(:,:,:) = 0.0d0
        Profile%Cld%AltExpDeriv(:,:,:) = 0.0d0

        ! Store parameters from input file
        DO n=1,Profile%Cld%nSpecies

          ! Set values from input file
          DO c=1,Profile%Cld%nSubPix
            Profile%Cld%ColumnOpticalDepth(c,n) = ProfOpt%CldColOD(n,c)
            Profile%Cld%AltMin(c,n)             = ProfOpt%CldZmin(n,c)
            Profile%Cld%AltMax(c,n)             = ProfOpt%CldZmax(n,c)
            Profile%Cld%AltPeak(c,n)            = ProfOpt%CldZpeak(n,c)
            Profile%Cld%AltSigma(c,n)           = ProfOpt%CldZwdth(n,c)
            Profile%Cld%AltExp(c,n)             = ProfOpt%CldZexp(n,c)
          ENDDO

          ! Set profile parameters based on input file
          IF( TRIM(ADJUSTL(ProfOpt%CldProfType(n))) .EQ. 'GDF' ) THEN
            Profile%Cld%TypeIndex(n) = 2

            ! Linearization parameters
            Profile%Cld%nProfilePar(n) = 3
            Profile%Cld%ProfileParName(n,1) = 'CldODCol'
            Profile%Cld%ProfileParName(n,2) = 'PeakHeight'
            Profile%Cld%ProfileParName(n,3) = 'PeakWidth'

          ELSEIF(TRIM(ADJUSTL(ProfOpt%CldProfType(n))) .EQ. 'EXP') THEN
            Profile%Cld%TypeIndex(n) = 3

            ! Linearization parameters
            Profile%Cld%nProfilePar(n) = 2
            Profile%Cld%ProfileParName(n,1) = 'CldODCol'
            Profile%Cld%ProfileParName(n,2) = 'RelaxHeight'

          ELSEIF(TRIM(ADJUSTL(ProfOpt%CldProfType(n))) .EQ. 'BOX') THEN
            Profile%Cld%TypeIndex(n) = 4
            Profile%Cld%nProfilePar(n) = 1
            Profile%Cld%ProfileParName(n,1) = 'CldODCol'
          ELSE
            STOP 'Incorrect profile parameter specified in input file (allowed options GDF,EXP, and BOX)'
          ENDIF

        ENDDO
          
        
      ENDIF

    ENDIF
    
    ! Gravity
    Profile%PlanetaryRadius = ProfOpt%PlanetaryRadius
    Profile%StandardGravity = ProfOpt%StandardGravity
    Profile%LinearizeHybridGrid = ProfOpt%LinearizeHybridGrid
    
    ! Placeholder values for surface quantities
    Profile%Surface%LandCoverFraction(:) = 0.0 ; Profile%Surface%LandCoverFraction(17) = 1.0 ! Barren land
    
  END SUBROUTINE InitProfile
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: SampleProfile
  ! 
  ! DESCRIPTION: Samples a profile based on supplied Level1 Geolocation
  
  SUBROUTINE SampleProfile(Level1,Profile,L1Opt,L2Prof,Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    
    TYPE(L1Type),          INTENT(IN)    :: Level1
    TYPE(ProfileType),     INTENT(INOUT) :: Profile
    TYPE(L1OptType),       INTENT(IN)    :: L1Opt
    TYPE(AuxProfType),     INTENT(IN)    :: L2Prof
    TYPE(ErrorType),       INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER                   :: rcode, vid, dimid, n, c, errstat, g
    CHARACTER(LEN=maxChar)    :: tmpchar, message, action
    TYPE(XYGridWtType)        :: GridWt
    REAL(KIND=8)              :: psurf
    LOGICAL                   :: do_derivatives, fail
    REAL(KIND=8)              :: zmin, zmax
    REAL(KIND=8)              :: aer_zedge(0:Profile%lmx)

    INTEGER                   :: lmx, gidx, l, s
    REAL(KIND=8), ALLOCATABLE :: tmpmid(:),tmpedge(:), ap(:), bp(:)
    INTEGER(KIND=2)           :: ClimGridTypeIndex
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SampleProfile'

    ! ============================================================
    ! SampleProfile starts here
    ! ============================================================

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Determine the sampling
    IF(Profile%SamplingIndex .EQ. 1) THEN
      CALL NearestNeighbourSampler(Profile%XYGrid, Level1%Geolocation, GridWt, Error)
    ELSEIF(Profile%SamplingIndex .EQ. 2) THEN
      CALL PointInPolygonSampler(Profile%XYGrid, Level1%Geolocation, GridWt, Error)
    ELSE
      CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,&
                            Message_in='Urecognized Profile sampling type',      &
                            Action_in='Change index in input file (Options=1..2)')
    ENDIF
    
    ! Allocate temporary arrays
    IF(Profile%UseL2Met) THEN
      lmx = Profile%ClimVertGrid%lmx
      ClimGridTypeIndex = Profile%ClimVertGrid%GridTypeIndex
    ELSE
      lmx = Profile%lmx
      ClimGridTypeIndex = Profile%GridTypeIndex
    ENDIF
    ALLOCATE(tmpmid(lmx),tmpedge(lmx+1))
    
    IF( ClimGridTypeIndex .EQ. 1 ) THEN

      ! Load Edge pressures
      CALL LoadProfileVar(Profile%ncid, 'pedge', GridWt,&
                          tmpedge, lmx+1,Error          )
      
    ELSEIF( ClimGridTypeIndex .EQ. 2 ) THEN

      ! Load Edge pressures
      CALL LoadSurfaceVar(Profile%ncid,'psurf', GridWt,&
                          psurf,Error                  )
      
    ELSE

      STOP 'Unrecognized vertical grid type index'
      
    ENDIF
    
    IF(Profile%UseL2Met) THEN

      ! Set the Profile Grid Based on the L2 Met
      Profile%Met%PressureEdge = L2Prof%PressureEdge
      Profile%Met%AP = L2Prof%AP
      Profile%Met%BP = L2Prof%BP

      ! Compute climatology grid pressure
      IF(ClimGridTypeIndex .EQ. 1) THEN
        Profile%ClimVertGrid%PressureEdge = tmpedge
        Profile%ClimVertGrid%AP = 0.0d0
        Profile%ClimVertGrid%BP = tmpedge/tmpedge(lmx+1)
      ELSEIF(ClimGridTypeIndex .EQ. 2) THEN
        Profile%ClimVertGrid%PressureEdge(:) = Profile%ClimVertGrid%AP(:) &
                                             + Profile%ClimVertGrid%BP(:)*psurf
      ENDIF
      Profile%ClimVertGrid%PressureMid(:) = 0.5d0*(                        &
           Profile%ClimVertGrid%PressureEdge(1:Profile%ClimVertGrid%lmx)   &
          +Profile%ClimVertGrid%PressureEdge(2:Profile%ClimVertGrid%lmx+1) )

      ! Compute the pressure grid using the L2 definition
      Profile%ClimVertGrid%PressureEdgeSigma(:) = Profile%ClimVertGrid%AP(:) &
          + Profile%ClimVertGrid%BP(:)*Profile%Met%PressureEdge(Profile%lmx+1)
      Profile%ClimVertGrid%PressureMidSigma(:) = 0.5d0*(                        &
           Profile%ClimVertGrid%PressureEdgeSigma(1:Profile%ClimVertGrid%lmx)   &
          +Profile%ClimVertGrid%PressureEdgeSigma(2:Profile%ClimVertGrid%lmx+1) )
      
      ! Now Compute the regrid Weights
      CALL VerticalRegridWeights( Profile%ClimVertGrid%lmx,              &
                                  Profile%ClimVertGrid%PressureEdge,     &
                                  Profile%lmx, Profile%Met%PressureEdge, &
                                  Profile%ClimVertGrid%RegridWt          )
      CALL VerticalRegridWeights( Profile%ClimVertGrid%lmx,              &
                                  Profile%ClimVertGrid%PressureEdgeSigma,&
                                  Profile%lmx, Profile%Met%PressureEdge, &
                                  Profile%ClimVertGrid%RegridWtSigma     )
    ELSE

      ! Set the simulated profile from the climatology
      IF(ClimGridTypeIndex .EQ. 1) THEN
        Profile%Met%PressureEdge = tmpedge
        Profile%Met%AP = 0.0d0
        Profile%Met%BP = tmpedge/tmpedge(lmx+1)
      ELSEIF(ClimGridTypeIndex .EQ. 2) THEN
        Profile%Met%PressureEdge(:) = Profile%Met%AP(:) &
                                    + Profile%Met%BP(:)*psurf
      ENDIF

    ENDIF
    
    ! Midpoint pressure
    DO n=1,Profile%lmx
      Profile%Met%PressureMid(n) = 0.5d0*(Profile%Met%PressureEdge(n)  &
                                         +Profile%Met%PressureEdge(n+1))
    ENDDO
    
    ! Load Edge temperatures
    IF(Profile%UseL2Met) THEN
      Profile%Met%TemperatureEdge(:) = L2Prof%TemperatureEdge(:)
    ELSE
      CALL LoadProfileVar(Profile%ncid, 'Tedge', GridWt,                  &
                          Profile%Met%TemperatureEdge, Profile%lmx+1,Error)
    ENDIF

    ! Interpolate temperature to midpoint pressure
    CALL BSPLINE( Profile%Met%PressureEdge(:), Profile%Met%TemperatureEdge(:), Profile%lmx+1,      &
                  Profile%Met%PressureMid(:),  Profile%Met%TemperatureMid(:),  Profile%lmx, errstat)
    
    ! -----------------
    ! Relative Humidity
    ! -----------------

    ! Set from L2 if using Meteorological profile
    IF(Profile%UseL2Met) THEN
      Profile%Met%RH(:) = L2Prof%RH(:)
      
    ELSE

      ! Set RH if present
      IF(ncdf_var_exists(Profile%ncid,'RH')) THEN
        CALL LoadProfileVar(Profile%ncid, 'RH', GridWt,      &
                            Profile%Met%RH, Profile%lmx,Error)
        
      ! Set default if field not in climatology
      ELSE
        Profile%Met%RH = 0.0d0
        
      ENDIF


    ENDIF
    
    ! ----------------
    ! Surface Altitude
    ! ----------------
    
    !  Get surface altitude
    Profile%Met%AltitudeEdge(:) = 0.0d0
    
    ! Load surface pressure from climatology
    CALL LoadXYVar(Profile%ncid, 'zsurf', GridWt,               &
                   Profile%Met%AltitudeEdge(Profile%lmx+1),Error)
    
    ! -----------
    ! Gas Species
    ! -----------
    DO n=1,Profile%Gas%nSpecies
      
      IF(Profile%Gas%UseL2Gas(n)) THEN

        ! Get index
        gidx = Profile%Gas%L2Idx(n)

        ! Set the profile from the L2 Auxiliary data
        Profile%Gas%MixingRatio(:,n) = L2Prof%GasMixingRatio(:,gidx)
        IF(Profile%Gas%SigmaAdjustType(n) .EQ. 3) THEN
          Profile%Gas%BackgroundMixingRatio(:,n) = &
            L2Prof%GasBackgroundMixingRatio(:,n)
        ENDIF

      ELSE

        IF(Profile%UseL2Met) THEN
          
          ! The data must be regridded (pass the regrid weights)
          CALL LoadProfileVar(Profile%ncid, TRIM(ADJUSTL(Profile%Gas%Name(n))), GridWt,&
                              Profile%Gas%MixingRatio(:,n), Profile%lmx,Error,         &
                              ClimVertGrid=Profile%ClimVertGrid                        )
          
          IF(Profile%Gas%SigmaAdjustType(n) .EQ. 3) THEN
            CALL LoadProfileVar(Profile%ncid, TRIM(ADJUSTL(Profile%Gas%Name(n)))//'_Background',&
                                GridWt,Profile%Gas%BackgroundMixingRatio(:,n),Profile%lmx,Error,&
                                ClimVertGrid=Profile%ClimVertGrid                               )
          ENDIF
          
        ELSE

          ! Directly load profile
          CALL LoadProfileVar(Profile%ncid, TRIM(ADJUSTL(Profile%Gas%Name(n))), GridWt,&
                              Profile%Gas%MixingRatio(:,n), Profile%lmx,Error)
          
          ! Load backround Profile if using the hybrid grid adjustment method
          IF(Profile%Gas%SigmaAdjustType(n) .EQ. 3) THEN
            CALL LoadProfileVar(Profile%ncid, TRIM(ADJUSTL(Profile%Gas%Name(n)))//'_Background',&
                                GridWt,Profile%Gas%BackgroundMixingRatio(:,n), Profile%lmx,Error)
          ENDIF
        ENDIF

      ENDIF
      
    ENDDO
    
    ! Compute derived quantities
    CALL ComputeProfileDerivedQuantities(Profile, Error)
    
    ! ==========================
    ! Aerosol Profiles
    ! ==========================
    CALL SetAerosolProfile(GridWt,Profile,L1Opt,L2Prof,Error)
    
    ! ==========================
    ! Cloud Profiles
    ! ==========================

    ! ------------------------------------------------------------------------
    ! (1) Auxiliary Lambertian Clouds
    ! ------------------------------------------------------------------------
    IF( Level1%Opt%AuxCloud%AuxLambertian%OverwriteProfileClim .AND. &
        Profile%Cld%DoLambertianCloud                                ) THEN
      
      IF( Level1%Auxiliary%CloudFraction .GT. TINY(0.0d0) ) THEN
        Profile%Cld%nSubPix = 1 
        Profile%Cld%TotalCloudFraction = Level1%Auxiliary%CloudFraction
        Profile%Cld%CloudFraction(1)   = Level1%Auxiliary%CloudFraction
        Profile%Cld%CloudPressure(1)   = Level1%Auxiliary%CloudPressure
      ELSE
        Profile%Cld%nSubPix = 0
        Profile%Cld%TotalCloudFraction = 0.0D0
        Profile%Cld%CloudFraction(1)   = 0.0D0
        Profile%Cld%CloudPressure(1)   = 0.0D0
      ENDIF

    ! ----------------------------------------------------------------------
    ! (2) Option where we are directly replacing clouds from L2 sampled file
    ! ----------------------------------------------------------------------
    ELSEIF(L1Opt%AuxCloud%UseL2MetCloud) THEN

      
      ! (*) Lambertian Case
      IF(Profile%Cld%DoLambertianCloud) THEN

        ! Check option exists
        IF(.NOT. L2Prof%LambertianCloudsPresent) THEN
          CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,  &
                 Message_in='L2 Lambertian clouds requested but not in profile dataset')
        ENDIF

        ! Reallocate field
        Profile%Cld%nSubPix = L2Prof%nCldPix
        IF(ALLOCATED(Profile%Cld%CloudFraction)) DEALLOCATE(Profile%Cld%CloudFraction)
        IF(ALLOCATED(Profile%Cld%CloudPressure)) DEALLOCATE(Profile%Cld%CloudPressure)
        ALLOCATE(Profile%Cld%CloudFraction(L2Prof%nCldPix))
        ALLOCATE(Profile%Cld%CloudPressure(L2Prof%nCldPix))
        IF(Profile%Cld%nSubPix .GT. 0) THEN
          Profile%Cld%CloudFraction = L2Prof%CldPixFraction
          Profile%Cld%CloudPressure = L2Prof%CldPixPressure
          Profile%Cld%TotalCloudFraction = SUM(Profile%Cld%CloudFraction)
          IF(Profile%Cld%TotalCloudFraction .GT. 1.0d0) Profile%Cld%TotalCloudFraction = 1.0d0
        ELSE
          Profile%Cld%TotalCloudFraction = 0.0d0
        ENDIF
        
      ! (*) Scattering Case
      ELSE

        ! Check option exists
        IF(.NOT. L2Prof%ScatteringCloudsPresent) THEN
          CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,  &
                 Message_in='L2 scattering clouds requested but not in profile dataset')
        ENDIF

        ! Reallocate fields
        Profile%Cld%nSubPix = L2Prof%nCldPix
        IF(ALLOCATED(Profile%Cld%CloudFraction)) DEALLOCATE(Profile%Cld%CloudFraction)
        ALLOCATE(Profile%Cld%CloudFraction(L2Prof%nCldPix))

        IF(Profile%Cld%nSubPix .GT. 0) THEN
        
          ! Reallocate arrays
          IF(ALLOCATED(Profile%Cld%LayerOpticalDepth)) DEALLOCATE(Profile%Cld%LayerOpticalDepth)
          IF(ALLOCATED(Profile%Cld%ColumnOpticalDepth)) DEALLOCATE(Profile%Cld%ColumnOpticalDepth)
          IF(ALLOCATED(Profile%Cld%ColumnOptDepthDeriv)) DEALLOCATE(Profile%Cld%ColumnOptDepthDeriv)
          ALLOCATE(Profile%Cld%LayerOpticalDepth(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%ColumnOpticalDepth(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%ColumnOptDepthDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
     
          ! Set arrays
          Profile%Cld%LayerOpticalDepth(:,:,:) = L2Prof%CldLayerOpticalDepth(:,1:Profile%Cld%nSubPix,:)
          Profile%Cld%ColumnOpticalDepth(:,:)  = L2Prof%CldColumnOpticalDepth(1:Profile%Cld%nSubPix,:)

          ! Set Column Optical Depth derivative (uniform scaling)
          Profile%Cld%ColumnOptDepthDeriv(:,:,:) = 0.0d0
          DO c=1,Profile%Cld%nSubPix
          DO s=1,Profile%Cld%nSpecies
            IF(Profile%Cld%ColumnOpticalDepth(c,s) .GT. TINY(0.0d0) ) THEN
              Profile%Cld%ColumnOptDepthDeriv(:,c,s) = Profile%Cld%LayerOpticalDepth(:,c,s) &
                                                     / Profile%Cld%ColumnOpticalDepth(c,s)
            ENDIF
          ENDDO
          ENDDO

        ENDIF

      ENDIF

    ! ----------------------------------------------------------------------------
    ! (3) Option Where we are taking clouds from the climatology (scattering only)
    ! ----------------------------------------------------------------------------
    ELSEIF(Profile%Cld%UseFileCloud) THEN 

      CALL LoadCloudVariables(GridWt,Profile,Error)

    ! ------------------------------------------------------------------
    ! (4) Option Where we are taking clouds from the input file
    ! ------------------------------------------------------------------
    ELSEIF(.NOT. Profile%Cld%UseFileCloud .AND. & 
           .NOT. Profile%Cld%DoLambertianCloud  ) THEN ! Using input file clouds

      ! Input aerosol profile
      aer_zedge(0:Profile%lmx) = Profile%Met%AltitudeEdge

      ! Zero Cloud optical properties
      Profile%Cld%LayerOpticalDepth(:,:,:) = 0.0d0
      Profile%Cld%ColumnOptDepthDeriv(:,:,:) = 0.0d0
      Profile%Cld%AltPeakDeriv(:,:,:) = 0.0d0
      Profile%Cld%AltSigmaDeriv(:,:,:) = 0.0d0
      Profile%Cld%AltExpDeriv(:,:,:) = 0.0d0

      ! Compute Parameterized profiles based on input 
      DO c=1,Profile%Cld%nSubPix
      DO n=1,Profile%Cld%nSpecies
        IF(Profile%Cld%ColumnOpticalDepth(c,n) .GT. TINY(0.0d0)) THEN
          CALL ComputeParameterizedProfile(Profile%Cld%TypeIndex(c),                      &
                Profile%lmx,aer_zedge,Profile%Cld%AltMin(c,n)+aer_zedge(Profile%lmx),     &
                Profile%Cld%AltMax(c,n)+aer_zedge(Profile%lmx),                           &
                Profile%Cld%ColumnOpticalDepth(c,n),                                      &
                Profile%Cld%AltPeak(c,n)+aer_zedge(Profile%lmx),Profile%Cld%AltSigma(c,n),&
                Profile%Cld%AltExp(c,n),Profile%Cld%LayerOpticalDepth(:,c,n),             &
                Profile%Cld%ColumnOptDepthDeriv(:,c,n),Profile%Cld%AltPeakDeriv(:,c,n),   &
                Profile%Cld%AltSigmaDeriv(:,c,n),Profile%Cld%AltExpDeriv(:,c,n), Error    )
        ENDIF
      ENDDO
      ENDDO

    ENDIF
    
    IF(Profile%Cld%DoLambertianCloud) THEN
      
      CALL ComputeLambertianCloudLevel(Profile,Error)

      ! print*,'Lambertian Cloud Levels::::::::'
      ! DO n=1,Profile%Cld%nSubPix
      !   print*,n,Profile%Cld%LambertianCldLevel(n),&
      !            Profile%Cld%LambertianCldLayerFrac(n)
      ! ENDDO

    ENDIF
    
    ! ------------------------------------------------------------------------
    ! Footprint variables (can be replaced by auxiliary information from L1-2)
    ! ------------------------------------------------------------------------

    ! Surface Winds
    IF(Level1%Opt%AuxSurfWind%OverwriteProfileClim) THEN
      Profile%Surface%WindSpeed = Level1%Auxiliary%WindSpeed
      Profile%Surface%WindDirection = Level1%Auxiliary%WindDirection
    ELSE
      CALL LoadSurfaceVar(Profile%ncid, 'wspd',GridWt,     &
                          Profile%Surface%WindSpeed , Error)
      CALL LoadSurfaceVar(Profile%ncid, 'wdir',GridWt,        &
                          Profile%Surface%WindDirection, Error)
    ENDIF

    ! Chlorophyll
    IF(Level1%Opt%AuxChlorophyll%OverwriteProfileClim) THEN
      Profile%Surface%Chlorophyll = Level1%Auxiliary%Chlorophyll
    ELSE
      CALL LoadSurfaceVar(Profile%ncid, 'chphyl',GridWt,      &
                          Profile%Surface%Chlorophyll  , Error)
    ENDIF

    ! Ocean Salinity
    IF(Level1%Opt%AuxSalinity%OverwriteProfileClim) THEN
      Profile%Surface%OceanSalinity = Level1%Auxiliary%OceanSalinity
    ELSE
      CALL LoadSurfaceVar(Profile%ncid, 'ocsal',GridWt,        &
                          Profile%Surface%OceanSalinity , Error)
    ENDIF

    ! Snow Fields
    IF(Level1%Opt%AuxSnow%OverwriteProfileClim) THEN
      Profile%Surface%SnowDepth = Level1%Auxiliary%SnowDepth
      Profile%Surface%SnowFraction = Level1%Auxiliary%SnowFraction
      Profile%Surface%SnowAge = Level1%Auxiliary%SnowAge
    ELSE
      CALL LoadSurfaceVar(Profile%ncid, 'snwdpth',GridWt,  &
                          Profile%Surface%SnowDepth , Error)
      CALL LoadSurfaceVar(Profile%ncid, 'snwfrc',GridWt,     &
                          Profile%Surface%SnowFraction, Error)
      CALL LoadSurfaceVar(Profile%ncid, 'snwage',GridWt,&
                          Profile%Surface%SnowAge, Error)
    ENDIF
    
    ! SIF @ 734 nm
    IF(ncdf_var_exists(Profile%ncid,'sif_734nm')) THEN
      CALL LoadSurfaceVar(Profile%ncid, 'sif_734nm',GridWt,&
                          Profile%Surface%SIF_734nm , Error)
    ELSE
      Profile%Surface%SIF_734nm = 0.0d0
    ENDIF
    
    ! land cover fraction
    Profile%Surface%LandCoverFraction(:)  = 0.0
    Profile%Surface%LandCoverFraction(17) = 1.0 ! Barren land
    
    ! -------------------------------------------------
    ! Load Uncertainties
    ! -------------------------------------------------
    IF(TRIM(ADJUSTL(Profile%CalculationMode)) .EQ. 'INVERSE') THEN
      
      ! Temperature
      CALL LoadVarUncert(Profile%ncid, 'Tmid',&
                         GridWt,Profile%Met%TempUncertainty, Profile, Error)
      
      ! Temperature shift
      CALL LoadVarUncert(Profile%ncid, 'Tshift',&
                         GridWt,Profile%Met%TempShiftUncertainty, Profile, Error)
      
      ! Surface Pressure
      CALL LoadVarUncert(Profile%ncid, 'psurf',&
                         GridWt,Profile%Met%SurfPresUncertainty, Profile, Error)

      ! Gas Profile (FIX-CCM usurp with auxiliary profile uncertaainties?)
      DO n=1,Profile%Gas%nSpecies

        IF(Profile%Gas%UseL2Gas(n)) THEN
          
          ! Get index
          gidx = Profile%Gas%L2Idx(n)
          
          IF(L2Prof%GasUncertainty(gidx)%CovarParType .GT. 0) THEN

            ! Set Profile uncertainty
            Profile%Gas%Uncertainty(n) =  L2Prof%GasUncertainty(gidx)

            ! The subcovariance matrix for the z-correlated parameterization
            ! must be computed here 
            IF(Profile%Gas%Uncertainty(n)%CovarParType .EQ. 2) THEN
              CALL ComputeZCorrCovar(Profile%lmx,                              &
                                     Profile%Gas%Uncertainty(n)%ProfilePar,    &
                                     Profile%Met%AltitudeMid,                  &
                                     Profile%Gas%Uncertainty(n)%ScalarPar(1),  &
                                     Profile%Gas%Uncertainty(n)%SubCovMatrix   )
            ENDIF
            
          ELSE

            ! Take uncertainty from climatology
            CALL LoadVarUncert(Profile%ncid, TRIM(ADJUSTL(Profile%Gas%Name(n))),&
                               GridWt,Profile%Gas%Uncertainty(n), Profile,Error )

          ENDIF

        ELSE

          ! Take uncertainty from climatology
          CALL LoadVarUncert(Profile%ncid, TRIM(ADJUSTL(Profile%Gas%Name(n))),&
                             GridWt,Profile%Gas%Uncertainty(n), Profile,Error )

        ENDIF
        
      ENDDO
      

      ! Aerosol Profile
      DO n=1,Profile%Aer%nSpecies
        CALL LoadVarUncert(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))),&
                           GridWt,Profile%Aer%Uncertainty(n), Profile,Error )
      ENDDO
      
    ENDIF

    ! Special overwrite case - If we are importing a terrain height from
    ! the L1 product then we need to adjust the pressure coordinates
    IF(Level1%Opt%AuxSurfAlt%OverwriteProfileClim) THEN
      
      CALL TerrainHeightAdjust(Level1%Geolocation%SurfaceAltitude,Profile,Error)
        
    ENDIF

    ! Archive the proxy species column mixing ratio
    g = Profile%Gas%ProxyNormIdx
    Profile%Gas%AprioriProxyNormMixingRatio = SUM(Profile%Gas%PartialColumn(:,g)) &
                                            / SUM(Profile%Met%DryAirPartialColumn)
    
  END SUBROUTINE SampleProfile
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: SetAerosolProfile
  ! 
  ! DESCRIPTION: Sets the profile aerosol fields from a supplied set of 
  !              interpolation X-Y grid weights
  
  SUBROUTINE SetAerosolProfile(GridWt,Profile,L1Opt,L2Prof,Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(XYGridWtType),         INTENT(IN)    :: GridWt  
    TYPE(ProfileType),          INTENT(INOUT) :: Profile
    TYPE(L1OptType),            INTENT(IN)    :: L1Opt
    TYPE(AuxProfType),          INTENT(IN)    :: L2Prof
    TYPE(ErrorTYpe),            INTENT(INOUT) :: Error
    
    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: n, idx
    REAL(KIND=8)              :: zmin, zmax, zhgt
    REAL(KIND=8), ALLOCATABLE :: aer_zedge(:)
    CHARACTER(LEN=maxChar)    :: tmpchar, message, action
    LOGICAL                   :: do_derivatives, fail
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SetAerosolProfile'

    ! ============================================================
    ! SetAerosolProfile starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Indexing needed to pass to aerosol routines
    ALLOCATE(aer_zedge(0:Profile%lmx))
    aer_zedge(0:Profile%lmx) = Profile%Met%AltitudeEdge(1:Profile%lmx+1)
    
    do_derivatives = .TRUE.
    
    DO n=1,Profile%Aer%nSpecies
      
      ! Don't overwrite input file 
      IF(Profile%Aer%UseFileAerosol .OR. Profile%Aer%UseL2Aerosol(n)) THEN
        Profile%Aer%AltPeak(n) = 0.0d0
        Profile%Aer%AltSigma(n) = 0.0d0
        Profile%Aer%AltExp(n) = 0.0d0
        Profile%Aer%ColumnOpticalDepth(n) = 0.0d0
      ENDIF
      Profile%Aer%LayerOpticalDepth(:,n) = 0.0d0
      Profile%Aer%ColumnOptDepthDeriv(:,n) = 0.0d0
      Profile%Aer%AltPeakDeriv(:,n) = 0.0d0
      Profile%Aer%AltSigmaDeriv(:,n) = 0.0d0
      Profile%Aer%AltExpDeriv(:,n) = 0.0d0


      ! Get the index in the L2 file
      IF(Profile%Aer%UseL2Aerosol(n)) idx = Profile%Aer%L2Idx(n)
      
      ! (1) Fully specified profile
      IF( Profile%Aer%TypeIndex(n) .EQ. 1 ) THEN

        IF(Profile%Aer%UseL2Aerosol(n)) THEN
          Profile%Aer%LayerOpticalDepth(:,n) = L2Prof%AerLayerOpticalDepth(:,idx)
        ELSE
          IF(Profile%UseL2Met) THEN
            CALL LoadProfileVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))), GridWt,&
                                Profile%Aer%LayerOpticalDepth(:,n), Profile%lmx,Error,   &
                                ClimVertGrid=Profile%ClimVertGrid, IsMixRatio_in=.FALSE. )
          ELSE
            CALL LoadProfileVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))), GridWt,&
                                Profile%Aer%LayerOpticalDepth(:,n), Profile%lmx,Error    )
          ENDIF
        ENDIF
        
        ! Set Column Optical depth (uniform scaling) derivative
        Profile%Aer%ColumnOpticalDepth(n) = SUM(Profile%Aer%LayerOpticalDepth(:,n))
        Profile%Aer%ColumnOptDepthDeriv(:,n) = Profile%Aer%LayerOpticalDepth(:,n) &
                                             / Profile%Aer%ColumnOpticalDepth(n)
        
      ELSE

        ! Overwrite altitude limits if using profile parameters
        IF(Profile%Aer%UseFileAerosol) THEN

          IF(Profile%Aer%UseL2Aerosol(n)) THEN
            Profile%Aer%ColumnOpticalDepth(n) = L2Prof%AerColumnOpticalDepth(idx)
            Profile%Aer%AltMin(n) = L2Prof%AerAltMin(idx)
            Profile%Aer%AltMax(n) = L2Prof%AerAltMax(idx)
          ELSE
            CALL LoadSurfaceVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))),&
                                GridWt, Profile%Aer%ColumnOpticalDepth(n), Error )
            CALL LoadSurfaceVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))) // '_zmin',&
                                GridWt, Profile%Aer%AltMin(n), Error )
            CALL LoadSurfaceVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))) // '_zmax',&
                                GridWt, Profile%Aer%AltMax(n), Error )
          ENDIF

        ENDIF

        IF(Profile%Aer%ColumnOpticalDepth(n) .LT. TINY(0.0d0)) CYCLE

        ! Local copy
        zmin = Profile%Aer%AltMin(n) + aer_zedge(Profile%lmx)
        zmax = Profile%Aer%AltMax(n) + aer_zedge(Profile%lmx)
        
        ! Check for errors
        IF (Profile%Aer%ColumnOpticalDepth(n) < 0.d0) THEN
          STOP "ERROR: Aerosol optical depth must be greater than 0"
        ENDIF
        IF (zmin >= zmax ) THEN
          STOP "ERROR: Aerosol bottom must be lower than aerosol top"
        ENDIF
        
        ! Bound altitudes
        IF(zmin .LT. aer_zedge(Profile%lmx)) zmin = aer_zedge(Profile%lmx)
        IF(zmax .GT. aer_zedge(0          )) zmax = aer_zedge(0          )

        ! (2) GDF
        IF( Profile%Aer%TypeIndex(n) .EQ. 2 ) THEN

          ! Load profile GDF Parameters
          IF( Profile%Aer%UseFileAerosol ) THEN
            IF(Profile%Aer%UseL2Aerosol(n)) THEN
              Profile%Aer%AltPeak(n)  = L2Prof%AerAltPeak(idx)
              Profile%Aer%AltSigma(n) = L2Prof%AerAltSigma(idx)
            ELSE
              CALL LoadSurfaceVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))) // '_pkh',&
                                  GridWt, Profile%Aer%AltPeak(n), Error )
              CALL LoadSurfaceVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))) // '_pkw',&
                                  GridWt, Profile%Aer%AltSigma(n), Error )
            ENDIF
          ENDIF

          ! Adjust Peak sigma based on local height
          zhgt = Profile%Aer%AltPeak(n) + aer_zedge(Profile%lmx)

          ! GDF specific error checks (shouldn't be relevant now referenced to surface ht)
          ! IF (Profile%Aer%AltPeak(n) .GT. zmax) THEN
          !   STOP "ERROR: Aerosol peak height should be <= upperlimit"
          ! ENDIF
          ! IF (Profile%Aer%AltPeak(n) .LT. zmin) THEN
          !   STOP "ERROR: Aerosol peak height should be >= lowerlimit"
          ! ENDIF
          
          ! Compute profile and derivatives
          CALL profiles_gdfone(Profile%lmx, Profile%lmx, aer_zedge, do_derivatives,&
                              zmax,  zhgt,zmin, Profile%Aer%AltSigma(n),          &
                              Profile%Aer%ColumnOpticalDepth(n),                  &
                              Profile%Aer%LayerOpticalDepth(:,n),                 &
                              Profile%Aer%ColumnOptDepthDeriv(:,n),               &
                              Profile%Aer%AltPeakDeriv(:,n),                      &
                              Profile%Aer%AltSigmaDeriv(:,n),                     &
                              fail, message, action                               ) 

        ! (3) EXP
        ELSEIF( Profile%Aer%TypeIndex(n) .EQ. 3 ) THEN

          ! Load profile EXP Parameters
          IF( Profile%Aer%UseFileAerosol ) THEN
            IF(Profile%Aer%UseL2Aerosol(n)) THEN
              Profile%Aer%AltExp(n) = L2Prof%AerAltExp(idx)
            ELSE
              CALL LoadSurfaceVar(Profile%ncid, TRIM(ADJUSTL(Profile%Aer%Name(n))) // '_rxh',&
                                  GridWt, Profile%Aer%AltExp(n), Error )
            ENDIF
          ENDIF
          
          CALL profiles_expone(Profile%lmx, Profile%lmx, aer_zedge, do_derivatives,      &
                              zmax, zmin,                                               &
                              Profile%Aer%AltExp(n), Profile%Aer%ColumnOpticalDepth(n), &
                              Profile%Aer%LayerOpticalDepth(:,n),                       &
                              Profile%Aer%ColumnOptDepthDeriv(:,n),                     &
                              Profile%Aer%AltExpDeriv(:,n),fail, message, action        )

        ! (4) Boxcar
        ELSEIF( Profile%Aer%TypeIndex(n) .EQ. 4 ) THEN
          
          ! Compute profile
          CALL profiles_uniform(Profile%lmx, Profile%lmx, aer_zedge, do_derivatives,&
                                zmax, zmin,                                         &
                                Profile%Aer%ColumnOpticalDepth(n),                  &
                                Profile%Aer%LayerOpticalDepth(:,n),                 &
                                Profile%Aer%ColumnOptDepthDeriv(:,n),               &
                                fail, message, action                               )
          
        ELSE
          
          STOP 'Unrecognized aerosol profile type'
          
        ENDIF
      ENDIF

    ENDDO
    
  END SUBROUTINE SetAerosolProfile
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: ComputeParameterizedProfile
  ! 
  ! DESCRIPTION: Computes a parameterized aerosol/cloud profile based 
  !              on a supplied profile type index
  
  SUBROUTINE ComputeParameterizedProfile(TypeIndex,lmx,aer_zedge,zmin_in,  &
                                         zmax_in,cod,pkh,pkw,rxh,layer_aod,&
                                         dcod,dpkh,dpkw,drxh, Error        )

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER(KIND=2), INTENT(IN)  :: TypeIndex
    INTEGER,         INTENT(IN)  :: lmx
    REAL(KIND=8),    INTENT(IN)  :: aer_zedge(0:lmx)
    REAL(KIND=8),    INTENT(IN)  :: zmin_in
    REAL(KIND=8),    INTENT(IN)  :: zmax_in
    REAL(KIND=8),    INTENT(IN)  :: cod
    REAL(KIND=8),    INTENT(IN)  :: pkh
    REAL(KIND=8),    INTENT(IN)  :: pkw
    REAL(KIND=8),    INTENT(IN)  :: rxh
    REAL(KIND=8),    INTENT(OUT) :: layer_aod(lmx)
    REAL(KIND=8),    INTENT(OUT) :: dcod(lmx)
    REAL(KIND=8),    INTENT(OUT) :: dpkh(lmx)
    REAL(KIND=8),    INTENT(OUT) :: dpkw(lmx)
    REAL(KIND=8),    INTENT(OUT) :: drxh(lmx)
    TYPE(ErrorType), INTENT(INOUT) :: Error

    ! ---------------
    ! Local Variables
    ! ---------------
    REAL(KIND=8) :: zmin,zmax
    LOGICAL :: do_derivatives, fail
    CHARACTER(LEN=maxChar)    :: tmpchar, message, action

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ComputeParameterizedProfile'

    ! ============================================================
    ! ComputeParameterizedProfile starts here
    ! ============================================================

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    do_derivatives = .TRUE.

    ! Zero outputs
    layer_aod(:) = 0.0d0 ; dcod(:) = 0.0d0 ; dpkh(:) = 0.0d0 
    dpkw(:) = 0.0d0 ; drxh(:) = 0.0d0

    ! Copy zmin/max
    zmin = zmin_in ; zmax = zmax_in

    ! If using parameterized profile check height is within altitude bounds

    ! Check for errors
    IF (cod < 0.d0) THEN
      STOP "ERROR: Aerosol optical depth must be greater than 0"
    ENDIF
    IF (zmin >= zmax ) THEN
      STOP "ERROR: Aerosol bottom must be lower than aerosol top"
    ENDIF
        
    ! Bound altitudes
    IF(zmin .LT. aer_zedge(lmx)) zmin = aer_zedge(lmx)
    IF(zmax .GT. aer_zedge(0  )) zmax = aer_zedge(0  )
    
    ! (2) GDF
    IF( TypeIndex .EQ. 2 ) THEN

      ! GDF specific error checks
      IF (pkh .GT. zmax) THEN
        STOP "ERROR: Aerosol peak height should be <= upperlimit"
      ENDIF
      IF (pkh .LT. zmin) THEN
        STOP "ERROR: Aerosol peak height should be >= lowerlimit"
      ENDIF
        
      ! Compute profile and derivatives
      CALL profiles_gdfone(lmx, lmx, aer_zedge, do_derivatives,&
                           zmax, pkh, zmin, pkw, cod,          &
                           layer_aod,dcod, dpkh, dpkw,         &
                           fail, message, action               ) 

      ! (3) EXP
      ELSEIF( TypeIndex .EQ. 3 ) THEN

        CALL profiles_expone(lmx, lmx, aer_zedge, do_derivatives,&
                             zmax, zmin,rxh, cod, layer_aod,     &
                             dcod, drxh,fail, message, action    )

      ! (4) Boxcar
      ELSEIF( TypeIndex .EQ. 4 ) THEN
        
        ! Compute profile
        CALL profiles_uniform(lmx, lmx, aer_zedge, do_derivatives,&
                              zmax, zmin, cod,layer_aod,  dcod,   &
                              fail, message, action               )
        
      ELSE
        
        STOP 'Unrecognized aerosol profile type'
        
      ENDIF

  END SUBROUTINE ComputeParameterizedProfile
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: TerrainHeightAdjust
  ! 
  ! DESCRIPTION: Adjusts the climatological profile based on an 
  !              updated surface terrain height
  
  SUBROUTINE TerrainHeightAdjust(NewSurfAltitude,Profile,Error)
      
    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL(KIND=8),           INTENT(IN)    :: NewSurfAltitude
    TYPE(ProfileType),      INTENT(INOUT) :: Profile
    TYPE(ErrorType),        INTENT(INOUT) :: Error
    
    
    ! ---------------
    ! Local Variables
    ! ---------------
    TYPE(ProfileType) :: Profile0
    INTEGER           :: l,ll,g,n
    REAL(KIND=8)      :: sf
    REAL(KIND=8)      :: NewSurfacePres, NewSurfTemp
    REAL(KIND=8)      :: VertWeights(Profile%lmx,Profile%lmx)
    REAL(KIND=8)      :: z0, T0, p0, g0, mr0, mrtot, Rspec, dp
    REAL(KIND=8)      :: p_eff, T_eff
    REAL(KIND=8)      :: p_int(Profile%lmx+2),T_int(Profile%lmx+2)
    INTEGER           :: n_int, errstat
    REAL(KIND=8)      :: tmp_prof(Profile%lmx), tot_mr(Profile%lmx)
    REAL(KIND=8)      :: aer_zedge(0:Profile%lmx)

    ! Assumptions for terrain height adjustment
    REAL(kind=8), PARAMETER :: Gamma = 6.5d0 ! Environmental lapse rate (K/km)
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'TerrainHeightAdjust'

    ! ============================================================
    ! TerrainHeightAdjust starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Save initial profile
    Profile0 = Profile
    
    ! Shorthand for original terrain height values
    z0 = Profile0%Met%AltitudeEdge(Profile0%lmx+1)
    T0 = Profile0%Met%TemperatureEdge(Profile0%lmx+1)
    p0 = Profile0%Met%PressureEdge(Profile0%lmx+1)
    g0 = Profile0%Met%Gravity(Profile0%lmx)
    mr0 = 0.0d0 ; mrtot = 0.0
    DO g=1,Profile0%Gas%nSpecies
      mr0 = mr0 + Profile0%Gas%MolecularWeight(g)*Profile0%Gas%MixingRatio(Profile0%lmx,g)
      mrtot = mrtot + Profile0%Gas%MixingRatio(Profile0%lmx,g)
    ENDDO
    mr0 = mr0 / mrtot
    Rspec = Constants_R / mr0
    
    ! Compute the new Surface pressure and temperature
    p_eff = p0 *(T0/(T0+Gamma*(z0-NewSurfAltitude)))**&
                (-1.0*g0/(Rspec*Gamma*1e-3))
    T_eff = T0 + Gamma*(z0-NewSurfAltitude)

    ! Compute the new pressure edges
    DO l=1,Profile%lmx+1
      Profile%Met%PressureEdge(l) = Profile%Met%AP(l) &
                                  + Profile%Met%BP(l)*p_eff
    ENDDO

    ! Compute new pressure midpoints
    DO l=1,Profile%lmx
      Profile%Met%PressureMid(l) = 0.5d0*(Profile%Met%PressureEdge(l)  &
                                         +Profile%Met%PressureEdge(l+1))
    ENDDO

    ! ------------------------
    ! Compute new Temperatures
    ! ------------------------
    IF(NewSurfAltitude .LT. Profile0%Met%AltitudeEdge(Profile0%lmx+1)) THEN
      n_int = Profile%lmx+2
      p_int(1) = LOG(p_eff) ; T_int(1) = T_eff
      p_int(2:n_int) = LOG(Profile0%Met%PressureEdge)
      T_int(2:n_int) = Profile0%Met%TemperatureEdge
    ELSE
      n_int = Profile%lmx+1
      p_int(1:n_int) = LOG(Profile0%Met%PressureEdge) ! T deriv is linear with height
      T_int(1:n_int) = Profile0%Met%TemperatureEdge
    ENDIF

    ! Edges
    CALL BSPLINE_EdgeFill(p_int,T_int,n_int,LOG(Profile%Met%PressureEdge),&
                          Profile%Met%TemperatureEdge,Profile%lmx+1,errstat)

    ! Midpoints
    CALL BSPLINE_EdgeFill(p_int,T_int,n_int,LOG(Profile%Met%PressureMid),&
                          Profile%Met%TemperatureMid,Profile%lmx+1,errstat)


    ! -------------------
    ! Adjust Gas Profiles 
    ! -------------------

    ! Keep track of total mixing ratio
    tot_mr(:) = 0.0d0

    DO g=1,Profile%Gas%nSpecies

      SELECT CASE (Profile%Gas%SigmaAdjustType(g))
        
        CASE(1) ! BACKGROUND (Preserve vertical mixing ratios)
          
          ! Interpolate Concentrations to pressure Midpoints
          CALL BSPLINE_EdgeFill(Profile0%Met%PressureMid,     &
                                Profile0%Gas%MixingRatio(:,g),&
                                Profile%lmx,                  &
                                Profile%Met%PressureMid,      &
                                Profile%Gas%MixingRatio(:,g), &
                                Profile%lmx, errstat          )
          
        CASE(2) ! SIGMA (Preserve layer column densities)
          
          ! Preserve Layer Densities (weight by pressure)
          DO l=1,Profile%lmx
            Profile%Gas%MixingRatio(l,g) = Profile0%Gas%MixingRatio(l,g)*&
                       (Profile%Met%PressureEdge(l+1)-Profile%Met%PressureEdge(l))/&
                      (Profile0%Met%PressureEdge(l+1)-Profile0%Met%PressureEdge(l))
          ENDDO
          
        CASE(3) ! HYBRID (SIGMA Adjustment on residual profile)
          
          ! Interpolate Background Concentrations to pressure Midpoints
          CALL BSPLINE_EdgeFill(Profile0%Met%PressureMid,               &
                                Profile0%Gas%BackgroundMixingRatio(:,g),&
                                Profile%lmx,                            &
                                Profile%Met%PressureMid,                &
                                Profile%Gas%BackgroundMixingRatio(:,g), &
                                Profile%lmx, errstat                    )

          ! Compute Residual Profile
          tmp_prof(:) = Profile0%Gas%MixingRatio(:,g) &
                      - Profile0%Gas%BackgroundMixingRatio(:,g)

          ! Sigma Adjust original profile and add background
          DO l=1,Profile%lmx
            Profile%Gas%MixingRatio(l,g) = tmp_prof(l)*                         &
                  (Profile%Met%PressureEdge(l+1)-Profile%Met%PressureEdge(l))/  &
                  (Profile0%Met%PressureEdge(l+1)-Profile0%Met%PressureEdge(l)) &
                 + Profile%Gas%BackgroundMixingRatio(l,g) 
          ENDDO

        CASE DEFAULT
          
          CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,&
                                Message_in='Urecognized Sigma adjust type'           )
          
      END SELECT

      ! Update total mixing ratio
      tot_mr(:) = tot_mr(:) + Profile%Gas%MixingRatio(:,g)

    ENDDO

    ! Ensure normalization of total mixing ratio
    DO g=1,Profile%Gas%nSpecies
      Profile%Gas%MixingRatio(:,g) = Profile%Gas%MixingRatio(:,g) / tot_mr(:)
    ENDDO

    ! Recompute the derived quantities
    CALL ComputeProfileDerivedQuantities(Profile, Error)
    
    ! ----------------
    ! Aerosol Profiles
    ! ----------------
    aer_zedge(0:Profile%lmx) = Profile%Met%AltitudeEdge(1:Profile%lmx+1)

    ! For now we assume all profiles undergo SIGMA Adjustment

    DO g=1,Profile%Aer%nSpecies

      ! Recompute profiles with new height if they are parameterized
      IF(Profile%Aer%TypeIndex(g) .GT. 1) THEN
          
        CALL ComputeParameterizedProfile(Profile%Aer%TypeIndex(g),              &
            Profile%lmx,aer_zedge,Profile%Aer%AltMin(g)+aer_zedge(Profile%lmx),   &
            Profile%Aer%AltMax(g)+aer_zedge(Profile%lmx),                         &
            Profile%Aer%ColumnOpticalDepth(g),                                    &
            Profile%Aer%AltPeak(g)+aer_zedge(Profile%lmx),Profile%Aer%AltSigma(g),&
            Profile%Aer%AltExp(g),Profile%Aer%LayerOpticalDepth(:,g),             &
            Profile%Aer%ColumnOptDepthDeriv(:,g),Profile%Aer%AltPeakDeriv(:,g),   &
            Profile%Aer%AltSigmaDeriv(:,g),Profile%Aer%AltExpDeriv(:,g), Error    )
          
      ENDIF

    ENDDO

    ! --------------
    ! Cloud Profiles
    ! --------------
    
    IF(Profile%Cld%DoLambertianCloud) THEN

      ! Recompute the Cloud heights on the new pressure grid
      CALL ComputeLambertianCloudLevel(Profile,Error)

    ELSE

      DO g=1,Profile%Cld%nSpecies

        
        IF(Profile%Cld%TypeIndex(g) .GT. 1) THEN

          DO n=1,Profile%Cld%nSubPix

            IF(Profile%Cld%ColumnOpticalDepth(n,g) .GT. 0.0d0) THEN
              CALL ComputeParameterizedProfile(Profile%Cld%TypeIndex(g),                  &
                Profile%lmx,aer_zedge,Profile%Cld%AltMin(n,g)+aer_zedge(Profile%lmx),     &
                Profile%Cld%AltMax(n,g)+aer_zedge(Profile%lmx),                           &
                Profile%Cld%ColumnOpticalDepth(n,g),                                      &
                Profile%Cld%AltPeak(n,g)+aer_zedge(Profile%lmx),Profile%Cld%AltSigma(n,g),&
                Profile%Cld%AltExp(n,g),Profile%Cld%LayerOpticalDepth(:,n,g),             &
                Profile%Cld%ColumnOptDepthDeriv(:,n,g),Profile%Cld%AltPeakDeriv(:,n,g),   &
                Profile%Cld%AltSigmaDeriv(:,n,g),Profile%Cld%AltExpDeriv(:,n,g), Error    )
            ENDIF

          ENDDO

        ENDIF

      ENDDO

    ENDIF
    
    ! Adjustments to uncertainties go here
    
    !STOP 'Testing TerrainHeightAdjust'
    
  END SUBROUTINE TerrainHeightAdjust
  
  SUBROUTINE LoadCloudVariables(GridWt,Profile,Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(XYGridWtType),         INTENT(IN)    :: GridWt  
    TYPE(ProfileType),          INTENT(INOUT) :: Profile
    TYPE(ErrorTYpe),            INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER                   :: n, m, i, j, l, t, ct, p
    
    TYPE Fld2DType
      INTEGER                      :: vmin(3,2), nv(3,2)
      TYPE(ArrSetType_3D)          :: fld
      ! INTEGER(KIND=2), ALLOCATABLE :: fld_i2(:,:,:)
      ! REAL(KIND=8),    ALLOCATABLE :: fld_r8(:,:,:)
      INTEGER, ALLOCATABLE         :: idx(:,:), jdx(:), tdx(:)
    ENDTYPE Fld2DType
    TYPE Fld3DType
      INTEGER                      :: vmin(4,2), nv(4,2)
      TYPE(ArrSetType_4D)          :: fld
      INTEGER, ALLOCATABLE         :: idx(:,:),jdx(:),ldx(:),tdx(:)
    ENDTYPE Fld3DType

    ! # Cloud pixels
    TYPE(Fld2DType)  :: ncpix
    TYPE(Fld2DType)  :: cf
    TYPE(Fld2DType)  :: pcld
    TYPE(Fld2DType)  :: cod
    TYPE(Fld2DType)  :: zmin
    TYPE(Fld2DType)  :: zmax
    TYPE(Fld2DType)  :: zpk
    TYPE(Fld2DType)  :: zwd
    TYPE(Fld2DType)  :: zexp
    TYPE(Fld3DType)  :: cldprof

    REAL(KIND=8) :: cod_in,zmin_in,zmax_in,zpk_in,zwd_in,zexp_in


    
    INTEGER(KIND=2)           :: this_ncpix
    REAL(KIND=8)              :: this_cf
    REAL(KIND=8), ALLOCATABLE :: aer_zedge(:)
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadCloudVariables'

    ! ============================================================
    ! LoadCloudVariables starts here
    ! ============================================================

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Indexing needed to pass to aerosol routines
    ALLOCATE(aer_zedge(0:Profile%lmx))
    aer_zedge(0:Profile%lmx) = Profile%Met%AltitudeEdge(1:Profile%lmx+1)

    IF(Profile%Cld%UseFileCloud) THEN

      CALL ReadSurfaceVarTmpArray(Profile%ncid,'ncpix',GridWt,ncpix%vmin,ncpix%nv,&
                                  ncpix%idx,ncpix%jdx,ncpix%tdx,ncpix%fld,Error,ReadI2_in=.TRUE.)

      ! CALL ReadSurfaceVarTmpArray_i2(Profile%ncid,'ncpix',GridWt,ncpix%vmin,ncpix%nv,&
      !                               ncpix%idx,ncpix%jdx,ncpix%tdx,ncpix%fld_i2,Error)
      
      CALL ReadSurfaceVarTmpArray(Profile%ncid,'cftot',GridWt,cf%vmin,cf%nv,&
                                  cf%idx,cf%jdx,cf%tdx,cf%fld,Error)
      
      IF(Profile%Cld%DoLambertianCloud) THEN
        
        CALL ReadSurfaceVarTmpArray(Profile%ncid,'pcld',GridWt,pcld%vmin,pcld%nv,&
                                  pcld%idx,pcld%jdx,pcld%tdx,pcld%fld,Error  )
        
      ENDIF
      
      ! Count the cloudy pixels
      ! Now Compute the averaged profile
      Profile%Cld%TotalCloudFraction = 0.0d0
      Profile%Cld%nSubPix = 0
      
      DO t=1,GridWt%nt
      DO j=1,Gridwt%nj
      DO i=1,Gridwt%ni
        IF(GridWt%Weight(i,j,t) .GT. TINY(0.0d0)) THEN
          this_cf    = cf%fld%data(1)%arr(cf%idx(i,1),cf%jdx(j),cf%tdx(t))
          this_ncpix = ncpix%fld%data(1)%arr_i2(ncpix%idx(i,1),ncpix%jdx(j),ncpix%tdx(t)) ! should be 0 or 1
          IF(this_ncpix .GT. 0) THEN
            Profile%Cld%TotalCloudFraction = Profile%Cld%TotalCloudFraction &
                                          + GridWt%Weight(i,j,t)*this_cf
            Profile%Cld%nSubPix = Profile%Cld%nSubPix + 1
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ! Include second subpixel in average if needed
      IF(GridWt%PeriodicOverlap) THEN
        DO t=1,GridWt%nt
        DO j=1,Gridwt%nj
        DO i=1,Gridwt%ni2
          IF(GridWt%Weight2(i,j,t) .GT. TINY(0.0d0)) THEN
            this_cf    = cf%fld%data(2)%arr(cf%idx(i,2),cf%jdx(j),cf%tdx(t))
            this_ncpix = ncpix%fld%data(2)%arr_i2(ncpix%idx(i,2),ncpix%jdx(j),ncpix%tdx(t)) ! should be 0 or 1
            IF(this_ncpix .GT. 0) THEN
              Profile%Cld%TotalCloudFraction = Profile%Cld%TotalCloudFraction &
                                            + GridWt%Weight2(i,j,t)*this_cf
              Profile%Cld%nSubPix = Profile%Cld%nSubPix + 1
            ENDIF
          ENDIF
        ENDDO
        ENDDO
        ENDDO
      ENDIF
      
      ! Allocate sub pixel arrays
      ! -------
      IF(Profile%Cld%nSubPix .GT. 0) THEN
        
        IF(ALLOCATED(Profile%Cld%CloudFraction)) DEALLOCATE(Profile%Cld%CloudFraction)
        ALLOCATE(Profile%Cld%CloudFraction(Profile%Cld%nSubPix))
        
        IF(Profile%Cld%DoLambertianCloud) THEN
          
          IF(ALLOCATED(Profile%Cld%CloudPressure)) DEALLOCATE(Profile%Cld%CloudPressure)
          ALLOCATE(Profile%Cld%CloudPressure(Profile%Cld%nSubPix))
          IF(ALLOCATED(Profile%Cld%LambertianCldLevel)) DEALLOCATE(Profile%Cld%LambertianCldLevel)
          ALLOCATE(Profile%Cld%LambertianCldLevel(Profile%Cld%nSubPix))
          IF(ALLOCATED(Profile%Cld%LambertianCldLayerFrac)) DEALLOCATE(Profile%Cld%LambertianCldLayerFrac)
          ALLOCATE(Profile%Cld%LambertianCldLayerFrac(Profile%Cld%nSubPix))

        ELSE
          
          ! Reallocate arrays
          IF(ALLOCATED(Profile%Cld%LayerOpticalDepth)) DEALLOCATE(Profile%Cld%LayerOpticalDepth)
          IF(ALLOCATED(Profile%Cld%ColumnOpticalDepth)) DEALLOCATE(Profile%Cld%ColumnOpticalDepth)
          IF(ALLOCATED(Profile%Cld%ColumnOptDepthDeriv)) DEALLOCATE(Profile%Cld%ColumnOptDepthDeriv)
          IF(ALLOCATED(Profile%Cld%AltMin)) DEALLOCATE(Profile%Cld%AltMin)
          IF(ALLOCATED(Profile%Cld%AltMax)) DEALLOCATE(Profile%Cld%AltMax)
          IF(ALLOCATED(Profile%Cld%AltPeak)) DEALLOCATE(Profile%Cld%AltPeak)
          IF(ALLOCATED(Profile%Cld%AltSigma)) DEALLOCATE(Profile%Cld%AltSigma)
          IF(ALLOCATED(Profile%Cld%AltExp)) DEALLOCATE(Profile%Cld%AltExp)
          IF(ALLOCATED(Profile%Cld%AltPeakDeriv)) DEALLOCATE(Profile%Cld%AltPeakDeriv)
          IF(ALLOCATED(Profile%Cld%AltSigmaDeriv)) DEALLOCATE(Profile%Cld%AltSigmaDeriv)
          IF(ALLOCATED(Profile%Cld%AltExpDeriv)) DEALLOCATE(Profile%Cld%AltExpDeriv)
          ALLOCATE(Profile%Cld%LayerOpticalDepth(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%ColumnOpticalDepth(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%ColumnOptDepthDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltMin(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltMax(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltPeak(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltSigma(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltExp(Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltPeakDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltSigmaDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          ALLOCATE(Profile%Cld%AltExpDeriv(Profile%lmx,Profile%Cld%nSubPix,Profile%Cld%nSpecies))
          
          ! Zero arrays
          Profile%Cld%LayerOpticalDepth(:,:,:) = 0.0d0
          Profile%Cld%ColumnOpticalDepth(:,:) = 0.0d0
          Profile%Cld%ColumnOptDepthDeriv(:,:,:) = 0.0d0
          Profile%Cld%AltMin(:,:) = 0.0d0
          Profile%Cld%AltMax(:,:) = 0.0d0
          Profile%Cld%AltPeak(:,:) = 0.0d0
          Profile%Cld%AltSigma(:,:) = 0.0d0
          Profile%Cld%AltExp(:,:) = 0.0d0
          Profile%Cld%AltPeakDeriv(:,:,:) = 0.0d0
          Profile%Cld%AltSigmaDeriv(:,:,:) = 0.0d0
          Profile%Cld%AltExpDeriv(:,:,:) = 0.0d0

        ENDIF
        
        ! Compute Subpixel fraction
        ct = 0

        DO t=1,GridWt%nt
        DO j=1,Gridwt%nj
        DO i=1,Gridwt%ni

          IF(GridWt%Weight(i,j,t) .GT. TINY(0.0d0)) THEN
            
            ! Get values
            this_cf    = cf%fld%Data(1)%arr(cf%idx(i,1),cf%jdx(j),cf%tdx(t))
            this_ncpix = ncpix%fld%Data(1)%arr_i2(ncpix%idx(i,1),ncpix%jdx(j),ncpix%tdx(t)) ! should be 0 or 1
            
            IF(this_ncpix .GT. 0) THEN
              
              ! Increment count
              ct = ct + 1
              
              ! Compute cloud fraction
              Profile%Cld%CloudFraction(ct) = GridWt%Weight(i,j,t)*this_cf
              
              ! -----------------------------------
              ! Fields needed for lambertian clouds
              ! -----------------------------------
              IF(Profile%Cld%DoLambertianCloud) THEN
                Profile%Cld%CloudPressure(ct) = GridWt%Weight(i,j,t)* &
                                                pcld%fld%data(1)%arr(pcld%idx(i,1),pcld%jdx(j),pcld%tdx(t))
              ENDIF

            ENDIF
            
          ENDIF
          
        ENDDO
        ENDDO
        ENDDO
        
        IF(GridWt%PeriodicOverlap) THEN

          DO t=1,GridWt%nt
          DO j=1,Gridwt%nj
          DO i=1,Gridwt%ni2

            IF(GridWt%Weight2(i,j,t) .GT. TINY(0.0d0)) THEN
              
              ! Get values
              this_cf    = cf%fld%data(2)%arr(cf%idx(i,2),cf%jdx(j),cf%tdx(t))
              this_ncpix = ncpix%fld%data(2)%arr_i2(ncpix%idx(i,2),ncpix%jdx(j),ncpix%tdx(t)) ! should be 0 or 1
              
              IF(this_ncpix .GT. 0) THEN
                
                ! Increment count
                ct = ct + 1
                
                ! Compute cloud fraction
                Profile%Cld%CloudFraction(ct) = GridWt%Weight2(i,j,t)*this_cf
                
                ! -----------------------------------
                ! Fields needed for lambertian clouds
                ! -----------------------------------
                IF(Profile%Cld%DoLambertianCloud) THEN
                  Profile%Cld%CloudPressure(ct) = GridWt%Weight(i,j,t)* &
                                                  pcld%fld%data(2)%arr(pcld%idx(i,2),pcld%jdx(j),pcld%tdx(t))
                ENDIF

              ENDIF
              
            ENDIF
            
          ENDDO
          ENDDO
          ENDDO

        ENDIF

        ! --------------------------------------------------------------------------
        ! Scattering cloud case
        ! --------------------------------------------------------------------------
        IF(.NOT. Profile%Cld%DoLambertianCloud) THEN
          
          ! Load Data for the profile if from disk
          DO n=1,Profile%Cld%nSpecies
            
            ! Specified Profiles
            ! ------------------
            IF( Profile%Cld%TypeIndex(n) .EQ. 1 ) THEN
              
              CALL ReadProfileVarTmpArray(Profile%ncid,&
                                          TRIM(ADJUSTL(Profile%Cld%Name(n))),&
                                          GridWt,Profile%lmx,cldprof%vmin,   &
                                          cldprof%nv,cldprof%idx,cldprof%jdx,&
                                          cldprof%ldx,cldprof%tdx,           &
                                          cldprof%fld,Error                  )
              
              ! Set Profile
              ct = 0
              DO t=1,GridWt%nt
              DO j=1,Gridwt%nj
              DO i=1,Gridwt%ni

                ! Increment Count
                ct = ct + 1

                ! Set Profile
                DO l=1,Profile%lmx
                  Profile%Cld%LayerOpticalDepth(l,ct,n) = &
                    cldprof%fld%data(1)%arr(cldprof%idx(i,1),cldprof%jdx(j),cldprof%ldx(l),cldprof%tdx(t))
                  Profile%Cld%ColumnOpticalDepth(ct,n) = Profile%Cld%ColumnOpticalDepth(ct,n) &
                                                       + Profile%Cld%LayerOpticalDepth(l,ct,n)
                ENDDO

              ENDDO
              ENDDO
              ENDDO
              
              IF(GridWt%PeriodicOverlap) THEN
                DO t=1,GridWt%nt
                DO j=1,Gridwt%nj
                DO i=1,Gridwt%ni2

                  ! Increment Count
                  ct = ct + 1

                  ! Set Profile
                  DO l=1,Profile%lmx
                    Profile%Cld%LayerOpticalDepth(l,ct,n) = &
                      cldprof%fld%data(2)%arr(cldprof%idx(i,2),cldprof%jdx(j),cldprof%ldx(l),cldprof%tdx(t))
                    Profile%Cld%ColumnOpticalDepth(ct,n) = Profile%Cld%ColumnOpticalDepth(ct,n) &
                                                        + Profile%Cld%LayerOpticalDepth(l,ct,n)
                  ENDDO

                ENDDO
                ENDDO
                ENDDO
              ENDIF

              ! Set Column Optical Depth derivative (uniform scaling)
              DO ct=1,Profile%Cld%nSubPix
                Profile%Cld%ColumnOptDepthDeriv(:,ct,n) = Profile%Cld%LayerOpticalDepth(:,ct,n) &
                                                       / Profile%Cld%ColumnOpticalDepth(ct,n)
              ENDDO
            
            ! Parameterized Profiles
            ! ----------------------
            ELSE

              ! Overwrite altitude limits if using profile parameters
              CALL ReadSurfaceVarTmpArray(Profile%ncid,                      &
                                          TRIM(ADJUSTL(Profile%Cld%Name(n))),&
                                          GridWt,cod%vmin,cod%nv,cod%idx,    &
                                          cod%jdx,cod%tdx,cod%fld,Error      )
              CALL ReadSurfaceVarTmpArray(Profile%ncid,                                &
                                          TRIM(ADJUSTL(Profile%Cld%Name(n)))// '_zmin',&
                                          GridWt,zmin%vmin,zmin%nv,zmin%idx,           &
                                          zmin%jdx,zmin%tdx,zmin%fld,Error             )
              CALL ReadSurfaceVarTmpArray(Profile%ncid,                                &
                                          TRIM(ADJUSTL(Profile%Cld%Name(n)))// '_zmax',&
                                          GridWt,zmax%vmin,zmax%nv,zmax%idx,           &
                                          zmax%jdx,zmax%tdx,zmax%fld,Error             )

              IF(Profile%Cld%TypeIndex(n) .EQ. 2) THEN
                CALL ReadSurfaceVarTmpArray(Profile%ncid,                             &
                                          TRIM(ADJUSTL(Profile%Cld%Name(n)))// '_pkh',&
                                          GridWt,zpk%vmin,zpk%nv,zpk%idx,             &
                                          zpk%jdx,zpk%tdx,zpk%fld,Error               )
                CALL ReadSurfaceVarTmpArray(Profile%ncid,                             &
                                          TRIM(ADJUSTL(Profile%Cld%Name(n)))// '_pkw',&
                                          GridWt,zwd%vmin,zwd%nv,zwd%idx,             &
                                          zwd%jdx,zwd%tdx,zwd%fld,Error               )

              ELSEIF(Profile%Cld%TypeIndex(n) .EQ. 3) THEN
                CALL ReadSurfaceVarTmpArray(Profile%ncid,                             &
                                          TRIM(ADJUSTL(Profile%Cld%Name(n)))// '_rxh',&
                                          GridWt,zexp%vmin,zexp%nv,zexp%idx,          &
                                          zexp%jdx,zexp%tdx,zexp%fld,Error            )
              ELSEIF(Profile%Cld%TypeIndex(n) .GT. 4) THEN
                  STOP 'Unrecognized aerosol profile type'
              ENDIF

              ! Compute the Cloud Profiles
              ct = 0 ; cod_in = 0.0d0 ; zmin_in = 0.0d0 ; zmax_in = 0.0d0 
              zpk_in = 0.0d0 ; zwd_in = 0.0d0 ; zexp_in = 0.0d0

              DO t=1,GridWt%nt
              DO j=1,Gridwt%nj
              DO i=1,Gridwt%ni

                ! Set inputs to cloud profile parameterization
                cod_in  = cod%fld%data(1)%arr(cod%idx(i,1),cod%jdx(j),cod%tdx(t))
                zmin_in = zmin%fld%data(1)%arr(zmin%idx(i,1),zmin%jdx(j),zmin%tdx(t))
                zmax_in = zmax%fld%data(1)%arr(zmax%idx(i,1),zmax%jdx(j),zmax%tdx(t))
                IF(Profile%Cld%TypeIndex(n) .EQ. 2) THEN
                  zpk_in = zpk%fld%data(1)%arr(zpk%idx(i,1),zpk%jdx(j),zpk%tdx(t))
                  zwd_in = zwd%fld%data(1)%arr(zwd%idx(i,1),zwd%jdx(j),zwd%tdx(t))
                ELSEIF(Profile%Cld%TypeIndex(n) .EQ. 3) THEN
                  zexp_in = zexp%fld%data(1)%arr(zexp%idx(i,1),zexp%jdx(j),zexp%tdx(t))
                ENDIF

                ! Increment Count
                ct = ct + 1

                ! Archive the Profile Parameters
                Profile%Cld%ColumnOpticalDepth(ct,n) = cod_in
                Profile%Cld%AltMin(ct,n)             = zmin_in
                Profile%Cld%AltMax(ct,n)             = zmax_in
                Profile%Cld%AltPeak(ct,n)            = zpk_in
                Profile%Cld%AltSigma(ct,n)           = zwd_in
                Profile%Cld%AltExp(ct,n)             = zexp_in
                
                ! Compute Profile for pixel
                IF(cod_in .GT. 0.0d0) THEN
                  CALL ComputeParameterizedProfile(Profile%Cld%TypeIndex(n),Profile%lmx,   &
                                                  aer_zedge,zmin_in+aer_zedge(Profile%lmx),&
                                                  zmax_in+aer_zedge(Profile%lmx),cod_in,   &
                                                  zpk_in+aer_zedge(Profile%lmx),zwd_in,    &
                                                  zexp_in,                                 &
                                                  Profile%Cld%LayerOpticalDepth(:,ct,n),   &
                                                  Profile%Cld%ColumnOptDepthDeriv(:,ct,n), &
                                                  Profile%Cld%AltPeakDeriv(:,ct,n),        &
                                                  Profile%Cld%AltSigmaDeriv(:,ct,n),       &
                                                  Profile%Cld%AltExpDeriv(:,ct,n), Error   )
                ENDIF

              ENDDO
              ENDDO
              ENDDO

              IF(GridWt%PeriodicOverlap) THEN
                DO t=1,GridWt%nt
                DO j=1,Gridwt%nj
                DO i=1,Gridwt%ni2

                  ! Set inputs to cloud profile parameterization
                  cod_in  = cod%fld%data(2)%arr(cod%idx(i,2),cod%jdx(j),cod%tdx(t))
                  zmin_in = zmin%fld%data(2)%arr(zmin%idx(i,2),zmin%jdx(j),zmin%tdx(t))
                  zmax_in = zmax%fld%data(2)%arr(zmax%idx(i,2),zmax%jdx(j),zmax%tdx(t))
                  IF(Profile%Cld%TypeIndex(n) .EQ. 2) THEN
                    zpk_in = zpk%fld%data(2)%arr(zpk%idx(i,2),zpk%jdx(j),zpk%tdx(t))
                    zwd_in = zwd%fld%data(2)%arr(zwd%idx(i,2),zwd%jdx(j),zwd%tdx(t))
                  ELSEIF(Profile%Cld%TypeIndex(n) .EQ. 3) THEN
                    zexp_in = zexp%fld%data(2)%arr(zexp%idx(i,1),zexp%jdx(j),zexp%tdx(t))
                  ENDIF

                  ! Increment Count
                  ct = ct + 1

                  ! Archive the Profile Parameters
                  Profile%Cld%ColumnOpticalDepth(ct,n) = cod_in
                  Profile%Cld%AltMin(ct,n)             = zmin_in
                  Profile%Cld%AltMax(ct,n)             = zmax_in
                  Profile%Cld%AltPeak(ct,n)            = zpk_in
                  Profile%Cld%AltSigma(ct,n)           = zwd_in
                  Profile%Cld%AltExp(ct,n)             = zexp_in
                  
                  ! Compute Profile for pixel
                  IF(cod_in .GT. 0.0d0) THEN
                    CALL ComputeParameterizedProfile(Profile%Cld%TypeIndex(n),Profile%lmx,   &
                                                    aer_zedge,zmin_in+aer_zedge(Profile%lmx),&
                                                    zmax_in+aer_zedge(Profile%lmx),cod_in,   &
                                                    zpk_in+aer_zedge(Profile%lmx),zwd_in,    &
                                                    zexp_in,                                 &
                                                    Profile%Cld%LayerOpticalDepth(:,ct,n),   &
                                                    Profile%Cld%ColumnOptDepthDeriv(:,ct,n), &
                                                    Profile%Cld%AltPeakDeriv(:,ct,n),        &
                                                    Profile%Cld%AltSigmaDeriv(:,ct,n),       &
                                                    Profile%Cld%AltExpDeriv(:,ct,n), Error   )
                  ENDIF

                ENDDO
                ENDDO
                ENDDO
              ENDIF

            ENDIF

          ENDDO ! nSpecies

        ENDIF ! DoLambertianCloud

      ENDIF ! CloudSubPix > 0
    
    ENDIF ! UseFileCloud

  END SUBROUTINE LoadCloudVariables
  
  SUBROUTINE ComputeLambertianCloudLevel(Profile,Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ProfileType), INTENT(INOUT) :: Profile
    TYPE(ErrorType),   INTENT(INOUT) :: Error

    ! ---------------
    ! Local Variables
    ! ---------------
    REAL(KIND=8) :: y2(Profile%lmx+1), plog(Profile%lmx+1)
    REAL(KIND=8) :: log_pcld(1), zcld(1)
    INTEGER      :: l,c

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ComputeLambertianCloudLevel'

    ! ============================================================
    ! ComputeLambertianCloudLevel starts here
    ! ============================================================

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Initialize arrays
    Profile%Cld%LambertianCldLevel(:)     = 0
    Profile%Cld%LambertianCldLayerFrac(:) = 0.0d0
    
    ! Compute height vs. logP
    plog = LOG10(Profile%Met%PressureEdge)

    ! If top layer is zero, compute height at pressure x1e-2 of the layer below
    ! This is the same assumption as applied to the Altitude
    IF(Profile%Met%PressureEdge(1) .LT. TINY(0.0d0)) THEN
      plog(1) = LOG10(Profile%Met%PressureEdge(2)*1e-2)
    ENDIF

    ! Compute splines
    CALL SPLINE1(plog,Profile%Met%AltitudeEdge,Profile%lmx+1,y2)
    
    ! Loop over cloud pixels
    DO c=1,Profile%Cld%nSubPix

      ! Compute log pcld
      log_pcld(1) = LOG10(Profile%Cld%CloudPressure(c))

      ! Convert layer to height coordinates
      CALL SPLINT1(plog,Profile%Met%AltitudeEdge,y2,Profile%lmx+1,&
                  log_pcld,zcld,1)
      
      ! Determine Layer
      DO l=1,Profile%lmx

        ! merge levels within 100m
        IF(ABS(zcld(1)-Profile%Met%AltitudeEdge(l+1)) < 1.0d-1) THEN

          ! Set the level
          zcld(1) = Profile%Met%AltitudeEdge(l+1)
          Profile%Cld%LambertianCldLevel(c) = l 

          ! Full Column
          Profile%Cld%LambertianCldLayerFrac(c) = 1.0d0

          ! Exit Loop
          EXIT
        
        ! Cloud is now above level
        ELSEIF(zcld(1) .GT. Profile%Met%AltitudeEdge(l+1)) THEN

          ! Set the level
          Profile%Cld%LambertianCldLevel(c) = l

          ! Compute Column Fraction
          Profile%Cld%LambertianCldLayerFrac(c) =                      &
            (Profile%Cld%CloudPressure(c)-Profile%Met%PressureEdge(l)) &
          / (Profile%Met%PressureEdge(l+1)-Profile%Met%PressureEdge(l))

          ! Exit Loop
          EXIT

        ENDIF

      ENDDO

      ! Check if Cloud is below the surface
       IF( zcld(1) < Profile%Met%AltitudeEdge(Profile%lmx+1) ) THEN
          
          Profile%Cld%LambertianCldLevel(c) = Profile%lmx
          Profile%Cld%LambertianCldLayerFrac(c) = 1.0d0
          print*,'Warning - Cloud below surface (setting cloud bottom at surface)'
          
       ENDIF
       
    ENDDO
    
  END SUBROUTINE ComputeLambertianCloudLevel


  SUBROUTINE LoadProfileVar(ncid , varname, GridWt, OutVar, lmx, Error, ClimVertGrid, IsMixRatio_in)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,                       INTENT(IN)    :: ncid
    CHARACTER(LEN=*),              INTENT(IN)    :: varname
    TYPE(XYGridWtType),            INTENT(IN)    :: GridWt  
    INTEGER,                       INTENT(IN)    :: lmx
    REAL(KIND=8),                  INTENT(OUT)   :: OutVar(lmx)
    TYPE(ErrorType),               INTENT(INOUT) :: Error
    TYPE(ClimVgridType), OPTIONAL, INTENT(IN)    :: ClimVertGrid
    LOGICAL,             OPTIONAL, INTENT(IN)    :: IsMixRatio_in

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER                   :: rcode, vid, lmx_in, errstat
    REAL(KIND=8), ALLOCATABLE :: profvar(:), dp(:), orig_val(:), bgrd(:), resid(:)
    REAL(KIND=8)              :: dp_out(lmx)
    INTEGER                   :: dimid(4), dims(4), n, m, i, j, l, t, p
    INTEGER                   :: vmin(4,2), nv(4,2)
    LOGICAL                   :: unit_dim(4), IsMixRatio
    TYPE(ArrSetType_4D)       :: tmpvar
    !REAL(KIND=8), ALLOCATABLE :: tmpvar(:,:,:,:)
    INTEGER, ALLOCATABLE      :: idx(:,:), jdx(:), ldx(:), tdx(:)
    CHARACTER(LEN=maxChar)    :: tmpchar
    INTEGER(KIND=2)           :: SigmaAdjustType

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadProfileVar'

    ! ============================================================
    ! LoadProfileVar starts here
    ! ============================================================

    ! Get Vertical grid type
    IF(PRESENT(ClimVertGrid)) THEN
      lmx_in = ClimVertGrid%lmx
    ELSE
      lmx_in = lmx
    ENDIF
    
    ALLOCATE(orig_val(lmx_in))
    ALLOCATE(bgrd(lmx_in))
    ALLOCATE(profvar(lmx_in))
    ALLOCATE(dp(lmx_in))
    ALLOCATE(resid(lmx_in))

    ! Default assumption for vertical regrid
    IsMixratio = .TRUE. ; IF(PRESENT(IsMixRatio_in)) IsMixratio = IsMixRatio_in

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    CALL ReadProfileVarTmpArray(ncid,varname,GridWt,lmx_in,vmin,&
                                nv,idx,jdx,ldx,tdx,tmpvar,Error,&
                                SigmaAdjustType)
    
    ! Now Compute the averaged profile
    profvar(:) = 0.0d0
    DO t=1,GridWt%nt
    DO l=1,lmx_in
    DO j=1,Gridwt%nj
    DO i=1,Gridwt%ni
      profvar(l) = profvar(l) &
                 + GridWt%Weight(i,j,t)*tmpvar%data(1)%arr(idx(i,1),jdx(j),ldx(l),tdx(t))
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    IF(GridWt%PeriodicOverlap) THEN

      DO t=1,GridWt%nt
      DO l=1,lmx_in
      DO j=1,Gridwt%nj
      DO i=1,Gridwt%ni2
        profvar(l) = profvar(l) &
                  + GridWt%Weight2(i,j,t)*tmpvar%data(2)%arr(idx(i,2),jdx(j),ldx(l),tdx(t))
      ENDDO
      ENDDO
      ENDDO
      ENDDO

    ENDIF
    
    IF(PRESENT(ClimVertGrid)) THEN

      ! We are performing a vertical regrid
      OutVar(:) = 0.0d0

      IF(IsMixRatio) THEN

        ! Store profile on original grid
        orig_val(:) = profvar(:)

        ! First we must perform a terrain height adjustment to the profile
        ! to match the output grid
        SELECT CASE(SigmaAdjustType)

          CASE(2) ! Sigma (Preserve layer column densities)
            
            DO l=1,lmx_in
              profvar(l) = orig_val(l)*&
                (ClimVertGrid%PressureEdgeSigma(l+1)-ClimVertGrid%PressureEdgeSigma(l))/&
                (ClimVertGrid%PressureEdge(l+1)-ClimVertGrid%PressureEdge(l))
            ENDDO

          CASE(3) ! Hybrid method
            
            ! We Need to also load the background Profile
            CALL ReadProfileVarTmpArray(ncid,varname//'_Background',GridWt,lmx_in,vmin,&
                                        nv,idx,jdx,ldx,tdx,tmpvar,Error,&
                                        SigmaAdjustType)
            bgrd = 0.0d0
            DO t=1,GridWt%nt
            DO l=1,lmx_in
            DO j=1,Gridwt%nj
            DO i=1,Gridwt%ni
              bgrd(l) = bgrd(l) &
                        + GridWt%Weight(i,j,t)*tmpvar%data(1)%arr(idx(i,1),jdx(j),ldx(l),tdx(t))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            IF(GridWt%PeriodicOverlap) THEN
              DO t=1,GridWt%nt
              DO l=1,lmx_in
              DO j=1,Gridwt%nj
              DO i=1,Gridwt%ni2
                profvar(l) = profvar(l) &
                          + GridWt%Weight2(i,j,t)*tmpvar%data(2)%arr(idx(i,2),jdx(j),ldx(l),tdx(t))
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDIF

            ! First Interpolate background component to grid
            CALL LinearInt_Edgefill(ClimVertGrid%PressureMid,bgrd,lmx_in,       &
                                   ClimVertGrid%PressureMidSigma,profvar,lmx_in,&
                                   errstat                                      )

            ! Compute residual profile
            resid = orig_val - bgrd

            ! Add the sigma-adjusted residual to the background
            DO l=1,lmx_in
              profvar(l) = profvar(l) + resid(l)*&
                (ClimVertGrid%PressureEdgeSigma(l+1)-ClimVertGrid%PressureEdgeSigma(l))/&
                (ClimVertGrid%PressureEdge(l+1)-ClimVertGrid%PressureEdge(l))
            ENDDO

          CASE DEFAULT ! Assume background concentration as default

            CALL LinearInt_Edgefill(ClimVertGrid%PressureMid(1:lmx_in),orig_val(1:lmx_in),lmx_in,&
                                   ClimVertGrid%PressureMidSigma(1:lmx_in),profvar(1:lmx_in),lmx_in, &
                                   errstat                                        )
        END SELECT
        
        ! Compute Pressure Difference
        dp = ClimVertGrid%PressureEdgeSigma(2:ClimVertGrid%lmx+1) &
           - ClimVertGrid%PressureEdgeSigma(1:ClimVertGrid%lmx) 
        
        dp_out = 0.0d0
        DO j=1,lmx
        DO i=1,lmx_in
          OutVar(j) = OutVar(j) + ClimVertGrid%RegridWtSigma(i,j)*profvar(i)*dp(i)
          dp_out(j) = dp_out(j) + ClimVertGrid%RegridWtSigma(i,j)*dp(i)
        ENDDO
        ENDDO
        
        ! Undo
        OutVar(:) = OutVar(:) / dp_out(:)

      ELSE

        ! Assume quantity is already partial layer
        DO j=1,lmx
        DO i=1,lmx_in
          OutVar(j) = OutVar(j) + ClimVertGrid%RegridWt(i,j)*profvar(i)
        ENDDO
        ENDDO

      ENDIF
    ELSE

      ! We are directly rading the profile
      OutVar(:) = profvar(:)

    ENDIF

    DEALLOCATE(orig_val,bgrd)
    DEALLOCATE(profvar,dp,resid)
  END SUBROUTINE LoadProfileVar
  
  SUBROUTINE LoadCovariance(ncid , varname, GridWt, OutVar, lmx, Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: varname
    TYPE(XYGridWtType),         INTENT(IN)    :: GridWt  
    INTEGER,                    INTENT(IN)    :: lmx
    REAL(KIND=8),               INTENT(OUT)   :: OutVar(lmx,lmx)
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    
    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: rcode, vid
    INTEGER                   :: dimid(5), dims(5), n, m, i, j, l, k, t, p
    INTEGER                   :: vmin(5), nv(5)
    LOGICAL                   :: unit_dim(5)
    REAL(KIND=8), ALLOCATABLE :: tmpvar(:,:,:,:,:)
    INTEGER, ALLOCATABLE      :: idx(:), jdx(:), ldx(:),kdx(:), tdx(:)
    CHARACTER(LEN=maxChar)    :: tmpchar 

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadCovariance'

    ! ============================================================
    ! LoadProfileVar starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Initialize local min/max indices
    vmin(1) = GridWt%imin ; nv(1) = GridWt%ni
    vmin(2) = GridWt%imin ; nv(2) = GridWt%nj
    vmin(3) = 1           ; nv(3) = lmx
    vmin(4) = 1           ; nv(4) = lmx
    vmin(5) = GridWt%tmin ; nv(5) = GridWt%nt

    ! Allocate index arrays
    ALLOCATE(idx(nv(1))) ; ALLOCATE(jdx(nv(2)))
    ALLOCATE(ldx(nv(3))) ; ALLOCATE(kdx(nv(4)))
    ALLOCATE(tdx(nv(5)))
    
    ! Initialize Indices
    DO i=1,nv(1)
      idx(i) = i
    ENDDO
    DO j=1,nv(2)
      jdx(j) = j
    ENDDO
    DO l=1,nv(3)
      ldx(l) = l
    ENDDO
    DO k=1,nv(4)
      kdx(k) = k
    ENDDO
    DO t=1,nv(4)
      tdx(t) = t
    ENDDO
    
    ! Get dimensions of variable
    rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname)), vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
          'nf_inq_varid:'//TRIM(ADJUSTL(varname))                     )
    rcode = nf_inq_vardimid(ncid, vid, dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
          'nf_inq_vardimid:'//TRIM(ADJUSTL(varname))                  )
    DO n=1,5

      ! Load dimension
      rcode = nf_inq_dim(ncid, dimid(n),tmpchar,dims(n))
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,    &
          'nf_inq_dim:'//TRIM(ADJUSTL(varname))//'-'//TRIM(ADJUSTL(tmpchar)))

      ! Check if dimension is unit
      unit_dim(n) = dims(n) .EQ. 1

      ! Update load dimensions
      IF(unit_dim(n)) THEN
        vmin(n) = 1 ;  nv(n) = 1
      ENDIF

    ENDDO
    
    ! Set unit fields to 1
    IF(unit_dim(1)) idx(:) = 1
    IF(unit_dim(2)) jdx(:) = 1
    IF(unit_dim(3)) ldx(:) = 1
    IF(unit_dim(4)) kdx(:) = 1
    IF(unit_dim(5)) tdx(:) = 1

    ! Allocate variable arrays
    ALLOCATE(tmpvar(nv(1),nv(2),nv(3),nv(4),nv(5))) ; tmpvar(:,:,:,:,:) = 0.0d0

    ! Load the data
    rcode = nf_get_vara_double( ncid, vid, vmin, nv, tmpvar(:,:,:,:,:) )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
          'nf_get_vara_double:'//TRIM(ADJUSTL(varname))               )  

    ! Now Compute the averaged profile
    OutVar(:,:) = 0.0d0
    
    DO t=1,GridWt%nt
    DO k=1,lmx
    DO l=1,lmx
    DO j=1,Gridwt%nj
    DO i=1,Gridwt%ni
      OutVar(l,k) = OutVar(l,k) + GridWt%Weight(i,j,t)*tmpvar(idx(i),jdx(j),ldx(l),kdx(k),tdx(t))
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
  END SUBROUTINE LoadCovariance
  
  SUBROUTINE LoadVarUncert(ncid , varname, GridWt, OutVarUncert, Profile, Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: varname
    TYPE(XYGridWtType),         INTENT(IN)    :: GridWt
    TYPE(UncertType),           INTENT(INOUT) :: OutVarUncert
    TYPE(ProfileType),          INTENT(IN)    :: Profile
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    
    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: rcode, vid, xmx, errstat
    INTEGER                   :: dimid(4), dims(4), n, m, i, j, l, k, t, p
    INTEGER                   :: vmin(4), nv(4)
    LOGICAL                   :: unit_dim(4)
    REAL(KIND=8), ALLOCATABLE :: tmpvar(:,:,:,:)
    INTEGER, ALLOCATABLE      :: idx(:), jdx(:), ldx(:), tdx(:), lmx
    CHARACTER(LEN=maxChar)    :: tmpchar
    REAL(KIND=8)              :: dz
    TYPE(UncertType)          :: VarUncert_Clim
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadVarUncert'

    ! ============================================================
    ! LoadVarUncert starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    IF(Profile%UseL2Met) THEN
      lmx = Profile%ClimVertGrid%lmx
    ELSE
      lmx = Profile%lmx
    ENDIF

    ! Read State Type
    rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname))//'_StateType', vid)
    rcode = nf_inq_vardimid(ncid, vid, dimid)
    rcode = nf_get_vara_int(ncid, vid, (/1/), (/1/), OutVarUncert%StateType)
    
    ! Read Error Covariance Parameterization Type
    rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname))//'_ErrorCovarType', vid)
    rcode = nf_inq_vardimid(ncid, vid, dimid)
    rcode = nf_get_vara_int(ncid, vid, (/1/), (/1/), OutVarUncert%CovarParType)
    
    ! Get the state vector dimension
    xmx = GetOutputCovarDim(OutVarUncert,lmx)

    ! Check allocation
    IF(ALLOCATED(OutVarUncert%ScalarPar))    DEALLOCATE(OutVarUncert%ScalarPar)
    IF(ALLOCATED(OutVarUncert%ProfilePar))   DEALLOCATE(OutVarUncert%ProfilePar)
    IF(ALLOCATED(OutVarUncert%SubCovMatrix)) DEALLOCATE(OutVarUncert%SubCovMatrix)
    
    ! Allocate/Zero Output sub-covariance matrix
    ALLOCATE(OutVarUncert%SubCovMatrix(xmx,xmx))
    OutVarUncert%SubCovMatrix(:,:) = 0.0d0
    
    ! Load the uncertainty
    ! --------------------
    ! Has not been specified ==> Set negative value
    IF(OutVarUncert%CovarParType .LT. 0) THEN
      
      ! Load the profile parameter 
      ALLOCATE(OutVarUncert%ProfilePar(xmx))

      DO l=1,xmx
          OutVarUncert%ProfilePar(l) = -1.0d0
          OutVarUncert%SubCovMatrix(l,l) = -1.0d0
      ENDDO

    ! Covariance Matrix for footprint variables
    ELSEIF(OutVarUncert%CovarParType .EQ. 0) THEN

      ! Check its not a profile variable
      IF(OutVarUncert%StateType .GT. 0) THEN
        CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,      &
                              Message_in='Error: trying to assign a scalar covariance '  &
                                        //'to profile variable:'//TRIM(ADJUSTL(varname)),&
                              Action_in=''                                               )
      ENDIF

      ! Load scalar uncertainty
      CALL LoadSurfaceVar(ncid , TRIM(ADJUSTL(varname))//'_ErrorCovar', &
                          GridWt,OutVarUncert%SubCovMatrix(1,1), Error  )
    
    ! Fully Specified Covariance Matrix
    ELSEIF( OutVarUncert%CovarParType .EQ. 3 ) THEN
      
      ! Directly read sub covariance matrix from disk
      CALL LoadCovariance(ncid, TRIM(ADJUSTL(varname))//'_ErrorCovar',&
                          GridWt,OutVarUncert%SubCovMatrix(:,:), xmx, Error)


      ! CCM FIX - IF we are vertically regridding this needs to be reweighted!!!!
      ! Stop code for this condition for now
      IF(Profile%UseL2Met) THEN
        STOP 'Fully specfied error covar matrix must be reweighted ' &
             //'when using L1 met (not currently implemented'

      ENDIF

    ! State-type dependent parameterizations
    ELSE
      
      ! Load the profile parameter 
      ALLOCATE(OutVarUncert%ProfilePar(xmx))
      
      ! -----------------
      ! Profile Variables
      ! -----------------
      IF(OutVarUncert%StateType .GE. 1 .AND. OutVarUncert%StateType .LE. 4) THEN
        
        CALL LoadProfileVar(ncid , TRIM(ADJUSTL(varname))//'_ErrorCovar', GridWt,&
                            OutVarUncert%ProfilePar, xmx, Error)
        
        ! Diagonal Matrix parameterization
        IF( OutVarUncert%CovarParType .EQ. 1) THEN
          
          DO l=1,xmx
            OutVarUncert%SubCovMatrix(l,l) = OutVarUncert%ProfilePar(l)**2
          ENDDO
        
        ! Altitude correlated prior parameterization
        ELSEIF(OutVarUncert%CovarParType .EQ. 2) THEN
          
          ! Load correlation length scale
          ALLOCATE(OutVarUncert%ScalarPar(1))
          CALL LoadSurfaceVar(ncid , TRIM(ADJUSTL(varname))//'_CorrelationLength', GridWt, &
                              OutVarUncert%ScalarPar(1), Error)
          
          IF(xmx .NE. Profile%lmx) THEN
            CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,   &
                                  Message_in='Attempting to compute ZCORRELATED Profile ' &
                                           //'incorrectly dimensioned matrix',            &
                                  Action_in=''                                            )
          ENDIF

          ! This is computed after interpolating the parameters to the new grid
          IF(.NOT. Profile%UseL2Met) THEN
            CALL ComputeZCorrCovar(xmx,OutVarUncert%ProfilePar,Profile%Met%AltitudeMid,&
                                  OutVarUncert%ScalarPar(1),OutVarUncert%SubCovMatrix)
          ENDIF
        ENDIF
      

      ELSE
        
        print*,'Unknown state element type for ',TRIM(ADJUSTL(varname))
        STOP 1
        
      ENDIF
      
    ENDIF

    ! We now need to check if we need to do a vertical regrid
    ! StateType > 0 => Profile
    IF((OutVarUncert%StateType .GE. 1 .AND. OutVarUncert%StateType .LE. 4) .AND. Profile%UseL2Met) THEN
      
      ! Store Copy of out variable uncertainty to regrid
      VarUncert_Clim = OutVarUncert

      ! Reallocate to L2 Met Dimension
      IF(ALLOCATED(OutVarUncert%ProfilePar)  ) DEALLOCATE(OutVarUncert%ProfilePar)
      IF(ALLOCATED(OutVarUncert%SubCovMatrix)) DEALLOCATE(OutVarUncert%SubCovMatrix)
      ALLOCATE(OutVarUncert%ProfilePar(Profile%lmx)) ; OutVarUncert%ProfilePar = 0.0d0
      ALLOCATE(OutVarUncert%SubCovMatrix(Profile%lmx,Profile%lmx)) ; OutVarUncert%SubCovMatrix = 0.0d0

      ! Interpolate uncertainty Profile
      CALL BSPLINE_Edgefill(Profile%ClimVertGrid%PressureMid,VarUncert_Clim%ProfilePar,lmx,&
                            Profile%Met%PressureMid,OutVarUncert%ProfilePar,Profile%lmx,   &
                            errstat                                                        )

      ! Diagonal Matrix parameterization
      IF(OutVarUncert%CovarParType .LT. 1) THEN
        DO l=1,Profile%lmx
          OutVarUncert%SubCovMatrix(l,l) = -1.0d0
        ENDDO

      ELSEIF( OutVarUncert%CovarParType .EQ. 1) THEN
        
        DO l=1,Profile%lmx
          OutVarUncert%SubCovMatrix(l,l) = OutVarUncert%ProfilePar(l)**2
        ENDDO

      ! Altitude correlated prior parameterization
      ELSEIF(OutVarUncert%CovarParType .EQ. 2) THEN

        CALL ComputeZCorrCovar(Profile%lmx,OutVarUncert%ProfilePar,Profile%Met%AltitudeMid,&
                               OutVarUncert%ScalarPar(1),OutVarUncert%SubCovMatrix)

      ENDIF

    ENDIF
    
  END SUBROUTINE LoadVarUncert
  
  SUBROUTINE ReadSurfaceVarTmpArray(ncid,varname,GridWt,vmin,nv,idx,jdx,tdx,tmpvar,Error,ReadI2_in)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,                     INTENT(IN)    :: ncid
    CHARACTER(LEN=*),            INTENT(IN)    :: varname
    TYPE(XYGridWtType),          INTENT(IN)    :: GridWt
    INTEGER,                     INTENT(INOUT) :: vmin(3,2)
    INTEGER,                     INTENT(INOUT) :: nv(3,2)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: idx(:,:)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: jdx(:)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: tdx(:)
    TYPE(ArrSetType_3D),         INTENT(INOUT) :: tmpvar 
    TYPE(ErrorType),             INTENT(INOUT) :: Error
    LOGICAL, OPTIONAL,           INTENT(IN)    :: ReadI2_in

    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: vmin_vec(3), nv_vec(3)
    INTEGER                   :: rcode, vid, imx
    INTEGER                   :: dimid(3), dims(3), n, m, i, j, t
    LOGICAL                   :: unit_dim(3)
    CHARACTER(LEN=maxChar)    :: tmpchar
    LOGICAL                   :: ReadI2
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ReadSurfaceVarTmpArray'

    ! ============================================================
    ! ReadSurfaceVarTmpArray starts here
    ! ============================================================

    ! Set default (F=> Read Double)
    ReadI2 = .FALSE. ; IF(PRESENT(ReadI2_in)) ReadI2 = ReadI2_in

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Deallocate memory for temporary array
    IF(ALLOCATED(tmpvar%Data)) DEALLOCATE(tmpvar%Data)
    IF(ALLOCATED(tmpvar%imx)) DEALLOCATE(tmpvar%imx)
    IF(ALLOCATED(tmpvar%jmx)) DEALLOCATE(tmpvar%jmx)
    IF(ALLOCATED(tmpvar%kmx)) DEALLOCATE(tmpvar%kmx)

    ! Determine Array Allocation
    tmpvar%n = 1 ; IF(GridWt%PeriodicOverlap) tmpvar%n = 2
    ALLOCATE(tmpvar%imx(tmpvar%n)) ; ALLOCATE(tmpvar%jmx(tmpvar%n)) ; 
    ALLOCATE(tmpvar%kmx(tmpvar%n)) ; ALLOCATE(tmpvar%Data(tmpvar%n))

    ! Initialize local min/max indices
    vmin(:,:) = 0 ; nv(:,:) = 0
    vmin(1,1) = GridWt%imin ; nv(1,1) = GridWt%ni
    vmin(2,1) = GridWt%jmin ; nv(2,1) = GridWt%nj
    vmin(3,1) = GridWt%tmin ; nv(3,1) = GridWt%nt
    IF(GridWt%PeriodicOverlap) THEN
      vmin(1,2) = GridWt%imin2 ; nv(1,2) = GridWt%ni2
      vmin(2,2) = GridWt%jmin ; nv(2,2) = GridWt%nj
      vmin(3,2) = GridWt%tmin ; nv(3,2) = GridWt%nt
    ENDIF

    ! Check allocation
    IF(ALLOCATED(idx)) DEALLOCATE(idx)
    IF(ALLOCATED(jdx)) DEALLOCATE(jdx)
    IF(ALLOCATED(tdx)) DEALLOCATE(tdx)

    ! Maximimum i dimension
    imx = MAXVAL(nv(1,:))
    
    ! Allocate index arrays
    ALLOCATE(idx(imx,tmpvar%n))
    ALLOCATE(jdx(nv(2,1))) ; ALLOCATE(tdx(nv(3,1)))
    
    ! Initialize Indices
    DO i=1,imx
      idx(i,:) = i
    ENDDO
    DO j=1,nv(2,1)
      jdx(j) = j
    ENDDO
    DO t=1,nv(3,1)
      tdx(t) = t
    ENDDO
    
    ! Get dimensions of variable
    rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname)), vid)
    rcode = nf_inq_vardimid(ncid, vid, dimid)

    DO n=1,3

      ! Load dimension
      rcode = nf_inq_dim(ncid, dimid(n),tmpchar,dims(n))
      
      ! Check if dimension is unit
      unit_dim(n) = dims(n) .EQ. 1

      ! Update load dimensions
      IF(unit_dim(n)) THEN
        vmin(n,:) = 1 ;  nv(n,:) = 1
      ENDIF

    ENDDO
    
    ! Set unit fields to 1
    IF(unit_dim(1)) idx(:,:) = 1
    IF(unit_dim(2)) jdx(:)   = 1
    IF(unit_dim(3)) tdx(:)   = 1

    ! Load data
    DO n=1,tmpvar%n

      ! Block to load
      vmin_vec(:) = vmin(:,n) ; nv_vec(:) = nv(:,n)

      

      IF(ReadI2) THEN

        ! Allocate variable arrays
        ALLOCATE( tmpvar%data(n)%arr_i2(nv_vec(1),nv_vec(2),nv_vec(3)) )
        
        ! Zero array
        tmpvar%data(n)%arr_i2(:,:,:) = 0.0d0

        ! Load the data
        rcode = nf_get_vara_double( ncid, vid, vmin_vec, nv_vec, tmpvar%data(n)%arr_i2(:,:,:) )

      ELSE

        ! Allocate variable arrays
        ALLOCATE(tmpvar%data(n)%arr(nv_vec(1),nv_vec(2),nv_vec(3)))
        
        ! Zero array
        tmpvar%data(n)%arr(:,:,:) = 0.0d0

        ! Load the data
        rcode = nf_get_vara_double( ncid, vid, vmin_vec, nv_vec, tmpvar%data(n)%arr(:,:,:) )

      ENDIF

    ENDDO

  END SUBROUTINE ReadSurfaceVarTmpArray
  
  
  SUBROUTINE ReadProfileVarTmpArray(ncid,varname,GridWt,lmx,vmin,nv,idx,jdx,ldx,tdx,tmpvar,Error,SigmaAdjustType)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,                     INTENT(IN)    :: ncid
    CHARACTER(LEN=*),            INTENT(IN)    :: varname
    TYPE(XYGridWtType),          INTENT(IN)    :: GridWt
    INTEGER,                     INTENT(IN)    :: lmx
    INTEGER,                     INTENT(INOUT) :: vmin(4,2)
    INTEGER,                     INTENT(INOUT) :: nv(4,2)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: idx(:,:)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: jdx(:)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: ldx(:)
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: tdx(:)
    TYPE(ArrSetType_4D),         INTENT(INOUT) :: tmpvar
    !REAL(KIND=8),   ALLOCATABLE, INTENT(INOUT) :: tmpvar(:,:,:,:)
    
    TYPE(ErrorType),             INTENT(INOUT) :: Error 
    INTEGER(KIND=2), OPTIONAL,   INTENT(OUT)   :: SigmaAdjustType
    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: rcode, vid, imx
    INTEGER                   :: dimid(4), dims(4), n, m, i, j, l, t
    LOGICAL                   :: unit_dim(4)
    CHARACTER(LEN=maxChar)    :: tmpchar 
    INTEGER                   :: vmin_vec(4), nv_vec(4)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ReadProfileVarTmpArray'

    ! ============================================================
    ! ReadProfileVarTmpArraystarts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Deallocate memory for temporary array
    IF(ALLOCATED(tmpvar%Data)) DEALLOCATE(tmpvar%Data)
    IF(ALLOCATED(tmpvar%imx)) DEALLOCATE(tmpvar%imx)
    IF(ALLOCATED(tmpvar%jmx)) DEALLOCATE(tmpvar%jmx)
    IF(ALLOCATED(tmpvar%kmx)) DEALLOCATE(tmpvar%kmx)
    IF(ALLOCATED(tmpvar%lmx)) DEALLOCATE(tmpvar%lmx)

    ! Determine Array Allocation
    tmpvar%n = 1 ; IF(GridWt%PeriodicOverlap) tmpvar%n = 2
    ALLOCATE(tmpvar%imx(tmpvar%n)) ; ALLOCATE(tmpvar%jmx(tmpvar%n)) ; 
    ALLOCATE(tmpvar%kmx(tmpvar%n)) ; ALLOCATE(tmpvar%lmx(tmpvar%n))
    ALLOCATE(tmpvar%Data(tmpvar%n))

    ! Initialize local min/max indices
    vmin(:,:) = 0 ; nv(:,:) = 0
    vmin(1,1) = GridWt%imin ; nv(1,1) = GridWt%ni
    vmin(2,1) = GridWt%jmin ; nv(2,1) = GridWt%nj
    vmin(3,1) =           1 ; nv(3,1) = lmx
    vmin(4,1) = GridWt%tmin ; nv(4,1) = GridWt%nt
    IF(GridWt%PeriodicOverlap) THEN
      vmin(1,2) = GridWt%imin2 ; nv(1,2) = GridWt%ni2
      vmin(2,2) = GridWt%jmin  ; nv(2,2) = GridWt%nj
      vmin(3,2) =           1  ; nv(3,2) = lmx
      vmin(4,2) = GridWt%tmin  ; nv(4,2) = GridWt%nt
    ENDIF

    ! Check allocation
    IF(ALLOCATED(idx)) DEALLOCATE(idx)
    IF(ALLOCATED(jdx)) DEALLOCATE(jdx)
    IF(ALLOCATED(ldx)) DEALLOCATE(ldx)
    IF(ALLOCATED(tdx)) DEALLOCATE(tdx)

    ! Maximimum i dimension
    imx = MAXVAL(nv(1,:))
    
    ! Allocate index arrays
    ALLOCATE(idx(imx,tmpvar%n)) ; ALLOCATE(jdx(nv(2,1))) 
    ALLOCATE(ldx(nv(3,1)))      ; ALLOCATE(tdx(nv(3,1)))

    ! Initialize Indices
    DO i=1,imx
      idx(i,:) = i
    ENDDO
    DO j=1,nv(2,1)
      jdx(j) = j
    ENDDO
    DO l=1,nv(3,1)
      ldx(l) = l
    ENDDO
    DO t=1,nv(4,1)
      tdx(t) = t
    ENDDO

    ! Get dimensions of variable
    rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname)), vid)
    CALL netcdf_handle_error('Error: NetCDF could not get dim for '//TRIM(ADJUSTL(varname)),rcode)
    rcode = nf_inq_vardimid(ncid, vid, dimid)
    CALL netcdf_handle_error('Error: NetCDF could not get dim for '//TRIM(ADJUSTL(varname)),rcode)
    
    DO n=1,4

      ! Load dimension
      rcode = nf_inq_dim(ncid, dimid(n),tmpchar,dims(n))
      CALL netcdf_handle_error('Error: NetCDF could not get dim for '//TRIM(ADJUSTL(varname)),rcode)
      
      ! Check if dimension is unit
      unit_dim(n) = dims(n) .EQ. 1

      ! Update load dimensions
      IF(unit_dim(n)) THEN
        vmin(n,:) = 1 ;  nv(n,:) = 1
      ENDIF

    ENDDO
    
    ! Set unit fields to 1
    IF(unit_dim(1)) idx(:,:) = 1
    IF(unit_dim(2)) jdx(:)   = 1
    IF(unit_dim(3)) ldx(:)   = 1
    IF(unit_dim(4)) tdx(:)   = 1

    DO n=1,tmpvar%n

      ! Block to load
      vmin_vec(:) = vmin(:,n) ; nv_vec(:) = nv(:,n)

      ! Allocate variable arrays
      ALLOCATE(tmpvar%data(n)%arr(nv_vec(1),nv_vec(2),nv_vec(3),nv_vec(4) ))
      
      ! Zero array
      tmpvar%data(n)%arr(:,:,:,:) = 0.0d0

      ! Load the data
      rcode = nf_get_vara_double( ncid, vid, vmin_vec, nv_vec, tmpvar%data(n)%arr(:,:,:,:) )
      CALL netcdf_handle_error('Error: NetCDF could not read '//TRIM(ADJUSTL(varname)),rcode)

    ENDDO
    
    ! Output Sigma adjust type if requested
    IF(PRESENT(SigmaAdjustType)) THEN
      SigmaAdjustType = -1
      rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname))//'_SigmaAdjustType', vid)

      ! If variable has a sigma adjust type return it
      IF(rcode .EQ. NF_NOERR ) THEN
        rcode = nf_get_vara_double( ncid, vid, (/1/), (/1/),SigmaAdjustType)
      ENDIF

    ENDIF

  END SUBROUTINE ReadProfileVarTmpArray

  SUBROUTINE LoadSurfaceVar(ncid , varname, GridWt, OutVar, Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,               INTENT(IN)    :: ncid
    CHARACTER(LEN=*),      INTENT(IN)    :: varname
    TYPE(XYGridWtType),    INTENT(IN)    :: GridWt  
    REAL(KIND=8),          INTENT(OUT)   :: OutVar
    TYPE(ErrorType),       INTENT(INOUT) :: Error

    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: n, m, i, j, t
    INTEGER                   :: vmin(3,2), nv(3,2)
    !REAL(KIND=8), ALLOCATABLE :: tmpvar(:,:,:)
    TYPE(ArrSetType_3D)       :: tmpvar 
    INTEGER, ALLOCATABLE      :: idx(:,:), jdx(:), tdx(:)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadSurfaceVar'

    ! ============================================================
    ! LoadSurfaceVar starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
    
    CALL ReadSurfaceVarTmpArray(ncid,varname,GridWt,vmin,nv,idx,jdx,tdx,tmpvar,Error)
    
    ! Now Compute the averaged profile
    OutVar= 0.0d0
    DO t=1,GridWt%nt
    DO j=1,Gridwt%nj
    DO i=1,Gridwt%ni
      OutVar = OutVar + GridWt%Weight(i,j,t)*tmpvar%data(1)%arr(idx(i,1),jdx(j),tdx(t))
    ENDDO
    ENDDO
    ENDDO

    ! Include second subpixel if needed
    IF(GridWt%PeriodicOverlap) THEN

      DO t=1,GridWt%nt
      DO j=1,Gridwt%nj
      DO i=1,Gridwt%ni2
        OutVar = OutVar + GridWt%Weight2(i,j,t)*tmpvar%data(2)%arr(idx(i,2),jdx(j),tdx(t))
      ENDDO
      ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE LoadSurfaceVar
  
  SUBROUTINE LoadXYVar(ncid , varname, GridWt, OutVar, Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,               INTENT(IN)    :: ncid
    CHARACTER(LEN=*),      INTENT(IN)    :: varname
    TYPE(XYGridWtType),    INTENT(IN)    :: GridWt  
    REAL(KIND=8),          INTENT(OUT)   :: OutVar
    TYPE(ErrorType),       INTENT(INOUT) :: Error

    ! --------------------
    ! Local Variables
    ! --------------------
    INTEGER                   :: rcode, vid
    INTEGER                   :: dimid(2), dims(2), n, m, i, j, t
    INTEGER                   :: vmin(2), nv(2)
    LOGICAL                   :: unit_dim(2)
    REAL(KIND=8), ALLOCATABLE :: tmpvar(:,:)
    INTEGER, ALLOCATABLE      :: idx(:), jdx(:), tdx(:)
    CHARACTER(LEN=maxChar)    :: tmpchar 

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadXYVar'

    ! ============================================================
    ! LoadXYVar starts here
    ! ============================================================

    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Initialize local min/max indices
    vmin(1) = GridWt%imin ; nv(1) = GridWt%ni
    vmin(2) = GridWt%imin ; nv(2) = GridWt%nj
    
    ! Allocate index arrays
    ALLOCATE(idx(nv(1))) ; ALLOCATE(jdx(nv(2)))
    
    ! Initialize Indices
    DO i=1,nv(1)
      idx(i) = i
    ENDDO
    DO j=1,nv(2)
      jdx(j) = j
    ENDDO
    
    ! Get dimensions of variable
    rcode = nf_inq_varid(ncid, TRIM(ADJUSTL(varname)), vid)
    rcode = nf_inq_vardimid(ncid, vid, dimid)

    DO n=1,2

      ! Load dimension
      rcode = nf_inq_dim(ncid, dimid(n),tmpchar,dims(n))
      
      ! Check if dimension is unit
      unit_dim(n) = dims(n) .EQ. 1

      ! Update load dimensions
      IF(unit_dim(n)) THEN
        vmin(n) = 1 ;  nv(n) = 1
      ENDIF

    ENDDO
    
    ! Set unit fields to 1
    IF(unit_dim(1)) idx(:) = 1
    IF(unit_dim(2)) jdx(:) = 1

    ! Allocate variable arrays
    ALLOCATE(tmpvar(nv(1),nv(2))) ; tmpvar(:,:) = 0.0d0

    ! Load the data
    rcode = nf_get_vara_double( ncid, vid, vmin, nv, tmpvar(:,:) )
    
    ! Now Compute the averaged profile
    OutVar= 0.0d0
    DO t=1,GridWt%nt
    DO j=1,Gridwt%nj
    DO i=1,Gridwt%ni
      OutVar = OutVar + GridWt%Weight(i,j,t)*tmpvar(idx(i),jdx(j))
    ENDDO
    ENDDO
    ENDDO

    ! Check if the second subset needs to be included
    IF(GridWt%PeriodicOverlap) THEN

      ! Update x indices
      vmin(1) = GridWt%imin2 ; nv(1) = GridWt%ni2
      DEALLOCATE(idx) ; ALLOCATE(idx(nv(1)))
      IF(unit_dim(1)) THEN
        idx(:) = 1
      ELSE
        DO i=1,nv(1)
          idx(i) = i
        ENDDO
      ENDIF

      ! Allocate variable arrays
      DEALLOCATE(tmpvar) ; ALLOCATE(tmpvar(nv(1),nv(2))) ; tmpvar(:,:) = 0.0d0

      ! Load the data
      rcode = nf_get_vara_double( ncid, vid, vmin, nv, tmpvar(:,:) )

      ! Add contribution from second sub pixel to average
      DO t=1,GridWt%nt
      DO j=1,Gridwt%nj
      DO i=1,Gridwt%ni2
        OutVar = OutVar + GridWt%Weight2(i,j,t)*tmpvar(idx(i),jdx(j))
      ENDDO
      ENDDO
      ENDDO

    ENDIF

    ! Deallocae temporary variable
    DEALLOCATE(tmpvar)

  END SUBROUTINE LoadXYVar

  SUBROUTINE ComputeProfileDerivedQuantities( Profile, Error, RecomputeZ_in )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ProfileType),     INTENT(INOUT) :: Profile
    TYPE(ErrorType),       INTENT(INOUT) :: Error
    LOGICAL, OPTIONAl,     INTENT(IN)    :: RecomputeZ_in
    
    ! ---------------
    ! Local Variables
    ! ---------------
    REAL(KIND=8), ALLOCATABLE :: tmp_prof(:), tmp_edge(:)
    INTEGER                   :: n, nr, g
    REAL(KIND=8)              :: g0, dp, pl, dp_sub
    LOGICAL                   :: RecomputeZ
    INTEGER, PARAMETER        :: npart_sublay = 1000 ! For the partial column calculation

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ComputeProfileDerivedQuantities'

    ! ============================================================
    ! SampleProfile starts here
    ! ============================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Default for optional arguments
    RecomputeZ = .TRUE. ; IF(PRESENT(RecomputeZ_in)) RecomputeZ = RecomputeZ_in
    
    ! Allocate temporary profile array
    ALLOCATE(tmp_prof(Profile%lmx))
    ALLOCATE(tmp_edge(Profile%lmx+1))
    
    ! ------------------------------------------------------------
    ! Air Molecular Weight
    ! ------------------------------------------------------------
    
    ! Zero air molecular weight
    Profile%Met%AirMolecularWeight(:) = 0.0d0
    
    ! Array to hold summed mixing ratio
    tmp_prof(:) = 0.0d0
    
    DO n=1,Profile%Gas%nSpecies

      ! Mixing ratio weighted molecular weight
      Profile%Met%AirMolecularWeight(:) = Profile%Met%AirMolecularWeight(:) &
                                        + Profile%Gas%MixingRatio(:,n)*Profile%Gas%MolecularWeight(n)
      tmp_prof(:) = tmp_prof(:) + Profile%Gas%MixingRatio(:,n)
      
    ENDDO

    ! Complete average
    Profile%Met%AirMolecularWeight(:) = Profile%Met%AirMolecularWeight(:) / tmp_prof(:)
    
    ! ------------------------------------------------------------
    ! Altitude + Gravity + Layer Density
    ! ------------------------------------------------------------
    
    ! Gravity at lower layer
    g0 = Profile%StandardGravity*(Profile%PlanetaryRadius &
       /(Profile%PlanetaryRadius + Profile%Met%AltitudeEdge(Profile%lmx+1)))**2
    
    ! Compute altitude edges
    DO nr=1, Profile%lmx
      
      ! Loop from bottom to TOA
      n = Profile%lmx - nr + 1
      
      IF(RecomputeZ) THEN

        ! If top layer is zero, compute height at pressure x1e-2 of the layer below
        IF(n .EQ. 1 .AND. Profile%Met%PressureEdge(n) .LT. TINY(0.0d0)) THEN
          ! Compute height at edge
          Profile%Met%AltitudeEdge(n) = Profile%Met%AltitudeEdge(n+1) &
                                    + Constants_R*Profile%Met%TemperatureEdge(n+1) &
                                    /(Profile%Met%AirMolecularWeight(n)*g0)       &
                                    *DLOG( Profile%Met%PressureEdge(n+1)*1e2 )*1e-3

          ! Compute height at midpoint
          Profile%Met%AltitudeMid(n) = Profile%Met%AltitudeEdge(n+1) &
                                      + Constants_R*Profile%Met%TemperatureEdge(n+1) &
                                      /(Profile%Met%AirMolecularWeight(n)*g0)     &
                                      *DLOG( Profile%Met%PressureMid(n)/(Profile%Met%PressureEdge(n+1)*1e-2) )*1e-3

          ! Raise a warning 
          CALL RaiseWarning( Error, ErrorCode_Profile, ModuleName, SubroutineName,                 &
                             'TOA Pressure is 0: Using 0.01xpenultimate layer for TOA height calc',&
                             'Change climatology (or L2 Met) top layer pressure > 0'               )

        ELSE 

          ! Compute height at edge
          Profile%Met%AltitudeEdge(n) = Profile%Met%AltitudeEdge(n+1) &
                                      + Constants_R*Profile%Met%TemperatureEdge(n+1) &
                                      /(Profile%Met%AirMolecularWeight(n)*g0)       &
                                      *DLOG( Profile%Met%PressureEdge(n+1)/Profile%Met%PressureEdge(n) )*1e-3
          
          ! Compute height at midpoint
          Profile%Met%AltitudeMid(n) = Profile%Met%AltitudeEdge(n+1) &
                                      + Constants_R*Profile%Met%TemperatureEdge(n+1) &
                                      /(Profile%Met%AirMolecularWeight(n)*g0)     &
                                      *DLOG( Profile%Met%PressureMid(n)/Profile%Met%PressureEdge(n) )*1e-3

        ENDIF
        
        
        ! Compute Gravity at midpoint
        Profile%Met%Gravity(n) = Profile%StandardGravity*(Profile%PlanetaryRadius &
                              /(Profile%PlanetaryRadius + Profile%Met%AltitudeMid(n)))**2
      
      ENDIF

      ! Update g
      g0 = Profile%StandardGravity*(Profile%PlanetaryRadius &
         /(Profile%PlanetaryRadius + Profile%Met%AltitudeEdge(n)))**2
      
      ! Layer pressure difference in Pa
      dp = 1.0d2*(Profile%Met%PressureEdge(n+1)-Profile%Met%PressureEdge(n))
      
      ! Air Layer Density (molec/cm2)
      Profile%Met%AirPartialColumn(n) = dp*Constants_NA*1.0d-4 & 
                                        /(Profile%Met%AirMolecularWeight(n)*Profile%Met%Gravity(n))
      
      ! Compute Partial Column squared for CIA - Break into 100 sub columns to ensure 
      CALL ComputePartColSquared(Profile%Met%PressureEdge(n+1),Profile%Met%PressureEdge(n),&
                                 Profile%Met%AltitudeEdge(n+1),Profile%Met%AltitudeEdge(n),&
                                 Profile%Met%AirMolecularWeight(n),Profile%Met%Gravity(n), &
                                 npart_sublay,Profile%Met%AirPartialColumnSquared(n)       )
      
    ENDDO
    
    ! ------------------------------------------------------------
    ! Individual Gas Layer Densities
    ! ------------------------------------------------------------
    
    DO n=1,Profile%Gas%nSpecies
      Profile%Gas%PartialColumn(:,n) = Profile%Met%AirPartialColumn(:)*Profile%Gas%MixingRatio(:,n)
    ENDDO

    ! ------------------------------------------------------------
    ! Compute Dry Air Column
    ! ------------------------------------------------------------
    
    ! Initialize with the full column
    Profile%Met%DryAirPartialColumn = Profile%Met%AirPartialColumn

    ! Subtract wet columns
    DO n=1,Profile%Gas%nWetAirSpc
      g = Profile%Gas%WetAirIdx(n)
      Profile%Met%DryAirPartialColumn(:) = Profile%Met%DryAirPartialColumn(:) &
                                         - Profile%Gas%PartialColumn(:,g)
    ENDDO


    ! Deallocate help arrays
    DEALLOCATE(tmp_prof)
    DEALLOCATE(tmp_edge)
    
  END SUBROUTINE ComputeProfileDerivedQuantities
  
  SUBROUTINE ComputePartColSquared(p_bot,p_top,z_bot,z_top,air_mw,g,nsublay,pcol_sq)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL(KIND=8), INTENT(IN)  :: p_bot   ! Pressure at bottom of layer
    REAL(KIND=8), INTENT(IN)  :: p_top   ! Pressure at top of layer 
    REAL(KIND=8), INTENT(IN)  :: z_bot   ! Height at bottom of layer
    REAL(KIND=8), INTENT(IN)  :: z_top   ! Height at top of layer (z_top>z_bot)
    REAL(KIND=8), INTENT(IN)  :: air_mw  ! Molecular weight of air
    REAL(KIND=8), INTENT(IN)  :: g       ! Layer gravity (assumed constant)
    INTEGER,      INTENT(IN)  :: nsublay ! Number of sub layers for calculation
    REAL(KIND=8), INTENT(OUT) :: pcol_sq ! The integrated air density squared

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: dp, dair_sq, H, dz, pb, pt, ztot
    INTEGER      :: l

    ! =====================================================================
    ! ComputePartColSquared starts here
    ! =====================================================================

    ! Compute scale height [cm]
    H = (z_top-z_bot)/DLOG(p_bot/p_top)*1e5
    
    ! Partial pressure for each sub layer [hPa]
    dp = (p_bot-p_top)/REAL(nsublay,KIND=8)
    dair_sq = ( 1.0e2*dp*Constants_NA*1.0d-4/(air_mw*g) )**2

    ! Initialize bottom of layer height/pressure
    pb = p_bot

    ! Initialize partial column
    pcol_sq = 0.0d0
    ztot = 0.0
    DO l=1,nsublay
        
        ! Compute top of sublayer pressure
        pt = pb - dp
        
        ! Change in height in cm
        dz = H*DLOG(pb/pt)

        ! Update partial column
        pcol_sq = pcol_sq + dair_sq / dz

        ! Set next bottom layer
        pb = pt

    ENDDO

  END SUBROUTINE ComputePartColSquared

  SUBROUTINE InitProfileDiags(nc_dimlist, nc_ndimlist, ProfDiagOpt, ProfDiag, Error)
    
    ! Allocates necessary arrays to track diagnostic quantities
    ! currently a placeholder
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(NCDimType), DIMENSION(nc_ndimlist), INTENT(IN)    :: nc_dimlist
    INTEGER,                                 INTENT(IN)    :: nc_ndimlist
    TYPE(ProfDiagOptType),                   INTENT(IN)    :: ProfDiagOpt
    TYPE(ProfDiagType),                      INTENT(INOUT) :: ProfDiag
    TYPE(ErrorType),                         INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER, DIMENSION(1) :: DIM_ID, DIMS
    LOGICAL               :: YN_ERROR
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitProfileDiags'

    ! =====================================================================
    ! InitProfileDiags starts here
    ! =====================================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Copy options
    ProfDiag%Options = ProfDiagOpt
    
    ! Get vertical dimension from netCDF dimension list
    CALL match_names_in_dimlist( nc_dimlist, nc_ndimlist, (/'lmx'/), 1, DIM_ID, DIMS, Error,&
                                 VARNAME='InitProfileDiags'                                 )
    ProfDiag%lmx = DIMS(1)
    
  END SUBROUTINE InitProfileDiags
  
  SUBROUTINE ArchiveProfileDiagnostics(Profile, ProfDiag, Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(ProfileType),  INTENT(IN)    :: Profile
    TYPE(ProfDiagType), INTENT(IN)    :: ProfDiag
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    
    ! =====================================================================
    !  ArchiveProfileDiagnostics starts here
    ! =====================================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

  END SUBROUTINE ArchiveProfileDiagnostics
  
  SUBROUTINE WriteProfileDiagnostics(ncid, nc_dimlist, nc_ndimlist,NC_XID, NC_YID,&
                                     Profile, ProfDiag, ACTION, DO_XY, Error      )
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,            INTENT(IN)    :: ncid
    TYPE(NCDimType),    INTENT(IN)    :: nc_dimlist(nc_ndimlist)
    INTEGER,            INTENT(IN)    :: nc_ndimlist
    INTEGER,            INTENT(IN)    :: NC_XID
    INTEGER,            INTENT(IN)    :: NC_YID
    TYPE(ProfileType),  INTENT(IN)    :: Profile
    TYPE(ProfDiagType), INTENT(IN)    :: ProfDiag
    INTEGER,            INTENT(IN)    :: ACTION
    LOGICAL,            INTENT(IN)    :: DO_XY
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                              :: VAR_ID
    INTEGER, DIMENSION(1)                :: size_1d, tmp_1d
    CHARACTER(LEN=maxChar), DIMENSION(1) :: dim_1d
    INTEGER                              :: N, G, G_Prx
    CHARACTER(LEN=maxChar)               :: VARNAME
    REAL(KIND=8)                         :: tmpprof(ProfDiag%lmx)
    REAL(KIND=8)                         :: XGas(1)
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'WriteProfileDiagnostics'

    ! =====================================================================
    ! WriteProfileDiagnostics starts here
    ! =====================================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! ---------------------------------------------------------------------
    ! Meteorological Variables
    ! ---------------------------------------------------------------------
    
    ! Match Dimension Names for Edge variables
    dim_1d(1) = 'lmx_e'
    
    ! Get dimension
    CALL match_names_in_dimlist( nc_dimlist, nc_ndimlist, dim_1d, &
                                 1, tmp_1d, size_1d, Error        )
    
    ! Write Pressure Layer Edge
    IF( ProfDiag%Options%PressureEdge ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%PressureEdge(:), 'PressureEdge',           &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error )
      
    ENDIF
    
    ! Write Temperature Layer Edge
    IF( ProfDiag%Options%TemperatureEdge ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%TemperatureEdge(:), 'TemperatureEdge',    &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
      
    ENDIF
    
    ! Write Altitude Layer Edge
    IF( ProfDiag%Options%AltitudeEdge ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%AltitudeEdge(:), 'AltitudeEdge',          &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
      
    ENDIF
    
    
    ! Match Dimension Names for Midpoint variables
    dim_1d(1) = 'lmx'
      
    ! Get dimension
    CALL match_names_in_dimlist( nc_dimlist, nc_ndimlist, dim_1d, &
                                 1, tmp_1d, size_1d, Error        )
    
    ! Write Pressure Layer Midpoint
    IF(ProfDiag%Options%PressureMid) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%PressureMid(:), 'PressureMid',            &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
      
    ENDIF
    
    ! Write Temperature Layer Midpoint
    IF( ProfDiag%Options%TemperatureMid ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%TemperatureMid(:), 'TemperatureMid',      &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID,Error )
      
    ENDIF
    
    ! Write Altitude Layer Midpoint
    IF( ProfDiag%Options%AltitudeMid ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%AltitudeMid(:), 'AltitudeMid',            &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
      
    ENDIF
    
    IF( ProfDiag%Options%AirMolecularWeight ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%AirMolecularWeight(:), 'AirMolecularWeight',  &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error  )
      
    ENDIF
    
    IF( ProfDiag%Options%Gravity ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%Gravity(:), 'Gravity',                   &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID,Error)
      
    ENDIF
    
    IF( ProfDiag%Options%RH ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%RH(:), 'RH',  &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
      
    ENDIF
    
    IF( ProfDiag%Options%AirPartialColumn ) THEN
      
      CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Met%AirPartialColumn(:), 'AirPartialColumn',   &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error )
      
    ENDIF
    
    ! ---------------------------------------------------------------------
    ! Gas Variables
    ! ---------------------------------------------------------------------
    
    IF( ProfDiag%Options%GasPartialColumn%DoDiag ) THEN
      
      DO N=1,ProfDiag%Options%GasPartialColumn%nSpc
        
        G       = ProfDiag%Options%GasPartialColumn%Idx(N)
        VARNAME = TRIM(ADJUSTL(ProfDiag%Options%GasPartialColumn%Name(N))) // '_GasPartialColumn'
        
        CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Gas%PartialColumn(:,G), TRIM(ADJUSTL(VARNAME)), &
                       NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error   )
        
      ENDDO
      
    ENDIF
    
!     IF( ProfDiag%Options%GasUncertainty%DoDiag ) THEN
!       
!       DO N=1,ProfDiag%Options%GasUncertainty%nSpc
!         
!         G       = ProfDiag%Options%GasUncertainty%Idx(N)
!         VARNAME = TRIM(ADJUSTL(ProfDiag%Options%GasUncertainty%Name(N))) // '_GasUncertainty'
!         
!         CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Gas%Uncertainty(:,G), TRIM(ADJUSTL(VARNAME)), &
!                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID        )
!         
!       ENDDO
!       
!     ENDIF
    
    IF( ProfDiag%Options%GasMixingRatio%DoDiag ) THEN
      
      DO N=1,ProfDiag%Options%GasMixingRatio%nSpc
        
        G       = ProfDiag%Options%GasMixingRatio%Idx(N)
        VARNAME = TRIM(ADJUSTL(ProfDiag%Options%GasMixingRatio%Name(N))) // '_GasMixingRatio'
        
        CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Gas%MixingRatio(:,G), TRIM(ADJUSTL(VARNAME)), &
                       NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
        
      ENDDO
      
    ENDIF
    
    IF( ProfDiag%Options%DryGasMixingRatio%DoDiag ) THEN

      DO N=1,ProfDiag%Options%DryGasMixingRatio%nSpc
        
        G       = ProfDiag%Options%DryGasMixingRatio%Idx(N)
        VARNAME = TRIM(ADJUSTL(ProfDiag%Options%DryGasMixingRatio%Name(N))) // '_GasDryAirMixingRatio'
        
        ! Set Profile if time to write
        IF(ACTION .EQ. 2) THEN
          tmpprof(:) = Profile%Gas%PartialColumn(:,G) & 
                     / Profile%Met%DryAirPartialColumn
        ELSE
          tmpprof(:) = 0.0d0
        ENDIF

        CALL nc_fld_1d( dim_1d, size_1d(1), tmpprof, TRIM(ADJUSTL(VARNAME)),              &
                       NCID, ACTION, DO_XY, nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
        
      ENDDO

    ENDIF



    ! ---------------------------------------------------------------------
    ! Aerosol Variables
    ! ---------------------------------------------------------------------
    IF( ProfDiag%Options%TotalAOD ) THEN
      
      ! Zero temporary profile array
      tmpprof(:) = 0.0d0

      DO G=1,Profile%Aer%nSpecies
        tmpprof(:) = tmpprof(:) + Profile%Aer%LayerOpticalDepth(:,G)
      ENDDO

      ! Write field
      CALL nc_fld_1d( dim_1d, size_1d(1), tmpprof, 'TotalAOD', &
                      NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID,Error)
    ENDIF
    
    IF( ProfDiag%Options%LayerOpticalDepth%DoDiag ) THEN
      
      DO N=1,ProfDiag%Options%LayerOpticalDepth%nSpc
        
        G       = ProfDiag%Options%LayerOpticalDepth%Idx(N)
        VARNAME = TRIM(ADJUSTL(ProfDiag%Options%LayerOpticalDepth%Name(N))) // '_LayerAOD'
        
        CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%LayerOpticalDepth(:,G), TRIM(ADJUSTL(VARNAME)), &
                       NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist,  NC_XID, NC_YID, Error)
        
      ENDDO
      
    ENDIF

    ! Profile Parameter derivatives
    IF( ProfDiag%Options%ProfileParDeriv%DoDiag ) THEN
      
      DO N=1,ProfDiag%Options%ProfileParDeriv%nSpc
        
        ! Get the index
        G = ProfDiag%Options%ProfileParDeriv%Idx(N)

        ! Column AOD derivative
        IF(Profile%Aer%TypeIndex(G) .GT. 1) THEN ! check for param. profile
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfileParDeriv%Name(N))) // '_ColumnAODDeriv'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%ColumnOptDepthDeriv(:,G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error          )
        ENDIF

        ! Additional derivatives for GDF profile type
        IF( Profile%Aer%TypeIndex(G) .EQ. 2 ) THEN
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfileParDeriv%Name(N))) // '_AltPeakDeriv'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltPeakDeriv(:,G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error   )
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfileParDeriv%Name(N))) // '_AltWidthDeriv'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltSigmaDeriv(:,G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error    )
        ENDIF

        ! Additional derivatives for EXP profile type
        IF(Profile%Aer%TypeIndex(G) .EQ. 3) THEN ! check for param. profile
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfileParDeriv%Name(N))) // '_AltExpDeriv'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltExpDeriv(:,G), TRIM(ADJUSTL(VARNAME)),&
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error )
        ENDIF
        
      ENDDO
      
    ENDIF

    ! The following fields are scalar
    dim_1d(1) = 'one'
      
    ! Get dimension
    CALL match_names_in_dimlist( nc_dimlist, nc_ndimlist, dim_1d, &
                                 1, tmp_1d, size_1d, Error        )

    IF(ProfDiag%Options%ProfilePar%DoDiag) THEN


      DO N=1,ProfDiag%Options%ProfilePar%nSpc
        
        ! Get the index
        G = ProfDiag%Options%ProfilePar%Idx(N)
        
        ! Variables common to all par. profiles
        IF(Profile%Aer%TypeIndex(G) .GT. 1) THEN
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfilePar%Name(N))) // '_ColumnAOD'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%ColumnOpticalDepth(G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfilePar%Name(N))) // '_AltMin'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltMin(G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfilePar%Name(N))) // '_AltMax'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltMax(G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
          
        ENDIF

        ! GDF parameters
        IF( Profile%Aer%TypeIndex(G) .EQ. 2 ) THEN

          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfilePar%Name(N))) // '_AltPeak'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltPeak(G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfilePar%Name(N))) // '_AltWidth'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltSigma(G), TRIM(ADJUSTL(VARNAME)), &
                        NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)

        ENDIF

        ! EXP parameters
        IF( Profile%Aer%TypeIndex(G) .EQ. 3 ) THEN
          VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProfilePar%Name(N))) // '_AltExp'
          CALL nc_fld_1d( dim_1d, size_1d(1), Profile%Aer%AltExp(G), TRIM(ADJUSTL(VARNAME)), &
               +         NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)
        ENDIF
      ENDDO
      
    ENDIF

    ! Proxy Gas Column Mixing Ratio
    IF( ProfDiag%Options%ProxyColumnMixingRatio%DoDiag ) THEN

      G_Prx = Profile%Gas%ProxyNormIdx
      DO N=1,ProfDiag%Options%ProxyColumnMixingRatio%nSpc
        
        G       = ProfDiag%Options%ProxyColumnMixingRatio%Idx(N)
        VARNAME = TRIM(ADJUSTL(ProfDiag%Options%ProxyColumnMixingRatio%Name(N))) // '_ProxyMixingRatio'
        
        IF(ACTION .EQ. 2) THEN
          XGas(1) = SUM(Profile%Gas%PartialColumn(:,G))     &
                  / SUM(Profile%Gas%PartialColumn(:,G_Prx)) &
                  * Profile%Gas%AprioriProxyNormMixingRatio
        ELSE
          XGas(1) = 0.0d0
        ENDIF

        CALL nc_fld_1d(dim_1d, size_1d(1), XGas, TRIM(ADJUSTL(VARNAME)), &
                       NCID,   ACTION,     DO_XY , nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error)

      ENDDO

    ENDIF

  END SUBROUTINE WriteProfileDiagnostics

END MODULE profile_module
