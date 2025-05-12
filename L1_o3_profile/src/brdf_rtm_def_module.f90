MODULE rtm_def

  USE parameters_module
!#if defined( USE_VLIDORT2p7 ) || defined(USE_VLIDORTpca)
  USE VLIDORT_PARS
!#elif defined( USE_VLIDORT2p8 )
!  USE VLIDORT_PARS_M
!#endif

  IMPLICIT NONE
  
  ! Common Radiative Transfer Settings to ensure consistency between models
  LOGICAL, PARAMETER :: VL_DO_FOCORR              = .TRUE.  ! Do overall Single scatter 
  LOGICAL, PARAMETER :: VL_DO_FULLRAD_MODE        = .TRUE.  ! Do full Stokes vector calculation? (Error when false)
  LOGICAL, PARAMETER :: VL_DO_FOCORR_NADIR        = .FALSE. ! Do nadir single scatter correction?
  LOGICAL, PARAMETER :: VL_DO_FOCORR_OUTGOING     = .TRUE.  ! Do outgoing single scatter correction?
  LOGICAL, PARAMETER :: VL_DO_SOLAR_SOURCES       = .TRUE.  ! Add sun
  LOGICAL, PARAMETER :: VL_DO_DO_CHAPMAN_FUNCTION = .TRUE.  ! Do internal Chapman function calculation?

  ! These flags are useful for single contiguous cloud layer in an otherwise Rayleigh atmos.
  LOGICAL, PARAMETER :: VL_DO_SOLUTION_SAVING      = .FALSE.
  LOGICAL, PARAMETER :: VL_DO_BVP_TELESCOPING      = .FALSE.

  ! Multiple geometry triplets (T for sat obs)
  LOGICAL, PARAMETER :: VL_DO_OBSERVATION_GEOMETRY = .FALSE. ! Do Observation Geometry?

  

  ! Initialization options
  TYPE JacOptType
    LOGICAL              :: DoTraceGas
    LOGICAL              :: DoTraceGasStokes
    LOGICAL              :: DoTraceGasNorm
    LOGICAL              :: DoSurfaceAmplitude
    LOGICAL              :: DoSurfaceAmplitudeStokes
    LOGICAL              :: DoSurfaceAmplitudeNorm
    LOGICAL              :: DoSurfaceParameter
    LOGICAL              :: DoSurfaceParameterStokes
    LOGICAL              :: DoSurfaceParameterNorm
    LOGICAL              :: DoTemperature
    LOGICAL              :: DoTemperatureStokes
    LOGICAL              :: DoTemperatureNorm
    LOGICAL              :: DoAerOD
    LOGICAL              :: DoAerODStokes
    LOGICAL              :: DoAerODNorm
    INTEGER              :: nAer
    LOGICAL, ALLOCATABLE :: DoAerODSpec(:)
    LOGICAL              :: DoAerSSA
    LOGICAL              :: DoAerSSAStokes
    LOGICAL              :: DoAerSSANorm
    LOGICAL, ALLOCATABLE :: DoAerSSASpec(:)
    LOGICAL              :: DoRayleigh
    LOGICAL              :: DoRayleighStokes
    LOGICAL              :: DoSurfPressureNorm
    LOGICAL              :: LinearizeHybridGrid
    LOGICAL              :: DoISRFPar
    LOGICAL              :: DoISRFParStokes
    LOGICAL              :: DoWvlShift
    LOGICAL              :: DoWvlShiftStokes

    ! Not yet implemented
    LOGICAL :: DoCldOD
    LOGICAL :: DoCldODStokes
    LOGICAL :: DoCldODNorm
    LOGICAL :: DoCldSSA
    LOGICAL :: DoCldSSAStokes
    LOGICAL :: DoCldSSANorm
    LOGICAL :: DoCldFraction
    LOGICAL :: DoCldFractionStokes
    LOGICAL :: DoCldFractionNorm

  ENDTYPE JacOptType
  
  ! Settings for VLIDORT
  TYPE VLOptType
    LOGICAL          :: do_full_stokes
    INTEGER          :: nstokes
    REAL(KIND=8)     :: ThermalEmissionCutOff
    INTEGER          :: ViewDirectionIndex
    LOGICAL          :: do_upwelling_radiation
    LOGICAL          :: do_downwelling_radiation
    INTEGER          :: nstreams
    INTEGER          :: nmoments
    LOGICAL          :: do_debug_write
  ENDTYPE VLOptType
  
  ! Settings VLIDORT PCA
  TYPE PCAVLOptType
    CHARACTER(LEN=maxChar) :: WinControlFile
    LOGICAL                :: do_full_stokes
    INTEGER                :: nstokes
    REAL(KIND=8)           :: ThermalEmissionCutOff
    INTEGER                :: ViewDirectionIndex
    LOGICAL                :: do_upwelling_radiation
    LOGICAL                :: do_downwelling_radiation
    INTEGER                :: nstreams
    INTEGER                :: nmoments
    LOGICAL                :: do_debug_write
  ENDTYPE PCAVLOptType
  
  ! Settings For First Order Model
  TYPE FrstOrdOptType
    LOGICAL                :: do_full_stokes
    INTEGER                :: nstokes
    REAL(KIND=8)           :: ThermalEmissionCutOff
    INTEGER                :: ViewDirectionIndex
    LOGICAL                :: do_upwelling_radiation
    LOGICAL                :: do_downwelling_radiation
    INTEGER                :: nmoments
  ENDTYPE FrstOrdOptType

  ! Settings for two stream model
  TYPE TwoStrOptType
    LOGICAL                :: do_full_stokes
    INTEGER                :: nstokes
    REAL(KIND=8)           :: ThermalEmissionCutOff
    INTEGER                :: ViewDirectionIndex
    LOGICAL                :: do_upwelling_radiation
    LOGICAL                :: do_downwelling_radiation
    INTEGER                :: nmoments
  ENDTYPE TwoStrOptType

  TYPE AllSimOptType
    TYPE(VLOptType)      :: VLIDORT
    TYPE(PCAVLOptType)   :: VLIDORT_PCA
    TYPE(FrstOrdOptType) :: FIRST_ORDER
    TYPE(TwoStrOptType)  :: TWO_STREAM
  ENDTYPE
  
  ! The overarching RTM Options Read by the file
  TYPE RTMOptType
    INTEGER                             :: nSimType
    CHARACTER(LEN=maxChar), ALLOCATABLE :: SimName(:)
    CHARACTER(LEN=maxChar), ALLOCATABLE :: RTMName(:)
    TYPE(AllSimOptType),    ALLOCATABLE :: SimOpt(:)
    TYPE(JacOptType)                    :: Jacobian
    LOGICAL                             :: DoScatteringWeights
    LOGICAL                             :: DoFluxCalculation
  ENDTYPE RTMOptType
  
  ! For A specific Window
  TYPE WinRTMSettingsType
    CHARACTER(LEN=maxChar) :: RTMName
    TYPE(ALLSimOptType)    :: SimOpt
    TYPE(JacOptType)       :: Jacobian
    LOGICAL                :: DoScatteringWeights
    LOGICAL                :: DoFluxCalculation
  ENDTYPE WinRTMSettingsType
  
  TYPE SurfJacType
    INTEGER                                :: idx ! VLIDORT Output Index
    LOGICAL                                :: DoStokes
    INTEGER                                :: nJacobian
    CHARACTER(LEN=maxChar), ALLOCATABLE    :: Name(:)
    REAL(KIND=8), ALLOCATABLE              :: Output(:,:,:)
    REAL(KIND=8), ALLOCATABLE              :: OutputISRF(:,:,:)
    REAL(KIND=8), ALLOCATABLE              :: pres(:,:,:)
  ENDTYPE SurfJacType
  
  TYPE ProfJacType
    INTEGER                   :: idx ! VLIDORT Output Index
    LOGICAL                   :: DoStokes
    CHARACTER(LEN=maxChar)    :: Name
    REAL(KIND=8), ALLOCATABLE :: Output(:,:,:)
    REAL(KIND=8), ALLOCATABLE :: OutputISRF(:,:,:)
  ENDTYPE ProfJacType
  
  TYPE RTMOutputType
    INTEGER                        :: BandIndex
    INTEGER                        :: nL1_Wvl
    INTEGER                        :: L1_iw0
    INTEGER                        :: L1_iwf
    INTEGER                        :: nStokes
    REAL(KIND=8),      ALLOCATABLE :: Irradiance(:)
    REAL(KIND=8),      ALLOCATABLE :: IrradianceISRF(:)
    REAL(KIND=8),      ALLOCATABLE :: RadianceOffset(:,:)
    REAL(KIND=8),      ALLOCATABLE :: RadianceScale(:,:)
    REAL(KIND=8),      ALLOCATABLE :: ResidualEOF(:,:,:)
    REAL(KIND=8),      ALLOCATABLE :: Radiance(:,:)
    REAL(KIND=8),      ALLOCATABLE :: RadianceISRF(:,:)
    REAL(KIND=8),      ALLOCATABLE :: RadiantFlux(:,:)
    REAL(KIND=8),      ALLOCATABLE :: RadiantFluxISRF(:,:)
    REAL(KIND=8),      ALLOCATABLE :: DirectFlux(:,:)
    REAL(KIND=8),      ALLOCATABLE :: DirectFluxISRF(:,:)
    REAL(KIND=8),      ALLOCATABLE :: ScatteringWeights(:,:)
    REAL(KIND=8),      ALLOCATABLE :: ScatteringWeightsISRF(:,:)
    REAL(KIND=8),      ALLOCATABLE :: RamanProbability(:)
    REAL(KIND=8),      ALLOCATABLE :: RamanSourceNorm(:)
    TYPE(JacOptType)               :: JacOpt ! Store A copy of the output jacobian options
    TYPE(ProfJacType)              :: TraceGasJacobian
    TYPE(ProfJacType)              :: RayleighJacobian
    TYPE(ProfJacType), ALLOCATABLE :: AerODJacobian(:)
    TYPE(ProfJacType), ALLOCATABLE :: AerSSAJacobian(:)
    TYPE(SurfJacType)              :: BRDFJacobian
    TYPE(SurfJacType)              :: ISRFParJacobian
    TYPE(SurfJacType)              :: WvlShiftJacobian
  ENDTYPE RTMOutputType
  
  ! --------------------------------------------
  ! Diagnostics 
  ! --------------------------------------------

  TYPE RTMDiagOptType
    INTEGER          :: nStokes
    LOGICAL          :: Wavelength
    LOGICAL          :: Irradiance
    TYPE(DiagSpcOpt) :: Radiance
    TYPE(DiagSpcOpt) :: RadiantFlux
    TYPE(DiagSpcOpt) :: DirectFlux
    TYPE(DiagSpcOpt) :: TraceGasJacobian
    LOGICAL          :: ScatteringWeights
    TYPE(DiagSpcOpt) :: AirMassFactor
    TYPE(DiagSpcOpt) :: TemperatureJacobian
    TYPE(DiagSpcOpt) :: TemperatureShiftJacobian
    TYPE(DiagSpcOpt) :: RayleighJacobian
    TYPE(DiagSpcOpt) :: AerODJacobian
    TYPE(DiagSpcOpt) :: AerODProfileParJacobian
    TYPE(DiagSpcOpt) :: AerSSAJacobian
    TYPE(DiagSpcOpt) :: AerSSAProfileParJacobian
    TYPE(DiagSpcOpt) :: BRDFAmplitudeJacobian
    TYPE(DiagSpcOpt) :: BRDFParameterJacobian
    TYPE(DiagSpcOpt) :: ISRFParJacobian
    TYPE(DiagSpcOpt) :: WavelengthShiftJacobian

    ! Not yet implemented
    TYPE(DiagSpcOpt) :: CldODJacobian
    TYPE(DiagSpcOpt) :: CldODProfileParJacobian
    TYPE(DiagSpcOpt) :: CldSSAJacobian
    TYPE(DiagSpcOpt) :: CldSSAProfileParJacobian
    TYPE(DiagSpcOpt) :: CldFractionJacobian
  ENDTYPE RTMDiagOptType

  TYPE RTMDiagType
    TYPE(RTMDiagOptType)      :: Options
    INTEGER                   :: nBand
    INTEGER, ALLOCATABLE      :: wmx(:)
  ENDTYPE RTMDiagType
  
  
  ! Some convenience routines for output allocation/deallocation
  
  CONTAINS
  
  SUBROUTINE DeallocateRTMOutput(RTMOutput)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(RTMOutputType), INTENT(INOUT) :: RTMOutput
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: N
    
    ! ===============================================================
    ! DeallocateRTMOutput Starts here
    ! ===============================================================
    
    ! Radiance/flux variables
    IF( ALLOCATED(RTMOutput%Radiance) ) DEALLOCATE(RTMOutput%Radiance)
    IF( ALLOCATED(RTMOutput%RadianceISRF) ) DEALLOCATE(RTMOutput%RadianceISRF)
    IF( ALLOCATED(RTMOutput%RadiantFlux) ) DEALLOCATE(RTMOutput%RadiantFlux)
    IF( ALLOCATED(RTMOutput%RadiantFluxISRF) ) DEALLOCATE(RTMOutput%RadiantFluxISRF)
    IF( ALLOCATED(RTMOutput%DirectFlux) ) DEALLOCATE(RTMOutput%DirectFlux)
    IF( ALLOCATED(RTMOutput%DirectFluxISRF) ) DEALLOCATE(RTMOutput%DirectFluxISRF)
    IF( ALLOCATED(RTMOutput%RadianceOffset)) DEALLOCATE(RTMOutput%RadianceOffset)
    IF( ALLOCATED(RTMOutput%ResidualEOF)) DEALLOCATE(RTMOutput%ResidualEOF)
    IF( ALLOCATED(RTMOutput%RadianceScale)) DEALLOCATE(RTMOutput%RadianceScale)
    
    ! Trace Gas Jacobian
    IF( ALLOCATED(RTMOutput%TraceGasJacobian%Output) ) &
        DEALLOCATE(RTMOutput%TraceGasJacobian%Output)
    IF( ALLOCATED(RTMOutput%TraceGasJacobian%OutputISRF) ) &
        DEALLOCATE(RTMOutput%TraceGasJacobian%OutputISRF)
    IF( ALLOCATED(RTMOutput%ScatteringWeights) ) &
          DEALLOCATE(RTMOutput%ScatteringWeights)
    IF( ALLOCATED(RTMOutput%ScatteringWeightsISRF) ) &
          DEALLOCATE(RTMOutput%ScatteringWeightsISRF)
    
    ! Surface Pressure
    IF( ALLOCATED(RTMOutput%RayleighJacobian%Output) ) &
        DEALLOCATE(RTMOutput%RayleighJacobian%Output)
    IF( ALLOCATED(RTMOutput%RayleighJacobian%OutputISRF) ) &
        DEALLOCATE(RTMOutput%RayleighJacobian%OutputISRF)
    
    ! Aerosol OD+SSA
    IF(ALLOCATED(RTMOutput%AerODJacobian)) DEALLOCATE(RTMOutput%AerODJacobian)
    IF(ALLOCATED(RTMOutput%AerSSAJacobian)) DEALLOCATE(RTMOutput%AerSSAJacobian)
    
    ! Surface BRDF
    IF(ALLOCATED(RTMOutput%BRDFJacobian%Output)) &
        DEALLOCATE(RTMOutput%BRDFJacobian%Output)
    IF(ALLOCATED(RTMOutput%BRDFJacobian%OutputISRF)) &
        DEALLOCATE(RTMOutput%BRDFJacobian%OutputISRF)
    
  END SUBROUTINE DeallocateRTMOutput
  
  
  SUBROUTINE AllocateRTMOutput(RTMOutput,Opt,nRTM_Wvl,nL1_Wvl,nStokes,nLevels,nSurfaceJac)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(RTMOutputType),      INTENT(INOUT) :: RTMOutput
    TYPE(WinRTMSettingsType), INTENT(IN)    :: Opt
    INTEGER,                  INTENT(IN)    :: nRTM_Wvl, nL1_Wvl
    INTEGER,                  INTENT(IN)    :: nStokes
    INTEGER,                  INTENT(IN)    :: nLevels
    INTEGER,                  INTENT(IN)    :: nSurfaceJac
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: N
    
    ! ===============================================================
    ! AllocateRTMOutput Starts here
    ! ===============================================================
    
    ! Deallocate arrays
    CALL DeallocateRTMOutput(RTMOutput)
    
    ! ------------------------
    ! Radiance/ flux variables
    ! ------------------------
    ALLOCATE(RTMOutput%Radiance(nRTM_Wvl,nStokes))         ; RTMOutput%Radiance(:,:) = 0.0d0
    ALLOCATE(RTMOutput%RadianceISRF(nL1_Wvl,nStokes))      ; RTMOutput%RadianceISRF(:,:) = 0.0d0
    ALLOCATE(RTMOutput%RadiantFlux(nRTM_Wvl,nStokes))      ; RTMOutput%RadiantFlux(:,:) = 0.0d0
    ALLOCATE(RTMOutput%RadiantFluxISRF(nL1_Wvl,nStokes))   ; RTMOutput%RadiantFluxISRF(:,:) = 0.0d0
    ALLOCATE(RTMOutput%DirectFlux(nRTM_Wvl,nStokes))       ; RTMOutput%DirectFlux(:,:) = 0.0d0
    ALLOCATE(RTMOutput%DirectFluxISRF(nL1_Wvl,nStokes))    ; RTMOutput%DirectFluxISRF(:,:) = 0.0d0
    ALLOCATE(RTMOutput%RadianceOffset(nL1_Wvl,nStokes))   ; RTMOutput%RadianceOffset(:,:) = 0.0d0
    ALLOCATE(RTMOutput%RadianceScale(nL1_Wvl,nStokes))    ; RTMOutput%RadianceScale(:,:) = 1.0d0
    
    ! ------------------
    ! Trace Gas Jacobian
    ! ------------------
    RTMOutput%TraceGasJacobian%Idx      = 0
    RTMOutput%TraceGasJacobian%DoStokes = Opt%Jacobian%DoTraceGasStokes
    IF(Opt%Jacobian%DoTraceGas) THEN
      ALLOCATE(RTMOutput%TraceGasJacobian%Output(nRTM_Wvl,nLevels,nStokes))
      RTMOutput%TraceGasJacobian%Output(:,:,:) = 0.0d0
      ALLOCATE(RTMOutput%TraceGasJacobian%OutputISRF(nL1_Wvl,nLevels,nStokes))
      RTMOutput%TraceGasJacobian%OutputISRF(:,:,:) = 0.0d0
      IF(Opt%DoScatteringWeights) THEN
        ALLOCATE(RTMOutput%ScatteringWeights(nRTM_Wvl,nLevels))
        RTMOutput%ScatteringWeights(:,:) = 0.0d0
        ALLOCATE(RTMOutput%ScatteringWeightsISRF(nL1_Wvl,nLevels))
        RTMOutput%ScatteringWeightsISRF(:,:) = 0.0d0
      ENDIF
    ENDIF
    
    ! -----------------
    ! Pressure Jacobian
    ! -----------------
    RTMOutput%RayleighJacobian%Idx      = 0
    RTMOutput%RayleighJacobian%DoStokes = Opt%Jacobian%DoRayleighStokes
    IF(Opt%Jacobian%DoRayleigh) THEN
      ALLOCATE(RTMOutput%RayleighJacobian%Output(nRTM_Wvl,nLevels,nStokes))
      RTMOutput%RayleighJacobian%Output(:,:,:) = 0.0d0
      ALLOCATE(RTMOutput%RayleighJacobian%OutputISRF(nL1_Wvl,nLevels,nStokes))
      RTMOutput%RayleighJacobian%OutputISRF(:,:,:) = 0.0d0
    ENDIF
    
    ! ---------------------
    ! Aerosol Optical Depth
    ! ---------------------
    
    ! Allocate # of species
    ALLOCATE(RTMOutput%AerODJacobian(Opt%Jacobian%nAer))
    
    ! Initialize Aerosol indices
    DO N=1,Opt%Jacobian%nAer
      RTMOutput%AerODJacobian(N)%Idx       = 0
      RTMOutput%AerODJacobian(N)%DoStokes  = Opt%Jacobian%DoAerODStokes
    ENDDO
    
    ! Allocate arrays
    IF(Opt%Jacobian%DoAerOD) THEN
      DO N=1,Opt%Jacobian%nAer
        IF(Opt%Jacobian%DoAerODSpec(N)) THEN
          ALLOCATE(RTMOutput%AerODJacobian(N)%Output(nRTM_Wvl,nLevels,nStokes))
          RTMOutput%AerODJacobian(N)%Output(:,:,:) = 0.0d0
          ALLOCATE(RTMOutput%AerODJacobian(N)%OutputISRF(nL1_Wvl,nLevels,nStokes))
          RTMOutput%AerODJacobian(N)%OutputISRF(:,:,:) = 0.0d0
        ENDIF
      ENDDO
    ENDIF
    
    ! --------------------------------
    ! Aerosol Single Scattering Albedo
    ! --------------------------------
    
    ! Allocate # of species
    ALLOCATE(RTMOutput%AerSSAJacobian(Opt%Jacobian%nAer))

    ! Initialize Aerosol indices
    DO N=1,Opt%Jacobian%nAer
      RTMOutput%AerSSAJacobian(N)%Idx       = 0 
      RTMOutput%AerSSAJacobian(N)%DoStokes  = Opt%Jacobian%DoAerSSAStokes
    ENDDO
    
    ! Allocate Arrays
    IF( Opt%Jacobian%DoAerSSA ) THEN
      DO N=1,Opt%Jacobian%nAer
        IF(Opt%Jacobian%DoAerSSASpec(N)) THEN
          ALLOCATE(RTMOutput%AerSSAJacobian(N)%Output(nRTM_Wvl,nLevels,nStokes))
          RTMOutput%AerSSAJacobian(N)%Output(:,:,:) = 0.0d0
          ALLOCATE(RTMOutput%AerSSAJacobian(N)%OutputISRF(nL1_Wvl,nLevels,nStokes))
          RTMOutput%AerSSAJacobian(N)%OutputISRF(:,:,:) = 0.0d0
        ENDIF
      ENDDO
    ENDIF
    
    ! ----------------------
    ! Surface Linearizations
    ! ----------------------
    
    ! BRDF
    RTMOutput%BRDFJacobian%Idx           = 0
    RTMOutput%BRDFJacobian%DoStokes      = Opt%Jacobian%DoSurfaceAmplitudeStokes .OR. &
                                            Opt%Jacobian%DoSurfaceParameterStokes
    RTMOutput%BRDFJacobian%nJacobian     = 0
    
    ! Allocate arrays
    IF(nSurfaceJac .GT. 0) THEN
      
      ! Save number of jacobians
      RTMOutput%BRDFJacobian%nJacobian = nSurfaceJac

      ! Allocate output arrays
      ALLOCATE(RTMOutput%BRDFJacobian%Output(nRTM_Wvl,nStokes,nSurfaceJac))
      RTMOutput%BRDFJacobian%Output(:,:,:) = 0.0d0
      ALLOCATE(RTMOutput%BRDFJacobian%OutputISRF(nL1_Wvl,nStokes,nSurfaceJac))
      RTMOutput%BRDFJacobian%Output(:,:,:) = 0.0d0
      
    ENDIF
    
  END SUBROUTINE AllocateRTMOutput
  
  SUBROUTINE SetJacobianIndices(JacOpt,nAer,AerName,RTMOutput,nTotalWFS,VL_JacName)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(JacOptType),         INTENT(IN)    :: JacOpt
    INTEGER,                  INTENT(IN)    :: nAer          ! OptProp%AerSca%Options%nSpecies
    CHARACTER(LEN=maxChar),   INTENT(IN)    :: AerName(nAer) ! OptProp%AerSca%Options%Name(N)
    TYPE(RTMOutputType),      INTENT(INOUT) :: RTMOutput
    INTEGER,                  INTENT(INOUT) :: nTotalWFS
    CHARACTER(LEN=*),OPTIONAL,INTENT(INOUT) :: VL_JacName(MAX_ATMOSWFS)

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: N

    ! ===============================================================
    ! SetJacobianIndices Starts here
    ! ===============================================================
    
    ! Initialize Jacobian Counts
    nTotalWFS = 0

    ! Store the Jacobian Options
    RTMOutput%JacOpt = JacOpt!Window%RTMSettings%Jacobian

    ! -----------------------------------------------------------------
    ! Trace gas
    ! -----------------------------------------------------------------
    
    ! Set Trace Gas Jacobians
    IF( JacOpt%DoTraceGas ) THEN
      
      ! Increment Profile Jacobian Count
      nTotalWFS = nTotalWFS + 1
      
      ! Store Index 
      RTMOutput%TraceGasJacobian%Idx = nTotalWFS
      
      ! Name the jacobian
      IF(PRESENT(VL_JacName)) &
        VL_JacName(nTotalWFS) = '-Trace Gas Volume Mixing Ratio-'
      
    ENDIF
    
    ! -----------------------------------------------------------------
    ! Surface pressure
    ! -----------------------------------------------------------------
    
    IF( JacOpt%DoRayleigh ) THEN

      ! Increment Profile Jacobian Count
      nTotalWFS = nTotalWFS + 1
      
      ! Store Index 
      RTMOutput%RayleighJacobian%Idx = nTotalWFS
      
      ! Name the jacobian
      IF(PRESENT(VL_JacName)) &
        VL_JacName(nTotalWFS) = '-Rayleigh Scatter-'
      
    ENDIF

    ! -----------------------------------------------------------------
    ! Aerosol Optical Depth
    ! -----------------------------------------------------------------
    
    IF( JacOpt%DoAerOD ) THEN

      ! jacobians organized by aer profile index
      DO N=1,nAer
        IF(JacOpt%DoAerODSpec(N)) THEN

          ! Increment Profile Jacobian Count
          nTotalWFS = nTotalWFS + 1
          
          ! Store index
          RTMOutput%AerODJacobian(N)%Idx = nTotalWFS

          ! Name Jacobian
          IF(PRESENT(VL_JacName)) &
            VL_JacName(nTotalWFS) = 'AerOD_'//TRIM(ADJUSTL(AerName(N)))
          
        ENDIF
      ENDDO

    ENDIF

    ! -----------------------------------------------------------------
    ! Aerosol SSA
    ! -----------------------------------------------------------------
    IF( JacOpt%DoAerSSA ) THEN
      
      ! jacobians organized by aer profile index
      DO N=1,nAer
        IF(JacOpt%DoAerSSASpec(N)) THEN
        
          ! Increment Profile Jacobian Count
          nTotalWFS = nTotalWFS + 1
          
          ! Store index
          RTMOutput%AerSSAJacobian(N)%Idx = nTotalWFS

          ! Name Jacobian
          IF(PRESENT(VL_JacName)) &
            VL_JacName(nTotalWFS) ='AerSSA_'//TRIM(ADJUSTL(AerName(N)))
          
        
        ENDIF
      ENDDO

    ENDIF


  END SUBROUTINE SetJacobianIndices

END MODULE rtm_def
