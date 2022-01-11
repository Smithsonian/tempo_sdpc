MODULE surface_module
  
  
!#if defined( USE_VLIDORT2p7 ) || defined( USE_VLIDORTpca )
  ! VLIDORT
  USE VLIDORT_PARS
  USE VLIDORT_IO_DEFS
  USE VLIDORT_LININPUTS_DEF
  USE VLIDORT_LINSUP_INOUT_DEF
  USE VLIDORT_LINOUTPUTS_DEF
  USE VLIDORT_OUTPUTS_DEF
  
  ! VLIDORT SUPPLEMENT
  USE VBRDF_SUP_INPUTS_DEF
  USE VBRDF_LINSUP_INPUTS_DEF
  USE VBRDF_SUP_OUTPUTS_DEF
  USE VBRDF_LINSUP_OUTPUTS_DEF
!#elif defined( USE_VLIDORT2p8 )
!  ! VLIDORT
!  USE VLIDORT_PARS_M
!  USE VLIDORT_IO_DEFS_M
!  USE VLIDORT_LININPUTS_DEF_M
!  USE VLIDORT_LINSUP_INOUT_DEF_M
!  USE VLIDORT_LINOUTPUTS_DEF_M
!  USE VLIDORT_OUTPUTS_DEF_M
  
  ! VLIDORT SUPPLEMENT
!  USE VBRDF_SUP_INPUTS_DEF_M
!  USE VBRDF_LINSUP_INPUTS_DEF_M
!  USE VBRDF_SUP_OUTPUTS_DEF_M
!  USE VBRDF_LINSUP_OUTPUTS_DEF_M
!#endif
  
  USE parameters_module
  USE error_module,         ONLY : ErrorType, CheckError
  USE window_module,        ONLY : WinType,GetWindowNumStokes,GetWindowNumStreams
  USE level1_def,           ONLY : GeolocationType, AuxSurfType, L1Type
  USE profile_module,       ONLY : SurfProfType
  USE spatial_fa_module,    ONLY : SpatialFAType, InitSpatialFA, SampleSpatialFA
  USE netcdf_module,        ONLY : CheckNetCDFErrorStatus
  USE chebyshev_module,     ONLY : EvaluateChebyshev
  USE interpolation_module, ONLY : BSPLINE_EdgeFill
  USE sif_module,           ONLY : SIFType, SIFOptType

  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  
  ! Surface Options
  TYPE SurfOptType
    CHARACTER(LEN=maxChar) :: RootDataDir
    INTEGER                :: OptionIndex
    REAL(KIND=8)           :: Option1_FixedAlbedo
    CHARACTER(LEN=maxChar) :: Option2_Infile
    CHARACTER(LEN=maxChar) :: Option3_Infile
    LOGICAL                :: Option3_UseConstantWvl
    REAL(KIND=8)           :: Option3_ConstWvl
    CHARACTER(LEN=maxChar) :: Option4_Infile
    CHARACTER(LEN=maxChar) :: Option4_ClimDir
    LOGICAL                :: Option4_DoIsotropic
    INTEGER                :: Option4_WhichAlbedo
    LOGICAL                :: Option4_DoOceanGlint
    CHARACTER(LEN=maxChar) :: Option5_KernName(3)
    INTEGER                :: Option5_KernIdx(3)
    REAL(KIND=8)           :: Option5_KernAmp(3)
    INTEGER                :: Option5_nKernPar(3)
    REAL(KIND=8)           :: Option5_KernPar(5,3)
    CHARACTER(LEN=maxChar) :: Option6_Infile
    LOGICAL                :: DoPlantFluorescence
    TYPE(SIFOptType)       :: SIF
    LOGICAL                :: DoAmplitudeLinearization
    LOGICAL                :: DoParameterLinearization
    INTEGER                :: EmissivityOptIndex
    CHARACTER(LEN=maxChar) :: EmissivityInfile
    REAL(KIND=8)           :: ConstantEmissivity
  ENDTYPE SurfOptType
  
  ! BRDF arrays To be passed to VLIDORT
  TYPE VLBRDFInpType
    REAL(KIND=8), ALLOCATABLE :: BRDF_F_0(:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: BRDF_F(:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: USER_BRDF_F_0(:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: USER_BRDF_F(:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: DBOUNCE_BRDFUNC(:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: EMISSIVITY(:,:)
    REAL(KIND=8), ALLOCATABLE :: USER_EMISSIVITY(:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_BRDF_F_0(:,:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_BRDF_F(:,:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_USER_BRDF_F_0(:,:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_USER_BRDF_F(:,:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_DBOUNCE_BRDFUNC(:,:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_USER_EMISSIVITY(:,:,:)
    REAL(KIND=8), ALLOCATABLE :: LS_EMISSIVITY(:,:,:)
  ENDTYPE VLBRDFInpType
  
  ! Case 1: Fixed lambertian Albedo
  
  ! Case 2: Fixed Lambertian Spectrum
  TYPE Surf2Type
    INTEGER                   :: nWvl
    REAL(KIND=8), ALLOCATABLE :: Wavelength(:)
    REAL(KIND=8), ALLOCATABLE :: Albedo(:)
    REAL(KIND=8), ALLOCATABLE :: AlbedoSP(:)
  ENDTYPE Surf2Type
  
  ! Case 3: LER Climatology
  TYPE Surf3Type
    INTEGER                   :: ncid
    INTEGER                   :: imx
    INTEGER                   :: jmx
    INTEGER                   :: wmx
    REAL(KIND=8), ALLOCATABLE :: Wavelength(:)
    REAL(KIND=8), ALLOCATABLE :: Longitude(:)
    REAL(KIND=8), ALLOCATABLE :: Latitude(:)
    INTEGER                   :: i
    INTEGER                   :: j
    INTEGER                   :: Month
    REAL(KIND=8), ALLOCATABLE :: Albedo(:)
    REAL(KIND=8), ALLOCATABLE :: AlbedoSP(:)
    REAL(KIND=8)              :: ConstWvl(1)
  ENDTYPE Surf3Type
  
  ! Case 4: SCIAMACHY-USGS FA Model

  ! Case 5: Fixed BRDF Kernels
  
  ! Case 6: BRDF Climatology
  TYPE Surf6Type
    INTEGER                   :: ncid
    INTEGER                   :: imx
    INTEGER                   :: jmx
    INTEGER                   :: wmx
    INTEGER                   :: pmx ! Max pars
    REAL(KIND=8), ALLOCATABLE :: Wavelength(:)
    REAL(KIND=8), ALLOCATABLE :: Longitude(:)
    REAL(KIND=8), ALLOCATABLE :: Latitude(:)
    INTEGER                   :: i
    INTEGER                   :: j
    INTEGER,      ALLOCATABLE :: KernIdx(:)
    REAL(KIND=8), ALLOCATABLE :: Amplitude(:,:)
    REAL(KIND=8), ALLOCATABLE :: AmplitudeSP(:,:)
    REAL(KIND=8), ALLOCATABLE :: Parameters(:,:,:)
    REAL(KIND=8), ALLOCATABLE :: ParametersSP(:,:,:)
  ENDTYPE Surf6Type
  
  TYPE SurfaceEmissType
    INTEGER                   :: OptionIndex
    CHARACTER(LEN=maxChar)    :: Filename
    REAL(KIND=8)              :: ConstantEmissivity
    REAL(KIND=8), ALLOCATABLE :: Value(:) ! Emissivity at RTM wavelengths
    
    ! Kernels
  ENDTYPE SurfaceEmissType
  
  TYPE SurfaceType
    LOGICAL                             :: UseL2Reflectance
    INTEGER                             :: OptionIndex
    INTEGER                             :: nstokes
    INTEGER                             :: nstreams
    INTEGER                             :: nkern ! Number of BRDF Kernels
    INTEGER,                ALLOCATABLE :: npar(:) ! Number of parameters for each BRDF kernel
    CHARACTER(LEN=maxChar), ALLOCATABLE :: kern_name(:) ! Names of BRDF Kernels
    INTEGER,                ALLOCATABLE :: kern_idx(:) ! VLIDORT kernel indices
    REAL(KIND=8),           ALLOCATABLE :: kern_amp(:,:) ! Kernel Amplitudes
    REAL(KIND=8),           ALLOCATABLE :: kern_par(:,:,:) ! Kernel Parameters
    LOGICAL                             :: DoLambertian ! Kernel model is purely lambertian
    REAL(KIND=8),           ALLOCATABLE :: WavelengthKernelAmplitude(:)
    REAL(KIND=8),           ALLOCATABLE :: WavelengthKernelParameters(:,:)
    LOGICAL                             :: fixed_par ! Kernel model has non-wavelength dependent parameters
    LOGICAL                             :: fixed_amp ! Kernel model has non-wavelength dependent amplitudes
    TYPE(Surf2Type)                     :: Option2
    TYPE(Surf3Type)                     :: Option3
    TYPE(SpatialFAType)                 :: Option4
    TYPE(Surf6Type)                     :: Option6
    TYPE(VLBRDFInpType)                 :: VLBRDF
    LOGICAL                             :: DoAmplitudeLinearization
    LOGICAL                             :: DoParameterLinearization
    INTEGER                             :: nJac
    CHARACTER(LEN=maxChar), ALLOCATABLE :: JacName(:)
    INTEGER,                ALLOCATABLE :: JacKernIdx(:)
    INTEGER,                ALLOCATABLE :: JacParIdx(:)
    TYPE(SurfaceEmissType)              :: Emissivity
    LOGICAL                             :: DoPlantFluorescence
    TYPE(SIFType)                       :: SIF

    ! New Kernel Model to avoid multiple VBRDF supplement calls
    LOGICAL                             :: UseFixedParKernels
    TYPE(VLBRDFInpType), ALLOCATABLE    :: VBRDF_Kernels(:)

    ! Functions
    CONTAINS
      PROCEDURE :: Init                           => InitSurface
      PROCEDURE :: SamplePixel                    => SampleSurfaceProperties
      PROCEDURE :: SetOpticalProperties_VL        => SetSurfaceOpticalProperties
      PROCEDURE :: SetOpticalProperties_2S        => SetSurfaceOpticalProperties_2S
      PROCEDURE :: GetDirectBounceReflectance     => GetDirectBounceReflectance
      PROCEDURE :: SetLinearization               => SetSurfaceLinearization

  ENDTYPE SurfaceType
  
  ! These are used locally by the surface module
  ! Include here so we don't allocate/deallocate each call
  TYPE(VBRDF_Sup_Inputs),                PRIVATE  :: VBRDF_Sup_In_surfmod
  TYPE(VBRDF_LinSup_Inputs),             PRIVATE  :: VBRDF_LinSup_In_surfmod
  TYPE(VBRDF_Input_Exception_Handling),  PRIVATE  :: VBRDF_Sup_InputStatus_surfmod 
  TYPE(VBRDF_Sup_Outputs),               PRIVATE  :: VBRDF_SupOut_surfmod
  TYPE(VBRDF_LinSup_Outputs),            PRIVATE  :: VBRDF_LinSupOut_surfmod
  TYPE(VBRDF_Output_Exception_Handling), PRIVATE  :: VBRDF_Sup_OutputStatus_surfmod
  


  ! While Testing new BRDF calculation option to avoid recomputation of kernels
  LOGICAL, PRIVATE, PARAMETER :: UseNewFixedParBRDF = .TRUE.

  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'surface_module'
  PRIVATE :: ModuleName

  CONTAINS
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: InitSurface
  ! 
  ! DESCRIPTION: Initialization of surface reflectance
  
  SUBROUTINE InitSurface(self, SurfOpt, Window, L2Surface, Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    CLASS(SurfaceType),INTENT(INOUT)   :: self
    TYPE(SurfOptType), INTENT(IN)      :: SurfOpt
    TYPE(WinType),     INTENT(IN)      :: Window
    TYPE(AuxSurfType), INTENT(IN)      :: L2Surface
    TYPE(ErrorType),   INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: n
    
    ! =====================================================================
    ! InitSurface starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Check if we are using the L2 Product
    self%UseL2Reflectance = L2Surface%OverwriteSurface

    ! Save Option Index
    self%OptionIndex = SurfOpt%OptionIndex
    
    ! Get stokes/streams from window
    self%nstokes  = GetWindowNumStokes(Window,Error)
    self%nstreams = GetWindowNumStreams(Window,Error)
    
    ! Save linearization options
    self%DoAmplitudeLinearization = SurfOpt%DoAmplitudeLinearization
    self%DoParameterLinearization = SurfOpt%DoParameterLinearization
    
    ! Initialize the selected surface model
    IF(self%UseL2Reflectance) THEN
      CALL InitSurfaceOption0_L2(L2Surface,Window,self,Error)
    ELSEIF(SurfOpt%OptionIndex .EQ. 1) THEN
      CALL InitSurfaceOption1_FixedAlb(SurfOpt,Window,self,Error)
    ELSEIF(SurfOpt%OptionIndex .EQ. 2) THEN
      CALL InitSurfaceOption2_AlbSpec(SurfOpt,self,Window,Error)
    ELSEIF(SurfOpt%OptionIndex .EQ. 3) THEN
      CALL InitSurfaceOption3_LERClim(SurfOpt,self,Window,Error)
    ELSEIF(SurfOpt%OptionIndex .EQ. 4) THEN
      CALL InitSurfaceOption4_SciaFA(SurfOpt,self,Window,Error)
    ELSEIF(SurfOpt%OptionIndex .EQ. 5) THEN
      CALL InitSurfaceOption5_FixedBRDF(SurfOpt,self,Window,Error)
    ELSEIF(SurfOpt%OptionIndex .EQ. 6) THEN
      CALL InitSurfaceOption6_ClimBRDF(SurfOpt,self,Window,Error)
    ELSE
      STOP 'Surface Option must range from 1..6'
    ENDIF
    
    ! Add option 0 for kirchoffs law?
    
    ! Emissivity Options
    self%Emissivity%OptionIndex        = SurfOpt%EmissivityOptIndex
    self%Emissivity%Filename           = SurfOpt%EmissivityInfile
    self%Emissivity%ConstantEmissivity = SurfOpt%ConstantEmissivity
    
    ! We need to skip the linearizations to prepare diagnostic outputs
    CALL self%SetLinearization(SurfOpt%DoAmplitudeLinearization, &
                               SurfOpt%DoAmplitudeLinearization, &
                               Error                             )
    
    
    ! Save SIF
    self%DoPlantFluorescence = SurfOpt%DoPlantFluorescence
    
    ! Initialize solar induced fluorescence
    IF(self%DoPlantFluorescence) THEN 
      CALL self%SIF%Initialize(SurfOpt%SIF,Error)
    ENDIF

    ! Check if we can precompute the BRDF Kernels
    IF(UseNewFixedParBRDF .AND. self%fixed_par .AND. .NOT. self%DoLambertian) THEN

      ! Do the allocation
      ALLOCATE(self%VBRDF_Kernels(self%nkern))
      DO n=1,self%nkern
        CALL AllocateVLBRDF(self%VBRDF_Kernels(n),Error)
      ENDDO

    ENDIF

  END SUBROUTINE InitSurface
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: SampleSurfaceProperties
  ! 
  ! DESCRIPTION: Computes BRDF spectrum for a given geolocation
  
  SUBROUTINE SampleSurfaceProperties( self, Window, Geolocation, SurfProf, L2Surface, Error )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    CLASS(SurfaceType),    INTENT(INOUT)   :: self
    TYPE(WinType),         INTENT(IN)      :: Window
    TYPE(GeolocationType), INTENT(IN)      :: Geolocation
    TYPE(SurfProfType),    INTENT(IN)      :: SurfProf
    TYPE(AuxSurfType),     INTENT(IN)      :: L2Surface
    ! TYPE(SurfaceType),     INTENT(INOUT)   :: Surface
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    ! =====================================================================
    ! SampleSurfaceProperties starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Sample the surface reflectance for the wavelength/geo-dependent cases
    IF(self%UseL2Reflectance) THEN
      CALL Surface0_RTMAlbedo(Window,L2Surface,self,Error)
    ELSEIF(self%OptionIndex .EQ. 1) THEN
      ! Do nothing
    ELSEIF(self%OptionIndex .EQ. 2) THEN
      CALL Surface2_RTMAlbedo(Window, self, Error)
    ELSEIF(self%OptionIndex .EQ. 3) THEN
      CALL Surface3_RTMAlbedo(Window, Geolocation, self, Error)
    ELSEIF(self%OptionIndex .EQ. 4) THEN
      CALL Surface4_RTMAlbedo(Window, Geolocation, SurfProf, self, Error)
    ELSEIF(self%OptionIndex .EQ. 5) THEN
      ! Do nothing - its fixed
      !CALL Surface5_RTMAlbedo(Window, Geolocation, Surface, Error)
    ELSEIF(self%OptionIndex .EQ. 6) THEN
      CALL Surface6_RTMAlbedo(Window, Geolocation, self, Error)
    ELSE
      STOP 'Surface Option must range from 1..6'
    ENDIF ; IF(CheckError(Error)) RETURN

    ! Update SIF
    IF(self%DoPlantFluorescence) THEN
      CALL self%SIF%ComputeSpectrum(Window%nRTM_wvl ,&
                      Window%RTM_wvl, SurfProf%SIF_734nm)
    ENDIF

    ! Emissivity Options
    IF(ALLOCATED(self%Emissivity%Value)) DEALLOCATE(self%Emissivity%Value)
    ALLOCATE(self%Emissivity%Value(Window%nRTM_Wvl))
    
    ! Constant value case
    IF(self%Emissivity%OptionIndex .EQ. 1) THEN
      self%Emissivity%Value(:) = self%Emissivity%ConstantEmissivity
    ELSEIF(self%Emissivity%OptionIndex .EQ. 2) THEN
      ! CALL Surface2_RTMEmiss
      STOP 'Surface Emissivity Option 2 has not yet been implemented'
    ELSEIF(self%Emissivity%OptionIndex .EQ. 3) THEN
      ! CALL Surface3_RTMEmiss
      STOP 'Surface Emissivity Option 3 has not yet been implemented'
    ELSE
      STOP 'Surface Emissivity Option must range from 1..3'
    ENDIF

    ! Check if we can precompute the BRDF Kernels
    IF(UseNewFixedParBRDF .AND. self%fixed_par .AND. .NOT. self%DoLambertian) THEN

      !
    ENDIF
    
    ! print*,'self%fixed_par',self%fixed_par
    ! STOP 'SampleSurfaceProperties'

    
  END SUBROUTINE SampleSurfaceProperties
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: SetSurfaceOpticalProperties
  ! 
  ! DESCRIPTION: Sets the surface optical properties for a particular
  !              wavelength from the VLIDORT RTM
  
  SUBROUTINE SetSurfaceOpticalProperties(self, w, Wavelength, Geolocation, I0_ConversionFactor, SurfProf,   &
                                         VLIDORT_FixIn, VLIDORT_LinModIn, VLIDORT_Sup, VLIDORT_LinSup, Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(SurfaceType),               INTENT(INOUT) :: self
    INTEGER,                          INTENT(IN)    :: w
    REAL(KIND=8),                     INTENT(IN)    :: Wavelength
    TYPE(GeolocationType),            INTENT(IN)    :: Geolocation
    REAL(KIND=8),                     INTENT(IN)    :: I0_conversionFactor
    TYPE(SurfProfType),               INTENT(IN)    :: SurfProf
    TYPE(VLIDORT_Fixed_Inputs),       INTENT(INOUT) :: VLIDORT_FixIn
    TYPE(VLIDORT_Modified_LinInputs), INTENT(INOUT) :: VLIDORT_LinModIn
    TYPE(VLIDORT_Sup_InOut),          INTENT(INOUT) :: VLIDORT_Sup ! VLIDORT supplements i/o structure
    TYPE(VLIDORT_LinSup_InOut),       INTENT(INOUT) :: VLIDORT_LinSup 
    TYPE(ErrorType),                  INTENT(INOUT) :: Error

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: n, j, w_a, w_p
    
    ! =====================================================================
    ! SetSurfaceOpticalProperties starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set Kernel Amplitude
    w_a = w ; IF(self%fixed_amp) w_a = 1
    DO n=1,self%nkern
      self%WavelengthKernelAmplitude(n) = self%kern_amp(w_a,n)
    ENDDO
    
    ! Added stuff for BRDF
    IF( .NOT. self%DoLambertian ) THEN

      ! First Compute BRDF arrays for cases where kernel varies by wavelength
      ! Its probably fine to precompute the kernels for the fixed par case - future change

      ! Set the kernel Parameters
      w_p = w ; IF(self%fixed_par) w_p = 1
      DO n=1,self%nkern
      DO j=1,self%npar(n)
        self%WavelengthKernelParameters(j,n) = self%kern_par(w_p,j,n)
      ENDDO
      ENDDO
      
      ! Recompute every time for safety (some kernels have wvl-dep pars)
      CALL ComputeBRDFKernels(Wavelength,Geolocation%SZA, Geolocation%VZA,&
                              Geolocation%AZA, self, Error, SurfProf     )

    ENDIF

    ! -------------------------------------------------------------------
    ! Transfer to VLIDORT Input Arrays
    ! -------------------------------------------------------------------

    ! Set the Lambertian Albedo Flag (can get reset if using lambertian clouds)
    VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = self%DoLambertian
    
    ! Lambertian Case
    IF( self%DoLambertian ) THEN
      
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = self%WavelengthKernelAmplitude(1)
      
    ! BRDF Case
    ELSE
      
      ! Set the VLIDORT Kernels
      VLIDORT_Sup%BRDF%TS_BRDF_F_0        = self%VLBRDF%BRDF_F_0
      VLIDORT_Sup%BRDF%TS_BRDF_F          = self%VLBRDF%BRDF_F
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = self%VLBRDF%USER_BRDF_F_0
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = self%VLBRDF%USER_BRDF_F
      VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = self%VLBRDF%DBOUNCE_BRDFUNC
      VLIDORT_Sup%BRDF%TS_EMISSIVITY      = self%VLBRDF%EMISSIVITY
      VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY = self%VLBRDF%USER_EMISSIVITY

      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = self%VLBRDF%LS_BRDF_F_0
      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = self%VLBRDF%LS_BRDF_F
      VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = self%VLBRDF%LS_USER_BRDF_F_0
      VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = self%VLBRDF%LS_USER_BRDF_F
      VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = self%VLBRDF%LS_DBOUNCE_BRDFUNC
      VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = self%VLBRDF%LS_USER_EMISSIVITY
      VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = self%VLBRDF%LS_EMISSIVITY
      
    
    ENDIF
    
    ! -------------------------------------------------------------------
    ! Solar Induced Plant Fluorescence
    ! -------------------------------------------------------------------
    IF(self%DoPlantFluorescence) THEN
    
      ! Switch on surface leaving
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING = .TRUE.
    
      ! Assume isotropic emission for now
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC = .TRUE.
    
      ! First zero isotropic leaving term
      VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(:,:) = 0.0
    
      ! Add leaving radiance
      VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1,1) = self%SIF%SIFSpectrum(w)*I0_ConversionFactor

    
    ENDIF
    
    ! Emissivity? - Need to input some sort of planck function but not sure of units
    ! Could possible just add as an isotropic source of radiation
    
    ! Based on lidort manual should be normalized to Solar Flux (same units)
    
!     IF(VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION) THEN
!       VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT = self%Emissivity%Value(w)
!     ENDIF
    
  END SUBROUTINE SetSurfaceOpticalProperties
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: SetSurfaceOpticalProperties
  ! 
  ! DESCRIPTION: Sets the surface optical properties for a particular
  !              wavelength from the 2-Stream RTM
  
  SUBROUTINE SetSurfaceOpticalProperties_2S(self, w, Wavelength, Geolocation, SurfProf, DO_BRDF_SURFACE,ALBEDO,&
                                            BRDF_F_0, BRDF_F, UBRDF_F, LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(SurfaceType),               INTENT(INOUT) :: self
    INTEGER,                          INTENT(IN)    :: w
    REAL(KIND=8),                     INTENT(IN)    :: Wavelength
    TYPE(GeolocationType),            INTENT(IN)    :: Geolocation
    TYPE(SurfProfType),               INTENT(IN)    :: SurfProf
    LOGICAL,                          INTENT(OUT)   :: DO_BRDF_SURFACE
    REAL(KIND=8),                     INTENT(OUT)   :: ALBEDO
    REAL(KIND=8),                     INTENT(OUT)   :: BRDF_F_0(0:1, MAX_GEOMETRIES)
    REAL(KIND=8),                     INTENT(OUT)   :: BRDF_F(0:1)
    REAL(KIND=8),                     INTENT(OUT)   :: UBRDF_F(0:1, MAX_GEOMETRIES)
    REAL(KIND=8),                     INTENT(OUT)   :: LS_BRDF_F_0(MAX_SURFACEWFS, 0:1, MAX_GEOMETRIES)
    REAL(KIND=8),                     INTENT(OUT)   :: LS_BRDF_F(MAX_SURFACEWFS, 0:1)
    REAL(KIND=8),                     INTENT(OUT)   :: LS_UBRDF_F(MAX_SURFACEWFS, 0:1, MAX_GEOMETRIES)
    TYPE(ErrorType),                  INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: n, j, w_a, w_p, v
    
    ! =====================================================================
    ! SetSurfaceOpticalProperties starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set Kernel Amplitude
    w_a = w ; IF(self%fixed_amp) w_a = 1
    DO n=1,self%nkern
      self%WavelengthKernelAmplitude(n) = self%kern_amp(w_a,n)
    ENDDO
    
    ! Initialize logical for doing BRDF
    DO_BRDF_SURFACE = .FALSE.

    ! Added stuff for BRDF
    IF( .NOT. self%DoLambertian ) THEN

      ! First Compute BRDF arrays for cases where kernel varies by wavelength
      ! Its probably fine to precompute the kernels for the fixed par case - future change
      
      DO_BRDF_SURFACE = .TRUE. 

      ! Set the kernel Parameters
      w_p = w ; IF(self%fixed_par) w_p = 1
      DO n=1,self%nkern
      DO j=1,self%npar(n)
        self%WavelengthKernelParameters(j,n) = self%kern_par(w_p,j,n)
      ENDDO
      ENDDO
      
      ! Recompute every time for safety (some kernels have wvl-dep pars)
      CALL ComputeBRDFKernels( Wavelength,Geolocation%SZA, Geolocation%VZA,&
                                 Geolocation%AZA, self, Error, SurfProf,   &
                                 nstreams=1, nstokes=1                     )

    ENDIF

    ! -------------------------------------------------------------------
    ! Transfer to LIDORT Input Arrays
    ! -------------------------------------------------------------------
    
    ! Lambertian Case
    IF( self%DoLambertian ) THEN
      ALBEDO = self%WavelengthKernelAmplitude(1)

    ! BRDF Case
    ELSE
      
      ! Set the VLIDORT Kernels (1 Geometry for now)
      BRDF_F_0(0:1,1) = self%VLBRDF%BRDF_F_0(0:1,1,1,1)
      BRDF_F(0:1) = self%VLBRDF%BRDF_F(0:1,1,1,1)
      UBRDF_F(0:1, 1) = self%VLBRDF%USER_BRDF_F(0:1,1,1,1)
      LS_BRDF_F_0(1:MAX_SURFACEWFS, 0:1, 1) = self%VLBRDF%LS_BRDF_F_0(1:MAX_SURFACEWFS,0:1,1,1,1)
      LS_BRDF_F(1:MAX_SURFACEWFS, 0:1) = self%VLBRDF%LS_BRDF_F(1:MAX_SURFACEWFS,0:1,1,1,1)
      LS_UBRDF_F(1:MAX_SURFACEWFS, 0:1, 1) = self%VLBRDF%LS_USER_BRDF_F(1:MAX_SURFACEWFS,0:1,1,1,1)
    
    ENDIF
    
    ! Emissivity? - Need to input some sort of planck function but not sure of units
    ! Could possible just add as an isotropic source of radiation
    
    ! Based on lidort manual should be normalized to Solar Flux (same units)
    
!     IF(VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION) THEN
!       VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT = Surface%Emissivity%Value(w)
!     ENDIF
    
  END SUBROUTINE SetSurfaceOpticalProperties_2S

  SUBROUTINE GetDirectBounceReflectance(self, w, Wavelength, Geolocation, SurfProf,&
                                        DBOUNCE,L_DBOUNCE, Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(SurfaceType),               INTENT(INOUT) :: self
    INTEGER,                          INTENT(IN)    :: w
    REAL(KIND=8),                     INTENT(IN)    :: Wavelength
    TYPE(GeolocationType),            INTENT(IN)    :: Geolocation
    TYPE(SurfProfType),               INTENT(IN)    :: SurfProf
    ! TYPE(SurfaceType),                INTENT(INOUT) :: Surface
    REAL(KIND=8),                     INTENT(OUT)   :: DBOUNCE(4,4,MAX_GEOMETRIES)
    REAL(KIND=8),                     INTENT(OUT)   :: L_DBOUNCE(4,4,MAX_GEOMETRIES,MAX_SURFACEWFS)
    TYPE(ErrorType),                  INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: n, i, j, w_a, w_p, k
    
    ! =====================================================================
    ! GetDirectBounceReflectance starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Zero the return
    DBOUNCE(:,:,:) = 0.0d0 ; L_DBOUNCE(:,:,:,:) = 0.0d0
    
    ! Set Kernel Amplitude
    w_a = w ; IF(self%fixed_amp) w_a = 1
    DO n=1,self%nkern
      self%WavelengthKernelAmplitude(n) = self%kern_amp(w_a,n)
    ENDDO
    
    ! Return the reflectance
    IF( self%DoLambertian ) THEN
      
      ! The kernel amplitude is the lambertian albedo
      DBOUNCE(1,1,:) = self%WavelengthKernelAmplitude(1)
      L_DBOUNCE(1,1,:,1) = 1.0d0
      
    ELSE

      ! First Compute BRDF arrays for cases where kernel varies by wavelength
      ! Its probably fine to precompute the kernels for the fixed par case - future change
      

      ! Set the kernel Parameters
      w_p = w ; IF(self%fixed_par) w_p = 1
      DO n=1,self%nkern
      DO j=1,self%npar(n)
        self%WavelengthKernelParameters(j,n) = self%kern_par(w_p,j,n)
      ENDDO
      ENDDO
      
      ! Recompute every time for safety (some kernels have wvl-dep pars)
      CALL ComputeBRDFKernels( Wavelength,Geolocation%SZA, Geolocation%VZA,&
                               Geolocation%AZA, self, Error, SurfProf      )
      
      ! Return the direct bounce
      ! CCM FIX - May not work for BPDF (I think the 16->4x4 unwraps this way 
      ! but need to check with Rob
      N = 1
      DO I=1,4
      DO J=1,4
        DBOUNCE(I,J,:) = self%VLBRDF%DBOUNCE_BRDFUNC(N,1,1,1)
        DO K=1,self%nJac
          L_DBOUNCE(I,J,:,K) = self%VLBRDF%LS_DBOUNCE_BRDFUNC(K,N,1,1,1)
        ENDDO
        N = N+1
      ENDDO
      ENDDO
      
    ENDIF
    
  END SUBROUTINE GetDirectBounceReflectance
  
  SUBROUTINE InitSurfaceOption0_L2(L2Surface, Window, Surface, Error)

    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(AuxSurfType),  INTENT(IN)    :: L2Surface
    TYPE(WinType),      INTENT(IN)    :: Window
    TYPE(SurfaceType),  INTENT(INOUT) :: Surface
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: k, idx

    ! =====================================================================
    ! InitSurfaceOption0_L2 starts here
    ! =====================================================================

    ! Set surface flags
    Surface%DoLambertian  = L2Surface%isLambertian
    IF(L2Surface%wmx .EQ. 1) THEN
      Surface%fixed_amp = .TRUE.
      Surface%fixed_par = .TRUE.
    ELSE
      Surface%fixed_amp = .FALSE.
      Surface%fixed_par = .FALSE.
    ENDIF
  
    ! Set VLIDORT Lambertian kernel
    Surface%nkern = L2Surface%kmx

    ! Set the variables associated with the kernel
    ALLOCATE(Surface%npar(Surface%nkern)) ; Surface%npar(1) = 0
    ALLOCATE(Surface%kern_name(Surface%nkern))    ; Surface%kern_name(1)    = BRDF_CHECK_NAMES(1)
    ALLOCATE(Surface%kern_idx(Surface%nkern))     ; Surface%kern_idx(1)     = 1
    
    DO k=1,L2Surface%kmx

      ! Kernel Idx
      idx = L2Surface%kern_idx(k)

      ! Set Kernel properties
      Surface%npar(k)      = vl_brdf_npar(idx)
      Surface%kern_name(k) = BRDF_CHECK_NAMES(idx)
      Surface%kern_idx(k)  = idx

    ENDDO

  END SUBROUTINE InitSurfaceOption0_L2
  
  SUBROUTINE Surface0_RTMAlbedo(Window,L2Surface,Surface,Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(WinType),         INTENT(IN)      :: Window
    TYPE(AuxSurfType),     INTENT(IN)      :: L2Surface
    TYPE(SurfaceType),     INTENT(INOUT)   :: Surface
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: k, p, idx, wmx, errstat

    ! =====================================================================
    ! Surface0_RTMAlbedo starts here
    ! =====================================================================

    ! fixed_amp = fixed_par (both on same grid)
    IF(Surface%fixed_amp) THEN 
      wmx = 1
    ELSE
      wmx = Window%nRTM_Wvl
    ENDIF

    ! Allocate Kernel Inputs for window
    CALL AllocateKernelInputs(wmx,L2Surface%kmx,L2Surface%pmx,Surface,Error)

    DO k=1,Surface%nkern

      ! Kernel Idx
      idx = L2Surface%kern_idx(k)

      ! These need to get reset (fix later)
      Surface%npar(k)      = vl_brdf_npar(idx)
      Surface%kern_name(k) = BRDF_CHECK_NAMES(idx)
      Surface%kern_idx(k)  = idx

    ENDDO

    ! Fixed Amplitudes/Parameters
    IF(Surface%fixed_amp) THEN

      ! If fixed then just copy the values
      DO k=1,Surface%nkern
        Surface%kern_amp(1,k)   = L2Surface%kern_amp(1,k)
        IF(L2Surface%pmx .GT. 0) THEN
          DO p=1,Surface%npar(k)
            Surface%kern_par(1,k,p) = L2Surface%kern_par(1,k,p)
          ENDDO
        ENDIF
      ENDDO
    
    ! Wavelength Dependent Amplitudes/Parameters
    ELSE

      IF(L2Surface%IsChebyshevParam) THEN

        ! Evaluate Chebyshev on the L1 Grid
        DO k=1,Surface%nkern

          ! Kernel Amplitude
          CALL EvaluateChebyshev(Window%nRTM_Wvl,L2Surface%WvlMin,L2Surface%WvlMax,   &
                                 Window%RTM_Wvl,L2Surface%wmx,L2Surface%kern_amp(:,k),&
                                 Surface%kern_amp(:,k),Error                          )

          ! Kernel Parameters
          IF(L2Surface%pmx .GT. 0) THEN
            DO p=1,Surface%npar(k)
              CALL EvaluateChebyshev(Window%nRTM_Wvl,L2Surface%WvlMin,L2Surface%WvlMax,    &
                                    Window%RTM_Wvl,L2Surface%wmx,L2Surface%kern_par(:,k,p),&
                                    Surface%kern_par(:,k,p), Error                         )
            ENDDO
          ENDIF
          
        ENDDO

      ELSE

        ! Interpolate to grid
        DO k=1,Surface%nkern

          ! Kernel Amplitude
          CALL BSPLINE_EdgeFill(L2Surface%wvl,L2Surface%kern_amp(:,k),L2Surface%wmx,&
                     Window%RTM_Wvl, Surface%kern_amp(:,k), Window%nRTM_Wvl, errstat)

          ! Kernel Parameters
          IF(L2Surface%pmx .GT. 0) THEN
            DO p=1,Surface%npar(k)
              CALL BSPLINE_EdgeFill(L2Surface%wvl,L2Surface%kern_par(:,k,p),L2Surface%wmx,&
                         Window%RTM_Wvl, Surface%kern_par(:,k,p), Window%nRTM_Wvl, errstat)
            ENDDO
          ENDIF
          
        ENDDO

      ENDIF

    ENDIF
    
  END SUBROUTINE Surface0_RTMAlbedo

  SUBROUTINE InitSurfaceOption1_FixedAlb(SurfOpt, Window, Surface, Error)
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(SurfOptType),  INTENT(IN)    :: SurfOpt
    TYPE(WinType),      INTENT(IN)    :: Window
    TYPE(SurfaceType),  INTENT(INOUT) :: Surface
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    ! =====================================================================
    ! InitSurfaceOption1_FixedAlb starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set surface flags
    Surface%DoLambertian  = .TRUE.
    Surface%fixed_par     = .TRUE.
    Surface%fixed_amp     = .TRUE.
    
    ! Set VLIDORT Lambertian kernel
    Surface%nkern = 1
    
    ! Set the fixed lambertian albedo
    ALLOCATE(Surface%npar(1))         ; Surface%npar(1) = 0
    ALLOCATE(Surface%kern_name(1))    ; Surface%kern_name(1)    = BRDF_CHECK_NAMES(1)
    ALLOCATE(Surface%kern_idx(1))     ; Surface%kern_idx(1)     = 1
    ALLOCATE(Surface%kern_amp(1,1))   ; Surface%kern_amp(1,1)   = SurfOpt%Option1_FixedAlbedo
    ALLOCATE(Surface%kern_par(1,1,1)) ; Surface%kern_par(1,1,1) = 0.0
    ALLOCATE(Surface%WavelengthKernelAmplitude(1))
    
  END SUBROUTINE InitSurfaceOption1_FixedAlb
  
  SUBROUTINE InitSurfaceOption2_AlbSpec(SurfOpt,Surface,Window,Error)
    
    ! Adapated from the GEOCAPE-TOOL
    USE interpolation_module, ONLY : SPLINE1
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(SurfOptType),  INTENT(IN)    :: SurfOpt
    TYPE(WinType),      INTENT(IN)    :: Window
    TYPE(SurfaceType),  INTENT(INOUT) :: Surface
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: I
    INTEGER :: nAlb
    INTEGER :: iw0, iwf
    
    REAL(KIND=8), ALLOCATABLE :: Wavelength(:), AlbSpec(:)
    
    ! =====================================================================
    ! InitSurfaceOption2_AlbSpec starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set surface flags
    Surface%DoLambertian  = .TRUE.
    Surface%fixed_par     = .TRUE.
    Surface%fixed_amp     = .FALSE.
    
    ! Set VLIDORT Lambertian kernel
    Surface%nkern = 1
    
    ! Open File
    OPEN(srfunit, FILE=TRIM(ADJUSTL(SurfOpt%Option2_Infile)), STATUS='old')
    READ(srfunit, *) nAlb
    
    ! Allocate Wavelength/Albedo spectra
    ALLOCATE(Wavelength(nAlb))
    ALLOCATE(AlbSpec(nAlb))
    
    ! Read Data
    DO I = 1, nAlb
       READ(srfunit, *) Wavelength(I), AlbSpec(I)
    ENDDO
    
    ! Rescale Values
    Wavelength = Wavelength * 1.0d3 ! microns -> nm
    AlbSpec    = AlbSpec * 1e-2     ! % -> frac
    
    ! Close file
    CLOSE(srfunit)
    
    ! Get index range covering window
    iw0 = MAXLOC( Wavelength, DIM=1, MASK= Wavelength .LT. Window%Settings%StartWvl )
    iwf = MINLOC( Wavelength, DIM=1, MASK= Wavelength .GT. Window%Settings%EndWvl )
    
    ! Allocate output arrays
    Surface%Option2%nWvl = iwf - iw0 + 1
    
    ! Expand to at least 5 points
    IF( nAlb .GE. 5 ) THEN
      DO WHILE(Surface%Option2%nWvl .LT. 5 )
        
        IF(iw0 .GT. 1   ) iw0 = iw0 - 1
        IF(iwf .LT. nAlb) iwf = iwf + 1
        
        Surface%Option2%nWvl = iwf - iw0 + 1
        
      ENDDO
    ENDIF
    
    ALLOCATE(Surface%Option2%Wavelength(Surface%Option2%nWvl))
    ALLOCATE(Surface%Option2%Albedo(Surface%Option2%nWvl))
    ALLOCATE(Surface%Option2%AlbedoSP(Surface%Option2%nWvl))
    
    ! Store spectrum
    Surface%Option2%Wavelength(:) = Wavelength(iw0:iwf)
    Surface%Option2%Albedo(:) = AlbSpec(iw0:iwf)
    
    ! Compute Basis Spline Coefficients
    CALL SPLINE1(Wavelength(iw0:iwf),  AlbSpec(iw0:iwf),       &
                 Surface%Option2%nWvl, Surface%Option2%AlbedoSP)
    
  END SUBROUTINE InitSurfaceOption2_AlbSpec
  
  SUBROUTINE Surface2_RTMAlbedo(Window, Surface, Error)
    
    USE interpolation_module, ONLY : SPLINT1
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(WinType),         INTENT(IN)      :: Window
    TYPE(SurfaceType),     INTENT(INOUT)   :: Surface
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    ! =====================================================================
    ! Surface2_RTMAlbedo starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Deallocate arrays if necessary
    IF(ALLOCATED(Surface%npar))      DEALLOCATE(Surface%npar)
    IF(ALLOCATED(Surface%kern_name)) DEALLOCATE(Surface%kern_name)
    IF(ALLOCATED(Surface%kern_idx))  DEALLOCATE(Surface%kern_idx)
    IF(ALLOCATED(Surface%kern_amp))  DEALLOCATE(Surface%kern_amp)
    IF(ALLOCATED(Surface%kern_par))  DEALLOCATE(Surface%kern_par)
    IF(ALLOCATED(Surface%WavelengthKernelAmplitude)) DEALLOCATE(Surface%WavelengthKernelAmplitude)
    
    ! Set the fixed lambertian albedo
    ALLOCATE(Surface%npar(1))         ; Surface%npar(1) = 0
    ALLOCATE(Surface%kern_name(1))    ; Surface%kern_name(1)    = BRDF_CHECK_NAMES(1)
    ALLOCATE(Surface%kern_idx(1))     ; Surface%kern_idx(1)     = 1
    ALLOCATE(Surface%kern_par(1,1,1)) ; Surface%kern_par(1,1,1) = 0.0
    ALLOCATE(Surface%WavelengthKernelAmplitude(1))
    
    ! Set the kernel amplitudes via spline interpolation
    ALLOCATE(Surface%kern_amp(Window%nRTM_Wvl,1))
    
    CALL SPLINT1(Surface%Option2%Wavelength,           &
                 Surface%Option2%Albedo,               &
                 Surface%Option2%AlbedoSP,             &
                 Surface%Option2%nWvl,                 &
                 Window%RTM_Wvl,                       &
                 Surface%kern_amp(1:Window%nRTM_Wvl,1),&
                 Window%nRTM_Wvl)
    
  END SUBROUTINE Surface2_RTMAlbedo
  
  SUBROUTINE InitSurfaceOption3_LERClim(SurfOpt,Surface,Window,Error)
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(SurfOptType),  INTENT(IN)    :: SurfOpt
    TYPE(WinType),      INTENT(IN)    :: Window
    TYPE(SurfaceType),  INTENT(INOUT) :: Surface
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                :: rcode, vid, N, dimid, maxdim
    CHARACTER(LEN=maxChar) :: tmpchar

    REAL(KIND=4), ALLOCATABLE :: tmp_floatarr(:)

    ! =====================================================================
    ! InitSurfaceOption3_LERClim starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set surface flags
    Surface%DoLambertian  = .TRUE.
    Surface%fixed_par     = .TRUE.
    Surface%fixed_amp     = SurfOpt%Option3_UseConstantWvl
    
    ! Set VLIDORT Lambertian kernel
    Surface%nkern = 1

    ! Initialize i,j
    Surface%Option3%i = 0
    Surface%Option3%j = 0
    
    ! Store the fixed wavelength
    Surface%Option3%ConstWvl(1) = SurfOpt%Option3_ConstWvl
    
    ! Open File
    rcode = nf_open(TRIM(ADJUSTL(SurfOpt%Option3_Infile)),&
                    nf_Share,Surface%Option3%ncid         )
    
    ! Get Dimensions
    ! --------------

    ! Wavelength
    rcode = nf_inq_varid(Surface%Option3%ncid, 'Wavelength', vid)
    rcode = nf_inq_vardimid(Surface%Option3%ncid, vid, dimid)
    rcode = nf_inq_dim(Surface%Option3%ncid, dimid, tmpchar, &
                       Surface%Option3%wmx)
    maxdim = Surface%Option3%wmx


    ! Longitude
    rcode = nf_inq_varid(Surface%Option3%ncid, 'xmid', vid)
    rcode = nf_inq_vardimid(Surface%Option3%ncid, vid, dimid)
    rcode = nf_inq_dim(Surface%Option3%ncid, dimid, tmpchar, &
                       Surface%Option3%imx)
    maxdim = MAX(maxdim, Surface%Option3%imx)

    ! Latitude
    rcode = nf_inq_varid(Surface%Option3%ncid, 'ymid', vid)
    rcode = nf_inq_vardimid(Surface%Option3%ncid, vid, dimid)
    rcode = nf_inq_dim(Surface%Option3%ncid, dimid, tmpchar, &
                       Surface%Option3%jmx)
    maxdim = MAX(maxdim, Surface%Option3%jmx)

    ! Allocate Arrays
    ! ---------------
    ALLOCATE(Surface%Option3%Wavelength(Surface%Option3%wmx))
    ALLOCATE(Surface%Option3%Albedo(Surface%Option3%wmx))
    ALLOCATE(Surface%Option3%AlbedoSP(Surface%Option3%wmx))
    ALLOCATE(Surface%Option3%Longitude(Surface%Option3%imx))
    ALLOCATE(Surface%Option3%Latitude(Surface%Option3%jmx))

    ! Allocate temporary float array for reading
    ALLOCATE(tmp_floatarr(maxdim))

    ! Read Coordinates
    ! ----------------

    ! Wavelength
    rcode = nf_inq_varid(Surface%Option3%ncid, 'Wavelength', vid)
    rcode = nf_get_vara_real(Surface%Option3%ncid, vid,    &
                              (/1/),(/Surface%Option3%wmx/),&
                              tmp_floatarr(1:Surface%Option3%wmx))
    Surface%Option3%Wavelength = REAL(tmp_floatarr(1:Surface%Option3%wmx),KIND=8)

    ! Longitude
    rcode = nf_inq_varid(Surface%Option3%ncid, 'xmid', vid)
    rcode = nf_get_vara_real(Surface%Option3%ncid, vid,    &
                              (/1/),(/Surface%Option3%imx/),&
                              tmp_floatarr(1:Surface%Option3%imx))
    Surface%Option3%Longitude = REAL(tmp_floatarr(1:Surface%Option3%imx),KIND=8)

    ! Latitude
    rcode = nf_inq_varid(Surface%Option3%ncid, 'ymid', vid)
    rcode = nf_get_vara_real(Surface%Option3%ncid, vid,    &
                              (/1/),(/Surface%Option3%jmx/),&
                              tmp_floatarr(1:Surface%Option3%jmx))
    Surface%Option3%Latitude = REAL(tmp_floatarr(1:Surface%Option3%jmx),KIND=8)

    ! Deallocate temporary read array
    DEALLOCATE(tmp_floatarr)

  END SUBROUTINE InitSurfaceOption3_LERClim
  
  SUBROUTINE Surface3_RTMAlbedo(Window, Geolocation, Surface, Error)
    
    USE interpolation_module, ONLY : SPLINT1, SPLINE1, BSPLINE
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(WinType),         INTENT(IN)      :: Window
    TYPE(GeolocationType), INTENT(IN)      :: Geolocation
    TYPE(SurfaceType),     INTENT(INOUT)   :: Surface
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: rcode, vid, N, nwvl, errstat
    INTEGER :: idx(1), idy(1)
    REAL(KIND=4), ALLOCATABLE :: tmp_floatarr(:)
    REAL(KIND=8), ALLOCATABLE :: TruncWavelength(:)
    REAL(KIND=8)              :: ConstAlb(1)
    
    ! =====================================================================
    ! Surface3_RTMAlbedo starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Check which case
    IF( Surface%fixed_amp ) THEN
      nwvl =  1
    ELSE
      nwvl =  Window%nRTM_wvl
    ENDIF
    
    ! Allocate kernel inputs
    CALL AllocateKernelInputs( nwvl, 1, 1, Surface, Error )
    
    ! Set the fixed lambertian albedo
    Surface%npar(1)         = 0
    Surface%kern_name(1)    = BRDF_CHECK_NAMES(1)
    Surface%kern_idx(1)     = 1
    Surface%kern_par(1,1,1) = 0.0

    ! Find Nearest Index
    idx = MINLOC(ABS(Surface%Option3%Longitude - Geolocation%Longitude))
    idy = MINLOC(ABS(Surface%Option3%Latitude - Geolocation%Latitude))

    ! Check if we need to load a new spectrum
    IF( idx(1) .NE. Surface%Option3%i .OR.              &
        idy(1) .NE. Surface%Option3%j .OR.              &
        Surface%Option3%Month .NE. Geolocation%Time%Month ) THEN

      ! Allocate the temporary array
      ALLOCATE(tmp_floatarr(Surface%Option3%wmx))

      ! Read 
      rcode = nf_inq_varid(Surface%Option3%ncid, 'LER', vid)
      rcode = nf_get_vara_real(Surface%Option3%ncid, vid,    &
                                (/idx(1),idy(1),1,Geolocation%Time%Month/),&
                                (/1,1,Surface%Option3%wmx,1/),&
                                tmp_floatarr(1:Surface%Option3%wmx))
      Surface%Option3%Albedo = REAL(tmp_floatarr(1:Surface%Option3%wmx),KIND=8)
      
      IF( Surface%fixed_amp ) THEN
        CALL BSPLINE(Surface%Option3%Wavelength,Surface%Option3%Albedo,Surface%Option3%wmx,&
                     Surface%Option3%ConstWvl, ConstAlb, 1, errstat)
        Surface%kern_amp(1,1)   = ConstAlb(1)
      ELSE
      
        ! Compute Basis Spline Coefficients
        CALL SPLINE1(Surface%Option3%Wavelength,Surface%Option3%Albedo,&
                     Surface%Option3%wmx, Surface%Option3%AlbedoSP     )
      
      ENDIF
      
      ! Deallocate temporary read array
      DEALLOCATE(tmp_floatarr)
      
      ! Update i,j
      Surface%Option3%i = idx(1) 
      Surface%Option3%j = idy(1) 
    ENDIF
    
    ! Interpolate to RTM Grid for non-fixed wavelength case
    IF( .NOT. Surface%fixed_amp ) THEN
      ALLOCATE(TruncWavelength(Window%nRTM_Wvl)) ; TruncWavelength =  Window%RTM_Wvl
      DO N=1,Window%nRTM_Wvl
        IF(TruncWavelength(N) .LT. Surface%Option3%Wavelength(1)) THEN
          TruncWavelength(N) = Surface%Option3%Wavelength(1)
        ENDIF
        IF( TruncWavelength(N) .GT. Surface%Option3%Wavelength(Surface%Option3%wmx) ) THEN
          TruncWavelength(N) = Surface%Option3%Wavelength(Surface%Option3%wmx) 
        ENDIF
      ENDDO

      ! Interpolate to RTM Grid
      CALL SPLINT1(Surface%Option3%Wavelength,           &
                  Surface%Option3%Albedo,                &
                  Surface%Option3%AlbedoSP,              &
                  Surface%Option3%wmx,                   &
                  TruncWavelength,                       &
                  Surface%kern_amp(1:Window%nRTM_Wvl,1), &
                  Window%nRTM_Wvl)
      
      ! Deallocate array
      DEALLOCATE(TruncWavelength)
      
    ENDIF

  END SUBROUTINE Surface3_RTMAlbedo
  
  
  SUBROUTINE InitSurfaceOption4_SciaFA(SurfOpt,Surface,Window,Error)
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(SurfOptType),   INTENT(IN)    :: SurfOpt
    TYPE(WinType),      INTENT(IN)    :: Window
    TYPE(SurfaceType),  INTENT(INOUT) :: Surface
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    ! =====================================================================
    ! InitSurfaceOption4_SciaFA starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! -----------------------------
    ! Set Options for brdf sampling
    ! -----------------------------
    Surface%DoLambertian  = SurfOpt%Option4_DoIsotropic
    Surface%fixed_par     = .TRUE.
    Surface%fixed_amp     = .FALSE.
    
    ! Set Options for spatial FA Module
    Surface%Option4%Infile       = SurfOpt%Option4_Infile
    Surface%Option4%RootDataDir  = SurfOpt%RootDataDir
    Surface%Option4%DoIsotropic  = SurfOpt%Option4_DoIsotropic
    Surface%Option4%WhichAlbedo  = SurfOpt%Option4_WhichAlbedo
    Surface%Option4%DoOceanGlint = SurfOpt%Option4_DoOceanGlint
    Surface%Option4%RootBRDFDir  = TRIM(ADJUSTL(SurfOpt%RootDataDir))//SurfOpt%Option4_ClimDir
 
    ! Initialize Spatial FA
    CALL InitSpatialFA( Surface%Option4, Error )
    
    ! Allocate parameters
    IF(ALLOCATED(Surface%npar)) DEALLOCATE(Surface%npar)
    ALLOCATE(Surface%npar(Surface%Option4%nkern))
    
    ! Allocate BRDF arrays for VLIDORT
    IF(.NOT. Surface%DoLambertian) THEN
      CALL AllocateVLBRDF(Surface%VLBRDF, Error)
    ENDIF

    
    ! Copy surface options to main
    Surface%nkern = Surface%Option4%nkern
    IF(ALLOCATED(Surface%npar)) DEALLOCATE(Surface%npar)
    ALLOCATE(Surface%nPar(Surface%nkern))
    Surface%npar(:)  = Surface%Option4%npar(:)
    IF(ALLOCATED(Surface%kern_idx)) DEALLOCATE(Surface%kern_idx)
    ALLOCATE(Surface%kern_idx(Surface%nkern))
    Surface%kern_idx(:) = Surface%Option4%KernIdx(:)
    
  END SUBROUTINE InitSurfaceOption4_SciaFA
  
  SUBROUTINE Surface4_RTMAlbedo(Window, Geolocation, SurfProf, Surface, Error)
    
    USE interpolation_module, ONLY : SPLINT1
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(WinType),         INTENT(IN)      :: Window
    TYPE(GeolocationType), INTENT(IN)      :: Geolocation
    TYPE(SurfProfType),    INTENT(IN)      :: SurfProf
    TYPE(SurfaceType),     INTENT(INOUT)   :: Surface
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    REAL(KIND=8), ALLOCATABLE :: TruncWavelength(:)
    INTEGER                   :: N, J
    ! =====================================================================
    ! Surface4_RTMAlbedo starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! -----------------------------------------------------------
    ! Compute the kernels on the FA Wavelength Grid
    ! -----------------------------------------------------------
    CALL SampleSpatialFA( Geolocation, SurfProf, Surface%Option4, Error )
    IF(CheckError(Error)) RETURN

    ! -----------------------------------------------------------
    ! Truncated wavelength grid
    ! -----------------------------------------------------------

    ALLOCATE(TruncWavelength(Window%nRTM_Wvl)) ; TruncWavelength =  Window%RTM_Wvl
    DO N=1,Window%nRTM_Wvl
      IF(TruncWavelength(N) .LT. Surface%Option4%Wvl(1)) THEN
        TruncWavelength(N) = Surface%Option4%Wvl(1)
      ENDIF
      IF( TruncWavelength(N) .GT. Surface%Option4%Wvl(Surface%Option4%wmx) ) THEN
        TruncWavelength(N) = Surface%Option4%Wvl(Surface%Option4%wmx) 
      ENDIF
    ENDDO

    ! -----------------------------------------------------------
    ! Set Kernel properties
    ! -----------------------------------------------------------
    
    ! Allocate Kernel Inputs for window
    CALL AllocateKernelInputs(Window%nRTM_Wvl,       &
                              Surface%Option4%nkern, &
                              Surface%Option4%maxpar,&
                              Surface, Error         )
    
    ! Spline Amplitudes to the window grid
    DO n=1,Surface%Option4%nkern
      
      ! Kernel Index
      Surface%kern_idx(n) = Surface%Option4%KernIdx(n)
      
      ! Kernel name
      Surface%kern_name(n) = BRDF_CHECK_NAMES(Surface%kern_idx(n))
      
      ! Number of kernel parameters
      DO j=1,Surface%Option4%nkern
        Surface%npar(j) = Surface%Option4%npar(j)
      ENDDO
      
      ! Kernel Parameters
      DO j=1,Surface%npar(n)
        Surface%kern_par(:,j,n) = Surface%Option4%KernPar(j,n)
      ENDDO

      ! Kernel Amplitudes
      CALL SPLINT1(Surface%Option4%Wvl,                    &
                   Surface%Option4%KernelAmplitudes(:,n),  &
                   Surface%Option4%KernelAmplitudesSP(:,n),&
                   Surface%Option4%wmx,                    &
                   TruncWavelength,                        &
                   Surface%kern_amp(1:Window%nRTM_wvl,n),  &
                   Window%nRTM_wvl                         )
      
    ENDDO

    ! ! Check the amplitudes
    ! DO j=1, Window%nRTM_wvl
    !     WRITE(992,*) Window%RTM_Wvl(j),Surface%kern_amp(j,1),&
    !                  Surface%kern_amp(j,2),Surface%kern_amp(j,3),&
    !                  Surface%kern_amp(j,4)
    ! ENDDO
    ! STOP 'Testing SpatialFA'
    
    DEALLOCATE(TruncWavelength)
    
  END SUBROUTINE Surface4_RTMAlbedo
  
  SUBROUTINE InitSurfaceOption5_FixedBRDF(SurfOpt,Surface,Window,Error)
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(SurfOptType),  INTENT(IN)    :: SurfOpt
    TYPE(WinType),      INTENT(IN)    :: Window
    TYPE(SurfaceType),  INTENT(INOUT) :: Surface
    TYPE(ErrorType),    INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER      :: max_par
    INTEGER      :: N, J
    REAL(KIND=8) :: Wavelength
    
    ! =====================================================================
    ! InitSurfaceOption5_FixedBRDF starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set Flags
    Surface%fixed_amp    = .TRUE.
    Surface%fixed_par    = .TRUE.
    Surface%DoLambertian = .FALSE.

    ! Maximum number of parameters
    max_par = MAXVAL( SurfOpt%Option5_nKernPar )
    max_par = MAX(max_par,1)

    CALL AllocateKernelInputs( 1, 3, max_par, Surface, Error )

    ! Set Kernel values
    Surface%nkern = 3
    Surface%npar = SurfOpt%Option5_nKernPar
    Surface%kern_idx = SurfOpt%Option5_KernIdx
    DO N=1,Surface%nkern
      Surface%kern_name(N) = BRDF_CHECK_NAMES(Surface%kern_idx(N))
      DO J=1,Surface%npar(N)
        Surface%kern_par(1,:,N) = SurfOpt%Option5_KernPar(J,N)
      ENDDO
    ENDDO
    Surface%kern_amp(1,:) = SurfOpt%Option5_KernAmp(:)
    Surface%WavelengthKernelAmplitude(:) = SurfOpt%Option5_KernAmp(:)
    Surface%WavelengthKernelParameters(1:max_par,:) = SurfOpt%Option5_KernPar(1:max_par,:)
    
    ! Allocate BRDF arrays for VLIDORT
    CALL AllocateVLBRDF(Surface%VLBRDF, Error)
    
  END SUBROUTINE InitSurfaceOption5_FixedBRDF
  
  SUBROUTINE InitSurfaceOption6_ClimBRDF(SurfOpt,Surface,Window,Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(SurfOptType),  INTENT(IN)      :: SurfOpt
    TYPE(WinType),      INTENT(IN)      :: Window
    TYPE(SurfaceType),  INTENT(INOUT)   :: Surface
    TYPE(ErrorType),    INTENT(INOUT)   :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, vid, dimid, dimid5(5)
    CHARACTER(LEN=maxChar) :: tmpchar
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitSurfaceOption6_ClimBRDF'

    ! =====================================================================
    ! InitSurfaceOption6_ClimBRDF starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! -----------------------------
    ! Set Options for brdf sampling
    ! -----------------------------
    Surface%DoLambertian  = .FALSE.
    Surface%fixed_par     = .FALSE.
    Surface%fixed_amp     = .FALSE.
    
    ! Initialize Coordinates
    Surface%Option6%i = 0
    Surface%Option6%j = 0
    
    ! Allocate BRDF arrays for VLIDORT
    CALL AllocateVLBRDF(Surface%VLBRDF, Error)
    
    ! -----------
    ! Attach file
    ! -----------
    
    rcode = nf_open(TRIM(ADJUSTL(SurfOpt%Option6_Infile)),&
                    nf_Share,Surface%Option6%ncid         )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_open')

    ! Climatology wavelength dimension
    rcode = nf_inq_varid(Surface%Option6%ncid, 'clim_wvl', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:clim_wvl')
    rcode = nf_inq_vardimid(Surface%Option6%ncid,   vid,  dimid  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dimid:clim_wvl')
    rcode = nf_inq_dim(Surface%Option6%ncid, dimid,tmpchar, &
                       Surface%Option6%wmx                  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:clim_wvl')

    ! Climatology Longitude Dimension
    rcode = nf_inq_varid(Surface%Option6%ncid, 'clim_lon', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:clim_lon')
    rcode = nf_inq_vardimid(Surface%Option6%ncid,   vid,  dimid  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dimid:clim_lon')
    rcode = nf_inq_dim(Surface%Option6%ncid, dimid,tmpchar, &
                       Surface%Option6%imx                  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:clim_lon')

    ! Climatology latitude dimension
    rcode = nf_inq_varid(Surface%Option6%ncid, 'clim_lat', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:clim_lat')
    rcode = nf_inq_vardimid(Surface%Option6%ncid,   vid,  dimid  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dimid:clim_lat')
    rcode =      nf_inq_dim(Surface%Option6%ncid, dimid,tmpchar, &
                            Surface%Option6%jmx                  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:clim_lat')

    ! Climatology number of kernels
    rcode = nf_inq_varid(Surface%Option6%ncid, 'fkern_idx', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:fkern_idx')
    rcode = nf_inq_vardimid(Surface%Option6%ncid,   vid,  dimid   )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dimid:fkern_idx')
    rcode = nf_inq_dim(Surface%Option6%ncid, dimid,tmpchar,  &
                       Surface%nkern                         )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:fkern_idx')
    
    ! Climatology - max number of parameters
    rcode = nf_inq_varid(Surface%Option6%ncid, 'fkern_par', vid  )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:fkern_par')
    rcode = nf_inq_vardimid(Surface%Option6%ncid,   vid,  dimid5    )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dimid:fkern_par')
    rcode = nf_inq_dim(Surface%Option6%ncid, dimid5(5),tmpchar,&
                       Surface%Option6%pmx                     )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:fkern_par')
    
    ! ---------------
    ! Allocate arrays
    ! ---------------
    ALLOCATE(Surface%Option6%Wavelength(Surface%Option6%wmx))
    ALLOCATE(Surface%Option6%Longitude(Surface%Option6%imx))
    ALLOCATE(Surface%Option6%Latitude(Surface%Option6%jmx))
    ALLOCATE(Surface%Option6%KernIdx(Surface%nkern))
    ALLOCATE(Surface%Option6%Amplitude(Surface%Option6%wmx,Surface%nkern))
    ALLOCATE(Surface%Option6%AmplitudeSP(Surface%Option6%wmx,Surface%nkern))
    ALLOCATE(Surface%Option6%Parameters(Surface%Option6%wmx,Surface%Option6%pmx,Surface%nkern))
    ALLOCATE(Surface%Option6%ParametersSP(Surface%Option6%wmx,Surface%Option6%pmx,Surface%nkern))

    ! ---------
    ! Read Grid
    ! ---------
    
    ! Wavelength
    rcode = nf_inq_varid(Surface%Option6%ncid, 'clim_wvl', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:clim_wvl')
    rcode = nf_get_var_double(Surface%Option6%ncid, vid, &
                              Surface%Option6%Wavelength )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:clim_wvl')

    ! Longitude
    rcode = nf_inq_varid(Surface%Option6%ncid, 'clim_lon', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:clim_lon')
    rcode = nf_get_var_double(Surface%Option6%ncid, vid, &
                              Surface%Option6%Longitude )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:clim_lon')

    ! Latitude
    rcode = nf_inq_varid(Surface%Option6%ncid, 'clim_lat', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:clim_lat')
    rcode = nf_get_var_double(Surface%Option6%ncid, vid, &
                              Surface%Option6%Latitude )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:clim_lat')
    
    ! Kernel Indices
    rcode = nf_inq_varid(Surface%Option6%ncid, 'fkern_idx', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:fkern_idx')
    rcode = nf_get_var_int(Surface%Option6%ncid, vid, &
                           Surface%Option6%KernIdx )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:fkern_idx')

    ! We must figure out the appropriate dimensions now for figuring out the output 
    ! linearizations for diagnostic initialization
    CALL SetBRDFNamesAndDimensions(Surface,Surface%nkern,Surface%Option6%KernIdx)
    
  END SUBROUTINE InitSurfaceOption6_ClimBRDF
  
  SUBROUTINE Surface6_RTMAlbedo(Window, Geolocation, Surface, Error)
    
    USE interpolation_module, ONLY : SPLINT1, SPLINE1
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(WinType),         INTENT(IN)      :: Window
    TYPE(GeolocationType), INTENT(IN)      :: Geolocation
    TYPE(SurfaceType),     INTENT(INOUT)   :: Surface
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: rcode, vid, N,I,J, nwvl, errstat, nw, nk, np
    INTEGER :: idx(1), idy(1)
    REAL(KIND=4), ALLOCATABLE :: clim_amp_r4(:,:), clim_par_r4(:,:,:)
    REAL(KIND=8), ALLOCATABLE :: TruncWavelength(:)
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'Surface6_RTMAlbedo'

    ! =====================================================================
    ! Surface6_RTMAlbedo starts here 
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! For convenience
    nw = Surface%Option6%wmx
    nk = Surface%nkern
    np = Surface%Option6%pmx
    
    ! Allocate kernel inputs
    CALL AllocateKernelInputs( Window%nRTM_wvl, Surface%nkern, Surface%Option6%pmx, Surface, Error )

    ! Ensure kernel index set
    Surface%kern_idx = Surface%Option6%KernIdx
    
    ! Currently just nearest neighbour sample
    idx = MINLOC(ABS(Surface%Option6%Longitude - Geolocation%Longitude))
    idy = MINLOC(ABS(Surface%Option6%Latitude - Geolocation%Latitude))
    
    ! Check if we need to load a new spectrum
    IF( idx(1) .NE. Surface%Option6%i .OR.  &
        idy(1) .NE. Surface%Option6%j       ) THEN

      ! Allocate the temporary array
      ALLOCATE(clim_amp_r4(Surface%Option6%wmx,Surface%nkern))
      ALLOCATE(clim_par_r4(Surface%Option6%wmx,Surface%Option6%pmx,Surface%nkern))
      
      ! Read kernel amplitudes
      rcode = nf_inq_varid(Surface%Option6%ncid, 'fkern_amp', vid)
      rcode = nf_get_vara_real( Surface%Option6%ncid, vid,&
                             (/idx(1),idy(1), 1, 1/),     &
                             (/1,1,nw,nk/),               &
                              clim_amp_r4                 )
      Surface%Option6%Amplitude(:,:) = REAL(clim_amp_r4,KIND=8)
      
      ! Read kernel parameters
      rcode = nf_inq_varid(Surface%Option6%ncid, 'fkern_par', vid)
      rcode = nf_get_vara_real( Surface%Option6%ncid, vid,&
                             (/idx(1),idy(1), 1, 1, 1/),  &
                             (/      1,    1,nw,nk,np/),  &
                              clim_par_r4                 )
      Surface%Option6%Parameters(:,:,:) = REAL(clim_par_r4,KIND=8)
      
      ! Create spline objects
      DO I=1,Surface%nkern
        
        ! Compute Basis Spline Coefficients for Kernel amplitudes
        IF(nw .GT. 1) THEN
          CALL SPLINE1(Surface%Option6%Wavelength,Surface%Option6%Amplitude(:,I),&
                      Surface%Option6%wmx, Surface%Option6%AmplitudeSP(:,I)     )
          
          ! Loop over parameters
          DO J=1,vl_brdf_npar(Surface%Option6%KernIdx(I))
            
            ! Compute Basis Spline Coefficients for Kernel parameters
            CALL SPLINE1(Surface%Option6%Wavelength,Surface%Option6%Parameters(:,J,I),&
                        Surface%Option6%wmx, Surface%Option6%ParametersSP(:,J,I)     )
          ENDDO
        ENDIF
        
      ENDDO
      
      ! Update i,j
      Surface%Option6%i = idx(1) 
      Surface%Option6%j = idy(1) 
      
      ! Deallocate temporary reading arrays
      DEALLOCATE(clim_amp_r4) ; DEALLOCATE(clim_par_r4)
      
    ENDIF
    
    ! Set up wavelength grid to avoid extrapolation
    ALLOCATE(TruncWavelength(Window%nRTM_Wvl)) ; TruncWavelength =  Window%RTM_Wvl
    DO N=1,Window%nRTM_Wvl
      IF(TruncWavelength(N) .LT. Surface%Option6%Wavelength(1)) THEN
        TruncWavelength(N) = Surface%Option6%Wavelength(1)
      ENDIF
      IF( TruncWavelength(N) .GT. Surface%Option6%Wavelength(Surface%Option6%wmx) ) THEN
        TruncWavelength(N) = Surface%Option6%Wavelength(Surface%Option6%wmx) 
      ENDIF
    ENDDO
    
    DO I=1,Surface%nkern
      
      ! Set Kernel Index
      Surface%kern_idx(I) = Surface%Option6%KernIdx(I)
      
      ! Set Kernel Name
      Surface%kern_name(I)    = BRDF_CHECK_NAMES(I)
        
      ! Get Number of parameters
      Surface%npar(I) = vl_brdf_npar(Surface%kern_idx(I))
      
      IF(nw .GT. 1) THEN
        ! Interpolate Amplitude to RTM Wavelength Grid
        CALL SPLINT1(Surface%Option6%Wavelength,          &
                    Surface%Option6%Amplitude(:,I),       &
                    Surface%Option6%AmplitudeSP(:,I),     &
                    Surface%Option6%wmx,                  &
                    TruncWavelength,                      &
                    Surface%kern_amp(1:Window%nRTM_Wvl,I),&
                    Window%nRTM_Wvl                       )
        
        
        ! Interpolate Parameters to RTM Wavelength Grid
        DO J=1,Surface%npar(I)
          CALL SPLINT1(Surface%Option6%Wavelength,             &
                      Surface%Option6%Parameters(:,J,I),      &
                      Surface%Option6%ParametersSP(:,J,I),    &
                      Surface%Option6%wmx,                    &
                      TruncWavelength,                        &
                      Surface%kern_par(1:Window%nRTM_Wvl,J,I),&
                      Window%nRTM_Wvl                         )
        ENDDO
      ELSE

        ! Kernel amplitude
        Surface%kern_amp(1:Window%nRTM_Wvl,I) = Surface%Option6%Amplitude(1,I)

        ! Parameters
        DO J=1,Surface%npar(I)
          Surface%kern_par(1:Window%nRTM_Wvl,J,I) = Surface%Option6%Parameters(1,J,I)
        ENDDO

      ENDIF
      
      ! Deallocate array
      DEALLOCATE(TruncWavelength)
      
    ENDDO

  END SUBROUTINE Surface6_RTMAlbedo
  
  
  SUBROUTINE AllocateKernelInputs( nwvl, nkern, npar, Surface, Error )

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,            INTENT(IN)      :: nwvl
    INTEGER,            INTENT(IN)      :: nkern
    INTEGER,            INTENT(IN)      :: npar
    TYPE(SurfaceType),  INTENT(INOUT)   :: Surface
    TYPE(ErrorType),    INTENT(INOUT)   :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: idx, idy

    ! =====================================================================
    ! AllocateKernelInputs starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Deallocate arrays if necessary
    IF(ALLOCATED(Surface%npar))      DEALLOCATE(Surface%npar)
    IF(ALLOCATED(Surface%kern_name)) DEALLOCATE(Surface%kern_name)
    IF(ALLOCATED(Surface%kern_idx))  DEALLOCATE(Surface%kern_idx)
    IF(ALLOCATED(Surface%kern_amp))  DEALLOCATE(Surface%kern_amp)
    IF(ALLOCATED(Surface%kern_par))  DEALLOCATE(Surface%kern_par)
    IF(ALLOCATED(Surface%WavelengthKernelAmplitude)) DEALLOCATE(Surface%WavelengthKernelAmplitude)
    IF(ALLOCATED(Surface%WavelengthKernelParameters)) DEALLOCATE(Surface%WavelengthKernelParameters)
    
    ! Set the fixed lambertian albedo
    ALLOCATE(Surface%npar(nkern))      ; Surface%npar(:) = 0
    ALLOCATE(Surface%kern_name(nkern)) ; Surface%kern_name(:)    = ''
    ALLOCATE(Surface%kern_idx(nkern))  ; Surface%kern_idx(:)     = 0
    ALLOCATE(Surface%kern_par(nwvl,npar,nkern)) ; Surface%kern_par(:,:,:) = 0.0d0
    ALLOCATE(Surface%WavelengthKernelAmplitude(nkern)) ; Surface%WavelengthKernelAmplitude(:) = 0.0d0
    ALLOCATE(Surface%WavelengthKernelParameters(npar,nkern)) ; Surface%WavelengthKernelParameters(:,:) = 0.0d0
    ALLOCATE(Surface%kern_amp(nwvl,nkern)) ; Surface%kern_amp(:,:) = 0.0d0
    
  END SUBROUTINE AllocateKernelInputs

  
  SUBROUTINE SetBRDFNamesAndDimensions(Surface,nkern,fkern_idx)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(SurfaceType),      INTENT(INOUT) :: Surface
    INTEGER,                INTENT(IN)    :: nkern
    INTEGER,                INTENT(IN)    :: fkern_idx(nkern)

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: k


    ! =====================================================================
    ! SetBRDFNamesAndDimensions starts here
    ! =====================================================================
    
    ! Set # kernels
    Surface%nkern = nkern

    ! Set the kernel indices
    IF(ALLOCATED(Surface%kern_idx)) DEALLOCATE(Surface%kern_idx)
    ALLOCATE(Surface%kern_idx(nkern))
    Surface%kern_idx(:) = fkern_idx(:)
    
    ! Allocate arrays for storing parameter dimensions/brdf names
    IF(ALLOCATED(Surface%npar)) DEALLOCATE(Surface%npar)
    ALLOCATE(Surface%nPar(Surface%nkern))
    IF(ALLOCATED(Surface%kern_name)) DEALLOCATE(Surface%kern_name)
    ALLOCATE(Surface%kern_name(Surface%nkern))


    ! Set Names and dimensions
    DO k=1,nkern

      ! Get BRDF Name
      Surface%kern_name(k) = BRDF_CHECK_NAMES(fkern_idx(k))

      ! Get BRDF Kernel Dimension
      Surface%npar(k) = vl_brdf_npar(fkern_idx(k))

    ENDDO
    
  END SUBROUTINE SetBRDFNamesAndDimensions

  SUBROUTINE AllocateVLBRDF(VLBRDF,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(VLBRDFInpType), INTENT(INOUT) :: VLBRDF
    TYPE(ErrorType),     INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------


    ! =====================================================================
    ! AllocateVLIDORT_BRDF starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Later apply logic to only allocate what is needed
    ALLOCATE(VLBRDF%BRDF_F_0(0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS,  MAXBEAMS))
    ALLOCATE(VLBRDF%BRDF_F(0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS))
    ALLOCATE(VLBRDF%USER_BRDF_F_0(0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS))
    ALLOCATE(VLBRDF%USER_BRDF_F(0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS))
    ALLOCATE(VLBRDF%DBOUNCE_BRDFUNC(MAXSTOKES_SQ,MAX_USER_STREAMS,MAX_USER_RELAZMS, MAXBEAMS))
    ALLOCATE(VLBRDF%EMISSIVITY(MAXSTOKES, MAXSTREAMS))
    ALLOCATE(VLBRDF%USER_EMISSIVITY(MAXSTOKES, MAX_USER_STREAMS))
    ALLOCATE(VLBRDF%LS_BRDF_F_0(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS))
    ALLOCATE(VLBRDF%LS_BRDF_F(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS))
    ALLOCATE(VLBRDF%LS_USER_BRDF_F_0(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS))
    ALLOCATE(VLBRDF%LS_USER_BRDF_F(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS ))
    ALLOCATE(VLBRDF%LS_DBOUNCE_BRDFUNC(MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS))
    ALLOCATE(VLBRDF%LS_USER_EMISSIVITY(MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS))
    ALLOCATE(VLBRDF%LS_EMISSIVITY(MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS))

    
  END SUBROUTINE AllocateVLBRDF

  SUBROUTINE InitializeBRDFInputs( Surface, VBRDF_Sup_In, VBRDF_LinSup_In, VBRDF_Sup_InputStatus, Error )
    
    ! --------------------
    ! subroutine arguments
    ! ---------- ----------
    TYPE(SurfaceType),                    INTENT(INOUT) :: Surface
    TYPE(VBRDF_Sup_Inputs),               INTENT(INOUT) :: VBRDF_Sup_In
    TYPE(VBRDF_LinSup_Inputs),            INTENT(INOUT) :: VBRDF_LinSup_In
    TYPE(VBRDF_Input_Exception_Handling), INTENT(INOUT) :: VBRDF_Sup_InputStatus
    TYPE(ErrorType),                      INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------

    ! =====================================================================
    ! InitializeBRDFInputs starts here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! These inputs are just from the old VBRDF supplement the GEOCAPE-TOOL was using
    ! At some point they should be updated 
    
    ! Control Inputs
    VBRDF_Sup_In%BS_DO_USER_STREAMS     = .TRUE.
    VBRDF_Sup_In%BS_DO_BRDF_SURFACE     = .TRUE.
    VBRDF_Sup_In%BS_DO_SURFACE_EMISSION = .FALSE.
    VBRDF_Sup_In%BS_DO_SOLAR_SOURCES    = .TRUE.
    VBRDF_Sup_In%BS_DO_USER_OBSGEOMS    = .FALSE.

    ! Geometry Results
    VBRDF_Sup_In%BS_NSTOKES              = Surface%nstokes
    VBRDF_Sup_In%BS_NSTREAMS             = Surface%nstreams

    ! BRDF Inputs
    VBRDF_Sup_In%BS_N_BRDF_KERNELS              = 3
    VBRDF_Sup_In%BS_BRDF_NAMES(1:3)             = (/'Lambertian','Ross-thick','Li-sparse '/)
    VBRDF_Sup_In%BS_WHICH_BRDF(1:3)             = (/1,3,4/)
    VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1:3)      = (/0,0,2/)
    VBRDF_Sup_In%BS_BRDF_PARAMETERS(:,:)        = 0.0
    VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG(1:3) = (/.TRUE.,.FALSE.,.FALSE./)
    VBRDF_Sup_In%BS_BRDF_FACTORS(:)             = 1.0d-5
    VBRDF_Sup_In%BS_NSTREAMS_BRDF               = 100 ! Azimuth angles
    VBRDF_Sup_In%BS_DO_SHADOW_EFFECT            = .FALSE.
    VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY        = .FALSE.

    VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR           = .FALSE.
    VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY    = .FALSE.
    VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER        = 0
    VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD      = 0
    VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD     = 0

    VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT = .FALSE.
    VBRDF_Sup_In%BS_DO_WSA_SCALING = .FALSE.
    VBRDF_Sup_In%BS_DO_BSA_SCALING = .FALSE.
    VBRDF_Sup_In%BS_WSA_VALUE      = 0.0
    VBRDF_Sup_In%BS_BSA_VALUE      = 0.0


    ! NewCM options
    VBRDF_Sup_In%BS_DO_NewCMGLINT = .FALSE.
    VBRDF_Sup_In%BS_SALINITY      = 0.0d0
    VBRDF_Sup_In%BS_WAVELENGTH    = 0.0d0
    VBRDF_Sup_In%BS_WINDSPEED     = 0.0d0
    VBRDF_Sup_In%BS_WINDDIR       = 0.0d0
    VBRDF_Sup_In%BS_DO_GlintShadow   = .FALSE.
    VBRDF_Sup_In%BS_DO_FoamOption    = .FALSE.
    VBRDF_Sup_In%BS_DO_FacetIsotropy = .FALSE.
!     
!     ! Linearizations
    VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(:) = .FALSE.
    VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS    = .FALSE.
    VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(:)    = .FALSE.
    VBRDF_LinSup_In%BS_N_SURFACE_WFS           = 0
    VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS     = 0
    VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS     = 0
    VBRDF_LinSup_In%BS_DO_WSAVALUE_WF          = .FALSE.         ! New Version 2.7
    VBRDF_LinSup_In%BS_DO_BSAVALUE_WF          = .FALSE.        ! New Version 2.7
    VBRDF_LinSup_In%BS_DO_WINDSPEED_WF         = .FALSE.       ! New Version 2.7

    ! Exception handling
    VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD              = 0
    VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES                = 0
    VBRDF_Sup_InputStatus%BS_INPUTMESSAGES(:)   = ' '
    VBRDF_Sup_InputStatus%BS_INPUTACTIONS(:)    = ' '
    VBRDF_Sup_InputStatus%BS_INPUTMESSAGES(0)              = 'Successful Read of VLIDORT Input file'
    VBRDF_Sup_InputStatus%BS_INPUTACTIONS                  = 'No Action required for this Task'

  END SUBROUTINE InitializeBRDFInputs
  
  SUBROUTINE SetSurfaceLinearization(self, DoAmp, DoPar, Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(SurfaceType),INTENT(INOUT) :: self
    LOGICAL,           INTENT(IN)    :: DoAmp
    LOGICAL,           INTENT(IN)    :: DoPar
    ! TYPE(SurfaceType), INTENT(INOUT) :: Surface
    TYPE(ErrorType),   INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER            :: I, J, ct
    CHARACTER(LEN=500) :: numstr

    ! =====================================================================
    ! SetSurfaceLinearization starts here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Save the linearization options
    self%DoAmplitudeLinearization = DoAmp
    self%DoParameterLinearization = DoPar
    self%nJac = 0
    
    ! Count the jacobians
    ! -------------------

    IF(self%DoLambertian) THEN

      IF( self%DoAmplitudeLinearization ) THEN

        ! 1 Kernel
        self%nJac = 1
        
        ! Check allocation
        IF( SIZE(self%JacName) .NE. self%nJac .OR. .NOT. ALLOCATED(self%JacName) ) THEN

          ! Allocate array to hold Jacobian Names
          IF(ALLOCATED(self%JacName)) DEALLOCATE(self%JacName)
          ALLOCATE(self%JacName(self%nJac))

          ! Kernel/Parameter indices
          IF(ALLOCATED(self%JacKernIdx)) DEALLOCATE(self%JacKernIdx)
          IF(ALLOCATED(self%JacParIdx)) DEALLOCATE(self%JacParIdx)
          ALLOCATE(self%JacKernIdx(self%nJac))
          ALLOCATE(self%JacParIdx(self%nJac))

        ENDIF
        
        ! Name the kernel
        self%JacName(1) = 'f_' // short_kern_name(1)
        self%JacKernIdx(1) =  1
        self%JacParIdx(1)  = -1
        
      ENDIF

    ELSE
      
      ! Kernel Amplitudes (factors)
      IF(self%DoAmplitudeLinearization) THEN
        self%nJac = self%nJac + self%nkern
      ENDIF

      ! Kernel parameters
      IF(self%DoParameterLinearization) THEN
        self%nJac = self%nJac + SUM(self%npar(1:self%nkern))
      ENDIF
      

      ! Check allocation
      IF( SIZE(self%JacName) .NE. self%nJac .OR. .NOT. ALLOCATED(self%JacName) ) THEN

        ! Allocate array to hold Jacobian Names
        IF(ALLOCATED(self%JacName)) DEALLOCATE(self%JacName)
        ALLOCATE(self%JacName(self%nJac))

        ! Kernel/Parameter indices
        IF(ALLOCATED(self%JacKernIdx)) DEALLOCATE(self%JacKernIdx)
        IF(ALLOCATED(self%JacParIdx)) DEALLOCATE(self%JacParIdx)
        ALLOCATE(self%JacKernIdx(self%nJac))
        ALLOCATE(self%JacParIdx(self%nJac))

      ENDIF

      ! Loop over kernels in order from VLIDORT to get kernel names
      ct = 0
      DO I=1,self%nkern

        ! Add factor jacobian
        IF( self%DoAmplitudeLinearization ) THEN
          ct = ct + 1
          self%JacName(ct)    = 'f_' // short_kern_name(self%kern_idx(I))
          self%JacKernIdx(ct) =  I
          self%JacParIdx(ct)  = -1
        ENDIF

        ! Add parameter jacobians
        IF(self%DoParameterLinearization) THEN

          ! Set parameter jacobians 
          VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_PARAMS_WFS(I,1:self%npar(I)) = .TRUE.
          VBRDF_LinSup_In_surfmod%BS_DO_KPARAMS_DERIVS(I) = .TRUE.

          ! Get Names
          DO J=1,self%npar(I)
            ct = ct + 1

            ! Write string
            WRITE(numstr,'(I500)') j

            self%JacName(ct) = 'p' // TRIM(ADJUSTL(numstr)) // '_' &
                              // short_kern_name(self%kern_idx(I))
            self%JacKernIdx(ct) =  I
            self%JacParIdx(ct)  =  J

          ENDDO

        ENDIF

      ENDDO

    ENDIF

  END SUBROUTINE SetSurfaceLinearization

  
  SUBROUTINE ComputeBRDFKernels(Wavelength, SZA, VZA, AZA, Surface, Error, SurfAux, nstreams, nstokes,DoIndividualKern)

     USE VBRDF_LINSUP_MASTERS_M

    ! --------------------
    ! subroutine arguments
    ! ---------------------
    REAL(KIND=8),                 INTENT(IN)    :: Wavelength
    REAL(KIND=8),                 INTENT(IN)    :: SZA
    REAL(KIND=8),                 INTENT(IN)    :: VZA
    REAL(KIND=8),                 INTENT(IN)    :: AZA
    TYPE(SurfaceType),            INTENT(INOUT) :: Surface
    TYPE(ErrorType),              INTENT(INOUT) :: Error
    TYPE(SurfProfType), OPTIONAL, INTENT(IN)    :: SurfAux
    INTEGER,            OPTIONAL, INTENT(IN)    :: nstreams, nstokes
    LOGICAL,            OPTIONAL, INTENT(IN)    :: DoIndividualKern

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER            :: i,j,k
    LOGICAL            :: DO_DEBUG_RESTORATION
    INTEGER            :: BS_NMOMENTS_INPUT
    CHARACTER(LEN=500) :: numstr
    LOGICAL            :: DoIndividualKern_in
    
    ! =====================================================================
    ! ComputeBRDFKernels starts here
    ! =====================================================================
    
    DoIndividualKern_in = .FALSE.
    IF(PRESENT(DoIndividualKern)) DoIndividualKern_in = DoIndividualKern

    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Intialize BRDF supplement inputs
    CALL InitializeBRDFInputs( Surface, VBRDF_Sup_In_surfmod, &
                               VBRDF_LinSup_In_surfmod,       &
                               VBRDF_Sup_InputStatus_surfmod, &
                               Error                          )
    
    ! Overwrite nstokes/nstreams (needed for two stream)
    IF(PRESENT(nstokes))  VBRDF_Sup_In_surfmod%BS_NSTOKES = nstokes
    IF(PRESENT(nstreams)) VBRDF_Sup_In_surfmod%BS_NSTREAMS = nstreams
    
    ! Set Geometry
    VBRDF_Sup_In_surfmod%BS_NBEAMS               = 1 ! Single geometry observations for now
    VBRDF_Sup_In_surfmod%BS_N_USER_STREAMS       = 1
    VBRDF_Sup_In_surfmod%BS_N_USER_RELAZMS       = 1
    VBRDF_Sup_In_surfmod%BS_BEAM_SZAS(1)         = SZA
    VBRDF_Sup_In_surfmod%BS_USER_ANGLES_INPUT(1) = VZA
    VBRDF_Sup_In_surfmod%BS_USER_RELAZMS(1)      = AZA
    
    VBRDF_Sup_In_surfmod%BS_WAVELENGTH = Wavelength*1.0d-3 ! supplement takes wavelength in microns
    
    ! Default glint options
    VBRDF_Sup_In_surfmod%BS_DO_GlintShadow   = .TRUE.
    VBRDF_Sup_In_surfmod%BS_DO_FoamOption    = .TRUE.
    VBRDF_Sup_In_surfmod%BS_DO_FacetIsotropy = .FALSE.
    
    ! Set additional surface properties if passed
    IF(PRESENT(SurfAux)) THEN
      VBRDF_Sup_In_surfmod%BS_SALINITY         = SurfAux%OceanSalinity
      VBRDF_Sup_In_surfmod%BS_WINDSPEED        = SurfAux%WindSpeed
      VBRDF_Sup_In_surfmod%BS_WINDDIR(1)       = SurfAux%WindDirection
    ENDIF
    
    ! Set number of kernels
    VBRDF_Sup_In_surfmod%BS_N_BRDF_KERNELS = Surface%nkern

    ! Initialize Lambertian kernel flag
    VBRDF_Sup_In_surfmod%BS_LAMBERTIAN_KERNEL_FLAG(:) = .FALSE.

    ! Set VBRDF supplement inputs
    DO I=1,Surface%nkern
      VBRDF_Sup_In_surfmod%BS_BRDF_NAMES(I)             = TRIM(ADJUSTL(Surface%kern_name(I)))
      VBRDF_Sup_In_surfmod%BS_WHICH_BRDF(I)             = Surface%kern_idx(I)
      VBRDF_Sup_In_surfmod%BS_N_BRDF_PARAMETERS(I)      = Surface%npar(I)
      IF(Surface%kern_idx(I) .EQ. 1) VBRDF_Sup_In_surfmod%BS_LAMBERTIAN_KERNEL_FLAG(I) = .TRUE.
      VBRDF_Sup_In_surfmod%BS_BRDF_FACTORS(I)           = Surface%WavelengthKernelAmplitude(I)
      DO J=1,Surface%npar(I)
        VBRDF_Sup_In_surfmod%BS_BRDF_PARAMETERS(I,J) = Surface%WavelengthKernelParameters(J,I) 
      ENDDO
    ENDDO
    
    ! Do BRDF surface
    VBRDF_Sup_In_surfmod%BS_DO_BRDF_SURFACE = .TRUE.

    ! Do not want debug restoration
    DO_DEBUG_RESTORATION = .false.

    ! A normal calculation will require
    BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In_surfmod%BS_NSTREAMS - 1
    
    ! Zero linearizations
    VBRDF_LinSup_In_surfmod%BS_N_SURFACE_WFS       = 0
    VBRDF_LinSup_In_surfmod%BS_N_KERNEL_FACTOR_WFS = 0 
    VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_FACTOR_WFS(:) = .FALSE.
    
    IF(DoIndividualKern_in) THEN

      ! Factor = 1.0d0
      VBRDF_Sup_In_surfmod%BS_BRDF_FACTORS(:) = 0.0d0
      VBRDF_Sup_In_surfmod%BS_BRDF_FACTORS(1) = 1.0d0

      ! Kernel Amplitudes (factors)
      IF(Surface%DoAmplitudeLinearization) THEN

        ! Update Jacobian number
        VBRDF_LinSup_In_surfmod%BS_N_KERNEL_FACTOR_WFS = 1

        ! Turn kernel factor linearization on
        VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.

      ENDIF

      ! Loop over kernels
      DO k=1,Surface%nkern

        ! Zero Parameter derivatives
        VBRDF_LinSup_In_surfmod%BS_N_KERNEL_PARAMS_WFS = 0
        VBRDF_LinSup_In_surfmod%BS_DO_KPARAMS_DERIVS(:) = .FALSE.
        VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_PARAMS_WFS(:,:) = .FALSE.

        ! Kernel parameters
        IF(Surface%DoParameterLinearization) THEN
          
          ! Count number of parameters
          VBRDF_LinSup_In_surfmod%BS_N_KERNEL_PARAMS_WFS = Surface%npar(k)
          
          ! Set parameter jacobians 
          VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_PARAMS_WFS(1,1:Surface%npar(I)) = .TRUE.
          VBRDF_LinSup_In_surfmod%BS_DO_KPARAMS_DERIVS(I) = .TRUE.

        ENDIF

        ! Total Number of kernels
        VBRDF_LinSup_In_surfmod%BS_N_SURFACE_WFS =              &
                VBRDF_LinSup_In_surfmod%BS_N_KERNEL_FACTOR_WFS  &
              + VBRDF_LinSup_In_surfmod%BS_N_KERNEL_PARAMS_WFS

        ! Call supplement to compute BRDF
        CALL VBRDF_LIN_MAINMASTER (           &
              DO_DEBUG_RESTORATION,           & ! Inputs
              BS_NMOMENTS_INPUT,              & ! Inputs
              VBRDF_Sup_In_surfmod,           & ! Inputs
              VBRDF_LinSup_In_surfmod,        & ! Inputs
              VBRDF_SupOut_surfmod,           & ! Outputs
              VBRDF_LinSupOut_surfmod,        & ! Outputs
              VBRDF_Sup_OutputStatus_surfmod  ) ! Output Status

        ! Store values in arrays
        Surface%VBRDF_Kernels(K)%BRDF_F_0        = VBRDF_SupOut_surfmod%BS_BRDF_F_0
        Surface%VBRDF_Kernels(K)%BRDF_F          = VBRDF_SupOut_surfmod%BS_BRDF_F
        Surface%VBRDF_Kernels(K)%USER_BRDF_F_0   = VBRDF_SupOut_surfmod%BS_USER_BRDF_F_0
        Surface%VBRDF_Kernels(K)%USER_BRDF_F     = VBRDF_SupOut_surfmod%BS_USER_BRDF_F
        Surface%VBRDF_Kernels(K)%DBOUNCE_BRDFUNC = VBRDF_SupOut_surfmod%BS_DBOUNCE_BRDFUNC
        Surface%VBRDF_Kernels(K)%EMISSIVITY      = VBRDF_SupOut_surfmod%BS_EMISSIVITY
        Surface%VBRDF_Kernels(K)%USER_EMISSIVITY = VBRDF_SupOut_surfmod%BS_USER_EMISSIVITY

        Surface%VBRDF_Kernels(K)%LS_BRDF_F_0        = VBRDF_LinSupOut_surfmod%BS_LS_BRDF_F_0
        Surface%VBRDF_Kernels(K)%LS_BRDF_F          = VBRDF_LinSupOut_surfmod%BS_LS_BRDF_F
        Surface%VBRDF_Kernels(K)%LS_USER_BRDF_F_0   = VBRDF_LinSupOut_surfmod%BS_LS_USER_BRDF_F_0
        Surface%VBRDF_Kernels(K)%LS_USER_BRDF_F     = VBRDF_LinSupOut_surfmod%BS_LS_USER_BRDF_F
        Surface%VBRDF_Kernels(K)%LS_DBOUNCE_BRDFUNC = VBRDF_LinSupOut_surfmod%BS_LS_DBOUNCE_BRDFUNC
        Surface%VBRDF_Kernels(K)%LS_USER_EMISSIVITY = VBRDF_LinSupOut_surfmod%BS_LS_USER_EMISSIVITY
        Surface%VBRDF_Kernels(K)%LS_EMISSIVITY      = VBRDF_LinSupOut_surfmod%BS_LS_EMISSIVITY

      ENDDO

      
      



      

    ELSE
      
      ! Kernel Amplitudes (factors)
      IF(Surface%DoAmplitudeLinearization) THEN

        ! Update Jacobian number
        VBRDF_LinSup_In_surfmod%BS_N_KERNEL_FACTOR_WFS = Surface%nkern

        ! Turn kernel factor linearization on
        VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_FACTOR_WFS(1:Surface%nkern) = .TRUE.

      ENDIF

      ! Kernel parameters
      IF(Surface%DoParameterLinearization) THEN
        
        ! Count number of parameters
        VBRDF_LinSup_In_surfmod%BS_N_KERNEL_PARAMS_WFS = SUM(Surface%npar(1:Surface%nkern))
        
      ENDIF

      ! Total number of weighting functions
      VBRDF_LinSup_In_surfmod%BS_N_SURFACE_WFS = Surface%nJac

      ! Loop over kernels in order from VLIDORT to get kernel names
      IF(Surface%DoParameterLinearization) THEN
        DO I=1,Surface%nkern
          ! Set parameter jacobians 
          VBRDF_LinSup_In_surfmod%BS_DO_KERNEL_PARAMS_WFS(I,1:Surface%npar(I)) = .TRUE.
          VBRDF_LinSup_In_surfmod%BS_DO_KPARAMS_DERIVS(I) = .TRUE.
        ENDDO
      ENDIF
      
      ! Call supplement to compute BRDF
      CALL VBRDF_LIN_MAINMASTER (           &
            DO_DEBUG_RESTORATION,           & ! Inputs
            BS_NMOMENTS_INPUT,              & ! Inputs
            VBRDF_Sup_In_surfmod,           & ! Inputs
            VBRDF_LinSup_In_surfmod,        & ! Inputs
            VBRDF_SupOut_surfmod,           & ! Outputs
            VBRDF_LinSupOut_surfmod,        & ! Outputs
            VBRDF_Sup_OutputStatus_surfmod  ) ! Output Status

      !  Exception handling
      IF ( VBRDF_Sup_OutputStatus_surfmod%BS_STATUS_OUTPUT .EQ. VLIDORT_SERIOUS ) THEN
        WRITE(*,*)'SPLAT: program failed, VBRDF calculation aborted'
        WRITE(*,*)'Here are the error messages from the VBRDF supplement : - '
        WRITE(*,*)' - Number of error messages = ',VBRDF_Sup_OutputStatus_surfmod%BS_NOUTPUTMESSAGES
        DO I = 1, VBRDF_Sup_OutputStatus_surfmod%BS_NOUTPUTMESSAGES
            WRITE(*,*) '  * ',adjustl(Trim(VBRDF_Sup_OutputStatus_surfmod%BS_OUTPUTMESSAGES(I)))
        ENDDO
      ENDIF
      
      ! Transfer to VLIDORT Input Arrays
      Surface%VLBRDF%BRDF_F_0        = VBRDF_SupOut_surfmod%BS_BRDF_F_0
      Surface%VLBRDF%BRDF_F          = VBRDF_SupOut_surfmod%BS_BRDF_F
      Surface%VLBRDF%USER_BRDF_F_0   = VBRDF_SupOut_surfmod%BS_USER_BRDF_F_0
      Surface%VLBRDF%USER_BRDF_F     = VBRDF_SupOut_surfmod%BS_USER_BRDF_F
      Surface%VLBRDF%DBOUNCE_BRDFUNC = VBRDF_SupOut_surfmod%BS_DBOUNCE_BRDFUNC
      Surface%VLBRDF%EMISSIVITY      = VBRDF_SupOut_surfmod%BS_EMISSIVITY
      Surface%VLBRDF%USER_EMISSIVITY = VBRDF_SupOut_surfmod%BS_USER_EMISSIVITY

      Surface%VLBRDF%LS_BRDF_F_0        = VBRDF_LinSupOut_surfmod%BS_LS_BRDF_F_0
      Surface%VLBRDF%LS_BRDF_F          = VBRDF_LinSupOut_surfmod%BS_LS_BRDF_F
      Surface%VLBRDF%LS_USER_BRDF_F_0   = VBRDF_LinSupOut_surfmod%BS_LS_USER_BRDF_F_0
      Surface%VLBRDF%LS_USER_BRDF_F     = VBRDF_LinSupOut_surfmod%BS_LS_USER_BRDF_F
      Surface%VLBRDF%LS_DBOUNCE_BRDFUNC = VBRDF_LinSupOut_surfmod%BS_LS_DBOUNCE_BRDFUNC
      Surface%VLBRDF%LS_USER_EMISSIVITY = VBRDF_LinSupOut_surfmod%BS_LS_USER_EMISSIVITY
      Surface%VLBRDF%LS_EMISSIVITY      = VBRDF_LinSupOut_surfmod%BS_LS_EMISSIVITY
    ENDIF

  END SUBROUTINE ComputeBRDFKernels

END MODULE surface_module
