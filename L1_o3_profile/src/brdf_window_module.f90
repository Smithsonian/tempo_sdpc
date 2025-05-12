MODULE window_module
  
  USE parameters_module
  USE error_module,         ONLY : ErrorType, RaiseFatalError, CheckError
  USE level1_def,           ONLY : L1Type
  USE interpolation_module, ONLY : SPLINE1, SPLINT1, LinearInterpolationWeights,&
                                   BSPLINE_EdgeFill
  USE rtm_def,              ONLY : WinRTMSettingsType, RTMOptType
  USE netcdf_module,        ONLY : CheckNetCDFErrorStatus
  USE pca_bin_module,       ONLY : PCABinType, CreatePCABins
  USE isrf_module,          ONLY : ISRF_FunctionType
  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  
  TYPE ISRFOptType
    CHARACTER(LEN=maxChar) :: Name
    LOGICAL                :: PixDepPar
    CHARACTER(LEN=maxChar) :: Infile
    INTEGER                :: Infile_BandIndex
    INTEGER                :: OptionIndex
    
    ! For reading from the input file
    INTEGER                   :: nFixedPar
    REAL(KIND=8), ALLOCATABLE :: FixedPar(:)
  END TYPE ISRFOptType

  ! Input Options
  TYPE WinOptType
    INTEGER                :: BandIndex
    REAL(KIND=8)           :: StartWvl
    REAL(KIND=8)           :: EndWvl
    REAL(KIND=8)           :: BufferWvl
    CHARACTER(LEN=maxChar) :: ForwardModelName
    CHARACTER(LEN=maxChar) :: RefSolarInfile
    LOGICAL                :: DO_RTM_AT_L1
    LOGICAL                :: DoSolarI0Correction
    REAL(KIND=8)           :: Convol_dWvl
    TYPE(ISRFOptType)      :: ISRF
    REAL(KIND=8)           :: ConvolWidth_HW1E
  ENDTYPE WinOptType

  TYPE WinType
    TYPE(WinOptType)          :: Settings
    TYPE(WinRTMSettingsType)  :: RTMSettings
    CHARACTER(LEN=maxChar)    :: CalculationMode
    
    ! Radiance Output Unit
    CHARACTER(LEN=maxChar)    :: RadianceUnit

    ! Grid for RTM Calculations
    INTEGER                   :: nRTM_wvl
    REAL(KIND=8), ALLOCATABLE :: RTM_wvl(:)
    INTEGER                   :: nConvol_Wvl
    REAL(KIND=8), ALLOCATABLE :: Convol_Wvl(:)

    ! Convolved IO Reference
    REAL(KIND=8)              :: I0_ScaleFactor
    REAL(KIND=8), ALLOCATABLE :: I0_ConvolGrid(:)
    REAL(KIND=8), ALLOCATABLE :: I0_ConvolGrid_ISRF(:)
    REAL(KIND=8), ALLOCATABLE :: I0_ConvolGrid_ISRF_SP(:)
    REAL(KIND=8), ALLOCATABLE :: I0_L1Grid_ISRF(:)

    ! ISRF
    LOGICAL                   :: UseInfileISRF
    REAL(KIND=8), ALLOCATABLE :: ISRFPar_ConvolGrid(:,:)
    REAL(KIND=8), ALLOCATABLE :: ISRFPar_L1Grid(:,:)
    REAL(KIND=8), ALLOCATABLE :: ISRFPar_Fixed(:)
    REAL(KIND=8)              :: ISRF_hw1e_Fixed
    INTEGER                   :: nT_lut, nw_lut, nx_lut, np_lut
    REAL(KIND=8), ALLOCATABLE :: ISRF_LUT_Par(:,:,:,:)
    REAL(KIND=8), ALLOCATABLE :: ISRF_LUT_T(:)
    LOGICAL                   :: ISRF_LUT_Tdep
    LOGICAL                   :: ISRF_LUT_wdep
    REAL(KIND=8), ALLOCATABLE :: ISRF_hw1e_ConvolGrid(:)
    REAL(KIND=8), ALLOCATABLE :: ISRF_hw1e_L1Grid(:) ! To figure out convolution window
    TYPE(ISRF_FunctionType)   :: ISRF_Function

    ! Level1 Grid within window
    INTEGER                   :: nL1_Wvl
    INTEGER                   :: L1_iw0
    INTEGER                   :: L1_iwf
    REAL(KIND=8), ALLOCATABLE :: L1_Wvl(:)

    ! Solar Reference
    INTEGER                   :: nRefWvl
    REAL(KIND=8), ALLOCATABLE :: RefSolarWavelength(:)
    REAL(KIND=8), ALLOCATABLE :: RefSolarIrradiance(:)
    REAL(KIND=8), ALLOCATABLE :: RefSolarIrradianceSP(:)

    ! PCA Binning Strategy for window
    TYPE(PCABinType)          :: PCABin

  ENDTYPE WinType
  
  PUBLIC :: InitWindow, &
            SetWindow, &
            ApplyISRF,&
            GetWindowNumStokes,&
            GetWindowNumStreams,&
            GetWindowNumMoments
            
  ! PRIVATE :: super_gaussian_sf, &
  !            signdp

  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'window_module'
  PRIVATE :: ModuleName

  CONTAINS
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: InitWindow
  ! 
  ! DESCRIPTION: Initializes a window for simulation/retrieval
  !                - Sets convolution grid
  !                - Initializes ISRF 
  !                - Prepares Solar irradiance on convolution grid

  SUBROUTINE InitWindow( WinOpt, RTMOpt, CalculationMode, Window, Error, AMFWavelength )
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(WinOptType),       INTENT(IN)    :: WinOpt
    TYPE(RTMOptType),       INTENT(IN)    :: RTMOpt
    CHARACTER(LEN=maxChar), INTENT(IN)    :: CalculationMode
    TYPE(WinType),          INTENT(INOUT) :: Window
    TYPE(ErrorType),        INTENT(INOUT) :: Error
    REAL(KIND=8), OPTIONAL, INTENT(IN)    :: AMFWavelength
    
    ! ---------------
    ! Local Variables
    ! ---------------
    REAL(KIND=8)              :: wvl_0,wvl_f
    INTEGER                   :: w, ncid, rcode, vid, dimid
    REAL(KIND=8), ALLOCATABLE :: wvl_in(:), sol_in(:)
    CHARACTER(LEN=maxChar)    :: tmpchar
    INTEGER                   :: nwvl, n_below, i, n, nISRFPar

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitWindow'

    ! =====================================================================
    ! InitWindow starts here 
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Store Input Settings
    Window%Settings = WinOpt ; Window%CalculationMode = CalculationMode
    
    ! Hardwire Convolution width for now
    Window%Settings%ConvolWidth_HW1E =10.0d0

    ! Select the RTM/Save settings for this window
    Window%RTMSettings%RTMName = ''
    DO N=1,RTMOpt%nSimType
      
      IF( TRIM(ADJUSTL(RTMOpt%SimName(N))) .EQ.  &
          TRIM(ADJUSTL(WinOpt%ForwardModelName)) ) THEN
        
        ! Store the settings for the window
        Window%RTMSettings%RTMName             = RTMOpt%RTMName(N)
        Window%RTMSettings%SimOpt              = RTMOpt%SimOpt(N)
        Window%RTMSettings%Jacobian            = RTMOpt%Jacobian
        Window%RTMSettings%DoScatteringWeights = RTMOpt%DoScatteringWeights
        Window%RTMSettings%DoFluxcalculation   = RTMOpt%DoFluxCalculation
        
        ! Exit loop
        EXIT
        
      ENDIF
      
    ENDDO
    
    ! Check if simulation has been matched
    IF(TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. '') THEN
      print*,TRIM(ADJUSTL(Window%RTMSettings%RTMName)) // &
             ' could not be matched to any RTM sim settings'
      STOP 1
    ENDIF
    
    ! Determine if we need to load the ISRF from file
    ! The input module flags this by setting OptionIndex = 0
    ! OptionIndex gets set after loading file
    Window%UseInfileISRF = .FALSE.
    IF(Window%Settings%ISRF%OptionIndex .LE. 0) THEN

      ! Load the ISRF
      CALL LoadInputFileISRF(Window,Error)

      ! Set flag to use input file parameter ISRF
      !  0 => LUT ISRF 
      ! >0 => Table of ISRF Parameters
      Window%UseInfileISRF = .TRUE.

    ELSE 

      ! Transfer Options
      ! CCM-Fix - Probably don't need two copies of this...
      nISRFPar = Window%Settings%ISRF%nFixedPar
      ALLOCATE(Window%ISRFPar_Fixed(nISRFPar))
      Window%ISRFPar_Fixed(1:nISRFPar) = Window%Settings%ISRF%FixedPar

    ENDIF
    
    ! Initialize the ISRF parameterization
    IF(Window%Settings%ISRF%OptionIndex .EQ. 0) THEN
      CALL Window%ISRF_Function%Initialize(Window%Settings%ISRF%OptionIndex,Error,   &
                                    LookupTableFile=Window%Settings%ISRF%Infile,     &
                                    BandIndex=Window%Settings%ISRF%Infile_BandIndex, &
                                    Convol_dWvl=Window%Settings%Convol_dWvl,         &
                                    ConvolWidth_HW1E=Window%Settings%ConvolWidth_HW1E)
      
      ! In this case ISRF is pixel dependent - Allocate 
      
    ELSE
      CALL Window%ISRF_Function%Initialize(Window%Settings%ISRF%OptionIndex,Error)
    ENDIF

    IF(TRIM(ADJUSTL(CalculationMode)) .EQ. 'AMF' ) THEN

      ! Check if AMF wavelength present
      IF(.NOT. PRESENT(AMFWavelength)) THEN
        print*,'Calculation is in AMF Mode...'
        print*,'ERROR:Could not initialize spectral window as AMFWavelength was not provided'
        STOP 1
      ENDIF

      ! Define window
      Window%Settings%StartWvl = AMFWavelength ; Window%Settings%EndWvl = AMFWavelength
      
      IF(Window%Settings%ISRF%OptionIndex .EQ. 0) THEN
        print*,'Cannot use instrument LUT ISRF in AMF Calculation mode'
        print*,'Use one of the fixed options in the input file'
        STOP 1
      ENDIF
    ENDIF

    ! Add 2x buffer to the wavelengths (Force overlap with RTM Grid)
    wvl_0 = Window%Settings%StartWvl - Window%Settings%BufferWvl
    wvl_f = Window%Settings%EndWvl   + Window%Settings%BufferWvl
      
    ! Optical Property grid - Number of wavelengths
    Window%nConvol_Wvl = NINT ((wvl_f - wvl_0)/Window%Settings%Convol_dWvl) + 1
    
    ! Allocate Grid
    ALLOCATE(Window%Convol_Wvl(Window%nConvol_Wvl))
    ALLOCATE(Window%I0_ConvolGrid(Window%nConvol_Wvl))
    ALLOCATE(Window%I0_ConvolGrid_ISRF(Window%nConvol_Wvl))
    ALLOCATE(Window%I0_ConvolGrid_ISRF_SP(Window%nConvol_Wvl))
    
    ! Compute Grid
    Window%Convol_Wvl(1) = wvl_0
    DO w=2,Window%nConvol_Wvl
      Window%Convol_Wvl(w) = Window%Convol_Wvl(w-1) + Window%Settings%Convol_dWvl
    ENDDO
    
    ! Wavelength is constant - Can set the window on initialization
    IF(TRIM(ADJUSTL(CalculationMode)) .EQ. 'AMF' ) THEN

      ! Define "L1 Grid" as the AMF Wavelength
      Window%L1_iw0  = 1 ; Window%L1_iwf = 1 ; Window%nL1_Wvl = 1
      ALLOCATE(Window%L1_Wvl(1)) ; Window%L1_Wvl(1) = AMFWavelength

      ! Set up RTM Grid
      IF( Window%Settings%DO_RTM_AT_L1 ) THEN
        
        ! RTM = L1 Grid 
        Window%nRTM_wvl = 1
        ALLOCATE(Window%RTM_wvl(1)) ; Window%RTM_wvl(1) = AMFWavelength
        
      ELSE

        ! The RTM grid -> Convolution Grid
        Window%nRTM_Wvl = Window%nConvol_Wvl
        ALLOCATE(Window%RTM_Wvl(Window%nRTM_Wvl))
        Window%RTM_Wvl = Window%Convol_Wvl(:)
        
      ENDIF
      
    ELSE

      ! Set up RTM Grid
      IF( Window%Settings%DO_RTM_AT_L1 ) THEN
        
        ! RTM Must Be set to L1 Grid (Varies by position)
        Window%nRTM_wvl = 0
        
      ELSE

        ! The RTM grid -> Convolution Grid
        Window%nRTM_Wvl = Window%nConvol_Wvl
        ALLOCATE(Window%RTM_Wvl(Window%nRTM_Wvl))
        Window%RTM_Wvl = Window%Convol_Wvl(:)
        
      ENDIF

    ENDIF

    ! Additional Flags for Fixed ISRFs
    IF(.NOT. Window%Settings%ISRF%PixDepPar) THEN

      ! Flag Instances where we do not need to perform convolutions
      ! -----------------------------------------------------------

      ! Supergaussian case (width = 0.0)
      IF(Window%Settings%ISRF%OptionIndex .EQ. 1 .AND.   &
         Window%ISRFPar_Fixed(1) .LT. TINY(0.0d0)        ) THEN
          Window%ISRF_Function%IsDeltaFunction = .TRUE.
          Window%ISRF_hw1e_Fixed = 0.0d0
      ENDIF

      ! Determine Convolution hw1e
      IF(.NOT. Window%ISRF_Function%IsDeltaFunction) THEN
        CALL Window%ISRF_Function%DetermineHW1E(                        &
                      Window%ISRFPar_Fixed, Window%Settings%Convol_dWvl,&
                      Window%nConvol_Wvl, Window%ISRF_hw1e_Fixed        )
      ENDIF

    ENDIF

    ! Check for Pixel dependent ISRF
    IF(Window%Settings%ISRF%PixDepPar .OR. Window%UseInfileISRF) THEN

      ! Allocate ISRF array on convolution grid
      nISRFPar = Window%ISRF_Function%nPar
      ALLOCATE(Window%ISRFPar_ConvolGrid(Window%nConvol_Wvl,nISRFPar))
      ALLOCATE(Window%ISRF_hw1e_ConvolGrid(Window%nConvol_Wvl))

    ENDIF

    ! -------------------------------------
    ! Load Solar Reference
    ! -------------------------------------

    ! Open file
    rcode = nf_open(TRIM(ADJUSTL(Window%Settings%RefSolarInfile)), NF_SHARE, ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
             'nf_open:'//TRIM(ADJUSTL(Window%Settings%RefSolarInfile)))

    ! Find the dimension of the spectrum
    rcode = nf_inq_varid( ncid, 'Wavelength', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_var:solar_wvl')
    rcode = nf_inq_vardimid(ncid,vid,dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_vardim:solar_wvl')
    rcode = NF_INQ_DIM(ncid,dimid,tmpchar, nwvl)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_dim:solar_wvl')
    
    ! read the whole wavelength array
    ALLOCATE(wvl_in(nwvl))
    rcode = nf_get_var_double(ncid, vid, wvl_in)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_get_var:solar_wvl')

    ! Read spectrum
    ALLOCATE( sol_in(nwvl) )
    rcode = nf_inq_varid(ncid, 'Irradiance', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_var:solar_irr')
    rcode = nf_get_var_double( ncid, vid, sol_in)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_get_var:solar_irr')
  
    ! Close file
    rcode = nf_close(ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_close:solar')

    ! Load wavelengths (use wide buffer in case HITRAN cross section option set)
    wvl_0 = Window%Convol_Wvl(1                 )-10.0d0*Window%Settings%BufferWvl
    wvl_f = Window%Convol_Wvl(Window%nConvol_Wvl)+10.0d0*Window%Settings%BufferWvl
    
    ! Determine indices within range
    n_below = 0
    Window%nRefWvl = 0
    DO i=1,nwvl
      IF( wvl_in(i) .LT. wvl_0 ) THEN
        n_below = n_below + 1
      ELSEIF( wvl_in(i) .GT. wvl_f ) THEN
        EXIT
      ELSE
        Window%nRefWvl = Window%nRefWvl + 1
      ENDIF 
    ENDDO

    ! Allocate output arrays
    ALLOCATE(Window%RefSolarWavelength(Window%nRefWvl))
    ALLOCATE(Window%RefSolarIrradiance(Window%nRefWvl))
    ALLOCATE(Window%RefSolarIrradianceSP(Window%nRefWvl))

    ! Get subset
    Window%RefSolarWavelength = wvl_in(n_below+1:n_below+Window%nRefWvl)
    Window%RefSolarIrradiance = sol_in(n_below+1:n_below+Window%nRefWvl)

    ! Compute basis spline coefficients
    CALL SPLINE1(Window%RefSolarWavelength,Window%RefSolarIrradiance,&
                 Window%nRefWvl,Window%RefSolarIrradianceSP          )
    
    ! Spline Solar reference to I0 Grid
    CALL SPLINT1 (Window%RefSolarWavelength,Window%RefSolarIrradiance,       &
                  Window%RefSolarIrradianceSP, Window%nRefWvl,               &
                  Window%Convol_Wvl, Window%I0_ConvolGrid, Window%nConvol_Wvl)
    
    ! Normalize I0 Spectrum to ~ 1
    Window%I0_ScaleFactor = PhotonScalingUnit
    Window%I0_ConvolGrid = Window%I0_ConvolGrid / Window%I0_ScaleFactor
    
    ! Set up PCA Binning if needed
    IF(TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-PCA' ) THEN

      CALL Window%PCABin%InitializePCABins(                             &
              Window%RTMSettings%SimOpt%VLIDORT_PCA%WinControlFile,Error)

    ENDIF
    
    ! Can preconvolve the solar spectrum if using the input file ISRF
    IF(.NOT. Window%Settings%ISRF%PixDepPar) THEN
      CALL ApplyISRF_HRSP(Window, Window%I0_ConvolGrid, Window%I0_ConvolGrid_ISRF,&
                          Window%I0_ConvolGrid_ISRF_SP, Error                     )
    ENDIF

  END SUBROUTINE InitWindow
  
  SUBROUTINE SetWindow( Level1, Window, Error )
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    TYPE(L1Type),    INTENT(IN)    :: Level1
    TYPE(WinType),   INTENT(INOUT) :: Window
    TYPE(ErrorType), INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: nwvl, b, w, errstat
    LOGICAL :: BelowWindow
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SetWindow'

    ! =====================================================================
    ! SetWindow starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    IF(Error%VerbosityLevel .GT. 1) print*,'         - Setting Window'

    ! Store Radiance Unit
    Window%RadianceUnit = Level1%Spectrum(Window%Settings%BandIndex)%RadianceUnit
    
    ! Store Cross Track (needed by ISRF)
    Window%ISRF_Function%xtrk_idx = Level1%ix
    
    ! We dont need to do anything for AMF Mode as window is set on initialization
    IF(TRIM(ADJUSTL(Window%CalculationMode)) .EQ. 'AMF' ) THEN

      ! Load ISRF On convolution grid if needed
      IF(Window%UseInfileISRF) THEN
        CALL GetInputFileISRF(Window,Level1,Error)
      ENDIF

    ELSE

      ! Scan to first index
      BelowWindow = .TRUE.
      
      ! Get band Index
      b = Window%Settings%BandIndex
      
      ! Get Wavelength Limits
      Window%L1_iw0  = MINLOC(Level1%Spectrum(b)%Wavelength, DIM=1, &
        MASK = Level1%Spectrum(b)%Wavelength  .GE. Window%Settings%StartWvl)
      Window%L1_iwf = MAXLOC(Level1%Spectrum(b)%Wavelength, DIM=1, &
                  MASK = Level1%Spectrum(b)%Wavelength  .LE. Window%Settings%EndWvl)
      Window%nL1_Wvl = Window%L1_iwf - Window%L1_iw0 + 1

      ! Save the L1 wavelength grid
      ! Reallocate L1 Wavelength array if necessary
      IF(ALLOCATED(Window%L1_Wvl)) DEALLOCATE(Window%L1_Wvl)
      ALLOCATE(Window%L1_Wvl(Window%nL1_Wvl))
      Window%L1_Wvl(:)=Level1%Spectrum(b)%Wavelength(Window%L1_iw0:Window%L1_iwf)

      ! Set RTM Grid if we are doing at native satellite resolution
      IF( Window%Settings%DO_RTM_AT_L1 ) THEN 
        
        ! Reallocate RTM Wavelength array if necessary
        IF(Window%nRTM_wvl .NE. Window%nL1_Wvl .OR. &
          .NOT. ALLOCATED(Window%RTM_wvl)          ) THEN
          IF(ALLOCATED( Window%RTM_wvl )) DEALLOCATE( Window%RTM_wvl )
          Window%nRTM_wvl = Window%nL1_Wvl
          ALLOCATE( Window%RTM_wvl(Window%nRTM_wvl ) )
        ENDIF

        ! Copy L1 Wavelength Grid
        Window%RTM_wvl(:) = Window%L1_Wvl(:)
        
      ENDIF
      
      ! Load ISRF (On Convolution or fixed Grids)
      IF(Window%UseInfileISRF) THEN
        CALL GetInputFileISRF(Window,Level1,Error)
      ELSEIF(Window%Settings%ISRF%OptionIndex .EQ. 0) THEN
        

      ENDIF

    ENDIF

    ! Determine HW1E For Convolution Limits
    IF(Window%Settings%ISRF%PixDepPar) THEN
      
      ! Reallocate HW1E array stored on the L1 Grid
      IF(ALLOCATED(Window%ISRF_hw1e_L1Grid)) DEALLOCATE(Window%ISRF_hw1e_L1Grid)
      ALLOCATE(Window%ISRF_hw1e_L1Grid(Window%nL1_Wvl))
      
      ! Compute HW1E on L1 Grid
      DO w=1,Window%nL1_Wvl
        CALL Window%ISRF_Function%DetermineHW1E(                           &
                    Window%ISRFPar_L1Grid(w,:),Window%Settings%Convol_dWvl,&
                    Window%nConvol_Wvl,Window%ISRF_hw1e_L1Grid(w)          )
      ENDDO

      ! Now on the convolution Grid
      DO w=1,Window%nConvol_Wvl
        CALL Window%ISRF_Function%DetermineHW1E(                                      &
                           Window%ISRFPar_ConvolGrid(w,:),Window%Settings%Convol_dWvl,&
                           Window%nConvol_Wvl,Window%ISRF_hw1e_ConvolGrid(w)          )
      ENDDO

    ELSE
      
      ! Compute HW1E using fixed parameters
      CALL Window%ISRF_Function%DetermineHW1E(                               &
                         Window%ISRFPar_Fixed,Window%Settings%Convol_dWvl,   &
                         Window%nConvol_Wvl,Window%ISRF_hw1e_Fixed           )
      
    ENDIF
    
    ! Update solar spectrum if not using variable ISRF parameters (can change spectrum to spectrum)
    IF(Window%Settings%ISRF%PixDepPar .AND. .NOT. Window%UseInfileISRF) THEN
      CALL ApplyISRF_HRSP(Window, Window%I0_ConvolGrid, Window%I0_ConvolGrid_ISRF,&
                          Window%I0_ConvolGrid_ISRF_SP, Error                     )
    ENDIF
    
    ! Interpolate to L1 Grid if doing an AMF Calculation
    IF(ALLOCATED(Window%I0_L1Grid_ISRF)) DEALLOCATE(Window%I0_L1Grid_ISRF)
    ALLOCATE(Window%I0_L1Grid_ISRF(Window%nL1_Wvl))
    CALL SPLINT1 (Window%Convol_Wvl,Window%I0_ConvolGrid_ISRF,        &
                  Window%I0_ConvolGrid_ISRF_SP, Window%nConvol_Wvl,   &
                  Window%L1_Wvl, Window%I0_L1Grid_ISRF, Window%nL1_Wvl)
    
    
  END SUBROUTINE SetWindow

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: GetWindowNumStokes
  ! 
  ! DESCRIPTION: Helper function to get the number of stokes elements 
  !              simulated by the RTM associated with the window

  INTEGER(KIND=4) FUNCTION GetWindowNumStokes(Window,Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    TYPE(ErrorType), INTENT(INOUT) :: Error
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'GetWindowNumStokes'

    ! =====================================================================
    ! GetWindowNumStokes starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Use the window RTM Settings
    IF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-LBL' ) THEN
      GetWindowNumStokes = Window%RTMSettings%SimOpt%VLIDORT%nstokes
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-PCA' ) THEN
      GetWindowNumStokes = Window%RTMSettings%SimOpt%VLIDORT_PCA%nstokes
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'FIRSTORD' ) THEN
      GetWindowNumStokes = Window%RTMSettings%SimOpt%FIRST_ORDER%nstokes
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. '2STRM' ) THEN
      GetWindowNumStokes = Window%RTMSettings%SimOpt%TWO_STREAM%nstokes
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'TRANSMISSION' ) THEN
      GetWindowNumStokes = 1
    ELSE
      print*,'In Window Model'
      print*,'Unrecognized Radiative Transfer Model'
      STOP 1
    ENDIF
    
  END FUNCTION GetWindowNumStokes
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: GetWindowNumMoments
  ! 
  ! DESCRIPTION: Helper function to get the number of phase func. 
  !              expansion coefficients used by the RTM that is 
  !              associated with the window

  INTEGER(KIND=4) FUNCTION GetWindowNumMoments(Window,Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    TYPE(ErrorType), INTENT(INOUT) :: Error
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'GetWindowNumMoments'

    ! =====================================================================
    ! GetWindowNumMoments starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Use the window RTM Settings
    IF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-LBL' ) THEN
      GetWindowNumMoments = Window%RTMSettings%SimOpt%VLIDORT%nmoments
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-PCA' ) THEN
      GetWindowNumMoments = Window%RTMSettings%SimOpt%VLIDORT_PCA%nmoments
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'FIRSTORD' ) THEN
      GetWindowNumMoments = Window%RTMSettings%SimOpt%FIRST_ORDER%nmoments
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. '2STRM' ) THEN
      GetWindowNumMoments = Window%RTMSettings%SimOpt%TWO_STREAM%nmoments
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'TRANSMISSION' ) THEN
      GetWindowNumMoments = 2
    ELSE
      print*,'In Window Model'
      print*,'Unrecognized Radiative Transfer Model'
      STOP 1
    ENDIF
    
  END FUNCTION GetWindowNumMoments
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: ForceISRFPixDepPar
  ! 
  ! DESCRIPTION: Makes sure that the pixel-dependent ISRF convolution 
  !              is used. This is needed when a pixel-dependent 
  !              ISRF parameter is optimized in INVERSE mode

  SUBROUTINE ForceISRFPixDepPar(Window,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(INOUT) :: Window
    TYPE(ErrorType), INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: P, nP

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ForceISRFPixDepPar'

    ! =====================================================================
    ! ForceISRFPixDepPar starts here 
    ! =====================================================================

    ! Check if we need to do anything
    IF(Window%Settings%ISRF%PixDepPar) RETURN

    ! Get the number of parameters
    nP = Window%ISRF_Function%nPar

    ! Allocate ISRF parameters
    IF(ALLOCATED(Window%ISRFPar_ConvolGrid)) DEALLOCATE(Window%ISRFPar_ConvolGrid)
    ALLOCATE(Window%ISRFPar_ConvolGrid(Window%nConvol_Wvl,nP))
    IF(ALLOCATED(Window%ISRFPar_L1Grid)) DEALLOCATE(Window%ISRFPar_L1Grid)
    ALLOCATE(Window%ISRFPar_L1Grid(Window%nL1_Wvl,nP))
    IF(ALLOCATED(Window%ISRF_hw1e_L1Grid)) DEALLOCATE(Window%ISRF_hw1e_L1Grid)
    ALLOCATE(Window%ISRF_hw1e_L1Grid(Window%nL1_Wvl))
    IF(ALLOCATED(Window%ISRF_hw1e_ConvolGrid)) DEALLOCATE(Window%ISRF_hw1e_ConvolGrid)
    ALLOCATE(Window%ISRF_hw1e_ConvolGrid(Window%nConvol_Wvl))

    ! Set the wvl-dep parameters from the fixed ones
    DO P=1,nP
      Window%ISRFPar_ConvolGrid(:,P) = Window%ISRFPar_Fixed(P)
      Window%ISRFPar_L1Grid(:,P) = Window%ISRFPar_Fixed(P)
    ENDDO

    ! Set hw1e
    Window%ISRF_hw1e_L1Grid(:)     = Window%ISRF_hw1e_Fixed
    Window%ISRF_hw1e_ConvolGrid(:) = Window%ISRF_hw1e_Fixed


  END SUBROUTINE ForceISRFPixDepPar

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: UpdateISRFParameter
  ! 
  ! DESCRIPTION: Updates the value(s) of an ISRF parameter for an
  !              instance of the window type

  SUBROUTINE UpdateISRFParameter(Window,ParIdx,nParVal,ParVal,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(INOUT) :: Window
    INTEGER,         INTENT(IN)    :: ParIdx
    INTEGER,         INTENT(IN)    :: nParVal
    REAL(KIND=8),    INTENT(IN)    :: ParVal(nParVal)
    TYPE(ErrorType), INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: nISRFPar, errstat, w

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'UpdateISRFParameter'

    ! =====================================================================
    ! UpdateISRFParameter starts here 
    ! =====================================================================

    ! Get the number of ISRF Parameters
    nISRFPar = Window%ISRF_Function%nPar

    ! Raise Exception if number of ISRF parameters is exceeded
    IF(ParIdx .LT. 1 .OR. ParIdx .Gt. nISRFPar) THEN
      CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName,&
                            Message_in='invalid value of ParIdx'                  )
    ENDIF

    ! Update Parameter based on window settings
    IF(Window%Settings%ISRF%PixDepPar) THEN

      IF(nParVal .EQ. 1) THEN

        ! Assume value is constant across grid
        Window%ISRFPar_L1Grid(:,ParIdx) = ParVal(1)
        Window%ISRFPar_ConvolGrid(:,ParIdx) = ParVal(1)

      ELSEIF(nParVal .EQ. Window%nL1_Wvl) THEN

        ! Set value for L1 Grid
        Window%ISRFPar_L1Grid(:,ParIdx) = ParVal

        ! Interpolate over convolution grid
        CALL BSPLINE_EdgeFill(Window%L1_Wvl, ParVal, nParVal,Window%Convol_Wvl,&
                              Window%ISRFPar_ConvolGrid(:,ParIdx),             &
                              Window%nConvol_Wvl, errstat                      )

      ELSE
        CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName, &
               Message_in='# of ISRF pix-dep par values must = 1 or # L1 wavelengths')
      ENDIF

      ! Compute HW1E on L1 Grid
      DO w=1,Window%nL1_Wvl
        CALL Window%ISRF_Function%DetermineHW1E(                                  &
                           Window%ISRFPar_L1Grid(w,:),Window%Settings%Convol_dWvl,&
                           Window%nConvol_Wvl,Window%ISRF_hw1e_L1Grid(w)          )
      ENDDO

      ! Now on the convolution Grid
      DO w=1,Window%nConvol_Wvl
        CALL Window%ISRF_Function%DetermineHW1E(                                      &
                           Window%ISRFPar_ConvolGrid(w,:),Window%Settings%Convol_dWvl,&
                           Window%nConvol_Wvl,Window%ISRF_hw1e_ConvolGrid(w)          )
      ENDDO

    ELSE
      
      IF(nParVal .EQ. 1) THEN
        Window%ISRFPar_Fixed(ParIdx) = ParVal(1)
      ELSE
        CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName, &
                          Message_in='Inputted more than 1 value for fixed parameter')
      ENDIF
      
      ! Update HW1E
      CALL Window%ISRF_Function%DetermineHW1E(                               &
                         Window%ISRFPar_Fixed,Window%Settings%Convol_dWvl,   &
                         Window%nConvol_Wvl,Window%ISRF_hw1e_Fixed           )

    ENDIF
    
  END SUBROUTINE UpdateISRFParameter

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: LoadInputFileISRF
  ! 
  ! DESCRIPTION: Loads an ISRF lookup table supplied for a given 
  !              instrument

  SUBROUTINE LoadInputFileISRF(Window,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),     INTENT(INOUT) :: Window
    TYPE(ErrorType),   INTENT(INOUT) :: Error  

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                 :: rcode, ncid, gid, dimid(4)
    INTEGER                 :: nPar, vid
    CHARACTER(LEN=maxChar)  :: tmpchar
    INTEGER(KIND=2)         :: Infile_nBand
    CHARACTER(LEN=100)      :: nBandStr
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadInputFileISRF'

    ! =====================================================================
    ! LoadInputFileISRF starts here 
    ! =====================================================================

    ! Check for error
    IF(CheckError(Error)) RETURN

    ! Open the file
    rcode = nf_open(TRIM(ADJUSTL(Window%Settings%ISRF%Infile)), NF_SHARE, ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
               'nf_open:'//TRIM(ADJUSTL(Window%Settings%ISRF%Infile)) )
    
    ! Get the band string
    WRITE(nBandStr,'(I100)') Window%Settings%ISRF%Infile_BandIndex

    ! Attach Band
    rcode = nf_inq_ncid(ncid,'Band' // TRIM(ADJUSTL(nBandStr)),gid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                        'inq_ncid:'//'Band' // TRIM(ADJUSTL(nBandStr)))

    ! Read the ISRF Type Index
    rcode = nf_inq_varid(gid, 'ISRFTypeIndex', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                               'nf_inq_varid(ISRF_LUT):ISRFTypeIndex' )
    rcode = nf_get_vara_int2(gid, vid, (/1/), (/1/), Window%Settings%ISRF%OptionIndex)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                               'nf_get_var(ISRF_LUT):ISRFTypeIndex'   )

    ! OptionIndex > 0 Implies parameterized ISRF
    IF(Window%Settings%ISRF%OptionIndex .GT. 0) THEN
      
      ! Get dimensions
      rcode = nf_inq_varid(gid,'ISRFParameters', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:ISRFParameters')
      rcode = nf_inq_vardimid(gid, vid, dimid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_vardimid:ISRFParameters')
      rcode = nf_inq_dim(gid,dimid(1),tmpchar,Window%nw_lut)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_inq_dim:ISRFParameters')
      rcode = nf_inq_dim(gid,dimid(2),tmpchar,Window%nx_lut)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_inq_dim:ISRFParameters')
      rcode = nf_inq_dim(gid,dimid(3),tmpchar,Window%nT_lut)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_inq_dim:ISRFParameters')
      rcode = nf_inq_dim(gid,dimid(4),tmpchar,Window%np_lut)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_inq_dim:ISRFParameters')

      ! Allocate Arrays
      ALLOCATE(Window%ISRF_LUT_Par(Window%nw_lut,Window%nx_lut,Window%nT_lut,Window%np_lut))
      ALLOCATE(Window%ISRF_LUT_T(Window%nT_lut))
      
      ! Determine grid dependence 
      Window%ISRF_LUT_wdep = Window%nw_lut .GT. 1
      Window%ISRF_LUT_Tdep = Window%nT_lut .GT. 1
      
      ! Set Wavelength grid dependence 
      Window%Settings%ISRF%PixDepPar = Window%ISRF_LUT_wdep .OR. &
                                       Window%nx_lut .GT. 1 ! Must also update param. with cross track
      
      ! Allocate array for fixed option
      IF(.NOT. ALLOCATED(Window%ISRFPar_Fixed)) &
        ALLOCATE(Window%ISRFPar_Fixed(Window%np_lut))
      
      ! Read ISRF Parameters
      rcode = nf_get_vara_double(gid, vid, (/1,1,1,1/), (/Window%nw_lut,Window%nx_lut,Window%nT_lut,Window%np_lut/),&
                                Window%ISRF_LUT_Par)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:ISRFParameters')
      
      ! Read Temperature grid
      rcode = nf_inq_varid(gid,'Temperature', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:Temperature')
      rcode = nf_get_vara_double(gid, vid, (/1/), (/Window%nT_lut/), Window%ISRF_LUT_T)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:ISRFParameters')
      
    ENDIF

    ! We are done with the file
    rcode = nf_close(ncid)

  END SUBROUTINE LoadInputFileISRF

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: GetInputFileISRF
  ! 
  ! DESCRIPTION: Sets the ISRF by interpolating over the supplied 
  !              lookup table

  SUBROUTINE GetInputFileISRF(Window,Level1,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),     INTENT(INOUT) :: Window
    TYPE(L1Type),      INTENT(IN)    :: Level1
    TYPE(ErrorType),   INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: idT0,idx0, p, bidx, i, xtr
    REAL(KIND=8) :: Twt0, Twt1, xwt0, xwt1
    REAL(KIND=8) :: OptBenchTemp
    REAL(KIND=8) :: ISRFPar_band(Window%nw_lut,Window%np_lut)
    REAL(KIND=8), ALLOCATABLE :: iGrid_band(:)
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'GetInputFileISRF'
    
    ! =====================================================================
    ! GetInputFileISRF starts here 
    ! =====================================================================

    ! OptionIndex > 0 Implies parameterized ISRF: Now determine parameters
    IF(Window%ISRF_Function%IsLUT_ISRF) THEN

      ! Allocate ISRF parameters
      IF(ALLOCATED(Window%ISRFPar_L1Grid)) DEALLOCATE(Window%ISRFPar_L1Grid)
      ALLOCATE(Window%ISRFPar_L1Grid(Window%nL1_Wvl,Window%ISRF_Function%nPar))

      ! Shift Init to 0, Squeeze Factor to 1.0
      Window%ISRFPar_L1Grid(:,1) = 0.0d0
      Window%ISRFPar_L1Grid(:,2) = 1.0d0

      ! Same for the convolution grid
      Window%ISRFPar_ConvolGrid(:,1) = 0.0d0
      Window%ISRFPar_ConvolGrid(:,2) = 1.0d0 

    ELSE
      ! Cross track index
      xtr = Window%ISRF_Function%xtrk_idx

      ! Get Optical Bench Temperature
      IF(Level1%Geolocation%OpticalBenchTemperature .GT. 0.0) THEN
        OptBenchTemp = Level1%Geolocation%OpticalBenchTemperature
      ELSE
        OptBenchTemp = Window%ISRF_LUT_T(1) ! Set to coldest
      ENDIF

      ! Allocate wvl-dep ISRF on L1 grid needed
      IF(Window%ISRF_LUT_wdep) THEN

        ! Allocate ISRF parameters
        IF(ALLOCATED(Window%ISRFPar_L1Grid)) DEALLOCATE(Window%ISRFPar_L1Grid)
        ALLOCATE(Window%ISRFPar_L1Grid(Window%nL1_Wvl,Window%np_lut))

        ! Get Level 1 Band index for window
        bidx = Window%Settings%BandIndex

      ENDIF

      ! ------------------------
      ! Variable Wavelength Case
      ! ------------------------
      IF(Window%ISRF_LUT_wdep) THEN

        IF(Window%ISRF_LUT_Tdep) THEN

          ! Compute the interpolation weights for the temperature grid
          CALL LinearInterpolationWeights(Window%nT_lut,Window%ISRF_LUT_T,&
                                          OptBenchTemp,idT0,Twt0,Twt1)
          
          ! Compute Band ISRF
          DO p=1,Window%np_lut
            ISRFPar_band(:,p) = Window%ISRF_LUT_Par(:, xtr,idT0,  p)*Twt0 &
                              + Window%ISRF_LUT_Par(:, xtr,idT0+1,p)*Twt1
          ENDDO
          
          ! Get the ISRF Parameters on the L1 Grid
          Window%ISRFPar_L1Grid(:,:) = ISRFPar_band(Window%L1_iw0:Window%L1_iwf,:)

          ! Get ISRF Parameters on convolution grid
          DO i=1,Window%nConvol_Wvl

            ! Get fractional pixel index
            CALL LinearInterpolationWeights(Level1%Spectrum(bidx)%nwvl,&
                                            Level1%Spectrum(bidx)%Wavelength,&
                                            Window%Convol_Wvl(i),idx0,xwt0,xwt1)

            ! Compute ISRF
            Window%ISRFPar_ConvolGrid(i,:) = ISRFPar_band(idx0  ,:)*xwt0 &
                                           + ISRFPar_band(idx0+1,:)*xwt1

          ENDDO
          
        ELSE

          ! Store ISRF Parameters on L1 Grid
          DO p=1,Window%np_lut
            Window%ISRFPar_L1Grid(:,p) = &
              Window%ISRF_LUT_Par(Window%L1_iw0:Window%L1_iwf,xtr,1,p)
          ENDDO
          
          ! Interpolate to convolution Grid
          DO i=1,Window%nConvol_Wvl

            ! Get fractional pixel index
            CALL LinearInterpolationWeights(Level1%Spectrum(bidx)%nwvl,        &
                                            Level1%Spectrum(bidx)%Wavelength,  &
                                            Window%Convol_Wvl(i),idx0,xwt0,xwt1)

            ! Compute ISRF
            Window%ISRFPar_ConvolGrid(i,:) = Window%ISRF_LUT_Par(idx0  ,xtr,1,:)*xwt0 &
                                          + Window%ISRF_LUT_Par(idx0+1,xtr,1,:)*xwt1
            
          ENDDO

        ENDIF

      ! ---------------------
      ! Fixed Wavelength Case
      ! ---------------------
      ELSE

        IF(Window%ISRF_LUT_Tdep) THEN

          ! Compute the interpolation weights
          CALL LinearInterpolationWeights(Window%nT_lut,Window%ISRF_LUT_T,&
                                          OptBenchTemp,idT0,Twt0,Twt1)

          ! Set the fixed ISRF
          Window%ISRFPar_Fixed(:) = Window%ISRF_LUT_Par(1, xtr, idT0  ,:)*Twt0 &
                                  + Window%ISRF_LUT_Par(1, xtr, idT0+1,:)*Twt1

        ELSE

          ! Set the fixed ISRF
          Window%ISRFPar_Fixed(:) = Window%ISRF_LUT_Par(1,xtr,1,:)

        ENDIF

      ENDIF

    ENDIF

  END SUBROUTINE GetInputFileISRF

  ! Helper Functions for dimensions
  INTEGER(KIND=4) FUNCTION GetWindowNumStreams(Window,Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    TYPE(ErrorType), INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'GetWindowNumStreams'

    ! =====================================================================
    ! GetWindowNumStreams starts here 
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Use the window RTM Settings
    IF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-LBL' ) THEN
      GetWindowNumStreams = Window%RTMSettings%SimOpt%VLIDORT%nstreams
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'VL-PCA' ) THEN
      GetWindowNumStreams = Window%RTMSettings%SimOpt%VLIDORT_PCA%nstreams
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'FIRSTORD' ) THEN
      GetWindowNumStreams = 0
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. '2STRM' ) THEN
      GetWindowNumStreams = 2
    ELSEIF( TRIM(ADJUSTL(Window%RTMSettings%RTMName)) .EQ. 'TRANSMISSION' ) THEN
      GetWindowNumStreams = 0
    ELSE
      print*,'In Window Model'
      print*,'Unrecognized Radiative Transfer Model'
      STOP 1
    ENDIF
    
  END FUNCTION GetWindowNumStreams
  
  SUBROUTINE ApplyISRF(Window, SpectrumHR, SpectrumLR, Error, DoI0Corr_in, DerivativeIndex)
    
    USE interpolation_module, ONLY : BSPLINE
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),     INTENT(IN)    :: Window
    REAL(KIND=8),      INTENT(IN)    :: SpectrumHR(Window%nConvol_Wvl)
    REAL(KIND=8),      INTENT(OUT)   :: SpectrumLR(Window%nL1_Wvl)
    TYPE(ErrorType),   INTENT(INOUT) :: Error
    LOGICAL, OPTIONAL, INTENT(IN)    :: DoI0Corr_in
    INTEGER,OPTIONAL,  INTENT(IN)    :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: SpectrumConvoled(Window%nConvol_Wvl)
    REAL(KIND=8) :: SpectrumSP(Window%nConvol_Wvl)
    REAL(KIND=8) :: RadTmp(Window%nConvol_Wvl)
    REAL(KIND=8) :: RadTmpLR(Window%nConvol_Wvl),RadTmpLRSP(Window%nConvol_Wvl)
    INTEGER      :: errstat_sp, w, deriv_idx
    LOGICAL      :: DoI0Corr
    REAL(KIND=8) :: Ns, xs_max
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ApplyISRF'

    ! =====================================================================
    ! ApplyISRF begins here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex

    ! Flag for I0 Correction
    DoI0Corr = .FALSE. ; IF(PRESENT(DoI0Corr_in)) DoI0Corr = DoI0Corr_in
    
    ! Get the max absolute spectrum value
    xs_max = MAXVAL(ABS(SpectrumHR))
    
    IF(xs_max .GT. TINY(0.0d0)) THEN

      ! If its pixel dependent it is likely faster to compute on the L1 grid
      IF(Window%Settings%ISRF%PixDepPar) THEN

        IF(DoI0Corr) THEN
          STOP 'Not implemented'
        ELSE
          CALL ApplyISRF_LR(Window, SpectrumHR, SpectrumLR, Error, &
                            DerivativeIndex=deriv_idx              )
        ENDIF

      ELSE

        IF(DoI0Corr) THEN
          
          ! Apply ISRF with I0 - output on convolution grid
          CALL ApplyISRF_HRSP_I0(Window, SpectrumHR, SpectrumConvoled, SpectrumSP,&
                                 Error, DerivativeIndex=deriv_idx                 )
          
        ELSE
        
          ! Apply ISRF - output to convolution grid
          CALL ApplyISRF_HRSP(Window, SpectrumHR, SpectrumConvoled, SpectrumSP,&
                              Error, DerivativeIndex=deriv_idx                 )
          
        ENDIF
        
        ! Interpolate to L1 Grid
        CALL SPLINT1(Window%Convol_Wvl, SpectrumConvoled, SpectrumSP, Window%nConvol_Wvl,&
                    Window%L1_Wvl,     SpectrumLR,                   Window%nL1_Wvl      )
      ENDIF

    ELSE
      
      SpectrumLR(:) = 0.0d0
      
    ENDIF
     
  END SUBROUTINE ApplyISRF

  ! Subroutine for applying the ISRF for each output grid position
  SUBROUTINE ApplyISRF_LR(Window, SpectrumHR, SpectrumConvol, Error, DerivativeIndex)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    REAL(KIND=8),    INTENT(IN)    :: SpectrumHR(Window%nConvol_Wvl)
    REAL(KIND=8),    INTENT(OUT)   :: SpectrumConvol(Window%nL1_Wvl)
    TYPE(ErrorType), INTENT(INOUT) :: Error
    INTEGER,OPTIONAL,INTENT(IN)    :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: p, deriv_idx
    REAL(KIND=8),   ALLOCATABLE :: ISRFPar(:,:)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ApplyISRF_LR'
    
    ! =====================================================================
    ! ApplyISRF_LR begins here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex
    
    ! Do the convolution
    IF(Window%Settings%ISRF%PixDepPar) THEN
      CALL Window%ISRF_Function%Convolve_Wdep(Window%nConvol_Wvl, Window%nL1_Wvl,&
                                 Window%Convol_Wvl, SpectrumHR,                  &
                                 Window%ISRFPar_L1Grid, Window%ISRF_hw1e_L1Grid, & 
                                 Window%L1_Wvl, SpectrumConvol,                  &
                                 DerivativeIndex=deriv_idx                       )

      ! CALL Convolution_WvlDepPar(Window%ISRF_Function,                          &
      !                            Window%nConvol_Wvl, Window%nL1_Wvl,            &
      !                            Window%Convol_Wvl, SpectrumHR,                 &
      !                            Window%ISRFPar_L1Grid, Window%ISRF_hw1e_L1Grid,& 
      !                            Window%L1_Wvl, SpectrumConvol,                 &
      !                            DerivativeIndex=deriv_idx                      )
    ELSE

      ! Kludge - Use the fixed parameters
      ALLOCATE(ISRFPar(Window%nL1_Wvl,Window%ISRF_Function%nPar))
      DO p=1,Window%ISRF_Function%nPar
        ISRFPar(:,p) = Window%ISRFPar_Fixed(p)
      ENDDO

      ! Do the Convolution
      CALL Window%ISRF_Function%Convolve_Wdep(Window%nConvol_Wvl, Window%nL1_Wvl,   &
                                              Window%Convol_Wvl, SpectrumHR,ISRFPar,&
                                              Window%ISRF_hw1e_L1Grid,              &
                                              Window%L1_Wvl,SpectrumConvol,         &
                                              DerivativeIndex=deriv_idx             )

      ! CALL Convolution_WvlDepPar(Window%ISRF_Function,                 &
      !                            Window%nConvol_Wvl, Window%nL1_Wvl,   &
      !                            Window%Convol_Wvl, SpectrumHR,ISRFPar,&
      !                            Window%ISRF_hw1e_L1Grid,              &
      !                            Window%L1_Wvl,SpectrumConvol,         &
      !                            DerivativeIndex=deriv_idx             )
      
    ENDIF

  END SUBROUTINE ApplyISRF_LR
  
  SUBROUTINE ApplyISRF_LR_I0(Window, SpectrumHR, SpectrumConvol, Error, DerivativeIndex)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    REAL(KIND=8),    INTENT(IN)    :: SpectrumHR(Window%nConvol_Wvl)
    REAL(KIND=8),    INTENT(OUT)   :: SpectrumConvol(Window%nL1_Wvl)
    TYPE(ErrorType), INTENT(INOUT) :: Error
    INTEGER,OPTIONAL,INTENT(IN)    :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: p, w, deriv_idx
    REAL(KIND=8) :: xs_max, Ns, RadTmp(Window%nConvol_Wvl)
    REAL(KIND=8) :: RadTmp_LR(Window%nL1_Wvl)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ApplyISRF_LR'
    

    ! =====================================================================
    ! ApplyISRF_LR_I0 begins here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex

    ! Get the max spectrum value
    xs_max = MAXVAL(SpectrumHR)
    
    IF(xs_max .GT. TINY(0.0d0)) THEN
      
      ! I0 "weak absorber" correction - Make peak optical depth = 1e-4
      ! Approximation is valid for weak absorptions 
      Ns = 1e-4/xs_max
        
      ! Compute HR "Radiance" spectrum
      DO w=1,Window%nConvol_Wvl
        RadTmp(w) = Window%I0_ConvolGrid(w)*EXP(-1.0*Ns*SpectrumHR(w))
      ENDDO
      
      ! Convolve RadTmp
      CALL ApplyISRF_LR(Window,RadTmp,RadTmp_LR,Error,DerivativeIndex=deriv_idx)
      
      DO w=1,Window%nL1_Wvl
        SpectrumConvol(w) = LOG(Window%I0_L1Grid_ISRF(w)/RadTmp_LR(w))/Ns
      ENDDO
      
    ELSE
      SpectrumConvol(:) = 0.0d0
    ENDIF
    

  END SUBROUTINE ApplyISRF_LR_I0

  ! Subroutine for applying the ISRF for each convolution grid position
  SUBROUTINE ApplyISRF_HRSP(Window, SpectrumHR, SpectrumConvol, SpectrumConvSP, Error, DerivativeIndex)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    REAL(KIND=8),    INTENT(IN)    :: SpectrumHR(Window%nConvol_Wvl)
    REAL(KIND=8),    INTENT(OUT)   :: SpectrumConvol(Window%nConvol_Wvl)
    REAL(KIND=8),    INTENT(OUT)   :: SpectrumConvSP(Window%nConvol_Wvl)
    TYPE(ErrorType), INTENT(INOUT) :: Error
    INTEGER,OPTIONAL,INTENT(IN)    :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------

    ! For error checking
    INTEGER :: deriv_idx
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ApplyISRF_HRSP'

    ! =====================================================================
    ! ApplyISRF_HRSP begins here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex

    IF( Window%Settings%ISRF%PixDepPar ) THEN
      CALL Window%ISRF_Function%Convolve_Wdep(Window%nConvol_Wvl,   &
                                 Window%nConvol_Wvl,                &
                                 Window%Convol_Wvl, SpectrumHR,     &
                                 Window%ISRFPar_ConvolGrid,         &
                                 Window%ISRF_hw1e_ConvolGrid,       &
                                 Window%Convol_Wvl, SpectrumConvol, &
                                 DerivativeIndex=deriv_idx          )
    ELSE
      CALL Window%ISRF_Function%Convolve_Fixed(Window%nConvol_Wvl,      &
                                 Window%Convol_Wvl,                     & 
                                 SpectrumHR, Window%ISRFPar_Fixed,      &
                                 Window%ISRF_hw1e_Fixed, SpectrumConvol,&
                                 DerivativeIndex=deriv_idx              )
    ENDIF
    
    ! Compute basis spline expansion
    CALL SPLINE1(Window%Convol_Wvl,  SpectrumConvol,  &
                 Window%nConvol_Wvl, SpectrumConvSP   )
    
  END SUBROUTINE ApplyISRF_HRSP

  SUBROUTINE ApplyISRF_HRSP_I0(Window, SpectrumHR, SpectrumConvol, SpectrumConvSP, Error,DerivativeIndex)
  
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(WinType),   INTENT(IN)    :: Window
    REAL(KIND=8),    INTENT(IN)    :: SpectrumHR(Window%nConvol_Wvl)
    REAL(KIND=8),    INTENT(OUT)   :: SpectrumConvol(Window%nConvol_Wvl)
    REAL(KIND=8),    INTENT(OUT)   :: SpectrumConvSP(Window%nConvol_Wvl)
    TYPE(ErrorType), INTENT(INOUT) :: Error
    INTEGER,OPTIONAL,INTENT(IN)    :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: RadTmp(Window%nConvol_Wvl)
    REAL(KIND=8) :: RadTmpLR(Window%nConvol_Wvl)
    REAL(KIND=8) :: RadTmpLRSP(Window%nConvol_Wvl)
    REAL(KIND=8) :: Ns, xs_max
    INTEGER      :: w, deriv_idx
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ApplyISRF_HRSP_I0'

    ! =====================================================================
    ! ApplyISRF_HRSP_I0 begins here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex

    ! Get the max spectrum value
    xs_max = MAXVAL(SpectrumHR)
    
    IF(xs_max .GT. TINY(0.0d0)) THEN
      
      ! I0 "weak absorber" correction - Make peak optical depth = 1e-4
      ! Approximation is valid for weak absorptions 
      Ns = 1e-4/xs_max
      
      ! Compute HR "Radiance" spectrum
      DO w=1,Window%nConvol_Wvl
        RadTmp(w) = Window%I0_ConvolGrid(w)*EXP(-1.0*Ns*SpectrumHR(w))
      ENDDO
      
      ! Convolve RadTmp
      CALL ApplyISRF_HRSP(Window, RadTmp, RadTmpLR, RadTmpLRSP, Error, &
                          DerivativeIndex=deriv_idx                    )
      
      ! Compute effective cross section
      DO w=1,Window%nConvol_Wvl
        SpectrumConvol(w) = DLOG(Window%I0_ConvolGrid_ISRF(w)/RadTmpLR(w))/Ns
      ENDDO
      
    ELSE
      SpectrumConvol(:) = 0.0d0
    ENDIF
    
    ! Compute basis spline expansion
    CALL SPLINE1(Window%Convol_Wvl,  SpectrumConvol,&
                 Window%nConvol_Wvl, SpectrumConvSP )
    
  END SUBROUTINE ApplyISRF_HRSP_I0
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: InitISRFOpt
  ! 
  ! DESCRIPTION: Initialize empty ISRF Option array

  SUBROUTINE InitISRFOpt(ISRFOpt)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(ISRFOptType),    INTENT(INOUT) :: ISRFOpt

    ! ---------------
    ! local variables
    ! ---------------
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitISRFOpt'

    ! ============================================================
    ! InitISRFOpt starts here
    ! ============================================================
    
    ISRFOpt%Name               = ''
    ISRFOpt%PixDepPar          = .FALSE.
    ISRFOpt%OptionIndex        = -1
    ISRFOpt%Infile             = ''
    ISRFOpt%Infile_BandIndex   = -1  

    ! New
    ISRFOpt%nFixedPar = 0
    IF(ALLOCATED(ISRFOpt%FixedPar)) DEALLOCATE(ISRFOpt%FixedPar)
    
  END SUBROUTINE InitISRFOpt
  
END MODULE window_module

