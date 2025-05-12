MODULE isrf_module

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  !
  ! MODULE: isrf_module
  ! 
  ! DESCRIPTION: The ISRF Module contains the routines needed to perform 
  !              convolutions between spectra and instrument spectral 
  !              response functions and their derivatives w.r.t. ISRF 
  !              parameters
  !
  !              Current ISRFs:
  !                - Supergaussian: super_gaussian_sf (fixed parameters)
  !                                 SuperGaussISRF_WvlDep (wvl-dep parameters)
  !                - TROPOMI-Like: tropomi_isrf (fixed parameters)
  !                                tropomi_isrf_wvldep (wvl-dep parameters)
  ! 
  ! The above routines assume an evenly-spaced wavelength grid
  !###################################################################

  USE parameters_module
  USE error_module,         ONLY : ErrorType, RaiseFatalError
  USE interpolation_module, ONLY : BSPLINE, BSPLINE_EdgeFill,&
                                   SPLINE1, SPLINT1,         &
                                   LinearInterpolationWeights
  USE netcdf_module,        ONLY : CheckNetCDFErrorStatus

  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'

  ! Holds Routines
  TYPE ISRF_FunctionType

    ! Things Needed for the ISRF Parameterization
    PROCEDURE(ISRF_FunctionInterface), NOPASS, POINTER :: ParameterizedISRF
    LOGICAL                                            :: IsDeltaFunction
    INTEGER                                            :: nPar
    REAL(KIND=8)                                       :: NumConvolutionHW1E

    ! Things Needed for the LUT Approach
    LOGICAL                   :: IsLUT_ISRF
    REAL(KIND=8), ALLOCATABLE :: LUT_ISRF(:,:,:) ! w x w_convol
    REAL(KIND=8), ALLOCATABLE :: LUT_ISRF_SP(:,:,:) ! w x w_convol
    REAL(KIND=8), ALLOCATABLE :: LUT_CentralWvl(:)
    REAL(KIND=8), ALLOCATABLE :: LUT_DeltaWvl(:)
    REAL(KIND=8)              :: LUT_DeltaMin
    REAL(KIND=8)              :: LUT_DeltaMax
    INTEGER                   :: LUT_iw0 ! Start Index For LUT
    INTEGER                   :: LUT_dwmx
    INTEGER                   :: LUT_wmx
    INTEGER                   :: LUT_xmx
    INTEGER                   :: xtrk_idx
    
    ! For the convolution
    CONTAINS
      PROCEDURE :: Initialize => InitISRFFunctionType
      PROCEDURE :: EVAL => EvaluateISRF
      PROCEDURE :: Convolve_Wdep => Convolution_WvlDepPar
      PROCEDURE :: Convolve_Fixed => Convolution_FixedPar
      PROCEDURE :: DetermineHW1E => DetermineHW1E
      

  ENDTYPE ISRF_FunctionType

  ! Define interface to ISRF Function Call
  INTERFACE
    SUBROUTINE ISRF_FunctionInterface(nx,x,np,p,S,ScaleBydW_in)
      INTEGER ,                 INTENT(IN)  :: nx
      REAL(KIND=8),             INTENT(IN)  :: x(nx) ! Assumes even spacing
      INTEGER,                  INTENT(IN)  :: np
      REAL(KIND=8),             INTENT(IN)  :: p(np)
      REAL(KIND=8),             INTENT(OUT) :: S(nx)
      LOGICAL, OPTIONAL,        INTENT(IN)  :: ScaleBydW_in ! Default=TRUE
    END SUBROUTINE ISRF_FunctionInterface
  END INTERFACE

  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'isrf_module'
  PRIVATE :: ModuleName

  CONTAINS
  
  SUBROUTINE InitISRFFunctionType(self,OptionIndex,Error,     &
                                  LookupTableFile,BandIndex,  &
                                  Convol_dWvl,ConvolWidth_HW1E)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(ISRF_FunctionType),   INTENT(INOUT) :: self
    INTEGER,                    INTENT(IN)    :: OptionIndex
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: LookupTableFile
    INTEGER,          OPTIONAL, INTENT(IN)    :: BandIndex
    REAL(KIND=8),     OPTIONAL, INTENT(IN)    :: Convol_dWvl
    REAL(KIND=8),     OPTIONAL, INTENT(IN)    :: ConvolWidth_HW1E

    ! ---------------
    ! local variables
    ! ---------------
    CHARACTER(LEN=10)      :: IdxStr
    CHARACTER(LEN=maxChar) :: tmpchar
    INTEGER                :: ncid, gid, rcode, vid, dimid(3), errstat
    INTEGER                :: dwmx, w, i, j, n_below, n_above
    REAL(KIND=8)           :: NormFac
    
    ! Local variables to store the delta grid
    REAL(KIND=8), ALLOCATABLE :: DeltaWvl(:), DeltaWvl_Convol(:)
    REAL(KIND=8), ALLOCATABLE :: ISRF(:,:,:)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitISRFFunctionType'

    ! =====================================================================
    ! InitISRFType begins here
    ! =====================================================================

    ! Initialize Values
    self%ParameterizedISRF  => null()
    self%IsDeltaFunction = .FALSE.
    self%nPar = 0
    self%IsLUT_ISRF = .FALSE.
    self%LUT_iw0 = 0
    self%LUT_wmx = 0
    self%LUT_xmx = 0
    self%LUT_dwmx = 0

    ! Set Convolution window if passed
    self%NumConvolutionHW1E = 10.0d0
    IF(PRESENT(ConvolWidth_HW1E)) self%NumConvolutionHW1E = ConvolWidth_HW1E
    
    ! ---------------------------------------------
    ! Read Lookup table ISRF if input file provided
    ! ---------------------------------------------
    IF(PRESENT(LookupTableFile)) THEN
      
      ! Band Index Must be present
      IF(.NOT. PRESENT(BandIndex)) THEN
        CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName, &
                              Message_in='BandIndex Must be present when using LUT'  )
      ENDIF

      IF(.NOT. PRESENT(Convol_dWvl)) THEN
        CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName, &
                              Message_in='Convol_dWvl Must be present when using LUT'  )
      ENDIF

      ! Set Flag for LUT
      self%IsLUT_ISRF = .TRUE.

      ! Open File
      rcode = nf_open(TRIM(ADJUSTL(LookupTableFile)), NF_SHARE, ncid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
              'nf_open:'//TRIM(ADJUSTL(LookupTableFile)))

      ! Get the band string
      WRITE(IdxStr,'(I10)') BandIndex

      ! Attach Band
      rcode = nf_inq_ncid(ncid,'Band' // TRIM(ADJUSTL(IdxStr)),gid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                            'inq_ncid:'//'Band' // TRIM(ADJUSTL(IdxStr)))

      ! Get Dimensions of LUT
      rcode = nf_inq_varid(gid,'ISRF', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_varid:ISRF')
      rcode = nf_inq_vardimid(gid, vid, dimid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_vardimid:ISRF')
      rcode = nf_inq_dim(gid,dimid(1),tmpchar,self%LUT_dwmx) ! Dispersed Grid
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:ISRF')
      rcode = nf_inq_dim(gid,dimid(2),tmpchar,self%LUT_wmx) ! Central Wavelength
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:ISRF')
      rcode = nf_inq_dim(gid,dimid(3),tmpchar,self%LUT_xmx) ! Cross Track 
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_dim:ISRF')

      ! Load Central Wavelength
      ALLOCATE(self%LUT_CentralWvl(self%LUT_wmx))
      rcode = nf_inq_varid( gid, 'CenterWavelength', vid )
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_var:CenterWavelength')
      rcode = nf_get_var_double(gid, vid, self%LUT_CentralWvl)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:CentralWavelength')

      ! Load the delta wavelength Grid
      ALLOCATE(self%LUT_DeltaWvl(self%LUT_dwmx))
      rcode = nf_inq_varid( gid, 'DeltaWavelength', vid )
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_var:DeltaWavelength')
      rcode = nf_get_var_double(gid, vid, self%LUT_DeltaWvl)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:DeltaWavelength')

      ! Also save Min/Max
      self%LUT_DeltaMin = self%LUT_DeltaWvl(1)
      self%LUT_DeltaMax = self%LUT_DeltaWvl(self%LUT_dwmx)
      
      ! Load the ISRF
      ALLOCATE(self%LUT_ISRF(self%LUT_dwmx,self%LUT_wmx,self%LUT_xmx))
      rcode = nf_inq_varid(gid, 'ISRF', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'inq_var:ISRF')
      rcode = nf_get_var_double(gid, vid, self%LUT_ISRF)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'get_var:ISRF')
      
      ! Get the Basis Spline expansion coefficients
      ALLOCATE(self%LUT_ISRF_SP(self%LUT_dwmx,self%LUT_wmx,self%LUT_xmx))
      DO j=1,self%LUT_xmx
      DO i=1,self%LUT_wmx
        CALL SPLINE1(self%LUT_DeltaWvl,      &
                     self%LUT_ISRF(:,i,j),   &
                     self%LUT_dwmx,          &
                     self%LUT_ISRF_SP(:,i,j),&
                     SubroutineName          )
      ENDDO
      ENDDO
    
  ENDIF


  ! -------------------------------------------------------
  ! Set Parameter Dimensions / Point to Parameterized ISRFs
  ! -------------------------------------------------------
      
  ! (0) Lookup Table
  IF(OptionIndex .EQ. 0) THEN

    ! Shift/Squeeze
    self%nPar = 2

  ! (1) SuperGaussian
  ELSEIF(OptionIndex .EQ. 1) THEN
    self%ParameterizedISRF => SuperGauss_ISRF_Func
    self%nPar = 3

  ! (2) TROPOMI 7-Parameter (~Asymm Voigt)
  ELSEIF(OptionIndex .EQ.2 ) THEN
    self%ParameterizedISRF => TROPOMI_ISRF_Func
    self%nPar = 7

  ELSE
    WRITE(IdxStr,'(I10)') OptionIndex
    CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName, &
                          Message_in='ISRF Option Index is not valid'            )

  ENDIF
  
  END SUBROUTINE InitISRFFunctionType

  ! =======================================================================
  ! =======================================================================
  !                ISRF Function Options For referencing
  ! =======================================================================
  ! =======================================================================

  SUBROUTINE SuperGauss_ISRF_Func(nx,x,np,p,S,ScaleBydW_in)

    ! Parameter Definitions
    ! p(1): hw1e
    ! p(2): e_asym
    ! p(3): g_shap
    
    ! =========================================================================
    !
    ! The asymetric Gaussian g(x) is defined as
    !                   _                                            _
    !                  |   |            x^2                  |^g_shap |
    !      g(x) =  EXP | - |---------------------------------|        |
    !                  |_  | (hw1e * (1 + SIGN(x)*e_asym))   |       _|
    !
    ! g(x) becomes symmetric for E_ASYM = 0.
    !
    ! =========================================================================

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER ,          INTENT(IN)  :: nx 
    REAL(KIND=8),      INTENT(IN)  :: x(nx) ! Assumes even spacing
    INTEGER,           INTENT(IN)  :: np
    REAL(KIND=8),      INTENT(IN)  :: p(np)
    REAL(KIND=8),      INTENT(OUT) :: S(nx)
    LOGICAL, OPTIONAL, INTENT(IN)  :: ScaleBydW_in ! Default=TRUE

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i
    REAL(KIND=8) :: sg_nfac
    LOGICAL      :: ScaleBydW
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SuperGauss_ISRF_Func'

    ! =====================================================================
    ! SuperGauss_Func begins here
    ! =====================================================================

    ! Default Scaling 
    ScaleBydW = .TRUE. ; IF(PRESENT(ScaleBydW_in)) ScaleBydW = ScaleBydW_in

    ! Compute Normalization Factor (asymm in width is equiv to normal SG)
    sg_nfac = p(3)/(2.0d0*p(1)*GAMMA(1.0d0/p(3)))

    DO i=1,nx
      S(i) = sg_nfac*EXP( -1.0d0*(ABS( x(i)/(p(1) + signdp(x(i))*p(2)) ) )**p(3) )
    ENDDO

    ! Check if we are returning S(x)*dx
    IF(ScaleBydW) THEN
      S(:) = S(:)*(x(2)-x(1))
    ENDIF

  END SUBROUTINE SuperGauss_ISRF_Func

  SUBROUTINE TROPOMI_ISRF_Func(nx,x,np,p,S,ScaleBydW_in)


    ! Parameter Definitions
    ! p(1): x0
    ! p(2): d
    ! p(3): s 
    ! p(4): w
    ! p(5): eta
    ! p(6): gamma 
    ! p(7): m

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER ,          INTENT(IN)  :: nx 
    REAL(KIND=8),      INTENT(IN)  :: x(nx) ! Assumes even spacing
    INTEGER,           INTENT(IN)  :: np
    REAL(KIND=8),      INTENT(IN)  :: p(np)
    REAL(KIND=8),      INTENT(OUT) :: S(nx)
    LOGICAL, OPTIONAL, INTENT(IN)  :: ScaleBydW_in ! Default=TRUE

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i
    REAL(KIND=8) :: P7_nfac, S_nfac, al_p, al_m,delt, delta, al_pre, sr2, m_exp, dw
    LOGICAL      :: ScaleBydW
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'TROPOMI_ISRF_FUNC'

    ! =====================================================================
    ! TROPOMI_ISRF_FUNC begins here
    ! =====================================================================

    ! Default Scaling 
    ScaleBydW = .TRUE. ; IF(PRESENT(ScaleBydW_in)) ScaleBydW = ScaleBydW_in

    ! Helper variables
    sr2 = SQRT(2.0d0)
    delta = sr2*p(3)/SQRT(Constants_pi*(1.0d0+p(3)*p(3)) )
    al_pre = SQRT(1.0d0-delta*delta)/p(2)
    m_exp = -1.0d0*p(7)
    
    ! Normalization Factors
    P7_nfac = p(5)*GAMMA(p(7)) &
            / (p(6)*SQRT(Constants_pi)*GAMMA(p(7)-0.5d0) )
    S_nfac = (1.0d0-p(5))/p(4)

    ! Compute Function for each pixel
    DO i=1,nx

      ! Lorentz Component (tails)
      S(i) = P7_nfac*(1.0d0+(x(i)-p(1))**2/p(6)/p(6))**m_exp

      ! Component from asym gauss
      al_p = al_pre*(x(i)-p(1)+0.5d0*p(4))+delta
      al_m = al_pre*(x(i)-p(1)-0.5d0*p(4))+delta
      S(i) = S(i) + S_nfac*( 0.5d0*( ERF(al_p/sr2)-ERF(al_m/sr2) ) &
                            -2.0d0*( tha(al_p,1.0d0,p(3),1.0d0)-tha(al_m,1.0d0,p(3),1.0d0) ) )
    ENDDO
    
    IF(ScaleBydW) THEN
      IF(nx > 1) THEN
        S(:) = S(:) * (x(2)-x(1))
      ELSE
        print*,'Cannot scale by dwvl as only 1 point',nx
      ENDIF
    ENDIF

  END SUBROUTINE TROPOMI_ISRF_FUNC

  SUBROUTINE EvaluateISRF(self,wvl0,nx,x,p,S,ScaleBydW_in)

    ! Parameter Definitions for LUT
    ! p(1): Wavelength grid shift
    ! p(2): Wavelength grid squeeze

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(ISRF_FunctionType), INTENT(IN)  :: self
    REAL(KIND=8),             INTENT(IN)  :: wvl0 ! Required by LUT Shape Only
    INTEGER ,                 INTENT(IN)  :: nx
    REAL(KIND=8),             INTENT(IN)  :: x(nx) ! Assumes even spacing
    REAL(KIND=8),             INTENT(IN)  :: p(self%nPar)
    REAL(KIND=8),             INTENT(OUT) :: S(nx)
    LOGICAL, OPTIONAL,        INTENT(IN)  :: ScaleBydW_in ! Default=TRUE

    ! These must be input for this case
    

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: x_in(nx), Stmp0(nx), Stmp1(nx), dx
    REAL(KIND=8) :: NormFac, xwt0, xwt1, x_test
    INTEGER      :: i_in(nx), w, n_in, idx0
    LOGICAL      :: ScaleBydW
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'EvaluateISRF'

    ! =====================================================================
    ! EvaluateISRF begins here
    ! =====================================================================
    
    ! Default Scaling 
    ScaleBydW = .TRUE. ; IF(PRESENT(ScaleBydW_in)) ScaleBydW = ScaleBydW_in

    IF(self%IsLUT_ISRF) THEN
      
      ! Initialize Return values
      S(:) = 0.0d0
        
      ! Grid Spacing
      dx = x(2) - x(1)

      ! First Interpolate Slit to regular grid
      n_in = 0
      DO w=1, nx

        ! Apply shift/squeeze to grid
        x_test = x(w)*p(2) + p(1)
        IF( x_test .GE. self%LUT_DeltaMin .OR. &
          x_test .LE. self%LUT_DeltaMax      ) THEN
          n_in = n_in + 1
          x_in(n_in) = x_test
          i_in(n_in) = w
        ENDIF
          
      ENDDO

      ! Linearly interpolate to wavelength
      CALL LinearInterpolationWeights(self%LUT_wmx,       &
                                      self%LUT_CentralWvl,&
                                      wvl0,idx0,xwt0,xwt1 )
      
      ! Compute weighted ISRF for interpolation
      CALL SPLINT1 (self%LUT_DeltaWvl,                     &
                    self%LUT_ISRF(:,idx0,self%xtrk_idx),   &
                    self%LUT_ISRF_SP(:,idx0,self%xtrk_idx),&
                    self%LUT_dwmx,                         &
                    x_in(1:n_in), Stmp0(1:n_in), n_in,     &
                    SubroutineName                         )
      
      CALL SPLINT1 (self%LUT_DeltaWvl,                       &
                    self%LUT_ISRF(:,idx0,self%xtrk_idx),     &
                    self%LUT_ISRF_SP(:,idx0+1,self%xtrk_idx),&
                    self%LUT_dwmx,                           &
                    x_in(1:n_in), Stmp1(1:n_in), n_in,       &
                    SubroutineName                           )
      
      ! Add to output ISRF
      DO w=1,n_in
        S(i_in(w)) = xwt0*Stmp0(w) + xwt1*Stmp1(w)
      ENDDO
      
      ! ! Now Check the normalization
      ! IF(self%LUT_dwmx .EQ. 2) THEN
      !   NormFac = 0.5d0*dx*SUM(S(1:2))
      ! ELSE
      !   NormFac = 0.5d0*dx*( S(1) + S(nx) + 2.0d0*SUM(S(2:nx-1)) )
      ! ENDIF
      ! print*,NormFac

      ! Initialize Normalization Factor
      NormFac = 1.0d0

      ! Must Renormalize if shifting the wavelength grid
      IF( ABS(p(2)-1.0d0) .GT. TINY(0.0d0) .OR. &
          ABS(p(1)) .GT. TINY(0.0d0)            ) THEN
        STOP 'Need to implemented renormalization before using shift/squeeze in ISRF'
      ENDIF  

      ! Normfac must be recomputed for shift/squeeze grid - For now not worrying
      IF(ScaleBydW) THEN
        NormFac = NormFac / dx
      ENDIF

      ! Apply Normalization
      S(:) = S(:) / NormFac
      
    ELSE

      ! Compute parameterized slit
      CALL self%ParameterizedISRF(nx,x,self%nPar,p,S,ScaleBydW_in=ScaleBydW)

    ENDIF
    
  END SUBROUTINE EvaluateISRF

  ! =======================================================================
  ! =======================================================================
  !                          Convolution Routines
  ! =======================================================================
  ! =======================================================================

  ! Returns Convolution on the L1 Grid for Wavelength dependent parameters
  SUBROUTINE Convolution_WvlDepPar(self,                               &
                                   nhr, nlr, wvl_hr, spc_hr, p,        &
                                   hw1e,wvl_lr, spc_lr, DerivativeIndex)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(ISRF_FunctionType), INTENT(IN)  :: self
    INTEGER,                  INTENT(IN)  :: nhr, nlr
    REAL(KIND=8),             INTENT(IN)  :: wvl_hr(nhr), spc_hr(nhr)
    REAL(KIND=8),             INTENT(IN)  :: p(nlr,self%nPar)
    REAL(KIND=8),             INTENT(IN)  :: hw1e(nlr)
    REAL(KIND=8),             INTENT(IN)  :: wvl_lr(nlr)
    REAL(KIND=8),             INTENT(OUT) :: spc_lr(nlr)
    INTEGER, OPTIONAL,        INTENT(IN)  :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i, iw0, iwf, nw, deriv_idx, errstat, idx0
    REAL(KIND=8) :: dw, wvl0, wvlf, dp, xwt0, xwt1
    REAL(KIND=8) :: S_p(nhr), S_m(nhr),  S_conv(nhr), x(nhr)
    REAL(KIND=8) :: p_pix(self%nPar), p0_pix(self%nPar)
    REAL (KIND=8), PARAMETER :: dhalf = 2.0d0
    REAL(KIND=8)             :: dhalf_var

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'Convolution_WvlDepPar'

    ! =====================================================================
    ! Convolution_WvlDepPar begins here
    ! =====================================================================
    
    ! Get the wavelength increment
    dw = wvl_hr(2)-wvl_hr(1)

    ! Set default derivative index; 0 => Regular Convolution
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex
    
    ! First Check if function has been flagged as a delta impulse
    IF(self%IsDeltaFunction) THEN

      ! CCM Fix later with proper error handling
      ! Raise an error if derivative index is > 0 (not defined)
      IF(deriv_idx > 0) THEN
        STOP 'No derivatives for delta-function'
      ENDIF

      ! Interpolate to the lower resolution Grid
      CALL BSPLINE(wvl_hr, spc_hr, nhr, wvl_lr, spc_lr, nlr, errstat)

      RETURN

    ENDIF
    
    ! =======================================================
    ! Work out the grid perturbation
    ! =======================================================
    dp = 0.0d0
    IF(deriv_idx .GT. 0) THEN

      !==========================================
      ! FIXME
      ! Check the appropriate step sizes later
      ! I suspect its not important
      ! =========================================
      dp = 1e-4

    ENDIF
    
    ! ! ===================================
    ! ! Do the Convolution with the ISRF
    ! ! ===================================
    DO i=1,nlr
        
      ! Convolve out to 10 hw1e
      dhalf_var = self%NumConvolutionHW1E*hw1e(i)

      ! Get range to perform integral
      wvl0 = wvl_lr(i) - dhalf_var
      wvlf = wvl_lr(i) + dhalf_var

      ! Find indices
      iw0 = MAX( FLOOR((wvl0-wvl_hr(1))/dw) ,1)
      iwf = MIN( CEILING((wvlf-wvl_hr(1))/dw),nhr)
      nw  = iwf - iw0 + 1

      ! Compute ISRF on grid
      x(1:nw) = wvl_hr(iw0:iwf) - wvl_lr(i)
      
      ! Set Pixel Parameters
      p0_pix(:) = p(i,:)
      
      IF(deriv_idx .GT. 0) THEN
        
        ! Positive Perturbation
        p_pix = p0_pix ; p_pix(deriv_idx) = p_pix(deriv_idx) + dp
        CALL self%EVAL(wvl_lr(i),nw,x,p_pix,S_p(1:nw))

        ! Negative Perturbation
        p_pix = p0_pix ; p_pix(deriv_idx) = p_pix(deriv_idx) - dp
        CALL self%EVAL(wvl_lr(i),nw,x,p_pix,S_m(1:nw))
          

        ! Compute finite difference derivative
        S_conv = (S_p - S_m) / ( 2.0*dp )

      ELSE
        CALL self%EVAL(wvl_lr(i),nw,x,p0_pix,S_conv(1:nw))
          
      ENDIF
      
      ! Compute integral
      spc_lr(i) = DOT_PRODUCT(S_conv(1:nw),spc_hr(iw0:iwf))
      
    ENDDO

  END SUBROUTINE Convolution_WvlDepPar

  ! Returns convolution on the input grid for fixed parameters
  SUBROUTINE Convolution_FixedPar(self,                          &
                                  nhr, wvl_hr, spc_hr, p,        &
                                  hw1e, spc_conv, DerivativeIndex)

    ! Copies method from OMI/OMPS operational code

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(ISRF_FunctionType), INTENT(IN)  :: self
    INTEGER,                  INTENT(IN)  :: nhr
    REAL(KIND=8),             INTENT(IN)  :: wvl_hr(nhr), spc_hr(nhr)
    REAL(KIND=8),             INTENT(IN)  :: p(self%nPar)
    REAL(KIND=8),             INTENT(IN)  :: hw1e
    REAL(KIND=8),             INTENT(OUT) :: spc_conv(nhr)
    INTEGER, OPTIONAL,        INTENT(IN)  :: DerivativeIndex

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: spc_temp(3*nhr) 
    REAL(KIND=8) :: x(nhr)
    REAL(KIND=8) :: p_pix(self%nPar)
    REAL(KIND=8) :: sf_val(3*nhr)
    REAL(KIND=8) :: sf_val_p(3*nhr),sf_val_m(3*nhr)
    REAL(KIND=8) :: delwvl, dp
    INTEGER      :: sslit, i, nslit, nhalf, deriv_idx

    ! 
    REAL(KIND=8), PARAMETER :: dhalf = 2.0d0
    REAL(KIND=8) :: dhalf_var

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'Convolution_FixedPar'

    ! =====================================================================
    ! Convolution_FixedPar begins here
    ! =====================================================================

    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex

    ! First Check if function has been flagged as a delta impulse
    IF(self%IsDeltaFunction) THEN
      spc_conv(1:nhr) = spc_hr(1:nhr)
      RETURN
    ENDIF

    ! Skip Convolution if hw1e is 0
    IF(hw1e .EQ. 0.0d0 .AND. DerivativeIndex .EQ. 0) THEN
      spc_conv(1:nhr) = spc_hr(1:nhr)
      RETURN
    ENDIF

    ! -----------------------------------------------------------------------
    ! One temporary variable is SPC_TEMP, which is three times the size of
    ! SPEC. For the convolution routine to work in each and every case, we 
    ! reflect the spectrum at its end points to always have a fully filled 
    ! slit function. But this causes some real index headaches when the slit
    ! function wraps around at the ends. Performing the mirror imaging before
    ! we get to the convolution helps to keep things a little more simple.
    ! -----------------------------------------------------------------------
    spc_temp(nhr+1:2*nhr) = spc_hr(1:nhr)
    DO i = 1, nhr
       spc_temp(nhr+1-i) = spc_hr(i)
       spc_temp(2*nhr+i) = spc_hr(nhr+1-i)
    END DO

    ! Convolve out to 10 hw1e
    dhalf_var = self%NumConvolutionHW1E*hw1e ! MAX(10.0d0*hw1e,dhalf)

    ! ---------------------------------------------------------------------------
    ! Because all reference cross sections are provided in regular grids, equally
    ! spaced we only need to work out the super Gaussian slit functin once.
    ! ---------------------------------------------------------------------------
    delwvl = wvl_hr(2) - wvl_hr(1)
    nhalf = CEILING(dhalf_var/delwvl)
    nslit = nhalf*2+1

    ! CCM Fix - Add proper error handling later
    IF(nslit > nhr) THEN
      print*,'Warning: Convolution slit width is greater than window'
      print*,'nslit,nwindow:',nslit,nhr,hw1e,delwvl
      nhalf = FLOOR(REAL(nhr-1,KIND=8)*0.5d0)
      nslit = nhr 
    ENDIF

    ! Compute slit function grid
    DO i = 1, nslit
      x(i) = delwvl * REAL(i-1,KIND=8) - dhalf_var
    ENDDO

    IF(deriv_idx .GT. 0) THEN
      
      ! FIXME - Not sure what a reasonable number is
      dp = 1.0e-5

      ! Positive Perturbation
      p_pix = p ; p_pix(deriv_idx) = p_pix(deriv_idx) + dp
      CALL self%EVAL(0.0d0,nslit,x(1:nslit),p_pix,sf_val_p(1:nslit))

      ! Negative Perturbation
      p_pix = p ; p_pix(deriv_idx) = p_pix(deriv_idx) - dp
      CALL self%EVAL(0.0d0,nslit,x(1:nslit),p_pix,sf_val_m(1:nslit))

      ! Compute finite difference derivative
      sf_val = (sf_val_p - sf_val_m) / ( 2.0*dp )
      
    ELSE

      ! Regular ISRF
      CALL self%EVAL(0.0d0,nslit,x(1:nslit),p,sf_val(1:nslit))
    ENDIF
    
    ! ------------------------
    ! Proceed with convolution
    ! ------------------------
    DO i = 1, nhr

      ! -----------------------------------------------------
      ! Starting index of spectra contributing to convolution
      ! -----------------------------------------------------
      sslit = nhr+1-nhalf+i
      
      ! ----------------
      ! Safe convolution
      ! ----------------
      spc_conv(i) = DOT_PRODUCT(sf_val(1:nslit), spc_temp(sslit:sslit+nslit-1))

    ENDDO

  END SUBROUTINE Convolution_FixedPar

  SUBROUTINE DetermineHW1E(self,p,dw,wmx,hw1e)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(ISRF_FunctionType), INTENT(IN)  :: self
    REAL(KIND=8),             INTENT(IN)  :: p(self%nPar)
    REAL(KIND=8),             INTENT(IN)  :: dw
    INTEGER,                  INTENT(IN)  :: wmx
    REAL(KIND=8),             INTENT(OUT) :: hw1e

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: x(1), S(1), S_thresh, Sp, Sm, x_abs, S_peak, wvl0
    INTEGER      :: iw

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'DetermineHW1E'

    ! =====================================================================
    ! DetermineHW1E Starts here
    ! =====================================================================
    
    IF(self%IsLUT_ISRF) THEN
      
      ! Approximate value
      hw1e = MAXVAL(ABS(self%LUT_DeltaWvl))/self%NumConvolutionHW1E

    ELSE

      ! THis is not used in the parameterizations but is part of the EVAL interface
      wvl0 = 0.0

      ! Compute Value at centre
      x(1) = 0.0d0 ; CALL self%EVAL(wvl0,1,x,p,S,ScaleBydW_in=.FALSE.) ; S_peak = ABS(S(1))
      
      ! Threshold Value
      S_thresh = S_peak*EXP(-1.0d0)

      ! Search to find minimum
      iw = 0 ; hw1e = -1.0d0 ; x_abs = 0.0d0

      DO iw=1,wmx
        
        ! Increment Counter
        x_abs = x_abs + dw 

        ! Postive Case
        x(1) = x_abs
        CALL self%EVAL(wvl0,1,x,p,S,ScaleBydW_in=.FALSE.) ; Sp = S(1)
        IF(ABS(Sp) .GT. S_peak) THEN
          S_peak = ABS(Sp) ; S_thresh = S_peak*EXP(-1.0d0)
        ENDIF

        ! Negative Case
        x(1) = -1.0d0*x_abs
        CALL self%EVAL(wvl0,1,x,p,S,ScaleBydW_in=.FALSE.) ; Sm = S(1)
        IF(ABS(Sm) .GT. S_peak) THEN
          S_peak = ABS(Sm) ; S_thresh = S_peak*EXP(-1.0d0)
        ENDIF
        
        ! Check if we have found hw1e (accounting for either side of slit)
        IF(ABS(Sp) .LT. S_thresh .AND. ABS(Sm) .LT. S_thresh) THEN
          hw1e = x_abs
          EXIT
        ENDIF
        
      ENDDO

    ENDIF

  END SUBROUTINE DetermineHW1E

  ! =======================================================================
  ! =======================================================================
  !                  Utility Routines for various ISRFs
  ! =======================================================================
  ! =======================================================================

  ! Supergaussian sign func.
  REAL(KIND=8) FUNCTION signdp ( x )
    
    REAL(KIND=8), INTENT (IN) :: x

    signdp = 0.0d0
    IF ( x < 0.0d0 ) THEN
      signdp = -1.0d0
    ELSE
      signdp = +1.0d0
    ENDIF

    RETURN
    
  END FUNCTION signdp

! ============================================================================
! The following are for computing Owen's T function needed by the TROPOMI ISRF
! ============================================================================

REAL(KIND=8) FUNCTION alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none
  
  ! Arguments
  real ( kind = 8 ), intent(in) :: x
  logical,           intent(in) :: upper

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
END FUNCTION alnorm

subroutine owen_values ( n_data, h, a, t )

!*****************************************************************************80
!
!! OWEN_VALUES returns some values of Owen's T function.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        Integral ( 0 <= X <= A ) e^(-H^2*(1+X^2)/2) / (1+X^2) dX
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) H, a parameter.
!
!    Output, real ( kind = 8 ) A, the upper limit of the integral.
!
!    Output, real ( kind = 8 ) T, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.2500000000000000D+00, &
    0.4375000000000000D+00, &
    0.9687500000000000D+00, &
    0.0625000000000000D+00, &
    0.5000000000000000D+00, &
    0.9999975000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03 /)
  real ( kind = 8 ) h
  real ( kind = 8 ), save, dimension ( n_max ) :: h_vec = (/ &
    0.0625000000000000D+00, &
    6.5000000000000000D+00, &
    7.0000000000000000D+00, &
    4.7812500000000000D+00, &
    2.0000000000000000D+00, &
    1.0000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( n_max ) :: t_vec = (/ &
    3.8911930234701366D-02, &
    2.0005773048508315D-11, &
    6.3990627193898685D-13, &
    1.0632974804687463D-07, &
    8.6250779855215071D-03, &
    6.6741808978228592D-02, &
    0.4306469112078537D-01, &
    0.6674188216570097D-01, &
    0.7846818699308410D-01, &
    0.7929950474887259D-01, &
    0.6448860284750376D-01, &
    0.1066710629614485D+00, &
    0.1415806036539784D+00, &
    0.1510840430760184D+00, &
    0.7134663382271778D-01, &
    0.1201285306350883D+00, &
    0.1666128410939293D+00, &
    0.1847501847929859D+00, &
    0.7317273327500385D-01, &
    0.1237630544953746D+00, &
    0.1737438887583106D+00, &
    0.1951190307092811D+00, &
    0.7378938035365546D-01, &
    0.1249951430754052D+00, &
    0.1761984774738108D+00, &
    0.1987772386442824D+00, &
    0.2340886964802671D+00, &
    0.2479460829231492D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    h = 0.0D+00
    a = 0.0D+00
    t = 0.0D+00
  else
    h = h_vec(n_data)
    a = a_vec(n_data)
    t = t_vec(n_data)
  end if

  return
end subroutine owen_values

real(kind=8) function tfn ( x, fx )

!*****************************************************************************80
!
!! TFN calculates the T-function of Owen.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by JC Young, Christoph Minder.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    MA Porter, DJ Winstanley,
!    Remark AS R30:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, page 113.
!
!    JC Young, Christoph Minder,
!    Algorithm AS 76:
!    An Algorithm Useful in Calculating Non-Central T and
!    Bivariate Normal Distributions,
!    Applied Statistics,
!    Volume 23, Number 3, 1974, pages 455-457.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, FX, the parameters of the function.
!
!    Output, real ( kind = 8 ) TFN, the value of the T-function.
!
  implicit none

  integer ( kind = 4 ), parameter :: ng = 5

  real ( kind = 8 ) fx
  real ( kind = 8 ) fxs
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( ng ) :: r = (/ &
    0.1477621D+00, &
    0.1346334D+00, &
    0.1095432D+00, &
    0.0747257D+00, &
    0.0333357D+00 /)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rt
  !real ( kind = 8 ) tfn
  real ( kind = 8 ), parameter :: tp = 0.159155D+00
  real ( kind = 8 ), parameter :: tv1 = 1.0D-35
  real ( kind = 8 ), parameter :: tv2 = 15.0D+00
  real ( kind = 8 ), parameter :: tv3 = 15.0D+00
  real ( kind = 8 ), parameter :: tv4 = 1.0D-05
  real ( kind = 8 ), dimension ( ng ) :: u = (/ &
    0.0744372D+00, &
    0.2166977D+00, &
    0.3397048D+00, &
    0.4325317D+00, &
    0.4869533D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xs
!
!  Test for X near zero.
!
  if ( abs ( x ) < tv1 ) then
    tfn = tp * atan ( fx )
    return
  end if
!
!  Test for large values of abs(X).
!
  if ( tv2 < abs ( x ) ) then
    tfn = 0.0D+00
    return
  end if
!
!  Test for FX near zero.
!
  if ( abs ( fx ) < tv1 ) then
    tfn = 0.0D+00
    return
  end if
!
!  Test whether abs ( FX ) is so large that it must be truncated.
!
  xs = - 0.5D+00 * x * x
  x2 = fx
  fxs = fx * fx
!
!  Computation of truncation point by Newton iteration.
!
  if ( tv3 <= log ( 1.0D+00 + fxs ) - xs * fxs ) then

    x1 = 0.5D+00 * fx
    fxs = 0.25D+00 * fxs

    do

      rt = fxs + 1.0D+00

      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) &
      / ( 2.0D+00 * x1 * ( 1.0D+00 / rt - xs ) )

      fxs = x2 * x2

      if ( abs ( x2 - x1 ) < tv4 ) then
        exit
      end if

      x1 = x2

    end do

  end if
!
!  Gaussian quadrature.
!
  rt = 0.0D+00

  do i = 1, ng

    r1 = 1.0D+00 + fxs * ( 0.5D+00 + u(i) )**2
    r2 = 1.0D+00 + fxs * ( 0.5D+00 - u(i) )**2

    rt = rt + r(i) * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 )

  end do

  tfn = rt * x2 * tp

  return
end function tfn

function tha ( h1, h2, a1, a2 )

!*****************************************************************************80
!
!! THA computes Owen's T function.
!
!  Discussion:
!
!    This function computes T(H1/H2, A1/A2) for any real numbers H1, H2,
!    A1 and A2.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by JC Young, Christoph Minder.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Boys,
!    Remark AS R80:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 38, Number 3, 1989, pages 580-582.
!
!    Youn-Min Chou,
!    Remark AS R55:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 34, Number 1, 1985, pages 100-101.
!
!    PW Goedhart, MJW Jansen,
!    Remark AS R89:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 41, Number 2, 1992, pages 496-497.
!
!    JC Young, Christoph Minder,
!    Algorithm AS 76:
!    An Algorithm Useful in Calculating Noncentral T and
!    Bivariate Normal Distributions,
!    Applied Statistics,
!    Volume 23, Number 3, 1974, pages 455-457.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H1, H2, A1, A2, define the arguments
!    of the T function.
!
!    Output, real ( kind = 8 ) THA, the value of Owen's T function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) absa
  real ( kind = 8 ) ah
  !real ( kind = 8 ) alnorm
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) ex
  real ( kind = 8 ) g
  real ( kind = 8 ) gah
  real ( kind = 8 ) gh
  real ( kind = 8 ) h
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) lam
  !real ( kind = 8 ) tfn
  real ( kind = 8 ) tha
  real ( kind = 8 ), parameter :: twopi = 6.2831853071795864769D+00

  if ( h2 == 0.0D+00 ) then
    tha = 0.0D+00
    return
  end if

  h = h1 / h2

  if ( a2 == 0.0D+00 ) then

    g = alnorm ( h, .false. )

    if ( h < 0.0D+00 ) then
      tha = g / 2.0D+00
    else
      tha = ( 1.0D+00 - g ) / 2.0D+00
    end if

    if ( a1 < 0.0D+00 ) then
      tha = - tha
    end if

    return
  end if

  a = a1 / a2

  if ( abs ( h ) < 0.3D+00 .and. 7.0D+00 < abs ( a ) ) then

    lam = abs ( a * h )
    ex = exp ( - lam * lam / 2.0D+00 )
    g = alnorm ( lam, .false. )
    c1 = ( ex / lam + sqrt ( twopi ) * ( g - 0.5D+00 ) ) / twopi
    c2 = ( ( lam * lam + 2.0D+00 ) * ex / lam**3 &
    + sqrt ( twopi ) * ( g - 0.5D+00 ) ) / ( 6.0D+00 * twopi )
    ah = abs ( h )
    tha = 0.25D+00 - c1 * ah + c2 * ah**3
    tha = sign ( tha, a )

  else
!
!  Correction AS R89
!
    absa = abs ( a )

    if ( absa <= 1.0D+00 ) then
      tha = tfn ( h, a )
      return
    end if

    ah = absa * h
    gh = alnorm ( h, .false. )
    gah = alnorm ( ah, .false. )
    tha = 0.5D+00 * ( gh + gah ) - gh * gah &
    - tfn ( ah, 1.0D+00 / absa )

    if ( a < 0.0D+00 ) then
      tha = - tha
    end if

  end if

  return
end function tha 

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp

! ====================================================================================================
! ====================================================================================================
! ====================================================================================================
! OLD ROUTINES
! ====================================================================================================
! ====================================================================================================
! ====================================================================================================

SUBROUTINE super_gaussian_sf(npoints, hw1e, e_asym, g_shap, wvlarr, specarr, specmod, DerivativeIndex)

    ! =========================================================================
    !
    ! Convolves input spectrum with an asymmetric Gaussian slit function of
    ! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
    !
    ! The asymetric Gaussian g(x) is defined as
    !                   _                                            _
    !                  |   |            x^2                  |^g_shap |
    !      g(x) =  EXP | - |---------------------------------|        |
    !                  |_  | (hw1e * (1 + SIGN(x)*e_asym))   |       _|
    !
    ! g(x) becomes symmetric for E_ASYM = 0.
    !
    ! =========================================================================


    IMPLICIT NONE

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,           INTENT (IN)  :: npoints
    REAL(KIND=8),      INTENT (IN)  :: hw1e, e_asym, g_shap
    REAL(KIND=8),      INTENT (IN)  :: wvlarr(npoints), specarr(npoints)
    REAL (KIND=8),     INTENT (OUT) :: specmod(npoints)
    INTEGER, OPTIONAL, INTENT(IN)   :: DerivativeIndex

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=8)  :: delwvl, dp
    INTEGER        :: i, nslit, sslit, nhalf, deriv_idx
    REAL (KIND=8)  :: slitsum, wvl
    REAL (KIND=8)  :: slitsum_p, slitsum_m
    REAL (KIND=8)  :: spc_temp(3*npoints) 
    REAL (KIND=8)  :: sf_val(3*npoints)
    REAL (KIND=8)  :: sf_val_p(3*npoints),sf_val_m(3*npoints)
!     REAL (KIND=8)  :: signdp
    
    
    REAL (KIND=8), PARAMETER :: dhalf = 2.0d0
    REAL(KIND=8) :: dhalf_var
    
!     EXTERNAL signdp
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'super_gaussian_sf'

    ! --------------------------------------------------------
    ! Initialize output variable (default for "no convolution"
    ! --------------------------------------------------------
    specmod(1:npoints) = specarr(1:npoints)
    
    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex

    ! -----------------------------------------------
    ! No Gaussian convolution if Halfwidth @ 1/e is 0
    ! -----------------------------------------------
    IF ( hw1e == 0.0d0 .AND. deriv_idx .EQ. 0) RETURN

    ! -----------------------------------------------------------------------
    ! One temporary variable is SPC_TEMP, which is three times the size of
    ! SPEC. For the convolution routine to work in each and every case, we 
    ! reflect the spectrum at its end points to always have a fully filled 
    ! slit function. But this causes some real index headaches when the slit
    ! function wraps around at the ends. Performing the mirror imaging before
    ! we get to the convolution helps to keep things a little more simple.
    ! -----------------------------------------------------------------------
    spc_temp(npoints+1:2*npoints) = specarr(1:npoints)
    DO i = 1, npoints
       spc_temp(npoints+1-i) = specarr(i)
       spc_temp(2*npoints+i) = specarr(npoints+1-i)
    END DO
    
    ! Convolve out to 10 hw1e
    dhalf_var = MAX(10.0d0*hw1e,dhalf)
    
    ! ---------------------------------------------------------------------------
    ! Because all reference cross sections are provided in regular grids, equally
    ! spaced we only need to work out the super Gaussian slit functin once.
    ! ---------------------------------------------------------------------------
    delwvl = wvlarr(2) - wvlarr(1)
    nhalf = CEILING (dhalf_var/delwvl)
    nslit = nhalf*2+1
    IF(deriv_idx .GT. 0) THEN
      dp = 1.0d-3
      IF(deriv_idx .EQ. 1) THEN
        CALL getslit_sgf(delwvl,dhalf_var,hw1e+dp,e_asym,g_shap,nslit,sf_val_p(1:nslit),slitsum_p)
        CALL getslit_sgf(delwvl,dhalf_var,hw1e-dp,e_asym,g_shap,nslit,sf_val_m(1:nslit),slitsum_m)
      ELSEIF(deriv_idx .EQ. 2) THEN
        CALL getslit_sgf(delwvl,dhalf_var,hw1e,e_asym+dp,g_shap,nslit,sf_val_p(1:nslit),slitsum_p)
        CALL getslit_sgf(delwvl,dhalf_var,hw1e,e_asym-dp,g_shap,nslit,sf_val_m(1:nslit),slitsum_m)
      ELSEIF(deriv_idx .EQ. 3) THEN
        CALL getslit_sgf(delwvl,dhalf_var,hw1e,e_asym,g_shap+dp,nslit,sf_val_p(1:nslit),slitsum_p)
        CALL getslit_sgf(delwvl,dhalf_var,hw1e,e_asym,g_shap-dp,nslit,sf_val_m(1:nslit),slitsum_m)
      ELSE
        print*,'Unknown derivative Index (super_gaussian_sf)',deriv_idx ! RaiseFatalError
        STOP 1
      ENDIF
      
      ! To complete the derivative
      sf_val(1:nslit) = sf_val_p(1:nslit)/slitsum_p - sf_val_m(1:nslit)/slitsum_m
      slitsum = 2.0d0*dp
      
    ELSE
      CALL getslit_sgf(delwvl,dhalf_var,hw1e,e_asym,g_shap,nslit,sf_val(1:nslit),slitsum)
    ENDIF
    ! ------------------------
    ! Proceed with convolution
    ! ------------------------
    DO i = 1, npoints
       ! -----------------------------------------------------
       ! Starting index of spectra contributing to convolution
       ! -----------------------------------------------------
       sslit = npoints+1-nhalf+i
       
       ! ----------------
       ! Safe convolution
       ! ----------------
       specmod(i) = DOT_PRODUCT(sf_val(1:nslit), spc_temp(sslit:sslit+nslit-1))/slitsum

    END DO
    
    RETURN
    
  END SUBROUTINE super_gaussian_sf
  
  ! For convenience
  SUBROUTINE getslit_sgf(delwvl,dhalf_var,hw1e,e_asym,g_shap,nslit,sf_val,slitsum)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL(KIND=8), INTENT(IN)  :: delwvl,dhalf_var,hw1e,e_asym,g_shap
    INTEGER,      INTENT(IN)  :: nslit
    REAL(KIND=8), INTENT(OUT) :: sf_val(nslit),slitsum

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: wvl
    INTEGER      :: i

    ! =====================================================================
    ! getslit_sgf begins here
    ! =====================================================================

    DO i = 1, nslit
      wvl = delwvl * REAL(i-1,KIND=8) - dhalf_var
      sf_val(i) = EXP(-(ABS(wvl / ( hw1e + signdp(wvl)*e_asym ) ) )**g_shap )
    END DO
    slitsum = SUM(sf_val(1:nslit))

  END SUBROUTINE getslit_sgf

  SUBROUTINE SuperGauss_Func(nx,x,hw1e,e_asym,g_shap,S,ScaleBydW_in)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER ,          INTENT(IN)  :: nx 
    REAL(KIND=8),      INTENT(IN)  :: x(nx) ! Assumes even spacing
    REAL(KIND=8),      INTENT(IN)  :: hw1e,e_asym,g_shap
    REAL(KIND=8),      INTENT(OUT) :: S(nx)
    LOGICAL, OPTIONAL, INTENT(IN)  :: ScaleBydW_in ! Default=TRUE

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: i
    LOGICAL :: ScaleBydW
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SuperGauss_Func'

    ! =====================================================================
    ! SuperGauss_Func begins here
    ! =====================================================================

    ! Default Scaling 
    ScaleBydW = .TRUE. ; IF(PRESENT(ScaleBydW_in)) ScaleBydW = ScaleBydW_in


    DO i=1,nx
      S(i) = EXP( -(ABS(x(i)/(hw1e + signdp(x(i))*e_asym)) )**g_shap )
    ENDDO

    ! Normalize
    IF(ScaleBydW) THEN

      ! Return S(x)*dx
      S = S / SUM(S)

    ELSE

      ! Return S(x)
      S = S / SUM(S) / (x(2)-x(1))

    ENDIF
  END SUBROUTINE SuperGauss_Func

  SUBROUTINE SuperGaussISRF_WvlDep( nhr, nlr, wvl_hr, spc_hr, hw1e, e_asym, g_shap,&
                                    wvl_lr, spc_lr,DerivativeIndex                 )

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,          INTENT(IN)  :: nhr, nlr
    REAL(KIND=8),     INTENT(IN)  :: wvl_hr(nhr), spc_hr(nhr)
    REAL(KIND=8),     INTENT(IN)  :: hw1e(nlr), e_asym(nlr), g_shap(nlr)
    REAL(KIND=8),     INTENT(IN)  :: wvl_lr(nlr)
    REAL(KIND=8),     INTENT(OUT) :: spc_lr(nlr)
    INTEGER, OPTIONAL,INTENT(IN)  :: DerivativeIndex
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i, iw0, iwf, nw, deriv_idx
    REAL(KIND=8) :: dw, wvl0, wvlf
    REAL(KIND=8) :: dhw1e, dasym, dshap, dp
    REAL(KIND=8) :: S_p(nhr), S_m(nhr),  S_conv(nhr), x(nhr)

    REAL (KIND=8), PARAMETER :: dhalf = 2.0d0
    REAL(KIND=8)             :: dhalf_var

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SuperGaussISRF_WvlDep'

    ! =====================================================================
    ! SuperGaussISRF_WvlDep begins here
    ! =====================================================================

    ! Get the wavelength increment
    dw = wvl_hr(2)-wvl_hr(1)

    ! Set default derivative index
    deriv_idx = 0 ; IF(PRESENT(DerivativeIndex)) deriv_idx = DerivativeIndex
    
    ! =======================================================
    ! Work out the grid perturbation
    ! =======================================================
    IF(deriv_idx .GT. 0) THEN

      ! Zero the perturbations
      dhw1e = 0.0d0 ;  dasym = 0.0d0 ; dshap = 0.0d0

      !==========================================
      ! FIXME
      ! Check the appropriate step sizes later
      ! I suspect its not important
      ! =========================================
      IF(deriv_idx .EQ. 1) THEN
        dhw1e = 1.0d-4 ; dp = dhw1e
      ELSEIF(deriv_idx .EQ. 2) THEN
        dasym = 1.0d-4 ; dp = dasym
      ELSEIF(deriv_idx .EQ. 3) THEN
        dshap = 1.0d-4 ; dp = dshap
      ELSE
        print*,'Unknown derivative Index (SuperGaussISRF_WvlDep)',deriv_idx ! RaiseFatalError
        STOP 1
      ENDIF
    
    ENDIF

    ! ===================================
    ! Do the Convolution with the ISRF
    ! ===================================
    DO i=1,nlr
      
      ! Convolve out to 10 hw1e
      dhalf_var = MAX(10.0d0*hw1e(i),dhalf)

      ! Get range to perform integral
      wvl0 = wvl_lr(i) - dhalf_var
      wvlf = wvl_lr(i) + dhalf_var

      ! Find indices
      iw0 = MAX( FLOOR((wvl0-wvl_hr(1))/dw) ,1)
      iwf = MIN( CEILING((wvlf-wvl_hr(1))/dw),nhr)
      nw  = iwf - iw0 + 1

      ! Compute Supergaussian on grid
      x(1:nw) = wvl_hr(iw0:iwf) - wvl_lr(i)

      IF(deriv_idx .GT. 0) THEN
        CALL SuperGauss_Func(nw,x(1:nw),hw1e(i)+dhw1e,e_asym(i)+dasym,g_shap(i)+dshap,S_p(1:nw))
        CALL SuperGauss_Func(nw,x(1:nw),hw1e(i)-dhw1e,e_asym(i)-dasym,g_shap(i)-dshap,S_m(1:nw))
        S_conv = (S_p - S_m) / ( 2.0*dp )
      ELSE
        CALL SuperGauss_Func(nw,x(1:nw),hw1e(i),e_asym(i),g_shap(i),S_conv(1:nw))
      ENDIF

      ! Compute integral
      spc_lr(i) = DOT_PRODUCT(S_conv(1:nw),spc_hr(iw0:iwf))

    ENDDO

  END SUBROUTINE SuperGaussISRF_WvlDep




END MODULE isrf_module
