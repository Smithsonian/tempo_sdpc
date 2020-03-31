MODULE sif_module
  
  USE parameters_module
  USE error_module,         ONLY : ErrorType, CheckError, RaiseFatalError
  USE netcdf_module,        ONLY : CheckNetCDFErrorStatus
  USE interpolation_module, ONLY : SPLINE1, SPLINT1

  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'

  TYPE SIFOptType
    CHARACTER(LEN=maxChar) :: NormalizedSpecInfile
  ENDTYPE SIFOptType

  TYPE SIFType
    TYPE(SIFOptType)          :: Opt
    INTEGER(KIND=2)           :: TypeIndex
    INTEGER                   :: nNormWvl
    REAL(KIND=8), ALLOCATABLE :: NormSpecWvl(:)
    REAL(KIND=8)              :: NormSpecMinWvl
    REAL(KIND=8)              :: NormSpecMaxWvl
    REAL(KIND=8), ALLOCATABLE :: NormSpecSIF(:)
    REAL(KIND=8), ALLOCATABLE :: NormSpecSIFSP(:)

    INTEGER                   :: nSpecWvl
    REAL(KIND=8), ALLOCATABLE :: SIFSpectrum(:)
    ! Functions
    CONTAINS
      PROCEDURE :: Initialize      => InitSIF
      PROCEDURE :: ComputeSpectrum => ComputeSIF

  ENDTYPE SIFType
  
  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'sif_module'
  PRIVATE :: ModuleName

  CONTAINS

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: InitSIF
  ! 
  ! DESCRIPTION: Initialization of solar induced fluorescence 

  SUBROUTINE InitSIF(self,SIFOpt,Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(SIFType),   INTENT(INOUT) :: self
    TYPE(SIFOptType), INTENT(IN)    :: SIFOpt
    TYPE(ErrorType),  INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, ncid, vid, dimid
    CHARACTER(LEN=maxChar) :: tmpchar

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitSIF'

    ! =====================================================================
    ! InitSIF begins here
    ! =====================================================================

    IF(CheckError(Error)) RETURN

    ! Save SIF Options
    self%Opt = SIFOpt

    ! Open Normalized spectrum file
    rcode = nf_open(TRIM(ADJUSTL(SIFOpt%NormalizedSpecInfile)), NF_SHARE, ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                'nf_open:'//TRIM(ADJUSTL(SIFOpt%NormalizedSpecInfile)))

    ! Read Type
    rcode = nf_inq_varid( ncid, 'SIF_Type_Index', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_var:SIF_Type_Index')
    rcode = nf_get_vara_int2(ncid, vid, (/1/), (/1/), self%TypeIndex)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    
    ! (1) Single Normalized Spectrum
    IF(self%TypeIndex .EQ. 1) THEN

      ! Get dimension
      rcode = nf_inq_varid(ncid, 'Wavelength', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_var:Wavelength')
      rcode = nf_inq_vardimid(ncid,vid,dimid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_vardim:Wavelength')
      rcode = NF_INQ_DIM(ncid,dimid,tmpchar,self%nNormWvl)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_dim:Wavelength')

      ! Allocate arrays
      ALLOCATE(self%NormSpecWvl(self%nNormWvl))
      ALLOCATE(self%NormSpecSIF(self%nNormWvl))
      ALLOCATE(self%NormSpecSIFSP(self%nNormWvl))

      ! Load Spectrum
      rcode = nf_get_var_double(ncid, vid, self%NormSpecWvl)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_get_var:Wavelength')
      rcode = nf_inq_varid(ncid, 'SIF_Radiance', vid)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_inq_var:SIF_Radiance')
      rcode = nf_get_var_double(ncid, vid, self%NormSpecSIF)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'nf_get_var:SIF_Radiance')

      ! Store max/min wvl
      self%NormSpecMinWvl = MINVAL(self%NormSpecWvl)
      self%NormSpecMaxWvl = MAXVAL(self%NormSpecWvl)

      ! Compute basis spline coefficients
      CALL SPLINE1(self%NormSpecWvl,self%NormSpecSIF,self%nNormWvl,self%NormSpecSIFSP)
      
    ELSE
      CALL RaiseFatalError( Error, ErrorCode_Profile, ModuleName, SubroutineName,                 &
                            Message_in='Urecognized SIF Type',                                    &
                            Action_in='Check sif file'//TRIM(ADJUSTL(SIFOpt%NormalizedSpecInfile)))
    ENDIF

    ! Intiialize interpolated spectrum dimension
    self%nSpecWvl = 0

  END SUBROUTINE InitSIF
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: ComputeSIF
  ! 
  ! DESCRIPTION: Returns plant induced fluorescense for reference SIF

  SUBROUTINE ComputeSIF(self,nWvl,Wvl,SIF_734nm)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(SIFType),   INTENT(INOUT) :: self
    INTEGER,          INTENT(IN)    :: nWvl
    REAL(KIND=8),     INTENT(IN)    :: Wvl(nWvl)
    REAL(KIND=8),     INTENT(IN)    :: SIF_734nm
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: nWvl_int, N
    REAL(KIND=8) :: unit_conv
    REAL(KIND=8) :: Wvl_int(nWvl), Spc_int(nWvl)
    INTEGER      :: Idx_int(nWvl)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ComputeSIF'

    ! =====================================================================
    ! ComputeSIF begins here
    ! =====================================================================
    
    ! Unit Conversion (mW/m2/Sr/nm -> photons/cm2/nm/s/Sr)
    unit_conv = 1.0d-16 * 734.0d0 / ( 2.99792458d8*6.62607004d-34 )

    ! Check SIF Spectrum allocation
    IF(self%nSpecWvl .NE. nWvl) THEN
      IF(ALLOCATED(self%SIFSpectrum)) DEALLOCATE(self%SIFSpectrum)
      ALLOCATE(self%SIFSpectrum(nWvl))
      self%nSpecWvl = nWvl
    ENDIF

    ! Zero SIF spectrum
    self%SIFSpectrum(:) = 0.0d0

    ! (1) Single Normalized Spectrum
    IF(self%TypeIndex .EQ. 1) THEN

      ! Find wavelength within range
      nWvl_int = 0
      DO N=1,nWvl
        IF(Wvl(N) .GE. self%NormSpecMinWvl .AND. &
           Wvl(N) .LE. self%NormSpecMaxWvl       ) THEN
          nWvl_int = nWvl_int + 1
          Wvl_int(nWvl_int) = Wvl(N)
          Idx_int(nWvl_int) = N
        ENDIF
      ENDDO

      
      IF(nWvl_int .GT. 0) THEN

        ! Do the interpolation
        CALL SPLINT1(self%NormSpecWvl, self%NormSpecSIF,               &
                     self%NormSpecSIFSP, self%nNormWvl,                &
                     Wvl_int(1:nWvl_int), Spc_int(1:nWvl_int), nWvl_int)

        ! Update SIF spectrum (photons/cm2/nm/s/Sr)
        DO N=1,nWvl_int
          self%SIFSpectrum(Idx_int(N)) = SIF_734nm*Spc_int(N)*unit_conv*Wvl_int(N)/734.0d0
        ENDDO

      ENDIF

      ! Apply unit conversion
    ENDIF

  END SUBROUTINE ComputeSIF
  

END MODULE sif_module