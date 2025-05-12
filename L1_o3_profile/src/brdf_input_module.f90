MODULE input_module

  USE parameters_module
  USE error_module,        ONLY : ErrorType, RaiseFatalError
  USE string_utils_module, ONLY : SPLIT_ONE_LINE, READ_ONE_LINE,&
                                  STRREPL, SkipFileLines
  USE surface_module,      ONLY : SurfOptType
  IMPLICIT NONE

  ! =====================================================================
  ! Module Variables
  ! =====================================================================
  
  
  TYPE InpOptType
    INTEGER                             :: funit
    CHARACTER(LEN=maxChar)              :: CalculationMode
    CHARACTER(LEN=maxChar)              :: l2_outfile
    TYPE(SurfOptType)                   :: Surface
  ENDTYPE InpOptType
  
  PUBLIC  :: ReadInputFile
             
  PRIVATE :: READ_SURFACE_OPTIONS
  
  CONTAINS
    
    
    !###################################################################
    !#                              SPLAT                              #
    !###################################################################
    
    ! SUBROUTINE: ReadInputFile
    ! 
    ! DESCRIPTION: Reads the entire SPLAT input file. Transfers common options
    !              that need to be consistent between modules to their input 
    !              structures. Also determines the jacobian flags that need to 
    !              be set in the forward model based on the input options.
    
    SUBROUTINE ReadInputFile( input_filename, SurfOpt, error )

      ! Subroutine Arguments
      CHARACTER(LEN=*), INTENT(IN) :: input_filename
      TYPE(SurfOptType), INTENT(INOUT)    :: SurfOpt
      TYPE(ErrorType),  INTENT(INOUT)    :: error

      ! ---------------
      ! Local Variables
      ! ---------------
      TYPE(InpOptType) :: InpOpt
      CHARACTER(LEN=1)       :: TAB   = ACHAR(9)
      CHARACTER(LEN=1)       :: SPACE = ' '
      CHARACTER(LEN=maxChar) :: LINE
      INTEGER                :: ios
      LOGICAL                :: EOF
      
      ! =====================================================================
      ! readInputFile starts here 
      ! =====================================================================
      
      ! Open control file
      InpOpt%funit=ctrunit
      OPEN (UNIT=ctrunit, FILE=TRIM(ADJUSTL(input_filename)), &
            ACTION='READ', STATUS='OLD', IOSTAT=ios                )
      IF(ios /= 0) THEN
        WRITE(*,*) 'BRDF: cannot read'//TRIM(ADJUSTL(input_filename))
      ENDIF
      
      ! Loop over files
      DO
        
        ! Read a line from the file, exit if EOF
        LINE = READ_ONE_LINE( ctrunit, EOF )
        IF ( EOF ) EXIT
          
        ! Replace tab characters in LINE (if any) w/ spaces
        CALL STRREPL( LINE, TAB, SPACE )
          
        IF(INDEX( LINE, '%%% SURFACE OPTIONS %%%'     ) > 0 ) THEN
          CALL READ_SURFACE_OPTIONS(InpOpt, Error)
        ELSEIF(INDEX( LINE, 'END OF FILE'                 ) > 0 ) THEN 
          EXIT
        ENDIF
        
      ENDDO
      
      ! Close input file
      CLOSE( InpOpt%funit )
      SurfOpt = InpOpt%Surface
      !InpOpt%Surface%DoAmplitudeLinearization = InpOpt%Diag%RTM%BRDFAmplitudeJacobian%DoDiag
      !InpOpt%Surface%DoParameterLinearization = InpOpt%Diag%RTM%BRDFParameterJacobian%DoDiag

      END SUBROUTINE ReadInputFile
      
       SUBROUTINE READ_SURFACE_OPTIONS(InpOpt, Error)

      ! --------------------
      ! Subroutine Arguments
      ! -------------------- 
      TYPE(InpOptType), INTENT(INOUT) :: InpOpt
      TYPE(ErrorType),  INTENT(INOUT) :: error

      ! ---------------
      ! Local Variables
      ! ---------------
      INTEGER                :: N, I, J
      CHARACTER(LEN=maxChar) :: SUBSTRS(maxChar)

      ! =====================================================================
      ! READ_SURFACE_OPTIONS starts here 
      ! =====================================================================

      ! Initialize linearization flags
      InpOpt%Surface%DoAmplitudeLinearization = .FALSE.
      InpOpt%Surface%DoParameterLinearization = .FALSE.

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%OptionIndex
      ! Skip line
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )
      
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%Option1_FixedAlbedo

      ! Skip line
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), '(A)' ) InpOpt%Surface%Option2_Infile

      ! Skip line
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), '(A)' ) InpOpt%Surface%Option3_Infile

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%Option3_UseConstantWvl

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%Option3_ConstWvl

      ! Skip line
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), '(A)' ) InpOpt%Surface%Option4_Infile

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), '(A)' ) InpOpt%Surface%Option4_ClimDir

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%Option4_DoIsotropic

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%Option4_WhichAlbedo

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%Option4_DoOceanGlint

      ! Skip 2 lines
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )

      ! Zero Kernel Parameters
      InpOpt%Surface%Option5_KernPar(:,:) = 0.0d0
       
      ! Read Fixed Kernels
      DO I=1,3
         CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )
         READ(SUBSTRS(1),'(A)') InpOpt%Surface%Option5_KernName(I)
         READ(SUBSTRS(2),*) InpOpt%Surface%Option5_KernIdx(I)
         READ(SUBSTRS(3),*) InpOpt%Surface%Option5_KernAmp(I)
         READ(SUBSTRS(4),*) InpOpt%Surface%Option5_nKernPar(I)
         DO J=1,InpOpt%Surface%Option5_nKernPar(I)
           READ(SUBSTRS(4+J),*) InpOpt%Surface%Option5_KernPar(J,I)
         ENDDO
      ENDDO

      ! Skip line
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )

      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, -1, 'read_surface_options:1' )
      READ( SUBSTRS(1), '(A)' ) InpOpt%Surface%Option6_Infile

      ! Read Emissivity Options
      ! -----------------------------------------------------------------------
      
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1), * ) InpOpt%Surface%EmissivityOptIndex
      
      ! Option 1 - Constant
      CALL SkipFileLines(InpOpt%funit,1)
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1), * ) InpOpt%Surface%ConstantEmissivity
      
      ! Option 2 - Emissivity Spectrum
      CALL SkipFileLines(InpOpt%funit,1)
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      IF(InpOpt%Surface%EmissivityOptIndex .EQ. 2) &
        READ( SUBSTRS(1),'(A)')  InpOpt%Surface%EmissivityInfile
      
      ! Option 3 - Emissivity Climatology
      CALL SkipFileLines(InpOpt%funit,1)
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      IF(InpOpt%Surface%EmissivityOptIndex .EQ. 3) &
        READ( SUBSTRS(1),'(A)')  InpOpt%Surface%EmissivityInfile
        
      ! Stope code for now if we are using features not yet implemented
      IF(InpOpt%Surface%EmissivityOptIndex > 1) THEN
        STOP 'Emissivity options > 1 have not yet been implemented :( '
      ENDIF
      
      ! Read Options for additional surface sources
      ! -----------------------------------------------------------------------
      
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1:N), * ) InpOpt%Surface%DoPlantFluorescence

      RETURN
      CALL SPLIT_ONE_LINE( InpOpt%funit, SUBSTRS, N, 1, 'read_surface_options:1' )
      READ( SUBSTRS(1),'(A)')  InpOpt%Surface%SIF%NormalizedSpecInfile

    
    END SUBROUTINE READ_SURFACE_OPTIONS
  
    
    SUBROUTINE SkipPropertiesList(funit,EndOfListStr_in)
      
      
      ! --------------------
      ! subroutine arguments
      ! --------------------
      INTEGER,                    INTENT(IN) :: funit
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: EndOfListStr_in
      ! ---------------
      ! Local variables
      ! ---------------
      LOGICAL                :: NOT_END_OF_LIST
      INTEGER                :: N
      CHARACTER(LEN=maxChar) :: SUBSTRS(maxChar), EndOFListStr
      
      ! ============================================================
      ! SkipPropertiesList starts here
      ! ============================================================
      
      ! Set default
      EndOfListStr = '###END_OF_LIST###'
      IF(PRESENT(EndOfListStr_in)) EndOfListStr = EndOfListStr_in

      NOT_END_OF_LIST = .TRUE.
      DO WHILE(NOT_END_OF_LIST)
        
        ! Read one line
        CALL SPLIT_ONE_LINE( funit, SUBSTRS, N, -1, 'read_diagnostic_options:2' )
        
        ! Check for end of list
        IF(TRIM(ADJUSTL(SUBSTRS(1))) .EQ. TRIM(ADJUSTL(EndOfListStr))) THEN
          NOT_END_OF_LIST = .FALSE.
        ENDIF
      
      ENDDO
    END SUBROUTINE SkipPropertiesList
    

    
    !###################################################################
    !#                              SPLAT                              #
    !###################################################################
    
    ! FUNCTION: DO_ALL_SPC
    ! 
    ! DESCRIPTION: Checks if the species diagnostic name "all" has been set
    
    LOGICAL FUNCTION DO_ALL_SPC( spcname )
      
      USE string_utils_module, ONLY : TRANUC
      
      ! --------------------
      ! Subroutine arguments
      ! --------------------
      CHARACTER(LEN=maxChar), INTENT(IN) :: spcname
      
      ! ---------------
      ! Local variables
      ! ---------------
      CHARACTER(LEN=maxChar)           :: upper_str
      
      ! ============================================================
      ! DO_ALL_SPC starts here
      ! ============================================================
      
      ! Init
      DO_ALL_SPC = .FALSE.
      
      ! Check first element is not "ALL"
      upper_str = spcname
      CALL TRANUC( upper_str )
      
      IF( (TRIM(ADJUSTL(upper_str)).EQ. 'ALL' ) ) DO_ALL_SPC = .TRUE.
      
      RETURN
      
    END FUNCTION DO_ALL_SPC
    
END MODULE input_module
