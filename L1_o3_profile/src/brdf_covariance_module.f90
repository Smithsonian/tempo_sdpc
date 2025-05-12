MODULE covariance_module

  ! Contains routines for setting/manipulating covariance matrices

  USE parameters_module
  USE error_module,      ONLY : ErrorType, RaiseFatalError, CheckError
  IMPLICIT NONE

  TYPE CovarOptType
    CHARACTER(LEN=maxChar)    :: ParType
    INTEGER                   :: nPar
    REAL(KIND=8), ALLOCATABLE :: Par(:)
  ENDTYPE CovarOptType

  ! Error uncertainties
  TYPE UncertType
    INTEGER(KIND=2)           :: StateType ! 0-Scalar, 1-Profile, 2-GDF, 3-EXP, 4-BOX 
    INTEGER(KIND=2)           :: CovarParType ! 0-Scalar, 1-Diagnonal, 2-ZCorrelated, 3-FullySpecified
    INTEGER(KIND=4)           :: nScalarPar
    REAL(KIND=8), ALLOCATABLE :: ScalarPar(:)
    REAL(KIND=8), ALLOCATABLE :: ProfilePar(:)
    REAL(KIND=8), ALLOCATABLE :: SubCovMatrix(:,:)
  ENDTYPE UncertType

  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'level1_module'
  PRIVATE :: ModuleName

  CONTAINS
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: GetOutputCovarianceDim
  ! 
  ! DESCRIPTION: Returns the dimension of the sub covariance matrix 
  !              for a given UncertType

  INTEGER FUNCTION GetOutputCovarDim(Uncertainty,lmx)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(UncertType), INTENT(IN) :: Uncertainty
    INTEGER,          INTENT(IN) :: lmx ! Vertical grid dimension

    ! ---------------
    ! local variables
    ! ---------------

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'GetOutputCovarDim'

    ! ============================================================
    ! GetOutputCovarDim starts here
    ! ============================================================

    GetOutputCovarDim = 1
    IF(Uncertainty%StateType .EQ. 1) THEN ! Full Profile
      GetOutputCovarDim =  lmx
    ELSEIF(Uncertainty%StateType .EQ. 2) THEN ! GDF
      GetOutputCovarDim = 3
    ELSEIF(Uncertainty%StateType .EQ. 3) THEN ! EXP
      GetOutputCovarDim = 2
    ENDIF

  END FUNCTION GetOutputCovarDim

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
  
  ! SUBROUTINE: SetSubStateCovariance
  ! 
  ! DESCRIPTION: ---
  
  SUBROUTINE SetSubStateCovariance(nY,nX,IsScaleFactor,StateType,InpFileCovar,&
                                   OverwriteClimatology,ClimUncertainty,      &
                                   lmx, zmid, Sa, Error, V, ProfPar, LayerCol,&
                                   ProfParDeriv                               )

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,               INTENT(IN)    :: nY,nX
    LOGICAL,               INTENT(IN)    :: IsScaleFactor
    CHARACTER(LEN=*),      INTENT(IN)    :: StateType ! 'Scalar','Column','Profile','ProfilePar' 
    TYPE(CovarOptType),    INTENT(IN)    :: InpFileCovar ! The new options
    LOGICAL,               INTENT(IN)    :: OverwriteClimatology
    TYPE(UncertType),      INTENT(IN)    :: ClimUncertainty
    INTEGER,               INTENT(IN)    :: lmx
    REAL(KIND=8),          INTENT(IN)    :: zmid(lmx)
    REAL(KIND=8),          INTENT(INOUT) :: Sa(nX,nX)
    TYPE(ErrorType),       INTENT(INOUT) :: Error
    REAL(KIND=8), OPTIONAL,INTENT(IN)    :: V(nX) ! For scale factor
    REAL(KIND=8), OPTIONAL,INTENT(IN)    :: ProfPar(MaxProfPar)
    REAL(KIND=8), OPTIONAL,INTENT(IN)    :: LayerCol(lmx)
    REAL(KIND=8), OPTIONAL,INTENT(IN)    :: ProfParDeriv(lmx,MaxProfPar)
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: l, nXc, P
    REAL(KIND=8) :: TotalCol, H(lmx,1), SaHT(lmx,1)
    REAL(KIND=8), ALLOCATABLE :: K(:,:), ScKT(:,:)
    REAL(KIND=8) :: Sprof(lmx,lmx)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SetSubStateCovariance'

    ! ============================================================
    ! SetSubStateCovariance starts here
    ! ============================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Check V is present
    IF(.NOT. PRESENT(V) .AND. IsScaleFactor) THEN
      STOP 'Must use optional argument V for scale factor case'
    ENDIF
    
    IF(StateType .EQ. 'Scalar' .AND. nX .NE. 1) THEN
      STOP 'Statetype set to scalar but X substate is a vector!'
    ENDIF

    IF(StateType .EQ. 'Column' .AND. .NOT. PRESENT(LayerCol)) THEN
      STOP 'When Column Type must set LayerCol'
    ENDIF


    IF(OverwriteClimatology) THEN  ! Compute parameterization from CovarOpt

      ! Scalar Type
      IF(StateType .EQ. 'Scalar') THEN
        Sa(1,1) = InpFileCovar%Par(1)**2

      ! Profile Type
      ELSE
        CALL ProfCovarParameterization(InpFileCovar,lmx,zmid,nX,Sa,Error)

      ENDIF

    ELSE
      
      ! Dimension of climatology output covariance matrix
      nXc = GetOutputCovarDim(ClimUncertainty,lmx)

      ! Compute Profile Jacobian if needed
      IF(StateType .EQ. 'Profile' .OR. StateType .EQ. 'Column') THEN 
        
        ! If using climatology copy the error covariance matrix
        IF(ClimUncertainty%StateType .EQ. 1) THEN
          Sprof = ClimUncertainty%SubCovMatrix
        ELSE

          ! Allocate  Jacobian
          ALLOCATE(K(lmx,nXc)) ; K(:,:) = 0.0D0

          ! Compute Jacobian
          DO L=1,lmx
            IF(LayerCol(L) .GT. TINY(0.0D0)) THEN
              DO P=1,nXc
                K(L,P) = ProfParDeriv(L,P)*ProfPar(P)/LayerCol(L)
              ENDDO
            ENDIF
          ENDDO

          ! Compute the covariance matrix K*S_a*K^T
          ALLOCATE(ScKT(nXc,lmx)) ; ScKT(:,:) = 0.0d0
          CALL DGEMM('N','T',nXc,lmx,nXc,1.0D0,&
                     ClimUncertainty%SubCovMatrix,&
                     nXc,K,lmx,0.0D0,ScKT,nXc)
          CALL DGEMM('N','N',lmx,lmx,nXc,1.0d0, &
                     K,lmx,ScKT,nXc,0.0D0,Sprof,lmx)

        ENDIF

        ! Set the covariance matrix
        IF(StateType .EQ. 'Profile') THEN

          ! Can just copy the profile coviance
          Sa = Sprof

        ELSEIF(StateType .EQ. 'Column') THEN

          ! Compute column operator
          TotalCol = SUM(LayerCol)
          H(:,1) = LayerCol/TotalCol
          
          ! Compute the covariance of the column h*Sa*h^T
          SaHT(:,:) = 0.0d0
          CALL DGEMM('N','N',lmx,1,lmx,1.0D0,      &
                     Sprof,lmx,H,lmx,0.0D0,SaHT,lmx)
          CALL DGEMM('T','N',1,1,lmx,1.0D0,   &
                     H,lmx,SaHT,lmx,0.0d0,Sa,1)

        ENDIF

      ELSEIF(StateType .EQ. 'ProfilePar') THEN
        
        ! Check if statetypes match (based on dimension)
        IF(nX .EQ. nXc) THEN
          Sa = ClimUncertainty%SubCovMatrix
        ELSE
          print*,'Only Parameterized covariance of the same ' &
               // 'profile type may be used in ProfilePar case'
          STOP 1
        ENDIF
      
      ELSEIF(StateType .EQ. 'Scalar') THEN
        Sa(:,:) = ClimUncertainty%SubCovMatrix(:,:) ! 1x1
      ENDIF

    ENDIF
    
    ! Rescale to scale factors if set
    IF(IsScaleFactor) THEN
      CALL ScaleFactorCovarianceMatrix(Nx,V,Sa,Error)
    ENDIF

  END SUBROUTINE SetSubStateCovariance

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: ProfCovarParameterization
  
  ! DESCRIPTION: Routine for computing the various covariance 
  !              parameterization options
  
  SUBROUTINE ProfCovarParameterization(InpFileCovar,lmx,zmid,nX,Sa,Error,Sa0)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(CovarOptType),     INTENT(IN)    :: InpFileCovar
    INTEGER,                INTENT(IN)    :: lmx
    REAL(KIND=8),           INTENT(IN)    :: zmid(lmx)
    INTEGER,                INTENT(IN)    :: nX
    REAL(KIND=8),           INTENT(INOUT) :: Sa(nX,nX)
    TYPE(ErrorType),        INTENT(INOUT) :: Error
    REAL(KIND=8), OPTIONAL, INTENT(IN)    :: Sa0(nX,nX)
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: I,J
    REAL(KIND=8) :: dz

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ProfCovarParameterization'

    ! ============================================================
    ! CovarianceParameterization starts here
    ! ============================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Zero covariance matrix
    Sa(:,:) = 0.0d0

    IF(TRIM(ADJUSTL(InpFileCovar%ParType)) .EQ. 'SCALAR') THEN

      DO I=1,nX
        Sa(I,I) = InpFileCovar%Par(1)**2
      ENDDO

    ELSEIF(TRIM(ADJUSTL(InpFileCovar%ParType)) .EQ. 'DIAGONAL') THEN

      DO I=1,nX
        Sa(I,I) = InpFileCovar%Par(I)**2
      ENDDO

    ELSEIF(TRIM(ADJUSTL(InpFileCovar%ParType)) .EQ. 'ZCORRELATED') THEN

      IF(nX .NE. lmx) THEN
        print*,'FATAL: Cannot use ZCORRELATED Covariance Param. when nX != lmx'
        STOP 1
      ENDIF

      DO I=1,lmx
      DO J=1,lmx
        dz = ABS(zmid(I)-zmid(J))
        Sa(I,J) = EXP(-1.0d0*dz/InpFileCovar%Par(2))*InpFileCovar%Par(1)**2
      ENDDO
      ENDDO

    ELSEIF(TRIM(ADJUSTL(InpFileCovar%ParType)) .EQ. 'CLIMSCALING') THEN

      IF(PRESENT(Sa0)) THEN
        Sa(:,:) = Sa0(:,:)*InpFileCovar%Par(1)
      ELSE
        STOP 'Sa0 Must be present to apply covariance scaling'
      ENDIF

    ELSE

      print*,InpFileCovar%ParType,': Unrecognized covariance parameterization type'

    ENDIF

  END SUBROUTINE ProfCovarParameterization

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: ComputeZCorrCovar
  ! 
  ! DESCRIPTION: Routine for computing a covariance matrix with a 
  !              correlation length scale 

  SUBROUTINE ComputeZCorrCovar(lmx,Sdiag,z,zcorr,Sa)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,             INTENT(IN)    :: lmx
    REAL(KIND=8),        INTENT(IN)    :: Sdiag(lmx)
    REAL(KIND=8),        INTENT(IN)    :: z(lmx)
    REAL(KIND=8),        INTENT(IN)    :: zcorr
    REAL(KIND=8),        INTENT(OUT)   :: Sa(lmx,lmx)
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: l,k
    REAL(KIND=8) :: dz

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ComputeZCorrCovar'

    ! ============================================================
    ! CovarParameterization starts here
    ! ============================================================

    DO k=1,lmx
    DO l=1,lmx
      dz = z(k)-z(l)
      Sa(l,k) = Sdiag(l)*Sdiag(k)*EXP(-1.0d0*ABS(dz)/zcorr)
    ENDDO
    ENDDO

  END SUBROUTINE ComputeZCorrCovar

  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: CovarParameterization
  ! 
  ! DESCRIPTION: Routine for computing the various covariance 
  !              parameterization options for non-profile variables
  
  SUBROUTINE CovarParameterization(CovarOpt,nX,Sa,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(CovarOptType),  INTENT(IN)    :: CovarOpt
    INTEGER,             INTENT(IN)    :: nX
    REAL(KIND=8),        INTENT(INOUT) :: Sa(nX,nX)
    TYPE(ErrorType),     INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: I

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'CovarParameterization'

    ! ============================================================
    ! CovarParameterization starts here
    ! ============================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Zero covariance matrix
    Sa(:,:) = 0.0d0
    
    IF(TRIM(ADJUSTL(CovarOpt%ParType)) .EQ. 'DIAGONAL') THEN
      
      DO I=1,nX
        Sa(I,I) = CovarOpt%Par(I)**2
      ENDDO
    
    ELSEIF(TRIM(ADJUSTL(CovarOpt%ParType)) .EQ. 'SCALAR') THEN
      
      DO I=1,nX
        Sa(I,I) = CovarOpt%Par(1)**2
      ENDDO
      
    ELSE
      print*,CovarOpt%ParType,': Unrecognized covariance parameterization type'
      
    ENDIF

  END SUBROUTINE CovarParameterization
  
  !###################################################################
  !#                              SPLAT                              #
  !###################################################################
    
  ! SUBROUTINE: ScaleFactorCovarianceMatrix
  ! 
  ! DESCRIPTION: Convert a prior covariance matrix from absolute values
  !              for a scale factor type state subvector

  SUBROUTINE ScaleFactorCovarianceMatrix(lmx,v,Sa,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,         INTENT(IN)    :: lmx
    REAL(KIND=8),    INTENT(IN)    :: v(lmx)
    REAL(KIND=8),    INTENT(INOUT) :: Sa(lmx,lmx)
    TYPE(ErrorType), INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i,j
    REAL(KIND=8) :: Sa_in(lmx,lmx), zero

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ScaleFactorCovarianceMatrix'

    ! ============================================================
    ! ScaleFactorCovarianceMatrix starts here
    ! ============================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Save copy of input
    Sa_in = Sa

    ! Zero Sa
    Sa(:,:) = 0.0d0

    ! "zero"
    zero = TINY(0.0d0)
    DO j=1,lmx
    DO i=1,lmx

      IF(v(i) .GT. zero .AND. v(j) .GT. zero) THEN
        Sa(i,j) = Sa_in(i,j)/v(i)/v(j)
      ENDIF

    ENDDO
    ENDDO

  END SUBROUTINE ScaleFactorCovarianceMatrix

END MODULE covariance_module
