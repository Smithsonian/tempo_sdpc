MODULE pca_bin_module

  USE parameters_module
  USE error_module,        ONLY : ErrorType, CheckError, RaiseFatalError
  USE string_utils_module, ONLY : SPLIT_ONE_LINE, READ_ONE_LINE,&
                                  STRREPL, SkipFileLines
  
  IMPLICIT NONE

  TYPE PCABinInpType

    ! Bin Mode type
    INTEGER                   :: BinModeIndex
    INTEGER                    :: nbin
    INTEGER                    :: nbin0
    REAL(KIND=8),  ALLOCATABLE :: INITIAL_BINLIMS(:)
    INTEGER                    :: THRESH_COUNT

    ! New Strategy 
    LOGICAL :: alb_pcainclude
    LOGICAL :: do_4M_correction


  ENDTYPE PCABinInpType
    
  ! Interface to PCA Binning Function
  INTERFACE
    SUBROUTINE PcaBinInterface(PCABinInp,MaxWav,MaxLayers,MaxBins,NLAY,NDAT,&
                               gasdat,taudp,omega,nbin,ncnt_new,index_new,bin_new,binlims)
      IMPORT
      TYPE(PCABinInpType), INTENT(IN)    :: PCABinInp
      INTEGER,             INTENT(IN)    :: MaxWav
      INTEGER,             INTENT(IN)    :: MaxLayers
      INTEGER,             INTENT(IN)    :: MaxBins
      INTEGER,             INTENT(IN)    :: NLAY, NDAT
      REAL(KIND=4),        INTENT(IN)    :: gasdat(Maxlayers,Maxwav)
      REAL(KIND=4),        INTENT(IN)    :: taudp(Maxlayers,MaxWav)
      REAL(KIND=4),        INTENT(IN)    :: omega(Maxlayers,MaxWav)
      INTEGER,             INTENT(OUT)   :: NBIN
      INTEGER,             INTENT(OUT)   :: ncnt_new(0:MaxBins)
      INTEGER,             INTENT(OUT)   :: index_new(MaxWav)
      INTEGER,             INTENT(OUT)   :: bin_new(NDAT)
      REAL(KIND=8),        INTENT(OUT)   :: binlims(0:MaxBins)
    END SUBROUTINE PcaBinInterface
  END INTERFACE

  ! PCA Binning
  TYPE PCABinType
    TYPE(PCABinInpType)       :: Settings
        
    ! These are the names of the variables in the MethaneSAT PCA Driver
    INTEGER,      ALLOCATABLE :: BIN_NEW(:)
    INTEGER,      ALLOCATABLE :: NCT_NEW(:)
    INTEGER,      ALLOCATABLE :: INDEX_NEW(:)
    REAL(KIND=8), ALLOCATABLE :: BINLIMS(:)
    
    ! Pointer to Bin creation function
    PROCEDURE(PcaBinInterface), NOPASS, POINTER :: BinningSubroutine

    ! Functions
    CONTAINS
      PROCEDURE :: InitializePCABins => InitializePCABins
      PROCEDURE :: CreatePCABins     => CreatePCABins

  ENDTYPE PCABinType

  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'pca_bin_module'
  PRIVATE :: ModuleName

  CONTAINS

  SUBROUTINE InitializePCABins(self,ControlFile,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(PCABinType),   INTENT(INOUT) :: self
    CHARACTER(LEN=*),    INTENT(IN)    :: ControlFile
    TYPE(ErrorType),     INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: LINE
    CHARACTER(LEN=1)       :: TAB   = ACHAR(9)
    CHARACTER(LEN=1)       :: SPACE = ' '
    LOGICAL                :: EOF, PCA_HAS_BEEN_READ, OPTSTRAT_HAS_BEEN_READ

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitializePCABins'

    ! =====================================================================
    ! InitializePCABins begins here
    ! =====================================================================

    IF(CheckError(Error)) RETURN

    ! Open Control File
    OPEN(UNIT=pcaunit,FILE=TRIM(ADJUSTL(ControlFile)),ACTION='READ', STATUS='OLD')

    ! Check control file has been read
    PCA_HAS_BEEN_READ = .FALSE.
    OPTSTRAT_HAS_BEEN_READ = .FALSE.
    DO

      ! Read a line from the file, exit if EOF
      LINE = READ_ONE_LINE( pcaunit, EOF )
      IF ( EOF ) EXIT

      ! Replace tab characters in LINE (if any) w/ spaces
      CALL STRREPL( LINE, TAB, SPACE )

      IF(INDEX( LINE, '%%% PCA BIN OPTIONS  %%%'    ) > 0 ) THEN
        CALL ReadPCABinOptions(self%Settings,Error)
        PCA_HAS_BEEN_READ = .TRUE.
      ELSEIF(INDEX( LINE, '%%% CALCULATION OPTIONS %%%') > 0 ) THEN
        CALL ReadPCAOptPropStrat(self%Settings,Error)
        OPTSTRAT_HAS_BEEN_READ = .TRUE.
      ELSEIF(INDEX( LINE, 'END OF FILE'             ) > 0 ) THEN 
        EXIT
      ENDIF

    ENDDO

    ! Close file
    CLOSE(pcaunit)

    ! Check file has been read
    IF(.NOT. PCA_HAS_BEEN_READ) THEN
      CALL RaiseFatalError( Error, ErrorCode_Input, ModuleName, SubroutineName,&
                            Message_in='Could Not read PCA control file (Bin)' )
    ENDIF
    
    IF(.NOT. OPTSTRAT_HAS_BEEN_READ) THEN
      CALL RaiseFatalError( Error, ErrorCode_Input, ModuleName, SubroutineName,   &
                            Message_in='Could Not read PCA control file(Optstrat)')
    ENDIF

    ! Point to function for selected optical property strategy
    self%BinningSubroutine =>  null()
    IF(self%Settings%BinModeIndex .EQ. 1) THEN
      STOP 'PCA Bin Type 1 has not been implemented yet'
    ELSEIF(self%Settings%BinModeIndex .EQ. 2) THEN
      STOP 'PCA Bin Type 2 has not been implemented yet'
    ELSEIF(self%Settings%BinModeIndex .EQ. 3) THEN
      self%BinningSubroutine => BinningSubroutine_Type3
    ELSE

    ENDIF

  END SUBROUTINE InitializePCABins

  SUBROUTINE ReadPCABinOptions(PCABinInp,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(PCABinInpType),   INTENT(INOUT) :: PCABinInp
    TYPE(ErrorType),       INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: SUBSTRS(maxChar)
    INTEGER                :: I,N

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ReadPCABinOptions'

    ! =====================================================================
    ! ReadPCABinOptions begins here
    ! =====================================================================

    IF(CheckError(Error)) RETURN
    
    ! Get the bin mode index
    CALL SPLIT_ONE_LINE( pcaunit, SUBSTRS, N, -1, 'pca_bin_options:1' )
    READ( SUBSTRS(1), * ) PCABinInp%BinModeIndex

    ! Skip line (Header for Type 1)
    CALL SkipFileLines(pcaunit,1)

    ! Skip header for Type 2
    CALL SkipFileLines(pcaunit,1)

    ! Skip header for Type 3
    CALL SkipFileLines(pcaunit,1)

    ! Read Options if type 3
    IF(PCABinInp%BinModeIndex .EQ. 3) THEN

      ! Bin Optical Dephs
      CALL SPLIT_ONE_LINE( pcaunit, SUBSTRS, N, -1, 'pca_bin_options:3' )
      PCABinInp%NBIN = N
      ALLOCATE(PCABinInp%INITIAL_BINLIMS(N))
      DO I=1,PCABinInp%NBIN
        READ(SUBSTRS(I),*) PCABinInp%INITIAL_BINLIMS(I)
      ENDDO 

      ! Count threshold
      CALL SPLIT_ONE_LINE( pcaunit, SUBSTRS, N, -1, 'pca_bin_options:3' )
      READ(SUBSTRS(1),*) PCABinInp%THRESH_COUNT

    ELSE
      CALL SkipFileLines(pcaunit,2)
    ENDIF

  END SUBROUTINE ReadPCABinOptions

  SUBROUTINE ReadPCAOptPropStrat(PCABinInp,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(PCABinInpType),   INTENT(INOUT) :: PCABinInp
    TYPE(ErrorType),       INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: SUBSTRS(maxChar)
    INTEGER                :: I,N,M

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ReadPCAOptPropStrat'

    ! =====================================================================
    ! ReadPCAOptPropStrat begins here
    ! =====================================================================

    IF(CheckError(Error)) RETURN

    ! Use Albedo in PCA
    CALL SPLIT_ONE_LINE( pcaunit, SUBSTRS, N, -1, 'pca_optprop_options:1' )
    READ( SUBSTRS(1), * ) PCABinInp%alb_pcainclude

    ! 4m correction
    CALL SPLIT_ONE_LINE( pcaunit, SUBSTRS, N, -1, 'pca_optprop_options:1' )
    READ( SUBSTRS(1), * ) PCABinInp%do_4M_correction


  END SUBROUTINE ReadPCAOptPropStrat

  SUBROUTINE CreatePCABins(self,MaxWav,MaxLayers,MaxBins,NLAY,NDAT,&
                          gasdat,taudp,omega,nbin,ncnt_new,index_new,bin_new,binlims)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    CLASS(PCABinType),   INTENT(INOUT)   :: self
    INTEGER,             INTENT(IN)    :: MaxWav
    INTEGER,             INTENT(IN)    :: MaxLayers
    INTEGER,             INTENT(IN)    :: MaxBins
    INTEGER,             INTENT(IN)    :: NLAY, NDAT
    REAL(KIND=4),        INTENT(IN)    :: gasdat(Maxlayers,Maxwav)
    REAL(KIND=4),        INTENT(IN)    :: taudp(Maxlayers,MaxWav)
    REAL(KIND=4),        INTENT(IN)    :: omega(Maxlayers,MaxWav)
    INTEGER,             INTENT(OUT)   :: nbin
    INTEGER,             INTENT(OUT)   :: ncnt_new(0:MaxBins)
    INTEGER,             INTENT(OUT)   :: index_new(MaxWav)
    INTEGER,             INTENT(OUT)   :: bin_new(NDAT)
    REAL(KIND=8),        INTENT(OUT)   :: binlims(0:MaxBins)
    ! =====================================================================
    ! CreatePCABins begins here
    ! =====================================================================

    ! This is just a wrapper for the various binning schemes
    CALL self%BinningSubroutine(self%Settings,MaxWav,MaxLayers,MaxBins,NLAY,NDAT,&
                            gasdat,taudp,omega,nbin,ncnt_new,index_new,bin_new,binlims)

  END SUBROUTINE CreatePCABins

  SUBROUTINE BinningSubroutine_Type3(PCABinInp,MaxWav,MaxLayers,MaxBins,NLAY,NDAT,&
                                     gasdat,taudp,omega,nbin,ncnt_new,index_new,bin_new,binlims)

    ! --------------------
    ! subroutine arguments
    ! --------------------

    TYPE(PCABinInpType), INTENT(IN)    :: PCABinInp
    INTEGER,             INTENT(IN)    :: MaxWav
    INTEGER,             INTENT(IN)    :: MaxLayers
    INTEGER,             INTENT(IN)    :: MaxBins
    INTEGER,             INTENT(IN)    :: NLAY, NDAT
    REAL(KIND=4),        INTENT(IN)    :: gasdat(Maxlayers,Maxwav)
    REAL(KIND=4),        INTENT(IN)    :: taudp(Maxlayers,MaxWav)
    REAL(KIND=4),        INTENT(IN)    :: omega(Maxlayers,MaxWav)
    INTEGER,             INTENT(OUT)   :: NBIN
    INTEGER,             INTENT(OUT)   :: ncnt_new(0:MaxBins)
    INTEGER,             INTENT(OUT)   :: index_new(MaxWav)
    INTEGER,             INTENT(OUT)   :: bin_new(NDAT)
    REAL(KIND=8),        INTENT(OUT)   :: binlims(0:MaxBins)

    ! ---------------
    ! local variables
    ! ---------------

    !  precision
    integer, parameter :: sp = SELECTED_REAL_KIND(6)
    integer, parameter :: dp = SELECTED_REAL_KIND(15)

    ! JBak alteration, 8/2/18. Add nbin0, which_jbak, STEP_INDEX
    logical :: iterating
    INTEGER :: IBIN, W, COUNT, N, K, NBIN_OLD, NBIN_NEW, THRESH_COUNT, STEP_INDEX

    INTEGER, ALLOCATABLE :: NCNT_OLD(:), INDEX_OLD(:)
    INTEGER, ALLOCATABLE :: BIN_OLD(:), BINDEX_OLD(:), BINDEX_NEW(:)
    ! INTEGER :: NCNT_OLD(0:Maxbins), INDEX_OLD(MaxWav)
    ! INTEGER :: BIN_OLD(NDAT), BINDEX_OLD(NDAT), BINDEX_NEW(NDAT)
    REAL(kind=dp) :: TAUTOT, DTAU
    REAL(kind=dp), ALLOCATABLE ::  REAL_INITIAL_BINLIMS(:), INITIAL_BINLIMS(:)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'CreatePCABins_Type3'

    ! =====================================================================
    ! CreatePCABins_Type3 starts here
    ! =====================================================================
    
    ! Initialize Input values from settings
    NBIN = PCABinInp%NBIN
    THRESH_COUNT = PCABinInp%THRESH_COUNT
    
    ! Set initial bin limits
    ALLOCATE(INITIAL_BINLIMS(NBIN))
    DO IBIN=1,NBIN
      INITIAL_BINLIMS(IBIN) = LOG(PCABinInp%INITIAL_BINLIMS(IBIN))
    ENDDO

    ! Zero Output
    NCNT_NEW = 0 ; BINLIMS = 0.0_dp ; INDEX_NEW = 0 ; BIN_NEW = 0

    ! Allocate/Zero "OLD" values
    ALLOCATE(NCNT_OLD(0:MaxBins)) ; NBIN_OLD = NBIN ; NCNT_OLD = 0
    ALLOCATE(BIN_OLD(NDAT))     ; BIN_OLD = 0
    ALLOCATE(INDEX_OLD(MaxWav)) ; INDEX_OLD = 0
    ALLOCATE(BINDEX_OLD(NDAT))  ; BINDEX_OLD = 0
    ALLOCATE(BINDEX_NEW(NDAT))  ; BINDEX_NEW = 0

    print*,'NBIN:',NBIN
    print*,'THRESH_COUNT:',THRESH_COUNT
    print*,INITIAL_BINLIMS
    ! start wavelength loop
    DO W = 1, NDAT

      ! Compute log(total profile OD)
      TAUTOT = REAL(SUM(GASDAT(1:nlay,W)),dp); DTAU = LOG(TAUTOT)

      ! Assign bin
      DO N=1,NBIN

        ! Last bin index
        IF(N .EQ. NBIN) THEN
          IF(DTAU .LT. INITIAL_BINLIMS(N)) THEN
            IBIN = N-1 ; BINLIMS(IBIN) = INITIAL_BINLIMS(N)
          ELSE
            IBIN = N ; BINLIMS(IBIN) = 5000.0 ! Change based on window limit
          ENDIF
          EXIT

        ! Somewhere in the middle/below limit
        ELSE
          IF(DTAU .LT. INITIAL_BINLIMS(N)) THEN
            IBIN = N-1 ; BINLIMS(IBIN) = INITIAL_BINLIMS(N)
            EXIT
          ENDIF
        ENDIF

      ENDDO

      ! Update helper indices
      NCNT_OLD(IBIN) = NCNT_OLD(IBIN)+1
      BIN_OLD(W)     = IBIN
      BINDEX_OLD(W)  = NCNT_OLD(IBIN)

    ENDDO

    !  get the indices to change between bin space and wavenumber space
    DO W = 1, NDAT
      IF (BIN_OLD(W) .EQ. 0) THEN
        COUNT = BINDEX_OLD(W)
      ELSE
        COUNT = SUM(NCNT_OLD(0:BIN_OLD(W)-1)) + BINDEX_OLD(W)
      ENDIF
      INDEX_OLD(COUNT) = W
    ENDDO

    
    DO !  Now reduce the number of bins by 1 (from the bottom), fills up next bin in sequence
      
      K = NBIN_OLD ; ITERATING = .FALSE.
      IF( NCNT_OLD(0) .lt. THRESH_COUNT ) THEN
        ITERATING = .TRUE.
        NBIN_NEW = NBIN_OLD - 1
        NCNT_NEW(1:K-2) = NCNT_OLD(2:K-1)
        NCNT_NEW(0) = NCNT_OLD(0)
      ENDIF

      
      ! Check if it has reduced
      IF(.NOT. ITERATING) THEN
        
        ! copy old to new
        NBIN_NEW = NBIN_OLD
        NCNT_NEW = NCNT_OLD
        DO W = 1, NDAT
          BINDEX_NEW(W) = BINDEX_OLD(W)
          BIN_NEW(W)    = BIN_OLD(W)
          IF (BIN_NEW(W) .EQ. 0) THEN
            COUNT = BINDEX_NEW(W)
          ELSE
            COUNT = SUM(NCNT_NEW(0:BIN_NEW(W)-1)) + BINDEX_NEW(W)
          ENDIF
          INDEX_NEW(COUNT) = W
        enddo
        NBIN = NBIN_NEW

        ! Exit Loop 
        EXIT

      ENDIF

      !  Fill up BIN_NEW
      DO W = 1, NDAT
        IF ( BIN_OLD(W).eq. 1 ) then
          IBIN = 0
          BIN_NEW(W) = IBIN
          NCNT_NEW(IBIN) = NCNT_NEW(IBIN) + 1
          BINDEX_NEW(W)  = NCNT_NEW(IBIN)
        ELSE
          BINDEX_NEW(W) = BINDEX_OLD(W)
          BIN_NEW(W)    = BIN_OLD(W) - 1
        ENDIF
      ENDDO

      !  Re-assign output
      NBIN = NBIN_NEW
      DO IBIN = 0, NBIN
        BINLIMS(IBIN) = BINLIMS(IBIN+1) 
      ENDDO
      BINLIMS(NBIN) = 0.0_dp 

      ! Re-set and return to iteration
      NBIN_OLD = NBIN_NEW
      NCNT_OLD = 0
      NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
      BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
      BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

    ENDDO

    !  Re-set
    NBIN_OLD = NBIN_NEW
    NCNT_OLD = 0
    NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
    BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
    BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

    !  Now reduce the number of bins by 1 (from the top), fills up penumltimate bin
    DO 

      K = NBIN_OLD-1 ; ITERATING = .false.
      IF ( NCNT_OLD(K) .lt. THRESH_COUNT ) THEN
        ITERATING = .true.
        NBIN_NEW = NBIN_OLD - 1
        NCNT_NEW(0:K-2) = NCNT_OLD(0:K-2)
        NCNT_NEW(K-1) = NCNT_OLD(K-1) !+ NCNT_OLD(K)
      ENDIF

      !  Not iterating : copy Old to new and exit
      IF( .not. ITERATING ) then
        NBIN_NEW = NBIN_OLD
        NCNT_NEW = NCNT_OLD
        DO W = 1, NDAT
          BINDEX_NEW(W) = BINDEX_OLD(W)
          BIN_NEW(W)    = BIN_OLD(W)
          IF (BIN_NEW(W) .EQ. 0) THEN
            COUNT = BINDEX_NEW(W)
          ELSE
            COUNT = SUM(NCNT_NEW(0:BIN_NEW(W)-1)) + BINDEX_NEW(W)
          ENDIF
          INDEX_NEW(COUNT) = W
        enddo
        NBIN = NBIN_NEW

        ! Exit Loop
        EXIT
      ENDIF

      !  Fill up BIN_NEW
      DO W = 1, NDAT
        IF ( BIN_OLD(W).eq. NBIN_OLD - 1) then
          IBIN = NBIN_NEW - 1
          BIN_NEW(W) = IBIN
          NCNT_NEW(IBIN) = NCNT_NEW(IBIN) + 1
          BINDEX_NEW(W)  = NCNT_NEW(IBIN)
        ELSE
          BINDEX_NEW(W) = BINDEX_OLD(W)
          BIN_NEW(W)    = BIN_OLD(W)
        ENDIF
      ENDDO

      !  Re-assign output
      NBIN = NBIN_NEW
      BINLIMS(NBIN_NEW) = BINLIMS(NBIN_OLD) 
      BINLIMS(NBIN_OLD) = 0.0_dp 

      !  Re-set and return to iteration
      NBIN_OLD = NBIN_NEW
      NCNT_OLD = 0
      NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
      BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
      BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

    ENDDO

    ! Cleanup
    DEALLOCATE(INITIAL_BINLIMS,NCNT_OLD,BIN_OLD,INDEX_OLD,BINDEX_OLD,BINDEX_NEW)

  END SUBROUTINE BinningSubroutine_Type3


  ! ###############################################################
  ! OLD STUFF
  ! ###############################################################



  ! SUBROUTINE CreatePCABins(nwav,wav,PCABin,Error)

  !   ! --------------------
  !   ! subroutine arguments
  !   ! --------------------
  !   INTEGER,             INTENT(IN)    :: nwav
  !   REAL(KIND=8),        INTENT(IN)    :: wav(nwav)
  !   TYPE(PCABinType),    INTENT(INOUT) :: PCABin
  !   TYPE(ErrorType),     INTENT(INOUT) :: Error

  !   ! ---------------
  !   ! local variables
  !   ! ---------------

  !   ! For error checking
  !   CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'CreatePCABins'

  !   ! =====================================================================
  !   ! CreatePCABins begins here
  !   ! =====================================================================
    
  !   IF(PCABin%Settings%BinModeIndex .EQ. 1) THEN
  !     CALL CreatePCABins_Type1(nwav,wav,PCABin,Error)
  !   ELSEIF(PCABin%Settings%BinModeIndex .EQ. 2) THEN
  !     ! CALL CreatePCABins_Type1(nwav,nlvl,wav,tautot,omega,taugas,PCABin,Error)
  !   ELSE
  !     print*,PCABin%Settings%BinModeIndex
  !     CALL RaiseFatalError( Error, ErrorCode_OptProp , ModuleName, SubroutineName,&
  !                           Message_in='Unknwon PCA BinModeIndex'                 )
  !   ENDIF


  ! END SUBROUTINE CreatePCABins

  ! SUBROUTINE CreatePCABins_Type1(nwav,wav,PCABin,Error)

  !   ! --------------------
  !   ! subroutine arguments
  !   ! --------------------
  !   INTEGER,             INTENT(IN)    :: nwav
  !   REAL(KIND=8),        INTENT(IN)    :: wav(nwav)
  !   TYPE(PCABinType),    INTENT(INOUT) :: PCABin
  !   TYPE(ErrorType),     INTENT(INOUT) :: Error

  !   ! ---------------
  !   ! local variables
  !   ! ---------------
  !   INTEGER      :: k, it, w, cc, IBIN
  !   INTEGER      :: STEP_INDEX, COUNT1
  !   REAL(KIND=8) :: THIS_STEP_BIN, min_tau, dtau
  !   LOGICAL      :: Appending_SubBins

  !   INTEGER, ALLOCATABLE :: BINDEX(:), TAU(:)

  !   ! For error checking
  !   CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'CreatePCABins_Type1'

  !   ! =====================================================================
  !   ! CreatePCABins_Type1 begins here
  !   ! =====================================================================

  !   IF(CheckError(Error)) RETURN

  !   ! Initialize bin number
  !   PCABin%NBIN = 0

  !   ! Deallocate helper arrays
  !   IF(ALLOCATED(PCABin%STEP_BINS)) DEALLOCATE(PCABin%STEP_BINS)
  !   IF(ALLOCATED(PCABin%STEP_NEOF)) DEALLOCATE(PCABin%STEP_NEOF)

  !   ! Initialize minimum wavelength
  !   min_tau = PCABin%Settings%min_tau

  !   ! Append lower limit
  !   CALL AppendStepBin(PCABin,Error,min_tau)

  !   k = 1
  !   DO it=1,PCABin%Settings%niter

  !     ! Initialize flag for appending sub bins
  !     Appending_SubBins = .TRUE.

  !     w = 0
  !     DO WHILE(Appending_SubBins)
  !       IF(k .EQ. 1) THEN
  !         THIS_STEP_BIN = PCABin%Settings%iter_maxtau(it)
  !       ELSE
  !         THIS_STEP_BIN = REAL(w,KIND=8)*PCABin%Settings%iter_dtau(it) &
  !                       + PCABin%Settings%min_tau
  !       ENDIF

  !       ! Get numberof points within range
  !       cc = COUNT(mask = (wav(1:nwav) .GE. PCABin%STEP_BINS(K-1) .and. &
  !                          wav(1:nwav) .LT. THIS_STEP_BIN))
        
  !       IF(THIS_STEP_BIN .GE. PCABin%Settings%iter_maxtau(it)) THEN
  !         Appending_SubBins = .FALSE.
  !       ELSE

  !         ! Append the sub bin
  !         IF(CC .GE. PCABin%Settings%THRESH_COUNT) THEN
  !           CALL AppendStepBin(PCABin,Error,THIS_STEP_BIN,&
  !                              PCABin%Settings%iter_teof(it))
  !           PCABin%NBIN = PCABin%NBIN + 1
  !           k = k + 1
  !         ENDIF

  !       ENDIF
  !       w = w + 1
  !     ENDDO

  !     ! Update starting wavelength
  !     min_tau = PCABin%STEP_BINS(K-1)
  !     IF (cc > PCABin%Settings%THRESH_COUNT) THEN
  !       CALL AppendStepBin(PCABin,Error,THIS_STEP_BIN,&
  !                          PCABin%Settings%iter_teof(it))
  !       PCABin%NBIN = PCABin%NBIN + 1
  !       k = k + 1
  !     ENDIF

  !   ENDDO

  !   ! Reallocate Indexing arrays
  !   IF(ALLOCATED(PCABin%BIN)    ) DEALLOCATE(PCABin%BIN)
  !   IF(ALLOCATED(PCABin%INDEX)  ) DEALLOCATE(PCABin%INDEX)
  !   IF(ALLOCATED(PCABin%NCNT)   ) DEALLOCATE(PCABin%NCNT)
  !   IF(ALLOCATED(PCABin%NEOF)   ) DEALLOCATE(PCABin%NEOF)
  !   IF(ALLOCATED(PCABin%BINLIMS)) DEALLOCATE(PCABin%BINLIMS)
  !   ALLOCATE(PCABin%BIN(nwav))              ; PCABin%BIN(:) = 0
  !   ALLOCATE(PCABin%INDEX(nwav))            ; PCABin%INDEX(:) = 0
  !   ALLOCATE(PCABin%NCNT(0:PCABin%NBIN))    ; PCABin%NCNT(:) = 0
  !   ALLOCATE(PCABin%NEOF(0:PCABin%NBIN))    ; PCABin%NEOF(:) = 0
  !   ALLOCATE(PCABin%BINLIMS(0:PCABin%NBIN)) ; PCABin%BINLIMS(:) = 0.0
  !   ALLOCATE(BINDEX(nwav)) ; BINDEX(:) = 0

  !   ! Binning process
  !   IBIN = 0
  !   DO W = 1, nwav
  !     dtau = wav(w)
  !     STEP_INDEX = MINVAL(MAXLOC( PCABin%STEP_BINS(1:PCABin%NBIN),&
  !                           MASK=(PCABin%STEP_BINS(1:PCABin%NBIN) < dtau)))
  !     IF (dtau .LE. PCABin%STEP_BINS(1)) THEN
  !       IBIN = 0 ; STEP_INDEX = 1 !@bottom bin
  !     ELSEIF (dtau > PCABin%STEP_BINS(PCABin%NBIN)) THEN 
  !       IBIN = PCABin%NBIN
  !       STEP_INDEX = PCABin%NBIN
  !       PCABin%STEP_BINS(STEP_INDEX) = dtau !@top bin
  !     ELSEIF (dtau .GT. PCABin%STEP_BINS(1)      .AND. &
  !             dtau .LE. PCABin%STEP_BINS(PCABin%NBIN)  ) THEN
  !       IBIN = STEP_INDEX
  !       STEP_INDEX = STEP_INDEX + 1
  !     ENDIF
    
  !     PCABin%NCNT(IBIN) = PCABin%NCNT(IBIN)+1
  !     PCABin%NEOF(IBIN) = PCABin%STEP_NEOF(STEP_INDEX)
  !     PCABin%BINLIMS(IBIN)  = PCABin%STEP_BINS(STEP_INDEX)
  !     PCABin%BIN(W)     = IBIN
  !     BINDEX(W)  = PCABin%NCNT(IBIN)

  !     print*,dtau,  IBIN, PCABin%BINLIMS(IBIN), &
  !            PCABin%NCNT(IBIN)
  !   ENDDO

  !   !  get the indices to change between bin space and wavenumber space
  !   DO W = 1, nwav
  !     IF (PCABin%BIN(W) .EQ. 0) THEN
  !       COUNT1 = BINDEX(W)
  !     ELSE
  !       COUNT1 = SUM(PCABin%NCNT(0:PCABin%BIN(W)-1)) + BINDEX(W)
  !     ENDIF
  !     PCABin%INDEX(COUNT1) = W
  !   ENDDO
    
  !   ! Finish with deallocation
  !   deallocate(BINDEX)

  ! END SUBROUTINE CreatePCABins_Type1
  
  ! SUBROUTINE AppendStepBin(PCABin,Error,STEP_BIN,STEP_NEOF)

  !   ! --------------------
  !   ! subroutine arguments
  !   ! --------------------
  !   TYPE(PCABinType),       INTENT(INOUT) :: PCABin
  !   TYPE(ErrorType),        INTENT(INOUT) :: Error
  !   REAL(KIND=8),           INTENT(IN)    :: STEP_BIN
  !   INTEGER, OPTIONAL,      INTENT(IN)    :: STEP_NEOF

  !   ! ---------------
  !   ! local variables
  !   ! ---------------
  !   REAL(KIND=8), ALLOCATABLE :: STEP_BIN_TMP(:)
  !   INTEGER,      ALLOCATABLE :: STEP_NEOF_TMP(:)
  !   INTEGER                   :: NDIM

  !   ! For error checking
  !   CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'AppendStepBin'

  !   ! =====================================================================
  !   ! AppendStepBin begins here
  !   ! =====================================================================

  !   IF(CheckError(Error)) RETURN
    
    

  !   IF(ALLOCATED(PCABin%STEP_BINS)) THEN

  !     ! Get original dimension
  !     NDIM = SIZE(PCABin%STEP_BINS)

  !     ! Store original array
  !     ALLOCATE(STEP_BIN_TMP(0:NDIM-1))
  !     STEP_BIN_TMP = PCABin%STEP_BINS

  !     ! Increase dimension by one
  !     DEALLOCATE(PCABin%STEP_BINS)
  !     ALLOCATE(PCABin%STEP_BINS(0:NDIM))

  !     ! Populate array
  !     PCABin%STEP_BINS(0:NDIM-1) = STEP_BIN_TMP
  !     PCABin%STEP_BINS(NDIM) = STEP_BIN
  !     print*,'NDIM',NDIM

  !     IF(PRESENT(STEP_NEOF)) THEN

  !       IF(NDIM .EQ. 1) THEN
  !         ALLOCATE(PCABin%STEP_NEOF(1))
  !         PCABin%STEP_NEOF(1) = STEP_NEOF
  !         print*,'~~>',NDIM
  !       ELSE

  !         print*,'-->',NDIM

  !         ! Store original array
  !         ALLOCATE(STEP_NEOF_TMP(NDIM-1))
  !         STEP_NEOF_TMP = PCABin%STEP_NEOF

  !         ! Increase dimension by one
  !         DEALLOCATE(PCABin%STEP_NEOF)
  !         ALLOCATE(PCABin%STEP_NEOF(NDIM))

  !         ! Populate array
  !         PCABin%STEP_NEOF(1:NDIM-1) = STEP_NEOF_TMP
  !         PCABin%STEP_NEOF(NDIM) = STEP_NEOF

  !       ENDIF

  !     ENDIF

  !   ELSE
  !     ALLOCATE(PCABin%STEP_BINS(0:0))
  !     PCABin%STEP_BINS(0) = STEP_BIN
  !   ENDIF

  ! END SUBROUTINE AppendStepBin
  
END MODULE pca_bin_module