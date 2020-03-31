module GEMSTOOL_createbins_m

!  Modules in GEMSTOOL_sourcecode/structures

   USE GEMSTOOL_pars_m, Only : GTPK, GTZERO, GT_Maxbins, GT_maxlayers, GT_maxwav
   REAL    :: the_bin_size 
  INTEGER :: the_eof_num, which_binvar
  LOGICAL :: do_debug_bin = .false.
!  private  midval, GEMSTOOL_CreateBins_Rob_8bin_V3, GEMSTOOL_CreateBins_V3_HighEndOnly
!  public   GEMSTOOL_CreateBins_V0, GEMSTOOL_CreateBins_V1, GEMSTOOL_CreateBins_V2, GEMSTOOL_CreateBins_V3,&
! GEMSTOOL_CreateBIns_V4, 
  PUBLIC :: GEMSTOOL_CreateBins_V5

contains

!------------------------------------
! Jbak, 2020-01
! @ input variables
! which_win : put anyvalue excpet zero
! ndat      : number of wavelengths
! hgascum(ndat)   : H2O + O2 ODS
! gascum (ndat)    : total gas ODS
! wav (ndat)       : wavelengths
! @ Output variables
! index( ndat): wavindex as a functino  of bin
! bin (ndat) : binindex as a function of wav
! nbin : number of bin
! ncnt : number of sampling w.r.t each bin
! neof : number of EOFs w.r.t each bin
! binlimes : boundaries of each bins
!                 
!------------------------------------
SUBROUTINE GEMSTOOL_CreateBins_V5 ( which_win, ndat, hgascum, gascum, wav, amf,& !INPUT
             index, bin, nbin, ncnt, neof, binlims) ! OUTPUT
   USE GEMSTOOL_binset
   IMPLICIT NONE
!-----------------------------------------------
!  INPUT
!  ---------------------------------------------
   INTEGER, INTENT (IN)    :: NDAT, which_win
   INTEGER, INTENT (OUT)   :: NBIN
   REAL (GTPK), INTENT(IN) :: amf
   REAL (GTPK), DIMENSION (ndat), INTENT(IN) :: WAV ! used for binning
   REAL (GTPK), DIMENSION (ndat), INTENT(IN) :: hgascum,gascum ! used for binning
!-----------------------------------------------
!  OUTPUT
!  ---------------------------------------------
   INTEGER , INTENT(OUT), DIMENSION (ndat) :: BIN
   INTEGER , INTENT(OUT), DIMENSION (ndat) :: INDEX
   INTEGER , INTENT(OUT)    :: NCNT(0:GT_Maxbins), NEOF(0:GT_Maxbins)
   REAL(GTPK), INTENT(OUT)  :: BINLIMS (0:GT_Maxbins)
!----------------------------------------------
!  local variables
!  -------------------------------------------
   INTEGER :: i,j, fidx, lidx,IBIN, W, ipca, COUNT, K,STEP_INDEX, & 
             cc, it, count1, W1, W2, off_bin,ntau
   REAL    :: dtau, max_tau, min_tau, mtau,sd,dw
   INTEGER,DIMENSION(:), ALLOCATABLE :: BINDEX  !(ndat)
   REAL,DIMENSION (:), ALLOCATABLE   :: TAU
   REAL,DIMENSION (:), ALLOCATABLE   :: TAU_sub !(ndat)

   ! PCA windows
   INTEGER :: mpca, npca,nuv,nvis
   REAL (GTPK),DIMENSION (0:50):: pca_wav ! wavlimt for pcawin
   INTEGER :: PCA_wpix(50,2)  ! index of for each pcawin
   INTEGER :: PCA_nbin(50) ! N of bins for each pcawin
   REAL,DIMENSION (50,0:GT_Maxbins) :: PCA_BINS ! binlimt for each pcawin
   INTEGER,    DIMENSION (50,GT_Maxbins) :: PCA_NEOFS !Neof for each pcawin
   TYPE (binset) :: iter
   !----------------------
   ! Initialize
   !----------------------
   allocate(tau(ndat), bindex(ndat))
   !----------------------
   ! binning setting for multiple windows and binning parameters w.r.t UV or VIS
   !------------------------
   IF (which_win == 0) THEN ! testing 
    mpca=1; pca_wav(0:mpca) = (/minval(wav), maxval(wav)/)
    npca=1; pca_wpix(1,1:2) = (/1, ndat/)
    IF (wav(1) < 500) THEN 
       TAU(1:ndat) = -alog(REAL(gascum(1:ndat), KIND=4))
    ELSE
       TAU(1:ndat) = REAL(hgascum(1:ndat), KIND=4)
    ENDIF
   ELSE
     IF (wav(1) < 500) THEN 
      TAU(1:ndat) = -alog(REAL(gascum(1:ndat), KIND=4))
      mpca = 8
      pca_wav(0:mpca) = (/250.,340.,350.0,360.0,390.,420.,460.,490.,500./)
     ELSE
      TAU(1:ndat) = hgascum(1:ndat) ! gascum for absorption lines here they are H2o + O2
      mpca  = 13
      pca_wav(0:mpca) = (/500.,520., 540.,550.,565.,595.,625.,650.,686.75, 690., 705.,716.0,735.,752.0/)
     ENDIF
   ENDIF
   ! extract PCA windows limited for the selected wavelengths (wav)
   IF (which_win /=0) THEN 
    w1 = 1 ; npca=0
    !print * , 'nwav:', ndat
    DO ipca = 1, mpca
      w2 = MINVAL(MAXLOC( wav(1:ndat),MASK=(wav(1:ndat) < pca_wav(ipca))))
     ! WRITE(*,'(A, 2f8.2, 3i6, 3f8.2)') "loop of im:",wav(w1),wav(w2), w1, w2, w2-w1+1
      cc = w2-w1+1
      IF (cc > 1) THEN 
        npca=npca + 1
        PCA_wpix(npca,1) = w1 ;PCA_wpix(npca,2) = w2
        w1 = w2 +1   
      ELSE IF (npca > 1) THEN 
        PCA_wpix(npca, 2) = w2
      ENDIF
     ENDDO
   ENDIF

   ! first loop : 1 to mwav(number of pca_wav)
   ! second loop : 1 to nbin (number of bin for each pcawav)
   w1 = 1
   pca_nbin(:) = 0 
   DO ipca = 1, npca
     ! set the wavelength range for each sub window
     w1 = PCA_wpix(ipca,1)  ;w2=PCA_wpix(ipca,2)
     ntau = w2-w1 + 1
     allocate (tau_sub(ntau))
     tau_sub = tau(w1:w2)
     min_tau = minval(tau_sub) ; max_tau=maxval(tau_sub)+1E-20
     mtau    = sum(tau_sub)/ntau
     sd      = sum((tau_sub - mtau)*(tau_sub-mtau))/ntau
     sd      = sqrt(sd)
     IF (which_win .ne. 0) THEN 
     call get_bin(real(wav(w1),kind=4), real(wav(w2), kind=4),& 
          min_tau,max_tau,mtau,sd,real(amf, kind=4),iter)
     ELSE
      CALL bin_UV (real(min_tau, kind=4),real(max_tau, kind=4),& 
      real(the_bin_size, kind=4), the_eof_num,iter) ! o3 dominant
      CALL bin_UV (real(min_tau, kind=4),real(max_tau, kind=4),& 
      real(the_bin_size, kind=4), the_eof_num,iter) ! o3 dominant
     ENDIF
     k = 1; nbin=0
     PCA_bins(ipca,:) = min_tau
     
     !WRITE(*,'(A, 2f8.2, 3i5, 3f8.2)') "loop of im:",& 
     !wav(w1),wav(w2), w1, w2, ntau, min_tau, max_tau,mtau
     DO it = 1, iter%nbin
       DO w = 1, GT_Maxbins
         PCA_bins(ipca,k)  = iter%dtau(it) + min_tau
         IF (PCA_bins(ipca,k) > iter%maxtau(it) ) THEN 
             PCA_bins(ipca,k) = iter%maxtau(it)
         ENDIF
         min_tau           = PCA_bins(ipca,k)
         PCA_NEOFS(ipca,k) = iter%neof(it) 
         cc = COUNT(mask = (tau_sub(1:ntau) >= PCA_BINS(ipca,K-1) .and. &
                            tau_sub(1:ntau) < PCA_BINS(ipca,k) ))
         !IF (cc> 0) WRITE(*,'(4i3, i5, 199f8.2)') & 
         !   ipca,it,w, k, cc, PCA_BINS(ipca,k-1), PCA_BINS(ipca,k), &
         !   iter%dtau(it), min_tau, iter%maxtau(it)
         IF (cc <= iter%cc(it) ) THEN 
           IF (cc >= 1 ) THEN
             IF (cc ==1) THEN 
                 PCA_NEOFS(ipca,k) = 0
                 K = K + 1
             ELSE IF (cc > 1 .and. cc <= 3) THEN 
             !    PCA_NEOFS(ipca,k) = 2
             !    K = K + 1
             ENDIF
           ELSE IF  (cc == 0) THEN
             k = k 
           ENDIF
         ELSE 
           k = k +1
         ENDIF
         IF (min_tau >= iter%maxtau(it)) EXIT
         IF (min_tau >= maxval(tau_sub)) EXIT
       ENDDO ! loop of bin for each iteration       
     ENDDO ! loop of iteration for each window
     nbin = k -1
     IF (nbin ==0) THEN 
        PRINT *, 'nbin ==0' , nbin, cc,wav(w1), wav(w2)
     ENDIF
     PCA_nbin (ipca) = nbin
     IF (PCA_BINS(ipca,0) < minval(tau_sub)) PCA_BINS(ipca,0) = minval(tau_sub)
     IF (PCA_BINS(ipca,nbin) < max_tau) PCA_BINS(ipca,nbin) = max_tau
     deallocate (tau_sub)
     IF (nbin > GT_Maxbins) THEN 
       WRITE(*,*) 'nbin > GT_Maxbins' ; stop
     ENDIF
   ENDDO ! loop of win   
   ! Binning process
   IBIN = 0  ; NCNT = 0 ; NEOF = 0 ; BINLIMS = 0.0 ; BIN = 0 ; BINDEX = 0
   OFF_BIN= 0
   DO ipca = 1, npca
     w1 = PCA_wpix(ipca,1) ; w2=PCA_wpix(ipca,2) ; nbin = PCA_nbin(ipca)
     IF (nbin == 0) cycle
     !WRITE(*,'(A,2f8.2, 4i5)') 'w1:w2=' , wav(w1), wav(w2),w2-w1+1,w1,w2, nbin
     DO w = w1, w2
       dtau = tau(w)
       STEP_INDEX = MINVAL(MAXLOC( PCA_BINS(ipca,1:nbin),MASK=(PCA_BINS(ipca,1:nbin) < dtau)))
       IF (dtau <= PCA_BINS(ipca,1)) THEN
         IBIN = OFF_BIN
         STEP_INDEX = 1 !@bottom bin
       ELSE IF (dtau > PCA_BINS(ipca,nbin)) THEN 
         IBIN = NBIN + OFF_BIN
         STEP_INDEX = NBIN !@Top bin
         PCA_BINS(STEP_INDEX, ipca) = DTAU
       ELSE IF (dtau > PCA_BINS(ipca,1) .and. dtau <= PCA_BINS(ipca,nbin)) THEN
         IBIN = STEP_INDEX+ OFF_BIN
         STEP_INDEX = STEP_INDEX + 1
       ENDIF
  
       NCNT(IBIN) = NCNT(IBIN)+1
       NEOF(IBIN) = PCA_NEOFS(ipca,STEP_INDEX)
       BINLIMS(IBIN)  = PCA_BINS(ipca,STEP_INDEX)
       BIN(W)     = IBIN
       BINDEX(W)  = NCNT(IBIN)
      ! if (w == w1 .or. w == w2) print *, bin(w), bindex(w), step_index, nbin
     ENDDO     
     off_bin = off_bin + nbin
   ENDDO
   NBIN = OFF_BIN

   IF (do_debug_bin) THEN     
     PRINT * , 'tau:', minval(tau), maxval(tau),NDAT
     DO IBIN = 0, NBIN-1
       PRINT * , IBIN, NCNT(IBIN),SUM(NCNT(0:IBIN)), NEOF(IBIN), BINLIMS(IBIN)
     ENDDO
   ENDIF
   !  get the indices to change between bin space and wavenumber space
   DO W = 1, NDAT
     IF (BIN(W) .EQ. 0) THEN
         COUNT1 = BINDEX(W)
     ELSE
         COUNT1 = SUM(NCNT(0:BIN(W)-1)) + BINDEX(W)
     ENDIF
     IF (COUNT1 == 0) THEN 
        WRITE(*,*) 'Errors in check bin_v6' ; STOP
     ENDIF
     INDEX(COUNT1) = W
   ENDDO

   IF (sum(NCNT(0:nbin-1)) /= ndat) THEN 
     WRITE(*,*) 'ncnt=ndat', sum(ncnt(0:nbin-1)), ndat ; stop
   ENDIF
   ! Finish with deallocation
   deallocate(tau,bindex)
   RETURN

END SUBROUTINE GEMSTOOL_CreateBins_V5

subroutine GEMSTOOL_CreateBins_V3 &
            ( ndat, nlay, nbin, binlims, gasdat, &
              ncnt_new, index_new, bin_new )

   IMPLICIT NONE

!  inputs
!  ------

!  Numbers
!     1/11/16 Flexible number of bins, NBIN isIntent(inout)

   integer, intent(in)    :: NLAY, NDAT
   integer, intent(inout) :: NBIN

!  Optical data

   real(GTPK), intent(in) :: gasdat (GT_Maxlayers,GT_MaxWav)

!  outputs
!  -------

!  1/11/16. Changed dimensioning on BINLIMS to symbolic 

   INTEGER   , intent(inout) :: BIN_NEW(NDAT)
   INTEGER   , intent(out)   :: NCNT_NEW(0:GT_Maxbins), INDEX_NEW(GT_MaxWav)
   real(GTPK), intent(out)   :: BINLIMS(0:GT_Maxbins)

!  local variables (Dynamic memory here!!)
!  ---------------

! JBak alteration, 8/2/18. Add nbin0, which_jbak, STEP_INDEX

      logical :: iterating
      INTEGER :: IBIN, W, COUNT, K, NBIN_OLD, NBIN_NEW, THRESH_COUNT, STEP_INDEX
      INTEGER :: NCNT_OLD(0:GT_Maxbins), INDEX_OLD(GT_MaxWav)
      INTEGER :: BIN_OLD(NDAT), BINDEX_OLD(NDAT), BINDEX_NEW(NDAT)
      REAL(GTPK) :: TAUTOT, DTAU
      INTEGER, PARAMETER :: nbin0=8, which_jbak = 2
      real(GTPK) :: REAL_INITIAL_BINLIMS(nbin0), INITIAL_BINLIMS(nbin0)
      Data REAL_INITIAL_BINLIMS / 0.01_GTPK, 0.1_GTPK, 0.213_GTPK, 0.467_GTPK, 1.0_GTPK, 2.13_GTPK, 4.67_GTPK, 10.0_GTPK /

!  Threshholds

!      Data THRESH_COUNT / 30 /     ! Set 1, 12 January 2016
      Data THRESH_COUNT / 30 /      ! Set 2, 12 January 2016

!  Initialize

      DO IBIN = 1, nbin0
         INITIAL_BINLIMS(IBIN) = LOG(REAL_INITIAL_BINLIMS(IBIN))
         BINLIMS(IBIN-1) = INITIAL_BINLIMS(IBIN)
      ENDDO
      NCNT_NEW = 0 ; BINLIMS = GTZERO ; INDEX_NEW = 0 ; BIN_NEW = 0

!  initial value of NBIN should always be 9

      NBIN_OLD = NBIN ; NCNT_OLD = 0
! JBak alteration, 8/2/18. Change condition which which_jbak = 0
!      if ( NBIN .ne.9 ) stop'bad boy, change input of NBIN, Brian'
      if ( NBIN0+1 .ne.9 .and. which_jbak == 0) stop 'bad boy, change input of NBIN, Brian'

!  NOTES. 29 May 2015
!  mick fix - changed argument passing on GASDAT from ":" to "1:nlay" at several locations
!             [e.g., GASDAT(:,1) to GASDAT(1:nlay,1)] 

!  start wavelength loop
      IBIN = 0
      DO W = 1, NDAT

!  basic input, each wavelength
! JBak alteration, 8/2/18. Use Step_index, which_jbak and nbin0

        TAUTOT = SUM(GASDAT(1:nlay,W)) ; DTAU = LOG(TAUTOT)
        STEP_INDEX = MINVAL(MAXLOC( INITIAL_BINLIMS(1:nbin0),MASK=(INITIAL_BINLIMS(1:nbin0) < DTAU)))
        IF (which_jbak == 1 ) THEN 
          IF (STEP_INDEX < nbin0) THEN 
             IBIN = STEP_INDEX; BINLIMS(IBIN) = INITIAL_BINLIMS(IBIN+1)
          ELSE IF (STEP_INDEX == nbin0) THEN 
             IBIN = nbin0;      BINLIMS(IBIN) = 5000.0_GTPK
          ENDIF
        ELSE
         IF (DTAU .LT. INITIAL_BINLIMS(1)) THEN
           IBIN = 0 ; BINLIMS(IBIN) = INITIAL_BINLIMS(1)
         ELSE IF (DTAU .LT. INITIAL_BINLIMS(2)) THEN
           IBIN = 1 ; BINLIMS(IBIN) = INITIAL_BINLIMS(2)
         ELSE IF ( DTAU .LT. INITIAL_BINLIMS(3)) THEN
           IBIN = 2 ; BINLIMS(IBIN) = INITIAL_BINLIMS(3)
         ELSE IF (DTAU .LT. INITIAL_BINLIMS(4)) THEN
           IBIN = 3 ; BINLIMS(IBIN) = INITIAL_BINLIMS(4)
         ELSE IF (DTAU .LT. INITIAL_BINLIMS(5)) THEN
           IBIN = 4 ; BINLIMS(IBIN) = INITIAL_BINLIMS(5)
         ELSE IF (DTAU .LT. INITIAL_BINLIMS(6)) THEN
           IBIN = 5 ; BINLIMS(IBIN) = INITIAL_BINLIMS(6)
         ELSE IF (DTAU .LT. INITIAL_BINLIMS(7)) THEN
           IBIN = 6 ; BINLIMS(IBIN) = INITIAL_BINLIMS(7)
         ELSE IF (DTAU .LT. INITIAL_BINLIMS(8)) THEN
           IBIN = 7 ; BINLIMS(IBIN) = INITIAL_BINLIMS(8)
         ELSE 
           IBIN = 8 
           BINLIMS(IBIN) = 5000.0_GTPK
         ENDIF
        ENDIF
        NCNT_OLD(IBIN) = NCNT_OLD(IBIN)+1
        BIN_OLD(W) = IBIN
        BINDEX_OLD(W) = NCNT_OLD(IBIN)
      ENDDO
   ! DO W = 0, nbin-1
   !  print * , W, exp(BINLIMS(W)), NCNT_OLD(W)
   ! ENDDO
!    STOP
!write(*,'(i2,2x,10i5)')NBIN_OLD, NCNT_OLD(0:8), SUM(NCNT_OLD(0:8))
!write(*,'(9F12.4)')BINLIMS(0:8)

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
        IF (BIN_OLD(W) .EQ. 0) THEN
          COUNT = BINDEX_OLD(W)
        ELSE
          COUNT = SUM(NCNT_OLD(0:BIN_OLD(W)-1)) + BINDEX_OLD(W)
        ENDIF
        INDEX_OLD(COUNT) = W
      ENDDO

!  debug before
!      write(*,*)'OLD NCNT, Sum/Details = ',SUM(NCNT_OLD(0:NBIN_OLD-1)),NCNT_OLD(0:NBIN_OLD-1)
!      do w = 1, ndat
!         write(77,*)W, LOG(SUM(GASDAT(1:nlay,W))), BINDEX_OLD(W), BIN_OLD(W), INDEX_OLD(W)
!      enddo

!  BIN REASSIGNMENT - BOTTOM END
!  =============================

!  continuation point

568   continue

!  Now reduce the number of bins by 1 (from the bottom), fills up next bin in sequence
      
      K = NBIN_OLD ; iterating = .false.
      if ( NCNT_OLD(0) .lt. THRESH_COUNT ) THEN
        ITERATING = .true.
        NBIN_NEW = NBIN_OLD - 1
        NCNT_NEW(1:K-2) = NCNT_OLD(2:K-1)
        NCNT_NEW(0) = NCNT_OLD(0)
      endif

!  Not iterating : copy old to new and exit

      if ( .not. iterating ) then
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
        GO TO 571
      endif

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

!  get the indices to change between bin space and wavenumber space

!      DO W = 1, NDAT
!        IF (BIN_NEW(W) .EQ. 0) THEN
!          COUNT = BINDEX_NEW(W)
!        ELSE
!          COUNT = SUM(NCNT_NEW(0:BIN_NEW(W)-1)) + BINDEX_NEW(W)
!        ENDIF
!        INDEX_NEW(COUNT) = W
!      ENDDO

!  Re-assign output

      NBIN = NBIN_NEW
      DO IBIN = 0, NBIN
        BINLIMS(IBIN) = BINLIMS(IBIN+1) 
      ENDDO
      BINLIMS(NBIN) = GTZERO 

!  debug after

!      write(*,*)'NEW NCNT, Sum/Details = ',SUM(NCNT_NEW(0:NBIN_NEW-1)),NCNT_NEW(0:NBIN_NEW-1)
!      do w = 1, ndat
!         write(78,*)W, BINDEX_NEW(W), BIN_NEW(W), LOG(SUM(GASDAT(1:nlay,W)))
!      enddo

!  Re-set and return to iteration

      NBIN_OLD = NBIN_NEW
      NCNT_OLD = 0
      NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
      BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
      BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

      go to 568
      
!  BIN REASSIGNMENT - TOP END
!  ==========================

571   continue

!  Re-set

      NBIN_OLD = NBIN_NEW
      NCNT_OLD = 0
      NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
      BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
      BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

!  Debug after reassignment at the bottom end
!      do w = 1, ndat
!         write(99,*)W, LOG(SUM(GASDAT(1:nlay,W))), BINDEX_OLD(W), BIN_OLD(W), INDEX_NEW(W)
!      enddo

!  continuation point

567   continue

!  Now reduce the number of bins by 1 (from the top), fills up penumltimate bin
!    LAST BIN ONLY
      
      K = NBIN_OLD-1 ; iterating = .false.
      if ( NCNT_OLD(K) .lt. THRESH_COUNT ) THEN
        ITERATING = .true.
        NBIN_NEW = NBIN_OLD - 1
        NCNT_NEW(0:K-2) = NCNT_OLD(0:K-2)
        NCNT_NEW(K-1) = NCNT_OLD(K-1) !+ NCNT_OLD(K)
      endif

!  Not iterating : copy Old to new and exit

      if ( .not. iterating ) then
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
        RETURN
      endif

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

!  get the indices to change between bin space and wavenumber space

!      DO W = 1, NDAT
!        IF (BIN_NEW(W) .EQ. 0) THEN
!          COUNT = BINDEX_NEW(W)
!        ELSE
!          COUNT = SUM(NCNT_NEW(0:BIN_NEW(W)-1)) + BINDEX_NEW(W)
!        ENDIF
!        INDEX_NEW(COUNT) = W
!      ENDDO

!  Re-assign output

      NBIN = NBIN_NEW
      BINLIMS(NBIN_NEW) = BINLIMS(NBIN_OLD) 
      BINLIMS(NBIN_OLD) = GTZERO 

!  debug after

!      write(*,*)'NEW NCNT, Sum/Details = ',SUM(NCNT_NEW(0:NBIN_NEW-1)),NCNT_NEW(0:NBIN_NEW-1)
!      do w = 1, ndat
!         write(78,*)W, BINDEX_NEW(W), BIN_NEW(W), LOG(SUM(GASDAT(1:nlay,W)))
!      enddo

!  Re-set and return to iteration

      NBIN_OLD = NBIN_NEW
      NCNT_OLD = 0
      NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
      BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
      BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

      go to 567
      
!  finish

   return
end subroutine GEMSTOOL_CreateBins_V3

!------------------------------------
! Jbak, 2018-07-20
! New binnings with size of wavelength
!------------------------------------
SUBROUTINE GEMSTOOL_CreateBins_V4 ( which_win, ndat, wav, &
             index, bin, &
             nbin, ncnt, neof, binlims)

   IMPLICIT NONE
!-----------------------------------------------
!  INPUT
!  ---------------------------------------------
   INTEGER, INTENT (IN)    :: NDAT, which_win
   INTEGER, INTENT (OUT)   :: NBIN
   REAL (GTPK), DIMENSION (ndat), INTENT(IN) :: wav ! used for binning
!-----------------------------------------------
!  OUTPUT
!  ---------------------------------------------
   INTEGER , INTENT(OUT), DIMENSION (ndat) :: BIN
   INTEGER , INTENT(OUT), DIMENSION (ndat) :: INDEX
   INTEGER , INTENT(OUT)    :: NCNT(0:GT_Maxbins), NEOF(0:GT_Maxbins)
   REAL(GTPK), INTENT(OUT)  :: BINLIMS (0:GT_Maxbins)
!----------------------------------------------
!  local variables
!  -------------------------------------------
   INTEGER :: IBIN, W, COUNT, K,STEP_INDEX, cc, it, count1
   INTEGER :: THRESH_COUNT, NITER
   REAL (GTPK) :: TAUTOT, DTAU,  max_tau, min_tau
   INTEGER,    DIMENSION(:), ALLOCATABLE :: BINDEX !(ndat)
   REAL (GTPK),DIMENSION(:), ALLOCATABLE :: TAU !(ndat)
   REAL (GTPK),DIMENSION (0:GT_Maxbins) :: STEP_BINS 
   INTEGER,    DIMENSION (GT_Maxbins) :: STEP_NEOFS 
   INTEGER,    PARAMETER :: miter  = 10
   INTEGER,    DIMENSION (miter) :: iter_teof
   REAL (GTPK),DIMENSION (miter) :: iter_maxtau ,iter_dtau
   !----------------------
   ! Initialize
   !----------------------
   THRESH_COUNT = 20
   allocate(tau(ndat), bindex(ndat))
   !-----------------------
   ! decide binning steps
   !------------------------
   
   niter = 10 ; min_tau = 260
   iter_maxtau(1:niter) = (/270., 280. ,290. ,295. ,301. ,313. ,323. ,335.,340., 360./)
   iter_dtau (1:niter) = (/  20.,  10.,  10., 5.0,  3.0,   2.,  2.5,  2.,  5., 5.0/) !182
   iter_teof (1:niter) = (/  1,     1,   1,    1,    1,    1,    3,   2 ,  2,  2/)
   
   IF (which_win == 2) THEN  ! visible channel
    niter = 2 ; min_tau = 530
    iter_maxtau(1:niter) = (/540.0, 660.0/)
    iter_dtau (1:niter) = 10
    iter_teof (1:niter) = 4
   ENDIF
   IF (which_win == 0) THEN 
    niter = 2 ; min_tau = minval(wav)-10
    iter_maxtau(1:niter) = (/minval(wav), maxval(wav)/)
    iter_dtau (1:niter) = the_bin_size
    iter_teof (1:niter) = the_eof_num
   ENDIF
   STEP_BINS(:) = min_tau 
   k = 1 
   DO it = 1, niter
     DO w = 1, GT_Maxbins
       IF (k == 1) THEN
         STEP_BINS(k) = iter_maxtau(it)
       ELSE 
         STEP_BINS(k) =  w*iter_dtau(it) + min_tau
       ENDIF
       STEP_NEOFS(k) = iter_teof(it) 
       nbin = k  
       cc = COUNT(mask = (wav(1:ndat) >= STEP_BINS(K-1) .and. wav(1:ndat) < STEP_BINS(K)))
       !WRITE(*,'(i3, i5, 10f8.5)')  k, cc, STEP_BINS(K-1), STEP_BINS(K) , STEP_BINS(k) -STEP_BINS(K-1)
 
       IF (STEP_BINS(k) >= iter_maxtau(it)) exit
       IF (CC < THRESH_COUNT ) THEN 
         k = k 
       ELSE 
         k = k +1
       ENDIF
       
     ENDDO
     min_tau = step_bins(k) 
     IF (cc > THRESH_COUNT) THEN 
       k = k + 1
     ENDIF
   ENDDO

   IF (CC < THRESH_COUNT) THEN 
     step_bins(nbin-1) = MAXVAL(wav)
     nbin = nbin -1
   ENDIF

   IF (nbin > GT_Maxbins) THEN 
     WRITE(*,*) 'nbin > GT_Maxbins' ; stop
   ENDIF

   ! Binning process
   IBIN = 0  ; NCNT = 0 ; NEOF = 0 ; BINLIMS = 0.0 ; BIN = 0 ; BINDEX = 0
   DO W = 1, NDAT
     dtau = wav(w)
     STEP_INDEX = MINVAL(MAXLOC( STEP_BINS(1:NBIN),MASK=(STEP_BINS(1:NBIN) < dtau)))
     IF (dtau <= STEP_BINS(1)) THEN
       IBIN = 0 ; STEP_INDEX = 1 !@bottom bin
     ELSE IF (dtau > STEP_BINS(NBIN)) THEN 
       IBIN = NBIN
       STEP_INDEX = NBIN
       STEP_BINS(STEP_INDEX) = DTAU !@top bin
     ELSE IF (dtau > STEP_BINS(1) .and. dtau <= STEP_BINS(NBIN)) THEN
       IBIN = STEP_INDEX
       STEP_INDEX = STEP_INDEX + 1
     ENDIF
  
     NCNT(IBIN) = NCNT(IBIN)+1
     NEOF(IBIN) = STEP_NEOFS(STEP_INDEX)
     BINLIMS(IBIN)  = STEP_BINS(STEP_INDEX)
     BIN(W)     = IBIN
     BINDEX(W)  = NCNT(IBIN)
    !print * ,dtau,  IBIN, BINLIMS(IBIN), NCNT(IBIN)
   ENDDO

   IF (do_debug_bin) THEN    
     DO IBIN = 0, NBIN-1
       PRINT * , IBIN, NCNT(IBIN), NEOF(IBIN), BINLIMS(IBIN)
     ENDDO

   ENDIF
   !  get the indices to change between bin space and wavenumber space
   DO W = 1, NDAT
     IF (BIN(W) .EQ. 0) THEN
         COUNT1 = BINDEX(W)
     ELSE
         COUNT1 = SUM(NCNT(0:BIN(W)-1)) + BINDEX(W)
     ENDIF
     INDEX(COUNT1) = W
   ENDDO
   ! Finish with deallocation
   deallocate(tau, bindex)
   RETURN

END SUBROUTINE GEMSTOOL_CreateBins_V4

!  finish module

End Module GEMSTOOL_createbins_m

