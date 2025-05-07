module GEMSTOOL_PCACaller_m

!  Modules in GEMSTOOL_sourcecode/structures

   USE GEMSTOOL_pars_m
   USE GEMSTOOL_Input_types_m
   USE GEMSTOOL_Geophys_types_m

   USE GEMSTOOL_PCAPROJ_type_m            ! Output

!  This is the Regular PCA Caller for the Radiance Applicaton (No linearizations)

!  History of PCA developments for GODFIT
!  --------------------------------------

!   Version 1  : With original    EOF Scheme, 22 November 2011
!   Version 2  : With Alternative EOF Scheme, 26 January  2012
!   Version 2a : With Alternative EOF Scheme + CP-Linearization, 06 February  2012
!   Version 3  : With Alternative Second-Bin Scheme, 24 APril  2012 (Version 3)
!   Version 4  : With Option to do PCA on albedos for Internal Closure, 24 May 2012
!   Version 5  : With Option to do PCA on albedos for Internal Closure, 29 May 2012
!   Version 6  : With Profile linearization fully working, 05 November 2012
!   Version 7  : Robustness and Recoding, 08 November 2012

!  History of PCA developments for CALTECH
!  --------------------------------------

!   Mark 1  : With original    EOF Scheme, 22 November 2011
!   Mark 2  : With Alternative EOF Scheme, 26 January  2012
!   Mark 3  : With Alternative Second-Bin Scheme, 04 June 2012
!   Mark4/5 : (absent, as for Mark 3)
!   Mark 6  : Re-engineered to follow GODFIT Version 7, 05 April 2013
!               - Allows for albedo in the PCA, but not interpolation scheme
!               - clouds removed, changed indexing on AERMOMS
!               - Multiple aerosol types, 18 September 2013 (V. Natraj)
!               - Phase function projections introduced and aerosol moments 
!                 projected only up to 2*nstreams, 21 October 2013 (V. Natraj)


!  ---------------------------------------------------------------
!  This Version for GEMSTOOL, assembled Septebmer 2017 by R. Spurr
!  ---------------------------------------------------------------

!  Auxiliary routines

use pca_auxiliaries_m, only : PCA_LINTP2, PCA_Ranker

!  PCA routines
!    1/20/16, Introduce new solvers with Solar spectrum included.

use pca_svdcmpsolver_m
use pca_eigensolvers_m
use pca_eigensolvers_IncSun_m

public

contains

subroutine GEMSTOOL_PCACaller &
      ( do_aerosol, do_Sun_Normalized, use_hitran,        & ! Input Flags  
        do_svd_cmp, alb_pcainclude, B2S_index, k, n_eofs, & ! Input PCA Control  
        nlayers, nlayers2, nlayers21, nlayers22, nstr2,   & ! Input numbers
        nmuller, ndat, npoints, istart, index, Geophys,   & ! Binning and Optical Input
        PCAProj, PrinComps,                               & ! outputs
        pca_fail, pca_message, pca_trace_1, pca_trace_2 )   ! Exception handling

!  1/20/16. Introduce control for including solar spectrum in the PCA

   IMPLICIT NONE

!  inputs
!  ======

!  Flags
!  -----

!  Cloud and aerosol control

   logical, intent(in)     :: do_aerosol

!  Whether to use HITRAN or not

   logical, intent(in)     :: use_hitran
!  1/20/16 Sun normalization condition
!          If not set, will include solar spectrum in PCA process

   logical, intent(in)     :: Do_Sun_Normalized

!  PCA Control
!  -----------

!  Choice for the PCA solver

   logical, intent(in)     :: do_svd_cmp

!  Choice to include albedo in the PCA solver

   logical, intent(in)     :: alb_pcainclude

!  Second-bin strategy, PCA control number

   integer, intent(in)     :: B2S_index
   INTEGER, intent(in)     :: k, n_eofs

!  Choice to interpolate albedo using projected PCA results. Not included 1/20/16
!   logical, intent(in)     :: alb_pcainterp

!  Actual numbers
!  --------------
!  --> 1/20/16  revision to include nlayers22

   INTEGER, intent(in)     :: nlayers, nlayers2, nlayers21, nlayers22
   integer, intent(in)     :: ndat, npoints, nstr2, nmuller

!  offset and index for Given Bin of data

   integer      , intent(in) :: index(GT_MaxWav) 
   integer      , intent(in) :: istart

!  GEMSTOOL Type structure Geophysics (Input)

   TYPE(GEMSTOOL_Geophys_Setups), intent(in) :: Geophys

!  OUTPUTS
!  =======

!  PCA Projections. Type structure must now be Intent(InOut)

   TYPE(GEMSTOOL_PCAProj_Optical), intent(inout) :: PCAProj

!  Principal components (This should be Pre-allocated)
!  Global array, must now be Intent(InOut). 06 May 2015

   Real(GTPK), intent(inout) :: PrinComps(GT_MaxEofs,ndat)     !  Valid for Both solvers

!  Exception handling

   CHARACTER*(*), intent(out) :: PCA_MESSAGE, PCA_TRACE_1, PCA_TRACE_2
   LOGICAL      , intent(out) :: PCA_FAIL

!  LOCAL VARIABLES
!  ---------------

!  Local NoSun flag

   logical       :: do_NoSun

!  Number of 2 *EOFS + 1

   integer       :: n_eofs_2p1

!  Choice to include aerosols in the PCA solver; set based on B2S_index

   logical       :: aer_pcainclude

!  Local arrays in this section were formerly allocatable
!   -------Using dynamic memory allocation now.

   real(GTPK) :: tau_bin     (nlayers, npoints)
   real(GTPK) :: second_bin  (nlayers, npoints)
   real(GTPK) :: fr_bin      (nlayers, npoints)
   real(GTPK) :: reflec_bin  (npoints)
   real(GTPK) :: depol_bin   (npoints)
   real(GTPK) :: lambdas_bin (npoints)

   real(GTPK) :: fa_bin      (nlayers, npoints)
   real(GTPK) :: aersca_bin  (nlayers, npoints) 
   real(GTPK) :: aermoms_bin (6,0:GT_maxaermoms,npoints)

!  1/20/16. Introduce control for including solar spectrum in the PCA

   real(GTPK) :: solar_bin   (npoints)

!  Averaged properties in the old code

   Real(GTPK) :: FR_AVG, FA_AVG
   Real(GTPK) :: DEPOL_AVG, ALBEDO_AVG, LAMBDAS_AVG
   Real(GTPK) :: AERSCA_AVG(nlayers)

!  Second bin output

   real(GTPK) :: SEC_PRJ(nlayers,GT_MaxEofs2p1)

!  Aerosol bin output

   real(GTPK) :: AOD_PRJ(GT_MaxEofs2p1)
   real(GTPK) :: AER_PRJ(nlayers,GT_MaxEofs2p1)

!  Interpolation arrays
!   real(GTPK) :: scadep(npoints),scadep1(npoints),scalam(npoints)
!   integer       :: order(npoints)
!   real(GTPK) :: scaprj(GT_MaxEofs2p1), lamprj(GT_MaxEofs2p1)

!  Outputs from EIGENSOLVERS: pca_eigensolver, pca_eigensolver_alb, pca_eigensolver_aer

!  1/20/16. Introduce control for including solar spectrum in the PCA
!           Additional arrays to deal with solar inputs

   Real(GTPK) :: Atmosmean(nlayers,2)            ! Both Solvers
   Real(GTPK) :: Albmean                         ! +Albedo solver only
   Real(GTPK) :: Aodmean                         ! +Aerosol solver only
   Real(GTPK) :: SolarMean                       ! Solvers with Solar

   Real(GTPK) :: Eofs    (GT_MaxEofs,nlayers21)    ! Regular solver with SUN
   Real(GTPK) :: Eofs_Alb(GT_MaxEofs,nlayers22)    ! Albedo  solver with SUN
   Real(GTPK) :: Eofs_NS    (GT_MaxEofs,nlayers2)  ! Regular solver No SUN
   Real(GTPK) :: Eofs_NS_Alb(GT_MaxEofs,nlayers21) ! Albedo  solver No SUN
   Real(GTPK) :: Eofs_Aer(GT_MaxEofs,nlayers2+5)   ! Reg/NoSun+Aerosol solver

!  *** Conversion using Local PrinComps (06 May 2015)

   Real(GTPK) :: PrinComps_Local (n_eofs,npoints)

!  Outputs  from the SVDCMP_SOLVER

   Real(GTPK) :: Eofs_SVD      (nlayers2,nlayers2)
   Real(GTPK) :: PrinComps_SVD (nlayers2,npoints)

!  local variables

   real(GTPK) :: DDIM, TOTSCA, RAY2_AVG, AER_AVG
   real(GTPK) :: aod_bin  (npoints), AERSCA_TOT, ssafinal
   INTEGER    :: nf, N, M, MM, MP, N1, L, indexl, K1, local_nmoms, mmoms

!   INTEGER    :: L1, L2, indexl, Lstar
!   real(GTPK) :: s2s1, ss1, ww1, w2w1

!  2m + 1

   n_eofs_2p1 = 2 * n_eofs + 1

!  GTZERO indexing for incoming bin

   K1 = K + 1

!  Initialize exception handling

   pca_fail = .false. ; pca_message = ' '
   pca_trace_1 = ' '  ; pca_trace_2 = ' '

!  Initialize flag

   AER_PCAINCLUDE = (B2S_index .eq. 3)

!  Initialize output. Now done outside the routine

!   PrinComps = GTZERO
!   OPD_PRJ   = GTZERO   ; SSA_PRJ = GTZERO     ; FR_PRJ = GTZERO      ; FA_PRJ = GTZERO 
!   RAY2MOM_PRJ = GTZERO ; DEP_PRJ = GTZERO     ; AERMOMS_PRJ = GTZERO ; ALBEDO_PRJ = GTZERO
!   PF_PRJ = GTZERO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        Mapping Full grid of data to binned indexed grid
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Former Versions (1-6) used allocatable arrays

!   allocate ( tau_bin(nlayers,npoints), second_bin(nlayers,npoints), reflec_bin(npoints),   &
!              fr_bin(nlayers,npoints),  depol_bin(npoints) )
!   if ( do_aerosol ) then
!      allocate(fa_bin(nlayers,npoints), aersca_bin(nlayers,npoints), aercoeffs_bin(0:nmom_input,nlayers,6,npoints) )
!   endif
!   if ( do_cloud2 ) then
!      allocate(fc_bin(nlayers,npoints), cldsca_bin(nlayers,npoints), cldmoms_bin(0:2*nstr,npoints) )
!   endif

!  Bin common stuff

   do l = 1, npoints
      indexl = index(istart+l)
      tau_bin(1:nlayers,l)   = Geophys%TotalODs%taudp(1:nlayers,indexl)
      reflec_bin(l)          = Geophys%Surface%albedo(indexl)
      fr_bin(1:nlayers,l)    = Geophys%TotalODs%fr(1:nlayers,indexl)
      depol_bin(l)           = Geophys%Xsecs%Rayleigh_depol(indexl)
      if (.not. use_hitran) then
         lambdas_bin(l)      = Geophys%WavGrids%Wav(indexl)
      else
         lambdas_bin(l)      = Geophys%WavGrids%Wav(ndat-indexl+1)
      endif
      !print * , indexl, tau_bin(nlayers, l), reflec_bin(l), fr_bin(nlayers, l),depol_bin(l), lambdas_bin(l)
   enddo
!  1/20/16. Solar-bin inclusion

   if (.not.do_Sun_Normalized ) then
      do l = 1, npoints
         indexl = index(istart+l)
         solar_bin(l) = Geophys%SolarSpec%SunSpec(indexl)
      enddo
   else
      solar_bin = 1.0d0
   endif

!  Set NoSun flag

   Do_NoSun = do_Sun_Normalized

!  separate treatment of aerosols
!   - Moments, for now include them all, as we have to average later on...

   if ( do_aerosol ) then
      mmoms = -99
      do l = 1, npoints
         indexl = index(istart+l)
         fa_bin(1:nlayers,l)    = Geophys%TotalODs%fa(1:nlayers,indexl)
         mmoms = max(mmoms,Geophys%Aerosols%aerosol_nscatmoms(indexl))
         local_nmoms = GT_maxaermoms
         do nf = 1, nmuller
           aermoms_bin(nf,0:local_nmoms,l) = Geophys%Aerosols%aerosol_scatmoms(nf,0:local_nmoms,indexl)
         enddo
      enddo
   endif

!  Determine second bin according to strategy
!    Strategy 1, Second bin property is just single-scattering albedo
   if ( B2S_index .eq. 1 ) then
      do l = 1, npoints
         indexl = index(istart+l)
         Second_bin(1:nlayers,l) =  Geophys%TotalODs%omega(1:nlayers,indexl)
      enddo
   endif

!  Determine second bin according to strategy
!    Strategy 2, Second bin property is Rayleigh scattering optical depth
!                Assign also the aerosol scattering OD (if flagged)

   if ( B2S_index.eq.2 ) then
      do l = 1, npoints
         indexl = index(istart+l)
         do n = 1, nlayers
            Second_bin(n,l) =  Geophys%TotalODs%omega(n,indexl)* Geophys%TotalODs%fr(n,indexl)* Geophys%TotalODs%taudp(n,indexl)
         enddo
      enddo
      if ( do_aerosol ) then
         do l = 1, npoints
            indexl = index(istart+l)
            do n = 1, nlayers
               aersca_bin(n,l) = Geophys%TotalODs%omega(n,indexl)* Geophys%TotalODs%fa(n,indexl)*Geophys%TotalODs%taudp(n,indexl)
            enddo
         enddo
      endif
   endif

!  Determine second and third bins according to strategy
!    Strategy 3, Second bin property is Rayleigh scattering optical depth
!                Assign also the aerosol scattering OD (if flagged)
!                Third bin property is column aerosol OD for each aerosol type

   if ( B2S_index.eq.3 ) then
     do l = 1, npoints
       indexl = index(istart+l)
       do n = 1, nlayers
         Second_bin(n,l) = Geophys%TotalODs%omega(n,indexl)* Geophys%TotalODs%fr(n,indexl)*Geophys%TotalODs%taudp(n,indexl)
       enddo
     enddo
     if ( do_aerosol ) then
       do l = 1, npoints
         indexl = index(istart+l)
         do n = 1, nlayers
           aersca_bin(n,l) = Geophys%TotalODs%omega(n,indexl)* Geophys%TotalODs%fa(n,indexl)*Geophys%TotalODs%taudp(n,indexl)
         enddo
         aod_bin(l) = sum(aersca_bin(1:nlayers,l))
       enddo
     endif
   endif

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!             PRINCIPAL COMPONENT ANALYSIS Section
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Original "SVD" solution to the PCA
!  ----------------------------------

!    *** (Only applicable if there is no albedo analysis)
!    *** This routine is wasteful, it finds ALL the Eofs
!    *** We have to copy output for the first few Eofs
!    *** Convert to Global Princomps array (06 May 2015)

   if ( DO_SVD_CMP ) THEN
      IF ( .not. ALB_PCAINCLUDE ) then
         CALL pca_svdcmpsolver &
           ( npoints, nlayers, TAU_BIN, SECOND_BIN,     & ! Inputs
             ATMOSMEAN, Eofs_SVD, PrinComps_SVD)    ! Outputs
         DO m = 1, n_eofs
           Eofs(m,1:nlayers)  = Eofs_SVD(M,1:nlayers)
         enddo
         do l = 1, npoints
            indexl = index(istart+l)
            PrinComps(1:n_eofs,indexl) = PrinComps_SVD(1:n_eofs,l)
         ENDDO
      else
         pca_fail = .true.
         pca_message = 'Not allowed to use Original SVD EOF routine with albedos'
         pca_trace_1 = 'Action: Turn off DO_SVD_PCA if you want albedos included in the PCA'
         pca_trace_2 = 'pca_svdcmpsolver cannot be used here in pca_driver_V7'
         go to 679
      endif
   endif

!  "Eigensolution" solution to the PCA (New Default)
!  ------------------------------------------------

!    *** Convert to Global Princomps array (06 May 2015)

   if ( .not. DO_SVD_CMP ) then

      if ( ALB_PCAINCLUDE ) then
        if ( do_NoSun ) then
          CALL pca_eigensolver_alb &
           ( GT_MaxEofs, npoints, nlayers, nlayers21,            & ! Input dimensions (allocated here)
             n_Eofs, npoints, nlayers, nlayers2, nlayers21,    & ! Input control
             tau_bin, second_bin, reflec_bin,                  & ! Input optical properties
             Atmosmean, Albmean, Eofs_NS_Alb, PrinComps_Local, & ! Outputs
             pca_fail, pca_message, pca_trace_1 )                ! Exception handling
        else
          CALL pca_eigensolver_p1_alb &
           ( GT_MaxEofs, npoints, nlayers, nlayers22,                    & ! Input dimensions (allocated here)
             n_Eofs, npoints, nlayers, nlayers2, nlayers21, nlayers22, & ! Input control
             tau_bin, second_bin, reflec_bin, solar_bin,               & ! Input optical properties
             Atmosmean, Solarmean, Albmean, Eofs_Alb, PrinComps_Local, & ! Outputs
             pca_fail, pca_message, pca_trace_1 )                        ! Exception handling
        endif
      else

        if ( AER_PCAINCLUDE ) then

          CALL pca_eigensolver_aer &
           ( GT_MaxEofs, npoints, nlayers, nlayers21,       & ! Input dimensions (allocated here)
             n_Eofs, npoints, nlayers, nlayers21,           & ! Input control
             tau_bin, second_bin, aod_bin,                  & ! Input optical properties
             Atmosmean, Aodmean, Eofs_Aer, PrinComps_Local, & ! Outputs
             pca_fail, pca_message, pca_trace_1 )             ! Exception handling

        else

          if ( do_NoSun ) then
            CALL pca_eigensolver &
             ( GT_MaxEofs, npoints, nlayers, nlayers2, & ! Input dimensions (allocated here)
               n_Eofs, npoints, nlayers, nlayers2,     & ! Input control
               tau_bin, second_bin,                    & ! Input optical properties
               Atmosmean, Eofs_NS, PrinComps_Local,    & ! Outputs
               pca_fail, pca_message, pca_trace_1 )     ! Exception handling
          else
            CALL pca_eigensolver_p1 &
             ( GT_MaxEofs, npoints, nlayers, nlayers21,          & ! Input dimensions (allocated here)
               n_Eofs, npoints, nlayers, nlayers2, nlayers21,  & ! Input control
               tau_bin, second_bin, solar_bin,                 & ! Input optical properties
               Atmosmean, Solarmean, Eofs, PrinComps_Local,    & ! Outputs
               pca_fail, pca_message, pca_trace_1 )              ! Exception handling
          endif

        endif
      endif

!  Exception handling check, finish if set
!     Following options rewritten, 1/20/16

      if ( pca_fail ) then
        if ( ALB_PCAINCLUDE ) then
          if ( do_NoSun ) then
            pca_trace_2 = 'pca_eigensolver_alb (NoSun) failed in FastV2p7_Eofpc_PCAdriver'
          else
            pca_trace_2 = 'pca_eigensolver_p1_alb failed in FastV2p7_Eofpc_PCAdriver'
          endif
        else
          if ( AER_PCAINCLUDE ) then
            pca_trace_2 = 'pca_eigensolver_aer failed in FastV2p7_Eofpc_PCAdriver'
          else
            if ( do_NoSun ) then
              pca_trace_2 = 'pca_eigensolver (NoSun) failed in FastV2p7_Eofpc_PCAdriver'
            else
              pca_trace_2 = 'pca_eigensolver_p1 failed in FastV2p7_Eofpc_PCAdriver'
            endif
          endif
        endif
        go to 679
      endif

!    *** Convert to Global Princomps array

      do l = 1, npoints
         indexl = index(istart+l)
         PrinComps(1:n_eofs,indexl) = PrinComps_LOCAL(1:n_eofs,l)
      enddo
   endif

!   pause 'eof'

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     NEW TREATMENT: EOF-derived versus Averaged ALBEDOS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Set Projections from the PCA analysis
!  -------------------------------------

!  Mean atmospheric optical properties (same for ALL SOLVERS)

   DO N = 1,nlayers
      PCAProj%OPD_PRJ(N,K1,1) = EXP(Atmosmean(N,1))
      SEC_PRJ(N,1)    = EXP(Atmosmean(N,2))
   ENDDO

!  Atmospheric EOF +/- properties reconstructed from EOF perturbations
!    Separate EOF arrays for the albedo and aerosol alternatives !!!!
!    1/15/16. revision for No_Sun option.

   DO M = 1, n_eofs
      MP = 2*M ; MM = MP + 1
      if ( alb_pcainclude ) then
         if ( do_NoSun ) then
            DO N = 1, nlayers
               N1 = N + nlayers
               PCAProj%OPD_PRJ(N,K1,MP) = EXP(Atmosmean(N,1)+Eofs_NS_Alb(M,N))
               SEC_PRJ(N,MP)            = EXP(Atmosmean(N,2)+Eofs_NS_Alb(M,N1))
               PCAProj%OPD_PRJ(N,K1,MM) = EXP(Atmosmean(N,1)-Eofs_NS_Alb(M,N))
               SEC_PRJ(N,MM)            = EXP(Atmosmean(N,2)-Eofs_NS_Alb(M,N1))
            ENDDO
         else
            DO N = 1, nlayers
               N1 = N + nlayers
               PCAProj%OPD_PRJ(N,K1,MP) = EXP(Atmosmean(N,1)+Eofs_Alb(M,N))
               SEC_PRJ(N,MP)            = EXP(Atmosmean(N,2)+Eofs_Alb(M,N1))
               PCAProj%OPD_PRJ(N,K1,MM) = EXP(Atmosmean(N,1)-Eofs_Alb(M,N))
               SEC_PRJ(N,MM)            = EXP(Atmosmean(N,2)-Eofs_Alb(M,N1))
            ENDDO
         endif
      else
         if ( aer_pcainclude ) then
            DO N = 1, nlayers
               N1 = N + nlayers
               PCAProj%OPD_PRJ(N,K1,MP) = EXP(Atmosmean(N,1)+Eofs_Aer(M,N))
               SEC_PRJ(N,MP)            = EXP(Atmosmean(N,2)+Eofs_Aer(M,N1))
               PCAProj%OPD_PRJ(N,K1,MM) = EXP(Atmosmean(N,1)-Eofs_Aer(M,N))
               SEC_PRJ(N,MM)            = EXP(Atmosmean(N,2)-Eofs_Aer(M,N1))
            ENDDO
         else
            if ( do_NoSun ) then
               DO N = 1, nlayers
                  N1 = N + nlayers
                  PCAProj%OPD_PRJ(N,K1,MP) = EXP(Atmosmean(N,1)+Eofs_NS(M,N))
                  SEC_PRJ(N,MP)            = EXP(Atmosmean(N,2)+Eofs_NS(M,N1))
                  PCAProj%OPD_PRJ(N,K1,MM) = EXP(Atmosmean(N,1)-Eofs_NS(M,N))
                  SEC_PRJ(N,MM)            = EXP(Atmosmean(N,2)-Eofs_NS(M,N1))
               ENDDO
            else
               DO N = 1, nlayers
                  N1 = N + nlayers
                  PCAProj%OPD_PRJ(N,K1,MP) = EXP(Atmosmean(N,1)+Eofs(M,N))
                  SEC_PRJ(N,MP)            = EXP(Atmosmean(N,2)+Eofs(M,N1))
                  PCAProj%OPD_PRJ(N,K1,MM) = EXP(Atmosmean(N,1)-Eofs(M,N))
                  SEC_PRJ(N,MM)            = EXP(Atmosmean(N,2)-Eofs(M,N1))
               ENDDO
            endif
         endif
      endif
   ENDDO
 
!  1/20/16.  Set Solar Projection

    if ( .not. do_NoSUN ) then
      PCAProj%SOLAR_PRJ(K1,1) = EXP(Solarmean)
      if ( alb_pcainclude ) then
        DO M = 1, n_eofs
          MP = 2*M ; MM = MP + 1
          PCAProj%SOLAR_PRJ(K1,MP) = EXP(Solarmean+Eofs_Alb(M,nlayers21))
          PCAProj%SOLAR_PRJ(K1,MM) = EXP(Solarmean-Eofs_Alb(M,nlayers21))
        ENDDO
      else
        DO M = 1, n_eofs
          MP = 2*M ; MM = MP + 1
          PCAProj%SOLAR_PRJ(K1,MP) = EXP(Solarmean+Eofs(M,nlayers21))
          PCAProj%SOLAR_PRJ(K1,MM) = EXP(Solarmean-Eofs(M,nlayers21))
        ENDDO
      endif
    endif

!  Form average quantities
!  -----------------------

   ddim = 1.0_GTPK/dble(npoints)

!  Averaging for LAMBDAS

   LAMBDAS_AVG = SUM(LAMBDAS_BIN(:))*ddim

!  NOT GOING TO DO THIS FOR NOW.
!   if ( do_aerosol ) then
!      FAC2 = (LAMBDAS_AVG - LAMBDAS(1)) / (LAMBDAS(NDAT) - LAMBDAS(1))
!      FAC1 = 1.0d0 - FAC2
!        DO L = 0, 2*nstr
!          AERCOEFFS_PRJ(L,K1,1:6) = FAC1 * AERCOEFFS(1,L,NF,1:6) + FAC2 * AERCOEFFS(2,L,NF,1:6)
!        ENDDO
!      ENDDO
!      DO V = 1, NGEOMS
!        DO NF = 1, 5
!          PF_PRJ(NF,K1,V) = FAC1 * PHASFUNC(1,NF,V) + FAC2 * PHASFUNC(2,NF,V)
!        ENDDO
!      ENDDO
!   endif

!  Averaged coefficients

   if ( do_aerosol ) then
      PCAProj%N_AERCOEFFS_PRJ(K1,1:n_eofs_2p1) = mmoms
      DO L = 0, mmoms
        DO nf = 1, nmuller
          AER_AVG = SUM(aermoms_bin(nf,L,:))*ddim
          PCAProj%AERCOEFFS_PRJ(nf,L,K1,1:n_eofs_2p1) = AER_AVG
        ENDDO
      ENDDO
   endif

!  Averaging for FR and FA projections for strategy 1

   IF ( B2S_index .eq. 1 ) then
      DO N = 1, nlayers
         FR_AVG = SUM(FR_BIN(N,:))*ddim
         DO M = 1, n_eofs_2p1
            PCAProj%FR_PRJ(N,K1,M) = FR_AVG
         ENDDO
      ENDDO
      if ( do_aerosol ) then
         DO N = 1, nlayers
            FA_AVG = SUM(FA_BIN(N,:))*ddim
            DO M = 1, n_eofs_2p1
               PCAProj%FA_PRJ(N,K1,M) = FA_AVG
            ENDDO
         ENDDO
      endif
   ENDIF

!  For strategy 2 and 3, Averaged FR and FA projections come later, for now set AERSCA_AVG

   IF ( B2S_index .eq. 2 .or. B2S_index .eq. 3 ) then
      if ( do_aerosol ) then
         DO N = 1, nlayers
            AERSCA_AVG(N) = SUM(AERSCA_BIN(N,:))*ddim
         ENDDO
      ENDIF
   endif

!  Output optical properties
!  -------------------------

!  STRATEGY 1 : Just copy the second projection into SSA, also toggle ssa

   IF ( B2S_index .eq. 1 ) then
      DO M = 1, n_eofs_2p1
         DO N = 1, nlayers
            ssafinal = SEC_PRJ(N,M)
            IF (ssafinal  .GT. 0.999999_GTPK) ssafinal = 0.999999_GTPK ! consistent setting across code
            IF (ssafinal  .LT. 1.0e-6_GTPK)   ssafinal = 1.0e-6_GTPK   ! consistent setting across code
            PCAProj%SSA_PRJ(N,K1,M) = ssafinal
         ENDDO
      ENDDO
   ENDIF

!  STRATEGY 2: Must additionally define projections for FR and FA

   IF ( B2S_index .eq. 2 ) then
      if ( do_aerosol  ) then
         DO M = 1, n_eofs_2p1
            DO N = 1, nlayers
               TOTSCA = SEC_PRJ(N,M)
               TOTSCA = TOTSCA + AERSCA_AVG(N)
               PCAProj%FR_PRJ(N,K1,M)  = SEC_PRJ(N,M)  / TOTSCA
               PCAProj%FA_PRJ(N,K1,M)  = AERSCA_AVG(N) / TOTSCA
               PCAProj%SSA_PRJ(N,K1,M) = TOTSCA / PCAProj%OPD_PRJ(N,K1,M)
            ENDDO
         ENDDO
      else
         DO M = 1, n_eofs_2p1
            DO N = 1, nlayers         
               PCAProj%SSA_PRJ(N,K1,M) = SEC_PRJ(N,M) / PCAProj%OPD_PRJ(N,K1,M)
               PCAProj%FR_PRJ(N,K1,M)  = GTONE
            ENDDO
         ENDDO
      ENDIF
      DO M = 1, n_eofs_2p1
         DO N = 1, nlayers         
            ssafinal = PCAProj%SSA_PRJ(N,K1,M)
            IF (ssafinal  .GT. 0.999999_GTPK) ssafinal = 0.999999_GTPK ! consistent setting across code
            IF (ssafinal  .LT. 1.0e-6_GTPK)   ssafinal = 1.0e-6_GTPK   ! consistent setting across code
            PCAProj%SSA_PRJ(N,K1,M) = ssafinal
         ENDDO
      ENDDO
   endif

!  STRATEGY 3: Must additionally define projections for FR and FA
!  Aerosol projections are constructed as usual from AODMEAN and EOFs

   IF ( B2S_index .eq. 3 ) then
      if ( do_aerosol  ) then
         AOD_PRJ(1) = EXP(Aodmean)
         DO M = 1, n_eofs
            MP = 2*M ; MM = MP + 1
            AOD_PRJ(MP) = EXP(Aodmean+Eofs_Aer(M,nlayers21))
            AOD_PRJ(MM) = EXP(Aodmean-Eofs_Aer(M,nlayers21))
         ENDDO
         DO M = 1, n_eofs_2p1
            DO N = 1, nlayers
               TOTSCA = SEC_PRJ(N,M)
               AERSCA_TOT = SUM(AERSCA_AVG(:))
               AER_PRJ(N,M) = (AOD_PRJ(M)/AERSCA_TOT)*AERSCA_AVG(N)
               TOTSCA  = TOTSCA + AER_PRJ(N,M)
               PCAProj%FR_PRJ(N,K1,M)  = SEC_PRJ(N,M) / TOTSCA
               PCAProj%FA_PRJ(N,K1,M)  = AER_PRJ(N,M) / TOTSCA
               PCAProj%SSA_PRJ(N,K1,M) = TOTSCA / PCAProj%OPD_PRJ(N,K1,M)
            ENDDO
         ENDDO
      else
         DO M = 1, n_eofs_2p1
            DO N = 1, nlayers         
               PCAProj%SSA_PRJ(N,K1,M) = SEC_PRJ(N,M) / PCAProj%OPD_PRJ(N,K1,M)
               PCAProj%FR_PRJ(N,K1,M)  = GTONE
            ENDDO
         ENDDO
      ENDIF
      DO M = 1, n_eofs_2p1
         DO N = 1, nlayers         
            ssafinal = PCAProj%SSA_PRJ(N,K1,M)
            IF (ssafinal  .GT. 0.999999_GTPK) ssafinal = 0.999999_GTPK ! consistent setting across code
            IF (ssafinal  .LT. 1.0e-6_GTPK)   ssafinal = 1.0e-6_GTPK   ! consistent setting across code
            PCAProj%SSA_PRJ(N,K1,M) = ssafinal
         ENDDO
      ENDDO
   endif

!  Rayleigh moment
!   - For now, assign all values to average
!            (otherwise linearization gets complex)

   DEPOL_AVG   = SUM(DEPOL_BIN)*ddim
   RAY2_AVG    = ( GTONE - DEPOL_AVG ) / (2.0_GTPK + DEPOL_AVG  )
   DO M = 1, n_eofs_2p1
      PCAProj%RAY2MOM_PRJ(K1,M) = RAY2_AVG
      PCAProj%DEP_PRJ(K1,M) = DEPOL_AVG
   ENDDO

!  here is the code that could be used.......
!         call PCA_LINTP2(npoints,scalam,depol_bin,n_eofs_2p1,lamprj,dep_prj)
!         do m = 1, n_eofs_2p1
!            ray2mom_prj(m) = ( 1.0_GTPK - dep_prj(m) ) / (2.0_GTPK + dep_prj(m) )
!         enddo

!  Alternative scheme "FRFA_PCAINTERP"
!    develop Interpolated FR and FA projections corresponding to projected scattering optical depth

!    if ( frfa_pcainterp ) then
!      do l = 1, npoints
!         indexl = index(istart+l) ; scadep1(l) = dot_product(taudp(1:nlayers,indexl),omega(1:nlayers,indexl))
!      enddo
!      do m = 1, n_eofs_2p1 
!        scaprj(m) = dot_product(opd_prj(:,m),ssa_prj(:,m))
!      enddo
!      order = 0 ; call PCA_Ranker(npoints,scadep1,order)
!      do l = 1, npoints
!          scadep(l) = scadep1(order(l))
!      enddo
!      do m = 1, n_eofs_2p1 
!         l = 1
!         do while (scaprj(m).gt.scadep(l))
!           l = l + 1
!         enddo
!         l1 = l - 1; l2 = l1 + 1
!         s2s1 = scadep(l2) - scadep(l1) ; ss1 = scaprj(m) - scadep(l1)
!         DO N = 1, nlayers
!            do l = 1, npoints
!               Lstar = index(istart+order(l))
!               scafr(l) = fr(n,Lstar)
!               scafa(l) = fa(n,Lstar)
!            enddo
!            if ( FA_PRJ(N,M) .ne. 0.0d0 ) then
!               w2w1 = scafa(l2) - scafa(l1) ; ww1 = w2w1 * ss1 / s2s1
!               fa_prj(n,m) = ww1 + scafa(l1)
!               fr_prj(n,m) = 1.0d0 - fa_prj(n,m)        !alternative restriction, makes NO DIFFERENCE
!            endif
!            if ( FR_PRJ(N,M) .ne. 1.0d0 ) then
!               w2w1 = scafr(l2) - scafr(l1) ; ww1 = w2w1 * ss1 / s2s1
!               fr_prj(n,m) = ww1 + scafr(l1)
!            endif
!            write(68,*)istart,npoints,m,n,FR_PRJ(N,M),FA_PRJ(N,M),L1,L2
!         enddo
!      enddo
!   endif

!  Albedo schemes
!  ==============

!  If the Albedo is coupled to the atmospheric optical properties, then
!   Albedo projections are constructed as usual from ALBMEAN and EOFs
!    revision 1/15/16 to include Sun Option

   if ( alb_pcainclude ) then
      ALBEDO_AVG  = GTZERO
      PCAProj%ALBEDO_PRJ(K1,1) = EXP(Albmean)
      if ( do_NoSun ) then
        DO M = 1, n_eofs
          MP = 2*M ; MM = MP + 1
          PCAProj%ALBEDO_PRJ(K1,MP) = EXP(Albmean+Eofs_NS_Alb(M,nlayers21))
          PCAProj%ALBEDO_PRJ(K1,MM) = EXP(Albmean-Eofs_NS_Alb(M,nlayers21))
        ENDDO
      else
        DO M = 1, n_eofs
          MP = 2*M ; MM = MP + 1
          PCAProj%ALBEDO_PRJ(K1,MP) = EXP(Albmean+Eofs_Alb(M,nlayers22))
          PCAProj%ALBEDO_PRJ(K1,MM) = EXP(Albmean-Eofs_Alb(M,nlayers22))
        ENDDO
      endif
   endif

!  One Alternative is "ALBEDO AVERAGING". The only one for now,

   if ( .not.alb_pcainclude ) then
      ALBEDO_AVG  = SUM(REFLEC_BIN)*ddim
      DO M = 1, n_eofs_2p1
         PCAProj%ALBEDO_PRJ(K1,M) = ALBEDO_AVG
      ENDDO
   endif

!   COMMENTED OUT HERE ----------------------------------------------- 05  April 2013
!  Another Alternative scheme "ALBEDO_PCAINTERP"
!    develop Interpolated albedo projections corresponding to projected scattering optical depth

!    if ( .not.alb_pcainclude .and. alb_pcainterp ) then
!      do l = 1, npoints
!         indexl = index(istart+l)
!         scadep1(l) = dot_product(taudp(1:nlayers,indexl),omega(1:nlayers,indexl))
!      enddo
!      do m = 1, n_eofs_2p1 
!        scaprj(m) = dot_product(opd_prj(:,m),ssa_prj(:,m))
!      enddo
!      call PCA_Ranker(npoints,scadep1,order)
!      do l = 1, npoints
!          Lstar = index(istart+order(l))
!          scadep(l) = scadep1(order(l))
!          scalam(l) = lambdas(Lstar)
!      enddo
!      do m = 1, n_eofs_2p1 
!         l = 1
!         do while (scaprj(m).gt.scadep(l))
!           l = l + 1
!         enddo
!         l1 = l - 1; l2 = l1 + 1
!         s2s1 = scadep(l2) - scadep(l1) ; ss1 = scaprj(m) - scadep(l1) 
!         w2w1 = scalam(l2) - scalam(l1) ; ww1 = w2w1 * ss1 / s2s1
!         lamprj(m) = ww1 + scalam(l1)
!         albedo_prj(m) =  albcoeffs(1)
!         mult = GTONE ; fac = GTONE - (lamprj(m)/reflam)
!         do q = 2, nalbcoeffs
!            mult = mult * fac
!            albedo_prj(m) = albedo_prj(m) + mult*albcoeffs(q)
!         enddo
!!         write(*,'(i2,1pe24.12)')m,albedo_prj(m) 
!      enddo
!   endif

!  One Alternative is "ALBEDO AVERAGING". The only one for now,

!   if ( .not.alb_pcainclude .and..not. alb_pcainterp ) then
!      ALBEDO_AVG  = SUM(REFLEC_BIN)*ddim
!      DO M = 1, n_eofs_2p1
!         ALBEDO_PRJ(M) = ALBEDO_AVG
!      ENDDO
!   endif

!  Final section
!  -------------

!  Continuation point for failure output

679 continue

!  finish

   return
end subroutine GEMSTOOL_PCACaller

end module GEMSTOOL_PCACaller_m

