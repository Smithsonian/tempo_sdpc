module VLIDORT_eofpc_module

!  History.

!  V1: used VLIDORT for the Fo code. Too slow.

!  V2: used VLIDORT_NOFO and Stand-alone Optimized/Leveraged FO code
!  V2: used Optimized 2stream code.    

!  Juseon Bak fix : 8/2/18. V4 Binning Routine

!  Rob Fix 8/3/18. MultiGeometry capability properly installed
!                  Use of new 2S code.
!  Modules in GEMSTOOL_sourcecode/structures

   USE GEMSTOOL_pars_m
   !USE GEMSTOOL_Input_types_m
   !USE GEMSTOOL_Geophys_types_m
!   USE GEMSTOOL_Result_types_m
   USE GEMSTOOL_PCAPROJ_type_m

   !USE GEMSTOOL_L_Geophys_types_m
   USE GEMSTOOL_L_Result_types_m

!  VLIDORT (Version 2.7) Modules

   use VLIDORT_PARS
   use VLIDORT_IO_DEFS
   USE VLIDORT_LIN_IO_DEFS

   use VLIDORT_masters
   use VLIDORT_LPS_masters

!  Two-stream (Version 2.4) modules. Optimized Version

   use twostream_pars_m
   USE twostream_master_m
   USE twostream_lps_master_m

!  First Order (Version 1.4) modules. Optimized, Leveraged Version

   use FO_SSGeometry_Master_m
   use FO_VectorSS_spherfuncs_Optimized_m
   use FO_VectorSS_PhasMat_m
   use FO_VectorSS_PhasMat_Plus_m
   use FO_VectorSS_RTCalcs_I_Optimized_m
   use FO_VectorSS_RTCalcs_ILPS_Optimized_m

!  PCA-based modules

   USE pca_correction_m
   USE pca_correction_Plus_m

!  Initializer routine for VLIDORT
   USE GEMSTOOL_RTinitialize_m
!    - PCA Binning routine. Replaced 1/20/16.
!    - PCA Caller  routine (New version with polarization, 20 January 2016)
   use GEMSTOOL_PCACaller_m
   USE GEMSTOOL_createbins_pm
   ! JCH: because GEMSTOOL_createbins_m is private to vlidort, I'm changing
   !      the name of this module to GEMSTOOL_createbins_pm to indicate that
   !      this is a private/custom version

!  Jbak : shared variables between PCA master and O3P algorithm
   USE vlidort_data_module, ONLY: & 
       Inputs,& ! control input
       Geophys, L_Geophys,& ! optical input 
       layer_vary_number_cc, layer_vary_flag_cc, &
                                  which_win, & ! 1=uv or 2=vis
                                  RO, RI
   USE ozprof_data_module, ONLY: num_iter 
!  routine is public
   INTEGER, PARAMETER :: dir =1 ! 1=upwelling 2 = downwelling
   INTEGER, DIMENSION(:), ALLOCATABLE :: Assigned_bins
   REAL(KIND=dp), DIMENSION(GT_maxwav) :: absline, habsline
   INTEGER :: n_call_vlidort 
   public

contains

   subroutine VLIDORT_Eofpc_Master &
    ( Monitor_CPU, do_VLIDORT_initialize, do_Jacobians,      & ! Flags 
      fail, Nmessages, messages )      ! Status & Timings

   implicit none 

!  Inputs
!  ======

   logical, intent(in) :: Monitor_CPU ! monitoring flag
   logical, intent(inout) :: do_VLIDORT_initialize
   logical, intent(in) :: do_Jacobians

!  Stokes-vector and Jacobian output results
!  =========================================
 
!  Timing and status
!  =================

!  Exception handling variables

   logical, intent(inout)       :: Fail
   integer, intent(inout)       :: Nmessages
   character*(*), intent(inout) :: messages (GT_maxmessages)

!  Local variables
!  ===============

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      VLIDORT  ARGUMENTS (V2.7 code)
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), SAVE       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs),SAVE     :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut), SAVE          :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn

!  VLIDORT linearized supplements i/o structure

      TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup

!  VLIDORT linearized output structure

      TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut

!  Stokes-vector and Jacobian output results

!   TYPE(GEMSTOOL_Eofpc_IC)      :: Results_Eofpc
   TYPE(GEMSTOOL_LP_Eofpc_IC)   :: Results_LP_Eofpc
   TYPE(GEMSTOOL_LS_Eofpc_IC)   :: Results_LS_Eofpc
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          2STREAM ARGUMENTS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Version 2.4 code. Dimensioning uses that from VLIDORT
!    Optimized Code is single geometry, Upwelling, One Level

!  Rob Fix 8/3/18. MultiGeometry capability properly installed in 2S code.

!  Use separate names for the 2S Dimensioning

   integer, parameter :: S2_MAXLAYERS   = GT_maxlayers
   integer, parameter :: S2_MAXGEOMS    = GT_MaxGeometries   ! New 8/3/18
   integer, parameter :: S2_MAXTOTAL    = 2 * S2_MAXLAYERS
   integer, parameter :: S2_MAXMessages = GT_Maxmessages

!  Jacobian dimensions

   INTEGER, PARAMETER :: S2_max_atmoswfs   = GT_maxatmoswfs
   INTEGER, PARAMETER :: S2_max_surfacewfs = GT_maxsurfacewfs
   INTEGER, PARAMETER :: S2_max_sleavewfs  = GT_maxsurfacewfs    ! Currently no separate ariable for this.

!  New for  Version 2.1: Observational Geometry
!   INTEGER, PARAMETER :: max_user_obsgeoms_2s = max_geometries     !@@

!  Layering and geometry (New 8/3/18) indices

   integer              :: S2_NLAYERS, S2_NTOTAL, S2_NGEOMS

!  2S Flags: plane-parallel, deltam-2stream scaling

      LOGICAL           :: DO_PLANE_PARALLEL, DO_D2S_SCALING

!  Level output

      INTEGER           :: TOPMOST_LEVEL

!  2S flags: SURFACE and SLEAVE options (1/23/14, Version 2.3)

      LOGICAL           :: DO_BRDF_SURFACE
      LOGICAL           :: DO_SURFACE_LEAVING
      LOGICAL           :: DO_SL_ISOTROPIC

!  2Stream BVP control --- New 6/25/14, Version 2.3 and higher
!     2STREAM BVP INdex (0 = LAPACK,1 = PentaDiag) and inverse flag
!  * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!  BVP control --- New 6/25/14, Version 2.3 and higher
!  * BVP Scale Factor. Debug only. Set this to 1.0 on input

      logical, parameter    :: DO_PDINVERSE = .false.
      INTEGER, parameter    :: BVPINDEX     = 1
      REAL(GTPK)            :: BVPSCALEFACTOR

!  Linearization flags

      LOGICAL        :: DO_PROFILE_WFS
      LOGICAL        :: DO_SURFACE_WFS
      LOGICAL        :: DO_SLEAVE_WFS

!  Jacobian (linearization) control

      LOGICAL        :: LAYER_VARY_FLAG   ( S2_MAXLAYERS )
      INTEGER        :: LAYER_VARY_NUMBER ( S2_MAXLAYERS )
      INTEGER        :: N_SURFACE_WFS
      INTEGER        :: N_SLEAVE_WFS

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  2S Flags: Observational Geometry, Levelout , Flux options. Not in optimized Version
!      LOGICAL           :: DO_USER_OBSGEOMS,DO_2S_LEVELOUT,: DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT
!  2S flags: Sources control, including thermal. Not in optimized Version
!      LOGICAL           :: DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_SOLAR_SOURCES
!  Order of Taylor series (including terms up to EPS^n). Version 2.4. Not in optimized Version
!      INTEGER           :: TAYLOR_ORDER
!      REAL(gtpk)     :: TAYLOR_SMALL
!  Thermal Cutoff (actually a layer optical thickness minimum). Not in optimized Version
!     Rob, introduced 14 May 2015, Version 2.4, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers
!      REAL(gtpk)     :: TCUTOFF
!  2Stream Geometry. Lattice-geometry input. Angles MUST be in DEGREES
!      N_GEOMETRIES = NBEAMS * N_USER_ANGLES * N_USER_RELAZMS (Lattice value)
!      N_GEOMETRIES = N_USER_OBSGEOMS                         (ObsGeom value)
!      integer         :: nbeams, n_user_angles, n_user_relazms, n_geometries
!      REAL(gtpk)   :: BEAM_SZAS    ( MAXBEAMS )
!      REAL(gtpk)   :: USER_ANGLES  ( MAX_USER_STREAMS )
!      REAL(gtpk)   :: USER_RELAZMS ( MAX_USER_RELAZMS )
!   2Stream Geometry. Observational Geometry Input (New for  Version 2.1)
!     If set, this will override Lattice-geometry input
!      INTEGER       :: N_USER_OBSGEOMS                      !@@
!      REAL(gtpk) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS_2S,3) !@@
!  2Stream Thermal inputs. Not in optimized Version
!      REAL(gtpk)   :: SURFBB, THERMAL_BB_INPUT  ( 0:MAXLAYERS )
!  2Stream Emissivity. NOt required
!      REAL(gtpk)   :: EMISSIVITY
!  2Stream height and earth radius
!      REAL(gtpk)   :: EARTH_RADIUS
!      REAL(gtpk)   :: HEIGHT_GRID ( 0:MAXLAYERS )
!  Intensities, TOA/BOA
!      REAL(gtpk)   :: INTENSITY_TOA(MAX_GEOMETRIES)
!      REAL(gtpk)   :: INTENSITY_BOA(MAX_GEOMETRIES)
!  2stream Flux output
!     REAL(gtpk)    :: FLUXES_TOA(MAXBEAMS,2)
!     REAL(gtpk)    :: FLUXES_BOA(MAXBEAMS,2)
!  2stream output solutions at ALL levels
!      REAL(gtpk)    :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
!      REAL(gtpk)    :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2Stream Stream value

      REAL(gtpk)   :: STREAMVAL, STREAMINV

!  Geometry input   ! Geometry dimension, New 8/3/18
!     GEOMETRIES(v,1,1:3) = SZA angle, cosine, secant, sine
!     GEOMETRIES(v,2,1:3) = SZA angle, cosine, secant, sine
!     GEOMETRIES(v,3,1:3) = SZA angle, cosine, secant, sine

      REAL(gtpk)  :: GEOMETRIES(S2_MAXGEOMS,3,3)
      REAL(gtpk)  :: AUX_GEOMS(S2_MAXGEOMS,7)

!  Flux factor

      REAL(gtpk) :: FLUX_FACTOR

!  Chapman Factors  ! Geometry dimension, New 8/3/18

      real(gtpk)  :: CHAPFACS ( S2_MAXLAYERS, S2_MAXLAYERS, S2_MAXGEOMS )

!  2Stream Atmospheric Optical properties

      REAL(gtpk)   :: DELTAU_INPUT(S2_MAXLAYERS)
      REAL(gtpk)   :: OMEGA_INPUT (S2_MAXLAYERS)
      REAL(gtpk)   :: ASYMM_INPUT (S2_MAXLAYERS)
      REAL(gtpk)   :: D2S_SCALING (S2_MAXLAYERS)

!  2Stream Atmospheric Linearized optical properties

      REAL(gtpk)  :: L_DELTAU_INPUT(S2_MAXLAYERS, S2_MAX_ATMOSWFS)
      REAL(gtpk)  :: L_OMEGA_INPUT (S2_MAXLAYERS, S2_MAX_ATMOSWFS)
      REAL(gtpk)  :: L_ASYMM_INPUT (S2_MAXLAYERS, S2_MAX_ATMOSWFS)
      REAL(gtpk)  :: L_D2S_SCALING (S2_MAXLAYERS, S2_MAX_ATMOSWFS)

!  Albedo and BRDF fourier components ! Geometry dimension, New 8/3/18
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams    !  NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(gtpk)   :: ALBEDO
      REAL(gtpk)   :: BRDF_F_0  ( 0:1, S2_MAXGEOMS )
      REAL(gtpk)   :: BRDF_F    ( 0:1 )
!      REAL(gtpk)  :: UBRDF_F_0 ( 0:1, S2_MAXGEOMS )
      REAL(gtpk)   :: UBRDF_F   ( 0:1, S2_MAXGEOMS )

!  Linearized BRDF fourier components ! Geometry dimension, New 8/3/18
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident quadrature stream, reflected user streams

      REAL(gtpk) :: LS_BRDF_F_0  ( S2_MAX_SURFACEWFS, 0:1, S2_MAXGEOMS  )
      REAL(gtpk) :: LS_BRDF_F    ( S2_MAX_SURFACEWFS, 0:1)
      REAL(gtpk) :: LS_UBRDF_F   ( S2_MAX_SURFACEWFS, 0:1, S2_MAXGEOMS  )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  Do not require any first-order inputs (exact or Fourier)
!  Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms  ! Geometry dimension, New 8/3/18

      REAL(gtpk)   ::  SLTERM_ISOTROPIC ( S2_MAXGEOMS )
      REAL(gtpk)   ::  SLTERM_F_0 ( 0:1, S2_MAXGEOMS )

!  Linearized surface leaving ! Geometry dimension, New 8/3/18
!   @@@ Addition of SLEAVE WF inputs, R. Spurr, 23 January 2014 @@@@@@@@@

      REAL(gtpk)   :: LSSL_ISOTROPIC  ( S2_MAX_SLEAVEWFS, S2_MAXGEOMS )
      REAL(gtpk)   :: LSSL_F_0        ( S2_MAX_SLEAVEWFS, 0:1, S2_MAXGEOMS )

!  2Stream Results output.
!  ----------------------

!  Single-level intensity  ! Geometry dimension, New 8/3/18

      REAL(gtpk)   :: INTENSITY_2S ( S2_MAXGEOMS )

!  Jacobian output ! Geometry dimension, New 8/3/18

      REAL(gtpk)   :: PROFILEWFS_2S(S2_MAXGEOMS, S2_MAXLAYERS, S2_MAX_ATMOSWFS)
      REAL(gtpk)   :: SURFACEWFS_2S(S2_MAXGEOMS, S2_MAX_SURFACEWFS)

!  2Stream Exception handling

!    1. Check Messages and actions

      INTEGER          :: STATUS_INPUTCHECK
      INTEGER          :: C_NMessages
      CHARACTER*100    :: C_Messages(0:S2_MAXMessages)
      CHARACTER*100    :: C_ACTIONS (0:S2_MAXMessages)

!    2. Execution Message and 2 Traces

      INTEGER          :: STATUS_EXECUTION
      CHARACTER*100    :: S2_Message, S2_TRACE_1, S2_TRACE_2

!   2Stream  ARGUMENTS (Not I/O to the model)

!      LOGICAL, parameter :: DO_FULLQUADRATURE = .false.
      LOGICAL, parameter :: DO_FULLQUADRATURE = .true.

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          FIRST ORDER VARIABLES
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  FIFTH REVISION: Complete relabeling
!     --- Keep distinction of variables

!  Dimensions - Use Fine Layering, One output level

   integer, parameter :: FO_MaxGeometries  = GT_MaxGeometries        ! VLIDORT value, GEMSTOOL value
   integer, parameter :: FO_maxlayers      = GT_maxlayers            ! VLIDORT_value, GEMSTOOL value
   integer, parameter :: FO_maxmoments     = GT_MaxMoments           ! VLIDORT value, GEMSTOOL value
   integer, parameter :: FO_maxatmoswfs    = GT_maxatmoswfs          ! VLIDORT value, GEMSTOOL value
   integer, parameter :: FO_maxsurfacewfs  = GT_maxsurfacewfs        ! VLIDORT value, GEMSTOOL value

!  Now for the Optimized, Leveraged Code, this is a unique dimension for FO

   integer, parameter   :: FO_maxfine       = 15

!  Critical attenuation, set as parameter here

   real(gtpk), parameter :: FO_Acrit = 1.0e-10_gtpk

!  Constants (Copy VLIDORT values)

   real(gtpk) :: FO_Pie, FO_dtr

!  2. Control Variables
!  --------------------

!  Derived flags, optical settings

   logical    :: FO_do_regular_ps
   logical    :: FO_do_enhanced_ps
   logical    :: FO_do_planpar
   logical    :: FO_do_deltam
   logical    :: FO_do_lambertian
   logical    :: FO_do_obsgeoms
   logical    :: FO_do_Chapman
   logical    :: FO_Do_Sunlight
   logical    :: FO_Do_Sleave

!  Jacobian Flags. Renamed surface flag, added sleave flag (6/25/18)

   LOGICAL :: FO_do_reflecwfs
   LOGICAL :: FO_do_sleavewfs
!   LOGICAL :: FO_do_surfacewfs
   LOGICAL :: FO_do_profilewfs

!  Layer control. Finelayer input

   integer    :: FO_nlayers
   integer    :: FO_nfineinput

!  user level

   integer    :: FO_aclevel

!  Number of geometries, stokes

   integer    :: FO_ngeoms, FO_nszas, FO_nvzas, FO_nazms
   integer    :: FO_nstokes

!  Jacobian control. Distinguish relfec and sleave numbers. (6/25/18)
!    n_surfacewfs = SUM OF n_reflecwfs + n_sleavewfs. [Not declared here]

   LOGICAL :: FO_Lvaryflags(FO_maxlayers)
   INTEGER :: FO_Lvarynums (FO_maxlayers)
   LOGICAL :: FO_Lvarymoms (FO_maxlayers,FO_maxatmoswfs)

   INTEGER :: FO_n_reflecwfs
   INTEGER :: FO_n_sleavewfs
!   INTEGER, Intent(in) :: FO_n_surfacewfs

!  Radius + heights

   real(gtpk)  :: FO_eradius, FO_heights (0:FO_maxlayers)
   real(gtpk)  :: FO_diffgrid(FO_maxlayers)

!  input angles (Degrees), VSIGN = +1 (Up); -1(Down)

   real(gtpk)  :: FO_vsign
   real(gtpk)  :: FO_alpha_boa(FO_MaxGeometries), FO_theta_boa(FO_MaxGeometries), FO_phi_boa(FO_MaxGeometries)
   real(gtpk)  :: FO_obsgeom_boa(FO_MaxGeometries,3)

!  Critical adjustment for cloud layers

   logical     :: FO_doCrit
   real(gtpk)  :: FO_extincs(FO_maxlayers)

!  Flag for the Nadir case.

   logical     :: FO_doNadir(FO_MaxGeometries)
  
!  Alphas,  Cotangents, Radii, Ray constant.

   real(gtpk)  :: FO_radii    (0:FO_maxlayers)
   real(gtpk)  :: FO_Raycon   (FO_MaxGeometries)
   real(gtpk)  :: FO_alpha    (0:FO_maxlayers,FO_MaxGeometries)
   real(gtpk)  :: FO_cota     (0:FO_maxlayers,FO_MaxGeometries)

!  LOS Quadratures for Enhanced PS

   integer     :: FO_nfinedivs(FO_maxlayers,FO_MaxGeometries)
   real(gtpk)  :: FO_xfine    (FO_maxlayers,FO_maxfine,FO_MaxGeometries)
   real(gtpk)  :: FO_wfine    (FO_maxlayers,FO_maxfine,FO_MaxGeometries)
   real(gtpk)  :: FO_csqfine  (FO_maxlayers,FO_maxfine,FO_MaxGeometries)
   real(gtpk)  :: FO_cotfine  (FO_maxlayers,FO_maxfine,FO_MaxGeometries)

!  Fine layering output

   real(gtpk)  :: FO_alphafine (FO_maxlayers,FO_maxfine,FO_MaxGeometries)
   real(gtpk)  :: FO_radiifine (FO_maxlayers,FO_maxfine,FO_MaxGeometries)

!  Critical layer

   integer    :: FO_Ncrit(FO_MaxGeometries)
   real(gtpk) :: FO_RadCrit(FO_MaxGeometries), FO_CotCrit(FO_MaxGeometries)

!  solar paths, Intent out 

   integer    :: FO_ntraverse      (0:FO_maxlayers,FO_MaxGeometries)
   real(gtpk) :: FO_sunpaths       (0:FO_maxlayers,FO_maxlayers,FO_MaxGeometries)

   integer    :: FO_ntraverse_fine (FO_maxlayers,FO_maxfine,FO_MaxGeometries)
   real(gtpk) :: FO_sunpaths_fine  (FO_maxlayers,FO_maxlayers,FO_maxfine,FO_MaxGeometries)
   real(gtpk) :: FO_Chapfacs       (FO_maxlayers,FO_maxlayers,FO_MaxGeometries)

!  Cosine scattering angle and other cosines

   real(gtpk) :: FO_cosscat_up(FO_MaxGeometries) 
   real(gtpk) :: FO_Mu1       (FO_MaxGeometries)
   real(gtpk) :: FO_Mu0       (FO_MaxGeometries)

!  Spherical Functions: Inputs. [Starter flag may be re-set]

   logICAL   :: FO_STARTER
   intEGER   :: FO_NMOMENTS

!  Outputs from the spherical Functions routine

   REAL(gtpk) :: FO_ROTATIONS(4,FO_MaxGeometries)
   REAL(gtpk) :: FO_GENSPHER(0:FO_MaxMoments,4,FO_MaxGeometries)
   REAL(gtpk) :: FO_GSHELP(7,FO_MaxMoments)

!  4. Optical properties
!  ---------------------

!  optical inputs Atmosphere

   real(gtpk) :: FO_Extinction   ( FO_MAXLAYERS )
   real(gtpk) :: FO_DELTAUS      ( FO_MAXLAYERS )
   real(gtpk) :: FO_OMEGAS       ( FO_MAXLAYERS )
   real(gtpk) :: FO_TRUNCFAC     ( FO_MAXLAYERS )
   real(gtpk) :: FO_EXACTSCAT_UP ( FO_MAXLAYERS, 4, 4, FO_MaxGeometries )
   real(gtpk) :: FO_GREEKMAT     ( 0:FO_MaxMoments, FO_Maxlayers, 16 )

!  Solar Flux and Surface reflectivity (Could be the albedo)
!  surface-leaving term, 6/25/18 for AVIRIS.

   real(gtpk) :: FO_REFLEC(4,4,FO_MaxGeometries), FO_FLUX,  FO_FluxVec(4)
   real(gtpk) :: FO_SLTERM(4,FO_MaxGeometries)

!  Linearized optical inputs

   real(gtpk) :: FO_L_EXTINCTION   ( FO_MAXLAYERS, FO_maxatmoswfs )
   real(gtpk) :: FO_L_DELTAUS      ( FO_MAXLAYERS, FO_maxatmoswfs )
   real(gtpk) :: FO_L_OMEGAS       ( FO_MAXLAYERS, FO_maxatmoswfs )
   real(gtpk) :: FO_L_TRUNCFAC     ( FO_MAXLAYERS, FO_maxatmoswfs )
   real(gtpk) :: FO_L_exactscat_up ( FO_MAXLAYERS,4,4,FO_MaxGeometries,FO_maxatmoswfs)
   real(gtpk) :: FO_L_GREEKMAT     ( 0:FO_MaxMoments, FO_Maxlayers, 16, FO_maxatmoswfs )

   real(gtpk) :: FO_LS_REFLEC     ( 4,4,FO_MaxGeometries, FO_maxsurfacewfs)
   real(gtpk) :: FO_LSSL_SLTERM   ( 4,FO_MaxGeometries,   FO_maxsurfacewfs)

!  5. Output Variables and arrays
!  ------------------------------

!  First-order Stokes vectors and Direct-beam reflectance

   real(gtpk) :: FO_Stokes_up     ( 4, FO_MaxGeometries )
   real(gtpk) :: FO_Stokes_db     ( 4, FO_MaxGeometries )

!  First order output dummy

   real(gtpk) :: FO_cumsource_up  ( 0:FO_maxlayers, 4, FO_MaxGeometries)


!  First order Jacobian output

   real(gtpk) :: FO_LP_Jacobians_up  ( 4, FO_MaxGeometries, FO_maxlayers, FO_maxatmoswfs )
   real(gtpk) :: FO_LP_Jacobians_db  ( 4, FO_MaxGeometries, FO_maxlayers, FO_maxatmoswfs )
   real(gtpk) :: FO_LS_Jacobians_db  ( 4, FO_MaxGeometries, FO_maxsurfacewfs )

!  Exception handling on Geometry routine

   logical         :: FO_fail
   character*100   :: FO_message
   character*100   :: FO_trace

!  Local output variables
!  Eofpc_Stokes
   real (gtpk), dimension (:,:,:),allocatable :: Stokes_Eofpc, Stokes_Corrfacs
!  !Polarization
   real (gtpk), dimension (:,:), allocatable :: DOLP_Eofpc, DOCP_Eofpc
!  Profile Jacobians
   real (gtpk), dimension (:,:,:), allocatable :: LS_Jacobians_Eofpc
!  

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          PCAL2S TOOL LOCAL ARRAYS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Existence flags (avoids any numerical problems)

   LOGICAL, allocatable :: LP_Jacs_Exist ( :,: )
   LOGICAL, allocatable :: LS_Jacs_Exist ( : )

!  SVD CMP (parameter here)

   logical, parameter :: do_svd_cmp = .false.

!  Variables for creatGT_bin
!    NCNT and INDEX are always produced
!    BIN or BINLIMS not always produced (one or the other)

   integer :: PCA_neofs (0:GT_maxbins), PCA_nbins
   integer    :: PCA_ncnt(0:GT_maxbins),index(GT_MaxWav)
   real(gtpk) :: PCA_bins(0:GT_maxbins)

!  Variables for pca_performer

   integer :: ndim

!  PCA Performer Output (VLIDORT optical properties for each EOF)
!  -------------------------------------------------------------

!  PCA Projections. Type structure must now be Intent(InOut)

   TYPE(GEMSTOOL_PCAProj_Optical) :: PCAProj

!  Principal components. Local and Global arrays are allocatable

   real(gtpk), allocatable :: PrinComps(:,:)
   real(gtpk), allocatable :: Princomps_Local(:,:)

!  Correction factors --> Local allocatable arrays, Output to  PCA Correction routines. 
!   Q array added 5/20/16, U array 9/23/16. Merged QU, 12/28/17.

   real(gtpk), allocatable :: Intensity_Corrfacs_Local(:,:)
   real(gtpk), allocatable :: LP_IJacs_Corrfacs_Local(:,:,:,:)
   real(gtpk), allocatable :: LS_IJacs_Corrfacs_Local(:,:)

   real(gtpk), allocatable :: StokesQU_Corrfacs_Local(:,:,:)
   real(gtpk), allocatable :: LP_QUJacs_Corrfacs_Local(:,:,:,:,:)
   real(gtpk), allocatable :: LS_QUJacs_Corrfacs_Local(:,:,:)

!  Exception handling for PCA driver

   CHARACTER*100 :: PCA_MessAGE, PCA_trace_1, PCA_trace_2
   LOGICAL       :: PCA_FAIL

!  Intensity/StokesQU results for PCA-projections
!    - Global arrays, must now be Intent(InOut) to PCA routines. 06 May 2015
!    - StokesQU stuff added 5/20/16 (Q), 9/23/16 (U), Merged, 12/28/17
!    - Indices 1/2 = Q/U (First). 1 = LD, 2 = FO (Second)

   real(gtpk) :: intensity_Bin  (0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,3)  ! 1 = LD, 2 = 2S, 3 = FO
   !! real(gtpk) :: LP_IJacs_Bin   (GT_maxatmoswfs,GT_maxlayers,0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,3)
   real(gtpk), allocatable :: LP_IJacs_Bin(:,:,:,:,:,:) !   (GT_maxatmoswfs,GT_maxlayers,0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,3)
   real(gtpk) :: LS_IJacs_Bin   (0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,3)

   real(gtpk) :: StokesQU_Bin   (2,0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,2)
   real(gtpk) :: LP_QUJacs_Bin  (2,GT_maxatmoswfs,GT_maxlayers,0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,2)
   real(gtpk) :: LS_QUJacs_Bin  (2,0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,2)

!  Intensity/StokesQU results for PCA-projections
!     - Local arrays, Inputs to PCA Correction routines.06 May 2015. Made allocatable, 12/28/17

   real(gtpk), allocatable :: intensity_Bin_Local(:,:,:)
   real(gtpk), allocatable :: StokesQU_Bin_Local(:,:,:,:)
   real(gtpk), allocatable :: LP_IJacs_Bin_Local(:,:,:,:,:)
   real(gtpk), allocatable :: LS_IJacs_Bin_Local(:,:,:)
   real(gtpk), allocatable :: LP_QUJacs_Bin_Local(:,:,:,:,:,:)
   real(gtpk), allocatable :: LS_QUJacs_Bin_Local(:,:,:,:)

!  Fast Intensity and stokes results

   real(gtpk) :: Intensity_2S_Fast(GT_maxGeometries)
   real(gtpk) :: LP_Intensity_2S_Fast(GT_maxatmoswfs,GT_maxlayers,GT_maxGeometries)
   real(gtpk) :: LS_Intensity_2S_Fast(GT_maxSurfacewfs,GT_maxGeometries)

   real(gtpk) :: Stokes_FO_Fast(GT_maxGeometries,4)
   real(gtpk) :: LP_Stokes_FO_Fast(GT_maxatmoswfs,GT_maxlayers,GT_maxGeometries,4)
   real(gtpk) :: LS_Stokes_FO_Fast(GT_maxSurfacewfs,GT_maxGeometries,4)

!  Proxy GEMSTOOL variables
!  ========================

!  Aerosols/Clouds

   logical :: do_Aerosols, do_Clouds

!  Enhanced sphericity flag, Optional 2S flag.  Controlled by the Config-file input "do_fast_calculation"

   logical :: do_enhanced_ps
   logical :: do_optional_2stream

!    --> Number of LIDORT discrete ordinates, number of stokes components (3=Vector, 1 = Scalar)
!    --> Number of Geometries

   integer :: ngeoms, nstreams, nstokes, nlayers, nlevels, nwav, nstr2, ndat, nmuller

!    --> PCA inputs

   integer :: nlayers1, nlayers2, nlayers21, nlayers22, nspars, npars
   integer :: PCA_Strategy_index, PCA_Binning_index
   logical :: do_fast_calculation, alb_pcainclude, do_3M_correction
   REAL (KIND=8) :: amf
!  Serial/Timing Variables
!  =======================

!  Serial tests (timing)

   REAL  :: ser_e1, ser_e2
   real  :: e1,e2,e3 
   real  :: SereofpcTime

!  Other Variables
!  ===============

!  timing variables

   real :: Eofpctimes(20), EofpcRunTime

!  OD help
   INTEGER, PARAMETER :: maxscatter=3, maxgkmatc=8, maxgksec=6
   logical, DIMENSION(GT_maxlayers)  :: Aerflag, CldFlag
   real(kind=GTPK) :: Aer_WT, Ray_WT, CLD_WT, L_Aer_WT, L_Ray_WT, L_Cld_wt, &
                      LD_taudp, LD_omega, fac, beta2, dnm1, L_fact1
   real(kind=GTPK) :: SK, SQ2U2, SSRad(4), SV2, DEPOL, Int, jac, mom, momv, fact1, extinc, six, mrtsix
   real(kind=GTPK) :: depol_prj, MOMA, MOMR, MOM2M, RAYWT, AERWT, sumtaug
   real(kind=GTPK) :: A1Cofs(0:GT_maxmoments), A2Cofs(0:GT_maxmoments), RYCofs(0:2)
   real(kind=GTPK) :: MSRad(4), MSJac(4), FORad(4), FOJac(4), Rad(4), Rad1, Ifast(3), Jfast(3), ISurf(4)
   Real(kind=GTPK), DIMENSION(0:GT_maxmoments, 1:maxgksec, maxscatter) :: phasmoms_input
   real(kind=GTPK), DIMENSION(0:GT_maxmoments, 1:maxgksec) :: phasmoms_total_input
   real(kind=GTPK), DIMENSION(0:GT_maxmoments, 1:maxgksec) :: l_phasmoms_total_input
   real(kind=GTPK), DIMENSION(gt_maxatmoswfs, 0:GT_maxmoments, 1:MAXSTOKES_SQ) :: l_greekmat_total_input
   real(kind=GTPK), ALLOCATABLE, DIMENSION(:) :: taugcum, tauhgcum, taucum
  
!  Help variables

   logical        :: Fail1, do_Continuum, Use_Hitran, DO_Vlidort_Inpdebug, LFVary(GT_maxlayers)
   integer        :: s,q,l,n,nf,ns,w,m,v,k,kk,k1,i,mm,ii
   integer        :: LM, local_nmoms, nscatter, ngkmatc, ncoeffs , n_totalatmos_wfs !, n_surface_wfs
   integer, DIMENSION (GT_maxlayers) ::  LNVary
   INTEGER        :: VK, GK, RK, CK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)
   character*5    :: c5
   character*3    :: c3
   character*4    :: c4
   character*100  :: Message1

!  Linearized optical projections

   real(GTPK) :: L_OPD_PRJ(GT_MaxAtmoswfs,GT_Maxlayers,GT_Maxbins,GT_MaxEofs2p1)
   real(GTPK) :: L_SSA_PRJ(GT_MaxAtmoswfs,GT_Maxlayers,GT_Maxbins,GT_MaxEofs2p1)
   real(GTPK) :: L_FR_PRJ (GT_MaxAtmoswfs,GT_Maxlayers,GT_Maxbins,GT_MaxEofs2p1)
   real(GTPK) :: L_FA_PRJ (GT_MaxAtmoswfs,GT_Maxlayers,GT_Maxbins,GT_MaxEofs2p1)

!  PCA bookkeeping help

   integer       :: istart,istart_save(0:GT_maxbins-1),iend,indexl,ndiff_neofp1,PCA_neof
   integer       :: irt_prj,irt_map(GT_maxbins*GT_maxeofs2p1,2)

!  Debug output. FORT 8811.

   logical, parameter :: do_debug_input    = .false.
   logical, parameter :: do_debug_output   = .false.
   logical, parameter :: do_Sun_Normalized = .true.
   logical, save :: first = .true.
!  Mask Data
!  =========

   GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)     ! Greek Matrix indices
   RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)         ! Convention for Rayleigh

!   CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)         ! Convention for TEMPO aerosol output
!   SMASK = (/  1, 1, 1, 1, 1, -1, 1, 1 /)        ! Sign Mask  for TEMPO aerosol output

   CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)        ! Convention for RTS MIE aerosol output
   SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)     ! Sign Mask  for RTS MIE aerosol output

!  START CODE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   !write(*,*) '    ---- Doing PCA Eofpc-plus RT  (with Jacobians)'

!  Initialize time

!   eofpcTimes   = 0.0 ! pre-initialized now, 2/18/16
!   eofpcRunTime = 0.0 ! pre-initialized now, 2/18/16

    SereofpcTime = 0.0

!  Initialize output - no longer needed, pre-initialized now, 2/18/16
!   Intensity_eofpc    = zero
!   Intensity_Corrfacs = zero

!  Initialize Exception handling

   fail = .false.
   Nmessages = 0
   messages  = ' '

!  Initial CPU call

   if ( Monitor_CPU ) call cpu_time(e1)

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  1. INITIAL SECTION. RTM SETUP AND GEOMETRY
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  start CPU time

   if ( Monitor_CPU ) call cpu_time(e2)

!  1a. VLIDORT control Read input, abort if failed
!  ===============================================

!  Proxies from Inputs
   nstokes  = Inputs%RTMcontrol%NVlidort_nstokes  ; ns = nstokes
   nstreams = Inputs%RTMcontrol%NVlidort_nstreams ; nstr2 = 2 * nstreams ; dnm1 = real(2*nstr2+1,gtpk)
   ngeoms   = Inputs%Geometry%N_GEMS_geometries

   do_aerosols = Inputs%Atmosph%do_aerosols
   do_clouds   = Inputs%Atmosph%do_clouds
   do_optional_2stream = Inputs%PCAcontrol%do_fast_Calculation
!  Proxies from Geophys

   nwav    = Geophys%WavGrids%nwav ; ndat = nwav
   nlayers = Geophys%Atmos%nlayers ; nlevels = nlayers + 1
   local_nmoms = Geophys%Atmos%ngreek_moments_input
!  Layers

   NLAYERS1 = NLAYERS + 1
   NLAYERS2  = 2 * nlayers
   NLAYERS21 = NLAYERS2 + 1
   NLAYERS22 = NLAYERS2 + 2

!  Fixtures

   use_hitran   = .false.
   do_Continuum = .false.

!  constants

   six    = 6.0_gtpk
   mrtsix = -sqrt(six)

! Optical Depth Properties
   n_totalatmos_wfs = RI%n_totalatmos_wfs
   n_surface_wfs =    RI%n_surface_wfs
! allocating global variables
  IF (allocated(Assigned_bins)) deallocate(Assigned_bins)
  allocate (Assigned_bins(ndat))
  allocate ( Stokes_Eofpc (ngeoms, ns, nwav))

! allocating local variables
  allocate ( Stokes_Corrfacs (ngeoms, ns, nwav))
  IF (ns > 1 ) THEN 
    allocate ( DOLP_Eofpc (ngeoms, nwav), DOCP_Eofpc(ngeoms, nwav))
  ENDIF
  IF (do_Jacobians) THEN 
  !  allocate (LP_Jacobians_Eofpc(nlayers, n_totalatmos_wfs, ns, nwav))
  !  allocate (LP_Jacobians_Corrfacs(nlayers, n_totalatmos_wfs,ns, nwav))
  ENDIF
  !print * , nlayers, n_totalatmos_wfs, ngeoms, ns, nwav
   LFvary (1:nlayers) =layer_vary_flag_cc(1:nlayers)
   LNvary (1:nlayers) =layer_vary_number_cc(1:nlayers)
   !Geophys%totalods%taudp(1:nlayers, 1:nwav)=RI%taudp(1:nlayers, 1:nwav)
   !Geophys%totalods%taug(1:nlayers, 1:nwav)=RI%taug(1:nlayers, 1:nwav)
   !Geophys%totalods%omega(1:nlayers,1:nwav)=RI%omega(1:nlayers, 1:nwav)
   !Geophys%TotalODs%fr(1:nlayers, 1:nwav) = RI%fr(1:nlayers,1:nwav)
   !Geophys%TotalODs%fa(1:nlayers, 1:nwav) = RI%fa(1:nlayers,1:nwav)
   !L_Geophys%L_TotalODs%l_taudp(1:n_totalatmos_wfs, 1:nlayers,1:nwav) = &
   !                   RI%l_taudp(1:n_totalatmos_wfs,1:nlayers,1:nwav)
   !L_Geophys%L_TotalODs%l_omega(1:n_totalatmos_wfs, 1:nlayers,1:nwav) = &
   !                   RI%l_omega(1:n_totalatmos_wfs,1:nlayers,1:nwav)

!  PCA Proxies
   PCA_Strategy_index = Inputs%PCAControl%PCA_Strategy_index 
   PCA_Binning_index  = Inputs%PCAControl%PCA_Binning_index
   IF (PCA_Binning_index /= 4) THEN 
     PCA_nbins          = Inputs%PCAControl%PCA_nbins
     PCA_neofs = 0 ; PCA_neofs (0:PCA_nbins-1) = Inputs%PCAControl%PCA_neofs (0:PCA_nbins-1)
   ENDIF
   do_fast_calculation = .true.
   alb_pcainclude      = Inputs%PCAControl%alb_pcainclude
   do_3M_correction    = Inputs%PCAControl%do_3M_correction
!  VLIDORT. Initialize the Input structures.
   !IF (num_iter == 0) THEN
   !IF ( first ) THEN 
        call GEMSTOOL_VLIDORT_Initialize &
       ( Inputs, nlayers, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, fail, Message1 )
   !ENDIF

   VLIDORT_FixIn%Cont%TS_NLAYERS          = nlayers
   IF ( VLIDORT_FixIn%Bool%TS_do_dnwelling ) then
      VLIDORT_ModIn%MUserVal%TS_user_levels(1) = Real(nlayers,fpk)
   ENDIF

     VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = 1
     VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1)   = Inputs%Geometry%GEMS_azms(1)
     VLIDORT_ModIn%MSunrays%TS_N_SZANGLES       = 1
     VLIDORT_ModIn%MSunrays%TS_SZANGLES(1)      = Inputs%Geometry%GEMS_szas(1)
     VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES  = 1
     VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1)= Inputs%Geometry%GEMS_vzas(1)
     VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS  = 1
     VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1)   = Inputs%Geometry%GEMS_szas(1)
     VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2)   = Inputs%Geometry%GEMS_vzas(1)
     VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3)   = Inputs%Geometry%GEMS_azms(1)
     VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT    = local_nmoms
   !ENDIF
!  Fix heights

   VLIDORT_FixIn%Chapman%TS_height_grid(0:nlayers) = Geophys%Atmos%Level_heights(0:nlayers)
   VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT   = Geophys%Atmos%Level_heights(nlayers)
!  linearization control
!mick fix 7/30/2018 - initialize DO_ATMOS_LBBF, DO_SURFACE_LBBF, COLUMNWF_NAMES, & PROFILEWF_NAMES 

  ! IF (num_iter == 0) THEN 
   VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION      = .false.
   VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION       = .false.
   VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION      = .false.
   VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY            = .true.
   VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS   = .false.
   VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF   = .false.
   VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF = .false.

   VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS       = 0
   VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS   = 0
   VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS  = 0
   VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS        = 0
   VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES      = ''
   VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES     = ''

!  Initial setting, just one profile Jacobian
!   Rob 8/3/18. Add Surface albedo Jacobian control

   if ( do_Jacobians ) then
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                 = n_surface_Wfs
      IF (n_surface_wfs > 0) THEN 
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION     = .true.
      ENDIF
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION     = .true.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION       = .true.
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION             = .true.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:nlayers)    =  LFVary(1:nlayers)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:nlayers)  =  LNVary(1:nlayers)
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY           = .false.
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS            = n_totalatmos_wfs
   endif
  !ENDIF
!  Optical property zeroing
!  ------------------------

!  Entries in Greekmat

   if ( nstokes .eq. 1 ) then
      ngkmatc = 1 ; nmuller = 1
   else if ( nstokes .eq. 3 ) then
      ngkmatc = 5 ; nmuller = 6
   else if ( nstokes .eq. 4 ) then
      ngkmatc = 8 ; nmuller = 6
   endif

!  Preliminary zeroing
   
   VLIDORT_FixIn%Optical%TS_Greekmat_total_input(:,:,:) = GTZero
   VLIDORT_FixIn%Optical%TS_deltau_vert_input (:)       = GTZero
   VLIDORT_ModIn%MOptical%TS_omega_total_input(:)       = GTZero

   VLIDORT_LinFixIn%Optical%TS_L_Greekmat_total_input = GTZero
   VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = GTZero
   VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = GTZero

   ASYMM_INPUT    = GTZero ; D2S_SCALING   = GTZero
   deltau_input   = GTZero ; omega_input   = GTZero
   L_ASYMM_INPUT  = GTZero ; L_D2S_SCALING = GTZero
   L_deltau_input = GTZero ; L_omega_input = GTZero

!  Flags
!mick fix 7/30/2018 - added IF conditions

   !do n = 1, nlayers
   !  AERFLAG(n) = DO_AEROSOLS  .AND. Geophys%Aerosols%AEROSOL_LAYERFLAGS(N)
   !  CLDFLAG(n) = DO_CLOUDS    .AND. Geophys%Clouds%CLOUD_LAYERFLAGS(N)
   !enddo
   if ( do_aerosols ) then
     AERFLAG(1:nlayers) = Geophys%Aerosols%AEROSOL_LAYERFLAGS(1:nlayers)
   else
     AERFLAG(1:nlayers) = .false.
   endif
   if ( do_clouds ) then
     CLDFLAG(1:nlayers) = Geophys%Clouds%CLOUD_LAYERFLAGS(1:nlayers)
   else
     CLDFLAG(1:nlayers) = .false.
   endif

!  Exception handling
   if ( fail ) then
      NMessages = 1 ; fail = .true.
      Messages(1) = Trim(Message1) ; return
   endif

   VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT    = local_nmoms
!  Fix heights

!  1b. FO Settings
!  ===============

!  Flux vector

   FO_do_Sunlight  = .true.
   FO_FluxVec(1)   = GTone
   FO_FluxVec(2:4) = GTZero

!  Nmoments comes from the largest Geophys values

   FO_nmoments = 2
   if ( do_aerosols ) FO_nmoments = MaxVal(Geophys%Aerosols%aerosol_nscatmoms(1:nwav))
   if ( do_clouds )   FO_nmoments = Max(FO_nmoments,Geophys%Clouds%cloud_nscatmoms)

!  Set Local sphericity for FO. Control comes from GEMSTOOLinputs

   FO_do_enhanced_ps = .not.Inputs%RTMControl%do_regular_ps
   FO_do_regular_ps  = Inputs%RTMControl%do_regular_ps
   FO_do_planpar     = .false.
   FO_do_Chapman     = .true.

!  Number of Stokes components

   FO_Nstokes  = Inputs%RTMControl%NVlidort_nstokes

!  constants (copy VLIDORT values)

   FO_Dtr = DEG_TO_RAD
   FO_PIE = PIE

!  Fine layer control

   FO_nfineinput = 6

!  General flags

   FO_starter           = .true.   ! Should always be true  before first Spherfuncs call
   FO_do_lambertian     = .not.Inputs%BRDF%do_Gemstool_BRDF
   FO_do_deltam         = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING

!  Toa output

   FO_aclevel  = 0

!  surface leaving

   FO_do_sleave = .false.
   FO_slterm         = GTZero
   FO_LSSL_slterm    = GTZero

!  Linearization Control. Initial setting 7/9/18.
!mick fix 7/30/2018 - initialize "FO_n_reflecwfs" & "FO_n_sleavewfs"

   FO_do_reflecwfs  = .false.
   FO_do_sleavewfs  = .false.
   FO_do_profilewfs = VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION

   FO_Lvaryflags(1:nlayers) = LFVary(1:nlayers)
   FO_LVarynums (1:nlayers) = LNVary(1:nlayers)
   FO_Lvarymoms = .false.                          ! Initial setting.

   FO_n_reflecwfs = 0
   FO_n_sleavewfs = 0

!  Rob Fix 8/3/18. Add surface Jacobian control

   if ( Inputs%LinControl%do_Surface_Jacobians ) then
      FO_do_reflecwfs  = .true.
      FO_n_reflecwfs   = 1
   endif

!  Local Geomemtry. Single valued. Now Multi-valued 8/3/18.

   FO_do_obsgeoms = .true.
   FO_nszas = ngeoms ; FO_nvzas = ngeoms ; FO_nazms = ngeoms ; FO_ngeoms  = ngeoms
   do n = 1, FO_ngeoms
      FO_obsgeom_boa(n,1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(N,1:3)
      FO_theta_boa(n) = FO_obsgeom_boa(n,1)
      FO_alpha_boa(n) = FO_obsgeom_boa(n,2)
      FO_phi_boa(n)   = FO_obsgeom_boa(n,3)
   enddo
   FO_vsign = GTOne   ! scattering angle control, Upwelling only

!mick fix 2/26/2016 - initialize FO_doNadir flag

   FO_doNadir = .false.

!  Attenuation control

   FO_doCrit          = .false.        ! Criticality not set (fast option)

!   if ( do_cloud2.and.FO_do_enhanced_ps )  FO_doCrit = .true.

!  Copy the earth radius, and height grid

   FO_nlayers = nlayers
   FO_eradius = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS
   FO_heights(0:nlayers) = Geophys%Atmos%Level_heights(0:nlayers)

!  Find the maximum extinction for criticality(Clouds as layers)
!    - Corrected Code 4/22/17, Takes proper care with Delta-M scaling ---------------NOT USED HERE.
!    - Extinc here is only for criticality, so will be a worst-case value.

   FO_extincs = GTZero
   do n = 1, FO_nlayers
      FO_diffgrid(n)   = FO_heights(n-1) - FO_heights(n)
      do w = 1, nwav
         extinc = Geophys%TotalODs%taudp(n,w)/FO_diffgrid(n)
         FO_extincs(n) = max(FO_extincs(n),extinc)
      enddo
   enddo

!  Initialize

   FO_Greekmat   = GTZero
   FO_L_Greekmat = GTZero
   FO_L_Truncfac = GTZero

!  1c. 2stream settings
!  ====================

   BVPSCALEFACTOR = GTOne

!  Set ntotal and nlayers, and ngeoms (new,8/3/18)

   S2_NLAYERS = nlayers
   S2_NTOTAL  = 2 * nlayers
   S2_NGEOMS  = ngeoms
   IF (ngeoms /= 1) THEN 
      WRITE(*,*) 'check if ngeoms == 1 or not' !!! ; stop
   ENDIF
!  Stream Value

   if ( DO_FULLQUADRATURE ) then
      STREAMVAL = GTOne / sqrt(3.0_gtpk)
   else
      STREAMVAL = 0.5_gtpk
   endif
   STREAMINV = GTOne / STREAMVAL

!  output level is TOA

   TOPMOST_LEVEL = 0

!  Linearization control (Remainder comes from input)

   DO_PROFILE_WFS = VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION
   DO_SURFACE_WFS = VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION
   DO_SLEAVE_WFS  =  VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS
   N_SURFACE_WFS  = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
   N_SLEAVE_WFS   = VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS

   LAYER_VARY_FLAG(1:NLAYERS)   = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:nlayers)
   LAYER_VARY_NUMBER(1:NLAYERS) = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:nlayers)

!  Geometry settings 
!  Originally, single geometry for optimized 2S code
!  Rob Fix 8/3/18. Multiple geometry for 2stream introduced.

   do v = 1, ngeoms
      GEOMETRIES(v,1,1)  = Inputs%Geometry%GEMS_szas(v)
      GEOMETRIES(v,2,1)  = Inputs%Geometry%GEMS_vzas(v)
      GEOMETRIES(v,3,1)  = Inputs%Geometry%GEMS_azms(v)
      DO K = 1, 3
         GEOMETRIES(v,K,2) = COS ( DEG_TO_RAD * GEOMETRIES(v,K,1) ) !  Cosines
         GEOMETRIES(v,K,3) = GTOne / GEOMETRIES(v,K,2)              !  Secants
      ENDDO
   enddo

!  Surface Flags

   DO_BRDF_SURFACE     = Inputs%BRDF%do_Gemstool_BRDF
   DO_SURFACE_LEAVING  = Inputs%Sleave%do_WaterLeaving
   DO_SL_ISOTROPIC     = .true.

!  Deltam-scaling flag MUST BE SET, regardless of the VLIDORT setting
!    Otherwise, results are unphysical. 2/5/15

   DO_D2S_SCALING      = .true.
   DO_PLANE_PARALLEL   = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL

!  Earth radius and heights (Copy FO values)
!   HEIGHT_GRID(0:nlayers)  = FO_heights(0:nlayers) 
!   EARTH_RADIUS            = FO_eradius

!  5/6/15 upgrade - initialize Surface BRDF inputs

   BRDF_F_0 = GTZero ; LS_BRDF_F_0 = GTZero
   BRDF_F   = GTZero ; LS_BRDF_F   = GTZero
   UBRDF_F  = GTZero ; LS_UBRDF_F  = GTZero

!  SLEAVE Terms

   SLTERM_ISOTROPIC    = GTZero ; LSSL_ISOTROPIC = GTZero
   SLTERM_F_0          = GTZero ; LSSL_F_0       = GTZero

!  Other linearization settings for PCA.

   npars  = 0 ; nspars = 0
   if ( do_Profile_Wfs ) then
     do n = 1, nlayers
       if ( LFVary(n) ) npars  = max(npars,LNvary(n))
     enddo
     allocate(LP_Jacs_Exist(npars,nlayers)) ; LP_Jacs_Exist(1:npars,1:nlayers) = .true.
   endif
   if ( do_Surface_Wfs ) then
      nspars = 1 ; allocate(LS_Jacs_Exist(nspars)) ; LS_Jacs_Exist(1) = .true.
   endif

!  Monitor for the above settings for VLIDORT and FO and S

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Eofpctimes(1) = e3 - e2  ! setuptime VL/FO/2S initialize.
   endif

!  1d. First-Order Preliminary Geometry calculations
!  =================================================

!  Monitoring

   if ( Monitor_CPU ) call CPU_time(e2)

!  Call to geometry master

   call FO_SSGeometry_Master &
      ( FO_MaxGeometries, FO_MaxGeometries, FO_MaxGeometries, FO_MaxGeometries, FO_maxlayers, FO_maxfine, & ! Input Dimensions
        FO_do_obsgeoms, FO_do_Chapman, FO_do_planpar, FO_do_enhanced_ps,                              & ! Input flags
        FO_ngeoms, FO_nszas, FO_nvzas, FO_nazms, FO_nlayers, FO_nfineinput, FO_dtr, FO_Pie, FO_vsign, & ! Input control and constants
        FO_eradius, FO_heights, FO_obsgeom_boa, FO_alpha_boa, FO_theta_boa, FO_phi_boa,               & ! Input geometry/heights
        FO_doNadir, FO_doCrit, FO_Acrit, FO_extincs, FO_Raycon, FO_radii, FO_alpha, FO_cota,          & ! Input/Output(level)
        FO_nfinedivs, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine, FO_alphafine, FO_radiifine,         & ! Output(Fine)
        FO_NCrit, FO_RadCrit, FO_CotCrit, FO_Mu0, FO_Mu1, FO_cosscat_up, FO_chapfacs,                 & ! Output(Crit/scat)
        FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,                               & ! Output(Sunpaths)
        FO_fail, FO_message, FO_trace )                                                                 ! Output(Status)

!  Exception handling. Fatal error at this stage = return
               
   if ( FO_fail ) then
      messages(nmessages+1) = Adjustl(Trim(FO_message))
      messages(nmessages+1) = Adjustl(Trim(FO_trace))
      nmessages = nmessages+2
      fail = .true. ; return
   endif

!  1e. First-Order Spherical Function calculations
!  ===============================================

!  Call optimized spherical functions routine

   Call FO_VectorSS_spherfuncs_Optimized &
        ( FO_STARTER, FO_Maxmoments, FO_MaxGeometries, FO_Do_Sunlight, & ! Inputs
          FO_nmoments, FO_ngeoms, FO_nstokes, FO_dtr, FO_vsign,        & ! Inputs
          FO_theta_boa, FO_alpha_boa, FO_phi_boa, FO_cosscat_up,       & ! Inputs
          FO_ROTATIONS, FO_GSHELP, FO_GENSPHER )                         ! Outputs &

!  1f. 2S preliminary geometry calculations
!  ========================================

!  1/30/18. New section. Get necessary Chapman factors for 2stream
!    call TWOSTREAM_BEAM_GEOMETRY_PREPARE &
!         ( S2_MAXLAYERS, NLAYERS, DO_PLANE_PARALLEL,    & ! Input
!           GEOMETRIES(1,1), EARTH_RADIUS, HEIGHT_GRID,  & ! Input
!           CHAPFACS, SZA_LEVEL_OUTPUT )               ! In/Out

!  Copy FO chapman factors. Reverse layering in 2S
!  Rob Fix 8/3/18. Multiple geometry for 2stream introduced.

   do v = 1, ngeoms
      do n = 1, nlayers
         CHAPFACS(1:nlayers,n,v) = FO_Chapfacs(n,1:nlayers,v)
      enddo
   enddo

!  Auxiliary Geometry
!  Rob Fix 8/3/18. Multiple geometry for 2stream introduced.

!   MU0_2S         = GEOMETRIES(1,2)
!   USER_STREAM_2S = GEOMETRIES(2,2)
!   CALL TWOSTREAM_AUXGEOM_PREPARE  ( GEOMETRIES(1,2), GEOMETRIES(2,2), STREAMVAL, AUX_GEOMS )
   CALL TWOSTREAM_AUXGEOM_PREPARE &
         ( S2_MAXGEOMS, S2_NGEOMS, GEOMETRIES(1:S2_MAXGEOMS,1,2), GEOMETRIES(1:S2_MAXGEOMS,2,2), STREAMVAL, AUX_GEOMS )

!  Monitor for the above settings for FO Geometry and Spherical functions, and 2S geometry

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Eofpctimes(2) = e3 - e2  ! setuptime FO Geometry/SpherFuncs, 2Stream geometry
   endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  2.     BIN CREATION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Create bins

   if ( Monitor_CPU ) call cpu_time(e2)

!  =================== UV Binning ===============================
!   Binning Index 0 : 1 or 2 bins only, gas absorption only  (325-335 nm option). Continuum option. 8/6/13
!   Binning Index 1 : 11 bins 0-->10, gas absorption only
!   Binning Index 2 : 22 bins 0-->21, gas absorption + halfvalue total scat
!    27 October 2011, Versions 1 and 2
!    06 August  2013, Versions 0. Use Taudp for the Continuum case
!  =================== UV Binning ================================

!  =================== TIR Binning ===============================
!   Binning Index 2, May/June 2015   : 6-10 bins, incremental and dynamic, with bisection
!   Binning Index 3, 11 January 2016 : 9 bin input, reduction and merging of sub-size bins
!  =================== TIR Binning ===============================

!  =================== XOVER Binning =============================
!   Binning Index 2, May/June 2015   : 6-10 bins, incremental and dynamic, with bisection
!   Binning Index 3, 11 January 2016 : 9 bin input, reduction and merging of sub-size bins
!  =================== XOVER Binning =============================

!  =================== JBAK Binning July 2018 =============================
!   Binning Index 4, Wavelength groupings
!  =================== JBAK Binning July 2018 =============================

   allocate (taugcum(nwav), tauhgcum(nwav), taucum(nwav))
   taugcum(1:nwav) = GTzero
   DO n = 1, nwav
      taugcum(n) = sum( Geophys%TotalODs%taug(1:nlayers, n) )
      tauhgcum(n) = sum( Geophys%TotalODs%tauhg(1:nlayers, n) )
      taucum(n) = sum( Geophys%TotalODs%taudp(1:nlayers, n) )
   ENDDO
   which_win = 1 ! UV
   IF (Geophys%WavGrids%wav(1) >= 500) THEN 
      which_win = 2 !VIS
   ENDIF
!   do_debug_bin = .false.
!   if ( PCA_Binning_index .eq. 0 ) then
!      if ( do_Continuum ) then
!         call GEMSTOOL_CreateBins_V0 &
!             ( ndat, nlayers, PCA_nbins, bins, Geophys%TotalODs%taudp, &
!               PCA_ncnt,index )
!      else
!         call GEMSTOOL_CreateBins_V0 &
!             ( ndat, nlayers, PCA_nbins, bins, Geophys%TotalODs%taug, &
!               PCA_ncnt,index )
!      endif
!   else if ( PCA_Binning_index .eq. 1 ) then
!      call GEMSTOOL_CreateBins_V1 &
!             ( ndat, nlayers, PCA_nbins, bins, Geophys%TotalODs%taug, Geophys%TotalODs%omega, &
!               PCA_ncnt,index )
!   else if ( PCA_Binning_index .eq. 2 ) then
!      call GEMSTOOL_CreateBins_V2 &
!             ( ndat, nlayers, PCA_nbins, Geophys%TotalODs%taug, Geophys%TotalODs%taudp, Geophys%TotalODs%omega, &
!               PCA_ncnt,index,Assigned_bins )
   if ( PCA_Binning_index .eq. 3 .and. which_win == 1 ) then
   !   CALL GEMSTOOL_CreateBins_V3 &
   !          ( ndat, nlayers, PCA_nbins, PCA_bins, Geophys%TotalODs%taug, &
   !            PCA_ncnt,index,Assigned_bins )
   !if ( PCA_Binning_index .eq. 4 ) then
   else if ( PCA_Binning_index .eq. 4 .and. which_win == 1) then
   !  CALL GEMSTOOL_CreateBins_V4 &
   !       (which_win, ndat, Geophys%WavGrids%wav(1:ndat), &
   !       index(1:ndat),Assigned_bins(1:ndat),&
   !       PCA_nbins, PCA_ncnt,PCA_neofs, PCA_bins)
   ELSE IF (PCA_Binning_index .eq. 5) THEN
         amf = minval([cos(Inputs%Geometry%GEMS_szas(1)*3.14/180.),cos(Inputs%Geometry%GEMS_vzas(1)*3.14/180.)])
         amf = maxval([(Inputs%Geometry%GEMS_szas(1)),(Inputs%Geometry%GEMS_vzas(1))])
          !print * , taugcum(1:ndat)
    
     CALL GEMSTOOL_CreateBins_V5 (which_win, ndat, tauhgcum(1:ndat),taugcum(1:ndat), &
          Geophys%WavGrids%Wav(1:ndat), amf, &
          index(1:ndat), Assigned_bins(1:ndat), &
          PCA_nbins, PCA_ncnt, PCA_neofs, PCA_bins)
     
   ELSE IF (PCA_Binning_index .eq. 6) THEN !Add new PCA binning option from Sunny
      amf = maxval([Inputs%Geometry%GEMS_szas(1), Inputs%Geometry%GEMS_vzas(1)])
      bintest=.false.
      
      !CALL GEMSTOOL_CreateBins_jbak  (ndat,&
      !   habsline(1:ndat), absline(1:ndat),Geophys%WavGrids%Wav(1:ndat),amf, &
      !   index(1:ndat), Results_Eofpc%Assigned_bins(1:ndat), &
      !   PCA_nbins, PCA_ncnt, PCA_neofs, PCA_bins, pca_fail, PCA_message)

      CALL GEMSTOOL_CreateBins_jbak  (ndat,&
         habsline(1:ndat), absline(1:ndat), Geophys%WavGrids%Wav(1:ndat),amf, &
         index(1:ndat), Assigned_bins(1:ndat), &
         PCA_nbins, PCA_ncnt, PCA_neofs, PCA_bins, pca_fail, PCA_message)
     
   ELSE 
     print * , 'wrong binning index', PCA_binning_index
     stop 1
   ENDIF
   deallocate (taugcum,tauhgcum, taucum)
   n_call_vlidort = SUM(PCA_neofs(1:PCA_nbins)*2)+ PCA_nbins
!  CPU timing
   
   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Eofpctimes(3) = e3 - e2 ! CreateBintime
   endif

!  Dump of bins
   !if ( do_debug_output ) then
   !  open(99,file='BINDUMP_Generic',status='unknown')
   !  write(99,'(A,I2,10(4x,i5,1x))')'Bin Numbers : ',PCA_nbins, PCA_ncnt(0:PCA_nbins-1)
   !  write(99,'(A,2x,10f10.4)')'Bin Limits  : ', bins(0:PCA_nbins-1)
   !  do w = 1, ndat
   !    write(99,'(f8.2, i5, i5, f8.2)') Geophys%wavgrids%wav(w),index(w),Assigned_bins(w), -log(NNtaugcum(w))
   !  enddo
   !  close(99)
   !endif

!  3. FIRST BIN LOOP, PCA Execution
!  ================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  3. PCA Execution. First loop over bins.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Allocation and zeroing

   call CPU_time(e2)

!  Global allocation of Principal Components
   allocate(PrinComps(GT_Maxeofs,ndat))
!  Initialize counting

   istart  = 0
   iend    = 0
   irt_prj = 0

!  Initialize output for this section. Now done outside the PCA Driver routines

   PrinComps = GTZERO

   PCAProj%OPD_PRJ   = GTZERO   ; PCAProj%SSA_PRJ    = GTZERO
   PCAProj%FR_PRJ    = GTZERO   ; PCAProj%FA_PRJ     = GTZERO
   PCAProj%SOLAR_PRJ = GTZERO   ; PCAProj%ALBEDO_PRJ = GTZERO 

   PCAProj%RAY2MOM_PRJ   = GTZERO ; PCAProj%DEP_PRJ = GTZERO
   PCAProj%AERCOEFFS_PRJ = GTZERO ; PCAProj%PF_PRJ  = GTZERO
   PCAProj%N_AERCOEFFS_PRJ = 0

!  timing

   call CPU_time(e3) ; EofpcTimes(4) = e3 - e2 !  Time for Allocation and Zeroing.

!  Start Bin loop
  !print * , geophys%wavgrids%wav(1), geophys%wavgrids%wav(nwav)
   do k = 0, PCA_nbins - 1

!  Progress

      !write(*,*)'Doing PCA for Bin # ',k, PCA_ncnt(k)

!  3a. PCA analysis
!  ----------------

!  Timing
 
      if ( Monitor_CPU ) call cpu_time(e2)

!  set offsets and counts

      ndim   = PCA_ncnt(k) ; k1 = k + 1
      istart = iend; istart_save(k) = istart
      iend   = iend+ndim
  !    IF (which_win == 2) WRITE(*,'(f8.2, f8.2, i5, f8.2)') geophys%wavgrids%wav(istart+1), geophys%wavgrids%wav(iend) ,PCA_ncnt(k), geophys%wavgrids%wav(iend)-geophys%wavgrids%wav(istart+1)
!  Bookkeeping - Number of actual RT binned calculations

!      ndiff_neofp1 = 1 + 2 * PCA_neof           ! THis is the regular value
!      ndiff_neofp1 = 1 + 4 * PCA_neof           ! USE THIS FOR Third-Order Differencing

       PCA_neof = PCA_neofs(k) 
       ndiff_neofp1 = 1 + 2 * PCA_neof        ! Regular value only
       do m = 1, ndiff_neofp1
          irt_prj = irt_prj + 1
          irt_map(irt_prj,1) = k
          irt_map(irt_prj,2) = m
       enddo
!  Strategy 1 : Perform PCA on taudp, omega
!  Strategy 2 : Perform PCA on taudp, rayop

      Call GEMSTOOL_PCACaller &
        ( do_aerosols, do_Sun_Normalized, use_hitran,     & ! Input Flags  
          do_svd_cmp, alb_pcainclude, PCA_Strategy_index, k, PCA_neof,   & ! Input PCA Control  
          nlayers, nlayers2, nlayers21, nlayers22, nstr2, & ! Input numbers
          nmuller, ndat, ndim, istart, index, Geophys,    & ! Binning and Optical Input
          PCAProj, PrinComps,                             & ! outputs
          pca_fail, pca_message, pca_trace_1, pca_trace_2 )   ! Exception handling
!  Linearization. Q = 1 only

      if ( .not. do_aerosols ) then
        do i = 1, ndiff_neofp1
          do n = 1, nlayers
            PCAProj%gas_prj(n,k1,i) = PCAProj%opd_prj(n,k1,i) * ( one - PCAProj%ssa_prj(n,k1,i) )
            if ( LFvary(n) ) then
              do q = 1, LNvary(n) ! xdtau/dx
            ! print * , TRIM(ADJUSTL(profilewf_names_cc(q)))
                !IF (profilewf_names_cc(q) == 'ozone volume mixing ratio------' ) THEN
                    L_opd_prj(q,n,k1,i) = PCAProj%gas_prj(n,k1,i)
                    L_ssa_prj(q,n,k1,i) = - PCAProj%ssa_prj(n,k1,i) * PCAProj%gas_prj(n,k1,i) / PCAProj%opd_prj(n,k1,i)
                    L_fr_prj(q,n,k1,i) = zero ; L_fa_prj(q,n,k1,i) = zero
!             print * , l_opd_prj(q, n, k1, i), '((((((((((((('
                !ENDIF 
              enddo
            endif
          enddo
        enddo
      else
!  PLACEHOLDER WITH AEROSOLS
      endif

!  timing

      if ( Monitor_CPU ) then
         call cpu_time(e3) ; Eofpctimes(5) = Eofpctimes(5) + e3 - e2  ! PCA timing
      endif

!  Exception handling. (If failed, go to 69, where allocation is removed and module returns)

      if ( pca_fail ) then
         messages(Nmessages+1) = '(PCA_Message) - '//Adjustl(TRIM(pca_Message))
         messages(Nmessages+2) = '(PCA_Trace1)  - '//Adjustl(TRIM(pca_trace_1))
         messages(Nmessages+3) = '(PCA_Trace2)  - '//Adjustl(TRIM(pca_trace_2))
         messages(Nmessages+4) = '(PCA_Driver)  - '//'Driver error'
         Nmessages = Nmessages + 4 ; Fail = .true. ; return
      endif

!  continuation point
68   continue

!  End Bin loop, first occasion

   enddo
!  Failure at this point - return

   if ( fail ) return

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  4. RTM (FO/2S/LD) Calculations based on PCA-derived profiles 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!mick mod 4/6/2015 - added extra timing for better comparison with OpenMP
!                    version of code

    if ( Monitor_CPU ) call cpu_time(ser_e1)

!######################################## START PARALLEL REGION ##########

!  Initialize Intensity results (PRJ calculations). Add Q/U Zeroing. 9/23/16
!   12/28/17 Use simplified array declarations.

    if (.not.allocated(lp_ijacs_bin)) then
      allocate (LP_IJacs_Bin (GT_maxatmoswfs,GT_maxlayers,0:GT_maxbins-1,GT_Maxeofs2p1,GT_maxGeometries,3))
    endif

   INTENSITY_BIN = GTZero ; STOKESQU_BIN  = GTZero 
   LP_IJACS_BIN  = GTZero ; LP_QUJACS_BIN = GTZero 
   LS_IJACS_BIN  = GTZero ; LS_QUJACS_BIN = GTZero 

!  ENTER PARALLEL REGION HERE
!  @@@@@@@@@@@@@@@@@@@@@@@@@@

!  Start RT Loop
   !write(*,*) 'Doing RT PRJ calculation'
   do mm = 1, irt_prj

      k = irt_map(mm,1) ;  k1 = k + 1 !    Bin number (Map is zero-based)
      i = irt_map(mm,2)               !    EOF or Exact number within bin

!  Progress

   !   write(*,'(1X,A,I3,A,2I3)')'Doing RT PRJ calculation # ',MM,' ; Bin/EOF = ',K,I

!  4a. Compute (PCA-binned) input optical properties for VLIDORT 
!  -------------------------------------------------------------

      if ( Monitor_CPU ) call cpu_time(e2)

!  "Problem Ray"

      beta2  = PCAProj%ray2mom_prj(k1,i)
      phasmoms_input = GTzero
      phasmoms_input(0,1,1) = GTOne
      phasmoms_input(2,1,1) = beta2
      phasmoms_input(2,2,1) = six    * beta2
      phasmoms_input(2,5,1) = mrtsix * beta2
      if ( nstokes.eq.4) then
         depol_prj = ( GTOne - 2.0d0 * beta2 ) / ( GTOne + beta2 )
         phasmoms_input(1,4,1) =  3.D0 * ( GTOne - 2.D0*depol_prj ) / (GTOne + depol_prj )
      endif

!  set local nmoms for VLIDORT assignations (avoids needless copying if no aerosols)
     LM = local_nmoms

!  Albedo for VLIDORT

      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = PCAProj%albedo_prj(k1,i)
     ! IF (k == 0) print * ,'alb',k1, i, PCAProj%albedo_prj(k1, i)
!  VLIDORT Solar Flux. Now set directly (Optional default to sun-normalized).

      if ( do_Sun_Normalized ) then
         VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = ONE
      else
         VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = PCAProj%solar_prj(k1,i)
      endif

!  VLIDORT optical properties (Use Local results where possible)
!    15 December 2015, Vector implementation

      do n = 1, nlayers
         raywt = PCAProj%fr_prj(n,k1,i)
         VLIDORT_FixIn%Optical%TS_deltau_vert_input(n)  = PCAProj%opd_prj(n,k1,i)
         VLIDORT_ModIn%MOptical%TS_omega_total_input(n) = PCAProj%ssa_prj(n,k1,i)
         if ( aerflag(n) ) then
            aerwt = PCAProj%fa_prj(n,k1,i)
            DO VK = 1, ngkmatc
               GK = gmask(vk) ; rk = rmask(vk) ; ck = cmask(vk) ; sk  = real(smask(vk),dp)
               DO L = 0, local_nmoms
                 MOM = AERWT * SK * PCAProj%AERCOEFFS_PRJ(ck,L,k1,i)
                 IF ( L.LT.3) MOM = MOM + RAYWT*phasmoms_input(L, rk, 1)
                 VLIDORT_FixIn%Optical%TS_Greekmat_total_input(L,N,GK) = MOM
               ENDDO
            enddo
         else
            DO VK = 1, ngkmatc
               GK = gmask(vk); rk = rmask(vk)
               DO L = 0, 2
                  VLIDORT_FixIn%Optical%TS_Greekmat_total_input(L,N,GK) = phasmoms_input(L, rk,1)
               ENDDO
            ENDDO
         endif
         VLIDORT_FixIn%Optical%TS_Greekmat_total_input(0,N,1) = ONE
      enddo

!  VLIDORT Linearized properties. 
!   NEEDS WORK FOR INPUTS with AEROSOLS.....
!     No linearized setup (yet) for aerosols -   PLACEHOLDER

      do n = 1, nlayers
         if ( LFvary(n) ) then
           do q = 1, LNvary(n)
             VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(q,n) = L_opd_prj(q,n,k1,i) / PCAProj%opd_prj(n,k1,i)
             VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(q,n) = L_ssa_prj(q,n,k1,i) / PCAProj%ssa_prj(n,k1,i)
           enddo
         endif
      enddo
      
!  time

      if (monitor_CPU) then
        call cpu_time(e3) ; Eofpctimes(6) = Eofpctimes(6) + e3 - e2  !  EofpcOpTime_PCAProj_VLIDORT
      endif

!  4b. Compute (PCA-binned) optical properties for FO Model 
!  --------------------------------------------------------

!  initialize CPU monitoring

      if ( Monitor_CPU ) call cpu_time(e2)

!  FO albedos. Copy VLIDORT

      FO_reflec(1,1,1:FO_ngeoms) = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
      if ( FO_do_reflecwfs ) FO_LS_reflec(1,1,1,1:FO_ngeoms) = one

!  FO Solar Flux. Now set directly (Optional default to sun-normalized).Sufficient to copy VLIDORT input

      FO_FLUX  = 0.25_gtpk * VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR  / FO_Pie

!  FO layer Bulk optical properties (New Style). Use VLIDORT input.
!     Corrected Code 5/2/17, Takes proper care with Delta-M scaling 
 
      do n = 1, FO_nlayers
         FO_omegas(n)  = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
         if ( FO_do_deltam ) then
            FO_truncfac(n) = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(nstr2,n,1)/dnm1
            fact1 = GTOne - FO_truncfac(n) * FO_omegas(n)
            FO_deltaus(n) = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n) * fact1
         else
            FO_truncfac(n) = GTZero
            FO_deltaus(n)  = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n)
         endif
         FO_extinction(n) = FO_deltaus(n) / FO_diffgrid(n)
         FO_Greekmat(0:local_nmoms,n,:) = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:local_nmoms,n,:)
      enddo

!  Linearized bulk properties. (everything pre-zeroed)
!     Corrected Code 5/2/17, Takes proper care with Delta-M scaling 

      do n = 1, FO_nlayers
         if ( FO_do_deltam ) then
            fact1 = one - FO_truncfac(n) * FO_omegas(n)
            if ( LFvary(n) ) then
               do q = 1, LNvary(n)
                  FO_L_truncfac(n,q) = FO_L_Greekmat(nstr2,n,1,q)/dnm1
                  L_fact1 = - FO_truncfac(n) * L_opd_prj(q,n,k1,i) - FO_L_truncfac(n,q) * FO_omegas(n)
                  FO_L_deltaus(n,q)    = L_opd_prj(q,n,k1,i) * fact1 + PCAProj%opd_prj(n,k1,i) * L_fact1
                  FO_L_omegas(n,q)     = L_ssa_prj(q,n,k1,i)
                  FO_L_extinction(n,q) = FO_L_deltaus(n,q) / FO_diffgrid(n)
               enddo
            endif
         else
            if ( LFvary(n) ) then
               do q = 1, LNvary(n)
                  FO_L_deltaus(n,q)    = L_opd_prj(q,n,k1,i)
                  FO_L_omegas(n,q)     = L_ssa_prj(q,n,k1,i)
                  FO_L_extinction(n,q) = FO_L_deltaus(n,q) / FO_diffgrid(n)
               enddo
            endif
         endif
       enddo

!   1/20/16. Solar Flux. Now set directly (Optional default to sun-normalized).
!                        Sufficient to copy VLIDORT input

       FO_FLUX  = 0.25_gtpk * VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR  / FO_Pie

!  Phase Matrix calculation. Use the VLIDORT Greekmat array.

       if ( do_Jacobians ) then
         Call FO_VectorSS_PhasMat_Plus &
         ( FO_MaxGeometries, FO_maxlayers, FO_maxmoments, FO_maxatmoswfs, FO_do_profilewfs, & ! Dimensions/Flag
           FO_do_sunlight, FO_do_deltam, FO_Lvaryflags, FO_Lvarynums, FO_Lvarymoms,         & ! Flags/Linearization
           FO_nstokes, FO_ngeoms, FO_nlayers, local_nmoms, FO_aclevel,                      & ! Numbers
           FO_omegas, FO_truncfac, FO_Greekmat, FO_L_omegas, FO_L_truncfac, FO_L_Greekmat,  & ! optical
           FO_GenSpher, FO_Rotations, FO_ExactScat_up, FO_L_ExactScat_up )                    ! Functions/Output
       else
         Call FO_VectorSS_PhasMat &
         ( FO_MaxGeometries, FO_maxlayers, FO_maxmoments, FO_do_sunlight, FO_do_deltam,   & ! Dimensions/Flags
           FO_nstokes, FO_ngeoms, FO_nlayers, local_nmoms, FO_aclevel,                    & ! Numbers
           FO_omegas, FO_truncfac, FO_Greekmat,                                           & ! Optical input
           FO_GenSpher, FO_Rotations, FO_ExactScat_up )                                     ! Functions/Output
       endif

!do n = 1, 74
!   write(*,*)FO_ExactScat_up(n,1,1,1),FO_L_ExactScat_up(n,1,1,1,1)!,FO_omegas(n), FO_truncfac(n), FO_Greekmat(0:2,n,1)
!enddo
!stop'FO Exactcat'

!  time

      if (monitor_CPU) then
         call cpu_time(e3) ; Eofpctimes(7) = Eofpctimes(7) + e3 - e2  !  EofpcOpTime_PCAProj_FO
      endif

!  4c. Compute (PCA-binned) 2stream optical properties
!  ---------------------------------------------------

!  CPU monitoring

      if (monitor_CPU) call cpu_time(e2)

!  2stream optical properties (mostly, copy the VLIDORT values)

      do n = 1, nlayers
         deltau_input(n) = VLIDORT_FixIn%Optical%TS_deltau_vert_input(n)
         omega_input(n)  = VLIDORT_ModIn%MOptical%TS_omega_total_input(n)
         ASYMM_INPUT(n)  = VLIDORT_FixIn%Optical%TS_greekmat_total_input(1,n,1)/3.d0 
         D2S_SCALING(n)  = VLIDORT_FixIn%Optical%TS_greekmat_total_input(2,n,1)/5.d0
      enddo

!  Linearized properties

      do n = 1, nlayers
         if ( LFvary(n) ) then
            do q = 1, LNvary(n)
               L_deltau_input(n,q) = deltau_input(n) * VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(q,n) 
               L_omega_input(n,q)  = omega_input(n)  * VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(q,n)
               L_ASYMM_INPUT(n,q)  = ASYMM_INPUT(n)  * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,1,n,1)  ! Zero for now
               L_D2S_SCALING(n,q)  = D2S_SCALING(n)  * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,2,n,1)  ! Zero for now
            enddo
         endif
      enddo

!  surface and solar

      ALBEDO      = VLIDORT_FixIn%Optical%TS_lambertian_albedo
      FLUX_FACTOR = VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR
!  time

      if (monitor_CPU) then
         call cpu_time(e3) ; Eofpctimes(8) = Eofpctimes(8) + e3 - e2  !  EofpcOpTime_PCAProj_2S
      endif

!  4d. VLIDORT call
!  ----------------

!  First-Order optical properties Now part of VLIDORT
!    ****** REQUIRED IF YOU ARE DOING 3-Intensity correction

      if ( .not. do_3M_Correction  ) stop 'Must do 3M for Stokes Q/U'

!  time

      if (monitor_CPU) call cpu_time(e2)

!  Debug

      DO_Vlidort_Inpdebug = .false. ! ; IF ( w.eq.40 ) DO_Vlidort_Inpdebug = .true.

!  Call

      If ( do_Jacobians ) then
         CALL VLIDORT_LPS_master ( do_debug_input,& 
           VLIDORT_FixIn,    & ! INPUTS
           VLIDORT_ModIn,    & ! INPUTS (possibly modified)
           VLIDORT_Sup,      & ! INPUTS/OUTPUTS
           VLIDORT_Out,      & ! OUTPUTS
           VLIDORT_LinFixIn, & ! INPUTS
           VLIDORT_LinModIn, & ! INPUTS (possibly modified)
           VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
           VLIDORT_LinOut )    ! OUTPUTS
      else
         CALL VLIDORT_master ( do_debug_input, & 
           VLIDORT_FixIn,    & ! INPUTS
           VLIDORT_ModIn,    & ! INPUTS (possibly modified)
           VLIDORT_Sup,      & ! INPUTS/OUTPUTS
           VLIDORT_Out )       ! OUTPUTS
      endif

!  Exception handling for Input checks. Check on indexing now added

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_SERIOUS ) then
         DO M = 1, VLIDORT_Out%Status%TS_NCHECKmessages
           messages(Nmessages+2*M-1) = '(VLIDORT_2p7 Message) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_CHECKmessages(M)))
           messages(Nmessages+2*M  ) = '(VLIDORT_2p7_Action ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_ACTIONS(M)))
         ENDDO
         Nmessages = Nmessages + 2*VLIDORT_Out%Status%TS_NCHECKmessages
         write(C3,'(I3)')mm ; messages(Nmessages+1) = '(VLIDORT Input Check, RT_PRJ # = '//C3
         Nmessages = Nmessages + 1 ; Fail = .true.;  go to 69
      endif

!  Exception handling for Calculation. 5/6/15 upgrade, added wavelength number

      if ( VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         messages(Nmessages+1) = '(VLIDORT_2p7 Message) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_MessAGE))
         messages(Nmessages+2) = '(VLIDORT_2p7_Trace ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_trace_1))
         messages(Nmessages+3) = '(VLIDORT_2p7_Trace ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_trace_2))
         messages(Nmessages+4) = '(VLIDORT_2p7_Trace ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_trace_3))
         Nmessages = Nmessages + 4
         write(C3,'(I3)')mm ; messages(Nmessages+1) = '(VLIDORT Execution, RT_PRJ # = '//C3
         Nmessages = Nmessages + 1 ; Fail = .true.;  go to 69
      endif

!  High bin (VLIDORT) results. Add Q/U components (2016), upgraded, 12/28/17

      do v = 1, ngeoms
         intensity_Bin(k,i,v,1) = VLIDORT_Out%Main%TS_Stokes(1,v,1,dir) 
         if ( nstokes.gt.1 ) then
            StokesQU_Bin(1:2,k,i,v,1)   = VLIDORT_Out%Main%TS_Stokes(1,v,2:3,dir)
         endif
      enddo

!  Intensities, Stokes QU results. OLD CODE......

!      do v = 1, ngeoms
!         Isurf(1:ns) = 0.0d0 ; if ( dir .eq. UPIDX ) Isurf =  VLIDORT_Sup%SS%TS_STOKES_DB(1,v,1:ns)
!         Intensity_bin(k,i,v,3) = VLIDORT_Sup%SS%TS_Stokes_SS(1,v,1,dir) + ISurf(1)
!         intensity_Bin(k,i,v,1) = VLIDORT_Out%Main%TS_Stokes(1,v,1,dir) - Intensity_bin(k,i,v,3)
!         if ( nstokes.gt.1 ) then
!            StokesQU_bin(1:2,k,i,v,2)   = VLIDORT_Sup%SS%TS_Stokes_SS(1,v,2:3,dir) + ISurf(2:3)
!            StokesQU_Bin(1:2,k,i,v,1)   = VLIDORT_Out%Main%TS_Stokes(1,v,2:3,dir) - StokesQU_bin(1:2,k,i,v,2)
!         endif
!      enddo

!  Linearized Profile Jacobian results

      if ( do_profile_WFs ) then
        do v = 1, ngeoms
          do n = 1, nlayers
            if ( LFvary(n) ) then
              do q = 1, LNvary(n)
                MSJac(1:ns) = VLIDORT_LinOut%Prof%TS_ProfileWF(q,n,1,v,1:ns,DIR)
                !print *, n,q, k, i,  MSjac(1)
                LP_IJacs_Bin(q,n,k,i,v,1) = MSJac(1)
                if ( nstokes.gt.1 ) then
                  LP_QUJacs_Bin(1:2,q,n,k,i,v,1) = MSJac(2:3)    ! LD result    
                endif
              enddo
            endif
          enddo
        enddo
      endif

     
!  Linearized Surface Jacobian results

      if ( do_surface_WFs ) then
        do v = 1, ngeoms
          q = 1 ; MSJac(1:ns) = VLIDORT_LinOut%Surf%TS_Surfacewf(q,1,v,1:ns,DIR)
          !print *, q, k, i,  MSjac(1)
          LS_IJacs_Bin(k,i,v,1) = MSJac(1)
          if ( nstokes.gt.1 ) then
            LS_QUJacs_Bin(1:2,k,i,v,1) = MSJac(2:3) 
          endif
        enddo
      endif

!  Timing

      if (monitor_CPU) then
         call cpu_time(e3) ; Eofpctimes(9) = Eofpctimes(9) + e3 - e2  ! binvlidorttime_4d
      endif

!  4e. FO CALL
!  -----------

!  Time

      if (monitor_CPU) call cpu_time(e2)

!  Call

      if ( do_Jacobians ) then
        Call SSV_Integral_ILPS_UP_Optimized &
          ( FO_maxgeometries, FO_maxlayers, FO_maxfine, FO_maxatmoswfs, FO_maxsurfacewfs,             & ! Inputs (dimensioning)
           FO_do_planpar, FO_do_regular_ps, FO_do_enhanced_ps, FO_doNadir, FO_do_sleave,              & ! Inputs (Flags)
           FO_do_sunlight, FO_do_lambertian, FO_do_profilewfs, FO_do_reflecwfs, FO_do_sleavewfs,      & ! Inputs (control, Jacobian)
           FO_Lvaryflags, FO_Lvarynums, FO_n_reflecwfs, FO_n_sleavewfs,                               & ! Inputs (control, Jacobian)
           FO_nstokes, FO_ngeoms, FO_nlayers, FO_nfinedivs, FO_AcLevel,                               & ! Inputs (control output)
           FO_reflec, FO_slterm, FO_extinction, FO_deltaus, FO_exactscat_up, FO_flux, FO_fluxvec,     & ! Inputs (Optical)
           FO_LS_reflec, FO_LSSL_SLterm, FO_L_extinction, FO_L_deltaus, FO_L_exactscat_up,            & ! Inputs (Optical - Lin)
           FO_Mu0, FO_Mu1, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,                      & ! Inputs (Geometry)
           FO_Raycon, FO_cota, FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,        & ! Inputs (Geometry)
           FO_stokes_up, FO_stokes_db, FO_LP_Jacobians_up, FO_LP_Jacobians_db, FO_LS_Jacobians_db )     ! Output
      else
        Call SSV_Integral_I_UP_Optimized &
         ( FO_maxgeometries, FO_maxlayers, FO_maxfine,                                                    & ! Inputs (dimension)
           FO_do_planpar, FO_do_regular_ps, FO_do_enhanced_ps, FO_doNadir, FO_do_sleave,                  & ! Inputs (Flags)
           FO_do_sunlight, FO_do_lambertian, FO_nstokes, FO_ngeoms, FO_nlayers, FO_nfinedivs, FO_AcLevel, & ! Inputs (ctrl output)
           FO_reflec, FO_slterm, FO_extinction, FO_deltaus, FO_exactscat_up, FO_flux, FO_fluxvec,         & ! Inputs (Optical)
           FO_Mu0, FO_Mu1, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,                          & ! Inputs (Geometry)
           FO_Raycon, FO_cota, FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,            & ! Inputs (Geometry)
           FO_stokes_up, FO_stokes_db, FO_cumsource_up )                                                    ! Outputs
      endif

!  Saved results

      do v = 1, ngeoms
         FORad(1:ns)= FO_stokes_up(1:ns,v)
         if ( dir .eq. UPIDX ) FORad(1:ns) =  FORad(1:ns) + FO_stokes_db(1:ns,v)
         Intensity_bin(k,i,v,3) = FORad(1)
         if ( nstokes.gt.1 ) then
            StokesQU_bin(1:2,k,i,v,2)   = FORad(2:3)
         endif
      enddo

!  Linearized Profile Jacobian results

      if ( do_profile_WFs ) then
        do v = 1, ngeoms
          do n = 1, nlayers
            if ( LFvary(n) ) then
              do q = 1, LNvary(n)
                FOJac(1:ns)  = FO_LP_Jacobians_up(1:ns,v,n,q)
                IF ( DIR.eq.UpIdx ) FOJac(1:ns) = FOJac(1:ns) + FO_LP_Jacobians_db(1:ns,v,n,q)
                LP_IJacs_Bin(q,n,k,i,v,3) = FOJac(1)
                if ( nstokes.gt.1 ) then
                  LP_QUJacs_Bin(1:2,q,n,k,i,v,2) = FOJac(2:3)    ! FO result
                endif
              enddo
            endif
          enddo
        enddo
      endif

!  Linearized Surface Jacobian results.

      if ( do_surface_WFs ) then
        do v = 1, ngeoms
          q = 1 ; FOJac(1:ns)  = GTzero ; IF ( DIR.eq.UpIdx ) FOJac(1:ns) = FO_LS_Jacobians_db(1:ns,v,q)
          LS_IJacs_Bin(k,i,v,3) = FOJac(1)
          if ( nstokes.gt.1 ) then
             LS_QUJacs_Bin(1:2,k,i,v,2) = FOJac(2:3)
          endif
        enddo
      endif

!  Timing

      if (monitor_CPU) then
         call cpu_time(e3) ; Eofpctimes(10) = Eofpctimes(10) + e3 - e2  ! binFOtime_4e
      endif

!  4f.  2stream binned calculation
!  -------------------------------

!  timing

      if (monitor_CPU) call cpu_time(e2)

!  Call
!  Rob Fix 8/3/18. New Calling statements for multiple geometries
      if ( do_Jacobians ) then
         CALL TWOSTREAM_LPS_MASTER &
          ( S2_MAXLAYERS, S2_MAXTOTAL, S2_MAXGEOMS, S2_MAXMESSAGES, S2_NLAYERS, S2_NTOTAL, S2_NGEOMS, & ! Inputs
            S2_MAX_ATMOSWFS, S2_MAX_SURFACEWFS, S2_MAX_SLEAVEWFS,             & ! Dimensions
            DO_PLANE_PARALLEL, DO_D2S_SCALING, TOPMOST_LEVEL,                 & ! Inputs
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,             & ! Inputs
            DO_PDINVERSE, BVPINDEX, BVPSCALEFACTOR, FLUX_FACTOR,              & ! Inputs
            DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs 
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,  & ! Inputs
            STREAMVAL, STREAMINV, GEOMETRIES, AUX_GEOMS, CHAPFACS,            & ! Inputs
            DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, ALBEDO,      & ! Inputs
            BRDF_F_0, BRDF_F, UBRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,          & ! Inputs
            L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,      & ! Inputs
            LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LSSL_ISOTROPIC, LSSL_F_0,     & ! Inputs
            INTENSITY_2S, PROFILEWFS_2S, SURFACEWFS_2S,                       & ! Outputs
            STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,            & ! Exception handling
            STATUS_EXECUTION, S2_MESSAGE, S2_TRACE_1, S2_TRACE_2 )              ! Exception handling
      else
         Call TWOSTREAM_MASTER &
          ( S2_MAXLAYERS, S2_MAXTOTAL, S2_MAXGEOMS, S2_MAXMESSAGES, S2_NLAYERS, S2_NTOTAL, S2_NGEOMS, & ! Inputs
            DO_PLANE_PARALLEL, DO_D2S_SCALING, TOPMOST_LEVEL,                 & ! Inputs
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,             & ! Inputs
            DO_PDINVERSE, BVPINDEX, BVPSCALEFACTOR, FLUX_FACTOR,              & ! Inputs
            STREAMVAL, STREAMINV, GEOMETRIES, AUX_GEOMS, CHAPFACS,            & ! Inputs
            DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, ALBEDO,      & ! Inputs
            BRDF_F_0, BRDF_F, UBRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,          & ! Inputs
            INTENSITY_2S, STATUS_INPUTCHECK, STATUS_EXECUTION,                & ! Outputs 
            C_NMESSAGES, C_MESSAGES, C_ACTIONS, S2_MESSAGE, S2_TRACE_1, S2_TRACE_2 )
      endif
!  Exception handling for Input checks. 5/6/15 upgrade, added wavelength number

      if ( STATUS_INPUTCHECK .eq. 1 ) then
         DO M = 1, C_NMessages
            Messages(NMessages+2*M-1) = '(2STREAM_2p4 Message) '//Adjustl(TRIM(C_Messages(M-1)))
            Messages(NMessages+2*M)   = '(2STREAM_2p4 action ) '//Adjustl(TRIM(C_ACTIONS(M-1)))
         ENDDO
         NMessages = NMessages + 2*C_NMessages
         write(C5,'(I5)')w ; Messages(NMessages+1) = '(Optional 2S Input  Check, RT Exact # = '//C5//')'
         NMessages = NMessages + 1 ; fail = .true.;  go to 69
      endif

!stop'after first 2S'

!  Exception handling for Calculation

      IF ( STATUS_EXECUTION .eq. 1 ) THEN
         messages(Nmessages+1)   = '(2STREAM_2p4 Message) '//Adjustl(TRIM(S2_Message))
         messages(Nmessages+2)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(S2_Trace_1))
         messages(Nmessages+3)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(S2_Trace_2))
         Nmessages = Nmessages + 3
         write(C3,'(I3)')mm ; messages(Nmessages+1) = '(Twostream Execution, RT_PRJ # = '//C3
         Nmessages = Nmessages + 1 ; Fail = .true.;  go to 69
      endif

!  Lo bin (2stream) results. Up only, 1 geometry
!   Rob Fix 8/3/18. Add q = 1 for surface WFs. Index should be k not q
!   Rob Fix 8/3/18. Allowed multiple geometries

      do v = 1, ngeoms
        intensity_Bin(k,i,v,2) = intensity_2S(v)
       !  print *, k,mm, intensity_bin(k, i, v, 1:3) 
        if ( do_profile_WFs ) then
          do n = 1, nlayers
            if ( LFvary(n) ) then
              do q = 1, LNvary(n)
                LP_IJacs_Bin(q,n,k,i,v,2) = PROFILEWFS_2S(v,n,q)
              enddo
            endif
          enddo
        endif
        if ( do_surface_Wfs ) then
          q = 1 ; LS_IJacs_Bin(k,i,v,2)  = SURFACEWFS_2S(v,q)
        endif
      enddo

!  Timing

      if (monitor_CPU) then
         call cpu_time(e3) ; Eofpctimes(11) = Eofpctimes(11) + e3 - e2  ! bin2stime_4f
      endif

!  Continuation point  for errors

69    continue
!  End of RT Calculation loop  (EOF projections)
   enddo

!pause'End of Step 1'

!  END PARALLEL REGION HERE, RT based on PCA-derived profiles
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!mick fix 4/6/2015 - moved from down below and modified
!  Time spent in equivalent parallel region #1

   if (monitor_CPU) then
      call cpu_time(ser_e2) ; SerEofpcTime = ser_e2 - ser_e1
   endif
!  Return if failure. This must be done outside parallel region and bin loop
   if ( Fail ) Return

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  5. Correction Factors (third Bin Loop)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialize

   istart = 0
   iend   = 0

!  Start Bin Loop
   !write(*,*) 'Doing PCA corrections'
   do k = 0, PCA_nbins - 1

!  Progress

      !write(*,'(1X,A,I3)')'Doing PCA corrections, Bin = ',K

!  set offsets

      ndim   = PCA_ncnt(k)
      istart = iend
      iend   = iend+ndim

!  Special case, PLACEHOLDER
!      if (SPECIAL_FLAG(k).eqv..true.) go to 555

!  Number of EOFs

      PCA_neof = PCA_neofs(k) 
      ndiff_neofp1 = 1 + 2 * PCA_neof

!  Local allocate.

      allocate(Princomps_Local      (PCA_neof,ndim))

      allocate(Intensity_Corrfacs_Local(ndim,ngeoms))
      allocate(Intensity_Bin_Local      (ndiff_neofp1,ngeoms,3))
      if (nstokes.gt.1 ) then
        allocate(StokesQU_Bin_Local       (2,ndiff_neofp1,ngeoms,2))
        allocate(StokesQU_Corrfacs_Local  (2,ndim,ngeoms))
      endif

      if ( do_profile_wfs ) then
        allocate(LP_IJacs_bin_Local(npars,nlayers,ndiff_neofp1,ngeoms,3))
        allocate(LP_IJacs_Corrfacs_Local(npars,nlayers,ndim,ngeoms))
        if (nstokes.gt.1 ) then
!           allocate(LP_QUJacs_Corrfacs_Local(npars,nlayers,ndim,ngeoms,2))
           allocate(LP_QUJacs_Corrfacs_Local(2,npars,nlayers,ndim,ngeoms))
           allocate(LP_QUJacs_bin_Local(2,npars,nlayers,ndiff_neofp1,ngeoms,2))
        endif
      endif
      if ( do_surface_wfs ) then
        allocate(LS_IJacs_bin_Local(ndiff_neofp1,ngeoms,3))
        allocate(LS_IJacs_Corrfacs_Local(ndim,ngeoms))
        if (nstokes.gt.1 ) then
           allocate(LS_QUJacs_bin_Local(2,ndiff_neofp1,ngeoms,2))
           allocate(LS_QUJacs_Corrfacs_Local(2,ndim,ngeoms))
!           allocate(LS_QUJacs_Corrfacs_Local(ndim,ngeoms,2))
        endif
      endif

!  5a. Calculating ratios, logarithms and differences (2010 paper, Equations 2-7)

!  timing

      if (monitor_CPU) call cpu_time(e2)     

!  Assign Local Princomps and Intensities/StokesQU. Simple copying, no real overhead.....

      do m = 1, PCA_neof
        do i = 1, ndim
          ii = index(istart+i)
          PrinComps_Local(m,i) = PrinComps(m,ii) 
         enddo
      enddo

      do m = 1, ndiff_neofp1
        do v = 1, ngeoms
          Intensity_Bin_Local(m,v,1:3) = Intensity_Bin(k,m,v,1:3)             ! I-component, LD/2S/FO
          !IF (m ==1) print * , intensity_bin_local(m, v, 1:3)
          if (nstokes.gt.1 ) then
            StokesQU_Bin_Local(1,m,v,1:2) = StokesQU_bin(1,k,m,v,1:2)         ! Q-component, LD/FO
            StokesQU_Bin_Local(2,m,v,1:2) = StokesQU_bin(2,k,m,v,1:2)         ! U-component, LD/FO
          endif
        enddo
      enddo

!  Assign Profile linearization local arrays

      if ( do_Profile_wfs ) then
        do m = 1, ndiff_neofp1
          do v = 1, ngeoms
            do n = 1, nlayers
              if ( LFvary(n) ) then
                do q = 1, LNvary(n)
                  LP_IJacs_bin_Local(q,n,m,v,1:3) = LP_IJacs_bin(q,n,k,m,v,1:3)
!IF(m ==1 .and. q ==1)  print * , LP_IJacs_bin_Local(q,n,m,v,1:3) , q, n

                  if (nstokes.gt.1 ) then
                    LP_QUJacs_bin_Local(1,q,n,m,v,1:2) = LP_QUJacs_bin(1,q,n,k,m,v,1:2)  ! Stokes-Q
                    LP_QUJacs_bin_Local(2,q,n,m,v,1:2) = LP_QUJacs_bin(2,q,n,k,m,v,1:2)  ! Stokes-U
                  endif
!if ( n.eq.1)write(68,*)k,m,q,n,v,LP_QUJacs_bin(1,q,n,k,m,v,1:2),LP_QUJacs_bin(2,q,n,k,m,v,1:2)
                enddo
              endif
            enddo
          enddo
        enddo
      endif

!  Suggested scattering weight code
!                  if ( q.eq.1 ) then
!                    LP_Scatwts_LD_bin_Local(n,m,v) = LP_Scatwts_LD_bin(n,k,m,v)
!                    LP_Scatwts_2S_bin_Local(n,m,v) = LP_Scatwts_2S_bin(n,k,m,v)
!                    LP_Scatwts_FO_bin_Local(n,m,v) = LP_Scatwts_FO_bin(n,k,m,v)
!                  endif

!  Assign Surface linearization local arrays. Only 1 Jacobian

      if ( do_Surface_wfs ) then
        do m = 1, ndiff_neofp1
          do v = 1, ngeoms
            LS_IJacs_bin_Local(m,v,1:3) = LS_IJacs_bin(k,m,v,1:3)              ! LD/2S/FO
            if (nstokes.gt.1 ) then
              LS_QUJacs_bin_Local(1:2,m,v,1) = LS_QUJacs_bin(1:2,k,m,v,1)      ! LD 
              LS_QUJacs_bin_Local(1:2,m,v,2) = LS_QUJacs_bin(1:2,k,m,v,2)      ! FO
            endif
          enddo
        enddo
      endif

!  Calculate Correction Factors
!    Correction for Q/U, added 5/20/16 (Q), 9/23/16(U). Merged, December 2017

      if ( do_3M_correction ) then

!  Stokes-I and Jacobians, Stokes QU

        if ( do_Profile_wfs ) then
          if ( do_surface_wfs ) then
            call pca_3M_correction_LPS &
            ( PCA_neof, ndiff_neofp1, ndim, ngeoms, nlayers, npars, nspars, PrinComps_Local,             & ! Inputs (control, PCs)
              Intensity_bin_Local, LP_IJacs_bin_Local, LS_IJacs_bin_Local, LP_Jacs_Exist, LS_Jacs_Exist, & ! Inputs (BinI/BinJacs)
              Intensity_Corrfacs_Local, LP_IJacs_Corrfacs_Local, LS_IJacs_Corrfacs_Local )                 ! Output (Corrections)
            if ( nstokes.gt.1 ) then
              Call pca_QU_correction_LPS &
               ( PCA_neof, ndiff_neofp1, ndim, ngeoms, nlayers, npars, nspars, PrinComps_Local,              & ! In (ctrl, PCs)
                 StokesQU_bin_Local, LP_QUJacs_bin_Local, LS_QUJacs_bin_Local, LP_Jacs_Exist, LS_Jacs_Exist, & ! In (BinQU/BinJacs)
                 StokesQU_Corrfacs_Local, LP_QUJacs_Corrfacs_Local, LS_QUJacs_Corrfacs_Local )                 ! Out (Corrections)
            endif
          else
            call pca_3M_correction_LP &
            ( PCA_neof, ndiff_neofp1, ndim, ngeoms, nlayers, npars, PrinComps_Local,  & ! Inputs (control, PCs)
              Intensity_bin_Local, LP_IJacs_bin_Local, LP_Jacs_Exist,                 & ! Inputs (BinI/BinJacobians)
              Intensity_Corrfacs_Local, LP_IJacs_Corrfacs_Local )                       ! Output (Corrections)
            if ( nstokes.gt.1 ) then
              Call pca_QU_correction_LP &
               ( PCA_neof, ndiff_neofp1, ndim, ngeoms, nlayers, npars, PrinComps_Local, & ! Inputs (control, PCs)
                 StokesQU_bin_Local, LP_QUJacs_bin_Local, LP_Jacs_Exist,                & ! Inputs (BinQU/BinJacobians)
                 StokesQU_Corrfacs_Local, LP_QUJacs_Corrfacs_Local )                      ! Output (Corrections)
            endif
          endif
        else
          call  pca_3M_correction &
          ( PCA_neof, ndiff_neofp1, ndim, ngeoms, PrinComps_Local, & ! Inputs (control, PCs)
            Intensity_Bin_Local, Intensity_Corrfacs_Local )          ! Input/Output (BinI/Corrections)
          if ( nstokes.gt.1 ) then
            call  pca_QU_correction &
              ( PCA_neof, ndiff_neofp1, ndim, ngeoms, PrinComps_Local,  & ! Inputs (control, PCs)
                StokesQU_bin_Local, StokesQU_Corrfacs_Local )             ! Input/Output (BinQU/Corrections)
          endif
        endif
      endif

!  2D Corrections. PLACEHOLDER
!      if ( .not. do_3M_correction ) then
!        call  pca_2M_correction &
!          ( PCA_neof, ndiff_neofp1, ndim, ngeoms, PrinComps_Local, & ! Inputs (control, PCs)
!            Intensity_Bin_Local, Intensity_Corrfacs_Local )          ! Input/Output (BinI/Corrections)
!      endif

!  Use the following for 3M, third-order differencing. WITH CAUTION.
!      if ( .not. do_SwitchCorr .and. do_3M_correction .and. do_3orderdiff ) then
!         call  pca_3M_correction_3OD &
!          ( PCA_neof, neofs_4p1, ndim, ngeoms, PrinComps_Local, & ! Inputs (control, PCs)
!            Intensity_Bin_Local, Intensity_Corrfacs_Local )       ! Input/Output (BinI/Corrections)
!      endif

!  Save Correction factor results. Add Q 5/20/16, Add U component, 9/23/16.
      do i = 1, ndim
        ii = index(istart+i)
        do v = 1, ngeoms
          Stokes_Corrfacs(v,1,ii) = intensity_corrfacs_Local(i,v)
          !IF (k == 0) print * ,i, ii, Stokes_Corrfacs(v,1,ii) 
        enddo
      enddo
      if ( do_3M_correction .and. nstokes.gt.1 ) then
        do i = 1, ndim
          ii = index(istart+i)
          do v = 1, ngeoms
            Stokes_Corrfacs(v,2:3,ii) = StokesQU_corrfacs_Local(1:2,i,v) 
          enddo
        enddo
      endif

!  Save Correction factor results. Profile Jacobians

      if ( do_profile_wfs ) then
        do i = 1, ndim
          ii = index(istart+i)
          do v = 1, ngeoms
            do n = 1, nlayers
              if ( LFvary(n) ) then
                do q = 1, LNvary(n)
                  Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,1,ii) = LP_IJacs_Corrfacs_Local(q,n,i,v)
                  if ( nstokes.gt.1 ) then
                    Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,2:3,ii) = LP_QUJacs_Corrfacs_Local(1:2,q,n,i,v)
!                    Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,2:3,ii) = LP_QUJacs_Corrfacs_Local(q,n,i,v,1:2)
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      endif

!  Save Correction factor results. surface Jacobians (Only 1 of them

      if ( do_surface_wfs ) then
        q = 1
        do i = 1, ndim
          ii = index(istart+i)
          do v = 1, ngeoms
            Results_LS_Eofpc%LS_Jacobians_Corrfacs(q,v,1,ii) = LS_IJacs_Corrfacs_Local(i,v)
            !print * , ii, LS_IJacs_corrfacs_local(i,v)
            if ( nstokes.gt.1 ) then
              Results_LS_Eofpc%LS_Jacobians_Corrfacs(q,v,2:3,ii) = LS_QUJacs_Corrfacs_Local(1:2,i,v)
            endif
          enddo
        enddo
      endif

!  Time

      if (monitor_CPU) then
         call cpu_time(e3) ; Eofpctimes(12) = Eofpctimes(12) + e3 - e2  ! eofcortime_5a
      endif

!  Local de-allocate.

      deallocate(PrinComps_Local,Intensity_CorrFacs_Local,Intensity_Bin_Local)
      if ( do_Profile_wfs ) deallocate(LP_IJacs_CorrFacs_Local,LP_IJacs_Bin_Local)
      if ( do_Surface_wfs ) deallocate(LS_IJacs_CorrFacs_Local,LS_IJacs_Bin_Local)
      if ( nstokes.gt.1) then
         deallocate(StokesQU_CorrFacs_Local,StokesQU_Bin_Local)
         if ( do_Profile_wfs ) deallocate(LP_QUJacs_CorrFacs_Local,LP_QUJacs_Bin_Local)
         if ( do_Surface_wfs ) deallocate(LS_QUJacs_CorrFacs_Local,LS_QUJacs_Bin_Local)
      endif
!  End third bin loop

   enddo

!  Global de-allocate

   deallocate(PrinComps)
   if ( do_Profile_WFs ) deallocate(LP_Jacs_Exist)
   if ( do_Surface_WFs ) deallocate(LS_Jacs_Exist)

!  return if Fail

   if ( Fail ) return

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  6. Fast 2S/FO RTM Calculations (4th Bin loop)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!mick fix 4/6/2015 - moved from up above

   if ( Monitor_CPU ) call cpu_time(ser_e1)

!  Start Bin Loop
   !write(*,*) 'Doing Fast RT Calculations'
   do k = 0, PCA_nbins - 1

!  Progress

!      write(*,'(1X,A,I3)')'Doing Fast RT Calculations, Bin = ',K

!  set offsets

      ndim   = PCA_ncnt(k)
      istart = istart_save(k)

!  Calculation NOT DONE if special_flag set. or if fast calculation not required
!      if ( Special_Flag(k) ) go to 556

      if ( .not. do_fast_calculation ) go to 556

!@@@@@@@@@@@@@@@@@ ENTER PARALLEL REGION HERE
!
!  Full 2-stream calculations (All wavelengths)
!  ==========================

!  Start loop over all points within bin
     !print * ,k,ndim, Geophys%WavGrids%wav(index(istart+1)),geophys%wavgrids%wav(index(istart+ndim))
      do i = 1, ndim

!  Progress

         !if ( mod(i,20).eq.0) write(*,*)' -- Fast Calculation #', i

!  6a. Assign full-monochromatic (non-binned) properties, FO Model 
!  ---------------------------------------------------------------

!  initialize CPU monitoring

         if ( Monitor_CPU ) call cpu_time(e2)

!  Index in the original data set

         ii = index(istart+i)
         !print * , ii, Geophys%WavGrids%wav(ii)
!  "Problem Ray"

         depol = Geophys%Xsecs%Rayleigh_depol(ii)
         beta2  = ( GTOne - depol ) / ( 2.0_gtpk + depol ) 

         phasmoms_input = GTZero
         phasmoms_input(0, 1, 1) =  GTOne
         phasmoms_input(2, 1, 1) =  beta2
         phasmoms_input(2, 2, 1) =  six * beta2
         phasmoms_input(2, 5, 1) =  mrtsix * beta2
         phasmoms_input(1, 4, 1) =  3.D0 * ( GTOne - 2.D0*DEPOL ) / (GTOne + DEPOL)

!  set local nmoms for VLIDORT assignations (avoids needless copying if no aerosols)
         
         LM = local_nmoms

!  Start layer loop

         do n = 1, FO_nlayers

!  weightings

            AER_WT = Geophys%TotalODs%fa(n,ii) ; RAY_WT = Geophys%TotalODs%fr(n,ii)
 
!  Make FO Greekmat. No clouds here.....

            IF ( AERFLAG(n) ) THEN
              DO KK = 1, ngkmatc
                GK = gmask(kk) ; sk = dble(smask(kk)) ; ck = cmask(kk) ; rk = rmask(kk)
                RYCofs(0:2)     = RAY_WT * phasmoms_input(0:2, rk,1)
                A1Cofs(0:local_nmoms) = AER_WT * SK * Geophys%Aerosols%AEROSOL_SCATMOMS(ck,0:local_nmoms,ii)
                FO_Greekmat(0:2,n,gk)     = A1Cofs(0:2)    + RYCofs(0:2)
                FO_Greekmat(3:local_nmoms,n,gk) = A1Cofs(3:local_nmoms) 
              ENDDO
            ELSE
              DO KK = 1, ngkmatc
                GK = gmask(kk)  ; rk = rmask(kk) ; RYCofs(0:2) = RAY_WT * phasmoms_input(0:2, rk,1)
                FO_Greekmat(0:2,n,gk) =  RYCofs(0:2)
              ENDDO
            ENDIF

!  First Order Bulk optical properties (New Style).
!     Corrected Code 5/2/17, Takes proper care with Delta-M scaling 
 
            FO_omegas(n)  = Geophys%TotalODs%omega(n,ii)
            if ( FO_do_deltam .and. AERFLAG(n) ) then
              MOM2M = AER_WT * Geophys%Aerosols%AEROSOL_SCATMOMS(1,nstr2,ii)
              FO_truncfac(n) = MOM2M/dnm1
              fact1 = GTOne - FO_truncfac(n) * FO_omegas(n)
              FO_deltaus(n) = Geophys%TotalODs%taudp(n,ii) * fact1
            else
              FO_truncfac(n) = GTZero
              FO_deltaus(n)  = Geophys%TotalODs%taudp(n,ii)
            endif
            FO_extinction(n) = FO_deltaus(n) / FO_diffgrid(n)
 
!  Linearization of FO_Greekmat

            IF ( AERFLAG(n) ) THEN
              if ( LFvary(n) ) then
                do q = 1, LNvary(n)
                  L_Ray_wt = L_Geophys%L_TotalODs%L_fr(q,n,ii) 
                  L_Aer_wt = L_Geophys%L_TotalODs%L_fr(q,n,ii) 
                  DO kk = 1, ngkmatc
                    GK = gmask(kk) ; sk = dble(smask(kk)) ; ck = cmask(kk); rk = rmask(kk)
                    RYCofs(0:2)  = L_RAY_WT * phasmoms_input(0:2, rk,1)
                    A1Cofs(0:local_nmoms) = L_AER_WT * SK * Geophys%Aerosols%AEROSOL_SCATMOMS(ck,0:local_nmoms,ii)
                    FO_L_Greekmat(0:2,n,gk,q)            = A1Cofs(0:2) + RYCofs(0:2)
                    FO_L_Greekmat(3:local_nmoms,n,gk,q)  = A1Cofs(3:local_nmoms) 
                  ENDDO
                enddo
              endif
            endif

!  Linearized bulk properties. (everything pre-zeroed)
!     Corrected Code 5/2/17, Takes proper care with Delta-M scaling 

            if ( FO_do_deltam ) then
               fact1 = one - FO_truncfac(n) * FO_omegas(n)
               if ( LFvary(n) ) then
                  do q = 1, LNvary(n)
                    FO_L_truncfac(n,q) = FO_L_Greekmat(nstr2,n,1,q)/dnm1
                    L_fact1 = - FO_truncfac(n) * L_Geophys%L_TotalODs%L_omega(q,n,ii) - FO_L_truncfac(n,q) * FO_omegas(n)
                    FO_L_deltaus(n,q)    = L_Geophys%L_TotalODs%L_taudp(q,n,ii) * fact1 + Geophys%TotalODs%taudp(n,ii) * L_fact1
                    FO_L_omegas(n,q)     = L_Geophys%L_TotalODs%L_omega(q,n,ii)
                    FO_L_extinction(n,q) = FO_L_deltaus(n,q) / FO_diffgrid(n)
                  enddo
               endif
            else
               if ( LFvary(n) ) then
                  do q = 1, LNvary(n)
                    FO_L_deltaus(n,q)    = L_Geophys%L_TotalODs%L_taudp(q,n,ii)
                    FO_L_omegas(n,q)     = L_Geophys%L_TotalODs%L_omega(q,n,ii)
                    FO_L_extinction(n,q) = FO_L_deltaus(n,q) / FO_diffgrid(n)
                  enddo
               endif
            endif

!  End layer loop

         enddo

!  First Order albedos
         ALBEDO = Geophys%Surface%albedo(ii) ! Jbak 2019/05/03
         FO_reflec(1,1,1:FO_ngeoms) = ALBEDO
       
!   Solar Flux. Copy 2S value

         FO_FLUX  = 0.25_gtpk * FLUX_FACTOR  / FO_Pie

!  Phase Matrix calculation. Use the Local Greekmat array.
         if ( do_Jacobians ) then
           Call FO_VectorSS_PhasMat_Plus &
           ( FO_MaxGeometries, FO_maxlayers, FO_maxmoments, FO_maxatmoswfs, FO_do_profilewfs, & ! Dimensions/Flag
             FO_do_sunlight, FO_do_deltam, FO_Lvaryflags, FO_Lvarynums, FO_Lvarymoms,         & ! Flags/Linearization
             FO_nstokes, FO_ngeoms, FO_nlayers, local_nmoms, FO_aclevel,                      & ! Numbers
             FO_omegas, FO_truncfac, FO_Greekmat, FO_L_omegas, FO_L_truncfac, FO_L_Greekmat,  & ! optical
             FO_GenSpher, FO_Rotations, FO_ExactScat_up, FO_L_ExactScat_up )                    ! Functions/Output
         else
           Call FO_VectorSS_PhasMat &
           ( FO_MaxGeometries, FO_maxlayers, FO_maxmoments, FO_do_sunlight, FO_do_deltam,   & ! Dimensions/Flags
             FO_nstokes, FO_ngeoms, FO_nlayers, local_nmoms, FO_aclevel,                    & ! Numbers
             FO_omegas, FO_truncfac, FO_Greekmat,                                           & ! Optical input
             FO_GenSpher, FO_Rotations, FO_ExactScat_up )                                     ! Functions/Output
         endif

!do n = 1, 74
!   write(*,*)FO_ExactScat_up(n,1,1,1),FO_L_ExactScat_up(n,1,1,1,1)!,FO_omegas(n), FO_truncfac(n), FO_Greekmat(0:2,n,1)
!enddo
!stop'FO Exactcat fast'

!  time

         if (monitor_CPU) then
            call cpu_time(e3) ; Eofpctimes(13) = Eofpctimes(13) + e3 - e2  ! eofpc_FastOptime_FO
         endif

!  6b. Assign full-monochromatic (non-binned) properties, 2Stream
!      ----------------------------------------------------------

!  Monitoring

        if (monitor_CPU) call cpu_time(e2)

!  Preliminary zeroing

        ASYMM_INPUT  = zero ; D2S_SCALING = zero
        deltau_input = zero ; omega_input = zero

!  Main properties

         do n = 1, nlayers
            deltau_input(n) = Geophys%TotalODs%taudp(n,ii)
            omega_input(n)  = Geophys%TotalODs%omega(n,ii)
            if (omega_input(n) .gt. 0.999999d0) omega_input(n) = 0.999999d0 ! consistent setting across code
!            if (omega_input(n) .gt. 0.999999d0) omega_input(n) = 0.999999d0
            if (omega_input(n) .lt. 1.0d-6) omega_input(n) = 1.0d-6 ! consistent setting across code
!            if (omega_input(n) .lt. 0.000001d0) omega_input(n) = 0.000001d0
            if ( Aerflag(n) ) then
               ASYMM_INPUT(n) = FO_Greekmat(1,n,1) / 3.0_gtpk
               D2S_SCALING(n) = FO_Greekmat(2,n,1) / 5.0_gtpk
            else
               D2S_SCALING(n) = beta2 / 5.0_gtpk
            endif
         enddo

!  Linearized properties. No linearized moments.

         do n = 1, nlayers
           if ( LFvary(n) ) then
             do q = 1, LNvary(n)
               L_deltau_input(n,q) = L_Geophys%L_TotalODs%L_taudp(q,n,ii)
               L_omega_input(n,q)  = L_Geophys%L_TotalODs%L_omega(q,n,ii)
               if ( Aerflag(n) ) then
                 L_ASYMM_INPUT(n,q) = FO_L_Greekmat(1,n,1,q) / 3.0_gtpk
                 L_D2S_SCALING(n,q) = FO_L_Greekmat(2,n,1,q) / 5.0_gtpk
               endif
             enddo
           endif
         enddo

!  Albedo
         ALBEDO = Geophys%Surface%albedo(ii)
!         FO_reflec (1, 1, 1) =  ALBEDO
!   1/15/16. Solar Flux. Now set directly (Optional default to sun-normalized).

         if ( do_Sun_Normalized ) then
            FLUX_FACTOR = GTONE
         else
            FLUX_FACTOR = Geophys%SolarSpec%SunSpec(ii)
         endif

!  Time

         if (monitor_CPU) then
            call cpu_time(e3) ; Eofpctimes(14) = Eofpctimes(14) + e3 - e2  ! eofpc_FastOptime_2s 
         endif


!  6c. Full-wavelength FO calculation
!      ------------------------------

         if ( Monitor_CPU ) call cpu_time(e2)

!  Call

         if ( do_Jacobians ) then
           Call SSV_Integral_ILPS_UP_Optimized &
            ( FO_maxgeometries, FO_maxlayers, FO_maxfine, FO_maxatmoswfs, FO_maxsurfacewfs,              & ! Inputs (dimensioning)
              FO_do_planpar, FO_do_regular_ps, FO_do_enhanced_ps, FO_doNadir, FO_do_sleave,              & ! Inputs (Flags)
              FO_do_sunlight, FO_do_lambertian, FO_do_profilewfs, FO_do_reflecwfs, FO_do_sleavewfs,      & ! Inputs (control, Jac )
              FO_Lvaryflags, FO_Lvarynums, FO_n_reflecwfs, FO_n_sleavewfs,                               & ! Inputs (control, Jac )
              FO_nstokes, FO_ngeoms, FO_nlayers, FO_nfinedivs, FO_AcLevel,                               & ! Inputs (control output)
              FO_reflec, FO_slterm, FO_extinction, FO_deltaus, FO_exactscat_up, FO_flux, FO_fluxvec,     & ! Inputs (Optical)
              FO_LS_reflec, FO_LSSL_SLterm, FO_L_extinction, FO_L_deltaus, FO_L_exactscat_up,            & ! Inputs (Optical - Lin)
              FO_Mu0, FO_Mu1, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,                      & ! Inputs (Geometry)
              FO_Raycon, FO_cota, FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,        & ! Inputs (Geometry)
              FO_stokes_up, FO_stokes_db, FO_LP_Jacobians_up, FO_LP_Jacobians_db, FO_LS_Jacobians_db )     ! Output
         else
           Call SSV_Integral_I_UP_Optimized &
            ( FO_maxgeometries, FO_maxlayers, FO_maxfine,                                                    & ! Inputs (dimension)
              FO_do_planpar, FO_do_regular_ps, FO_do_enhanced_ps, FO_doNadir, FO_do_sleave,                  & ! Inputs (Flags)
              FO_do_sunlight, FO_do_lambertian, FO_nstokes, FO_ngeoms, FO_nlayers, FO_nfinedivs, FO_AcLevel, & ! Inputs (ctrl out)
              FO_reflec, FO_slterm, FO_extinction, FO_deltaus, FO_exactscat_up, FO_flux, FO_fluxvec,         & ! Inputs (Optical)
              FO_Mu0, FO_Mu1, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,                          & ! Inputs (Geometry)
              FO_Raycon, FO_cota, FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,            & ! Inputs (Geometry)
              FO_stokes_up, FO_stokes_db, FO_cumsource_up )                                                    ! Outputs
         endif


!  FO fast results

         if ( dir.eq.UPIDX ) then
           do v = 1, ngeoms
             Stokes_FO_Fast(v,1:ns) = FO_stokes_up(1:ns,v) + FO_stokes_db(1:ns,v)
           enddo
         endif
         if ( do_profile_WFs ) then
           do v = 1, ngeoms
             do n = 1, nlayers
               if ( LFvary(n) ) then
                 do q = 1, LNvary(n)
                   FOJac(1:ns)  = FO_LP_Jacobians_up(1:ns,v,n,q)
                   IF ( DIR.eq.UpIdx ) FOJac(1:ns) = FOJac(1:ns) + FO_LP_Jacobians_db(1:ns,v,n,q)
                   LP_Stokes_FO_Fast(q,n,v,1:ns) = FOJac(1:ns)
                 enddo
               endif
             enddo
           enddo
         endif
      
         if ( do_Surface_WFs ) then
           do v = 1, ngeoms
             q = 1 ; FOJac(1:ns)  = GTzero ; IF ( DIR.eq.UpIdx ) FOJac(1:ns) = FO_LS_Jacobians_db(1:ns,v,q)
             LS_Stokes_FO_Fast(q,v,1:ns) = FOJac(1:ns)
             !print * ,ii,  FOJAC(1) 
           enddo
         endif

!  Timing

         if (monitor_CPU) then
            call cpu_time(e3) ; Eofpctimes(15) = Eofpctimes(15) + e3 - e2  ! eofpc_FastRTM_FO
         endif

!  6e. Full-wavelength 2S calculation
!      ------------------------------

!  timing

        if (monitor_CPU) call cpu_time(e2)

!  Call
!  Rob Fix 8/3/18. New Calling statements for multiple geometries

        if ( do_Jacobians ) then
          CALL TWOSTREAM_LPS_MASTER &
          ( S2_MAXLAYERS, S2_MAXTOTAL, S2_MAXGEOMS, S2_MAXMESSAGES, S2_NLAYERS, S2_NTOTAL, S2_NGEOMS, & ! Inputs
            S2_MAX_ATMOSWFS, S2_MAX_SURFACEWFS, S2_MAX_SLEAVEWFS,             & ! Dimensions
            DO_PLANE_PARALLEL, DO_D2S_SCALING, TOPMOST_LEVEL,                 & ! Inputs
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,             & ! Inputs
            DO_PDINVERSE, BVPINDEX, BVPSCALEFACTOR, FLUX_FACTOR,              & ! Inputs
            DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs 
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,  & ! Inputs
            STREAMVAL, STREAMINV, GEOMETRIES, AUX_GEOMS, CHAPFACS,            & ! Inputs
            DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, ALBEDO,      & ! Inputs
            BRDF_F_0, BRDF_F, UBRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,          & ! Inputs
            L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,      & ! Inputs
            LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LSSL_ISOTROPIC, LSSL_F_0,     & ! Inputs
            INTENSITY_2S, PROFILEWFS_2S, SURFACEWFS_2S,                       & ! Outputs
            STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,            & ! Exception handling
            STATUS_EXECUTION, S2_MESSAGE, S2_TRACE_1, S2_TRACE_2 )              ! Exception handling
        else
          Call TWOSTREAM_MASTER &
          ( S2_MAXLAYERS, S2_MAXTOTAL, S2_MAXGEOMS, S2_MAXMESSAGES, S2_NLAYERS, S2_NTOTAL, S2_NGEOMS, & ! Inputs
            DO_PLANE_PARALLEL, DO_D2S_SCALING, TOPMOST_LEVEL,                 & ! Inputs
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,             & ! Inputs
            DO_PDINVERSE, BVPINDEX, BVPSCALEFACTOR, FLUX_FACTOR,              & ! Inputs
            STREAMVAL, STREAMINV, GEOMETRIES, AUX_GEOMS, CHAPFACS,            & ! Inputs
            DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, ALBEDO,      & ! Inputs
            BRDF_F_0, BRDF_F, UBRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,          & ! Inputs
            INTENSITY_2S, STATUS_INPUTCHECK, STATUS_EXECUTION,                & ! Outputs 
            C_NMESSAGES, C_MESSAGES, C_ACTIONS, S2_MESSAGE, S2_TRACE_1, S2_TRACE_2 )
        endif

!  Exception handling for Input checks !mick fix 2/14/2015 - adjusted C_messages indexing

         if ( STATUS_INPUTCHECK .eq. 1 ) then
            DO M = 1, C_Nmessages
               messages(Nmessages+2*M-1) = '(2STREAM_2p4 Message) '//Adjustl(TRIM(C_messages(M-1)))
               messages(Nmessages+2*M)   = '(2STREAM_2p4 action ) '//Adjustl(TRIM(C_ACTIONS(M-1)))
            ENDDO
            Nmessages = Nmessages + 2*C_Nmessages
            write(C4,'(I4)')ii ; messages(Nmessages+1) = '(Twostream Input Check, RT_FULL # = '//C4
            Nmessages = Nmessages + 1 ; Fail = .true.;  go to 70
         endif

!  Exception handling for Calculation

         IF ( STATUS_EXECUTION .eq. 1 ) THEN
            messages(Nmessages+1)   = '(2STREAM_2p4 Message) '//Adjustl(TRIM(S2_Message))
            messages(Nmessages+2)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(S2_Trace_1))
            messages(Nmessages+3)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(S2_Trace_2))
            Nmessages = Nmessages + 3
            write(C4,'(I4)')ii ; messages(Nmessages+1) = '(Twostream Execution, RT_FULL # = '//C4
            Nmessages = Nmessages + 1 ; Fail = .true. ; go to 70
         endif

!  2stream results (upwelling only)
!  Rob Fix 8/3/18. Allowed multiple geometries
         if ( dir.eq.UPIDX ) Intensity_2S_Fast(1:ngeoms)  = INTENSITY_2S(1:ngeoms)
         if ( do_profile_WFs ) then
           do v = 1, ngeoms
             do n = 1, nlayers
               if ( LFvary(n) ) then
                 do q = 1, LNvary(n)
                   LP_Intensity_2S_fast(q,n,v)  = PROFILEWFS_2S(v,n,q)

!  debug check on fast output. Compare with fort.5555 output in exact tool
!write(7777,*)n,q,v,ii,Geophys%WavGrids%wav(ii),LP_Intensity_2S_fast(q,n,v),LP_Stokes_FO_Fast(q,n,v,1)

                 enddo
               endif
             enddo
           enddo
         endif

         if ( do_surface_WFs ) then
           do v = 1, ngeoms
             q = 1 ; LS_Intensity_2S_fast(q,v)  = SURFACEWFS_2S(v,q)
            !print * , geophys%wavgrids%wav(ii),ii
           enddo
         endif

!  fast     results debugging


!  Timing

         if (monitor_CPU) then
            call cpu_time(e3) ; Eofpctimes(16) = Eofpctimes(16) + e3 - e2  ! eofpc_FastRTM_2S
         endif

!  6f. Backmapping to get PCA-Estimated intensity
!  ----------------------------------------------

 !  timing

        if (monitor_CPU) call cpu_time(e2)

!  Apply correction factors to Fast calculation. Stokes-Q 5/20/16, Stokes-U, 9/23/16. Merged, 12/27/17

        if ( do_3M_correction ) then
           do v = 1, ngeoms
             Ifast(1)= Intensity_2S_Fast(v) + Stokes_FO_Fast(v,1)
             Stokes_Eofpc(v,1,ii) = Ifast(1) * Stokes_Corrfacs(v,1,ii)
             if ( nstokes.gt.1 ) then
               Ifast(2:ns) = 2.0d0 * Stokes_FO_Fast(v,2:ns)
               Stokes_Eofpc(v,2:ns,ii)  = Ifast(2:ns) + Stokes_Corrfacs(v,2:ns,ii)
             endif
           enddo
         else
           do v = 1, ngeoms
             Stokes_Eofpc(v,1,ii) = Stokes_Corrfacs(v,1,ii)*Intensity_2S_Fast(v)
           enddo
         endif

!  Apply correction factors to Profile Jacobian fast calculations

         if ( do_3M_correction .and. do_Profile_WFs ) then
           do v = 1, ngeoms
             do n = 1, nlayers
               if ( LFvary(n) ) then
                 do q = 1, LNvary(n)
                   Jfast(1) = LP_Intensity_2S_fast(q,n,v) + LP_Stokes_FO_Fast(q,n,v,1)
                   Results_LP_Eofpc%LP_Jacobians_Eofpc(n,q,v,1,ii) = &
                                 Jfast(1) * Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,1,ii)
                   if ( nstokes.gt.1 ) then
                     Jfast(2:ns) = 2.0d0 * LP_Stokes_FO_Fast(q,n,v,2:ns)
                     Results_LP_Eofpc%LP_Jacobians_Eofpc(n,q,v,2:ns,ii) = &
                                 Jfast(2:ns) + Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,2:ns,ii)
!if ( n.eq.1)write(67,*)k,ii,v,n,q,Jfast(2:ns),Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,2:ns,ii), &
!Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,1,ii)
                   endif
                 enddo
               endif
             enddo
           enddo
         endif

!  Apply correction factors to Surface Jacobian fast calculations

         if ( do_3M_correction .and. do_Surface_Wfs ) then
           q = 1
           do v = 1, ngeoms
             Jfast(1) = LS_Intensity_2S_fast(q,v) + LS_Stokes_FO_Fast(q,v,1)
             Results_LS_Eofpc%LS_Jacobians_Eofpc(q,v,1,ii) = &
                                 Jfast(1) * Results_LS_Eofpc%LS_Jacobians_Corrfacs(q,v,1,ii)
             if ( nstokes.gt.1 ) then
               Jfast(2:ns) = 2.0d0 * LS_Stokes_FO_Fast(q,v,2:ns)
               Results_LS_Eofpc%LS_Jacobians_Eofpc(q,v,2:ns,ii) = &
                                 Jfast(2:ns) + Results_LS_Eofpc%LS_Jacobians_Corrfacs(q,v,2:ns,ii)
             endif
           enddo
         endif

!  debug check on fast output, to FORT.8820+v, FORT.9920+v
!   Rob Fix, MultiGoemetry 8/3/18

         if ( do_debug_output ) then
           do v = 1, ngeoms
             if ( nstokes.eq.1 ) then
               write(8820+v,251)Assigned_bins(ii),i,ii,Geophys%WavGrids%wav(ii),Intensity_2S_Fast(v),&
                       Stokes_FO_Fast(v,1)
             else
               write(9920+v,251)Assigned_bins(ii),i,ii,Geophys%WavGrids%wav(ii),Intensity_2S_Fast(v),&
                       Stokes_FO_Fast(v,1:nstokes)
             endif
           enddo
251        format(i2,2i4,1p5e20.10)
         endif

!  Set DOLP to zero if no polarization (Nstokes = 1), and skip DOLP calculation
!         degree of polarization (only for NSTOKES = 3 or 4)
!mick fix 7/30/2018 - changed array dims below from using "w" to "ii"

         IF ( NSTOKES.gt.1 ) THEN
           do v = 1, ngeoms
             !SQ2U2 = SQRT ( Stokes_Eofpc(v,2,w) * Stokes_Eofpc(v,2,w) &
             !             + Stokes_Eofpc(v,3,w) * Stokes_Eofpc(v,3,w) )
             !DOLP_Eofpc(V,W) = SQ2U2 / Stokes_Eofpc(v,1,w)
             !if ( nstokes .eq. 4 ) then
             !  SV2 = SQRT ( Stokes_Eofpc(v,4,w) * Stokes_Eofpc(v,4,w) )
             !  DOCP_Eofpc(V,W) = SV2 / Stokes_Eofpc(v,1,w)
             !endif
             SQ2U2 = SQRT ( Stokes_Eofpc(v,2,ii) * Stokes_Eofpc(v,2,ii) &
                          + Stokes_Eofpc(v,3,ii) * Stokes_Eofpc(v,3,ii) )
             DOLP_Eofpc(V,ii) = SQ2U2 / Stokes_Eofpc(v,1,ii)
             if ( nstokes .eq. 4 ) then
               SV2 = SQRT ( Stokes_Eofpc(v,4,ii) * Stokes_Eofpc(v,4,ii) )
               DOCP_Eofpc(V,ii) = SV2 / Stokes_Eofpc(v,1,ii)
             endif
           ENDDO
         ENDIF

!  Profile Jacobians results

         if ( do_profile_WFs ) then
          do v = 1, ngeoms
            do n = 1, nlayers
              if ( LFvary(n) ) then
                do q = 1, LNvary(n)
                  Rad1     = LP_Stokes_FO_Fast(q,n,v,1) + LP_Intensity_2S_Fast(q,n,v)
                  Results_LP_Eofpc%LP_Jacobians_Eofpc(n,q,v,1,ii) = &
                              Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,1,ii) * Rad1
!      IF(q ==1) print * ,ii, Results_LP_Eofpc%LP_jacobians_corrfacs(n, q, v, 1,ii), lp_stokes_fo_fast(q, n, v, 1) , lp_intensity_2s_fast(q, n, v)
                  if ( nstokes.gt.1 ) then
                    do s = 2, nstokes
                      Rad1  = 2.0d0 * LP_Stokes_FO_Fast(q,n,v,s)
                      Results_LP_Eofpc%LP_Jacobians_Eofpc(n,q,v,s,ii) = &
                              Results_LP_Eofpc%LP_Jacobians_Corrfacs(n,q,v,s,ii) + Rad1
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        endif
	 
!  Surface Jacobian results

        if ( do_Surface_WFs ) then
          do v = 1, ngeoms
            q = 1 ; Rad(1)     = LS_Stokes_FO_Fast(q,v,1) + LS_Intensity_2S_Fast(q,v)
            Results_LS_Eofpc%LS_Jacobians_Eofpc(q,v,1,ii) = &
                            Results_LS_Eofpc%LS_Jacobians_Corrfacs(q,v,1,ii)*Rad(1)
           
            if ( nstokes.gt.1 ) then
              Rad(2:ns)  = 2.0d0 * LS_Stokes_FO_Fast(q,v,2:ns)
              Results_LS_Eofpc%LS_Jacobians_Eofpc(q,v,2:ns,ii) = &
                            Results_LS_Eofpc%LS_Jacobians_Corrfacs(q,v,2:ns,ii) + Rad(2:ns)
            endif
          enddo
        endif
!  Timing

         if (monitor_CPU) then
            call cpu_time(e3) ; Eofpctimes(17) = Eofpctimes(17) + e3 - e2  ! Eofpctime (Final correction)
         endif
 
!  End points-loop within one bin
       
      enddo

!   END PARALLEL REGION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Continuation point for error handling (Fast calculation)
!  If Fail, Return

70    continue
      if ( Fail )  return

!  Continuation point for avoiding Fast Full-2S calculation
!    ( Continuation point at end of Fourth Bin loop )

556   continue

!      write(*,*)' -- Done Fast RT calculations for all points --'

!  Go to end of Bin loop
!      goto 633
!  Now do the special case if flagged
!  Full Lidort and FO calculations (All wavelengths)
!        PLACEHOLDER
!  Continuation point for avoiding fast calculation
!633   continue

!  End fourth and final bin loop

   enddo

!  print * , Results_LP_Eofpc%LP_Jacobians_Eofpc(1,1,1,1,80:90)
!          print  *  , '------------'
!               print * , RO%LS_Jacobians(1,1,1,80:90)
!mick fix 4/6/2015 - moved from down below and modified
!  Time spent in equivalent OpenMP parallel region #2 and parallel regions overall

   if (monitor_CPU) then
      call cpu_time(ser_e2) ; SerEofpcTime = SerEofpcTime + (ser_e2 - ser_e1)
   endif

!  Time spent in equivalent OpenMP parallel region
!    Have to subtract off time spent in 2S+FO setup, FO Exec, and 2S Exec

   if (monitor_CPU) then
      SerEofpcTime = sum(EofpcTimes(1:12)) + sum(EofpcTimes(13:17))
      EofpcTimes(18) = SerEofpcTime
   endif

!  Final CPU call

   if (monitor_CPU) then
      call cpu_time(e3) ; EofpcRuntime = EofpcRuntime + e3-e1
   endif

!  Output Variables
   IF (allocated (RO%stokes) ) THEN 
      deallocate (RO%stokes, RO%PJAC, RO%SJAC)
   ENDIF
   allocate (RO%STOKES(4, nwav))
   allocate (RO%PJAC (nlayers,n_totalatmos_wfs, 4, nwav))
   allocate (RO%SJAC (n_surface_wfs,  4, nwav))
   RO%Stokes(1:ns, 1:nwav)   = Stokes_Eofpc(1, 1:ns, 1:nwav)
   RO%PJAC(1:nlayers, 1:n_totalatmos_wfs, 1:4, 1:nwav) = & 
           Results_LP_EofPc%LP_Jacobians_EofPc(1:nlayers, 1:n_totalatmos_wfs, 1, :, 1:nwav)
   RO%SJAC(1:n_surface_wfs,1:4, 1:nwav)  = Results_LS_EofPc%LS_Jacobians_EofPc(1:n_surface_wfs, 1, :, 1:nwav)
!  deallocating local variables
   deallocate ( Stokes_Eofpc, Stokes_Corrfacs)
   IF (ns > 1) THEN 
    deallocate ( DOLP_Eofpc, DOCP_Eofpc)
   ENDIF
   IF (do_Jacobians) THEN
   !  deallocate (LP_Jacobians_Eofpc, LP_Jacobians_Corrfacs)
     !deallocate (LS_Jacobians_Eofpc, LS_Jacobians_Corrfacs)
   ENDIF
   first = .false.
!  normal Finish
   RETURN
!  Bookkeeping

END SUBROUTINE VLIDORT_Eofpc_Master

!  End module

end module VLIDORT_eofpc_module

