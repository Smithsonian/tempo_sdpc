module VLIDORT_exact_module

!  History.

!  V1: used VLIDORT for the Fo code. Too slow.

!  V2: used VLIDORT_NOFO and Stand-alone Optimized/Leveraged FO code
!  V2: used Optimized 2stream code.    

!  Rob Fix 8/3/18. MultiGeometry capability properly installed
!                  Use of new 2S code.
!   USE LIDORT_DATA_MODULE
!  Modules in GEMSTOOL_sourcecode/structures
   USE ozprof_data_module, ONLY: num_iter
   USE GEMSTOOL_pars_m
   USE GEMSTOOL_Input_types_m
   USE GEMSTOOL_Geophys_types_m
   USE GEMSTOOL_Result_types_m

   USE GEMSTOOL_L_Geophys_types_m
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

!  Initializer routine for VLIDORT

   USE GEMSTOOL_RTinitialize_m
   use vlidort_data_module, ONLY: & 
       Geophys, L_Geophys, Inputs, & 
       layer_vary_flag_cc, layer_vary_number_cc, profilewf_names_cc, &
       RO, RI
   
   INTEGER, PARAMETER :: dir = 1
   public

contains

   subroutine VLIDORT_Exact_Master &
    ( Monitor_CPU, do_VLIDORT_initialize, do_Jacobians,        & ! Flags 
      fail, NMessages, Messages )      ! Status 

   implicit none 
!  Inputs
!  ======

!  Monitoring flag

   logical, intent(in) :: Monitor_CPU

!  VLIDORT initialize flag

   logical, intent(inout) :: DO_VLidort_initialize

!  Jacobians flag

   logical, intent(in) :: do_Jacobians

!  Timing and status
!  =================

!  Exception handling variables

   logical, intent(inout)       :: fail
   integer, intent(inout)       :: NMessages
   character*(*), intent(inout) :: Messages (GT_maxMessages)

!  Local variables
!  ===============

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      VLIDORT  ARGUMENTS (V2.7 code)
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), SAVE       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), SAVE    :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut), SAVE          :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn

!  VLIDORT lnearized supplements itt structure

      TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup

!  VLIDORT linearized output structure

      TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut

!  Stokes-vector and Jacobian output results

   TYPE(GEMSTOOL_Exact_IC)       :: Results_Exact
   TYPE(GEMSTOOL_LP_Exact_IC)    :: Results_LP_Exact
   TYPE(GEMSTOOL_LS_Exact_IC)    :: Results_LS_Exact
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

!  
      LOGICAL              :: DO_RAYLEIGH_ONLY 
 
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
!      INTEGER        :: N_SURFACE_WFS
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
   real(gtpk) :: FO_L_exactscat_up ( FO_MAXLAYERS,4,4,FO_MaxGeometries, FO_maxatmoswfs)
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

!   Proxy GEMSTOOL variables
!   ========================

!  Aerosols/Clouds

   logical :: do_Aerosols, do_Clouds

!  Enhanced sphericity flag, Optional 2S flag.  Controlled by the Config-file input "do_fast_calculation"

   logical :: do_enhanced_ps
   logical :: do_optional_2stream

!    --> Number of LIDORT discrete ordinates, number of stokes components (3=Vector, 1 = Scalar)
!    --> Number of Geometries

   integer :: ngeoms, nstreams, nstokes, nlayers, nlevels, nwav, nstr2

!  Serial/Timing Variables
!  =======================

!  Serial tests (timing)

   REAL  :: ser_e1, ser_e2
   real  :: e1,e2,e3 
   real  :: SerExactTime

!  Other variables
!  ===============

!  timing variables

   real :: Exacttimes(20), ExactRunTime

!  OD help
   INTEGER, PARAMETER :: maxscatter=3, maxgkmatc=8, maxgksec=6
   logical, DIMENSION(GT_maxlayers)  :: Aerflag, CldFlag
   real(kind=GTPK) :: Aer_WT, Ray_WT, CLD_WT, L_Aer_WT, L_Ray_WT, L_Cld_wt, &
                      LD_taudp, LD_omega, fac, beta2, dnm1, L_fact1
   real(kind=GTPK) ::  SK, SQ2U2, SSRad(4), SV2, DEPOL, Int, jac, mom, momv, fact1, extinc, six, mrtsix
   real(kind=GTPK) :: A1Cofs(0:GT_maxmoments), A2Cofs(0:GT_maxmoments), RYCofs(0:2)
   real(kind=GTPK) :: MSRad(4), MSJac(4), FORad(4), FOJac(4), AllRad(4)
   real(kind=GTPK), DIMENSION(0:GT_maxmoments, 1:maxgksec, maxscatter) :: phasmoms_input
   real(kind=GTPK), DIMENSION(0:GT_maxmoments, 1:maxgksec) :: phasmoms_total_input
   real(kind=GTPK), DIMENSION(0:GT_maxmoments, 1:maxgksec) :: l_phasmoms_total_input
   real(kind=GTPK), DIMENSION(gt_maxatmoswfs, 0:GT_maxmoments, 1:MAXSTOKES_SQ) :: l_greekmat_total_input
   CHARACTER (LEN=31), DIMENSION(gt_maxatmoswfs) :: profilewf_names_cc

!  Help variables
   logical        :: Fail1, DO_Vlidort_Inpdebug, LFVary(GT_maxlayers)
   integer        :: q,l,n,s,nf,ns,w,m,v
   integer        :: LM, local_nmoms, nscatter, ngkmatc, ncoeffs , &
                     n_totalatmos_wfs, n_surface_wfs
   INTEGER, DIMENSION(GT_maxlayers) :: LNVary
   INTEGER        :: K, GK, RK, CK
   INTEGER, DIMENSION (maxgkmatc) ::  CMASK, GMASK, SMASK, RMASK
   character*5    :: c5
   character*100  :: Message1

!  Debug output. FORT 8810.
   logical, parameter :: do_debug_input     = .false.
   logical, parameter :: do_debug_output     = .false.
   logical, parameter :: do_print = .false.
   logical, parameter :: do_Sun_Normalized = .true.
   logical :: first=.true.
!  Mask Data
!  =========
   GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)     ! Greek Matrix indices
   RMASK = (/  1, 5, 5, 2, 3,  6,  6,  4 /)         ! Convention for Rayleigh

!   CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)         ! Convention for TEMPO aerosol output
!   SMASK = (/  1, 1, 1, 1, 1, -1, 1, 1 /)        ! Sign Mask  for TEMPO aerosol output

   CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)        ! Convention for RTS MIE aerosol output
   SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)     ! Sign Mask  for RTS MIE aerosol output

!  START CODE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   !write(*,*) '    ---- Doing PCA Exact-Plus RT (with Jacobians)'

!  Initialize time
!   ExactTimes   = 0.0 ! pre-initialized now, 2/18/16
!   ExactRunTime = 0.0 ! pre-initialized now, 2/18/16

    SerExactTime = 0.0

!  Initialize output - no longer needed, pre-initialized now, 2/18/16
!   Intensity_Exact     = zero
!   Intensity_LD_Exact  = zero
!   Intensity_FO_Exact  = zero

!  Initialize Exception handling

   fail = .false.
   NMessages = 0
   Messages  = ' '

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
   nwav        = Geophys%WavGrids%nwav
   nlayers     = Geophys%Atmos%nlayers ; nlevels = nlayers + 1
   local_nmoms = Geophys%Atmos%ngreek_moments_input
   !print * , do_aerosols, do_clouds, do_optional_2stream, nwav, nlayers,local_nmoms, n_totalatmos_wfs, n_surface_wfs
   IF (first .and. nwav > 50) THEN 
      WRITE(*,'(A,i3)') 'N call in exact:', nwav
       first = .false.
   ENDIF
!  constants

   six    = 6.0_gtpk
   mrtsix = -sqrt(six)

! Optical Depth Properties
   n_totalatmos_wfs = RI%n_totalatmos_wfs
   n_surface_wfs =    RI%n_surface_wfs
   LFvary (1:nlayers) = layer_vary_flag_cc(1:nlayers)
   LNvary (1:nlayers) = layer_vary_number_cc(1:nlayers)
   !profilewf_names_cc(1:n_totalatmos_wfs) = RI%profilewf_names_cc(1:n_totalatmos_wfs)

   !Geophys%totalods%taudp(1:nlayers, 1:nwav)=RI%taudp(1:nlayers, 1:nwav)
   !Geophys%totalods%taug(1:nlayers, 1:nwav)=RI%taug(1:nlayers, 1:nwav)
   !Geophys%totalods%omega(1:nlayers,1:nwav)=RI%omega(1:nlayers, 1:nwav)
   !Geophys%TotalODs%fr(1:nlayers, 1:nwav) = RI%fr(1:nlayers,1:nwav)
   !Geophys%TotalODs%fa(1:nlayers, 1:nwav) = RI%fa(1:nlayers,1:nwav)
   !L_Geophys%L_TotalODs%l_taudp(1:n_totalatmos_wfs, 1:nlayers,1:nwav) = &
   !                   RI%l_taudp(1:n_totalatmos_wfs,1:nlayers,1:nwav)
   !L_Geophys%L_TotalODs%l_omega(1:n_totalatmos_wfs, 1:nlayers,1:nwav) = &
   !                   RI%l_omega(1:n_totalatmos_wfs,1:nlayers,1:nwav)

!  VLIDORT. Initialize the Input structures.
   IF (first) THEN 
   call GEMSTOOL_VLIDORT_Initialize &
       ( Inputs, nlayers, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, fail1, Message1 )   
   ENDIF
    VLIDORT_FixIn%Cont%TS_NLAYERS          = nlayers
   IF ( VLIDORT_FixIn%Bool%TS_do_dnwelling ) then
      VLIDORT_ModIn%MUserVal%TS_user_levels(1) = Real(nlayers,fpk)
   ENDIF
    
   IF (num_iter == 0) THEN
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
   ENDIF
   !WRITE(*,*) 'Initialize', nwav, nlayers

   VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT    = local_nmoms
!  Fix heights

   VLIDORT_FixIn%Chapman%TS_height_grid(0:nlayers) = Geophys%Atmos%Level_heights(0:nlayers)
   VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT   = Geophys%Atmos%Level_heights(nlayers)
   
!  linearization control
!mick fix 7/30/2018 - initialize DO_ATMOS_LBBF, DO_SURFACE_LBBF, COLUMNWF_NAMES, & PROFILEWF_NAMES

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
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                 = n_surface_wfs
      IF (n_surface_wfs > 0) VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION     = .true.
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION     = .true.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION       = .true.
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION             = .true.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:nlayers)    =  LFVary(1:nlayers)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:nlayers)  =  LNVary(1:nlayers)
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY           = .false.
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS            =  n_totalatmos_wfs
   endif
   do_rayleigh_only = VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
!  Optical property zeroing
!  ------------------------

!  Entries in Greekmat

   if ( nstokes .eq. 1 ) then
      ngkmatc = 1 ; ncoeffs = 1
   else if ( nstokes .eq. 3 ) then
      ngkmatc = 5 ; ncoeffs = 6
   else if ( nstokes .eq. 4 ) then
      ngkmatc = 8 ; ncoeffs = 6
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

   if ( do_aerosols ) then
     AERFLAG(1:nlayers) = Geophys%Aerosols%AEROSOL_LAYERFLAGS(1:nlayers)
   else
     AERFLAG(1:nlayers) = .false.
   endif
   if ( do_clouds ) then
     !CLDFLAG(1:nlayers) = Geophys%Clouds%CLOUD_LAYERFLAGS(N)
   else
     CLDFLAG(1:nlayers) = .false.
   endif

!  Exception handling

   if ( fail1 ) then
      NMessages = 1 ; fail = .true.
      Messages(1) = Trim(Message1) ; return
   endif

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

!  Monitor for the above settings for VLIDORT and FO

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Exacttimes(1) = e3 - e2  ! setuptime VL and FO initialize.
   endif

!  1c. 2stream settings (Optional). This is not CPU-monitored
!  ==========================================================

!  Skip if not flagged
!    do_optional_2stream = .false.
   if ( do_optional_2stream ) then

!  Set BVP Scale factor and other control

      BVPSCALEFACTOR = GTOne

!  Set ntotal and nlayers, and ngeoms (new,8/3/18)

      S2_NLAYERS = nlayers
      S2_NTOTAL  = 2 * nlayers
      S2_NGEOMS  = ngeoms

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
      DO_SLEAVE_WFS  = VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS
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
      DO_SL_ISOTROPIC     = .false.  ! Until we enable it

!  Deltam-scaling flag MUST BE SET, regardless of the VLIDORT setting
!    Otherwise, results are unphysical. 2/5/15

      DO_D2S_SCALING      = .true.
      DO_PLANE_PARALLEL   = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL

!  Earth radius and heights (Copy FO values)
!      HEIGHT_GRID(0:nlayers)  = FO_heights(0:nlayers) 
!      EARTH_RADIUS            = FO_eradius

!  5/6/15 upgrade - initialize Surface BRDF inputs

      BRDF_F_0 = GTZero ; LS_BRDF_F_0 = GTZero
      BRDF_F   = GTZero ; LS_BRDF_F   = GTZero
      UBRDF_F  = GTZero ; LS_UBRDF_F  = GTZero

!  SLEAVE Terms

      SLTERM_ISOTROPIC    = GTZero ; LSSL_ISOTROPIC = GTZero
      SLTERM_F_0          = GTZero ; LSSL_F_0       = GTZero

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

!  Exception handling
               
   if ( FO_fail ) then
      messages(nmessages+1) = Adjustl(Trim(FO_message))
      messages(nmessages+1) = Adjustl(Trim(FO_trace))
      nmessages = nmessages+2
      go to 69
   endif

!  1e. First-Order Spherical Function calculations
!  ===============================================

!  Call optimized spherical functions routine

   Call FO_VectorSS_spherfuncs_Optimized &
        ( FO_STARTER, FO_Maxmoments, FO_MaxGeometries, FO_Do_Sunlight, & ! Inputs
          FO_nmoments, FO_ngeoms, FO_nstokes, FO_dtr, FO_vsign,        & ! Inputs
          FO_theta_boa, FO_alpha_boa, FO_phi_boa, FO_cosscat_up,       & ! Inputs
          FO_ROTATIONS, FO_GSHELP, FO_GENSPHER )                         ! Outputs &

!  Monitor for the above settings for FO Geometry and Spherical functions

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Exacttimes(2) = e3 - e2  ! setuptime FO Geometry/SpherFuncs
   endif

!  1f. 2S preliminary geometry calculations (Optional)
!  ===================================================

!  This is not monitored for CPU

   IF ( do_Optional_2Stream ) then

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

!      MU0_2S         = GEOMETRIES(1,2)
!      USER_STREAM_2S = GEOMETRIES(2,2)
!      CALL TWOSTREAM_AUXGEOM_PREPARE  ( GEOMETRIES(1,2), GEOMETRIES(2,2), STREAMVAL, AUX_GEOMS )
      CALL TWOSTREAM_AUXGEOM_PREPARE &
         ( S2_MAXGEOMS, S2_NGEOMS, GEOMETRIES(1:S2_MAXGEOMS,1,2), GEOMETRIES(1:S2_MAXGEOMS,2,2), STREAMVAL, AUX_GEOMS )

   ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  2. MAIN LOOP OVER Full WAVELENGTHS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!mick mod 4/6/2015 - added extra timing for better comparison with OpenMP
!                    version of code

    if ( Monitor_CPU ) call cpu_time(ser_e1)

!######################################## START PARALLEL REGION ##########

!  Start "point" loop
!   write(*,'(10(a, i5))') 'nwav = ',nwav, 'Nlayer=',nlayers, 'Nstreams=',nstreams,'nstokes=', nstokes
!  Start wavelength Loop
   do w = 1, nwav

!  2a. Compute input optical properties for VLIDORT 
!  ------------------------------------------------

!  "Problem Ray"  (Same for each layer)

      depol = Geophys%Xsecs%Rayleigh_depol(w)
      beta2  = ( GTOne - depol ) / ( 2.0_gtpk + depol ) 

      phasmoms_input = GTZero
      phasmoms_input(0, 1, 1) =  GTOne
      phasmoms_input(2, 1, 1) =  beta2
      phasmoms_input(2, 2, 1) =  six * beta2
      phasmoms_input(2, 5, 1) =  mrtsix * beta2
      phasmoms_input(1, 4, 1) =  3.D0 * ( GTOne - 2.D0*DEPOL ) / (GTOne + DEPOL )
      
!  set local nmoms for VLIDORT assignations (avoids needless copying if no aerosols)

      LM = local_nmoms

!  Albedo for VLIDORT

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = Geophys%Surface%albedo(w)

!  VLIDORT Solar Flux. Now set directly (Optional default to sun-normalized).

      if ( do_Sun_Normalized ) then
         VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = GTOne
      else
         VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = Geophys%SolarSpec%SunSpec(w)
      endif

!  Start layer loop
      do n = 1, nlayers


!  Omega Toggling (Already done)

        LD_taudp = Geophys%TotalODs%taudp(n,w)
        LD_omega = Geophys%TotalODs%omega(n,w)
        if ( LD_omega .gt. 0.999999_gtpk) LD_omega = 0.999999_gtpk ! consistent setting across code
        if ( LD_omega .lt. 1.e-6_gtpk)    LD_omega = 1.e-6_gtpk    ! consistent setting across code

!  VLIDORT bulks

        VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n)    = LD_taudp
        VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)   = LD_omega

! Aerosol Phase function
        nscatter = 1
! Cloud Phase Function
  
!  VLIDORT  Linearized properties. 
!   NEEDS WORK FOR INPUTS with AEROSOLS.....
       if ( LFvary(n) ) then
          do q = 1, LNvary(n)
            VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(q,n) = L_Geophys%L_TotalODs%L_taudp(q,n,w) / LD_taudp
            VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(q,n) = L_Geophys%L_TotalODs%L_omega(q,n,w) / LD_omega
          enddo
       endif

!  weightings
        CLD_WT = Geophys%TotalODs%fc(n,w) ; AER_WT = Geophys%TotalODs%fa(n,w) ; RAY_WT = Geophys%TotalODs%fr(n,w) 

!  Phase function
        ! sum up phase moments as required in LIDORT
        phasmoms_total_input(0:LM,1:ncoeffs) = phasmoms_input(0:LM,1:ncoeffs, 1)*RAY_WT
        IF (AER_WT /= 0.0 ) phasmoms_total_input(0:LM, 1:ncoeffs) = & 
           phasmoms_input(0:LM,1:ncoeffs,2)*AER_WT+phasmoms_total_input(0:LM, 1:ncoeffs)
        IF (CLD_WT /= 0.0 ) phasmoms_total_input(0:LM, 1:ncoeffs) = & 
           phasmoms_input(0:LM,1:ncoeffs,3)*CLD_WT+phasmoms_total_input(0:LM, 1:ncoeffs)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:LM,n,gmask(1:ngkmatc)) = & 
        phasmoms_total_input(0:LM,rmask(1:ngkmatc))
        
!  Linearization 
        
      l_greekmat_total_input(:, : , :) = GTZERO
       do q = 1, n_totalatmos_wfs
          IF ( profilewf_names_cc(q) == 'rayleigh optical thickness-----' ) THEN
           DO l = 0, local_nmoms
             DO k = 1, ncoeffs
                IF (phasmoms_total_input(l, k) /= 0.0) THEN
                   l_phasmoms_total_input(l, k) = & 
                   ( phasmoms_input(l, k, 1) - phasmoms_total_input(l, k) ) & 
                   *RAY_WT/phasmoms_total_input(l,k)
                ELSE
                   l_phasmoms_total_input(l, k) = 0.0
                ENDIF
             ENDDO
           ENDDO
           l_greekmat_total_input(q, 0:LM, gmask(1:ngkmatc)) = &
                                          l_phasmoms_total_input(0:LM, rmask(1:ngkmatc))
           IF ( ngkmatc > 1 )  l_greekmat_total_input(q, 0:LM, 15) &
                                  = - l_greekmat_total_input(q, 0:LM,15)
           ELSE IF ( profilewf_names_cc(q) == 'rayleigh scattering coefficient' ) THEN
           !  Still need to consider the variation in phase function
           l_greekmat_total_input(q, : , :) = GTZERO
           !  w.r.t aerosol extinction coefficient / aerosol optical thickness
           !  aerosol scattering albedo does not change
           !  xliu: April 13, 2007 (consider the variation in phase function)         
           ENDIF
       enddo
       VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:n_totalatmos_wfs,0:LM, n,:) = & 
           l_greekmat_total_input(1:n_totalatmos_wfs, 0:LM, :) 

!  Phase function normalization

        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
        
!  End layer loop

      enddo

!  time

      if (monitor_CPU) then
        call cpu_time(e3) ; Exacttimes(3) = Exacttimes(3) + e3 - e2  !  ExactOpTime_VLIDORT
      endif

!  2b. Compute layer input optical properties for FO Model 
!  -------------------------------------------------------

!  initialize CPU monitoring

      if ( Monitor_CPU ) call cpu_time(e2)

!  First Order albedos

      FO_reflec(1,1,1:FO_ngeoms) = Geophys%Surface%albedo(w)
      if ( FO_do_reflecwfs ) FO_LS_reflec(1,1,1,1:FO_ngeoms) = one

!  First Order optical properties (New Style). Use VLIDORT input where possible
!     Corrected Code 5/2/17, Takes proper care with Delta-M scaling 
 
      do n = 1, FO_nlayers

!  bulk properties and trauncation factors.

        FO_omegas(n)  = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
        if ( FO_do_deltam ) then
          FO_truncfac(n) = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(nstr2,n,1)/dnm1
          fact1 = GTOne - FO_truncfac(n) * FO_omegas(n)
          FO_deltaus(n) = Geophys%TotalODs%taudp(n,w) * fact1
        else
          FO_truncfac(n) = GTZero
          FO_deltaus(n)  = Geophys%TotalODs%taudp(n,w)
        endif
        FO_extinction(n) = FO_deltaus(n) / FO_diffgrid(n)

!  Greekmat and its linearization (not required for Rayleigh-only)

        FO_Greekmat(0:local_nmoms,n,:) = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:local_nmoms,n,:)
        if ( LFvary(n).and. (AERFLAG(n).or.cldflag(n)) ) then
          do q = 1, LNvary(n)
            Do L = 0, local_nmoms
              do gk  = 1, 16
                FO_L_Greekmat(L,n,gk,q) = FO_Greekmat(L,n,gk) * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,L,n,gk)
              enddo
            enddo
          enddo
        endif

!  End layer loop

      enddo

!  Linearized bulk properties. (everything pre-zeroed)
!     Corrected Code 5/2/17, Takes proper care with Delta-M scaling 

      do n = 1, FO_nlayers
         if ( FO_do_deltam ) then
            fact1 = one - FO_truncfac(n) * FO_omegas(n)
            if ( LFvary(n) ) then
               do q = 1, LNvary(n)
                  FO_L_truncfac(n,q) = FO_L_Greekmat(nstr2,n,1,q)/dnm1
                  L_fact1 = - FO_truncfac(n) * L_Geophys%L_TotalODs%L_omega(q,n,w) - FO_L_truncfac(n,q) * FO_omegas(n)
                  FO_L_deltaus(n,q)    = L_Geophys%L_TotalODs%L_taudp(q,n,w) * fact1 + Geophys%TotalODs%taudp(n,w) * L_fact1
                  FO_L_omegas(n,q)     = L_Geophys%L_TotalODs%L_omega(q,n,w)
                  FO_L_extinction(n,q) = FO_L_deltaus(n,q) / FO_diffgrid(n)
               enddo
            endif
         else
            if ( LFvary(n) ) then
               do q = 1, LNvary(n)
                  FO_L_deltaus(n,q)    = L_Geophys%L_TotalODs%L_taudp(q,n,w)
                  FO_L_omegas(n,q)     = L_Geophys%L_TotalODs%L_omega(q,n,w)
                  FO_L_extinction(n,q) = FO_L_deltaus(n,q) / FO_diffgrid(n)
               enddo
            endif
         endif
      enddo

!   1/20/16. Solar Flux. Now set directly (Optional default to sun-normalized).
!                        Sufficient to copy VLIDORT input

      FO_FLUX  = 0.25_gtpk * VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR  / FO_Pie

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

!do n = 1, FO_nlayers
!   write(*,*)FO_omegas(n), FO_ExactScat_up(n,1,1,1),FO_L_ExactScat_up(n,1,1,1,1)!,FO_truncfac(n), FO_Greekmat(0:2,n,1),&
!enddo
!stop'FO Exactcat'

!  time

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(4) = Exacttimes(4) + e3 - e2  !  ExactOpTime_FO
      endif

!  2c. Compute Optional 2stream optical properties
!  -----------------------------------------------

!  2-stream optical properties (Copy VLIDORT where possible)

      if ( DO_OPTIONAL_2STREAM ) then

!  CPU monitoring

        if (monitor_CPU) call cpu_time(e2)

!  Layer optical depths

        do n = 1, nlayers
          deltau_input(n) = Geophys%TotalODs%taudp(n,w)
          omega_input(n)  = Geophys%TotalODs%omega(n,w)
          ASYMM_INPUT(n)  = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(1,n,1)/3.d0
          D2S_SCALING(n)  = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(2,n,1)/5.d0
        enddo

!  Linearized properties

        do n = 1, nlayers
          if ( LFvary(n) ) then
            do q = 1, LNvary(n)
               L_deltau_input(n,q) = deltau_input(n) * VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(q,n) 
               L_omega_input(n,q)  = omega_input(n)  * VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(q,n)
!     if (w ==1)    print * , omega_input(n),vlidort_linfixin%optical%TS_l_omega_total_input(1, n)
               L_ASYMM_INPUT(n,q)  = ASYMM_INPUT(n)  * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,1,n,1)
               L_D2S_SCALING(n,q)  = D2S_SCALING(n)  * VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,2,n,1)
            enddo
          endif
        enddo

!  surface and solar

        ALBEDO      = VLIDORT_FixIn%Optical%TS_lambertian_albedo
        FLUX_FACTOR = VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  time

        if (monitor_CPU) then
           call cpu_time(e3) ; Exacttimes(5) = Exacttimes(5) + e3 - e2  !  ExactOpTime_2S
        endif

!  end optional 2S clause

      endif

!  2d. VLIDORT call
!  ----------------

!  Time

      if (monitor_CPU) call cpu_time(e2)

!  debug flag

      DO_Vlidort_Inpdebug = .false. ; IF ( w.eq.1 ) DO_Vlidort_Inpdebug = .true.

!  Call 
      If ( do_Jacobians ) then
           CALL VLIDORT_LPS_master ( do_debug_input, &
           VLIDORT_FixIn,    & ! INPUTS
           VLIDORT_ModIn,    & ! INPUTS (possibly modified)
           VLIDORT_Sup,      & ! INPUTS/OUTPUTS
           VLIDORT_Out,      & ! OUTPUTS
           VLIDORT_LinFixIn, & ! INPUTS
           VLIDORT_LinModIn, & ! INPUTS (possibly modified)
           VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
           VLIDORT_LinOut )    ! OUTPUTS
      else
         CALL VLIDORT_master ( &
           do_debug_input, &
           VLIDORT_FixIn,    & ! INPUTS
           VLIDORT_ModIn,    & ! INPUTS (possibly modified)
           VLIDORT_Sup,      & ! INPUTS/OUTPUTS
           VLIDORT_Out )       ! OUTPUTS
      endif
!  stop 'after first VLIDORT call'

!  debug progress

      !if (mod(w,20).eq.0)write(*,*)'Exact, Progress counter = ', w,VLIDORT_Out%Main%TS_FOURIER_SAVED(1)

!  Exception handling for Input checks. 5/6/15 upgrade, added wavelength number

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_SERIOUS ) then
         DO M = 1, VLIDORT_Out%Status%TS_NCHECKMessages
           Messages(NMessages+2*M-1) = '(VLIDORT_2p7 Message) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_CHECKMessages(M)))
           Messages(NMessages+2*M  ) = '(VLIDORT_2p7_Action ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_ACTIONS(M)))
         ENDDO
         NMessages = NMessages + 2*VLIDORT_Out%Status%TS_NCHECKMessages
         write(C5,'(I5)')w ; Messages(NMessages+1) = '(VLIDORT Input Check, RT Exact # = '//C5//')'
         NMessages = NMessages + 1 ; fail = .true.;  go to 69
      endif

!  Exception handling for Calculation. 5/6/15 upgrade, added wavelength number

      if ( VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Messages(NMessages+1) = '(VLIDORT_2p7 Message) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_MessAGE))
         Messages(NMessages+2) = '(VLIDORT_2p7_Trace ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_TRACE_1))
         Messages(NMessages+3) = '(VLIDORT_2p7_Trace ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_TRACE_2))
         Messages(NMessages+4) = '(VLIDORT_2p7_Trace ) '//Adjustl(TRIM(VLIDORT_Out%Status%TS_TRACE_3))
         NMessages = NMessages + 4
         write(C5,'(I5)')w ; Messages(NMessages+1) = '(VLIDORT Execution, RT Exact # = '//C5//')'
         NMessages = NMessages + 1 ; fail = .true.;  go to 69
      endif

!  Saved results,  Jacobians are normalized, scattering weights not.

      do v = 1, ngeoms
        MSRad(1:ns) = VLIDORT_Out%Main%TS_Stokes (1,v,1:ns,DIR)
        Results_Exact%Stokes_LD_Exact(v,1:ns,w) = MSRad(1:ns)
        if ( do_profile_WFs ) then
          do k = 1, nlayers
            if ( LFvary(k) ) then
              do q = 1, LNvary(k)
                MSJac(1:ns) = VLIDORT_LinOut%Prof%TS_ProfileWF(q,k,1,v,1:ns,DIR)
                Results_LP_Exact%LP_Jacobians_LD_Exact(k,q,v,1:ns,w) = MSJac(1:ns)

!if ( w.eq.1.and.k.eq.1)write(*,*)'VL',k,q,v,Results_LP_Exact%LP_Jacobians_LD_Exact(k,q,v,1,w)

                if (q.eq.1) then
                  fac = - one / Geophys%TotalODs%taug(k,w)
                  Results_LP_Exact%Scatwts_LD_Exact(k,v,w) = MSJac(1) * fac / MsRad(1)
                endif
              enddo
            endif
          enddo
        endif
        if ( DO_SURFACE_WFS ) then
          q = 1 ; MSJac(1:ns) = VLIDORT_LinOut%Surf%TS_Surfacewf(1,1,v,1:ns,DIR)
          Results_LS_Exact%LS_Jacobians_LD_Exact(q,v,1:ns,w) = MSJac(1:ns)
        endif
      enddo

!  Old Code for FO output
!      do v = 1, ngeoms
!        FORad(1:ns)  = VLIDORT_Sup%SS%TS_Stokes_SS(1,v,1:ns,DIR)
!        IF ( DIR.eq.UpIdx ) FORad(1:ns) = FORad(1:ns) + VLIDORT_Sup%SS%TS_STOKES_DB(1,v,1:ns)
!        Results_Exact%Stokes_FO_Exact(v,1:ns,w) = FORAD(1:ns)
!        Results_Exact%Stokes_LD_Exact(v,1:ns,w) = VLIDORT_Out%Main%TS_Stokes(1,v,1:ns,DIR) - Results_Exact%Stokes_FO_Exact(v,1:ns,w)
!        if ( do_profile_WFs ) then
!          do k = 1, nlayers
!            if ( LFvary(k) ) then
!              do q = 1, LNvary(k)
!                FORad(1:ns)  = VLIDORT_LinSup%SS%Prof%TS_ProfileWF_SS(q,k,1,v,1:ns,DIR)
!                IF ( DIR.eq.UpIdx ) FORad(1:ns) = FORad(1:ns) + VLIDORT_LinSup%SS%Prof%TS_ProfileWF_DB(q,k,1,v,1:ns)
!                Results_LP_Exact%LP_Jacobians_FO_Exact(k,q,v,1:ns,w) = FORad(1:ns)
!                if (q.eq.1) then
!                  fac = - one / Geophys%TotalODs%taug(k,w)
!                  Results_LP_Exact%Scatwts_FO_Exact(k,v,w) = FORad(1) * fac / Results_Exact%Stokes_FO_Exact(v,1,w)
!                endif
!              enddo
!            endif
!          enddo
!        endif
!        if ( DO_SURFACE_WFS ) then
!          FORad(1:ns)  = GTzero
!          IF ( DIR.eq.UpIdx ) FORad(1:ns) = VLIDORT_LinSup%SS%Surf%TS_Surfacewf_DB(1,1,v,1:ns)
!          Results_LS_Exact%LS_Jacobians_FO_Exact(q,v,1:ns,w) = FORad(1:ns)
!        endif
!      enddo

!  Timing

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(6) = Exacttimes(6) + e3 - e2  ! ExactRTMTime (VLIDORT)
      endif

!  2e. FO CALL
!  -----------

!  Time

      if (monitor_CPU) call cpu_time(e2)

!  Call

      if ( do_Jacobians ) then
        Call SSV_Integral_ILPS_UP_Optimized &
          ( FO_maxgeometries, FO_maxlayers, FO_maxfine, FO_maxatmoswfs, FO_maxsurfacewfs,             & ! Inputs (dimensioning)
           FO_do_planpar, FO_do_regular_ps, FO_do_enhanced_ps, FO_doNadir, FO_do_sleave,              & ! Inputs (Flags)
           FO_do_sunlight, FO_do_lambertian, FO_do_profilewfs, FO_do_reflecwfs, FO_do_sleavewfs,      & ! Inputs (control, Jacobian )
           FO_Lvaryflags, FO_Lvarynums, FO_n_reflecwfs, FO_n_sleavewfs,                               & ! Inputs (control, Jacobian )
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
           FO_do_sunlight, FO_do_lambertian, FO_nstokes, FO_ngeoms, FO_nlayers, FO_nfinedivs, FO_AcLevel, & ! Inputs (control output)
           FO_reflec, FO_slterm, FO_extinction, FO_deltaus, FO_exactscat_up, FO_flux, FO_fluxvec,         & ! Inputs (Optical)
           FO_Mu0, FO_Mu1, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,                          & ! Inputs (Geometry)
           FO_Raycon, FO_cota, FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,            & ! Inputs (Geometry)
           FO_stokes_up, FO_stokes_db, FO_cumsource_up )                                                    ! Outputs
      endif

!  Saved results. Upwelling only (DIR = 1)

      do v = 1, ngeoms
        FORad(1:ns)= FO_stokes_up(1:ns,v)
        Results_Exact%Stokes_FO_Exact(v,1:ns,w) =  FORad(1:ns)
        IF ( DIR.eq.UpIdx ) Results_Exact%Stokes_FO_Exact(v,1:ns,w) = FORad(1:ns) + FO_stokes_db(1:ns,v)
        if ( do_profile_WFs ) then
          do k = 1, nlayers
            if ( LFvary(k) ) then
              do q = 1, LNvary(k)
                FOJac(1:ns)  = FO_LP_Jacobians_up(1:ns,v,k,q)
                IF ( DIR.eq.UpIdx ) FOJac(1:ns) = FOJac(1:ns) + FO_LP_Jacobians_db(1:ns,v,k,q)
                Results_LP_Exact%LP_Jacobians_FO_Exact(k,q,v,1:ns,w) = FOJac(1:ns)
!if ( w.eq.1.and.k.lt.10)write(*,*)'FO',k,q,v,Results_LP_Exact%LP_Jacobians_FO_Exact(k,q,v,1,w)
                if (q.eq.1) then
                  fac = - one / Geophys%TotalODs%taug(k,w)
                  Results_LP_Exact%Scatwts_FO_Exact(k,v,w) = FOJac(1) * fac / FORad(1)
                endif
              enddo
            endif
          enddo
        endif
        if ( DO_SURFACE_WFS ) then
          q = 1 ; FOJac(1:ns)  = GTzero ; IF ( DIR.eq.UpIdx ) FOJac(1:ns) = FO_LS_Jacobians_db(1:ns,v,q)
          Results_LS_Exact%LS_Jacobians_FO_Exact(q,v,1:ns,w) = FOJac(1:ns)
        endif
      enddo

!  Stop for first-order
!do k = 1, nlayers
!write(*,*)k,FO_LP_Jacobians_up(1,1,k,1),FO_LP_Jacobians_db(1,1,k,1)
!enddo
!stop'FO Jacs'

!  Timing

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(7) = Exacttimes(7) + e3 - e2  ! ExactRTMTime (FO)
      endif

!  2f. 2STREAM call (OPTIONAL)
!  ---------------------------

      if ( DO_OPTIONAL_2STREAM ) then

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

!  Exception handling for Calculation. 5/6/15 upgrade, added wavelength number

        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          Messages(NMessages+1)   = '(2STREAM_2p4 Message) '//Adjustl(TRIM(S2_Message))
          Messages(NMessages+2)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(S2_TRACE_1))
          Messages(NMessages+3)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(S2_TRACE_2))
          NMessages = NMessages + 3
          write(C5,'(I5)')w ; Messages(NMessages+1) = '(Optional 2S Execution, RT Exact # = '//C5//')'
          NMessages = NMessages + 1 ; fail = .true.;  go to 69
        endif

!  Save exact 2stream results. Upwelling only.....
!  Rob Fix 8/3/18. Allowed multiple geometries

        do v = 1, ngeoms
          Int = intensity_2S(v)
          Results_Exact%intensity_2S_Exact(v,w) = Int
          if ( do_profile_WFs ) then
            do k = 1, nlayers
              if ( LFvary(k) ) then
                do q = 1, LNvary(k)
                  Jac = PROFILEWFS_2S(v,k,q) ; Results_LP_Exact%LP_Jacobians_2S_Exact(k,q,v,w) = JAC
                  if (q.eq.1) Results_LP_Exact%ScatWts_2S_Exact(k,v,w) = - Jac / Int / Geophys%TotalODs%taug(k,w)
                enddo
              endif
            enddo
          endif
          if ( do_surface_Wfs ) then
            q = 1 ; Results_LS_Exact%LS_Jacobians_2S_Exact(q,v,w)  = SURFACEWFS_2S(v,q)
          endif
        enddo

!  Timing

        if (monitor_CPU) then
           call cpu_time(e3) ; Exacttimes(8) = Exacttimes(8) + e3 - e2 ! ExactRTMtime (Optional 2S) 
        endif

!  End 2-stream if block

      end if

!  2g. Final Calculations and debug
!  --------------------------------

!  Timing

      if (monitor_CPU) call cpu_time(e2)

!  Final output (Stokes  vector)
      do v = 1, ngeoms
         Results_Exact%Stokes_Exact(v,1:ns,w) = &
                  Results_Exact%Stokes_LD_Exact(v,1:ns,w) + Results_Exact%Stokes_FO_Exact(v,1:ns,w)
     
      !if (w ==138) print * , Results_Exact%Stokes_Exact(v,1,w), Results_Exact%Stokes_LD_Exact(v,1,w), Results_Exact%Stokes_FO_Exact(v,1,w)
!       print * , results_exact%stokes_exact(v, 1:ns, w)     
      enddo

!  develop Fluxes (VLIDORT only)

      if ( Inputs%RTMControl%do_SphericalAlbedo ) then
         do v = 1, ngeoms
            Results_Exact%Fluxes_Exact(v,1:ns,1,w) = VLIDORT_Out%Main%TS_MEAN_STOKES(1,v,1:ns,DIR)  ! Actinic
            Results_Exact%Fluxes_Exact(v,1:ns,2,w) = VLIDORT_Out%Main%TS_FLUX_STOKES(1,v,1:ns,DIR)  ! Flux
         enddo
      endif

!  Develop Polarization: scalar case --> skip DOLP/DOCP calculation (only for NSTOKES = 3 or 4)

      IF ( NSTOKES.gt.1 ) THEN
         do v = 1, ngeoms
            SQ2U2 = SQRT ( Results_Exact%Stokes_Exact(v,2,w) * Results_Exact%Stokes_Exact(v,2,w) &
                         + Results_Exact%Stokes_Exact(v,2,w) * Results_Exact%Stokes_Exact(v,2,w) )
            Results_Exact%DOLP_Exact(V,W) = SQ2U2 / Results_Exact%Stokes_Exact(v,1,w)
            if ( nstokes .eq. 4 ) then
               SV2 = SQRT ( Results_Exact%Stokes_Exact(v,4,w) * Results_Exact%Stokes_Exact(v,4,w) )
               Results_Exact%DOCP_Exact(V,W) = SV2 / Results_Exact%Stokes_Exact(v,1,w)
            endif
         ENDDO
      ENDIF

!  Linearization final output
      do v = 1, ngeoms
        if ( do_profile_WFs ) then
          do k = 1, nlayers
            if ( LFvary(k) ) then
              do q = 1, LNvary(k)

                Results_LP_exact%LP_Jacobians_Exact(k,q,v,1:ns,w) = &
                  Results_LP_exact%LP_Jacobians_LD_exact(k,q,v,1:ns,w) + Results_LP_exact%LP_Jacobians_FO_Exact(k,q,v,1:ns,w)
!if ( k.lt.7)write(443,*)w,Geophys%WavGrids%wav(w), Results_LP_exact%LP_Jacobians_Exact(k,q,v,1:ns,w) 
!!Results_LP_exact%LP_Jacobians_LD_exact(k,q,v,1:ns,w), &
!!                                   Results_LP_exact%LP_Jacobians_FO_Exact(k,q,v,1:ns,w)

              enddo
              fac = - one / Geophys%TotalODs%taug(k,w)
              Results_LP_exact%ScatWts_Exact (k,v,w) = fac * &
                Results_LP_exact%LP_Jacobians_Exact(k,1,v,1,w) / Results_Exact%Stokes_Exact(v,1,w) 

!write(*,*)k,q,Results_LP_exact%ScatWts_Exact(k,w,v),&
!    - ( Results_LP_exact%LP_Jacobians_2S_Exact(k,1,v,w) + Results_LP_exact%LP_Jacobians_FO_Exact(k,1,v,1,w) ) &
! /( Results_Exact%intensity_2S_Exact(v,w) + Results_Exact%Stokes_FO_Exact(v,1,w)  )  / Geophys%TotalODs%taug(k,w)

            endif
          enddo
        endif
        if ( DO_SURFACE_WFS ) THEN
          q = 1 ;  Results_LS_exact%LS_Jacobians_Exact(q,v,1:ns,w) = &
               Results_LS_exact%LS_Jacobians_LD_exact(q,v,1:ns,w) +  Results_LS_exact%LS_Jacobians_FO_Exact(q,v,1:ns,w)
         
        endif
      enddo

!  FD Testing: Debug First Order output. FORT.75/76/77
!
!      if ( w.lt.11 .and.do_debug_output ) then
!         write(*,202)w,intensity_FO_Exact(w,1),FO_intensity_ss(1,1), FO_intensity_db(1,1)
!      endif
!202   format(i4,1p3e20.10)

! Debug output to FORT.8810+v, FORT.9910+v
!   Rob Fix, MultiGoemetry 8/3/18

      if ( do_debug_output ) then
        do v = 1, ngeoms
          if ( .not. do_optional_2stream) Results_Exact%intensity_2S_Exact(v,w) = GTZERO
          if ( nstokes.eq.1 ) then
            write(8810+v,250)w,Geophys%WavGrids%wav(w),&
             Results_Exact%Stokes_LD_Exact(v,1,w),Results_Exact%intensity_2S_Exact(v,w),Results_Exact%Stokes_FO_Exact(v,1,w)
          else
            write(9910+v,250)w,Geophys%WavGrids%wav(w),&
             Results_Exact%Stokes_LD_Exact(v,1:nstokes,w),Results_Exact%intensity_2S_Exact(v,w),&
             Results_Exact%Stokes_FO_Exact(v,1:nstokes,w)
          endif
        enddo
      endif
250   format(i4,1p10e20.10)

!  Timing

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(9) = Exacttimes(9) + e3 - e2  ! ExactRTMtime (Final)
      endif

!  End of main wavelength loop
!if ( w.eq.1)stop'tempo 1 wavelength'

   enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  END PARALLEL REGION @@@@@@@@@@@@@@@@@@@@@@

!  Time spent in equivalent OpenMP parallel region
!    Have to subtract off time spent in 2S+FO setup, FO Exec, and 2S Exec

   if (monitor_CPU) then
      call cpu_time(ser_e2) ; SerExactTime = ser_e2 - ser_e1 - sum(ExactTimes(4:5))- sum(ExactTimes(7:8))
      ExactTimes(10) = SerExactTime
   endif

!  Final CPU call

   if (monitor_CPU) then
      call cpu_time(e3) ; ExactRuntime = ExactRuntime + e3-e1
   endif

!  Output Variables
   IF (allocated (RO%stokes) ) THEN 
      deallocate (RO%stokes, RO%PJAC, RO%SJAC)
   ENDIF

   allocate (RO%STOKES(4, nwav))
   allocate (RO%PJAC (nlayers,n_totalatmos_wfs, 4, nwav))
   allocate (RO%SJAC (n_surface_wfs,  4, nwav))
   RO%Stokes(1:4, 1:nwav)   = Results_exact%Stokes_Exact(1, :, 1:nwav)
   RO%PJAC(1:nlayers, 1:n_totalatmos_wfs, 1:4, 1:nwav) = & 
           Results_LP_exact%LP_Jacobians_Exact(1:nlayers, 1:n_totalatmos_wfs, 1, :, 1:nwav)
   RO%SJAC(1:n_surface_wfs,1:4, 1:nwav)  = Results_LS_exact%LS_Jacobians_Exact(1:n_surface_wfs, 1, :, 1:nwav)
!  Normal Finish
   RETURN

!  Error Finish

69 continue

   RETURN
END SUBROUTINE VLIDORT_Exact_master

!  End module

end module VLIDORT_exact_module

