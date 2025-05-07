!>
MODULE OMSAO_omidata_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, n_rad_winwav, maxwin, &
       max_fit_pts, max_ring_pts, mrefl
  USE OMSAO_indices_module,    ONLY: n_max_fitpars, max_rs_idx, &
       max_calfit_idx, o3_t1_idx, o3_t3_idx, spc_idx
 ! USE OMSAO_variables_module, ONLY: cali_group, irrad_group, rad_group, refl_group, ring_group, geo_group
  IMPLICIT NONE

  ! max number of channels (UV-1 and UV-2) reqd for OMI O3 profile retrievals
  INTEGER, PARAMETER :: mswath  = 2
  INTEGER, PARAMETER :: zoom_p1 = 16, zoom_p2 = 45
  LOGICAL            :: zoom_mode
  ! ---------------------------------
  ! Maximum OMI data/swath dimensions
  ! ---------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       ntimes_max     = 1644, nxtrack_max  =  60, nwavel_max  = 750, &
       nwavelcoef_max =    5, nlines_max = 100

  ! ---------------------------------
  ! Dimensions for measurement swaths
  ! ---------------------------------
  ! nfxtrack: is the number of pixels after coadding UV-2
  INTEGER  :: nfxtrack

  ! ---------------------------------
  ! OMI swath names (UV-1 and UV-2)
  ! --------------------------------
  CHARACTER (LEN=maxchlen), DIMENSION(mswath), PARAMETER :: & 
     omi_irradiance_swathname=(/'Sun Volumne UV-1 Swath','Sun Volumne UV-2 Swath'/), &
     omi_radiance_swathname=(/'Earth UV-1 Swath','Earth UV-2 Swath'/)

  ! ------------------------------------------------------------
  ! Boundary wavelengths (approximate) for UV-2 and VIS channels
  ! ------------------------------------------------------------
  REAL (KIND=r4), DIMENSION(mswath),  PARAMETER :: &
       upper_wvls = (/310.0, 387.0/), lower_wvls = (/260.0, 310.0/)
  REAL (KIND=dp), PARAMETER :: lower_spec = 0.0, upper_spec = 4.0E14

  ! ---------------------------------------
  ! Minimum OMI spectral resolution (in nm)
  ! ---------------------------------------
  REAL (KIND=r8), PARAMETER :: omi_min_specres = 0.5_r8
  INTEGER (KIND=i4)         :: omisol_version
  ! --------------------------------------------------
  ! Parameters defined by the NISE snow cover approach
  ! --------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: NISE_snowfree =   0, NISE_allsnow = 100, &
                                  NISE_permice = 101, NISE_drysnow = 103, &
                                   NISE_ocean    = 104, NISE_suspect = 125, NISE_error   = 127

  ! --------------------------------------------------------------------
  ! A diagnostic array that shows how the AMF was computed, with values
  ! that indicate missing cloud products, glint, and geometric or no AMF
  ! --------------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       omi_cfr_addmiss = 1000, omi_ctp_addmiss = 2000, omi_glint_add = 10000, &
       omi_geo_amf = -1, omi_oobview_amf = -2
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1) :: amf_diagnostic

  ! ---------------------------------------
  ! Swath attributes for measurement swaths
  ! ---------------------------------------
  INTEGER (KIND=i4), DIMENSION (nxtrack_max):: n_omi_database_wvl
  INTEGER (KIND=i2), DIMENSION (nxtrack_max)                         :: &
       omi_solcal_itnum, omi_radcal_itnum, omi_radref_itnum,            &
       omi_solcal_xflag, omi_radcal_xflag, omi_radref_xflag
  REAL (KIND=r8), DIMENSION (max_calfit_idx, nxtrack_max)         :: &
       omi_solcal_pars,  omi_radcal_pars,  omi_radref_pars
  REAL (KIND=r8), DIMENSION ( nwavel_max, nxtrack_max) :: omi_database_wvl
  REAL (KIND=r8), DIMENSION (nxtrack_max) :: omi_sol_wav_avg
  REAL (KIND=r8), DIMENSION (nxtrack_max) :: &
       omi_solcal_chisq, omi_radcal_chisq, omi_radref_chisq, &
       omi_radref_col,   omi_radref_dcol,  omi_radref_rms
  REAL (KIND=r8), DIMENSION (2,nxtrack_max,0:nlines_max-1) :: omi_wavwin_rad, &
       omi_fitwin_rad
  REAL (KIND=r8), DIMENSION (2,nxtrack_max) :: omi_wavwin_sol, omi_fitwin_sol

  INTEGER (KIND=i4), DIMENSION (nxtrack_max) :: omi_solfit_xflag
  REAL    (KIND=dp), DIMENSION (max_calfit_idx, nxtrack_max) :: omi_solfit_pars
  INTEGER (KIND=i4), DIMENSION (nxtrack_max, n_rad_winwav) :: &
       omi_rad_winwav_idx


  INTEGER (KIND=i4), PARAMETER :: qa_percent_passed = 80, &
       qa_percent_suspect = 50
 
  !TYPE (cali_group) :: omi_cali
  !TYPE (irrad_group):: omi_irrad
  !TYPE (rad_group)  :: omi_rad
  !TYPE (ring_group) :: omi_ring
  !TYPE (refl_group) :: omi_refl
  !TYPE (geo_group) :: omi_geo
  ! OUTPUT variables
  CHARACTER(len=*), PARAMETER :: l2_swathname = 'OMI Vertical Ozone Profile'

END MODULE OMSAO_omidata_module
