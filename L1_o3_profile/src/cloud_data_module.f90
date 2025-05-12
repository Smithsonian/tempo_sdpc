MODULE cloud_data_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, mreflcld, max_fit_pts


  IMPLICIT NONE

  INTEGER, PARAMETER    :: maxalb = 4, maxwfc = 4, maxvars = 10, maxgases = 3
  INTEGER, PARAMETER    :: maxcldvar = maxalb + maxwfc + maxvars + maxgases
  INTEGER, PARAMETER    :: maxit = 10,  maxnfr=11001
  REAL(KIND=8), PARAMETER                :: ATM1 = 1013.25d0
  INTEGER                                :: ringOrder, shiOrder
  INTEGER                                :: the_solwin
  INTEGER                                :: nalb, nwfc
  INTEGER                                :: n_fitvar_cld, nfitwav_cld
  REAL(KIND=dp)                          :: totozone, the_grref, the_clref, tozap

  REAL(KIND=dp), DIMENSION (max_fit_pts) :: currwvls, currweights, currspec
  REAL(KIND=dp), DIMENSION(mreflcld)     :: fitwvls
  REAL(KIND=dp), DIMENSION(mreflcld)     :: fitspec_cld, fitres_cld, weights, simrad_cld

  REAL(KIND=dp), DIMENSION(maxcldvar)    :: fitvar_cld, fitvar_cld_saved, fitvar_cld_aperror, &
       fitvar_cld_init, fitvar_cld_apriori, fitvar_cld_std, fitvar_cld_nstd, fitvar_cld_dfs

  REAL(KIND=dp), DIMENSION(maxcldvar)    :: alldfs
  REAL(KIND=dp), DIMENSION(6)            :: tozrrs
  REAL(KIND=dp), DIMENSION(maxcldvar, maxcldvar)   :: covar, ncovar
  INTEGER         :: ctp_idx, ps_idx, cfrac_idx,  alb_idx, &
       poly_fidx, rc_idx, toz_idx, shi_idx!, radd_idx

  ! Overall control flags
  LOGICAL         ::  ret_cfrac, ret_ctp, ret_ps, ret_alb, ret_ringc, &
       update_toz, ret_shi
  LOGICAL         :: use_o2o2, use_rrs, wrt_retcld, use_retctp, use_retalb, &
       wrt_rescld, debug, do_subwin
  REAL(KIND=dp)   :: cld_swav, cld_ewav, toz_swav, toz_ewav, cfrac_thread, the_cfrac1

  ! Fitting control flags
  LOGICAL         :: use_tomsv8_clima
  LOGICAL         :: has_ps, has_cfrac, has_alb, has_ctp, has_toz, two_steps
  LOGICAL         :: has_ringc, ret_toz, has_shi , has_ringc_toz, has_shi_toz!, has_radd
  LOGICAL         :: do_cloud_only, do_kai_method, softcali_forcld
  REAL(KIND=dp)   :: o2o2swav, o2o2ewav,  rrsswav, rrsewav
  CHARACTER(LEN=6), DIMENSION(maxcldvar)           :: cldvarname

  ! Trace gases
  INTEGER                                       :: ngases, gas_idx    
  CHARACTER(LEN=4), DIMENSION(maxgases)         :: which_gases        ! OCLO, BRO   
  LOGICAL                                       :: ret_gas, has_gas
  !CHARACTER(LEN=maxchlen),dimension(maxgases)   :: gastbl_fname
  REAL(KIND=dp), DIMENSION (maxgases)              :: gas_init_col
  REAL(KIND=dp), DIMENSION (max_fit_pts, maxgases) :: tracegas_xsecs

  ! Ocean Raman
  LOGICAL                   :: ret_vraman, ret_chl, has_ocean, inc_chl_only
  CHARACTER(LEN=maxchlen)   :: vraman_fname
  INTEGER                   :: vraman_idx
  REAL(KIND=dp)             :: chlorophyll

  ! File unit
  INTEGER, PARAMETER    :: oc_raman_unit = 40,  cldfit_ctrl_unit = 41

  ! CLDO2 and CLDRRS data
  LOGICAL          :: get_opt_cld
  REAL(KIND=dp)    :: o2o2_ctp, o2o2_ps0, o2o2_cfrac, rrs_ctp, rrs_ps0, rrs_cfrac

  ! Final results
  INTEGER          :: niter_cld, exval_cld, n_newalb, n_newwfc
  REAL(KIND=dp)    :: avgwav_cld, rms_cld, avgres_cld
  REAL(KIND=dp)    :: new_ps0, new_ctp, new_toz, ps0_kai, ctp_kai, cfrac_kai, AIb_kai
  !character(LEN=4), dimmension(10), parameter   :: all_cldvarname
  REAL(KIND=dp), DIMENSION(maxalb) :: new_alb
  REAL(KIND=dp), DIMENSION(maxwfc) :: new_cfrac
  ! Saving 
  INTEGER  :: numwin_saved, nradwvl_saved, n_refwvl_saved, n_refwvl_sav_saved

  ! Final results, include all possible parameter
  ! 1-6: a-priori, a-priori error, ret value, ret. error, Rnd. error, dfs
  !INTEGER                                 :: num_retcld
  !INTEGER      , DIMENSION(maxcldvar)     :: retcld_idx
  !REAL(KIND=dp), DIMENSION(maxcldvar, 6)  :: retcld_results
  !CHARACTER(LEN=6), DIMENSION(maxcldvar)  :: retcld_varname

  ! All possible variables, note that rnc and shi are fixed
  ! ectp
  ! sprs
  ! alb0
  ! alb1
  ! alb2
  ! cfr0
  ! cfr1
  ! cfr2
  ! rnc0
  ! rnc1
  ! tocl
  ! shi0
  ! vram



END MODULE cloud_data_module
