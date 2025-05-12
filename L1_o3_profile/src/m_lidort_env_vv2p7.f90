MODULE m_LIDORT_ENV_vv2p7

  USE VLIDORT_PARS
  USE VLIDORT_IO_DEFS
  USE VLIDORT_LIN_IO_DEFS

  USE VLIDORT_AUX
  USE VLIDORT_L_INPUTS
  USE VLIDORT_LPS_MASTERS
 

  TYPE(VLIDORT_Fixed_Inputs), SAVE       :: VLIDORT_FixIn
  TYPE(VLIDORT_Modified_Inputs), SAVE    :: VLIDORT_ModIn
  TYPE(VLIDORT_Fixed_LinInputs), SAVE    :: VLIDORT_LinFixIn
  TYPE(VLIDORT_Modified_LinInputs), SAVE :: VLIDORT_LinModIn

  INTEGER, SAVE :: nstreams, nstokes, nlayers, n_atmos_wfs, n_sfc_wfs
 ! VLIDORT INPUT
  INTEGER, DIMENSION (maxlayers), SAVE :: layer_vary_number
  LOGICAL, DIMENSION (maxlayers), SAVE :: layer_vary_flag
  CHARACTER (LEN=31), DIMENSION (max_atmoswfs), SAVE :: profilewf_names

  ! dummy for splat
  REAL(KIND=8), PARAMETER    :: I0_conversionFactor = 1.0

  PUBLIC lidort_prof_env 
  PRIVATE
! ******************************************************************************************
! When use_effcrs is set to TRUE, outputs (radiance and weighting function) correspond to 
! input wavelengths. 
! When use_effcrs is set to FALSE, i.e., use high resolution cross 
! section in the radiative transfer calculation and convolve with slot functions after
! radiative transfer calculation, waves do not correspond to the output rad and weighting 
! function, waves(1:nw) store the wavelengths (waves(1:ncalcp) where radiance calculation 
! are done (waves(ncalp+1:nw) = 0.0. And outputs corresponds to measurents (1:ns).
! Here nw is max(ncalcp, ns).
! ******************************************************************************************
 CONTAINS
SUBROUTINE LIDORT_PROF_ENV(do_ozwf, do_albwf, do_tmpwf, do_o3shi, &
     do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
     do_sprswf, do_so2zwf, do_pslwf, ozvary,&
     nw, waves, nalbwf, albs, wfcs,nos, o3shi, nl, ozprof, tprof, nostk, rto, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY  : du2mol, rearth
  USE OMSAO_indices_module,   ONLY  : so2_idx, o2o2_idx, &
                                      o2_idx, h2o_idx !so2v_idx
  USE OMSAO_variables_module, ONLY  : scnwrt, ctrdbdir, &
       sza=>the_sza_atm, vza=>the_vza_atm, aza=>the_aza_atm, lat=>the_lat, currloop, &
       numwin, winlim, band_selectors, wcenter_uvvis, &
       nradpix, n_rad_wvl, fitwavs, radwvl_sav, n_radwvl_sav, & 
       mask_fitvar_rad, &
       fitvar_rad, fitvar_rad_apriori, &
       refidx, refspec_norm, the_surfalt, database_pslwf, npsl, tabdir
  USE ozprof_data_module,     ONLY : num_iter, do_debug_o3p, & 
       do_tracewf, use_effcrs, &
       ncalcp, radcwav, &
       mpolcorr, npolcorr,polcorr, polcorr_idxs, the_str, do_radinter, & 
       mflay, nfsfc, nsfc, ntp, nt_fit, nflay,nup2p, &
       atmosprof, nup2p, fts, fps, fzs, fozs, frhos,&
       nallgas, ngas, mgasprof, tracegas, fgasidxs, &
       do_lambcld, lambcld_refl, has_clouds, aerosol, strat_aerosol,&
       nmom, aerwavs,actawin,useasy, maxgksec, ngksec, &
       the_cod, the_cfrac, the_cbeta, ncbp, nctp, nwfc, &
       gaext, gasca, gaasy, gamoms, gcq, gcasy, gcmoms, &
       tropsca, tropaod, tropwaer, strataod, stratsca, &
       taodind, saodind, twaeind, sprsind, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       so2idx, o4idx, o2idx, h2oidx, broidx, hchoidx, no2idx, &
       ccrs, hcrs, dadp, dadt, dads,  & 
       rtm_outputs, set_rtmvar, do_brdf
  USE OMSAO_errstat_module
  USE m_get_xcrs, ONLY: get_alb_ozcrs_ray , &
                        get_hres_gascrs_ray, get_effres_gascrs_ray
  USE m_lidort_util, ONLY: set_polcorr, polcorr_online, polcorr_online_with_lut, &
                           radwf_interpol, hres_radwf_inter_convol,&
                           get_tracegas_wf, debug_rtm, debug_taug, set_polcorr_new
  USE m_ezspline_interpolation, ONLY: bspline
  USE m_set_brdf, ONLY: error, surface, SurfProf, window, Geolocation
  USE Surface_module, ONLY: SetSurfaceLinearization,SetSurfaceOpticalProperties

  implicit none
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)  :: nw, nl, nos, nostk, nalbwf
  LOGICAL, INTENT(IN)  :: do_ozwf,   do_albwf,  do_tmpwf, do_o3shi,&
                          do_taodwf, do_twaewf, do_saodwf, &
                          do_cfracwf,do_codwf,  do_ctpwf, &
                          do_sprswf, do_so2zwf, do_pslwf
  LOGICAL, DIMENSION(nl), INTENT(IN)           :: ozvary
  REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)   :: waves, albs, wfcs
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)   :: o3shi
  REAL (KIND=dp), DIMENSION(nl),  INTENT(IN)           :: ozprof, tprof
  TYPE (rtm_outputs), INTENT(OUT) :: rto
  INTEGER, INTENT(OUT) :: errstat
  ! =======================
  ! Local variables
  ! =======================
  INTEGER, PARAMETER :: nlowstr = 2
  LOGICAL :: problems, do_clouds, do_fozwf, do_faerwf, do_fraywf, do_abs
  INTEGER :: ic, iw, ialbwf, i, j, k, low, hgh, fidx, lidx, iwf, idx, &
             nstep, istk, npolmod, ipol, nsprs, nz1, nfgas, nw0
  INTEGER :: ozwfidx, aodwfidx, twaewfidx, codwfidx, sprswfidx, raywfidx
  INTEGER, DIMENSION (nallgas)                  :: gasin   
  INTEGER, DIMENSION (5)                        :: tmp_gasidxs
  LOGICAL, DIMENSION(nflay):: cldmsk, varyprof!, aermsk
  REAL (KIND=dp)  :: lamda, xg, frac, toz, temp, aodscl, waerscl,lambertian_albedo
  REAL (KIND=dp), DIMENSION(0:nflay)           :: ozs, delps
  REAL (KIND=dp), DIMENSION(nflay)             :: cldsca, cldext, cldasy, cldext0, aersca, aerext, aerasy
  REAL (KIND=dp), DIMENSION(0:nmom,maxgksec,nflay) :: cldmoms, aermoms
  REAL (KIND=dp)                               :: tmpalb
  REAL (KIND=dp), DIMENSION(nw)                :: tmpalbs
  ! rad/wf variables
  REAL (KIND=dp), DIMENSION(2, nostk)          :: radclrcld
  REAL (KIND=dp), DIMENSION(nw, nostk)         :: rad, &
       cfracwf, o3shiwf, codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(nw, nalbwf, nostk) :: albwf
  REAL (KIND=dp), DIMENSION(nw, nl, nostk)     :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(nw, nflay       )  :: tauwf
  REAL (KIND=dp), DIMENSION(nw, nflay, nostk)  :: fozwf, &
       faerwf, faerswf, fcodwf, fsprswf, fraywf

  ! polarization correction variables
  ! @ lut-correction
  LOGICAL :: do_plutcorr_after
  REAL (KIND=dp) :: albclrcld(nw, 2)
  REAL (KIND=dp), DIMENSION(nw, 2, nostk)      :: polerr
   CHARACTER(LEN=100)                          :: VLDLUTdir
  ! @ on-correction
  REAL (KIND=dp), DIMENSION(mpolcorr, nflay)   :: ptauwf,pfozwf, &
        pfaerwf, pfaerswf, pfcodwf, pfsprswf, pfraywf
  REAL (KIND=dp), DIMENSION(mpolcorr)          :: prad, pctpwf, pcfracwf
  REAL (KIND=dp), DIMENSION(mpolcorr, nalbwf)  :: palbwf

  ! cross section & optical depth variables
  REAL (KIND=dp), DIMENSION(nw)                 :: raycof, depol
  REAL (KIND=dp), DIMENSION(nw, nallgas, nflay) :: allcrs
  REAL (KIND=dp), DIMENSION(nw, nflay)          :: abscrs
  REAL (KIND=dp), DIMENSION(nallgas, nflay)     :: alleta, allcol  
  REAL (KIND=dp), DIMENSION(nw, nflay) :: deltau, delsca, delray, delabs
  ! Others
  REAL (KIND=dp) :: e_loop1, e_loop2, e_s, e_n ,e1, e2, e_vlidort, e_pol, e_inter
  ! Save variables
  LOGICAL, SAVE                                 :: first = .TRUE.
  INTEGER, SAVE                                 :: nz, faer_lvl, nradcal
  REAL (KIND=dp), DIMENSION(0:mflay),      SAVE :: ts, ps 
  LOGICAL, DIMENSION(mflay),               SAVE :: aermsk
  LOGICAL, DIMENSION(:),ALLOCATABLE,       SAVE :: do_radcals, do_polcorrs
  INTEGER, DIMENSION(:),ALLOCATABLE,       SAVE :: radcal_idxs, polidx

  ! Will become an input parameter later in ozprof.inp
  LOGICAL :: do_ssfullb295 = .FALSE.
  LOGICAL, PARAMETER :: do_debug_input = .FALSE., do_debug_rtm = .FALSE.

  ! ========================================================================
  ! VLIDORT variables
  ! ========================================================================

!  VLIDORT file inputs status structure
   TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus
!  VLIDORT output structure
   TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out
   TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut
!  VLIDORT supplements i/o structure
   TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup
   TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup
   CHARACTER (LEN=15) :: vldlut = 'vldlut'
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=15), PARAMETER :: modulename = 'LIDORT_PROF_ENV'
  
  IF (do_debug_o3p) WRITE(www_lun, '(i3,A, i4)') num_iter, ': lidort_prof_env (start)', nw
  CALL cpu_time(e_s)
  do_fozwf = .FALSE.
  do_faerwf = .FALSE.
  do_fraywf = .FALSE.
  errstat = pge_errstat_ok

  IF (first) THEN
    ! ======================= Read LIDORT Control Input ==========================
    !xliu: 03/07/2011, switch VLIDORT from vv2p4 to vv2p4RTC
    !CALL VLIDORT_L_INPUT_MASTER ('INP/vlidort_control.inp', &
    !     'o3prof_lidort_error', status_inputread)
    WRITE(*,*) 'VLIDORT VV2P7, polcorr=', polcorr
    CALL VLIDORT_L_INPUT_MASTER (& 
       TRIM(ADJUSTL(ctrdbdir))//'vlidort_control_vv2p7.inp', & !Input
       VLIDORT_FixIn,      & ! Outputs
       VLIDORT_ModIn,      & ! Outputs
       VLIDORT_LinFixIn,   & ! Outputs
       VLIDORT_LinModIn,   & ! Outputs
       VLIDORT_InputStatus ) ! Outputs
   
    IF (VLIDORT_InputStatus%TS_STATUS_INPUTREAD .NE. 0 ) THEN 
       WRITE(*,*) 'Errors in VLIDORT_L_INPUT_MASTER' 
       open(1,file = 'V2p7_Profiling_ReadInput.log',status = 'unknown')
       WRITE(1,*)' FATAL:   Wrong input from VLIDORT input file-read'
       WRITE(1,*)'  ------ Here are the messages and actions '
       write(1,'(A,I3)')'    ** Number of messages = ',VLIDORT_InputStatus%TS_NINPUTMESSAGES
       DO i = 1, VLIDORT_InputStatus%TS_NINPUTMESSAGES
        write(1,'(A,I3,A,A)')'Message # ',i,' : ',Trim(VLIDORT_InputStatus%TS_INPUTMESSAGES(i))
        write(1,'(A,I3,A,A)')'Action  # ',i,' : ',Trim(VLIDORT_InputStatus%TS_INPUTACTIONS(i))
       ENDDO
      close(1)
      STOP 1
    ENDIF

    VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS = 0 
    VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF = 0 
    VLIDORT_ModIn%MCont%TS_ngreek_moments_input  = nmom  
    VLIDORT_ModIn%MChapman%TS_earth_radius = rearth

    IF (.NOT. aerosol .AND. (.NOT. has_clouds .OR. do_lambcld)) THEN
      VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .FALSE. ! DO_SSCORR_TRUNCATION = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.     ! DO_DELTAM_SCALING = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY = .TRUE.         ! DO_RAYLEIGH_ONLY = .TRUE.
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .FALSE.      ! DO_SOLUTION_SAVING = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .FALSE.   
    ENDIF

    first = .FALSE.
  ENDIF
  ! ============= Overridden some control and atmospheric variables ============== 
  IF (num_iter == 0 ) THEN

     !**geometries related variables
     VLIDORT_ModIn%MSunRays%TS_N_SZANGLES = 1
     VLIDORT_ModIn%MSunRays%TS_SZANGLES = sza
     IF (sza >= 90.0 .OR. sza < 0) THEN
        WRITE(*, *) modulename, ' : SZA is >= 90 or < 0 !!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 1
     VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1)= vza
     VLIDORT_ModIn%MUserVal%TS_n_user_relazms= 1
     VLIDORT_ModIn%MUserVal%TS_user_relazms(1) = aza
     VLIDORT_ModIn%MUserVal%TS_n_user_obsgeoms = 1
     VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1, 1:3)= (/sza,vza,aza/)

     ! ** geophysical related variables

     nz = nflay;   nlayers = nz

     IF (nz > maxlayers) THEN
        WRITE(*, *) modulename, ' : # of layers exceeded allowed !!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
     IF (nl > nlayers) THEN
        WRITE(*, *) modulename, ' : Coarse grids cannot be finer than fine grids!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     ! set the first aerosol layer (from TOA down) 1:faer_lvl-1 (without aerosols)
     faer_lvl = 1;   IF (.NOT. strat_aerosol) faer_lvl = nup2p(ntp) + 1

     IF (aerosol) THEN 
        aermsk = .TRUE.
        IF (.NOT. strat_aerosol) aermsk(1:faer_lvl-1) = .FALSE.
     ELSE
        aermsk = .FALSE.
     ENDIF

     ts(1:nz) = (fts(1:nz) + fts(0:nz-1)) / 2.0
     ps(1:nz) = exp((log(fps(1:nz)) + log(fps(0:nz-1))) / 2.0)
     ! Set height grid for doing Chapman Function Calculation 
     VLIDORT_FixIn%Chapman%TS_FINEGRID         = ZERO
     VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:nz)    = fps(0:nz)
     VLIDORT_FixIn%Chapman%TS_height_grid(0:nz)      = fzs(0:nz)
     VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(1:nz) = ts(1:nz)
     VLIDORT_ModIn%MUserVal%TS_geometry_specheight  = the_surfalt
      
     VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC  = ZERO
     VLIDORT_Sup%BRDF%TS_BRDF_F_0         = ZERO
     VLIDORT_Sup%BRDF%TS_BRDF_F           = ZERO
     VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0    = ZERO
     VLIDORT_Sup%BRDF%TS_USER_BRDF_F      = ZERO
     VLIDORT_Sup%BRDF%TS_EMISSIVITY       = ZERO
     VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY  = ZERO

     VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
     VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES = ZERO
     VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0        = ZERO
     VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

     VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
     VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
     VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
     VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
     VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO
     VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
     VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

     VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS = .FALSE.
     VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS   = 0

     VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
     VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
     VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
     VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

 ENDIF

  ! =================== Interpolate Ozone, T to fine grids =======================
  IF (nz /= nl) THEN
     DO i = 1, nl
        varyprof(nup2p(i-1)+1:nup2p(i)) = ozvary(i)
        ozs(nup2p(i-1)+1:nup2p(i)) = fozs(nup2p(i-1)+1:nup2p(i)) * ozprof(i) / &
             SUM(fozs(nup2p(i-1)+1:nup2p(i)))
     ENDDO
     IF (nt_fit > 0) THEN

        DO i = 1, nl
           atmosprof(3, i) = tprof(i) * 2.0 - atmosprof(3,  i-1)
        ENDDO

        CALL BSPLINE(atmosprof(2, 0:nl), atmosprof(3, 0:nl), &
             nl+1, fzs(0:nz), ts(0:nz), nz+1, errstat)
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ' : BSPLINE error, errstat = ', errstat; RETURN
        ENDIF
        ts(1:nz) = (ts(1:nz) + ts(0:nz-1)) / 2.0
        ps(1:nz) = exp((log(ps(1:nz)) + log(ps(0:nz-1))) / 2.0)
     ENDIF
  ELSE
     ozs(1:nz) = ozprof(1:nl);    varyprof(1:nz) = ozvary(1:nl)
  ENDIF

  nz1 = nfsfc - 1 
  ! Update aerosol fields: first AOD
  IF (do_taodwf .AND. num_iter > 0) THEN
     aodscl = fitvar_rad(taodind) / tropaod(actawin)
     tropaod(1:actawin) = tropaod(1:actawin) * aodscl
     gaext(1:actawin, nup2p(ntp)+1:nz1) = gaext(1:actawin, nup2p(ntp)+1:nz1) * aodscl

     IF (.NOT. do_twaewf) THEN ! Single scattering albedo does not change
        tropsca(1:actawin) = tropsca(1:actawin) * aodscl
        gasca(1:actawin, nup2p(ntp)+1:nz1) = gasca(1:actawin, nup2p(ntp)+1:nz1) * aodscl
     ENDIF
  ELSE
     aodscl = 1.0
  ENDIF

  IF (do_twaewf .AND. num_iter > 0) THEN
     waerscl = fitvar_rad(twaeind) / tropwaer(actawin)   ! Scale single scattering albedo
     tropwaer(1:actawin) = tropwaer(1:actawin) * waerscl
     tropsca(1:actawin)  = tropsca(1:actawin) * waerscl * aodscl
     gasca(1:actawin, nup2p(ntp)+1:nz1) = gasca(1:actawin, nup2p(ntp)+1:nz1) * waerscl * aodscl
  ENDIF

  IF (do_saodwf .AND. num_iter > 0) THEN
     aodscl = fitvar_rad(saodind) / strataod(actawin)
     strataod(1:actawin) = strataod(1:actawin) * aodscl
     stratsca(1:actawin) = stratsca(1:actawin) * aodscl
     gaext(1:actawin, 1:nup2p(ntp)) = gaext(1:actawin, 1:nup2p(ntp)) * aodscl
     gasca(1:actawin, 1:nup2p(ntp)) = gasca(1:actawin, 1:nup2p(ntp)) * aodscl
  ENDIF

  IF (do_sprswf .AND. num_iter > 0) THEN
     temp = (fitvar_rad(sprsind) - fps(nup2p(nsfc-1)))/ (fps(nz1) - fps(nup2p(nsfc-1)))
     frhos(nup2p(nsfc-1)+1:nz1) = frhos(nup2p(nsfc-1)+1:nz1) * temp
     delps(nup2p(nsfc-1)+1:nz1) = (fps(nup2p(nsfc-1)+1:nz1)-fps(nup2p(nsfc-1):nz1-1))  * temp
     DO i = nup2p(nsfc-1)+1, nz1
        fps(i) = fps(i-1) + delps(i)
     ENDDO
  ENDIF

  ! Determinine atmospheric weighting functions to be calculated
  IF (do_ozwf .OR.  do_tmpwf .OR. do_o3shi .OR. do_taodwf .OR. &
      do_twaewf .OR. do_saodwf .OR. do_codwf .OR. do_sprswf )  & 
      VLIDORT_LinModIn%MCont%TS_do_atmos_linearization = .TRUE.

  VLIDORT_LinModIn%MCont%TS_do_surface_linearization = do_albwf    
  IF (VLIDORT_LinModIn%MCont%TS_do_atmos_linearization .OR. do_albwf) THEN
     VLIDORT_LinModIn%MCont%TS_do_simulation_only = .FALSE.
     VLIDORT_LinModIn%MCont%TS_do_linearization   = .TRUE.
  ELSE 
     VLIDORT_LinModIn%MCont%TS_do_simulation_only = .TRUE.
     VLIDORT_LinModIn%MCont%TS_do_linearization = .FALSE.
  ENDIF

  i = 0;  ozwfidx = 0; aodwfidx = 0; twaewfidx = 0; codwfidx = 0; sprswfidx = 0
  layer_vary_flag(1:nz) = .FALSE.
  layer_vary_number(1:nz) = 0
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi) THEN
     i = i + 1;  ozwfidx = i
     profilewf_names(i) = 'ozone volume mixing ratio------'

     layer_vary_flag(1:nz) = varyprof(1:nz)
     WHERE (varyprof(1:nz) )
        layer_vary_number(1:nz) = & 
        layer_vary_number(1:nz) + 1
     ENDWHERE
  ENDIF

  IF ( do_taodwf .OR. do_saodwf) THEN
     i = i + 1; aodwfidx = i
     profilewf_names(i) = 'aerosol extinction coefficient-'

     IF (do_taodwf) THEN
        layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(ntp)+1:nz1) = & 
        layer_vary_number(nup2p(ntp)+1:nz1) + 1
     ENDIF

     IF (do_saodwf) THEN
        layer_vary_flag(1:nup2p(ntp)) = .TRUE.
        layer_vary_number(1:nup2p(ntp)) = & 
        layer_vary_number(1:nup2p(ntp)) + 1
     ENDIF
  ENDIF
  IF ( do_twaewf) THEN
     i = i + 1; twaewfidx = i
     profilewf_names(i) = 'aerosol scattering coefficient-'
     layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
     layer_vary_number(nup2p(ntp)+1:nz1) = & 
     layer_vary_number(nup2p(ntp)+1:nz1) + 1
  ENDIF
  IF ( do_codwf) THEN
     i = i + 1; codwfidx = i
     profilewf_names(i) = 'cloud extinction coefficient---'
     layer_vary_flag(nctp:ncbp) = .TRUE.
     layer_vary_number(nctp:ncbp) = & 
     layer_vary_number(nctp:ncbp) + 1
  ENDIF

  do_fraywf = .FALSE.
  IF ( do_sprswf  .OR. (.NOT. use_effcrs .AND. nw > 1) ) THEN
     ! Need to use jacobians wrt rayleigh OD to perform interpolation

     i = i + 1; sprswfidx = i; raywfidx = i
     profilewf_names(i) = 'rayleigh optical thickness-----'
     IF (.NOT. use_effcrs) THEN
        do_fraywf = .TRUE.
        layer_vary_flag(1:nz1) = .TRUE.
        layer_vary_number(1:nz1) = & 
        layer_vary_number(1:nz1) + 1
     ELSE
        layer_vary_flag(nup2p(nsfc-1)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(nsfc-1)+1:nz1) = & 
        layer_vary_number(nup2p(nsfc-1)+1:nz1) + 1
     ENDIF
  ENDIF

  n_atmos_wfs = i ;  n_sfc_wfs =1
 
  IF ( n_atmos_wfs   > 1) THEN
     WHERE(layer_vary_number(1:nz) > 0) 
        layer_vary_number(1:nz) = & 
        n_atmos_wfs 
     ENDWHERE
  ENDIF

  do_fozwf = .FALSE.
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi .OR. (.NOT. use_effcrs .AND. nw > 1)) do_fozwf = .TRUE.
  do_faerwf = .FALSE.
  IF (do_taodwf .OR. do_saodwf) do_faerwf = .TRUE.

  ! copy to vlidort variables
  
  VLIDORT_LinFixIn%Cont%TS_n_totalprofile_wfs = n_atmos_wfs 
  VLIDORT_LinFixIn%Cont%TS_n_surface_wfs = n_sfc_wfs
  VLIDORT_LinFixIn%Cont%TS_layer_vary_flag(1:nz) = layer_vary_flag(1:nz)
  VLIDORT_LinFixIn%Cont%TS_layer_vary_number(1:nz) = layer_vary_number(1:nz)
  VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES(1:n_atmos_wfs) = profilewf_names(1:n_atmos_wfs)
  ! ==================== Get Ozone Absorption Cross Section ====================   
  IF (do_debug_o3p) WRITE(www_lun, '(i3,A)') num_iter, ': lidort_prof_env (start of get xcrs)'

  IF (num_iter == 0) THEN
    ! identify wavelengths for calc
    IF (allocated(do_radcals)) deallocate(do_radcals, radcal_idxs)
    allocate (do_radcals(nw), radcal_idxs(nw))
  ENDIF

  IF (nw > 1) THEN
     allcol(1, 1:nz1) = ozs(1:nz1) * du2mol
     nfgas = 1; gasin(1) = 1

     DO k = 1, ngas
        IF (fgasidxs(k) > 0) THEN
           nfgas = nfgas + 1; gasin(nfgas) = nfgas
           ! molecules cm^-2, but normalized by refspec_norm(gasidxs(k))
           ! allcol / refspec_norm will be molecules cm^-2
           allcol(nfgas, 1:nz1) = mgasprof(k, 1:nz1) * tracegas(k, 4) / mgasprof(k, nz+1)
        ENDIF
     ENDDO
 
     IF (use_effcrs) THEN
        CALL get_effres_gascrs_ray( num_iter, n_radwvl_sav, radwvl_sav, nw, &
             nz1, ts(1:nz1), ps(1:nz1), nfgas,allcol(1:nfgas,1:nz1),frhos(1:nz1), &
             do_o3shi, o3shi, do_tmpwf,do_pslwf, allcrs(1:nw,1:nfgas, 1:nz1), &
             raycof(1:nw), depol(1:nw), problems)
        abscrs(1:nw, 1:nz1) = allcrs(1:nw, 1, 1:nz1)
        IF (num_iter == 0) THEN 
          IF (do_radinter) THEN 
            do_radcals(1:nw) = .FALSE.
            fidx = 1; k = nz1 / 2

            DO iw = 1, numwin
              lidx = fidx + nradpix(iw) - 1

              ! Do radiative transfer calculations at local maxima/minima
              ! always do radiative transfer calculation at end points
              IF (band_selectors(iw) == 2) THEN
                DO i = fidx + 1, lidx - 1
                  IF (abscrs(i, k) > abscrs(i-1, k) .AND. abscrs(i, k) > abscrs(i+1, k)) &
                     do_radcals(i) = .TRUE.
                  IF (abscrs(i, k) < abscrs(i-1, k) .AND. abscrs(i, k) < abscrs(i+1, k)) &
                     do_radcals(i) = .TRUE.
                  !IF (abscrs(i, nz1) > abscrs(i-1, nz1) .AND. abscrs(i, nz1) > abscrs(i+1, nz1)  &
                  !    do_radcals(i) = .TRUE.
                  !IF (abscrs(i, nz1) < abscrs(i-1, nz1) .AND. abscrs(i, nz1) < abscrs(i+1, nz1)) &
                  !    do_radcals(i) = .TRUE.
                ENDDO

                nstep = 5; do_radcals(fidx) = .TRUE.; do_radcals(lidx-1:lidx) = .TRUE.                
              ELSE IF (band_selectors(iw) == 1) THEN
                nstep = 3; do_radcals(fidx) = .TRUE.; do_radcals(lidx) = .TRUE.
              ELSE
                nstep = 1
              ENDIF

              DO i = fidx + 1, lidx - 1, nstep
                do_radcals(i) = .TRUE.
              ENDDO

              fidx = lidx + 1        
            ENDDO

            nradcal = 0
            DO i = 1, nw
              IF (do_radcals(i)) THEN
                nradcal = nradcal + 1; radcal_idxs(nradcal) = i
              ENDIF
            ENDDO
          ELSE
            do_radcals(1:nw) = .TRUE.; nradcal = nw
            radcal_idxs(1:nw) = (/(i, i=1, nw)/)
          ENDIF
        ENDIF
     ELSE   ! .NOT. use_effcrs
       ! O3/SO2 (use_so2dtcrs=.TRUE.) cross section: if do_tmpwf = .FALSE. and do_o3shi is FALSE, 
       ! just need to get once for each retrieval
       ! Other trace gas cross section: just need to get it once for all the retrievals if no shifts 
     
       CALL GET_HRES_GASCRS_RAY(num_iter, nw, waves, & 
             nz1, ts(1:nz1), ps(1:nz1),nfgas, allcol(1:nfgas, 1:nz1),frhos(1:nz1), &
             do_o3shi, o3shi, do_tmpwf, allcrs(1:nw, 1:nfgas, 1:nz1),& 
             raycof(1:nw), depol(1:nw), problems)
      
       abscrs(1:ncalcp, 1:nz1) = allcrs(1:ncalcp, 1, 1:nz1)
       IF (num_iter == 0) THEN
         do_radcals(1:nw) = .FALSE. ; nradcal = ncalcp
         do_radcals(1:ncalcp) = .TRUE.
         radcal_idxs(1:ncalcp) = (/(i, i=1, ncalcp)/)
       ENDIF
     ENDIF
     IF (num_iter == 0 ) THEN 
       IF (use_so2dtcrs) THEN
           fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(so2idx)))= mgasprof(so2idx, nflay+1) * refspec_norm(so2_idx)
       ENDIF
       IF (use_o4dtcrs) THEN 
           fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(o4idx)))= mgasprof(o4idx, nflay+1) * refspec_norm(o2o2_idx)
       ENDIF
       IF (use_o2dptcrs) THEN 
           fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(o2idx)))= mgasprof(o2idx, nflay+1) * refspec_norm(o2_idx)
       ENDIF
       IF (use_h2odptcrs) THEN 
           fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(h2oidx)))= mgasprof(h2oidx, nflay+1) * refspec_norm(h2o_idx)
       ENDIF
     ENDIF
   ELSE
     do_radcals(1) = .TRUE.; nradcal = 1; radcal_idxs(1) = 1

     ! o3 absorption coefficient at 370.2 nm with TOMS FWHM
     IF (allocated(dads%o3)) deallocate (dads%o3)
     IF (allocated(dadt%o3)) deallocate (dadt%o3)
     allocate(dads%o3(1, nz1), dadt%o3(1, nz1))
     dads%o3(1, 1:nz1) = 0.0; dadt%o3(1, 1:nz1) = 0.0

     ! Weighted by solar flux
     !abscrs(1, 1:nz) = 9.1231787D-24 + (ts(1:nz) - 273.15) * 1.9005502D-25 + &
     !     (ts(1:nz) - 273.15)**2.0 * 1.2275286D-27     
     !raycof(1) = 2.3184501D-26; depol(1) = 0.030247913D0

     nfgas = 7
     CALL GET_ALB_OZCRS_RAY(nz1, ts(1:nz1), nfgas, allcrs(1, 1:nfgas, 1:nz1), raycof(1), depol(1), problems)
     IF (problems) THEN
        WRITE(*, *) modulename, ' : Problems in reading O3 XSec for determining Fc!!!'
        errstat = pge_errstat_error; STOP 1
     ENDIF

     !CALL GET_ALL_RAYCOF_DEPOL(1, 360.D0, raycof(1), depol(1))
     !print *, raycof(1), depol(1)

     ! Ozone
     nfgas = 6;  gasin(1) = 1
     allcol(1, 1:nz1) = ozs(1:nz1) * du2mol

     ! O3, No2, So2, bro, Hco, O4, oclo, f0
     !tmp_gasidxs = (/1, 4, 5, 7, 3/)
     tmp_gasidxs = (/no2idx, so2idx, broidx, hchoidx,o4idx/)
     DO i = 1, 4
        allcol(i+1, 1:nz1) = mgasprof(tmp_gasidxs(i), 1:nz1)
        gasin(i+1) = i + 1
     ENDDO

     ! O2-O2 Concentration (Oxygen: 20.95%)
     allcol(6, 1:nz1) = ( frhos(1:nz1) * 0.2095 ) ** 2.0 / (fzs(0:nz1-1) - fzs(1:nz1)) / 1.0D5
     gasin(6) = 6
  ENDIF

  IF (do_debug_o3p) WRITE(www_lun, '(i3,A, I3)') num_iter, ': lidort_prof_env (end of get xcrs)', errstat
  IF (problems) THEN
    WRITE(*, *) modulename, ' : Problems in reading O3 XSec !!!'
    errstat = pge_errstat_error; RETURN
  ENDIF

  ! Initialize aerosol/cloud property profiles
  IF (num_iter == 0) THEN
     aersca(1:nz1) = 0.0; aerext(1:nz1) = 0.0
     aerasy(1:nz1) = 0.0;  aermoms(0:nmom, 1:ngksec, 1:nz1) = 0.0
  ENDIF

  cldmsk = .FALSE.; IF (nctp /= 0) cldmsk(nctp:ncbp)=.TRUE.
  cldsca=0.0; cldext=0.0; cldasy=0.0;  cldmoms = 0.0

  ! Initialize output variables  
  rad = 0.0
  cfracwf = 0.0
  IF ( do_fozwf )  fozwf   = 0.0
  IF ( do_albwf)   albwf   = 0.0
  IF ( do_codwf )  fcodwf  = 0.0
  IF ( do_sprswf ) fsprswf  = 0.0
  IF ( do_fraywf ) fraywf  = 0.0
  IF ( do_ctpwf )  ctpwf   = 0.0
  IF ( do_faerwf)  faerwf  = 0.0
  IF ( do_twaewf ) faerswf = 0.0

  ! senstitivity of absorption cross section to temperature, used 
  ! for calculating temperature wf directly with LIDORT
  !eta = 0.0             ! dummy variable here
  alleta = 0.0

 !******************** polarization correction setting  *********************
  ! Determine wavelengths where exact polarization correction (NSTOKES: 4 vs 1) is calculated
  ! In UV1 (or between 270 and 310 nm): ~292 nm, ~298 nm, ~300 nm, ~302 nm, ~304 nm, ~306 nm, last wavelength 
  ! In UV2 (or between 310 and 340 nm): first, 1/4, middle and last wavelength
  ! So exact vector LIDORT calculation is done at 11 wavelengths.
  ! This option works when radiance interpolation option is turned on

  IF (num_iter == 0) THEN
    IF (allocated(do_polcorrs)) deallocate (do_polcorrs, polidx)
    allocate (do_polcorrs(nw), polidx(nw))
    npolcorr=0;do_polcorrs(1:nw) = .FALSE.;polidx(1:nw) = 0; polcorr_idxs(:) = 0
    IF ( (polcorr >= 3 .AND. polcorr <= 5) .AND. nw > 1 ) THEN
      IF (use_effcrs) THEN
       CALL set_polcorr_new(numwin, winlim(1:numwin,:),nw, waves, do_radcals(1:nw), &
                       npolcorr,  do_polcorrs(1:nw), polidx(1:nw), polcorr_idxs)
      ELSE
       CALL set_polcorr_new(numwin, winlim(1:numwin,:),ncalcp, radcwav, &
           do_radcals(1:ncalcp),npolcorr,do_polcorrs(1:ncalcp),polidx(1:ncalcp), &
           polcorr_idxs)
    ENDIF
   ENDIF
  ENDIF

  polerr = 1.0
  e_vlidort = 0.0
  CALL cpu_time(e_loop1)  
  ! ====================== Call LIDORT and Do Post Processing =====================
  !print *, 'nw, nw0, ncalcp: ', nw, nw0, ncalcp
  DO iw = 1, nw 
     
     IF (.NOT. do_radcals(iw) ) CYCLE   ! Radiances and weighting functions will be interpolated 
     lamda = waves(iw) 

     IF (nwfc > 0) the_cfrac = wfcs(iw)

     ! NSTOKES = 4 when
     ! (1) Vector LIDORT or 
     ! (2) on-line polarization correction (polcorr=3) at cloud wavelength or wavelengths 
     !     where exact radiances are calcualted or
     ! (3) on-line polarization correction (polcorr=4) at cloud wavelength or wavelengths 
     !     where exact radiances are calcualted in the first iteration
     IF (polcorr == 0 .OR. nw==1 .OR.((polcorr == 3 .OR. polcorr == 5) .AND. (nw == 1 .OR. do_polcorrs(iw)))  &
          .OR. (polcorr == 4 .AND. (nw == 1 .OR. (do_polcorrs(iw) .AND. num_iter == 0))) ) THEN
       NSTOKES = 3 ;    NSTREAMS = the_str
       IF (nw ==1) NSTREAMS = 4
     ELSE
       NSTOKES = 1 ;    NSTREAMS = the_str
     ENDIF
     IF ( lamda < 295.0 .AND. NSTOKES == 1 .AND. nw > 1 .AND. do_ssfullb295) THEN
       VLIDORT_FixIn%Bool%TS_DO_SSFULL = .TRUE.;  
       VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .FALSE.
       VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.
     ELSE
       VLIDORT_FixIn%Bool%TS_DO_SSFULL = .FALSE.
       IF (aerosol .OR. (has_clouds .AND. .NOT. do_lambcld)) THEN
         VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .TRUE.
         VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = .TRUE.
       ENDIF
     ENDIF

     IF ( aerosol ) THEN  !.AND. num_iter == 0 ) THEN
       ! Interpolate/extrapoalte for aerosol properties
       hgh = actawin
       DO i = 1, actawin
         IF (lamda < aerwavs(i)) THEN
           hgh = i; EXIT
         ENDIF
       ENDDO
       IF (hgh==1) hgh = 2    ! extrapolation
       low = hgh - 1
       xg = (lamda - aerwavs(low)) / (aerwavs(hgh) - aerwavs(low))

       aersca(faer_lvl:nfsfc-1) = gasca(low, faer_lvl:nfsfc-1) * (1.0 - xg) +  gasca(hgh, faer_lvl:nfsfc-1) * xg
       aerext(faer_lvl:nfsfc-1) = gaext(low, faer_lvl:nfsfc-1) * (1.0 - xg) +  gaext(hgh, faer_lvl:nfsfc-1) * xg

       !IF (num_iter == 0) THEN  ! Don't need to be updated
       IF (useasy) THEN
         aerasy(faer_lvl:nfsfc-1) = gaasy(low, faer_lvl:nfsfc-1) * (1.0 - xg) + &
           gaasy(hgh, faer_lvl:nfsfc-1) * xg
       ELSE        
         DO i = 0, nmom
           DO j = 1, ngksec
             aermoms(i, j, faer_lvl:nfsfc-1) = &
               gamoms(low, faer_lvl:nfsfc-1, i, j) * (1.0 - xg) +  &
               gamoms(hgh, faer_lvl:nfsfc-1, i, j) * xg
            ENDDO
          ENDDO
       ENDIF
        !ENDIF
     ENDIF

     IF (has_clouds .AND. .NOT. do_lambcld) THEN
       cldext(nctp:ncbp) = (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) * &
             (fzs(nctp-1 : ncbp-1) - fzs(nctp:ncbp)) * the_cbeta
       cldsca(nctp:ncbp) = cldext(nctp:ncbp)  ! w0 = 1.0
       IF (iw == 1) cldext0(nctp:ncbp) = cldext(nctp:ncbp) * gcq(actawin) / (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) 

       IF (useasy) THEN
         cldasy(nctp:ncbp) = gcasy(low) * (1.0 - xg) + gcasy(hgh) * xg
       ELSE
         DO i = 0, nmom
           DO j = 1, ngksec
             cldmoms(i, j, nctp:ncbp) = gcmoms(low, i, j) * (1.0 - xg) + gcmoms(hgh, i, j) * xg 
           ENDDO
         ENDDO
       ENDIF
     ENDIF

     IF ( do_polcorrs(iw) .AND. ((polcorr == 3 .OR. polcorr == 5) &
          .OR. (polcorr == 4 .AND. num_iter == 0) ) ) THEN
        npolmod = 2   ! Twice, one vector and one scalar          
     ELSE 
        npolmod = 1   ! Only once either scalar or vector
     ENDIF

     ! When polcorr = 5 is selected, only calculate weighting function for iteration 1
     ! iteration 0 if 1st pixel of a x-track position is being retrieved
     IF (polcorr == 5 .AND. nw > 1 .AND. (num_iter > 1 .OR. (num_iter == 0 .AND. currloop /= 0))) THEN
        VLIDORT_LinModIn%MCont%TS_do_simulation_only = .TRUE.
        VLIDORT_LinModIn%MCont%TS_do_linearization   = .FALSE.
     ENDIF
        
     DO ipol = 1, npolmod
       IF (ipol == 2) THEN
         NSTOKES = 1 ;   NSTREAMS = the_str ! Always scalar for second mode

         IF ( lamda < 295.0 .AND. nw > 1 .AND. do_ssfullb295) THEN
           VLIDORT_FixIn%Bool%TS_DO_SSFULL = .TRUE.
           VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .FALSE.
           VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.
         ELSE
           VLIDORT_FixIn%Bool%TS_DO_SSFULL = .FALSE.
           VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .TRUE.
           VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = .TRUE.
          !DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
         ENDIF
         IF (.NOT. aerosol) THEN
           VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .FALSE.
           VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.
         ENDIF
              
        ! Save VECTOR LIDORT results
        prad(polidx(iw)) = rad(iw, 1)
        rad(iw, 1:nostk) = 0.d0
        IF (do_cfracwf) pcfracwf(polidx(iw)) = cfracwf(iw, 1)
        IF (do_cfracwf) cfracwf(iw, 1:nostk) =  0.d0

        IF ( VLIDORT_LinModIn%MCont%TS_do_linearization ) THEN
          IF (do_albwf)  palbwf(polidx(iw), 1:nalbwf)   = albwf(iw,1:nalbwf, 1)
          IF (do_fozwf)  pfozwf(polidx(iw), :)          = fozwf(iw, :, 1)
          IF (do_fozwf)  ptauwf(polidx(iw), :)          = tauwf(iw, :)
          IF (do_codwf)  pfcodwf(polidx(iw), nctp:ncbp) = fcodwf(iw, nctp:ncbp, 1)
          IF (do_sprswf) pfsprswf(polidx(iw), nup2p(nsfc-1)+1:nfsfc-1) = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1)
          IF (do_fraywf) pfraywf(polidx(iw), :)         = fraywf(iw, :, 1)
          IF (do_ctpwf)  pctpwf(polidx(iw))             = ctpwf(iw, 1)
          IF (do_faerwf) pfaerwf(polidx(iw), faer_lvl:nfsfc-1)  = faerwf(iw, faer_lvl:nfsfc-1, 1)           
          IF (do_twaewf) pfaerswf(polidx(iw), faer_lvl:nfsfc-1) = faerwf(iw, faer_lvl:nfsfc-1, 1)

          ! Initalize those variables to zero again          
          IF (do_albwf) albwf(iw, :, 1:nostk) = 0.d0
          IF (do_fozwf) fozwf(iw, :, 1:nostk) = 0.d0 
          IF (do_codwf) fcodwf(iw, nctp:ncbp, 1:nostk) = 0.d0
          IF (do_sprswf) fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1:nostk) = 0.d0
          IF (do_fraywf) fraywf(iw, :, 1:nostk) = 0.d0
          IF (do_ctpwf)  ctpwf(iw, 1:nostk)       = 0.d0
          IF (do_faerwf) faerwf(iw, faer_lvl:nfsfc-1, 1:nostk)  = 0.d0
          IF (do_twaewf) faerswf(iw, faer_lvl:nfsfc-1, 1:nostk) = 0.d0
        ENDIF
      ENDIF

      radclrcld = 0.0
      DO ic = 1, 2 ! for clear and cloud
        IF (ic == 1) THEN
          do_clouds = .FALSE.; frac = 1.0 - the_cfrac
        ELSE
          do_clouds = .TRUE. ; frac = the_cfrac
        ENDIF
        IF (frac == 0.0) CYCLE  ! No clear/cloudy part
        ! Note convert ozone from DU to molecule/cm^2  here

        IF ((ic == 1) .OR. (.NOT. do_lambcld))  THEN
          ! Reset up number of layers since Lambertian 
          ! cloudy scene got different layers
          nlayers = nfsfc - 1; nz1 = nlayers   
          ! zero O3 weighting function below surface (This is necessary)           
          lambertian_albedo = albs(iw)
        ELSE
          nlayers = nctp-1  ! from cloud top to TOA
          nz1 = nlayers;    do_clouds = .FALSE.
          IF (the_cfrac == 1.0 .AND. nw /= 1) THEN
            lambertian_albedo = albs(iw)
            lambcld_refl = albs(iw)
          ELSE
            lambertian_albedo = lambcld_refl 
            ! use 80% (could be adjusted when the_cfrac gt 0.90)
          ENDIF
        ENDIF
        IF (polcorr==2) THEN
          tmpalb = lambertian_albedo
          IF (tmpalb <0.001) tmpalb = 0.001
          IF (tmpalb >0.999) tmpalb = 0.999
          albclrcld(iw, ic) = tmpalb
        ENDIF

        IF (ipol == 1) THEN
           problems = .FALSE.
            IF (ic == 1 .or. (ic == 2 .and. frac == 1.0)) THEN 
            CALL LIDORT_PROF_PREP(lamda, raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
              varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
              alleta(1:nfgas, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
              aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
              do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
              cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems, &
              deltau(iw, 1:nz1), delsca(iw, 1:nz1), delray(iw, 1:nz1))
            ENDIF
            IF (problems) THEN
              WRITE(*, *) modulename, ' : Problems encountered in lidort preparation!!!'
              errstat = pge_errstat_error; RETURN
            END IF
        ENDIF
                  
        VLIDORT_FixIn%Cont%TS_NSTOKES  = NSTOKES
        VLIDORT_FixIn%Cont%TS_NSTREAMS = NSTREAMS
        VLIDORT_FixIn%Cont%TS_NLAYERS  = NLAYERS
        VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo

        IF (ic == 1 .and.  lamda > wcenter_uvvis .and. do_brdf .and. the_cfrac < 1) THEN
          Surface%nstokes = NSTOKES ;  Surface%nstreams = NSTREAMS
          CALL Surface%SetOpticalProperties_VL(iw, Window%RTM_Wvl(iw), &
            Geolocation, I0_conversionFactor, SurfProf, &
            VLIDORT_FixIn, VLIDORT_LinModIn, VLIDORT_Sup, VLIDORT_LinSup,Error)
          CALL Surface%SetLinearization (.TRUE., .FALSE., Error)
          VLIDORT_FIxIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .FALSE.
          n_sfc_wfs = Surface%njac
        ELSE
          VLIDORT_FIxIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .TRUE.
          n_sfc_wfs = 1
        ENDIF

        VLIDORT_LinFixIn%Cont%TS_n_surface_wfs = n_sfc_wfs
        VLIDORT_LinOut%Prof%TS_profilewf(1:n_atmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0
        CALL cpu_time (e1)
        CALL VLIDORT_LPS_master ( &
          do_debug_input, &
          VLIDORT_FixIn,    & ! INPUTS
          VLIDORT_ModIn,    & ! INPUTS (possibly modified)
          VLIDORT_Sup,      & ! INPUTS/OUTPUTS
          VLIDORT_Out,      & ! OUTPUTS
          VLIDORT_LinFixIn, & ! INPUTS
          VLIDORT_LinModIn, & ! INPUTS (possibly modified)
          VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
          VLIDORT_LinOut )    ! OUTPUTS

        VLIDORT_LinOut%Prof%TS_profilewf(1:n_atmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0
        IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_SERIOUS ) THEN
          WRITE(*,*)'VLIDORT input abort, PROFILEWF calculation', lamda, iw ; STOP 1
        ELSE IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .ne. VLIDORT_SERIOUS .and. &
              VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
              print * ,  VLIDORT_Out%Status%TS_STATUS_INPUTCHECK , &
              VLIDORT_Out%Status%TS_STATUS_CALCULATION 
          WRITE(*,*)'VLIDORT calculation abort, PROFILEWF calculation',lamda,iw!; STOP 1 
        ENDIF
        CALL cpu_time(e2)
        e_vlidort = e_vlidort + e2 - e1

        ! Pixel-independent approximation
        radclrcld(ic, 1:nostk) = VLIDORT_Out%Main%TS_stokes(1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk)
        rad(iw, 1:nostk)       = rad(iw, 1:nostk) + radclrcld(ic, 1:nostk) * frac
        IF (ic == 1 .and. do_brdf .and. lamda > 649 )  THEN 
          WRITE(*,'(f7.2, 20e14.6)')  lamda, rad(iw, 1),&
           vlidort_sup%brdf%ts_exactdb_brdfunc (1,1, 1, 1), &
           vlidort_sup%brdf%ts_brdf_f_0(1,1, 1,1), &
           vlidort_sup%brdf%ts_brdf_f (0,1,1, 1), &
           vlidort_sup%brdf%ts_user_brdf_f_0 (0,1,1, 1), &
           vlidort_sup%brdf%ts_user_brdf_f (0,1,1, 1), &
          VLIDORT_Out%Main%TS_stokes(1, 1, 1, 1),&
          VLIDORT_LinOut%Surf%TS_surfacewf(1:nalbwf, 1, 1, 1, 1) !, surface%kern_amp(iw, 1:4)
        ENDIF
        IF (VLIDORT_LinModIn%MCont%TS_do_linearization ) THEN
          ! weighting function per Dobson Unit
          IF ( do_ozwf ) THEN
            DO istk = 1, nostk
              fozwf(iw, 1:nfsfc-1, istk) = fozwf(iw, 1:nfsfc-1, istk) + &
                VLIDORT_LinOut%Prof%TS_profilewf(ozwfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                / ozs(1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
            ENDDO
          ENDIF
          ! Weighting function with respect to aerosol/cloud optical depth at the last wavelength
          ! so as to keep the scaling since we are fitting the aod at that wavelength
          IF ( do_taodwf .OR. do_saodwf ) THEN
            DO istk = 1, nostk
              faerwf(iw, faer_lvl:nfsfc-1, istk)  = faerwf(iw, faer_lvl:nfsfc-1, istk) + &
                VLIDORT_LinOut%Prof%TS_profilewf(aodwfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                gaext(actawin, faer_lvl:nfsfc-1) * polerr(iw, ic, istk) *  frac
            ENDDO
          ENDIF
          IF ( do_twaewf ) THEN
            DO istk = 1, nostk
              faerswf(iw, faer_lvl:nfsfc-1, istk)  = faerswf(iw, faer_lvl:nfsfc-1, istk) + &
                VLIDORT_LinOut%Prof%TS_profilewf(twaewfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                gasca(actawin, faer_lvl:nfsfc-1) * gaext(actawin, faer_lvl:nfsfc-1) * polerr(iw, ic, istk) *  frac
            ENDDO
          ENDIF
          IF ( do_codwf ) THEN
            DO istk = 1, nostk
              fcodwf(iw, nctp:ncbp, istk)  = fcodwf(iw, nctp:ncbp, istk) + &
                VLIDORT_LinOut%Prof%TS_profilewf(codwfidx, nctp:ncbp, 1, 1, istk, 1) / &
                cldext0(nctp:ncbp) * polerr(iw, ic, istk) *  frac
            ENDDO
          ENDIF
          IF ( do_sprswf ) THEN
            DO istk = 1, nostk 
              fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk)  = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk) + &
                VLIDORT_LinOut%Prof%TS_profilewf(sprswfidx, nup2p(nsfc-1)+1:nfsfc-1, 1, 1, istk, 1) / &
                (fps(nup2p(nsfc-1)+1:nfsfc-1) - fps(nup2p(nsfc-1):nfsfc-2)) * polerr(iw, ic, istk) *  frac
            ENDDO
          ENDIF
          IF ( do_fraywf ) THEN
            DO istk = 1, nostk
              fraywf(iw, 1:nfsfc-1, istk) = fraywf(iw, 1:nfsfc-1, istk) + &
                VLIDORT_LinOut%Prof%TS_profilewf(raywfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                / delray(iw, 1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
            ENDDO
          ENDIF

          ! Non-lambertian clouds: albedo wf is from both clear and cloud
          ! Lamberitan clouds:     if the_cfrac  < 1, albwf from clear only
          !                        if the_cfrac == 1, albwf from cloud only
          ! xliu, 03/08/11: the surface albedo weighting functions from v2p4RTC are un-normalized.
          ! Also, we need to initialize n_surface_wfs = 1 for calculating lambertian surface weighting function
          IF (do_albwf) THEN 
            DO ialbwf = 1, VLIDORT_LinFixIn%Cont%TS_n_surface_wfs
              IF (.NOT. do_lambcld) THEN
                albwf(iw, ialbwf, 1:nostk) = albwf(iw, ialbwf, 1:nostk) + &
                         VLIDORT_LinOut%Surf%TS_surfacewf(ialbwf, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac 
              ELSE IF (the_cfrac < 1.0) THEN
                IF (ic == 1) albwf(iw, ialbwf, 1:nostk) = albwf(iw, ialbwf, 1:nostk) + &
                         VLIDORT_LinOut%Surf%TS_surfacewf(ialbwf, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac 
              ELSE IF (the_cfrac == 1.0) THEN
                IF (ic == 2) albwf(iw, ialbwf, 1:nostk) = albwf(iw, ialbwf, 1:nostk) + &
                         VLIDORT_LinOut%Surf%TS_surfacewf(ialbwf, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac
              ENDIF
              IF (.NOT. VLIDORT_FIxIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN 
                ! convert dI/da => dI/dc
                albwf(iw, ialbwf, 1:nostk) = albwf(iw, ialbwf, 1:nostk) *surface%kern_amp(iw, ialbwf)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO ! end clear/cloudy scene loop
      IF (do_cfracwf) cfracwf(iw, 1:nostk) = radclrcld(2, 1:nostk) - radclrcld(1, 1:nostk)
    ENDDO    ! end scalar/vector modes
  ENDDO       ! end wavelength loop 
  CALL cpu_time (e_loop2)

  IF (do_debug_o3p) WRITE(www_lun, '(i3,A, I3)') num_iter, ': lidort_prof_env (do iw = 1, nw)'
  !----------------------------------
  ! Post-pocessing after RTM simulations
  !-----------------------------------

  IF ( the_cfrac == 1.0 .AND. do_lambcld) THEN
     nz1 = nctp - 1
  ELSE
     nz1 = nfsfc - 1
  ENDIF
  nsprs = nup2p(nsfc-1)+1
  IF (nw == 1 .or. use_effcrs) THEN
    nw0 = nw 
  ELSE
    nw0 = ncalcp
  ENDIF

  IF (num_iter >=0 .and.  nw > 1 .and. do_debug_rtm ) THEN ! before polcorr
    rto%nw = nw0 ; rto%nl = nz1 ; rto%nos = nostk ; rto%nalbwf = nalbwf
    CALL set_rtmvar(.TRUE., rto)
    IF (use_effcrs) THEN 
      rto%wav(1:nw0)   = fitwavs(1:nw0)
    ELSE
      rto%wav(1:nw0)   = radcwav(1:nw0)
    ENDIF   
    rto%alb(1:nw0)   = albs(1:nw0)
    rto%rad(1:nw0,1:nostk)   = rad(1:nw0, 1:nostk)
    rto%albwf(1:nw0,1:nalbwf,1:nostk) = albwf(1:nw0, 1:nalbwf,1:nostk)
    rto%ozwf(1:nw0,1:nz1,1:nostk)    = fozwf(1:nw0,1:nz1,1:nostk)
    IF (use_effcrs) THEN
      CALL debug_rtm (10,rto, ccrs)
    ELSE
      CALL debug_rtm (11,rto, hcrs)
    ENDIF
    CALL set_rtmvar(.FALSE., rto)
  ENDIF

  !----------------------------------
  ! On-line polarization correction
  !-----------------------------------
  do_plutcorr_after = .false.
  IF (polcorr == 2 .and. nw > 1 .and. .NOT. do_plutcorr_after) THEN
    fidx = 1 ; lidx = nw0
    !do_ssfullb295 = .true.
    !IF (do_ssfullb295) fidx = MINVAL(MINLOC(waves(1:nw0), MASK=waves(1:nw0) >
    !280))
    delabs(fidx:lidx, 1:nz1) = deltau(fidx:lidx, 1:nz1) - delsca(fidx:lidx,1:nz1)
    DO i = 1, nz1
     tauwf (fidx:lidx, i) = fozwf(fidx:lidx, i, 1)/abscrs(fidx:lidx,i)/rad(fidx:lidx,1)/du2mol
    ENDDO

    toz = SUM(ozs(1:nfsfc-1))
    ! Junsung: revised the ots library for VLDLUTdir
    ! Junsung: chnage the nfsfc to nsprs
    ! VLDLUTdir='/home/jbak/data/GEMSTOOL/lutdatav2.8-r/LUT-48/'
    ! VLDLUTdir='/home/jbak/OzoneFit/tbl/vldlut/'
    VLDLUTdir= TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(vldlut))
    CALL polcorr_online_with_lut(num_iter,VLDLUTdir,lidx-fidx+1, nz1, nctp,nsprs, nalbwf,&
         do_albwf, do_cfracwf, the_cfrac, albclrcld(fidx:lidx, 1:2),&
         sza, vza, aza, lat, toz, fps(0:nz1), ps(1:nz1), &
         waves(fidx:lidx),  rad(fidx:lidx, 1),&
         tauwf(fidx:lidx, 1:nz1),delabs(fidx:lidx, 1:nz1), delsca(fidx:lidx,1:nz1), &
         albwf(fidx:lidx, 1:nalbwf,1),fozwf(fidx:lidx, 1:nz1, 1),  &
         cfracwf(fidx:lidx, 1))
  ENDIF
  IF ( (polcorr >= 3 .AND. polcorr <= 5 ) .AND. nw > 1) THEN
    delabs(1:nw0, 1:nz1) = deltau(1:nw0, 1:nz1) - delsca(1:nw0, 1:nz1)
    DO i = 1, nz1
       tauwf (1:nw0, i) = fozwf(1:nw0, i, 1)/abscrs(1:nw0, i)/rad(1:nw0,1)/du2mol
       ptauwf (1:npolcorr, i) = pfozwf(1:npolcorr, i)/abscrs(polcorr_idxs(1:npolcorr), i) & 
                                /prad(1:npolcorr)/du2mol 
    ENDDO
  
    ! jbak, dimension for p-variables changed from nw to npolcorr
    CALL polcorr_online(num_iter, polcorr, nw0, nz1, nctp, ncbp, nsprs, nalbwf,faer_lvl,&
       npolcorr, polcorr_idxs(1:npolcorr), do_fozwf, do_albwf, do_faerwf,            &
       do_twaewf, do_codwf, do_sprswf, do_fraywf, do_cfracwf, & 
       waves(1:nw0),  rad(1:nw0, 1), waves(polcorr_idxs(1:npolcorr)), prad(1:npolcorr), & 
       tauwf(1:nw0, 1:nz1), ptauwf(1:npolcorr, 1:nz1), delabs(1:nw0, 1:nz1),& 
       albwf(1:nw0,1, 1:nalbwf),palbwf(1:npolcorr, 1:nalbwf),fozwf(1:nw0, 1:nz1, 1),pfozwf(1:npolcorr,   1:nz1), & 
       faerwf(1:nw0, 1:nz1, 1), pfaerwf(1:npolcorr,1:nz1), faerswf(1:nw0, 1:nz1, 1),pfaerswf(1:npolcorr, 1:nz1),&
       fcodwf(1:nw0, 1:nz1, 1), pfcodwf(1:npolcorr, 1:nz1),fsprswf(1:nw0, 1:nz1, 1),pfsprswf(1:npolcorr, 1:nz1), & 
       fraywf(1:nw0, 1:nz1, 1), pfraywf(1:npolcorr, 1:nz1),cfracwf(1:nw0, 1),       pcfracwf(1:npolcorr))

  ENDIF
  IF (do_debug_o3p) WRITE(www_lun, '(i3,A, I3)') num_iter, ': lidort_prof_env (end of on-line pol corr)'
  CALL cpu_time(e_pol)

  IF (  nw > 1 .and. do_debug_rtm) THEN ! after polcorr
    rto%nw = nw0 ; rto%nl = nz1 ; rto%nos = nostk ;rto%nalbwf=nalbwf
    CALL set_rtmvar(.TRUE., rto)
    IF (use_effcrs) THEN 
      rto%wav(1:nw0)   = fitwavs(1:nw0)
    ELSE
      rto%wav(1:nw0)   = radcwav(1:nw0)
    ENDIF
    rto%alb(1:nw0)   = albs(1:nw0)
    rto%rad(1:nw0,1:nostk)   = rad(1:nw0, 1:nostk)
    rto%albwf(1:nw0,1:nalbwf, 1:nostk) = albwf(1:nw0, 1:nalbwf, 1:nostk)
    rto%ozwf(1:nw0,1:nz1,1:nostk)    = fozwf(1:nw0,1:nz1,1:nostk)
    IF (use_effcrs) THEN   
      CALL debug_rtm (10,rto, ccrs)
      CALL debug_taug(100, nw0, nz1, nfgas, rto%wav(1:nw0), allcol(1:nfgas,1:nz1),allcrs(1:nw0, 1:nfgas, 1:nz1))
    ELSE
      CALL debug_rtm (11,rto, hcrs)
    ENDIF
    CALL set_rtmvar(.FALSE., rto)
  ENDIF

  IF (do_brdf .and. do_albwf .and. nw > 1) THEN  
    ! convert dI/da => dI/dc
    idx = MINVAL(MAXLOC(radcwav(1:ncalcp), MASK=(radcwav(1:ncalcp) < wcenter_uvvis)))+1
    DO i = 1, nostk
    DO iwf = 1, nalbwf          
        albwf(idx:nw0, iwf, i) = albwf(idx:nw0, iwf, i) *surface%kern_amp(idx:nw0, iwf)
    ENDDO
    ENDDO
  ENDIF

  !-----------------------------------------------------------------
  ! convert simulations to higher resolution into instrument resoution
  !-----------------------------------------------------------------
  IF (do_debug_o3p) WRITE(www_lun, '(i3,A, I3)') num_iter, ': lidort_prof_env (radwf_interpol)'
  IF (nw > 1 .AND. do_radinter ) THEN
     CALL radwf_interpol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, do_radcals(1:nw),             &
          do_fozwf, do_albwf, do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf,         &
          waves, abscrs(1:nw, 1:nz1),ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1),        &
          albwf(1:nw, 1, 1), cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
          fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), errstat)
     IF (errstat == pge_errstat_error) RETURN
  ENDIF

  IF (nw > 1 .AND. .NOT. use_effcrs) THEN
    IF (do_fraywf .AND. do_fozwf) THEN 
      do_abs = .false.
      if (polcorr==2 .and. do_plutcorr_after) do_abs = .true.
      delabs(1:nw, 1:nz1) = deltau(1:nw, 1:nz1) - delsca(1:nw, 1:nz1)
      CALL hres_radwf_inter_convol(nw, nz1, nctp, ncbp, nsprs, nalbwf, faer_lvl,  & 
      do_albwf, do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf, do_tracewf,&
      do_o3shi, do_tmpwf,do_pslwf, waves, ozs(1:nz1), do_abs, delabs(1:nw,1:nz1), &
      rad(1:nw, 1),fozwf(1:nw, 1:nz1, 1), albwf(1:nw,1:nalbwf, 1), cfracwf(1:nw, 1), &
      faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
      fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), fraywf(1:nw,1:nz1, 1), errstat)
      IF (errstat == pge_errstat_error) RETURN
      IF (polcorr == 2  .and. do_plutcorr_after) THEN
        fidx = 1 ; lidx = n_radwvl_sav
        DO i = 1, nz1
         tauwf (fidx:lidx, i) = fozwf(fidx:lidx, i, 1)/ccrs%o3(fidx:lidx,i)/rad(fidx:lidx,1)/du2mol
        ENDDO
        ! Junsung: revised the ots library for VLDLUTdir
        ! Junsung: chnage the nfsfc to nsprs
        ! VLDLUTdir='/home/jbak/data/GEMSTOOL/lutdatav2.8/LUT-conv/'
        ! VLDLUTdir='/home/jbak/OzoneFit/tbl/vldlut/'
        VLDLUTdir= TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(vldlut))
        CALL polcorr_online_with_lut(num_iter,VLDLUTdir,lidx-fidx+1, nz1, nctp,nsprs, nalbwf,&
             do_albwf, do_cfracwf, the_cfrac, albclrcld(fidx:lidx, 1:2),&
             sza, vza, aza, lat, toz, fps(0:nz1), ps(1:nz1), &
             radwvl_sav(fidx:lidx),  rad(fidx:lidx, 1),&
             tauwf(fidx:lidx, 1:nz1),delabs(fidx:lidx, 1:nz1), delsca(fidx:lidx,1:nz1),&
             albwf(fidx:lidx, 1:nalbwf,1),fozwf(fidx:lidx, 1:nz1, 1),  &
             cfracwf(fidx:lidx, 1))
      ENDIF
    ELSE
      WRITE(*, *) 'Must have O3/rayleigh weighting function to do radiance interpolation!!!'
      errstat = pge_errstat_error; RETURN
    ENDIF
  ENDIF

  CALL cpu_time(e_inter)
  !-----------------------------------------------------------------
  ! Calculate desired weighting functions
  ! 1. convert profilewf at fps to at umkp
  ! 2. derive tracegas weighting function
  !----------------------------------------------------------------
  IF (nw > 1 .AND. .NOT. use_effcrs) THEN
    nw0 = n_rad_wvl
  ELSE
    nw0 = nw
  ENDIF

  IF ( do_ozwf ) THEN  
    DO i = 1, nl
      fidx = nup2p(i - 1) + 1; lidx = nup2p(i)         
      DO iw = 1, nw0
        DO istk = 1, nostk
          ozwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx)) / &
                              SUM(ozs(fidx:lidx))
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF ( do_tmpwf ) THEN  
    DO i = 1, nl
      fidx = nup2p(i - 1) + 1; lidx = nup2p(i)         
      DO iw = 1, nw0
        DO istk = 1, nostk
          tmpwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx) &
            * dadt%o3(iw, fidx:lidx) * (fzs(fidx-1:lidx-1)-fzs(fidx:lidx))) / &
              (fzs(fidx-1) - fzs(lidx))
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (do_o3shi) THEN
    DO iw = 1, nw0
      DO istk = 1, nostk
        o3shiwf(iw, istk) = SUM(fozwf(iw, 1:nz1, istk) * ozs(1:nz1) * dads%o3(iw, 1:nz1))
      ENDDO
    ENDDO
  ENDIF


  IF (do_taodwf) THEN
    DO iw = 1, nw0
      DO istk = 1, nostk
        taodwf(iw, istk) = SUM( faerwf(iw, nup2p(ntp)+1:nz1, istk) * &
          gaext(actawin, nup2p(ntp)+1:nz1) ) / SUM (gaext(actawin, nup2p(ntp)+1:nz1))
      ENDDO
    ENDDO
  ENDIF

  IF (do_twaewf) THEN
    DO iw = 1, nw0
      DO istk = 1, nostk
        twaewf(iw, istk) = SUM( faerswf(iw, nup2p(ntp)+1:nz1, istk) * &
          gasca(actawin, nup2p(ntp)+1:nz1) / gaext(actawin, nup2p(ntp)+1:nz1) ) &
          / SUM (gasca(actawin, nup2p(ntp)+1:nz1) / gaext(actawin, nup2p(ntp)+1:nz1))
      ENDDO
    ENDDO
  ENDIF

  IF (do_saodwf) THEN
    DO iw = 1, nw0
      DO istk = 1, nostk
        saodwf(iw, istk) = SUM( faerwf(iw, 1:nup2p(ntp), istk) * &
          gaext(actawin, 1:nup2p(ntp)) ) / SUM (gaext(actawin, 1:nup2p(ntp)))
      ENDDO
    ENDDO
  ENDIF

  IF (do_codwf) THEN
    DO iw = 1, nw0
      DO istk = 1, nostk
        codwf(iw, istk) = SUM( fcodwf(iw, nctp:ncbp, istk) * &
          cldext0(nctp:ncbp) ) / SUM (cldext0(nctp:ncbp))
      ENDDO
    ENDDO
  ENDIF

  IF (do_sprswf) THEN
    DO iw = 1, nw0
      DO istk = 1, nostk
        sprswf(iw, istk) = SUM( fsprswf(iw, nup2p(nsfc-1)+1:nz1, istk) * &
                (fps(nup2p(nsfc-1)+1:nz1) - fps(nup2p(nsfc-1):nz1-1))) / (fps(nz1)-fps(nup2p(nsfc-1)))
      ENDDO
    ENDDO
  ENDIF

  IF (nw > 1 .AND. do_ozwf .AND. do_tracewf ) THEN
    CALL GET_TRACEGAS_WF (nw0, nz, nz1, rad(1:nw0, 1), & 
          fozwf(1:nw0, 1:nz, 1), do_so2zwf, so2zwf(1:nw0, 1))
  ENDIF

  IF (do_pslwf .and. use_effcrs) THEN
    DO k = 1, npsl
    DO iw = 1, nw0
      database_pslwf(refidx(iw), k) = &
      SUM(fozwf(iw, 1:nz1, 1) * ozs(1:nz1) * dadp%o3(iw, 1:nz1, k))

      !database_pslwf(refidx(iw), K) = database_pslwf(refidx(i), k)  + &
      ! SUM(trace_profwf(o4idx,iw, 1:nz1) *refspec_norm(o2o2_idx)* dadp%o4(iw, 1:nz1, k))
      database_pslwf(refidx(iw),k) = database_pslwf(refidx(iw),k )/rad(iw,1)
        !      print * ,  dadp%o3(iw, 1:nz1, k), 'here'
    ENDDO
    ENDDO
  ENDIF

  !------------------------------------
  ! copying local variables to output variables
  !------------------------------------
  rto%nw = nw0 ; rto%nl = nl ; rto%nos = nostk; rto%nalbwf=nalbwf
  CALL set_rtmvar(.TRUE., rto)
  rto%rad   = rad(1:nw0, :)
  rto%albwf   = albwf(1:nw0,:, :)
  rto%o3shiwf = o3shiwf(1:nw0,:)
  rto%cfracwf = cfracwf(1:nw0,:)
  rto%codwf   = codwf(1:nw0,:)
  rto%ctpwf   = ctpwf(1:nw0,:)
  rto%taodwf  = taodwf(1:nw0,:)
  rto%twaewf  = twaewf(1:nw0,:)
  rto%saodwf  = saodwf(1:nw0,:)
  rto%sprswf  = sprswf(1:nw0,:)
  rto%so2zwf  = so2zwf(1:nw0,:)
  rto%ozwf    = ozwf(1:nw0,1:nl,:)
  rto%tmpwf   = tmpwf(1:nw0,1:nl,:)
!  RETURN
  !--------------------------------
  ! Check running time in RT simulation
  !--------------------------------
  IF (do_debug_rtm .and. scnwrt .and. nw > 1 ) THEN
  call cpu_time(e_n)
  PRINT * , 'N of simulation wavelength-------:', ncalcp
  PRINT * , 'N of OMI wavelength--------------:', nw0
  print * , ' START----------------END (EXACT):', e_n     - e_s
  print * , ' CALL read_config, call read_xcrs:', e_loop1 - e_s
  print * , ' RTM LOOP iw = 1, nw-------------:', e_loop2-e_loop1
  PRINT * , '     + CALL VLIDORT--------------:', e_vlidort
  PRINT * , '     - CALL VLIDORT--------------:', e_loop2-e_loop1- e_vlidort
  print * , ' After RTM calculation-----------:', e_n-e_loop2
  print * , '     *polcorr_online ------------:', e_pol-e_loop2
  print * , '     *interpolation/convolution--:', e_inter-e_pol
  PRINT * , '++++++++++++++++++++++++++++++++++++++++++++++++++'
  ENDIF
  RETURN
END SUBROUTINE LIDORT_PROF_ENV

! ====================================================================
!  Author: Xiong Liu
!  Date: Jan. 30, 2004
!  Purpose: Prepare input (tau, single scattering albedo, phase
!		moments, their variation, surface BDRF, spherical attenutaion
!		for LIDORT calculation
! ====================================================================
! Reminder:
! Lambertian Surface model is implemented
! Refractive atmosphere is not implemented
! Modication History:
! ====================================================================

! Description of Auguments
! lamda:   wavelength 
! zsgrid:  altitude at each level from TOA to BOS in km, nlayers+1 levels
! airgrid: air column density for each layer, molecules / cm^2 
! varyprof:arrays of linearization flags for each layer, 0:no 1:yes

! Currently, number of gases to be allowed is one, i.e., O3
! Number of gases possibly to be considered later includes
! 1: O3 2: NO2 3: O4 4: BrO 5: SO2 6: HCHO  7: OCLO  8: O2  9:H2O
! These species will be added when needed, just need to add corresponding
! cross section database and read them
 
! ngas  :  number of gases
! gasin :  pointer to gases that are used
! abscrs:  Input/output
!		   If get_crs is set, then it refers to absorption cross section 
!		   at each layer and for each species for each molecule 
!          On return, it gives the absorption od for each species at each layer
! gascol : column density for each species at each layer, molecules/cm^2

! useasy    :  use asymmetric factor/phase moments for clouds/aerosols
! nmoms     :  number of phase moments
! do_aerosols: include aerosols
! aersca    :  aerosol scattering coefficients at each layer
! aerext    :  aerosol extinction coefficients at each layer
! aerasy    :  aerosol asymmetric factor  at each layer
! aermsk    :  Indicator of aerosols for each layer, 1: with aerosol, 0: no aerosols
! aermoms   :  Aerosol moments at layers with aerosols

! do_clouds : include clouds 
! cldsca    :  cloud scattering coefficients at each layer
! cldext    :  cloud extinction coefficients at each layer
! cldasy    :  cloud asymmetric factor  at each layer
! cldmsk    :  Indicator of clouds for each layer, 1: with cloud, 0: no clouds
! cldmoms   :  cloud moments at layers with clouds

SUBROUTINE LIDORT_PROF_PREP (lamda, raycof, depol, &
     zsgrid, airgrid,varyprof, &
     ngas, gasin, abscrs, gascol, eta, useasy, nmoms, &          
     do_aerosols, aersca, aerext, aerasy, aermoms, aermsk, &
     do_clouds, cldsca, cldext, cldasy, cldmoms, cldmsk, problems, &
     deltau, delsca,  delray)

  USE OMSAO_precision_module
  USE ozprof_data_module, ONLY : maxgksec, maxgkmatc, ngksec, ngkmatc
  USE OMSAO_errstat_module
  IMPLICIT NONE  
   
  !===============================  Define Variables ===========================
  ! Include files of dimensions and numbers
  !INCLUDE 'VLIDORT.PARS'  
  ! Include files of input variables
  !INCLUDE 'VLIDORT_INPUTS.VARS'
  !INCLUDE 'VLIDORT_SETUPS.VARS'
  !INCLUDE 'VLIDORT_L_INPUTS.VARS'
  !INCLUDE 'VLIDORT_BOOKKEEP.VARS'

  ! Input variables
  INTEGER, INTENT(IN)                     :: ngas, nmoms
  INTEGER, INTENT(IN), DIMENSION(ngas)    :: gasin
  LOGICAL, INTENT(IN), DIMENSION(nlayers) :: cldmsk, aermsk, varyprof
  LOGICAL, INTENT(IN)                     :: useasy
  REAL (KIND=dp), INTENT(IN)              :: raycof, depol, lamda
  REAL (KIND=dp), DIMENSION(0:nlayers), INTENT(IN) :: zsgrid
  REAL (KIND=dp), DIMENSION(nlayers), INTENT(IN)   :: airgrid,  &
       aersca, aerext, aerasy, cldsca, cldext, cldasy
  REAL (KIND=dp), DIMENSION(0:nmoms, maxgksec, nlayers), INTENT(IN) :: aermoms
  REAL (KIND=dp), DIMENSION(0:nmoms, maxgksec, nlayers), INTENT(IN) :: cldmoms
  REAL (KIND=dp), DIMENSION(ngas, nlayers), INTENT(INOUT)  :: abscrs
  REAL (KIND=dp), DIMENSION(ngas, nlayers), INTENT(IN)  :: gascol, eta
  ! Optional output
  REAL (KIND=dp), DIMENSION(nlayers), INTENT(OUT) :: deltau, delsca, delray
  
  
  ! Output variables
  LOGICAL, INTENT(OUT)   :: problems

  ! Modified variables
  LOGICAL, INTENT(INOUT) :: do_aerosols, do_clouds

  ! Local variables
  INTEGER, PARAMETER     :: maxngas = 7, maxscatter=3, allngas = 9
  INTEGER, DIMENSION(maxgkmatc), PARAMETER :: &
       greekmat_idxs = (/1, 2, 5, 6, 11, 12, 15, 16/), phasmoms_idxs = (/1, 5, 5, 2, 3, 6, 6, 4/)

  INTEGER  :: ui, i, j, k, q, nscatter, idx, cldidx, aeridx, nactgksec, nactgkmatc, ig
  INTEGER, DIMENSION(allngas) :: absin
  REAL (KIND=dp)     :: scaco_r, absco_r, omega, &  ! raycof, depol, 
     extco_r, extco, scaco, pvar, extco_a, scaco_a, extco_c, scaco_c, j0, j1
  REAL (KIND=dp), DIMENSION(maxscatter)              :: scaco_input
  REAL (KIND=dp), DIMENSION(nlayers)                 :: extconf
  REAL (KIND=dp), DIMENSION(ngas, nlayers)           :: absod
  REAL (KIND=dp), DIMENSION(0:maxmoments_input, 1:maxgksec, maxscatter),     SAVE :: phasmoms_input
  REAL (KIND=dp), DIMENSION(0:maxmoments_input, 1:maxgksec),                 SAVE :: phasmoms_total_input
  REAL (KIND=dp), DIMENSION(max_atmoswfs, 0:maxmoments_input, 1:maxgksec), SAVE :: l_phasmoms_total_input
  LOGICAL, SAVE                                                             :: first = .TRUE.
  
  ! ========================== Check for Input ==================================
  problems = .FALSE.
  !IF (MAXVAL(gasin) > allngas) THEN
  !   WRITE(www_lun, *) 'Not weighting functions are implemented for all gases!!!'
  !   problems = .TRUE.; RETURN
  !ENDIF

  IF (first) THEN
     IF (.NOT. useasy .AND. nmoms > maxmoments_input) THEN
        WRITE(www_lun, *) 'Need to increase maxmoments_input for aerosols/clouds!!!'
        problems = .TRUE.; RETURN
     ENDIF
     
     ! This only needs to be initialized once
     phasmoms_input        = ZERO
     phasmoms_total_input  = ZERO
     l_phasmoms_total_input= ZERO
     VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input   = ZERO
     VLIDORT_LinFixIn%Optical%TS_l_omega_total_input   = ZERO
     VLIDORT_FixIn%Optical%TS_greekmat_total_input  = ZERO
     VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input= ZERO 
     
     first =.FALSE.
  ENDIF
  
  IF (NSTOKES == 1) THEN
     nactgksec = 1;  nactgkmatc = 1
  ELSE IF (NSTOKES == 3) THEN 
     !nactgksec = ngksec; nactgkmatc = ngkmatc
     nactgksec = 5; nactgkmatc = ngkmatc
  ELSE IF (NSTOKES == 4) THEN 
     nactgksec = 8; nactgkmatc = ngkmatc
  eLSE IF (NSTOKES == 2) THEN 
      WRITE(*,*) 'nstokes=2' ; STOP 1
  ENDIF
   

  ! 1: O3 2: NO2  3:O2  4: O4 5: BrO 6: H2O 7 SO2 8: HCHO  9: OCLO
  absin(:) = 0
  DO i =1, ngas
    absin(gasin(i)) = i
  ENDDO

  ! Disable clouds and aerosols if for Rayleigh scattering atmosph:q!ere
  !IF (do_rayleigh_only) THEN
  !   do_clouds = .FALSE.; do_aerosols = .FALSE.
  !ENDIF

  ! Enable delta-M scaling for clouds or aerosols
  !IF (do_clouds .OR. do_aerosols) do_deltam_scaling = .TRUE.

  ! Start layer loop

  ! Get rayleigh scattering phase function moments (Same for each layer)
  ! unassigned elements have already initialized to zero
  phasmoms_input(0, 1, 1) = ONE
  phasmoms_input(2, 1, 1) = (ONE - depol) / (TWO + depol)  
  IF (nactgksec /= 1) THEN
     phasmoms_input(2, 2, 1) = 6.0D0 * phasmoms_input(2, 1, 1)
     phasmoms_input(2, 5, 1) = -SQRT(6.0D0) * phasmoms_input(2, 1, 1)
     phasmoms_input(1, 4, 1) = 3.0D0 * (ONE - 2.0D0 * depol) / (TWO + depol)
  ENDIF
!  DO i = 1, 3 
!     print * , phasmoms_input(i-1, 1:6,1)
!  ENDDO 
  DO i = 1, nlayers   
     ! Rayleigh scattering
     scaco_r = raycof * airgrid(i)
     delray(i) = scaco_r

     IF (any(abscrs(1:ngas, i) < 0)) THEN
        DO ig = 1, ngas
           IF (abscrs(ig, i) < 0) abscrs(ig,i) = 0.0
        ENDDO
     ENDIF

     ! Gas absorption
     absod(1:ngas, i) = abscrs(1:ngas, i) * gascol(1:ngas, i)
     absco_r = SUM(absod(1:ngas, i))         
!    WRITE(*,'(f8.3, i4, 3e15.7)') lamda,i, absod(2, i), abscrs(2, i), gascol(2, i)

     IF (absco_r <= 0.0) THEN
        problems  = .TRUE.
        print *, 'Negative total absorption: ', lamda, i, ngas
        print *, abscrs(1:ngas, i)
        print *, gascol(1:ngas, i)
        RETURN
     ENDIF

     extco_r = absco_r + scaco_r
     scaco_input(1) = scaco_r
     
     ! Aerosols and clouds
     extco = extco_r
     nscatter = 1; cldidx = 0; aeridx = 0
     extco_a = ZERO; scaco_a = ZERO; extco_c = ZERO; scaco_c = ZERO; 

     !IF (.NOT. do_rayleigh_only) THEN        
        IF (do_clouds .AND. cldmsk(i)) THEN
           nscatter = nscatter + 1
           extco = extco + cldext(i)
           extco_c = cldext(i);   scaco_c = cldsca(i)
           scaco_input(nscatter) = cldsca(i)
           cldidx = nscatter
           
           ! get phase moments for clouds
           IF (.NOT. useasy) THEN
              phasmoms_input(0:nmoms, 1:nactgksec, nscatter) = cldmoms(0:nmoms, 1:nactgksec, i)                         
           ELSE ! use H_G function
              phasmoms_input(0, 2, nscatter) = ONE
              j0 = ONE             
              DO j = 1, nmoms
                 j1 = REAL(2*j+1, KIND=dp)
                 phasmoms_input(j, 2, nscatter) = (j1/j0) * cldasy(i) * phasmoms_input(j-1, 2, nscatter)
                 j0 = j1
              ENDDO
           ENDIF
           
        ENDIF  ! end clouds
        
        IF (do_aerosols .AND. aermsk(i)) THEN
           nscatter = nscatter + 1
           extco = extco + aerext(i)
           scaco_input(nscatter) = aersca(i)
           extco_a = aerext(i);  scaco_a = aersca(i)
           aeridx = nscatter

           ! get phase moments for aerosols
           IF (.NOT. useasy) THEN  
              phasmoms_input(0:nmoms, 1:nactgksec, nscatter) = aermoms(0:nmoms, 1:nactgksec, i)  
           ELSE ! use H_G function
              phasmoms_input(0, 2, nscatter) = ONE
              j0 = ONE             
              DO j = 1, nmoms
                 j1 = REAL(2*j+1, KIND=dp)
                 phasmoms_input(j, 2, nscatter) = (j1/j0) * aerasy(i) * phasmoms_input(j-1, 2, nscatter)

                 j0 = j1
              ENDDO
           ENDIF
        ENDIF  ! end aerosols
     !ENDIF     ! end non-rayleigh
     
     ! setup LIDORT input for tau and omega
     scaco = SUM(scaco_input(1:nscatter))
     omega = scaco / extco
     
     IF (omega < OMEGA_SMALLNUM) omega = OMEGA_SMALLNUM 
     IF (omega > 1.0 - OMEGA_SMALLNUM) omega = 1.0 - OMEGA_SMALLNUM
     VLIDORT_ModIn%MOptical%TS_omega_total_input(i) = omega
     VLIDORT_FixIn%Optical%TS_deltau_vert_input(i) = extco
  
     ! sum up phase moments as required in LIDORT
     DO j = 0, nmoms
        DO k = 1, nactgksec
           phasmoms_total_input(j, k) = SUM(phasmoms_input(j, k, 1:nscatter) &
                * scaco_input(1:nscatter)) / scaco
        ENDDO
     ENDDO
     !phasmoms_total_input(nmoms+1:maxmoments, 1:maxgksec) = 0.0  
     
     !extconf(i) = extco / (zsgrid(i-1) - zsgrid(i))   ! extinction coefficients
     ! Set up greek scattering matrix for each moment 
     !greekmat_total_input(0:nmoms, i, 1)  = phasmoms_total_input(0:nmoms, 1)
     !greekmat_total_input(0:nmoms, i, 2)  = phasmoms_total_input(0:nmoms, 5)
     !greekmat_total_input(0:nmoms, i, 5)  = phasmoms_total_input(0:nmoms, 5)
     !greekmat_total_input(0:nmoms, i, 6)  = phasmoms_total_input(0:nmoms, 2)
     !greekmat_total_input(0:nmoms, i, 11) = phasmoms_total_input(0:nmoms, 3)
     !greekmat_total_input(0:nmoms, i, 12) = phasmoms_total_input(0:nmoms, 6)
     !greekmat_total_input(0:nmoms, i, 15) = -phasmoms_total_input(0:nmoms, 6)
     !greekmat_total_input(0:nmoms, i, 16) = phasmoms_total_input(0:nmoms, 4)
     !greekmat_total_input(nmoms+1:maxmoments, i, 1:MAXSTOKES_SQ) = 0.0
     VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
          phasmoms_total_input(0:nmoms, phasmoms_idxs(1:nactgkmatc))
     IF ( nactgkmatc > 1 ) VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmoms, i, 15) &
          = -VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmoms, i, 15)

     ! This should always be 1, but may be slightly different due to numerical truncation
     VLIDORT_FixIn%Optical%TS_greekmat_total_input(0, i, 1) = 1.0  

     !IF  (i == 26) THEN
     !WRITE (91, '(I5, 7D24.12, I5)') i, extco, scaco, scaco_a, absco_r, scaco_r, &
     !     omega_total_input(i), greekmat_total_input(0, i, 1), nscatter
     !DO k = 1, 16
     !   WRITE (91, '(1000D24.12)') (greekmat_total_input(j, i, k), j=0, nmoms)
     !ENDDO
     !ENDIF

     !IF (do_simulation_only .OR. .NOT. varyprof(i)) THEN   ! no linearition
     !   ! zero out quantity for safety
     !   layer_vary_flag(i) = .FALSE.
     !   layer_vary_number(i) = 0	
     !   !VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(:, i) =   ZERO
     !   !VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(:, i) = ZERO
     !   !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(:, : , :, i) = ZERO        
     !
     !ELSE
     !   layer_vary_flag(i) = .TRUE.
     !   layer_vary_number(i) = VLIDORT_LinFixIn%Cont%TS_n_totalprofile_wfs
     !The above part has been taken care of in routine: lidort_prof_env.f90
        
     DO q = 1, VLIDORT_LinFixIn%Cont%TS_n_totalprofile_wfs        
        !  w.r.t ozone volume mixing ratio: 1
        !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF ( profilewf_names(q) == 'ozone volume mixing ratio------' ) THEN
           idx = absin(1)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:maxmoments , i, 1:16) = ZERO
           
           !  w.r.t NO2 volume mixing ratio: 2
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'NO2 volume mixing ratio------' ) THEN
           idx = absin(2)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O2 volume mixing ratio: 8
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O2 volume mixing ratio------' ) THEN
           idx = absin(8)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O4 volume mixing ratio: 3
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O4 volume mixing ratio------' ) THEN
           idx = absin(3)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.;  RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t BrO volume mixing ratio: 4
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'BrO volume mixing ratio------' ) THEN
           idx = absin(4)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t H2O volume mixing ratio: 9
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'H2O volume mixing ratio------' ) THEN
           idx = absin(9)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t SO2 volume mixing ratio: 5
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'SO2 volume mixing ratio------' ) THEN
           idx = absin(5)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t HCHO volume mixing ratio: 6
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'HCHO volume mixing ratio------' ) THEN
           idx = absin(6)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t OCLO volume mixing ratio: 7
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'OCLO volume mixing ratio------' ) THEN
           idx = absin(7)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - absod(idx, i) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i )= + absod(idx, i) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average temperature of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
           !  Assume no effects on air density              
        ELSE IF ( profilewf_names(q) == 'average temperature of layer---' ) THEN
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = - SUM(absod(:, i) * eta(:, i)) / extco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q, i) = + SUM(absod(:, i) * eta(:, i)) / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average pressure of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
        ELSE IF ( profilewf_names(q) == 'average pressure of layer------' ) THEN
           
           pvar = extco_r/ extco
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = ((ONE - pvar) * scaco_input(1) - &
                pvar * (scaco - scaco_input(1))) / scaco
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i) = extco_r / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t rayleigh soptical thickness
           ! xliu: August 12, 2008
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        ELSE IF ( profilewf_names(q) == 'rayleigh optical thickness-----' ) THEN
           pvar = scaco_r / extco  
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i) = (1.0 - omega) & 
                                      * scaco_r / extco /omega
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i) = pvar
           VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
           DO j = 0, nmoms
             DO k = 1, nactgksec
                IF (phasmoms_total_input(j, k) /= 0.0) THEN
                   l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, 1) - phasmoms_total_input(j, k) ) &
                        / phasmoms_total_input(j, k) * scaco_a / scaco
                ELSE
                   l_phasmoms_total_input(q, j, k) = 0.0
                ENDIF
             ENDDO
          ENDDO
          VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
               l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
          IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15) &
               = - VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15)          
        ELSE IF ( profilewf_names(q) == 'rayleigh scattering coefficient' ) THEN
           !  xliu: April 13, 2007 
           !  Still need to consider the variation in phase function  
           pvar = scaco_r / extco
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q, i) = ((ONE - pvar) * scaco_input(1) - &
                pvar * (scaco - scaco_input(1)) ) / scaco
           VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i) = pvar
           VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i)     = scaco_r / extco
           !VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t aerosol extinction coefficient / aerosol optical thickness
           !  aerosol scattering albedo does not change   
           !  xliu: April 13, 2007 (consider the variation in phase function)
        ELSE IF ( profilewf_names(q) == 'aerosol extinction coefficient-' ) THEN
           IF (aeridx > 0) THEN
              VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i) = + extco_a / extco
              VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i) = &
                  (scaco_a  / extco_a - omega )/omega  * extco_a / extco
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
              DO j = 0, nmoms
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aeridx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15) &
                   = - VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15)
              !print *, maxval(VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input), minval(VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input)
              !print *, i, VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i), VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i)
              !WRITE(www_lun, '(6D14.6)') (l_phasmoms_total_input(1, j, 1:nactgksec), j = 0, nmoms)
              !STOP 1
           ENDIF
           !  w.r.t  aerosol scattering coefficient / single scattering albedo
           !  aerosol optical thickness will not change
        ELSE IF ( profilewf_names(q) == 'aerosol scattering coefficient-' ) THEN
           IF (aeridx > 0) THEN
              VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i) = ZERO
              VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i) = scaco_a  / omega / extco
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
              DO j = 0, nmoms
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aeridx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15) &
                   = - VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15)
           ENDIF
           !   w.r.t cloud extinction coefficient / optical thickness
        ELSE IF ( profilewf_names(q) == 'cloud extinction coefficient---' ) THEN
           IF (cldidx > 0) THEN
              VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i) = + extco_c / extco
              VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i) = (scaco_c  / extco_c - omega ) &
                   / omega  * extco_c / extco
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
              DO j = 0, nmoms
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldidx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 1:nmoms, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15) &
                   = - VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15)
           ENDIF
           
           !  w.r.t clouds scattering coefficient
        ELSE IF ( profilewf_names(q) == 'cloud scattering coefficient---' ) THEN
           IF (cldidx > 0) THEN
              VLIDORT_LinFixIn%Optical%TS_l_deltau_vert_input(q,i) = ZERO
              VLIDORT_LinFixIn%Optical%TS_l_omega_total_input(q,i) = scaco_c  / omega / extco
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
              DO j = 0, nmoms
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldidx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15) &
                   = - VLIDORT_LinFixIn%Optical%TS_l_greekmat_total_input(q, 0:nmoms, i, 15)
           ENDIF
        ENDIF        ! end selection of weighting function 
     ENDDO           ! VLIDORT_LinFixIn%Cont%TS_n_totalprofile_wfs loop
  !ENDIF             ! end of do_linearization    
  ENDDO              ! layer loop
  
  deltau(1:nlayers) =  VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)
  delsca(1:nlayers) =  deltau(1:nlayers)*VLIDORT_ModIn%Moptical%TS_omega_total_input(1:nlayers)
    
    
!  ! Prepare for surface albedo  (unnecessary)
!  IF (do_lambertian_surface) THEN
!     DO i = 1, nstreams
!        bireflec_0 (0, i, 1) = ONE
!        emissivity (i) = ONE - lambertian_albedo
!        DO j = 1, nstreams
!           bireflec (0, i, j) = ONE
!        ENDDO
!     ENDDO
!     DO ui = 1, n_user_streams
!        user_bireflec_0 (0, ui, 1) = ONE
!        user_emissivity (ui) = ONE - lambertian_albedo
!        DO j = 1, nstreams
!           user_bireflec (0, ui, j) = ONE
!        ENDDO
!     ENDDO
!  ENDIF
!   
!  ! Prepare for spherical attenuation
!  IF ( .NOT. do_chapman_function ) THEN
!     CALL PREPARE_SPHERICAL (nlayers, do_plane_parallel, beam_szas(1), earth_radius, &
!          extconf, zsgrid, deltau_slant_input(1:nlayers, 1:nlayers, 1), &
!          sza_local_input(1:nlayers, 1))
!  ENDIF
  
END SUBROUTINE LIDORT_PROF_PREP

! Compute slant optical thickness using Chapman function in LIDORT
SUBROUTINE GET_SLANT_TAU(nz, zs, tauin, sza, tauout)
  USE OMSAO_parameters_module, ONLY  : rearth
  USE OMSAO_precision_module
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                          :: nz
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, tauin
  REAL (KIND=dp), INTENT(IN)                   :: sza
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(out) :: tauout
  ! =======================
  ! Local variables
  ! =======================
  INTEGER                   :: i, j

  VLIDORT_ModIn%MSunRays%TS_N_SZANGLES = 1
  VLIDORT_ModIn%MSunRays%TS_SZANGLES = sza
  IF (sza >= 90.0 .OR. sza < 0) THEN
     STOP 'GET_SLANT_TAU: SZA is >= 90 or < 0!!!'
  ENDIF

  nlayers = nz
  !! FIX ME taugrid_input(0:nz) = tauin

  VLIDORT_FixIn%Chapman%TS_height_grid(0:nz) = zs
  VLIDORT_ModIn%MChapman%TS_earth_radius = rearth
  IF (nz > maxlayers) THEN
     STOP 'LIDORT_PROF_ENV: # of layers exceeded allowed !!!'
  ENDIF

  !CALL VLIDORT_CHAPMAN(fail, message, trace)
  PRINT * , 'not implemented'; STOP 1
  tauout = 0.0
  DO i = 1, nz
     DO j = 1, i
        !! FIX ME tauout(i) = tauout(i) + VLIDORT_Work_Miscellanous%deltau_slant(i, j, 1)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE GET_SLANT_TAU

END MODULE m_lidort_env_vv2p7

