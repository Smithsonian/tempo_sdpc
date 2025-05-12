MODULE m_pseudo_model

  ! *********************************************************************************
  ! Author: Xingu Liu
  ! Date:   July 24, 2003
  ! Purpose: calculated simulated reflectance and weighting function/first derivative
  ! for all the fitting variables. Basically, call the lidort to calculate radiance
  ! and ozone, albedo weighting functions, where species other than ozone are not
  ! taken into account. Then calibrate the measured radiance/solar spectra to obtain
  ! the measured reflectance and use the finite difference to obtain the first
  ! derivative for all the other variables.
  ! xliu, 08/10/2010
  ! Current VLIDORT calculation is based on single surface albedo (per channel).
  ! The wavelength dependence of surface albedo on radiance is corrected through
  ! weighting function. Howver, there are not accounted for in the calculation of
  ! weighting functions.
  ! Jbak (2017 July to 2020~)
  !   adding slit function fitting variables (do_pslwf)
  !   adding PCA RTM simulation scheme
  !   adding visible simluation with hres cross section
  !   adding negative value correction for all trace gases
  !   improving w-cloud fraction fitting from Cai
  !   improving cloud fitting from Cai
  !   adding BRDF and EOF-based spectrum from chris's work
  !   ** Need to do more about the albedo and polarization reltaed fitting
  !   ** Need to do more about negative albedo treatment with vary_sfcalb
  ! *********************************************************************************

  USE m_lidort_env_vv2p7,     ONLY: lidort_prof_env
  USE m_lidort_env_vv2p7_pca, ONLY: lidort_prof_env_pca
  USE m_set_brdf, ONLY: surface
  INTEGER, PARAMETER :: nalbwf = 1 ! should increase if do_brdf is turn on
  PUBLIC  pseudo_model 
  PRIVATE HRES_RADCALC_ENV

CONTAINS
SUBROUTINE pseudo_model (num_iter, refl_only, ns, nf, fitvar, fitvarap, dyda, gspec,    &
     fitres, fitspec, fitqres, fitq, chisq, relrms, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : maxlay, maxwin
  USE OMSAO_indices_module,   ONLY : ring_idx,  &
       maxalb, maxoth, bro_idx, so2_idx,   &
       no2_t1_idx, hcho_idx, shift_offset, us1_idx, us2_idx, maxwfc, so2v_idx, &
       o2o2_idx, h2o_idx, o2_idx, &
       com_idx, com1_idx, fsl_idx, rsl_idx, com2_idx, com3_idx
  USE OMSAO_variables_module, ONLY : scnwrt, numwin, nviswin, band_selectors, &
       refidx, refspec_norm, npix_fitted, nradpix, &
       fitwavs, fitweights, actspec_rad, &
       fitvar_rad, fitvar_rad_apriori, fitvar_rad_str, &
       mask_fitvar_rad, rmask_fitvar_rad,  &
       database_indices, database, database_shiwf, slwf, &
       npsl, psl_fpos, database_pslwf, database_cmwf, &
       stokfrac, stokwaves
  USE ozprof_data_module,     ONLY : nlay, use_lograd, do_debug_o3p,         &
       ozf_fidx => ozfit_start_index, ozf_lidx => ozfit_end_index,           &
       ozp_fidx=>ozprof_start_index, ozp_lidx => ozprof_end_index,           &
       stlay => start_layer, endlay => end_layer, albfidx, nalb, nfalb,      &
       albidx, albfpix, alblpix, t_fidx, t_lidx, tf_fidx,    &
       tf_lidx, nt_fit, do_subfit, saa_flag, fgasidxs, tracegas, ngas,       &
       gasidxs, osind, osfind, slind, slfind, shind, shfind, rnind, rnfind,  &
       dcind, dcfind, isind, isfind, irind, irfind, oswins, slwins, shwins,  &
       rnwins, dcwins, iswins, irwins, nos, nsl, nsh, nrn, ndc, nis, nir,    &
       p1ind, p1find, p1wins, np1, p2ind, p2find, p2wins, np2, &
       p3ind, p3find, p3wins, np3, cmfind, cmind, cmwins, ncm, &
       nothgrp, the_cfrac, polcorr, radcalwrt, do_simu, ecfrfind, ecodfind,  &
       ectpfind, taodfind, twaefind, saodfind, ecfrind, &
       sprsfind, wfcfidx, nwfc, nfwfc,   &
       wfcidx, wfcfpix, wfclpix, so2zfind, fit_atanring, &
       use_effcrs, do_simu_rmring, nsaa_spike, the_ai, &
       do_alb_longwav, alb_swav, alb_ewav,  which_cld, use_prefitalb, pf2ba0, pf2ba1, pf2fc0, pf2fc1, &
       do_sy_diagonal, merr_covar, &
       is_albspcvar, albspcs, use_albspc, use_albeofs, nactalbspc, sfcalbs, &
       do_rtm_pca, do_brdf, vary_sfcalb, the_snowice, rtm_outputs, set_rtmvar, &
       allrms, allradrms
  USE cloud_data_module, ONLY : new_cfrac, new_alb, n_newalb, &
       n_newwfc, use_retalb, avgwav_cld
  USE OMSAO_gome_data_module, ONLY: n_gome_q, gome_q, stkidx, &
      use_origin_q ,n_stokfrac
  USE m_gsvd_o3prof_utilities, ONLY: sq_matrix_invert
  USE m_fitting_util, ONLY: uv1_spike_detect, spike_detect_correct
  USE m_spectra_reflectance, ONLY: spectra_reflectance
  USE m_lidort_util, ONLY: debug_rtm, hres_stkwf_inter_convol2
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns, nf, num_iter
  INTEGER, INTENT(OUT)                            :: errstat
  LOGICAL, INTENT(INOUT)                          :: refl_only
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (nf)   :: fitvar, fitvarap
  REAL (KIND=dp), INTENT(OUT), DIMENSION (ns)     :: gspec, fitres, fitspec
  REAL (KIND=dp), INTENT(OUT), DIMENSION (ns)     :: fitq, fitqres
  REAL (KIND=dp), INTENT(OUT), DIMENSION (ns, nf) :: dyda
  REAL (KIND=dp), INTENT(OUT)                     :: chisq, relrms

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER :: mstks = 4 ! jbak it is changed from 4 to 1 to reduce memory because Q/U is not used here
  TYPE (rtm_outputs) :: rio
  INTEGER :: n0alb, idx0alb(maxalb), n0wfc, i, j, k, iw, ridx, sidx, fidx, lidx, &
       albord, min_ssa_iter, swin, ewin, ig, nord, nostk, wfcord, &
       albsidx, albeidx, m, slit_idx, idx330 
  INTEGER, DIMENSION(maxalb)               :: albpmax, albpmin
  INTEGER, DIMENSION(maxwfc)               :: wfcpmax, wfcpmin
  INTEGER, DIMENSION(maxwin, maxoth)       :: tmpind, tmpfind
  INTEGER, DIMENSION(maxoth, 2)            :: tmpwins
  REAL (KIND=dp), DIMENSION(maxalb)        :: albarr 
  REAL (KIND=dp), DIMENSION(ns)            :: albadj
  REAL (KIND=dp), DIMENSION(maxwfc)        :: wfcarr 
  REAL (KIND=dp), DIMENSION (ns)           :: delpos, waves, meas1, meas2, sim1, sim2, &
       simrad, simrad1, fitspec1, temporwf, simq, Qsimrad, direc !, tempo, simqcor
  REAL (KIND=dp), DIMENSION(ns,nlay,mstks) :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(ns,nlay)       :: q_ozwf ! Q/I wf w.r.t. ozone,
  REAL (KIND=dp), DIMENSION(ns, 0:maxalb)  :: albothwf, wfcothwf
  REAL (KIND=dp), DIMENSION(ns, nalbwf, mstks) :: albwf
  REAL (KIND=dp), DIMENSION(ns, mstks)     :: o3shiwf, cfracwf, fsimrad, &
       ctpwf, codwf, saodwf, taodwf, twaewf, sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(numwin, maxoth):: o3shi
  REAL (KIND=dp), DIMENSION(nlay)          :: tprof, ozprof, ozadj, ozaprof
  REAL (KIND=dp), DIMENSION(nf)            :: fitvar_saved
  REAL (KIND=dp)                           :: rms, radrms, wavavg, newoz
  REAL (KIND=dp)                           :: so2adj, so2vadj, newso2, newbro, &
       newhcho, newno2, newo4, newo2, newh2o, broadj, hchoadj, no2adj, o4adj, o2adj, h2oadj, newalb 
  REAL (KIND=dp), DIMENSION (numwin)       :: allchisq
  LOGICAL :: do_ozwf, do_albwf, do_o3shi, do_tmpwf, do_shiwf, do_taodwf, do_twaewf, &
       do_saodwf, do_cfracwf, do_ctpwf, do_codwf, negval, do_sprswf, do_so2zwf, do_pslwf, &
       so2negval, so2vnegval, hchonegval, bronegval, no2negval, o4negval, o2negval, h2onegval, albnegval
  LOGICAL, DIMENSION (nlay)     :: ozvary
  REAL (KIND=dp), DIMENSION(ns) :: walb0s, wfc0s 
  CHARACTER(LEN=1)  :: ordchar
  LOGICAL           :: no2alb0
  ! measurement error covariance Random + Systematic
  REAL(KIND=dp), DIMENSION(ns, ns)  :: Sy, Sy_inv
  REAL(KIND=dp), DIMENSION(ns, 1)   :: y1
  REAL(KIND=dp), DIMENSION(1, 1)    :: chi
  REAL(KIND=dp), ALLOCATABLE        :: y1tmp(:, :), Sy_invtmp(:, :)

  LOGICAL, SAVE :: first
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  INTEGER :: funit
  LOGICAL :: do_debug_out = .false.
  CHARACTER (LEN=12), PARAMETER :: modulename = 'pseudo_model'

  !============================================================================
  ! debuging part, it might be done in read_ozprof_input. Re-done for safe
  !============================================================================
  IF (first) THEN
    IF (use_albspc) vary_sfcalb = .TRUE.
    IF (do_rtm_pca .and. nfwfc > 0) THEN
       WRITE(*,'(A)') modulename//'nfwfc is not implemnted in PCA simulation' ; stop 1
    ENDIF
    IF (do_rtm_pca .and. do_brdf) THEN
       WRITE(*,'(A)') modulename//'BRDF is not implemnted in PCA simulation' ; stop 1
    ENDIF
    IF (do_brdf) THEN
       IF (nalbwf ==1) THEN
          WRITE(*,'(A)') modulename//'nalbwf need to increase in BRDF mode' ; stop 1
       ENDIF
       WRITE(*,'(A)') modulename//'BRDF need to be improved'
    ENDIF
    IF (use_prefitalb) THEN
       WRITE(*,'(A)') modulename//'prefitalb is not really implemented'; stop 1
    ENDIF
    IF (do_alb_longwav) THEN
       IF (nviswin > 1) THEN
           WRITE(*,'(A)') modulename//'do_alb_longwav should be only for UV channel'
       ENDIF
       WRITE(*,'(A)') modulename//'do_alb_longwav need to be improved'
    ENDIF
    IF (n_gome_q > 0 .and. mstks == 1) THEN
       WRITE(*,'(A)') modulename//'mstks need to increase 3 for  GOME polarization from Q/U' ; stop 1
    ENDIF
    first = .true.
  ENDIF

  IF (do_debug_o3p) WRITE(www_lun, '(A,i3)') modulename//':  (start)', num_iter

  ! Initialize
  errstat = pge_errstat_ok
  ! ================ Determine flags for linearization ======================
  do_albwf  = .TRUE.; do_ozwf = .TRUE.; do_o3shi = .TRUE.; do_tmpwf = .TRUE.
  do_cfracwf = .TRUE.; do_pslwf = .TRUE.
  IF (nfalb <= 0) do_albwf = .FALSE.
  IF (nfwfc <= 0) do_cfracwf = .FALSE.
  IF (refl_only .AND. .NOT. use_effcrs) THEN
     do_ozwf = .FALSE.; do_o3shi = .FALSE.; do_tmpwf = .FALSE.
  END IF
  IF (nos <= 0)    do_o3shi = .FALSE.
  IF (npsl <= 0)    do_pslwf  = .FALSE.
  IF (nt_fit <= 0) do_tmpwf = .FALSE. 
  ozvary = .FALSE.; ozvary(stlay:endlay) = .TRUE.
  do_ctpwf  = .FALSE.; do_codwf  = .FALSE.
  do_taodwf = .FALSE.; do_twaewf = .FALSE.; do_saodwf = .FALSE.
  do_sprswf = .FALSE.; do_so2zwf = .FALSE. 
  IF (.NOT. refl_only) THEN
     IF (ecfrfind > 0) do_cfracwf = .TRUE.
     IF (ecodfind > 0) do_codwf   = .TRUE.
     IF (ectpfind > 0) do_ctpwf   = .TRUE.
     IF (taodfind > 0) do_taodwf  = .TRUE.
     IF (saodfind > 0) do_saodwf  = .TRUE.
     IF (twaefind > 0) do_twaewf  = .TRUE.
     IF (sprsfind > 0) do_sprswf  = .TRUE.
     IF (so2zfind > 0) do_so2zwf  = .TRUE.
  ENDIF

  ! Update cloud fraction or disable fitting cloud fraction when fc<=1.0E-2 or fc>=0.99
  IF (ecfrfind > 1) THEN !RECHECK
     the_cfrac = fitvar_rad(ecfrind)
     IF (the_cfrac <= 0.2) THEN
        do_cfracwf = .FALSE.
     ENDIF
     IF (the_cfrac <= 0.01) THEN
        the_cfrac = 0.0D0;  fitvar(ecfrfind) = 0.0D0; fitvar_rad(ecfrind) = 0.0D0
        fitvarap(ecfrfind) = 0.0D0;    fitvar_rad_apriori(ecfrind) = 0.0D0
        do_cfracwf = .FALSE.      
     ELSE IF (the_cfrac >= 0.99) THEN
        the_cfrac = 1.0D0;  fitvar(ecfrfind) = 1.0D0; fitvar_rad(ecfrind) = 1.0D0
        fitvarap(ecfrfind) = 1.0D0;    fitvar_rad_apriori(ecfrind) = 1.0D0
        do_cfracwf = .FALSE. 
        ! When using Lambertian cloud model, surface albedo is already fitted, so nothing needs to be done
        ! But when using scattering cloud model, need to fit an optical thickness and a surface albedo
        ! To be implemented !!!
     ENDIF
  ENDIF

  ! Disable fitting surface albedo when 0.99>=fc>=0.01
  ! JBAK For safe, as we have set sa(i, i) = 0.0 in specfit_ozprof 
  !      when fitting both w-dependent cloud fraction and surface albedo,
  !      these two state vectors are compete eath others, therefore
  !      if one is dominant, other is turn off.
  IF (nfwfc > 0 ) THEN ! ecfraind should be zero
     DO i = wfcidx + nwfc - 1, wfcidx, -1
        IF (fitvar_rad_str(i)(4:4) == '0') EXIT
     ENDDO
     the_cfrac = fitvar_rad(i)   ! Use UV2/last channel cloud fraction
     DO i = albidx + nalb - 1, albidx, -1
        IF (fitvar_rad_str(i)(4:4) == '0') EXIT
     ENDDO
     IF (the_cfrac <= 1.0E-2) THEN
        the_cfrac = 0.0D0; do_cfracwf = .FALSE.
        DO i = 1, nwfc
           j = wfcidx - 1 + i
           fitvar_rad(j) = 0.0D0;  fitvar_rad_apriori(j) = 0.0D0
        ENDDO
        DO i = wfcfidx, wfcfidx + nfwfc - 1
           fitvar(i) = 0.0D0;  fitvarap(i) = 0.0D0
        ENDDO
     ELSE IF (the_cfrac >= 0.99) THEN
        the_cfrac = 1.0D0; do_cfracwf = .FALSE.
        DO i = 1, nwfc
           j = wfcidx - 1 + i
           IF (fitvar_rad_str(j)(4:4) == '0') THEN
              fitvar_rad(j) = 1.0D0;  fitvar_rad_apriori(j) = 1.0D0
           ELSE
              fitvar_rad(j) = 0.0D0;  fitvar_rad_apriori(j) = 0.0D0
           ENDIF
        ENDDO
        DO i = wfcfidx, wfcfidx + nfwfc - 1           
           IF (fitvar_rad_str(mask_fitvar_rad(i))(4:4) == '0') THEN
              fitvar(i) = 1.0D0;  fitvarap(i) = 1.0D0
           ELSE
              fitvar(i) = 0.0D0;  fitvarap(i) = 0.0D0
           ENDIF
        ENDDO
     ENDIF
     IF ( the_cfrac <= 0.2 .or. (the_snowice >= 1 .and. the_snowice < 104)) THEN 
          do_cfracwf = .false.
          DO i = 1, nwfc
           j = wfcidx - 1 + i
           IF (fitvar_rad_str(j)(4:4) /= '0') THEN
              fitvar_rad(j) = 0.0D0;  fitvar_rad_apriori(j) = 0.0D0
           ENDIF
          ENDDO
          DO i = wfcfidx, wfcfidx + nfwfc - 1
              IF (fitvar_rad_str(mask_fitvar_rad(i))(4:4) /= '0') THEN
                 fitvar(i) = 0.0D0;  fitvarap(i) = 0.0D0
              ENDIF
          ENDDO
     ELSE IF (the_cfrac >=0.8 .and. (the_snowice <1 .and. the_snowice >= 104)) THEN
           do_albwf = .false.
     ELSE IF ( the_cfrac > 0.2 .and. do_cfracwf) THEN  
            do_albwf = .false.
            DO i = 1, nalb
              j = albidx -1 + i 
              IF (fitvar_rad_str(j)(4:4) /= '0') THEN
                  fitvar_rad(i) = 0.0D0;  fitvar_rad_apriori(i) = 0.0D0
              ELSE
                 IF (fitvar_rad(i) < 0.0 ) THEN
                    fitvar_rad(i) = 0.0D0;  fitvar_rad_apriori(i) = 0.0D0
                 ENDIF
              ENDIF
            ENDDO

            DO i = albfidx , albfidx + nfalb -1 
              IF (fitvar_rad_str(mask_fitvar_rad(i))(4:4) /= '0') THEN
                  fitvar(i) = 0.0D0;  fitvarap(i) = 0.0D0
              ELSE
                 IF (fitvar(i) < 0.0) THEN 
                    fitvar(i) = 0.0D0;  fitvarap(i) = 0.0D0
                 ENDIF
              ENDIF
            ENDDO
     ENDIF
  ENDIF

  ! ======= Set up ozone, temperature, trace gases, albedo, lamda for LIDORT ============
  ozprof(1:nlay)  = fitvar_rad (ozp_fidx:ozp_lidx)
  ozaprof(1:nlay) = fitvar_rad_apriori(ozp_fidx:ozp_lidx)
  tprof(1:nlay)   = fitvar_rad(t_fidx:t_lidx)

  IF (do_subfit) THEN
     DO i = 1, maxoth
        o3shi(1:numwin, i) = fitvar_rad(osind(1:numwin, i))
     ENDDO
  ELSE
     o3shi(1, 1:maxoth)    = fitvar_rad(osind(1, 1:maxoth))
  ENDIF

  ! set up for how many stokes we need
  nostk = 1
  IF (polcorr >= 3 .AND. polcorr < 5) THEN
    nostk = 3
  ELSE
    IF (polcorr == 0 .OR. polcorr ==6) nostk = 3
  ENDIF
  IF (mstks == 1) nostk = 1

  !xliu (02/01/2007): adjust ozone profile for negative ozone values
  !                   radiances will be corrected using ozone weighting function
  !@ Set up flag for negative ozone value
  ozadj(1:nlay) = 0.0; negval = .FALSE.
  DO i = 1, nlay
     IF (ozprof(i) <= 0.0) THEN
        newoz  = MIN(0.5, ozaprof(i))
        negval = .TRUE. ; ozadj(i)  = newoz - ozprof(i); ozprof(i) = newoz
        do_ozwf = .TRUE.; ozvary(i) = .TRUE.
     ENDIF
  ENDDO

  !@ Set up flag for other tracegases
  so2negval = .FALSE.; so2vnegval = .FALSE. ; hchonegval= .FALSE.; bronegval  = .FALSE.
  no2negval = .FALSE.; o4negval   = .FALSE. ; o2negval = .FALSE. ; h2onegval = .FALSE.
  DO k = 1, ngas
   i = fgasidxs(k)
   IF (i > 0) THEN
      j = mask_fitvar_rad(i)
      tracegas(k, 4) = fitvar_rad(j) !/ refspec_norm(gasidxs(k)) ! trace gas column in molecules cm-2

      IF (gasidxs(k) == so2_idx .AND. tracegas(k, 4) < 0.d0)  THEN
         so2negval = .TRUE.
         newso2 = 1.0E15 * refspec_norm(gasidxs(k))
         so2adj = newso2 - tracegas(k, 4)
         tracegas(k, 4) = newso2
      ENDIF
      IF (gasidxs(k) == so2v_idx .AND. tracegas(k, 4) < 0.d0)  THEN
         so2vnegval = .TRUE.
         newso2 = 1.0E15 * refspec_norm(gasidxs(k))
         so2vadj = newso2 - tracegas(k, 4)
         tracegas(k, 4) = newso2
      ENDIF
      IF (gasidxs(k) == no2_t1_idx .AND. tracegas(k, 4) < 0.d0)  THEN
         no2negval = .TRUE.
         newno2 = 1.0E15 * refspec_norm(gasidxs(k))
         no2adj = newno2 - tracegas(k, 4)
         tracegas(k, 4) = newno2
      ENDIF
      IF (gasidxs(k) == hcho_idx .AND. tracegas(k, 4) < 0.d0)  THEN
         hchonegval = .TRUE.
         newhcho = 1.0E15 * refspec_norm(gasidxs(k))
         hchoadj = newhcho- tracegas(k, 4)
         tracegas(k, 4) = newhcho
      ENDIF
      IF (gasidxs(k) == bro_idx .AND. tracegas(k, 4) < 0.d0)  THEN
         bronegval = .TRUE.
         newbro = 1.0E12 * refspec_norm(gasidxs(k))
         broadj = newbro - tracegas(k, 4)
         tracegas(k, 4) = newbro
      ENDIF
      IF (gasidxs(k) == o2o2_idx .AND. tracegas(k,4) < 0.d0) THEN 
         o4negval = .TRUE.
         newo4 = 1.3E33 * refspec_norm(gasidxs(k))
         o4adj = newo4 - tracegas(k, 4)
         tracegas(k, 4) = newo4
      ENDIF

      IF (gasidxs(k) == o2_idx .AND. tracegas(k,4) < 0.d0) THEN
         o2negval = .TRUE.
         newo2 = 4.5E24 * refspec_norm(gasidxs(k))
         o2adj = newo2 - tracegas(k, 4)
         tracegas(k, 4) = newo2
      ENDIF

      IF (gasidxs(k) == h2o_idx .AND. tracegas(k,4) < 0.d0) THEN
         h2onegval = .TRUE.
         newh2o = 1.0E22 * refspec_norm(gasidxs(k))
         h2oadj = newh2o - tracegas(k, 4)
         tracegas(k, 4) = newh2o
      ENDIF
    ENDIF
  ENDDO

  !@ Set up W-dependent surface albedo
  waves = fitwavs(1:ns) 
  n0alb = 0
  idx0alb(1:nalb) = 0
  no2alb0 = .false.
  DO i = 1, nalb
     j = albidx - 1 + i
     READ(fitvar_rad_str(j)(4:5), '(I2)') albord
     fidx=albfpix(i); lidx=alblpix(i)
     wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
     IF (.NOT. is_albspcvar(i)) THEN 
       IF (albord == 0) THEN
          n0alb = n0alb + 1
          !IF (fitvar_rad(j) > 1.0) THEN
          !   fitvar_rad(j) = 1.0
          !   IF (rmask_fitvar_rad(j) > 0) fitvar(rmask_fitvar_rad(j)) = 1.0
          !ENDIF
          albarr(n0alb)  = fitvar_rad(j)
          albpmax(n0alb) = lidx; albpmin(n0alb) = fidx
          walb0s(fidx:lidx) = fitvar_rad(j)
          idx0alb(n0alb) = j
          ! Borrow 1st order albedo from longer wavelength, assuming that for 312-355 nm the wavelenght dependent 
          ! parameter are same for each sub window 312-330 and 330-350 nm, but deriving 0-order, independently.
          IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j)(1:4) == '2ba0' .AND. rmask_fitvar_rad(j) <= 0) THEN
            DO k = i, nalb
              m = albidx -1 + k
              IF (fitvar_rad_str(m)(1:4) == '3ba0'.AND. rmask_fitvar_rad(m) > 0) EXIT 
            ENDDO
            IF (k > nalb) THEN
              WRITE(*, *)  modulename, ': ERROR 2ba0 is not filled'
              errstat = pge_errstat_error; RETURN   
            ENDIF
            albarr(n0alb)     = fitvar_rad(m)
            walb0s(fidx:lidx) = albarr(n0alb)
            no2alb0 = .TRUE.
           !ELSE IF (do_alb_longwav.AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
           !   fidx=albfpix(i); lidx=alblpix(i)
           !   walb0s(fidx:lidx) = new_alb(1)
          ELSE IF (use_prefitalb .AND. fitvar_rad_str(j)(1:4) == '2ba0' .AND. rmask_fitvar_rad(j) <= 0) THEN
            albarr(n0alb) = pf2ba0
            walb0s(fidx:lidx) = pf2ba0
          ENDIF
       ELSE
          IF (vary_sfcalb) THEN
             walb0s(fidx:lidx) = walb0s(fidx:lidx) + fitvar_rad(j) * (waves(fidx:lidx) - wavavg)**albord
             IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j)(1:3) == '2ba' &
                .AND. rmask_fitvar_rad(j) <= 0)  THEN
              WRITE(ordchar, '(I1)') albord 

              DO k = i, nalb
                 m = albidx -1 + k
                 IF (fitvar_rad_str(m) == '3ba'//ordchar.AND. rmask_fitvar_rad(m) > 0) EXIT 
              ENDDO
              IF (k > nalb) THEN
                 WRITE(*, *)  modulename, ': ERROR 2ba1 is not filled'
                 errstat = pge_errstat_error; RETURN   
              ENDIF
              IF (no2alb0) THEN ! if UV2 the zero fit value is not retrieved
                 fidx = albfpix(k)
                 lidx = alblpix(k)
                 wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
              ENDIF
              fidx = albfpix(i)
              lidx = alblpix(i)
              walb0s(fidx:lidx) = walb0s(fidx:lidx) + fitvar_rad(m) * (waves(fidx:lidx) - wavavg)**albord               
             ELSE IF (do_alb_longwav.AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
              fidx = albfpix(i)
              lidx = alblpix(i)
              m = albord + 1
              IF (m > n_newalb) THEN
                 WRITE(*, *)  modulename, ': which_cld==5, alb variables not consist'
              ENDIF
              walb0s(fidx:lidx) = walb0s(fidx:lidx) + new_alb(m) * (waves(fidx:lidx) - avgwav_cld)**albord
           
              !walb0s(fidx:lidx) = walb0s(fidx:lidx) + new_alb(m) * (waves(fidx:lidx) - wavavg)**albord
             ELSE IF (use_prefitalb .AND. fitvar_rad_str(j)(1:3) == '2ba'  &
                .AND. rmask_fitvar_rad(j) <= 0)  THEN
              fidx = albfpix(i); lidx = alblpix(i)
              walb0s(fidx:lidx) = walb0s(fidx:lidx) + pf2ba1 * (waves(fidx:lidx) - wavavg)**albord
             ENDIF
          ENDIF
       ENDIF
    ELSE 
        IF (do_brdf .and. use_effcrs) THEN
        ELSE IF (use_albeofs) THEN
           IF (nactalbspc == 1 ) THEN ! Snow/ice/water
              IF (albord == 1) THEN
                 walb0s(fidx:lidx) = albspcs(fidx:lidx, 0) * fitvar_rad(j)
              ELSE
                 walb0s(fidx:lidx) = walb0s(fidx:lidx) + (waves(fidx:lidx) - wavavg)**(albord-1) &
                      * albspcs(fidx:lidx, 0) * fitvar_rad(j)
              ENDIF
           ELSE
              IF (albord == 1)  walb0s(fidx:lidx) = albspcs(fidx:lidx, 0)
              walb0s(fidx:lidx) =  walb0s(fidx:lidx) + albspcs(fidx:lidx,albord) * fitvar_rad(j)
           ENDIF
        ELSE
           IF (albord == 0) THEN
              walb0s(fidx:lidx) = albspcs(fidx:lidx, 0) * fitvar_rad(j)
           ELSE
              walb0s(fidx:lidx) = walb0s(fidx:lidx) + (waves(fidx:lidx) - wavavg)**albord &
                   * albspcs(fidx:lidx, 0) * fitvar_rad(j)
           ENDIF
        ENDIF
     ENDIF
    ! WRITE(*,'(4i4, 2f8.2, 4e15.5)') i, albord, fidx, lidx,waves(fidx), waves(lidx), walb0s(fidx),walb0s(lidx), albspcs(fidx, albord), albspcs(lidx, albord)
  ENDDO
  IF (num_iter == 0) sfcalbs(1:ns, 1) = walb0s(1:ns) 
  !IF (nactalbspc > 2 .and. the_cfrac > 0)   walb0s(fidx:lidx) =
  !walb0s(fidx:lidx)*0.1
  !print * , 'avgalb:', sum(walb0s(fidx:lidx))/(lidx-fidx+1),nactalbspc,the_cfrac

  n0wfc = 0; no2alb0 = .FALSE.
  IF (nwfc > 0) THEN  
   DO i = 1, nwfc
     j = wfcidx - 1 + i
     READ(fitvar_rad_str(j)(4:4), '(I1)') wfcord
     IF (wfcord == 0) THEN
        n0wfc = n0wfc + 1
        wfcarr(n0wfc) = fitvar_rad(j)
        wfcpmax(n0wfc) = wfclpix(i); wfcpmin(n0wfc) = wfcfpix(i)
        IF (vary_sfcalb) wfc0s(wfcpmin(n0wfc):wfcpmax(n0wfc)) = wfcarr(n0wfc)
        IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j) == '2fc0' &
             .AND. rmask_fitvar_rad(j) <= 0)  THEN
           DO k = i, nwfc
              m = wfcidx -1 + k
              IF (fitvar_rad_str(m) == '3fc0'.AND. rmask_fitvar_rad(m) > 0) EXIT 
           ENDDO
           IF (k > nwfc) THEN
              WRITE(*, *)  modulename, ': ERROR 2fc0 is not filled'
              errstat = pge_errstat_error; RETURN   
           ENDIF
           fidx=wfcfpix(i); lidx=wfclpix(i)
           wfcarr(n0wfc) = fitvar_rad(m)
           IF (vary_sfcalb)  wfc0s(fidx:lidx) = fitvar_rad(m)      
           no2alb0 = .TRUE.
           !ELSE IF (do_alb_longwav.AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
           !   fidx=wfcfpix(i); lidx=wfclpix(i)
           !   wfc0s(fidx:lidx) = new_cfrac(1)
        ELSE IF (use_prefitalb .AND. fitvar_rad_str(j) == '2fc0' &
             .AND. rmask_fitvar_rad(j) <= 0 .AND. vary_sfcalb) THEN
           wfcarr(n0wfc) = pf2fc0
           fidx=wfcfpix(i); lidx=wfclpix(i)
           IF (vary_sfcalb) wfc0s(wfcpmin(n0wfc):wfcpmax(n0wfc)) =  pf2fc0
        ENDIF
     ELSE
        IF (vary_sfcalb) THEN
           fidx = wfcfpix(i); lidx = wfclpix(i)   
           wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
           wfc0s(fidx:lidx) =  wfc0s(fidx:lidx) + fitvar_rad(j) * (waves(fidx:lidx) - wavavg)**wfcord
           IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j)(1:3) == '2fc' &
                .AND. rmask_fitvar_rad(j) <= 0)  THEN
              WRITE(ordchar, '(I1)') wfcord 

              DO k = i, nwfc
                 m = wfcidx -1 + k
                 IF (fitvar_rad_str(m) == '3fc'//ordchar.AND. rmask_fitvar_rad(m) >= 0) EXIT 
              ENDDO
              IF (k > nwfc) THEN
                 WRITE(*, *)  modulename, ': ERROR 2fc1 is not filled'
                 errstat = pge_errstat_error; RETURN   
              ENDIF
              IF (no2alb0) THEN
                 fidx = wfcfpix(k)
                 lidx = wfclpix(k)
                 wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
              ENDIF
              fidx = wfcfpix(i)
              lidx = wfclpix(i)
              wfc0s(fidx:lidx) = wfc0s(fidx:lidx) + fitvar_rad(m) * (waves(fidx:lidx) - wavavg)**wfcord
           ELSE IF (do_alb_longwav .AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
              fidx = wfcfpix(i)
              lidx = wfclpix(i)
              m = wfcord + 1
              IF (m > n_newwfc) THEN
                 WRITE(*, *)  modulename, ': which_cld==5, cfrac variables not consist'
              ENDIF
              wfc0s(fidx:lidx) = wfc0s(fidx:lidx) + new_cfrac(m) * (waves(fidx:lidx) - avgwav_cld)**wfcord
           ELSE IF (use_prefitalb .AND. fitvar_rad_str(j)(1:3) == '2fc'  &
                .AND. rmask_fitvar_rad(j) <= 0)  THEN
              fidx = wfcfpix(i)
              lidx = wfclpix(i)
              wfc0s(fidx:lidx) = wfc0s(fidx:lidx) + pf2fc1 * (waves(fidx:lidx) - wavavg)**wfcord
           ENDIF
        ENDIF
     ENDIF
   ENDDO
  ENDIF

  ! Negative albedo treatment
  albadj(1:ns) = 0.0
  albnegval = .false.
  IF (any(albarr(1:n0alb) < 0) .or. any(walb0s(1:ns) <0)) THEN
    albnegval = .TRUE.
    albadj = 0.0
    newalb = 0.01 !Junsung: this should be used for omler
    !newalb = 0.05 !Junsung: 0.01 has been changed to 0.05 (12/03/2025) for gler
    IF (.NOT. vary_sfcalb) THEN
      IF (any(albarr(1:n0alb) <0)) THEN 
        DO i = 1, n0alb 
          IF (albarr(i) < 0) THEN 
            albadj(i) = newalb-albarr(i)
            albarr(i) = newalb
            walb0s(albpmin(i):albpmax(i)) = newalb 
          ENDIF
        ENDDO
      ELSE
        albnegval = .false.
      ENDIF
    ELSE
      DO i = 1, nalb
         j = albidx - 1 + i
         READ(fitvar_rad_str(j)(4:5), '(I2)') albord
         fidx=albfpix(i); lidx=alblpix(i)
         IF (any(walb0s(fidx:lidx) < 0)) THEN
           albarr(i) = newalb
           IF (all(walb0s (fidx:lidx) /= newalb)) THEN
              albadj(fidx:lidx) = newalb - walb0s(fidx:lidx)
           ENDIF
            IF (scnwrt) WRITE(*,*) modulename//': alb < 0', &
            sum(albadj(fidx:lidx))/(lidx-fidx+1),sum(walb0s(fidx:lidx))/(lidx-fidx+1)
            walb0s(fidx:lidx) = newalb
         ENDIF
      ENDDO
    ENDIF
  ENDIF
  ! === Call LIDORT, polarization correction, and additional wf calc =====
  IF (use_effcrs) THEN
     CALL LIDORT_PROF_ENV(do_ozwf, do_albwf, do_tmpwf, do_o3shi,&
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf,&
          do_sprswf, do_so2zwf, do_pslwf, ozvary,&
          ns, waves, nalbwf,walb0s, wfc0s, & 
          maxoth, o3shi, nlay, ozprof, tprof,&
          nostk,rio, errstat)
  ELSE
     CALL HRES_RADCALC_ENV(do_ozwf, do_albwf, do_tmpwf, do_o3shi,&
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
          do_sprswf,do_so2zwf, do_pslwf, ozvary,& 
          ns, waves, n0alb, albarr, albpmin, albpmax, nalbwf,walb0s, albnegval,newalb, & 
                     n0wfc, wfcarr, wfcpmin, wfcpmax, wfc0s,  &
          maxoth, o3shi, nlay, ozprof, tprof,&
          nostk,rio, errstat)
  ENDIF
  IF (errstat == pge_errstat_error) THEN
    IF (use_effcrs) THEN
      WRITE(*, *) modulename, ': Errors in calling LIDORT_PROF_ENV!!!'
    ELSE
       WRITE(*, *) modulename, ': Errors in calling HRES_RADCALC_ENV!!!'
    ENDIF
    CALL set_rtmvar (.false., rio)
    RETURN
  ENDIF
  fsimrad(1:ns, 1:nostk)     = rio%rad(1:ns, 1:nostk)
  albwf(1:ns, 1:nalbwf, 1:nostk) = rio%albwf(1:ns, 1:nalbwf,1:nostk)
  o3shiwf(1:ns, 1:nostk)     = rio%o3shiwf(1:ns, 1:nostk)
  cfracwf(1:ns, 1:nostk)     = rio%cfracwf(1:ns, 1:nostk)
  codwf(1:ns, 1:nostk)       = rio%codwf(1:ns, 1:nostk)
  ctpwf(1:ns, 1:nostk)       = rio%ctpwf(1:ns, 1:nostk)
  taodwf(1:ns, 1:nostk)      = rio%taodwf(1:ns, 1:nostk)
  twaewf(1:ns, 1:nostk)      = rio%twaewf(1:ns, 1:nostk)
  saodwf(1:ns, 1:nostk)      = rio%saodwf(1:ns, 1:nostk)
  sprswf(1:ns, 1:nostk)      = rio%sprswf(1:ns, 1:nostk)
  so2zwf(1:ns, 1:nostk)      = rio%so2zwf(1:ns, 1:nostk)
  ozwf(1:ns, 1:nlay, 1:nostk)  = rio%ozwf(1:ns, 1:nlay, 1:nostk)
  tmpwf(1:ns, 1:nlay, 1:nostk) = rio%tmpwf(1:ns, 1:nlay, 1:nostk)
  simq(1:ns)                   = rio%radq(1:ns)
  q_ozwf(1:ns, 1:nlay)         = rio%qozwf(1:ns, 1:nlay)
  rio%wav(1:ns)                = waves(1:ns)
  rio%alb(1:ns)                = walb0s(1:ns)

  IF (do_debug_out) THEN
     funit = 11
     IF (.NOT. do_rtm_pca) THEN
        IF (.NOT. use_effcrs) THEN
          funit = 12
        ELSE
          funit = 13
        ENDIF
     ENDIF
     !CALL debug_rtm (funit, rio, ccrs)
  ENDIF
  CALL set_rtmvar (.false., rio)

  !xliu (02/01/2007): correct radiances based on ozone weighting function to deal with negative ozone values
  IF (negval) THEN
     IF (scnwrt) print * , 'pseudo_mode: negaval' 
     DO i = 1, nlay 
        IF (ozadj(i) > 0) THEN
           fsimrad(1:ns, 1:nostk) = fsimrad(1:ns, 1:nostk) - ozadj(i) * ozwf(1:ns, i, 1:nostk) 
        ENDIF
     ENDDO
  ENDIF

 !xliu (12/11/2014): correct radiances based on SO2 to deal with negative ozone values
  DO k = 1, ngas
    i = fgasidxs(k)
    IF (i > 0) THEN
      j = mask_fitvar_rad(i)
      tracegas(k, 4) = fitvar_rad(j)
      IF ( gasidxs(k) == so2_idx .AND. so2negval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * (1.0d0 + database(so2_idx,refidx(1:ns)) * so2adj)
         tracegas(k, 4) = newso2 - so2adj; fitvar_rad(j) = tracegas(k, 4)
      ENDIF
       IF ( gasidxs(k) == so2v_idx .AND. so2vnegval) THEN
        fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * (1.0d0 + database(so2v_idx,refidx(1:ns)) * so2vadj)
        tracegas(k, 4) = newso2 - so2vadj; fitvar_rad(j) = tracegas(k, 4)
      ENDIF
       IF ( gasidxs(k) == no2_t1_idx .AND. no2negval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * ( 1.0d0 + database(no2_t1_idx,refidx(1:ns)) * no2adj )
         tracegas(k, 4) = newno2 - no2adj; fitvar_rad(j) = tracegas(k, 4)
      ENDIF
       IF ( gasidxs(k) == hcho_idx .AND. hchonegval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * (1.0d0 + database(hcho_idx,refidx(1:ns)) * hchoadj)
         tracegas(k, 4) = newhcho - hchoadj; fitvar_rad(j) = tracegas(k, 4)
      ENDIF
       IF ( gasidxs(k) == bro_idx .AND. bronegval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * ( 1.d0 + database(bro_idx,refidx(1:ns)) * broadj)
         tracegas(k, 4) = newbro - broadj; fitvar_rad(j) = tracegas(k, 4)
      ENDIF
      IF ( gasidxs(k) == o2o2_idx .AND. o4negval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * ( 1.d0 + database(o2o2_idx,refidx(1:ns)) * o4adj)
         tracegas(k, 4) = newo4 - o4adj; fitvar_rad(j) = tracegas(k, 4)
         IF (scnwrt) WRITE(*,*) 'pseudo_model: negative o4'
      ENDIF
      IF ( gasidxs(k) == o2_idx .AND. o2negval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * ( 1.d0 + database(o2_idx,refidx(1:ns)) * o2adj)
         tracegas(k, 4) = newo2 - o2adj; fitvar_rad(j) = tracegas(k, 4)
         IF (scnwrt) WRITE(*,*) 'pseudo_model: negative o2'
      ENDIF
      IF ( gasidxs(k) == h2o_idx .AND. h2onegval) THEN
         fsimrad(1:ns, 1) = fsimrad(1:ns, 1) * ( 1.d0 + database(h2o_idx,refidx(1:ns)) * h2oadj)
         tracegas(k, 4) = newh2o - h2oadj; fitvar_rad(j) = tracegas(k, 4)
         IF (scnwrt) WRITE(*,*) 'pseudo_model: negative h2o'
      ENDIF
    ENDIF
  ENDDO
  IF (albnegval) THEN
    IF (.NOT. vary_sfcalb) THEN
      DO i = 1, n0alb
        IF (albadj(i) > 0 ) THEN
          !fitspec(albpmin(i):albpmax(i)) = fitspec(albpmin(i):albpmax(i)) &
          ! + albadj(i) * albwf(albpmin(i):albpmax(i), 1, 1)
          fsimrad(albpmin(i):albpmax(i),1) = fsimrad(albpmin(i):albpmax(i),1) &
           - albadj(i) * albwf(albpmin(i):albpmax(i), 1, 1)
          fitvar_rad(idx0alb(i)) = newalb !- albadj(i) !Junsung: deleted (03/12/2025) for gler, this should be used for omler
          fitvar_rad(idx0alb(i)) = albadj(i)!newalb - albadj(i) !Junsung: added (03/12/2025) for gler
          !fitvar_rad(idx0alb(i)) = newalb - albadj(i) !Junsung: added (03/12/2025) for gler
          !fitvar(28) = albadj(1) !Junsung: added (03/12/2025) for gler
        ENDIF
      ENDDO
    ELSE
      DO i = 1, nalb
         j = albidx - 1 + i
         READ(fitvar_rad_str(j)(4:5), '(I2)') albord
         fidx=albfpix(i); lidx=alblpix(i)
         IF (any(albadj(fidx:lidx) > 0)) THEN
            ! correction spectrum
            IF ( (.NOT. is_albspcvar(i) .and. albord == 0) .or. &
                (is_albspcvar(i) .and. use_albeofs .and. albord == 1) .or. &
                (is_albspcvar(i) .and. .NOT. use_albeofs .and. albord == 0)) THEN
            ! fitspec(fidx:lidx) = fitspec(fidx:lidx) + albadj(fidx:lidx)
            ! *albwf(fidx:lidx,1,1)
             fsimrad(fidx:lidx,1) = fsimrad(fidx:lidx,1) - albadj(fidx:lidx)*albwf(fidx:lidx,1,1)
            ENDIF
            ! reset state vector
            IF (.NOT. is_albspcvar(i) .and. albord == 0) THEN
                fitvar_rad(j) = newalb
            ELSE
                fitvar_rad(j) = fitvar_rad_apriori(j)
            ENDIF
         ENDIF
      ENDDO
    ENDIF
  ENDIF
  simrad(1:ns) = fsimrad(1:ns, 1)
  IF (nostk > 1) THEN
    Qsimrad (1:ns) = fsimrad(1:ns, 2)
    direc=1.0
    WHERE (Qsimrad(1:ns) < 0_dp)
       direc = -1.0
    END WHERE
  ENDIF
  !IF (radcalwrt .AND. do_simu) THEN
  !   albwf(1:ns, 1) = 0.0
  !   o3shiwf(1:ns, 1) = 0.0
  !   ozwf(1:ns, :, 1) = 0.0
  !   database(so2_idx, refidx(1:ns)) = 0.0
  !   database(bro_idx, refidx(1:ns)) = 0.0
  !   database(bro2_idx, refidx(1:ns)) = 0.0
  !   database(hcho_idx, refidx(1:ns)) = 0.0
  !   database(no2_t1_idx, refidx(1:ns)) = 0.0    
  !ENDIF

  ! get dlnI/dx
  IF (use_lograd) THEN
     DO j=1, nostk
     IF (do_o3shi)   o3shiwf(1:ns, j) = o3shiwf(1:ns, j) / fsimrad(1:ns, j)
     IF (do_codwf)   codwf(1:ns, j)   = codwf(1:ns, j)   / fsimrad(1:ns, j)
     IF (do_sprswf)  sprswf(1:ns, j)  = sprswf(1:ns, j)  / fsimrad(1:ns, j)
     IF (do_so2zwf)  so2zwf(1:ns, j)  = so2zwf(1:ns, j)  / fsimrad(1:ns, j)
     IF (do_ctpwf)   ctpwf(1:ns, j)   = ctpwf(1:ns, j)   / fsimrad(1:ns, j)
     IF (do_cfracwf) cfracwf(1:ns, j) = cfracwf(1:ns, j) / fsimrad(1:ns, j)
     IF (do_taodwf)  taodwf(1:ns, j)  = taodwf(1:ns, j)  / fsimrad(1:ns, j)
     IF (do_saodwf)  saodwf(1:ns, j)  = saodwf(1:ns, j)  / fsimrad(1:ns, j)
     IF (do_twaewf)  twaewf(1:ns, j)  = twaewf(1:ns, j)  / fsimrad(1:ns, j)

     IF (do_albwf)THEN   
       DO i = 1 , nalbwf
         albwf(1:ns, i, j)   = albwf(1:ns, i, j)   / fsimrad(1:ns, j)   
       ENDDO
     ENDIF
     IF (do_ozwf) THEN
       DO i = stlay, endlay
         ozwf (:, i, j) = ozwf(:, i, j) / fsimrad(1:ns, j)
         IF (nt_fit > 0) tmpwf(:, i, j) = tmpwf(:, i, j) / fsimrad(1:ns, j)
       END DO
     ENDIF
     ENDDO
     IF (nostk > 1) Qsimrad(1:ns)= LOG(ABS(Qsimrad(1:ns)))
     simrad(1:ns) = LOG(simrad(1:ns))           ! get dlnI
  END IF

  IF (.NOT. do_albwf) albwf(:,:, :) = 0.0D0
  IF (.NOT. do_cfracwf) cfracwf(:,:) = 0.0D0

  ! correct for linear/quardratic wavelength dependent in albedo
  albothwf = 0.0
  DO i = 1, nalb
    j = albidx + i - 1
    READ(fitvar_rad_str(j)(4:5), '(I2)') albord
    IF (albord == 0 .AND. .NOT. is_albspcvar(i)) CYCLE
    fidx=albfpix(i); lidx=alblpix(i)
    delpos(fidx:lidx) =waves(fidx:lidx) -  SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
    IF (.NOT. is_albspcvar(i)) THEN
      albothwf(fidx:lidx,albord) = albwf(fidx:lidx, 1,1)*(delpos(fidx:lidx))**albord
      IF (.NOT. vary_sfcalb) &
         simrad(fidx:lidx) = simrad(fidx:lidx) +  albothwf(fidx:lidx, albord) *fitvar_rad(j)
      IF (nostk > 1) Qsimrad(fidx:lidx)= Qsimrad(fidx:lidx)+ &
         (albwf(fidx:lidx,1, 2)*(delpos(fidx:lidx))**albord)*fitvar_rad(j)
    ELSE ! vary_sfcalb is turn on and hence not correct simrad
      IF (do_brdf) THEN
         ! re-check
         albothwf(fidx:lidx, albord) = albwf(fidx:lidx,albord,1)
      ELSE IF (use_albeofs) THEN
        IF (nactalbspc == 1) THEN ! snow/ice/water
          albothwf(fidx:lidx, albord) = albwf(fidx:lidx, 1, 1) *albspcs(fidx:lidx, 0) * (delpos(fidx:lidx))**(albord-1)
        ELSE
          albothwf(fidx:lidx, albord) = albwf(fidx:lidx, 1, 1) *albspcs(fidx:lidx, albord)
        ENDIF
      ELSE
        albothwf(fidx:lidx, albord) = albwf(fidx:lidx, 1, 1) *albspcs(fidx:lidx, 0) * (delpos(fidx:lidx))**albord
      ENDIF
    ENDIF
  ENDDO

  IF (nwfc > 0) THEN
     wfcothwf = 0.0
     DO i = 1, nwfc
        j = wfcidx + i - 1
        READ(fitvar_rad_str(j)(4:4), '(I1)') wfcord
        IF (wfcord == 0) CYCLE
        fidx=wfcfpix(i); lidx=wfclpix(i)
        wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
        wfcothwf(fidx:lidx,wfcord) = cfracwf(fidx:lidx, 1)*(waves(fidx:lidx) - wavavg)**wfcord 
        IF (.NOT. vary_sfcalb) & 
        simrad(fidx:lidx) = simrad(fidx:lidx) +  wfcothwf(fidx:lidx, wfcord) * fitvar_rad(j)        
     ENDDO
  ENDIF

  IF (radcalwrt .AND. do_simu .AND. .NOT. do_simu_rmring) THEN
     fitspec = actspec_rad(1:ns)
     IF ( use_lograd ) fitspec = LOG( fitspec )
     fitres = fitspec - simrad
     chisq  = SUM((fitres / fitweights(1:ns))**2.0)
     RETURN
  ENDIF

  ! get calibrated reflectance, correcting for ring, undersampling, and trace gases
  ! calculate weighting function for shift parameter if required, just need to do it once
  IF (num_iter == 0) THEN 
     do_shiwf = .TRUE.
  ELSE
     do_shiwf = .FALSE.
  ENDIF
  CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, simrad, fitspec, errstat)
  IF (errstat == pge_errstat_error) THEN
     WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'; RETURN
  ENDIF
  
  ! commented out by zcai@2016m0630 
  IF (num_iter == 0) THEN
     CALL UV1_SPIKE_DETECT(ns, fitspec, simrad, nsaa_spike)
  ENDIF
  
  ! get residual between measured and simulated reflectance
  fitres = fitspec - simrad
 ! get residual between measured and simulated stokes fraction
  ! *as we didn't model other trace gas, so correction might be taken later
  IF (n_gome_q > 0) THEN
     IF (use_origin_q) THEN
        fitq(1:n_gome_q) = gome_q(2,stkidx(1:n_gome_q))
     ELSE
        fitq = stokfrac(1:n_stokfrac)
        CALL hres_stkwf_inter_convol2(n_stokfrac, stokwaves, fitq, errstat)
     ENDIF
  ENDIF
  IF (polcorr == 0) THEN
     !fitq = stokfrac(1:ns)
     !CALL hres_stkwf_inter_convol2(ns, waves, fitq, errstat)
     CALL hres_stkwf_inter_convol2(ns, waves, simq, errstat)
     !simqcor = EXP(Qsimrad)/EXP(simrad)*direc(:)
     !IF(ANY(simqcor <-1 .OR. simqcor >1) .OR. ANY(simq <-1 .OR. simq >1)) THEN
     !   write(70,*) sza,vza
     !   do i=1, ns
     !     write(70,'(4d16.7)') simq(i), simqcor(i),
     !     EXP(Qsimrad(i)),EXP(simrad(i))
     !   enddo
     !ENDIF
     !fitqres = fitq-simqcor
     fitqres= fitq-simq
  ENDIF
  IF (polcorr == 6 .AND. n_gome_q > 0) THEN
     !fitq(1:n_gome_q) = gome_q(2,stkidx(1:n_gome_q))
     fitqres(1:n_gome_q) = fitq(1:n_gome_q) - simq(1:n_gome_q)
  ENDIF

  ! compute chi-square difference
  IF (.NOT.do_sy_diagonal) THEN 
     Sy = merr_covar(1:ns, 1:ns) 
     CALL sq_matrix_invert(Sy, ns, Sy_inv)
     y1(1:ns, 1) = fitres
     chi = MATMUL(MATMUL(TRANSPOSE(y1), sy_inv),y1)
     chisq = chi(1, 1)
  ELSE
     chisq  = SUM((fitres / fitweights(1:ns))**2.0)
  ENDIF

  rms =    SQRT(chisq  / REAL(ns, KIND=dp))
  relrms = 100.D0 * SQRT(SUM(ABS((simrad-fitspec) / fitspec)**2.0) &
       / REAL(ns, KIND=dp))
  IF (scnwrt) THEN
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        IF (.NOT.do_sy_diagonal) THEN
           j = lidx-fidx+1
           ALLOCATE(y1tmp(j, 1), sy_invtmp(j, j)) 
           y1tmp(:, 1)    = y1(fidx:lidx, 1)
           sy_invtmp(:, :)= sy_inv(fidx:lidx, fidx:lidx)
    
           chi = MATMUL(MATMUL(TRANSPOSE(y1tmp),sy_invtmp),y1tmp)
           allchisq(i) = chi(1, 1)
           DEALLOCATE(y1tmp, sy_invtmp)
        ELSE
           allchisq(i) = SUM((fitres(fidx:lidx) / fitweights(fidx:lidx))**2.0)
        ENDIF
        allrms(i)   = SQRT(allchisq(i) / nradpix(i))
        allradrms(i) = 100.D0 * SQRT(SUM(ABS((simrad(fidx:lidx)-fitspec(fidx:lidx)) &
             / fitspec(fidx:lidx))**2) / nradpix(i)) 
        fidx = lidx + 1
     ENDDO

     ! Relative rms difference between calculated and simulated log-radiances
     ! /radiances depending on the flag use_lograd
     IF (use_lograd) THEN     
        simrad1 = EXP(simrad); fitspec1 = EXP(fitspec)
        radrms = 100.D0 * SQRT(SUM(ABS((simrad1-fitspec1) / fitspec1)**2) &
             / REAL(ns, KIND=dp))   

        fidx = 1
        DO i = 1, numwin
           lidx = fidx + nradpix(i) - 1
           allradrms(i) = 100.D0 * SQRT(SUM(ABS((simrad1(fidx:lidx)-fitspec1(fidx:lidx)) &
                / fitspec1(fidx:lidx))**2) / nradpix(i))
           fidx = lidx + 1
        ENDDO
     ELSE 
        radrms = relrms
     ENDIF

    WRITE (*, '(I5, 4(A10,1pd11.3))') num_iter, ' Chi = ', chisq, ' rms = ', rms, &
          ' rms(%) = ', relrms, ' Irms(%) = ', radrms
    fidx=1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) -1
      WRITE (*, '(A, I2, 4(A, 1pd11.3))') 'Win ', i, ': allrms = ', &
        allrms(i), 'allIrms(%) = ', allradrms(i),&
        'Is:',sum(simrad1(fidx:lidx))/nradpix(i),'Im:',sum(fitspec1(fidx:lidx))/nradpix(i)
      fidx = lidx + 1
    ENDDO
  ENDIF

  IF (npix_fitted == 0) THEN
     min_ssa_iter = 2
  ELSE
     min_ssa_iter = 1
  ENDIF

  ! South Atlantic Anomaly Correction, use 10 channel 2 pixels as reference because they are
  ! less subject to South Atlantic Anomaly
  IF (num_iter == min_ssa_iter .AND. num_iter <=min_ssa_iter+2 .AND. .NOT. refl_only .AND. &
       band_selectors(1) == 1 .AND. saa_flag  ) THEN 
     CALL spike_detect_correct(nradpix(1)+10, fitspec(1:nradpix(1)+10), simrad(1:nradpix(1)+10))  
  ENDIF
  fitres = fitspec - simrad

  IF (.NOT. refl_only) THEN 
     dyda = 0.0

     ! albedo weighting function
     DO i = 1, nfalb
        j = albfidx + i - 1
        k = mask_fitvar_rad(j) 
        fidx = albfpix(k - albidx + 1); lidx=alblpix(k - albidx + 1)
           
        READ(fitvar_rad_str(k)(4:5), '(I1)') albord
        IF ( do_alb_longwav .AND. numwin == 2 .AND. .NOT. (which_cld == 5 .AND. use_retalb)) THEN
           ! Find the index of alb_swav
           albeidx = MINVAL(MINLOC(waves(fidx:lidx), mask=(waves(fidx:lidx) >=alb_ewav)))  
           albeidx = albeidx + lidx 
           albsidx = MINVAL(MAXLOC(waves(fidx:lidx), mask=(waves(fidx:lidx) <=alb_swav)))
           albsidx = albsidx + fidx 
           IF (albsidx < fidx) albsidx = fidx
           IF (albeidx > lidx) albeidx = lidx
        ENDIF
        IF (albord == 0) THEN
           !IF ( do_alb_longwav .AND. fitvar_rad_str(k) == '2ba0' ) THEN
           !   dyda(albsidx:lidx, j) = albwf(albsidx:lidx, 1)
           !ELSE 
           dyda(fidx:lidx, j)=albwf(fidx:lidx, 1, 1)
           !ENDIF
           !dyda(albsidx:lidx, j) = albwf(albsidx:lidx, 1)
        ELSE
           IF (do_alb_longwav .AND. numwin == 2 .AND. .NOT. (which_cld == 5 .AND. use_retalb) &
                .AND. (fitvar_rad_str(k) == '2ba1' .OR. fitvar_rad_str(k) == '2ba2') .AND.     &
                albsidx > 0 .AND. albeidx > 0 .AND.  albsidx >= fidx .AND. albeidx <= lidx) THEN
              WRITE(*, *) 'Using longer wavelength dyda only '
              !IF (albeidx > 0 .AND. albeidx < lidx) dyda(albeidx:lidx, j) = albothwf(albeidx:lidx, albord)
              !IF (albsidx > 0 .AND. albsidx > fidx) dyda(fidx:albsidx, j) = albothwf(fidx:albsidx, albord)
              dyda(albsidx:albeidx, j) = albothwf(albsidx:albeidx, albord) 
              !dyda(fidx:lidx, j)=albothwf(fidx:lidx, albord)
           ELSE 
              dyda(fidx:lidx, j)=albothwf(fidx:lidx, albord)
           ENDIF
           !dyda(albsidx:lidx, j) = albothwf(albsidx:lidx, albord)
        ENDIF

        ! xliu: 07/01/2010, compute aerosol index (defined with relative to 20 nm distance)
        ! defined similar to TOMS aerosol index
        IF (albord == 1) THEN
           the_ai = (dyda(lidx, j) - dyda(fidx, j)) * fitvar(j) * 100. / (waves(lidx)-waves(fidx)) * 20.
           IF ( radcalwrt .AND. do_simu) THEN !JBAK
               the_ai = (dyda(lidx, j) - dyda(fidx, j))*100 /(waves(lidx)-waves(fidx)) * 20.
           ENDIF
        ENDIF
     ENDDO

     ! wavelength-dependent cloud fraction weighting function
     DO i = 1, nfwfc
        j = wfcfidx + i - 1
        k = mask_fitvar_rad(j) 
        fidx = wfcfpix(k -wfcidx + 1); lidx=wfclpix(k - wfcidx + 1)

        READ(fitvar_rad_str(k)(4:4), '(I1)') wfcord

        IF ( do_alb_longwav .AND. numwin == 2 .AND. .NOT. (which_cld == 5 .AND. use_retalb)) THEN
           ! Find the index of alb_swav
           albeidx = MINVAL(MINLOC(waves(fidx:lidx), mask=(waves(fidx:lidx) >=alb_ewav)))
           albeidx = albeidx + lidx 
           albsidx = MINVAL(MAXLOC(waves(fidx:lidx), mask=(waves(fidx:lidx) <=alb_swav)))
           albsidx = albsidx + fidx 
           IF (albsidx < fidx) albsidx = fidx
           IF (albeidx > lidx) albeidx = lidx      
        ENDIF

        IF (wfcord == 0) THEN
           !IF ( do_alb_longwav .AND. fitvar_rad_str(k) == '2fc0' ) THEN
           !   dyda(fidx:lidx, j) = cfracwf(albsidx:lidx, 1)
           !ELSE
           dyda(fidx:lidx, j) = cfracwf(fidx:lidx, 1)
           !ENDIF
        ELSE 
           IF (do_alb_longwav .AND. numwin == 2  .AND. .NOT. (which_cld == 5 .AND. use_retalb) &
                .AND.  (fitvar_rad_str(k) == '2fc1' .OR. fitvar_rad_str(k) == '2fc2') .AND. &
                albsidx > 0 .AND. albeidx > 0 .AND.  albsidx >= fidx .AND. albeidx <= lidx) THEN
             ! WRITE(*, *) 'Using longer wavelength dyda only '          
             ! IF (albeidx > 0 .AND. albeidx < lidx) dyda(albeidx:lidx, j) = wfcothwf(albeidx:lidx, albord)
             ! IF (albsidx > 0 .AND. albsidx > fidx) dyda(fidx:albsidx, j) = wfcothwf(fidx:albsidx, albord)
                   dyda(albsidx:albeidx, j) = wfcothwf(albsidx:  albeidx, albord) 
           ELSE 
              dyda(fidx:lidx, j) = wfcothwf(fidx:lidx, wfcord)
           ENDIF
        ENDIF
     ENDDO

     ! ozone profile and temperature weighting function
     dyda(:, ozf_fidx:ozf_lidx) = ozwf(:, stlay:endlay, 1) 
     ! do not use longer waves for ozone, this make no sense 
     !IF (do_alb_longwav .AND. waves(ns) > 330.0 ) THEN
    !    lidx = MINVAL(MINLOC(waves(1:ns), mask=(waves(1:ns) >=330.0)))
    !    IF (lidx > 0 .AND. lidx <= ns) dyda(lidx:ns, ozf_fidx:ozf_lidx) = 0.0D0
     !ENDIF

     IF (nt_fit > 0) THEN
        fidx = stlay + mask_fitvar_rad(tf_fidx)-maxlay-ozp_fidx
        lidx = endlay + mask_fitvar_rad(tf_lidx)-maxlay-ozp_lidx
        dyda(:, tf_fidx:tf_lidx) =tmpwf(1:ns, fidx:lidx, 1) 
     ENDIF

     ! Aerosol/cloud parameters 
     IF (taodfind > 0) dyda(:, taodfind) = taodwf(:, 1)
     IF (saodfind > 0) dyda(:, saodfind) = saodwf(:, 1)
     IF (twaefind > 0) dyda(:, twaefind) = twaewf(:, 1)
     IF (ecfrfind > 0) dyda(:, ecfrfind) = cfracwf(:, 1)
     IF (ecodfind > 0) dyda(:, ecodfind) = codwf(:, 1)
     IF (sprsfind > 0) dyda(:, sprsfind) = sprswf(:, 1)
     IF (so2zfind > 0) dyda(:, so2zfind) = so2zwf(:, 1)
     IF (ectpfind > 0) dyda(:, ectpfind) = ctpwf(:, 1)     

     ! get 1st derivative for calibration and reference parameters
     ! wfs for undersampling, common mode, other gases and shift terms
     fitvar_saved(1:nf) = fitvar(1:nf); do_shiwf = .FALSE.
     DO i = 1, ozf_fidx - 1     
        ! shift parameter, already calculate dR/dS
        IF (mask_fitvar_rad(i) > shift_offset) CYCLE 
        ! check for shift indices
        ridx = database_indices(i) ! indices in the database
        sidx = shift_offset + ridx ! shift indice in fitvar_rad

        !IF (ridx == so2_idx .OR. alb_ewav, alb_swavno2_t1_idx .OR. &
        !   ridx == bro_idx .OR. ridx == bro2_idx .OR. ridx == so2v_idx .OR. ridx == o2o2_idx) THEN
        IF (ridx /= us1_idx .AND. ridx /= us2_idx .AND. ridx /= com_idx .AND. ridx /= com1_idx &
             .AND. ridx /= com2_idx .AND. ridx /= com3_idx .AND. ridx /= fsl_idx .AND. ridx /= rsl_idx) THEN
           ! for trace gases, wfs are just the cross sections (or amf * cross sections)
           dyda(:, i) = -database(ridx, refidx(1:ns))
           IF (.NOT. use_lograd)  dyda(:, i) =  dyda(:, i) * simrad
           ! xliu, 11/01/2011, the following for undersampling is incorrect
           !ELSE  IF (ridx == us1_idx .OR. ridx == us2_idx) THEN
           !   dyda(:, i) = - EXP(fitspec(1:ns)) * database(ridx, refidx(1:ns))  
        ELSE
           sim1 = simrad
           fitvar(i) = fitvar_saved(i) * 1.001

           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim1, meas1, errstat)
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF
           stop

           sim2 = simrad
           fitvar(i) = fitvar_saved(i) * 0.999
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim2, meas2, errstat)  
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF
           dyda(:, i) = -(meas1 - sim1 - meas2 + sim2) / (0.002 * fitvar_saved(i))
           fitvar(i) = fitvar_saved(i)
           
        ENDIF
        ! check for shifting
        IF (rmask_fitvar_rad(sidx) > 0) THEN
           dyda(:, rmask_fitvar_rad(sidx)) = dyda(:, i) * database_shiwf(ridx, refidx(1:ns))  
        ENDIF
     ENDDO

     ! Weighting function when use fit_atanring
     ! 1. analytical  2. finite difference (use 2 to validate 1)
     IF (fit_atanring) THEN
        do_shiwf = .FALSE.
        DO i = rnfind(1, 1), rnfind(1, 3)
           sim1 = simrad
           fitvar(i) = fitvar_saved(i) * 1.001           
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim1, meas1, errstat)
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF
           sim2 = simrad
           fitvar(i) = fitvar_saved(i) * 0.999
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim2, meas2, errstat)  
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF
           dyda(:, i) = -(meas1 - sim1 - meas2 + sim2) / (0.002 * fitvar_saved(i))
           fitvar(i) = fitvar_saved(i)
        ENDDO

        !! Analytical weighting function
        !delpos = (fitwavs-fitvar(rnfind(1, 2))) / fitvar(rnfind(1, 3))
        !dyda(:, rnfind(1, 1)) = database(ring_idx, refidx(1:ns)) * (atan(delpos)+1.54223)
        !dyda(:, rnfind(1, 2)) = -database(ring_idx, refidx(1:ns)) / (1.0 + delpos**2) * &
        !     fitvar(rnfind(1, 1)) / fitvar(rnfind(1, 3))
        !dyda(:, rnfind(1, 3)) = dyda(:, rnfind(1, 2)) * delpos             
     ENDIF

     ! wfs for ozcrs shift, slit_shift, wavelength shift, Ring effect, degradation correction,
     ! internal scattering in irradiance, internal scattering in radiance
     DO ig = 1, nothgrp
        nord = 0
        IF (ig == 1) THEN            ! ozone cross section  
           nord = nos; tmpind = osind; tmpfind = osfind; tmpwins = oswins
           temporwf = -o3shiwf(1:ns, 1)
        ELSE IF ( ig == 2) THEN      ! Radiance/Irradiance Slit Difference
           nord = nsl; tmpind = slind; tmpfind = slfind; tmpwins = slwins
           temporwf = -slwf(1:ns)
        ELSE IF (ig == 3)  THEN      ! Radiance/irradince wavelength shift
           nord = nsh; tmpind = shind; tmpfind = shfind; tmpwins = shwins
        ELSE IF (ig == 4)  THEN      ! Ring effect, default: fit_atanring =false
           IF (fit_atanring) CYCLE
           nord = nrn; tmpind = rnind; tmpfind = rnfind; tmpwins = rnwins
           temporwf = -database(ring_idx, refidx(1:ns))
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
        ELSE IF (ig == 5 ) THEN      ! Degradation correction
           nord = ndc; tmpind = dcind; tmpfind = dcfind; tmpwins = dcwins
        ELSE IF (ig == 6)  THEN      ! Internal scattering in irradinace
           nord = nis; tmpind = isind; tmpfind = isfind; tmpwins = iswins
        ELSE IF (ig == 7) THEN       ! Internal scattering in radiance
           nord = nir; tmpind = irind; tmpfind = irfind; tmpwins = irwins
        ELSE IF (ig == 8 ) THEN                                                         
           nord = np1; tmpind = p1ind; tmpfind = p1find; tmpwins = p1wins ; slit_idx=psl_fpos(1)
           temporwf(1:ns) = - database_pslwf(refidx(1:ns), 1)
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
        ELSE IF (ig == 9 ) THEN 
           nord = np2; tmpind = p2ind; tmpfind = p2find; tmpwins = p2wins ; slit_idx=psl_fpos(2)
           temporwf(1:ns) = - database_pslwf(refidx(1:ns), 2)
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
        ELSE IF (ig == 10 ) THEN 
           nord = np3; tmpind = p3ind; tmpfind = p3find; tmpwins = p3wins ; slit_idx=psl_fpos(3)
           temporwf(1:ns) = - database_pslwf(refidx(1:ns), 3)
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
        ELSE IF (ig == 11 ) THEN 
           nord = ncm; tmpind = cmind; tmpfind = cmfind; tmpwins = cmwins
           temporwf(1:ns) = database_cmwf(refidx(1:ns)) 
          !if (num_iter > 0 ) print *, num_iter, fitvar_rad(cmind(1,1)),fitvar_rad(cmind(2,1))
           !print * , temporwf(ns-10:ns)
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
           !if (ig ==10  ) print *, fitvar(tmpfind(1:2, 1))
        ENDIF
        ! Use finite differences to get the zero-order weighting functions
        IF (nord > 0 .AND. ig >= 3 .AND. ig <= nothgrp .AND. ig /= 4 .AND. ( ig < 8 .or. ig > 11) ) THEN
           temporwf = 0.0
           IF (do_subfit) THEN
              swin = tmpwins(1, 1); ewin = tmpwins(1, 2)
           ELSE
              swin = 1; ewin = 1
           ENDIF

           IF (ig == 3) THEN  ! irradiance/radiance wavelength shift
              sim1 = simrad
              fitvar(tmpfind(swin:ewin, 1)) = 0.001 
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim1, meas1, errstat) 
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              sim2 = simrad
              fitvar(tmpfind(swin:ewin, 1)) = -0.001 
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim2, meas2, errstat)
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              temporwf = -(meas1 -  sim1 - meas2 + sim2) / 0.002
              fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1))
           ELSE
              sim1 = simrad
              fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)) * 1.001 
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim1, meas1, errstat) 
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              sim2 = simrad
              fitvar(tmpfind(swin:ewin, 1))  = fitvar_saved(tmpfind(swin:ewin, 1)) * 0.999
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, sim2, meas2, errstat)
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              fitvar(tmpfind(swin:ewin, 1))  = fitvar_saved(tmpfind(swin:ewin, 1))
              IF (swin == ewin) THEN
                 IF (fitvar_saved(tmpfind(swin, 1)) /= 0.0) THEN
                    temporwf = -(meas1 - sim1 - meas2 + sim2) / & 
                                (0.002 * fitvar_saved(tmpfind(swin, 1)))
                 ELSE
                    temporwf = 0.0d0
                 ENDIF
              ELSE
                 IF (swin == 1) THEN
                    fidx = 1
                 ELSE
                    fidx = SUM(nradpix(1: swin-1)) + 1
                 ENDIF

                 DO iw = swin, ewin
                    lidx = fidx + nradpix(iw) - 1
                    IF (fitvar_saved(tmpfind(iw, 1)) /= 0.0) THEN
                       temporwf(fidx:lidx) = -(meas1(fidx:lidx) - sim1(fidx:lidx) - &
                             meas2(fidx:lidx) + sim2(fidx:lidx)) / &
                            (0.002 * fitvar_saved(tmpfind(iw, 1)))
                    ELSE
                       temporwf(fidx:lidx) = 0.d0
                    ENDIF
                    fidx = lidx + 1
                 ENDDO
              ENDIF
           ENDIF
        ENDIF

        IF (nord > 0) THEN
           IF (do_subfit) THEN
              fidx = 1
              DO j = 1, numwin
                 lidx = fidx + nradpix(j) - 1
                 IF ( tmpfind(j, 1) > 0 ) dyda(fidx:lidx, tmpfind(j, 1)) = temporwf(fidx:lidx)

                 delpos(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 0.5
                 DO i = 2, nord
                    IF ( tmpfind(j, i) > 0 ) THEN
                       dyda(fidx:lidx, tmpfind(j, i)) = dyda(fidx:lidx, tmpfind(j, i-1)) * delpos(fidx:lidx)
                    ENDIF
                 ENDDO
                 fidx = lidx + 1
              ENDDO
           ELSE 
              IF (tmpwins(1, 1) == 1) THEN
                 fidx = 1
              ELSE
                 fidx = SUM(nradpix(1: tmpwins(1, 1) - 1)) + 1
              ENDIF
              lidx = SUM(nradpix(1: tmpwins(1, 2)))

              dyda(fidx:lidx, tmpfind(1, 1)) = temporwf(fidx:lidx)
              delpos(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 0.5

              DO i = 2, nord
                 IF ( tmpfind(1, i) > 0 ) THEN
                    dyda(fidx:lidx, tmpfind(1, i)) = dyda(fidx:lidx, tmpfind(1, i-1)) * delpos(fidx:lidx)
                 ENDIF
              ENDDO
           ENDIF
        ENDIF
     ENDDO

     ! Remove longer wavelengths wfs except for alb and cloud
     idx330 = maxval(minloc(waves(1:ns), mask=(waves(1:ns)>=330.0)))
     !IF (idx330  > 0 .AND. idx330  <= ns) dyda(idx330:ns, ozf_fidx:ozf_lidx) = 0.0D0
     IF ( do_alb_longwav) then 
      !do i = 1, nf
        !if ( (i < albfidx .or. i > albfidx + nfalb - 1) .and. i .ne. ecfrfind ) then  
        !if ( (i < albfidx .or. i > albfidx + nfalb - 1) )then 
        !if ( (i < albfidx +nfalb-1 .or. i > albfidx + nfalb - 1) )then 
        !   dyda(idx330:ns,i) = 0.0d0 
        !endif
      !enddo 
      ! do not use longer waves for ozone, this make no sense 
      IF (idx330  > 0 .AND. idx330  <= ns) dyda(idx330:ns, ozf_fidx:ozf_lidx) = 0.0D0
      !IF (idx330  > 0 .AND. idx330  <= ns) dyda(idx330:ns, :) = 0.0D0
     ENDIF

     WRITE(92, *) ns, nf
     DO i = 1, ns
        WRITE(92, '(f10.4, 80d14.6)') fitwavs(i), fitspec(i), simrad(i), walb0s(i), dyda(i, 1:nf)
     ENDDO
     CLOSE(92) !; stop 1
     !stop

     DO i = 1, nf  
        dyda(:, i) = dyda(:, i) / fitweights(1:ns)
     END DO 
     ! finnally obtain the new spectrum to be fitted in the GSVD
     gspec(1:ns) = fitres(1:ns) / fitweights(1:ns)
     ! Restore the unperturbated fitting variables
     fitvar(1:nf) = fitvar_saved(1:nf)
     fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar(1:nf)  
  END IF

  RETURN

END SUBROUTINE pseudo_model

SUBROUTINE HRES_RADCALC_ENV (do_ozwf, do_albwf, do_tmpwf, do_o3shi, &
     do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, & 
     do_sprswf, do_so2zwf, do_pslwf, ozvary,&
     nw, waves,&
     n0alb, albarr, albpmin, albpmax, nalbwf, walb0s, albnegval, newalb, &
     n0wfc, wfcarr, wfcpmin, wfcpmax, wfc0s, &
     nos, o3shi, nl, ozprof, tprof, &
     nostk, rio, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : numwin, &
       fitvar_rad_str, fitvar_rad, rmask_fitvar_rad
  USE ozprof_data_module,     ONLY : radcwav, ncalcp,  &
       nalb, albidx, albfpix, alblpix, albfpix_r, alblpix_r, &
       nwfc, wfcidx, wfcmin, wfcmax, wfcfpix, wfclpix, &
       which_cld, do_alb_longwav, use_prefitalb, pf2ba0,&
       pf2ba1, pf2fc0, pf2fc1, & 
       is_albspcvar, use_albeofs, nactalbspc, albspcs_hres,  &
       do_rtm_pca, rtm_outputs, vary_sfcalb, do_brdf
  USE cloud_data_module, ONLY : new_cfrac, new_alb, n_newalb, &
       n_newwfc, use_retalb, avgwav_cld
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN) :: nw, nl, nos, n0alb, nostk, n0wfc, nalbwf
  LOGICAL, INTENT(IN) :: do_ozwf, do_albwf, do_tmpwf, do_o3shi, do_taodwf, &
       do_twaewf, do_saodwf, do_cfracwf, do_codwf, do_ctpwf, do_sprswf, do_so2zwf, do_pslwf
  INTEGER, DIMENSION(n0alb), INTENT(IN)                :: albpmax, albpmin ! for OMI grids
  INTEGER, DIMENSION(n0wfc), INTENT(IN)                :: wfcpmax, wfcpmin
  LOGICAL, DIMENSION(nl), INTENT(IN)                   :: ozvary
  REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)           :: waves, walb0s, wfc0s
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)   :: o3shi
  REAL (KIND=dp), DIMENSION(n0alb), INTENT(IN)         :: albarr 
  REAL (KIND=dp), DIMENSION(n0wfc), INTENT(IN)         :: wfcarr 
  REAL (KIND=dp), DIMENSION(nl),    INTENT(IN)         :: ozprof, tprof
  LOGICAL,                          INTENT(IN)         :: albnegval
  REAL (KIND=dp), INTENT(IN)                           :: newalb
  TYPE (rtm_outputs),               INTENT(OUT)        :: rio
  INTEGER, INTENT(OUT)                                 :: errstat
  ! =======================
  ! Local variables
  ! =======================
  LOGICAL  :: no2alb0
  INTEGER :: nw0, i, j, k, m, n, fidx, lidx, fidx0, lidx0, albord, wfcord
  INTEGER, DIMENSION(n0wfc)                   :: wfcpmax0, wfcpmin0
  REAL (KIND=dp)                              :: wavavg
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE   :: waves0, walb0s0, wfc0s0
  CHARACTER(LEN=1)  :: ordchar
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  !CHARACTER(16), PARAMETER :: modulename = 'HRES_RADCALC_ENV'
  errstat = pge_errstat_ok
  nw0 = MAXVAL([nw, ncalcp])
  IF (allocated(waves0)) deallocate (waves0, walb0s0, wfc0s0)
  allocate (waves0(nw0), walb0s0(nw0), wfc0s0(nw0))
  waves0 = 0.0 ; walb0s0 = 0.0
  waves0(1:ncalcp) = radcwav(1:ncalcp)

  ! Need to find indices of boundaries for using different surface albedos/cloud fractions
  ! no2alb0 this make no change 
  no2alb0 = .FALSE.
  k = 0; fidx0 = 1; lidx0=0
  DO i = 1, nalb
     j = albidx - 1 + i
     READ(fitvar_rad_str(j)(4:5), '(I2)') albord
     fidx=albfpix(i); lidx=alblpix(i)
     fidx0=albfpix_r(i); lidx0 = alblpix_r(i)
     wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx)) 
     IF (.NOT. is_albspcvar(i)) THEN
       IF (albord == 0) THEN
         k = k + 1
         walb0s0(fidx0:lidx0) = albarr(k)
         IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j)(1:4) == '2ba0' &
             .AND. rmask_fitvar_rad(j) <= 0) THEN
           DO n = i, nalb
              m = albidx -1 + n
              IF (fitvar_rad_str(m) == '3ba0'.AND. rmask_fitvar_rad(m) > 0) EXIT 
           ENDDO
           IF (n > nalb) THEN
              WRITE(*, *)  'HRES_RADCALC_ENV: ERROR 2ba0 is not filled'
              errstat = pge_errstat_error; RETURN   
           ENDIF
           walb0s0(fidx0:lidx0) = fitvar_rad(m)
           no2alb0 = .TRUE.
           !ELSE IF (do_alb_longwav.AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
           !   walb0s0(fidx0:lidx0) = new_alb(1)
         ELSE IF (use_prefitalb .AND. fitvar_rad_str(j)(1:4) == '2ba0' &
             .AND. rmask_fitvar_rad(j) <= 0 .AND. vary_sfcalb) THEN
           walb0s0(fidx0:lidx0) = pf2ba0
         ENDIF
       ELSE
         IF (vary_sfcalb) THEN
           ! Note use the exact average wavelength as retrieval grid
           ! Get surface albedo for each radiance calculation wavelength
           walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + fitvar_rad(j) * (waves0(fidx0:lidx0) - wavavg)**albord
           IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j)(1:3) == '2ba' &
                .AND. rmask_fitvar_rad(j) <= 0)  THEN
              WRITE(ordchar, '(I1)') albord 
              DO n = i, nalb
                 m = albidx -1 + n
                 IF (fitvar_rad_str(m)(1:4) == '3ba'//ordchar.AND. rmask_fitvar_rad(m) > 0) EXIT 
              ENDDO
              IF (n > nalb) THEN
                 WRITE(*, *)  'HRES_RADCALC_ENV: ERROR 2ba1 is not filled'
                 errstat = pge_errstat_error; RETURN   
              ENDIF
              IF (no2alb0) THEN
                 fidx = albfpix(n)
                 lidx = alblpix(n)
                 wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
              ENDIF
              walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + fitvar_rad(m) * (waves0(fidx0:lidx0) - wavavg)**albord
           ELSE IF (do_alb_longwav.AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
              m = albord + 1
              IF (m > n_newalb) THEN
                 WRITE(*, *)  'HRES_RADCALC_ENV: which_cld==5, alb variables not consist'
              ENDIF
              walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + new_alb(m) * (waves0(fidx0:lidx0) - avgwav_cld)**albord
           ELSE IF (use_prefitalb .AND. fitvar_rad_str(j)(1:3) == '2ba'  &
                .AND. rmask_fitvar_rad(j) <= 0)  THEN
              walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + pf2ba1 * (waves0(fidx0:lidx0) - wavavg)**albord
           ENDIF
         ENDIF ! IF (vary_sfcalb)
       ENDIF ! noth  
     ELSE
      IF (do_brdf) THEN     
        IF (albord > 0) THEN 
          print *, "Using kern_amp here and in m_lidort_env_vv2p7.f90 is wrong, but do_brdf is not working anyway!"
          print *, "fidx0,lidx0 are for UV+vis, while kern_amp is for vis only"
          Surface%kern_amp(fidx0:lidx0, albord) = Surface%kern_amp(fidx0:lidx0,albord)*fitvar_rad(j)
        ENDIF
      ELSE IF (use_albeofs) THEN
        IF (nactalbspc == 1 ) THEN ! Snow/ice/water
          IF (albord == 1) THEN
            walb0s0(fidx0:lidx0) = albspcs_hres(fidx0:lidx0, 0) * fitvar_rad(j)
          ELSE
            walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + (waves0(fidx0:lidx0) - wavavg)**(albord-1) &
            * albspcs_hres(fidx0:lidx0, 0) * fitvar_rad(j)
          ENDIF
        ELSE
          IF (albord == 1)  walb0s0(fidx0:lidx0) = albspcs_hres(fidx0:lidx0, 0)
          walb0s0(fidx0:lidx0) =  walb0s0(fidx0:lidx0) + albspcs_hres(fidx0:lidx0,albord) * fitvar_rad(j)
        ENDIF
      ELSE
        IF (albord == 0) THEN
          walb0s0(fidx0:lidx0) = albspcs_hres(fidx0:lidx0, 0) * fitvar_rad(j)
        ELSE
          walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + (waves0(fidx0:lidx0) - wavavg)**albord &
          * albspcs_hres(fidx0:lidx0, 0) * fitvar_rad(j)
        ENDIF
      ENDIF
    ENDIF
   !WRITE(*,'(5i4,15f8.2)')fidx, lidx, fidx0, lidx0, albord, waves0(fidx0), waves0(lidx0), waves(fidx),waves(lidx),&
  !   walb0s0(fidx0), walb0s0(lidx0), walb0s(fidx),walb0s(lidx), & 
  !   albspcs0(fidx0,albord), albspcs0(lidx0,albord), albspcs(fidx, albord), albspcs(lidx,albord), albmin(i), albmax(i)
  ENDDO

  IF (albnegval) THEN
    DO i = 1, nalb
       j = albidx - 1 + i
       fidx=albfpix(i); lidx=alblpix(i)
       fidx0=albfpix_r(i); lidx0=alblpix_r(i)
       IF (all(walb0s(fidx:lidx) == newalb)) THEN
             walb0s0(fidx0:lidx0) = newalb
       ENDIF
    ENDDO
  ENDIF

  IF (nwfc > 0) THEN 
  no2alb0 = .FALSE.
  k = 0;  wfcpmin0(1) = 1; wfcpmax0(n0wfc) = ncalcp 
  DO i = 1, nwfc
    j = wfcidx - 1 + i
    READ(fitvar_rad_str(j)(4:4), '(I1)') wfcord
    IF (wfcord == 0) THEN
      k = k + 1
      IF (k > 1) wfcpmin0(k)= wfcpmax0(k-1) + 1 
      IF (k < n0wfc) wfcpmax0(k)= MINVAL(MAXLOC(waves0(1:ncalcp), MASK=(waves0(1:ncalcp) &
          >= wfcmin(i) .AND. waves0(1:ncalcp) < wfcmax(i))))  
      fidx0 = wfcpmin0(k); lidx0 = wfcpmax0(k)
      wfc0s0(wfcpmin0(k):wfcpmax0(k)) = wfcarr(k)
      ! If use longwave albedo/cfrac
      IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j) == '2fc0' &
         .AND. rmask_fitvar_rad(j) <= 0 .AND. vary_sfcalb)  THEN
        DO n = i, nwfc
          m = wfcidx -1 + n
          IF (fitvar_rad_str(m) == '3fc0'.AND. rmask_fitvar_rad(m) > 0) EXIT 
        ENDDO
        IF (n > nwfc) THEN
           WRITE(*, *)  'HRES_RADCALC_ENV: ERROR 2fc0 is not filled'
           errstat = pge_errstat_error; RETURN   
        ENDIF
           wfc0s0(fidx0:lidx0) = fitvar_rad(m)
           no2alb0 = .TRUE.
           !ELSE IF (do_alb_longwav.AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
           !   wfc0s0(fidx0:lidx0) = new_cfrac(1)
      ELSE IF (use_prefitalb .AND. fitvar_rad_str(j) == '2fc0' &
           .AND. rmask_fitvar_rad(j) <= 0 .AND. vary_sfcalb) THEN
        wfc0s0(wfcpmin0(k):wfcpmax0(k)) = pf2fc0
      ENDIF
    ELSE
      IF (vary_sfcalb) THEN
        fidx=wfcfpix(i); lidx=wfclpix(i)
        wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx)) 
        wfc0s0(fidx0:lidx0) = wfc0s0(fidx0:lidx0) + fitvar_rad(j) * (waves0(fidx0:lidx0) - wavavg)**wfcord
        ! If use longwave albedo/cfrac
        IF (do_alb_longwav .AND. numwin == 3 .AND. fitvar_rad_str(j)(1:3) == '2fc' &
          .AND. rmask_fitvar_rad(j) <= 0)  THEN
          WRITE(ordchar, '(I1)') wfcord 
          DO n = i, nwfc
            m = wfcidx -1 + n
            IF (fitvar_rad_str(m) == '3fc'//ordchar.AND. rmask_fitvar_rad(m) >= 0) EXIT 
          ENDDO
          IF (n > nwfc) THEN
            WRITE(*, *)  'HRES_RADCALC_ENV: ERROR 2fc1 is not filled'
            errstat = pge_errstat_error; RETURN   
          ENDIF
          IF (no2alb0) THEN
            fidx = wfcfpix(n) ;  lidx = wfclpix(n)
            wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
          ENDIF
          wfc0s0(fidx0:lidx0) = wfc0s0(fidx0:lidx0) + fitvar_rad(m) * (waves0(fidx0:lidx0) - wavavg)**wfcord
        ELSE IF (do_alb_longwav .AND. numwin == 2 .AND. which_cld == 5 .AND. use_retalb) THEN
          m = wfcord + 1
          IF (m > n_newwfc) THEN
            WRITE(*, *) 'HRES_RADCALC_ENV: which_cld==5, cfrac variables not consist'
          ENDIF
            wfc0s0(fidx0:lidx0) = wfc0s0(fidx0:lidx0) + new_cfrac(m) * (waves0(fidx0:lidx0) - avgwav_cld)**wfcord
          ELSE IF (use_prefitalb .AND. fitvar_rad_str(j)(1:3) == '2fc'  &
            .AND. rmask_fitvar_rad(j) <= 0)  THEN
            wfc0s0(fidx0:lidx0) = wfc0s0(fidx0:lidx0) + pf2fc1 * (waves0(fidx0:lidx0) - wavavg)**wfcord
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  ENDIF
  ! Call LIDORT_PROF_ENV on fine wavelength grids
  ! Return raidances and weighting functions on required resolution wavelength grid
  IF (do_rtm_pca) THEN 
    CALL LIDORT_PROF_ENV_PCA (do_ozwf, do_albwf, do_tmpwf, do_o3shi,do_taodwf,      &
       do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, & 
       do_sprswf, do_so2zwf, do_pslwf, ozvary,    &
       nw0, waves0, nalbwf, walb0s0, wfc0s0, nos, o3shi, nl, ozprof, tprof, & 
       nostk, rio,errstat)
  ELSE
    CALL LIDORT_PROF_ENV (do_ozwf, do_albwf, do_tmpwf, do_o3shi,do_taodwf,      &
       do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, & 
       do_sprswf, do_so2zwf, do_pslwf, ozvary,    &
       nw0, waves0, nalbwf, walb0s0, wfc0s0, nos, o3shi, nl, ozprof, tprof, & 
       nostk, rio,errstat)
  ENDIF

  deallocate (waves0, walb0s0, wfc0s0)
  RETURN
END SUBROUTINE HRES_RADCALC_ENV

END MODULE m_pseudo_model
