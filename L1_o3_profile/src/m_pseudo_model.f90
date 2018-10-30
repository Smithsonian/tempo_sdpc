!
module m_pseudo_model
  use m_lidort_utilities, ONLY:get_hres_radcal_waves
  USE m_lidort_master, ONLY:lidort_prof_env
  use m_spectra_reflectance
  use adj_measurement_data, only: uv1_spike_detect, spike_detect_correct

  public pseudo_model
  private

contains

  ! **************************************************************************
  ! Author: Xingu Liu
  ! Date:   July 24, 2003
  ! Purpose: calculated simulated reflectance and weighting function/first 
  ! derivative for all the fitting variables. Basically, call the lidort 
  ! to calculate radiance and ozone, albedo weighting functions, where 
  ! species other than ozone are not taken into account. Then calibrate 
  ! the measured radiance/solar spectra to obtain the measured reflectance 
  ! and use the finite difference to obtain the first derivative for all 
  ! the other variables.
  ! **************************************************************************
  ! Need to do more about the albedo 

SUBROUTINE pseudo_model (num_iter, refl_only, ns, nf, fitvar, fitvarap, dyda, gspec,    &
     fitres, fitspec, chisq, relrms, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : maxlay
  USE OMSAO_indices_module,   ONLY : ring_idx, &
        shift_offset, maxalb, maxoth, maxwfc,&
         bro_idx, bro2_idx, so2_idx, so2v_idx, no2_t1_idx, hcho_idx, o2o2_idx,&
         ring_idx, us1_idx, us2_idx, com_idx, com1_idx, com2_idx, com3_idx, &
         fsl_idx, rsl_idx, sdc_idx
  USE OMSAO_variables_module, ONLY : fitwavs, fitweights, sza => the_sza_atm, &
       vza => the_vza_atm, aza => the_aza_atm, sca=>the_sca_atm, &
       fitvar_rad, mask_fitvar_rad, database_indices,           &
       database_shiwf, slwf, npix_fitted, database, fitvar_rad_str, & 
       numwin, nradpix, band_selectors, refidx,  scnwrt, & 
       refspec_norm, rmask_fitvar_rad, fitvar_rad_apriori, actspec_rad, &
       n_slitvar, mask_slitvar, database_pslwf
  USE ozprof_data_module,     ONLY : nlay, use_lograd,                       &
       ozf_fidx => ozfit_start_index, ozf_lidx => ozfit_end_index,           &
       ozp_fidx=>ozprof_start_index, ozp_lidx => ozprof_end_index,           &
       stlay => start_layer, endlay => end_layer, & 
       nt_fit, t_fidx, t_lidx, tf_fidx, tf_lidx,  & 
       nwfc, wfcidx, nfwfc, wfcfidx, wfcfpix, wfclpix, &
       nalb, albidx,nfalb, albfidx,  albfpix, alblpix, &
       ngas, gasidxs, fgasidxs, tracegas, the_ai, the_cfrac, the_snowice, &
       ecfrind, ecfrfind, ecodfind, ectpfind, taodfind, twaefind, saodfind,  &
       sprsfind,so2zfind,  &
       nothgrp, osind, osfind, slind, slfind, shind, shfind, rnind, rnfind,  &
       dcfind, dcind, isind, isfind, irind, irfind, oswins, slwins, shwins,  &
       rnwins, dcwins, iswins, irwins, nos, nsl, nsh, nrn, ndc, nis, nir,   &
       p1ind, p1find, p1wins, np1, p2ind, p2find, p2wins, np2, &
       cmfind, cmind, cmwins, ncm, &
       do_subfit, radcalwrt, do_simu, do_simu_rmring, fit_atanring, & 
       use_effcrs, ncalcp, saa_flag, nsaa_spike,  vary_sfcalb, &
       pos_alb, which_cld, & 
       maxpol, npol, nfpol, polmin, polmax, polfpix, pollpix, polidx, polfidx, &
       polwf, polcorr, polcc, &
       is_albspcvar, albspcs, sfcalbs, use_albspc, use_albeofs, nactalbspc

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
  REAL (KIND=dp), INTENT(OUT), DIMENSION (ns, nf) :: dyda
  REAL (KIND=dp), INTENT(OUT)                     :: chisq, relrms

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER :: MSTKS = 4
  INTEGER :: n0alb, n0wfc, i, j, k, iw, ReturnStatus, ridx, sidx, fidx, lidx, &
       idx, albord, min_ssa_iter, swin, ewin, ig, nord, nostk, wfcord, ntmp,  &
       albsidx, albeidx , m, slit_idx, idx330, n0, ord
  INTEGER, DIMENSION(maxalb)               :: albpmax, albpmin
  INTEGER, DIMENSION(maxwfc)               :: wfcpmax, wfcpmin
  INTEGER, DIMENSION(maxpol)               :: polpmax, polpmin
  INTEGER, DIMENSION(numwin, maxoth)       :: tmpind, tmpfind
  INTEGER, DIMENSION(maxoth, 2)            :: tmpwins
  REAL (KIND=dp), DIMENSION(maxalb)        :: albarr 
  REAL (KIND=dp), DIMENSION(maxwfc)        :: wfcarr 
  REAL (KIND=dp), DIMENSION (ns)           :: delpos, waves, meas1, meas2, &
       sim1, sim2, simrad, simrad1, fitspec1, temporwf, corr
  REAL (KIND=dp), DIMENSION(ns,nlay,MSTKS) :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(ns, 4)         :: albothwf, wfcothwf, polothwf
  REAL (KIND=dp), DIMENSION(ns, MSTKS)     :: o3shiwf, cfracwf, albwf, fsimrad, &
       ctpwf, codwf, saodwf, taodwf, twaewf,  sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(numwin, maxoth):: o3shi
  REAL (KIND=dp), DIMENSION(nlay)          :: tprof, ozprof, ozadj, ozaprof
  REAL (KIND=dp), DIMENSION (nf)           :: fitvar_saved
  REAL (KIND=dp)                           :: rms, radrms, wavavg, cfrac, the_salb
  REAL (KIND=dp)                           :: newoz, newso2, newbro,   newhcho, newno2, newo4
  REAL (KIND=dp)                           :: so2adj, so2vadj, broadj, hchoadj, no2adj, o4adj
  REAL (KIND=dp), DIMENSION (numwin)       :: allrms, allchisq, allradrms

  LOGICAL :: do_ozwf, do_albwf, do_o3shi, do_tmpwf, do_shiwf, do_taodwf, do_twaewf, &
              do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf, do_so2zwf, do_pslwf, do_polwf 
  LOGICAL :: negval, so2negval, so2vnegval, hchonegval, bronegval, no2negval, o4negval
  LOGICAL, DIMENSION (nlay)     :: ozvary
  REAL (KIND=dp), DIMENSION(ns) :: walb0s, wfc0s 
  CHARACTER(LEN=1)  :: ordchar
  CHARACTER(LEN=3)  ::tmpc

  ! measurement error covariance Random + Systematic
  REAL(KIND=dp), DIMENSION(ns, ns)  :: Sy, Sy_inv
  REAL(KIND=dp), DIMENSION(ns, 1)   :: y1
  REAL(KIND=dp), DIMENSION(1, 1)    :: chi
  REAL(KIND=dp), ALLOCATABLE        :: y1tmp(:, :), Sy_invtmp(:, :)

  ! xliu, 08/10/2010
  ! Current VLIDORT calculation is based on single surface albedo (per channel). 
  ! The wavelength dependence of surface albedo on radiance is corrected through
  ! weighting function. Howver, there are not accounted for in the calculation of
  ! weighting functions. 
  ! LOGICAL :: vary_sfcalb ! global variables                

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=12), PARAMETER :: modulename = 'pseudo_model'

  ! To compute fine and radiance calculaiton wavelngth grids when using high resolution wavelength grids
  LOGICAL, SAVE             :: first = .TRUE.

  errstat = pge_errstat_ok
  IF (first .AND. .NOT. use_effcrs) THEN
     CALL get_hres_radcal_waves(errstat)
     IF (errstat == pge_errstat_error) THEN
        WRITE(*, *) modulename, ': Errors in getting fine & radiance calculation wavelength grids!!!'
        RETURN
     ENDIF
     first = .FALSE.
  ENDIF
  ! ================ Determine flags for linearization ======================
  do_albwf  = .TRUE.; do_ozwf = .TRUE.; do_o3shi = .TRUE.; do_tmpwf = .TRUE.; do_cfracwf = .TRUE.; do_pslwf = .TRUE.
  do_polwf=.TRUE.
  IF (nfalb <= 0) do_albwf = .FALSE. 
  IF (nfwfc <= 0) do_cfracwf = .FALSE.
  IF (refl_only .AND. .NOT. use_effcrs) THEN
     do_ozwf = .FALSE.; do_o3shi = .FALSE.; do_tmpwf = .FALSE.
  END IF
  IF (nos <= 0)    do_o3shi = .FALSE.
  IF (np1+np2 <= 0)    do_pslwf  = .FALSE. ; IF (nfpol <=0 ) do_polwf = .FALSE.
  IF (nt_fit <= 0) do_tmpwf = .FALSE. 
  ozvary = .FALSE.; ozvary(stlay:endlay) = .TRUE.
  do_ctpwf  = .FALSE.;  do_codwf  = .FALSE.
  do_taodwf  = .FALSE.; do_twaewf = .FALSE.;  do_saodwf = .FALSE.
  do_sprswf = .FALSE. ; do_so2zwf = .FALSE. 
  IF (.NOT. refl_only) THEN
     IF (ecfrfind > 0) do_cfracwf = .TRUE.
     IF (ecodfind > 0) do_codwf   = .TRUE.
     IF (ectpfind > 0) do_codwf   = .TRUE.
     IF (taodfind > 0) do_taodwf  = .TRUE.
     IF (saodfind > 0) do_saodwf  = .TRUE.
     IF (twaefind > 0) do_twaewf  = .TRUE.
     IF (sprsfind > 0) do_sprswf  = .TRUE.
     IF (so2zfind > 0) do_so2zwf  = .TRUE.
  ENDIF

  ! Update cloud fraction or disable fitting cloud fraction when fc<=1.0E-3 or fc>=0.999
  IF ( ecfrfind > 0 ) THEN
     the_cfrac = fitvar_rad(ecfrind)
     IF (the_cfrac <= 0.2) THEN
        do_cfracwf = .FALSE.
     ENDIF
     IF (the_cfrac <= 1.0E-2 ) THEN
        the_cfrac = 0.0D0;  fitvar(ecfrfind) = 0.0D0; fitvar_rad(ecfrind) = 0.0D0
        fitvarap(ecfrfind) = 0.0D0;    fitvar_rad_apriori(ecfrind) = 0.0D0
        do_cfracwf = .FALSE.      
     ELSE IF (the_cfrac >= 0.99 ) THEN
        the_cfrac = 1.0D0;  fitvar(ecfrfind) = 1.0D0; fitvar_rad(ecfrind) = 1.0D0
        fitvarap(ecfrfind) = 1.0D0;    fitvar_rad_apriori(ecfrind) = 1.0D0
        do_cfracwf = .FALSE. 

        ! When using Lambertian cloud model, surface albedo is already fitted, so nothing needs to be done
        ! But when using scattering cloud model, need to fit an optical thickness and a surface albedo
        ! To be implemented !!!
     ENDIF
  ENDIF
  
  ! Disable fitting surface albedo when 0.99>=fc>=0.01
  ! For safe, as we have set sa(i, i) = 0.0 in specfit_ozprof 
  !IF ( do_alb_longwav .AND. nfwfc > 0) THEN
  IF (nfwfc > 0) THEN
     DO i = wfcidx + nwfc - 1, wfcidx, -1
        IF (fitvar_rad_str(i)(4:4) == '0') EXIT
     ENDDO
     the_cfrac = fitvar_rad(i)   ! Use UV2/last channel cloud fraction
     DO i = albidx + nalb - 1, albidx, -1
        IF (fitvar_rad_str(i)(4:4) == '0') EXIT
     ENDDO
     the_salb = fitvar_rad(i)   ! Use UV2/last channel cloud fraction

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
      ELSE IF (the_cfrac >=0.8 .and. (the_snowice <1 .and. the_snowice >= 104) ) THEN
           do_albwf = .false.
      ELSE IF ( the_cfrac > 0.2 .and. do_cfracwf == .true.) THEN  
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
  tmpc = 'F,F'
  if ( do_albwf)    tmpc(3:3) = 'T'
  if ( do_cfracwf)  tmpc(1:1) = 'T'
 ! WRITE(*,'(i3, a3, 2f8.2, 2e15.7)') nfwfc, tmpc, the_cfrac,the_salb,    fitvar(wfcfidx+nfwfc-1), fitvar(albfidx+nfalb-1)
!  print * , fitvar(wfcfidx:wfcfidx+nfwfc-1)
!  print * , fitvar(albfidx:albfidx+nfalb-1)
  ! ======= Set up ozone, temperature, trace gases, albedo, lamda for LIDORT ============
  ozprof(1:nlay)  = fitvar_rad (ozp_fidx:ozp_lidx)
  ozaprof(1:nlay) = fitvar_rad_apriori(ozp_fidx:ozp_lidx)
  tprof(1:nlay)   = fitvar_rad(t_fidx:t_lidx)

  !WRITE(*, '(A)') 'Initial Guess (pseudo_model-1): '
  !WRITE(*, '(12F8.3)') fitvar_rad(ozp_fidx:ozp_lidx), SUM(fitvar_rad(ozp_fidx:ozp_lidx))
  !WRITE(*, '(A)') 'Initial Guess (pseudo_model-2): '
  !WRITE(*, '(12F8.3)') fitvar(ozf_fidx:ozf_lidx), SUM(fitvar(ozf_fidx:ozf_lidx))

  !xliu (02/01/2007): adjust ozone profile for negative ozone values
  !                   radiances will be corrected using ozone weighting function
  ozadj(1:nlay) = 0.0; negval = .FALSE.
  DO i = 1, nlay
     IF (ozprof(i) <= 0.0) THEN
        newoz  = MIN(0.5, ozaprof(i))
        negval = .TRUE. ; ozadj(i)  = newoz - ozprof(i); ozprof(i) = newoz
        do_ozwf = .TRUE.; ozvary(i) = .TRUE.
     ENDIF
  ENDDO

  so2negval = .FALSE.; so2vnegval = .FALSE.
  hchonegval = .FALSE.; bronegval = .FALSE.
  no2negval = .FALSE. ; o4negval = .FALSE.

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
           newo4 = 1.0E12 * refspec_norm(gasidxs(k))
           o4adj = newo4 - tracegas(k, 4)
           tracegas(k, 4) = newo4
        ENDIF
     ENDIF
  ENDDO

  waves = fitwavs(1:ns) 
  walb0s = 0.0
  n0alb =  0
  DO i = 1, nalb
     j = albidx - 1 + i
     !xliu, 02/08/2012, add albord and **albord 
     READ(fitvar_rad_str(j)(4:4), '(I1)') albord
     fidx=albfpix(i); lidx=alblpix(i)
     wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
     IF (albord == 0) THEN
        n0alb = n0alb + 1
        albarr(n0alb)  = fitvar_rad(j)
        albpmax(n0alb) = alblpix(i); albpmin(n0alb) = albfpix(i)
     ENDIF
     IF (vary_sfcalb) THEN 
        walb0s(fidx:lidx) = walb0s(fidx:lidx) + fitvar_rad(j) * (waves(fidx:lidx) - wavavg)**albord
     ENDIF
  ENDDO

  n0wfc = 0
  wfc0s (:) = 0.0
  DO i = 1, nwfc
     j = wfcidx - 1 + i
     READ(fitvar_rad_str(j)(4:4), '(I1)') wfcord
     fidx = wfcfpix(i); lidx = wfclpix(i)   
     wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
     IF (wfcord == 0) THEN
        n0wfc = n0wfc + 1
        wfcarr(n0wfc) = fitvar_rad(j)
        wfcpmax(n0wfc) = wfclpix(i); wfcpmin(n0wfc) = wfcfpix(i)
     ENDIF
     IF (vary_sfcalb) THEN
         wfc0s(fidx:lidx) =  wfc0s(fidx:lidx) + fitvar_rad(j) * (waves(fidx:lidx) - wavavg)**wfcord
     ENDIF
  ENDDO

  n0 = 0
  polcc(:) = 1.0
  IF (do_polwf) THEN 
  DO i = 1, npol
      j = polidx - 1 + i
      READ (fitvar_rad_str(j)(4:4), '(I1)')  ord
      IF (ord == 0) THEN
          n0 = n0 + 1
          polpmax(n0) = pollpix(i); polpmin(n0) = polfpix(i)
          polcc(polpmin(n0):polpmax(n0)) = fitvar_rad(j)
      ELSE
          fidx = polfpix(i) ; lidx = pollpix(i)
          wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
          polcc(fidx:lidx) = polcc(fidx:lidx) + fitvar_rad(j)*(waves(fidx:lidx) - wavavg)**ord
      ENDIF
      print * , i,nfpol, fitvar_rad(j), polcorr
  ENDDO
  ENDIF
 
  IF (do_subfit) THEN
     DO i = 1, maxoth
        o3shi(1:numwin, i) = fitvar_rad(osind(1:numwin, i))
     ENDDO
  ELSE
     o3shi(1, 1:maxoth)    = fitvar_rad(osind(1, 1:maxoth))
  ENDIF

  waves = fitwavs(1:ns)

  !Iraddaince/radiance shift is done when interpolating solar reference to wavelength grid, not here
  !IF (nsh > 0) THEN
  !   IF (do_subfit) THEN
  !      fidx = 1
  !      DO j = 1, numwin
  !         lidx = fidx + nradpix(j) - 1
  !         delpos(fidx:lidx) =  waves(fidx:lidx) - (waves(fidx) + waves(lidx)) / 2.0
  !         IF (shfind(j, 1) > 0)  waves(fidx:lidx) = waves(fidx:lidx) + fitvar_rad(shind(j, 1)) 
  !
  !         DO i = 2, nsh            
  !            IF (shfind(j, i) > 0)  waves(fidx:lidx) = waves(fidx:lidx) + &
  !                 fitvar_rad(shind(j, i)) * (delpos(fidx:lidx) ** (i-1))
  !         ENDDO
  !         fidx = lidx + 1
  !      ENDDO
  !   ELSE
  !      IF (shwins(1, 1) == 1) THEN
  !         fidx = 1
  !      ELSE
  !         fidx = SUM(nradpix(1: shwins(1, 1)-1)) + 1
  !      ENDIF
  !      lidx = SUM(nradpix(1: shwins(1, 2)))
  !      delpos(fidx:lidx) =  waves(fidx:lidx) - (waves(fidx) + waves(lidx)) / 2.0
  !      IF (shfind(1, 1) > 0) waves(fidx:lidx)  = waves(fidx:lidx) + fitvar_rad(shind(1, 1))
  !
  !      DO i = 2, nsh  
  !         IF (shfind(1, i) > 0) THEN
  !            waves(fidx:lidx)  = waves(fidx:lidx) + fitvar_rad(shind(1, i)) * (delpos(fidx:lidx) ** (i-1))
  !         ENDIF
  !      ENDDO
  !   ENDIF
  !ENDIF
  ! === Call LIDORT, polarization correction, and additional wf calc =====
  !the_cfrac = 0.
  nostk = 1
  IF (use_effcrs) THEN
     CALL LIDORT_PROF_ENV(do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary,   &
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
          do_sprswf, do_so2zwf, do_pslwf,  ns, waves, maxoth, o3shi, sza, vza, aza, &
          nlay, ozprof, tprof, n0alb, albarr, albpmin,albpmax, vary_sfcalb, walb0s, &
          n0wfc, wfcarr, wfcpmin, wfcpmax, wfc0s, &
          nostk, albwf(1:ns, 1:nostk), ozwf(1:ns, 1:nlay, 1:nostk), &
          tmpwf(1:ns, 1:nlay, 1:nostk),o3shiwf(1:ns, 1:nostk), &
          cfracwf(1:ns, 1:nostk), codwf(1:ns, 1:nostk), ctpwf(1:ns, 1:nostk), &
          taodwf(1:ns, 1:nostk), twaewf(1:ns, 1:nostk), saodwf(1:ns, 1:nostk),&
          sprswf(1:ns, 1:nostk), so2zwf(1:ns, 1:nostk), fsimrad(1:ns, 1:nostk), errstat)
     IF (errstat == pge_errstat_error) &
          WRITE(*, *) modulename, ': Errors in calling LIDORT_PROF_ENV!!!'
  ELSE
     ntmp = MAX(ncalcp, ns)
     CALL HRES_RADCALC_ENV(ntmp, do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary,            &
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
          do_sprswf, do_so2zwf, do_pslwf, & 
          ns, waves, maxoth, o3shi, sza, vza, aza, nlay, ozprof, tprof, n0alb,  &
          albarr, albpmin, albpmax, vary_sfcalb, walb0s, & 
          n0wfc, wfcarr, wfcpmin, wfcpmax, wfc0s, &
          nostk, albwf(1:ns, 1:nostk), ozwf(1:ns, 1:nlay, 1:nostk), &
          tmpwf(1:ns, 1:nlay, 1:nostk), o3shiwf(1:ns, 1:nostk), &
          cfracwf(1:ns, 1:nostk), codwf(1:ns, 1:nostk), ctpwf(1:ns, 1:nostk), &
          taodwf(1:ns, 1:nostk), twaewf(1:ns, 1:nostk), saodwf(1:ns, 1:nostk), &
          sprswf(1:ns, 1:nostk), so2zwf(1:ns, 1:nostk), fsimrad(1:ns, 1:nostk), errstat)
     IF (errstat == pge_errstat_error) &
          WRITE(*, *) modulename, ': Errors in calling HRES_RADCALC_ENV!!!'
  ENDIF
  IF (errstat == pge_errstat_error) RETURN
  !xliu (02/01/2007): correct radiances based on ozone weighting function to deal with negative ozone values
  IF (negval) THEN
     DO i = 1, nlay 
        IF (ozadj(i) > 0) THEN
           fsimrad(1:ns, 1:nostk) = fsimrad(1:ns, 1:nostk) - ozadj(i) * ozwf(1:ns, i, 1:nostk) 
        ENDIF
     ENDDO
  ENDIF
  IF (do_albwf == .false.) albwf(:,:) = 0.0D0
  IF (do_cfracwf == .false.) cfracwf(:,:) = 0.0D0
  !!xliu (12/11/2014): correct radiances based on SO2 to deal with negative
  !ozone values
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
        ENDIF

     ENDIF
  ENDDO

  simrad = fsimrad(1:ns, 1)

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

  !PRINT *, polcorr, ns, nlay
  !print *, sza, aza, vza
  !DO i = 1, ns
  ! WRITE(77, '(f12.5, 43d14.6)') fitwavs(i), simrad(i), ozwf(i, 1:nlay, 1), albwf(i, 1) &!, taodwf(i, 1) !&
  !        , o3shiwf(i, 1) !* ozprof(1:nlay)
  !ENDDO
  !CLOSE(77)
  !STOP


  ! get dlnI/dx
  IF (use_lograd) THEN
     IF (do_albwf)   albwf(1:ns, 1)   = albwf(1:ns, 1)   / simrad     
     IF (do_o3shi)   o3shiwf(1:ns, 1) = o3shiwf(1:ns, 1) / simrad
     IF (do_codwf)   codwf(1:ns, 1)   = codwf(1:ns, 1)   / simrad
     IF (do_sprswf)  sprswf(1:ns, 1)  = sprswf(1:ns, 1)  / simrad
     IF (do_so2zwf)  so2zwf(1:ns, 1)  = so2zwf(1:ns, 1)  / simrad
     IF (do_ctpwf)   ctpwf(1:ns, 1)   = ctpwf(1:ns, 1)   / simrad
     IF (do_cfracwf) cfracwf(1:ns, 1) = cfracwf(1:ns, 1) / simrad
     IF (do_taodwf)  taodwf(1:ns, 1)  = taodwf(1:ns, 1)  / simrad
     IF (do_saodwf)  saodwf(1:ns, 1)  = saodwf(1:ns, 1)  / simrad
     IF (do_twaewf)  twaewf(1:ns, 1)  = twaewf(1:ns, 1)  / simrad
     IF (do_polwf)   polwf(1:ns, 1)   = polwf(1:ns, 1)   / simrad

     IF (do_ozwf) THEN
        DO i = stlay, endlay
           ozwf (:, i, 1) = ozwf(:, i, 1) / simrad
           IF (nt_fit > 0) tmpwf(:, i, 1) = tmpwf(:, i, 1) / simrad
        END DO
     ENDIF

     simrad = LOG(simrad)           ! get dlnI     
  END IF
  !WRITE(77, *) 'Ozone weighting function D(lnI)/D(lnx)'
  !DO i = 1, ns
  !  WRITE(77, '(f8.4, 43d14.6)') fitwavs(i), ozwf(i, 1:nlay) !* ozprof(1:nlay)
  !ENDDO

  IF (do_albwf   == .false.) albwf(:,:) = 0.0D0
  IF (do_cfracwf == .false.) cfracwf(:,:) = 0.0D0
  IF (do_polwf   == .false.) polwf(:,:) = 0.0D0
  ! correct for linear/quardratic wavelength dependent in albedo
  albothwf = 0.0
  DO i = 1, nalb
     j = albidx + i - 1

     READ(fitvar_rad_str(j)(4:4), '(I1)') albord
     IF (albord == 0) CYCLE

     fidx=albfpix(i); lidx=alblpix(i)
     wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
     !albothwf(fidx:lidx,albord) = albwf(fidx:lidx, 1)*(waves(fidx:lidx) / wavavg)**albord
     albothwf(fidx:lidx,albord) = albwf(fidx:lidx, 1)*(waves(fidx:lidx) - wavavg)**albord ! much better than above
     IF (.NOT. vary_sfcalb) simrad(fidx:lidx) = simrad(fidx:lidx) +  albothwf(fidx:lidx, albord) * fitvar_rad(j)        
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
        IF (.NOT. vary_sfcalb) simrad(fidx:lidx) = simrad(fidx:lidx) +  wfcothwf(fidx:lidx, wfcord) * fitvar_rad(j)        
     ENDDO
  ENDIF

  polothwf = 0.0
  IF (do_polwf) THEN 
  DO i = 1, npol
     j = polidx + i - 1
     READ (fitvar_rad_str(j)(4:4), '(I1)') ord
     IF (ord == 0) CYCLE
     fidx = polfpix(i); lidx = pollpix(i)
     wavavg = SUM(waves(fidx:lidx)/(1.0 + lidx - fidx))
     polothwf(fidx:lidx, ord) = polwf(fidx:lidx, 1) * (waves(fidx:lidx) -wavavg)**ord
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

  IF (ncm > 0) THEN 
   print * , 'please update common mode spectra for new instrumnet'
   !call get_sdc_spec (ns, fitwavs, corr)
   !      database (sdc_idx, refidx(1:ns)) = corr(1:ns) !*exp(simrad(1:ns))
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
  ! compute chi-square difference
 
     chisq  = SUM((fitres / fitweights(1:ns))**2.0)
 

  rms =    SQRT(chisq  / REAL(ns, KIND=dp))
  relrms = 100.D0 * SQRT(SUM(ABS((simrad-fitspec) / fitspec)**2.0) &
       / REAL(ns, KIND=dp))
  
 
  IF (scnwrt) THEN
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1

           allchisq(i) = SUM((fitres(fidx:lidx) / fitweights(fidx:lidx))**2.0)

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
    DO i = 1, numwin
        WRITE (*, '(A13, I2, 2(A14, 1pd11.3))') 'Win ', i, ': allrms = ', &
             allrms(i), 'allIrms(%) = ', allradrms(i)
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
        fidx = albfpix(k -albidx + 1); lidx=alblpix(k - albidx + 1)
           
        READ(fitvar_rad_str(k)(4:4), '(I1)') albord
        IF (albord == 0) THEN
           !IF ( do_alb_longwav .AND. fitvar_rad_str(k) == '2ba0' ) THEN
           !   dyda(albsidx:lidx, j) = albwf(albsidx:lidx, 1)
           !ELSE 
           dyda(fidx:lidx, j)=albwf(fidx:lidx, 1)
           !ENDIF
           !dyda(albsidx:lidx, j) = albwf(albsidx:lidx, 1)
        ELSE

              dyda(fidx:lidx, j)=albothwf(fidx:lidx, albord)      
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

        IF (wfcord == 0) THEN
           !IF ( do_alb_longwav .AND. fitvar_rad_str(k) == '2fc0' ) THEN
           !   dyda(fidx:lidx, j) = cfracwf(albsidx:lidx, 1)
           !ELSE
           dyda(fidx:lidx, j) = cfracwf(fidx:lidx, 1)
           !ENDIF
        ELSE 
      
           dyda(fidx:lidx, j) = wfcothwf(fidx:lidx, wfcord)

        ENDIF
     ENDDO

     DO i = 1, nfpol
        j = polfidx + i -1
        k = mask_fitvar_rad(j)
        fidx = polfpix(k - polidx + 1); lidx = pollpix(k - polidx + 1)
        READ (fitvar_rad_str(k)(4:4), '(I1)') ord
        IF (ord == 0) THEN
           !print  * , i, j, fidx, lidx, polidx
           dyda(fidx:lidx, j) = polwf(fidx:lidx, 1)
        ELSE
           dyda(fidx:lidx, j) = polothwf(fidx:lidx, ord)
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
        IF (ridx /= us1_idx .AND. ridx /= us2_idx .AND. &
            ridx /= com_idx .AND. ridx /= com1_idx .AND. ridx /= com2_idx .AND. ridx /= com3_idx .AND. &
            ridx /= fsl_idx .AND. ridx /= rsl_idx) THEN
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
           fitvar(i) = fitvar_saved(i) * 0.999
           sim2 = simrad
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
           nord =0
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
        ELSE IF (ig == 8 .or. ig == 9) THEN 
           IF (ig  == 8 ) THEN                                                         
             nord = np1; tmpind = p1ind; tmpfind = p1find; tmpwins = p1wins ; slit_idx=mask_slitvar(1)
           ELSE IF (ig == 9 ) THEN 
             nord = np2; tmpind = p2ind; tmpfind = p2find; tmpwins = p2wins ; slit_idx=mask_slitvar(2)
           ENDIF
           temporwf(1:ns) = - database_pslwf(slit_idx, 1:ns)
           !IF (use_lograd) temporwf(1:ns) = temporwf /simrad           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
        ELSE IF (ig == 10 ) THEN 
           nord = ncm; tmpind = cmind; tmpfind = cmfind; tmpwins = cmwins
           temporwf(1:ns) = database(sdc_idx, refidx(1:ns)) 
           !if (num_iter > 0 ) print *, num_iter, fitvar_rad(cmind(1,1)),fitvar_rad(cmind(2,1))
           !print * , temporwf(ns-10:ns)
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
           !if (ig ==10  ) print *, fitvar(tmpfind(1:2, 1))
        ENDIF
        ! Use finite differences to get the zero-order weighting functions
        IF (nord > 0 .AND. ig >= 3 .AND. ig <= nothgrp .AND. ig /= 4 .AND. ( ig < 8 .or. ig > 10) ) THEN
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
              temporwf = -(meas1 - sim1 - meas2 + sim2) / 0.002
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
              if (ig == 10) THEN 
              !     print * , meas1(ns-10:ns) - meas2(ns-10:ns)
              ENDIF
              IF (swin == ewin) THEN
                 IF (fitvar_saved(tmpfind(swin, 1)) /= 0.0) THEN
                    temporwf = -(meas1 - sim1 - meas2 + sim2) / (0.002 * fitvar_saved(tmpfind(swin, 1)))
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

     !WRITE(92, *) ns, nf
     !DO i = 1, ns
     !   WRITE(92, '(f10.4, 80d14.6)') fitwavs(i), dyda(i, 1:nf)
     !ENDDO
     !errstat = pge_errstat_error; RETURN
     !
     !CLOSE(92)
     !STOP

     ! Remove longer wavelengths wfs except for alb and cloud
     !idx330 = maxval(minloc(waves(1:ns), mask=(waves(1:ns)>=330.0)))
     !IF (idx330  > 0 .AND. idx330  <= ns) dyda(idx330:ns, ozf_fidx:ozf_lidx) = 0.0D0
     !IF ( do_alb_longwav) then 
      !do i = 1, nf
        !if ( (i < albfidx .or. i > albfidx + nfalb - 1) .and. i .ne. ecfrfind ) then  
        !if ( (i < albfidx .or. i > albfidx + nfalb - 1) )then 
        !if ( (i < albfidx +nfalb-1 .or. i > albfidx + nfalb - 1) )then 
        !   dyda(idx330:ns,i) = 0.0d0 
        !endif
      !enddo 
      ! do not use longer waves for ozone, this make no sense 
     ! IF (idx330  > 0 .AND. idx330  <= ns) dyda(idx330:ns, ozf_fidx:ozf_lidx) = 0.0D0
     !ENDIF

     DO i = 1, nf  
        dyda(:, i) = dyda(:, i) / fitweights(1:ns)
     !print * , i, mask_fitvar_rad(i),  fitvar_rad_str(mask_fitvar_rad(i)), dyda(1:10, i)
     END DO
     ! finnally obtain the new spectrum to be fitted in the GSVD
     gspec(1:ns) = fitres(1:ns) / fitweights(1:ns)
     ! Restore the unperturbated fitting variables
     fitvar(1:nf) = fitvar_saved(1:nf)
     fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar(1:nf)  
  END IF
  RETURN

END SUBROUTINE pseudo_model

SUBROUTINE HRES_RADCALC_ENV (nw0, do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary,    &
     do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf,    &
     do_so2zwf, do_pslwf, nw, waves, nos, o3shi, sza, vza, aza, nl, ozprof, tprof, n0alb,     &
     albarr, albpmin, albpmax, vary_sfcalb, walb0s, n0wfc, wfcarr, wfcpmin, wfcpmax, wfc0s, &
     nostk, albwf, ozwf, tmpwf, o3shiwf, cfracwf, codwf, ctpwf, taodwf, twaewf,     &
     saodwf, sprswf, so2zwf, rad, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : fitvar_rad_str, fitwavs, numwin,  &
       fitvar_rad, mask_fitvar_rad, rmask_fitvar_rad
  USE ozprof_data_module,     ONLY : use_effcrs, radcwav, ncalcp,      &
       albfidx, nalb, nfalb, albidx, albmin, albmax, albfpix, alblpix, &
       wfcfidx, nwfc, nfwfc, wfcidx, wfcmin, wfcmax, wfcfpix, wfclpix, &
       hreswav, which_cld
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN) :: nw0, nw, nl, nos, n0alb, nostk, n0wfc
  LOGICAL, INTENT(IN) :: do_ozwf, do_albwf, do_tmpwf, do_o3shi, do_taodwf, vary_sfcalb, &
       do_twaewf, do_saodwf, do_cfracwf, do_codwf, do_ctpwf, do_sprswf, do_so2zwf, do_pslwf
  INTEGER, INTENT(OUT)                                 :: errstat
  INTEGER, DIMENSION(n0alb), INTENT(IN)                :: albpmax, albpmin
  INTEGER, DIMENSION(n0wfc), INTENT(IN)                :: wfcpmax, wfcpmin
  LOGICAL, DIMENSION(nl), INTENT(IN)                   :: ozvary
  REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)           :: waves, walb0s, wfc0s
  REAL (KIND=dp), DIMENSION(nw, nostk),  INTENT(OUT)   :: rad, albwf, cfracwf, o3shiwf, &
       codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)   :: o3shi
  REAL (KIND=dp), DIMENSION(nl),  INTENT(IN)           :: ozprof, tprof
  REAL (KIND=dp), DIMENSION(nw, nl, nostk),INTENT(OUT) :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(n0alb), INTENT(IN)         :: albarr 
  REAL (KIND=dp), DIMENSION(n0wfc), INTENT(IN)         :: wfcarr 
  REAL (KIND=dp), INTENT(IN)                           :: sza, vza, aza

  ! =======================
  ! Local variables
  ! =======================
  INTEGER :: i, j, k, m, n, fidx, lidx, fidx0, lidx0, albord, wfcord
  REAL (KIND=dp), DIMENSION(nw0)            :: waves0, walb0s0, wfc0s0
  REAL (KIND=dp), DIMENSION(nw0, nostk)     :: rad0, albwf0, cfracwf0, o3shiwf0, &
       codwf0, ctpwf0, taodwf0, twaewf0, saodwf0, sprswf0, so2zwf0
  REAL (KIND=dp), DIMENSION(nw0, nl, nostk) :: ozwf0, tmpwf0
  REAL (KIND=dp)                            :: wavavg
  INTEGER, DIMENSION(n0alb)                 :: albpmax0, albpmin0
  INTEGER, DIMENSION(n0wfc)                 :: wfcpmax0, wfcpmin0
  CHARACTER(LEN=1)  :: ordchar
  errstat = pge_errstat_ok

  ! Note nw0 is MAX(ncalcp, nw)
  waves0 = 0.0
  waves0(1:ncalcp) = radcwav(1:ncalcp)

  ! Need to find indices of boundaries for using different surface albedos/cloud fractions
   
  k = 0; albpmin0(1) = 1; albpmax0(n0alb) = ncalcp 
  DO i = 1, nalb
     j = albidx - 1 + i

     READ(fitvar_rad_str(j)(4:4), '(I1)') albord
     IF (albord == 0) THEN
        k = k + 1

        IF (k > 1) albpmin0(k)= albpmax0(k-1) + 1 
        IF (k < n0alb) albpmax0(k)= MINVAL(MAXLOC(waves0(1:ncalcp), MASK=(waves0(1:ncalcp) &
             >= albmin(i) .AND. waves0(1:ncalcp) < albmax(i))))  

        fidx0 = albpmin0(k); lidx0 = albpmax0(k)
        IF (vary_sfcalb) walb0s0(albpmin0(k):albpmax0(k)) = albarr(k)

  
     ELSE
        IF (vary_sfcalb) THEN
           ! Note use the exact average wavelength as retrieval grid
           fidx=albfpix(i); lidx=alblpix(i)
           wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx)) 

           ! Get surface albedo for each radiance calculation wavelength
           walb0s0(fidx0:lidx0) = walb0s0(fidx0:lidx0) + fitvar_rad(j) * (waves0(fidx0:lidx0) - wavavg)**albord


        ENDIF
     ENDIF
  ENDDO

  IF (n0wfc > 0) THEN
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
           IF (vary_sfcalb) wfc0s0(wfcpmin0(k):wfcpmax0(k)) = wfcarr(k)
      
        ELSE
           IF (vary_sfcalb) THEN
              fidx=wfcfpix(i); lidx=wfclpix(i)
              wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx)) 
              wfc0s0(fidx0:lidx0) = wfc0s0(fidx0:lidx0) + fitvar_rad(j) * (waves0(fidx0:lidx0) - wavavg)**wfcord
         
           ENDIF
        ENDIF
     ENDDO
  ENDIF
  !do i = 1, ncalcp
  !   write(89, '(F9.4, 1p30D16.7)') waves0(i), walb0s0(i), wfc0s0(i)
  !enddo
  !pause
  !print *, n0alb, albpmin(1:n0alb), albpmax(1:n0alb)
  !print *, waves(albpmin(1:n0alb)), waves(albpmax(1:n0alb))
  !print *, nw0, ncalcp, nw
  !print *, albpmin0(1:n0alb), albpmax0(1:n0alb)
  !print *, waves0(albpmin0(1:n0alb)), waves0(albpmax0(1:n0alb))
  !STOP

  ! Call LIDORT_PROF_ENV on fine wavelength grids
  ! Return raidances and weighting functions on required resolution wavelength grid
  CALL LIDORT_PROF_ENV (do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary, do_taodwf,      &
       do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf, do_so2zwf, do_pslwf,     &
       nw0, waves0, nos, o3shi, sza, vza, aza, nl, ozprof, tprof, n0alb, albarr,       &
       albpmin0, albpmax0, vary_sfcalb, walb0s0, n0wfc, wfcarr, wfcpmin0, wfcpmax0,    &
       wfc0s0, nostk, albwf0, ozwf0, tmpwf0, o3shiwf0, cfracwf0, codwf0, ctpwf0, taodwf0,     &
       twaewf0, saodwf0, sprswf0, so2zwf0, rad0, errstat)
  rad(1:nw, 1:nostk)         = rad0(1:nw, 1:nostk)
  albwf(1:nw, 1:nostk)       = albwf0(1:nw, 1:nostk)
  o3shiwf(1:nw, 1:nostk)     = o3shiwf0(1:nw, 1:nostk)
  cfracwf(1:nw, 1:nostk)     = cfracwf0(1:nw, 1:nostk)
  codwf(1:nw, 1:nostk)       = codwf0(1:nw, 1:nostk)
  ctpwf(1:nw, 1:nostk)       = ctpwf0(1:nw, 1:nostk)
  taodwf(1:nw, 1:nostk)      = taodwf0(1:nw, 1:nostk)
  twaewf(1:nw, 1:nostk)      = twaewf0(1:nw, 1:nostk)
  saodwf(1:nw, 1:nostk)      = saodwf0(1:nw, 1:nostk)
  sprswf(1:nw, 1:nostk)      = sprswf0(1:nw, 1:nostk)
  so2zwf(1:nw, 1:nostk)      = so2zwf0(1:nw, 1:nostk)
  ozwf(1:nw, 1:nl, 1:nostk)  = ozwf0(1:nw, 1:nl, 1:nostk)
  tmpwf(1:nw, 1:nl, 1:nostk) = tmpwf0(1:nw, 1:nl, 1:nostk)
  RETURN
END SUBROUTINE HRES_RADCALC_ENV
end module m_pseudo_model
