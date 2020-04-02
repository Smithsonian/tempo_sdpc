MODULE m_lidort_env_vv2p7_pca
  USE GEMSTOOL_PARS_m
  USE vlidort_data_module, ONLY: &
      lidort_read_config, allocate_od, deallocate_od, &
      n_totalatmos_wfs, n_surface_wfs, &
      Inputs, Geophys, L_Geophys, &
      layer_vary_flag_cc, layer_vary_number_cc, profilewf_names_cc, &
      deltau_vert_input, deltaug_vert_input, omega_total_input,&
      fr,fa,fc, &
      l_deltau_vert_input, l_omega_total_input, & 
      l_fr, l_fa, l_fc, &
      which_win, RO, RI
  USE VLIDORT_eofpc_module, ONLY: vlidort_eofpc_master, & 
      n_call_vlidort, Assigned_bins
  USE VLIDORT_exact_module, ONLY: vlidort_exact_master
  USE VLIDORT_PARS
  USE m_ezspline_interpolation, ONLY: bspline

  IMPLICIT NONE
  INTEGER, SAVE :: nstreams, nstokes, nlayers !, n_totalatmos_wfs, n_surface_wfs
! **********************************************************
! JBAK 
! - vlidort source code : vv2p7  
!
 CONTAINS

SUBROUTINE LIDORT_PROF_ENV_PCA (do_ozwf, do_albwf, do_tmpwf, do_o3shi, &
     do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
     do_sprswf,  do_so2zwf,do_pslwf, ozvary,&
     nw, waves, nalbwf, albs, wfcs,nos, o3shi, nl, ozprof, tprof, nostk, rto, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY  : maxlay, du2mol, rearth
  USE OMSAO_indices_module,   ONLY  : so2_idx, so2v_idx, o2o2_idx, &
                                      hcho_idx, o2_idx, h2o_idx
  USE OMSAO_variables_module, ONLY  : scnwrt, currloop, &
       sza=>the_sza_atm, vza=>the_vza_atm, aza=>the_aza_atm, lat=>the_lat,&
       numwin, winlim,&
       n_rad_wvl, fitwavs, radwvl_sav, n_radwvl_sav, &
       rmask_fitvar_rad, mask_fitvar_rad, &
       fitvar_rad, fitvar_rad_init, fitvar_rad_apriori, &
       refidx, refspec_norm, database_pslwf, ctrdbdir
  USE ozprof_data_module,     ONLY : num_iter,&
       do_simu, radcalwrt, do_tracewf, use_effcrs, do_multi_vza,&
       ncalcp, radcwav, osfind, oswins,&
       mpolcorr, npolcorr,polcorr, polcorr_idxs, the_str, do_radinter, &
       mflay, nfsfc, nsfc, ntp, nt_fit, nflay,nup2p, &
       atmosprof, nup2p, fts, fps, fzs, fozs, frhos,&
       nallgas, ngas, mgasprof, tracegas, fgasidxs, fgaspos, gasidxs, nhgas, hgaspos, &
       do_lambcld, lambcld_refl, has_clouds, aerosol, strat_aerosol,&
       nmom, aerwavs,actawin,useasy, maxgksec, ngksec,&
       the_cfrac, the_cbeta, ncbp, nctp, nwfc, &
       gaext, gasca, gaasy, gamoms, gcq, gcw, gcasy, gcmoms, &
       tropsca, tropaod, tropwaer, strataod, stratsca, &
       taodind, saodind, saodfind, twaeind, sprsind, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       so2idx, so2vidx, o4idx, o2idx, h2oidx, broidx, hchoidx, no2idx, &
       ccrs, hcrs, dadp, dadt, dads, &
       rtm_outputs, set_rtmvar, npca, npcapix, which_pcabin
  USE OMSAO_errstat_module
  USE m_get_xcrs, ONLY: get_alb_ozcrs_ray , &
                        get_hres_gascrs_ray, get_effres_gascrs_ray
  USE m_lidort_util, ONLY: set_polcorr, polcorr_online, polcorr_online_with_lut,&
                           radwf_inter_convol,hres_radwf_inter_convol, &
                           get_tracegas_wf, debug_rtm, debug_taug

  IMPLICIT NONE
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
  LOGICAL :: problems, do_clouds, do_prep, do_fozwf, do_faerwf, do_fraywf, do_abs
  INTEGER :: ic, iw, i, j, k, ii, kk, jj, jk, low, hgh, fidx, lidx, &
             istk, ipol, nsprs, nz1, nfgas, nw0, &
             iwf,  zz, ip,  w1, w2, nalb, ipca, woff, nw_pca
  INTEGER :: ozwfidx, aodwfidx, twaewfidx, codwfidx, sprswfidx, raywfidx
  INTEGER, DIMENSION (nallgas)                 :: gasin   
  INTEGER, DIMENSION (5)                       :: tmp_gasidxs
  LOGICAL, DIMENSION(nflay)                    :: cldmsk, varyprof!, aermsk
  REAL (KIND=dp)                               :: lambertian_albedo, lamda, xg, frac, toz, temp, aodscl, waerscl
  REAL (KIND=dp), DIMENSION(0:nflay)           :: ozs, delps
  REAL (KIND=dp), DIMENSION(nflay)             :: cldsca, cldext, cldasy, eta, cldext0, aersca, aerext, aerasy
  REAL (KIND=dp), DIMENSION(0:nmom,maxgksec,nflay) :: cldmoms, aermoms
  REAL (KIND=dp), DIMENSION(nw)                :: tmpalbs, wfrac
  REAL (KIND=dp), DIMENSION(nallgas, nflay)    :: absod
! rad/wf variables
  REAL (KIND=dp), DIMENSION(2,nw, nostk)       :: radclrclds
  REAL (KIND=dp), DIMENSION(nw, nostk)         :: rad, &
       cfracwf, o3shiwf, &
       codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(nw, nalbwf, nostk) :: albwf
  REAL (KIND=dp), DIMENSION(nw, nl, nostk)     :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(nw, nflay       )  :: tauwf
  REAL (KIND=dp), DIMENSION(nw, nflay, nostk)  :: fozwf, &
       faerwf, faerswf, fcodwf, fsprswf, fraywf
  ! polarization correction variables
  ! @ lut-correction
  CHARACTER (LEN=maxchlen) :: VLDLUTdir
  LOGICAL :: do_plutcorr_after
  REAL (KIND=dp) :: albclrcld(nw, 2)
  REAL (KIND=dp), DIMENSION(nw, 2, nostk)      :: polerr
  ! @ on-correction  
  REAL (KIND=dp), DIMENSION(mpolcorr, nflay)   :: ptauwf,pfozwf, &
                                                  pfaerwf, pfaerswf, pfcodwf, pfsprswf, pfraywf
  REAL (KIND=dp), DIMENSION(mpolcorr)          :: prad, pctpwf, pcfracwf
  REAL (KIND=dp), DIMENSION(mpolcorr, nalbwf)  :: palbwf

  !cross section & optical depth variables
  REAL (KIND=dp), DIMENSION(nw)                 :: raycof, depol
  REAL (KIND=dp), DIMENSION(nw, nallgas, nflay) :: allcrs
  REAL (KIND=dp), DIMENSION(nw, nflay)          :: abscrs
  REAL (KIND=dp), DIMENSION(nallgas, nflay)     :: alleta, allcol
  REAL (KIND=dp), DIMENSION(nw, nflay) :: deltau, delsca, delray, delabs
  ! Others
  REAL (KIND=dp) :: e_loop1, e_loop2, e_s, e_n ,e1, e2,e_vlidort,e_pol, e_inter
  ! Save variables
  LOGICAL, SAVE                                 :: first = .TRUE.
  INTEGER, SAVE                                 :: nz, faer_lvl, nradcal, npolmod
  REAL (KIND=dp), DIMENSION(0:mflay),      SAVE :: ts, ps
  LOGICAL, DIMENSION(mflay),               SAVE :: aermsk
  LOGICAL, DIMENSION(:),ALLOCATABLE,       SAVE :: do_radcals, do_polcorrs
  INTEGER, DIMENSION(:),ALLOCATABLE,       SAVE :: radcal_idxs, polidx
  LOGICAL :: do_ssfullb295 = .FALSE.

  ! Added variables for VLIDORT-PCA
  CHARACTER (LEN=maxchlen) :: message
  CHARACTER (LEN=maxchlen), DIMENSION (GT_MAXMESSAGES) :: messages
  INTEGER :: nmessages
  LOGICAL :: Monitor_CPU = .false., do_vlidort_initialize = .false., do_jacobians=.TRUE.
  LOGICAL, PARAMETER :: do_debug_rtm=.false.
  LOGICAL :: do_atmos_linearization, do_surface_linearization, do_linearization

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=15), PARAMETER :: modulename = 'LIDORT_PROF_PCA'

  call cpu_time(e_s)
  IF (nw <= 1 .or. use_effcrs ) THEN     
    WRITE(*, *) modulename//'is it PCA-RT simulation ?'
  ENDIF
  errstat = pge_errstat_ok
  do_atmos_linearization = .false.
  do_surface_linearization = .false.
  do_fraywf = .false. ; fraywf=0.0D0

  IF (first) THEN
     ! ======================= Read LIDORT Control Input ==========================
     WRITE(*,'(A,A,3I2)') 'RTM=PCA', 'polcorr=', polcorr, the_str, nflay
     CALL LIDORT_Read_Config (ADJUSTL(TRIM(ctrdbdir))//'vlidort_control_vv2p7_pca.inp', & 
     aerosol, has_clouds, Inputs, problems, message)
     IF (problems) THEN 
         WRITE(*,*) 'Errors in lidort_read_config'
         print *, message
         STOP 1
     ENDIF 
     first = .FALSE.
  ENDIF

  ! ============= Overridden some control and atmospheric variables ============== 
  IF (num_iter == 0 ) THEN

     !**geometries related variables
     Inputs%Geometry%n_gems_geometries = 1
     Inputs%Geometry%gems_szas(1) = sza

     IF (sza >= 90.0 .OR. sza < 0) THEN
        WRITE(*, *) modulename, ' : SZA is >= 90 or < 0 !!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     Inputs%Geometry%gems_vzas(1) = vza
     Inputs%Geometry%gems_azms(1) = aza

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
     ps(1:nz) = exp((log(fps(1:nz)) + log(fps(0:nz-1))) / 2.0)  
     ts(1:nz) = (fts(1:nz) + fts(0:nz-1)) / 2.0

     IF (nwfc > 0) THEN 
       PRINT *, 'not implemented for nwfc > 0 in PCA' ; STOP 1
     ENDIF
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
      CALL BSPLINE(atmosprof(2, 0:nl), atmosprof(1, 0:nl), &
          nl+1, fzs(0:nz), ps(0:nz), nz+1, errstat)
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
      do_atmos_linearization = .TRUE.

  do_surface_linearization = do_albwf    
  IF (do_atmos_linearization .OR. do_albwf) THEN
     do_linearization = .TRUE.
  ELSE 
     do_linearization = .FALSE.
  ENDIF

  !call allocate_rt (n_totalatmos_wfs, nlayers, ncalcp)
  i = 0;  ozwfidx = 0; aodwfidx = 0; twaewfidx = 0; codwfidx = 0; sprswfidx = 0
  layer_vary_flag_cc(1:nz) = .FALSE.
  layer_vary_number_cc(1:nz) = 0

  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi) THEN
     i = i + 1;  ozwfidx = i
     profilewf_names_cc(i) = 'ozone volume mixing ratio------'

     layer_vary_flag_cc(1:nz) = varyprof(1:nz)
     WHERE (varyprof(1:nz) )
        layer_vary_number_cc(1:nz) = & 
        layer_vary_number_cc(1:nz) + 1
     ENDWHERE
  ENDIF
  IF ( do_taodwf .OR. do_saodwf) THEN
     i = i + 1; aodwfidx = i
     profilewf_names_cc(i) = 'aerosol extinction coefficient-'

     IF (do_taodwf) THEN
        layer_vary_flag_cc(nup2p(ntp)+1:nz1) = .TRUE.
        layer_vary_number_cc(nup2p(ntp)+1:nz1) = & 
        layer_vary_number_cc(nup2p(ntp)+1:nz1) + 1
     ENDIF

     IF (do_saodwf) THEN
        layer_vary_flag_cc(1:nup2p(ntp)) = .TRUE.
        layer_vary_number_cc(1:nup2p(ntp)) = & 
        layer_vary_number_cc(1:nup2p(ntp)) + 1
     ENDIF
  ENDIF
  IF ( do_twaewf) THEN
     i = i + 1; twaewfidx = i
     profilewf_names_cc(i) = 'aerosol scattering coefficient-'
     layer_vary_flag_cc(nup2p(ntp)+1:nz1) = .TRUE.
     layer_vary_number_cc(nup2p(ntp)+1:nz1) = & 
     layer_vary_number_cc(nup2p(ntp)+1:nz1) + 1
  ENDIF
  IF ( do_codwf) THEN
     i = i + 1; codwfidx = i
     profilewf_names_cc(i) = 'cloud extinction coefficient---'
     layer_vary_flag_cc(nctp:ncbp) = .TRUE.
     layer_vary_number_cc(nctp:ncbp) = & 
     layer_vary_number_cc(nctp:ncbp) + 1
  ENDIF

  IF ( do_sprswf ) THEN ! FIX ME
     ! Need to use jacobians wrt rayleigh OD to perform interpolation

     i = i + 1; sprswfidx = i; raywfidx = i
     profilewf_names_cc(i) = 'rayleigh optical thickness-----'
     layer_vary_flag_cc(nup2p(nsfc-1)+1:nz1) = .TRUE.
     layer_vary_number_cc(nup2p(nsfc-1)+1:nz1) = & 
     layer_vary_number_cc(nup2p(nsfc-1)+1:nz1) + 1
  ENDIF

  n_totalatmos_wfs = i
  n_surface_wfs = nalbwf   

  IF ( n_totalatmos_wfs   > 1) THEN
     WHERE(layer_vary_number_cc(1:nz) > 0) 
        layer_vary_number_cc(1:nz) = n_totalatmos_wfs 
     ENDWHERE
  ENDIF

  do_fozwf = .FALSE.
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi ) do_fozwf = .TRUE.
  do_faerwf = .FALSE.
  IF (do_taodwf .OR. do_saodwf) do_faerwf = .TRUE.

  e_vlidort = 0.0
  ! ==================== Get Ozone Absorption Cross Section ====================   
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

      ! O3/SO2 (use_so2dtcrs=.TRUE.) cross section: if do_tmpwf = .FALSE. and do_o3shi is false, 
      ! just need to get once for each retrieval
      ! Other trace gas cross section: just need to get it once for all the retrievals if no shifts 
      CALL GET_HRES_GASCRS_RAY(num_iter, nw, waves, nz1, ts(1:nz1), ts(1:nz1),nfgas,&
          allcol(1:nfgas, 1:nz1), frhos(1:nz1), &
          do_o3shi, o3shi, do_tmpwf,  &
          allcrs(1:nw, 1:nfgas, 1:nz1),raycof(1:nw), depol(1:nw), problems)

      abscrs(1:ncalcp, 1:nz1) = allcrs(1:ncalcp, 1, 1:nz1)

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

     ! o3 absorption coefficient at 370.2 nm with TOMS FWHM
     dads%o3(1, 1:nz1) = 0.0; dadt%o3(1, 1:nz1) = 0.0

     ! Weighted by solar flux
     !abscrs(1, 1:nz) = 9.1231787D-24 + (ts(1:nz) - 273.15) * 1.9005502D-25 + &
     !     (ts(1:nz) - 273.15)**2.0 * 1.2275286D-27     
     !raycof(1) = 2.3184501D-26; depol(1) = 0.030247913D0

     nfgas = 7
     CALL GET_ALB_OZCRS_RAY(nz1, ts(1:nz1), nfgas, allcrs(1, 1:nfgas, 1:nz1), raycof(1), depol(1), problems)

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

  IF (problems) THEN
        WRITE(*, *) modulename, ' : Problems in reading O3 XSec !!!'
        errstat = pge_errstat_error; RETURN
  ENDIF
  ! Initialize aerosol/cloud property profiles
  IF (num_iter == 0) THEN
     aersca(1:nz1) = 0.0; aerext(1:nz1) = 0.0
     aerasy(1:nz1) = 0.0;  aermoms(0:nmom, 1:ngksec, 1:nz1) = 0.0
  ENDIF

  IF (do_lambcld) THEN
    IF (nctp /= 0) cldmsk(nctp:ncbp)=.TRUE.
  ELSE
    cldmsk = .FALSE.
  ENDIF

  cldsca=0.0; cldext=0.0; cldasy=0.0;  cldmoms = 0.0

  ! senstitivity of absorption cross section to temperature, used 
  ! for calculating temperature wf directly with LIDORT
  !eta = 0.0             ! dummy variable here
  alleta = 0.0

  !******************** polarization correction setting  *********************
  IF (num_iter == 0) THEN 
    ! identify wavelengths for calc
    IF (allocated(do_radcals)) deallocate(do_radcals, radcal_idxs)    
    allocate (do_radcals(nw), radcal_idxs(nw))
    IF (nw>1) THEN 
      do_radcals(1:nw) = .FALSE.
      nradcal = ncalcp
      do_radcals(1:ncalcp) = .TRUE.
      radcal_idxs(1:ncalcp) = (/(i, i=1, ncalcp)/)
    ELSE
      do_radcals(1) = .TRUE.; nradcal = 1; radcal_idxs(1) = 1
    ENDIF

    ! identify npolmod
    IF (nw > 1) THEN
      IF ( polcorr == 3 .OR. (polcorr == 4 .AND. num_iter == 0) )  THEN
        npolmod = 2   ! Twice, one vector and one scalar
      ELSE
        npolmod = 1   ! Only once either scalar or vector
      ENDIF
    ELSE
      npolmod = 1
    ENDIF

   ! Determine wavelengths where exact polarization correction (NSTOKES: 4 vs 1) is calculated
   ! In UV1 (or between 270 and 310 nm): ~292 nm, ~298 nm, ~300 nm, ~302 nm, ~304 nm, ~306 nm, last wavelength 
   ! In UV2 (or between 310 and 340 nm): first, 1/4, middle and last wavelength
   ! So exact vector LIDORT calculation is done at 11 wavelengths.
   ! This option works when radiance interpolation option is turned on
   IF (allocated(do_polcorrs)) deallocate (do_polcorrs, polidx)
   allocate (do_polcorrs(nw), polidx(nw))
   npolcorr=0; do_polcorrs(1:nw) = .FALSE. ; polidx(1:nw) = 0; polcorr_idxs(:) = 0 
   IF ( (polcorr >= 3 .AND. polcorr <= 5) .AND. nw > 1 ) THEN
     call set_polcorr(numwin, winlim(1:numwin, :),nw, waves, do_radcals(1:nw), &
     npolcorr,do_polcorrs, polidx, polcorr_idxs)
   ENDIF
  ENDIF

  polerr = 1.0
  e_vlidort = 0.0
  call cpu_time(e_loop1)  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 1. simulation at everywavelengths
  ! 2. simulation at polcorr wavelengths if required
  
  IF (allocated(fr)) THEN 
   CALL deallocate_od
  ENDIF
  CALL allocate_od(nz, n_totalatmos_wfs)

  IF (npolmod == 2 .or. nw == 1) THEN
    nstokes = 3 ; nstreams = the_str
    ! Initialize output variables
    prad = 0.0 ; pcfracwf = 0.0
    IF ( do_linearization ) THEN
     IF ( do_fozwf )  pfozwf   = 0.0
     IF ( do_albwf)   palbwf   = 0.0
     IF ( do_codwf )  pfcodwf  = 0.0
     IF ( do_sprswf ) pfsprswf  = 0.0
     IF ( do_ctpwf )  pctpwf   = 0.0
     IF ( do_faerwf)  pfaerwf  = 0.0
     IF ( do_twaewf ) pfaerswf = 0.0
    ENDIF

    radclrclds(:,:,:) = 0.0
    DO ic = 1, 2
      wfrac (:) = 0.0
      DO ip = 1, npolcorr
        iw = polcorr_idxs(ip)
        lamda = waves(iw)
        IF (nwfc > 0 ) the_cfrac = wfcs(iw)  
        IF (ic == 1) THEN 
          do_clouds = .FALSE. ; frac = 1.0 - the_cfrac ! clear sky simulation
        ELSE
          do_clouds = .TRUE.  ; frac = the_cfrac ! cloud pixel simulation
        ENDIF       

        wfrac (ip) = frac
        IF (frac == 0.0) CYCLE
        do_prep = .false.

        IF ((ic == 1)  .OR. ( .NOT. do_lambcld ) ) THEN ! clear sky or treat cloud as a particle
          lambertian_albedo = albs(iw)
          nz1 = nfsfc-1 ; nlayers = nz1
          do_prep = .TRUE.
        ELSE ! lambertian cloud simulation
          nz1 = nctp -1 ; nlayers = nz1
          do_clouds = .FALSE.
          IF (the_cfrac == 1.0 .AND. nw /=1 ) THEN 
            lambertian_albedo = albs(iw) 
            lambcld_refl = lambertian_albedo
          ELSE
            lambertian_albedo = lambcld_refl
          ENDIF
          IF (.NOT. do_prep ) THEN 
            do_prep = .TRUE.
          ENDIF
        ENDIF
        Geophys%Surface%albedo(ip)      = lambertian_albedo

        IF (do_prep) THEN 
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
          IF (ip == 1) cldext0(nctp:ncbp) = cldext(nctp:ncbp) * gcq(actawin) / (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) 
          IF (useasy) THEN
            cldasy(nctp:ncbp) = gcasy(low) * (1.0 - xg) + gcasy(hgh) * xg
          ELSE
            DO i = 0, nmom
              DO j = 1, ngksec
                 cldmoms(i, j, nctp:ncbp) = gcmoms(low, i, j) &
                      * (1.0 - xg) + gcmoms(hgh, i, j) * xg 
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        CALL LIDORT_PROF_PREP(iw,lamda ,nstokes, nz1, n_totalatmos_wfs, &
             raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
             varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
             alleta(1:nfgas, 1:nz1), useasy, nmom, profilewf_names_cc(1:n_totalatmos_wfs), & 
             aerosol, aersca(1:nz1), aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
             do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
             cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1),absod(1:nfgas,1:nz1), problems)
        IF (problems) THEN 
          WRITE(*, *) ADJUSTL(TRIM(modulename))//': Errors in vlidort preparation'
          errstat = pge_errstat_error 
          RETURN 
        ENDIF
      
        deltau(ip, 1:nz1) = deltau_vert_input(1:nz1)
        delsca(ip, 1:nz1) = deltau_vert_input(1:nz1)*omega_total_input(1:nz1)
        delray(ip, 1:nz1) = delsca(ip,1:nz1)*fr(1:nz1)
        Geophys%Xsecs%Rayleigh_depol(ip)= depol(iw)
        Geophys%totalods%taudp(1:nz1,ip)= deltau_vert_input(1:nz1)
        Geophys%totalods%omega(1:nz1,ip)= omega_total_input(1:nz1)
        Geophys%totalods%taug(1:nz1, ip)= deltaug_vert_input(1:nz1)
        Geophys%totalods%fr(1:nz1, ip)  = fr(1:nz1) !rayod(1:nz1)/scaod(1:nz1)
        Geophys%totalods%fa(1:nz1, ip)  = fa(1:nz1) !aerod(1:nz1)/scaod(1:nz1)
        l_geophys%l_totalods%l_taudp(1:n_totalatmos_wfs, 1:nz1, ip)   = l_deltau_vert_input(1:n_totalatmos_wfs, 1:nz1)
        l_geophys%l_totalods%l_omega(1:n_totalatmos_wfs, 1:nz1, ip)   = l_omega_total_input(1:n_totalatmos_wfs, 1:nz1)
        ENDIF
        wfrac(ip) = frac
        Geophys%WavGrids%nwav           = ip
        Geophys%WavGrids%wav(ip)        = waves(iw) !@@@ Need to change
        Geophys%Xsecs%Rayleigh_depol(ip)= depol(iw)
        
    ENDDO ! loop of wavelength
        nw0 = ip-1
        IF (nw0 ==0) cycle  
        Inputs%RTMcontrol%NVLIDORT_nstreams = nstreams
        Inputs%RTMcontrol%NVlidort_nstokes  = nstokes
        Inputs%Atmosph%do_clouds   = do_clouds
        Inputs%Atmosph%do_aerosols = aerosol
        Geophys%Atmos%ngreek_moments_input   = nmom
        Geophys%Atmos%nlayers   = nz1
        Geophys%Atmos%Level_heights(0:nz1) = fzs(0:nz1)
        !Geophys%Aerosols%AEROSOL_LAYERFLAGS(0:nz1) =aermsk(1:nz1)
        !Geophys%Clouds%CLOUD_LAYERFLAGS(0:nz1) = cldmsk(1:nz1)
        RI%n_totalatmos_wfs = n_totalatmos_wfs
        RI%n_surface_wfs    = n_surface_wfs
        ! initialize output variable
        !RO%Stokes = 0.0
        !RO%LP_jacobians = 0.0
        !RO%LS_jacobians = 0.0

        CALL cpu_time(e1)
        CALL VLIDORT_Exact_master (Monitor_CPU, do_VLIDORT_initialize, do_Jacobians,  & !flag
           problems, NMessages , Messages) 
        CALL cpu_time(e2)
        e_vlidort = e_vlidort + e2 - e1
        IF (problems) THEN 
           WRITE(*, *) ADJUSTL(TRIM(modulename))//': Errors in vlidort calculation'
           errstat = pge_errstat_error
           RETURN
        ENDIF
        ! radiance
        Do istk = 1, 1 !nostk
          radclrclds(ic,1:nw0, istk) =  RO%Stokes(istk,1:nw0)*polerr(1:nw0, ic, istk)
          prad(1:nw0) = prad(1:nw0) + wfrac(1:nw0) * radclrclds(ic,1:nw0, istk)
        ENDDO
        IF (do_linearization ) THEN
         PostProcess0: DO istk = 1, 1 !nostk
        DO iwf = 1, n_totalatmos_wfs
          IF (iwf == ozwfidx .and. do_ozwf) THEN
            DO zz = 1, nz1 !nfsfc-1            
              pfozwf(1:nw0,zz) = pfozwf(1:nw0,zz) + &
              RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)/ozs(zz)
            ENDDO          
          ENDIF

          ! Weighting function with respect to aerosol/cloud optical depth at the last wavelength
          ! so as to keep the scaling since we are fitting the aod at that wavelength
          IF (iwf == aodwfidx .and. (do_taodwf .or. do_saodwf)) THEN
            DO zz = faer_lvl, nfsfc-1
              pfaerwf(1:nw0, zz)  = pfaerwf(1:nw0, zz) + &
              RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(1:nw0, ic, istk) / gaext(actawin,zz)
            ENDDO
          ENDIF                
          IF (iwf == twaewfidx .and. do_twaewf ) THEN
            DO zz = faer_lvl, nfsfc-1
              pfaerswf(1:nw0,zz)  = pfaerswf(1:nw0,zz) + &
              RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(1:nw0, ic, istk)/ &
                                          gasca(actawin,zz) * gaext(actawin,zz)
            ENDDO
          ENDIF
          IF (iwf == codwfidx .and. do_codwf) THEN
            DO zz = nctp, ncbp
              pfcodwf(1:nw0, zz)  = pfcodwf(1:nw0, zz) + &
              RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(1:nw0,ic, istk) /cldext0(zz) 
            ENDDO
          ENDIF
          IF (iwf == sprswfidx .and. do_sprswf ) THEN
            DO zz = nup2p(nsfc-1)+1, nfsfc-1
              pfsprswf(1:nw0, zz)  = pfsprswf(1:nw0,zz) + &
              RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(1:nw0,ic,istk) /(fps(zz) -fps(zz-1))
            ENDDO
          ENDIF
        ENDDO ! loop of iwf
        
        DO iwf = 1, n_surface_wfs
             IF (do_albwf) THEN 
                 IF (.NOT. do_lambcld) THEN
                   palbwf(1:nw0, iwf) = palbwf(1:nw0, iwf) + RO%SJac(iwf, istk, 1:nw0) * wfrac(1:nw0) 
                 ELSE IF (the_cfrac < 1.0) THEN
                   IF (ic == 1) palbwf(1:nw0, iwf) = palbwf(1:nw0, iwf) + RO%SJac(iwf, istk, 1:nw0) *  wfrac(1:nw0)
                 ELSE IF (the_cfrac == 1.0) THEN
                   IF (ic == 2) palbwf(1:nw0, iwf) = palbwf(1:nw0, iwf) + RO%SJac(iwf, istk, 1:nw0) *  wfrac(1:nw0)
                 ENDIF
             ENDIF
        ENDDO
        ENDDO PostProcess0 ! loop of istk
       ENDIF ! do_linearization
    ENDDO ! loop of cloud
    IF (do_cfracwf) THEN 
        pcfracwf(1:nw0) = radclrclds(2, 1:nw0, 1) - radclrclds(1, 1:nw0, 1)
    ENDIF
    IF (nw == 1) THEN
    rad (1:nw0,1) = prad(1:nw0)
    albwf(1:nw0,1:nalbwf, 1) = palbwf(1:nw0, 1:nalbwf)
    cfracwf(1:nw0,1) = pcfracwf(1:nw0)
    ENDIF
  ENDIF


  ! (2) PCA mode
  ! Initialize output variables
  rad = 0.0
  cfracwf = 0.0
  IF ( do_linearization) THEN 
    IF ( do_fozwf )  fozwf   = 0.0
    IF ( do_albwf)   albwf   = 0.0
    IF ( do_codwf )  fcodwf  = 0.0
    IF ( do_sprswf ) fsprswf  = 0.0
    IF ( do_ctpwf )  ctpwf   = 0.0
    IF ( do_faerwf)  faerwf  = 0.0
    IF ( do_twaewf ) faerswf = 0.0
  ENDIF
  
  nstokes = 1 ; nstreams = the_str
  IF (polcorr ==0 ) nstokes = 3
  
  radclrclds(:,:,:) = 0.0
  DO ic = 1, 2
    IF (ic == 1) THEN 
      do_clouds = .FALSE. ; frac = 1.0 - the_cfrac
    ELSE
      do_clouds = .TRUE.  ; frac = the_cfrac
    ENDIF

    wfrac (:) = frac
    IF (frac == 0.0) CYCLE

    woff = 0
    nw_pca = 0
    DO ipca = 1, npca
      w1 = npcapix(ipca,1)
      w2 = npcapix(ipca,2)
      which_win = which_pcabin(ipca)
      ip = 1
      DO iw = w1, w2
        do_prep = .FALSE.
        IF ((ic == 1)  .OR. ( .NOT. do_lambcld ) ) THEN
          lambertian_albedo = albs(iw)
          nz1 = nfsfc-1 ; nlayers = nz1 
          do_prep = .TRUE.
        ELSE
          nz1 = nctp -1 ; nlayers = nz1
          do_clouds = .FALSE.
          IF (the_cfrac == 1.0 .AND. nw /=1 ) THEN 
            lambertian_albedo = albs(iw) 
            lambcld_refl = lambertian_albedo
          ELSE
            lambertian_albedo = lambcld_refl
          ENDIF
          IF (.NOT. do_prep) THEN
            do_prep = .TRUE.
          ENDIF  
        ENDIF 
        IF ( polcorr==2) THEN 
          tmpalbs(ip) = lambertian_albedo
          IF (tmpalbs(ip) <0.001) tmpalbs(ip) = 0.001
          IF (tmpalbs(ip) >0.999) tmpalbs(ip) = 0.999
        ENDIF

        Geophys%Surface%albedo(ip)      = lambertian_albedo
        IF (do_prep) THEN
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

          !IF (num_iter == 0) THEN  ! Don't need to be updated, aerasy, aermoms should be save variables
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
          IF (ip == 1) cldext0(nctp:ncbp) = cldext(nctp:ncbp) * gcq(actawin) / (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) 
          IF (useasy) THEN
            cldasy(nctp:ncbp) = gcasy(low) * (1.0 - xg) + gcasy(hgh) * xg
          ELSE
            DO i = 0, nmom
              DO j = 1, ngksec
                 cldmoms(i, j, nctp:ncbp) = gcmoms(low, i, j) &
                      * (1.0 - xg) + gcmoms(hgh, i, j) * xg 
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        CALL LIDORT_PROF_PREP(iw,lamda ,nstokes, nz1, n_totalatmos_wfs, raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
                alleta(1:nfgas, 1:nz1), useasy, nmom, profilewf_names_cc(1:n_totalatmos_wfs), aerosol, aersca(1:nz1),      &
                aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), absod(1:nfgas, 1:nz1), problems)
        IF (problems) THEN 
             print * , 'problem in lidort prof prep'
             STOP 1
        ENDIF
        deltau(ip+woff, 1:nz1) = deltau_vert_input(1:nz1)
        delsca(ip+woff, 1:nz1) = deltau_vert_input(1:nz1)*omega_total_input(1:nz1) 
        delray(ip+woff, 1:nz1) = delsca(ip,1:nz1)*fr(1:nz1)
        Geophys%Xsecs%Rayleigh_depol(ip)= depol(iw)
        Geophys%totalods%taudp(1:nz1,ip)= deltau_vert_input(1:nz1)
        Geophys%totalods%omega(1:nz1,ip)= omega_total_input(1:nz1) 
        Geophys%totalods%taug(1:nz1, ip)= deltaug_vert_input(1:nz1)
        IF (hgaspos(1) > 0) Geophys%totalods%tauhg(1:nz1, ip)= absod(hgaspos(1),1:nz1) 
        IF (hgaspos(2) > 0) Geophys%totalods%tauhg(1:nz1, ip)= Geophys%totalods%tauhg(1:nz1, ip)+ absod(hgaspos(1),1:nz1) 
        Geophys%totalods%fr(1:nz1, ip)  = fr(1:nz1) !rayod(1:nz1)/scaod(1:nz1)
        Geophys%totalods%fa(1:nz1, ip)  = fa(1:nz1) !aerod(1:nz1)/scaod(1:nz1)
        l_geophys%l_totalods%l_taudp(1:n_totalatmos_wfs, 1:nz1, ip)   = l_deltau_vert_input(1:n_totalatmos_wfs, 1:nz1)
        l_geophys%l_totalods%l_omega(1:n_totalatmos_wfs, 1:nz1, ip)   = l_omega_total_input(1:n_totalatmos_wfs, 1:nz1)
        ENDIF  
        Geophys%WavGrids%nwav     = ip
        Geophys%WavGrids%wav(ip)  = waves(iw) !@@@ Need to change
        ip = ip + 1
    ENDDO ! loop of wavelength
        nw0 = ip-1
        IF (nw0 == 0) cycle 
        Inputs%RTMcontrol%NVLIDORT_nstreams = nstreams
        Inputs%RTMcontrol%NVlidort_nstokes = nstokes
        Inputs%Atmosph%do_clouds   = do_clouds
        Inputs%Atmosph%do_aerosols = aerosol
        Geophys%Atmos%ngreek_moments_input   = nmom
        Geophys%Atmos%nlayers   = nz1
        Geophys%Atmos%Level_heights(0:nz1) = fzs(0:nz1)
        !Geophys%Aerosols%AEROSOL_LAYERFLAGS(0:nz1) =aermsk(1:nz1)
        !Geophys%Clouds%CLOUD_LAYERFLAGS(0:nz1) = cldmsk(1:nz1)
        RI%n_totalatmos_wfs = n_totalatmos_wfs
        RI%n_surface_wfs    = n_surface_wfs
        ! initialize output variable
        !RO%Stokes = 0.0
        !RO%LP_jacobians = 0.0
        !RO%LS_jacobians = 0.0
               
        woff = woff + nw0   
        call cpu_time(e1)
        CALL VLIDORT_eofpc_master (Monitor_CPU, do_VLIDORT_initialize, do_Jacobians,  & !flag
           problems, NMessages , Messages) 
        call cpu_time(e2)
        nw_pca = nw_pca + n_call_vlidort
        e_vlidort = e_vlidort + e2 - e1
        IF (problems) THEN 
          WRITE(*, '(A, i3)') 'Problems in EOF VLIDORT maseter', Nmessages
          print *, messages(1)
          errstat = pge_errstat_error
          return
        ENDIF
        albclrcld(w1:w2,ic) = tmpalbs(1:nw0)
        ! radiance
        Do istk = 1, nostk
          radclrclds(ic,w1:w2, istk) =  RO%Stokes(istk,1:nw0)*polerr(w1:w2, ic, istk)
          rad(w1:w2, istk) = rad(w1:w2,istk) + wfrac(1:nw0) * radclrclds(ic,w1:w2, istk)
        ENDDO
        !WRITE(*,'(i3,i3,i5,f8.2, 2f8.2)') ic, ipca,nw0, sum(ro%stokes(1, 1:nw0)),-log(taug)
       IF (do_linearization ) THEN
         PostProcess: DO istk = 1, nostk
         DO iwf = 1, n_totalatmos_wfs
           IF (iwf == ozwfidx .and. do_ozwf) THEN
             DO zz = 1, nz1 !nfsfc-1            
               fozwf(w1:w2,zz,istk) = fozwf(w1:w2,zz, istk) + &
                  RO%PJac(zz, iwf, istk,1:nw0)*wfrac(1:nw0)*polerr(w1:w2, ic, istk)/ozs(zz)
               ! write(*,'(i3,2e15.7, 2f8.2)') zz,fozwf(1, zz, 1), profilewf(zz,1), frac, ozs(zz)
             ENDDO
           ENDIF
          ! Weighting function with respect to aerosol/cloud optical depth at the last wavelength
           ! so as to keep the scaling since we are fitting the aod at that wavelength
           IF (iwf == aodwfidx .and. (do_taodwf .or. do_saodwf)) THEN
             DO zz = faer_lvl, nfsfc-1
               faerwf(w1:w2, zz,istk)  = faerwf(w1:w2, zz, istk) + &
                RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(w1:w2, ic, istk) / gaext(actawin,zz)
             ENDDO
           ENDIF
                
           IF (iwf == twaewfidx .and. do_twaewf ) THEN
             DO zz = faer_lvl, nfsfc-1
               faerswf(w1:w2,zz, istk)  = faerswf(w1:w2,zz, istk) + &
                  RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(w1:w2, ic, istk)/ &
                  gasca(actawin,zz) * gaext(actawin,zz)
             ENDDO
           ENDIF
 
           IF (iwf == codwfidx .and. do_codwf) THEN
             DO zz = nctp, ncbp
               fcodwf(w1:w2, zz, istk)  = fcodwf(w1:w2, zz, istk) + &
                  RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(w1:w2,ic, istk) /cldext0(zz) 
             ENDDO
           ENDIF

           IF (iwf == sprswfidx .and. do_sprswf ) THEN
             DO zz = nup2p(nsfc-1)+1, nfsfc-1
                fsprswf(w1:w2, zz, istk)  = fsprswf(w1:w2,zz, istk) + &
                   RO%PJac(zz, iwf, istk, 1:nw0)*wfrac(1:nw0)*polerr(w1:w2,ic,istk) / &
                  (fps(zz) -fps(zz-1))
             ENDDO
           ENDIF
         ENDDO ! loop of iwf

         DO iwf = 1, n_surface_wfs
           IF (do_albwf) THEN 
             IF (.NOT. do_lambcld) THEN
               albwf(w1:w2, iwf, istk) = albwf(w1:w2, iwf, istk) + RO%SJac(iwf, istk, 1:nw0) * wfrac(1:nw0) 
             ELSE IF (the_cfrac < 1.0) THEN
               IF (ic == 1) albwf(w1:w2, iwf, istk) = albwf(w1:w2, iwf, istk) + RO%SJac(iwf, istk, 1:nw0) *  wfrac(1:nw0)
             ELSE IF (the_cfrac == 1.0) THEN
               IF (ic == 2) albwf(w1:w2, iwf, istk) = albwf(w1:w2, iwf, istk) + RO%SJac(iwf, istk, 1:nw0) *  wfrac(1:nw0)
             ENDIF
           ENDIF
         ENDDO
         ENDDO PostProcess ! loop of istk
       ENDIF ! do_linearization
    ENDDO ! loop of ipca
  ENDDO ! loop of cloud
  IF (do_cfracwf) THEN 
    nw0 = ncalcp
    cfracwf(1:nw0, 1:nostk) = radclrclds(2, 1:nw0, 1:nostk) - radclrclds(1, 1:nw0, 1:nostk)
  ENDIF

  CALL cpu_time (e_loop2)
  CALL deallocate_od
  !----------------------------------
  ! Post-pocessing after RTM simulations
  !-----------------------------------
  IF ( the_cfrac == 1.0 .AND. do_lambcld) THEN
     nz1 = nctp - 1
  ELSE
     nz1 = nfsfc - 1
  ENDIF
  nsprs = nup2p(nsfc-1)+1
  nw0 = ncalcp
  
  !----------------------------------
  ! polarization correction
  !-----------------------------------
  do_plutcorr_after = .false.
  IF (polcorr == 2 .and. nw > 1 .and. .NOT. do_plutcorr_after) THEN
    fidx = 1 ; lidx = nw0
    !do_ssfullb295 = .true.
    !IF (do_ssfullb295) fidx = MINVAL(MINLOC(waves(1:nw0), MASK=waves(1:nw0) > 280))
    delabs(fidx:lidx, 1:nz1) = deltau(fidx:lidx, 1:nz1) - delsca(fidx:lidx, 1:nz1)
    DO i = 1, nz1
     tauwf (fidx:lidx, i) = fozwf(fidx:lidx, i, 1)/abscrs(fidx:lidx, i)/rad(fidx:lidx,1)/du2mol
    ENDDO

    toz = SUM(ozs(1:nfsfc-1))
    VLDLUTdir='/home/jbak/data/GEMSTOOL/lutdatav2.8-r/LUT-48/'
    CALL polcorr_online_with_lut(num_iter,VLDLUTdir,lidx-fidx+1, nz1, nctp,nfsfc, nalbwf,&
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
    CALL polcorr_online(num_iter, polcorr, nw0, nz1, nctp, ncbp, nsprs, nalbwf, faer_lvl, &
       npolcorr, polcorr_idxs(1:npolcorr), do_fozwf, do_albwf, do_faerwf, &
       do_twaewf, do_codwf, do_sprswf, do_fraywf, do_cfracwf, &
       waves(1:nw0),  rad(1:nw0, 1), waves(polcorr_idxs(1:npolcorr)), prad(1:npolcorr), &
       tauwf(1:nw0, 1:nz1), ptauwf(1:npolcorr, 1:nz1), delabs(1:nw0, 1:nz1), &
       albwf(1:nw0, 1:nalbwf,1), palbwf(1:npolcorr, 1:nalbwf), fozwf(1:nw0, 1:nz1, 1), pfozwf(1:npolcorr, 1:nz1), &
       faerwf(1:nw0, 1:nz1, 1), pfaerwf(1:npolcorr, 1:nz1),faerswf(1:nw0, 1:nz1,1),pfaerswf(1:npolcorr, 1:nz1),&
       fcodwf(1:nw0, 1:nz1, 1), pfcodwf(1:npolcorr, 1:nz1),fsprswf(1:nw0, 1:nz1,1),pfsprswf(1:npolcorr, 1:nz1), &
       fraywf(1:nw0, 1:nz1, 1), pfraywf(1:npolcorr, 1:nz1),cfracwf(1:nw0, 1),      pcfracwf(1:npolcorr))
  ENDIF
  CALL cpu_time(e_pol) 

  IF (nw > 1 .and. do_debug_rtm) THEN
   close(21)
   rto%nw = nw0 ; rto%nl = nz1 ; rto%nos = nostk ; rto%nalbwf = nalbwf
   CALL set_rtmvar(.TRUE., rto)
   rto%wav   = radcwav(1:nw0)
   rto%alb   = albs(1:nw0)
   rto%rad   = rad(1:nw0, :)
   rto%albwf = albwf(1:nw0,1:nalbwf, :)
   rto%ozwf    = fozwf(1:nw0,1:nz1,:)
   CALL debug_rtm (21,rto, hcrs)
   CALL set_rtmvar(.false., rto)
   !CALL debug_taug(121, nw0, nz1, nfgas, radcwav(1:nw0),&
   !allcol(1:nfgas,1:nz1),allcrs(1:nw0, 1:nfgas, 1:nz1))
  ENDIF
  !-----------------------------------------------------------------
  ! convert simulations to higher resolution into instrument resoution
  !-----------------------------------------------------------------
  IF (nw > 1) THEN  ! hres==>radwf
   do_abs = .false.
   if (polcorr==2 .and. do_plutcorr_after) do_abs = .true.
   delabs(1:nw, 1:nz1) = deltau(1:nw, 1:nz1) - delsca(1:nw, 1:nz1)
   print * , 'pca interpolation scheme ! hres_radwf_inter_convol'
   CALL hres_radwf_inter_convol(nw, nz1, nctp, ncbp, nsprs, nalbwf, faer_lvl, & 
      do_albwf, do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf, do_tracewf,&
      do_o3shi, do_tmpwf,do_pslwf, waves, ozs(1:nz1), do_abs, delabs(1:nw,1:nz1), &
      rad(1:nw, 1),fozwf(1:nw, 1:nz1, 1), albwf(1:nw,1:nalbwf, 1), cfracwf(1:nw, 1), &
      faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
      fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), fraywf(1:nw,1:nz1, 1), errstat)
    IF (errstat == pge_errstat_error) RETURN
 
    IF (polcorr == 2  .and. do_plutcorr_after) THEN
    fidx = 1 ; lidx = n_radwvl_sav
    DO i = 1, nz1
     tauwf (fidx:lidx, i) = fozwf(fidx:lidx, i, 1)/ccrs%o3(fidx:lidx, i)/rad(fidx:lidx,1)/du2mol
    ENDDO
    VLDLUTdir='/home/jbak/data/GEMSTOOL/lutdatav2.8/LUT-W/'
    CALL polcorr_online_with_lut(num_iter, VLDLUTdir, lidx-fidx+1, nz1, nctp,nfsfc, nalbwf,&
         do_albwf, do_cfracwf, the_cfrac, albclrcld(fidx:lidx, 1:2),&
         sza, vza, aza, lat, toz, fps(0:nz1), ps(1:nz1), &
         radwvl_sav(fidx:lidx),  rad(fidx:lidx, 1),&
         tauwf(fidx:lidx, 1:nz1),delabs(fidx:lidx, 1:nz1), delsca(fidx:lidx,1:nz1), &
         albwf(fidx:lidx, 1:nalbwf,1),fozwf(fidx:lidx, 1:nz1, 1),  &
         cfracwf(fidx:lidx, 1))
    ENDIF
  ENDIF
  CALL cpu_time(e_inter)
  
  !-----------------------------------------------------------------
  ! Calculate desired weighting functions 
  ! 1. convert profilewf at fps to at umkp
  ! 2. derive tracegas weighting function 
  !----------------------------------------------------------------
  nw0 = n_radwvl_sav

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

  !------------------------------------
  ! copying local variables to output variables
  !------------------------------------
  rto%nw = nw0 ; rto%nl = nl ; rto%nos = nostk ; rto%nalbwf= nalbwf
  CALL set_rtmvar(.TRUE., rto)
  rto%rad   = rad(1:nw0, :)
  rto%albwf   = albwf(1:nw0,:,:)
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

  !--------------------------------
  ! Check running time in RT simulation
  !--------------------------------
  RETURN
  IF (scnwrt .and. nw > 1 ) THEN
  call cpu_time(e_n)
  PRINT * , 'N of VLIDORT CALL----------------:', nw_pca
  PRINT * , 'N of simulation wavelength-------:', nw
  PRINT * , 'N of OMI wavelength--------------:', nw0
  print * , ' START------------------END (PCA):', e_n     - e_s
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
END SUBROUTINE LIDORT_PROF_ENV_PCA

SUBROUTINE LIDORT_PROF_PREP (iw, lamda, nstokes, nlayers, n_totalatmos_wfs, raycof, depol, zsgrid, airgrid,  varyprof, &
     ngas, gasin, abscrs, gascol, eta, useasy, nmoms, profilewf_names, &          
     do_aerosols, aersca, aerext, aerasy, aermoms, aermsk, &
     do_clouds, cldsca, cldext, cldasy, cldmoms, cldmsk, absod, problems)
  
  USE OMSAO_precision_module
  USE ozprof_data_module, ONLY : maxgksec, maxgkmatc, ngksec, ngkmatc
  USE OMSAO_errstat_module
  USE vlidort_PARS, ONLY: omega_smallnum

  IMPLICIT NONE  
   
  !===============================  Define Variables ===========================
  ! Input variables
  INTEGER, INTENT(IN)                     :: ngas, nmoms, iw, nstokes, nlayers, n_totalatmos_wfs
  INTEGER, INTENT(IN), DIMENSION(ngas)    :: gasin
  LOGICAL, INTENT(IN), DIMENSION(nlayers) :: cldmsk, aermsk, varyprof
  LOGICAL, INTENT(IN)                     :: useasy
  REAL (KIND=dp), INTENT(IN)              :: raycof, depol, lamda
  REAL (KIND=dp), DIMENSION(0:nlayers), INTENT(IN) :: zsgrid
  REAL (KIND=dp), DIMENSION(nlayers), INTENT(IN)   :: airgrid,  &
       aersca, aerext, aerasy, cldsca, cldext, cldasy
  REAL (KIND=dp), DIMENSION(0:nmoms, maxgksec, nlayers), INTENT(IN) :: aermoms
  REAL (KIND=dp), DIMENSION(0:nmoms, maxgksec, nlayers), INTENT(IN) :: cldmoms
  REAL (KIND=dp), DIMENSION(ngas, nlayers),              INTENT(INOUT) :: abscrs
  REAL (KIND=dp), DIMENSION(ngas, nlayers),              INTENT(IN) :: gascol, eta
  CHARACTER(len=31), DIMENSION(n_totalatmos_wfs),        INTENT(IN) :: profilewf_names
  REAL (KIND=dp), DIMENSION(ngas, nlayers),              INTENT(OUT):: absod

  ! Output variables
  LOGICAL, INTENT(OUT)   :: problems

  ! Modified variables
  LOGICAL, INTENT(INOUT) :: do_aerosols, do_clouds

  ! Local variables
  INTEGER, PARAMETER     :: maxngas = 7, maxscatter=3, allngas = 9
  INTEGER, DIMENSION(maxgkmatc), PARAMETER :: &
       greekmat_idxs = (/1, 2, 5, 6, 11, 12, 15, 16/), &
       phasmoms_idxs = (/1, 5, 5, 2, 3, 6, 6, 4/)
  INTEGER  :: ui, i, j, k, q, ig, nscatter, idx, cldidx, aeridx, nactgksec, nactgkmatc
  INTEGER, DIMENSION(allngas)                     :: absin
  REAL (KIND=dp) :: scaco_r, absco_r, omega, extco_r, extco, scaco, &
                    extco_a, scaco_a, extco_c, scaco_c, pvar, j0, j1
  REAL (KIND=dp), DIMENSION(maxscatter)           :: scaco_input
  REAL (KIND=dp), DIMENSION(nlayers)              :: extconf

!  scattering phase function realted variables are calculated in!  vlidort_eof_master
   !Optical Properity
!  REAL(fpk), dimension ( 0:nmoms, nlayers, MAXSTOKES_SQ ) :: GREEKMAT_TOTAL_INPUT
!  REAL(fpk), dimension(n_totalatmos_wfs,nlayers) :: L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT
!  REAL(fpk), dimension(n_totalatmos_wfs, 0:nmoms, nlayers, MAXSTOKES_SQ) :: L_GREEKMAT_TOTAL_INPUT
!  REAL (KIND=dp), DIMENSION(0:nmoms, 1:maxgksec, maxscatter):: phasmoms_input
!  REAL (KIND=dp), DIMENSION(0:nmoms, 1:maxgksec) :: phasmoms_total_input
!  REAL (KIND=dp), DIMENSION(n_totalatmos_wfs, 0:nmoms, 1:maxgksec) :: l_phasmoms_total_input
  LOGICAL, SAVE :: first = .TRUE.
  
  ! ========================== Check for Input ==================================
  problems = .FALSE.
  !IF (MAXVAL(gasin) > allngas) THEN
  !   WRITE(www_lun, *) 'Not weighting functions are implemented for all gases!!!'
  !   problems = .TRUE.; RETURN
  !ENDIF

  IF (first) THEN
     IF (.NOT. useasy .AND. nmoms > gt_maxmoments) THEN
        WRITE(www_lun, *) 'Need to increase maxmoments_input for aerosols/clouds!!!'
        problems = .TRUE.; RETURN
     ENDIF
     
     ! This only needs to be initialized once
!     phasmoms_input        = ZERO
!     phasmoms_total_input  = ZERO
!     greekmat_total_input  = ZERO
     l_deltau_vert_input   = ZERO
     l_omega_total_input   = ZERO
!     l_greekmat_total_input= ZERO 
!     l_phasmoms_total_input= ZERO
     
     first =.FALSE.
  ENDIF
  
  IF (NSTOKES == 1) THEN
     nactgksec = 1;  nactgkmatc = 1
  ELSE
     nactgksec = ngksec; nactgkmatc = ngkmatc
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
 ! phasmoms_input(0, 1, 1) = ONE
 ! phasmoms_input(2, 1, 1) = (ONE - depol) / (GTTWO + depol)  
 ! IF (nactgksec == 6) THEN
 !    phasmoms_input(2, 2, 1) = 6.0D0 * phasmoms_input(2, 1, 1)
 !    phasmoms_input(2, 5, 1) = -SQRT(6.0D0) * phasmoms_input(2, 1, 1)
 !    phasmoms_input(1, 4, 1) = 3.0D0 * (ONE - 2.0D0 * depol) / (GTTWO + depol)
 ! ENDIF
 
  DO i = 1, nlayers   
     ! Rayleigh scattering
     scaco_r = raycof * airgrid(i)
     IF (any(abscrs(1:ngas, i) < 0)) THEN 
        DO ig = 1, ngas 
           IF (abscrs(ig, i) < 0) abscrs(ig,i) = 0.0
        ENDDO
     ENDIF

     ! Gas absorption
     absod(1:ngas, i) = abscrs(1:ngas, i) * gascol(1:ngas, i)
     absco_r = SUM(absod(1:ngas, i))         
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
          ! IF (.NOT. useasy) THEN
          !    phasmoms_input(0:nmoms, 1:nactgksec, nscatter) = cldmoms(0:nmoms, 1:nactgksec, i)                         
          ! ELSE ! use H_G function
          !    phasmoms_input(0, 2, nscatter) = ONE
          !    j0 = ONE             
          !    DO j = 1, nmoms
          !       j1 = REAL(2*j+1, KIND=dp)
          !       phasmoms_input(j, 2, nscatter) = (j1/j0) * cldasy(i) * phasmoms_input(j-1, 2, nscatter)
          !       j0 = j1
          !    ENDDO
          ! ENDIF
           
        ENDIF  ! end clouds
        
        IF (do_aerosols .AND. aermsk(i)) THEN
           nscatter = nscatter + 1
           extco = extco + aerext(i)
           scaco_input(nscatter) = aersca(i)
           extco_a = aerext(i);  scaco_a = aersca(i)
           aeridx = nscatter

           ! get phase moments for aerosols
          ! IF (.NOT. useasy) THEN  
          !    phasmoms_input(0:nmoms, 1:nactgksec, nscatter) = aermoms(0:nmoms, 1:nactgksec, i)  
          ! ELSE ! use H_G function
          !    phasmoms_input(0, 2, nscatter) = ONE
          !    j0 = ONE             
          !    DO j = 1, nmoms
          !       j1 = REAL(2*j+1, KIND=dp)
          !       phasmoms_input(j, 2, nscatter) = (j1/j0) * aerasy(i) * phasmoms_input(j-1, 2, nscatter)
          !       j0 = j1
          !    ENDDO
          ! ENDIF
        ENDIF  ! end aerosols
     !ENDIF     ! end non-rayleigh
     
     ! setup LIDORT input for tau and omega
     scaco = SUM(scaco_input(1:nscatter))
     omega = scaco / extco

     IF (omega < OMEGA_SMALLNUM) omega = OMEGA_SMALLNUM 
     IF (omega > 1.0 - OMEGA_SMALLNUM) omega = 1.0 - OMEGA_SMALLNUM

     deltau_vert_input(i) = extco
     deltaug_vert_input(i) = absco_r
     omega_total_input(i) = omega

     ! sum up phase moments as required in LIDORT
     !DO j = 0, nmoms
     !   DO k = 1, nactgksec
     !      phasmoms_total_input(j, k) = SUM(phasmoms_input(j, k, 1:nscatter) &
     !           * scaco_input(1:nscatter)) / scaco
     !   ENDDO
     !ENDDO
     !phasmoms_total_input(ngreek_moments_input+1:maxmoments, 1:maxgksec) = 0.0  
     
     ! Set up greek scattering matrix for each moment 

     !greekmat_total_input(0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
     !     phasmoms_total_input(0:nmoms, phasmoms_idxs(1:nactgkmatc))
     !IF ( nactgkmatc > 1 ) greekmat_total_input(0:nmoms, i, 15) &
     !     = -greekmat_total_input(0:nmoms, i, 15)!

     ! This should always be 1, but may be slightly different due to numerical truncation
     !greekmat_total_input(0, i, 1) = 1.0  


     !IF (do_simulation_only .OR. .NOT. varyprof(i)) THEN   ! no linearition
     !   ! zero out quantity for safety
     !   layer_vary_flag(i) = .FALSE.
     !   layer_vary_number(i) = 0	
     !   !l_deltau_vert_input(:, i) =   ZERO
     !   !l_omega_total_input(:, i) = ZERO
     !   !l_greekmat_total_input(:, : , :, i) = ZERO        
     !
     !ELSE
     !   layer_vary_flag(i) = .TRUE.
     !   layer_vary_number(i) = n_totalatmos_wfs
     !The above part has been taken care of in routine: lidort_prof_env.f90
         
     DO q = 1, n_totalatmos_wfs
        IF ( profilewf_names(q) == 'totaltauabs--------------------' ) THEN
           l_deltau_vert_input(q, i) = + absco_r !/ extco
           l_omega_total_input(q, i) = - absco_r *omega/ extco

        !  w.r.t ozone volume mixing ratio: 1
        !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'ozone volume mixing ratio------' ) THEN
           idx = absin(1)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, 0:maxmoments , i, 1:16) = ZERO
           !  w.r.t NO2 volume mixing ratio: 2
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'NO2 volume mixing ratio------' ) THEN
           idx = absin(2)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O2 volume mixing ratio: 8
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O2 volume mixing ratio------' ) THEN
           idx = absin(8)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O4 volume mixing ratio: 3
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O4 volume mixing ratio------' ) THEN
           idx = absin(3)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.;  RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t BrO volume mixing ratio: 4
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'BrO volume mixing ratio------' ) THEN
           idx = absin(4)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco

           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t H2O volume mixing ratio: 9
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'H2O volume mixing ratio------' ) THEN
           idx = absin(9)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t SO2 volume mixing ratio: 5
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'SO2 volume mixing ratio------' ) THEN
           idx = absin(5)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t HCHO volume mixing ratio: 6
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'HCHO volume mixing ratio------' ) THEN
           idx = absin(6)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t OCLO volume mixing ratio: 7
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'OCLO volume mixing ratio------' ) THEN
           idx = absin(7)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_deltau_vert_input(q, i) = + absod(idx, i) !/ extco
           l_omega_total_input(q, i) = - absod(idx, i) *omega/ extco

           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average temperature of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
           !  Assume no effects on air density              
        ELSE IF ( profilewf_names(q) == 'average temperature of layer---' ) THEN
           PRINT * , 'NEED WORK MORE'
           STOP 1
           l_omega_total_input(q, i) = - SUM(absod(:, i) * eta(:, i)) / extco
           l_deltau_vert_input(q, i) = + SUM(absod(:, i) * eta(:, i)) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average pressure of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
        ELSE IF ( profilewf_names(q) == 'average pressure of layer------' ) THEN
           PRINT * , 'NEED WORK MORE'
           STOP 1
           pvar = extco_r/ extco
           l_omega_total_input(q, i) = ((ONE - pvar) * scaco_input(1) - &
                pvar * (scaco - scaco_input(1))) / scaco
           l_deltau_vert_input(q,i) = extco_r / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t rayleigh soptical thickness
           ! xliu: August 12, 2008
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        ELSE IF ( profilewf_names(q) == 'rayleigh optical thickness-----' ) THEN
           pvar = scaco_r / extco    
           !l_omega_total_input(q,i) = (1.0 - omega) * scaco_r / extco / omega
           !l_deltau_vert_input(q,i) = pvar

           l_deltau_vert_input(q,i) = scaco_r
           l_omega_total_input(q,i) = (1.0 - omega) * pvar !/ omega

           !l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
           !DO j = 0, nmoms
           !  DO k = 1, nactgksec
           !     IF (phasmoms_total_input(j, k) /= 0.0) THEN
!          !        l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, 1) - phasmoms_total_input(j, k) ) &
!          !              / phasmoms_total_input(j, k) * scaco_r / scaco
           !        l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, 1) - phasmoms_total_input(j, k) ) *pvar
           !     ELSE
           !        l_phasmoms_total_input(q, j, k) = 0.0
           !     ENDIF
           !     !if ( i == 10 .and. iw == 1) print * ,k, l_phasmoms_total_input(q,0:2, k) 
           !  ENDDO            
          !ENDDO
          !l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
          !     l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
          !IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:nmoms, i, 15) &
          !     = - l_greekmat_total_input(q, 0:nmoms, i, 15)          
        ELSE IF ( profilewf_names(q) == 'rayleigh scattering coefficient' ) THEN
           !  xliu: April 13, 2007 
           !  Still need to consider the variation in phase function  
           pvar = scaco_r / extco
           l_omega_total_input(q, i) = ((ONE - pvar) * scaco_input(1) - &
                pvar * (scaco - scaco_input(1)) ) / scaco
           l_omega_total_input(q,i) = pvar
           l_deltau_vert_input(q,i)     = scaco_r / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t aerosol extinction coefficient / aerosol optical thickness
           !  aerosol scattering albedo does not change   
           !  xliu: April 13, 2007 (consider the variation in phase function)
        ELSE IF ( profilewf_names(q) == 'aerosol extinction coefficient-' ) THEN
           IF (aeridx > 0) THEN
              l_deltau_vert_input(q,i) = + extco_a / extco
              l_omega_total_input(q,i) = (scaco_a  / extco_a - omega ) &
                   / omega  * extco_a / extco
             ! l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
             ! DO j = 0, nmoms
             !    DO k = 1, nactgksec
             !       IF (phasmoms_total_input(j, k) /= 0.0) THEN
             !          l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aeridx) &
             !               - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
             !               * scaco_a / scaco
             !       ELSE
             !          l_phasmoms_total_input(q, j, k) = 0.0
             !       ENDIF
             !    ENDDO
             ! ENDDO
             ! l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
             !      l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
             ! IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:nmoms, i, 15) &
             !      = - l_greekmat_total_input(q, 0:nmoms, i, 15)
              !print *, maxval(l_greekmat_total_input), minval(l_greekmat_total_input)
              !print *, i, l_deltau_vert_input(q,i), l_omega_total_input(q,i)
              !WRITE(www_lun, '(6D14.6)') (l_phasmoms_total_input(1, j, 1:nactgksec), j = 0, ngreek_moments_input)
              !STOP 1
           ENDIF
           !  w.r.t  aerosol scattering coefficient / single scattering albedo
           !  aerosol optical thickness will not change
        ELSE IF ( profilewf_names(q) == 'aerosol scattering coefficient-' ) THEN
           IF (aeridx > 0) THEN
              l_deltau_vert_input(q,i) = ZERO
              l_omega_total_input(q,i) = scaco_a  / omega / extco
             ! l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
             ! DO j = 0, nmoms
             !    DO k = 1, nactgksec
             !       IF (phasmoms_total_input(j, k) /= 0.0) THEN
             !          l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aeridx) &
             !               - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
             !               * scaco_a / scaco
             !       ELSE
             !          l_phasmoms_total_input(q, j, k) = 0.0
             !       ENDIF
             !    ENDDO
             ! ENDDO
             ! l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
             !      l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
             ! IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:nmoms, i, 15) &
             !      = - l_greekmat_total_input(q, 0:nmoms, i, 15)
           ENDIF
           !   w.r.t cloud extinction coefficient / optical thickness
        ELSE IF ( profilewf_names(q) == 'cloud extinction coefficient---' ) THEN
           IF (cldidx > 0) THEN
              l_deltau_vert_input(q,i) = + extco_c / extco
              l_omega_total_input(q,i) = (scaco_c  / extco_c - omega ) &
                   / omega  * extco_c / extco
             ! l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
             ! DO j = 0, nmoms
             !    DO k = 1, nactgksec
             !       IF (phasmoms_total_input(j, k) /= 0.0) THEN
             !          l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldidx) &
             !               - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
             !               * scaco_a / scaco
             !       ELSE
             !          l_phasmoms_total_input(q, j, k) = 0.0
             !       ENDIF
             !    ENDDO
             ! ENDDO
             ! l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
             !      l_phasmoms_total_input(0, 1:nmoms, phasmoms_idxs(1:nactgkmatc))
             ! IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:nmoms, i, 15) &
             !      = - l_greekmat_total_input(q, 0:nmoms, i, 15)
           ENDIF
           
           !  w.r.t clouds scattering coefficient
        ELSE IF ( profilewf_names(q) == 'cloud scattering coefficient---' ) THEN
           IF (cldidx > 0) THEN
              l_deltau_vert_input(q,i) = ZERO
              l_omega_total_input(q,i) = scaco_c  / omega / extco
             ! l_greekmat_total_input(q, 0:nmoms, i, :) = ZERO
             ! DO j = 0, nmoms
             !    DO k = 1, nactgksec
             !!       IF (phasmoms_total_input(j, k) /= 0.0) THEN
              !         l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldidx) &
              !              - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
              !              * scaco_a / scaco
              !!      ELSE
               !        l_phasmoms_total_input(q, j, k) = 0.0
               !     ENDIF
               !  ENDDO
             ! ENDDO
             ! l_greekmat_total_input(q, 0:nmoms, i, greekmat_idxs(1:nactgkmatc)) = &
             !      l_phasmoms_total_input(q, 0:nmoms, phasmoms_idxs(1:nactgkmatc))
             ! IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:nmoms, i, 15) &
             !      = - l_greekmat_total_input(q, 0:nmoms, i, 15)
           ENDIF
        ENDIF     ! end selection of weighting function 
        l_fr(q,i) = ZERO
        l_fc(q,i) = ZERO
        l_fa(q,i) = ZERO                
     ENDDO           ! n_totalatmos_wfs loop
  !ENDIF             ! end of do_linearization    
     fr(i) = scaco_r/scaco     
     fa(i) = scaco_a/scaco
     fc(i) = scaco_c/scaco
  ENDDO              ! layer loop
  
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
  !USE OMSAO_parameters_module, ONLY  : maxchlen, rearth
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
  !LOGICAL                   :: fail
  !CHARACTER (len=maxchlen)  :: message, trace

!  VLIDORT_ModIn%MSunRays%TS_N_SZANGLES = 1
!  VLIDORT_ModIn%MSunRays%TS_SZANGLES = sza
  IF (sza >= 90.0 .OR. sza < 0) THEN
     STOP 'GET_SLANT_TAU: SZA is >= 90 or < 0!!!'
  ENDIF

  nlayers = nz
  !! FIX ME taugrid_input(0:nz) = tauin

!  VLIDORT_FixIn%Chapman%TS_height_grid(0:nz) = zs
!  VLIDORT_ModIn%MChapman%TS_earth_radius = rearth
  IF (nz > maxlayers) THEN
     STOP 'LIDORT_PROF_ENV: # of layers exceeded allowed !!!'
  ENDIF
  STOP ' VLIDORT_CAPMAN NOT implemented, which is from v2p4 !!!'
  !CALL VLIDORT_CHAPMAN(fail, message, trace)
  tauout = 0.0
  DO i = 1, nz
     DO j = 1, i
        !! FIX ME tauout(i) = tauout(i) + VLIDORT_Work_Miscellanous%deltau_slant(i, j, 1)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE GET_SLANT_TAU

END MODULE m_lidort_env_vv2p7_pca

