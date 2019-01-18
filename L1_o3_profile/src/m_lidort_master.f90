!
module m_lidort_master
  
! ******************************************************************************************
! use_effcrs 
!  1) true, outputs (radiance and weighting function) correspond to input wavelengths. 
!     geto3_crs ==> abcrs *normc
!  2) false, i.e., use high resolution cross  section in the radiative transfer calculation and 
!          convolve with slot functions after radiative transfer calculation, 
!          waves do not correspond to the output rad and weighting function, 
!          waves(1:nw) store the wavelengths (waves(1:ncalcp) where radiance calculation 
!          are done (waves(ncalp+1:nw) = 0.0. And outputs corresponds to measurents (1:ns).
!          Here nw is max(ncalcp, ns).
!    
! ******************************************************************************************

  PUBLIC  lidort_prof_env
  PRIVATE lidort_prof_prep!, get_slant_tau !, prepare_spherical
  CONTAINS

  SUBROUTINE LIDORT_PROF_ENV (do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary, &
  do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
          do_sprswf,  do_so2zwf,do_pslwf, &
          nw, waves, nos, o3shi, sza, vza, aza, nl, ozprof, tprof, &
          n0alb, albarr, albpmin, albpmax,  vary_sfcalb, walb0s, &
          n0wfc, wfcarr, wfcpmin, wfcpmax, wfc0s, nostk,  &
          albwf, ozwf,  tmpwf, o3shiwf, cfracwf, codwf, ctpwf, &
          taodwf, twaewf, saodwf, sprswf, so2zwf, rad, errstat)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,ONLY  : maxlay, max_fit_pts, du2mol, rearth, max_spec_pts
    USE OMSAO_indices_module,   ONLY  : so2_idx, so2v_idx, o2o2_idx, o2_idx, o2t2_idx, &
                                        h2o_idx, h2ot2_idx, refspec_strings
    USE OMSAO_variables_module, ONLY  : currloop, the_surfalt, band_selectors, &
       fitvar_rad, mask_fitvar_rad, fitvar_rad_apriori, &
       refspec_norm, &
       numwin, winlim, nradpix, radwvl_sav, n_radwvl_sav, &
       refidx, &
       database_pslwf, max_psl, npsl, tabdir, rtmdbg, rtm_unit
    USE ozprof_data_module,     ONLY : num_iter, VlidortNstream, &
       nflay, mflay, ntp, nfsfc, nsfc, nup2p, nt_fit, &
       atmosprof, fts, fps, fzs, fozs, frhos, &
       ngas, nallgas, gasidxs, fgasidxs,fgaspos, gassidxs, & 
       do_tracewf, tracegas, mgasprof,      &
       do_lambcld, lambcld_refl, has_clouds, aerosol, useasy, strat_aerosol, &
       nwfc, ncbp, nctp, the_cbeta,  nmom, maxgksec, ngksec, actawin, aerwavs, &
       the_cod, the_cfrac, taodind, saodind, twaeind, sprsind, &
       gaext, gasca, gaasy, gamoms, gcq, gcasy, gcmoms, &
       raycof, depol, tropsca, tropaod, tropwaer, strataod, stratsca, &
       np1, np2, &
       do_simu, radcalwrt, &
       use_effcrs, ncalcp, do_radinter, polcorr, & 
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       no2idx, broidx, hchoidx, so2idx, o4idx, o2idx, h2oidx, &
       so2crsidx, o2crsidx, o4crsidx, o2crsidx, h2ocrsidx

    USE m_lidort_util, ONLY:  hres_radwf_inter_convol, &
        get_tracegas_wf, radwf_interpol, polcorr_online
    USE m_ezspline_interpolation, only: bspline, bspline1
    USE m_get_xcrs, only: get_alb_ozcrs_ray,  get_hres_gascrs_ray, &
        get_effres_gascrs_ray
!    USE m_vlidort90_include
    USE OMSAO_errstat_module

    IMPLICIT NONE

    !===============================  Define Variables ===========================
    ! Include files of dimensions and numbers
    INCLUDE 'VLIDORT.PARS'

    ! Include files of input variables
    INCLUDE 'VLIDORT_INPUTS.VARS'
    INCLUDE 'VLIDORT_SETUPS.VARS'
    INCLUDE 'VLIDORT_L_INPUTS.VARS'
    INCLUDE 'VLIDORT_BOOKKEEP.VARS'

    ! Include files of result variables
    INCLUDE 'VLIDORT_RESULTS.VARS'
    INCLUDE 'VLIDORT_L_RESULTS.VARS'

    !INTEGER, PARAMETER :: ngas=1 

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN)  :: nw, nl, nos, n0alb, nostk, n0wfc
    LOGICAL, INTENT(IN)  :: do_ozwf, do_albwf, do_tmpwf, do_o3shi, do_taodwf, vary_sfcalb, &
         do_twaewf, do_saodwf, do_cfracwf, do_codwf, do_ctpwf, do_sprswf, do_so2zwf, do_pslwf
    INTEGER, INTENT(OUT) :: errstat
    INTEGER, DIMENSION(n0alb), INTENT(IN)        :: albpmax, albpmin
    INTEGER, DIMENSION(n0wfc), INTENT(IN)        :: wfcpmax, wfcpmin
    LOGICAL, DIMENSION(nl), INTENT(IN)           :: ozvary
    REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)   :: waves, walb0s, wfc0s
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
    LOGICAL :: problems,do_clouds, do_fozwf, do_faerwf, do_fraywf,  &
               do_o3hwe, do_o3spk
    INTEGER :: ic, iw, i, j, k, kk, jj, jk, idum, low, hgh, fidx, lidx, &
         nz1,nfgas1, nstep, istk, npolmod, ipol, nsprs, nw0
    INTEGER :: ozwfidx, aodwfidx, twaewfidx, codwfidx, sprswfidx, raywfidx
    LOGICAL, DIMENSION(nflay)                    :: cldmsk, varyprof!, aermsk
    REAL (KIND=dp)                               :: lamda, xg, frac, toz, temp, aodscl, waerscl
    REAL (KIND=dp), DIMENSION(0:nflay)           :: ozs, delps
    REAL (KIND=dp), DIMENSION(nflay)             :: cldsca, cldext, cldasy, cldext0, aersca, aerext, aerasy
    REAL (KIND=dp), DIMENSION(0:nmom,maxgksec,nflay) :: cldmoms, aermoms
    REAL (KIND=dp), DIMENSION(nw)                :: albs, tmpalbs, wfcs
    REAL (KIND=dp), DIMENSION(nw, 2, nostk)      :: polerr
    REAL (KIND=dp), DIMENSION(nw, nostk)         :: radclr, radcld
    !REAL (KIND=dp), DIMENSION(n_radwvl_sav)      :: swaves
    REAL (KIND=dp), DIMENSION(2, nostk)          :: radclrcld
    REAL (KIND=dp), DIMENSION(nw, nflay, nostk)  :: fozwf, faerwf, faerswf, fcodwf, fsprswf, fraywf
    REAL (KIND=dp), DIMENSION(nw, nflay)         :: pfozwf, pfaerwf, pfaerswf, pfcodwf, pfsprswf, pfraywf
    REAL (KIND=dp), DIMENSION(nw)        :: prad, palbwf, pctpwf, pcfracwf
    !REAL (KIND=dp), DIMENSION(nw)        :: tmprefwav, tmprefspec

    INTEGER, DIMENSION (nallgas)                  :: gasin   
    INTEGER, DIMENSION (5)                        :: tmp_gaspos ! used to sort gases when nw=1
    REAL (KIND=dp), DIMENSION(nw, nallgas, nflay) :: allcrs
    REAL (KIND=dp), DIMENSION(nallgas, nflay)     :: alleta, allcol  
    REAL (KIND=dp), DIMENSION(nw, nflay) :: deltau, delsca, delray, delo3abs, delabs

    !  Exception handling for VLIDORT Model Calculation. New code, 13 October 2010
    LOGICAL :: openfileflag
    INTEGER :: status_inputcheck, status_inputread, status_calculation, & 
             ncheckmessages, nreadmessages
    CHARACTER(LEN=maxchlen), DIMENSION(0:MAX_MESSAGES) :: checkmessages, checkactions, &
                                                        readmessages, readactions
    CHARACTER(LEN=maxchlen) :: message, trace_1, trace_2, trace_3

    LOGICAL, SAVE                                 :: first = .TRUE.
    INTEGER, SAVE                                 :: nz, faer_lvl, npolcorr, nradcal
    REAL (KIND=dp), DIMENSION(0:mflay),                   SAVE :: ts, ps
    !REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),       SAVE :: aersca, aerext, aerasy 
    !REAL (KIND=dp), DIMENSION(max_fit_pts,0:maxmom,maxgksec,mflay), SAVE :: aermoms
    !REAL (KIND=dp), DIMENSION(3, max_fit_pts),            SAVE :: abscrs_qtdepen
    REAL (KIND=dp), DIMENSION(max_fit_pts, mflay,max_psl),SAVE :: dadp
    LOGICAL, DIMENSION(mflay),                            SAVE :: aermsk
    LOGICAL, DIMENSION(max_fit_pts),                      SAVE :: do_radcals, do_polcorrs
    INTEGER, DIMENSION(max_fit_pts),                      SAVE :: polcorr_idxs, radcal_idxs
    REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),        SAVE :: abscrs, dads, dadt
    REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),        SAVE :: so2crs, o4crs,o2crs, h2ocrs

    ! Will become an input parameter later in ozprof.inp
    LOGICAL :: do_ssfullb295 = .FALSE.
    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=15), PARAMETER :: modulename = 'LIDORT_PROF_ENV'

    !do_polwf = .false. ;IF (nfpol > 0) do_polwf = .true.
    errstat = pge_errstat_ok
    status_inputcheck  = vlidort_success
    status_calculation = vlidort_success
    status_inputread   = vlidort_success
    IF (first) THEN
     ! ======================= Read LIDORT Control Input ==========================
      ! FIXME file path below should NOT be hard coded !
      !CALL VLIDORT_L_INPUT_MASTER ( &
      !     TRIM(tabdir)//'../control/vlidort_control.inp', &
      !     'o3prof_lidort_error', status_inputread)
       openfileflag   = .FALSE.
       CALL VLIDORT_L_INPUT_MASTER ( &
         TRIM(tabdir)//'../control/vlidort_control_vv2p4RTC.inp', &
            status_inputread, nreadmessages, readmessages, readactions )

       IF (status_inputread .NE. vlidort_success) THEN
           WRITE(www_lun, *) modulename, &
               ': Problems encountered with VLIDORT input read!!!'
               errstat = pge_errstat_error
           !GO TO 9999
           RETURN     
       ENDIF

       ngreek_moments_input = nmom
       earth_radius = rearth

       IF (.NOT. aerosol .AND. (.NOT. has_clouds .OR. do_lambcld)) THEN
          DO_SSCORR_TRUNCATION = .FALSE.
          DO_DELTAM_SCALING = .FALSE.
          DO_RAYLEIGH_ONLY = .TRUE.
          DO_SOLUTION_SAVING = .FALSE.
          DO_BVP_TELESCOPING = .FALSE.
       ENDIF
       first = .FALSE.
    ENDIF
  ! ============= Overridden some control and atmospheric variables ============== 
  IF (num_iter == 0 ) THEN
     n_szangles = 1
     szangles(1) = sza 
     IF (sza >= 90.0 .OR. sza < 0) THEN
        WRITE(*, *) modulename, ' : SZA is >= 90 or < 0 !!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     n_user_vzangles = 1
     n_user_relazms = 1
     user_vzangles(1) = vza
     user_relazms(1) = aza
     !IF (.NOT. do_multi_vza) THEN
     !ELSE
     !   n_user_streams = 3; n_user_relazms = 3
     !   user_angles_input(1:3) = ABS(vza_atm); user_relazms(1:3) = aza_atm
     !   ! vza will be sorted automatically in LIDORT, need to get order index
     !   vind(1) = MINVAL(MINLOC(ABS(vza_atm)))
     !   vind(3) = MINVAL(MAXLOC(ABS(vza_atm)))
     !   vind(2) = 6 - vind(1) - vind(3)
     !   IF (vind(2) < 1 .OR. vind(2) > 3) THEN
     !      WRITE(*, *) modulename, ' : Sth. wrong in view geometry!!!'
     !      errstat = pge_errstat_error; RETURN
     !   ENDIF
     !
     !   vind((/vind(1), vind(2), vind(3)/)) = (/1, 2, 3/)
     !   naza =3; nvza = 3
     !ENDIF
     nz = nflay;   nlayers = nz

     IF (nz > maxlayers) THEN
        WRITE(*, *) modulename, ' : # of layers exceeded allowed !!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
     IF (nl > nlayers) THEN
        WRITE(*, *) modulename, ' : Coarse grids cannot be finer than fine grids!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     ! Set height grid for doing Chapman Function Calculation 
     height_grid(0:nlayers) = fzs(0:nz)
     geometry_specheight  = the_surfalt

     ! set the first aerosol layer (from TOA down) 1:faer_lvl-1 (without aerosols)
     faer_lvl = 1;   IF (.NOT. strat_aerosol) faer_lvl = nup2p(ntp) + 1

     IF (aerosol) THEN 
        aermsk = .TRUE.
        IF (.NOT. strat_aerosol) aermsk(1:faer_lvl-1) = .FALSE.
     ELSE
        aermsk = .FALSE.
     ENDIF

     ts(1:nz) = (fts(1:nz) + fts(0:nz-1)) / 2.0
     ps(1:nz) = exp((log(fps(1:nz)) + log(fps(0:nz-1)))/2.0)
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

  ! albedo/cloud array
  IF (.NOT. vary_sfcalb) THEN
     DO i = 1, n0alb
        albs(albpmin(i):albpmax(i)) = albarr(i)
     ENDDO
     DO i = 1, n0wfc
        wfcs(wfcpmin(i):wfcpmax(i)) = wfcarr(i)
     ENDDO

  ELSE
     albs(1:nw) = walb0s(1:nw)
     wfcs(1:nw) = wfc0s(1:nw)
  ENDIF


  ! Determinine atmospheric weighting functions to be calculated
  IF (do_ozwf .OR.  do_tmpwf .OR. do_o3shi .OR. do_taodwf .OR. &
       do_twaewf .OR. do_saodwf .OR. do_codwf .OR. do_sprswf )  do_atmos_linearization = .TRUE.
  do_surface_linearization = do_albwf    
  IF (do_atmos_linearization .OR. do_albwf) THEN
     do_simulation_only =    .FALSE.;    do_linearization = .TRUE.
  ELSE 
     do_simulation_only =    .TRUE.;     do_linearization = .FALSE.
  ENDIF
  
  do_o3hwe = .false. ; do_o3spk = .false.
  IF (do_pslwf) THEN 
     IF (np1 > 0 ) do_o3hwe = .true.
     IF (np2 > 0 ) do_o3spk = .true.
  ENDIF

  i = 0;  ozwfidx = 0; aodwfidx = 0; twaewfidx = 0; codwfidx = 0; sprswfidx = 0
  layer_vary_flag(1:nz) = .FALSE.
  layer_vary_number(1:nz) = 0
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi) THEN
     i = i + 1;  ozwfidx = i
     profilewf_names(i) = 'ozone volume mixing ratio------'

     layer_vary_flag(1:nz) = varyprof(1:nz)
     WHERE (varyprof(1:nz) )
        layer_vary_number(1:nz) = layer_vary_number(1:nz) + 1
     ENDWHERE
  ENDIF
  IF ( do_taodwf .OR. do_saodwf) THEN
     i = i + 1; aodwfidx = i
     profilewf_names(i) = 'aerosol extinction coefficient-'

     IF (do_taodwf) THEN
        layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(ntp)+1:nz1) = layer_vary_number(nup2p(ntp)+1:nz1) + 1
     ENDIF

     IF (do_saodwf) THEN
        layer_vary_flag(1:nup2p(ntp)) = .TRUE.
        layer_vary_number(1:nup2p(ntp)) = layer_vary_number(1:nup2p(ntp)) + 1
     ENDIF
  ENDIF
  IF ( do_twaewf) THEN
     i = i + 1; twaewfidx = i
     profilewf_names(i) = 'aerosol scattering coefficient-'
     layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
     layer_vary_number(nup2p(ntp)+1:nz1) = layer_vary_number(nup2p(ntp)+1:nz1) + 1
  ENDIF
  IF ( do_codwf) THEN
     i = i + 1; codwfidx = i
     profilewf_names(i) = 'cloud extinction coefficient---'
     layer_vary_flag(nctp:ncbp) = .TRUE.
     layer_vary_number(nctp:ncbp) = layer_vary_number(nctp:ncbp) + 1
  ENDIF
  do_fraywf = .FALSE.
  IF ( do_sprswf  .OR. (.NOT. use_effcrs .AND. nw > 1) ) THEN
     ! Need to use jacobians wrt rayleigh OD to perform interpolation

     i = i + 1; sprswfidx = i; raywfidx = i
     profilewf_names(i) = 'rayleigh optical thickness-----'
     IF (.NOT. use_effcrs) THEN
        do_fraywf = .TRUE.
        layer_vary_flag(1:nz1) = .TRUE.
        layer_vary_number(1:nz1) = layer_vary_number(1:nz1) + 1
     ELSE
        layer_vary_flag(nup2p(nsfc-1)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(nsfc-1)+1:nz1) = layer_vary_number(nup2p(nsfc-1)+1:nz1) + 1
     ENDIF
  ENDIF

  n_totalatmos_wfs = i
  ! xliu, 03/8/11, has to be initialized in v2p4RTC
  n_surface_wfs = 1   

  IF (n_totalatmos_wfs > 1) THEN
     WHERE(layer_vary_number(1:nz) > 0) 
        layer_vary_number(1:nz) = n_totalatmos_wfs
     ENDWHERE
  ENDIF

  do_fozwf = .FALSE.
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi .OR. (.NOT. use_effcrs .AND. nw > 1)) do_fozwf = .TRUE.
  do_faerwf = .FALSE.
  IF (do_taodwf .OR. do_saodwf) do_faerwf = .TRUE.

  !print *, n_totalatmos_wfs, nz
  !print *, layer_vary_flag(1:nz)
  !print *, layer_vary_number(1:nz)
  nw0 = nw 
  ! ==================== Get Ozone Absorption Cross Section ====================   
  IF (nw > 1) THEN
    allcol(1, 1:nz1) = ozs(1:nz1) * du2mol
    nfgas1 = 1; gasin(1) = 1
    DO k = 1, ngas
      IF (fgasidxs(k) > 0) THEN
        nfgas1 = nfgas1 + 1; gasin(nfgas1) = nfgas1
        ! molecules cm^-2, but normalized by refspec_norm(gasidxs(k))
        ! allcol / refspec_norm will be molecules cm^-2
        !WRITE(*,'(2i3, a4, 10e15.7)')  k, nz+1,refspec_strings(gasidxs(k)), mgasprof(k,nz1),&
        ! tracegas(k,4), mgasprof(k, nz + 1)
        allcol(nfgas1, 1:nz1) = mgasprof(k, 1:nz1) * tracegas(k, 4) / mgasprof(k, nz+1)
      ENDIF
    ENDDO
    IF (use_effcrs) THEN
    ! Get temperature-dependent ozone cross section at instrumental spectral resolution
      dadp = 0.0
      CALL get_effres_gascrs_ray (num_iter, n_radwvl_sav, radwvl_sav, nw,&
          nz1,ts(1:nz1), ps(1:nz1),  nfgas1, allcol(1:nfgas1, 1:nz1), frhos(1:nz1),&
          do_o3shi, o3shi, do_tmpwf, do_o3hwe,do_o3spk, &
          allcrs(1:nw, 1:nfgas1, 1:nz1), &
          dads(1:nw, 1:nz1),   dadt(1:nw, 1:nz1),& 
          dadp(1:nw, 1:nz1,1), dadp(1:nw, 1:nz1, 2), &
          raycof(1:nw), depol(1:nw), problems)
      abscrs(1:nw, 1:nz1) = allcrs(1:nw, 1, 1:nz1)
      
      IF (problems) THEN 
        WRITE(www_lun, *) modulename//':Errors in get_effres_gascrs_ray'
        errstat = pge_errstat_error; RETURN
      ENDIF

      IF (do_radinter) THEN
        IF (num_iter == 0) THEN
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
                       !IF (abscrs(i, nz1) > abscrs(i-1, nz1) .AND. abscrs(i, nz1) > abscrs(i+1, nz1)) &
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
           ENDIF
        ELSE
           IF (num_iter == 0) THEN
              do_radcals(1:nw) = .TRUE.; nradcal = nw
              DO i = 1, nw
                 radcal_idxs(i) = i
              ENDDO
           ENDIF
        ENDIF

    ELSE   ! .NOT. use_effcrs
        ! O3/SO2 (use_so2dtcrs=.TRUE.) cross section: if do_tmpwf = .FALSE. and do_o3shi is false, 
        ! just need to get once for each retrieval
        ! Other trace gas cross section: just need to get it once for all the retrievals if no shifts 
        CALL GET_HRES_GASCRS_RAY(num_iter, nw, waves, nz1, ts(1:nz1),ps(1:nz1),&
             nfgas1, allcol(1:nfgas1, 1:nz1), frhos(1:nz1), &  
             do_o3shi,  o3shi, do_tmpwf,&
             allcrs(1:nw, 1:nfgas1, 1:nz1),raycof(1:nw), depol(1:nw), errstat)
        nw0 = ncalcp

    ENDIF 
    
    abscrs(1:nw0, 1:nz1) = allcrs(1:nw0, 1, 1:nz1)
    IF (use_so2dtcrs) so2crs(1:nw0, 1:nz1) = allcrs(1:nw0, so2crsidx, 1:nz1)
    IF (use_o4dtcrs)  o4crs(1:nw0, 1:nz1)  = allcrs(1:nw0, o4crsidx, 1:nz1)
    IF (use_o2dptcrs) o2crs(1:nw0, 1:nz1)  = allcrs(1:nw0, o2crsidx, 1:nz1)
    IF (use_h2odptcrs) h2ocrs(1:nw0, 1:nz1)= allcrs(1:nw0, h2ocrsidx, 1:nz1)        
        
    IF (num_iter == 0) THEN
        do_radcals(1:nw) = .FALSE.
        nradcal = nw0
        do_radcals(1:nw0) = .TRUE.
        radcal_idxs(1:nw0) = (/(i, i=1,nw0)/)
        IF (use_so2dtcrs) THEN
            fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(so2idx)))= mgasprof(so2idx, nflay+1) * refspec_norm(so2_idx)
        ENDIF
        IF (use_o4dtcrs) THEN 
            fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(o4idx)))  = mgasprof(o4idx, nflay+1) * refspec_norm(o2o2_idx)
        ENDIF
        IF (use_o2dptcrs) THEN 
            fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(o2idx))) = mgasprof(o2idx, nflay+1) * refspec_norm(o2_idx)
        ENDIF
        IF (use_h2odptcrs) THEN 
            fitvar_rad_apriori (mask_fitvar_rad(fgasidxs(h2oidx))) = mgasprof(h2oidx, nflay+1) * refspec_norm(h2o_idx)
        ENDIF
     ENDIF
  ELSE
     do_radcals(1) = .TRUE.; nradcal = 1; radcal_idxs(1) = 1

     ! o3 absorption coefficient at 370.2 nm with TOMS FWHM
     dads(1, 1:nz1) = 0.0; dadt(1, 1:nz1) = 0.0

     ! Weighted by solar flux
     !abscrs(1, 1:nz) = 9.1231787D-24 + (ts(1:nz) - 273.15) * 1.9005502D-25 + &
     !     (ts(1:nz) - 273.15)**2.0 * 1.2275286D-27     
     !raycof(1) = 2.3184501D-26; depol(1) = 0.030247913D0

     nfgas1 = 7
     CALL GET_ALB_OZCRS_RAY(nz1, ts(1:nz1), nfgas1, allcrs(1, 1:nfgas1, 1:nz1), raycof(1), depol(1), problems)
     IF (problems) THEN
        WRITE(*, *) modulename, ' : Problems in reading O3 XSec for determining Fc!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     !CALL GET_ALL_RAYCOF_DEPOL(1, 360.D0, raycof(1), depol(1))
     !print *, raycof(1), depol(1)

     ! Ozone
     nfgas1 = 6;  gasin(1) = 1
     allcol(1, 1:nz1) = ozs(1:nz1) * du2mol
   
     tmp_gaspos =([no2idx,so2idx,broidx, hchoidx,o4idx])
     DO i = 2, 5
        allcol(i, 1:nz1) = mgasprof(tmp_gaspos(i-1), 1:nz1)
        gasin(i) = i
     ENDDO
    
     ! O2-O2 Conc !NO2, SO2, BRO, HCHO, O4 (Oxygen: 20.95%)
     allcol(6, 1:nz1) = ( frhos(1:nz1) * 0.2095 ) ** 2.0 / (fzs(0:nz1-1) - fzs(1:nz1)) / 1.0D5
     gasin(6) = 6
  ENDIF

  WRITE(www_lun,*) modulename//'==> Finishing get allcol, allcrs'
  ! Initialize aerosol/cloud property profiles
  IF (num_iter == 0) THEN
     aersca(1:nz1) = 0.0; aerext(1:nz1) = 0.0
     aerasy(1:nz1) = 0.0;  aermoms(0:nmom, 1:ngksec, 1:nz1) = 0.0
  ENDIF

  cldmsk = .FALSE.; IF (nctp /= 0) cldmsk(nctp:ncbp)=.TRUE.
  cldsca=0.0; cldext=0.0; cldasy=0.0;  cldmoms = 0.0

  ! Initialize output variables  
  rad = 0.0; radclr = 0.0; radcld = 0.0
  cfracwf = 0.0
  IF ( do_fozwf )  fozwf   = 0.0
  IF ( do_albwf)   albwf   = 0.0
  IF ( do_codwf )  fcodwf  = 0.0
  IF ( do_sprswf ) fsprswf  = 0.0
  IF ( do_fraywf ) fraywf  = 0.0
  IF ( do_ctpwf )  ctpwf   = 0.0
  IF ( do_faerwf)  faerwf  = 0.0
  IF ( do_twaewf ) faerswf = 0.0
  IF ( do_pslwf) database_pslwf = 0.0

  ! senstitivity of absorption cross section to temperature, used 
  ! for calculating temperature wf directly with LIDORT
  !eta = 0.0             ! dummy variable here
  alleta = 0.0

  polerr = 1.0
  IF ( polcorr == 2) THEN 
     IF (nwfc > 0) THEN
        STOP 'Polarization correction + variable fc: not implemented!!!'
     ENDIF
     tmpalbs = albs 
     ! clear-sky
     IF (the_cfrac /= 1.0) THEN
        toz = SUM(ozs(1:nz1))

        WHERE(tmpalbs >= 1.0)
           tmpalbs = 0.999
        ENDWHERE

        WHERE(tmpalbs < 0.0)
           tmpalbs = 0.001
        ENDWHERE

        !        CALL lup_polerror(toz, fps(nz1), nw, waves, tmpalbs, polerr(1:nw, 1, 1), errstat)
        IF (errstat == pge_errstat_error) RETURN
        polerr(1:nw, 1, 1) = 1.0 + polerr(1:nw, 1, 1) / 100.0
     ENDIF

     ! cloudy part
     IF (the_cfrac /= 0.0) THEN
        toz = SUM(ozs(1:nctp-1))
        IF (do_lambcld) THEN
           tmpalbs = lambcld_refl
        ELSE
           ! approxmation of cod as cloud albedo for polarization correction
           ! Follows Kokhanovsky et al., 2003
           tmpalbs = 1.0 - 1.0 / (1.072 + 0.1125 * the_cod)
        ENDIF

        WHERE(tmpalbs >= 1.0)
           tmpalbs = 0.999
        ENDWHERE

        WHERE(tmpalbs < 0.0)
           tmpalbs = 0.001
        ENDWHERE
        !        CALL lup_polerror(toz, fps(nctp-1), nw, waves, tmpalbs, polerr(1:nw, 2, 1), errstat)
        IF (errstat == pge_errstat_error) RETURN
        polerr(1:nw, 2, 1) = 1.0 + polerr(1:nw, 2, 1) / 100.0
     ENDIF
  ENDIF

  ! Determine wavelengths where exact polarization correction (NSTOKES: 4 vs 1) is calculated
  ! In UV1 (or between 270 and 310 nm): ~292 nm, ~298 nm, ~300 nm, ~302 nm, ~304 nm, ~306 nm, last wavelength 
  ! In UV2 (or between 310 and 340 nm): first, 1/4, middle and last wavelength
  ! So exact vector LIDORT calculation is done at 11 wavelengths.
  ! This option works when radiance interpolation option is turned on
  IF (num_iter == 0) THEN
     do_polcorrs(1:nw) = .FALSE.; npolcorr = 0
     IF ( (polcorr >= 3 .AND. polcorr <= 5) .AND. nw > 1 ) THEN
        fidx = 1; idum = 0
        DO iw = 1, numwin
           IF (iw == numwin) THEN
              temp = winlim(iw, 2)
           ELSE
              temp = (winlim(iw, 2) + winlim(iw + 1, 1)) / 2.
           ENDIF
           lidx = MINVAL(MAXLOC(waves(1:nw), MASK=(waves(1:nw) < temp .AND. waves(1:nw) > 0))) 
           !lidx = fidx + nradpix(iw) - 1
           IF ( waves(lidx) <= 312.0 ) THEN
              ! Error from using single scattering is about 0.2% at 270 nm, it needs o be corrected
              !IF (do_ssfullb295) do_polcorrs(fidx) = .TRUE. 
              IF (idum == 0) idum = 1
              DO i = fidx + 1, lidx - 1
                 IF (do_radcals(i) ) THEN
                    IF ( (waves(idum) < 290. .AND. waves(i) >= 290.0) .OR.  &
                         (waves(idum) < 295. .AND. waves(i) >= 295.0) .OR.  &
                         (waves(idum) < 299. .AND. waves(i) >= 299.0) .OR.  &
                         (waves(idum) < 301. .AND. waves(i) >= 301.0) .OR.  &
                         (waves(idum) < 303. .AND. waves(i) >= 303.) .OR.  &
                                !(waves(idum) < 304. .AND. waves(i) >= 304.0) .OR.  &
                         (waves(idum) < 305. .AND. waves(i) >= 305.0) )  do_polcorrs(i) = .TRUE.
                    idum = i
                 ENDIF
              ENDDO
              do_polcorrs(lidx) = .TRUE.
           ELSE
              IF (idum == 0) idum = 1
              do_polcorrs(fidx) = .TRUE.; do_polcorrs(lidx) = .TRUE.     
              j = (fidx + lidx ) / 2; k = fidx + (j - fidx) / 3  

              jj = (j + lidx )  / 2; kk = fidx + (k - fidx) / 3
              jk = (j + k ) / 2

              DO i = fidx + 1, lidx - 1
                 IF (do_radcals(i) ) THEN
                    !IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                    !     (waves(idum) < waves(k) .AND. waves(i) >= waves(k)) .OR. &
                    !     (waves(idum) < waves(jj) .AND. waves(i) >= waves(jj)) .OR. &
                    !     (waves(idum) < waves(kk) .AND. waves(i) >= waves(kk)) .OR. &
                    !     (waves(idum) < waves(jk) .AND. waves(i) >= waves(jk)))  do_polcorrs(i) = .TRUE.
                    IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                         (waves(idum) < waves(k) .AND. waves(i) >= waves(k)))  do_polcorrs(i) = .TRUE.
                    idum = i
                 ENDIF
              ENDDO
           ENDIF

           fidx = lidx + 1
        ENDDO
        DO i = 1, nw
           IF ( do_polcorrs(i) ) THEN
              !IF (do_ssfullb295 .AND. npolcorr == 2 .AND. waves(1) < 290.0) THEN
              !   ! Use single scattering (ss), scalar below 295 nm
              !   ! Use ss + multiple scattering, scalar above 295 nm
              !   ! Correction: 270 (ss + scalr), 290  (ss + scalar),
              !   !             first lamda before 295 (ss + scalar)
              !   !             first lamda after  295 (scalar)
              !   do_polcorrs(i-1) = .TRUE.; npolcorr = npolcorr + 2
              !   polcorr_idxs(3) = i-1; polcorr_idxs(4) = i
              !ELSE
              npolcorr = npolcorr + 1; polcorr_idxs(npolcorr) = i
              !ENDIF
           ENDIF
        ENDDO

        !print *, nradcal, npolcorr
        !print *, polcorr_idxs(1:npolcorr)
        !print *, do_radcals(polcorr_idxs(1:npolcorr))
        !print *, waves(polcorr_idxs(1:npolcorr))
        !STOP
     ENDIF
  ENDIF
  ! ====================== Call LIDORT and Do Post Processing =====================
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
     IF (polcorr == 0 .OR. ((polcorr == 3 .OR. polcorr == 5) .AND. (nw == 1 .OR. do_polcorrs(iw)))  &
          .OR. (polcorr == 4 .AND. (nw == 1 .OR. (do_polcorrs(iw) .AND. num_iter == 0))) ) THEN
        NSTOKES = 3 ; NSTREAMS = 4
     ELSE
        NSTOKES = 1 ; NSTREAMS = 4
     ENDIF
     VlidortNstream = NSTREAMS
     
     IF ( lamda < 295.0 .AND. NSTOKES == 1 .AND. nw > 1 .AND. do_ssfullb295) THEN
        DO_SSFULL = .TRUE.;  DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
     ELSE
        DO_SSFULL = .FALSE.
        IF (aerosol .OR. (has_clouds .AND. .NOT. do_lambcld)) THEN
           DO_SSCORR_TRUNCATION = .TRUE.; DO_DELTAM_SCALING = .TRUE.
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
                 cldmoms(i, j, nctp:ncbp) = gcmoms(low, i, j) &
                      * (1.0 - xg) + gcmoms(hgh, i, j) * xg 
              ENDDO
           ENDDO
        ENDIF
     ENDIF

     IF ( do_polcorrs(iw) .AND. &
         ((polcorr == 3 .OR. polcorr == 5) .OR. (polcorr == 4 .AND. &
                                                 (num_iter == 0 .or. num_iter == 2 ) )) ) THEN
        npolmod = 2   ! Twice, one vector and one scalar
     ELSE 
        npolmod = 1   ! Only once either scalar or vector
     ENDIF

     ! When polcorr = 5 is selected, only calculate weighting function for iteration 1
     ! iteration 0 if 1st pixel of a x-track position is being retrieved
     IF (polcorr == 5 .AND. nw > 1 .AND. (num_iter > 1 .OR. (num_iter == 0 .AND. currloop /= 0))) THEN
        do_simulation_only = .TRUE.
        do_linearization   = .FALSE.
     ENDIF

     DO ipol = 1, npolmod
        IF (ipol == 2) THEN
           NSTOKES = 1 ; NSTREAMS = 4 ! Always scalar for second mode

           IF ( lamda < 295.0 .AND. nw > 1 .AND. do_ssfullb295) THEN
              DO_SSFULL = .TRUE.;  DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
           ELSE
              DO_SSFULL = .FALSE.; DO_SSCORR_TRUNCATION = .TRUE.; DO_DELTAM_SCALING = .TRUE.
              !DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
           ENDIF
           IF (.NOT. aerosol) THEN
              DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
           ENDIF
              
           ! Save VECTOR LIDORT results
           prad(iw) = rad(iw, 1)
           rad(iw, 1:nostk) = 0.d0
           IF (do_cfracwf) pcfracwf(iw) = cfracwf(iw, 1)
           IF (do_cfracwf) cfracwf(iw, 1:nostk) =  0.d0

           IF ( do_linearization ) THEN
              IF (do_albwf)  palbwf(iw)             = albwf(iw, 1)
              IF (do_fozwf)  pfozwf(iw, 1:nz)          = fozwf(iw, 1:nz, 1)
              IF (do_codwf)  pfcodwf(iw, nctp:ncbp) = fcodwf(iw, nctp:ncbp, 1)
              IF (do_sprswf) pfsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1) = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1)
              IF (do_fraywf) pfraywf(iw, 1:nz)         = fraywf(iw, 1:nz, 1)
              IF (do_ctpwf)  pctpwf(iw)             = ctpwf(iw, 1)
              IF (do_faerwf) pfaerwf(iw, faer_lvl:nfsfc-1)  = faerwf(iw, faer_lvl:nfsfc-1, 1)           
              IF (do_twaewf) pfaerswf(iw, faer_lvl:nfsfc-1) = faerwf(iw, faer_lvl:nfsfc-1, 1)

              ! Initalize those variables to zero again          
              IF (do_albwf) albwf(iw, 1:nostk) = 0.d0
              IF (do_fozwf) fozwf(iw, :, 1:nostk) = 0.d0 
              IF (do_codwf) fcodwf(iw, nctp:ncbp, 1:nostk) = 0.d0
              IF (do_sprswf) fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1:nostk) = 0.d0
              IF (do_fraywf) fraywf(iw, :, 1:nostk) = 0.d0
              IF (do_ctpwf)  ctpwf(iw, 1:nostk)       = 0.d0
              IF (do_faerwf) faerwf(iw, faer_lvl:nfsfc-1, 1:nostk)  = 0.d0
              IF (do_twaewf) faerswf(iw, faer_lvl:nfsfc-1, 1:nostk) = 0.d0

           ENDIF
        ENDIF

        !print *, iw, ipol, do_simulation_only, do_linearization, do_atmos_linearization, do_surface_linearization

        radclrcld = 0.0
        DO ic = 1, 2               ! for clear and cloud

           IF (ic == 1) THEN
              do_clouds = .FALSE.; frac = 1.0 - the_cfrac
           ELSE
              do_clouds = .TRUE. ; frac = the_cfrac
           ENDIF
           IF (frac == 0.0) CYCLE  ! No clear/cloudy part

           ! Note convert ozone from DU to molecule/cm^2  here
           IF ((ic == 1) .OR. (.NOT. do_lambcld))  THEN
              lambertian_albedo = albs(iw)

              ! Reset up number of layers since Lambertian 
              ! cloudy scene got different layers
              nlayers = nfsfc - 1; nz1 = nlayers   

              ! zero O3 weighting function below surface (This is necessary)           
              profilewf(1:n_totalatmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0

              IF (ipol == 1) THEN
                 !WRITE(91, '(4F10.4,3I5)') sza, vza, aza, lambertian_albedo, 1, nlayers, ngreek_moments_input
                 !WRITE(91, '(40D14.6)') height_grid(0:nlayers)
                 !WRITE(91, '(2D14.6)') lamda, depol(iw)
                 CALL LIDORT_PROF_PREP(lamda, raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                      varyprof(1:nz1), nfgas1, gasin(1:nfgas1), allcrs(iw, 1:nfgas1, 1:nz1), allcol(1:nfgas1, 1:nz1), &
                      alleta(1:nfgas1, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
                      aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                      do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                      cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems, &
                      deltau(iw, 1:nz1), delsca(iw, 1:nz1), delo3abs(iw, 1:nz1), delray(iw, 1:nz1))
                 IF (problems) THEN
                    WRITE(www_lun, *) modulename, ' : Problems encountered in lidort preparation!!!' 
                    errstat = pge_errstat_error; RETURN
                 END IF
              ENDIF
           ELSE                 ! lambertian clouds
              nlayers = nctp-1  ! from cloud top to TOA
              nz1 = nlayers;    do_clouds = .FALSE.

              IF (the_cfrac == 1.0 .AND. nw /= 1) THEN
                 lambertian_albedo = albs(iw); lambcld_refl = lambertian_albedo
              ELSE
                 lambertian_albedo = lambcld_refl ! use 80% (could be adjusted when the_cfrac gt 0.90)
              ENDIF

              ! zero O3 weighting function (This is necessary)           
              profilewf(1:n_totalatmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0

              IF (frac == 1.0 .AND. ipol == 1) THEN
                 CALL LIDORT_PROF_PREP(lamda, raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                      varyprof(1:nz1), nfgas1, gasin(1:nfgas1), allcrs(iw, 1:nfgas1, 1:nz1), allcol(1:nfgas1, 1:nz1), &
                      alleta(1:nfgas1, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
                      aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                      do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                      cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems , &
                      deltau(iw, 1:nz1), delsca(iw, 1:nz1), delo3abs(iw, 1:nz1), delray(iw, 1:nz1))

                 IF (problems) THEN
                    WRITE(*, *) modulename, ' : Problems encountered in lidort preparation!!!'
                    errstat = pge_errstat_error; RETURN
                 END IF
              ENDIF
           ENDIF
         
           ! No Need to reset viewing geometry for each LIDORT call
           !geometry_specheight = height_grid(nlayers)
           !szangles(1) = sza; user_vzangles(1) = vza; user_relazms(1) = aza  

           !xliu: 03/07/2011, switch VLIDORT version (from vv2p4 to vv2p4RTC)
           ! Used in v2p4
           !CALL VLIDORT_L_MASTER_LAMBERTIAN (status_inputcheck, status_calculation) 
           !CALL VLIDORT_STATUS ( status_inputcheck, status_calculation )

           ! used in v2p4RTC
           CALL VLIDORT_L_MASTER ( status_inputcheck, ncheckmessages, checkmessages, &
                checkactions, status_calculation, message, trace_1, trace_2, trace_3)
           CALL VLIDORT_STATUS ( 'o3prof_lidort_error', vlidort_errunit, openfileflag, &
                status_inputcheck,  ncheckmessages, checkmessages, checkactions, &
                status_calculation, message, trace_1, trace_2, trace_3 )

           IF (status_inputcheck   /= vlidort_success) THEN
              WRITE(*, *) modulename, ' : Problems encountered in lidort input check!!!'
              errstat = pge_errstat_error; RETURN
           ENDIF
           IF (status_calculation  /= vlidort_success) THEN
              WRITE(*, *) modulename, ' : Problems encountered in lidort calculation!!!'
              errstat = pge_errstat_error; RETURN
           ENDIF

           ! Pixel-independent approximation
           radclrcld(ic, 1:nostk) = stokes(1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk)
           rad(iw, 1:nostk)       = rad(iw, 1:nostk) + radclrcld(ic, 1:nostk) * frac

!           print * ,iw, lamda,  deltau(iw, nz1), delsca(iw, nz1), delo3abs(iw, nz1), rad(iw, 1)
           IF (do_linearization ) THEN

              ! weighting function per Dobson Unit
              IF ( do_ozwf ) THEN
                 DO istk = 1, nostk
                    fozwf(iw, 1:nfsfc-1, istk) = fozwf(iw, 1:nfsfc-1, istk) + profilewf(ozwfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                         / ozs(1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
                 ENDDO
              ENDIF

              ! Weighting function with respect to aerosol/cloud optical depth at the last wavelength
              ! so as to keep the scaling since we are fitting the aod at that wavelength
              IF ( do_taodwf .OR. do_saodwf ) THEN
                 !print *, ipol, faer_lvl, nfsfc-1, frac
                 DO istk = 1, nostk
                    faerwf(iw, faer_lvl:nfsfc-1, istk)  = faerwf(iw, faer_lvl:nfsfc-1, istk) + &
                         profilewf(aodwfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                         gaext(actawin, faer_lvl:nfsfc-1) * polerr(iw, ic, istk) *  frac
                    !print *, faerwf(iw, faer_lvl:nfsfc-1, 1)
                 ENDDO
              ENDIF
              IF ( do_twaewf ) THEN
                 DO istk = 1, nostk
                    faerswf(iw, faer_lvl:nfsfc-1, istk)  = faerswf(iw, faer_lvl:nfsfc-1, istk) + &
                         profilewf(twaewfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                         gasca(actawin, faer_lvl:nfsfc-1) * gaext(actawin, faer_lvl:nfsfc-1) * &
                         polerr(iw, ic, istk) *  frac
                 ENDDO
              ENDIF
              IF ( do_codwf ) THEN
                 DO istk = 1, nostk
                    fcodwf(iw, nctp:ncbp, istk)  = fcodwf(iw, nctp:ncbp, istk) + &
                         profilewf(codwfidx, nctp:ncbp, 1, 1, istk, 1) / &
                         cldext0(nctp:ncbp) * polerr(iw, ic, istk) *  frac
                 ENDDO
              ENDIF
              IF ( do_sprswf ) THEN
                 DO istk = 1, nostk 
                    fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk)  = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk) + &
                         profilewf(sprswfidx, nup2p(nsfc-1)+1:nfsfc-1, 1, 1, istk, 1) / &
                         (fps(nup2p(nsfc-1)+1:nfsfc-1) - fps(nup2p(nsfc-1):nfsfc-2)) * polerr(iw, ic, istk) *  frac
                 ENDDO
              ENDIF

              IF ( do_fraywf ) THEN
                 DO istk = 1, nostk
                    fraywf(iw, 1:nfsfc-1, istk) = fraywf(iw, 1:nfsfc-1, istk) + profilewf(raywfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                         / delray(iw, 1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
                 ENDDO
              ENDIF

              ! Non-lambertian clouds: albedo wf is from both clear and cloud
              ! Lamberitan clouds:     if the_cfrac  < 1, albwf from clear only
              !                        if the_cfrac == 1, albwf from cloud only


              ! xliu, 03/08/11: the surface albedo weighting functions from v2p4RTC are un-normalized.
              ! Also, we need to initialize n_surface_wfs = 1 for calculating lambertian surface weighting function
              IF (do_albwf) THEN 
                 IF (.NOT. do_lambcld) THEN
                    albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + surfacewf(1, 1, 1, 1:nostk, 1) &
                         * polerr(iw, ic, 1:nostk) * frac !/ lambertian_albedo
                 ELSE IF (the_cfrac < 1.0) THEN
                    IF (ic == 1) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
                         surfacewf(1, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac ! / lambertian_albedo 
                 ELSE IF (the_cfrac == 1.0) THEN
                    IF (ic == 2) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
                         surfacewf(1, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac !/ lambertian_albedo
                 ENDIF
              ENDIF
              !IF (do_albwf) THEN 
              !  IF (.NOT. do_lambcld) THEN
              !    albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
              !       surfacewf(1, 1, 1, 1:nostk, 1) &
              !       * polerr(iw, ic, 1:nostk) * frac / lambertian_albedo
              !  ELSE IF (the_cfrac < 1.0) THEN
              !  IF (ic == 1) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
              !       surfacewf(1, 1, 1, 1:nostk, 1) * &
              !       polerr(iw, ic, 1:nostk) * frac / lambertian_albedo
              !  ELSE IF (the_cfrac == 1.0) THEN
              !  IF (ic == 2) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
              !       surfacewf(1, 1, 1, 1:nostk, 1) * &
              !       polerr(iw, ic, 1:nostk) * frac / lambertian_albedo
              !  ENDIF
             !ENDIF
           ENDIF

        ENDDO ! end clear/cloudy scene loop
        radclr(iw, 1:nostk) = radclrcld(1, 1:nostk); radcld(iw, 1:nostk) = radclrcld(2, 1:nostk)
        !write(*,'(i4, f8.2, 2e15.5)') i, lamda, radclr(iw,1), radcld(iw,1)
        IF (do_cfracwf) cfracwf(iw, 1:nostk) = radcld(iw, 1:nostk) - radclr(iw, 1:nostk)
     ENDDO    ! end scalar/vector modes
  ENDDO       ! end wavelength loop 
  ! On-line polarization correction
  IF ( the_cfrac == 1.0 .AND. do_lambcld) THEN
     nz1 = nctp - 1
  ELSE
     nz1 = nfsfc - 1
  ENDIF
  IF (rtmdbg) THEN
     IF ( nw ==1 .and. num_iter == 0) THEN 
        WRITE(rtm_unit,*) 'nlayer=', nlayers
        WRITE(rtm_unit,*) 'i fts fps fozs frhos==================================='
        DO i = 1, nflay
           WRITE(rtm_unit,'(i3, 3f8.2, e15.5)') i,fts(i), fps(i),fozs(i), frhos(i)
        ENDDO
        WRITE(rtm_unit,*)
     ENDIF 

     WRITE(rtm_unit, '(A,i3, A, i4,A)') '# iter=', num_iter, '++++++    nw=',nw0,''
     IF (nw==1) THEN
        DO i = 2, nfgas1
         WRITE(rtm_unit,'(i3,a5,500e15.5)') i, &
           refspec_strings(gasidxs(tmp_gaspos(i-1))) , allcrs(1:nw, i, 10), allcol(i, 10)
        ENDDO
     ELSE
       WRITE(rtm_unit, '(a6,100a12)') 'Lamda', &
         'O3',refspec_strings(gasidxs(fgaspos(1:nfgas1-1))), &
         'o3',refspec_strings(gasidxs(fgaspos(1:nfgas1-1)))
       DO i = 1, nw0 , 20
         k = radcal_idxs(i)
         WRITE(rtm_unit,'(f6.2, 500e12.3)')  waves(k), allcrs(k,1:nfgas1,10), allcol(1:nfgas1, 10)
       ENDDO
     ENDIF
      WRITE(rtm_unit, '(A, 2L)') 'waves radclr radcld cfracwf albwf, fozwf',do_cfracwf, do_albwf
     DO i = 1, nw0, 20 
       k = radcal_idxs(i)
       WRITE(rtm_unit, '(f8.2,500e15.5)') waves(k),radclr(k, 1), radcld(k, 1),cfracwf(k,1), albwf(k,1),fozwf(k,10,1)
     ENDDO
  ENDIF

  IF ( (polcorr >= 3 .AND. polcorr <= 5 ) .AND. nw > 1) THEN
     IF (use_effcrs) THEN
        nw0 = nw
     ELSE
        nw0 = ncalcp
     ENDIF
     nsprs = nup2p(nsfc-1)+1
     delabs(1:nw0, 1:nz1) = deltau(1:nw0, 1:nz1) - delsca(1:nw0, 1:nz1)
     !xliu, 11/02/2011, add abscrs in the variables
     CALL polcorr_online(num_iter, polcorr, nw0, nz1, nctp, ncbp, nsprs, & 
          faer_lvl,  npolcorr, polcorr_idxs(1:npolcorr), &
          do_fozwf, do_albwf, do_faerwf,do_twaewf, do_codwf, do_sprswf, do_fraywf, do_cfracwf, &
          waves(1:nw0), rad(1:nw0, 1),  prad(1:nw0), delabs(1:nw0, 1:nz1), abscrs(1:nw0, 1:nz1), ozs(1:nz1), &
          albwf(1:nw0, 1), palbwf(1:nw0), fozwf(1:nw0, 1:nz1, 1), pfozwf(1:nw0, 1:nz1), &
          faerwf(1:nw0, 1:nz1, 1),   pfaerwf(1:nw0, 1:nz1), faerswf(1:nw0, 1:nz1, 1), pfaerswf(1:nw0, 1:nz1),&
          fcodwf(1:nw0, 1:nz1, 1), pfcodwf(1:nw0, 1:nz1), fsprswf(1:nw0, 1:nz1, 1),pfsprswf(1:nw0, 1:nz1), &
          fraywf(1:nw0, 1:nz1, 1), pfraywf(1:nw0, 1:nz1), cfracwf(1:nw0, 1), pcfracwf(1:nw0))

  ENDIF
 ! IF (nw > 1) then 
    !print * ,  nw0,allcrs(nw0,1:nfgas1, 10), allcol(1:nfgas1, 10)
    !print * , allcrs(1:nw0, 3, 10)
 ! ENDIF
  ! Radiance Interpolation
  IF (nw > 1 .AND. do_radinter ) THEN
     nsprs = nup2p(nsfc-1)+1
     CALL radwf_interpol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, do_radcals(1:nw),             &
          do_fozwf, do_albwf, do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf,         &
          waves, abscrs(1:nw, 1:nz1), ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1),       &
          albwf(1:nw, 1), cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
          fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), errstat)
     IF (errstat == pge_errstat_error) RETURN
  ENDIF
  
  IF (nw > 1 .AND. .NOT. use_effcrs) THEN
     IF (do_fraywf .AND. do_fozwf) THEN 
        nsprs = nup2p(nsfc-1)+1
        CALL hres_radwf_inter_convol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, do_albwf, &
             do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf, do_tracewf, &
             do_o3shi, do_tmpwf,do_pslwf, waves, ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1), &
             albwf(1:nw, 1), cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
             fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), fraywf(1:nw, 1:nz1, 1), &
             dads(1:nw, 1:nz1), dadt(1:nw, 1:nz1), abscrs(1:nw, 1:nz1), &
             so2crs(1:nw, 1:nz1), o4crs(1:nw, 1:nz1), o2crs(1:nw, 1:nz1), h2ocrs(1:nw, 1:nz1), errstat)
        IF (errstat == pge_errstat_error) RETURN
     ELSE
        WRITE(*, *) 'Must have O3/rayleigh weighting function to do radiance interpolation!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ENDIF

  !WRITE(71, *) nw, nflay
  !DO iw = 1, nw
  !   IF (do_radcals(iw) ) THEN
  !      write(71, *) '1 '
  !   ELSE
  !      write(71, *) '0 '
  !   ENDIF
  !ENDDO
  !    
  !DO iw = 1, nw 
  !   WRITE(71, '(1000D16.7)') waves(iw), rad(iw, 1), fozwf(iw, 1:nflay, 1), albwf(iw, 1), &
  !        deltau(iw, 1:nflay), delsca(iw, 1:nflay), delo3abs(iw, 1:nflay)
  !ENDDO
  !WRITE(71, '(10D16.7)') ozs(1:nflay)
  !print *, fozwf(100, 1:nflay, 1)
  !
  !STOP
 
   nw0=nw

  ! Calculate desired weighting functions at the end after applying all the correction
  IF ( do_ozwf ) THEN  
     DO i = 1, nl
        fidx = nup2p(i - 1) + 1; lidx = nup2p(i)         
        DO iw = 1, nw0
           DO istk = 1, nostk
              ozwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx)) / SUM(ozs(fidx:lidx))
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
                   * dadt(iw, fidx:lidx) * (fzs(fidx-1:lidx-1)-fzs(fidx:lidx))) / (fzs(fidx-1) - fzs(lidx))
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  IF (do_o3shi) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           o3shiwf(iw, istk) = SUM(fozwf(iw, 1:nz1, istk) * ozs(1:nz1) * dads(iw, 1:nz1))
        ENDDO
     ENDDO
  ENDIF

  IF (do_pslwf .and. use_effcrs) THEN 
    DO k = 1, npsl
      DO iw = 1, nw0
          database_pslwf(k, refidx(iw)) = &
                     SUM(fozwf(iw, 1:nz1, 1) * ozs(1:nz1) * dadp(iw, 1:nz1, k))   
          database_pslwf(k, refidx(iw)) = database_pslwf(k,refidx(iw) ) /rad(iw,1) 
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
  !IF (nw > 1) THEN
  !   DO iw = 1, nw
  !      WRITE(78, '(4D16.7)') waves(iw), albwf(iw, 1), sprswf(iw, 1), cfracwf(iw, 1)
  !   ENDDO
  !ENDIF

  IF (nw > 1 .AND. do_ozwf .AND. do_tracewf ) THEN 
      CALL get_tracegas_wf  (nw, nz, nz1, rad(1:nw,1), &
                             fozwf(1:nw, 1:nz, 1), abscrs(1:nw, 1:nz), &
                             use_so2dtcrs, so2crs(1:nw, 1:nz), &
                             use_o4dtcrs,  o4crs(1:nw,  1:nz), &
                             use_o2dptcrs, o2crs(1:nw,  1:nz), &
                             use_h2odptcrs,h2ocrs(1:nw, 1:nz), &
                             do_so2zwf,    so2zwf(1:nw, 1))
  ENDIF

  IF (nw > 1 .AND. do_simu .AND. .NOT. radcalwrt) THEN
     WRITE(78, *) 2, nz1
     WRITE(78, *) 'Profile: Z, P, T, O3, SO2'
     DO i = 1, nz1
        WRITE(78, '(F10.4, 4D16.7)') fzs(i), fps(i), fts(i), fozs(i), allcol(3, i)/refspec_norm(9)
     ENDDO
     WRITE(78, '(A)') 'Wavelength, radiance, albedo wf, ozone wf, aerosol wf, ozcrs, so2crs'
     DO i = 1, nw !160, 369, 209
        WRITE(78, '(F10.4, 500D16.7)') waves(i), rad(i, 1), albwf(i, 1), fozwf(i, 1:nz1, 1), & !, faerwf(i, 1:nz1, 1), &
             allcrs(i, 1, 1:nz1), allcrs(i, 3, nz1)*refspec_norm(9)
     ENDDO
     STOP
  ENDIF
     !IF (wrtozcrs) THEN
     !   WRITE(91, *) nw
     !   DO iw = 1, nw 
     !      WRITE(91, '(F10.4, 6D14.6)') waves(iw), abscrs_qtdepen(1:3, iw), raycof(iw), depol(iw) !,  &
     !           !database(6, refidx(iw))*refspec_norm(6), database(9, refidx(iw))*refspec_norm(9), &
     !           !database(10, refidx(iw))*refspec_norm(10), database(12, refidx(iw))*refspec_norm(12)
     !   ENDDO
     !   
     !   WRITE(91, *) nz1
     !   WRITE(91, '(2D14.6)') fps(0), fzs(0)
     !   DO i = 1, nz1 
     !      WRITE(91, '(10D14.6)') fps(i), fzs(i), ts(i), frhos(i), allcol(1, i), allcol(2, i)/refspec_norm(6), &
     !           allcol(3, i)/refspec_norm(9), allcol(4, i)/refspec_norm(10), allcol(5, i)/refspec_norm(12) 
     !   ENDDO
     !ENDIF

  RETURN

    !9999 CONTINUE 
    !  OPEN(vlidort_errunit, file='o3prof_lidort_input_error.log', status = 'unknown')
    !  WRITE(vlidort_errunit, *)' FATAL:   Wrong input from VLIDORT input file-read'
    !  WRITE(vlidort_errunit, *)'  ------ Here are the messages and actions '
    !  WRITE(vlidort_errunit,'(A,I3)')'    ** Number of messages = ', nreadmessages
    !  DO i = 1, nreadmessages
    !       nf = LEN(readmessages(i))
    !       na = LEN(readactions(i))
    !       WRITE(vlidort_errunit,'(A,I3,A,A)')'Message # ', i,' : ',readmessages(i)(1:nf)
    !       WRITE(vlidort_errunit,'(A,I3,A,A)')'Action  # ', i,' : ',readactions(i)(1:na)
    !  ENDDO
    !  CLOSE(vlidort_errunit)
    ! ENDIF
 
  END SUBROUTINE lidort_prof_env

  SUBROUTINE LIDORT_PROF_PREP (lamda, raycof, depol, zsgrid, airgrid,  varyprof, &
     ngas, gasin, abscrs, gascol, eta, useasy, nmoms, &          
     do_aerosols, aersca, aerext, aerasy, aermoms, aermsk, &
     do_clouds, cldsca, cldext, cldasy, cldmoms, cldmsk, problems, &
     deltau, delsca, delo3abs, delray)

  USE OMSAO_precision_module
  USE ozprof_data_module, ONLY : maxgksec, maxgkmatc, ngksec, ngkmatc
  USE OMSAO_errstat_module
  !USE m_vlidort90_include
  IMPLICIT NONE  
   
    !===============================  Define Variables ========================
    ! Include files of dimensions and numbers
    INCLUDE 'VLIDORT.PARS'
    ! Include files of input variables
    INCLUDE 'VLIDORT_INPUTS.VARS'
    INCLUDE 'VLIDORT_SETUPS.VARS'
    !  INCLUDE 'VLIDORT_REFLECTANCE.VARS'
    INCLUDE 'VLIDORT_L_INPUTS.VARS'
    INCLUDE 'VLIDORT_BOOKKEEP.VARS'

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
  REAL (KIND=dp), DIMENSION(ngas, nlayers), INTENT(IN)  :: abscrs, gascol, eta

  ! Optional output
  REAL (KIND=dp), DIMENSION(nlayers), INTENT(OUT) :: deltau, delsca, delo3abs, delray
  
  
  ! Output variables
  LOGICAL, INTENT(OUT)   :: problems

  ! Modified variables
  LOGICAL, INTENT(INOUT) :: do_aerosols, do_clouds

  ! Local variables
  INTEGER, PARAMETER     :: maxscatter=3, allngas = 9
  INTEGER, DIMENSION(maxgkmatc), PARAMETER :: &
       greekmat_idxs = (/1, 2, 5, 6, 11, 12, 15, 16/), phasmoms_idxs = (/1, 5, 5, 2, 3, 6, 6, 4/)

  INTEGER :: i, j, k, q, nscatter, idx, cldidx, aeridx, nactgksec, nactgkmatc
  INTEGER, DIMENSION(allngas) :: absin
  REAL (KIND=dp)  :: scaco_r, absco_r, &  ! raycof, depol, 
             extco_r, extco, scaco, pvar, extco_a, scaco_a, extco_c, scaco_c, j0, j1
  REAL (KIND=dp), DIMENSION(maxscatter)              :: scaco_input
  REAL (KIND=dp), DIMENSION(ngas, nlayers)           :: absod
  REAL (KIND=dp), DIMENSION(0:maxmoments_input, 1:maxgksec, maxscatter),     SAVE :: phasmoms_input
  REAL (KIND=dp), DIMENSION(0:maxmoments_input, 1:maxgksec),                 SAVE :: phasmoms_total_input
  REAL (KIND=dp), DIMENSION(0:max_atmoswfs, 0:maxmoments_input, 1:maxgksec), SAVE :: l_phasmoms_total_input
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
     
     IF (.NOT. useasy .AND. nmoms < ngreek_moments_input) THEN
        WRITE(www_lun, *) 'Need to increase input moments for aerosols/clouds!!!'
        problems = .TRUE.; RETURN
     ENDIF

     ! This only needs to be initialized once
     phasmoms_input        = ZERO
     phasmoms_total_input  = ZERO
     greekmat_total_input  = ZERO
				   
     l_deltau_vert_input   = ZERO
     l_omega_total_input   = ZERO
     l_greekmat_total_input= ZERO 
     l_phasmoms_total_input= ZERO     
     first =.FALSE.
  ENDIF
  
  IF (NSTOKES == 1) THEN
     nactgksec = 1;  nactgkmatc = 1
  ELSE
     nactgksec = ngksec; nactgkmatc = ngkmatc
  ENDIF

  !WRITE(www_lun, *) nmoms, nmoments, ngreek_moments_input, maxmoments 
  !WRITE(www_lun, *) SUM(gascol), SUM(airgrid)
  !WRITE(www_lun, *) varyprof(1), varyprof(nlayers)
  !WRITE(www_lun, *) abscrs(1, 1), abscrs(1, 30)

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
  taugrid_input(0) = ZERO 

  ! Get rayleigh scattering phase function moments (Same for each layer)
  ! unassigned elements have already initialized to zero
  phasmoms_input(0, 1, 1) = ONE
  phasmoms_input(2, 1, 1) = (ONE - depol) / (TWO + depol)  
  IF (nactgksec == 6) THEN
     phasmoms_input(2, 2, 1) = 6.0D0 * phasmoms_input(2, 1, 1)
     phasmoms_input(2, 5, 1) = -SQRT(6.0D0) * phasmoms_input(2, 1, 1)
     phasmoms_input(1, 4, 1) = 3.0D0 * (ONE - 2.0D0 * depol) / (TWO + depol)
  ENDIF
 
  DO i = 1, nlayers   
     ! Rayleigh scattering
     scaco_r = raycof * airgrid(i)
     delray(i) = scaco_r

     ! Gas absorption
     absod(1:ngas, i) = abscrs(1:ngas, i) * gascol(1:ngas, i)
     absco_r = SUM(absod(1:ngas, i))         
     !if ( i == nlayers) print *, absco_r, absod(1:ngas, i)
     extco_r = absco_r + scaco_r
     scaco_input(1) = scaco_r
   
     !IF (absco_r <= 0.0) THEN
     !   print *, 'Negative total absorption (lamda,layer, ngas)', lamda, i, ngas
     !   write(*,'(100e15.7)') lamda,abscrs(1:ngas, i), gascol(1:ngas, i)
     !   problems  = .TRUE.
     !   RETURN
     !ENDIF

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
              DO j = 1, ngreek_moments_input
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
              DO j = 1, ngreek_moments_input
                 j1 = REAL(2*j+1, KIND=dp)
                 phasmoms_input(j, 2, nscatter) = (j1/j0) * aerasy(i) * phasmoms_input(j-1, 2, nscatter)
                 j0 = j1
              ENDDO
           ENDIF
        ENDIF  ! end aerosols
     !ENDIF     ! end non-rayleigh
     
     ! setup LIDORT input for tau and omega
     scaco = SUM(scaco_input(1:nscatter))
     omega_total_input(i) = scaco / extco
     
     IF (omega_total_input(i) < OMEGA_SMALLNUM) omega_total_input(i) = OMEGA_SMALLNUM 
     IF (omega_total_input(i) > 1.0 - OMEGA_SMALLNUM)  &
          omega_total_input(i) = 1.0 - OMEGA_SMALLNUM
     taugrid_input(i) = taugrid_input(i-1) + extco
     deltau_vert_input(i) = extco
     !extconf(i) = extco / (zsgrid(i-1) - zsgrid(i))   ! extinction coefficients
     !IF (i > 25) THEN
     !   print *, i, aermsk(i), taugrid_input(i), omega_total_input(i)
     !   print *, extco, absco_r, scaco_r, extco_a, scaco_a
     !ENDIF
  
     ! sum up phase moments as required in LIDORT
     DO j = 0, ngreek_moments_input
        DO k = 1, nactgksec
           phasmoms_total_input(j, k) = SUM(phasmoms_input(j, k, 1:nscatter) &
                * scaco_input(1:nscatter)) / scaco
        ENDDO
     ENDDO
     !phasmoms_total_input(ngreek_moments_input+1:maxmoments, 1:maxgksec) = 0.0  
     
     ! Set up greek scattering matrix for each moment 
     !greekmat_total_input(0:ngreek_moments_input, i, 1)  = phasmoms_total_input(0:ngreek_moments_input, 1)
     !greekmat_total_input(0:ngreek_moments_input, i, 2)  = phasmoms_total_input(0:ngreek_moments_input, 5)
     !greekmat_total_input(0:ngreek_moments_input, i, 5)  = phasmoms_total_input(0:ngreek_moments_input, 5)
     !greekmat_total_input(0:ngreek_moments_input, i, 6)  = phasmoms_total_input(0:ngreek_moments_input, 2)
     !greekmat_total_input(0:ngreek_moments_input, i, 11) = phasmoms_total_input(0:ngreek_moments_input, 3)
     !greekmat_total_input(0:ngreek_moments_input, i, 12) = phasmoms_total_input(0:ngreek_moments_input, 6)
     !greekmat_total_input(0:ngreek_moments_input, i, 15) = -phasmoms_total_input(0:ngreek_moments_input, 6)
     !greekmat_total_input(0:ngreek_moments_input, i, 16) = phasmoms_total_input(0:ngreek_moments_input, 4)
     !greekmat_total_input(ngreek_moments_input+1:maxmoments, i, 1:MAXSTOKES_SQ) = 0.0
     greekmat_total_input(0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
          phasmoms_total_input(0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
     IF ( nactgkmatc > 1 ) greekmat_total_input(0:ngreek_moments_input, i, 15) &
          = -greekmat_total_input(0:ngreek_moments_input, i, 15)

     ! This should always be 1, but may be slightly different due to numerical truncation
     greekmat_total_input(0, i, 1) = 1.0  

     !IF  (i == 26) THEN
     !WRITE (91, '(I5, 7D24.12, I5)') i, extco, scaco, scaco_a, absco_r, scaco_r, &
     !     omega_total_input(i), greekmat_total_input(0, i, 1), nscatter
     !DO k = 1, 16
     !   WRITE (91, '(1000D24.12)') (greekmat_total_input(j, i, k), j=0, ngreek_moments_input)
     !ENDDO
     !ENDIF

     !IF (do_simulation_only .OR. .NOT. varyprof(i)) THEN   ! no linearition
     !   ! zero out quantity for safety
     !   layer_vary_flag(i) = .FALSE.
     !   layer_vary_number(i) = 0	
	 !   	   
     !   !l_deltau_vert_input(:, i) =   ZERO
     !   !l_omega_total_input(:, i) = ZERO
     !   !l_greekmat_total_input(:, : , :, i) = ZERO        
     !
     !ELSE
     !   layer_vary_flag(i) = .TRUE.
     !   layer_vary_number(i) = n_totalatmos_wfs
     !The above part has been taken care of in routine: lidort_prof_env.f90
        
     DO q = 1, n_totalatmos_wfs        
        !  w.r.t ozone volume mixing ratio: 1
        !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF ( profilewf_names(q) == 'ozone volume mixing ratio------' ) THEN
           idx = absin(1)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, 0:maxmoments , i, 1:16) = ZERO
           
           !  w.r.t NO2 volume mixing ratio: 2
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'NO2 volume mixing ratio------' ) THEN
           idx = absin(2)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O2 volume mixing ratio: 8
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O2 volume mixing ratio------' ) THEN
           idx = absin(8)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O4 volume mixing ratio: 3
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O4 volume mixing ratio------' ) THEN
           idx = absin(3)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.;  RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t BrO volume mixing ratio: 4
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'BrO volume mixing ratio------' ) THEN
           idx = absin(4)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t H2O volume mixing ratio: 9
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'H2O volume mixing ratio------' ) THEN
           idx = absin(9)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t SO2 volume mixing ratio: 5
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'SO2 volume mixing ratio------' ) THEN
           idx = absin(5)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t HCHO volume mixing ratio: 6
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'HCHO volume mixing ratio------' ) THEN
           idx = absin(6)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t OCLO volume mixing ratio: 7
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'OCLO volume mixing ratio------' ) THEN
           idx = absin(7)
           IF (idx < 1) THEN
              WRITE(www_lun, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average temperature of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
           !  Assume no effects on air density              
        ELSE IF ( profilewf_names(q) == 'average temperature of layer---' ) THEN
           l_omega_total_input(q, i) = - SUM(absod(:, i) * eta(:, i)) / extco
           l_deltau_vert_input(q, i) = + SUM(absod(:, i) * eta(:, i)) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average pressure of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
        ELSE IF ( profilewf_names(q) == 'average pressure of layer------' ) THEN
           
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
           l_omega_total_input(q,i) = (1.0 - omega_total_input(i)) * scaco_r / extco / omega_total_input(i)
           l_deltau_vert_input(q,i) = pvar
           l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
           DO j = 0, ngreek_moments_input
             DO k = 1, nactgksec
                IF (phasmoms_total_input(j, k) /= 0.0) THEN
                   l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, 1) - phasmoms_total_input(j, k) ) &
                        / phasmoms_total_input(j, k) * scaco_a / scaco
                ELSE
                   l_phasmoms_total_input(q, j, k) = 0.0
                ENDIF
             ENDDO
          ENDDO
          l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
               l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
          IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
               = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)          
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
              l_omega_total_input(q,i) = (scaco_a  / extco_a - omega_total_input(i) ) &
                   / omega_total_input(i)  * extco_a / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
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
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
              !print *, maxval(l_greekmat_total_input), minval(l_greekmat_total_input)
              !print *, i, l_deltau_vert_input(q,i), l_omega_total_input(q,i)
              !WRITE(www_lun, '(6D14.6)') (l_phasmoms_total_input(1, j, 1:nactgksec), j = 0, ngreek_moments_input)
              !STOP
           ENDIF
           !  w.r.t  aerosol scattering coefficient / single scattering albedo
           !  aerosol optical thickness will not change
        ELSE IF ( profilewf_names(q) == 'aerosol scattering coefficient-' ) THEN
           IF (aeridx > 0) THEN
              l_deltau_vert_input(q,i) = ZERO
              l_omega_total_input(q,i) = scaco_a  / omega_total_input(i) / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
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
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
           ENDIF
           !   w.r.t cloud extinction coefficient / optical thickness
        ELSE IF ( profilewf_names(q) == 'cloud extinction coefficient---' ) THEN
           IF (cldidx > 0) THEN
              l_deltau_vert_input(q,i) = + extco_c / extco
              l_omega_total_input(q,i) = (scaco_c  / extco_c - omega_total_input(i) ) &
                   / omega_total_input(i)  * extco_c / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
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
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(0, 1:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
           ENDIF
           
           !  w.r.t clouds scattering coefficient
        ELSE IF ( profilewf_names(q) == 'cloud scattering coefficient---' ) THEN
           IF (cldidx > 0) THEN
              l_deltau_vert_input(q,i) = ZERO
              l_omega_total_input(q,i) = scaco_c  / omega_total_input(i) / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
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
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
           ENDIF
        ENDIF        ! end selection of weighting function 
     ENDDO           ! n_totalatmos_wfs loop
  !ENDIF             ! end of do_linearization    
  ENDDO              ! layer loop
  
  deltau(1:nlayers) = deltau_vert_input(1:nlayers)
  delsca(1:nlayers) = deltau_vert_input(1:nlayers) * omega_total_input(1:nlayers)
  delo3abs(1:nlayers) = absod(1, 1:nlayers)
    
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

  !SUBROUTINE GET_SLANT_TAU(nz, zs, tauin, sza, tauout)
  !USE OMSAO_parameters_module, ONLY  : maxchlen, rearth
  !USE OMSAO_precision_module
  !
  !IMPLICIT NONE
  !
  !!===============================  Define Variables ===========================
  !! Include files of dimensions and numbers
  !INCLUDE 'VLIDORT.PARS'
  !INCLUDE 'VLIDORT_INPUTS.VARS'
  !INCLUDE 'VLIDORT_SETUPS.VARS'
  !
  !! =======================
  !! Input/Output variables
  !! =======================
  !INTEGER, INTENT(IN)                          :: nz
  !REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, tauin
  !REAL (KIND=dp), INTENT(IN)                   :: sza
  !REAL (KIND=dp), DIMENSION(0:nz), INTENT(out) :: tauout
  !
  !! =======================
  !! Local variables
  !! ======================= 
  !INTEGER                   :: i, j
  !LOGICAL                   :: fail
  !CHARACTER (len=maxchlen)  :: message, trace
  !
  !n_szangles = 1; szangles(1) = sza
  !IF (sza >= 90.0 .OR. sza < 0) THEN
  !   STOP 'GET_SLANT_TAU: SZA is >= 90 or < 0!!!'
  !ENDIF
  !
  !nlayers = nz; taugrid_input(0:nz) = tauin
  !height_grid(0:nz) = zs; earth_radius = rearth
  !IF (nz > maxlayers) THEN
  !   STOP 'LIDORT_PROF_ENV: # of layers exceeded allowed !!!'
  !ENDIF 
  !
  !CALL VLIDORT_CHAPMAN(fail, message, trace)
  !tauout = 0.0
  !DO i = 1, nz
  !   DO j = 1, i
  !      tauout(i) = tauout(i) + deltau_slant(i, j, 1)
  !   ENDDO
  !ENDDO
  !
  !RETURN
  !END SUBROUTINE GET_SLANT_TAU



  !  Unused
  !
  !  ! ===============================================================
  !  ! Modified from provided routine in LODORT V23 by Rob (F77->F90)
  !  ! This routine is not consistent with CHAPMAN_FUNCTION, maybe sth.
  !  ! is wrong
  !  ! ===============================================================
  !
  !  SUBROUTINE PREPARE_SPHERICAL (nlayers, do_plane_parallel,  input_sunzen, &
  !       re, ext, z_grid, tauthick_input, sunlocal_input )
  !
  !    USE OMSAO_precision_module
  !    USE OMSAO_parameters_module, ONLY : deg2rad
  !    IMPLICIT NONE
  !
  !    !  Input arguments
  !    INTEGER, INTENT(IN)	:: nlayers
  !    LOGICAL, INTENT(IN)	:: do_plane_parallel
  !    REAL (KIND=dp), INTENT(IN)   :: input_sunzen, re
  !    REAL (KIND=dp), DIMENSION (0:nlayers), INTENT(IN) :: z_grid
  !    REAL (KIND=dp), DIMENSION (nlayers), INTENT(IN)   :: ext
  !
  !    !  Output
  !    REAL (KIND=dp), DIMENSION (nlayers, nlayers), INTENT(OUT) :: tauthick_input
  !    REAL (KIND=dp), DIMENSION (nlayers), INTENT(OUT)          :: sunlocal_input
  !
  !    !  local variables
  !    INTEGER        :: n, j, m
  !    REAL (KIND=dp) :: gm_toa, th_toa, th0, th1, gm0, gm1
  !    REAL (KIND=dp) :: h(0:nlayers), delz, taup, mu_toa
  !    REAL (KIND=dp) :: z_dIff, z_0, z, const0
  !    REAL (KIND=dp) :: x, xd, cumdep, s2, hf, deltm, delt
  !
  !
  !    !  get spherical optical depths
  !    !  ----------------------------
  !    !  prepare spherical attenuation (shell geometry)
  !
  !    IF ( .NOT. do_plane_parallel ) THEN
  !
  !      mu_toa = COS ( input_sunzen * deg2rad )
  !      gm_toa = SQRT ( 1.0d0 - mu_toa * mu_toa )
  !      th_toa = ASIN (gm_toa)
  !      h(0:nlayers) = z_grid(0:nlayers) + re
  !      const0 = gm_toa / h(0)
  !      cumdep = 0.0d0
  !
  !      DO n = 1, nlayers
  !        delz = z_grid(n-1)-z_grid(n)
  !        z_diff   = delz
  !        z_0 = z_grid(n-1)
  !        z  = z_0
  !        DO j = 1, 1
  !          z = z - z_dIFf
  !          x = z_0 - z
  !          xd = x + cumdep
  !          hf = h(0) - xd
  !          gm0 = const0 * hf
  !          th0 = ASIN ( gm0 )
  !          taup = 0.0d0
  !          DO m = 1, n-1
  !            gm1 = h(m-1) * gm0 / h(m)
  !            th1 = ASIN(gm1)
  !            s2 = h(m-1) * SIN(th1-th0) / gm1
  !            IF ( j == 1 ) tauthick_input(n,m) = ext(m) * s2
  !            taup = taup + ext(m) * s2
  !            th0 = th1
  !            gm0 = gm1
  !          ENDDO
  !          s2 = h(n-1)* SIN(th_toa-th0) / gm_toa
  !          taup = taup + ext(n) * s2
  !          IF ( j== 1 ) tauthick_input(n,n) = ext(n) * s2
  !        ENDDO
  !        cumdep = xd
  !        sunlocal_input(n) = mu_toa
  !      ENDDO
  !    ELSE
  !
  !      mu_toa = COS( input_sunzen * deg2rad )
  !      DO n = 1, nlayers
  !        sunlocal_input(n) = mu_toa
  !        delz = z_grid(n-1)-z_grid(n)
  !        delt = ext(n) * delz
  !        tauthick_input(n,n) = delt/mu_toa
  !        DO m = 1, n-1
  !          deltm = ext(m) * delz
  !          tauthick_input(n,m) = deltm / mu_toa
  !        ENDDO
  !      ENDDO
  !
  !    ENDIF
  !
  !  END SUBROUTINE PREPARE_SPHERICAL


END MODULE m_lidort_master

